suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})

# 02_build_statescope_inputs.R
# Build the raw gene-count matrix and purity table consumed by StateScope.

star_root <- "FOLDER/STAR_out"
gtf_path <- "REFgencode/GRCh38/gencode.v46.annotation.gtf"
purity_tsv <- "/CNVkit_output/purity_calls/purity_summary.tsv"
out_dir <- "/output/folder/"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

##Functions
read_star_tab <- function(f) {
  x <- suppressMessages(readr::read_tsv(f, col_names = FALSE, show_col_types = FALSE))
  colnames(x) <- c("GeneID","Unstranded","Forward","Reverse")
  x %>% filter(!str_starts(GeneID, "N_"))
}
files <- list.files(star_root, pattern = "ReadsPerGene\\.out\\.tab$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0) stop("No ReadsPerGene.out.tab files found under: ", star_root)
samples <- basename(dirname(files))

pick_col_per_sample <- function(f) {
  x <- read_star_tab(f)
  sums <- colSums(x[, c("Forward","Reverse")], na.rm = TRUE)
  if (is.na(sums[1]) || is.na(sums[2])) return(4L)
  if (sums["Forward"] > sums["Reverse"]) 3L else 4L
}
cols_voted <- vapply(files, pick_col_per_sample, integer(1))
selected_col <- as.integer(names(sort(table(cols_voted), decreasing = TRUE)[1]))
message("Selected STAR counts column = ", selected_col,
        " (2=unstranded, 3=forward, 4=reverse). Votes: ",
        paste0(names(table(cols_voted)), "→", as.integer(table(cols_voted)), collapse = ", "))

read_one_counts <- function(f, sel_col) {
  x <- read_star_tab(f) %>% transmute(
    gene_id_clean = sub("\\.\\d+$","", GeneID),
    Count = as.numeric(.data[[c("Unstranded","Forward","Reverse")][sel_col - 1]]))
  tibble(gene_id_clean = x$gene_id_clean, Count = x$Count)
}

gtf <- suppressMessages(readr::read_tsv(
  gtf_path, comment = "#",
  col_names = c("seqname","source","feature","start","end","score","strand","frame","attribute"),
  col_types = "ccccccccc"))
gene_attr <- gtf %>% filter(feature == "gene") %>% transmute(attribute)
extract_attr <- function(attr, key) {
  m <- stringr::str_match(attr, paste0(key, " \\\"([^\\\"]+)\\\""))
  m[,2]
}
gene_map <- gene_attr %>%
  mutate(gene_id = extract_attr(attribute, "gene_id"), gene_name = extract_attr(attribute, "gene_name")) %>%
  select(gene_id, gene_name) %>% distinct() %>%
  mutate(gene_id_clean = sub("\\.\\d+$","", gene_id)) %>%
  select(gene_id_clean, gene_name) %>% distinct()
if (nrow(gene_map) == 0) stop("Could not parse gene_id/gene_name from GTF: ", gtf_path)

counts_list <- map(files, ~ read_one_counts(.x, selected_col)); names(counts_list) <- samples
counts <- counts_list %>%
  imap(~ rename(.x, !!.y := Count)) %>% reduce(full_join, by = "gene_id_clean") %>%
  left_join(gene_map, by = "gene_id_clean") %>%
  mutate(Gene = ifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id_clean)) %>%
  select(-gene_name, -gene_id_clean) %>% group_by(Gene) %>%
  summarise(across(everything(), ~ as.numeric(replace_na(., 0))), .groups = "drop")
readr::write_csv(counts, file.path(out_dir, "gene_counts_matrix.csv"))

pur <- suppressMessages(readr::read_tsv(purity_tsv, show_col_types = FALSE))
nms <- names(pur)
sample_col <- nms[str_detect(tolower(nms), "sample|^id$|sample_id|tumor|name")][1]
purity_candidates <- nms[str_detect(tolower(nms), "purity")]
pur_col <- if (length(purity_candidates) == 0) NA_character_ else {
  abs_idx <- which(str_detect(tolower(purity_candidates), "absolute"))
  if (length(abs_idx)) purity_candidates[abs_idx[1]] else purity_candidates[1]
}
if (is.na(sample_col) || is.na(pur_col)) stop("Could not identify sample/purity columns")
pur_out <- pur %>% transmute(Sample = as.character(.data[[sample_col]]), `ABSOLUTE Purity` = as.numeric(.data[[pur_col]])) %>% distinct() %>% arrange(Sample)
readr::write_csv(pur_out, file.path(out_dir, "tumor_purity_absolute.csv"))
