#install.packages("purrr")
library(readr)
library(dplyr)
library(ggplot2)
library(purrr)

setwd("WORKING/DIR")
files <- list.files(path = "Part1BAMfiles", pattern = "dup_metrics\\.txt$", full.names = TRUE)

# Function to read the histogram table
read_histogram_table <- function(file_path) {
  lines <- readLines(file_path)
  
  # Find the line that says "## HISTOGRAM	java.lang.Double"
  hist_line <- grep("^## HISTOGRAM", lines)
  
  if (length(hist_line) == 0) {
    message("Skipping (no histogram found): ", basename(file_path))
    return(tibble())
  }
  
  # Read the table from the next line onward
  hist_df <- read_tsv(
    file_path,
    skip = hist_line,    # skip includes the "## HISTOGRAM" line itself
    comment = "#",        # ignore any additional comments
    show_col_types = FALSE
  )
  
  hist_df$File <- basename(file_path)
  return(hist_df)
}

# Read all histograms into one dataframe 
all_histograms <- map_dfr(files, read_histogram_table)
sampleNames <- unique(all_histograms$File)

coverage_summary <- data.frame(Sample = character(), AverageCoverage = numeric(), stringsAsFactors = FALSE)

# Loop through each sample
for (i in 1:length(sampleNames)) {
  tempDF <- all_histograms[all_histograms$File == sampleNames[i],]
  
  # Compute average (or max if you prefer) coverage
  averageCoverage <- mean(tempDF$CoverageMult, na.rm = TRUE)
  # Or use max:
  # averageCoverage <- max(tempDF$CoverageMult, na.rm = TRUE)
  
  # Save to summary table
  coverage_summary <- rbind(coverage_summary, data.frame(Sample = sampleNames[i], AverageCoverage = averageCoverage))
  
  # Create the plot
  plotTitle <- paste0("Average genome coverage = ", round(averageCoverage, 2))
  plot <- ggplot(data = tempDF, aes(x = BIN, y = all_sets)) + 
    geom_bar(stat = "identity") + 
    labs(title = plotTitle) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Clean up filename
  cleanName <- gsub("[[:space:]]", "_", sampleNames[i])
  plotFileName <- paste0("Extra_information/genome_coverage_", cleanName, ".png")
  
  # Save the plot
  ggsave(plotFileName, plot = plot, width = 5, height = 2)
}

write.csv(coverage_summary, "Extra_information/average_coverage_per_sample.csv", row.names = FALSE)

