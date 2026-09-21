# Computational workflow

## RNA sequencing
Paired-end RNA-seq FASTQs were aligned to GRCh38 with STAR v2.7.11b in two-pass mode using the GENCODE v46 GRCh38 reference FASTA and matching GTF. STAR `ReadsPerGene.out.tab` files were combined in R. Forward versus reverse stranded counts were compared per sample and the final stranded count column was chosen by majority vote across samples. Ensembl version suffixes were removed adn ID's were mapped to GENCODE v46 gene symbols. When duplicate symbols were rpesent, these were summed, and missing values were set to zero. No additional normalization was applied before creation of the raw count matrix as StateScope uses the raw count input.

## DNA sequencing and CNVkit
Paired-end DNA reads were aligned to hg19 using BWA-MEM, coordinate-sorted with samtools, supplied with read groups and duplicate-marked using Picard. CNVkit v0.9.12 was run in WGS mode. Candidate purity values from 0.30 to 1.00 in increments of 0.05 were evaluated. For each candidate, expected purity-adjusted log2 ratios were compared with observed segment log2 ratios by segment-weighted RMSE; the minimum-RMSE purity was selected.

## StateScope
StateScope v1.0.6 was initialized with `TumorType="PDAC"` and `Ncores=16`. CNVkit-derived purity was supplied only to `Epithelial` and capped at 0.99. Deconvolution was followed by two sequential `Refinement()` calls. State discovery used the verified v1.0.6 defaults: Omega weighting, automatic K testing from 2 to 9, 10 initial cNMF restarts per K, minimum cophenetic coefficient 0.9, and 100 final restarts.

