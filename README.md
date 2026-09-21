# PDAC StateScope deconvolution

Publication-facing code and computational provenance for the StateScope-based deconvolution of bulk RNA-seq data from pancreatic ductal adenocarcinoma (PDAC).

## Workflow

```text
RNA-seq FASTQ
  -> 01_align_rnaseq_star.sh
  -> 02_build_statescope_inputs.R ----------------------+
                                                        |                                                                                            
DNA-seq FASTQ                                           |
  -> 03_align_dna_fastq.sh   
  -> 04_run_purity_grid.sh                              |
  -> 05_run_cnvkit_wgs.sh                               |
  -> 06_select_best_purity.py --------------------------+
                                                        |
                                                        v
                                             07_run_statescope.py
                                               deconvolution
                                               refinement x2
                                               state discovery
                                                        |
                                                        v
                                          08_export_state_fractions.py
```

## Key software
- STAR 2.7.11b
- GENCODE v46 / GRCh38
- CNVkit 0.9.12
- scikit-learn 1.7.1 (exploratory GMM purity script)
- StateScope 1.0.6

## StateScope settings
- `TumorType="PDAC"`
- `Ncores=16`
- CNVkit-derived prior supplied only to `Epithelial`
- prior capped at 0.99
- two sequential `Refinement()` calls
- `StateDiscovery()` defaults from the original v1.0.6 environment: Omega weighting, K=2–9, 10 initial restarts/K, cophenetic threshold 0.9, 100 final restarts

## Biological state labels
Biological annotation was performed separately and is described in the manuscript rather than encoded in this repository.
