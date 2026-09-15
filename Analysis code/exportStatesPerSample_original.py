import os

# Output folder
outdir = "PATH/TO/OUTDIR"
os.makedirs(outdir, exist_ok=True)

# Loop over all cell types in StateScores
for celltype, scores_df in ss.StateScores.items():
    if not hasattr(scores_df, "div"):
        print(f"Skipping {celltype} (not a DataFrame)")
        continue

    # Fractions for this cell type
    if celltype not in ss.Fractions.columns:
        print(f"Skipping {celltype}, no matching fraction found")
        continue
    frac = ss.Fractions[celltype]

    # Normalize scores (avoid negatives and  scale rows to 1)
    scores_norm = scores_df.clip(lower=0)
    scores_norm = scores_norm.div(scores_norm.sum(axis=1), axis=0)

    # Multiply normalized scores × cell fraction
    state_fractions = scores_norm.mul(frac, axis=0)

    # Save to CSV
    outpath = os.path.join(outdir, f"{celltype}_StateFractions.csv")
    state_fractions.to_csv(outpath)
    print(f"Saved {celltype} state fractions -> {outpath}")

