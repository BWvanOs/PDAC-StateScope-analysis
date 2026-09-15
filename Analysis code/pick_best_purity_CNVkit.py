#!/usr/bin/env python3
import os, glob, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

indir = "/INPUT/DIR"
purdir = os.path.join(indir, "purity_calls")
plotdir = os.path.join(purdir, "purity_curves")
os.makedirs(plotdir, exist_ok=True)

def score(original_cns, called_cns, p):
    """Compute weighted RMSE between observed log2 ratios and expected given purity p."""
    base = pd.read_table(original_cns)
    called = pd.read_table(called_cns)

    Lobs = base['log2'].values
    w = base['weight'].values
    cn = called['cn'].astype(float).values

    expected_cn = p * cn + (1 - p) * 2.0
    eps = 1e-6
    expected_cn = np.clip(expected_cn, eps, None)
    Ltheory = np.log2(expected_cn / 2.0)

    if Lobs.shape[0] != Ltheory.shape[0]:
        raise ValueError(f"Length mismatch for {original_cns} vs {called_cns}")

    rmse = math.sqrt(np.average((Lobs - Ltheory) ** 2, weights=w))
    return rmse

records = []
for cns in sorted(glob.glob(os.path.join(indir, "*.final.cns"))): 
    sample = os.path.basename(cns).replace(".final.cns","")
    purity_vals, rmse_vals = [], []
    for cand in sorted(glob.glob(os.path.join(purdir, f"{sample}.final.purity*.cns"))):
        part = cand.split("purity")[-1].replace(".cns","")
        try:
            p = float(part)
        except ValueError:
            continue
        try:
            rmse = score(cns, cand, p)
            purity_vals.append(p)
            rmse_vals.append(rmse)
        except Exception as e:
            print(f"Skipping {cand}: {e}")
            continue
    if purity_vals:
        idx = int(np.argmin(rmse_vals))
        best_p, best_rmse = purity_vals[idx], rmse_vals[idx]
        records.append((sample, best_p, best_rmse))

        # Plot RMSE curve
        plt.figure()
        plt.plot(purity_vals, rmse_vals, marker="o", color="blue")
        plt.axvline(best_p, color="red", linestyle="--", label=f"Best={best_p:.2f}")
        plt.xlabel("Purity")
        plt.ylabel("RMSE (lower = better fit)")
        plt.title(f"Purity fit — {sample}")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(plotdir, f"{sample}.rmse_curve.png"), dpi=150)
        plt.close()

df = pd.DataFrame(records, columns=["Sample","BestPurity","RMSE"])
outf = os.path.join(purdir, "purity_summary.tsv")
df.to_csv(outf, sep="\t", index=False)

print(f"Wrote summary: {outf} with {len(df)} samples")
print(f"Individual RMSE curve plots saved in {plotdir}/")

