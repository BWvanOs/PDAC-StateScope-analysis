#!/usr/bin/env python3
import argparse, math, numpy as np, pandas as pd
from sklearn.mixture import GaussianMixture

def load_cns(path):
    df = pd.read_table(path)
    # CNVkit .cns typically has: chromosome, start, end, gene, log2, probes, weight
    # Normalize column naming we’ll use
    colMap = {c.lower(): c for c in df.columns}
    # Ensure required columns
    req = ['chromosome', 'start', 'end', 'log2']
    for r in req:
        if r not in colMap:
            raise ValueError(f"Missing column '{r}' in {path}")
    df.rename(columns={colMap['chromosome']:'Chromosome',
                       colMap['start']:'Start',
                       colMap['end']:'End',
                       colMap['log2']:'Log2'}, inplace=True)
    if 'weight' in df.columns:
        df.rename(columns={colMap['weight']:'Weight'}, inplace=True)
    else:
        df['Weight'] = (df['End'] - df['Start']).clip(lower=1)
    # Filter sex chr if desired, keep autosomes by default
    df = df[~df['Chromosome'].astype(str).str.upper().isin(['X','Y','MT','M'])]
    # Remove tiny segments and extreme outliers
    df = df[df['Weight'] > 10]  # tweak if needed
    df = df[np.isfinite(df['Log2'])].copy()
    return df

def purity_from_peak(log2center, direction):
    # direction: 'loss' (expects negative) or 'gain' (expects positive)
    L = log2center
    if direction == 'loss' and L < 0:
        p = 2*(1 - 2**L)
    elif direction == 'gain' and L > 0:
        p = 2*((2**L) - 1)
    else:
        return np.nan
    return float(np.clip(p, 0.0, 1.0))

def fit_peaks(df):
    # Weight by segment size, fit 3-component GMM to capture loss/neutral/gain
    X = df['Log2'].values.reshape(-1,1)
    W = df['Weight'].values
    k = 3 if len(df) >= 30 else 2
    gmm = GaussianMixture(n_components=k, covariance_type='full', random_state=1)
    gmm.fit(X)
    means = gmm.means_.ravel()
    # Sort components by mean
    means.sort()
    # Heuristic: lowest = loss-like, highest = gain-like, middle ~ neutral
    lossCenter = means[0]
    gainCenter = means[-1]
    return lossCenter, gainCenter

def combine_purities(pLoss, pGain):
    vals = [v for v in [pLoss, pGain] if np.isfinite(v)]
    if not vals:
        return np.nan, np.nan, np.nan
    # median + a rough 95% CI using IQR heuristic
    median = float(np.median(vals))
    lo = float(np.percentile(vals, 25))
    hi = float(np.percentile(vals, 75))
    # widen a bit to be conservative
    spread = (hi - lo)
    return median, max(0.0, median - 1.57*spread), min(1.0, median + 1.57*spread)

def main():
    ap = argparse.ArgumentParser(description="Estimate tumor purity from CNVkit .cns")
    ap.add_argument("cns", help=".cns file")
    ap.add_argument("--print-peaks", action="store_true", help="Print peak centers")
    args = ap.parse_args()

    df = load_cns(args.cns)
    if df.empty:
        raise SystemExit("No segments after filtering.")
    lossCenter, gainCenter = fit_peaks(df)

    pLoss = purity_from_peak(lossCenter, 'loss')
    pGain = purity_from_peak(gainCenter, 'gain')
    pHat, pLo, pHi = combine_purities(pLoss, pGain)

    if args.print_peaks:
        print(f"Loss peak ≈ {lossCenter:.3f}, Gain peak ≈ {gainCenter:.3f}")
    print(f"Purity ≈ {pHat:.3f} (rough CI: {pLo:.3f}–{pHi:.3f})")
    if np.isfinite(pLoss): print(f"  from losses: {pLoss:.3f}")
    if np.isfinite(pGain): print(f"  from gains : {pGain:.3f}")

if __name__ == "__main__":
    main()

