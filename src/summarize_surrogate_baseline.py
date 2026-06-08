#!/usr/bin/env python3
"""src/summarize_surrogate_baseline.py — Step 2 summary (numbers only, no interpretation).

Reads results/step2/surrogate_vs_measured.csv and computes, split by surrogate scheme:
  1. Regression on the NUMERIC rows: Spearman rho in log space (predicted Aff vs measured IC50).
     Spearman is rank-based, so log10 is monotone-invariant -> identical rho; computed on log10
     to honour the "log space" instruction. n is small (pilot) and reported as such.
  2. Binary binder/non-binder classification using the assay's OWN dichotomy:
       truth: binder = measured IC50 within assay range (numeric); non-binder = censored (>40,000 nM).
       predicted-binder threshold (stated, not optimised): EL %rank <= 2.0 (primary) and Aff <= 500 nM.
     Reports the full 2x2 + accuracy, and the 186 censored rows specifically (TN/FP).
  3. The 135 rows where naive-surrogate == native ("assume the modification does nothing"):
     residual distribution on the numeric subset + how often it predicts binding for a true
     non-binder on the censored subset.

Stdlib only. Writes results/surrogate_baseline.csv (tidy metrics) + results/surrogate_baseline_summary.txt.
"""
import csv, os, math, statistics as st

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC  = os.path.join(REPO, "results", "step2", "surrogate_vs_measured.csv")
OUTC = os.path.join(REPO, "results", "surrogate_baseline.csv")
OUTT = os.path.join(REPO, "results", "surrogate_baseline_summary.txt")

EL_RANK_BINDER = 2.0     # canonical NetMHCpan binder threshold (%rank)
AFF_BINDER_NM  = 500.0   # standard high-affinity IC50 threshold

def _avg_ranks(a):
    n = len(a); order = sorted(range(n), key=lambda i: a[i]); r = [0.0]*n; i = 0
    while i < n:
        j = i
        while j+1 < n and a[order[j+1]] == a[order[i]]: j += 1
        avg = (i + j)/2.0 + 1.0
        for k in range(i, j+1): r[order[k]] = avg
        i = j+1
    return r

def _pearson(x, y):
    n = len(x); mx = sum(x)/n; my = sum(y)/n
    cov = sum((a-mx)*(b-my) for a, b in zip(x, y))
    sx = math.sqrt(sum((a-mx)**2 for a in x)); sy = math.sqrt(sum((b-my)**2 for b in y))
    return cov/(sx*sy) if sx > 0 and sy > 0 else float("nan")

def spearman(x, y):
    return _pearson(_avg_ranks(x), _avg_ranks(y))

rows = list(csv.DictReader(open(SRC)))
schemes = ["naive_parent_map", "charge_aware"]
metrics = []   # (analysis, scheme, n, metric, value)
def add(a, s, n, m, v): metrics.append((a, s, n, m, v))

lines = ["STEP 2 — SURROGATE-ERROR SUMMARY (numbers only; no interpretation)",
         "=" * 70,
         "predicted = NetMHCpan-4.1 of the surrogate peptide; measured = modified-peptide IC50.",
         "log10_resid = log10(pred_Aff_nM) - log10(measured_IC50_nM)  [+ => pred weaker; - => stronger]",
         ""]

# ---------- 1. Regression (numeric rows) ----------
lines.append("1) REGRESSION (numeric measured rows) — Spearman rho in log space")
for s in schemes:
    sub = [r for r in rows if r["scheme"] == s and r["measured_censored"] == "False"]
    lp = [math.log10(float(r["pred_Aff_nM"])) for r in sub]
    lm = [math.log10(float(r["measured_modified_IC50_nM"])) for r in sub]
    rho = spearman(lp, lm)
    add("regression_logspace", s, len(sub), "spearman_rho", round(rho, 3))
    lines.append(f"   {s:16}  n={len(sub):3}  Spearman rho = {rho:+.3f}   [PILOT: n small]")
lines.append("")

# ---------- 2. Classification (truth = assay dichotomy: numeric=binder, censored=non-binder) ----------
lines.append("2) CLASSIFICATION — truth: binder=measured numeric, non-binder=measured >40,000 (censored)")
for thr_name, is_pred_binder in (
        (f"EL%rank<= {EL_RANK_BINDER}", lambda r: float(r["pred_EL_rank"]) <= EL_RANK_BINDER),
        (f"Aff<= {AFF_BINDER_NM:.0f}nM", lambda r: float(r["pred_Aff_nM"]) <= AFF_BINDER_NM)):
    lines.append(f"  predicted-binder rule: {thr_name}")
    for s in schemes:
        sub = [r for r in rows if r["scheme"] == s]
        TP = FP = TN = FN = 0
        for r in sub:
            truth_binder = (r["measured_censored"] == "False")
            pb = is_pred_binder(r)
            if   truth_binder and pb:         TP += 1
            elif truth_binder and not pb:     FN += 1
            elif (not truth_binder) and pb:   FP += 1
            else:                             TN += 1
        n = len(sub); acc = (TP+TN)/n
        ncens = FP + TN  # actual non-binders (censored)
        add("classification", s, n, f"{thr_name}|accuracy", round(acc, 3))
        add("classification", s, n, f"{thr_name}|TP", TP); add("classification", s, n, f"{thr_name}|FP", FP)
        add("classification", s, n, f"{thr_name}|TN", TN); add("classification", s, n, f"{thr_name}|FN", FN)
        lines.append(f"     {s:16}  acc={acc:.3f}  [TP={TP} FN={FN} | FP={FP} TN={TN}]  "
                     f"censored(non-binders) n={ncens}: TN={TN} correct, FP={FP} called binder")
    lines.append("")

# ---------- 3. "assume the modification does nothing" : naive==native rows ----------
lines.append("3) ASSUME-NO-EFFECT (rows where naive-surrogate == native; naive scheme)")
nn = [r for r in rows if r["scheme"] == "naive_parent_map" and r["surrogate_eq_native"] == "True"]
nn_num = [r for r in nn if r["measured_censored"] == "False"]
nn_cen = [r for r in nn if r["measured_censored"] == "True"]
resid = [abs(float(r["log10_resid_pred_minus_meas"])) for r in nn_num]
lines.append(f"   rows where naive==native: {len(nn)}  ({len(nn_num)} numeric, {len(nn_cen)} censored)")
if resid:
    le03 = sum(1 for x in resid if x <= 0.3); l03_1 = sum(1 for x in resid if 0.3 < x <= 1.0)
    gt1 = sum(1 for x in resid if x > 1.0)
    lines.append(f"   |log10_resid| on the {len(nn_num)} numeric: median={st.median(resid):.2f}  "
                 f"max={max(resid):.2f}")
    lines.append(f"     within 2-fold (|resid|<=0.3): {le03}   2-10x (0.3-1.0): {l03_1}   "
                 f">10-fold off (|resid|>1.0): {gt1}")
    add("assume_no_effect", "naive_parent_map", len(nn_num), "median_abs_log10_resid", round(st.median(resid), 3))
    add("assume_no_effect", "naive_parent_map", len(nn_num), "n_gt_10fold_off", gt1)
fp_cen = sum(1 for r in nn_cen if float(r["pred_EL_rank"]) <= EL_RANK_BINDER)
fp_cen_aff = sum(1 for r in nn_cen if float(r["pred_Aff_nM"]) <= AFF_BINDER_NM)
lines.append(f"   censored subset ({len(nn_cen)} true non-binders): naive-surrogate predicts BINDER for "
             f"{fp_cen} (EL%rank<=2) / {fp_cen_aff} (Aff<=500nM)")
add("assume_no_effect", "naive_parent_map", len(nn_cen), "censored_predicted_binder_ELrank", fp_cen)

# ---------- write ----------
with open(OUTC, "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n"); w.writerow(["analysis", "scheme", "n", "metric", "value"])
    w.writerows(metrics)
with open(OUTT, "w", newline="\n") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
print(f"\nwrote {OUTC}  and  {OUTT}")
