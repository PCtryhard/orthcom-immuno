#!/usr/bin/env python3
"""src/parse_surrogate_baseline.py — Step 2: join NetMHCpan predictions -> manifest -> measured.

For each testable source row and each surrogate scheme (naive_parent_map / charge_aware) it joins
the predicted binding of THAT scheme's surrogate peptide to the MEASURED modified IC50, and records
a per-row surrogate error (log10 residual) for numeric rows.

JOIN INTEGRITY (the two failure modes called out):
  (a) Join is keyed on the EXACT peptide sequence (netMHCpan 'Peptide' column, NOT 'Core'), via a
      dict lookup. A sequence has one length, so 9-mers and 10-mers can never collide. Each source
      row keeps its own output row even when two source rows share a surrogate sequence (one
      prediction, both rows preserved) — never collapsed, never cross-joined.
  (b) Censored measured values (modified peptide is a non-binder, IC50 > 40,000 nM) are carried as
      the literal flag '>40000'. They are NEVER coerced to 40000 or any number; log10_resid is blank.

NO summary statistic is computed here (no Spearman, no classification metric) — table only.
"""
import csv, os, sys, math
from collections import defaultdict

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S2   = os.path.join(REPO, "results", "step2")
RAW  = os.path.join(S2, "pred_netmhcpan_el_ba_raw.txt")
MAN  = os.path.join(S2, "input_manifest.csv")
OUT  = os.path.join(S2, "surrogate_vs_measured.csv")

# ---- parse predictions: peptide -> fields (key on EXACT Peptide column f[2]) ----
preds = {}
with open(RAW) as fh:
    for line in fh:
        f = line.split()
        if len(f) >= 16 and f[0].isdigit() and f[1] == "HLA-A*02:01":   # data row
            pep = f[2]
            if pep in preds:
                sys.exit(f"FATAL: duplicate prediction for peptide {pep}")
            preds[pep] = {"Score_EL": f[11], "Rnk_EL": f[12],
                          "Score_BA": f[13], "Rnk_BA": f[14], "Aff_nM": f[15]}
print(f"predictions parsed: {len(preds)} unique peptides")

# ---- read manifest, group by source row ----
man = list(csv.DictReader(open(MAN)))
by_src = defaultdict(dict)
for r in man:
    by_src[r["source_row"]][r["variant"]] = r
print(f"manifest rows: {len(man)}  |  source rows: {len(by_src)}")

SCHEME = {"naive_surrogate": "naive_parent_map", "charge_surrogate": "charge_aware"}
rows_out, missing = [], []
for src, v in by_src.items():
    assert len(v) == 3, f"{src}: expected 3 variants, got {sorted(v)}"
    nat = v["native"]
    measured_raw = nat["measured_modified_IC50_nM"].strip()   # measured MODIFIED IC50 (same across variants)
    censored = (measured_raw == "")
    native_pred = preds.get(nat["peptide"])
    for var in ("naive_surrogate", "charge_surrogate"):
        m = v[var]; pep = m["peptide"]; p = preds.get(pep)
        if p is None:
            missing.append((src, var, pep)); continue
        aff = float(p["Aff_nM"])
        resid = "" if censored else round(math.log10(aff) - math.log10(float(measured_raw)), 3)
        rows_out.append({
            "source_row": src, "modification_type": m["modification_type"],
            "pair_status": m["pair_status"], "native_seq": nat["peptide"],
            "modified_seq": m["modified_sequence"], "scheme": SCHEME[var],
            "surrogate_peptide": pep, "surrogate_eq_native": pep == nat["peptide"],
            "pred_Aff_nM": p["Aff_nM"], "pred_EL_score": p["Score_EL"], "pred_EL_rank": p["Rnk_EL"],
            "pred_native_Aff_nM": native_pred["Aff_nM"] if native_pred else "",
            "measured_modified_IC50_nM": "" if censored else measured_raw,
            "measured_modified": ">40000" if censored else measured_raw,
            "measured_censored": censored,
            "log10_resid_pred_minus_meas": resid,
        })

if missing:
    sys.exit(f"FATAL: {len(missing)} surrogate peptides missing from predictions: {missing[:5]}")

cols = ["source_row","modification_type","pair_status","native_seq","modified_seq","scheme",
        "surrogate_peptide","surrogate_eq_native","pred_Aff_nM","pred_EL_score","pred_EL_rank",
        "pred_native_Aff_nM","measured_modified_IC50_nM","measured_modified","measured_censored",
        "log10_resid_pred_minus_meas"]
with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n"); w.writeheader(); w.writerows(rows_out)

# ---- join reconciliation ----
n_src = len(by_src)
print("\n===== JOIN RECONCILIATION =====")
print(f"source rows IN (testable, 1 measured-IC50 each): {n_src}")
print(f"comparison rows OUT: {len(rows_out)}  (= {n_src} source x 2 schemes; expect {n_src*2})")
print(f"every surrogate peptide found in predictions: {len(missing)==0}  (missing={len(missing)})")
uniq = {r['surrogate_peptide'] for r in rows_out}
print(f"distinct surrogate peptides used: {len(uniq)} of {len(preds)} predicted")
src_per_pep = defaultdict(set)
for r in rows_out: src_per_pep[r['surrogate_peptide']].add(r['source_row'])
many = {p:s for p,s in src_per_pep.items() if len(s) > 1}
print(f"many-to-one: {len(many)} surrogate peptides shared by >1 source row (BY DESIGN: identical "
      f"surrogate sequence -> identical prediction; all rows preserved, none collapsed). Examples:")
for p, s in list(many.items())[:4]:
    print(f"    {p}  <- source rows {sorted(s)}")
print(f"prediction-key lengths present: {sorted({len(p) for p in preds})} "
      f"(each sequence has one length -> no 9/10-mer key collision)")
nnum = sum(1 for r in rows_out if not r['measured_censored'])
print(f"measured: {nnum} numeric rows (log10_resid computed) | {len(rows_out)-nnum} censored '>40000' (resid blank)")
print("counts by pair_status x scheme:")
for ps in ("complete","native_nonbinder","modified_nonbinder","both_nonbinder"):
    a = sum(1 for r in rows_out if r['pair_status']==ps and r['scheme']=='naive_parent_map')
    b = sum(1 for r in rows_out if r['pair_status']==ps and r['scheme']=='charge_aware')
    print(f"    {ps:18}  naive={a:3}  charge={b:3}")
coincide = sum(1 for s in by_src if by_src[s]['naive_surrogate']['peptide']==by_src[s]['charge_surrogate']['peptide'])
print(f"source rows where naive==charge surrogate (no U/B; schemes coincide): {coincide}/{n_src} "
      f"(so naive vs charge only differs for the {n_src-coincide} U/B rows)")
print(f"\nwrote {OUT}  ({len(rows_out)} rows)")
