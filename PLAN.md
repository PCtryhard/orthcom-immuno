# PLAN.md — OrthCom-Immuno

**Project:** Structure-based assessment of how non-canonical amino acids perturb MHC-I
presentation, benchmarked against the field's canonical-surrogate shortcut.

**Status legend:** `TODO` → `PROMPTED` (sent to Claude Code, awaiting result) → `DONE`
(verified by the human). Only one step should be `PROMPTED` at a time. Do not advance a step
to `DONE` on the agent's say-so alone — verify against actual files/outputs first (this is
inherited institutional discipline from OrthCom: docs and agents have been wrong before).

---

## The thesis (one sentence)

The field's current method for judging whether a non-canonical amino acid (ncAA) creates an
immunological liability is to **replace it with the nearest canonical residue and run a
standard tool** (Mattei et al. 2022). This project tests **how wrong that shortcut is** by
modeling the modified residue's actual chemistry in the HLA-A\*02:01 groove, and produces a
**position-resolved map** of presentation liability the surrogate misses.

## Why it's novel (and the exact boundary)

- **Already done — do NOT claim:**
  - Mattei et al. 2022 (*Front. Drug Discov.*): the canonical-surrogate **workflow** for UAA
    immunogenicity (e.g. proxy 1-naphthyl-Ala as **Trp**). Explicitly a "vision for future
    enhancements" — i.e. unvalidated.
  - Jiang, Su, Meiler et al. 2025 (*PLOS One*): a single-allele (A\*02:01) **ML
    binding-affinity regressor** for ncAA-peptides (R²≈0.48). The naive "train a model"
    headline is taken.
- **Open — the defensible claim:** nobody has (a) scored the **actual** modified residue
  structurally in the MHC-I groove, (b) **quantified the error** of the canonical-surrogate
  shortcut, or (c) produced a **position-resolved liability map** for designed/bioorthogonal
  ncAAs. The contribution is the **structure-based correction to, and validation of, the
  existing sequence-surrogate paradigm.**
- **Both outcomes publish:** surrogate-is-wrong → "field's shortcut mis-ranks ncAA
  immunogenicity + here's the fix" (high impact). Surrogate-is-right → "first structural
  validation that the shortcut is justified, and the boundary where it breaks" (a real
  methods/validation result). The **only** worthless outcome is *assuming* either — which is
  the current state of the art.

## The critical de-risking correction (read before starting)

There is **no** open-literature dataset of *aromatic synthetic ncAA* (naphthyl-Ala,
fluoro-Phe) on A\*02:01 with measured binding/T-cell data. So the project is staged:

1. **Validate** the structure-based method against data that **does** exist — PTM-modified
   A\*02:01 peptides with measured IC50s (citrullination, oxidation, deamidation,
   chlorination on insulin-B; Sidney/Sette 2018, BMC Immunol, PMC5863483). These are
   chemically-defined non-canonical residues with quantitative ground truth on the target
   allele.
2. **Extend predictively** to the synthetic aromatic ncAAs (2-Nal, fluoro-Phe) where no data
   exists — framed as falsifiable predictions for a future wet-lab collaboration.

This gives a **validated method** (what reviewers demand) + a **novel predictive extension**
(what's interesting). The project no longer bets the paper on unmeasurable molecules.

## Dataset statistical shape (Step 1 result — constrains Steps 2–5)

`data/validation/validation_ptm_a0201.csv` (181 HLA-A\*02:01 modified-peptide rows from
Sidney/Sette 2018 MOESM3) is **mostly censored**, which bounds what later steps can claim.
Honest accounting (`pair_status` column):

- **`complete` = 50** (both native & modified IC50 numeric). **Continuous predicted-vs-measured
  IC50 regression runs on these only.** Adequate for a **pilot**; **thin for a strong
  correlation claim** — do not oversell Spearman ρ at n≈50.
- **One-sided = 41** (`native_nonbinder` 17 + `modified_nonbinder` 24): **censored at IC50
  >40,000 nM.** Treat as **binary threshold-crossing classification** (binder ↔ non-binder),
  **never** regression on the censored bound (that would fabricate an unmeasured value).
- **`both_nonbinder` = 90**: neither peptide binds → **no quantitative signal**; true-negative
  context only.
- **Highest-value subset: the 17 `native_nonbinder` rows** — the modification *creates* binding
  (native >40,000 → modified measurable). The cleanest **measured ground truth for the thesis**
  ("a non-canonical modification changes MHC-I presentation"). **Step 4's position analysis
  should target these threshold-crossing cases specifically.**
- **Exact-by-construction exclusion (the real surrogate-error test set).** A row can test surrogate
  error only if its modified sequence still contains a non-canonical code (`U/J/B/O/Z/X`); when the
  modification yields an all-canonical sequence (oxidation `F→Y`; straight deamidation `N→D`/`Q→E`;
  Phe-only pan-oxidation) the surrogate **equals** the modified peptide → error ≡ 0 by construction.
  **27 of 181 rows are exact-by-construction** (Ox-of-F 14, straight deamidation 8, Phe-only
  pan-ox 5). Net test sets: **regression on 48 of the 50 `complete` pairs**; **binary classification
  on 32 of the 41 censored one-sided rows** (80 measured rows with a genuine surrogate approximation;
  the 90 `both_nonbinder` stay excluded). Counts deterministic + self-checked.

## Assets carried over from OrthCom (one-way import; do not import the rest)

- **Physics stack (now the scientific core, not a deferred tail):** `generate_conformer.py`,
  `run_psi4_dft.py`, and the RESP/AmberTools charge→scoring approach.
- **Verified ligands:** `ligand_2NAL.pdb`, `ligand_L4F.pdb`, `ligand_PHE.pdb` — clean,
  hash-checked side-chain anchors usable as residues to graft into a pMHC peptide.
- **Institutional memory:** the "What didn't work" section and the operating discipline
  (verify-don't-trust-docs; gate before spend; deprecate-by-rename).
- **Explicitly NOT imported:** RFdiffusion(AA), LigandMPNN, the AF2-multimer instrument, the
  Lock/Key framing, the cloud package. OrthCom stays frozen in its own repo as a shelved
  sub-project.

---

## Step-by-step plan

### Step 0 — Repo creation + restructure  `DONE`
Create a clean new repo `orthcom-immuno`. Import only the carry-over assets above. Write a
fresh `README.md` and `CLAUDE.md` for the new thesis. Seed `PROGRESS_LOG.md` with the project
origin, the OrthCom→Immuno pivot rationale, and the carried-over "What didn't work" history.
**Gate:** show the new tree + both docs before any analysis. Do NOT import or reference the DL
pipeline. (This is the step the first handoff prompt covers.)

### Step 1 — Pull the validation dataset  `DONE`
Extract the modified-vs-native IC50 pairs from Sidney/Sette 2018 (PMC5863483): Table 3 + the
supplementary files (Additional files referenced as MOESM). Need, per peptide: native
sequence, modified position, modification type (citrullination / deamidation / oxidation /
chlorination), native IC50 (nM), modified IC50 (nM), restricting allele (filter to A\*02:01).
Build `data/validation_ptm_a0201.csv`. **Human must spot-check** ~5 rows against the PDF — do
not trust auto-extraction of numeric tables. **Note:** these exact IC50s were NOT pre-supplied
by the advisor; they must be read from the paper. Flag any row you're unsure of.

### Step 2 — Reproduce the surrogate baseline  `PROMPTED`
For each validation peptide, apply the canonical-surrogate method (replace the modified
residue with its nearest natural amino acid per Mattei's logic; for the PTM cases use the
sensible canonical analog, e.g. citrulline→Gln/Arg debate documented explicitly) and predict
A\*02:01 binding with NetMHCpan-4.1. Record predicted vs. measured IC50. Compute baseline
correlation (Spearman ρ, and rank-error on the "did the modification increase binding?"
call). **This is the bar to beat.** Output: `results/surrogate_baseline.csv` + a short
written summary of the surrogate's error.

**Statistical-shape constraint (see "Dataset statistical shape" above):** run the continuous
predicted-vs-measured IC50 regression on the **~50 `complete` pairs only** (report ρ as a
**pilot**, n≈50); treat the **~41 one-sided rows as binary binder/non-binder classification**
(censored >40,000 nM — never regress the censored bound); **exclude the 90 `both_nonbinder`
rows** from quantitative scoring. **Gate:** the canonical-surrogate residue mapping is reviewed
by the human *before* any NetMHCpan run (the mapping is the thing under test — a sloppy mapping
would strawman the field's method).
**Surrogate labeling (no pre-ranking — test the method as practiced):** for the two debated
residues run **both** and label honestly — citrulline `U`: *naive parent-map* = Arg, *charge-aware*
= Gln; cysteic acid `B`: *naive parent-map* = Cys, *charge-aware* = Asp. The **naive parent-map is
the field-as-practiced baseline**; report both, never only the charge-aware version (that would
flatter the surrogate). Unambiguous single maps: `J→Y`, `X→H`, `O→D`, `Z→E`.

### Step 3 — Structure-based pilot (the real go/no-go)  `TODO`
Pick ONE well-characterized PTM case. Build the modified pMHC-A\*02:01 complex from a solved
A\*02:01 crystal template (e.g. structures in PMC4317018 and the standard 2.4–2.7 Å A\*02:01
templates), graft the modified residue, run the conformer→charge→physics score adapted from
OrthCom. Compare structure-based score vs. measured IC50 vs. the surrogate prediction for that
case. **Decision:** if structure-based scoring cannot beat the surrogate on even one clean
case, the core premise is wrong → STOP and write it up as a validated negative. If it wins →
proceed. **Gate:** report the single-case result before scaling.

### Step 4 — Scale + position sweep  `TODO`
Run the structure-based method across the full PTM validation set. Then run the **position
sweep**: place the modification at anchor (P2, P9) vs. non-anchor positions and map the
presentation consequence. Compare against the 2021 prior (anchor mods tend to *abolish*
binding). Output: the surrogate-error quantification + the position-resolved liability map.
**Target the 17 `native_nonbinder` threshold-crossing cases (modification *creates* binding)
specifically — the cleanest measured ground truth for the thesis.**

### Step 5 — Predictive extension to synthetic ncAAs  `TODO`
Apply the validated method to 2-Nal and the fluoro-Phe gradient (carried-over ligands) in
A\*02:01 peptides. Produce falsifiable predictions: where does naphthyl create/destroy a
presentation liability, and does the Trp-surrogate get it wrong? Frame as wet-lab-testable.

### Step 6 — Write-up + validation collaboration  `TODO`
Target tier: *Front. Immunol.* / *PLOS Comp Biol* / *Bioinformatics* (where the comparators
live). Reach tier with a wet collaborator: *Briefings in Bioinformatics* / NAR web-server /
*Cell Reports Methods*. Use Teichmann-lab contacts to scope a small wet validation of the
Step-5 predictions — that's what elevates this from methods note to validated method.

---

## Gates (non-negotiable, inherited from OrthCom)
- No heavy spend (cloud, multi-GB downloads, GPU jobs) without an explicit human go. This
  project is CPU/local-feasible through Step 4; flag the moment that stops being true.
- Verify against actual files/outputs, not docs or agent claims.
- Only one step `PROMPTED` at a time; human marks `DONE` after verification.
- Deprecate by renaming, never delete. OrthCom repo stays frozen, not deleted.
