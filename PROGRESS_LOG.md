# PROGRESS_LOG.md — orthcom-immuno

Authoritative decision log. Every decision as **what / why / effect**; **newest at the bottom.**
If this log disagrees with `PLAN.md` or `README.md`, the log wins — fix the other. `PLAN.md` holds
the step plan and status; this file holds the reasoning and the receipts.

---

## Project origin & the OrthCom → Immuno pivot (2026-06-07)

**Origin.** This project descends from **OrthCom**, which asked whether a computationally designed
protein pocket could **discriminate a non-canonical amino acid (ncAA) from canonical residues by
an affinity gap.** OrthCom built a ligand-aware design pipeline (RFdiffusionAA → LigandMPNN →
AF2-multimer co-fold → physics) around verified ncAA anchors (2-Nal, 4-fluoro-Phe, Phe).

**Why pivot.** The core capability OrthCom was chasing — designing pockets that discriminate
ncAAs — turned out to be **already published** (Baker/DeGrado-lab work, 2024). The pocket-design
framing is therefore **shelved, not deleted**: OrthCom stays frozen in its own repo (its HEAD
`ed46c3d` records this same Strategic Hold on the OrthCom side).

**The new question — open, and adjacent to the human's prior skill.** Instead of *designing
pockets for* ncAAs, ask whether incorporating an ncAA into a peptide creates an **immunological
liability** — does it change **MHC-I presentation**? The field's current method (Mattei et al.
2022) replaces the ncAA with its nearest canonical residue and runs a standard, NetMHCpan-adjacent
tool; **nobody has validated whether that shortcut is correct.** We will, using **structure-based**
scoring against **PTM-modified HLA-A\*02:01 peptides with measured IC50s** (Sidney/Sette 2018,
*BMC Immunol*, PMC5863483) as ground truth, then extend predictively to synthetic aromatic ncAAs.
Both outcomes publish; the only failure is assuming the answer. Full plan in `PLAN.md`.

**What carried over (one-way copy):** the physics stack (`generate_conformer.py`,
`run_psi4_dft.py`, and the RESP/AmberTools charge → scoring approach), the three verified ligand
anchors, and the institutional memory below. **What did not:** the entire DL pocket-design
pipeline and the Lock/Key framing (see `CLAUDE.md`).

---

## Inherited from OrthCom — "What didn't work (and why)"

> Copied **verbatim** from OrthCom `PROGRESS_LOG.md` (HEAD `ed46c3d`). It refers throughout to the
> now-shelved pocket-design pipeline (RFdiffusionAA, LigandMPNN, partial diffusion, the Lock/Key
> framing); it is preserved as durable institutional memory, **not** as current strategy. The
> lesson that carries directly into *this* project is item 1 — a single fluorine is not a
> discriminable feature for sequence design, so discrimination must come from **shape and/or
> physics** — together with the operating discipline (verify-don't-trust; gate before spend;
> deprecate-by-rename).

Institutional memory — keep it. Newest first. (Detailed evidence is in the dated decision-log
entries above; this is the durable index, mirrored in `README.md`.)

1. **A single fluorine is not a discriminable feature (for sequence design).** LigandMPNN with
   ligand atom-context ON moved its sequence choice between L4F and PHE by Jensen–Shannon
   divergence **mean ~0.0034 bits** F-proximal (max 0.0099) vs up to **~0.60 bits** for
   whole-ligand context — ~50–100× smaller. → MPNN *sequence* design cannot discriminate one
   para-F; discrimination must come from pocket **shape** (RFdiffusionAA) and/or **physics**.
   Hence **L4F → negative control** and the primary target → **2-Nal**.
2. **Helical Key in a barrel = frustrated interface (era 2, partial diffusion).** Shallow,
   one-sided dock (44/82 Key atoms buried, COM ~13 Å off centre); AF2 co-folds buried it more
   but **deformed the helix** (helical windows 12→2–8), introduced hard clashes (to 0.5 Å),
   Key pLDDT stuck ~30 while the Lock folded fine. Helix-vs-barrel-loop mismatch; `partial_T`
   cannot change loop length or relieve a hard clash. → abandoned partial diffusion.
3. **The original anchor was chemically malformed.** `graft_residue.py` auto-named atoms by
   element count (`O1/O2/C3…C9`), kept a spurious 2nd carboxylate O and explicit H — wrong for
   a mid-peptide residue, unusable as a clean ligand. → rebuilt by `prep_ligand.py` from
   verified CCD coords (side-chain-centric, 0 free O).
4. **De-novo hallucination = "pool-noodle" helices (era 1).** Hallucinating 100–200 residues
   around the fixed Key gave long α-helices, not pockets (~1000 designs to fish for a
   concavity). → motivated the scaffold pivot (era 2), then the ligand-aware approach (current).

*Smaller gotchas:* DL tools crash on `L4F` (`KeyError`) → PHE "Trojan-Horse" stand-in for DL
steps; Tyr rejected (hydroxyl polar penalty repels the swap-back F); AF2's BFC allocator
segfaults at recycle 2 → JAX `platform` allocator; synthetic proteins have no MSA →
`--msa-mode single_sequence`.

---

## Recorded asset hashes (post-import; verify before use)

The three ligand anchors are byte-identical to OrthCom's verified copies. Confirmed three ways:
(a) `git hash-object` matches OrthCom's recorded values; (b) byte-identical regeneration from the
committed builders (`prep_ligand.py`, `prep_ligand_2nal.py`); (c) the builders' own geometry and
SHA-1 self-checks pass. The `git hash-object` values `ee44ce78` (L4F) and `fc4bdc78` (PHE) match
the values OrthCom's handoff recorded; 2-Nal's builder pins raw-SHA-1 `ac9f19f7…`.

> **The two 2-Nal hashes are different hash types over different objects — not a contradiction.**
> `f1f6ede6…` (the value in the table) is the **`git hash-object` of the file blob** — git-SHA-1
> over the PDB's bytes including git's blob header — whereas `ac9f19f7…` is the **builder's
> internal raw-SHA-1 self-check over its generated coordinate content** — plain SHA-1, no git
> header. Both are expected, both pass, and neither is in conflict with the other.

| file (in `data/ligands/`) | `git hash-object` (git blob SHA-1)         | SHA-256                                                            |
| ------------------------- | ------------------------------------------ | ----------------------------------------------------------------- |
| `ligand_2NAL.pdb`         | `f1f6ede6c33d3edb0b2aa1116d6989aa8ef93ee9` | `25aa6c87be1d836e04e8bcb0740df542cf85915ee6f37892bc59c09998e83e8c` |
| `ligand_L4F.pdb`          | `ee44ce78f53399e70b56945c7656192d061eaf2f` | `3782376b30706ad48ce339023e972876146b2c4d230c89e8fc47c059cfc4636a` |
| `ligand_PHE.pdb`          | `fc4bdc78cbfae8da665b9d6070ffc831ac99a358` | `8f9a94a7d36a63c9aa8e84cbdd81529c526d45291f31476accb6a017a14993c3` |

---

## 2026-06-07 — Step 0: repo creation + carry-over import

- **What.** Created `orthcom-immuno` (fresh `git init`, new history, branch `main`) with layout
  `src/physics/  data/ligands/  data/validation/  results/  docs/`. Imported, one-way copy:
  - `src/physics/generate_conformer.py`   ← OrthCom `src/generate_conformer.py`
  - `src/physics/run_psi4_dft.py`         ← OrthCom `src/run_psi4_dft.py`
  - `data/ligands/ligand_2NAL.pdb`, `ligand_L4F.pdb`, `ligand_PHE.pdb` ← OrthCom `data/raw/`
  Wrote fresh `README.md`, `CLAUDE.md`, and this log; copied `PLAN.md` in unchanged (byte-identical
  to the supplied source).
- **Why.** Establish a clean home for the immunogenicity thesis, free of the shelved DL pipeline's
  gravity, while preserving the physics core, the verified ligands, and the institutional memory
  (handoff Step 0).
- **Effect / receipts.**
  - OrthCom was confirmed **clean at HEAD `ed46c3d`** before any copy — no copy from a dirty tree.
    (At session start its Task B was uncommitted; it was committed externally to `ed46c3d` during
    the session — by the user / a separate Claude session — not by this agent.)
  - The earlier ligand "hash mismatch" was resolved: OrthCom's recorded `ee44ce78` / `fc4bdc78`
    are `git hash-object` values, not SHA-256; under the correct function the on-disk files match
    exactly. The two physics scripts import only RDKit/Psi4 + stdlib — no internal helper needed,
    so the two files are the complete physics carry-over.
  - Ligand hashes recorded above; byte-identity to OrthCom verified by `cmp`.
  - **NOT imported:** RFdiffusion(AA), LigandMPNN, ProteinMPNN, localcolabfold, the AF2-multimer
    instrument (`prep_multimer.py` / `run_multimer.sh` / `parse_interface.py`), the cloud package,
    the Lock/Key framing — shelved with OrthCom.
- **Gate.** Step 0 stops here for human verification; flipping `PLAN.md` Step 0 → `DONE` is the
  human's call. No Step 1 work (no NetMHCpan, no Psi4 runs, no dataset extraction) was performed.
- **Known follow-up.** GitHub repo creation + push is pending: `gh` is not installed on this
  machine. The local repo + initial commit are complete; remote creation/push awaits a decision
  (install `gh`, or create the remote another way).

---

## 2026-06-07 — Step 0 closed (pushed)
- Scaffold commit `8c382f7` + doc-clarification commit `ffd793e` (the 2-Nal hash-types note)
  pushed to **private** GitHub `PCtryhard/orthcom-immuno` (`main`, via authed `gh`; no auth
  prompt). User confirmed; `PLAN.md` Step 0 → `DONE`.

## 2026-06-07 — Step 1: PTM validation dataset (build only; awaiting spot-check)
- **What.** Built `data/validation/validation_ptm_a0201.csv` — the modified-vs-native peptide
  IC50 dataset for **HLA-A\*02:01** — from Sidney & Sette et al. 2018 (*BMC Immunol*,
  PMC5863483, DOI 10.1186/s12865-018-0250-3), **Additional file 3 (MOESM3)** "Measured binding
  affinities." 181 modified-peptide rows (the complete A\*02:01 PTM set); 16 columns incl.
  per-row provenance (MOESM3 table-row index). CSV sha256 `54407a90…b557a95107`.
- **Source + method (reproducible; no fabrication).**
  - MOESM3 `.docx` archived at `data/validation/sources/` (sha256 `ffbebafe…a0a27f`, CC-BY) —
    downloaded from Springer static-content (PMC `/bin/` URL returned a stub). See
    `sources/SOURCES.md`.
  - Extracted **deterministically** from the Word table XML by
    `data/validation/build_validation_ptm_a0201.py` (stdlib only) — no model, OCR, or
    hand-transcription. IC50 columns mapped: "Modified" → modified-peptide IC50, "WT" → native
    IC50; orientation **verified** by reproducing the paper's own "4 of 9 reported epitopes
    improved ≥2-fold" (ALWMRLLPL, ALWMRLLPLL, SHLVEALYLV, HLVEALYLV).
  - Modification **position** derived by diffing native vs modified sequence (1-based peptide
    coords; ";"-joined for multi-site oxidation).
- **`-` handling (critical).** Per Methods, a dash = **IC50 > 40,000 nM (non-binder)**, not
  "missing." Recorded as a **blank** numeric field (never guessed/interpolated); `pair_status`
  = `*_nonbinder`; the raw `-` is preserved in `*_IC50_raw`.
- **Counts.** `pair_status`: **complete 50** (both numeric), modified_nonbinder 24,
  native_nonbinder 17, both_nonbinder 90 (= 181). By modification: Citrullination 70,
  Pan oxidation 38, Chlorination 29, Deamidation 16, Oxidation of C 14, Oxidation of F 14.
- **Flags / uncertainties (for spot-check).**
  - `X` (oxidized His in "Pan oxidation" rows) appears in modified sequences but is **not** in
    the MOESM3 letter legend (defines only B/J/O/Z/U); meaning inferred from Methods.
  - The `-` = ">40,000 nM" threshold was read from the Methods via web fetch — **confirm vs the PDF**.
  - No value was guessed: MOESM3 is digital text (not a scanned figure), so every numeric IC50 is
    read exactly; no figure-derived values were used.
- **Scope.** Dataset assembly only — **no NetMHCpan, no Psi4, no analysis** (Steps 2–3).
- **Gate.** STOP for human spot-check (≥5 CSV rows vs the PDF/MOESM3) before `PLAN.md` Step 1 →
  `DONE`. *(Resolved: spot-check passed; committed `154b5c1` + PLAN `3ef0247`, pushed; Step 1 `DONE`.)*

## 2026-06-07 — Step 2 prep: surrogate mapping approved + test-set sizing (no run yet)
- **Approved canonical-surrogate mapping** (the method under test; reviewed before any run):
  - Citrulline `U`: **naive parent-map → Arg**, **charge-aware → Gln** — run **both**, label honestly,
    **do NOT pre-rank** (the field practices the naive parent-map; Gln is the biochemical upgrade).
  - Cysteic acid `B`: **naive parent-map → Cys**, **charge-aware → Asp** (same labeling logic).
  - Unambiguous single maps: chlorotyrosine `J`→Tyr, oxo-His `X`→His, iso-Asp `O`→Asp, iso-Glu `Z`→Glu.
  - Oxidized-Met → Met (parent) — **not instantiated** (no Met in any modified peptide here).
- **Exact-by-construction exclusion (deterministic, self-checked).** A row tests surrogate error only
  if its modified sequence still holds a non-canonical code (`U/J/B/O/Z/X`); an all-canonical modified
  sequence ⇒ surrogate ≡ modified ⇒ error ≡ 0. **27 of 181 exact** (Ox-of-F 14, straight deamidation
  8, Phe-only pan-ox 5). Real test sets: **regression on 48/50 complete pairs**; **binary
  classification on 32/41 censored one-sided rows** (80 measured rows total carry a real surrogate
  approximation; 90 `both_nonbinder` excluded). Verified: 27 exact + 154 testable = 181; per-status
  sums reconcile (complete 2+48=50, native_nb 4+13=17, modified_nb 5+19=24, both_nb 16+74=90).
- **NetMHCpan-4.1 NOT installed** — absent from git-bash PATH, WSL PATH, and common dirs. **Gate:
  STOP** — academic license/download flagged to the human (DTU Health Tech; Linux build for WSL). No
  NetMHCpan run, no Psi4 performed.
- **Effect.** PLAN "Dataset statistical shape" gains the exact-by-construction test-set sizes; Step 2
  records the naive-parent-map vs charge-aware labeling. Awaiting NetMHCpan + human review of the run
  command before any execution.

## 2026-06-07 — Step 2: predictor route (IEDB standalone) + hard input guard
- **DTU NetMHCpan-4.1 blocked** (academic form inaccessible). **Pivot to the IEDB MHC-I standalone**,
  which bundles **`netmhcpan_el` = NetMHCpan-4.1 EL** as a local method (no DTU form, runs offline).
  Route **confirmed viable** (research + WSL probe):
  - Package: **`IEDB_MHC_I-3.1.7.tar.gz`** (latest; **~360 MB**), free academic license (Academic
    Software License + Non-Profit OSL v3.0). URL:
    `https://downloads.iedb.org/tools/mhci/3.1.7/IEDB_MHC_I-3.1.7.tar.gz`.
  - Requirements: Python 3.6+ (**WSL has 3.12.3 ✓**), **tcsh** (MISSING), **gawk** (✓), Linux 64-bit
    (**WSL Ubuntu 24.04 ✓**). Install: unpack → run `./configure`; path containing `/mhc_i` must be
    < 57 chars.
  - Run syntax (from README): `./src/predict_binding.py netmhcpan_el HLA-A*02:01 <len> <input.fasta>`;
    input = FASTA.
  - **Planned install path:** `/home/max/mhc_i` (17 chars, < 57 ✓).
- **Blocker (human action — `sudo` needs a password here, I can't apt-install):**
  `sudo apt install -y tcsh make unzip`. **tcsh is the hard requirement** (the netMHCpan wrapper is a
  tcsh script). The ~360 MB download also awaits an explicit go.
- **Hard input guard added — `src/seq_guard.py`** (committed): HALTS if any peptide sent to the
  predictor contains a non-canonical code (`U/J/B/O/Z/X`). NetMHCpan silently maps unknown letters to
  `X`, so a leftover code would yield confident-looking garbage; the guard makes that loud and stops
  the run. **Self-test passes** (canonical pass; all six codes halt with position+meaning; `guard_all`
  halts the batch on first offender). Surrogate substitution must run upstream; `guard_all()` is to be
  called immediately before every predictor invocation.
- **No predictions, no Psi4.** Awaiting: (a) human installs tcsh; (b) approval of the ~360 MB download;
  then I install + configure and show the exact run command + a 3-row formatted-input sample for
  review before any prediction.

## 2026-06-07 — Step 2 RUN: surrogate predictions generated (NetMHCpan-4.1 EL+BA)
- **Install completed.** IEDB MHC-I standalone unpacked to **`/home/max/mhc_i`** (WSL-native, NOT
  vendored in git). `./configure` first failed (Ubuntu 24.04 has `python3` but no `python`, and the
  configure shebang is `#!/usr/bin/env python`) → ran `python3 src/configure.py` = "successfully
  installed". `netmhcpan-4.1-executable` + EL & BA percentile-data present.
- **Input built** by `src/build_surrogate_input.py`: 154 testable rows → per-row `native` /
  `naive_surrogate` (U→Arg, B→Cys) / `charge_surrogate` (U→Gln, B→Asp), each passed through
  `src/seq_guard.py` (HALT on any leftover U/J/B/O/Z/X; all 462 passed). Outputs in `results/step2/`:
  `input_{9,10}mers.fasta`, **`input_peptides.txt` (209 unique peptides)**, `input_manifest.csv`
  (peptide → row/variant/measured-IC50).
- **CORE FRAMING (as requested).** For **135 of the 154 testable rows the naive parent-map surrogate
  reconstructs the *native* peptide** (citrulline→Arg, chlorotyrosine→Tyr, cysteic acid→Cys,
  oxo-His→His all revert to the native residue). For those rows the field's naive surrogate is
  **mathematically equivalent to predicting "the modification has no effect"** — the predictor is
  handed the native sequence and is structurally blind to the PTM. **That is the core thing Step 2
  measures.** The **19 exceptions are partial-modeling cases**: deamidation (O/Z → the canonical
  product Asp/Glu, ≠ native Asn/Gln) and pan-oxidation peptides whose Phe→Tyr (`F→Y` is canonical and
  kept, so the surrogate models that part but not the C→cysteic-acid / His→oxo-His parts).
- **Predictor HANG — debugged.** The IEDB **`predict_binding.py` Python wrapper hangs (>1 h)** and —
  confirmed via `ps` — **never spawns the netMHCpan binary** (Python blocked, load 0.00). Ruled out:
  drvfs (`/mnt/c` read = 12 ms; hang reproduces fully WSL-native), peptide count (3 peptides also
  hang), CPU contention / stale lock (system idle, no leftover procs; netMHCpan uses per-run
  `mkdtemp`). The **netMHCpan-4.1 binary itself is fine** — invoked directly in `-p` (peptide-list)
  mode it did **all 209 peptides in ~1 s**.
- **FIX = call the binary directly** via `src/run_netmhcpan.sh`:
  `netMHCpan -p -BA -a HLA-A02:01 -f input_peptides.txt`. Identical NetMHCpan-4.1 EL+BA model; one
  model load; mixed 9/10-mers in one list. `-BA` puts **both heads in one file** →
  `results/step2/pred_netmhcpan_el_ba_raw.txt` (209 rows; 4 strong + 11 weak binders). Columns:
  `Score_EL Rnk_EL` (presentation head) and `Score_BA Rnk_BA Aff(nM)` — **`Aff(nM)` = predicted IC50**
  for the predicted-vs-measured regression. (Sanity: citrulline naive surrogate `ALWMRLLPL` (=native)
  → 12.79 nM; charge-aware `ALWMQLLPL` → 11.44 nM.)
- **Note vs the originally-named wrapper.** Same NetMHCpan-4.1 model; the binary reports its own native
  %rank (≈ the IEDB percentile-data %rank). The wrapper is simply unusable (hangs).
- **Gotcha logged.** Invoking the WSL command from git-bash needs `MSYS_NO_PATHCONV=1`, else git-bash
  rewrites `/home/...` arguments into `C:\Program Files\Git\...`.
- **STOP — no parsing/analysis yet, no Psi4.** Next gate: join predictions to `input_manifest.csv`
  (predicted `Aff(nM)` / EL-rank vs measured IC50) for the **48 regression + 32 classification** rows,
  reported under the naive-parent-map vs charge-aware labels.

## 2026-06-07 — Step 2 PARSE: predicted-vs-measured joined table (no stats yet)
- **What.** `src/parse_surrogate_baseline.py` joins `pred_netmhcpan_el_ba_raw.txt` → `input_manifest.csv`
  → measured modified IC50, emitting **`results/step2/surrogate_vs_measured.csv` (308 rows)**: one row
  per (testable source row × surrogate scheme {naive_parent_map, charge_aware}) with the scheme's
  surrogate peptide, its predicted `Aff(nM)`/EL-score/EL-rank, the predicted native `Aff(nM)` (context),
  the measured modified IC50, and a per-row surrogate error
  `log10_resid = log10(pred_Aff) − log10(measured)`.
- **Why.** The per-row residual (in log space) IS the surrogate error; split by scheme so the field's
  *naive parent-map* baseline and the *charge-aware* refinement are reported separately.
- **Join integrity (both failure modes guarded).**
  - Keyed on the **exact peptide sequence** (netMHCpan `Peptide` column, not `Core`), via a dict — so
    9-mers and 10-mers cannot collide (key lengths present = {9,10}).
  - **154 source rows IN → 308 comparison rows OUT** (×2 schemes); every surrogate peptide found in the
    209 predictions (0 missing); all 209 predictions used.
  - **Many-to-one is by design and does NOT collapse rows:** 43 surrogate sequences are shared by >1
    source row (e.g. `FVNQHLCGS` ← r14, r18 — two different modifications of the same native), each
    attached to the same prediction while keeping both source rows distinct. No cross-join.
- **Censoring (special handling).** 186 of 308 rows have a non-binder modified peptide (IC50 > 40,000
  nM) → `measured_modified` carries the literal flag **`>40000`**, never coerced to a number;
  `log10_resid` left **blank**. 122 rows have a numeric measured modified IC50 (complete 48×2 +
  native_nonbinder 13×2) → residual computed.
- **Scheme divergence.** naive vs charge differ only for the **107** rows containing U (citrulline) or
  B (cysteic acid); they **coincide for 47** rows (deamidation O/Z, chlorination J, oxo-His X — same
  map in both schemes), where the two scheme rows are identical apart from the label.
- **No summary statistic computed** (no Spearman, no classification) — table only. **STOP for human
  review of the joined table.** No Psi4.
