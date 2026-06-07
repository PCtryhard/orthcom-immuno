# CLAUDE.md — orthcom-immuno

Orientation for AI agents (Claude Code) working in this repo. **Read this first.**

## What this project is

A **methods/validation project** in computational immunology. One question: when a non-canonical
amino acid (ncAA) is incorporated into a peptide, does it change **MHC-I (HLA-A\*02:01) antigen
presentation** — and **how wrong is the field's shortcut** of replacing the ncAA with its nearest
canonical residue before running a standard predictor (Mattei et al. 2022)?

The deliverable is a **structure-based** score of the *actual* modified residue in the A\*02:01
groove, validated against PTM-modified peptides with measured IC50s, then extended predictively
to synthetic aromatic ncAAs. Both outcomes (shortcut-is-wrong / shortcut-is-right) are
publishable; the only failure is *assuming* the answer. See `README.md` for the full narrative.

## Sources of truth (read in this order)

1. **`PLAN.md`** — the step-by-step plan (Steps 0–6) and the **authoritative status** of each
   step. Status flows `TODO → PROMPTED → DONE`. **Only one step is `PROMPTED` at a time**, and a
   step becomes `DONE` only after a **human** has verified it against real files/outputs — never
   on an agent's say-so.
2. **`PROGRESS_LOG.md`** — the decision log: every decision as *what / why / effect*, newest at
   the bottom. **Append to it whenever you make a decision.** If this log ever disagrees with
   `PLAN.md` or `README.md`, the **log wins** — and fix the other file.
3. **`README.md`** — the human-facing narrative.

## Carried over from OrthCom (use these)

This repo pivots off **OrthCom**, a shelved pocket-design sub-project that stays frozen in its own
repo. One-way **copies** were imported — there is no submodule and no live link:

- `src/physics/generate_conformer.py`, `src/physics/run_psi4_dft.py` — the conformer → DFT
  scoring path (RDKit ETKDG/MMFF → Psi4 B3LYP/6-31G\*). The **scientific core** of this project.
  Currently L4F-specific; generalize as the work requires.
- `data/ligands/ligand_2NAL.pdb`, `ligand_L4F.pdb`, `ligand_PHE.pdb` — verified, hash-checked
  side-chain anchors (heavy-atom CCD-ideal coords). Treat as **read-only reference**; their
  post-import hashes (`git hash-object` + SHA-256) are recorded in `PROGRESS_LOG.md`. Re-verify
  before use.
- The **"What didn't work (and why)"** institutional memory and the operating discipline (below).

## Do NOT import or reintroduce — the shelved DL pipeline

The OrthCom **deep-learning pocket-design pipeline is shelved and must not be wired back in:**
- RFdiffusion / RFdiffusionAA, LigandMPNN, ProteinMPNN, localcolabfold.
- The AF2-multimer instrument (`prep_multimer.py`, `run_multimer.sh`, `parse_interface.py`).
- The cloud (Apptainer / GPU) package.
- The **Lock / Key** framing, and any "pocket discrimination" or "platform / skyscraper" framing.

If you find yourself reaching for any of these, **stop** — that is the old project's gravity, and
escaping it is exactly why this repo exists. The pivot rationale is recorded in `PROGRESS_LOG.md`.

## Gates (non-negotiable)

- **No heavy spend without an explicit human go.** Cloud runs, GPU jobs, multi-GB downloads:
  show the exact command and get a "yes" first. This project is CPU/local-feasible through
  Step 4 — flag the moment that stops being true.
- **Verify, don't trust.** Check claims against actual files/outputs, not docs or an agent's
  summary. Docs and agents have been wrong here before — that discipline is inherited, not
  hypothetical.
- **One step at a time.** Do only the step that is `PROMPTED` in `PLAN.md`; stop at its gate and
  report for human verification before touching the next.
- **Deprecate by renaming, never deleting.** Preserve history; OrthCom stays frozen, not deleted.
- **Small, reviewable commits; show diffs; never force-push.**

## Environment

CPU/local through Step 4. The physics stack needs **RDKit + Psi4** (and, for partial charges, the
**RESP / AmberTools** approach); these live in a conda environment that is **not** on the default
PATH — activate it explicitly or call tools by absolute path. Record environment specifics in
`PROGRESS_LOG.md` as you establish them.
