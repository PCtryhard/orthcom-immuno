# orthcom-immuno

**Structure-based assessment of how non-canonical amino acids perturb MHC-I antigen
presentation — benchmarked against the field's canonical-surrogate shortcut.**

## Thesis (one sentence)

The field's current way to judge whether a non-canonical amino acid (ncAA) creates an
immunological liability is to *replace it with the nearest canonical residue and run a standard
tool* (Mattei et al. 2022); this project measures **how wrong that shortcut is** by modeling the
modified residue's *actual* chemistry in the HLA-A\*02:01 groove, and produces a
**position-resolved map** of the presentation liability the surrogate misses.

## Why this is novel (and the exact boundary)

**Already published — we do NOT claim it:**
- *Mattei et al. 2022 (Front. Drug Discov.)* — the canonical-surrogate **workflow** for UAA
  immunogenicity (e.g. proxy 1-naphthyl-Ala as Trp). Explicitly framed as a "vision for future
  enhancements," i.e. unvalidated.
- *Jiang, Su, Meiler et al. 2025 (PLOS One)* — a single-allele (A\*02:01) **ML binding-affinity
  regressor** for ncAA-peptides (R²≈0.48). The "train a model" headline is taken.

**Open — the defensible contribution:** nobody has (a) scored the **actual** modified residue
structurally in the MHC-I groove, (b) **quantified the error** of the canonical-surrogate
shortcut, or (c) produced a **position-resolved liability map** for designed/bioorthogonal ncAAs.
The contribution is a **structure-based correction to, and validation of, the existing
sequence-surrogate paradigm.**

**Both outcomes publish.** Surrogate-is-wrong → "the field's shortcut mis-ranks ncAA
immunogenicity, and here's the fix" (high impact). Surrogate-is-right → "first structural
validation that the shortcut is justified, and the boundary where it breaks" (a real
methods/validation result). The **only** worthless outcome is *assuming* either — which is the
current state of the art. This is a **methods/validation project**: honest measurement of which
of those two it turns out to be is the entire point.

## Approach (staged, deliberately de-risked)

There is **no** open-literature dataset of *aromatic synthetic ncAAs* (naphthyl-Ala, fluoro-Phe)
on A\*02:01 with measured binding/T-cell data. So the project is staged so it never bets the
result on an unmeasurable molecule:

1. **Validate** the structure-based method against data that **does** exist — PTM-modified
   A\*02:01 peptides with measured IC50s (citrullination, oxidation, deamidation, chlorination on
   insulin-B; Sidney/Sette 2018, *BMC Immunol*, PMC5863483). These are chemically-defined
   non-canonical residues with quantitative ground truth on the target allele.
2. **Extend predictively** to the synthetic aromatic ncAAs (2-Nal, fluoro-Phe) where no data
   exists — framed as falsifiable predictions for a future wet-lab collaboration.

A **validated method** (what reviewers demand) plus a **novel predictive extension** (what's
interesting). The full step-by-step plan and its live status are in [`PLAN.md`](PLAN.md); the
decision log and receipts are in [`PROGRESS_LOG.md`](PROGRESS_LOG.md).

## Carried over from OrthCom (one-way import — copied, not linked)

This project pivots off **OrthCom**, a now-shelved sub-project that asked whether a *designed
protein pocket* could discriminate an ncAA from canonical residues by an affinity gap. That
capability turned out to be already published, so the pocket-design framing was retired. OrthCom
stays **frozen in its own repo**; only these assets crossed over:

- **Physics scripts** → [`src/physics/`](src/physics/): `generate_conformer.py` (SMILES → 3D
  conformer via RDKit ETKDG + MMFF) and `run_psi4_dft.py` (Psi4 B3LYP/6-31G\* single-point). The
  conformer → charge → DFT scoring path is now the **scientific core**, not a deferred tail.
  (Both are currently L4F-specific and will be generalized as the validation work requires.)
- **Verified ligands** → [`data/ligands/`](data/ligands/): `ligand_2NAL.pdb`, `ligand_L4F.pdb`,
  `ligand_PHE.pdb` — clean, hash-checked side-chain anchors (heavy-atom, CCD-ideal coordinates)
  usable as residues to graft into a pMHC peptide. Post-import hashes are recorded in
  `PROGRESS_LOG.md`.
- **Institutional memory** → the "What didn't work (and why)" history and the operating
  discipline (verify-don't-trust-docs; gate before spend; deprecate-by-rename), seeded into
  `PROGRESS_LOG.md`.

**Explicitly NOT carried over:** the deep-learning pocket-design pipeline — RFdiffusion(AA),
LigandMPNN, ProteinMPNN, localcolabfold, the AF2-multimer instrument, the cloud package, and the
Lock/Key framing. That pipeline is shelved with OrthCom and must not be reintroduced here.

## Repo map

```
orthcom-immuno/
├── README.md         # this file — the human-facing narrative
├── PLAN.md           # source of truth for the step-by-step plan (Steps 0–6) and status
├── PROGRESS_LOG.md   # authoritative decision log (what / why / effect); newest at bottom
├── CLAUDE.md         # orientation for AI agents working in this repo
├── src/
│   └── physics/      # carried-over conformer + DFT scripts (the scientific core)
├── data/
│   ├── ligands/      # the 3 verified ncAA side-chain anchors (tracked, byte-preserved)
│   └── validation/   # PTM-modified A*02:01 IC50 dataset (Step 1 output)
├── results/          # surrogate baseline + structure-based scores (Steps 2–4 output)
└── docs/             # notes, references, write-up drafts
```

## Status & discipline

Work proceeds **one step at a time** per [`PLAN.md`](PLAN.md): a step is `DONE` only after it is
**verified against actual files/outputs by a human**, never on an agent's say-so. The project is
**CPU/local-feasible through Step 4**; any heavy spend (GPU, multi-GB downloads, cloud) requires
an explicit go. Current state: **Step 0 — repo creation + restructure.**
