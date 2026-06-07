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
