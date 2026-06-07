# Sources — validation dataset (Step 1)

## Primary source
- **Paper:** Sidney J, Vela JL, Friedrich D, Kolla R, von Herrath M, Wesley JD, Sette A.
  *"Low HLA binding of diabetes-associated CD8+ T-cell epitopes is increased by post
  translational modifications."* **BMC Immunology** 2018.
- **DOI:** 10.1186/s12865-018-0250-3 · **PMC:** PMC5863483 · **License:** CC-BY 4.0 (open access).
- **Accessed:** 2026-06-07.

## File archived here
- `12865_2018_250_MOESM3_ESM.docx` — *Additional file 3: "Measured binding affinities."* The
  complete table of modified and native insulin-derived peptides and their measured
  **HLA-A\*02:01** binding (IC50, nM).
  - Downloaded from Springer static content:
    `https://static-content.springer.com/esm/art%3A10.1186%2Fs12865-018-0250-3/MediaObjects/12865_2018_250_MOESM3_ESM.docx`
    (the PMC `/articles/instance/5863483/bin/…` URL returned a non-docx stub).
  - **sha256:** `ffbebafed1737fc06cab1124e1c0043fa09780fb95c8e68f60d3b47bd3a0a27f`

## How the CSV was built (deterministic; no fabrication)
```
python ../build_validation_ptm_a0201.py 12865_2018_250_MOESM3_ESM.docx ../validation_ptm_a0201.csv
```
Reads the Word table XML directly (Python stdlib), maps columns, and derives the modified
peptide position by diffing native vs modified sequence. No model, OCR, or hand-transcription.

## Source conventions (read before using the CSV)
- **`-` in an IC50 column = non-binder, IC50 > 40,000 nM** (paper Methods). In the CSV the numeric
  field is left **blank** (never a guessed number); `pair_status` marks it `*_nonbinder`; the raw
  `-` is preserved in `*_IC50_raw`. *(Exact threshold value to confirm on human spot-check.)*
- Modified-residue letter codes (MOESM3 legend): `U` citrulline (←Arg), `J` chlorotyrosine
  (←Tyr), `O` iso-aspartate / `Z` iso-glutamate (deamidation of Asn/Gln), `B` cysteic acid
  (←Cys; "Oxidation of C"). **`X`** (oxidized His, used in "Pan oxidation" rows) appears in
  sequences but is **not** defined in the MOESM3 legend — inferred from the Methods.
- Assay: classical radiolabeled-ligand competition binding to purified HLA-A\*02:01; IC50 in nM.
