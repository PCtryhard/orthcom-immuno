# orthcom-immuno

Measures how wrong it is to predict a chemically modified peptide's MHC-I binding by
swapping the modified residue for its nearest canonical residue and running a standard
predictor.

## The question

That swap is how the field currently judges whether a non-canonical residue creates an
immunological liability (Mattei et al. 2022, Front. Drug Discov., where it is proposed
rather than tested). It has never been checked against measured binding. So: how large
is its error, and is it the kind of error that changes a decision?

The answer could be that the shortcut is fine. That would also be worth reporting, since
right now it is assumed rather than known. What I can say so far is that the shortcut
preserves rank order well and is quantitatively unreliable, and that in most rows it is
not really predicting anything.

Prior work that I am not claiming: the surrogate workflow itself is Mattei 2022, and a
single-allele ML affinity regressor for ncAA peptides already exists (Jiang, Su, Meiler
et al. 2025, PLOS One, R^2 about 0.48). What is open is quantifying the surrogate's
error against measured IC50s, and scoring the modified residue structurally instead.

## Method

**Data.** 181 modified and native peptide pairs on HLA-A\*02:01 with measured IC50s,
from Sidney and Sette 2018 (BMC Immunol, PMC5863483), derived from insulin B segments.
After pairing and filtering, 154 rows are scored per scheme: 61 with a numeric measured
IC50, and 93 censored at above 40,000 nM. I treat the numeric rows as binders and the
censored rows as non-binders. Modifications are citrullination (70), pan oxidation (33),
chlorination (29), oxidation of cysteine (14) and deamidation (8).

**Predictor.** NetMHCpan-4.1 standalone, both EL and BA outputs.

**Surrogate schemes.** Two, so the result is not an artefact of one mapping choice.
`naive_parent_map` replaces the modified residue with its parent canonical residue.
`charge_aware` picks the canonical residue that better matches the modified side chain's
charge state, so citrulline maps to Gln rather than back to Arg.

**Readouts.** Spearman correlation between predicted and measured values in log space on
the 61 numeric rows, and binary classification at the two standard thresholds,
EL %rank <= 2.0 and predicted affinity <= 500 nM.

## Results

Produced by Step 2. `PLAN.md` still has Step 2 as `PROMPTED`, not `DONE`, because this
repo only marks a step done after I verify it against the files by hand.

**Rank order survives the swap.** Spearman rho in log space, n = 61:
`naive_parent_map` +0.807, `charge_aware` +0.759. If all you need is to sort candidate
peptides, the shortcut is serviceable.

**Calibration does not.** At EL %rank <= 2.0, `naive_parent_map` gives accuracy 0.669
with TP 10, FN 51, FP 0, TN 93. So it recovers 10 of 61 measured binders and
misclassifies none of the 93 non-binders. `charge_aware` is within noise of this at
0.675, and 0.682 under the 500 nM affinity rule. Part of that false-negative count is a
threshold mismatch rather than surrogate error, and I say how much in Limitations.

**Most of the time the shortcut predicts nothing at all.** In 135 of 154 rows the naive
surrogate peptide is identical to the native peptide, which means the prediction asserts
that the modification changes binding not at all. On the 58 of those rows with a numeric
measured IC50:

| error of assuming no effect | rows |
|---|---|
| within 2-fold (abs log10 resid <= 0.3) | 28 |
| 2 to 10-fold | 20 |
| more than 10-fold | 10 |

Median absolute log10 residual 0.32, so about 2.1-fold. Maximum 2.04, so about 110-fold.
Assuming a post-translational modification is chemically invisible is therefore within
2-fold about half the time and off by more than tenfold in 10 of 58 cases.

`results/step2/surrogate_failure_cases.csv` holds 23 cases worth looking at directly.
The worst is a chlorinated tyrosine at position 8 of SHLVEALYLV: predicted 9.79 nM,
measured 1066 nM.

**No structure-based score exists yet.** Step 3 is the actual go/no-go and has not run,
so the repo's eventual claim, that scoring the real modified residue does better than
the swap, is at this point unmeasured.

## Limitations

The false-negative count overstates the surrogate's failure. "Binder" here means any
measurable IC50 up to 40,000 nM, while both prediction thresholds correspond to roughly
500 nM, so a sizeable number of rows counted as missed binders are genuinely weak
binders that NetMHCpan is not wrong to rank low. The assume-no-effect table is the
cleaner result because it compares a prediction against a measurement directly instead
of against a threshold.

n = 61 numeric rows is small, and the summary file flags both regression figures as a
pilot. Every peptide comes from insulin B segments on a single allele, so none of this
generalises to MHC-I as a class without more alleles and more source proteins.

The modifications measured here are post-translational, not the synthetic aromatic ncAAs
this project is ultimately aimed at. There is no open dataset of naphthylalanine or
fluoro-phenylalanine peptides on A\*02:01 with measured binding, which is why the
project validates on PTMs first. Extending to synthetic ncAAs is Step 5 and will be a
falsifiable prediction, not a result.

## How to run

Step 2 needs only the Python standard library plus NetMHCpan. NetMHCpan-4.1 standalone
has to be requested from DTU under its own licence and is not vendored here.

```bash
python data/validation/build_validation_ptm_a0201.py   # -> data/validation/validation_ptm_a0201.csv
python src/build_surrogate_input.py                    # -> results/step2/input_*.fasta, input_manifest.csv
bash   src/run_netmhcpan.sh                            # -> results/step2/pred_netmhcpan_el_ba_raw.txt
python src/parse_surrogate_baseline.py                 # -> results/step2/surrogate_vs_measured.csv
python src/summarize_surrogate_baseline.py             # -> results/surrogate_baseline_summary.txt
```

`src/seq_guard.py` rejects any peptide containing a non-canonical letter before it
reaches the predictor, which is what stops a modified sequence being sent to a tool that
cannot represent it.

The physics scripts in `src/physics/` need `rdkit` and `psi4`, both conda only:

```bash
conda install -c conda-forge rdkit psi4
```

## Relationship to OrthCom

[OrthCom](https://github.com/PCtryhard/OrthCom) is a separate, active project that
designs a protein pocket to discriminate L-2-naphthylalanine from canonical aromatics.
This repo shares two things with it: the verified ncAA side-chain anchors in
`data/ligands/` (`ligand_2NAL.pdb`, `ligand_L4F.pdb`, `ligand_PHE.pdb`, byte-preserved
with hashes recorded in `PROGRESS_LOG.md`) and the conformer and DFT scripts in
`src/physics/`.

The two ask different questions, so the pocket-design pipeline stays out of this repo by
scope: no RFdiffusionAA, no LigandMPNN, no AF2-multimer instrument here. That is a
boundary between two live projects, not a statement that either is abandoned.

## Where the reasoning is written down

`PLAN.md` is the step plan and its status. `PROGRESS_LOG.md` is the decision log,
including what failed. A step is marked done only after I check it against the files.
