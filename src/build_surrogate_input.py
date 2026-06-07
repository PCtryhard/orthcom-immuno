#!/usr/bin/env python3
"""src/build_surrogate_input.py — Step 2: build the predictor input from the surrogate mapping.

For each TESTABLE row of data/validation/validation_ptm_a0201.csv (modified sequence still
contains a non-canonical code U/J/B/O/Z/X) it emits three canonical variants:
  - native           : the unmodified peptide
  - naive_surrogate  : naive parent-map   U->Arg, B->Cys ; J->Tyr, X->His, O->Asp, Z->Glu
  - charge_surrogate : charge-aware        U->Gln, B->Asp ; J->Tyr, X->His, O->Asp, Z->Glu
(The two schemes differ ONLY at U and B — the residues the field genuinely debates.)

Every emitted sequence is passed through seq_guard.assert_canonical (HALT on any leftover
non-canonical letter) BEFORE it can reach the predictor. Writes one FASTA per peptide length
plus a manifest linking each peptide back to its row/variant/measured IC50. Runs NO predictor.
"""
import csv, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from seq_guard import assert_canonical

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CSV  = os.path.join(REPO, "data", "validation", "validation_ptm_a0201.csv")
OUT  = os.path.join(REPO, "results", "step2")

NONCANON = set("UJBOZX")
NAIVE  = {"U": "R", "B": "C", "J": "Y", "X": "H", "O": "D", "Z": "E"}  # naive parent-map
CHARGE = {"U": "Q", "B": "D", "J": "Y", "X": "H", "O": "D", "Z": "E"}  # charge-aware

def sub(scheme, seq): return "".join(scheme.get(c, c) for c in seq)
def testable(seq):    return any(c in NONCANON for c in seq)

def main():
    os.makedirs(OUT, exist_ok=True)
    rows = list(csv.DictReader(open(CSV)))
    trows = [r for r in rows if testable(r["modified_sequence"])]

    manifest, naive_eq_native = [], 0
    for r in trows:
        mod, nat = r["modified_sequence"], r["native_sequence"]
        nsur, csur = sub(NAIVE, mod), sub(CHARGE, mod)
        for s, v in ((nat, "native"), (nsur, "naive_surrogate"), (csur, "charge_surrogate")):
            assert_canonical(s, label=f"{v}:{r['provenance']}")   # HALT if anything slipped through
        naive_eq_native += (nsur == nat)
        for pep, v in ((nat, "native"), (nsur, "naive_surrogate"), (csur, "charge_surrogate")):
            manifest.append({
                "peptide": pep, "length": len(pep), "variant": v,
                "modification_type": r["modification_type"], "pair_status": r["pair_status"],
                "modified_sequence": mod,
                "measured_modified_IC50_nM": r["modified_IC50_nM"],
                "measured_native_IC50_nM": r["native_IC50_nM"],
                "source_row": r["provenance"].split("tbl-row ")[-1],
            })

    for L in (9, 10):
        peps = sorted({m["peptide"] for m in manifest if m["length"] == L})
        with open(os.path.join(OUT, f"input_{L}mers.fasta"), "w", newline="\n") as fh:
            for i, p in enumerate(peps, 1):
                fh.write(f">pep{L}_{i}\n{p}\n")
        print(f"length {L}: {len(peps)} unique peptides -> results/step2/input_{L}mers.fasta")

    # Unified peptide list for the netMHCpan -p (peptide-list) run: ALL unique peptides, both
    # lengths, one per line. This is what the predictor actually consumes (see src/run_netmhcpan.sh);
    # the per-length FASTA files are kept for reference. LF newlines (netMHCpan chokes on CR).
    allpeps = sorted({m["peptide"] for m in manifest})
    with open(os.path.join(OUT, "input_peptides.txt"), "w", newline="\n") as fh:
        fh.write("\n".join(allpeps) + "\n")
    print(f"unified peptide list: {len(allpeps)} peptides -> results/step2/input_peptides.txt")

    cols = ["peptide", "length", "variant", "modification_type", "pair_status",
            "modified_sequence", "measured_modified_IC50_nM", "measured_native_IC50_nM", "source_row"]
    with open(os.path.join(OUT, "input_manifest.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n"); w.writeheader(); w.writerows(manifest)

    uniq = len({(m["peptide"], m["length"]) for m in manifest})
    print(f"testable rows: {len(trows)} | manifest entries: {len(manifest)} (rows x3 variants) | "
          f"unique peptides: {uniq}")
    print(f"naive_surrogate == native: {naive_eq_native}/{len(trows)} "
          f"(naive != native only for deamidation, where both schemes use the canonical product Asp/Glu)")

    print("\n--- 3-row sample (one per modification family) ---")
    for want in ("Citrullination", "Oxidation of C", "Deamidation"):
        r = next((x for x in trows if x["modification_type"] == want), None)
        if not r: continue
        mod = r["modified_sequence"]
        meas = r["modified_IC50_nM"] or "-(>40k nonbinder)"
        print(f"[{want}] row {r['source_row'] if 'source_row' in r else r['provenance'].split('tbl-row ')[-1]} "
              f"({r['pair_status']}) modified={mod}  measured_modified_IC50={meas} nM")
        print(f"    native           : {r['native_sequence']}")
        print(f"    naive_surrogate  : {sub(NAIVE, mod)}   (U->Arg, B->Cys)")
        print(f"    charge_surrogate : {sub(CHARGE, mod)}   (U->Gln, B->Asp)")

if __name__ == "__main__":
    main()
