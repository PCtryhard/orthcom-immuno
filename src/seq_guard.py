#!/usr/bin/env python3
"""src/seq_guard.py — hard guard for predictor input (Step 2+).

WHY THIS EXISTS
---------------
NetMHCpan (and the IEDB MHC-I wrapper around it) silently maps any UNKNOWN residue
letter to 'X' (its wildcard). So a leftover non-canonical code from our PTM dataset
(U/J/B/O/Z/X) does NOT raise — it yields a confident-looking but meaningless score.
That is exactly the failure mode that would silently corrupt the surrogate baseline.

This guard makes that failure LOUD: any sequence containing a letter outside the 20
canonical amino acids HALTS the run. The surrogate substitution MUST happen upstream
(U->Arg/Gln, J->Tyr, B->Cys/Asp, O->Asp, Z->Glu, X->His); nothing non-canonical may
reach the predictor. Import `assert_canonical` / `guard_all` and call them on EVERY
peptide immediately before the predictor invocation.
"""
from __future__ import annotations

# The 20 standard amino acids. Deliberately EXCLUDES B, J, O, U, X, Z so every
# non-canonical code in data/validation/ is rejected (X included: in our data it is
# oxidized-His and must be substituted, and it is also NetMHCpan's silent wildcard).
CANONICAL_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")

# For a helpful halt message (maps each non-canonical code to its meaning + surrogate).
KNOWN_NONCANONICAL = {
    "U": "citrulline (from Arg; surrogate Arg/Gln)",
    "J": "chlorotyrosine (from Tyr; surrogate Tyr)",
    "B": "cysteic acid (from Cys; surrogate Cys/Asp)",
    "O": "iso-aspartate (from Asn; surrogate Asp)",
    "Z": "iso-glutamate (from Gln; surrogate Glu)",
    "X": "oxidized-His / NetMHCpan wildcard (surrogate His)",
}


class NonCanonicalSequenceError(ValueError):
    """Raised when a sequence bound for the predictor contains a non-canonical residue."""


def assert_canonical(seq: str, *, label: str | None = None) -> str:
    """Return `seq` unchanged iff every residue is canonical; otherwise HALT (raise)."""
    offenders = [(i + 1, c) for i, c in enumerate(seq) if c not in CANONICAL_AA]
    if offenders:
        where = ", ".join(
            f"pos{i}:{c}" + (f" [{KNOWN_NONCANONICAL[c]}]" if c in KNOWN_NONCANONICAL else "")
            for i, c in offenders
        )
        tag = f" [{label}]" if label else ""
        raise NonCanonicalSequenceError(
            f"non-canonical residue(s) would be silently mangled to 'X' by the predictor"
            f"{tag}: {seq!r} -> {where}. Surrogate substitution must happen UPSTREAM. HALTING."
        )
    return seq


def guard_all(peptides, *, labels=None):
    """Validate every peptide; raise on the FIRST offender so no partial/garbage run starts."""
    peptides = list(peptides)
    labels = list(labels) if labels is not None else [None] * len(peptides)
    for pep, lab in zip(peptides, labels):
        assert_canonical(pep, label=lab)
    return peptides


def _selftest() -> None:
    # 1) clean canonical peptides pass through untouched
    for ok in ("ALWMRLLPL", "SHLVEALYL", "MALWMRLLP", "HLVEALYLV"):
        assert assert_canonical(ok) == ok, ok
    # 2) every non-canonical code in our dataset must HALT
    bad_cases = {
        "U": "ALWMULLPL", "J": "SHLVEALJL", "B": "HLBGSHLVEA",
        "OZ": "FVOZHLCGS", "X": "XLBGSXLVEA",
    }
    for code, bad in bad_cases.items():
        try:
            assert_canonical(bad, label=f"test:{code}")
        except NonCanonicalSequenceError as e:
            print(f"  HALT (correct) [{code}] {bad}: {e}")
        else:
            raise SystemExit(f"GUARD FAILED: {bad} passed but must halt")
    # 3) guard_all halts on the first offender within a batch
    try:
        guard_all(["ALWMRLLPL", "ALWMULLPL", "SHLVEALYL"])
    except NonCanonicalSequenceError:
        print("  guard_all halts the batch on first offender (correct)")
    else:
        raise SystemExit("GUARD FAILED: guard_all did not halt the batch")
    print("seq_guard self-test PASS: canonical pass; all of U/J/B/O/Z/X halt.")


if __name__ == "__main__":
    _selftest()
