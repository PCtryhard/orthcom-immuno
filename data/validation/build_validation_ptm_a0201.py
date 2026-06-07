#!/usr/bin/env python3
"""Extract A*02:01 modified-vs-native IC50 pairs from MOESM3 .docx -> CSV.

Deterministic: reads the actual Word table XML. A '-' cell means the peptide is a
NON-BINDER (IC50 above the assay ceiling; the paper's Methods state '-' = IC50
>40,000 nM). It becomes a BLANK numeric field (never a guessed/interpolated number)
and the pair is labelled *_nonbinder. Peptide modification position is derived by
diffing native vs modified sequence.
Source: Sidney & Sette et al. 2018, BMC Immunol, PMC5863483, Additional file 3.
"""
import sys, zipfile, csv, os
import xml.etree.ElementTree as ET
from collections import Counter

W = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
docx, out_csv = sys.argv[1], sys.argv[2]

with zipfile.ZipFile(docx) as z:
    root = ET.fromstring(z.read('word/document.xml'))

def para_text(p): return ''.join(n.text or '' for n in p.iter(f'{W}t'))
def cell_text(tc): return ' '.join(para_text(p) for p in tc.findall(f'{W}p')).strip()

rows = root.find(f'.//{W}tbl').findall(f'{W}tr')

DASHES = {'', '-', '–', '—'}
def num_or_blank(s):
    s = s.strip()
    if s in DASHES: return ''
    try:
        float(s.replace(',', '')); return s.replace(',', '')
    except ValueError:
        return ''  # anything non-numeric -> blank (flagged via raw column)

def diff_positions(native, modified):
    if len(native) != len(modified):
        return '', 'len_mismatch'
    pos = [str(i + 1) for i, (a, b) in enumerate(zip(native, modified)) if a != b]
    return ';'.join(pos), ('' if pos else 'no_diff')

records, legend, skips = [], '', []
for ri, tr in enumerate(rows):
    cells = [cell_text(tc) for tc in tr.findall(f'{W}tc')]
    if ri in (0, 1):            # group-header + column-header rows
        continue
    if len(cells) == 1:         # legend row
        legend = cells[0]; continue
    if len(cells) < 9:
        skips.append(f'r{ri}: {len(cells)} cells -> {cells}'); continue
    start, segment, native, length, reported, modified, modtype, mod_raw, wt_raw = cells[:9]
    modpos, posflag = diff_positions(native, modified)
    nat_ic50, mod_ic50 = num_or_blank(wt_raw), num_or_blank(mod_raw)
    if nat_ic50 and mod_ic50:           pair = 'complete'
    elif not nat_ic50 and not mod_ic50: pair = 'both_nonbinder'
    elif not mod_ic50:                  pair = 'modified_nonbinder'
    else:                               pair = 'native_nonbinder'
    notes = [x for x in (posflag, ('multi_site' if modtype.lower().startswith('pan') else '')) if x]
    records.append({
        'native_sequence': native, 'modified_sequence': modified,
        'modified_position': modpos, 'modification_type': modtype,
        'native_IC50_nM': nat_ic50, 'modified_IC50_nM': mod_ic50,
        'allele': 'HLA-A*02:01', 'length': length, 'insulin_segment': segment,
        'insulin_start': start, 'reported_A0201_epitope': reported,
        'pair_status': pair, 'native_IC50_raw': wt_raw, 'modified_IC50_raw': mod_raw,
        'provenance': f'Sidney&Sette2018 PMC5863483 Add.file3(MOESM3) tbl-row r{ri}',
        'notes': ';'.join(notes),
    })

cols = ['native_sequence','modified_sequence','modified_position','modification_type',
        'native_IC50_nM','modified_IC50_nM','allele','length','insulin_segment',
        'insulin_start','reported_A0201_epitope','pair_status',
        'native_IC50_raw','modified_IC50_raw','provenance','notes']
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
with open(out_csv, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols, lineterminator='\n')
    w.writeheader(); w.writerows(records)

print('LEGEND:', legend)
print('rows written:', len(records))
print('pair_status:', dict(Counter(r['pair_status'] for r in records)))
print('modification_type:', dict(Counter(r['modification_type'] for r in records)))
if skips: print('SKIPS:', skips)
print('\n--- COMPLETE pairs (both native & modified IC50 numeric) ---')
n = 0
for r in records:
    if r['pair_status'] == 'complete':
        n += 1
        print(f"{r['native_sequence']:<11} {r['native_IC50_nM']:>6} -> "
              f"{r['modified_sequence']:<11} {r['modified_IC50_nM']:>6}  "
              f"{r['modification_type']:<16} P{r['modified_position']}")
print('COMPLETE pairs total:', n)
