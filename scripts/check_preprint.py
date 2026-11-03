"""Fail publication when manuscript tables drift from the released result table."""
import csv
from pathlib import Path
import sys

root = Path(__file__).resolve().parents[1]
paper = (root/'publications/MiAnemone_preprint.md').read_text()
tables = paper.split('## Table ')
with (root/'all_results.csv').open() as file:
    pairs = list(csv.DictReader(file))
fields = [('gene','name','F','R','span'),
          ('gene','name','evaluable','covered','exact','WC','GT','WC3','GT2')]
for table, keys in zip(tables[1:], fields, strict=True):
    actual = [line.strip('|').strip().split(' | ') for line in table.splitlines() if line.startswith('|')][2:]
    assert actual == [[row[key] for key in keys] for row in pairs], 'Update the manuscript result tables'
fasta = (root/'primers.fasta').read_text().splitlines()
expected = {f"{r['gene']}_{r['name']}_{end}":r[end] for r in pairs for end in 'FR'}
assert dict(zip([k[1:] for k in fasta[::2]],fasta[1::2],strict=True)) == expected
if len(sys.argv)>1:
    text = Path(sys.argv[1]).read_text()
    assert all(s in text for s in expected.values()), 'Missing primer sequence in rendered PDF'
