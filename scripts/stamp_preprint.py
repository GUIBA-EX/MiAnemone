"""Stamp a new version before the PDF is built."""
import datetime
from pathlib import Path
import re
import sys

version = sys.argv[1]
assert re.fullmatch(r'1\.\d+',version)
day = datetime.date.today()
english = f'{day.day} {day:%B %Y}'
for filename in ['README.md','README.en.md','CITATION.cff','publications/MiAnemone_preprint.md','publications/build_preprint.py']:
    path = Path(filename)
    text = path.read_text()
    text = re.sub(r'preprint-v1\.\d+',f'preprint-v{version}',text)
    text = re.sub(r'(Preprint v|预印本 v|version: "|version )1\.\d+',lambda m:m[1]+version,text)
    text = re.sub(r'date-released: "[\d-]+"',f'date-released: "{day}"',text)
    text = re.sub(r'^  year: \d+$',f'  year: {day.year}',text,flags=re.M)
    text = re.sub(r'(Computational preprint · |Preprint v[\d.]+ · )\d+ [A-Za-z]+ \d{4}',lambda m:m[1]+english,text)
    text = re.sub(r'(预印本 v[\d.]+ · )\d+ 年 \d+ 月 \d+ 日',lambda m:m[1]+f'{day.year} 年 {day.month} 月 {day.day} 日',text)
    path.write_text(text)
