"""Render a versioned manuscript with Pandoc and Chrome; set CHROME if needed."""
import json,os,re,shutil,subprocess,tempfile
from pathlib import Path
from pypdf import PdfReader,PdfWriter
out=Path(__file__).resolve().parent;c=json.loads((out/'release.json').read_text());name=c['project']
article=subprocess.check_output(['pandoc','--from=gfm','--to=html5',str(out/f'{name}_preprint.md')],text=True)
article=re.sub(r'href="(?!https?:|#)([^"]+)"',lambda m:'href="https://github.com/GUIBA-EX/'+name+'/blob/preprint-v'+c['version']+'/'+str((out/m[1]).resolve().relative_to(out.parent))+'"',article)
css=(out/'preprint.css').read_text().replace('PROJECT',name)
chrome=os.environ.get('CHROME') or shutil.which('chromium') or shutil.which('google-chrome') or '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
with tempfile.TemporaryDirectory() as tmp:
 t=Path(tmp);doc=t/'paper.html';raw=t/'paper.pdf'
 doc.write_text(f'<!doctype html><meta charset=utf-8><base href="{out.as_uri()}/"><title>{c["title"]}</title><style>{css}</style>{article}')
 subprocess.run([chrome,'--headless','--disable-gpu','--allow-file-access-from-files','--no-pdf-header-footer',f'--user-data-dir={t}/profile',f'--print-to-pdf={raw}',doc.as_uri()],check=True)
 writer=PdfWriter(clone_from=PdfReader(raw));writer.add_metadata({'/Author':'Fei Xia','/Title':c['title'],'/Subject':f'Computational preprint v{c["version"]}; not peer reviewed'})
 with (out/f'{name}_preprint.pdf').open('wb') as f:writer.write(f)
