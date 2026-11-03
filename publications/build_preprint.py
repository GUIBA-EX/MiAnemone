"""Build the preprint with Python 3, Biopython and Chrome/Chromium (set CHROME if needed)."""
import csv
import html
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'publications'
def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))
def section(md,heading):
    return md.split('## '+heading+'\n',1)[1].split('\n## ',1)[0].strip()
def paragraphs(s):
    return ''.join('<p>'+html.escape(p).replace('*','')+'</p>' for p in s.split('\n\n'))
def table(headers,data):
    line=lambda tag,row:'<tr>'+''.join(f'<{tag}>{html.escape(str(c))}</{tag}>' for c in row)+'</tr>'
    return '<table>'+line('th',headers)+''.join(line('td',r) for r in data)+'</table>'
def fig(path):
    return '<img src="'+(ROOT/'figures'/path).as_uri()+'">'
def heading(x):return '<h2>'+x+'</h2>'
en=(ROOT/'README.en.md').read_text()
zh=(ROOT/'README.md').read_text()
pairs=read_csv(ROOT/'all_results.csv')
title=re.search(r'^# (.+)',en,re.M)[1]
pages=[]
pages.append('<h1>'+title+'</h1><p class="by">Fei Xia</p><small>Computational preprint · 4 November 2026</small>'+heading('Abstract')+paragraphs(section(en,'Abstract'))+'<p><b>Keywords</b> · Actiniaria, environmental DNA, mitochondrial rRNA, primer degeneracy, metabarcoding</p>'+heading('中文摘要')+paragraphs(section(zh,'摘要'))+heading('1. Introduction')+paragraphs(section(en,'Introduction'))+fig('MiAnemone_overview.png'))
methods=section(en,'Materials and methods').split('\n\n')
pages.append(heading('2. Materials and methods')+'<h3>2.1 Reference panels and primer design</h3>'+paragraphs(methods[0])+'<h3>2.2 Target coverage and non-target matching</h3>'+paragraphs(methods[1])+paragraphs(methods[2].split(' Each table cell')[0])+'''<h3>2.3 Base-wise compatibility</h3><p>For each primer end, we retained evaluable species with an unambiguous, contiguous aligned binding site of the reference length. Insertions and deletions within that end were excluded. We reverse-complemented reverse sites so all primer bases are read in synthesis orientation. Each site percentage is the fraction of retained species whose template base belongs to the primer’s IUPAC base set. This site-specific denominator can differ between ends. Paired coverage was evaluated separately by searching both strands of the raw sequences.</p>
<h3>2.4 Reference insert identity</h3><p>We removed both primer sites from the homologous interval and counted unambiguous inserts unique to one species within each window panel. A unique reference insert measures discrimination by sequence identity in the sampled panel.</p><h3>2.5 Independent verification</h3><p>Independent paired-hit checks on target and background sequences agreed across 8,940 record–model comparisons for 12S and 45,780 for 16S. Model counts use records in the corresponding marker panel.</p>''')
for i,g in enumerate(('12S','16S'),1):
    rr=[r for r in pairs if r['gene']==g]
    stat=[[r['name'],r['evaluable'],r['covered'],r['exact'],r['WC']+' / '+r['GT'],r['WC3']+' / '+r['GT2']] for r in rr]
    p=heading(f'3.{i} Results · mitochondrial {g}')+('<p>The 12S panel comprises 250 target species and 322 non-target cnidarian species. MD covered 224/242 species with 208 exact matches. Both M and MD hit nine background records under strict WC pairing.</p>' if i==1 else '<p>The 16S panel comprises 265 target species and 1,775 non-target records. LD covered 231/258 species with 180 exact matches. H and HD covered 146/165 and 154/165 at the original site, but both exceeded the relaxed background thresholds.</p>')
    p+=fig(g+'_binding.svg')+f'<p class="cap">Figure {i}. {g} reference binding sites on NC_000933.1. Arrows point in the synthesis direction. Spans include both primers.</p>'
    p+=f'<p class="cap">Table {i}. Predicted paired matches among evaluable target species (N). Background cells show WC / G–T record counts for each mismatch rule.</p>'+table(['Pair','N','Covered','Exact','Strict background','Relaxed background'],stat)
    pages.append(p)
for i,g in enumerate(('12S','16S'),3):
    pages.append(heading(f'3.{i} Results · {g} base matching')+fig(g+'_base_matching.svg')+f'<p class="cap">Figure {i}. Each cell gives the synthesis base and rounded percentage of compatible reference species. Rows are ordinary and D forward ends followed by ordinary and D reverse ends. Colors mark &lt;50%, 50–&lt;90%, 90–&lt;100% and 100% compatibility. The 3′ triplet is outlined. Counts for the two ends may differ. Exact numerators and denominators are provided with the figure data.</p>')
identity=[]
for g,file in [('12S','barcode_identity'),('16S','insert_resolution')]:
    for r in read_csv(ROOT/g/(file+'.csv')):
        identity.append([g,r['window'],r['unambiguous'],r.get('unique',r.get('identity_unique'))])
pages.append(heading('4. Discussion')+paragraphs(section(en,'Discussion'))+'<p>The 12S XL window had unique inserts in 100/239 species with unambiguous inserts, compared with 35/244 for XS. In 16S, the H window had 87 unique inserts among 159 species. Different taxa contribute to each denominator, so these counts cannot isolate an effect of length alone.</p><p class="cap">Table 3. Primer-free reference insert identity by window.</p>'+table(['Marker','Window','Unambiguous inserts','Unique to one species'],identity))
pages.append(heading('5. Primer sequences')+'<p>Both ends are written in the 5′ to 3′ synthesis direction. Each D pair shares the four reference endpoints of its ordinary counterpart. Degeneracy counts the number of oligo variants per end.</p><p class="cap">Table 4. All 22 primer pairs in the released FASTA.</p>'+table(['Marker','Pair','Forward (5′–3′)','Reverse (5′–3′)','Deg F/R'],[[r['gene'],r['name'],r['F'],r['R'],r['deg_F']+'/'+r['deg_R']] for r in pairs])+heading('Data availability')+'<p>Accession panels, alignments, oligos, per-species predictions and the scripts for verification and figure regeneration are archived in the <a href="https://github.com/GUIBA-EX/MiAnemone">MiAnemone repository</a>.</p>'+heading('References')+''.join('<p class="ref">'+re.sub(r'\[DOI\]\((.*?)\)',r'<a href="\1">DOI</a>',line.replace('*',''))+'</p>' for line in section(en,'References').splitlines()))
css="""@page{size:A4;margin:0}*{box-sizing:border-box}body{margin:0;font:11.3pt/1.44 Georgia,'Times New Roman',serif;color:#182722}.page{height:297mm;padding:16mm 17mm 18mm;position:relative;break-after:page}.page:last-child{break-after:auto}h1{font-size:22pt;line-height:1.16;text-align:center;margin:0 0 5mm}h2{font-size:16pt;margin:4mm 0 2mm}h3{font-size:12pt;margin:4mm 0 1mm}p{margin:2mm 0 3mm}.by{font-weight:bold}small,.cap,.ref{font-size:9.2pt}.cap{color:#345149}.page img{width:100%;display:block;margin:4mm auto 2mm}.page:first-child img{width:77%}table{width:100%;border-collapse:collapse;font:9pt/1.35 Arial,sans-serif;margin:2mm 0}th,td{padding:1.45mm;border-bottom:1px solid #dae4df;text-align:left}th{background:#e7efed}tr:nth-child(even){background:#f7faf9}footer{position:absolute;bottom:8mm;left:17mm;right:17mm;border-top:.6pt solid #b2c3ba;padding-top:2mm;font:8pt Arial;color:#53635d}footer span{float:right}a{color:#235954}.page:last-child table{font-size:8pt}.page:last-child td{padding:1.1mm}.page:last-child td:nth-child(3),.page:last-child td:nth-child(4){font-family:monospace}"""
body=''.join(f'<section class="page">{p}<footer>MiAnemone · Fei Xia · Computational preprint · 4 November 2026 <span>{j}</span></footer></section>' for j,p in enumerate(pages,1))
chrome=os.environ.get('CHROME') or shutil.which('chromium') or '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome'
with tempfile.TemporaryDirectory() as temp:
    doc=Path(temp)/'preprint.html'
    doc.write_text(f'<!doctype html><meta charset="utf-8"><title>{title}</title><style>{css}</style>{body}')
    subprocess.run([chrome,'--headless','--disable-gpu','--no-pdf-header-footer','--allow-file-access-from-files',f'--user-data-dir={temp}/profile',f'--print-to-pdf={OUT}/MiAnemone_preprint.pdf',doc.as_uri()],check=True)
