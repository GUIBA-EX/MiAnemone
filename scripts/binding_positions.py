"""Draw paired primer windows on the gene-relative NC_000933.1 coordinate scale."""
import csv
from pathlib import Path
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
COLORS = {'XS': '#288064', 'S': '#13aaba', 'M': '#c88312', 'L': '#2779ef',
          'XL': '#922bde', 'H': '#697885'}
KEYS = ('F_start', 'F_end', 'R_start', 'R_end')

def element(parent, tag, **attrs):
    return ET.SubElement(parent, tag, {k.replace('_', '-'): str(v) for k, v in attrs.items()})

def draw(gene, rows):
    groups = {}
    for row in rows:
        coords = tuple(int(row[k]) for k in KEYS)
        assert coords == groups.setdefault(row['group'], coords)
        a, b, c, d = coords
        assert d-a+1 == int(row['span']) and b-a+1 == len(row['F']) and d-c+1 == len(row['R'])
    order = [k for k in COLORS if k in groups]
    lo = min(v[0] for v in groups.values()) // 50 * 50
    hi = (max(v[3] for v in groups.values()) + 49) // 50 * 50
    x = lambda n: 155 + 830 * (n-lo) / (hi-lo)
    height = 160 + len(order)*65
    svg = ET.Element('svg', xmlns='http://www.w3.org/2000/svg', width='1100', height=str(height),
                     viewBox=f'0 0 1100 {height}', role='img')
    element(svg, 'title').text = f'MiAnemone {gene} primer binding positions'
    element(svg, 'rect', width='100%', height='100%', rx=18, fill='#ffffff')
    g = element(svg, 'g', font_family='Avenir Next, Arial, sans-serif', fill='#173647')
    def label(px, y, s, size, **attrs):
        element(g, 'text', x=px, y=y, font_size=size, **attrs).text = s
    label(30, 40, f'MiAnemone  |  {gene} binding positions', 27, font_weight='bold')
    reference = 'Halcampoides purpureus  ·  MT676790.1' if gene == '28S' else 'Metridium senile  ·  NC_000933.1'
    label(30, 69, reference + '  ·  gene-relative coordinates (1-based, inclusive)', 17)
    for n in range(lo, hi+1, 50):
        element(g, 'line', x1=x(n), x2=x(n), y1=102, y2=height-40, stroke='#dce7ef')
        label(x(n), 97, str(n), 13, text_anchor='middle')
    for i, key in enumerate(order):
        y = 134 + i*65
        a,b,c,d = groups[key]
        color = COLORS[key]
        label(24, y+5, key+' / '+key+'D', 17, font_weight='bold')
        element(g, 'line', x1=x(a), x2=x(d), y1=y, y2=y, stroke=color, stroke_opacity='.35', stroke_width=5)
        for start,end in [(a,b),(d,c)]:
            tail, tip = x(start), x(end)
            neck = tip + (1 if tail>tip else -1)*min(7,abs(tip-tail)/2)
            element(g,'polygon',points=f'{tail},{y-4} {neck},{y-4} {neck},{y-8} {tip},{y} {neck},{y+8} {neck},{y+4} {tail},{y+4}',fill=color)
        label((x(a)+x(d))/2, y-10, f'{d-a+1} bp', 15, text_anchor='middle')
        for mid, s in [((a+b)/2,f'F {a}–{b}'), ((c+d)/2,f'R {c}–{d}')]:
            label(x(mid), y+23, s, 14, text_anchor='middle')
    footer = 'Arrows show primer orientation (5′ to 3′). Spans include both primers.'
    label(24,height-17,footer,14)
    ET.ElementTree(svg).write(ROOT/'figures'/f'{gene}_binding.svg', encoding='utf-8', xml_declaration=True)

if __name__=='__main__':
    with (ROOT/'all_results.csv').open() as f:
        rows=list(csv.DictReader(f))
    for gene in ('12S','16S','28S'):
        draw(gene,[r for r in rows if r['gene']==gene])
