"""Calculate column-wise primer compatibility on homologous, oriented 16S and 12S alignments."""
import csv
from pathlib import Path
from xml.etree import ElementTree as ET

from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Seq import Seq
from verify import fasta, rows

ROOT = Path(__file__).resolve().parents[1]
KEYS = ('F_start', 'F_end', 'R_start', 'R_end')
COLORS = ['#fff0cf', '#cde7e4', '#37877f', '#14524e']


def sites(alignment, columns, start, end, accessions, reverse):
    # Retain contiguous, unambiguous sites with no insertion or deletion relative to the reference.
    indices = columns[start-1:end]
    words = []
    for accession in accessions:
        s = alignment[accession]
        segment = s[indices[0]:indices[-1]+1].replace('-', '')
        word = ''.join(s[i] for i in indices)
        if word == segment and set(word) <= set('ACGT'):
            words.append(str(Seq(word).reverse_complement()) if reverse else word)
    return words


def element(group, tag, **attrs):
    return ET.SubElement(group, tag, {k.replace('_', '-'):str(v) for k,v in attrs.items()})


def text(g,x,y,label,size,**attrs):
    element(g,'text',x=x,y=y,font_size=size,**attrs).text = str(label)


def build(gene):
    folder = ROOT/gene
    align = {k.removeprefix('_R_'):v for k,v in fasta(folder/'aligned.fasta').items()}
    columns = [i for i,c in enumerate(align['NC_000933.1']) if c!='-']
    panel = rows(folder/'species_results.csv')
    pairs = {r['name']:r for r in rows(folder/'candidate_results.csv')}
    windows = ['XS','S','M','L','XL'] + (['H'] if gene=='16S' else [])
    svg = ET.Element('svg', xmlns='http://www.w3.org/2000/svg', viewBox=f'0 0 850 {66+len(windows)*138}')
    g = element(svg,'g',font_family='Arial, sans-serif',fill='#183a3a')
    text(g,5,22,f'{gene}  |  Primer–template compatibility by nucleotide',20,font_weight='bold')
    text(g,5,44,'5′ → 3′ in each row. Cells show the primer base and percent compatible; the final three bases are outlined.',13)
    detail = []
    for index,window in enumerate(windows):
        y = 66+138*index
        text(g,5,y+9,window+' / '+window+'D',16,font_weight='bold')
        a,b,c,d = (int(pairs[window][k]) for k in KEYS)
        text(g,145,y+9,f'F {a}–{b}     R {c}–{d}',13)
        for rowIndex,(name,side) in enumerate([(window,'F'),(window+'D','F'),(window,'R'),(window+'D','R')]):
            pair = pairs[name]
            assert all(pair[k]==pairs[window][k] for k in KEYS)
            pool = [r['accession'] for r in panel if r['name']==name and r['evaluable']=='1']
            q = pair[side]
            words = sites(align,columns,int(pair[side+'_start']),int(pair[side+'_end']),pool,side=='R')
            assert words and all(len(w)==len(q) for w in words)
            top = y+16+rowIndex*26
            text(g,5,top+17,f'{name} {side}',13)
            text(g,64,top+17,f'n={len(words)}',12)
            for j,base in enumerate(q):
                count=sum(w[j] in ambiguous_dna_values[base] for w in words)
                pct = 100*count/len(words)
                color = COLORS[3 if count==len(words) else 2 if pct>=90 else 1 if pct>=50 else 0]
                x = 132+j*31
                element(g,'rect',x=x,y=top,width=30,height=25,fill=color)
                ink='white' if pct>=90 else '#183a3a'
                text(g,x+15,top+11,base,12,text_anchor='middle',fill=ink,font_weight='bold')
                text(g,x+15,top+22,f'{pct:.0f}',10,text_anchor='middle',fill=ink)
                detail.append((gene,name,side,j+1,base,(int(pair[side+'_end'])-j if side=='R' else int(pair[side+'_start'])+j),len(words),count))
            element(g,'rect',x=132+(len(q)-3)*31,y=top,width=92,height=25,fill='none',stroke='#242c38',stroke_width=1.4)
    ET.ElementTree(svg).write(ROOT/'figures'/f'{gene}_base_matching.svg',encoding='utf-8',xml_declaration=True)
    return detail


if __name__=='__main__':
    with (ROOT/'publications/base_matching_counts.csv').open('w',newline='') as f:
        writer=csv.writer(f,lineterminator='\n')
        writer.writerow(('gene','pair','side','position_5to3','primer_base','reference_coordinate','site_species','compatible_species'))
        for gene in ('12S','16S'):
            writer.writerows(build(gene))
