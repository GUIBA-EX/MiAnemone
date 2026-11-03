"""Check published pair counts on both strands of each underlying sequence."""
import csv
import json
from pathlib import Path
from independent import flags

ROOT = Path(__file__).resolve().parents[1]

def fasta(path):
    sequences = {}
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            key = line[1:].split()[0]
            assert key not in sequences
            sequences[key] = ''
        else:
            sequences[key] += line.strip().upper()
    return sequences

def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f))

def verify(gene):
    folder = ROOT / gene
    targets = fasta(folder/'target.fasta')
    background = fasta(folder/'background.fasta')
    aligned = {k.removeprefix('_R_'):v for k,v in fasta(folder/'aligned.fasta').items()}
    cols = [i for i,c in enumerate(aligned['NC_000933.1']) if c!='-']
    details = {(r['name'],r['accession']):r for r in rows(folder/'species_results.csv')}
    hits = {(r['name'],r['accession'],r['model']):r for r in rows(folder/'background_hits.csv')}
    pairs = rows(folder/'candidate_results.csv')
    checks = 0
    for r in pairs:
        name = r['name']
        a,b,c,d = [cols[int(r[k])-1] for k in ('F_start','F_end','R_start','R_end')]
        counts = dict.fromkeys(('evaluable','covered','exact','WC','GT','WC3','GT2'),0)
        for acc,seq in targets.items():
            aln = aligned[acc]
            ok = bool(aln[:a+1].strip('-') and aln[d:].strip('-') and set(aln[a:b+1]+aln[c:d+1])<=set('ACGT-'))
            exact, covered, _ = flags(seq,r['F'],r['R'],True,False)
            detail = details[name,acc]
            for key,value in [('evaluable',ok),('covered',ok and covered),('exact',ok and exact)]:
                assert int(detail[key])==value,(name,acc,key)
                counts[key]+=value
            checks+=1
        for acc,seq in background.items():
            for model,loose in [('WC','WC3'),('GT','GT2')]:
                _, strict, broad = flags(seq,r['F'],r['R'],False,model=='GT')
                detail = hits.get((name,acc,model),{})
                assert (int(detail.get('strict',0)),int(detail.get('relaxed',0)))==(strict,broad)
                counts[model]+=strict
                counts[loose]+=broad
                checks+=1
        assert all(counts[k]==int(r[k]) for k in counts),name
    validation=dict(record_model_checks=checks,disagreements=0,background_records=len(background),target_species=len(targets))
    (folder/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(gene,checks,'checks passed',flush=True)

if __name__=='__main__':
    for gene in ('12S','16S'):
        verify(gene)
