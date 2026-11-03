![MiAnemone 12S and 16S eDNA primers for sea anemones. Classification in progress](figures/MiAnemone_overview.png)

# MiAnemone · Mitochondrial 12S and 16S eDNA primers for sea anemones

[中文](README.md) | English

MiAnemone contains 22 computationally evaluated primer pairs for Actiniaria. Each marker has five size windows (XS–XL), each with an ordinary and a degenerate (D) pair sharing both binding sites. The 16S panel also includes H and HD at the original design site. Oligos are reported in synthesis orientation (5′→3′) and identified by marker, pair and direction. [Primer FASTA](primers.fasta) · [Sequences, coordinates and results](all_results.csv)

## Primer selection

**12S MD** covered 224 of 242 evaluable species (92.6%) with 208 exact matches (86.0%). Ordinary M covered 220 with 198 exact matches. Both hit nine background species under the strict Watson–Crick model. For a shorter product, SD covers 195/246 at 197 bp with zero hits across the four background models.

**16S LD** covered 231/258 species (89.5%) with 180 exact matches (69.8%), compared with 138 and 106 for ordinary L. LD hit zero and three background records under the two strict models, increasing to 46 and ten under the relaxed models.

High coverage can coincide with few exact matches. For example, 12S XSD covered 230/245 species with zero exact matches. Ordinary and D oligos were selected independently at the same coordinates.

## Methods and table key

Each target panel contains one GenBank sequence per named species. All coordinates are 1-based, inclusive positions in the respective ungapped mitochondrial rRNA gene of *Metridium senile* (NC_000933.1). Ordinary primers were selected from observed templates and consensus motifs following the conserved-flank approach of MiFish. D candidates incorporate lineage variation following the Anth-28S approach, with at most eight variants per oligo and ambiguity confined to its 3′ half.

Evaluable references span both sites and have no ambiguous nucleotides in their aligned binding regions. Both strands were searched for 100–500 bp products with at least 50 internal bases. Coverage allows two mismatches per oligo with compatibility at each of its final three bases. Exact matches require compatibility at every base on both ends.

Background scans allow 0–1,000 internal bases. Strict WC and G–T models allow two mismatches per oligo and require compatibility across its last three bases. Relaxed WC and G–T models allow three and two mismatches respectively without a terminal constraint. The selection thresholds were at most ten hits per strict model and 50 per relaxed model. Each table cell lists **ordinary / D** and each background entry lists **WC,G–T**. Each window has its own evaluable denominator.

## 12S

![12S binding positions](figures/12S_binding.svg)

The panel contains 250 target species and 322 other cnidarian species.

| Window | Span (bp) | Evaluable species | Covered | Exact | Strict background | Relaxed background |
|---|---:|---:|---:|---:|---:|---:|
| XS | 161 | 245 | 205 / 230 | 2 / 0 | 9,9 / 9,9 | 9,9 / 9,9 |
| S | 197 | 246 | 173 / 195 | 0 / 26 | 0,0 / 0,0 | 0,0 / 0,0 |
| M | 230 | 242 | 220 / 224 | 198 / 208 | 9,9 / 9,9 | 9,9 / 9,9 |
| L | 269 | 246 | 226 / 229 | 1 / 1 | 9,9 / 9,9 | 9,9 / 9,9 |
| XL | 385 | 243 | 185 / 225 | 46 / 52 | 0,0 / 0,0 | 1,0 / 18,0 |

[Per-species results](12S/species_results.csv) · [Background hits](12S/background_hits.csv) · [Coordinates](12S/coordinates.csv) · [Insert uniqueness](12S/barcode_identity.csv)

## 16S

![16S binding positions including H and HD](figures/16S_binding.svg)

The panel contains 265 target species and 1,775 other cnidarian records. Only 76–77 species are evaluable for the S, M and XL windows.

| Window | Span (bp) | Evaluable species | Covered | Exact | Strict background | Relaxed background |
|---|---:|---:|---:|---:|---:|---:|
| XS | 152 | 262 | 165 / 223 | 132 / 1 | 0,0 / 0,0 | 0,0 / 1,0 |
| S | 186 | 77 | 63 / 70 | 0 / 56 | 0,0 / 0,0 | 0,0 / 0,0 |
| M | 248 | 76 | 61 / 61 | 51 / 52 | 0,0 / 0,0 | 0,0 / 0,0 |
| L | 334 | 258 | 138 / 231 | 106 / 180 | 0,0 / 0,3 | 7,4 / 46,10 |
| XL | 404 | 76 | 41 / 61 | 36 / 52 | 0,0 / 0,0 | 0,0 / 0,0 |
| H / HD | 334 | 165 | 146 / 154 | 122 / 136 | 0,0 / 0,1 | 486,500 / 976,501 |

H and HD occupy the original design site (F 1497–1518, R 1809–1830). HD covers eight more species and exactly matches 14 more than H. Both exceed the relaxed background threshold.

[Per-species results](16S/species_results.csv) · [Background hits](16S/background_hits.csv) · [Coordinates](16S/coordinate_check.csv) · [Insert uniqueness](16S/insert_resolution.csv) · [16S sequences including H and HD](16S/primers.fasta)

## Reproduction and citation

Each marker directory includes the reference sequences, target alignment and accession panels. From the repository root, run `python3 scripts/binding_positions.py` to redraw the figures. Run `python3 scripts/verify.py` with Biopython installed to check per-record outcomes against both panels. There were 8,940 checks for 12S and 45,780 for 16S. [12S validation](12S/validation.json) · [16S validation](16S/validation.json)

Fei Xia. [Citation](CITATION.cff) · [Copyright and permissions](COPYRIGHT.md)

## References

1. Miya et al. (2015). MiFish, a set of universal PCR primers for metabarcoding environmental DNA from fishes. [Royal Society Open Science](https://doi.org/10.1098/rsos.150088)
2. Anth-28S primer design for coral metabarcoding. [PeerJ](https://doi.org/10.7717/peerj.18607)
