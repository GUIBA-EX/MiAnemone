![MiAnemone 12S and 16S eDNA primers for sea anemones. Taxonomic framework under development](figures/MiAnemone_overview.png)

# MiAnemone, a set of universal PCR primers for metabarcoding environmental DNA from sea anemones

**Fei Xia**

[中文](README.md) | English

## Abstract

We designed and evaluated 22 primer pairs in silico for mitochondrial 12S and 16S rRNA genes in Actiniaria. In a 12S reference library of 250 sea anemone species, MD covered 224 of 242 evaluable species and matched both binding sites exactly in 208. In a 16S library of 265 species, LD covered 231 of 258 evaluable species, with 180 exact matches. Under strict Watson–Crick pairing, MD and LD matched nine and zero records in their respective background panels.

## Introduction

Primer design for environmental DNA metabarcoding must account for amplicon length, variation at binding sites among target species and matches to related non-target taxa. We drew on the conserved-flank design of MiFish and the use of degenerate sites in Anth-28S (Miya et al., 2015, McCartin et al., 2024). We compared ordinary and degenerate primers at identical coordinates in five size windows (XS–XL) per mitochondrial marker to determine how gains in coverage affect exact matching and background hits.

## Materials and methods

Each target panel contains one GenBank sequence per named species. All coordinates are 1-based, inclusive positions in the respective ungapped mitochondrial rRNA gene of *Metridium senile* (NC_000933.1). Ordinary primers were selected from observed templates and consensus motifs following the conserved-flank approach of MiFish. D candidates incorporate lineage variation following the Anth-28S approach, with at most eight variants per oligo and ambiguity confined to its 3′ half.

Evaluable references span both sites and have no ambiguous nucleotides in their aligned binding regions. Both strands were searched for 100–500 bp products with at least 50 internal bases. Coverage allows two mismatches per oligo with compatibility at each of its final three bases. Exact matches require compatibility at every base on both ends.

Background scans allow 0–1,000 internal bases. Strict WC and G–T models allow two mismatches per oligo and require compatibility across its last three bases. Relaxed WC and G–T models allow three and two mismatches respectively without a terminal constraint. The selection thresholds were at most ten hits per strict model and 50 per relaxed model. Each table cell lists **ordinary / D** and each background entry lists **WC,G–T**. Each window has its own evaluable denominator.

## Results

### 12S

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

### 16S

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

## Discussion

12S MD yielded the most exact matches among the 12S pairs. Compared with ordinary M, it covered four additional species and provided ten additional exact matches with the same strict background count. Where a shorter amplicon and fewer background hits are priorities, the 197 bp SD pair covered 195/246 species with zero hits in all four background models. For 16S, LD covered 93 more species than ordinary L, while relaxed WC hits increased from seven to 46.

Shorter windows illustrate the distinction between coverage and exact matching. 12S XSD covered 230/245 species with zero exact matches. For 16S, XSD covered 58 more species than XS, but exact matches fell from 132 to one. Ordinary and D sequences were selected independently at each site, so both metrics matter when choosing a pair. Comparisons across windows also depend on how many references span both sites, particularly for the sparsely represented 16S windows. These results describe predicted template compatibility. Amplification efficiency in environmental samples remains to be measured.

## Data availability and reproducibility

Oligos are reported in synthesis orientation (5′→3′) with identifiers for gene, pair and direction, for example `12S_MD_F`. The full results table includes sequences, binding coordinates and background counts. [Primer FASTA](primers.fasta) · [Complete results](all_results.csv)

Each marker directory contains target and background references, accession panels and the target alignment. With Biopython installed, run `python3 scripts/verify.py` from the repository root to verify hits, and `python3 scripts/binding_positions.py` to redraw the binding maps. All 8,940 checks for 12S and 45,780 checks for 16S agreed with the independent implementation. [12S validation](12S/validation.json) · [16S validation](16S/validation.json)

[Citation](CITATION.cff) · [Copyright and permissions](COPYRIGHT.md)

## References

1. Miya M et al. (2015). MiFish, a set of universal PCR primers for metabarcoding environmental DNA from fishes: detection of more than 230 subtropical marine species. *Royal Society Open Science*, 2, 150088. [DOI](https://doi.org/10.1098/rsos.150088)
2. McCartin L et al. (2024). Nuclear eDNA metabarcoding primers for anthozoan coral biodiversity assessment. *PeerJ*, 12, e18607. [DOI](https://doi.org/10.7717/peerj.18607)
