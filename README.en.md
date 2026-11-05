![MiAnemone overview](figures/MiAnemone_overview.png)

# MiAnemone: in silico design and evaluation of mitochondrial and nuclear rRNA primers for sea anemone (Actiniaria) eDNA metabarcoding

**Fei Xia**

<!-- publication:start -->
**Preprint v1.5 · 5 November 2026 · Not peer reviewed.**

[PDF](publications/MiAnemone_preprint.pdf) · [Manuscript](publications/MiAnemone_preprint.md) · [Versioned release](https://github.com/GUIBA-EX/MiAnemone/releases/tag/preprint-v1.5) · [CITATION.cff](CITATION.cff) · [BibTeX](CITATION.bib) · [Build and updates](publications/README.md)
<!-- publication:end -->


[中文](README.md) · [Revision PDF](publications/MiAnemone_preprint.pdf) · [Response to review](publications/response_to_review.md) · [Reproduction guide](revision/README.md)

## Start here

| Intended use | Priority |
| --- | --- |
| Short fragments | 12S XSD |
| Matching and reference discrimination | 12S MD |
| Broader target coverage | 16S LD |
| Complement mitochondrial barcode discrimination | 28S LD; background amplification needs evaluation |

[Primer sequences](primers.fasta) · [Primer diagrams](#primers-and-base-matching-diagrams) · [PDF](publications/MiAnemone_preprint.pdf)

## Abstract

Sea-anemone environmental DNA surveys require primers that accommodate target variation while retaining informative barcode regions. We evaluated 22 ordinary and degenerate primer pairs across five mitochondrial 12S and 16S rRNA windows and an additional manually designed 16S site. Reference panels contained 250 and 265 named species, respectively. Paired matching, common-species comparisons, multi-record analyses and genus-blocked reselection at fixed windows assessed compatibility and reference discrimination. The adopted 12S XSD covered 227 of 245 evaluable species and exactly matched 212, with nine strict background hits. Factorial controls showed that fixed-base substitutions, rather than degeneracy itself, explained the preceding candidate's loss of exact matching. Common coverage panels contained 240 species for 12S and 67 for 16S, highlighting differences in reference availability among windows. Reference uniqueness depended on barcode location as well as length. The results prioritize 12S XSD for short-fragment studies, 12S MD for matching and reference discrimination, and 16S LD for broader target coverage where relaxed background matching is acceptable. These computational comparisons provide a basis for selecting assays according to target taxa and required resolution.

Eight nuclear 28S pairs bring the catalogue to 30 pairs. 28S LD covers 166/172 species with 158/172 exact matches. Of 67 mitochondrial-conflict pairs evaluable at 28S L, 64 have distinct inserts and LD covers both species in 63. Strict background hits are zero, but relaxed WC/GT hits reach 340/339 species; specificity therefore depends strongly on terminal constraints.

**Keywords:** Actiniaria; environmental DNA; mitochondrial and nuclear rRNA; degenerate primers; reference databases; in silico PCR

## Introduction

Actiniaria comprises a diverse assemblage of benthic anthozoans whose detection and identification can be difficult from fragmented, cryptic or inaccessible specimens. Rodríguez et al. (2014) recognized approximately 1,200 described species and reorganized the order into Anenthemonae and Enthemonae using molecular evidence that included mitochondrial 12S and 16S. Environmental DNA offers a complementary route to surveying this diversity, but a candidate assay must solve two distinct problems: compatible primer-binding sites and a taxonomically informative barcode.

MiFish illustrates the use of conserved flanks around an informative mitochondrial region (Miya et al., 2015). Anth-28S instead uses degenerate nuclear-ribosomal primers to accommodate anthozoan variation (McCartin et al., 2024). Neither principle guarantees that an analogous mitochondrial assay will distinguish sea-anemone species. Slow mitochondrial sequence evolution in Anthozoa is well documented (Shearer et al., 2002; Hellberg, 2006), and shared mitochondrial sequences can limit assignment even when primer compatibility is high.

Mismatch-tolerant coverage also differs from exact matching. Experimental metabarcoding studies show that primer mismatches can influence amplification bias (Piñol et al., 2015), but their effects depend on position and reaction conditions. We use mismatch reduction to quantify sequence compatibility. Following the comparative approach illustrated by Bae (2025), this study evaluates binding-site compatibility, background matching and barcode discrimination separately. We ask how these assessments change after controlling the species panel, separating degeneracy from fixed-base substitutions, and testing additional records and withheld genera.

## Materials and methods

### Data acquisition and reference curation

The archived GenBank searches retrieved 1,360 records for 12S and 2,527 for 16S. The complete search expressions and accession-level decisions are distributed with the frozen input archives. Searches combined Actiniaria with marker names and mitochondrial-genome terms. Annotated rRNA features were preferred; restricted title-based extraction accommodated unannotated marker fragments. The initial eligible pools contained 786 and 1,245 records, respectively. GenBank organism names and lineages were retained. Records labelled sp., cf. or other non-binomial names were excluded from the named-species design panels. These names were not independently reconciled against WoRMS, and voucher identifications were not revalidated; named-species counts consequently describe the archived GenBank labels.

One reference was retained per named species. For 12S, ranking used the number of called A/C/G/T bases, rRNA annotation, RefSeq status and accession as successive tie-breakers. Complete-genome status was not used in the 12S ranking. For 16S, ranking prioritized genome-supported completeness, called reference columns, called sequence length, low ambiguity, feature-based extraction, RefSeq status and accession. The resulting design panels contained 250 12S species and 265 16S species. Relative to the approximate diversity reported in 2014, these are about one fifth of described Actiniaria, not comprehensive inventories of the order.

The original background panels comprised 322 other cnidarian species for 12S and 1,775 other cnidarian records for 16S. Their numerical units differ; we retain explicit record and species labels and do not compare raw hit counts across markers as specificity estimates. Background results describe these cnidarian panels.

### Alignment and reference coordinates

MAFFT (Katoh & Standley, 2013) was run with `--adjustdirection --retree 2 --maxiterate 0 --thread 4`. The alignment used nucleotide sequence similarity. Binding-site coordinates are one-based and inclusive on the ungapped genes of *Metridium senile* NC_000933.1. The reference 12S feature occupies mitochondrial positions 12,201–13,282 and the 16S feature positions 72–2,260, both on the positive strand. We use 12S and 16S throughout for mitochondrial small- and large-subunit rRNA, respectively.

A sequence was evaluable at a window when its alignment spanned both binding regions and those regions contained no ambiguous nucleotide codes; alignment gaps were allowed. This criterion determines the denominator. Paired raw-sequence scanning then determines whether a product is predicted. A complete homologous binding motif is a stricter subset: after removing alignment gaps, each binding segment must have the primer length and contain only A/C/G/T. Thus evaluability and membership in the motif-frequency calculation need not coincide.

### Initial candidate generation and selection

The original search examined reference start positions at three-base increments and oligo lengths of 20 or 22 nt. At least ten species, or 10% of the target panel if greater, had to provide complete motifs. Candidate ordinary oligos comprised the position-wise consensus and frequent observed motifs. Degenerate alternatives added bases supported by at least 4% of complete motifs in the 3′ half, with at most eight sequence variants per end. Oligos were ranked by compatible-species count, then exact-species count, then lower degeneracy. The search retained a small ranked subset at each start position; it was a heuristic search, not an exhaustive proof of optimality.

Candidate pairs spanned 150–470 reference bases. Coordinate bins of 25 bp retained spatially diverse candidates, and successive size bands defined XS, S, M, L and XL. Raw target and background screening determined the selected windows. Within each selected window, observed motifs and a consensus were reconsidered. Degenerate candidates were built from up to 12 high-ranking ordinary seeds by merging compatible motifs sharing the same 5′ half for up to three rounds. Up to six ordinary and six degenerate candidates per end entered paired screening. Ordinary and D pairs were selected separately at identical endpoints. H denotes the hand-designed 16S site; H and HD are revised sequences at that site, retained as comparators rather than preferred assays.

Terminal filters rejected oligos capable of a homopolymer of five bases or a homogeneous terminal triplet, required at least one assured G/C among the final three positions, and limited possible G/C content in the final five positions. There was no additional global GC-percentage gate. A degenerate terminal base was permitted when all other constraints were met; compatibility means that at least one concrete variant supplies a matching terminus. This permits terminal variants while diluting the concentration of each matching oligo.

Thermodynamic filtering used Primer3 (Untergasser et al., 2012) through primer3-py 2.4.0 with 50 mM monovalent salt, 1.5 mM divalent salt and 0.6 mM dNTP. Nearest-neighbor Tm calculations used the SantaLucia parameter and salt-correction defaults (SantaLucia, 1998). Scenarios were 50 nM per concrete variant and 200 nM per direction divided equally among variants. The latter gives 25 nM per variant at degeneracy eight. Minimum variant Tm was 45°C, maximum modeled secondary-structure Tm was below 55°C, and the initial allowed F/R difference was 5°C. These thresholds defined candidate screening. Oligo properties and both concentration scenarios are reported in the supplementary table.

The initial ranking and thermodynamic gates can retain non-consensus bases. In particular, candidates were truncated before all paired combinations were examined, and a fixed temperature-difference threshold can exclude a consensus pair. This explains why shared conserved sites can have different selected bases; it is not evidence that each selected sequence is a universal consensus. The revision therefore audits every position, reports rejected alternatives, and tests training-only reselection rather than interpreting the initial winner as a global optimum.

### Matching models and continuous mismatch diagnostics

Coverage denotes a paired hit with at most two mismatches per end and zero mismatches in each terminal three nucleotides. Exact matching denotes zero mismatches at both ends; for D oligos, any concrete forward variant and any concrete reverse variant may supply that match. A species contributes once even if it has multiple compatible products. Both DNA orientations are scanned. Uncalled bases at a binding site do not count as matches.

In the original target analysis, a product had to span 100–500 bp with at least 50 internal bases. Original background screening allowed 0–1,000 internal bases to search more broadly for potential off-target products. Strict WC and GT models used two mismatches plus the terminal-triplet rule. Relaxed WC allowed three mismatches without the terminal restriction; relaxed GT allowed two. The GT model additionally accepts primer G against template A and primer T against template C in primer-oriented sequences, representing the corresponding wobble pairings on the opposite strand. These accepted pairings are also exempt from the terminal mismatch count in the GT model.

For direct sensitivity–background comparisons, the revision applies identical 100–500 bp and ≥50 internal-base geometry to target and background, WC and GT compatibility, zero to four mismatches per end, and zero, three or five protected terminal bases. Mismatch positions are recorded as distance from the terminal base, with distance zero at the 3′ end. An exploratory positional penalty sums `1 + 4 exp(-d/3)` over mismatches at both ends. This exploratory score summarizes proximity to the terminus. Accession-level scores are supplied in `revision/positional_penalties.csv`.

Background selection limits of ten strict and 50 relaxed hits were heuristic tolerances on the original panels. They do not represent a false-positive rate. H/HD exceeded the relaxed limits and are retained to show why strict matches alone are insufficient for assay choice. No 500-hit truncation was applied: the original values of 500 and 501 are actual record counts under that model and geometry.

### Consensus revisions and factorial controls

Before the XSD replacement, every one of the 912 positions across all 44 mitochondrial oligos, including D oligos, was audited against complete binding motifs. Replacements were tested when an excluded template base was more frequent than every base already represented at that position. Single substitutions and combinations within pairs produced 23 alternatives after the adopted XS G10A change. We tracked gained and lost species, all four original background models, terminal filters and thermodynamics.

The exception allowing a temperature difference up to 6°C for G10A was formulated after inspecting that candidate.  Its criteria were a nonterminal substitution, new-base frequency ≥90%, exact-match gain ≥10 percentage points, retention of every previously covered species and no added background record. We now expose sensitivity to temperature limits of 4, 5, 6 and 7°C and frequency limits of 80%, 85%, 90% and 95%. The threshold grid quantifies the sensitivity of candidate selection.

Three ordinary/D pairs in the preceding candidate set were non-nested. For each, we constructed a 2 × 2 control: ordinary sequence, fixed-base changes alone, degeneracy changes alone, and both changes. The degeneracy-only state nests the ordinary oligos. Controls were evaluated separately against terminal and thermodynamic filters. Fixed endpoints alone isolate location, not the biochemical effect of degeneracy when fixed bases also change.

### Common panels, additional records and genus-blocked reselection

Coverage and exact matching were recomputed on the intersection of evaluable species across all windows of each marker, including H for 16S. Wilson (1927) intervals describe binomial uncertainty conditional on these panels; GenBank species are not a random sample of Actiniaria. Ordinary-only and D-only species counts are reported in `revision/nesting.csv`.

The frozen archives supplied 694 named-species 12S records and 1,116 named-species 16S records. We aligned each expanded pool and reported any-record, strict-majority-record and all-record matching among evaluable records for each species. Complete records with different binding motifs flag observed within-species variability. Expanded alignments can change evaluability, so this analysis combines record and alignment sensitivity rather than estimating a pure replicate-record effect.

Five genus blocks were assigned deterministically by SHA-256 of the genus name modulo five. In each fold, ordinary and D candidates were regenerated using only training genera at the fixed published endpoints, with the same motif-generation and thermodynamic gates. Candidate selection maximized training coverage, exact matching and lower degeneracy. The background panel was not used to tune these diagnostic candidates. Withheld genera were scored by raw paired scanning. This tests sequence reselection conditional on endpoints and an alignment derived from the complete panel.

### Reference discrimination and assay compatibility

Primer-binding regions were removed before barcode comparison. Own-panel reference uniqueness and common-panel uniqueness were calculated separately. For the latter, all windows used the same species with complete unambiguous internal sequences. A joint metric counted species both predicted to amplify and unique within the corresponding reference panel. It describes the sampled reference library, not assured species identification in nature.

Nearest-neighbor differences were calculated as global unit-cost edit distance, including substitutions and indels. We report distance classes zero, one, two and at least three, and whether identical inserts are confined to a genus or family. Leave-one-reference nearest-neighbor assignments are counted correct only when all equally nearest references share the query's taxon. With one record per species, leave-one-reference species classification is not estimable. The expanded multi-record analysis instead reports observed maximum within-species distance, minimum between-species distance, and strict leave-one-record species assignment. A positive observed gap requires the former to be smaller than the latter and at least two records for that species.

Cross-assay heterodimers were screened across concrete variants of all pairs within each marker at 50 nM per variant under the same salt assumptions. Product-length summaries include each distinct accession–product-length combination and report potential ≥20-bp overlap for nominal paired reads of 150, 250 or 300 bases. These estimates precede quality and adapter trimming. Actual adapter sequences remain unspecified, so no adapter-bearing structure result is reported.

### Independent implementation checks

A Python seed-and-extension scanner and an independent C++ exhaustive all-window scanner checked the original frozen-panel outputs. The original 8,940 and 45,780 comparisons combine one target model with two background model families; each background family checks its strict and relaxed flags together. The new equal-geometry comparison checked 57,584 evaluable-record/model combinations at k = 2 and terminal protection of three bases, with no discrepancies. Agreement supports computational reproducibility on the frozen inputs.

### Nuclear 28S references and selection

The 28S search retrieved 11,964 records. Annotated nuclear rRNA and unambiguous single-locus titles were curated, excluding non-binomial names, sequences shorter than 150 nt and more than 2% ambiguous bases. One representative per species maximized concrete length, yielding 221 Actiniaria and 448 other hexacorallian species. MAFFT aligned direction-corrected targets and added backgrounds at fixed alignment length. Coordinates refer to *Halcampoides purpureus* MT676790.1. Anth-28S-inspired 95% position-wise consensus and observed motifs guided a bounded search; each D oligo has at most eight variants, without restricting ambiguity to its 3′ half. Ordinary/D endpoints are identical; expanded variants pass terminal and thermal filters.

Both target and background scans use 100–500 bp products and at least 50 internal bases. Strict matching allows at most two mismatches per end and requires compatible terminal three bases; relaxed WC3/GT2 removes terminal protection. Independent C++ and Python scans agreed on 90,984 checks. Cross-marker comparisons link species names and compare primer-trimmed inserts. Mitochondrial genus-holdout results do not validate 28S. See [28S methods and evidence](28S/README.md).

## Results

### Equal-denominator comparisons change the interpretation of window coverage

The common coverage panels contained 240 12S species and only 67 16S species. The latter is much smaller than the 258–262 species evaluable at the 16S L and XS windows. Percentages calculated on the original window-specific panels therefore should not be ranked as if they represented the same taxa. On the common 16S panel, SD covered 62/67, LD 59/67, and MD and XLD each 58/67. HD covered 65/67 but has substantial relaxed background matching. Table 2 presents the original denominators together with common-panel counts and insert discrimination.

![Common-panel matching](publications/figures/review_common_panel.svg)

**Figure 1.** Covered and exact species on the same evaluable panel within each marker. Numbers above bars are counts. Whiskers are conditional Wilson 95% intervals. H is the manually designed site comparator, placed last. The displayed uncertainty does not account for taxonomic sampling bias or candidate selection.

![Binding-site variation](publications/figures/review_entropy.svg)

**Figure 2.** Reference-position variation, sampling support and primer geometry. Each marker has three panels sharing a gene-relative, 1-based coordinate axis for NC_000933.1: Shannon entropy of called A/C/G/T bases (top), the number of species contributing a called base (middle), and forward/reverse binding intervals (bottom). Gaps and ambiguous calls are excluded; entropy is not plotted where fewer than two species contribute. Low entropy at poorly sampled positions does not establish conservation. Each arrow spans the recorded binding interval and points toward its 3′ end; connecting lines show the intervening region. Ordinary and D primers share endpoints and occupy one row. The 12S and 16S panels have separate coordinate scales. Shared sites include the 12S XS/S/L forward region and M reverse region, and the 16S M/XL forward region.

### Factorial controls support replacement of 12S XSD

Following replacement of 12S XSD, nine of eleven D pairs include their corresponding ordinary sequences. The remaining exceptions are 12S S/SD and 16S XS/XSD. The factorial audit retains the previous 12S XSD for comparison. The nested 12S XS control passed both terminal and thermal filters: maximum F/R difference was 4.84°C, minimum variant Tm 58.99°C, and maximum secondary-structure Tm 49.74°C. For 12S XS, adding only the degenerate positions raised coverage from 205 to 227/245 and exact matches from 192 to 212. Changing the fixed bases alone preserved 205 covered species but reduced exact matching to zero. Combining both changes gave the previous XSD result, 230 covered and zero exact. We replaced that pair with the nested control as the current XSD: coverage decreased by three species, exact matching increased by 212 species, and all nine strict background hits were retained. The 16S nested XS control similarly covered 223/262 with 167 exact matches, compared with 223 covered and one exact for released XSD. Its terminal filter fails, so this numerical improvement is not sufficient for adoption. The full four-state decomposition is given in Table 3.

For nested pairs, 16S LD increased original-panel coverage from 138 to 231/258 and exact matching from 106 to 180. That change is attributable to adding compatible variants at fixed positions, subject to concentration and structure assumptions. It also increases relaxed background matching. Neither result supports a general claim that higher degeneracy is always better.

### Consensus substitutions are sensitive to the selection rule

At XS forward position 10, 228 of 230 complete 12S motifs contained A and two contained G. Fifteen additional evaluable species were absent from this complete-motif calculation: twelve had projected binding-site gaps and three had other binding-segment length changes. They remained eligible for raw-sequence scanning. G10A increased paired exact matches from 2 to 192/245, with 192 gains and two losses of exact matches. All 205 covered species and all four background hit sets were retained. The predicted F/R Tm differences were 5.18°C at 50 nM and 5.33°C at 200 nM per ordinary oligo. Thus the original 5°C gate would reject the change; the post hoc 6°C exception admitted it.

For 12S S, forward C9A retained 173 covered species and increased exact matching from zero to four, with maximum temperature difference 5.95°C. It did not meet the post hoc ten-percentage-point exact-gain criterion. Combining it with reverse T16C produced 113 exact matches but lost eight previously covered species, despite the smaller net coverage decrease. The full 23-candidate audit, including temperature and background changes, is included in Table S1 and the supplementary sensitivity grid.

At 12S L, forward T1C increased exact matching from one to 202 while two *Lebrunia* species no longer met the two-mismatch rule. Changes to the previous XSD additionally affected *Phyllodiscus semoni*. All three are Aliciidae in the archived taxonomy. These are threshold-dependent predictions, not evidence that PCR would necessarily exclude the family. The distance-to-3′ profiles distinguish these 5′ changes from terminal incompatibilities.

![Mismatch positions](publications/figures/review_mismatch_profiles.svg)

**Figure 3.** Mismatch frequency at called homologous positions, ordered by distance from the 3′ terminus. The line separates the protected terminal triplet from other positions. Gaps and ambiguous calls are excluded from each position's denominator. White cells indicate positions outside an oligo or without called observations. These aligned-site diagnostics complement, rather than replace, raw paired-product scanning.

### The background cost persists under a common product geometry

Using identical product geometry for target and background retained the nine strict 12S hits at XS, XSD, M, MD, L and LD. They are the same nine Antipatharia accessions in the original background panel. Geometry reduced relaxed WC hits for 16S H and HD from 486 and 976 to 452 and 939, respectively, and for LD from 46 to 45. H/HD therefore remain unfavorable background comparators even after length is controlled. In the original LD relaxed-WC set, most hits were Antipatharia, linking the trade-off to a specific co-occurring anthozoan group.

![Common-model trade-offs](publications/figures/review_unified_models.svg)

**Figure 4.** Target coverage against background record counts under identical product geometry, separately for WC and GT. The horizontal axis uses a symmetric-log scale with a linear region through one record. Each curve connects zero through four allowed mismatches per end with the terminal three bases protected. Solid lines denote ordinary pairs and dashed lines D pairs. Target percentages use each window's evaluable denominator; these are model-sensitivity curves, not experimentally measured ROC curves.

![All background groups](publications/figures/review_background_all.svg)

**Figure 5.** Original-panel background hits for all 22 pairs. Columns within each taxonomic group are strict WC, strict GT, relaxed WC and relaxed GT; nonzero cells display record counts. Both panels share one count color scale (symlog, linear from zero to one). W and G denote strict WC and GT; W+ and G+ denote their relaxed models. Marker-specific panel sizes and original background geometry apply; cross-marker raw counts are not directly comparable.

### Additional records and withheld genera expose different uncertainties

Across the 694 named 12S records, XS covered at least one evaluable record in 208 species, a strict majority in 206, and every evaluable record in 203; 247 species were evaluable in the expanded alignment. For 16S LD, the corresponding counts were 233, 232 and 225 of 259. Binding motifs varied within 10 and 23 species, respectively. An assay's species count therefore depends on whether one representative, any available record or every available record is required to match.

The holdout analysis included 230 species with complete, equal-length binding motifs at both ends, a subset of the 245 raw-scan-evaluable species. Genus-blocked fixed-window reselection yielded 205/230 covered but zero exactly matched held-out motifs for ordinary 12S XS under the original 5°C gate. D reselection yielded 226 covered and 212 exact. The released G10A ordinary primer cannot be treated as validated by this diagnostic because it uses a different post hoc temperature gate. For 12S M, ordinary reselection covered 219/226 with 197 exact; D reselection covered 219 with 205 exact. For 16S L, ordinary reselection covered 137/246 with 106 exact; D reselection covered 213 with 186 exact. These candidates were regenerated within each training fold and are not the released primer pair tested repeatedly. All fold-specific sequences, denominators and unavailable-candidate flags are provided.

### Reference discrimination depends on location, sampling and amplification

In the original one-record panels, unique 12S inserts occurred in 35/244 species for XS, 44/244 for S, 86/239 for M, 68/243 for L and 100/239 for XL. The decrease from M to L contradicts a monotonic length explanation. M and XL include a different variable region from L. The reference entropy profile provides positional context for this pattern.

The common complete-insert analysis and joint covered-and-unique counts are shown in Figure 6. This metric excludes species whose barcode is unique but whose primer pair is not predicted to amplify. At 16S, the original H and L inserts had different proportions unique (87/159 and 86/255), despite nearly identical numerators; the window-specific species panels explain why the denominators differ.

![Joint coverage and discrimination](publications/figures/review_joint_resolution.svg)

**Figure 6.** Percentage of species both covered and unique in the same complete-insert species panel across all windows of a marker. Numbers above bars are species counts. This is sampled-reference uniqueness, not a validated assignment success rate.

In the multi-record analysis, a positive observed barcode gap was found for only 11 of 111 multi-record 12S XS species, compared with 28 of 110 at XL. Strict leave-one-record assignment correctly recovered 52 of 658 XS records and 179 of 646 XL records; ties spanning species were not accepted as correct species identifications. At 16S L, 18 of 130 multi-record species had a positive observed gap. Sparse and uneven within-species sampling limits extrapolation from these values. Nearest-neighbor edit-distance distributions and genus/family diagnostics are supplied at accession level.

### Assays should be evaluated separately before any multiplex combination

The cross-assay screen identified strong complementarities between the 12S XS forward and M/MD reverse oligos, and between SD forward and M/MD reverse oligos. Their maximum heterodimer Tm exceeded the 55°C screening threshold; S forward with an MD reverse variant also exceeded it. The shared 574-region explains these interactions. They should not be combined as a multiplex on the basis of the present single-pair design. Other combinations, including 16S oligos with shorter terminal complementarities, are tabulated but have not been experimentally cleared for multiplexing. Actual product lengths are reported in Table 2; nominal read-overlap counts are supplied in `revision/length_summary.csv`. Adapter-bearing structures require the final adapter sequences.

### Nuclear 28S complementary discrimination

The nuclear 28S panel adds 221 sea-anemone species and 448 other hexacorallian species. Eight ordinary/D pairs share coordinates within XS, S, M and L; degeneracy is at most eight per oligo. The 95% position-wise consensus guided a finite search, not a claim of 95% paired coverage. Strict background hits are zero, but relaxed hits remain high. Of 112 species pairs sharing all comparable 12S and 16S inserts, 67 have evaluable 28S L inserts; 64 differ, and LD covers both species in 63 pairs. Names link reference species across markers, not individual specimens.

## Discussion

The comparison identifies three complementary priorities: 12S XSD for short fragments, 12S MD for compatibility and reference discrimination, and 16S LD for broader target coverage. Their trade-offs follow from binding-site variation and the information carried by the intervening sequence. G10A improves exact compatibility with many complete motifs, but the unchanged two-mismatch coverage and its internal position prevent translating that gain into a measured detection improvement. Similarly, a 5′ substitution can move a species across an arbitrary mismatch boundary without establishing that amplification is biologically impossible.

Fixed coordinates support a controlled comparison of location. They isolate degeneracy only when the ordinary oligos are included among the D variants. The factorial controls resolve the previously misleading XSD contrasts and reveal candidate sequences worth evaluating under a prospective set of constraints. They also show why matching counts cannot be the sole acceptance criterion: the promising nested 16S XS control violates the existing terminal filter. The next selection round should declare its thermodynamic and terminal tolerances before scoring withheld taxa, retain consensus candidates through paired screening, and report a Pareto set instead of treating a single heuristic winner as an optimum.

The reference library is another part of the assay. Common panels remove one source of confounding but reduce the available 16S evidence to 67 species. Multi-record analysis exposes within-species variability, while conditional genus-blocked reselection addresses dependence on the sampled genera without validating window discovery. Broader taxonomic sampling and independent window selection would extend the test to additional sea-anemone lineages.

Slow anthozoan mitochondrial evolution provides biological context for shared inserts (Shearer et al., 2002; Hellberg, 2006). Neither greater amplicon length nor higher coverage guarantees species assignment. A short assay can support detection or broader taxonomic assignment while requiring a complementary locus for closely related species. The observed barcode-gap and nearest-neighbor results justify this limitation directly rather than inferring it from primer matching.

### Priorities for assay evaluation

| Intended use | Computational priority | Reason and qualification |
| --- | --- | --- |
| Short-fragment sea-anemone detection | 12S XSD | Passes terminal and thermal filters; 227 covered and 212 exact, compared with 205/192 for ordinary XS and 230/0 for the previous XSD. Fish and crustacean co-amplification was not evaluated. |
| Matching plus greater reference discrimination | 12S MD | Strong exact matching and more unique inserts than XS; black-coral backgrounds remain relevant. |
| Broad 16S target compatibility | 16S LD | Large coverage gain over L; relaxed background hits include Antipatharia. |
| Habitats rich in black corals | Compare 12S S/SD with 16S alternatives on local references | Original background screening favors some windows, but target completeness and regional species composition constrain interpretation. |
| Historical hand-designed site | H/HD as comparators | High relaxed background matching argues against prioritizing them for mixed-community assays. |
| Multiplex reactions | No released combination endorsed | Cross-assay complementarity and final adapter sequences require separate assessment. |

The background panels comprise other cnidarians. Co-amplification of fishes, crustaceans and other non-cnidarian organisms was not evaluated, an important consideration for the 12S candidates. This computational study evaluates binding compatibility and reference discrimination. Environmental mixtures introduce template competition and reaction-dependent amplification, which require assay-level evaluation.


## Primers and base-matching diagrams

All 30 pairs are organized by 12S, 16S and 28S. Sequences run 5′→3′; D denotes a degenerate version. Binding maps use gene-specific reference coordinates; arrows indicate orientation and spans include primers. Coverage/exact counts use window-specific evaluable species and are not equal-denominator rankings.

[FASTA](primers.fasta) · [CSV](all_results.csv)

### 12S

![12S binding positions](figures/12S_binding.svg)

| Pair | F (5′–3′) | R (5′–3′) | Degeneracy F/R | Covered | Exact | Base matching |
|---|---|---|---|---|---|---|
| XS | CGAAAGGTAAAATCTTCCGCCT | CGACGGCCATGCAATACCGG | 1/1 | 205/245 | 192/245 | [F / R](figures/template_matching/12S_XS.svg) |
| XSD | CGAAAGGTAAAATCTTCCGCCT | CGACGGCCATGCARTACCDG | 1/6 | 227/245 | 212/245 | [F / R](figures/template_matching/12S_XSD.svg) |
| S | CGAAAGGTCAAATCTTCCGC | CGGGTTGGGTCCCTTTTAGC | 1/1 | 173/246 | 0/246 | [F / R](figures/template_matching/12S_S.svg) |
| SD | CGAAAGGTAAAATCTTCCGC | CGGATTGGGTCCCTTYTAGY | 1/4 | 195/246 | 26/246 | [F / R](figures/template_matching/12S_SD.svg) |
| M | AAGGCGTAAAGGATGTGTAG | GCGGAAGATTTTACCTTTCG | 1/1 | 220/242 | 198/242 | [F / R](figures/template_matching/12S_M.svg) |
| MD | AAGGCGTAAAGGATGTGTVG | GCGGAAGATTTDACCTTTCR | 3/6 | 224/242 | 208/242 | [F / R](figures/template_matching/12S_MD.svg) |
| L | TGAAAGGTAAAATCTTCCGC | TATAAGAAAATGTAGCGCAT | 1/1 | 226/246 | 1/246 | [F / R](figures/template_matching/12S_L.svg) |
| LD | TGAAAGGTAAAATCTTCCGC | TATAAGAAAATGTAGCGCAY | 1/2 | 229/246 | 1/246 | [F / R](figures/template_matching/12S_LD.svg) |
| XL | GCGTAAAGGATGTGTAGGCC | CTTAACACACGGACTTACGACG | 1/1 | 185/243 | 46/243 | [F / R](figures/template_matching/12S_XL.svg) |
| XLD | GCGTAAAGGATGTGTARGYY | CTTAACACACGRACTTACGACG | 8/2 | 225/243 | 52/243 | [F / R](figures/template_matching/12S_XLD.svg) |

<details>
<summary>12S MD — base-matching preview</summary>

![12S MD](figures/template_matching/12S_MD_preview.svg)

[Full species](figures/template_matching/12S_MD.svg)

</details>

### 16S

![16S binding positions](figures/16S_binding.svg)

| Pair | F (5′–3′) | R (5′–3′) | Degeneracy F/R | Covered | Exact | Base matching |
|---|---|---|---|---|---|---|
| H | AAGACCCCATTGAGCTTTACTA | CCCCATGGTAGCTTTTATTCGT | 1/1 | 146/165 | 122/165 | [F / R](figures/template_matching/16S_H.svg) |
| HD | AAGACCCCATTGAGCTTTACYW | CCCCATGGTAGCTTTTATYCRT | 4/4 | 154/165 | 136/165 | [F / R](figures/template_matching/16S_HD.svg) |
| XS | AGTTTAGATAATGTGGGATC | TAACTACGGACCTTAAGCTC | 1/1 | 165/262 | 132/262 | [F / R](figures/template_matching/16S_XS.svg) |
| XSD | AGTTTAGATAATGTGGARBC | TAATTACGGACCTTAAGCTC | 6/1 | 223/262 | 1/262 | [F / R](figures/template_matching/16S_XSD.svg) |
| S | CGTAGTTAATAGCAAGTGGCCT | AGGTCGCAAACATCGTTGTCGA | 1/1 | 63/77 | 0/77 | [F / R](figures/template_matching/16S_S.svg) |
| SD | CGTAGTTAATAGCHRGTGGCCT | AGGTCGCAAACATCGKYGTCRA | 6/8 | 70/77 | 56/77 | [F / R](figures/template_matching/16S_SD.svg) |
| M | AAGCTAAATGCTGAAGTATG | TTCGTCTTACAATATGGACG | 1/1 | 61/76 | 51/76 | [F / R](figures/template_matching/16S_M.svg) |
| MD | AAGCTAAATGCTGAAGTATG | TTCGTCTTACAATWTGDACG | 1/6 | 61/76 | 52/76 | [F / R](figures/template_matching/16S_MD.svg) |
| L | GGCCGGTATGAATGGCGTCACG | CCCCTTGCAGTCAGGCCACCGG | 1/1 | 138/258 | 106/258 | [F / R](figures/template_matching/16S_L.svg) |
| LD | GGCCGGTATGAATGGCRYYACG | CCCCTTGCAGTCAGGCCACYDG | 8/6 | 231/258 | 180/258 | [F / R](figures/template_matching/16S_LD.svg) |
| XL | AAGCTAAATGCTGAAGTATG | AACCAAACTATCCAACTTAC | 1/1 | 41/76 | 36/76 | [F / R](figures/template_matching/16S_XL.svg) |
| XLD | AAGCTAAATGCTGAAGTATG | AACCAAACTATCYAACTYRC | 1/8 | 61/76 | 52/76 | [F / R](figures/template_matching/16S_XLD.svg) |

<details>
<summary>16S LD — base-matching preview</summary>

![16S LD](figures/template_matching/16S_LD_preview.svg)

[Full species](figures/template_matching/16S_LD.svg)

</details>

### 28S

![28S binding positions](figures/28S_binding.svg)

| Pair | F (5′–3′) | R (5′–3′) | Degeneracy F/R | Covered | Exact | Base matching |
|---|---|---|---|---|---|---|
| XS | GGGTTGTTTGGGAATGCA | TTGACTCTCTTTTCAAAGTG | 1/1 | 179/182 | 7/182 | [F / R](figures/template_matching/28S_XS.svg) |
| XSD | GGGTTGTTYKGGAATGCA | TKRACTCTCTTTTCAAAGTG | 4/4 | 179/182 | 177/182 | [F / R](figures/template_matching/28S_XSD.svg) |
| S | CGTCATAGAGGGTGACAA | TTAACTCTCTTTTCAAAGTG | 1/1 | 163/182 | 134/182 | [F / R](figures/template_matching/28S_S.svg) |
| SD | CGTCRTAGAGGGTGASAA | TTRACTCTCTTTTCAAAGTG | 4/2 | 172/182 | 148/182 | [F / R](figures/template_matching/28S_SD.svg) |
| M | AAGTTGCTTGGAACAGCA | TTTTAACTCTCTTTTCAAAGTG | 1/1 | 169/180 | 139/180 | [F / R](figures/template_matching/28S_M.svg) |
| MD | AAGTTGCTTGGAACAGCA | KTTTAACTCTCTTTTCAAAGTG | 1/2 | 169/180 | 154/180 | [F / R](figures/template_matching/28S_MD.svg) |
| L | TGAGGGAAAGGTGAAAAGCA | ACTGCAAGCAGTCCACCG | 1/1 | 94/172 | 0/172 | [F / R](figures/template_matching/28S_L.svg) |
| LD | TGAGGGAAAGRTGAAAAGCA | ACTGCAHGCAGTCCACCR | 2/6 | 166/172 | 158/172 | [F / R](figures/template_matching/28S_LD.svg) |

<details>
<summary>28S XS — base-matching preview</summary>

![28S XS](figures/template_matching/28S_XS_preview.svg)

[Full species](figures/template_matching/28S_XS.svg)

</details>

<details>
<summary>28S XSD — base-matching preview</summary>

![28S XSD](figures/template_matching/28S_XSD_preview.svg)

[Full species](figures/template_matching/28S_XSD.svg)

</details>

<details>
<summary>28S S — base-matching preview</summary>

![28S S](figures/template_matching/28S_S_preview.svg)

[Full species](figures/template_matching/28S_S.svg)

</details>

<details>
<summary>28S SD — base-matching preview</summary>

![28S SD](figures/template_matching/28S_SD_preview.svg)

[Full species](figures/template_matching/28S_SD.svg)

</details>

<details>
<summary>28S M — base-matching preview</summary>

![28S M](figures/template_matching/28S_M_preview.svg)

[Full species](figures/template_matching/28S_M.svg)

</details>

<details>
<summary>28S MD — base-matching preview</summary>

![28S MD](figures/template_matching/28S_MD_preview.svg)

[Full species](figures/template_matching/28S_MD.svg)

</details>

<details>
<summary>28S L — base-matching preview</summary>

![28S L](figures/template_matching/28S_L_preview.svg)

[Full species](figures/template_matching/28S_L.svg)

</details>

<details>
<summary>28S LD — base-matching preview</summary>

![28S LD](figures/template_matching/28S_LD_preview.svg)

[Full species](figures/template_matching/28S_LD.svg)

</details>

## Data and reproducibility

[12S/16S analysis](revision/README.md) · [28S analysis](28S/README.md) · [Preprint and build](publications/README.md) · [Supplement](publications/MiAnemone_supplement.md)

## Citation

Xia, F. (2026). *MiAnemone: in silico design and evaluation of mitochondrial and nuclear rRNA primers for sea anemone (Actiniaria) eDNA metabarcoding*. Computational preprint v1.5. [GitHub Release](https://github.com/GUIBA-EX/MiAnemone/releases/tag/preprint-v1.5).

[Copyright and permissions](COPYRIGHT.md)

## References

Bae S. 2025. Validation of improved cytochrome c oxidase I (COI) primers for comprehensive biodiversity assessment of ascidians. *PeerJ* 13:e19671. https://doi.org/10.7717/peerj.19671

Hellberg ME. 2006. No variation and low synonymous substitution rates in coral mtDNA despite high nuclear variation. *BMC Evolutionary Biology* 6:24. https://doi.org/10.1186/1471-2148-6-24

Katoh K, Standley DM. 2013. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30:772–780. https://doi.org/10.1093/molbev/mst010

McCartin L et al. 2024. Nuclear eDNA metabarcoding primers for anthozoan coral biodiversity assessment. *PeerJ* 12:e18607. https://doi.org/10.7717/peerj.18607

Miya M et al. 2015. MiFish, a set of universal PCR primers for metabarcoding environmental DNA from fishes: detection of more than 230 subtropical marine species. *Royal Society Open Science* 2:150088. https://doi.org/10.1098/rsos.150088

Piñol J, Mir G, Gomez-Polo P, Agustí N. 2015. Universal and blocking primer mismatches limit the use of high-throughput DNA sequencing for the quantitative metabarcoding of arthropods. *Molecular Ecology Resources* 15:819–830. https://doi.org/10.1111/1755-0998.12355

Rodríguez E et al. 2014. Hidden among sea anemones: the first comprehensive phylogenetic reconstruction of the order Actiniaria (Cnidaria, Anthozoa, Hexacorallia) reveals a novel group of hexacorals. *PLOS ONE* 9:e96998. https://doi.org/10.1371/journal.pone.0096998

SantaLucia J Jr. 1998. A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *Proceedings of the National Academy of Sciences USA* 95:1460–1465. https://doi.org/10.1073/pnas.95.4.1460

Shearer TL, Van Oppen MJH, Romano SL, Wörheide G. 2002. Slow mitochondrial DNA sequence evolution in the Anthozoa (Cnidaria). *Molecular Ecology* 11:2475–2487. https://doi.org/10.1046/j.1365-294X.2002.01652.x

Untergasser A et al. 2012. Primer3—new capabilities and interfaces. *Nucleic Acids Research* 40:e115. https://doi.org/10.1093/nar/gks596

Wilson EB. 1927. Probable inference, the law of succession, and statistical inference. *Journal of the American Statistical Association* 22:209–212. https://doi.org/10.2307/2276774
