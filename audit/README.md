# Consensus-site review · MiAnemone

All 22 pairs were reviewed on the released target and background panels. [Observed bases](position_support.csv) lists A, C, G and T counts for all 912 oligo positions, in synthesis orientation, among window-evaluable species with complete, ungapped sites on the indicated end. [Candidate comparisons](consensus_trials.csv) contains the 23 additional trials after correcting 12S XS. Sites were flagged when the most frequent observed nucleotide was excluded by the primer and occurred more often than any included nucleotide. We tested the single replacements and all combinations of flagged replacements within each pair. This review tests consensus substitutions. It does not search new windows or exhaustively add new degenerate alleles.

The corrected 12S XS forward sequence is `CGAAAGGTAAAATCTTCCGCCT`. G10A matches 228 rather than 2 of 230 complete forward sites. Paired exact matches increased from 2 to 192 among 245 evaluable species, retaining the same 205 covered species and the same nine background accessions in each model. Predicted maximum paired ΔTm was 5.33°C. An exhaustive scanner independently agreed with the seed-based scanner across all 894 target and background record–model comparisons for this pair.

Nine further pairs contained at least one flagged position. Below are the largest exact-match gains with a loss of at most three covered species. Each comparison uses the same denominator and no new background records were hit in any of the four models.

| Pair | Change | Covered, released → trial | Exact, released → trial | Species losing predicted coverage |
|---|---|---:|---:|---|
| 12S XSD | F1 T→C and R1 A→C | 230 → 227 | 0 → 212 | *Phyllodiscus semoni*, *Lebrunia neglecta*, *L. coralligens* |
| 12S L | F1 T→C | 226 → 224 | 1 → 202 | *L. neglecta*, *L. coralligens* |
| 12S LD | F1 T→C | 229 → 227 | 1 → 204 | *L. neglecta*, *L. coralligens* |
| 12S SD | R4 A→G | 195 → 194 | 26 → 117 | Two previously covered species lost and one gained (see accession-level changes below). |

These are candidates for a different coverage trade-off. The 16S XSD F17 A→G trial would increase coverage from 223 to 224 and exact matches from one to 16 without a lost species. Its degenerate oligo can encode a prohibited five-base homopolymer, so it fails the existing sequence filter. H and HD had no site meeting the flagging criterion. Among the 23 additional trials, none met the full acceptance rules while retaining every previously covered species.

[Accession-level coverage changes](coverage_changes.csv) records gains and losses for every tested alternative.
