# MiAnemone preprint

The current version is recorded in [release.json](release.json). The manuscript, bilingual READMEs, sequences and seven figures describe the same 30 pairs. The 12S XSD replacement is documented in Table 3; prior sequences are confined to the factorial audit.

## Build

From the repository root, install the Python dependencies and make Pandoc and Chrome available. Figure and PDF rendering require an installed Avenir Next font; the published PDF embeds it.

```sh
python publications/plot_revision.py
python publications/unify_figure_style.py
python scripts/figures_28S.py
python scripts/stamp_preprint.py 1.5
python publications/build_preprint.py
python scripts/verify.py
python scripts/check_figures.py
python scripts/check_28S.py
pdftotext -layout publications/MiAnemone_preprint.pdf /tmp/paper.txt
python scripts/check_preprint.py /tmp/paper.txt
python scripts/check_publication.py
python scripts/preprint_integrity.py
```

Computational inputs, methods and rerun commands are indexed in [the reproduction guide](../revision/README.md). [Table S1](MiAnemone_supplement.md) reports the preceding candidate audit, with its relationship to the adopted XSD made explicit. Checksums cover the PDF and primer FASTA release assets. GitHub Actions validates the committed artifacts without substituting fonts or rebuilding the PDF. Release tags freeze the complete source and input snapshot.

Author: Fei Xia. Copyright and permissions follow [COPYRIGHT.md](../COPYRIGHT.md).

The current manuscript and READMEs include the nuclear 28S extension. The older assembly scripts regenerate the mitochondrial revision only; edit the integrated documents directly when revising this release. Nuclear inputs and reproduction steps are in [28S](../28S/README.md).
