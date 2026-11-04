# MiAnemone preprint

**Fei Xia** · Computational preprint; not peer reviewed.

[PDF](MiAnemone_preprint.pdf) · [Manuscript](MiAnemone_preprint.md) · [Versioned releases](https://github.com/GUIBA-EX/MiAnemone/releases) · [Citation](../CITATION.cff)

## Build and validate

Install Python 3, Pandoc, Poppler and Google Chrome or Chromium. Set `CHROME` to the browser executable when needed. Run from the repository root:

```sh
python3 -m pip install -r publications/requirements-preprint.txt biopython
python3 publications/build_preprint.py
python3 scripts/verify.py
pdftotext -layout publications/MiAnemone_preprint.pdf /tmp/paper.txt
python3 scripts/check_preprint.py /tmp/paper.txt
python3 scripts/check_publication.py
python3 scripts/preprint_integrity.py
```

## Automatic updates

Pushes to `main` and manual workflow runs check the committed research results and publication files. Changes to README editions, citations, data, figures, publication sources or release scripts trigger a new minor preprint version. The workflow stamps all editions using Asia/Shanghai time, rebuilds the PDF, checks scientific tables and primer sequences, validates publication metadata and links, and verifies checksums before creating a tagged GitHub Release. A run with no relevant changes validates the existing version without creating a new one.

The release contains the PDF, matching primer FASTA and SHA-256 checksums. Existing public releases remain unchanged. Relative links in the PDF point to the matching source tag. Scientific analyses and their interpretations are maintained in the source files; the workflow rebuilds the paper from those committed results.

## Consistency checks

`release.json` holds the title, version, date, README editions and primer FASTA path. `scripts/stamp_preprint.py` synchronizes the manuscript and citation formats. `scripts/check_publication.py` verifies metadata, local links, visible figure captions and every oligo in the PDF and manuscript. Project-specific checks separately verify scientific tables against their source CSV files. `publication_check.json` records the checked version, page and oligo counts, and output hashes.
