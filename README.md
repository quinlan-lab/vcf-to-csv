# Annotated VCF to CSV

This utility reads an annotated biallelic VCF or BCF with `cyvcf2` and writes
one CSV row per variant. Multiallelic records are rejected and should
be split before export.

Each row includes:

- `chrom`, `pos`, `ref`, `alt`, `alt_index`, `rsid`, `qual`, and `filter`
- strict call-rate counts (`call_rate`, `n_called`, `n_missing`)
- allele-specific `het_samples` and `homalt_samples` lists and counts
- HET and HOMALT counts within samples of interest and all remaining samples
- a JSON `vep_annotations` cell containing selected transcript annotations
- any additional VEP or top-level INFO fields defined in TOML

Sample lists are comma-separated, missing values are empty, and VCF floating
point values are written with seven significant digits.

## Install

Python 3.9 or newer is required.

```bash
python3 -m venv .venv
.venv/bin/pip install -e '.[dev]'
```

## Run

```bash
.venv/bin/vcf-to-csv \
  input.annotated.vcf.gz \
  variants.csv \
  --config config.example.toml \
  --samples-of-interest samples.txt
```

`python -m vcf_to_csv` is also available after installation.

The command refuses to replace an existing CSV unless `--force` is supplied.
Use `-` as the output path to write CSV to standard output.

### Run the public test fixture

The repository includes a four-sample, two-variant subset of the public 1000
Genomes Phase 3 callset. Run it directly and write the CSV to standard output:

```bash
.venv/bin/vcf-to-csv \
  tests/data/1000genomes-phase3-subset.vcf \
  - \
  --config tests/data/1000genomes.toml
```

Selected columns from the two output rows are:

| chrom | pos | rsid | global_af | het_samples | homalt_samples |
| --- | ---: | --- | ---: | --- | --- |
| 20 | 10000117 | rs4816203 | 0.605431 | HG00096,HG00097,NA12878 | |
| 20 | 10000598 | rs6057087 | 0.810503 | HG00096,HG00097,HG00099 | NA12878 |

### Samples of interest

The optional `--samples-of-interest` file is plain text with one VCF sample ID
per line and no header. For example:

```text
# Case samples
HG00096
NA12878
```

IDs must exactly match the VCF `#CHROM` header. Blank lines and `#` comments are
ignored; unknown, duplicate, or whitespace-containing IDs are errors. The
`interest_*` columns count carriers in this file, while `other_*` columns count
all remaining carriers. Without the option, all carriers are counted as other.

## Configure annotations

Copy `config.example.toml` and update it to match the VCF header. A custom
column can come from either:

- `source = "vep"`: a named subfield inside the configured `CSQ` field
- `source = "info"`: a top-level VCF INFO field

Supported aggregation modes are:

- `first`: first populated value
- `join`: distinct values in input order
- `max` or `min`: numeric aggregation, useful for population frequencies

For example:

```toml
[[columns]]
name = "gnomad_exomes_af"
source = "vep"
field = "gnomADe_AF"
aggregate = "max"
required = false
```

Set `required = true` when a missing header field should stop the export. This
helps catch spelling mistakes and annotation-pipeline drift.

## Genotype and multiallelic rules

Partially missing calls such as `0/.` and `./1` are treated as missing. Call
rate is the number of fully called samples divided by the number of VCF
samples.

Multiallelic records are rejected with an error. Split or normalize them into
biallelic records before running this exporter.

The `vep_annotations` column is JSON rather than a lossy display string. This
preserves multiple transcripts for a future web interface while still keeping
the current output to a single CSV column.

## Test

```bash
.venv/bin/pytest
.venv/bin/ruff check .
```
