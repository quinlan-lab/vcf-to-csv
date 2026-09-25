[![Test](https://github.com/quinlan-lab/vcf-to-csv/actions/workflows/test.yml/badge.svg)](https://github.com/quinlan-lab/vcf-to-csv/actions/workflows/test.yml)

# Annotated VCF to CSV

This utility reads an annotated biallelic VCF or BCF with `cyvcf2` and writes
one CSV row per variant. Multiallelic records are rejected and should
be split before export.

Each row includes:

- `chrom`, `pos`, `ref`, `alt`, `alt_index`, `rsid`, `qual`, and `filter`
- strict call-rate counts (`call_rate`, `n_called`, `n_missing`)
- allele-specific `het_samples` and `homalt_samples` lists and counts
- HET and HOMALT counts within samples of interest and all remaining samples
- flat, delimiter-separated VEP columns defined in TOML
- any additional top-level INFO fields defined in TOML

Sample lists are comma-separated, missing values are empty, and VCF floating
point values are written with seven significant digits.

## Install

Python 3.9 or newer is required.

```bash
python3 -m venv .venv
.venv/bin/pip install -e '.[dev]'
```

On a FIPS-enabled cluster, the prebuilt `cyvcf2` wheel may fail during import
with `FATAL FIPS SELFTEST FAILURE` because it bundles OpenSSL through libcurl.
Rebuild `cyvcf2` without remote-URL support:

```bash
CYVCF2_HTSLIB_CONFIGURE_OPTIONS=--disable-libcurl \
uv pip install \
  --python .venv/bin/python \
  --reinstall \
  --no-binary cyvcf2 \
  'cyvcf2>=0.34,<1'
```

## Run

```bash
.venv/bin/vcf-to-csv \
  input.annotated.vcf.gz \
  variants.csv \
  --config config.example.toml \
  --samples-of-interest samples.txt \
  --genes genes.txt
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
  --config tests/data/1000genomes.toml \
  --genes tests/data/genes.txt
```

Selected columns from the two output rows are:

| chrom | pos | rsid | global_af | genes | consequence | het_samples | homalt_samples |
| --- | ---: | --- | ---: | --- | --- | --- | --- |
| 20 | 10000117 | rs4816203 | 0.605431 | ANKEF1;SNAP25-AS1 | intron_variant&non_coding_transcript_variant;downstream_gene_variant | HG00096,HG00097,NA12878 | |
| 20 | 10000598 | rs6057087 | 0.810503 | ANKEF1;SNAP25-AS1 | intron_variant&non_coding_transcript_variant;downstream_gene_variant | HG00096,HG00097,HG00099 | NA12878 |

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

### Gene filtering

Use `--genes FILE` (or `--gene-list FILE`) to emit only variants having at
least one matching VEP `CSQ` annotation. The file contains one gene per line,
with the same blank-line and `#` comment handling as the samples-of-interest
file. Values are matched exactly and case-sensitively against both the VEP
`SYMBOL` and `Gene` subfields, so either gene symbols or Ensembl gene IDs can
be used. An empty gene file emits only the CSV header.

## Configure annotations

Copy `config.example.toml` and update it to match the VCF header. A custom
column can come from either:

- `source = "vep"`: a named subfield inside the configured `CSQ` field
- `source = "info"`: a top-level VCF INFO field

Supported aggregation modes are:

- `first`: first populated value
- `join`: values in input order (distinct values for top-level INFO fields)
- `max` or `min`: numeric aggregation, useful for population frequencies

For CSV-friendly VEP output, define one column per useful subfield. Repeated
and empty values are preserved in annotation order, so columns such as symbol,
consequence, and transcript remain positionally aligned:

```toml
[[columns]]
name = "genes"
source = "vep"
field = "SYMBOL"
aggregate = "join"
separator = ";"

[[columns]]
name = "consequences"
source = "vep"
field = "Consequence"
aggregate = "join"
separator = ";"
```

For example:

```toml
[[columns]]
name = "gnomad_exomes_af"
source = "vep"
field = "gnomADe_AF"
aggregate = "max"
required = false
```

Structured INFO fields with a `Format:` declaration can expose one named
subfield as a column. For example, given
`Format: ALLELE|ALL_AF|ALL_AC|ALL_AN`, this extracts only `ALL_AF`:

```toml
[[columns]]
name = "gnomad_all_af"
source = "info"
field = "FV_GNOMAD"
subfield = "ALL_AF"
aggregate = "first"
```

Subfields are pipe-delimited by default. Set `subfield_sep` when a header uses
a different delimiter:

```toml
subfield_sep = "~"
```

Set `required = true` when a missing header field should stop the export. This
helps catch spelling mistakes and annotation-pipeline drift.

## Genotype and multiallelic rules

Partially missing calls such as `0/.` and `./1` are treated as missing. Call
rate is the number of fully called samples divided by the number of VCF
samples.

Multiallelic records are rejected with an error. Split or normalize them into
biallelic records before running this exporter.

The lossless transcript-level JSON is optional. To include it in addition to
the flat columns, set `output_column` and optionally select its fields:

```toml
[vep]
info_field = "CSQ"
output_column = "vep_annotations"
include_fields = ["SYMBOL", "Consequence", "Feature", "HGVSc", "HGVSp"]
```

## Test

```bash
.venv/bin/pytest
.venv/bin/ruff check .
```
