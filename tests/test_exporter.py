import csv
from dataclasses import replace
from pathlib import Path

import pytest

from vcf_to_csv.config import ColumnConfig, config_from_mapping, load_config
from vcf_to_csv.exporter import _format_call_rate, export_vcf

ROOT = Path(__file__).parents[1]


def test_call_rate_uses_up_to_three_decimal_places() -> None:
    assert _format_call_rate(1, 3) == "0.333"
    assert _format_call_rate(2, 3) == "0.667"
    assert _format_call_rate(3, 4) == "0.75"
    assert _format_call_rate(0, 3) == "0"
    assert _format_call_rate(0, 0) == ""


def test_export_biallelic_with_interest_counts(tmp_path: Path) -> None:
    output = tmp_path / "variants.csv"
    config = load_config(ROOT / "config.example.toml")
    config = replace(
        config,
        columns=(
            *config.columns,
            ColumnConfig(
                name="test_af",
                source="info",
                field="TEST_AF",
                aggregate="first",
            ),
            ColumnConfig(
                name="allele_label",
                source="info",
                field="ALLELE_LABEL",
            ),
            ColumnConfig(
                name="tags",
                source="info",
                field="TAGS",
                aggregate="join",
                separator=";",
            ),
        ),
    )
    stats = export_vcf(
        ROOT / "tests" / "data" / "biallelic.vcf",
        output,
        config,
        samples_of_interest=("S1", "S3"),
    )

    assert stats.variants_read == 2
    assert stats.rows_written == 2
    assert stats.samples == 4
    assert stats.rows_without_vep == 0

    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    first = rows[0]
    assert first["chrom"] == "1"
    assert first["rsid"] == "rsOne"
    assert first["call_rate"] == "0.75"
    assert first["filter"] == "PASS"
    assert first["het_samples"] == "S1"
    assert first["homalt_samples"] == "S2"
    assert first["interest_het_count"] == "1"
    assert first["interest_homalt_count"] == "0"
    assert first["other_het_count"] == "0"
    assert first["other_homalt_count"] == "1"
    assert first["gnomad_exomes_af"] == "0.0001"
    assert first["clinvar_significance"] == "Pathogenic"
    assert first["allele_label"] == "alternate"
    assert first["tags"] == "a;b"
    assert first["genes"] == "ATM"
    assert first["consequences"] == "missense_variant"
    assert first["transcripts"] == "ENST1"
    assert first["hgvsc"] == "ENST1:c.1A>G"
    assert first["hgvsp"] == "ENSP1:p.Lys1Arg"
    assert "vep_annotations" not in first

    alt_t = rows[1]
    assert alt_t["alt"] == "T"
    assert alt_t["filter"] == ""
    assert alt_t["het_samples"] == "S1,S3"
    assert alt_t["homalt_samples"] == "S2"
    assert alt_t["test_af"] == "0.01"
    assert alt_t["allele_label"] == "alternate-2"
    assert alt_t["interest_het_count"] == "2"
    assert alt_t["interest_homalt_count"] == "0"
    assert alt_t["other_het_count"] == "0"
    assert alt_t["other_homalt_count"] == "1"


def test_export_public_1000_genomes_multisample_subset(tmp_path: Path) -> None:
    output = tmp_path / "1000genomes.csv"
    config = load_config(ROOT / "tests" / "data" / "1000genomes.toml")

    stats = export_vcf(
        ROOT / "tests" / "data" / "1000genomes-phase3-subset.vcf",
        output,
        config,
        samples_of_interest=("HG00096", "NA12878"),
    )

    assert stats.variants_read == 2
    assert stats.rows_written == 2
    assert stats.samples == 4

    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert rows[0]["rsid"] == "rs4816203"
    assert rows[0]["global_af"] == "0.605431"
    assert rows[0]["call_rate"] == "1"
    assert rows[0]["het_samples"] == "HG00096,HG00097,NA12878"
    assert rows[0]["homalt_samples"] == ""
    assert rows[0]["interest_het_count"] == "2"
    assert rows[0]["other_het_count"] == "1"
    assert rows[0]["genes"] == "ANKEF1;SNAP25-AS1"
    assert rows[0]["consequence"] == (
        "intron_variant&non_coding_transcript_variant;downstream_gene_variant"
    )
    assert "vep_annotations" not in rows[0]

    assert rows[1]["rsid"] == "rs6057087"
    assert rows[1]["global_af"] == "0.810503"
    assert rows[1]["het_samples"] == "HG00096,HG00097,HG00099"
    assert rows[1]["homalt_samples"] == "NA12878"
    assert rows[1]["interest_het_count"] == "1"
    assert rows[1]["interest_homalt_count"] == "1"
    assert rows[1]["other_het_count"] == "2"
    assert rows[1]["other_homalt_count"] == "0"


def test_rejects_multiallelic_variant(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Multiallelic variant at 1:200"):
        export_vcf(
            ROOT / "tests" / "data" / "example.vcf",
            tmp_path / "variants.csv",
            load_config(ROOT / "config.example.toml"),
        )


def test_rejects_unknown_sample_of_interest(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="UNKNOWN"):
        export_vcf(
            ROOT / "tests" / "data" / "biallelic.vcf",
            tmp_path / "variants.csv",
            load_config(ROOT / "config.example.toml"),
            samples_of_interest=("UNKNOWN",),
        )


def test_gene_filter_matches_vep_gene_id(tmp_path: Path) -> None:
    output = tmp_path / "variants.csv"

    stats = export_vcf(
        ROOT / "tests" / "data" / "biallelic.vcf",
        output,
        load_config(ROOT / "config.example.toml"),
        genes=("ENSG2",),
    )

    assert stats.variants_read == 2
    assert stats.rows_written == 1
    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["rsid"] for row in rows] == ["rsTwo"]


def test_exports_site_only_vcf(tmp_path: Path) -> None:
    input_vcf = tmp_path / "site-only.vcf"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t1\t.\tA\tG\t.\tPASS\t.\n",
        encoding="utf-8",
    )
    output = tmp_path / "variants.csv"

    stats = export_vcf(
        input_vcf,
        output,
        config_from_mapping({"vep": {"required": False}}),
        samples_of_interest=(),
    )

    assert stats.samples == 0
    with output.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    assert row["n_samples"] == "0"
    assert row["interest_het_count"] == "0"
    assert row["other_het_count"] == "0"


def test_extracts_formatted_info_subfields(tmp_path: Path) -> None:
    input_vcf = tmp_path / "formatted-info.vcf"
    input_vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=1>\n"
        '##INFO=<ID=FV_GNOMAD,Number=.,Type=String,Description="Format: '
        'ALLELE|ALL_AF|ALL_AC">\n'
        '##INFO=<ID=CUSTOM,Number=1,Type=String,Description="Format: '
        'ALLELE~LABEL">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t1\t.\tA\tG\t.\tPASS\t"
        "FV_GNOMAD=G|0.012|7;CUSTOM=G~selected\n",
        encoding="utf-8",
    )
    output = tmp_path / "variants.csv"
    config = config_from_mapping(
        {
            "vep": {"required": False},
            "columns": [
                {
                    "name": "gnomad_all_af",
                    "source": "info",
                    "field": "FV_GNOMAD",
                    "subfield": "ALL_AF",
                },
                {
                    "name": "custom_label",
                    "source": "info",
                    "field": "CUSTOM",
                    "subfield": "LABEL",
                    "subfield_sep": "~",
                },
            ],
        }
    )

    export_vcf(input_vcf, output, config)

    with output.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    assert row["gnomad_all_af"] == "0.012"
    assert row["custom_label"] == "selected"


def test_rejects_unknown_required_info_subfield(tmp_path: Path) -> None:
    config = config_from_mapping(
        {
            "columns": [
                {
                    "name": "missing",
                    "source": "info",
                    "field": "TAGS",
                    "subfield": "MISSING",
                }
            ]
        }
    )

    with pytest.raises(ValueError, match="does not contain a 'Format:'"):
        export_vcf(
            ROOT / "tests" / "data" / "biallelic.vcf",
            tmp_path / "variants.csv",
            config,
        )


def test_rejects_duplicate_columns_from_direct_config(tmp_path: Path) -> None:
    config = load_config(ROOT / "config.example.toml")
    config = replace(config, vep=replace(config.vep, output_column="chrom"))

    with pytest.raises(ValueError, match="Duplicate output column"):
        export_vcf(
            ROOT / "tests" / "data" / "biallelic.vcf",
            tmp_path / "variants.csv",
            config,
        )


def test_json_vep_annotations_are_opt_in(tmp_path: Path) -> None:
    output = tmp_path / "variants.csv"
    config = config_from_mapping(
        {
            "vep": {
                "output_column": "vep_annotations",
                "include_fields": ["SYMBOL", "Consequence"],
            }
        }
    )

    export_vcf(ROOT / "tests" / "data" / "biallelic.vcf", output, config)

    with output.open(newline="", encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    assert row["vep_annotations"] == (
        '[{"SYMBOL":"ATM","Consequence":"missense_variant"}]'
    )


def test_rejects_same_input_and_output_path(tmp_path: Path) -> None:
    path = tmp_path / "variants.vcf"
    path.write_text(
        (ROOT / "tests" / "data" / "biallelic.vcf").read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="paths must differ"):
        export_vcf(path, path, load_config(ROOT / "config.example.toml"))
