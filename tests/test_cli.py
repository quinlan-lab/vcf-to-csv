import csv
from pathlib import Path

from vcf_to_csv.cli import main

ROOT = Path(__file__).parents[1]


def test_cli_loads_samples_of_interest_file(tmp_path: Path) -> None:
    sample_file = tmp_path / "samples.txt"
    sample_file.write_text("# selected samples\nS1\n\nS3\n", encoding="utf-8")
    output = tmp_path / "variants.csv"

    result = main(
        [
            str(ROOT / "tests" / "data" / "biallelic.vcf"),
            str(output),
            "--config",
            str(ROOT / "config.example.toml"),
            "--samples-of-interest",
            str(sample_file),
        ]
    )

    assert result == 0
    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[1]["interest_het_count"] == "2"
    assert rows[1]["other_homalt_count"] == "1"


def test_cli_accepts_empty_samples_of_interest_file(tmp_path: Path) -> None:
    sample_file = tmp_path / "samples.txt"
    sample_file.write_text("", encoding="utf-8")
    output = tmp_path / "variants.csv"

    result = main(
        [
            str(ROOT / "tests" / "data" / "biallelic.vcf"),
            str(output),
            "--config",
            str(ROOT / "config.example.toml"),
            "--samples-of-interest",
            str(sample_file),
        ]
    )

    assert result == 0
    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert rows[1]["interest_het_count"] == "0"
    assert rows[1]["other_het_count"] == "2"


def test_cli_filters_variants_by_gene_list(tmp_path: Path) -> None:
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("# genes of interest\nATM\n", encoding="utf-8")
    output = tmp_path / "variants.csv"

    result = main(
        [
            str(ROOT / "tests" / "data" / "biallelic.vcf"),
            str(output),
            "--config",
            str(ROOT / "config.example.toml"),
            "--genes",
            str(gene_file),
        ]
    )

    assert result == 0
    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [row["rsid"] for row in rows] == ["rsOne"]
