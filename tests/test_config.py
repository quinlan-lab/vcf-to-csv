from pathlib import Path

import pytest

from vcf_to_csv.config import config_from_mapping, load_config


def test_load_example_config() -> None:
    path = Path(__file__).parents[1] / "config.example.toml"
    config = load_config(path)

    assert config.vep.info_field == "CSQ"
    assert config.vep.output_column is None
    assert [column.name for column in config.columns] == [
        "genes",
        "consequences",
        "transcripts",
        "hgvsc",
        "hgvsp",
        "gnomad_exomes_af",
        "gnomad_genomes_af",
        "clinvar_significance",
    ]


def test_rejects_reserved_column_name() -> None:
    with pytest.raises(ValueError, match="reserved"):
        config_from_mapping(
            {"columns": [{"name": "chrom", "source": "info", "field": "CHROM"}]}
        )


def test_rejects_reserved_impact_column_name() -> None:
    with pytest.raises(ValueError, match="reserved"):
        config_from_mapping(
            {"columns": [{"name": "impact", "source": "vep", "field": "IMPACT"}]}
        )


def test_rejects_reserved_vep_output_column_name() -> None:
    with pytest.raises(ValueError, match="reserved"):
        config_from_mapping({"vep": {"output_column": "chrom"}})


def test_include_fields_requires_json_output_column() -> None:
    with pytest.raises(ValueError, match="requires vep.output_column"):
        config_from_mapping({"vep": {"include_fields": ["SYMBOL"]}})


def test_info_subfield_separator_defaults_to_pipe() -> None:
    config = config_from_mapping(
        {
            "columns": [
                {
                    "name": "all_af",
                    "source": "info",
                    "field": "FV_GNOMAD",
                    "subfield": "ALL_AF",
                }
            ]
        }
    )

    assert config.columns[0].subfield_sep == "|"


def test_rejects_subfield_on_vep_column() -> None:
    with pytest.raises(ValueError, match="only valid for source='info'"):
        config_from_mapping(
            {
                "columns": [
                    {
                        "name": "symbol",
                        "source": "vep",
                        "field": "SYMBOL",
                        "subfield": "value",
                    }
                ]
            }
        )


def test_rejects_removed_output_settings() -> None:
    with pytest.raises(ValueError, match="Unknown top-level"):
        config_from_mapping({"output": {"sample_separator": ";"}})
