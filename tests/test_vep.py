from vcf_to_csv.vep import VepParser, parse_info_format, parse_vep_format


def test_parse_vep_format() -> None:
    description = (
        "Consequence annotations from Ensembl VEP. "
        "Format: Allele|Consequence|ALLELE_NUM"
    )
    fields = parse_vep_format(description)
    assert fields == ("Allele", "Consequence", "ALLELE_NUM")


def test_parse_info_format_with_custom_separator() -> None:
    fields = parse_info_format("Annotations. Format: ALLELE~SCORE", "~")
    assert fields == ("ALLELE", "SCORE")


def test_parses_and_projects_records() -> None:
    parser = VepParser.from_header("Format: Allele|Consequence|ALLELE_NUM", None)
    records = parser.parse("T|synonymous_variant|1")

    assert parser.project(records) == [
        {
            "Allele": "T",
            "Consequence": "synonymous_variant",
            "ALLELE_NUM": "1",
        }
    ]
