"""VCF-to-CSV export implementation."""

from __future__ import annotations

import csv
import json
import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO

from cyvcf2 import VCF

from .config import CORE_COLUMNS, ColumnConfig, ExportConfig
from .vep import VepParser

SAMPLE_SEPARATOR = ","
MISSING_VALUE = ""
FLOAT_SIGNIFICANT_DIGITS = 7


@dataclass(frozen=True)
class ExportStats:
    """Counts reported after a completed export."""

    variants_read: int
    rows_written: int
    samples: int
    rows_without_vep: int


@dataclass(frozen=True)
class HeaderField:
    """VCF header metadata needed for INFO extraction."""

    number: str
    value_type: str


def _header_record(vcf: VCF, field_name: str) -> Mapping[str, str] | None:
    try:
        return vcf.get_header_type(field_name)
    except KeyError:
        return None


def _build_vep_parser(vcf: VCF, config: ExportConfig) -> VepParser | None:
    header = _header_record(vcf, config.vep.info_field)
    if header is None:
        if config.vep.required:
            raise ValueError(
                f"Required VEP INFO field {config.vep.info_field!r} is absent "
                "from the VCF header"
            )
        return None
    description = header.get("Description", "")
    return VepParser.from_header(description, config.vep.include_fields)


def _validate_custom_columns(
    vcf: VCF,
    config: ExportConfig,
    vep_parser: VepParser | None,
) -> dict[str, HeaderField]:
    info_headers: dict[str, HeaderField] = {}
    for column in config.columns:
        if column.source == "vep":
            field_exists = vep_parser is not None and column.field in vep_parser.fields
            if column.required and not field_exists:
                raise ValueError(
                    f"Required VEP field {column.field!r} for column {column.name!r} "
                    "is absent from the VEP header"
                )
            continue

        header = _header_record(vcf, column.field)
        if header is None:
            if column.required:
                raise ValueError(
                    f"Required INFO field {column.field!r} for column {column.name!r} "
                    "is absent from the VCF header"
                )
            continue
        info_headers[column.field] = HeaderField(
            number=str(header.get("Number", ".")),
            value_type=str(header.get("Type", "String")),
        )
    return info_headers


def _called_and_carrier_samples(
    sample_names: Sequence[str],
    genotypes: Sequence[Sequence[object]],
) -> tuple[int, list[str], list[str]]:
    """Classify carriers for a biallelic variant using strict genotypes."""

    called = 0
    het_samples: list[str] = []
    homalt_samples: list[str] = []

    for sample_name, genotype in zip(sample_names, genotypes):
        # cyvcf2 appends a phasing boolean to the allele indexes.
        alleles = genotype[:-1]
        if not alleles or any(allele is None or int(allele) < 0 for allele in alleles):
            continue

        called += 1
        allele_indexes = [int(allele) for allele in alleles]
        if 1 not in allele_indexes:
            continue
        if all(allele == 1 for allele in allele_indexes):
            homalt_samples.append(sample_name)
        else:
            het_samples.append(sample_name)

    return called, het_samples, homalt_samples


def _info_value(
    value: object,
    header: HeaderField,
) -> object:
    """Normalize an INFO value and select its biallelic ALT component."""

    if (
        isinstance(value, str)
        and header.value_type in {"String", "Character"}
        and header.number != "1"
    ):
        value = value.split(",")

    if header.number not in {"A", "R"}:
        return value

    values = list(value) if isinstance(value, (list, tuple)) else [value]
    index = 0 if header.number == "A" else 1
    return values[index] if index < len(values) else None


def _iter_atomic_values(values: Iterable[object]) -> Iterable[object]:
    for value in values:
        if value is None:
            continue
        if isinstance(value, (list, tuple)):
            yield from _iter_atomic_values(value)
        else:
            yield value


def _clean_scalar(value: object) -> object | None:
    if value is None or value == "." or value == "":
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


def _aggregate(values: Iterable[object], column: ColumnConfig) -> object | None:
    cleaned = [
        item
        for item in (_clean_scalar(value) for value in _iter_atomic_values(values))
        if item is not None
    ]
    if not cleaned:
        return None
    if column.aggregate == "first":
        return cleaned[0]
    if column.aggregate == "join":
        unique = []
        seen = set()
        for value in cleaned:
            text = str(value)
            if text not in seen:
                seen.add(text)
                unique.append(text)
        return column.separator.join(unique)

    numeric_values = []
    for value in cleaned:
        parts = str(value).replace("&", ",").split(",")
        try:
            numeric_values.extend(
                float(part) for part in parts if part not in {"", "."}
            )
        except ValueError as error:
            raise ValueError(
                f"Column {column.name!r} uses aggregate={column.aggregate!r}, "
                f"but value {value!r} is not numeric"
            ) from error
    if not numeric_values:
        return None
    if column.aggregate == "max":
        return max(numeric_values)
    if column.aggregate == "min":
        return min(numeric_values)
    raise ValueError(f"Unsupported aggregation mode: {column.aggregate!r}")


def _custom_column_value(
    variant: object,
    records: Sequence[Mapping[str, str]],
    column: ColumnConfig,
    info_headers: Mapping[str, HeaderField],
) -> object | None:
    if column.source == "vep":
        return _aggregate((record.get(column.field) for record in records), column)

    header = info_headers.get(column.field)
    if header is None:
        return None
    raw_value = variant.INFO.get(column.field)
    selected = _info_value(raw_value, header)
    return _aggregate([selected], column)


def _display_value(value: object | None) -> object:
    cleaned = _clean_scalar(value)
    if cleaned is None:
        return MISSING_VALUE
    if isinstance(cleaned, float):
        return format(cleaned, f".{FLOAT_SIGNIFICANT_DIGITS}g")
    return cleaned


def _format_call_rate(called: int, total: int) -> str:
    if not total:
        return MISSING_VALUE
    value = format(called / total, ".3f")
    return value.rstrip("0").rstrip(".") if "." in value else value


def _open_output(path: Path | str) -> tuple[TextIO, bool]:
    if str(path) == "-":
        import sys

        return sys.stdout, False
    # The caller owns this handle and closes it in a finally block. A context
    # manager cannot also wrap stdout, which must stay open after export.
    handle = Path(path).open("w", encoding="utf-8", newline="")  # noqa: SIM115
    return handle, True


def export_vcf(
    input_path: Path | str,
    output_path: Path | str,
    config: ExportConfig,
    samples_of_interest: Iterable[str] | None = None,
) -> ExportStats:
    """Export an annotated biallelic VCF to CSV."""

    if (
        str(output_path) != "-"
        and Path(input_path).resolve() == Path(output_path).resolve()
    ):
        raise ValueError("Input VCF and output CSV paths must differ")

    variants_read = 0
    rows_written = 0
    rows_without_vep = 0
    vcf = VCF(str(input_path), gts012=True, strict_gt=True)
    try:
        sample_names = tuple(vcf.samples)
        interest_samples = frozenset(samples_of_interest or ())
        unknown_samples = sorted(interest_samples - set(sample_names))
        if unknown_samples:
            raise ValueError(
                "Sample(s) of interest are absent from the VCF header: "
                + ", ".join(unknown_samples)
            )
        vep_parser = _build_vep_parser(vcf, config)
        info_headers = _validate_custom_columns(vcf, config, vep_parser)
        fieldnames = (
            *CORE_COLUMNS,
            config.vep.output_column,
            *(column.name for column in config.columns),
        )
        duplicate_fieldnames = sorted(
            name for name in set(fieldnames) if fieldnames.count(name) > 1
        )
        if duplicate_fieldnames:
            raise ValueError(
                "Duplicate output column name(s): " + ", ".join(duplicate_fieldnames)
            )

        output_handle, should_close = _open_output(output_path)
        try:
            writer = csv.DictWriter(output_handle, fieldnames=fieldnames)
            writer.writeheader()

            for variant in vcf:
                variants_read += 1
                if len(variant.ALT) > 1:
                    raise ValueError(
                        f"Multiallelic variant at {variant.CHROM}:{variant.POS} has "
                        f"{len(variant.ALT)} ALT alleles; split multiallelic records "
                        "before export"
                    )
                if not variant.ALT or variant.ALT[0] is None:
                    raise ValueError(
                        f"Variant at {variant.CHROM}:{variant.POS} has no ALT allele"
                    )

                raw_vep = (
                    variant.INFO.get(config.vep.info_field)
                    if vep_parser is not None
                    else None
                )
                vep_records = vep_parser.parse(raw_vep) if vep_parser else []
                try:
                    genotypes = variant.genotypes or ()
                except Exception as error:
                    raise ValueError(
                        f"Could not read GT at {variant.CHROM}:{variant.POS}"
                    ) from error
                n_samples = len(sample_names)
                called, het_samples, homalt_samples = _called_and_carrier_samples(
                    sample_names, genotypes
                )
                projected_records = (
                    vep_parser.project(vep_records) if vep_parser else []
                )
                interest_het_count = len(interest_samples.intersection(het_samples))
                interest_homalt_count = len(
                    interest_samples.intersection(homalt_samples)
                )
                if vep_parser is not None and not vep_records:
                    rows_without_vep += 1

                row = {
                    "chrom": variant.CHROM,
                    "pos": variant.POS,
                    "ref": variant.REF,
                    "alt": variant.ALT[0],
                    "alt_index": 1,
                    "rsid": _display_value(variant.ID),
                    "qual": _display_value(variant.QUAL),
                    "filter": ";".join(variant.FILTERS) or MISSING_VALUE,
                    "call_rate": _format_call_rate(called, n_samples),
                    "n_samples": n_samples,
                    "n_called": called,
                    "n_missing": n_samples - called,
                    "het_count": len(het_samples),
                    "homalt_count": len(homalt_samples),
                    "interest_het_count": interest_het_count,
                    "interest_homalt_count": interest_homalt_count,
                    "other_het_count": len(het_samples) - interest_het_count,
                    "other_homalt_count": len(homalt_samples) - interest_homalt_count,
                    "het_samples": SAMPLE_SEPARATOR.join(het_samples),
                    "homalt_samples": SAMPLE_SEPARATOR.join(homalt_samples),
                    config.vep.output_column: json.dumps(
                        projected_records,
                        ensure_ascii=False,
                        separators=(",", ":"),
                    ),
                }
                for column in config.columns:
                    value = _custom_column_value(
                        variant,
                        vep_records,
                        column,
                        info_headers,
                    )
                    row[column.name] = _display_value(value)

                writer.writerow(row)
                rows_written += 1
        finally:
            if should_close:
                output_handle.close()
    finally:
        vcf.close()

    return ExportStats(
        variants_read=variants_read,
        rows_written=rows_written,
        samples=len(sample_names),
        rows_without_vep=rows_without_vep,
    )
