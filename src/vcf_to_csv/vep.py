"""Parsing helpers for VEP-compatible CSQ annotations."""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

DEFAULT_VEP_FIELDS = (
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "Existing_variation",
    "CANONICAL",
    "MANE_SELECT",
    "ALLELE_NUM",
)


def parse_vep_format(description: str) -> tuple[str, ...]:
    """Extract ordered subfield names from a VEP INFO header description."""

    match = re.search(r"\bFormat:\s*([^\"\r\n]+)", description, flags=re.IGNORECASE)
    if not match:
        raise ValueError("VEP INFO header does not contain a 'Format:' declaration")
    fields = tuple(part.strip() for part in match.group(1).strip().split("|"))
    if not fields or any(not item for item in fields):
        raise ValueError("VEP INFO header contains an invalid field list")
    return fields


@dataclass(frozen=True)
class VepParser:
    """Parse and project transcript annotations from a VEP-compatible INFO field."""

    fields: tuple[str, ...]
    output_fields: tuple[str, ...]

    @classmethod
    def from_header(
        cls,
        description: str,
        requested_output_fields: Sequence[str] | None,
    ) -> VepParser:
        fields = parse_vep_format(description)
        if requested_output_fields is None:
            output_fields = tuple(
                field for field in DEFAULT_VEP_FIELDS if field in fields
            )
            if not output_fields:
                output_fields = fields
        else:
            missing = sorted(set(requested_output_fields) - set(fields))
            if missing:
                raise ValueError(
                    "Configured VEP output field(s) are absent from the header: "
                    + ", ".join(missing)
                )
            output_fields = tuple(requested_output_fields)
        return cls(fields=fields, output_fields=output_fields)

    def parse(self, raw_value: object) -> list[dict[str, str]]:
        """Decode a CSQ value into one mapping per transcript annotation."""

        if raw_value is None:
            return []
        if isinstance(raw_value, (list, tuple)):
            text = ",".join(str(item) for item in raw_value)
        else:
            text = str(raw_value)
        if not text or text == ".":
            return []

        records = []
        for raw_record in text.split(","):
            values = raw_record.split("|")
            if len(values) < len(self.fields):
                values.extend([""] * (len(self.fields) - len(values)))
            records.append(dict(zip(self.fields, values)))
        return records

    def project(self, records: Iterable[Mapping[str, str]]) -> list[dict[str, str]]:
        """Keep only configured fields while preserving transcript order."""

        return [
            {field: record.get(field, "") for field in self.output_fields}
            for record in records
        ]
