"""Configuration models and TOML loading for the VCF exporter."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.9/3.10 only
    import tomli as tomllib  # type: ignore[no-redef]


CORE_COLUMNS = (
    "chrom",
    "pos",
    "ref",
    "alt",
    "alt_index",
    "rsid",
    "qual",
    "filter",
    "call_rate",
    "n_samples",
    "n_called",
    "n_missing",
    "het_count",
    "homalt_count",
    "interest_het_count",
    "interest_homalt_count",
    "other_het_count",
    "other_homalt_count",
    "het_samples",
    "homalt_samples",
)
VALID_SOURCES = frozenset({"info", "vep"})
VALID_AGGREGATES = frozenset({"first", "join", "max", "min"})


@dataclass(frozen=True)
class VepConfig:
    """Settings for the VEP-compatible INFO field."""

    info_field: str = "CSQ"
    output_column: str | None = None
    include_fields: tuple[str, ...] | None = None
    required: bool = True


@dataclass(frozen=True)
class ColumnConfig:
    """One additional CSV column sourced from INFO or a VEP subfield."""

    name: str
    source: str
    field: str
    aggregate: str = "first"
    separator: str = ";"
    required: bool = True
    subfield: str | None = None
    subfield_sep: str = "|"


@dataclass(frozen=True)
class ExportConfig:
    """Complete exporter configuration."""

    vep: VepConfig = field(default_factory=VepConfig)
    columns: tuple[ColumnConfig, ...] = ()


def _check_keys(mapping: Mapping[str, Any], allowed: set[str], context: str) -> None:
    unknown = sorted(set(mapping) - allowed)
    if unknown:
        raise ValueError(f"Unknown {context} setting(s): {', '.join(unknown)}")


def _as_nonempty_string(value: Any, context: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{context} must be a non-empty string")
    return value


def _as_bool(value: Any, context: str) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{context} must be true or false")
    return value


def config_from_mapping(data: Mapping[str, Any]) -> ExportConfig:
    """Build and validate :class:`ExportConfig` from decoded TOML data."""

    _check_keys(data, {"vep", "columns"}, "top-level")

    vep_data = data.get("vep", {})
    if not isinstance(vep_data, Mapping):
        raise ValueError("[vep] must be a TOML table")
    _check_keys(
        vep_data,
        {"info_field", "output_column", "include_fields", "required"},
        "[vep]",
    )
    include_fields_value = vep_data.get("include_fields")
    include_fields: tuple[str, ...] | None
    if include_fields_value is None:
        include_fields = None
    elif isinstance(include_fields_value, list) and all(
        isinstance(item, str) and item for item in include_fields_value
    ):
        include_fields = tuple(include_fields_value)
    else:
        raise ValueError("vep.include_fields must be an array of non-empty strings")

    output_column_value = vep_data.get("output_column")
    output_column = (
        None
        if output_column_value is None
        else _as_nonempty_string(output_column_value, "vep.output_column")
    )
    if output_column is None and include_fields is not None:
        raise ValueError("vep.include_fields requires vep.output_column")

    vep = VepConfig(
        info_field=_as_nonempty_string(
            vep_data.get("info_field", "CSQ"), "vep.info_field"
        ),
        output_column=output_column,
        include_fields=include_fields,
        required=_as_bool(vep_data.get("required", True), "vep.required"),
    )
    if vep.output_column is not None and vep.output_column in CORE_COLUMNS:
        raise ValueError(
            f"Duplicate or reserved output column name: {vep.output_column}"
        )

    raw_columns = data.get("columns", [])
    if not isinstance(raw_columns, list):
        raise ValueError("[[columns]] entries must be TOML tables")

    columns = []
    seen_names = set(CORE_COLUMNS)
    if vep.output_column is not None:
        seen_names.add(vep.output_column)
    for index, raw_column in enumerate(raw_columns, start=1):
        context = f"[[columns]] entry {index}"
        if not isinstance(raw_column, Mapping):
            raise ValueError(f"{context} must be a TOML table")
        _check_keys(
            raw_column,
            {
                "name",
                "source",
                "field",
                "aggregate",
                "separator",
                "required",
                "subfield",
                "subfield_sep",
            },
            context,
        )
        name = _as_nonempty_string(raw_column.get("name"), f"{context}.name")
        source = _as_nonempty_string(
            raw_column.get("source"), f"{context}.source"
        ).lower()
        field_name = _as_nonempty_string(raw_column.get("field"), f"{context}.field")
        aggregate = str(raw_column.get("aggregate", "first")).lower()
        separator = _as_nonempty_string(
            raw_column.get("separator", ";"), f"{context}.separator"
        )
        subfield_value = raw_column.get("subfield")
        subfield = (
            None
            if subfield_value is None
            else _as_nonempty_string(subfield_value, f"{context}.subfield")
        )
        subfield_sep = _as_nonempty_string(
            raw_column.get("subfield_sep", "|"), f"{context}.subfield_sep"
        )

        if name in seen_names:
            raise ValueError(f"Duplicate or reserved output column name: {name}")
        if source not in VALID_SOURCES:
            raise ValueError(
                f"{context}.source must be one of: {', '.join(sorted(VALID_SOURCES))}"
            )
        if aggregate not in VALID_AGGREGATES:
            raise ValueError(
                f"{context}.aggregate must be one of: "
                f"{', '.join(sorted(VALID_AGGREGATES))}"
            )
        if source != "info" and subfield is not None:
            raise ValueError(f"{context}.subfield is only valid for source='info'")
        if subfield is None and "subfield_sep" in raw_column:
            raise ValueError(f"{context}.subfield_sep requires subfield")

        seen_names.add(name)
        columns.append(
            ColumnConfig(
                name=name,
                source=source,
                field=field_name,
                aggregate=aggregate,
                separator=separator,
                required=_as_bool(
                    raw_column.get("required", True), f"{context}.required"
                ),
                subfield=subfield,
                subfield_sep=subfield_sep,
            )
        )

    return ExportConfig(vep=vep, columns=tuple(columns))


def load_config(path: Path | None = None) -> ExportConfig:
    """Load exporter settings from TOML, or return defaults when no path is given."""

    if path is None:
        return ExportConfig()
    with path.open("rb") as handle:
        data = tomllib.load(handle)
    return config_from_mapping(data)
