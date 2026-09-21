"""Command-line entry point."""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .config import load_config
from .exporter import export_vcf


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="vcf-to-csv",
        description=(
            "Export annotated VCF records to CSV with allele-specific HET and "
            "HOMALT sample lists."
        ),
    )
    parser.add_argument("input_vcf", type=Path, help="Input .vcf, .vcf.gz, or .bcf")
    parser.add_argument(
        "output_csv",
        help="Output CSV path, or '-' to write CSV to standard output",
    )
    parser.add_argument(
        "--config",
        type=Path,
        help="Optional TOML file defining VEP and additional INFO columns",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace an existing output CSV",
    )
    parser.add_argument(
        "--samples-of-interest",
        type=Path,
        metavar="FILE",
        help="Possibly empty text file with one sample ID of interest per line",
    )
    parser.add_argument("--version", action="version", version=__version__)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.input_vcf.is_file():
        print(f"error: input VCF does not exist: {args.input_vcf}", file=sys.stderr)
        return 2

    if args.output_csv != "-":
        output_path = Path(args.output_csv)
        if output_path.exists() and not args.force:
            message = (
                f"error: output already exists: {output_path} "
                "(use --force to replace it)"
            )
            print(
                message,
                file=sys.stderr,
            )
            return 2
        output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        config = load_config(args.config)
        samples_of_interest = None
        if args.samples_of_interest is not None:
            samples_of_interest = _load_samples_of_interest(args.samples_of_interest)
        stats = export_vcf(
            args.input_vcf,
            args.output_csv,
            config,
            samples_of_interest=samples_of_interest,
        )
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1

    if args.output_csv != "-":
        print(
            f"Wrote {stats.rows_written} CSV row(s) from "
            f"{stats.variants_read} variant record(s) and {stats.samples} sample(s).",
            file=sys.stderr,
        )
        if stats.rows_without_vep:
            message = (
                f"Warning: {stats.rows_without_vep} row(s) had no matching "
                "VEP annotation."
            )
            print(
                message,
                file=sys.stderr,
            )
    return 0


def _load_samples_of_interest(path: Path) -> frozenset[str]:
    samples: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            sample = raw_line.strip()
            if not sample or sample.startswith("#"):
                continue
            if any(character.isspace() for character in sample):
                raise ValueError(
                    f"Invalid sample ID at {path}:{line_number}; "
                    "expected one ID per line"
                )
            if sample in samples:
                raise ValueError(
                    f"Duplicate sample ID {sample!r} at {path}:{line_number}"
                )
            samples.add(sample)
    return frozenset(samples)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
