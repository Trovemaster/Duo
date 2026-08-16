#!/usr/bin/env python3
"""Compare numerical results in two Duo output files.

Only the Python standard library is used.  The parser currently understands:

* rovibronic eigenvalues on lines containing ``||``;
* contracted vibrational energies on lines containing ``[ state v ]``.

Exit status 0 means that the benchmark passed, 1 means that differences were
found, and 2 means that an input file could not be parsed.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"

FINAL_RE = re.compile(
    rf"^\s*(?P<j>{NUMBER})\s+(?P<number>\d+)\s+(?P<energy>{NUMBER})\s+"
    rf"(?P<state>\d+)\s+(?P<v>\d+)\s+(?P<lambda>[+-]?\d+)\s+"
    rf"(?P<spin>{NUMBER})\s+(?P<sigma>{NUMBER})\s+"
    rf"(?P<omega>{NUMBER})\s+(?P<parity>\S+)\s*\|\|\s*(?P<name>.+?)\s*$"
)

VIB_RE = re.compile(
    rf"^\s*\d+\s+(?P<energy>{NUMBER})\s+"
    rf"\[\s*(?P<state>\d+)\s+(?P<v>\d+)\s*\]\s*(?P<name>.+?)\s*$"
)

CSV_COLUMNS = (
    "kind",
    "energy",
    "j",
    "number",
    "state",
    "v",
    "lambda",
    "spin",
    "sigma",
    "omega",
    "parity",
    "name",
)


def as_float(text: str) -> float:
    """Accept both the E and D exponent notation used by Fortran."""
    return float(text.replace("D", "E").replace("d", "e"))


def normal_number(text: str) -> str:
    """Make numerically equivalent quantum-number spellings compare equally."""
    value = as_float(text)
    if value == 0.0:
        value = 0.0  # Avoid different keys for -0.0 and 0.0.
    return f"{value:g}"


@dataclass(frozen=True)
class Record:
    kind: str
    energy: float
    fields: tuple[str, ...]

    @property
    def key(self) -> tuple[str, ...]:
        if self.kind == "final":
            # A final energy is rigorously identified by its position within
            # a J/parity block.  The assignment quantum numbers are diagnostic.
            j, number, _state, _v, _lam, _spin, _sigma, _omega, parity, _name = (
                self.fields
            )
            return (self.kind, j, parity, number)
        return (self.kind,) + self.fields

    @property
    def assignment(self) -> tuple[str, ...]:
        if self.kind != "final":
            return ()
        (_j, _number, state, v, lam, spin, sigma, omega, _parity, name) = self.fields
        return (state, v, lam, spin, sigma, omega, name)

    def label(self) -> str:
        if self.kind == "final":
            j, number, state, v, lam, spin, sigma, omega, parity, name = self.fields
            return (
                f"final J={j}, parity={parity}, i={number}, "
                f"state={state}, v={v}, lambda={lam}, "
                f"spin={spin}, sigma={sigma}, omega={omega}, "
                f"name={name}"
            )
        state, v, name = self.fields
        return f"vibrational state={state}, v={v}, name={name}"


@dataclass
class ParsedOutput:
    records: list[Record]
    final_count: int
    vibrational_count: int


@dataclass
class Difference:
    label: str
    reference: float
    actual: float
    delta: float
    limit: float


@dataclass
class AssignmentChange:
    label: str
    reference: tuple[str, ...]
    actual: tuple[str, ...]


@dataclass
class Comparison:
    reference: ParsedOutput
    actual: ParsedOutput
    differences: list[Difference]
    assignment_changes: list[AssignmentChange]
    missing: list[str]
    extra: list[str]
    compared: int
    maximum_error: float
    squared_error_sum: float

    @property
    def passed(self) -> bool:
        return not (self.differences or self.missing or self.extra)

    @property
    def rms_error(self) -> float:
        if not self.compared:
            return 0.0
        return math.sqrt(self.squared_error_sum / self.compared)


def final_record(match: re.Match[str]) -> Record:
    fields = (
        normal_number(match["j"]),
        match["number"],
        match["state"],
        match["v"],
        match["lambda"],
        normal_number(match["spin"]),
        normal_number(match["sigma"]),
        normal_number(match["omega"]),
        match["parity"],
        match["name"].strip(),
    )
    return Record("final", as_float(match["energy"]), fields)


def vibrational_record(match: re.Match[str]) -> Record:
    fields = (match["state"], match["v"], match["name"].strip())
    return Record("vibrational", as_float(match["energy"]), fields)


def parse_out(path: Path) -> ParsedOutput:
    records: list[Record] = []
    final_count = 0
    vibrational_count = 0

    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line in stream:
            if "||" in line:
                match = FINAL_RE.match(line)
                if match:
                    records.append(final_record(match))
                    final_count += 1
                continue

            if "[" in line and "]" in line:
                match = VIB_RE.match(line)
                if match:
                    records.append(vibrational_record(match))
                    vibrational_count += 1

    if final_count == 0 and vibrational_count == 0:
        raise ValueError(f"no supported Duo energy records found in {path}")
    return ParsedOutput(records, final_count, vibrational_count)


def record_to_row(record: Record) -> dict[str, str]:
    row = dict.fromkeys(CSV_COLUMNS, "")
    row["kind"] = record.kind
    row["energy"] = f"{record.energy:.15g}"
    if record.kind == "final":
        (
            row["j"],
            row["number"],
            row["state"],
            row["v"],
            row["lambda"],
            row["spin"],
            row["sigma"],
            row["omega"],
            row["parity"],
            row["name"],
        ) = record.fields
    else:
        row["state"], row["v"], row["name"] = record.fields
    return row


def write_snapshot(output: ParsedOutput, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for record in output.records:
            writer.writerow(record_to_row(record))


def parse_snapshot(path: Path) -> ParsedOutput:
    records: list[Record] = []
    counts = {"final": 0, "vibrational": 0}
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != list(CSV_COLUMNS):
            raise ValueError(f"{path} is not a Duo benchmark snapshot")
        for line_number, row in enumerate(reader, start=2):
            kind = row["kind"]
            if kind == "final":
                fields = tuple(row[name] for name in CSV_COLUMNS[2:])
            elif kind == "vibrational":
                fields = (row["state"], row["v"], row["name"])
            else:
                raise ValueError(f"unknown record kind on {path}:{line_number}: {kind!r}")
            records.append(Record(kind, as_float(row["energy"]), fields))
            counts[kind] += 1

    if not records:
        raise ValueError(f"empty Duo benchmark snapshot: {path}")
    return ParsedOutput(records, counts["final"], counts["vibrational"])


def parse_reference(path: Path) -> ParsedOutput:
    return parse_snapshot(path) if path.suffix.lower() == ".csv" else parse_out(path)


def group_records(records: Iterable[Record]) -> dict[tuple[str, ...], list[Record]]:
    grouped: dict[tuple[str, ...], list[Record]] = defaultdict(list)
    for record in records:
        grouped[record.key].append(record)
    for values in grouped.values():
        values.sort(key=lambda record: record.energy)
    return grouped


def compare(
    reference: ParsedOutput,
    actual: ParsedOutput,
    final_atol: float,
    vibrational_atol: float,
    rtol: float,
) -> Comparison:
    expected = group_records(reference.records)
    observed = group_records(actual.records)
    differences: list[Difference] = []
    assignment_changes: list[AssignmentChange] = []
    missing: list[str] = []
    extra: list[str] = []
    compared = 0
    maximum_error = 0.0
    squared_error_sum = 0.0

    for key in sorted(expected.keys() | observed.keys()):
        expected_values = expected.get(key, [])
        observed_values = observed.get(key, [])
        common = min(len(expected_values), len(observed_values))

        for index in range(common):
            expected_record = expected_values[index]
            observed_record = observed_values[index]
            delta = observed_record.energy - expected_record.energy
            error = abs(delta)
            atol = final_atol if expected_record.kind == "final" else vibrational_atol
            limit = atol + rtol * abs(expected_record.energy)
            compared += 1
            maximum_error = max(maximum_error, error)
            squared_error_sum += error * error
            if expected_record.assignment != observed_record.assignment:
                assignment_changes.append(
                    AssignmentChange(
                        # Show only the rigorous identity in the heading.
                        label=(
                            f"J={expected_record.fields[0]}, "
                            f"parity={expected_record.fields[8]}, "
                            f"i={expected_record.fields[1]}"
                        ),
                        reference=expected_record.assignment,
                        actual=observed_record.assignment,
                    )
                )
            if error > limit:
                differences.append(
                    Difference(
                        expected_record.label(),
                        expected_record.energy,
                        observed_record.energy,
                        delta,
                        limit,
                    )
                )

        missing.extend(record.label() for record in expected_values[common:])
        extra.extend(record.label() for record in observed_values[common:])

    differences.sort(key=lambda item: abs(item.delta), reverse=True)
    return Comparison(
        reference,
        actual,
        differences,
        assignment_changes,
        missing,
        extra,
        compared,
        maximum_error,
        squared_error_sum,
    )


def print_report(result: Comparison, max_report: int) -> None:
    status = "PASS" if result.passed else "FAIL"
    print(f"Duo benchmark: {status}")
    print(
        "Records: "
        f"reference={len(result.reference.records)} "
        f"(final={result.reference.final_count}, "
        f"vibrational={result.reference.vibrational_count}), "
        f"actual={len(result.actual.records)} "
        f"(final={result.actual.final_count}, "
        f"vibrational={result.actual.vibrational_count})"
    )
    print(
        f"Matched={result.compared}, changed={len(result.differences)}, "
        f"missing={len(result.missing)}, extra={len(result.extra)}, "
        f"assignment comments={len(result.assignment_changes)}"
    )
    print(
        f"Maximum absolute error={result.maximum_error:.6g} cm-1, "
        f"RMS error={result.rms_error:.6g} cm-1"
    )

    shown = 0
    for difference in result.differences:
        if shown >= max_report:
            break
        print(
            f"CHANGED: {difference.label}\n"
            f"  reference={difference.reference:.9f}, "
            f"actual={difference.actual:.9f}, delta={difference.delta:+.9f}, "
            f"limit={difference.limit:.6g} cm-1"
        )
        shown += 1
    for prefix, items in (("MISSING", result.missing), ("EXTRA", result.extra)):
        for label in items:
            if shown >= max_report:
                break
            print(f"{prefix}: {label}")
            shown += 1

    total_problems = len(result.differences) + len(result.missing) + len(result.extra)
    if total_problems > shown:
        print(f"... {total_problems - shown} more problem(s) not shown")

    if result.assignment_changes:
        field_names = ("state", "v", "lambda", "spin", "sigma", "omega", "name")
        print("Quantum-number assignment comments (do not affect PASS/FAIL):")
        for change in result.assignment_changes[:max_report]:
            changed_fields = []
            for name, reference, actual in zip(
                field_names, change.reference, change.actual
            ):
                if reference != actual:
                    changed_fields.append(f"{name}: {reference} -> {actual}")
            print(f"COMMENT: {change.label}: " + ", ".join(changed_fields))
        if len(result.assignment_changes) > max_report:
            print(
                f"... {len(result.assignment_changes) - max_report} "
                "more assignment comment(s) not shown"
            )


def non_negative(text: str) -> float:
    value = float(text)
    if value < 0.0:
        raise argparse.ArgumentTypeError("must not be negative")
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Benchmark Duo energies using a reference .out file or CSV snapshot."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    compare_parser = subparsers.add_parser("compare", help="compare a new Duo run")
    compare_parser.add_argument("reference", type=Path, help="reference .out or snapshot .csv")
    compare_parser.add_argument("actual", type=Path, help="new Duo .out file")
    compare_parser.add_argument(
        "--final-atol",
        type=non_negative,
        default=1.0e-3,
        help="absolute tolerance for || energies in cm-1 (default: 0.001)",
    )
    compare_parser.add_argument(
        "--vibrational-atol",
        type=non_negative,
        default=1.0e-3,
        help="absolute tolerance for [state v] energies in cm-1 (default: 0.001)",
    )
    compare_parser.add_argument(
        "--rtol",
        type=non_negative,
        default=0.0,
        help="relative tolerance added to the absolute tolerance (default: 0)",
    )
    compare_parser.add_argument(
        "--max-report",
        type=int,
        default=20,
        help="maximum number of individual problems to print (default: 20)",
    )

    snapshot_parser = subparsers.add_parser(
        "snapshot", help="extract a portable CSV reference from a Duo .out file"
    )
    snapshot_parser.add_argument("source", type=Path, help="trusted Duo .out file")
    snapshot_parser.add_argument("destination", type=Path, help="new snapshot .csv")
    return parser


def run(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "snapshot":
            parsed = parse_out(args.source)
            write_snapshot(parsed, args.destination)
            print(
                f"Wrote {len(parsed.records)} records "
                f"(final={parsed.final_count}, vibrational={parsed.vibrational_count}) "
                f"to {args.destination}"
            )
            return 0

        reference = parse_reference(args.reference)
        actual = parse_out(args.actual)
        result = compare(
            reference,
            actual,
            final_atol=args.final_atol,
            vibrational_atol=args.vibrational_atol,
            rtol=args.rtol,
        )
        print_report(result, max(0, args.max_report))
        return 0 if result.passed else 1
    except (OSError, ValueError, csv.Error) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(run())
