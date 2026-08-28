#!/usr/bin/env python3
"""Compare two Duo intensity-fitting ``.int`` files.

Only the Python standard library is used.  Records are grouped by their
rigorous transition labels (direction, upper/lower J and symmetry, and branch),
so output line order is irrelevant.  Assignment quantum numbers help pair
records inside a group, but assignment changes are comments and do not affect
PASS/FAIL.

Exit status 0 means PASS, 1 means benchmark differences were found, and 2
means that a file could not be read or parsed (or there was insufficient
memory).
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Sequence


@dataclass(frozen=True)
class IntRecord:
    line_number: int
    direction: str
    upper_two_j: int
    upper_symmetry: str
    lower_two_j: int
    lower_symmetry: str
    branch: str
    upper_energy: float
    lower_energy: float
    frequency: float
    line_strength: float
    einstein_a: float
    calculated_intensity: float
    reference_intensity: float
    assignment: tuple[int, ...]

    @property
    def rigorous_key(self) -> tuple[object, ...]:
        return (
            self.direction,
            self.upper_two_j,
            self.upper_symmetry,
            self.lower_two_j,
            self.lower_symmetry,
            self.branch,
        )

    @property
    def energy_sort_key(self) -> tuple[float, ...]:
        # Used only as a fallback if assignment quantum numbers have changed.
        return (
            self.upper_energy,
            self.lower_energy,
            self.frequency,
            self.reference_intensity,
        )


@dataclass
class QuantityComparison:
    changed: int = 0
    maximum_absolute_error: float = 0.0
    maximum_relative_error: float = 0.0

    def add(
        self,
        reference: float,
        actual: float,
        absolute_tolerance: float,
        relative_tolerance: float,
    ) -> tuple[bool, float, float, float]:
        absolute_error = abs(actual - reference)
        if reference != 0.0:
            relative_error = absolute_error / abs(reference)
        else:
            relative_error = 0.0 if actual == 0.0 else math.inf
        limit = absolute_tolerance + relative_tolerance * abs(reference)
        passed = absolute_error <= limit
        if not passed:
            self.changed += 1
        self.maximum_absolute_error = max(self.maximum_absolute_error, absolute_error)
        self.maximum_relative_error = max(
            self.maximum_relative_error, relative_error
        )
        return passed, absolute_error, relative_error, limit


@dataclass
class ComparisonResult:
    reference_count: int = 0
    actual_count: int = 0
    matched: int = 0
    missing: int = 0
    extra: int = 0
    assignment_comments: int = 0
    upper_energy: QuantityComparison = field(default_factory=QuantityComparison)
    lower_energy: QuantityComparison = field(default_factory=QuantityComparison)
    frequency: QuantityComparison = field(default_factory=QuantityComparison)
    line_strength: QuantityComparison = field(default_factory=QuantityComparison)
    einstein_a: QuantityComparison = field(default_factory=QuantityComparison)
    calculated_intensity: QuantityComparison = field(default_factory=QuantityComparison)
    reference_intensity: QuantityComparison = field(default_factory=QuantityComparison)

    @property
    def passed(self) -> bool:
        quantities = (
            self.upper_energy,
            self.lower_energy,
            self.frequency,
            self.line_strength,
            self.einstein_a,
            self.calculated_intensity,
            self.reference_intensity,
        )
        return not (self.missing or self.extra or any(q.changed for q in quantities))


def fortran_float(text: str) -> float:
    """Parse ordinary or Fortran D-exponent floating-point text."""
    return float(text.replace("D", "E").replace("d", "e"))


def twice_quantum_number(text: str, name: str, path: Path, line_number: int) -> int:
    value = fortran_float(text)
    twice_value = round(2.0 * value)
    if not math.isclose(2.0 * value, twice_value, rel_tol=0.0, abs_tol=1.0e-8):
        raise ValueError(
            f"{name} is not an integer or half-integer at "
            f"{path}:{line_number}: {text}"
        )
    return int(twice_value)


def parse_record(fields: list[str], path: Path, line_number: int) -> IntRecord:
    if len(fields) != 23:
        raise ValueError(
            f"expected 23 .int columns at {path}:{line_number}, found {len(fields)}"
        )
    if fields[0] not in ("<-", "->"):
        raise ValueError(
            f"expected '<-' or '->' in column 1 at {path}:{line_number}"
        )

    try:
        numerical_values = [fortran_float(value) for value in fields[6:13]]
        upper_state, upper_v, upper_lambda = map(int, fields[13:16])
        lower_state, lower_v, lower_lambda = map(int, fields[18:21])
        upper_two_sigma = twice_quantum_number(
            fields[16], "upper sigma", path, line_number
        )
        upper_two_omega = twice_quantum_number(
            fields[17], "upper omega", path, line_number
        )
        lower_two_sigma = twice_quantum_number(
            fields[21], "lower sigma", path, line_number
        )
        lower_two_omega = twice_quantum_number(
            fields[22], "lower omega", path, line_number
        )
        record = IntRecord(
            line_number=line_number,
            direction=fields[0],
            upper_two_j=twice_quantum_number(
                fields[1], "upper J", path, line_number
            ),
            upper_symmetry=fields[2],
            lower_two_j=twice_quantum_number(
                fields[3], "lower J", path, line_number
            ),
            lower_symmetry=fields[4],
            branch=fields[5],
            upper_energy=numerical_values[0],
            lower_energy=numerical_values[1],
            frequency=numerical_values[2],
            line_strength=numerical_values[3],
            einstein_a=numerical_values[4],
            calculated_intensity=numerical_values[5],
            reference_intensity=numerical_values[6],
            assignment=(
                upper_state,
                upper_v,
                upper_lambda,
                upper_two_sigma,
                upper_two_omega,
                lower_state,
                lower_v,
                lower_lambda,
                lower_two_sigma,
                lower_two_omega,
            ),
        )
    except ValueError as error:
        raise ValueError(f"invalid .int record at {path}:{line_number}: {error}") from error

    values = (
        record.upper_energy,
        record.lower_energy,
        record.frequency,
        record.line_strength,
        record.einstein_a,
        record.calculated_intensity,
        record.reference_intensity,
    )
    if not all(map(math.isfinite, values)):
        raise ValueError(f"non-finite value at {path}:{line_number}")
    return record


def parse_int(path: Path) -> list[IntRecord]:
    records: list[IntRecord] = []
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line_number, line in enumerate(stream, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith(("#", "!")):
                continue
            records.append(parse_record(stripped.split(), path, line_number))
    if not records:
        raise ValueError(f"no .int records found in {path}")
    return records


def format_group(key: tuple[object, ...]) -> str:
    direction, upper_two_j, upper_symmetry, lower_two_j, lower_symmetry, branch = key
    return (
        f"{direction} J'={upper_two_j / 2:g} {upper_symmetry}, "
        f"J''={lower_two_j / 2:g} {lower_symmetry}, branch={branch}"
    )


def format_assignment(assignment: tuple[int, ...]) -> str:
    us, uv, ul, uts, uto, ls, lv, ll, lts, lto = assignment
    return (
        f"upper=({us},{uv},{ul},{uts / 2:g},{uto / 2:g}), "
        f"lower=({ls},{lv},{ll},{lts / 2:g},{lto / 2:g})"
    )


def energy_distance(reference: IntRecord, actual: IntRecord) -> float:
    return (
        abs(actual.upper_energy - reference.upper_energy)
        + abs(actual.lower_energy - reference.lower_energy)
        + abs(actual.frequency - reference.frequency)
    )


def pair_group(
    reference: list[IntRecord], actual: list[IntRecord]
) -> tuple[list[tuple[IntRecord, IntRecord]], list[IntRecord], list[IntRecord]]:
    """Pair one rigorous group independently of its output line order.

    Exact assignment matches are used first.  If assignments changed, the
    remaining records are paired by their sorted energy positions.  Therefore
    assignment changes can be reported without making the benchmark fail.
    """
    unmatched_reference = list(reference)
    unmatched_actual = list(actual)
    pairs: list[tuple[IntRecord, IntRecord]] = []

    for observed in sorted(actual, key=lambda record: record.energy_sort_key):
        candidates = [
            expected
            for expected in unmatched_reference
            if expected.assignment == observed.assignment
        ]
        if not candidates:
            continue
        expected = min(candidates, key=lambda item: energy_distance(item, observed))
        pairs.append((expected, observed))
        unmatched_reference.remove(expected)
        unmatched_actual.remove(observed)

    unmatched_reference.sort(key=lambda record: record.energy_sort_key)
    unmatched_actual.sort(key=lambda record: record.energy_sort_key)
    fallback_count = min(len(unmatched_reference), len(unmatched_actual))
    pairs.extend(
        zip(
            unmatched_reference[:fallback_count],
            unmatched_actual[:fallback_count],
        )
    )
    return (
        pairs,
        unmatched_reference[fallback_count:],
        unmatched_actual[fallback_count:],
    )


def add_message(messages: list[str], maximum: int, message: str) -> None:
    if len(messages) < maximum:
        messages.append(message)


def compare_quantity(
    result: QuantityComparison,
    name: str,
    reference_value: float,
    actual_value: float,
    absolute_tolerance: float,
    relative_tolerance: float,
    identity: str,
    problems: list[str],
    maximum_report: int,
) -> None:
    passed, _absolute_error, relative_error, limit = result.add(
        reference_value,
        actual_value,
        absolute_tolerance,
        relative_tolerance,
    )
    if not passed:
        add_message(
            problems,
            maximum_report,
            f"{name}: {identity}, reference={reference_value:.9e}, "
            f"actual={actual_value:.9e}, delta={actual_value-reference_value:+.9e}, "
            f"relative error={relative_error:.6g}, limit={limit:.6g}",
        )


def compare_files(
    reference_path: Path,
    actual_path: Path,
    energy_atol: float = 1.0e-3,
    frequency_atol: float = 1.0e-3,
    value_atol: float = 0.0,
    value_rtol: float = 1.0e-6,
    maximum_report: int = 20,
) -> tuple[ComparisonResult, list[str], list[str]]:
    reference_records = parse_int(reference_path)
    actual_records = parse_int(actual_path)
    result = ComparisonResult(
        reference_count=len(reference_records), actual_count=len(actual_records)
    )
    problems: list[str] = []
    comments: list[str] = []

    reference_groups: dict[tuple[object, ...], list[IntRecord]] = defaultdict(list)
    actual_groups: dict[tuple[object, ...], list[IntRecord]] = defaultdict(list)
    for record in reference_records:
        reference_groups[record.rigorous_key].append(record)
    for record in actual_records:
        actual_groups[record.rigorous_key].append(record)

    for key in sorted(reference_groups.keys() | actual_groups.keys(), key=str):
        pairs, missing, extra = pair_group(
            reference_groups.get(key, []), actual_groups.get(key, [])
        )
        result.matched += len(pairs)
        result.missing += len(missing)
        result.extra += len(extra)
        group_text = format_group(key)

        for record in missing:
            add_message(
                problems,
                maximum_report,
                f"MISSING RECORD: {group_text}, reference line={record.line_number}, "
                f"frequency={record.frequency:.6f}",
            )
        for record in extra:
            add_message(
                problems,
                maximum_report,
                f"EXTRA RECORD: {group_text}, actual line={record.line_number}, "
                f"frequency={record.frequency:.6f}",
            )

        for expected, observed in pairs:
            identity = (
                f"{group_text}, reference line={expected.line_number}, "
                f"actual line={observed.line_number}"
            )
            compare_quantity(
                result.upper_energy, "UPPER ENERGY", expected.upper_energy,
                observed.upper_energy, energy_atol, 0.0, identity,
                problems, maximum_report,
            )
            compare_quantity(
                result.lower_energy, "LOWER ENERGY", expected.lower_energy,
                observed.lower_energy, energy_atol, 0.0, identity,
                problems, maximum_report,
            )
            compare_quantity(
                result.frequency, "FREQUENCY", expected.frequency,
                observed.frequency, frequency_atol, 0.0, identity,
                problems, maximum_report,
            )
            compare_quantity(
                result.line_strength, "LINE STRENGTH", expected.line_strength,
                observed.line_strength, value_atol, value_rtol, identity,
                problems, maximum_report,
            )
            compare_quantity(
                result.einstein_a, "EINSTEIN A", expected.einstein_a,
                observed.einstein_a, value_atol, value_rtol, identity,
                problems, maximum_report,
            )
            compare_quantity(
                result.calculated_intensity, "CALCULATED INTENSITY",
                expected.calculated_intensity, observed.calculated_intensity,
                value_atol, value_rtol, identity, problems, maximum_report,
            )
            # This is copied from the fitting data rather than calculated.  An
            # exact numerical comparison catches accidentally using a different
            # reference dataset.
            compare_quantity(
                result.reference_intensity, "REFERENCE INTENSITY",
                expected.reference_intensity, observed.reference_intensity,
                0.0, 0.0, identity, problems, maximum_report,
            )

            if expected.assignment != observed.assignment:
                result.assignment_comments += 1
                add_message(
                    comments,
                    maximum_report,
                    f"ASSIGNMENT COMMENT: {identity}; "
                    f"reference {format_assignment(expected.assignment)}; "
                    f"actual {format_assignment(observed.assignment)}",
                )

    return result, problems, comments


def print_quantity(name: str, quantity: QuantityComparison, units: str = "") -> None:
    unit_text = f" {units}" if units else ""
    print(
        f"{name}: changed={quantity.changed}, "
        f"maximum absolute error={quantity.maximum_absolute_error:.6g}{unit_text}, "
        f"maximum relative error={quantity.maximum_relative_error:.6g}"
    )


def print_report(
    result: ComparisonResult, problems: list[str], comments: list[str]
) -> None:
    print(f"Duo intensity-fitting benchmark: {'PASS' if result.passed else 'FAIL'}")
    print(
        f"Records: reference={result.reference_count}, actual={result.actual_count}, "
        f"matched={result.matched}, missing={result.missing}, extra={result.extra}, "
        f"assignment comments={result.assignment_comments}"
    )
    print_quantity("Upper energy", result.upper_energy, "cm-1")
    print_quantity("Lower energy", result.lower_energy, "cm-1")
    print_quantity("Frequency", result.frequency, "cm-1")
    print_quantity("Line strength", result.line_strength)
    print_quantity("Einstein A", result.einstein_a, "s-1")
    print_quantity("Calculated intensity", result.calculated_intensity)
    print_quantity("Reference intensity", result.reference_intensity)
    for message in problems:
        print(message)
    if comments:
        print("Informative assignment comments (do not affect PASS/FAIL):")
        for message in comments:
            print(message)


def non_negative(text: str) -> float:
    value = float(text)
    if value < 0.0:
        raise argparse.ArgumentTypeError("must not be negative")
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Benchmark two Duo intensity-fitting .int files."
    )
    parser.add_argument("reference_int", type=Path)
    parser.add_argument("actual_int", type=Path)
    parser.add_argument(
        "--energy-atol", type=non_negative, default=1.0e-3,
        help="absolute upper/lower energy tolerance in cm-1 (default: 0.001)",
    )
    parser.add_argument(
        "--frequency-atol", type=non_negative, default=1.0e-3,
        help="absolute transition-frequency tolerance in cm-1 (default: 0.001)",
    )
    parser.add_argument(
        "--value-rtol", type=non_negative, default=1.0e-6,
        help="relative tolerance for strength, A, and calculated intensity (default: 1e-6)",
    )
    parser.add_argument(
        "--value-atol", type=non_negative, default=0.0,
        help="absolute tolerance added for strength, A, and intensity (default: 0)",
    )
    parser.add_argument(
        "--max-report", type=int, default=20,
        help="maximum stored problem and comment messages (default: 20)",
    )
    return parser


def run(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        result, problems, comments = compare_files(
            args.reference_int,
            args.actual_int,
            energy_atol=args.energy_atol,
            frequency_atol=args.frequency_atol,
            value_atol=args.value_atol,
            value_rtol=args.value_rtol,
            maximum_report=max(0, args.max_report),
        )
        print_report(result, problems, comments)
        return 0 if result.passed else 1
    except MemoryError:
        print(
            "ERROR: insufficient memory while reading the .int files. "
            "Reduce the fitting dataset or use a machine with more memory.",
            file=sys.stderr,
        )
        return 2
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(run())
