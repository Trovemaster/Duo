#!/usr/bin/env python3
"""Benchmark Duo electric-dipole .states and .trans line-list files.

The implementation uses only the Python standard library.  State records are
matched by state ID.  Transition records are matched by (upper ID, lower ID),
so the order of lines in a .trans file is irrelevant.

Exit status 0 means PASS, 1 means benchmark differences were found, and 2
means that a file could not be read/parsed or there was insufficient memory.
"""

from __future__ import annotations

import argparse
import itertools
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence


@dataclass(frozen=True)
class StateRecord:
    state_id: int
    energy: float
    degeneracy: int
    two_j: int
    parity: str
    rotationless_parity: str
    assignment: tuple[str, ...]


@dataclass
class StateComparison:
    reference_count: int = 0
    actual_count: int = 0
    matched: int = 0
    energy_changed: int = 0
    rigorous_changed: int = 0
    missing: int = 0
    extra: int = 0
    assignment_comments: int = 0
    maximum_energy_error: float = 0.0
    squared_energy_error: float = 0.0

    @property
    def passed(self) -> bool:
        return not (
            self.energy_changed
            or self.rigorous_changed
            or self.missing
            or self.extra
        )

    @property
    def rms_energy_error(self) -> float:
        if not self.matched:
            return 0.0
        return math.sqrt(self.squared_energy_error / self.matched)


@dataclass
class TransitionComparison:
    reference_count: int = 0
    actual_count: int = 0
    matched: int = 0
    a_changed: int = 0
    frequency_changed: int = 0
    missing: int = 0
    extra: int = 0
    maximum_a_relative_error: float = 0.0
    maximum_frequency_error: float = 0.0
    squared_frequency_error: float = 0.0

    @property
    def passed(self) -> bool:
        return not (
            self.a_changed
            or self.frequency_changed
            or self.missing
            or self.extra
        )

    @property
    def rms_frequency_error(self) -> float:
        if not self.matched:
            return 0.0
        return math.sqrt(self.squared_frequency_error / self.matched)


def fortran_float(text: str) -> float:
    """Parse ordinary or Fortran D-exponent floating-point text."""
    return float(text.replace("D", "E").replace("d", "e"))


def parse_two_j(text: str, path: Path, line_number: int) -> int:
    value = fortran_float(text)
    two_j = round(2.0 * value)
    if not math.isclose(2.0 * value, two_j, rel_tol=0.0, abs_tol=1.0e-8):
        raise ValueError(f"non-integer/half-integer J at {path}:{line_number}: {text}")
    return int(two_j)


def parse_states(path: Path) -> dict[int, StateRecord]:
    states: dict[int, StateRecord] = {}
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line_number, line in enumerate(stream, start=1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) < 6:
                raise ValueError(
                    f"expected at least 6 .states columns at {path}:{line_number}"
                )
            try:
                record = StateRecord(
                    state_id=int(fields[0]),
                    energy=fortran_float(fields[1]),
                    degeneracy=int(fields[2]),
                    two_j=parse_two_j(fields[3], path, line_number),
                    parity=fields[4],
                    rotationless_parity=fields[5],
                    assignment=tuple(fields[6:]),
                )
            except ValueError as error:
                raise ValueError(f"invalid .states record at {path}:{line_number}: {error}") from error
            if record.state_id in states:
                raise ValueError(
                    f"duplicate state ID {record.state_id} at {path}:{line_number}"
                )
            states[record.state_id] = record
    if not states:
        raise ValueError(f"no state records found in {path}")
    return states


def transition_key(upper_id: int, lower_id: int) -> int:
    """Pack two non-negative 32-bit IDs into one memory-efficient integer."""
    if upper_id < 0 or lower_id < 0 or upper_id > 0xFFFFFFFF or lower_id > 0xFFFFFFFF:
        raise ValueError("transition IDs must be in the range 0..4294967295")
    return (upper_id << 32) | lower_id


def unpack_transition_key(key: int) -> tuple[int, int]:
    return key >> 32, key & 0xFFFFFFFF


def iter_transitions(path: Path) -> Iterator[tuple[int, int, float, float, int]]:
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for line_number, line in enumerate(stream, start=1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) < 4:
                raise ValueError(
                    f"expected 4 .trans columns at {path}:{line_number}"
                )
            try:
                upper_id = int(fields[0])
                lower_id = int(fields[1])
                einstein_a = fortran_float(fields[2])
                frequency = fortran_float(fields[3])
                transition_key(upper_id, lower_id)  # Validate ID range.
            except ValueError as error:
                raise ValueError(f"invalid .trans record at {path}:{line_number}: {error}") from error
            if not all(map(math.isfinite, (einstein_a, frequency))):
                raise ValueError(f"non-finite value at {path}:{line_number}")
            yield upper_id, lower_id, einstein_a, frequency, line_number


def load_reference_transitions(path: Path) -> dict[int, tuple[float, float]]:
    transitions: dict[int, tuple[float, float]] = {}
    for upper_id, lower_id, einstein_a, frequency, line_number in iter_transitions(path):
        key = transition_key(upper_id, lower_id)
        if key in transitions:
            raise ValueError(
                f"duplicate transition ({upper_id}, {lower_id}) at {path}:{line_number}"
            )
        transitions[key] = (einstein_a, frequency)
    if not transitions:
        raise ValueError(f"no transition records found in {path}")
    return transitions


def within(reference: float, actual: float, atol: float, rtol: float) -> tuple[bool, float]:
    limit = atol + rtol * abs(reference)
    return abs(actual - reference) <= limit, limit


def add_problem(problems: list[str], maximum: int, message: str) -> None:
    if len(problems) < maximum:
        problems.append(message)


def compare_states(
    reference_path: Path,
    actual_path: Path,
    energy_atol: float,
    energy_rtol: float,
    maximum_report: int,
) -> tuple[StateComparison, list[str], list[str]]:
    reference = parse_states(reference_path)
    actual = parse_states(actual_path)
    result = StateComparison(reference_count=len(reference), actual_count=len(actual))
    problems: list[str] = []
    comments: list[str] = []

    for state_id in sorted(reference.keys() | actual.keys()):
        expected = reference.get(state_id)
        observed = actual.get(state_id)
        if expected is None:
            result.extra += 1
            add_problem(problems, maximum_report, f"EXTRA STATE: ID={state_id}")
            continue
        if observed is None:
            result.missing += 1
            add_problem(problems, maximum_report, f"MISSING STATE: ID={state_id}")
            continue

        result.matched += 1
        error = abs(observed.energy - expected.energy)
        result.maximum_energy_error = max(result.maximum_energy_error, error)
        result.squared_energy_error += error * error
        energy_ok, limit = within(expected.energy, observed.energy, energy_atol, energy_rtol)
        if not energy_ok:
            result.energy_changed += 1
            add_problem(
                problems,
                maximum_report,
                f"STATE ENERGY: ID={state_id}, reference={expected.energy:.9f}, "
                f"actual={observed.energy:.9f}, delta={observed.energy-expected.energy:+.9f}, "
                f"limit={limit:.6g} cm-1",
            )

        expected_rigorous = (
            expected.degeneracy,
            expected.two_j,
            expected.parity,
            expected.rotationless_parity,
        )
        observed_rigorous = (
            observed.degeneracy,
            observed.two_j,
            observed.parity,
            observed.rotationless_parity,
        )
        if expected_rigorous != observed_rigorous:
            result.rigorous_changed += 1
            add_problem(
                problems,
                maximum_report,
                f"STATE COLUMNS 3-6: ID={state_id}, "
                f"reference=(g={expected.degeneracy}, J={expected.two_j/2:g}, "
                f"parity={expected.parity}, e/f={expected.rotationless_parity}), "
                f"actual=(g={observed.degeneracy}, J={observed.two_j/2:g}, "
                f"parity={observed.parity}, e/f={observed.rotationless_parity})",
            )

        if expected.assignment != observed.assignment:
            result.assignment_comments += 1
            add_problem(
                comments,
                maximum_report,
                f"STATE COMMENT: ID={state_id}, reference={' '.join(expected.assignment)!r}, "
                f"actual={' '.join(observed.assignment)!r}",
            )

    return result, problems, comments


def compare_transitions(
    reference_path: Path,
    actual_path: Path,
    a_atol: float,
    a_rtol: float,
    frequency_atol: float,
    frequency_rtol: float,
    maximum_report: int,
) -> tuple[TransitionComparison, list[str]]:
    expected = load_reference_transitions(reference_path)
    result = TransitionComparison(reference_count=len(expected))
    problems: list[str] = []

    for upper_id, lower_id, actual_a, actual_frequency, _line_number in iter_transitions(actual_path):
        result.actual_count += 1
        key = transition_key(upper_id, lower_id)
        reference_values = expected.pop(key, None)
        if reference_values is None:
            result.extra += 1
            add_problem(
                problems,
                maximum_report,
                f"EXTRA TRANSITION: upper={upper_id}, lower={lower_id}",
            )
            continue

        reference_a, reference_frequency = reference_values
        result.matched += 1

        a_ok, a_limit = within(reference_a, actual_a, a_atol, a_rtol)
        if reference_a != 0.0:
            a_relative_error = abs(actual_a - reference_a) / abs(reference_a)
        else:
            a_relative_error = 0.0 if actual_a == 0.0 else math.inf
        result.maximum_a_relative_error = max(
            result.maximum_a_relative_error, a_relative_error
        )
        if not a_ok:
            result.a_changed += 1
            add_problem(
                problems,
                maximum_report,
                f"EINSTEIN A: upper={upper_id}, lower={lower_id}, "
                f"reference={reference_a:.9e}, actual={actual_a:.9e}, "
                f"relative error={a_relative_error:.6g}, limit={a_limit:.6g} s-1",
            )

        frequency_error = abs(actual_frequency - reference_frequency)
        result.maximum_frequency_error = max(
            result.maximum_frequency_error, frequency_error
        )
        result.squared_frequency_error += frequency_error * frequency_error
        frequency_ok, frequency_limit = within(
            reference_frequency, actual_frequency, frequency_atol, frequency_rtol
        )
        if not frequency_ok:
            result.frequency_changed += 1
            add_problem(
                problems,
                maximum_report,
                f"FREQUENCY: upper={upper_id}, lower={lower_id}, "
                f"reference={reference_frequency:.9f}, actual={actual_frequency:.9f}, "
                f"delta={actual_frequency-reference_frequency:+.9f}, "
                f"limit={frequency_limit:.6g} cm-1",
            )

    result.missing = len(expected)
    remaining_report_slots = max(0, maximum_report - len(problems))
    for key in itertools.islice(expected, remaining_report_slots):
        upper_id, lower_id = unpack_transition_key(key)
        add_problem(
            problems,
            maximum_report,
            f"MISSING TRANSITION: upper={upper_id}, lower={lower_id}",
        )
    return result, problems


def print_report(
    states: StateComparison,
    state_problems: list[str],
    state_comments: list[str],
    transitions: TransitionComparison,
    transition_problems: list[str],
) -> None:
    passed = states.passed and transitions.passed
    print(f"Duo intensity benchmark: {'PASS' if passed else 'FAIL'}")
    print(
        "States: "
        f"reference={states.reference_count}, actual={states.actual_count}, "
        f"matched={states.matched}, energy changed={states.energy_changed}, "
        f"columns 3-6 changed={states.rigorous_changed}, missing={states.missing}, "
        f"extra={states.extra}, assignment comments={states.assignment_comments}"
    )
    print(
        f"State energy maximum error={states.maximum_energy_error:.6g} cm-1, "
        f"RMS error={states.rms_energy_error:.6g} cm-1"
    )
    print(
        "Transitions: "
        f"reference={transitions.reference_count}, actual={transitions.actual_count}, "
        f"matched={transitions.matched}, A changed={transitions.a_changed}, "
        f"frequency changed={transitions.frequency_changed}, "
        f"missing={transitions.missing}, extra={transitions.extra}"
    )
    print(
        f"Einstein-A maximum relative error={transitions.maximum_a_relative_error:.6g}; "
        f"frequency maximum error={transitions.maximum_frequency_error:.6g} cm-1, "
        f"RMS error={transitions.rms_frequency_error:.6g} cm-1"
    )
    for message in state_problems + transition_problems:
        print(message)
    if state_comments:
        print("Informative state-assignment comments (do not affect PASS/FAIL):")
        for message in state_comments:
            print(message)


def non_negative(text: str) -> float:
    value = float(text)
    if value < 0.0:
        raise argparse.ArgumentTypeError("must not be negative")
    return value


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Benchmark Duo electric-dipole .states and .trans files."
    )
    parser.add_argument("reference_states", type=Path)
    parser.add_argument("reference_trans", type=Path)
    parser.add_argument("actual_states", type=Path)
    parser.add_argument("actual_trans", type=Path)
    parser.add_argument(
        "--energy-atol", type=non_negative, default=1.0e-3,
        help="absolute .states energy tolerance in cm-1 (default: 0.001)",
    )
    parser.add_argument(
        "--energy-rtol", type=non_negative, default=0.0,
        help="relative .states energy tolerance (default: 0)",
    )
    parser.add_argument(
        "--a-rtol", type=non_negative, default=1.0e-6,
        help="relative Einstein-A tolerance, approximately six significant figures (default: 1e-6)",
    )
    parser.add_argument(
        "--a-atol", type=non_negative, default=0.0,
        help="absolute Einstein-A tolerance in s-1 (default: 0)",
    )
    parser.add_argument(
        "--frequency-atol", type=non_negative, default=1.0e-3,
        help="absolute frequency tolerance in cm-1 (default: 0.001)",
    )
    parser.add_argument(
        "--frequency-rtol", type=non_negative, default=0.0,
        help="relative frequency tolerance (default: 0)",
    )
    parser.add_argument(
        "--max-report", type=int, default=20,
        help="maximum stored diagnostic messages per category (default: 20)",
    )
    return parser


def run(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        maximum_report = max(0, args.max_report)
        states, state_problems, state_comments = compare_states(
            args.reference_states,
            args.actual_states,
            args.energy_atol,
            args.energy_rtol,
            maximum_report,
        )
        transitions, transition_problems = compare_transitions(
            args.reference_trans,
            args.actual_trans,
            args.a_atol,
            args.a_rtol,
            args.frequency_atol,
            args.frequency_rtol,
            maximum_report,
        )
        print_report(
            states,
            state_problems,
            state_comments,
            transitions,
            transition_problems,
        )
        return 0 if states.passed and transitions.passed else 1
    except MemoryError:
        print(
            "ERROR: insufficient memory while loading the reference .trans file. "
            "Reduce the transition dataset (for example with tighter intensity, "
            "frequency, energy, or J thresholds) or use a machine with more memory.",
            file=sys.stderr,
        )
        return 2
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(run())
