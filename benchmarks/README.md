# Duo output benchmarks

`duo_benchmark.py` compares the numerical energies in a new Duo `.out` file
with a trusted reference result.  It deliberately ignores banners, timings,
memory reports, echoed input, and other text that varies between machines.

The program uses only the Python standard library and works with Python 3.9 or
newer.  It checks:

- rovibronic eigenvalues on lines containing `||`;
- contracted vibrational energies on lines containing `[ state v ]`.

Final rovibronic records are matched by the rigorous `(J, parity, i)` identity,
where `i` is Duo's running number within the J/parity block.  Assignment quantum
numbers (`state`, `v`, `lambda`, `spin`, `sigma`, `omega`, and the state name)
are not used for matching and do not determine PASS/FAIL.  If an assignment
changes, the program prints a `COMMENT` showing the changed fields.

Contracted vibrational records are matched by electronic state, `v`, and state
name.  Missing and additional identities cause a clear failure rather than
shifting every subsequent comparison.

## Quick start

Run Duo in the usual way, then compare the new output directly with a trusted
sample output:

```text
duo.exe < project.inp > new.out
python duo_benchmark.py compare reference.out new.out
```

This repository includes a CSV snapshot made from the supplied PS output:

```text
python duo_benchmark.py compare benchmarks/PS/reference.csv new.out
```

A passing comparison returns exit status 0.  A numerical or structural
difference returns 1, which makes the command suitable for CI.  A file or
parsing error returns 2.

## Tolerances

The default absolute tolerance is `0.001 cm-1` for both kinds of energy.  This
corresponds to comparing approximately three digits after the decimal point,
without rounding the printed numbers first.  The limits can be changed:

```text
python duo_benchmark.py compare reference.out new.out --final-atol 0.0001 --vibrational-atol 0.001
```

The test for every energy is:

```text
absolute_error <= absolute_tolerance + rtol * abs(reference_energy)
```

For spectroscopic energies, an absolute tolerance has a more direct physical
meaning than a relative tolerance, so `--rtol` defaults to zero.

## Creating or updating a benchmark

Keep the full trusted `.out` file if convenient, or extract just the supported
results into a readable CSV file:

```text
python duo_benchmark.py snapshot trusted.out benchmarks/molecule/reference.csv
```

The CSV retains Duo's running number as the `number` column.  Review a changed
reference before committing it.  Updating the snapshot merely
because a test failed would hide a real numerical regression.

## Running the Python tests

No test framework needs to be installed:

```text
python -m unittest discover -v
```

## Extending the parser

The two regular expressions near the top of `duo_benchmark.py` describe the
supported Duo lines.  `Record.key` defines the identity used to match a level,
while `Record.assignment` defines non-rigorous diagnostic fields.  A further
output section can be supported by adding a parser that
creates another `Record` kind, then assigning its tolerance in `compare()`.
