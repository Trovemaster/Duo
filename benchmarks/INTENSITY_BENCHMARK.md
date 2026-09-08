# Duo intensity benchmarks

`duo_intensity_benchmark.py` compares the `.states` and `.trans` files produced
by a Duo `INTENSITY` calculation.  It uses only the Python standard library.

## Running a comparison

```text
python duo_intensity_benchmark.py reference.states reference.trans new.states new.trans
```

On this Windows installation, Python bundled with Codex can be used directly:

```powershell
& "C:\Users\Sergey\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" `
  .\duo_intensity_benchmark.py `
  "C:\path\reference.states" "C:\path\reference.trans" `
  "C:\path\new.states" "C:\path\new.trans"
```

Exit status `0` means PASS, `1` means a benchmark difference, and `2` means a
file/parsing error or insufficient memory.

## What is compared

For `.states`, records are matched by the state ID in column 1.  The following
are pass/fail quantities:

- energy (column 2), with an absolute tolerance of `0.001 cm-1` by default;
- degeneracy (column 3), exactly;
- J (column 4), exactly as an integer or half-integer;
- `+/-` and `e/f` parities (columns 5 and 6), exactly.

The remaining assignment columns are informative.  Differences are printed as
comments but do not affect PASS/FAIL.

For `.trans`, records are matched by `(upper state ID, lower state ID)`.  Line
order is irrelevant, so neither file needs to be sorted.  Einstein A is tested
with relative tolerance `1e-6` by default, and frequency is tested with absolute
tolerance `0.001 cm-1`.

The test for a numerical value is:

```text
absolute_error <= absolute_tolerance + relative_tolerance * abs(reference)
```

Tolerances can be changed, for example:

```text
python duo_intensity_benchmark.py reference.states reference.trans new.states new.trans --a-rtol 1e-5 --frequency-atol 0.0001
```

The supplied CH `.trans` file prints Einstein A with five digits in its
mantissa.  Consequently, the sixth significant digit is not present in that
file: the default `1e-6` relative tolerance effectively requires the displayed
five-digit values to remain unchanged.  To test six or more actual significant
digits, Duo must first print A with greater precision.

## Large transition files

The reference transitions are loaded into a dictionary keyed by the two state
IDs.  The new `.trans` file is streamed one line at a time, and matched reference
entries are removed immediately.  This is independent of line order and avoids
holding two large transition files in memory.

If Python cannot allocate enough memory, the program stops with an explanatory
message.  Reduce the line-list size using tighter intensity, frequency, energy,
or J thresholds, or run the benchmark on a machine with more memory.
