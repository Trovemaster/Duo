# Duo intensity-fitting `.int` benchmark

`duo_intensity_fit_benchmark.py` compares a new Duo intensity-fitting `.int`
file with a trusted reference.  It uses only the Python standard library and
works with Python 3.9 or newer.

## Run Duo and compare the result

In the input, the `FITTING` block contains both `intensity` and an `output`
name.  Duo writes the selected calculated and reference intensities to:

```text
<output-name>.int
```

The attached input currently contains:

```text
output O2-66_quad_f27d_J60_analytic_m03_morph_damp_i4_qzz_analyt_all-states
```

It therefore creates
`O2-66_quad_f27d_J60_analytic_m03_morph_damp_i4_qzz_analyt_all-states.int`.
This does not overwrite the trusted `O2_quad_int-fit_check1.int` reference.

From Windows PowerShell, run Duo as follows. Calling `setvars.bat` makes the
Intel runtime DLLs available to an ifx-built executable:

```powershell
cd C:\sergei\programs\duo-inputs\O2
cmd /c 'call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat" intel64 >nul && C:\sergei\programs\Duo\build\ifx\duo.exe < O2_quad_int-fit_check1.inp > O2_quad_int-fit_check1_current.out'
```

Then compare the trusted and newly generated files with:

```powershell
python "C:\Users\Sergey\Documents\ChatGPT\Duo\duo_intensity_fit_benchmark.py" `
  "O2_quad_int-fit_check1.int" `
  "O2-66_quad_f27d_J60_analytic_m03_morph_damp_i4_qzz_analyt_all-states.int"
```

If `python` is not installed or not on `PATH`, the Python runtime bundled with
Codex can be used directly:

```powershell
& "C:\Users\Sergey\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe" `
  "C:\Users\Sergey\Documents\ChatGPT\Duo\duo_intensity_fit_benchmark.py" `
  "O2_quad_int-fit_check1.int" `
  "O2-66_quad_f27d_J60_analytic_m03_morph_damp_i4_qzz_analyt_all-states.int"
```

Alternatively, change `output` in a working copy of the input to a shorter
name such as `O2_quad_int-fit_check1_current`; the comparison script accepts
any filenames.

Exit status `0` means PASS, `1` means a benchmark difference, and `2` means a
file/parsing error or insufficient memory.

## What the 23 columns mean

The script reads the `.int` records written by Duo as:

1. transition direction;
2. upper J and symmetry;
3. lower J and symmetry;
4. branch label;
5. upper energy, lower energy, and transition frequency;
6. calculated line strength, Einstein A, and absorption intensity;
7. the reference intensity copied from the `FITTING` data;
8. upper and lower assignment quantum numbers.

Records are grouped by direction, upper/lower J and symmetry, and branch.  Line
order is irrelevant.  Exact assignment matches are used to pair records within
each group.  If assignments change, unmatched records are paired by energy and
the assignment differences are printed as informative comments; they do not
affect PASS/FAIL.

The pass/fail comparisons are:

- rigorous transition labels and number of records: exact;
- upper/lower energies: `0.001 cm-1` absolute tolerance;
- transition frequency: `0.001 cm-1` absolute tolerance;
- line strength, Einstein A, and calculated intensity: `1e-6` relative
  tolerance;
- reference intensity from the input data: exact numerical match.

The numerical rule for calculated intensity-like values is:

```text
absolute_error <= value_atol + value_rtol * abs(reference)
```

Tolerances can be changed without editing the script:

```powershell
python duo_intensity_fit_benchmark.py reference.int current.int `
  --energy-atol 0.0001 --frequency-atol 0.0001 --value-rtol 1e-5
```

For reference values equal to zero, the default absolute value tolerance is
zero.  If tiny numerical noise around zero should be allowed, supply, for
example, `--value-atol 1e-40`.

## Maintaining the reference

Keep `O2_quad_int-fit_check1.int` unchanged as the trusted reference.  When an
intentional physics change affects the result, inspect the reported transitions
before replacing the reference.  Do not update the reference solely to make a
failed test pass.

The script loads both `.int` files into memory because records must be paired
independently of output order.  If memory allocation fails, it stops with a
message suggesting a smaller fitting dataset or a machine with more memory.
