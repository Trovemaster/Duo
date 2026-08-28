.. _duo-benchmarks:

Benchmarking Duo calculations
=============================

Duo includes lightweight Python programs for checking numerical results after
changes to the Fortran source, compiler, libraries, or build options.  The
programs are in the `benchmarks directory
<https://github.com/Trovemaster/Duo/tree/master/benchmarks>`_ and use only the
Python standard library.  Python 3.9 or newer is recommended.

A benchmark compares a new calculation with a trusted reference calculation.
It deliberately checks spectroscopic results rather than the complete textual
output: compiler banners, timings, memory reports, echoed input, and similar
machine-dependent text are ignored.

Three kinds of calculation are currently supported:

.. list-table:: Benchmark programs
   :header-rows: 1
   :widths: 23 25 25 27

   * - Calculation
     - Program
     - Files compared
     - Principal identity
   * - Rovibronic energies
     - ``duo_benchmark.py``
     - ``.out`` or reference ``.csv`` against ``.out``
     - :math:`J`, parity, and running number within the block
   * - Line-list intensities
     - ``duo_intensity_benchmark.py``
     - ``.states`` and ``.trans``
     - State ID, or upper and lower state IDs
   * - Selected fitting intensities
     - ``duo_intensity_fit_benchmark.py``
     - ``.int``
     - Direction, upper/lower :math:`J` and symmetry, and branch

Run the commands below from the root of a Duo source checkout.  All comparison
programs return exit status ``0`` for PASS, ``1`` for a numerical or structural
difference, and ``2`` for a file, parsing, or memory error.  They can therefore
be used directly in continuous-integration scripts.

Numerical tolerances
--------------------

Unless stated otherwise, numerical values are accepted when

.. math::

   |x_{\rm new}-x_{\rm ref}| \leq
   \epsilon_{\rm abs}+\epsilon_{\rm rel}|x_{\rm ref}|.

Absolute energy tolerances are expressed in cm\ :sup:`-1`.  A relative
tolerance is more appropriate for intensities and Einstein coefficients,
which can span many orders of magnitude.  Every tolerance can be changed with
a command-line option; use ``python <program> --help`` to list the available
options.

Energy benchmarks
-----------------

Run Duo in the normal way and keep the new main output separately from the
trusted output:

.. code-block:: console

   duo.exe < project.inp > current.out
   python benchmarks/duo_benchmark.py compare reference.out current.out

The comparison extracts two kinds of record from the main output:

* final rovibronic energies on lines containing ``||``;
* contracted vibrational energies on lines containing ``[ state v ]``.

Final rovibronic levels are matched by :math:`(J,p,i)`, where ``i`` is Duo's
running number within the same :math:`J` and parity block.  This running number
is used instead of approximate state assignments because assignments can swap
after very small numerical changes.  Changes to electronic-state, vibrational,
:math:`\Lambda`, spin, :math:`\Sigma`, :math:`\Omega`, or state-name assignments
are reported as comments and do not determine PASS or FAIL.

Contracted vibrational energies are matched by electronic state, vibrational
quantum number, and state name.  Missing or additional records always fail the
benchmark.  The default absolute tolerance is ``0.001 cm-1`` for both final
and contracted energies, and the default relative tolerance is zero.

A compact CSV reference can be created when keeping a complete ``.out`` file
is inconvenient:

.. code-block:: console

   python benchmarks/duo_benchmark.py snapshot trusted.out reference.csv
   python benchmarks/duo_benchmark.py compare reference.csv current.out

The CSV snapshot retains the running number used in the final-state identity.
For example, the repository contains PS references under ``benchmarks/PS``.

Line-list intensity benchmarks
------------------------------

An ``INTENSITY`` section containing ``linelist <name>`` produces
``<name>.states`` and ``<name>.trans``.  Compare a new pair with a trusted pair
using

.. code-block:: console

   python benchmarks/duo_intensity_benchmark.py \
       reference.states reference.trans \
       current.states current.trans

The format of these files is independent of whether the transition operator
is electric dipole, magnetic dipole, or electric quadrupole, so the same
comparison procedure can be used for E1, M1, and E2 line lists.

States file
^^^^^^^^^^^

Records in ``.states`` are matched by the state ID in column 1.  The following
columns determine PASS or FAIL:

* energy in column 2, with default absolute tolerance ``0.001 cm-1``;
* total degeneracy in column 3, exactly;
* :math:`J` in column 4, exactly as an integer or half-integer;
* total parity and rotationless ``e/f`` parity in columns 5 and 6, exactly.

The remaining state-assignment columns are informative.  Differences in them
are printed as comments but do not fail the benchmark.

Transitions file
^^^^^^^^^^^^^^^^

Records in ``.trans`` are matched by ``(upper state ID, lower state ID)``.
Their order in either file is irrelevant.  Einstein :math:`A` is compared with
relative tolerance ``1e-6`` by default, while the transition wavenumber uses
an absolute tolerance of ``0.001 cm-1``.

The reference transitions are stored in a dictionary and the new file is read
one line at a time.  This avoids holding two large transition files in memory.
If the reference is too large for available memory, the program exits with an
explanation; reduce the calculation using tighter intensity, frequency,
energy, or :math:`J` limits, or use a machine with more memory.

The repository contains a CH E1 example under ``benchmarks/CH``.

Intensity calculations selected by ``FITTING``
-----------------------------------------------

Duo can combine an ``INTENSITY`` calculation with a ``FITTING`` transition
list to calculate intensities only for selected reference transitions.  Put
the ``intensity`` keyword in the ``FITTING`` section and specify an ``output``
name.  Setting ``itmax 0`` performs a straight-through calculation without
changing fitted parameters.  Duo writes ``<output>.int``.

Compare this file with a trusted reference using

.. code-block:: console

   python benchmarks/duo_intensity_fit_benchmark.py \
       reference.int current.int

Each ``.int`` record contains direction, upper and lower :math:`J` and
symmetry, branch, upper and lower energies, transition wavenumber, line
strength, Einstein :math:`A`, calculated absorption intensity, the reference
intensity copied from the fitting data, and upper/lower assignment quantum
numbers.

Records are first grouped by the rigorous direction, :math:`J`, symmetry, and
branch labels.  Exact assignment matches are used to pair records inside a
group.  If assignments change, the remaining records are paired by energy and
the assignment changes are reported only as comments.  Output order therefore
does not affect the result.

The default comparisons are:

* upper and lower energies: absolute tolerance ``0.001 cm-1``;
* transition wavenumber: absolute tolerance ``0.001 cm-1``;
* line strength, Einstein :math:`A`, and calculated intensity: relative
  tolerance ``1e-6``;
* reference intensity copied from the input: exact numerical match.

For reference values equal to zero, relative tolerance alone permits no
nonzero result.  A small absolute floor can be supplied explicitly, for
example

.. code-block:: console

   python benchmarks/duo_intensity_fit_benchmark.py \
       reference.int current.int --value-atol 1e-40

Changing tolerances
-------------------

The following examples make an energy comparison stricter and an intensity
comparison less strict, respectively:

.. code-block:: console

   python benchmarks/duo_benchmark.py compare reference.out current.out \
       --final-atol 0.0001 --vibrational-atol 0.001

   python benchmarks/duo_intensity_benchmark.py \
       reference.states reference.trans current.states current.trans \
       --a-rtol 1e-5 --frequency-atol 0.0001

On Windows PowerShell, use the backtick instead of a backslash to continue a
command onto the next line, or put the complete command on one line.

Running Duo on Windows
----------------------

Input redirection with ``<`` is supported by Command Prompt.  From PowerShell,
invoke it through ``cmd /c``:

.. code-block:: powershell

   cmd /c ".\duo.exe < project.inp > current.out"

Use ``python`` or ``python3`` according to the local Python installation.  A
full path to ``python.exe`` can also be used when Python is not on ``PATH``.

Maintaining trusted references
------------------------------

Reference results should be produced with a reviewed Duo version and retained
unchanged during ordinary regression testing.  When a deliberate physics or
format change alters a benchmark:

#. inspect every reported structural and numerical difference;
#. confirm that the change is expected on more than one build when practical;
#. create the replacement reference from the reviewed executable;
#. commit the implementation, benchmark input, reference, and explanation
   together.

Do not replace a reference merely because a test failed: doing so can conceal
a numerical regression.  Compiler- or architecture-specific numerical noise
should normally be handled with a justified tolerance rather than separate
unreviewed references.
