Checkpointing (wavefunctions, moments, and reduced density)
==========================================================
.. _checkpointing:

Duo can *checkpoint* (save) and later *reuse* the key numerical objects required for post-processing:
rovibronic eigenvalues and eigenvectors, vibrational basis functions on the radial grid, and (optionally)
matrix elements of transition-moment operators (electric dipole, magnetic dipole, electric quadrupole, etc.).
This enables intensity calculations **without recomputing** the eigenproblem, basis functions, or (previously
computed) transition-moment integrals.

The checkpointing workflow is controlled by the ``checkpoint`` input block.

Overview of the ``checkpoint`` block
------------------------------------

Basic structure:
::

  checkpoint
    eigenfunc <save|read>
    dipole   <save|read|calc>
    filename <prefix>
  end

The keyword ``filename`` defines the checkpoint prefix (for example ``YO_01``), which is used to name the
produced files.

What gets written: checkpoint files
-----------------------------------

When checkpointing is enabled, Duo creates up to four files:

* ``<prefix>_values.chk`` (ASCII)
    Rovibronic eigenvalues together with descriptors (e.g., quantum numbers, state labels, properties).
    The **first ~10 lines** form a *fingerprint* of the calculation (molecule, masses, basis size, grid, state list, etc.).
    When reading checkpoints, Duo compares this fingerprint with the current input to avoid mixing incompatible files.

* ``<prefix>_vectors.chk`` (binary)
    Rovibronic eigenvectors for all computed solutions.

* ``<prefix>_vib.chk`` (ASCII)
    Vibrational (contracted) basis functions tabulated on the radial grid.

* ``external.chk`` (binary)
    Stored matrix elements of transition moments (electric/magnetic dipole, quadrupole, etc.), written when
    ``dipole save`` is used.  (Despite the generic name, this file belongs to the checkpoint dataset and should be kept
    together with the other ``<prefix>_*`` files.)

.. note::
   The exact set of operators stored in ``external.chk`` depends on what transition-moment curves are present in the input
   and enabled in the calculation (e.g., electric dipole, magnetic dipole, quadrupole).

Writing checkpoints for later intensity work
--------------------------------------------

A typical setup to compute energies/eigenvectors (and optionally moment matrix elements) and write them to disk is:
::

  checkpoint
    eigenfunc save
    dipole   save
    filename YO_01
  end

This writes ``YO_01_values.chk``, ``YO_01_vectors.chk``, ``YO_01_vib.chk``, and ``external.chk``.

Controlling the J-range that is saved
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The range of :math:`J` values to be computed (and thus stored) can be defined either in the ``intensity`` block or via
the global ``jrot`` line.

Example using the ``intensity`` block:
::

  intensity
    absorption
    states-only
    thresh_coeff  1e-60
    thresh_dipole 1e-9
    temperature   300.0
    J,  0.5, 6.5
    freq-window  -0.001,  20000.0
    energy low   -0.001,  3000.00, upper   -0.00, 23000.0
  end

Here, ``states-only`` tells Duo *not* to generate intensities during this run. The calculation still produces and stores
all rovibronic states for :math:`J = 0.5 \dots 6.5` (as requested), which can then be reused later.

Alternatively, you can define the rotational range globally:
::

  jrot 0.5  -  6.5

.. tip::
   Using ``states-only`` together with checkpointing is a convenient way to generate and store states once, and perform
   intensity calculations later (potentially split into multiple runs).

Reading checkpoints to compute intensities (no recomputation)
------------------------------------------------------------

Once the checkpoint data exist, Duo can compute intensities **directly from the saved information**.

Use:
::

  checkpoint
    eigenfunc read
    dipole   read
    filename YO_01
  end

and provide an ``intensity`` block for the transitions you want:
::

  intensity
    absorption
    thresh_coeff  1e-60
    thresh_dipole 1e-9
    temperature   300.0
    Qstat  337.0515
    J,  5.5, 6.5
    freq-window  -0.001,  20000.0
    energy low   -0.001,  3000.00, upper   -0.00, 23000.0
  end

In this example, Duo reads the stored eigenvalues/eigenvectors/basis and transition-moment matrix elements and computes
intensities only for :math:`J = 5.5 \dots 6.5`.

Main advantage: split and parallelise intensity production
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because intensities can be computed for *subsets* of :math:`J` (or different frequency/energy windows) from the same
checkpoint dataset, large spectra can be split into independent jobs, e.g.

* compute and store states once for a broad :math:`J` range;
* run multiple intensity jobs in parallel, each handling different :math:`J` intervals or spectral windows;
* run “test” intensity calculations with different thresholds/windows without re-solving the Schrödinger equation.

Recomputing moment matrix elements with new dipole curves (``dipole calc``)
--------------------------------------------------------------------------

Another common use case is to keep eigenfunctions fixed while changing the transition-moment curves
(e.g. refining an electric dipole moment function to experimental intensities).

In this mode Duo reads the stored wavefunctions/basis from checkpoint files, **recomputes the transition-moment matrix
elements**, and writes them (so that intensities can then be generated consistently with the updated moments):
::

  checkpoint
    eigenfunc read
    dipole   calc
    filename YO_01
  end

You would then run an intensity calculation using the updated matrix elements (either in the same run, or in a
subsequent ``dipole read`` run, depending on your workflow).

Example: structure of ``*_values.chk`` (eigenvalue checkpoint file)
-------------------------------------------------------------------

Below is an excerpt from ``KH_01_values.chk`` (KH project):
::

  Molecule = K               H
  masses   =      38.963706486430      1.007825032230
  Nroots   =       20
  Nbasis   =      160
  Nestates =        4
  Npoints   =      501
  range   =      0.9000000    16.0000000
  nJs     =       11
  Jrange  =   0.00000 10.00000
  X1Sigma+, A1Sigma+, B1Pi, C1Sigma+,    <- States
           140  <- Nlevels
           140  <- Ndimen
             1    492.2960321790      0.0    1    1   X1Sigma+       0  0      0.0      0.0      0.0    1       1  T     2.20122 0.11673E-29
             2   1448.1762331140      0.0    1    1   X1Sigma+       1  0      0.0      0.0      0.0    1       1  T     2.26446 0.97793E-30
             3   2374.1127556956      0.0    1    1   X1Sigma+       2  0      0.0      0.0      0.0    1       1  T     2.32626 0.32747E-29
             4   3271.0297066643      0.0    1    1   X1Sigma+       3  0      0.0      0.0      0.0    1       1  T  ...
  ...
  End of eigenvalues

The header (roughly the first 10 lines) is used as a **fingerprint** of the project. Duo checks it against the current
input to ensure that checkpoint files are not accidentally reused across different models, parameter sets, grids, or
molecules.

Reduced density checkpoint (optional)
-------------------------------------

In addition to wavefunction/moment checkpointing, Duo can compute and store the vibrational reduced density on the radial
grid. Enable it via:
::

  checkpoint
    density save
    filename <prefix>
  end

The reduced density :math:`\rho(r)` is computed as:
.. math::

   \rho(r) =
   \sum_{v,v'} \sum_{\mathrm{State},\Lambda,\Sigma}
   \left[C_{v,\mathrm{State},\Lambda,\Sigma}^{(i)}\right]^*
   C_{v',\mathrm{State},\Lambda,\Sigma}^{(i)}
   \,\phi_v(r)^* \phi_{v'}(r)\,\Delta r .

A typical density checkpoint record looks like:
::

  0.545190480438E-08 ||      1.5  0       1
  0.286121234769E-07 ||      1.5  0       1
  0.134835397210E-06 ||      1.5  0       1
  ...

The first column is the density value at grid point :math:`r_i`, followed by delimiter ``||``, total angular momentum
:math:`J`, parity :math:`\tau`, and the state number (as in the Duo output).

Notes on terminology and legacy keywords
----------------------------------------

Older inputs and parts of the manual may refer to ``eigenvectors save`` or use default prefixes such as ``eigen_*``.
In the current syntax, the recommended form is:

* ``eigenfunc save`` / ``eigenfunc read`` for wavefunction checkpointing
* ``dipole save`` / ``dipole read`` / ``dipole calc`` for moment matrix elements
* ``filename <prefix>`` to control output names

.. rubric:: Footnote

.. [#1] Strictly speaking, :math:`\mathbf{J} = \mathbf{R} + \mathbf{L} + \mathbf{S}`
   is the sum of the rotational and total electronic angular momenta; it is the total angular momentum only if
   the nuclear angular momentum :math:`\mathbf{I}` is zero (or neglected).
