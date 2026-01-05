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
