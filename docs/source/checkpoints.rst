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
