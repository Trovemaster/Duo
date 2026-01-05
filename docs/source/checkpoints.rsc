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

