.. _omega-representation:

Omega representation
====================

Duo traditionally works in a Hund’s case (a) electronic basis, :math:`|\mathrm{state},\Lambda,S,\Sigma\rangle`, i.e. the ``Lambda–S`` (or ``vib``) contraction. Duo can now also work in an :math:`\Omega`-based contracted representation, where the electronic + spin–orbit problem is diagonalised at each bond length and the resulting :math:`\Omega`-labelled channels are used for vibrational contraction and for the final rovibronic basis.

The :math:`\Lambda S \rightarrow \Omega` transformation (often called the *state-interacting* method) has been widely used to simplify the treatment of spin–orbit coupling in rovibronic calculations; see, e.g., :cite:p:21BeJoLia.SO`, :cite:p:`19PaEvAn.ai`, :cite:p:`11YuBixx.SO`. The key idea is to diagonalise the electronic Hamiltonian together with the Breit–Pauli spin–orbit Hamiltonian to obtain effective potentials for each spin–orbit component. While this can make the problem look "single-state-like", strict equivalence with the original :math:`\Lambda S` formulation requires
transforming the *full* nuclear-motion Hamiltonian, which introduces spin–orbit-induced non-adiabatic couplings (NACs); see, e.g., :cite:p:`10TaKlKr.KCs`, :cite:p:`24BrDrYu`, :cite:p:`25Brady`.

Transforming to the :math:`\Omega` representation
-------------------------------------------------

A general workflow for building the :math:`\Omega` representation is:

1. Solve the electronic Schrödinger equation to obtain electronic wavefunctions, and construct potential energy curves and spin–orbit coupling curves for the electronic states of interest.
2. Build, at each bond length :math:`r`, the electronic + spin–orbit Hamiltonian matrix  :math:`\mathbf{H}_\Omega(r) = \mathbf{V}(r) + \mathbf{H}_{\rm SO}(r)` in the chosen electronic basis.
3. Diagonalise :math:`\mathbf{H}_\Omega(r)` at each :math:`r` to obtain spin–orbit-decoupled channels and effective potentials labelled by :math:`\Omega`.
4. To achieve exact equivalence with the original :math:`\Lambda S` representation, apply the same *r-dependent* unitary transformation to the remaining parts of the rovibronic Hamiltonian. In particular, transforming radial  derivative operators generates NAC terms that must be included for accurate energies, wavefunctions, and intensities.

The diatom in the :math:`\Lambda S` and :math:`\Omega` representations
----------------------------------------------------------------------

In a Hund’s case (a) basis, the coupled rovibronic Schrödinger equation can be written as
(see also the standard Duo theory in the manual):

.. math::
   :label: eq-diatomic-schrodinger-lambdas

   \left[
   \frac{\hbar^2}{2\mu}\left(-\frac{d^2}{dr^2} + \frac{1}{r^2}\hat{\mathbf{R}}^2\right)
   + \mathbf{V}(r) + \mathbf{H}_{\rm SO}(r)
   \right]\vec{\chi}(r) = E_i\,\vec{\chi}(r),

where :math:`\mu` is the reduced mass, :math:`r` is the internuclear separation, :math:`\hat{\mathbf{R}}` is the nuclear rotational angular momentum operator, :math:`\mathbf{V}` contains diagonal Born–Oppenheimer PECs :cite:`27BoOpxx.diabat`, and :math:`\mathbf{H}_{\rm SO}` contains spin–orbit matrix elements (e.g. from *ab initio* electronic-structure calculations).

State-interacting diagonalisation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The state-interacting :math:`\Omega` transformation is defined by diagonalising :math:`\mathbf{V}+\mathbf{H}_{\rm SO}` at each :math:`r`:

.. math::
   :label: eq-omega-diagonalisation

   \mathbf{U}^\dagger(r)\left(\mathbf{V}(r)+\mathbf{H}_{\rm SO}(r)\right)\mathbf{U}(r)
   = \mathbf{V}_\Omega(r),

where :math:`\mathbf{U}(r)` is an :math:`r`-dependent unitary matrix and :math:`\mathbf{V}_\Omega(r)` is diagonal. This defines the transformation of the electronic basis

.. math::
   :label: eq-omega-basis

   |\mathrm{state},\Lambda,S,\Sigma\rangle \;\rightarrow\; |\mathrm{state},\Omega\rangle.

It is tempting to assume that :math:`\mathbf{V}_\Omega(r)` yields fully decoupled single-channel rovibronic problems. However, to remain formally consistent, the kinetic-energy operators in Eq. :eq:`eq-diatomic-schrodinger-lambdas` must also be transformed. The transformation of radial derivatives introduces NAC terms.

Spin–orbit-induced NACs from the vibrational kinetic energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Upon transforming the radial kinetic-energy operator, one obtains additional terms of the standard adiabatic/NAC form (see diabatisation theory, e.g. :cite:`00Baer.diabat`, :cite:`00BaAlxx.diabat`, :cite:p:`24BrDrYu`, :cite:p:`25Brady`):

.. math::
   :label: eq-so-nacs

   -\frac{\hbar^2}{2\mu}\mathbf{U}^\dagger \frac{d^2}{dr^2}\mathbf{U}
   =
   -\frac{\hbar^2}{2\mu}\left[
     \frac{d^2}{dr^2}
     + \left(\mathbf{W}\frac{d}{dr} + \frac{d}{dr}\mathbf{W}\right)
     + \mathbf{W}^2
   \right],

where

.. math::

   \mathbf{W}(r) = \mathbf{U}^\dagger(r)\,\frac{d\mathbf{U}(r)}{dr}, \qquad W_{ij}(r) = \langle i | \frac{d}{dr} | j \rangle,

is a skew-Hermitian (real antisymmetric) matrix of derivative couplings (the columns of :math:`\mathbf{U}` are the
:math:`\Omega`-eigenvectors, see Eq. :eq:`eq-omega-diagonalisation`). The terms :math:`\mathbf{W}\frac{d}{dr} + \frac{d}{dr}\mathbf{W}`
mediate non-adiabatic transitions between channels of the same :math:`\Omega`, while
:math:`-\frac{\hbar^2}{2\mu}\mathbf{W}^2 = \frac{\hbar^2}{2\mu}\mathbf{W}\mathbf{W}^\dagger` act as additional, DBOC-like corrections:
diagonal for two channels with the same :math:`\Omega`, and with off-diagonal elements when there are three or more channels.
All of these terms are included in Duo.

.. note::
   The relative phases of the :math:`\Omega`-eigenvectors are arbitrary at each :math:`r`. Duo fixes them by requiring
   a positive diagonal component at the first grid point and a positive overlap between the eigenvectors at neighbouring
   grid points, so that :math:`\mathbf{U}(r)` is a smooth function of :math:`r`. The same :math:`\mathbf{U}(r)` is used
   for :math:`\mathbf{W}(r)` and for all other transformed couplings (spin, dipoles etc.), which is required
   for the equivalence with the :math:`\Lambda S` representation.

.. note::
   In the :math:`\Omega` representation, “simplifying” the potential by diagonalisation relocates the physics into  induced non-adiabatic terms. Neglecting these terms can lead to substantial errors in energies and wavefunctions, and therefore also in intensities and lifetimes (see, e.g. :cite:p:`10TaKlKr.KC`, :cite:p:`24BrDrYu`, :cite:p:`26BrYu`).

How to use the Omega representation in Duo
------------------------------------------

Switching from the standard ``Lambda–S`` contraction to the :math:`\Omega` contraction requires only a change in the
``contraction`` block:

::

  contraction
    omega
    vmax 20 40 40
  end

Here the usual ``vib`` (``Lambda–S``) option is replaced by ``omega``.

High-level algorithm used by Duo (omega mode)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``omega`` is selected, Duo performs (conceptually) the following steps:

1. At each radial grid point :math:`r`, build an extended electronic matrix and diagonalise it to obtain :math:`\mathbf{U}(r)` and :math:`\mathbf{V}_\Omega(r)`. Duo then uses :math:`\mathbf{U}(r)` to transform relevant curve-based operators to the :math:`\Omega` representation and to generate the induced NAC terms arising from the  :math:`r`-dependence of :math:`\mathbf{U}(r)`.

   .. note::
      The “extended matrix” used for diagonalisation includes spin–orbit together with other *J*-independent electronic terms available in the model. Any terms not included at this stage are still handled in the full rovibronic Hamiltonian, but the precise partitioning of terms between “diagonalisation” and “final Hamiltonian” affects the  structure of the induced NACs and should be treated consistently.

2. Solve the pure vibrational problems for each :math:`\Omega` channel using the chosen DVR method (sinc DVR, 5-point finite differences or Lobatto DVR, see below), obtaining vibrational functions :math:`|\Omega, v\rangle`.

3. Use these vibrational functions to compute vibrational matrix elements of all required operators and curves (PECs, DMCs, couplings) in the :math:`\Omega` representation.

4. Construct the final rovibronic basis and Hamiltonian. Schematically,

   .. math::

      \Phi_{J,\Omega,v} \;=\; |J,\Omega\rangle\,|\Omega,v\rangle,

   which are then symmetrised (Wang functions) and used to build and diagonalise the rovibronic Hamiltonian.

5. If requested, compute rovibronic intensities and/or line lists using the same workflows as in the ``Lambda–S`` mode.

.. note::
   The vibrational basis functions :math:`|\Omega, v\rangle` are solutions of the single-channel problems with the
   :math:`\Omega`-PECs only, while the NAC and DBOC-like terms are added to the rovibronic Hamiltonian. For very sharp
   avoided crossings (strongly localised NACs) the convergence with respect to ``vmax`` can therefore be slow and should be checked.


Grids and DVR methods: sinc and Lobatto
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :math:`\Omega` representation can be used with all solution methods of Duo, including the Lobatto DVR
(Gauss–Lobatto quadrature grid, ``type 6``), e.g.

::

   SOLUTIONMETHOD LOBATTO

   grid
     npoints 1000
     range   0.80 4.00
     type    6
   end

   contraction
     omega
     vmax 80
   end

The derivative couplings :math:`\mathbf{W}(r)=\mathbf{U}^\dagger d\mathbf{U}/dr` require :math:`d\mathbf{U}/dr` on the grid.
By default it is evaluated with 9-point finite-difference formulas using the Fornberg weights for the actual
(in general non-uniformly spaced, e.g. Lobatto) grid points; for a uniform grid these are the standard central
(one-sided at the grid edges) formulas. For the Lobatto grid, the spectral derivative of the Lagrange interpolating polynomial
through all Lobatto points (the same polynomial representation as used by the Lobatto DVR) can be requested instead:

::

   contraction
     omega
     vmax 80
     derivatives spectral
   end

The spectral derivative converges very fast for smooth, well resolved :math:`\mathbf{U}(r)`, but, being global, it spreads the
errors of any under-resolved sharp features of :math:`\mathbf{U}(r)` (e.g. a near crossing of two weakly coupled states of the same
:math:`\Omega`) over the whole grid; the local finite-difference formulas (``derivatives FD``, default) are more robust.

The matrix of :math:`d/dr` used for the NAC terms is the antisymmetrised Lobatto-DVR derivative matrix
:math:`\langle u_i|d/dr|u_j\rangle` for ``LOBATTO`` and the sinc-DVR derivative matrix for ``SINC``.
As for the ``Lambda–S`` contraction, ``LOBATTO`` requires the grid ``type 6`` and an even number of grid points.

Omega-specific output
---------------------

The output is similar to the standard ``Lambda–S`` mode, but includes additional diagnostics for curves and operators transformed to the :math:`\Omega` representation. For example, Duo can print PECs labelled by :math:`\Omega`:

::

   PECS in the Omega representation
   #  #    Omega Lambda  Sigma  State
   1  1     -1.0    0     -1.0 a3Sigma-
   2  1      0.0    0      0.0 X1Sig+
   ...

and transformed multi-component operators such as the electronic spin operator (schematically shown here):

::

   S = (S+,S-) in the Omega representation
   #   i Omega State Lambda Sigma <-> j Omega State Lambda Sigma
   1   1  -1.0   2    0    -1.0      1   0.0   1    0     0.0
   ...

Similarly, dipole moment curves can be printed in the :math:`\Omega` representation:

::

   Dipole in the Omega representation
   #    Omega State Lambda Sigma   Omega  State Lambda Sigma
   1      0.0   1    0     0.0       0.0   1    0     0.0
   ...

Warning and words of caution
----------------------------

.. warning::
   The :math:`\Omega` representation in Duo is currently a work in progress. It has been tested on a limited set of    systems and couplings (see :cite:p:`26BrYu`) and may produce incorrect results for models involving additional  couplings or corrections not yet fully validated in ``omega`` mode. If you use ``omega`` mode beyond the tested scope, it is strongly recommended to validate against the standard ``Lambda–S`` calculation.

Even when technically available, the fully decoupled (single-state) :math:`\Omega` approximation should not be used blindly. Transforming to the :math:`\Omega` representation unavoidably introduces induced NAC terms which may be essential for accurate energies and transition properties, especially for spin-forbidden bands. The broader lesson is that unitary “simplifications” must be accompanied by an accounting of the physics moved elsewhere in the Hamiltonian (see discussion and examples in :cite:p:`26BrYu`).

Keywords
--------

.. glossary::
   :sorted:

   omega
      Option in the ``contraction`` block that requests contraction and rovibronic calculations in the       :math:`\Omega` representation. This triggers an :math:`r`-dependent diagonalisation to define :math:`\Omega`  channels and includes the induced non-adiabatic couplings required for formal consistency.

      Example:
      ::

         contraction
           omega
           vmax 20 40 40
         end

   vib
      Standard option in the ``contraction`` block for the Hund’s case (a) ``Lambda–S`` representation  (the default behaviour in Duo). In this manual it may also be referred to as the ``Lambda–S`` contraction.

   derivatives
      Option in the ``contraction`` block (``omega`` only) defining how :math:`d\mathbf{U}/dr` for the NACs is computed:
      ``FD`` (default, 9-point finite differences on the actual grid) or ``SPECTRAL`` (Lagrange-interpolation derivative,
      ``SOLUTIONMETHOD LOBATTO`` only).

      Example:
      ::

         contraction
           vib
           vmax 20 40 40
         end

References
----------

.. bibliography:: 

