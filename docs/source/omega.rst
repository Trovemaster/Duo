.. _omega-representation:

Omega representation
====================

Duo traditionally works in a Hund’s case (a) electronic basis, :math:`|\mathrm{state},\Lambda,S,\Sigma\rangle`, i.e. the ``Lambda–S`` (or ``vib``) contraction. Duo can now also work in an :math:`\Omega`-based contracted representation, where the electronic + spin–orbit problem is diagonalised at each bond length and the resulting :math:`\Omega`-labelled channels are used for vibrational contraction and for the final rovibronic basis.

The :math:`\Lambda S \rightarrow \Omega` transformation (often called the *state-interacting* method) has been widely used to simplify the treatment of spin–orbit coupling in rovibronic calculations; see, e.g., :cite:p:`21BeJoLia.SO`, :cite:p:`19PaEvAn.ai`, :cite:p:`11YuBixx.SO`. The key idea is to diagonalise the electronic Hamiltonian together with the Breit–Pauli spin–orbit Hamiltonian to obtain effective potentials for each spin–orbit component. While this can make the problem look "single-state-like", strict equivalence with the original :math:`\Lambda S` formulation requires
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
   The relative phases of the :math:`\Omega`-eigenvectors are arbitrary at each :math:`r`; Duo fixes them so that
   :math:`\mathbf{U}(r)` is a smooth function of :math:`r` and consistent with the parity of the rovibronic states,
   see :ref:`omega-phases`. The same :math:`\mathbf{U}(r)` is used for :math:`\mathbf{W}(r)` and for all other
   transformed couplings (spin, dipoles etc.), which is required for the equivalence with the :math:`\Lambda S`
   representation.

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

Here the usual ``vib`` (``Lambda–S``) option is replaced by ``omega``. All other sections of the input (states, grid,
potentials, couplings, dipoles, intensities, checkpoints) are the same as for the ``Lambda–S`` contraction: the curves and
couplings are always given in the :math:`\Lambda S` representation and are transformed to the :math:`\Omega` representation
internally.

Quick start
^^^^^^^^^^^

A complete input for a toy :math:`X\,^3\Sigma^-`/:math:`b\,^1\Sigma^+` model coupled by spin–orbit, using the Lobatto DVR and
the coupled :math:`\Omega` vibrational basis (see :ref:`omega-coupled`), is:

::

   atoms O S

   states X b

   jrot 0 - 5

   symmetry Cs(M)

   SOLUTIONMETHOD LOBATTO      (optional, sinc DVR by default)

   grid
     npoints 1000              (even number of points for LOBATTO)
     range   0.80 4.00
     type    6                 (Lobatto grid; type 0 for the sinc DVR)
   end

   contraction
     omega coupled             (or just omega for the uncoupled basis)
     vmax 30
   end

   poten X
   name "X3Sigma-"
   symmetry -
   lambda 0
   mult 3
   type   Morse
   values
     VE    0.0
     RE    1.21
     B0    1.90
     DE    42056.0
   end

   poten b
   name "b1Sigma+"
   symmetry +
   lambda 0
   mult 1
   type   Morse
   values
     VE    5000.0
     RE    1.24
     B0    2.3
     DE    27000.0
   end

   spin-orbit X b
   name "<X|LSZ|b>"
   spin  1 0
   sigma 0 0
   lambda 0 0
   factor i
   type  polynom
   values
     V0  400.0
     RE  1.5
     B1  300.0
   end

Replacing ``omega coupled`` by ``vib`` gives the same model in the ``Lambda–S`` representation. With a complete vibrational basis
(``vmax`` equal to ``npoints``) the two representations give the same energies, which is a useful check of a new model.
In this example the two potential curves cross near 1.63 Å (about 12750 cm\ :sup:`-1`), where the NACs between the two
:math:`\Omega=0` channels are localised: with ``vmax 30`` the ``coupled`` basis reproduces the converged energies below
20000 cm\ :sup:`-1` to :math:`3\times 10^{-5}` cm\ :sup:`-1`, while the ``uncoupled`` basis (``omega``) with the same ``vmax``
is in error by up to 9 cm\ :sup:`-1`.

High-level algorithm used by Duo (omega mode)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``omega`` is selected, Duo performs (conceptually) the following steps:

1. At each radial grid point :math:`r`, build an extended electronic matrix and diagonalise it to obtain :math:`\mathbf{U}(r)` and :math:`\mathbf{V}_\Omega(r)`. Duo then uses :math:`\mathbf{U}(r)` to transform relevant curve-based operators to the :math:`\Omega` representation and to generate the induced NAC terms arising from the  :math:`r`-dependence of :math:`\mathbf{U}(r)`.

   .. note::
      The “extended matrix” used for diagonalisation includes spin–orbit together with other *J*-independent electronic terms available in the model. Any terms not included at this stage are still handled in the full rovibronic Hamiltonian, but the precise partitioning of terms between “diagonalisation” and “final Hamiltonian” affects the  structure of the induced NACs and should be treated consistently.

2. Solve the pure vibrational (:math:`J=0`) problems for each :math:`\Omega` using the chosen DVR method (sinc DVR, 5-point finite differences or Lobatto DVR, see below), obtaining vibrational functions :math:`|\Omega, v\rangle`. By default (``uncoupled``) these are solutions of the single-channel problems for each :math:`\Omega` channel; with the ``coupled`` option they are solutions of the coupled-channel problem for each :math:`\Omega`, see :ref:`omega-coupled` below.

3. Use these vibrational functions to compute vibrational matrix elements of all required operators and curves (PECs, DMCs, couplings) in the :math:`\Omega` representation.

4. Construct the final rovibronic basis and Hamiltonian. Schematically,

   .. math::

      \Phi_{J,\Omega,v} \;=\; |J,\Omega\rangle\,|\Omega,v\rangle,

   which are then symmetrised (Wang functions) and used to build and diagonalise the rovibronic Hamiltonian.

5. If requested, compute rovibronic intensities and/or line lists using the same workflows as in the ``Lambda–S`` mode.

.. _omega-coupled:

Uncoupled and coupled vibrational basis sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Two choices of the vibrational (:math:`J=0`) basis set are available in the :math:`\Omega` representation.

``uncoupled`` (default)
   The vibrational functions :math:`|\Omega,i,v\rangle` are solutions of the single-channel problems with the
   adiabatic :math:`\Omega`-PECs :math:`V^{\Omega}_i(r)` only. The NAC and DBOC-like terms are added to the rovibronic
   Hamiltonian together with all other couplings. This basis is compact and well labelled when the :math:`\Omega`
   representation is the natural one (smooth, weak NACs), but for sharp avoided crossings (strongly localised NACs, cusps of
   the adiabatic PECs) the convergence with respect to ``vmax`` can be very slow.

``coupled``
   For each :math:`\Omega` the :math:`J`-independent coupled-channel problem

   .. math::

      \hat{H}_\Omega = -\frac{\hbar^2}{2\mu}\frac{d^2}{dr^2} + \mathbf{V}^{\Omega}(r)
      -\frac{\hbar^2}{2\mu}\left(\mathbf{W}\frac{d}{dr}+\frac{d}{dr}\mathbf{W}\right)+\frac{\hbar^2}{2\mu}\mathbf{W}\mathbf{W}^{T}

   is solved on the grid for all channels of this :math:`\Omega` simultaneously (a dense matrix of dimension
   :math:`N_\Omega\times N_{\rm grid}`, where :math:`N_\Omega` is the number of channels). The NACs are represented
   with the same grid derivative matrix as used for the rovibronic Hamiltonian (sinc, 5-point or Lobatto DVR).
   The lowest eigenvectors, i.e. multi-channel vibronic functions :math:`|\Omega,k\rangle`, form the vibrational basis set;
   their number for each :math:`\Omega` is the sum of ``vmax`` of its channels (the same basis size as for ``uncoupled``),
   optionally limited by ``enermax``. The NAC and DBOC-like terms are then exactly included at :math:`J=0` within each
   :math:`\Omega` block and are not added to the rovibronic Hamiltonian again; the remaining (:math:`J`-dependent) couplings
   between the :math:`\Omega` blocks are evaluated with the multi-channel functions. Each function is labelled by its dominant
   channel (``state``, :math:`\Lambda`, :math:`\Sigma`) and by the vibrational quantum number counted within this channel.
   This option is recommended for sharp avoided crossings, long-range and predissociation problems.

::

   contraction
     omega coupled
     vmax 20
   end

(``coupled`` or ``uncoupled`` can also be given on a separate line of the ``contraction`` block.)

As an example, for the :math:`X\,^3\Sigma^-`/:math:`b\,^1\Sigma^+` model of SO from :cite:p:`26BrYu`, with a sharp crossing at 1.63 Å
(:math:`\approx` 12800 cm\ :sup:`-1`, adiabatic gap 153 cm\ :sup:`-1`), the errors of the :math:`J=0` energies relative to the
converged ``Lambda-S`` values on the same grid (sinc, 1001 points) are:

============  ==================  ====================  ====================  ============================
Level/cm-1    ``vib``, vmax=20    ``omega``, vmax=20    ``omega``, vmax=80    ``omega coupled``, vmax=20
============  ==================  ====================  ====================  ============================
9170          2e-4                0.19                  0.06                  < 1e-6
10024         7e-4                1.2                   0.40                  < 1e-6
11910         4e-4                22                    8.1                   1e-5
12606         8e-3                61                    23                    6e-5
============  ==================  ====================  ====================  ============================

**Choosing the basis and checking the convergence.**

* For models where the :math:`\Omega` channels are only weakly coupled by smooth NACs, the ``uncoupled`` basis is compact and
  gives well defined :math:`(\Omega, v)` labels.
* For sharp avoided crossings (e.g. a weakly coupled crossing of two potentials of the same :math:`\Omega`), long-range and
  predissociation problems, use ``coupled``.
* In both cases check the convergence by increasing ``vmax`` and, for a new model, by comparing with the ``Lambda–S``
  (``vib``) calculation on the same grid.
* The cost of ``coupled`` is one dense diagonalisation of dimension :math:`N_\Omega\times N_{\rm grid}` for each :math:`\Omega`
  (only the lowest eigenvalues are computed). For an :math:`\Omega=0` block with channels of both :math:`\sigma_v` parities
  (see :ref:`omega-phases`), the two parities are solved separately.

**Output.** With ``verbose 4`` or higher, the contracted energies of the ``coupled`` basis set are printed as
::

   Vibrational (contracted) energies, coupled Omega basis set:
       i        Energy/cm    Channel  v  Omega   Weight  State

       1          0.000000 [    1   0     0.0 ]   1.0000  X3Sigma-
       ...

where ``Channel`` is the dominant :math:`\Omega` channel, ``v`` is the vibrational quantum number counted within this channel and
``Weight`` is the weight of the dominant channel in the multi-channel function (1 for a pure channel). The final rovibronic states
are assigned using the same labels.

**Checkpoints.** With ``eigenfunc save`` in the ``checkpoint`` section, each contracted function of the ``coupled`` basis set is
written to the ``*_vib.chk`` file as :math:`N_{\rm channels}\times N_{\rm points}` values (the channel components one after
another, zeros for the missing channels of a given :math:`\Omega`), and ``Nchannels`` is given in the ``Npoints`` line, see
:ref:`vib-chk`. The reduced densities (``density save`` and the bound/unbound filters, see :ref:`unbound-states`) are summed over the
channels.


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

.. _omega-phases:

Phases and parity
^^^^^^^^^^^^^^^^^

The eigenvectors (columns of :math:`\mathbf{U}(r)`) are defined up to a sign at each grid point. Duo fixes the signs as follows:

* At the first grid point the diagonal component of each eigenvector is made positive (or its largest component, if the diagonal
  one vanishes); at all other grid points the sign is chosen to give a positive overlap with the eigenvector at the previous
  point. :math:`\mathbf{U}(r)` is therefore smooth also when its components change sign, e.g. when a spin–orbit coupling curve
  crosses zero.

* The parity-adapted (Wang) combinations of the :math:`+\Omega` and :math:`-\Omega` functions,
  :math:`\frac{1}{\sqrt{2}}\left(|J,\Omega,i\rangle\pm(-1)^{J-\Omega}|J,-\Omega,i\rangle\right)`, require the channel :math:`i` of
  the :math:`-\Omega` block to be the :math:`\sigma_v` image of the channel :math:`i` of the :math:`+\Omega` block. Duo therefore
  defines the eigenvectors of each :math:`-\Omega` block as

  .. math::

     U^{(-\Omega)}_{n',i} = (-1)^{e_n}\, U^{(+\Omega)}_{n,i},\qquad e_n=-\Lambda+S-\Sigma\;(+1 \mbox{ for } \Sigma^- \mbox{ states}),

  where :math:`n'=\sigma_v(n)` is the :math:`(\mathrm{state},-\Lambda,-\Sigma)` partner of the :math:`\Lambda S` basis function
  :math:`n`. This is the same phase convention as used for the parity of the :math:`\Lambda S` basis functions and for the
  :math:`\sigma_v` images of the couplings. The NACs and all other :math:`J`-independent couplings are then identical in the
  :math:`\pm\Omega` blocks, and the parities agree with the ``Lambda–S`` calculation for all states (e.g. also for
  :math:`^2\Sigma^-`, :math:`^3\Pi` and :math:`^4\Pi` states).

* The :math:`\Omega=0` block is mapped onto itself by :math:`\sigma_v`. It is diagonalised in the :math:`\sigma_v`-adapted basis
  :math:`(|n\rangle\pm(-1)^{e_n}|n'\rangle)/\sqrt{2}` (:math:`|n\rangle` for :math:`\Lambda=\Sigma=0`), separately for the two
  :math:`\sigma_v` parities :math:`p=\pm1`, so that each :math:`\Omega=0` channel has a well defined parity and the parity of
  :math:`|J,0,i\rangle` is :math:`(-1)^J p_i`. This is needed e.g. for :math:`^3\Pi_{0^\pm}`, :math:`^1\Sigma^-` and the
  :math:`\Sigma=0` component of :math:`^3\Sigma^+`. The :math:`\Omega=0` channels with :math:`p=+1` are listed first in the output.

Duo prints a warning if an electronic :math:`-\Omega` matrix is not the :math:`\sigma_v` image of the :math:`+\Omega` one, or if the
:math:`\Omega=0` block is not :math:`\sigma_v`-symmetric (this would indicate inconsistent phases of the input couplings).


Supported couplings and current limitations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

=======================================================  ===========================================================
Term                                                     Status in the ``omega`` contraction
=======================================================  ===========================================================
PECs, spin-orbit (diagonal and non-diagonal)             included in the :math:`\Omega` diagonalisation; validated
NACs and DBOC-like terms (:math:`r`-dependence of U)     included (``uncoupled`` or ``coupled``); validated
S-uncoupling (rotation-spin)                             included; validated
L-uncoupling (``Lx``/``L+``)                             included; validated
Spin-rotation                                            included; validated
BOB-rot                                                  included; validated
Electric dipole intensities                              included; validated
Bound/unbound filters, reduced densities, checkpoints    included
:math:`L^2`                                              included (as for ``vib``), not separately validated
Spin-spin, BOB-vib                                       included in the diagonalisation, not validated
:math:`\Lambda`-doubling (``p2q``, ``q``)                not implemented (Duo stops)
Magnetic dipole and quadrupole intensities, fitting      not validated
=======================================================  ===========================================================

"Validated" means that the :math:`\Omega` results agree with the ``Lambda–S`` results on the same grid for test models with these
terms (to about :math:`10^{-5}` cm\ :sup:`-1` with complete vibrational basis sets). The ``overlap`` option is not available in the
``omega`` contraction. With ``SOLUTIONMETHOD 5POINTDIFFERENCES`` the first-derivative matrix used for the NACs is only of second order,
which leaves small differences from the ``Lambda–S`` results; the sinc and Lobatto DVRs are recommended.

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
           omega
           vmax 80
           derivatives spectral
         end

   coupled
      Option in the ``contraction`` block (``omega`` only): the vibrational basis set is obtained from the :math:`J=0`
      coupled-channel problem for each :math:`\Omega` (adiabatic :math:`\Omega`-PECs + NACs + DBOC), see
      :ref:`omega-coupled`. It can be given as ``omega coupled`` or on a separate line.

      Example:
      ::

         contraction
           omega coupled
           vmax 20
         end

   uncoupled
      Option in the ``contraction`` block (``omega`` only, default): the vibrational basis set is obtained from the
      single-channel problems for each :math:`\Omega` channel, see :ref:`omega-coupled`.

References
----------

.. bibliography:: 

