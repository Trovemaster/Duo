Nuclear hyperfine structure
===========================

Provided nuclear hyperfine interaction curves,
Duo calculates field-free hyperfine structure of diatomic molecules
for the following three cases:

1. only one of the nuclei possesses non-zero spin, e.g., 14N16O and 24MgH.
2. two identical nuclei with non-zero spin, e.g., H2 and D2. 
3. two different nuclei, both with non-zero spin, e.g., 14N17O and AlF.

To determine which case and coupling scheme to use, 
the nuclei and their isotopes must be specified using keyword ``ATOMS``. 

.. Nuclear electric quadrupole interaction
   and nuclear magnetic dipole interactions
   including Fermi-contact, nuclear spin-electron spin dipole-dipole,
   nuclear spin-orbit, nuclear spin-rotation, nuclear spin-nuclear spin dipole-dipole,
   can be involved in the calculation.

Input files
-------------------

The hyperfine calculation for case 1 or 2 is turned on 
with the following section:
::
    
    hyperfine 
        I 1
    end

where the value after the keyword ``I``  specifies 
the spin of the only non-zero-spin nucleus for case 1, 
or the common spin of both identical nuclei for case 2. 

The hyperfine calculation for case 3 can be turned on 
using either of the following two input formats:
::

    hyperfine
        I 1 0.5
    end

or

::

    hyperfine
        I1 1
        I2 0.5
    end

where the keywords ``I1`` and ``I2`` refer to 
the first and second nuclei, respectively,
and the two values specify the nuclear spins in the order given by ``ATOMS``. 

The following rules apply to nuclear-spin inputs:

* Nuclear spins must be finite, non-negative integers or half-integers.
* An isotope with zero nuclear spin in the database cannot be assigned non-zero value.
* For case 1 or 2, supply only single value or none. 
* For case 3, supply either both spins or neither.
  To neglect one nuclear spin, explicitly set its value to ``0``.
* If no nuclear spins are supplied, Duo uses the values from the internal isotope database.

The artificial model with nuclear spins different from the database values 
may lead to unexpected results.
Such calculations should be interpreted with care.

Hyperfine interaction curves must be supplied, 
which are typically from quantum chemistry calculations or fits to experimental data. 
For details on specifying these curves,
see :ref:`hyperfine-coupling-curves`.

Currently, Duo does not support refinement of hyperfine curves.
Thus, ``fitting`` and ``hyperfine`` sections cannot be activated simultaneously.

You may turn off the hyperfine calculation by:
::

    hyperfine none
        I 1
    end

Total angular momentum
^^^^^^^^^^^^^^^^^^^^^^^^

The global setup of ``J``,
::

    jrot 0.5 - 3.5

affects the maximum of ``F``.
:math:`F_\mathrm{max} = J_\mathrm{max}-I_\mathrm{max}`,
where :math:`I_\mathrm{max}` is:

* the spin of the active nucleus for case 1;
* twice the common nuclear spin for case 2;
* :math:`I_1 + I_2` for case 3.

The minimum of ``F`` is always 0 or 1/2.

Transition moment curves
^^^^^^^^^^^^^^^^^^^^^^^^

For hyperfine-resolved transitions, 
the electric dipole (E1) transition is supported for all three cases,
while the electric quadrupole (E2) transition is supported only for cases 2 and 3. 
The magnetic dipole (M1) transition is not yet supported.

The electric dipole or electric quadrupole moment curves 
for transition strength calculations
should be defined using the same formats
described in `dipole and dipole-x <fields.html#dipole-and-dipole-x>`_
and :ref:`quadrupole curves`, respectively.

Output files
------------

Two output files are generated after calculation.
The one named ``hyperfine.states`` contains hyperfine-resolved states, 
with 12 columns: 

1. State ID.
2. Energy [cm-1].
3. Degeneracy, :math:`2F+1`.
4. Total angular momentum, :math:`F`.
5. Symmetry label: parity (+/-) for cases 1 and 3, and parity with g/u symmetry for case 2.
6. Nuclear spin :math:`I` (case 1), total nuclear spin :math:`I` (case 2), or :math:`F_1` (case 3).
7. Angular momentum excluding nuclear spin, :math:`J`.
8. Electronic state.
9. Vibrational quantum number, :math:`v`.
10. :math:`\Lambda`.
11. :math:`\Sigma`.
12. :math:`\Omega`.

The other named  ``hyperfine.trans`` contains the hyperfine-resolved transitions, 
with 4 columns:

1. State ID of the upper state.
2. State ID of the lower state.
3. Einstein A coefficient [s-1].
4. Transition wavenumber [cm-1].

.. 5.line strength.
   Line strengths in this file are calculated by

   :math:`S = |\langle \psi(m,\tau, F)||T^{(1)}(\mu)||\psi_(m', \tau', F') \rangle|^2`

   where :math:`\psi(m,\tau, F)` and :math:`\psi_(m', \tau', F')` are the wavefunctions
   of the hyperfine states. :math:`T^{(1)}(\mu)` is the tensor of transition
   electric dipole moment.
   Line strengths have the unit of :math:`Debye^2`.

.. An electric dipole moment curve should be defined in a ``dipole`` section in the first place.
   ::

       dipole  1 1
       name "<X,2Pi|DMC|X,2Pi>"
       spin   0.5 0.5
       lambda  1  1
       factor   1   (0, 1 or i)
       type polynom
       values
           A0           1
       end


.. _hyperfine-coupling-curves:

Hyperfine coupling curves
-------------------------

The nuclear hyperfine interaction curves are introduced with the following ten keywords:

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Keyword
     - Interaction
   * - ``hfcc-a``
     - Nuclear spin-orbit
   * - ``hfcc-bF``
     - Fermi-contact
   * - ``hfcc-c``
     - Electron spin-nuclear spin dipole-dipole, :math:`c(R)`
   * - ``hfcc-d``
     - Electron spin-nuclear spin dipole-dipole, :math:`d(R)`
   * - ``hfcc-ci``
     - Nuclear spin-rotation, :math:`c_I(R)`
   * - ``hfcc-cip``
     - Off-diagonal nuclear spin-rotation, :math:`c_I'(R)`
   * - ``hfcc-eqq0``
     - Nuclear electric quadrupole, :math:`eQq_0(R)`
   * - ``hfcc-eqq2``
     - Nuclear electric quadrupole, :math:`eQq_2(R)`
   * - ``hfcc-ii``
     - Dipolar nuclear spin-spin
   * - ``hfcc-iis``
     - Scalar nuclear spin-spin

For the first eight keywords, the suffix ``-1`` is optional: 
for example, ``hfcc-bF`` and ``hfcc-bF-1`` are equivalent. 

For case 3, use the suffix ``-2`` for the second nucleus, 
e.g., ``hfcc-bF-2``, following the order specified in ``ATOMS``.

For case 1, use the keywords without suffix or with ``-1``,
even if the nucleus with non-zero spin is listed second in ``ATOMS``.
e.g., for ``ATOMS Mg-24 H-1``, 
use ``hfcc-bF`` or ``hfcc-bF-1`` for H, not ``hfcc-bF-2``.

The default units are ``angstrom`` and ``cm-1``.
Hyperfine couplings between electronic states are not allowed at present.


``hfcc-a``  
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Nuclear spin - orbit interaction curve is defined by an ``hfcc-a`` section.
This term is important when the electron orbital angular momentum is non-zero.
::

    hfcc-a 1 1
    name "<X2Pi|NSO|X2Pi>" (Nuclear spin-orbit)
    spin 0.5
    factor 1.0
    type polynom
    values 
        A0           0.1
    end


``hfcc-bF``  
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Fermi-contact interaction curve :math:`b_F(R)` is defined by an ``hfcc-bF`` section.
This term is important when the electron spin angular momentum is non-zero.
The :math:`b` constant defined by Frosch and Foley can be 
calculated by 
:math:`b = b_F - c/3`.
::

    hfcc-bf 1 1
    name "<X2Pi|FC|X2Pi>" (Fermi-contact)
    spin 0.5
    factor 1.0
    type polynom
    factor 1 
    values 
        A0           0.1
    end


``hfcc-c`` and   ``hfcc-d``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The electron spin - nuclear spin dipole-dipole
interaction curves :math:`c(R)` and :math:`d(R)` are defined by ``hfcc-c`` and ``hfcc-d`` sections.
::

    hfcc-c 1 1
    name "<X2Pi|SDND_C|X2Pi>" (Electron spin - nuclear spin dipole-dipole, c)
    spin 0.5
    factor 1.0
    type polynom
    values 
        A0           0.1
    end

    hfcc-d 1 1
    name "<X2Pi|SDND_D|X2Pi>" (Electron spin - nuclear spin dipole-dipole, d)
    spin 0.5
    factor 1.0
    type polynom
    values 
        A0           0.1
    end


``hfcc-eqq0`` and ``hfcc-eqq2``  
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The nuclear electric quadrupole
interaction curves :math:`eQq_0(R)` and :math:`eQq_2(R)` are defined by ``hfcc-eqq0`` and ``hfcc-eqq2`` sections.
These terms are active when 
the nuclear spin is not less than 1.
::

    hfcc-eqq0 1 1
    name "<X2Pi|eQq0|X2Pi>" (Electric quadrupole eQq0)
    spin 0.5
    factor 1.0
    type polynom
    values 
        A0           0.1
    end

    hfcc-eqq2 1 1
    name "<X2Pi|eQq2|X2Pi>" (Electric quadrupole eQq2)
    spin 0.5
    factor 1.0
    type polynom
    values 
        A0           0.1
    end


``hfcc-ci`` and ``hfcc-cip``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The nuclear spin - rotation curves 
:math:`c_I(R)` and :math:`c_I'(R)` are defined by ``hfcc-ci`` and ``hfcc-cip`` sections.
These two terms are usually negligible compared with other nuclear hyperfine interactions.
Nevertheless,
when all the other hyperfine couplings are inactive,
this term becomes important, e.g. 
for :math:`^1\Sigma` state of nuclear spin 1/2 molecule.
::

    hfcc-ci 1 1
    name "<X2Pi|NSR|X2Pi>" (Nuclear spin - rotation)
    spin 0.5
    factor 1.0
    type polynom
    values 
        A0           0.1
    end

    hfcc-cip 1 1
    name "<X2Pi|NSR-prime|X2Pi>" (Off-diagonal nuclear spin-rotation)
    spin 0.5
    factor 1.0
    type polynom
    values
        A0           0.1
    end

``hfcc-ii``  
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The nuclear spin-nuclear spin dipole-dipole interaction curve is defined by an ``hfcc-ii`` section, 
with aliases: ``hfcc-ii-d`` and ``hfcc-ii-dipole``. 
This term is available only for cases 2 and 3.
::

    hfcc-ii 1
    name "<X1Sigma+|II|X1Sigma+>" (Nuclear spin-nuclear spin dipole-dipole)
    factor 1.0
    type polynom
    values 
        A0           0.1
    end 

``hfcc-iis``
^^^^^^^^^^^^
The scalar nuclear spin-spin interaction curve is defined
by an ``hfcc-iis`` section, with alias ``hfcc-ii-scalar``.
This term is available only for cases 2 and 3.
::

    hfcc-iis 1 1
    name "<X1Sigma+|IIS|X1Sigma+>" (Scalar nuclear spin-spin)
    spin 0
    lambda 0
    factor 1.0
    type polynom
    values
        A0           0.1
    end
