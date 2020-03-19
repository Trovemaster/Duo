Specification of curves and couplings (Duo objects)
***************************************************

Once the main global parameters have been specified as described in the
previous sections, it is necessary to introduce the PECs and the various coupling
curves defining the Hamiltonian. Dipole moment curves (DMCs), which are necessary for
calculating spectral line intensities, are also discussed in this section, as well
as some special objects which are used for fitting.
Each object specification consists in a first part in which
keywords are given and a second one (starting from the
``values`` keyword) in which numerical values are given;
the order of the keywords is not important, except for ``values``.
Each object specification is terminated by the ``end`` keyword.

Objects of type ``poten`` (i.e., PECs, discussed in more detail below)
begin with a line of the kind ``poten N``
where N is an integer index number counting over potentials and identifying them.
It is recommended that PECs are numbered progressively as 1,2,3,...,
although this only restriction is that the total number Nmax of PECs
should be not less than the total number of states specified by the keywork ``nstates``.

Most other objects (e.g., ``spin-orbit``) are assumed to be matrix elements
of some operator between electronic wave functions and after
the keyword identifying their type require two integer numbers
specifying the two indexes of the two electronic states involved (bra and ket).
The indexes are the numbers specified after the \texttt{poten} keyword.

Currently Duo supports the following types of objects: ``potential``, ``spinorbit``, ``L2``, ``Lx``, ``spinspin``, ``spinspino``, ``bobrot``, 
``spinrot``, ``diabatic``, ``lambdaopq``, ``lambdap2q``, ``lambdaq``, ``abinitio``, ``brot``, ``dipoletm``.


``poten`` (alias: potential) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Objects of type ``poten`` represent potential energy curves (PECs) and are
the most fundamental objects underlying each calculation.
From the point of view of theory each PEC is the solution of the electronic
Schoedinger equation with clamped nuclei, possibly complemented with the
scalar-relativistic correction and with the 
Born-Oppenheimer Diagonal correction
(also known as adiabatic correction). Approximate PECs can be obtained with
well-known quantum chemistry methods such as Hartree-Fock, coupled cluster theory etc.
Objects of type ``poten`` should always appear before
all other objects as they are used to assign to each electronic states its quantum numbers.
Here is an example for a PEC showing the general structure:   
::

      poten 1
      name "a 3Piu"
      symmetry u
      type  EMO
      lambda 1
      mult   3
      values
      V0          0.82956283449835E+03
      RE          0.13544137530870E+01
      DE          0.50061051451709E+05
      RREF       -0.10000000000000E+01
      PL          0.40000000000000E+01
      PR          0.40000000000000E+01
      NL          0.20000000000000E+01
      NR          0.20000000000000E+01
      B0          0.20320375686486E+01
      B1         -0.92543284427290E-02
      B2          0.00000000000000E+00
      end



L2:  (alias: ``L**2``)
^^^^^^^^^^^^^^^^^^^^^^

These objects represent matrix elements between electronic states of the molecule-fixed
  angular momentum operator :math:`\hat{L}^2 = \hat{L}_x^2 + \hat{L}_y^2 +\hat{L}_z^2`.


L+:   (aliases: ``Lplus``, ``LxLy`` and  ``Lx``) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


It represent matrix elements between electronic states of the molecule-fixed
  angular momentum operator :math:`\hat{L}_+ = \hat{L}_x + i \hat{L}_y` and
  :math:`\hat{L}_x` in the :math:`\Lambda`- and Cartesian-representations, respectively.



``spin-orbit`` and ``spin-orbit-x`` 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

These objects are matrix elements of the Breit-Pauli spin-orbit Hamiltonian
in the :math:`\Lambda`- and Cartesian-representations, respectively.

Example:
::

    spin-orbit  1 3
    name "<0,S=0 (X1Sigma+)|LSY|+1 (a3Pi),S=1> SO1"
    spin   0.0 1.0
    lambda 0 -1
    sigma 0.0 -1.0
    type   grid
    factor sqrt(2)  (1 or i)
    units bohr  cm-1
    values
      2.80     17.500000
      2.90     15.159900
      3.00     12.347700
      3.10      9.050780
      3.20      5.391190
      3.30      1.256660
      3.40     -3.304040
      3.50     -8.104950
      3.60    -12.848400
      3.70    -17.229100
      3.80    -21.049000
      3.90    -24.250400
      4.00    -26.876900
      4.10    -29.014700
      4.20    -30.756100
      4.30    -32.181900
      4.50    -34.335500
      5.00    -37.348300
    end

For the ``spin-orbit-x`` case (:math:`\Lambda`-representation), the value of the matrix elements of the
 :math:`\hat{L}_z` operator nust be specified using the ``<x|Lz|y>`` keyword. 
 This representation is designed to work with e.g., the MOLPRO outputs. 
 For :math:`\Lambda\ne 0`, the diagonal SO-matrix element (e.g. between to :math:`\Pi`-components of :math:`\Lambda=1`) 
 should be specified using the :math:`\langle \Pi_x|LSZ |\Pi_y \rangle` component 
 (e.g. :math:`\langle 1.2 |{\rm LSZ} |1.3 \rangle`).



``spin-spin-p`` and ``spin-spin-o`` 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Parametrised phenomenological spin-spin operator (diagonal and off-diagonal. 

``spin-rot`` 
^^^^^^^^^^^^^^

Matrix elements of the spin-rotational operator .

``bob-rot``   
^^^^^^^^^^^

Alias: ``bobrot``. Specifies the rotational :math:`g` factor (rotational Born-Oppenheimer breakdown term),
which can be interpreted as a position-dependent modification to the rotational mass.

``diabatic``
^^^^^^^^^^^^

Alias: ``diabat``. Non-diagonal coupling of potential energy functions in the diabatic 
representation. 

``lambda-opq``, ``lambda-p2q``, and ``lambda-q``  
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
These objects are three Lambda-doubling objects which correspond to 
  :math:`o^{\rm LD }+p^{\rm LD }+q^{\rm LD }`, :math:`p^{\rm LD }+2q^{\rm LD }`, and :math:`q^{\rm LD }` couplings.

Example:
::

     lambda-p2q  1 1
     name "<X,2Pi|lambda-p2q|X,2Pi>"
     lambda     1 1
     spin   0.5 0.5
     type  BOBLEROY
     factor    1.0
     values
       RE           0.16200000000000E+01
       RREF        -0.10000000000000E+01
       P            0.10000000000000E+01
       NT           0.20000000000000E+01
       B0           0.98500969657331E-01
       B1           0.00000000000000E+00
       B2           0.00000000000000E+00
       BINF         0.00000000000000E+00
     end


``abinitio`` 
^^^^^^^^^^^^
  
Objects of type ``abinitio`` (aliases: ``reference``, ``anchor``) are reference, ``abinitio`` curves which may be specified
during fitting. When they are used they constrain the fit so that the fitted function differs as little as possible from the
`ab initio` (reference). The reference curve is typically obtained by `ab initio` methods.
For any Duo object one can specify a corresponding reference curve as in the following example:
::

     abinitio spin-orbit 1 2
     name "<3.1,S=0,0 (B1pSigma)|LSX|+1 (d3Pig),S=1,1>"
     spin   0.0 1.0
     type   grid
     units bohr cm-1
     values
      2.3        -3.207178925    13.0
      2.4        -3.668814404    24.0
      2.5        -4.010985122    35.0
      2.6        -4.271163495    46.0
      2.7        -4.445721312    47.0
      2.8        -4.468083270    48.0
     end


``dipole``  and ``dipole-x``  
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  
``Dipole`` (aliases: ``dipole-moment``, ``TM``):  Diagonal or transition dipole moment curves (DMCs),  necessary for computing 
(dipole-allowed) transition line intensities and related quantities (Einstein :math:`A` coefficients etc.). 

``dipole-x`` is related to the Cartesian-representation.

At the moment Duo cannot compute electric-quadrupole or magnetic dipole transition line intensities.

.. _quadrupole-curves:

``quadrupole``
^^^^^^^^^^^^^^

The keyword ``quadrupole`` is used to specify transition quadrupole moment curves, which are necessary for computing electric-quadrupole transition line intensities and related quantities. The actual calculation of line strengths requires the ``quadrupole`` keyword in the ``intensity`` section also (:ref:`see here <computing-spectra>`).

Keywords used in the specification of objects 
=============================================

Name and quantum numbers
^^^^^^^^^^^^^^^^^^^^^^^^

This is a list of keywords used to specify various parameters of Duo objects. 

* ``name``: object name.

``name`` is a text label which can be assigned to any object for reference in the output. The string must appear within quotation marks. 
Examples:
::

    name "X 1Sigma+"
    name "<X1Sigma\|HSO\|A3Pi>"


* ``lambda``: The quantum number(s) :math:`\Lambda`.

``Lambda`` specifies the quantum number(s) :math:`\Lambda`, 
i.e. projections of the electronic angular momentum onto the molecular axis, either for one (PECs) or two states (couplings).
It must be an integral number and is allowed to be either positive or negative.
The sign of :math:`\Lambda` is relevant when specifying couplings between degenerate states in the spherical representaion (e.g. ``spin-orbit``)
Examples:
::

   lambda 1
   lambda 0 -1

The last example is relative to a coupling-type object and the two numbers refer to the bra and ket states.

* ``sigma``: Spin-projection.


``sigma`` specifies the quantum number(s) :math:`\Sigma`, i.e. the  projections of the total spin onto the molecular axis, 
either for one (diagonal) or two  states (couplings). These values should be real (:math:`-S\le \Sigma \le S`) and can be half-integral,
where :math:`S` is the total spin. ``sigma`` is currently required for the spin-orbit couplings only.

Example:
::

   sigma 0.5 1.5

where two numbers refer to the bra and ket states.

* ``mult`` (alias: ``multiplicity``): Multiplicity


``mult`` specifies the multiplicity of the electronic state(s), given by :math:`(2S + 1)`, where :math:`S` is the total spin.
It must be an integer number and is an alternative to the ``spin`` keyword. 

Examples:
::


   mult 3
   mult 1 3

The last example is relative to a coupling-type object and the two numbers refer to the bra and ket states.

* ``spin``: Total spin.

The total ``spin`` of the electronic state(s), an integer or half-integer number.
Example:
::

   spin 1.0
   spin 0.5 1.5

The last example is relative to a  coupling-type object and the two numbers refer to the bra and ket states.

* ``symmetry``: State symmetry


This keyword tells Duo if the electronic state has gerade ``g`` or ungerade ``u`` symmetry (only for homonuclear diatomics)
and whether it has positive (``+``) or negative ``-`` parity (only for
:math:`\Sigma` states, i.e. states with :math:`\Lambda=0`, for which it is mandatory).

Examples:
::

    symmetry +

::

    symmetry + u

::

    symmetry g

The keywords ``g``/``u`` or ``+``/``-`` can appear in any order.


Other control keys
^^^^^^^^^^^^^^^^^^


* ``type``: Type of the functional representaion. 

``Type`` defines if the object is given on a grid ``type grid`` or
selects the parametrised analytical function  used for representing the objects
or selects the interpolation type to be used. The function types supported by Duo
are listed in :ref:`functions`.

Examples: 
::

   type grid
   type polynomial
   type morse

In the examples above ``grid`` selects numerical interpolation of values given on a grid,
``polynomial`` selects a polynomial expansion and ``morse`` selects a polynomial expansion in the Morse variable.
See :ref:`functions` for details.


* ``Interpolationtype``: Grid interpolation 


is used only for ``type grid`` and specifies
the method used for the numerical interpolation of the numerical values.
The currently implemented interpolation methods are ``Cubicsplines`` and ``Quinticsplines`` (default).

Example:
::

    Interpolationtype Cubicsplines
    Interpolationtype Quinticsplines


* ``factor``: Scaling factor  

This optional keyword permits to rescale any object by
an arbitrary multiplication factor. At the moment the accepted values are any real number,
the imaginary unit :math:`i`, the square root of two, written as ``sqrt(2)``, or products
of these quantities. To write a product simply leave a space between the factors, but do not
use the ``*`` sign. All factor can have a :math:`\pm` sign.
The default value for ``factor`` is 1. This keyword is useful, for example,
to temporarily zero a certain object without removing it from the input file.

Examples:
::

   factor 1.5
   factor -sqrt(2)
   factor sqrt(2)
   factor 5 i
   factor -2 sqrt(2) i


In the last example the factor is read in as :math:`-2 \sqrt{2} i`.
Note that imaginary factors make sense only in some cases for some coupling terms (in particular, spin-orbit) 
in the Cartesian-representation, see Section~\ref{s:representations}.


* ``units``

This keyword selects the units of measure used for the 
the object in question. Supported units are: ``angstroms``
(default) and ``bohr`` for the bond lengths; ``cm-1`` (default),
``hartree`` (aliases are ``au``, ``a.u.``, and ``Eh``), and ``eV`` (electronvolts)
for energies; ``debye`` (default) and ``ea0`` (i.e., atomic units) for dipoles; units can appear in any order. 

Example:
::

    units angstrom cm-1 (default for poten, spin-orbit, lambda-doubling etc)
    units bohr cm-1
    units debye  (default)
    units ae0 bohr


* ``<x|Lz|y>``, ``<z|Lz|xy>`` (aliases ``<a|Lz|b>`` and ``<1|Lz|2>``)  

This keyword is sometimes needed when specifying coupling curves between electronic states
with :math:`|\Lambda| > 0` in order to resolve ambiguities in the definition of the
degenerate components of each electronic state, see:ref:`representations`.

This keyword specifies the matrix element of the :math:`\hat{L}_z` operator between the degenerate components
of the electronic wave function. 

Examples:
::

    <x|Lz|y>   i  -i
    <z|Lz|xy> -2i  i

These matrix elements are pure imaginary number in the form :math:`\pm |\Lambda | i`.
It is the overall :math:`\pm` sign which Duo needs and cannot be otherwise guessed.
As shown in the examples above, each factor should be written in the form :math:`\pm |\Lambda | i` without any
space or ``*`` sign.



* ``Molpro`` is a single, stand-alone keywrd to trigger the molpro even for `non-x` fields.

Example:


    molpro


* ``morphing`` This keyword is used for fitting and switches on the morphing method. 

* ``ZPE``: Zero-point-energy 

``ZPE`` allows to explicitly input the zero-point energy (ZPE) of the molecule (in cm\ :sup:`-1`). This affects the value printed, as by default
Duo  prints energy of rovibronic levels by subtracting the ZPE. If not specified, the lowest energy of the first :math:`J`-block 
(independent of parity) will be used as appear on the line ``Jlist``.

* ``fit_factor`` 

This factor (:math:`d_{\lambda}`) is used as a part of the reference *ab initio* curves of the ``abinitio`` type which (when given) 
is applied to the corresponding weights assigned to the corresponding values of this object. 
It is different from ``fit_factor`` defined within in :ref:`fitting`.

Example:
::

    abinitio poten 1
    name "A 1Pi"
    type   grid
    lambda 1
    mult   1
    units bohr cm-1
    fit_factor  1e1
    values
      2.00	32841.37010	0.01
      2.20	17837.88960	0.10
      2.40	8785.33147	0.70
      2.60	3648.35520	1.00
      2.70	2107.10737	1.00
      2.80	1073.95670	1.00
      2.90	442.52180	1.00
      3.00	114.94960	1.00
      3.10	0.00000	    1.00
      3.20	48.46120	1.00
      3.30	213.34240	1.00
      3.40	455.16980	1.00
      3.50	739.61170	1.00
      3.60	1038.82620	1.00
      3.70	1332.46170	1.00
      4.00	2059.31119	1.00
      4.50	2619.19233	0.30
      5.00	2682.84741	0.30
      6.00	2554.34992	0.30
      8.00	2524.31106	0.30
      10.00	2561.48269	1.00
      12.00	2575.09861	1.00
    end




Definition of the function or a grid 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


* ``values``  

This keyword starts the subsection containing the numerical
values defining the object. 
For one of the ``type``s corresponding to an analytical function (see :ref:`functions`),
the input between ``values`` and ``end`` contains the values of the parameters of the function.
The input consists in two columns separated by spaces containing (i) a string label
identifying the parameter and (ii) the value of the parameter (a real number).

In case of ``fitting`` (see :ref:`fitting`) a third column should
also be provided; the parameters which are permitted to vary during fitting
must have in the third column the string ``fit`` or, alternatively, the letter ``f``
or the number 1. Any other string or number (for example, the string ``nofit`` or the number 0)
implies the parameter should be kept at its initial value.
In the case of fitting, the keyword ``link``
can be also appear at the end of each the line; this keyword permits to
cross-reference values from different objects and is explained
below in this section.

In the case of objects of type ``grid``, the third column can be also used to specify if the grid point needs to vary. 
The first columns contains the bond length :math:`r_i` and a second with the value of the object.
In the case of object of the ``abinitio`` (``reference``) type and specified as ``grid``
a third column can be used to specify the fitting weights (see :ref:`fitting`).


* ``link``  

This special keyword is used in fitting
to force a set of parameters 
(which may be relative to a different object) to have the same value.
For example, in a typical situation one may want to fit a set of PECs and to constrain their
dissociation (asymptotic) energy to the same value (because they are expected from theory to share the same
dissociation channel).


After the keyword ``link`` one should provide three numbers :math:`i_1`, :math:`i_2`, :math:`i_3` defining the parameter ID, where
:math:`i_1` identifies the object type (e.g. ``poten``, ``spin-orbit``, ``spin-rot`` etc.), 
:math:`i_2` is the object number within the type :math:`i_1` and :math:`i_3` is the parameter number as it appears after ``values``. The ID numbers :math:`i_1, i_2, i_3` 
are specified in the fitting outputs in the form `[i,j,k]`. 

Example of the input:
::

    DE     0.50960000000000E+05   fit     link   1   1   3

Example of the corresponding output
::

    DE     0.50960000000000E+05   [ 1   1   3 ]




.. _representations:

Using ab initio couplings in Duo: Representations of the electronic wave functions
==================================================================================


Quantum chemistry programs generally use real-valued electronic wave functions which transform according to the irreducible representations
of the C:sub:`2v` point group (for heteronuclear diatomics) or of D:math:`2h` (for homonuclear diatomics).
On the other hand Duo internally assumes the electronic wave functions are eigenfunctions of the :math:`\hat{L}_z`
operator, which implies they must be complex valued for :math:`|\Lambda| > 0`. Converting from one representation to the other is simple, as

:math:`|\Lambda\rangle =\frac{1}{\sqrt{2}}\left[\mp |1\rangle - i|2\rangle \right].`

where :math:`1\rangle` and :math:`2\rangle` are two Cartesian components of the electronic wave functions in a quantum chemistry program. 
Duo uses the matrix elements of the :math:`\hat{L}_z` to reconstruct the transformation between two representations: 


The keyword ``<x|Lz|y>`` and ``<z|Lz|xy>`` (aliases ``<a|Lz|b>`` and ``<1|Lz|2>``) is required when specifying coupling curves between electronic states
in the ``MOLPRO`` representation (``spin-orbit-x``, ``Lx`` and ``dipole-x``)  with :math:`|\Lambda| > 0`
in order to resolve ambiguities in the definition of the   degenerate components of each electronic state.
This is also the value of the matrix element of the :math:`\hat{L}_z` operator computed for
the two component spherical harmonic, degenerate functions :math:`|x\rangle` and :math:`|y\rangle` for the :math:`\Pi` states or 
:math:`|z\rangle` and :math:`|xy\rangle` for the :math:`\Delta` states etc. 
The corresponding `<x|Lz|y>` values for both coupled states must be provided.

Examples:
::

     <x|Lz|y>   i  -i

::
     
     <z|Lz|xy> -2i  i

This keyword is required for the couplings of the following types: ``spin-orbit-x``, ``Lx`` and ``dipole-x``. 
The suffix ``-x`` indicates that Duo expects the ``x``-component (non-zero) of the corresponding coupling. 
 
This keyword should appear anywhere in the object section, before the ``values`` keyword. 
::

    spin-orbit-x 1 1
    name "X-X SO term"  
    spin 1.0 1.0
    lambda 2 2
    sigma 1.0 1.0
    units angstrom cm-1
    type polynomial
    factor i
    *<x|Lz|y>  2i 2i*
    values
      f 101.2157	  
    end 



These matrix elements are pure imaginary number in the form :math:`\pm |\Lambda | i`.
It is the overall :math:`\pm` sign which \Duo\ needs and cannot be otherwise guessed.
As shown in the examples above, each factor should be written in the form :math:`\pm |\Lambda | i` without any
space or `*` sign.


