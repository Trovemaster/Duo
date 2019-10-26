Keywords used in the specification of objects 
---------------------------------------------


This is a list of keywords used to specify various parameters of Duo objects. 

* `name`: object name.

`name` is a text label which can be assigned to any object for reference in the output. The string must appear within quotation marks. 
Examples:
::

    name "X 1Sigma+"
    name "<X1Sigma\|HSO\|A3Pi>"


* `lambda`: The quantum number(s) :math:`\Lambda`.

`Lambda` specifies the quantum number(s) :math:`\Lambda`, 
i.e. projections of the electronic angular momentum onto the molecular axis, either for one (PECs) or two states (couplings).
It must be an integral number and is allowed to be either positive or negative.
The sign of :math:`\Lambda` is relevant when specifying couplings between degenerate states in the spherical representaion (e.g. `spin-orbit`)
Examples:
::

   lambda 1
   lambda 0 -1

The last example is relative to a coupling-type object and the two numbers refer to the bra and ket states.

* `sigma`: Spin-projection.


`sigma` specifies the quantum number(s) :math:`\Sigma`, i.e. the  projections of the total spin onto the molecular axis, 
either for one (diagonal) or two  states (couplings). These values should be real (:math:`-S\le \Sigma \le S`) and can be half-integral,
where :math:`S` is the total spin. `sigma` is currently required for the spin-orbit couplings only.

Example:
::

   sigma 0.5 1.5

where two numbers refer to the bra and ket states.

* `mult` (alias: `multiplicity`): Multiplicity


`mult` specifies the multiplicity of the electronic state(s), given by :math:`(2S + 1)`, where :math:`S` is the total spin.
It must be an integer number and is an alternative to the ``spin`` keyword. 

Examples:
::


   mult 3
   mult 1 3

The last example is relative to a coupling-type object and the two numbers refer to the bra and ket states.

* `spin`: Total spin.

The total `spin` of the electronic state(s), an integer or half-integer number.
Example:
::

   spin 1.0
   spin 0.5 1.5

The last example is relative to a  coupling-type object and the two numbers refer to the bra and ket states.

* `symmetry`: State symmetry


This keyword tells Duo if the electronic state has gerade `g` or ungerade `u` symmetry (only for homonuclear diatomics)
and whether it has positive (`+`) or negative `-` parity (only for
:math:`\Sigma` states, i.e. states with :math:`\Lambda=0`, for which it is mandatory).

Examples:
::

    symmetry +

::

    symmetry + u

::

    symmetry g

The `g`/`u` or `+`/`-` can appear in any order.

* `type`: Type of the functional representaion. 

`Type` defines if the object is given on a grid `type grid` or
selects the parametrised analytical function  used for representing the objects
or selects the interpolation type to be used. The function types supported by Duo
are listed in :ref:`functions`.

Examples: 
::

   type grid
   type polynomial
   type morse

In the examples above `grid` selects numerical interpolation of values given on a grid,
`polynomial` selects a polynomial expansion and `morse` selects a polynomial expansion in the Morse variable.
See :ref:`functions` for details.


* `Interpolationtype`: Grid interpolation 


is used only for ``type grid`` and specifies
the method used for the numerical interpolation of the numerical values.
The currently implemented interpolation methods are `Cubicsplines` and `Quinticsplines` (default).

Example:
::

    Interpolationtype Cubicsplines
    Interpolationtype Quinticsplines


* `factor`: Scaling factor  

This optional keyword permits to rescale any object by
an arbitrary multiplication factor. At the moment the accepted values are any real number,
the imaginary unit :math:`i`, the square root of two, written as ``sqrt(2)``, or products
of these quantities. To write a product simply leave a space between the factors, but do not
use the `*` sign. All factor can have a :math:`\pm` sign.
The default value for `factor` is 1. This keyword is useful, for example,
to temporarily zero a certain object without removing it from the input file.

Examples:
::

   factor 1.5
   factor -sqrt(2)
   factor -2 sqrt(2) i

In the last example the factor is read in as :math:`-2 \sqrt{2} i`.
Note that imaginary factors make sense only in some cases for some coupling terms (in particular, spin-orbit) 
in the Cartesian-representation, see Section~\ref{s:representations}.


* `units`

This keyword selects the units of measure used for the 
the object in question. Supported units are: `angstroms`
(default) and `bohr` for the bond lengths; `cm-1` (default),
`hartree` (aliases are `au`, `a.u.`, and `Eh`), and `eV` (electronvolts)
for energies; `debye` (default) and `ea0` (i.e., atomic units) for dipoles; units can appear in any order. 

Example:
::

    units angstrom cm-1 (default for poten, spin-orbit, lambda-doubling etc)
    units bohr cm-1
    units debye  (default)
    units ae0 bohr



* `values`  

This keyword starts the subsection containing the numerical
values defining the object. 
For one of the `type`s corresponding to an analytical function (see :ref:`functions`),
the input between `values` and `end` contains the values of the parameters of the function.
The input consists in two columns separated by spaces containing (i) a string label
identifying the parameter and (ii) the value of the parameter (a real number).

In case of `fitting` (see :ref:`fitting`) a third column should
also be provided; the parameters which are permitted to vary during fitting
must have in the third column the string `fit` or, alternatively, the letter `f`
or the number 1. Any other string or number (for example, the string `nofit` or the number 0)
implies the parameter should be kept at its initial value.
In the case of fitting, the keyword `link`
can be also appear at the end of each the line; this keyword permits to
cross-reference values from different objects and is explained
below in this section.

In the case of objects of type `grid`, the third column can be also used to specify if the grid point needs to vary. 
The first columns contains the bond length :math:`r_i` and a second with the value of the object.
In the case of object of the `abinitio` (`reference`) type and specified as `grid`
a third column can be used to specify the fitting weights (see :ref:`fitting`).

* `<x|Lz|y>`, `<z|Lz|xy>` (aliases `<a|Lz|b>` and `<1|Lz|2>`)  

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
space or `*` sign.


* `link`  

This special keyword is used in fitting
to force a set of parameters 
(which may be relative to a different object) to have the same value.
For example, in a typical situation one may want to fit a set of PECs and to constrain their
dissociation (asymptotic) energy to the same value (because they are expected from theory to share the same
dissociation channel).


After the keyword `link` one should provide three numbers :math:`i_1`, :math:`i_2`, :math:`i_3` defining the parameter ID, where
:math:`i_1` identifies the object type (e.g. `poten`, `spin-orbit`, `spin-rot` etc.), 
:math:`i_2` is the object number within the type :math:`i_1` and :math:`i_3` is the parameter number as it appears after `values`. The ID numbers :math:`i_1, i_2, i_3` 
are specified in the fitting outputs in the form `[i,j,k]`. 

Example of the input:
::

    DE     0.50960000000000E+05   fit     link   1   1   3

Example of the corresponding output
::

    DE     0.50960000000000E+05   [ 1   1   3 ]



* `morphing` This keyword is used for fitting and switches on the morphing method. 

* `ZPE`: Zero-point-energy 

`ZPE` allows to explicitly input the zero-point energy (ZPE) of the molecule (in cm:sup:`-1`). This affects the value printed, as by default
Duo  prints energy of rovibronic levels by subtracting the ZPE. If not specified, the lowest energy of the first :math:`J`-block 
(independent of parity) will be used as appear on the line `Jlist`.

* `fit_factor` 

This factor (:math:`d_{\lambda}`) is used as a part of the reference *ab initio* curves of the `abinitio` type which (when given) 
is applied to the corresponding weights assigned to the corresponding values of this object. 
It is different from `fit_factor` defined within in :ref:`fitting`.

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



