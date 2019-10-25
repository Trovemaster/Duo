Computing energy levels
=======================

In the following we present all keywords and options relevant to the calculations of energy levels.

Calculation (input) setup
^^^^^^^^^^^^^^^^^^^^^^^^^

* Atoms:  defines the chemical symbols of the two atoms.

Example:
::

    atoms Na-23 H-2

specifies the 23NaD diatomic. Duo includes an extensive database of atomic properties (atomic masses,
nuclear spins, isotopic abundances and other quantities) and will use the appropriate values
when required. The database should cover all naturally-occurring nuclei as well as all radioactive ones
with a half-life greater than one day and is based on the AME2012 and
NUBASE2012 databases. Each atom should be specified by its chemical symbol,
a hyphen (minus sign) and the atomic mass number, like in the example above.
Atomic masses will be used, which is generally the most appropriate choice unless one is explicitely including
non-adiabatic corrections. The hydrogen isotopes deuterium and tritium can also
be optionally specified by the symbols `D` and `T`.


The atomic mass number can be omitted, like in the following example:
::

     atoms Li F

In this case Duo will use the most-abundant isotopes (7Li and 19F in the example above) or, for
radioactive nuclei not naturally found, the longest lived one. For example 
::

     atoms Tc H

selects for technetium the isotope 97Tc, which is the longest-lived one. A few nuclides in the database are nuclear metastable isomers,
i.e. long-lived excited states of nuclei; these can be specified with a notation of the kind
::

     atoms Sb-120m H


In the example above the radioactive isotope of antimony 120mSb is specified (and hydrogen).
Another example
::

     atoms Sc-44m3 H

specifies the scandium radioactive isotope 44m3Sc (and hydrogen).



* masses:  

This is an optional keyword which specifies explicitely 
the masses of the two atoms (in Daltons, i.e. unified atomic mass units),
overriding the values from the internal database if the keyword \texttt{atoms} is also specified.

For example, the masses for the CaO molecule would be:
::

     masses 39.9625906 15.99491463
     
     
The masses may be atomic masses (the recommended choice if one does not include adiabatic or non-adiabatic corrections), nuclear masses.
An up-to-date reference of atomic masses is provided by the AME2012 catalogue (Chin. Phys. C, 36:1603–2014, 2012.)
Duo can also make use of position-dependent masses  (which is a practical way to account for non-adiabatic effects), see Bob-Rot.


* nstates:  is the number of potential energy curves (PECs) included in the calculation.
  
For example, if the ground state and four excited states of a molecule are to be included:
::

    nstates 5

Note that if `nstates` is set to a number different from the actual number of PECs  included in the 
input file no error message is issued; if more than nstates PECs are included in the input file then the PECs 
with `state` > `nstates`  will be ignored.

Note also that, consistently with the way Duo works internally, `nstates` is the number of unique PECs in absence of spin-orbit couplings.


* jrot: specifies the set of total angular momentum quantum numbers to be computed.

These must be integers or half-integers, depending on whether there is an even or odd number of electrons.
One can directly specify the values (separated by spaces or commas), specify a range of values (a minimum
and a maximum values separated by a hyphen; note than the hyphen must be surrounded
by at least by one space on each side). The values do not have to appear in ascending order.
For example, the following line
:: 
      
      jrot 2.5, 0.5, 10.5 - 12.5,  20.5
     

specifies the set `J = 0.5, 2.5, 10.5, 11.5, 12.5, 20.5`.

The first `J` in the `jrot` list will be used to define the reference zero-point-energy (ZPE) value for the
run.

Note that in the optional sections specifying calculation of spectra (See Intensity) or
specifying fitting (section \ref{sec:fitting}) is necessary to specify again a list of
J values   by `J` and `Jlist` respectively, which are completely independent from the `jrot` value
specified for energy level calculation.


* symmetry:  (or `Symgroup`) is an optional keywork specifying the molecular permutation-inversion symmetry group.

`Cs(M)` is for heteronuclear diatomics and `C2v(M)` is for homonuclear diatomics. 

For example: 
::

    symmetry Cs(M)

Instead of `C2v(M)` one can write equivalently `C2h(M)` or `G4(M)`, as these groups are isomorphic; 
the only difference will be in the labels used for the energy levels.
The short-hand notations `Cs`, `C2v`, `C2h` and `G4` can also be used and are equivalent to the ones with `(M)`.
The energy calculations are done using `Cs(M)`, which is also the default, while for the intensities the `C2v(M)` group can be also used.
Note that this keyword refers to the symmetry of the exact total (electronic, vibrational and rotational)
Hamiltonian and `not` to the :math:`C_{\infty v}` or the  :math:`D_{\infty h}` point groups, which are relative 
to the clamped-nuclei electronic Hamiltonian.

* DO_NOT_SHIFT_PECS: suppresses shifting the PEC to the minimum of potential 1 (assumed to the lowest). 

The default is to do the shift of the PECs to the minimum of  `poten 1`. In order to suppress shifting  energies to ZPE, use
::

   ZPE 0.0

see also the description of the keyword `ZPE`.



* SOLUTIONMETHOD * defines the DVR basis set and thus the DVR solution method for the vibrational problem. 

Possible methods include `5POINTDIFFERENCES`.
:: 
     
     SOLUTIONMETHOD  5POINTDIFFERENCES
     

for the 5 points stencil finite differences to derive the kinetic energy operator. A more efficient method is Sinc DVR (default), which is switched on with
::

    SOLUTIONMETHOD  SINC

Since Sinc is also currently the default method, this does not have to be specified.



Defining the grid
^^^^^^^^^^^^^^^^^

* grid: specifies an input section with the specifications of the grid of points.

It is used for the solution of the vibrational problem.

Example:
::

     grid
       npoints 501
       unit angstroms
       range 1.48 , 2.65
       type 0
     end


Keywords 


* npoints: is the number of grid points Np. 
 
Typical runs use 100 to 500 points. 
::

   npoints 301

* units: is optional and specifies the unit of measure of the grid specifications; possible values are `angstroms` (default)
                    or `bohrs`.

* range:   specifies the range of the grid in terms of rmin and rmax the lower and upper bond lengths. 

rmin should be strictly greater than zero and rmax strictly greater than rmin. As elsewhere in the program, the value
may be separated by a space or a comma. 


* type: Grid type

Is an integer number :math:`\ge 0` which specifies the type of grid.

Duo support not only uniformely spaced grids (default), which correspond to `type 0`,
but also various kind on non-uniformely spaced ones, which are particularly useful for near-dissociation, 
very weakly bound states.
Example: 
::                   

     type 0

In the case of uniformely-spaced grids the mesh points `rj`, j=0, Np-1 are given by
:math:`r_j = r_\mathrm{min} + \Delta j \qquad \mbox{where} \qquad \Delta = \frac{ r_\mathrm{max}-r_\mathrm{min}} {N_p - 1}`


Non-uniformely spaced grids are based on a change of variables from `r` to `z=f(r)`; it is then
the transformed variable `z` that is uniformely sampled. The transformed variables `z` are 
parametrised by two parameters, :math:`r_e` and :math:`\alpha`, which have to be specified
for the grid types > 0 (see below).

Transformed variable currently implemented are (Phys. Rev. A, 78:052510, 2008):\\

:math:`\quad z=\exp( -e^{-\alpha (r-r_e)} )`

:math:`\quad z=1 - \left(1+ e^{\alpha (r-r_e)} \right)^{-1}`

:math:`\quad z= \arctan( \alpha (r-r_e) )`

:math:`\quad z= (y-1)/(y+1)` with :math:`y=(r/r_e)^\alpha`

All the transformed grids have the property of decreasing the density of points for large `r`, so that one does not `waste` too many points in
regions where the potential is almost constant and the corresponding vibrational wave function slowly varying.

* re:  (alias `ref`)  Reference bond length used for `type` > 0 (see above). 

* alpha:  Parameter :math:`\alpha` for `type >0` (see above). 


