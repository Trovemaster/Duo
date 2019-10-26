Specification of curves and couplings
=====================================

Once the main global parameters have been specified as described in the
previous sections, it is necessary to introduce the PECs and the various coupling
curves defining the Hamiltonian. Dipole moment curves (DMCs), which are necessary for
calculating spectral line intensities, are also discussed in this section, as well
as some special objects which are used for fitting.
Each object specification consists in a first part in which
keywords are given and a second one (starting from the
`values` keyword) in which numerical values are given;
the order of the keywords is not important, except for `values`.
Each object specification is terminated by the `end` keyword.

Objects of type `poten` (i.e., PECs, discussed in more detail below)
begin with a line of the kind `poten N`
where N is an integer index number counting over potentials and identifying them.
It is recommended that PECs are numbered progressively as 1,2,3,...,
although this only restriction is that the total number Nmax of PECs
should be not less than the total number of states specified by the keywork `nstates`.

Most other objects (e.g., `spin-orbit`) are assumed to be matrix elements
of some operator between electronic wave functions and after
the keyword identifying their type require two integer numbers
specifying the two indexes of the two electronic states involved (bra and ket).
The indexes are the numbers specified after the \texttt{poten} keyword.

Currently Duo supports the following types of objects:


* `poten` (alias: potential) 

Objects of type `poten` represent potential energy curves (PECs) and are
the most fundamental objects underlying each calculation.
From the point of view of theory each PEC is the solution of the electronic
Schoedinger equation with clamped nuclei, possibly complemented with the
scalar-relativistic correction and with the 
Born-Oppenheimer Diagonal correction
(also known as adiabatic correction). Approximate PECs can be obtained with
well-known quantum chemistry methods such as Hartree-Fock, coupled cluster theory etc.
Objects of type `poten` should always appear before
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



* L2:  (alias: `L**2`)

These objects represent matrix elements between electronic states of the molecule-fixed
  angular momentum operator :math:`\hat{L}^2 = \hat{L}_x^2 + \hat{L}_y^2 +\hat{L}_z^2`.


* L+:   (aliases: `Lplus`, `LxLy` and  `Lx`) 


It represent matrix elements between electronic states of the molecule-fixed
  angular momentum operator :math:`\hat{L}_+ = \hat{L}_x + i \hat{L}_y` and
  :math:`\hat{L}_x` in the :math:`\Lambda`- and Cartesian-representations, respectively.



* `spin-orbit` and `spin-orbit-x` 

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

For the `spin-orbit-x` case (:math:`\Lambda`-representation), the value of the matrix elements of the
 :math:`\hat{L}_z` operator nust be specified using the `<x|Lz|y>` keyword. 
 This representation is designed to work with e.g., the MOLPRO outputs. 
 For :math:`\Lambda\ne 0`, the diagonal SO-matrix element (e.g. between to :math:`\Pi`-components of :math:`\Lambda=1`) 
 should be specified using the :math:`\langle \Pi_x|LSZ |\Pi_y \rangle` component 
 (e.g. :math:`\langle 1.2 |{\rm LSZ} |1.3 \rangle`).

 %Matrix element of the spin-orbit Hamiltonian.


* `spin-spin-p` and `spin-spin-o` Parametrised phenomenological spin-spin 
operator (diagonal and off-diagonal. 

* `spin-rot` Matrix elements of the spin-rotational operator .

* `bob-rot`  (alias: `bobrot`) 
  
  Specifies the rotational :math:`g` factor (rotational Born-Oppenheimer breakdown term),
  which can be interpreted as a position-dependent modification to the rotational mass.

* `diabatic` (alias: `diabat`) 

Non-diagonal coupling of potential energy functions in the diabatic 
representation. 

* `lambda-opq`, `lambda-p2q`, and `lambda-q`  
  
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


* `abinitio` 
  
Objects of type `abinitio` (aliases: `reference`, `anchor`) are reference, `abinitio' curves which may be specified
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


* `dipole` (aliases: `dipole-moment`, `TM`) and `dipole-x`  
  
  Diagonal or transition dipole moment curves (DMCs),  necessary for computing 
  (dipole-allowed) transition line intensities and related quantities (Einstein :math:`A` coefficients etc.). 
  `dipole-x` is related to the Cartesian-representation.

  At the moment Duo cannot compute electric-quadrupole or magnetic dipole transition line intensities.



