.. _representations:

Using ab initio couplings in Duo: Representations of the electronic wave functions
==================================================================================


Quantum chemistry programs generally use real-valued electronic wave functions which transform according to the irreducible representations
of the C:sub:`2v` point group (for heteronuclear diatomics) or of D:math:`2h` (for homonuclear diatomics).
On the other hand Duo internally assumes the electronic wave functions are eigenfunctions of the :math:`\hat{L}_z`
operator, which implies they must be complex valued for :math:`|\Lambda| > 0`. Converting from one representation to the other is simple, as

:math:` |\Lambda\rangle =\frac{1}{\sqrt{2}}\left[\mp |1\rangle - i|2\rangle \right].`

where :math:`1\rangle` and :math:`2\rangle` are two Cartesian components of the electronic wave functions in a quantum chemistry program. 
Duo uses the matrix elements of the :math:`\hat{L}_z` to reconstruct the transformation between two representations: 


The keyword `<x|Lz|y>! and `<z|Lz|xy>` (aliases `<a|Lz|b>` and `<1|Lz|2>`) is required when specifying coupling curves between electronic states
in the `MOLPRO` representation (`spin-orbit-x`, `Lx` and 'dipole-x')  with :math:`|\Lambda| > 0` 
in order to resolve ambiguities in the definition of the   degenerate components of each electronic state.
This is also the value of the matrix element of the :math:`\hat{L}_z` operator computed for
the two component spherical harmonic, degenerate functions :math:`|x\rangle` and :math:`|y\rangle' for the :math:`\Pi` states or 
:math:`|z\rangle` and :math:`|xy\rangle` for the :math:`\Delta` states etc. 
The corresponding `<x|Lz|y>` values for both coupled states must be provided.

Examples:
::

     <x|Lz|y>   i  -i
::
     
     <z|Lz|xy> -2i  i

This keyword is required for the couplings of the following types: `spin-orbit-x`, `Lx` and `dipole-x`. 
The suffix `-x` indicates that Duo expects the `x`-component (non-zero) of the corresponding coupling. 
 
This keyword should appear anywhere in the object section, before the `values` keyword. 
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
It is the overall :math:`\pm' sign which \Duo\ needs and cannot be otherwise guessed.
As shown in the examples above, each factor should be written in the form :math:`\pm |\Lambda | i` without any
space or `*` sign.


