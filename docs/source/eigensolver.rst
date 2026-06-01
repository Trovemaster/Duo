.. _Eigensolver:

Eigensolver: Specifying the eigen-solution of the Hamiltonian 
=============================================================

Eigensolver
^^^^^^^^^^^

The input section ``EigenSolver`` (aliases: ``FinalStates``, ``diagonalizer``, ``FinalStates``) 
specifies 
various options relative to the `J>0` and/or the coupled problem; it also specifies
the LAPACK routine which should be used for matrix diagonalization (both for the solution of the
vibrational problem and for the solution of the coupled problem).
Example:
::

    Eigensolver
      enermax 25000.0
      nroots 500
      ZPE 1200.0
      SYEVR
    END


**Keywords**


* nroots  is the number of energy levels of the coupled problem to be computed (for any of the specified values of ``jrot``).

Example: 
::

    nroots 500

* enermax:  is an energy threshold  

``Enermax`` (aliases: ``uplimit``, ``enercut``) is to select the energy levels of the coupled problem to be computed (cm\ :math:`-1`). 

For example:
::

    enermax 15000.


If both ``nroots`` and ``enermax`` are  specified then only levels satisfying both criteria are selected.
Note that the present ``enermax`` threshold 
is distinct from the homonymous one in the ``vibrationalbasis`` input section, as the latter refers to the solution of the `J=0` 
uncoupled problem while the one being discussed
at present refers to the solution of the full (rotationally excited and/or coupled) problem.

* ZPE: Zero point energy 

``ZPE`` allows to explicitly input the zero-point energy (ZPE) of the molecule (in cm\ :math:`-1`). This affects the value printed, as
Duo always prints energy of rovibronic levels by subtracting the ZPE. Example:
::

     ZPE 931.418890
     
If ``ZPE`` is not included Duo will define the ZPE value as the lowest computed energy for the first value of `J` 
listed next to the ``jrot`` keyword (``jlist``), from the positive parity block. Currently it is not possible to 
take an automatic ZPE from the negative parity block (it is however possible in the intensity and fitting parts of the output).  
Thus ZPE does not necessarily have to be from the ground electronic state. 
This ZPE taken from the ``eigensolver``/``diagonalizer`` section changes the energies in the main, standard  Duo output.

The ZPE shift can be suppressed by setting the ZPE value to zero. This should be done either in the ``Diagonalizer``, ``Fitting`` or ``Intensity`` 
sections, depending on the current task:
::

    ZPE 0.0

* SYEVR or SYEV: LAPACK Eigesolvers 

This optional keywords permits to specify which routine from the LAPACK library should be used for
matrix diagonalization. At the moment only the two options quoted are implemented.
Example:
::

     SYEV

The SYEV routine (default) first reduces the matrix to diagonalize to tridiagonal form
using orthogonal similarity transformations, and then the QR algorithm is applied to the
tridiagonal matrix to compute the eigenvalues and the eigenvectors.
The SYEVR routine also reduces the matrix to diagonalize to tridiagonal form
using orthogonal similarity transformations but then, whenever possible,
computes the eigenspectrum using Multiple Relatively Robust Representations (MR).
SYEVR might give better performance, although exact timings are system- and case-dependent.

