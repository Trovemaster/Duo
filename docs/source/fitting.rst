.. _fitting:

Duo Fitting
===========

Duo allows the user to modify (`refine`)
the potential energy curves and other coupling curves
by least-squares-fit to `experimental` energy term values or wavenumbers.

Fitting is, by far, the trickiest part of Duo, both on the part of the
program itself and on the part of the user. While the calculation of energy levels
and spectra from given PECs, couplings and dipole curves is relatively straighforward
(the most critical point being the consistency of the
phases specified for the various coupling terms), fitting is often more difficult
and may require a trial-and-error approach.
Fitting is also the part of Duo where most improvements are to be expected in
future new versions.

Example of a fitting section:
::

  FITTING
  JLIST 2.5,0.5, 1.5 - 11.5, 22.5 - 112.5
  itmax 30
  fit_factor  1e6
  output alo_01
  fit_type dgelss
  lock      5.0
  robust      0.001
  energies   (J parity NN  energy ) (e-state v ilambda isigma omega  weight)
   0.5   +   1     0.0000  1  0  0   0.5   0.5    100.000
   0.5   +   2   965.4519  1  1  0   0.5   0.5      7.071
   0.5   +   3  1916.8596  1  2  0   0.5   0.5      5.774
   0.5   +   4  2854.2366  1  3  0   0.5   0.5      5.000
   0.5   +   5  3777.5016  1  4  0   0.5   0.5      4.472
   0.5   +   6  4686.7136  1  5  0   0.5   0.5      4.082
   0.5   +   7  5346.1146  2  0  1  -0.5   0.5    100.000
  end



Keywords
^^^^^^^^



* ``FITTING``  

This keyword marks the beginning of the fitting input section. The whole section
can be deactivated by putting ``none`` or ``off`` next to the keyword ``FITTING``. This is useful to disable the fitting
without removing the input block from the input file.

* ``jlist`` (aliases are ``jrot`` and ``J``)

This keyword allows the user to specify the values of the :math:`J` quantum number to be used in the fit.
It superseedes the corresponding ``jrot`` keyword specified in the general setup
Individual values of :math:`J` can be separated by spaces or commas, while ranges are specified by two values separated by a hyphen
(hyphens should be surrounded by spaces). The first :math:`J` value is used to determine ZPE. For example
::

    JLIST  1.5, 5.5, 15.5 - 25.5, 112.5

selects the values 1.5, 5.5, all values from 15.5 to 25.5 and the value 112.5.


* ``itmax`` (alias ``itermax``)  An integer defining the maximum number of fitting iterations.

Setting ``itmax`` to zero implies that no fit will be performed (straight-through  calculation); however, the differences between
the computed energy levels (or frequences) and the reference (experimental) ones will be printed.
Example:
::

    itmax 15

* ``fit_factor``  This factor is used when reference curves of the
    ``abinitio`` type are included in the fit and used to define the importance of the energy/frequency data relative to the reference ``abinitio`` data. This factor is applied to all energy (frequencies) weight factors :math:`w_i^{\rm en}`.


When the factor is very large (e.g.  :math:`10^6`, like in the example above) the penalty for
differing for the reference curve is very small, so that only the `obs. - calc.` for energy levels
matter. Vice versa, if the factor is very small (e.g.  :math:`10^{-6}`) the fit is constrained so that the fitted
curves stay very close to the reference (``abinitio``) ones. When this number is extremely small (smaller than :math:`10^{-16}`)
the experimental data are completely ignored and the fit is performed to the \ai\ values only. Thus this feature also allows one to use the ``FITTING`` section for building analytical representations (see ``type``-s currently available) of different objects by fitting to the corresponding \ai\ or reference data provided in the ``abinitio``-sections of the input.

Example:
::
 
 
    fit_factor 1e2


* ``lock``   


``Lock`` denotes the threshold (cm\ :sup:`-1`) for which the quantum numbers are locked. 
The quantum numbers defining ``state``, :math:`v`, :math:`|\lambda|`, :math:`|\sigma|` and :math:`|\Omega|` 
will be used to identify and lock the energy value in place of the row number within the :math:`J`/parity block. 
When negative, the match is reconstructed based solely on the closest value within the lock-threshold given. 
If the match within the lock-region is not found, the row :math:`J`/parity number is used to match the theoretical 
and experimental energies. For example to match and lock to the calculated energy to the `experimental` one based 
on the quantum numbers within 20 cm\ :sup:`-1` use:
::
 
   lock 20.0

* ``thresh_obs-calc``  This keywords triggers switching off states from the fit if the obs.-calc. residuals become larger than the threshold specified. This feature is useful in case of multiple  swapping of the states during the fits and even the lock ``option`` does not help. The default value is zero (the feature is off).

* ``robust``   This keyword allows the user to switch on

Watson's robust fitting procedure: ``0`` is `off`, any other positive value
is `on` and defines the target accuracy of the fit as given by the weighted  root-mean-square
error. The ``robust``-value  is the targeted accuracy (obs.-calc.) of the fit. 
Example:
::
 
    robust 0.01

* ``target_rms`` 


This is to define the convergence threshold   (cm\ :sup:`-1`) for the total, not-weighted root-mean-squares (rms) fitting error. 
Example:
::

    target_rms 0.1 

* ``output``

This is the `filename`  for the files `name`.en, `name`.freq and `name`.pot, containing
detailed information on the fitting, including the fitting residuals for each iteration.
Example:
::

   output NaH_fit 


* ``energies``   

This keyword starts the section with the
energy levels to be fit to (e.g., obtained from an analysis of the experimental
line positions). Energy levels are written as in the following example:
::

  energies
     0.5   +    1     0.0000 1  0 0 0.5  0.5  1.00
     0.5   +    2   965.4519 1  1 0 0.5  0.5  0.90
     0.5   +    3  1916.8596 1  2 0 0.5  0.5  0.80
  end

where the meaning of the various quantities is as follows; col.1 is the total angular momentum quantum number :math:`J`;

col. 2  either the total parity :math:`\tau = \pm` or the :math:`e/f` parity;

col. 3  is a running number :math:`N` couting levels in ascending order of the energy within a :math:`(J,\tau)` symmetry block;

col. 4  is the energy term value :math:`\tilde{E}`, in cm\ :sup:`-1`;

col. 5  is the electronic state index `state`, as numbered in the ``poten`` sections;

col. 6  is the vibrational quantum number :math:`v`;

col. 7  is the projection of the electronic angular momentum :math:`\Lambda` for the state in question (an integer);

col. 8  is the projection of the total electronic spin :math:`\Sigma` (integer of half integer);

col. 9  is the projection of the total angular momentum :math:`\Omega` (integer of half integer);

col. 10 is the weight :math:`W` of the experimental energy in question (a real and positive number  usually given by :math:`\sigma^{-2}`, where :math:`\sigma` is the uncertainty of the energy level).


* ``frequency``  (aliases are ``frequencies`` and ``wavenumbers``)

This keyword works similarly to the ``energies``  keyword above but starts the section specifying the wavenumbers (i.e., line positions) to be fitted to.

Example:
::

    frequencies
      0.0  +   2 0.0 +  1   720.0000   2  0   1  -1.0   0.5    1  0   0   0.0   0.0  1.00
      2.0  +  17 3.0 -  2  5638.1376   4  0   0   1.0   1.0    2  0  -1  -1.0  -2.0  1.00
      4.0  +  17 5.0 -  2  5627.5270   4  0   0   1.0   1.0    2  0  -1  -1.0  -2.0  1.00
      4.0  +  18 7.0 -  2  5616.7976   4  0   0   0.0   0.0    2  0  -1  -1.0  -2.0  1.00
    end


The meaning of the quantities in each line are the following (see the keyword ``energies`` 
above for an explanation of the symbols. The prime/double prime symbol correspond to upper/lower level):  
:math:`J'`, :math:`\tau'`, :math:`N'`, :math:`J''`, :math:`\tau''`, :math:`N''`; frequency (cm\ :sup:`-1`); 
state\ :math:`'`, :math:`v'`, :math:`\Lambda'`, :math:`\Sigma'`, :math:`\Omega'`; state\ :math:`''`, :math:`v''`, 
:math:`\Lambda''`, :math:`\Sigma''`, :math:`\Omega''`; weight.

* ``off``, ``none``  are  used to switch off ``Fitting``, ``Intensity`` or ``Overlap``, when put next to these keywords.



Structure of the fitting output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

During fitting Duo will print for each iterations the fitting residuals using the following structure
(the first line with numbers 1 to 20 is not part of the output but serves as a legend): 
::

     1  2     3   4           5           6          7         8   9    10   11    12    13    14  15  16    17    18    19    20
     
     1  1    0.5  +      0.0000      0.0000     0.0000  0.60E-02   1     0    1  -0.5   0.5   0.5   1   0     1  -0.5   0.5   0.5
     2  2    0.5  +   1970.2743   1970.3983    -0.1240  0.59E-02   1     1    1  -0.5   0.5   0.5   1   1     1  -0.5   0.5   0.5
     3  3    0.5  +   3869.6639   3869.7934    -0.1295  0.30E-02   1     2    1  -0.5   0.5   0.5   1   2     1  -0.5   0.5   0.5
     4  4    0.5  +   5698.7392   5699.2951    -0.5559  0.20E-02   1     3    1  -0.5   0.5   0.5   1   3     1  -0.5   0.5   0.5
     5  1    0.5  -      0.1001      0.0000     0.1001  0.60E-02   1     0   -1   0.5  -0.5   0.5   1   0    -1   0.5  -0.5   0.5
     6  2    0.5  -   1970.4156   1970.3983     0.0173  0.59E-02   1     1   -1   0.5  -0.5   0.5   1   1    -1   0.5  -0.5   0.5


The meaning of the quantities in the various columns is as follows; 

col.1 is a simple line counter :math:`i` counting over all lines; 

col.2 is a counter :math:`N` counting lines within each :math:`J, \tau` symmetry block;

col. 3 is :math:`J`; col. 4 is the parity :math:`\tau`; 

col.5,6 are, respetively, the reference (`Observed`) and the calculated value of the line position; 

col.7 is the difference between observed and computed line positions;

col. 8 is the weight assigned to the transition in the fit; 

col. 9 to 14 are the quantum numbers of the lower state: `state`, :math:`v`, :math:`\Lambda`, :math:`\Sigma`, :math:`\Omega` and :math:`S`; 

col. 15 to 20 are the quantum numbers for the upper state (same definition as for columns 9 to 14).


* The auxiliary files .en, .freq, .pot


The files `name`.en contains all computed term values together with the theoretical quantum numbers, compared to the experimental
values, when available, along with the `experimental` quantum numbers as specified in the
``fitting`` section, for all iterations of the least-squares fit. Here ``name``
is the file name as speficied by the ``output`` keyword. The output is in the same format as in the
standard output file (see above) with the difference that it contains all calculated
values (subject of the ``nroots`` keyword, see Section \ref{s:diagonaliser}). An
asterisk ``*`` at the end of the line indicates that either the theoretical and
``experimental`` assignments don't agree or a residuals obs.-calc. is too large (large than
the ``lock`` parameter).

The frequency file `name`.freq with the keyword
``frequencies``. It has a similar structure as the standard output, with the
difference that for each transition from the ``frequency`` section the program will
estimate additional transition frequencies involving energies (both lower and
upper) which are within ``lock`` cm\ :sup:`-1` of the corresponding input values. This is done
to facilitate the search for possible miss-assignment, which is typical for transitions.
This is printed out for all iterations.

The file `{name`.pot (``potential``) contains the
residuals between the fitted and the reference
curve (if specified by an ``abinitio`` object).
The file is overwritten at each iteration. 
