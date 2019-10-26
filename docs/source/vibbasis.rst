Vibrational basis set
=====================

The keyword ``vibrationalbasis`` (aliases: ``vibrations``, ``contraction``) 
specifies the size of the vibrational basis set.

Duo uses a ``contraction`` scheme to construct the rovibronic basis set used for the solution
of the coupled problem. As a first step the `J=0` vibration problem is solved for each electronic state, in which the
corresponding Schroedinger equation is solved in the grid representation
of ``npoints``. Then a certain number  of the resulted
vibrational eigenfunctions :math:`|v\rangle` with :math:`0 \le v\le` vmax and :math:`\tilde{E} \le` ``EnerMax``  is selected to
form the vibrational part of the basis set

:math:`| J \Omega S \Sigma \Lambda v \rangle = | J \Omega \rangle | S \Sigma \rangle | \Lambda \rangle | v \rangle`

where :math:`| J \Omega \rangle`  and :math:`| S \Sigma \rangle`  are the rigid rotor functions and :math:`| \Lambda \rangle`  are the
electronic wavefunctions implicitly taken from the ab initio calculations.
Example: 
:: 


     vibrationalbasis
       vmax 30
       enermax 25000
     end


Keywords
^^^^^^^^


* vmax

(alias: ``vibmax``) specifies the value of the maximum vibrational quantum number to be computed and kept for
the solution of the coupled problem. For example
::

    vmax 15

specifies to compute for each PEC the lowest-energy 15 vibrational levels; it is also possible 
to specify different values of \texttt{vmax} for each PEC, in which case the values must be given as a list; for example
::

    vmax 10 15 8


specifies that for the PEC identified as ``poten 1`` Duo should use ``vmax=10``, for
``poten 2``, ``vmax=15`` and for ``poten 3``, ``vmax=8``.
If there are more PEC (``poten 4`` etc.) they will use for ``vmax`` the last value specified (``vmax=8`` in this example).

* enermax

Alternatively or complementary to ``vmax`` one can select the vibrational energy levels to compute
by specifying an upper energy threshold (in cm\ :sup:`-1`). Similarly to ``vmax``, one can specify a different value of ``enermax``
for each PEC by writing a list of values; for example
::

      enermax 30000.0 25000.0
      
      
selects a threshold of 30000 cm\ :sup:`-1`  for ``poten 1`` and one of 25000 cm\ :sup:`-1` for ``poten 2`` and any other potential present.
Note that by default Duo will shift the PECs so that the lowest point of the lowest-lying PEC has zero energy, and that the energy
used for the ``enermax`` threshold are ``total`` vibrational energies including the zero point energy.
One can prevent Duo from shifting the PECs by writing in the input (anywhere but not within an input section)
the option ``do_not_shift_pecs``.

If both ``enermax`` and ``vmax`` are specified only levels which satisfy both criteria are kept for the solution of the coupled problem.
If neither of them is specified (or the ``vibrationalbasis`` input section is missing altogether) then ``vmax``
is taken equal to ``npoints`` for all PECs and there is a hard-coded limit of 10\   :sup:`8` cm\ :sup:`-1` for ``enermax``.

\end{itemize}
