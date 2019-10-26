Computing spectra
*****************

Absorption or emission spectra as well as line lists, partition functions and other
related quantities can be computed by adding in the input file an
``intensity`` section.
Here is an example of its general structure:
::


  intensity
    absorption
    thresh_intens  1e-15
    thresh_coeff   1e-15
    temperature   300.0
    qstat         10.0
    gns           1.0 1.0
    ZPE  931.418890
    selection (rules) 1 1
    J,  0.5, 1.5
    freq-window  -0.001,  25000.0
    energy low   -0.001, 6000.00, upper   -0.00, 30000.0
  end

If the keyword ``intensity`` is followed by ``none`` or ``off`` then the calculation of
intensities is disabled and the section is ignored. This is useful to temporarily
avoid the intensity calculation without removing or commenting out
the relative input section from the input file.
The meaning of the keywords is explained in the following.

Keywords
^^^^^^^^

* ``absorption``, ``emission``, ``partfunc``

These keywords define the type of the spectra
(absorption or emission) or whether Duo should only compute the partition function.
This keyword should appear immediately after ``intensity``. 

Example:
::

   absorption
   emission
   partfunc



* ``J`` (aliases  ``Jrot``, ``Jlist``) 

This keywords defines the range of rotational angular momentum quantum numbers for which line transitions should be computed. 
Note that this parameter is independent from ``jrot`` specified in the general setup.

Example:
::

   J  0,10

.. note::
   Using the ``J`` keyword the intensity production can be split into independent  
   ``J`` :math:`J_{\rm min},J_{\rm max}` ranges. In order to prevent overlaps, the range :math:`J_{\rm min},J_{\rm max}` 
   does not include transitions :math:`J_{\rm min} \leftrightarrow  J_{\rm min}`, except for :math:`J_{\rm min} = 0.5`, 
   where the transitions :math:`0.5 \leftrightarrow 0.5` are included. Transitions :math:`0 \leftrightarrow 0` are forbidden.


* ``Lande`` Compute Lande :math:`g` factors and write to the .states file.

* ``energy low`` and ``upper``   


These keywords to restrict the calculation to transitions
between levels satisfying the specified lower and upper energy thresholds (in cm\ :sup:`-1`): 
In the following we select transitions for which the lower state is between 0 and 6000 cm\ :sup:`-1` and the upper state is between 10000 and 30000 cm\ :sup:`-1`:
::

   energy low 0.0, 6000.00, upper 10000., 30000.0


Note that in this context level energies are measured by setting the energy of the lowest energy level to zero,
i.e. they do not include the zero-point energy, in contrast with
the threshold ``enermax`` specified in the general setup.


* ``freq-window`` specifies a frequency window for line positions (in cm\ :sup:`-1`). 

Example:
::

   freq-window 0.001, 25000.0


* ``gns`` nuclear statistical weight

``gns`` specifies the nuclear statistical weight, which for heteronuclear diatomics
is given by :math:`g_{ns} = (2 I_1+1)(2I_2+1)`, where :math:`I_1` and :math:`I_2` are the spins of the two nuclei.
In the case of homonuclear diatomics four numbers are expected, one for each symmetry species of the
`C_{2v}`(M) or :math:`C_{2h}(M)` symmetry groups.
Example:
::

   GNS 3.0


For the :math:`C_{2v}`(M) or :math:`C_{2h}`(M) symmetries associated with the homonuclear molecules the :math:`g_{\rm ns}` values must be specified for all of the four irreducible representation in the order :math:`A_1`, :math:`A_2`, :math:`B_1`, :math:`B_2` and :math:`A_g`, :math:`A_u`, :math:`B_g`, :math:`B_u`, respectively.
::

    GNS 1.0 1.0 0.0 0.0


* ``overlap`` allows for printing vibrational overlap integral, aka Franck-Condon factors. 

The default is not to print (``off``). One can also explicitly switch the overlaps off by  adding ``off`` next to ``overlap``:
::

    overlap off

The format is
::

    < i,   v| i',   v'> = value

where ``i`` and ``i'`` are the electronic state numbers, ``v`` and ``v'`` are the vibrational labels and ``value`` is the overlap:
`` \langle i,v | i',v' \rangle.
`` 
* ``vib-dipole`` prints  out vibrational transition moments :math:`\langle i,v | \mu(r) | i',v' \rangle`. By default these values are print out whenever the ``intensity`` is invoked. In order to switch this option off write ``off`` next to ``vib-dipole``:
::

    vib-dipole`` off

The format is
::

    < i,   v| <State | mu | State'> i',   v'> = value

where ``i`` and ``i'`` are the electronic state numbers, ``v`` and ``v'`` are the vibrational labels, ``State`` is the electronic state label and ``value`` is the transition dipole moment.

* ``Temperature`` specifies the temperature (in Kelvin) to be used for the calculation of line intensities.

It can be considered as a reference temperature since the Einstein coefficients as the main computational product and are temperature independent. The partition function associated with this {``Temperature`` should be also specified.
Example:
::

   temperature  298.0

* ``qstat`` (aliases: ``part-func`` and ``Q``). 

This keyword is
    to specify the value of the partition function :math:`Q` for the reference temperature defined by {``Temperature``.
    If not given, :math:`Q` is computed by Duo.

Example:
::

   qstat 10.0


* ``ZPE``

This keyword defines the zero point energy (cm\ :sup:`-1`) used for the calculation of line intensities, overriding
the value specified by the same keyword in the ``EigenSolver`` input section.
It is important to explicitly specify ``ZPE`` when the ground rovibronic state (whose energy defined the ZPE)
is not included in the calculation. Omitting
this keyword corresponds to using as ZPE the energy of the lowest-lying level used in the calculation. 

Example:
::
   
   ZPE 931.418890



* ``Thresh-intes`` specifies a minimum intensity threshold (in cm/molecule) for printing the transition into the
    output file as well as into the line list. 
    
Example:
::

    Thresh-intes  1e-35


* ``Thresh-Einstein`` 

specifies a threshold for the Einstein coefficient (in 1/s) for printing out the
transition into the output file as well as into the line list.

Example:
::

  Thresh-Einstein  1e-50

* ``linelist`` specifies a file name for writing a line list in the ExoMol format.

Example:
::

    linelist ScH

In the example above two files will be written, ``ScH.states``, containing a list of energy levels,
and ``ScH.trans``, containing the line transition data (line positions and Einstein :math:`A` coefficients).
 


