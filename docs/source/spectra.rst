.. _computing-spectra:

Computing spectra (intensities and line lists)
**********************************************

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


It is possible to run intensities by batches, J1-J2: 0-1,1-2,2-3.... In order to prevern any overlaps, the 
transitions between J1-J1 are always omitted from the intensity calculations. 
The line list will contain states and transitions for the given batch only and the states with J1 will be 
exluded from the .states file, except for J1=0 or J1=0.5. The complete line list is the produced by stitchig  
together the individual  .statets and .trans files. The states numbering will be internally consistent. 

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


* ``quadrupole``

This keyword is used to trigger the calculation of electric quadrupole transitions (see also :ref:`quadrupole curves`).

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


* ``Richmol`` (`matelem`) is to generate the **RichMol** checkpoints. 

The **Richmol** checkpoint files ``.rchk`` are to be used for laser-driven molecular dynamics, see [18OwYaxx].

:cite:`18OwYa`.

* ``Raman`` or ``Polarizability`` is to generate the matrix elements for the **Raman** intensity calculations with Richmol. 

* ``Overlap`` is to trigger on/off the vibrational overlap integrals. 

The default is on, for switching of: 
::

   overlap off 

* ``VIB-DIPOLE`` is to generate transition dipole moments, and similarly ``VIB-QUADRUPOLE`` to generate transition quadrupole moments. 

Example
::

   VIB-DIPOLE 
   VIB-DIPOLE off 


* ``freq-window`` specifies a frequency window for line positions (in cm\ :sup:`-1`). 

Example:
::

   freq-window 0.001, 25000.0


* ``gns``  nuclear statistical weight (depriciated). **It is now automatically estimated using nuclear spin values
from the internal atomic and nuclear database.** 

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
 
 
* ``Nspin``  Nuclear spins of both atoms (**depriciated**). 

**The nuclear spins are now provided in the internal atomic and nuclear databases and is not required 
to be specified anymore.** 

The nuclear spin values are used to define the nuclear degeneracy factors as follows. Example
::

    nspin 0.0 0.5

::      
    nspin 0.0 0.0  

The nuclear degeneracy factors :math:`g_ns` are defined as follows. For the heteronuclear molecules:

:math:`g_{ns} = (2 I_1+1)(2I_2+1)`

For a homonuclear diatomic, it is given by 

:math:`g_{ns}^{A} = \frac{1}{2} ((2 I+1)^2+(2 I +1))`

and 

:math:`g_{ns}^{B} = \frac{1}{2} ((2 I+1)^2-(2 I +1))`

where :math:`I_1, I_2`  and `I` are the nuclear spins and `A` and `B` are the two irreps of the D2h symmetry group. 

 
* ``Gns`` is an alternative to ``nspin`` defining the nuclear spin degeneracy explicitly. 
 
Example: 
::

       GNS 3.0 3.0

::

       GNS 1.0 1.0 0.0 0.0 
 
 
Thresholds 
^^^^^^^^^^


** ``THRESH_LINE`` line strength  threshold (Debye:sup:`2`)

** ``THRESH_EINSTEIN`` Einstein A coefficient threshold (1/s).
 
** ``thresh_intes`` intensity (TM) threshold (cm/molecule)

** ``THRESH_DIPOLE`` transition dipole threshold (debye)


* ``states-only``, ``states_only``: to switch off the transition intensity when building the line list. When this option is given in 
the INTENSITY block, only a .states file is generated.  
    
Example:
::

    INTENSITY
    ....
    linelist CH_A-X
    states-only
    .....
    END


** ``bound`` used in INTENSITY calculations to produce bound-only spectra or linelists. For ``bound``, Duo 
identifies bound wavefunctions corresponding to the upper state and uses them to compute bound transition intensities. See also :ref:`_unboud states`. 

** ``unbound`` (oposite of ``bound``) used in INTENSITY calculations to produce unbound-only spectra or linelists. For ``unbound``, Duo 
identifies unbound wavefunctions corresponding to the upper state and uses them to compute unbound transition intensities. See also :ref:`_unboud states`. 



Example: Intensities of BeH
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Here we use the potential energy function of BeH from the example :ref:`energy_BeH`. 

For intensity calculations one needs an electric dipole moment curve, which we take from the spectroscopic model used in 
the ExoMol-I_ paper by Yadin et. al (2011)

.. _ExoMol-I: http://exomol.com/db/BeH/9Be-1H/Yadin/9Be-1H__Yadin__LEVEL8.0.inp

::
    
    dipole  1 1
    name "<2Sigma+|DMZ|2Sigma+>"
    spin   0.5 0.5
    lambda  0  0
    type   grid
    values 
       0.400     -0.4166624920
       0.500     -0.0241871531
       0.600      0.2217732500
       0.700      0.3386323420
       0.800      0.3661076190
       0.900      0.3311512400
       1.000      0.2513061130
       1.100      0.1379591390
       1.200     -0.0012406430
       1.300     -0.1588361650
       1.320     -0.1920270000
       1.340     -0.2256736540
       1.350     -0.2426539090
       1.360     -0.2597311920
       1.400     -0.3288944440
       1.500     -0.5056369720
       1.600     -0.6824442480
       1.700     -0.8513506410
       1.800     -1.0025214800
       1.900     -1.1238133700
       1.950     -1.1687609400
       2.000     -1.2005094800
       2.020     -1.2089972000
       2.050     -1.2166847200
       2.070     -1.2181089800
       2.100     -1.2136337000
       2.300     -1.0182994100
       2.400     -0.8538885220
       2.500     -0.6736179730
       2.600     -0.5046631750
       2.700     -0.3634556350
       2.800     -0.2548814520
       2.900     -0.1758884440
       3.000     -0.1201861300
       3.100     -0.0815224742
       3.200     -0.0549121655
       3.300     -0.0367099205
       3.400     -0.0243335573
       3.500     -0.0159701097
       3.600     -0.0103484461
       3.700     -0.0065800412
       3.800     -0.0040495078
       3.900     -0.0023383813
       4.000     -0.0011684378
       4.200      0.0002034367
       4.400      0.0008546009
       4.600      0.0011177434
       4.800      0.0011645509
       5.000      0.0011023829
       6.000      0.0005429083
       8.000     -0.0000033249
      10.000     -0.0000085504
    end
        
    INTENSITY
     absorption
     thresh_intes  1e-30
     thresh_line   1e-30
     temperature   300.0
     nspin         1.5  0.5 (see Wikipedia isotope Be)
     selection (rules) 1 1
     linelist   BeH
     J,  0.5, 10.5
     freq-window   0.0,  7000.0
     energy low   -0.001, 5000.00, upper   -0.00, 12000.0
    END
    

This will produce a line list for BeH in ExoMol format in two files .states and .trans, 
which can be processed using ExoCross_, see also ExoCross-tutorial_. 




.. _ExoCross: https://github.com/Trovemaster/exocross

.. _ExoCross-tutorial: https://github.com/Trovemaster/exocross/wiki/Configuring-the-ExoCross-session


.. bibliography:: references.bib



