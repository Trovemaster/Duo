Duo Input file: general structure
=================================



The input file is organized in self-contained input lines (e.g., 
::

     masses  1.00000  1.00000

specifies the masses of the two atoms in Daltons or in input sections beginning with
a specific keyword (e.g., `grid` and ending with the keyword `end`.

.. note:: The position of the keywords is not important. The input is not key-sensitive, 
          so `masses`, `MASSES`, `Masses` or any other combinations of uppercase and lowercase letters work 
          in exactly the same way.


A comma, a space or a hyphen (minus sign) can all be used as delimiters, so, e.g., one can also
:: 

     masses  1.00000, 1.00000

Sometimes keywords have several aliases, which are all equivalent.
Lines delimited by parentheses (i.e., round brackets) are ignored and can be used for comments.
If in the input there is a line with one of the keyword `END`, `STOP` or `FINISH` all lines after it are ignored.

Here is an example of a Duo inout to compute rovibrational energies of BeH in its ground electronic state using 
a grid-type potential energy curve by of Jacek Koput, JCP  135, 244308 (2011), Table III (see paper_).

.. _paper: http://dx.doi.org/10.1063/1.3671610


::
    
     atoms Be H
     (Total number of states taken into account)
     nstates 16
     
     (Total angular momentum quantum  - a value or an interval)
     jrot 0.5 - 2.5 
     
     (Defining the integration grid)
     grid
       npoints 501
       range   0.4 8.0
       type 0 
     end
     
     CONTRACTION
      vib
      vmax  30
     END

    
     poten 1
     units cm-1 angstroms
     name 'X2Sigma+'
     lambda 0
     symmetry +
     mult   2
     type grid
     values   
     0.60     105169.63
     0.65      77543.34
     0.70      55670.88
     0.75      38357.64
     0.80      24675.42
     0.85      13896.77
     0.90       5447.96
     0.95      -1125.87
     1.00      -6186.94
     1.05     -10024.96
     1.10     -12872.63
     1.15     -14917.62
     1.20     -16311.92
     1.25     -17179.13
     1.30     -17620.16
     1.32     -17696.29
     1.33     -17715.26
     1.34     -17722.22
     1.35     -17717.69
     1.36     -17702.19
     1.37     -17676.19
     1.38     -17640.16
     1.40     -17539.76
     1.45     -17142.53
     1.50     -16572.59
     1.55     -15868.72
     1.60     -15063.34
     1.65     -14183.71
     1.70     -13252.86
     1.80       -11313.
     1.90      -9369.74
     2.00      -7518.32
     2.10      -5832.29
     2.20      -4366.71
     2.30      -3155.94
     2.40      -2208.98
     2.50      -1507.72
     2.60      -1013.23
     2.80       -456.87
     3.00       -221.85
     3.50        -72.13
     4.00        -41.65
     4.50         -24.9
     5.00        -14.32
     6.00         -4.74
     8.00         -0.75
     10.00        -0.19
     20.00          0.0
    end
    

Control keys
------------

The following keys can appear anywere in the input file but outsides any sections. 


* ``Print_PECs_and_Couplings_to_File``

This keyword will tell Duo to print out all curves to a separate, auxiliary file. 

* ``Print_Vibrational_Energies_to_File``

This keyword is to print out all vibrational energies into a separate, auxiliary file. 

* ``Print_Rovibronic_Energies_To_File``

This keyword is to print out all rovibronic energies into a separate, auxiliary file.  

* DO_NOT_ECHO_INPUT is switch off the printing the inout file at the beginning of the output. 

* ``Do_not_Shift_PECs``

By default the PECs are shifted such that the minimum of the lowest PEC is at zero. This leads to Zero-Point-Energy (ZPE) to be 
defined relative to this zero.
All rovibronic energies are by default defined relative to the ZPE. This keyword will suppress shifting PECs so that ZPE is on the absolute scale. 

* ``DO_NOT_INCLUDE_JS_COUPLING``

This option is to switch the JS coupling in the Hamiltonian, can be used for debugging purposes. 

* ``ASSIGN_V_BY_COUNT``

This keyword will switch off the default assigning method  (based on the largest basis set contribution) of the vibrational to simple 
counting of the states, starting from :math:`v=0` within the same rotational-electronic configuration :math:`|{\rm State}, \Lambda, \Sigma, \Omega \rangle`. 
The default method  to assign the states with vibrational quantum numbers is known to fail at high excitations. 

* ``Legacy`` 

Aliases: ``Old-Version``, ``Version xxxx`` (xxxx is the year). This keyword to switch to the original, older version of the molpro function, 
which was modified in 2019 (bugs fixed and restructured). This keyword should help to reproduce the results published 
with the old version of the code. 
    

* ``L2Convention``

There are two conventions to include the electronic angular momentum :math:`\hat{L}_z^2` components: 
it can be defined either as part of the kinetic energy operator (``SPECIFY_L^2``, ``SPECIFY_L**2``,``Default``)
as :math:`\Lambda^2` or as part of the :math:`\hat{L}^2`  operator (``SPECIFY_LX^2_PLUS_LY^2``,``SPECIFY_LX**2_PLUS_LY**2``).

Example:
::

    L2Convention SPECIFY_LX^2_PLUS_LY^2


* ``Mem,``, ``Memory``: defines the maximal memory (RAM) available for the calculations.  

The program will stop with an error if the memory will be acceded before attempting to allocate a new array. The memory can be specified in 
``B``, ``Kb``, ``Mb``, ``Gb`` or ``Tb``. Example: 
::

    64 Gb 
    
    
