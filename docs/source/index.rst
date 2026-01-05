Welcome to Duo's documentation!
###############################

.. toctree::
   :glob:
   :maxdepth: 4
   :caption: Contents:

   duo
   compile
   input
   energies
   spectra
   contraction
   fields
   functions
   fitting
   grid
   eigenfunc
   checkpoints
   unbound
   keywords
   hyperfine
   license


News
====

- 12.02.2023: The diabatic and adiabatic options now give the same results (after some bugs fixed in the NAC section)! 
- 07.02.2023: Reduced radial density is implemented as part of the checkpoint section. Can be activated using density save. It is cool! 
- 18.02.2023: The dipole intensity calculations have been optimised; the half_linestrength subroutine openmp parallelised (most intensive part) and the print out of the Einstein A coefficients (.trans) is now done in batches, not lin-by-line. 
- 16.02.2023: A problem in assign_by_count associated with the v=0 wrongly assigned at high J has been fixed.  
- 12.02.2023: The definition of the nuclear spin (used for g_ns in intensity calculations) is now fully automatic using the nuclear spin data from the internal database. Specifying "nspin" or "gns" are not required anymore.
- 06.02.2023: A new keyword states_only in INTENSITY introduced used to generate .states file and skip .trans. 
- 04.02.2023: New forms of PECs and DMCs are intriduced: poten_EHH and dipole_medvedev_sing
- 02.01.2023: State labels (strings) can be used as state reference numbers. e.g. "potential A" or "spin-orbit A C"
- 02.01.2023: Finding crossing point between two diabatic states and using it as an expansion centre. 
- 02.01.2023: relation between abinitio fields and their parent fields has been refactored to make it more stable when using in the fitting  
- 17.07.2022: bound/continuum state treatment is introduced: with "bound" and "unbound" in INTENSITY compute intensities for bound and continuum spectra. 

 
Indices and tables
==================

* :ref:`genindex`
* :ref:`search`
