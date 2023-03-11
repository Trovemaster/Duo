.. _grids:

Defining the grid
=================


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


**Keywords**


* npoints: is the number of grid points Np. 
 
Typical runs use 100 to 500 points. 
::

   npoints 301

* units: is optional and specifies the unit of measure of the grid specifications; possible values are ``angstroms`` (default)
                    or ``bohrs``.

* range:   specifies the range of the grid in terms of rmin and rmax the lower and upper bond lengths. 

rmin should be strictly greater than zero and rmax strictly greater than rmin. As elsewhere in the program, the value
may be separated by a space or a comma. 


* type: Grid type

Is an integer number :math:`\ge 0` which specifies the type of grid.

Duo support not only uniformely spaced grids (default), which correspond to ``type 0``,
but also various kind on non-uniformely spaced ones, which are particularly useful for near-dissociation, 
very weakly bound states.
Example: 
::                   

     type 0

*(Default)* In the case of uniformely-spaced grids the mesh points `rj`, j=0, Np-1 are given by
:math:`r_j = r_\mathrm{min} + \Delta j \qquad \mbox{where} \qquad \Delta = \frac{ r_\mathrm{max}-r_\mathrm{min}} {N_p - 1}`


Non-uniformely spaced grids are based on a change of variables from `r` to `z=f(r)`; it is then
the transformed variable `z` that is uniformely sampled. The transformed variables `z` are 
parametrised by two parameters, :math:`r_e` and :math:`\alpha`, which have to be specified
for the grid types > 0 (see below).

Transformed variable currently implemented are (Phys. Rev. A, 78:052510, 2008):

:math:`\quad z=\exp( -e^{-\alpha (r-r_e)} )`

:math:`\quad z=1 - \left(1+ e^{\alpha (r-r_e)} \right)^{-1}`

:math:`\quad z= \arctan( \alpha (r-r_e) )`

:math:`\quad z= (y-1)/(y+1)` with :math:`y=(r/r_e)^\alpha`

All the transformed grids have the property of decreasing the density of points for large `r`, so that one does not `waste` too many points in
regions where the potential is almost constant and the corresponding vibrational wave function slowly varying.

* re:  (alias ``ref``)  Reference bond length used for ``type`` > 0 (see above). 

* alpha:  Parameter :math:`\alpha` for ``type >0`` (see above). 
