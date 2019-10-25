Duo keywords
============

.. note:: The position of the keywords is not important. The input is not key-sensitive.

* GNS: Nuclear statistical weight `gns 2.0`


* Temperature:  Temperature (K):

::

    Temperature  3000.0


`Aliases`: Temp, Temperature

+ Temperature can be split into rotational temperature and vibrational temperature, which in input is given as 
       

* `grid` is to define a multi-grid with different resolutions in different sub-grids as follows 

Example 
::     
    
    grid
      Range   0    100   Npoints 10000 offset 25.
      Range 100   1000  Npoints 1000 
      Range 1000 10000  Npoints 100
    end
     
