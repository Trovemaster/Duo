.. _unboud states

Treating unbound states
=======================

Identifying unbound states using density
----------------------------------------

Unbound states appear above the state dissociations. The Duo is developed to treat bound state problems
with an effective  boundary condition for the rovibronic eigenfunctions to vanish at the borders of the simulation grid.
However it is a typical problem when an electronic system contans unbound region between states  or above  their dissociations,
where unboud solutions are possible. Moreover, some of these unbound eigenfunctions are exactly zero at
:math:`r= r_{\rm max}` (other side is automatically zero duo to the steep repulsive wall).
Duo methodology can be applied also to compute the unbound spectra, here we show how to remove the spurious unbound states from the spectra (line lists)
calculations and thus to produce bound a line list. The unbound wavefunctions :math:`\psi_{\lambda}(r)` can be identified based on their asymptotic properties
via non-zero density in the small region of :math:`\delta` at the right border :math:`r= r_{\rm max}`:

.. math::
       
       \epsilon = \int_{r_{\rm max - \delta}}^{r_{\rm max}} |\psi_{\lambda}(r)|^2 dr > \epsilon_{\rm thr}
       
where :math:`\epsilon \sim 10^{-8}` is a small threshold value. The default threshold value is  :math:`\epsilon \sim 10^{-8}` is chosen as :math:`\sqrt{\epsilon(1.0d0)} \sim 1.5 \times 10^{-8}`.
The threshold value can be specified in the input file using the keyword `thresh_bound`.

The default value of :math:`\delta`  is 0.5 :math:`\AA`, which can be changed using the keyword ``THRESH_DELTA_R`` (also in :math:`\AA`).


.. sidebar::

   .. figure:: img/AlH_density.jpg
       :alt: AlH density

       Reduced densities of AlH together with an integration box used to disentangle (quasi-)bound and continuum states, :math:`\delta = 40\,\AA`.

Identifying unbound states using expectation value of the bond length
---------------------------------------------------------------------

Here we use the expectation values of :math:`r` as the measure of the unbound/bound character of the wavefunction. To this end, we calculate:

.. math::
       
       \bar{r} = \int_{r_{\rm min}}^{r_{\rm max}} \psi_{\lambda}(r) r \psi_{\lambda}^{*}(r)  dr
        
and compare it to a threshold value :math:`r_{\rm max}` ``thresh_bound_rmax``. The state is labeled "unbound" if :math:`r>r_{\rm thresh}`. 

This mechanism can be used in conjunction with the density criteria. First, Duo will check against the density criterium and then apply the expectation value criterium. 

The .states file in the case of the ``unbound`` calculations also provides the value the value of :math:`\bar{r}` used to decide on the unbound/bound property (last column), as well as a label ``u``/``b`` defining if it is unbound/bound (second from last), for example 

::
     
           1     0.000000      4     0.5 + e X2Pi         0  1    -0.5     0.5 b   1.145131
           2  2652.362427      4     0.5 + e X2Pi         1  1    -0.5     0.5 b   1.188210
           3  5200.319353      4     0.5 + e X2Pi         2  1    -0.5     0.5 b   1.230171
           4  7642.107722      4     0.5 + e X2Pi         3  1    -0.5     0.5 b   1.275067
           5  9977.930743      4     0.5 + e X2Pi         4  1    -0.5     0.5 b   1.314642
           6 25367.675041      4     0.5 + e B2sigma-     0  0     0.5     0.5 b   1.229628
     


Excluding  unbound states
-------------------------

Once the unbound states are identified they can be excluded from the intensity or line list calculations using the `bound` keyword in the INTENSITY section,
which tells Duo to compute bound-bound spectra only, e.g.:

Example (using the density criterium):
::

  intensity
    absorption
    bound
    thresh_intens 1e-50
    thresh_bound  1e-6
    thresh_delta_R   1 (Angstrom)
    temperature 3000 (K)
    linelist AlCl-37_61_J160
    J 0, 20
    freq-window  0.0  48000.0
    energy low 0.0 30000, upper 0.0,48000
  end


Here, the integrated density :math:`\epsilon` of the given state over the region of :math:`\delta= 1 {\rm \AA}` will be compared to the threshold value of :math:`\epsilon_{\rm thr} = 10^{-6}`. The state in question is considered as unbound if :math:`\epsilon>\epsilon_{\rm thr}`.

This selection is very dependent on the geometry of the box and the integrated region and thus will be different for different systems. The selection criterium   for the quasi-bound states is not very well defined and must be chosen on the case-by-case basis. For example, it is sometimes useful to plot the corresponding reduced densities, check the lifetimes or even compare to the experiment to see what should be considered as 'bound' states.

An alternative way, less dependent on the specific geometry is to specify the threshold using the average density over the specific integration region:

:math:`\bar\epsilon =  \frac{\epsilon}{\delta} > \bar\epsilon_{\rm thr}`

defined using the ``THRESH_AVERAGE_DENSITY`` keyword, for example:
::
     
    INTENSITY
      absorption
      bound
      THRESH_BOUND  0.1
      THRESH_DELTA_R   1 (Angstrom)
      THRESH_AVERAGE_density  1e-4
      THRESH_DELTA_R  4
      THRESH_INTENS  1e-40
      THRESH_LINE   1e-40
      thresh_dipole 1e-7
      TEMPERATURE   750
      linelist  AlH_446_A-X_L60.695_J10
      J,  0.0, 1
      freq-window    0.0,   30000.0
      energy low   -0.001, 30000.00, upper   -0.00, 30000.0
    END
    
    
The default value of :math:`\bar\epsilon_{\rm thr}` is :math:`\sim 10^{-8}`.



Excluding  bound upper states
-----------------------------

Sometimes only the transitions to the unbound state are needed. In this case we exclude tansitions to the upper bound states with a keyword `unbound` placed anywhere in the
INTENSITY section.

Example:
::

  intensity
    absorption
    unbound
    thresh_intens 1e-50
    thresh_bound  1e-6
    temperature 3000 (K)
    linelist AlCl-37_61_J160
    J 0, 20
    freq-window  0.0  48000.0
    energy low 0.0 30000, upper 0.0,48000
  end



Here is an example excluding  bound upper states using the criterium for the bond-expectation value:
::

  intensity
    absorption
    bound
    thresh_intens 1e-50
    thresh_bound  1e-6
    thresh_bound_rmax  2
    temperature 3000 (K)
    linelist AlCl-37_61_J160
    J 0, 20
    freq-window  0.0  48000.0
    energy low 0.0 30000, upper 0.0,48000
  end


where ``thresh_bound_rmax`` defines the value of :math:`r_{\rm max}` in the equation above. 

