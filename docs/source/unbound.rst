.. _unboud states

Treating unbound states
=======================

Identifying unbound states
^^^^^^^^^^^^^^^^^^^^^^^^^^

Unbound states appear above the state dissociations. The Duo is developed to treat bound state problems
with an effective  boundary condition for the rovibronic eigenfunctions to vanish at the borders of the simulation grid.
However it is a typical problem when an electronic system contans unbound region between states  or above  their dissociations,
where unboud solutions are possible. Moreover, some of these unbound eigenfunctions are exactly zero at
:math:`r= r_{\rm max}` (other side is automatically zero duo to the steep repulsive wall).
Duo methodology can be applied also to compute the unbound spectra, here we show how to remove the spurious unbound states from the spectra (line lists)
calculations and thus to produce bound a line list. The unbound wavefunctions :math:`\psi_{\lambda}(r)` can be identified based on their asymptotic properties
via non-zero density in the small region of :math:`\delta` at the right border :math:`r= r_{\rm max}`:

:math:`\epsilon = \int_{r_{\rm max - \delta}}^{r_{\rm max}} |\psi_{\lambda}(r)|^2 dr > \epsilon_{\rm thr}`

:math:`\epsilon \sim 10^{-8}` is a small threshold value. The default threshold value is  :math:`\epsilon \sim 10^{-8}` is chosen as :math:`\sqrt{\epsilon(1.0d0)} \sim 1.5 \times 10^{-8}`.
The threshold value can be specified in the input file using the keyword `thresh_bound`.

The default value of :math:`\delta`  is 0.5 :math:`\AA`, which can be changed using the keyword ``THRESH_DELTA_R`` (also in :math:`\AA`).


.. sidebar::

   .. figure:: img/AlH_density.jpg
       :alt: AlH density

       Reduced densities of AlH together with an integration box used to disentangle (quasi-)bound and continuum states, :math:`\delta = 40\,\AA`. 




Excluding  unbound states
^^^^^^^^^^^^^^^^^^^^^^^^^^

Once the unbound states are identified they can be excluded from the intensity or line list calculations using the `bound` keyword in the INTENSITY section,
which tells Duo to compute bound-bound spectra only, e.g.:

Example:
::

  intensity
    absorption
    bound
    thresh_intes 1e-50
    thresh_bound  1e-6
    thresh_delta_R   1 (Angstrom)
    temperature 3000 (K)
    nspin 2.5  1.5
    linelist AlCl-37_61_J160
    J 0, 20
    freq-window  0.0  48000.0
    energy low 0.0 30000, upper 0.0,48000
  end


Her, the integrated density :math:`\epsilon` of the given state over the region of :math:`\delta= 1 \AA` will be compared to the threshold value of :math:`\epsilon_{\rm thr} = 1e-6`. The state in question is considered as unbound if :math:`\epsilon>\epsilon_{\rm thr}`. 

This selection is very dependent on the geometry of the box and the integrated region and thus will be different for different systems. The selection criterium   for the quasi-bound states is not very well defined and must be chosen on the case-by-case basis. For example, it is sometimes useful to plot the corresponding reduced densities, check the lifetimes or even compare to the experiment to see what should be considered as 'bound' states. 

An alternative way, less dependent on the specific geometry is to specify the threshold using the average density over the specific integration region:

:math:`\bar\epsilon =  \frac{\epsilon}{\delta} > \bar\epsilon_{\rm thr}`

defined using the ``THRESH_AVERAGE_DENSITY`` keyword, for example:

    INTENSITY
      absorption
      bound
      THRESH_BOUND  0.1
      THRESH_DELTA_R   1 (Angstrom)
      THRESH_AVERAGE_density  1e-4
      THRESH_DELTA_R  4
      THRESH_INTES  1e-40
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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes only the transitions to the unbound state are needed. In this case we exclude tansitions to the upper bound states with a keyword `unbound` placed anywhere in the
INTENSITY section.

Example:
::

  intensity
    absorption
    unbound
    thresh_intes 1e-50
    thresh_bound  1e-6
    temperature 3000 (K)
    nspin 2.5  1.5
    linelist AlCl-37_61_J160
    J 0, 20
    freq-window  0.0  48000.0
    energy low 0.0 30000, upper 0.0,48000
  end

