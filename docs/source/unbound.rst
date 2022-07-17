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

:math:`\int_{r_{\rm max - \delta}}^{r_{\rm max}} |\psi_{\lambda}(r)|^2 dr > \epsilon` 

:math:`\epsilon \sim 10^{-8}` is a small thrshold value. The threshold is uotomaticall chosen as :math:`\sqrt{\epsilon(1.0)\sim 1.5 \times 10^{-8}`. 
The threshold value can be specified in the input file using the keyword `thresh_bound`. 


Excluding  unbound states
^^^^^^^^^^^^^^^^^^^^^^^^^^

Once the unboud states are identified they can be excluded from the intensity or line list calculations using the `bound` keyword in the INTENSITY section, 
which tells Duo to compute boud-bound spectra only. 


Excluding  bound upper states
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes only the transitions to the unbound state are needed. In this case we exclude tansitions to the upper bound states with a keyword `unbound` placed anywehere in the 
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
    


