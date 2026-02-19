.. _unbound-states:

Treating unbound states
=======================

Duo is primarily designed for bound-state rovibronic problems, using an effective boundary condition that wavefunctions vanish at the borders of the radial grid (:math:`r_{\rm min}` and :math:`r_{\rm max}`). In practice, however, many diatomic models include electronic states with dissociation limits inside the energy range of interest, or with continua above dissociation. In such cases Duo may produce solutions that are formally “eigenstates” on the finite grid but correspond to continuum-like (unbound) behaviour. These states can contaminate spectra and line lists if not treated carefully.

This section describes two practical criteria implemented in Duo to identify (un)bound character and how to use them to include or exclude states in intensity/line-list calculations. It also explains how to separate bound–continuum and continuum–bound contributions in intensity calculations.

Identifying unbound states using a boundary density test
--------------------------------------------------------

A continuum-like (unbound) wavefunction typically carries non-negligible probability density close to the outer grid boundary. Duo can therefore flag a state :math:`\psi_\lambda(r)` as unbound by integrating the density in a small region of width :math:`\delta` at the right boundary:

.. math::

   \epsilon = \int_{r_{\rm max}-\delta}^{r_{\rm max}} \left|\psi_\lambda(r)\right|^2 \, dr.

If :math:`\epsilon` exceeds a threshold :math:`\epsilon_{\rm thr}`, the state is treated as unbound:

.. math::

   \epsilon > \epsilon_{\rm thr}.

The threshold :math:`\epsilon_{\rm thr}` is controlled by the input keyword ``thresh_bound``. The integration width :math:`\delta` (in Å) is controlled by ``thresh_delta_r``.

Typical defaults are:

* :math:`\epsilon_{\rm thr} \sim 10^{-8}` (``thresh_bound``)
* :math:`\delta = 0.5\,\AA` (``thresh_delta_r``)

.. sidebar::

   .. figure:: img/AlH_density.jpg
      :alt: AlH density

      Reduced densities of AlH together with an integration box used to disentangle (quasi-)bound and continuum states (example shown with a large :math:`\delta` for illustration).

.. note::
   This criterion depends on the choice of :math:`r_{\rm max}` and :math:`\delta`. For quasi-bound (resonance) states the separation between “bound” and “continuum-like” behaviour is not universal; the thresholds often need to be tuned for a given molecule and energy range.

Identifying unbound states using :math:`\langle r \rangle`
----------------------------------------------------------

A complementary indicator of continuum-like behaviour is the expectation value of the bond length:

.. math::

   \langle r \rangle = \int_{r_{\rm min}}^{r_{\rm max}} \psi_\lambda(r)\, r \,
   \psi_\lambda^{*}(r)\, dr.

Duo can flag a state as unbound if :math:`\langle r \rangle` exceeds a user-defined threshold :math:`r_{\rm thresh}`. This threshold is specified via ``thresh_bound_rmax``.


.. note::
       When both criteria are enabled, to be a bound state, it must satisfy both criteria, for  :math:`\langle r \rangle < r_{\rm thresh}` and  :math:`\epsilon < \epsilon_{\rm thr}`. That is, a state is unbound if any of these two conditions is satisfied. 

The ``.states`` file can include the bound/unbound label (``b``/``u``) and the value of :math:`\langle r \rangle` used for classification (last column), e.g.
::

   1     0.000000      4     0.5 + e X2Pi         0  1    -0.5     0.5 b   1.145131
   2  2652.362427      4     0.5 + e X2Pi         1  1    -0.5     0.5 b   1.188210
   3  5200.319353      4     0.5 + e X2Pi         2  1    -0.5     0.5 b   1.230171
   4  7642.107722      4     0.5 + e X2Pi         3  1    -0.5     0.5 b   1.275067
   5  9977.930743      4     0.5 + e X2Pi         4  1    -0.5     0.5 b   1.314642
   6 25367.675041      4     0.5 + e B2Sigma-     0  0     0.5     0.5 b   1.229628

Excluding unbound states (bound–bound spectra)
----------------------------------------------

Once unbound states are identified, they can be excluded from intensity or line-list calculations by adding the keyword ``bound`` to the ``intensity`` section. This instructs Duo to compute **bound–bound** transitions only.

Example (using the density criterion):
::

  intensity
    absorption
    bound
    thresh_intens   1e-50
    thresh_bound    1e-6
    thresh_delta_r  1.0
    temperature     3000.0
    linelist        AlCl-37_61_J160
    J               0, 20
    freq-window     0.0, 48000.0
    energy low      0.0, 30000.0, upper 0.0, 48000.0
  end

Here, the integrated density :math:`\epsilon` over the region :math:`[r_{\rm max}-\delta, r_{\rm max}]` with :math:`\delta=1\,\AA` is compared to :math:`\epsilon_{\rm thr}=10^{-6}`. The state is considered unbound if :math:`\epsilon>\epsilon_{\rm thr}`.

Average-density option
^^^^^^^^^^^^^^^^^^^^^^

Because :math:`\epsilon` scales with :math:`\delta`, it can be convenient to use the **average density** over the integration region:

.. math::

   \bar{\epsilon} = \frac{\epsilon}{\delta}.

This criterion is enabled using ``thresh_average_density``, and the state is treated as unbound if :math:`\bar{\epsilon} > \bar{\epsilon}_{\rm thr}`.

Example:
::

  intensity
    absorption
    bound
    thresh_bound            0.1
    thresh_delta_r          1.0
    thresh_average_density  1e-4
    thresh_intens           1e-40
    thresh_line             1e-40
    thresh_dipole           1e-7
    temperature             750.0
    linelist                AlH_446_A-X_L60.695_J10
    J                       0.0, 1.0
    freq-window             0.0, 30000.0
    energy low             -0.001, 30000.0, upper -0.0, 30000.0
  end

The default value of :math:`\bar{\epsilon}_{\rm thr}` is typically :math:`\sim 10^{-8}`.

Selecting transitions involving unbound states
----------------------------------------------

Duo can include transitions involving continuum-like states using the keyword ``unbound`` inside the ``intensity`` section. This is useful for modelling photodissociation continua, predissociation features, or any spectra where one side of a transition lies above dissociation.

By default, a single keyword ``unbound`` enables **both** of the following classes of transitions:

* **bound → unbound** (unbound upper states), sometimes referred to as *bound–continuum*;
* **unbound → bound** (unbound lower states), sometimes referred to as *continuum–bound*.

To process these two contributions independently, Duo accepts an optional second selector after ``unbound``:

* ``unbound upper``: include only transitions whose **upper** state is classified as unbound;
* ``unbound lower``: include only transitions whose **lower** state is classified as unbound.

If the selector is omitted, ``unbound`` is equivalent to requesting both ``upper`` and ``lower`` contributions.

Examples
^^^^^^^^

Bound → unbound (unbound *upper* states only):
::

  intensity
    absorption
    unbound upper
    ...
  end

Unbound → bound (unbound *lower* states only):
::

  intensity
    absorption
    unbound lower
    ...
  end

Both classes simultaneously:
::

  intensity
    absorption
    unbound
    ...
  end

Using :math:`\langle r \rangle` thresholds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :math:`\langle r \rangle` criterion is enabled with ``thresh_bound_rmax``:
::

  intensity
    absorption
    unbound upper
    thresh_intens      1e-50
    thresh_bound       1e-6
    thresh_bound_rmax  3.0
    temperature        3000.0
    linelist           AlCl-37_61_J160
    J                  0, 20
    freq-window        0.0, 48000.0
    energy low         0.0, 30000.0, upper 0.0, 48000.0
  end

State-dependent thresholds can be provided by listing one value per electronic state, in the same order as the states are defined (e.g. by the ``states`` keyword):
::

  intensity
    absorption
    unbound upper
    ...
    thresh_bound_rmax  3.0  5.5  4.0
    ...
  end

.. note::
   When both criteria are requested, the :math:`\langle r \rangle` threshold can be used as a *boundness safeguard*: a state is considered **bound** if :math:`\langle r \rangle < r_{\rm thresh}` even if the boundary density test is marginal. A state is considered **unbound** if it is not **bound**.

Printing unbound diagnostics in the ``.states`` file
----------------------------------------------------

By default, unbound handling adds a bound/unbound label (``b``/``u``) to the ``.states`` file. When ``thresh_bound_rmax`` is used, the corresponding value of :math:`\langle r \rangle` is printed next to this label. The boundary-density integral :math:`\epsilon` can also be printed by adding ``print_bound_density``:

::

  intensity
    absorption
    unbound
    thresh_bound       1e-4
    thresh_delta_r     1.0
    thresh_bound_rmax  2.0
    print_bound_density
    temperature        3000.0
    linelist           CH
    J                  0.5, 10.5
    freq-window        0.0, 50000.0
    energy low        -0.001, 20000.0, upper -0.0, 50000.0
  end

Example excerpt from the resulting ``.states`` file:
::

   1     0.000000      4     0.5 + e X2Pi         0  1    -0.5     0.5 b   1.13934 0.58E-21
   2  2720.584833      4     0.5 + e X2Pi         1  1    -0.5     0.5 b   1.17893 0.17E-21
   3  5314.682367      4     0.5 + e X2Pi         2  1    -0.5     0.5 b   1.22026 0.47E-21
   4  7784.710908      4     0.5 + e X2Pi         3  1    -0.5     0.5 b   1.26339 0.97E-22
   5 10132.366441      4     0.5 + e X2Pi         4  1    -0.5     0.5 b   1.30890 0.31E-21
   6 25771.515741      4     0.5 + e B2Sigma-     0  0     0.5     0.5 b   1.21463 0.83E-21
   7 27530.492176      4     0.5 + e B2Sigma-     1  0     0.5     0.5 b   1.32459 0.28E-19
   8 27877.572537      4     0.5 + e B2Sigma-     2  0     0.5     0.5 u   4.15693 0.17E-01
   9 27905.414170      4     0.5 + e B2Sigma-     3  0     0.5     0.5 u   4.96461 0.18
  10 27947.636255      4     0.5 + e B2Sigma-     4  0     0.5     0.5 u   4.79088 0.28
   ...

Keywords
--------

.. glossary::
   :sorted:

   bound
      Used inside the ``intensity`` block to compute **bound–bound** spectra/line lists only.
      Duo first classifies rovibronic states as bound/unbound using the criteria controlled by
      :term:`thresh_bound`, :term:`thresh_delta_r`, :term:`thresh_average_density`, and/or
      :term:`thresh_bound_rmax`, and then keeps only transitions whose **lower and upper**
      levels are classified as bound.

   unbound
      Used inside the ``intensity`` block to include transitions involving unbound states.

      * ``unbound`` (without a selector) enables both classes simultaneously: bound→unbound and unbound→bound.
      * ``unbound upper`` keeps only transitions whose **upper** level is classified as unbound (bound→unbound).
      * ``unbound lower`` keeps only transitions whose **lower** level is classified as unbound (unbound→bound).

      Examples:
      ::

         intensity
           absorption
           unbound upper
           ...
         end

      ::

         intensity
           absorption
           unbound lower
           ...
         end

   upper
      Selector used only as a second keyword after :term:`unbound` in the ``intensity`` block.
      ``unbound upper`` keeps transitions with **unbound upper** levels (bound→unbound).

   lower
      Selector used only as a second keyword after :term:`unbound` in the ``intensity`` block.
      ``unbound lower`` keeps transitions with **unbound lower** levels (unbound→bound).

   thresh_bound
      Threshold :math:`\epsilon_{\rm thr}` used in the boundary-density criterion. Duo evaluates

      .. math::

         \epsilon = \int_{r_{\rm max}-\delta}^{r_{\rm max}} |\psi_\lambda(r)|^2 \, dr,

      where :math:`\delta` is controlled by :term:`thresh_delta_r`. A state is treated as unbound
      if :math:`\epsilon > \epsilon_{\rm thr}`.

   thresh_delta_r
      Width :math:`\delta` (in Å) of the outer-boundary integration region used in the
      boundary-density criterion for :term:`thresh_bound` (and, if enabled, for the average-density
      criterion :term:`thresh_average_density`).

   thresh_average_density
      Enables an average-density variant of the boundary criterion using

      .. math::

         \bar{\epsilon} = \epsilon/\delta.

      A state is treated as unbound if :math:`\bar{\epsilon} > \bar{\epsilon}_{\rm thr}`.

   thresh_bound_rmax
      Threshold :math:`r_{\rm thresh}` (in Å) used in the :math:`\langle r \rangle` criterion.
      Duo computes

      .. math::

         \langle r \rangle = \int_{r_{\rm min}}^{r_{\rm max}} \psi_\lambda(r)\, r\,
         \psi_\lambda^{*}(r)\, dr,

      and flags a state as unbound if :math:`\langle r \rangle > r_{\rm thresh}`.

      The keyword can accept either a single value (applied to all states), or one value per
      electronic state (applied in the same order as the electronic states are defined in the input).

   print_bound_density
      Requests printing of the boundary-density integral :math:`\epsilon` in the ``.states`` file
      (typically as an extra final column), alongside the bound/unbound label and (if requested)
      :math:`\langle r \rangle`.
