.. _functions:

Duo Functions
=============

This section shows examples of the definitions of the analytical functions supported in Duo.


Extended Morse Oscillator ``EMO`` 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:math:`V(r)=T_{\rm e} + (A_{\rm e}-T_{\rm e})\left( 1 - \exp\left\{-\beta_{\rm EMO}(r) (r-r_{\rm e})\right\} \right)^2`,

which has the form of a Morse potential with a exponential tail and the distance-dependent exponent coefficient

:math:`\beta_{\rm EMO}(r) =  \sum_{i=0} a_i y_p^{\rm eq}(r)^i`,

expressed as a simple power series in the reduced variable:

:math:`y_p^{\rm e}(r) = \frac{r^p-r_{\rm e}^p}{r^p+r_{\rm e}^p}`

with :math:`p` as a parameter. This form guarantees the correct dissociation limit and allows
for extra flexibility in the degree of the polynomial on the left or on the right sides
of a reference position :math:`R_{\rm ref}` which we take at :math:`R_{\rm ref} = r_{\rm e}`. This is
specified by the parameters :math:`N=` :math:`N_{l}` (:math:`N_{r}`) and  :math:`p=` :math:`p_{l}` (:math:`p_{r}`),
respectively.

Example:
::

    poten 2
    name "a 3Piu"
    symmetry u
    type  EMO
    lambda 1
    mult   3
    values
    Te          0.81769829519421E+03
    Re          0.13115676812526E+01
    Ae          0.50960000000000E+05
    RREF       -0.10000000000000E+01
    PL          4
    PR          4
    NL          2
    NR          3
    a0          0.21868146887665E+01
    a1          0.88875855351916E-01
    a2          0.84932592800179E-01
    a3          0.23343175838290E+00
    end



**Taylor expansion around** :math:`r_0`:

:math:`V(r) = T_{\rm e} + (A_{\rm e} - T_{\rm e}) a_0^2 (r-r_0)^2 + (A_{\rm e} - T_{\rm e}) \left( p \frac{a_0 a_1}{r_e} - a_0^3 \right) (r-r_0)^3 + \cdots`


Morse Long-Range (MLR) function ``MLR``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


:math:`V(r) = T_{\rm e}+ (A_{\rm e}-T_{\rm e}) \left(1 - \frac{u_{\textrm{LR}}(r)} {u_{\textrm{LR}}(r_e)} \exp\left\{ -\beta_{\rm MLR}(r) y_p^{\textrm{eq}}(r)\right\}\right)^2,`

where the radial variable :math:`y_p^{\rm eq}` in the exponent, the long-range potential :math:`u_{\textrm{LR}}(r)` by
:math:`V(r)\to u_{\rm LR}(r) = \sum_{n} \frac{C_n}{r^n}` while the exponent coefficient function

:math:`\beta_{\rm MLR}(r) = y_p^{\rm{ref}}(r)\, \beta_{\infty}  +  \left[1 -y_p^{\textrm{ref}}(r)\right] \sum_{i=0} a_i[y_q^{\textrm{ref}}(r)]^i`

is defined in terms of two radial variables which are similar to :math:`y_p^{\rm eq}`, 
but are defined with respect to a different expansion center
`r_\textrm{ref}`, and involve two different powers, :math:`p` and :math:`q`. The above
definition of the function :math:`\beta_{\rm MLR}(r)` means that:

:math:`\beta_{\rm MLR}(r\to\infty)  \equiv  \beta_{\infty}  =  \ln[2D_{\rm e}/u_{\textrm{LR}}(r_{\rm e})].`


Example:
::

   poten 6
   name "d 3Pig"
   symmetry g
   lambda 1
   mult   3
   units bohr cm-1
   type  MLR
   values
     Te        0.20151357236994E+05
     RE        0.12398935933004E+01
     AE        0.50960000000000E+05        link   1   1   3
     RREF     -0.10000000000000E+01
     P         0.40000000000000E+01
     NL        0.20000000000000E+01
     NR        0.80000000000000E+01
     b0        0.30652655627150E+01
     b1       -0.93393246763924E+00
     b2        0.45686541184906E+01
     b3       -0.37637923145046E+01
     b4       -0.41028177891391E+01
     b5        0.00000000000000E+00
     b6        0.00000000000000E+00
     b7        0.00000000000000E+00
     b8        0.00000000000000E+00
     Binf      1.000000000000000000
     a1        0.00000000000000E+00
     a2        0.00000000000000E+00
     a3        0.00000000000000E+00
     a4        0.00000000000000E+00
     a5        0.00000000000000E+00
     a6        192774.
     a7        0.00000000000000E+00
     a8        0.00000000000000E+00
   end



Surkus-polynomial expansion ``Surkus``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(alias ``BobLeroy``)

:math:`V(r) = T_{\rm e} + (1-y_p^{\textrm{eq}}) \sum_{i\ge 0} a_i [y_p^{\textrm{eq}}]^i + y_p^{\textrm{eq}} a_{\rm inf},`


where :math:`y_p^{\textrm{eq}}` is the \v{S}urkus variable (\ref{eq:ypEQ}) and
`a_{\rm inf}` is the asymptote of the potential at :math:`r\to \infty`.


Example:
::

    spin-orbit  2 2
    name "<Lambda=+1,S=1 (a2Pi)|LSZ|+1 (a2Pi),S=1>"
    spin   0.5 0.5
    lambda 1 1
    sigma  0.5 0.5
    type  BOBLEROY
    units  cm-1
    factor    1.0   (0, 1 or i)
    values
     re         0.17700000000000E+01
     rref      -0.10000000000000E+01
     P          0.20000000000000E+01
     NT         0.30000000000000E+01
     a0        -0.63452015232176E+02
     a1        -0.20566444179565E+01
     a2        -0.13784613913938E+02
     a3         0.00000000000000E+00
     ainf      -0.56030500000000E+02
    end



``Surkus-damp`` (alias ``BobLeroy_damp``) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Surkus-polynomial expansion with a damping  function:

:math:`V(r) =  T_{\rm e} + \left[ (1-y_p^{\textrm{eq}}) \sum_{i\ge 0} a_i [y_p^{\textrm{eq}}]^i + y_p^{\textrm{eq}} a_{\rm inf}\right] f^{\rm damp} + t^{\rm damp} (1- f^{\rm damp}),`

where the damping function is defined by
:math:`f^{\rm damp} = 1-\tanh[\alpha(r-r_0)]`, and  :math:`t^{\rm damp}`, :math:`r_0` and :math:`\alpha` are parameters.


Example:
::

    spin-orbit  2 2
    name "<Lambda=+1,S=1 (a2Pi)|LSZ|+1 (a2Pi),S=1>"
    spin   0.5 0.5
    lambda 1 1
    sigma  0.5 0.5
    type  BOBLEROY
    units  cm-1
    factor    1.0   (0, 1 or i)
    values
    re         0.17700000000000E+01
    rref      -0.10000000000000E+01
    P          0.20000000000000E+01
    NT         0.30000000000000E+01
    a0        -0.63452015232176E+02
    a1        -0.20566444179565E+01
    a2        -0.13784613913938E+02
    a3         0.00000000000000E+00
    ainf      -0.56030500000000E+02
    tdamp      0.00000000000000E+00
    r0         0.10000000000000E+01
    alpha      0.30000000000000E+01
    end


Morse oscillator ``Morse`` 
^^^^^^^^^^^^^^^^^^^^^^^^^^

A polynomial expansion in the  Morse variable :math:`y_{\rm M}=1-e^{-a(r-r_0)}` is used

:math:`V(r)=T_{\rm e}+ (A_{\rm e}-T_{\rm e})  y_{\rm M}^2 +  \sum_{i=1}^N a_i \, y_{\rm M}^{i+2}.`



Example
::

    poten 1
    name "X 1Sigmag+"
    symmetry g +
    type   MORSE
    lambda 0
    mult   1
    units bohr cm-1
    values
    TE             0.00000000000000E+00
    RE             0.12423216077595E+01
    a              0.20372796052933E+01
    AE             0.73955889175514E+05
    A1            -0.62744302960091E+04
    A2            -0.57683579529693E+04
    end




``Morse_damp``
^^^^^^^^^^^^^^

:math:`V_(r)=T_{\rm e}+ (A_{\rm e}-T_{\rm e})  y_{\rm M}^2  + e^{-d_{\rm damp} (r-r_{\rm e})^4} \sum_{i=1} a_i  \left( \frac{r-r_{\rm e}}{r+r_{\rm e}} \right)^{i+2}.`

Example:
::

    poten 6
    name "d 3Pig"
    symmetry g
    lambda 1
    mult   3
    units bohr cm-1
    type  Morse_damp
    values
     Te      20121.09769
     re      0.12545760270976E+01
     Ae      0.50937907750000E+05        link   1   1   3
     a0      0.30398932686950E+01
     DAMP    0.10000000000000E-02
     a1      0.11437702960146E+05
     a2     -0.36585731834570E+03
     a3     -0.20920472718062E+05
     a4      0.90487097982036E-03
     a5      0.00000000000000E+00
     a6      0.00000000000000E+00
     a7      0.00000000000000E+00
     a8      0.00000000000000E+00
    end



``Modified-Morse``
^^^^^^^^^^^^^^^^^^

Alias ``MMorse``

:math:`V_(r)=T_{\rm e}+ (A_{\rm e}-T_{\rm e}) \frac{ \left[ 1-\exp\left(-\sum_{i=0} a_i \xi^{i+1}\right)  \right]^2}{\left[ 1-\exp\left(-\sum_{i=0} a_i \right) \right]^2},`

where  :math:`\xi = (r-r_{\rm e})/(r+r_{\rm e})`.

Example:
::

    poten 8
    name "Bp 1Sigmag+"
    symmetry g +
    lambda 0
    mult   1
    type  MMorse
    values
    Te            1.5408840263E+04
    rE            1.3778208709E+00
    Ae            5.0937907750E+04               link   1   1   3
    a0            6.2733066935E+00
    a1            1.4954972843E+01
    a2            4.5160872659E+01
    end

where the value :math:`A_{\rm e}` is `linked` to the corresponding value of ``poten 1``.

``Polynomial`` 
^^^^^^^^^^^^^^

This keyword selects a polynomial expansion in the variable :math:`y=(r-r_0)`

:math:`V(r) = T_{\rm e} + a_1 y + a_2 y^2 + \cdots`


Example:
::

    spin-orbit  2 2
    name "<+1,S=1 (a3Pi)|LSZ|+1  (a3Pi),S=1>"
    spin   1.0 1.0
    sigma  1.0 1.0
    lambda 1 1
    type  polynom
    factor   1
    values
    a0           14.97
    re           1.3
    a1           0.0
    end


**Taylor expansion around** :math:`r_0`:
:math:`V(r) = T_{\rm e} + a_1 (r-r_0)^2 + a_2 (r-r_0)^2 + a_3 (r-r_0)^3  + \cdots` 

``Dunham`` expansion 

``Dunham`` selects a polynomial expansion in the Dunham variable  :math:`y=(r-r_0)/r_0` 

:math:`V(r) = T_{\rm e}+ a_0 y^2 \left( 1 + a_1 y + a_2 y^2 + \cdots \right)`

Example:
::

    poten 1
    name "X 2 Delta"
    lambda 2
    mult   2 type   Dunham values
    Te              0.00000
    Re              1.4399282269779912
    a0         123727.20496894409      (= omega**2 / 4 B)
    a2             -2.31
    a3              3.80
    a4             -6.00
    a5              5.00
    end


.. 
   As a function form ``Dunham`` is equivalent to a ``Polynomial`` object with the linear term absent and 
   a redefinition of the expansion coefficients; the comments given for ``Polynomial`` also apply to ``Dunham``.

**Taylor expansion around** :math:`r_0`:
:math:`V(r) = T_{\rm e} + \frac{a_0}{r_0^2} (r-r_0)^2 + \frac{a_0 a_1}{r_0^3} (r-r_0)^3 + \cdots`

Simons, Parr and Finlan ``SPF``  
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``SPF`` selects a polynomial expansion in the the so-called Simons, Parr and Finlan variable :math:`y=(r-r_0)/r` 

:math:`V(r) = T_{\rm e} + a_0 y^2 \left( 1 + a_1 y + a_2 y^2 + \cdots \right)`


Example:
::

    poten 1
    name "X 2Sigma+"
    symmetry +
    type   SPF
    lambda 0
    mult   2
    values
    Te         0.00000000000000E+00
    RE         0.16292698613903E+01
    a1         0.37922070444743E+06
    a2         0.00000000000000E+00
    a3        -0.53314483965665E+01
    a4         0.00000000000000E+00
    a5         0.19407192336518E+02
    a4         0.00000000000000E+00
    a5        -0.17800496953835E+03
    end


**Taylor expansion around** :math:`r_0`:
:math:`V(r) = T_{\rm e} + \frac{a_0}{r_0^2} (r-r_0)^2 + \frac{a_0 a_1 - 2 a_0}{r_0^3} (r-r_0)^3 + \cdots`

.. 
  Behaviour for :math:`r \to +\infty`:

:math:`V(r) = T_{\rm e} + a_0 \left(1+\sum_{i=1}^N a_i\right) - \frac{a_0 r_0}{r} \left( 2+\sum_{i=1}^N (i+2)a_i \right ) + \cdots`

where :math:`N` is the maximum exponent included in the expansion.
For long :math:`r` the potential goes to a constant value; convergence to the constant
is of the :math:`1/r` type (correct for ions but too slow for neutral molecules).

**Behaviour for** :math:`r \to 0`

:math:`V(r) = a_0 a_N \left(\frac{r_0}{r}\right)^{N+2} + \cdots`

The coefficient :math:`a_0` is definitely positive, but :math:`a_N` can be positive and negative, 
so that :math:`V(r)` can go to :math:`\pm \infty` for short :math:`r`.

Murrell-Sorbie ``M-S``
^^^^^^^^^^^^^^^^^^^^^^^^

:math:`V(r)=A_{\rm e}- (A_{\rm e}-T_{\rm e})\left( 1 + a_1 \rho + a_2 \rho^2 + a_3 \rho^3 + \ldots  \right) e^{-a_1 \rho},`
`
where  :math:`\rho = r-r_{\rm e}`.

Example:
::

   poten 4
   name "B 2Sigma"
   symmetry -
   type  M-S  (Murrell-Sorbie)
   lambda 0
   mult   2
   values
   V0            21000.0
   RE            1.6
   DE            25653.27131
   a1            2.81468
   a2            1.68719
   a3            0.757787
   a4            -0.5963168
   a5            -0.54596343
   a6            0.20611664
   end



**Taylor expansion around** :math:`r_0`:
:math:`V(r) = T_{\rm e} + \frac{A_\mathrm{e}-T_\mathrm{e}}{2} (a_1^2 - 2a_2) (r-r_0)^2 + \frac{A_\mathrm{e}-T_\mathrm{e}}{3} (a_1^3 -3a_1 a_2+3 a_3) (r-r_0)^3 + \cdots`


**Behaviour for** :math:`r \to +\infty`:
:math:`V(r) = A_{\rm e} - a_N (A_\mathrm{e}-T_\mathrm{e}) (r-r_e)^N e^{-a_1 (r-r_e)} + \cdots`
`
where :math:`N` is the maximum exponent included in the expansion.
For long :math:`r` the potential goes to the constant value :math:`A_\mathrm{e}`, and the aymptotic behavior is
determined by the coefficients of the term with the highest exponent.

``Chebyshev`` 
^^^^^^^^^^^^^

This keyword selects an expansion in Chebyshev polynomials in the variable 
:math:`y= [r-(b+a)/2]/[(b-a)/2]`. The scaled variable :math:`y` ranges from :math:`-1` to 1 for :math:`r`     
in :math:`[a,b]`. The expansion is  

:math:`V(r) = a_0 + a_1 T_1(y) + a_2 T_2(y) + \cdots`

Example:
::

    spin-orbit  2 2
    name "<+1,S=1 (a3Pi)|LSZ|+1  (a3Pi),S=1>"
    spin   1.0 1.0
    type  chebyshev
    factor   1
    values
       a               0.80000000000000E+00
       b               0.26500000000000E+01
       A0             -0.25881057805341E+02
       A1              0.82258425882627E+01
       A2              0.52391700137878E+00
       A3              0.28483394288286E+01
       A4             -0.15136422837793E+00
       A5              0.97553692867070E-01
       A6             -0.25825811071417E+00
       A7             -0.69105144347567E-01
       A8             -0.44700771508442E-01
       A9              0.11793957297111E-01
       A10             0.16403055376257E-01
       A11             0.92509900186428E-02
       A12             0.50789943150707E-02
       A13            -0.39439903216016E-03
    end


    
``COSH-POLY`` 
^^^^^^^^^^^^^

This function can be used as a coupling for a diabatic representation of potentials characterised by
an avoiding crossing and is given by:
:math:`F(r)= F_0 + \frac{ \sum_{i=0}^N a_i \, (r-r_{\rm ref})^{i}.}{\cosh\beta(r-r_{\rm ref})} .`


Example
::

    diabatic  1 8
    name "<X1Sigmag+|D|Bp 1Sigmag+>"
    spin   0.0 0.0
    lambda  0  0
    type  COSH-poly
    units  cm-1
    factor    i   (0, 1 or i)
    values
    v0            0.0000
    beta          5.62133
    RE            1.610505
    B0           -0.307997
    B1            0.0000000000E+00
    B2            0.0000000000E+00
    BINF          0.0000000000E+00
    end




``POLYNOM_DECAY_24`` 
^^^^^^^^^^^^^^^^^^^^

This function is similar to ``Surkus`` expansion
:math:`F(r)=\sum^{N}_{k=0}B_{k}\, z^{k} (1-\xi_p) + \xi_p\, B_{\infty},`

where :math:`z` is either taken as the damped-coordinate given by:

:math:`z = (r-r_{\rm ref})\, e^{-\beta_2 (r-r_{\rm ref})^2-\beta_4 (r - r_{\rm ref})^4},`

Here :math:`r_{\rm ref}` is a reference position equal to :math:`r_{\rm e}` by default and 
:math:`\beta_2` and :math:`\beta_4` are damping factors. 
When used for morphing, the parameter :math:`B_{\infty}` is usually fixed to 1.


``POLYNOM_DIMENSIONLESS`` 
^^^^^^^^^^^^^^^^^^^^^^^^^

This function is a polynomial 
:math:`F(r)=\sum^{N}_{k=0} a_{k}\, y^{k} ,`
in terms of the dimensionless variable 
:math:`y = \frac{r-r_{\rm e}}{r_{\rm e}}.`

The order of the parameters in the input is as follows :math:`a_0,r_{\rm e}, a_1,a_2, \ldots`

Example
::

  dipole 1 1 
  name "L_2015"
  type POLYNOM_DIMENSIONLESS 
  spin   0.0 0.0
  lambda  0  0
  values
   re   1.12832252847d0     
   a0   -0.1229099d0
   a1    3.604742d0
   a2   -0.23716d0
   a3   -3.67326d0
   a4    1.4892d0 
   a5    1.8293d0 
   a6   -4.342d0  
  end



``CO_X_UBOS`` 
^^^^^^^^^^^^^

This CO PEC was used in `Meshkov et. al, JQSRT, 217, 262 (2017) <https://doi.org/10.1016/j.jqsrt.2018.06.001>`_ to compute energies 
of CO in its ground electronic state.  All parameters are predefined internally.  




Implementation guide  
^^^^^^^^^^^^^^^^^^^^

All these analytical functions are programmed as Fortran double precision functions 
in the module ``functions.f90``. 

Below is an example of a function for the `EMO` potential energy function. 
::

    function poten_EMO(r,parameters) result(f)
      !
      real(rk),intent(in)    :: r             ! geometry (Ang)
      real(rk),intent(in)    :: parameters(:) ! potential parameters
      real(rk)               :: y,v0,r0,de,f,rref,z,phi
      integer(ik)            :: k,N,p
      !
      v0 = parameters(1)
      r0 = parameters(2)
      ! Note that the De is relative the absolute minimum of the ground state
      De = parameters(3)-v0
      !
      rref = parameters(4)
      !
      if (rref<=0.0_rk) rref = r0
      !
      if (r<=rref) then 
        p = nint(parameters(5))
        N = parameters(7)
      else
        p = nint(parameters(6))
        N = parameters(8)
      endif 
      !
      if (size(parameters)/=8+max(parameters(7),parameters(8))+1) then 
        write(out,"('poten_EMO: Illegal number of parameters in EMO, check NS and NL, must be max(NS,NL)+9')")
        print*,parameters(:)
        stop 'poten_EMO: Illegal number of parameters, check NS and NL'
      endif 
      !
      z = (r**p-rref**p)/(r**p+rref**p)
      !
      phi = 0
      do k=0,N
       phi = phi + parameters(k+9)*z**k
      enddo
      !
      y  = 1.0_rk-exp(-phi*(r-r0))
      !
      f = de*y**2+v0
      !
    end function poten_EMO


To define a new functional form, apart from the actual function, a new reference ``case`` identifying this calculation 
options needs to be added as part of the ``case select`` section in the ``subroutine define_analytical_field``, for example:
::

    case("EMO") ! "Expanded MorseOscillator"
      !
      fanalytical_field => poten_EMO


