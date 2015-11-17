module functions
  !
  use accuracy
  use timer
  !
  implicit none
  !
  ! Different 1D functions (potential, spin-orbit etc) 
  !
  integer(ik),parameter :: verbose=5
  !
  public analytical_fieldT
  !
  !
  abstract interface
    !
    function analytical_fieldT(r,parameters)
      use accuracy
      !
      real(rk)                      :: analytical_fieldT !NB: NAG Fortran 6.0 doesn't like in-line declaration
      real(rk),intent(in)           :: r             ! geometry (Ang)
      real(rk),intent(in)           :: parameters(:) ! potential parameters
      !
    end function analytical_fieldT
    !
  end interface 
  !
  !
  contains
  !
  !
  subroutine define_analytical_field(ftype,fanalytical_field)
    !
    character(len=cl),intent(in)      :: ftype
                                            ! NB: NAG Fortran 6.0 doesn't like intent(in) and initial nullification
    procedure(analytical_fieldT),pointer :: fanalytical_field !=> null()
    !
    select case(ftype)
      !
    case("MORSE")
      !
      fanalytical_field => poten_morse
      !
    case("MORSE_DAMP")
      !
      fanalytical_field => poten_morse_damp
      !
    case("MODIFIED-MORSE","MODIFIED_MORSE","MMORSE")
      !
      fanalytical_field => poten_morse_modified
      !
    case("EMO") ! "Expanded MorseOscillator"
      !
      fanalytical_field => poten_EMO
      !
    case("MLR") ! "Morse/Long-Range"
      !
      fanalytical_field => poten_MLR
      !
    case("MARQUARDT") ! "Marquardt"
      !
      fanalytical_field => poten_Marquardt
      !
    case("COSH-POLY") ! "Diabatic coupling as a polynom/cosh"
      !
      fanalytical_field => poten_cosh_polynom
      !
    case("BOBLEROY","BOB","SURKUS") ! "BOB expansion"
      !
      fanalytical_field => poten_BOBLeRoy
      !
    case("BOBLEROY_DAMP","BOB_DAMP","SURKUS_DAMP") ! "BOB expansion with damping to zero at r=0"
      !
      fanalytical_field => poten_BOBLeRoy_damp
      !
    case("DUNHAM")
      !
      fanalytical_field => poten_dunham
      !
    case("SPF")
      !
      fanalytical_field => poten_spf
      !
    case("SPF_1")
      !
      fanalytical_field => poten_spf_1
      !
    case("SPF_H2")
      !
      fanalytical_field => poten_spf_h2
      !
    case("HH","HULBERT-HIRSCHFELDER")
      !
      fanalytical_field => poten_Hulbert_Hirschfelder
      !
    case("CHEBYSHEV")
      !
      fanalytical_field => poten_cheb
      !
    case("POLYNOM", "POLYNOMIAL")
      !
      fanalytical_field => poten_polynom
      !
    case("EXPONENTIAL")
      !
      fanalytical_field => poten_exponential
      !
    case("POTEN_LAURA")
      !
      fanalytical_field => poten_laura
      !
    case("POTEN_FERMI_T")
      !
      fanalytical_field => poten_fermi_t
      !
    case("M-S")
      !
      fanalytical_field => poten_Murrell_Sorbie
      !
    case("NONE")
      !
      write(out,'("Analytical: Some fields are not properly defined and produce function type ",a)') trim(ftype)
      stop "Analytical: Unknown function type "
      !
    case default
      !
      write(out,'("Analytical: Unknown function type ",a)') trim(ftype)
      stop "Analytical: Unknown field type "
      !
    end select
    !
  end subroutine define_analytical_field
  !
  !
  function poten_morse(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0,a0,f,de
    integer(ik)            :: k,N
    !
    N = size(parameters)
    !
    v0 = parameters(1)
    r0 = parameters(2)
    a0 = parameters(3)
    !
    y = 1.0_rk-exp(-a0*(r-r0))
    ! Note that the De is relative the absolute minimum of the ground state
    !
    f = v0
    !
    de = parameters(4)-v0
    !
    f = f + de*y**2
    !
    do k=5,N 
      f = f + parameters(k)*y**(k-2)
    enddo
    !
  end function poten_morse
  !
  function poten_morse_damp(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0,a0,f,de,vdamp,vshort,vlong,damp,dr,z
    integer(ik)            :: k,N,Npot,Nstruct
    !
    N = size(parameters)
    !
    v0 = parameters(1)
    r0 = parameters(2)
    de = parameters(3)-v0
    a0 = parameters(4)
    damp = parameters(5)
    !
    ! number of structural parameters
    !
    Nstruct = 5
    !
    Npot = N-Nstruct
    !
    dr = (r-r0)
    !
    y = 1.0_rk-exp(-a0*dr)
    !
    ! Note that the De is relative the absolute minimum of the ground state
    !
    vlong = de*y**2
    !
    vdamp = exp(-damp*dr**4)
    !
    vshort = 0
    !
    z = (r-r0)/(r+r0)
    !
    do k=1,Npot
      vshort = vshort + parameters(Nstruct+k)*y**(k+2)
    enddo
    !
    f = v0+vdamp*vshort+vlong
    !
  end function poten_morse_damp
  !
  function poten_morse_modified(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0,de,f,f0,f1,beta,beta_inf,z,y_inf
    integer(ik)            :: k,N
    !
    N = size(parameters)
    !
    v0 = parameters(1)
    r0 = parameters(2)
    !
    ! Note that the De is relative the absolute minimum of the ground state
    De = parameters(3)-v0
    !
    z = (r-r0)/(r+r0)
    !
    f0 = 0
    f1 = 0
    do k=4,N
     f0 = f0 + parameters(k)*z**(k-4)
     f1 = f1 + parameters(k) ! *(1.0_rk)**(k-4)
    enddo
    !
    beta = z*f0
    beta_inf = f1
    !
    y     = 1.0_rk-exp(-beta)
    y_inf = 1.0_rk-exp(-beta_inf)
    !
    f = de*y**2/y_inf**2+v0
    !
  end function poten_morse_modified
  !
  function poten_cosh_polynom(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,r0,f,f0,beta,z,v0
    integer(ik)            :: k,N
    !
    N = size(parameters)
    !
    v0 = parameters(1)
    !
    beta = parameters(2)
    r0 = parameters(3)
    !
    z = (r-r0)
    !
    f0 = 0
    do k=4,N
     f0 = f0 + parameters(k)*z**(k-4)
    enddo
    !
    y     = cosh(beta*z)
    !
    f = v0+f0/y
    !
  end function poten_cosh_polynom
  !
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
  !
  ! Morse/Long-Range, see Le Roy manuals
  !
  function poten_MLR(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0,de,f,rref,z,phi,phiinf,uLR,uLR0
    integer(ik)            :: k,N,p,M,Nstruc,Ntot,Npot
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
    p = nint(parameters(5))
    !
    if (r<=rref) then 
      N = parameters(6)
    else
      N = parameters(7)
    endif
    !
    ! Number of structural parameters 
    !
    Nstruc = 7
    !
    ! total number of parameters 
    !
    Ntot = size(parameters(:),dim=1)
    !
    ! Number of pot-parameters
    !
    Npot = max(parameters(6),parameters(7))+1
    !
    ! number of long range parameters 
    !
    M = Ntot-Npot-Nstruc
    !
    !if (size(parameters)/=8+max(parameters(7),parameters(8))+1) then 
    !  write(out,"('poten_EMO: Illegal number of parameters in EMO, check NS and NL, must be max(NS,NL)+9')")
    !  stop 'poten_EMO: Illegal number of parameters, check NS and NL'
    !endif 
    !
    z = (r**p-rref**p)/(r**p+rref**p)
    !
    phi = 0
    do k=0,N
     phi = phi + parameters(k+Nstruc+1)*z**k
    enddo
    !
    ! double check uLR parameters
    !
    uLR = sum(parameters(1+Npot+Nstruc:M+Npot+Nstruc)**2)
    !
    if (uLR<small_) then
      write(out,"('poten_MLR: At least one uLR should be non-zero')")
      stop 'poten_MLR: At least one uLR should be non-zer'
    endif
    !
    ! long-range part
    !
    uLR = 0
    do k=1,M
     uLR = uLR + parameters(k+Npot+Nstruc)/r**k
    enddo
    !
    ! at R=Re
    uLR0 = 0
    do k=1,M
     uLR0 = uLR0 + parameters(k+Npot+Nstruc)/r0**k
    enddo
    !
    phiinf = log(2.0_rk*de/uLR0)
    !
    phi = (1.0_rk-z)*phi+z*phiinf
    !
    phi = uLR/uLR0*exp(-phi*z)
    ! 
    y  = 1.0_rk-phi
    !
    f = de*y**2+v0
    !
  end function poten_MLR
  !
  function poten_Marquardt(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0,de,f,rref,z,phi,eps6,eps8,rs,damp
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
    eps6 = parameters(9)
    eps8 = parameters(10)
    rs = parameters(11)
    !
    z = (r**p-rref**p)/(r**p+rref**p)
    !
    phi = 0
    do k=0,N
     phi = phi + parameters(k+12)*z**k
    enddo
    !
    y  = 1.0_rk-exp(-phi*(r-r0))
    !
    damp = ( 1.0_rk+eps6*(-(rs/r)**6) )*( 1.0_rk+eps8*(-(rs/r)**8) ) 
    y = y*damp
    !
    f = de*y**2+v0
    !
  end function poten_Marquardt
  !
  function poten_BOBLeRoy(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: r0,f,rref,z,t,tinf
    integer(ik)            :: k,N,p
    !
    r0 = parameters(1)
    !
    rref = parameters(2)
    !
    p = nint(parameters(3))
    N = size(parameters)-4-2
    !
    !if (size(parameters)/=4+N+2) then 
    !  write(out,"('poten_BOBLeRoy: Illegal number of parameters, check N')")
    !  stop 'poten_BOBLeRoy: Illegal number of parameters, check N'
    !endif 
    !
    z = (r**p-r0**p)/(r**p+r0**p)
    !
    t = 0
    do k=0,N
     t = t + parameters(k+5)*z**k
    enddo
    !
    tinf = parameters(N+6)
    !
    f = (1.0_rk-z)*t+z*tinf
    !
  end function poten_BOBLeRoy
  !
  function poten_BOBLeRoy_damp(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: r0,f,rref,z,t,tinf,f_damp,a_d,r_d,t_0
    integer(ik)            :: k,N,p
    !
    r0 = parameters(1)
    !
    rref = parameters(2)
    !
    p = nint(parameters(3))
    N = size(parameters)-4-2-3
    !
    !if (size(parameters)/=4+N+2) then 
    !  write(out,"('poten_BOBLeRoy: Illegal number of parameters, check N')")
    !  stop 'poten_BOBLeRoy: Illegal number of parameters, check N'
    !endif 
    !
    z = (r**p-r0**p)/(r**p+r0**p)
    !
    t = 0
    do k=0,N
     t = t + parameters(k+5)*z**k
    enddo
    !
    tinf = parameters(N+6)
    !
    t_0 = parameters(N+7)
    r_d = parameters(N+8)
    a_d = parameters(N+9)
    !
    !f = (1.0_rk-z)*t+z*tinf
    !
    f_damp = ( 1.0_rk+tanh(a_d*(r-r_d)) )*0.5_rk
    !
    f = ( (1.0_rk-z)*t+z*tinf )*f_damp+t_0*(1.0_rk-f_damp)
    !
  end function poten_BOBLeRoy_damp
  !
  function poten_dunham(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0,a0,f
    integer(ik)            :: k,N
    !
    N = size(parameters)
    !
    v0 = parameters(1)
    r0 = parameters(2)
    a0 = parameters(3)
    !
    y = (r-r0)/r0
    !
    f = 1.0_rk
    do k=4,N 
      f = f + parameters(k)*y**(k-3)
    enddo
    !
    f = f*a0*y**2 + v0
    !
  end function poten_dunham
  !
  function poten_spf(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0,b0,f
    integer(ik)            :: k,N
    !
    N = size(parameters)
    !
    v0 = parameters(1)
    r0 = parameters(2)
    b0 = parameters(3)
    !
    y = (r-r0)/r
    !
    f = 1.0_rk
    do k=4,N 
      f = f + parameters(k)*y**(k-3)
    enddo
    !
    f = f*b0*y**2 + v0
    !
  end function poten_spf
  !
  function poten_spf_1(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0, f !, b0
    integer(ik)            :: k,N
    !
    N = size(parameters)
    !
    v0 = parameters(1)
    r0 = parameters(2)
    !
    y = (r-r0)/r
    !
    f = 0.0_rk
    do k=3,N 
      f = f + parameters(k)*y**(k-2)
    enddo
    !
    f = f + v0
    !
  end function poten_spf_1
  !
  function poten_spf_h2(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0,b0,f,vh2
    integer(ik)            :: k,N
    !
    N = size(parameters)
    !
    v0 = parameters(1)
    r0 = parameters(2)
    b0 = parameters(3)
    !
    y = (r-r0)/r
    !
    f = 1.0_rk
    do k=4,N-4 
      f = f + parameters(k)*y**(k-3)
    enddo
    !
    f = f*b0*y**2 + v0
    !
    vh2 = 0
    !
    do k=1,4
       vh2 = vh2+parameters(N-4+k)*y**k
    enddo
    !
    f = f + vh2
    !
  end function poten_spf_h2
  !
  !
  ! SO functional form  
  !
  function poten_Hulbert_Hirschfelder(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: x,De,Aatom,re,a,b,c,f,UHH
    !
    if (size(parameters)/=6) stop 'Illegal number of parameters in poten_Hulbert_Hirschfelder (must be 6)'
    !
    ! Aatom is the atomic limit
    ! De is ASO(eq)-A(atom)
    ! re is the "equilibrium"
    ! a,b,c are the expansion parameters
    !
    De    = parameters(1)
    Aatom = parameters(2)
    re   = parameters(3)
    a    = parameters(4)
    b    = parameters(5)
    c    = parameters(6)
    !
    x = a*(r/re-1.0_rk)
    !
    ! Note that the De is difference between SO at equilibrium and SO(atomic limit)
    !
    UHH = De*( ( 1.0_rk-exp(-x) )**2+c*exp( -2.0_rk*x )*x**3*( 1.0_rk+b*x ) )
    !
    f = Aatom-De+UHH
    !
  end function poten_Hulbert_Hirschfelder


  !
  function poten_Murrell_Sorbie(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0,a0,f,de,f_t,f_exp
    integer(ik)            :: k,N
    !
    N = size(parameters)
    !
    v0 = parameters(1)
    r0 = parameters(2)
    a0 = parameters(4)
    !
    y = (r-r0)
    ! Note that the De is relative the absolute minimum of the ground state
    !
    de = parameters(3)-v0
    !
    f_exp = exp(-a0*y)
    !
    f_t = 1.0_rk
    !
    do k=4,N
      f_t = f_t + parameters(k)*y**(k-3)
    enddo
    !
    f = v0-de*f_t*f_exp
    !
  end function poten_Murrell_Sorbie

  !
  function poten_polynom(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0,f
    integer(ik)            :: k,N
    !
    N = size(parameters)
    !
    v0 = parameters(1)
    r0 = parameters(2)
    !
    y = (r-r0)
    !
    f = v0
    do k=3,N 
      f = f + parameters(k)*y**(k-2)
    enddo
    !
  end function poten_polynom
  !
  ! The Chebyshev polynomial is evaluated at a point y=[x-(b + a)/2]/[(b - a)/2],
  ! and the result is returned as the function value.
  !
  real(rk) function poten_cheb(x,c)
    real(rk),intent(in) :: x,c(:)
    !
    integer(ik) :: m
    real(rk) :: a,b
    !
    integer(ik) ::  j
    real(rk)    ::  d,dd,sv,y,y2
      !
      m = size(c)
      !
      a = c(1) ; b = c(2)
      if ((x-a)*(x-b).gt.0.) stop 'x not in range in chebev'
      d=0
      dd=0
      y=(2.0_rk*x-a-b)/(b-a) ! Change of variable.
      y2=2.0_rk*y
      do j=m,4,-1  ! Clenshaw's recurrence.
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
      enddo
      !
      poten_cheb=y*d-dd+0.5_rk*c(3) ! Last step is different.
      !
  end function poten_cheb
  !
 function poten_exponential(r, parameters) result(fun)
  ! expansion in powers of y = e^(-a(r-re)), 
  !   F = sum_{n=0}^m c_n y^n = c(1) + c(3)*exp( -a*r ) + ...
  !   where a=c(2) and re=c(3)
  ! One should provide either one parameter or at least 3.
    real(rk),intent(in) :: r, parameters(:)
    real(rk)            :: fun
    !
    integer(ik) :: m, j
    real(rk)    :: a
    real(rk)    :: y
      !
      m = size(parameters)
      !
      ! do some checks
      if( m < 1) stop 'Too few parameters in poten_exponential: m < 1 '
      fun = parameters(1)
      if( m ==1) return ! simple case of a constant 
      if( m ==2) stop 'Too few parameters in poten_exponential: m ==2 '
      a  = parameters(2)

      y=exp( -a*r )

      do j=1,m-2
        fun = fun + parameters(j+2) * y**j
      enddo
      !
  end function poten_exponential
  !
  !
  function poten_laura(r,parameters) result(fun)
    implicit none
    real(rk),intent(in) :: r,parameters(:)
    integer(ik)  :: m
    real(rk)            :: fun,v0,r0,lambda,qm,j,q
    m=size(parameters)
    !
    if(m/=5) stop 'poten_laura: illegal number of parameters parameters (/=5)'
    !
    v0     = parameters(1)
    r0     = parameters(2)
    lambda = parameters(3)
    qm     = parameters(4)
    j      = parameters(5)
    !
    q = (r-r0)/qm
    !
    fun = 0.25_rk*q**2+J*( 1.0_rk-sqrt(1.0_rk+( 0.5_rk*lambda/J*q )**2 ) )
    !
  end function poten_laura
  !
  function poten_fermi_T(r,parameters) result(fun)
    implicit none
    real(rk),intent(in) :: r,parameters(:)
    integer(ik)  :: m
    real(rk)            :: fun
    m=size(parameters)
    if(m/=7) stop 'poten_fermi_T: illegal number of parameters parameters (/=7)'
    fun = parameters(1) +parameters(2)/(parameters(3)+&
          exp((-parameters(4)+parameters(8))/(parameters(5)*r**2+parameters(6)*r+parameters(7))))
  end function poten_fermi_T
  !
  function poten_fermi_deriv(r,parameters) result(fun)
    implicit none
    real(rk), intent(in) :: r,parameters(:)
    integer(ik) :: m
    real(rk)            :: fun
    !
    m=size(parameters)
    if(m<4) stop 'poten_fermi_deriv: Too few parameters'
    if(m>4) stop 'poten_fermi_deriv: Too many parameters'
    fun = -parameters(1)*exp((-parameters(2)+r)/parameters(3))/(parameters(3)*(parameters(4)+&
          exp((-parameters(2)+r)/parameters(3)))**2)
    !
  end function poten_fermi_deriv
  !
end module functions
