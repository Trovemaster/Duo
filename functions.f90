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
      real(rk),intent(in)           :: r                 ! geometry (Ang)
      real(rk),intent(in)           :: parameters(:)     ! potential parameters
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
    case("MLR_DS") ! "Morse/Long-Range with Douketis-damping"
      !
      fanalytical_field => poten_MLR_Douketis
      !
    case("MLR_DS_DARBY") ! "Morse/Long-Range with Douketis-damping"
      !
      fanalytical_field => poten_MLR_Douketis_Darby
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
    case("BOBNA") ! "BOB-NA expansion"
      !
      fanalytical_field => poten_BOBna
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
    case("LAURA_SO","ATAN_SO")
      !
      fanalytical_field => SO_arctan
      !
    case("POTEN_FERMI_T")
      !
      fanalytical_field => poten_fermi_t
      !
    case("M-S")
      !
      fanalytical_field => poten_Murrell_Sorbie
      !
    case("PADE_GOODISMAN2","PADE2")
      !
      fanalytical_field => poten_Pade_Goodisman_2
      !
    case("POLYNOM_DECAY")
      !
      fanalytical_field => dipole_polynom_exp
      !
    case("POLYNOM_DECAY_DAMP")
      !
      fanalytical_field => dipole_polynom_exp_damp
      !
    case("POLYNOM_DECAY_24")
      !
      fanalytical_field => dipole_polynom_exp24
      !
    case("DOUBLEEXP2")
      !
      fanalytical_field => dipole_doubleexp
      !
    case("DIABATIC_MU_DIAG","AVOIDEDCROSSING_DIAG_MU")
      !
      fanalytical_field => dipole_avoidedcrossing_diag_mu
      !
    case("TWO_COUPLED_EMOS")
      !
      fanalytical_field => poten_two_coupled_EMOs
      !
    case("TWO_COUPLED_BOBS")
      !
      fanalytical_field => poten_two_coupled_BOBs
      !
    case("COUPLED_EMO_REPULSIVE")
      !
      fanalytical_field => poten_two_coupled_EMO_repulsive
      !
    case("REPULSIVE")
      !
      fanalytical_field => poten_repulsive
      !
    case("LORENTZ","LORENTZIAN")
      !
      fanalytical_field => poten_lorentzian_polynom
      !
    case("POLYNOM_DIMENSIONLESS","POLYNOMIAL_DIMENSIONLESS")
      !
      fanalytical_field => polynomial_dimensionless
      !
    case("CO_X_UBO")
      !
      fanalytical_field => potential_stolyarov_CO_X_UBO
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
 
 
   ! Morse/Long-Range with Douketis damping, see Le Roy manuals
  !
  function poten_MLR_Douketis(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0,de,f,rref,z,beta,betainf,betaN,yq,yp,uLR,uLR0,rho,b,c,s,damp
    integer(ik)            :: k,N,p,M,Nstruc,Ntot,Npot,q
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
    q = nint(parameters(6))
    N = parameters(7)
    rho = parameters(8)
    !
    ! Number of structural parameters 
    !
    Nstruc = 8
    !
    ! total number of parameters 
    !
    Ntot = size(parameters(:),dim=1)
    !
    ! Number of pot-parameters
    !
    Npot = N+1
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
    z = (r**p-r0**p)/(r**p+r0**p)
    yp = (r**p-rref**p)/(r**p+rref**p)
    yq = (r**q-rref**q)/(r**q+rref**q)
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
    ! For the Damping part   
    ! the values of s, b,c are as suggested by LeRoy 2011 (MLR paper)
    !
    s = -1.0_rk
    b = 3.3_rk
    c = 0.423_rk

    !s = -2.0_rk
    !b = 2.50_rk
    !c = 0.468_rk


    !
    ! long-range part
    !
    uLR = 0
    do k=1,M
     !
     ! Douketis damping function    
     Damp = ( 1.0_rk-exp( -b*rho*r/real(k,rk)-c*(rho*r)**2/sqrt(real(k,rk)) ) )**(real(k,rk)+s)
     !
     uLR = uLR + Damp*parameters(k+Npot+Nstruc)/r**k
     !
    enddo
    !
    !
    ! at R=Re
    uLR0 = 0
    do k=1,M
     ! Douketis damping function
     Damp = ( 1.0_rk-exp( -b*rho*r0/real(k,rk)-c*(rho*r0)**2/sqrt(real(k,rk)) ) )**(real(k)+s)
     !
     uLR0 = uLR0 + Damp*parameters(k+Npot+Nstruc)/r0**k
     !
    enddo


    !
    betaN = 0
    do k=0,N
     betaN = betaN + parameters(k+Nstruc+1)*yq**k
    enddo
    !
    betainf = log(2.0_rk*de/uLR0)
    !
    beta = (1.0_rk-yp)*betaN+yp*betainf
    ! 
    y  = 1.0_rk-uLR/uLR0*exp(-beta*z)
    !
    f = de*y**2+v0
    !
  end function poten_MLR_Douketis
 
 
 
   function poten_MLR_Douketis_Darby(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0,de,f,rref,z,beta,betainf,betaN,yq,yp,uLR,uLR0,rho,b,c,s,damp,u,uinf,ma,mb
    integer(ik)            :: k,N,p,M,Nstruc,Npot,q,NUa,NUb
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
    q = nint(parameters(6))
    N = parameters(7)
    rho = parameters(8)
    !
    ! Number of structural parameters 
    !
    Nstruc = 8
    !
    ! total number of parameters 
    !
    ! Number of pot-parameters
    !
    Npot = N+1
    !
    ! number of long range parameters 
    !
    M = parameters(Nstruc+Npot+1)
    !
    !if (size(parameters)/=8+max(parameters(7),parameters(8))+1) then 
    !  write(out,"('poten_EMO: Illegal number of parameters in EMO, check NS and NL, must be max(NS,NL)+9')")
    !  stop 'poten_EMO: Illegal number of parameters, check NS and NL'
    !endif 
    !
    z = (r**p-r0**p)/(r**p+r0**p)
    yp = (r**p-rref**p)/(r**p+rref**p)
    yq = (r**q-rref**q)/(r**q+rref**q)
    !
    ! double check uLR parameters
    !
    uLR = sum(parameters(2+Npot+Nstruc:M+Npot+Nstruc)**2)
    !
    if (uLR<small_) then
      write(out,"('poten_MLR: At least one uLR should be non-zero')")
      stop 'poten_MLR: At least one uLR should be non-zer'
    endif
    ! 
    ! For the Damping part   
    ! the values of s, b,c are as suggested by LeRoy 2011 (MLR paper)
    !
    s = -1.0_rk
    b = 3.3_rk
    c = 0.423_rk

    !s = -2.0_rk
    !b = 2.50_rk
    !c = 0.468_rk


    !
    ! long-range part
    !
    uLR = 0
    do k=1,M
     !
     ! Douketis damping function
     Damp = ( 1.0_rk-exp( -b*rho*r/real(k,rk)-c*(rho*r)**2/sqrt(real(k,rk)) ) )**(real(k,rk)+s)
     !
     uLR = uLR + Damp*parameters(1+k+Npot+Nstruc)/r**k
     !
    enddo
    !
    !
    ! at R=Re
    uLR0 = 0
    do k=1,M
     ! Douketis damping function
     Damp = ( 1.0_rk-exp( -b*rho*r0/real(k,rk)-c*(rho*r0)**2/sqrt(real(k,rk)) ) )**(real(k)+s)
     !
     uLR0 = uLR0 + Damp*parameters(1+k+Npot+Nstruc)/r0**k
     !
    enddo
    !
    betaN = 0
    do k=0,N
     betaN = betaN + parameters(k+Nstruc+1)*yq**k
    enddo
    !
    betainf = log(2.0_rk*de/uLR0)
    !
    beta = (1.0_rk-yp)*betaN+yp*betainf
    ! 
    y  = 1.0_rk-uLR/uLR0*exp(-beta*z)
    !



    !
    ma = 1.0_rk - parameters(2+M+Npot+Nstruc)/parameters(3+M+Npot+Nstruc)
    mb = 1.0_rk - parameters(4+M+Npot+Nstruc)/parameters(5+M+Npot+Nstruc)
    !
    
!    write(out,*) parameters(2+M+Npot+Nstruc), parameters(3+M+Npot+Nstruc)
!    write(out,*) parameters(4+M+Npot+Nstruc), parameters(5+M+Npot+Nstruc)


    p = nint(parameters(6+M+Npot+Nstruc))
    q = nint(parameters(7+M+Npot+Nstruc))
    !
    NUa = nint(parameters(8+M+Npot+Nstruc))
    !
    NUb = nint(parameters(9+M+Npot+Nstruc))
    !
    yp = (r**p-r0**p)/(r**p+r0**p)
    yq = (r**q-r0**q)/(r**q+r0**q)
    !
    u = 0
    do k=0,NUa
     u = u + ma*parameters(k+10+M+Npot+Nstruc)*yq**k
    enddo
    do k=0,NUb
     u = u + mb*parameters(k+12+NUa+M+Npot+Nstruc)*yq**k
    enddo
    !
    uinf = parameters(11+NUa+M+Npot+Nstruc)*ma + parameters(13+NUa+NUb+M+Npot+Nstruc)*mb



    f = de*y**2+v0 + ( (1.0_rk-yp)*u + uinf*yp )
    !
  end function poten_MLR_Douketis_Darby

 
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
     t = t + parameters(k+4)*z**k
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
  function poten_BOBna(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: r0,f,yq,t,tinf,ma,mb,q
    integer(ik)            :: k,NTa,NTb
    !
    r0 = parameters(1)
    !
    ma = parameters(2)/parameters(3)
    mb = parameters(4)/parameters(5)
!    write(out,*) r, ma, parameters(2), parameters(3), mb, parameters(4), parameters(5)
    !
    q = parameters(6)
    !
    NTa = nint(parameters(7))
    !
    NTb = nint(parameters(8))
    !
    yq = (r**q-r0**q)/(r**q+r0**q)
    !
    t = 0
    do k=0,NTa
     t = t + ma*parameters(k+9)*yq**k
    enddo
    do k=0,NTb
     t = t + mb*parameters(k+11+NTa)*yq**k
    enddo
    !
    tinf = parameters(10+NTa)*ma + parameters(12+NTa+NTb)*mb
    !
!    write(out,*) r, ( (1.0_rk-yq)*t + tinf*yq ),( (1.0_rk-yq)*t + tinf*yq )+1.0_rk
    f = ( (1.0_rk-yq)*t + tinf*yq )
    !
  end function poten_BOBna
  !
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
    !
    f = v0
    !
    if (N==1) return 
    !
    r0 = parameters(2)
    !
    y = (r-r0)
    do k=3,N 
      f = f + parameters(k)*y**(k-2)
    enddo
    !
  end function poten_polynom
  !
  ! The Chebyshev polynomial is evaluated at a point y=[x-(b + a)/2]/[(b - a)/2],
  ! and the result is returned as the function value = \sum_{k=0} a_k P_k  - a0/2
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
  ! Dipole moment Pade-2 by J. Goodisman, Dipole-moment function for diatomic molecules, 
  ! J. Chem. Phys. 38 (1963) 2597-2599. doi: 10.1063/1.1733557.
  ! based on Chebyshev polynomials
  !
  function poten_Pade_Goodisman_2(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: z,y,r0,f
    integer(ik)            :: N
    real(rk)               :: a(100)        ! temporal array of fixed size
    !
    N = size(parameters)
    !
    if (N+1>100) stop 'too many parameters requested in poten_Pade2'
    !
    r0 = parameters(1)
    !
    if (abs(r0)<small_) stop 'It is illegal to set re to zero'
    !
    z = r/r0
    y = (z-1.0_rk)/(z+1.0_rk)
    !
    a = 0
    !
    a(1) = -1.0_rk
    a(2) =  1.0_rk
    !
    a(3:N+1) = parameters(2:N)
    !
    f = poten_cheb(y,a(1:N+1))+a(3)*0.5_rk
    !
    f = z**3/(1.0_rk+z**7)*f
    !
  end function poten_Pade_Goodisman_2
  !
  !
  ! polynomial with exp decay
  !
  function dipole_polynom_exp(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: beta,z,f,r0
    integer(ik)            :: N,k
    !
    N = size(parameters)
    !
    r0 = parameters(1)
    beta = parameters(2)
    !
    if (abs(r0)<small_) stop 'It is illegal to set re to zero'
    !
    z = (r-r0)*exp(-beta*(r-r0)**2)
    !
    f =  0
    do k=0,N-3
      f = f + parameters(k+3)*z**k
    enddo
    !
  end function dipole_polynom_exp
  !
  function dipole_polynom_exp24(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: beta,z,f,r0,gamma,y,tinf,p
    integer(ik)            :: N,k
    !
    N = size(parameters)
    !
    r0 = parameters(1)
    beta = parameters(2)
    gamma = parameters(3)
    p = parameters(4)
    !
    if (abs(r0)<small_) stop 'It is illegal to set re to zero'
    !
    z = (r-r0)*exp(-beta*(r-r0)**2-gamma*(r-r0)**4)
    !
    f =  0
    do k=0,N-5-1
      f = f + parameters(k+5)*z**k
    enddo
    !
    y = (r**p-r0**p)/(r**p+r0**p)
    !
    tinf = parameters(N)
    !
    f = (1.0_rk-y)*f+y*tinf
    !
  end function dipole_polynom_exp24


  !
  ! polynomial with exp decay
  !
  function dipole_polynom_exp_damp(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: beta,z,f,r0,s,b,c,damp,tinf,p,y
    integer(ik)            :: N,k
    !
    N = size(parameters)
    !
    r0   = parameters(1)
    beta = parameters(2)
    !
    p = parameters(3)
    s = parameters(4)
    b = parameters(5)
    c = parameters(6)
    !
    tinf = parameters(N)
    !
    if (abs(r0)<small_) stop 'It is illegal to set re to zero'
    !
    z = (r-r0)*exp(-beta*(r-r0)**2)
    !
    f =  0
    do k=7,N-1
      f = f + parameters(k)*z**(k-7)
    enddo
    !
    ! Douketis damping function    
    Damp = ( 1.0_rk-exp( -b*r-c*r**2 ) )**s
    !
    y = (r**p-r0**p)/(r**p+r0**p)
    !
    f = ((1.0_rk-y)*f+y*tinf)*damp
    !
  end function dipole_polynom_exp_damp



  !
  function SO_arctan(r,parameters) result(fun)
     implicit none
     real(rk),intent(in) :: r,parameters(:)
     integer(ik)  :: m
     real(rk)            :: fun,r0,SO0,SOinf,k,q,f
     !
     m=size(parameters)
     !
     if(m/=5) stop 'SO_arctan: illegal number of parameters parameters (/=5)'
     !
     r0     = parameters(1)
     SO0    = parameters(2)
     SOinf  = parameters(3)
     k      = parameters(4)
     f      = parameters(5)
     !
     q = (r-r0)*k
     !
     fun = f*(0.5_rk*(SO0+SOinf)+(1.0_rk/pi)*(SO0-SOinf)*atan(q))
     !
   end function SO_arctan
   !
  !
  function dipole_doubleexp(r,parameters) result(fun)
    implicit none
    real(rk),intent(in) :: r,parameters(:)
    integer(ik)  :: m
    real(rk)            :: fun,c1,c2,a1,a2,m1,m2
    !
    m=size(parameters)
    !
    if(m/=6) stop 'dipole_doubleexp: illegal number of parameters parameters (/=6)'
    !
    c1=parameters(1)
    a1=parameters(2)
    m1=parameters(3)
    c2=parameters(4)
    a2=parameters(5)
    m2=parameters(6)
    !
    fun = c1*exp(-a1*r**m1)+c2*exp(-a2*r**m2)
    !
  end function dipole_doubleexp
  !
  function dipole_avoidedcrossing_diag_mu(r,parameters) result(fun)
    implicit none
    real(rk),intent(in) :: r,parameters(:)
    integer(ik)  :: m
    real(rk)            :: fun,j0,m0,d0,k0,t0,l1,l2,lambda,inter,muionic
    !
    m=size(parameters)
    !
    if(m/=7) stop 'dipole_avoidedcrossing_diag_mu: illegal number of parameters parameters (/=7)'
    !
    j0=parameters(1)
    m0=parameters(2)
    d0=parameters(3)
    k0=parameters(4)
    t0=parameters(5)
    l1=parameters(6)
    l2=parameters(7)
    !
    lambda=m0*(r-d0)*(r-k0)
    inter=(sqrt(4*j0**2+lambda**2)+lambda)**2
    muionic=t0*r+l1/r+l2/(r**3)
    !    
    fun = inter*muionic/(inter+4.0_rk*j0**2)
    !
  end function dipole_avoidedcrossing_diag_mu


  function poten_two_coupled_EMOs(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: f1,f2,a,e(2),f,discr
    integer(ik)            :: nparams1,nparams2,nparams3,icomponent
    !
    nparams1 = parameters(8)+9
    nparams2 = parameters(nparams1+8)+9
    nparams3 = size(parameters)-(nparams1+nparams2)-1 ! last parameter is the adiabatic component 
    icomponent = parameters(nparams1+nparams2+nparams3+1)
    !
    f1 = poten_EMO(r,parameters(1:nparams1))
    f2 = poten_EMO(r,parameters(nparams1+1:nparams1+nparams2))
    a  = poten_cosh_polynom(r,parameters(nparams1+nparams2+1:nparams1+nparams2+nparams3))
    !
    discr = f1**2-2.0_rk*f1*f2+f2**2+4.0_rk*a**2
    !
    if (discr<-small_) then
      write(out,"('poten_two_coupled_EMOs: discriminant is negative')")
      stop 'poten_two_coupled_EMOs: discriminant is negative'
    endif
    !
    e(1)=0.5_rk*(f1+f2)-0.5_rk*sqrt(discr) 
    e(2)=0.5_rk*(f1+f2)+0.5_rk*sqrt(discr) 
    !
    f = e(icomponent)
    !
  end function poten_two_coupled_EMOs


  function poten_two_coupled_BOBs(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: f1,f2,f,f_switch,r_s,a_s
    integer(ik)            :: nparams1,nparams2,nparams3,icomponent,n
    !
    nparams1 = parameters(4)+6
    nparams2 = parameters(nparams1+4)+6
    nparams3 = 2
    icomponent = parameters(nparams1+nparams2+nparams3+1)
    !
    N = nparams1+nparams2
    !
    f1 = poten_BOBLeRoy(r,parameters(1:nparams1))
    f2 = poten_BOBLeRoy(r,parameters(nparams1+1:nparams1+nparams2))
    !
    r_s = parameters(N+1)
    a_s = parameters(N+2)
    !
    f_switch = ( 1.0_rk+tanh(a_s*(r-r_s)) )*0.5_rk
    !
    f = 0 
    !
    if (icomponent==1) then 
      f = f_switch*f2+f1*(1.0_rk-f_switch)
    else
      f = f_switch*f1+f2*(1.0_rk-f_switch)
    endif
    !
  end function poten_two_coupled_BOBs
  !
  ! Morse/Long-Range, see Le Roy manuals
  !
  function poten_repulsive(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: v0,f,uLR
    integer(ik)            :: k,N
    !
    N = parameters(1)
    v0 = parameters(2)
    !
    if (N/=size(parameters)-1) then
      write(out,"('Repulsive: Npar inconsistent with total number of parameters:',2i9)") N,size(parameters)+1
      stop 'Repulsive: illegal number of parameters'
    endif 
    !
    ! long-range part
    !
    uLR = V0
    do k=2,N
     uLR = uLR + parameters(k+1)/r**(k-1)
    enddo
    !
    f = uLR
    !
  end function poten_repulsive


  function poten_two_coupled_EMO_repulsive(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: f1,f2,a,e(2),f,discr
    integer(ik)            :: nparams1,nparams2,nparams3,icomponent
    !
    nparams1 = parameters(8)+9
    nparams2 = parameters(nparams1+1)+1
    nparams3 = size(parameters)-(nparams1+nparams2)-1 ! last parameter is the adiabatic component 
    icomponent = parameters(nparams1+nparams2+nparams3+1)
    !
    f1 = poten_EMO(r,parameters(1:nparams1))
    f2 = poten_repulsive(r,parameters(nparams1+1:nparams1+nparams2))
    a  = poten_cosh_polynom(r,parameters(nparams1+nparams2+1:nparams1+nparams2+nparams3))
    !
    discr = f1**2-2.0_rk*f1*f2+f2**2+4*a**2
    !
    if (discr<-small_) then
      write(out,"('poten_two_coupled_EMOs: discriminant is negative')")
      stop 'poten_two_coupled_EMOs: discriminant is negative'
    endif
    !
    e(1)=0.5_rk*(f1+f2)-0.5_rk*sqrt(discr) 
    e(2)=0.5_rk*(f1+f2)+0.5_rk*sqrt(discr) 
    !
    f = e(icomponent)
    !
  end function poten_two_coupled_EMO_repulsive

  !
  ! A lorentzian function for the couplings between diabatic curves 
  !
  function poten_lorentzian_polynom(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y0,r0,w,a,z,f0,f
    integer(ik)            :: k,N
    !
    N = size(parameters)
    !
    y0 = parameters(1)
    r0 = parameters(2)
    w = parameters(3)
    a = parameters(4)
    !
    z = (r-r0)
    !
    f0 = a
    do k=5,N
     f0 = f0 + parameters(k)*z**(k-4)
    enddo
    !
    f = y0+2.0_rk*f0/pi*( w/( 4.0_rk*(r-r0)**2+w**2 ) )
    !
  end function poten_lorentzian_polynom

                                                                
  function polynomial_dimensionless(r,parameters) result(f)
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
    y = (r-r0)/r0
    !
    f = parameters(N)
    do k=N-1,3,-1
      !f = f + parameters(k)*y**(k-2)
      f = f*y + parameters(k)
    enddo
    !
    f = f + v0
    !
  end function polynomial_dimensionless

  function potential_stolyarov_CO_X_UBO(r,parameters) result(f)
  
   !     CO X BO potential. R in angstrom, CO_X_U in cm^-1
   !     JQSRT 217 (2018) 262–273
   !
   implicit none
   real(rk),intent(in)    :: r             ! geometry (Ang)
   real(rk),intent(in)    :: parameters(:) ! potential parameters
   integer(ik),parameter ::          na = 14
   real(rk) :: T, K, E0, d, d1, d2, c1, c2, C5, C6, C8, f , a(na)
   real(rk) :: z, y0, yinf, t1, t2, t3
   integer(ik)  :: i 
     !
     data a/  -1.72162379120d+00,&
              -1.57663684255d+00,&
               1.62435769080d+00,&
               6.48364603740d-01,&
               2.95603114170d-01,&
               6.11299378400d-02,&
              -1.19431556940d-01,&
              -9.00931765500d-02,&
              -6.48242550300d-02,&
              -3.62863354500d-02,&
              -1.30105004800d-02,&
              -6.07793671000d-03,&
              -1.07070127000d-03,&
              -4.00250370000d-04/
     !
     T    =  parameters(1)
     K    =  0.55747667156e+07_rk
     E0   = -0.38846898e+08_rk
     !    
     d    =  6.441_rk
     d1   = 16.829_rk
     d2   =  0.5363_rk
     c1   =  0.5704_rk
     c2   =  1.2668_rk
     !
     C5   =  0.102093e+06_rk
     C6   =  0.655440e+05_rk
     C8   =  0.496639e+06_rk
     y0   = (K/r + K*d + E0 + d*(0.5_rk*K*d + e0)*r)*exp(-d*r)
     !    
     yinf = (DD(r,5,d1,d2)*C5 + DD(r,6,d1,d2)*C6/r + DD(r,8,d1,d2)*C8/r**3)/r**5
     !       
     z    = dtanh(c1*r - c2/r)
     !       
     t1 = 0._rk 
     t2 = 0._rk
     do i = na, 2, -1
         t3 = t1
         t1 = 2.0_rk*z*t1 - t2 + a(i)
         t2 = t3
     enddo
     !
     t3 = z*t1 - t2 + a(1)
     !     
     f  = (y0 + yinf)*t3 + T
     !
     contains  
      !
      function DD(r, n, d1, d2) result(f)
      real(rk) :: r, d1, d2, f
      integer(ik) ::   n
      !
      f = (1.0_rk - dexp(-d1*r/n - d2*r*r/dsqrt(dble(n))))**(n + 2)
      !
      end function DD
      !
  end function potential_stolyarov_CO_X_UBO
  !
  !
end module functions
