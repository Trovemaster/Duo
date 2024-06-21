module functions
  !
  use accuracy
  use timer
  !
  implicit none
  !
  public fanalytic_fieldT
  !
  ! Different 1D functions (potential, spin-orbit etc)
  !
  integer(ik),parameter :: verbose=5
  !
  !logical,save :: defined_complex_field = .false.
  integer(ik),save  :: N1,N2,N3
  type array_of_funcitonsT
    procedure(fanalytic_fieldT) ,pointer ,nopass :: f =>null()
  end type array_of_funcitonsT
  !
  procedure (fanalytic_fieldT),pointer :: function_V1 => null()
  procedure (fanalytic_fieldT),pointer :: function_V2 => null()
  procedure (fanalytic_fieldT),pointer :: function_VD => null()
  procedure (fanalytic_fieldT),pointer :: function_beta => null()
  !
  type(array_of_funcitonsT),pointer :: function_multi(:)
  !
  integer(ik) :: N_multi_subfunctions
  integer(ik),allocatable :: nsize_multi(:)
  !
  abstract interface
    !
    function fanalytic_fieldT(r,parameters)
      use accuracy
      !
      real(rk)                      :: fanalytic_fieldT !NB: NAG Fortran 6.0 doesn't like in-line declaration
      real(rk),intent(in)           :: r                 ! geometry (Ang)
      real(rk),intent(in)           :: parameters(:)     ! potential parameters
      !
    end function fanalytic_fieldT
    !
  end interface
  !
  abstract interface
    !
    function fanalytic_complex_fieldT(r,parameters,sub_type,Nsub_terms)
      use accuracy
      !
      real(rk)                      :: fanalytic_complex_fieldT !NB: NAG Fortran 6.0 doesn't like in-line declaration
      real(rk),intent(in)           :: r                 ! geometry (Ang)
      real(rk),intent(in)           :: parameters(:)     ! potential parameters
      character(len=cl),intent(in)  :: sub_type(:)
      integer(ik),intent(in)        :: Nsub_terms(:) ! Number of terms oin sub-functions, used for coupled representations 
      !
    end function fanalytic_complex_fieldT
    !
  end interface
  !
  contains
  !
  !
  subroutine define_fanalytic_field(ftype,fanalytic_field)
    !
    character(len=cl),intent(in)           :: ftype
                                           ! NB: NAG Fortran 6.0 doesn't like intent(in) and initial nullification
    procedure(fanalytic_fieldT),pointer :: fanalytic_field !=> null()
    !
    select case(ftype)
      !
    case("MORSE")
      !
      fanalytic_field => poten_morse
      !
    case("MORSE_DAMP")
      !
      fanalytic_field => poten_morse_damp
      !
    case("MODIFIED-MORSE","MODIFIED_MORSE","MMORSE")
      !
      fanalytic_field => poten_morse_modified
      !
    case("EMO") ! "Expanded MorseOscillator"
      !
      fanalytic_field => poten_EMO
      !
    case("EMO-BOB") ! "Expanded MorseOscillator with BOB correction"
      !
      fanalytic_field => poten_EMO_BOB
      !
    case("MLR") ! "Morse/Long-Range"
      !
      fanalytic_field => poten_MLR
      !
    case("MLR_DS") ! "Morse/Long-Range with Douketis-damping"
      !
      fanalytic_field => poten_MLR_Douketis
      !
    case("DELR") ! "Double-exponential-long-range"
      !
      fanalytic_field => poten_DELR
      !
    case("MLR_DS_DARBY") ! "Morse/Long-Range with Douketis-damping"
      !
      fanalytic_field => poten_MLR_Douketis_Darby
    case("MLR3") ! "MLR3 with Douketis damping (see Coxon & Hajigeorgiou 2010)
      !
      fanalytic_field => poten_MLR3
    case("MARQUARDT") ! "Marquardt"
      !
      fanalytic_field => poten_Marquardt
      !
    case("COSH-POLY") ! "Diabatic coupling as a polynom/cosh"
      !
      fanalytic_field => poten_cosh_polynom
      !
    case("BOBLEROY","BOB","SURKUS") ! "BOB expansion"
      !
      fanalytic_field => poten_BOBLeRoy
      !
    case("BOBNA") ! "BOB-NA expansion"
      !
      fanalytic_field => poten_BOBna
      !
    case("BOBLEROY_DAMP","BOB_DAMP","SURKUS_DAMP") ! "BOB expansion with damping to zero at r=0"
      !
      fanalytic_field => poten_BOBLeRoy_damp
      !
    case("DUNHAM")
      !
      fanalytic_field => poten_dunham
      !
    case("SPF")
      !
      fanalytic_field => poten_spf
      !
    case("SPF_1")
      !
      fanalytic_field => poten_spf_1
      !
    case("SPF_H2")
      !
      fanalytic_field => poten_spf_h2
      !
    case("HH","HULBERT-HIRSCHFELDER")
      !
      fanalytic_field => poten_Hulbert_Hirschfelder
      !
    case("CHEBYSHEV")
      !
      fanalytic_field => poten_cheb
      !
    case("IRREG_CHEBYSHEV")
      !
      fanalytic_field => irreg_chebyshev_DMC
      !
    case("POLYNOM", "POLYNOMIAL")
      !
      fanalytic_field => poten_polynom
      !
    case("EXPONENTIAL")
      !
      fanalytic_field => poten_exponential
      !
    case("LAURA_SO","ATAN_SO")
      !
      fanalytic_field => SO_arctan
      !
    case("POTEN_FERMI_T")
      !
      fanalytic_field => poten_fermi_t
      !
    case("M-S")
      !
      fanalytic_field => poten_Murrell_Sorbie
      !
    case("PADE_GOODISMAN2","PADE2")
      !
      fanalytic_field => poten_Pade_Goodisman_2
      !
    case("POLYNOM_DECAY")
      !
      fanalytic_field => dipole_polynom_exp
      !
    case("POLYNOM_DECAY_DAMP")
      !
      fanalytic_field => dipole_polynom_exp_damp
      !
    case("POLYNOM_DECAY_24")
      !
      fanalytic_field => dipole_polynom_exp24
      !
    case("DOUBLEEXP2")
      !
      fanalytic_field => dipole_doubleexp
      !
    case("SIGMOID")
      !
      fanalytic_field => poten_sigmoid
      !
    case("DIABATIC_MU_DIAG","AVOIDEDCROSSING_DIAG_MU")
      !
      fanalytic_field => dipole_avoidedcrossing_diag_mu
      !
    case("TWO_COUPLED_EMOS")
      !
      fanalytic_field => poten_two_coupled_EMOs
      !
    case("TWO_COUPLED_EMOS_LORENTZ")
      !
      fanalytic_field => poten_two_coupled_EMOs_Lorentz
      !
    case("DIABATIC_LORENTZ_TWO_EMOS")
      !
      fanalytic_field => poten_diabatic_Lorentz_TWO_EMOS
      !
    case("TWO_COUPLED_EMOS_SQRTLORENTZ")
      !
      fanalytic_field => poten_two_coupled_EMOs_SqrtLorentz
      !
    case("TWO_COUPLED_BOBS")
      !
      fanalytic_field => poten_two_coupled_BOBs
      !
    case("COUPLED_EMO_REPULSIVE")
      !
      fanalytic_field => poten_two_coupled_EMO_repulsive
      !
    case("COUPLED_EMOS_WITH_EMO")
      !
      fanalytic_field => poten_two_coupled_EMOs_with_EMO
      !
    case("REPULSIVE")
      !
      fanalytic_field => poten_repulsive
      !
    case("LORENTZ","LORENTZIAN")
      !
      fanalytic_field => poten_lorentzian_polynom
      !
    case("LORENTZ-SURKUS","LORENTZIAN-SURKUS")
      !
      fanalytic_field => poten_lorentzian_surkus_polynom
      !
    case("SQRT(LORENTZ)","SQRT(LORENTZIAN)")
      !
      fanalytic_field => poten_sqrt_lorentzian_polynom
      !
    case("BETA_LORENTZ","BETA_LORENTZIAN")
      !
      fanalytic_field => beta_diabatic_lorentzian
      !
    case("BETA_LAPLACIAN","BETA_LAPLACE")
      !
      fanalytic_field => beta_diabatic_laplacian
      !
    case("BETA_LORENTZ_LAPLACE")
      !
      fanalytic_field => beta_diabatic_lorentzian_laplacian
      !
    case("POLYNOM_DIMENSIONLESS","POLYNOMIAL_DIMENSIONLESS")
      !
      fanalytic_field => polynomial_dimensionless
      !
    case("CO_X_UBO")
      !
      fanalytic_field => potential_stolyarov_CO_X_UBO
      !
    case("EHH") !  Extended Hulburt-Hirschfelde
      !
      fanalytic_field => poten_EHH
      !
    case("MEDVEDEV_SING2","SING2") ! Irregular DMF by Medvedev Opt. spectrosc. 130, 1334 (2022)
      !
      fanalytic_field => dipole_medvedev_sing
      !
    case("NONE")
      !
      write(out,'(//"Analytic: Some fields are not properly defined and produce function type ",a)') trim(ftype)
      stop "Analytic: Unknown function type "
      !
    case("COUPLED-PEC","GRID")
      !
      !fanalytic_complex_field => polynomial_coupled_PECs!
      !
      fanalytic_field => function_dummy
      !
    case("COUPLED-PEC-BETA","COUPLED-DIABATIC","COUPLED-TRANSIT-BETA")
      !
      fanalytic_field => function_dummy
      !
    case default
      !
      write(out,'(//"Analytic: Unknown function type ",a)') trim(ftype)
      stop "Analytic: Unknown field type "
      !
    end select
    !
  end subroutine define_fanalytic_field
  !
  !
  subroutine define_complex_analytic_field(ftype,fanalytic_field,sub_type,Nsub_terms)
    !
    character(len=cl),intent(in)           :: ftype
                                           ! NB: NAG Fortran 6.0 doesn't like intent(in) and initial nullification
    procedure(fanalytic_fieldT),pointer :: fanalytic_field !=> null()
    character(len=cl),intent(in)  :: sub_type(:)
    integer(ik),intent(in)        :: Nsub_terms(:) ! Number of terms oin sub-functions, used for coupled representations 
    integer(ik) :: i,alloc
    !
    N_multi_subfunctions = size(Nsub_terms)
    !
    select case(ftype)
      !
    case("COUPLED-PEC")
      !
      if (N_multi_subfunctions==3) then 
        !
        call define_fanalytic_field(sub_type(1),function_V1)
        call define_fanalytic_field(sub_type(2),function_V2)
        call define_fanalytic_field(sub_type(3),function_VD)
        !
        N1 =Nsub_terms(1)
        N2 =Nsub_terms(2)
        N3 =Nsub_terms(3)
        !
        fanalytic_field => polynomial_coupled_PECs
        !
      else
        !
        if (associated(function_multi)) deallocate(function_multi)
        if (allocated(nsize_multi)) deallocate(nsize_multi)
        !
        allocate(function_multi(N_multi_subfunctions),stat=alloc)
        if (alloc/=0) then
          write (out,"(' Error define_complex_analytic_field ',i0,' initializing function_multi')") alloc
          stop 'Error define_complex_analytic_field initializing function_multi - alloc'
        end if
        !
        allocate(nsize_multi(N_multi_subfunctions),stat=alloc)
        if (alloc/=0) stop 'Error nsize_multi - alloc'
        ! 
        do i = 1,N_multi_subfunctions
          !
          call define_fanalytic_field(sub_type(i),function_multi(i)%f)
          nsize_multi(i) = Nsub_terms(i)
          !
        enddo
        !
        fanalytic_field => polynomial_multi_coupled_PECs
        !
      endif
      !
    case("COUPLED-PEC-BETA")
      !
      if (N_multi_subfunctions/=3) then 
        write(out,"('Illegal number of sub-funcitons for COUPLED-PEC-BETA',i3,' only 3 has been immplemented')") &
              N_multi_subfunctions
        stop 'Illegal number of sub-funcitons for COUPLED-PEC-BETA/=3'
      endif
      !
      call define_fanalytic_field(sub_type(1),function_V1)
      call define_fanalytic_field(sub_type(2),function_V2)
      call define_fanalytic_field(sub_type(3),function_beta)
      !
      N1 =Nsub_terms(1)
      N2 =Nsub_terms(2)
      N3 =Nsub_terms(3)
      !
      fanalytic_field => polynomial_coupled_PECs_beta
      !
    case("COUPLED-TRANSIT-BETA")
      !
      call define_fanalytic_field(sub_type(1),function_V1)
      call define_fanalytic_field(sub_type(2),function_V2)
      call define_fanalytic_field(sub_type(3),function_beta)
      !
      N1 =Nsub_terms(1)
      N2 =Nsub_terms(2)
      N3 =Nsub_terms(3)
      !
      fanalytic_field => coupled_transition_curves_beta
      !
    case default
      !
      write(out,'(//"Complex Analytic: Unknown function type ",a)') trim(ftype)
      stop "Complex Analytic: Unknown field type "
      !
    end select
    !
  end subroutine define_complex_analytic_field
  !  
  !
  subroutine define_sub_terms_complex_analytic_field(Nsub_terms)
    !
    integer(ik),intent(in)        :: Nsub_terms(:) ! Number of terms oin sub-functions, used for coupled representations 
    integer(ik) :: alloc,i
    !
    N_multi_subfunctions = size(Nsub_terms)
    !
    if (N_multi_subfunctions==3) then 
      !
      N1 =Nsub_terms(1)
      N2 =Nsub_terms(2)
      N3 =Nsub_terms(3)
      !
    else
      !
      if (allocated(nsize_multi)) deallocate(nsize_multi)
      !
      allocate(nsize_multi(N_multi_subfunctions),stat=alloc)
      if (alloc/=0) stop 'Error nsize_multi - alloc'
      ! 
      do i = 1,N_multi_subfunctions
        !
        nsize_multi(i) = Nsub_terms(i)
        !
      enddo
      !
    endif    
    !
  end subroutine define_sub_terms_complex_analytic_field
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
  function poten_sigmoid(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,r0,f,z,v0,D0,phi
    integer(ik)            :: k,N,Nbeta,p,N0
    !
    N = size(parameters)
    !
    Nbeta = N - 5
    !
    v0 = parameters(1)
    r0 = parameters(2)
    ! Note that the De is relative the absolute minimum of the ground state
    D0 = parameters(3)-v0
    p = nint(parameters(4))
    !
    N0 = 5
    !
    z = (r**p-r0**p)/(r**p+r0**p)
    !
    phi = 0
    do k=0,Nbeta
     phi = phi + parameters(k+N0)*z**k
    enddo
    !
    y  = 1.0_rk+exp(-phi*(r-r0))
    !
    f = d0/y+v0
    !
  end function poten_sigmoid


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
  !
  function poten_EMO_BOB(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0,de,f,rref,z,phi, ma, mb, yp, yq, u, uinf
    integer(ik)            :: k,N,p, NUa, NUb, q
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
    NUa = nint(parameters(16+N))
    NUb = nint(parameters(17+N))
    !
    if (size(parameters)/=8+max(parameters(7),parameters(8))+1+12+NUa+NUb) then 
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
    ma = 1.0_rk - parameters(10+N)/parameters(11+N)
    mb = 1.0_rk - parameters(12+N)/parameters(13+N)
    !
    p = nint(parameters(14+N))
    q = nint(parameters(15+N))
    !
    !
    yp = (r**p-r0**p)/(r**p+r0**p)
    yq = (r**q-r0**q)/(r**q+r0**q)
    !
    u = 0
    do k=0,NUa
     u = u + ma*parameters(k+18+N)*yq**k
    enddo
    do k=0,NUb
     u = u + mb*parameters(k+20+NUa+N)*yq**k
    enddo
    !
    uinf = parameters(19+NUa+N)*ma + parameters(21+NUa+NUb+N)*mb
    !
    f = de*y**2+v0+ ( (1.0_rk-yp)*u + uinf*yp )
    !
  end function poten_EMO_BOB


  
  
  
  !
  ! Morse/Long-Range, see Le Roy manuals
  !
  function poten_MLR(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0,de,f,rref,z,phi,phiinf,uLR,uLR0,phi0
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
    if (r<rref) then
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
    Npot = max(parameters(6),parameters(7))
    !
    ! number of long range parameters
    !
    M = Ntot-Npot-Nstruc-1
    !
    !if (size(parameters)/=8+max(parameters(7),parameters(8))+1) then
    !  write(out,"('poten_EMO: Illegal number of parameters in EMO, check NS and NL, must be max(NS,NL)+9')")
    !  stop 'poten_EMO: Illegal number of parameters, check NS and NL'
    !endif
    !
    z = (r**p-rref**p)/(r**p+rref**p)
    !
    phi0 = 0
    do k=0,N
     phi0 = phi0 + parameters(k+Nstruc+1)*z**k
    enddo
    !
    !phiinf = parameters(Npot+Nstruc+2)
    !
    ! double check uLR parameters
    !
    uLR = sum(parameters(Npot+Nstruc+2:M+Npot+Nstruc+1)**2)
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
     uLR = uLR + parameters(k+Npot+Nstruc+1)/r**k
    enddo
    !
    ! at R=Re
    uLR0 = 0
    do k=1,M
     uLR0 = uLR0 + parameters(k+Npot+Nstruc+1)/r0**k
    enddo
    !
    phiinf = log(2.0_rk*de/uLR0)
    !phiinf = -8.8481839714D-01
    !
    phi = (1.0_rk-z)*phi0+z*phiinf
    !
    z = (r**p-r0**p)/(r**p+r0**p)
    !
    y  = 1.0_rk-uLR/uLR0*exp(-phi*z)
    !y  = 1.0_rk-uLR/uLR0*exp(-phi*(r-r0))
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
    !
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


  ! Double-exponential-long-range, see Le Roy http://dx.doi.org/10.1063/1.1607313
  !
  function poten_DELR(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: v0,r0,de,f,rref,beta_,yp,uLR,uLR0,uLR1,B,A,D_LR
    real(rk)               :: beta(0:100),C(1:100)
    integer(ik)            :: k,p,Nstruc,Ntot,NLR,NL,NS
    !
    v0 = parameters(1)
    r0 = parameters(2)
    ! Note that the De is relative the absolute minimum of the ground state
    De   = parameters(3)-V0
    D_LR = parameters(4)-V0
    !
    rref = parameters(5)
    !
    if (rref<=0.0_rk) rref = r0
    !
    p = nint(parameters(6))
    !
    NS = parameters(7)
    NL = parameters(8)
    !
    Nstruc = 8
    !
    if (NL>100) stop 'poten_DELR_Douketis: illegally large NL'
    !
    beta = 0 
    beta(0:NL) = parameters(Nstruc+1:Nstruc+1+NL+1)
    !
    ! total number of parameters
    !
    Ntot = size(parameters(:),dim=1)
    !
    NLR = Ntot-(Nstruc+1+NL)
    !
    if (NLR>100) stop 'poten_DELR_Douketis: illegally large LR expansion'
    !
    C = 0 
    !
    C(1:NLR) = parameters(Nstruc+1+NL+1:Nstruc+1+NL+NLR)
    !
    yp = (r**p-rref**p)/(r**p+rref**p)
    !
    ! double check uLR parameters
    !
    uLR = sum(C(1:NLR)**2)
    !
    if (uLR<small_) then
      write(out,"('poten_DELR_Douketis: At least one uLR should be non-zero')")
      stop 'poten_DELR_Douketis: At least one uLR should be non-zer'
    endif
    !
    !s = -1.0_rk
    !b = 3.97_rk
    !c0 = 0.39_rk
    !
    ! For the Damping part
    ! the values of s, b,c are as suggested by LeRoy 2011 (MLR paper)
    !
    !s = -1.0_rk
    !b = 3.3_rk
    !c = 0.423_rk
    !
    !s = -2.0_rk
    !b = 2.50_rk
    !c = 0.468_rk
    !
    ! long-range part
    !
    uLR  = D_LR
    uLR0 = D_LR
    uLR0 = 0
    do k=1,NLR
     !
     ! Douketis damping function
     !Damp = ( 1.0_rk-exp( -b*rho*r/real(k,rk)-c*(rho*r)**2/sqrt(real(k,rk)) ) )**(real(k,rk)+s)
     !
     uLR  = uLR  + C(k)/r**k
     uLR0 = uLR0 + C(k)/r0**k
     uLR1 = uLR1 + real(-k,rk)*C(k)/r0**(k+1)
    enddo
    !
    beta_ = 0
    do k=0,NL
     beta_ = beta_ + beta(k)*yp**k
    enddo
    !
    A = De + uLR0+uLR1/beta(0)
    B = 2.0_rk*(De + uLR0)+uLR1/beta(0)
    !
    f = v0+A*exp(-2.0_rk*beta_*(r-r0))-B*exp(-beta_*(r-r0))+De+uLR
    !
  end function poten_DELR





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
    !
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
    ma = 1.0_rk - parameters(2+M+Npot+Nstruc)/parameters(3+M+Npot+Nstruc)
    mb = 1.0_rk - parameters(4+M+Npot+Nstruc)/parameters(5+M+Npot+Nstruc)
    !
!    write(out,*) parameters(2+M+Npot+Nstruc), parameters(3+M+Npot+Nstruc)
!    write(out,*) parameters(4+M+Npot+Nstruc), parameters(5+M+Npot+Nstruc)
    !
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
    !
    f = de*y**2+v0 + ( (1.0_rk-yp)*u + uinf*yp )
    !
  end function poten_MLR_Douketis_Darby

  function poten_MLR3(r, parameters) result(f)
    !
    real(rk), intent(in)  :: r ! geometry (Angstroms)
    real(rk), intent(in)  :: parameters(:) ! vector of potential parameters
    real(rk)              :: v0, r0, De, rRef, a, s, rho, b, c
    real(rk), allocatable :: coefs(:), phis(:), pwrs(:)
    integer(ik)           :: nPwrs, nPhis, n, i, j, p, m, q
    real(rk)              :: f, uLR, uLRe, Dn, Dne, phiInf, ym, yq, ypa, sumPhis

    v0 = parameters(1) ! potential minimum
    r0 = parameters(2) ! equilibrium bond
    De = parameters(3) - v0 ! well depth relative to minimum
    ! MLR3 parameters
    rRef = parameters(4) ! reference r for y_m y_q variables in phi function
    p    = nint(parameters(5)) ! power of r for y_pa variable in MLR exponent
    m    = nint(parameters(6)) ! power of r for y_m variable in phi function
    q    = nint(parameters(7)) ! power of r for y_q variable in phi function
    a    = parameters(8) ! denominator r_e factor in y_pa variable
    ! Douketis damping parameters
    s      = parameters(9) ! specifies member of Douketis functions
    rho    = parameters(10) ! related to ionisation potential of atoms
    b      = parameters(11) ! b,c constants for Douketis damping
    c      = parameters(12)
    ! number & order of inverse powers can vary (in theory) so assign dynamically
    nPwrs = nint(parameters(13)) ! number of inverse power terms
    allocate ( coefs(nPwrs) )
    allocate ( pwrs(nPwrs) )
    pwrs  = nint(parameters(14:13+nPwrs)) ! the order of each power term
    coefs = parameters(14+nPwrs:13+2*nPwrs) ! the coefficients for each power term
    ! so can the number of phi terms
    nPhis = nint(parameters(14+2*nPwrs)) ! number of phi coefficients
    allocate ( phis(nPhis) )
    phis  = parameters(15+2*nPwrs:) ! the phi coefficients for MLR3 (arb. number)

    ! check the no. of power term and phi coefficients given matches the number specified
    if ((size(phis) .NE. nPhis) .OR. (size(pwrs) .NE. nPwrs)) then
      write(out, "('poten_MLR3: Number of power terms or phi coefficients&
        & specified does not match nPwrs or nPhis')")
      stop 'poten_MLR3: Number of power terms or phi coefficients&
        & specified does not match nPwrs or nPhis'
    endif
    ! check that some phis and power terms have been given
    if ((size(phis) .EQ. 0) .OR. (size(pwrs) .EQ. 0)) then
      write(out, "('poten_MLR3: Requires more than zero phi coefficients&
        & and more than zero power terms')")
      stop 'poten_MLR3: Requires more than zero phi coefficients&
        & and more than zero power terms'
    endif

    ! if rRef not given, equals equilibrium bond length
    if (rRef .LE. 0.0_rk) then
      write(out, "('poten_MLR3: rRef equal to zero or less, using rRef = r0')")
      rRef = r0
    endif

    ! calculate U_LR (Long-range potential) by summing over power terms
    uLR = 0.0_rk
    uLRe = 0.0_rk
    do j = 1, nPwrs
      n = pwrs(j)
      Dn  = (1.0_rk - exp( -b*rho*r/real(n,rk) - c*(rho*r)**2/real(n,rk)**.5_rk ))**(n+s)
      uLR = uLR + (Dn * coefs(j)/r**n)
      ! also uLR at r_e (this is the same for every r, but I can't be bothered
      ! to think of an implementation that means it only needs computing once
      Dne = (1.0_rk - exp( -b*rho*r0/real(n,rk) - c*(rho*r0)**2/real(n,rk)**.5_rk ))**(n+s)
      uLRe= uLRe + (Dne * coefs(j)/r0**n)
    enddo

    ! calculate phi_MLR3(r) - sum of
    phiInf = log(2.0_rk*De/uLRe)
    ym     = (r**m - rRef**m)/(r**m + rRef**m)
    yq     = (r**q - rRef**q)/(r**q + rRef**q)
    sumPhis = 0.0_rk
    do i = 1, nPhis
      sumPhis = sumPhis + phis(i) * yq**(i-1)
    enddo
    sumPhis = sumPhis*(1.0_rk - ym) + phiInf*ym

    ! put parts together to give total potential
    ypa = (r**p - r0**p)/(r**p + a*r0**p)
    f = De*( 1.0_rk - (uLR/uLRe) * exp(-sumPhis*ypa) )**2
  end function poten_MLR3
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
  function irreg_chebyshev_DMC(x,c) result(f)        ! based on eq.(3) of https://doi.org/10.1016/j.jqsrt.2022.108255
  real(rk),intent(in) :: x,c(:)                      ! x = bond length, c = functional parameters  
  integer(ik) ::  m,j                                ! Chebyshev expansion order, iterable 
  real(rk)    ::  d,dd,sv,z,z2,cheb_expansion,chi,f  ! d,dd,sv = appear in Clenshaw's recurrence, z = damped polynomial coordinate,
                                                     ! cheb_expansion = sum_k=0^size
    
    ! c_1,c_2,...,c_7 are the chebyshev expansion coefficients, c_9,...,c_(13) are additional parameters for chi below
    ! for chi and c_8 enters the damped polynomial coordinate z.
    m = size(c)-6 ! expansion order
    !
    ! check if number of parameters is correct
    if (size(c).gt.13) stop 'Too many parameters, should have 7 expansion coeficients and 6 fitting parameters.'
    if (size(c).lt.13) stop 'Too few parameters, should have 7 expansion coeficients and 6 fitting parameters.'
    !
    ! Change of variable, the damped polynomial coordinate, z, which maps the r E [0,infty] -> z E [-1,+1] interval
    z = 1-2*exp(-c(8)*x) 
    !
    ! check to see if the c0 is > 0
    if (c(8).lt.0.) stop 'c1 is unphysical. c1 should never be < 0!'
    !
    ! Initialise parameters for the Clenshaw reccurance formaular: see http://www.foo.be/docs-free/Numerical_Recipe_In_C/c5-8.pdf
    d=0
    dd=0
    z2=2.0_rk*z
    do j=7,2,-1  ! Clenshaw's recurrence.
      sv=d
      d=z2*d-dd+c(j)
      dd=sv
    enddo
    !
    ! Compute the Chebyshev expansion as a fucntion of z
    cheb_expansion=z*d-dd+c(1) ! not correct so far
    ! Compute the r-dependent, empirical term chi that multiplies the Chebyshev expansion
    chi=(1-exp(-c(9)*x))**3/(sqrt((x**2-c(10)**2)**2+c(11)**2)*sqrt((x**2-c(12)**2)**2+c(13)**2)) ! correct
    !
    ! Compute the irregular Chebysehv expansion DMC
    f=chi*cheb_expansion
    !
  end function irreg_chebyshev_DMC
  !
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


  function poten_two_coupled_EMOs_Lorentz(r,parameters) result(f)
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
    !
    a  = poten_lorentzian_polynom(r,parameters(nparams1+nparams2+1:nparams1+nparams2+nparams3))
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
  end function poten_two_coupled_EMOs_Lorentz


  recursive function poten_diabatic_Lorentz_TWO_EMOS(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: f1,f2,gamma,r0,f
    integer(ik)            :: nparams1,nparams2,nparams0
    !
    nparams0 = 2
    nparams1 = parameters(2+8)+9
    nparams2 = parameters(2+nparams1+8)+9
    !
    f1 = poten_EMO(r,parameters(2+1:2+nparams1))
    f2 = poten_EMO(r,parameters(2+nparams1+1:2+nparams1+nparams2))
    !
    gamma = parameters(1)
    r0 =  parameters(2)
    !
    if (abs(r-r0)>sqrt(small_)) then 
      !
      f = gamma/(r-r0)*0.5_rk*(f2-f1)
    else
      f=0.5_rk*( poten_diabatic_Lorentz_TWO_EMOS(r+0.001_rk,parameters)+poten_diabatic_Lorentz_TWO_EMOS(r-0.001_rk,parameters) )
    endif
    !
  end function poten_diabatic_Lorentz_TWO_EMOS


  function poten_two_coupled_EMOs_SqrtLorentz(r,parameters) result(f)
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
    !
    a  = poten_sqrt_lorentzian_polynom(r,parameters(nparams1+nparams2+1:nparams1+nparams2+nparams3))
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
  end function poten_two_coupled_EMOs_SqrtLorentz



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

  function poten_two_coupled_EMOs_with_EMO(r,parameters) result(f)
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
    !
    a  = poten_EMO(r,parameters(nparams1+nparams2+1:nparams1+nparams2+nparams3))
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
  end function poten_two_coupled_EMOs_with_EMO


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
  !
  !
  function poten_lorentzian_surkus_polynom(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: r0,gamma,z,f0,f
    integer(ik)            :: k,p,N
    !
    N = size(parameters)
    !
    gamma = parameters(1)
    r0 = parameters(2)
    z = 0.0_rk
    !
    if (N>3) then 
       p = int(parameters(3),ik)
       z = (r**p-r0**p)/(r**p+r0**p)
    endif
    !
    f0 = 1.0
    do k=4,N
     f0 = f0 + parameters(k)*z**(k-3)
    enddo
    !
    f = 0.5_rk*gamma*f0/( (r-r0)**2+gamma**2 )
    !
  end function poten_lorentzian_surkus_polynom

  !
  ! A sqrt(lorentzian) function for the couplings between diabatic curves
  !
  function poten_sqrt_lorentzian_polynom(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y0,r0,w,a,z,f0,f
    integer(ik)            :: k,N,p
    !
    N = size(parameters)
    !
    y0 = parameters(1)
    r0 = parameters(2)
    p  = parameters(3)
    !
    w = parameters(4)
    a = parameters(5)
    !
    z = (r**p-r0**p)/(r**p+r0**p)
    !
    f0 = a
    do k=6,N
     f0 = f0 + parameters(k)*z**(k-5)
    enddo
    !
    f = y0+f0*sqrt(2.0_rk*pi*( w/( 4.0_rk*(r-r0)**2+w**2 ) ))
    !
  end function poten_sqrt_lorentzian_polynom


  !
  ! An angle used for a unitary transformation to the diabatic representaion 
  ! for the lorentzian case
  !
  function beta_diabatic_lorentzian(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: gamma0,r0,f,t
    !
    gamma0 = parameters(1)
    r0 = parameters(2)
    !
    t = (r-r0)/gamma0
    !
    f=0.5_rk*atan(t)+pi/4.0_rk;
    !
  end function beta_diabatic_lorentzian


  !
  ! An angle used for a unitary transformation to the diabatic representaion 
  ! for the Laplacian case
  !
  function beta_diabatic_laplacian(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: gamma0,r0,f
    !
    gamma0 = parameters(1)
    r0 = parameters(2)
    !
    if (r<r0-small_) then 
      f= Pi/4.0_rk*exp( (r-r0)/gamma0 )
    elseif (r>r0+small_) then 
      f= Pi/2.0_rk-Pi/4.0_rk*exp(-(r-r0)/gamma0 )
    else
      f= Pi/4.0_rk
    endif
    !
  end function beta_diabatic_laplacian


  !
  ! An angle used for a unitary transformation to the diabatic representaion 
  ! for the Laplacian case
  !
  function beta_diabatic_lorentzian_laplacian(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: r0,f,tLor,tLap,gammaLap,gammaLor,betaLor,betaLap
    !
    gammaLor = parameters(1)
    r0 = parameters(2)
    gammaLap = parameters(3)
    !
    ! special case of the option ratio between gammas; 
    ! it is used to define gammaLap for gammaLap=0
    if (gammaLap<small_) gammaLap = 1.397_rk*gammaLor
    !
    tLor = (r-r0)/gammaLor
    tLap = (r-r0)/gammaLap
    !
    betaLor=0.5_rk*atan(tLor)+pi/4.0_rk
    !
    if (r<r0-small_) then 
      betaLap= Pi/4.0_rk*exp( tLap )
    elseif (r>r0+small_) then 
      betaLap= Pi/2.0_rk-Pi/4.0_rk*exp(-tLap )
    else
      betaLap= Pi/4.0_rk
    endif
    !
    ! beta is a geometric average of betaLor and betaLap
    !
    f = 0.5_rk*asin(sqrt(sin(2.0_rk*betaLor)*sin(2.0_rk*betaLap)))
    !
  end function beta_diabatic_lorentzian_laplacian


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
     z    = tanh(c1*r - c2/r)
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
      f = (1.0_rk - exp(-d1*r/n - d2*r*r/sqrt(real(n,rk))))**(n + 2)
      !
      end function DD
      !
  end function potential_stolyarov_CO_X_UBO
  !
  !
  function poten_EHH(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: y,v0,r0,de,f,alpha,q,phi,c
    integer(ik)            :: k,N
    !
    v0 = parameters(1)
    r0 = parameters(2)
    ! Note that the De is relative the absolute minimum of the ground state
    De = parameters(3)-v0
    !
    alpha = parameters(4)
    c     = parameters(5)
    !
    N     = size(parameters) - 5
    !
    if (size(parameters)<5) then
      write(out,"('poten_EHH: Illegal number of parameters in EHH, no beta -present')")
      print*,parameters(:)
      stop 'poten_EHH: Illegal number of parameters, no beta'
    endif
    !
    q = alpha*(r-r0)
    !
    phi = 1.0_rk
    do k=1,N
     phi = phi + parameters(k+5)*q**k
    enddo
    !
    y  = exp(-q)
    !
    f = de*((1.0_rk-y)**2+c*q**3*phi*y*y)+v0
    !
  end function poten_EHH
  !
  !
  ! Irregular DMF by Medvedev Opt. spectrosc. 130, 1334 (2022)
  !
  function dipole_medvedev_sing(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: alpha,beta,z,y,f,r1,r2,b1,b2,s
    integer(ik)            :: N,k,i
    !
    alpha = parameters(1)
    beta  = parameters(2)
    r1    = parameters(3)
    b1    = parameters(4)
    r2    = parameters(5)
    b2    = parameters(6)
    n     = int(parameters(7))
    !
    z=1.0_rk-2.0_rk*exp(-r*beta)
    y=1.0_rk-exp(-r*alpha)
    !
    k = size(parameters)-8
    !    
    s = parameters(8)
    do i = 1, k
        s = s+parameters(i+8)*z**i
    end do
    !
    f = s*y**5/&
        sqrt( (r**2-r1**2)**2+b1**2 )/&
        sqrt( (r**2-r2**2)**2+b2**2 )
    !
  end function dipole_medvedev_sing
  !
  ! does not do anything 
  !
  function function_dummy(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    real(rk)               :: f
    !
    f = 0
    !
  end function function_dummy
  !
  !
  function polynomial_coupled_PECs(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    !
    real(rk)    :: V1,V2,VD,Discr,E(2)
    integer(ik) :: icomponent,Ntot
    !
    real(rk)   :: f
    !
    f = 0
    !
    V1 = function_V1(r,parameters(1:N1))
    V2 = function_V2(r,parameters(N1+1:N1+N2))
    VD = function_VD(r,parameters(N1+N2+1:N1+N2+N3))
    !
    ! to diagobalize the 2x2 matrix 
    !/     \
    ! V1 VD |
    ! VD V2 |
    !\     /
    ! we solve a quadratic equation with the discriminant 
    !
    Discr = V1**2-2.0_rk*V1*V2+V2**2+4.0_rk*VD**2
    !
    if (discr<-small_) then
      write(out,"('COUPLED: discriminant is negative for polynomial_coupled_PECs, r= ',f17.8)") r
      stop 'COUPLED: discriminant is negative'
    endif
    !
    E(1)=0.5_rk*(V1+V2)-0.5_rk*sqrt(Discr)
    E(2)=0.5_rk*(V1+V2)+0.5_rk*sqrt(Discr)
    !
    ! we use the component as specified in the last parameter:
    !
    Ntot = size(parameters)
    !
    icomponent = parameters(Ntot)
    !
    f = E(icomponent)              
    !
  end function polynomial_coupled_PECs
  !
  !
  ! A coupled PEC with an arbitrary number of states 
  !
  function polynomial_multi_coupled_PECs(r,parameters) result(f)
    !
    use lapack,only : lapack_syev
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    !
    real(rk),allocatable   :: h(:,:),e(:)
    integer(ik) :: icomponent,Ntot,Ndim,N,i,i1,i2,alloc
    !
    real(rk)   :: f
    !
    f = 0
    !
    Ndim = ( nint( sqrt(1.0_rk+8_rk*real(N_multi_subfunctions,rk)))-1 )/2
    !
    allocate(h(Ndim,Ndim),e(Ndim),stat=alloc)
    if (alloc/=0) then
      write (out,"('Error alloac polynomial_multi_coupled_PECs:  h and e ',i0)") alloc
      stop 'Error alloac polynomial_multi_coupled_PECs - alloc'
    end if
    !
    i = 0
    Ntot = 0 ! accumulated size
    !
    do i1 =1,Ndim
       i = i + 1
       N = nsize_multi(i)
       h(i1,i1) = function_multi(i)%f(r,parameters(Ntot+1:Ntot+N))
       Ntot = Ntot + N
    enddo
    !
    do i1 =1,Ndim
       do i2 =i1+1,Ndim
         i = i + 1
         N = nsize_multi(i)
         h(i2,i1) = function_multi(i)%f(r,parameters(Ntot+1:Ntot+N))
         h(i1,i2) = h(i2,i1)
         Ntot = Ntot + N
       enddo
    enddo
    !
    ! to diagobalize the NxN matrix 
    !
    call lapack_syev(h,e)
    !
    ! we use the component as specified in the last parameter:
    !
    Ntot = size(parameters)
    !
    icomponent = parameters(Ntot)
    !
    f = e(icomponent)  
    !
    deallocate(h,e)            
    !
  end function polynomial_multi_coupled_PECs
  !
  recursive function polynomial_coupled_PECs_beta(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    !
    real(rk)    :: V1,V2,VD,Discr,E(2),beta
    integer(ik) :: icomponent,Ntot
    !
    real(rk)   :: f
    !
    f = 0
    !
    V1 = function_V1(r,parameters(1:N1))
    V2 = function_V2(r,parameters(N1+1:N1+N2))
    !
    beta = function_beta(r,parameters(N1+N2+1:N1+N2+N3))
    !
    ! Define the diabatic coupling: VD = 0.5*tan(2*gamma)*(V1-V2)
    !
    if (abs(beta-pi*0.25_rk)>sqrt(small_)) then 
      !
      VD = 0.5_rk*tan(2.0_rk*beta)*(V2-V1)
      !
    else
      !
      VD = polynomial_coupled_PECs_beta(r-sqrt(small_),parameters) 
      !
    endif
    !
    ! to diagobalize the 2x2 matrix 
    !/     \
    ! V1 VD |
    ! VD V2 |
    !\     /
    ! we solve a quadratic equation with the discriminant 
    !
    Discr = V1**2-2.0_rk*V1*V2+V2**2+4.0_rk*VD**2
    !
    if (discr<-small_) then
      write(out,"('COUPLED: discriminant is negative in polynomial_coupled_PECs_beta r=',f15.8)") r
      stop 'COUPLED: discriminant is negative in polynomial_coupled_PECs_beta'
    endif
    !
    E(1)=0.5_rk*(V1+V2)-0.5_rk*sqrt(Discr)
    E(2)=0.5_rk*(V1+V2)+0.5_rk*sqrt(Discr)
    !
    ! we use the component as specified in the last parameter:
    !
    Ntot = size(parameters)
    !
    icomponent = parameters(Ntot)
    !
    f = E(icomponent)              
    !
  end function polynomial_coupled_PECs_beta


  !
  recursive function coupled_transition_curves_beta(r,parameters) result(f)
    !
    real(rk),intent(in)    :: r             ! geometry (Ang)
    real(rk),intent(in)    :: parameters(:) ! potential parameters
    !
    real(rk)    :: f1,f2,g(2),beta
    integer(ik) :: icomponent,Ntot
    !
    real(rk)   :: f
    !
    f = 0
    !
    f1 = function_V1(r,parameters(1:N1))
    f2 = function_V2(r,parameters(N1+1:N1+N2))
    !
    beta = function_beta(r,parameters(N1+N2+1:N1+N2+N3))
    !
    ! unitary transforamtion of one state in the transition property (e.g. a dipole between a diabatised and single states)
    !
    g(1) = cos(beta)*f1-sin(beta)*f2
    g(2) = sin(beta)*f1+cos(beta)*f2
    !
    ! we use the component as specified in the last parameter:
    !
    Ntot = size(parameters)
    !
    icomponent = parameters(Ntot)
    !
    f = g(icomponent)              
    !
  end function coupled_transition_curves_beta

  !
end module functions
