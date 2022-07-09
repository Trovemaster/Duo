module diatom_module
  !
  use accuracy
  use timer
  use functions, only : analytical_fieldT
  use symmetry,  only : sym,SymmetryInitialize
  use Lobatto,   only : LobattoAbsWeights,derLobattoMat
  use me_numer,  only : ME_numerov
  !
  implicit none
  !                     by Lorenzo Lodi
  !                     Specifies variables needed for storing a temporary (scratch) copy
  !                     of the input (useful for jumping around).
  !
  integer             :: ierr
  character(len=wl)   :: line_buffer
  integer,parameter   :: max_input_lines=500000  ! maximum length (in lines) of input. 500,000 lines is plenty..
  !                                              ! to avoid feeding in GB of data by mistake.
  !
  ! main method used for solving the J=0 problem, i.e. to solve the one-dimentional Schrodinger equation.
  ! character(len=wl)   :: solution_method = "5POINTDIFFERENCES"  ! kinetic energy approximated by 5-point finite-diff. formula
  character(len=wl)   :: solution_method = "SINC"  ! kinetic energy approximated by sinc DVR
  !
  ! Type to describe different terms from the hamiltonian, e.g. potential energy, spin-orbit, <L^2>, <Lx>, <Ly> functions.
  !
  integer(ik),parameter   :: verbose=5
  integer(ik),parameter   :: Nobjects = 32  ! number of different terms of the Hamiltonian 
  !                                          (poten,spinorbit,L2,LxLy,spinspin,spinspino,bobrot,spinrot,diabatic,
  !                                             lambda-opq,lambda-p2q)
  !
  ! In order to add a new field:
  ! 1. Change Nobjects
  ! 2. Add the name of the field to the CASE line where all fields are listed: case("SPIN-ORBIT","SPIN-ORBIT-X"....
  ! 3. Add the CASE-section describing the input of the new field 
  ! 4. Add a similar CASE-section describing the input of the new field to the CASE("ABINITIO") part
  ! 5. Introduce a new name for the new field (object). 
  ! 6. Define new array in type(fieldT),pointer .....
  ! 7. Add a new case to all cases where the fieds appear: at the end of the input subroutine, 
  !    two times in in map_fields_onto_grid, in duo_j0, in sf_fitting
  ! 8. Add corresposponding check_and_print_coupling in map_fields_onto_grid
  ! 9. Add the corresponding name of the field to  "use diatom_module,only" in "refinement"
  !
  ! Current list of fields:
  ! 
  !        case (1) poten(iterm)
  !        case (2) spinorbit(iterm)
  !        case (3) l2(iterm)
  !        case (4) lxly(iterm)
  !        case (5) spinspin(iterm)
  !        case (6) spinspino(iterm)
  !        case (7) bobrot(iterm)
  !        case (8) spinrot(iterm)
  !        case (9) diabatic(iterm)
  !        case (10) lambdaopq(iterm)
  !        case (11) lambdap2q(iterm)
  !        case (12) lambdaq(iterm)
  !        case (14 - 20) reserved
  !        case (21) hfcc1(1) for Fermi contact, bF
  !        case (22) hfcc1(2) for nuclear spin - orbit, a
  !        case (23) hfcc1(3) for nuclear dipole - spin dipole, c
  !        case (24) hfcc1(4) for nuclear dipole - spind dipole, d
  !        case (25) hfcc1(5) for nuclear spin - rotaton, cI
  !        case (26) hfcc1(6) for electric dipole, eQq0
  !        case (27) hfcc1(7) for electric diople, eQq2
  !        case (28) reserved
  !        case(Nobjects-3) quadrupoletm(iterm)
  !        case(Nobjects-2) abinitio(iterm)
  !        case(Nobjects-1) brot(iterm)
  !        case(Nobjects) dipoletm(iterm)
  !
  !
  ! Lorenzo Lodi: it is necessary to repeat the character length specification in the array contructor for Fortran2003 conformance
  character(len=wl),parameter :: CLASSNAMES(1:Nobjects)  = (/ character(len=wl):: "POTEN","SPINORBIT","L2", "L+","SPIN-SPIN",&
                                                             "SPIN-SPIN-O","BOBROT","SPIN-ROT","DIABATIC","LAMBDAOPQ", &
                                                             "LAMBDAP2Q","LAMBDAQ", &
                                                             "NAC", "", "", "", "", "", "", "", & ! reserved
                                                             "HFCC-BF-1", "HFCC-A-1", &
                                                             "HFCC-C-1", "HFCC-D-1", &
                                                             "HFCC-CI-1", &
                                                             "HFCC-EQQ0-1", "HFCC-EQQ2-1", &
                                                             "", & ! reserved
                                                             "ABINITIO","BROT","DIPOLE","QUADRUPOLE"/)
  !
  ! Lorenzo Lodi
  ! What follows is a bunch of temporary variables needed for computation of derivatives of pecs and other curves
  real(kind=rk) :: fmmm,fmm,fm,f0,fp,fpp,fppp ! value of PEC on grid stencil
  real(kind=rk) :: der1, der2, der3, der4     ! derivatives
  real(kind=rk), parameter :: thresh_rel_err_min_x = 1e-13_rk ! threshold (relative error) for convergence of minimum of PEC
  integer(kind=ik), parameter       :: max_iter_min_search=20  ! maximum number of iterations for finding minimum of PEC
  real(kind=rk), parameter :: h=2e-3_rk ! `h' is the step size in ang. for computing numerical derivatives of pecs
  real(kind=rk) :: x0, x1 ! generic tmp variables

  ! Lorenzo Lodi -- text variables containing the symbols of the two atoms. The following formats are supported:
  !              1) Chemical symbol. E.g.  C, O, Na, Sc, Fe ....
  !              2) Chemical symbol - atomic number, e.g.: C-14, O-16, Na-23
  character(len=15) :: symbol1="Undefined", symbol2="Undefined"
  !
  !
  type symmetryT
    !
    integer(ik) :: gu = 0
    integer(ik) :: pm = 0
    !
  end type symmetryT
  !
  ! type describing the parameter geneology
  !
  type linkT
    !
    integer(ik) :: iobject
    integer(ik) :: ifield
    integer(ik) :: iparam
    !
  end type linkT
  !
  type braketT
    integer(ik) :: ilambda = 0.0_rk
    integer(ik) :: jlambda = 0.0_rk
    real(rk) :: sigmai = 0.0_rk
    real(rk) :: sigmaj = 0.0_rk
    real(rk) :: value = 0.0_rk
  end type braketT
  !
  ! files with the eigenvectors 
  !
  type  eigenfileT
    !
    character(len=cl)  :: dscr       ! file with fingeprints and descriptions of each levels + energy values
    character(len=cl)  :: primitives ! file with the primitive quantum numbres   
    character(len=cl)  :: vectors    ! eigenvectors stored here 
    !
  end type  eigenfileT 
  !
  type weightT
    character(cl) :: wtype = 'GRID'
    real(rk) :: alpha = 0.001_rk
    real(rk) :: Vtop = 10000.0_rk
  end type weightT
  !
  integer(ik), parameter :: bad_value=-10000    ! default `impossible' value for some quantum numbers
  !
  type fieldT
    !
    character(len=cl)    :: name='(unnamed)' ! Identifying name of the function (default to avoid garbled outputs)
    character(len=cl)    :: type='NONE'  ! Identifying type of the function
    character(len=cl)    :: class='NONE' ! Identifying class of the function (poten, spinorbit,dipole,abinitio etc)
    !
    ! variable used for GRID curves to specify interpolation and extrapolation  options
    ! 
    character(len=cl)    :: interpolation_type='CUBICSPLINES'
    !
    !character(len=cl)    :: interpolation_type='QUINTICSPLINES'
    !
    integer(ik)          :: iref         ! reference number of the term as given in input (bra in case of the coupling)
    integer(ik)          :: jref         ! reference number of the coupling term as given in input (ket in case of the coupling)
    integer(ik)          :: istate       ! the actual state number (bra in case of the coupling)
    integer(ik)          :: jstate       ! the actual state number (ket in case of the coupling)
    integer(ik)          :: Nterms       ! Number of terms or grid points
    integer(ik)          :: Lambda =bad_value     ! identification of the electronic state Lambda
    integer(ik)          :: Lambdaj=bad_value     ! identification of the electronic state Lambda (ket)
    integer(ik)          :: multi        ! identification of the electronic spin (bra for the coupling)
    integer(ik)          :: jmulti       ! identification of the ket-multiplicity
    real(rk)             :: sigmai=real(bad_value,rk)    ! the bra-projection of the spin
    real(rk)             :: sigmaj=real(bad_value,rk)    ! the ket-projection of the spin
    real(rk)             :: spini        ! electronic spin (bra component for couplings)
    real(rk)             :: spinj        ! electronic spin of the ket vector for couplings
    real(rk)             :: omegai       ! projection of J
    real(rk)             :: omegaj       ! projection of J
    integer(ik)          :: iomega       ! additional QN
    integer(ik)          :: jomega       ! additional QN
    integer(ik)          :: ilevel       ! additional QN
    integer(ik)          :: jlevel       ! additional QN
    complex(rk)          :: complex_f=(1._rk,0._rk)  ! defines if the term is imaginary or real
    real(rk)             :: factor=1.0_rk      ! defines if the term is imaginary or real
    real(rk)             :: fit_factor=1.0     ! defines if the term is imaginary or real
    real(rk),pointer     :: value(:)=>null()   ! Expansion parameter or grid values from the input
    type(symmetryT)      :: parity       ! parity of the electronic state as defined by the molecular inversion (g,u), 
     !                                        or laboratory inversion (+,-)
    real(rk),pointer     :: gridvalue(:) ! Expansion parameter or a grid value on the grid used inside the program
    real(rk),pointer     :: weight(:)    ! fit (1) or no fit (0)
    real(rk),pointer     :: grid(:)      ! grid value
    real(rk),pointer     :: matelem(:,:) ! matrix elements
    real(rk)             :: refvalue = 0 ! reference value will be used as a shift to be applied to the ab initio function used 
    !                                         for the fitting constraints
    character(len=cl),pointer :: forcename(:) ! The parameter name
    integer(ik)          :: nbrakets=0   ! total number of different combinations of lambda and sigma in matrix elements (maximum 4)
    type(braketT)        :: braket(4)   ! here all possible combinations of values <+/-lambda,+/-sigma|A|+/-lambdaj,+/-sigma> 
    !                                                                                                                  can be listed
    integer(ik)          :: imin         ! grid point index at which a potential energy curve is minimum
    real(rk)             :: Vimin        ! V(imin), minimum of potential energy curve on the grid in cm-1
    real(rk)             :: rimin        ! r_imin, r of the minimum of potential energy curve on the grid
    integer(ik)          :: imax         ! grid point index at which a potential energy curve is maximum
    real(rk)             :: Vimax        ! V(imax), maximum of potential energy curve on the grid in cm-1

    logical              :: zHasMinimum =.true.  ! whether a potential energy curve has a minimum
    real(rk)             :: re   ! pecs only: minimum for given pec, in angstroms.
    real(rk)             :: V0   ! pecs only: V(re)    , in cm-1
    real(rk)             :: V1   ! pecs only: V'(re)   , in cm-1 / angstroms (should be zero).
    real(rk)             :: V2   ! pecs only: V''(re)  , in cm-1  / angstroms**2.
    real(rk)             :: V3   ! pecs only: V'''(re) , in cm-1  / angstroms**3.
    real(rk)             :: V4   ! pecs only: V''''(re), in cm-1  / angstroms**4.
    real(rk)             :: we   ! pecs only: harmonic vibrational freq, cm-1.
    real(rk)             :: B0   ! pecs only: v=0 vibrational constant vibrational, cm-1.
    real(rk)             :: xe   ! pecs only: anharmonicity constant (pure number)
    real(rk)             :: alphae   ! pecs only: Coriolis coupling constant, cm-1.
    real(rk)             :: Debar    ! pecs only: Centrifugal distortion constant, cm-1.
    real(rk)             :: Y00      ! pecs only: ZPE correction, cm-1.
    real(rk)             :: Omega_min  !  pecs only: minimum physically possible value for |Omega|=|Lambda+Sigma|
    real(rk)             :: approxEJ0    !  pecs only: approximate J=0, v=0 energy (no couplings)
    real(rk)             :: approxEJmin  !  pecs only: approximate J=Jmin, v=0 energy (no couplings)


    procedure (analytical_fieldT),pointer, nopass :: analytical_field => null()
    type(linkT),pointer   :: link(:)       ! address to link with the fitting parameter in a different object in the fit
    logical               :: morphing = .false.    ! When morphing is the field is used to morph the ab initio couterpart
    !                                                towards the final object
    logical               :: molpro = .false.      ! The object is given in the molpro representaion (|x> and |y>) 
    integer(ik)           :: ix_lz_y = 1000        ! The value of the matrix element (-I)<x|Lz|y> for i-state, 1000 is for undefined 
    integer(ik)           :: jx_lz_y = 1000        ! The value of the matrix element (-I)<x|Lz|y> for j-state, 1000 is for undefined 
    type(weightT)         :: weighting             ! When morphing is the field is used to morph the ab initio couterpart
    character(len=cl)     :: integration_method='DVR-SINC' ! Identifying the type of the integration method used specifically for this field, DVR-SINC default
    !
    real(rk)              :: adjust_val = 0.0_rk
    real(rk)              :: asymptote = 0.0_rk      ! reference asymptote energy used e.g. in the renormalisation of the continuum wavefunctions to sin(kr)
    logical               :: adjust = .false.        ! Add constant adjust_val to all fields
  end type fieldT
  !
  type FieldListT
    TYPE(fieldT), POINTER :: field(:)
    INTEGER :: num_field = 0
    PROCEDURE(), POINTER, NOPASS :: hfcc_matrix_element => NULL()
  end type FieldListT

  type jobT
      logical             :: select_gamma(4) ! the diagonalization will be done only for selected gamma's
      integer(ik)         :: nroots(4)=1e9_rk ! number of the roots to be found in variational diagonalization with syevr
      integer(ik)         :: maxiter = 1000  ! maximal number of iterations in arpack
      real(rk)            :: tolerance = 0.0_rk   ! tolerance for arpack diagonalization, 0 means the machine accuracy
      real(rk)            :: upper_ener = 1e9_rk  ! upper energy limit for the eigenvalues found by diagonalization with syevr
      real(rk)            :: thresh = -1e-18_rk   ! thresh of general use
      real(rk)            :: zpe=0.0_rk             ! zero-point-energy
      logical             :: shift_to_zpe = .true. ! 
      character(len=cl)   :: diagonalizer = 'SYEV'
      character(len=cl)   :: molecule = 'H2'
      character(len=cl)   :: contraction = 'VIB' ! contraction
      real(rk),pointer    :: vibenermax(:)    !     contraction parameter: energy
      integer(ik),pointer :: vibmax(:)           !     contraction parameter: vibration quantum
      real(rk)            :: potmin              ! absolute minimum of the PEC with the lowest rotational-vibrational state
      logical             :: zShiftPECsToZero = .true.
      logical             :: zEchoInput = .true.
      logical             :: zExclude_JS_coupling =.false. ! If set to true will disable J.S coupling (aka S-uncoupling)
      integer(ik)         :: total_parameters =0  !  total number of parameters used to define different hamiltonian fields
      real(rk)            :: degen_threshold = 1e-6_rk
      real(rk),pointer    :: j_list(:)     ! J values processed
      integer(ik)         :: nJ = 1        ! Number of J values processed 
      character(len=cl)   :: IO_eigen = 'NONE'   ! we can either SAVE to or READ from the eigenfunctions from an external file
      character(len=cl)   :: IO_dipole = 'NONE'  ! we can either SAVE to or READ from an external file the dipole moment 
      character(len=cl)   :: basis_set = 'NONE'   ! we can keep the vibrational basis functions for further usage 
      !                                                matrix elements on the contr. basis 
      type(eigenfileT)    :: eigenfile
      character(len=cl)   :: symmetry = 'CS(M)'    ! molecular symmetry
      real(rk)   :: diag_L2_fact = 1._rk    ! specifies the convention used for the diagonal contribution
                                            !  due to L^2 = LxLx+LyLy+LzLz
      logical,pointer     :: isym_do(:)     ! process or not the symmetry in question
      logical             :: intensity      ! switch on the intensity calculations
      logical             :: print_vibrational_energies_to_file = .false. ! if .true. prints to file
      logical             :: print_rovibronic_energies_to_file = .false. ! if .true. prints to file
      logical             :: print_pecs_and_couplings_to_file = .false. ! if .true. prints to file
      logical             :: assign_v_by_count = .false.
      logical             :: legacy_version = .false.
      !
  end type jobT
  !
  type gridT
      integer(ik)   :: npoints = 1000       ! grid size
      real(rk)      :: rmin = 1.0,rmax=3.00 ! range of the grid
      real(rk)      :: step = 1e-2          ! step size
      real(rk)      :: alpha = 0.2          ! grid parameter
      real(rk)      :: re = 1.0             ! grid parameter
      integer(ik)   :: nsub = 0             ! grid type parameter (0=uniformly spaced)
      real(rk),pointer :: r(:)=>null()      ! the molecular geometry at the grid point
  end type gridT

  type quantaT
    real(rk) :: I1 ! nuclear spin 1
    real(rk) :: F1  ! \hat{F1} = \hat{I1} + \hat{J}
    INTEGER(ik) :: index_F1 ! index of the F1 in F1_list
    real(rk)     :: Jrot       ! J - real
    integer(ik)  :: irot       ! index of the J value in J_list
    integer(ik)  :: istate     ! e-state
    integer(ik)  :: imulti     ! imulti = 1,..,multiplicity
    real(rk)     :: sigma      ! spin-projection = -spin,-spin+1,..,spin-1,spin
    real(rk)     :: omega      ! sigma+lambda
    real(rk)     :: spin       ! spin
    integer(ik)  :: ilambda    ! ilambda
    integer(ik)  :: v  = 0     ! vibrational quantum
    integer(ik)  :: ivib = 1   ! vibrational quantum number counting all vib-contracted states
    integer(ik)  :: ilevel = 1  ! primitive quantum
    integer(ik)  :: iroot       ! running number
    integer(ik)  :: iF1_ID       ! running number within the same F and parity
    integer(ik)  :: iJ_ID       ! running number within the same J
    integer(ik)  :: iparity = 0
    integer(ik)  :: igamma = 1
    integer(ik)  :: iomega = 1   ! countig number of omega
    character(len=cl) :: name    ! Identifying name of the  function
    logical :: bound = .true.    ! is this state bound or unbound 
  end type quantaT
  !
  type eigenT
      real(rk),pointer :: vect(:,:)      ! the eigenvector in the J=0 contracted representaion
      real(rk),pointer :: val(:)         ! the eigenvalue
      type(quantaT),pointer ::quanta(:)  ! the quantum numbers
      integer(ik)      :: Nlevels
      integer(ik)      :: Ndimen
  end type eigenT
  !
  type basisT
      type(quantaT),pointer :: icontr(:)  ! the quantum numbers
      integer(ik)           :: Ndimen
  end type basisT
  !
  type actionT
     !
     logical :: fitting       = .false.
     logical :: intensity     = .false.
     logical :: frequency     = .false.
     logical :: matelem       = .false.
     logical :: raman         = .false.
     logical :: quadrupole    = .false.
     logical :: RWF           = .false.
     logical :: hyperfine     = .false.
     logical :: save_eigen_J  = .false.
     !
  end type actionT
  !
  type obsT
     !
     real(rk)      :: Jrot
     real(rk)      :: Jrot_      ! J - real (lower)
     integer(ik)   :: irot       ! index of the J value in J_list
     integer(ik)   :: irot_
     integer(ik)   :: symmetry
     integer(ik)   :: N
     integer(ik)   :: N_         ! N lower
     integer(ik)   :: iparity
     integer(ik)   :: iparity_   ! parity +/- (0,1) of the lower state lower
     real(rk)      :: energy
     real(rk)      :: frequency
     real(rk)      :: weight
     type(quantaT) :: quanta
     type(quantaT) :: quanta_  ! quantum numbers of the lower state
     !
  end type obsT
  type thresholdsT
     real(rk) :: intensity    = -1e0    ! threshold defining the output intensities
     real(rk) :: linestrength = -1e0    ! threshold defining the output linestrength
     real(rk) :: dipole       = -1e0    ! threshold defining the output linestrength
     real(rk) :: coeff        = -1e0    ! threshold defining the eigenfunction coefficients
                                        ! taken into account in the matrix elements evaluation.
  end type thresholdsT
  !
  type landeT
    real, allocatable :: ar(:)
  end type landeT
  !
  type IntensityT
     logical             :: do = .false.     ! process (.true.) or not (.false.) the intensity (or TM) calculations 
     character(cl)       :: action           ! type of the intensity calculations:
                                             ! absorption, emission, tm (transition moments),
                                             !  raman, and so on. 
     real(rk)            :: temperature      ! temperature in K
     real(rk)            :: part_func=0      ! partition function 
     real(rk)            :: ZPE=0            ! zero point energy
     type(thresholdsT)   :: threshold        ! different thresholds
     real(rk),pointer    :: gns(:)           ! nuclear stat. weights
     integer(ik),pointer :: isym_pairs(:)    ! numbers defining symmetry pairs with allowed transitions, analogous to gns
     real(rk)            :: freq_window(1:2) ! frequency window (1/cm)
     real(rk)            :: erange_low(1:2)  ! energy range for the lower state
     real(rk)            :: erange_upp(1:2)  ! energy range for the upper state
     real(rk)            :: J(1:2)           ! range of J-values, from..to; in order to avoid double counting of transitions
                                             ! in the calculations it is always assumed that 
                                             ! J1<=J_lower<=J2 and J1<=J_upper<J2;
                                             !
     type(quantaT) :: lower                   ! lower state range of the quantun numbers employed 
     type(quantaT) :: upper                   ! upper state range of the quantun numbers employed 
                                             ! in intensity calculations; (imode,1:2), 
                                             ! where 1 stands for the beginning and 2 for the end. 
     !
     integer(ik)         :: swap_size    = 0       ! the number of vectors to keep in memory
     character(cl)       :: swap = "NONE"          ! whether save the compacted vectors or read
     character(cl)       :: swap_file  ="compress" ! where swap the compacted eigenvectors to
     character(cl)       :: linelist_file="NONE"   ! filename for the line list (filename.states and filename.trans)
     integer(ik)         :: int_increm = 1e9       ! used to print out the lower energies needed to select int_increm intensities
     real(rk)            :: factor = 1.0d0         ! factor <1 to be applied the maxsize of the vector adn thus to be shrunk 
     logical             :: matelem =.false.       ! switch for the line-strenth-type matelems  (matrix elements of the dipole moment)
     logical             :: lande_calc = .false.   ! checks whether calculation for Lande should be conducted
     logical             :: overlap = .false.      ! print out overlap integrals (Franck-Condon)
     logical             :: tdm      = .true.      ! print out dipole transition moments 
     logical             :: tqm      = .true.      ! print out quadrupole transition moments
     integer(ik)         :: Npoints = -1           ! used for cross sections grids 
     real(rk)            :: gamma = 0.05_rk        ! Lorentzian FWHM, needed for cross-sections
     integer(ik)         :: N_RWF_order  = 1       ! Expansion order of the matrix fraction needed for RWF 
     character(cl)       :: RWF_type="GAUSSIAN"    ! Type of RWH
     logical             :: renorm = .false.       ! renormalize the continuum/unbound wavefunctions to sin(kr) for r -> infty
     logical             :: bound = .false.        ! filter bound states
     logical             :: unbound = .false.       ! filter and process unbound upper states only 
     !
 end type IntensityT
  !
  type matrixT
     !
     real(rk),pointer    :: matrix(:,:)
     integer(ik),pointer :: irec(:)
     !
  end type matrixT
  !
  type fittingT
     !
     logical              :: run
     integer(ik)          :: nJ = 1        ! Number of J values processed 
     real(rk),pointer     :: j_list(:)     ! J-values processed in the fit
     integer(ik)          :: iparam(1:2) = (/1,100000/)
     integer(ik)          :: itermax = 500
     integer(ik)          :: Nenergies = 1
     integer(ik)          :: parmax =0          ! total number of all parameters used
     real(rk)             :: factor = 1.0_rk
     real(rk)             :: target_rms = 1e-8
     real(rk)             :: robust = 0
     character(len=cl)    :: geom_file = 'pot.fit'
     character(len=cl)    :: output_file = 'fitting'
     character(len=cl)    :: fit_type = 'LINUR'      ! to switch between fitting methods.
     real(rk)             :: threshold_coeff    = -1e-18
     real(rk)             :: threshold_lock     = -1e-18
     real(rk)             :: threshold_obs_calc  = -1e-16
     real(rk)             :: zpe=0
     logical              :: shift_to_zpe = .true.   ! 
     real(rk)             :: fit_scaling=1.0_rk         ! scaling the fitting correction with this factor >0 and <1
     integer(ik)          :: linear_search = 0 ! use linear scaling to achieve better convergence with the Armijo condition
     type(obsT),pointer   :: obs(:)           ! experimental data
     !
     !type(paramT),pointer :: param(:)         ! fitting parameters
     !
  end type fittingT
  !
  type fieldmapT
    integer(ik)          :: Nfields
  end type fieldmapT
  !
  ! This type is designed for the Omega-state representations
  type Omega_gridT
      integer(ik)      :: Nstates    
      real(rk)         :: omega
      real(rk),pointer :: energy(:,:)=>null()      ! the Omega-potential energy at each grid point 
      real(rk),pointer :: vector(:,:,:)=>null()      ! the Omega-vector at each grid point in the lambda-sigma space 
      type(quantaT),pointer ::  QN(:)=>null()  ! quantum numbers 
      type(quantaT),pointer ::  basis(:)=>null()  ! basis set
  end type Omega_gridT
  !
  ! This type is  for the vibronic contacted eigensolutions in Omega represetation
  type contract_solT
      integer(ik)      :: Ndimen    
      real(rk),pointer :: Energy(:)
      real(rk),pointer :: vector(:,:)=>null()         ! the Omega-vector at each grid point in the lambda-sigma space 
      integer(ik),pointer :: ilevel(:)=>null()      ! level index
  end type contract_solT

  type F1_hyperfine_steup_T
      real(rk) :: I1
      real(rk) :: fit_bound_ratio
      character(cl) :: fit_optimization_algorithm
  end type F1_hyperfine_steup_T
  !
  integer, parameter :: trk        = selected_real_kind(12)
  integer,parameter  :: jlist_max = 500
  type(fieldT),pointer :: poten(:),spinorbit(:),l2(:),lxly(:),abinitio(:),dipoletm(:)=>null(),&
                          spinspin(:),spinspino(:),bobrot(:),spinrot(:),diabatic(:),lambdaopq(:),lambdap2q(:),lambdaq(:),nac(:)
  type(fieldT),pointer :: brot(:),quadrupoletm(:)
  !
  ! Fields in the Omega representation
  !
  integer(ik) :: Nspins,NLplus_omega,NSplus_omega,NSR_omega,NBob_omega,Nomegas,Np2q_omega,Nq_omega,NKin_omega,NBRot_omega
  type(Omega_gridT),allocatable :: omega_grid(:)
  integer(ik),allocatable :: iLplus_omega(:,:,:,:),iSplus_omega(:,:,:,:),iSR_omega(:,:,:,:),iBOB_omega(:,:,:)
  integer(ik),allocatable :: iP2Q_omega(:,:,:,:),iQ_omega(:,:,:,:),iKin_omega(:,:,:),iBRot_omega(:,:,:)
  !
  type(fieldT),pointer :: l_omega_obj(:)=>null(),s_omega_obj(:)=>null(),sr_omega_obj(:)=>null(),bob_omega_obj(:)=>null()
  type(fieldT),pointer :: p2q_omega_obj(:)=>null(),q_omega_obj(:)=>null(),kin_omega_obj(:)=>null(),brot_omega_obj(:)=>null()
  !
  type(jobT)   :: job
  type(gridT)  :: grid
  type(quantaT),allocatable :: quanta(:)
  integer(ik),allocatable   :: iquanta2ilevel(:,:,:)
  real(rk),allocatable      :: r(:), d2dr(:), r2sc(:),z(:)
  type(actionT)             :: action   ! defines different actions to perform
  type(fittingT)            :: fitting
  type(IntensityT)            :: Intensity
  !type(symmetryT)             :: sym
  !
  integer(ik)   :: nestates,Nspinorbits,Ndipoles,Nlxly,Nl2,Nabi,Ntotalfields=0,Nss,Nsso,Nbobrot,Nsr,Ndiabatic,&
                   Nlambdaopq,Nlambdap2q,Nlambdaq,Nnac,vmax,nQuadrupoles,NBrot,nrefstates = 1
  real(rk)      :: m1=-1._rk,m2=-1._rk ! impossible, negative initial values for the atom masses
  real(rk)      :: jmin,jmax,amass,hstep,Nspin1,Nspin2
  real(rk)      :: jmin_global
  !
  !type(fieldT),pointer :: refined(:)
  type(fieldmapT) :: fieldmap(Nobjects)
  type(eigenT),allocatable :: eigen(:,:)
  !
  type(basisT),allocatable :: basis(:)
  !
  logical :: gridvalue_allocated  = .false.
  logical :: fields_allocated  = .false.
  real(rk),parameter :: enermax = safe_max  ! largest energy allowed 

  real(rk), allocatable :: vibrational_contrfunc(:,:)  
  type(quantaT), allocatable :: vibrational_quantum_number(:)
  integer(ik) :: vibrational_totalroots
  type(F1_hyperfine_steup_T) :: F1_hyperfine_setup
  integer(ik), parameter :: GLOBAL_NUM_HFCC_OBJECT = 7
  type(FieldListT) :: hfcc1(GLOBAL_NUM_HFCC_OBJECT)
  !
  public ReadInput,poten,spinorbit,l2,lxly,abinitio,brot,map_fields_onto_grid,fitting,&
         jmin,jmax,vmax,fieldmap,Intensity,eigen,basis,Ndipoles,dipoletm,linkT,three_j,quadrupoletm,&
         l_omega_obj,s_omega_obj,sr_omega_obj,brot_omega_obj,p2q_omega_obj,q_omega_obj,kin_omega_obj
  !
  save grid, Intensity, fitting, action, job, gridvalue_allocated, fields_allocated, hfcc1
  !
  contains
  !
  subroutine ReadInput
    !
    use  input
    !
    integer(ik)  :: iobject(Nobjects)
    integer(ik)  :: ipot=0,iso=0,ncouples=0,il2=0,ilxly=0,iabi=0,idip=0,iss=0,isso=0,ibobrot=0,isr=0,idiab=0,iquad=0
    integer(ik)  :: Nparam,alloc,iparam,i,j,iobs,i_t,iref,jref,istate,jstate,istate_,jstate_,item_,ibraket,iabi_,iterm,iobj
    integer(ik)  :: Nparam_check    !number of parameters as determined automatically by duo (Nparam is specified in input).
    logical      :: zNparam_defined ! true if Nparam is in the input, false otherwise..
    integer(ik)  :: itau,lambda_,x_lz_y_,iobject_
    logical      :: integer_spin = .false., matchfound
    real(rk)     :: unit_field = 1.0_rk,unit_adjust = 1.0_rk, unit_r = 1.0_rk,spin_,jrot2,gns_a,gns_b
    real(rk)     :: f_t,jrot,j_list_(1:jlist_max)=-1.0_rk,omega_,sigma_,hstep = -1.0_rk
    !
    character(len=cl) :: w,ioname
    character(len=wl) :: large_fmt
    !
    integer(ik)       :: iut !  iut is a unit number. 
    !
    type(fieldT),pointer      :: field
    logical :: eof,include_state,allgrids
    logical :: symmetry_defined=.false.
    integer :: ic,ierr
    !
    ! -----------------------------------------------------------
    !
    !
    ! read the general input
    ! by Lorenzo Lodi
    ! read everything from std input and write to a temporary (scratch) file.
    !
    write(ioname, '(a, i4)') 'write to a temporary (scratch) file.'
    call IOstart(trim(ioname), iut)
    !
    open(unit=iut, status='scratch', action='readwrite')
    write(large_fmt, '(A,i0,A)') '(A', wl, ')'
    trans_loop: do i=1, max_input_lines
      read(unit=*,fmt=large_fmt,iostat=ierr) line_buffer
      if(ierr /=0) exit trans_loop
      write(iut, '(a)') trim(line_buffer)

      ! This is a hack; I need to know if to echo the input or not before processing it
      ! The option 'do_not_echo_input' is dealt with here as a special case
      line_buffer = adjustl(line_buffer) ! remove leading spaces

      do j=1, len(trim(line_buffer)) ! convert to uppercase
       ic = ichar( line_buffer(j:j))
       if( ic >= 97) line_buffer(j:j) = achar(ic-32)
      enddo

      if( trim(line_buffer) == 'DO_NOT_ECHO_INPUT' ) job%zEchoInput = .false.
    enddo trans_loop
    rewind(iut)

    !
    !
    ! default constants
    !
    jmin = 0 ; jmax = 0
    !
    ! To count objects
    iobject = 0
    !
    if( job%zEchoInput) then
      write(out,"('(Transcript of the input --->)')")
      call input_options(echo_lines=.true.,error_flag=1)
    else
      call input_options(echo_lines=.false.,error_flag=1)
    endif

    do
        zNparam_defined = .false. ! For each input block, set number of points/params to undefined
        call read_line(eof,iut) ; if (eof) exit
        call readu(w)
        select case(w)
          !
        case("STOP","FINISH","END")
          !
          exit
          !
        case("HYPERFINE")
          ! skip if hyperfine NONE
          if (Nitems>1) then
            call readu(w)
            if (trim(w)=="NONE".or.trim(w)=="OFF") then
              action%hyperfine = .false.
              do while (trim(w)/="".and.trim(w)/="END")
                call read_line(eof,iut) ; if (eof) exit
                call readu(w)
              enddo
              cycle
            endif
          endif
          !
          action%save_eigen_J = .true.
          action%hyperfine = .True.
          !
          job%basis_set = 'KEEP'
          !
          do while (trim(w)/="".and.trim(w)/="END")
            call read_line(eof,iut) ; if (eof) exit
            call readu(w)
            select case(w)
            case("I")
              call readf(F1_hyperfine_setup%I1)
            case default
              call report ("Unrecognized hyperfine keyword: "//trim(w),.true.)
            end select

            call read_line(eof,iut) ; if (eof) exit
            call readu(w)
         end do
         !
        case("PRINT_PECS_AND_COUPLINGS_TO_FILE")
          job%print_pecs_and_couplings_to_file = .true.
          !
        case("PRINT_VIBRATIONAL_ENERGIES_TO_FILE")
          job%print_vibrational_energies_to_file = .true.
          !
        case("PRINT_ROVIBRONIC_ENERGIES_TO_FILE")
          !
          job%print_rovibronic_energies_to_file = .true.
        case("DO_NOT_ECHO_INPUT") 
          !
          job%zEchoInput = .false.
        case("DO_NOT_SHIFT_PECS") 
          ! 
          job%zShiftPECsToZero = .false.
        case("DO_NOT_INCLUDE_JS_COUPLING") 
          ! 
          job%zExclude_JS_coupling = .true.
          !
        case("") ! do nothing in case of blank lines
          !
        case("ASSIGN_V_BY_COUNT")
          !
          job%assign_v_by_count = .true.
          !
        case("LEGACY","OLD-VERSION","VERSION")
          !
          job%legacy_version = .true.
          !
          if (item>1) then 
            call readi(item_)
            if (item_>2018) job%legacy_version = .false.
          endif 
          !
        case ("SOLUTIONMETHOD")
          !
          call readu(w)
          solution_method = trim(w)
          ! 
        case ("L2CONVENTION")
          !
          call readu(w)
          !
          select case(w)
             case ("SPECIFY_L^2","SPECIFY_L**2" ,"DEFAULT")
                job%diag_L2_fact = 1._rk  !default case
             case ("SPECIFY_LX^2_PLUS_LY^2", "SPECIFY_LX**2_PLUS_LY**2" )
                job%diag_L2_fact = 0._rk  !Lorenzo's choice
             case default
               call report ("Unrecognized L2CONVENTION "//trim(w)// ". " // &
                           "Implemented: DEFAULT, SPECIFY_L**2, SPECIFY_L^2, SPECIFY_LX^2_PLUS_LY^2" // &
                            ", SPECIFY_LX**2_PLUS_LY**2")
          end select
          !
        case ("MEM","MEMORY")
          !
          call readf(memory_limit)
          !
          if (nitems<2) call report("Please specify the memory units",.true.)
          !
          call readu(w)
          !
          select case(w)
              !
            case default 
              !
              call report("Unexpected argument in MEMORY",.true.)
              !
            case("TB","T")
              !
              memory_limit = memory_limit*1024.0_rk
              !
            case("GB","G")
              !
              memory_limit = memory_limit
              !
            case("MB","M")
              !
              memory_limit = memory_limit/1024.0_rk
              !
            case("KB","K")
              !
              memory_limit = memory_limit/1024.0_rk**2
              !
            case("B")
              !
              memory_limit = memory_limit/1024.0_rk**3
              !
          end select
          !
        case ("ATOMS") ! chemical symbols of the atoms (possibly with atomic numbers)
          !
          call reada(symbol1)
          call reada(symbol2)
          !
        case ("MASSES","MASS")
          !
          call readf(m1)
          call readf(m2)
          !
        case ("MOLECULE","MOL")
          !
          call readu(w)
          !
          select case(w)
          !
          case ("C2","ALO","X2","XY","CO","CAO","NIH","MGH")
            !
            job%molecule = trim(w)
            !
          case default
            !
            ! write (out,"('  I see this molecule for the first time.')")  
            !
            !call report ("Unrecognized unit name "//trim(w)//"implemented: (C2,ALO,X2, XY)",.true.)
            !
          end select
          !
        case('J_LIST','JLIST','JROT','J')
          !
          jmin =  1e6
          jmax = -1.0
          integer_spin = .false.
          i = 0
          do while (item<Nitems.and.i<jlist_max)
             !
             i = i + 1
             !
             call readu(w)
             !
             if (trim(w)/="-") then
               !
               read(w,*) jrot
               !
               j_list_(i) = jrot
               !
             else
               !
               call readf(jrot2)
               !
               do while (item<=Nitems.and.nint(2.0*jrot)<nint(2.0*jrot2))
                 !
                 jrot = jrot + 1.0
                 j_list_(i) = jrot
                 i = i + 1
                 !
               enddo
               i = i - 1
               !
             endif
             !
             if (i==1.and.mod(nint(2.0_rk*jrot+1.0_rk),2)==1) integer_spin = .true.
             !
             if (i>1.and.mod(nint(2.0_rk*jrot+1.0_rk),2)==0.and.integer_spin) then
               !
               call report("The multiplicities of J-s in J-list/Jrot are inconsistent",.true.)
               !
             endif
             !
             jmin = min(jmin,j_list_(i))
             jmax = max(jmax,j_list_(i))
             !
          enddo
          !
          job%nJ = i
          !
          allocate(job%j_list(i),stat=alloc)
          !
          job%J_list(1:i) = J_list_(1:i)
          !
        !case ("JROT")
        !  !
        !  call readf(jmin)
        !  !
        !  if (nitems>2) then
        !    !
        !    call readf(jmax)
        !    !
        !  else
        !    !
        !    jmax = jmin
        !    !
        !  endif
        !  !
        !  ! check the multiplicity
        !  !
        !  if (mod(nint(2.0_rk*jmin+1.0_rk),2)==1) integer_spin = .true.
        !  !
        !  if (mod(nint(2.0_rk*jmax+1.0_rk),2)==0.and.integer_spin) then
        !    !
        !    call report("The multiplicities of jmin and jmax are inconsistent",.true.)
        !    !
        !  endif
        !  !
        case ("NSTATES","NESTATES")
          !
          call readi(nestates)
          !
          if (nestates<1) call report("nestates cannot be 0 or negative",.true.)
          !
          ! the maximum number of couplings possible, assuming each state may be doubly degenerate
          !
          ncouples = 2*nestates*(2*nestates-1)/2
          !
          ! allocate all fields representing the hamiltonian: PEC, spin-orbit, L2, LxLy etc.
          !
          allocate(poten(nestates),spinorbit(ncouples),l2(ncouples),lxly(ncouples),spinspin(nestates),spinspino(nestates), &
                   bobrot(nestates),spinrot(nestates),job%vibmax(nestates),job%vibenermax(nestates),diabatic(ncouples),&
                   lambdaopq(nestates),lambdap2q(nestates),lambdaq(nestates),nac(nestates),quadrupoletm(ncouples),stat=alloc)

          do i = 1, GLOBAL_NUM_HFCC_OBJECT
            allocate(hfcc1(i)%field(nestates), stat=alloc)
          end do   
          !
          ! initializing the fields
          !
          job%vibmax = 1e8
          job%vibenermax = enermax
          !
          allocate(abinitio(nestates*Nobjects+4*ncouples),stat=alloc)
          !
        case ("NREFSTATES")
          !
          call readi(nrefstates)
          !
          if (nrefstates>nestates) call report("nrefstates cannot be larger than nestates",.true.)
          !
       case ("GRID")
         !
         call read_line(eof,iut) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
             !
           case ("NPOINTS")
             !
             call readi(grid%npoints)
             !
           case ("STEP")
             !
             call readf(hstep)
             !
           case ("NSUB","TYPE")
             !
             call readi(grid%nsub)
             !
           case ("RE","REF")
             !
             call readf(grid%re)
             !
           case ("ALPHA")
             !
             call readf(grid%alpha)
             !
           case("RANGE")
             !
             call readf(grid%rmin)
             call readf(grid%rmax)
             !
           case default
             !
             call report ("Unrecognized keyword in GRID: "//trim(w),.true.)
             !
           end select
           !
           if (hstep>0.and.grid%npoints/=0) then 
             write(out,"('Illegal grid-input: npoints and step should not appear together')")
             stop "Illegal grid-input: npoints and step should not appear together"
           endif
           !
           if (hstep>0) then 
             grid%npoints = (grid%rmax - grid%rmin)/hstep+1
           endif
           !
           call read_line(eof,iut) ; if (eof) exit
           call readu(w)
           !
         enddo
         !
         if (trim(w)/="".and.trim(w)/="END") then
            !
            call report ("Unrecognized keyword in GRID: "//trim(w),.true.)
            !
         endif
         !
       case ("CONTRACTION","VIBRATIONS","VIBRATIONALBASIS")
         !
         call read_line(eof,iut) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
             !
           case('VIB','ROT','OMEGA')
             !
             job%contraction = trim(w)
             !
           case ("VMAX","VIBMAX","NMAX")
             !
             vmax = 1
             !
             if (Nitems<2) job%vibmax = grid%npoints-1
             !
             do i = 1,min(Nitems-1,Nestates)
               call readi(job%vibmax(i))
               job%vibmax(i+1:Nestates) = job%vibmax(i)
               vmax = max(vmax,job%vibmax(i))
             enddo
             !
           case("ENERMAX")
             !
             do i = 1,min(Nitems-1,Nestates)
               call readf(job%vibenermax(i))
               job%vibenermax(i+1:Nestates) = job%vibenermax(i)
             enddo
             !
           case default
             !
             call report ("Unrecognized keyword in CONTRACTION: "//trim(w),.true.)
             !
           end select
           !
           call read_line(eof,iut) ; if (eof) exit
           call readu(w)
           !
         enddo
         !
         if (trim(w)/="".and.trim(w)/="END") then
            !
            write (out,"('input: wrong last line in CONTRACTION =',a)") trim(w)
            stop 'input - illegal last line in CONTRACTION'
            !
         endif
         !
       case("CHECK_POINT","CHECKPOINT","CHECKPOINTS")
         !
         job%eigenfile%vectors  = 'eigen'
         !
         call readu(w)
         !
         call read_line(eof,iut) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
           !
           case('EIGENFUNC','EIGENVECT','EIGENVECTORS')
             !
             call readu(w)
             !
             job%IO_eigen = trim(w)
             !
             if (all(trim(w)/=(/'READ','SAVE','NONE'/))) then 
               call report('ReadInput: illegal key in CHECK_POINT '//trim(w),.true.)
             endif 
             !
           case('VECTOR-FILENAME','VECTOR','FILENAME')
             !
             call reada(w)
             !
             job%eigenfile%vectors = trim(w)
             !
           case('BASIS_SET','BASIS-SET') 
             !
             call readu(w)
             !
             if (all(trim(w)/=(/'KEEP','NONE'/))) then 
               call report('ReadInput: illegal key in CHECK_POINT /= KEEP or NONE '//trim(w),.true.)
             endif 
             !
             job%basis_set  = trim(w)
             !
           case('TM','DIPOLE')
             !
             call readu(job%IO_dipole)
             !
             if (all(trim(w)/=(/'READ','SAVE','NONE'/))) then 
               call report('ReadInput: illegal key in CHECK_POINT '//trim(w),.true.)
             endif 
             !
           case default 
             !
             call report('ReadInput: illegal key in CHECK_POINT '//trim(w),.true.)
             !
           end select 
           !
           call read_line(eof,iut) ; if (eof) exit
           call readu(w)
           !
         enddo 
         !
         if (trim(w)/="".and.trim(w)/="END") then 
            call report('ReadInput: wrong last line in CHECK_POINTS ='//trim(w),.true.)
         endif 
         !
       case ("DIAGONALIZER","EIGENSOLVER","FINALSTATES")
         !
         call read_line(eof,iut) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
           !
           case('SYEV','SYEVR')
             !
             job%diagonalizer = trim(w)
             !
           case ("GAMMA")
             !
             i = 1
             job%select_gamma = .false.
             call readi(i_t)
             !
             do while (i_t/=0.and.i<=size(job%select_gamma))
                !
                i = i+1
                !
                job%select_gamma(i_t) = .true.
                !
                call readi(i_t)
                !
             enddo
             !
           case("NROOTS")
             !
             !call readi(job%nroots)
             !
             if (nitems-1==1) then
                !
                call readi(i)
                job%nroots(:) = i
                !
             else
               !
               if (nitems-1>20) then
                  !
                  write (out,"('input: too many entries in roots (>20): ',i8)") nitems-1
                  stop 'input - illigal number of entries in nroots'
                  !
               endif
               !
               do i =1,nitems-1
                  !
                  call readi(job%nroots(i))
                  !
               end do
               !
             endif
             !
           case("MAXITER")
             !
             call readi(job%maxiter)
             !
           case("CONTRACTION")
             !
             call readu(job%contraction)
             !
           case("TOLERANCE","TOL")
             !
             call readf(job%tolerance)
             !
           case("UPLIMIT","ENERMAX","ENERCUT")
             !
             call readf(job%upper_ener)
             !
           case("THRESHOLD","THRESH")
             !
             call readf(job%thresh)
             !
           case("ZPE")
             !
             call readf(job%zpe)
             job%shift_to_zpe = .false.
             !
           case default
             !
             call report ("Unrecognized keyword in DIAGONALIZER: "//trim(w),.true.)
             !
           end select
           !
           call read_line(eof,iut) ; if (eof) exit
           call readu(w)
           !
         enddo
         !
         if (trim(w)/="".and.trim(w)/="END") then
            !
            write (out,"('input: wrong last line in DIAGONALIZER =',a)") trim(w)
            stop 'input - illigal last line in DIAGONALIZER'
            !
         endif
         !
       case("FITTING")
         !
         action%fitting = .true.
         !
         ! skip if fitting NONE
         if (Nitems>1) then
            call readu(w)
            if (trim(w)=="NONE".or.trim(w)=="OFF") then
               action%fitting = .false.
               do while (trim(w)/="".and.trim(w)/="END")
                 call read_line(eof,iut) ; if (eof) exit
                 call readu(w)
               enddo
               cycle
            endif
         endif
         !
         call read_line(eof,iut) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
              !
           case('J_LIST','JLIST','J','JROT')
             !
             i = 0
             !do while (item<Nitems.and.i<jlist_max)
             !   !
             !   i = i + 1
             !   !
             !   call readf(j_list_(i))
             !   !
             !enddo
             !
             do while (item<Nitems.and.i<jlist_max)
                !
                i = i + 1
                !
                call readu(w)
                !
                if (trim(w)/="-") then
                  !
                  read(w,*) jrot
                  !
                  j_list_(i) = jrot
                  !
                else
                  !
                  call readf(jrot2)
                  !
                  do while (item<=Nitems.and.nint(2.0*jrot)<nint(2.0*jrot2))
                    !
                    jrot = jrot + 1.0
                    j_list_(i) = jrot
                    i = i + 1
                    !
                  enddo
                  !
                  i = i - 1
                  !
                endif
                !
             enddo
             !
             fitting%nJ = i
             allocate(fitting%j_list(i),stat=alloc)
             fitting%J_list(1:i) = J_list_(1:i)
             !
           case('ITMAX','ITERMAX','ITER')
             !
             call readi(fitting%itermax)
             !
           case('ROBUST')
             !
             call readf(fitting%robust)
             !
           case('TARGET_RMS')
             !
             call readf(fitting%target_rms)
             !
           case('FIT_TYPE')
             !
             call readu(fitting%fit_type)
             !
           case('THRESH_ASSIGN','THRESH_REASSIGN','THRESH_LOCK','LOCK','LOCK_QUANTA')
             !
             call readf(fitting%threshold_lock)
             !
           case('THRESH_OBS-CALC') 
             ! switch off weights for residuals larger than THRESH_OBS-CALC
             !
             call readf(fitting%threshold_obs_calc)
             !
           case("IPARAM")
             !
             call readi(fitting%iparam(1))
             call readi(fitting%iparam(2))
             !
           case('GEOMETRIES')
             !
             call readl(fitting%geom_file)
             !
           case('OUTPUT')
             !
             call reada(fitting%output_file)
             !
           case('ZPE')
             !
             call readf(fitting%zpe)
             fitting%shift_to_zpe = .false.
             !
           case('FIT_FACTOR')
             !
             call readf(fitting%factor)
             !
           case('FIT_SCALE')
             !
             call readf(fitting%fit_scaling)
             !
           case('LINEAR_SEARCH','LINEAR-SEARCH')
             !
             call readi(fitting%linear_search)
             !
           case('ABINIIO')
             !
             ! ignore the experiment and fit to the ab initio curves only
             !
             fitting%factor = small_
             !
           case('ENER','ENERGIES','FREQ','FREQUENCY','FREQUENCIES','WAVENUMBERS')
             !
             if (w(1:4)=='FREQ'.or.w(1:4)=='WAVE') then
               action%frequency = .true.
             endif
             !
             Nparam_check = 0
             !
             call input_options(echo_lines=.false.)
             !
             do while (trim(w)/="END")
                !
                call read_line(eof,iut) ; if (eof) exit
                !
                call readu(w)
                !
                Nparam_check = Nparam_check+1
                !
             enddo
             !
             Nparam_check = Nparam_check-1
             !
             if( job%zEchoInput) call input_options(echo_lines=.true.)
             !
             if (trim(w) /= "END") then
                 call report ("ERROR: Cannot find `END' statement)",.true.)
             endif
             !
             ! go back to beginning of VALUES block and reset `w' to original value
             do i=1, Nparam_check+1
               backspace(unit=iut)
             enddo
             !
             fitting%Nenergies = Nparam_check
             !
             !call readi(fitting%Nenergies)
             !
             allocate (fitting%obs(1:fitting%Nenergies),stat=alloc)
             if (alloc/=0) then
               write (out,"(' Error ',i0,' initializing obs. energy related arrays')") alloc
               stop 'obs. energy arrays - alloc'
             end if
             !
             iobs = 0
             !
             call read_line(eof,iut) ; if (eof) exit
             call readu(w)
             !
             do while (trim(w)/="END".and.iobs<fitting%Nenergies)
                !
                iobs = iobs + 1
                !
                if (.not.action%frequency.and.nitems<9) then
                   call report('input: wrong number of records in obs_fitting_energies (maybe using old input' &
                                                                            // ' and need to add omega)',.true.)
                endif
                !
                if (action%frequency.and.nitems<16) then
                   call report('input: wrong number of records in obs_fitting_frequency (maybe line is too long ' &
                                                           // ' (>300) or old input and need to add omega)',.true.)
                endif
                !
                CALL reread(-1)
                !
                call readf(fitting%obs(iobs)%Jrot)
                !
                !read(w,"(f9.1)") fitting%obs(iobs)%Jrot
                !
                ! get the index number of the current Jrot in J_list
                !
                i = 0
                matchfound = .false.
                if (.not.associated(fitting%J_list)) then 
                   fitting%nJ = 1
                   allocate(fitting%j_list(1),stat=alloc)
                   fitting%J_list(1) = 0
                   if (.not.integer_spin) fitting%J_list(1) = 0.5
                endif 
                do while( i<fitting%nJ.and..not.matchfound )
                  !
                  i = i + 1
                  if (fitting%J_list(i)/=fitting%obs(iobs)%Jrot) cycle
                  matchfound = .true.
                  !
                enddo
                !
                ! skip current line if these Jrot-s are not processed
                if (.not.matchfound) then
                  iobs = iobs-1
                  call read_line(eof,iut) ; if (eof) exit
                  call readu(w)
                  cycle
                endif
                !
                fitting%obs(iobs)%irot = i
                !
                ! parity:
                !
                call readu(w)
                !
                select case(w)
                  !
                case ('E')
                  !
                  if ( mod( nint( 2.0*fitting%obs(iobs)%Jrot ),2 )==1 ) then
                    itau = mod( nint( fitting%obs(iobs)%Jrot-0.5 ),2 )
                  else
                    itau = mod( nint( fitting%obs(iobs)%Jrot ),2 )
                  endif
                  !
                case ('F')
                  !
                  if ( mod( nint( 2.0*fitting%obs(iobs)%Jrot ),2 )==1 ) then
                    itau = mod( nint( fitting%obs(iobs)%Jrot-0.5 )+1,2 )
                  else
                    itau = mod( nint( fitting%obs(iobs)%Jrot )+1,2 )
                  endif
                  !
                case ('+','+1','1')
                  !
                  itau = 0
                  !
                case ('-','-1')
                  !
                  itau = 1
                  !
                case default
                  !
                  read(w,"(i2)") itau
                  !
                end select
                !
                fitting%obs(iobs)%iparity = itau
                !
                call readi(fitting%obs(iobs)%N)
                !
                if (action%frequency) then
                  call readf(fitting%obs(iobs)%Jrot_)
                  !
                  i = 0 
                  matchfound = .false.
                  do while( i<fitting%nJ.and..not.matchfound )
                    !
                    i = i + 1
                    if (fitting%J_list(i)/=fitting%obs(iobs)%Jrot_) cycle
                    matchfound = .true.
                    !
                  enddo
                  !
                  ! skip current line if these Jrot-s are not processed
                  if (.not.matchfound) then
                    iobs = iobs-1
                    call read_line(eof,iut) ; if (eof) exit
                    call readu(w)
                    cycle
                  endif
                  !
                  fitting%obs(iobs)%irot_ = i
                  !
                  ! parity:
                  !
                  call readu(w)
                  !
                  select case(w)
                    !
                  case ('E')
                    !
                    if ( mod( nint( 2.0*fitting%obs(iobs)%Jrot_ ),2 )==1 ) then
                      itau = mod( nint( fitting%obs(iobs)%Jrot_-0.5 ),2 )
                    else
                      itau = mod( nint( fitting%obs(iobs)%Jrot_ ),2 )
                    endif
                    !
                  case ('F')
                    !
                    if ( mod( nint( 2.0*fitting%obs(iobs)%Jrot_ ),2 )==1 ) then
                      itau = mod( nint( fitting%obs(iobs)%Jrot_-0.5 )+1,2 )
                    else
                      itau = mod( nint( fitting%obs(iobs)%Jrot_ )+1,2 )
                    endif
                    !
                  case ('+','+1','1')
                    !
                    itau = 0
                    !
                  case ('-','-1')
                    !
                    itau = 1
                    !
                  case default
                    !
                    read(w,"(i2)") itau
                    !
                  end select
                  !
                  fitting%obs(iobs)%iparity_ = itau
                  !
                  call readi(fitting%obs(iobs)%N_)
                  !
                  call readf(fitting%obs(iobs)%energy)
                  !
                else
                  call readf(fitting%obs(iobs)%energy)
                endif
                !
                call readi(fitting%obs(iobs)%quanta%istate)
                !
                ! skip current line if this state is not processed
                !
                if (fitting%obs(iobs)%quanta%istate>Nestates) then
                  iobs = iobs-1
                  call read_line(eof,iut) ; if (eof) exit
                  call readu(w)
                  cycle
                endif
                !
                call readi(fitting%obs(iobs)%quanta%v)
                call readi(fitting%obs(iobs)%quanta%ilambda)
                call readf(fitting%obs(iobs)%quanta%sigma)
                !
                if (.not.action%frequency.and.nitems==9) then
                  ! old input where omega was not present
                  !
                  fitting%obs(iobs)%quanta%omega = fitting%obs(iobs)%quanta%sigma + real(fitting%obs(iobs)%quanta%ilambda,rk)
                  !
                else
                  !
                  call readf(fitting%obs(iobs)%quanta%omega)
                  !
                endif
                !
                fitting%obs(iobs)%quanta%spin = poten(fitting%obs(iobs)%quanta%istate)%spini
                !
                if (action%frequency) then
                  !
                  call readi(fitting%obs(iobs)%quanta_%istate)
                  !
                  ! skip current line if this state is not processed
                  !
                  if (fitting%obs(iobs)%quanta_%istate>Nestates) then
                    iobs = iobs-1
                    call read_line(eof,iut) ; if (eof) exit
                    call readu(w)
                    cycle
                  endif
                  !
                  call readi(fitting%obs(iobs)%quanta_%v)
                  call readi(fitting%obs(iobs)%quanta_%ilambda)
                  call readf(fitting%obs(iobs)%quanta_%sigma)
                  call readf(fitting%obs(iobs)%quanta_%omega)
                  !
                  fitting%obs(iobs)%quanta_%spin = poten(fitting%obs(iobs)%quanta_%istate)%spini
                  !
                  !fitting%obs(iobs)%quanta_%omega = fitting%obs(iobs)%quanta_%sigma + real(fitting%obs(iobs)%quanta_%ilambda,rk)
                  !
                endif
                !
                if (fitting%obs(iobs)%quanta%istate>Nestates.or.fitting%obs(iobs)%quanta%istate<1) then
                   call report('input: illegal state',.true.)
                endif
                !
                !call readf(fitting%obs(i)%quanta%omega)
                !
                !if (abs(fitting%obs(i)%quanta%sigma)>fitting%obs(i)%quanta%spin) then
                !   call report('input: sigma is large than spin',.true.)
                !endif
                !
                call readf(fitting%obs(iobs)%weight)
                !
                !if (fitting%obs(iobs)%Jrot<jmin.or.fitting%obs(iobs)%Jrot>jmax.or.&
                !   (action%frequency.and.fitting%obs(iobs)%Jrot_<jmin.or.fitting%obs(iobs)%Jrot_>jmax) ) then
                !    iobs = iobs-1
                !endif
                !
                call read_line(eof,iut) ; if (eof) exit
                call readu(w)
                !
             enddo
             !
           case default
             !
             call report ("Unrecognized keyword name (error 01): "//trim(w),.true.)
             !
           end select
           !
           call read_line(eof,iut) ; if (eof) exit
           call readu(w)
           !
         enddo
         !
         fitting%Nenergies = iobs
         !
         if (trim(w)/="".and.trim(w)/="END") then
            call report ('wrong last line in FITTING ',.false.)
         endif
          !
       case("SPIN-ORBIT","SPIN-ORBIT-X","POTEN","POTENTIAL","L2","L**2","LXLY","LYLX","ABINITIO",&
            "LPLUS","L+","L_+","LX","DIPOLE","TM","DIPOLE-MOMENT","DIPOLE-X",&
            "SPIN-SPIN","SPIN-SPIN-O","BOBROT","BOB-ROT","SPIN-ROT","SPIN-ROTATION","DIABATIC","DIABAT",&
            "LAMBDA-OPQ","LAMBDA-P2Q","LAMBDA-Q","LAMBDAOPQ","LAMBDAP2Q","LAMBDAQ","NAC",&
            "QUADRUPOLE", &
            "HFCC-BF", "HFCC-A", "HFCC-C", "HFCC-D", "HFCC-CI", "HFCC-EQQ0", "HFCC-EQQ2") 
          !
          ibraket = 0
          !
          ! initializing units
          unit_field = 1 ; unit_r = 1
          !
          select case (w)
             !
          case("DIPOLE","TM","DIPOLE-MOMENT","DIPOLE-X")
             !
             if (idip==0) then 
                allocate(dipoletm(ncouples),stat=alloc)
             endif
             !
             idip = idip + 1
             !
             call readi(iref)
             call readi(jref)
             !
             include_state = .false.
             loop_istated : do istate=1,Nestates
               do jstate=1,Nestates
                 !
                 if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                   include_state = .true.
                   istate_ = istate
                   jstate_ = jstate
                   exit loop_istated
                 endif
                 !
               enddo
             enddo loop_istated
             !
             if (.not.include_state) then
                 !write(out,"('The interaction ',2i8,' is skipped')") iref,jref
                 idip = idip - 1
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 cycle
             endif
             !
             if (idip>ncouples) then
                 write(out,"(2a,i4,a,i6)") trim(w),": Number of dipoles = ",idip," exceeds the maximal allowed value",ncouples
                 call report ("Too many couplings given in the input for"//trim(w),.true.)
             endif
             !
             field => dipoletm(idip)
             !
             call set_field_refs(field,iref,jref,istate_,jstate_)
             !
             field%class = "DIPOLE"
             !
             if (trim(w)=='DIPOLE-X') then
               field%molpro = .true.
             endif
             !
          case("SPIN-ORBIT","SPIN-ORBIT-X")
             !
             iobject(2) = iobject(2) + 1
             !
             call readi(iref)
             call readi(jref)
             !
             include_state = .false.
             loop_istate : do istate=1,Nestates
               do jstate=1,Nestates
                 !
                 if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                   include_state = .true.
                   istate_ = istate
                   jstate_ = jstate
                   exit loop_istate
                 endif
                 !
               enddo
             enddo loop_istate
             !
             if (.not.include_state) then
                 !write(out,"('The interaction ',2i8,' is skipped')") iref,jref
                 iobject(2) = iobject(2) - 1
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 cycle
             endif
             !
             iso = iobject(2)
             !
             if (iso>ncouples) then
                 write(out, "(2a,i4,a,i6)") trim(w),": Number of couplings = ",iso," exceeds the maximal allowed value",ncouples
                 call report ("Too many couplings given in the input for"//trim(w),.true.)
             endif
             !
             field => spinorbit(iso)
             !
             !call set_field_refs(field,iref,jref,istate_,jstate_)
             !
             field%iref = iref
             field%jref = jref
             field%istate = istate_
             field%jstate = jstate_
             !
             if (action%fitting) call report ("SPIN-ORBIT cannot appear after FITTING",.true.)
             field%class = "SPINORBIT"
             !
             if (trim(w)=='SPIN-ORBIT-X') then
               field%molpro = .true.
             endif
             !
          case("LXLY","LYLX","L+","L_+","LX")
             !
             iobject(4) = iobject(4) + 1
             !
             call readi(iref)
             call readi(jref)
             !
             include_state = .false.
             loop_istatex : do istate=1,Nestates
               do jstate=1,Nestates
                 !
                 if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                   include_state = .true.
                   istate_ = istate
                   jstate_ = jstate
                   exit loop_istatex
                 endif
                 !
               enddo
             enddo loop_istatex
             !
             if (.not.include_state) then
                 !write(out,"('The interaction ',2i8,' is skipped')") iref,jref
                 iobject(4) = iobject(4) - 1
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 cycle
             endif
             !
             ilxly = iobject(4)
             !
             if (ilxly>ncouples) then
                 print "(2a,i4,a,i6)",trim(w),": Number of L+ couplings = ",ilxly," exceeds the maximal allowed value",ncouples
                 call report ("Too many L+ couplings given in the input for"//trim(w),.true.)
             endif
             !
             field => lxly(ilxly)
             !
             call set_field_refs(field,iref,jref,istate_,jstate_)
             field%class = "L+"
             !
             if (action%fitting) call report ("LXLY (L+) cannot appear after FITTING",.true.)
             !
             if (trim(w)=='LX') then
               field%molpro = .true.
             endif
             !
          case("POTEN","POTENTIAL")
             !
             iobject(1) = iobject(1) + 1
             !
             if (iobject(1)>nestates) then
                 print "(a,i4,a,i6)","The state # ",iobject(1)," is not included for the total number of states",nestates
                 !call report ("Too many potentials given in the input",.true.)
                 iobject(1) = iobject(1) - 1
                 !
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 !
                 !call read_line(eof,iut) ; if (eof) exit
                 !call readu(w)
                 cycle
             endif
             !
             ipot = iobject(1)
             !
             field => poten(iobject(1))
             field%istate = iobject(1)
             field%jstate = iobject(1)
             !
             call readi(field%iref)
             field%jref = field%iref
             field%class = "POTEN"
             !
             ! Check if it was defined before 
             do istate=1,iobject(1)-1
                if (field%iref==poten(istate)%iref) then
                  call report ("poten object is repeated",.true.)
                endif
             enddo
             !
             if (action%fitting) call report ("POTEN cannot appear after FITTING",.true.)
             !
          case("L2","L**2", "L^2")
             !
             iobject(3) = iobject(3) + 1
             !
             call readi(iref) ; jref = iref
             !
             ! for nondiagonal L2 terms
             if (nitems>2) call readi(jref)
             !
             ! find the corresponding potential
             !
             include_state = .false.
             loop_istate_l2 : do istate=1,Nestates
               do jstate=1,Nestates
                 if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                   include_state = .true.
                   istate_ = istate
                   jstate_ = jstate
                   exit loop_istate_l2
                 endif
               enddo
             enddo loop_istate_l2
             !
             ! Check if it was defined before 
             do istate=1,iobject(3)-1
                if (iref==l2(istate)%iref.and.jref==l2(istate)%jref) then
                  call report ("L2 object is repeated",.true.)
                endif
             enddo
             !
             if (.not.include_state) then
                 !write(out,"('The L2 term ',2i8,' is skipped')") iref,jref
                 iobject(3) = iobject(3) - 1
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 cycle
             endif
             !
             il2 = iobject(3)
             !
             field => l2(il2)
             field%class = trim(classnames(3))
             !
             call set_field_refs(field,iref,jref,istate_,jstate_)
             !
             if (action%fitting) call report ("L2 cannot appear after FITTING",.true.)
             !
             !
          case("BOB-ROT","BOBROT")
             !
             iobject(7) = iobject(7) + 1
             !
             call readi(iref) ; jref = iref
             !
             ! find the corresponding potential
             !
             include_state = .false.
             loop_istate_bobrot : do istate=1,Nestates
                 if (iref==poten(istate)%iref) then
                   include_state = .true.
                   istate_ = istate
                   exit loop_istate_bobrot
                 endif
             enddo loop_istate_bobrot
             !
             ! Check if it was defined before 
             do istate=1,iobject(7)-1
                if (iref==bobrot(istate)%iref.and.jref==bobrot(istate)%jref) then
                  call report ("BROT object is repeated",.true.)
                endif
             enddo
             !
             if (.not.include_state) then
                 !write(out,"('The BOB-ROT term ',1i8,' is skipped')") iref
                 iobject(7) = iobject(7) - 1
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 cycle
             endif
             !
             ibobrot = iobject(7)
             !
             field => bobrot(ibobrot)
             !
             call set_field_refs(field,iref,jref,istate_,istate_)
             !
             field%class = trim(classnames(7))
             !
             if (action%fitting) call report ("BOBrot cannot appear after FITTING",.true.)
             !
          case("SPIN-SPIN")
             !
             iobject(5) = iobject(5) + 1
             !
             call readi(iref) ; jref = iref
             !
             if (nitems>2) call readi(jref)
             !
             ! find the corresponding potential
             !
             include_state = .false.
             loop_istate_ss : do istate=1,Nestates
               do jstate=1,Nestates
                 if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                   include_state = .true.
                   istate_ = istate
                   jstate_ = jstate
                   exit loop_istate_ss
                 endif
               enddo
             enddo loop_istate_ss
             !
             ! Check if it was defined before 
             do istate=1,iobject(5)-1
                if (iref==spinspin(istate)%iref.and.jref==spinspin(istate)%jref) then
                  call report ("Spin-spin object is repeated",.true.)
                endif
             enddo
             !
             if (.not.include_state) then
                 !write(out,"('The SS term ',2i8,' is skipped')") iref,jref
                 iobject(5) = iobject(5) - 1
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 cycle
             endif
             !
             iss = iobject(5)
             !
             field => spinspin(iss)
             !
             call set_field_refs(field,iref,jref,istate_,jstate_)
             !
             field%class = trim(classnames(5))
             !
             if (action%fitting) call report ("Spin-spin cannot appear after FITTING",.true.)
             !
             ! non-diagonal spin-spin term 
             !
          case("SPIN-SPIN-O")
             !
             iobject(6) = iobject(6) + 1
             !
             call readi(iref) ; jref = iref
             !
             if (nitems>2) call readi(jref)
             !
             ! find the corresponding potential
             !
             include_state = .false.
             loop_istate_sso : do istate=1,Nestates
               do jstate=1,Nestates
                 if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                   include_state = .true.
                   istate_ = istate
                   jstate_ = jstate
                   exit loop_istate_sso
                 endif
               enddo
             enddo loop_istate_sso
             !
             ! Check if it was defined before 
             do istate=1,iobject(6)-1
                if (iref==spinspino(istate)%iref.and.jref==spinspino(istate)%jref) then
                  call report ("SS-o object is repeated",.true.)
                endif
             enddo
             !
             if (.not.include_state) then
                 !write(out,"('The SS-o term ',2i8,' is skipped')") iref,jref
                 iobject(6) = iobject(6) - 1
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 cycle
             endif
             !
             isso = iobject(6)
             !
             field => spinspino(isso)
             !
             call set_field_refs(field,iref,jref,istate_,jstate_)
             field%class = trim(classnames(6))
             !
             if (action%fitting) call report ("Spin-spin-o cannot appear after FITTING",.true.)
             !
          case("SPIN-ROT","SPIN-ROTATION")
             !
             ! spin-rotation (gammma) term 
             !
             iobject(8) = iobject(8) + 1
             !
             call readi(iref) ; jref = iref
             !
             ! find the corresponding potential
             !
             include_state = .false.
             loop_istate_sr : do istate=1,Nestates
               do jstate=1,Nestates
                 if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                   include_state = .true.
                   istate_ = istate
                   jstate_ = jstate
                   exit loop_istate_sr
                 endif
               enddo
             enddo loop_istate_sr
             !
             ! Check if it was defined before 
             do istate=1,iobject(8)-1
                if (iref==spinrot(istate)%iref.and.jref==spinrot(istate)%jref) then
                  call report ("SR object is repeated",.true.)
                endif
             enddo
             !
             if (.not.include_state) then
                 !write(out,"('The SR term ',2i8,' is skipped')") iref,jref
                 iobject(8) = iobject(8) - 1
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 cycle
             endif
             !
             isr = iobject(8)
             !
             field => spinrot(isr)
             !
             call set_field_refs(field,iref,jref,istate_,jstate_)
             !
             field%class = trim(classnames(8))
             !
             if (action%fitting) call report ("Spin-rot cannot appear after FITTING",.true.)
             !
          case("DIABAT","DIABATIC")
             !
             iobject(9) = iobject(9) + 1
             !
             call readi(iref) ; jref = iref
             !
             ! for nondiagonal terms
             if (nitems>2) call readi(jref)
             !
             ! find the corresponding potential
             !
             include_state = .false.
             loop_istate_diab : do istate=1,Nestates
               do jstate=1,Nestates
                 if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                   include_state = .true.
                   istate_ = istate
                   jstate_ = jstate
                   exit loop_istate_diab
                 endif
               enddo
             enddo loop_istate_diab
             !
             ! Check if it was defined before 
             do istate=1,iobject(9)-1
                if (iref==diabatic(istate)%iref.and.jref==diabatic(istate)%jref) then
                  call report ("diabatic object is repeated",.true.)
                endif
             enddo
             !
             if (.not.include_state) then
                 !write(out,"('The DIABATIC term ',2i8,' is skipped')") iref,jref
                 iobject(9) = iobject(9) - 1
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 cycle
             endif
             !
             idiab = iobject(9)
             !
             field => diabatic(idiab)
             !
             call set_field_refs(field,iref,jref,istate_,jstate_)
             field%class = trim(classnames(9))
             !
             if (action%fitting) call report ("DIABATIC cannot appear after FITTING",.true.)
             !
          case("LAMBDA-OPQ","LAMBDAOPQ")  ! o+p+q
             !
             iobject(10) = iobject(10) + 1
             !
             call readi(iref) ; jref = iref
             !
             ! for nondiagonal terms
             if (nitems>2) call readi(jref)
             !
             ! find the corresponding potential
             !
             include_state = .false.
             loop_istate_10 : do istate=1,Nestates
               do jstate=1,Nestates
                 if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                   include_state = .true.
                   istate_ = istate
                   jstate_ = jstate
                   exit loop_istate_10
                 endif
               enddo
             enddo loop_istate_10
             !
             ! Check if it was defined before 
             do istate=1,iobject(10)-1
                if (iref==lambdaopq(istate)%iref.and.jref==lambdaopq(istate)%jref) then
                  call report ("lambdaopq object is repeated",.true.)
                endif
             enddo
             !
             if (.not.include_state) then
                 !write(out,"('The LAMBDA-O term ',2i8,' is skipped')") iref,jref
                 iobject(10) = iobject(10) - 1
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 cycle
             endif
             !
             field => lambdaopq(iobject(10))
             !
             call set_field_refs(field,iref,jref,istate_,jstate_)
             !
             field%class = trim(CLASSNAMES(10))
             !
             if (action%fitting) call report (trim(field%class)//" cannot appear after FITTING",.true.)
             !
             ! -(p+2q)
             !
          case("LAMBDA-P2Q","LAMBDAP2Q")
             !
             iobject(11) = iobject(11) + 1
             !
             call readi(iref) ; jref = iref
             !
             ! for nondiagonal terms
             if (nitems>2) call readi(jref)
             !
             ! find the corresponding potential
             !
             include_state = .false.
             loop_istate_11 : do istate=1,Nestates
               do jstate=1,Nestates
                 if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                   include_state = .true.
                   istate_ = istate
                   jstate_ = jstate
                   exit loop_istate_11
                 endif
               enddo
             enddo loop_istate_11
             !
             ! Check if it was defined before 
             do istate=1,iobject(11)-1
                if (iref==lambdap2q(istate)%iref.and.jref==lambdap2q(istate)%jref) then
                  call report ("lambdap2q object is repeated",.true.)
                endif
             enddo
             !
             if (.not.include_state) then
                 !write(out,"('The LAMBDA-P term ',2i8,' is skipped')") iref,jref
                 iobject(11) = iobject(11) - 1
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 cycle
             endif
             !
             field => lambdap2q(iobject(11))
             !
             call set_field_refs(field,iref,jref,istate_,jstate_)
             !
             field%class = trim(CLASSNAMES(11))
             !
             if (action%fitting) call report (trim(field%class)//" cannot appear after FITTING",.true.)
             !
          case("LAMBDA-Q","LAMBDAQ")
             !
             iobject(12) = iobject(12) + 1
             !
             call readi(iref) ; jref = iref
             !
             ! for nondiagonal terms
             if (nitems>2) call readi(jref)
             !
             ! find the corresponding potential
             !
             include_state = .false.
             loop_istate_12 : do istate=1,Nestates
               do jstate=1,Nestates
                 if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                   include_state = .true.
                   istate_ = istate
                   jstate_ = jstate
                   exit loop_istate_12
                 endif
               enddo
             enddo loop_istate_12
             !
             ! Check if it was defined before 
             do istate=1,iobject(12)-1
                if (iref==lambdaq(istate)%iref.and.jref==lambdaq(istate)%jref) then
                  call report ("lambdaq object is repeated",.true.)
                endif
             enddo
             !
             if (.not.include_state) then
                 !write(out,"('The LAMBDA-Q term ',2i8,' is skipped')") iref,jref
                 iobject(12) = iobject(12) - 1
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 cycle
             endif
             !
             field => lambdaq(iobject(12))
             !
             call set_field_refs(field,iref,jref,istate_,jstate_)
             !
             field%class = trim(CLASSNAMES(12))
             !
             if (action%fitting) call report (trim(field%class)//" cannot appear after FITTING",.true.)
             !
          case("NAC")
             !
             call input_non_diagonal_field(Nobjects,13,iobject(13),nac,ierr)
             !
             field => nac(iobject(13))
             !
             if (ierr>0) cycle
             !
          case("QUADRUPOLE")
             !
             if (iquad==0) then 
                allocate(quadrupoletm(ncouples),stat=alloc)
             endif
             !
             iquad = iquad + 1
             !
             call readi(iref)
             call readi(jref)
             !
             include_state = .false.
             loop_quad : do istate=1,Nestates
               do jstate=1,Nestates
                 !
                 if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                   include_state = .true.
                   istate_ = istate
                   jstate_ = jstate
                   exit loop_quad
                 endif
                 !
               enddo
             enddo loop_quad
             !
             if (.not.include_state) then
                 iquad = iquad - 1
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 cycle
             endif
             !
             if (iquad>ncouples) then
                 write(out,"(2a,i4,a,i6)") trim(w),": Number of quadrupoles = ",iquad," exceeds the maximal allowed value",ncouples
                 call report ("Too many couplings given in the input for"//trim(w),.true.)
             endif
             !
             field => quadrupoletm(iquad)
             !
             call set_field_refs(field,iref,jref,istate_,jstate_)
             !
             field%class = "QUADRUPOLE"
             !
          case("HFCC-BF")
            hfcc1(1)%num_field = hfcc1(1)%num_field + 1
            iobject(21) = iobject(21) + 1
            call readi(iref); jref = iref
            include_state = .false.
            !
            ! find the corresponding potential
            include_state = .false.
            loop_istate_hfcc_bf : do istate=1, Nestates
              do jstate=1, Nestates
                if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                  include_state = .true.
                  istate_ = istate
                  jstate_ = jstate
                  exit loop_istate_hfcc_bf
                endif
              enddo
            enddo loop_istate_hfcc_bf
            !
            ! Check if it was defined before 
            do istate=1, hfcc1(1)%num_field - 1
              if (iref==hfcc1(1)%field(istate)%iref.and.jref==hfcc1(1)%field(istate)%jref) then
                call report ("Fermi contact object is repeated",.true.)
              endif
            enddo
            !
            if (.not.include_state) then
              !write(out,"('The Fermi-contact term ',2i8,' is skipped')") iref,jref
              hfcc1(1)%num_field = hfcc1(1)%num_field - 1
              do while (trim(w)/="".and.trim(w)/="END")
                call read_line(eof,iut) ; if (eof) exit
                call readu(w)
              enddo
              cycle
            endif
            !
            field => hfcc1(1)%field(hfcc1(1)%num_field) 
            !
            call set_field_refs(field, iref, jref, istate_, jstate_)
            field%class = "HFCC-BF-1"
            if (action%fitting) call report ("Fermi-contact object cannot appear after FITTING",.true.)
            !
          case("HFCC-A")
            hfcc1(2)%num_field = hfcc1(2)%num_field + 1
            iobject(22) = iobject(22) + 1
            call readi(iref); jref = iref
            include_state = .false.
            !
            ! find the corresponding potential
            include_state = .false.
            loop_istate_hfcc_a : do istate=1, Nestates
              do jstate=1, Nestates
                if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                  include_state = .true.
                  istate_ = istate
                  jstate_ = jstate
                  exit loop_istate_hfcc_a
                endif
              enddo
            enddo loop_istate_hfcc_a
            !
            ! Check if it was defined before 
            do istate=1, hfcc1(2)%num_field - 1
              if (iref==hfcc1(2)%field(istate)%iref.and.jref==hfcc1(2)%field(istate)%jref) then
                call report ("Nuclear spin -- orbit object is repeated",.true.)
              endif
            enddo
            !
            if (.not.include_state) then
              !write(out,"('The Nuclear spin -- orbit ',2i8,' is skipped')") iref,jref
              hfcc1(2)%num_field = hfcc1(2)%num_field - 1
              do while (trim(w)/="".and.trim(w)/="END")
                call read_line(eof,iut) ; if (eof) exit
                call readu(w)
              enddo
              cycle
            endif
            !
            field => hfcc1(2)%field(hfcc1(2)%num_field)
            !
            call set_field_refs(field, iref, jref, istate_, jstate_)
            field%class = "HFCC-A-1"
            if (action%fitting) call report ("Nuclear spin -- orbit object cannot appear after FITTING",.true.)
            !
          case("HFCC-C")
            hfcc1(3)%num_field = hfcc1(3)%num_field + 1
            iobject(23) = iobject(23) + 1
            call readi(iref); jref = iref
            include_state = .false.
            !
            ! find the corresponding potential
            include_state = .false.
            loop_istate_hfcc_c : do istate=1, Nestates
              do jstate=1, Nestates
                if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                  include_state = .true.
                  istate_ = istate
                  jstate_ = jstate
                  exit loop_istate_hfcc_c
                endif
              enddo
            enddo loop_istate_hfcc_c
            !
            ! Check if it was defined before 
            do istate=1, hfcc1(3)%num_field - 1
              if (iref==hfcc1(3)%field(istate)%iref.and.jref==hfcc1(3)%field(istate)%jref) then
                call report ("Diagonal nuclear spin - electron spin dipole-dipole object is repeated",.true.)
              endif
            enddo
            !
            if (.not.include_state) then
              !write(out,"('Diagonal nuclear spin - electron spin dipole-dipole term ',2i8,' is skipped')") iref,jref
              hfcc1(3)%num_field = hfcc1(3)%num_field - 1
              do while (trim(w)/="".and.trim(w)/="END")
                call read_line(eof,iut) ; if (eof) exit
                call readu(w)
              enddo
              cycle
            endif
            !
            field => hfcc1(3)%field(hfcc1(3)%num_field)
            !
            call set_field_refs(field, iref, jref, istate_, jstate_)
            field%class = "HFCC-C-1"
            if (action%fitting) then
              call report ("Diagonal nuclear spin - electron spin dipole-dipole object cannot appear after FITTING",.true.)
            endif
            !
          case("HFCC-D")
            hfcc1(4)%num_field = hfcc1(4)%num_field + 1
            iobject(24) = iobject(24) + 1
            call readi(iref); jref = iref
            include_state = .false.
            !
            ! find the corresponding potential
            include_state = .false.
            loop_istate_hfcc_d : do istate=1, Nestates
              do jstate=1, Nestates
                if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                  include_state = .true.
                  istate_ = istate
                  jstate_ = jstate
                  exit loop_istate_hfcc_d
                endif
              enddo
            enddo loop_istate_hfcc_d
            !
            ! Check if it was defined before 
            do istate=1, hfcc1(4)%num_field - 1
              if (iref==hfcc1(4)%field(istate)%iref.and.jref==hfcc1(4)%field(istate)%jref) then
                call report ("Off-diagonal nuclear spin - electron spin dipole-dipole object  is repeated",.true.)
              endif
            enddo
            !
            if (.not.include_state) then
              !write(out,"('The Off-diagonal nuclear spin - electron spin dipole-dipole object ',2i8,' is skipped')") iref,jref
              hfcc1(4)%num_field = hfcc1(4)%num_field - 1
              do while (trim(w)/="".and.trim(w)/="END")
                call read_line(eof,iut) ; if (eof) exit
                call readu(w)
              enddo
              cycle
            endif
            !
            field => hfcc1(4)%field(hfcc1(4)%num_field)
            !
            call set_field_refs(field, iref, jref, istate_, jstate_)
            field%class = "HFCC-D-1"
            if (action%fitting) then
              call report ("Off-diagonal nuclear spin - electron spin dipole-dipole object cannot appear after FITTING",.true.)
            endif
          case("HFCC-CI")
            !
            hfcc1(5)%num_field = hfcc1(5)%num_field + 1
            iobject(25) = iobject(25) + 1
            call readi(iref); jref = iref
            include_state = .false.
            !
            ! find the corresponding potential
            include_state = .false.
            loop_istate_hfcc_ci : do istate=1, Nestates
              do jstate=1, Nestates
                if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                  include_state = .true.
                  istate_ = istate
                  jstate_ = jstate
                  exit loop_istate_hfcc_ci
                endif
              enddo
            enddo loop_istate_hfcc_ci
            !
            ! Check if it was defined before 
            do istate=1, hfcc1(5)%num_field - 1
              if (iref==hfcc1(5)%field(istate)%iref.and.jref==hfcc1(5)%field(istate)%jref) then
                call report ("Nuclear spin -- rotation object is repeated",.true.)
              endif
            enddo
            !
            if (.not.include_state) then
              !write(out,"('The Fermi-contact term ',2i8,' is skipped')") iref,jref
              hfcc1(5)%num_field = hfcc1(5)%num_field - 1
              do while (trim(w)/="".and.trim(w)/="END")
                call read_line(eof,iut) ; if (eof) exit
                call readu(w)
              enddo
              cycle
            endif
            !
            field => hfcc1(5)%field(hfcc1(5)%num_field)
            !
            call set_field_refs(field, iref, jref, istate_, jstate_)
            field%class = "HFCC-CI-1"
            if (action%fitting) call report ("Nuclear spin -- rotation object cannot appear after FITTING",.true.)
            !
          case("HFCC-EQQ0")
            hfcc1(6)%num_field = hfcc1(6)%num_field + 1
            iobject(26) = iobject(26) + 1
            call readi(iref); jref = iref
            include_state = .false.
            !
            ! find the corresponding potential
            include_state = .false.
            loop_istate_hfcc_eqq0 : do istate=1, Nestates
              do jstate=1, Nestates
                if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                  include_state = .true.
                  istate_ = istate
                  jstate_ = jstate
                  exit loop_istate_hfcc_eqq0
                endif
              enddo
            enddo loop_istate_hfcc_eqq0
            !
            ! Check if it was defined before 
            do istate=1, hfcc1(6)%num_field - 1
              if (iref==hfcc1(6)%field(istate)%iref.and.jref==hfcc1(6)%field(istate)%jref) then
                call report ("Diagonal nuclear electric quadrupole object is repeated",.true.)
              endif
            enddo
            !
            if (.not.include_state) then
              !write(out,"('Diagonal nuclear electric quadrupole term ',2i8,' is skipped')") iref,jref
              hfcc1(6)%num_field = hfcc1(6)%num_field - 1
              do while (trim(w)/="".and.trim(w)/="END")
                call read_line(eof,iut) ; if (eof) exit
                call readu(w)
              enddo
              cycle
            endif
            !
            field => hfcc1(6)%field(hfcc1(6)%num_field)
            !
            call set_field_refs(field, iref, jref, istate_, jstate_)
            field%class = "HFCC-EQQ0-1"
            if (action%fitting) call report ("Diagonal nuclear electric quadrupole object cannot appear after FITTING",.true.)
            !
          case("HFCC-EQQ2")
            hfcc1(7)%num_field = hfcc1(7)%num_field + 1
            iobject(27) = iobject(27) + 1
            call readi(iref); jref = iref
            include_state = .false.
            !
            ! find the corresponding potential
            include_state = .false.
            loop_istate_hfcc_eqq2 : do istate=1, Nestates
              do jstate=1, Nestates
                if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                  include_state = .true.
                  istate_ = istate
                  jstate_ = jstate
                  exit loop_istate_hfcc_eqq2
                endif
              enddo
            enddo loop_istate_hfcc_eqq2
            !
            ! Check if it was defined before 
            do istate=1, hfcc1(7)%num_field - 1
              if (iref==hfcc1(7)%field(istate)%iref.and.jref==hfcc1(7)%field(istate)%jref) then
                call report ("Off-diagonal nuclear electric quadrupole object is repeated",.true.)
              endif
            enddo
            !
            if (.not.include_state) then
              !write(out,"('Off-diagonal nuclear electric quadrupole term ',2i8,' is skipped')") iref,jref
              hfcc1(7)%num_field = hfcc1(7)%num_field - 1
              do while (trim(w)/="".and.trim(w)/="END")
                call read_line(eof,iut) ; if (eof) exit
                call readu(w)
              enddo
              cycle
            endif
            !
            field => hfcc1(7)%field(hfcc1(7)%num_field)
            !
            call set_field_refs(field, iref, jref, istate_, jstate_)
            field%class = "HFCC-EQQ2-1"
            if (action%fitting) call report ("Off-diagonal nuclear electric quadrupole cannot appear after FITTING",.true.)
             !
          case("ABINITIO")
             !
             iabi = iabi + 1
             !
             call readu(w)
             !
             jref = 0
             jstate_ = 0
             iabi_ = 0
             !
             select case (w)
               !
             case ("POTEN","POTENTIAL")
               !
               ! find the corresponding potential
               !
               call readi(iref)
               !
               include_state = .false.
               loop_istate_abpot : do istate=1,Nestates
                   if (iref==poten(istate)%iref) then
                     include_state = .true.
                     !istate_ = istate
                     iabi_ = istate
                     exit loop_istate_abpot
                   endif
               enddo loop_istate_abpot
               !
             case("L2","L**2")
               !
               ! find the corresponding L2
               !
               call readi(iref) ; jref = iref
               if (nitems>2) call readi(jref)
               !
               include_state = .false.
               loop_istate_abl2 : do i=1,NL2
                   if (iref==l2(i)%iref.and.jref==l2(i)%jref) then
                     include_state = .true.
                     !istate_ = istate
                     !
                     iabi_ = Nestates + iso + i
                     !
                     exit loop_istate_abl2
                   endif
               enddo loop_istate_abl2

             case("BOB-ROT","BOBROT")
               !
               ! find the corresponding BB
               !
               call readi(iref) ; jref = iref
               !
               include_state = .false.
               loop_istate_abbobr : do i=1,NL2
                   if (iref==bobrot(i)%iref) then
                     include_state = .true.
                     !istate_ = istate
                     !
                     iabi_ = Nestates + iso + il2 + ilxly + iss+isso+i
                     !
                     exit loop_istate_abbobr
                   endif
               enddo loop_istate_abbobr
               !
             case("SPIN-ORBIT","SPIN-ORBIT-X")
               !
               call readi(iref)
               call readi(jref)
               !
               include_state = .false.
               loop_istate_abiso : do i=1,iso
                   !
                   if (iref==spinorbit(i)%iref.and.jref==spinorbit(i)%jref) then
                     include_state = .true.
                     !
                     iabi_ = Nestates + i
                     !
                     if (trim(w)=='SPIN-ORBIT-X') then
                       abinitio(iabi_)%molpro = .true.
                     endif
                     !
                     exit loop_istate_abiso
                   endif
                   !
               enddo loop_istate_abiso
               !
               w = "SPINORBIT"
               !
             case("LXLY","LYLX","L+","L_+","LX")
               !
               call readi(iref)
               call readi(jref)
               !
               include_state = .false.
               loop_istatex_abi : do i=1,ilxly
                   !
                   if (iref==lxly(i)%iref.and.jref==lxly(i)%jref) then
                     include_state = .true.
                     !
                     iabi_ = Nestates + iso + il2 + i
                     !
                     exit loop_istatex_abi
                   endif
                   !
               enddo loop_istatex_abi
               !
               if (trim(w)=='LX'.and.iabi_>0) then
                 abinitio(iabi_)%molpro = .true.
               endif
               !
               w = "LX"
               !
             case("SPIN-SPIN")
               !
               ! find the corresponding SS
               !
               call readi(iref) ; jref = iref
               if (nitems>2) call readi(jref)
               !
               include_state = .false.
               loop_istate_abss : do i=1,iss
                   if (iref==spinspin(i)%iref.and.jref==spinspin(i)%jref) then
                     include_state = .true.
                     !istate_ = istate
                     !
                     iabi_ = Nestates + iso + il2 + ilxly + i
                     !
                     exit loop_istate_abss
                   endif
               enddo loop_istate_abss
               !
             case("SPIN-SPIN-O")
               !
               ! find the corresponding SSO
               !
               call readi(iref) ; jref = iref
               if (nitems>2) call readi(jref)
               !
               include_state = .false.
               loop_istate_absso : do i=1,isso
                   if (iref==spinspino(i)%iref.and.jref==spinspino(i)%jref) then
                     include_state = .true.
                     !istate_ = istate
                     !
                     iabi_ = Nestates + iso + il2 + ilxly + iss+i
                     !
                     exit loop_istate_absso
                   endif
               enddo loop_istate_absso
               !
             case("SPIN-ROT","SPIN-ROTATION")
               !
               ! find the corresponding object
               !
               call readi(iref) ; jref = iref
               if (nitems>2) call readi(jref)
               !
               include_state = .false.
               loop_istate_absr : do i=1,isr
                   if (iref==spinrot(i)%iref.and.jref==spinrot(i)%jref) then
                     include_state = .true.
                     !istate_ = istate
                     !
                     iabi_ = Nestates + iso + il2 + ilxly + iss + isso + ibobrot + i
                     !
                     exit loop_istate_absr
                   endif
               enddo loop_istate_absr
               !
             case("DIABATIC","DIABAT")
               !
               ! find the corresponding object
               !
               call readi(iref) ; jref = iref
               if (nitems>2) call readi(jref)
               !
               include_state = .false.
               loop_istate_abdia : do i=1,idiab
                   if (iref==diabatic(i)%iref.and.jref==diabatic(i)%jref) then
                     include_state = .true.
                     !istate_ = istate
                     !
                     iabi_ = Nestates + iso + il2 + ilxly + iss + isso + ibobrot + isr + i
                     !
                     exit loop_istate_abdia
                   endif
               enddo loop_istate_abdia
               !
             case("LAMBDA-OPQ","LAMBDAOPQ")
               !
               ! find the corresponding object
               !
               call readi(iref) ; jref = iref
               if (nitems>2) call readi(jref)
               !
               iobject_ = 10
               !
               include_state = .false.
               loop_istate_ab10 : do i=1,iobject(iobject_)
                   if (iref==lambdaopq(i)%iref.and.jref==lambdaopq(i)%jref) then
                     include_state = .true.
                     !
                     iabi_ = sum(iobject(1:iobject_-1)) + i
                     !
                     exit loop_istate_ab10
                   endif
               enddo loop_istate_ab10
               !
             case("LAMBDA-P2Q","LAMBDAP2Q")
               !
               ! find the corresponding object
               !
               call readi(iref) ; jref = iref
               if (nitems>2) call readi(jref)
               !
               iobject_ = 11
               !
               include_state = .false.
               loop_istate_ab11 : do i=1,iobject(iobject_)
                   if (iref==lambdap2q(i)%iref.and.jref==lambdap2q(i)%jref) then
                     include_state = .true.
                     !
                     iabi_ = sum(iobject(1:iobject_-1)) + i
                     !
                     exit loop_istate_ab11
                   endif
               enddo loop_istate_ab11
               !
             case("LAMBDA-Q","LAMBDAQ")
               !
               ! find the corresponding object
               !
               call readi(iref) ; jref = iref
               if (nitems>2) call readi(jref)
               !
               iobject_ = 12
               !
               include_state = .false.
               loop_istate_ab12 : do i=1,iobject(iobject_)
                   if (iref==lambdaq(i)%iref.and.jref==lambdaq(i)%jref) then
                     include_state = .true.
                     !
                     iabi_ = sum(iobject(1:iobject_-1)) + i
                     !
                     exit loop_istate_ab12
                   endif
               enddo loop_istate_ab12
               !
             case("DIPOLE","DIPOLE-X")
               !
               call readi(iref)
               call readi(jref)
               !
               iobject_ = Nobjects
               !
               include_state = .false.
               loop_istate_abdip : do i=1,idip
                   !
                   if (iref==dipoletm(i)%iref.and.jref==dipoletm(i)%jref) then
                     include_state = .true.
                     !
                     iabi_ = sum(iobject(1:Nobjects-2)) + i
                     !
                     exit loop_istate_abdip
                   endif
                   !
               enddo loop_istate_abdip
               !
               if (trim(w)=='DIPOLE-X') then
                 abinitio(iabi_)%molpro = .true.
               endif
               !
             case("QUADRUPOLE")
               !
               call report ("Ab initio field is crrently not working with QUADRUPOLE",.true.)
               !
               call readi(iref)
               call readi(jref)
               !
               iobject_ = Nobjects-3
               !
               include_state = .false.
               loop_istate_abquad : do i=1,iquad
                   !
                   if (iref==quadrupoletm(i)%iref.and.jref==quadrupoletm(i)%jref) then
                     include_state = .true.
                     !
                     iabi_ = sum(iobject(1:Nobjects-3)) + i
                     !
                     exit loop_istate_abquad
                   endif
                   !
               enddo loop_istate_abquad
               !
             case("HFCC-BF")
                call readi(iref) ; jref = iref
                if (nitems>2) call readi(jref)
                !
                include_state = .false.
                do i = 1, hfcc1(1)%num_field
                    if (iref == hfcc1(1)%field(i)%iref.and.jref == hfcc1(1)%field(i)%jref) then
                      include_state = .true.
                      iabi_ = sum(iobject(1:21-1)) + i
                      exit
                    endif
                enddo
                !
             case("HFCC-A")
                call readi(iref) ; jref = iref
                if (nitems>2) call readi(jref)
                !
                include_state = .false.
                do i = 1, hfcc1(2)%num_field
                    if (iref == hfcc1(2)%field(i)%iref.and.jref == hfcc1(2)%field(i)%jref) then
                      include_state = .true.
                      iabi_ = sum(iobject(1:22-1)) + i
                      exit
                    endif
                enddo
                !
             case("HFCC-C")
                call readi(iref) ; jref = iref
                if (nitems>2) call readi(jref)
                !
                include_state = .false.
                do i = 1, hfcc1(3)%num_field
                    if (iref == hfcc1(3)%field(i)%iref.and.jref == hfcc1(3)%field(i)%jref) then
                      include_state = .true.
                      iabi_ = sum(iobject(1:23-1)) + i
                      exit
                    endif
                enddo
                !             
             case("HFCC-D")
                call readi(iref) ; jref = iref
                if (nitems>2) call readi(jref)
                !
                include_state = .false.
                do i = 1, hfcc1(4)%num_field
                    if (iref == hfcc1(4)%field(i)%iref.and.jref == hfcc1(4)%field(i)%jref) then
                      include_state = .true.
                      iabi_ = sum(iobject(1:24-1)) + i
                      exit
                    endif
                enddo
                !
             case("HFCC-CI")
                call readi(iref) ; jref = iref
                if (nitems>2) call readi(jref)
                !
                include_state = .false.
                do i = 1, hfcc1(5)%num_field
                    if (iref == hfcc1(5)%field(i)%iref.and.jref == hfcc1(5)%field(i)%jref) then
                      include_state = .true.
                      iabi_ = sum(iobject(1:25-1)) + i
                      exit
                    endif
                enddo
                !
             case("HFCC-EQQ0")
                call readi(iref) ; jref = iref
                if (nitems>2) call readi(jref)
                !
                include_state = .false.
                do i = 1, hfcc1(6)%num_field
                    if (iref == hfcc1(6)%field(i)%iref.and.jref == hfcc1(6)%field(i)%jref) then
                      include_state = .true.
                      iabi_ = sum(iobject(1:26-1)) + i
                      exit
                    endif
                enddo
                !
             case("HFCC-EQQ2")
                call readi(iref) ; jref = iref
                if (nitems>2) call readi(jref)
                !
                include_state = .false.
                do i = 1, hfcc1(7)%num_field
                    if (iref == hfcc1(7)%field(i)%iref.and.jref == hfcc1(7)%field(i)%jref) then
                      include_state = .true.
                      iabi_ = sum(iobject(1:27-1)) + i
                      exit
                    endif
                enddo
                !
             end select
             !
             if (.not.include_state) then
                 !write(out,"('The ab initio potential  ',i8,' is skipped')") iref
                 iabi = iabi - 1
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 cycle
             endif
             !
             field => abinitio(iabi_)
             !
             field%iref = iref
             field%jref = jref
             !
             field%refvalue = 0
             !
             !field%istate = istate_
             !field%jstate = jstate_
             !
             if (.not.include_state) then
                 !write(out,"('The ab initio potential  ',i8,' is skipped')") iref
                 iabi = iabi - 1
                 do while (trim(w)/="".and.trim(w)/="END")
                   call read_line(eof,iut) ; if (eof) exit
                   call readu(w)
                 enddo
                 cycle
             endif
             !
             field => abinitio(iabi_)
             !
             field%iref = iref
             field%jref = jref
             field%class = "ABINITIO-"//trim(w)
             !
             field%refvalue = 0
             !
             loop_istate_ai : do istate=1,Nestates
               do jstate=1,Nestates
                 !
                 if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
                   field%istate = istate
                   field%jstate = jstate
                   exit loop_istate_ai
                 endif
                 !
               enddo
             enddo loop_istate_ai
             !
          case default
             call report ("Unrecognized keyword (error 02): "//trim(w),.true.)
          end select
          !
          ! refnumbers of the states to couple
          !
          call read_line(eof,iut) ; if (eof) exit
          call readu(w)
          !
          do while (trim(w)/="".and.trim(w)/="END") !we read until the END the object input block
            !
            select case(w)
            !
            case("INTERPOLATIONTYPE")
              !
              call readu(w)
              field%interpolation_type = trim(w)
              !
            case("TYPE")
              !
              call readu(w)
              !
              field%type = trim(w)
              !
            case("NAME")
              !
              call reada(w)
              !
              field%name = trim(w)
              !
            case("LAMBDA")
              !
              call readi(field%lambda)
              field%lambdaj = field%lambda
              if (nitems>2) call readi(field%lambdaj)
              !
            case("SHIFT","REFVALUE","REF","F1","V0","VE")
              !
              call readf(field%refvalue)
              !
            case("<X|LZ|Y>")
              !
              if (nitems<=2) call report ("Too few entries in "//trim(w),.true.)
              !
              item_ = 1
              do while (trim(w)/="".and.trim(w)/="END".and.item_<nitems)
                !
                call readu(w)
                !
                item_ = item_ + 1
                !
                select case (trim(w))
                  !
                case('0')
                  x_lz_y_ = 0
                case('I')
                  x_lz_y_ = 1
                case('-I')
                  x_lz_y_ = -1
                case('2*I','2I','I*2')
                  x_lz_y_ = 2
                case('-2*I','-2I','-I*2')
                  x_lz_y_ = -2
                case('3*I','3I','I*3')
                  x_lz_y_ = 3
                case('-3*I','-3I','-I*3')
                  x_lz_y_ = -3
                case('4*I','4I','I*4')
                  x_lz_y_ = 4
                case('-4*I','-4I','-I*4')
                  x_lz_y_ = -4
                case('5*I','5I','I*5')
                  x_lz_y_ = 5
                case('-5*I','-5I','-I*5')
                  x_lz_y_ = -5
                case('6*I','6I','I*6')
                  x_lz_y_ = 6
                case('-6*I','-6I','-I*6')
                  x_lz_y_ = -6
                !case('(-I)*N')
                !  call readi(x_lz_y_)
                case default
                  call report ("Illegal input field"//trim(w)//"; SHOULD BE IN THE FORM I, -I, 2*I, -2*I...,(-i)*N",.true.)
                end select
                !
                if (item_==2) then
                  !
                  field%ix_lz_y = x_lz_y_
                  !
                  if (poten(field%istate)%ix_lz_y==1000) poten(field%istate)%ix_lz_y = field%ix_lz_y
                  !
                  if (poten(field%istate)%ix_lz_y/=field%ix_lz_y) then
                    !
                    write(out,"('input: ',a,2i4,' <x|lz|y> disagrees with the poten-value ',2i8)") & 
                            trim(field%class),field%iref,field%jref,poten(field%istate)%ix_lz_y/=field%ix_lz_y
                    call report (" <x|lz|y> disagrees with the poten-value",.true.)
                    !
                  endif
                  !
                elseif(item_==3) then
                  !
                  field%jx_lz_y = x_lz_y_
                  !
                  ! poten can have only one ix_lz_y, i.e. diagonal, its jx_lz_y should not be used
                  !
                  if (poten(field%jstate)%ix_lz_y==1000) poten(field%jstate)%ix_lz_y = field%jx_lz_y
                  !
                  if (poten(field%jstate)%ix_lz_y/=field%jx_lz_y) then
                    !
                    write(out,"('input: ',a,2i4,' <x|lz|y> disagrees with previsouly given <x|Lz|z> value ',2i8)") & 
                            trim(field%class),field%iref,field%jref,poten(field%istate)%ix_lz_y/=field%jx_lz_y
                    call report (" <x|lz|y> disagrees with the previsouly given <x|Lz|z>-value",.true.)
                    !
                  endif
                endif
                !
              enddo                
              !
            case("UNITS")
              !
              item_ = 1
              do while (trim(w)/="".and.trim(w)/="END".and.item_<nitems)
                !
                call readu(w)
                !
                item_ = item_ + 1
                !
                select case(w)
                  !
                case ('BOHR', 'BOHRS')
                  !
                  unit_r = bohr
                  !
                case ('ANG','ANGSTROM','ANGSTROMS')
                  !
                  unit_r = 1.0_rk
                  !
                case ('CM-1')
                  !
                  unit_field = 1.0_rk
                  !
                case ('EV')
                  !
                  unit_field = ev
                  !
                case ('HARTREE','EH','A.U.', 'AU')
                  !
                  unit_field = hartree
                  !
                  if (trim(field%class)=="DIPOLE") unit_field = todebye
                  if (trim(field%class)=="QUADRUPOLE") unit_field = 1.0_rk
                  !
                case ('EA0')
                  !
                  unit_field = todebye
                  !
                case ('DEBYE')
                  !
                  unit_field = 1.0_rk
                  !
                case default
                  !
                  call report ("Illegal input field"//trim(w),.true.)
                  !
                end select
                !
              enddo
              !
            case("SYM","SYMMETRY")
              !
              item_ = 1
              do while (trim(w)/="".and.trim(w)/="END".and.item_<nitems)
                !
                call readu(w)
                !
                item_ = item_ + 1
                !
                select case(w)
                  !
                case ('G')
                  !
                  field%parity%gu = 1
                  !
                case ('U')
                  !
                  field%parity%gu = -1
                  !
                case ('+')
                  !
                  field%parity%pm = 1
                  !
                case ('-')
                  !
                  field%parity%pm = -1
                  !
                case default
                  !
                  call report ("Illegal input field"//trim(w),.true.)
                  !
                end select
                !
              enddo
              !
            case("WEIGHTING")
              !
              call readu(field%weighting%wtype)
              !
              if (trim(field%weighting%wtype)=="PS1997") then 
                !
                call readf(field%weighting%alpha)
                call readf(field%weighting%Vtop)
                !
              endif 
              !
            case("MORPHING","MORPH")
              !
              field%morphing = .true.
              !
            case("ADJUST")
              ! W. Somogyi : add specified constant value to all field values
              !              e.g for shifting PECs when there is a known error
              !
              call readf(field%adjust_val)
              !
              if (nitems>2) then
                call readu(w)
                !
                select case (trim(w))
                  !
                case ('CM-1')
                  !
                  unit_adjust = 1.0_rk
                  !
                case ('EV')
                  !
                  unit_adjust = ev
                  !
                case ('HARTREE','EH','A.U.', 'AU')
                  !
                  unit_adjust = hartree
                  !
                  if (trim(field%class)=="DIPOLE") unit_adjust = todebye
                  if (trim(field%class)=="QUADRUPOLE") unit_adjust = 1.0_rk
                  !
                case ('EA0')
                  !
                  unit_adjust = todebye
                  !
                case ('DEBYE')
                  !
                  unit_adjust = 1.0_rk
                  !
                case default
                  !
                  call report ("Illegal input field"//trim(w),.true.)
                  !
                end select
              else
                unit_adjust = 1.0_rk
              endif
              !
              field%adjust_val = field%adjust_val * unit_adjust
              field%adjust = .true.
              !
            case("MOLPRO")
              !
              field%molpro = .true.
              !
            case("ASYMPTOTE")
              !
              call readf(field%asymptote)
              !
            case("INTEGRATION")
              !
              call readu(field%integration_method)
              !
            case("FACTOR")
              !
              field%complex_f = cmplx(1.0_rk,0.0_rk, kind=rk)
              field%factor = 1.0_rk
              !
              do while (item<min(Nitems,3))
                !
                call readu(w)
                !
                select case (trim(w))
                  !
                case('I')
                  field%complex_f = cmplx(0.0_rk,1.0_rk, kind=rk)
                case('-I')
                  field%complex_f = cmplx(0.0_rk,-1.0_rk, kind=rk)
                case('SQRT(2)')
                  field%factor =  field%factor*sqrt(2.0_rk)
                case('-SQRT(2)')
                  field%factor = -field%factor*sqrt(2.0_rk)
                case default
                  read(w,*,iostat=ierr) x0
                  if(ierr ==0) then
                    field%factor = field%factor*x0
                  else
                    write(out, '(A)') 'Fatal error reading FACTOR. Trying to read as a number: ' // trim(w)
                    write(out, '(A)') 'Stopping now. Please check the input.'
                    stop
                  endif 
                end select
                !
             enddo
             !
             !
            case("FIT_FACTOR")
              !
              call readf(field%fit_factor)
              !
            case("COMPLEX")
              !
              call readu(w)
              !
              select case (trim(w))
                !
              case('0')
                field%complex_f = cmplx(0.0_rk,0.0_rk, kind=rk)
              case('I')
                field%complex_f = cmplx(0.0_rk,1.0_rk, kind=rk)
              case('-I')
                field%complex_f = cmplx(0.0_rk,-1.0_rk, kind=rk)
              case('1')
                field%complex_f = cmplx(1.0_rk,0.0_rk, kind=rk)
              case('SQRT(2)')
                field%complex_f = cmplx(sqrt(2.0_rk),0.0_rk, kind=rk)
              case('-SQRT(2)')
                field%complex_f = cmplx(-sqrt(2.0_rk),0.0_rk, kind=rk)
              case('I*SQRT(2)')
                field%complex_f = cmplx(0.0_rk,sqrt(2.0_rk), kind=rk)
              case('-I*SQRT(2)')
                field%complex_f = cmplx(0.0_rk,-sqrt(2.0_rk), kind=rk)
              case default
                call report ("Illegal input field"//trim(w),.true.)
              end select
              !
            case("<")
              !
              call report ("Braket input structure < | | > is obsolete ",.true.)
              !
              if (Nitems/=9) call report ("Illegal number of characters in the bra-ket field, has to be 9",.true.)
              !
              field%nbrakets = field%nbrakets + 1
              ibraket = ibraket + 1
              !
              if (ibraket>4) then
                 write(out,'("Two many bra-kets = ",i4," (maximum 4) ")') ibraket
                 call report ("Two many bra-kets ",.true.)
              endif
              !
              call readi(field%braket(ibraket)%ilambda)
              call readf(field%braket(ibraket)%sigmai)
              !
              call readu(w)
              if (trim(w)/="|") call report ("Illegal character between bra and ket, has to be a bar |"//trim(w),.true.)
              !
              call readi(field%braket(ibraket)%jlambda)
              call readf(field%braket(ibraket)%sigmaj)
              call readu(w)
              if (trim(w)/=">") call report ("Illegal character after the bra-ket, has to be >"//trim(w),.true.)
              call readu(w)
              if (trim(w)/="=") call report ("Illegal character before the value of the bra-ket, has to be ="//trim(w),.true.)
              !
              call readf(field%braket(ibraket)%value)
              !
              istate = field%istate
              jstate = field%jstate
              !
              if (abs(field%braket(ibraket)%ilambda)/=poten(istate)%lambda.or. &
                  abs(field%braket(ibraket)%jlambda)/=poten(jstate)%lambda) then
                write(out,'("bra-ket lambdas (",2i4,") do not agree with the state lambdas (",2i4,") ")') &
                           field%braket(ibraket)%ilambda, & 
                           field%braket(ibraket)%jlambda,poten(istate)%lambda,poten(jstate)%lambda
                stop "Illegal bra-ket lambdas"
              endif
              !
              if ( nint( abs( 2.0*field%braket(ibraket)%sigmai ) )>nint( 2.0*field%spini ) .or. &
                   nint( abs( 2.0*field%braket(ibraket)%sigmaj ) )>nint( 2.0*field%spinj ) ) then
                write(out,'("bra-ket sigmai (",2f8.1,") greater than the field spini (",2f8.1,") ")') &
                            field%braket(ibraket)%sigmai,field%braket(ibraket)%sigmaj,field%spini,field%spinj
                stop "Illegal bra-ket sigmai"
              endif
              !
            case("SIGMA")
              !
              call readf(field%sigmai)
              field%sigmaj = field%sigmai
              if (nitems>2) call readf(field%sigmaj)
              !
            case("SPIN")
              !
              call readf(field%spini)
              field%spinj = field%spini
              if (nitems>2) call readf(field%spinj)
              !
              if (mod(nint(2.0_rk*field%spini+1.0_rk),2)==0.and.integer_spin) then
                call report("The spin of the field is inconsistent with the multiplicity of J-s in J_list/Jrot (top)",.true.)
              endif
              !
            case("MULT","MULTIPLICITY")
              !
              call readi(field%multi)
              !
              if ( field%multi < 1 ) then ! check that multiplicity >=1
                !
                write(out,'(A,i4,A,i4,A)') "The multiplicity ",field%multi, " of field ", field%iref, " is less than one"
                call report("Error: multiplicity should be an integer greater or equal to one")
              endif
              !
              if (mod(field%multi,2)==0.and.integer_spin) then
                !
                write(out,'(A,i4," of the field ",i4," is inconsistent with multiplicity of jmin/jmax = ",2f8.1)') &
                            "The multiplicity ", field%multi,field%iref,jmin,jmax
                write(out,'("Please check that Jrot at the top of input is integer/half-integer.")') 
                call report("The multiplicity of the field is inconsistent with Jrot/Jlist")
                !
              endif
              !
              !if (mod(field%multi,2)==0) integer_spin = .true.
              !
              field%spini = real(field%multi-1,rk)*0.5_rk
              field%spinj = field%spini ! set the `ket' spin to the same value as the bra by default
              !
              field%jmulti = field%multi
              if (nitems>2) then
                call readi(field%jmulti)
                field%spinj = real(field%jmulti-1,rk)*0.5_rk
              endif
              !
            case("NPARAM","N","NPOINTS")
              !
              ! Obsolete
              !
              call readi(Nparam)
              !
              if (Nparam<0) then
                  call report ("The size of the potential is illegar (<1)",.true.)
              endif
              !
              field%Nterms = Nparam
              !
            case("VALUES")
              !
              ! by Lorenzo Lodi
              ! find the number of points/parameters in input
              !
              Nparam_check = 0
              !write(my_fmt, '(A,I0,A)') '(A',cl,')'
              !
              call input_options(echo_lines=.false.)
              !
              if (field%molpro.and.(field%ix_lz_y==1000.or.field%jx_lz_y==1000)) then
                write(out,"('For MOLPRO-X representaion please define <x|lx|y>',a,2i4)") trim(field%class),field%iref,field%jref
                stop '<x|lz|y> is undefined in dipole-x or spin-orbit-x'
              endif
              !
              do while (trim(w)/="END")
                 !
                 call read_line(eof,iut) ; if (eof) exit
                 !
                 call readu(w)
                 !
                 Nparam_check = Nparam_check+1
                 !
              enddo
              !
              Nparam_check = Nparam_check-1
              !
              if( job%zEchoInput) call input_options(echo_lines=.true.)
              !
              if (trim(w) /= "END") then
                  call report ("ERROR: Cannot find `END' statement)",.true.)
              endif
              !
              ! go back to beginning of VALUES block and reset `w' to original value
              do i=1, Nparam_check+1
                backspace(unit=iut)
              enddo
              !
              w = "VALUES"
              !
              Nparam = Nparam_check
              !
              if (Nparam <= 0) then
                  call report ("ERROR: Number of points or parameters <= 0 )",.true.)
              endif
              !
              field%Nterms = Nparam
              !
              ! Allocation of the pot. parameters
              !
              allocate(field%value(Nparam),field%forcename(Nparam),field%grid(Nparam),field%weight(Nparam),stat=alloc)
              call ArrayStart(trim(field%type),alloc,Nparam,kind(field%value))
              call ArrayStart(trim(field%type),alloc,Nparam,kind(field%grid))
              call ArrayStart(trim(field%type),alloc,Nparam,kind(field%weight))
              !
              allocate(field%link(Nparam),stat=alloc)
              call ArrayStart(trim(field%type),alloc,3*Nparam,ik)
              !
              field%value = 0
              field%forcename = 'dummy'
              field%weight = 0
              !
              iparam = 0
              !
              do while (trim(w)/="".and.iparam<Nparam.and.trim(w)/="END")
                 !
                 call read_line(eof,iut) ; if (eof) exit
                 !
                 iparam = iparam+1
                 !
                 select case(trim(field%type))
                 !
                 case("GRID")
                   !
                   call readu(w)
                   !
                   if (trim(w)=="END") call report("Two many grid-entries in the field "//trim(field%name)// &
                                                                                " (or N is too small)",.true.)
                   !
                   !call readf(f_t)
                   read(w,*) f_t
                   field%grid(iparam) = f_t*unit_r
                   !
                   call readf(f_t)
                   field%value(iparam) = f_t*unit_field
                   !
                   ! these fields are to link to other parameters in the refinement
                   !
                   field%link(iparam)%iobject = 0
                   field%link(iparam)%ifield = 0
                   field%link(iparam)%iparam = 0
                   !
                   write(field%forcename(iparam),"(f18.6)") field%grid(iparam)
                   !
                   field%weight(iparam) = 0
                   !
                   if(nitems>=3) then
                     !
                     call readu(w)
                     !
                     if (trim(w(1:1))=="F") then 
                       !
                       field%weight(iparam) = 1.0_rk
                       !
                       if (nitems>=4) call readu(w)
                       !
                     elseif(trim(field%class(1:4))=="ABIN") then 
                       !
                       read(w,*) field%weight(iparam)
                       !
                     endif
                     !
                     if(trim(w(1:1))=="L".and.nitems>5) then
                       !
                       call readi(field%link(iparam)%iobject)
                       call readi(field%link(iparam)%ifield)
                       call readi(field%link(iparam)%iparam)
                       !
                       ! set the weight of the linked parameter to zero
                       !
                       field%weight(iparam) = 0
                       !
                     endif
                     !
                   endif 
                   !
                   job%total_parameters = job%total_parameters + 1
                   !
                 case ("NONE")
                   !
                   call report ("The field type (e.g. GRID) is undefined for the current field "//trim(w),.true.)
                   !
                 case default
                   !
                   if (nitems<2) then
                      !
                      write(out,"(a,i4)") "wrong number of records for an analytical field-type," // &
                                      "must be two at least (name value)",nitems
                      call report ("illegal number of records (<2) in the current field-line "//trim(w),.true.)
                      !
                   endif
                   !
                   ! these fields are to link to other parameters in the refinement
                   !
                   field%link(iparam)%iobject = 0
                   field%link(iparam)%ifield = 0
                   field%link(iparam)%iparam = 0
                   !
                   call readu(field%forcename(iparam))
                   call readf(f_t)
                   !
                   field%value(iparam) = f_t*unit_field
                   field%weight(iparam) = 0
                   !
                   if(nitems>=3) then
                     !
                     call readu(w)
                     !
                     if (trim(w(1:1))=="F") then 
                       !
                       field%weight(iparam) = 1.0_rk
                       !
                       if (nitems>=4) call readu(w)
                       !
                     elseif(trim(w(1:1))/="L") then
                       !
                       ! old input? 
                       !
                       field%weight(iparam) = f_t
                       read(w,*) field%value(iparam)
                       !
                       if (nitems>=4) call readu(w)
                       !
                     endif
                     !
                     if(trim(w(1:1))=="L".and.nitems>5) then
                       !
                       call readi(field%link(iparam)%iobject)
                       call readi(field%link(iparam)%ifield)
                       call readi(field%link(iparam)%iparam)
                       !
                       ! set the weight of the linked parameter to zero
                       !
                       field%weight(iparam) = 0
                       !
                     endif
                     !
                   endif 
                   !
                   job%total_parameters = job%total_parameters + 1
                   !
                 end select
                 !
              enddo
              !
              call readu(w)
              !
              if (iparam/=Nparam.or.(trim(w)/="".and.trim(w)/="END")) then
                 !
                 print "(2a,2i6)","wrong number of rows in section: ",trim(field%name),iparam,Nparam
                 call report("illegal number of rows in a section: "//trim(w),.true.)
                 !
              endif
                 !
            case default
                 !
                 call report ("Unrecognized keyword (error 03): "//trim(w),.true.)
                 !
            end select
            !
            call read_line(eof,iut) ; if (eof) exit
            call readu(w)
          enddo
         !
       case("SYMGROUP","SYMMETRY","SYMM","SYM","SYM_GROUP") 
         !
         if (symmetry_defined) then 
            !
            write (out,"('input: Symmetry is already defined, SYMGROUP should not apear after INTENSITY')") 
            stop 'input - SYMGROUP should not apear after INTENSITY'
            !
         endif 
         !
         call readu(w)
         !
         job%symmetry = trim(w)
         !
         if (trim(job%symmetry)=="CS") job%symmetry = "CS(M)"
         if (trim(job%symmetry)=="C2V") job%symmetry = "C2V(M)"
         !
         if (trim(job%symmetry)/="CS(M)".and.trim(job%symmetry)/="C2V(M)") then 
           call report ("SYMGROUP: ONLY CS(M) or C2V(M) are allowed, not "//trim(w),.true.)
         endif 
         !
         ! Initialize the group symmetry 
         !
         call SymmetryInitialize(job%symmetry)
         !
         symmetry_defined = .true.
         !
         allocate(job%isym_do(sym%Nrepresen),stat=alloc)
         if (alloc/=0)  stop 'input, isym_do - out of memory'
         !
         job%isym_do = .true.
         !
       case("INTENSITY")
         !
         ! skip if intensity NONE
         !
         if (Nitems>1) then
           call readu(w)
           if (trim(w)=="NONE".or.trim(w)=="OFF") then
             do while (trim(w)/="".and.trim(w)/="END")
               call read_line(eof,iut) ; if (eof) exit
               call readu(w)
             enddo
             cycle
           endif
         endif
         !
         if (Nestates==0) then 
            !
            write (out,"('input: INTENSITY cannot appear before anypoten entries')") 
            stop 'input - INTENSITY defined before all poten entries'
            !
         endif 
         !
         if (.not.symmetry_defined) then 
           !
           ! Initialize symmetry if it has not been done before 
           !
           if ((m1<small_.or.m2<small_).and.(trim(symbol1)=="Undefined".or.trim(symbol2)=="Undefined")) then 
               write(out,"('at least one of MASSES or ATOMS should be defined before INTENSITY')")
               call report("either masses or atom should be defined before INTENSITY",.true.) 
           endif 
           !
           if ( ( m1>0.and.abs(m1-m2)<small_ ).or.( trim(symbol1)/="Undefined".and.trim(symbol1)==trim(symbol2) ) ) then
             !
             job%symmetry = "C2V(M)"
             !
           else
             !
             job%symmetry = "CS(M)"
             !
           endif
           !
           ! Initialize the group symmetry 
           !
           call SymmetryInitialize(job%symmetry)
           !
           symmetry_defined = .true.
           !
           allocate(job%isym_do(sym%Nrepresen),stat=alloc)
           if (alloc/=0)  stop 'input, isym_do - out of memory'
           !
           job%isym_do = .true.
           !
         endif 
         !
         allocate(intensity%gns(sym%Nrepresen),intensity%isym_pairs(sym%Nrepresen),stat=alloc)
         if (alloc/=0) stop 'input, intensity-arrays - out of memory'
         !
         if (sym%Nrepresen==4) then
           !
           intensity%isym_pairs(1:2) = 1
           intensity%isym_pairs(3:4) = 2
           !
         else
           !
           intensity%isym_pairs(:) = 1
           !
         endif
         !
         ! defauls values
         !
         intensity%gns = 1
         forall(i=1:sym%Nrepresen) intensity%isym_pairs(i) = 1
         !
         call read_line(eof,iut) ; if (eof) exit
         call readu(w)
         !
         do while (trim(w)/="".and.trim(w)/="END")
           !
           select case(w)
           !
           case('NONE','ABSORPTION','EMISSION','TM','DIPOLE-TM','PARTFUNC')
             !
             intensity%action = trim(w)
             !
             if (trim(intensity%action)=='DIPOLE-TM') intensity%action = 'TM'
             !
             ! L Lodi: character length must be specified in the array constructor for Fortran 2003 conformance 
             if (any(trim(intensity%action)==(/ character(len=wl) :: 'TM','ABSORPTION','EMISSION','PARTFUNC'/))) then
               !
               action%intensity = .true.
               intensity%do = .true.
               !
             endif
             !
           case ('RWF')
             !
             action%RWF = .true.
             !
             if (nitems>1) call readi(intensity%N_RWF_Order) 
             !
             ! we will need the vibrational basis for RWF
             !
             job%basis_set  = 'KEEP'
             !
           case ('LANDE')
             !
             intensity%lande_calc= .true.
             !
           case('LINELIST')
             !
             call reada (intensity%linelist_file)
             !
           case('MATELEM','RICHMOL')
             !
             intensity%matelem = .true.
             !
           case('RAMAN','POLARIZABILITY')
             !
             action%raman = .true.
             !
           case('QUADRUPOLE')
             !
             action%quadrupole = .true.
             !action%dipole     = .false.
             !
           case('OVERLAP')
             !
             intensity%overlap = .true.
             if (nitems>1) call readu(w) 
             if (trim(w)=="OFF") intensity%overlap = .false.
             !
           case('RENORM','RENORMALIZE')
             !
             intensity%renorm = .true.
             job%basis_set='KEEP'
             !
           case('BOUND')
             !
             intensity%bound = .true.
             job%basis_set='KEEP'
             !
           case('UNBOUND')
             !
             intensity%unbound = .true.
             job%basis_set='KEEP'
             !
           case('VIB-DIPOLE','MU')
             !
             intensity%tdm = .true.
             !
             if (nitems>1) call readu(w) 
             if (trim(w)=="OFF") intensity%tdm = .false.
             !
           case('VIB-QUADRUPOLE')
             !
             intensity%tqm = .true.
             !
             if (nitems>1) call readu(w)
             if (trim(w)=="OFF") intensity%tqm = .false.
             !
           case('THRESH_INTES','THRESH_TM','THRESH-INTES','THRESH_INTENS','THRESH-INTENS')
             !
             call readf(intensity%threshold%intensity)
             !
           case('THRESH_LINE','THRESH_LINESTRENGHT','THRESH_EINSTEIN','THRESH-EINSTEIN')
             !
             call readf(intensity%threshold%linestrength)
             !
           case('THRESH_DIPOLE')
             !
             call readf(intensity%threshold%dipole)
             !
           case('THRESH_COEFF','THRESH_COEFFICIENTS')
             !
             call readf(intensity%threshold%coeff)
             !
           case('TEMPERATURE')
             !
             call readf(intensity%temperature)
             !
           case('QSTAT','PARTITION','PART_FUNC','Q','PART-FUNC')
             !
             call readf(intensity%part_func)
             !
           case ("NSPIN","NUCLEAR-SPIN","NSPINS")
             !
             if ((m1<small_.or.m2<small_).and.(trim(symbol1)=="Undefined".or.trim(symbol2)=="Undefined")) & 
                 call report("masses or atom should be defined before nuclear-spins",.true.)
             !
             call readf(Nspin1)
             call readf(Nspin2)
             !
             if ( ( m1>0.and.abs(m1-m2)<small_ ).or.( trim(symbol1)/="Undefined".and.trim(symbol1)==trim(symbol2) ) ) then
               !
               gns_a = 0.5_rk*((2.0_rk*Nspin1+1.0_rk)**2+(2.0_rk*Nspin1+1.0_rk))
               gns_b = 0.5_rk*((2.0_rk*Nspin1+1.0_rk)**2-(2.0_rk*Nspin1+1.0_rk))
               !
               if (mod(nint(2.0*Nspin1),2)==1) then 
                 !
                 intensity%gns(1:2) = gns_b
                 intensity%gns(3:4) = gns_a
                 !
               else
                 !
                 intensity%gns(1:2) = gns_a
                 intensity%gns(3:4) = gns_b
                 !
               endif
               !
             else
               !
               intensity%gns(1:2) = (2.0_rk*Nspin1+1.0_rk)*(2.0_rk*Nspin2+1.0_rk)
               !
             endif
             !
           case('GNS')
             !
             i = 0
             !
             do while (item<Nitems.and.i<sym%Nrepresen)
               !
               i = i + 1
               call readf(intensity%gns(i))
               !
               ! currently with Cs(M) we have to allow for all symmetries 
               !
               !if (intensity%gns(i)<small_) job%isym_do(i) = .false.
               !
             enddo
             !
             if (i/=sym%Nrepresen.and.sym%Nrepresen==2) then 
               !
               intensity%gns(2) = intensity%gns(1)
               !
             elseif (i/=sym%Nrepresen.and.sym%Nrepresen==2) then
               !
               write (out,"('input: illegal number entries in gns: ',i8,' /= ',i8)") i,sym%Nrepresen
               stop 'input - illegal number entries in gns'
               !
             endif 
             !
           case('SELECTION','SELECTION_RULES','SELECT','PAIRS')
             !
             i = 0
             !
             do while (trim(w)/="".and.i<sym%Nrepresen)
               !
               i = i + 1
               !
               call readi(intensity%isym_pairs(i))
               !
             enddo
             !
             if (i/=sym%Nrepresen) then 
               !
               write (out,"('input: illegal number entries in SELECTION: ',i8,' /= ',i8)") i,sym%Nrepresen
               stop 'input - illegal number entries in SELECTION'
               !
             endif 
             !
           case('ZPE')
             !
             call readf(intensity%zpe)
             job%zpe = intensity%zpe
             job%shift_to_zpe = .false.
             !
           case('J','JLIST','JROT')
             !
             call readf(intensity%j(1))
             !
             if (nitems>3) then 
               call readu(w)
               if (trim(w)/="-") call report ("Unrecognized delimeter, can be comma or dash "//trim(w),.true.)
             endif
             !
             call readf(intensity%j(2))
             !
           case('FREQ-WINDOW','FREQ','NU','FREQUENCY')
             !
             call readf(intensity%freq_window(1))
             call readf(intensity%freq_window(2))
             !
             if (intensity%freq_window(1)<small_) intensity%freq_window(1) = -small_ 
             !
           case('NPOINTS')
             !
             call readi(intensity%npoints)
             if (mod(intensity%npoints,2)==2) intensity%npoints = intensity%npoints + 1
             !
           case('GAMMA','FWHM')
             !
             call readf(intensity%gamma)
             !
           case('LORENTZIAN','GAUSSIAN')
             !
             intensity%RWF_type  = trim(w)
             !
             call readf(intensity%gamma)
             !
           case('ENERGY')
             !
             call readu(w)
             !
             do while (trim(w)/="")
                !
                select case(w)
                !
                case("LOWER","LOW","L")
                  !
                  call readf(intensity%erange_low(1))
                  call readf(intensity%erange_low(2))
                  !
                  if (intensity%erange_low(1)<small_) intensity%erange_low(1)= -small_ 
                  !
                case("UPPER","UPP","UP","U")
                  !
                  call readf(intensity%erange_upp(1))
                  call readf(intensity%erange_upp(2))
                  !
                  if (intensity%erange_upp(1)<small_) intensity%erange_upp(1)=-small_ 
                  !
                end select 
                !
                call readu(w)
                !
             enddo 
             !
           !case('LOWER','LOW','L')
           !  !
           !  call readu(w)
           !  !
           !  call input_quanta(w,intensity%lower)
           !  !
           !case('UPPER','UP','U')
           !  !
           !  call readu(w)
           !  !
           !  call input_quanta(w,intensity%upper)
           !  !
           case default
             !
             call report ("Unrecognized keyword (error 04): "//trim(w),.true.)
             !
           end select 
           !
           call read_line(eof,iut) ; if (eof) exit
           call readu(w)
           !
         enddo
         !
         if (trim(intensity%action) == 'ABSORPTION'.or.trim(intensity%action) == 'EMISSION') then 
           !
           ! define the selection pairs by gns if not yet defined
           !
           i_t = 0
           !
           do i = 1,sym%Nrepresen
             !
             if (intensity%isym_pairs(i)/=0) cycle
             !
             do j = 1,sym%Nrepresen
               !
               if (i/=j.and.intensity%isym_pairs(j)==0.and.intensity%gns(i)==intensity%gns(j)) then 
                 !
                 i_t = i_t + 1
                 !
                 intensity%isym_pairs(i) = i_t
                 intensity%isym_pairs(j) = i_t
                 !
               endif 
               !
             enddo
             !
           enddo
           !
         endif 
         !
       case default
         call report ("Principal keyword "//trim(w)//" not recognized",.true.)
       end select
       !
    end do
    !
    Nestates = iobject(1)
    !
    ! make sure job%vibmax is not larger than npoints 
    do i = 1,Nestates
      if ( job%vibmax(i) == 1e8 ) job%vibmax(i) = grid%npoints-1
    enddo
    !
    if (Nestates<1) call report ("At least one POTEN object must be present (abinitio poten does not count)",.true.)
    !
    Nspinorbits = iso
    Ndipoles = idip
    Nlxly = ilxly
    Nl2   = il2
    Nabi  = 0
    Nss   = iss
    Nsso  = isso
    Nbobrot  = ibobrot
    Nsr = isr
    Ndiabatic = idiab
    Nlambdaopq = iobject(10)
    Nlambdap2q = iobject(11)
    Nlambdaq = iobject(12)
    Nnac = iobject(13)    
    nQuadrupoles = iquad
    !
    ! create a map with field distribution
    !
    do i = 1,Nobjects-4
      fieldmap(i)%Nfields = iobject(i)
    enddo
    !
    !fieldmap(1)%Nfields = Nestates
    !fieldmap(2)%Nfields = Nspinorbits
    !fieldmap(3)%Nfields = Nl2
    !fieldmap(4)%Nfields = Nlxly
    !fieldmap(5)%Nfields = Nss
    !fieldmap(6)%Nfields = Nsso
    !fieldmap(7)%Nfields = Nbobrot
    !fieldmap(8)%Nfields = Nsr
    !fieldmap(9)%Nfields = Ndiabatic
    !fieldmap(10)%Nfields = iobject(10)
    !
    fieldmap(Nobjects-3)%Nfields = nQuadrupoles
    fieldmap(Nobjects-2)%Nfields = Nabi
    fieldmap(Nobjects-1)%Nfields = 1  ! Brot
    fieldmap(Nobjects)%Nfields = Ndipoles
    !
    !Ntotalfields = Nestates+Nspinorbits+NL2+NLxLy+Nss+Nsso+Nbobrot+Nsr+Ndiabatic+iobject(10)
    !
    Ntotalfields = sum(iobject(1:Nobjects-4))
    !
    ! check if all abinitio fields are initialized. If not we need to make dummy abinitio fields;
    ! we also check whether not all fields are given on a grid and thus can be varied.
    !
    !if (action%fitting .eqv. .true.) then
      !
      Nabi = Ntotalfields
      fieldmap(Nobjects-2)%Nfields = Nabi
      !
      ! we also check whether not all fields are given on a grid and thus can be varied.
      !
      allgrids = .true.
      !
      iabi = 0
      !
      do iobj = 1,Nobjects-3
        !
        do iterm = 1,fieldmap(iobj)%Nfields
          !
          select case (iobj)
          case (1)
            field => poten(iterm)
          case (2)
            field => spinorbit(iterm)
          case (3)
            field => l2(iterm)
          case (4)
            field => lxly(iterm)
          case (5)
            field => spinspin(iterm)
          case (6)
            field => spinspino(iterm)
          case (7)
            field => bobrot(iterm)
          case (8)
            field => spinrot(iterm)
          case (9)
            field => diabatic(iterm)
          case (10)
            field => lambdaopq(iterm)
          case (11)
            field => lambdap2q(iterm)
          case (12)
            field => lambdaq(iterm)
          case (13)
            field => nac(iterm)
          case (21, 22, 23, 24, 25, 26, 27)
            field => hfcc1(iobj - 20)%field(iterm)
          case (Nobjects-3)
            field => quadrupoletm(iterm)
          case (Nobjects-2)
            field => abinitio(iterm)
          case default
             print "(a,i0)", "iobject = ",iobj
             stop "illegal iobject  "
          end select
          !
          iabi = iabi + 1
          !
          if (trim(field%type)/="GRID") allgrids = .false.
          !
          !field => abinitio(iabi)
          !
          if (.not.associated(abinitio(iabi)%value)) then
            !
            Nparam = 1 ; abinitio(iabi)%Nterms = 0
            !
            allocate(abinitio(iabi)%value(Nparam),abinitio(iabi)%forcename(Nparam),abinitio(iabi)%grid(Nparam), & 
                     abinitio(iabi)%weight(Nparam),stat=alloc)
            call ArrayStart(trim(abinitio(iabi)%type),alloc,Nparam,kind(abinitio(iabi)%value))
            call ArrayStart(trim(abinitio(iabi)%type),alloc,Nparam,kind(abinitio(iabi)%grid))
            call ArrayStart(trim(abinitio(iabi)%type),alloc,Nparam,kind(abinitio(iabi)%weight))
            !
            abinitio(iabi)%value = 0
            abinitio(iabi)%grid = 1.0_rk
            abinitio(iabi)%weight = 0
            abinitio(iabi)%type = 'DUMMY'  ! dummy field
            abinitio(iabi)%name    = field%name
            abinitio(iabi)%spini   = field%spini
            abinitio(iabi)%spinj   = field%spinj
            abinitio(iabi)%sigmai  = field%sigmai
            abinitio(iabi)%sigmaj  = field%sigmaj
            abinitio(iabi)%multi   = field%multi
            abinitio(iabi)%lambda  = field%lambda
            abinitio(iabi)%lambdaj = field%lambdaj
            !
          endif
          !
        enddo
      enddo
      !
      !if (allgrids.and.action%fitting) then
      !  call report ("Fitting is not possible: No field of not the GRID-type!",.true.)
      !endif 
      !
    !endif
    !
    !if (Nabi>Ntotalfields) then
    !    print "(2a,i4,a,i6)",trim(w),": Number of ab initio fields ",iabi," exceeds the total number of fields ",Ntotalfields
    !    call report ("Too many ab initio fields given in the input for"//trim(w),.true.)
    !endif
    !
    if (.not.symmetry_defined) then
         !
         ! Initialize the group symmetry  (NB: use capital letters)
         !
         job%symmetry = "CS(M)"
         !
         call SymmetryInitialize(job%symmetry)
         !
         symmetry_defined = .true.
         !
         allocate(job%isym_do(sym%Nrepresen),stat=alloc)
         if (alloc/=0)  stop 'input, isym_do - out of memory'
         !
         job%isym_do = .true.
    !
    write(out,"('Symmetry was not specified. ',a,' is assumed based on the masses/atoms', /)") trim(job%symmetry)
         !
    endif 
    ! !I think the following message should be outputed only if line intensity
    ! are computed, not all the time.
    !write(out,"('A dipole threshold of ',e18.8,' will be used')") intensity%threshold%dipole
    !
    if (iobject(1)/=nestates) then
      write(out,'("The number of states required ",i8," is inconcistent (smaller) with the number of PECs ",i8," included")') & 
                 nestates,iobject(1)
      stop "Illegal number of states: ipo/=nestates"
    endif
    !
    if ((trim(solution_method)=="LOBATTO".and.grid%nsub /= 6).or.(trim(solution_method)/="LOBATTO".and.grid%nsub == 6)) then 
        write(out,"('Error: The grid 6 can be used without LOBATTO  (key word SOLUTIONMETHOD)')")
        write(out,"('solution_method, grid = ',a,i8)") trim(solution_method),grid%nsub 
    endif
    !
    ! find lowest Jmin allowed by the spin and lambda
    ! this is done by computing all possible |Omegas| = | Lambda + Sigma |
    ! we assume Jmin=min( Omegas )
    !
    !lambda_ = 100
    !spin_ = lambda_

    ! for each term, find the minimum possible Omega
    do istate=1,Nestates
      omega_ = safe_max
      spin_   = poten(istate)%spini
      lambda_ = poten(istate)%lambda
      do i=0, nint( 2._rk * spin_ )
       sigma_=-spin_ + real(i, rk)
       omega_=min(omega_,abs(real(lambda_, rk) +sigma_) )
      end do
      poten(istate)%Omega_min = omega_
    enddo
    !
    ! find globally minimum Omega
    omega_= safe_max
    do istate=1,Nestates
      if( poten(istate)%Omega_min < omega_) omega_ = poten(istate)%Omega_min
    enddo
    !
    !jmin_ = abs( real(lambda_) ) ; 
    !if (.not.integer_spin) jmin_ = abs( real(lambda_)-0.5_rk )
    !
    !jmin_ = abs( real(lambda_)-abs(spin_) ) ! ; if (.not.integer_spin) jmin_ = jmin_-0.5_rk
    !
    !if (.not.integer_spin) jmin_ = abs(jmin_- 0.5_rk)
    !
    jmin = omega_
    if (jmax<jmin) jmax = jmin
    !
    jmin_global = Jmin
    !
    ! check the L2 terms:
    !
    !if (Nl2>nestates) then
    !  write(out,'("The number of L2 components ",i8," is inconsistent with the number of electronic states ",i8," included")') & 
    !              Nl2,nestates
    !  stop "Illegal number of L2 components: Nl2>nestates"
    !endif
    !
    !do istate = 1,Nl2
    !  if (L2(il2)%iref/=poten(L2(il2)%istate)%iref) then
    !    write(out,'("For the state ",i4," the reference number of L2  = ",i8," is inconsistent with the ref-number of the" &
    !                 // " electronic states = ",i8)') istate,L2(istate)%iref,poten(istate)%iref
    !    stop "Illegal ref-number of L2 components: il2-ref/=ipot-ref"
    !  endif
    !enddo
    !

    if( job%zEchoInput) write(out,"('(<--- End of the input.)'/)")

  contains
    !
    subroutine input_quanta(w,field)
      !
      character(len=cl),intent(in) :: w
      type(quantaT),intent(inout)  :: field
      
      do while (trim(w)/="")
        !      
        select case (w)
          !
        case("J")
          !
          call readf(field%jrot)
          !
        case("OMEGA")
          !
          call readf(field%omega)
          !
        case("V")
          !
          call readi(field%v)
          !
        case("LAMBDA")
          !
          call readi(field%ilambda)
          !if (nitems>2) call readi(field%ilambdaj)
          !
        case("STATE")
          !
          call readi(field%istate)
          !
        case("SIGMA")
          !
          call readf(field%sigma)
          !if (nitems>2) call readf(field%sigmaj)
          !
        case("SPIN")
          !
          call readf(field%spin)
          !
          !if (nitems>2) call readf(field%spinj)
          !
          !if (mod(nint(2.0_rk*field%spini+1.0_rk),2)==0.and.integer_spin) then
          !  call report("The multiplicity of j-s in J_lits are inconcistent with the SPIN of the field defined at:",.true.)
          !endif
          !
        end select
        !
        call readu(w)
        !
     enddo 
     !
    end subroutine input_quanta

    !
    ! set some default values of diffferent fields based on the descritpion of the corresponding poten-s
    !
    subroutine set_field_refs(field,iref,jref,istate,jstate)
      !
      integer(ik),intent(in) :: iref,jref,istate,jstate
      type(fieldT),pointer,intent(in)  :: field
         !
         field%iref = iref
         field%jref = jref
         field%istate = istate
         field%jstate = jstate
         if ( field%sigmai<bad_value+1 ) field%sigmai  = field%spini ! `sigma' should be irrelevant; but we define it to `spin'
         if ( field%sigmaj<bad_value+1 ) field%sigmaj  = field%spinj
         !
    end subroutine set_field_refs
    !
    !
    subroutine input_non_diagonal_field(Nobjects,iType,iobject,fields,ierr)
        !
        integer(ik),intent(in) :: iType,Nobjects
        integer(ik),intent(inout)  :: iobject
        integer(ik),intent(out) :: ierr
        type(FieldT),pointer :: fields(:)
        type(FieldT),pointer :: field
        !
        integer(ik) :: iref,jref,istare,jstate,istate_,jstate_
        !
        ierr = 0 
        !
        iobject = iobject + 1
        !
        field => fields(iobject)
        !
        call readi(iref) ; jref = iref
        !
        ! for nondiagonal terms
        if (nitems>2) call readi(jref)
        !
        ! find the corresponding potential
        !
        include_state = .false.
        loop_istate_nd : do istate=1,Nestates
          do jstate=1,Nestates
            if (iref==poten(istate)%iref.and.jref==poten(jstate)%iref) then
              include_state = .true.
              istate_ = istate
              jstate_ = jstate
              exit loop_istate_nd
            endif
          enddo
        enddo loop_istate_nd
        !
        ! Check if it was defined before 
        do istate=1,iobject-1
           if ( iref==fields(istate)%iref.and.jref==fields(istate)%jref ) then
             call report (trim(CLASSNAMES(iType))//" is repeated",.true.)
           endif
        enddo
        !
        if (.not.include_state) then
            !write(out,"('The LAMBDA-Q term ',2i8,' is skipped')") iref,jref
            iobject = iobject - 1
            do while (trim(w)/="".and.trim(w)/="END")
              call read_line(eof,iut) ; if (eof) exit
              call readu(w)
            enddo
            ierr = 1
            return
        endif
        !
        call set_field_refs(field,iref,jref,istate_,jstate_)
        !
        field%class = trim(CLASSNAMES(iType))
        !
        if (action%fitting) call report (trim(field%class)//" cannot appear after FITTING",.true.)
        !
    end subroutine input_non_diagonal_field    
    !
  end subroutine ReadInput




  !
  ! This basis set composition is for contr=VIB, State-Lambda-Sigma
  !
  subroutine define_quanta_bookkeeping(iverbose,jval,Nestates,Nlambdasigmas)
    !
    integer(ik),intent(in) :: iverbose
    real(rk),intent(in)    :: jval
    integer(ik),intent(in) :: nestates
    integer(ik),intent(out) :: Nlambdasigmas ! to count states with different lambda/sigma
    integer(ik) :: ilevel,itau,ilambda,nlevels,multi_max,imulti,istate,multi,alloc,taumax
    real(rk)    :: sigma,omega
    !
    if (iverbose>=4) call TimerStart('Define quanta')
    !
    ilevel = 0
    multi_max = 1
    !
    ! count states
    !
    do istate = 1,nestates
      !
      multi_max = max(poten(istate)%multi,multi_max)
      multi = poten(istate)%multi
      !
      taumax = 1
      if (poten(istate)%lambda==0) taumax = 0
      !
      do itau = 0,taumax
        ilambda = (-1)**itau*poten(istate)%lambda
        !
        sigma = -poten(istate)%spini
        do imulti = 1,multi
          !
          omega = real(ilambda,rk)+sigma
          !
          if (nint(2.0_rk*abs(omega))<=nint(2.0_rk*jval)) then
            ilevel = ilevel + 1
          endif
          sigma = sigma + 1.0_rk
          !
        enddo
        !
      enddo
      !
    enddo
    !
    nlevels = ilevel
    !
    if (allocated(iquanta2ilevel)) then
       deallocate(quanta)
       deallocate(iquanta2ilevel)
       call ArrayStop('quanta')
    endif
    !
    allocate(quanta(nlevels),stat=alloc)
    allocate(iquanta2ilevel(Nestates,0:1,multi_max),stat=alloc)
    call ArrayStart('quanta',alloc,size(iquanta2ilevel),kind(iquanta2ilevel))
    iquanta2ilevel = 1e4
    !
    ! the total number of the spin-lambda-electronic states
    !
    Nlambdasigmas = nlevels
    !
    if (iverbose>=4) write(out,'("The total number sigma/lambda states (size of the sigma-lambda submatrix) = ",i0)') Nlambdasigmas
    !
    ! assign quanta (rotation-spin-sigma)
    !
    ilevel = 0
    !
    if (iverbose>=4) write(out,'(/"Sigma-Lambda basis set:")')
    if (iverbose>=4) write(out,'("     i     jrot  state   spin    sigma lambda   omega")')
    !
    do istate = 1,nestates
      multi = poten(istate)%multi
      !
      taumax = 1
      if (poten(istate)%lambda==0) taumax = 0
      do itau = 0,taumax
        ilambda = (-1)**itau*poten(istate)%lambda
        !
        sigma = -poten(istate)%spini          !set sigma to -S (its most negative possible value
        if( sigma == -0.0_rk) sigma = +0.0_rk  !if sigma=0 use `positive signed' zero (for consistency across compilers)
        !
        do imulti = 1,multi
          !
          omega = real(ilambda,rk)+sigma
          !
          if (nint(2.0_rk*abs(omega))<=nint(2.0_rk*jval)) then
            !
            ilevel = ilevel + 1
            !
            quanta(ilevel)%spin = poten(istate)%spini
            quanta(ilevel)%istate = istate
            quanta(ilevel)%sigma = sigma
            quanta(ilevel)%imulti = multi
            quanta(ilevel)%ilambda = ilambda
            quanta(ilevel)%omega = real(ilambda,rk)+sigma
            iquanta2ilevel(istate,itau,imulti) = ilevel
            !
            ! print out quanta
            !
            if (iverbose>=4) write(out,'(i6,1x,f8.1,1x,i4,1x,f8.1,1x,f8.1,1x,i4,1x,f8.1,3x,a)') & 
                             ilevel,jval,istate,quanta(ilevel)%spin,sigma,ilambda,omega,trim(poten(istate)%name)
            !
          endif
          !
          sigma = sigma + 1.0_rk
          !
        enddo
        !
      enddo
    enddo
    !
    if (iverbose>=4) call TimerStop('Define quanta')
    !
  end subroutine define_quanta_bookkeeping
  !



  !
  !
  ! Book keeping for a reduced basis |Omega>=|Lambda>|Sigma> with Omega = Lambda + Sigma as main QN
  subroutine define_quanta_bookkeeping_Omega(iverbose,omega,Nestates,Nlambdasigmas)
    !
    integer(ik),intent(in) :: iverbose
    real(rk),intent(in)    :: omega
    integer(ik),intent(in) :: nestates
    integer(ik),intent(out) :: Nlambdasigmas ! to count states with different lambda/sigma
    integer(ik) :: ilevel,itau,ilambda,nlevels,multi_max,imulti,istate,multi,alloc,taumax,lambda_max,lambda_min
    real(rk)    :: sigma,omega_
    !
    if (iverbose>=4) call TimerStart('Define quanta')
    !
    ilevel = 0
    multi_max = 1
    lambda_max = 0
    lambda_min = 10000
    !
    do istate = 1,nestates
      !
      multi_max = max(poten(istate)%multi,multi_max)
      lambda_max = max(poten(istate)%lambda,lambda_max)
      lambda_min = min(poten(istate)%lambda,lambda_min)
      !
    enddo
    !
    ! count states
    !
    do istate = 1,nestates
      !
      multi_max = max(poten(istate)%multi,multi_max)
      multi = poten(istate)%multi
      !
      taumax = 1
      if (poten(istate)%lambda==0) taumax = 0
      !
      do itau = 0,taumax
        ilambda = (-1)**itau*poten(istate)%lambda
        !
        sigma = -poten(istate)%spini-1.0_rk
        do imulti = 1,multi
          !
          sigma = sigma + 1.0_rk
          omega_ = real(ilambda,rk)+sigma
          !
          if (nint(omega_-omega)/=0) cycle
          !
          ilevel = ilevel + 1
          !
        enddo
        !
      enddo
      !
    enddo
    !
    nlevels = ilevel
    !
    if (allocated(iquanta2ilevel)) then
       deallocate(quanta)
       deallocate(iquanta2ilevel)
       call ArrayStop('quanta')
    endif
    !
    allocate(quanta(nlevels),stat=alloc)
    allocate(iquanta2ilevel(Nestates,0:1,multi_max),stat=alloc)
    call ArrayStart('quanta',alloc,size(iquanta2ilevel),kind(iquanta2ilevel))
    iquanta2ilevel = 1e4
    !
    ! the total number of the spin-lambda-electronic states
    !
    Nlambdasigmas = nlevels
    !
    if (iverbose>=6) write(out,'("The total number sigma/lambda states (size of the sigma-lambda submatrix) = ",i0)') Nlambdasigmas
    !
    ! assign quanta (rotation-spin-sigma)
    !
    ilevel = 0
    !
    if (iverbose>=6) write(out,'(/"Sigma-Omega basis set:")')
    !
    do istate = 1,nestates
      multi = poten(istate)%multi
      !
      taumax = 1
      if (poten(istate)%lambda==0) taumax = 0
      do itau = 0,taumax
        ilambda = (-1)**itau*poten(istate)%lambda
        !
        sigma = -poten(istate)%spini          !set sigma to -S (its most negative possible value
        if( sigma == -0.0_rk) sigma = +0.0_rk  !if sigma=0 use `positive signed' zero (for consistency across compilers)
        !
        sigma = sigma-1.0_rk
        !
        do imulti = 1,multi
          !
          sigma = sigma + 1.0_rk
          !
          omega_ = real(ilambda,rk)+sigma
          !
          !if (taumax==0.and.poten(istate)%parity%pm==-1) omega_ = - omega_
          !
          if (nint(omega_-omega)/=0) cycle
          !
          ilevel = ilevel + 1
          !
          quanta(ilevel)%spin = poten(istate)%spini
          quanta(ilevel)%istate = istate
          quanta(ilevel)%sigma = sigma
          quanta(ilevel)%imulti = multi
          quanta(ilevel)%ilambda = ilambda
          quanta(ilevel)%omega = real(ilambda,rk)+sigma
          quanta(ilevel)%name = trim(poten(istate)%name)
          iquanta2ilevel(istate,itau,imulti) = ilevel
          !
          ! print out quanta
          !
          if (iverbose>=5) write(out,'(i6,1x,i4,1x,f8.1,1x,f8.1,1x,i4,1x,f8.1,3x,a)') & 
                           ilevel,istate,quanta(ilevel)%spin,sigma,ilambda,omega,trim(poten(istate)%name)
          !
        enddo
        !
      enddo
    enddo
    !
    if (iverbose>=4) call TimerStop('Define quanta')
    !
  end subroutine define_quanta_bookkeeping_Omega


! if chemical symbols of the atoms are given in ATOMS gets from those masses and spins.
! if masses or spins are given explicitely, uses those values
subroutine check_and_set_atomic_data(iverbose)
   use atomic_and_nuclear_data, only: print_atomic_and_nuclear_info, get_z, get_a, get_m, atomic_mass, &
                                      get_name_from_mass !, nuclear_spin
   integer,intent(in) :: iverbose
   integer :: zi, ai, m
   real(kind=rk) m1_ref, m2_ref
   integer :: iselect ! flag use to identify four cases:
                      ! iselect =1 => mass is given and atom symbol is given. The given mass is used, info on the atom printed.
                      ! iselect =2 => mass is given and atom symbol NOT given. The given mass is used, info on the atom printed.
                      ! iselect =3 => mass is NOT given and atom symbol given. Mass taken from database, info on the atom printed.
                      ! iselect =4 => mass is NOT given and atom symbol NOT given. Error message printed.

    iselect = 0
    if( m1 >  0._rk .and. symbol1(1:9) /="Undefined") iselect = 1
    if( m1 >  0._rk .and. symbol1(1:9) =="Undefined") iselect = 2
    if( m1 <= 0._rk .and. symbol1(1:9) /="Undefined") iselect = 3
    if( m1 <= 0._rk .and. symbol1(1:9) =="Undefined") iselect = 4

    select case(iselect)
      case default ! it should never happen...
        if (iverbose>=4) write(out,'(a)') ' Error in subroutine check_and_set_atomic_data (m1): iselect = 0'

      case(1) !mass is given and atom symbol is given
        if (iverbose>=4) write(out,'(A)') 'Atom #1 specified as: '  // trim(symbol1)
        if (iverbose>=4) write(out,'(A)') 'Information from the internal database:'
        call print_atomic_and_nuclear_info(symbol1,iverbose,out)
    !    write(out,'(A)')
        if (iverbose>=4) write(out,'(A, F25.12,/)') 'Using the mass explicitely specified in the input, m1=', m1

      case(2) !mass is given and atom symbol NOT given
        call get_name_from_mass(m1, symbol1)
        if( symbol1(1:7) =='Unknown') then
          write(out,'(a)') 'The given mass could not be matched to any atom in the internal database.' &
                        // ' Are you sure it is right?'
         else
            if (iverbose>=4) write(out,'(a)') 'The given mass is best-matched in the internal database by:'
            call print_atomic_and_nuclear_info(symbol1,iverbose)
        !    write(out,'(a)')
        endif
        if (iverbose>=4) write(out,'(A, F25.12,/)') 'Using the mass explicitely specified in the input, m1=', m1


      case(3) !mass is NOT given and atom symbol given
       if (iverbose>=4) write(out,'(A)') 'Atom #1 specified as: '  // trim(symbol1)
       if (iverbose>=4) write(out,'(A)') 'Information from the internal database:'
       call print_atomic_and_nuclear_info(symbol1,iverbose,out)
       if (iverbose>=4) write(out,'(A)')
       zi = get_z(symbol1) ; ai = get_a(symbol1) ; m  = get_m(symbol1)
       if( zi <= 0) then
         write(out,'(A)') 'Error: cannot use mass from internal database because: element unrecognized.'
         return
       endif

       if( ai >0   .and. m >= 0 ) m1_ref = atomic_mass(zi,ai,m)
       if( ai >0   .and. m <  0 ) m1_ref = atomic_mass(zi,ai)
       if( ai <= 0 .and. m >= 0 ) m1_ref = atomic_mass(zi,m=m)
       if( ai <= 0 .and. m  < 0 ) m1_ref = atomic_mass(zi)

       if( m1_ref <= 0._rk) then
           write(out,'(A)') 'Error: cannot use mass from internal database because: element not found.'
!          write(out,'(A)') 'Possibly radioactive with a half-life < 1 hour.'
           return
       endif
       m1 = m1_ref
       if (iverbose>=4) write(out,'(A,F25.12,/)') 'Using the following mass m1 = ', m1

      case(4) ! mass is NOT given and atom symbol NOT given
        return ! do nothing, error message triggered later
    end select

    ! same thing for atom2 (duplicated code, for now).
    iselect = 0
    if( m2 >  0._rk .and. symbol2(1:9) /="Undefined") iselect = 1
    if( m2 >  0._rk .and. symbol2(1:9) =="Undefined") iselect = 2
    if( m2 <= 0._rk .and. symbol2(1:9) /="Undefined") iselect = 3
    if( m2 <= 0._rk .and. symbol2(1:9) =="Undefined") iselect = 4

    select case(iselect)
      case default ! it should never happen...
        if (iverbose>=4) write(out,'(a)') ' Error in subroutine check_and_set_atomic_data (m2): iselect = 0'

      case(1) !mass is given and atom symbol is given
        if (iverbose>=4) write(out,'(A)') 'Atom #2 specified as: '  // trim(symbol2)
        if (iverbose>=4) write(out,'(A)') 'Information from the internal database:'
        call print_atomic_and_nuclear_info(symbol2,iverbose,out)
      !  write(out,'(A)')
        if (iverbose>=4) write(out,'(A, F25.12,/)') 'Using the mass explicitely specified in the input, m2=', m2

      case(2) !mass is given and atom symbol NOT given
        call get_name_from_mass(m2, symbol2)
        if( symbol2(1:7) =='Unknown') then
          if (iverbose>=4) write(out,'(a)') 'The given mass could not be matched to any atom in the internal database.' &
                        // ' Are you sure it is right?'
         else
            if (iverbose>=4) write(out,'(a)') 'The given mass is best-matched in the internal database by:'
            call print_atomic_and_nuclear_info(symbol2,iverbose)
 !           write(out,'(a)')
        endif
        if (iverbose>=4) write(out,'(A, F25.12,/)') 'Using the mass explicitely specified in the input, m2=', m2


      case(3) !mass is NOT given and atom symbol given
       if (iverbose>=4) write(out,'(A)') 'Atom #2 specified as: '  // trim(symbol2)
       if (iverbose>=4) write(out,'(A)') 'Information from the internal database:'
       call print_atomic_and_nuclear_info(symbol2,iverbose,out)
       if (iverbose>=4) write(out,'(A)')
       zi = get_z(symbol2) ; ai = get_a(symbol2) ; m  = get_m(symbol2)
       if( zi <= 0) then
         write(out,'(A)') 'Error: cannot use mass from internal database because: element unrecognized.'
         return
       endif

       if( ai >0   .and. m >= 0 ) m2_ref = atomic_mass(zi,ai,m)
       if( ai >0   .and. m <  0 ) m2_ref = atomic_mass(zi,ai)
       if( ai <= 0 .and. m >= 0 ) m2_ref = atomic_mass(zi,m=m)
       if( ai <= 0 .and. m  < 0 ) m2_ref = atomic_mass(zi)

       if( m2_ref <= 0._rk) then
           if (iverbose>=4) write(out,'(A)') 'Error: cannot use mass from internal database because: element not found.'
!          write(out,'(A)') 'Possibly radioactive with a half-life < 1 hour.'
           return
       endif
       m2 = m2_ref
       if (iverbose>=4) write(out,'(A,F25.12,/)') 'Using the following mass m2 = ', m2

      case(4) ! mass is NOT given and atom symbol NOT given
        return ! do nothing, error message triggered later
    end select
    !
end subroutine  check_and_set_atomic_data

subroutine map_fields_onto_grid(iverbose)
     !
     use functions,only : define_analytical_field
     !
     character(len=130)     :: my_fmt  !text variable containing formats for reads/writes
     ! 
     integer(ik),intent(in) :: iverbose
     !
     integer(ik)             :: ngrid,alloc,j,nsub,Nmax,iterm,nterms,i,ipotmin,istate,jstate,itotal
     integer(ik)             :: ifterm,iobject,ifield
     real(rk)                :: rmin, rmax, re, alpha, h,sc,h12,scale,check_ai
     real(rk),allocatable    :: f(:)
     !
     integer             :: np   ! tmp variable for extrapolation
     real(kind=rk)       :: x1, x2, y1, y2, aa, bb ! tmp variables used for extrapolation
     !
     real(rk),allocatable    :: spline_wk_vec(:) ! working vector needed for spline interpolation
     real(rk),allocatable    :: spline_wk_vec_B(:) ! working vector needed for spline interpolation
     real(rk),allocatable    :: spline_wk_vec_C(:) ! working vector needed for spline interpolation
     real(rk),allocatable    :: spline_wk_vec_D(:) ! working vector needed for spline interpolation
     real(rk),allocatable    :: spline_wk_vec_E(:) ! working vector needed for spline interpolation
     real(rk),allocatable    :: spline_wk_vec_F(:) ! working vector needed for spline interpolation
     real(rk),allocatable    :: xx(:), yy(:), ww(:)! tmp vectors for extrapolation
     type(linkT),allocatable :: link_(:) ! tmp link values for extrapolation
     real(rk) :: yp1, ypn, Vtop, beta, vmin, DeltaR
     integer(ik) :: imin
     !
     type(fieldT),pointer      :: field
     !
     ngrid = grid%npoints
     !
     if (allocated(r)) then
       deallocate(r,z,d2dr,r2sc)
       call ArrayStop('r-field')
     endif
     !
     !
     if (associated(grid%r)) then
       deallocate(grid%r)
       call ArrayStop('grid%r')
     endif
     !
     allocate(r(ngrid),z(ngrid),f(ngrid),d2dr(ngrid),r2sc(ngrid),grid%r(ngrid),stat=alloc)
     call ArrayStart('grid%r',alloc,ngrid,kind(grid%r))
     call ArrayStart('r-field',alloc,ngrid,kind(r))
     call ArrayStart('f-field',alloc,ngrid,kind(f))
     call ArrayStart('r-field',alloc,ngrid,kind(z))
     call ArrayStart('r-field',alloc,ngrid,kind(d2dr))
     call ArrayStart('r-field',alloc,ngrid,kind(r2sc))
     !
     rmin = grid%rmin
     rmax = grid%rmax
     !
     nsub = grid%nsub
     alpha = grid%alpha
     re = grid%re
     !
     ! mapping grid
     !
     call check_and_set_atomic_data(iverbose)
     ! check masses make sense
     if( m1 <= 0._rk .or. m2 <= 0._rk) then
       write(out, '(A,2F20.8)') 'Illegal value for masses! m1, m2= ', m1, m2
       write(out, '(A)') 'Masses for the two atoms should be specified with a line of the type: '
       write(out, '(A)') 'MASSES 12.00000000 15.99491461957'
       write(out, '(A)') 'Alternatively, specify the atoms with a line of the type: '
       write(out, '(A)') 'ATOMS C-12 O-16'
       stop
     endif
     !
     ! reduced mass in amu (i.e., Daltons)
     !
     amass = m1*m2/(m1+m2)
     !
     if (iverbose>=4) write(out, '(a,G25.15,5x,a)'  ) 'Reduced mass is ', amass,  'atomic mass units (Daltons)'
     if (iverbose>=4) write(out, '(a,G25.15,5x,a,/)') 'Reduced mass is ', amass*umatoau, 'atomic units (electron masses)'
     !
     scale = amass/aston
     !
     if (iverbose>=4) call TimerStart('Define grid')
     !
     call gridred(nsub,alpha,re,ngrid,rmin,rmax,h,r,z,f,iverbose)
     !
     grid%r = r
     !
     if (iverbose>=4) call TimerStop('Define grid')
     !
     hstep = h
     !
     h12 = 12.0_rk*h**2
     sc  = h12*scale
     !
     ! For uniformly spaced grid z(j)= 1 ; f(j) = 0
     do j = 1, ngrid
        r2sc(j)= scale*r(j)*r(j)
        d2dr(j) = 30.0_rk*z(j)*z(j) - h12*f(j)
     enddo
     !
     deallocate(f)
     call ArrayStop('f-field')
     !
     ! generate grid-representaion for the all type of of the hamiltonian fields
     !
     if (iverbose>=4) call TimerStart('Grid representaions')
     !
     if (iverbose>=3) write(out,'("Generate a grid representation for all Hamiltonian terms")')
     !
     ! in case of fitting we will need this object to store the parameters to refine
     !
     if (action%fitting) then
       !
       itotal = 0
       ifterm = 0
       !
     endif
     !
     object_loop: do iobject = 1,Nobjects
        !
        Nmax = fieldmap(iobject)%Nfields
        !
        ! each field type consists of Nmax terms
        !
        do iterm = 1,Nmax
          !
          select case (iobject)
          case (1)
            field => poten(iterm)
          case (2)
            field => spinorbit(iterm)
          case (3)
            field => l2(iterm)
          case (4)
            field => lxly(iterm)
          case (5)
            field => spinspin(iterm)
          case (6)
            field => spinspino(iterm)
          case (7)
            field => bobrot(iterm)
          case (8)
            field => spinrot(iterm)
          case (9)
            field => diabatic(iterm)
          case (10)
            field => lambdaopq(iterm)
          case (11)
            field => lambdap2q(iterm)
          case (12)
            field => lambdaq(iterm)
          case (13)
            field => nac(iterm)
          case (21, 22, 23, 24, 25, 26, 27)
            field => hfcc1(iobject - 20)%field(iterm)
          case (Nobjects-3)
            field => quadrupoletm(iterm)
          case (Nobjects-2)
            field => abinitio(iterm)
          case (Nobjects-1)
            cycle
          case (Nobjects)
            field => dipoletm(iterm)
          case default
             print "(a,i0)", "iobject = ",iobject
             stop "illegal iobject  "
          end select
          !
          if (.not.gridvalue_allocated) then 
            !
            allocate(field%gridvalue(ngrid),stat=alloc)
            call ArrayStart(trim(field%type),alloc,ngrid,kind(field%gridvalue))
            !
          endif
          !
          select case(trim(field%type))
          !
          case("GRID")
            !
            nterms = field%Nterms
            !
            ! Lorenzo Lodi, 13 February 2014
            ! This section will add extrapolated points at short and long bond lengths
            ! for tabulated `GRID' functions
            ! and only for fitting 
            if( field%grid(1) > rmin ) then
               if (iverbose>=4) write(out, '(/,A)') 'Extrapolating at short bond length curve ' // trim(field%name) // &
                                 ' (class ' // trim(field%class) // ')'
         
         
                 np = 20              ! I always add `np' short bond length guide points
                 ! I go one step `beyond' rmin to minimize interpolating artifacts at the end of the interval
         
                 x1=field%grid(1)  ; y1 =field%value(1)
                 x2=field%grid(2)  ; y2 =field%value(2)
                 allocate( xx(nterms+np), yy(nterms+np),ww(nterms+np),stat=alloc)
                 call ArrayStart('extrap_tmp',alloc,size(xx),kind(xx))
                 call ArrayStart('extrap_tmp',alloc,size(ww),kind(ww))
                 call ArrayStart('extrap_tmp',alloc,size(yy),kind(yy))
                 allocate( link_(nterms+np),stat=alloc)
                 if (alloc/=0) then 
                   write(out,"('allocation problem: extrap_tmp')")
                   stop 'allocation problem'
                 endif
                 xx=0._rk
                 yy=0._rk
                 do i=1, np
                    xx(i) = rmin + (field%grid(1)-rmin)*real(i-2,rk) / real( np-1, rk)
                 enddo
                 !
               select case(field%class)
                 !
               case ("POTEN")
                 !
                 if (iverbose>=4) write(out, '(A, I20)') 'Using A + B/r extrapolation; points added = ', np
                 bb = -x1*x2*(y1-y2)/(x1-x2)
                 aa = y1 - bb/x1
                 do i=1, np
                    yy(i) = aa + bb/xx(i)
                 enddo
                 !
               case ("DIPOLE")
                 if (iverbose>=4) write(out, '(A, I20)') 'Using A*r + B*r^2 extrapolation; points added = ', np
                 bb = (x2*y1- x1*y2) / (x1*x2 *(x1-x2) )
                 aa = (-bb*x1**2 + y1 ) / x1
                 do i=1, np
                    yy(i) = aa*xx(i) + bb*xx(i)**2
                 enddo
                 !
               case default  ! linear extrapolation using the first two points to r = rmin - 1/(np-1)
                 if (iverbose>=4) write(out, '(A, I20)') 'Using linear extrapolation; points added = ', np
                 do i=1, np
                      yy(i) = y1 + (xx(i)-x1) * (y1-y2)/(x1-x2)
                 enddo
               end select
               !
               do i=np+1, np+nterms
                  xx(i) =  field%grid(i-np)
                  yy(i) =  field%value(i-np)
               enddo
               !
               ww(np+1:) = field%weight(1:) ; ww(1:np) = 0
               !
               link_(1:np)%iobject = 0
               link_(1:np)%ifield = 0
               link_(1:np)%iparam = 0
               link_(np+1:) = field%link(1:)
               !
               nterms = np+nterms
               field%Nterms = nterms
               !
               call ArrayMinus(trim(field%type),size(field%value),kind(field%value))
               call ArrayMinus(trim(field%type),size(field%grid),kind(field%grid))
               call ArrayMinus(trim(field%type),size(field%weight),kind(field%weight))
               !
               deallocate(field%grid, field%value, field%weight, field%forcename,stat=alloc)
               !
               deallocate(field%link)
               !
               allocate(field%value(nterms),field%forcename(nterms),field%grid(nterms),field%weight(nterms),stat=alloc)
               allocate(field%link(nterms),stat=alloc)
               !
               call ArrayStart(trim(field%type),alloc,nterms,kind(field%value))
               call ArrayStart(trim(field%type),alloc,nterms,kind(field%grid))
               call ArrayStart(trim(field%type),alloc,nterms,kind(field%weight))
               !
               field%link(:) = link_(:)
               field%grid  = xx
               field%value = yy
               field%forcename= 'dummy'
               field%weight = ww
               deallocate(xx, yy, ww,link_)
               call ArrayStop('extrap_tmp')
               !
            endif
            !
            !*************end of short bond length extrapolation **********************************************
            nterms = field%Nterms
            if( field%grid(nterms) < rmax) then
               if (iverbose>=4) write(out, '(/,A)') 'Extrapolating at long bond length curve ' // trim(field%name) // &
                                ' (class ' // trim(field%class) // ')'
               !
               np = 20    ! I always add `np' long bond length guide points
               ! I go one step `beyond' rmax to minimize interpolating artifacts at the end of the interval
               !
               x1=field%grid(nterms-1)  ; y1 =field%value(nterms-1)
               x2=field%grid(nterms)    ; y2 =field%value(nterms)
               allocate( xx(nterms+np), yy(nterms+np),ww(nterms+np),stat=alloc)
               call ArrayStart('extrap_tmp',alloc,size(xx),kind(xx))
               call ArrayStart('extrap_tmp',alloc,size(ww),kind(ww))
               call ArrayStart('extrap_tmp',alloc,size(yy),kind(yy))
               allocate( link_(nterms+np),stat=alloc)
               if (alloc/=0) then 
                 write(out,"('allocation problem: extrap_tmp')")
                 stop 'allocation problem'
               endif
               !
               xx=0._rk
               yy=0._rk
               do i=1, nterms
                  xx(i) =  field%grid(i)
                  yy(i) =  field%value(i)
               enddo
               !
               do i=nterms+1, nterms+np
                  xx(i) = field%grid(nterms) + (rmax-field%grid(nterms))*real(i-nterms,rk) / real( np-1, rk)
               enddo
               !
               select case(field%class)
         
               case ("POTEN")
                 if (iverbose>=4) write(out, '(A, I20)') 'Using De + A/r^6 extrapolation; points added = ', np
                 bb = (x1*x2)**6 * (y2-y1) / (x1**6 - x2**6)
                 aa = y1 - bb/x1**6
                 !       write(out, '(A, 2F20.4)') 'Dissociation set to = ', aa !, bb
                 do i=nterms+1,nterms+np
                    yy(i) = aa + bb / xx(i)**6
                 enddo
                 !
               case ("DIPOLE")
                 if (iverbose>=4) write(out, '(A, I20)') 'Using A/r^2 + B/r^3 extrapolation; points added = ', np
                 bb = x1*x2*(x1**2 * y1 - x2**2 * y2)/(x2-x1)
                 aa = (y1*x1**3 - bb)/x1
                 do i=nterms+1,nterms+np
                    yy(i) = aa/xx(i)**2 + bb / xx(i)**3
                 enddo
                 !
               case default  ! extrapolation using the last two points
                 if (iverbose>=4) write(out, '(A, I20)') 'Using linear extrapolation; points added = ', np
                 do i=nterms+1,nterms+np
                      yy(i) = y1 + (xx(i)-x1) * (y1-y2)/(x1-x2)
                  enddo
               end select
               !
               ww(1:nterms) = field%weight(1:nterms) ; ww(nterms+1:) = 0
               !
               link_(1:np)%iobject = 0
               link_(1:np)%ifield = 0
               link_(1:np)%iparam = 0
               link_(np+1:) = field%link(1:)
               !
               nterms = np+nterms
               field%Nterms = nterms
               !
               call ArrayMinus(trim(field%type),size(field%value),kind(field%value))
               call ArrayMinus(trim(field%type),size(field%grid),kind(field%grid))
               call ArrayMinus(trim(field%type),size(field%weight),kind(field%weight))
               !
               deallocate(field%grid, field%value, field%weight, field%forcename,stat=alloc)
               !
               deallocate(field%link)
               !
               allocate(field%value(nterms),field%forcename(nterms),field%grid(nterms),field%weight(nterms),stat=alloc)
               allocate(field%link(nterms),stat=alloc)
               !
               call ArrayStart(trim(field%type),alloc,nterms,kind(field%value))
               call ArrayStart(trim(field%type),alloc,nterms,kind(field%grid))
               call ArrayStart(trim(field%type),alloc,nterms,kind(field%weight))
               !
               field%link(:) = link_(:)
               field%grid  = xx
               field%value = yy
               field%forcename= 'dummy'
               field%weight = ww
               deallocate(xx, yy, ww,link_)
               call ArrayStop('extrap_tmp')
               !
            endif
            !
            !*************end of long bond length extrapolation **********************************************
            !
            ! build spline interpolant
            select case (field%interpolation_type)
              !
            case default
               write(out,"(a)") "Unrecognized interpolation type ", field%interpolation_type
               stop "illegal iobject"
               ! -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- .
            case("CUBICSPLINES")  ! cubic natural spline
               !
               if (iverbose>=6) write(out, '(A, I20)') 'Interpolating with cubic natural splines'
               allocate(spline_wk_vec(nterms),stat=alloc)
               call ArrayStart('spline_wk_vec-field',alloc,ngrid,kind(spline_wk_vec))
               !
               yp1= 0._rk ; ypn =0._rk  ! 1nd derivatives at the first and last point (ignored)
               call spline(field%grid,field%value,field%Nterms,yp1,ypn,spline_wk_vec)
               !
               !$omp parallel do private(i) schedule(guided)
               do i=1,ngrid
                ! evaluate spline interpolant
                call splint(field%grid,field%value,spline_wk_vec,field%Nterms,r(i),field%gridvalue(i))
               enddo
               !$omp end parallel do
               !
               deallocate(spline_wk_vec)
               call ArrayStop('spline_wk_vec-field')
               !
               ! -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- .
            case("QUINTICSPLINES")  ! quintic spline
               ! Lorenzo Lodi, 7 April 2016
               ! Tests show that the present quintic spline interpolation does not give better interpolation results
               ! than cubic splines. Furthermore the error seems to decrease as npoints**(-4)
               ! instead of npoints**(-6). The reason for this behaviour is yet unknown, but possibly related 
               ! to a bad choice of the boundary conditions. 
               !
               ! FOR NOW USE OF QUINTIC SPLINE INTERPOLATION IS NOT RECOMMENDED.
               ! USE CUBIC SPLINES INSTEAD
               !
               if (iverbose>=6) write(out, '(A, I20)') 'Interpolating with quintic splines'
               nterms=field%Nterms
               allocate(spline_wk_vec_B(nterms),stat=alloc); call ArrayStart('spline_wk_vec_B',alloc,nterms,kind(spline_wk_vec))
               allocate(spline_wk_vec_C(nterms),stat=alloc); call ArrayStart('spline_wk_vec_C',alloc,nterms,kind(spline_wk_vec))
               allocate(spline_wk_vec_D(nterms),stat=alloc); call ArrayStart('spline_wk_vec_D',alloc,nterms,kind(spline_wk_vec))
               allocate(spline_wk_vec_E(nterms),stat=alloc); call ArrayStart('spline_wk_vec_E',alloc,nterms,kind(spline_wk_vec))
               allocate(spline_wk_vec_F(nterms),stat=alloc); call ArrayStart('spline_wk_vec_F',alloc,nterms,kind(spline_wk_vec))
               !
               call QUINAT(nterms,field%grid,field%value, spline_wk_vec_B, spline_wk_vec_C, &
                                            spline_wk_vec_D, spline_wk_vec_E, spline_wk_vec_F)
               !
               !$omp parallel do private(i) schedule(guided)
               do i=1,ngrid    ! evaluate spline interpolant
!                 call splint(field%grid,field%value,spline_wk_vec,field%Nterms,r(i),field%gridvalue(i))

                  call splint_quint(field%grid,field%value,nterms,r(i),field%gridvalue(i), spline_wk_vec_B, &
                                     spline_wk_vec_C, spline_wk_vec_D, spline_wk_vec_E, spline_wk_vec_F)
               enddo
               !$omp end parallel do
               !
               deallocate(spline_wk_vec_B) ; call ArrayStop('spline_wk_vec_B')
               deallocate(spline_wk_vec_C) ; call ArrayStop('spline_wk_vec_C')
               deallocate(spline_wk_vec_D) ; call ArrayStop('spline_wk_vec_D')
               deallocate(spline_wk_vec_E) ; call ArrayStop('spline_wk_vec_E')
               deallocate(spline_wk_vec_F) ; call ArrayStop('spline_wk_vec_F')
               ! -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- .
               !
            end select
            !
            ! for dummy fields not used in fittings
            !
          case("DUMMY")
            !
            nterms = field%Nterms
            field%gridvalue = 0._rk
            !
          case default
            !
            call define_analytical_field(field%type,field%analytical_field)
            !
            !$omp parallel do private(i) schedule(guided)
            do i=1,ngrid
              !
              field%gridvalue(i) = field%analytical_field(r(i),field%value)
              !
            enddo
            !$omp end parallel do
            !
            ! counting the total number of parameters when the fitting is on
            !
            if (action%fitting) then
              itotal = itotal + field%Nterms
            endif
            !
          end select
          !
          field%gridvalue =  field%gridvalue*field%factor
          !
        enddo
        !
     enddo object_loop
     !
     ! change the status of  gridvalue_allocated to true (allocated)
     ! to prevent reallocation next time this subroutine is called (e.g. from the refinement)
     gridvalue_allocated = .true.
     !
     ! Now morph the objects by applying the morphing function to the corresponding ab initio field if neccessary
     !
     ifield = 0 
     
     object_loop2: do iobject = 1,Nobjects
        !
        Nmax = fieldmap(iobject)%Nfields
        !
        ! each field type consists of Nmax terms
        !
        do iterm = 1,Nmax
          !
          select case (iobject)
          case (1)
            field => poten(iterm)
          case (2)
            field => spinorbit(iterm)
          case (3)
            field => l2(iterm)
          case (4)
            field => lxly(iterm)
          case (5)
            field => spinspin(iterm)
          case (6)
            field => spinspino(iterm)
          case (7)
            field => bobrot(iterm)
          case (8)
            field => spinrot(iterm)
          case (9)
            field => diabatic(iterm)
          case (10)
            field => lambdaopq(iterm)
          case (11)
            field => lambdap2q(iterm)
          case (12)
            field => lambdaq(iterm)
          case (13)
            field => nac(iterm)
          case (21, 22, 23, 24, 25, 26, 27)
            field => hfcc1(iobject - 20)%field(iterm)
          case (Nobjects-3)
            field => quadrupoletm(iterm)
          case (Nobjects-2)
            field => abinitio(iterm)
          case (Nobjects-1)
            !
            ! no morphing for Brot 
            cycle
            !
          case (Nobjects)
            field => dipoletm(iterm)
          case default
             print "(a,i0)", "iobject = ",iobject
             stop "illegal iobject  "
          end select
          !
          ifield = ifield + 1
          !
          ! Shift values up or down by specified constant value
          if (field%adjust) then
            field%gridvalue = field%gridvalue + field%adjust_val
          endif
          !
          ! Introduce morphing (not for abinitio)
          !
          if (field%morphing.and.iobject/=Nobjects-2) then
            !
            ! check if ai field was defined 
            !
            check_ai = sum((abinitio(ifield)%gridvalue)**2)
            !
            if (check_ai<small_) then 
              !
              write(out,"('Error: Corresponding ab initio field is undefined or zero when using MORPHING ',a)") trim(field%name)
              stop 'ab initio field is undefined while using MORPHING'
              !
            endif
            ! 
            field%gridvalue = field%gridvalue*abinitio(ifield)%gridvalue
            !
          endif
          !
          ! Generate weights if an analytical expression is given
          !
          if (trim(field%weighting%wtype)=="PS1997") then 
            !
            istate = field%iref
            beta = field%weighting%alpha
            Vtop = field%weighting%Vtop
            imin = (field%grid(1)-rmin) / real( np-1, rk)
            Vmin = minval(poten(istate)%gridvalue(:))
            !
            DeltaR = (rmax-rmin)/real(ngrid-1,rk)
            !
            do i=1,field%Nterms
              !
              np = nint(( field%grid(i) - rmin )/ DeltaR)
              !
              np = max(min(ngrid,np),1)
              !
              field%weight(i) = ( tanh(-beta*( poten(istate)%gridvalue(np) - Vmin - Vtop ) ) & 
                                   +1.000020000200002_rk )/2.000020000200002_rk
              !
            enddo
            !
          endif  
          !
          !
          ! transform from the MOLRPO to Duo representation
          !
          if (field%molpro) then
            !
            if (iverbose>=4) then 
                write(out,'(/,a)') 'Transforming '//trim(field%class)//' '//trim(field%name)//&
                ' from MOLPRO (Cartesian) to Duo (Spherical)'            
                write(out,'(a)') '  Cartesian (complex):'
                write(out,'(a,t20,i4,1x,t30,i4)') '  |Lambda|:',field%lambda,field%lambdaj
                if (iobject==2) write(out,'(a,t18,f8.1,t28,f8.1)') '  Sigma:',field%sigmai,field%sigmaj 
                write(out,'(a,t20,i4,a,t30,i4,a)') '  <a|Lz|b>:',field%ix_lz_y,'i',field%jx_lz_y,'i'
                write(out,'(a,t18,f8.1,a,f8.1,a)') '  factor',real(field%complex_f),'+',aimag(field%complex_f),'i'
            endif
            !
            if (job%legacy_version) then  
              !
              call molpro_duo_old_2018(field)
              !
            else 
              !
              call molpro_duo(field)
              !
            endif
            !
            if (iverbose>=4) then 
                write(out,'(a)') '  Spherical (real):'
                write(out,'(a,t20,i4,1x,t30,i4)') '  Lambda:',field%lambda,field%lambdaj
                if (iobject==2) write(out,'(a,t18,f8.1,t28,f8.1)') '  Sigma:',field%sigmai,field%sigmaj 
            endif
          endif
          !
        enddo
        !
     enddo object_loop2
     !
     if (action%fitting) then
       fitting%parmax = itotal
     endif
     !
     if (iverbose>=3) write(out,'("...done!"/)')
     !
     if (iverbose>=4) call TimerStop('Grid representaions')
     !
     !
     ! Lorenzo Lodi, 8 April 2015
     ! Find the minima of all potential functions (real mimima, not on the grid) and compute other
     ! equilibrium quantities
     !
     loop_pecs: do istate=1,Nestates ! loop over potential energy curves
       ! find minimum and maximum on the grid 
       poten(istate)%imin  = minloc(poten(istate)%gridvalue,dim=1)
       poten(istate)%Vimin = poten(istate)%gridvalue( poten(istate)%imin )
       poten(istate)%imax  = maxloc(poten(istate)%gridvalue,dim=1)
       poten(istate)%Vimax = poten(istate)%gridvalue( poten(istate)%imax )
       poten(istate)%rimin = r(poten(istate)%imin)
       ! if minimum is at extremes of the grid and/or minimum==maximum then the potential curve has no minimum
       if( poten(istate)%imin ==1                            ) poten(istate)%zHasMinimum = .false.
       if( poten(istate)%imin ==size(poten(istate)%gridvalue)) poten(istate)%zHasMinimum = .false.
       ! if minimum is same as maximum potential is constant
       if( poten(istate)%Vimin == poten(istate)%Vimax        ) poten(istate)%zHasMinimum = .false.

       poten(istate)%V0  = poten(istate)%Vimin ! for now set the minimum to the grid-minimum
       !
       if( poten(istate)%zHasMinimum .eqv. .false. ) cycle loop_pecs ! exit if pec has no minimum

       ! For not if we have a spline interpolant we don't go on... (to be fixed)
       if( poten(istate)%type == 'GRID') then
          !
          if( poten(istate)%imin-4<1 .or. poten(istate)%imin+4>grid%npoints) cycle loop_pecs ! exit if the minimum is too close to the border
          !
          x0  = r(poten(istate)%imin)
          fmmm = poten(istate)%gridvalue(poten(istate)%imin-3)
          fmm  = poten(istate)%gridvalue(poten(istate)%imin-2)
          fm   = poten(istate)%gridvalue(poten(istate)%imin-1)
          f0   = poten(istate)%gridvalue(poten(istate)%imin)
          fp   = poten(istate)%gridvalue(poten(istate)%imin+1)
          fpp  = poten(istate)%gridvalue(poten(istate)%imin+2)
          fppp = poten(istate)%gridvalue(poten(istate)%imin+3)
      
          der1 = (-fmmm + 9.0_rk*fmm-45.0_rk*fm+45._rk*fp-9.0_rk*fpp+fppp) /(60._rk*h) ! 6-point, error h^6
          der2 = (2.0_rk*(fmmm+fppp)-27._rk*(fmm+fpp)+270._rk*(fm+fp)-490._rk*f0 )/(2.0_rk*90._rk*h**2) ! 7-point, error h^6
      
          !  der1 = (fmm-8._rk*fm+8._rk*fp-fpp) /(6._rk*h)                          ! 4-point, error h^4
          !  der2 = (-fmm +16_rk*fm -30._rk*f0 +16._rk*fp - fpp ) / (3._rk * h**2)  ! 5-point, error h^4
          der3 = (fmmm-8._rk*fmm+13._rk*fm -13._rk*fp+8._rk*fpp-fppp) / (h**3*8.0_rk)        ! 6-point, error h^4
          der4 = (-fmmm+12._rk*fmm-39._rk*fm+56._rk*f0-39._rk*fp+12._rk*fpp-fppp )/( 6._rk* h**4)   !7-point, error h^4
          !
       else 
          ! Find miminum of each PEC using Newton-type search
          ! Note: computing derivatives by finite differences is an intrinsically ill-conditioned problem
          ! because of cancellation errors. It gets worse with higher order, so
          ! we use a 6th order formula for the 1st derivative but lower order for higher
          ! derivatives. The step size `h' and the derivative formulas were chosen on the
          ! basis of tests using morse-type potentials typical for molecules.
          ! A value of the step size of about 2e-3 or larger should work best.
          !
          x0  = r(poten(istate)%imin)
          find_minimum: do i=1, max_iter_min_search
      
               fmmm = poten(istate)%analytical_field( x0-1.5_rk*h, poten(istate)%value )
               fmm  = poten(istate)%analytical_field( x0-h       , poten(istate)%value )
               fm   = poten(istate)%analytical_field( x0-h/2._rk , poten(istate)%value )
               f0   = poten(istate)%analytical_field( x0         , poten(istate)%value )
               fp   = poten(istate)%analytical_field( x0+h/2._rk , poten(istate)%value )
               fpp  = poten(istate)%analytical_field( x0+h       , poten(istate)%value )
               fppp = poten(istate)%analytical_field( x0+1.5_rk*h, poten(istate)%value )
      
               der1 = (-fmmm + 9._rk*fmm-45._rk*fm+45._rk*fp-9.0_rk*fpp+fppp) /(30._rk*h) ! 6-point, error h^6
               der2 = (2._rk*(fmmm+fppp)-27._rk*(fmm+fpp)+270._rk*(fm+fp)-490._rk*f0 )/(45._rk*h**2) ! 7-point, error h^6
      
               !  der1 = (fmm-8._rk*fm+8._rk*fp-fpp) /(6._rk*h)                          ! 4-point, error h^4
               !  der2 = (-fmm +16_rk*fm -30._rk*f0 +16._rk*fp - fpp ) / (3._rk * h**2)  ! 5-point, error h^4
               der3 = (fmmm-8._rk*fmm+13._rk*fm -13._rk*fp+8._rk*fpp-fppp) / h**3        ! 6-point, error h^4
               der4 = (-8._rk*fmmm+96._rk*fmm-312._rk*fm+448._rk*f0-312._rk*fp+ &
                                                96._rk*fpp-8._rk*fppp )/( 3._rk* h**4)   !7-point, error h^4
                !  der1 = ( fp - fm ) / h                                    ! 2-point, error h^2
                !  der2 = 4._rk*( fm - 2_rk*f0 + fp ) / h**2                 ! 3-point, error h^2
                !  der3 = 4._rk*(-fmm + 2._rk*fm -2._rk*fp + fpp ) / h**3    ! 4-point, error h^2
                !  der4 = 16._rk*(fmm-4._rk*fm+6._rk*f0-4._rk*fp+fpp) / h**4 ! 5-point, error h^2
      
                if(der2 ==0._rk)   der2 = tiny(1._rk) ! because we are going to divide by der2 
                x1 = x0 - der1/der2
      
               ! check convergence
               if( abs(x1 - x0) / max( abs(x0), tiny(1._rk) ) <= thresh_rel_err_min_x) exit find_minimum
               x0 = x1
      
          enddo find_minimum
          !
          if( i >= max_iter_min_search) poten(istate)%zHasMinimum = .false.
          !
       endif
       !
       if( der2 < 0._rk ) poten(istate)%zHasMinimum = .false.
       poten(istate)%re  = x0
       poten(istate)%V0  = f0
       poten(istate)%V1  = der1
       poten(istate)%V2  = der2
       poten(istate)%V3  = der3
       poten(istate)%V4  = der4
       ! From formulas 4.67 to 4.72 (pag. 155-156) of I. N. Levine, Molecular Spectroscopy, John Wiley & Sons (1975)
       poten(istate)%we  = bohr*sqrt(abs(der2)*hartree / (amass*umatoau) )
       poten(istate)%B0  = 0.5_rk*bohr**2*hartree /( amass*umatoau * x0**2 )
       poten(istate)%xe  = (  (poten(istate)%B0**2*x0**4) / (12._rk*poten(istate)%we**5) )  & 
                             * (10._rk*poten(istate)%B0*der3**2*x0**2-3._rk*der4*poten(istate)%we**2)
       poten(istate)%alphae  = -( 2._rk*poten(istate)%B0**2 / ( poten(istate)%we ) ) &
                               * (3._rk + (2._rk*poten(istate)%B0*der3*x0**3)/poten(istate)%we**2 )
       poten(istate)%Debar  = (4._rk*poten(istate)%B0**3) / ( poten(istate)%we**2 )
       poten(istate)%Y00  = ( 9._rk*der4*poten(istate)%we**2 - 14._rk*poten(istate)%B0*der3**2*x0**2 ) &
                               *((poten(istate)%B0**2 * x0**4) / ( 144._rk*poten(istate)%we**4) )

       ! This is an approximation to the energy of the v=0, J=0 level assuming Spin=0,
       !  Lambda=0 (1 Sigma term) and neglecting all couplings (spin-orbit, non-adiabatic,...).
       poten(istate)%approxEJ0 = poten(istate)%Y00 + poten(istate)%V1 &
                                                   + poten(istate)%we*0.5_rk &
                                                   - poten(istate)%we*poten(istate)%xe*0.5_rk**2
       ! This is an approximation to the energy of the physically-possible v=0
       ! level, still desregarding couplings to other curves (spin-orbit etc).
       ! NB: TO BE FIXED! the expression below is incorrect
       poten(istate)%approxEJmin = poten(istate)%approxEJ0 - poten(istate)%alphae*0.5_rk*poten(istate)%Omega_min &
                                          - poten(istate)%Debar* ( poten(istate)%Omega_min*(poten(istate)%Omega_min+1) )**2

     enddo loop_pecs

     ipotmin    = 1 ! find state with lowest-lying pec
     job%potmin = poten(ipotmin)%V0
     do istate=1,Nestates
        if( poten(istate)%V0 <= job%potmin ) then
            ipotmin = istate
            job%potmin = poten(istate)%V0
        endif
     enddo
     job%potmin = poten(ipotmin)%V0
     !
     !
     if(iverbose >= 6 )  write(out, '(A, F30.6, A)')   'The lowest lying energy point has energy ', job%potmin, ' cm^-1'
     if( job%zShiftPECsToZero ) then ! shift all PECs relative to the lowest point in all electronic curves
      if(iverbose >= 6)  write(out, '(A, F30.6, A)') 'Shifting all potential energy curves by  ', job%potmin, ' cm^-1'
       do istate=1,Nestates
          poten(istate)%gridvalue(:) = poten(istate)%gridvalue(:)-job%potmin
       enddo
     else
       if(iverbose >= 6) write(out, '(A, F30.6, A)') 'Potential energy curves not shifted.'
     endif
     !
     call check_and_print_field(Nestates,iverbose,poten,"Potential functions:")
     !
     call check_and_print_coupling(Nspinorbits,iverbose,spinorbit,"Spin-Orbit:")
     call check_and_print_coupling(Nlxly,      iverbose,lxly,     "<L+> functions:")
     call check_and_print_coupling(NL2,        iverbose,l2,       "<L**2> functions:")
     call check_and_print_coupling(Nss,        iverbose,spinspin, "Spin-spin functions:")
     call check_and_print_coupling(Nsso,       iverbose,spinspino,"Spin-spin-o (non-diagonal) functions:")
     call check_and_print_coupling(Nsr,        iverbose,spinrot,  "Spin-rotation functions:")
     call check_and_print_coupling(Nbobrot,    iverbose,bobrot,   "Bob-Rot centrifugal functions:")
     call check_and_print_coupling(Ndiabatic,  iverbose,diabatic, "Diabatic functions:")
     call check_and_print_coupling(Nlambdaopq, iverbose,lambdaopq,"Lambda-opq:")
     call check_and_print_coupling(Nlambdap2q, iverbose,lambdap2q,"Lambda-p2q:")
     call check_and_print_coupling(Nlambdaq,   iverbose,lambdaq,  "Lambda-q:")
     call check_and_print_coupling(Nnac,       iverbose,nac,  "NACouplings")
     if(associated(dipoletm)) call check_and_print_coupling(Ndipoles,   iverbose,dipoletm, "Dipole moment functions:")
     if(associated(quadrupoletm)) call check_and_print_coupling(nQuadrupoles,iverbose,quadrupoletm, "Quadrupole moment functions:")
     !
   contains 
     !
     subroutine molpro_duo(field)
        !
        use lapack,only : lapack_heev     
        !
        type(fieldT),intent(inout) :: field
        integer(ik) :: ix_lz_y,jx_lz_y,iroot,jroot,il_temp,ngrid
        complex(rk) :: a(2,2),b(2,2),coupling(2,2),f_t(2,2),b0(2,2),c
        real(rk)    :: lambda_i(2),lambda_j(2)
            !
            ngrid = grid%npoints
            !
            ix_lz_y = field%ix_lz_y
            jx_lz_y = field%jx_lz_y
            !
            a = 0 ; b = 0
            !
            a(1,1) =  cmplx(1.0_rk,0.0_rk, kind=rk)
            b(1,1) =  cmplx(1.0_rk,0.0_rk, kind=rk)
            a(2,2) =  cmplx(1.0_rk,0.0_rk, kind=rk)
            b(2,2) =  cmplx(1.0_rk,0.0_rk, kind=rk)
            !
            iroot = 1
            jroot = 1
            !
            if (ix_lz_y/=0) then
              a = 0 
              a(1,2) = cmplx(0.0_rk,ix_lz_y , kind=rk)
              a(2,1) = cmplx(0.0_rk,-ix_lz_y, kind=rk)
              !
              call lapack_heev(a,lambda_i)
              !
              ! swap to have the first root positive 
              !
              f_t = a
              a(:,1) = f_t(:,2)
              a(:,2) = f_t(:,1)
              !
              il_temp = lambda_i(2)
              lambda_i(2) = lambda_i(1)
              lambda_i(1) = il_temp
              !
              field%lambda = nint(lambda_i(1))
              !
              a = a*cmplx(0.0_rk,1.0_rk, kind=rk)
              !
            elseif (poten(field%istate)%parity%pm==-1) then 
              !
              a(1,1) =  cmplx(0.0_rk,1.0_rk, kind=rk)
              !
            endif
            !
            if (jx_lz_y/=0) then
              b = 0 
              b(1,2) = cmplx(0.0_rk,jx_lz_y , kind=rk)
              b(2,1) = cmplx(0.0_rk,-jx_lz_y, kind=rk)
              !
              b0 = b
              !
              call lapack_heev(b,lambda_j)
              !
              ! swap to have the first root positive 
              !
              f_t = b
              b(:,1) = f_t(:,2)
              b(:,2) = f_t(:,1)
              !
              il_temp = lambda_j(2)
              lambda_j(2) = lambda_j(1)
              lambda_j(1) = il_temp
              !
              field%lambdaj = nint(lambda_j(1))
              !
              b = b*cmplx(0.0_rk,1.0_rk, kind=rk)
              !
              f_t = matmul( conjg(transpose(b)),matmul(b0,(b)) )
              !
            elseif (poten(field%jstate)%parity%pm==-1) then 
              !
              b(1,1) =  cmplx(0.0_rk,1.0_rk, kind=rk)
              !
            endif
            !
            ! Check the selection rules
            select case(trim(field%class))
              !
            case('SPINORBIT','ABINITIO-SPINORBIT')
              !
              if ((nint(field%sigmaj-field%sigmai))/=(field%lambda-field%lambdaj)) then
                !
                ! try to select the root#2 for the i-state first 
                !
                if (field%lambda/=0.and.(nint(field%sigmaj-field%sigmai))==( lambda_i(2)-field%lambdaj ) ) then
                  !
                  il_temp = lambda_i(2)
                  lambda_i(2) = lambda_i(1)
                  lambda_i(1) = il_temp
                  !
                  field%lambda = nint(lambda_i(1))
                  !
                  f_t = a
                  a(:,1) = f_t(:,2)
                  a(:,2) = f_t(:,1)
                  !
                elseif( field%lambdaj/=0.and.(nint(field%sigmaj-field%sigmai))==( field%lambda-lambda_j(2) ) ) then
                  !
                  il_temp = lambda_j(2)
                  lambda_j(2) = lambda_j(1)
                  lambda_j(1) = il_temp
                  !
                  jroot = 2
                  field%lambdaj = nint(lambda_j(1))
                  !
                  f_t = b
                  b(:,1) = f_t(:,2)
                  b(:,2) = f_t(:,1)
                  !
                elseif( field%lambdaj/=0.and.(nint(field%sigmaj-field%sigmai))==( lambda_i(2)-lambda_j(2) ) ) then
                  !
                  il_temp = lambda_i(2)
                  lambda_i(2) = lambda_i(1)
                  lambda_i(1) = il_temp
                  !
                  field%lambda = nint(lambda_i(1))
                  !
                  il_temp = lambda_j(2)
                  lambda_j(2) = lambda_j(1)
                  lambda_j(1) = il_temp
                  !
                  jroot = 2
                  field%lambdaj = nint(lambda_j(1))
                  !
                  f_t = a
                  a(:,1) = f_t(:,2)
                  a(:,2) = f_t(:,1)
                  !
                  f_t = b
                  b(:,1) = f_t(:,2)
                  b(:,2) = f_t(:,1)
                  !
                else
                  !
                  write(out,"(/'molpro_duo: cannot find selecion rule for',2i8,3x,a)") field%iref,field%jref,trim(field%class)
                  write(out,"(' sigma = ',2f8.1,' lambda (i) = ',i4,' lamda (j) = ',i4)") field%sigmai, field%sigmaj, &
                                                                                                  int(lambda_i(1)),int(lambda_j(1))
                  stop 'molpro_duo: cannot find the selecion rules'
                endif
                !
              endif
              ! 
            case ('L+','ABINITIO-LX')
              !
              !write(out,"('molpro_duo: this L+-part is not implemented')")
              !stop 'molpro_duo: this L+-part is not implemented'
              ! 
            case ('DIPOLE')
              !
              !write(out,"('molpro_duo: this Dipole-part is not implemented')")
              !stop 'molpro_duo: this Dipole-part is not implemented'
              !
            case default
              !
              write(out,"(/'molpro_duo: this part is not implemented:',a)") trim(field%class)
              stop 'molpro_duo: this part is not implemented'
              !
            end select
            !
            !omp parallel do private(i) schedule(guided)
            do i=1,ngrid
              !
              coupling = 0 
              if (ix_lz_y==0.and.jx_lz_y==0) then
                !
                coupling(1,1) = field%gridvalue(i)*field%complex_f 
                !
              elseif(ix_lz_y/=0.and.jx_lz_y==0) then
                !
                select case(trim(field%class))
                  !
                case('SPINORBIT','ABINITIO-SPINORBIT')
                  !
                  if (abs(nint(field%sigmaj-field%sigmai))/=abs(field%lambda-field%lambdaj)) then
                    !
                    write(out,"(/'molpro_duo: SO ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") &
                         field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                    stop 'molpro_duo: illegal selection rules '
                    !
                  endif
                  !
                  ! for SOX it is either <1.2|SOX (regular) or <1.3|SOY (iregular)
                  !
                  if (poten(field%jstate)%parity%pm==-1) then 
                    !
                    ! regular
                    coupling(1,1) = field%gridvalue(i)*field%complex_f
                    coupling(2,1) =-coupling(1,1)*conjg(a(1,2)/a(2,2))
                    !
                  else
                    !
                    ! irregular
                    coupling(2,1) = field%gridvalue(i)*field%complex_f
                    coupling(1,1) =-coupling(2,1)*conjg(a(2,2)/a(1,2))
                    !
                  endif
                  ! 
                case ('L+','ABINITIO-LX')
                  !
                  if (poten(field%jstate)%parity%pm==-1) then 
                    !
                    ! regular
                    coupling(1,1) = field%gridvalue(i)*field%complex_f
                    coupling(2,1) =-coupling(1,1)*conjg(a(1,2)/a(2,2))
                    !
                  else
                    !
                    ! iregular
                    coupling(2,1) = field%gridvalue(i)*field%complex_f
                    coupling(1,1) = -coupling(2,1)*conjg(a(2,2)/a(1,2))
                    !
                  endif
                  ! 
                case ('DIPOLE')
                  !
                  if (poten(field%jstate)%parity%pm==-1) then 
                    !
                    ! iregular
                    coupling(2,1) = field%gridvalue(i)*field%complex_f*(-sqrt(0.5_rk))
                    coupling(1,1) =-coupling(2,1)*conjg(a(2,2)/a(1,2))
                    !
                  else
                    !
                    ! regular
                    coupling(1,1) = field%gridvalue(i)*field%complex_f*(-sqrt(0.5_rk))
                    coupling(2,1) =-coupling(1,1)*conjg(a(1,2)/a(2,2))
                    !
                  endif
                  !
                case default
                  !
                  write(out,"(/'molpro_duo (lambdaj=0): for class ',a,' is not implemented ')") field%class
                  stop 'molpro_duo: not for this object'
                  !
                end select
                !
              elseif(ix_lz_y==0.and.jx_lz_y/=0) then
                !
                select case(trim(field%class)) 
                  !
                case('SPINORBIT','ABINITIO-SPINORBIT')
                  !
                  if (abs(nint(field%sigmaj-field%sigmai))/=abs(field%lambda-field%lambdaj)) then
                    !
                    write(out,"(/'molpro_duo: SO ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") & 
                          field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                    stop 'molpro_duo: illegal selection rules '
                    !
                  endif
                  !
                  ! The following relations are found using the condition <0,Sigma|SOX+SOY|+/-1,Sigma+/-1> = 0
                  ! Regualr case is for <0|SOX|Pix> and iregular is for <0|SOX|Piy>
                  !
                  if (poten(field%istate)%parity%pm==-1) then 
                    !
                    ! regular
                    coupling(1,1) = field%gridvalue(i)*field%complex_f
                    coupling(1,2) =-coupling(1,1)*b(1,2)/b(2,2)
                    !
                  else
                    !
                    ! iregular
                    coupling(1,2) = field%gridvalue(i)*field%complex_f
                    coupling(1,1) =-coupling(1,2)*b(2,2)/b(1,2)
                    !
                  endif
                  ! 
                case('L+','ABINITIO-LX')
                  !
                  ! The following relations are found using the condition <0|LX+iLY|+1> = 0
                  ! Regualr case is for <0|LX|Pix> and iregular is for <0|LX|Piy>
                  !
                  if (poten(field%istate)%parity%pm==-1) then 
                    !
                    ! regular
                    coupling(1,1) = field%gridvalue(i)*field%complex_f
                    coupling(1,2) =-coupling(1,1)*b(1,2)/b(2,2)
                    !
                  else
                    !
                    ! iregular
                    coupling(1,2) = field%gridvalue(i)*field%complex_f
                    coupling(1,1) =-coupling(1,2)*b(2,2)/b(1,2)
                    !
                  endif
                  ! 
                case ('DIPOLE')
                  !
                  ! The following relations are found using the condition <0|-mux/sqrt(2)+imuy/sqrt(2)|+1> = 0
                  ! Regualr case is for <0|muX|Pix> and iregular is for <0|mux|Piy>
                  !
                  if (poten(field%istate)%parity%pm==-1) then 
                    !
                    coupling(1,2) = field%gridvalue(i)*field%complex_f*(-sqrt(0.5_rk))
                    coupling(1,1) = -coupling(1,2)*b(2,2)/b(1,2)
                    !
                  else
                    !
                    coupling(1,1) = field%gridvalue(i)*field%complex_f*(-sqrt(0.5_rk))
                    coupling(1,2) = -coupling(1,1)*b(1,2)/b(2,2)
                    !
                  endif
                  !
                case default
                  !
                  write(out,"(/'molpro_duo (lambdai=0): for class ',a,' is not implemented ')") field%class
                  stop 'molpro_duo: not for this object'
                  !
                end select
                !
              elseif (abs(ix_lz_y)/=abs(jx_lz_y)) then
                !
                select case(trim(field%class))
                  !
                case('SPINORBIT','ABINITIO-SPINORBIT')
                  !
                  if (abs(nint(field%sigmaj-field%sigmai))/=abs(field%lambda-field%lambdaj)) then
                    !
                    write(out,"(/'molpro_duo: SO ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") &
                          field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                    stop 'molpro_duo: illegal selection rules '
                    !
                  endif
                  !
                  c = field%gridvalue(i)*field%complex_f
                  !
                  !coupling(1,1) =  0
                  !coupling(1,2) =  c
                  !coupling(2,1) =  -c*conjg(a(1,1))*b(2,2)/(conjg(a(2,1))*b(1,2))
                  !coupling(2,2) =  0
                  !
                  ! these cominations are found using the following three conditions 
                  ! <-Lamba,Sigma|SO|Lambda',Sigma'> = 0 
                  ! <Lamba,Sigma|SO|-Lambda',Sigma'> = 0 
                  ! <-Lamba,Sigma|SO|-Lambda',Sigma'> = 0 
                  ! where we assume 
                  ! c = <x|LSX|y> and  < Lamba,Sigma|SO|Lambda',Sigma'> /= 0 
                  !
                  coupling(1,1) =  -c*b(2,2)/b(1,2)
                  coupling(1,2) =   c
                  coupling(2,1) =   c*conjg(a(1,2)/a(2,2))*b(2,2)/b(1,2)
                  coupling(2,2) =  -c*conjg(a(1,2)/a(2,2))
                  !
                case('DIPOLE')
                  !
                  if (abs(field%lambda-field%lambdaj)/=1) then
                    !
                    write(out,"('molpro_duo: DIPOLE ',2i4,'; illegal selection rules for lambda = ',2i4,' not +/-1')") & 
                          field%iref,field%jref,field%lambda,field%lambdaj
                    stop 'molpro_duo: illegal selection rules for transition dipole'
                    !
                  endif
                  !
                  c = -field%gridvalue(i)*field%complex_f*sqrt(0.5_rk)
                  !
                  ! maple for Lx
                  !c = -conjugate(A[1,2])*a/conjugate(A[2,2]), b = -a*B[1,2]/B[2,2], 
                  !d = B[1,2]*conjugate(A[1,2])*a/(B[2,2]*conjugate(A[2,2]))
                  ! m+ = -1/sqrt(2)l+
                  !
                  ! these cominations are found using the following three conditions 
                  ! <-Lamba|r+|Lambda'> = 0 
                  ! <Lamba|r+|-Lambda'> = 0 
                  ! <-Lamba|r+|-Lambda'> = 0 
                  ! where we assume 
                  ! c = <x|r+|x> and < Lamba|L+|Lambda'> /= 0 
                  !
                  coupling(1,1) =  c
                  coupling(1,2) = -c*b(1,2)/b(2,2)
                  coupling(2,1) = -c*conjg(a(1,2)/a(2,2))
                  coupling(2,2) =  c*conjg(a(1,2)/a(2,2))*b(1,2)/b(2,2)
                  !
                  !coupling(1,1) =  -c*b(2,2)/b(1,2)
                  !coupling(1,2) =   c
                  !coupling(2,1) =   c*conjg(a(1,2)/a(2,2))*b(2,2)/b(1,2)
                  !coupling(2,2) =  -c*conjg(a(1,2)/a(2,2))
                  !
                case('L+','ABINITIO-LX')
                  !
                  if (abs(field%lambda-field%lambdaj)/=1) then
                    !
                    write(out,"('molpro_duo: L+ ',2i4,'; illegal selection rules for lambda = ',2i4,' not +/-1')") & 
                          field%iref,field%jref,field%lambda,field%lambdaj
                    stop 'molpro_duo: L+ illegal selection rules for transition dipole'
                    !
                  endif
                  !
                  ! The <x|Lx+iLy|y> element is <x|Lx|y>
                  !
                  c = field%gridvalue(i)*field%complex_f
                  !
                  ! maple:
                  !
                  !d = -b*conjugate(A[1,2])/conjugate(A[2,2]), c = B[2,2]*b*conjugate(A[1,2])/(B[1,2]*conjugate(A[2,2])), 
                  !a = -B[2,2]*b/B[1,2]
                  !
                  ! these cominations are found using the following three conditions 
                  ! <-Lamba|L+|Lambda'> = 0 
                  ! <Lamba|L+|-Lambda'> = 0 
                  ! <-Lamba|L+|-Lambda'> = 0 
                  ! where we assume 
                  ! c = <x|L+|y> and  < Lamba|L+|Lambda'> /= 0 
                  !
                  coupling(1,1) = -c*b(2,2)/b(1,2)
                  coupling(1,2) =  c
                  coupling(2,1) =  c*conjg(a(1,2)/a(2,2))*b(2,2)/b(1,2)
                  coupling(2,2) = -c*conjg(a(1,2)/a(2,2))
                  !
                case default
                  !
                  write(out,"('molpro_duo (lambdaj<>lambdai): for class ',a,' is not implemented ')") field%class
                  stop 'molpro_duo: not for this object'
                  !
                end select
                !
              else
                !
                select case(trim(field%class))
                  !
                case('SPINORBIT','ABINITIO-SPINORBIT')
                  !
                  if (nint(field%sigmaj-field%sigmai)/=0) then
                    !
                    write(out,"('molpro_duo: SOZ ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") & 
                          field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                    stop 'molpro_duo: SOZ illegal selection rules '
                    !
                  endif
                  !
                  if ( field%sigmai<0 ) then
                    !
                    !write(out,"('molpro_duo: SO ',2i4,'; illegal reference sigmai <0 ',2f8.1)") &
                    !      field%iref,field%jref,field%sigmai,field%sigmaj
                    !stop 'molpro_duo: illegal reference sigmai'
                    !
                  endif
                  !
                  c = field%gridvalue(i)*field%complex_f
                  !
                  !coupling(1,1) =  c
                  !coupling(1,2) =  0
                  !coupling(2,1) =  0
                  !coupling(2,2) = -c*b(1,2)/b(2,2)*conjg(a(1,1))/conjg(a(2,1))
                  !
                  coupling(1,1) =  0
                  coupling(1,2) =  c
                  coupling(2,1) = -conjg(a(1,1))*c*b(2,2)/(conjg(a(2,1))*b(1,2))
                  coupling(2,2) =  0
                  !
                case('DIPOLE')
                  !
                  if (abs(field%lambda-field%lambdaj)/=0) then
                    !
                    write(out,"('molpro_duo: DMZ ',2i4,'; illegal selection rules for lambda = ',2i4,' not 0')") &
                          field%iref,field%jref,field%lambda,field%lambdaj
                    stop 'molpro_duo: illegal selection rules for DMZ'
                    !
                  endif
                  !
                  c = field%gridvalue(i)*field%complex_f
                  !
                  !coupling(1,1) =  0
                  !coupling(1,2) =  c
                  !coupling(2,1) = -conjg(a(1,1))*c*b(2,2)/(conjg(a(2,1))*b(1,2))
                  !coupling(2,2) =  0
                  !
                  coupling(1,1) =  c
                  coupling(1,2) =  0
                  coupling(2,1) =  0
                  coupling(2,2) =  c
                  !
                case default
                  !
                  write(out,"('molpro_duo: for class ',a,' is not implemented ')") field%class
                  stop 'molpro_duo: not for this object'
                  !
                end select
                !
              endif
              !
              f_t = matmul( conjg(transpose(a)),matmul(coupling,(b)) )
              !
              field%gridvalue(i) = real(f_t(1,1))
              !
              if (any( abs( aimag( f_t ) )>small_ ) ) then
                !
                write(out,"(/'molpro_duo: ',a,' ',2i3,'; duo-complex values ',8f8.1)") trim(field%class),field%iref,&
                                                                                       field%jref,f_t(:,:)
                                                                                       
                write(out,"('Please check the MOLPRO output and use FACTOR I if required')")
                write(out,"('It is important that MOLPROS LSX (spin-orbit-x), LX and DMX (dipole-x) values are used')")
                write(out,"('If this is input from older publications (before 2019) please try keyword LEGACY anywhere in input')")
                write(out,"('This will use the previously standard for the molpro objects')")
                stop 'molpro_duo error: duo-complex values'
                !
              endif
              !
              if (abs( real( f_t(1,1) ) )<=sqrt(small_) .and.abs(field%gridvalue(i))>sqrt(small_)) then
                !
                write(out,"('molpro_duo: ',a,' ',2i3,'; duo-zero values ',8f8.1)") trim(field%class),field%iref,field%jref,f_t(:,:)
                stop 'molpro_duo: duo-zero values ?'
                !
              endif
              !
            enddo
            !omp end parallel do

     end subroutine molpro_duo


     subroutine molpro_duo_old_2018(field)
        !
        use lapack,only : lapack_heev     
        !
        type(fieldT),intent(inout) :: field
        integer(ik) :: ix_lz_y,jx_lz_y,iroot,jroot,il_temp,ngrid
        complex(rk) :: a(2,2),b(2,2),coupling(2,2),f_t(2,2),b0(2,2),c
        real(rk)    :: lambda_i(2),lambda_j(2)
            !
            ngrid = grid%npoints
            !
            ix_lz_y = field%ix_lz_y
            jx_lz_y = field%jx_lz_y
            !
            a = 0 ; b = 0
            !
            a(1,1) =  cmplx(1.0_rk,0.0_rk, kind=rk)
            b(1,1) =  cmplx(1.0_rk,0.0_rk, kind=rk)
            a(2,2) =  cmplx(1.0_rk,0.0_rk, kind=rk)
            b(2,2) =  cmplx(1.0_rk,0.0_rk, kind=rk)
            !
            iroot = 1
            jroot = 1
            !
            if (ix_lz_y/=0) then
              a = 0 
              a(1,2) = cmplx(0.0_rk,ix_lz_y , kind=rk)
              a(2,1) = cmplx(0.0_rk,-ix_lz_y, kind=rk)
              !
              call lapack_heev(a,lambda_i)
              !
              ! swap to have the first root positive 
              !
              f_t = a
              a(:,1) = f_t(:,2)
              a(:,2) = f_t(:,1)
              !
              il_temp = lambda_i(2)
              lambda_i(2) = lambda_i(1)
              lambda_i(1) = il_temp
              !
              field%lambda = nint(lambda_i(1))
              !
              a = a*cmplx(0.0_rk,1.0_rk, kind=rk)
              !
            elseif (poten(field%istate)%parity%pm==-1) then 
              !
              a(1,1) =  cmplx(0.0_rk,1.0_rk, kind=rk)
              !
            endif
            !
            if (jx_lz_y/=0) then
              b = 0 
              b(1,2) = cmplx(0.0_rk,jx_lz_y , kind=rk)
              b(2,1) = cmplx(0.0_rk,-jx_lz_y, kind=rk)
              !
              b0 = b
              !
              call lapack_heev(b,lambda_j)
              !
              ! swap to have the first root positive 
              !
              f_t = b
              b(:,1) = f_t(:,2)
              b(:,2) = f_t(:,1)
              !
              il_temp = lambda_j(2)
              lambda_j(2) = lambda_j(1)
              lambda_j(1) = il_temp
              !
              field%lambdaj = nint(lambda_j(1))
              !
              b = b*cmplx(0.0_rk,1.0_rk, kind=rk)
              !
              f_t = matmul( conjg(transpose(b)),matmul(b0,(b)) )
              !
            elseif (poten(field%jstate)%parity%pm==-1) then 
              !
              b(1,1) =  cmplx(0.0_rk,1.0_rk, kind=rk)
              !
            endif
            !
            ! Check the selection rules
            select case(trim(field%class))
              !
            case('SPINORBIT','ABINITIO-SPINORBIT')
              !
              if ((nint(field%sigmaj-field%sigmai))/=(field%lambda-field%lambdaj)) then
                !
                ! try to select the root#2 for the i-state first 
                !
                if (field%lambda/=0.and.(nint(field%sigmaj-field%sigmai))==( lambda_i(2)-field%lambdaj ) ) then
                  !
                  il_temp = lambda_i(2)
                  lambda_i(2) = lambda_i(1)
                  lambda_i(1) = il_temp
                  !
                  field%lambda = nint(lambda_i(1))
                  !
                  f_t = a
                  a(:,1) = f_t(:,2)
                  a(:,2) = f_t(:,1)
                  !
                elseif( field%lambdaj/=0.and.(nint(field%sigmaj-field%sigmai))==( field%lambda-lambda_j(2) ) ) then
                  !
                  il_temp = lambda_j(2)
                  lambda_j(2) = lambda_j(1)
                  lambda_j(1) = il_temp
                  !
                  jroot = 2
                  field%lambdaj = nint(lambda_j(1))
                  !
                  f_t = b
                  b(:,1) = f_t(:,2)
                  b(:,2) = f_t(:,1)
                  !
                elseif( field%lambdaj/=0.and.(nint(field%sigmaj-field%sigmai))==( lambda_i(2)-lambda_j(2) ) ) then
                  !
                  il_temp = lambda_i(2)
                  lambda_i(2) = lambda_i(1)
                  lambda_i(1) = il_temp
                  !
                  field%lambda = nint(lambda_i(1))
                  !
                  il_temp = lambda_j(2)
                  lambda_j(2) = lambda_j(1)
                  lambda_j(1) = il_temp
                  !
                  jroot = 2
                  field%lambdaj = nint(lambda_j(1))
                  !
                  f_t = a
                  a(:,1) = f_t(:,2)
                  a(:,2) = f_t(:,1)
                  !
                  f_t = b
                  b(:,1) = f_t(:,2)
                  b(:,2) = f_t(:,1)
                  !
                else
                  !
                  write(out,"(/'molpro_duo: cannot find selecion rule for',2i8,3x,a)") field%iref,field%jref,trim(field%class)
                  write(out,"(' sigma = ',2f8.1,' lambda (i) = ',i4,' lamda (j) = ',i4)") field%sigmai, field%sigmaj, &
                                                                                                  int(lambda_i(1)),int(lambda_j(1))
                  stop 'molpro_duo: cannot find the selecion rules'
                endif
                !
              endif
              ! 
            case ('L+','ABINITIO-LX')
              !
              !write(out,"('molpro_duo: this L+-part is not implemented')")
              !stop 'molpro_duo: this L+-part is not implemented'
              ! 
            case ('DIPOLE')
              !
              !write(out,"('molpro_duo: this Dipole-part is not implemented')")
              !stop 'molpro_duo: this Dipole-part is not implemented'
              !
            case default
              !
              write(out,"(/'molpro_duo: this part is not implemented:',a)") trim(field%class)
              stop 'molpro_duo: this part is not implemented'
              !
            end select
            !
            !omp parallel do private(i) schedule(guided)
            do i=1,ngrid
              !
              coupling = 0 
              if (ix_lz_y==0.and.jx_lz_y==0) then
                !
                coupling(1,1) = field%gridvalue(i)*field%complex_f 
                !
              elseif(ix_lz_y/=0.and.jx_lz_y==0) then
                !
                !if ( lambda_i(1)/=1_rk.or.lambda_i(2)/=-1.0_rk) then
                !  !
                !  write(out,"('molpro_duo: lambda_i ',2f8.1,' are not -1 and 1, coupling')") lambda_i,field%iref,field%jref
                !  !stop 'molpro_duo: not for this object'
                !  !
                !endif
                !
                select case(trim(field%class))
                  !
                case('SPINORBIT','ABINITIO-SPINORBIT')
                  !
                  if (abs(nint(field%sigmaj-field%sigmai))/=abs(field%lambda-field%lambdaj)) then
                    !
                    write(out,"(/'molpro_duo: SO ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") &
                         field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                    stop 'molpro_duo: illegal selection rules '
                    !
                  endif
                  !
                  if ( field%sigmai<0 ) then
                    !
                    !write(out,"('molpro_duo: SO ',2i4,'; illegal reference sigmai <0 ?',2f8.1)") &
                    !    field%iref,field%jref,field%sigmai,field%sigmaj
                    !stop 'molpro_duo: illegal reference sigmai'
                    !
                  endif
                  !
                  ! for SOX it is always <1.3| which is given, i.e. we need to solve for the 1st, <1.2| component:
                  !
                  if (field%lambda>0) then 
                    !
                    coupling(2,1) = field%gridvalue(i)*field%complex_f
                    coupling(1,1) =-field%gridvalue(i)*field%complex_f*conjg(a(2,2))/conjg(a(1,2))
                    !
                  else 
                    !
                    coupling(2,1) = field%gridvalue(i)*field%complex_f
                    coupling(1,1) =-field%gridvalue(i)*field%complex_f*conjg(a(2,1))/conjg(a(1,1))
                    !
                  endif 
                  ! 
                case ('L+','ABINITIO-LX')
                  !
                  ! eigen-vector 2 is for Lambda
                  !
                  coupling(1,1) = field%gridvalue(i)*field%complex_f*cmplx(0.0_rk,-1.0_rk,kind=rk)  
                  coupling(2,1) = field%gridvalue(i)*field%complex_f*cmplx(0.0_rk, 1.0_rk,kind=rk)*conjg(a(1,2))/conjg(a(2,2))
                  ! 
                case ('DIPOLE')
                  !
                  ! eigen-vector 1 is for -Lambda
                  !
                  coupling(1,1) = field%gridvalue(i)*field%complex_f*(-sqrt(0.5_rk))
                  coupling(2,1) = field%gridvalue(i)*field%complex_f*(conjg(a(1,2)/a(2,2))*sqrt(0.5_rk))
                  !
                case default
                  !
                  write(out,"(/'molpro_duo (lambdaj=0): for class ',a,' is not implemented ')") field%class
                  stop 'molpro_duo: not for this object'
                  !
                end select
                !
              elseif(ix_lz_y==0.and.jx_lz_y/=0) then
                !
                !if (lambda_j(1)/=1_rk.or.lambda_j(2)/=-1.0_rk) then
                !  !
                !  !write(out,"('molpro_duo: lambda_j ',2f8.1,' are not -1 and 1, coupling for states',2i)") lambda_j,
                !                                                                                field%iref,field%jref
                !  !stop 'molpro_duo: not for this object'
                !  !
                !endif
                !
                select case(trim(field%class)) 
                  !
                case('SPINORBIT','ABINITIO-SPINORBIT')
                  !
                  if (abs(nint(field%sigmaj-field%sigmai))/=abs(field%lambda-field%lambdaj)) then
                    !
                    write(out,"(/'molpro_duo: SO ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") & 
                          field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                    stop 'molpro_duo: illegal selection rules '
                    !
                  endif
                  !
                  !if ( field%sigmaj<0 ) then
                  !  !
                  !  write(out,"('molpro_duo: SO ',2i4,'; illegal reference sigmaj <0 ',2f8.1)") & 
                  !        field%iref,field%jref,field%sigmai,field%sigmaj
                  !  stop 'molpro_duo: illegal reference sigmaj'
                  !  !
                  !endif
                  !
                  ! eigen-vector 1 is for Lambda
                  !
                  ! for SOX the non-zero is for |y> vector, i.e. the second component of coupling
                  !
                  !coupling(1,1) = -field%gridvalue(i)*field%complex_f*b(1,2)/b(2,2)
                  !coupling(1,2) = field%gridvalue(i)*field%complex_f

                  ! for SOX it is always <1.3| which is given, i.e. we need to solve for the 1st, <1.2| component:
                  !
                  if (field%lambdaj>0) then 
                    !
                    coupling(1,2) = field%gridvalue(i)*field%complex_f
                    coupling(1,1) =-field%gridvalue(i)*field%complex_f*b(2,1)/b(1,1)
                    !
                  else 
                    !
                    coupling(1,2) = field%gridvalue(i)*field%complex_f
                    coupling(1,1) =-field%gridvalue(i)*field%complex_f*b(2,2)/b(1,2)
                    !
                  endif 
                  ! 
                case('L+','ABINITIO-LX')
                  !
                  ! eigen-vector 1 is for Lambda
                  !
                  coupling(1,1) = field%gridvalue(i)*field%complex_f*cmplx(0.0_rk, 1.0_rk,kind=rk)  
                  coupling(1,2) = field%gridvalue(i)*field%complex_f*cmplx(0.0_rk,-1.0_rk,kind=rk)*b(1,2)/b(2,2)
                  ! 
                case ('DIPOLE')
                  !
                  ! eigen-vector 1 is for Lambda
                  !
                  coupling(1,1) = field%gridvalue(i)*field%complex_f*(-sqrt(0.5_rk))
                  coupling(1,2) = field%gridvalue(i)*field%complex_f*(b(1,2)/b(2,2)*sqrt(0.5_rk))
                  !
                case default
                  !
                  write(out,"(/'molpro_duo (lambdai=0): for class ',a,' is not implemented ')") field%class
                  stop 'molpro_duo: not for this object'
                  !
                end select
                !
              elseif (abs(ix_lz_y)/=abs(jx_lz_y)) then
                !
                select case(trim(field%class))
                  !
                case('SPINORBIT','ABINITIO-SPINORBIT')
                  !
                  if (abs(nint(field%sigmaj-field%sigmai))/=abs(field%lambda-field%lambdaj)) then
                    !
                    write(out,"(/'molpro_duo: SO ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") &
                          field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                    stop 'molpro_duo: illegal selection rules '
                    !
                  endif
                  !
                  if ( field%sigmai<0 ) then
                    !
                    !write(out,"('molpro_duo: SO ',2i4,'; illegal reference sigmai <0 ?',2f8.1)") &
                    !      field%iref,field%jref,field%sigmai,field%sigmaj
                    !stop 'molpro_duo: illegal reference sigmai'
                    !
                  endif
                  !
                  c = field%gridvalue(i)*field%complex_f
                  !
                  coupling(1,1) =  0
                  coupling(1,2) =  c
                  coupling(2,1) =  -c*conjg(a(1,1))*b(2,2)/(conjg(a(2,1))*b(1,2))
                  coupling(2,2) =  0
                  !
                case('DIPOLE')
                  !
                  if (abs(field%lambda-field%lambdaj)/=1) then
                    !
                    write(out,"('molpro_duo: DIPOLE ',2i4,'; illegal selection rules for lambda = ',2i4,' not +/-1')") & 
                          field%iref,field%jref,field%lambda,field%lambdaj
                    stop 'molpro_duo: illegal selection rules for transition dipole'
                    !
                  endif
                  !
                  c = -field%gridvalue(i)*field%complex_f*sqrt(0.5_rk)
                  !
                  ! maple for Lx
                  !c = -conjugate(A[1,2])*a/conjugate(A[2,2]), b = -a*B[1,2]/B[2,2], 
                  !d = B[1,2]*conjugate(A[1,2])*a/(B[2,2]*conjugate(A[2,2]))
                  ! m+ = -1/sqrt(2)l+
                  !
                  coupling(1,1) =  c
                  coupling(1,2) = -c*b(1,2)/b(2,2)
                  coupling(2,1) = -c*conjg(a(1,2))/conjg(a(2,2))
                  coupling(2,2) =  c*b(1,2)*conjg(a(1,2))/(conjg(a(2,2))*b(2,2))
                  !
                case('L+','ABINITIO-LX')
                  !
                  if (abs(field%lambda-field%lambdaj)/=1) then
                    !
                    write(out,"('molpro_duo: L+ ',2i4,'; illegal selection rules for lambda = ',2i4,' not +/-1')") & 
                          field%iref,field%jref,field%lambda,field%lambdaj
                    stop 'molpro_duo: L+ illegal selection rules for transition dipole'
                    !
                  endif
                  !
                  ! The <x|Lx+iLy|y> element is <x|Lx|y>
                  !
                  c = field%gridvalue(i)*field%complex_f
                  !
                  ! maple:
                  !
                  !d = -b*conjugate(A[1,2])/conjugate(A[2,2]), c = B[2,2]*b*conjugate(A[1,2])/(B[1,2]*conjugate(A[2,2])), 
                  !a = -B[2,2]*b/B[1,2]
                  !
                  coupling(1,2) =  c
                  coupling(1,1) = -c*b(2,2)/b(1,2)
                  coupling(2,1) =  c*b(2,2)*conjg(a(1,2))/(conjg(a(2,2))*b(1,2))
                  coupling(2,2) = -c*conjg(a(1,2))/conjg(a(2,2))
                  !
                case default
                  !
                  write(out,"('molpro_duo (lambdaj<>lambdai): for class ',a,' is not implemented ')") field%class
                  stop 'molpro_duo: not for this object'
                  !
                end select
                !
              else
                !
                select case(trim(field%class))
                  !
                case('SPINORBIT','ABINITIO-SPINORBIT')
                  !
                  if (nint(field%sigmaj-field%sigmai)/=0) then
                    !
                    write(out,"('molpro_duo: SOZ ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") & 
                          field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                    stop 'molpro_duo: SOZ illegal selection rules '
                    !
                  endif
                  !
                  if ( field%sigmai<0 ) then
                    !
                    !write(out,"('molpro_duo: SO ',2i4,'; illegal reference sigmai <0 ',2f8.1)") &
                    !      field%iref,field%jref,field%sigmai,field%sigmaj
                    !stop 'molpro_duo: illegal reference sigmai'
                    !
                  endif
                  !
                  c = field%gridvalue(i)*field%complex_f
                  !
                  !coupling(1,1) =  c
                  !coupling(1,2) =  0
                  !coupling(2,1) =  0
                  !coupling(2,2) = -c*b(1,2)/b(2,2)*conjg(a(1,1))/conjg(a(2,1))
                  !
                  coupling(1,1) =  0
                  coupling(1,2) =  c
                  coupling(2,1) = -conjg(a(1,1))*c*b(2,2)/(conjg(a(2,1))*b(1,2))
                  coupling(2,2) =  0
                  !
                case('DIPOLE')
                  !
                  if (abs(field%lambda-field%lambdaj)/=0) then
                    !
                    write(out,"('molpro_duo: DMZ ',2i4,'; illegal selection rules for lambda = ',2i4,' not 0')") &
                          field%iref,field%jref,field%lambda,field%lambdaj
                    stop 'molpro_duo: illegal selection rules for DMZ'
                    !
                  endif
                  !
                  c = field%gridvalue(i)*field%complex_f
                  !
                  !coupling(1,1) =  0
                  !coupling(1,2) =  c
                  !coupling(2,1) = -conjg(a(1,1))*c*b(2,2)/(conjg(a(2,1))*b(1,2))
                  !coupling(2,2) =  0
                  !
                  coupling(1,1) =  c
                  coupling(1,2) =  0
                  coupling(2,1) =  0
                  coupling(2,2) =  c
                  !
                case default
                  !
                  write(out,"('molpro_duo: for class ',a,' is not implemented ')") field%class
                  stop 'molpro_duo: not for this object'
                  !
                end select
                !
              endif
              !
              f_t = matmul( conjg(transpose(a)),matmul(coupling,(b)) )
              !
              field%gridvalue(i) = real(f_t(1,1))
              !
              if (any( abs( aimag( f_t ) )>small_ ) ) then
                !
                write(out,"(/'molpro_duo: ',a,' ',2i3,'; duo-complex values ',8f8.1)") trim(field%class),field%iref,&
                                                                                       field%jref,f_t(:,:)
                stop 'molpro_duo: duo-complex values ?'
                !
              endif
              !
              if (abs( real( f_t(1,1) ) )<=sqrt(small_) .and.abs(field%gridvalue(i))>sqrt(small_)) then
                !
                write(out,"('molpro_duo: ',a,' ',2i3,'; duo-zero values ',8f8.1)") trim(field%class),field%iref,field%jref,f_t(:,:)
                stop 'molpro_duo: duo-zero values ?'
                !
              endif
              !
            enddo
            !omp end parallel do

     end subroutine molpro_duo_old_2018


     !
     subroutine check_and_print_coupling(N,iverbose,fl,name)
       !
       type(fieldT),intent(in) :: fl(:)
       integer(ik),intent(in)  :: N,iverbose
       character(len=*),intent(in) :: name
       integer(ik)             :: i,istate,ngrid
       character(len=200) :: filename
       integer :: u1
         !
         if (N<=0.or.iverbose<5) return
         !
         write(out,'(/a)') trim(name)
         !
         ngrid = grid%npoints
         !
         ! double check
         !
         do i=1,N
           !
           istate = fl(i)%istate
           jstate = fl(i)%jstate
           !
           if (trim(fl(i)%class)=='SPINORBIT'.and.(abs(fl(i)%sigmai)>fl(i)%spini.or.abs(fl(i)%sigmaj)>fl(i)%spinj)) then
              write(out,'("For N =",i4," one of sigmas (",2f8.1,") large than  spins or undefined (",2f8.1,")")') & 
                        i,fl(i)%sigmai,fl(i)%sigmaj,fl(i)%spini,fl(i)%spinj
              stop 'illegal sigma or spin'
           endif
           if (nint(2.0_rk*fl(i)%spini)+1/=poten(istate)%multi.or. &
               nint(2.0_rk*fl(i)%spinj)+1/=poten(jstate)%multi ) then
              write(out,'("For N =",i3," multi (",2i3,") dont agree with either of multi (",2i3,") of states ",i2," and ",i2)') &
                        i,nint(2.0*fl(i)%spini)+1,nint(2.0*fl(i)%spinj)+1,poten(istate)%multi, &
                        poten(jstate)%multi,istate,jstate
              stop 'illegal multi in map_fields_onto_grid'
           endif
           !
           if (  fl(i)%lambda<bad_value+1.or.fl(i)%lambdaj<bad_value+1 ) then
              write(out,'("For N =",i3," lambdas  are undefined of states ",i2," and ",i2," ",a)') &
                        i,istate,jstate,trim(name)
              stop 'lambdas are undefined: map_fields_onto_grid'
           endif
           !
           if (  abs(fl(i)%lambda)/=poten(istate)%lambda.or.abs(fl(i)%lambdaj)/=poten(jstate)%lambda ) then
              write(out,'("For N =",i3," lambdas (",2i3,") dont agree with either of lambdas (",2i3,") ' // &
                                                                     'of states ",i2," and ",i2)') &
                        i,fl(i)%lambda,fl(i)%lambdaj,poten(istate)%lambda, &
                        poten(jstate)%lambda,istate,jstate
              stop 'illegal lambdas in map_fields_onto_grid'
           endif
           !
           if (  trim(name)=="Spin-Orbit:".and.(abs(fl(i)%sigmai)>fl(i)%spini.or.abs(fl(i)%sigmaj)>fl(i)%spinj) ) then
              write(out,'("For N =",i3," sigmas (",2f9.2,") dont agree with their spins (",2f9.2,") of states ",i2," and ",i2)') &
                        i,fl(i)%sigmai,fl(i)%sigmai,fl(i)%spini, &
                        fl(i)%spinj,istate,jstate
              stop 'illegal sigmas in map_fields_onto_grid'
           endif
           !
           if ( trim(name)=="Spin-Orbit:".and.(fl(i)%sigmai<bad_value+1.0.or.fl(i)%sigmaj<bad_value+1) ) then
              write(out,'("For N =",i3," sigmas are undefined for states ",i2," and ",i2," ",a)') &
                        i,istate,jstate,trim(name)
              stop 'sigmas are undefined: map_fields_onto_grid'
           endif
           !
           if ( trim(name)=="Spin-Orbit:".and.(fl(i)%sigmai<bad_value+1.0.or.fl(i)%sigmaj<bad_value+1.0) ) then
              write(out,'("For N =",i3," sigmas are undefined for states ",i2," and ",i2," ",a)') &
                        i,istate,jstate,trim(name)
              stop 'sigmas are undefined: map_fields_onto_grid'
           endif
           !
           if ( trim(name)=="<L+> functions:".and.istate==jstate ) then
              write(out,'("For N =",i3," Lx/L+ are defined for the same state ",i2," and ",i2," ",a)') &
                        i,istate,jstate,trim(name)
              stop 'illegal  - diagonal  - L+/Lx coupling: map_fields_onto_grid'
           endif
           !
           if ( trim(name(1:6))=="Lambda".and.istate/=jstate ) then
              write(out,'("For N =",i3," Lambda-doubling must be defined for the same state, not ",i2," and ",i2," ",a)') &
                        i,istate,jstate,trim(name)
              stop 'illegal  - non-diagonal  - Lambda doubling coupling: map_fields_onto_grid'
           endif
           !
           if ( trim(name)=="NACouplings".and.istate==jstate ) then
              write(out,'("For N =",i3," NACouplings must be defined for the different states, not ",i2," and ",i2," ",a)') &
                        i,istate,jstate,trim(name)
              stop 'illegal  - diagonal  - NAcoupling: map_fields_onto_grid'
           endif
           !
         enddo
         !
         do istate=1,N
           write(out,'(i4,2x,a)') istate,trim(fl(istate)%name)
         enddo
         !
         !
         ! write to file if required
         if(job%print_pecs_and_couplings_to_file .eqv. .true.) then
           ! set up filename
           write(filename, '(a)') trim(name)
           i=len_trim(filename)
           if( filename(i:i) == ':') filename(i:i)=' ' !remove trailing colons
           do i=1, len_trim(filename) ! remove/change some character in the file name
             if( filename(i:i) == ' ') filename(i:i)='_' !spaces to underscores
             if( filename(i:i) == '*') filename(i:i)='_' !stars to underscores
             if( filename(i:i) == '<') filename(i:i)='_' !< to underscores
             if( filename(i:i) == '>') filename(i:i)='_' !> to underscores
             if( filename(i:i) == '(') filename(i:i)='_' !( to underscores
             if( filename(i:i) == ')') filename(i:i)='_' !) to underscores
             if( filename(i:i) == '+') filename(i:i)='p' !+ to `p'
           enddo
           filename=trim(filename) // '.dat'

           call IOstart("check and_print coupling",u1)
           open(unit=u1, file=trim(filename), status='unknown',action='write')
           write(my_fmt,'(A,I0,A)') '(f18.10,', N, '(f18.9))'
           do i=1,ngrid
            write(u1 ,my_fmt) r(i),(fl(istate)%gridvalue(i),istate=1,N)
           enddo
           close(u1)
           call IOstop("check and_print coupling")
         endif
         !
         write(my_fmt,'(A,I0,A)') '("        r(Ang)  ",2x,', N, '(i9,10x))'
         write(out,my_fmt) (istate,istate=1,N)
         !
         write(my_fmt,'(A,I0,A)') '(f18.8,', N, '(1x,f18.8))'
         do i=1,ngrid
            write(out,my_fmt) r(i),(fl(istate)%gridvalue(i),istate=1,N)
         enddo
         !
     end subroutine check_and_print_coupling
     !

     !
     subroutine check_and_print_field(N,iverbose,fl,name)
       !
       type(fieldT),intent(in) :: fl(:)
       integer(ik),intent(in)  :: N,iverbose
       character(len=*),intent(in) :: name
       integer(ik)             :: i,istate,ngrid
       character(len=100)   :: filename ! file containing tabulated field (for plotting etc)
       integer, parameter  :: u1=1001
           !
           if (N<=0.or.iverbose<5) return
           !
           write(out,'(/a)') trim(name)
           !
           ngrid = grid%npoints
           !
           do istate=1,N
              !
              write(out,'(i6,2x,a30,2x,a7,a20)') istate,trim(fl(istate)%name), "type = " , trim(fl(istate)%type)
              !
              ! double check
              if (fl(istate)%lambda == 0 .and.fl(istate)%parity%pm==0 ) then
                 write(out,'("Please define the +/- symmetry of the Lambda=0 state (",i4,")")') istate
                 write(out,'("It is important for the SO component")')
                 stop 'Parity (+/-) for Lambda=0 undefined'
              endif
              !
           enddo
           ! 
           ! write to file if required
           if(job%print_pecs_and_couplings_to_file .eqv. .true.) then
             ! set up filename for output
             write(filename, '(a)') trim(name)
             i=len_trim(filename)
             if( filename(i:i) == ':') filename(i:i)=' ' !remove trailing colons
             do i=1, len_trim(filename) ! remove/change some character in the file name
               if( filename(i:i) == ' ') filename(i:i)='_' !spaces to underscores
               if( filename(i:i) == '*') filename(i:i)='_' !stars to underscores
               if( filename(i:i) == '<') filename(i:i)='_' !< to underscores
               if( filename(i:i) == '>') filename(i:i)='_' !> to underscores
               if( filename(i:i) == '(') filename(i:i)='_' !( to underscores
               if( filename(i:i) == ')') filename(i:i)='_' !) to underscores
               if( filename(i:i) == '+') filename(i:i)='p' !+ to `p'
             enddo
             filename=trim(filename) // '.dat'
             open(unit=u1, file=trim(filename), status='unknown',action='write')
             write(my_fmt, '(A,I0,A)') '(f18.10,', N, '(es22.14))'
             do i=1,ngrid
               write(u1, my_fmt) r(i),(fl(istate)%gridvalue(i),istate=1,N)
             enddo
             close(u1)
           endif
           !
           write(my_fmt, '(A,I0,A)') '("            r(Ang)",', N, '(i22))'
           write(out,my_fmt) (istate,istate=1,N)
           !
           write(my_fmt, '(A,I0,A)') '(f18.8,', N, '(es22.8))'
           do i=1,ngrid
              write(out,my_fmt) r(i),(fl(istate)%gridvalue(i),istate=1,N)
           enddo
           !
           ! write equilibrium properties
           write(out,*)
           write(out,'(A)') 'Equilibrium properties'
           !
           write(out, '(A18)',advance='no') 'imin = '
           do istate=1,N; write(out,'(I22)',advance='no') fl(istate)%imin ; enddo ! grid index for which poten is minimum
           write(out,*)
           !
           write(out, '(A18)',advance='no') 'r_imin / ang = '
           do istate=1,N; write(out,'(F22.10)',advance='no') fl(istate)%rimin ; enddo ! minimum of potential on the grid
           write(out,*)
           !
           write(out, '(A18)',advance='no') 'V(r_imin) / cm-1 = '
           do istate=1,N; write(out,'(F22.8)',advance='no') fl(istate)%Vimin ; enddo ! minimum of potential on the grid
           write(out,*)
           !
           write(out, '(A18)',advance='no') 'Has a single min?'
           do istate=1,N;
              if(fl(istate)%zHasMinimum) then
                   write(out,'(A22)',advance='no') 'Yes'
              else
                   write(out,'(A22)',advance='no') ' No'
              endif
           enddo ! minimum of potential on the grid
           write(out,*)
           !
           write(out, '(A18)',advance='no') 'True r_e  / ang ='
           do istate=1,N;
              if(fl(istate)%zHasMinimum) then
                write(out,'(F22.10)',advance='no') fl(istate)%re
              else
                write(out, '(a22)', advance='no') 'N.A.'
              endif
           enddo
           write(out,*)
           !
           write(out,'(/,A)') 'Derivatives at true re'
           !
           write(out, '(A18)',advance='no') 'der0, cm-1       ='
           do istate=1,N
              if(fl(istate)%zHasMinimum) then
                write(out,'(F22.10)',advance='no') fl(istate)%V0
              else
                write(out, '(a22)', advance='no') 'N.A.'
              endif
           enddo
           write(out,*)
           !
           write(out, '(A18)',advance='no') "der1, cm-1/ang   ="
           do istate=1,N
              if(fl(istate)%zHasMinimum) then
                write(out,'(F22.10)',advance='no') fl(istate)%V1
              else
                write(out, '(a22)', advance='no') 'N.A.'
              endif
           enddo
           write(out,*)
           !
           write(out, '(A18)',advance='no') "der2, cm-1/ang^2 ="
           do istate=1,N
              if(fl(istate)%zHasMinimum) then
                write(out,'(F22.10)',advance='no') fl(istate)%V2
              else
                write(out, '(a22)', advance='no') 'N.A.'
              endif
           enddo
           write(out,*)
           !
           write(out, '(A18)',advance='no') "der3, cm-1/ang^3 ="
           do istate=1,N
              if(fl(istate)%zHasMinimum) then
                write(out,'(F22.10)',advance='no') fl(istate)%V3
              else
                write(out, '(a22)', advance='no') 'N.A.'
              endif
           enddo
           write(out,*)
           !
           write(out, '(A18)',advance='no') "der4, cm-1/ang^4 ="
           do istate=1,N
              if(fl(istate)%zHasMinimum) then
                write(out,'(F22.10)',advance='no') fl(istate)%V4
              else
                write(out, '(a22)', advance='no') 'N.A.'
              endif
           enddo
           write(out,*)
           !
           write(out,*)
           !
           write(out, '(A18)',advance='no') "Harmonic we, cm-1="
           do istate=1,N
              if(fl(istate)%zHasMinimum) then
                write(out,'(F22.10)',advance='no') fl(istate)%we
              else
                write(out, '(a22)', advance='no') 'N.A.'
              endif
           enddo
           write(out,*)
           !
           write(out, '(A18)',advance='no') "Rotat. B0, cm-1="
           do istate=1,N
              if(fl(istate)%zHasMinimum) then
                write(out,'(F22.10)',advance='no') fl(istate)%B0
              else
                write(out, '(a22)', advance='no') 'N.A.'
              endif
           enddo
           write(out,*)
           !
           write(out, '(A18)',advance='no') "Anharm. const. xe="
           do istate=1,N
              if(fl(istate)%zHasMinimum) then
                write(out,'(F22.10)',advance='no') fl(istate)%xe
              else
                write(out, '(a22)', advance='no') 'N.A.'
              endif
           enddo
           write(out,*)
           !
           write(out, '(A18)',advance='no') "Coriol. ae, cm-1="
           do istate=1,N
              if(fl(istate)%zHasMinimum) then
                write(out,'(F22.10)',advance='no') fl(istate)%alphae
              else
                write(out, '(a22)', advance='no') 'N.A.'
              endif
           enddo
           write(out,*)
           !
           write(out, '(A18)',advance='no') "Centr. De, cm-1="
           do istate=1,N
              if(fl(istate)%zHasMinimum) then
                write(out,'(F22.10)',advance='no') fl(istate)%Debar
              else
                write(out, '(a22)', advance='no') 'N.A.'
              endif
           enddo
           write(out,*)
           !
           write(out, '(A18)',advance='no') "Y00, cm-1="
           do istate=1,N
              if(fl(istate)%zHasMinimum) then
                write(out,'(F22.10)',advance='no') fl(istate)%Y00
              else
                write(out, '(a22)', advance='no') 'N.A.'
              endif
           enddo
           write(out,*)
           !
           write(out,*)
           write(out,'(A)') 'Approximate J=0 vibrational energy levels (no couplings) '
           write(out,'(A)') ' given by E(v, J=0) = V(re) + Y00 + we*(v+0.5) - we*xe*(v+0.5)^2'
           write(out,'(A18)') 'v'
           do i=0, 3
           write(out,'(I18)',advance='no') i
           do istate=1,N
              if(fl(istate)%zHasMinimum) then
                write(out,'(F22.10)',advance='no') fl(istate)%Y00+fl(istate)%V1+fl(istate)%we*(real(i,rk)+0.5_rk) - &
                                                   fl(istate)%we*fl(istate)%xe*(real(i,rk)+0.5_rk)**2
              else
                write(out, '(a22)', advance='no') 'N.A.'
              endif
           enddo
           write(out,*)
           enddo
           write(out,*)
           write(out,'(a18)', advance='no') 'Spin ='
           do istate=1,N
                write(out,'(F22.1)',advance='no') fl(istate)%spini
           enddo
           write(out,*)
           write(out,'(a18)', advance='no') '|Lambda| ='
           do istate=1,N
                write(out,'(F22.1)',advance='no') real( fl(istate)%Lambda, rk)
           enddo
           write(out,*)
           write(out,'(a18)', advance='no') 'Physical J_min ='
           do istate=1,N
                write(out,'(F22.1)',advance='no') fl(istate)%Omega_min
           enddo
           write(out,*)
           write(out,*)
           write(out,'(A)') 'Approximate J=J_min vibrational-rotational energy levels (no couplings)'
           write(out,'(A)') ' given by E(v, J) = E(v, J=0) B0*J*(J+1)- ae*(v+0.5)*J*(J+1) - De*(J(J+1))^2'
           write(out,'(A18)') 'v'
           do i=0, 3
           write(out,'(I18)',advance='no') i
           do istate=1,N
              if(fl(istate)%zHasMinimum) then
                write(out,'(F22.10)',advance='no')  fl(istate)%Y00+fl(istate)%V1+fl(istate)%we*(real(i,rk)+0.5_rk) &
                                                   -fl(istate)%we*fl(istate)%xe*(real(i,rk)+0.5_rk)**2             &
                                                   +fl(istate)%B0*fl(istate)%Omega_min*(fl(istate)%Omega_min+1) &
                                                   -fl(istate)%alphae*(real(i,rk)+0.5_rk)*fl(istate)%Omega_min* &
                                                                                       (fl(istate)%Omega_min+1) &
                                                   -fl(istate)%Debar*(fl(istate)%Omega_min*(fl(istate)%Omega_min+1))**2
              else
                write(out, '(a22)', advance='no') 'N.A.'
              endif
           enddo
           write(out,*)
           enddo
           !
           write(out,*)
           !
     end subroutine check_and_print_field
     !
end subroutine map_fields_onto_grid



  subroutine duo_j0(iverbose_,J_list_,enerout,quantaout,nenerout)
    !
    use accuracy
    use timer
    use input
    use lapack,only : lapack_syev,lapack_heev,lapack_syevr     
     !
     implicit none
     !
     integer, parameter :: u1=102 ! unit for output file with J=0 energies
     integer, parameter :: u2=103 ! unit for output file with rovibronic energies
     !
     integer(ik),intent(in),optional  :: iverbose_
     !
     real(rk),intent(in),optional  :: J_list_(:) ! range of J values
     !
     real(rk),intent(out),optional  :: enerout(:,:,:)
     type(quantaT),intent(out),optional  :: quantaout(:,:,:)
     integer(ik),intent(out),optional  :: nenerout(:,:)
     !
     real(rk)                :: scale,sigma,omega,omegai,omegaj,spini,spinj,sigmaj,sigmai,jval
     integer(ik)             :: alloc,Ntotal,nmax,iterm,Nlambdasigmas,iverbose
     integer(ik)             :: ngrid,j,i,igrid,jgrid,kgrid
     integer(ik)             :: ilevel,mlevel,istate,imulti,jmulti,ilambda,jlambda,iso,jstate,jlevel,iobject
     integer(ik)             :: mterm,Nroots,tau_lambdai,irot,ilxly,itau,isigmav,isigmav_max
     integer(ik)             :: ilambda_,jlambda_,ilambda_we,jlambda_we,iL2,iss,isso,ibobrot,idiab,totalroots,ivib,jvib,v,inac
     integer(ik)             :: ipermute,istate_,jstate_,nener_total
     real(rk)                :: sigmai_,sigmaj_,f_l2,zpe,spini_,spinj_,omegai_,omegaj_,f_ss,f_bobrot,f_diabatic,f_nac
     real(rk)                :: sc, h12,f_rot,b_rot,epot,erot
     real(rk)                :: f_t,f_grid,energy_,f_s,f_l,psipsi_t,f_sr
     real(rk)                :: three_j_ref, three_j_,q_we, sigmai_we, sigmaj_we, SO,f_s1,f_s2,f_lo,f_o2,f_o1
     integer(ik)             :: isigma2
     character(len=1)        :: rng,jobz,plusminus(2)=(/'+','-'/)
     character(cl)           :: printout_
     real(rk)                :: vrange(2),veci(2,2),vecj(2,2),pmat(2,2),smat(2,2),maxcontr
     integer(ik)             :: irange(2),Nsym(2),jsym,isym,Nlevels,jtau,Nsym_,nJ,k
     integer(ik)             :: total_roots,irrep,jrrep,isr,ild
     real(rk),allocatable    :: eigenval(:),hmat(:,:),vec(:),vibmat(:,:),vibener(:),hsym(:,:),kinmat(:,:)
     real(rk),allocatable    :: LobAbs(:),LobWeights(:),LobDerivs(:,:),vibTmat(:,:)
     real(rk),allocatable    :: contrfunc(:,:),contrenergy(:),tau(:),J_list(:),Utransform(:,:,:)
     integer(ik),allocatable :: iswap(:),Nirr(:,:),ilevel2i(:,:),ilevel2isym(:,:),QNs(:)
     integer(ik),allocatable :: vib_count(:)
     type(quantaT),allocatable :: icontrvib(:),icontr(:)
     real(rk),allocatable    :: psi_vib(:),vec_t(:),vec0(:)
     integer(ik),allocatable :: ilambdasigmas_v_icontr(:,:)
     character(len=250),allocatable :: printout(:)
     double precision,parameter :: alpha = 1.0d0,beta=0.0d0
     type(matrixT)              :: transform(2)
     type(fieldT),pointer       :: field
     !
     real(ark),allocatable      :: psipsi_ark(:)
     real(rk),allocatable       :: mu_rr(:)
     !real(rk),allocatable      :: contrfunc_rk(:,:),vibmat_rk(:,:),matelem_rk(:,:),grid_rk(:)
     !real(rk)                  :: f_rk
     real(ark)                  :: f_ark
     character(len=cl)          :: filename,ioname
     integer(ik)                :: iunit,vibunit,imaxcontr,i0,imaxcontr_,mterm_,iroot,jroot,iomega_,jomega_,k_
     !
     real(rk)                   :: psi1,psi2,amplit1,amplit2,amplit3,diff,sum_wv,rhonorm,energy_unbound_sqrsqr
     integer(ik)                :: npoints_last,icount_max
     !
     ! Lambda-Sigma-> State-Omega contraction
     integer(ik) :: lambda_max,multi_max,lambda_min,iomega,Nomega_states
     integer(ik) :: ilambdasigma,Nlambdasigmas_max
     integer(ik) :: Nspins,Ndimen,jomega,v_i,v_j
     real(rk)    :: omega_min,omega_max,spin_min
     !
     type(contract_solT),allocatable :: contracted(:)
     real(rk),allocatable            :: vect_i(:),vect_j(:)
     !
     !
     ! open file for later (if option is set)
     if (job%print_rovibronic_energies_to_file ) &
         open(unit=u2, file='rovibronic_energies.dat',status='replace',action='write')
     !
     !
     ! define verbose level
     !
     iverbose = verbose
     if (present(iverbose_)) iverbose = iverbose_
     !
     ! define the range of the angular momentum
     !
     if (present(J_list_)) then
        !
        nJ = size(J_list_)
        !
        allocate(J_list(nJ),stat=alloc)
        !
        J_list = J_list_
        !
     else
        !
        nJ = size(job%J_list)
        allocate(J_list(nJ),stat=alloc)
        J_list = job%J_list
        !
     endif
     !
     if (iverbose>=4) call TimerStart('Map on grid')
     !
     ! Here we map all fields onto the same grid
     call map_fields_onto_grid(iverbose)
     !
     if (iverbose>=4) call TimerStop('Map on grid')
     !
     ! mapping grid
     !
     scale = amass/aston
     !
     h12 = 12.0_rk*hstep**2
     sc  = h12*scale

     !
     b_rot = aston/amass
     !
     ngrid = grid%npoints
     !
     ! First solve and contract the J=0,Sigma=0,Spin=0 problem and then use
     ! the corresponding eigenfunctions as basis set for the main hamiltonian.
     ! For this we use the vibrational hamiltonian + the L2(R) part.
     !
     if (iverbose>=3) write(out,'(/"Construct the J=0 matrix")')
     if (iverbose>=3) write(out,"(a)") 'Solving one-dimentional Schrodinger equations using : ' // trim(solution_method) 
     !
     allocate(psipsi_ark(ngrid),stat=alloc)
     call ArrayStart('psipsi_ark',alloc,size(psipsi_ark),kind(psipsi_ark))
     !
     !do i = 1,ngrid
     !  grid_rk(i) = 0.7_rk+(3.0_rk-.7_rk)/real(ngrid-1,rk)*real(i-1,rk)
     !enddo
     !
     if (iverbose>=4) call TimerStart('Solve vibrational part')
     !
     ! this will count all vibrational energies obtained for different istates
     totalroots = 0
     !
     zpe = 0
     !
     select case (job%contraction)  
       !
     case default
       !
       write(out,"('Error: contraction',a,' not implemented')") trim(job%contraction)
       stop 'Error: contraction is not implemented'
       !
     case ("OMEGA") 
       !
       if (action%intensity) then 
         write(out,'(/a)') "Error: INTENSITIES have not been implemented for the OMEGA contraction"
         stop "Error: INTENSITIES have not been implemented for the OMEGA contraction"
         !
       endif
       !
       if (iverbose>=5) then 
          write(out,'(/a)') "The Omega (adiabatic) representaion is under development."
          write(out,'(a)') "It is currently working only with PES, SO and Lx, other fields are being added."
          write(out,'(a)') "Please reported any bugs to Sergey Yurchenko"
       endif
       !
       multi_max = 1
       lambda_max = 0
       lambda_min = 10000
       !
       do istate = 1,nestates
         !
         multi_max = max(poten(istate)%multi,multi_max)
         lambda_max = max(poten(istate)%lambda,lambda_max)
         lambda_min = min(poten(istate)%lambda,lambda_min)
         !
       enddo
       !
       iomega = 0
       omega_min = 10000
       omega_max = -10000
       !
       do jlambda = lambda_min,lambda_max
         !
         do itau = 0,min(jlambda,1)
           !
           ilambda = (-1)**itau*jlambda
           !
           sigma = -real(multi_max-1,rk)*0.5_rk    !set sigma to -Smax (its most negative possible value)
           if ( sigma == -0.0_rk) sigma = +0.0_rk  !if sigma=0 use `positive signed' zero (for consistency across compilers)
           !
           do imulti = 1,multi_max
             !
             iomega = iomega + 1
             omega = sigma+ilambda
             !
             omega_min = min(omega_min,omega)
             omega_max = max(omega_max,omega)
             !
             sigma = sigma + 1.0_rk
             !
           enddo
           !
         enddo
         !
       enddo
       !
       if (abs(omega_min+omega_max)>small_) stop 'Something wrong with omega_min+omega_max'
       !
       Nomegas = nint(2.0_rk*omega_max)+1
       !
       ! count omega states 
       Nomega_states = 0
       Nlambdasigmas_max = 0
       omega = omega_min - 1.0_rk
       do iomega=1,Nomegas
          !
          omega = omega + 1.0_rk
          !
          call define_quanta_bookkeeping_Omega(0,omega,Nestates,Nlambdasigmas)
          Nomega_states = Nomega_states + Nlambdasigmas
          Nlambdasigmas_max = max(Nlambdasigmas,Nlambdasigmas_max)
          !
       enddo
       !
       Nspins = nint(real(multi_max-1,rk)*0.5_rk)
       !
       allocate(Omega_grid(Nomegas),stat=alloc)
       if (alloc/=0) stop 'Omega_grid cannot be allocated'
       !
       ! counter iLplus_Omega 
       allocate(iLplus_omega(Nomegas,Nomegas,Nlambdasigmas_max,Nlambdasigmas_max),stat=alloc)
       call ArrayStart('iLplus_omega',alloc,size(iLplus_omega),kind(iLplus_omega))
       iLplus_omega = 0 
       !
       ! counter iSplus_omega 
       allocate(iSplus_omega(Nomegas,Nomegas,Nlambdasigmas_max,Nlambdasigmas_max),stat=alloc)
       call ArrayStart('iSplus_omega',alloc,size(iSplus_omega),kind(iSplus_omega))
       iSplus_omega = 0 
       !
       ! counter iSR_omega 
       allocate(iSR_omega(Nomegas,Nomegas,Nlambdasigmas_max,Nlambdasigmas_max),stat=alloc)
       call ArrayStart('iSR_omega',alloc,size(iSR_omega),kind(iSR_omega))
       iSR_omega = 0
       !
       ! counter iBob_omega 
       allocate(iBob_omega(Nomegas,Nlambdasigmas_max,Nlambdasigmas_max),stat=alloc)
       call ArrayStart('iBob_omega',alloc,size(iBob_omega),kind(iBob_omega))
       iBob_omega = 0
       !
       ! counter ip2q_omega 
       allocate(ip2q_omega(Nomegas,Nomegas,Nlambdasigmas_max,Nlambdasigmas_max),stat=alloc)
       call ArrayStart('ip2q_omega',alloc,size(ip2q_omega),kind(ip2q_omega))
       ip2q_omega = 0
       !
       ! counter iQ_omega 
       allocate(iQ_omega(Nomegas,Nomegas,Nlambdasigmas_max,Nlambdasigmas_max),stat=alloc)
       call ArrayStart('iQ_omega',alloc,size(iQ_omega),kind(iQ_omega))
       iQ_omega = 0
       !
       ! counter iBRot_omega 
       allocate(iBRot_omega(Nomegas,Nlambdasigmas_max,Nlambdasigmas_max),stat=alloc)
       call ArrayStart('iBRot_omega',alloc,size(iBRot_omega),kind(iBRot_omega))
       iBob_omega = 0
       !
       ! counter iKin_omega 
       allocate(iKin_omega(Nomegas,Nlambdasigmas_max,Nlambdasigmas_max),stat=alloc)
       call ArrayStart('iKin_omega',alloc,size(iKin_omega),kind(iKin_omega))
       iKin_omega = 0
       !
       Nomega_states = 0
       omega = omega_min-1.0_rk 
       !
       if (iverbose>=5) write(out,'(/"     i   state   spin    sigma lambda   omega")')
       !
       do iomega=1,Nomegas
          !
          omega = omega + 1.0_rk
          !
          call define_quanta_bookkeeping_Omega(iverbose,omega,Nestates,Nlambdasigmas)
          !
          allocate(Omega_grid(iomega)%energy(Nlambdasigmas,ngrid),&
            Omega_grid(iomega)%vector(Nlambdasigmas,Nlambdasigmas,ngrid),stat=alloc)
          call ArrayStart('omegamat_grid',alloc,size(Omega_grid(iomega)%energy),kind(Omega_grid(iomega)%energy))
          call ArrayStart('omegamat_grid',alloc,size(Omega_grid(iomega)%vector),kind(Omega_grid(iomega)%vector))
          allocate(Omega_grid(iomega)%qn(Nlambdasigmas),stat=alloc)
          if (alloc/=0) stop 'Omega_grid(iomega)%qn cannot be allocated'
          allocate(Omega_grid(iomega)%basis(Nlambdasigmas),stat=alloc)
          if (alloc/=0) stop 'Omega_grid(iomega)%basis cannot be allocated'
          !
          Omega_grid(iomega)%Nstates = Nlambdasigmas
          Omega_grid(iomega)%omega = omega
          !
          do i = 1,Nlambdasigmas
            !
            Nomega_states = Nomega_states + 1
            Omega_grid(iomega)%basis(i)%istate = quanta(i)%istate
            Omega_grid(iomega)%basis(i)%name = trim(quanta(i)%name)
            Omega_grid(iomega)%basis(i)%sigma = quanta(i)%sigma
            Omega_grid(iomega)%basis(i)%ilambda = quanta(i)%ilambda
            Omega_grid(iomega)%basis(i)%spin = quanta(i)%spin
            Omega_grid(iomega)%basis(i)%omega = omega
            Omega_grid(iomega)%basis(i)%ilevel = ilambdasigma
            !
          enddo
          !
       enddo
       !
       ! Grid independent Splus-spin matrix
       !
       Nspins = nint(real(multi_max-1,rk)*0.5_rk)
       !
       spin_min = 0
       if (mod(nint(2.0_rk*multi_max+1.0_rk),2)==1) spin_min = 0.5
       !
       ! count and creat Lplus_omega object from lxly
       !
       ! Lplus/Lminus
       !
       ! count coupling 
       call L_omega_create(NLplus_omega,onlycount=.true.)
       !
       if (NLplus_omega/=0.and..not.fields_allocated) then
         !
         allocate(L_omega_obj(NLplus_omega),stat=alloc)
         if (alloc/=0) stop 'L_omega_obj cannot be allocated'
         !
         do i = 1,NLplus_omega
           L_omega_obj(i)%type = "grid"
           L_omega_obj(i)%name = "L+ Omega obj"
           allocate(L_omega_obj(i)%gridvalue(ngrid),stat=alloc)
           call ArrayStart("L+ Omega obj",alloc,ngrid,kind(L_omega_obj(i)%gridvalue))
           L_omega_obj(i)%gridvalue = 0
         enddo
         !
         call L_omega_create(NLplus_omega,onlycount=.false.)
         !
       endif
       !
       ! count the number of different Splus objects
       !
       call S_omega_create(NSplus_omega,onlycount=.true.)
       !
       if (NSplus_omega/=0.and..not.fields_allocated) then 
         !
         allocate(S_omega_obj(NSplus_omega),stat=alloc)
         if (alloc/=0) stop 'S_omega_obj cannot be allocated'
         !
         do i = 1,NSplus_omega
           S_omega_obj(i)%type = "grid"
           S_omega_obj(i)%name = "S+ Omega obj"
           allocate(S_omega_obj(i)%gridvalue(ngrid),stat=alloc)
           call ArrayStart("S+ Omega obj",alloc,ngrid,kind(S_omega_obj(i)%gridvalue))
           S_omega_obj(i)%gridvalue = 0
         enddo
         !
         call S_omega_create(NSplus_omega,onlycount=.false.)
         !
       endif
       !
       ! count the number of different SR objects
       !
       call SR_omega_create(NSR_omega,onlycount=.true.)
       !
       if (NSR_omega/=0.and..not.fields_allocated) then 
         !
         allocate(SR_omega_obj(NSR_omega),stat=alloc)
         if (alloc/=0) stop 'S_omega_obj cannot be allocated'
         !
         do i = 1,NSR_omega
           SR_omega_obj(i)%type = "grid"
           SR_omega_obj(i)%name = "SR Omega obj"
           allocate(SR_omega_obj(i)%gridvalue(ngrid),stat=alloc)
           call ArrayStart("SR Omega obj",alloc,ngrid,kind(SR_omega_obj(i)%gridvalue))
           SR_omega_obj(i)%gridvalue = 0
         enddo
         !
         call SR_omega_create(NSR_omega,onlycount=.false.)
         !
       endif
       !
       ! count the number of different BOB objects
       !
       call BOB_omega_create(NBob_omega,onlycount=.true.)
       !
       if (NBob_omega/=0.and..not.fields_allocated) then 
         !
         allocate(BOB_omega_obj(NBOB_omega),stat=alloc)
         if (alloc/=0) stop 'bob_omega_obj cannot be allocated'
         !
         do i = 1,NBob_omega
           Bob_omega_obj(i)%type = "grid"
           Bob_omega_obj(i)%name = "BOB Omega obj"
           allocate(bob_omega_obj(i)%gridvalue(ngrid),stat=alloc)
           call ArrayStart("BOB Omega obj",alloc,ngrid,kind(bob_omega_obj(i)%gridvalue))
           bob_omega_obj(i)%gridvalue = 0
         enddo
         !
         call Bob_omega_create(NBob_omega,onlycount=.false.)
         !
       endif
       !
       ! count the number of different p2q objects
       !
       call P2Q_omega_create(Np2q_omega,onlycount=.true.)
       !
       if (Np2q_omega/=0.and..not.fields_allocated) then 
         !
         allocate(p2q_omega_obj(Np2q_omega),stat=alloc)
         if (alloc/=0) stop 'p2q_omega_obj cannot be allocated'
         !
         do i = 1,Np2q_omega
           p2q_omega_obj(i)%type = "grid"
           p2q_omega_obj(i)%name = "P2Q Omega obj"
           allocate(p2q_omega_obj(i)%gridvalue(ngrid),stat=alloc)
           call ArrayStart("P2Q Omega obj",alloc,ngrid,kind(p2q_omega_obj(i)%gridvalue))
           p2q_omega_obj(i)%gridvalue = 0
         enddo
         !
         call P2Q_omega_create(Np2q_omega,onlycount=.false.)
         !
       endif
       !
       ! count the number of different p2q objects
       !
       call Q_omega_create(Nq_omega,onlycount=.true.)
       !
       if (Nq_omega/=0.and..not.fields_allocated) then 
         !
         allocate(q_omega_obj(Nq_omega),stat=alloc)
         if (alloc/=0) stop 'q_omega_obj cannot be allocated'
         !
         do i = 1,Nq_omega
           q_omega_obj(i)%type = "grid"
           q_omega_obj(i)%name = "Q Omega obj"
           allocate(q_omega_obj(i)%gridvalue(ngrid),stat=alloc)
           call ArrayStart("Q Omega obj",alloc,ngrid,kind(q_omega_obj(i)%gridvalue))
           q_omega_obj(i)%gridvalue = 0
         enddo
         !
         call Q_omega_create(Nq_omega,onlycount=.false.)
         !
       endif
       ! count the number of different BRot objects
       !
       call Brot_omega_create(NBRot_omega,onlycount=.true.)
       !
       if (NBRot_omega/=0.and..not.fields_allocated) then 
         !
         allocate(BRot_omega_obj(NBRot_omega),stat=alloc)
         if (alloc/=0) stop 'BRot_omega_obj cannot be allocated'
         !
         do i = 1,NBRot_omega
           BRot_omega_obj(i)%type = "grid"
           BRot_omega_obj(i)%name = "BRot Omega obj"
           allocate(BRot_omega_obj(i)%gridvalue(ngrid),stat=alloc)
           call ArrayStart("BRot Omega obj",alloc,ngrid,kind(BRot_omega_obj(i)%gridvalue))
           BRot_omega_obj(i)%gridvalue = 0 
         enddo
         !
         call Brot_omega_create(NBRot_omega,onlycount=.false.)
         !
       endif
       !
       ! count the number of different Kin objects
       !
       call Kin_omega_create(NKin_omega,onlycount=.true.)
       !
       if (NKin_omega/=0.and..not.fields_allocated) then 
         !
         allocate(Kin_omega_obj(NKin_omega),stat=alloc)
         if (alloc/=0) stop 'Kin_omega_obj cannot be allocated'
         !
         do i = 1,NKin_omega
           Kin_omega_obj(i)%type = "grid"
           Kin_omega_obj(i)%name = "Kin Omega obj"
           allocate(Kin_omega_obj(i)%gridvalue(ngrid),stat=alloc)
           call ArrayStart("Kin Omega obj",alloc,ngrid,kind(Kin_omega_obj(i)%gridvalue))
           !
           allocate(Kin_omega_obj(i)%matelem(ngrid,ngrid),stat=alloc)
           call ArrayStart("Kin Omega obj",alloc,ngrid*ngrid,kind(Kin_omega_obj(i)%matelem))
           !
           Kin_omega_obj(i)%matelem = 0
           !
         enddo
         !
         call Kin_omega_create(NKin_omega,onlycount=.false.)
         !
       endif
       ! 
       ! Diagonalise the PECs+SOCs+couplings and transform all other curves to the Omega representation 
       !
       call Transfrorm_Sigma_Lambda_to_Omega_representation(iverbose,sc,Nlambdasigmas_max,Nomega_states)
       !
       deallocate(iLPlus_omega,iSPlus_omega,iSR_omega,iBob_omega,ip2q_omega,iQ_omega,iKin_omega,iBRot_omega)
       call ArrayStop('iLplus_omega')
       call ArrayStop('iSplus_omega')
       call ArrayStop('iSR_omega')
       call ArrayStop('iBob_omega')
       call ArrayStop('ip2q_omega')
       call ArrayStop('iQ_omega')
       call ArrayStop('iBRot_omega')
       call ArrayStop('iKin_omega')
       !
       ! print out some fields in the new Omega reprsentation 
       !
       call Print_fileds_in_Omega_representation(iverbose,Nomega_states)
       !
       ! Eigensolve the vibrational problem for individual omega states
       !
       allocate(contracted(Nomegas))
       !
       do iomega=1,Nomegas
          !
          Nlambdasigmas = omega_grid(iomega)%Nstates
          Ndimen = Nlambdasigmas*ngrid
          contracted(iomega)%Ndimen = Ndimen
          !
          allocate(contracted(iomega)%vector(Ndimen,Ndimen),contracted(iomega)%energy(Ndimen),&
                   contracted(iomega)%ilevel(Ndimen),stat=alloc)
          call ArrayStart('contracted%vector',alloc,size(contracted(iomega)%vector),kind(contracted(iomega)%vector))
          call ArrayStart('contracted%energy',alloc,size(contracted(iomega)%energy),kind(contracted(iomega)%energy))
          call ArrayStart('contracted%ilevel',alloc,size(contracted(iomega)%ilevel),kind(contracted(iomega)%ilevel))
          !
       enddo       
       !
       allocate(contrenergy(ngrid*Nomega_states),stat=alloc)
       call ArrayStart('contrenergy',alloc,size(contrenergy),kind(contrenergy))
       !
       !allocate(contrfunc(ngrid,ngrid*Nomega_states),stat=alloc)
       !call ArrayStart('contrfunc',alloc,size(contrfunc),kind(contrfunc))
       !
       allocate(icontrvib(ngrid*Nomega_states),stat=alloc)
       !
       call Solve_vibrational_problem_for_Omega_states(iverbose,ngrid,Nomega_states,sc,totalroots,icontrvib,contrenergy,contracted)
       !
       ! Now we need to compute all vibrational matrix elements of all field of the Hamiltonian, except for the potentials V,
       ! which together with the vibrational kinetic energy operator are diagonal on the contracted basis developed
       !
       ! allocate arrays for matrix elements for all hamiltonian fields
       !
       ! introducing a new field for the centrifugal matrix
       !
       if (.not.fields_allocated) then
         allocate(brot(1),stat=alloc)
       endif
       !
       allocate(vect_i(ngrid),vect_j(ngrid),stat=alloc)
       call ArrayStart('vect_ij',alloc,size(vect_i),kind(vect_i))
       call ArrayStart('vect_ij',alloc,size(vect_j),kind(vect_j))
       !
       do iobject = 1,7
          !
          select case (iobject)
            !
          case (1)
            Nmax = NLplus_omega
          case (2)
            Nmax = NSplus_omega
          case (3)
            Nmax = NSR_omega
          case (4)
            Nmax = NBob_omega
          case (5)
            Nmax = Np2q_omega
          case (6)
            Nmax = Nq_omega
          case (7)
            !
            Nmax = NBrot_omega
            !
            !field => brot(1)
            !field%name = 'BROT'
            !if (.not.fields_allocated) then 
            !  allocate(field%gridvalue(ngrid),stat=alloc)
            !  call ArrayStart(field%name,alloc,size(field%gridvalue),kind(field%gridvalue))
            !endif 
            !field%gridvalue(:) = b_rot/r(:)**2*sc
            !Nmax = 1
          case default 
            Nmax = 0
          end select
          !
          do iterm = 1,Nmax
             !
             isigmav_max = 1 ! permuations (account for -omega where omega is the component from the input, only needed for L and S)
             !
             select case (iobject)
               !
             case (1)
               field => L_omega_obj(iterm)
               field%gridvalue(:) = field%gridvalue(:)*b_rot/r(:)**2*sc
               isigmav_max = 1
             case (2)
               field => S_omega_obj(iterm)
               field%gridvalue(:) = field%gridvalue(:)*b_rot/r(:)**2*sc
               isigmav_max = 1
             case (3)
               field => SR_omega_obj(iterm)
               field%gridvalue(:) = field%gridvalue(:)*sc
             case (4)
               field => bob_omega_obj(iterm)
               field%gridvalue(:) = field%gridvalue(:)*b_rot/r(:)**2*sc
             case (5)
               field => p2q_omega_obj(iterm)
               field%gridvalue(:) = field%gridvalue(:)*sc
             case (6)
               field => q_omega_obj(iterm)
               field%gridvalue(:) = field%gridvalue(:)*sc
             case (7)
               !field => brot(1)
               field => brot_omega_obj(iterm)
               field%gridvalue(:) = field%gridvalue(:)*b_rot/r(:)**2*sc
             case default 
               stop 'illegal object in LambdaSigma-Omega'
             end select
             !
             ! check if the field wass not allocated at previous call to prevent multiple allocatetions 
             !
             if (.not.fields_allocated) then 
               allocate(field%matelem(totalroots,totalroots),stat=alloc)
               call ArrayStart(field%name,alloc,size(field%matelem),kind(field%matelem))
             endif
             !
             if (associated(field%matelem) .eqv. .false.) then
               allocate(field%matelem(totalroots,totalroots),stat=alloc)
               call ArrayStart(field%name,alloc,size(field%matelem),kind(field%matelem))
             endif
             !
             field%matelem = 0 
             !
             iomega  = field%iomega
             jomega  = field%jomega
             !
             ilevel  = field%ilevel
             jlevel  = field%jlevel
             !
             !Nlambdasigmas = omega_grid(iomega)%Nstates
             !Nlambdasigmas_= omega_grid(jomega)%Nstates
             !
             !omp parallel do private(ilevel,jlevel) schedule(guided)
             !
             do iroot = 1,totalroots
               !
               i = icontrvib(iroot)%v+1
               iomega_ = icontrvib(iroot)%iomega
               omegai_ = icontrvib(iroot)%omega
               !
               !if (iomega/=iomega_.and.jomega/=iomega_) cycle 
               !
               do jroot = 1,iroot
                 !
                 j = icontrvib(jroot)%v+1
                 jomega_ = icontrvib(jroot)%iomega
                 omegaj_ = icontrvib(jroot)%omega
                 !
                 do isigmav = 0,isigmav_max
                   !
                   ! the permutation is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
                   ! avoid the double counting.
                   if( isigmav==1.and. nint(abs( field%omegai ) + abs( field%omegaj ))==0 ) cycle
                   !           
                   ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                   omegai = field%omegai*(-1)**isigmav
                   omegaj = field%omegaj*(-1)**isigmav
                   !
                   ! proceed only if the quantum numbers of the field equal to the corresponding <i| and |j> quantum numbers:
                   !
                   if( nint(omegai_-omegai)/=0.or.nint(omegaj_-omegaj)/=0 ) cycle
                   !
                   !if (iomega/=jomega_.and.jomega/=jomega_) cycle
                   !
                   !ilevel_ = ilevel
                   !jlevel_ = jlevel
                   !
                   !if (isigmav==1) then
                   !  ilevel_ = jlevel
                   !  jlevel_ = ilevel
                   !endif
                   !
                   !if (ilevel_/=ilevel.or.jlevel_/=ilevel.or.nint(omegai_-omega)/=0) cycle
                   !
                   vect_i(1:Ngrid) = contracted(iomega_)%vector((ilevel-1)*ngrid+1:ilevel*ngrid,i)
                   vect_j(1:Ngrid) = contracted(jomega_)%vector((jlevel-1)*ngrid+1:jlevel*ngrid,j)
                   !
                   ! in the grid representation of the vibrational basis set
                   ! the matrix elements are evaluated simply by a summation of over the grid points
                   !
                   field%matelem(iroot,jroot)  = sum(vect_i(:)*(field%gridvalue(:))*vect_j(:))
                   !
                   ! If intensity%threshold%dipole is given and TM is smaller than this threshold set the TM-value to zero
                   ! is applied to the dipole (iobject=Nobjects) and quadrupole (iobject=Nobjects-3) moments 
                   !if (iobject==Nobjects-3.or.iobject==Nobjects) then
                   !  if (abs(field%matelem(ilevel,jlevel))<intensity%threshold%dipole) field%matelem(ilevel,jlevel) = 0 
                   !endif
                   !
                   field%matelem(jroot,iroot) = field%matelem(iroot,jroot)
                   !
                   !
                 enddo
                 !
               enddo
             enddo
             !
          enddo
          !
       enddo
       !
       deallocate(vect_i,vect_j)
       call ArrayStop('vect_ij')
       !
       fields_allocated = .true.
       !
       if (allocated(psipsi_ark)) then 
         deallocate(psipsi_ark)
         call ArrayStop('psipsi_ark')
       endif
       !
       if (present(nenerout)) nenerout = 0
       !
       if (present(enerout)) then
         enerout = 0
       endif
       !
     case ("VIB") 
       !
       allocate(icontrvib(ngrid*Nestates),stat=alloc)
       allocate(contrfunc(ngrid,ngrid*Nestates),stat=alloc)
       call ArrayStart('contrfunc',alloc,size(contrfunc),kind(contrfunc))
       call ArrayStart('icontrvib',alloc,ngrid*Nestates*(14+6*2),ik)
       !
       allocate(vibmat(ngrid,ngrid),vibener(ngrid),contrenergy(ngrid*Nestates),vec(ngrid),stat=alloc)
       call ArrayStart('vibmat',alloc,size(vibmat),kind(vibmat))
       call ArrayStart('vibener',alloc,size(vibener),kind(vibmat))
       call ArrayStart('contrenergy',alloc,size(contrenergy),kind(contrenergy))
       call ArrayStart('vec',alloc,size(vec),kind(vec))
       allocate(kinmat(ngrid,ngrid),stat=alloc)
       call ArrayStart('kinmat',alloc,size(kinmat),kind(kinmat))
       !
       if (trim(solution_method)=="LOBATTO") then 
         allocate(LobAbs(ngrid),LobWeights(ngrid),LobDerivs(ngrid,ngrid),vibTmat(ngrid,ngrid),stat=alloc)
         call ArrayStart('LobAbs',alloc,size(LobAbs),kind(LobAbs))
         call ArrayStart('LobWeights',alloc,size(LobWeights),kind(LobWeights))
         call ArrayStart('LobDerivs',alloc,size(LobDerivs),kind(LobDerivs))
         call ArrayStart('vibTmat',alloc,size(vibTmat),kind(vibTmat))
       endif
       !
       ! Lobbato grids, abcissas and weighs 
       if (grid%nsub == 6) then
         call LobattoAbsWeights(LobAbs,LobWeights,ngrid,grid%rmin,grid%rmax)
         call derLobattoMat(LobDerivs,ngrid-2,LobAbs,LobWeights)
       endif
       !
       !
       call kinetic_energy_grid_points(ngrid,kinmat,vibTmat,LobWeights,LobDerivs)
       !
       do istate = 1,Nestates
         !
         !vibmat = 0
         !
         if (iverbose>=6) write(out,'("istate = ",i0)') istate
         !
         ! reconstruct quanta for the bra-state
         !
         imulti = poten(istate)%multi
         ilambda = poten(istate)%lambda
         spini = poten(istate)%spini
         !
         if (iverbose>=4) call TimerStart('Build vibrational Hamiltonian')
         !
         vibmat = kinmat
         !
         !$omp parallel do private(igrid,f_rot,epot,f_l2,iL2,erot) shared(vibmat) schedule(guided)
         do igrid =1, ngrid
           !
           if (iverbose>=6) write(out,'("igrid = ",i0)') igrid
           !
           ! the centrifugal factor will be needed for the L**2 term
           !
           f_rot=b_rot/r(igrid)**2*sc
           !
           !
           ! the diagonal term with the potential function
           !
           epot = poten(istate)%gridvalue(igrid)*sc
           !
           !
           ! Another diagonal term:
           ! The L^2 term (diagonal): (1) L2(R) is used if provided otherwise
           ! an approximate value Lambda*(Lamda+1) is assumed.
           !
           f_l2 = 0 ! real(ilambda*(ilambda+1),rk)*f_rot
           do iL2 = 1,Nl2
             if (L2(iL2)%istate==istate.and.L2(iL2)%jstate==istate) then
               f_l2 = f_rot*L2(iL2)%gridvalue(igrid)
               exit
             endif
           enddo
           !
           erot = f_l2
           !
           ! the diagonal matrix element will include PEC +L**2 as well as the vibrational kinetic contributions.
           vibmat(igrid,igrid) = vibmat(igrid,igrid) + epot + erot
           !
         enddo
         !$omp end parallel do
         !
         if (iverbose>=4) call TimerStop('Build vibrational Hamiltonian')
         !
         select case (trim(poten(istate)%integration_method))
           !
         case ('NUMEROV') 
           !
           ! we need only these many roots
           !
           nroots = min(job%vibmax(istate),Ngrid)
           allocate(mu_rr(1:ngrid),stat=alloc)
           call ArrayStart('mu_rr',alloc,size(mu_rr),kind(mu_rr))
           !
           mu_rr = 2.0_rk*b_rot
           !
           !write(out,"('Error: ME_numerov is not implemented yet')")
           !stop 'Error: ME_numerov is not implemented yet'
           !
           !iverbose
           !
           nroots = min(nroots,ngrid-2)
           !
           call ME_numerov(nroots,(/grid%rmin,grid%rmax/),ngrid-1,ngrid-1,r,poten(istate)%gridvalue,mu_rr,1,0,&
                           job%vibenermax(istate),iverbose,vibener(1:nroots+1),vibmat)
           deallocate(mu_rr)
           call ArrayStop('mu_rr')
           !
           vibener = vibener*sc
           !
           ! or as many as below job%upper_ener if required by the input
           if ((job%vibenermax(istate))*sc<safe_max) then
             nroots = maxloc(vibener(:)-vibener(1),dim=1,mask=vibener(:).le.job%vibenermax(istate)*sc)
           endif
           !
           if (istate==1) zpe = vibener(1)
           !
         case ('NONE','RAW')
           !
           vibener = poten(istate)%gridvalue*sc
           !
           nroots = Ngrid
           !
           vibmat = 0
           do i=1,nroots
             vibmat(i,i) = 1.0_ark !/sqrt(hstep)
           enddo
           !
           if (istate==1) zpe = minval(vibener(:))
           !
         case default
           !
           if (job%vibmax(istate)>ngrid/2) then
              !
              call lapack_syev(vibmat,vibener)
              !
              ! we need only these many roots
              Nroots = min(ngrid,job%vibmax(istate))
              !
              ! or as many as below job%upper_ener if required by the input
              if ((job%vibenermax(istate))*sc<safe_max) then
                nroots = maxloc(vibener(:)-vibener(1),dim=1,mask=vibener(:).le.job%vibenermax(istate)*sc)
              endif
              !
            else
              !
              ! some diagonalizers needs the following parameters to be defined
              !
              ! diagonalize the vibrational hamiltonian using the DSYEVR routine from LAPACK
              ! DSYEVR computes selected eigenvalues and, optionally, eigenvectors of a real n by n symmetric matrix A.
              ! The matrix is first reduced to tridiagonal form, using orthogonal similarity transformations.
              ! Then whenever possible, DSYEVR computes the eigenspectrum using Multiple Relatively Robust Representations (MR).
              !
              jobz = 'V'
              vrange(1) = -0.0_rk ; vrange(2) = (job%vibenermax(istate))*sc
              if (.not.job%zShiftPECsToZero) vrange(1) = -safe_max
              irange(1) = 1 ; irange(2) = min(job%vibmax(istate),Ngrid)
              nroots = Ngrid
              rng = 'A'
              !
              if (job%vibmax(istate)/=1e8) then
                 rng = 'I'
              elseif (job%vibenermax(istate)<1e8) then
                 rng = 'V'
              endif
              !
              call lapack_syevr(vibmat,vibener,rng=rng,jobz=jobz,iroots=nroots,vrange=real(vrange,kind=8),irange=irange)
              !
              !call lapack_syev(vibmat,vibener)
              !
           endif
           !
           ! ZPE is obatined only from the lowest state
           !
           if (istate==1) zpe = vibener(1)
           !
         end select
         !
         if (nroots<1) then
           nroots = 1
           vibener = 0
           vibmat = 0
           vibmat(1,1) = 1.0_rk
         endif
         !
         ! write the pure vibrational energies and the corresponding eigenfunctions into global matrices
         contrfunc(:,totalroots+1:totalroots+nroots) = vibmat(:,1:nroots)
         contrenergy(totalroots+1:totalroots+nroots) = vibener(1:nroots)
         !
         !vibmat_rk = vibmat
         !
         !call schmidt_orthogonalization(ngrid,nroots,vibmat_rk)
         !
         !contrfunc_rk(:,totalroots+1:totalroots+nroots) = vibmat_rk(:,1:nroots)
         !
         ! assign the eigenstates with quanta
         do i=1,nroots
           icontrvib(totalroots + i)%istate =  istate
           icontrvib(totalroots + i)%v = i-1
         enddo
         !
         ! increment the global counter of the vibrational states
         !
         totalroots = totalroots + nroots
         !
       enddo
       !
       ! sorting basis states (energies, basis functions and quantum numbers) from different
       ! states all together according with their energies
       !
       do ilevel = 1,totalroots
         !
         energy_ = contrenergy(ilevel)
         !
         ! skip sorting for DVR-raw representaions 
         if (icontrvib(ilevel)%istate>Nrefstates) cycle
         !
         do jlevel=ilevel+1,totalroots
           !
           if (icontrvib(jlevel)%istate>Nrefstates) cycle
           !
           if ( energy_>contrenergy(jlevel) ) then
             !
             ! energy
             !
             energy_=contrenergy(jlevel)
             contrenergy(jlevel) = contrenergy(ilevel)
             contrenergy(ilevel) = energy_
             !
             ! basis function
             !
             vec(:) = contrfunc(:,jlevel)
             contrfunc(:,jlevel) = contrfunc(:,ilevel)
             contrfunc(:,ilevel) = vec(:)
             !
             ! qunatum numbers
             !
             istate = icontrvib(jlevel)%istate
             icontrvib(jlevel)%istate = icontrvib(ilevel)%istate
             icontrvib(ilevel)%istate = istate
             !
             i = icontrvib(jlevel)%v
             icontrvib(jlevel)%v = icontrvib(ilevel)%v
             icontrvib(ilevel)%v = i
             !
           endif
           !
         enddo
         !
       enddo
       !
       ! print out the vibrational fields in the J=0 representaion
       if (iverbose>=4) then
          write(out,'(/"Vibrational (contracted) energies: ")')
          write(out,'("    i        Energy/cm    State v"/)')
          do i = 1,totalroots
            istate = icontrvib(i)%istate
            write(out,'(i5,f18.6," [ ",2i4," ] ",a)') i,(contrenergy(i)-contrenergy(1))/sc,istate,icontrvib(i)%v, &
                                                      trim(poten(istate)%name)
          enddo
       endif
       !
       if (job%print_vibrational_energies_to_file ) then
          open(unit=u1, file='J0_vibrational_energies.dat',status='replace',action='write')
          write(u1,'(/"Vibrational (contracted) energies: ")')
          write(u1,'("    i        Energy/cm    State v"/)')
          do i = 1,totalroots
            istate = icontrvib(i)%istate
            write(u1,'(i5,f18.6," [ ",2i4," ] ",a)') i,(contrenergy(i))/sc,istate,icontrvib(i)%v, &
                                                      trim(poten(istate)%name)
          enddo
          close(u1)
       endif
       !
       !
       ! check the orthogonality of the basis
       !
       if (iverbose>=3) then
         !
         if (iverbose>=6) write(out,'(/"Check the contracted basis for ortho-normality")')
         !
         if (action%intensity.and.intensity%overlap) then 
           !
           write(out,'(/"Vibrational overlap integrals: ")')
           ! write(out,'("    State-i    <i|j>   State-j"/)')
           write(out,"(1x,a7,1x,a7,6x,a10)") 'State-i','State-j', '<i|j>'
           !
         endif
         !
         !omp parallel do private(ilevel,jlevel,psipsi_t) schedule(guided)
         do ilevel = 1,totalroots
           do jlevel = 1,ilevel
             !
             !
             psipsi_t  = sum(contrfunc(:,ilevel)*contrfunc(:,jlevel))
             !
             if (iverbose>=6) then
                if (icontrvib(ilevel)%istate/=icontrvib(jlevel)%istate.and.ilevel/=jlevel.and.abs(psipsi_t)>sqrt(small_)) then
                   write(out,"('orthogonality is brocken : <',i4,'|',i4,'> (',f16.6,')')") ilevel,jlevel,psipsi_t
                   stop 'Brocken orthogonality'
                endif
                !
                if (ilevel==jlevel.and.abs(psipsi_t-1.0_rk)>sqrt(small_)) then
                   write(out,"('normalization is brocken:  <',i4,'|',i4,'> (',f16.6,')')") ilevel,jlevel,psipsi_t
                   stop 'Brocken normalization'
                endif
             endif
             !
             ! Reporting the quality of the matrix elemenst
             !
             if (action%intensity.and.intensity%overlap.and.&
                 icontrvib(ilevel)%istate/=icontrvib(jlevel)%istate) then
                !
                write(out,'("<",i2,",",i4,"|",i2,",",i4,"> = ",es18.8)') icontrvib(ilevel)%istate,    &
                                                                         icontrvib(ilevel)%v,            &
                                                                         icontrvib(jlevel)%istate,       &
                                                                         icontrvib(jlevel)%v,            &
                                                                         psipsi_t 
             endif
             !
             if (iverbose>=6) then
               if (ilevel/=jlevel) then
                 write(out,"('<',i4,'|',i4,'> = ',e16.2,'<-',8x,'0.0')") ilevel,jlevel,psipsi_t
               else
                 write(out,"('<',i4,'|',i4,'> = ',f16.2,'<-',8x,'1.0')") ilevel,jlevel,psipsi_t
               endif
             endif
             !
           enddo
         enddo
         !omp end parallel do
         !
       endif
       !
       !deallocate(vibmat_rk)
       !call ArrayStop('vibmat_rk')
       !
       !allocate(matelem_rk(totalroots,totalroots),stat=alloc)
       !call ArrayStart('matelem_rk',alloc,size(matelem_rk),kind(matelem_rk))
       !
       deallocate(vec)
       call ArrayStop('vec')
       !
       ! save contrfunc(:, :) 
       vibrational_totalroots = totalroots
       !
       ! keep vibrational basis functions for other applicaitons 
       if (trim(job%basis_set)=='KEEP') then 
         !
         if (.not.allocated(vibrational_contrfunc)) then 
            !
            allocate(vibrational_quantum_number(ngrid*Nestates),stat=alloc)
            allocate(vibrational_contrfunc(ngrid,ngrid*Nestates),stat=alloc)
            call ArrayStart('vibrational_contrfunc',alloc,size(vibrational_contrfunc),kind(vibrational_contrfunc))
            call ArrayStart('vibrational_quantum_number',alloc,ngrid*Nestates*(14+6*2),ik)
            !
         endif
         !
         vibrational_contrfunc = contrfunc(:, 1:totalroots)
         vibrational_quantum_number = icontrvib(1:totalroots)
         !
       endif
       !
       if (iverbose>=4) call TimerStop('Solve vibrational part')
       !
       ! Now we need to compute all vibrational matrix elements of all field of the Hamiltonian, except for the potentials V,
       ! which together with the vibrational kinetic energy operator are diagonal on the contracted basis developed
       !
       ! allocate arrays for matrix elements for all hamiltonian fields
       !
       ! introducing a new field for the centrifugal matrix
       !
       if (.not.fields_allocated) then
         !
         allocate(brot(1),stat=alloc)
         allocate(brot(1)%matelem(totalroots,totalroots),stat=alloc)
         call ArrayStart('brot',alloc,size(brot(1)%matelem),kind(brot(1)%matelem))
         !
       endif
       !
       if (Nnac>0) then 
           !
           vibmat = 0 
           !
           ! f'(0) = [ f(h) - f(-h) ] / (2 h) 
           !kinetic factor is  12*h**2/(2*h) = 6*h 
           !
           do igrid =1, ngrid
             if (igrid>1) then
               vibmat(igrid,igrid-1) = -z(igrid-1)*hstep*6.0_rk
               vibmat(igrid-1,igrid) = -vibmat(igrid,igrid-1)
             endif
           enddo
       endif
       !
       do iobject = 1,Nobjects
          !
          if (iobject==Nobjects-2) cycle
          !
          if ( action%intensity.and.(iobject==Nobjects.and.iverbose>=3.and.(intensity%tdm.or.intensity%tqm)) ) then 
          !if ( iobject==Nobjects-3.and.iverbose>=3.and.action%intensity.and.intensity%tqm) then 
             !
             write(out,'(/"Vibrational transition moments: ")')
!              write(out,'("    State    TM   State"/)')
             write(out,"(A8,A20,25X,A8,A19)") 'State', 'TM', 'State', 'Value'
             !
          endif
          !
          Nmax = fieldmap(iobject)%Nfields
          !
          ! each field type constits of Nmax terms
          !
          do iterm = 1,Nmax
            !
            select case (iobject)
              !
            case (1)
              field => poten(iterm)
            case (2)
              field => spinorbit(iterm)
            case (3)
              field => l2(iterm)
            case (4)
              field => lxly(iterm)
              field%gridvalue(:) = field%gridvalue(:)*b_rot/r(:)**2*sc
            case (5)
              field => spinspin(iterm)
            case (6)
              field => spinspino(iterm)
            case (7)
              field => bobrot(iterm)
              field%gridvalue(:) = field%gridvalue(:)*b_rot/r(:)**2*sc
            case (8)
              field => spinrot(iterm)
            case (9)
              field => diabatic(iterm)
            case (10)
              field => lambdaopq(iterm)
            case (11)
              field => lambdap2q(iterm)
            case (12)
              field => lambdaq(iterm)
            case (13)
              ! A special case of NAC couplings with 1st derivatives wrt r
              field => nac(iterm)
            case (Nobjects-3)
              field => quadrupoletm(iterm)
            case (Nobjects-2)
              field => abinitio(iterm)
            case (Nobjects-1)
              field => brot(iterm)
              field%name = 'BROT'
              if (.not.fields_allocated) then 
                allocate(field%gridvalue(ngrid),stat=alloc)
                call ArrayStart(field%name,alloc,size(field%gridvalue),kind(field%gridvalue))
              endif 
              field%gridvalue(:) = b_rot/r(:)**2*sc
            case (Nobjects)
              !if (.not.action%intensity) cycle 
              field => dipoletm(iterm)
            end select
            !
            ! check if the field wass not allocated at previous call to prevent multiple allocatetions 
            !
            if (.not.fields_allocated) then
              allocate(field%matelem(totalroots,totalroots),stat=alloc)
              call ArrayStart(field%name,alloc,size(field%matelem),kind(field%matelem))
            endif
            !
            field%matelem = 0
            !
            !$omp parallel do private(ilevel,jlevel) schedule(guided)
            do ilevel = 1,totalroots
              do jlevel = 1,ilevel
                 !
                 ! in the grid representation of the vibrational basis set
                 ! the matrix elements are evaluated simply by a summation of over the grid points
                 !
                 !psipsi_ark = real(contrfunc(:,ilevel)*(field%gridvalue(:))*contrfunc(:,jlevel),kind=ark)
                 !
                 !f_ark = simpsonintegral_ark(ngrid-1,psipsi_ark)
                 !
                 !field%matelem(ilevel,jlevel) = f_ark
                 !
                 field%matelem(ilevel,jlevel)  = sum(contrfunc(:,ilevel)*(field%gridvalue(:))*contrfunc(:,jlevel))
                 !
                 ! A special case of the non-diagonal integration of NAC
                 !
                 ! If intensity%threshold%dipole is given and TM is smaller than this threshold set the TM-value to zero
                 ! is applied to the dipole (iobject=Nobjects) and quadrupole (iobject=Nobjects-3) moments 
                 if (iobject==Nobjects-3.or.iobject==Nobjects) then
                   if (abs(field%matelem(ilevel,jlevel))<intensity%threshold%dipole) field%matelem(ilevel,jlevel) = 0 
                 endif
                 !
                 field%matelem(jlevel,ilevel) = field%matelem(ilevel,jlevel)
                 !
                 if (iobject==13) then 
                   !
                   vibener =  matmul(vibmat,contrfunc(1:,jlevel))
                   field%matelem(ilevel,jlevel)  = sum(contrfunc(:,ilevel)*(field%gridvalue(:))*vibener(:))
                   field%matelem(jlevel,ilevel) = -field%matelem(ilevel,jlevel)
                   !
                 endif
                 !
                 !matelem_rk(ilevel,jlevel)  = sum(contrfunc_rk(:,ilevel)*real(field%gridvalue(:),rk)*contrfunc_rk(:,jlevel))
                 !
                 !matelem_rk(ilevel,jlevel)  = sum(contrfunc_rk(:,ilevel)*( (grid_rk(:)-2.24_rk )*0.6_rk )* & 
                 !                                                                        contrfunc_rk(:,jlevel))
                 !
                 !psipsi_rk = contrfunc_rk(:,ilevel)*real(field%gridvalue(:),rk)*contrfunc_rk(:,jlevel)
                 !
                 !psipsi_rk = contrfunc_rk(:,ilevel)*( (grid_rk(:)-2.24_rk )*0.6_rk )*contrfunc_rk(:,jlevel)
                 !
                 !f_rk = simpsonintegral_rk(ngrid-1,psipsi_rk)
                 !
                 !matelem_rk(ilevel,jlevel) = f_rk
                 !
                 !matelem_rk(jlevel,ilevel) = matelem_rk(ilevel,jlevel)
                 !
                 !
              enddo
            enddo
            !$omp end parallel do
            !
            ! printing out transition moments 
            !
            if (  ( iobject==Nobjects  .and.action%intensity.and.intensity%tdm ).or.&
                  ( iobject==Nobjects-3.and.action%intensity.and.intensity%tqm ) ) then
                !
                !write(out,'(/"Vibrational transition moments: ")')
                !write(out,'("    State    TM   State"/)')
                !
                do ilevel = 1,totalroots
                  do jlevel = 1,totalroots
                    !
                    istate = icontrvib(ilevel)%istate
                    jstate = icontrvib(jlevel)%istate
                    !
                    ! dipole selection rules
                    !
                    !if (nint(field%spini-field%spinj)==0.and.abs(field%lambda-field%lambdaj)<=1) then 
                      !
                      !field%matelem(ilevel,jlevel) = field%matelem(ilevel,jlevel)*field%factor
                      !
                      !  if ( iverbose>=4.and.abs(field%matelem(ilevel,jlevel))>sqrt(small_).and.istate==field%istate.and.&
                      if ( iverbose>=4 .and. istate==field%istate.and.&   ! remove the check on magnitude --- print all
                           jstate==field%jstate ) then 
                        !                        NB:   hard limit 40 characters to field name, may lead to truncation!!!
                        write(out,'("<",i2,",",i4,"|",a40,5x,"|",i2,",",i4,"> = ",2es18.8)') icontrvib(ilevel)%istate, &
                                                                                            icontrvib(ilevel)%v,       &
                                                                                            trim(field%name),          &
                                                                                            icontrvib(jlevel)%istate,  &
                                                                                            icontrvib(jlevel)%v,       &
                                                                                            field%matelem(ilevel,jlevel)  ! & 
                                                                                            !,matelem_rk(ilevel,jlevel)
                        !
                      endif
                      !
                    !endif
                    !
                  enddo
                enddo
                !
            endif 
            !
          enddo
          !
       enddo
       !
       fields_allocated = .true.
       !
       !deallocate(contrfunc_rk)
       !call ArrayStop('contrfunc_rk')
       !
       !deallocate(matelem_rk)
       !call ArrayStop('matelem_rk')
       !
       if (allocated(psipsi_ark)) then 
         deallocate(psipsi_ark)
         call ArrayStop('psipsi_ark')
       endif
       !
       ! dealocate some objects
       !
       deallocate(vibmat,vibener)
       call ArrayStop('vibmat')
       call ArrayStop('vibener')
       !
       !deallocate(grid_rk)
       !call ArrayStop('grid_rk')
       !
       ! checkpoint the matrix elements of dipoles if required 
       !
       !if (trim(job%IO_dipole=='SAVE')) then 
       !!    call check_point_dipoles('SAVE',iverbose,totalroots) 
       !endif
       !
     end select 
     !
     ! First we start a loop over J - the total angular momentum quantum number
     !
     if (action%save_eigen_J) then
       ! 
       allocate(eigen(nJ,sym%NrepresCs),basis(nJ),stat=alloc)
       if (alloc/=0) stop 'problem allocating eigen'
       !
       ! initialize the following fields
       do irot = 1,nJ
         do irrep = 1,sym%NrepresCs
           eigen(irot,irrep)%Nlevels = 0
           eigen(irot,irrep)%Ndimen = 0
         enddo
       enddo
       !
     endif
     !
     if (job%IO_eigen=='SAVE') then 
        !
        filename =  trim(job%eigenfile%vectors)//'_vectors.chk'
        write(ioname, '(a, i4)') 'Eigenvectors file '
        call IOstart(trim(ioname),iunit)
        open(unit = iunit, action = 'write',status='replace' , file = filename)
        !
        write(iunit,'(" Molecule = ",a,1x,a)' ) symbol1,symbol2
        write(iunit,'(" masses   = ",2f20.12)') m1,m2
        write(iunit,'(" Nroots   = ",i8)') Nroots
        write(iunit,'(" Nbasis   = ",i8)') totalroots
        write(iunit,'(" Nestates = ",i8)') nestates
        write(iunit,'(" Npoints   = ",i8)') grid%npoints
        write(iunit,'(" range   = ",2f14.7)') grid%rmin,grid%rmax
        !
        do k =1,nestates
          write(iunit, '(a,", ")',advance='no') trim(poten(k)%name)
        enddo
        !
        write(iunit,'(a)') '   <- States'
        !
        write(iunit,"(6x,'|   # |    J | p |           Coeff.   | St vib Lambda Spin     Sigma    Omega ivib|')")
        !
        filename =  trim(job%eigenfile%vectors)//'_vib.chk'
        write(ioname, '(a, i4)') 'Contracted vib basis set on the grid'
        call IOstart(trim(ioname),vibunit)
        open(unit = vibunit, action = 'write',status='replace' , file = filename)
        !
        do i = 1,totalroots
         istate = icontrvib(i)%istate
          write(vibunit,'(i5,f18.6,3x,2i4,3x,a)') i,(contrenergy(i)-contrenergy(1))/sc,icontrvib(i)%istate,icontrvib(i)%v, &
                                                     trim(poten(istate)%name)
          do k = 1,grid%npoints
            write(vibunit,'(e20.12)') contrfunc(k,i)
          enddo
          !
        enddo
        !
        write(vibunit,"('End of contracted basis')")
        if (trim(solution_method)=="LOBATTO") then 
          write(vibunit,"('Start of Lobatto Weights - needed for internal testing in RmatReact code')")
           do k=1, grid%npoints
              write(vibunit,'(e20.12)') LobWeights(k)
          enddo
          write(vibunit,"('End of Lobatto Weights')")
        endif
        !
        close(unit = vibunit, status='keep')
        !
     endif
     !
     if (trim(solution_method)=="LOBATTO") then 
       deallocate(LobDerivs,LobAbs,LobWeights,vibTmat)
       call ArrayStop('LobDerivs')
       call ArrayStop('LobAbs')
       call ArrayStop('LobWeights')
       call ArrayStop('vibTmat')
     endif
     !
     if (allocated(contrfunc)) then 
       deallocate(contrfunc)
       call ArrayStop('contrfunc')
     endif
     !
     if (present(nenerout)) nenerout = 0
     !
     if (present(enerout)) then
       enerout = 0
     endif
     !
     loop_jval : do irot = 1,nJ
       !
       jval = J_list(irot)
       !
       if (jval<jmin) cycle
       !
       if (iverbose>=4) write(out,'(/"j = ",f9.1/)') jval
       !
       ! define the bookkeeping of the quantum numbers for the sigma-lambda basis set
       !
       select case(trim(job%contraction))
         !
       case("OMEGA")
         !
         ! define the bookkeeping of the quantum numbers for the sigma-lambda basis set
         !
         if (iverbose>=3) write(out,'("Define the quanta book-keeping")')
         !
         ! Now we combine together the vibrational-Omega basis functions  and |J,Omega> (as product)
         ! and the corresponding quantum numbers to form our final contracted basis set as well as
         ! the numbering of the contratced basis functions using only one index i.
         !
         ! first count the contracted states of the product basis set
         !
         i = 0
         do ivib =1, totalroots
             !
             omega = icontrvib(ivib)%omega
             !
             ! link the states in Omega-vibrarional and |J,Omega> basis functions
             !
             if (abs(omega)>jval) cycle
             !
             i = i + 1
             !
         enddo
         !
         ! this how many states we get in total after the product of the
         ! vibrational and sigma-lambda basis sets:
         Ntotal = i
         !
         if (Ntotal==0) then
           write(out,'("The size of the rovibronic basis set is zero. Check the CONTRACTION parameters.")')
           stop "The size of the rovibronic basis set is zero"
         endif
         !
         ! allocate the book keeping array to manage the mapping between
         ! the running index i and the vibrational ivib and lamda-sigma ilevel quantum numbers
         allocate(icontr(Ntotal),stat=alloc)
         printout = ''
         !
         if (iverbose>=4) write(out,'(/"Contracted basis set:")')
         if (iverbose>=4) write(out,'("     i     jrot ilevel ivib state v     spin    sigma lambda   omega   Name")')
         !
         ! build the bookkeeping: the object icontr will store this informtion
         !
         i = 0
         do ivib = 1,totalroots
           !
           omega = icontrvib(ivib)%omega
           !
           if (abs(omega)>jval) cycle
           !
           i = i + 1
           !
           ilevel =  icontrvib(ivib)%ilevel 
           iomega = icontrvib(ivib)%iomega
           tau_lambdai = 0 ; if (omega<0) tau_lambdai = 1
    
           istate      =  Omega_grid(iomega)%qn(ilevel)%istate
           sigma       =  Omega_grid(iomega)%qn(ilevel)%sigma 
           ilambda     =  Omega_grid(iomega)%qn(ilevel)%ilambda
           spini       =  Omega_grid(iomega)%qn(ilevel)%spin
           !
           icontr(i) = omega_grid(iomega)%qn(ilevel)
           icontr(i)%ivib = ivib
           icontr(i)%ilevel = ilevel
           icontr(i)%v = icontrvib(ivib)%v
           icontr(i)%ilevel = ilevel
           icontr(i)%iomega = iomega
           icontr(i)%omega = omega
           icontr(i)%spin = Omega_grid(iomega)%qn(ilevel)%spin
           !
           ! print the quantum numbers
           if (iverbose>=4) then
               write(out,'(i6,1x,f8.1,1x,i4,1x,i4,1x,i4,1x,i4,1x,f8.1,1x,f8.1,1x,i4,1x,f8.1,3x,a)') &
                       i,jval,ilevel,ivib,ilevel,&
                       icontr(i)%v,spini,sigma,ilambda,omega,trim(poten(istate)%name)
           endif
           !
         enddo
         !
         ! allocate the hamiltonian matrix and an array for the energies of this size Ntotal
         allocate(hmat(Ntotal,Ntotal),stat=alloc)
         call ArrayStart('hmat',alloc,size(hmat),kind(hmat))
         !
         if (iverbose>=4) call MemoryReport
         !
         call Compute_rovibronic_Hamiltonian_in_omega_vib_representation(iverbose,jval,ngrid,Ntotal,Nomega_states,&
                                                                         sc,icontr,contrenergy,hmat)
         !
         ! Transformation to the symmetrized basis set
         !
         ! |v,Lambda,Sigma,J,Omega,tau> = 1/sqrt(2) [ |v,Lambda,Sigma,J,Omega>+(-1)^tau |v,-Lambda,-Sigma,J,-Omega> ]
         !
         allocate(iswap(Ntotal),vec(Ntotal),Nirr(Ntotal,2),ilevel2i(Ntotal,2),tau(Ntotal),ilevel2isym(Ntotal,2),stat=alloc)
         call ArrayStart('iswap-vec',alloc,size(iswap),kind(iswap))
         call ArrayStart('iswap-vec',alloc,size(vec),kind(vec))
         call ArrayStart('Nirr',alloc,size(Nirr),kind(Nirr))
         !
         iswap = 0
         Nsym = 0
         Nlevels = 0
         ilevel2i = 0
         ilevel2isym = 0
         Nirr = 0
         !
         !omp parallel do private(istate,sigmai,ilambda,spini,omegai,ibib,j,jstate,sigmaj,jlambda,omegaj,spinj,jvib) & 
         !                        shared(iswap,vec) schedule(guided)
         do i = 1,Ntotal
           !
           istate  = icontr(i)%istate
           sigmai  = icontr(i)%sigma
           ilambda = icontr(i)%ilambda
           spini   = icontr(i)%spin
           omegai  = icontr(i)%omega
           ilevel  = icontr(i)%ilevel
           ivib    = icontr(i)%v
           !
           if (iswap(i)/=0) cycle
           !
           if (nint(omegai)==0) then
             !
             Nlevels = Nlevels + 1
             !
             itau = nint(Jval)
             !
             itau = abs(itau)
             !
             tau(Nlevels) = (-1.0_rk)**itau
             !
             if (mod(itau,2)==1) then
               Nsym(2) = Nsym(2) + 1
               ilevel2isym(Nlevels,2) = nsym(2)
               ilevel2i(Nlevels,2) = i
               Nirr(Nlevels,2) = 1
             else
               Nsym(1) = Nsym(1) + 1
               ilevel2isym(Nlevels,1) = nsym(1)
               ilevel2i(Nlevels,1) = i
               Nirr(Nlevels,1) = 1
             endif
             !
             iswap(i) = i
             !ilevel2i(Nlevels,1) = i
             !
           else
             !
             do  j = 1,Ntotal
               !
               jstate  = icontr(j)%istate
               sigmaj  = icontr(j)%sigma
               jlambda = icontr(j)%ilambda
               omegaj  = icontr(j)%omega
               jlevel  = icontr(j)%ilevel
               spinj   = icontr(j)%spin
               jvib    = icontr(j)%v
               !
               if (nint(omegai+omegaj)==0.and. ilevel==jlevel.and.ivib==jvib) then
                 !
                 ! total number of irreps
                 !
                 Nsym(:) = Nsym(:) + 1
                 !
                 ! Levels of each irrep
                 !
                 Nlevels = Nlevels + 1
                 !
                 Nirr(Nlevels,:) = 1
                 !
                 !ilevel2i(Nlevels,1) = i
                 !ilevel2i(Nlevels,2) = j
                 !
                 if (omegaj>omegai) then
                 !if (mod(nint(omegai-omegaj),2)==0) then
                   !
                   ilevel2i(Nlevels,1) = i
                   ilevel2i(Nlevels,2) = j
                   !
                 else
                   !
                   ilevel2i(Nlevels,1) = j
                   ilevel2i(Nlevels,2) = i
                   !
                 endif
                 !
                 ilevel2isym(Nlevels,1:2) = nsym(1:2)
                 !
                 itau = nint(Jval)
                 !
                 if (mod(nint(2.0_rk*Jval),2)/=0) itau = nint(Jval+0.5_rk+2.0_rk*omegai)
                 !
                 tau(Nlevels) = (-1.0_rk)**itau
                 !
                 iswap(i) = j*(-1)**itau
                 iswap(j) = i*(-1)**itau
                 !
                 exit
                 !
               endif
               !
             enddo
             !
           endif
           !
         enddo
         !omp end parallel do
         !
         ! Nlevels is the number of states disregarding the degeneracy 
         ! Nroots is the total number of roots including the degenerate states 
         !
         allocate(transform(1)%matrix(max(1,Nsym(1)),max(1,Nsym(1))),stat=alloc)
         allocate(transform(2)%matrix(max(1,Nsym(2)),max(1,Nsym(2))),stat=alloc)
         allocate(transform(1)%irec( max( 1,Nsym(1) ) ),stat=alloc)
         allocate(transform(2)%irec( max( 1,Nsym(2) ) ),stat=alloc)
         !
         call ArrayStart('transform',alloc,size(transform(1)%matrix),kind(transform(1)%matrix))
         call ArrayStart('transform',alloc,size(transform(2)%matrix),kind(transform(2)%matrix))
         call ArrayStart('transform',alloc,size(transform(1)%irec),kind(transform(1)%irec))
         call ArrayStart('transform',alloc,size(transform(2)%irec),kind(transform(2)%irec))
         !
         allocate(Utransform(Nlevels,2,2),stat=alloc)
         call ArrayStart('Utransform',alloc,size(Utransform),kind(Utransform))
         !
         ! Building the transformation to the symmetrized representaion (irrep)
         !
         do ilevel = 1,Nlevels
           !
           veci = 0
           !
           if (any(Nirr(ilevel,:)==0)) then
             !
             do irrep = 1,sym%NrepresCs
               do itau = 1,Nirr(ilevel,irrep)
                 !
                 i = ilevel2i(ilevel,irrep)
                 veci(irrep,irrep) = 1.0_rk
                 isym = ilevel2isym(ilevel,irrep)
                 transform(irrep)%irec(isym) = ilevel
                 !
               enddo
             enddo
             !
             !
             !if (tau(ilevel)<0) then
             !  veci(1,2) = 1.0_rk
             !  isym = ilevel2isym(ilevel,2)
             !  transform(2)%irec(isym) = ilevel
             !else
             !  veci(1,1) = 1.0_rk
             !  isym = ilevel2isym(ilevel,1)
             !  transform(1)%irec(isym) = ilevel
             !endif
             !
           else
             !
             veci(1,1) = sqrt(0.5_rk)
             veci(2,1) = sqrt(0.5_rk)*tau(ilevel)
             veci(1,2) = sqrt(0.5_rk)
             veci(2,2) =-sqrt(0.5_rk)*tau(ilevel)
             !
             do irrep = 1,sym%NrepresCs
                isym = ilevel2isym(ilevel,irrep)
                transform(irrep)%irec(isym) = ilevel
             enddo
             !
           endif
           !
           Utransform(ilevel,:,:) = veci(:,:)
           !
           do jlevel = 1,Nlevels
             !
             vecj = 0
             !
             if (any(Nirr(jlevel,:)==0)) then
                !
                do jrrep = 1,sym%NrepresCs
                  do jtau = 1,Nirr(jlevel,jrrep)
                    vecj(jrrep,jrrep) = 1.0_rk
                  enddo
                enddo
                !
             else
                !
                vecj(1,1) = sqrt(0.5_rk)
                vecj(2,1) = sqrt(0.5_rk)*tau(jlevel)
                vecj(1,2) = sqrt(0.5_rk)
                vecj(2,2) =-sqrt(0.5_rk)*tau(jlevel)
                !
             endif
             !
             pmat = 0
             !
             do isym = 1,2
                do itau = 1,Nirr(ilevel,isym)
                   i = ilevel2i(ilevel,isym)
                   do jsym = 1,2
                      do jtau = 1,Nirr(jlevel,jsym)
                         j = ilevel2i(jlevel,jsym)
                         !
                         if (i<=j) then
                           pmat(isym,jsym) = hmat(i,j)
                         else
                           pmat(isym,jsym) = hmat(j,i)
                         endif
                         !
                      enddo
                   enddo
                enddo
             enddo
             !
             smat = matmul(transpose(veci),matmul(pmat,vecj))
             !
             !smat = matmul((veci),matmul(pmat,transpose(vecj)))
             !
             do irrep = 1,sym%NrepresCs
                do itau = 1,Nirr(ilevel,irrep)
                   !i = ilevel2i(ilevel,isym)
                   !
                   isym = ilevel2isym(ilevel,irrep)
                   !
                   do jrrep = 1,sym%NrepresCs
                      do jtau = 1,Nirr(jlevel,jrrep)
                         !j = ilevel2i(jlevel,jsym)
                         jsym = ilevel2isym(jlevel,jrrep)
                         !
                         if (irrep==jrrep) then
                           !
                           transform(irrep)%matrix(isym,jsym) = smat(irrep,irrep)
                           !
                         else
                           !
                           if (abs(smat(irrep,jrrep))>sqrt(small_)) then
                             !
                             i = ilevel2i(ilevel,itau)
                             j = ilevel2i(jlevel,jtau)
                             !
                             istate = icontr(i)%istate
                             sigmai  = icontr(i)%sigma
                             ilambda = icontr(i)%ilambda
                             spini   = icontr(i)%spin
                             omegai  = icontr(i)%omega
                             ivib    = icontr(i)%v
                             !
                             jstate  = icontr(j)%istate
                             sigmaj  = icontr(j)%sigma
                             jlambda = icontr(j)%ilambda
                             omegaj  = icontr(j)%omega
                             spinj   = icontr(j)%spin
                             jvib    = icontr(j)%v
                             !
                             write(out,'(/"Problem with symmetry: The non-diagonal matrix element is not zero:")')
                             write(out,'(/"i,j = ",2i8," irrep,jrrep = ",2i8," isym,jsym = ",2i8," ilevel,jlevel = ", &
                                        & 2i3," , matelem =  ",g16.9," with zero = ",g16.9)') &
                                        i,j,irrep,jrrep,isym,jsym,ilevel,jlevel,smat(itau,jtau),sqrt(small_)
                             write(out,'(/"<State   v  lambda spin   sigma  omega |H(sym)| State   v  lambda spin   '//&
                               'sigma  omega>")')
                             write(out,'("<",i3,2x,2i4,3f8.1," |H(sym)| ",i3,2x,2i4,3f8.1,"> /= 0")') &
                                         istate,ivib,ilambda,spini,sigmai,omegai,jstate,jvib,jlambda,spinj,sigmaj,omegaj
                             write(out,'("<",a10,"|H(sym)|",a10,"> /= 0")') trim(poten(istate)%name),trim(poten(jstate)%name)
                             !
                             stop 'Problem with symmetry: The non-diagonal matrix element is not zero'
                             !
                           endif
                           !
                         endif
                         !
                      enddo
                   enddo
                enddo
             enddo
    
           enddo
           !
         enddo
         !
       case('VIB')
         !
         if (iverbose>=3) write(out,'("Define the quanta book-keeping")')
         !
         call define_quanta_bookkeeping(iverbose,jval,Nestates,Nlambdasigmas)
         !
         if (iverbose>=3) write(out,'("...done!")')
         !
         ! Now we combine together the vibrational and sigma-lambda basis functions (as product)
         ! and the corresponding quantum numbers to form our final contracted basis set as well as
         ! the numbering of the contratced basis functions using only one index i.
         !
         ! first count the contracted states of the product basis set
         !
         i = 0
         do ilevel = 1,Nlambdasigmas
           do ivib =1, totalroots
             !
             if (quanta(ilevel)%istate/=icontrvib(ivib)%istate) cycle
             i = i + 1
           enddo
         enddo
         !
         ! this how many states we get in total after the product of the
         ! vibrational and sigma-lambda basis sets:
         Ntotal = i
         !
         if (Ntotal==0) then
           write(out,'("The size of the rovibronic basis set is zero. Check the CONTRACTION parameters.")')
           stop "The size of the rovibronic basis set is zero"
         endif
         !
         ! allocate the book keeping array to manage the mapping between
         ! the running index i and the vibrational ivib and lamda-sigma ilevel quantum numbers
         allocate(icontr(Ntotal),printout(Nlambdasigmas),stat=alloc)
         printout = ''
         !
         allocate(ilambdasigmas_v_icontr(totalroots,Nlambdasigmas),stat=alloc)
         call ArrayStart('ilambdasigmas_v_icontr',alloc,size(ilambdasigmas_v_icontr),kind(ilambdasigmas_v_icontr))
         !
         ilambdasigmas_v_icontr = 0 
         !
         if (iverbose>=4) write(out,'(/"Contracted basis set:")')
         if (iverbose>=4) write(out,'("     i     jrot ilevel ivib state v     spin    sigma lambda   omega   Name")')
         !
         ! build the bookkeeping: the object icontr will store this informtion
         !
         i = 0
         do ilevel = 1,Nlambdasigmas
           !
           istate = quanta(ilevel)%istate
           sigma = quanta(ilevel)%sigma
           imulti = quanta(ilevel)%imulti
           ilambda = quanta(ilevel)%ilambda
           omega = quanta(ilevel)%omega
           spini = quanta(ilevel)%spin
           tau_lambdai = 0 ; if (ilambda<0) tau_lambdai = 1
           !
           do ivib =1,totalroots
             !
             ! link the states in vibrarional and spin-rot basis components for the VIB contraction
             !
             if (quanta(ilevel)%istate/=icontrvib(ivib)%istate) cycle
             !
             i = i + 1
             !
             icontr(i) = quanta(ilevel)
             icontr(i)%ivib = ivib
             icontr(i)%ilevel = ilevel
             icontr(i)%v = icontrvib(ivib)%v
             !
             ilambdasigmas_v_icontr(ivib,ilevel) = i
             !
             ! print the quantum numbers
             if (iverbose>=4) then
                 write(out,'(i6,1x,f8.1,1x,i4,1x,i4,1x,i4,1x,i4,1x,f8.1,1x,f8.1,1x,i4,1x,f8.1,3x,a)') &
                         i,jval,ilevel,ivib,istate,&
                         icontr(i)%v,spini,sigma,ilambda,omega,trim(poten(istate)%name)
             endif
             !
           enddo
         enddo
         !
         ! allocate the hamiltonian matrix and an array for the energies of this size Ntotal
         allocate(hmat(Ntotal,Ntotal),stat=alloc)
         call ArrayStart('hmat',alloc,size(hmat),kind(hmat))
         !
         if (iverbose>=4) call MemoryReport
         !
         hmat = 0
         !
         if (iverbose>=4) call TimerStart('Construct the hamiltonian')
         !
         if (iverbose>=3) write(out,'(/"Construct the hamiltonian matrix")')
         !
         !omp parallel do private(i,ivib,ilevel,istate,sigmai,imulti,ilambda,omegai,spini,v_i,jvib,jlevel,jstate,sigmaj, & 
         !                        jmulti,jlambda,omegaj,spinj,v_j,f_rot,erot,iL2,field,f_l2,f_s,f_t,iso,ibraket,ipermute,&
         !                        istate_,ilambda_,sigmai_,spini_,jstate_,jlambda_,sigmaj_,spinj_,isigmav,omegai_,       &
         !                        omegaj_,itau,ilxly,f_grid,f_l,f_ss) shared(hmat) schedule(guided)
         do i = 1,Ntotal
           !
           ivib = icontr(i)%ivib
           ilevel = icontr(i)%ilevel
           !
           istate = icontr(i)%istate
           sigmai = icontr(i)%sigma
           imulti = icontr(i)%imulti
           ilambda = icontr(i)%ilambda
           omegai = icontr(i)%omega
           spini = icontr(i)%spin
           v_i = icontr(i)%v
           !
           ! the diagonal contribution is the energy from the contracted vibrational solution
           !
           hmat(i,i) = contrenergy(ivib)
           !
           if (trim(poten(istate)%integration_method)=="NUMEROV") cycle
           !
           do j =i,Ntotal
              !
              jvib = icontr(j)%ivib
              jlevel = icontr(j)%ilevel
              jstate = icontr(j)%istate
              sigmaj = icontr(j)%sigma
              jmulti = icontr(j)%imulti
              jlambda = icontr(j)%ilambda
              omegaj = icontr(j)%omega
              spinj = icontr(j)%spin
              v_j = icontr(j)%v
              !
              if (iverbose>=6) write(out,'("ilevel,ivib = ",2(i0,2x) )') ilevel,ivib
              !
              ! the centrifugal factor will be needed for different terms
              !
              f_rot=brot(1)%matelem(ivib,jvib)
              !
              ! BOB centrifugal (rotational) term, i.e. a correction to f_rot
              !
              do ibobrot = 1,Nbobrot
                if (bobrot(ibobrot)%istate==istate.and.bobrot(ibobrot)%jstate==jstate.and.istate==jstate) then
                  field => bobrot(ibobrot)
                  f_bobrot = field%matelem(ivib,jvib)
                  f_rot = f_rot + f_bobrot
                  exit
                endif
              enddo
              !
              ! For the raw non-integrated basis add the kinetic energy matrix elemenent which otherwsie is not included in the 
              ! corresponding contracted enegy values 
              !
              if (istate==jstate.and.ilevel==jlevel) then
                !
                select case (trim(poten(istate)%integration_method))
                !
                case ('NONE','RAW')
                  hmat(i,j) = hmat(i,j) + kinmat(v_i+1,v_j+1)
                end select 
                !
              endif
              !
              if (action%RWF) then
                if (istate>Nrefstates.or.jstate>Nrefstates) cycle
              endif
              !
              if (trim(poten(jstate)%integration_method)=="NUMEROV") cycle
              !
              ! diagonal elements
              !
              if (ilevel==jlevel) then
                !                                             ! L Lodi -job%diag_L2_fact is either zero or one
                erot = f_rot*( Jval*(Jval+1.0_rk) - omegai**2 -job%diag_L2_fact*real(ilambda**2,rk)  & 
                         +   spini*(spini+1.0_rk) - sigmai**2 )
                !
                ! add the diagonal matrix element to the local spin-rotational matrix hmat
                hmat(i,j) = hmat(i,j) + erot
                !
                ! Diagonal spin-spin term
                !
                do iss = 1,Nss
                  if (spinspin(iss)%istate==istate.and.spinspin(iss)%jstate==jstate.and.istate==jstate) then
                    field => spinspin(iss)
                    f_ss = field%matelem(ivib,jvib)*(3.0_rk*sigmai**2-spini*(spini+1.0_rk))*sc
                    hmat(i,j) = hmat(i,j) + f_ss
                    exit
                  endif
                enddo
                !
                ! Diagonal spin-rotation term
                !
                do isr = 1,Nsr
                  if (spinrot(isr)%istate==istate.and.spinrot(isr)%jstate==jstate.and.istate==jstate) then
                    field => spinrot(isr)
                    f_sr = field%matelem(ivib,jvib)*(sigmai**2-spini*(spini+1.0_rk))*sc
                    hmat(i,j) = hmat(i,j) + f_sr
                    exit
                  endif
                enddo
                !
                ! print out the internal matrix at the first grid point
                if (iverbose>=4.and.abs(hmat(i,j)) >small_) then
                    write(printout(ilevel),'(A, F15.3,A)') "RV=", hmat(i,j)/sc, "; "
                endif
                !
              endif
              !
              ! Non diagonal L2 term
              !
              do iL2 = 1,Nl2
                if (L2(iL2)%istate==istate.and.L2(iL2)%jstate==jstate.and.istate/=jstate) then
                  field => L2(iL2)
                  f_l2 = f_rot*field%matelem(ivib,jvib)
                  hmat(i,j) = hmat(i,j) + f_l2
                  exit
                endif
              enddo
              !
              ! Diabatic non-diagonal contribution  term
              !
              do idiab = 1,Ndiabatic
                !
                if (diabatic(idiab)%istate==istate.and.diabatic(idiab)%jstate==jstate.and.&
                    abs(nint(sigmaj-sigmai))==0.and.(ilambda==jlambda).and.nint(spini-spinj)==0 ) then
                  field => diabatic(idiab)
                  f_diabatic = field%matelem(ivib,jvib)*sc
                  hmat(i,j) = hmat(i,j) + f_diabatic
                  exit
                endif
              enddo
              !
              ! NAC non-diagonal contribution  term
              !
              do iNAC = 1,Nnac
                if (nac(iNAC)%istate==istate.and.nac(iNAC)%jstate==jstate.and.&
                    abs(nint(sigmaj-sigmai))==0.and.(ilambda==jlambda).and.nint(spini-spinj)==0 ) then
                  field => nac(iNAC) 
                  f_nac = (field%matelem(ivib,jvib)-field%matelem(jvib,ivib))*sc
                  hmat(i,j) = hmat(i,j) + f_nac
                  exit
                endif
              enddo
              !
              !  Non-diagonal spin-spin term
              !
              loop_iss : do iss = 1,Nss
                !
                if ( spinspin(iss)%istate/=istate.or.spinspin(iss)%jstate/=jstate.or.&
                   spinspin(iss)%istate==spinspin(iss)%jstate ) cycle
                !
                field => spinspin(iss)
              
                ! The selection rules are (Lefebvre-Brion and Field, Eq. (3.4.50)): 
                ! Delta J = 0 ; Delta Omega  = 0 ; g<-/->u; e<->f; Sigma+<->Sigma-;
                ! Delta S = 0 or Delta S = 1 ; Delta Lambda = Delta Sigma = 0 or Delta Lambda = - Delta Sigma = +/- 1
                !
                if (nint(omegai-omegaj)/=0.or.nint(spini-spinj)>2 ) cycle
                !if ( ilambda==0.and.jlambda==0.and.poten(istate)%parity%pm==poten(jstate)%parity%pm ) cycle
                !if ( poten(istate)%parity%gu/=0.and.poten(istate)%parity%gu/=poten(jstate)%parity%gu ) cycle
                !
                do ipermute  = 0,1
                  !
                  if (ipermute==0) then
                    !
                    istate_ = field%istate ; ilambda_we = field%lambda  ; sigmai_we = field%sigmai ; spini_ = field%spini
                    jstate_ = field%jstate ; jlambda_we = field%lambdaj ; sigmaj_we = field%sigmaj ; spinj_ = field%spinj
                    !
                  else  ! permute
                    !
                    jstate_ = field%istate ; jlambda_we = field%lambda  ; sigmaj_we = field%sigmai ; spinj_ = field%spini
                    istate_ = field%jstate ; ilambda_we = field%lambdaj ; sigmai_we = field%sigmaj ; spini_ = field%spinj
                    !
                  endif
                  ! proceed only if the spins of the field equal the corresponding <i| and |j> spins of the current matrix elements. 
                  ! otherwise skip it:
                  if ( nint(spini_-spini)/=0.or.nint(spinj_-spinj)/=0 ) cycle
                  !
                  ! however the permutation makes sense only when for non diagonal <State,Lambda,Spin|F|State',Lambda',Spin'>
                  ! otherwise it will cause a double counting:
                  !
                  if (ipermute==1.and.istate_==jstate_.and.ilambda_we==jlambda_we.and.nint(sigmai_we-sigmaj_we)==0.and. & 
                      nint(spini_-spinj_)==0) cycle
                  !
                  ! check if we are at the right electronic states
                  if( istate/=istate_.or.jstate/=jstate_ ) cycle
                  !
                  ! We apply the Wigner-Eckart theorem to reconstruct all combinations of <Lamba Sigma |HSS|Lamba Sigma' > 
                  ! connected with the reference (input) <Lamba Sigma_r |HSS|Lamba Sigma_r' > by this theorem. 
                  ! Basically, we loop over Sigma (Sigma = -S..S).  The following 3j-symbol for the reference states will 
                  ! be conidered:
                  ! / Si      k  Sj     \    k  = 2
                  ! \ -Sigmai q  Sigmaj /    q  = Sigmai - Sigmaj
                  !
                  ! reference q from Wigner-Eckart
                  q_we = sigmai_we-sigmaj_we
                  !
                  ! We should consider also a permutation <State',Lambda',Spin'|F|State,Lambda,Spin> if this makes a change.
                  ! This could be imortant providing that we constrain the i,j indexes to be i<=j (or i>=j).
                  ! We also assume that the matrix elements are real!
                  !
                  ! First of all we can check if the input values are not unphysical and consistent with Wigner-Eckart:
                  ! the corresponding three_j should be non-zero:
                  three_j_ref = three_j(spini_, 2.0_rk, spinj_, -sigmai_we, q_we, sigmaj_we)
                  !
                  if (abs(three_j_ref)<small_) then 
                    !
                    write(out,"('The Spin-orbit field ',2i3,' is incorrect according to Wigner-Eckart, three_j = 0 ')") & 
                          field%istate,field%jstate
                    write(out,"('Check S_i, S_j, Sigma_i, Sigma_j =  ',4f9.2)") spini_,spinj_,sigmai_we,sigmaj_we
                    stop "The S_i, S_j, Sigma_i, Sigma_j are inconsistent"
                    !
                  end if 
                  !
                  ! Also check the that the SO is consistent with the selection rules for SS
                  !
                  if ( ilambda_we-jlambda_we+nint(sigmai_we-sigmaj_we)/=0.or.nint(spini_-spinj_)>2.or.&
                     ( ilambda_we==0.and.jlambda_we==0.and.poten(field%istate)%parity%pm/=poten(field%jstate)%parity%pm ).or.&
                     ( (ilambda_we-jlambda_we)/=-nint(sigmai_we-sigmaj_we) ).or.&
                        abs(ilambda_we-jlambda_we)>2.or.abs(nint(sigmai_we-sigmaj_we))>2.or.&
                     ( poten(field%istate)%parity%gu/=0.and.poten(field%istate)%parity%gu/=poten(field%jstate)%parity%gu ) ) then
                     !
                     write(out,"('The quantum numbers of the spin-spin field ',2i3,' are inconsistent" // &
                                     " with SO selection rules: ')") field%istate,field%jstate
                     write(out,"('Delta J = 0 ; Delta Omega  = 0 ; g<-/->u; e<-/->f; Sigma+<->Sigma-; " // &
                        "Delta S = 0 or Delta S = 1,2 ; Delta Lambda = Delta Sigma = 0 or Delta Lambda = - Delta Sigma = +/- 2')")
                     write(out,"('Check S_i, S_j, Sigma_i, Sigma_j, lambdai, lambdaj =  ',4f9.2,2i4)") &
                                                                        spini_,spinj_,sigmai_we,sigmaj_we,ilambda_we,jlambda_we
                     stop "The S_i, S_j, Sigma_i, Sigma_j lambdai, lambdaj are inconsistent with selection rules"
                     !
                  endif
                  !
                  do isigma2 = -nint(2.0*spini_),nint(2.0*spini_),2
                    !
                    ! Sigmas from Wigner-Eckart
                    sigmai_ = real(isigma2,rk)*0.5 
                    sigmaj_ = sigmai_ - q_we
                    !
                    ! three_j for current Sigmas
                    three_j_ = three_j(spini_, 2.0_rk, spinj_, -sigmai_, q_we, sigmaj_)
                    !
                    ! current value of the SO-matrix element from Wigner-Eckart
                    SO = (-1.0_rk)**(sigmai_-sigmai_we)*three_j_/three_j_ref*field%matelem(ivib,jvib)
                    !
                    ! We should also take into account that Lambda and Sigma can change sign
                    ! since in the input we give only a unique combination of matrix elements, for example
                    ! < 0 0 |  1  1 > is given, but < 0 0 | -1 -1 > is not, assuming that the program will generate the missing
                    ! combinations.
                    !
                    ! In order to recover other combinations we apply the symmetry transformation
                    ! laboratory fixed inversion which is equivalent to the sigmav operation 
                    !                    (sigmav= 0 correspond to the unitary transformation)
                    do isigmav = 0,1
                      !
                      ! sigmav is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
                      ! avoid the double counting.
                      if( isigmav==1.and.nint( abs( 2.0*sigmai_ )+ abs( 2.0*sigmaj_ ) )+abs( ilambda_we )+abs( jlambda_we )==0 ) &
                                                                                                                             cycle
                      !
                      ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                      ilambda_ = ilambda_we*(-1)**isigmav
                      jlambda_ = jlambda_we*(-1)**isigmav
                      sigmai_ = sigmai_*(-1.0_rk)**isigmav
                      sigmaj_ = sigmaj_*(-1.0_rk)**isigmav
                      !
                      omegai_ = sigmai_+real(ilambda_)
                      omegaj_ = sigmaj_+real(jlambda_)
                      !
                      ! Check So selection rules
                      if ( ( ilambda_-jlambda_)/=-nint(sigmai_-sigmaj_).or.abs(sigmai_-sigmaj_)>1.or.omegai_/=omegaj_ ) cycle
                      !
                      ! proceed only if the quantum numbers of the field equal
                      ! to the corresponding <i| and |j> quantum numbers of the basis set. otherwise skip it:
                      if ( nint(sigmai_-sigmai)/=0.or.nint(sigmaj_-sigmaj)/=0.or.ilambda_/=ilambda.or.jlambda_/=jlambda ) cycle
                      !
                      f_ss = SO*sc
                      !
                      ! the result of the symmetry transformtaion applied to the <Lambda,Sigma|HSO|Lambda',Sigma'> only
                      if (isigmav==1) then
                        !
                        ! still not everything is clear here: CHECK!
                        !
                        itau = -ilambda_-jlambda_ +nint(spini_-sigmai_)+nint(spinj_-sigmaj_) !+nint(jval-omegai)+(jval-omegaj)
                        !
                        !itau = nint(spini_-sigmai_)+nint(spinj_-sigmaj_) ! +nint(jval-omegai)+(jval-omegaj)
                        !
                        !itau = 0
                        !
                        if (ilambda_==0.and.poten(istate)%parity%pm==-1) itau = itau+1
                        if (jlambda_==0.and.poten(jstate)%parity%pm==-1) itau = itau+1
                        !
                        f_ss = f_ss*(-1.0_rk)**(itau)
                        !
                      endif
                      !
                      ! double check
                      if ( nint(omegai-omegai_)/=0 .or. nint(omegaj-omegaj_)/=0 ) then
                        write(out,'(A,f8.1," or omegaj ",f8.1," do not agree with stored values ",f8.1,1x,f8.1)') &
                                   "SO: reconsrtucted omegai", omegai_,omegaj_,omegai,omegaj
                        stop 'SO: wrongly reconsrtucted omegai or omegaj'
                      endif
                      !
                      ! we might end up in eilther parts of the matrix (upper or lower),
                      ! so it is safer to be general here and
                      ! don't restrict to lower part as we have done above
                      !
                      hmat(i,j) = hmat(i,j) + f_ss
                      !
                      !hmat(j,i) = hmat(i,j)
                      !
                      ! print out the internal matrix at the first grid point
                      if (iverbose>=4.and.abs(hmat(i,j))>small_) then
                          !
                          write(printout_,'(i3,"-SS",2i3)') iss,ilevel,jlevel
                          printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                          write(printout_,'(g12.4)') f_ss/sc
                          printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                         !
                      endif
                      !
                      cycle loop_iss
                      !
                    enddo
                  enddo
                enddo              
                !
              enddo loop_iss
              !
              do isso = 1,Nsso
                !
                if (spinspino(isso)%istate==istate.and.spinspino(isso)%jstate==jstate.and.&!istate==jstate.and.&
                    abs(nint(sigmaj-sigmai))==1.and.(ilambda-jlambda)==nint(sigmaj-sigmai)) then
                   !
                   field => spinspino(isso)
                   !
                   f_s = sigmaj-sigmai
                   !
                   !f_t = sqrt( (spinj-f_s*sigmaj )*( spinj + f_s*sigmaj+1.0_rk ) )*&
                   !      sqrt( (jval -f_s*omegaj )*( jval  + f_s*omegaj+1.0_rk ) )
                   !
                   f_t = sqrt( spini*(spini+1.0_rk)-(sigmai+0.5_rk*f_s)*(sigmai    ) )*&
                         sqrt( spini*(spini+1.0_rk)-(sigmai+0.5_rk*f_s)*(sigmai+f_s) )
                   !
                   f_ss = field%matelem(ivib,jvib)*f_t*sc
                   !
                   hmat(i,j) = hmat(i,j) + f_ss
                   !
                   ! print out the internal matrix at the first grid point
                   if (iverbose>=4.and.abs(hmat(i,j))>sqrt(small_)) then
                      write(printout_,'("    SS-o",2i3)') ilevel,jlevel
                      printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                      if (abs(hmat(i,j))>sqrt(small_)) then
                        write(printout_,'(g12.4)') hmat(i,j)/sc
                        printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                      endif
                   endif
                   !
                endif
                ! 
              enddo
              !
              ! Non-diagonal spin-rotaion term
              !
              do isr = 1,Nsr
                !
                field => spinrot(isr)
                !
                ! Two options are possible: 
                ! 1. <Sigma,Omega,Lambda|HSR|Sigma+/-1,Omega+/-1,Lambda>
                ! 2. <Sigma,Omega,Lambda|HSR|Sigma+/-1,Omega,Lambda-/+>
                !
                ! 1. <Sigma,Omega,Lambda|HSR|Sigma+/-1,Omega+/-1,Lambda>
                if (spinrot(isr)%istate==istate.and.spinrot(isr)%jstate==jstate.and.istate==jstate.and.&
                    abs(nint(sigmaj-sigmai))==1.and.(ilambda==jlambda).and.nint(spini-spinj)==0) then
                   !
                   f_s = sigmaj-sigmai
                   !
                   f_t = sqrt( jval* (jval +1.0_rk)-omegai*(omegai+f_s) )*&
                         sqrt( spini*(spini+1.0_rk)-sigmai*(sigmai+f_s) )
                   !
                   f_sr = field%matelem(ivib,jvib)*f_t*sc
                   !
                   hmat(i,j) = hmat(i,j) + f_sr*0.5_rk
                   !
                   ! print out the internal matrix at the first grid point
                   if (iverbose>=4.and.abs(hmat(i,j))>sqrt(small_)) then
                      write(printout_,'("    SR",2i3)') ilevel,jlevel
                      printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                      if (abs(hmat(i,j))>sqrt(small_)) then
                        write(printout_,'(g12.4)') hmat(i,j)/sc
                        printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                      endif
                   endif
                   !
                endif
                !
                ! 2. <Sigma,Omega,Lambda|HSR|Sigma+/-1,Omega,Lambda-/+>
                ! with the effective parameter gamma_v including the matrix element <Lambda|L+/-|lambda-/+1>
                if (spinrot(isr)%istate==istate.and.spinrot(isr)%jstate==jstate.and.&
                    abs(nint(sigmaj-sigmai))==1.and.abs(ilambda-jlambda)==1.and.nint(spini-spinj)==0) then
                    !
                    do ipermute  = 0,1
                      !
                      if (ipermute==0) then
                        !
                        istate_ = field%istate ; ilambda_ = field%lambda  
                        jstate_ = field%jstate ; jlambda_ = field%lambdaj 
                        !
                      else  ! permute
                        !
                        jstate_ = field%istate ; jlambda_ = field%lambda 
                        istate_ = field%jstate ; ilambda_ = field%lambdaj
                        !
                      endif
                      !
                      ! however the permutation makes sense only when for non diagonal <State,Lambda,Spin|F|State',Lambda',Spin'>
                      ! otherwise it will cause a double counting:
                      !
                      if (ipermute==1.and.istate_==jstate_.and.ilambda_==jlambda_) cycle
                      !
                      ! check if we at the right electronic states
                      if( istate/=istate_.or.jstate/=jstate_ ) cycle
                      !
                      ! We should also take into account that Lambda can change sign (only Lambda>0 is given in input)
                      ! In order to recover other combinations we apply the symmetry transformation
                      ! laboratory fixed inversion which is equivalent to the sigmav operation 
                      !                    (sigmav= 0 correspond to the unitary transformation)
                      do isigmav = 0,1
                        !
                        ! the permutation is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
                        ! avoid the double counting.
                        if( isigmav==1.and. abs( field%lambda ) + abs( field%lambdaj )==0 ) cycle
                 
                        ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                        ilambda_ = ilambda_*(-1)**isigmav
                        jlambda_ = jlambda_*(-1)**isigmav
                        !
                        ! proceed only if the quantum numbers of the field equal to the corresponding <i| and |j> quantum numbers:
                        if (ilambda_/=ilambda.or.jlambda_/=jlambda) cycle
                        !
                        !
                        ! double check
                        if (spini/=poten(istate)%spini.or.spinj/=poten(jstate)%spini) then
                         write(out,'("SR: reconstructed spini ",f8.1," or spinj ",f8.1," do not agree with stored values ", & 
                                    & f8.1,1x,f8.1)') spini,spinj,poten(istate)%spini,poten(jstate)%spini
                          stop 'SR: wrongly reconsrtucted spini or spinj'
                        endif
                        !
                        f_grid  = field%matelem(ivib,jvib)
                        !
                        ! <Lx> and <Ly> don't depend on Sigma
                        !
                        ! L*S part of the spin-rotation 
                        !
                        ! the selection rules are Delta Sigma = - Delta Lambda (Delta Spin = 0)
                        !
                        ! factor to switch between <Sigma+1|S+|Sigma> and <Sigma-1|S-|Sigma>:
                        f_s = real(ilambda-jlambda,rk)
                        !
                        ! the bra-component of Sigma (i.e. sigmaj):
                        sigmaj_ = sigmai+f_s
                        !
                        ! make sure that this sigmaj_ is consistent with the current ket-sigmaj
                        if (nint(2.0_rk*sigmaj_)==nint(2.0*sigmaj)) then
                          !
                          f_t = f_grid
                          !
                          ! the result of the symmetry transformation:
                          if (isigmav==1) then
                            !
                            ! we assume that
                            ! sigmav <Lamba|L+|Lambda'> => <-Lamba|L-|-Lambda'> == <Lamba|L+|Lambda'>(-1)^(Lamba+Lambda')
                            ! and <Lamba|L+|Lambda'> is an unique quantity given in the input
                            ! however we don't apply the sigmav transformation to sigma or omega
                            ! since we only need to know how <Lamba|L+/-|Lambda'> transforms in order to relate it to the
                            ! value given in input.
                            !
                            itau = 0 
                            !
                            if (ilambda_==0.and.poten(istate)%parity%pm==-1) itau = itau+1
                            if (jlambda_==0.and.poten(jstate)%parity%pm==-1) itau = itau+1
                            !
                            f_t = f_t*(-1.0_rk)**(itau)
                            !
                          endif
                           !
                          ! the matrix element <Sigmai| S+/- |Sigmai+/-1>
                          !
                          f_t = sqrt( (spini-f_s*sigmai)*(spini+f_s*sigmai+1.0_rk) )*f_t
                          !
                          !f_t = sqrt( spini*(spini+1.0_rk)-sigmai*(sigmai+f_s)  )*f_t
                          !
                          hmat(i,j) = hmat(i,j) - f_t*0.5_rk
                          !
                          ! print out the internal matrix at the first grid point
                          if (iverbose>=4.and.abs(hmat(i,j))>sqrt(small_)) then
                             write(printout_,'(i3,"-SR",2i3)') isr,ilevel,jlevel
                             printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                             write(printout_,'(g12.4)') hmat(i,j)/sc
                             printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                          endif
                          !
                        endif
                        !
                      enddo
                      !
                    enddo
                    !
                endif
                ! 
              enddo
              !
              ! J*S part (S-uncoupling)
              !
              ! The selection rules are Delta Spin=0, Delta Lambda = 0, Delta Sigma = +/- 1
              ! and (CHECK!!) istate==jstate (?)
              !
              if (abs(nint(sigmaj-sigmai))==1.and.ilambda==jlambda.and.nint(spini-spinj)==0.and.istate==jstate) then
                !
                !if (nint(sigmai-sigmaj)/=nint(omegai-omegaj)) cycle
                !
                f_s = sigmaj-sigmai
                !
                !f_t = sqrt( (spinj-f_s*sigmaj )*( spinj + f_s*sigmaj+1.0_rk ) )*&
                !      sqrt( (jval -f_s*omegaj )*( jval  + f_s*omegaj+1.0_rk ) )
                !
                f_t = sqrt( spini*(spini+1.0_rk)-sigmai*(sigmai+f_s) )*&  !== sqrt( ... sigmai*sigmaj)
                      sqrt( jval* (jval +1.0_rk)-omegai*(omegai+f_s) )    !== sqrt( ... omegai*omegaj)
                !
                if( job%zExclude_JS_coupling .eqv. .true.)  f_t = 0.0_rk ! do not include J.S coupling if flag is set
                !
                hmat(i,j) = hmat(i,j) - f_t*f_rot
                !
                !hmat(j,i) = hmat(i,j)
                !
                if ((nint(omegai-omegaj))/=nint(sigmai-sigmaj)) then
                  write(out,'("J*S: omegai-omegaj/=sigmai-sigmaj ",2f8.1,2x,2f8.1," for i,j=",2(i0,2x))') omegai,omegaj, &
                                                                                                         sigmai,sigmaj,i,j
                  stop 'J*S: omegai/=omegaj+/-1 '
                endif
                !
                ! print out the internal matrix at the first grid point
                if (iverbose>=4.and.abs(hmat(i,j))>sqrt(small_)) then
                   write(printout_,'("  J-S(",2i3,")=")') ilevel,jlevel
                   printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                   if (abs(hmat(i,j))>sqrt(small_)) then
                     write(printout_,'(F12.4, A)') hmat(i,j)/sc, " ;"
                     printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                   endif
                endif
                !
              endif
              !
              ! spin-orbit part:
              loop_iso : do iso =1,Nspinorbits
                !
                field => spinorbit(iso)
                !
                ! The selection rules are (Lefebvre-Brion and Field, Eq. (3.4.6)): 
                ! Delta J = 0 ; Delta Omega  = 0 ; g<-/->u; e<->f; Sigma+<->Sigma-;
                ! Delta S = 0 or Delta S = 1 ; Delta Lambda = Delta Sigma = 0 or Delta Lambda = - Delta Sigma = +/- 1
                !
                if (nint(omegai-omegaj)/=0.or.nint(spini-spinj)>1 ) cycle
                if ( ilambda==0.and.jlambda==0.and.poten(istate)%parity%pm==poten(jstate)%parity%pm ) cycle
                if ( poten(istate)%parity%gu/=0.and.poten(istate)%parity%gu/=poten(jstate)%parity%gu ) cycle
                !
                do ipermute  = 0,1
                  !
                  if (ipermute==0) then
                    !
                    istate_ = field%istate ; ilambda_we = field%lambda  ; sigmai_we = field%sigmai ; spini_ = field%spini
                    jstate_ = field%jstate ; jlambda_we = field%lambdaj ; sigmaj_we = field%sigmaj ; spinj_ = field%spinj
                    !
                  else  ! permute
                    !
                    jstate_ = field%istate ; jlambda_we = field%lambda  ; sigmaj_we = field%sigmai ; spinj_ = field%spini
                    istate_ = field%jstate ; ilambda_we = field%lambdaj ; sigmai_we = field%sigmaj ; spini_ = field%spinj
                    !
                  endif
                  ! proceed only if the spins of the field equal the corresponding <i| and |j> spins of the current matrix elements. 
                  ! otherwise skip it:
                  if ( nint(spini_-spini)/=0.or.nint(spinj_-spinj)/=0 ) cycle
                  !
                  ! however the permutation makes sense only when for non diagonal <State,Lambda,Spin|F|State',Lambda',Spin'>
                  ! otherwise it will cause a double counting:
                  !
                  if (ipermute==1.and.istate_==jstate_.and.ilambda_we==jlambda_we.and.nint(sigmai_we-sigmaj_we)==0.and. & 
                      nint(spini_-spinj_)==0) cycle
                  !
                  ! check if we at the right electronic states
                  if( istate/=istate_.or.jstate/=jstate_ ) cycle
                  !
                  ! We apply the Wigner-Eckart theorem to reconstruct all combinations of <Lamba Sigma |HSO|Lamba Sigma' > 
                  ! connected with the reference (input) <Lamba Sigma_r |HSO|Lamba Sigma_r' > by this theorem. 
                  ! Basically, we loop over Sigma (Sigma = -S..S).  The following 3j-symbol for the reference states will 
                  ! be conidered:
                  ! / Si      k  Sj     \    k  = 1
                  ! \ -Sigmai q  Sigmaj /    q  = Sigmai - Sigmaj
                  !
                  ! reference q from Wigner-Eckart
                  q_we = sigmai_we-sigmaj_we
                  !
                  ! We should consider also a permutation <State',Lambda',Spin'|F|State,Lambda,Spin> if this makes a change.
                  ! This could be imortant providing that we constrain the i,j indexes to be i<=j (or i>=j).
                  ! We also assume that the matrix elements are real!
                  !
                  ! First of all we can check if the input values are not unphysical and consistent with Wigner-Eckart:
                  ! the corresponding three_j should be non-zero:
                  three_j_ref = three_j(spini_, 1.0_rk, spinj_, -sigmai_we, q_we, sigmaj_we)
                  !
                  if (abs(three_j_ref)<small_) then 
                    !
                    write(out,"('The Spin-orbit field ',2i3,' is incorrect according to Wigner-Eckart, three_j = 0 ')") & 
                          field%istate,field%jstate
                    write(out,"('Check S_i, S_j, Sigma_i, Sigma_j =  ',4f9.2)") spini_,spinj_,sigmai_we,sigmaj_we
                    stop "The S_i, S_j, Sigma_i, Sigma_j are inconsistent"
                    !
                  end if 
                  !
                  ! Also check the that the SO is consistent with the selection rules for SO
                  !
                  if ( ilambda_we-jlambda_we+nint(sigmai_we-sigmaj_we)/=0.or.nint(spini_-spinj_)>1.or.&
                     ( ilambda_we==0.and.jlambda_we==0.and.poten(field%istate)%parity%pm==poten(field%jstate)%parity%pm ).or.&
                     ( (ilambda_we-jlambda_we)/=-nint(sigmai_we-sigmaj_we) ).or.&
                        abs(ilambda_we-jlambda_we)>1.or.abs(nint(sigmai_we-sigmaj_we))>1.or.&
                     ( poten(field%istate)%parity%gu/=0.and.poten(field%istate)%parity%gu/=poten(field%jstate)%parity%gu ) ) then
                     !
                     write(out,"('The quantum numbers of the spin-orbit field ',2i3,' are inconsistent" // &
                                     " with SO selection rules: ')") field%istate,field%jstate
                     write(out,"('Delta J = 0 ; Delta Omega  = 0 ; g<-/->u; e<-/->f; Sigma+<->Sigma-; " // &
                        "Delta S = 0 or Delta S = 1 ; Delta Lambda = Delta Sigma = 0 or Delta Lambda = - Delta Sigma = +/- 1')")
                     write(out,"('Check S_i, S_j, Sigma_i, Sigma_j, lambdai, lambdaj =  ',4f9.2,2i4)") &
                                                                        spini_,spinj_,sigmai_we,sigmaj_we,ilambda_we,jlambda_we
                     stop "The S_i, S_j, Sigma_i, Sigma_j lambdai, lambdaj are inconsistent with selection rules"
                     !
                  endif
                  !
                  do isigma2 = -nint(2.0*spini_),nint(2.0*spini_),2
                    !
                    ! Sigmas from Wigner-Eckart
                    sigmai_ = real(isigma2,rk)*0.5 
                    sigmaj_ = sigmai_ - q_we
                    !
                    ! three_j for current Sigmas
                    three_j_ = three_j(spini_, 1.0_rk, spinj_, -sigmai_, q_we, sigmaj_)
                    !
                    ! current value of the SO-matrix element from Wigner-Eckart
                    SO = (-1.0_rk)**(sigmai_-sigmai_we)*three_j_/three_j_ref*field%matelem(ivib,jvib)
                    !
                    ! We should also take into account that Lambda and Sigma can change sign
                    ! since in the input we give only a unique combination of matrix elements, for example
                    ! < 0 0 |  1  1 > is given, but < 0 0 | -1 -1 > is not, assuming that the program will generate the missing
                    ! combinations.
                    !
                    ! In order to recover other combinations we apply the symmetry transformation
                    ! laboratory fixed inversion which is equivalent to the sigmav operation 
                    !                    (sigmav= 0 correspond to the unitary transformation)
                    do isigmav = 0,1
                      !
                      ! sigmav is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
                      ! avoid the double counting.
                      if( isigmav==1.and.nint( abs( 2.0*sigmai_ )+ abs( 2.0*sigmaj_ ) )+abs( ilambda_we )+abs( jlambda_we )==0 ) &
                                                                                                                             cycle
                      !
                      ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                      ilambda_ = ilambda_we*(-1)**isigmav
                      jlambda_ = jlambda_we*(-1)**isigmav
                      sigmai_ = sigmai_*(-1.0_rk)**isigmav
                      sigmaj_ = sigmaj_*(-1.0_rk)**isigmav
                      !
                      omegai_ = sigmai_+real(ilambda_)
                      omegaj_ = sigmaj_+real(jlambda_)
                      !
                      ! Check So selection rules
                      if ( ( ilambda_-jlambda_)/=-nint(sigmai_-sigmaj_).or.abs(sigmai_-sigmaj_)>1.or.omegai_/=omegaj_ ) cycle
                      !
                      ! proceed only if the quantum numbers of the field equal
                      ! to the corresponding <i| and |j> quantum numbers of the basis set. otherwise skip it:
                      if ( nint(sigmai_-sigmai)/=0.or.nint(sigmaj_-sigmaj)/=0.or.ilambda_/=ilambda.or.jlambda_/=jlambda ) cycle
                      !
                      f_t = SO*sc
                      !
                      ! the result of the symmetry transformtaion applied to the <Lambda,Sigma|HSO|Lambda',Sigma'> only
                      if (isigmav==1) then
                        !
                        ! still not everything is clear here: CHECK!
                        !
                        itau = -ilambda_-jlambda_ +nint(spini_-sigmai_)+nint(spinj_-sigmaj_) !+nint(jval-omegai)+(jval-omegaj)
                        !
                        !itau = nint(spini_-sigmai_)+nint(spinj_-sigmaj_) ! +nint(jval-omegai)+(jval-omegaj)
                        !
                        !itau = 0
                        !
                        if (ilambda_==0.and.poten(istate)%parity%pm==-1) itau = itau+1
                        if (jlambda_==0.and.poten(jstate)%parity%pm==-1) itau = itau+1
                        !
                        f_t = f_t*(-1.0_rk)**(itau)
                        !
                      endif
                      !
                      ! double check
                      if ( nint(omegai-omegai_)/=0 .or. nint(omegaj-omegaj_)/=0 ) then
                        write(out,'(A,f8.1," or omegaj ",f8.1," do not agree with stored values ",f8.1,1x,f8.1)') &
                                   "SO: reconsrtucted omegai", omegai_,omegaj_,omegai,omegaj
                        stop 'SO: wrongly reconsrtucted omegai or omegaj'
                      endif
                      !
                      ! we might end up in eilther parts of the matrix (upper or lower),
                      ! so it is safer to be general here and
                      ! don't restrict to lower part as we have done above
                      !
                      hmat(i,j) = hmat(i,j) + f_t
                      !
                      !hmat(j,i) = hmat(i,j)
                      !
                      ! print out the internal matrix at the first grid point
                      if (iverbose>=4.and.abs(hmat(i,j))>small_) then
                          !
                          write(printout_,'(i3,"-SO",2i3)') iso,ilevel,jlevel
                          printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                          write(printout_,'(g12.4)') f_t/sc
                          printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                         !
                      endif
                      !
                      cycle loop_iso
                      !
                    enddo
                  enddo
                enddo
              enddo  loop_iso
              !
              !
              ! L*S and J*L parts
              !
              loop_ilxly : do ilxly =1,Nlxly
                !
                field => lxly(ilxly)
                !
                ! Also check that L+ is consistent with the selection rules
                !
                if ( field%istate==field%jstate .or.abs(field%lambda-field%lambdaj)/=1 ) then
                   !
                   write(out,"('The quantum numbers of the L+/Lx field ',2i3,' are inconsistent" // &
                                   " with L+selection rules: ')") field%istate,field%jstate
                   write(out,"('Delta Lamda = +/-1')")
                   stop "Lx/L+ input is inconsistent with selection rules"
                   !
                endif
                !
                ! the field entry in the input gives only one combination of the quantum numbers for
                ! the matrix element <State,Lambda,Spin|F|State',Lambda',Spin'>
                ! LxLy  4 6 ;  lambda  0 1 ; spin   1.0 1.0
                ! we should consider also a permutation <State',Lambda',Spin'|F|State,Lambda,Spin> if this makes a change.
                ! This could be imortant providing that we constrain the i,j indexes to be i<=j (or i>=j)
                ! We also assume that the matrix elements are real!
                !
                do ipermute  = 0,1
                  !
                  if (ipermute==0) then
                    !
                    istate_ = field%istate ; ilambda_ = field%lambda  ; spini_ = field%spini
                    jstate_ = field%jstate ; jlambda_ = field%lambdaj ; spinj_ = field%spinj
                    !
                  else  ! permute
                    !
                    jstate_ = field%istate ; jlambda_ = field%lambda  ; spinj_ = field%spini
                    istate_ = field%jstate ; ilambda_ = field%lambdaj ; spini_ = field%spinj
                    !
                  endif
                  !
                  ! however the permutation makes sense only when for non diagonal <State,Lambda,Spin|F|State',Lambda',Spin'>
                  ! otherwise it will cause a double counting:
                  !
                  if (ipermute==1.and.istate_==jstate_.and.ilambda_==jlambda_.and.nint(spini_-spinj_)==0) cycle
                  !
                  ! check if we at the right electronic states
                  if( istate/=istate_.or.jstate/=jstate_ ) cycle
                  !
                  ! We should also take into account that Lambda can change sign (only Lambda>0 is given in input)
                  ! In order to recover other combinations we apply the symmetry transformation
                  ! laboratory fixed inversion which is equivalent to the sigmav operation 
                  !                    (sigmav= 0 correspond to the unitary transformation)
                  do isigmav = 0,1
                    !
                    ! the permutation is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
                    ! avoid the double counting.
                    if( isigmav==1.and. abs( field%lambda ) + abs( field%lambdaj )==0 ) cycle
       
                    ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                    ilambda_ = ilambda_*(-1)**isigmav
                    jlambda_ = jlambda_*(-1)**isigmav
                    !
                    ! proceed only if the quantum numbers of the field equal to the corresponding <i| and |j> quantum numbers:
                    if (istate/=istate_.or.jstate_/=jstate.or.ilambda_/=ilambda.or.jlambda_/=jlambda) cycle
                    !
                    ! check the selection rule Delta Lambda = +/1
                    if (abs(ilambda-jlambda)/=1) cycle
                    !
                    ! double check
                    if (spini/=poten(istate)%spini.or.spinj/=poten(jstate)%spini) then
                     write(out,'("LJ: reconstructed spini ",f8.1," or spinj ",f8.1," do not agree with stored values ", & 
                                & f8.1,1x,f8.1)') spini,spinj,poten(istate)%spini,poten(jstate)%spini
                      stop 'LJ: wrongly reconsrtucted spini or spinj'
                    endif
                    !
                    f_grid  = field%matelem(ivib,jvib)
                    !
                    ! <Lx> and <Ly> don't depend on Sigma
                    !
                    ! L*S part (spin-electronic coupling)
                    !
                    ! the selection rules are Delta Sigma = - Delta Lambda (Delta Spin = 0)
                    !
                    ! factor to switch between <Sigma+1|S+|Sigma> and <Sigma-1|S-|Sigma>:
                    f_s = real(ilambda-jlambda,rk)
                    !
                    ! the bra-component of Sigma (i.e. sigmaj):
                    sigmaj_ = sigmai+f_s
                    !
                    ! make sure that this sigmaj_ is consistent with the current ket-sigmaj
                    if (nint(2.0_rk*sigmaj_)==nint(2.0*sigmaj)) then
                      !
                      f_t = f_grid
                      !
                      ! the result of the symmetry transformation:
                      if (isigmav==1) then
                        !
                        ! we assume that
                        ! sigmav <Lamba|L+|Lambda'> => <-Lamba|L-|-Lambda'> == <Lamba|L+|Lambda'>(-1)^(Lamba+Lambda')
                        ! and <Lamba|L+|Lambda'> is an unique quantity given in the input
                        ! however we don't apply the sigmav transformation to sigma or omega
                        ! since we only need to know how <Lamba|L+/-|Lambda'> transforms in order to relate it to the
                        ! value given in input.
                        !
                        !itau = ilambda-jlambda+nint(spini-sigmai)+nint(spinj-sigmaj)!-nint(omegai+omegaj)
                        !
                        ! we try to remove also lambda from the sigmav transformation!!
                        !
                        itau = 0 !!!! ilambda-jlambda
                        !
                        !itau = ilambda-jlambda
                        !
                        if (ilambda_==0.and.poten(istate)%parity%pm==-1) itau = itau+1
                        if (jlambda_==0.and.poten(jstate)%parity%pm==-1) itau = itau+1
                        !
                        f_t = f_t*(-1.0_rk)**(itau)
                        !
                        ! change the sign of <|L+|> to get sigmav*<|L+|>=-<|L-|> if necessary
                        !! f_t = sign(1.0_rk,f_s)*f_t
                        !
                      endif
                       !
                      ! the matrix element <Sigmai| S+/- |Sigmai+/-1>
                      !
                      f_t = sqrt( (spini-f_s*sigmai)*(spini+f_s*sigmai+1.0_rk) )*f_t
                      !
                      !f_t = sqrt( spini*(spini+1.0_rk)-sigmai*(sigmai+f_s)  )*f_t
                      !
                      hmat(i,j) = hmat(i,j) + f_t
                      !hmat(j,i) = hmat(i,j)
                      !
                      ! print out the internal matrix at the first grid point
                      if (iverbose>=4.and.abs(hmat(i,j))>sqrt(small_)) then
                         write(printout_,'(i3,"-LS",2i3)') ilxly,ilevel,jlevel
                         printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                         write(printout_,'(g12.4)') hmat(i,j)/sc
                         printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                      endif
                      !
                    endif
                    !
                    ! L*J part (L-uncoupling)
                    !
                    ! The selection rule is simple: Delta Sigma = 0, Delta Spin = 0,
                    ! i.e. bra and ket sigma are equal:
                    !
                    if (nint(2.0_rk*sigmai)==nint(2.0*sigmaj)) then
                      !
                      ! however omega should change by 1 (via J+/-) exactly as lambda changes (with L-/+):
                      ! (f_l will be needed to switch between J+ and J-)
                      f_l = real(jlambda-ilambda,rk)
                      !
                      ! we should obtain  omegaj = omega+f_l
                      ! double check
                      if ( nint( 2.0_rk*omegaj )/=nint( 2.0_rk*(omegai+f_l) ) ) then
                         write(out,'("L*J omegaj ",f8.1," does agree with assumed ",f8.1," value omegai+/-1")') omegaj,omegai+f_l
                         stop 'wrongly reconsrtucted omegaj'
                      endif
                      !
                      f_t = f_grid
                      !
                      ! the result of the symmetry transformation sigmav:
                      !
                      if (isigmav==1) then
                        !
                        ! we assume that
                        ! sigmav <Lamba|L+|Lambda'> => <-Lamba|L-|-Lambda'> == <Lamba|L+|Lambda'>(-1)^(Lamba+Lambda')
                        ! and <Lamba|L+|Lambda'> is an unique quantity given in the input
                        ! (see alos above)
                        !
                        !itau = ilambda-jlambda+nint(spini-sigmai)+nint(spinj-sigmaj)-nint(omegai+omegaj)
                        !
                        ! we try now removing lambda from sigmav transormation!!!
                        !
                        itau = 0 !! ilambda-jlambda
                        !
                        !itau = ilambda-jlambda
                        !
                        if (ilambda==0.and.poten(istate)%parity%pm==-1) itau = itau+1
                        if (jlambda==0.and.poten(jstate)%parity%pm==-1) itau = itau+1
                        !
                        f_t = f_t*(-1.0_rk)**(itau)
                        !
                      endif
                      !
                      ! For <Omega|<Lambda|L+ J- |Lambda+1>|Omega+1> f_l = 1
                      ! For <Omega|<Lambda|L- J+ |Lambda-1>|Omega-1> f_l =-1
                      !
                      f_t = sqrt( (jval-f_l*omegai)*(jval+f_l*omegai+1.0_rk) )*f_t
                      !
                      !f_t = sqrt( jval*(jval+1.0_rk)-omegai*(omega+f_l) )*f_t
                      !
                      hmat(i,j) = hmat(i,j) - f_t
                      !hmat(j,i) = hmat(i,j)
                      !
                      ! print out the internal matrix at the first grid point
                      if (iverbose>=4.and.abs(hmat(i,j))>small_) then
                         write(printout_,'(i3,"-LJ",2i3)') ilxly,ilevel,jlevel
                         printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                         write(printout_,'(g12.4)') hmat(i,j)/sc
                         printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                      endif
                      !
                    endif
                    !
                  enddo
                  !
                enddo
                !
              enddo loop_ilxly
              !
              !
              ! Non-diagonal lambda-opq doubling
              !
              do ild = 1,Nlambdaopq
                !
                field => lambdaopq(ild)
                !
                ! 1. <Sigma,Omega,Lambda|Lambda-O|Sigma+/-2,Omega,-Lambda>
                if (lambdaopq(ild)%istate==istate.and.lambdaopq(ild)%jstate==jstate.and.istate==jstate.and.&
                    abs(ilambda)==1.and.(ilambda-jlambda)==nint(sigmaj-sigmai).and.abs(nint(sigmaj-sigmai))==2 &
                                   .and.(ilambda==-jlambda).and.nint(spini-spinj)==0.and.nint(omegai-omegaj)==0) then
                   !
                   f_s2 = sigmai-sigmaj
                   f_s1 = sign(1.0_rk,f_s2)
                   !
                   f_t = sqrt( spini*(spini+1.0_rk)-(sigmaj     )*(sigmaj+f_s1) )*&
                         sqrt( spini*(spini+1.0_rk)-(sigmaj+f_s1)*(sigmaj+f_s2) )
                   !
                   f_lo = field%matelem(ivib,jvib)*f_t*sc
                   !
                   hmat(i,j) = hmat(i,j) + f_lo*0.5_rk
                   !
                   ! print out the internal matrix at the first grid point
                   if (iverbose>=4.and.abs(hmat(i,j))>sqrt(small_)) then
                      write(printout_,'("    LO",2i3)') ilevel,jlevel
                      printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                      if (abs(hmat(i,j))>sqrt(small_)) then
                        write(printout_,'(g12.4)') hmat(i,j)/sc
                        printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                      endif
                   endif
                   !
                endif
                ! 
              enddo
              !
              ! Non-diagonal lambda-p2q doubling
              !
              do ild = 1,Nlambdap2q
                !
                field => lambdap2q(ild)
                !
                ! 1. <Sigma,Omega,Lambda|Lambda-p|Sigma+/-1,Omega-/+1,Lambda>
                if (lambdap2q(ild)%istate==istate.and.lambdap2q(ild)%jstate==jstate.and.istate==jstate.and.&
                   abs(ilambda)==1.and.abs(nint(sigmaj-sigmai))==1.and.(ilambda==-jlambda).and.nint(spini-spinj)==0.and.&
                   abs(nint(omegai-omegaj))==1.and.nint(sigmaj-sigmai)==nint(omegai-omegaj)) then
                   !
                   f_s = sigmai-sigmaj
                   !
                   f_t = sqrt( spini*(spini+1.0_rk)-sigmaj*(sigmaj+f_s) )*&
                         sqrt( jval* (jval +1.0_rk)-omegaj*(omegaj-f_s) )  
                   !
                   f_lo = field%matelem(ivib,jvib)*f_t*sc
                   !
                   hmat(i,j) = hmat(i,j) - f_lo*0.5_rk
                   !
                   ! print out the internal matrix at the first grid point
                   if (iverbose>=4.and.abs(hmat(i,j))>sqrt(small_)) then
                      write(printout_,'("    LP",2i3)') ilevel,jlevel
                      printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                      if (abs(hmat(i,j))>sqrt(small_)) then
                        write(printout_,'(g12.4)') hmat(i,j)/sc
                        printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                      endif
                   endif
                   !
                endif
                ! 
              enddo
              !
              ! Non-diagonal lambda-q doubling
              !
              do ild = 1,Nlambdaq
                !
                field => lambdaq(ild)
                !
                ! 1. <Sigma,Omega,Lambda|Lambda-O|Sigma+/-2,Omega,-Lambda>
                if (lambdaq(ild)%istate==istate.and.lambdaq(ild)%jstate==jstate.and.istate==jstate.and.&
                    abs(ilambda)==1.and.(ilambda-jlambda)==nint(omegai-omegaj).and.abs(nint(sigmaj-sigmai))==0.and.&
                       (ilambda==-jlambda).and.nint(spini-spinj)==0.and.nint(omegai-omegaj)==2) then
                   !
                   f_o2 = omegaj-omegai
                   f_o1 = sign(1.0_rk,f_o2)
                   !
                   f_t = sqrt( jval*(jval+1.0_rk)-(omegaj     )*(omegaj-f_o1) )*&
                         sqrt( jval*(jval+1.0_rk)-(omegaj-f_o1)*(omegaj-f_o2) )
                   !
                   f_lo = field%matelem(ivib,jvib)*f_t*sc
                   !
                   hmat(i,j) = hmat(i,j) + f_lo*0.5_rk
                   !
                   ! print out the internal matrix at the first grid point
                   if (iverbose>=4.and.abs(hmat(i,j))>sqrt(small_)) then
                      write(printout_,'("    LO",2i3)') ilevel,jlevel
                      printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                      if (abs(hmat(i,j))>sqrt(small_)) then
                        write(printout_,'(g12.4)') hmat(i,j)/sc
                        printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                      endif
                   endif
                   !
                endif
                ! 
              enddo
              !
           enddo  ! j
         enddo  ! i
         !omp end parallel do
         !
         if (iverbose>=3) write(out,'("...done!")')
         !
         if (iverbose>=5) then
            ! print out the structure of the submatrix
            !
            write(out,'(/"Non-zero matrix elements of the coupled Sigma-Lambda matrix:")')
            write(out, '(A, ES10.2,A)') 'Threshold for printing is ', small_, ' cm^-1'
            write(out,'(A)') 'RV == Rotational-vibrational'
            write(out,'(A)') 'SO == Spin-Orbit interaction'
            write(out,'(A)') 'SS == Spin-Spin interaction'
            write(out,'(A)') 'JS == J.S interaction (aka S-uncoupling)'
            write(out,'(A)') 'LS == L.J interaction (aka L-uncoupling)'
            write(out,'(A)') 'LS == L.S interaction (spin-electronic)'
            !
            do ilevel = 1,Nlambdasigmas
              write(out,'(a)') trim( printout(ilevel) )
            enddo
            !
            write(out,'(" "/)')
            !
         endif
         !
         if (iverbose>=4) call TimerStop('Construct the hamiltonian')
         !
         ! Transformation to the symmetrized basis set
         !
         ! |v,Lambda,Sigma,J,Omega,tau> = 1/sqrt(2) [ |v,Lambda,Sigma,J,Omega>+(-1)^tau |v,-Lambda,-Sigma,J,-Omega> ]
         !
         allocate(iswap(Ntotal),vec(Ntotal),Nirr(Ntotal,2),ilevel2i(Ntotal,2),tau(Ntotal),ilevel2isym(Ntotal,2),stat=alloc)
         call ArrayStart('iswap-vec',alloc,size(iswap),kind(iswap))
         call ArrayStart('iswap-vec',alloc,size(vec),kind(vec))
         call ArrayStart('Nirr',alloc,size(Nirr),kind(Nirr))
         !
         iswap = 0
         Nsym = 0
         Nlevels = 0
         ilevel2i = 0
         ilevel2isym = 0
         Nirr = 0
         !
         !omp parallel do private(istate,sigmai,ilambda,spini,omegai,ibib,j,jstate,sigmaj,jlambda,omegaj,spinj,jvib) & 
         !                        shared(iswap,vec) schedule(guided)
         do i = 1,Ntotal
           !
           istate = icontr(i)%istate
           sigmai = icontr(i)%sigma
           ilambda = icontr(i)%ilambda
           spini = icontr(i)%spin
           omegai = real(ilambda,rk)+sigmai
           ivib    = icontr(i)%v
           !
           if (iswap(i)/=0) cycle
           !
           if (ilambda==0.and.nint(sigmai)==0) then
             !
             Nlevels = Nlevels + 1
       
             itau = nint(spini+Jval)
             if (poten(istate)%parity%pm==-1) itau = itau+1
             !
             itau = abs(itau)
             !
             tau(Nlevels) = (-1.0_rk)**itau
             !
             if (mod(itau,2)==1) then
               Nsym(2) = Nsym(2) + 1
               ilevel2isym(Nlevels,2) = nsym(2)
               ilevel2i(Nlevels,2) = i
               Nirr(Nlevels,2) = 1
             else
               Nsym(1) = Nsym(1) + 1
               ilevel2isym(Nlevels,1) = nsym(1)
               ilevel2i(Nlevels,1) = i
               Nirr(Nlevels,1) = 1
             endif
             !
             iswap(i) = i
             !ilevel2i(Nlevels,1) = i
             !
           else
             !
             do  j = 1,Ntotal
               !
               jstate  = icontr(j)%istate
               sigmaj  = icontr(j)%sigma
               jlambda = icontr(j)%ilambda
               omegaj  = real(jlambda,rk)+sigmaj
               spinj   = icontr(j)%spin
               jvib    = icontr(j)%v
               !
               if (ilambda==-jlambda.and.nint(2.0*sigmai)==-nint(2.0*sigmaj).and. &
                   istate==jstate.and.ivib==jvib.and.nint(2.0*spini)==nint(2.0*spinj)) then
                 !
                 Nsym(:) = Nsym(:) + 1
                 Nlevels = Nlevels + 1
                 Nirr(Nlevels,:) = 1
                 !
                 if (ilambda>jlambda.or.sigmaj<sigmai) then
                   !
                   ilevel2i(Nlevels,1) = i
                   ilevel2i(Nlevels,2) = j
                   !
                 else
                   !
                   ilevel2i(Nlevels,1) = j
                   ilevel2i(Nlevels,2) = i
                   !
                 endif
                 !
                 ilevel2isym(Nlevels,1:2) = nsym(1:2)
                 !
                 itau = -ilambda+nint(spini-sigmai)+nint(Jval-omegai)
                 !
                 if (ilambda==0.and.poten(istate)%parity%pm==-1) itau = itau+1
                 !
                 tau(Nlevels) = (-1.0_rk)**itau
                 !
                 iswap(i) = j*(-1)**itau
                 iswap(j) = i*(-1)**itau
                 !
                 exit
                 !
               endif
               !
             enddo
             !
           endif
           !
         enddo
         !omp end parallel do
         !
         ! Nlevels is the number of states disregarding the degeneracy 
         ! Nroots is the total number of roots including the degenerate states 
         !
         allocate(transform(1)%matrix(max(1,Nsym(1)),max(1,Nsym(1))),stat=alloc)
         allocate(transform(2)%matrix(max(1,Nsym(2)),max(1,Nsym(2))),stat=alloc)
         allocate(transform(1)%irec( max( 1,Nsym(1) ) ),stat=alloc)
         allocate(transform(2)%irec( max( 1,Nsym(2) ) ),stat=alloc)
         !
         call ArrayStart('transform',alloc,size(transform(1)%matrix),kind(transform(1)%matrix))
         call ArrayStart('transform',alloc,size(transform(2)%matrix),kind(transform(2)%matrix))
         call ArrayStart('transform',alloc,size(transform(1)%irec),kind(transform(1)%irec))
         call ArrayStart('transform',alloc,size(transform(2)%irec),kind(transform(2)%irec))
         !
         allocate(Utransform(Nlevels,2,2),stat=alloc)
         call ArrayStart('Utransform',alloc,size(Utransform),kind(Utransform))
         !
         ! Building the transformation to the symmetrized representaion 
         !
         do ilevel = 1,Nlevels
           !
           veci = 0
           !
           if (any(Nirr(ilevel,:)==0)) then
             !
             do irrep = 1,sym%NrepresCs
               do itau = 1,Nirr(ilevel,irrep)
                 !
                 i = ilevel2i(ilevel,irrep)
                 istate = icontr(i)%istate
                 veci(irrep,irrep) = 1.0_rk
                 isym = ilevel2isym(ilevel,irrep)
                 transform(irrep)%irec(isym) = ilevel
                 !
               enddo
             enddo
             !
             !
             !if (tau(ilevel)<0) then
             !  veci(1,2) = 1.0_rk
             !  isym = ilevel2isym(ilevel,2)
             !  transform(2)%irec(isym) = ilevel
             !else
             !  veci(1,1) = 1.0_rk
             !  isym = ilevel2isym(ilevel,1)
             !  transform(1)%irec(isym) = ilevel
             !endif
             !
           else
             !
             veci(1,1) = sqrt(0.5_rk)
             veci(2,1) = sqrt(0.5_rk)*tau(ilevel)
             veci(1,2) = sqrt(0.5_rk)
             veci(2,2) =-sqrt(0.5_rk)*tau(ilevel)
             !
             do irrep = 1,sym%NrepresCs
                isym = ilevel2isym(ilevel,irrep)
                transform(irrep)%irec(isym) = ilevel
             enddo
             !
             !do itau=1,Nirr(ilevel)
             !  !
             !  isym = ilevel2isym(ilevel,itau)
             !  transform(itau)%irec(isym) = ilevel ! ilevel2i(ilevel,itau)
             !  !
             !enddo
             !
           endif
           !
           Utransform(ilevel,:,:) = veci(:,:)
           !
           do jlevel = 1,Nlevels
             !
             vecj = 0
             !
             !if (Nirr(jlevel)==1) then
             !  !
             !  j = ilevel2i(jlevel,1)
             !  if (j==0) j = ilevel2i(jlevel,2)
             !  !
             !  jstate = icontr(j)%istate
             !  !
             !  if (tau(jlevel)<0) then
             !    vecj(1,2) = 1.0_rk
             !    jsym = ilevel2isym(jlevel,2)
             !    !transform(2)%irec(jsym) = i
             !  else
             !    vecj(1,1) = 1.0_rk
             !    jsym = ilevel2isym(jlevel,1)
             !    !transform(1)%irec(jsym) = i
             !  endif
             !
             if (any(Nirr(jlevel,:)==0)) then
                !
                do jrrep = 1,sym%NrepresCs
                  do jtau = 1,Nirr(jlevel,jrrep)
                    vecj(jrrep,jrrep) = 1.0_rk
                  enddo
                enddo
                !
             else
                !
                vecj(1,1) = sqrt(0.5_rk)
                vecj(2,1) = sqrt(0.5_rk)*tau(jlevel)
                vecj(1,2) = sqrt(0.5_rk)
                vecj(2,2) =-sqrt(0.5_rk)*tau(jlevel)
                !
             endif
             !
             pmat = 0
             !
             do isym = 1,2
                do itau = 1,Nirr(ilevel,isym)
                   i = ilevel2i(ilevel,isym)
                   do jsym = 1,2
                      do jtau = 1,Nirr(jlevel,jsym)
                         j = ilevel2i(jlevel,jsym)
                         !
                         if (i<=j) then
                           pmat(isym,jsym) = hmat(i,j)
                         else
                           pmat(isym,jsym) = hmat(j,i)
                         endif
                         !
                      enddo
                   enddo
                enddo
             enddo
       
             !do isym=1,Nirr(ilevel)
             !  i = ilevel2i(ilevel,isym)
             !  do jsym=1,Nirr(jlevel)
             !     j = ilevel2i(jlevel,jsym)
             !     !
             !     if (i<=j) then
             !       pmat(isym,jsym) = hmat(i,j)
             !     else
             !       pmat(isym,jsym) = hmat(j,i)
             !     endif
             !     !
             !  enddo
             !enddo
             !
             smat = matmul(transpose(veci),matmul(pmat,vecj))
             !
             !smat = matmul((veci),matmul(pmat,transpose(vecj)))
             !
             do irrep = 1,sym%NrepresCs
                do itau = 1,Nirr(ilevel,irrep)
                   !i = ilevel2i(ilevel,isym)
                   !
                   isym = ilevel2isym(ilevel,irrep)
                   !
                   do jrrep = 1,sym%NrepresCs
                      do jtau = 1,Nirr(jlevel,jrrep)
                         !j = ilevel2i(jlevel,jsym)
                         jsym = ilevel2isym(jlevel,jrrep)
                         !
                         if (irrep==jrrep) then
                           !
                           transform(irrep)%matrix(isym,jsym) = smat(irrep,irrep)
                           !
                         else
                           !
                           if (abs(smat(irrep,jrrep))>sqrt(small_)) then
                             !
                             i = ilevel2i(ilevel,itau)
                             j = ilevel2i(jlevel,jtau)
                             !
                             istate = icontr(i)%istate
                             sigmai = icontr(i)%sigma
                             ilambda = icontr(i)%ilambda
                             spini = icontr(i)%spin
                             omegai = real(ilambda,rk)+sigmai
                             ivib    = icontr(i)%v
                             !
                             jstate  = icontr(j)%istate
                             sigmaj  = icontr(j)%sigma
                             jlambda = icontr(j)%ilambda
                             omegaj  = real(jlambda,rk)+sigmaj
                             spinj   = icontr(j)%spin
                             jvib    = icontr(j)%v
                             !
                             write(out,'(/"Problem with symmetry: The non-diagonal matrix element is not zero:")')
                             write(out,'(/"i,j = ",2i8," irrep,jrrep = ",2i8," isym,jsym = ",2i8," ilevel,jlevel = ", &
                                        & 2i3," , matelem =  ",g16.9," with zero = ",g16.9)') &
                                        i,j,irrep,jrrep,isym,jsym,ilevel,jlevel,smat(itau,jtau),sqrt(small_)
                             write(out,'(/"<State   v  lambda spin   sigma  omega |H(sym)| State   v  lambda spin' //&
                              '   sigma  omega>")')
                             write(out,'("<",i3,2x,2i4,3f8.1," |H(sym)| ",i3,2x,2i4,3f8.1,"> /= 0")') &
                                         istate,ivib,ilambda,spini,sigmai,omegai,jstate,jvib,jlambda,spinj,sigmaj,omegaj
                             write(out,'("<",a10,"|H(sym)|",a10,"> /= 0")') trim(poten(istate)%name),trim(poten(jstate)%name)
                             !
                             stop 'Problem with symmetry: The non-diagonal matrix element is not zero'
                             !
                           endif
                           !
                         endif
                         !
                      enddo
                   enddo
                enddo
             enddo
       
           enddo
           !
         enddo
         !
       case default
         !
         write (out,"('error: illegal contraction',a)") trim(job%contraction)
         stop 'error - illegal CONTRACTION'
         !
       end select
       !
       ! Now we diagonalize the two matrices contructed one by one 
       !
       if (iverbose>=2) write(out,'(/"Eigenvalues for J = ",f8.1)') jval
       !
       ! the loop over the two parities
       do irrep = 1,sym%NrepresCs
          !
          nener_total = 0
          !
          Nsym_ = Nsym(irrep)
          !
          if (Nsym_<1) cycle
          !
          if (iverbose>=3) write(out,'(/"       J      i        Energy/cm  State   v  lambda spin   sigma   omega  parity")')
          !
          if (iverbose>=4) call TimerStart('Diagonalization')
          !
          ! Prepare the Hamiltonian matrix in the symmetrized representaion
          !
          allocate(eigenval(Nsym_),hsym(Nsym_,Nsym_),stat=alloc)
          call ArrayStart('eigenval',alloc,size(eigenval),kind(eigenval))
          call ArrayStart('hsym',alloc,size(hsym),kind(eigenval))
          !
          hsym = transform(irrep)%matrix
          !
          ! Diagonalization of the hamiltonian matrix
          !
          if (iverbose>=6) write(out,'(/"Diagonalization of the hamiltonian matrix")')
          !
          select case (trim(job%diagonalizer))
            !
          case('SYEV')
            !
            call lapack_syev(hsym,eigenval)
            !
            ! we need only these many roots
            Nroots = min(job%nroots(1),Nsym_)
            !
            ! or as many as below job%upper_ener if required by the input
            if (job%upper_ener<1e8) then
              nroots = maxloc(eigenval(:)-eigenval(1),dim=1,mask=eigenval(:).le.job%upper_ener*sc)
            endif
            !
          case('SYEVR')
            !
            ! some diagonalizers needs the following parameters to be defined
            !
            zpe = job%ZPE
            !
            jobz = 'V'
            vrange(1) = -0.0_rk ; vrange(2) = (job%upper_ener+zpe)*sc
            !
            if (.not.job%zShiftPECsToZero) vrange(1) = -safe_max
            !
            irange(1) = 1 ; irange(2) = min(job%nroots(1),Ntotal)
            nroots = job%nroots(1)
            !
            rng = 'A'
            !
            if (irange(2)==Nsym_) then
               rng = 'A'
            elseif (irange(2)<Nsym_) then
               rng = 'I'
            elseif (job%upper_ener<1e15.and.job%upper_ener>small_) then
               rng = 'V'
            endif
            !
            call lapack_syevr(hsym,eigenval,rng=rng,jobz=jobz,iroots=nroots,vrange=real(vrange,kind=8),irange=irange)
            !
            !
          case default
            !
            print "('Unrecognized diagonalizer ',a)", trim(job%diagonalizer)
            stop  'Unrecognized diagonalizer'
            !
          end select
          !
          if (iverbose>=6) write(out,'("...done!")')
          !
          if (iverbose>=4) call TimerStop('Diagonalization')
          !
          ! The ZPE value can be obtained only for the first J in the J_list tau=1.
          ! Otherwise the value from input must be used.
          if (irot==1.and.job%shift_to_zpe.and.iverbose/=0) then
            !
            if (irrep==1) then
              job%ZPE = eigenval(1)/sc
            endif
            !
            if (action%intensity) intensity%ZPE = job%ZPE
            !
          endif
          !
          eigenval(:) = eigenval(:)/sc
          !
          if (iverbose>=4) call TimerStart('Assignment')
          !
          ! Assign the eigevalues with quantum numbers and print them out
          !
          allocate(QNs(Nroots),stat=alloc)
          call ArrayStart('QNs',alloc,size(QNs),kind(QNs))
          !
          allocate(vib_count(Nroots),stat=alloc)
          call ArrayStart('vib_count',alloc,size(vib_count),kind(vib_count))
          vib_count = 0
          !
          !omp parallel do private(i,mterm,f_t,plusminus) shared(maxTerm) schedule(dynamic)
          do i=1,Nroots
            !
            ! to get the assignement we find the term with the largest contribution
            !
            imaxcontr = maxloc(hsym(:,i)**2,dim=1,mask=hsym(:,i)**2.ge.small_)
            !
            ! check if this contribution has already been used and take the second largest coefficient in this case
            !
            maxcontr = hsym(imaxcontr,i)**2
            !
            i0 = 0 
            !
            loop_check : do
              !
              loop_ilevel : do ilevel = 1,i-1
               !
               if (imaxcontr==QNs(ilevel)) then
                 !
                 j = maxloc(hsym(:,i)**2,dim=1,mask=hsym(:,i)**2.ge.small_.and.hsym(:,i)**2.lt.maxcontr)
                 imaxcontr = j
                 maxcontr = hsym(imaxcontr,i)**2
                 i0 = i0 + 1 
                 !
                 if (i0>10) then 
                   imaxcontr = maxloc(hsym(:,i)**2,dim=1,mask=hsym(:,i)**2.ge.small_)
                   exit loop_check
                 endif 
                 !
                 cycle loop_check 
                 !
               endif
               !
              enddo loop_ilevel
              exit loop_check
            enddo loop_check
            !
            QNs(i) = imaxcontr
            !
            mlevel = transform(irrep)%irec(imaxcontr)
            mterm = ilevel2i(mlevel,irrep)
            !
            istate = icontr(mterm)%istate
            sigma = icontr(mterm)%sigma
            imulti = icontr(mterm)%imulti
            ilambda = icontr(mterm)%ilambda
            omega = icontr(mterm)%omega
            spini = icontr(mterm)%spin
            v = icontr(mterm)%v
            vib_count(i) = v
            !
            ! assign vibrational QN v based on the increasing energy for the same State-Sigma-Lambda 
            if (job%assign_v_by_count) then
              !
              loop_ilevel2 : do ilevel = i-1,1,-1
                !
                imaxcontr_ = QNs(ilevel)
                mlevel = transform(irrep)%irec(imaxcontr_)
                mterm_ = ilevel2i(mlevel,irrep)
                istate_ = icontr(mterm_)%istate
                sigmaj = icontr(mterm_)%sigma
                ilambda_ = icontr(mterm_)%ilambda
                !
                if (istate==istate_.and.ilambda==ilambda_.and.nint(sigmaj-sigma)==0) then
                  v = vib_count(ilevel)+1
                  vib_count(i) = v
                  exit loop_ilevel2
                endif
                !
              enddo loop_ilevel2
              !
            endif
            !
            if (iverbose>=3) write(out,'(2x,f8.1,i5,f18.6,1x,i3,2x,2i4,3f8.1,3x,a1,4x,"||",a)') & 
                             jval,i,eigenval(i)-job%ZPE,istate,v,ilambda,spini,sigma,omega,plusminus(irrep), &
                                                  trim(poten(istate)%name)

            if (job%print_rovibronic_energies_to_file ) then
               !  open(unit=u2, file='rovibronic_energies.dat',status='replace',action='write')
                 write(u2,'(2x,f8.1,i5,f20.12,1x,i3,2x,2i4,3f8.1,3x,a1,4x,"||",a)') & 
                     jval,i,eigenval(i)-job%ZPE,istate,v,ilambda,spini,sigma,omega,plusminus(irrep), &
                                                  trim(poten(istate)%name)
               !  close(u2)
             endif
            !
            ! do not count degenerate solutions
            !
            if (i>1.and.abs( eigenval(i)-eigenval(max(i-1,1)) )<job%degen_threshold) cycle
            !
            nener_total = nener_total + 1
            !
            ! for fitting we will need the energies and quantum numbers as a result of this subroutine
            !
            if (present(enerout)) then
               if (nener_total<=size(enerout,dim=3)) then
                  enerout(irot,irrep,nener_total) = eigenval(i)
               endif
            endif
            !
            if (present(quantaout)) then
               if (nener_total<=size(enerout,dim=3)) then
                  !
                  quantaout(irot,irrep,nener_total)%Jrot = Jval
                  quantaout(irot,irrep,nener_total)%irot = irot
                  quantaout(irot,irrep,nener_total)%istate = istate
                  quantaout(irot,irrep,nener_total)%sigma = sigma
                  quantaout(irot,irrep,nener_total)%imulti = imulti
                  quantaout(irot,irrep,nener_total)%ilambda = ilambda
                  quantaout(irot,irrep,nener_total)%omega = omega
                  quantaout(irot,irrep,nener_total)%spin  = spini
                  quantaout(irot,irrep,nener_total)%v = v
                  quantaout(irot,irrep,nener_total)%iparity = irrep-1
                  !
               endif
            endif
          enddo
          !
          if (iverbose>=4) call TimerStop('Assignment')
          !
          if (action%save_eigen_J.or.job%IO_eigen=='SAVE') then
            !
            ! total number of levels for given J,gamma selected for the intensity calculations
            total_roots = 0
            !
            if (iverbose>=4) call TimerStart('Prepare_eigenfuncs_for_intens')
            !
            do i=1,Nroots
              !
              !call Energy_filter(Jval,eigenval(i),irrep,passed)
              !
              if (job%isym_do(irrep)) total_roots = total_roots + 1
              !
            enddo
            !
            total_roots = max(total_roots,1)
            !
            if (action%save_eigen_J) then 
               !
               allocate(eigen(irot,irrep)%vect(Ntotal,total_roots),eigen(irot,irrep)%val(total_roots), & 
                                  eigen(irot,irrep)%quanta(total_roots),stat=alloc)
               call ArrayStart('eigens',alloc,size(eigen(irot,irrep)%vect),kind(eigen(irot,irrep)%vect))
               call ArrayStart('eigens',alloc,size(eigen(irot,irrep)%val),kind(eigen(irot,irrep)%val))
               !
               allocate(basis(irot)%icontr(Ntotal),stat=alloc)
               !
               do i=1,Ntotal
                 !
                 basis(irot)%icontr(i)%istate = icontr(i)%istate
                 basis(irot)%icontr(i)%sigma  = icontr(i)%sigma
                 basis(irot)%icontr(i)%ilambda= icontr(i)%ilambda
                 basis(irot)%icontr(i)%spin   = icontr(i)%spin
                 basis(irot)%icontr(i)%omega  = real(icontr(i)%ilambda,rk)+icontr(i)%sigma
                 basis(irot)%icontr(i)%ivib   = icontr(i)%ivib
                 basis(irot)%icontr(i)%v   = icontr(i)%v
                 !
               enddo
               !
               eigen(irot,irrep)%Nlevels = 0
               eigen(irot,irrep)%Ndimen = Ntotal
               !
               basis(irot)%Ndimen = Ntotal
               !
            endif
            !
            if (intensity%renorm.or.intensity%bound.or.intensity%unbound) then
               allocate(psi_vib(ngrid),vec_t(ngrid),vec0(Ntotal),stat=alloc)
               call ArrayStart('psi_vib',alloc,size(psi_vib),kind(psi_vib))
               call ArrayStart('psi_vib',alloc,size(vec_t),kind(vec_t))
               call ArrayStart('psi_vib',alloc,size(vec0),kind(vec0))               
               psi_vib = 0
            endif
            !
            if (intensity%renorm) then
               !
               rhonorm = sqrt(sqrt(8.0_rk*vellgt*amass*uma/planck))
               !
               write(out,"(/'  Renormalization of unbound states (listing non-converged to sin(kr) at large r)...')")
               write(out,"(6x,'|   # |    J | p | last 3 coeffs. | St vib Lambda Spin     Sigma    Omega ivib|')")
               !
            endif
            !
            total_roots = 0
            !
            if (job%assign_v_by_count) vib_count = 0
            !
            do i=1,Nroots
              !
              !call Energy_filter(Jval,eigenval(i),irrep,passed)
              !
              if (job%isym_do(irrep)) then  
                !
                vec = 0
                !
                do isym = 1,Nsym_
                  !
                  ilevel = transform(irrep)%irec(isym)
                  !
                  do jrrep = 1,sym%NrepresCs
                    !
                    do itau=1,Nirr(ilevel,jrrep)
                      !
                      k = ilevel2i(ilevel,jrrep)
                      !
                      jsym = ilevel2isym(ilevel,irrep)
                      !
                      vec(k) = vec(k) + hsym(isym,i)*Utransform(ilevel,jrrep,irrep)
                      !
                    enddo
                    !
                  enddo
                  !
                enddo
                !
                total_roots = total_roots + 1
                !
                ! to get the assignement we find the term with the largest contribution
                !
                imaxcontr = maxloc(hsym(:,i)**2,dim=1,mask=hsym(:,i)**2.ge.small_)
                !
                ! check if this contribution has already been used and take the second largest coefficient in this case
                !
                maxcontr = hsym(imaxcontr,i)**2
                !
                i0 = 0 
                !
                loop_check_i : do
                  loop_ilevel_i : do ilevel = 1,i-1
                   if (imaxcontr==QNs(ilevel)) then
                     !
                     j = maxloc(hsym(:,i)**2,dim=1,mask=hsym(:,i)**2.ge.small_.and.hsym(:,i)**2.lt.maxcontr)
                     imaxcontr = j
                     maxcontr = hsym(imaxcontr,i)**2
                     i0 = i0 + 1 
                     !
                     if (i0>10) then 
                       imaxcontr = maxloc(hsym(:,i)**2,dim=1,mask=hsym(:,i)**2.ge.small_)
                       exit loop_check_i
                     endif 
                     !
                     cycle loop_check_i
                   endif
                  enddo loop_ilevel_i
                  exit loop_check_i
                enddo loop_check_i
                !
                j = imaxcontr
                QNs(i) = imaxcontr
                !
                mlevel = transform(irrep)%irec(j)
                mterm = ilevel2i(mlevel,irrep)
                !
                istate = icontr(mterm)%istate
                sigma = icontr(mterm)%sigma
                imulti = icontr(mterm)%imulti
                ilambda = icontr(mterm)%ilambda
                omega = icontr(mterm)%omega
                spini = icontr(mterm)%spin
                !v = icontr(mterm)%v
                v  = vib_count(i)
                !vib_count(i) = v
                !
                ! assign vibrational QN v based on the increasing energy for the same State-Sigma-Lambda 
                if (job%assign_v_by_count) then
                  !
                  loop_ilevel3 : do ilevel = i-1,1,-1
                    !
                    imaxcontr_ = QNs(ilevel)
                    mlevel = transform(irrep)%irec(imaxcontr_)
                    mterm_ = ilevel2i(mlevel,irrep)
                    istate_ = icontr(mterm_)%istate
                    sigmaj = icontr(mterm_)%sigma
                    ilambda_ = icontr(mterm_)%ilambda
                    !
                    if (istate==istate_.and.ilambda==ilambda_.and.nint(sigmaj-sigma)==0) then
                      v = vib_count(ilevel)+1
                      vib_count(i) = v
                      exit loop_ilevel3
                    endif
                    !
                  enddo loop_ilevel3
                  !
                endif
                !
                if (action%save_eigen_J) then 
                   !
                   eigen(irot,irrep)%vect(:,total_roots) = vec(:)
                   eigen(irot,irrep)%val(total_roots) = eigenval(i)
                   eigen(irot,irrep)%quanta(total_roots)%istate = istate
                   eigen(irot,irrep)%quanta(total_roots)%sigma = sigma
                   eigen(irot,irrep)%quanta(total_roots)%imulti = imulti
                   eigen(irot,irrep)%quanta(total_roots)%ilambda = ilambda
                   eigen(irot,irrep)%quanta(total_roots)%omega = omega
                   eigen(irot,irrep)%quanta(total_roots)%spin = spini
                   eigen(irot,irrep)%quanta(total_roots)%iparity = irrep-1
                   eigen(irot,irrep)%quanta(total_roots)%igamma = irrep
                   eigen(irot,irrep)%quanta(total_roots)%v = v
                   eigen(irot,irrep)%quanta(total_roots)%name = trim(poten(istate)%name)
                   !
                   !
                   ! This is a special case of the Raman-Wave-Function for which we use set all the 
                   ! eigenvectors back to the basis set
                   !
                   !if (action%RWF.and.istate>Nrefstates) then
                   !    eigen(irot,irrep)%vect(:,total_roots) = 0
                   !    eigen(irot,irrep)%vect(i,total_roots) = 1.0_rk
                   !endif

                   if (intensity%bound.or.intensity%unbound) then
                      !
                      ! funding unboud states
                      !
                      if (iverbose>=4) call TimerStart('Find unbound states')
                      !
                      npoints_last = max(10,grid%npoints/50) 
                      if (npoints_last>=grid%npoints) then
                        write(out,"('wavefunciton unboud check error: too few grid points = ',i7,' use at least 50')") grid%npoints
                        stop 'wavefunciton unboud check error: too few grid points'
                      endif
                      !
                      !$omp parallel do private(k) shared(psi_vib) schedule(guided)
                      do k= grid%npoints-npoints_last+1,grid%npoints
                        psi_vib(k) = vibrational_reduced_density(k,Ntotal,totalroots,Nlambdasigmas,ilambdasigmas_v_icontr,0,&
                        vec,psi_vib)
                      enddo
                      !$omp end parallel do
                      !
                      sum_wv = sum(psi_vib(grid%npoints-npoints_last+1:grid%npoints))
                      !
                      ! condition for the unbound state
                      !
                      if (sum_wv>sqrt(small_)) then 
                         !
                         eigen(irot,irrep)%quanta(total_roots)%bound = .false.
                         !
                      endif
                      !
                      if (iverbose>=4) call TimerStop('Find unbound states')
                      !
                   endif 
                   !
                   if (intensity%renorm) then
                      !
                      ! in order to renormalize the wavefuncitons to sin(kr) at r->infty
                      ! we first need to average the eigenfunction over other degrees of freedom
                      !
                      !psi_vib = 0
                      !vec0(0) = 0 
                      !vec0(1:Ntotal) = vec(1:Ntotal)
                      !
                      if (iverbose>=4) call TimerStart('Reduced vibrational density')
                      !
                      !omp parallel do private(igrid,k,k_,ivib,jvib) shared(psi_vib) schedule(guided)
                      !do igrid=1,grid%npoints
                      !do k = 1,Ntotal
                      !  ivib = icontr(k)%ivib
                      !  do k_ = 1,Ntotal
                      !     jvib = icontr(k_)%ivib
                      !     !
                      !     if (icontr(k)%ilevel==icontr(k_)%ilevel) then 
                      !         !.and.&
                      !         !icontr(k)%istate==icontr(k_)%istate.and.icontr(k)%ilambda==icontr(k_)%ilambda.and.&
                      !         !icontr(k)%sigma==icontr(k_)%sigma) then
                      !         !
                      !         psi_vib(:) = psi_vib(:) + vec(k )*vibrational_contrfunc(:,ivib)*&
                      !                                   vec(k_)*vibrational_contrfunc(:,jvib)
                      !         !
                      !     endif
                      !     !
                      !  enddo
                      !enddo
                      !enddo
                      !omp end parallel do
                      !
                      !omp parallel do private(igrid,ilevel,ivib,k,vec_t,jvib,k_) shared(psi_vib) schedule(guided)
                      !do igrid=1,grid%npoints
                      !do ilevel = 1,Nlambdasigmas
                      !  !
                      !  do ivib =1,totalroots
                      !    !
                      !    k = ilambdasigmas_v_icontr(ivib,ilevel)
                      !    !
                      !    !if (k==0) cycle
                      !    !
                      !    vec_t(:) = vec0(k)*vibrational_contrfunc(:,ivib)
                      !    !
                      !    do jvib =1,totalroots
                      !       !
                      !       k_ = ilambdasigmas_v_icontr(jvib,ilevel)
                      !       !
                      !       !if (k_==0) cycle
                      !       !
                      !       psi_vib(:) = psi_vib(:) + vec_t(:)*vec0(k_)*vibrational_contrfunc(:,jvib)
                      !       !
                      !    enddo
                      !    !
                      !  enddo
                      !enddo
                      !
                      !enddo
                      !omp end parallel do
                      !
                      if (iverbose>=4) call TimerStop('Reduced vibrational density')
                      !
                      ! is this wavefunciton unboud? Check a few last points
                      !
                      npoints_last = max(10,grid%npoints/50) 
                      if (npoints_last>=grid%npoints) then
                        write(out,"('wavefunciton unboud check error: too few grid points = ',i7,' use at least 50')") grid%npoints
                        stop 'wavefunciton unboud check error: too few grid points'
                      endif
                      !
                      !$omp parallel do private(k) shared(psi_vib) schedule(guided)
                      do k= grid%npoints-npoints_last+1,grid%npoints
                        psi_vib(k) = vibrational_reduced_density(k,Ntotal,totalroots,Nlambdasigmas,ilambdasigmas_v_icontr,0,&
                                                                 vec,psi_vib)
                      enddo
                      !$omp end parallel do
                      !
                      sum_wv = sum(psi_vib(grid%npoints-npoints_last+1:grid%npoints))
                      !
                      ! condition for the unbound state
                      !
                      if (sum_wv>sqrt(small_)) then 
                         !
                         if (iverbose>=4) call TimerStart('Find aplitudes of unbound wavefuncs')
                         !
                         ! energy of the unbound state aboove the asympote energy
                         !
                         energy_unbound_sqrsqr = sqrt(sqrt(max(eigenval(i) - poten(istate)%asymptote,0.0_rk)))
                         !
                         ! inspect maxima of the |wavefucntion(r)|^2 at large distances and identify unboud states
                         !
                         psi1 = 0
                         psi2 = psi_vib(grid%npoints)
                         amplit1 = 0 
                         amplit2 = 0
                         amplit3 = 0
                         diff = 0
                         !
                         icount_max = 0
                         !
                         !omp parallel do private(k,psi1,psi2,amplit1,amplit2,amplit3,diff) schedule(guided)
                         loop_gid_dens : do k=grid%npoints-2,max(3,grid%npoints/2),-1
                            !
                            psi1=psi2
                            !
                            psi2= vibrational_reduced_density(k,Ntotal,totalroots,Nlambdasigmas,ilambdasigmas_v_icontr,&
                                                              npoints_last,vec,psi_vib)
                            !psi_vib(k+1)
                            !
                            if ( psi2>100.0*small_.and.psi1<psi2.and.psi2>psi_vib(k) ) then
                               !
                               icount_max = icount_max + 1
                               !
                               amplit1 = amplit2
                               amplit2 = amplit3
                               amplit3 = psi2
                               !
                               diff = abs(amplit2-amplit1)
                               !
                               if (icount_max>2) exit  loop_gid_dens
                               !
                            endif
                         enddo loop_gid_dens
                         !omp end parallel do
                         !
                         if (iverbose>=4) call TimerStop('Find aplitudes of unbound wavefuncs')
                         !
                         if (all((/amplit1,amplit2,amplit3/)>1000*small_)) then
                           !
                           ! now we renormalize wavefunctions that oscilate at large r to 1 at the last amplitude
                           ! and to the density states, see Le Roy J. Chem. Phys. 65, 1485 (1976)
                           !
                           vec(:) = vec(:)*sqrt(amplit3)*rhonorm/energy_unbound_sqrsqr
                           !
                           eigen(irot,irrep)%vect(:,total_roots) = vec(:)
                           !
                           if (diff>1e-3) then 
                             !
                             write(out,'(2x,i8,1x,f8.1,1x,i2,1x,3e12.5,1x,i3,1x,i3,1x,i3,1x,f8.1,1x,f8.1,1x,f8.1,1x,i4)') &
                                        total_roots,J_list(irot),irrep-1,amplit1,amplit2,amplit3,icontr(k)%istate,&
                                        icontr(k)%v,icontr(k)%ilambda,&
                                        icontr(k)%spin,icontr(k)%sigma,icontr(k)%omega,icontr(k)%ivib
                             !
                           endif
                           !
                         endif
                         !
                      endif
                      !
                   endif
                   !
                endif
                !
                if (present(quantaout)) then
                   if (nener_total<=size(enerout,dim=3)) then
                      !
                      ! in order to use the number of nodes as for the vibrational assignement 
                      ! we first need to average the eigenfunction over other degrees of freedom
                      !
                      !psi_vib = 0
                      !!
                      !do k = 1,Ntotal
                      !  ivib = icontr(k)%ivib
                      !  do k_ = 1,Ntotal
                      !     jvib = icontr(k_)%ivib
                      !     !
                      !     if (icontr(k)%istate==icontr(k_)%istate.and.icontr(k)%ilambda==icontr(k_)%ilambda.and.icontr(k)%sigma==icontr(k_)%sigma) then
                      !       !
                      !       psi_vib(:) = psi_vib(:) + vec(k)*vec(k_)*contrfunc(:,ivib)*contrfunc(:,jvib)
                      !       !
                      !     endif
                      !     !
                      !  enddo
                      !enddo
                      !
                      ! Now we can count zeros
                      !
                      !numnod=0
                      !do k=2,ngrid-1
                      !   if (psi_vib(k)<-small_) then
                      !      write(out,"('psi_vib is negative = ',i6,g12.4)") k,psi_vib(k)
                      !      stop 'psi_vib negative' 
                      !   endif
                      !   if ( psi_vib(k)<1e-4.and.&
                      !       (psi_vib(k-1)>sqrt(small_).and.psi_vib(k-1)>sqrt(small_)).and.&
                      !       (psi_vib(k)<psi_vib(k-1).and.psi_vib(k)<psi_vib(k+1)) ) then 
                      !      numnod=numnod+1
                      !   endif
                      !enddo
                      !
                      !v = quantaout(irot,irrep,nener_total)%v 
                      !if (numnod/=v) then 
                      !  write(out,"('total_roots,v,numnod = ',3i8)") total_roots,v,numnod
                      !endif
                      !
                   endif
                endif
                !
                if (job%IO_eigen=='SAVE') then 
                  !
                  do k = 1,Ntotal
                    write(iunit,'(2x,i8,1x,f8.1,1x,i2,1x,e20.12,1x,i3,1x,i3,1x,i3,1x,f8.1,1x,f8.1,1x,f8.1,1x,i4)') &
                          total_roots,J_list(irot),irrep-1,vec(k),icontr(k)%istate,icontr(k)%v,icontr(k)%ilambda,&
                          icontr(k)%spin,icontr(k)%sigma,icontr(k)%omega,icontr(k)%ivib
                  enddo
                  !
                endif
                !
              endif
              !
            enddo
            !
            if (allocated(psi_vib)) then 
               deallocate(psi_vib,vec_t,vec0)
               call ArrayStop('psi_vib')
            endif
            !
            if (action%save_eigen_J) then
              !
              eigen(irot,irrep)%Nlevels = total_roots
              !
            endif
            !
            if (iverbose>=4) call TimerStop('Prepare_eigenfuncs_for_intens')
            !
          endif
          !
          deallocate(QNs)
          call ArrayStop('QNs')
          !
          deallocate(vib_count)
          call ArrayStop('vib_count')
          !
          deallocate(hsym,eigenval)
          call ArrayStop('hsym')
          call ArrayStop('eigenval')
          !
          if (present(nenerout)) nenerout(irot,irrep) = nener_total
          !
       enddo
       !
       if (iverbose>=2) write(out,'(/"Zero point energy (ZPE) = ",f20.12)') job%ZPE
       !
       !omp end parallel do
       !
       !jval = jval + 1.0_rk ; irot = irot + 1
       !
       deallocate(iswap,vec,Nirr)
       call ArrayStop('iswap-vec')
       call ArrayStop('Nirr')
       !
       deallocate(ilevel2i,tau,ilevel2isym)
       !
       deallocate(transform(1)%matrix,transform(2)%matrix,transform(1)%irec,transform(2)%irec)
       call ArrayStop('transform')
       !
       deallocate(Utransform)
       call ArrayStop('Utransform')
       !
       deallocate(hmat)
       call ArrayStop('hmat')
       !
       if (allocated(printout)) deallocate(printout)
       !
       deallocate(icontr)
       !
       if (allocated(ilambdasigmas_v_icontr)) then 
          deallocate(ilambdasigmas_v_icontr)
          call ArrayStop('ilambdasigmas_v_icontr')
       endif
       !
     enddo loop_jval
     !
     if (job%IO_eigen=='SAVE') then 
       !
       write(iunit,"('End of eigenvector')")
       !
       close(unit = iunit, status='keep')
       !  
     endif
     !
     deallocate(kinmat)
     call ArrayStop('kinmat')
     !
     if (allocated(contrenergy)) then 
       deallocate(contrenergy)
       call ArrayStop('contrenergy')
     endif
     !
     if (allocated(icontrvib)) then 
       deallocate(icontrvib)
       call ArrayStop('icontrvib')
     endif
     !
     deallocate(J_list)
     !
     if (allocated(Omega_grid)) then
       !
       do iomega=1,Nomegas
          !
          deallocate(Omega_grid(iomega)%energy,Omega_grid(iomega)%vector,Omega_grid(iomega)%qn,stat=alloc)
          if (alloc/=0) stop 'Omega_grid(iomega)%fields cannot be deallocated'
          deallocate(Omega_grid(iomega)%basis,stat=alloc)
          if (alloc/=0) stop 'Omega_grid(iomega)%basis cannot be deallocated'
          !
       enddo
       !
       call ArrayStop('omegamat_grid')
       !
       deallocate(Omega_grid,stat=alloc)
       if (alloc/=0) stop 'Omega_grid cannot be deallocated'
       !
     endif
     !
     !if (associated(L_omega_obj)) then 
     !  !
     !  do i = 1,NLplus_omega
     !    deallocate(L_omega_obj(i)%gridvalue,stat=alloc)
     !    if (alloc/=0) stop 'L_omega_obj-grid cannot be deallocated'
     !  enddo
     !  deallocate(L_omega_obj,stat=alloc)
     !  if (alloc/=0) stop 'L_omega_obj cannot be deallocated'
     !  call ArrayStop("L+ Omega obj")
     !  !
     !endif
     !
     !if (associated(S_omega_obj)) then 
     !  !
     !  do i = 1,NSplus_omega
     !    deallocate(S_omega_obj(i)%gridvalue,stat=alloc)
     !    if (alloc/=0) stop 'S_omega_obj-grid cannot be deallocated'
     !  enddo
     !  deallocate(S_omega_obj,stat=alloc)
     !  if (alloc/=0) stop 'S_omega_obj cannot be deallocated'
     !  call ArrayStop("S+ Omega obj")
     !endif
     !
     !if (associated(SR_omega_obj)) then 
     !  !
     !  do i = 1,NSR_omega
     !    deallocate(SR_omega_obj(i)%gridvalue,stat=alloc)
     !    if (alloc/=0) stop 'SR_omega_obj-grid cannot be deallocated'
     !  enddo
     !  deallocate(SR_omega_obj,stat=alloc)
     !  if (alloc/=0) stop 'SR_omega_obj cannot be deallocated'
     !  call ArrayStop("SR Omega obj")
     !endif
     !!
     !if (associated(bob_omega_obj)) then 
     !  do i = 1,NBob_omega
     !    deallocate(bob_omega_obj(i)%gridvalue,stat=alloc)
     !    if (alloc/=0) stop 'bob_omega_obj-grid cannot be deallocated'
     !  enddo
     !  deallocate(bob_omega_obj,stat=alloc)
     !  if (alloc/=0) stop 'bob_omega_obj cannot be deallocated'
     !  call ArrayStop("BOB Omega obj")
     !endif
     !!
     !if (associated(p2q_omega_obj)) then 
     !  do i = 1,Np2q_omega
     !    deallocate(p2q_omega_obj(i)%gridvalue,stat=alloc)
     !    if (alloc/=0) stop 'p2q_omega_obj-grid cannot be deallocated'
     !  enddo
     !  deallocate(p2q_omega_obj,stat=alloc)
     !  if (alloc/=0) stop 'p2q_omega_obj cannot be deallocated'
     !  call ArrayStop("P2Q Omega obj")
     !endif
     !
     !if (associated(q_omega_obj)) then 
     !  do i = 1,Nq_omega
     !    deallocate(q_omega_obj(i)%gridvalue,stat=alloc)
     !    if (alloc/=0) stop 'q_omega_obj-grid cannot be deallocated'
     !  enddo
     !  deallocate(q_omega_obj,stat=alloc)
     !  if (alloc/=0) stop 'q_omega_obj cannot be deallocated'
     !  call ArrayStop("Q Omega obj")
     !endif
     !!
     !if (associated(brot)) then 
     !  deallocate(brot(1)%gridvalue,stat=alloc)
     !  if (alloc/=0) stop 'brot-grid cannot be deallocated'
     !  call ArrayStop('BROT')
     !endif
     !
     if (job%print_rovibronic_energies_to_file ) close(u2)
     !
  end subroutine duo_j0


  !
  !  vibrational reduced density of a rovibronic eigenstate at the igrid point
  
  
  function vibrational_reduced_density(igrid,Ntotal,Nvib,Nlambdasigmas,ilambdasigmas_v_icontr,npoints_last,vec,psi_in) &
           result (psi_vib)
     !  
     integer(ik),intent(in) :: igrid,Ntotal,Nvib,Nlambdasigmas,npoints_last,ilambdasigmas_v_icontr(Nvib,Nlambdasigmas)
     real(rk),intent(in)    :: vec(Ntotal),psi_in(grid%npoints)
     real(rk)    :: psi_vib,vec_t
     integer(ik) :: k,k_,ivib,jvib,ilevel
  
        !
        if (igrid>grid%npoints-npoints_last) then 
          !
          psi_vib = psi_in(igrid)
          return 
          !
        endif
        !
        psi_vib = 0 
        do ilevel = 1,Nlambdasigmas
          !
          do ivib =1,Nvib
            !
            k = ilambdasigmas_v_icontr(ivib,ilevel)
            !
            if (k==0) cycle
            !
            vec_t = vec(k)*vibrational_contrfunc(igrid,ivib)
            !
            do jvib =1,Nvib
               !
               k_ = ilambdasigmas_v_icontr(jvib,ilevel)
               !
               if (k_==0) cycle
               !
               psi_vib = psi_vib + vec_t*vec(k_)*vibrational_contrfunc(igrid,jvib)
               !
            enddo
            !
          enddo
        enddo
    
    end function vibrational_reduced_density
    !
    !

     subroutine kinetic_energy_grid_points(ngrid,kinmat,vibTmat,LobWeights,LobDerivs)
        !
        integer(ik),intent(in)  :: ngrid
        real(rk),intent(inout)  :: kinmat(ngrid,ngrid)
        real(rk),intent(inout),optional  :: vibTmat(ngrid,ngrid)
        real(rk),intent(in),optional  :: LobWeights(ngrid),LobDerivs(ngrid,ngrid)
        integer(ik)   :: jgrid,kgrid,igrid
        !
        kinmat = 0
        !
        do igrid =1, ngrid
          !
          select case(solution_method)
          case ("5POINTDIFFERENCES")
           !
           kinmat(igrid,igrid) = kinmat(igrid,igrid) + d2dr(igrid)
           !
           ! The nondiagonal matrix elemenets are:
           ! The vibrational kinetic energy operator will connect only the
           ! neighbouring grid points igrid+/1 and igrid+/2.
           !
           ! Comment by Lorenzo Lodi
           ! The following method corresponds to approximating the second derivative of the wave function
           ! psi''  by the 5-point finite difference formula:
           !
           ! f''(0) = [-f(-2h) +16*f(-h) - 30*f(0) +16*f(h) - f(2h) ] / (12 h^2)  + O( h^4 )
           !
           if (igrid>1) then
             kinmat(igrid,igrid-1) = -16.0_rk*z(igrid-1)*z(igrid)
             kinmat(igrid-1,igrid) = kinmat(igrid,igrid-1)
           endif
           !
           if (igrid>2) then
             kinmat(igrid,igrid-2) = z(igrid-2)*z(igrid)
             kinmat(igrid-2,igrid) = kinmat(igrid,igrid-2)
           endif
           !
           case("SINC")   ! Colbert Miller sinc DVR (works only for uniform grids at the moment)
                          ! This is the `simple' sinc DVR version, derived for the range (-infty, +infty).
             if( grid%nsub /=0) then
               write(out, '(A)') 'Sorry, at the moment only uniformely-spaced grids (type 0) can be used with the SINC method.'
               write(out, '(A)') 'Use 5PointDifferences as solution method for non uniformely-spaced grids.'
               stop
             endif
             kinmat(igrid,igrid) = kinmat(igrid,igrid) +(12._rk)* pi**2 / 3._rk
             !
             do jgrid =igrid+1, ngrid
               kinmat(igrid,jgrid) = +(12._rk)*2._rk* real( (-1)**(igrid+jgrid), rk) / real(igrid - jgrid, rk)**2
               kinmat(jgrid,igrid) = kinmat(igrid,jgrid)
             enddo
             !
           case("LOBATTO") ! Implements a DVR method based on Lobatto quadrature
                           ! Requires the Lobatto nonuniform grid to work
             if(grid%nsub /= 6) then
               write(out, '(A)') 'The Lobatto DVR method only works with the'
               write(out, '(A)') 'Lobatto grid (grid type 6).'
               stop
             endif
             !
             if (.not.present(vibTmat)) stop 'vibTmat is not present'
             !
             vibTmat = 0 
             !
             do jgrid = igrid,ngrid
                do kgrid=1,ngrid
                   vibTmat(igrid,jgrid) = vibTmat(igrid,jgrid) + (12._rk)*(hstep**2)*(LobWeights(kgrid))*& 
                                          LobDerivs(igrid,kgrid)*LobDerivs(jgrid,kgrid)
                enddo
                vibTmat(jgrid,igrid) = vibTmat(igrid,jgrid)
                kinmat(igrid,jgrid) = kinmat(igrid,jgrid) + vibTmat(igrid,jgrid)
                kinmat(jgrid,igrid) = kinmat(igrid,jgrid)
             enddo
             !
           case default
             write(out, '(A)') 'Error: unrecognized solution method' // trim(solution_method)
             write(out, '(A)') 'Possible options are: '
             write(out, '(A)') '                      5POINTDIFFERENCES'
             write(out, '(A)') '                      SINC'
             write(out, '(A)') '                      LOBATTO'
           end select
           !
        enddo
        !     
     end subroutine kinetic_energy_grid_points





   subroutine Transfrorm_Sigma_Lambda_to_Omega_representation(iverbose,sc,Nlambdasigmas_max,Nomega_states)
     !
     use lapack,only : lapack_syev,lapack_heev,lapack_syevr
     !
     implicit none 
     !
     integer(ik),intent(in) ::iverbose,Nlambdasigmas_max,Nomega_states
     real(rk),intent(in)    :: sc
     !
     integer(ik) :: omega_min,omega_max,iomega,jomega
     integer(ik) :: igrid,jgrid,istate,jstate,imulti,jmulti,ilambda,jlambda,iL2,ieq,ispin,nspins
     integer(ik) :: i,j,idiab,ipermute,istate_,jstate_,ilambda_we,jlambda_we,isigma2,isigmav,itau,N_i,N_j
     integer(ik) :: alloc,ngrid,Nlambdasigmas,iso,ibob,ilambda_,jlambda_,ilxly,iLplus_omega_,iomega_count
     integer(ik) :: multi_max,iSplus_omega_,iSR_omega_,ibob_omega_,ip2q_omega_,iq_omega_,iKin_omega_,iBRot_omega_
     integer(ik) :: imaxcontr,isr,iss,isso,ild,ip2q,iq,ibobrot
     !
     real(rk)    :: f_rot,omegai,omegaj,sigmai,sigmaj,spini,spinj,epot,f_l2,sigmai_we,sigmaj_we,spini_,spinj_,q_we,f_centrif
     real(rk)    :: three_j_ref,three_j_,SO,omegai_,omegaj_,f_grid,f_s,b_rot,erot,f_diabatic
     real(rk)    :: sigmai_,sigmaj_,f_t,spin_min,f_sr,f_ss,f_s1,f_s2,f_lo
     !
     type(fieldT),pointer      :: field
     !
     real(rk), allocatable :: mat_1(:,:),mat_2(:,:),vect(:,:,:)
     real(rk), allocatable :: omegamat(:,:),omegaenergy(:)
     real(rk), allocatable :: L_LambdaSigma(:,:)
     integer(ik),allocatable :: iomega_state(:,:),imax_contr(:,:)
     real(rk),allocatable    :: vibmat(:,:),kinmat(:,:)
       !
       ngrid = grid%npoints
       !
       b_rot = aston/amass
       !
       iomega = 0
       omega_min = 10000
       omega_max = -10000
       multi_max = 1
       !
       do istate = 1,nestates
         multi_max = max(poten(istate)%multi,multi_max)
       enddo
       !
       spin_min = 0
       if (mod(nint(2.0_rk*multi_max+1.0_rk),2)==1) spin_min = 0.5
       !
       Nspins = nint(real(multi_max-1,rk)*0.5_rk)
       !
       allocate(L_LambdaSigma(Nomega_states,Nomega_states),stat=alloc)
       call ArrayStart('L_LambdaSigma',alloc,size(L_LambdaSigma),kind(L_LambdaSigma))
       !
       allocate(iomega_state(Nomegas,Nlambdasigmas_max),stat=alloc)
       call ArrayStart('iomega_state',alloc,size(iomega_state),kind(iomega_state))
       !
       iomega_state = 0
       iomega_count = 0
       !
       do iomega=1,Nomegas
          !
          Nlambdasigmas = Omega_grid(iomega)%Nstates
          !
          do i = 1,Nlambdasigmas
            !
            iomega_count = iomega_count + 1
            iomega_state(iomega,i) = iomega_count
            !
          enddo
          !
       enddo
       !
       allocate(mat_1(Nlambdasigmas_max,Nlambdasigmas_max),stat=alloc)
       call ArrayStart('mat_1',alloc,size(mat_1),kind(mat_1))
       allocate(mat_2(Nlambdasigmas_max,Nlambdasigmas_max),stat=alloc)
       call ArrayStart('mat_2',alloc,size(mat_2),kind(mat_2))
       !
       ! Vector at the previous step to keep track of the phase of the wavefunctions
       allocate(vect(Nomegas,Nlambdasigmas_max,Nlambdasigmas_max),stat=alloc)
       call ArrayStart('vect',alloc,size(vect),kind(vect))
       !
       ! index to keep track of the state of the wavefunctions
       allocate(imax_contr(Nomegas,Nlambdasigmas_max),stat=alloc)
       call ArrayStart('imax_contr',alloc,size(imax_contr),kind(imax_contr))
       !
       ! For each grid point diagonalise the Sigma-Lambda PECs + SOs and transform to the Omega-represenation
       do igrid =1, ngrid
          !
          i = (istate-1)*ngrid + igrid
          !
          if (iverbose>=6) write(out,'("igrid = ",i0)') igrid
          !
          ! the centrifugal factor will be needed for the L**2 term
          !
          f_centrif=b_rot/r(igrid)**2*sc
          !
          L_LambdaSigma = 0 
          !
          do iomega=1,Nomegas
             !
             omegai = Omega_grid(iomega)%omega
             !
             Nlambdasigmas = Omega_grid(iomega)%Nstates
             !
             if (Nlambdasigmas==0) cycle
             !
             allocate(omegamat(Nlambdasigmas,Nlambdasigmas),omegaenergy(Nlambdasigmas),stat=alloc)
             if (alloc/=0) stop 'Cannot allocate omegamat and omegaenergy'
             omegamat = 0
             !
             do i = 1,Nlambdasigmas
               !
               istate  = Omega_grid(iomega)%basis(i)%istate
               sigmai  = Omega_grid(iomega)%basis(i)%sigma
               imulti  = Omega_grid(iomega)%basis(i)%imulti
               ilambda = Omega_grid(iomega)%basis(i)%ilambda
               spini   = Omega_grid(iomega)%basis(i)%spin
               !
               if (omegai/=Omega_grid(iomega)%basis(i)%omega) then
                 print*,'something wrong with omegas',omegai,quanta(i)%omega
                 stop 'something wrong with omega'
               endif
               !
               epot = poten(istate)%gridvalue(igrid)*sc
               !
               ! Another diagonal term:
               ! The L^2 term (diagonal): (1) L2(R) is used if provided otherwise
               ! an approximate value Lambda*(Lamda+1) is assumed.
               !
               f_l2 = 0 ! real(ilambda*(ilambda+1),rk)*f_rot
               do iL2 = 1,Nl2
                 if (L2(iL2)%istate==istate.and.L2(iL2)%jstate==istate) then
                   f_l2 = f_rot*L2(iL2)%gridvalue(igrid)
                   exit
                 endif
               enddo
               !
               erot = f_l2
               !
               f_rot=f_centrif
               !
               ! BOB centrifugal (rotational) term, i.e. a correction to f_rot
               !
               do ibobrot = 1,Nbobrot
                 if (bobrot(ibobrot)%istate==istate.and.bobrot(ibobrot)%jstate==jstate.and.istate==jstate) then
                   field => bobrot(ibobrot)
                   f_rot = f_rot + field%gridvalue(igrid)*sc
                   exit
                 endif
               enddo
               !
               ! rotational diagonal element 
               !                                             ! L Lodi -job%diag_L2_fact is either zero or one
               !erot = erot + f_rot*(  - omegai**2 -job%diag_L2_fact*real(ilambda**2,rk)  & 
               !         +   spini*(spini+1.0_rk) - sigmai**2 )
               !
               erot = erot + f_rot*(  -job%diag_L2_fact*real(ilambda**2,rk)  & 
                        +   spini*(spini+1.0_rk) - sigmai**2 )
               !
               ! add the diagonal matrix element to the local spin-rotational matrix hmat
               !
               omegamat(i,i) = epot+erot
               !
               !
               ! Diagonal spin-rotation term
               !
               do isr = 1,Nsr
                 if (spinrot(isr)%istate==istate.and.spinrot(isr)%jstate==istate) then
                   field => spinrot(isr)
                   f_sr = field%gridvalue(igrid)*(sigmai**2-spini*(spini+1.0_rk))*sc
                   omegamat(i,i) = omegamat(i,i) + f_sr
                   exit
                 endif
               enddo
               !
               ! Diagonal spin-spin term
               !
               do iss = 1,Nss
                 if (spinspin(iss)%istate==istate.and.spinspin(iss)%jstate==istate) then
                   field => spinspin(iss)
                   f_ss = field%gridvalue(igrid)*(3.0_rk*sigmai**2-spini*(spini+1.0_rk))*sc
                   omegamat(i,i) = omegamat(i,i) + f_ss
                   exit
                 endif
               enddo
               !
               do j = i,Nlambdasigmas
                   !
                   jstate  = Omega_grid(iomega)%basis(j)%istate
                   sigmaj  = Omega_grid(iomega)%basis(j)%sigma
                   jmulti  = Omega_grid(iomega)%basis(j)%imulti
                   jlambda = Omega_grid(iomega)%basis(j)%ilambda
                   spinj   = Omega_grid(iomega)%basis(j)%spin
                   !
                   if (omegai/=Omega_grid(iomega)%basis(i)%omega) stop 'something wrong with omega'
                   !
                   if (iverbose>=6) write(out,'("i,j = ",2(i0,2x) )') i,j
                   !
                   ! Diabatic non-diagonal contribution  term
                   !
                   do idiab = 1,Ndiabatic
                     if (diabatic(idiab)%istate==istate.and.diabatic(idiab)%jstate==jstate.and.&
                         abs(nint(sigmaj-sigmai))==0.and.(ilambda==jlambda).and.nint(spini-spinj)==0 ) then
                       field => diabatic(idiab) 
                       f_diabatic = field%gridvalue(igrid)*sc
                       omegamat(i,j) = omegamat(i,j) + f_diabatic
                       omegamat(j,i) = omegamat(i,j)
                       exit
                     endif
                   enddo
                   !
                   ! spin-orbit part:
                   loop_iso_omega : do iso =1,Nspinorbits
                     !
                     field => spinorbit(iso)
                     !
                     ! The selection rules are (Lefebvre-Brion and Field, Eq. (3.4.6)): 
                     ! Delta J = 0 ; Delta Omega  = 0 ; g<-/->u; e<->f; Sigma+<->Sigma-;
                     ! Delta S = 0 or Delta S = 1 ; Delta Lambda = Delta Sigma = 0 or Delta Lambda = - Delta Sigma = +/- 1
                     !
                     if (nint(spini-spinj)>1 ) cycle
                     if ( ilambda==0.and.jlambda==0.and.poten(istate)%parity%pm==poten(jstate)%parity%pm ) cycle
                     if ( poten(istate)%parity%gu/=0.and.poten(istate)%parity%gu/=poten(jstate)%parity%gu ) cycle
                     !
                     do ipermute  = 0,1
                       !
                       if (ipermute==0) then
                         !
                         istate_ = field%istate ; ilambda_we = field%lambda  ; sigmai_we = field%sigmai ; spini_ = field%spini
                         jstate_ = field%jstate ; jlambda_we = field%lambdaj ; sigmaj_we = field%sigmaj ; spinj_ = field%spinj
                         !
                       else  ! permute
                         !
                         jstate_ = field%istate ; jlambda_we = field%lambda  ; sigmaj_we = field%sigmai ; spinj_ = field%spini
                         istate_ = field%jstate ; ilambda_we = field%lambdaj ; sigmai_we = field%sigmaj ; spini_ = field%spinj
                         !
                       endif
                       !
                       ! proceed only if the spins of the field equal the corresponding <i| and |j> spins of the current matrix elements. 
                       ! otherwise skip it:
                       if ( nint(spini_-spini)/=0.or.nint(spinj_-spinj)/=0 ) cycle
                       !
                       ! however the permutation makes sense only when for non diagonal <State,Lambda,Spin|F|State',Lambda',Spin'>
                       ! otherwise it will cause a double counting:
                       !
                       if (ipermute==1.and.istate_==jstate_.and.ilambda_we==jlambda_we.and.nint(sigmai_we-sigmaj_we)==0.and. & 
                           nint(spini_-spinj_)==0) cycle
                       !
                       ! check if we at the right electronic states
                       if( istate/=istate_.or.jstate/=jstate_ ) cycle
                       !
                       ! We apply the Wigner-Eckart theorem to reconstruct all combinations of <Lamba Sigma |HSO|Lamba Sigma' > 
                       ! connected with the reference (input) <Lamba Sigma_r |HSO|Lamba Sigma_r' > by this theorem. 
                       ! Basically, we loop over Sigma (Sigma = -S..S).  The following 3j-symbol for the reference states will 
                       ! be conidered:
                       ! / Si      k  Sj     \    k  = 1
                       ! \ -Sigmai q  Sigmaj /    q  = Sigmai - Sigmaj
                       !
                       ! reference q from Wigner-Eckart
                       q_we = sigmai_we-sigmaj_we
                       !
                       ! We should consider also a permutation <State',Lambda',Spin'|F|State,Lambda,Spin> if this makes a change.
                       ! This could be imortant providing that we constrain the i,j indexes to be i<=j (or i>=j).
                       ! We also assume that the matrix elements are real!
                       !
                       ! First of all we can check if the input values are not unphysical and consistent with Wigner-Eckart:
                       ! the corresponding three_j should be non-zero:
                       three_j_ref = three_j(spini_, 1.0_rk, spinj_, -sigmai_we, q_we, sigmaj_we)
                       !
                       if (abs(three_j_ref)<small_) then 
                         !
                         write(out,"('The Spin-orbit field (J=0) ',2i3,' is incorrect according to Wigner-Eckart, three_j = 0 ')") & 
                               field%istate,field%jstate
                         write(out,"('Check S_i, S_j, Sigma_i, Sigma_j =  ',4f9.2)") spini_,spinj_,sigmai_we,sigmaj_we
                         stop "The S_i, S_j, Sigma_i, Sigma_j are inconsistent"
                         !
                       end if 
                       !
                       ! Also check the that the SO is consistent with the selection rules for SO
                       !
                       if ( ilambda_we-jlambda_we+nint(sigmai_we-sigmaj_we)/=0.or.nint(spini_-spinj_)>1.or.&
                          ( ilambda_we==0 .and. jlambda_we==0 .and. &
                            poten(field%istate)%parity%pm==poten(field%jstate)%parity%pm ).or. &
                          ( (ilambda_we-jlambda_we)/=-nint(sigmai_we-sigmaj_we) ).or.&
                             abs(ilambda_we-jlambda_we)>1.or.abs(nint(sigmai_we-sigmaj_we))>1.or. &
                          ( poten(field%istate)%parity%gu/=0.and.&
                          poten(field%istate)%parity%gu/=poten(field%jstate)%parity%gu ) ) then
                          !
                          write(out,"('The quantum numbers of the spin-orbit field (J=0) ',2i3,' are inconsistent" // &
                                          " with SO selection rules: ')") field%istate,field%jstate
                          write(out,"('Delta J = 0 ; Delta Omega  = 0 ; g<-/->u; e<-/->f; Sigma+<->Sigma-; " // &
                             "Delta S = 0 or Delta S = 1 ; Delta Lambda = Delta Sigma = 0 " // &
                                "or Delta Lambda = - Delta Sigma = +/- 1')")
                          write(out,"('Check S_i, S_j, Sigma_i, Sigma_j, lambdai, lambdaj =  ',4f9.2,2i4)") &
                                                                spini_,spinj_,sigmai_we,sigmaj_we,ilambda_we,jlambda_we
                          stop "The S_i, S_j, Sigma_i, Sigma_j lambdai, lambdaj are inconsistent with selection rules"
                          !
                       endif
                       !
                       do isigma2 = -nint(2.0*spini_),nint(2.0*spini_),2
                         !
                         ! Sigmas from Wigner-Eckart
                         sigmai_ = real(isigma2,rk)*0.5 
                         sigmaj_ = sigmai_ - q_we
                         !
                         ! three_j for current Sigmas
                         three_j_ = three_j(spini_, 1.0_rk, spinj_, -sigmai_, q_we, sigmaj_)
                         !
                         ! current value of the SO-matrix element from Wigner-Eckart
                         SO = (-1.0_rk)**(sigmai_-sigmai_we)*three_j_/three_j_ref*field%gridvalue(igrid)
                         !
                         ! We should also take into account that Lambda and Sigma can change sign
                         ! since in the input we give only a unique combination of matrix elements, for example
                         ! < 0 0 |  1  1 > is given, but < 0 0 | -1 -1 > is not, assuming that the program will generate the missing
                         ! combinations.
                         !
                         ! In order to recover other combinations we apply the symmetry transformation
                         ! laboratory fixed inversion which is equivalent to the sigmav operation 
                         !                    (sigmav= 0 correspond to the unitary transformation)
                         do isigmav = 0,1
                           !
                           ! sigmav is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
                           ! avoid the double counting.
                           if( isigmav==1.and.&
                             nint( abs( 2.0*sigmai_ )+ abs( 2.0*sigmaj_ ) )+abs( ilambda_we )+abs( jlambda_we )==0 &
                           )cycle
                           !
                           ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                           ilambda_ = ilambda_we*(-1)**isigmav
                           jlambda_ = jlambda_we*(-1)**isigmav
                           sigmai_ = sigmai_*(-1.0_rk)**isigmav
                           sigmaj_ = sigmaj_*(-1.0_rk)**isigmav
                           !
                           omegai_ = sigmai_+real(ilambda_)
                           omegaj_ = sigmaj_+real(jlambda_)
                           !
                           ! Check So selection rules
                           if ( ( ilambda_-jlambda_)/=-nint(sigmai_-sigmaj_).or. &
                             abs(sigmai_-sigmaj_)>1.or.omegai_/=omegaj_ ) cycle
                           !
                           ! proceed only if the quantum numbers of the field equal
                           ! to the corresponding <i| and |j> quantum numbers of the basis set. otherwise skip it:
                           if ( nint(sigmai_-sigmai)/=0.or.nint(sigmaj_-sigmaj)/=0 &
                             .or.ilambda_/=ilambda.or.jlambda_/=jlambda ) cycle
                           !
                           f_t = SO*sc
                           !
                           ! the result of the symmetry transformtaion applied to the <Lambda,Sigma|HSO|Lambda',Sigma'> only
                           if (isigmav==1) then
                             !
                             ! still not everything is clear here: CHECK!
                             !
                             itau = -ilambda_-jlambda_ +nint(spini_-sigmai_)+nint(spinj_-sigmaj_) !+nint(jval-omegai)+(jval-omegaj)
                             !
                             !itau = nint(spini_-sigmai_)+nint(spinj_-sigmaj_) ! +nint(jval-omegai)+(jval-omegaj)
                             !
                             !itau = 0
                             !
                             if (ilambda_==0.and.poten(istate)%parity%pm==-1) itau = itau+1
                             if (jlambda_==0.and.poten(jstate)%parity%pm==-1) itau = itau+1
                             !
                             f_t = f_t*(-1.0_rk)**(itau)
                             !
                           endif
                           !
                           ! double check
                           if ( nint(omegai-omegai_)/=0 .or. nint(omegai-omegaj_)/=0 ) then
                             write(out,'(A,f8.1," or omegaj ",f8.1," do not agree with stored values ",f8.1,1x,f8.1)') &
                                        "SO: reconsrtucted omegai", omegai_,omegaj_,omegai
                             stop 'SO: wrongly reconsrtucted omegai or omegaj'
                           endif
                           !
                           ! we might end up in eilther parts of the matrix (upper or lower),
                           ! so it is safer to be general here and
                           ! don't restrict to lower part as we have done above
                           !
                           omegamat(i,j) = omegamat(i,j) + f_t
                           !
                           omegamat(j,i) = omegamat(i,j)
                           !
                           cycle loop_iso_omega
                           !
                         enddo
                       enddo
                     enddo
                   enddo  loop_iso_omega
                   !
                   ! L*S
                   !
                   loop_ilxly_omega : do ilxly =1,Nlxly
                     !
                     field => lxly(ilxly)
                     !
                     ! Also check that L+ is consistent with the selection rules
                     !
                     if ( field%istate==field%jstate .or.abs(field%lambda-field%lambdaj)/=1 ) then
                        !
                        write(out,"('The quantum numbers of the L+/Lx field ',2i3,' are inconsistent" // &
                                        " with L+selection rules: ')") field%istate,field%jstate
                        write(out,"('Delta Lamda = +/-1')")
                        stop "Lx/L+ input is inconsistent with selection rules"
                        !
                     endif
                     !
                     ! the field entry in the input gives only one combination of the quantum numbers for
                     ! the matrix element <State,Lambda,Spin|F|State',Lambda',Spin'>
                     ! LxLy  4 6 ;  lambda  0 1 ; spin   1.0 1.0
                     ! we should consider also a permutation <State',Lambda',Spin'|F|State,Lambda,Spin> if this makes a change.
                     ! This could be imortant providing that we constrain the i,j indexes to be i<=j (or i>=j)
                     ! We also assume that the matrix elements are real!
                     !
                     do ipermute  = 0,1
                       !
                       if (ipermute==0) then
                         !
                         istate_ = field%istate ; ilambda_ = field%lambda  ; spini_ = field%spini
                         jstate_ = field%jstate ; jlambda_ = field%lambdaj ; spinj_ = field%spinj
                         !
                       else  ! permute
                         !
                         jstate_ = field%istate ; jlambda_ = field%lambda  ; spinj_ = field%spini
                         istate_ = field%jstate ; ilambda_ = field%lambdaj ; spini_ = field%spinj
                         !
                       endif
                       !
                       ! however the permutation makes sense only when for non diagonal <State,Lambda,Spin|F|State',Lambda',Spin'>
                       ! otherwise it will cause a double counting:
                       !
                       if (ipermute==1.and.istate_==jstate_.and.ilambda_==jlambda_.and.nint(spini_-spinj_)==0) cycle
                       !
                       ! check if we at the right electronic states
                       if( istate/=istate_.or.jstate/=jstate_ ) cycle
                       !
                       ! We should also take into account that Lambda can change sign (only Lambda>0 is given in input)
                       ! In order to recover other combinations we apply the symmetry transformation
                       ! laboratory fixed inversion which is equivalent to the sigmav operation 
                       !                    (sigmav= 0 correspond to the unitary transformation)
                       do isigmav = 0,1
                         !
                         ! the permutation is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
                         ! avoid the double counting.
                         if( isigmav==1.and. abs( field%lambda ) + abs( field%lambdaj )==0 ) cycle
           
                         ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                         ilambda_ = ilambda_*(-1)**isigmav
                         jlambda_ = jlambda_*(-1)**isigmav
                         !
                         ! proceed only if the quantum numbers of the field equal to the corresponding <i| and |j> quantum numbers:
                         if (istate/=istate_.or.jstate_/=jstate.or.ilambda_/=ilambda.or.jlambda_/=jlambda) cycle
                         !
                         ! check the selection rule Delta Lambda = +/1
                         if (abs(ilambda-jlambda)/=1) cycle
                         !
                         ! double check
                         if (spini/=poten(istate)%spini.or.spinj/=poten(jstate)%spini) then
                          write(out,&
                            '("LS: reconstructed spini ",f8.1," or spinj ",f8.1," do not agree with stored values ", & 
                                     & f8.1,1x,f8.1)') spini,spinj,poten(istate)%spini,poten(jstate)%spini
                           stop 'LS: wrongly reconsrtucted spini or spinj'
                         endif
                         !
                         f_grid  = field%gridvalue(igrid)
                         !
                         ! <Lx> and <Ly> don't depend on Sigma
                         !
                         ! L*S part (spin-electronic coupling)
                         !
                         ! the selection rules are Delta Sigma = - Delta Lambda (Delta Spin = 0)
                         !
                         ! factor to switch between <Sigma+1|S+|Sigma> and <Sigma-1|S-|Sigma>:
                         f_s = real(ilambda-jlambda,rk)
                         !
                         ! the bra-component of Sigma (i.e. sigmaj):
                         sigmaj_ = sigmai+f_s
                         !
                         ! make sure that this sigmaj_ is consistent with the current ket-sigmaj
                         if (nint(2.0_rk*sigmaj_)==nint(2.0*sigmaj)) then
                           !
                           f_t = f_grid*f_rot
                           !
                           ! the result of the symmetry transformation:
                           if (isigmav==1) then
                             !
                             ! we assume that
                             ! sigmav <Lamba|L+|Lambda'> => <-Lamba|L-|-Lambda'> == <Lamba|L+|Lambda'>(-1)^(Lamba+Lambda')
                             ! and <Lamba|L+|Lambda'> is an unique quantity given in the input
                             ! however we don't apply the sigmav transformation to sigma or omega
                             ! since we only need to know how <Lamba|L+/-|Lambda'> transforms in order to relate it to the
                             ! value given in input.
                             !
                             itau = 0 
                             !
                             if (ilambda_==0.and.poten(istate)%parity%pm==-1) itau = itau+1
                             if (jlambda_==0.and.poten(jstate)%parity%pm==-1) itau = itau+1
                             !
                             f_t = f_t*(-1.0_rk)**(itau)
                             !
                           endif
                           !
                           ! the matrix element <Sigmai| S+/- |Sigmai+/-1>
                           !
                           f_t = sqrt( (spini-f_s*sigmai)*(spini+f_s*sigmai+1.0_rk) )*f_t
                           !
                           omegamat(i,j) = omegamat(i,j) + f_t
                           omegamat(j,i) = omegamat(i,j)
                           !
                         endif
                         !
                       enddo
                       !
                     enddo
                     !
                   enddo loop_ilxly_omega
                   !
                   !
                   ! Non-diagonal spin-rotaion term (J-independent part only)
                   !
                   do isr = 1,Nsr
                     !
                     ! non-diagonal part
                     if (i==j) cycle
                     !
                     field => spinrot(isr)
                     !
                     ! Only the part which does not depend on J: 
                     !
                     ! 2. <Sigma,Omega,Lambda|HSR|Sigma+/-1,Omega,Lambda-/+>
                     ! with the effective parameter gamma_v including the matrix element <Lambda|L+/-|lambda-/+1>
                     if (spinrot(isr)%istate==istate.and.spinrot(isr)%jstate==jstate.and.&
                         abs(nint(sigmaj-sigmai))==1.and.abs(ilambda-jlambda)==1.and.nint(spini-spinj)==0) then
                         !
                         do ipermute  = 0,1
                           !
                           if (ipermute==0) then
                             !
                             istate_ = field%istate ; ilambda_ = field%lambda  
                             jstate_ = field%jstate ; jlambda_ = field%lambdaj 
                             !
                           else  ! permute
                             !
                             jstate_ = field%istate ; jlambda_ = field%lambda 
                             istate_ = field%jstate ; ilambda_ = field%lambdaj
                             !
                           endif
                           !
                           ! however the permutation makes sense only when for non diagonal <State,Lambda,Spin|F|State',Lambda',Spin'>
                           ! otherwise it will cause a double counting:
                           !
                           if (ipermute==1.and.istate_==jstate_.and.ilambda_==jlambda_) cycle
                           !
                           ! check if we at the right electronic states
                           if( istate/=istate_.or.jstate/=jstate_ ) cycle
                           !
                           ! We should also take into account that Lambda can change sign (only Lambda>0 is given in input)
                           ! In order to recover other combinations we apply the symmetry transformation
                           ! laboratory fixed inversion which is equivalent to the sigmav operation 
                           !                    (sigmav= 0 correspond to the unitary transformation)
                           do isigmav = 0,1
                             !
                             ! the permutation is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
                             ! avoid the double counting.
                             if( isigmav==1.and. abs( field%lambda ) + abs( field%lambdaj )==0 ) cycle
                      
                             ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                             ilambda_ = ilambda_*(-1)**isigmav
                             jlambda_ = jlambda_*(-1)**isigmav
                             !
                             ! proceed only if the quantum numbers of the field equal to the corresponding <i| and |j> quantum numbers:
                             if (ilambda_/=ilambda.or.jlambda_/=jlambda) cycle
                             !
                             !
                             ! double check
                             if (spini/=poten(istate)%spini.or.spinj/=poten(jstate)%spini) then
                              write(out,'("SR: reconstructed spini ",f8.1," or spinj ",f8.1," do not agree with stored values ", & 
                                         & f8.1,1x,f8.1)') spini,spinj,poten(istate)%spini,poten(jstate)%spini
                               stop 'SR: wrongly reconsrtucted spini or spinj'
                             endif
                             !
                             f_grid  = field%gridvalue(igrid)*sc
                             !
                             ! <Lx> and <Ly> don't depend on Sigma
                             !
                             ! L*S part of the spin-rotation 
                             !
                             ! the selection rules are Delta Sigma = - Delta Lambda (Delta Spin = 0)
                             !
                             ! factor to switch between <Sigma+1|S+|Sigma> and <Sigma-1|S-|Sigma>:
                             f_s = real(ilambda-jlambda,rk)
                             !
                             ! the bra-component of Sigma (i.e. sigmaj):
                             sigmaj_ = sigmai+f_s
                             !
                             ! make sure that this sigmaj_ is consistent with the current ket-sigmaj
                             if (nint(2.0_rk*sigmaj_)==nint(2.0*sigmaj)) then
                               !
                               f_t = f_grid
                               !
                               ! the result of the symmetry transformation:
                               if (isigmav==1) then
                                 !
                                 ! we assume that
                                 ! sigmav <Lamba|L+|Lambda'> => <-Lamba|L-|-Lambda'> == <Lamba|L+|Lambda'>(-1)^(Lamba+Lambda')
                                 ! and <Lamba|L+|Lambda'> is an unique quantity given in the input
                                 ! however we don't apply the sigmav transformation to sigma or omega
                                 ! since we only need to know how <Lamba|L+/-|Lambda'> transforms in order to relate it to the
                                 ! value given in input.
                                 !
                                 itau = 0 
                                 !
                                 if (ilambda_==0.and.poten(istate)%parity%pm==-1) itau = itau+1
                                 if (jlambda_==0.and.poten(jstate)%parity%pm==-1) itau = itau+1
                                 !
                                 f_t = f_t*(-1.0_rk)**(itau)
                                 !
                               endif
                                !
                               ! the matrix element <Sigmai| S+/- |Sigmai+/-1>
                               !
                               f_t = sqrt( (spini-f_s*sigmai)*(spini+f_s*sigmai+1.0_rk) )*f_t
                               !
                               !f_t = sqrt( spini*(spini+1.0_rk)-sigmai*(sigmai+f_s)  )*f_t
                               !
                               omegamat(i,j) = omegamat(i,j) + f_t
                               omegamat(j,i) = omegamat(i,j)
                               !
                             endif
                             !
                           enddo
                           !
                         enddo
                         !
                     endif
                     ! 
                   enddo
                   !
                   ! Non-diagonal spin-spin term
                   !
                   do isso = 1,Nsso
                     !
                     if (spinspino(isso)%istate==istate.and.spinspino(isso)%jstate==jstate.and.istate==jstate.and.&
                         abs(nint(sigmaj-sigmai))==1.and.(ilambda-jlambda)==nint(sigmaj-sigmai)) then
                        !
                        field => spinspino(isso)
                        !
                        f_s = sigmaj-sigmai
                        !
                        f_t = sqrt( spini*(spini+1.0_rk)-(sigmai+0.5_rk*f_s)*(sigmai    ) )*&
                              sqrt( spini*(spini+1.0_rk)-(sigmai+0.5_rk*f_s)*(sigmai+f_s) )
                        !
                        f_ss = field%gridvalue(igrid)*f_t*sc
                        !
                        omegamat(i,j) = omegamat(i,j) + f_ss
                        omegamat(j,i) = omegamat(i,j)
                        !
                     endif
                     ! 
                   enddo ! S-S
                   !
                   !
                   ! Non-diagonal lambda-opq doubling, J-independent part 
                   !
                   do ild = 1,Nlambdaopq
                     !
                     field => lambdaopq(ild)
                     !
                     ! 1. <Sigma,Omega,Lambda|Lambda-O|Sigma+/-2,Omega,-Lambda>
                     if (lambdaopq(ild)%istate==istate.and.lambdaopq(ild)%jstate==jstate.and.istate==jstate.and.&
                         abs(ilambda)==1.and.(ilambda-jlambda)==nint(sigmaj-sigmai).and.abs(nint(sigmaj-sigmai))==2 &
                                    .and.(ilambda==-jlambda).and.nint(spini-spinj)==0.and.nint(omegai-omegaj)==0) then
                        !
                        f_s2 = sigmai-sigmaj
                        f_s1 = sign(1.0_rk,f_s2)
                        !
                        f_t = sqrt( spini*(spini+1.0_rk)-(sigmaj     )*(sigmaj+f_s1) )*&
                              sqrt( spini*(spini+1.0_rk)-(sigmaj+f_s1)*(sigmaj+f_s2) )
                        !
                        f_lo = field%gridvalue(igrid)*f_t*sc
                        !
                        omegamat(i,j) = omegamat(i,j) + f_lo*0.5_rk
                        omegamat(j,i) = omegamat(i,j)
                        !
                     endif
                     ! 
                   enddo
                   !
               enddo  ! j
             enddo  ! i
             !
             !Solve the Sigma-Lambda-State Hamiltonian at each grid point 
             !
             call lapack_syev(omegamat,omegaenergy)
             !
             N_i = Nlambdasigmas
             !
             omega_grid(iomega)%energy(1:N_i,igrid) = omegaenergy/sc
             omega_grid(iomega)%vector(1:N_i,1:N_i,igrid) = omegamat
             omega_grid(iomega)%Nstates = Nlambdasigmas
             !
             deallocate(omegamat,omegaenergy)
             !
             if (igrid==1) then 
                vect(iomega,1:N_i,1:N_i) = omega_grid(iomega)%vector(1:N_i,1:N_i,igrid)
             else
               !
               mat_1(1:N_i,1:N_i) = &
                 matmul(transpose(vect(iomega,1:N_i,1:N_i)),omega_grid(iomega)%vector(1:N_i,1:N_i,igrid))
               !
               ! scalar product of vect with the previous step
               do i = 1,N_i
                 !
                 imaxcontr = maxloc(mat_1(1:N_i,i)**2,dim=1,mask=mat_1(1:N_i,i)**2.ge.small_)
                 !
                 if (mat_1(imaxcontr,i)<-small_) then 
                     omega_grid(iomega)%vector(1:N_i,i,igrid) = -omega_grid(iomega)%vector(1:N_i,i,igrid)
                 endif
                 !
                 vect(iomega,1:N_i,1:N_i) = omega_grid(iomega)%vector(1:N_i,1:N_i,igrid)
                 !
               enddo
               !
             endif
             !
          enddo
          !
          ! Now we can use the eigenfunctions to  unitary trnsform different objects to the Omega representation
          !
          ! Compute the L+ matrix elements in the primitive Lambda-Sigma representation
          !
          if (NLPlus_omega/=0) then 
            !
            do iomega=1,Nomegas
              !
              omegai = Omega_grid(iomega)%omega
              N_i = Omega_grid(iomega)%Nstates
              !
              do jomega=1,Nomegas
                  !
                  omegaj = Omega_grid(jomega)%omega
                  N_j = Omega_grid(jomega)%Nstates
                  !
                  if (all(iLplus_omega(iomega,jomega,1:N_i,1:N_j)==0)) cycle
                  !
                  L_LambdaSigma = 0
                  !
                  do i = 1,N_i
                    !
                    istate  = Omega_grid(iomega)%basis(i)%istate
                    ilambda = Omega_grid(iomega)%basis(i)%ilambda
                    spini   = Omega_grid(iomega)%basis(i)%spin
                    sigmai   = Omega_grid(iomega)%basis(i)%sigma
                    !
                    do j = 1,N_j
                      !
                      jstate  = Omega_grid(jomega)%basis(j)%istate
                      jlambda = Omega_grid(jomega)%basis(j)%ilambda
                      spinj   = Omega_grid(jomega)%basis(j)%spin
                      sigmaj  = Omega_grid(jomega)%basis(j)%sigma
                      !
                      ! Lx is always diagonal in sigma and we also restrict to L+ (omegai = omegaj+1)
                      !
                      if( nint(sigmaj-sigmai)/=0.or.(omegai<omegaj) ) cycle
                      !
                      ! Lplus/Lminus
                      !
                      do ilxly =1,Nlxly
                        !
                        field => lxly(ilxly)
                        !
                        ! Also check that L+ is consistent with the selection rules
                        !
                        if ( field%istate==field%jstate .or.abs(field%lambda-field%lambdaj)/=1 ) then
                           !
                           write(out,"('The quantum numbers of the L+/Lx field ',2i3,' are inconsistent" // &
                                           " with L+selection rules: ')") field%istate,field%jstate
                           write(out,"('Delta Lamda = +/-1')")
                           stop "Lx/L+ input is inconsistent with selection rules"
                           !
                        endif
                        !
                        ! the field entry in the input gives only one combination of the quantum numbers for
                        ! the matrix element <State,Lambda,Spin|F|State',Lambda',Spin'>
                        ! LxLy  4 6 ;  lambda  0 1 ; spin   1.0 1.0
                        ! we should consider also a permutation <State',Lambda',Spin'|F|State,Lambda,Spin> if this makes a change.
                        !
                        do ipermute  = 0,1
                          !
                          if (ipermute==0) then
                            !
                            istate_ = field%istate ; ilambda_ = field%lambda  ; spini_ = field%spini
                            jstate_ = field%jstate ; jlambda_ = field%lambdaj ; spinj_ = field%spinj
                            !
                          else  ! permute
                            !
                            jstate_ = field%istate ; jlambda_ = field%lambda  ; spinj_ = field%spini
                            istate_ = field%jstate ; ilambda_ = field%lambdaj ; spini_ = field%spinj
                            !
                          endif
                          !
                          ! check if we at the right electronic states
                          if( istate/=istate_.or.jstate/=jstate_ ) cycle
                          !
                          ! We should also take into account that Lambda can change sign (only Lambda>0 is given in input)
                          ! In order to recover other combinations we apply the symmetry transformation
                          ! laboratory fixed inversion which is equivalent to the sigmav operation 
                          !                    (sigmav= 0 correspond to the unitary transformation)
                          do isigmav = 0,1
                            !
                            ! the permutation is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
                            ! avoid the double counting.
                            if( isigmav==1.and. abs( field%lambda ) + abs( field%lambdaj )==0 ) cycle
                       
                            ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                            ilambda_ = ilambda_*(-1)**isigmav
                            jlambda_ = jlambda_*(-1)**isigmav
                            !
                            ! proceed only if the quantum numbers of the field equal to the corresponding <i| and |j> quantum numbers:
                            if (istate/=istate_.or.jstate_/=jstate.or.ilambda_/=ilambda.or.jlambda_/=jlambda) cycle
                            !
                            ! check the selection rule Delta Lambda = +/1
                            if (abs(ilambda-jlambda)/=1) cycle
                            !
                            f_grid  = field%gridvalue(igrid)
                            !
                            ! <Lx> and <Ly> don't depend on Sigma
                            !
                            f_s = real(ilambda-jlambda,rk)
                            !
                            f_t = f_grid
                            !
                            ! the result of the symmetry transformation:
                            if (isigmav==1) then
                              !
                              ! we assume that
                              ! sigmav <Lamba|L+|Lambda'> => <-Lamba|L-|-Lambda'> == <Lamba|L+|Lambda'>(-1)^(Lamba+Lambda')
                              ! and <Lamba|L+|Lambda'> is an unique quantity given in the input
                              ! however we don't apply the sigmav transformation to sigma or omega
                              ! since we only need to know how <Lamba|L+/-|Lambda'> transforms in order to relate it to the
                              ! value given in input.
                              !
                              !itau = ilambda-jlambda+nint(spini-sigmai)+nint(spinj-sigmaj)!-nint(omegai+omegaj)
                              !
                              ! we try to remove also lambda from the sigmav transformation!!
                              !
                              itau = 0 
                              !
                              if (ilambda_==0.and.poten(istate)%parity%pm==-1) itau = itau+1
                              if (jlambda_==0.and.poten(jstate)%parity%pm==-1) itau = itau+1
                              !
                              f_t = f_t*(-1.0_rk)**(itau)
                              !
                            endif
                            !
                            L_LambdaSigma(i,j) = f_t
                            !
                          enddo
                          !
                        enddo
                        !
                      enddo
                      !
                    enddo
                    !
                  enddo
                  !
                  mat_1(1:N_i,1:N_j) =matmul(L_LambdaSigma(1:N_i,1:N_j),omega_grid(jomega)%vector(1:N_j,1:N_j,igrid))
                  mat_2(1:N_i,1:N_j) =matmul(transpose(omega_grid(iomega)%vector(1:N_i,1:N_i,igrid)),mat_1(1:N_i,1:N_j))
                  !
                  do i = 1,N_i
                     do j = 1,N_j
                        !
                        iLplus_omega_ = iLplus_omega(iomega,jomega,i,j)
                        !
                        if (iLplus_omega_==0) cycle
                        !
                        l_omega_obj(iLplus_omega_)%gridvalue(igrid) = mat_2(i,j)
                        !
                     enddo
                  enddo
                  !
               enddo
            enddo
            !
          endif
          !
          ! Compute the S+ matrix elements in the primitive Lambda-Sigma representation
          !
          if (NSPlus_omega/=0) then 
            !
            do iomega=1,Nomegas
              !
              omegai = Omega_grid(iomega)%omega
              N_i = Omega_grid(iomega)%Nstates
              !
              do jomega=1,Nomegas
                  !
                  omegaj = Omega_grid(jomega)%omega
                  N_j = Omega_grid(jomega)%Nstates
                  !
                  if (all(iSplus_omega(iomega,jomega,1:N_i,1:N_j)==0)) cycle
                  !
                  L_LambdaSigma = 0 
                  !
                  do i = 1,N_i
                    !
                    istate  = Omega_grid(iomega)%basis(i)%istate
                    ilambda = Omega_grid(iomega)%basis(i)%ilambda
                    sigmai  = Omega_grid(iomega)%basis(i)%sigma
                    spini   = Omega_grid(iomega)%basis(i)%spin
                    !
                    do j = 1,N_j
                      !
                      jstate  = Omega_grid(jomega)%basis(j)%istate
                      jlambda = Omega_grid(jomega)%basis(j)%ilambda
                      sigmaj  = Omega_grid(jomega)%basis(j)%sigma
                      spinj   = Omega_grid(jomega)%basis(j)%spin
                      !
                      if (spini/=spini.or.(ilambda/=jlambda).or.nint(sigmai-sigmaj)/=1) cycle
                      !
                      do ispin = 1,Nspins
                         !
                         spini_ = real(ispin-1,rk)+spin_min
                         !
                         spinj = spini
                         !
                         L_LambdaSigma(i,j) = sqrt( (spini-sigmaj)*(spini+sigmaj+1.0_rk) )
                         !
                      enddo
                      !
                    enddo
                  enddo
                  !
                  mat_1(1:N_i,1:N_j) =matmul(L_LambdaSigma(1:N_i,1:N_j),omega_grid(jomega)%vector(1:N_j,1:N_j,igrid))
                  mat_2(1:N_i,1:N_j) =matmul(transpose(omega_grid(iomega)%vector(1:N_i,1:N_i,igrid)),mat_1(1:N_i,1:N_j))
                  !
                  do i = 1,N_i
                     do j = 1,N_j
                        !
                        iSplus_omega_ = iSplus_omega(iomega,jomega,i,j)
                        !
                        if (iSplus_omega_==0) cycle
                        !
                        S_omega_obj(iSplus_omega_)%gridvalue(igrid) = mat_2(i,j)
                        !                            
                     enddo
                  enddo
                  !
               enddo
            enddo
            !
          endif
          !
          ! Compute the SR (spin-rotation) matrix elements in the primitive Lambda-Sigma representation
          !
          if (NSR_omega/=0) then 
            !
            do iomega=1,Nomegas
              !
              omegai = Omega_grid(iomega)%omega
              N_i = Omega_grid(iomega)%Nstates
              !
              do jomega=1,Nomegas
                  !
                  omegaj = Omega_grid(jomega)%omega
                  N_j = Omega_grid(jomega)%Nstates
                  !
                  if (all(iSR_omega(iomega,jomega,1:N_i,1:N_j)==0)) cycle
                  !
                  L_LambdaSigma = 0
                  !
                  do i = 1,N_i
                    !
                    istate  = Omega_grid(iomega)%basis(i)%istate
                    ilambda = Omega_grid(iomega)%basis(i)%ilambda
                    spini   = Omega_grid(iomega)%basis(i)%spin
                    sigmai   = Omega_grid(iomega)%basis(i)%sigma
                    !
                    do j = 1,N_j
                      !
                      jstate  = Omega_grid(jomega)%basis(j)%istate
                      jlambda = Omega_grid(jomega)%basis(j)%ilambda
                      spinj   = Omega_grid(jomega)%basis(j)%spin
                      sigmaj  = Omega_grid(jomega)%basis(j)%sigma
                      !
                      ! Sx,Sy are always diagonal in Lambda and omegai = omegaj+/-1)
                      !
                      if( (ilambda-jlambda)/=0.or.abs(nint(omegai-omegaj))/=1 ) cycle
                      !
                      do iSR =1,NSR
                        !
                        field => spinrot(isr)
                        !
                        ! 1. <Sigma,Omega,Lambda|HSR|Sigma+/-1,Omega+/-1,Lambda>
                        if (spinrot(isr)%istate==istate.and.spinrot(isr)%jstate==jstate.and.istate==jstate.and.&
                            abs(nint(sigmaj-sigmai))==1.and.(ilambda==jlambda).and.nint(spini-spinj)==0) then
                           !
                           f_s = sigmaj-sigmai
                           !
                           f_t = sqrt( spini*(spini+1.0_rk)-sigmai*(sigmai+f_s) )
                           !
                           f_grid  = field%gridvalue(igrid)*f_t
                           !
                           L_LambdaSigma(i,j) = f_grid
                           !
                           !
                        endif
                        !
                      enddo
                      !
                    enddo
                    !
                  enddo
                  !
                  mat_1(1:N_i,1:N_j) =matmul(L_LambdaSigma(1:N_i,1:N_j),omega_grid(jomega)%vector(1:N_j,1:N_j,igrid))
                  mat_2(1:N_i,1:N_j) =matmul(transpose(omega_grid(iomega)%vector(1:N_i,1:N_i,igrid)),mat_1(1:N_i,1:N_j))
                  !
                  do i = 1,N_i
                     do j = 1,N_j
                        !
                        iSR_omega_ = iSR_omega(iomega,jomega,i,j)
                        !
                        if (iSR_omega_==0) cycle
                        !
                        sr_omega_obj(iSR_omega_)%gridvalue(igrid) = mat_2(i,j)
                        !
                     enddo
                  enddo
                  !
               enddo
            enddo
            !
          endif
          !
          !
          ! Compute the BobRot (BOB) matrix elements in the primitive Lambda-Sigma representation
          !
          if (NBOB_omega/=0) then 
            !
            do iomega=1,Nomegas
              !
              omegai = Omega_grid(iomega)%omega
              N_i = Omega_grid(iomega)%Nstates
              !
              if (all(iBOB_omega(iomega,1:N_i,1:N_i)==0)) cycle
              !
              L_LambdaSigma = 0
              !
              do i = 1,N_i
                !
                istate  = Omega_grid(iomega)%basis(i)%istate
                ilambda = Omega_grid(iomega)%basis(i)%ilambda
                spini   = Omega_grid(iomega)%basis(i)%spin
                sigmai  = Omega_grid(iomega)%basis(i)%sigma
                !
                do ibob =1,Nbobrot
                  !
                  field => bobrot(ibob)
                  !
                  if (field%istate==istate.and.field%jstate==istate) then
                     !
                     f_grid  = field%gridvalue(igrid)*f_t
                     !
                     L_LambdaSigma(i,i) = f_grid
                     !
                     !
                  endif
                  !
                enddo
                !
              enddo
              !
              mat_1(1:N_i,1:N_i) = matmul(L_LambdaSigma(1:N_i,1:N_i),omega_grid(iomega)%vector(1:N_i,1:N_i,igrid))
              mat_2(1:N_i,1:N_i) = matmul(transpose(omega_grid(iomega)%vector(1:N_i,1:N_i,igrid)),mat_1(1:N_i,1:N_i))
              !
              do i = 1,N_i
                 do j = 1,N_i
                    !
                    iBob_omega_ = iBob_omega(iomega,i,j)
                    !
                    if (iBob_omega_==0) cycle
                    !
                    BOB_omega_obj(iBob_omega_)%gridvalue(igrid) = mat_2(i,j)
                    !
                 enddo
              enddo
              !
            enddo
            !
          endif
          !
          ! Compute the BobRot (BRot) matrix elements in the primitive Lambda-Sigma representation
          !
          if (NBRot_omega/=0) then 
            !
            do iomega=1,Nomegas
              !
              omegai = Omega_grid(iomega)%omega
              N_i = Omega_grid(iomega)%Nstates
              !
              if (all(iBRot_omega(iomega,1:N_i,1:N_i)==0)) cycle
              !
              L_LambdaSigma = 0
              !
              do i = 1,N_i
                !
                L_LambdaSigma(i,i) = 1.0_rk
                !
              enddo
              !
              mat_1(1:N_i,1:N_i) = matmul(L_LambdaSigma(1:N_i,1:N_i),omega_grid(iomega)%vector(1:N_i,1:N_i,igrid))
              mat_2(1:N_i,1:N_i) = matmul(transpose(omega_grid(iomega)%vector(1:N_i,1:N_i,igrid)),mat_1(1:N_i,1:N_i))
              !
              do i = 1,N_i
                 do j = 1,N_i
                    !
                    iBRot_omega_ = iBRot_omega(iomega,i,j)
                    !
                    if (iBRot_omega_==0) cycle
                    !
                    BRot_omega_obj(iBRot_omega_)%gridvalue(igrid) = mat_2(i,j)
                    !
                 enddo
              enddo
              !
            enddo
            !
          endif
          !
          !
          ! Compute the p2q (lambda-doubling) matrix elements in the primitive Lambda-Sigma representation
          !
          if (Np2q_omega/=0) then 
            !
            do iomega=1,Nomegas
              !
              omegai = Omega_grid(iomega)%omega
              N_i = Omega_grid(iomega)%Nstates
              !
              do jomega=1,Nomegas
                  !
                  omegaj = Omega_grid(jomega)%omega
                  N_j = Omega_grid(jomega)%Nstates
                  !
                  if (all(ip2q_omega(iomega,jomega,1:N_i,1:N_j)==0)) cycle
                  !
                  L_LambdaSigma = 0
                  !
                  do i = 1,N_i
                    !
                    istate   = Omega_grid(iomega)%basis(i)%istate
                    ilambda  = Omega_grid(iomega)%basis(i)%ilambda
                    spini    = Omega_grid(iomega)%basis(i)%spin
                    sigmai   = Omega_grid(iomega)%basis(i)%sigma
                    !
                    do j = 1,N_j
                      !
                      jstate  = Omega_grid(jomega)%basis(j)%istate
                      jlambda = Omega_grid(jomega)%basis(j)%ilambda
                      spinj   = Omega_grid(jomega)%basis(j)%spin
                      sigmaj  = Omega_grid(jomega)%basis(j)%sigma
                      !
                      if( abs(nint(omegai-omegaj))/=1.or.abs(nint(sigmai-sigmaj))/=1 ) cycle
                      if( abs(ilambda)/=1.or.abs(jlambda)/=1.or.abs(ilambda-jlambda)/=2 ) cycle
                      if (istate/=jstate.or.nint(spini-spinj)/=0.or.nint(sigmaj-sigmai)/=nint(omegai-omegaj) ) cycle
                      !
                      do ip2q =1,Nlambdap2q
                        !
                        field => lambdap2q(ip2q)
                        !
                        ! <Sigma+/-1,Omega-/+1,Lambda=-/+1|Hp2q|Sigma,Omega,-Lambda>
                        if (field%istate==istate.and.field%jstate==jstate.and.abs(field%lambda)==1) then
                           !
                           f_s = sigmai-sigmaj
                           !
                           f_t = sqrt( spini*(spini+1.0_rk)-sigmaj*(sigmaj+f_s) )
                           !
                           f_grid = field%gridvalue(igrid)*f_t
                           !
                           L_LambdaSigma(i,j) = f_grid
                           !
                           !
                        endif
                        !
                      enddo
                      !
                    enddo
                    !
                  enddo
                  !
                  mat_1(1:N_i,1:N_j) =matmul(L_LambdaSigma(1:N_i,1:N_j),omega_grid(jomega)%vector(1:N_j,1:N_j,igrid))
                  mat_2(1:N_i,1:N_j) =matmul(transpose(omega_grid(iomega)%vector(1:N_i,1:N_i,igrid)),mat_1(1:N_i,1:N_j))
                  !
                  do i = 1,N_i
                     do j = 1,N_j
                        !
                        ip2q_omega_ = ip2q_omega(iomega,jomega,i,j)
                        !
                        if (ip2q_omega_==0) cycle
                        !
                        p2q_omega_obj(ip2q_omega_)%gridvalue(igrid) = mat_2(i,j)
                        !
                     enddo
                  enddo
                  !
               enddo
            enddo
            !
          endif
          !
          ! Compute the q (lambda-doubling) matrix elements in the primitive Lambda-Sigma representation
          !
          if (Nq_omega/=0) then 
            !
            do iomega=1,Nomegas
              !
              omegai = Omega_grid(iomega)%omega
              N_i = Omega_grid(iomega)%Nstates
              !
              do jomega=1,Nomegas
                  !
                  omegaj = Omega_grid(jomega)%omega
                  N_j = Omega_grid(jomega)%Nstates
                  !
                  if (all(iq_omega(iomega,jomega,1:N_i,1:N_j)==0)) cycle
                  !
                  L_LambdaSigma = 0
                  !
                  do i = 1,N_i
                    !
                    istate  = Omega_grid(iomega)%basis(i)%istate
                    ilambda = Omega_grid(iomega)%basis(i)%ilambda
                    spini   = Omega_grid(iomega)%basis(i)%spin
                    sigmai   = Omega_grid(iomega)%basis(i)%sigma
                    !
                    do j = 1,N_j
                      !
                      jstate  = Omega_grid(jomega)%basis(j)%istate
                      jlambda = Omega_grid(jomega)%basis(j)%ilambda
                      spinj   = Omega_grid(jomega)%basis(j)%spin
                      sigmaj  = Omega_grid(jomega)%basis(j)%sigma
                      !
                      if( abs(nint(omegai-omegaj))/=2.or.nint(sigmai-sigmaj)/=0 ) cycle
                      if( abs(ilambda)/=1.or.abs(jlambda)/=1.or.abs(ilambda-jlambda)/=2 ) cycle
                      if (istate/=jstate.or.(ilambda-jlambda)/=nint(omegai-omegaj).or.nint(spini-spinj)/=0) cycle
                      !
                      do iq =1,Nlambdaq
                        !
                        field => lambdaq(iq)
                        !
                        ! 1. <Sigma,Omega,Lambda|Lambda-O|Sigma+/-2,Omega,-Lambda>
                        if (field%istate==istate.and.field%jstate==jstate.and.abs(field%lambda)==1) then
                           !
                           f_grid = field%gridvalue(igrid)
                           !
                           L_LambdaSigma(i,j) = f_grid
                           !
                           !
                        endif
                        !
                      enddo
                      !
                    enddo
                    !
                  enddo
                  !
                  mat_1(1:N_i,1:N_j) =matmul(L_LambdaSigma(1:N_i,1:N_j),omega_grid(jomega)%vector(1:N_j,1:N_j,igrid))
                  mat_2(1:N_i,1:N_j) =matmul(transpose(omega_grid(iomega)%vector(1:N_i,1:N_i,igrid)),mat_1(1:N_i,1:N_j))
                  !
                  do i = 1,N_i
                     do j = 1,N_j
                        !
                        iq_omega_ = iq_omega(iomega,jomega,i,j)
                        !
                        if (iq_omega_==0) cycle
                        !
                        q_omega_obj(iq_omega_)%gridvalue(igrid) = mat_2(i,j)
                        !
                     enddo
                  enddo
                  !
               enddo
            enddo
            !
          endif
          !
       enddo
       !
       !
       allocate(vibmat(ngrid,ngrid),stat=alloc)
       call ArrayStart('vibmat-omega',alloc,size(vibmat),kind(vibmat))
       allocate(kinmat(ngrid,ngrid),stat=alloc)
       call ArrayStart('kinmat-omega',alloc,size(kinmat),kind(kinmat))
       !
       !
       ! Kinetic energy part 
       !
       call kinetic_energy_grid_points(ngrid,kinmat,vibmat)
       !
       ! For each grid point diagonalise the Sigma-Lambda PECs + SOs and transform to the Omega-represenation
       do igrid =1, ngrid
          do jgrid =1,ngrid
            !
            ! Compute the Kinetic (Kin) matrix elements in the primitive Lambda-Sigma representation
            !
            do iomega=1,Nomegas
              !
              omegai = Omega_grid(iomega)%omega
              N_i = Omega_grid(iomega)%Nstates
              !
              L_LambdaSigma = 0
              !
              do i = 1,N_i
                !
                istate  = Omega_grid(iomega)%basis(i)%istate
                ilambda = Omega_grid(iomega)%basis(i)%ilambda
                spini   = Omega_grid(iomega)%basis(i)%spin
                sigmai  = Omega_grid(iomega)%basis(i)%sigma
                !
                L_LambdaSigma(i,i) = vibmat(igrid,jgrid)
                !
              enddo
              !
              mat_1(1:N_i,1:N_i) = matmul(L_LambdaSigma(1:N_i,1:N_i),omega_grid(iomega)%vector(1:N_i,1:N_i,igrid))
              mat_2(1:N_i,1:N_i) = matmul(transpose(omega_grid(iomega)%vector(1:N_i,1:N_i,jgrid)),mat_1(1:N_i,1:N_i))
              !
              do i = 1,N_i
                 do j = 1,N_i
                    !
                    iKin_omega_ = iKin_omega(iomega,i,j)
                    !
                    Kin_omega_obj(iKin_omega_)%matelem(igrid,jgrid) = mat_2(i,j)
                    !
                 enddo
              enddo
              !
            enddo
            !
         enddo
       enddo
       !
       deallocate(vibmat)
       call ArrayStop('vibmat-omega')
       deallocate(kinmat)
       call ArrayStop('kinmat-omega')
       !
       deallocate(L_LambdaSigma)
       call ArrayStop('L_LambdaSigma')
       !
       deallocate(mat_1,mat_2)
       call ArrayStop('mat_1')
       call ArrayStop('mat_2')
       !
       deallocate(vect)
       call ArrayStop('vect')
       !
       deallocate(imax_contr)
       call ArrayStop('imax_contr')
       !
       deallocate(iomega_state)
       call ArrayStop('iomega_state')
       !
       do iomega=1,Nomegas
          !
          omegai = Omega_grid(iomega)%omega
          !
          N_i = omega_grid(iomega)%Nstates
          !
          do i = 1,N_i
            !
            ieq =minloc(omega_grid(iomega)%energy(i,:),dim=1)
            !
            imaxcontr = maxloc(omega_grid(iomega)%vector(1:N_i,i,ieq)**2,&
                        dim=1,mask=omega_grid(iomega)%vector(1:N_i,i,ieq)**2.ge.small_)
            !
            Omega_grid(iomega)%qn(i)%istate  = Omega_grid(iomega)%basis(imaxcontr)%istate
            Omega_grid(iomega)%qn(i)%name    = trim(Omega_grid(iomega)%basis(imaxcontr)%name)
            Omega_grid(iomega)%qn(i)%sigma   = Omega_grid(iomega)%basis(imaxcontr)%sigma
            Omega_grid(iomega)%qn(i)%ilambda = Omega_grid(iomega)%basis(imaxcontr)%ilambda
            Omega_grid(iomega)%qn(i)%spin    = Omega_grid(iomega)%basis(imaxcontr)%spin
            Omega_grid(iomega)%qn(i)%omega   = omegai
            Omega_grid(iomega)%qn(i)%ilevel  = Omega_grid(iomega)%basis(imaxcontr)%ilevel
            !
          enddo
          !
       enddo
       !
       !
 end subroutine Transfrorm_Sigma_Lambda_to_Omega_representation


 subroutine print_fileds_in_Omega_representation(iverbose,Nomega_states)
      !
      implicit none
      !
      integer(ik),intent(in) :: iverbose,Nomega_states
      !
      integer(ik) :: i,j,igrid,iterm,iomega,jomega,ngrid
       !
       ngrid = grid%npoints
       !
       ! Print PECs 
       !
       if (iverbose>=3) then 
          !
          write(out,'(/7x,a)') "PECS in the Omega representation"
          write(out,'(7x,a)') "#  #    Omega Lambda  Sigma  State"
          !
          do iomega=1,Nomegas
            do i=1,Omega_grid(iomega)%Nstates
             write(out,'(4x,i4,1x,i2,1x,f8.1,1x,1x,i3,1x,f8.1,1x,a)') iomega,i,&
                       Omega_grid(iomega)%omega,Omega_grid(iomega)%qn(i)%ilambda,&
                       Omega_grid(iomega)%qn(i)%sigma,trim(Omega_grid(iomega)%qn(i)%name)
            enddo
          enddo
          !
          write(my_fmt, '(A,I0,A)') '("            r(Ang)",', Nomega_states, '(i22))'
          write(out,my_fmt) (i,i=1,Nomega_states)
          !
          write(my_fmt, '(A,I0,A)') '(f18.8,', Nomega_states, '(es22.8))'
          do igrid=1,ngrid
             write(out,my_fmt) &
               r(igrid),((omega_grid(iomega)%energy(j,igrid),j=1,Omega_grid(iomega)%Nstates),iomega=1,Nomegas)
          enddo
          !
          ! Print EAMC L+ in Omega
          !
          if (NLplus_omega>0) then
            !
            write(out,'(/9x,a)') "L+ in the Omega representation"
            write(out,'(9x,a)') "#    Omega State Lambda Sigma <->Omega State Lambda Sigma"
            !
            do iterm=1,NLplus_omega
               !
               i = L_omega_obj(iterm)%ilevel
               j = L_omega_obj(iterm)%jlevel
               iomega = L_omega_obj(iterm)%iomega
               jomega = L_omega_obj(iterm)%jomega
               !
               write(out,'(7x,i4,1x,f8.1,1x,i2,1x,i4,1x,f8.1,3x,f8.1,1x,i2,1x,i4,1x,f8.1)') iterm,&
                 L_omega_obj(iterm)%omegai,&
                 Omega_grid(iomega)%qn(i)%istate,Omega_grid(iomega)%qn(i)%ilambda,Omega_grid(iomega)%qn(i)%sigma ,&
                 L_omega_obj(iterm)%omegaj,&
                 Omega_grid(jomega)%qn(j)%istate,Omega_grid(jomega)%qn(j)%ilambda,Omega_grid(jomega)%qn(j)%sigma
            enddo
            !
            write(my_fmt, '(A,I0,A)') '("            r(Ang)",', NLplus_omega, '(i22))'
            write(out,my_fmt) (i,i=1,NLplus_omega)
            !
            write(my_fmt, '(A,I0,A)') '(f18.8,', NLplus_omega, '(es22.8))'
            do igrid=1,ngrid
               write(out,my_fmt) r(igrid),(L_omega_obj(iterm)%gridvalue(igrid),iterm=1,NLplus_omega)
            enddo
            !
          endif
          !
          ! Print Spin S+ in Omega
          !
          if (NSplus_omega>0) then
            !
            write(out,'(/10x,a)') "S+ in the Omega representation"
            write(out,'(10x,a)') "#    Omega State Lambda Sigma <->Omega State Lambda Sigma"
         
            do iterm=1,NSplus_omega
               !
               i = S_omega_obj(iterm)%ilevel
               j = S_omega_obj(iterm)%jlevel
               iomega = S_omega_obj(iterm)%iomega
               jomega = S_omega_obj(iterm)%jomega
               !
               write(out,'(7x,i4,1x,f8.1,1x,i2,1x,i4,1x,f8.1,3x,f8.1,1x,i2,1x,i4,1x,f8.1)') iterm,&
                 S_omega_obj(iterm)%omegai,&
                 Omega_grid(iomega)%qn(i)%istate,Omega_grid(iomega)%qn(i)%ilambda,Omega_grid(iomega)%qn(i)%sigma ,&
                 S_omega_obj(iterm)%omegaj,&
                 Omega_grid(jomega)%qn(j)%istate,Omega_grid(jomega)%qn(j)%ilambda,Omega_grid(jomega)%qn(j)%sigma
            enddo
            !
            write(my_fmt, '(A,I0,A)') '("            r(Ang)",', NSplus_omega, '(i22))'
            write(out,my_fmt) (i,i=1,NSplus_omega)
            !
            write(my_fmt, '(A,I0,A)') '(f18.8,', NSplus_omega, '(es22.8))'
            do igrid=1,ngrid
               write(out,my_fmt) r(igrid),(S_omega_obj(iterm)%gridvalue(igrid),iterm=1,NSplus_omega)
            enddo
            !
          endif
          !
          ! Print spint-rotation in Omega
          !
          if (NSR_omega>0) then
            !
            write(out,'(/10x,a)') "Spin-rotation in the Omega representation"
            write(out,'(10x,a)') "#    Omega State Lambda Sigma <->Omega State Lambda Sigma"
         
            do iterm=1,NSR_omega
               !
               i = sr_omega_obj(iterm)%ilevel
               j = sr_omega_obj(iterm)%jlevel
               iomega = sr_omega_obj(iterm)%iomega
               jomega = sr_omega_obj(iterm)%jomega
               !
               write(out,'(7x,i4,1x,f8.1,1x,i2,1x,i4,1x,f8.1,3x,f8.1,1x,i2,1x,i4,1x,f8.1)') iterm,&
                 sr_omega_obj(iterm)%omegai,&
                 Omega_grid(iomega)%qn(i)%istate,Omega_grid(iomega)%qn(i)%ilambda,Omega_grid(iomega)%qn(i)%sigma,&
                 sr_omega_obj(iterm)%omegaj,&
                 Omega_grid(jomega)%qn(j)%istate,Omega_grid(jomega)%qn(j)%ilambda,Omega_grid(jomega)%qn(j)%sigma
            enddo
            !
            write(my_fmt, '(A,I0,A)') '("            r(Ang)",', NSR_omega, '(i22))'
            write(out,my_fmt) (i,i=1,NSR_omega)
            !
            write(my_fmt, '(A,I0,A)') '(f18.8,', NSR_omega, '(es22.8))'
            do igrid=1,ngrid
               write(out,my_fmt) r(igrid),(sr_omega_obj(iterm)%gridvalue(igrid),iterm=1,NSR_omega)
            enddo
            !
          endif
          !
          ! Print Bob-rot in Omega
          !
          if (NBob_omega>0) then
            !
            write(out,'(/10x,a)') "Bob-rotation in the Omega representation"
            write(out,'(10x,a)') "#    Omega State Lambda Sigma State Lambda Sigma"
         
            do iterm=1,NBob_omega
               !
               i = bob_omega_obj(iterm)%ilevel
               j = bob_omega_obj(iterm)%jlevel
               iomega = bob_omega_obj(iterm)%iomega
               !
               write(out,'(7x,i4,1x,f8.1,1x,i2,1x,i4,1x,f8.1,3x,i2,1x,i4,1x,f8.1,1x)') iterm,&
                 bob_omega_obj(iterm)%omegai,&
                 Omega_grid(iomega)%qn(i)%istate,Omega_grid(iomega)%qn(i)%ilambda,Omega_grid(iomega)%qn(i)%sigma ,&
                 Omega_grid(jomega)%qn(j)%istate,Omega_grid(jomega)%qn(j)%ilambda,Omega_grid(iomega)%qn(j)%sigma
            enddo
            !
            write(my_fmt, '(A,I0,A)') '("            r(Ang)",', NBob_omega, '(i22))'
            write(out,my_fmt) (i,i=1,NBob_omega)
            !
            write(my_fmt, '(A,I0,A)') '(f18.8,', NBob_omega, '(es22.8))'
            do igrid=1,ngrid
               write(out,my_fmt) r(igrid),(bob_omega_obj(iterm)%gridvalue(igrid),iterm=1,NBob_omega)
            enddo
            !
          endif
          !
          !
          ! Print lambda-doubling-p2q in Omega
          !
          if (Np2q_omega>0) then
            !
            write(out,'(/10x,a)') "Lambda-doubling for Pi p2q in the Omega representation"
            write(out,'(10x,a)') "#    Omega State Lambda Sigma <->Omega State Lambda Sigma"
            !
            do iterm=1,Np2q_omega
               !
               i = p2q_omega_obj(iterm)%ilevel
               j = p2q_omega_obj(iterm)%jlevel
               iomega = p2q_omega_obj(iterm)%iomega
               jomega = p2q_omega_obj(iterm)%jomega
               !
               write(out,'(7x,i4,1x,f8.1,1x,i2,1x,i4,1x,f8.1,3x,f8.1,1x,i2,1x,i4,1x,f8.1)') iterm,&
                 p2q_omega_obj(iterm)%omegai,&
                 Omega_grid(iomega)%qn(i)%istate,Omega_grid(iomega)%qn(i)%ilambda,Omega_grid(iomega)%qn(i)%sigma,&
                 p2q_omega_obj(iterm)%omegaj,&
                 Omega_grid(jomega)%qn(j)%istate,Omega_grid(jomega)%qn(j)%ilambda,Omega_grid(jomega)%qn(j)%sigma
            enddo
            !
            write(my_fmt, '(A,I0,A)') '("            r(Ang)",', Np2q_omega, '(i22))'
            write(out,my_fmt) (i,i=1,Np2q_omega)
            !
            write(my_fmt, '(A,I0,A)') '(f18.8,', Np2q_omega, '(es22.8))'
            do igrid=1,ngrid
               write(out,my_fmt) r(igrid),(p2q_omega_obj(iterm)%gridvalue(igrid),iterm=1,Np2q_omega)
            enddo
            !
          endif
          !
          ! Print lambda-doubling-q in Omega
          !
          if (Nq_omega>0) then
            !
            write(out,'(/10x,a)') "Lambda-doubling for Pi q in the Omega representation"
            write(out,'(10x,a)') "#    Omega State Lambda Sigma <->Omega State Lambda Sigma"
         
            do iterm=1,Nq_omega
               !
               i = q_omega_obj(iterm)%ilevel
               j = q_omega_obj(iterm)%jlevel
               iomega = q_omega_obj(iterm)%iomega
               jomega = q_omega_obj(iterm)%jomega
               !
               write(out,'(7x,i4,1x,f8.1,1x,i2,1x,i4,1x,f8.1,3x,f8.1,1x,i2,1x,i4,1x,f8.1)') iterm,&
                 q_omega_obj(iterm)%omegai,&
                 Omega_grid(iomega)%qn(i)%istate,Omega_grid(iomega)%qn(i)%ilambda,Omega_grid(iomega)%qn(i)%sigma,&
                 q_omega_obj(iterm)%omegaj,&
                 Omega_grid(jomega)%qn(j)%istate,Omega_grid(jomega)%qn(j)%ilambda,Omega_grid(jomega)%qn(j)%sigma
            enddo
            !
            write(my_fmt, '(A,I0,A)') '("            r(Ang)",', Nq_omega, '(i22))'
            write(out,my_fmt) (i,i=1,Nq_omega)
            !
            write(my_fmt, '(A,I0,A)') '(f18.8,', Nq_omega, '(es22.8))'
            do igrid=1,ngrid
               write(out,my_fmt) r(igrid),(q_omega_obj(iterm)%gridvalue(igrid),iterm=1,Nq_omega)
            enddo
            !
          endif
          !
       endif
       !
       if (iverbose>=3) write(out,'("...done!")')
       !
  end subroutine print_fileds_in_Omega_representation



  subroutine Solve_vibrational_problem_for_Omega_states(iverbose,ngrid,Nomega_states,sc,totalroots,&
                               icontrvib,contrenergy,contracted)
      !
      use lapack
      !
      implicit none
      !
      integer(ik),intent(in) :: iverbose,Nomega_states,ngrid
      !type(Omega_gridT),intent(in) :: omega_grid(Nomegas)
      real(rk),intent(in)       :: sc
      integer(ik),intent(out)   :: totalroots
      type(contract_solT),intent(inout) :: contracted(Nomegas)
      type(quantaT),intent(out) :: icontrvib(ngrid*Nomega_states)
      real(rk),intent(out)      :: contrenergy(ngrid*Nomega_states)
      !
      integer(ik) :: alloc,iomega,jomega,ilevel,jlevel,Nlambdasigmas,igrid,Nroots,i,j,istate,u1,iKin,ilevel_,jlevel_
      integer(ik) :: Ndimen,ielem,imaxcontr
      real(rk)    :: omega,b_rot,epot,zpe,energy_
      real(rk)    :: psipsi_t,omegai_,energy_j
      real(rk),allocatable    :: vibmat(:,:),vibener(:),Kinmat(:,:)
      character(len=1)        :: rng,jobz
      real(rk)                :: vrange(2)
      integer(ik)             :: irange(2)
      type(quantaT)           :: icontrvib_
      type(fieldT),pointer    :: field
       !
       !
       if (trim(solution_method)=="LOBATTO".or.grid%nsub == 6) then 
         write(out,"('The LOBATTO is currently not working for the OMEGA representation')")
         stop 'The LOBATTO is currently not working for the OMEGA representation'
       endif
       !
       b_rot = aston/amass
       !
       totalroots = 0
       !
       allocate(kinmat(Ngrid,Ngrid),stat=alloc)
       call ArrayStart('kinmat',alloc,size(kinmat),kind(kinmat))
       !
       do iomega=1,Nomegas
          !
          if (iverbose>=4) call TimerStart('Build vibrational Hamiltonian')
          !
          Nlambdasigmas = omega_grid(iomega)%Nstates
          Ndimen = Nlambdasigmas*ngrid
          !
          allocate(vibmat(Ndimen,Ndimen),vibener(Ndimen),stat=alloc)
          call ArrayStart('vibmat',alloc,size(vibmat),kind(vibmat))
          call ArrayStart('vibener',alloc,size(vibener),kind(vibmat))
          !
          vibmat = 0
          !
          do ilevel = 1,Nlambdasigmas
            !
            if (iverbose>=6) write(out,'("ilevel = ",i0)') ilevel
            !
            omega  = Omega_grid(iomega)%basis(ilevel)%omega
            !
            kinmat = 0 
            !
            do jlevel = 1,Nlambdasigmas
               !
               ! Kinetic energy values 
               !
               do iKin = 1,NKin_omega
                  !
                  field => Kin_omega_obj(iKin)
                  !
                  omegai_ = field%omegai
                  ilevel_ = field%ilevel
                  jlevel_ = field%jlevel
                  !
                  if (ilevel_/=ilevel.or.jlevel_/=ilevel.or.nint(omegai_-omega)/=0) cycle
                  !
                  kinmat = Kin_omega_obj(iKin)%matelem
                  !
                  exit
                  !
               enddo
               !
               vibmat((ilevel-1)*ngrid+1:ilevel*ngrid,(jlevel-1)*ngrid+1:jlevel*ngrid) = kinmat
               !
            enddo
            !
          enddo
          !
          do ilevel = 1,Nlambdasigmas
            !
            omega  = Omega_grid(iomega)%basis(ilevel)%omega
            !
            if (iverbose>=6) write(out,'("ilevel = ",i0)') ilevel
            !
            !$omp parallel do private(igrid,epot,ielem) shared(vibmat) schedule(guided)
            do igrid =1, ngrid
              !
              if (iverbose>=6) write(out,'("igrid = ",i0)') igrid
              !
              ! the centrifugal factor will be needed for the L**2 term
              !
              !f_rot=b_rot/r(igrid)**2*sc
              !
              ! the diagonal term with the potential function
              !
              epot = omega_grid(iomega)%energy(ilevel,igrid)*sc
              !
              !call kinetic_energy_grid_points(ngrid,igrid,vibmat)
              !
              ielem = (ilevel-1)*ngrid+igrid
              !
              ! the diagonal matrix element will include PEC +L**2 as well as the vibrational kinetic contributions.
              vibmat(ielem,ielem) = vibmat(ielem,ielem) + epot
              !
            enddo
            !$omp end parallel do
            !
          enddo
          !
          if (iverbose>=4) call TimerStop('Build vibrational Hamiltonian')
          !
          istate = Omega_grid(iomega)%basis(1)%istate
          !
          if (job%vibmax(istate)>Ndimen/2) then
             !
             call lapack_syev(vibmat,vibener)
             !
             ! we need only these many roots
             Nroots = min(Ndimen,job%vibmax(istate))
             !
             ! or as many as below job%upper_ener if required by the input
             if ((job%vibenermax(istate))*sc<safe_max) then
               nroots = maxloc(vibener(:)-vibener(1),dim=1,mask=vibener(:).le.job%vibenermax(istate)*sc)
             endif
             !
           else
             !
             ! some diagonalizers needs the following parameters to be defined
             !
             ! diagonalize the vibrational hamiltonian using the DSYEVR routine from LAPACK
             ! DSYEVR computes selected eigenvalues and, optionally, eigenvectors of a real n by n symmetric matrix A.
             ! The matrix is first reduced to tridiagonal form, using orthogonal similarity transformations.
             ! Then whenever possible, DSYEVR computes the eigenspectrum using Multiple Relatively Robust Representations (MR).
             !
             jobz = 'V'
             vrange(1) = -0.0_rk ; vrange(2) = (job%vibenermax(istate))*sc
             if (.not.job%zShiftPECsToZero) vrange(1) = -safe_max
             irange(1) = 1 ; irange(2) = min(job%vibmax(istate),Ndimen)
             nroots = Ndimen
             rng = 'A'
             !
             if (job%vibmax(istate)/=1e8) then
                rng = 'I'
             elseif (job%vibenermax(istate)<1e8) then
                rng = 'V'
             endif
             !
             call lapack_syevr(vibmat,vibener,rng=rng,jobz=jobz,iroots=nroots,vrange=real(vrange,kind=8),irange=irange)
             !
          endif
          !
          if (nroots<1) then
            nroots = 1
            vibener = 0
            vibmat = 0
            vibmat(1,1) = 1.0_rk
          endif
          !
          ! ZPE is obatined only from the lowest state
          !
          if (ilevel==1) zpe = vibener(1)
          !
          ! write the pure vibrational energies and the corresponding eigenfunctions into global matrices
          contracted(iomega)%vector(:,1:nroots) = vibmat(:,1:nroots)
          contracted(iomega)%energy(1:nroots)   = vibener(1:nroots)
          contrenergy(totalroots+1:totalroots+nroots) = vibener(1:nroots)
          !
          !vibmat_rk = vibmat
          !
          !call schmidt_orthogonalization(ngrid,nroots,vibmat_rk)
          !
          !contrfunc_rk(:,totalroots+1:totalroots+nroots) = vibmat_rk(:,1:nroots)
          !
          ! assign the eigenstates with quanta
          do i=1,nroots
            !
            contracted(iomega)%ilevel(i) = totalroots+i
            !
            imaxcontr = maxloc(vibmat(:,i)**2,dim=1,mask=vibmat(:,i)**2.ge.small_)
            !
            ilevel = ceiling(real(imaxcontr,rk)/real(ngrid,rk))
            !
            icontrvib(totalroots + i)%ilevel =  ilevel
            icontrvib(totalroots + i)%omega  =  omega
            icontrvib(totalroots + i)%iomega =  iomega
            icontrvib(totalroots + i)%istate =  Omega_grid(iomega)%qn(ilevel)%istate
            icontrvib(totalroots + i)%ilambda=  Omega_grid(iomega)%qn(ilevel)%ilambda
            icontrvib(totalroots + i)%sigma  =  Omega_grid(iomega)%qn(ilevel)%sigma
            icontrvib(totalroots + i)%v = i-1
            !
          enddo
          !
          ! increment the global counter of the vibrational states
          !
          totalroots = totalroots + nroots
          !
          ! dealocate some objects
          !
          deallocate(vibmat,vibener)
          call ArrayStop('vibmat')
          call ArrayStop('vibener')
          !
          !deallocate(vibmat_rk)
          !call ArrayStop('vibmat_rk')
          !
          !allocate(matelem_rk(totalroots,totalroots),stat=alloc)
          !call ArrayStart('matelem_rk',alloc,size(matelem_rk),kind(matelem_rk))
          !
          !deallocate(vec)
          !call ArrayStop('vec')
          !          
       enddo
       !
       ! sorting basis states (energies, basis functions and quantum numbers) from different
       ! states all together according with their energies
       !
       do i = 1,totalroots
         !
         energy_ = contrenergy(i)
         !
         do j=i+1,totalroots
           !
           energy_j = contrenergy(j)
           !
           if ( energy_>energy_j ) then
             !
             ! energy
             !
             energy_=energy_j
             contrenergy(j) = contrenergy(i)
             contrenergy(i) = energy_
             !
             ! basis function
             !
             !vec(:) = contrfunc(:,j)
             !contrfunc(:,j) = contrfunc(:,i)
             !contrfunc(:,i) = vec(:)
             !
             ! qunatum numbers
             !
             icontrvib_ = icontrvib(j)
             icontrvib(j) = icontrvib(i)
             icontrvib(i) = icontrvib_
             !
           endif
           !
         enddo
         !
       enddo
       !
       ! print out the vibrational fields in the J=0 representaion
       if (iverbose>=4) then
          write(out,'(/"Vibrational (contracted) energies: ")')
          write(out,'("    i        Energy/cm    State v"/)')
          do i = 1,totalroots
            ilevel = icontrvib(i)%ilevel
            write(out,'(i5,f18.6," [ ",2i4,f8.1," ] ",a)') i,( contrenergy(i)-contrenergy(1))/sc,ilevel,icontrvib(i)%v,&
                                                   icontrvib(i)%omega,trim(poten(icontrvib(i)%istate)%name)
          enddo
       endif
       !
       if (job%print_vibrational_energies_to_file ) then
          !
          call IOstart("print vibrational energies to file",u1)
          !
          open(unit=u1, file='J0_vibrational_energies.dat',status='replace',action='write')
          write(u1,'(/"Vibrational (contracted) energies: ")')
          write(u1,'("    i        Energy/cm    State v"/)')
          do i = 1,totalroots
            ilevel = icontrvib(i)%istate
            write(u1,'(i5,f18.6," [ ",2i4,f8.1," ] ",a)') i,(contrenergy(i))/sc,ilevel,icontrvib(i)%v,&
                                                   icontrvib(i)%omega,trim(poten(icontrvib(i)%istate)%name)
          enddo
          close(u1)
          call IOstop("print vibrational energies to file")
       endif
       !
       !
       ! check the orthogonality of the basis
       !
       if (iverbose>=3) then
         !
         if (iverbose>=6) write(out,'(/"Check the contracted basis for ortho-normality")')
         !
         if (action%intensity.and.intensity%overlap) then 
           !
           write(out,'(/"Vibrational overlap integrals: ")')
           ! write(out,'("    State-i    <i|j>   State-j"/)')
           write(out,"(1x,a7,1x,a7,6x,a10)") 'State-i','State-j', '<i|j>'
           !
         endif
         !
         !omp parallel do private(ilevel,jlevel,psipsi_t) schedule(guided)
         do ilevel = 1,totalroots
           iomega = icontrvib(ilevel)%iomega
           i = icontrvib(ilevel)%v+1
           do jlevel = 1,ilevel

             jomega = icontrvib(jlevel)%iomega
             j = icontrvib(jlevel)%v+1
             !
             if (iomega/=jomega) cycle
             !
             psipsi_t  = sum(contracted(iomega)%vector(:,i)*contracted(jomega)%vector(:,j))
             !
             if (iverbose>=6) then
                if (icontrvib(ilevel)%ilevel/=icontrvib(jlevel)%ilevel&
                  .and.ilevel/=jlevel.and.abs(psipsi_t)>sqrt(small_)) then
                   write(out,"('orthogonality is brocken : <',i4,'|',i4,'> (',f16.6,')')") ilevel,jlevel,psipsi_t
                   stop 'Brocken orthogonality'
                endif
                !
                if (ilevel==jlevel.and.abs(psipsi_t-1.0_rk)>sqrt(small_)) then
                   write(out,"('normalization is brocken:  <',i4,'|',i4,'> (',f16.6,')')") ilevel,jlevel,psipsi_t
                   stop 'Brocken normalization'
                endif
             endif
             !
             ! Reporting the quality of the matrix elemenst
             !
             if (action%intensity.and.intensity%overlap.and.&
                 icontrvib(ilevel)%ilevel/=icontrvib(jlevel)%ilevel) then
                !
                write(out,'("<",i2,",",i4,"|",i2,",",i4,"> = ",es18.8)') icontrvib(ilevel)%ilevel,    &
                                                                         icontrvib(ilevel)%v,            &
                                                                         icontrvib(jlevel)%ilevel,       &
                                                                         icontrvib(jlevel)%v,            &
                                                                         psipsi_t 
             endif
             !
             if (iverbose>=6) then
               if (ilevel/=jlevel) then
                 write(out,"('<',i4,'|',i4,'> = ',e16.2,'<-',8x,'0.0')") ilevel,jlevel,psipsi_t
               else
                 write(out,"('<',i4,'|',i4,'> = ',f16.2,'<-',8x,'1.0')") ilevel,jlevel,psipsi_t
               endif
             endif
             !
           enddo
         enddo
         !omp end parallel do
         !
       endif
       !
       deallocate(Kinmat)
       call ArrayStop('kinmat')
       !
       if (iverbose>=4) call TimerStop('Solve vibrational part')
       !
  end subroutine Solve_vibrational_problem_for_Omega_states




  subroutine Compute_rovibronic_Hamiltonian_in_omega_vib_representation(iverbose,jval,ngrid,Ntotal,Nomega_states,&
                                                                        sc,icontr,contrenergy,hmat)
     !
     implicit none
     !
     integer(ik),intent(in)   :: iverbose,Nomega_states,ngrid,Ntotal
     real(rk),intent(in)     :: sc,jval
     type(quantaT),intent(in) :: icontr(Ntotal)
     real(rk),intent(in)      :: contrenergy(ngrid*Nomega_states)
     real(rk),intent(out)     :: hmat(Ntotal,Ntotal)
     !
     integer(ik) :: i,j,ivib,ilevel,jlevel,istate,jstate,ilambda,jlambda,imulti,jmulti,iomega,jomega,jvib,alloc
     integer(ik) :: iSplus_omega_,iLplus_omega_,iSR_omega_,isigmav,ilevel_,jlevel_,ibobrot,ip2q_omega_,iq_omega_
     integer(ik) :: ibrot_omega
     real(rk)  :: sigmai,sigmaj,omegai,omegaj,spini,spinj,f_rot,erot,omegai_,omegaj_,f_w,f_t,f_o2,f_o1,f_lo
     character(len=250),allocatable :: printout(:)
     character(cl)         :: printout_
     type(fieldT),pointer  :: field
      !
      allocate(printout(Nomega_states),stat=alloc) ; if (alloc/=0) stop 'cannot allocate printout'
      printout = ''

      hmat = 0
      !
      if (iverbose>=4) call TimerStart('Construct the hamiltonian')
      !
      if (iverbose>=3) write(out,'(/"Construct the hamiltonian matrix")')
      !
      !omp parallel do private(i,ivib,ilevel,istate,sigmai,imulti,ilambda,omegai,spini,jvib,jlevel,jstate,sigmaj,  & 
      !                        jmulti,jlambda,omegaj,spinj,f_rot,erot,iL2,field,f_l2,f_s,f_t,iso,ibraket,ipermute, &
      !                        istate_,ilambda_,sigmai_,spini_,jstate_,jlambda_,sigmaj_,spinj_,isigmav,omegai_,    &
      !                        omegaj_,itau,ilxly,f_grid,f_l,f_ss) shared(hmat) schedule(guided)
      do i = 1,Ntotal
        !
        ivib    = icontr(i)%ivib
        ilevel  = icontr(i)%ilevel
        !
        istate  = icontr(i)%istate
        sigmai  = icontr(i)%sigma
        imulti  = icontr(i)%imulti
        ilambda = icontr(i)%ilambda
        omegai  = icontr(i)%omega
        spini   = icontr(i)%spin
        ilevel  = icontr(i)%ilevel
        iomega  = icontr(i)%iomega
        omegai  = icontr(i)%omega
        !
        ! the diagonal contribution is the energy from the contracted vibrational solution
        !
        hmat(i,i) = contrenergy(ivib)
        !
        if (iverbose>=6) write(out,'("ilevel,ivib = ",2(i0,2x) )') ilevel,ivib
        !
        do j =i,Ntotal
           !
           jvib    = icontr(j)%ivib
           jlevel  = icontr(j)%ilevel
           jstate  = icontr(j)%istate
           sigmaj  = icontr(j)%sigma
           jmulti  = icontr(j)%imulti
           jlambda = icontr(j)%ilambda
           omegaj  = icontr(j)%omega
           spinj   = icontr(j)%spin
           jlevel = icontr(j)%ilevel
           jomega  = icontr(j)%iomega
           omegaj  = icontr(j)%omega
           !
           ! the centrifugal factor will be needed for different terms
           ! BRot centrifugal (rotational) term, i.e. a correction to f_rot
           !
           if (nint(omegai-omegaj)==0) then
             !
             do ibrot_omega = 1,NBRot_omega
                 !
                 field => BRot_omega_obj(ibrot_omega)
                 !
                 omegai_ = field%omegai
                 ilevel_ = field%ilevel
                 jlevel_ = field%jlevel
                 !
                 if ( nint(omegai_-omegai)/=0.or.ilevel_/=ilevel.or.jlevel_/=jlevel) cycle
                 !
                 f_rot = field%matelem(ivib,jvib)
                 !
                 erot = f_rot*( Jval*(Jval+1.0_rk) - omegai**2)
                 !
                 hmat(i,j) = hmat(i,j) + erot
                 !
             enddo
             !
             ! print out the internal matrix at the first grid point
             if (iverbose>=4.and.abs(hmat(i,j)) >small_) then
                 write(printout(ilevel),'(A, F15.3,A)') "RV=", hmat(i,j)/sc, "; "
             endif
             !
           endif
           !
           !f_rot=brot(1)%matelem(ivib,jvib)
           !
           ! diagonal elements
           !
           !if (nint(omegai-omegaj)==0.and.ilevel==jlevel) then
           !  !
           !  erot = f_rot*( Jval*(Jval+1.0_rk) - omegai**2)
           !  !
           !  ! add the diagonal matrix element to the local spin-rotational matrix hmat
           !  hmat(i,j) = hmat(i,j) + erot
           !  !
           !  ! print out the internal matrix at the first grid point
           !  if (iverbose>=4.and.abs(hmat(i,j)) >small_) then
           !      write(printout(ilevel),'(A, F15.3,A)') "RV=", hmat(i,j)/sc, "; "
           !  endif
           !  !
           !endif
           !
           ! BOB centrifugal (rotational) term, i.e. a correction to f_rot
           !
           if (nint(omegai-omegaj)==0) then
             !
             do ibobrot = 1,NBob_omega
                 !
                 field => Bob_omega_obj(ibobrot)
                 !
                 omegai_ = field%omegai
                 ilevel_ = field%ilevel
                 jlevel_ = field%jlevel
                 !
                 if ( nint(omegai_-omegai)/=0.or.ilevel_/=ilevel.or.jlevel_/=jlevel) cycle
                 !
                 f_rot = field%matelem(ivib,jvib)
                 !
                 erot = f_rot*( Jval*(Jval+1.0_rk) - omegai**2)
                 !
                 hmat(i,j) = hmat(i,j) + erot
                 !
             enddo
             !
           endif
           !
           ! J*S part (S-uncoupling)
           !
           if (abs(nint(omegaj-omegai))==1) then
             !
             do iSplus_omega_ = 1,NSplus_omega 
               !
               field => S_omega_obj(iSplus_omega_)
               !
               omegai_ = field%omegai
               omegaj_ = field%omegaj
               ilevel_ = field%ilevel
               jlevel_ = field%jlevel
               !
               if (ilevel_/=ilevel.or.jlevel_/=jlevel) cycle
               !
               do isigmav = 0,1
                 !
                 ! the permutation is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
                 ! avoid the double counting.
                 if( isigmav==1.and.nint( abs( field%omegai ) + abs( field%omegaj ) )==0 ) cycle
                 !           
                 ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                 omegai_ = omegai_*(-1)**isigmav
                 omegaj_ = omegaj_*(-1)**isigmav
                 !
                 ! proceed only if the quantum numbers of the field equal to the corresponding <i| and |j> quantum numbers:
                 !
                 if( nint(omegai_-omegai)/=0.or.nint(omegaj_-omegaj)/=0 ) cycle
                 !
                 f_w = nint(omegai-omegaj)
                 !
                 f_t = sqrt( jval* (jval +1.0_rk)-omegai*(omegai-f_w) )*field%matelem(ivib,jvib)
                 !
                 hmat(i,j) = hmat(i,j) - f_t
                 !
               enddo
               !
               ! print out the internal matrix at the first grid point
               if (iverbose>=4.and.abs(f_t)>sqrt(small_)) then
                  write(printout_,'("  J-S(",2i3,")=")') ilevel,jlevel
                  printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                  if (abs(hmat(i,j))>sqrt(small_)) then
                    write(printout_,'(F12.4, A)') -f_t/sc, " ;"
                    printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                  endif
               endif
               !
             enddo
             !
           endif
           !
           ! J*L parts
           !
           ! Also check that L+ is consistent with the selection rules
           !
           if (abs(nint(omegaj-omegai))==1) then
             !
             loop_iLomega : do iLplus_omega_ =1,NLplus_omega
               !
               field => L_omega_obj(iLplus_omega_)
               !
               omegai_ = field%omegai
               omegaj_ = field%omegaj
               ilevel_ = field%ilevel
               jlevel_ = field%jlevel
               !
               !
               f_t = field%matelem(ivib,jvib)
               !
               if (ilevel_/=ilevel.or.jlevel_/=jlevel) cycle
               !
               do isigmav = 0,1
                 !
                 ! the permutation is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
                 ! avoid the double counting.
                 if( isigmav==1.and. nint(abs( field%omegai ) + abs( field%omegaj ))==0 ) cycle
                 !           
                 ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                 omegai_ = omegai_*(-1)**isigmav
                 omegaj_ = omegaj_*(-1)**isigmav
                 !
                 ! proceed only if the quantum numbers of the field equal to the corresponding <i| and |j> quantum numbers:
                 !
                 if( nint(omegai_-omegai)/=0.or.nint(omegaj_-omegaj)/=0 ) cycle
                 !
                 ! L*J part (L-uncoupling)
                 !
                 !  omega should change by 1 (via J+/-) exactly as l_omega
                 f_w = nint(omegai-omegaj)
                 !
                 !f_t = field%matelem(ivib,jvib)
                 !f_t = sqrt( (jval-f_w*omegai)*(jval+f_w*omegai-1.0_rk) )*f_t
                 !
                 f_t = sqrt( jval* (jval +1.0_rk)-omegai*(omegai-f_w) )*field%matelem(ivib,jvib)
                 !
                 hmat(i,j) = hmat(i,j) - f_t
                 !
               enddo
               !
               ! print out the internal matrix at the first grid point
               if (iverbose>=4.and.abs(f_t)>small_) then
                  write(printout_,'(i3,"-LJ",2i3)') iLplus_omega_,ilevel,jlevel
                  printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                  write(printout_,'(g12.4)') -f_t/sc
                  printout(ilevel) = trim(printout(ilevel))//trim(printout_)
               endif
               !
             enddo loop_iLomega
             !
           endif
           !
           ! Spin-rotation 1st type (J*S)
           !
           ! selection rules
           !
           if (abs(nint(omegaj-omegai))==1) then
              !
              ! Non-diagonal spin-rotaion term
              !
              do iSR_omega_ = 1,NSR_omega
                !
                field => SR_omega_obj(iSR_omega_)
                !
                omegai_ = field%omegai
                omegaj_ = field%omegaj
                ilevel_ = field%ilevel
                jlevel_ = field%jlevel
                !
                if (ilevel_/=ilevel.or.jlevel_/=jlevel) cycle
                !
                ! Only one option is possible in the Omega representaion
                ! 1. <Omega|HSR|Omega+/-1>
                !
                !do isigmav = 0,1
                !  !
                !  ! the permutation is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
                !  ! avoid the double counting.
                !  if( isigmav==1.and. abs( field%iomega ) + abs( field%jomega )==0 ) cycle
                !!  !           
                !  ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                !  omegai_ = omegai_*(-1)**isigmav
                !  omegaj_ = omegaj_*(-1)**isigmav
                  !
                  ! proceed only if the quantum numbers of the field equal to the corresponding <i| and |j> quantum numbers:
                  !
                  if( nint(omegai_-omegai)/=0.or.nint(omegaj_-omegaj)/=0 ) cycle
                  !
                  f_w = nint(omegai-omegaj)
                  !
                  f_t = sqrt( jval* (jval +1.0_rk)-omegai*(omegai-f_w) )*field%matelem(ivib,jvib)
                  !
                  hmat(i,j) = hmat(i,j) + f_t*0.5_rk
                  !
                !enddo
                !
                ! print out the internal matrix at the first grid point
                if (iverbose>=4.and.abs(hmat(i,j))>sqrt(small_)) then
                   write(printout_,'("    SR",2i3)') ilevel,jlevel
                   printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                   if (abs(hmat(i,j))>sqrt(small_)) then
                     write(printout_,'(g12.4)') hmat(i,j)/sc
                     printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                   endif
                endif
                ! 
              enddo
              !
           endif
           !
           ! lambda doubling p2q 
           !
           if (abs(nint(omegaj-omegai))==1) then
             !
             loop_p2q_omega : do ip2q_omega_ =1,Np2q_omega
               !
               field => p2q_omega_obj(ip2q_omega_)
               !
               omegai_ = field%omegai
               omegaj_ = field%omegaj
               ilevel_ = field%ilevel
               jlevel_ = field%jlevel
               !
               f_t = field%matelem(ivib,jvib)
               !
               if (ilevel_/=ilevel.or.jlevel_/=jlevel) cycle
               !
               ! proceed only if the quantum numbers of the field equal to the corresponding <i| and |j> quantum numbers:
               !
               if( nint(omegai_-omegai)/=0.or.nint(omegaj_-omegaj)/=0 ) cycle
               !
               ! L*J part of p2q
               !
               !  omega should change by 1 (via J+/-) exactly as l_omega
               f_w = nint(omegai-omegaj)
               !
               f_t = sqrt( jval* (jval +1.0_rk)-omegaj*(omegaj+f_w) )*field%matelem(ivib,jvib)
               !
               hmat(i,j) = hmat(i,j) - f_t*0.5_rk
               !
               ! print out the internal matrix at the first grid point
               if (iverbose>=4.and.abs(f_t)>small_) then
                  write(printout_,'(i3,"-LJ",2i3)') ip2q_omega_,ilevel,jlevel
                  printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                  write(printout_,'(g12.4)') -f_t/sc
                  printout(ilevel) = trim(printout(ilevel))//trim(printout_)
               endif
               !
             enddo loop_p2q_omega
             !
           endif
           !
           ! lambda-doubling q
           !
           ! selection rules
           !
           if (abs(nint(omegaj-omegai))==2) then
              !
              do iq_omega_ = 1,Nq_omega
                !
                field => q_omega_obj(iq_omega_)
                !
                omegai_ = field%omegai
                omegaj_ = field%omegaj
                ilevel_ = field%ilevel
                jlevel_ = field%jlevel
                !
                if (ilevel_/=ilevel.or.jlevel_/=jlevel) cycle
                !
                ! Only one option is possible in the Omega representaion
                ! 1. <Omega|Hq|Omega+/-2>
                !
                ! proceed only if the quantum numbers of the field equal to the corresponding <i| and |j> quantum numbers:
                !
                if( nint(omegai_-omegai)/=0.or.nint(omegaj_-omegaj)/=0 ) cycle
                !
                f_o2 = omegaj-omegai
                f_o1 = sign(1.0_rk,f_o2)
                !
                f_t = sqrt( jval*(jval+1.0_rk)-(omegaj     )*(omegaj-f_o1) )*&
                      sqrt( jval*(jval+1.0_rk)-(omegaj-f_o1)*(omegaj-f_o2) )
                !
                f_lo = field%matelem(ivib,jvib)*f_t
                !
                hmat(i,j) = hmat(i,j) + f_lo*0.5_rk
                !
                ! print out the internal matrix at the first grid point
                if (iverbose>=4.and.abs(hmat(i,j))>sqrt(small_)) then
                   write(printout_,'("    q",2i3)') ilevel,jlevel
                   printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                   if (abs(hmat(i,j))>sqrt(small_)) then
                     write(printout_,'(g12.4)') hmat(i,j)/sc
                     printout(ilevel) = trim(printout(ilevel))//trim(printout_)
                   endif
                endif
                ! 
              enddo
              !
           endif
           !
        enddo  ! j
      enddo  ! i
      !omp end parallel do
      !
      if (iverbose>=4) call TimerStop('Construct the hamiltonian')
      !
      deallocate(printout)
      !

  end subroutine Compute_rovibronic_Hamiltonian_in_omega_vib_representation





  subroutine L_omega_create(NLplus_omega,onlycount)
    !
    implicit none
    !
    logical,intent(in) :: onlycount
    integer(ik),intent(inout) :: NLplus_omega
    !
    integer(ik)  :: ipermute,ilxly
    integer(ik)  :: i,j,istate_,jstate_,iLxy,imulti,jmulti,multi,multj
    real(rk)     :: spini_,spinj_,omegai,omegaj,omegai_,omegaj_
    real(rk)     :: sigmai,sigmaj
    integer(ik)  :: ilambda_,jlambda_,iomega,jomega,N_i,N_j
    type(fieldT),pointer       :: field
    !
    NLplus_omega = 0 
    !
    if (Nlxly==0) return
    !
    do iomega=1,Nomegas
      !
      omegai = Omega_grid(iomega)%omega
      !
      N_i = Omega_grid(iomega)%Nstates
      !
      do i = 1,N_i
        !
        do jomega=1,Nomegas
            !
            omegaj = Omega_grid(jomega)%omega
            !
            N_j = Omega_grid(jomega)%Nstates
            !
            do j = 1,N_j
              !
              ! assuming we will need L+ only 
              !
              if (nint(omegai-omegaj)/=1) cycle
              !
              ! check if any Lxy objects satisfy this criteria 
              !
              iLxy = 0 
              !
              loop_Lxy : do ilxly =1,Nlxly
                !
                field => lxly(ilxly)
                !
                ! Also check that L+ is consistent with the selection rules
                !
                if ( field%istate==field%jstate .or.abs(field%lambda-field%lambdaj)/=1 ) then
                   !
                   write(out,"('The quantum numbers of the L+/Lx field ',2i3,' are inconsistent" // &
                                   " with L+selection rules: ')") field%istate,field%jstate
                   write(out,"('Delta Lamda = +/-1')")
                   stop "Lx/L+ input is inconsistent with selection rules"
                   !
                endif
                !
                !
                do ipermute  = 0,1
                  !
                  if (ipermute==0) then
                    !
                    istate_ = field%istate ; ilambda_ = field%lambda  ; spini_ = field%spini
                    jstate_ = field%jstate ; jlambda_ = field%lambdaj ; spinj_ = field%spinj
                    !
                  else  ! permute
                    !
                    jstate_ = field%istate ; jlambda_ = field%lambda  ; spinj_ = field%spini
                    istate_ = field%jstate ; ilambda_ = field%lambdaj ; spini_ = field%spinj
                    !
                  endif
                  !
                  ! ipermute only makes sense for different states
                  !
                  if (ipermute==1.and.istate_==jstate_.and.ilambda_==jlambda_) cycle
                  !
                  multi = nint(2.0_rk*spini_+1.0_rk)
                  !
                  sigmai = -spini_-1.0_rk
                  do imulti = 1,multi
                    !
                    sigmai = sigmai + 1.0_rk
                    omegai_ = real(ilambda_,rk)+sigmai
                    !
                    multj = nint(2.0_rk*spinj_+1.0_rk)
                    sigmaj = -spinj_-1.0_rk
                    do jmulti = 1,multj
                      !
                      sigmaj = sigmaj + 1.0_rk
                      omegaj_ = real(jlambda_,rk)+sigmaj
                      !
                      if (nint(omegai_-omegai)==1.and.nint(omegaj_-omegaj)==1) then
                        iLxy = iLxy + 1
                        exit loop_Lxy 
                      endif
                      !
                    enddo
                    !
                  enddo
                  !
                enddo
                !
              enddo  loop_Lxy
              !
              if (iLxy==0) cycle
              !
              NLplus_omega = NLplus_omega + 1
              !
              iLplus_omega(iomega,jomega,i,j) = NLplus_omega
              !
              if (.not.onlycount) then
                 !
                 L_omega_obj(NLplus_omega)%ilevel = i
                 L_omega_obj(NLplus_omega)%omegai  = omegai
                 L_omega_obj(NLplus_omega)%iomega  = iomega
                 !
                 L_omega_obj(NLplus_omega)%jlevel = j
                 L_omega_obj(NLplus_omega)%omegaj  = omegaj
                 L_omega_obj(NLplus_omega)%jomega  = jomega
                 !
              endif
              !
           enddo
           !
        enddo  
        !
      enddo  
      !           
    enddo
    !
  end subroutine L_omega_create
  !
  !
  subroutine S_omega_create(NSplus_omega,onlycount)
    !
    implicit none
    !
    integer(ik),intent(inout) :: NSplus_omega
    !
    logical,intent(in) :: onlycount
    real(rk)     :: sigmai,sigmaj,spini
    real(rk)     :: omegai,omegaj,spinj
    integer(ik)  :: iomega,jomega,N_i,N_j,i,j
    integer(ik)  :: istate,jstate
    integer(ik)  :: ilambda,jlambda
    !
    NSplus_omega = 0 
    !
    do iomega=1,Nomegas
      !
      omegai = Omega_grid(iomega)%omega
      !
      N_i = Omega_grid(iomega)%Nstates
      !
      do i = 1,N_i
        !
        istate  = Omega_grid(iomega)%basis(i)%istate
        ilambda = Omega_grid(iomega)%basis(i)%ilambda
        sigmai  =  Omega_grid(iomega)%basis(i)%sigma
        spini   = Omega_grid(iomega)%basis(i)%spin
        !
        do jomega=1,Nomegas
            !
            omegaj = Omega_grid(jomega)%omega
            !
            N_j = Omega_grid(jomega)%Nstates
            !
            do j = 1,N_j
              !
              jstate  = Omega_grid(jomega)%basis(j)%istate
              jlambda = Omega_grid(jomega)%basis(j)%ilambda
              sigmaj  = Omega_grid(jomega)%basis(j)%sigma
              spinj   = Omega_grid(jomega)%basis(j)%spin
              !
              if (nint(spini-spinj)/=0) cycle
              !
              ! apply selection rules for S+ 
              !
              if (nint(sigmai-sigmaj)/=1.or.nint(spini-spinj)/=0.or.(ilambda/=jlambda)) cycle
              !
              NSplus_omega = NSplus_omega + 1
              !
              iSplus_omega(iomega,jomega,i,j) = NSplus_omega
              !
              if (.not.onlycount) then
                 !
                 S_omega_obj(NSplus_omega)%ilevel = i
                 S_omega_obj(NSplus_omega)%omegai  = omegai
                 S_omega_obj(NSplus_omega)%iomega  = iomega
                 !
                 S_omega_obj(NSplus_omega)%jlevel = j
                 S_omega_obj(NSplus_omega)%omegaj  = omegaj
                 S_omega_obj(NSplus_omega)%jomega  = jomega
                 !
              endif
              !
            enddo
            !
          enddo
          !
       enddo
       !
    enddo  
    !
  end subroutine S_omega_create
  !
  !
  subroutine SR_omega_create(NSR_omega,onlycount)
    !
    implicit none
    !
    integer(ik),intent(inout) :: NSR_omega
    !
    logical,intent(in) :: onlycount
    integer(ik)  :: ipermute,iSR
    integer(ik)  :: i,j,istate_,jstate_,iSR_
    real(rk)     :: spini_,spinj_,omegai,omegaj
    real(rk)     :: sigmai,sigmaj
    integer(ik)  :: ilambda_,jlambda_,iomega,jomega,N_i,N_j,isigmav
    type(fieldT),pointer       :: field
    !
    NSR_omega = 0 
    !
    if (Nsr==0) return
    !
    do iomega=1,Nomegas
      !
      omegai = Omega_grid(iomega)%omega
      !
      N_i = Omega_grid(iomega)%Nstates
      !
      do i = 1,N_i
        !
        do jomega=1,Nomegas
            !
            omegaj = Omega_grid(jomega)%omega
            !
            N_j = Omega_grid(jomega)%Nstates
            !
            do j = 1,N_j
              !
              if (abs(nint(omegai-omegaj))/=1) cycle
              !
              iSR_ = 0 
              !
              loop_SR : do isr = 1,Nsr
                !
                field => spinrot(isr)
                !
                ! Using the J-free part of the J*S term 
                ! 1. <Sigma,Omega,Lambda|HSR|Sigma+/-1,Omega+/-1,Lambda>
                ! For this term Lambda does not change
                !
                do ipermute  = 0,1
                   !
                   if (ipermute==0) then
                     !
                     istate_ = field%istate ; ilambda_ = field%lambda  ; spini_ = field%spini
                     jstate_ = field%jstate ; jlambda_ = field%lambdaj ; spinj_ = field%spinj
                     !
                   else  ! permute
                     !
                     jstate_ = field%istate ; jlambda_ = field%lambda  ; spinj_ = field%spini
                     istate_ = field%jstate ; ilambda_ = field%lambdaj ; spini_ = field%spini
                     !
                   endif
                   !
                   if (ilambda_/=jlambda_) cycle
                   !
                   ! however the permutation makes sense only when for non diagonal <State,Lambda,Spin|F|State',Lambda',Spin'>
                   ! otherwise it will cause a double counting:
                   !
                   if (ipermute==1.and.istate_==jstate_.and.ilambda_==jlambda_) cycle
                   !
                   do isigmav = 0,1
                     !
                     ! the permutation is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
                     ! avoid the double counting.
                     if( isigmav==1.and. abs( field%lambda ) + abs( field%lambdaj )==0 ) cycle
                     !                 
                     ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                     ilambda_ = ilambda_*(-1)**isigmav
                     jlambda_ = jlambda_*(-1)**isigmav
                     !
                     sigmai = omegai - real(ilambda_,rk)
                     sigmaj = omegaj - real(jlambda_,rk)
                     !
                     if (abs(sigmai)>spini_.or.abs(sigmaj)>spinj_) cycle
                     !
                     NSR_omega = NSR_omega + 1
                     !
                     iSR_omega(iomega,jomega,i,j) = NSR_omega
                     !
                     if (.not.onlycount) then
                        !
                        SR_omega_obj(NSR_omega)%ilevel = i
                        SR_omega_obj(NSR_omega)%omegai  = omegai
                        SR_omega_obj(NSR_omega)%iomega  = iomega
                        !
                        SR_omega_obj(NSR_omega)%jlevel = j
                        SR_omega_obj(NSR_omega)%omegaj  = omegaj
                        SR_omega_obj(NSR_omega)%jomega  = jomega
                        !
                     endif
                     !
                   enddo
                enddo
                ! 
              enddo loop_SR
              !
           enddo
           !
        enddo  
        !
      enddo  
      !           
    enddo
    !
  end subroutine SR_omega_create


  !
  subroutine Bob_omega_create(NBob_omega,onlycount)
    !
    implicit none
    !
    integer(ik),intent(inout) :: NBob_omega
    !
    logical,intent(in) :: onlycount
    integer(ik)  :: iBob
    integer(ik)  :: i,j,istate_,jstate_,iBob_
    real(rk)     :: spini_,spinj_,omegai
    real(rk)     :: sigmai
    integer(ik)  :: ilambda_,jlambda_,iomega,N_i,isigmav
    type(fieldT),pointer       :: field
    !
    NBob_omega = 0 
    !
    if (Nbobrot==0) return
    !
    do iomega=1,Nomegas
      !
      omegai = Omega_grid(iomega)%omega
      !
      N_i = Omega_grid(iomega)%Nstates
      !
      do i = 1,N_i
         !
         do j = 1,N_i
           !
           iBob_ = 0 
           !
           loop_Bob : do iBob = 1,Nbobrot
             !
             field => bobrot(iBob)
             !
             ! Bob is diagonal in everything
             !
             istate_ = field%istate ; ilambda_ = field%lambda  ; spini_ = field%spini
             jstate_ = field%jstate ; jlambda_ = field%lambdaj ; spinj_ = field%spinj
             !
             do isigmav = 0,1
               !
               ! the permutation is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
               ! avoid the double counting.
               if( isigmav==1.and.field%lambda==0 ) cycle
               !
               ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
               ilambda_ = ilambda_*(-1)**isigmav
               !
               sigmai = omegai - real(ilambda_,rk)
               !
               if (abs(sigmai)>spini_) cycle
               !
               NBob_omega = NBob_omega + 1
               !
               iBob_omega(iomega,i,j) = NBob_omega
               !
               if (.not.onlycount) then
                  !
                  Bob_omega_obj(NBob_omega)%ilevel = i
                  Bob_omega_obj(NBob_omega)%omegai  = omegai
                  Bob_omega_obj(NBob_omega)%iomega  = iomega
                  !
                  Bob_omega_obj(NBob_omega)%jlevel = j
                  Bob_omega_obj(NBob_omega)%omegaj  = omegai
                  Bob_omega_obj(NBob_omega)%jomega  = iomega
                  !
               endif
               !
             enddo
             ! 
           enddo loop_Bob
           !
        enddo  
        !
      enddo  
      !           
    enddo
    !
  end subroutine Bob_omega_create
  !
  !
  subroutine P2Q_omega_create(NP2Q_omega,onlycount)
    !
    implicit none
    !
    integer(ik),intent(inout) :: NP2Q_omega
    !
    logical,intent(in) :: onlycount
    integer(ik)  :: iP2Q
    integer(ik)  :: i,j,istate_,jstate_,iP2Q_
    real(rk)     :: spini_,spinj_,omegai,omegaj
    real(rk)     :: sigmai,sigmaj
    integer(ik)  :: ilambda_,jlambda_,iomega,jomega,N_i,N_j,isigmav
    type(fieldT),pointer       :: field
    !
    NP2Q_omega = 0 
    !
    if (Nlambdap2q==0) return
    !
    do iomega=1,Nomegas
      !
      omegai = Omega_grid(iomega)%omega
      !
      N_i = Omega_grid(iomega)%Nstates
      !
      do i = 1,N_i
        !
        do jomega=1,Nomegas
            !
            omegaj = Omega_grid(jomega)%omega
            !
            N_j = Omega_grid(jomega)%Nstates
            !
            if ( abs(nint(omegai-omegaj))/=1 ) cycle
            !
            do j = 1,N_j
              !
              iP2Q_ = 0 
              !
              loop_P2Q : do ip2q = 1,Nlambdap2q
                !
                field => lambdap2q(ip2q)
                !
                ! Using the J-free part of the J*S term of p2q
                ! <Sigma+/-1,Omega-/+1,Lambda=-/+1|HP2Q|SigmaOmega,Lambda=+/-1>
                !
                ! p2q is only for Lambda=+/-1 and only for the same electronic state
                !
                istate_ = field%istate ; ilambda_ = field%lambda  ; spini_ = field%spini
                jstate_ = field%jstate ; jlambda_ = field%lambdaj ; spinj_ = field%spinj
                !
                do isigmav = 0,1
                  !
                  ! the permutation is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
                  ! avoid the double counting. Although this should not happen for opq 
                  if( isigmav==1.and.field%lambda==0 ) cycle
                  !
                  ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                  ilambda_ = ilambda_*(-1)**isigmav
                  !
                  ! lambda should change sign in p2q
                  jlambda_ = -ilambda_
                  !
                  sigmai = omegai - real(ilambda_,rk)
                  sigmaj = omegaj - real(jlambda_,rk)
                  !
                  if (abs(sigmai)>spini_.or.abs(sigmaj)>spini_.or.abs(nint(sigmaj-sigmai))/=1) cycle 
                  if (nint(sigmaj-sigmai)/=nint(omegai-omegaj)) cycle
                  !
                  NP2Q_omega = NP2Q_omega + 1
                  !
                  iP2Q_omega(iomega,jomega,i,j) = NP2Q_omega
                  !
                  if (.not.onlycount) then
                     !
                     P2Q_omega_obj(NP2Q_omega)%ilevel  = i
                     P2Q_omega_obj(NP2Q_omega)%omegai  = omegai
                     P2Q_omega_obj(NP2Q_omega)%iomega  = iomega
                     !
                     P2Q_omega_obj(NP2Q_omega)%jlevel  = j
                     P2Q_omega_obj(NP2Q_omega)%omegaj  = omegaj
                     P2Q_omega_obj(NP2Q_omega)%jomega  = jomega
                     !
                  endif
                enddo
                ! 
              enddo loop_P2Q
              !
           enddo
           !
        enddo  
        !
      enddo  
      !           
    enddo
    !
  end subroutine P2Q_omega_create  
  !
  !
  subroutine Q_omega_create(NQ_omega,onlycount)
    !
    implicit none
    !
    integer(ik),intent(inout) :: NQ_omega
    !
    logical,intent(in) :: onlycount
    integer(ik)  :: iQ
    integer(ik)  :: i,j,istate_,jstate_,iQ_
    real(rk)     :: spini_,spinj_,omegai,omegaj
    real(rk)     :: sigmai,sigmaj
    integer(ik)  :: ilambda_,jlambda_,iomega,jomega,N_i,N_j,isigmav
    type(fieldT),pointer       :: field
    !
    NQ_omega = 0 
    !
    if (Nlambdaq==0) return
    !
    do iomega=1,Nomegas
      !
      omegai = Omega_grid(iomega)%omega
      !
      N_i = Omega_grid(iomega)%Nstates
      !
      do i = 1,N_i
        !
        do jomega=1,Nomegas
            !
            omegaj = Omega_grid(jomega)%omega
            !
            N_j = Omega_grid(jomega)%Nstates
            !
            if ( abs(nint(omegai-omegaj))/=2 ) cycle
            !
            do j = 1,N_j
              !
              iQ_ = 0 
              !
              loop_Q : do iq = 1,Nlambdaq
                !
                field => lambdaq(iq)
                ! Q-lambd-doubling
                ! <Sigma,Omega-/+2,Lambda=-/+1|Lambda-O|Sigma,Omega,-Lambda>
                !
                ! q is only for Lambda=+/-1 and only for the same electronic state
                !
                istate_ = field%istate ; ilambda_ = field%lambda  ; spini_ = field%spini
                jstate_ = field%jstate ; jlambda_ = field%lambdaj ; spinj_ = field%spinj
                !
                do isigmav = 0,1
                  !
                  ! the permutation is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
                  ! avoid the double counting. Although this should not happen for opq 
                  if( isigmav==1.and.field%lambda==0 ) cycle
                  !
                  ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                  ilambda_ = ilambda_*(-1)**isigmav
                  !
                  ! lambda should change sign in q
                  jlambda_ = -ilambda_
                  !
                  sigmai = omegai - real(ilambda_,rk)
                  sigmaj = omegaj - real(jlambda_,rk)
                  !
                  if (abs(sigmai)>spini_.or.abs(sigmaj)>spini_.or.abs(nint(sigmaj-sigmai))/=0) cycle 
                  if ((ilambda_-jlambda_)/=nint(omegai-omegaj)) cycle
                  !
                  NQ_omega = NQ_omega + 1
                  !
                  iQ_omega(iomega,jomega,i,j) = NQ_omega
                  !
                  if (.not.onlycount) then
                     !
                     Q_omega_obj(NQ_omega)%ilevel  = i
                     Q_omega_obj(NQ_omega)%omegai  = omegai
                     Q_omega_obj(NQ_omega)%iomega  = iomega
                     !
                     Q_omega_obj(NQ_omega)%jlevel  = j
                     Q_omega_obj(NQ_omega)%omegaj  = omegaj
                     Q_omega_obj(NQ_omega)%jomega  = jomega
                     !
                  endif
                enddo
                ! 
              enddo loop_Q
              !
           enddo
           !
        enddo  
        !
      enddo  
      !           
    enddo
    !
  end subroutine Q_omega_create    
  !
  !
  subroutine Kin_omega_create(NKin_omega,onlycount)
    !
    implicit none
    !
    integer(ik),intent(inout) :: NKin_omega
    !
    logical,intent(in) :: onlycount
    integer(ik)  :: iKin
    integer(ik)  :: i,j,istate_,jstate_,iKin_
    real(rk)     :: spini_,spinj_,omegai
    real(rk)     :: sigmai
    integer(ik)  :: ilambda_,jlambda_,iomega,N_i,isigmav
    type(fieldT),pointer       :: field
    !
    NKin_omega = 0 
    !
    if (nestates==0) return
    !
    do iomega=1,Nomegas
      !
      omegai = Omega_grid(iomega)%omega
      !
      N_i = Omega_grid(iomega)%Nstates
      !
      do i = 1,N_i
         !
         do j = 1,N_i
           !
           iKin_ = 0 
           !
           loop_Kin : do iKin = 1,nestates
             !
             field => poten(iKin)
             !
             ! Kin is diagonal in everything
             !
             istate_ = field%istate ; ilambda_ = field%lambda  ; spini_ = field%spini
             jstate_ = field%jstate ; jlambda_ = field%lambdaj ; spinj_ = field%spinj
             !
             do isigmav = 0,1
               !
               ! the permutation is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
               ! avoid the double counting.
               if( isigmav==1.and.field%lambda==0 ) cycle
               !
               ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
               ilambda_ = ilambda_*(-1)**isigmav
               !
               sigmai = omegai - real(ilambda_,rk)
               !
               if (abs(sigmai)>spini_) cycle
               !
               NKin_omega = NKin_omega + 1
               !
               iKin_omega(iomega,i,j) = NKin_omega
               !
               if (.not.onlycount) then
                  !
                  Kin_omega_obj(NKin_omega)%ilevel = i
                  Kin_omega_obj(NKin_omega)%omegai  = omegai
                  Kin_omega_obj(NKin_omega)%iomega  = iomega
                  !
                  Kin_omega_obj(NKin_omega)%jlevel = j
                  Kin_omega_obj(NKin_omega)%omegaj  = omegai
                  Kin_omega_obj(NKin_omega)%jomega  = iomega
                  !
               endif
               !
             enddo
             ! 
           enddo loop_Kin
           !
        enddo  
        !
      enddo  
      !           
    enddo
    !
  end subroutine Kin_omega_create


  subroutine Brot_omega_create(NBrot_omega,onlycount)
    !
    implicit none
    !
    integer(ik),intent(inout) :: NBrot_omega
    !
    logical,intent(in) :: onlycount
    integer(ik)  :: iKin
    integer(ik)  :: i,j,istate_,jstate_,iKin_
    real(rk)     :: spini_,spinj_,omegai
    real(rk)     :: sigmai
    integer(ik)  :: ilambda_,jlambda_,iomega,N_i,isigmav
    type(fieldT),pointer       :: field
    !
    NBrot_omega = 0 
    !
    if (nestates==0) return
    !
    do iomega=1,Nomegas
      !
      omegai = Omega_grid(iomega)%omega
      !
      N_i = Omega_grid(iomega)%Nstates
      !
      do i = 1,N_i
         !
         do j = 1,N_i
           !
           iKin_ = 0 
           !
           loop_Kin : do iKin = 1,nestates
             !
             field => poten(iKin)
             !
             ! Kin is diagonal in everything
             !
             istate_ = field%istate ; ilambda_ = field%lambda  ; spini_ = field%spini
             jstate_ = field%jstate ; jlambda_ = field%lambdaj ; spinj_ = field%spinj
             !
             do isigmav = 0,1
               !
               ! the permutation is only needed if at least some of the quanta is not zero. otherwise it should be skipped to
               ! avoid the double counting.
               if( isigmav==1.and.field%lambda==0 ) cycle
               !
               ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
               ilambda_ = ilambda_*(-1)**isigmav
               !
               sigmai = omegai - real(ilambda_,rk)
               !
               if (abs(sigmai)>spini_) cycle
               !
               NBrot_omega = NBrot_omega + 1
               !
               iBrot_omega(iomega,i,j) = NBrot_omega
               !
               if (.not.onlycount) then
                  !
                  Brot_omega_obj(NBrot_omega)%ilevel = i
                  Brot_omega_obj(NBrot_omega)%omegai  = omegai
                  Brot_omega_obj(NBrot_omega)%iomega  = iomega
                  !
                  Brot_omega_obj(NBrot_omega)%jlevel = j
                  Brot_omega_obj(NBrot_omega)%omegaj  = omegai
                  Brot_omega_obj(NBrot_omega)%jomega  = iomega
                  !
               endif
               !
             enddo
             ! 
           enddo loop_Kin
           !
        enddo  
        !
      enddo  
      !           
    enddo
    !
  end subroutine Brot_omega_create

  !
  subroutine schmidt_orthogonalization(dimen,nroots,mat)
      !
      integer(ik),intent(in) :: dimen,nroots
      real(rk),intent(inout) :: mat(dimen,dimen)
      integer(ik) :: ielem,jelem
      real(rk)   :: cross_prod,factor
      !
      do ielem =  1,nroots
        !
        cross_prod = sum(mat(:,ielem)**2)
        !
        factor = 1.0_rk/sqrt(cross_prod)
        !
        mat(:,ielem) = mat(:,ielem)*factor
        !
        !$omp parallel do private(jelem,cross_prod) shared(mat) schedule(dynamic)
        do jelem = ielem+1,nroots
          !
          cross_prod = sum(mat(:,ielem)*mat(:,jelem))
          !
          cross_prod = -cross_prod
          !
          mat(:,jelem) = mat(:,jelem)+mat(:,ielem)*cross_prod
          ! 
        enddo 
        !$omp end parallel do
        !
      enddo 
      !
  end subroutine schmidt_orthogonalization


!
! integration with Simpson rules 
!                                      
  function simpsonintegral_rk(npoints,f) result (si) 
    integer(ik),intent(in) :: npoints
    !
    real(rk),intent(in) :: f(0:npoints)
    !
    real(rk) :: si
    !
    integer(ik) :: i
    !
    real(rk) ::  feven,fodd,f0,fmax,h
      !
      h = 1.0_rk 
      !xmax/real(Npoints,kind=rk)  !   integration step, it is already inlcuded
      feven=0         
      fodd =0
      f0   =f(0)
      fmax =f(Npoints)

     !
     !  sum of odd and even contributions 
     !
     do i = 1,npoints-2,2
        fodd   = fodd  + f(i  )
        feven  = feven + f(i+1)
     enddo
     !
     fodd   = fodd  + f(npoints-1)
     !
     si =  h/3.0_rk*( 4.0_rk*fodd + 2.0_rk*feven + f0 + fmax)

  end function  simpsonintegral_rk


!
! integration with Simpson rules 
!                                      
  function simpsonintegral_ark(npoints,f) result (si) 
    integer(ik),intent(in) :: npoints
    !
    real(ark),intent(in) :: f(0:npoints)
    !
    real(rk) :: si
    !
    integer(ik) :: i
    !
    real(ark) ::  feven,fodd,f0,fmax,h
      !
      h = 1.0_ark 
      !xmax/real(Npoints,kind=rk)  !   integration step, it is already inlcuded
      feven=0         
      fodd =0
      f0   =f(0)
      fmax =f(Npoints)

     !
     !  sum of odd and even contributions 
     !
     do i = 1,npoints-2,2
        fodd   = fodd  + f(i  )
        feven  = feven + f(i+1)
     enddo
     !
     fodd   = fodd  + f(npoints-1)
     !
     si =  h/3.0_ark*( 4.0_rk*fodd + 2.0_ark*feven + f0 + fmax)

  end function  simpsonintegral_ark

 !
 
  function three_j(a,b,c,al,be,ga)

      real(rk) :: three_j
      real(rk),intent(in) :: a,b,c,al,be,ga
      !
      integer(ik):: newmin,newmax,new,iphase
      real(rk)   :: delta,clebsh,minus
      real(rk)   :: term,term1,term2,term3,summ,dnew,term4,term5,term6,delta_log,term16,termlog


      three_j=0
!
!     (j1+j2).ge.j and j.ge.abs(a-b)    -m=m1+m2    j1,j2,j.ge.0
!     abs(m1).le.j1    abs(m2).le.j2   abs(m).le.j
!
      if(c.gt.a+b) return
      if(c.lt.abs(a-b)) return
      if(a.lt.0.or.b.lt.0.or.c.lt.0) return
      if(a.lt.abs(al).or.b.lt.abs(be).or.c.lt.abs(ga)) return
      if(-1.0_rk*ga.ne.al+be) return
!
!
!     compute delta(abc)
!
!     delta=sqrt(fakt(a+b-c)*fakt(a+c-b)*fakt(b+c-a)/fakt(a+b+c+1.0_rk))
      delta_log = faclog(a+b-c)+faclog(a+c-b)+faclog(b+c-a)-faclog(a+b+c+1.0_rk)
      !
      delta=sqrt(exp(delta_log)) 
!
!
      !term1=fakt(a+al)*fakt(a-al)
      !term2=fakt(b-be)*fakt(b+be)
      !term3=fakt(c+ga)*fakt(c-ga)
      !
      !term=sqrt( (2.0_rk*c+1.0_rk)*term1*term2*term3 )
      !
      !
      term1=faclog(a+al)+faclog(a-al)
      term2=faclog(b-be)+faclog(b+be)
      term3=faclog(c+ga)+faclog(c-ga)
      !
      termlog = ( term1+term2+term3+delta_log )*0.5_rk
 
      term=sqrt( (2.0_rk*c+1.0_rk) )
!
!
!     now compute summation term
!
!     sum to get summation in eq(2.34) of brink and satchler.  sum until
!     a term inside factorial goes negative.  new is index for summation
!     .  now find what the range of new is.
!
!
      newmin=nint(max((a+be-c),(b-c-al),0.0_rk))
      newmax=nint(min((a-al),(b+be),(a+b-c)))
!
!
      summ=0
!
!
      do new=newmin,newmax
        !
        dnew=real(new,rk)
        !
        term4=faclog(a-al-dnew)+faclog(c-b+al+dnew)
        term5=faclog(b+be-dnew)+faclog(c-a-be+dnew)
        term6=faclog(dnew)+faclog(a+b-c-dnew)
        !
        term16=termlog-(term4+term5+term6)
        !
        summ=summ+(-1.0_rk)**new*exp(term16)
        !
      enddo
!
!     so clebsch-gordon <j1j2m1m2ljm> is clebsh
!
      clebsh=term*summ ! /sqrt(10.0_rk)
!
!     convert clebsch-gordon to three_j
!
      iphase=nint(a-b-ga)
      minus = -1.0_rk
      if (mod(iphase,2).eq.0) minus = 1.0_rk
      three_j=minus*clebsh/term

!     threej=(-1.0_rk)**(iphase)*clebsh/sqrt(2.0_rk*c+1.0_rk)
!
!
   end function three_j
  !
  !
  ! calculate factorial by log function 
  ! 
  function faclog(a)   result (v)
    real(rk),intent(in) ::  a
    real(rk)            :: v 
    integer(ik) j,k

    v=0
    k=nint(a)
    if(k>=2) then 
      do j=2,k
         v=v+log(real(j,rk))
      enddo 
    endif 
    
  end function faclog
  !
end module diatom_module



