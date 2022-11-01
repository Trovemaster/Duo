module accuracy
  implicit none
  private
  public sik, ik, hik, rk, ark, out, inp, safe_max,safe_min,max_exp, pi, twopi, cl, wl
  public accuracyInitialize, print_physical_constants
  public planck,avogno,vellgt,boltz,bohr,todebye, umatoau,uma, g_s
  public epsil,small_,sqrt2,sqrt3,rad,fititermax,aston,hartree,ev,my_fmt
  !
  integer, parameter :: sik         = selected_int_kind(4)       ! Small integers
  integer, parameter :: ik          = selected_int_kind(8)       ! "Normal" integers. This must map on
                                                                 ! C "int" type, or the DX interface won't
                                                                 ! work correctly.
  integer, parameter :: hik         = selected_int_kind(8)       ! "Pointer" integers - sufficient to store
                                                                 ! memory address
  integer, parameter :: drk         = selected_real_kind(12,25)  ! double    precision floating-point numbers 
  integer, parameter :: rk          = selected_real_kind(12,25)  ! single/double/quadruple  precision floating-point numbers
  integer, parameter :: ark         = selected_real_kind(25,32)  ! quadruple precision floating-point numbers
!  integer, parameter :: ark         = selected_real_kind(25,32)  ! double precision floating-point numbers
  integer, parameter :: inp         = 5                          ! Output I/O channel
  integer, parameter :: out         = 6                          ! Output I/O channel
  integer, parameter :: nfilelegendre = 101                      ! Damp-output channel for eigenfunction 


  real(rk), parameter  :: safe_max   = exp(log(huge(1.0_rk))*0.25_rk)   ! Largest number we want to work with
  real(rk), parameter  :: safe_min   = exp(log(tiny(1.0_rk))*0.25_rk)   ! Smallest number we want to work with
                                                                        ! (The somewhat convoluted syntax is used to standard F95 
                                                                        !  which prohibits non-integer exponents in this context)
  real(rk), parameter  :: max_exp    =log(safe_max)                     ! Largest number OK for exponentiating
  real(rk), parameter  :: small_     =epsilon(1.0_rk)                   ! a positive model number that is almost 
                                                                        ! negligible compared to unity in the current model
                                                                 ! epsil - antisymmetric tensor (Levi-Civita symbol)
  real(rk), parameter  :: epsil(3,3,3)=reshape( (/0,0,0,0,0,-1,0,1,0,0,0,1,0,0,0,& 
                                        -1,0,0,0,-1,0,1,0,0,0,0,0/), (/3,3,3/) )
  real(ark), parameter :: pi    = 4.0_ark * atan2(1.0_ark,1.0_ark)   ! PI=3.14...
  real(ark), parameter :: twopi = 2.0_ark * pi                     ! 2*PI=6.28...
  real(rk), parameter  :: sqrt2 = sqrt(2._rk)                     ! \sqrt{2}
  real(rk), parameter  :: sqrt3 = sqrt(3._rk)                     ! \sqrt{3}
  real(rk), parameter  :: rad   = 180._rk/pi                      ! radian = 180/pi
  integer, parameter   :: cl          = 80                        ! Max character string length
  integer, parameter   :: wl          = 500                       ! Very large max character string length 
  character(len=cl)    :: my_fmt                                  !text variable containing formats for reads/writes
  integer, parameter   :: fititermax  = 200                       ! Max number of iterations in different fittings 


  ! physical constants -- All constants updated 24 April 2015 from the NIST
  ! website http://physics.nist.gov/cuu/Constants/index.html
  ! They are the CODATA 2010 values
  real(drk), parameter :: planck     =  6.62606957e-27_rk         ! Planck constant in (non-SI) erg*second
  real(drk), parameter :: vellgt     =  2.99792458e+10_rk         ! Speed of light constant in (non-SI) cm/second
  real(drk), parameter :: avogno     =  6.02214129e+23_rk         ! Avogadro constant
  real(drk), parameter :: boltz      =  1.3806488e-16_rk          ! Boltzmann constant in (non-SI) erg/Kelvin
  real(drk), parameter :: bohr       =  0.52917721092_rk          ! bohr constant in Angstroms
  real(drk), parameter :: hartree    =  219474.6313708_rk         ! hartree in cm-1
  real(drk), parameter :: uma        =  1.660538921e-24_rk        ! unified atomic mass unit [=mass(C-12)/12 ] in grams
  real(drk), parameter :: e_charge   =  1.602176565e-19_rk        ! electron charge in Coulombs
  real(drk), parameter :: aston      =  planck/(8._rk*PI**2*vellgt*uma*1e-16_rk)  !rotational factor in cm-1 amu Ang^2
  real(drk), parameter :: umatoau    =  bohr**2 *hartree /(2._rk*aston)   ! uma (Daltons) to a.u. (electron masses)=~1822.888...
  real(drk), parameter :: me         =  uma/umatoau                       ! electron mass in grams
  real(drk), parameter :: alpha      =  2._rk*pi*bohr*hartree*1e-8_rk     ! fine structure alpha =~ 1/137.035999074
  real(drk), parameter :: todebye    =  e_charge*planck*1e17_rk/(2._rk*pi*me*alpha)   ! 1 a.u. of dipole in debye
  real(drk), parameter :: ev         =  1e7_rk*e_charge/(planck*vellgt)           ! ev in cm-1
  real(drk), parameter :: kJoule     =  1.e10_rk / (planck*vellgt*avogno)         ! 1kJ/mol in cm-1
  real(drk), parameter :: kCal       =  kJoule*4.184_rk   !1 Cal (thermochemical)= 4.184 Joules
  real(drk), parameter :: THz        =  1e12_rk / vellgt
  real(drk), parameter :: hc         =  vellgt*planck
  real(drk), parameter :: g_s        =  2.00231930436256_rk       ! spin Lande g-factor
  
  contains

    subroutine accuracyInitialize  ! This initialization routine is now superfluous but kept for compatibility
    !
    end subroutine accuracyInitialize
    !
    subroutine print_physical_constants ! (to be completed)
      write(out,'(A)') 'Values of physical constants used by DUO:'
      write(out,'(A40,ES20.12,2x,a12)') 'Planck constant h = ', planck, 'erg*second'
      write(out,'(A40,F20.4,2x,a12)')   'Speed of light c = ', vellgt, 'cm/second'
      write(out,'(A40,F20.13,2x,a12)')  'Bohr radius a0 = ', bohr, 'angstroms'
      write(out,'(A40,F20.10,2x,a12)')  'Hartree energy Eh = ', hartree, 'cm^-1'
      write(out,'(A40,ES20.12,2x,a12)') 'Unified atomic mass unit u = ', uma, 'grams'
      write(out,'(A40,F20.13,2x,a12)')  'Unified atomic mass unit u = ', umatoau, 'me'
      write(out,'(A40,ES20.12,2x,a12)') 'Electron charge e = ', e_charge, 'Coulombs'
      write(out,'(A40,ES20.12,2x,a12)') 'Boltzmann constant kB = ', boltz, 'erg/Kelvin'
      write(out,'(A40,F20.16,2x,a12)') 'Boltzmann constant kB = ', boltz/hc, 'cm^-1/Kelvin'
      write(out,'(A40,ES20.12,2x,a12)') "Avogadro constant = ", avogno, 'mol^-1'

      write(out,'(A40,F20.14)')   'Fine structure constant 1/alpha = ', 1._rk / alpha
      write(out,'(A40,ES20.12,2x,a12)') 'Electron mass me = ', me, 'grams'
      write(out,'(A40,F20.12,2x,a12)') 'a.u. of dipole e*a0 = ', todebye, 'debyes'
      write(out,'(A40,ES20.12,2x,a12)') "1 erg      => ",1e-7_rk , 'Joule'
      write(out,'(A40,ES20.12,2x,a12)') "1 erg      => ", 1._rk/hc, 'cm^-1'
      write(out,'(A40,F20.10,2x,a12)') "1 eV       => ", ev, 'cm^-1'
      write(out,'(A40,F20.10,2x,a12)') "1 kCal/mol => ", kcal, 'cm^-1'
      write(out,'(A40,F20.10,2x,a12)') "1 kJ/mol   => ", kjoule, 'cm^-1'
      write(out,'(A40,F20.10,2x,a12)') "1 THz      => ", thz, 'cm^-1'
      write(out,'(a)')

    end subroutine print_physical_constants
end module accuracy
