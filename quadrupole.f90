module quadrupole

use accuracy,      only : hik, ik, rk, ark, cl, out,&
                          vellgt, planck, avogno, boltz, pi, small_
use diatom_module, only : job, Intensity, quantaT, eigen, basis,&
                          nQuadrupoles, quadrupoletm, duo_j0, fieldT, poten,&
                          three_j, jmin_global,check_point_eigenfunc,check_point_dipoles,&
                          check_point_basis_set,compute_vib_integrals,fitting
use timer,         only : IOstart,IOstop, Arraystart, Arraystop, ArrayMinus,&
                          Timerstart, Timerstop, MemoryReport, &
                          TimerReport, memory_limit, memory_now
use symmetry,      only : sym, correlate_to_Cs
use dipole,        only : transitions_filter_from_fitting

private
public qm_tranint

type ElevelT
  integer(ik)  :: jind  ! J index
  integer(ik)  :: igamma
  integer(ik)  :: ilevel
end type ElevelT

type(ElevelT), allocatable  :: Elevel(:)

integer(ik)    :: nEigenLevels

! Wilfrid Somogyi: Routine for evaluating the quadrupole moment matrix elements
!                  and linestrength equations. Adapted from dipole.f90 for the
!                  case of a symmetric tensor of rank two.

contains

  subroutine qm_tranint
    ! computes the partition function and calls qm_intensity to calculate
    ! the electric quadrupole transition moments, linestrengths and
    ! intensities if required.

    implicit none

    ! counter, indexes, etc.
    real(rk)              :: jMin, jMax, jValMin, jVal_
    integer(ik)           :: nJ, jInd
    integer(ik)           :: indLevel, indGamma, indState, guParity, &
                             indSym

    ! calculation variables
    real(rk)              :: energy, beta, expEn, part
    integer(ik)           :: irrep,totalroots,Ndimen_vib,NVibBasis

    real(rk), allocatable :: jVal(:), qPart(:,:)

    ! assorted
    integer(ik)   :: info, iVerbose = 4

    ! find the min and max J values to compute transitions for
    jMax = maxval(Intensity%J(1:2))
    jMin = minval(Intensity%J(1:2))
    jMin = max(jmin_global, jMin)

    ! now find lowest possible J value of system to make sure all
    ! energies J > 0 are included in line list file regardless of the
    ! value Intensity%J, otherwise energy list numbering will not be
    ! consistent with the actual J values.
    jValMin = 0
    if  ( mod(nint(2.0_rk * job%j_list(1)), 2) /= 0 ) jValMin = 0.5_rk

    jValMin = max(jmin_global, jValMin)

    ! now count number of J states
    jVal_ =  jValMin
    nJ = 1
    do while ( jVal_ < jMax )
      jVal_ = jVal_ + 1.0_rk
      nJ = nJ +  1
    enddo

    ! prepare arrays to store J values and partition functions
    allocate(jVal(nJ), stat=info)
    if ( info /= 0 ) then
      stop 'qm_tranint allocation error: jVal - out of memory'
    endif

    allocate(qPart(20, nJ), stat=info)
    if ( info /= 0 ) then
      stop 'qm_tranint allocation error: qPart - out of memory'
    endif

    ! add J values to array
    jVal_ = jValMin
    jInd  = 1
    Jval(jind) = Jval_

    do while (Jval_<Jmax)
       jind = jind + 1
       Jval_ = Jval_ + 1.0_rk
       Jval(jind) = Jval_
    end do
    !
    if (trim(job%IO_eigen)=='READ') then
       !
       call check_point_eigenfunc('READ',iverbose,nJ,Jval,totalroots)
       !
       if (trim(job%IO_dipole)=='READ') then 
          call check_point_dipoles('READ',iverbose,Ndimen_vib)
       elseif(trim(job%IO_dipole)=='CALC') then
          call check_point_basis_set('READ',iverbose,NVibBasis)
          call compute_vib_integrals(NVibBasis)
       else
          stop 'IO_eigen = read with IO_dipole undefined is illegal'
       endif
       !
    else
       call duo_j0(iverbose,Jval)
    endif
    !
    !call duo_j0(iverbose, jVal)
    !
    if  ( job%shift_to_zpe ) then

      do jInd = 1, nJ

        jVal_ = jVal(jInd)

        do indGamma = 1, sym%NrepresCs
          !
          irrep = correlate_to_Cs(indGamma, guParity)
          !
          do indLevel = 1, eigen(jInd, indGamma)%Nlevels
            !            
            ! for homonuclear Nrepres = 4 and the irrep can be
            ! reconstructed from parity and g/u
            !
            indState = eigen(jInd, indGamma)%quanta(indLevel)%istate
            guParity = poten(indState)%parity%gu
            indSym   = correlate_to_Cs(indGamma, guParity)
            !
            if (Intensity%gns(indSym)<small_) cycle
            ! if energy of state < ZPE then adjust ZPE
            energy = eigen(jInd, indGamma)%val(indLevel)
            Intensity%ZPE = min(Intensity%ZPE, energy)

          enddo
        enddo
      enddo

      ! write out result
      if ( iverbose >= 2 ) then
        write(out, '(/"Zero point energy (ZPE) = ",f18.6,"&
          & (global zero, used for intensities)")') Intensity%ZPE
      endif

      if ( iverbose >= 4 ) then
        write(out, "(/'Partition function = ',f18.4,'&
          & T = ',f12.2)") Intensity%part_func, Intensity%temperature
      endif

    endif

    select case ( trim(Intensity%action) )

      case('ABSORPTION', 'EMISSION', 'TM')

        call Sort_levels(iVerbose, nJ, jVal(1:nJ))

        beta = planck * vellgt / (boltz * Intensity%temperature)

        if ( Intensity%part_func < small_ ) then

          Intensity%part_func = 0

          do jInd = 1, nJ

            jVal_ = jVal(jInd)

            do indGamma = 1, sym%NrepresCs

              do indLevel = 1, eigen(jInd, indGamma)%Nlevels

                energy =  eigen(jInd, indGamma)%val(indLevel)
                irrep = eigen(jInd, indGamma)%quanta(indLevel)%igamma

                ! for homonuclear Nrepres = 4 and the irrep can be
                ! reconstructed from parity and g/u
                indState = eigen(jInd, indGamma)%quanta(indLevel)%istate
                guParity = poten(indState)%parity%gu
                indSym   = correlate_to_Cs(indGamma, guParity)

                ! calculate the Boltzmann exponent of the energy
                expEn = exp( -(energy - Intensity%ZPE) * beta)

                ! add to partition function
                Intensity%part_func = Intensity%part_func &
                  & + Intensity%gns(indSym) &
                  & * (2.0_rk*jVal_ + 1.0_rk) * expEn

              enddo
            enddo
          enddo

          ! write out result
          if ( iVerbose >= 4 ) then
            write(out, &
              "(/'Partition function = ',f18.4,' T = ',f12.2)" &
              ) intensity%part_func, intensity%temperature
          endif

        endif

        ! if absorption, emission of tm then calculate linestrengths
        call  qm_intensity(jVal, iVerbose)
        write(out, '(/a)') 'done'

      case('PARTFUNC')

        write(out, '(/a)') 'compute partition function'

        qPart = 0

        do jInd = 1, nJ

          jVal_ = jVal(jInd)

          do indGamma = 1, sym%NrepresCs

            do indLevel = 1, eigen(jInd, indGamma)%Nlevels

              energy =  eigen(jInd, indGamma)%val(indLevel)
              irrep = eigen(jInd, indGamma)%quanta(indLevel)%igamma

              ! for homonuclear Nrepres = 4 and the irrep can be
              ! reconstructed from parity and g/u
              indState = eigen(jInd, indGamma)%quanta(indLevel)%istate
              guParity = poten(indState)%parity%gu
              indSym   = correlate_to_Cs(indGamma, guParity)

              ! calculate the Boltzmann exponent of the energy
              expEn = exp( -(energy - Intensity%ZPE) * beta)

              ! add to partition function
              qPart(irrep, jInd) = qPart(irrep, jInd) &
                & + Intensity%gns(indSym) &
                & * (2.0_rk*jVal_ + 1.0_rk) * expEn

            enddo
          enddo
        enddo

        ! sum of partition function
        part = sum(qPart(:,:))

        do jInd = 1, nJ
          do irrep = 1, sym%NrepresCs
            write(out, '(i4,1x,f18.1,1x,es16.8)') &
              irrep, jVal(jInd), qPart(irrep, jInd)
          enddo
        enddo

        ! write out result
        write(out, '(es16.8)') part

    end select

    call MemoryReport
    call TimerReport
  end subroutine qm_tranint

  subroutine qm_intensity(Jval, iVerbose)
    ! performs the actual transition moment, linestrength and intensity
    ! calculations when called by qm_tranint

    implicit none

    ! I/O variables
    real(rk), intent(in)      :: jVal(:)
    integer(ik), intent(in)   :: iVerbose

    ! filenames, identifiers, etc.
    character(len=cl)         :: filename, ioname
    integer(ik)               :: enunit, transUnit,intunit, info
    integer(ik), allocatable  :: richUnit(:, :)
    character(len=130)        :: myFmt
    character(len=12)         :: char_Jf, char_Ji, char_LF
    character(len=1)          :: letLFa, letLFb
    character(len=2)          :: letLF, dir
    integer(ik)               :: nDecimals, alloc_p

    ! indexes and counters
    integer(ik)               :: nJ, Jmax_, indJ, indI, indF, IDj
    integer(ik)               :: iLFa, iLFb, iLF, iflag_rich
    integer(ik)               :: nTrans, indTrans, nLower
    integer(ik)               :: indGamma, guParity, indTau, &
                                 indGammaI, indGammaF, indSymI, indSymF
    integer(ik)               :: nLevels, nLevelsI, nLevelsF,&
                                 indLevelI, indLevelF
    integer(ik)               :: nLevelsG(sym%Nrepresen)
    integer(ik)               :: indRoot, k, k_

    ! (pseudo) quantum numbers
    real(rk), allocatable     :: vecI(:), vecF(:)
    type(quantaT), pointer    :: quantaI, quantaF
    real(rk)                  :: j_, jI, jF
    integer(ik)               :: ivibI, ivibF, vI, vF, ilambdaI, ilambdaF,parityI
    integer(ik)               :: vF_, ilambdaF_,iEntry_fitting
    real(rk)                  :: spinI, spinF, sigmaI, sigmaF, &
                                 omegaI, omegaF
    real(rk)                  :: spinF_, sigmaF_, iomegaF_
    character(len=10)         :: statename
    character(len=1)          :: ef, pm, branch

    ! calculation variables
    real(rk)                  :: energyI, energyF
    integer(ik)               :: istateI, istateF
    integer(ik)               :: nRepresen, dimenMax, dimenI, dimenF, &
                                 iGammaPair(sym%nRepresen)
    real(hik)                 :: matSize
    real(rk)                  :: beta, inten_cm_mol, emcoef, A_coef_s_1, &
                                 A_einst, boltz_fc, absorption_int, unitConv, vacPerm
    real(rk)                  :: lande, nu_if, lineStr, linestr2, tm, ddot
    real(rk), allocatable     :: halfLineStr(:)
     
    ! cash transitons
    real(rk),allocatable :: acoef_RAM(:),nu_ram(:)  
    integer(ik),allocatable :: indexi_RAM(:),indexf_RAM(:)
    ! logicals
    logical                   :: intSpin = .true.
    logical                   :: passed, passed_
    !
    ! bound/unbound,integer/non-integer labels
    character(len=1) :: label_bound
    character(len=7) :: print_J,print_omega,print_sigma

    call TimerStart('Intensity calculations')

    ! define values of some constants
    beta         = planck * vellgt / (boltz * Intensity%temperature)
    inten_cm_mol = 8.0d-36*pi**3 / (3.0_rk * planck * vellgt)
    emcoef       = planck*vellgt/(4.0_rk*pi)

    !
    ! vacuum permittivity (NIST 2018) - needs to be
    ! properly programmed later
    vacPerm = 8.8541878128d-12

    ! conversion factor for Q[a.u] -> Q[S.I],
    ! h[erg.s] -> h[J.s] and nu_if[/cm] -> nu_if[/m]
    unitConv = 2.012914458d-62

    ! calculate the common factor for the Einstein coefficient
    A_coef_s_1 = unitConv*(8.0_rk * pi**5)/(5.0_rk * vacPerm * planck)
    !

    if ( sym%maxdegen > 2) then
      write(out, "('qm_intensity: this procedure has not been tested&
            & for the symmetries with degeneracies higher than 2...&
            & In fact, this procedure has not been tested at all!')" &
            )
    endif

    nRepresen = sym%NrepresCs

    ! define number of J states from input
    nJ = size(Jval)

    if  ( trim(Intensity%linelist_file) /= 'NONE' ) then

      ! prepare and open the .states file
      filename = trim(Intensity%linelist_file)//'.states'
      write(ioname, '(a, i4)') 'Energy file'
      call IOStart(trim(ioname), enunit)
      open(unit = enunit, action='write', &
           status='replace', file=filename)

      ! prepare and open the .trans file
      filename = trim(Intensity%linelist_file)//'.trans'
      write(ioname, '(a, i4)') 'Transition file'
      call IOStart(trim(ioname), transUnit)
      open(unit = transUnit, action='write', &
           status='replace', file=filename)

      ! calculate matrix elements
      if ( Intensity%matelem ) then

        ! define maximum J value
        Jmax_ = nint( maxval(Jval(:)) )

        allocate( richUnit(nJ, nJ) )

        ! loop over initial J states
        do indI = 1, nJ

          ! assign value of initial J and write as string to char_Ji
          jI = Jval(indI)
          write(char_Ji, '(i12)') nint(jI)

          do indF = 1, nJ

            ! assign value of final J
            ! only loop over jF > jI, permute to calculate others
            jF = Jval(indF)
            if ( jF < jI) cycle

            ! J selection rules
            if (     nint(abs(jI - jF)) > 2 &
                .or. nint(jI + jF) < 2 &
                ) cycle

            ! write to file final J value as string to char_Jf
            write(char_Jf, '(i12)') nint(jF)

            ! New RichMol format - one file for all components
            filename = 'matelem_Q'//&
                       '_j'//trim(adjustl(char_Ji))//&
                       '_j'//trim(adjustl(char_Jf))//&
                       '_'//trim(Intensity%linelist_file)//'.rchm'

            ! open RichMol file for matrix elements between jI, jF
            call IOstart(trim(filename), richUnit(indI, indF))
            open(unit = richUnit(indI, indF), action='write', &
                 status = 'replace', file=filename)

            ! add headers to RichMol file
            write(richUnit(indI, indF), "('Start richmol format')")
            write(richUnit(indI, indF), "('Q ','   2','    9')")
            write(richUnit(indI, indF), "('M-tensor')")

            ! nine cartesian LF-components
            do iLFa = 1, 3

              letLFa = 'x'
              if ( iLFa > 3) letLFa = 'y'
              if ( iLFa > 6) letLFa = 'z'

              do iLFb = 1, 3

                letLFb = 'x'
                if ( iLFb == 2) letLFb = 'y'
                if ( iLFb == 3) letLFb = 'z'

                ! concatenate first and second index and convert iLFa
                ! and iLFb to a single running index over all 9 comps.
                ! 1=xx, 2=xy, 3=xz, 4=yx, 5=yy, 6=yz, 7=zx, 8=zy, 9=zz
                letLF = letLFa//letLFb
                iLF = 3*(iLFa - 1) + iLFb

                ! flag for imaginary components
                iflag_rich = -1
                if (     iLFa == 2 .and. iLFb /= 2 &
                    .or. iLFa /= 2 .and. iLFb == 2 &
                    ) iflag_rich  = 0

                write(char_LF, '(i12)') iLF
                write(richUnit(indI, indF), "('alpha',i5,i3,1x,a2)") &
                      iLF, iflag_rich, letLF

                ! calculate the part of the matrix element equation
                ! that corresponds to the molecule -> lab transformation
                call do_LF_matrix_elements(iLFa, iLFb, &
                                           richUnit(indI, indF), &
                                           jI, jF)
              enddo

              ! heading for next section
              write(richUnit(indI, indF), "('K-tensor')")

            enddo
          enddo
        enddo
      endif
    endif
    !
    if (intensity%use_fitting) then
       !
       filename =  trim(fitting%output_file)//'.int'
       write(ioname, '(a)') 'Selected intensities'
       call IOstart(trim(ioname),intunit)
       open(unit = intunit, action = 'write',status='replace' , file = filename)
       !
    endif
    !
    ! estimate the maximum size of the basis set
    dimenMax = 0
    do indJ = 1, nJ
      do indGamma = 1, nRepresen
        dimenMax = max(dimenmax, eigen(indJ, indGamma)%Ndimen)
      enddo
    enddo

    ! count the number of transitions that need to be calculated first
    ! this will help keep track of the calculation progress, also count
    ! the number of lower levels from which transitions occur
    nTrans = 0
    nLower = 0

    indRoot = 0

    ! number of initial states
    nLevels = nEigenLevels

    ! for a given symmetry, iGamma, with some g_ns(iGamma) we find
    ! it's counterpart jGamma /= iGamma with the same g_ns(iGamma). It
    ! is assumed that there is only one such pair in the case of
    ! absorption and emission calculations
    call find_igamma_pair(iGammaPair)

    call TimerStart('Intens_Filter-1')

    do indI =  1, nJ

      jI = JVal(indI)
      !
      do indGammaI = 1, nRepresen

        nLevelsI = eigen(indI, indGammaI)%Nlevels

        do indLevelI = 1, nLevelsI

          ! obtain the energy and quanta of the initial state
          energyI = eigen(indI, indGammaI)%val(indLevelI)

          ! obtain the symmetry of the initial state
          istateI = eigen(indI, indGammaI)%quanta(indLevelI)%iState
          guParity = poten(istateI)%parity%gu
          indSymI = correlate_to_Cs(indGammaI, guParity)
          quantaI => eigen(indI,indGammaI)%quanta(indLevelI)

          ! check energy of lower state is in range and > ZPE
          call energy_filter_ul(jI, energyI, passed, 'lower')
          !
          ! skipping unbound lower states if the bound filter is on
          if ( intensity%bound_filter.and.intensity%bound .and. .not. quantaI%bound ) passed = .false.

          if ( .not. passed ) cycle

          nLower = nLower + 1

          do indF = 1, nJ

            jF = jVal(indF)

            do indGammaF = 1, nRepresen

              nLevelsF = eigen(indF, indGammaF)%Nlevels

              do indLevelF = 1, nLevelsF

                ! obtain the energy and quanta of the final state
                energyF = eigen(indF, indGammaF)%val(indLevelF)

                ! obtain the symmetry of the final state
                istateF = eigen(indF, indGammaF)&
                            %quanta(indLevelF)%iState
                guParity = poten(istateF)%parity%gu
                indSymF = correlate_to_Cs(indGammaF, guParity)
                quantaF => eigen(indF,indGammaF)%quanta(indLevelF)

                ! check the Intensity of the transition passes filter
                call intens_filter(jI, jF, energyI, energyF, &
                                   indSymI, indSymF, iGammaPair, &
                                   passed)
                !
                ! skipping unbound upper states if the bound filter is on
                !if ( intensity%bound_filter.and.intensity%bound .and. .not.quantaF%bound ) passed = .false.
                !
                ! skip if both the upper and lower states are bound states if the unbound filter is on
                !if ( intensity%bound_filter.and.intensity%unbound.and.(quantaF%bound.and.quantaI%bound) ) passed = .false.
                !
                if ( intensity%bound_filter) then
                   ! skipping unbound upper states if the bound filter is on                  
                   if( intensity%bound .and. .not.quantaF%bound ) passed = .false.
                   !
                   ! skip if both the upper and lower states are bound states if the unbound filter is on
                   if ( intensity%unbound.and.(quantaF%bound.and.quantaI%bound) ) passed = .false.
                   !
                   ! skip if the upper states is bound states for unbound_upper as true 
                   if ( intensity%unbound_upper.and.quantaF%bound ) passed = .false.
                   !
                endif
                !
                if (intensity%use_fitting) then
                  !
                  quantaF => eigen(indF,indGammaF)%quanta(indLevelF)
                  !
                  call transitions_filter_from_fitting(jI,jF,indSymI,indSymF,indLevelI,indLevelF,&
                       energyI,energyF,quantaI,quantaF,passed)

                  !
                endif
                !
                if ( Intensity%matelem ) then
                  call matelem_filter(jI, jF, energyI, energyF, &
                                      indSymI, indSymF, iGammaPair, &
                                      passed)
                endif

                if ( passed ) then

                  nTrans = nTrans + 1

                endif
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo

    call TimerStop('Intens_Filter-1')

    ! we now count the number of states with a given symmetry
    nLevelsG = 0

    if ( mod(eigen(1,1)%quanta(1)%imulti, 2) == 0) intSpin = .false.

    allocate(vecI(dimenMax), stat = info)
    call ArrayStart('intensity-vecI', info, size(vecI), kind(vecI))

    do indI = 1, nJ

      jI = jVal(indI)
      j_ = -1 ! for output in RichMol format (w/ mat. elems.)

      do indGammaI = 1, nRepresen

        ! number of initial levels
        nLevelsI = eigen(indI, indGammaI)%Nlevels
        dimenI   = eigen(indI, indGammaI)%Ndimen

        do indLevelI = 1, nLevelsI

          ! obtain the energy of the initial state
          energyI = eigen(indI, indGammaI)%val(indLevelI)
          !
          ! obtain the symmetry of the state
          quantaI => eigen(indI,indGammaI)%quanta(indLevelI)
          istateI = quantaI%istate
          guParity = poten(istateI)%parity%gu
          indSymI = correlate_to_Cs(indGammaI, guParity)

          ! ignore states with zero nuclear spin statistical weighting
          if  ( Intensity%gns(indSymI) < small_ ) cycle
          !
          ! skipping unbound lower states if the bound filter is on
          !if (intensity%bound_filter.and.intensity%bound.and..not.quantaI%bound) cycle
          !
          if ( intensity%bound_filter) then
             !
             ! skipping unbound lower states if the bound filter is on                
             if( intensity%bound .and. .not.quantaI%bound ) passed = .false.
             !
             ! skip if the lower state is unbound  but the filter unbound_lower is false
             if ( .not.intensity%unbound_lower.and..not.quantaI%bound ) passed = .false.
             !
          endif
          !
          ! indRoot is a running number over states
          indRoot = indRoot + 1
          eigen(indI, indGammaI)%quanta(indLevelI)%iroot = indRoot
          !
          call energy_filter_ul(jI, energyI, passed, 'upper')
          call energy_filter_ul(jI, energyI, passed_, 'lower')
          !
          if (      .not. passed &
              .and. .not. passed_) cycle

          if ( trim(Intensity%linelist_file) /= 'NONE') then

            ! assign quantum numbers for initial state
            quantaI   => eigen(indI, indGammaI)%quanta(indLevelI)
            ivibI      = quantaI%ivib
            vI        = quantaI%v
            spinI     = quantaI%spin
            sigmaI    = quantaI%sigma
            ilambdaI   = quantaI%ilambda
            omegaI    = quantaI%omega
            parityI   = quantaI%iparity
            statename = trim(quantaI%name)

            ! reconstruct +/- and e/f parities

            if ( parityI == 1) then
              pm = '-'
            else
              pm = '+'
            endif

            if ( mod( nint(2.0_rk*jI), 2 ) == 1 ) then
              indTau = mod( nint(jI - 0.5), 2 )
            else
              indTau = mod( nint(jI), 2 )
            endif

            if ( indTau == parityI ) then
              ef = 'e'
            else
              ef = 'f'
            endif

            ! the variable nDecimals determines the number of decimal
            ! places to which we print the energy levels. By default we
            ! use 6 decimals for energies up to 100,000 /cm, sacrificing
            ! more for higher energies in order that the value fits
            ! within the 12 allocated character spaces. The present
            ! format works for energies -10,000 /cm < E < 1e11 /cm.
            nDecimals = 6 - max(0, &
              int(log10(abs(energyI - Intensity%ZPE) + 1.d-6) - 4))

            ! 
            label_bound = "b"
            if (.not.quantaI%bound) label_bound = "u"
            !
            if (intSpin) then 
              write(print_J,"(i7)") nint(jI)
              write(print_sigma,"(i7)") nint(sigmaI)
              write(print_omega,"(i7)") nint(omegaI)
            else
              write(print_J,"(f7.1)") jI
              write(print_sigma,"(f7.1)") sigmaI
              write(print_omega,"(f7.1)") omegaI
            endif

            ! if requested, calculate and print the Lande g-factor for
            ! the selected eigenstate
            if ( Intensity%lande_calc ) then

              lande = 0

              vecI(1:dimenI) = eigen(indI, indGammaI)%&
                                vect(1:dimenI, indLevelI)

              if ( jI > 0 ) then

                do k = 1, dimenI

                  spinF   = basis(indI)%icontr(k)%spin
                  sigmaF  = basis(indI)%icontr(k)%sigma
                  ilambdaF = basis(indI)%icontr(k)%ilambda
                  omegaF  = basis(indI)%icontr(k)%omega
                  vF    = basis(indI)%icontr(k)%ivib

                  do k_ = 1, dimenI

                    spinF_   = basis(indI)%icontr(k)%spin
                    sigmaF_  = basis(indI)%icontr(k)%sigma
                    ilambdaF_ = basis(indI)%icontr(k)%ilambda
                    iomegaF_  = basis(indI)%icontr(k)%omega
                    vF_    = basis(indI)%icontr(k)%ivib
                    !!!!!
                    if (     ilambdaF /= ilambdaF_ &
                        .or. nint(spinF - spinF_) /= 0 &
                        .or. vF /= vF_ ) cycle

                    if ( k == k_ ) then

                      lande = lande + vecI(k)*vecI(k)*omegaF &
                              *(real(ilambdaF, rk) + 2.0023_rk*sigmaF)

                    elseif ( nint(abs(sigmaF_ - sigmaF)) == 1 ) then

                      lande = lande + vecI(k)*vecI(k_) &
                              *sqrt( &
                                spinF*(spinF + 1.0_rk) &
                                - sigmaF*(sigmaF + sigmaF_ - sigmaF) &
                              ) &
                              *sqrt( &
                                jI*(jI + 1.0_rk) &
                                - omegaF*(omegaF + iomegaF_ - omegaF) &
                              ) &
                              * (2.002319_rk/2.0_rk)
                    endif
                  enddo
                enddo

                lande = lande / ( jI*(jI + 1.0_rk) )

              endif

              write(myFmt,'(A,i0,a,a)') "(i12,1x,f12.",ndecimals,",1x,i6,1x,a7,1x,f13.6,1x,a1,1x,a1,",&
                            "1x,a10,1x,i3,1x,i2,1x,a7,1x,a7)"

              write(enunit,myFmt,advance="no") & 
                         indRoot,energyI-intensity%ZPE,nint(intensity%gns(indSymI)*( 2.0_rk*jI + 1.0_rk )),print_J,&
                         lande,pm,ef,statename,vI,ilambdaI,print_sigma,print_omega
                         !
              if (intensity%unbound) then                          
                 write(enunit,"(1x,a1)",advance="no") label_bound
                 if (intensity%use_bound_rmax) then
                   write(enunit,"(1x,f9.5)",advance="no") quantaI%r_exp
                 endif
                if (intensity%bound_eps_print) then
                   write(enunit,"(1x,g8.2)",advance="no") quantaI%epsilon
                endif
              endif
              !
              !
              write(enunit,"(a1)",advance="yes") ""
              !
            ! alternative format for RichMol matrix elements
            elseif ( Intensity%matelem ) then

              nDecimals = 6 - max(0, &
                int(log10(abs(energyI - Intensity%ZPE) + 1.d-6) - 4))

              if ( nint(2*jI) /= nint(2*j_) ) then
                j_ = jI
                IDj = 0
              endif

              IDj = IDj + 1
              quantaI%iJ_ID = IDj

              ! if integer spin, then integerise quantum numbers
              if ( intSpin ) then
                write(myFmt, '(a)') &
                  "(i6,1x,i8,1x,i2,1x,i2,3x,e21.14,5x,a4,i3,1x,a2,i4,&
                  &1x,a2,f8.4,1x,i6,1x,i6,1x,i4,1x,i6,1x,a1,1x,a10)"
                write(enunit, myFmt) &
                  nint(j_), IDj, parityI+1, 1, energyI-Intensity%ZPE, &
                  'tau:', parityI, 'j:', nint(j_), 'c', 1.000_rk, &
                  nint(omegaI), vI, ilambdaI, nint(sigmaI), pm, statename

              ! if not then write quantum numbers as reals
              else
                write(myFmt, '(A,i0,a)') &
                  "(i7,1x,i12,1x,i1,1x,i2,1x,f12.",ndecimals,&
                  ",1x,f7.1,1x,i6,1x,i4,1x,f7.1,1x,a1,1x,a10)"
                write(enunit, myFmt) &
                  nint(j_), IDj, parityI+1, 1, energyI-Intensity%ZPE, &
                  omegaI, vI, ilambdaI, sigmaI, pm, statename
              endif

            ! standard output format if matelem or lande not required
            else

              nDecimals = 6 - max(0, &
                int(log10(abs(energyI - Intensity%ZPE) + 1.d-6) - 4))

              ! if integer spin, then integerise quantum numbers


              write(myFmt,'(A,i0,a)') "(i12,1x,f12.",ndecimals,",1x,i6,1x,a7,1x,a1,1x,a1,1x,a10,1x,i3,1x,i2,1x,a7,1x,a7)"
              !
              write(enunit,myFmt,advance="no") & 
                        indRoot,energyI-intensity%ZPE,nint(intensity%gns(indSymI)*( 2.0_rk*jI + 1.0_rk )),print_J,&
                        pm,ef,statename,vI,(ilambdaI),print_sigma,print_omega
                        !
              if (intensity%unbound) then                          
                 write(enunit,"(1x,a1)",advance="no") label_bound
                 if (intensity%use_bound_rmax) then
                   write(enunit,"(1x,f9.5)",advance="no") quantaI%r_exp
                 endif
                if (intensity%bound_eps_print) then
                   write(enunit,"(1x,g8.2)",advance="no") quantaI%epsilon
                endif
              endif
              !
              write(enunit,"(a1)",advance="yes") ""
              !
            endif
            !
          endif

          istateI = eigen(indI, indGammaI)%quanta(indLevelI)%istate
          guParity = poten(istateI)%parity%gu
          indSymI = correlate_to_Cs(indGammaI, guParity)

          nLevelsG(indSymI) = nLevelsG(indSymI) + 1

        enddo
      enddo
    enddo

    deallocate(vecI)
    call ArrayStop('intensity-vecI')

    if ( trim(Intensity%linelist_file) /= "NONE") then
      close(enunit, status='keep')
    endif

    write(myFmt, '(a,i0,a)') &
      "('number of states for each sym = ',", sym%Nrepresen, "i8)"
    write(out, myFmt) nLevelsG(:)

    matSize = int( sum(nLevelsG(:)), hik )

    if ( iVerbose >= 4 ) then
      write(out, "(/'Quadrupole moment integration (i)...')")
    endif
    !
    if (intensity%states_only) then 
         write(out,"('The transition intensities are not requested (states_only option)')") 
         return 
    endif 
    !
    if ( nTrans == 0 ) then
      write(out, "('qm_intensity: the transition filters are too tight:&
            & no entry')")
      stop 'qm_intensity: the filters are too tight'
    endif

    write(out, "(/'...done!')")

    allocate(vecI(dimenMax), vecF(dimenMax), stat=info)
    call ArrayStart('intensity-vectors', info, size(vecI), kind(vecI))
    call ArrayStart('intensity-vectors', info, size(vecF), kind(vecF))

    !!! why is this duplicated from above?
    write(myFmt, '(a,i0,a)') &
      "('number of states for each sym = ',", sym%Nrepresen, "i8)"
    write(out, myFmt) nLevelsG(:)

    if ( iVerbose >= 0 ) then
      write(out, "('Total number of lower states = ',i8)") nLower
      write(out, "('Total number of transitions  = ',i8)") nTrans
    endif

    if ( iVerbose >= 0 ) then
      write(out, "(/'Statistical weight g_ns =',4f12.1)") &
        Intensity%gns(1:)
    endif

    ! To speed up line strength evaluation, we perform the calculation:
    !   S_{if} = | <i|a|f> |^2  = | \sum_{nm} C_in C_fm <n|a|m> |^2
    ! in three steps:
    !   1. Evaluate the expansion of the initial state:
    !         s_{im} = sum_{n} C_in <n|a|m>
    !   2. Evaluate the product with the expansion of the final state:
    !         s_{if} = sum_{m} C_fm s_{im}
    !   3. Square the result to obtain S_{if}
    !         S_{if} = s_{if}^2
    !
    ! The transitory object s_{im} we refer to as 'half linestrength',
    ! with corresponding variable "half_linestr".

    ! initialise array to store eigenvectors
    allocate(halfLineStr(dimenMax), stat=info)
    call ArrayStart('halfLineStr', info, &
      size(halfLineStr), kind(halfLineStr))

    if ( iVerbose >= 5) call MemoryReport

    !
    if (trim(intensity%linelist_file)/="NONE") then 
       write(out,"(/'This is a line list production only, intensity print-out is swtitched off')")
       write(out,"('To see intensities in the standard output remove the keyword LINELIST from INTENSITY'/)")
    else 
    
      ! prepare the table header
      write(out, "(/a,a,a,a)") &
        'Linestrength S(f<-i) [Debye**2],', &
        'Transition moments [Debye],', &
        'Einstein coefficient A(if) [1/s],', &
        'and Intensities [cm/mol]'
   
      ! depending on the case we have different file formats
      select case ( trim(intensity%action) )
   
        ! absorption lines
        case('ABSORPTION')
          write(out, &
            ! Fixed width output
            "(/t5,'J',t7,'Gamma <-',t18,'J',t21,'Gamma',t27,'Typ',t37,&
            &'Ei',t44,'<-',t52,'Ef',t64,'nu_if',8x,'S(f<-i)',10x,'A(if)',&
            &12x,'I(f<-i)',7x,'State v lambda sigma  omega <- State v &
            &lambda sigma  omega ')" &
            !
            ! CSV output
            ! "('dir, J_i, Gamma_i, J_f, Gamma_f, Branch,&
            ! & E_i, Ef, nu_if,&
            ! & S_fi, A_if, I_fi,&
            ! & state_f, v_f, lambda_f, sigma_f, omega_f,&
            ! & state_i, v_i, lambda_i, sigma_i, omega_i')"&
          )
          dir = '<-'
   
        ! emission lines
        case('EMISSION')
          write(out, &
            "(/t5,'J',t7,'Gamma ->',t18,'J',t21,'Gamma',t27,'Typ',t37,&
            &'Ei',t44,'->',t52,'Ef',t64,'nu_if',8x,'S(i->f)',10x,'A(if)',&
            &12x,'I(i->f)',7x,'State v lambda sigma  omega -> State v &
            &lambda sigma  omega ')" &
          )
          dir = '->'
   
        ! Transition moments
        case('TM')
          write(out, &
            "(/t4,'J',t6,'Gamma <-',t17,'J',t19,'Gamma',t25,'Typ',t35,&
            &'Ei',t42,'<-',t52,'Ef',t65,'nu_if',10x,'TM(f->i)')" &
          )
   
      end select

    endif

    deallocate(vecF)
    
    ! ------------------------------------------------
    ! now begin the actual line intensity calculations
    ! ------------------------------------------------

    ! counter for the no. transitions
    indTrans = 0

    ! loop over initial J states and assign corresponding J value
    do indI = 1, nJ
      jI = jVal(indI)
      !
      if (trim(intensity%linelist_file)/="NONE".and.iverbose>=4) write(out,"('J = ',f9.1)") jI
      !
      ! loop over symmetries
      do indGammaI = 1, nRepresen

        ! no. of levels in, and dimension of, basis for initial state
        nLevelsI = eigen(indI, indGammaI)%nLevels
        dimenI   = eigen(indI, indGammaI)%nDimen

        ! loop over final J states and assign corresponding J value
        do indF = 1, nJ
          jF = jVal(indF)

          ! J selection rules for quadrupole operator
          if (     abs(nint(jI - jF)) > 2 &
              .or. abs(nint(jI + jF)) < 2) cycle

          ! loop over symmetries
          do indGammaF = 1, nRepresen

            ! no. of levels in, and dimension of, basis for final state
            nLevelsF = eigen(indF, indGammaF)%nLevels
            dimenF   = eigen(indF, indGammaF)%nDimen

            ! loop over levels in the initial state
            loopLevelsI : do indLevelI = 1, nLevelsI

              ! energy and quantum numbers of initial state
              energyI =  eigen(indI, indGammaI)%val(indLevelI)
              quantaI => eigen(indI, indGammaI)%quanta(indLevelI)

              istateI  = quantaI%istate  ! electronic state
              ivibI    = quantaI%ivib    ! vibrational (contracted)
              vI      = quantaI%v       ! vibrational
              spinI   = quantaI%spin    ! electron spin
              sigmaI  = quantaI%sigma   ! spin projection
              ilambdaI = quantaI%ilambda ! e- orb. ang. mom. projection
              omegaI  = quantaI%omega   ! tot. ang. mom. proj. mol. ax.

              ! reconstruct symmetry for C2v case
              guParity = poten(istateI)%parity%gu
              indSymI = correlate_to_Cs(indGammaI, guParity)

              ! apply energy filter to initial (lower) state
              call energy_filter_ul(jI, energyI, passed, 'lower')
              !
              ! skipping unbound lower states if the bound filter is on
              !if (intensity%bound_filter.and.intensity%bound.and..not.quantaI%bound) passed = .false.
              !
              if ( intensity%bound_filter) then
                 !
                 ! skipping unbound lower states if the bound filter is on                
                 if( intensity%bound .and. .not.quantaI%bound ) passed = .false.
                 !
                 ! skip if the lower state is unbound  but the filter unbound_lower is false
                 if ( .not.intensity%unbound_lower.and..not.quantaI%bound ) passed = .false.
                 !
              endif
              !
              if ( .not. passed ) cycle loopLevelsI

              ! vector of basis state coefficients for inital state
              vecI(1:dimenI) = eigen(indI, indGammaI)&
                                %vect(1:dimenI, indLevelI)
              !
              halfLineStr = 0
              !
              ! -----
              ! Before the actual calculations we check if there are any
              ! allowed transitions from current level of jI to any
              ! levels of jF, if not skip initial J state level.

              ! loop over levels in the final state to check for trans.
              do indLevelF = 1, nLevelsF

                energyF =  eigen(indF, indGammaF)%val(indLevelF)
                quantaF => eigen(indF, indGammaF)%quanta(indLevelF)

                istateF  = quantaF%istate  ! electronic state
                ivibF    = quantaF%ivib    ! vibrational (contracted)
                vF      = quantaF%v       ! vibrational
                spinF   = quantaF%spin    ! electron spin
                sigmaF  = quantaF%sigma   ! spin projection
                ilambdaF = quantaF%ilambda ! e- orb. ang. mom. projection
                omegaF  = quantaF%omega   ! tot. ang. mom. proj. mol. ax

                ! reconstruct symmetry for C2v case
                guParity = poten(istateF)%parity%gu
                indSymF  = correlate_to_Cs(indGammaF, guParity)

                ! apply transition intensity filter, result of which is
                ! overidden by mat. elem. filter if we want mat. elems.
                call intens_filter(jI, jF, energyI, energyF, &
                  indSymI, indSymF, iGammaPair, passed)
                !
                ! skipping unbound upper states if the bound filter is on
                !if (intensity%bound_filter.and.intensity%bound.and..not.quantaF%bound) passed = .false.
                !
                if ( intensity%bound_filter) then
                   !
                   ! skipping unbound upper states if the bound filter is on                  
                   if( intensity%bound .and. .not.quantaF%bound ) passed = .false.
                   !
                   ! skip if the upper states is bound states for unbound_upper as true 
                   if ( intensity%unbound_upper.and.quantaF%bound ) passed = .false.
                   !
                endif
                !
                if (.not.passed) cycle 
                !
                if (intensity%use_fitting) then
                  !
                  quantaF => eigen(indF,indGammaF)%quanta(indLevelF)
                  !
                  call transitions_filter_from_fitting(jI,jF,indSymI,indSymF,indLevelI,indLevelF,&
                       energyI,energyF,quantaI,quantaF,passed)

                  !
                endif

                if ( Intensity%matelem ) then
                  call matelem_filter(jI, jF, energyI, energyF, &
                    indSymI, indSymF, iGammaPair, passed)
                endif

                if ( passed ) exit

              enddo

              ! if no transition then move to next initial J state level
              if ( .not. passed ) cycle

              ! end of transition checking, move to calculation
              ! -----

              select case ( trim(intensity%action) )

                case('ABSORPTION', 'EMISSION')

                  if (indSymF /= iGammaPair(indSymI)) cycle

                  if (      Intensity%J(1) + Intensity%J(2) > 0 &
                      .and. abs(nint(jI - jF)) <= 2 &
                      .and. nint(jI + jF) >= 2 &
                      ) then

                    call do_1st_half_linestrength(jI, jF, indI, indF, &
                      dimenI, dimenF, vecI(1:dimenI), halfLineStr)
                  endif

                case('TM')

                  stop 'TM is not yet coded'

              end select

              !
              ! allocating all different arrays to keep the data in RAM in the exomol and matelem format
              if (trim(intensity%linelist_file)/="NONE") then
                !
                allocate(acoef_RAM(nlevelsF),stat=info)
                call ArrayStart('swap:acoef_RAM',info,size(acoef_RAM),kind(acoef_RAM))
                allocate(nu_RAM(nlevelsF),stat=info)
                call ArrayStart('swap:nu_RAM',info,size(nu_RAM),kind(nu_RAM))
                allocate(indexf_RAM(nlevelsF),indexi_RAM(nlevelsF),stat=info)
                call ArrayStart('swap:indexf_RAM',info,size(indexf_RAM),kind(indexf_RAM))
                call ArrayStart('swap:indexi_RAM',info,size(indexi_RAM),kind(indexi_RAM))
                !
                acoef_RAM = 0
                !
              endif
              
              !$omp parallel private(vecF, alloc_p) shared(nu_ram,acoef_RAM,indexi_RAM,indexf_RAM)
              allocate(vecF(dimenMax), stat=alloc_p)

              if (alloc_p /= 0) then
                write(out, &
                  "('quadrupole: ',i9,' &
                  &trying to allocate array -vecF')") alloc_p
                stop 'quadrupole-vecF - out of memory'
              endif

              !$omp  do private(indLevelF, energyF, quantaF, istateF, ivibF, vF, spinF, sigmaF, ilambdaF, omegaF, &
              !$omp&  dimenF, guParity, indSymF, passed, branch, nu_if, lineStr, linestr2, A_einst, boltz_fc, absorption_int, tm) &
              !$omp&  schedule(static) reduction(+:indTrans)
              ! loop over levels in the final state
              loopLevelsF : do indLevelF = 1, nLevelsF

                ! energy and quantum numbers of final state
                energyF =  eigen(indF, indGammaF)%val(indLevelF)
                quantaF => eigen(indF, indGammaF)%quanta(indLevelF)

                istateF  = quantaF%istate  ! electronic state
                ivibF    = quantaF%ivib    ! vibrational (contracted)
                vF      = quantaF%v       ! vibrational
                spinF   = quantaF%spin    ! electron spin
                sigmaF  = quantaF%sigma   ! spin projection
                ilambdaF = quantaF%ilambda ! e- orb. ang. mom. projection
                omegaF  = quantaF%omega   ! tot. ang. mom. proj. mol. ax

                dimenF  = eigen(indF, indgammaF)%Ndimen

                ! reconstruct symmetry for C2v case
                guParity = poten(istateF)%parity%gu
                indSymF  = correlate_to_Cs(indGammaF, guParity)

                ! apply energy filter to final (upper) state
                call energy_filter_ul(jF, energyF, passed, 'upper')
                !
                ! skip if both the upper and lower states are bound states if the unbound filter is on
                !if (intensity%bound_filter.and.intensity%unbound.and.(quantaF%bound.and.quantaI%bound)) passed = .false.
                !
                ! skipping unbound upper states if the bound filter is on
                !if (intensity%bound_filter.and.intensity%bound.and..not.quantaF%bound) passed = .false.
                !
                if ( intensity%bound_filter) then
                   !
                   if ( intensity%unbound.and.(quantaF%bound.and.quantaI%bound) ) passed = .false. 
                   !
                   ! skipping unbound upper states if the bound filter is on                  
                   if( intensity%bound .and. .not.quantaF%bound ) passed = .false.
                   !
                   ! skip if the upper states is bound states for unbound_upper as true 
                   if ( intensity%unbound_upper.and.quantaF%bound ) passed = .false.
                   !
                endif
                !
                if ( .not. passed ) cycle loopLevelsF

                ! apply transition intensity filter, result of which is
                ! overidden by mat. elem. filter if we want mat elems.
                call intens_filter(jI, jF, energyI, energyF, &
                  indSymI, indSymF, iGammaPair, passed)
                !
                if (intensity%use_fitting) then
                  !
                  quantaF => eigen(indF,indGammaF)%quanta(indLevelF)
                  !
                  call transitions_filter_from_fitting(jI,jF,indSymI,indSymF,indLevelI,indLevelF,&
                       energyI,energyF,quantaI,quantaF,passed,iEntry_fitting)

                  !
                endif


                if ( Intensity%matelem ) then
                  call matelem_filter(jI, jF, energyI, energyF, &
                    indSymI, indSymF, iGammaPair, passed)
                endif

                if ( .not. passed ) cycle loopLevelsF

                ! check which OQS branch this transition belongs to
                branch = OQS_branch(jI, jF)

                ! transitions should have energy change > 0
                nu_if = energyF - energyI
                if ( nu_if < small_) cycle

                ! increment transition counter for valid transition
                indTrans = indTrans + 1

                ! vector of basis state coefficients for final state
                vecF(1:dimenF) = eigen(indF, indGammaF)&
                                  %vect(1:dimenF, indLevelF)

                select case ( trim(intensity%action) )

                  case default
                    stop 'only ABSORPTION and TM are properly coded&
                      & at this time'

                  case('ABSORPTION', 'EMISSION')
                    lineStr = ddot(dimenF, halfLineStr, 1, vecF, 1)
                    linestr2 = lineStr**2

                    ! calculate the Einstein A coefficient
                    !A_einst = unitConv * (2.0_rk*jI + 1.0_rk) * linestr2 &
                    !  * (8.0_rk * pi**5 * abs(nu_if)**5) / (5.0_rk * vacPerm * planck)

                    ! calculate the Einstein A coefficient
                    A_einst =  A_coef_s_1 * (2.0_rk*jI + 1.0_rk) * linestr2*abs(nu_if)**5

                    linestr2 = linestr2 * Intensity%gns(indSymI) &
                    ! linestrength times transition degeneracy
                      * (2.0_rk*jI + 1.0_rk) * (2.0_rk*jF + 1.0_rk)

                    if ( trim(Intensity%action) == 'ABSORPTION') then

                      boltz_fc = exp(-(energyI - Intensity%ZPE) * beta)&
                        * (1.0_rk - exp(-abs(nu_if) * beta)) &
                        / (Intensity%part_func * nu_if**2)

                      ! intensity in cm/mol
                      absorption_int = 1.0_rk / (8.0_rk * pi * vellgt) &
                        * Intensity%gns(indSymF)*(2.0_rk*jF + 1.0_rk) &
                        * A_einst * boltz_fc

                    else
                      stop 'EMISSION has not yet been coded for the &
                        &quadrupole operator'

                      !!! this is still dipole equation
                      ! emissivity in ergs/mol/sr
                      boltz_fc = exp(-(energyF - Intensity%ZPE) * beta)

                      absorption_int = emcoef * A_einst * boltz_fc * nu_if &
                        * Intensity%gns(indSymI)*(2.0_rk*jF + 1.0_rk) &
                        / Intensity%part_func

                    endif
                    !
                    if ( linestr2 >= Intensity%threshold%linestrength &
                        .and. absorption_int >= Intensity%threshold%intensity &
                        ) then
                       !
                       if (trim(intensity%linelist_file)/="NONE") then
                          !
                          if ( intensity%matelem ) then 
                            !
                            acoef_RAM(indLevelF) = lineStr 
                            nu_ram(indLevelF) = nu_if
                            indexi_RAM(indLevelF) = quantaI%iJ_ID
                            indexf_RAM(indLevelF) = quantaF%iJ_ID
                            !
                          else ! exomol
                            !
                            acoef_RAM(indLevelF) = A_einst 
                            nu_ram(indLevelF) = nu_if
                            indexi_RAM(indLevelF) = quantaI%iroot
                            indexf_RAM(indLevelF) = quantaF%iroot
                            !
                          endif
                          !
                        else
                          !
                          !$omp critical
                          write(out, &
                            ! Fixed width output
                            "( (f5.1, 1x, a4, 3x),a2, (f5.1, 1x, a4, 3x),&
                            &a1,(2x, f11.4,1x),a2,(1x, f11.4,1x),f11.4,2x,&
                            & 3(1x, es16.8),&
                            & ' ( ',i2,1x,i3,1x,i2,2f8.1,' )',a2,&
                            &'( ',i2,1x,i3,1x,i2,2f8.1,' )')") &
                            jF, sym%label(indSymF), dir, &
                            jI, sym%label(indSymI), branch, &
                            energyF - Intensity%ZPE, dir, &
                            energyI - Intensity%ZPE, nu_if, &
                            linestr2, A_einst, absorption_int, &
                            istateF, vF, ilambdaF, sigmaF, omegaF, dir, &
                            istateI, vI, ilambdaI, sigmaI, omegaI
                            !
                            ! CSV output
                            ! "(a2,',',f5.1,',',a2,',',f5.1,',',a2,',',a1,',',&
                            ! &f11.4,',',f11.4,',',f11.4,',',&
                            ! &es16.8,',',es16.8,',',es16.8,',',&
                            ! &i2,',',i3,',',i2,',',f8.1,',',f8.1,',',&
                            ! &i2,',',i3,',',i2,',',f8.1,',',f8.1,',')")&
                            ! dir, jF, sym%label(indSymF), jI, sym%label(indSymI), branch, &
                            ! energyF - Intensity%ZPE, energyI - Intensity%ZPE, nu_if, &
                            ! linestr2, A_einst, absorption_int, &
                            ! istateF, vF, ilambdaF, sigmaF, omegaF, &
                            ! istateI, vI, ilambdaI, sigmaI, omegaI
                            !
                            if ( trim(Intensity%linelist_file) /= 'NONE') then
                        
                              if ( Intensity%matelem ) then
                        
                                write(richUnit(indI, indF), &
                                  "(i8,i8,2i3,4x,e24.14)") &
                                  quantaI%iJ_ID, quantaF%iJ_ID, 1, 1, lineStr
                        
                              else
                        
                                write(transUnit, &
                                  "(i12,1x,i12,1x,es10.4,4x,f16.6)") &
                                  quantaF%iroot, quantaI%iroot, A_einst, nu_if
                        
                              endif
                            endif
                            !$omp end critical
                       endif
                       !
                       if (intensity%use_fitting .and. iEntry_fitting>0) then
                           !
                           !$omp critical
                           write(intunit, &
                           ! Fixed width output
                           "(a2,1x,f5.1,1x,a2,1x,f5.1,1x,a2,1x,a1,1x,&
                           &f11.4,1x,f11.4,1x,f11.4,1x,&
                           &es16.8,1x,es16.8,1x,es16.8,1x,es16.8,1x,&
                           &i2,1x,i3,1x,i2,1x,f8.1,1x,f8.1,1x,&
                           &i2,1x,i3,1x,i2,1x,f8.1,1x,f8.1,1x)")&
                           !
                           dir, jF, sym%label(indSymF), jI, sym%label(indSymI), branch, &
                           energyF - Intensity%ZPE, energyI - Intensity%ZPE, nu_if, &
                           linestr2, A_einst, absorption_int,fitting%obs(iEntry_fitting)%weight,&
                           istateF, vF, ilambdaF, sigmaF, omegaF, &
                           istateI, vI, ilambdaI, sigmaI, omegaI
                           !
                           !$omp end critical
                           !                     
                       endif                       
                       !
                    endif
                    !
                  case('TM')

                    tm = &
                      dot_product(halfLineStr(1:dimenF), vecF(1:dimenF))

                    lineStr = tm

                    if ( lineStr >= Intensity%threshold%intensity) then
                      !$omp critical
                      write(out, &
                        "( (i4, 1x, a3, 3x),'->', (i4, 1x, a3, 3x),a1, &
                        &(2x, f13.6,1x),'->',(1x, f13.6,1x),f12.6, &
                        &f15.8)") &
                        jI, sym%label(indSymI), jF, sym%label(indSymF),&
                        branch, lineStr, indTrans, tm
                      !$omp end critical
                    endif
                end select
              enddo loopLevelsF
              !$omp enddo

              deallocate(vecF)
              !$omp end parallel
              !
              ! generate the line list (Transition file)
              !
              if (trim(intensity%linelist_file)/="NONE") then
                 !
                 do indLevelF=1,nlevelsF
                   !
                   if (acoef_RAM(indLevelF)>0.0_rk) then 
                      !
                      if ( intensity%matelem ) then 
                        !
                        write(richunit(indI,indF),"(i8,i8,2i3,4x,e24.14)") & 
                                          indexi_RAM(indLevelF),indexf_RAM(indLevelF),1,1,acoef_RAM(indLevelF)
                        !
                      else
                        !
                        write(transunit,"(i12,1x,i12,1x,es10.4,4x,f16.6)") & 
                                  indexf_RAM(indLevelF),indexi_RAM(indLevelF),acoef_RAM(indLevelF),nu_ram(indLevelF)
                      endif
                   endif
                   !
                 enddo
                 !
                 if (allocated(acoef_RAM)) then 
                    deallocate(acoef_RAM)
                    call ArrayStop('swap:acoef_RAM')
                 endif 
                 if (allocated(nu_RAM)) then 
                    deallocate(nu_RAM)
                    call ArrayStop('swap:nu_RAM')
                 endif 
                 if (allocated(indexf_RAM)) then 
                    deallocate(indexf_RAM)
                    call ArrayStop('swap:indexf_RAM')
                 endif 
                 if (allocated(indexi_RAM)) then 
                    deallocate(indexi_RAM)
                    call ArrayStop('swap:indexi_RAM')
                 endif 
                 !
              endif
              !
              if ( iVerbose >= 5 ) call TimerReport
              !
            enddo loopLevelsI
          enddo
        enddo
      enddo

      ! close some RichMol J files
      if ( Intensity%matelem ) then

        do indF = max(1, indI - 1), min(nJ, indI + 1)

          jF = jVal(indF)

          if (     nint(abs(jI - jF)) > 2 &
              .or. nint(jI + jF) == 0) cycle

          if (jI > jF) cycle

          write(richUnit(indI, indF), "('End richmol format')")
          close(richUnit(indI, indF))

        enddo
      endif
    enddo

    deallocate(vecI)
    call ArrayStop('intensity-vectors')

    deallocate(halfLineStr)
    call ArrayStop('halfLineStr')

    if (  trim(intensity%linelist_file) /= 'NONE' ) then
      close(transUnit, status='keep')
    endif
    !
    if (intensity%use_fitting) then
       call IOstop('Selected intensities')
       close(unit = intunit,status="keep")
    endif
    !
    call TimerStop('Intensity calculations')

  end subroutine qm_intensity

  subroutine do_LF_matrix_elements(iLFa, iLFb, iunit, jI, jF, icount)
    ! calculates the components of the matrix element due to the
    ! molecule to laboratory frame transformation (i.e rotational part)
    ! when called by qm_intensity

    implicit none

    ! I/O variables
    real(rk), intent(in)                :: jI, jF
    integer(ik), intent(in)             :: iLFa, iLFb, iunit
    integer(ik), intent(out), optional  :: icount

    ! quantum numbers and indexes
    integer(ik)                         :: indJf, indJi, indMf, indMi
    integer(ik)                         :: Mf_, Mi_
    real(rk)                            :: Mf, Mi

    ! calculation variables
    real(rk)                            :: T(3, 3, 0:2, -2:2)
    integer(ik)                         :: m
    real(rk)                            :: matElem, f3j

    ! assorted
    integer(ik)                         :: icount_

    ! The matrix T gives the transformation from the spherical,
    ! irreducible representation to the Cartesian representation.
    ! Entries corresponding to rank zero and one components are zero
    ! as the quadrupole moment is zero for these in diatomic molecules.
    !
    ! Imaginary values are made real here as this does not change the
    ! value of the line strength, but we must indicate when there is an
    ! additional phase factor in the RichMol file.

    ! real
    T(1, 1, 2, 2) = 0.5_rk
    T(1, 1, 2, 0) = -1.0_rk/sqrt(6.0_rk)
    T(1, 1, 2, -2) = 0.5_rk

    ! real
    T(2, 2, 2, 2) = -0.5_rk
    T(2, 2, 2, 0) = -1.0_rk/sqrt(6.0_rk)
    T(2, 2, 2, -2) = -0.5_rk

    ! real
    T(3, 3, 2, 0) = 2.0_rk/sqrt(6.0_rk)

    ! imaginary
    T(1, 2, 2, 2) = -0.5_rk
    T(1, 2, 2, -2) = 0.5_rk

    ! imaginary
    T(2, 1, 2, 2) = -0.5_rk
    T(2, 1, 2, -2) = 0.5_rk

    ! real
    T(1, 3, 2, 1) = -0.5_rk
    T(1, 3, 2, -1) = 0.5_rk

    ! real
    T(3, 1, 2, 1) = -0.5_rk
    T(3, 1, 2, -1) = 0.5_rk

    ! imaginary
    T(2, 3, 2, 1) = 0.5_rk
    T(2, 3, 2, -1) = 0.5_rk

    ! imaginary
    T(3, 2, 2, 1) = 0.5_rk
    T(3, 2, 2, -1) = 0.5_rk

    ! initialise variables
    icount_ = 0
    matElem = 0
    T = 0
    indJf = nint(2*jF) ! for indexing loops (jF/i may be 1/2 integer)
    indJi = nint(2*jI)

    call TimerStart('do_LF_matrix_elements')

    loop_F :  do indMf = -indJf, indJf, 2

      Mf = 0.5_rk * indMf

      ! remove imaginary factor due to (-)^Mf when jF is half-integer

      ! this method requires only +i phase factor later
      if ( mod(indJf + 1, 2) == 0 ) then
        Mf_ = nint(0.5_rk*indMf - 0.5_rk)
      else
        Mf_ = nint(0.5_rk*indMf)
      endif

      loop_I : do indMi = -indJi, indJi, 2

        Mi = 0.5_rk * indMi

        ! remove imaginary phase due to (-)^Mi when jI is half-integer

        ! this method requires only +i phase factor later
        if ( mod(indJi + 1, 2) == 0 ) then
          Mi_ = nint(0.5_rk*indMi - 0.5_rk)
        else
          Mi_ = nint(0.5_rk*indMi)
        endif

        ! M quantum number selection rule (from 3-j symbol)
        m = nint(Mi - Mf)
        if  ( abs(m) > 2 )  cycle

        ! calculate value of 3j symbol and then LF matrix elem
        f3j = three_j(jI, 2.0_rk, jF, Mi, real(m, rk), -Mf)
        matElem = T(iLFa,iLFb, 2, m) * (-1.0_rk)**Mf * f3j *&
                  sqrt( (2.0_rk*jI + 1.0_rk)*(2.0_rk*jF +  1.0_rk) )

        if ( abs(matElem) > small_ ) then
          icount_ = icount_ +  1
          if (.not. present(icount)) then
            ! write integerised Mi/f and LF matrix elem to file.=
            write (iunit, "(2(i6),3x,e24.14)") Mi_, Mf_, matElem
          endif
        endif
      enddo loop_I
    enddo loop_F

    if ( present(icount) ) icount = icount_

    call TimerStop('do_LF_matrix_elements')

  end subroutine do_LF_matrix_elements

  subroutine do_1st_half_linestrength(jI, jF, indI, indF, dimenI,&
                                      dimenF, vector, half_ls)
    ! calculates the first half of the linestrength when called by
    ! qm_intensity

    implicit none

    ! I/O  variables
    real(rk), intent(in)    :: jI, jF
    integer(ik), intent(in) :: indI, indF, dimenI, dimenF
    real(rk), intent(in)    :: vector(:)
    real(rk), intent(out)   :: half_ls(:)

    ! quantum numbers
    integer(ik)             :: istateI, istateF, ivibI, ivibF, &
                                ilambdaI, ilambdaF
    real(rk)                :: omegaI, omegaF, spinI, spinF, &
                                sigmaI, sigmaF

    ! psuedo quantum numbers and indexes
    integer(ik)             :: iomegaI_, iomegaF_, istateI_, istateF_, &
                                ilambdaI_, ilambdaF_
    real(rk)                :: spinI_, spinF_
    integer(ik)             :: icontrI, icontrF, indQuad, indPermute, &
                                indSigmaV

    ! calculation variables
    real(rk)                :: f3j, vibME, ls
    integer(ik)             :: itau
    !
    call TimerStart('do_1st_half_linestr')
    !
    half_ls = 0
    !
    !$omp parallel do private(&
    !$omp & icontrF,istateF,ivibF,ilambdaF,omegaF,spinF,sigmaF,iomegaF_,&
    !$omp & icontrI,istateI,ivibI,ilambdaI,omegaI,spinI,sigmaI,iomegaI_,f3j,ls,indQuad,indPermute,&
    !$omp & istateI_,ilambdaI_,spinI_,istateF_,ilambdaF_,spinF_,indSigmaV,vibME,itau) shared(half_ls) schedule(guided)
    loop_F : do icontrF = 1, dimenF

      istateF  = basis(indF)%icontr(icontrF)%istate
      ivibF    = basis(indF)%icontr(icontrF)%ivib
      ilambdaF = basis(indF)%icontr(icontrF)%ilambda
      omegaF  = basis(indF)%icontr(icontrF)%omega
      spinF   = basis(indF)%icontr(icontrF)%spin
      sigmaF  = basis(indF)%icontr(icontrF)%sigma

      ! remove imaginary factor (-1)^Omega when J is half-int
      iomegaF_ = nint(omegaF)
      if (mod(nint(2.0_rk*omegaF+1.0_rk),2)==0 ) iomegaF_ = nint((2.0_rk*omegaF-1.0_rk)*0.5_rk)         

      loop_I : do icontrI = 1, dimenI

        istateI  = basis(indI)%icontr(icontrI)%istate
        ivibI    = basis(indI)%icontr(icontrI)%ivib
        ilambdaI = basis(indI)%icontr(icontrI)%ilambda
        omegaI  = basis(indI)%icontr(icontrI)%omega
        spinI   = basis(indI)%icontr(icontrI)%spin
        sigmaI  = basis(indI)%icontr(icontrI)%sigma

        if (abs(nint(omegaF - omegaI))>2.or.nint(spinI-spinF)/=0.or.nint(sigmaI-sigmaF)/=0) cycle loop_I
        if (abs(nint(omegaF - omegaI))/=abs(ilambdaI-ilambdaF)) cycle loop_I
        
        iomegaI_ = int(omegaI)
        if (mod(nint(2.0_rk*omegaI+1.0_rk),2)==0 ) iomegaI_ = nint((2.0_rk*omegaI-1.0_rk)*0.5_rk)

        !remove imaginary factor (-1)^Omega when J is half-int
        iomegaI_ = int(omegaI)
        if (mod(nint(2.0_rk*omegaI+1.0_rk),2)==0 ) iomegaI_ = nint((2.0_rk*omegaI-1.0_rk)*0.5_rk)

        f3j = three_j(jI, 2.0_rk, jF, omegaI, omegaF - omegaI, -omegaF)
        if ( abs(f3j) < Intensity%threshold%coeff ) cycle loop_I

        ls = 0

        loop_quadpole :  do indQuad = 1, nQuadrupoles

          !field  => quadrupoletm(indQuad)

          ! we can calculate opposite matrix elements at same time
          do indPermute =  0, 1

            if  ( indPermute == 0 ) then
              istateI_  = quadrupoletm(indQuad)%istate  ; istateF_  = quadrupoletm(indQuad)%jstate
              ilambdaI_ = quadrupoletm(indQuad)%lambda  ; ilambdaF_ = quadrupoletm(indQuad)%lambdaj
              spinI_   = quadrupoletm(indQuad)%spini   ; spinF_   = quadrupoletm(indQuad)%spinj
            else
              istateI_  = quadrupoletm(indQuad)%jstate  ; istateF_  = quadrupoletm(indQuad)%istate
              ilambdaI_ = quadrupoletm(indQuad)%lambdaj ; ilambdaF_ = quadrupoletm(indQuad)%lambda
              spinI_   = quadrupoletm(indQuad)%spinj   ; spinF_   = quadrupoletm(indQuad)%spini
            endif

            ! do not double count diagonal elements
            if (      indPermute == 1 &
                .and. istateI_  == istateF_ &
                .and. ilambdaI_ == ilambdaF_ &
                .and. nint(spinI_ - spinF_) == 0 &
                ) cycle

            ! check the field I and F states match the basis loop
            if (     istateI_ /= istateI &
                .or. istateF_ /= istateF &
                ) cycle

            ! since only one of Lambda > or < 0 is given as input, we
            ! must also account for negative values. This is done using
            ! a lab-fixed inversion, equivalent to sigmav operation.
            do indSigmaV = 0, 1

              ! skip if both lambda = 0, else double counting
              if ( indSigmaV == 1 &
                  .and. abs(ilambdaF_) + abs(ilambdaI_) == 0 &
                  ) cycle

              ! apply sigmav transformation to build Lambda < 0 states
              ilambdaI_ = ilambdaI_ * (-1)**indSigmaV
              ilambdaF_ = ilambdaF_ * (-1)**indSigmaV

              ! check field states match basis states
              if (     ilambdaI_ /= ilambdaI &
                  .or. ilambdaF_ /= ilambdaF &
                  ) cycle

              ! reinforce lambda selection rule
              if  ( abs(ilambdaF - ilambdaI) > 2 ) cycle

              ! vibrational matrix element
              vibME = quadrupoletm(indQuad)%matelem(ivibI, ivibF)

              ! the Sigma electronic states (with lambda = 0) have
              ! definite +/- parities - thus, if one of the mat.
              ! elem. states is the Sigma^-, lambda = 0 state, the
              ! sign of the matrix element will change under the sigmaV
              ! inversion. This cancels out if both states are ungerade.
              if ( indSigmaV == 1 ) then

                itau = 0

                if (      ilambdaI_ == 0 &
                    .and. poten(istateI)%parity%pm == -1 &
                    ) itau = itau + 1

                if (      ilambdaF_ == 0 &
                    .and. poten(istateF)%parity%pm == -1 &
                    ) itau = itau + 1

                vibME = vibME * (-1.0_rk)**itau
              endif

              ! calculate linestrength and add to sum
              ls =  vibME * f3j * vector(icontrI)
              half_ls(icontrF) = half_ls(icontrF) + (-1.0_rk)**(iomegaI_) * ls

            enddo
          enddo
        enddo loop_quadpole
      enddo  loop_I
    enddo loop_F
    !$omp end parallel do

    call TimerStop('do_1st_half_linestr')

  end subroutine do_1st_half_linestrength



  function OQS_branch(jI,jF) result (X)

    real(rk),intent(in)  :: jI,jF
    character(len=1)        :: X

    select case(trim(intensity%action))

      case default

        X = 'D'

      case ('ABSORPTION','EMISSION')

        if ( nint(jF-jI)==-2 ) then
          X = 'O'
        elseif ( nint(jF-jI)==-1) then
          X = 'P'
        elseif( nint(jF-jI)==0 ) then
          X = 'Q'
        elseif( nint(jF-jI)==1) then
          X = 'R'
        elseif( nint(jF-jI)==2) then
          X = 'S'
        else
          X = 'X'
        endif

    end select

  end function OQS_branch


  subroutine find_igamma_pair(igamma_pair)

    integer(ik),intent(out) :: igamma_pair(:)
    integer(ik)     :: igammaI,ngamma

    if ( trim(intensity%action) == 'TM' ) then
      igamma_pair = 1
      return
    endif

    do igammaI = 1,sym%Nrepresen

      ! count number of hits
      ngamma = 0
      igamma_pair(igammaI) = igammaI

      ! do igammaF = 1,sym%Nrepresen
      !   !
      !   if (igammaI/=igammaF.and.intensity%isym_pairs(igammaI)==intensity%isym_pairs(igammaF)) then 
      !     !
      !     igamma_pair(igammaI) = igammaF
      !     !
      !     ngamma = ngamma + 1 
      !     !
      !     if (ngamma>1) then 
      !       !
      !       write(out,"('qm_intensity: Assumption that selection rules come in pairs is wrong!')")
      !       stop 'qm_intensity: Assumption that all selection rules work in pairs is wrong!'
      !       !
      !     endif   
      !     !
      !   endif
      !   !
      ! enddo

      if ( nint(intensity%gns(igammaI) &
          - intensity%gns(igamma_pair(igammaI))) /= 0 &
         ) then
        write(out, &
          "('qm_intensity: selection rules do not agree with Gns')")
        stop 'qm_intensity: selection rules do not agree with Gns!'
      endif
    enddo

  end subroutine find_igamma_pair


  subroutine energy_filter_ul(J, energy, passed, upp_low)

    real(rk), intent(in)         :: J
    real(rk), intent(in)         :: energy
    character(len=5), intent(in) :: upp_low
    logical, intent(out)         :: passed
    real(rk)                     :: bottom, top

    passed = .false.

    select case (upp_low)
      case ('lower')
        bottom = Intensity%erange_low(1)
        top = Intensity%erange_low(2)
      case ('upper')
        bottom = Intensity%erange_upp(1)
        top = Intensity%erange_upp(2)
    end select

    if (      J >= intensity%J(1) &
        .and. J <= intensity%J(2) &
        .and. energy-intensity%ZPE >= bottom &
        .and. energy-intensity%ZPE <= top &
       ) then

      passed = .true.

    endif
  end subroutine energy_filter_ul


  subroutine intens_filter(jI, jF, energyI, energyF, isymI, isymF, &
                           igamma_pair, passed)

    real(rk), intent(in)    :: jI, jF
    integer(ik), intent(in) :: isymI, isymF
    real(rk), intent(in)    :: energyI, energyF
    integer(ik), intent(in) :: igamma_pair(sym%Nrepresen)
    real(rk)                :: nu_if
    logical, intent(out)    :: passed
    real(rk), dimension(4)  :: local_sym

    passed = .false.
    local_sym = Intensity%isym_pairs
    nu_if = energyF - energyI

    if ( &
        ! nuclear stat.weight:
        intensity%gns(isymI) > small_ &

        ! absorption/emission go only in one direction
        .and. nu_if > small_ &

        ! spectroscopic window
        .and. nu_if >= intensity%freq_window(1) &
        .and. nu_if <= intensity%freq_window(2) &

        .and. jI >= intensity%J(1) &
        .and. jI <= intensity%J(2) &

        .and. jF >= intensity%J(1) &
        .and. jF <= intensity%J(2) &

        .and. energyI - intensity%ZPE >= intensity%erange_low(1) &
        .and. energyI - intensity%ZPE <= intensity%erange_low(2) &

        .and. energyF - intensity%ZPE >= intensity%erange_upp(1) &
        .and. energyF - intensity%ZPE <= intensity%erange_upp(2) &

    ) then
      passed = .true.
    endif
    !
    if (     trim(intensity%action) == 'ABSORPTION' &
        .or. trim(intensity%action) == 'EMISSION' &
    ) then

      ! In order to avoid double counting of transitions we exclude
      ! jI=jF==intensity%J(2), i.e. Q branch for the highest J is
      ! never considered:
      
      ! set passed
      passed = passed &
        .and. &
        ( nint(jF - intensity%J(1)) /= 0 &
            .or. &
          nint(jI - intensity%J(1)) /= 0 &
            .or. &
          nint(jI + jF) == 1 &  ! check ????
        ) &
        .and. &
        (intensity%J(1) + intensity%J(2) > 0) &
        .and. &

        ! selection rules:
        intensity%isym_pairs(isymI) == intensity%isym_pairs(isymF) &
        .and. &
        igamma_pair(isymI) == isymF &
        .and. &

        ! selection rules from the 3j-symbols
        abs(nint(jI - jF)) <= 2 &
        .and. &
        nint(jI + jF) >= 2
      ! end set passed
    endif
  end subroutine intens_filter


  subroutine matelem_filter(jI,jF,energyI,energyF,isymI,isymF,igamma_pair,passed)

    real(rk), intent(in)    :: jI, jF
    integer(ik), intent(in) :: isymI, isymF
    real(rk), intent(in)    :: energyI, energyF
    integer(ik), intent(in) :: igamma_pair(sym%Nrepresen)
    real(rk)                :: nu_if
    logical, intent(out)    :: passed

    passed = .false.

    nu_if = energyF - energyI

    if ( &   ! nuclear stat.weight:
             intensity%gns(isymI) > small_ &
       .and. &
             ! absorption/emission go only in one direction
             (Jf-Ji) > -small_ &
       .and. &
             ! spectroscopic window
             nu_if >=intensity%freq_window(1) &
       .and. nu_if <=intensity%freq_window(2) &
       .and. jI >= intensity%J(1) &
       .and. jI <= intensity%J(2) &
       .and. jF >= intensity%J(1) &
       .and. jF <= intensity%J(2) &
       .and. energyI - intensity%ZPE >= intensity%erange_low(1) &
       .and. energyI - intensity%ZPE <= intensity%erange_low(2) &
       .and. energyF - intensity%ZPE >= intensity%erange_upp(1) &
       .and. energyF - intensity%ZPE <= intensity%erange_upp(2) &
    ) then
      passed = .true.
    endif

    if (     trim(intensity%action) == 'ABSORPTION' &
        .or. trim(intensity%action) == 'EMISSION') then

      ! In order to avoid double counting of transitions we exclude
      ! jI=jF==intensity%J(2), i.e. Q branch for the highest J is never considered:

      passed = passed &
          .and. &
        (nint(jF - intensity%J(1)) /= 0 &
            .or. &
         nint(jI - intensity%J(1)) /= 0 &
            .or. &
         nint(jI + jF) == 1 &
        ) &
          .and. &
        intensity%J(1) + intensity%J(2) > 0 &
          .and. &

        ! selection rules:
        intensity%isym_pairs(isymI) == intensity%isym_pairs(isymF) &
          .and. &
        igamma_pair(isymI) == isymF &
          .and. &
        ! selection rules from the 3j-symbols
        abs(nint(jI - jF)) <= 2 &
          .and. &
        nint(jI+jF) >= 2
    endif
  end subroutine matelem_filter


  subroutine Sort_levels(iverbose,njval, jval)

    integer(ik), intent(in) :: iverbose,njval
    real(rk), intent(in)    :: jval(njval)
    integer(ik)             :: jind, nroots, nlevels,  iroot, ilevel, &
                               jlevel, info, itau, jtau
    real(rk)                :: energy
    logical                 :: passed
    integer(ik)             :: iind, jroot, igamma

    if (iverbose >= 2) then
      write(out,"(/'Read and sort eigenvalues in increasing order...')")
    endif

    call TimerStart('Sort eigenvalues')

    if ( .not. allocated(basis) ) then
      stop 'Sort_levels error: associated(basis) = .false.'
    endif

    ! In practice we do not need all stored rootes, but only those that
    ! pass the filters given in the input. Here we perform a
    ! pre-selection of the states, which are actually required in the
    ! calculations. Towards this we read the unit twice: a) to count and
    ! b) to read. After reading all states will be combined together and
    ! sorted wrt energy increasing to produce just one list of states.

    ! Estimate the maximal number of energy records
    Nlevels = 0
    do jind = 1, njval
      do itau = 1,sym%NrepresCs
         Nlevels = max(eigen(jind,itau)%Nlevels,Nlevels)
      enddo
    enddo

    nroots  = 0

    do jind = 1, njval
      do itau = 1,sym%NrepresCs

        Nlevels = eigen(jind,itau)%Nlevels

        do ilevel = 1,Nlevels
          energy = eigen(jind,itau)%val(ilevel)
          igamma =  eigen(jind,itau)%quanta(ilevel)%igamma

          if (      job%ZPE < 0 &
              .and. jind == 1 &
              .and. ilevel == 1 &
          ) then
            job%zpe = energy
          endif

          call Energy_filter(Jval(jind),energy,igamma,passed)

          if (passed) then
            nroots = nroots + 1
          endif
        enddo
      enddo
    enddo

    if (nroots==0) then
      write(out,"('Sort_levels: the filters are too tight: no entry')")
    endif

    ! total number of levels to be processed, a global variable
    Neigenlevels = nroots

    return

    allocate(Elevel(Neigenlevels),stat = info)

    if (info /= 0) then
      stop 'Sort_levels allocation error: eigen - out of memory'
    endif

    if (iverbose >= 4) then
      write(out, &
        "('   number of selected eigenvectors: ',i8)") Neigenlevels
    end if

    ! Now we can read and store the energies and their description.
    iroot  = 0
    do jind = 1, njval
      do itau = 1,sym%NrepresCs
        Nlevels = eigen(jind,itau)%Nlevels
        do ilevel = 1,Nlevels
          iroot  = iroot + 1
          Elevel(iroot)%jind   = jind
          Elevel(iroot)%igamma = itau
          Elevel(iroot)%ilevel = ilevel
        enddo
      enddo
    enddo

    if (iroot/=Neigenlevels) then
      write(out,&
        '("The number of selected levels ",i0," does&
        & not agree with Neigenlevels =  ",i0)') ilevel, Neigenlevels
      stop 'iroot/=Neigenlevels'
    endif

    ! Sort according the energy increasing
    do iroot =1,Neigenlevels

      ilevel = Elevel(iroot)%ilevel
      iind = Elevel(iroot)%jind
      itau = Elevel(iroot)%igamma

      energy = eigen(iind,itau)%val(ilevel)

      do jroot =iroot+1,Neigenlevels

        jind   = Elevel(jroot)%jind
        jtau   = Elevel(jroot)%igamma
        jlevel = Elevel(jroot)%ilevel

        if (eigen(jind,jtau)%val(jlevel) < energy) then

          Elevel(jroot)%jind = iind
          Elevel(jroot)%ilevel = ilevel
          Elevel(jroot)%igamma = itau

          Elevel(iroot)%jind= jind
          Elevel(iroot)%ilevel = jlevel
          Elevel(jroot)%igamma = jtau

          ilevel = Elevel(iroot)%ilevel
          iind   = Elevel(iroot)%jind
          itau   = Elevel(iroot)%igamma
          energy = eigen(iind,itau)%val(ilevel)

        endif
      enddo
    enddo

    if (iverbose >= 2) write(out,"('...done!')")

    call TimerStop('Sort eigenvalues')
    !
  end subroutine Sort_levels


  subroutine Energy_filter(Jval,energy,igamma,passed)

    real(rk), intent(in)    :: Jval, energy
    integer(ik), intent(in) :: igamma
    logical, intent(out)    :: passed
    real(rk)                :: erange(2)

    passed = .false.
    erange(1) = min(intensity%erange_low(1),intensity%erange_upp(1))
    erange(2) = max(intensity%erange_low(2),intensity%erange_upp(2))

    if (      job%isym_do(igamma) &
        .and. energy - job%ZPE >= erange(1) &
        .and. Jval >= intensity%J(1) &
        .and. Jval <= intensity%J(2) &
        .and. energy - job%ZPE <= erange(2) &
    ) then

      passed = .true.

    endif

  end subroutine Energy_filter


end module quadrupole
