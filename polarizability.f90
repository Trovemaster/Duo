module polarizability
  !
  use accuracy,      only : hik, ik, rk, ark, cl, out, vellgt, planck, &
                            avogno, boltz, pi, small_
  use diatom_module, only : job, Intensity, quantaT, eigen, basis, Ndipoles, Nss, &
                            dipoletm, duo_j0, fieldT, poten, three_j, &
                            jmin_global,spinspin
  use timer,         only : IOstart, Arraystart, Arraystop, ArrayMinus, &
                            Timerstart, Timerstop, MemoryReport, &
                            TimerReport, memory_limit, memory_now
  use symmetry,      only : sym,correlate_to_Cs
  
  private
  public pol_tranint

  real(rk), allocatable, save :: polarisability_me(:, :)

  type DmatT  
    real(rk), pointer   :: mat(:, :, :)
  end type DmatT
  
  type DimatT
    integer(ik),pointer :: imat(:,:,:)
  end type DimatT

  type DkmatT
    integer(ik),pointer :: kmat(:,:)
  end type DkmatT
  
  type ElevelT
    integer(ik)  :: jind  ! J index
    integer(ik)  :: igamma
    integer(ik)  :: ilevel
  end type ElevelT
  
  type dipoleT                      
    real(rk),pointer    :: rot(:,:,:,:,:)
  end type dipoleT

  integer(ik)           :: Neigenlevels
  
  type(ElevelT), allocatable  :: Elevel(:)
  
  type(dipoleT), allocatable  :: wigner(:,:) ! rot. component of matrix elems

  contains
  
  subroutine pol_tranint

    ! Obtain the partition function, polarisability transition moments,
    !   linestrengths, and intensities.
    ! Transition moments and intensities are computed by the call to 
    !   pol_intensity subroutine.

    integer(ik)             :: info
    integer(ik)             :: nJ, jind
    real(rk), allocatable   :: Jval(:), q_part(:,:)
    real(rk)                :: Jmin, Jmax, Jval_min, Jval_, energy, &
                               beta, exp_en, part
    integer(ik)             :: ilevel, igamma, irrep, istate, parity_gu, isym
    integer(ik)             :: iverbose = 4
    
    ! initialise array of J values
    ! J = 0 always present, regardless of whether it is used in calculations
    !
    Jmin = minval(intensity%J(1:2))
    Jmax = maxval(intensity%J(1:2))
    !
    Jmin = max(jmin_global, Jmin)
    !
    ! we now count J states, starting from the lowest J = 0 or 0.5, to make 
    !   sure all J > 0 are included in the line list for any intensity%J to
    !   guarantee consistent energy list numbering
    !
    ! find lowest possible J value
    Jval_min = 0
    if (mod(nint(2.0_rk*job%j_list(1)),2) /= 0) Jval_min = 0.5_rk
    Jval_min = max(jmin_global, Jval_min)
    !
    Jval_ = Jval_min
    nJ = 1
    !
    ! count the number of J states
    do while (Jval_ < Jmax)
      Jval_ = Jval_ + 1.0_rk
      nJ = nJ + 1
    enddo
    !
    allocate(Jval(nJ), stat = info)
    if (info /= 0) stop 'pol_tranint allocation error: Jval - out of memory'
    !
    allocate(q_part(20, nJ), stat = info)
    if (info /= 0) stop 'pol_tranint allocation error: q_part - out of memory'
    !
    Jval_ = Jval_min
    jind = 1
    Jval(jind) = Jval_
    do while (Jval_ < Jmax)
      jind = jind + 1
      Jval_ = Jval_ + 1.0_rk
      Jval(jind) = Jval_
    enddo
    !
    call duo_j0(iverbose, Jval)
    !
    ! include ZPE if not provided as part of intensity input
    !
    if (job%shift_to_zpe) then
      !
      do jind = 1, nJ
        !
        Jval_ = Jval(jind)
        !
        do igamma = 1, sym%NrepresCs
          !
          do ilevel = 1, eigen(jind, igamma)%Nlevels
            !
            energy = eigen(jind, igamma)%val(ilevel)
            !
            intensity%zpe = min(intensity%zpe, energy)
            !
          enddo
          !
        enddo
      enddo
      !
      if (iverbose >= 4) &
        write (out,"(/'Partition function = ',f18.4,' T = ',f12.2)") intensity%part_func,intensity%temperature
      !
    endif
    !
    select case (trim(intensity%action))
    !
      case('ABSORPTION', 'EMISSION', 'TM')
        !
        call Sort_levels(iverbose, nJ, Jval(1:nJ))
        !
        ! run intensity simulations
        !
        beta = planck * vellgt / (boltz * intensity%temperature)
        !
        ! compute partition function if not given
        if (intensity%part_func < small_) then
          !
          intensity%part_func = 0
          !
          do jind = 1, nJ
            !
            Jval_ = Jval(jind)
            !
            do igamma = 1, sym%NrepresCs
              !
              do ilevel = 1, eigen(jind, igamma)%Nlevels
                !
                energy = eigen(jind, igamma)%val(ilevel)
                irrep = eigen(jind, igamma)%quanta(ilevel)%igamma
                !
                ! for homonuclear symmetries Nrepres = 4 and the irrep can be reconstucted baed on the parity and g/u: 
                istate  = eigen(jind,igamma)%quanta(ilevel)%istate
                parity_gu = poten(istate)%parity%gu
                isym = correlate_to_Cs(igamma,parity_gu)
                !
                exp_en = exp(-(energy - intensity%ZPE) * beta)
                !
                intensity%part_func = intensity%part_func + intensity%gns(isym)*(2.0_rk*Jval_+1.0_rk) * exp_en
                !
              enddo
              !
            enddo
            !
          enddo
          !
          if (iverbose >= 4) write(out, "(/'Partition function = ', f18.4,'T = ', f12.2)") &
                                   intensity%part_func, intensity%temperature
        endif
        !
        call pol_intensity(Jval, iverbose)
        !
        write(out, '(/a)') 'done'
        !
      ! compute parition function
      case('PARTFUNC')
        !
        write(out, '(/a)') 'compute parition function'
        !
        beta = planck * vellgt / (boltz * intensity%temperature)
        !
        q_part = 0
        !
        ! loop over J values
        do jind = 1, nJ
          !
          Jval_ = Jval(jind)
          !
          do igamma = 1, sym%NrepresCs
            !
            do ilevel = 1, eigen(jind, igamma)%Nlevels
              !
              energy = eigen(jind, igamma)%val(ilevel)
              irrep = eigen(jind, igamma)%quanta(ilevel)%igamma
              !
              ! for homonuclear, symmetries Nrepres = 4, and the irrep can be 
              !   reconstructed based on the parity and g/u:
              istate = eigen(jind, igamma)%quanta(ilevel)%istate
              parity_gu = poten(istate)%parity%gu
              isym = correlate_to_Cs(igamma, parity_gu)
              !
              exp_en = exp(-(energy - intensity%ZPE) * beta)
              !
              q_part(irrep, jind) = q_part(irrep, jind)&
                  + intensity%gns(isym)*(2.0_rk*Jval_+1.0_rk) * exp_en
              !
            enddo
            !
          enddo
        enddo
        !
        part = sum(q_part(:,:))
        !
        do jind = 1, nJ
          do irrep = 1, sym%NrepresCs
            write(out, '(i4,1x,f18.1,1x,es16.8)') &
                irrep, Jval(jind), q_part(irrep, jind)
          enddo
        enddo
        !
        write(out, '(es16.8)') part
        !
    end select
    !
    call MemoryReport
    !
    call TimerReport
    !
  end subroutine pol_tranint
  !
  !
  !
  function cg(j0, k0, dj, dk)
    !
    ! Returns values for some Clebsch-Gordan coefficients
    ! 
    ! Inputs:
    !   j0, k0 are j" and k" (initial state)
    !   dj, dk such that j0 + dj and k0 + dk are j' and k' (final state)
    ! Returns:
    !   cg = <j0 k0, 1 dk| j0 + dj j0 + dk>
    !
    integer(ik), intent(in) :: j0, k0, dj, dk
    real(rk)                :: cg
    real(rk)                :: j, m

    j = real(j0, kind = rk)
    m = real(k0 + dk, kind = rk)

    cg = 0.0_rk

    if (dj ==  0 .and. dk ==  0) cg = m / sqrt(j * (j + 1.0_rk))

    if (dj ==  0 .and. dk ==  1) cg = -sqrt((j + m) * (j - m + 1.0_rk) / (2.0_rk * j * (j + 1.0_rk)))

    if (dj ==  0 .and. dk == -1) cg = sqrt((j - m) * (j + m + 1.0_rk) / (2.0_rk * j * (j + 1.0_rk)))

    if (dj ==  1 .and. dk ==  0) cg = sqrt((j - m + 1.0_rk) * (j + m + 1.0_rk) / ((2.0_rk * j + 1.0_rk) * (j + 1.0_rk)))

    if (dj ==  1 .and. dk ==  1) cg = sqrt((j + m) * (j + m + 1.0_rk) / ((2.0_rk * j + 1.0_rk) * (2.0_rk * j + 2.0_rk)))

    if (dj ==  1 .and. dk == -1) cg = sqrt((j - m) * (j - m + 1.0_rk) / ((2.0_rk * j + 1.0_rk) * (2.0_rk * j + 2.0_rk)))

    if (dj == -1 .and. dk ==  0) cg =-sqrt((j - m) * (j + m) / (j * (2.0_rk * j + 1.0_rk)))

    if (dj == -1 .and. dk ==  1) cg = sqrt((j - m) * (j - m + 1.0_rk) / (2.0_rk * j * (2.0_rk * j + 1.0_rk)))

    if (dj == -1 .and. dk == -1) cg = sqrt((j + m + 1.0_rk) * (j + m) / (2.0_rk * j * (2.0_rk * j + 1.0_rk)))

    return

  end function cg
  !
  !
  !
  subroutine pol_intensity(Jval, iverbose)
    ! Calculates the polarisability transition moments and intensities
    !
    implicit none
    !
    real(rk),intent(in)  :: Jval(:)
    integer(ik),intent(in)   :: iverbose
    !
    integer(ik)    :: nJ, dimenmax
    integer(ik)    :: ilevelI, ilevelF
    integer(ik)    :: nlevelsG(sym%Nrepresen)
    integer(ik)    :: info, indI, indF, itransit, Ntransit, Nrepresen
    integer(ik)    :: igammaI, igammaF
    integer(ik)    :: dimenI, dimenF, nmax, parity_gu, isymI, isymF
    real(rk)       :: energyI, energyF, nu_if, linestr, linestr_zero, &
                      ener_, linestr2, linestr2_zero
    real(rk)       :: tm, jI, jF, ddot
    logical        :: passed, passed_
    !
    real(rk), allocatable     :: vecI(:), vecF(:)
    real(rk), allocatable     :: half_linestr(:), half_linestr_zero(:)
    !
    integer(ik)  :: jind, nlevels, j1, j2
    !
    integer(ik)  :: iroot, NlevelsI, NlevelsF, nlower, k, k_, iLF, iflag_rich
    !
    integer(ik)  :: igamma_pair(sym%Nrepresen), igamma, istateI, istateF,&
                    ivibI, ivibF, ivI, ivF, ilambdaI, ilambdaF, iparityI, itau
    integer(ik)  :: ivF_, ilambdaF_
    real(rk)     :: spinI, spinF, omegaI, omegaF, sigmaI, sigmaF,& 
                    sigmaF_, omegaF_, spinF_
    integer(hik) :: matsize
    !
    character(len=1)  :: branch, ef, pm
    character(len=2)  :: dir
    character(len=10) :: statename
    !
    type(quantaT),pointer  :: quantaI, quantaF
    !
    real(rk)     :: boltz_fc, beta, intens_cm_mol, emcoef,&
                    A_coef_s_1, A_einst, absorption_int, lande
    !
    character(len=130) :: my_fmt !format for I/O specification
    !
    integer           :: ndecimals
    integer(ik)       :: enunit, transunit
    character(len=cl) :: filename, ioname
    !
    logical     :: integer_spin = .true.
    !
    integer(ik) :: alloc_p
    !
    integer(ik) :: Jmax_, ID_J, nMs
    real(rk)    :: J_
    !
    character(len=12)        :: char_Jf,char_Ji,char_LF
    integer(ik), allocatable :: richunit(:,:)
    character(len=2)         :: let_LF ! richmol letters x,y,z
    !
    call TimerStart('Intensity calculations')
    !
    if (sym%maxdegen > 2) then
      !
      write(out, "('pol_intensity: this procedure has not been tested for symmetries with degeneracy higher than 2')")
      !
    endif
    !
    Nrepresen = sym%NrepresCs
    !
    ! constants that will be needed
    beta = planck * vellgt / (boltz * intensity%temperature)
    intens_cm_mol = 8.0d-36 * pi**3 / (3.0_rk * planck * vellgt)
    emcoef = planck*vellgt / (4.0_rk*pi)
    A_coef_s_1 = 64.0d-36 * pi**4 / (3.0_rk * planck)
    !
    nJ = size(Jval)
    !
    ! open units for the line list in the Exomol format
    ! the next section establishes the files to which the line intensities 
    !   and energies will be written
    if (trim(intensity%linelist_file)/="NONE") then
      !
      ! open the '.states' file, which contains the energy levels themselves
      filename = trim(intensity%linelist_file)//'.states'
      write(ioname, '(a, i4)') 'Energy file' ! heading
      call IOstart(trim(ioname), enunit)
      open(unit = enunit, action = 'write', &
          status = 'replace', file = filename)
      !
      ! open the '.tans' file, which contains the positions and A coeffs. 
      !   of the spectral lines
      filename = trim(intensity%linelist_file)//'.trans'
      write(ioname, '(a, i4)') 'Transition file' ! heading
      call IOstart(trim(ioname), transunit)
      open(unit = transunit, action = 'write', &
          status = 'replace', file = filename)
      !
      if (intensity%matelem) then
        !
        Jmax_ = nint(maxval(Jval(:)))
        !
        allocate(richunit(nJ, nJ)) ! used to mark I/O unit
        !
        ! we now begin writing the matrix elements to the file
        do indI = 1, nJ
          !
          jI = Jval(indI)
          !
          write(char_jI, '(i12)') nint(jI)
          !
          do indF = 1, nJ
            !
            jF = Jval(indF)
            !
            if (Jf<Ji) cycle 
            !
            ! selection rules for SECOND rank spherical component
            if(nint(abs(jI - jF)) > 2) cycle
            !
            write(char_jF, '(i12)') nint(jF)
            !
            ! New Richmol format - one file for all components
            !
            ! open the file for the matrix elements of the i-th lab component
            filename = &
            "matelem_ALPHA"//"_j"//trim(adjustl(char_jI))//"_j"//trim(adjustl(char_jF))//"_"//trim(intensity%linelist_file)//".rchm"
            !
            call IOstart(trim(filename), richunit(indI, indF))
            open(unit = richunit(indI, indF), action = 'write', &
                 status='replace', file=filename)
            !
            write(richunit(indI, indF), "('Start richmol format')")
            !
            ! 2nd rank tensor, with 6 non-zero components
            write(richunit(indI, indF), "('ALPHA','   2','   6')")
            write(richunit(indI, indF), "('M-tensor')")
            !
            ! Nine lab-frame polarisability tensor components
            do iLF = 1,6
              !
              let_LF = "xx"
              if (iLF == 2) let_LF = "xy" 
              if (iLF == 3) let_LF = "xz"
              if (iLF == 4) let_LF = "yy"
              if (iLF == 5) let_LF = "yz"
              if (iLF == 6) let_LF = "zz"
              if (iLF == 7) let_LF = "yx"
              if (iLF == 8) let_LF = "zx"
              if (iLF == 9) let_LF = "zy"
              !
              ! flag is -1 when sum over alpha^k_m has prefactor of i
              iflag_rich = 0
              !if (iLF == 2 .or. iLF == 5) iflag_rich = -1
              if (iLF == 2 .or. iLF == 5 .or. iLF == 7 .or. iLF == 9) iflag_rich = -1
              !
              write(char_LF, '(i12)') iLF
              !
              write(richunit(indI, indF), "('alpha',i5,i3,1x,a2)") &
                iLF, iflag_rich, let_LF
              !
              ! call subroutine to calculate LF portion of matrix elements
              call do_LF_matrix_elements(iLF, richunit(indI, indF), jI, jF)
            enddo
            !
            write(richunit(indI, indF), "('K-tensor')")
            !
          enddo
        enddo
      endif
      !
    endif
    !
    ! maximal size of basis functions
    !
    dimenmax = 0
    !
    ! loop over J values
    !
    do jind = 1, nJ
      do igamma = 1, Nrepresen
        !
        ! Estimate the maximal size of basis functions
        dimenmax = max(dimenmax, eigen(jind, igamma)%Ndimen)
        !
      enddo
    enddo
    !
    Nmax = nint(Jval(nJ)) + 1
    !
    ! We now count the number of non-zero transitions, this will help keep
    !   keep track of the calculation process.
    !
    Ntransit = 0
    !
    ! number of inital states
    nlevels = Neigenlevels
    !
    ! For a given symmetry (igamma) with some gns(igamma), we find it's 
    !  counterpart jgamma /= igamma having the same gns(jgamma). We assume 
    !  that there is only one such pair in the case of absorption or emission
    !
    call find_igamma_pair(igamma_pair)
    !
    call TimerStart('Intens_Filter-1')
    !
    ener_ = 0
    nlower = 0
    iroot = 0
    !
    ! loop over initial states
    do indI = 1, nJ
      !
      ! rotational quantum number
      jI = Jval(indI)
      !
      ! loop for each symmetry
      do igammaI = 1, Nrepresen
        !
        nlevelsI = eigen(indI, igammaI)%Nlevels
        !
        ! loop over energy levels in initial J state
        do ilevelI = 1, nlevelsI
          !
          ! energy and quanta of initial state
          energyI = eigen(indI, igammaI)%val(ilevelI)
          istateI = eigen(indI, igammaI)%quanta(ilevelI)%istate
          parity_gu = poten(istateI)%parity%gu
          isymI = correlate_to_Cs(igammaI, parity_gu)
          !
          call energy_filter_lower(jI, energyI, passed)
          !
          if (.not. passed) cycle
          !
          nlower = nlower + 1
          !
          ! loop over final states
          do indF = 1, nJ
            !
            jF = jval(indF)
            !
            do igammaF = 1, Nrepresen
              !
              nlevelsF = eigen(indF, igammaF)%Nlevels
              !
              ! loop over energy levels in final J state
              do ilevelF = 1, nlevelsF
                !
                ! energy and quanta of final state
                energyF = eigen(indF, igammaF)%val(ilevelF)
                istateF = eigen(indF, igammaF)%quanta(ilevelF)%istate
                parity_gu = poten(istateF)%parity%gu
                isymF = correlate_to_Cs(igammaF, parity_gu)
                !
                call intens_filter(jI, jF, energyI, energyF, isymI, &
                                   isymF, igamma_pair, passed)
                !
                if (intensity%matelem) call matelem_filter(jI, jF, energyI,&
                    energyF, isymI, isymF, igamma_pair, passed)
                !
                ! if transition is non-zero then increment counter
                if (passed) then
                  !
                  Ntransit = Ntransit + 1
                  !
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo ! end of transition counting
    !
    call TimerStop('Intens_Filter-1')
    !
    ! loop over final states -> count states for each symmetry
    !
    nlevelsG = 0
    !
    ! guess if half-integer spin case
    if (mod(eigen(1, 1)%quanta(1)%imulti, 2) == 0) integer_spin = .false.
    !
    allocate(vecI(dimenmax), stat = info)
    call ArrayStart('intensity-vecI', info, size(vecI), kind(vecI))
    !
    do indI = 1, nJ
      !
      ! rotational quantum number
      jI = jval(indI)
      J_ = -1 ! for Richmol count
      !
      do igammaI = 1, Nrepresen
        !
        nlevelsI = eigen(indI, igammaI)%Nlevels
        dimenI = eigen(indI, igammaI)%Ndimen
        !
        do ilevelI = 1, nlevelsI
          !
          ! energy and symmetry of the state
          energyI = eigen(indI, igammaI)%val(ilevelI)
          !
          istateI = eigen(indI, igammaI)%quanta(ilevelI)%istate
          parity_gu = poten(istateI)%parity%gu 
          ! C2v/Cs symmetry
          isymI = correlate_to_Cs(igammaI, parity_gu)
          !
          ! ignore states with zero nuclear weight
          if (intensity%gns(isymI) < small_) cycle
          !
          iroot = iroot + 1
          eigen(indI, igammaI)%quanta(ilevelI)%iroot = iroot
          !
          if (trim(intensity%linelist_file) /= "NONE") then
            !
            ! dimension of the basis for initial states
            !
            ! energy, quanta, and degeneracy order of the initial state
            quantaI => eigen(indI, igammaI)%quanta(ilevelI)
            ivibI     = quantaI%ivib
            ivI       = quantaI%v
            sigmaI    = quantaI%sigma 
            spinI     = quantaI%spin 
            ilambdaI  = quantaI%ilambda
            omegaI    = quantaI%omega
            iparityI  = quantaI%iparity 
            statename = trim(quantaI%name)
            !
            ! reconstruct +/- and e/f parities
            pm = "+" ; if (iparityI == 1) pm = "-"
            ef = "e"
            !
            if (mod(nint(2.0 * jI), 2) == 1) then
              itau = mod(nint(jI - 0.5), 2)
            else
              itau = mod(nint(jI), 2)
            endif
            !
            if (itau /= iparityI) then
              ef = "f"
            endif
            !
            ! the variable 'ndecimals' gives the number of decimal digits to
            !   print the values of the energy levels to.
            ! we use 6 decimals for energy levels up to 100,000 cm-1 before
            !   sacrificing decimals to fit larger values in 12 spaces.
            ! this format works for energy levels larger than -10,000 cm-1 and
            !   less than 1e11 cm-1 - Lorenzo Lodi
            !
            ndecimals = 6  - max(0, int(log10(abs(energyI - intensity%ZPE) + 1.d-6) - 4))
            !
            if ( intensity%matelem ) then
              !
              ndecimals=6-max(0, int( log10(abs(energyI-intensity%ZPE)+1.d-6)-4) )
              !
              if (nint(2*Ji)/=nint(2*J_)) then
                !
                J_ = Ji
                ID_J = 0
                !
              endif
              !
              ID_J = ID_J + 1
              quantaI%iJ_ID = ID_J
              !
              if (integer_spin) then 
                !
                write(my_fmt,'(a)') "(i6,1x,i8,1x,i2,1x,i2,3x,e21.14,5x,a4,i3,1x,a2,i4,1x,a2,f8.4,1x,i6,1x,i6,1x,i4,1x,i6,1x,a1,1x,a10)"
                write(enunit,my_fmt) & 
                          nint(J_),ID_J,iparityI+1,1,energyI-intensity%ZPE,'tau:',iparityI,'j:',nint(J_),'c',1.000_rk,nint((omegaI)),&
                          ivI,(ilambdaI),nint((sigmaI)),pm,statename
                !
              else
                !
                !stop 'not tested'
                !
                write(my_fmt,'(A,i0,a)') "(i7,1x,i12,1x,i1,1x,i2,1x,f12.",ndecimals,",1x,f7.1,1x,i6,1x,i4,1x,f7.1,1x,a1,1x,a10)"
                write(enunit,my_fmt) & 
                          int(J_),ID_J,iparityI+1,1,energyI-intensity%ZPE,omegaI,&
                          ivI,(ilambdaI),sigmaI,pm,statename
                          !
              endif
              !
            else
              !
              ndecimals=6-max(0, int( log10(abs(energyI-intensity%ZPE)+1.d-6)-4) )
              if (integer_spin) then 
                !
                write(my_fmt,'(A,i0,a)') "(i12,1x,f12.",ndecimals,",1x,i6,1x,i7,1x,a1,1x,a1,1x,a10,1x,i3,1x,i2,2i8)"
                write(enunit,my_fmt) & 
                          iroot,energyI-intensity%ZPE,nint(intensity%gns(isymI)*( 2.0_rk*jI + 1.0_rk )),nint(jI),&
                          pm,ef,statename,ivI,(ilambdaI),nint((sigmaI)),nint((omegaI))
                !
              else
                !
                write(my_fmt,'(A,i0,a)') "(i12,1x,f12.",ndecimals,",1x,i6,1x,f7.1,1x,a1,1x,a1,1x,a10,1x,i3,1x,i2,2f8.1)"
                write(enunit,my_fmt) & 
                          iroot,energyI-intensity%ZPE,nint(intensity%gns(isymI)*( 2.0_rk*jI + 1.0_rk )),jI,&
                          pm,ef,statename,ivI,(ilambdaI),(sigmaI),(omegaI)
                          !
              endif
              !
            endif
            !
          endif
          !
          call energy_filter_upper(jI, energyI, passed)
          !
          call energy_filter_lower(jI, energyI, passed_)
          !
          if (.not. passed .and. .not. passed_) cycle
          !
          istateI = eigen(indI, igammaI)%quanta(ilevelI)%istate 
          parity_gu = poten(istateI)%parity%gu 
          isymI = correlate_to_Cs(igammaI, parity_gu)
          !
          nlevelsG(isymI) = nlevelsG(isymI) + 1
          !
        enddo
      enddo
    enddo
    !
    deallocate(vecI)
    call ArrayStop('intensity-vecI')
    !
    if (trim(intensity%linelist_file) /= "NONE") close(enunit, status='keep')
    !
    write(my_fmt,'(A,I0,A)') &
      "('Number of states for each symm = ',", sym%Nrepresen, "i8)"
    !
    write(out,my_fmt) nlevelsG(:)
    !
    matsize = int(sum(nlevelsG(:)), hik)
    !
    if (iverbose >= 4) write(out,"(/'Dipole moment integration (i)...')")
    !
    if (Ntransit == 0) then
      write(out,"('pol_intensity: the transition filters are too tight: no entry')")
      !
      stop 'pol_intensity: the filters are too tight' 
    endif
    !
    write(out, "(/'...done!')")
    !
    allocate(vecI(dimenmax), vecF(dimenmax), stat=info)
    !
    call ArrayStart('intensity-vectors', info, size(vecI), kind(vecI))
    call ArrayStart('intensity-vectors', info, size(vecF), kind(vecF))
    !
    ! loop over final states -> count states for each symmetry
    !
    write(my_fmt,'(A,I0,A)') &
      "('Number of states for each symm = ',", sym%Nrepresen, "i8)"
    !
    write(out,my_fmt) nlevelsG(:)
    !
    if (iverbose >= 0) then
      write(out, "(' Total number of lower states = ',i8)") nlower
      write(out, "(' Total number of transitions  = ',i8)") Ntransit
    end if
    !
    if (iverbose >= 0) then
      write(out, "(/' Statistical weights gns = ',4f12.1)") intensity%gns(1:)
    end if
    !
    ! To speed up the line strength evaluation, we perform the calculation:
    !   S_{if} = | <i|a|f> |^2  = | \sum_{nm} C_in C_fm <n|a|m> |^2
    ! in three steps:
    !   1. Evaluate the expansion of the initial state:
    !         s_{im} = sum_{n} C_in <n|a|m>
    !   2. Evaluate the product with the expansion of the final state:
    !         s_{if} = sum_{m} C_fm s_{im}
    !   3. Square the result to obtain S_{if}
    !         S_{if} = s_{if}^2
    !
    ! The temporary object s_{im} is referred to as the 'half linestrength',
    !   with corresponding variable "half_linestr".
    !
    allocate(half_linestr(dimenmax), stat = info)
    allocate(half_linestr_zero(dimenmax), stat = info)
    !
    call ArrayStart('half_linestr', info, size(half_linestr), &
                    kind(half_linestr))
    !
    if (iverbose >= 5) call MemoryReport
    !
    write(out,"(/a,a,a,a)") &
          'Linestrength S(f<-i) [Debye**2],', ' Transition moments [Debye],', &
          ' Einstein coefficient A(if) [1/s],', ' and Intensities [cm/mol]'
    !
    ! Prepare the table header
    !
    select case (trim(intensity%action))
      !
      case('ABSORPTION')
        !
        write(out, &
              "(/t5,'J',t7,'Gamma <-',t18,'J',t21,'Gamma',t27,'Typ',t37,'Ei',&
              &t44,'<-',t52,'Ef',t64,'nu_if',8x,'S(f<-i)',10x,'A(if)',12x,&
              &'I(f<-i)',7x,'State v lambda sigma  omega <- State v lambda&
              & sigma  omega ')")
        !
        dir = '<-'
        !
      case('EMISSION')
        !
        write(out, &
              "(/t5,'J',t7,'Gamma ->',t18,'J',t21,'Gamma',t27,'Typ',t37,'Ei',&
              &t44,'->',t52,'Ef',t64,'nu_if',8x,'S(i->f)',10x,'A(if)',12x,&
              &'I(i->f)',7x,'State v lambda sigma  omega -> State v lambda&
              & sigma  omega ')")
        !
        dir = '->'
        !
      case('TM')
        !
        write(out,"(/t4,'J',t6,'Gamma <-',t17,'J',t19,'Gamma',t25,'Typ',t35,'Ei',&
              &t42,'<-',t52,'Ef',t65,'nu_if',10x,'TM(f->i)')")
        !
      !
    end select
    !
    deallocate(vecF)
    !
    ! ---------------------------------
    ! The actual intensity calculations
    ! ---------------------------------
    !
    itransit = 0
    !
    ! loop over initial states
    do indI = 1, nJ
      !
      jI = jval(indI)
      !
      ! loop over initial state symmetries
      do igammaI = 1, Nrepresen
        !
        ! number of sublevels and dimension of basis for initial state
        nlevelsI = eigen(indI, igammaI)%Nlevels
        dimenI = eigen(indI, igammaI)%Ndimen
        !
        ! loop over final states
        do indF = 1, nJ
          !
          jF = jval(indF)
          !
          ! selection rules for irreducible second & zeroth rank 
          !   spherical polarisability tensor components
          if (abs(nint(jI - jF)) > 2 ) cycle
          !
          ! loop over final state symmetries
          do igammaF = 1, Nrepresen
            !
            ! number of sublevels and dimension of basis for final state
            nlevelsF = eigen(indF, igammaF)%Nlevels
            dimenF = eigen(indF, igammaF)%Ndimen
            !
            Ilevels_loop : do ilevelI = 1, nlevelsI
              !
              ! energy, quanta and degeneracy of the initial state
              energyI = eigen(indI, igammaI)%val(ilevelI)
              !
              quantaI   => eigen(indI, igammaI)%quanta(ilevelI)
              istateI   = quantaI%istate
              ivibI     = quantaI%ivib
              ivI       = quantaI%v
              sigmaI    = quantaI%sigma
              spinI     = quantaI%spin
              ilambdaI  = quantaI%ilambda
              omegaI    = quantaI%omega
              !
              ! reconstruct symmetry for C2v which is different to Cs case
              parity_gu = poten(istateI)%parity%gu
              isymI = correlate_to_Cs(igammaI, parity_gu)
              !
              call energy_filter_lower(jI, energyI, passed)
              !
              if (.not. passed) cycle
              !
              vecI(1:dimenI) = eigen(indI, igammaI)%vect(1:dimenI, ilevelI)
              !
              ! Compute the half-linestrength
              !
              ! if no transitions from ilevelI -> jF exist then skip 
              !   calculation, so we check for this condition
              passed = .false.
              !
              ! loop over final states
              do ilevelF = 1, nlevelsF 
                !
                ! energy and quanta of final state
                energyF = eigen(indF,igammaF)%val(ilevelF)
                !
                quantaF   => eigen(indF,igammaF)%quanta(ilevelF)
                istateF   = quantaF%istate
                ivibF     = quantaF%ivib
                ivF       = quantaF%v
                sigmaF    = quantaF%sigma
                spinF     = quantaF%spin
                ilambdaF  = quantaF%ilambda
                omegaF    = quantaF%omega
                !
                ! reconstruct symmetry for C2v which is different to Cs case
                parity_gu = poten(istateF)%parity%gu
                isymF = correlate_to_Cs(igammaF, parity_gu)
                !
                call intens_filter(jI, jF, energyI, energyF, isymI, isymF, &
                                   igamma_pair, passed)
                !
                ! prevents diagonal matrix elements from being skipped
                if (intensity%matelem) &
                  call matelem_filter(jI, jF, energyI, energyF, isymI, isymF,&
                                      igamma_pair,passed)
                !
                ! stop checking when a transition passed the filter
                if (passed) exit
                !
              enddo
              !
              ! if no transition pass the filter, cycle to next state
              if (.not. passed) cycle
              !
              select case (trim(intensity%action))
                !
                case('ABSORPTION', 'EMISSION')
                  !
                  if (isymF /= igamma_pair(isymI)) cycle
                  !
                  if ((intensity%J(1) + intensity%J(2) > 0) &
                      .and. abs(nint(jI - jF)) <= 2 ) then
                    !
                    call do_1st_half_linestrength(jI, jF, indI, indF, dimenI, &
                                                  dimenF, vecI(1:dimenI), &
                                                  half_linestr, &
                                                  half_linestr_zero)
                    !
                  endif
                  !
                case('TM')
                  !
                  call do_1st_half_tm(indI, indF, dimenI, dimenF, &
                                      vecI(1:dimenI), half_linestr)
                  !
              end select
              !
              ! loop over final states
              allocate(vecF(dimenmax), stat = alloc_p)
              !
              if (alloc_p /= 0) then
                write(out, &
                      "(' dipole: ',i9,' trying to allocate array -vecF')") &
                  alloc_p
                !
                stop 'dipole-vecF - out of memory'
                !
              end if
              !
              Flevels_loop : do ilevelF = 1, nlevelsF 
                !
                ! energy and quanta of final state
                energyF = eigen(indF, igammaF)%val(ilevelF)
                !
                ! dimension of bases for final state
                dimenF = eigen(indF, igammaF)%Ndimen
                !
                quantaF => eigen(indF, igammaF)%quanta(ilevelF)
                !
                istateF  = quantaF%istate
                ivibF    = quantaF%ivib
                ivF      = quantaF%v
                sigmaF   = quantaF%sigma
                spinF    = quantaF%spin
                ilambdaF = quantaF%ilambda
                omegaF   = quantaF%omega
                !
                call energy_filter_upper(jF, energyF, passed)
                !
                if (.not.passed) cycle Flevels_loop
                !
                parity_gu = poten(istateF)%parity%gu
                isymF = correlate_to_Cs(igammaF, parity_gu)
                !
                call intens_filter(jI, jF, energyI, energyF, isymI, isymF,&
                                   igamma_pair, passed)
                !
                if (intensity%matelem) &
                  call matelem_filter(jI, jF, energyI, energyF, isymI, &
                                      isymF, igamma_pair, passed)
                !
                if (.not. passed) cycle Flevels_loop
                !
                ! Find PQR branch for transition
                branch = PQR_branch(jI, jF)
                !
                nu_if = energyF - energyI
                !
                ! filer out zero-frequency transitions
                if (nu_if < small_) cycle 
                !
                ! Count processed transitions
                itransit = itransit + 1
                !
                vecF(1:dimenF) = eigen(indF, igammaF)%vect(1:dimenF, ilevelF)
                !
                select case (trim(intensity%action))
                  !
                  case default
                    !
                    stop 'only ABSORPTION and TM are properly coded'
                    !
                  case('ABSORPTION', 'EMISSION')
                    !
                    linestr =      ddot(dimenF, half_linestr, 1, vecF, 1)
                    linestr_zero = ddot(dimenF, half_linestr_zero, 1, vecF, 1)
                    !
                    linestr2 = linestr**2
                    linestr2_zero = linestr_zero**2
                    !
                    !
                    ! calculate intensity
                    A_einst = A_coef_s_1 * (2.0_rk * jI + 1.0_rk) &
                              * linestr2 * abs(nu_if)**3
                    !
                    linestr2 = linestr2 * intensity%gns(isymI) &
                              * (2.0_rk * jI + 1.0_rk) &
                              * (2.0_rk * jF + 1.0_rk)
                    !
                    linestr2_zero = linestr2_zero * intensity%gns(isymI) &
                              * (2.0_rk * jI + 1.0_rk) &
                              * (2.0_rk * jF + 1.0_rk)
                    !
                    if (trim(intensity%action) == 'ABSORPTION') then
                      !
                      ! Boltzmann factor
                      boltz_fc = exp( - (energyI - intensity%ZPE) * beta) &
                                * (1.0_rk - exp( - abs(nu_if) * beta)) &
                                / (intensity%part_func * nu_if**2)
                      !
                      ! absorption intensity in cm/mol
                      absorption_int = 1.0_rk / (8.0_rk * pi * vellgt) &
                                      * intensity%gns(isymF) &
                                      * (2.0_rk * jF + 1.0_rk) &
                                      * A_einst * boltz_fc
                    else
                      !
                      ! emissivity in Ergs/mol/Sr
                      boltz_fc = exp( - (energyF - intensity%ZPE) * beta)
                      !
                      ! emission intensity in cm/mol (despite label)
                      absorption_int = emcoef * A_einst * boltz_fc &
                                      * intensity%gns(isymI) &
                                      * (2.0_rk * jF + 1.0_rk) &
                                      * nu_if / intensity%part_func
                      !
                    endif
                    !
                    if (absorption_int >= intensity%threshold%intensity .and.&
                        linestr2 + linestr2_zero >= &
                        intensity%threshold%linestrength ) then
                      !
                      write(out, &
                            "( (f5.1, 1x, a4, 3x),a2, &
                            &(f5.1, 1x, a4, 3x),a1,&
                            &(2x, f11.4,1x),a2,&
                            &(1x, f11.4,1x),f11.4,2x, 3(1x, es16.8),&
                            &' ( ',i2,1x,i3,1x,i2,2f8.1,' )',a2,&
                            &'( ',i2,1x,i3,1x,i2,2f8.1,' )')") &
                        jF, sym%label(isymF), dir, &
                        jI, sym%label(isymI), branch, &
                        energyF - intensity%ZPE, dir, &
                        energyI - intensity%ZPE, nu_if, &
                        linestr2, A_einst, absorption_int, &
                        istateF, ivF, ilambdaF, sigmaF, omegaF, dir, &
                        istateI, ivI, ilambdaI, & sigmaI, omegaI
                        !
                      !
                      ! generate the line list (transition file)
                      if (trim(intensity%linelist_file) /= "NONE") then
                        !
                        if (intensity%matelem) then
                          !
                          write(richunit(indI, indF), &
                                "(i8,1x,i8,1x,i4,1x,i4,3x,e23.16,3x,e23.16)") & 
                                quantaI%iJ_ID, quantaF%iJ_ID, 1, 1,linestr_zero,linestr
                          !
                        else
                          !
                          write(transunit,"(i12,1x,i12,2x,es10.4,4x,f16.6)") &
                            quantaF%iroot, quantaI%iroot, A_einst, nu_if
                        endif
                      endif
                    endif
                    !
                  case('TM')
                    !
                    tm = dot_product(half_linestr(1:dimenF), vecF(1:dimenF))
                    !
                    linestr = tm 
                    !
                    if (linestr >= intensity%threshold%intensity) then
                      !
                      write(out, &
                          "( (i4, 1x, a3, 3x),'->', (i4, 1x, a3, 3x),a1, &
                          &(2x, f13.6,1x),'->',(1x, f13.6,1x),f12.6, f15.8)") &
                        jI, sym%label(isymI), jF, sym%label(isymF), branch, &
                        linestr, itransit, tm
                      !
                    endif
                  !
                end select
                !
              enddo Flevels_loop
              !
              deallocate(vecF)
              !
              if (iverbose >= 5) call TimerReport
              !
            enddo Ilevels_loop
          enddo
        enddo
      enddo
      !
      ! close some J-files for Richmol:
      if (intensity%matelem) then
        !
        do indF = max(1, indI - 2), min(nJ, indI + 2)
          !
          jF = Jval(indF)
          if (nint(abs(jI - jF)) > 2) cycle
          !
          if (jI > jF) cycle
          !
          write(richunit(indI,indF),"('End richmol format')")
          close(richunit(indI,indF))
          !
        enddo
      endif
    enddo
    !
    deallocate(vecI)
    call ArrayStop('intensity-vectors')
    !
    deallocate(half_linestr)
    call ArrayStop('half_linestr')
    !
    if (trim(intensity%linelist_file) /= "NONE") &
      close(transunit, status="keep")
    !
    call TimerStop('Intensity calculations')
    !
  end subroutine pol_intensity
  !
  !
  !
  function PQR_branch(jI, jF) result (X)
    !
    real(rk),intent(in)  :: jI,jF
    character(len=1)        :: X
    !
    select case(trim(intensity%action))
        !
      case default 
        !
        X = 'D'
        !
      case ('ABSORPTION','EMISSION')
        !
        if ( jI>jF ) then
          X = 'P'
        elseif(nint(jI - jF) /=0) then
          X = 'R'
        else
          X = 'Q'
        endif
        !
      !case ('EMISSION')
      !  !
      !  if (jI>jF) then
      !    X = 'R'
      !  elseif(jI/=jF) then
      !    X = 'P'
      !  else
      !    X = 'Q'
      !  endif
        !
    end select
    !
  end function PQR_branch
  !
  !
  !
  subroutine find_igamma_pair(igamma_pair)
    !
    integer(ik),intent(out) :: igamma_pair(:)
    integer(ik)     :: igammaI,igammaF,ngamma
    !
    if (trim(intensity%action)=='TM') then 
      !
      igamma_pair = 1
      return 
      ! 
    endif
    !
    do igammaI = 1, sym%Nrepresen
      !
      ! count number of hits
      !
      ngamma = 0
      igamma_pair(igammaI) = igammaI
      !
      do igammaF = 1,sym%Nrepresen
        !
        if (igammaI == igammaF .and. &
          intensity%isym_pairs(igammaI) == intensity%isym_pairs(igammaF)) then 
          !
          igamma_pair(igammaI) = igammaF
          !
          ngamma = ngamma + 1 
          !
          if (ngamma > 1) then 
            !
            write(out,"('pol_intensity: Assumption that selection rules come in pairs is wrong!')")
            stop 'pol_intensity: Assumption that all selection rules work in pairs is wrong!'
            !
          endif   
          !
        endif
        !
      enddo
      !
      if ( intensity%gns(igammaI)/=intensity%gns(igamma_pair(igammaI)) ) then 
        !
        write(out,"('pol_intensity: selection rules do not agree with Gns')")
        stop 'pol_intensity: selection rules do not agree with Gns!'
        !
      endif   
      !
    enddo 
    !
  end subroutine find_igamma_pair
  !
  !
  !
  subroutine energy_filter_lower(J,energy,passed)
    !
    real(rk),intent(in) :: J
    real(rk),intent(in)    :: energy
    !type(quantaT),intent(in) :: quanta
    logical,intent(out)    :: passed
    !
    ! passed = .true.
    !
    ! if (.not.intensity%do) return
    !
    passed = .false.
    !
    if (                                                             &
        ! nuclear stat.weight: 
        !
        J>=intensity%J(1).and.                                       &
        J<=intensity%J(2).and.                                       &
        !
        energy-intensity%ZPE>=intensity%erange_low(1).and.           &
        energy-intensity%ZPE<=intensity%erange_low(2)    ) then 
        !
        !quanta%istate >= intensity%lower%istate.and.   &
        !quanta%v >= intensity%lower%v.and.             &
        !quanta%omega >= intensity%lower%omega  
        !
        passed = .true.
        !
    endif
    !
  end subroutine energy_filter_lower
  !
  !
  !
  subroutine energy_filter_upper(J,energy,passed)
    !
    real(rk),intent(in) :: J
    real(rk),intent(in)    :: energy
    !type(quantaT),intent(in) :: quanta
    logical,intent(out)    :: passed
    !
    ! passed = .true.
    !
    ! if (.not.intensity%do) return
    !
    passed = .false.
    !
    if (                                                             &
        ! nuclear stat. weight: 
        !
        J>=intensity%J(1).and.                                       &
        J<=intensity%J(2).and.                                       &
        !
        energy-intensity%ZPE>=intensity%erange_upp(1).and.           &
        energy-intensity%ZPE<=intensity%erange_upp(2)   ) then 
        !
        !quanta%istate >= intensity%upper%istate.and.   &
        !quanta%v >= intensity%upper%v.and.             &
        !quanta%omega >= intensity%upper%omega
        !
        passed = .true.
        !
    endif 
    !
  end subroutine energy_filter_upper
  !
  !
  !
  subroutine intens_filter(jI, jF, energyI, energyF, isymI, isymF, &
                          igamma_pair, passed)
    !
    real(rk),intent(in) :: jI,jF
    integer(ik),intent(in) :: isymI,isymF
    real(rk),intent(in)    :: energyI,energyF
    integer(ik),intent(in) :: igamma_pair(sym%Nrepresen)
    real(rk)               :: nu_if
    logical,intent(out)    :: passed

    passed = .false.
    !
    nu_if = energyF - energyI
    !
    !if (trim(intensity%action)=='EMISSION') nu_if = -nu_if 
    !
    if (                                                             &
        ! nuclear stat.weight: 
        !
        intensity%gns(isymI)>small_.and.                             &
        !
        ! absorption/emission go only in one direction
        !
        nu_if>small_.and.                                            &
        !
        ! spectroscopic window
        !
        nu_if>=intensity%freq_window(1).and.                         &
        nu_if<=intensity%freq_window(2).and.                         &
        !
        jI>=intensity%J(1).and.                                      &
        jI<=intensity%J(2).and.                                      &
        !
        jF>=intensity%J(1).and.                                      &
        jF<=intensity%J(2).and.                                      &
        !
        energyI-intensity%ZPE>=intensity%erange_low(1).and.          &
        energyI-intensity%ZPE<=intensity%erange_low(2).and.          &
        !
        energyF-intensity%ZPE>=intensity%erange_upp(1).and.          &
        energyF-intensity%ZPE<=intensity%erange_upp(2)   ) then 
        !

        !quantaI%istate >= intensity%lower%istate.and.   &
        !quantaI%v >= intensity%lower%v.and.             &
        !quantaI%omega >= intensity%lower%omega.and.     &
        !
        !quantaF%istate >= intensity%upper%istate.and.   &
        !quantaF%v >= intensity%upper%v.and.             &
        !quantaF%omega >= intensity%upper%omega  ) then 
        !
        passed = .true.
        !
    endif 
    !
    if (trim(intensity%action)=='ABSORPTION' .or. &
        trim(intensity%action)=='EMISSION') then 
      !
      ! In order to avoid double counting of transitions we exclude jI = jF 
      !  == intensity%J(2), i.e. Q branch for highest J is never considered
      !
      passed = passed .and.                                              &
      ! ????? 
      (jF/=intensity%J(1).or.jI/=intensity%J(1).or.nint(jI+jF)==1).and.  &
      !
      !((nint(jF - intensity%J(1)) /=0 .or.                          &
      !  nint(jI - intensity%J(1)) /= 0) .and.                       &
      !  intensity%J(1)>0 ) .and.                                    &
      ( intensity%J(1)+intensity%J(2)>0 ).and.                       &
      !
      ! selection rules: 
      !
      intensity%isym_pairs(isymI)==intensity%isym_pairs(isymF).and.  &
      !
      igamma_pair(isymI)==isymF.and.                                 &
      !
      ! selection rules from the 3j-symbols
      !
      abs(nint(jI-jF))<=2.and.nint(jI+jF)>=2
      !
    endif
    !
  end subroutine intens_filter
  !
  !
  !
  subroutine matelem_filter(jI, jF, energyI, energyF, isymI, isymF, &
                            igamma_pair, passed)
    !
    real(rk),intent(in) :: jI,jF
    integer(ik),intent(in) :: isymI,isymF
    real(rk),intent(in)    :: energyI,energyF
    integer(ik),intent(in) :: igamma_pair(sym%Nrepresen)
    real(rk)               :: nu_if
    logical,intent(out)    :: passed

    passed = .false.
    !
    nu_if = energyF - energyI
    !
    !if (trim(intensity%action)=='EMISSION') nu_if = -nu_if 
    !
    if (                                                             &
        ! nuclear stat.weight: 
        !
        intensity%gns(isymI)>small_.and.                           &
        !
        ! absorption/emission go only in one direction
        !
        (Jf-Ji)>-small_.and.                                        &
        !
        ! spectroscopic window
        !
        nu_if>=intensity%freq_window(1).and.                         &
        nu_if<=intensity%freq_window(2).and.                         &
        !
        jI>=intensity%J(1).and.                                      &
        jI<=intensity%J(2).and.                                      &
        !
        jF>=intensity%J(1).and.                                      &
        jF<=intensity%J(2).and.                                      &
        !
        energyI-intensity%ZPE>=intensity%erange_low(1).and.          &
        energyI-intensity%ZPE<=intensity%erange_low(2).and.          &
        !
        energyF-intensity%ZPE>=intensity%erange_upp(1).and.          &
        energyF-intensity%ZPE<=intensity%erange_upp(2)   ) then 
        !
        passed = .true.
        !
    endif 
    !
    if (trim(intensity%action)=='ABSORPTION'.or.&
        trim(intensity%action)=='EMISSION') then 
      !
      ! In order to avoid double counting of transitions we exclude jI = jF 
      ! == intensity%J(2), i.e. Q branch for highest J is never considered:
      !
      passed = passed.and.                                              &
      !
      (jF/=intensity%J(1).or.jI/=intensity%J(1).or.nint(jI+jF)==1).and. &
      !
      !( ( nint(jF-intensity%J(1))/=0 .or.                            &
      ! nint(jI-intensity%J(1))/=0 ).and.                             &
      ! intensity%J(1)>0 ) .and.                                      &
      ( intensity%J(1)+intensity%J(2)>0 ).and. &
      !
      ! selection rules: 
      !
      intensity%isym_pairs(isymI)==intensity%isym_pairs(isymF).and.   &
      !
      igamma_pair(isymI)==isymF.and.                                  &
      !
      ! selection rules from the 3j-symbols
      !
      abs(nint(jI-jF))<=2.and.nint(jI+jF)>=2
      !
    endif
    !
  end subroutine matelem_filter
  !
  !
  !
  subroutine Jgamma_filter(jI, jF, isymI, isymF, igamma_pair, passed)
    !
    real(rk),intent(in) :: jI,jF
    integer(ik),intent(in) :: isymI,isymF
    integer(ik),intent(in) :: igamma_pair(sym%Nrepresen)
    logical,intent(out)    :: passed
    !
    passed = .true.
    !
    if (trim(intensity%action)=='ABSORPTION'.or. &
        trim(intensity%action)=='EMISSION') then 
      !
      ! In order to avoid double counting of transitions we exclude jI = jF 
      ! == intensity%J(2), i.e. Q branch for highest J is never considered:
      !
      passed = passed.and.                                            &
      !
      !((nint(jF-intensity%J(1))/=0 .or.                              &
      !  nint(jI-intensity%J(1))/=0 ).and.                            &
      !  intensity%J(1)+intensity%J(2)>0 ).and.                       &
      ( intensity%J(1)+intensity%J(2)>0 ).and.                        &
      !
      ! selection rules: 
      !
      intensity%isym_pairs(isymI)==intensity%isym_pairs(isymF).and.   &
      !
      igamma_pair(isymI)==isymF.and.                                  &
      !
      ! selection rules from the 3j-symbols
      !
      abs(nint(jI-jF))<=1.and.nint(jI+jF)>=1
      !
    endif
    !
  end subroutine Jgamma_filter
  !
  !
  !
  subroutine Sort_levels(iverbose,njval, jval)
    !
    integer(ik), intent(in) :: iverbose,njval
    real(rk), intent(in)    :: jval(njval)

    integer(ik)             :: jind, nroots, nlevels,  iroot, ilevel, jlevel, &
                               info,itau,jtau

    real(rk)                :: energy

    !
    logical                 :: passed
    integer(ik)             :: iind,jroot,igamma
    !
    if (iverbose>=2) &
      write(out,"(/'Read and sort eigenvalues in increasing order...')")
    !
    call TimerStart('Sort eigenvalues')

    if (.not. allocated(basis)) &
      stop 'Sort_levels error: associated(basis) = .false.'
    !
    ! In practice we do not need all stored roots, only those that pass the 
    ! filters given in the input. Here we perform a pre-selection of the 
    ! states which are actually required in the calculations. Towards this we ! read the unit twice:
    ! a) to count and 
    ! b) to read.  
    ! After reading all states will  be combined together and sorted wrt 
    ! energy increasing to produce just one list of states. 
    !
    ! Estimate the maximal number of energy records
    !
    Nlevels = 0
    do jind = 1, njval
      do itau = 1,sym%NrepresCs
         Nlevels = max(eigen(jind,itau)%Nlevels,Nlevels)
      enddo
    enddo 
    !
    nroots  = 0
    !
    do jind = 1, njval
      !
      !irot = jind
      do itau = 1,sym%NrepresCs
        !
        Nlevels = eigen(jind,itau)%Nlevels
        !
        do ilevel = 1,Nlevels
          !
          energy = eigen(jind,itau)%val(ilevel)
          !
          !state => eigen(iJ)%quanta
          !
          igamma =  eigen(jind,itau)%quanta(ilevel)%igamma
          !
          if (job%ZPE<0.and.jind==1.and.ilevel==1) job%zpe = energy
          !
          call Energy_filter(Jval(jind),energy,igamma,passed)
          !
          if (passed) then 
            !
            nroots = nroots + 1
            !
          endif 
          !
        enddo
      enddo
    enddo
    !
    if (nroots==0) then 
      write(out,"('Sort_levels: the filters are too tight: no entry')") 
      !stop 'Sort_levels: the filters are too tight' 
    endif 
    !
    ! total number of levels to be processed, a global variable 
    !
    Neigenlevels = nroots
    !
    return
    !
    allocate(Elevel(Neigenlevels),stat = info)
    if (info /= 0) stop 'Sort_levels allocation error: eigen - out of memory'
    !
    if (iverbose>=4) then 
      write(out,"('   Number of selected eigenvectors: ',i8)") Neigenlevels
    end if
    !
    ! Now we can actually read and store the energies and their description. 
    !
    iroot  = 0
    !
    do jind = 1, njval
      do itau = 1,sym%NrepresCs
        !
        Nlevels = eigen(jind,itau)%Nlevels
        !
        do ilevel = 1,Nlevels
          !
          iroot  = iroot + 1
          !
          Elevel(iroot)%jind   = jind
          Elevel(iroot)%igamma = itau
          Elevel(iroot)%ilevel = ilevel
          !
        enddo
      enddo
    enddo
    !
    if (iroot/=Neigenlevels) then 
      write(out, '("The number of selected levels ",i0,&
            &" does not agree with Neigenlevels =  ",i0)') &
            ilevel, Neigenlevels
      stop 'iroot/=Neigenlevels'
    endif
    !
    ! Sort according the energy increasing 
    !
    do iroot =1,Neigenlevels
      !
      ilevel = Elevel(iroot)%ilevel
      iind = Elevel(iroot)%jind
      itau = Elevel(iroot)%igamma
      !
      energy = eigen(iind,itau)%val(ilevel)
      !
      do jroot =iroot+1,Neigenlevels
        !
        jind   = Elevel(jroot)%jind
        jtau   = Elevel(jroot)%igamma
        jlevel = Elevel(jroot)%ilevel
        !
        if (eigen(jind,jtau)%val(jlevel)<energy) then 
          !
          Elevel(jroot)%jind= iind ; Elevel(jroot)%ilevel = ilevel ; Elevel(jroot)%igamma = itau
          Elevel(iroot)%jind= jind ; Elevel(iroot)%ilevel = jlevel ; Elevel(jroot)%igamma = jtau
          ilevel = Elevel(iroot)%ilevel
          iind   = Elevel(iroot)%jind
          itau   = Elevel(iroot)%igamma
          energy = eigen(iind,itau)%val(ilevel)
          !
        endif 
        !
      enddo
      !
    enddo
    !
    if (iverbose>=2) write(out,"('...done!')")
    !
    call TimerStop('Sort eigenvalues')
    !
  end subroutine Sort_levels
  !
  !
  !
  !`-------------------------------------------------------------------------
  !                 Matrix Element/Line Strength Evaluation                  
  ! -------------------------------------------------------------------------
    !   To speed up the line strength evaluation, we perform the calculation: 
    !       S_{if} = | <i|a|f> |^2  = | \sum_{nm} C_in C_fm <n|a|m> |^2       
    !   in three steps:                                                       
    !     1. Evaluate the expansion of the initial state:                     
    !             s_{im} = sum_{n} C_in <n|a|m>                               
    !     2. Evaluate the product with the expansion of the final state:      
    !             s_{if} = sum_{m} C_fm s_{im}                                
    !     3. Square the result to obtain S_{if}                               
    !             S_{if} = s_{if}^2                                           
    ! -----------------------------------------------------------------------
  !
  !
  subroutine do_1st_half_linestrength(jI, jF, indI, indF, &
                                      dimenI, dimenF, vector, &
                                      half_ls, half_ls_zero)
    !
    implicit none
    ! I/O variables
    real(rk), intent(in)    :: jI, jF
    integer(ik), intent(in) :: indI, indF, dimenI, dimenF
    real(rk), intent(in)    :: vector(:)
    real(rk), intent(out)   :: half_ls(:), half_ls_zero(:)
    !
    ! states
    integer(ik)             :: ivibF, ivibI, istateF, istateI, &
                               ilambdaF, ilambdaI
    real(rk)                :: omegaF, omegaI, sigmaF, sigmaI, spinF, spinI
    ! faux states
    integer(ik)             :: omegaF_, omegaI_, istateF_, istateI_, &
                               ilambdaF_, ilambdaI_
    real(rk)                :: spinF_, spinI_
    !
    ! counters
    integer(ik)             :: icontrF, icontrI, ipermute, isigmav, itau, idip
    !
    ! other
    real(rk)                :: f3j, ls, f_t, rank
    type(fieldT), pointer   :: field
    !
    call TimerStart('do_1st_half_linestr')		
    !
    half_ls = 0
    half_ls_zero = 0
    !
    ! loop over final states
    loop_F : do icontrF = 1, dimenF
      !
      ivibF = basis(indF)%icontr(icontrF)%ivib
      istateF = basis(indF)%icontr(icontrF)%istate
      omegaF = basis(indF)%icontr(icontrF)%omega
      sigmaF = basis(indF)%icontr(icontrF)%sigma
      spinF = basis(indF)%icontr(icontrF)%spin
      ilambdaF = basis(indF)%icontr(icontrF)%ilambda
      !
      ! ----- Note -----
      ! 1/2 integer spin cases introduce an imaginary factor due to the 
      ! (-1)^omega term. To simplify the code, we treat these cases as real by
      ! setting omega equal to the next lowest integer, such that the 
      ! multiplicative pre-factor in the matrix element equation is real.
      ! Note that 1/2 integer values of omega are still used in the
      ! calculation of the  3-j symbols.
      !
      ! check for 1/2 integer spin cases in final state
      omegaF_ = nint(omegaF)
      if (mod(nint(2.0_rk*omegaF + 1.0_rk), 2) == 0) &
        omegaF_ = nint(omegaF - 0.5_rk)
      !
      ! loop over initial states
      loop_I : do icontrI = 1, dimenI
        !
        ivibI = basis(indI)%icontr(icontrI)%ivib
        istateI = basis(indI)%icontr(icontrI)%istate
        omegaI = basis(indI)%icontr(icontrI)%omega
        sigmaI = basis(indI)%icontr(icontrI)%sigma
        spinI = basis(indI)%icontr(icontrI)%ilambda
        ilambdaI = basis(indI)%icontr(icontrF)%ilambda
        !
        ! selection rules
        if (abs(nint(omegaF - omegaI)) /= 0 &
            .or. nint(spinI - spinF)   /= 0 &
            .or. nint(sigmaI - sigmaF) /= 0)&
          cycle loop_I
        !
        if (abs(ilambdaI - ilambdaF) /= abs(nint(omegaF - omegaI))) &
          cycle loop_I
        !
        ! check for half integer spin in inital state
        omegaI_ = nint(omegaI)
        if (mod(nint(2.0_rk * omegaI + 1.0_rk), 2) == 0) &
          omegaI_ = nint(omegaI - 0.5_rk)
        !
        ls = 0
        !
        ! idip = Ndipoles + 1 is the zeroth rank component
        loop_idipole : do idip = 1, Ndipoles + Nss
          !
          if (idip <= Ndipoles) then
            field => dipoletm(idip)
            rank = 2.0_rk
          else
            ! use spinspin as field pointer for zeroth rank
            field => spinspin(idip-Ndipoles)
            rank = 0.0_rk
          endif
          !
          ! calculate 3-j symbol and check selection rule
          f3j = three_j(jI, jF, rank, -omegaI, omegaF, omegaI - omegaF)
          if (abs(f3j) < intensity%threshold%coeff) cycle loop_idipole
          !
          do ipermute = 0, 1
            !
            if (ipermute == 0) then
              !
              istateI_ = field%istate
              ilambdaI_ = field%lambda
              spinI_ = field%spini
              !
              istateF_ = field%jstate 
              ilambdaF_ = field%lambdaj
              spinF_ = field%spinj
              !
            else
              !
              istateF_ = field%istate
              ilambdaF_ = field%lambda
              spinF_ = field%spini
              !
              istateI_ = field%jstate
              ilambdaI_ = field%lambdaj
              spinI_ = field%spinj
              !
            endif
            !
            ! no permutation on diagonal entries, else double counting
            if (ipermute == 1 .and. &
                istateI_ == istateF_ .and. &
                ilambdaI_ == ilambdaF_ .and. &
                nint(spinI_ - spinF_) == 0) cycle
            !
            ! check we are at the right electronic state
            if (istateI /= istateI_ .or. istateF /= istateF_) cycle
            !
            ! also need to account for a change in sign of Lambda as the input 
            ! gives only lambda > 0. This is equivalent to a laboratory fixed 
            ! inversion, i.e sigmav operation (sigmav = unitary transformation)
            do isigmav = 0, 1
              !
              ! skip inversion if no quanta in state
              if (isigmav == 1 .and. &
                  abs(field%lambda) + abs(field%lambdaj) == 0) cycle
              !
              ! do transformation
              ilambdaI_ = ilambdaI_ * (-1)**isigmav
              ilambdaF_  = ilambdaF_ * (-1)**isigmav
              !
              ! proceed only if field quantum numbers are equal to 
              ! corresponding <i| and |j> quantum numbers
              if (ilambdaI_ /= ilambdaI .or. ilambdaF_ /= ilambdaF) cycle
              !
              ! delta Lambda = +/- 2 selection rule
              if ( abs(ilambdaI - ilambdaF) > 2 .or.abs(ilambdaI - ilambdaF) == 1 ) cycle
              !
              f_t = field%matelem(ivibI, ivibF)
              !
              ! apply result of symmetry transformation
              if (isigmav == 1) then
                !
                itau = 0
                !
                if (ilambdaI_ == 0 .and. poten(istateI)%parity%pm == -1) &
                  itau = itau + 1
                !
                if (ilambdaF_ == 0 .and. poten(istateI)%parity%pm == -1) &
                  itau = itau + 1
                !
                f_t = f_t*(-1.0_rk)**(itau)
                !
              endif
              !
              ! product of 3-j symbol parity, mat. elem and vector coefficient
              ls = (-1.0_rk)**(omegaI_) * f_t * f3j * vector(icontrI)
              !
              if (idip <= Ndipoles) then
                ! add to line strength for second rank
                half_ls(icontrF) = half_ls(icontrF) + ls
              else
                ! add to line strength for zeroth rank
                half_ls_zero(icontrF) = half_ls_zero(icontrF) + ls
                !
                ! forced selection rule
                !if (jI /= jF .and. nint(omegaI - omegaF) /= 0) half_ls_zero(icontrF) = 0
                !
              endif
              !
            enddo
          enddo
        enddo loop_idipole
      enddo loop_I
    enddo loop_F
    !
    call TimerStop('do_1st_half_linestr')
    !
  end subroutine do_1st_half_linestrength
  !
  !
  !
  subroutine do_1st_half_tm(indI, indF, dimenI, dimenF, vector, half_tm)
    !
    integer(ik),intent(in)  :: indI,indF,dimenI,dimenF
    real(rk),intent(in)     :: vector(:)
    real(rk),intent(out)    :: half_tm(:)
    integer(ik)             :: icontrF, icontrI, ivibF, ivibI, idip, &
                               istateI,istateF
    real(rk)                :: f
    !
    half_tm = 0
    !
    ! loop over final state basis components
    loop_F : do icontrF = 1, dimenF
      !
      ivibF = basis(indF)%icontr(icontrF)%ivib
      istateF = basis(indF)%icontr(icontrF)%istate
      !isigmaF = basis(indF)%icontr(icontrF)%sigma
      !
      loop_I : do icontrI = 1, dimenI
        !
        ivibI = basis(indI)%icontr(icontrI)%ivib
        istateI = basis(indI)%icontr(icontrI)%istate
        !isigmaI = basis(indI)%icontr(icontrI)%sigma
        !
        !compute TM
        !
        do idip = 1, Ndipoles
          !
          if (dipoletm(idip)%istate /= istateI .or. &
              dipoletm(idip)%jstate/=istateF) cycle
          !
          f = dipoletm(idip)%matelem(ivibI,ivibF)
          !
          half_tm(icontrF) = half_tm(icontrF) + f*vector(icontrI)
          !
        enddo
        !
      enddo loop_I
      !
    enddo loop_F
  !
  end subroutine do_1st_half_tm
  !
  !
  !
  subroutine do_LF_matrix_elements(iLF, iunit, jI, jF, icount)
    implicit none
    real(rk), intent(in)                :: jI, jF
    integer(ik), intent(in)             :: iLF, iunit
    integer(ik), intent(out), optional  :: icount
    !
    integer(ik)                         :: icount_, sphComp
    integer(ik)                         :: Mfx2, Mix2, Mi_, Mf_
    real(rk)                            :: Mf, Mi
    !
    real(rk)                            :: f3j, M_2, M_0
    real(rk)                            :: T(9,0:2,-2:2)
    ! 
    call TimerStart('do_LF_matrix_elements')
    !
    ! define transformation matrix from cartesian to irreducible spherical
    ! tensor components. Transformation matrix elements involving y-component 
    ! are imaginary - we make real here but indicate in RichMol.
    !
    ! transformation matrix elements = T(iLF, k, m), where 1=<iLF=<9 is the
    ! Cartesian laboratory component (defined in pol_intensity subroutine),
    ! k is the irreducible spherical tensor rank and m the component of that 
    ! tensor.
    !
    T = 0 
    !
    ! ONLY THE SECOND RANK SPHERICAL TENSOR IS CURRENTLY PROGRAMMED
    !
    ! Transformation matrix elements:
    ! spherical components of alpha_xx
      T(1,0,0) = -1.0_rk/sqrt(3.0_rk)
      T(1,2,2) =  0.5_rk
      T(1,2,0) = -1.0_rk/sqrt(6.0_rk)
      T(1,2,-2) = 0.5_rk
    !
    ! spherical components of alpha_xy (imaginary)
      T(2,1,0) = -1.0_rk/sqrt(2.0_rk)
      T(2,2,2) = -0.5_rk
      T(2,2,-2) = 0.5_rk
    !
    ! spherical components of alpha_xz
      T(3,1,1)  = -0.5_rk
      T(3,1,-1) =  0.5_rk
      T(3,2,1)  = -0.5_rk
      T(3,2,-1) = 0.5_rk
    !
    ! spherical components of alpha_yy
      T(4,0,0) = -1.0_rk/sqrt(3.0_rk)
      T(4,2,2) = -0.5_rk
      T(4,2,0) = -1.0_rk/sqrt(6.0_rk)
      T(4,2,-2) = -0.5_rk
    !
    ! spherical components of alpha_yz (imaginary)
      T(5,1,1)  =  0.5_rk
      T(5,1,-1) =  0.5_rk
      T(5,2,1)  = -0.5_rk
      T(5,2,-1) = -0.5_rk
    !
    ! spherical components of alpha_zz
      T(6,0,0) = -1.0_rk/sqrt(3.0_rk)
      T(6,2,0) =  2.0_rk/sqrt(6.0_rk)
    !
    ! spherical components of alpha_yx (imaginary)
      T(7,1,0) = 1.0_rk/sqrt(2.0_rk)
      T(7,2,2) = -0.5_rk
      T(7,2,-2) = 0.5_rk
    !
    ! spherical components of alpha_zx
      T(8,1,1) = 0.5_rk
      T(8,1,-1) = 0.5_rk
      T(8,2,1) = -0.5_rk
      T(8,2,-1) = 0.5_rk
    !
    ! spherical components of alpha_zy (imaginary)
      T(9,1,1) = -0.5_rk
      T(9,1,-1) = 0.5_rk
      T(9,2,1) = 0.5_rk
      T(9,2,-1) = 0.5_rk
    !
    !
    M_2 = 0
    M_0 = 0
    !
    icount_ = 0
    !
    ! loop final ang. mom. projection on lab axis
    loop_F : do Mfx2 = -nint(2.0_rk * Jf), nint(2.0_rk * Jf), 2
      !
      Mf = 0.5_rk * real(Mfx2,rk)
      Mf_ = nint(Mf)
      if (mod(nint(2.0_rk * Jf + 1.0_rk), 2) == 0) Mf_ = int(Mf)
      !
      ! loop initial ang. mom. projection on lab axis
      loop_I : do Mix2 = -nint(2.0_rk * Ji), nint(2.0_rk * Ji), 2
        !
        Mi = 0.5_rk * real(Mix2,rk)
        Mi_ = nint(Mi)
        if (mod(nint(2.0_rk * Ji + 1.0_rk), 2) == 0) Mi_ = int(Mi)
        !
        ! sum only over valid rank 2 spherical tensor components
        ! selection rule doubled to account exactly for 1/2-integer cases
        !
        sphComp = Mfx2 - Mix2
        !
        if(     abs(sphComp) /= 0 &
          .and. abs(sphComp) /= 2 &
          .and. abs(sphComp) /= 4) cycle
        !
        sphComp = nint(0.5_rk * sphComp)
        !
        ! calculation of line strength contribution from second rank component
        f3j = three_j(Ji, Jf, 2.0_rk, -Mi, Mf, Mi - Mf)
        !
        M_2 = (-1.0_rk)**Mi_ * T(iLF, 2, sphComp) * f3j &
            * sqrt((2.0_rk * Ji + 1.0_rk) * (2.0_rk * Jf + 1.0_rk))
        !
        ! if Jf = Ji and Mf = Mi then also count contribution from zeroth rank
        if (Jf == Ji .and. sphComp == 0) then
          !
          f3j = three_j(Ji, Jf, 0.0_rk, -Mi, Mf, Mi - Mf)
          !
          M_0 = (-1.0_rk)**Mi_ * T(iLF, 0, sphComp) * f3j &
              * sqrt((2.0_rk * Ji + 1.0_rk) * (2.0_rk * Jf + 1.0_rk))
        endif
        !
        !
        if (abs(M_2) > small_ .or. abs(M_0) > small_) then
          icount_ = icount_ + 1
          if (.not. present(icount)) write(iunit, "(i4,1x,i4,3x,e24.14,2x,e24.14)") &
            Mi_, Mf_, M_0, M_2
        endif
      enddo loop_I
    enddo loop_F
    !
    if (present(icount)) icount = icount_
    !
    call TimerStop('do_LF_matrix_elements')
    !
  end subroutine do_LF_matrix_elements
  !
  !
  !
  subroutine Energy_filter(Jval,energy,igamma,passed)
    !
    real(rk),intent(in)    :: Jval,energy
    integer(ik),intent(in) :: igamma
    logical,intent(out)    :: passed
    real(rk)               :: erange(2)
    !
    ! passed = .true.
    !
    ! if (.not.intensity%do) return
    !
    passed = .false.
    erange(1) = min(intensity%erange_low(1),intensity%erange_upp(1))
    erange(2) = max(intensity%erange_low(2),intensity%erange_upp(2))
    !
    if (job%isym_do(igamma).and.energy-job%ZPE>=erange(1) .and. &
        Jval>=intensity%J(1).and.Jval<=intensity%J(2).and. &
        energy-job%ZPE<=erange(2)) then 
      !
      passed = .true.
      !
    endif 
    !
  end subroutine Energy_filter
  !
  !
  !
  function three_j0(a,b,c,al,be,ga)
    !
    real(rk) :: three_j0
    real(rk),intent(in) :: a,b,c,al,be,ga
    !
    integer(ik):: newmin,newmax,new,iphase
    real(rk)   :: delta,clebsh,minus
    real(rk)   :: term,term1,term2,term3,summ,dnew,term4,term5,term6

    three_j0=0
    !
    !     (j1+j2).ge.j and j.ge.abs(a-b)    -m=m1+m2    j1,j2,j.ge.0
    !     abs(m1).le.j1    abs(m2).le.j2   abs(m).le.j
    !
    if(c.gt.a+b) return
    if(c.lt.abs(a-b)) return
    if(a.lt.0.or.b.lt.0.or.c.lt.0) return
    if(a.lt.abs(al).or.b.lt.abs(be).or.c.lt.abs(ga)) return
    if(-ga.ne.al+be) return
    !
    !
    ! compute delta(abc)
    !
    delta = sqrt( fakt(a+b-c)*fakt(a+c-b)*fakt(b+c-a)/fakt(a+b+c+1.0_rk) )
    !
    !
    term1=fakt(a+al)*fakt(a-al)
    term2=fakt(b-be)*fakt(b+be)
    term3=fakt(c+ga)*fakt(c-ga)
    term=sqrt( (2.0_rk*c+1.0_rk)*term1*term2*term3 )
    !
    !
    ! now compute summation term
    !
    ! sum to get summation in eq(2.34) of brink and satchler.  
    ! sum until a term inside factorial goes negative.  
    ! new is index for summation
    ! now find what the range of new is.
    !
    newmin=idnint(max((a+be-c),(b-c-al),0.0_rk))
    newmax=idnint(min((a-al),(b+be),(a+b-c)))
    !
    summ=0
    !
    !
    do new=newmin,newmax
      dnew=real(new,rk)
      term4=fakt(a-al-dnew)*fakt(c-b+al+dnew)
      term5=fakt(b+be-dnew)*fakt(c-a-be+dnew)
      term6=fakt(dnew)*fakt(a+b-c-dnew)
      summ=summ+(-1.0_rk)**new/(term4*term5*term6)
    enddo
    !
    ! so clebsch-gordon <j1j2m1m2ljm> is clebsh
    !
    clebsh=delta*term*summ/sqrt(10.0_rk)
    !
    ! convert clebsch-gordon to three_j
    !
    iphase=idnint(a-b-ga)
    minus = -1.0_rk
    if (mod(iphase,2).eq.0) minus = 1.0_rk
    three_j0 = minus * clebsh / sqrt(2.0_rk * c + 1.0_rk)
    !
    !threej=(-1.d0)**(iphase)*clebsh/sqrt(2.0_rk*c+1.d0)
    !
  end function three_j0
  !
  !
  !
  function fakt(a) result (f)

    real(rk),intent(in) :: a
    real(rk)            :: ax,f
    integer(ik)         :: i,ic
    !
    ax=a
    f=1.0_rk
    if(abs(ax)<1.d-24) return
    f=.1_rk
    if(ax.lt.0.d0) then 
       write (*,"(' fkt.err  negative argument for functi on fakt. argument = ',e12.5)") ax
       stop 'fkt.err  negative argument'
    endif 
    !
    ic=nint(ax)
    ax=ax/10.0_rk
    f=ax
    do i = 1, ic - 1
      f=f*(ax-real(i,rk)*0.1_rk)
    enddo
    !
  end function fakt
!
!
!         |     |
!          \   /
!           \_/
!      __   /^\   __
!     '  `. \_/ ,'  `
!          \/ \/     
!     _,--./| |\.--._  
!  _,'   _.-\_/-._  `._
!       |   / \   |
!       |  /   \  |
!      /   |   |   \
!    -'    \___/    `-
!
!
end module polarizability