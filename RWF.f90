module RWF

 use accuracy,     only : hik, ik, rk, ark, cl, out, vellgt, planck, avogno, boltz, pi, small_, aston
 use diatom_module,only : job,Intensity,quantaT,eigen,basis,Ndipoles,dipoletm,duo_j0,fieldT,poten,three_j,jmin_global,&
                          quadrupoletm,nQuadrupoles,grid,kinetic_energy_grid_points,Nestates,Nrefstates,brot,&
                          matrixT,amass,hstep,vibrational_totalroots, vibrational_contrfunc, vibrational_quantum_number
 use timer,        only : IOstart,Arraystart,Arraystop,ArrayMinus,Timerstart,Timerstop,MemoryReport, &
                          TimerReport,memory_limit,memory_now
 use symmetry,     only : sym,correlate_to_Cs
 use lapack,only : lapack_zgesv,lapack_gelss


 implicit none
 !
 private
 public Raman_wavefunction

 real(rk),allocatable,save  :: dipole_me(:, :)


 type DmatT
      real(rk),pointer   :: mat(:,:,:) 
 end type DmatT

 type DimatT
      integer(ik),pointer   :: imat(:,:,:)
 end type DimatT

 type DkmatT
      integer(ik),pointer   :: kmat(:,:)
 end type DkmatT
 !
 type ElevelT
      integer(ik)  :: jind  ! J index
      integer(ik)  :: igamma
      integer(ik)  :: ilevel
 end type ElevelT
 ! 
 type dipoleT                      
    real(rk),pointer    :: rot(:,:,:,:,:)
 end type dipoleT

 type Mat2DT
      real(rk),pointer   :: matelem(:,:) 
 end type Mat2DT


 integer(ik) :: Neigenlevels
 !
 type(ElevelT),allocatable  :: Elevel(:)  

 !
 type(dipoleT),allocatable     :: wigner(:,:) ! Rotational component of the dipole moment matrix elements 

 !
 type(quantaT),allocatable :: quanta_RWF(:)
 !
 !
 !type(IntensityT),save :: intensity

contains
  
 subroutine Raman_wavefunction

 ! compute  partition function, electric dipole transition moments, linestrength, and intensities

    integer(ik)          :: info

    integer(ik)              :: nJ, jind
    real(rk), allocatable :: Jval(:),q_part(:,:)

    real(rk)             :: Jval_,Jval_min,Jmin, Jmax,exp_en, part, beta, energy

    integer(ik)          :: ilevel, irrep,igamma,isym,istate,parity_gu
    integer(ik)          :: iverbose = 4

    ! initialize array of J values
    ! J=0 is always present no matter if we use it in intensity calcs or not. 


    !call define_jlist(Jval)
    !
    Jmin = minval(intensity%J(1:2))
    Jmax = maxval(intensity%J(1:2))
    !
    Jmin = max(jmin_global,Jmin)
    !
    ! This is the lowest possible J, 0 or 0.5 to make sure all computed energies J>0 are included into 
    ! the line list for any intensity%J in order to guarantee consistent Energy list numbering 
    Jval_min = 0 ; if (mod(nint(2.0_rk*job%j_list(1)),2)/=0) Jval_min = 0.5
    Jval_min = max(jmin_global,Jval_min)
    !
    Jval_ = Jval_min
    nJ = 1
    do while (Jval_<Jmax)
      Jval_ = Jval_ + 1.0_rk
      nJ = nJ + 1
    enddo
    !
    allocate(Jval(nJ), stat = info)
    if (info /= 0) stop 'Raman_wavefunction allocation error: Jval - out of memory'
    !
    allocate(q_part(20,nJ), stat = info)
    if (info /= 0) stop 'Raman_wavefunction allocation error: q_part - out of memory'
    !
    Jval_ = Jval_min
    jind = 1
    Jval(jind) = Jval_
    do while (Jval_<Jmax)
       jind = jind + 1
       Jval_ = Jval_ + 1.0_rk
       Jval(jind) = Jval_
    end do
    !
    call duo_j0(iverbose,Jval)
    !
    !call Sort_levels(iverbose,nJ, Jval(1:nJ))
    !
    ! update ZPE if not provided as part of the intensity input
    !
    if (job%shift_to_zpe) then
      !
      do jind = 1,nJ
        !
        Jval_ = Jval(jind)
        !
        do igamma = 1,sym%NrepresCs
          !
          do ilevel = 1,eigen(jind,igamma)%Nlevels
              !
              energy = eigen(jind,igamma)%val(ilevel)
              !
              intensity%zpe = min(intensity%zpe,energy)
              !
           enddo
           !
        enddo
      end do
      !
      if (iverbose>=2) write(out,'(/"Zero point energy (ZPE) = ",f18.6," (global zero, used for intensities)")') intensity%ZPE
      if (iverbose>=4) write(out,"(/'Partition funciton = ',f18.4,' T = ',f12.2)") intensity%part_func,intensity%temperature
      !
    endif
    !
    select case (trim(intensity%action))
    !
    case('ABSORPTION', 'EMISSION', 'TM')
       !
       ! read eigenvalues and their labeling, i.e. description;
       ! initialize file-units for reading eigenvectors
       !
       call Sort_levels(iverbose,nJ, Jval(1:nJ))
       !
       !restore vibrational contracted matrix elements
       !for all  dipole moment vector components
       !
       !call check_point_dipoles('READ',iverbose,totalroots)
       !
       ! Run intensity simulations 
       !
       beta = planck * vellgt / (boltz * intensity%temperature)
       !
       ! Compute part-func if not given
       !
       if (intensity%part_func<small_) then 
          !
          intensity%part_func  = 0 
          !
          do jind = 1,nJ
            !
            Jval_ = Jval(jind)
            !
            do igamma = 1,sym%NrepresCs
              !
              do ilevel = 1,eigen(jind,igamma)%Nlevels
                  !
                  energy = eigen(jind,igamma)%val(ilevel)
                  irrep  = eigen(jind,igamma)%quanta(ilevel)%igamma
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
          end do
          !
          if (iverbose>=4) write(out,"(/'Partition funciton = ',f18.4,' T = ',f12.2)") intensity%part_func,intensity%temperature
          !
       endif
       !
       !call rwf_intensity(Jval,iverbose)
       !
       call rwf_dvr_intensity(Jval,iverbose)
       !
       write(out, '(/a)') 'done'
       !
    !compute partition function
    !
    case('PARTFUNC')
       !
       write(out, '(/a)') 'compute partition function'
       !
       beta = planck * vellgt / (boltz * intensity%temperature)
       !
       !loop ove J values
       !
       q_part  = 0 
       !
       do jind = 1,nJ
         !
         Jval_ = Jval(jind)
         !
         do igamma = 1,sym%NrepresCs
           !
           do ilevel = 1,eigen(jind,igamma)%Nlevels
               !
               energy = eigen(jind,igamma)%val(ilevel)
               irrep  = eigen(jind,igamma)%quanta(ilevel)%igamma
               !
               ! for homonuclear symmetries Nrepres = 4 and the irrep can be reconstucted baed on the parity and g/u: 
               istate  = eigen(jind,igamma)%quanta(ilevel)%istate
               parity_gu = poten(istate)%parity%gu
               isym = correlate_to_Cs(igamma,parity_gu)
               !
               exp_en = exp(-(energy - intensity%ZPE) * beta)
               !
               q_part(irrep,jind) = q_part(irrep,jind) + intensity%gns(isym)*(2.0_rk*Jval_+1.0_rk) * exp_en
               !
            enddo
            !
         enddo
       end do
       !
       part = sum(q_part(:,:))
       !
       do jind=1,nJ
         do irrep = 1,sym%NrepresCs
            write(out, '(i4,1x,f18.1,1x,es16.8)') irrep,Jval(jind),q_part(irrep,jind)
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
 end subroutine Raman_wavefunction
 !

 !
 function cg(j0, k0, dj, dk)

 ! return values for some CG coeeficients
 ! j0      and k0      are j" and k" (initial state)
 ! j0 + dj and k0 + dk are j' and k' (final state)
 ! compute <j0 k0, 1 dk|j0+dj k0+dk>

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
 ! Electric dipole Raman-Wave-Function intensities and cross section calculations 
 !
 subroutine rwf_intensity(Jval,iverbose)
    !
    implicit none
    !
    real(rk),intent(in)  :: Jval(:)
    integer(ik),intent(in)   :: iverbose

    integer(ik)    :: nJ,dimenmax
    integer(ik)    :: ilevelI, ilevelF,ilevelR
    integer(ik)    :: nlevelsG(sym%Nrepresen)
    integer(ik)    :: info,indI,indF,itransit,Ntransit,Nrepresen
    integer(ik)    :: igammaI,igammaF
    integer(ik)    :: dimenI,dimenF,nmax,parity_gu,isymI,isymF
    real(rk)       :: energyI, energyF,energyR,nu_if,linestr,ener_
    real(rk)       :: jI,jF,ddot
    logical        :: passed,passed_

    real(rk),allocatable :: vecI(:), vecF(:)
    real(rk),allocatable :: half_linestr(:),half_pecme(:)
    !
    integer(ik)  :: jind,nlevels
    !
    integer(ik)  :: iroot,NlevelsI,NlevelsF,nlower
    !
    integer(ik)  :: igamma_pair(sym%Nrepresen),igamma,istateI,istateF,ivibI,ivibF,ivI,ivF,ilambdaI,ilambdaF
    real(rk)     :: spinI,spinF,omegaI,omegaF,sigmaI,sigmaF
    integer(hik) :: matsize
    !
    character(len=1) :: branch
    !
    type(quantaT),pointer  :: quantaI,quantaF
    !
    real(rk)     :: boltz_fc, beta, intens_cm_mol, emcoef, A_coef_s_1
    !
    character(len=130) :: my_fmt !format for I/O specification
    integer(ik)  :: transunit
    character(len=cl) :: filename,ioname
    !
    logical     :: integer_spin = .true.
    !
    integer(ik) :: alloc_p
    !
    integer(ik) :: inu
    real(rk) :: J_
    real(rk) :: dnu, nu,RWF2,intens_cm_molecule
    complex(rk),allocatable :: Amat(:,:),B(:,:)
    !
    real(rk),allocatable :: crosssections(:)
    !
    type(Mat2DT) :: mu,pec
    !
    call TimerStart('Intensity calculations')
    !
    if (sym%maxdegen>2) then 
      !
      write(out,"('dm_intensity: this procedure has not been tested for the symmetries with degeneracies higher than 2')")
      !stop 'dm_intensity was not tested for present symmetry'
      !
    endif
    !
    Nrepresen = sym%NrepresCs
    !
    beta = planck * vellgt / (boltz * intensity%temperature)
    intens_cm_mol  = 8.0d-36*pi**3 / (3.0_rk * planck * vellgt)
    emcoef = planck*vellgt/(4.0_rk*pi)
    A_coef_s_1     =64.0d-36 * pi**4  / (3.0_rk * planck)
    !
    nJ = size(Jval)
    !
    ! Prepare the list of units with the stored eigenvectors
    !
    !allocate(Jeigenvec_unit(nJ), stat = info)
    !if (info /= 0) stop 'dm_tranint allocation error: Jeigenvec_unit - out of memory'
    !
    !do jind=1,nJ 
    !  Jeigenvec_unit(jind) = TReigenvec_unit(jind,Jval)
    !enddo
    !
    ! open units for the line list in the exomol format
    if (trim(intensity%linelist_file)/="NONE") then
      !
      filename =  trim(intensity%linelist_file)//'.xsec'
      write(ioname, '(a, i4)') 'cross sections '
      call IOstart(trim(ioname),transunit)
      open(unit = transunit, action = 'write',status='replace' , file = filename)
      !
    endif
    !
    ! maximal size of basis functions 
    !
    dimenmax = 0 
    !
    !loop over J quantities
    !
    do jind = 1, nJ
      do igamma=1,Nrepresen
        !
        ! Estimate the maximal size of the basis functions. 
        !
        dimenmax = max(dimenmax,eigen(jind,igamma)%Ndimen)
        !
      enddo
    enddo 
    !
    Nmax = nint(Jval(nJ))+1
    !
    ! in the case of the proper rotatonal symmetry the matrix elements of the 
    ! wigner matrix have to be computed using the correspoding symmetrization 
    ! transformation derived here. they will be used together with the 
    ! vibrational matrix elements of the dipole moment when evaluating 
    ! line strengthes using eigenvectos.  
    !
    !if (job%rotsym_do) call conctraced_rotational_dipole(nJ, jval, jmax, threej)
    !
    !allocate arrays for eigenvectors and coefficients
    !
    ! First we count transitions to be calculated. It will help us to keep track of the 
    ! calculation progress.
    !
    Ntransit = 0 
    !
    !number of initial states
    !
    nlevels = Neigenlevels
    !
    ! For a given symmetry igamma with some gns(igamma) we find its counterpart jgamma/=igamma
    ! having the same gns(jgamma). We assume that there is only one such pair 
    ! in case of absorption or emission calcs. 
    !
    call find_igamma_pair(igamma_pair)
    !
    call TimerStart('Intens_Filter-1')
    !
    ! loop over initial states
    !
    ener_ = 0
    !
    nlower = 0
    iroot = 0
    !
    do indI = 1, nJ
       !
       ! rotational quantum number 
       !
       jI = Jval(indI)
       !
       do igammaI=1,Nrepresen
         !
         nlevelsI = eigen(indI,igammaI)%Nlevels
         !
         !omp parallel do private(ilevelI,jI,energyI,igammaI,quantaI,ilevelF,jF,energyF,igammaF,quantaF,passed) & 
         !                      & schedule(guided) reduction(+:Ntransit,nlevelI)
         do ilevelI = 1, nlevelsI
           !
           !energy energy and and quanta of the initial state
           !
           energyI = eigen(indI,igammaI)%val(ilevelI)
           !
           istateI  = eigen(indI,igammaI)%quanta(ilevelI)%istate
           parity_gu = poten(istateI)%parity%gu
           isymI = correlate_to_Cs(igammaI,parity_gu)
           !
           call energy_filter_lower(jI,energyI,passed)
           !
           if (.not.passed) cycle
           !
           nlower = nlower + 1
           !
           do indF = 1, nJ
              !
              ! rotational quantum number 
              !
              jF = jval(indF)
              !
              do igammaF=1,Nrepresen
                !
                nlevelsF = eigen(indF,igammaF)%Nlevels
                !
                !call Jgamma_filter(jI,jF,igammaI,igammaF,igamma_pair,passed)
                !
                !if (.not.passed) cycle
                !
                !omp parallel do private(ilevelI,jI,energyI,igammaI,quantaI,ilevelF,jF,energyF,igammaF,quantaF,passed) & 
                !                        & schedule(guided) reduction(+:Ntransit,nlevelI)
                do ilevelF = 1, nlevelsF
                  !
                  !energy and and quanta of the final state
                  !
                  energyF = eigen(indF,igammaF)%val(ilevelF)
                  !
                  istateF  = eigen(indF,igammaF)%quanta(ilevelF)%istate
                  parity_gu = poten(istateF)%parity%gu
                  isymF = correlate_to_Cs(igammaF,parity_gu)
                  !
                  call intens_filter(jI,jF,energyI,energyF,isymI,isymF,igamma_pair,passed)
                  !
                  if (passed) then 
                    !
                    Ntransit = Ntransit + 1
                    !
                  endif 
                  !
                enddo
              enddo
           enddo
           !
         enddo
         !
       enddo
    enddo
    !omp end parallel do
    !
    call TimerStop('Intens_Filter-1')
    !
    !
    !loop over final states -> count states for each symmetry
    !
    nlevelsG = 0
    !
    ! guessing if this is a half-integer spin case
    if (mod(eigen(1,1)%quanta(1)%imulti,2)==0) integer_spin = .false.
    !
    allocate(vecI(dimenmax), stat = info)
    call ArrayStart('intensity-vecI',info,size(vecI),kind(vecI))
    !
    do indI = 1, nJ
       !
       ! rotational quantum number 
       !
       jI = jval(indI)
       !
       J_ = -1 
       !
       do igammaI=1,Nrepresen
         !
         nlevelsI = eigen(indI,igammaI)%Nlevels
         !
         dimenI = eigen(indI,igammaI)%Ndimen
         !
         !
         !omp parallel do private(ilevelI,jI,energyI,igammaI,quantaI,ilevelF,jF,energyF,igammaF,quantaF,passed) & 
         !                      & schedule(guided) reduction(+:Ntransit,nlevelI)
         do ilevelI = 1, nlevelsI
           !
           !energy energy and and the symmetry of the state
           !
           energyI = eigen(indI,igammaI)%val(ilevelI)
           !
           istateI  = eigen(indI,igammaI)%quanta(ilevelI)%istate
           parity_gu = poten(istateI)%parity%gu
           ! Obtain the C2v/Cs symmetry 
           isymI = correlate_to_Cs(igammaI,parity_gu)
           !
           ! ignore states with zero nuclear weight 
           if (intensity%gns(isymI)<small_) cycle 
           !
           iroot = iroot + 1
           eigen(indI,igammaI)%quanta(ilevelI)%iroot = iroot
           !
           call energy_filter_upper(jI,energyI,passed)
           !
           call energy_filter_lower(jI,energyI,passed_)
           !
           if (.not.passed.and..not.passed_) cycle
           !
           istateI  = eigen(indI,igammaI)%quanta(ilevelI)%istate
           parity_gu = poten(istateI)%parity%gu
           isymI = correlate_to_Cs(igammaI,parity_gu)
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
    write(my_fmt,'(A,I0,A)') "('Number of states for each symm = ',", sym%Nrepresen, "i8)"
    write(out,my_fmt) nlevelsG(:)
    !
    matsize = int(sum( nlevelsG(:) ),hik)
    !
    if (iverbose>=4) write(out,"(/'Dipole moment integration (i)...')")
    !
    ! prescreen all eigenfunctions, compact and store on the disk
    !
    !call TimerStart('Dipole moment integration (i)')
    !
    !allocate(vecI(sym%Maxdegen,dimenmax_swap),vecF(sym%Maxdegen,dimenmax), stat = info)
    !if (info/=0)  stop 'vecI,vecF,icoeffF - out of memory'
    !
    !allocate(icoeffF(sym%Maxdegen,dimenmax), stat = info)
    !
    if (Ntransit==0) then 
         write(out,"('dm_intensity: the transition filters are too tight: no entry')") 
         stop 'dm_intensity: the filters are too tight' 
    endif 
    !
    !call TimerStop('Dipole moment integration (i)')
    !
    write(out,"(/'...done!')")
    !
    allocate(vecI(dimenmax),vecF(dimenmax), stat = info)
    !
    call ArrayStart('intensity-vectors',info,size(vecI),kind(vecI))
    call ArrayStart('intensity-vectors',info,size(vecF),kind(vecF))
    !
    ! loop over final states -> count states for each symmetry
    !
    write(my_fmt,'(A,I0,A)') "('Number of states for each symm = ',", sym%Nrepresen, "i8)"
    write(out,my_fmt) nlevelsG(:)
    !
    if (iverbose >= 0) then
       write(out,"(' Total number of lower states = ',i8)") nlower
       write(out,"(' Total number of transitions  = ',i8)") Ntransit
    end if
    !
    if (iverbose >= 0) then
       write(out,"(/' Statistical weights gns = ',4f12.1)") intensity%gns(1:)
    end if
    !
    ! In order to speed up the line strength evaluation, 
    ! S_if = | <i|mu|f> |^2  = | \sum_{nm} C_in C_fm <n|mu|m> |^2
    ! in three steps:
    ! 1. Evaluating the expansion of the i-state 
    !    s_{im} =  \sum_{n} C_in  <n|mu|m> 
    ! 2. performing the expansion of the final state f:
    !    s_{if} =  \sum_{m} C_fm  s_{im}. 
    ! 3. Building S_{if}
    !    S_{if} = s_{if}^2
    !
    !  The temporaly object s_{im} will be referted to as 
    !  a half-linestrength "half_linestr"
    !
    allocate(half_linestr(dimenmax),half_pecme(dimenmax),stat=info)
    !
    call ArrayStart('half_linestr',info,size(half_linestr),kind(half_linestr))
    call ArrayStart('half_pecme',info,size(half_pecme),kind(half_pecme))
    !
    !  The matrix where some of the eigenvectors will be stored
    !
    if (iverbose>=5) call MemoryReport
    !
    deallocate(vecF)
    !
    !tid = omp_get_thread_num()
    !
    !omp parallel private(tid)
    !  if (tid==0) then
    !     nprocs_ = omp_get_num_threads()
    !  endif
    !omp end parallel
    !
    ! ---------------------------------
    ! Dipole moment matrix elements
    ! --------------------------------- 
    !
    itransit = 0
    !
    if (intensity%npoints==-1) then
      write(out,"('RWF error: cross sectons Npoints is undefined. Define it using NPOINTS keyword in INTENSITY')")
      stop 'RWF error: cross sectons Npoints is undefined'
    endif
    !
    dnu = (intensity%freq_window(2)-intensity%freq_window(1))/real(intensity%Npoints-1,rk)
    !
    allocate(crosssections(intensity%npoints),stat = info)
    call ArrayStart('crosssections',info,size(crosssections),kind(crosssections))
    crosssections = 0
    !
    intens_cm_molecule  = 8.0d-36 * pi**3/ (3.0_rk * planck * vellgt)
    !
    ! loop over initial states
    !
    do indI = 1, nJ
       !
       jI = jval(indI)
       !
       if (iverbose>=3) write(out,"('J = ',f5.1)") jI
       !
       do igammaI=1,Nrepresen
         !
         nlevelsI = eigen(indI,igammaI)%Nlevels
         dimenI = eigen(indI,igammaI)%Ndimen
         !
         if (nlevelsI==0) cycle 
         parity_gu = poten(istateI)%parity%gu
         isymI = correlate_to_Cs(igammaI,parity_gu)
         !
         do indF = 1, nJ
           !
           jF = jval(indF)
           if (abs(nint(jI-jF))>1.or.abs(nint(jI+jF))==0) cycle 
           !
           do igammaF=1,Nrepresen
              !
              !call Jgamma_filter(jI,jF,igammaI,igammaF,igamma_pair,passed)
              !if (.not.passed) cycle
              !
              nlevelsF = eigen(indF,igammaF)%Nlevels
              dimenF = eigen(indF,igammaF)%Ndimen
              !
              !igammaF = igamma_pair(igammaI)
              !
              parity_gu = poten(istateF)%parity%gu
              isymF = correlate_to_Cs(igammaF,parity_gu)
              !
              if (isymF /= igamma_pair(isymI)) cycle
              !
              nlevelsF = eigen(indF,igammaF)%Nlevels
              !
              if (nlevelsF==0) cycle 
              !
              allocate(mu%matelem(nlevelsF,nlevelsI),stat = info)
              call ArrayStart('mu%matelem',info,size(mu%matelem),kind(mu%matelem))
              mu%matelem = 0
              !
              allocate(pec%matelem(nlevelsF,nlevelsF),stat = info)
              call ArrayStart('pec%matelem',info,size(pec%matelem),kind(pec%matelem))
              pec%matelem = 0
              !
              allocate(Amat(nlevelsF,nlevelsF),B(nlevelsF,1),stat = info)
              call ArrayStart('RWF:Amat',info,size(Amat),kind(Amat))
              call ArrayStart('RWF:Amat',info,size(B),kind(B))
              !
              !
              do ilevelF = 1, nlevelsF
                !
                !energy and and quanta of the final state
                !
                energyF = eigen(indF,igammaF)%val(ilevelF)
                !
                call energy_filter_upper(jF,energyF,passed)
                !
                if (.not.passed) cycle
                !
                vecI(1:dimenF) = eigen(indF,igammaF)%vect(1:dimenF,ilevelF)
                !
                half_pecme = 0
                !
                call do_matelem_pec(jF,jF,indF,indF,dimenF,dimenF,&
                                                 vecI(1:dimenF),&
                                                 half_pecme)
                !
                do ilevelR = 1,nlevelsF
                  !
                  energyR = eigen(indF,igammaF)%val(ilevelR)
                  !
                  call energy_filter_upper(jF,energyR,passed)
                  !
                  if (.not.passed) cycle
                  !
                  vecI(1:dimenF) = eigen(indF,igammaF)%vect(1:dimenF,ilevelF)
                  !
                  if (nint(jF-jI)==0)  then
                    pec%matelem(ilevelF,ilevelR) = ddot(dimenF,half_pecme,1,vecI,1)
                  endif
                  !
                enddo
                !
              enddo
              !
              !omp parallel do private(ilevelI,jI,energyI,igammaI,quantaI,ilevelF,jF,energyF,igammaF,quantaF,passed) & 
              !                        & schedule(guided) reduction(+:Ntransit,nlevelI)
              Ilevels_loop : do ilevelI = 1, nlevelsI
                !
                !energy and and quanta of the final state
                !
                energyI = eigen(indI,igammaI)%val(ilevelI)
                !
                !dimension of the bases for the initial states
                !
                !energy, quanta, and gedeneracy order of the initial state
                quantaI => eigen(indI,igammaI)%quanta(ilevelI)
                istateI  = quantaI%istate
                ivibI    = quantaI%ivib
                ivI      = quantaI%v
                sigmaI   = quantaI%sigma
                spinI    = quantaI%spin
                ilambdaI = quantaI%ilambda
                omegaI   = quantaI%omega
                !
                ! reconstruct the symmetry for the C2v case which is different from Cs
                parity_gu = poten(istateI)%parity%gu
                isymI = correlate_to_Cs(igammaI,parity_gu)
                !
                call energy_filter_lower(jI,energyI,passed)
                !
                if (.not.passed) cycle
                !
                vecI(1:dimenI) = eigen(indI,igammaI)%vect(1:dimenI,ilevelI)
                !
                ! Compute the half-linestrength
                !
                half_linestr = 0
                !
                if (isymF /= igamma_pair(isymI)) cycle
                !
                if (( intensity%J(1)+intensity%J(2)>0 )&
                    .and. abs(nint(jI-jF))<=1.and.nint(jI+jF)>=1) then 
                   !
                   call do_1st_half_linestrength(jI,jF,indI,indF,dimenI,dimenF,&
                                                 vecI(1:dimenI),&
                                                 half_linestr)
                   !
                endif
                !
                !loop over final states
                !
                !omp parallel private(vecF,alloc_p)
                allocate(vecF(dimenmax),stat = alloc_p)
                if (alloc_p/=0) then
                    write (out,"(' dipole: ',i9,' trying to allocate array -vecF')") alloc_p
                    stop 'dipole-vecF - out of memory'
                end if
                !
                !omp do private(ilevelF,energyF,dimenF,quantaF,istateF,ivibF,ivF,sigmaF,spinF,ilambdaF,omegaF,passed,&
                !omp& parity_gu,isymF,branch,nu_if,linestr,linestr2,A_einst,boltz_fc,absorption_int,tm) schedule(static) &
                !omp                                                                             & reduction(+:itransit)
                Flevels_loop: do ilevelF = 1,nlevelsF
                   !
                   !energy and and quanta of the final state
                   !
                   energyF = eigen(indF,igammaF)%val(ilevelF)
                   !
                   !dimension of the bases for the final state
                   !
                   dimenF = eigen(indF,igammaF)%Ndimen
                   !
                   quantaF => eigen(indF,igammaF)%quanta(ilevelF)
                   !
                   istateF  = quantaF%istate
                   ivibF    = quantaF%ivib
                   ivF      = quantaF%v
                   sigmaF   = quantaF%sigma
                   spinF    = quantaF%spin
                   ilambdaF = quantaF%ilambda
                   omegaF   = quantaF%omega
                   !
                   call energy_filter_upper(jF,energyF,passed)
                   !
                   if (.not.passed) cycle Flevels_loop
                   !
                   parity_gu = poten(istateF)%parity%gu
                   isymF = correlate_to_Cs(igammaF,parity_gu)
                   !
                   !call TimerStart('Intens_Filter-3')
                   !
                   call intens_filter(jI,jF,energyI,energyF,isymI,isymF,igamma_pair,passed)
                   !
                   !call TimerStop('Intens_Filter-3')
                   !
                   if (.not.passed) cycle Flevels_loop
                   !
                   ! Which PQR branch this transition belong to ->
                   ! 
                   branch = PQR_branch(jI,jF)
                   !
                   nu_if = energyF - energyI 
                   !
                   ! no zero-frequency transitions should be produced 
                   if ( nu_if < small_) cycle
                   !if (trim(intensity%action)=='EMISSION') nu_if = -nu_if 
                   !
                   ! Count the processed transitions 
                   !
                   itransit = itransit + 1
                   !
                   vecF(1:dimenF) = eigen(indF,igammaF)%vect(1:dimenF,ilevelF)
                   !
                   linestr = ddot(dimenF,half_linestr,1,vecF,1)
                   !
                   mu%matelem(ilevelF,ilevelI) = linestr
                   !
                end do Flevels_loop
                !omp enddo
                !
                deallocate(vecF)
                !omp end parallel
                !
                if (iverbose>=5) call TimerReport
                !
              enddo Ilevels_loop
              !
              !
              ! Wavenumber grid 
              !
              do inu = 1,intensity%npoints
                 !
                 nu = intensity%freq_window(1)+dnu*real(inu,rk)
                 !
                 if (iverbose>=4.and.mod(inu,50)==0) write(out,"('nu = ',f9.2)") nu
                 !
                 do ilevelI = 1, nlevelsI
                   !
                   energyI = eigen(indI,igammaI)%val(ilevelI)
                   !
                   Amat = 0 
                   B = 0
                   !
                   do ilevelF = 1, nlevelsF
                     !
                     !energy and and quanta of the final state
                     !
                     energyF = eigen(indF,igammaF)%val(ilevelF)
                     !
                     call energy_filter_upper(jF,energyF,passed)
                     !
                     if (.not.passed) cycle
                     !
                     B(ilevelF,1) = cmplx(0.0_rk,mu%matelem(ilevelF,ilevelI))
                     !
                     do ilevelR = 1,nlevelsF
                       !
                       energyR = eigen(indF,igammaF)%val(ilevelR)
                       !
                       call energy_filter_upper(jF,energyR,passed)
                       !
                       if (.not.passed) cycle
                       !
                       Amat(ilevelF,ilevelR) = -pec%matelem(ilevelF,ilevelR)
                       !
                       if (ilevelF==ilevelR) then
                         !
                         Amat(ilevelF,ilevelR) = Amat(ilevelF,ilevelR) + nu + energyI - energyR &
                            + cmplx(0.0_rk,intensity%gamma,kind=rk) 
                         !
                       endif
                       !
                       !endif
                       !
                     enddo
                     !
                   enddo
                   !
                   call lapack_gelss(Amat,b)
                   !
                   RWF2 = sum(conjg(b)*b)
                   boltz_fc = intensity%gns(isymI)*real( (2*jI + 1)*(2 * jF + 1),rk )*nu *&
                        exp(-(energyI-intensity%ZPE) * beta) * (1.0_rk - exp(-nu * beta))/intensity%part_func
                   !
                   crosssections(inu) = crosssections(inu) + intens_cm_molecule*boltz_fc*RWF2
                   !
                 enddo
                 !
              enddo
              !
              deallocate(mu%matelem,stat = info)
              call ArrayStop('mu%matelem')
              !
              deallocate(pec%matelem,stat = info)
              call ArrayStop('pec%matelem')
              !
              deallocate(Amat,B)
              call Arraystop('RWF:Amat')
              !
           enddo
           !
         enddo
         !
       enddo
       !
    enddo
    !
    do inu = 1,intensity%npoints
       !
       nu = intensity%freq_window(1)+dnu*real(inu,rk)
       !
       write(transunit,"(f12.5,2x,e12.5)") nu,crosssections(inu)
       !
    enddo
    !
    deallocate(crosssections)
    call ArrayStop('crosssections')
    !
    deallocate(vecI)
    call ArrayStop('intensity-vectors')
    !
    deallocate(half_linestr)
    call ArrayStop('half_linestr')
    !
    deallocate(half_pecme)
    call ArrayStop('half_pecme')
    !
    if (trim(intensity%linelist_file)/="NONE") close(transunit,status="keep")
    !
    call TimerStop('Intensity calculations')
    !
  end subroutine rwf_intensity
  !


 !
 ! Electric dipole Raman-Wave-Function intensities and cross section calculations 
 !
 subroutine rwf_dvr_intensity(Jval,iverbose)
    !
    implicit none
    !
    real(rk),intent(in)  :: Jval(:)
    integer(ik),intent(in)   :: iverbose

    integer(ik)    :: nJ,dimenmax
    integer(ik)    :: ilevelI, ilevelF,ilevelR
    integer(ik)    :: nlevelsG(sym%Nrepresen)
    integer(ik)    :: info,indI,indF,itransit,Ntransit,Nrepresen
    integer(ik)    :: igammaI,igammaF
    integer(ik)    :: dimenI,dimenF,nmax,parity_gu,isymI
    real(rk)       :: energyI
    real(rk)       :: jI,jF
    logical        :: passed,passed_

    real(rk),allocatable :: vecI(:), vecF(:)
    real(rk),allocatable :: half_linestr(:)
    !
    integer(ik)  :: jind,nlevels
    !
    integer(ik)  :: iroot,NlevelsI,NlevelsF,nlower
    !
    integer(ik)  :: igamma_pair(sym%Nrepresen),igamma,istateI,ivibI,ivI,ilambdaI
    integer(ik)  :: ngrid
    real(rk)     :: spinI,omegaI,sigmaI
    integer(hik) :: matsize
    !
    type(quantaT),pointer  :: quantaI
    !
    real(rk)     :: boltz_fc, beta, intens_cm_mol, emcoef, A_coef_s_1
    !
    character(len=130) :: my_fmt !format for I/O specification
    integer(ik)  :: enunit,transunit
    character(len=cl) :: filename,ioname
    !
    logical     :: integer_spin = .true.
    !
    integer(ik) :: alloc_p
    !
    integer(ik) :: inu,Nlambdasigmas,i,ilevel,ivib,Ntotal
    real(rk) :: J_,delta
    real(rk) :: dnu, nu,RWF2,intens_cm_molecule,sc,scale,h12
    complex(16),allocatable :: Amat(:,:),B(:),C(:)
    complex(16),parameter :: alpha_ = (1.0d0,0.0d0),beta_ = (0.0d0,0.0d0)
    double precision,parameter :: dalpha = 1.0d0 , dbeta = 0.d0
    !
    real(rk),allocatable :: crosssections(:)
    !
    real(rk),allocatable :: kinmat(:,:),hmat(:,:),dipole_mat(:,:),hmat_n(:,:,:)
    !
    integer(ik) :: istate,imulti,ilambda,igrid,j,jvib,jlevel,jstate,jmulti,jlambda,Nsym(2),Nvib_l
    real(rk) :: sigma,omega,f_rot,sigmaj,spinj,omegaj,erot
    type(quantaT),allocatable :: icontr(:)
    integer(ik),allocatable :: Nirr(:,:),ilevel2i(:,:)
    complex(16) :: zdotc
    !double precision :: ddot
    !
    call TimerStart('Intensity calculations')
    !
    if (sym%maxdegen>2) then 
      !
      write(out,"('dm_intensity: this procedure has not been tested for the symmetries with degeneracies higher than 2')")
      !stop 'dm_intensity was not tested for present symmetry'
      !
    endif
    !
    Nrepresen = sym%NrepresCs
    !
    beta = planck * vellgt / (boltz * intensity%temperature)
    intens_cm_mol  = 8.0d-36*pi**3 / (3.0_rk * planck * vellgt)
    emcoef = planck*vellgt/(4.0_rk*pi)
    A_coef_s_1     =64.0d-36 * pi**4  / (3.0_rk * planck)
    !
    scale = amass/aston
    !
    h12 = 12.0_rk*hstep**2
    sc  = h12*scale
    !
    nJ = size(Jval)
    !
    ! Prepare the list of units with the stored eigenvectors
    !
    !allocate(Jeigenvec_unit(nJ), stat = info)
    !if (info /= 0) stop 'dm_tranint allocation error: Jeigenvec_unit - out of memory'
    !
    !do jind=1,nJ 
    !  Jeigenvec_unit(jind) = TReigenvec_unit(jind,Jval)
    !enddo
    !
    ! open units for the line list in the exomol format
    if (trim(intensity%linelist_file)/="NONE") then
      !
      filename =  trim(intensity%linelist_file)//'.xsec'
      write(ioname, '(a, i4)') 'cross sections '
      call IOstart(trim(ioname),transunit)
      open(unit = transunit, action = 'write',status='replace' , file = filename)
      !
    endif
    !
    ! maximal size of basis functions 
    !
    dimenmax = 0 
    !
    !loop over J quantities
    !
    do jind = 1, nJ
      do igamma=1,Nrepresen
        !
        ! Estimate the maximal size of the basis functions. 
        !
        dimenmax = max(dimenmax,eigen(jind,igamma)%Ndimen)
        !
      enddo
    enddo 
    !
    Nmax = nint(Jval(nJ))+1
    !
    ! in the case of the proper rotatonal symmetry the matrix elements of the 
    ! wigner matrix have to be computed using the correspoding symmetrization 
    ! transformation derived here. they will be used together with the 
    ! vibrational matrix elements of the dipole moment when evaluating 
    ! line strengthes using eigenvectos.  
    !
    !if (job%rotsym_do) call conctraced_rotational_dipole(nJ, jval, jmax, threej)
    !
    !allocate arrays for eigenvectors and coefficients
    !
    ! First we count transitions to be calculated. It will help us to keep track of the 
    ! calculation progress.
    !
    Ntransit = 0 
    !
    !number of initial states
    !
    nlevels = Neigenlevels
    !
    ! For a given symmetry igamma with some gns(igamma) we find its counterpart jgamma/=igamma
    ! having the same gns(jgamma). We assume that there is only one such pair 
    ! in case of absorption or emission calcs. 
    !
    call find_igamma_pair(igamma_pair)
    !
    !loop over final states -> count states for each symmetry
    !
    nlevelsG = 0
    !
    ! guessing if this is a half-integer spin case
    if (mod(eigen(1,1)%quanta(1)%imulti,2)==0) integer_spin = .false.
    !
    allocate(vecI(dimenmax), stat = info)
    call ArrayStart('intensity-vecI',info,size(vecI),kind(vecI))
    !
    do indI = 1, nJ
       !
       ! rotational quantum number 
       !
       jI = jval(indI)
       !
       J_ = -1 ! For RichMol count
       !
       do igammaI=1,Nrepresen
         !
         nlevelsI = eigen(indI,igammaI)%Nlevels
         !
         dimenI = eigen(indI,igammaI)%Ndimen
         !
         !
         !omp parallel do private(ilevelI,jI,energyI,igammaI,quantaI,ilevelF,jF,energyF,igammaF,quantaF,passed) & 
         !                      & schedule(guided) reduction(+:Ntransit,nlevelI)
         do ilevelI = 1, nlevelsI
           !
           !energy energy and and the symmetry of the state
           !
           energyI = eigen(indI,igammaI)%val(ilevelI)
           !
           istateI  = eigen(indI,igammaI)%quanta(ilevelI)%istate
           parity_gu = poten(istateI)%parity%gu
           ! Obtain the C2v/Cs symmetry 
           isymI = correlate_to_Cs(igammaI,parity_gu)
           !
           ! ignore states with zero nuclear weight 
           if (intensity%gns(isymI)<small_) cycle 
           !
           iroot = iroot + 1
           eigen(indI,igammaI)%quanta(ilevelI)%iroot = iroot
           !
           call energy_filter_upper(jI,energyI,passed)
           !
           call energy_filter_lower(jI,energyI,passed_)
           !
           if (.not.passed.and..not.passed_) cycle
           !
           istateI  = eigen(indI,igammaI)%quanta(ilevelI)%istate
           parity_gu = poten(istateI)%parity%gu
           isymI = correlate_to_Cs(igammaI,parity_gu)
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
    if (trim(intensity%linelist_file)/="NONE") close(enunit,status='keep')
    !
    write(my_fmt,'(A,I0,A)') "('Number of states for each symm = ',", sym%Nrepresen, "i8)"
    write(out,my_fmt) nlevelsG(:)
    !
    matsize = int(sum( nlevelsG(:) ),hik)
    !
    if (iverbose>=4) write(out,"(/'Dipole moment integration (i)...')")
    !
    ! prescreen all eigenfunctions, compact and store on the disk
    !
    !call TimerStart('Dipole moment integration (i)')
    !
    !allocate(vecI(sym%Maxdegen,dimenmax_swap),vecF(sym%Maxdegen,dimenmax), stat = info)
    !if (info/=0)  stop 'vecI,vecF,icoeffF - out of memory'
    !
    !allocate(icoeffF(sym%Maxdegen,dimenmax), stat = info)
    !
    !call TimerStop('Dipole moment integration (i)')
    !
    write(out,"(/'...done!')")
    !
    allocate(vecI(dimenmax),vecF(dimenmax), stat = info)
    !
    call ArrayStart('intensity-vectors',info,size(vecI),kind(vecI))
    call ArrayStart('intensity-vectors',info,size(vecF),kind(vecF))
    !
    ! loop over final states -> count states for each symmetry
    !
    write(my_fmt,'(A,I0,A)') "('Number of states for each symm = ',", sym%Nrepresen, "i8)"
    write(out,my_fmt) nlevelsG(:)
    !
    if (iverbose >= 0) then
       write(out,"(' Total number of lower states = ',i8)") nlevelsI
       write(out,"(' Total number of transitions  = ',i8)") Ntransit
    end if
    !
    if (iverbose >= 0) then
       write(out,"(/' Statistical weights gns = ',4f12.1)") intensity%gns(1:)
    end if
    !
    ! In order to speed up the line strength evaluation, 
    ! S_if = | <i|mu|f> |^2  = | \sum_{nm} C_in C_fm <n|mu|m> |^2
    ! in three steps:
    ! 1. Evaluating the expansion of the i-state 
    !    s_{im} =  \sum_{n} C_in  <n|mu|m> 
    ! 2. performing the expansion of the final state f:
    !    s_{if} =  \sum_{m} C_fm  s_{im}. 
    ! 3. Building S_{if}
    !    S_{if} = s_{if}^2
    !
    !  The temporaly object s_{im} will be referted to as 
    !  a half-linestrength "half_linestr"
    !
    !  The matrix where some of the eigenvectors will be stored
    !
    if (iverbose>=5) call MemoryReport
    !
    deallocate(vecF)
    !
    ngrid = grid%npoints
    !
    allocate(kinmat(ngrid,ngrid),stat=info)
    call ArrayStart('kinmat',info,size(kinmat),kind(kinmat))
    !
    call kinetic_energy_grid_points(ngrid,kinmat)
    !
    ! ---------------------------------
    ! Dipole moment matrix elements
    ! --------------------------------- 
    !
    itransit = 0
    !
    if (intensity%npoints==-1) then
      write(out,"('RWF error: cross sectons Npoints is undefined. Define it using NPOINTS keyword in INTENSITY')")
      stop 'RWF error: cross sectons Npoints is undefined'
    endif
    !
    dnu = (intensity%freq_window(2)-intensity%freq_window(1))/real(intensity%Npoints-1,rk)
    !
    allocate(crosssections(intensity%npoints),stat = info)
    call ArrayStart('crosssections',info,size(crosssections),kind(crosssections))
    crosssections = 0
    !
    intens_cm_molecule  = 8.0d-36 * pi**3/ (3.0_rk * planck * vellgt)
    !
    ! loop over initial states
    !
    do indF = 1, nJ
       !
       jF = jval(indF)
       !
       if (iverbose>=3) write(out,"('J = ',f5.1)") jF
       !
       ! Primitive basis set 
       !
       if (iverbose>=3) write(out,'("Define the quanta book-keeping")')
       !
       call define_quanta_bookkeeping_RWF(iverbose,jval(indF),Nestates,Nrefstates,Nlambdasigmas)
       !
       if (iverbose>=3) write(out,'("...done!")')
       !
       ! Now we combine together the vibrational and sigma-lambda basis functions (as product)
       ! and the corresponding quantum numbers to form our final contracted basis set as well as
       ! the numbering of the contratced basis functions using only one index i.
       !
       !
       ! this how many states we get in total after the product of the
       ! vibrational and sigma-lambda basis sets:
       Ntotal = ngrid*Nlambdasigmas
       !
       if (Ntotal==0) cycle
       !
       allocate(icontr(Ntotal),stat=info)
       !
       allocate(half_linestr(Ntotal),stat=info)
       !
       call ArrayStart('half_linestr',info,size(half_linestr),kind(half_linestr))
       !
       if (iverbose>=4) write(out,'(/"Contracted basis set:")')
       if (iverbose>=4) write(out,'("     i     jrot ilevel ivib state v     spin    sigma lambda   omega   Name")')
       !
       ! How many vibrational state below Nrefstates?
       !
       Nvib_l = 0 
       !
       do ilevel = 1, vibrational_totalroots
          !
          istate  = vibrational_quantum_number(ilevel)%istate
          !
          if ( istate<=Nrefstates )  Nvib_l = Nvib_l + 1
          !
       enddo
       !
       ! build the bookkeeping: the object icontr will store this informtion
       !
       i = 0
       do ilevel = 1,Nlambdasigmas
         !
         istate = quanta_RWF(ilevel)%istate
         sigma = quanta_RWF(ilevel)%sigma
         imulti = quanta_RWF(ilevel)%imulti
         ilambda = quanta_RWF(ilevel)%ilambda
         omega = quanta_RWF(ilevel)%omega
         spini = quanta_RWF(ilevel)%spin
         !
         do igrid =1,ngrid
           !
           ! link the states in vibrarional and spin-rot basis components for the VIB-igrid basis
           !
           if (poten(istate)%gridvalue(igrid)>job%vibenermax(istate)) cycle
           !
           i = i + 1
           !
           icontr(i) = quanta_RWF(ilevel)
           icontr(i)%ivib = igrid
           icontr(i)%v = Nvib_l + igrid
           icontr(i)%ilevel = ilevel
           !
           ! print the quantum numbers
           if (iverbose>=4.and.indF==0) then
               write(out,'(i6,1x,f8.1,1x,i4,1x,i4,1x,i4,1x,f8.1,1x,f8.1,1x,i4,1x,f8.1,3x,a)') &
                       i,jval,ilevel,igrid,istate,&
                       spini,sigma,ilambda,omega,trim(poten(istate)%name)
           endif
           !
         enddo
       enddo
       !
       nlevelsF = i
       !
       ! allocate the hamiltonian matrix and an array for the energies of this size Ntotal
       allocate(hmat(nlevelsF,nlevelsF),stat=info)
       call ArrayStart('hmat',info,size(hmat),kind(hmat))
       !
       if (intensity%N_RWF_order>1) then 
          allocate(hmat_n(nlevelsF,nlevelsF,2:intensity%N_RWF_order),stat=info)
         call ArrayStart('hmat_n',info,size(hmat_n),kind(hmat_n))
       endif 
       !
       if (iverbose>=4) call MemoryReport
       !
       hmat = 0
       !
       !
       if (iverbose>=4) call TimerStart('Construct the hamiltonian')
       !
       if (iverbose>=3) write(out,'(/"Construct the hamiltonian matrix")')
       !
       !omp parallel do private(i,ivib,ilevel,istate,sigmai,imulti,ilambda,omegai,spini,v_i,jvib,jlevel,jstate,sigmaj, & 
       !                        jmulti,jlambda,omegaj,spinj,v_j,f_rot,erot,iL2,field,f_l2,f_s,f_t,iso,ibraket,ipermute,&
       !                        istate_,ilambda_,sigmai_,spini_,jstate_,jlambda_,sigmaj_,spinj_,isigmav,omegai_,       &
       !                        omegaj_,itau,ilxly,f_grid,f_l,f_ss) shared(hmat) schedule(guided)
       do i = 1,nlevelsF
         !
         ivib = icontr(i)%ivib ! igrid
         ilevel = icontr(i)%ilevel
         !
         istate = icontr(i)%istate
         sigmai = icontr(i)%sigma
         imulti = icontr(i)%imulti
         ilambda = icontr(i)%ilambda
         omegai = icontr(i)%omega
         spini = icontr(i)%spin
         !
         ! the diagonal contribution is the energy from the contracted vibrational solution
         !
         hmat(i,i) = poten(istate)%gridvalue(ivib)
         !
         f_rot=brot(1)%gridvalue(ivib)/sc
         erot = f_rot*( jF*(jF+1.0_rk) - omegai**2 -job%diag_L2_fact*real(ilambda**2,rk)  & 
                +   spini*(spini+1.0_rk) - sigmai**2 )
         !
         ! add the diagonal matrix element to the local spin-rotational matrix hmat
         hmat(i,i) = hmat(i,i) + erot
         !
         do j =i,nlevelsF
            !
            jvib = icontr(j)%ivib
            jlevel = icontr(j)%ilevel
            jstate = icontr(j)%istate
            sigmaj = icontr(j)%sigma
            jmulti = icontr(j)%imulti
            jlambda = icontr(j)%ilambda
            omegaj = icontr(j)%omega
            spinj = icontr(j)%spin
            !
            if (iverbose>=6) write(out,'("ilevel,ivib = ",2(i0,2x) )') ilevel,ivib
            !
            ! For the raw non-integrated basis add the kinetic energy matrix elemenent which otherwsie is not included in the 
            ! corresponding contracted enegy values 
            !
            if (istate==jstate.and.ilevel==jlevel) then
              !
              select case (trim(poten(istate)%integration_method))
              !
              case ('NONE','RAW')
                hmat(i,j) = hmat(i,j) + kinmat(ivib,jvib)/sc
              end select 
              !
            endif
            !
            hmat(j,i) =  hmat(i,j)
            !
         enddo  ! j
       enddo  ! i
       !omp end parallel do
       !
       if (iverbose>=3) write(out,'("...done!")')
       !
       if (iverbose>=4) call TimerStop('Construct the hamiltonian')
       !
       ! Transformation to the symmetrized basis set
       !
       ! |v,Lambda,Sigma,J,Omega,tau> = 1/sqrt(2) [ |v,Lambda,Sigma,J,Omega>+(-1)^tau |v,-Lambda,-Sigma,J,-Omega> ]
       allocate(Nirr(Ntotal,2),stat=info)
       call ArrayStart('Nirr',info,size(Nirr),kind(Nirr))
       !
       allocate(ilevel2i(Ntotal,2),stat=info)
       call ArrayStart('ilevel2i',info,size(ilevel2i),kind(ilevel2i))
       !
       !call transform_hmat_to_symmety_addapted_matrices(Ntotal,icontr,hmat,Nirr,Nsym,ilevel2i,transform)
       !
       ! transformaion to the symmetrised basis 
       !
       do indI = 1, nJ
         !
         jI = jval(indI)
         !
         do igammaI=1,Nrepresen
           !
           nlevelsI = eigen(indI,igammaI)%Nlevels
           dimenI = eigen(indI,igammaI)%Ndimen
           !
           if (nlevelsI==0) cycle
           !
           !parity_gu = poten(istateI)%parity%gu
           !isymI = correlate_to_Cs(igammaI,parity_gu)
           !
           if (abs(nint(jI-jF))>1.or.abs(nint(jI+jF))==0) cycle 
           !
           !do igammaF=1,Nrepresen
           !
           !nlevelsF =Nsym(igammaF)
           dimenF = Ntotal
           !
           !igammaF = igamma_pair(igammaI)
           !
           !parity_gu = poten(istateF)%parity%gu
           !isymF = correlate_to_Cs(igammaF,parity_gu)
           !
           !if (trim(job%symmetry) == "CS(M)".and.igammaF==igammaI) cycle
           !if (isymF /= igamma_pair(isymI)) cycle
           !
           allocate(dipole_mat(dimenF,nlevelsI),stat = info)
           call ArrayStart('dipole_mat',info,size(dipole_mat),kind(dipole_mat))
           dipole_mat = 0
           !
           Ilevels_loop : do ilevelI = 1, nlevelsI
             !
             !energy and and quanta of the final state
             !
             energyI = eigen(indI,igammaI)%val(ilevelI)
             !
             istateI  = eigen(indI,igammaI)%quanta(ilevelI)%istate
             !
             !dimension of the bases for the initial states
             !
             !energy, quanta, and gedeneracy order of the initial state
             quantaI => eigen(indI,igammaI)%quanta(ilevelI)
             istateI  = quantaI%istate
             ivibI    = quantaI%ivib
             ivI      = quantaI%v
             sigmaI   = quantaI%sigma
             spinI    = quantaI%spin
             ilambdaI = quantaI%ilambda
             omegaI   = quantaI%omega
             !
             ! reconstruct the symmetry for the C2v case which is different from Cs
             parity_gu = poten(istateI)%parity%gu
             isymI = correlate_to_Cs(igammaI,parity_gu)
             !
             call energy_filter_lower(jI,energyI,passed)
             !
             if (.not.passed) cycle
             !
             vecI(1:dimenI) = eigen(indI,igammaI)%vect(1:dimenI,ilevelI)
             !
             ! Compute the half-linestrength
             !
             half_linestr = 0
             !
             !if (isymF /= igamma_pair(isymI)) cycle
             !
             if (( intensity%J(1)+intensity%J(2)>0 )&
                 .and. abs(nint(jI-jF))<=1.and.nint(jI+jF)>=1) then 
                !
                call do_1st_half_linestrength_DVR(jI,jF,indI,indF,dimenI,dimenF,&
                                              vecI(1:dimenI),icontr,&
                                              half_linestr)
                !
             endif
             !
             !loop over final states
             !
             !
             !omp do private(ilevelF,energyF,dimenF,quantaF,istateF,ivibF,ivF,sigmaF,spinF,ilambdaF,omegaF,passed,&
             !omp& parity_gu,isymF,branch,nu_if,linestr,linestr2,A_einst,boltz_fc,absorption_int,tm) schedule(static) &
             !omp                                                                             & reduction(+:itransit)
             Flevels_loop: do ilevelF = 1,nlevelsF
                !
                !j = ilevel2i(ilevelF,isymF)
                !istateF = icontr(j)%istate
                !
                !parity_gu = poten(istateF)%parity%gu
                !isymF = correlate_to_Cs(igammaF,parity_gu)
                !
                !call TimerStart('Intens_Filter-3')
                !
                !call intens_filter_sym(jI,jF,isymI,isymF,igamma_pair,passed)
                !
                !call TimerStop('Intens_Filter-3')
                !
                !if (.not.passed) cycle Flevels_loop
                !
                !dipole_mat(ilevelF,ilevelI) = ddot(half_linestr,vibrational_contrfunc(:,ilevelF))
                !
                dipole_mat(ilevelF,ilevelI) = half_linestr(ilevelF)
                !
             end do Flevels_loop
             !omp enddo
             !
             if (iverbose>=5) call TimerReport
             !
           enddo Ilevels_loop
           !
           if (intensity%N_RWF_order>1) then 
             !
             call dgemm('N','N',nlevelsF,nlevelsF,nlevelsF,dalpha,hmat,nlevelsF,hmat,nlevelsF,dbeta,hmat_n(:,:,2),nlevelsF)
             !
           endif
           !
           if (intensity%N_RWF_order>3) then 
             !
             call dgemm('N','N',nlevelsF,nlevelsF,nlevelsF,dalpha,hmat_n(:,:,2),nlevelsF,hmat,nlevelsF,dbeta,hmat_n(:,:,3),nlevelsF)
             !
             call dgemm('N','N',nlevelsF,nlevelsF,nlevelsF,dalpha,hmat_n(:,:,3),nlevelsF,hmat,nlevelsF,dbeta,hmat_n(:,:,4),nlevelsF)
             !
           endif
           !
           !$omp parallel private(Amat,B,C,alloc_p) shared(crosssections) 
           allocate(Amat(nlevelsF,nlevelsF),B(nlevelsF),C(nlevelsF),stat = alloc_p)
           if (alloc_p/=0) then
               write (out,"(' RWF: ',i9,' trying to allocate arrays A, B, C')") alloc_p
               stop 'RWF A, B, C - out of memory'
           end if
           !call ArrayStart('RWF:Amat',info,size(Amat),kind(Amat))
           !call ArrayStart('RWF:Amat',info,size(B),kind(B))
           !
           ! Gaussuan parameter 
           !
           delta = log(2.0_rk)/intensity%gamma**2*0.5_rk
           !
           ! Wavenumber grid 
           !
           !$omp do private(inu,nu,ilevelI,istateI,parity_gu,isymI,energyI,ilevelF,ilevelR,RWF2,boltz_fc) schedule(static) 
           do inu = 1,intensity%npoints
              !
              nu = intensity%freq_window(1)+dnu*real(inu,rk)
              !
              if (iverbose>=4.and.mod(inu,100)==0) write(out,"(2f8.1,1x,i2,1x,i7,1x,'nu = ',f9.2)") Jf,Ji,igammaI,inu,nu
              !
              do ilevelI = 1, nlevelsI
                !
                istateI  = eigen(indI,igammaI)%quanta(ilevelI)%istate
                !
                ! reconstruct the symmetry for the C2v case which is different from Cs
                parity_gu = poten(istateI)%parity%gu
                isymI = correlate_to_Cs(igammaI,parity_gu)
                !
                energyI = eigen(indI,igammaI)%val(ilevelI)
                !
                Amat = 0 
                B = 0
                C = 0
                !
                do ilevelF = 1, nlevelsF
                  !
                  B(ilevelF) = dipole_mat(ilevelF,ilevelI)
                  !
                  do ilevelR = 1,nlevelsF
                    !
                    !
                    select case(trim(intensity%RWF_type)) 
                      !
                    case default
                      !
                      stop 'unknown RWF_type'
                      !
                    case ('GAUSSIAN')
                       !
                       Amat(ilevelF,ilevelR) = delta*( -hmat_n(ilevelF,ilevelR,2)+2.0_rk*hmat(ilevelF,ilevelR)*(nu + energyI) )
                       !
                       if (ilevelF==ilevelR) then
                         !
                         Amat(ilevelF,ilevelR) = Amat(ilevelF,ilevelR) + 1.0_rk  &
                         - delta*(nu + energyI)**2 
                         !
                       endif
                       !
                       if (intensity%N_RWF_order>3) then 
                          !
                          Amat(ilevelF,ilevelR) = Amat(ilevelF,ilevelR) + delta**2* &
                                           ( 0.5_rk*hmat_n(ilevelF,ilevelR,4) - 2.0_rk*hmat_n(ilevelF,ilevelR,3)*(nu + energyI)+&
                                            3.0_rk*hmat_n(ilevelF,ilevelR,2)*(nu + energyI)**2-&
                                            2.0_rk*hmat(ilevelF,ilevelR)*(nu + energyI)**3 )
                          !
                          if (ilevelF==ilevelR) then
                            !
                            Amat(ilevelF,ilevelR) = Amat(ilevelF,ilevelR) + delta**2*0.5_rk*(nu + energyI)**4 
                            !
                          endif
                          !
                       endif
                       !
                    case ('LORENTZIAN')
                       !
                       !Amat(ilevelF,ilevelR) = -transform(isymI)%matrix(ilevelF,ilevelR)/intensity%gamma**2
                       Amat(ilevelF,ilevelR) = hmat(ilevelF,ilevelR)*cmplx(0.0_rk,-1.0_rk/intensity%gamma**2)
                       !
                       if (ilevelF==ilevelR) then
                         !
                         Amat(ilevelF,ilevelR) = Amat(ilevelF,ilevelR) + (nu + energyI)*cmplx(0.0_rk,1.0_rk/intensity%gamma**2) &
                          +1.0_rk/intensity%gamma
                         !
                       endif
                       !
                       if (intensity%N_RWF_order>1) then 
                          !
                          Amat(ilevelF,ilevelR) = Amat(ilevelF,ilevelR) - &
                                                hmat_n(ilevelF,ilevelR,2)/intensity%gamma**3&
                                               +hmat(ilevelF,ilevelR)*2.0_rk*(nu + energyI)/intensity%gamma**3
                          !
                          if (ilevelF==ilevelR) then
                            !
                            Amat(ilevelF,ilevelR) = Amat(ilevelF,ilevelR) - &
                                                   (nu + energyI)**2/intensity%gamma**3
                            !
                          endif
                          !
                       endif
                       !
                    end select 
                    !
                  enddo
                  !
                enddo
                !
                !call lapack_gelss(Amat,b)
                !call  lapack_zgesv(Amat,b)
                !
                !C = matmul(Amat,B)
                !
                call zgemv('N',nlevelsF,nlevelsF,alpha_,Amat,nlevelsF,B,1,beta_,C,1)
                !
                !RWF2 = sum(conjg(C)*C)
                !
                RWF2 = zdotc(nlevelsF,C,1,C,1)
                !
                boltz_fc = intensity%gns(isymI)*real( (2*jI + 1)*(2 * jF + 1),rk )*nu *&
                     exp(-(energyI-intensity%ZPE) * beta) * (1.0_rk - exp(-nu * beta))/intensity%part_func
                !
                crosssections(inu) = crosssections(inu) + intens_cm_molecule*boltz_fc*RWF2
                !
              enddo
              !
           enddo
           !$omp enddo
           !
           deallocate(Amat,B,C)
           !$omp end parallel
           !call Arraystop('RWF:Amat')
           !
           deallocate(dipole_mat,stat = info)
           call ArrayStop('dipole_mat')
           !
         enddo
         !
       enddo
       !
       deallocate(hmat)
       call ArrayStop('hmat')
       !
       if (allocated(hmat_n)) then 
          deallocate(hmat_n)
          call ArrayStop('hmat_n')
       endif 
       !
       if (allocated(icontr))  deallocate(icontr)
       !
       !if (associated(transform(1)%matrix)) then 
       !   deallocate(transform(1)%matrix)
       !   call ArrayStop('transform')
       !endif
       !if (associated(transform(2)%matrix)) deallocate(transform(2)%matrix)
       !if (associated(transform(1)%irec)) deallocate(transform(1)%irec)
       !if (associated(transform(2)%irec)) deallocate(transform(2)%irec)
       !
       if (allocated(Nirr)) then 
          deallocate(Nirr)
          call ArrayStop('Nirr')
       endif
       !
       if (allocated(ilevel2i)) then 
          deallocate(ilevel2i)
          call ArrayStop('ilevel2i')
       endif
       !
       deallocate(half_linestr)
       call ArrayStop('half_linestr')
       !
    enddo
    !
    do inu = 1,intensity%npoints
       !
       nu = intensity%freq_window(1)+dnu*real(inu,rk)
       !
       write(transunit,"(f12.5,2x,e12.5)") nu,crosssections(inu)
       !
    enddo
    !
    deallocate(crosssections)
    call ArrayStop('crosssections')
    !
    deallocate(vecI)
    call ArrayStop('intensity-vectors')
    !
    deallocate(kinmat)
    call ArrayStop('kinmat')
    !
    if (trim(intensity%linelist_file)/="NONE") close(transunit,status="keep")
    !
    call TimerStop('Intensity calculations')
    !
  end subroutine rwf_dvr_intensity


  subroutine transform_hmat_to_symmety_addapted_matrices(Ntotal,icontr,hmat,Nirr,Nsym,ilevel2i,transform)
  
    implicit none
    integer(ik),intent(in)  :: Ntotal
    type(quantaT),intent(in) :: icontr(Ntotal)
    real(rk) :: hmat(Ntotal,Ntotal)
    integer(ik),intent(out) :: Nirr(Ntotal,2),Nsym(2),ilevel2i(Ntotal,2)
    type(matrixT),intent(out) :: transform(2)
    
    integer(ik),allocatable :: iswap(:),ilevel2isym(:,:)
    real(rk),allocatable :: vec(:),tau(:),Utransform(:,:,:)
    real(rk) :: vecti(2,2),vectj(2,2),pmat(2,2),smat(2,2)

    integer(ik) :: ilambda,ilevel,irrep,istate,isym,itau,ivib,jlambda,jlevel,jrrep,jstate,&
                   jsym,jtau,jvib,nlevels,i,j,info
    real(rk) :: sigmai,sigmaj,Ji,omegai,omegaj,spini,spinj

  
       allocate(iswap(Ntotal),vec(Ntotal),tau(Ntotal),ilevel2isym(Ntotal,2),stat=info)
       call ArrayStart('iswap-vec',info,size(iswap),kind(iswap))
       call ArrayStart('iswap-vec',info,size(vec),kind(vec))
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
         ivib    = icontr(i)%ivib
         !
         if (iswap(i)/=0) cycle
         !
         if (ilambda==0.and.nint(sigmai)==0) then
           !
           Nlevels = Nlevels + 1
       
           itau = nint(spini+jI)
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
             jvib    = icontr(j)%ivib
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
               itau = -ilambda+nint(spini-sigmai)+nint(jI-omegai)
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
       allocate(transform(1)%matrix(max(1,Nsym(1)),max(1,Nsym(1))),stat=info)
       allocate(transform(2)%matrix(max(1,Nsym(2)),max(1,Nsym(2))),stat=info)
       allocate(transform(1)%irec( max( 1,Nsym(1) ) ),stat=info)
       allocate(transform(2)%irec( max( 1,Nsym(2) ) ),stat=info)
       !
       call ArrayStart('transform',info,size(transform(1)%matrix),kind(transform(1)%matrix))
       call ArrayStart('transform',info,size(transform(2)%matrix),kind(transform(2)%matrix))
       call ArrayStart('transform',info,size(transform(1)%irec),kind(transform(1)%irec))
       call ArrayStart('transform',info,size(transform(2)%irec),kind(transform(2)%irec))
       !
       allocate(Utransform(Nlevels,2,2),stat=info)
       call ArrayStart('Utransform',info,size(Utransform),kind(Utransform))
       !
       ! Building the transformation to the symmetrized representaion 
       !
       do ilevel = 1,Nlevels
         !
         vecti = 0
         !
         if (any(Nirr(ilevel,:)==0)) then
           !
           do irrep = 1,sym%NrepresCs
             do itau = 1,Nirr(ilevel,irrep)
               !
               i = ilevel2i(ilevel,irrep)
               istate = icontr(i)%istate
               vecti(irrep,irrep) = 1.0_rk
               isym = ilevel2isym(ilevel,irrep)
               transform(irrep)%irec(isym) = ilevel
               !
             enddo
           enddo
           !
         else
           !
           vecti(1,1) = sqrt(0.5_rk)
           vecti(2,1) = sqrt(0.5_rk)*tau(ilevel)
           vecti(1,2) = sqrt(0.5_rk)
           vecti(2,2) =-sqrt(0.5_rk)*tau(ilevel)
           !
           do irrep = 1,sym%NrepresCs
              isym = ilevel2isym(ilevel,irrep)
              transform(irrep)%irec(isym) = ilevel
           enddo
           !
         endif
         !
         Utransform(ilevel,:,:) = vecti(:,:)
         !
         do jlevel = 1,Nlevels
           !
           vectj = 0
           !
           if (any(Nirr(jlevel,:)==0)) then
              !
              do jrrep = 1,sym%NrepresCs
                do jtau = 1,Nirr(jlevel,jrrep)
                  vectj(jrrep,jrrep) = 1.0_rk
                enddo
              enddo
              !
           else
              !
              vectj(1,1) = sqrt(0.5_rk)
              vectj(2,1) = sqrt(0.5_rk)*tau(jlevel)
              vectj(1,2) = sqrt(0.5_rk)
              vectj(2,2) =-sqrt(0.5_rk)*tau(jlevel)
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
           smat = matmul(transpose(vecti),matmul(pmat,vectj))
           !
           do irrep = 1,sym%NrepresCs
              do itau = 1,Nirr(ilevel,irrep)
                 !
                 isym = ilevel2isym(ilevel,irrep)
                 !
                 do jrrep = 1,sym%NrepresCs
                    do jtau = 1,Nirr(jlevel,jrrep)
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
       deallocate(iswap,vec)
       call ArrayStop('iswap-vec')
       !
       deallocate(tau,ilevel2isym)
       !
       deallocate(Utransform)
       call ArrayStop('Utransform')


  end subroutine transform_hmat_to_symmety_addapted_matrices
  


  !
  !
  function PQR_branch(jI,jF) result (X)
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
        elseif( nint(jI-jF)/=0 ) then
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
      do igammaI = 1,sym%Nrepresen
        !
        ! count number of hits
        !
        ngamma = 0
        igamma_pair(igammaI) = igammaI
        !
        do igammaF = 1,sym%Nrepresen
          !
          if (igammaI/=igammaF.and.intensity%isym_pairs(igammaI)==intensity%isym_pairs(igammaF)) then 
            !
            igamma_pair(igammaI) = igammaF
            !
            ngamma = ngamma + 1 
            !
            if (ngamma>1) then 
              !
              write(out,"('dm_intensity: Assumption that selection rules come in pairs is wrong!')")
              stop 'dm_intensity: Assumption that all selection rules work in pairs is wrong!'
              !
            endif   
            !
          endif
          !
        enddo
        !
        if ( nint(intensity%gns(igammaI)-intensity%gns(igamma_pair(igammaI)))/=0 ) then 
          !
          write(out,"('dm_intensity: selection rules do not agree with Gns')")
          stop 'dm_intensity: selection rules do not agree with Gns!'
          !
        endif   
        !
      enddo 
      !
     end subroutine find_igamma_pair


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

     end subroutine energy_filter_lower



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
             ! nuclear stat.weight: 
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

     end subroutine energy_filter_upper



     subroutine intens_filter(jI,jF,energyI,energyF,isymI,isymF,igamma_pair,passed)
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
              nu_if>small_.and.                                           &
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
          if (trim(intensity%action)=='ABSORPTION'.or.trim(intensity%action)=='EMISSION') then 
             !
             ! In order to avoid double counting of transitions
             ! we exclude jI=jF==intensity%J(2), i.e. Q branch for the highest J is never considered:
             !
             passed = passed.and.                                              &
             !
             (nint(jF-intensity%J(1))/=0.or.nint(jI-intensity%J(1))/=0.or.nint(jI+jF)==1).and.  &
             !
             !( ( nint(jF-intensity%J(1))/=0.or.nint(jI-intensity%J(1))/=0 ).and.intensity%J(1)>0 ).and.   &
             ( intensity%J(1)+intensity%J(2)>0 ).and. &
             !
             ! selection rules: 
             !
             intensity%isym_pairs(isymI)==intensity%isym_pairs(isymF).and.  &
             !
             igamma_pair(isymI)==isymF.and.                                 &
             !
             ! selection rules from the 3j-symbols
             !
             abs(nint(jI-jF))<=1.and.nint(jI+jF)>=1
             !
          endif
          !
     end subroutine intens_filter




     subroutine intens_filter_sym(jI,jF,isymI,isymF,igamma_pair,passed)
        !
        real(rk),intent(in) :: jI,jF
        integer(ik),intent(in) :: isymI,isymF
        integer(ik),intent(in) :: igamma_pair(sym%Nrepresen)
        logical,intent(out)    :: passed

          passed = .false.
          !
          if (                                                             &
              ! nuclear stat.weight: 
              !
              intensity%gns(isymI)>small_.and.                             &
              !
              jI>=intensity%J(1).and.                                      &
              jI<=intensity%J(2).and.                                      &
              !
              jF>=intensity%J(1).and.                                      &
              jF<=intensity%J(2)      ) then
              !
              passed = .true.
              !
          endif 
          !
          if (trim(intensity%action)=='ABSORPTION'.or.trim(intensity%action)=='EMISSION') then 
             !
             ! In order to avoid double counting of transitions
             ! we exclude jI=jF==intensity%J(2), i.e. Q branch for the highest J is never considered:
             !
             passed = passed.and.                                              &
             !
             (nint(jF-intensity%J(1))/=0.or.nint(jI-intensity%J(1))/=0.or.nint(jI+jF)==1).and.  &
             !
             !( ( nint(jF-intensity%J(1))/=0.or.nint(jI-intensity%J(1))/=0 ).and.intensity%J(1)>0 ).and.   &
             ( intensity%J(1)+intensity%J(2)>0 ).and. &
             !
             ! selection rules: 
             !
             intensity%isym_pairs(isymI)==intensity%isym_pairs(isymF).and.  &
             !
             igamma_pair(isymI)==isymF.and.                                 &
             !
             ! selection rules from the 3j-symbols
             !
             abs(nint(jI-jF))<=1.and.nint(jI+jF)>=1
             !
          endif
          !
     end subroutine intens_filter_sym


     subroutine matelem_filter(jI,jF,energyI,energyF,isymI,isymF,igamma_pair,passed)
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
          if (trim(intensity%action)=='ABSORPTION'.or.trim(intensity%action)=='EMISSION') then 
             !
             ! In order to avoid double counting of transitions
             ! we exclude jI=jF==intensity%J(2), i.e. Q branch for the highest J is never considered:
             !
             passed = passed.and.                                              &
             !
             (nint(jF-intensity%J(1))/=0.or.nint(jI-intensity%J(1))/=0.or.nint(jI+jF)==1).and.                    &
             !
             !( ( nint(jF-intensity%J(1))/=0.or.nint(jI-intensity%J(1))/=0 ).and.intensity%J(1)>0 ).and.   &
             ( intensity%J(1)+intensity%J(2)>0 ).and. &
             !
             ! selection rules: 
             !
             intensity%isym_pairs(isymI)==intensity%isym_pairs(isymF).and.  &
             !
             igamma_pair(isymI)==isymF.and.                                 &
             !
             ! selection rules from the 3j-symbols
             !
             abs(nint(jI-jF))<=1.and.nint(jI+jF)>=1
             !
          endif
          !
     end subroutine matelem_filter




     subroutine Jgamma_filter(jI,jF,isymI,isymF,igamma_pair,passed)
        !
        real(rk),intent(in) :: jI,jF
        integer(ik),intent(in) :: isymI,isymF
        integer(ik),intent(in) :: igamma_pair(sym%Nrepresen)
        logical,intent(out)    :: passed
          !
          passed = .true.
          !
          if (trim(intensity%action)=='ABSORPTION'.or.trim(intensity%action)=='EMISSION') then 
             !
             ! In order to avoid double counting of transitions
             ! we exclude jI=jF==intensity%J(2), i.e. Q branch for the highest J is never considered:
             !
             passed = passed.and.                                              &
             !
             !( ( nint(jF-intensity%J(1))/=0.or.nint(jI-intensity%J(1))/=0 ).and.intensity%J(1)+intensity%J(2)>0 ).and. &
             ( intensity%J(1)+intensity%J(2)>0 ).and. &
             !
             ! selection rules: 
             !
             intensity%isym_pairs(isymI)==intensity%isym_pairs(isymF).and.  &
             !
             igamma_pair(isymI)==isymF.and.                                 &
             !
             ! selection rules from the 3j-symbols
             !
             abs(nint(jI-jF))<=1.and.nint(jI+jF)>=1
             !
          endif
          !
     end subroutine Jgamma_filter





 ! 
 !read eigenvalues and their assignment 
 !
 subroutine Sort_levels(iverbose,njval, jval)
    integer(ik), intent(in) :: iverbose,njval
    real(rk), intent(in) :: jval(njval)

    integer(ik)             :: jind, nroots, nlevels,  iroot, ilevel, jlevel, &
                               info,itau,jtau

    real(rk)                :: energy

    !
    logical                 :: passed
    integer(ik)             :: iind,jroot,igamma
    !
    if (iverbose>=2) write(out,"(/'Read and sort eigenvalues in increasing order...')")
    !
    call TimerStart('Sort eigenvalues')

    if (.not. allocated(basis)) stop 'Sort_levels error: associated(basis) = .false.'
    !
    ! In practice we do not need all stored rootes, but only those that pass the 
    ! filters given in the input. Here we perform a pre-selection of the states, 
    ! which are actually required in the calculations. Towards this we read the unit twice:
    ! a) to count and b) to read.  
    ! After reading all states will  be combined together and sorted wrt energy increasing 
    ! to produce just one list of states. 
    !
    ! Astimate the maximal number of energy records
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
      write(out,'("The number of selected levels ",i0," does not agree with Neigenlevels =  ",i0)') ilevel,Neigenlevels
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
      ! In order to speed up the line strength evaluation, 
      ! S_if = | <i|mu|f> |^2  = | \sum_{nm} C_in C_fm <n|mu|m> |^2
      ! in three steps:
      ! 1. Evaluating the expansion of the i-state 
      !    s_{im} =  \sum_{n} C_in  <n|mu|m> 
      ! 2. performing the expansion of the final state f:
      !    s_{if} =  \sum_{m} C_fm  s_{im}. 
      ! 3. Building S_{if}
      !    S_{if} = s_{if}^2
      ! 
      !  This routine performs the first half of <i|mu|f>. 
      !

      subroutine do_1st_half_linestrength(jI,jF,indI,indF,dimenI,dimenF,vector,half_ls)

        real(rk),intent(in)     :: jI,jF
        integer(ik),intent(in)  :: indI,indF,dimenI,dimenF
        real(rk),intent(in)     :: vector(:)
        real(rk),intent(out)    :: half_ls(:)
        integer(ik)             :: icontrF,icontrI, & 
                                   ivibF,ivibI,idip,istateI,istateF,ilambdaF,ilambdaI
        integer(ik)             :: ipermute,istateI_,ilambdaI_,ilambdaF_,isigmav,iomegaI_,istateF_,itau,iomegaF_
        real(rk)                :: ls, f3j, omegaI,omegaF,sigmaF,sigmaI,spinF,spinI
        real(rk)                :: spinI_,spinF_,f_t
        type(fieldT),pointer    :: field
          !
          !dms_tmp = dipole_me
          !
          call TimerStart('do_1st_half_linestr')
          !
          half_ls    = 0
          !
          !loop over final state basis components
          !
          !omp parallel do private(irootF,icontrF,ktau,kF,tauF,cirootI,irootI,icontrI,tauI,sigmaI,sigmaF,kI, & 
          !                   &    irow,icol,cind,f3j,ls) shared(half_ls) schedule(guided)
          loop_F : do icontrF = 1, dimenF
               !
               ivibF = basis(indF)%icontr(icontrF)%ivib
               istateF = basis(indF)%icontr(icontrF)%istate
               omegaF = basis(indF)%icontr(icontrF)%omega
               sigmaF = basis(indF)%icontr(icontrF)%sigma
               spinF = basis(indF)%icontr(icontrF)%spin
               ilambdaF = basis(indF)%icontr(icontrF)%ilambda
               !
               iomegaF_ = nint(omegaF)
               if (mod(nint(2.0_rk*omegaF+1.0_rk),2)==0 ) iomegaF_ = nint((2.0_rk*omegaF-1.0_rk)*0.5_rk)
               !
               loop_I : do icontrI = 1, dimenI
                  !
                  ivibI   = basis(indI)%icontr(icontrI)%ivib
                  istateI = basis(indI)%icontr(icontrI)%istate
                  omegaI  = basis(indI)%icontr(icontrI)%omega
                  sigmaI  = basis(indI)%icontr(icontrI)%sigma
                  spinI   = basis(indI)%icontr(icontrI)%spin
                  ilambdaI= basis(indI)%icontr(icontrI)%ilambda
                  !
                  if (abs(nint(omegaF - omegaI))>1.or.nint(spinI-spinF)/=0.or.nint(sigmaI-sigmaF)/=0) cycle loop_I
                  if (abs(nint(omegaF - omegaI))==0.and.ilambdaI/=ilambdaF) cycle loop_I
                  if (abs(nint(omegaF - omegaI))==1.and.abs(ilambdaI-ilambdaF)/=1) cycle loop_I
                  !
                  iomegaI_ = int(omegaI)
                  if (mod(nint(2.0_rk*omegaI+1.0_rk),2)==0 ) iomegaI_ = nint((2.0_rk*omegaI-1.0_rk)*0.5_rk)
                  !
                  f3j = three_j(jI, 1.0_rk, jF, omegaI, omegaF - omegaI, -omegaF)
                  !f3j = three_j0(jI, 1.0_rk, jF, omegaI, omegaF - omegaI, -omegaF)
                  ! 
                  ! 3j-symbol selection rule
                  !
                  if (abs(f3j)<intensity%threshold%coeff) cycle loop_I
                  !
                  !index of the corresponding vibrational contracted matrix element (cind)
                  !compute line strength
                  !
                  ls = 0 
                  !
                  loop_idipole : do idip =1,Ndipoles
                    !
                    field => dipoletm(idip)
                    !
                    do ipermute  = 0,1
                      !
                      if (ipermute==0) then
                        !
                        istateI_ = field%istate ; ilambdaI_ = field%lambda  ; spinI_ = field%spini
                        istateF_ = field%jstate ; ilambdaF_ = field%lambdaj ; spinF_ = field%spinj
                        !
                      else  ! permute
                        !
                        istateF_ = field%istate ; ilambdaF_ = field%lambda  ; spinF_ = field%spini
                        istateI_ = field%jstate ; ilambdaI_ = field%lambdaj ; spinI_ = field%spinj
                        !
                      endif
                      !
                      ! however the permutation makes sense only when for non diagonal <State,Lambda,Spin|F|State',Lambda',Spin'>
                      ! otherwise it will cause a double counting:
                      !
                      if (ipermute==1.and.istateI_==istateF_.and.ilambdaI_==ilambdaF_.and.nint(spinI_-spinF_)==0) cycle
                      !
                      ! check if we at the right electronic states
                      if( istateI/=istateI_.or.istateF/=istateF_ ) cycle
                      !
                      ! We should also take into account that Lambda can change sign (only Lambda>0 is given in input)
                      ! In order to recover other combinations we apply the symmetry transformation
                      ! laboratory fixed inversion which is equivalent to the sigmav operation 
                      !                    (sigmav= 0 correspond to the unitary transformation)
                      do isigmav = 0,1
                        !
                        ! the permutation is only needed if at least some of the quanta is not zero. 
                        ! otherwise it should be skipped to avoid the double counting.
                        if( isigmav==1.and. abs( field%lambda ) + abs( field%lambdaj )==0 ) cycle
                
                        ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                        ilambdaI_ = ilambdaI_*(-1)**isigmav
                        ilambdaF_ = ilambdaF_*(-1)**isigmav
                        !
                        ! proceed only if the quantum numbers of the field equal to the corresponding <i| and |j> quantum numbers:
                        if (ilambdaI_/=ilambdaI.or.ilambdaF_/=ilambdaF) cycle
                        !
                        ! check the selection rule Delta Lambda = +/1
                        if (abs(ilambdaI-ilambdaF)>1) cycle
                        !
                        ! double check
                        !if (spini/=poten(istate)%spini.or.spinj/=poten(jstate)%spini) then
                        !  write(out,'("dipole_intens: reconsrtucted spini ",f8.1," or spinj ",f8.1, & 
                        !            & " do not agree with stored values ",f8.1,x,f8.1)') &
                        !        spini,spinj,poten(istate)%spini,poten(jstate)%spini
                        !  stop 'dipole_intens: wrongly reconsrtucted spini or spinj'
                        !endif
                        !
                        !f_grid  = field%matelem(ivib,jvib)
                        !
                        f_t = field%matelem(ivibI,ivibF)
                        !
                        ! the result of the symmetry transformation:
                        if (isigmav==1) then
                          !
                          itau = 0
                          !
                          if (ilambdaI_==0.and.poten(istateI)%parity%pm==-1) itau = itau+1
                          if (ilambdaF_==0.and.poten(istateF)%parity%pm==-1) itau = itau+1
                          !
                          f_t = f_t*(-1.0_rk)**(itau)
                          !
                        endif
                        !
                        ls  =  f_t*f3j*vector(icontrI)
                        !
                        half_ls(icontrF) = half_ls(icontrF) + (-1.0_rk)**(iomegaI_)*ls
                        !
                      enddo
                      !
                    enddo
                    !
                  enddo loop_idipole
                  !
               end do  loop_I
               !
            end do   loop_F
            !omp end parallel do
            !
            call TimerStop('do_1st_half_linestr')
            !
      end subroutine do_1st_half_linestrength
      !


      subroutine do_1st_half_linestrength_DVR(jI,jF,indI,indF,dimenI,dimenF,vector,icontr,half_ls)

        real(rk),intent(in)      :: jI,jF
        integer(ik),intent(in)   :: indI,indF,dimenI,dimenF
        real(rk),intent(in)      :: vector(:)
        type(quantaT),intent(in) :: icontr(dimenF)
        real(rk),intent(out)     :: half_ls(:)
        integer(ik)              :: icontrF,icontrI, & 
                                    ivibF,ivibI,idip,istateI,istateF,ilambdaF,ilambdaI,vF
        integer(ik)              :: ipermute,istateI_,ilambdaI_,ilambdaF_,isigmav,iomegaI_,istateF_,itau,iomegaF_
        real(rk)                 :: ls, f3j, omegaI,omegaF,sigmaF,sigmaI,spinF,spinI
        real(rk)                 :: spinI_,spinF_,f_t
        type(fieldT),pointer     :: field
          !
          !dms_tmp = dipole_me
          !
          call TimerStart('do_1st_half_linestr')
          !
          half_ls    = 0
          !
          !loop over final state basis components
          !
          !omp parallel do private(irootF,icontrF,ktau,kF,tauF,cirootI,irootI,icontrI,tauI,sigmaI,sigmaF,kI, & 
          !                   &    irow,icol,cind,f3j,ls) shared(half_ls) schedule(guided)
          loop_F : do icontrF = 1, dimenF
               !
               ivibF = icontr(icontrF)%ivib
               istateF =icontr(icontrF)%istate
               omegaF =icontr(icontrF)%omega
               sigmaF =icontr(icontrF)%sigma
               spinF =icontr(icontrF)%spin
               ilambdaF =icontr(icontrF)%ilambda
               vF = icontr(icontrF)%v
               !
               ! shift to the original vibrational basis which we assume to start from the count of the lower (Nrefstate) state
               !
               iomegaF_ = nint(omegaF)
               if (mod(nint(2.0_rk*omegaF+1.0_rk),2)==0 ) iomegaF_ = nint((2.0_rk*omegaF-1.0_rk)*0.5_rk)
               !
               loop_I : do icontrI = 1, dimenI
                  !
                  ivibI   = basis(indI)%icontr(icontrI)%ivib
                  istateI = basis(indI)%icontr(icontrI)%istate
                  omegaI  = basis(indI)%icontr(icontrI)%omega
                  sigmaI  = basis(indI)%icontr(icontrI)%sigma
                  spinI   = basis(indI)%icontr(icontrI)%spin
                  ilambdaI= basis(indI)%icontr(icontrI)%ilambda
                  !
                  if (abs(nint(omegaF - omegaI))>1.or.nint(spinI-spinF)/=0.or.nint(sigmaI-sigmaF)/=0) cycle loop_I
                  if (abs(nint(omegaF - omegaI))==0.and.ilambdaI/=ilambdaF) cycle loop_I
                  if (abs(nint(omegaF - omegaI))==1.and.abs(ilambdaI-ilambdaF)/=1) cycle loop_I
                  !
                  iomegaI_ = int(omegaI)
                  if (mod(nint(2.0_rk*omegaI+1.0_rk),2)==0 ) iomegaI_ = nint((2.0_rk*omegaI-1.0_rk)*0.5_rk)
                  !
                  f3j = three_j(jI, 1.0_rk, jF, omegaI, omegaF - omegaI, -omegaF)
                  !f3j = three_j0(jI, 1.0_rk, jF, omegaI, omegaF - omegaI, -omegaF)
                  ! 
                  ! 3j-symbol selection rule
                  !
                  if (abs(f3j)<intensity%threshold%coeff) cycle loop_I
                  !
                  !index of the corresponding vibrational contracted matrix element (cind)
                  !compute line strength
                  !
                  ls = 0 
                  !
                  loop_idipole : do idip =1,Ndipoles
                    !
                    field => dipoletm(idip)
                    !
                    do ipermute  = 0,1
                      !
                      if (ipermute==0) then
                        !
                        istateI_ = field%istate ; ilambdaI_ = field%lambda  ; spinI_ = field%spini
                        istateF_ = field%jstate ; ilambdaF_ = field%lambdaj ; spinF_ = field%spinj
                        !
                      else  ! permute
                        !
                        istateF_ = field%istate ; ilambdaF_ = field%lambda  ; spinF_ = field%spini
                        istateI_ = field%jstate ; ilambdaI_ = field%lambdaj ; spinI_ = field%spinj
                        !
                      endif
                      !
                      ! however the permutation makes sense only when for non diagonal <State,Lambda,Spin|F|State',Lambda',Spin'>
                      ! otherwise it will cause a double counting:
                      !
                      if (ipermute==1.and.istateI_==istateF_.and.ilambdaI_==ilambdaF_.and.nint(spinI_-spinF_)==0) cycle
                      !
                      ! check if we at the right electronic states
                      if( istateI/=istateI_.or.istateF/=istateF_ ) cycle
                      !
                      ! We should also take into account that Lambda can change sign (only Lambda>0 is given in input)
                      ! In order to recover other combinations we apply the symmetry transformation
                      ! laboratory fixed inversion which is equivalent to the sigmav operation 
                      !                    (sigmav= 0 correspond to the unitary transformation)
                      do isigmav = 0,1
                        !
                        ! the permutation is only needed if at least some of the quanta is not zero. 
                        ! otherwise it should be skipped to avoid the double counting.
                        if( isigmav==1.and. abs( field%lambda ) + abs( field%lambdaj )==0 ) cycle
                
                        ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                        ilambdaI_ = ilambdaI_*(-1)**isigmav
                        ilambdaF_ = ilambdaF_*(-1)**isigmav
                        !
                        ! proceed only if the quantum numbers of the field equal to the corresponding <i| and |j> quantum numbers:
                        if (ilambdaI_/=ilambdaI.or.ilambdaF_/=ilambdaF) cycle
                        !
                        ! check the selection rule Delta Lambda = +/1
                        if (abs(ilambdaI-ilambdaF)>1) cycle
                        !
                        ! double check
                        !if (spini/=poten(istate)%spini.or.spinj/=poten(jstate)%spini) then
                        !  write(out,'("dipole_intens: reconsrtucted spini ",f8.1," or spinj ",f8.1, & 
                        !            & " do not agree with stored values ",f8.1,x,f8.1)') &
                        !        spini,spinj,poten(istate)%spini,poten(jstate)%spini
                        !  stop 'dipole_intens: wrongly reconsrtucted spini or spinj'
                        !endif
                        !
                        !f_grid  = field%matelem(ivib,jvib)
                        !
                        f_t = field%matelem(ivibI,vF)
                        !
                        ! the result of the symmetry transformation:
                        if (isigmav==1) then
                          !
                          itau = 0
                          !
                          if (ilambdaI_==0.and.poten(istateI)%parity%pm==-1) itau = itau+1
                          if (ilambdaF_==0.and.poten(istateF)%parity%pm==-1) itau = itau+1
                          !
                          f_t = f_t*(-1.0_rk)**(itau)
                          !
                        endif
                        !
                        ls  =  f_t*f3j*vector(icontrI)
                        !
                        half_ls(icontrF) = half_ls(icontrF) + (-1.0_rk)**(iomegaI_)*ls
                        !
                      enddo
                      !
                    enddo
                    !
                  enddo loop_idipole
                  !
               end do  loop_I
               !
            end do   loop_F
            !omp end parallel do
            !
            call TimerStop('do_1st_half_linestr')
            !
      end subroutine do_1st_half_linestrength_DVR

      !


      subroutine do_1st_half_tm(indI,indF,dimenI,dimenF,vector,half_tm)

        integer(ik),intent(in)  :: indI,indF,dimenI,dimenF
        real(rk),intent(in)     :: vector(:)
        real(rk),intent(out)    :: half_tm(:)
        integer(ik)             :: icontrF,icontrI, ivibF,ivibI,idip,istateI,istateF
        real(rk)                :: f
          !
          half_tm    = 0
          !
          !loop over final state basis components
          !
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
                  do idip = 1,Ndipoles
                    !
                    if (dipoletm(idip)%istate/=istateI.or.dipoletm(idip)%jstate/=istateF) cycle
                    f = dipoletm(idip)%matelem(ivibI,ivibF)
                    !
                    half_tm(icontrF) = half_tm(icontrF) + f*vector(icontrI)
                    !
                  enddo
                  !
               end do loop_I
               !
          end do  loop_F
          !
      end subroutine do_1st_half_tm



      subroutine  do_LF_matrix_elements(iLF,iunit,jI,jF,icount)
        implicit none 
        real(rk),intent(in)     :: jI,jF
        integer(ik),intent(in)  :: iLF,iunit
        integer(ik),intent(out),optional :: icount
        !real(rk),intent(out)    :: M(-jmax:jmax,-jmax:jmax,3)
        integer(ik)             :: imI_,imF_,icount_,isigma
        real(rk)                :: f3j, mI,mF,M_
        real(rk)                :: f_t
        real(rk)                :: t(3,-1:1)
          !
          t = 0
          !
          !t(1,-1) = cmplx(-1.0_rk/sqrt(2.0_rk),0.0_rk)
          !t(1,1)  = cmplx( 1.0_rk/sqrt(2.0_rk),0.0_rk)
          !t(2,-1) = cmplx(0.0_rk,1.0_rk/sqrt(2.0_rk))
          !t(2,1)  = cmplx(0.0_rk,1.0_rk/sqrt(2.0_rk))
          !t(3,0)  = cmplx(1.0_rk,0.0_rk)
          !
          ! Only the y-component is imaginary, which we make real here but indicate in richmol file
          !
          t(1,-1) =  1.0_rk/sqrt(2.0_rk)
          t(1,1)  = -1.0_rk/sqrt(2.0_rk)
          t(2,-1) = 1.0_rk/sqrt(2.0_rk)
          t(2,1)  = 1.0_rk/sqrt(2.0_rk)
          t(3,0)  = 1.0_rk
          !
          !dms_tmp = dipole_me
          !
          call TimerStart('do_LF_matrix_elements')
          !
          M_ = 0
          !
          icount_ = 0
          !
          loop_F : do imF_ = -nint(Jf), nint(Jf)
             !
             mF = imF_
             if (mod(nint(2.0_rk*jF+1.0_rk),2)==0 ) mF = imF_+0.5_rk
             !
             loop_I : do imI_ = -nint(Ji), nint(Ji)
                  !
                  mI = imI_
                  if (mod(nint(2.0_rk*jI+1.0_rk),2)==0 ) mI = imI_+0.5_rk
                  !
                  isigma = nint(mF - mI)
                  !
                  if (abs(isigma)>1) cycle
                  !
                  f3j = three_j(jI, 1.0_rk, jF, mI, mF - mI, -mF)
                  !
                  f_t = 0
                  !
                  f_t = T(iLF,isigma)
                  !
                  !M_t(imI_,imF_,iLF) = f_t
                  !
                  M_ = (-1.0_rk)**(imI_)*f_t*f3j*sqrt(( 2.0_rk*jI + 1.0_rk )*( 2.0_rk * jF + 1.0_rk ))
                  !
                  if (abs(M_)>small_) then 
                    icount_ = icount_ + 1
                    if (.not.present(icount)) write(iunit,"(2(i6),3x,e24.14)") imI_,imF_,M_ 
                  endif
                  !
               end do  loop_I
               !
          end do   loop_F
          !
          if (present(icount)) icount = icount_ 
          !omp end parallel do
          !
          call TimerStop('do_LF_matrix_elements')
          !
      end subroutine do_LF_matrix_elements


      subroutine do_matelem_pec(jI,jF,indI,indF,dimenI,dimenF,vector,half_me)

        implicit none
        real(rk),intent(in)     :: jI,jF
        integer(ik),intent(in)  :: indI,indF,dimenI,dimenF
        real(rk),intent(in)     :: vector(:)
        real(rk),intent(out)    :: half_me(:)
        integer(ik)             :: icontrF,icontrI, & 
                                   ivibF,ivibI,ipec,istateI,istateF,ilambdaF,ilambdaI
        integer(ik)             :: ipermute,istateI_,ilambdaI_,ilambdaF_,isigmav,iomegaI_,istateF_,itau,iomegaF_
        real(rk)                :: omegaI,omegaF,sigmaF,sigmaI,spinF,spinI
        real(rk)                :: spinI_,spinF_,f_t
        type(fieldT),pointer    :: field
          !
          !dms_tmp = dipole_me
          !
          call TimerStart('do_matelem_pec')
          !
          half_me    = 0
          !
          if (nint(jF-jI)/=0) then 
            call TimerStop('do_matelem_pec')
            return 
          endif
          !
          !loop over final state basis components
          !
          !omp parallel do private(irootF,icontrF,ktau,kF,tauF,cirootI,irootI,icontrI,tauI,sigmaI,sigmaF,kI, & 
          !                   &    irow,icol,cind,f3j,me) shared(half_me) schedule(guided)
          loop_F : do icontrF = 1, dimenF
               !
               ivibF = basis(indF)%icontr(icontrF)%ivib
               istateF = basis(indF)%icontr(icontrF)%istate
               omegaF = basis(indF)%icontr(icontrF)%omega
               sigmaF = basis(indF)%icontr(icontrF)%sigma
               spinF = basis(indF)%icontr(icontrF)%spin
               ilambdaF = basis(indF)%icontr(icontrF)%ilambda
               !
               iomegaF_ = nint(omegaF)
               if (mod(nint(2.0_rk*omegaF+1.0_rk),2)==0 ) iomegaF_ = nint((2.0_rk*omegaF-1.0_rk)*0.5_rk)
               !
               loop_I : do icontrI = 1, dimenI
                  !
                  ivibI   = basis(indI)%icontr(icontrI)%ivib
                  istateI = basis(indI)%icontr(icontrI)%istate
                  omegaI  = basis(indI)%icontr(icontrI)%omega
                  sigmaI  = basis(indI)%icontr(icontrI)%sigma
                  spinI   = basis(indI)%icontr(icontrI)%spin
                  ilambdaI= basis(indI)%icontr(icontrI)%ilambda
                  !
                  if (abs(nint(omegaF - omegaI))>0.or.nint(spinI-spinF)/=0.or.nint(sigmaI-sigmaF)/=0) cycle loop_I
                  if (abs(nint(omegaF - omegaI))==0.and.ilambdaI/=ilambdaF) cycle loop_I
                  !
                  iomegaI_ = int(omegaI)
                  if (mod(nint(2.0_rk*omegaI+1.0_rk),2)==0 ) iomegaI_ = nint((2.0_rk*omegaI-1.0_rk)*0.5_rk)
                  !
                  !index of the corresponding vibrational contracted matrix element (cind)
                  !compute mat. element
                  !
                  loop_ipec : do ipec =1,nQuadrupoles
                    !
                    field => quadrupoletm(ipec)
                    !
                    do ipermute  = 0,1
                      !
                      if (ipermute==0) then
                        !
                        istateI_ = field%istate ; ilambdaI_ = field%lambda  ; spinI_ = field%spini
                        istateF_ = field%jstate ; ilambdaF_ = field%lambdaj ; spinF_ = field%spinj
                        !
                      else  ! permute
                        !
                        istateF_ = field%istate ; ilambdaF_ = field%lambda  ; spinF_ = field%spini
                        istateI_ = field%jstate ; ilambdaI_ = field%lambdaj ; spinI_ = field%spinj
                        !
                      endif
                      !
                      ! however the permutation makes sense only when for non diagonal <State,Lambda,Spin|F|State',Lambda',Spin'>
                      ! otherwise it will cause a double counting:
                      !
                      if (ipermute==1.and.istateI_==istateF_.and.ilambdaI_==ilambdaF_.and.nint(spinI_-spinF_)==0) cycle
                      !
                      ! check if we at the right electronic states
                      if( istateI/=istateI_.or.istateF/=istateF_ ) cycle
                      !
                      ! We should also take into account that Lambda can change sign (only Lambda>0 is given in input)
                      ! In order to recover other combinations we apply the symmetry transformation
                      ! laboratory fixed inversion which is equivalent to the sigmav operation 
                      !                    (sigmav= 0 correspond to the unitary transformation)
                      do isigmav = 0,1
                        !
                        ! the permutation is only needed if at least some of the quanta is not zero. 
                        ! otherwise it should be skipped to avoid the double counting.
                        if( isigmav==1.and. abs( field%lambda ) + abs( field%lambdaj )==0 ) cycle
                
                        ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                        ilambdaI_ = ilambdaI_*(-1)**isigmav
                        ilambdaF_ = ilambdaF_*(-1)**isigmav
                        !
                        ! proceed only if the quantum numbers of the field equal to the corresponding <i| and |j> quantum numbers:
                        if (ilambdaI_/=ilambdaI.or.ilambdaF_/=ilambdaF) cycle
                        !
                        ! check the selection rule Delta Lambda = +/1
                        if (abs(ilambdaI-ilambdaF)>1) cycle
                        !
                        f_t = field%matelem(ivibI,ivibF)
                        !
                        ! the result of the symmetry transformation:
                        if (isigmav==1) then
                          !
                          itau = 0
                          !
                          if (ilambdaI_==0.and.poten(istateI)%parity%pm==-1) itau = itau+1
                          if (ilambdaF_==0.and.poten(istateF)%parity%pm==-1) itau = itau+1
                          !
                          f_t = f_t*(-1.0_rk)**(itau)
                          !
                        endif
                        !
                        f_t  =  f_t*vector(icontrI)
                        !
                        half_me(icontrF) = half_me(icontrF) + f_t
                        !
                      enddo
                      !
                    enddo
                    !
                  enddo loop_ipec
                  !
               end do  loop_I
               !
            end do   loop_F
            !omp end parallel do
            !
            call TimerStop('do_matelem_pec')
            !
      end subroutine do_matelem_pec


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
     if (job%isym_do(igamma).and.energy-job%ZPE>=erange(1).and.  &
         Jval>=intensity%J(1).and.Jval<=intensity%J(2).and.&
         energy-job%ZPE<=erange(2)) then 
         !
         passed = .true.
         !
     endif 

 end subroutine Energy_filter


  function three_j0(a,b,c,al,be,ga)

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
      !if(-ga.ne.al+be) return
      if(nint(ga+al+be).ne.0) return
!
!
!     compute delta(abc)
!
      delta=sqrt( fakt(a+b-c)*fakt(a+c-b)*fakt(b+c-a)/fakt(a+b+c+1.0_rk) )
!
!
      term1=fakt(a+al)*fakt(a-al)
      term2=fakt(b-be)*fakt(b+be)
      term3=fakt(c+ga)*fakt(c-ga)
      term=sqrt( (2.0_rk*c+1.0_rk)*term1*term2*term3 )
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
        dnew=real(new,rk)
        term4=fakt(a-al-dnew)*fakt(c-b+al+dnew)
        term5=fakt(b+be-dnew)*fakt(c-a-be+dnew)
        term6=fakt(dnew)*fakt(a+b-c-dnew)
        summ=summ+(-1.0_rk)**new/(term4*term5*term6)
      enddo
!
!     so clebsch-gordon <j1j2m1m2ljm> is clebsh
!
      clebsh=delta*term*summ/sqrt(10.0_rk)
!
!     convert clebsch-gordon to three_j
!
      iphase=nint(a-b-ga)
      minus = -1.0_rk
      if (mod(iphase,2).eq.0) minus = 1.0_rk
      three_j0=minus*clebsh/sqrt(2.0_rk*c+1.0_rk)

!     threej=(-1.d0)**(iphase)*clebsh/sqrt(2.0_rk*c+1.d0)
!
!
   end function three_j0



   function fakt(a) result (f)

      real(rk),intent(in) :: a
      real(rk)            :: ax,f
      integer(ik)         :: i,ic
!

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
      do  i=1,ic-1
        f=f*(ax-real(i,rk)*0.1_rk)
      enddo

    end function fakt



  !
  ! This basis set composition is for contr=VIB, State-Lambda-Sigma
  !
  subroutine define_quanta_bookkeeping_RWF(iverbose,jval,Nestates,Nrefstates,Nlambdasigmas)
    !
    integer(ik),intent(in) :: iverbose
    real(rk),intent(in)    :: jval
    integer(ik),intent(in) :: nestates,Nrefstates
    integer(ik),intent(out) :: Nlambdasigmas ! to count states with different lambda/sigma
    integer(ik) :: ilevel,itau,ilambda,nlevels,multi_max,imulti,istate,multi,alloc,taumax
    real(rk)    :: sigma,omega
    integer(ik) :: iquanta2ilevel

    !
    if (iverbose>=4) call TimerStart('Define quanta')
    !
    ilevel = 0
    multi_max = 1
    !
    ! count states
    !
    do istate = Nrefstates+1,nestates
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
    if (allocated(quanta_RWF)) then
       deallocate(quanta_RWF)
       !deallocate(iquanta2ilevel)
       !call ArrayStop('quanta_RWF')
    endif
    !
    allocate(quanta_RWF(nlevels),stat=alloc)
    !allocate(iquanta2ilevel(Nestates,0:1,multi_max),stat=alloc)
    !call ArrayStart('quanta',alloc,size(iquanta2ilevel),kind(iquanta2ilevel))
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
    do istate = Nrefstates+1,nestates
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
            quanta_RWF(ilevel)%spin = poten(istate)%spini
            quanta_RWF(ilevel)%istate = istate
            quanta_RWF(ilevel)%sigma = sigma
            quanta_RWF(ilevel)%imulti = multi
            quanta_RWF(ilevel)%ilambda = ilambda
            quanta_RWF(ilevel)%omega = real(ilambda,rk)+sigma
            !iquanta2ilevel(istate,itau,imulti) = ilevel
            !
            ! print out quanta
            !
            if (iverbose>=4) write(out,'(i6,1x,f8.1,1x,i4,1x,f8.1,1x,f8.1,1x,i4,1x,f8.1,3x,a)') & 
                             ilevel,jval,istate,quanta_RWF(ilevel)%spin,sigma,ilambda,omega,trim(poten(istate)%name)
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
  end subroutine define_quanta_bookkeeping_RWF
  !






end module RWF
