module dipole

 use accuracy,     only : hik, ik, rk, ark, cl, out, vellgt, planck, avogno, boltz, pi, small_
 use diatom_module,only : job,Intensity,quantaT,eigen,basis,Ndipoles,dipoletm,duo_j0,fieldT,poten,three_j,jmin_global
 use timer,        only : IOstart,Arraystart,Arraystop,ArrayMinus,Timerstart,Timerstop,MemoryReport, &
                          TimerReport,memory_limit,memory_now
 use symmetry,     only : sym,correlate_to_Cs



 private
 public dm_tranint

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


 integer(ik) :: Neigenlevels
 !
 type(ElevelT),allocatable  :: Elevel(:)  

 !
 type(dipoleT),allocatable     :: wigner(:,:) ! Rotational component of the dipole moment matrix elements 


 !type(IntensityT),save :: intensity

contains
  
 subroutine dm_tranint

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
    if (info /= 0) stop 'dm_tranint allocation error: Jval - out of memory'
    !
    allocate(q_part(20,nJ), stat = info)
    if (info /= 0) stop 'dm_tranint allocation error: q_part - out of memory'
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
       ! read contraction indexes (icase,ilambda) for all J specified
       !
       !call read_contrind(nJ, Jval(1:nJ))
       !
       !if (allocated(bset_contr)) deallocate(bset_contr)
       !allocate(bset_contr(0:njval), stat = info)
       !if (info /= 0) stop 'read_contr_ind allocation error: bset_contr - out of memory'
       !
       !bset_contr(0)%Ndimen = Ntotalvib
       !bset_contr(0)%Ntotal = Ntotalvib
       !
       !allocate(bset_contr(0)%icontr(Ntotalvib),stat = info)
       !call ArrayStart('bset_contr',info,1,4)
       !
       !do iJ = 1,nJ
       !  !
       !  !bset_contr(iJ)%jval = jval(iJ)
       !  !
       !  Jrot = Jval(iJ)
       !  !
       !  ! CHECK SIZE OF ICONTR AND ALLOCATE IT???
       !  !
       !  Ndimen = eigen(iJ)%Ndimen
       !  Nlevels = eigen(iJ)%Nlevels
       !  !
       !  !allocate(bset_contr(iJ)%icontr(Ndimen),stat = info)
       !  !call ArrayStart('bset_contr',info,1,4)
       !  !
       !  !call check_point_eigenvectors('READ',iverbose,Jrot,Ntotal,bset_contr(iJ)%icontr)
       !  bset_contr(iJ)%Nlevels = Ntotal
       !  !
       !enddo
       !
       ! find correlation between indexes for J = 0 and J > 0
       ! We need this because the contracted matrix elements of the 
       ! dipole moment are computed on the J=0 conatracted basis functions. 
       ! When J/=0 the numbering (bookkeeping) of these basis set has changed and
       ! we need to find the correlation between the bookkeepings.  
       !
       !call index_correlation(nJ, Jval(1:nJ))
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
       call dm_intensity(Jval,iverbose)
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
 end subroutine dm_tranint
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
 ! Electric dipole moment intensities and transition moments calculations 
 !
 subroutine dm_intensity(Jval,iverbose)
    !
    implicit none
    !
    real(rk),intent(in)  :: Jval(:)
    integer(ik),intent(in)   :: iverbose

    integer(ik)    :: nJ,dimenmax
    integer(ik)    :: ilevelI, ilevelF
    integer(ik)    :: nlevelsG(sym%Nrepresen)
    integer(ik)    :: info,indI,indF,itransit,Ntransit,Nrepresen
    integer(ik)    :: igammaI,igammaF
    integer(ik)    :: dimenI,dimenF,nmax,parity_gu,isymI,isymF
    real(rk)       :: energyI, energyF, nu_if,linestr,ener_,linestr2
    real(rk)       :: tm,jI,jF,ddot
    logical        :: passed,passed_

    real(rk),    allocatable :: vecI(:), vecF(:)
    real(rk),allocatable     :: half_linestr(:)
    !
    integer(ik)  :: jind,nlevels
    !
    integer(ik)  :: iroot,NlevelsI,NlevelsF,nlower,k,k_,iLF,iflag_rich
    !
    integer(ik)  :: igamma_pair(sym%Nrepresen),igamma,istateI,istateF,ivibI,ivibF,ivI,ivF,ilambdaI,ilambdaF,iparityI,itau
    integer(ik)  :: ivF_,ilambdaF_
    real(rk)     :: spinI,spinF,omegaI,omegaF,sigmaI,sigmaF,sigmaF_,omegaF_,spinF_
    integer(hik) :: matsize
    !
    character(len=1) :: branch,ef,pm
    character(len=2) :: dir
    character(len=10) :: statename
    !
    type(quantaT),pointer  :: quantaI,quantaF
    !
    real(rk)     :: boltz_fc, beta, intens_cm_mol, emcoef, A_coef_s_1, A_einst, absorption_int,lande
    !
    character(len=130) :: my_fmt !format for I/O specification
    integer :: ndecimals
    integer(ik)  :: enunit,transunit
    character(len=cl) :: filename,ioname
    !
    logical     :: integer_spin = .true.
    !
    integer(ik) :: alloc_p
    !
    integer(ik) :: Jmax_,ID_J
    real(rk) :: J_
    character(len=12) :: char_Jf,char_Ji,char_LF
    integer(ik),allocatable :: richunit(:,:)
    character(1)  :: let_LF ! richmol letters x,y,z
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
      filename =  trim(intensity%linelist_file)//'.states'
      write(ioname, '(a, i4)') 'Energy file '
      call IOstart(trim(ioname),enunit)
      open(unit = enunit, action = 'write',status='replace' , file = filename)
      !
      filename =  trim(intensity%linelist_file)//'.trans'
      write(ioname, '(a, i4)') 'Transition file '
      call IOstart(trim(ioname),transunit)
      open(unit = transunit, action = 'write',status='replace' , file = filename)
      !
      if ( intensity%matelem ) then
        !
        Jmax_ = nint(maxval(Jval(:))) 
        !
        allocate(richunit(nJ,nJ))
        !
        do indI = 1, nJ
          !
          jI = Jval(indI)
          write(char_jI,'(i12)') nint(Ji)
          do indF = 1, nJ
            !
            jF = Jval(indF)
            !
            if (Jf<Ji) cycle 
            ! 
            ! selection rules: Delta J<=1
            !
            if (nint(abs(jI-jF))>1.or.nint(jI+jF)==0) cycle
            !
            write(char_Jf,'(i12)') nint(jF)
            !
            !  New RICHMOL format - one file for x,y,z
            !
            filename =  &
            "matelem_MU"//"_j"//trim(adjustl(char_jI))//"_j"//trim(adjustl(char_jF))//"_"//trim(intensity%linelist_file)//".rchm"
            !
            call IOstart(trim(filename),richunit(indI,indF))
            open(unit = richunit(indI,indF), action = 'write',status='replace' , file = filename)
            !
            write(richunit(indI,indF),"('Start richmol format')")
            !
            write(richunit(indI,indF),"('MU','   1','   3')")
            write(richunit(indI,indF),"('M-tensor')")
            !
            ! three cartesian LF-components 
            do iLF = 1,3
              !
              ! set let_LF FOR EACH COMPONENT
              !
              let_LF = "x" ; if (iLF==2) let_LF = "y" ; if (iLF==3) let_LF = "z"
              !
              iflag_rich = -1 ; if (iLF==2) iflag_rich = 0 ! NEW RICHMOL FORMAT
              !
              write(char_LF,'(i12)') iLF
              !
              ! this flag is 0 for real M-values, 1 for imaginary and -1 for i**2 (OLD FORMAT)
              ! iflag_rich = 0 ; if (iLF==2) iflag_rich = 1 ! OLD RICHMOL FORMAT
              write(richunit(indI,indF),"('alpha',i5,i3,1x,a1)") iLF,iflag_rich,let_LF
              !
              ! count the number of M-elements ! OLD RICHMOL FORMAT
              ! call do_LF_matrix_elements(iLF,richunit(indI,indF,iLF),jI,jF,nMs) ! OLD RICHMOL FORMAT
              !write(richunit(indI,indF,iLF),"(i7)") nMs ! OLD RICHMOL FORMAT
              !
              call do_LF_matrix_elements(iLF,richunit(indI,indF),jI,jF) ! write elements
              !
              !write(richunit(indI,indF,iLF),"('MF-block')") ! OLD RICHMOL FORMAT
              !
            enddo
            !
            write(richunit(indI,indF),"('K-tensor')") ! NEW RICHMOL FORMAT
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
                  if ( intensity%matelem ) call matelem_filter (jI,jF,energyI,energyF,isymI,isymF,igamma_pair,passed)
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
           if (trim(intensity%linelist_file)/="NONE") then
             !
             !dimension of the bases for the initial states
             !
             !energy, quanta, and gedeneracy order of the initial state
             
             quantaI => eigen(indI,igammaI)%quanta(ilevelI)
             ivibI    = quantaI%ivib
             ivI      = quantaI%v
             sigmaI   = quantaI%sigma
             spinI    = quantaI%spin
             ilambdaI = quantaI%ilambda
             omegaI   = quantaI%omega
             iparityI = quantaI%iparity
             statename = trim(quantaI%name)
             !
             ! reconstruct the +/- and e/f parities:
             !
             pm = "+" ; if (iparityI==1) pm = "-"
             ef = "e" ; 
             !
             if ( mod( nint( 2.0*jI ),2 )==1 ) then
               itau = mod( nint( jI-0.5 ),2 )
             else
               itau = mod( nint( jI ),2 )
             endif
             !
             if (itau==iparityI) then 
               ef = "e"
             else
               ef = "f"
             endif
             !
             ! ndecimals contains the number of decimal digits to print for energy levels
             ! we use 6 decimals for energy levels up to 100,000 cm-1, then sacrify more and more decimals 
             ! to make larger numbers fit in 12 spaces.
             ! The present format works for energy levels larger than -10000 cm-1 and less than 1e11 cm-1
             ! By Lorenzo Lodi
             ndecimals=6-max(0, int( log10(abs(energyI-intensity%ZPE)+1.d-6)-4) )
             !
             !Mikhail Semenov: Lande g-factor for the selected eigenstate
             !
             if (intensity%lande_calc) then
               !
               lande = 0
               !
               vecI(1:dimenI) = eigen(indI,igammaI)%vect(1:dimenI,ilevelI)
               !
               !do k = 1,dimenI
               !  !
               !  omegaF = basis(indI)%icontr(k)%omega
               !  sigmaF = basis(indI)%icontr(k)%sigma
               !  ilambdaF = basis(indI)%icontr(k)%ilambda
               !  !
               !  if ( Ji > 0) then
               !    !
               !    lande = lande + vecI(k)**2*( omegaF+sigmaF )*omegaF/real(Ji*(Ji + 1),rk)
               !    !
               !  endif
               !  !
               !enddo
               !
               ! This version of Lande factor takes into account the non-diagonal Sigma/Sigma+/-1 terms
               !
               if (Ji>0) then
                 !
                 do k = 1,dimenI
                   !
                   omegaF   = basis(indI)%icontr(k)%omega
                   sigmaF   = basis(indI)%icontr(k)%sigma
                   ilambdaF = basis(indI)%icontr(k)%ilambda
                   spinF    = basis(indI)%icontr(k)%spin
                   ivF      = basis(indI)%icontr(k)%ivib
                   !
                   do k_ = 1,dimenI
                     !
                     omegaF_   = basis(indI)%icontr(k_)%omega
                     sigmaF_   = basis(indI)%icontr(k_)%sigma
                     ilambdaF_ = basis(indI)%icontr(k_)%ilambda
                     spinF_    = basis(indI)%icontr(k_)%spin
                     ivF_      = basis(indI)%icontr(k_)%ivib
                     !
                     if (ilambdaF/=ilambdaF_.or.nint(spinF-spinF_)/=0.or.ivF/=ivF_) cycle
                     !
                     if ( k==k_ ) then
                       !
                       ! g_e = 2.0023
                       !
                       lande = lande + vecI(k)*vecI(k)*( real(ilambdaF,rk)+2.0023_rk*sigmaF )*omegaF
                       !
                     elseif (nint(abs(sigmaF_-sigmaF))==1) then
                       !
                       lande = lande + vecI(k)*vecI(k_)*&
                               sqrt( spinF*(spinF+1.0_rk)-sigmaF*(sigmaF+sigmaF_-sigmaF) )*&
                               sqrt( Ji*(Ji+1.0_rk)-omegaF*(omegaF+omegaF_-omegaF) )*2.002319_rk/2.0_rk
                       !
                     endif
                     !
                   enddo
                   !
                 enddo
                 !
                 lande = lande/( Ji*(Ji + 1.0_rk) )
                 !
               endif
               !
               if (integer_spin) then 
                 !
                 write(my_fmt,'(A,i0,a)') "(i12,1x,f12.",ndecimals,",1x,i6,1x,i7,1x,f13.6,1x,a1,1x,a1,1x,a10,1x,i3,1x,i2,2i8)"
                 write(enunit,my_fmt) & 
                           iroot,energyI-intensity%ZPE,nint(intensity%gns(isymI)*( 2.0_rk*jI + 1.0_rk )),nint(jI),&
                           lande,pm,ef,statename,ivI,(ilambdaI),nint((sigmaI)),nint((omegaI))
                 !
               else
                 !
                 write(my_fmt,'(A,i0,a)') "(i12,1x,f12.",ndecimals,",1x,i6,1x,f7.1,1x,f13.6,1x,a1,1x,a1,1x,a10,1x,i3,1x,i2,2f8.1)"
                 write(enunit,my_fmt) & 
                           iroot,energyI-intensity%ZPE,nint(intensity%gns(isymI)*( 2.0_rk*jI + 1.0_rk )),jI,&
                           lande,pm,ef,statename,ivI,(ilambdaI),(sigmaI),(omegaI)
                           !
               endif
               !
             elseif ( intensity%matelem ) then
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
                 write(my_fmt,'(a)') &
                       "(i6,1x,i8,1x,i2,1x,i2,3x,e21.14,5x,a4,i3,1x,a2,i4,1x,a2,f8.4,1x,i6,1x,i6,1x,i4,1x,i6,1x,a1,1x,a10)"
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
    allocate(half_linestr(dimenmax),stat=info)
    !
    call ArrayStart('half_linestr',info,size(half_linestr),kind(half_linestr))
    !
    !  The matrix where some of the eigenvectors will be stored
    !
    if (iverbose>=5) call MemoryReport
    !
    write(out,"(/a,a,a,a)") 'Linestrength S(f<-i) [Debye**2],',' Transition moments [Debye],'& 
                          &,'Einstein coefficient A(if) [1/s],','and Intensities [cm/mol]'
    !
    ! Prepare the table header
    !
    select case (trim(intensity%action))
      !
      case('ABSORPTION')
        !
        write(out,"(/t5,'J',t7,'Gamma <-',t18,'J',t21,'Gamma',t27,'Typ',t37,'Ei',t44,'<-',t52,'Ef',t64,'nu_if',&
                    &8x,'S(f<-i)',10x,'A(if)',12x,'I(f<-i)', &
                    &7x,'State v lambda sigma  omega <- State v lambda sigma  omega ')")
        dir = '<-'
        !
      case('EMISSION')
        !
        write(out,"(/t5,'J',t7,'Gamma ->',t18,'J',t21,'Gamma',t27,'Typ',t37,'Ei',t44,'->',t52,'Ef',t64,'nu_if',&
                    &8x,'S(i->f)',10x,'A(if)',12x,'I(i->f)', &
                    &7x,'State v lambda sigma  omega -> State v lambda sigma  omega ')")
        dir = '->'
        !
      case('TM')
        !
        write(out,"(/t4,'J',t6,'Gamma <-',t17,'J',t19,'Gamma',t25,'Typ',t35,'Ei',t42,'<-',t52,'Ef',t65,'nu_if',&
                    &10x,'TM(f->i)')")

       !
    end select
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
    ! the actual intensity calculations
    ! --------------------------------- 
    !
    itransit = 0
    !
    ! loop over initial states
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
                ! Check if it is really necessary to start the calculations for the given levelI -> jF, 
                ! i.e. one can skip the rest if no transitions will start from the given ilevelI and 
                ! finish anywehere at J= jF. 
                !
                passed = .false.
                !
                !loop over final states
                !
                do ilevelF = 1, nlevelsF
                  !
                  !energy and and quanta of the final state
                  !
                  energyF = eigen(indF,igammaF)%val(ilevelF)
                  !
                  quantaF => eigen(indF,igammaF)%quanta(ilevelF)
                  istateF  = quantaF%istate
                  ivibF    = quantaF%ivib
                  ivF      = quantaF%v
                  sigmaF   = quantaF%sigma
                  spinF    = quantaF%spin
                  ilambdaF = quantaF%ilambda
                  omegaF   = quantaF%omega
                  !
                  ! reconstruct the symmetry for the C2v case which is different from Cs
                  parity_gu = poten(istateF)%parity%gu
                  isymF = correlate_to_Cs(igammaF,parity_gu)
                  !
                  !call TimerStart('Intens_Filter-2')
                  !
                  call intens_filter(jI,jF,energyI,energyF,isymI,isymF,igamma_pair,passed)
                  if ( intensity%matelem ) call matelem_filter (jI,jF,energyI,energyF,isymI,isymF,igamma_pair,passed)
                  !
                  !call TimerStop('Intens_Filter-2')
                  !
                  if (passed) exit
                  !
                enddo
                !
                if (.not.passed) cycle
                !
                select case (trim(intensity%action))
                  !
                case('ABSORPTION','EMISSION')
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
                case('TM')
                  !
                  call do_1st_half_tm(indI,indF,dimenI,dimenF,&
                                      vecI(1:dimenI),half_linestr)
                  !
                end select
                !
                !loop over final states
                !
                !$omp parallel private(vecF,alloc_p)
                allocate(vecF(dimenmax),stat = alloc_p)
                if (alloc_p/=0) then
                    write (out,"(' dipole: ',i9,' trying to allocate array -vecF')") alloc_p
                    stop 'dipole-vecF - out of memory'
                end if
                !
                !$omp do private(ilevelF,energyF,dimenF,quantaF,istateF,ivibF,ivF,sigmaF,spinF,ilambdaF,omegaF,passed,&
                !$omp& parity_gu,isymF,branch,nu_if,linestr,linestr2,A_einst,boltz_fc,absorption_int,tm) schedule(static) &
                !$omp                                                                             & reduction(+:itransit)
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
                   if ( intensity%matelem ) call matelem_filter (jI,jF,energyI,energyF,isymI,isymF,igamma_pair,passed)
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
                   select case (trim(intensity%action))
                     !
                     case default
                       !
                       stop 'only ABSORPTION and TM are properly coded up to now'
                       !
                     case('ABSORPTION','EMISSION')
                       !
                       linestr = ddot(dimenF,half_linestr,1,vecF,1)
                       !
                       linestr2 = linestr**2
                       !
                       ! calculate the intensity 
                       !
                       A_einst = A_coef_s_1*(2.0_rk*jI+1.0_rk)*linestr2*abs(nu_if)**3
                       !
                       linestr2 = linestr2 * intensity%gns(isymI) * ( 2.0_rk*jI + 1.0_rk )*( 2.0_rk * jF + 1.0_rk )
                       !
                       if (trim(intensity%action)=='ABSORPTION') then 
                         !
                         boltz_fc = exp( -(energyI-intensity%ZPE)*beta )*(1.0_rk - exp( -abs(nu_if)*beta) ) &
                                    / ( intensity%part_func*nu_if**2 )
                         !
                         ! intensity in cm/mol
                         !
                         absorption_int = 1.0_rk/(8.0_rk*pi*vellgt)*intensity%gns(isymF)*(2.0_rk*jF+1.0_rk)*A_einst*boltz_fc
                         !
                       else
                         !
                         ! emissivity in Ergs/mol/Sr
                         !
                         boltz_fc=exp( -(energyF-intensity%ZPE)*beta )
                         !
                         absorption_int = emcoef*A_einst*boltz_fc*intensity%gns(isymI)*(2.0_rk*jF+1.0_rk)*nu_if & 
                                        &      /intensity%part_func
                         !
                       endif
                       !
                       !
                       if (absorption_int>=intensity%threshold%intensity.and.linestr2>=intensity%threshold%linestrength) then 
                         !
                         !$omp critical
                         write(out, "( (f5.1, 1x, a4, 3x),a2, (f5.1, 1x, a4, 3x),a1,&
                                      &(2x, f11.4,1x),a2,(1x, f11.4,1x),f11.4,2x,&
                                      & 3(1x, es16.8),&
                                      & ' ( ',i2,1x,i3,1x,i2,2f8.1,' )',a2,'( ',i2,1x,i3,1x,i2,2f8.1,' )')")  &
                                      jF,sym%label(isymF),dir,jI,sym%label(isymI),branch, &
                                      energyF-intensity%ZPE,dir,energyI-intensity%ZPE,nu_if,  &
                                      linestr2,A_einst,absorption_int,&
                                      istateF,ivF,ilambdaF,sigmaF,omegaF,dir,&
                                      istateI,ivI,ilambdaI,sigmaI,omegaI
                                      !
                         !
                         ! generate the line list (Transition file)
                         !
                         if (trim(intensity%linelist_file)/="NONE") then
                           !
                           if ( intensity%matelem ) then 
                             !
                             write(richunit(indI,indF),"(i8,i8,2i3,4x,e24.14)") & 
                                               quantaI%iJ_ID,quantaF%iJ_ID,1,1,linestr
                             !
                           else
                             !
                             write(transunit,"(i12,1x,i12,2x,es10.4,4x,f16.6)") & 
                                       quantaF%iroot,quantaI%iroot,A_einst,nu_if
                           endif 
                           !
                         endif
                         !
                         !$omp end critical
                         !
                       endif
                       !
                     case('TM')
                       !
                       tm = dot_product( half_linestr(1:dimenF),vecF(1:dimenF) )
                       !
                       linestr = tm
                       !
                       if (linestr>=intensity%threshold%intensity) then 
                         !
                         !$omp critical
                         write(out, "( (i4, 1x, a3, 3x),'->', (i4, 1x, a3, 3x),a1, &
                                      &(2x, f13.6,1x),'->',(1x, f13.6,1x),f12.6, &
                                      &f15.8)") &
                                      !
                                      jI,sym%label(isymI),jF,sym%label(isymF),branch, &
                                      linestr,itransit,tm
                         !$omp end critical
                       endif
                       !
                   end select
                   !
                end do Flevels_loop
                !$omp enddo
                !
                deallocate(vecF)
                !$omp end parallel
                !
                if (iverbose>=5) call TimerReport
                !
              enddo Ilevels_loop
              !
           enddo
           !
         enddo
         !
       enddo
       !
       ! close some J-files for Richmol:
       if ( intensity%matelem ) then 
         !
         do indF = max(1,indI-1),min(nJ,indI+1)
           !
           jF = Jval(indF)
           if (nint(abs(jI-jF))>1.or.nint(jI+jF)==0) cycle
           if (jI>jF) cycle
           !
           write(richunit(indI,indF),"('End richmol format')")
           close(richunit(indI,indF))
           !
         enddo
         !
       endif
       !
    enddo
    !
    deallocate(vecI)
    call ArrayStop('intensity-vectors')
    !
    deallocate(half_linestr)
    call ArrayStop('half_linestr')
    !
    if (trim(intensity%linelist_file)/="NONE") close(transunit,status="keep")
    !
    call TimerStop('Intensity calculations')
    !
  end subroutine dm_intensity
  !
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
      newmin=idnint(max((a+be-c),(b-c-al),0.0_rk))
      newmax=idnint(min((a-al),(b+be),(a+b-c)))
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
      iphase=idnint(a-b-ga)
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





end module dipole
