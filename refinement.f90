module refinement
  !
  use accuracy
  use timer
  !
  use functions,only : define_complex_analytic_field_subterms
  use diatom_module,only : verbose,fitting,Nobjects,Nestates,Nspinorbits,&
                           Ntotalfields,fieldT,poten,spinorbit,l2,lxly,NL2,NLxLy,Nbobrot,Ndiabatic,Nlambdaopq,&
                           Nlambdap2q,Nlambdaq,Nnac,&
                           grid,duo_j0,quantaT,fieldmap,Nabi,abinitio,quadrupoletm, &
                           action,spinspin,spinspino,spinrot,bobrot,diabatic,lambdaopq,lambdap2q,lambdaq,nac,linkT,vmax,&
                           l_omega_obj,s_omega_obj,rangeT
  !
  !
  implicit none
  !
  private 
   !
   real(rk)    :: stab_best=1e-12  ! best standard error and stability 
   integer(ik) :: maxiter_as = 3                  ! maximal number of iterations to find a match for assignement
   integer(ik) :: Nobjectmax = 14
   integer(ik) :: vmax_ = 100                     ! if vmax in input is undefined vmax_ will be used to predife size
   !                                              !  of the matrix with energies 
   !
   real(rk), allocatable,save :: Jval(:)
   integer(ik),allocatable,save :: Jindex(:)
   !integer(ik)  ::vmax = 100
   !
   type fit_indexT
     integer(ik)  :: i
     integer(ik)  :: iobject
     integer(ik)  :: iterm
     integer(ik)  :: ifield
   end type fit_indexT
   !
   type object_containerT
     type(fieldT), pointer :: field
   end type object_containerT
   !
  public sf_fitting,define_jlist
  !
  !
  contains
  !
 subroutine sf_fitting
      !
      integer(ik)       :: npts,en_npts,pot_npts,nused   ! number of data points: current, all obs. energies, 
      !                                                       all pot. points, and actually used. 
      integer(ik)       :: parmax                  ! total number of potential parameters 
      real(rk)          :: jrot,jrot_              ! current value of the rotational quantum number 
      integer(ik)       :: lwork,fititer
      integer           :: info, ierror      ! error state variables 
      real(rk)          :: wtsum, lock_factor

      real(rk),allocatable :: wt_bit(:),wtall(:) ! weight factors - only for the energies and total.
      real(rk),allocatable :: rjacob(:,:),eps(:),energy_(:,:,:),enercalc(:)
      real(rk),allocatable :: local(:,:),pot_values(:),wspace(:),Tsing(:)
      real(rk),allocatable :: al(:,:),bl(:),dx(:),ai(:,:),sterr(:),sigma(:),bi(:),solution(:)
      real(rk),allocatable :: am(:,:),bm(:)
      character(len=cl),allocatable :: nampar(:)    ! parameter names 
      real(ark),allocatable :: pot_terms(:)
      !
      integer(ik),allocatable :: ivar(:),ifitparam(:),iZPE(:),nenergies(:,:),nenergies_(:,:)
      !
      logical      :: still_run,do_deriv,deriv_recalc,do_Armijo,do_print
      real(rk)     :: stadev_old,stability,stadev,sum_sterr,conf_int
      real(rk)     :: ssq,rms,ssq1,ssq2,rms1,rms2,fit_factor
      real(rk)     :: a_wats = 1.0_rk,lambda = 0.01_rk,nu = 10.0_rk
      integer(ik)  :: i,numpar,itmax,j,jlistmax,rank,ncol,nroots,nroots_max,NumSVD
      integer(ik)  :: iener,jener,irow,icolumn,ndigits,nrot,irot,irot_
      integer(ik)  :: enunit,abinitunit,i0,i1,iter_th,k0,frequnit
      integer(ik)  :: nlevels,NArmijo,ialpha_Armijo
      real(rk)     :: alpha_Armijo
      character(len=cl) :: filename, ioname
      character(len=330) :: char_fmt
      !
      real(rk),allocatable  ::  deriv(:,:,:,:)   ! to store the calculated derivatives 
      real(rk),allocatable  :: r(:),ezero(:)
      real(rk),allocatable :: potparam(:)
      character(len=1),allocatable  :: mark(:)
      character(len=1),parameter    :: pm(1:2) = (/"+","-"/)
      character(len=1)  :: mark_
      character(len=3)  :: mark_f
      character(len=3) :: iTAG,iTAGC,iTAG_,iTAGC_
      !
      integer(ik)  :: ifield,jfield,ifield_,iobject,Nfields,iterm,ifitpar,istate,ilambda,ivib,istate_,ilambda_,ivib_,ngrid,Nterms
      integer(ik)  :: maxNfields,ipotmin,i_,nchar,k,k_
      real(rk)     :: JrotC,sigmaC,omegaC,spinC,JrotC_,sigmaC_,omegaC_,spinC_,rms0
      integer(ik)  :: istateC,istateC_,ilambdaC,ilambdaC_,ivibC,ivibC_,itau,itau_,itauC,itauC_
      !
      real(rk)     :: param_t,delta,sigmai,omega,spin,sigmai_,omega_,spin_,f,fR,fL,r_t,req,fshift
      type(fieldT),pointer      :: field
      type(linkT),pointer       :: flink
      type(rangeT),pointer      :: frange
      real(rk),allocatable      :: EnergyR(:,:,:),EnergyL(:,:,:),params_t(:)
      type(quantaT),allocatable   :: calc(:,:,:),calc_(:,:,:)
      type(fit_indexT),allocatable   :: fit_index(:)
      type(object_containerT), allocatable :: objects(:,:)
      type(fieldT),allocatable      :: object0(:,:)
      character(len=130)  :: my_fmt ! contains format specification for intput/output
      logical :: printed ! used to print frequencies
      real(rk),allocatable    :: spline_wk_vec(:) ! working vector needed for spline interpolation
      real(rk) :: yp1, ypn
       !
       if (verbose>=2) write(out,"(/'The least-squares fitting ...')")
       !
       ! Nobjectmax = Nobjects - 3
       Nobjectmax = 14
       !
       call TimerStart('Simultaneous Fitting')
       !
       en_npts = fitting%Nenergies
       !
       ngrid = grid%npoints
       !
       filename =  trim(fitting%output_file)//'.en'
       write(ioname, '(a, i4)') 'calc. energies from the fit '
       call IOstart(trim(ioname), enunit)
       open(unit = enunit, action = 'write',status='replace' , file = filename)
       !
       ! file (temp) where all computed potential energy corrections are printed out 
       !
       filename =  trim(fitting%output_file)//'.pot'
       write(ioname, '(a, i4)') 'pot. energies at the ab initio geometries'
       call IOstart(trim(ioname), abinitunit)
       open(unit = abinitunit, action = 'write',status='replace' , file = filename)
       !
       if (action%frequency) then
         !
         filename =  trim(fitting%output_file)//'.freq'
         write(ioname, '(a, i4)') 'calc. frequencies from the fit '
         call IOstart(trim(ioname), frequnit)
         open(unit = frequnit, action = 'write',status='replace' , file = filename)
         !
       endif
       !
       pot_npts = 0
       !
       ! count all grid points 
       !
       j = 0
       do ifield = 1,Nabi !Ntotalfields
         !
         field => abinitio(ifield)
         !
         do i = 1,field%Nterms
            j = j + 1
         enddo 
       enddo
       !
       pot_npts = j
       !
       write(out,"(/'Number of obs. energies: ',i9)") en_npts
       write(out,"(/'Number of potential energy data points: ',i9)") pot_npts
       !
       npts = en_npts + pot_npts
       !
       allocate (r(pot_npts),stat=info)
       call ArrayStart('coords-mat',info,size(local),kind(local))
       !
       ! count all parameters used to define different fields 
       !
       parmax = fitting%parmax
       !
       itmax = fitting%itermax
       !
       fit_factor = fitting%factor
       !
       jlistmax = size(Jval)
       !
       allocate (mark(en_npts),stat=info)
       allocate (wtall(npts),stat=info)
       call ArrayStart('weight-mat',info,size(wtall),kind(wtall))
       allocate (wt_bit(npts),stat=info)
       call ArrayStart('weight-mat',info,size(wt_bit),kind(wt_bit))
       !
       allocate (potparam(parmax),nampar(parmax),ivar(parmax),ifitparam(parmax),&
                 al(parmax,parmax),ai(parmax,parmax),bl(parmax),bi(parmax),dx(parmax),sterr(parmax),&
                 Tsing(parmax),pot_terms(parmax),fit_index(parmax),solution(parmax),stat=info)
       call ArrayStart('potparam-mat',info,size(potparam),kind(potparam))
       call ArrayStart('potparam-mat',info,size(ivar),kind(ivar))
       call ArrayStart('potparam-dx',info,size(dx),kind(dx))
       call ArrayStart('potparam-sterr',info,size(sterr),kind(sterr))
       call ArrayStart('potparam-Tsing',info,size(Tsing),kind(Tsing))
       call ArrayStart('a-b-mat',info,size(al),kind(al))
       call ArrayStart('a-b-mat',info,size(bl),kind(bl))
       call ArrayStart('pot_terms',info,size(pot_terms),kind(pot_terms))
       !
       allocate (am(parmax,parmax),bm(parmax),stat=info)
       call ArrayStart('a-b-mat',info,size(am),kind(am))
       call ArrayStart('a-b-mat',info,size(bm),kind(bm))
       !
       allocate (pot_values(pot_npts),stat=info)
       call ArrayStart('pot_values-mat',info,size(pot_values),kind(pot_values))
       !
       ! rough estimation of the total number of energies 
       !
       allocate (enercalc(en_npts),stat=info)
       call ArrayStart('enercalc',info,size(enercalc),kind(enercalc))
       allocate (ezero(Nestates),iZPE(Nestates),stat=info)
       call ArrayStart('ezero',info,size(ezero),kind(ezero))
       call ArrayStart('iZPE',info,size(iZPE),kind(iZPE))
       !
       ! total number of jval-ues
       nrot = size(fitting%J_list)
       !
       ! arrays to count and store energies within each J block
       allocate(nenergies(nrot,2),nenergies_(nrot,2),stat=info)
       call ArrayStart("nenergies",info,size(nenergies),kind(nenergies))
       !
       ! estimate maximal number of energies for each J
       !
       Nfields = fieldmap(1)%Nfields
       !
       sigmai = 0
       ilambda = 0
       do ifield =1,Nfields
          !
          sigmai = max(poten(ifield)%spini,sigmai)
          ilambda = max(poten(ifield)%lambda,ilambda)
          !
       enddo
       !
       ! maximal value of omega is
       omega = sigmai+real(ilambda)
       !
       ! we estimate the size of the basis as the product of (vmax+1) for vibration,
       !                           (2*sigma+1)x2 for spin-rotation, Nestates for electronic 
       ! all divided by two to account for two symmetry blocks
       !
       if (vmax>0) vmax_ = vmax
       !
       nroots_max = (vmax_+1)*Nestates*2*(2*nint(sigmai)+1)/2+1
       allocate (energy_(nrot,2,nroots_max),calc(nrot,2,nroots_max),stat=info)
       call ArrayStart('energy_',info,size(energy_),kind(energy_))
       !
       ! Allocate objects, that will be used for the fitting procedure:
       !
       allocate (rjacob(npts,parmax),eps(npts),deriv(nrot,2,nroots_max,parmax),stat=info)
       call ArrayStart('rjacob',info,size(rjacob),kind(rjacob))
       call ArrayStart('eps',info,size(eps),kind(eps))
       call ArrayStart('fit-ener-deriv',info,size(deriv),kind(deriv))
       !
       allocate (params_t(parmax),stat=info)
       call ArrayStart('params_t',info,size(params_t),kind(params_t))
       !
       if (fitting%robust>small_) then
         !
         allocate (sigma(npts),stat=info)
         call ArrayStart('sigma',info,size(sigma),kind(sigma))
         !
       endif 
       !
       ! The last object to allocate - the lapack related work array
       !
       lwork = 50*parmax
       !
       allocate (wspace(lwork),stat=info)
       call ArrayStart('wspace',info,size(wspace),kind(wspace))
       !
       wtall = 0 
       !
       if (verbose>=2)  write (out,"('...done!')") 
       !
       forall(i=1:en_npts) wtall(i) = fitting%obs(i)%weight
       wtsum = sum(fitting%obs(1:en_npts)%weight) ; wtsum = max(wtsum,sqrt(small_))
       wtall(1:en_npts) = wtall(1:en_npts)/max(wtsum,sqrt(small_))
       !
       ! weights for the ab initio fields
       !
       !j = 0
       !do ifield =1,Nabi
       !  !
       !  field => abinitio(ifield)
       !  do i=1,field%Nterms
       !    j = j + 1
       !    r(j) = field%grid(i)
       !    wtall(en_npts+j) = field%weight(i)
       !  enddo
       !  !
       !  ! normalize weights for each field separately
       !  !
       !  wtsum = sum(wtall(en_npts+j-field%Nterms+1:en_npts+j))
       !  wtall(en_npts+j-field%Nterms+1:en_npts+j) =  & 
       !       field%fit_factor*wtall(en_npts+j-field%Nterms+1:en_npts+j)/max(wtsum,sqrt(small_))
       !  !
       !enddo 
       !
       !wtsum = sum(wtall(en_npts+1:npts))
       !wtall(en_npts+1:npts) = wtall(en_npts+1:npts)/max(wtsum,sqrt(small_))
       !
       ! assigning mark with null
       !
       mark(:) = ' '
       !
       ! The fitting loop is about to start. 
       !
       ! fititer will count the iterations. 
       !
       fititer = 0
       !
       nlevels = fitting%Nenergies
       !
       ! find the maximal size of fields
       !
       maxNfields = 1
       do iobject =1,Nobjectmax
          maxNfields = max(fieldmap(iobject)%Nfields,maxNfields)
       enddo
       !
       allocate(objects(Nobjectmax,maxNfields),stat=info)
       allocate(object0(Nobjectmax,maxNfields),stat=info)
       !
       ! find eqilibrium
       !
       ipotmin = minloc(poten(1)%gridvalue,dim=1) ; req = grid%r(ipotmin)
       !
       ifield_ = 0
       j = 0
       !
       do iobject =1,Nobjectmax
          !
          Nfields = fieldmap(iobject)%Nfields
          !
          do ifield =1,Nfields
            !
            select case (iobject)
            case (1)
              objects(iobject,ifield)%field => poten(ifield)
            case (2)
              objects(iobject,ifield)%field => spinorbit(ifield)
            case (3)
              objects(iobject,ifield)%field => l2(ifield)
            case (4)
              objects(iobject,ifield)%field => lxly(ifield)
            case (5)
              objects(iobject,ifield)%field => spinspin(ifield)
            case (6)
              objects(iobject,ifield)%field => spinspino(ifield)
            case (7)
              objects(iobject,ifield)%field => bobrot(ifield)
            case (8)
              objects(iobject,ifield)%field => spinrot(ifield)
            case (9)
              objects(iobject,ifield)%field => diabatic(ifield)
            case (10)
              objects(iobject,ifield)%field => lambdaopq(ifield)
            case (11)
              objects(iobject,ifield)%field => lambdap2q(ifield)
            case (12)
              objects(iobject,ifield)%field => lambdaq(ifield)
            case (13)
              objects(iobject,ifield)%field => nac(ifield)
            case (14)
              objects(iobject,ifield)%field => quadrupoletm(ifield)
            end select
            !
            Nterms = objects(iobject,ifield)%field%Nterms
            allocate(object0(iobject,ifield)%value(Nterms))
            call ArrayStart("object0",info,Nterms,kind(object0(iobject,ifield)%value))
            !
            object0(iobject,ifield)%value = objects(iobject,ifield)%field%value
            !
            ! what is the number of the parent field in the abinitio counter?
            !
            ifield_  = objects(iobject,ifield)%field%iabi
            !
            if (ifield_==0) cycle 
            !
            ! change to the equilibrium
            !
            if (abinitio(ifield_)%Nterms==1) then
              !
              abinitio(ifield_)%grid(1) = req
              !
              !f = fanalytic_field(req,objects(iobject,ifield)%field%type,objects(iobject,ifield)%field%value)
              !
              abinitio(ifield_)%gridvalue(1) = objects(iobject,ifield)%field%gridvalue(ipotmin)
              !
            endif 
            !
            ! weights for the ab initio fields
            !
            do i = 1,abinitio(ifield_)%Nterms
              !
              j = j + 1
              !r(j) = field%grid(i)
              r(j)=  abinitio(ifield_)%grid(1)
              !wtall(en_npts+j) = field%weight(i)
              wtall(en_npts+j) = abinitio(ifield_)%weight(i)
            enddo
            !
            wtsum = sum(wtall(en_npts+j-field%Nterms+1:en_npts+j))
            wtall(en_npts+j-abinitio(ifield_)%Nterms+1:en_npts+j) =  & 
                 abinitio(ifield_)%fit_factor*wtall(en_npts+j-abinitio(ifield_)%Nterms+1:en_npts+j)/max(wtsum,sqrt(small_))
            !
         enddo
       enddo
       !
       ! 2. Factorizing the obs. weights by the factor "fit_factor":
       !
       wtall(1:en_npts) = wtall(1:en_npts)*fitting%factor
       !
       ! 3. And normilizing all weight factors. 
       !
       wtsum = sum(wtall(1:npts))
       wtall(1:npts) = wtall(1:npts)/max(wtsum,sqrt(small_))
       ! 
       ! Count how many data points actually will be fitted. 
       !
       nused=0
       wt_bit = 0 
       !
       do i=1,npts
         if (wtall(i)>small_) then 
           !
           nused=nused+1
           wt_bit(i) = 1.0_rk
           !
         endif
       enddo
       !
       write(out,"('Number of data points used in the fit: ',i9)") nused
       !
       ! sigma = exp. data precision for robust fit 
       !
       if (fitting%robust>small_) then
         !
         sigma = 1.0_rk
         !
         !$omp parallel do private(i) shared(sigma) schedule(dynamic)
         do i=1,npts
           if (wtall(i)>small_) sigma(i) = sigma(i)/sqrt(wtall(i))
         enddo
         !$omp end parallel do
         !
         wtsum = 1.0_rk ! sqrt(sum(sigma(1:en_npts)**2))
         !
         sigma(:) = sigma(:)*fitting%robust/wtsum
         !
       endif 
       !
       ! For "frequencies" reconstruct energies by combine differences
       !
       if (action%frequency) then
           !
           ! TODO
           !
           do iener = 1,en_npts
             !
             jrot = fitting%obs(iener)%jrot
             itau = fitting%obs(iener)%iparity+1
             !
             jrot_= fitting%obs(iener)%jrot_
             itau_= fitting%obs(iener)%iparity_+1
             !
           enddo
           !
       endif
       !
       !
       ! The outer fitting loop - allows to restart the fitting in case 
       ! we decide to remove some of the varying parameters. 
       ! This option is working together with fit_type ='linur'
       !
       still_run = .true.
       outer_loop: do while (still_run)  
          !
          ! Initial values for the standard error and  stability.
          !
          stadev_old = 1e10
          stability = 1e10
          stadev    = 1e10
          !
          ! relation between the fitting parameters and the corresponding parameter in each field
          !
          numpar = 0
          ifield_ = 0
          j = 0
          !
          do iobject =1,Nobjectmax
             !
             Nfields = fieldmap(iobject)%Nfields
             !
             do ifield =1,Nfields
               !
               do iterm=1,objects(iobject,ifield)%field%Nterms
                 !
                 j = j + 1
                 !
                 if (nint(objects(iobject,ifield)%field%weight(iterm)).eq.1) then 
                   !
                   numpar=numpar+1
                   !
                   fit_index(numpar)%i = j
                   fit_index(numpar)%iobject = iobject
                   fit_index(numpar)%ifield = ifield
                   fit_index(numpar)%iterm = iterm
                   !
                 endif
                 !
               enddo
               !
            enddo
          enddo
          !
          if (numpar==0.and.itmax/=0) then 
            !
            write (out,"('No varying paramters, check input!')") 
            stop  ! 'sf_fitting: No varying paramters'
            !
          endif 
          !
          rjacob = 0
          !
          ! NArmijo iterations for linear search
          do_Armijo = .false.
          if (fitting%linear_search>0) do_Armijo = .true.
          NArmijo = fitting%linear_search
          rms0 = 0
          ialpha_Armijo = 0
          !
          do_print = .true.
          !
          ! printing out format 
          !
          if (action%frequency) then 
            char_fmt = "(/1X,250('-'),/'|  ## |  N |    J  p |  N |    J  p |      Obs.     |     Calc.   &
             &|  Obs.-Calc. |   Weight |      Eup      |     Elow    |  State vib Lambda  Sigma Omega State vib Lambda Sigma &
                       &Omega    State vib Lambda  Sigma Omega  State vib Lambda Sigma Omega',/1X,250('-'))"
          else
            char_fmt = "(/1X,145('-'),/'| ## |  N |     J  p |      Obs.      |     Calc.   |  Obs.-Calc. | &
                      &  Weight | States  vib  Lambda Sigma  Omega  States  vib Lambda Sigma &
                      & Omega ',/1X,145('-'))"
          endif 
          !
          ! The loop starts here. 
          !
          loop_fititer : do while( fititer.le.itmax .and. stadev.ge.fitting%target_rms.and. stability.ge.stab_best)
            !
            fititer = fititer + 1
            !
            ! If fit_factor is set zero - there would be no fit to the experimental energies. 
            !
            if (fit_factor<0) rjacob(1:en_npts,:) = 0
            !
            !if (do_print) write(out,"(/'Iteration = ',i8)") fititer-1
            if (do_print) write(enunit,"(/'Iteration = ',i8)") fititer-1
            if (action%frequency.and.do_print) write(frequnit,"(/'Iteration = ',i8)") fititer-1
            if (do_print.and.fititer-1==0) write(out,"(/a)") &
              'Iteration = 0: Straight through calculations with initial parameters...'
            if (do_print.and.fititer-1==0.and.itmax>0) write(out,"(/a)") 'Generating derivatives for the least-squares fit...'
            !
            ! Reconstruct the potential expansion from the local to linearized coords.
            !
            pot_terms = potparam
            !
            eps = 0
            !
            ! Print out the calc. and obs.-calc., i.e. result of the fit. 
            ! Only fitted energies are printed. 
            !
            if (do_print) write(out ,char_fmt)
            if (do_print) write(enunit,"(/1X,145('-'),/'| ## |  N |     J  p |      Obs.      |     Calc.   |  Obs.-Calc. | &
                        &  Weight | States  vib  Lambda Sigma  Omega  States  vib Lambda Sigma &
                        & Omega ',/1X,145('-'))")
            !
            deriv_recalc = .true.
            !
            ! some parameters are linked in the fitting:
            !
            call update_linked_parameters
            !
            ! first compute energies
            !
            call TimerStart('Energy calculations')
            !
            call duo_j0(0_ik,fitting%J_list,energy_,calc,nenergies)
            !
            if (maxval(nenergies)>nroots_max) then
              !
              write(out,"('refinement: vmax_ is too small ',i0,' nroots_max = ',i0,' nenergies ',i0)") vmax_,nroots_max, & 
                                                                                                        maxval(nenergies)
              stop 'refinement: vmax_ is too small'
              !
            endif 
            !
            call TimerStop('Energy calculations')
            !
            ! compute derivatives
            !
            call TimerStart('Energy derivatives')
            !
            ! Check we if we need to computed the derivatives 
            !
            do_deriv = .false.
            if (itmax.ge.1.and.fit_factor>0) do_deriv = .true. 
            !
            do_print = .true.
            ! Switch off derivatives for linear search (Armijo) iteration > 0
            if (do_Armijo.and.ialpha_Armijo/=0) then 
               do_deriv = .false. 
               do_print = .false.
            endif
            !
            if (do_deriv) then 
              !
              deriv = 0
              !
              nroots = 0
              do itau = 1,2
                !
                nroots = max(nroots,maxval(nenergies(:,itau),dim=1))
                !
                !do irot = 1,nrot
                !   nenergies(irot,itau) = min(size(nenergies,dim=3))
                !enddo
                !
              enddo
              !
              allocate (calc_(nrot,2,nroots),stat=info)
              allocate(energyR(nrot,2,nroots),energyL(nrot,2,nroots),stat=info)
              call ArrayStart('EnergyLR',info,size(EnergyL),kind(EnergyL))
              call ArrayStart('EnergyLR',info,size(EnergyR),kind(EnergyR))
              !
              do ifitpar = 1,numpar
                !
                i = fit_index(ifitpar)%i
                iobject = fit_index(ifitpar)%iobject
                ifield = fit_index(ifitpar)%ifield
                iterm = fit_index(ifitpar)%iterm
                !
                param_t = objects(iobject,ifield)%field%value(iterm)
                !
                delta = 1e-4*abs(param_t)
                !
                if (delta < small_) delta = sqrt(small_)*10.0
                !
                objects(iobject,ifield)%field%value(iterm) = param_t + delta
                !
                call update_linked_parameters
                !
                call duo_j0(0,fitting%J_list,EnergyR,calc_,nenergies_)
                !
                objects(iobject,ifield)%field%value(iterm) = param_t - delta
                !
                call update_linked_parameters
                !
                call duo_j0(0,fitting%J_list,EnergyL,calc_,nenergies_)
                !
                !omp parallel do private(irot,nroots) shared(deriv) schedule(dynamic)
                do irot = 1,nrot
                  !
                  do itau = 1,2
                    !
                    nroots = min(nenergies_(irot,itau),nenergies(irot,itau))
                    !
                    deriv(irot,itau,1:nroots,ifitpar) = (EnergyR(irot,itau,1:nroots)-EnergyL(irot,itau,1:nroots))/(2.0_rk*delta)
                    !
                  enddo
                  !
                enddo
                !omp end parallel do
                !
                objects(iobject,ifield)%field%value(iterm) = param_t
                !
                call update_linked_parameters
                !
              enddo
              !
              deallocate(energyR,energyL,calc_)
              call ArrayStop('EnergyLR')
              !
            endif
            !
            ! Nentries
            !
            call TimerStop('Energy derivatives')
            !
            ! Zero point energy:
            ! obtain the lowest J=0 energies (ZPEs) for each electronic state (if possible)
            !
            ezero = fitting%zpe
            if (fitting%shift_to_zpe) then 
               ezero = energy_(1,1,1)
               if (energy_(1,2,1)<ezero(1).and.nenergies(1,2)>0) ezero = energy_(1,2,1)
            endif
            !
            ! find the lowest root for each electronic state  - will be used as a corresponding ZPE
            !
            loop_istates : do istate  =1,Nestates
              !
              do i = 1,nenergies(1,1)
                !
                if (istate == calc(1,1,i)%istate.and.nint(2.0*calc(1,1,i)%Jrot)==nint(2.0*fitting%J_list(1))) then 
                  !
                  ezero(istate) = fitting%zpe
                  if (fitting%shift_to_zpe) then 
                    ezero(istate) = energy_(1,1,i)
                    if (energy_(1,2,i)<ezero(istate).and.nenergies(1,2)>0) ezero(istate) = energy_(1,2,i)
                  endif
                  !
                  iZPE(istate) = i
                  !
                  cycle loop_istates 
                  !
                endif
                !
              enddo
            enddo loop_istates
            !
            ! if threshold_lock>0, correct the addresses of the energies in case of accidential swaps
            ! by comparing with energies within the "threshold_lock"-range and with "exp" quantum numbers.
            ! if threshold_lock<0, find closets matches with experimental energies 
            ! if threshold_lock=0, do nothing
            !
            mark = " "
            !
            if (action%frequency) then
              !
              ! check swaps 
              !
              if (abs(fitting%threshold_lock)>0) then
                !
                ! upper state 
                !
                !k0 = 1
                !
                if (fitting%threshold_lock<0)  maxiter_as = 1
                !
                jrot  = fitting%obs(1)%jrot
                itau  = fitting%obs(1)%iparity+1
                !
                jrot_ = fitting%obs(1)%jrot_
                itau_ = fitting%obs(1)%iparity_+1
                !
                do iener = 1,en_npts
                  !
                  mark(iener) = '*'
                  !
                  !if (jrot /=fitting%obs(iener)%jrot .or.itau  /= fitting%obs(iener)%iparity +1.or.&
                  !    jrot_/=fitting%obs(iener)%jrot_.or.itau_ /= fitting%obs(iener)%iparity_+1) k0 = 1
                  !
                  jrot = fitting%obs(iener)%jrot
                  itau = fitting%obs(iener)%iparity+1
                  !
                  jrot_= fitting%obs(iener)%jrot_
                  itau_= fitting%obs(iener)%iparity_+1
                  !
                  loop_thresh_u : do iter_th = 1,maxiter_as
                    !
                    !if (iter_th == 0) cycle
                    !
                    lock_factor = real(iter_th,8)
                    !if (iter_th<) lock_factor = 10**iter_th
                    !
                    !lock_factor = 10.0**iter_th
                    !
                    loop_jrot_u : do irot = 1,nrot
                      !
                      jrotC = fitting%J_list(irot)
                      !
                      if (nint( jrotC-fitting%obs(iener)%jrot )/=0) cycle
                      !
                      loop_jrot_l : do irot_ = 1,nrot
                        !
                        jrotC_ = fitting%J_list(irot_)
                        !
                        if (nint( jrotC_-fitting%obs(iener)%jrot_ )/=0) cycle
                        !
                        do itauC = 1,2
                          !
                          if (itauC-1/=fitting%obs(iener)%iparity) cycle
                          !
                          do itauC_ = 1,2
                            !
                            if (itauC_-1/=fitting%obs(iener)%iparity_) cycle
                            !
                            loop_i : do i0 = 1,nenergies(irot,itauC)
                              !
                              loop_j : do i1 = 1,nenergies(irot_,itauC_)
                                !
                                if (energy_(irot,itauC,i0)-energy_(irot_,itauC_,i1)<0) cycle
                                !
                                !do jener = max(iener-10,1),iener-1
                                !  !
                                !  if (nint(jrotC -fitting%obs(jener)%jrot )/=0.or.itauC -1/=fitting%obs(jener)%iparity.or.&
                                !      nint(jrotC_-fitting%obs(jener)%jrot_)/=0.or.itauC_-1/=fitting%obs(jener)%iparity_) cycle
                                !  !
                                !  if ( fitting%obs(jener)%N_ == i1.and.fitting%obs(jener)%N == i0 ) cycle loop_j
                                !  !
                                !enddo
                                !
                                if ( abs( fitting%obs(iener)%energy-( energy_(irot,itauC,i0)-energy_(irot_,itauC_,i1) )  )<= & 
                                                            & lock_factor*abs(fitting%threshold_lock).and.&
                                      ! threshold_lock > 0, search for a QN match
                                    ( ( calc(irot,itauC,i0)%istate==fitting%obs(iener)%quanta%istate.and.&
                                        abs(calc(irot,itauC,i0)%ilambda)==abs(fitting%obs(iener)%quanta%ilambda).and.&
                                        calc(irot,itauC,i0)%v==fitting%obs(iener)%quanta%v.and.&
                                        nint(abs(calc(irot,itauC,i0)%sigma)-abs(fitting%obs(iener)%quanta%sigma) )==0.and.&
                                        nint(abs(calc(irot,itauC,i0)%omega)-abs(fitting%obs(iener)%quanta%omega) )==0.and.&
                                        !
                                        calc(irot_,itauC_,i1)%istate==fitting%obs(iener)%quanta_%istate.and.&
                                        abs(calc(irot_,itauC_,i1)%ilambda)==abs(fitting%obs(iener)%quanta_%ilambda).and.&
                                        calc(irot_,itauC_,i1)%v==fitting%obs(iener)%quanta_%v.and.&
                                        nint(abs(calc(irot_,itauC_,i1)%sigma)-abs(fitting%obs(iener)%quanta_%sigma))==0 .and.&
                                        nint(abs(calc(irot_,itauC_,i1)%omega)-abs(fitting%obs(iener)%quanta_%omega))==0 ).or.&
                                      ! threshold_lock < 0, only frequency match 
                                      (  fitting%threshold_lock<0 ) ) ) then 
                                    !
                                    fitting%obs(iener)%N = i0
                                    fitting%obs(iener)%N_= i1
                                    mark(iener) = ' '
                                    exit loop_thresh_u
                                    !
                                endif 
                                !
                              enddo loop_j
                              !
                            enddo loop_i
                            !
                          enddo
                          !
                        enddo
                        !
                      enddo loop_jrot_l
                      !
                    enddo loop_jrot_u
                    !
                  enddo loop_thresh_u
                  !
                  !if (mark(iener) == ' ')  k0 = max(fitting%obs(iener)%N-10,1)
                  !
                enddo
                !
              endif
              !
              !
              ! frequency fit
              !
              do iener = 1,en_npts
                 !
                 Jrot    = fitting%obs(iener)%Jrot
                 irot    = fitting%obs(iener)%irot
                 istate  = fitting%obs(iener)%quanta%istate
                 sigmai  = fitting%obs(iener)%quanta%sigma
                 ilambda = fitting%obs(iener)%quanta%ilambda
                 omega   = sigmai + real(ilambda,rk)
                 spin    = fitting%obs(iener)%quanta%spin
                 ivib    = fitting%obs(iener)%quanta%v
                 iTAG    = fitting%obs(iener)%quanta%iTAG(1:3)
                 !
                 Jrot_   = fitting%obs(iener)%Jrot_
                 irot_   = fitting%obs(iener)%irot_
                 istate_ = fitting%obs(iener)%quanta_%istate
                 sigmai_ = fitting%obs(iener)%quanta_%sigma
                 ilambda_= fitting%obs(iener)%quanta_%ilambda
                 omega_  = fitting%obs(iener)%quanta_%omega
                 spin_   = fitting%obs(iener)%quanta_%spin
                 ivib_   = fitting%obs(iener)%quanta_%v
                 iTAG_   = fitting%obs(iener)%quanta_%iTAG(1:3)
                 !
                 itau  = fitting%obs(iener)%iparity +1 ;  if (itau <1.or.itau >2) cycle 
                 itau_ = fitting%obs(iener)%iparity_+1 ;  if (itau_<1.or.itau_>2) cycle 
                 !
                 i  = fitting%obs(iener)%N ;  if (i >nenergies(irot ,itau )) cycle 
                 i_ = fitting%obs(iener)%N_;  if (i_>nenergies(irot_,itau_)) cycle 
                 !
                 enercalc(iener) = energy_(irot,itau,i)-energy_(irot_,itau_,i_)
                 !
                 eps(iener) = fitting%obs(iener)%energy-enercalc(iener)
                 !
                 if (do_deriv) then  
                    rjacob(iener,1:numpar) =   deriv(irot,itau,i ,1:numpar) - deriv(irot_,itau_,i_,1:numpar)
                 endif
                 !
                 if (do_print) then
                    !
                    JrotC    = calc(irot,itau,i)%Jrot
                    istateC  = calc(irot,itau,i)%istate
                    sigmaC   = calc(irot,itau,i)%sigma
                    ilambdaC = calc(irot,itau,i)%ilambda
                    omegaC   = sigmaC+real(ilambdaC,rk)
                    spinC    = calc(irot,itau,i)%spin
                    ivibC    = calc(irot,itau,i)%v
                    iTAGC    = calc(irot,itau,i)%iTAG(1:3)
                    !
                    JrotC_   = calc(irot_,itau_,i_)%Jrot
                    istateC_ = calc(irot_,itau_,i_)%istate
                    sigmaC_  = calc(irot_,itau_,i_)%sigma
                    ilambdaC_= calc(irot_,itau_,i_)%ilambda
                    omegaC_  = sigmaC_+real(ilambdaC_,rk)
                    spinC_   = calc(irot_,itau_,i_)%spin
                    ivibC_   = calc(irot_,itau_,i_)%v
                    iTAGC_   = calc(irot_,itau_,i_)%iTAG(1:3)
                    
                    write (out,"(i5,2(i5,1x,f7.1,1x,a1),2x,' ',3f14.4,2x,e9.2,2x,2f14.4,2x," &
                               // " '(',1x,a3,1x,2i4,2f7.1,', <-',1x,a3,1x,2i4,2f7.1,')', " &
                               // " '(',1x,a3,1x,2i4,2f7.1,', <-',1x,a3,1x,2i4,2f7.1,')',a)") &
                           iener,i,Jrot,pm(itau),i_,Jrot_,pm(itau_),enercalc(iener)+eps(iener),enercalc(iener),eps(iener),&
                           wtall(iener),energy_(irot,itau,i)-ezero(1),energy_(irot_,itau_,i_)-ezero(1),&
                           iTAGC,ivibC,ilambdaC,sigmaC,omegaC,iTAGC_,ivibC_,ilambdaC_,sigmaC_,omegaC_,&
                           iTAG ,ivib ,ilambda ,sigmai,omega ,iTAG_ ,ivib_ ,ilambda_ ,sigmai_,omega_ ,&
                           mark(iener)
                endif
                !
              enddo ! --- i
              !
              ! print more frequencies in order to predict mis-assignements
              !
              if (do_print) then
                 !
                 do iener = 1,en_npts
                    !
                    Jrot    = fitting%obs(iener)%Jrot
                    irot    = fitting%obs(iener)%irot
                    istate  = fitting%obs(iener)%quanta%istate
                    sigmai  = fitting%obs(iener)%quanta%sigma
                    ilambda = fitting%obs(iener)%quanta%ilambda
                    omega   = fitting%obs(iener)%quanta%omega
                    spin    = fitting%obs(iener)%quanta%spin
                    ivib    = fitting%obs(iener)%quanta%v
                    iTAG    = fitting%obs(iener)%quanta%iTAG(1:3)
                    !
                    Jrot_   = fitting%obs(iener)%Jrot_
                    irot_   = fitting%obs(iener)%irot_
                    istate_ = fitting%obs(iener)%quanta_%istate
                    sigmai_ = fitting%obs(iener)%quanta_%sigma
                    ilambda_= fitting%obs(iener)%quanta_%ilambda
                    omega_  = fitting%obs(iener)%quanta_%omega
                    spin_   = fitting%obs(iener)%quanta_%spin
                    ivib_   = fitting%obs(iener)%quanta_%v
                    iTAG_   = fitting%obs(iener)%quanta_%iTAG(1:3)
                    !
                    itau  = fitting%obs(iener)%iparity +1 ;  if (itau <1.or.itau >2) cycle 
                    itau_ = fitting%obs(iener)%iparity_+1 ;  if (itau_<1.or.itau_>2) cycle 
                    !
                    i  = fitting%obs(iener)%N ;  if (i >nenergies(irot ,itau )) cycle 
                    i_ = fitting%obs(iener)%N_;  if (i_>nenergies(irot_,itau_)) cycle 
                    !
                    ! compare with theoretical frequencies around the main transition
                    !
                    mark_ = ' '
                    !
                    do k = max(i-5,1),min(i+5,nenergies(irot,itau))
                       !
                       do k_ = max(i_-3,1),min(i_+3,nenergies(irot_,itau_))
                          !
                          if (abs(fitting%obs(iener)%energy-(energy_(irot,itau,k)-energy_(irot_,itau_,k_)))& 
                              >real(maxiter_as,rk)*abs(fitting%threshold_lock)) cycle
                          !
                          if (i/=k.or.i_/=k_) mark_ = '*'
                          !
                          JrotC    = calc(irot,itau,k)%Jrot
                          istateC  = calc(irot,itau,k)%istate
                          sigmaC   = calc(irot,itau,k)%sigma
                          ilambdaC = calc(irot,itau,k)%ilambda
                          omegaC   = sigmaC+real(ilambdaC,rk)
                          spinC    = calc(irot,itau,k)%spin
                          ivibC    = calc(irot,itau,k)%v
                          iTAGC    = calc(irot,itau,i)%iTAG(1:3)
                          !
                          JrotC_   = calc(irot_,itau_,k_)%Jrot
                          istateC_ = calc(irot_,itau_,k_)%istate
                          sigmaC_  = calc(irot_,itau_,k_)%sigma
                          ilambdaC_= calc(irot_,itau_,k_)%ilambda
                          omegaC_  = sigmaC_+real(ilambdaC_,rk)
                          spinC_   = calc(irot_,itau_,k_)%spin
                          ivibC_   = calc(irot_,itau_,k_)%v
                          iTAGC_   = calc(irot_,itau_,i_)%iTAG
                          !
                          write (frequnit,"(i5,2(i5,1x,f7.1,1x,a1),2x,' ',3f14.4,2x,e9.2,2x,2f14.4,2x,&
                                 &'(',1x,a3,1x,2i4,2f7.1,' <-',1x,a3,1x,2i4,2f7.1,' )" &
                                 // "(',1x,a3,1x,2i4,2f7.1,' <-',1x,a3,1x,2i4,2f7.1,' )',a)") &
                                 iener,k,Jrot,pm(itau),k_,Jrot_,pm(itau_),&
                                 fitting%obs(iener)%energy,energy_(irot,itau,k)-energy_(irot_,itau_,k_),&
                                 fitting%obs(iener)%energy-(energy_(irot,itau,k)-energy_(irot_,itau_,k_)),&
                                 wtall(iener),energy_(irot,itau,k)-ezero(1),energy_(irot_,itau_,k_)-ezero(1),&
                                 iTAGC,ivibC,ilambdaC,sigmaC,omegaC,iTAGC_,ivibC_,ilambdaC_,sigmaC_,omegaC_,&
                                 iTAG ,ivib ,ilambda ,sigmai,omega ,iTAG_,ivib_ ,ilambda_ ,sigmai_,omega_ ,&
                                 mark_
                       enddo
                       !
                    enddo
                    !
                 enddo ! --- i
                 !
              endif
              !
              ! print all energies
              !
              if (do_print) then
                 !
                 do irot = 1,nrot
                    !
                    do itau=1,2
                        !
                        do i = 1,nenergies(irot,itau)
                           !
                           if (i==1.or.energy_(irot,itau,i)-ezero(1)>sqrt(small_)) then 
                              !
                              printed = .false.
                              !
                              do jener=1,en_npts
                                !
                                ! lower statet must be the gound state
                                !
                                if (fitting%obs(jener)%N ==i .and.fitting%obs(jener)%irot ==irot & 
                                    .and.fitting%obs(jener)%iparity+1==itau) then 
                                   !
                                   printed = .true.
                                   !
                                   Jrot   = fitting%obs(jener)%Jrot
                                   istate = fitting%obs(jener)%quanta%istate
                                   sigmai = fitting%obs(jener)%quanta%sigma
                                   ilambda= fitting%obs(jener)%quanta%ilambda
                                   omega  = fitting%obs(jener)%quanta%omega
                                   spin   = fitting%obs(jener)%quanta%spin
                                   ivib   = fitting%obs(jener)%quanta%v
                                   iTAG    = fitting%obs(jener)%quanta%iTAG(1:3)
                                   !
                                   Jrot_   = fitting%obs(jener)%Jrot_
                                   istate_ = fitting%obs(jener)%quanta_%istate
                                   sigmai_ = fitting%obs(jener)%quanta_%sigma
                                   ilambda_= fitting%obs(jener)%quanta_%ilambda
                                   omega_  = fitting%obs(jener)%quanta_%omega
                                   spin_   = fitting%obs(jener)%quanta_%spin
                                   ivib_   = fitting%obs(jener)%quanta_%v
                                   iTAG_   = fitting%obs(jener)%quanta_%iTAG(1:3)
                                   !
                                   i_ = fitting%obs(jener)%N_
                                   irot_ = fitting%obs(jener)%irot_
                                   itau_ = fitting%obs(jener)%iparity_+1
                                   !
                                   write(enunit,"(2i5,1x,f7.1,1x,a1,2x,' ',3f14.4,2x,e9.2,2x,'(',1x,a3,1x,2i4,2f7.1,' )', &
                                         &'(',1x,a3,1x,2i4,2f7.1,' )',a1,2x,a)") & 
                                      i,fitting%obs(jener)%N,Jrot,pm(itau),&
                                      enercalc(jener)+eps(jener)+energy_(irot_,itau_,i_)-ezero(1),&
                                      energy_(irot,itau,i)-ezero(1),&
                                      energy_(irot,itau,i)-(enercalc(jener)+eps(jener)+energy_(irot_,itau_,i_)),wtall(jener),&
                                      iTAG, ivib, ilambda ,sigmai ,omega ,&
                                      iTAG_,ivib_,ilambda_,sigmai_,omega_,&
                                      mark(jener),trim( poten(istate)%name )
                                   !
                                endif 
                                !
                              enddo
                              !
                              if (.not.printed) then 
                                   !
                                   Jrot     = calc(irot,itau,i)%Jrot
                                   istate_  = calc(irot,itau,i)%istate
                                   sigmai_  = calc(irot,itau,i)%sigma
                                   ilambda_ = calc(irot,itau,i)%ilambda
                                   omega_   = calc(irot,itau,i)%omega
                                   spin_    = calc(irot,itau,i)%spin
                                   ivib_    = calc(irot,itau,i)%v
                                   iTAG_    = calc(irot,itau,i)%iTAG(1:3)
                                   !
                                   write(enunit,"(2i5,1x,f7.1,1x,a1,2x,' ',3f14.4,2x,e9.2,2x,'(',1x,a3,1x,2i4,2f7.1,' )',2x,a)") &
                                      i,0,Jrot,pm(itau),0.0,energy_(irot,itau,i)-energy_(1,1,1),0.0,0.0,&
                                      iTAG_,ivib_,ilambda_,sigmai_,omega_,trim( poten(istate_)%name )
                                   !
                              endif
                              !
                           endif
                           !
                        enddo
                        !
                    enddo
                    !
                 enddo
                 !
              endif
              !
            else
              !
              if (abs(fitting%threshold_lock)>0) then 
                !
                k0 = 1
                !
                if (fitting%threshold_lock<0)  maxiter_as = 1
                !
                jrot_ = fitting%obs(1)%jrot
                itau_ = fitting%obs(1)%iparity+1
                !
                do iener = 1,en_npts
                  !
                  mark(iener) = '*'
                  !
                  if (jrot_/=fitting%obs(iener)%jrot.or.itau_ /= fitting%obs(iener)%iparity+1) k0 = 1
                  !
                  jrot_ = fitting%obs(iener)%jrot
                  itau_ = fitting%obs(iener)%iparity+1
                  !
                  loop_thresh : do iter_th = 1,maxiter_as
                    !
                    loop_jrot : do irot = 1,nrot
                      !
                      jrot = fitting%J_list(irot)
                      !
                      if (nint(jrot-fitting%obs(iener)%jrot)/=0) cycle
                      !
                      do itau = 1,2
                        !
                        if (itau-1/=fitting%obs(iener)%iparity) cycle
                        !
                        loop_i_l : do i0 = k0,nenergies(irot,itau)
                          !
                          do jener = max(iener-10,1),iener-1
                            !
                            if (nint(jrot-fitting%obs(jener)%jrot)/=0.or.itau-1/=fitting%obs(jener)%iparity) cycle
                            !
                            if (fitting%obs(jener)%N == i0) cycle loop_i_l
                            !
                          enddo
                          !
                          if (i0>size(energy_,dim=3)) then
                            write(out,"('sf_fitting error: the size of energy_ is too small ',i8,';')") size(energy_,dim=3)
                            write(out,"('                  J= ',f9.2,',  and state =  ',i8,'; increase vmax (',i8,')')") &
                                                                                                            jrot,i0,vmax
                            stop 'error: the size of energy_ is too small'
                          endif
                          !
                          if (abs( fitting%obs(iener)%energy-(energy_(irot,itau,i0)-ezero(1)) )<= & 
                                                      & real(iter_th,rk)*abs(fitting%threshold_lock).and.&
                                ! threshold_lock>0, QN match
                              ( ( abs(calc(irot,itau,i0)%ilambda)==abs(fitting%obs(iener)%quanta%ilambda).and.&
                                  calc(irot,itau,i0)%istate ==fitting%obs(iener)%quanta%istate.and.&
                                  calc(irot,itau,i0)%v      ==fitting%obs(iener)%quanta%v.and.&
                                  nint(abs(calc(irot,itau,i0)%omega)-abs(fitting%obs(iener)%quanta%omega) )==0.and.&
                                  nint(abs(calc(irot,itau,i0)%sigma)-abs(fitting%obs(iener)%quanta%sigma) )==0 ).or.& 
                                 ! energy match + state + vib QN match 
                                 ( fitting%threshold_lock<0.and. &
                                  calc(irot,itau,i0)%istate ==fitting%obs(iener)%quanta%istate.and.&
                                  calc(irot,itau,i0)%v      ==fitting%obs(iener)%quanta%v ) ) ) then 
                              !
                              fitting%obs(iener)%N = i0
                              mark(iener) = ' '
                              exit loop_thresh
                              !
                          endif 
                          !
                        enddo loop_i_l
                        !
                      enddo
                      !
                    enddo loop_jrot
                    !
                  enddo loop_thresh
                  !
                  if (mark(iener) == ' ')  k0 = max(fitting%obs(iener)%N-10,1)
                  !
                enddo
                !
              endif
              !
              ! energy fit
              !
              do iener = 1,en_npts
                 !
                 i  = fitting%obs(iener)%N
                 jrot  = fitting%obs(iener)%jrot 
                 irot  = fitting%obs(iener)%irot 
                 !
                 itau = fitting%obs(iener)%iparity +1 ;  if (itau <1.or.itau >2) cycle
                 !
                 if (i >nenergies(irot,itau )) cycle
                 !
                 enercalc(iener) = energy_(irot,itau,i)-ezero(1)
                 !
                 eps(iener) = fitting%obs(iener)%energy-enercalc(iener)
                 !
                 if (do_deriv) then  
                    rjacob(iener,1:numpar) =  deriv(irot,itau,i,1:numpar) - deriv(1,1,1,1:numpar)
                 endif
                 !
                 if (do_print) then 
                   !
                   istate = fitting%obs(iener)%quanta%istate
                   sigmai = fitting%obs(iener)%quanta%sigma
                   ilambda = fitting%obs(iener)%quanta%ilambda
                   omega   = fitting%obs(iener)%quanta%omega
                   spin = fitting%obs(iener)%quanta%spin
                   ivib = fitting%obs(iener)%quanta%v
                   iTAG = fitting%obs(iener)%quanta%iTAG(1:3)
                   !
                   JrotC   = calc(irot,itau,i)%Jrot
                   istateC = calc(irot,itau,i)%istate
                   sigmaC = calc(irot,itau,i)%sigma
                   ilambdaC= calc(irot,itau,i)%ilambda
                   omegaC  = sigmaC + real(ilambdaC,rk)
                   spinC   = calc(irot,itau,i)%spin
                   ivibC   = calc(irot,itau,i)%v
                   iTAGC = calc(irot,itau,i)%iTAG(1:3)
                   !
                   write (out,"(2i5,1x,f8.1,1x,a1,2x,' ',3f14.4,2x,e9.2,2x,&
                          &'(',1x,a3,1x,2i4,2f8.1,' )','(',1x,a3,1x,2i4,2f8.1,' )',a)") &
                          iener,i,Jrot,pm(itau),enercalc(iener)+eps(iener),enercalc(iener),eps(iener),&
                          wtall(iener),&
                          iTAGC,ivibC,ilambdaC,sigmaC,omegaC,&
                          iTAG ,ivib ,ilambda ,sigmai,omega ,&
                          mark(iener)
                 endif
              enddo
              !
              ! Printing all calculated term values. If the obs. counterpats exist, 
              ! the obs.-calc. are printed as well. 
              ! This list can be used to identify the obs-calc pairs, i.e. the number of 
              ! a term value as it appear in the list of calculated energies. 
              !
              do irot = 1,nrot
                 !
                 jrot = fitting%J_list(irot)
                 !
                 do itau = 1,2
                   !
                   do i = 1,nenergies(irot,itau)
                     !
                     if (i==1.or.energy_(irot,itau,i)-ezero(1)>sqrt(small_)) then 
                        !
                        jener = 1
                        !
                        do while ( jener/=en_npts+1)
                          !
                          if (fitting%obs(jener)%N==i.and.fitting%obs(jener)%irot==irot.and. & 
                              fitting%obs(jener)%iparity+1==itau) exit 
                          !
                          jener = jener+1
                          !
                        enddo
                        !
                        if (do_print) then
                          ! 
                          if (jener<=en_npts) then 
                             !
                             Jrot   = fitting%obs(jener)%Jrot
                             istate = fitting%obs(jener)%quanta%istate
                             sigmai = fitting%obs(jener)%quanta%sigma
                             ilambda= fitting%obs(jener)%quanta%ilambda
                             omega  = fitting%obs(jener)%quanta%omega
                             spin   = fitting%obs(jener)%quanta%spin
                             ivib   = fitting%obs(jener)%quanta%v
                             iTAG   = fitting%obs(jener)%quanta%iTAG(1:3)
                             !
                             istateC = calc(irot,itau,i)%istate
                             sigmaC  = calc(irot,itau,i)%sigma
                             ilambdaC= calc(irot,itau,i)%ilambda
                             omegaC  = sigmaC + real(ilambdaC,rk)
                             spinC   = calc(irot,itau,i)%spin
                             ivibC   = calc(irot,itau,i)%v
                             iTAGC   = calc(irot,itau,i)%iTAG(1:3)
                             !
                             write(enunit,"(2i5,1x,f8.1,1x,a1,2x,' ',3f14.4,2x,e9.2,2x,'(',1x,a3,1x,2i4,2f8.1,' )', &
                                          & '(',1x,a3,1x,2i4,2f8.1,' )',a)") &
                                i,fitting%obs(jener)%N,Jrot,pm(itau),&
                                enercalc(jener)+eps(jener),&
                                enercalc(jener),eps(jener),wtall(jener),&
                                iTAGC,ivibC,ilambdaC,sigmaC,omegaC,&
                                iTAG ,ivib  ,ilambda,sigmai,omega,&
                                mark(jener)
                            !
                          else
                            !
                            Jrot = calc(irot,itau,i)%Jrot
                            istate_ = calc(irot,itau,i)%istate
                            sigmai_ = calc(irot,itau,i)%sigma
                            ilambda_ = calc(irot,itau,i)%ilambda
                            omega_ = calc(irot,itau,i)%omega
                            spin_ = calc(irot,itau,i)%spin
                            ivib_ = calc(irot,itau,i)%v
                            iTAG_ = calc(irot,itau,i)%iTAG(1:3)
                            !
                            write(enunit,"(2i5,1x,f8.1,1x,a1,2x,' ',3f14.4,2x,e9.2,2x,'(',1x,a3,1x,2i4,2f8.1,' )')") &
                               i,0,Jrot,pm(itau),0.0,energy_(irot,itau,i)-ezero(1),0.0,0.0,&
                               iTAG_,ivib_,ilambda_,sigmai_,omega_
                               !
                          endif 
                          !
                        endif 
                        !
                     endif 
                     !
                   enddo
                   !
                 enddo
                 !
              enddo
              !
            endif 
            !
            ! switch off energies with large obs-calc assuming unwanted swapping 
            if (fitting%threshold_obs_calc>small_) then
               do iener = 1,en_npts
                  if (abs(eps(iener))>fitting%threshold_obs_calc.and.wtall(iener)>small_) then
                    wtall(iener) = 0
                    if (fitting%robust>small_) sigma(iener) = 0
                    eps(iener) = 0
                    nused = nused - 1
                  endif
               enddo
               nused = max(nused,1)
               wtsum = sum(wtall(1:npts))
               wtall(1:npts) = wtall(1:npts)/wtsum
            endif 
            !
            ! Here the potential energy section starts. 
            !
            call TimerStart('Grid points')
            !
            rewind(abinitunit)
            !
            write(abinitunit,"(' Grid points '/&
            &'  state        r                 ab initio             calc.           ab initio - calc        weight'/)")
            !
            ! this will count different ab initio fields one by one 
            ifield_ = 0
            !
            ! this will count different ab initio grid points one by one going from one field to another field 
            ! (remember that different ab initio fields might have different grids)
            j = 0
            !
            do iobject =1,Nobjectmax
               !
               Nfields = fieldmap(iobject)%Nfields
               !
               do ifield =1,Nfields
                 !
                 ifield_  = objects(iobject,ifield)%field%iabi
                 !
                 if (ifield_==0) cycle 
                 !
                 field => objects(iobject,ifield)%field
                 !
                 if (field%type/="DUMMY") write(abinitunit,'(a)') field%name
                 !
                 if (associated(field%sub_type)) call define_complex_analytic_field_subterms(field%type,&
                                                      field%sub_type,field%Nsub_terms)
                 !
                 do i = 1,abinitio(ifield_)%Nterms
                   !
                   j = j + 1
                   !
                   r_t = abinitio(ifield_)%grid(i)
                   !
                   ! dummy field 
                   !
                   if (objects(iobject,ifield)%field%type=='GRID') then
                     !
                     f = objects(iobject,ifield)%field%gridvalue(ipotmin)
                     !
                   else
                     !
                     f = objects(iobject,ifield)%field%fanalytic_field(r_t,objects(iobject,ifield)%field%value)
                     !
                   endif
                   !
                   f = objects(iobject,ifield)%field%factor*f
                   !
                   if (objects(iobject,ifield)%field%morphing) then 
                       f = f*abinitio(ifield_)%value(i)*abinitio(ifield_)%factor
                   endif 
                   !
                   ! eps - epsilon = ab initio energies - calculated pot. energies,
                   ! where we comntinue counting the fitting data points starting with en_npts - 
                   ! obs. data. We also apply a shift to the constraining field if necessary
                   !
                   eps(en_npts+j) = abinitio(ifield_)%value(i)*abinitio(ifield_)%factor-f
                   !
                   if (abs(abinitio(ifield_)%refvalue)>small_) then  
                      fshift = objects(iobject,ifield)%field%value(1)-abinitio(ifield_)%refvalue
                      eps(en_npts+j) = eps(en_npts+j)+fshift
                   endif
                   !
                   ! print out into .pot
                   !
                   write (abinitunit,"((2x,i4,2x,f16.9),3(2x,g19.8),2x,e12.4)") & 
                          ifield,r_t, &
                          abinitio(ifield_)%value(i)*abinitio(ifield_)%factor,f, &
                          eps(en_npts+j),wtall(en_npts+j) 
                   !
                   ! calculate the derivative of the potential function wrt the fitting parameters
                   !
                   if (itmax.ge.1.and.(.not.do_Armijo.or.ialpha_Armijo==0)) then 
                      !
                      select case (trim(objects(iobject,ifield)%field%type)) 
                        !
                      case('GRID')
                        !
                        allocate(spline_wk_vec(field%Nterms),stat=info)
                        call ArrayStart('spline_wk_vec-fit',info,ngrid,kind(spline_wk_vec))
                        !
                        yp1= 0._rk ; ypn =0._rk  ! 1nd derivatives at the first and last point (ignored)
                        !
                        do ifitpar = 1,numpar
                          !
                          rjacob(en_npts+j,ifitpar) = 0
                          !
                          jfield = fit_index(ifitpar)%ifield
                          iterm = fit_index(ifitpar)%iterm
                          !
                          ! when the constraint applied to the reference field the derivative wrt the first parameter is zero
                          if (abs(abinitio(ifield_)%refvalue)>small_.and.iterm==1) cycle
                          !
                          if (ifield/=jfield.or.iobject/=fit_index(ifitpar)%iobject)cycle 
                          !
                          param_t = field%value(iterm)
                          !
                          delta = 1e-3*abs(param_t)
                          !
                          if (delta < small_) delta = sqrt(small_)*10.0
                          !
                          field%value(iterm) = param_t + delta
                          !
                          call spline(field%grid,field%value,field%Nterms,yp1,ypn,spline_wk_vec)
                          !
                          ! evaluate spline interpolant
                          call splint(field%grid,field%value,spline_wk_vec,field%Nterms,r_t,FR)
                          !
                          field%value(iterm) = param_t - delta
                          !
                          call spline(field%grid,field%value,field%Nterms,yp1,ypn,spline_wk_vec)
                          !
                          ! evaluate spline interpolant
                          call splint(field%grid,field%value,spline_wk_vec,field%Nterms,r_t,FL)
                          !
                          rjacob(en_npts+j,ifitpar) = (fR-fL)/(2.0_rk*delta)*field%factor
                          !
                          field%value(iterm) = param_t
                          !
                        enddo
                        !
                        deallocate(spline_wk_vec)
                        call ArrayStop('spline_wk_vec-fit')
                        !
                      case default 
                        !
                        do ifitpar = 1,numpar
                          !
                          rjacob(en_npts+j,ifitpar) = 0
                          !
                          !if (objects(iobject,ifield)%field%type=='GRID') cycle
                          !
                          jfield = fit_index(ifitpar)%ifield
                          iterm = fit_index(ifitpar)%iterm
                          !
                          ! when the constraint applied to the reference field the derivative wrt the first parameter is zero
                          if (abs(abinitio(ifield_)%refvalue)>small_.and.iterm==1) cycle
                          !
                          if (ifield/=jfield.or.iobject/=fit_index(ifitpar)%iobject)cycle 
                          !
                          param_t = objects(iobject,ifield)%field%value(iterm)
                          !
                          delta = 1e-4*abs(param_t)
                          !
                          if (delta < small_) delta = sqrt(small_)*10.0
                          !
                          objects(iobject,ifield)%field%value(iterm) = param_t + delta
                          !
                          fR =  objects(iobject,ifield)%field%fanalytic_field(r_t,objects(iobject,ifield)%field%value)
                          !
                          objects(iobject,ifield)%field%value(iterm) = param_t - delta
                          !
                          fL =  objects(iobject,ifield)%field%fanalytic_field(r_t,objects(iobject,ifield)%field%value)
                          !
                          rjacob(en_npts+j,ifitpar) = (fR-fL)/(2.0_rk*delta)*objects(iobject,ifield)%field%factor
                          !
                          objects(iobject,ifield)%field%value(iterm) = param_t
                          !
                        enddo
                        !
                      end select 
                      !
                   endif 
                   !
                 enddo
                 !
               enddo
            enddo
            !
            call TimerStop('Grid points')
            !
            ! ssq  - weighted rms**2, rms  - root mean square deviation. 
            !
            ssq=sum(eps(1:npts)*eps(1:npts)*wtall(1:npts))
            rms=sqrt(sum(eps(1:npts)*eps(1:npts))/npts)
            !
            if (do_print.and.fititer==1) write(out,"(/a)") 'Refinement using the least-squared fitting ...'
            if (do_print) write(out,"(/'Iteration = ',i8)") fititer
            !
            ! Prepare the linear system a x = b as in the Newton fitting approach.  
            !
            if (itmax>=1.and.(.not.do_Armijo.or.ialpha_Armijo==0)) then
               !
               !----- form the a and b matrix ------c
               ! form A matrix 
               do irow=1,numpar       
                 do icolumn=1,irow    
                   al(irow,icolumn)=sum(rjacob(1:npts,icolumn)*rjacob(1:npts,irow)*wtall(1:npts))
                   al(icolumn,irow)=al(irow,icolumn)
                   !dec if (fit_debug > 2)
                   !  write (out,"('al (',i,',',i,')= ',es14.7)") irow,icolumn,al(irow,icolumn)
                   !dec end if
                 enddo
               enddo
               !
               ! form B matrix 
               do irow=1,numpar      
                 bl(irow)=sum(eps(1:npts)*rjacob(1:npts,irow)*wtall(1:npts))
                 !dec if (fit_debug > 2)
                 !  write (out,"('bl (',i,')= ',es14.7)") irow,bl(irow)
                 !dec end if
               enddo  
               !
               ! In case any diagonal matrix elements of al is zero, we can have remove this paramter 
               ! from the fit and set its value to zero. And start the iteration again. 
               do ncol=1,numpar 
                  !
                  if ( abs(al(ncol,ncol))<small_ ) then 
                      !
                      iobject = fit_index(ncol)%iobject
                      ifield = fit_index(ncol)%ifield
                      iterm = fit_index(ncol)%iterm
                      !
                      objects(iobject,ifield)%field%weight(iterm) = 0
                      objects(iobject,ifield)%field%value(iterm) = object0(iobject,ifield)%value(iterm)
                      !
                      write(out,"(i0,'-th is out - ',a8)") i,objects(iobject,ifield)%field%forcename(iterm)
                      !
                      cycle outer_loop
                      !
                  endif 
                  !
               enddo 
               !
               ! We use two different approaches to Least squares fit Ax=b
               ! 1. A(ij) = sum_k R(i,k) R(j,k), where R(i,k) = d Ei/d Ck and b(i) = sum_k R(i,k) [Eobs(i)-Ecalc(i)]
               ! 2. A(i,j) = R(i,j) and b(i) = [Eobs(i)-Ecalc(i)] <= SVD
               !
               select case (trim(fitting%fit_type)) 
               !
               case default
                 !
                 write (out,"('fit_type ',a,' unknown')") trim(fitting%fit_type)
                 stop 'fit_type unknown'
                 !
               case('LINUR','DGELSS') 
                  !
                  ! Using Marquardt's fitting method
                  !
                  ! solve the linear equatins for two values of lambda and lambda/10
                  !
                  ! Defining scaled (with covariance) A and b
                  ! 
                  ! form A matrix 
                  do irow=1,numpar       
                    do icolumn=1,irow    
                      Am(irow,icolumn)=al(irow,icolumn)/sqrt( al(irow,irow)*al(icolumn,icolumn) )
                      Am(icolumn,irow)=Am(irow,icolumn)
                    enddo
                    bm(irow) = bl(irow)/sqrt(al(irow,irow))
                  enddo
                  !
                  ! define shifted A as A =  A+lambda I
                  ! lambda is Marquard's scaling factor
                  !
                  do irow=1,numpar       
                      Am(irow,irow)=Am(irow,irow)*(1.0_rk+lambda)
                  enddo
                  !
                  ! Two types of the linear solver are availible: 
                  ! 1. linur (integrated into the program, from Ulenikov Oleg)
                  ! 2. dgelss - Lapack routine (recommended).
                  !
                  select case (trim(fitting%fit_type)) 
                  !
                  case default
                    !
                    write (out,"('fit_type ',a,' unknown')") trim(fitting%fit_type)
                    stop 'fit_type unknown'
                    !
                  case('LINUR') 
                    !
                    call MLlinur(numpar,numpar,am(1:numpar,1:numpar),bm(1:numpar),solution(1:numpar),ierror)
                    !
                    ! In case of dependent parameters  "linur" reports an error = ierror, 
                    ! which is a number of the dependent parameter. We can remove this paramter 
                    ! from the fit and set its value to zero. And start the iteration again. 
                    !
                    if (ierror.ne.0) then 
                      do ncol=1,numpar 
                         !
                         if ( ncol.eq.ierror ) then 
                             !
                             iobject = fit_index(ncol)%iobject
                             ifield = fit_index(ncol)%ifield
                             iterm = fit_index(ncol)%iterm
                             !
                             objects(iobject,ifield)%field%weight(iterm) = 0
                             objects(iobject,ifield)%field%value(iterm) = object0(iobject,ifield)%value(iterm)
                             !
                             write(out,"(i0,'-th is out - ',a8)") i,objects(iobject,ifield)%field%forcename(iterm)
                             !
                         endif 
                         !
                      enddo 
                      !
                      cycle outer_loop
                     endif 
                     !
                  case ('DGELSS')
                    !
                    ai = am 
                    bi = bm
                    call dgelss(numpar,numpar,1,ai(1:numpar,1:numpar),numpar,bi(1:numpar),numpar,Tsing,-1.D-12,rank,&
                                wspace,lwork,info)
                    !
                    if (info/=0) then
                      write(out,"('dgelss:error',i0)") info
                      stop 'dgelss'
                    endif
                    !
                    solution = bi ! *0.1
                    !
                 end select 
                 !
                 ! convert back from Marquardt's representation
                 !
                 do ncol=1,numpar
                    solution(ncol) =  solution(ncol)/sqrt(al(ncol,ncol))
                 enddo
                 !
               case ("SDD")
                 !
                 call lapack_sdd_pseudo_inverse(fitting%svd_tol,rjacob(1:npts,1:numpar),NumSVD)
                 !
                 solution(1:numpar) = matmul(eps(1:npts)*wtall(1:npts),rjacob(1:npts,1:numpar))
                 !
                 if (verbose>=4) write(out,"(/a,1x,i8,1x,a,g12.5)") 'Number of SDD roots = ',NumSVD,'tol = ',fitting%svd_tol
                 !
               case ("SVD")
                 !
                 call lapack_svd_pseudo_inverse(fitting%svd_tol,rjacob(1:npts,1:numpar),NumSVD)
                 !
                 solution(1:numpar) = matmul(eps(1:npts)*wtall(1:npts),rjacob(1:npts,1:numpar))
                 !
                 if (verbose>=4) write(out,"(/a,1x,i8,1x,a,g12.5)") 'Number of SVD roots = ',NumSVD,'tol = ',fitting%svd_tol
                 !
               end select
               !
               ! Scale the correction if required 
               !
               dx = solution*fitting%fit_scaling
               !
               !----- update the parameter values with a scaled correction------!
               !
               !$omp parallel do private(ncol,iobject,ifield,iterm,param_t) shared(params_t) schedule(dynamic)
               do ncol=1,numpar
                  !
                  iobject = fit_index(ncol)%iobject
                  ifield = fit_index(ncol)%ifield
                  iterm = fit_index(ncol)%iterm
                  !
                  params_t(ncol) = objects(iobject,ifield)%field%value(iterm)
                  !
                  param_t = params_t(ncol)+dx(ncol)
                  param_t = max(objects(iobject,ifield)%field%fit_range(iterm)%min,param_t)
                  param_t = min(objects(iobject,ifield)%field%fit_range(iterm)%max,param_t)
                  !
                  objects(iobject,ifield)%field%value(iterm)=param_t
                  !
               enddo
               !$omp end parallel do
               !
               !
               ! Robust fit: adjust the fitting weights
               !
               if (fitting%robust>small_) then
                 !
                 call robust_fit(a_wats,sigma(1:npts),eps(1:npts),wtall(1:npts))
                 !
                 ssq=sum(eps(1:npts)*eps(1:npts)*wtall(1:npts))
                 !
                 !fitting%robust = -1.0_rk
                 !
               endif 
               !
               ! Estimate standard deviation error. 
               !
               if ( nused.ne.numpar ) then 
                 stadev=sqrt(ssq/float(nused-numpar))
               else 
                 stadev=sqrt(ssq/nused)
               endif
               !               
               if (stadev<stadev_old) then
                 lambda = lambda/nu
               else 
                 lambda = min(lambda*nu,10000.0_rk)
               endif
               !
               ! Estimate the standard errors for each parameter using 
               ! the inverse matrix of a. 
               !
               call MLinvmat(al(1:numpar,1:numpar),ai(1:numpar,1:numpar),numpar,info)
               !
               sum_sterr=0.d0
               do ncol=1,numpar 
                  if (nused.eq.numpar) then  
                     sterr(ncol)=0
                  else
                     !
                     iobject = fit_index(ncol)%iobject
                     ifield = fit_index(ncol)%ifield
                     iterm = fit_index(ncol)%iterm
                     !                     
                     sterr(ncol)=sqrt(abs(ai(ncol,ncol)))*stadev
                     sum_sterr=sum_sterr+abs(sterr(ncol)/objects(iobject,ifield)%field%value(iterm))
                     !
                  endif
               enddo    
               !
               sum_sterr=sum_sterr/numpar 
               !
               ! This is how we define stability of the fit:
               ! as a relative change of stadev comparing with the step before. 
               !
               stability=abs( (stadev-stadev_old)/stadev )
               stadev_old=stadev
               !
               if(do_Armijo) then
                 !
                 alpha_Armijo = (1.0_ark-1.0_rk/real(NArmijo)*ialpha_Armijo)
                 dx = solution*alpha_Armijo
                 !
                 rms0=sum(eps(1:npts)*eps(1:npts)*wtall(1:npts))*0.5_rk
                 !
                 !----- update the parameter values with a scaled correction------!
                 !
                 !$omp parallel do private(ncol,iobject,ifield,iterm,param_t) schedule(dynamic)
                 do ncol=1,numpar
                    !
                    iobject = fit_index(ncol)%iobject
                    ifield = fit_index(ncol)%ifield
                    iterm = fit_index(ncol)%iterm
                    !
                    param_t = params_t(ncol)+dx(ncol)
                    param_t = max(objects(iobject,ifield)%field%fit_range(iterm)%min,param_t)
                    param_t = min(objects(iobject,ifield)%field%fit_range(iterm)%max,param_t)
                    !
                    objects(iobject,ifield)%field%value(iterm)=param_t
                    !
                 enddo
                 !$omp end parallel do 
                 !
                 fititer = fititer - 1
                 ialpha_Armijo = ialpha_Armijo + 1
                 !
                 do_print = .false.
                 !
                 cycle loop_fititer
                 !
               endif
               !
            else
               !
               if(do_Armijo.and.itmax>=1) then
                 !
                 ! Armijo condtion: rms2 must be < rms1 
                 !
                 rms1 = rms0 + alpha_Armijo*sum(bl(1:numpar)*solution(1:numpar))
                 !
                 rms2=sum(eps(1:npts)*eps(1:npts)*wtall(1:npts))*0.5_rk
                 !
                 alpha_Armijo = (1.0_rk - 1.0_rk/real(NArmijo)*ialpha_Armijo)
                 !
                 ialpha_Armijo = ialpha_Armijo + 1
                 !
                 dx = solution*alpha_Armijo
                 !
                 ! continue linear search if the condition is not fullfilled or too many iterations 
                 if (rms2>=rms1.and.ialpha_Armijo<NArmijo) then 
                    !
                    fititer = fititer - 1
                    !
                    !$omp parallel do private(ncol,iobject,ifield,iterm,param_t) schedule(dynamic)
                    do ncol=1,numpar
                       !
                       iobject = fit_index(ncol)%iobject
                       ifield = fit_index(ncol)%ifield
                       iterm = fit_index(ncol)%iterm
                       !
                       param_t = params_t(ncol)+dx(ncol)
                       param_t = max(objects(iobject,ifield)%field%fit_range(iterm)%min,param_t)
                       param_t = min(objects(iobject,ifield)%field%fit_range(iterm)%max,param_t)
                       !
                       objects(iobject,ifield)%field%value(iterm)=param_t
                       !
                    enddo
                    !$omp end parallel do 
                    !
                    cycle loop_fititer
                 endif
                 !
                 !----- update the parameter values with a scaled correction------!
                 !
                 if (verbose>=4.and.do_Armijo) write(out,"(/a,f15.8,1x,i5,' steps')") 'Armijo    alpha  parameter = ',&
                                               alpha_Armijo,ialpha_Armijo
                 !
                 ialpha_Armijo = 0
                 rms0 = 0 
                 do_print = .true.
                 !
               endif
               !
               stadev=sqrt(ssq/max(nused,1))
               !
            endif
            !
            ! Print the updated parameters. 
            !
            write(out,"(/'Parameters:')")
            !
            ifield_ = 0
            !
            do iobject =1,Nobjectmax
               !
               Nfields = fieldmap(iobject)%Nfields
               !
               do ifield =1,Nfields
                 !
                 ifield_ = ifield_ + 1
                 !
                 !if (abinitio(ifield_)%type=="DUMMY") cycle
                 !
                 nchar = len_trim(objects(iobject,ifield)%field%class)
                 write(my_fmt, '(A,I0,A)') "(/a",nchar,",4x,a,1x,a)"
                 write(out,my_fmt) trim(objects(iobject,ifield)%field%class),&
                                   trim(objects(iobject,ifield)%field%itag),&
                                   trim(objects(iobject,ifield)%field%jtag)
                 !
                 !write(out,'(/tl,a20,2i4)') adjustl(trim(objects(iobject,ifield)%field%class)), &
                 !   objects(iobject,ifield)%field%iref,objects(iobject,ifield)%field%jref
                 !
                 write(out,'(a4,4x,a)') "name",adjustl(trim(objects(iobject,ifield)%field%name))
                 write(out,'(a4,2x,a20)') "type",adjustl(trim(objects(iobject,ifield)%field%type))
                 !write(out,'("N",3x,i8)') objects(iobject,ifield)%field%Nterms
                 write(out,'(a6)') "Values"
                 !
                 do iterm = 1,objects(iobject,ifield)%field%Nterms
                   !
                   flink => objects(iobject,ifield)%field%link(iterm)
                   frange => objects(iobject,ifield)%field%fit_range(iterm)
                   !
                   mark_f = '' 
                   !
                   if (objects(iobject,ifield)%field%weight(iterm)>0) mark_f = 'fit'
                   !
                   if (flink%iobject/=0) then
                     !
                     write (out,"(a8,3x,es22.14,2x,a3,6x,'link',3(1x,i3))") objects(iobject,ifield)%field%forcename(iterm), & 
                           objects(iobject,ifield)%field%value(iterm),mark_f,&
                           flink%iobject,flink%ifield,flink%iparam
                     !
                   else
                     !
                     if (frange%set) then 
                       !
                       write (out,"(a8,3x,es22.14,2x,a3,2x,'range',2x,es15.7,',',es15.7)") &
                                 objects(iobject,ifield)%field%forcename(iterm),&
                                 objects(iobject,ifield)%field%value(iterm),mark_f,frange%min,frange%max
                       !
                     else
                       !
                       write (out,"(a8,3x,es22.14,2x,a3)") objects(iobject,ifield)%field%forcename(iterm), & 
                             objects(iobject,ifield)%field%value(iterm),mark_f
                           !
                     endif
                   endif
                   !
                 enddo
                 !
                 write(out,'(a3)') "end"
                 !
                 !
               enddo
            enddo
            !
            ! Output some statistics and results 
            !
            !  only if we are fitting:  
            !
            if (itmax.ne.0) then
              !               write (out,"(//'Potential parameters rounded in accord. with their standart errors'/)")
              ncol = 0 
              ifield_ = 0
              !
              write(out,"(/'Fitted parameters (rounded):')")
              !
              do iobject =1,Nobjectmax
                !
                Nfields = fieldmap(iobject)%Nfields
                !
                do ifield =1,Nfields
                  !
                  ifield_ = ifield_ + 1
                  !
                  !if (abinitio(ifield_)%type=="DUMMY") cycle
                  !
                  write(out,'(/)') 

                  nchar = len_trim(objects(iobject,ifield)%field%class)
                  write(my_fmt, '(A,I0,A)') "(/a",nchar,",4x,a,1x,a)"
                  write(out,my_fmt) trim(objects(iobject,ifield)%field%class),&
                                    trim(objects(iobject,ifield)%field%itag),&
                                    trim(objects(iobject,ifield)%field%jtag)
                  !
                  write(out,'(a4,4x,a)') "name",adjustl(trim(objects(iobject,ifield)%field%name))
                  write(out,'(a4,4x,a20)') "type",adjustl(trim(objects(iobject,ifield)%field%type))
                  write(out,'("N",3x,i8)') objects(iobject,ifield)%field%Nterms
                  write(out,'(a6)') "Values"
                  !
                  do iterm=1,objects(iobject,ifield)%field%Nterms
                    !
                    if (nint(objects(iobject,ifield)%field%weight(iterm)) > 0) then 
                      !
                      ncol = ncol + 1
                      !
                      ndigits = 0
                      conf_int = sterr(ncol)
                      do while (conf_int.le.10.0.and.ndigits<10)
                        ndigits = ndigits +1 
                        conf_int = conf_int*10
                      enddo
                      !
                      if (conf_int>1e8) conf_int = 0 
                      !
                      write(my_fmt, '(A,I0,A)') "(a8,1x,f8.1,2x,f22.",ndigits,"'(',i14,')',2x,'[',3i5,']')"
                      write (out,my_fmt) objects(iobject,ifield)%field%forcename(iterm),&
                                         objects(iobject,ifield)%field%weight(iterm),objects(iobject,ifield)%field%value(iterm),&
                                         nint(conf_int),iobject,ifield,iterm
                    else 
                        !
                        ndigits =2
                        if (abs(objects(iobject,ifield)%field%value(iterm))>small_) ndigits = 8
                        !
                        write(my_fmt, '(A,I0,A)') "(a8,1x,f8.1,2x,f22.",ndigits,",18x,'[',3i5,']')"
                        write (out,my_fmt) objects(iobject,ifield)%field%forcename(iterm), & 
                               objects(iobject,ifield)%field%weight(iterm),objects(iobject,ifield)%field%value(iterm),&
                               iobject,ifield,iterm
                        !
                    endif 
                  enddo
                enddo
              enddo
              !
            endif 
            !
            still_run = .false.
            !
            ! Print out the ssq for the rovib. energies and pot. data points separetely:
            !
            ssq1 = 0 ; ssq2 = 0 
            !
            wtsum = sum(wt_bit(1:en_npts))
            !
            if (wtsum>small_) ssq1 = sqrt( sum(eps(1:en_npts)**2*dble(wt_bit(1:en_npts)))/wtsum )
            !
            wtsum = sum(wt_bit(1+en_npts:npts))
            !
            if (wtsum>small_) ssq2 = sqrt( sum(eps(1+en_npts:npts)**2*dble(wt_bit(1+en_npts:npts)))/wtsum )
            !
            !rms1=sqrt(sum(eps(1:en_npts)**2)/en_npts)
            !rms2=sqrt(sum(eps(1+en_npts:npts)**2)/max(pot_npts,1))
            !
            write (out,6552) fititer,nused,numpar,stadev,ssq1,ssq2,stability
            !
            if (verbose>=4) write(out,"(/'Marquardt lambda parameter = ',g15.8)") lambda
            !
            if (verbose>=4) call TimerReport
            !
          enddo loop_fititer  ! --- fititer
          !
       enddo outer_loop
       !
       if (allocated(potparam)) deallocate(potparam)
       if (allocated(deriv)) deallocate(deriv )
       call ArrayStop('fit-ener-deriv')
       !
       if (allocated(params_t)) deallocate(params_t)
       call ArrayStop('params_t')
       !
       deallocate (nampar,ivar,ifitparam,al,ai,bl,bi,dx,sterr,Tsing,pot_terms,solution,am,bm)
       call ArrayStop('potparam-mat')
       call ArrayStop('potparam-dx')
       call ArrayStop('potparam-sterr')
       call ArrayStop('potparam-Tsing')
       call ArrayStop('a-b-mat')
       call ArrayStop('pot_terms')
       !
       deallocate (pot_values)
       call ArrayStop('pot_values-mat')
       !
       deallocate (rjacob,eps,wspace,energy_,calc,enercalc)
       call ArrayStop('rjacob')
       call ArrayStop('eps')
       call ArrayStop('wspace')
       !
       deallocate(object0)
       call ArrayStop('object0')
       !
       deallocate(nenergies)
       call ArrayStop("nenergies")
       !
       if (allocated(sigma)) then 
         deallocate(sigma)
         call ArrayStop("sigma")
       endif 
       !
       close(enunit,status="keep")
       close(abinitunit,status="keep")!
       if (action%frequency) close(frequnit,status="keep")
       !
       call TimerStop('Simultaneous Fitting')
       !
       !
6552   format(/3X,85('-')/'   |  Iter  | Points | Params |   w-rms(total)|',&
       '    rms_ener   |   rms_pot   | Stability |'/&
       3X,85('-')/,&
       '-->| ',i6,' | ',i6,' | ',i6,' |  ',ES12.5,' | ',ES12.5,'  |  ',&
            ES10.3,' |',ES10.3,' |',/3X,85('-')/)

!6553   format(/3X,85('-')/'   |  Iter  | Points | Params |   w-rms(total)|',&
!       '    rms_ener   |   rms_pot   | Stability |'/&
!       3X,85('-')/,&
!       '   | ',i6,' | ',i6,' | ',i6,' |  ',ES12.5,' | ',ES12.5,'  |  ',&
!            ES10.3,' |',ES10.3,' |',/3X,85('-')/)


  contains 



  subroutine parameters_copy_to_fields
    !
    integer(ik) :: ifield,Nterms,iterm
       !
       iterm = 0
       !
       ! potential functions
       do ifield = 1,Nestates
          Nterms = poten(ifield)%Nterms
          poten(ifield)%value(1:Nterms) = potparam(iterm+1:iterm+Nterms)
          iterm = iterm + Nterms
       enddo
       !
       ! spin-orbits
       do ifield = 1,Nspinorbits
          Nterms = spinorbit(ifield)%Nterms
          spinorbit(ifield)%value(1:Nterms) = potparam(iterm+1:iterm+Nterms)
          iterm = iterm + Nterms
       enddo
       !
       ! L**2
       do ifield = 1,NL2
          Nterms = L2(ifield)%Nterms
          L2(ifield)%value(1:Nterms) = potparam(iterm+1:iterm+Nterms)
          iterm = iterm + Nterms
       enddo
       !
       ! LxLy
       do ifield = 1,NLxLy
          Nterms = LxLy(ifield)%Nterms
          LxLy(ifield)%value(1:Nterms) = potparam(iterm+1:iterm+Nterms)
          iterm = iterm + Nterms
       enddo
    !
  end subroutine parameters_copy_to_fields


  subroutine update_linked_parameters
    !
    integer(ik) :: ifield,iterm,iobject,Nfields
    type(linkT),pointer :: flink
    !
    ! update linked parameters 
    !
    do iobject =1,Nobjectmax
       !
       Nfields = fieldmap(iobject)%Nfields
       !
       do ifield =1,Nfields
         !
         !ifield_ = ifield_ + 1
         !
         do iterm = 1,objects(iobject,ifield)%field%Nterms
           !
           flink => objects(iobject,ifield)%field%link(iterm)
           !
           if (flink%iobject/=0) then
             !
             objects(iobject,ifield)%field%value(iterm) = objects(flink%iobject,flink%ifield)%field%value(flink%iparam)
             !
           endif 
           !
         enddo
         !
       enddo
    enddo
    !
  end subroutine update_linked_parameters

  subroutine one_parameter_copy_to_fields(i,param)
    !
    integer(ik),intent(in) :: i
    real(rk),intent(in) :: param
    integer(ik) :: ifield,Nterms,iterm,jterm
       !
       iterm = 0
       !
       ! potential functions
       do ifield = 1,Nestates
          Nterms = poten(ifield)%Nterms
          do jterm = 1,Nterms
            if (i==iterm+jterm) then 
               poten(ifield)%value(jterm) = param
               return 
            endif 
          enddo
          iterm = iterm + Nterms
       enddo
       !
       ! spin-orbits
       do ifield = 1,Nspinorbits
          Nterms = spinorbit(ifield)%Nterms
          do jterm = 1,Nterms
            if (i==iterm+jterm) then 
               spinorbit(ifield)%value(jterm) = param
               return 
            endif 
          enddo
          iterm = iterm + Nterms
       enddo
       !
       ! L**2
       do ifield = 1,NL2
          Nterms = L2(ifield)%Nterms
          do jterm = 1,Nterms
            if (i==iterm+jterm) then 
               L2(ifield)%value(jterm) = param
               return 
            endif 
          enddo
       enddo
       !
       ! LxLy
       do ifield = 1,NLxLy
          Nterms = LxLy(ifield)%Nterms
          do jterm = 1,Nterms
            if (i==iterm+jterm) then 
               LxLy(ifield)%value(jterm) = param
               return 
            endif 
          enddo
       enddo
    !
  end subroutine one_parameter_copy_to_fields


  subroutine field_copy_to_parameters
    !
    integer(ik) :: ifield,Nterms,iterm
       !
       !potparam(:) = fitting%param(:)%value
       !ivar(:)     = fitting%param(:)%ifit
       !nampar(:)   = fitting%param(:)%name

       !
       iterm = 0
       !
       ! potential functions
       do ifield = 1,Nestates
          Nterms = poten(ifield)%Nterms
          potparam(iterm+1:iterm+Nterms) = poten(ifield)%value(1:Nterms)
          ivar(iterm+1:iterm+Nterms)     = nint(poten(ifield)%weight(1:Nterms))
          nampar(iterm+1:iterm+Nterms)   = poten(ifield)%name(1:Nterms)
          iterm = iterm + Nterms
       enddo
       !
       ! spin-orbits
       do ifield = 1,Nspinorbits
          Nterms = spinorbit(ifield)%Nterms
          potparam(iterm+1:iterm+Nterms) = spinorbit(ifield)%value(1:Nterms)
          ivar(iterm+1:iterm+Nterms)     = nint(spinorbit(ifield)%weight(1:Nterms))
          nampar(iterm+1:iterm+Nterms)   = spinorbit(ifield)%name(1:Nterms)
          iterm = iterm + Nterms
       enddo
       !
       ! L**2
       do ifield = 1,NL2
          Nterms = L2(ifield)%Nterms
          potparam(iterm+1:iterm+Nterms) = L2(ifield)%value(1:Nterms)
          ivar(iterm+1:iterm+Nterms)     = nint(L2(ifield)%weight(1:Nterms))
          nampar(iterm+1:iterm+Nterms)   = L2(ifield)%name(1:Nterms)
          iterm = iterm + Nterms
       enddo
       !
       ! LxLy
       do ifield = 1,NLxLy
          Nterms = LxLy(ifield)%Nterms
          potparam(iterm+1:iterm+Nterms) = LxLy(ifield)%value(1:Nterms)
          ivar(iterm+1:iterm+Nterms)     = nint(LxLy(ifield)%weight(1:Nterms))
          nampar(iterm+1:iterm+Nterms)   = LxLy(ifield)%name(1:Nterms)
          iterm = iterm + Nterms
       enddo
    !
  end subroutine field_copy_to_parameters

    !
 end subroutine sf_fitting

  subroutine MLlinur(dimen,npar,coeff,constant,solution,error)

  integer(ik),intent(in)  :: dimen,npar
  integer(ik),intent(out) :: error 
  real(rk),intent(in)  :: coeff(npar,npar),constant(npar)
  real(rk),intent(out) :: solution(npar)
  real(rk)          :: a0(npar,npar)
  real(rk)          :: c
  integer(ik)                   :: i1,i2,i,k8,k,l,k9

  !----- begin ----!
  
    do i1=1,dimen
    do i2=1,dimen 
       a0(i1,i2)=coeff(i1,i2)
    enddo
    enddo

    do i=1,dimen
      solution(i)=constant(i)
    enddo
    error=0
    do i=1,dimen
      c=0
      k8=i-1
      do k=1,k8
        c=c+a0(k,i)*a0(k,i)
      enddo

      if (c.ge.a0(i,i)) then
      !      write(6,*) '(',i,'-th adj. parameter is wrong )'
       error=i
       return
      endif

      a0(i,i)=sqrt(a0(i,i)-c)
      if (a0(i,i).eq.0) then
      !      write(6,*) '(',i,'-th adj. parameter is wrong )'
         error=i
         return
      endif
      k8=i+1
      do l=k8,dimen
         k9=i-1
         c=0
         do k=1,k9 
            c=c+a0(k,i)*a0(k,l)
         enddo
         a0(i,l)=(a0(l,i)-c)/a0(i,i)
      enddo
    enddo
    do i=1,dimen
      k8=i-1
      c=0
      do k=1,k8
         c=c+a0(k,i)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
    do i1=1,dimen
      i=1+dimen-i1
      k8=i+1
      c=0
      do k=k8,dimen
          c=c+a0(i,k)*solution(k)
      enddo
      solution(i)=(solution(i)-c)/a0(i,i)
    enddo
  return
  end subroutine MLlinur


  subroutine MLinvmat(al,ai,dimen,ierr)
  integer,intent(in)   :: dimen
  real(rk),intent(in)  :: al(dimen,dimen)
  real(rk),intent(out) :: ai(dimen,dimen)
  integer(ik),intent(out) :: ierr
  real(rk)             :: h(dimen),p,q
  integer(ik)          :: i1,k,i,j,k8,k9
      

    ierr = 0
    ai(1:dimen,1:dimen)=al(1:dimen,1:dimen)
 
    do i1=1,dimen
      k=dimen-i1+1
      p=ai(1,1)
      do i=2,dimen
        q=ai(i,1)
        !
        if (abs(p)<small_) then 
          !
          ierr = i
          !
          return
          !
        endif 
        !
        h(i)=q/p
        if(i.le.k) h(i)=-q/p
        do j=2,i
          k8=i-1
          k9=j-1
          ai(k8,k9)=ai(i,j)+q*h(j)
        enddo 
      enddo 
      ai(dimen,dimen)=1.0/p
      do i=2,dimen
        k8=i-1
        ai(dimen,k8)=h(i)
      enddo 
   end do 
   do i=1,dimen
     k8=i-1
     do j=1,k8
       ai(j,i)=ai(i,j)
     enddo 
   enddo 
   return
 end subroutine MLinvmat


    subroutine Robust_fit(a_wats,sigma,eps,wt)

      real(rk),intent(inout) :: a_wats
      real(rk),intent(in)    :: sigma(:),eps(:)
      real(rk),intent(inout) :: wt(:)
      !
      integer(ik)            :: npts,i,nrow,nused
      real(rk)               :: wtsum
      !
      if (verbose>=4) write(out,"(/'Robust fitting ...')")
      !
      npts = size(sigma)
      !
      nused = 0 
      do i=1,npts
        if (wt(i)>small_) nused=nused+1
      enddo
      !
      !Watson alpha-parameter
      ! 
      a_wats = 0.001_rk
      !
      if (verbose>=4) write(out,"('Watson parameter =',f18.8)") a_wats
      ! 
      !  adjusting the weights  
      ! 
      do nrow=1,npts
         if (wt(nrow)>small_) wt(nrow) = 1.0_rk/( sigma(nrow)**2 + a_wats*eps(nrow)**2 )
      enddo 
      ! 
      ! "re-normalizing" the weight factors
      !
      wtsum = sum(wt(1:npts))
      wt(1:npts) = wt(1:npts)/wtsum
      !
      if (verbose>=4) write(out,"('... done!')")
      !
    end subroutine Robust_fit
    !
    subroutine define_jlist
       !
       integer(ik)  :: nj,j,info,jind,iJmax,iJmin
       !
       !
       ! Count J-s
       !
       nj = 0 
       !
       do j = 1, size(fitting%J_list)
          !
          if (fitting%J_list(j)>-0.5_rk) nj = nj + 1
          !
       end do
       !
       iJmax = int(maxval(fitting%J_list(:)))
       iJmin = nint(minval(fitting%J_list(:),mask=fitting%J_list.gt.-1))
       !
       if (iJmin>0) nJ = nJ + 1
       !
       allocate(Jval(nJ),Jindex(0:ijmax),stat = info)
       if (info /= 0) stop 'refinement_by_fit allocation error: Jval - out of memory'
       !
       Jval = 0 
       jind = 0
       !
       do j = 1, size(fitting%J_list)
          !
          if (fitting%J_list(j)>=0) then
            !
            jind = jind + 1
            Jval(jind) = fitting%J_list(j)
            Jindex( int(Jval(jind)) ) = jind  ! round down half integer J, don't change integer J
            !
          endif 
          !
       end do
  end subroutine define_jlist

  subroutine lapack_svd_pseudo_inverse(tol,h,Nkeep)

    real(rk), intent(in)    :: tol 
    real(rk), intent(inout) :: h(:,:)  ! In:  matrix
    !                                          ! Out: pseudo-inverse matrix V Sigma- VT
    integer(ik), intent(out)  :: Nkeep
    character(len=1) :: jobu,jobvt,jobz
    !
    double precision,allocatable    :: work(:),u(:,:),vt(:,:),s(:),v(:,:),ut(:,:)
    integer           :: info
    integer           :: nh1, nh2,j,nu1,nu2,nvt1,nvt2,LDVT,LDU
    integer           :: lwork
    double precision  :: tol_
    double precision  :: alpha = 1.0d0,beta=0
    !
    jobu  = 'S'
    jobvt = 'A'
    jobz  = 'A'
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    lwork = 50*nh1
    !
    select case (jobu)
       case('A')
         LDU = nh1
         nu1 = nh1
         nu2 = nh1
       case('S')
         LDU = nh1
         nu1 = nh1
         nu2 = min(nh1,nh2)
    end select
    !
    select case (jobvt)
       case('A')
         LDVT = nh2
         nvt1 = nh2
         nvt2 = nh2
       case('S')
         LDVT = min(nh1,nh2)
         nvt1 = min(nh1,nh2)
         nvt2 = nh2
    end select
    !
    allocate(work(lwork),u(LDU,nu2),ut(nu2,LDU),vt(LDVT,nvt2),v(nvt2,LDVT),s(min(nh1,nh2)),stat=info)    !
    call ArrayStart('lapack_svd_pseudo_inverse',info,size(work),kind(work))
    call ArrayStart('lapack_svd_pseudo_inverse',info,size(u),kind(u))
    call ArrayStart('lapack_svd_pseudo_inverse',info,size(ut),kind(ut))
    call ArrayStart('lapack_svd_pseudo_inverse',info,size(vt),kind(vt))
    call ArrayStart('lapack_svd_pseudo_inverse',info,size(v),kind(v))
    call ArrayStart('lapack_svd_pseudo_inverse',info,size(s),kind(s))
    !
    call dgesvd(jobu,jobvt,nh1,nh2,h,nh1,s,u,LDU,vt,LDVT,work,-1,info)
    !call dgesdd(jobz,nh1,nh2,h,nh1,s,u,nh1,vt,nh2,work,-1,info)
    !
    if (int(work(1))>size(work)) then
      !
      lwork = int(work(1))
      !
      deallocate(work)
      !
      allocate(work(lwork))
      !
    endif
    !
    call dgesvd(jobu,jobvt,nh1,nh2,h,nh1,s,u,LDU,vt,LDVT,work,lwork,info)
    !call dgesdd(jobz,nh1,nh2,h,nh1,s,u,nh1,vt,nh2,work,lwork,info)
    !
    if (info/=0) then
      write (6,"(' dgesvd returned ',i8)") info
      stop 'lapack_dgesvd - dgesvd failed'
    end if
    !
    ! pseudo-inverse 
    !
    Nkeep = 0
    !
    tol_  = tol*maxval(s,dim=1)
    !
    Nkeep = minloc(s,dim=1,mask=s.ge.tol_)
    !
    Nkeep = min(Nkeep,nh1)
    !
    ut = 0 
    !
    do  j = 1,Nkeep
       ut(j,:) = U(:,j)/s(j)
    enddo
    !
    !v = transpose(vt)
    u = transpose(ut)
    !
    !ut(1:nh2,:) = matmul(v(1:nh2,1:nkeep),ut(1:Nkeep,:))
    !
    !h = 0
    !
    !h = transpose(ut)
    !
    !h(1:LDU,1:nvt2) = matmul(u(1:LDU,1:Nkeep),vt(1:nkeep,1:nvt2))
    !
    !call dgemm('T','N',nh2,nh1,nkeep,alpha,vt,nh2,ut,nh2,beta,ut,nh2)
    !h = transpose(ut)
    !
    call dgemm('T','N',LDU,nvt2,Nkeep,alpha,ut,nu2,vt,LDVT,beta,h,nh1)
    !
    deallocate(work,u,vt,ut,v,s)
    !
    call ArrayStop('lapack_svd_pseudo_inverse')
    !
  end subroutine lapack_svd_pseudo_inverse



  subroutine lapack_sdd_pseudo_inverse(tol,h,Nkeep)

    real(rk), intent(in)    :: tol 
    real(rk), intent(inout) :: h(:,:)  ! In:  matrix
    !                                          ! Out: pseudo-inverse matrix V Sigma- VT
    integer(ik), intent(out)  :: Nkeep
    character(len=1) :: jobu,jobvt,jobz
    !
    double precision,allocatable    :: work(:),u(:,:),vt(:,:),s(:),v(:,:),ut(:,:),A(:,:)
    integer,allocatable    :: iwork(:)
    integer           :: info
    integer           :: nh1, nh2,i,j,nu1,nu2,nvt1,nvt2,LDVT,LDU
    integer           :: lwork
    double precision  :: tol_
    !double precision  :: alpha = 1.0d0,beta=0
    !
    jobu  = 'S'
    jobvt = 'A'
    jobz  = 'A'
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    lwork = 50*nh1
    !
    select case (jobz)
       case('A')
         LDU = nh1
         nu1 = nh1
         nu2 = nh1
         LDVT = nh2
         nvt1 = nh2
         nvt2 = nh2
       case('S')
         LDU = nh1
         nu1 = nh1
         nu2 = min(nh1,nh2)
         LDVT = min(nh1,nh2)
         nvt1 = min(nh1,nh2)
         nvt2 = nh2
    end select
    !
    allocate(A(nh1,nh2),stat=info)
    call ArrayStart('lapack_sdd_pseudo_inverse',info,size(A),kind(A))
    !
    A = h
    !
    allocate(work(lwork),iwork(8*min(nh1,nh2)),u(LDU,nu2),ut(nu2,LDU),vt(LDVT,nvt2),v(nvt2,LDVT),s(min(nh1,nh2)),stat=info)
    call ArrayStart('lapack_sdd_pseudo_inverse',info,size(work),kind(work))
    call ArrayStart('lapack_sdd_pseudo_inverse',info,size(u),kind(u))
    call ArrayStart('lapack_sdd_pseudo_inverse',info,size(ut),kind(ut))
    call ArrayStart('lapack_sdd_pseudo_inverse',info,size(vt),kind(vt))
    call ArrayStart('lapack_sdd_pseudo_inverse',info,size(v),kind(v))
    call ArrayStart('lapack_sdd_pseudo_inverse',info,size(s),kind(s))
    !
    call dgesdd(jobz,nh1,nh2,A,nh1,s,u,LDU,vt,LDVT,work,-1,iwork,info)
    !
    if (int(work(1))>size(work)) then
      !
      lwork = int(work(1))
      !
      deallocate(work)
      !
      allocate(work(lwork))
      !
    endif
    !
    call dgesdd(jobz,nh1,nh2,A,nh1,s,u,LDU,vt,LDVT,work,lwork,iwork,info)
    !
    if (info/=0) then
      write (6,"(' dgesdd returned ',i8)") info
      stop 'lapack_dgesdd - dgesvd failed'
    end if
    !
    ! pseudo-inverse 
    !
    Nkeep = 0
    !
    tol_  = tol*maxval(s,dim=1)
    !
    Nkeep = minloc(s,dim=1,mask=s.ge.tol_)
    !
    Nkeep = min(Nkeep,nh1)
    !
    ut = 0 
    !
    do  i = 1,nu1
      do  j = 1,Nkeep
         ut(j,i) = U(i,j)/s(j)
      enddo
    enddo
    !
    !v = transpose(vt)
    u = transpose(ut)
    !
    !ut(1:nh2,:) = matmul(v(1:nh2,1:nkeep),ut(1:Nkeep,:))
    !
    !h = 0
    !
    !h = transpose(ut)
    !
    h(1:LDU,1:nvt2) = matmul(u(1:LDU,1:Nkeep),vt(1:nkeep,1:nvt2))
    !
    !call dgemm('T','N',nh2,nh1,nkeep,alpha,vt,nh2,ut,nh2,beta,ut,nh2)
    !h = transpose(ut)
    !
    !call dgemm('T','N',nh1,nkeep,nh2,alpha,ut,nh2,vt,nh2,beta,h,nh1)
    !
    deallocate(work,iwork,u,vt,ut,v,s,A)
    !
    call ArrayStop('lapack_sdd_pseudo_inverse')
    !
  end subroutine lapack_sdd_pseudo_inverse

  !
end module refinement
