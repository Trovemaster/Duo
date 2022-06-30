module me_numer
! 
! Numerical matrix elements based on the Numerov procedure
!
  use accuracy
  use lapack
  use timer
  implicit none
  private  

  public ik, rk, out
  public ME_numerov,simpsonintegral,simpsonintegral_rk,integral_rect_rk
  !
  !
  real(rk) :: enermax=  70000.0_rk   !  Upper limit for the energy search (1/cm)
  !
  real(rk), parameter :: enerdev=  0.1_rk        !  In the i-th solution has been already defined at 
                                                  !  previouse steps, we solve the eigen-problem anyway
                                                  !  but within +/- "enerdev" (cm-1) interval from 
                                                  !  the found solution
  real(rk), parameter :: enerstep = 100.0_rk    !  Energy interval needed for the step-wise search of the 
                                                  !  solution, when the iterative procedure does not help 
                                                  !
  real(rk), parameter :: enersmall= 0.000001_rk !  Small energy value to check whether the search interval not collapsed 
  !
  real(rk), parameter :: thrsh_int = 10.0_rk**(-rk) ! Threshold in Numerov-integration 
  real(rk), parameter :: thrsh_dif = 0.1e-2_rk       ! Threshold in numerical differentiation
  real(rk), parameter :: thrsh_upper  = 1.0e20_rk    ! Upper limit for fncts in out-numerov integr.
  !
  integer(ik), parameter :: itermax= 500000      ! Maximal possible number of iterations (tries) for Numerov 
  integer(ik), parameter :: maxslots = 10000       ! Maximal possible number of slots for found energies 
  integer(ik), parameter :: epoints = 40         ! N of points used for extrapolation in lsq_fit
  !
  real(rk) :: rhostep                           ! mesh size
  !
  real(rk),parameter :: wave_init=1.0_rk       ! initial arbitrary value at rho = rho_ref
  !
  !
  integer(ik), parameter :: iterlimit=5000       ! Iteration limit in numerov integration 
  !
  integer(ik) :: npoints                         ! number of grid points 
  real(rk)   :: rho_b(2)                        ! rhomin..rhomax
  integer(ik) :: iperiod                         ! the periodicity (can be negative for the reflecation)
  !
  integer(ik) :: verbose = 6                     ! Verbosity level
  !
  logical     :: periodic                        ! periodic boundary conditions: true or false
  !
  integer(ik) :: vmax,imin
  !
  integer(ik) :: io_slot                         ! unit numeber to store the numerov eigenvectors and their derivatives
  !
  integer(ik) :: imode                           ! the current mode
  !
  integer(ik) :: Nr = 4                          ! 2*Nr+1 is the number of interpolation points
  !
  integer(ik) :: iparity                         ! 
  !
  character(len=cl),parameter :: boundary_condition = 'UNBOUND'
  !
  character(len=cl),parameter :: deriv_method = 'ML_diffs' !  We use this method to estimate d pvi_v / d rho  
                                                         ! where phi_v is a numerical eigenfunction from numerov
                                                         ! deriv_method can be either 'd04aaf' of '5 points'
                                                         ! '5 points' 
  real(rk) ::  rho_switch  = .0174532925199432957692369_rk       ! the value of abcisse rho of the switch between regions (1 deg)

  integer(ik) :: iswitch                                 ! the grid point of switch


  !
  contains

  !
  ! Matrix elements with Numerov-eigenfunctions 
  !
  subroutine ME_numerov(vmax_,rho_b_,npoints_,numerpoints_,drho_,poten_,mu_rr_,icoord,iperiod_,&
                        enermax_,verbose_,energy,wavefunc)
   !
   integer(ik),intent(in)   :: vmax_,npoints_,numerpoints_,iperiod_
   real(rk),intent(out)    :: energy(0:vmax_),wavefunc(0:npoints_,0:npoints_)
   !
   real(rk),intent(in) :: rho_b_(2)
   real(rk),intent(in) :: poten_(0:npoints_),mu_rr_(0:npoints_),drho_(0:npoints_)
   integer(ik),intent(in) :: icoord ! coordinate number for which the numerov is employed
   real(rk),intent(in)    :: enermax_    !  Upper limit for the energy search (1/cm)
   integer(ik),intent(in) :: verbose_   ! Verbosity level
   !
   real(rk)            :: rho 
   real(rk)            :: h_t,sigma,sigma_t,rms,psipsi_t,characvalue,rhostep_,step_scale,fval,df_t
   !
   integer(ik) :: vl,vr,alloc,i,rec_len,i_,i1,i2
   !
   real(rk),allocatable :: phil(:),phir(:),dphil(:),dphir(:),phivphi(:),rho_kinet(:)
   real(rk),allocatable :: f(:),poten(:),mu_rr(:),d2fdr2(:),dfdr(:),rho_(:)
   character(len=cl)     :: unitfname 
   real(rk),allocatable :: enerslot(:),enerslot_(:)
   !
   integer,parameter :: NEND = 500,mxfs=5,mxprm=6,mxfsp=16001,mxisp=16001
   !
   real(rk) ::  BFCT,EFN,OVR,OVRCRT,RH,RMIN,VLIM
   integer(ik) :: ifs,IWR,JP,OTMF,TMFTYP(mxfs)
   real(rk) :: DER(0:mxprm-1),PSI(NEND),TMFPRM(0:mxprm-1,mxfs),VJ(mxfsp),z(mxisp)
   
   
    !
    if (verbose>=4) write (out,"(/'Numerov matrix elements calculations')")
     !
     ! global variables
     !
     vmax      = vmax_
     npoints   = numerpoints_
     rho_b = rho_b_
     iperiod = iperiod_
     verbose = verbose_
     imode=icoord
     enermax = enermax_
     !
     periodic = .false.
     if (iperiod>0) periodic = .true.
     !
     allocate(phil(0:npoints_),phir(0:npoints_),dphil(0:npoints_),dphir(0:npoints_), &
              phivphi(0:npoints_),rho_kinet(0:npoints_),enerslot(0:maxslots), &
              f(0:npoints),dfdr(0:npoints),d2fdr2(0:npoints),poten(0:npoints),mu_rr(0:npoints),stat=alloc)
     if (alloc/=0) then 
       write (out,"('phi - out of memory')")
       stop 'phi - out of memory'
     endif 
     !
     ! numerov step size 
     rhostep = (rho_b(2)-rho_b(1))/real(npoints,kind=rk)
     !
     ! real step size 
     rhostep_ = (rho_b(2)-rho_b(1))/real(npoints_,kind=rk)
     !
     ! define the rho-type coordinate 
     !
     rho_kinet(:) = drho_(:)
     !
     step_scale = rhostep/rhostep_
     !
     ! interpolation 
     !
     if (npoints==npoints_) then 
       !
       poten = poten_
       mu_rr = mu_rr_
       !
     else
       !
       allocate(rho_(0:npoints_),stat=alloc)
       if (alloc/=0) stop 'rho_ - out of memory'
       !
       forall(i_ = 0:npoints_) rho_(i_)  =  rho_b(1)+real(i_,kind=rk)*rhostep_
       !
       do i = 0,npoints 
          !
          rho =  rho_b(1)+real(i,kind=rk)*rhostep
          !
          i_ = int( real(i,rk)*step_scale )
          !
          i1 = max(0,i_-Nr) ; i2 = min(npoints_,i_+Nr)
          !
          call polintrk(rho_(i1:i2),poten_(i1:i2),rho,fval,df_t)
          poten(i) = fval
          !
          call polintrk(rho_(i1:i2),mu_rr_(i1:i2),rho,fval,df_t)
          mu_rr(i) = fval
          !
       enddo
       !
       deallocate(rho_)       !
     endif
     !
     ! Do some reporting
     !
     if (verbose>=3) then 
         write (out,"('vmax        = ',i8)") vmax
         write (out,"('icoord      = ',i4)") icoord
         write (out,"('rho_b (x)   = ',2f12.4)") rho_b(1:2) !*180.0_rk/pi
         write (out,"('rhostep (x) = ',2f12.4)") rhostep  !*180.0_rk/pi
     endif 
     !
     call diff_2d_4points_rk(npoints,rho_b,mu_rr,periodic,0_ik,dfdr,d2fdr2)
     !
     f(:) = (d2fdr2(:)/mu_rr(:)-dfdr(:)**2/mu_rr(:)**2*0.5_rk)*0.5_rk
     !
     ! Print out
     !
     if (verbose>=3) then 
        write(out,"('grid values (i,rho,rho_kinet,poten, mu_rr, f): ')") 
        do i_=0,npoints_,2
          i = int( real(i_,rk)/step_scale )
          rho = rho_b(1)+real(i_,kind=rk)*rhostep_
          write(out,"(i8,2f14.6,3g14.6)") i_,rho,rho_kinet(i_),poten(i),mu_rr(i),f(i)
        enddo 
     endif 
     !
     inquire(iolength=rec_len) phil(:),dphil(:)
     !
     write(unitfname,"('Numerov basis set # ',i6)") icoord
     call IOStart(trim(unitfname),io_slot)
     !
     open(unit=io_slot,status='scratch',access='direct',recl=rec_len)
     !
     ! solve 1D Schroedinger equation with Numerov algorithm
     !
     iparity = 0
     !
     iswitch = max(1,int( rho_switch/rhostep_))
     !
     if (iperiod/=0) vmax = vmax/2
     !
     call numerov(npoints,npoints_,step_scale,poten,mu_rr,f,enerslot)
     !
     psi = 0
     BFCT = 1.0
     der = 0 
     OTMF = -1
     !
     TMFTYP = 0
     !
     EFN = 48000.0
     BFCT = rhostep**2/mu_rr(1)
     VJ(1:npoints+1) = poten(0:npoints)*BFCT
     VLIM = poten(npoints)
     !MU(ISOT(iset))*RH*RH/16.85762908d0
     !
     !call OVRLAP(BFCT,DER,EFN,OVR,OVRCRT,PSI,RH,RMIN,TMFPRM,VJ,&
     !                   VLIM,z,ifs,IWR,JP,NEND,OTMF,TMFTYP)
     !
     if (iperiod/=0) then
       allocate(enerslot_(0:maxslots),stat=alloc)
       if (alloc/=0) then 
         write (out,"('phi - out of memory')")
         stop 'phi - out of memory'
       endif 
       ! 
       iparity = 1
       !
       call numerov(npoints,npoints_,step_scale,poten,mu_rr,f,enerslot_)
       !
       do vl = vmax,0,-1
         !
         enerslot(2*vl  ) = enerslot(vl)
         if (2*vl+1<=vmax_) enerslot(2*vl+1) = enerslot_(vl)
         !
       enddo
       !
       vmax = vmax_
       !
       deallocate(enerslot_)
       !
     endif
     !
     ! Matrix elements 
     !
     sigma = 0.0_rk 
     rms   = 0.0_rk 
     characvalue = maxval(enerslot(0:vmax))
     energy(0:vmax) = enerslot(0:vmax)
     if (trim(boundary_condition)/='UNBOUND') then
       energy(0:vmax) = energy(0:vmax)-energy(0)
     endif
     !
     do vl = 0,vmax
        !
        read (io_slot,rec=vl+1) (phil(i),i=0,npoints_),(dphil(i),i=0,npoints_)
        !
        wavefunc(:,vl) = phil(:)
        !
        do vr = vl,vmax
            !
            if (iperiod/=0.and.mod(abs(vl-vr),2)==1) cycle
            !
            if (vl==vr) then
                phir =  phil
               dphir = dphil
            else
               read (io_slot,rec=vr+1) (phir(i),i=0,npoints_),(dphir(i),i=0,npoints_)
            endif
            !
            ! Here we prepare integrals of the potential 
            ! <vl|poten|vr> and use to check the solution of the Schroedinger eq-n 
            ! obtained above by the Numerov
            !
            phivphi(:) = phil(:)*poten_(:)*phir(:)
            !
            h_t = simpsonintegral_rk(npoints_,rho_b(2)-rho_b(1),phivphi)
            !
            ! momenta-quadratic part 
            !
            phivphi(:) =-dphil(:)*mu_rr_(:)*dphir(:)
            !
            psipsi_t = simpsonintegral_rk(npoints_,rho_b(2)-rho_b(1),phivphi)
            !
            ! Add the diagonal kinetic part to the tested mat. elem-s
            !
            h_t = h_t - 0.5_rk*psipsi_t
            !
            !
            phivphi(:) = phil(:)*phir(:)
            !
            psipsi_t = simpsonintegral_rk(npoints_,rho_b(2)-rho_b(1),phivphi)
            !
            ! Count the error, as a maximal deviation sigma =  | <i|H|j>-E delta_ij |
            !
            sigma_t =  abs(h_t)
            if (vl==vr) sigma_t =  abs(h_t-enerslot(vl))
            !
            sigma = max(sigma,sigma_t)
            rms = rms + sigma_t**2
            !
            ! Now we test the h_t = <vl|h|vr> matrix elements and check if Numerov cracked
            ! the Schroedinger all right
            if (vl/=vr.and.abs(h_t)>sqrt(small_)*abs(characvalue)*1e4.and.trim(boundary_condition)/='UNBOUND') then 
               write(out,"('ME_numerov: wrong Numerovs solution for <',i4,'|H|',i4,'> = ',f20.10)") vl,vr,h_t
               stop 'ME_numerov: bad Numerov solution'
            endif 
            !
            if (vl==vr.and.abs(h_t-enerslot(vl))>sqrt(small_)*abs(characvalue)*1e4.and.trim(boundary_condition)/='UNBOUND') then 
               write(out,"('ME_numerov: wrong <',i4,'|H|',i4,'> (',f16.6,') =/= energy (',f16.6,')')") vl,vr,h_t,enerslot(vl)
               stop 'ME_numerov: bad Numerov solution'
            endif 
            !
            ! Reporting the quality of the matrix elemenst 
            !
            if (verbose>=3.and.trim(boundary_condition)/='UNBOUND') then 
              if (vl/=vr) then 
               write(out,"('<',i4,'|H|',i4,'> = ',e16.2,'<-',8x,'0.0',5x,'; <',i4,'|',i4,'> = ',e16.2,'<-',8x,'0.0')") & 
                                vl,vr,h_t,vl,vr,psipsi_t
              else
               write(out,"('<',i4,'|H|',i4,'> = ',f16.6,'<-',f16.6,'; <',i4,'|',i4,'> = ',f16.6)")& 
                              vl,vr,h_t,enerslot(vl),vl,vr,psipsi_t
              endif 
            endif 
            !
        enddo
     enddo
     !
     rms = sqrt(rms/real((vmax+1)*(vmax+2)/2,kind=rk))
     !
     if (verbose>=1) then 
        write(out,"('Maximal deviation sigma =  | <i|H|j>-E delta_ij | is ',f18.8)") sigma
        write(out,"('rms deviation is ',f18.8)") sqrt(rms)
     endif 
     !
     deallocate(phil,phir,dphil,dphir,phivphi,rho_kinet,enerslot,f,poten,mu_rr,d2fdr2,dfdr)
     !
     !
  end subroutine ME_numerov

!
!********************************************************************!
!                                                                    !
! numerov calculates bending functions and stores them               !
!                                                                    !
!                                                                    !
!********************************************************************!
  subroutine numerov(npoints,npoints_,step_scale,poten,mu_rr,f,enerslot) 
   !
   integer(ik),intent(in) :: npoints,npoints_
   real(rk),intent(in) :: step_scale
   real(rk),intent(in) :: poten(0:npoints),mu_rr(0:npoints),f(0:npoints)
   !
   real(rk),intent(out) :: enerslot(0:maxslots)

   integer(ik) :: v,i_
   integer(ik) :: ierr,numnod,ntries,vrec,k0,ipoint,mtries
   !
   real(rk) :: eold,tsum
   !
   real(rk) :: pcout
   !
   real(rk),allocatable :: pot_eff(:),i0(:),phi_f(:),phi_t(:),psi_t(:),phi_f_(:),phi_d_(:)
   !
   real(rk) :: eguess,enerlow,enerupp,enerdelta,enershift,enermid
   !
   !real(rk) :: simpsonintegral
   !
   logical :: notfound
   !
   real(rk) :: potmin,oldphi,newphi,v_t(-2:2),amplit,amplit0,diff
   !
   real(rk) :: df_t1,df_t2,deltaE = 1.0_rk
   !
   integer(ik) :: alloc,i,iter,nm2,ic,imin_l,imin_r
   integer(ik) :: icslots(0:maxslots),iright,ileft,ireflect
   !
   !
   allocate(pot_eff(0:npoints),i0(0:npoints),phi_f(0:npoints),phi_t(0:npoints),psi_t(-30:npoints+30),phi_f_(0:npoints),&
            phi_d_(0:npoints),stat=alloc)
   if (alloc/=0) then 
      write (6,"('numerov: allocation is faild - out of memory')")
      stop 'numerov: allocation is faild - out of memory'
   endif 
   !
   ! fcoef = planck*avogno*1.0d+16/( 4.0d+00*pi*pi*vellgt )
   !
   ! Check for the potential "1/0" probelm in mu_rr
   !
   do i=0,npoints
     if ( abs(mu_rr(i))<small_*100.0_rk ) then 
        write(out,"('numerov: mu_rr is zero at i = ',i8)") i
        stop 'numerov: mu_rr has a zero element'
     endif
   enddo 
   !
   ! before commencing the numerov-cooley numerical integration we
   ! store the following function in the f1 array:
   !
   !                       0                 
   !     w  = f  (p ) + 2 i  v  (p )   
   !      i    1   i          0   i
   !
   !      in the renner-teller version of morbid, this function is
   !      extended with a term originating in the non-zero
   !      electronic angular momentum.
   !
   !                  -1       0                 -1 
   ! v  is given in cm  , and i   is given in  cm
   !  0                         
   !
   ! calculate effective potential w_i
   !
   pot_eff(:) = f(:) + 2.0_rk*( poten(:) )/mu_rr(:)
   !
   ! Effective inertia moment
   !
   i0(:) = 2.0_rk/mu_rr(:)
   !
   ! find minimum : position 
   !
   potmin = safe_max
   !
   do i=0,npoints
      !
      if (pot_eff(i)<potmin) then 
         imin = i
         potmin = pot_eff(i)
      endif
      !
   enddo
   ! and value 
   !
   potmin = pot_eff(imin)
   !
   if(trim(boundary_condition)=='UNBOUND') then 
     potmin = pot_eff(npoints)
     imin = npoints
   endif
   !
   if (imin<0.or.imin>npoints) then 
       write(out,"('numerov: pot_eff has no minimum',i8)") 
       stop 'numerov: pot_eff has no minimum'
   endif 
   !
   imin_l = minloc(poten(0:imin-1),dim=1)
   imin_r = minloc(poten(imin+1:npoints),dim=1)+imin
   !
   ileft  = npoints
   iright = 0
   !
   ! Print out
   !
   if (verbose>=5) then 
      write(out,"('grid values (poten, mu_rr, poten_eff): ')") 
      do i=0,npoints,2
        write(out,"(i8,3f18.8)") i,poten(i),mu_rr(i),pot_eff(i)
      enddo 
      write(out,"('potmin(eff) =  ',f18.8,' at i = ',i8)") potmin,imin
   endif 
   !
   !
   ! Parameters for the Simson rule integration 
   !
   ! facodd=2.0_rk*rhostep/3.0_rk
   ! faceve=2.0_rk*facodd
   !
   ! Here we start solving Schroedinger equation for the v-th eigenvalue
   !
   select case (trim(boundary_condition))
     !
   case default
     !
     ! assign initial values 
     !
     enerslot(0:maxslots)=poten(imin)-sqrt(safe_max)
     icslots(0:maxslots) = npoints
     !
     enerlow=potmin/2.0_rk*mu_rr(imin)
     enerupp=enermax
     !
   case ('UNBOUND')
     !
     ! assign initial values 
     !
     !enermax = maxval(poten(:))
     !
     enerlow=potmin/2.0_rk*mu_rr(imin)
     !
     deltaE = (enermax-enerlow)/(npoints+1)
     !
     do i=0,min(npoints,maxslots)
        !
        enerslot(i) = enerlow+real(i,rk)*deltaE
        !
     enddo
     !
     icslots(0:maxslots) = npoints-2 !/2
     enerupp=enermax
     !
     !if (imin==npoints) imin = npoints/2
     !
   end select
   !
   do v=0,vmax+1
     !
     if (v/=0) enerlow = enerslot(v-1)
     !
     ! Determine the highest undefined energy slot 
     i = v
     do while (enerslot(i)<poten(imin).and.i<maxslots)
       i = i + 1
     enddo 
     !
     if (enerslot(i)>=poten(imin)) then
        enerupp=enerslot(i)
     else 
        enerupp=enermax
     endif 
     !
     if (enerslot(v)>poten(imin)) then
        !
        enerlow=max(enerslot(v)-enerdev,poten(imin))
        enerupp=min(enerslot(v)+enerdev,enermax)
        ic  = icslots(v)
        !
     endif
     !
     select case (trim(boundary_condition))
       !
     case ('UNBOUND')
       !
       ic  = icslots(v)
       enerupp=enerslot(v+1)
       enerlow=max(enerslot(v)-deltaE,poten(imin))
       enerlow = enerslot(v)
       !  
     end select 
     !
     enerdelta=  enerstep
     enershift = enerstep
     ntries=0
     mtries = 0
     !
     ! start itererative procedure for searching the current solution 
     !
     iter = 0
     !
     notfound = .true.
     !
     search_loop : do while(notfound.and.iter<itermax) 
       !
       iter = iter + 1
       !
       enermid=(enerlow+enerupp)/2.0
       !
       !enerdelta = enerdelta*0.5_rk
       !
       eguess=enermid
       !
       eold=eguess
       !
       ! Numerov procedure 
       !
       ic = icslots(v)
       !
       if (ic/=npoints.and.enerslot(v)>=poten(imin).and.iter==1) eguess = enerslot(v)
       !
       if (verbose>=5) write(out,"('eguess = ',e14.7)") eguess
       !
       call numcoo ( v, pot_eff, i0, eguess, enerlow, ic, phi_f, pcout, ierr)
       ! 
       if (ierr>1) then 
           !
           write(out,"('numerov: no solution found in numcoo, ierr = ',i8)") ierr
           !
           do i = 0,maxslots
             if (enerslot(i)>poten(imin)) write (out,"('        v,ener = ',i8,f20.10)") i,enerslot(i)
           enddo 
           !
           stop 'numerov: no solution found in numcoo'
           !
       endif 
       ! 
       !
       if (ierr == 0) then  
          !
          numnod=0
          nm2=npoints-2
          oldphi=phi_f(0)
          newphi=phi_f(2)
          !
          do i=2,nm2
             if ( oldphi*phi_f(i)<0.0 .or.  & 
                  ( phi_f(i)==0.0_rk  .and. newphi*oldphi<0.0_rk) ) numnod=numnod+1
             newphi=phi_f(i+2)
             oldphi=phi_f(i)
          enddo
          !
       endif 
       !
       if (ierr==0.and.iperiod<0) then 
           !
           phi_t(:) = phi_f(:)/sqrt(mu_rr(:))
           !
           tsum = simpsonintegral_rk(npoints,rho_b(2)-rho_b(1),phi_t(:)**2)
           !
           phi_t(:)=phi_t(:)/sqrt(tsum)
           !
           ireflect = 1
           if (iparity/=0) ireflect = -1
           !
           call diff_2d_4points_rk(Npoints,rho_b,phi_t,periodic,ireflect,psi_t(0:Npoints))
           !
           if ( (iparity==0.and.abs(psi_t(npoints  ))>sqrt(thrsh_int)).or.&
                (iparity==1.and.abs(phi_t(npoints  ))>sqrt(thrsh_int))) then 
              !
              if (verbose>=5) write(out,"(/'phi(N) and psi(N) = ',2g18.8,', energy = ',f12.4,', n = ',9i7)") &
                              phi_t(npoints),psi_t(npoints ),eguess,numnod
              !
              if ( numnod+1<maxslots.and.numnod>v) then 
                  enerslot(numnod)=eguess
                  icslots(numnod) = ic 
              endif 
              !
              numnod = maxslots+1
              !
           endif 
           !
       endif
       !
       if (verbose>=5) write(out,"('v = ',i5,'; numnod = ',i8,', efound = ',f16.7,', ierr = ',i9)") v,numnod,eguess,ierr
       !
       select case (trim(boundary_condition))
         !
       case default
         !
         if ( ierr==0.and.numnod+1<maxslots.and.numnod>=v ) then
            !
            enerslot(numnod)=eguess
            icslots(numnod) = ic 
            !
            if (verbose>=6) then 
               !
               write (out,"('v,numnod,ener = ',2i8,f20.10)") v,numnod,enerslot(numnod)
               !
               do i=0,npoints 
                  !
                  write(out,"(i8,2f18.8)") i,phi_f(i),exp((rho_b(1)+rhostep*real(i,rk)))
                  !
               enddo
               !
            endif 
         endif 
         !
         if ( ierr == 0.and.v==numnod ) notfound = .false.
         ! if it cannot find a solution by Numerov procedure and ierr = 1, we can still try to
         ! devide the searching interval and repeate
         ! the Numerov procedure with a closer initial guess
         ! We do the same thing if the found solution is not what we need (/=v)
         !
         if (abs(enerupp-enerlow) < enersmall.and.(ierr == 1 .or. numnod  /= v)) then 
            !
            ntries=0
            !
            !enerdelta = enerdelta*0.5_rk
            !
            ! Are we still at the begining?
            !
            if (v.eq.0) then
               !
               enerlow=poten(imin)
               !
            else
              !
              enerlow=enerslot(v-1)-0.001
              !
            endif
            !
            ! Determine the highest undefined energy slot 
            i = v+1
            do while (enerslot(i)<poten(imin).and.i<maxslots)
              i = i + 1
            enddo 
            !
            if (enerslot(i)>=poten(imin)) then
               enerupp=enerslot(i)
            else 
               enerupp=enermax
            endif 
            !
            if (enermid>=enerupp) then 
               enerdelta = enerdelta*0.5_rk
               enershift = 0.0_rk
            endif 
            !
            enerlow=enerlow+enershift
            enershift=enershift+enerdelta
            enerupp=enerlow
            !
            cycle search_loop
            !
         endif
         !
         if (ierr <=1 .and. numnod  /= v) then 
            !
            ! We need to do something if the seacrh interval is collapsed
            ! this must have happend at ierr=1, and no solution found 
            !
            if (mtries==0) then
               !
               mtries = 1
               ! 
               ! Are we still at the begining?
               !
               if (v.eq.0) then
                  enerlow=poten(imin)
               else
                  enerlow=enerslot(v-1) 
               endif
               !
               ! Determine the highest undefined energy slot 
               i = v+1
               do while (enerslot(i)<poten(imin).and.i<maxslots)
                 i = i + 1
               enddo 
               !
               if (enerslot(i)>=poten(imin)) then
                  enerupp=enerslot(i)
               else 
                  enerupp=enermax
               endif 
               !
            endif 
            !
            if (numnod<v.or.ierr==1) then
               enerlow=enermid
            else
               enerupp=enermid
            endif
            !
         endif
         !
         ! multiply the wavefunction with sqrt(irr) (eq. (6.4) of jensen)
         ! 
         phi_f(:) = phi_f(:)/sqrt(mu_rr(:))
         !
         phi_t(:) = phi_f(:)*phi_f(:)
         !
         !   numerical intagration with simpson's rule #2
         !
         tsum = simpsonintegral_rk(npoints,rho_b(2)-rho_b(1),phi_t)
         !
         phi_f(:)=phi_f(:)/sqrt(tsum)
         !
       case ('UNBOUND')
         !
         if (verbose>=6) then 
            !
            write (out,"('v,numnod,ener = ',2i8,f20.10)") v,numnod,enerslot(numnod)
            !
            do i=0,npoints 
               !
               write(out,"(i8,2f18.8)") i,phi_f(i),exp((rho_b(1)+rhostep*real(i,rk)))
               !
            enddo
            !
         endif 
         !
         if ( ierr == 0 ) notfound = .false.
         !
         ! multiply the wavefunction with sqrt(irr) (eq. (6.4) of jensen)
         ! 
         phi_f(:) = phi_f(:)/sqrt(mu_rr(:))
         !
         oldphi = phi_f(0)
         newphi = phi_f(1)
         amplit = 0 
         amplit0 = 0
         do i=2,npoints 
            oldphi=newphi
            newphi=phi_f(i-1)
            if ( oldphi<newphi.and.newphi<=phi_f(i) ) then
               diff = abs(amplit-amplit0)
               !
               amplit0 = amplit
               amplit = newphi
               !
            endif
         enddo

         !
         phi_t(:) = phi_f(:)*phi_f(:)
         !
         !   numerical intagration with simpson's rule #2
         !
         tsum = simpsonintegral_rk(npoints,rho_b(2)-rho_b(1),phi_t)
         !
         if (tsum>small_*100) then
           phi_f(:)=phi_f(:)/sqrt(tsum)
         endif
         !
       end select
       !
     enddo search_loop
     !
     if( notfound ) then
       !
        write(out,"('numerov: no solution found after ',i8,' iterations')") itermax
        write (out,"('        v = ',i8,', eguess= ',e20.10)") v,eguess
        !
        do i = 0,maxslots
          if (enerslot(i)>poten(imin)) write (out,"('        v,ener = ',i8,f20.10)") i,enerslot(i)
        enddo 
        !
        stop 'numerov: no solution found'
        !
     endif 
     !
     if (ierr==0.and.periodic) then 
        if (verbose>=5) write(out,"(/'phi_t(0)-phi_t(npoints  )    : ',3g18.8)") phi_t(0),phi_t(npoints  ),phi_t(0)-phi_t(npoints  )
        if (verbose>=5) write(out,"(/'dphi_t(0)-dphi_t(npoints)    : ',5g18.8)") phi_t(1),phi_t(npoints-1),df_t1,df_t2,df_t1-df_t2
     endif 
     !
     !if (iperiod<0) phi_f(:)=phi_f(:)/sqrt(2.0_rk)
     !
     !
     if (verbose>=5) then 
        write (out,"('v,ener = ',i8,f20.10)") v,enerslot(numnod)
     endif 
     !
     if (verbose>=6) then 
        write(out,"(f18.8)") phi_f
     endif 
     !
     ! direvatives of the wavefunctions method I
     !
     if(deriv_method == '5 points') then 
        !
        do i=0,npoints      ! --- np
           !
           do k0 = -2, 2
              !
              if (periodic) then 
                 ipoint = min(max(i+k0,0),npoints)
                 v_t(k0) = phi_f(ipoint)
              else
                 if (i+k0<0) then 
                    v_t(k0) = phi_f(-(i+k0))
                 elseif(i+k0>npoints) then 
                    v_t(k0) = phi_f(npoints)
                 else
                    v_t(k0) = phi_f(i+k0)
                 endif
              endif
              !
           enddo
           !
           ! 5-points expression
           !
           phi_t(i) = (-v_t( 2)/12.0_rk+2.0_rk/3.0_rk*v_t(1) & 
                       +v_t(-2)/12.0_rk-2.0_rk/3.0_rk*v_t(-1) )/rhostep
           !
           !phi_t(i) = (v_t(1)-v_t(-1) )/(rhostep*2.0_rk)
           !
        enddo 
        !
        ! direvatives of the wavefunctions method II
        !
     elseif(deriv_method == 'ML_diffs') then 
        !
        if (iperiod<0) then 
          !
          ireflect = 1
          if (iparity/=0) ireflect = -1
          !
          call diff_2d_4points_rk(Npoints,rho_b,phi_f,periodic,ireflect,phi_t)
          !
        else 
          !
          ireflect = 0
          if (iparity/=0) ireflect = -1
          !
          call diff_2d_4points_rk(Npoints,rho_b,phi_f,periodic,ireflect,phi_t)
          !
        endif
        !
        ! direvatives of the wavefunctions method III
        !
     else
         write(out,"('numerov: bad deriv_method ',a)") trim(deriv_method)
               stop 'numerov: bad deriv_method'

     endif
     !
     !
     if (verbose>=5) then 
        !
        do i=0,npoints 
           !
           write(out,"(i8,2f18.8)") i,phi_f(i),phi_t(i)
           !
        enddo
        !
     endif 
     !
     ! Find the outmost points where the function is not zero
     !
     i = 0 
     !
     do while (abs(phi_f(i))<100.0*sqrt(small_).and.i<npoints)
      i = i + 1
     enddo
     !
     ileft = min(ileft,i)
     !
     i = npoints 
     !
     do while (abs(phi_f(i))<100.0*sqrt(small_).and.i>0)
      i = i - 1
     enddo
     !
     iright = max(iright,i)
     !
     ! dump the eigenfunction
     !
     if (npoints_==npoints) then 
          phi_f_ = phi_f
          phi_d_ = phi_t
     else
       do i_ = 0,npoints_
          !
          i = nint( real(i_,rk)/step_scale )
          phi_f_(i_) = phi_f(i)
          phi_d_(i_) = phi_t(i)
          !
       enddo
     endif
     !
     vrec=v+1
     !
     if (iperiod/=0) then 
       !
       vrec = 2*v+1 
       !
       if (iparity/=0) vrec = 2*v+2
       !
     endif
     !
     write (io_slot,rec=vrec) (phi_f_(i_),i_=0,npoints_),(phi_d_(i_),i_=0,npoints_)
     !
   enddo
   !
   write(out,"(/' Outmost points of the nonvanishing wave-functions are: [',i6,'...',i6,'], i.e. [',f12.6,'...',f12.6,'] ')") &
                 ileft,iright,rho_b(1)+rhostep*real(ileft,rk),rho_b(1)+rhostep*real(iright,rk)
   !
   ! report the found Numerov energies 
   !
   write (out,"(/' Numerov-energies are:')") 
   !
   enerlow = minval(enerslot(0:vmax),dim=1) !mask=(enerslot>=poten(imin)
   !
   do v=0,vmax   
     write (out,"(i8,f18.8)") v,enerslot(v)-enerlow
   enddo
   !
   write (out,"('Zero-point-energy is :',f18.6)") enerlow
   !
   deallocate(pot_eff,i0,phi_f,phi_t,phi_f_,phi_d_)
   !
  end subroutine numerov



!     subroutine solves eqn 6.5 , for the given value of kqua ,
!     by the numerov - cooley technique.
!
!                  (                       0                  )
!      2           (                    2 i                   )
!     d            (                2      pp                 )
!     ---phi (p) = ( f (p) + f (p) k  + ----- ( v   (p) - e ) ) phi (p)
!       2   b      (  1       2         _2       eff          )    b
!     dp           (                    h                     )
!
!     for this the following defintions (equations 6.14 to 6.17) are
!     needed :-
!
!     ph  = phi  (p )
!       i      b   i
!
!             0         _2
!     i  = 2 i   (p ) / h
!      i      pp   i
!
!                             2      0                   _2
!     u  = f  (p ) + f  (p ) k  + 2 i   (p ) v    (p ) / h
!      i    1   i     2   i          pp   i   eff   i
!
!                 2
!     y  = ( 1 - h  ( u  - i  e ) ) ph
!      i         --    i    i         i
!                12
!
!     starting from an initial guess for the energy , given in wavfun
!     an iterative scheme is followed. this is performed by using the
!     recursion relation (eqn 6.18) :-
!
!                           2
!     y    + y    - 2 y  = h  ( u  - i  e ) ph
!      i+1    i-1      i         i    i       i
!
!     h is here the rhostep length used to generate the f1 , f2 , i0's
!
!     to calculate values of y(i) and phi_f(i) both from rho = 0
!     and from rho = rhomax. these two parts are called the outward
!     and inward integrations respectively.
!     the outward integration is carried out until a maximum is reached
!     in the wavefunction. this crossing point (ic , yc , pc) is then
!     the stopping point used for the inward integration.
!     both sets of wavefunctions are scaled so that pc(out) = pc(in) = 1
!     an error for the energy can be calculated from eqn 6.24 :-
!
!                                                         npoints
!               out           in       2                    \--      2
!     d(e) = (-y    + 2 y  - y    ) / h + u  - i  e ) y  /   \  i  ph
!               c-1      c    c+1          c    c      c     /   i   i
!                                                           /--
!                                                           i=1
!
!     and this is added to the eguess at each iteration , when d(e)
!     is less than thrsh3 the iteration has converged.
!

  subroutine numcoo( v, pot_eff, i0, eguess, enerlow, iref, phi_f, pcout, ierr )

     integer(ik),intent(in) :: v
     integer(ik),intent(inout) :: iref
     real(rk),intent(in)  ::  i0(0:npoints),pot_eff(0:npoints)
     real(rk),intent(in)  ::  enerlow
     real(rk),intent(inout)  :: eguess
     real(rk),intent(out) :: phi_f(0:npoints)
     real(rk),intent(out) :: pcout
     integer(ik),intent(out) :: ierr

     real(rk) :: hh,dx,sumout,sumin,tsum,ycm1,pcin,ycp1,yc,phi_t,k_coeff
     integer(ik) :: niter,ic,istart,i,id,i0_,jt,i1,i2
     real(rk) :: G1,G0,DI,SG0,SG1,SI,VV,Y1,Y2,S0,GI,SGI
     !
     !y(i)=(1.0_rk-hh*(pot_eff(i)-i0(i)*eguess)/12.0_rk)*phi_f(i)
     !
     hh=rhostep*rhostep
     !
     !
     niter=0
     dx = safe_max
     !
     pcout = 1.0_rk
     pcin  = 1.0_rk
     !
     if (mod(v,2)/=0)  pcout = -pcout
     !
     if ( .not.periodic.and.iperiod<0.and.iparity==0) pcout = 1.0_rk
     !
     if (periodic.and.iperiod>0) pcout = pcout*(-1.0_rk)**iparity
     !
     ierr = 0
     istart = 1
     !
     !
     do while ( ierr==0 .and. abs(dx)>thrsh_int  )   ! --- niter 
        !
        select case (trim(boundary_condition))
          !
        case('BOUND-0')
          !
          if (abs(rho_b(1))>0.01_rk) then 
            stop 'illegal combination of conditions for BOUND-0'
          endif
          !
          phi_f(1)  = sqrt(small_) ! small_
          phi_f(0)  = small_
          !
          phi_f(npoints-1)  = sqrt(small_) ! small_
          phi_f(npoints  )  = small_
          !
        case('PERIODIC')
          !
          if (pcin*pcout>small_.or.(mod(v,2)/=0.and.iperiod>0.and.iparity==1)) then 
            !
            phi_t = wave_init
            !
            !
            phi_f(npoints  )  = phi_t
            phi_f(npoints-1)  = phi_f(npoints)*(12.0_rk+5.0_rk*hh*(pot_eff(npoints)-i0(npoints)*eguess))/ &
                                               (12.0_rk-        hh*(pot_eff(npoints-1)-i0(npoints-1)*eguess))
            !
            phi_f(0)  = phi_t

            if ((mod(v,2)/=0.and.iperiod>0.and.iparity==1)) phi_f(0) = -phi_f(0)

            phi_f(1)  = phi_f(0)*(12.0_rk+5.0_rk*hh*(pot_eff(0)-i0(0)*eguess))/ &
                                 (12.0_rk-        hh*(pot_eff(1)-i0(1)*eguess))

           !
          else
            !
            phi_f(npoints-1)  = small_ !sqrt(small_) ! small_
            phi_f(npoints  )  = 0
            !
            phi_f(1)  = -small_
            phi_f(0)  = 0
            !
          endif 
          !
        case ('UNBOUND')
          !
          ic = npoints
          !
          ID= nint(0.2_rk/rhostep)
          !
          do i = id+1,npoints,1 !id   
            jt = i
            G1= pot_eff(i)-i0(i)*eguess
            if (G1 < small_) exit
            i0_ = i
            G0= G1
          enddo
          !
          !i1 = jt
          !jt = I0+(I1-I0)*G0/(G0-G1)
          !
          DI= 10.0_rk/(hh*(pot_eff(jt-1)-pot_eff(jt)))**(1.0_rk/3.0_rk)
          ID= nint(DI)
          I0_= max(0,jt-id)
          IF(I0_>npoints) stop 'Airy I0_>npoints'
          G0=  (pot_eff(I0_)-i0(I0_)*eguess)*hh
          !
          ! Adjust starting point outward to ensure integration scheme stability
          do while(G0>10.0_rk.and.i0_<npoints)
             i0_ = i0_+1
             G0= (pot_eff(i0_)-i0(i0_)*eguess)*hh
          enddo
          !
          IF(I0_>=npoints) stop 'Airy I0_>npoints'
          !
          i1= i0_+1
          i2= i1+1
          
          !** WKB starting condition for wave function    
          S0= 1.0_rk
          G0= (pot_eff(i0_)-i0(i0_)*eguess)*hh
          GI= (pot_eff(i1)-i0(i1)*eguess)*hh
          !
          if ((G0<=small_).or.(GI<=small_)) then
             VV= pot_eff(i1)*i0(i1)
             S0= 0
             SI= 1.0_rk
          else
             !
             SG0= sqrt(G0)
             SGI= sqrt(GI)
             SI= S0*sqrt(SG0/SGI)*exp((SG0+SGI)*0.5_ark)  
             if (SI<small_) S0=0
             Y1= S0*(1.0_rk-1.0_ark/12.0_ark*G0)
             Y2= SI*(1.0_rk-1.0_ark/12.0_ark*GI)
             !
          endif
          !
          istart = i0_+1
          !
          phi_f = 0
          !
          phi_f(istart-1)  = S0
          phi_f(istart)  = SI
          !
          if (i0(npoints)*eguess-pot_eff(npoints)<-small_) then 
            write(out,"('Error-Numerov: wrong usage of unbound integration for bound state ')")
            stop 'Error-Numerov: wrong usage of unbound integration for bound state '
          endif
          !
          k_coeff  = sqrt(i0(npoints)*eguess-pot_eff(npoints))
          !k_coeff  = sqrt(i0(npoints))*sqrt(eguess-pot_eff(imin)/i0(imin))
          !
          phi_f(npoints  )  = sin(k_coeff*rho_b(2))
          phi_f(npoints-1)  = sin(k_coeff*(rho_b(2)-rhostep))
          !
        case ('QUASI-BOUND')
          !
          stop 'QUASI-BOUND has not been implemented yet'
          !
        case default 
          !
          phi_f(0)  = 0.0_rk
          phi_f(1)  = small_
          phi_f(npoints-1)  = safe_min ! small_
          phi_f(npoints  )  = 0.0_rk
          !
        end select
        !
        ! Outer integration 
        !
        !
        call intout ( v, pot_eff, i0, eguess, phi_f, sumout, istart, iref,  ic, pcout,  ierr)
        !
        if (ierr/=0) return
        !
        ! Inner integration 
        !
        call intin ( pot_eff, i0, eguess, phi_f, sumin, ic, pcin)
        !
        tsum=0
        !
        if (sumin>small_) then 
           tsum=sumin/(pcin*pcin)
        endif
        !
        if (sumout>small_) then 
           tsum=tsum+sumout/(pcout*pcout)
        endif
        !
        !tsum=sumin/(pcin*pcin)+sumout/(pcout*pcout)
        !
        if (tsum > safe_max ) then  ! thrsh_upper
           ierr=2 ! 3 
           write (out,"(' numerov:  no solution, tsum > thrsh_upper:',2g18.8 )") tsum,thrsh_upper
           return
        endif 
        !
        yc=y(ic)/pcout
        ycm1=y(ic-1)/pcout
        !
        ycp1 = 0
        !
        if (pcin>small_) then 
          ycp1=y(ic+1)/pcin
        endif
        !
        if (trim(boundary_condition)=='UNBOUND') exit
        !
        dx=((-ycm1+yc+yc-ycp1)/hh+(pot_eff(ic)-i0(ic)*eguess))*yc/tsum
        !
        eguess=eguess+dx
        !
        niter=niter+1
        !
        if (niter==iterlimit) then 
           ierr=2
           write (out,"(' numerov:  iteration limit reached in numerical integration (niter=',i8)") iterlimit
           write (out,"('           energy difference =',d13.4)") dx
           write (out,"('           convergence threshold =',d13.4,/)") thrsh_int
           return
        endif
        !if (eguess<enerlow) then 
        !   ierr=1
        !endif 
        !
     enddo
     if (eguess<enerlow) then 
        ierr=1
     endif 
     !
     iref = ic
     !
     if (trim(boundary_condition)=='UNBOUND') then 
       !
       !phi_f(0:ic)=phi_f(0:ic)/pcout*pcin
       !phi_f(ic+1:npoints)=phi_f(ic+1:npoints)/pcin
       !
    else
       phi_f(0:ic)=phi_f(0:ic)/pcout
       phi_f(ic+1:npoints)=phi_f(ic+1:npoints)/pcin
     endif
     !
  contains 
     !
     function y(i) result (v)
     integer(ik),intent(in) :: i
     real(rk) :: v
        !
        v=(1.0_rk-hh*(pot_eff(i)-i0(i)*eguess)/12.0_rk)*phi_f(i)
        !
     end function y
     !
  end subroutine numcoo

!
!     subroutine performs the outward integration of numcoo ,
!     until the first maximum in the wave function is found.
!     the sum i0(i)*phi_f(i)**2 is saved for the outward integration
!     and will later br divided by pc(out)**2.

  subroutine intout ( v, pot_eff, i0, eguess, phi_f, sumout,istart, iref, ic, pcout, ierr)
      !
     integer(ik),intent(in ) :: v,istart,iref
     real(rk),intent(in ) :: eguess
     integer(ik),intent(inout) :: ierr
     integer(ik),intent(out) :: ic


     real(rk),intent(in)  ::  i0(0:npoints),pot_eff(0:npoints)
     real(rk),intent(inout) :: phi_f(0:npoints)
     real(rk),intent(out) :: sumout,pcout
     !
     integer(ik) :: i,iend,imin_ref
     logical :: notfound
     !
     real(rk) :: hh,const,tsum,redfac,yi,yim1,yip1,phi_t
     !
     hh=rhostep*rhostep
     !
     !istart=1
     !
     yim1=phi_f(istart-1)*(1.0_rk-hh*(pot_eff(istart-1)-i0(istart-1)*eguess)/12.0_rk)
     !
     yi=phi_f(istart)*(1.0_rk-hh*(pot_eff(istart)-i0(istart)*eguess)/12.0_rk)
     !
     const=hh*(pot_eff(istart)-i0(istart)*eguess)
     !
     yip1=const*phi_f(istart)+yi+yi-yim1
     !
     tsum=0.0_rk
     !
     notfound = .true.
     !
     i = istart-1
     iend = npoints-2 
     if (v==0) iend = imin
     !
     ! if minimum is at i= 0 choose iend in the middle 
     !
     !imin_ref = imin
     imin_ref = iref
     !
     do while(notfound.and.i<=iend-1)
        !
        i = i + 1
        !
        const=hh*(pot_eff(i)-i0(i)*eguess)
        !
        phi_t = yi/(1.0_rk-const/12.0_rk)
        !
        phi_f(i)=phi_t
        !
        yip1=const*phi_t+yi+yi-yim1
        !
        if (i>=imin_ref) then 
          !
          if ( sign( 1.0_rk,yi-yim1 )/=sign(1.0_rk,yip1-yi).and.i.ne.1) then 
          !if ( sign( 1.0_rk,yim1 )/=sign(1.0_rk,yip1 ).and.i.ne.1) then 
              notfound = .false.
              cycle 
          endif 
          !
          if ( i==iref ) then 
          !if ( sign( 1.0_rk,yim1 )/=sign(1.0_rk,yip1 ).and.i.ne.1) then 
              notfound = .false.
              cycle 
          endif 
          !
        endif 
        !
        yim1=yi
        !
        yi=yip1
        !
        tsum=tsum+i0(i)*phi_f(i)*phi_f(i)
        !
        ! renormalizing phi_f
        !
        if (tsum > thrsh_upper) then 
            redfac=sqrt(thrsh_upper)
            !
            phi_f(0:i)=phi_f(0:i)/redfac
            yi=yi/redfac
            yim1=yim1/redfac
            tsum=tsum/thrsh_upper
        endif 
        !
     enddo
     !
     sumout=tsum
     !
     if (notfound) then 
        if (v/=0) ierr=1
        ic=imin_ref
        if (verbose>=7) then 
          write (out,"('intout: no turning point found in outward integration')")
        endif 
     else
        sumout=tsum
        ic=i-1
        pcout=phi_f(ic)
     endif


     pcout=phi_f(ic)


  end subroutine intout


!     subroutine performs the inward integration from rho=rhomax
!     stopping at the maximum in the wavefunction determined in
!     intout. the recursion relations given in numcoo are used
!     in conjunction with the two starting values :-
!
!     phi(rhomax) = 0
!
!     phi(rhomax-rhostep) = small_
!
!
  subroutine intin ( pot_eff, i0, eguess, phi_f , sumin , ic , pcin )
      
     real(rk),intent(in ) :: eguess
     integer(ik),intent(in) :: ic

     real(rk),intent(in)  ::  i0(0:npoints),pot_eff(0:npoints)
     real(rk),intent(inout) :: phi_f(0:npoints)
     real(rk),intent(out) :: sumin,pcin
     !
     real(rk) :: hh,const,tsum,yi,yim1,yip1,redfac,phi_t
     integer(ik) :: i,ist,nm1,iend,kend,km1
     !
     hh=rhostep**2
     !
     ist=ic+1
     iend=npoints-1
     nm1=iend
     !
     yip1=phi_f(nm1+1)*(1.0_rk-hh*(pot_eff(nm1+1)-i0(nm1+1)*eguess)/12.0_rk)
     yi  =phi_f(nm1  )*(1.0_rk-hh*(pot_eff(nm1  )-i0(nm1  )*eguess)/12.0_rk)
     !
     if (periodic.and.abs(phi_f(0))>sqrt(small_).and..false.) then 
       !
       kend = 2
       km1 = kend
       !
       yip1=phi_f(km1+1)*(1.0_rk-hh*(pot_eff(km1+1)-i0(km1+1)*eguess)/12.0_rk)
       yi  =phi_f(km1  )*(1.0_rk-hh*(pot_eff(km1  )-i0(km1  )*eguess)/12.0_rk)
       !
      ! do i=kend,0,-1
      !   !
      !   const=hh*(pot_eff(i)-i0(i)*eguess)
      !   !
      !   phi_t=yi/(1.0_rk-const/12.0_rk)
      !   !
      !   phi_f(i)=phi_t
      !   !
      !   yim1=const*phi_t+yi+yi-yip1
      !   !
      !   yip1=yi
      !   yi=yim1
      !   !
      ! enddo
       !
       !iend=npoints-2
       !nm1=iend
       !
       !phi_f(npoints-2) = 2.0_rk*phi_f(npoints)-phi_f(2)
       !
       !yip1=phi_f(nm1+1)*(1.0_rk-hh*(pot_eff(nm1+1)-i0(nm1+1)*eguess)/12.0_rk)
       !yi  =phi_f(nm1  )*(1.0_rk-hh*(pot_eff(nm1  )-i0(nm1  )*eguess)/12.0_rk)
       !
     endif
     !
     !yip1=0.0_rk
     !
     !yi=phi_f(nm1)*(1.0_rk-hh*(pot_eff(nm1)-i0(nm1)*eguess)/12.0_rk)
     !
     tsum=0.0_rk
     !
     do i=iend,ist,-1
        !i=ist+iend-ii
        const=hh*(pot_eff(i)-i0(i)*eguess)
        !
        phi_t=yi/(1.0_rk-const/12.0_rk)
        !
        phi_f(i)=phi_t
        !
        yim1=const*phi_t+yi+yi-yip1
        !
        yip1=yi
        yi=yim1
        !
        tsum=tsum+i0(i)*phi_f(i)*phi_f(i)
        !
        ! renormalizing phi_f
        !
        if (tsum > thrsh_upper) then 
            redfac=sqrt(thrsh_upper)
            !
            phi_f(i:iend)=phi_f(i:iend)/redfac
            yi=yi/redfac
            yip1=yip1/redfac
            tsum=tsum/thrsh_upper
        endif 
        !
     enddo
     !
     sumin=tsum
     !
     pcin=yi/(1.0_rk-hh*(pot_eff(ic)-i0(ic)*eguess)/12.0_rk)
     !
  end subroutine intin

!
! integration with Simpson rules 
!                                      
  function simpsonintegral_rk(npoints,xmax,f) result (si) 
    integer(ik),intent(in) :: npoints
    !
    real(rk),intent(in) :: xmax,f(0:npoints)
    !
    real(rk) :: si
    !
    integer(ik) :: i
    !
    real(rk) ::  feven,fodd,f0,fmax,h
      !
      h = xmax/real(Npoints,kind=rk)  !   integration step   
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
  function simpsonintegral(npoints,xmax,f) result (si) 
    integer(ik),intent(in) :: npoints
    !
    real(rk),intent(in) :: xmax,f(0:npoints)
    !
    real(rk) :: si
    !
    integer(ik) :: i
    !
    real(rk) ::  feven,fodd,f0,fmax,h
      !
      h = xmax/real(Npoints,kind=rk)  !   integration step   
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

  end function  simpsonintegral


! integration with Simpson rules 
!                                      
  function integral_rect_rk(npoints,xmax,f) result (si) 
    integer(ik),intent(in) :: npoints
    !
    real(rk),intent(in) :: xmax,f(0:npoints)
    !
    real(rk) :: si
    !
    real(rk) ::  h
     !
     h = xmax/real(Npoints,kind=rk)  !   integration step   
     !
     !  sum of odd and even contributions 
     !
     si = sum(f)*h
     !
     !si  = simpsonintegral_rk(npoints,xmax,f)
     !
  end function  integral_rect_rk



      subroutine d04aaf(ifun, nder, hbase, der, erest, fun, ifail , maxpts)
!
!      integer(ik) :: gmatrix( k1,k2 ),mbasp1,leniw
!
!     ***   purpose   ***
!
!     a subroutine for numerical differentiation at a point.  it
!     returns a set of approximations to the j-th order derivative
!     (j=1,2,...14) of fun(x) evaluated at x = xval  and, for each
!     derivative, an error estimate ( which includes the effect of
!     amplification of round- off errors).
!
!
!     ***   input  parameters   ***
!
!     (1) xval   real.  the abscissa at which the set of
!     derivatives is required.
!     (2) nder   integer(ik) ::. the highest order derivative required.
!     if(nder.gt.0)  all derivatives up to min(nder,14) are
!     calculated.
!     if(nder.lt.0 and nder even)  even order derivatives
!     up to min(-nder,14) are calculated.
!     if(nder.lt.0 and nder odd )  odd  order derivatives
!     up to min(-nder,13) are calculated.
!     (3) hbase  real.  a step length.
!     (6) fun    the name of a real*8 function subprogramme,
!     which is required by the routine as a subprogramme
!     and which represents the function being differentiated
!     the routine requires 21 function evaluations fun(x)
!     located at    x = xval    and at
!     x = xval + (2*j-1)*hbase,  j=-9,-8, .... +9,+10.
!     the function value at  x = xval is disregarded when
!     only odd order derivatives are required.
!
!     ***   output parameters   ***
!
!     (4) der(j)   j=1,2,...14.  real. a vector of
!     approximations to the j-th derivative of fun(x)
!     evaluated at x = xval.
!     (5) erest(j) j=1,2,...14.  real. a vector of
!     estimates of the dabsolute accuracy of der(j).  these
!     are negative when erest(j).gt.abs(der(j)), or when,
!     for some other reason the routine is doubtful about
!     the validity of the result.
!
!     ***   important warning   ***
!
!     erest is an estimate of the overall error.  it is essential
!     for
!     proper use that the user checks each value of der
!     subsequently
!     used to see that it is accurate enough for his purposes.
!     failure to do this may result in the contamination of all
!     subsequent results.
!     it is to be expected that in nearly all cases  der(14) will
!     be
!     unusable.( 14 represents a limit in which for the easiest
!     function likely to be encountered this routine might just
!     obtain
!     an approximation of the correct sign.)
!
!     ***   note on calculation   ***
!
!     the calculation is based on the extended t-table (t sub
!     k,p,s)
!     described in lyness and moler, num. math., 8, 458-464,(1966).
!     references in comment cards to ttab(nk,np,ns) refer to
!     t-table
!     with nk=k+1,  np=p+1,  ns=s+1.   since only part of the
!     extended
!     t-table is needed at one time, that part is stored in
!     rtab(10,7)
!     and is subsequently overwritten.  here
!     rtab(nk,np-ns+1)  =  ttab(nk,np,ns).
!
!     nag copyright 1976
!     mrk 5 release
!     mrk 6b revised  ier-116 (mar 1978)
!     mrk 7 revised ier-139 (dec 1978)
!     mrk 8 revised. ier-219 (mar 1980)
!     mrk 8d revised. ier-271 (dec 1980).
      real(rk),intent(in) ::  hbase
      integer(ik) ::  nderp, nder, ndpar, ifail, n, nn1, nsq, npar, nsmax,j, &
                     neger, ns, npmin, np, nktop, npr, nk, ntop, nmax, nmin,nup, nlo, & 
                     ntopm2, jns, nst2,ii,ifun
      real(rk) ::  der(4), erest(4), fact, &
       acof(10),fz, h, hfact, hsq, a, b, trange(7), ranmin, temp, term, &
       rtab(10,7), tzerom(10), erprev, ernow, xjns, xmax, xmin,xnum, &
       thron2, zero, one, two, big, small, htest, ftest
      
      integer(ik) :: maxpts

      real(rk) ::  fun(maxpts)
      data zero, one, two, thron2 /0.0_rk,1.0_rk,2.0_rk,1.5_rk/
!
!     to alter the precision of the whole routine,alter the
!     precision
!     in the above declaration and data statements
!
!     quit on zero hbase and zero nder.  set control parameter
!     ndpar.
!     ndpar = +1  odd-order derivatives only.
!     ndpar =  0  all derivatives
!     ndpar = -1  even-order derivatives only.
      big=safe_max
      small=safe_min
      if (hbase.eq.zero) go to 20
      nderp = nder
      ndpar = 0
      if (nder) 40, 20, 60
20    ifail=1
      return
40    ndpar = 4*(nder/2) - 2*nder - 1
      nderp = -nder
60    continue
      if (nderp.gt.4) nderp = 4
!
!     next, evaluate 21 function values , set acof(n), and set
!     first
!     column of rtab for odd order derivatives and store in tzerom
!     the
!     first column of rtab for even order derivatives.
!     fz = fun(xval)
!
      fz=fun(ifun)
      nst2=maxpts+maxpts
!
      do 80 n=1,10
       nn1 = 2*n - 1
       nsq = nn1*nn1
       acof(n) = nsq
       xnum = nn1
       h = hbase*xnum
!
       ii=ifun+nn1
       if (ii.gt.maxpts) temp=zero
!
       if (ii.le.maxpts) temp=fun(ii)
!
       ii=ifun-nn1
       if (ii.eq.0) term=zero
!
       if (ii.lt.0) term=-fun(iabs(ii))
!
       if (ii.gt.0) term=fun(ii)
!
       rtab(n,1) = (temp-term)/(two*h)
       tzerom(n) = (temp-fz-fz+term)/(two*h**2)
80    continue
!
!     set up for odd-order derivatives.
      if (ndpar.eq.-1) go to 100
      npar = 1
      nsmax = (nderp+1)/2
      if (nsmax.le.0) go to 380
      go to 140
!
!     set up for even-order derivatives.
100   continue
      npar = -1
      if (ndpar.eq.+1) go to 460
      nsmax = nderp/2
      if (nsmax.le.0) go to 460
!
!     set the first column of the t-table for even order
!     derivatives.
      do 120 j=1,10
       rtab(j,1) = tzerom(j)
120   continue
140   continue
!
!     odd-order and even order paths join up here.
!
      neger = 0
      hsq = hbase*hbase
      fact = one
      hfact = one
      do 360 ns=1,nsmax
!
!     from here on to   statement 360 (near end)  we are dealing
!     with a
!     specified value of ns.
!
!     for each value of np we calculate range(np)  which is the
!     differe
!     between the greatest and the least of the elements
!     ttab(nk,np,ns)
!     as we go along we determine also the minimum of range(np)
!     which i
!     ranmin = range(npmin).
!     we retain nup and nlo which are the values  of nk giving the
!     extreme values of ttab(nk,npmin,ns).
!     this part of ns loop concludes at statement number 280.
!
!
!     first calculate elements of t-table  for current ns value.
       npmin = ns
       if (ns.eq.1) npmin = 2
       do 180 np=npmin,7
          nktop = 10 - np + 1
          npr = np - ns + 1
          do 160 nk=1,nktop
             j = np + nk - 1
             a = acof(j)
             b = acof(nk)
             term = zero
             temp = zero
             if (np.ne.ns) temp = a*rtab(nk,npr-1) -b*rtab(nk+1,npr-1)
             if (ns.ne.1) term = -rtab(nk,npr) + rtab(nk+1,npr)
             rtab(nk,npr) = (term+temp)/(a-b)
160         continue
180      continue
!
!     now calculate nup,nlo and ranmin.
       do 280 np=ns,7
          ntop = 11 - np
          npr = np - ns + 1
          xmax = rtab(1,npr)
          xmin = rtab(1,npr)
          nmax = 1
          nmin = 1
!
          do 200 nk=2,ntop
             temp = rtab(nk,npr)
             if (temp.gt.xmax) nmax = nk
             if (temp.gt.xmax) xmax = temp
             if (temp.lt.xmin) nmin = nk
             if (temp.lt.xmin) xmin = temp
200         continue
!
          trange(np) = xmax - xmin
          if (np.ne.ns) go to 220
          ranmin = trange(np)
          go to 240
220         continue
          if (trange(np).ge.ranmin) go to 260
          ranmin = trange(np)
240         continue
          npmin = np
          nup = nmax
          nlo = nmin
          if (nlo.eq.nup) nlo = nlo + 1
260         continue
280      continue
!
!     next we take the average of all except the extreme values of
!     ttab(nk,npmin,ns) as the result and ranmin as the error
!     estimate.
       term = zero
       ntop = 11 - npmin
       ntopm2 = ntop - 2
       xnum = ntopm2
       j = npmin - ns + 1
       do 320 nk=1,ntop
          if (nk.eq.nup) go to 300
          if (nk.eq.nlo) go to 300
          term = term + rtab(nk,j)
300         continue
320      continue
       term = term/xnum
! 
!     the above result and error estimate refer to taylor
!     coefficient.
!     next we do scaling to obtain derivative results instead.
!     jns and xjns are actual order of derivative being treated.
!     safety factors 1.5,1.5,2.0,2.0 and 2.0 are multiplied into
!     the
!     error estimates for j = 10,11,12,13 and 14 respectively.
!     these
!     have been chosen by inspection of performance statistics.
!     there
!     is no analytic justification for these particular factors.
       jns = 2*ns - 1
       if (npar.lt.0) jns = 2*ns
       xjns = jns
       fact = fact*xjns
!     test for underflow of hfact
       if (hfact.ge.1.0e0) go to 330
       htest = hfact*big
       ftest = term*fact
       if ((two*ranmin*fact).gt.ftest) ftest = two*ranmin*fact
       if (htest.gt.ftest) go to 330
       der(jns) = big
       erest(jns) = -der(jns)
       neger = neger + 1
       go to 340
!     end of code to handle zero hfact
330      continue
       der(jns) = term*fact/hfact
       if (jns.eq.10) ranmin = thron2*ranmin
       if (jns.eq.11) ranmin = thron2*ranmin
       if (jns.ge.12) ranmin = two*ranmin
       erest(jns) = ranmin*fact/(hfact)
!
!     set sign of erest.  erest negative either if it swamps der,
!     or if
!     two previous consecutive erests of this parity are negative.
!     it may also be set negative at end (750).
       if (neger.ge.2) erest(jns) = -erest(jns)
       if (neger.ge.2) go to 340
       if (term.lt.ranmin) erest(jns) = -erest(jns)
       if (-term.ge.ranmin) erest(jns) = -erest(jns)
       if (erest(jns).lt.zero) neger = neger + 1
       if (erest(jns).ge.zero) neger = 0
340      continue
!
!     second test for underflow of hfact
       if (hfact.ge.1.0e0) go to 346
       if (hfact.eq.0.0e0) go to 352
       if (hsq.ge.small/hfact) go to 346
       hfact = 0.0e0
       go to 352
346      hfact = hsq*hfact
352      fact = fact*(xjns+one)
360   continue
!
380   continue
      if (npar.gt.0) go to 100
!
!     set sign of erest negative if two previous consecutive erests
!     are negative.
      if (ndpar.ne.0) go to 460
      if (nderp.le.2) go to 460
      neger = 0
      erprev = erest(1)
      do 440 j=2,nderp
       ernow = erest(j)
       if (neger.eq.2) go to 400
       if (erprev.ge.zero) go to 420
       if (ernow.ge.zero) go to 420
400      neger = 2
       if (ernow.gt.zero) erest(j) = -ernow
420      erprev = ernow
440   continue
!
460   continue
      ifail = 0
!
!     do 500 i=1,nderp
!     erest(i)=erest(i)/der(i)
! 500 continue
!
      return
!     end of d04aaf   ***   ***   ***   ***   ***
      end subroutine d04aaf

!




     !
     ! extrapolation at borders
     !
     subroutine interapolate_at_center(N,V,dv)
     !
     integer,intent(in) :: N
     real(rk),intent(inout) ::dv,V(-N:N)
     !
     integer            :: i1,i2,fact
     real(rk)           :: x1,a(N,N),b(N,1)
        !
        fact = 2
        if (abs(v(-2)/v(2)+1.0_rk)<sqrt(small_)) then 
           dv = 0 
           return 
        endif 
        !
        do i1 = 1,N
           !
           x1 = i1*rhostep
           !x2 = i1*rhostep
           !
           !
           b(i1,1) = V( i1)
           !
           do i2 = 2,N
             !
             a(i1 ,i2) = x1**((i2-1)*fact)
             !
           enddo
        enddo
        !
        !b(1,1) = V(0)
        !a(1,1) = 1.0_rk
        a(:,1) = 1.0_rk
        !a(1,2:N) = 0
        !
        !  lapack_gelss 
        ! 
        call lapack_gelss(a(:,:),b(:,:))
        !
        dv = b(1,1)

        !
   end subroutine interapolate_at_center


  subroutine diff_2d_4points_rk(Npoints,rho_b,f,periodic,reflect,d1f,d2f)

   integer(ik) ,intent(in)  :: Npoints
   real(rk)   ,intent(in)  :: f(0:Npoints)
   real(rk)   ,intent(in)  :: rho_b(2)
   logical     ,intent(in)  :: periodic
   integer(ik),intent(in)   :: reflect
   !
   real(rk),   intent(out) ::  d1f(0:Npoints)
   real(rk),   optional    ::  d2f(0:Npoints)
   !
   real(rk)               :: rhostep,d1t,d2t,x,dy
   integer(ik)             :: ipoint,i,kl,kr
   integer(ik),parameter   :: Nextrap = 4
   integer(ik),parameter   :: Mextrap = 4
   real(rk)                :: v_t(-Nextrap:Nextrap),g_t(1:Nextrap)
     !
     !
     rhostep = (rho_b(2)-rho_b(1))/real(npoints,kind=rk)
     !
     do ipoint = 0,npoints
        !
        kl = max(0,ipoint-4) ;   kr = min(Npoints,ipoint+4)
        !
        call ML_diffs_rk(kl,kr,ipoint,f(kl:kr),rhostep,d1t,d2t) 
        !
        d1f(ipoint) = d1t
        !
        if (present(d2f)) then
           !
           d2f(ipoint) = d2t
           !
        endif
        !
     enddo
     !
     do i = 1,Nextrap
       !
       ipoint = Mextrap-1+i
       v_t(i) = rho_b(1)+real(ipoint,rk)*rhostep
       ipoint = Npoints-Mextrap-Nextrap+i
       g_t(i) = rho_b(1)+real(ipoint,rk)*rhostep
       !
     enddo

     do i = 0,Mextrap-1
       !
       ipoint = i
       x = rho_b(1)+real(ipoint,rk)*rhostep
       call polintrk(v_t(1:Nextrap), d1f(Mextrap:Mextrap+Nextrap-1), x, d1t, dy)
       d1f(ipoint) = d1t
       !
       if (present(d2f)) & 
          call polintrk(v_t(1:Nextrap), d2f(Mextrap:Mextrap+Nextrap-1), x, d2f(ipoint), dy)
       !
       ipoint = Npoints-Mextrap+1+i
       x = rho_b(1)+real(ipoint,rk)*rhostep
       !
       call polintrk(g_t(1:Nextrap), d1f(Npoints-Mextrap-Nextrap+1:Npoints-Mextrap), x, d1t, dy)
       !
       d1f(ipoint) = d1t
       !
       if (present(d2f)) & 
          call polintrk(g_t(1:Nextrap), d2f(Npoints-Mextrap-Nextrap+1:Npoints-Mextrap), x, d2f(ipoint), dy)
       !
     enddo
     !
     if (periodic) then
       !
       do ipoint = 0,3
         !
         v_t(0:4) = f(ipoint:ipoint+4)
         !
         do i = -4,-1
           !
           kl = mod(npoints+ipoint+i,npoints)
           !
           v_t(i) = f(kl)
           !
           if (reflect/=0.and.kl>4) v_t(i) = v_t(i)*real(reflect,rk)
           !
         enddo
         !
         call ML_diffs_rk(-4,4,0,v_t(-4:4),rhostep,d1t,d2t) 
         !
         d1f(ipoint) = d1t
         !
         if (present(d2f)) then
            !
            d2f(ipoint) = d2t
            !
         endif
         !
       enddo
       !
       do ipoint = npoints-3,npoints
         !
         v_t(-4:0) = f(ipoint-4:ipoint)
         !
         do i = 1,4
           !
           kl = mod(ipoint+i,npoints)
           !
           v_t(i) = f(kl)
           !
           if (reflect/=0.and.kl<npoints-4) v_t(i) = v_t(i)*real(reflect,rk)
           !
         enddo
         !
         call ML_diffs_rk(-4,4,0,v_t(-4:4),rhostep,d1t,d2t) 
         !
         d1f(ipoint) = d1t
         !
         if (present(d2f)) then
            !
            d2f(ipoint) = d2t
            !
         endif
         !
       enddo
       !
     elseif (reflect/=0) then
       !
       do ipoint = npoints-3,npoints
         !
         v_t(-4:0) = f(ipoint-4:ipoint)
         !
         do i = 1,4
           !
           if (ipoint+i<=npoints) then
             kl = ipoint+i
             v_t(i) = f(kl)
           else
             !
             kl = mod(ipoint+i,npoints)
             v_t(i) = f(npoints-kl)*real(reflect,rk)
             !
           endif 
           !
         enddo
         !
         call ML_diffs_rk(-4,4,0,v_t(-4:4),rhostep,d1t,d2t) 
         !
         d1f(ipoint) = d1t
         !
         if (present(d2f)) then
            !
            d2f(ipoint) = d2t
            !
         endif
         !
       enddo
       !
     endif 
     !
  end subroutine diff_2d_4points_rk


  subroutine ML_diffs_rk(n1,n2,n0,v,rhostep,d1,d2)
    !
    integer(ik),intent(in) :: n1,n2,n0
    real(rk),intent(in) :: v(n1:n2),rhostep
    real(rk),intent(out) :: d1,d2
    !
    real(rk)             :: v_t(-4:4)
    !
    !n1 = lbound(v,dim=1) ; n2 = ubound(v,dim=1)
    !
    if (n1>n0.or.n0>n2) then 
       !
       write (out,"('ML_diffs_rk: must be n1<=n0<=n2, you give: ',3i8)") n1,n0,n2
       stop 'ML_diffs_rk: wrong indexes n1,n0,n2'
       !
    endif 
    !
    if ( n0-n1==4 .and. n2-n0==4 ) then 
       !
       v_t(-4:4) = v(n1:n2)
       !
       d1 = (-v_t( 2)/12.0_rk+2.0_rk/3.0_rk*v_t( 1) & 
             +v_t(-2)/12.0_rk-2.0_rk/3.0_rk*v_t(-1) )/rhostep
       !
       d2 = ( v_t( 4)/144.0_rk-v_t( 3)/9.0_rk+v_t( 2)*4.0_rk/9.0_rk+v_t( 1)/9.0_rk  & 
           -v_t(0)*65.0_rk/72.0_rk  &
           +v_t(-4)/144.0_rk-v_t(-3)/9.0_rk+v_t(-2)*4.0_rk/9.0_rk+v_t(-1)/9.0_rk) &
           /rhostep**2
       !
    elseif ( n0-n1>=2 .and. n2-n0>=2 ) then
       !
       v_t(-2:2) = v(n0-2:n0+2)
       !
       d1 = (-v_t( 2)/12.0_rk+2.0_rk/3.0_rk*v_t( 1) & 
             +v_t(-2)/12.0_rk-2.0_rk/3.0_rk*v_t(-1) )/rhostep
       !
       d2 = (v_t( 2) + v_t(-2) - 2.0_rk*v_t(0) )/(rhostep**2*4.0_rk)
       !
    elseif ( n0-n1>=1 .and. n2-n0>=1 ) then
       !
       v_t(-1:1) = v(n0-1:n0+1)
       !
       d1 = ( v_t( 1) - v_t(-1) )/(rhostep*2.0_rk)
       !
       d2 = (v_t( 1) + v_t(-1) - 2.0_rk*v_t(0) )/(rhostep**2)
       !
    elseif ( n0-n1==0 .and. n2-n0>=2 ) then
       !
       v_t(0:2) = v(n0:n0+2)
       !
       d1 = ( v_t( 1) - v_t(0) )/rhostep
       !
       d2 = (v_t( 0) - 2.0_rk*v_t( 1) + v_t( 2) )/(rhostep**2)
       !
    elseif ( n0-n1>=2 .and. n2-n0==0 ) then
       !
       v_t(-2:0) = v(n0-2:n0)
       !
       d1 = ( v_t(0) - v_t(-1) )/rhostep
       !
       d2 = ( v_t( 0) - 2.0_rk*v_t(-1) + v_t(-2) )/(rhostep**2)
       !
    else
       !
       write (out,"('ML_diffs_rk: wrong n1,n0,n2: ',3i8)") n1,n0,n2
       stop 'ML_diffs_rk: wrong indexes n1,n0,n2'
       !
    endif
    !
  end subroutine ML_diffs_rk


recursive subroutine polintrk(xa, ya, x, y, dy)

  !use nrtype; use nrutil, only : assert_eq, iminloc, nrerror
  implicit none
  real(rk), dimension(:), intent(in) :: xa, ya
  real(rk), intent(in) :: x
  real(rk), intent(out) :: y, dy

  !
  ! given arrays xa and ya of length n, and given a value x, this routine
  ! returns a value y, and an error estimate dy. if p(x) is the polynomial
  ! of degree n - 1 such that p(x a_i) = y a_i, i = 1, ..., n, then the
  ! returned value y = p(x).

  integer(ik) :: m, n, ns
  real(rk), dimension(size(xa)) :: c, d, den, ho

  if (verbose>=6) write(out,"('polintrk...')")

  n = size(xa)
  if (size(ya)/=n) then 
      stop  'polint: wrong sizes'
  endif 

  !n = assert_eq(size(xa), size(ya), 'polint')
  c = ya                          ! initialize the tableau of c's and d's
  d = ya
  ho = xa - x
  ns = minloc(abs(x - xa),dim=1)  ! find index ns of closest table entry
  y = ya(ns)                      ! initial approximation to y.
  ns = ns - 1
  do m = 1, n - 1                        ! for each column of the tableau
     den(1:n-m) = ho(1:n-m) - ho(1+m:n)  ! we loop over c's and d's and
     !if (any(den(1:n-m) == 0.0)) &       ! update them
     !     call nrerror('polint: calculation failure')

     if (any(den(1:n-m) == 0.0)) then                   ! interpolating function
          write(out,"('failure in polint, has a pole here')")
          stop 'failure in polint'
          !call nrerror('failure in polint')         ! has a pole here
     endif 
     ! this error can occur only if two input xa's are (to within roundoff)
     ! identical.

     den(1:n - m) = (c(2:n-m+1) - d(1:n-m))/den(1:n-m)
     d(1:n-m) = ho(1+m:n) * den(1:n-m)   ! here c's and d's get updated
     c(1:n-m) = ho(1:n-m) * den(1:n-m)
     if (2 * ns < n-m) then       ! after each column in the tableau is
        dy=c(ns+1)                ! completed decide, which correction
     else                         ! c or d we add to y. we take the
        dy=d(ns)                  ! straightest line through the tableau
        ns=ns-1                   ! to its apex. the partial approximations
     end if                       ! are thus centred on x. the last dy
     y = y+dy                     ! is the measure of error.
  end do
  !
  if (verbose>=6) write(out,"('polintrk.')")
  !

end subroutine polintrk




      SUBROUTINE OVRLAP(BFCT,DER,EFN,OVR,OVRCRT,PSI,RH,RMIN,TMFPRM,VJ,&
                        VLIM,z,ifs,IWR,JP,NEND,OTMF,TMFTYP)
!=======================================================================
!        Routine by R.J. Le Roy;  Last Modified 19 August 2003
!=======================================================================
! Calculate overlap integral Franck-Condon Moment FCM(i) array 
! between the given bound state wave function PSI(i) (which is zero
! for  i > NEND) and the J' = JP continuum final state wave function
! (asymptotocally normalized to unit amplitude) at energy EFN on the
! effective potential VJ(I) with asymptote VLIM, with input array
! z(i).  z(i) is the input (radial) array whose moments are being taken.
! NOTE that z(i) = zin(i) = ztmf
!
! On entry, energy units for EFN and VLIM are (cm-1), while VJ(I)
! incorporates the factor BFCT (i.e., VJ/BFCT has units cm-1).
!
! Convergence of asymptotic wave function normalization defined by
! requirement that W.K.B. fits to 3 successive maxima must agree
! relatively to within OVRCRT. 
!-----------------------------------------------------------------------
!cc   INCLUDE 'arrsizes.h'
!-----------------------------------------------------------------------
!  Utility routine to summarize dimensioning of arrays
!-----------------------------------------------------------------------
      INTEGER mxdata,mxisp,mxfsp,mxnj,mxnp,mxntp,mxprm,mxv,mxfs,mxisot,&
              mxsets,mxfreq
      double precision :: CCM,PI
!-----------------------------------------------------------------------
!  mxdata - maximum number of input data points
!  mxisp  - maximum number of points for initial state potential array
!           (also used for number of points in transition moment array)
!  mxnj    - maxiumum value of j quantum number allowed
!  mxfsp  - maximum number of points for final state potential array
!  mxnp   - maximum number of parameters total
!  mxntp  - maximum number of turning points to be read in
!  mxprm  - maximum number of parameters for final state pot'l or TMF
!  mxv    - largest value for the v quantum number
!  mxfs   - maximum number of final states allowed
!  mxisot - maximum number of isotopomers allowed
!  mxsets - maximum number of data sets allowed
!  mxfreq - maximum number of data points allowed in a given set
!-----------------------------------------------------------------------
      PARAMETER (mxisp=16001)
      PARAMETER (mxnj=20)
      PARAMETER (mxfsp=16001)
      PARAMETER (mxntp=9999)
      PARAMETER (mxprm=6)
      PARAMETER (mxv=200)
      PARAMETER (mxfs=5)
      PARAMETER (mxisot=3)
      PARAMETER (mxsets=11)
      PARAMETER (mxfreq=501)
      PARAMETER (PI=3.141592653589793238d0)
      PARAMETER (mxnp=2*mxprm*mxfs+mxsets-1)
      PARAMETER (mxdata=mxfreq*mxsets)
      PARAMETER (CCM= 299792458d2)
      INTEGER i,ifs,MESH1,MESH2,MESH3,step,first,last,IWR,JP,TURNPT,m,&
              OTMF,NAMP,NEND,TMFTYP(mxfs)
      REAL(rk) ::  AMP1,AMP2,AMP3,AMP4,BFCT,DER(0:mxprm-1),DI,FCFACT,&
             EFN,ELIM,ER,FCM(0:MXPRM-1),EDIFF1,EDIFF2,EDIFFi,HALF,HARG,&
             NFACT,&
             RH,RMIN,OVR,OVRCRT,PSI(NEND),S0,S1,S2,SG1,SG2,SGi,Si,&
             SNARG,SQKINF,&
             THIRD,TMFPRM(0:mxprm-1,mxfs),VLIM,VJ(mxfsp),VV,XIITH,XX,&
             Y1,Y2,Y3,z(mxisp),&
             ZTST,ZZ0,ZZ1
!-----------------------------------------------------------------------
      HALF=  1.D0/2.D0
      XIITH= 1.D0/12.D0
      THIRD= 1.D0/3.D0
      ER= EFN*BFCT
      ELIM= VLIM*BFCT
      SQKINF= DSQRT(ER-ELIM)
      AMP1= 1.D0
      AMP2= 2.D0
      AMP3= 0.D0
      AMP4= 0.D0
      DO m= 0,OTMF
         DER(m)= 0.D0
         FCM(m)= 0.d0
      ENDDO
!-----------------------------------------------------------------------
!** Locate first turning point and use Airy function to estimate
!  appropriate integration starting point such that  PSI(1) .LE. 1.D-10
!-----------------------------------------------------------------------
      MESH1= 1
      EDIFF1= VJ(MESH1)-ER
      step= DINT(0.05d0/RH)
      IF(step.LT.1) step= 1
      first= step+1
      DO i= first,mxfsp,step
         TURNPT= i
         EDIFF2= VJ(i)-ER
         IF(EDIFF2.LE. 0.D0) GOTO 4
         MESH1= i
         EDIFF1= EDIFF2
      ENDDO     
      IF(IWR.NE.0) THEN
         WRITE(6,607) JP,EFN
         OVR= 0.d0
         RETURN
      ENDIF
    4 MESH2= TURNPT
      TURNPT= MESH1+(MESH2-MESH1)*EDIFF1/(EDIFF1-EDIFF2)
      IF(IABS(TURNPT-MESH2).LE.1) GOTO 6
      IF((TURNPT.LE.0).OR.(TURNPT.GT.mxfsp)) THEN
         IF(IWR.NE.0) WRITE(6,601) JP,EFN
         OVR= 0.d0
         RETURN
      ENDIF
      MESH1= MESH2
      EDIFF1= EDIFF2
      EDIFF2= VJ(TURNPT)-ER
      GOTO 4
    6 DI= 10.D0/(VJ(TURNPT-1)-VJ(TURNPT))**THIRD
      step= DINT(DI)
      MESH1= MAX0(1,TURNPT-step)
      IF(MESH1.GE.NEND) THEN
         OVR= 0.D0
         RETURN
      ENDIF
    8 EDIFF1= VJ(MESH1)-ER
      IF(EDIFF1.LT.10.D0) GOTO 10
!-----------------------------------------------------------------------
!** Adjust starting point outward to ensure integration scheme stability
!-----------------------------------------------------------------------
      MESH1= MESH1+1
      IF((MESH1-mxfsp).LT.0) GOTO 8
      IF((MESH1-mxfsp).GE.0) THEN
         OVR= 0.D0
         RETURN
      ENDIF
   10 MESH2= MESH1+1
      MESH3= MESH2+1
!-----------------------------------------------------------------------
!** WKB starting condition for wave function
!-----------------------------------------------------------------------
      S0= 1.D0
      EDIFF1= VJ(MESH1)-ER
      EDIFFi= VJ(MESH2)-ER
      IF((EDIFF1.GT. 0.D0).AND.(EDIFFi.GT. 0.D0)) THEN
         SG1= DSQRT(EDIFF1)
         SGi= DSQRT(EDIFFi)
         Si= S0*DSQRT(SG1/SGi)*DEXP((SG1+SGi)/2.D0)
         IF(Si.LE.S0) S0= 0.D0
      ELSE
         VV= VJ(MESH2)/BFCT
         IF(IWR.NE.0) WRITE(6,608) JP,EFN,VV,MESH2
         S0= 0.D0
         Si= 1.D0
      ENDIF   
!-----------------------------------------------------------------------
!  notationally speaking, all Yi's refer to values used in the Numerov
!  Algorithm for wavefunciont propagation.
!-----------------------------------------------------------------------
      Y1= S0*(1.D0-XIITH*EDIFF1)
      Y2= Si*(1.D0-XIITH*EDIFFi)
!-----------------------------------------------------------------------
!  Use trapezoid rule for numerical integration.  Initialize FCM(m)
!  values using first section of area.   
!-----------------------------------------------------------------------
      ZZ0= 1.D0
      ZZ1= 1.D0
      DO m= 0,OTMF
         FCM(m)= HALF*S0*PSI(MESH1)*ZZ0 + Si*PSI(MESH2)*ZZ1
         ZZ0= ZZ0*z(MESH1)
         ZZ1= ZZ1*z(MESH2)
      ENDDO
      S2= S0
!-----------------------------------------------------------------------
!** Integrate outward to first turning point.  NOTE that Airy-estimated
!  initialization minimizes need for renormalizations.
!-----------------------------------------------------------------------
      DO 16 i= MESH3,TURNPT
         Y3= Y2+Y2-Y1+EDIFFi*Si
         Y1= Y2
         Y2= Y3
         EDIFFi= VJ(I)-ER
         S1= S2
         S2= Si
         Si= Y3/(1.D0-XIITH*EDIFFi)
!-----------------------------------------------------------------------
!** If bound wavefx. non-negligible, accumulate overlap moments
!  NOTE that FCFACT is the Franck-Condon factor
!-----------------------------------------------------------------------
         IF(I.LE.NEND) THEN
            FCFACT= Si*PSI(i)
            DO m= 0,OTMF
               FCM(m)= FCM(m) + FCFACT
               FCFACT= FCFACT*z(i)
            ENDDO
         ENDIF
!-----------------------------------------------------------------------
!** If wavefuntion too large in forbidden region, renormalize it ...
!-----------------------------------------------------------------------
         IF((Si.GE.1.D32).OR.(i.EQ.TURNPT)) THEN
            NFACT= 1.D0/Si
            Si= 1.D0
            IF(S0.GT.1.D-30) S0= S0*NFACT
            DO m= 0,OTMF
               FCM(m)= FCM(m)*NFACT
            ENDDO
            Y1= Y1*NFACT
            Y2= Y2*NFACT
         ENDIF
   16 CONTINUE
      IF((IWR.NE.0).AND.(S0/SI.GT.1.D-8))&
                                 WRITE(6,602)JP,EFN,MESH1,S0/Si,TURNPT
      MESH2= TURNPT+1
!-----------------------------------------------------------------------
!** If turning point NOT past end of range for bound state wavefx., then
!   integrate from turning point to end of bound-state wave function
!-----------------------------------------------------------------------
      IF(TURNPT.LT.NEND) THEN
         DO i= MESH2,NEND
            Y3= Y2 + Y2 - Y1 + EDIFFi*Si
            Y1= Y2
            Y2= Y3
            EDIFFi= VJ(i)-ER
            S1= S2
            S2= Si
            Si= Y3/(1.D0 - XIITH*EDIFFi)
            FCFACT= Si*PSI(i)
            DO m= 0,OTMF
               FCM(m)= FCM(m) + FCFACT
               FCFACT= FCFACT*z(i)
            ENDDO
         ENDDO   
         MESH2= NEND+1
      ENDIF
!-----------------------------------------------------------------------
!** Continue wave function propagation until amplitude converges
!-----------------------------------------------------------------------
      NAMP= 0
      DO 30 i= MESH2,mxfsp
         Y3= Y2 + Y2 - Y1 + EDIFFi*Si
         Y1= Y2
         Y2= Y3
         EDIFF2= EDIFFi
         EDIFFi= VJ(i)-ER
         S1= S2
         S2= Si
         Si= Y3/(1.D0-XIITH*EDIFFi)
         IF((Si.GE.S2).OR.(S1.GT.S2)) GOTO 30
!-----------------------------------------------------------------------
!** At successive maxima, fit solution to W.K.B. form to determine
!  apparent asymptotic amplitude.
!-----------------------------------------------------------------------
         SG2= DSQRT(-EDIFF2)
         SGi= DSQRT(-EDIFFi)
         HARG= (SG2 + SGi)/2.D0
         SNARG= 1.D0/DSQRT(1.D0 + ((DSQRT(SG2/SGi)*S2/Si- DCOS(HARG))&
                /DSIN(HARG))**2)
         NAMP= NAMP+1
         AMP4= AMP3
         AMP3= AMP2
         AMP2= AMP1
         AMP1= Si*DSQRT(SGi/SQKINF)/SNARG
         XX= RMIN + (i-1)*RH
         IF(IWR.GE.3) WRITE(6,604) JP,EFN,XX,AMP1
         last= i
!-----------------------------------------------------------------------
!** Test successive amplitudes for convergence
!-----------------------------------------------------------------------
         ZTST= OVRCRT*AMP1
         IF((DABS(AMP1-AMP2).LT.ZTST).AND.(DABS(AMP2-AMP3).LT.ZTST))  GOTO 35
   30 CONTINUE
      IF(IWR.NE.0) WRITE(6,603) JP,EFN,AMP1,AMP2,AMP3,AMP4
   35 OVR= 0.d0
      DO m= 0,OTMF
          FCM(m)= FCM(m)*RH/AMP1
          OVR= OVR + TMFPRM(m,ifs)*FCM(m)
          ENDDO
      IF(IWR.GE.1) WRITE(6,605) JP,EFN,last,XX,(m,FCM(m),m=0,OTMF)
      IF(IWR.GE.2) WRITE(6,606) S0,AMP1,AMP2,AMP3,AMP4
      DO m= 0,OTMF
          DER(m)= 2.d0*OVR*FCM(m)
          ENDDO
      OVR= OVR*OVR
      RETURN
!-----------------------------------------------------------------------
  601 FORMAT(' *** OVRLAP BOMBed ***  For   J = ',I3,'   EFN = ',&
        F10.2,'   never got to first turning point')
  602 FORMAT(' ** WARNING **  For   J = ', I3 ,'   EFN = ',F10.2,&
       '  starting wavefunction is PSI(',I4,') = ',D10.3,' and  I(turn.pt.) = ',I4)
  603 FORMAT(' ** WARNING ** For J= ',I3,' EFN= ',F9.2,' amplitude not',&
       ' converged by end of range'/3x,'Last four values are ', 4(1PD14.6))
  604 FORMAT(' At  J=',I3,'   E=',F9.2,'   R=',F6.3,'  apparent asymptotic amplitude',1PD14.6)
  605 FORMAT(' At  J=',I3,'   E= ',F12.2,'   R(end)= R(',I5,')=',F7.4,'    FCM(',I1,')=',F12.8:/(4x,3(5x,'FCM(',I1,')=',F12.8:)))
  606 FORMAT(5X,'S0= ',1PD10.3,'   & last 4 amplitudes are',2D14.6/45x,2D14.6)
  607 FORMAT(' *** ERROR ***   At   J = ',I3,'   EFN = ',F10.2,'  have  V .GT. E  everywhere.' )
  608 FORMAT(' *** Caution ***  For   J = ',I3,'   (EFN= ',F10.2,') .GE. (V = ',F10.2,')   at  I = ',I4, &
           &' ,so initialize with a node.')
  
      END subroutine OVRLAP




 end module me_numer


