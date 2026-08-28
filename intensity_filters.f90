module intensity_filters

  use accuracy,      only : ik, rk, small_
  use diatom_module, only : quantaT, fitting

  implicit none
  private

  public :: transitions_filter_from_fitting

contains

  subroutine transitions_filter_from_fitting(jI,jF,isymI,isymF,ilevelI,ilevelF,EnergyI,EnergyF,&
                                               quantaI,quantaF,passed,iEntry_found)

    real(rk),intent(in) :: jI,jF,EnergyF,EnergyI
    integer(ik),intent(in) :: isymI,isymF,ilevelI,ilevelF
    type(quantaT),pointer,intent(in) :: quantaI,quantaF
    logical,intent(inout) :: passed
    integer(ik),optional :: iEntry_found
    integer(ik) :: ientry,Nentries,N,N_,itau,itau_,vI,vF
    integer(ik) :: istateI,istateF,ilambdaI,ilambdaF
    real(rk) :: omegaI,omegaF,sigmaI,sigmaF
    real(rk) :: jrot,jrot_,lock_factor,nu
    integer(ik) :: maxiter_as = 3,iter_th

    ! No need to check further if an earlier transition filter rejected the
    ! transition.  PASSED is therefore an input/output argument, rather than
    ! INTENT(OUT) as it was in dipole.f90.
    if (.not.passed) return

    if (present(iEntry_found)) iEntry_found = -1

    vI       = quantaI%v
    iStateI  = quantaI%iState
    sigmaI   = quantaI%sigma
    ilambdaI = quantaI%ilambda
    omegaI   = quantaI%omega

    vF       = quantaF%v
    iStateF  = quantaF%iState
    sigmaF   = quantaF%sigma
    ilambdaF = quantaF%ilambda
    omegaF   = quantaF%omega

    Nentries = fitting%Nenergies

    loop_freq : do ientry = 1,Nentries

      jrot = fitting%obs(ientry)%jrot
      itau = fitting%obs(ientry)%iparity+1
      nu = fitting%obs(ientry)%energy

      if (nint(jF-jrot)/=0) cycle
      if (isymF/=itau) cycle

      jrot_ = fitting%obs(ientry)%jrot_
      itau_ = fitting%obs(ientry)%iparity_+1

      if (nint(jI-jrot_)/=0) cycle
      if (isymI/=itau_) cycle

      N  = fitting%obs(ientry)%N
      N_ = fitting%obs(ientry)%N_

      if (abs(fitting%threshold_lock)<sqrt(small_)) then

        if (N==iLevelF .and. N_==iLevelI) then
          if (present(iEntry_found)) iEntry_found = ientry
          return
        endif

      else

        if (fitting%threshold_lock<0) maxiter_as = 1

        loop_thresh_u : do iter_th = 1,maxiter_as

          lock_factor = real(iter_th,rk)

          if (abs(nu-(EnergyF-EnergyI))<=abs(lock_factor*fitting%threshold_lock).and.&
              (((istateF==fitting%obs(ientry)%quanta%istate.and.&
                 abs(ilambdaF)==abs(fitting%obs(ientry)%quanta%ilambda).and.&
                 vF==fitting%obs(ientry)%quanta%v.and.&
                 nint(abs(sigmaF)-abs(fitting%obs(ientry)%quanta%sigma))==0.and.&
                 nint(abs(omegaF)-abs(fitting%obs(ientry)%quanta%omega))==0.and.&
                 istateI==fitting%obs(ientry)%quanta_%istate.and.&
                 abs(ilambdaI)==abs(fitting%obs(ientry)%quanta_%ilambda).and.&
                 vI==fitting%obs(ientry)%quanta_%v.and.&
                 nint(abs(sigmaI)-abs(fitting%obs(ientry)%quanta_%sigma))==0.and.&
                 nint(abs(omegaI)-abs(fitting%obs(ientry)%quanta_%omega))==0).or.&
                fitting%threshold_lock<0))) then

            fitting%obs(ientry)%N  = ilevelF
            fitting%obs(ientry)%N_ = ilevelI
            if (present(iEntry_found)) iEntry_found = ientry
            return

          endif

        enddo loop_thresh_u

      endif

    enddo loop_freq

    passed = .false.

  end subroutine transitions_filter_from_fitting

end module intensity_filters
