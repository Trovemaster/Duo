module magnetic_dipole

  ! Operator-specific magnetic-dipole matrix elements used by the shared
  ! rank-1 intensity pipeline in dipole.f90.

  use accuracy, only : ik, rk, g_s, safe_min
  use diatom_module, only : Intensity, basis, poten, three_j, &
                            nMagneticDipoles, magnetictm, &
                            nMagneticRotDipoles, magnetrot

  implicit none
  private

  public :: md_do_1st_half_linestrength

contains

  subroutine md_do_1st_half_linestrength(jI,jF,indI,indF,dimenI,dimenF,vector,half_ls)

    implicit none

    real(rk),intent(in)     :: jI,jF
    integer(ik),intent(in)  :: indI,indF,dimenI,dimenF
    real(rk),intent(in)     :: vector(:)
    real(rk),intent(out)    :: half_ls(:)
    integer(ik)             :: icontrF,icontrI,ivibF,ivibI,idip,istateI,istateF,ilambdaF,ilambdaI
    integer(ik)             :: ipermute,istateI_,ilambdaI_,ilambdaF_,isigmav,iomegaI_,istateF_,itau
    real(rk)                :: ls,f3j,omegaI,omegaF,sigmaF,sigmaI,spinF,spinI,f6j
    real(rk)                :: spinI_,spinF_,f_t,g00,g20

    half_ls = 0.0_rk

    ! Parallelism is owned by the shared outer initial-state loop. Keeping
    ! this basis loop serial avoids nested OpenMP regions and timer races.
    loop_F : do icontrF = 1,dimenF

      ivibF = basis(indF)%icontr(icontrF)%ivib
      istateF = basis(indF)%icontr(icontrF)%istate
      omegaF = basis(indF)%icontr(icontrF)%omega
      sigmaF = basis(indF)%icontr(icontrF)%sigma
      spinF = basis(indF)%icontr(icontrF)%spin
      ilambdaF = basis(indF)%icontr(icontrF)%ilambda

      loop_I : do icontrI = 1,dimenI

        ivibI = basis(indI)%icontr(icontrI)%ivib
        istateI = basis(indI)%icontr(icontrI)%istate
        omegaI = basis(indI)%icontr(icontrI)%omega
        sigmaI = basis(indI)%icontr(icontrI)%sigma
        spinI = basis(indI)%icontr(icontrI)%spin
        ilambdaI = basis(indI)%icontr(icontrI)%ilambda

        if (nint(spinI-spinF)/=0.or.abs(ilambdaI-ilambdaF)>1) cycle loop_I

        iomegaI_ = int(omegaI)
        if (mod(nint(2.0_rk*omegaI+1.0_rk),2)==0) &
          iomegaI_ = nint((2.0_rk*omegaI-1.0_rk)*0.5_rk)

        f3j = three_j(jI,1.0_rk,jF,omegaI,omegaF-omegaI,-omegaF)

        ! Spin and orbital magnetic moments have a rank-1 rotational factor.
        if (abs(f3j)>intensity%threshold%coeff) then

          f_t = 0.0_rk

          if (ilambdaI==ilambdaF.and.istateI==istateF.and.ivibI==ivibF) then
            if (nint(2.0_rk*sigmaF)==nint(2.0_rk*(sigmaI+1.0_rk))) then
              f_t = g_s/sqrt(2.0_rk)*sqrt(spinI*(spinI+1.0_rk)-sigmaI*(sigmaI+1.0_rk))
            elseif (nint(2.0_rk*sigmaF)==nint(2.0_rk*(sigmaI-1.0_rk))) then
              f_t =-g_s/sqrt(2.0_rk)*sqrt(spinI*(spinI+1.0_rk)-sigmaI*(sigmaI-1.0_rk))
            elseif (nint(2.0_rk*sigmaF)==nint(2.0_rk*sigmaI)) then
              f_t = real(ilambdaI,rk)+g_s*sigmaI
            else
              cycle loop_I
            endif

            ls = f_t*f3j*vector(icontrI)
            half_ls(icontrF) = half_ls(icontrF)+(-1.0_rk)**iomegaI_*ls
          endif

          loop_iorbdipole : do idip = 1,nMagneticDipoles
            do ipermute = 0,1

              if (ipermute==0) then
                istateI_ = magnetictm(idip)%istate
                ilambdaI_ = magnetictm(idip)%lambda
                spinI_ = magnetictm(idip)%spini
                istateF_ = magnetictm(idip)%jstate
                ilambdaF_ = magnetictm(idip)%lambdaj
                spinF_ = magnetictm(idip)%spinj
              else
                istateF_ = magnetictm(idip)%istate
                ilambdaF_ = magnetictm(idip)%lambda
                spinF_ = magnetictm(idip)%spini
                istateI_ = magnetictm(idip)%jstate
                ilambdaI_ = magnetictm(idip)%lambdaj
                spinI_ = magnetictm(idip)%spinj
              endif

              if (ipermute==1.and.istateI_==istateF_.and.ilambdaI_==ilambdaF_.and.&
                  nint(spinI_-spinF_)==0) cycle
              if (istateI/=istateI_.or.istateF/=istateF_) cycle

              do isigmav = 0,1
                if (isigmav==1.and.abs(magnetictm(idip)%lambda)+&
                    abs(magnetictm(idip)%lambdaj)==0) cycle

                ilambdaI_ = ilambdaI_*(-1)**isigmav
                ilambdaF_ = ilambdaF_*(-1)**isigmav

                if (ilambdaI_/=ilambdaI.or.ilambdaF_/=ilambdaF) cycle
                if (abs(ilambdaI-ilambdaF)/=1) cycle

                f_t = (ilambdaI_-ilambdaF_)/sqrt(2.0_rk)*magnetictm(idip)%matelem(ivibI,ivibF)

                if (isigmav==1) then
                  itau = 0
                  if (ilambdaI_==0.and.poten(istateI)%parity%pm==-1) itau = itau+1
                  if (ilambdaF_==0.and.poten(istateF)%parity%pm==-1) itau = itau+1
                  f_t = f_t*(-1.0_rk)**itau
                endif

                ls = f_t*f3j*vector(icontrI)
                half_ls(icontrF) = half_ls(icontrF)+(-1.0_rk)**iomegaI_*ls
              enddo
            enddo
          enddo loop_iorbdipole
        endif

        if (ilambdaI==ilambdaF.and.istateI==istateF.and.ivibI==ivibF) then
          loop_magnetrot : do idip = 1,nMagneticRotDipoles

            istateI_ = magnetrot(idip)%istate
            ilambdaI_ = magnetrot(idip)%lambda
            spinI_ = magnetrot(idip)%spini

            if (istateI/=istateI_.or.ilambdaI/=ilambdaI_) cycle

            f_t = magnetrot(idip)%matelem(ivibI,ivibF)
            ls = 0.0_rk

            if (nint(omegaF-omegaI)==0) then
              g00 = -2.0_rk/sqrt(3.0_rk)*f_t
              f6j = Wigner6j(0.0_rk,1.0_rk,1.0_rk,jF,jI,jF)
              f3j = three_j(jI,0.0_rk,jF,-omegaI,0.0_rk,omegaF)
              ls = g00*f6j*f3j

              g20 = -2.0_rk/sqrt(6.0_rk)*f_t
              f6j = Wigner6j(2.0_rk,1.0_rk,1.0_rk,jF,jI,jF)
              f3j = three_j(jI,2.0_rk,jF,-omegaI,0.0_rk,omegaF)
              ls = ls+sqrt(5.0_rk)*g20*f6j*f3j
            endif

            half_ls(icontrF) = half_ls(icontrF)-(-1.0_rk)**iomegaI_*ls
          enddo loop_magnetrot
        endif
      enddo loop_I
    enddo loop_F

  end subroutine md_do_1st_half_linestrength


  function Wigner6j(j1,j2,j3,J_1,J_2,J_3) result(w6j)

    implicit none

    real(rk),intent(in) :: j1,j2,j3,J_1,J_2,J_3
    real(rk) :: w6j,threshold,tri1,tri2,tri3,tri4
    real(rk) :: a1,a2,a3,a4,k1,k2,k3,upper_limit,lower_limit,t

    tri1 = triangle_coefficient(j1,j2,j3)
    tri2 = triangle_coefficient(j1,J_2,J_3)
    tri3 = triangle_coefficient(J_1,j2,J_3)
    tri4 = triangle_coefficient(J_1,J_2,j3)
    threshold = safe_min

    if (abs(tri1)<threshold.or.abs(tri2)<threshold.or.&
        abs(tri3)<threshold.or.abs(tri4)<threshold) then
      w6j = 0.0_rk
      return
    endif

    a1 = j1+j2+j3
    a2 = j1+J_2+J_3
    a3 = J_1+j2+J_3
    a4 = J_1+J_2+j3

    k1 = j1+j2+J_1+J_2
    k2 = j2+j3+J_2+J_3
    k3 = j3+j1+J_3+J_1

    lower_limit = max(a1,a2,a3,a4)
    upper_limit = min(k1,k3,k3)

    t = lower_limit
    threshold = 1.0e-6_rk
    w6j = 0.0_rk

    do while (t<upper_limit+threshold)
      w6j = w6j+get_sign(t)*exp(log_())
      t = t+1.0_rk
    enddo

    w6j = sqrt(tri1*tri2*tri3*tri4)*w6j

  contains

    real(rk) function log_()
      real(rk) :: log_gamma

      log_ = log_gamma(t+2.0_rk) &
        -log_gamma(t-j1-j2-j3+1.0_rk) &
        -log_gamma(t-j1-J_2-J_3+1.0_rk) &
        -log_gamma(t-J_1-j2-J_3+1.0_rk) &
        -log_gamma(t-J_1-J_2-j3+1.0_rk) &
        -log_gamma(j1+j2+J_1+J_2-t+1.0_rk) &
        -log_gamma(j2+j3+J_2+J_3-t+1.0_rk) &
        -log_gamma(j3+j1+J_3+J_1-t+1.0_rk)
    end function log_

    function triangle_coefficient(a,b,c) result(tri)
      implicit none

      real(rk),intent(in) :: a,b,c
      real(rk) :: tri
      integer(ik) :: xa

      tri = 0.0_rk
      if (a<0.or.b<0.or.c<0) return

      do xa = nint(abs(a-b)*2.0_rk),nint((a+b)*2.0_rk),2
        if (nint(c*2.0_rk)==xa) then
          tri = log_gamma(a+b-c+1.0_rk) &
            +log_gamma(a+c-b+1.0_rk) &
            +log_gamma(b+c-a+1.0_rk) &
            -log_gamma(a+b+c+2.0_rk)
          tri = exp(tri)
        endif
      enddo
    end function triangle_coefficient

  end function Wigner6j


  function get_sign(exponent) result(sign_)

    implicit none

    real(rk),intent(in) :: exponent
    real(rk) :: sign_

    if (mod(nint(exponent),2)==0) then
      sign_ = 1.0_rk
    else
      sign_ = -1.0_rk
    endif

  end function get_sign

end module magnetic_dipole
