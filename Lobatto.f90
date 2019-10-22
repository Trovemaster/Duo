module Lobatto
  !
  use accuracy,     only : hik, ik, rk, ark
  use timer
  !
  implicit none
  !
  private
  !
  public LobattoAbsWeights,derLobattoMat
  !
  !
  contains

 subroutine AllLobattoKE(KE,npnt2)
  implicit none
  integer(ik),intent(in) :: npnt2
  real(rk),intent(out) :: KE(1:npnt2,1:npnt2)
  integer(ik)::       bra,ket,Ntot,Ntotp2
  real(rk) :: a0,rmin
  real(rk),allocatable,dimension(:) ::      rpt, w
  real(rk), allocatable, dimension(:,:) :: derLobMat
  !
  Ntot=npnt2-1
  !
  Ntotp2=Ntot+2;
  !
  a0=10d0
  rmin = 0.0d0 !(SY)
  allocate(rpt(0:Ntot+1),w(0:Ntot+1))
  !
  call LobattoAbsWeights(rpt,w,Ntotp2,rmin,a0)
  !
  allocate(derLobMat(0:Ntot+1,0:Ntot+1))
  call derLobattoMat(derLobMat,Ntot,rpt,w)
  !
  ! sy; please dealocate when not needed and use Arraysstart/stop subroutnies 
  !
  do bra=1,Ntot+1
    do ket=1,Ntot+1
      !
      call KeLobattoDVR(KE(bra,ket),bra,ket,Ntot,rpt,w,derLobMat)
      !
    end do
  end do
  !
  return
  !
 end subroutine AllLobattoKE
 !
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 !
 !
  subroutine OlLobattoDVR(result,bra,ket)
    implicit none
    real(rk), intent(out)::   result
    integer(ik), intent(in)::   bra, ket
    !
    !
     if (bra == ket) then
        result=1d0
    else
        result=0d0
     end if
     !
     return
     !
  end subroutine OlLobattoDVR


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine KeLobattoDVR(result,bra,ket,Ntot,rpt,w,derLobMat)
    implicit none
    integer(ik), intent(in)::   bra,ket,Ntot
    real(rk), intent(out)::   result,rpt(0:Ntot+1),w(0:Ntot+1)
    real(rk), intent(in)::    derLobMat(0:Ntot+1,0:Ntot+1)
    real(rk)::                pf, intresult
    integer(ik)::             k
    !
    rpt = 0
    w = 0
    pf=1d0!/sqrt(w(bra)*w(ket))
    intresult=0d0;
    do k=0,Ntot+1
      intresult=intresult+w(k)*derLobMat(bra,k)*derLobMat(ket,k)
    end do


result=intresult*pf;

return

end subroutine KeLobattoDVR


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine rminus2LobattoDVR(result,bra,ket,Ntot,rpt)
    implicit none
    real(rk), intent(out)::   result
    integer(ik),intent(in)::    bra,ket,Ntot
    real(rk), intent(in)::    rpt(0:Ntot+1)

if (bra == ket) then
    result=1d0/(rpt(bra)**2)
else
    result=0d0
end if

return

end subroutine rminus2LobattoDVR


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


subroutine derLobattoMat(result,Ntot,rpt,w)
        implicit none
        integer(ik), intent(in)     ::      Ntot
        real(rk), intent(out)     ::      result(0:Ntot+1,0:Ntot+1)
        real(rk), intent(in)      ::      rpt(0:Ntot+1),w(0:Ntot+1)
        integer(ik)                 ::      n,eta

        do n=0,Ntot+1
          do eta=n,Ntot+1
             call derLobatto(result(n,eta),n,eta,Ntot,rpt,w)
             if(n.ne.eta) result(eta,n) = (-1.0d0) * (w(eta)/w(n)) * result(n,eta) 
             result(n,eta) = result(n,eta) * 1.0d0/sqrt(w(n))
             if(n.ne.eta) result(eta,n) = result(eta,n) * 1.0d0/sqrt(w(eta))
          end do
        end do
        return
end subroutine derLobattoMat


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine derLobatto(result,n,eta,Ntot,rpt,w)

!!!!!!!
!    Calculates u'_n(r_eta)
!!!!!!
        implicit none
    real(rk), intent(out)::   result
        integer(ik), intent(in)::   n, Ntot,eta
    real(rk), intent(in)::    w(0:Ntot+1), rpt(0:Ntot+1)
    real(rk)::                factor(0:Ntot+1),signrecorder
        integer(ik)::               j,Nhalf

    Nhalf=(Ntot+2)/2
    if (n==eta) then
            if (n==0) then
                    result=-1d0/(2d0*w(n))
            else if (n==Ntot+1) then
                    result=1d0/(2d0*w(n))
            else
                    result=0d0;
            end if
    else
            result=1d0/((rpt(n)-rpt(eta)))
        signrecorder = 1.0
            do j=0,Ntot+1
                    if ((j .ne. n) .and. (j .ne. eta)) then
               factor(j) = (rpt(eta)-rpt(j))/(rpt(n)-rpt(j))
               signrecorder = signrecorder * (factor(j) / abs(factor(j)))
               factor(j) = abs(factor(j))
           else
                factor(j) = 1.0d0
                    end if
            end do  
        call QuickSort(factor,1,Ntot+2)
        call Riffle(factor,Ntot+2,Nhalf)
        do j = 0,Ntot+1
           result = result * factor(j)
        end do
        result = result * signrecorder
    end if
    return
end subroutine derLobatto



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine LobattoAbsWeights(rpt,w,Ntotp2,rmin,a0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Calculates weights and abscissa of Gauss-Lobatto quadrature using subroutine from Manolopoulos paper
!
! This has been checked for Ntot=8 against the equivalent Mathematica routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        integer(ik), intent(in)::   Ntotp2
        real(rk), intent(in)::    rmin,a0
        real(rk), intent(out)::   w(1:Ntotp2),rpt(1:Ntotp2)
        real(rk)  :: gridlength,shift,scalelobatto,weight,z,p1,p2,p3,pi
        integer(ik) :: L,k,j,i

gridlength=a0-rmin
L=(Ntotp2+1)/2
shift=0.5d0*(a0+rmin)
scalelobatto=0.5d0*gridlength
weight=gridlength/(Ntotp2*(Ntotp2-1))
pi=dacos(-1d0)
rpt(1)=rmin
w(1)=weight
do k=2,L
        z=cos(pi*(4d0*k-3d0)/(4d0*Ntotp2-2d0))
        do i=1,50
                p2=0d0;
                p1=1d0;
                do j=1,Ntotp2-1
                        p3=p2
                        p2=p1
                        p1=((2d0*j-1d0)*z*p2-(j-1d0)*p3)/j
                end do
                p2=(Ntotp2-1d0)*(p2-z*p1)/(1d0-z*z)
                p3=(2*z*p2-Ntotp2*(Ntotp2-1)*p1)/(1d0-z*z)
                z=z-p2/p3
        end do
        rpt(k)=shift-scalelobatto*z
        rpt(ntotp2+1-k)=shift+scalelobatto*z
        w(k)=weight/(p1*p1)
        w(ntotp2+1-k)=w(k)
end do
rpt(ntotp2)=a0
w(ntotp2)=weight

return
end subroutine LobattoAbsWeights

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

recursive subroutine QuickSort(a, first, last)
    implicit none
    real(rk) a(*), x, t
    integer(ik) first, last
    integer(ik) i,j

    x = a( (first + last) / 2)
    i = first
    j = last
    do
        do while (a(i) < x)
           i = i+1
        end do
        do while (x < a(j))
        j = j-1
        end do
        if (i >= j) exit
        t = a(i); a(i) = a(j); a(j) = t
        i = i+1
        j = j-1
    end do
    if (first < i-1) call quicksort(a,first,i-1)
    if (j+1 < last) call quicksort(a,j+1,last)
end subroutine quicksort

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine Riffle(a,NAbsc,Nhalf)
    implicit none
    integer(ik),intent(in) :: NAbsc,Nhalf
    real(rk) ::  a(1:NAbsc),b(1:Nhalf),c(1:Nhalf)
    integer(ik) i,Nhalfplusone

    Nhalfplusone = Nhalf+1
    if(mod(NAbsc,2).eq.0) then
        do i=1,Nhalf
           b(i) = a(i)
           c(i) = a(NAbsc-i+1)
        enddo    
        do i = 1,Nhalf
           a(2*i) = c(i)
           a(2*i-1) = b(i)
        enddo
    else
        do i=1,Nhalfplusone
           b(i) = a(i)
           c(i) = a(NAbsc-i+1)
        enddo    
        do i = 1,Nhalfplusone
           if(i.ne.Nhalfplusone) a(2*i) = c(i)
           a(2*i-1) = b(i)
        enddo
    endif
end subroutine Riffle

end module Lobatto