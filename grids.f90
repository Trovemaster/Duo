    subroutine gridred(nsub,alpha,re,mes,rmin,rmax,h,r,z,add,iverbose)
    use accuracy, only : rk, ik, ark, out, bohr
    use Lobatto,   only : LobattoAbsWeights
    
    implicit none
!____________________________________________________________
!
     integer(ik)   :: nsub      ! type of grid (uniformely spaced or not)
     integer(ik)   :: mes       ! number of the mesh points
     real(kind=rk) :: rmin      ! left innermost point 0<rmin<re
     real(kind=rk) :: rmax      ! right outmost  point rmin<re<rmax
     real(kind=rk) :: alpha     ! alpha > 0
     real(kind=rk) :: re        ! re > 0
     integer(ik)   :: iverbose  ! verbose level
     intent(in)    :: alpha,re,rmin,rmax,mes,nsub ! input variables
!
!        | assumed radial coordinate belongs |
!        | to open interval ]0,+infty[       |
!____________________________________________________________
!
     real(kind=rk) :: h      ! uniform step of the z-coordinate treated
     real(kind=rk) :: r(mes) ! nonuniform grid points of the r-coordinate
     real(kind=rk) :: z(mes),zmin,zmax,y,s,t,w(mes)
     real(kind=rk) :: add(mes)
!      intent(out) zmin,zmax,h,r !output variables
! -----------------------------------------------------------------
!                intrisic scalars
!------------------------------------------------------------------
      integer       ::     j
!       real(kind=rk) ::     zmin, zmax, y, t, s ! intrisic scalars
      real(kind=rk), parameter :: zero = 0._rk, one = 1._rk, two = 2._rk
!
! --------------------------------------------------------------------
!
      if(rmax <= rmin) stop 'ERROR: rmax <= rmin'
      if(rmin <= zero) stop 'ERROR: rmin <= 0'
      if(alpha <= zero)stop 'ERROR: alpha <= 0'
!----------------------------------------------------------------------
      select case (nsub)
!----------------------------------------------------------------------

       case (1)         !  |  z = exp(-exp(-alpha*(r-re)))   |
!
         zmin = exp(-exp(-alpha*(rmin-re)))
         zmax = exp(-exp(-alpha*(rmax-re)))
         h    = ( zmax - zmin )/real( mes - 1 ,rk)
!----------------------------------------------------------------------
         s = (alpha/two)**2
         do j=1,mes
            y      = zmin + h*real(j-1,rk)
            r(j)   = re - (log (-log(y)) )/alpha
            z(j)   = -alpha*y*log(y)
            add(j) = s*(one + log(y)**2)
         end do
!----------------------------------------------------------------------
      case (2)         !  |  z = 1-[1+exp(alpha*(r-re))]^{-1}   |
!----------------------------------------------------------------------
         zmin = one - one/(one + exp(alpha*(rmin-re)))
         zmax = one - one/(one + exp(alpha*(rmax-re)))
         h    = ( zmax - zmin )/real( mes - 1,rk)
!----------------------------------------------------------------------
         s = (alpha/two)**2
         do j=1,mes
            y    =  zmin + h*real(j-1,rk)
            t    =  one - y
            r(j) =  re  + log(y/t)/alpha
            z(j) =  alpha*t*y
            add(j) =  s
         end do
!----------------------------------------------------------------------
      case (3)         !  |  z = arctan(alpha*(r-re))  |
!----------------------------------------------------------------------
         zmin = atan(alpha*(rmin-re))
         zmax = atan(alpha*(rmax-re))
         h    = ( zmax - zmin )/real( mes - 1,rk)
!----------------------------------------------------------------------
         do j=1,mes
            y    = zmin + h*real(j-1,rk)
            r(j) = re  + tan(y)/alpha
            z(j) =  alpha/(one+(alpha*(r(j)-re))**2)
            add(j) =  z(j)**2
         end do
!----------------------------------------------------------------------
      case (4)         !  |  z = (y-1)/(y+1); y = (r/re)^alpha  |
!----------------------------------------------------------------------
           y = (rmin/re)**alpha
           zmin = (y - one)/(y + one)
           y = (rmax/re)**alpha
           zmax = (y - one)/(y + one)
           h = ( zmax - zmin )/real( mes - 1,rk)
!----------------------------------------------------------------------
           s = (alpha**2-one)/(two**2)
           do j=1,mes
               y    = zmin + h*real(j-1, rk)
               r(j) = re * (( one + y)/( one - y))**(one/alpha)
               t=(r(j)/re)**alpha
               z(j) =  two*alpha*t/r(j)/(t+one)**2
               add(j) =  s/(r(j)**2)
          end do
!----------------------------------------------------------------------
      case (5)         !  |  z = 1-dexp(-alpha*r^re); |
!----------------------------------------------------------------------
           zmin = 1._rk - exp(-alpha*(rmin**re))
           zmax = 1._rk - exp(-alpha*(rmax**re))
           h    = ( zmax - zmin )/real( mes - 1, rk)
!----------------------------------------------------------------------
           write(out,*) zmin,zmax
           do j=1, mes, 1
              y       = zmin + h*real(j-1,rk)
              t       = -log(1._rk-y)
              r(j)    = (t/alpha)**(1._rk/re)
              z(j)    = r(j)/t/re/(1._rk-y)
              add(j)  = (1._rk-(re**2)*(1._rk+t**2))/(2._rk*r(j))**2
           enddo
!----------------------------------------------------------------------
      case (6)         !  |   Lobatto quadrature grid points |
!----------------------------------------------------------------------
           CALL LobattoAbsWeights(r,w,mes,rmin,rmax)
           zmin = rmin 
           zmax = rmax 
           h    = ( zmax - zmin )/real( mes - 1, rk)
!----------------------------------------------------------------------
     case (0)      !   uniform grid points: z=r |
!----------------------------------------------------------------------
           h = ( rmax - rmin )/real( mes - 1 ,rk)
           !
           if (iverbose>=4) then 
             write(out,"(a,I15)") 'Using a uniformly spaced grid with npoints = ', mes
             write(out,"(a,3F25.14)") 'rmin, rmax, step (in ang)   = ', rmin, rmax, h
             write(out,"(a,3F25.14)") 'rmin, rmax, step (in bohrs) = ', rmin/bohr, rmax/bohr, h/bohr
             write(out,"(a)")
           endif
           !
           do j=1,mes
               r(j) = rmin + h*real(j-1, rk)
               z(j) =  one
               add(j) = zero
           end do
           !
      case default
           !
           write(out,"(a,i15)") "illegal type of grid nsub = ",nsub
           stop "illegal type of grid"
           !
      end select
!==============================================================
      return
      end subroutine gridred

!==============================================================

!*****************************************************************************************
!*****************************************************************************************
!*****************************************************************************************
subroutine spline(x,y,n,yp1,ypn,y2)
! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
! x1 < x2 < ... < xN, and given values yp1 and ypn for the 1st derivative of the
! interpolating function at points 1 and n, respectively, this routine returns an array y2(1:n)
! of length n which contains the second derivatives of the interpolating function at
! the tabulated points xi. If yp1 and/or ypn are equal to 1e30 or larger, the routine is
! signaled to set the corresponding boundary condition for a natural spline, with zero second
! derivative on that boundary.
! From Numerical Recipes in Fortran 77
use accuracy, only : rk, out !,safe_max
!
integer            :: n,i,k
real(kind=rk)      :: yp1,ypn,x(n),y(n),y2(n)
real(kind=rk)      :: p,qn,sig,un,u(n)

! if (yp1.gt.safe_max) then   ! Lorenzo Lodi, 13 February 2014; modified so that natural splines are always used
  y2(1)=0._rk
  u(1)=0._rk
! else
!   y2(1)=-0.5_rk
!   u(1)=(3._rk/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
! endif

do i=2,n-1
  sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
  p=sig*y2(i-1)+2._rk
  y2(i)=(sig-1._rk)/p
  !
  if (abs(x(i+1)-x(i))<sqrt(epsilon(1.0_rk))) then
    write(out, "(A,2i5,' x(i),x(i+1) = ',2g16.8,' y(i),y(i+1) = ',2g16.8)") &
              'Error in spline: identical grid points: i,i+1 = ', i,i+1, x(i),x(i+1),y(i),y(i+1)
    stop 'Error in spline: identical grid points'
  endif
  !
  u(i)=(6._rk*( ( y(i+1)-y(i) )/( x(i+1)-x(i) )-( y(i)-y(i-1) )/( x(i)-x(i-1)) )/( x(i+1)-x(i-1) )-sig*u(i-1) )/p
enddo
! if (ypn.gt.safe_max) then   ! Lorenzo Lodi, 13 February 2014; modified so that natural splines are always used
  qn=0._rk
  un=0._rk
! else
!   qn=0.5_rk
!   un=(3._rk/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
! endif

y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1._rk)
do k=n-1,1,-1
  y2(k)=y2(k)*y2(k+1)+u(k)
enddo
return
end  subroutine spline
!*****************************************************************************************
subroutine splint(xa,ya,y2a,n,x,y)
!Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
! xai's in order), and given the array y2a(1:n), which is the output from spline above,
! and given a value of x, this routine returns a cubic-spline interpolated value y.
! From Numerical Recipes in Fortran 77
use accuracy, only : rk

integer :: n,klo,khi,k
real(kind=rk) :: x,y,xa(n),y2a(n),ya(n)
real(kind=rk) :: a,b,h
klo=1
khi=n

1 if (khi-klo.gt.1) then
  k=(khi+klo)/2
  if(xa(k).gt.x)then
    khi=k
  else
    klo=k
  endif
goto 1
endif
h=xa(khi)-xa(klo)
if (h.eq.0._rk) then 
    write(6,"('bad xa input in splint',f18.8)") x
    !stop 'bad xa input in splint'
    y = 0 
    return
endif
a=(xa(khi)-x)/h
b=(x-xa(klo))/h
y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6._rk
return
end subroutine splint
!*****************************************************************************************
!*****************************************************************************************
!*****************************************************************************************
SUBROUTINE QUINAT(N, X, Y, B, C, D, E, F)

use accuracy, only : rk, ik, ark
! Downloaded on 8-11-2007 from http://www.netlib.org/toms/600
! Modified by Lorenzo Lodi 11 Nov 2013 to make it fortran 95
! Note: this subroutine will NOT give the same results
! as MATHEMATICA v.9 using Interpolation[tbl, Method -> "Spline", InterpolationOrder -> 5]
! (although the interpolation is very close to what obtainable with mathematica).
! The difference is very probably due to different assumption for
! the derivatives at the end points.

!ALGORITHM 600, COLLECTED ALGORITHMS FROM ACM.
!ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.9, NO. 2,
!JUN., 1983, P. 258-259.
!QUINAT COMPUTES THE COEFFICIENTS OF A QUINTIC NATURAL QUINTIC SPLI
!S(X) WITH KNOTS X(I) INTERPOLATING THERE TO GIVEN FUNCTION VALUES:
!          S(X(I)) = Y(I)  FOR I = 1,2, ..., N.
!IN EACH INTERVAL (X(I),X(I+1)) THE SPLINE FUNCTION S(XX) IS A
!POLYNOMIAL OF FIFTH DEGREE:
!S(XX) = ((((F(I)*P+E(I))*P+D(I))*P+C(I))*P+B(I))*P+Y(I)    (*)
!      = ((((-F(I)*Q+E(I+1))*Q-D(I+1))*Q+C(I+1))*Q-B(I+1))*Q+Y(I+1)
!WHERE  P = XX - X(I)  AND  Q = X(I+1) - XX.
!(NOTE THE FIRST SUBSCRIPT IN THE SECOND EXPRESSION.)
!THE DIFFERENT POLYNOMIALS ARE PIECED TOGETHER SO THAT S(X) AND
!ITS DERIVATIVES UP TO S"" ARE CONTINUOUS.
!
!   INPUT:
!
!N          NUMBER OF DATA POINTS, (AT LEAST THREE, I.E. N > 2)
!X(1:N)     THE STRICTLY INCREASING OR DECREASING SEQUENCE OF
!           KNOTS.  THE SPACING MUST BE SUCH THAT THE FIFTH POWER
!           OF X(I+1) - X(I) CAN BE FORMED WITHOUT OVERFLOW OR
!           UNDERFLOW OF EXPONENTS.
!Y(1:N)     THE PRESCRIBED FUNCTION VALUES AT THE KNOTS.
!
!   OUTPUT:
!
!B,C,D,E,F  THE COMPUTED SPLINE COEFFICIENTS AS IN (*).
!    (1:N)  SPECIFICALLY
!           B(I) = S'(X(I)), C(I) = S"(X(I))/2, D(I) = S"'(X(I))/6,
!           E(I) = S""(X(I))/24,  F(I) = S""'(X(I))/120.
!           F(N) IS NEITHER USED NOR ALTERED.  THE FIVE ARRAYS
!           B,C,D,E,F MUST ALWAYS BE DISTINCT.
!
!   OPTION:
!
!IT IS POSSIBLE TO SPECIFY VALUES FOR THE FIRST AND SECOND
!DERIVATIVES OF THE SPLINE FUNCTION AT ARBITRARILY MANY KNOTS.
!THIS IS DONE BY RELAXING THE REQUIREMENT THAT THE SEQUENCE OF
!KNOTS BE STRICTLY INCREASING OR DECREASING.  SPECIFICALLY:
!
!IF X(J) = X(J+1) THEN S(X(J)) = Y(J) AND S'(X(J)) = Y(J+1),
!IF X(J) = X(J+1) = X(J+2) THEN IN ADDITION S"(X(J)) = Y(J+2).
!
!NOTE THAT S""(X) IS DISCONTINUOUS AT A DOUBLE KNOT AND, IN
!ADDITION, S"'(X) IS DISCONTINUOUS AT A TRIPLE KNOT.  THE
!SUBROUTINE ASSIGNS Y(I) TO Y(I+1) IN THESE CASES AND ALSO TO
!Y(I+2) AT A TRIPLE KNOT.  THE REPRESENTATION (*) REMAINS
!VALID IN EACH OPEN INTERVAL (X(I),X(I+1)).  AT A DOUBLE KNOT,
!X(J) = X(J+1), THE OUTPUT COEFFICIENTS HAVE THE FOLLOWING VALUES:
!  Y(J) = S(X(J))          = Y(J+1)
!  B(J) = S'(X(J))         = B(J+1)
!  C(J) = S"(X(J))/2       = C(J+1)
!  D(J) = S"'(X(J))/6      = D(J+1)
!  E(J) = S""(X(J)-0)/24     E(J+1) = S""(X(J)+0)/24
!  F(J) = S""'(X(J)-0)/120   F(J+1) = S""'(X(J)+0)/120
!AT A TRIPLE KNOT, X(J) = X(J+1) = X(J+2), THE OUTPUT
!COEFFICIENTS HAVE THE FOLLOWING VALUES:
!  Y(J) = S(X(J))         = Y(J+1)    = Y(J+2)
!  B(J) = S'(X(J))        = B(J+1)    = B(J+2)
!  C(J) = S"(X(J))/2      = C(J+1)    = C(J+2)
!  D(J) = S"'((X(J)-0)/6    D(J+1) = 0  D(J+2) = S"'(X(J)+0)/6
!  E(J) = S""(X(J)-0)/24    E(J+1) = 0  E(J+2) = S""(X(J)+0)/24
!  F(J) = S""'(X(J)-0)/120  F(J+1) = 0  F(J+2) = S""'(X(J)+0)/120
!
INTEGER(ik) :: N
real(rk) :: X(N), Y(N), B(N), C(N), D(N), E(N), F(N)
INTEGER(ik) :: I, M
real(rk) :: B1, P, PQ, PQQR, PR, P2, P3, Q, QR, Q2, Q3, R, R2, S, T, U, V

IF (N <= 2) return

!COEFFICIENTS OF A POSITIVE DEFINITE, PENTADIAGONAL MATRIX,
!STORED IN D,E,F FROM 2 TO N-2.

M = N - 2
Q = X(2) - X(1)
R = X(3) - X(2)
Q2 = Q*Q
R2 = R*R
QR = Q + R
D(1) = 0
E(1) = 0
D(2) = 0
IF (Q /= 0.0_rk) D(2) = 6.0_rk*Q*Q2/(QR*QR)

IF (M < 2) GO TO 40
   DO 30 I=2,M
     P = Q
     Q = R
     R = X(I+2) - X(I+1)
     P2 = Q2
     Q2 = R2
     R2 = R*R
     PQ = QR
     QR = Q + R
!     IF (Q) 20, 10, 20
     if(Q /= 0.0_rk) goto 20
     if(Q == 0.0_rk) goto 10

10   D(I+1) = 0
     E(I) = 0
     F(I-1) = 0
     GO TO 30
20   Q3 = Q2*Q
     PR = P*R
     PQQR = PQ*QR
     D(I+1) = 6.0_rk*Q3/(QR*QR)
     D(I) = D(I) + (Q+Q)*(15.0_rk*PR*PR+(P+R)*Q*(20.0_rk*PR+7.0_ark*Q2)+ &
  &   Q2*(8.0_rk*(P2+R2)+21.0_rk*PR+Q2+Q2))/(PQQR*PQQR)
     D(I-1) = D(I-1) + 6._ark*Q3/(PQ*PQ)
     E(I) = Q2*(P*QR+3.0_rk*PQ*(QR+R+R))/(PQQR*QR)
     E(I-1) = E(I-1) + Q2*(R*PQ+3.0_rk*QR*(PQ+P+P))/(PQQR*PQ)
     F(I-1) = Q3/PQQR
30 CONTINUE

40 IF (R.NE.0.0_rk) D(M) = D(M) + 6.0_rk*R*R2/(QR*QR)

!FIRST AND SECOND ORDER DIVIDED DIFFERENCES OF THE GIVEN FUNCTION
!VALUES, STORED IN B FROM 2 TO N AND IN C FROM 3 TO N
!RESPECTIVELY. CARE IS TAKEN OF DOUBLE AND TRIPLE KNOTS.

      DO 60 I=2,N
        IF (X(I).NE.X(I-1)) GO TO 50
        B(I) = Y(I)
        Y(I) = Y(I-1)
        GO TO 60
   50   B(I) = (Y(I)-Y(I-1))/(X(I)-X(I-1))
   60 CONTINUE

      DO 80 I=3,N
        IF (X(I).NE.X(I-2)) GO TO 70
        C(I) = B(I)*0.5D0
        B(I) = B(I-1)
        GO TO 80
   70   C(I) = (B(I)-B(I-1))/(X(I)-X(I-2))
   80 CONTINUE

!SOLVE THE LINEAR SYSTEM WITH C(I+2) - C(I+1) AS RIGHT-HAND SIDE.

      IF (M.LT.2) GO TO 100
      P = 0
      C(1) = 0
      E(M) = 0
      F(1) = 0
      F(M-1) = 0
      F(M) = 0
      C(2) = C(4) - C(3)
      D(2) = 1.0_rk/D(2)

      IF (M.LT.3) GO TO 100
      DO 90 I=3,M
        Q = D(I-1)*E(I-1)
        D(I) = 1.0_rk/(D(I)-P*F(I-2)-Q*E(I-1))
        E(I) = E(I) - Q*F(I-1)
        C(I) = C(I+2) - C(I+1) - P*C(I-2) - Q*C(I-1)
        P = D(I-1)*F(I-1)
   90 CONTINUE

  100 I = N - 1
      C(N-1) = 0.0_rk
      C(N) = 0.0_rk
      IF (N.LT.4) GO TO 120
      DO 110 M=4,N

!        I = N-2, ..., 2
        I = I - 1
        C(I) = (C(I)-E(I)*C(I+1)-F(I)*C(I+2))*D(I)
  110 CONTINUE

!     INTEGRATE THE THIRD DERIVATIVE OF S(X).

  120 M = N - 1
      Q = X(2) - X(1)
      R = X(3) - X(2)
      B1 = B(2)
      Q3 = Q*Q*Q
      QR = Q + R
!       IF (QR) 140, 130, 140
      if(QR /=0.0_rk) goto 140
      if(QR ==0.0_rk) goto 130

  130 V = 0.
      T = 0.
      GO TO 150
  140 V = C(2)/QR
      T = V
  150 F(1) = 0.0_rk
      IF (Q.NE.0.0_rk) F(1) = V/Q
      DO 180 I=2,M
        P = Q
        Q = R
        R = 0.0_rk
        IF (I.NE.M) R = X(I+2) - X(I+1)
        P3 = Q3
        Q3 = Q*Q*Q
        PQ = QR
        QR = Q + R
        S = T
        T = 0.0_rk
        IF (QR.NE.0.0_rk) T = (C(I+1)-C(I))/QR
        U = V
        V = T - S
!         IF (PQ) 170, 160, 170
         if(PQ /=0.0_rk) goto 170
         if(PQ ==0.0_rk) goto 160

  160   C(I) = C(I-1)
        D(I) = 0.0_rk
        E(I) = 0.0_rk
        F(I) = 0.0_rk
        GO TO 180
  170   F(I) = F(I-1)
        IF (Q.NE.0.0_rk) F(I) = V/Q
        E(I) = 5.0_rk*S
        D(I) = 10.0_rk*(C(I)-Q*S)
        C(I) = D(I)*(P-Q) + (B(I+1)-B(I)+(U-E(I))*P3-(V+E(I))*Q3)/PQ
        B(I) = (P*(B(I+1)-V*Q3)+Q*(B(I)-U*P3))/PQ - P*Q*(D(I)+E(I)*(Q-P))
  180 CONTINUE

!     END POINTS X(1) AND X(N).

      P = X(2) - X(1)
      S = F(1)*P*P*P
      E(1) = 0.0_rk
      D(1) = 0.0_rk
      C(1) = C(2) - 10.0_rk*S
      B(1) = B1 - (C(1)+S)*P

      Q = X(N) - X(N-1)
      T = F(N-1)*Q*Q*Q
      E(N) = 0.0_rk
      D(N) = 0.0_rk
      C(N) = C(N-1) + 10.0_rk*T
      B(N) = B(N) + (C(N)-T)*Q
      RETURN

END SUBROUTINE QUINAT

subroutine splint_quint(xa,ya,n,x,y, b,c,d,e,f)
!Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
! xai's in order), and given the arrays B,C,D,E,F (1:n), which is the output from spline above,
! and given a value of x, this routine returns a quintic-spline interpolated value y.
! Adapted from Numerical Recipes in Fortran 77
use accuracy, only : rk

integer, intent(in) :: n
real(kind=rk) , intent(in) :: xa(n),ya(n)
real(kind=rk) , intent(in)  :: x
real(kind=rk) , intent(out) :: y
real(kind=rk) ::  b(n),c(n),d(n),e(n), f(n)
real(kind=rk) ::  p, q
integer :: klo,khi,k, i
real(kind=rk) :: h
klo=1
khi=n

1 if (khi-klo > 1) then
  k=(khi+klo)/2
  if(xa(k) > x)then
    khi=k
  else
    klo=k
  endif
  goto 1
endif

h=xa(khi)-xa(klo)
if (h == 0._rk) stop 'bad xa input in splint'

i = klo
p = x -xa(klo)
q = xa(khi) - x
y = ((((F(I)*P+E(I))*P+D(I))*P+C(I))*P+B(I))*P+ya(I)

return
end subroutine splint_quint




recursive subroutine polint_rk(xa, ya, x, y, dy)

  use accuracy 
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

     if (any(abs(den(1:n-m)) < small_)) then                   ! interpolating function
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
end subroutine polint_rk
