module lapack

!
!  Simplistic type-agnostic LAPACK interface
!
  use accuracy
  use timer
  implicit none
  private verbose

  interface lapack_gelss
    module procedure lapack_cgelss
    module procedure lapack_zgelss
    module procedure lapack_qzgelss
    module procedure lapack_sgelss
    module procedure lapack_dgelss
    module procedure lapack_qgelss
  end interface ! lapack_gelss

  interface lapack_stev
    module procedure lapack_sstev
    module procedure lapack_dstev
  end interface ! lapack_stev

  interface lapack_heev
    module procedure lapack_cheev
    module procedure lapack_zheev
    module procedure lapack_qzheev
  end interface ! lapack_heev

  interface lapack_syev
    module procedure lapack_dsyev
    module procedure lapack_ssyev
    module procedure lapack_qsyev
  end interface ! lapack_syev

  interface lapack_syevd
    module procedure lapack_dsyevd
    module procedure lapack_ssyevd
  end interface ! lapack_syev

  interface lapack_syevr
    module procedure lapack_qsyevr
    module procedure lapack_dsyevr
    module procedure lapack_ssyevr
  end interface ! lapack_syevr

  interface lapack_gesvd
    module procedure lapack_dgesvd
    module procedure lapack_sgesvd
  end interface ! lapack_gesvd

  interface lapack_syevx
    module procedure lapack_dsyevx
  end interface ! lapack_syevx

  interface lapack_ginverse
    module procedure lapack_ginverse_real
    module procedure lapack_ginverse_double
  end interface ! lapack_ginverse


  interface lapack_gelsv
    module procedure lapack_zgesv
  end interface ! lapack_gelss


  integer,parameter:: verbose = 3
  real(rk),parameter :: singtol = -1.0d-12 ! 100.0d0*spacing(1.0d0)
  integer, parameter :: dp=kind(1.d0)
  !
  contains

  subroutine lapack_cgelss(a,b)
    complex, intent(inout) :: a(:,:)
    complex, intent(inout) :: b(:,:)

    external cgelss
    real                   :: s    (   min(size(a,dim=1),size(a,dim=2)))
    complex                :: work (50*max(size(a,dim=1),size(a,dim=2)))
    real                   :: rwork( 5*min(size(a,dim=1),size(a,dim=2)))
    integer                :: rank, info
    integer                :: na1, na2, nb1, nb2

    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call cgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0*spacing(1.0), rank, work, size(work), rwork, info)

    if (info/=0) then
      write (out,"(' cgelss returned ',i8)") info
      stop 'lapack_cgelss - cgelss failed'
    end if
  end subroutine lapack_cgelss

  subroutine lapack_zgelss(a,b)
    complex(dp), intent(inout) :: a(:,:)
    complex(dp), intent(inout) :: b(:,:)

    external zgelss
   
    double precision       :: s    (   min(size(a,dim=1),size(a,dim=2)))
    complex(dp)            :: work (50*max(size(a,dim=1),size(a,dim=2)))
    double precision       :: rwork( 5*min(size(a,dim=1),size(a,dim=2)))
    integer                :: rank, info
    integer                :: na1, na2, nb1, nb2

    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call zgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s,singtol, rank, work, size(work), rwork, info)

    if (info/=0) then
      write (out,"(' cgelss returned ',i8)") info
      stop 'lapack_cgelss - cgelss failed'
    end if
  end subroutine lapack_zgelss


  subroutine lapack_qzgelss(a_ark,b_ark)

    complex(ark), intent(inout) :: a_ark(:,:)
    complex(ark), intent(inout) :: b_ark(:,:)

    external zgelss

    complex(dp) :: a( size(a_ark,dim=1),size(a_ark,dim=2))
    complex(dp) :: b( size(b_ark,dim=1),size(b_ark,dim=2))
    
    double precision       :: s    (   min(size(a,dim=1),size(a,dim=2)))
    complex(dp)            :: work (50*max(size(a,dim=1),size(a,dim=2)))
    double precision       :: rwork( 5*min(size(a,dim=1),size(a,dim=2)))
    integer                :: rank, info
    integer                :: na1, na2, nb1, nb2

    a = a_ark
    b = b_ark

    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call zgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s,singtol, rank, work, size(work), rwork, info)

    if (info/=0) then
      write (out,"(' zgelss returned ',i8)") info
      stop 'lapack_qzgelss - zgelss failed'
    end if
    !
    a_ark = a
    b_ark = b
    !
  end subroutine lapack_qzgelss




  subroutine lapack_sgelss(a,b)
    real, intent(inout) :: a(:,:)
    real, intent(inout) :: b(:,:)

    external sgelss
    real                :: s    (   min(size(a,dim=1),size(a,dim=2)))
    real                :: work (50*max(size(a,dim=1),size(a,dim=2),size(b,dim=2)))
    integer             :: rank, info
    integer             :: na1, na2, nb1, nb2

    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call sgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0*spacing(1.0), rank, work, size(work), info)

    if (info/=0) then
      write (out,"(' sgelss returned ',i8)") info
      stop 'lapack_sgelss - sgelss failed'
    end if
  end subroutine lapack_sgelss

  subroutine lapack_dgelss(a,b,ierror)
    double precision, intent(inout) :: a(:,:)
    double precision, intent(inout) :: b(:,:)
    integer,intent(inout),optional  :: ierror

    external dgelss
    double precision    :: s    (   min(size(a,dim=1),size(a,dim=2))),tol
    double precision    :: work (70*max(size(a,dim=1),size(a,dim=2),size(b,dim=2)))
    double precision,allocatable  :: work_t (:),h(:,:)
    integer(ik)          :: rank,info1,info2,info
    integer(ik)          :: na1, na2, nb1, nb2, iw


    if (verbose>=4) write(out,"('lapack_dgelss...')")
    !
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    !
    if (verbose>=5) write(out,'(A,4i12)') 'lapack_dgelss...na1,nb1,na2,nb2 = ', na1,nb1,na2,nb2
    !
    allocate(h(na1,na2),stat=info1)
    !
    if (verbose>=5) write(out,"('lapack_dgelss...allocated')")
    !
    if (info1/=0) then
      !
      write (out,'(A,i8,A,2i0)') " allocation (h) error ", info1, " size = ", na1,na2
      stop 'lapack_dgelss - allocation of h failed'
      !
    end if
    !
    if (verbose>=5) write(out,"('lapack_dgelss: a=>h ')")
    !
    h = a
    !
     if (verbose>=5) write(out,"('lapack_dgelss: a=>h ... done')")
    !
    tol = singtol
    !
    call dgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s,singtol, rank, work, -1, info)
    !
    iw = int(work(1))
    !
    if (verbose>=5) write(out,"(A,i0)") '  size(work) = ', iw
    !
    allocate(work_t(iw),stat=info2)
    !
    if (verbose>=5) write(out,"('  work_t allocated? ')")
    !
    if (info2/=0) then
      !
      write (out,"(' allocation (of work) error ',i8,' iwork = ',i0)") info2,iw
      stop 'lapack_dgelss - allocation of work_t failed'
      !
    end if
    !
    call dgelss(na1,na2,nb2,h(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s,singtol, rank, work_t, iw, info)
    !
    a = h
    !
    deallocate(work_t,h)
    !
    if (info/=0.or.info1/=0.or.info2/=0) then
      !
      write (out,"(' dgelss returned ',i8)") info
      write (out,"(' dgelss: na1 = ',i8,' nb1 =  ',i8,' na2 = ',i8,' nb2 =  ',i8,' iwork = ',i0)") na1,nb1,na2,nb2,iw
      !
      if (present(ierror).and.info/=0) then
        !
        ierror = info
        !
      else
        !
        stop 'lapack_dgelss - dgelss failed'
        !
      endif
      !
    end if
    !
    if (rank/=na2.and.present(ierror)) then
      !
      ierror = 1
      !
    endif
    !
    if (verbose>=4) write(out,"('...lapack_dgelss done!')")
    !
  end subroutine lapack_dgelss


  subroutine lapack_qgelss(a_ark,b_ark,ierror)
    real(ark), intent(inout) :: a_ark(:,:)
    real(ark), intent(inout) :: b_ark(:,:)
    integer,intent(inout),optional  :: ierror

    double precision :: a( size(a_ark,dim=1),size(a_ark,dim=2))
    double precision :: b( size(b_ark,dim=1),size(b_ark,dim=2))

    external dgelss
    !
    double precision    :: s    (   min(size(a_ark,dim=1),size(a_ark,dim=2))),tol
    double precision    :: work (70*max(size(a_ark,dim=1),size(a_ark,dim=2),size(b_ark,dim=2)))
    double precision,allocatable  :: work_t (:),h(:,:)
    integer(ik)          :: rank,info1,info2,info
    integer(ik)          :: na1, na2, nb1, nb2, iw


    if (verbose>=4) write(out,"('lapack_qgelss...')")
    !
    a = real(a_ark,kind=dp)
    b = real(b_ark,kind=dp)
    !
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    !
    if (verbose>=5) write(out,'(A,4i12)') 'lapack_qgelss...na1,nb1,na2,nb2 = ', na1,nb1,na2,nb2
    !
    allocate(h(na1,na2),stat=info1)
    !
    if (verbose>=5) write(out,"('lapack_qgelss...allocated')")
    !
    if (info1/=0) then
      !
      write (out,'(A,i8,A,2i0)') " allocation (h) error ", info1, " size = ", na1,na2
      stop 'lapack_qgelss - allocation of h failed'
      !
    end if
    !
    if (verbose>=5) write(out,"('lapack_qgelss: a=>h ')")
    !
    h = a
    !
     if (verbose>=5) write(out,"('lapack_qgelss: a=>h ... done')")
    !
    tol = singtol
    !
    call dgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s,singtol, rank, work, -1, info)
    !
    iw = int(work(1))
    !
    if (verbose>=5) write(out,"(A,i0)") '  size(work) = ', iw
    !
    allocate(work_t(iw),stat=info2)
    !
    if (verbose>=5) write(out,"('  work_t allocated? ')")
    !
    if (info2/=0) then
      !
      write (out,"(' allocation (of work) error ',i8,' iwork = ',i0)") info2,iw
      stop 'lapack_qgelss - allocation of work_t failed'
      !
    end if
    !
    call dgelss(na1,na2,nb2,h(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s,singtol, rank, work_t, iw, info)
    !
    a = h
    !
    deallocate(work_t,h)
    !
    if (info/=0.or.info1/=0.or.info2/=0) then
      !
      write (out,"(' qgelss returned ',i8)") info
      write (out,"(' qgelss: na1 = ',i8,' nb1 =  ',i8,' na2 = ',i8,' nb2 =  ',i8,' iwork = ',i0)") na1,nb1,na2,nb2,iw
      !
      if (present(ierror).and.info/=0) then
        !
        ierror = info
        !
      else
        !
        stop 'lapack_dgelss - qgelss failed'
        !
      endif
      !
    end if
    !
    if (rank/=na2.and.present(ierror)) then
      !
      ierror = 1
      !
    endif
    !

    a_ark = real(a,kind=ark)
    b_ark = real(b,kind=ark)

    if (verbose>=4) write(out,"('...lapack_qgelss done!')")
    !
  end subroutine lapack_qgelss



  subroutine lapack_sstev(d,e,z)
    real, intent(inout) :: d(:)   ! In:  Diagonal elements of the matrix
                                  ! Out: Eigenvalues, ascending order
    real, intent(inout) :: e(:)   ! In:  Sub-/super-diagonal elements of the matrix
                                  ! Out: Destroyed
    real, intent(out)   :: z(:,:) ! Out: Eigenvectors

    real    :: work(max(1,2*size(d)-2))
    integer :: info
    integer :: nz1, nz2

    nz1 = size(z,dim=1) ; nz2 = size(z,dim=2)
    call sstev('V',size(d),d,e,z(1:nz1,1:nz2),nz1,work,info)

    if (info/=0) then
      write (out,"(' sstev returned ',i8)") info
      stop 'lapack_sstev - sstev failed'
    end if
  end subroutine lapack_sstev

  subroutine lapack_dstev(d,e,z)
    double precision, intent(inout) :: d(:)   ! In:  Diagonal elements of the matrix
                                              ! Out: Eigenvalues, ascending order
    double precision, intent(inout) :: e(:)   ! In:  Sub-/super-diagonal elements of the matrix
                                              ! Out: Destroyed
    double precision, intent(out)   :: z(:,:) ! Out: Eigenvectors

    double precision :: work(max(1,2*size(d)-2))
    integer          :: info
    integer          :: nz1, nz2

    nz1 = size(z,dim=1) ; nz2 = size(z,dim=2)
    call dstev('V',size(d),d,e,z(1:nz1,1:nz2),nz1,work,info)

    if (info/=0) then
      write (out,"(' dstev returned ',i8)") info
      stop 'lapack_dstev - dstev failed'
    end if
  end subroutine lapack_dstev

  subroutine lapack_cheev(h,e)
    complex, intent(inout) :: h(:,:)  ! In:  Hermitian matrix to be diagonalized
                                      ! Out: Eigenvectors
    real, intent(out)   :: e(:)       ! Out: Eigenvalues

    complex :: work(50*size(h,dim=1))
    real    :: rwork(3*size(h,dim=1))
    integer :: info
    integer :: nh1, nh2

    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call cheev('V','U',nh1,h(1:nh1,1:nh2),nh1,e(:),work,size(work),rwork,info)
    if (info/=0) then
      write (out,"(' cheev returned ',i8)") info
      stop 'lapack_cheev - cheev failed'
    end if
  end subroutine lapack_cheev

  subroutine lapack_zheev(h,e)
    integer, parameter :: dp=kind(1.d0)
    complex(dp), intent(inout) :: h(:,:)  ! In:  Hermitian matrix to be diagonalized
                                             ! Out: Eigenvectors
    double precision, intent(out)   :: e(:)  ! Out: Eigenvalues

    complex(dp)      :: work(50*size(h,dim=1))
    double precision :: rwork(3*size(h,dim=1))
    integer          :: info
    integer          :: nh1, nh2

    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call zheev('V','U',nh1,h(1:nh1,1:nh2),nh1,e(:),work,size(work),rwork,info)
    if (info/=0) then
      write (out,"(' zheev returned ',i8)") info
      stop 'lapack_zheev - zheev failed'
    end if
  end subroutine lapack_zheev



  subroutine lapack_qzheev(h_ark,e_ark)
    !
    complex(ark), intent(inout) :: h_ark(:,:)  ! In:  Hermitian matrix to be diagonalized
                                               ! Out: Eigenvectors
    real(ark), intent(out)      :: e_ark(:)    ! Out: Eigenvalues

    complex(dp) :: h(size(h_ark,dim=1),size(h_ark,dim=2)) 
    double precision   :: e(size(e_ark,dim=1)) 


    complex(dp)      :: work(50*size(h,dim=1))
    double precision :: rwork(3*size(h,dim=1))
    integer          :: info
    integer          :: nh1, nh2
    !
    h = h_ark
    !e = e_ark
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call zheev('V','U',nh1,h(1:nh1,1:nh2),nh1,e(:),work,size(work),rwork,info)
    if (info/=0) then
      write (out,"(' qzheev returned ',i8)") info
      stop 'lapack_qheev - zheev failed'
    end if
    !
    h_ark = h
    e_ark = e
    !
  end subroutine lapack_qzheev



  subroutine lapack_dgesvd(h)
    double precision, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    character(len=1) :: jobu,jobvt

    double precision,allocatable    :: work(:),u(:,:),vt(:,:),s(:)
    integer          :: info
    integer          :: nh1, nh2
    integer          :: lwork
    !
    if (verbose>=2) call TimerStart('lapack_dgesvd')
    !
    jobu = 'O'
    !
    jobvt = 'N'
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    lwork = 50*nh1
    !
    allocate(work(lwork),u(1,nh1),vt(1,nh2),s(min(nh1,nh2)),stat=info)
    !
    call ArrayStart('lapack_gesvd-arrays-work',info,size(u),kind(u))
    call ArrayStart('lapack_gesvd-arrays-work',info,size(vt),kind(vt))
    call ArrayStart('lapack_gesvd-arrays-work',info,size(s),kind(s))
    !
    call dgesvd(jobu,jobvt,nh1,nh2,h(1:nh1,1:nh2),nh1,s,u,nh1,vt,nh2,work,-1,info)
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
    call dgesvd(jobu,jobvt,nh1,nh2,h(1:nh1,1:nh2),nh1,s,u,nh1,vt,nh2,work,lwork,info)
    !
    deallocate(work,u,vt,s)
    !
    call ArrayStop('lapack_gesvd-arrays-work')
    !
    if (info/=0) then
      write (out,"(' dgesvd returned ',i8)") info
      stop 'lapack_dgesvd - dgesvd failed'
    end if
    !
    if (verbose>=2) call TimerStop('lapack_dgesvd')
    !
  end subroutine lapack_dgesvd



  subroutine lapack_sgesvd(h)
    real, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    character(len=1) :: jobu,jobvt

    real,allocatable    :: work(:),u(:,:),vt(:,:),s(:)
    integer          :: info
    integer          :: nh1, nh2
    integer          :: lwork
    !
    if (verbose>=2) call TimerStart('lapack_sgesvd: diagonalization')
    !
    jobu = 'O'
    !
    jobvt = 'N'
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    lwork = 50*nh1
    !
    allocate(work(lwork),u(nh1,nh1),vt(nh2,nh2),s(min(nh1,nh2)),stat=info)
    !
    call ArrayStart('lapack_gesvd-arrays-work',info,size(u),kind(u))
    call ArrayStart('lapack_gesvd-arrays-work',info,size(vt),kind(vt))
    call ArrayStart('lapack_gesvd-arrays-work',info,size(s),kind(s))
    !
    call sgesvd(jobu,jobvt,nh1,nh2,h(1:nh1,1:nh2),nh1,s,u,nh1,vt,nh2,work,-1,info)
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
    call sgesvd(jobu,jobvt,nh1,nh2,h(1:nh1,1:nh2),nh1,s,u,nh1,vt,nh2,work,lwork,info)
    !
    deallocate(work,u,vt,s)
    !
    call ArrayStop('lapack_gesvd-arrays-work')
    !
    if (info/=0) then
      write (out,"(' dgesvd returned ',i8)") info
      stop 'lapack_dgesvd - sgesvd failed'
    end if
    !
    !
    if (verbose>=2) call TimerStop('lapack_sgesvd: diagonalization')
    !
  end subroutine lapack_sgesvd




  subroutine lapack_dsyev(h,e,jobz)
    double precision, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    double precision, intent(out)   :: e(:)    ! Out: Eigenvalues
    character(len=1),intent(in),optional :: jobz
    character(len=1) :: jobz_

    double precision,allocatable    :: work(:)
    integer(ik)      :: info
    integer(ik)      :: nh1, nh2
    integer(ik)      :: lwork
    !
    !if (verbose>=2) call TimerStart('lapack_dsyev: diagonalization')
    !
    jobz_ = 'V'
    !
    if (present(jobz)) then
       jobz_ = jobz
    endif

    lwork = 50*size(h,dim=1)
    !
    allocate(work(lwork))
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call dsyev(jobz_,'U',nh1,h(1:nh1,1:nh2),nh1,e,work,-1,info)
    !
    if (int(work(1))>size(work)) then
      !
      lwork = int(work(1))
      !
      deallocate(work)
      !
      allocate(work(lwork))
      !
      call dsyev(jobz_,'U',nh1,h(1:nh1,1:nh2),nh1,e,work,lwork,info)
      !
    else
      !
      call dsyev(jobz_,'U',nh1,h(1:nh1,1:nh2),nh1,e,work,size(work),info)
      !
    endif
    !
    deallocate(work)
    !
    if (info/=0) then
      write (out,"(' dsyev returned ',i8)") info
      stop 'lapack_dsyev - dsyev failed'
    end if
    !
    !if (verbose>=2) call TimerStop('lapack_dsyev: diagonalization')
    !
  end subroutine lapack_dsyev


  subroutine lapack_qsyev(h_ark,e_ark,jobz)
    !
    real(ark), intent(inout) :: h_ark(:,:)  ! In:  Hermitian matrix to be diagonalized
                                             ! Out: Eigenvectors
    real(ark), intent(out)   :: e_ark(:)  ! Out: Eigenvalues

    real(dp) :: h(size(h_ark,dim=1),size(h_ark,dim=2)) 
    real(dp) :: e(size(e_ark,dim=1)) 
    !
    character(len=1),intent(in),optional :: jobz
    character(len=1) :: jobz_

    double precision,allocatable    :: work(:)
    integer(ik)      :: info,alloc
    integer(ik)      :: nh1, nh2
    integer(ik)      :: lwork
    !
    real(ark),allocatable ::  d_ark(:,:)
    !
    !if (verbose>=2) call TimerStart('lapack_qsyev: diagonalization')
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    allocate (d_ark(nh1,nh1),stat=alloc)
    call ArrayStart('lapack_qsyev-d',alloc,size(d_ark),kind(d_ark))
    !
    call MLdiag_ulen_ark(nh1,h_ark,e_ark,d_ark)
    !
    h_ark = d_ark
    !
    deallocate(d_ark)
    call ArrayStop('lapack_qsyev-d')
    !
    return 
    !
    h = real(h_ark,dp)
    e = real(e_ark,dp)
    !
    jobz_ = 'V'
    !
    if (present(jobz)) then
       jobz_ = jobz
    endif

    lwork = 50*size(h,dim=1)
    !
    allocate(work(lwork))
    !
    call dsyev(jobz_,'U',nh1,h(1:nh1,1:nh2),nh1,e,work,-1,info)
    !
    if (int(work(1))>size(work)) then
      !
      lwork = int(work(1))
      !
      deallocate(work)
      !
      allocate(work(lwork))
      !
      call dsyev(jobz_,'U',nh1,h(1:nh1,1:nh2),nh1,e,work,lwork,info)
      !
    else
      !
      call dsyev(jobz_,'U',nh1,h(1:nh1,1:nh2),nh1,e,work,size(work),info)
      !
    endif
    !
    deallocate(work)
    !
    if (info/=0) then
      write (out,"(' dsyev returned ',i8)") info
      stop 'lapack_qsyev - dsyev failed'
    end if
    !
    h_ark = real(h,ark)
    e_ark = real(e,ark)
    !
    !if (verbose>=2) call TimerStop('lapack_dsyev: diagonalization')
    !
  end subroutine lapack_qsyev



  subroutine lapack_ssyev(h,e,jobz)
    real, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                   ! Out: Eigenvectors
    real, intent(out)   :: e(:)    ! Out: Eigenvalues
    character(len=1),intent(in),optional :: jobz
    character(len=1) :: jobz_

    real,allocatable    :: work(:)
    integer          :: info
    integer          :: nh1, nh2
    integer          :: lwork
    !
    if (verbose>=2) call TimerStart('lapack_ssyev: diagonalization')
    !
    jobz_ = 'V'
    !
    if (present(jobz)) then
       jobz_ = jobz
    endif
    !
    lwork = 50*size(h,dim=1)
    !
    allocate(work(lwork))
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call ssyev(jobz_,'U',nh1,h(1:nh1,1:nh2),nh1,e,work,-1,info)
    !
    if (int(work(1))>size(work)) then
      !
      lwork = int(work(1))
      !
      deallocate(work)
      !
      allocate(work(lwork))
      !
      call ssyev(jobz_,'U',nh1,h(1:nh1,1:nh2),nh1,e,work,lwork,info)
      !
    else
      !
      call ssyev(jobz_,'U',nh1,h(1:nh1,1:nh2),nh1,e,work,size(work),info)
      !
    endif
    !
    deallocate(work)
    !
    if (info/=0) then
      write (out,"(' ssyev returned ',i8)") info
      stop 'lapack_ssyev - dsyev failed'
    end if
    !
    if (verbose>=2) call TimerStop('lapack_ssyev: diagonalization')
    !
  end subroutine lapack_ssyev



  subroutine lapack_dsyevd(h,e)
    double precision, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    double precision, intent(out)   :: e(:)    ! Out: Eigenvalues

    integer          :: info
    integer          :: nh1, nh2,liwork,lwork
    double precision,allocatable :: work(:)
    integer,allocatable :: iwork(:)


    lwork  = 50*size(h,dim=1)
    liwork = 50*size(h,dim=1)
    !
    allocate(work(lwork),iwork(liwork))
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call dsyevd('V','U',nh1,h(1:nh1,1:nh2),nh1,e,work,-1,iwork,-1,info)
    !
    if (int(work(1))>size(work).or.int(iwork(1))>size(iwork)) then
      !
      lwork = int(work(1))
      liwork = int(iwork(1))
      !
      deallocate(work,iwork)
      !
      allocate(work(lwork),iwork(liwork))
      !
      call dsyevd('V','U',nh1,h(1:nh1,1:nh2),nh1,e,work,lwork,iwork,liwork,info)
      !
    else
      !
      call dsyevd('V','U',nh1,h(1:nh1,1:nh2),nh1,e,work,lwork,iwork,liwork,info)
      !
    endif
    !
    deallocate(work,iwork)
    !
    if (info/=0) then
      write (out,"(' dsyevd returned ',i8)") info
      stop 'lapack_dsyevd - dsyev failed'
    end if
    !
  end subroutine lapack_dsyevd



  subroutine lapack_ssyevd(h,e)
    real, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                   ! Out: Eigenvectors
    real, intent(out)   :: e(:)    ! Out: Eigenvalues

    real     :: work(10*size(h,dim=1))
    integer  :: iwork(10*size(h,dim=1))
    integer  :: info
    integer  :: nh1, nh2,liwork,lwork
    real,allocatable :: twork(:)
    integer,allocatable :: tiwork(:)

    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    lwork = size(work)
    iwork = size(iwork)
    !
    call ssyevd('V','U',nh1,h(1:nh1,1:nh2),nh1,e,work,-1,iwork,-1,info)
    !
    lwork = int(work(1))
    liwork = int(iwork(1))
    !
    allocate(twork(lwork),tiwork(liwork))
    !
    call ssyevd('V','U',nh1,h(1:nh1,1:nh2),nh1,e,work,lwork,iwork,liwork,info)
    !
    if (info/=0) then
      write (out,"(' ssyev returned ',i8)") info
      stop 'lapack_ssyevd - ssyevd failed'
    end if
    !
    deallocate(twork,tiwork)
    !
  end subroutine lapack_ssyevd


  subroutine lapack_ssyevr(h,e,rng,jobz,iroots,vrange,irange)
    real, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    real, intent(out)   :: e(:)    ! Out: Eigenvalues
    character(len=1),intent(in),optional   :: rng,jobz
    real, intent(in),optional  :: vrange(2)  ! In:  energy range to be computed
    integer         , intent(in),optional  :: irange(2)  ! In:  index to be computed
    integer         , intent(out),optional :: iroots     ! Out: Number of roots found
    !
    real :: vl,vu,abstol
    integer          :: info
    integer          :: nh1, nh2, il, ir,n,k, alloc
    character(len=1) :: rng_,jobz_
    integer          :: ib,jb,niwork,nwork,ldz,msize
    real,allocatable :: a(:,:),work(:),twork(:)
    integer,allocatable :: iwork(:),isuppz(:),tiwork(:)
    !
    if (verbose>=2) call TimerStart('lapack_ssyevr: diagonalization')
    !
    rng_ = 'A'
    jobz_ = 'V'
    il = 1   ; ir = 100000
    vl = -.001 ; vu = 100000.0
    !
    if (present(rng)) then
       rng_ = rng
    endif
    !
    if (present(jobz)) then
       jobz_ = jobz
    endif
    !
    if (present(vrange).and.rng_=='V') then
       vl = vrange(1) ; vu = vrange(2)
    endif
    !
    if (present(irange).and.rng_=='I') then
       il = irange(1) ; ir = irange(2)
       if (irange(2)>=size(h,dim=1)) rng_ = 'A'
    endif
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    msize = nh2
    if (rng_=='I') msize = ir-il+1
    !
    if (nh1/=size(e)) then
      write (out,"(' ssyevr inconsistent sizes of h and e :',2i8)") nh1/=size(e)
      stop 'lapack_ssyevr inconsistent sizes of h and e'
    end if
    !
    abstol =(small_)
    !
    nwork = 50*nh1
    niwork = 20*nh1
    ldz = 2*nh1
    !
    if (verbose>=4) then
      write(out,"(//'MATRIX:')")
      do ib=1,nh1
        do jb=ib,nh1
          write(out,"(2i8,f19.8)") ib,jb,h(jb,ib)
        enddo
      enddo
    endif
    !
    allocate(work(nwork),stat=alloc)
    call ArrayStart('lapack_ssyevr-arrays-work',alloc,size(work),kind(work))
    allocate(iwork(niwork),stat=alloc)
    call ArrayStart('lapack_ssyevr-arrays',alloc,size(iwork),kind(iwork))
    allocate(isuppz(2*msize),stat=alloc)
    call ArrayStart('lapack_ssyevr-arrays',alloc,size(isuppz),kind(isuppz))
    !
    !if (verbose>=3) call MemoryReport
    !if (verbose>=4) write(out,"(/'range = ',a)") rng_
    !
    if (verbose>=4) write(out,"(/'abstol = ',e18.8)") abstol
    call ssyevr(jobz_,rng_,'U',nh1,h(1:nh1,1:nh2),nh1,vl,vu,il,ir,abstol,n,e(1:nh1),h(1:nh1,1:nh2),&
                                                          nh1,isuppz,work,-1,iwork,-1,info)
    !
    if (verbose>=4) write(out,"(/'nwork = ',i8,' niwork = ',i8,' changed to ',2i8)") nwork,niwork,int(work(1)),int(iwork(1))
    !
    nwork = int(work(1))
    niwork = int(iwork(1))
    !
    !if (verbose>=3) write(out,"(/'nwork = ',i8,' niwork = ',i8)") nwork,niwork
    !
    deallocate(work,iwork)
    call ArrayStop('lapack_ssyevr-arrays')
    call ArrayStop('lapack_ssyevr-arrays-work')
    !
    allocate(twork(nwork),stat=alloc)
    call ArrayStart('lapack_ssyevr-arrays-twork',alloc,size(twork),kind(twork))
    allocate(tiwork(niwork),stat=alloc)
    call ArrayStart('lapack_ssyevr-arrays',alloc,size(tiwork),kind(tiwork))
    !
    allocate (a(nh1,msize),stat=alloc)
    call ArrayStart('lapack_ssyevr-a',alloc,size(a),kind(a))
    !
    if (verbose>=3) call MemoryReport
    !
    call ssyevr(jobz_,rng_,'U',nh1,h(1:nh1,1:nh2),nh1,vl,vu,il,ir,abstol,n,e(1:nh1),a(1:nh1,1:msize),&
                nh1,isuppz,twork(1:nwork),nwork,tiwork(1:niwork),niwork,info)
    !
    if (verbose>=3) call TimerStop('lapack_ssyevr: diagonalization finished')
    !
    !$omp parallel do private(k) shared(h) schedule(dynamic)
    do k=1,msize
       h(:,k) = a(:,k)
    enddo
    !$omp end parallel do

    !h(1:nh1,1:nh2) = a(1:nh1,1:nh2)
    !
    deallocate(a,twork,tiwork,isuppz)
    call ArrayStop('lapack_ssyevr-a')
    call ArrayStop('lapack_ssyevr-arrays')
    call ArrayStop('lapack_ssyevr-arrays-twork')
    !
    if (verbose>=5) then
      write(out,"(//'ssyevr solution:')")
      do ib=1,n
         write(out,"(i8,f19.8)") ib,e(ib)
         do jb=1,nh1
           write(out,"(2i8,f19.8)") ib,jb,h(jb,ib)
         enddo
      enddo
      !
      call MemoryReport
      !
    endif
    if (present(iroots))  iroots = n
    !
    if (info/=0) then
      write (out,"(' ssyevr returned ',i8)") info
      stop 'lapack_ssyevr - ssyevr failed'
    end if
    !
    if (verbose>=2) call TimerStop('lapack_ssyevr: diagonalization')
    !
  end subroutine lapack_ssyevr


  subroutine lapack_dsyevr(h,e,rng,jobz,iroots,vrange,irange)
    double precision, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    double precision, intent(out)   :: e(:)    ! Out: Eigenvalues
    character(len=1),intent(in),optional   :: rng
    character(len=1),intent(in),optional   :: jobz
    double precision, intent(in),optional  :: vrange(2)  ! In:  energy range to be computed
    integer         , intent(in),optional  :: irange(2)  ! In:  index to be computed
    integer         , intent(out),optional :: iroots     ! Out: Number of roots found

    double precision :: vl,vu,abstol
    integer          :: info
    integer          :: nh1, nh2, il, ir,n=1,k, alloc
    character(len=1) :: rng_,jobz_
    integer          :: ib,jb,niwork,nwork,ldz,msize
    integer(hik)     :: hsize
    double precision,allocatable :: a(:,:),work(:),twork(:)
    integer,allocatable :: iwork(:),isuppz(:),tiwork(:)
    !
    if (verbose>=4) call TimerStart('lapack_dsyevr: diagonalization')
    !
    rng_ = 'A'
    jobz_ = 'V'
    il = 1   ; ir = 100000
    vl = -.001 ; vu = 100000.0
    !
    if (present(rng)) then
       rng_ = rng
    endif
    !
    if (present(jobz)) then
       jobz_ = jobz_
    endif
    !
    if (present(vrange).and.rng_=='V') then
       vl = vrange(1) ; vu = vrange(2)
    endif
    !
    if (present(irange).and.rng_=='I') then
       il = irange(1) ; ir = irange(2)
       if (irange(2)>=size(h,dim=1)) rng_ = 'A'
    endif
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    if (nh1<1) then
       iroots = 0
       return
    endif
    !
    msize = nh2
    if (rng_=='I') msize = ir-il+1
    !
    if (nh1/=size(e)) then
      write (out,"(' dsyevr inconsistent sizes of h and e :',2i8)") nh1/=size(e)
      stop 'lapack_dsyevr inconsistent sizes of h and e'
    end if
    !
    abstol =(small_)
    !
    nwork = 50*nh1
    niwork = 20*nh1
    ldz = 2*nh1
    !
    if (verbose>=5) then
      write(out,"(//'MATRIX:')")
      do ib=1,nh1
        do jb=ib,nh1
          write(out,"(2i8,f19.8)") ib,jb,h(jb,ib)
          !h(jb,ib) = h(ib,jb)
        enddo
      enddo
    endif
    !
    if (verbose>=4) write(out,"(/'nwork = ',i8,' niwork = ',i8)") nwork,niwork
    !
    allocate(work(nwork),stat=alloc)
    call ArrayStart('lapack_dsyevr-arrays-work',alloc,size(work),kind(work))
    allocate(iwork(niwork),stat=alloc)
    call ArrayStart('lapack_dsyevr-arrays',alloc,size(iwork),kind(iwork))
    allocate(isuppz(2*msize),stat=alloc)
    call ArrayStart('lapack_dsyevr-arrays',alloc,size(isuppz),kind(isuppz))
    !
    if (verbose>=4) call MemoryReport
    !if (verbose>=4) write(out,"(/'range = ',a)") rng_
    !
    hsize = int(nh1,hik)*int(msize,hik)
    !
!    allocate (a(nh1,msize),stat=alloc)
    allocate (a(nh1,nh2),stat=alloc) ! Lorenzo Lodi; should be right?
    call ArrayStart('lapack_dsyevr-a',alloc,1,kind(a),hsize)
    !
    if (verbose>=4) write(out,"(/'abstol = ',e18.8)") abstol
    call dsyevr(jobz_,rng_,'U',nh1,h(1:nh1,1:nh2),nh1,vl,vu,il,ir,abstol,n,e(1:nh1),a(1:nh1,1:nh2),&
                                                          nh1,isuppz,work,-1,iwork,-1,info)
    !
    if (verbose>=4) write(out,"(/'nwork = ',i8,' niwork = ',i8,' changed to ',2i8)") nwork,niwork,int(work(1)),int(iwork(1))
    !
    nwork = int(work(1))
    niwork = int(iwork(1))
    !
    deallocate(work,iwork)
    call ArrayStop('lapack_dsyevr-arrays')
    call ArrayStop('lapack_dsyevr-arrays-work')
    !
    allocate(twork(nwork),stat=alloc)
    call ArrayStart('lapack_dsyevr-arrays-twork',alloc,size(twork),kind(twork))
    allocate(tiwork(niwork),stat=alloc)
    call ArrayStart('lapack_dsyevr-arrays',alloc,size(tiwork),kind(tiwork))
    !
    if (verbose>=4) call MemoryReport
    !
    call dsyevr(jobz_,rng_,'U',nh1,h,nh1,vl,vu,il,ir,abstol,n,e,a,nh1,isuppz,twork,nwork,tiwork,niwork,info)
    !
    !$omp parallel do private(k) shared(h) schedule(dynamic)
    do k=1,msize
       h(:,k) = a(:,k)
    enddo
    !$omp end parallel do

    !h(1:nh1,1:nh2) = a(1:nh1,1:nh2)
    !
    deallocate(a,twork,tiwork,isuppz)
    call ArrayStop('lapack_dsyevr-a')
    call ArrayStop('lapack_dsyevr-arrays')
    call ArrayStop('lapack_dsyevr-arrays-twork')
    !
    if (verbose>=5) then
      write(out,"(//'dsyevr solution:')")
      do ib=1,n
         write(out,"(i8,f19.8)") ib,e(ib)
         do jb=1,nh1
           write(out,"(2i8,f19.8)") ib,jb,h(jb,ib)
         enddo
      enddo
      !
      if (verbose>=4) call MemoryReport
      !
    endif
    if (present(iroots))  iroots = n
    !
    if (info/=0) then
      write (out,"(' dsyevr returned ',i8)") info
      stop 'lapack_dsyevr - dsyevr failed'
    end if
    !
    if (verbose>=4) call TimerStop('lapack_dsyevr: diagonalization')
    !
  end subroutine lapack_dsyevr



  subroutine lapack_qsyevr(h_ark,e_ark,rng,jobz,iroots,vrange,irange)
    !
    real(ark), intent(inout) :: h_ark(:,:)  ! In:  Hermitian matrix to be diagonalized
                                             ! Out: Eigenvectors
    real(ark), intent(out)   :: e_ark(:)  ! Out: Eigenvalues

    double precision,allocatable   :: h(:,:) 
    double precision,allocatable   :: e(:) 

    character(len=1),intent(in),optional   :: rng
    character(len=1),intent(in),optional   :: jobz
    double precision, intent(in),optional  :: vrange(2)  ! In:  energy range to be computed
    integer         , intent(in),optional  :: irange(2)  ! In:  index to be computed
    integer         , intent(out),optional :: iroots     ! Out: Number of roots found
    integer(ik) :: nh1,nh2,ne,alloc
    !
    real(ark),allocatable ::  d_ark(:,:)
    real(ark)   :: energy_
    integer(ik) :: ilevel,jlevel
    !
    nh1 = size(h_ark,dim=1) ; nh2 = size(h_ark,dim=2)
    ne = size(e_ark,dim=1)
    !
    allocate (d_ark(nh1,nh1),stat=alloc)
    call ArrayStart('lapack_qsyevr-d',alloc,size(d_ark),kind(d_ark))
    !
    call MLdiag_ulen_ark(nh1,h_ark,e_ark,d_ark)
    !
    iroots = irange(2)
    !
    do ilevel = 1,iroots
      !
      energy_ = e_ark(ilevel)
      !
      do jlevel=ilevel+1,nh1
        !
        if ( energy_>e_ark(jlevel) ) then
          !
          ! energy
          !
          energy_=e_ark(jlevel)
          e_ark(jlevel) = e_ark(ilevel)
          e_ark(ilevel) = energy_
          !
          ! basis function
          !
          h_ark(:,1) = d_ark(:,jlevel)
          d_ark(:,jlevel) = d_ark(:,ilevel)
          d_ark(:,ilevel) = h_ark(:,1)
          !
        endif
        !
      enddo
      !
    enddo    
    !
    h_ark(:,1:iroots) = d_ark(:,1:iroots)
    !
    deallocate(d_ark)
    call ArrayStop('lapack_qsyevr-d')
    !
    return 
    !
    allocate (h(nh1,nh2),e(ne),stat=alloc)
    call ArrayStart('lapack_qsyevr-alloc',alloc,size(h),kind(h))
    call ArrayStart('lapack_qsyevr-alloc',alloc,size(e),kind(e))

    !
    h = real(h_ark,dp)
    e = real(e_ark,dp)
    !
    if (present(rng)) then 
      !
      if (present(jobz)) then 
        !
        if (present(iroots)) then 
          !
          if (present(vrange)) then 
            !
            if (present(vrange)) then 
              !
              if (present(vrange)) then 
                !
                call lapack_dsyevr(h,e,rng,jobz,iroots,vrange,irange)
                !
              else
                !
                call lapack_dsyevr(h,e,rng,jobz,iroots,vrange)
                !
              endif
              !
            else
              !
              call lapack_dsyevr(h,e,rng,jobz,iroots)
              !
            endif
            !
          else
            !
            call lapack_dsyevr(h,e,rng,jobz)
            !
          endif
          !
        else
          !
          call lapack_dsyevr(h,e,rng)
          !
        endif
        !
      else
        !
        call lapack_dsyevr(h,e)
        !
      endif
      !
    else
      !
      stop 'lapack_qsyevr not enoup parameters'
      !
    endif
    !
    h_ark = real(h,ark)
    e_ark = real(e,ark)
    !
    deallocate(h,e)
    !
    call ArrayStop('lapack_qsyevr-alloc')
    !
  end subroutine lapack_qsyevr



  subroutine lapack_dsyevx(h,e,iroots,vrange,irange,tol)
    double precision, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    double precision, intent(out)   :: e(:)    ! Out: Eigenvalues
    double precision, intent(in),optional  :: vrange(2)  ! In:  energy range to be computed
    integer         , intent(in),optional  :: irange(2)  ! In:  index to be computed
    integer         , intent(out),optional :: iroots     ! Out: Number of roots found
    double precision, intent(in),optional  :: tol    ! In: abs. tolerance

    double precision :: vl,vu,abstol
    integer          :: info
    integer          :: nh1, nh2, il, ir,n,k, alloc
    character(len=1) :: rng,jobz
    integer          :: ib,jb,niwork,nwork,msize
    double precision,allocatable :: a(:,:),work(:),twork(:)
    integer,allocatable :: iwork(:),ifail(:)
    !
    if (verbose>=2) call TimerStart('lapack_dsyevx: diagonalization')
    !
    rng = 'A'
    jobz = 'V'
    il = 1   ; ir = 100000
    vl = -.001 ; vu = 100000.0
    !
    if (present(vrange)) then
       rng = 'V'
       vl = vrange(1) ; vu = vrange(2)
    endif
    !
    abstol = sqrt(small_)
    if (present(tol))  abstol = tol
    !
    if (present(irange)) then
       rng = 'I'
       il = irange(1) ; ir = irange(2)
       if (irange(2)>=size(h,dim=1)) rng = 'A'
    endif
    !
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    !
    msize = nh2
    if (rng=='I') msize = ir-il+1
    !
    if (nh1/=size(e)) then
      write (out,"(' dsyevr inconsistent sizes of h and e :',2i8)") nh1/=size(e)
      stop 'lapack_dsyevx inconsistent sizes of h and e'
    end if
    !
    nwork = 50*nh1
    niwork = 5*nh1
    !
    if (verbose>=4) then
      write(out,"(//'MATRIX:')")
      do ib=1,nh1
        do jb=ib,nh1
          write(out,"(2i8,f19.8)") ib,jb,h(jb,ib)
        enddo
      enddo
    endif
    !
    if (verbose>=4) write(out,"(/'nwork = ',i8,' niwork = ',i8)") nwork,niwork
    !
    allocate(work(nwork),stat=alloc)
    call ArrayStart('lapack_dsyevx-arrays-work',alloc,size(work),kind(work))
    allocate(iwork(niwork),stat=alloc)
    call ArrayStart('lapack_dsyevx-arrays',alloc,size(iwork),kind(iwork))
    allocate(ifail(nh1),stat=alloc)
    call ArrayStart('lapack_dsyevx-arrays',alloc,size(ifail),kind(ifail))
    !
    if (verbose>=3) call MemoryReport
    !if (verbose>=3) write(out,"(/'range = ',a)") rng
    !
    if (verbose>=4) write(out,"(/'abstol = ',e18.8)") abstol
    call dsyevx('V',rng,'U',nh1,h(1:nh1,1:nh2),nh1,vl,vu,il,ir,abstol,n,e(1:nh1),h(1:nh1,1:nh2),&
                                                          nh1,work,-1,iwork,ifail,info)
    !
    nwork = int(work(1))
    !
    if (verbose>=4) write(out,"(/'nwork = ',i8)") nwork
    !
    deallocate(work,iwork)
    call ArrayStop('lapack_dsyevx-arrays')
    call ArrayStop('lapack_dsyevx-arrays-work')
    !
    allocate(twork(nwork),stat=alloc)
    call ArrayStart('lapack_dsyevx-arrays-twork',alloc,size(twork),kind(twork))
    !
    allocate (a(nh1,msize),stat=alloc)
    call ArrayStart('lapack_dsyevx-a',alloc,size(a),kind(a))
    !
    if (verbose>=3) call MemoryReport
    !
    call dsyevx('V',rng,'U',nh1,h(1:nh1,1:nh2),nh1,vl,vu,il,ir,abstol,n,e(1:nh1),a(1:nh1,1:msize),&
                nh1,twork(1:nwork),nwork,iwork(1:nh1),ifail,info)
    !
    !$omp parallel do private(k) shared(h) schedule(dynamic)
    do k=1,msize
       h(:,k) = a(:,k)
    enddo
    !$omp end parallel do

    !h(1:nh1,1:nh2) = a(1:nh1,1:nh2)
    !
    deallocate(a,twork,ifail)
    call ArrayStop('lapack_dsyevx-a')
    call ArrayStop('lapack_dsyevx-arrays')
    call ArrayStop('lapack_dsyevx-arrays-twork')
    !
    if (verbose>=5) then
      write(out,"(//'dsyevr solution:')")
      do ib=1,n
         write(out,"(i8,f19.8)") ib,e(ib)
         do jb=1,nh1
           write(out,"(2i8,f19.8)") ib,jb,h(jb,ib)
         enddo
      enddo
      !
      call MemoryReport
      !
    endif
    if (present(iroots))  iroots = n
    !
    if (info/=0) then
      write (out,"(' dsyevr returned ',i8)") info
      stop 'lapack_dsyevx - dsyevr failed'
    end if
    !
    if (verbose>=2) call TimerStop('lapack_dsyevx: diagonalization')
    !
  end subroutine lapack_dsyevx



  subroutine lapack_ginverse_real(amat)
    real, intent(inout) :: amat(:,:) ! In: matrix to invert
                                         ! Out: generalized inverse of the matrix

    real :: eval(size(amat,dim=1))
    real :: evec(size(amat,dim=1),size(amat,dim=1))
    real :: eps

    evec = amat
    call lapack_syev(evec(:,:),eval(:))
    eps = 100.0*spacing(maxval(eval))
    where (abs(eval)>eps)
      eval = 1.0 / sqrt(eval)
    elsewhere
      eval = 0.0
    end where
    evec = evec * spread(eval,dim=1,ncopies=size(evec,dim=1))
    amat = matmul(evec,transpose(evec))

  end subroutine lapack_ginverse_real

  subroutine lapack_ginverse_double(amat)
    double precision, intent(inout) :: amat(:,:) ! In: matrix to invert
                                         ! Out: generalized inverse of the matrix

    double precision :: eval(size(amat,dim=1))
    double precision :: evec(size(amat,dim=1),size(amat,dim=1))
    double precision :: eps

    evec = amat
    call lapack_syev(evec(:,:),eval(:))
    eps = 100.0d0*spacing(maxval(eval))
    !
    if (any(eval(:)<0)) then
      !
      write(out,"('lapack_ginverse_double: negative evals')")
      stop 'lapack_ginverse_double: negative evals'
      !
    endif
    !
    where (abs(eval)>eps)
      eval = 1.0d0 / sqrt(eval)
    elsewhere
      eval = 0.0d0
    end where
    evec = evec * spread(eval,dim=1,ncopies=size(evec,dim=1))
    amat = matmul(evec,transpose(evec))

  end subroutine lapack_ginverse_double



  subroutine lapack_zgesv(a,b)
    integer, parameter :: dp=kind(1.d0)
    complex(dp), intent(inout) :: a(:,:)
    complex(dp), intent(inout) :: b(:,:)
    integer(ik) :: ipiv(max(size(a,dim=1),1))

    external zgesv
    integer                :: info
    integer                :: na1, na2, nb1, nb2

    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call zgesv(na1,nb2,a(1:na1,1:na2),na1,ipiv,b(1:nb1,1:nb2),nb1,info)
    !
    if (info/=0) then
      write (out,"(' zgesv returned ',i8)") info
      stop 'lapack_zgesv - cgelss failed'
    end if
    !
  end subroutine lapack_zgesv


   subroutine MLdiag_ulen_ark(n,a,d,ve)
      !
      integer(ik)  ::  p,q,p1,q1,irot,i,n2,kp,kp1,ip1,j,iq1,kq1,ii,n
      real(ark)     ::  a(n,n),d(n),ve(n,n)
      real(ark)     ::  sm,tresh,g,h,s,c,t,tau,theta,ff1,ff2,ff3,ff4
      real(ark)     ::  err
      real(ark),allocatable  ::  b(:),z(:)
      !
      err = small_
      !
      allocate(b(n+10),z(n+10))


 101  format(5e14.5)
      do 10 p=1,n
      do 10 q=1,n
      ve(p,q)=0.0_ark
      if(p.eq.q) ve(p,q)=1.0_ark
  10  continue
      do 99 p=1,n
      z(p)=0.0_ark
      d(p)=a(p,p)
      b(p)=d(p)
 99   continue
      irot=0
      do 50 i=1,50
      sm=0.0_ark
      n2=n-1
      do 30 p=1,n2
      kp=p+1
      do 30 q=kp,n
      sm=sm+abs(a(p,q))
  30  continue
      if(sm.le.err) goto 50
      tresh=0.0_ark
      if(i-4) 3,4,4
  3   tresh=0.2_ark*sm/(n*n)
  4     do 33 p1=1,n2
        kp1=p1+1
      do 33 q1=kp1,n
      g=100*abs(a(p1,q1))
      ff1=abs(d(p1)+g)
      ff2=abs(d(p1))
      ff3=abs(d(q1)+g)
      ff4=abs(d(q1))
      if(i.gt.4.and.(ff1.eq.ff2).and.ff3.eq.ff4) goto 7
      ff1=abs(a(p1,q1))
        if(ff1.le.tresh) goto 33
      h=d(q1)-d(p1)
      ff1=abs(h)+g
      ff2=abs(h)
      if(ff1.ne.ff2) goto 13
      t=a(p1,q1)/h
      goto 6
  13    theta=0.5_ark*h/a(p1,q1)
        t=1._ark/(abs(theta)+sqrt(1._ark+theta*theta))
      if(theta) 5,6,6
  5     t=-t
  6     c=1._ark/sqrt(1._ark+t*t)
        s=t*c
      tau=s/(1._ark+c)
      h=t*a(p1,q1)
      z(p1)=z(p1)-h
      z(q1)=z(q1)+h
      d(p1)=d(p1)-h
      d(q1)=d(q1)+h
      a(p1,q1)=0.0_ark
      ip1=p1-1
        do 20 j=1,ip1
        g=a(j,p1)
        h=a(j,q1)
        a(j,p1)=g-s*(h+g*tau)
        a(j,q1)=h+s*(g-h*tau)
  20      continue
        iq1=q1-1
        do 21 j=kp1,iq1
        g=a(p1,j)
        h=a(j,q1)
        a(p1,j)=g-s*(h+g*tau)
        a(j,q1)=h+s*(g-h*tau)
  21    continue
        kq1=q1+1
        do 26 j=kq1,n
        g=a(p1,j)
        h=a(q1,j)
        a(p1,j)=g-s*(h+g*tau)
        a(q1,j)=h+s*(g-h*tau)
  26    continue
          do 29 j=1,n
        g=ve(j,p1)
        h=ve(j,q1)
        ve(j,p1)=g-s*(h+g*tau)
        ve(j,q1)=h+s*(g-h*tau)
  29      continue
        irot=irot+1
  7     a(p1,q1)=0.0_ark
  33    continue
        do 44 ii=1,n
      d(ii)=b(ii)+z(ii)
      b(ii)=d(ii)
      z(ii)=0.0_ark
  44  continue
  50  continue

  deallocate(b,z)

  end subroutine MLdiag_ulen_ark



   subroutine MLdiag_ulen_rk(n,a,d,ve)
      !
      integer(ik)  ::  p,q,p1,q1,irot,i,n2,kp,kp1,ip1,j,iq1,kq1,ii,n
      real(rk)     ::  a(n,n),d(n),ve(n,n)
      real(rk)     ::  sm,tresh,g,h,s,c,t,tau,theta,ff1,ff2,ff3,ff4
      real(rk)     ::  err
      real(rk),allocatable  ::  b(:),z(:)
      !
      err = small_
      !
      allocate(b(n+10),z(n+10))


 101  format(5e14.5)
      do 10 p=1,n
      do 10 q=1,n
      ve(p,q)=0.0_rk
      if(p.eq.q) ve(p,q)=1.0_rk
  10  continue
      do 99 p=1,n
      z(p)=0.0_rk
      d(p)=a(p,p)
      b(p)=d(p)
 99   continue
      irot=0
      do 50 i=1,50
      sm=0.0_rk
      n2=n-1
      do 30 p=1,n2
      kp=p+1
      do 30 q=kp,n
      sm=sm+abs(a(p,q))
  30  continue
      if(sm.le.err) goto 50
      tresh=0.0_rk
      if(i-4) 3,4,4
  3   tresh=0.2_rk*sm/(n*n)
  4     do 33 p1=1,n2
        kp1=p1+1
      do 33 q1=kp1,n
      g=100*abs(a(p1,q1))
      ff1=abs(d(p1)+g)
      ff2=abs(d(p1))
      ff3=abs(d(q1)+g)
      ff4=abs(d(q1))
      if(i.gt.4.and.(ff1.eq.ff2).and.ff3.eq.ff4) goto 7
      ff1=abs(a(p1,q1))
        if(ff1.le.tresh) goto 33
      h=d(q1)-d(p1)
      ff1=abs(h)+g
      ff2=abs(h)
      if(ff1.ne.ff2) goto 13
      t=a(p1,q1)/h
      goto 6
  13    theta=0.5_rk*h/a(p1,q1)
        t=1._rk/(abs(theta)+sqrt(1._rk+theta*theta))
      if(theta) 5,6,6
  5     t=-t
  6     c=1._rk/sqrt(1._rk+t*t)
        s=t*c
      tau=s/(1._rk+c)
      h=t*a(p1,q1)
      z(p1)=z(p1)-h
      z(q1)=z(q1)+h
      d(p1)=d(p1)-h
      d(q1)=d(q1)+h
      a(p1,q1)=0.0_rk
      ip1=p1-1
        do 20 j=1,ip1
        g=a(j,p1)
        h=a(j,q1)
        a(j,p1)=g-s*(h+g*tau)
        a(j,q1)=h+s*(g-h*tau)
  20      continue
        iq1=q1-1
        do 21 j=kp1,iq1
        g=a(p1,j)
        h=a(j,q1)
        a(p1,j)=g-s*(h+g*tau)
        a(j,q1)=h+s*(g-h*tau)
  21    continue
        kq1=q1+1
        do 26 j=kq1,n
        g=a(p1,j)
        h=a(q1,j)
        a(p1,j)=g-s*(h+g*tau)
        a(q1,j)=h+s*(g-h*tau)
  26    continue
          do 29 j=1,n
        g=ve(j,p1)
        h=ve(j,q1)
        ve(j,p1)=g-s*(h+g*tau)
        ve(j,q1)=h+s*(g-h*tau)
  29      continue
        irot=irot+1
  7     a(p1,q1)=0.0_rk
  33    continue
        do 44 ii=1,n
      d(ii)=b(ii)+z(ii)
      b(ii)=d(ii)
      z(ii)=0.0_rk
  44  continue
  50  continue



  deallocate(b,z)

  end subroutine MLdiag_ulen_rk


end module lapack
