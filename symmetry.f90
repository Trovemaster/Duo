!
module symmetry
  use accuracy

  implicit none
  public SymmetryInitialize,sym,correlate_to_Cs
 

  type  ScIIT
     integer(ik)          :: Noper     ! Number of operations in the CII operator
     integer(ik)          :: Nzeta     ! Number of symmetric elements taking into account the degeneracies 
     real(ark),pointer    :: ioper(:)  ! the operation number in the MS group
     integer(ik),pointer  :: coeff(:)  ! coefficients of the CII operator
     integer(ik),pointer  :: izeta(:)  ! symmetry indentification as a eigenvalues of the CII operator
  end type ScIIT


  type  SymmetryT
     character(len=cl)    :: group = 'C' ! The symmetry group 
     integer(ik)          :: NrepresCs = 2  ! Number of irreduc. represent. for Cs(M)
     integer(ik)          :: Nrepresen = 1  ! Number of irreduc. represent.
     integer(ik)          :: Noper     = 1  ! Number of operations
     integer(ik)          :: Nclasses  = 1  ! Number of classes
     integer(ik),pointer  :: Nelements(:)   ! Number of elements in a class
     integer(ik),pointer  :: characters(:,:)! Character table
     type(SrepresT),pointer :: irr(:,:)     ! irreducible representaion 
     integer(ik),pointer  :: degen(:)       ! degeneracy
     character(len=3),pointer  :: label(:)  ! The symmetry label 
     integer(ik)          :: Maxdegen  = 1  ! Maximal degeneracy order
     integer(ik),pointer  :: igenerator(:)  ! address of the class generator in the sym%Ngroup list
     type(ScIIT)          :: CII            ! the elements of the CII operator 
     real(ark),pointer    :: euler(:,:)     ! rotational angles equivalent to the group operations
  end type SymmetryT

  type  SrepresT
     real(ark),pointer  :: repres(:,:)      ! matrix representation of the group 
  end type SrepresT


  type(SymmetryT) , save  :: sym


contains 


  subroutine SymmetryInitialize(sym_group)
  character(len=cl),intent(in) :: sym_group
  integer(ik):: alloc,iclass,gamma,ioper,ielem
  real(ark)  :: o,p2,p3
  !   
  sym%group=sym_group
  !
  select case(trim(sym_group))

  case("C(M)","C")

    sym%Nrepresen=1
    sym%Noper=1
    sym%Nclasses=1
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters(1,1)=1
    sym%degen=(/1/)
    sym%Nelements=(/1/)
    sym%label=(/'A'/)

    do gamma = 1,sym%Nrepresen
      ioper = 0
      do iclass = 1,sym%Nclasses
        do ielem = 1,sym%Nelements(iclass)
          !
          ioper = ioper + 1
          !
          allocate(sym%irr(gamma,ioper)%repres(sym%degen(gamma),sym%degen(gamma)),stat=alloc)
          !
          if (sym%degen(gamma)==1) then
             sym%irr(gamma,ioper)%repres(1,1)= sym%characters(gamma,iclass)
          endif 
          !
        enddo 
      enddo 
    enddo 


  case("CS(M)","CS")

    sym%Nrepresen=2
    sym%Noper=2
    sym%Nclasses=2
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &
                                 (/ 1, 1,  & 
                                    1,-1  /),(/2,2/))
    sym%degen=(/1,1/)
    sym%Nelements=(/1,1/)
    sym%label=(/'A''','A"'/)

    call irr_allocation

  case("C2V(M)","C2V")

    sym%Nrepresen=4
    sym%Noper=4
    sym%Nclasses=4
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &
                                  (/ 1, 1, 1, 1, &   ! A1
                                     1, 1,-1,-1, &   ! A2
                                     1,-1,-1, 1, &   ! B1
                                     1,-1, 1,-1 /),(/4,4/)) ! B2
    sym%degen=(/1,1,1,1/)
    sym%Nelements=(/1,1,1,1/)
    sym%label=(/'A1','A2','B1','B2'/)
    !
    o  = 0.0_ark
    p2 = 0.5_ark*pi
    p3 = 1.5_ark*pi
    !
    sym%euler( 1,:) = 0
    sym%euler( 2,:) = (/pi,o,o /)
    sym%euler( 3,:) = (/o,pi,o/)
    sym%euler( 4,:) = (/p2,pi,p3/)
    !
    call irr_allocation


  case("C2H(M)","C2H")

    sym%Nrepresen=4
    sym%Noper=4
    sym%Nclasses=4
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &
                                  (/ 1, 1, 1, 1, &   ! Ag
                                     1, 1,-1,-1, &   ! Au
                                     1,-1, 1,-1, &   ! Bg
                                     1,-1,-1, 1 /),(/4,4/)) ! Bu
    sym%degen=(/1,1,1,1/)
    sym%Nelements=(/1,1,1,1/)
    sym%label=(/'Ag','Au','Bg','Bu'/)
    !
    call irr_allocation


  case("G4(M)","G4")

    sym%Nrepresen=4
    sym%Noper=4
    sym%Nclasses=4
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &       ! E (12)(34)  E*  (12)(34)*
                                  (/ 1, 1, 1, 1, &   ! A+
                                     1, 1,-1,-1, &   ! A-
                                     1,-1,-1, 1, &   ! B+
                                     1,-1, 1,-1 /),(/4,4/)) ! B-
    sym%degen=(/1,1,1,1/)
    sym%Nelements=(/1,1,1,1/)
    sym%label=(/'A+','A-','B+','B-'/)
    !
    call irr_allocation



  case("G4(EM)")

    sym%Nrepresen=8
    sym%Noper=8
    sym%Nclasses=8
    sym%CII%Noper = 0

    call simple_arrays_allocation

    sym%characters= reshape( &     ! E   a   b  ab   E' E'a E'b E'ab 
                                  (/ 1,  1,  1,  1,  1,  1,  1,  1, &   ! Ags
                                     1,  1, -1, -1,  1,  1, -1, -1, &   ! Aus
                                     1, -1, -1,  1,  1, -1, -1,  1, &   ! Bgs
                                     1, -1,  1, -1,  1, -1,  1, -1, &   ! Bus
                                     1,  1, -1, -1, -1, -1,  1,  1, &   ! Agd 
                                     1,  1,  1,  1, -1, -1, -1, -1, &   ! Aud
                                     1, -1,  1, -1, -1,  1, -1,  1, &   ! Bgd
                                     1, -1, -1,  1, -1,  1,  1, -1 /),(/8 ,8/)) ! Bud
                                     


    sym%degen=(/1,1,1,1,1,1,1,1/)
    sym%Nelements=(/1,1,1,1,1,1,1,1/)
    sym%label=(/'Ags','Aus','Bgs','Bus','Agd','Aud','Bgd','Bud'/)
    !
    call irr_allocation

  case default

    write(out,"('symmetry: undefined symmetry group ',a)") trim(sym_group)
    stop 'symmetry: undefined symmetry group '

  end select
  !
  if (20<sym%Nrepresen) then 
    !
    write(out,"('symmetry: number of elements in _select_gamma_ is too small: ',i8)") 20
    stop 'symmetry: size of _select_gamma_ is too small'
    !
  endif 
  !
  sym%maxdegen = maxval(sym%degen(:),dim=1)

  !
  ! store the address of the group generator from ioper = 1..Noper list 
  !
  ioper = 1
  !
  do iclass = 1,sym%Nclasses
    !
    sym%igenerator(iclass) = ioper
    ioper = ioper + sym%Nelements(iclass)
    !
  enddo
  
   


  contains

   subroutine simple_arrays_allocation

    integer(ik) :: alloc,nCII

    nCII = max(1,sym%CII%Noper)
    !
    allocate (sym%characters(sym%Nclasses,sym%Nrepresen),sym%irr(sym%Nrepresen,sym%Noper),&  
              sym%degen(sym%Nrepresen),sym%Nelements(sym%Nclasses),sym%label(sym%Nrepresen),&
              sym%igenerator(sym%Nclasses),&
              sym%CII%ioper(nCII),sym%CII%coeff(nCII),sym%euler(sym%Noper,3),stat=alloc)

    if (alloc/=0) stop 'simple_arrays_allocation - out of memory'
    !
    sym%CII%coeff = 0
    sym%CII%ioper = 1
    sym%euler = 0 
    !
   end subroutine simple_arrays_allocation



   subroutine irr_allocation

    integer(ik) :: gamma,ioper,iclass,ielem,alloc

    do gamma = 1,sym%Nrepresen
      ioper = 0
      do iclass = 1,sym%Nclasses
        do ielem = 1,sym%Nelements(iclass)
          !
          ioper = ioper + 1
          !
          allocate(sym%irr(gamma,ioper)%repres(sym%degen(gamma),sym%degen(gamma)),stat=alloc)
          !
          if (sym%degen(gamma)==1) then
             sym%irr(gamma,ioper)%repres(1,1)= sym%characters(gamma,iclass)
          endif 
          !
        enddo 
      enddo 
    enddo 


   if (alloc/=0) then
       write (out,"(' symmetryInitialize ',i9,' error trying to allocate symmetry')") alloc
      stop 'symmetryInitialize, symmetries - out of memory'
   end if


   end subroutine irr_allocation

 end subroutine symmetryInitialize


 function correlate_to_Cs(iparity,gu) result(isym)
   !
   implicit none 
   !
   integer(ik),intent(in) :: iparity,gu
   integer(ik) :: isym
     !
     if (sym%Nrepresen==sym%NrepresCs) then 
       !
       isym = iparity
       return
       !
     endif
     !
     if (iparity==1.and.gu==1) then 
       ! A1
       isym = 1
       !
     elseif (iparity==2.and.gu==1) then 
       ! B1
       isym = 3
       !
     elseif (iparity==1.and.gu==-1) then 
       ! B2
       isym = 4
       !
     else
       ! A2
       isym = 2
       !       
     endif
     ! 
  end function correlate_to_Cs


end module symmetry



