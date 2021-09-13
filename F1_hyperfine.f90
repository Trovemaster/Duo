module F1_hyperfine

    use accuracy
    use lapack
    use symmetry
    use diatom_module, only: basis, three_j, eigen, jmax, &
                            vibrational_totalroots, vibrational_contrfunc, vibrational_quantum_number, &
                            hfcc, I1, scaling_factor, &
                            r, d2dr, r2sc,z, grid, m1, m2, amass, hstep, action, &
                            GLOBAL_NUM_HFCC_OBJECT, abinitio, job, poten, bad_value, &
                            eigenT, basisT, quantaT, fieldT, linkT

    implicit none
    
    
    TYPE(basisT), ALLOCATABLE :: primitive_F1_basis(:)
    
    REAL(rk), ALLOCATABLE :: F1_list(:)
    
    INTEGER(ik) :: num_F1, num_primitive_F1_basis, num_represCs
    REAL(rk) :: F1_global_min, F1_global_max
    INTEGER(ik) :: alloc
    REAL(rk) :: J_global_min = 0.5_rk, J_global_max
    REAL(rk) :: J_min, J_max, MHz_to_wavenumber = 1.0E6_rk / vellgt, &
                wavenumber_to_MHz = vellgt / 1.0E6_rk

    TYPE(fieldT) :: FC_bF_field 
    TYPE(eigenT), allocatable :: eigen_all_F1(:,:)
    LOGICAL :: F1_hyperfine_gridvalue_allocated  = .false.
    
contains

    subroutine F1_hyperfine_structrure
        ! Top level 
        implicit none

        INTEGER(ik) :: index_F1, index_represCs
        INTEGER(ik) :: Ndimen_F1, Nlevels_F1
        REAL(rk), ALLOCATABLE :: eigen_value_F1(:), &
                                    parity_conserved_F1_matrix(:, :), &
                                    primitive_F1_hyperfine_matrix(:,:), &
                                    transformation_matrix(:,:)
        
        
        ! wavenumber_to_MHz = wavenumber_to_MHz

        J_global_max = jmax
        if (J_global_max < I1) then
            write (out, *) "J_global_max < I1"
            stop
        end if

        F1_global_max = J_global_max - I1

        if (mod(nint(F1_global_max*2.0_rk), 2) == 0) then
            F1_global_min = 0.0_rk
        else
            F1_global_min = 0.5_rk
        end if

        ! The basis of J always starts from the mininum of J
        if (mod(nint(J_global_max*2.0_rk), 2) == 0) then
            J_global_min = 0.0_rk
        else
            J_global_min = 0.5_rk
        end if

        num_F1 = nint(F1_global_max - F1_global_min) + 1
        
        
        CALL map_F1_hyperfine_fields_onto_grid(4)
        
        CALL construct_F1_hyperfine_constant_field_matrix        
        

        ALLOCATE(F1_list(num_F1))
        ALLOCATE(primitive_F1_basis(num_F1))

        CALL construct_primitive_F1_basis
        
        DEALLOCATE(F1_list)

        ALLOCATE(eigen_all_F1(num_F1, sym%NrepresCs))

        do index_F1 = 1, num_F1
            
            Ndimen_F1 = primitive_F1_basis(index_F1)%Ndimen
            ALLOCATE(primitive_F1_hyperfine_matrix(Ndimen_F1, Ndimen_F1))
            CALL construct_primitive_F1_hyperfine_matrix( &
                index_F1, primitive_F1_hyperfine_matrix)

            do index_represCs = 1, sym%NrepresCs

                Nlevels_F1 = Ndimen_F1 / sym%NrepresCs

                ALLOCATE(eigen_value_F1(Nlevels_F1))
                ALLOCATE(parity_conserved_F1_matrix(Nlevels_F1, Nlevels_F1))
                ALLOCATE(transformation_matrix(Ndimen_F1, Nlevels_F1))

                ALLOCATE(eigen_all_F1(index_F1, index_represCs)%quanta(Nlevels_F1))
                ALLOCATE(eigen_all_F1(index_F1, index_represCs)%val(Nlevels_F1))
                ALLOCATE(eigen_all_F1(index_F1, index_represCs)%vect(Ndimen_F1, Nlevels_F1))

                CALL construct_parity_conserved_F1_matrix( &
                    index_F1, index_represCs, primitive_F1_hyperfine_matrix, &
                    parity_conserved_F1_matrix, transformation_matrix)
                ! Now, parity_conserved_F1_matrix(:, :) is the sum of
                ! parity-conserved hyperfine matrix under the rovibronic basis set
                ! and the eigenvalues of rovibronic matrix.
                ! transformation_matrix(:, :) holds the rovibronic eigenvectors of F1
                ! viz Phi

                ! Diagonalize parity conserved F1 matrix
                CALL lapack_syev(parity_conserved_F1_matrix, eigen_value_F1)
                ! Now, parity_conserved_F1_matrix(:, :) is the eigenvector matrix
                ! under the rovibronic basis set, viz U.
                
                ! Psi = Phi * U
                transformation_matrix = matmul(transformation_matrix, &
                                               parity_conserved_F1_matrix)
                ! Now, transformation_matrix(:, :) is the eigenvector matrix
                ! under the Hund's case(a_beta) basis set.
                
                eigen_all_F1(index_F1, index_represCs)%Nlevels = Nlevels_F1
                eigen_all_F1(index_F1, index_represCs)%Ndimen = Ndimen_F1
                eigen_all_F1(index_F1, index_represCs)%val = eigen_value_F1
                eigen_all_F1(index_F1, index_represCs)%vect = transformation_matrix

                DEALLOCATE(transformation_matrix)
                DEALLOCATE(parity_conserved_F1_matrix)
                DEALLOCATE(eigen_value_F1)                
            end do

            DEALLOCATE(primitive_F1_hyperfine_matrix)
        end do
        
        ! DO index_F1 = 1, num_F1
        !     DEALLOCATE (primitive_F1_basis(index_F1)%icontr)
        ! end do
        
        ! DEALLOCATE(primitive_F1_basis)       

    end subroutine F1_hyperfine_structrure

    subroutine map_F1_hyperfine_fields_onto_grid(iverbose)
        !
        use functions,only : define_analytical_field
        !
        character(len=130)     :: my_fmt  !text variable containing formats for reads/writes
        ! 
        integer(ik),intent(in) :: iverbose
        !
        integer(ik)             :: ngrid,alloc,j,nsub,Nmax,iterm,nterms,i,ipotmin,istate,jstate,itotal
        integer(ik)             :: ifterm,iobject,ifield
        real(rk)                :: rmin, rmax, re, alpha, h,sc,h12,scale,check_ai
        real(rk),allocatable    :: f(:)
        !
        integer             :: np   ! tmp variable for extrapolation
        real(kind=rk)       :: x1, x2, y1, y2, aa, bb ! tmp variables used for extrapolation
        !
        real(rk),allocatable    :: spline_wk_vec(:) ! working vector needed for spline interpolation
        real(rk),allocatable    :: spline_wk_vec_B(:) ! working vector needed for spline interpolation
        real(rk),allocatable    :: spline_wk_vec_C(:) ! working vector needed for spline interpolation
        real(rk),allocatable    :: spline_wk_vec_D(:) ! working vector needed for spline interpolation
        real(rk),allocatable    :: spline_wk_vec_E(:) ! working vector needed for spline interpolation
        real(rk),allocatable    :: spline_wk_vec_F(:) ! working vector needed for spline interpolation
        real(rk),allocatable    :: xx(:), yy(:), ww(:)! tmp vectors for extrapolation
        type(linkT),allocatable :: link_(:) ! tmp link values for extrapolation
        real(rk) :: yp1, ypn, Vtop, beta, vmin, DeltaR
        integer(ik) :: imin
        !
        type(fieldT),pointer      :: field
        !
        ngrid = grid%npoints

        if (iverbose>=4) call TimerStart('Grid representaions')
        !
        if (iverbose>=3) write(out,'("Generate a grid representation for all Hamiltonian terms")')
        !
        ! in case of fitting we will need this object to store the parameters to refine
        !
        if (action%fitting) then          
            itotal = 0
            ifterm = 0          
        endif
        !
        object_loop: do iobject = 1, GLOBAL_NUM_HFCC_OBJECT
            !
            Nmax = hfcc(iobject)%num_field
            !
            ! each field type consists of Nmax terms
            !
            do iterm = 1,Nmax
                !
                field => hfcc(iobject)%field(iterm)
                !
                if (.not. f1_hyperfine_gridvalue_allocated) then 
                    !
                    allocate(field%gridvalue(ngrid),stat=alloc)
                    call ArrayStart(trim(field%type),alloc,ngrid,kind(field%gridvalue))
                    !
                endif
                !
                select case(trim(field%type))
                !
                case("GRID")
                    nterms = field%Nterms
                    !
                    ! Lorenzo Lodi, 13 February 2014
                    ! This section will add extrapolated points at short and long bond lengths
                    ! for tabulated `GRID' functions
                    ! and only for fitting 
                    if( field%grid(1) > rmin ) then
                        if (iverbose>=4) write(out, '(/,A)') 'Extrapolating at short bond length curve ' // trim(field%name) // &
                                            ' (class ' // trim(field%class) // ')'
                    
                    
                            np = 20              ! I always add `np' short bond length guide points
                            ! I go one step `beyond' rmin to minimize interpolating artifacts at the end of the interval
                    
                            x1=field%grid(1)  ; y1 =field%value(1)
                            x2=field%grid(2)  ; y2 =field%value(2)
                            allocate( xx(nterms+np), yy(nterms+np),ww(nterms+np),stat=alloc)
                            call ArrayStart('extrap_tmp',alloc,size(xx),kind(xx))
                            call ArrayStart('extrap_tmp',alloc,size(ww),kind(ww))
                            call ArrayStart('extrap_tmp',alloc,size(yy),kind(yy))
                            allocate( link_(nterms+np),stat=alloc)
                            if (alloc/=0) then 
                            write(out,"('allocation problem: extrap_tmp')")
                            stop 'allocation problem'
                            endif
                            xx=0._rk
                            yy=0._rk
                            do i=1, np
                            xx(i) = rmin + (field%grid(1)-rmin)*real(i-2,rk) / real( np-1, rk)
                            enddo
                            !
                        select case(field%class)
                            !
                        case ("POTEN")
                            !
                            if (iverbose>=4) write(out, '(A, I20)') 'Using A + B/r extrapolation; points added = ', np
                            bb = -x1*x2*(y1-y2)/(x1-x2)
                            aa = y1 - bb/x1
                            do i=1, np
                            yy(i) = aa + bb/xx(i)
                            enddo
                            !
                        case ("DIPOLE")
                            if (iverbose>=4) write(out, '(A, I20)') 'Using A*r + B*r^2 extrapolation; points added = ', np
                            bb = (x2*y1- x1*y2) / (x1*x2 *(x1-x2) )
                            aa = (-bb*x1**2 + y1 ) / x1
                            do i=1, np
                                yy(i) = aa*xx(i) + bb*xx(i)**2
                            enddo
                            !
                        case default  ! linear extrapolation using the first two points to r = rmin - 1/(np-1)
                            if (iverbose>=4) write(out, '(A, I20)') 'Using linear extrapolation; points added = ', np
                            do i=1, np
                                yy(i) = y1 + (xx(i)-x1) * (y1-y2)/(x1-x2)
                            enddo
                        end select
                        !
                        do i=np+1, np+nterms
                            xx(i) =  field%grid(i-np)
                            yy(i) =  field%value(i-np)
                        enddo
                        !
                        ww(np+1:) = field%weight(1:) ; ww(1:np) = 0
                        !
                        link_(1:np)%iobject = 0
                        link_(1:np)%ifield = 0
                        link_(1:np)%iparam = 0
                        link_(np+1:) = field%link(1:)
                        !
                        nterms = np+nterms
                        field%Nterms = nterms
                        !
                        call ArrayMinus(trim(field%type),size(field%value),kind(field%value))
                        call ArrayMinus(trim(field%type),size(field%grid),kind(field%grid))
                        call ArrayMinus(trim(field%type),size(field%weight),kind(field%weight))
                        !
                        deallocate(field%grid, field%value, field%weight, field%forcename,stat=alloc)
                        !
                        deallocate(field%link)
                        !
                        allocate(field%value(nterms),field%forcename(nterms),field%grid(nterms),field%weight(nterms),stat=alloc)
                        allocate(field%link(nterms),stat=alloc)
                        !
                        call ArrayStart(trim(field%type),alloc,nterms,kind(field%value))
                        call ArrayStart(trim(field%type),alloc,nterms,kind(field%grid))
                        call ArrayStart(trim(field%type),alloc,nterms,kind(field%weight))
                        !
                        field%link(:) = link_(:)
                        field%grid  = xx
                        field%value = yy
                        field%forcename= 'dummy'
                        field%weight = ww
                        deallocate(xx, yy, ww,link_)
                        call ArrayStop('extrap_tmp')
                        !
                    endif
                    !
                    !*************end of short bond length extrapolation **********************************************
                    nterms = field%Nterms
                    if( field%grid(nterms) < rmax) then
                        if (iverbose>=4) write(out, '(/,A)') 'Extrapolating at long bond length curve ' // trim(field%name) // &
                                        ' (class ' // trim(field%class) // ')'
                        !
                        np = 20    ! I always add `np' long bond length guide points
                        ! I go one step `beyond' rmax to minimize interpolating artifacts at the end of the interval
                        !
                        x1=field%grid(nterms-1)  ; y1 =field%value(nterms-1)
                        x2=field%grid(nterms)    ; y2 =field%value(nterms)
                        allocate( xx(nterms+np), yy(nterms+np),ww(nterms+np),stat=alloc)
                        call ArrayStart('extrap_tmp',alloc,size(xx),kind(xx))
                        call ArrayStart('extrap_tmp',alloc,size(ww),kind(ww))
                        call ArrayStart('extrap_tmp',alloc,size(yy),kind(yy))
                        allocate( link_(nterms+np),stat=alloc)
                        if (alloc/=0) then 
                            write(out,"('allocation problem: extrap_tmp')")
                            stop 'allocation problem'
                        endif
                        !
                        xx=0._rk
                        yy=0._rk
                        do i=1, nterms
                            xx(i) =  field%grid(i)
                            yy(i) =  field%value(i)
                        enddo
                        !
                        do i=nterms+1, nterms+np
                            xx(i) = field%grid(nterms) + (rmax-field%grid(nterms))*real(i-nterms,rk) / real( np-1, rk)
                        enddo
                        !
                        select case(field%class)
                    
                        case ("POTEN")
                            if (iverbose>=4) write(out, '(A, I20)') 'Using De + A/r^6 extrapolation; points added = ', np
                            bb = (x1*x2)**6 * (y2-y1) / (x1**6 - x2**6)
                            aa = y1 - bb/x1**6
                            !       write(out, '(A, 2F20.4)') 'Dissociation set to = ', aa !, bb
                            do i=nterms+1,nterms+np
                            yy(i) = aa + bb / xx(i)**6
                            enddo
                            !
                        case ("DIPOLE")
                            if (iverbose>=4) write(out, '(A, I20)') 'Using A/r^2 + B/r^3 extrapolation; points added = ', np
                            bb = x1*x2*(x1**2 * y1 - x2**2 * y2)/(x2-x1)
                            aa = (y1*x1**3 - bb)/x1
                            do i=nterms+1,nterms+np
                            yy(i) = aa/xx(i)**2 + bb / xx(i)**3
                            enddo
                            !
                        case default  ! extrapolation using the last two points
                            if (iverbose>=4) write(out, '(A, I20)') 'Using linear extrapolation; points added = ', np
                            do i=nterms+1,nterms+np
                                yy(i) = y1 + (xx(i)-x1) * (y1-y2)/(x1-x2)
                            enddo
                        end select
                        !
                        ww(1:nterms) = field%weight(1:nterms) ; ww(nterms+1:) = 0
                        !
                        link_(1:np)%iobject = 0
                        link_(1:np)%ifield = 0
                        link_(1:np)%iparam = 0
                        link_(np+1:) = field%link(1:)
                        !
                        nterms = np+nterms
                        field%Nterms = nterms
                        !
                        call ArrayMinus(trim(field%type),size(field%value),kind(field%value))
                        call ArrayMinus(trim(field%type),size(field%grid),kind(field%grid))
                        call ArrayMinus(trim(field%type),size(field%weight),kind(field%weight))
                        !
                        deallocate(field%grid, field%value, field%weight, field%forcename,stat=alloc)
                        !
                        deallocate(field%link)
                        !
                        allocate(field%value(nterms),field%forcename(nterms),field%grid(nterms),field%weight(nterms),stat=alloc)
                        allocate(field%link(nterms),stat=alloc)
                        !
                        call ArrayStart(trim(field%type),alloc,nterms,kind(field%value))
                        call ArrayStart(trim(field%type),alloc,nterms,kind(field%grid))
                        call ArrayStart(trim(field%type),alloc,nterms,kind(field%weight))
                        !
                        field%link(:) = link_(:)
                        field%grid  = xx
                        field%value = yy
                        field%forcename= 'dummy'
                        field%weight = ww
                        deallocate(xx, yy, ww,link_)
                        call ArrayStop('extrap_tmp')
                        !
                    endif
                    !
                    !*************end of long bond length extrapolation **********************************************
                    !
                    ! build spline interpolant
                    select case (field%interpolation_type)
                    case default
                        write(out,"(a)") "Unrecognized interpolation type ", field%interpolation_type
                        stop "illegal iobject"
                        ! -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- .
                    case("CUBICSPLINES")  ! cubic natural spline
                        !
                        if (iverbose>=6) write(out, '(A, I20)') 'Interpolating with cubic natural splines'
                        allocate(spline_wk_vec(nterms),stat=alloc)
                        call ArrayStart('spline_wk_vec-field',alloc,ngrid,kind(spline_wk_vec))
                        !
                        yp1= 0._rk ; ypn =0._rk  ! 1nd derivatives at the first and last point (ignored)
                        call spline(field%grid,field%value,field%Nterms,yp1,ypn,spline_wk_vec)
                        !
                        !$omp parallel do private(i) schedule(guided)
                        do i=1,ngrid
                        ! evaluate spline interpolant
                            call splint(field%grid,field%value,spline_wk_vec,field%Nterms,r(i),field%gridvalue(i))
                        enddo
                        !$omp end parallel do
                        !
                        deallocate(spline_wk_vec)
                        call ArrayStop('spline_wk_vec-field')
                        !
                        ! -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- .
                    case("QUINTICSPLINES")  ! quintic spline
                        ! Lorenzo Lodi, 7 April 2016
                        ! Tests show that the present quintic spline interpolation does not give better interpolation results
                        ! than cubic splines. Furthermore the error seems to decrease as npoints**(-4)
                        ! instead of npoints**(-6). The reason for this behaviour is yet unknown, but possibly related 
                        ! to a bad choice of the boundary conditions. 
                        !
                        ! FOR NOW USE OF QUINTIC SPLINE INTERPOLATION IS NOT RECOMMENDED.
                        ! USE CUBIC SPLINES INSTEAD
                        !
                        if (iverbose>=6) write(out, '(A, I20)') 'Interpolating with quintic splines'
                        nterms=field%Nterms
                        allocate(spline_wk_vec_B(nterms),stat=alloc); 
                        call ArrayStart('spline_wk_vec_B',alloc,nterms,kind(spline_wk_vec))
                        allocate(spline_wk_vec_C(nterms),stat=alloc); 
                        call ArrayStart('spline_wk_vec_C',alloc,nterms,kind(spline_wk_vec))
                        allocate(spline_wk_vec_D(nterms),stat=alloc); 
                        call ArrayStart('spline_wk_vec_D',alloc,nterms,kind(spline_wk_vec))
                        allocate(spline_wk_vec_E(nterms),stat=alloc); 
                        call ArrayStart('spline_wk_vec_E',alloc,nterms,kind(spline_wk_vec))
                        allocate(spline_wk_vec_F(nterms),stat=alloc); 
                        call ArrayStart('spline_wk_vec_F',alloc,nterms,kind(spline_wk_vec))
                        !
                        call QUINAT(nterms, field%grid,field%value, spline_wk_vec_B, spline_wk_vec_C, &
                                                    spline_wk_vec_D, spline_wk_vec_E, spline_wk_vec_F)
                        !
                        !$omp parallel do private(i) schedule(guided)
                        do i=1,ngrid    ! evaluate spline interpolant
                        !        call splint(field%grid,field%value,spline_wk_vec,field%Nterms,r(i),field%gridvalue(i))
        
                            call splint_quint(field%grid,field%value,nterms,r(i),field%gridvalue(i), spline_wk_vec_B, &
                                                spline_wk_vec_C, spline_wk_vec_D, spline_wk_vec_E, spline_wk_vec_F)
                        enddo
                        !$omp end parallel do
                        !
                        deallocate(spline_wk_vec_B) ; call ArrayStop('spline_wk_vec_B')
                        deallocate(spline_wk_vec_C) ; call ArrayStop('spline_wk_vec_C')
                        deallocate(spline_wk_vec_D) ; call ArrayStop('spline_wk_vec_D')
                        deallocate(spline_wk_vec_E) ; call ArrayStop('spline_wk_vec_E')
                        deallocate(spline_wk_vec_F) ; call ArrayStop('spline_wk_vec_F')
                        ! -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- . -- .
                        !
                    end select
                    !
                    ! for dummy fields not used in fittings
                    !
                case("DUMMY")
                    nterms = field%Nterms
                    field%gridvalue = 0._rk
                    !
                case default
                    !
                    call define_analytical_field(field%type,field%analytical_field)
                    !
                    !$omp parallel do private(i) schedule(guided)
                    do i=1,ngrid
                        !
                        field%gridvalue(i) = field%analytical_field(r(i),field%value)
                        !
                    enddo
                    !$omp end parallel do
                    !
                    ! counting the total number of parameters when the fitting is on
                    !
                    if (action%fitting) then
                        itotal = itotal + field%Nterms
                    endif
                    !
                end select
                    !
                field%gridvalue =  field%gridvalue*field%factor
            enddo
        enddo object_loop
        !
        ! change the status of  f1_hyperfine_gridvalue_allocated to true (allocated)
        ! to prevent reallocation next time this subroutine is called (e.g. from the refinement)
        f1_hyperfine_gridvalue_allocated = .true.
        !
        ! Now morph the objects by applying the morphing function to the corresponding ab initio field if neccessary
        !
        ifield = 0 
        
        object_loop2: do iobject = 1,GLOBAL_NUM_HFCC_OBJECT
            !
            Nmax = hfcc(iobject)%num_field
            !
            ! each field type consists of Nmax terms
            !
            do iterm = 1,Nmax
                !
                field => hfcc(iobject)%field(iterm)
                !
                ifield = ifield + 1
                !
                ! Shift values up or down by specified constant value
                if (field%adjust) then
                    field%gridvalue = field%gridvalue + field%adjust_val
                endif
                !
                ! Introduce morphing (not for abinitio)
                !
                if (field%morphing.and.iobject/=GLOBAL_NUM_HFCC_OBJECT-2) then
                    !
                    ! check if ai field was defined 
                    !
                    check_ai = sum((abinitio(ifield)%gridvalue)**2)
                    !
                    if (check_ai<small_) then 
                        !
                        write(out,"('Error: Corresponding ab initio field is undefined when using MORPHING for ',a)")&
                              trim(field%name)
                        stop 'ab initio field is undefined while using MORPHING'
                        !
                    endif
                    ! 
                    field%gridvalue = field%gridvalue*abinitio(ifield)%gridvalue
                    ! 
                endif
                !
                ! Generate weights if an analytical expression is given
                !
                if (trim(field%weighting%wtype)=="PS1997") then 
                !
                istate = field%iref
                beta = field%weighting%alpha
                Vtop = field%weighting%Vtop
                imin = (field%grid(1)-rmin) / real( np-1, rk)
                Vmin = minval(poten(istate)%gridvalue(:))
                !
                DeltaR = (rmax-rmin)/real(ngrid-1,rk)
                !
                do i=1,field%Nterms
                    !
                    np = nint(( field%grid(i) - rmin )/ DeltaR)
                    !
                    np = max(min(ngrid,np),1)
                    !
                    field%weight(i) = ( tanh(-beta*( poten(istate)%gridvalue(np) - Vmin - Vtop ) ) & 
                                        +1.000020000200002_rk )/2.000020000200002_rk
                    !
                enddo
                !
                endif  
                !
                !
                ! transform from the MOLRPO to Duo representation
                !
                if (field%molpro) then
                
                    if (iverbose>=4) then 
                        write(out,'(/,a)') 'Transforming '//trim(field%class)//' '//trim(field%name)//&
                        ' from MOLPRO (Cartesian) to Duo (Spherical)'            
                        write(out,'(a)') '  Cartesian (complex):'
                        write(out,'(a,t20,i4,1x,t30,i4)') '  |Lambda|:',field%lambda,field%lambdaj
                        if (iobject==2) write(out,'(a,t18,f8.1,t28,f8.1)') '  Sigma:',field%sigmai,field%sigmaj 
                        write(out,'(a,t20,i4,a,t30,i4,a)') '  <a|Lz|b>:',field%ix_lz_y,'i',field%jx_lz_y,'i'
                        write(out,'(a,t18,f8.1,a,f8.1,a)') '  factor',real(field%complex_f),'+',aimag(field%complex_f),'i'
                    endif
                
                    if (job%legacy_version) then  
                        !
                        call molpro_duo_old_2018(field)
                        !
                    else 
                        !
                        call molpro_duo(field)
                        !
                    endif
                
                    if (iverbose>=4) then 
                        write(out,'(a)') '  Spherical (real):'
                        write(out,'(a,t20,i4,1x,t30,i4)') '  Lambda:',field%lambda,field%lambdaj
                        if (iobject==2) write(out,'(a,t18,f8.1,t28,f8.1)') '  Sigma:',field%sigmai,field%sigmaj 
                    endif
                endif
             
            enddo
           
        enddo object_loop2
        !
        ! if (action%fitting) then
        !     fitting%parmax = itotal
        ! endif
        !
        if (iverbose>=3) write(out,'("...done!"/)')
        !
        if (iverbose>=4) call TimerStop('Grid representaions')

      contains 
        !
        subroutine molpro_duo(field)
           !
            use lapack,only : lapack_zheev     
            !
            type(fieldT),intent(inout) :: field
            integer(ik) :: ix_lz_y,jx_lz_y,iroot,jroot,il_temp,ngrid
            complex(rk) :: a(2,2),b(2,2),coupling(2,2),f_t(2,2),b0(2,2),c
            real(rk)    :: lambda_i(2),lambda_j(2)
                !
                ngrid = grid%npoints
                !
                ix_lz_y = field%ix_lz_y
                jx_lz_y = field%jx_lz_y
                !
                a = 0 ; b = 0
                !
                a(1,1) =  cmplx(1.0_rk,0.0_rk, kind=rk)
                b(1,1) =  cmplx(1.0_rk,0.0_rk, kind=rk)
                a(2,2) =  cmplx(1.0_rk,0.0_rk, kind=rk)
                b(2,2) =  cmplx(1.0_rk,0.0_rk, kind=rk)
                !
                iroot = 1
                jroot = 1
                !
                if (ix_lz_y/=0) then
                  a = 0 
                  a(1,2) = cmplx(0.0_rk,ix_lz_y , kind=rk)
                  a(2,1) = cmplx(0.0_rk,-ix_lz_y, kind=rk)
                  !
                  call lapack_zheev(a,lambda_i)
                  !
                  ! swap to have the first root positive 
                  !
                  f_t = a
                  a(:,1) = f_t(:,2)
                  a(:,2) = f_t(:,1)
                  !
                  il_temp = lambda_i(2)
                  lambda_i(2) = lambda_i(1)
                  lambda_i(1) = il_temp
                  !
                  field%lambda = nint(lambda_i(1))
                  !
                  a = a*cmplx(0.0_rk,1.0_rk, kind=rk)
                  !
                elseif (poten(field%istate)%parity%pm==-1) then 
                  !
                  a(1,1) =  cmplx(0.0_rk,1.0_rk, kind=rk)
                  !
                endif
                !
                if (jx_lz_y/=0) then
                  b = 0 
                  b(1,2) = cmplx(0.0_rk,jx_lz_y , kind=rk)
                  b(2,1) = cmplx(0.0_rk,-jx_lz_y, kind=rk)
                  !
                  b0 = b
                  !
                  call lapack_zheev(b,lambda_j)
                  !
                  ! swap to have the first root positive 
                  !
                  f_t = b
                  b(:,1) = f_t(:,2)
                  b(:,2) = f_t(:,1)
                  !
                  il_temp = lambda_j(2)
                  lambda_j(2) = lambda_j(1)
                  lambda_j(1) = il_temp
                  !
                  field%lambdaj = nint(lambda_j(1))
                  !
                  b = b*cmplx(0.0_rk,1.0_rk, kind=rk)
                  !
                  f_t = matmul( conjg(transpose(b)),matmul(b0,(b)) )
                  !
                elseif (poten(field%jstate)%parity%pm==-1) then 
                  !
                  b(1,1) =  cmplx(0.0_rk,1.0_rk, kind=rk)
                  !
                endif
                !
                ! Check the selection rules
                select case(trim(field%class))
                  !
                case('SPINORBIT','ABINITIO-SPINORBIT')
                  !
                  if ((nint(field%sigmaj-field%sigmai))/=(field%lambda-field%lambdaj)) then
                    !
                    ! try to select the root#2 for the i-state first 
                    !
                      if (field%lambda/=0.and.(nint(field%sigmaj-field%sigmai))==( lambda_i(2)-field%lambdaj ) ) then
                          !
                          il_temp = lambda_i(2)
                          lambda_i(2) = lambda_i(1)
                          lambda_i(1) = il_temp
                          !
                          field%lambda = nint(lambda_i(1))
                          !
                          f_t = a
                          a(:,1) = f_t(:,2)
                          a(:,2) = f_t(:,1)
                          !
                      elseif( field%lambdaj/=0.and.(nint(field%sigmaj-field%sigmai))==( field%lambda-lambda_j(2) ) ) then
                          !
                          il_temp = lambda_j(2)
                          lambda_j(2) = lambda_j(1)
                          lambda_j(1) = il_temp
                          !
                          jroot = 2
                          field%lambdaj = nint(lambda_j(1))
                          !
                          f_t = b
                          b(:,1) = f_t(:,2)
                          b(:,2) = f_t(:,1)
                          !
                      elseif( field%lambdaj/=0.and.(nint(field%sigmaj-field%sigmai))==( lambda_i(2)-lambda_j(2) ) ) then
                          !
                          il_temp = lambda_i(2)
                          lambda_i(2) = lambda_i(1)
                          lambda_i(1) = il_temp
                          !
                          field%lambda = nint(lambda_i(1))
                          !
                          il_temp = lambda_j(2)
                          lambda_j(2) = lambda_j(1)
                          lambda_j(1) = il_temp
                          !
                          jroot = 2
                          field%lambdaj = nint(lambda_j(1))
                          !
                          f_t = a
                          a(:,1) = f_t(:,2)
                          a(:,2) = f_t(:,1)
                          !
                          f_t = b
                          b(:,1) = f_t(:,2)
                          b(:,2) = f_t(:,1)
                          !
                      else
                          !
                          write(out,"(/'molpro_duo: cannot find selecion rule for',2i8,3x,a)") &
                                field%iref,field%jref,trim(field%class)
                          write(out,"(' sigma = ',2f8.1,' lambda (i) = ',i4,' lamda (j) = ',i4)") &
                              field%sigmai, field%sigmaj, int(lambda_i(1)),int(lambda_j(1))
                          stop 'molpro_duo: cannot find the selecion rules'
                      endif
                    !
                  endif
                  ! 
                case ('L+','ABINITIO-LX')
                  !
                  !write(out,"('molpro_duo: this L+-part is not implemented')")
                  !stop 'molpro_duo: this L+-part is not implemented'
                  ! 
                case ('DIPOLE')
                  !
                  !write(out,"('molpro_duo: this Dipole-part is not implemented')")
                  !stop 'molpro_duo: this Dipole-part is not implemented'
                  !
                case default
                  !
                  write(out,"(/'molpro_duo: this part is not implemented:',a)") trim(field%class)
                  stop 'molpro_duo: this part is not implemented'
                  !
                end select
                !
                !omp parallel do private(i) schedule(guided)
                do i=1,ngrid
                  !
                  coupling = 0 
                  if (ix_lz_y==0.and.jx_lz_y==0) then
                    !
                    coupling(1,1) = field%gridvalue(i)*field%complex_f 
                    !
                  elseif(ix_lz_y/=0.and.jx_lz_y==0) then
                    !
                    select case(trim(field%class))
                      !
                    case('SPINORBIT','ABINITIO-SPINORBIT')
                      !
                      if (abs(nint(field%sigmaj-field%sigmai))/=abs(field%lambda-field%lambdaj)) then
                        !
                        write(out,"(/'molpro_duo: SO ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") &
                              field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                        stop 'molpro_duo: illegal selection rules '
                        !
                      endif
                      !
                      ! for SOX it is either <1.2|SOX (regular) or <1.3|SOY (iregular)
                      !
                      if (poten(field%jstate)%parity%pm==-1) then 
                        !
                        ! regular
                        coupling(1,1) = field%gridvalue(i)*field%complex_f
                        coupling(2,1) =-coupling(1,1)*conjg(a(1,2)/a(2,2))
                        !
                      else
                        !
                        ! irregular
                        coupling(2,1) = field%gridvalue(i)*field%complex_f
                        coupling(1,1) =-coupling(2,1)*conjg(a(2,2)/a(1,2))
                        !
                      endif
                      ! 
                    case ('L+','ABINITIO-LX')
                      !
                      if (poten(field%jstate)%parity%pm==-1) then 
                        !
                        ! regular
                        coupling(1,1) = field%gridvalue(i)*field%complex_f
                        coupling(2,1) =-coupling(1,1)*conjg(a(1,2)/a(2,2))
                        !
                      else
                        !
                        ! iregular
                        coupling(2,1) = field%gridvalue(i)*field%complex_f
                        coupling(1,1) = -coupling(2,1)*conjg(a(2,2)/a(1,2))
                        !
                      endif
                      ! 
                    case ('DIPOLE')
                      !
                      if (poten(field%jstate)%parity%pm==-1) then 
                        !
                        ! iregular
                        coupling(2,1) = field%gridvalue(i)*field%complex_f*(-sqrt(0.5_rk))
                        coupling(1,1) =-coupling(2,1)*conjg(a(2,2)/a(1,2))
                        !
                      else
                        !
                        ! regular
                        coupling(1,1) = field%gridvalue(i)*field%complex_f*(-sqrt(0.5_rk))
                        coupling(2,1) =-coupling(1,1)*conjg(a(1,2)/a(2,2))
                        !
                      endif
                      !
                    case default
                      !
                      write(out,"(/'molpro_duo (lambdaj=0): for class ',a,' is not implemented ')") field%class
                      stop 'molpro_duo: not for this object'
                      !
                    end select
                    !
                  elseif(ix_lz_y==0.and.jx_lz_y/=0) then
                    !
                    select case(trim(field%class)) 
                      !
                    case('SPINORBIT','ABINITIO-SPINORBIT')
                      !
                      if (abs(nint(field%sigmaj-field%sigmai))/=abs(field%lambda-field%lambdaj)) then
                        !
                        write(out,"(/'molpro_duo: SO ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") & 
                              field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                        stop 'molpro_duo: illegal selection rules '
                        !
                      endif
                      !
                      ! The following relations are found using the condition <0,Sigma|SOX+SOY|+/-1,Sigma+/-1> = 0
                      ! Regualr case is for <0|SOX|Pix> and iregular is for <0|SOX|Piy>
                      !
                      if (poten(field%istate)%parity%pm==-1) then 
                        !
                        ! regular
                        coupling(1,1) = field%gridvalue(i)*field%complex_f
                        coupling(1,2) =-coupling(1,1)*b(1,2)/b(2,2)
                        !
                      else
                        !
                        ! iregular
                        coupling(1,2) = field%gridvalue(i)*field%complex_f
                        coupling(1,1) =-coupling(1,2)*b(2,2)/b(1,2)
                        !
                      endif
                      ! 
                    case('L+','ABINITIO-LX')
                      !
                      ! The following relations are found using the condition <0|LX+iLY|+1> = 0
                      ! Regualr case is for <0|LX|Pix> and iregular is for <0|LX|Piy>
                      !
                      if (poten(field%istate)%parity%pm==-1) then 
                        !
                        ! regular
                        coupling(1,1) = field%gridvalue(i)*field%complex_f
                        coupling(1,2) =-coupling(1,1)*b(1,2)/b(2,2)
                        !
                      else
                        !
                        ! iregular
                        coupling(1,2) = field%gridvalue(i)*field%complex_f
                        coupling(1,1) =-coupling(1,2)*b(2,2)/b(1,2)
                        !
                      endif
                      ! 
                    case ('DIPOLE')
                      !
                      ! The following relations are found using the condition <0|-mux/sqrt(2)+imuy/sqrt(2)|+1> = 0
                      ! Regualr case is for <0|muX|Pix> and iregular is for <0|mux|Piy>
                      !
                      if (poten(field%istate)%parity%pm==-1) then 
                        !
                        coupling(1,2) = field%gridvalue(i)*field%complex_f*(-sqrt(0.5_rk))
                        coupling(1,1) = -coupling(1,2)*b(2,2)/b(1,2)
                        !
                      else
                        !
                        coupling(1,1) = field%gridvalue(i)*field%complex_f*(-sqrt(0.5_rk))
                        coupling(1,2) = -coupling(1,1)*b(1,2)/b(2,2)
                        !
                      endif
                      !
                    case default
                      !
                      write(out,"(/'molpro_duo (lambdai=0): for class ',a,' is not implemented ')") field%class
                      stop 'molpro_duo: not for this object'
                      !
                    end select
                    !
                  elseif (abs(ix_lz_y)/=abs(jx_lz_y)) then
                    !
                    select case(trim(field%class))
                      !
                    case('SPINORBIT','ABINITIO-SPINORBIT')
                      !
                      if (abs(nint(field%sigmaj-field%sigmai))/=abs(field%lambda-field%lambdaj)) then
                        !
                        write(out,"(/'molpro_duo: SO ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") &
                              field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                        stop 'molpro_duo: illegal selection rules '
                        !
                      endif
                      !
                      c = field%gridvalue(i)*field%complex_f
                      !
                      !coupling(1,1) =  0
                      !coupling(1,2) =  c
                      !coupling(2,1) =  -c*conjg(a(1,1))*b(2,2)/(conjg(a(2,1))*b(1,2))
                      !coupling(2,2) =  0
                      !
                      ! these cominations are found using the following three conditions 
                      ! <-Lamba,Sigma|SO|Lambda',Sigma'> = 0 
                      ! <Lamba,Sigma|SO|-Lambda',Sigma'> = 0 
                      ! <-Lamba,Sigma|SO|-Lambda',Sigma'> = 0 
                      ! where we assume 
                      ! c = <x|LSX|y> and  < Lamba,Sigma|SO|Lambda',Sigma'> /= 0 
                      !
                      coupling(1,1) =  -c*b(2,2)/b(1,2)
                      coupling(1,2) =   c
                      coupling(2,1) =   c*conjg(a(1,2)/a(2,2))*b(2,2)/b(1,2)
                      coupling(2,2) =  -c*conjg(a(1,2)/a(2,2))
                      !
                    case('DIPOLE')
                      !
                      if (abs(field%lambda-field%lambdaj)/=1) then
                        !
                        write(out,"('molpro_duo: DIPOLE ',2i4,'; illegal selection rules for lambda = ',2i4,' not +/-1')") & 
                              field%iref,field%jref,field%lambda,field%lambdaj
                        stop 'molpro_duo: illegal selection rules for transition dipole'
                        !
                      endif
                      !
                      c = -field%gridvalue(i)*field%complex_f*sqrt(0.5_rk)
                      !
                      ! maple for Lx
                      !c = -conjugate(A[1,2])*a/conjugate(A[2,2]), b = -a*B[1,2]/B[2,2], 
                      !d = B[1,2]*conjugate(A[1,2])*a/(B[2,2]*conjugate(A[2,2]))
                      ! m+ = -1/sqrt(2)l+
                      !
                      ! these cominations are found using the following three conditions 
                      ! <-Lamba|r+|Lambda'> = 0 
                      ! <Lamba|r+|-Lambda'> = 0 
                      ! <-Lamba|r+|-Lambda'> = 0 
                      ! where we assume 
                      ! c = <x|r+|x> and < Lamba|L+|Lambda'> /= 0 
                      !
                      coupling(1,1) =  c
                      coupling(1,2) = -c*b(1,2)/b(2,2)
                      coupling(2,1) = -c*conjg(a(1,2)/a(2,2))
                      coupling(2,2) =  c*conjg(a(1,2)/a(2,2))*b(1,2)/b(2,2)
                      !
                      !coupling(1,1) =  -c*b(2,2)/b(1,2)
                      !coupling(1,2) =   c
                      !coupling(2,1) =   c*conjg(a(1,2)/a(2,2))*b(2,2)/b(1,2)
                      !coupling(2,2) =  -c*conjg(a(1,2)/a(2,2))
                      !
                    case('L+','ABINITIO-LX')
                      !
                      if (abs(field%lambda-field%lambdaj)/=1) then
                        !
                        write(out,"('molpro_duo: L+ ',2i4,'; illegal selection rules for lambda = ',2i4,' not +/-1')") & 
                              field%iref,field%jref,field%lambda,field%lambdaj
                        stop 'molpro_duo: L+ illegal selection rules for transition dipole'
                        !
                      endif
                      !
                      ! The <x|Lx+iLy|y> element is <x|Lx|y>
                      !
                      c = field%gridvalue(i)*field%complex_f
                      !
                      ! maple:
                      !
                      !d = -b*conjugate(A[1,2])/conjugate(A[2,2]), c = B[2,2]*b*conjugate(A[1,2])/(B[1,2]*conjugate(A[2,2])), 
                      !a = -B[2,2]*b/B[1,2]
                      !
                      ! these cominations are found using the following three conditions 
                      ! <-Lamba|L+|Lambda'> = 0 
                      ! <Lamba|L+|-Lambda'> = 0 
                      ! <-Lamba|L+|-Lambda'> = 0 
                      ! where we assume 
                      ! c = <x|L+|y> and  < Lamba|L+|Lambda'> /= 0 
                      !
                      coupling(1,1) = -c*b(2,2)/b(1,2)
                      coupling(1,2) =  c
                      coupling(2,1) =  c*conjg(a(1,2)/a(2,2))*b(2,2)/b(1,2)
                      coupling(2,2) = -c*conjg(a(1,2)/a(2,2))
                      !
                    case default
                      !
                      write(out,"('molpro_duo (lambdaj<>lambdai): for class ',a,' is not implemented ')") field%class
                      stop 'molpro_duo: not for this object'
                      !
                    end select
                    !
                  else
                    !
                    select case(trim(field%class))
                      !
                    case('SPINORBIT','ABINITIO-SPINORBIT')
                      !
                      if (nint(field%sigmaj-field%sigmai)/=0) then
                        !
                        write(out,"('molpro_duo: SOZ ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") & 
                              field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                        stop 'molpro_duo: SOZ illegal selection rules '
                        !
                      endif
                      !
                      if ( field%sigmai<0 ) then
                        !
                        !write(out,"('molpro_duo: SO ',2i4,'; illegal reference sigmai <0 ',2f8.1)") &
                        !      field%iref,field%jref,field%sigmai,field%sigmaj
                        !stop 'molpro_duo: illegal reference sigmai'
                        !
                      endif
                      !
                      c = field%gridvalue(i)*field%complex_f
                      !
                      !coupling(1,1) =  c
                      !coupling(1,2) =  0
                      !coupling(2,1) =  0
                      !coupling(2,2) = -c*b(1,2)/b(2,2)*conjg(a(1,1))/conjg(a(2,1))
                      !
                      coupling(1,1) =  0
                      coupling(1,2) =  c
                      coupling(2,1) = -conjg(a(1,1))*c*b(2,2)/(conjg(a(2,1))*b(1,2))
                      coupling(2,2) =  0
                      !
                    case('DIPOLE')
                      !
                      if (abs(field%lambda-field%lambdaj)/=0) then
                        !
                        write(out,"('molpro_duo: DMZ ',2i4,'; illegal selection rules for lambda = ',2i4,' not 0')") &
                              field%iref,field%jref,field%lambda,field%lambdaj
                        stop 'molpro_duo: illegal selection rules for DMZ'
                        !
                      endif
                      !
                      c = field%gridvalue(i)*field%complex_f
                      !
                      !coupling(1,1) =  0
                      !coupling(1,2) =  c
                      !coupling(2,1) = -conjg(a(1,1))*c*b(2,2)/(conjg(a(2,1))*b(1,2))
                      !coupling(2,2) =  0
                      !
                      coupling(1,1) =  c
                      coupling(1,2) =  0
                      coupling(2,1) =  0
                      coupling(2,2) =  c
                      !
                    case default
                      !
                      write(out,"('molpro_duo: for class ',a,' is not implemented ')") field%class
                      stop 'molpro_duo: not for this object'
                      !
                    end select
                    !
                  endif
                  !
                  f_t = matmul( conjg(transpose(a)),matmul(coupling,(b)) )
                  !
                  field%gridvalue(i) = real(f_t(1,1))
                  !
                  if (any( abs( aimag( f_t ) )>small_ ) ) then
                    !
                    write(out,"(/'molpro_duo: ',a,' ',2i3,'; duo-complex values ',8f8.1)") trim(field%class),field%iref,&
                                                                                            field%jref,f_t(:,:)
                                                                                            
                    write(out,"('Please check the MOLPRO output and use FACTOR I if required')")
                    write(out,"('It is important that MOLPROS LSX (spin-orbit-x), LX and DMX (dipole-x) values are used')")
                    write(out, &
                    "('If this is input from older publications (before 2019) please try keyword LEGACY anywhere in input')")
                    write(out,"('This will use the previously standard for the molpro objects')")
                    stop 'molpro_duo error: duo-complex values'
                    !
                  endif
                  !
                  if (abs( real( f_t(1,1) ) )<=sqrt(small_) .and.abs(field%gridvalue(i))>sqrt(small_)) then
                    !
                    write(out,"('molpro_duo: ',a,' ',2i3,'; duo-zero values ',8f8.1)") &
                          trim(field%class),field%iref,field%jref,f_t(:,:)
                    stop 'molpro_duo: duo-zero values ?'
                    !
                  endif
                  !
                enddo
                !omp end parallel do
   
        end subroutine molpro_duo
   
   
        subroutine molpro_duo_old_2018(field)
           !
           use lapack,only : lapack_zheev     
           !
           type(fieldT),intent(inout) :: field
           integer(ik) :: ix_lz_y,jx_lz_y,iroot,jroot,il_temp,ngrid
           complex(rk) :: a(2,2),b(2,2),coupling(2,2),f_t(2,2),b0(2,2),c
           real(rk)    :: lambda_i(2),lambda_j(2)
               !
               ngrid = grid%npoints
               !
               ix_lz_y = field%ix_lz_y
               jx_lz_y = field%jx_lz_y
               !
               a = 0 ; b = 0
               !
               a(1,1) =  cmplx(1.0_rk,0.0_rk, kind=rk)
               b(1,1) =  cmplx(1.0_rk,0.0_rk, kind=rk)
               a(2,2) =  cmplx(1.0_rk,0.0_rk, kind=rk)
               b(2,2) =  cmplx(1.0_rk,0.0_rk, kind=rk)
               !
               iroot = 1
               jroot = 1
               !
               if (ix_lz_y/=0) then
                 a = 0 
                 a(1,2) = cmplx(0.0_rk,ix_lz_y , kind=rk)
                 a(2,1) = cmplx(0.0_rk,-ix_lz_y, kind=rk)
                 !
                 call lapack_zheev(a,lambda_i)
                 !
                 ! swap to have the first root positive 
                 !
                 f_t = a
                 a(:,1) = f_t(:,2)
                 a(:,2) = f_t(:,1)
                 !
                 il_temp = lambda_i(2)
                 lambda_i(2) = lambda_i(1)
                 lambda_i(1) = il_temp
                 !
                 field%lambda = nint(lambda_i(1))
                 !
                 a = a*cmplx(0.0_rk,1.0_rk, kind=rk)
                 !
               elseif (poten(field%istate)%parity%pm==-1) then 
                 !
                 a(1,1) =  cmplx(0.0_rk,1.0_rk, kind=rk)
                 !
               endif
               !
               if (jx_lz_y/=0) then
                 b = 0 
                 b(1,2) = cmplx(0.0_rk,jx_lz_y , kind=rk)
                 b(2,1) = cmplx(0.0_rk,-jx_lz_y, kind=rk)
                 !
                 b0 = b
                 !
                 call lapack_zheev(b,lambda_j)
                 !
                 ! swap to have the first root positive 
                 !
                 f_t = b
                 b(:,1) = f_t(:,2)
                 b(:,2) = f_t(:,1)
                 !
                 il_temp = lambda_j(2)
                 lambda_j(2) = lambda_j(1)
                 lambda_j(1) = il_temp
                 !
                 field%lambdaj = nint(lambda_j(1))
                 !
                 b = b*cmplx(0.0_rk,1.0_rk, kind=rk)
                 !
                 f_t = matmul( conjg(transpose(b)),matmul(b0,(b)) )
                 !
               elseif (poten(field%jstate)%parity%pm==-1) then 
                 !
                 b(1,1) =  cmplx(0.0_rk,1.0_rk, kind=rk)
                 !
               endif
               !
               ! Check the selection rules
               select case(trim(field%class))
                 !
               case('SPINORBIT','ABINITIO-SPINORBIT')
                 !
                 if ((nint(field%sigmaj-field%sigmai))/=(field%lambda-field%lambdaj)) then
                   !
                   ! try to select the root#2 for the i-state first 
                   !
                   if (field%lambda/=0.and.(nint(field%sigmaj-field%sigmai))==( lambda_i(2)-field%lambdaj ) ) then
                     !
                     il_temp = lambda_i(2)
                     lambda_i(2) = lambda_i(1)
                     lambda_i(1) = il_temp
                     !
                     field%lambda = nint(lambda_i(1))
                     !
                     f_t = a
                     a(:,1) = f_t(:,2)
                     a(:,2) = f_t(:,1)
                     !
                   elseif( field%lambdaj/=0.and.(nint(field%sigmaj-field%sigmai))==( field%lambda-lambda_j(2) ) ) then
                     !
                     il_temp = lambda_j(2)
                     lambda_j(2) = lambda_j(1)
                     lambda_j(1) = il_temp
                     !
                     jroot = 2
                     field%lambdaj = nint(lambda_j(1))
                     !
                     f_t = b
                     b(:,1) = f_t(:,2)
                     b(:,2) = f_t(:,1)
                     !
                   elseif( field%lambdaj/=0.and.(nint(field%sigmaj-field%sigmai))==( lambda_i(2)-lambda_j(2) ) ) then
                     !
                     il_temp = lambda_i(2)
                     lambda_i(2) = lambda_i(1)
                     lambda_i(1) = il_temp
                     !
                     field%lambda = nint(lambda_i(1))
                     !
                     il_temp = lambda_j(2)
                     lambda_j(2) = lambda_j(1)
                     lambda_j(1) = il_temp
                     !
                     jroot = 2
                     field%lambdaj = nint(lambda_j(1))
                     !
                     f_t = a
                     a(:,1) = f_t(:,2)
                     a(:,2) = f_t(:,1)
                     !
                     f_t = b
                     b(:,1) = f_t(:,2)
                     b(:,2) = f_t(:,1)
                     !
                   else
                     !
                     write(out,"(/'molpro_duo: cannot find selecion rule for',2i8,3x,a)") field%iref,field%jref,trim(field%class)
                     write(out,"(' sigma = ',2f8.1,' lambda (i) = ',i4,' lamda (j) = ',i4)") &
                                field%sigmai, field%sigmaj, int(lambda_i(1)),int(lambda_j(1))
                     stop 'molpro_duo: cannot find the selecion rules'
                   endif
                   !
                 endif
                 ! 
               case ('L+','ABINITIO-LX')
                 !
                 !write(out,"('molpro_duo: this L+-part is not implemented')")
                 !stop 'molpro_duo: this L+-part is not implemented'
                 ! 
               case ('DIPOLE')
                 !
                 !write(out,"('molpro_duo: this Dipole-part is not implemented')")
                 !stop 'molpro_duo: this Dipole-part is not implemented'
                 !
               case default
                 !
                 write(out,"(/'molpro_duo: this part is not implemented:',a)") trim(field%class)
                 stop 'molpro_duo: this part is not implemented'
                 !
               end select
               !
               !omp parallel do private(i) schedule(guided)
               do i=1,ngrid
                 !
                 coupling = 0 
                 if (ix_lz_y==0.and.jx_lz_y==0) then
                   !
                   coupling(1,1) = field%gridvalue(i)*field%complex_f 
                   !
                 elseif(ix_lz_y/=0.and.jx_lz_y==0) then
                   !
                   !if ( lambda_i(1)/=1_rk.or.lambda_i(2)/=-1.0_rk) then
                   !  !
                   !  write(out,"('molpro_duo: lambda_i ',2f8.1,' are not -1 and 1, coupling')") lambda_i,field%iref,field%jref
                   !  !stop 'molpro_duo: not for this object'
                   !  !
                   !endif
                   !
                   select case(trim(field%class))
                     !
                   case('SPINORBIT','ABINITIO-SPINORBIT')
                     !
                     if (abs(nint(field%sigmaj-field%sigmai))/=abs(field%lambda-field%lambdaj)) then
                       !
                       write(out,"(/'molpro_duo: SO ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") &
                            field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                       stop 'molpro_duo: illegal selection rules '
                       !
                     endif
                     !
                     if ( field%sigmai<0 ) then
                       !
                       !write(out,"('molpro_duo: SO ',2i4,'; illegal reference sigmai <0 ?',2f8.1)") &
                       !    field%iref,field%jref,field%sigmai,field%sigmaj
                       !stop 'molpro_duo: illegal reference sigmai'
                       !
                     endif
                     !
                     ! for SOX it is always <1.3| which is given, i.e. we need to solve for the 1st, <1.2| component:
                     !
                     if (field%lambda>0) then 
                       !
                       coupling(2,1) = field%gridvalue(i)*field%complex_f
                       coupling(1,1) =-field%gridvalue(i)*field%complex_f*conjg(a(2,2))/conjg(a(1,2))
                       !
                     else 
                       !
                       coupling(2,1) = field%gridvalue(i)*field%complex_f
                       coupling(1,1) =-field%gridvalue(i)*field%complex_f*conjg(a(2,1))/conjg(a(1,1))
                       !
                     endif 
                     ! 
                   case ('L+','ABINITIO-LX')
                     !
                     ! eigen-vector 2 is for Lambda
                     !
                     coupling(1,1) = field%gridvalue(i)*field%complex_f*cmplx(0.0_rk,-1.0_rk,kind=rk)  
                     coupling(2,1) = field%gridvalue(i)*field%complex_f*cmplx(0.0_rk, 1.0_rk,kind=rk)*conjg(a(1,2))/conjg(a(2,2))
                     ! 
                   case ('DIPOLE')
                     !
                     ! eigen-vector 1 is for -Lambda
                     !
                     coupling(1,1) = field%gridvalue(i)*field%complex_f*(-sqrt(0.5_rk))
                     coupling(2,1) = field%gridvalue(i)*field%complex_f*(conjg(a(1,2)/a(2,2))*sqrt(0.5_rk))
                     !
                   case default
                     !
                     write(out,"(/'molpro_duo (lambdaj=0): for class ',a,' is not implemented ')") field%class
                     stop 'molpro_duo: not for this object'
                     !
                   end select
                   !
                 elseif(ix_lz_y==0.and.jx_lz_y/=0) then
                   !
                   !if (lambda_j(1)/=1_rk.or.lambda_j(2)/=-1.0_rk) then
                   !  !
                   !  !write(out,"('molpro_duo: lambda_j ',2f8.1,' are not -1 and 1, coupling for states',2i)") lambda_j,
                   !                                                                                field%iref,field%jref
                   !  !stop 'molpro_duo: not for this object'
                   !  !
                   !endif
                   !
                   select case(trim(field%class)) 
                     !
                   case('SPINORBIT','ABINITIO-SPINORBIT')
                     !
                     if (abs(nint(field%sigmaj-field%sigmai))/=abs(field%lambda-field%lambdaj)) then
                       !
                       write(out,"(/'molpro_duo: SO ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") & 
                             field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                       stop 'molpro_duo: illegal selection rules '
                       !
                     endif
                     !
                     !if ( field%sigmaj<0 ) then
                     !  !
                     !  write(out,"('molpro_duo: SO ',2i4,'; illegal reference sigmaj <0 ',2f8.1)") & 
                     !        field%iref,field%jref,field%sigmai,field%sigmaj
                     !  stop 'molpro_duo: illegal reference sigmaj'
                     !  !
                     !endif
                     !
                     ! eigen-vector 1 is for Lambda
                     !
                     ! for SOX the non-zero is for |y> vector, i.e. the second component of coupling
                     !
                     !coupling(1,1) = -field%gridvalue(i)*field%complex_f*b(1,2)/b(2,2)
                     !coupling(1,2) = field%gridvalue(i)*field%complex_f
   
                     ! for SOX it is always <1.3| which is given, i.e. we need to solve for the 1st, <1.2| component:
                     !
                     if (field%lambdaj>0) then 
                       !
                       coupling(1,2) = field%gridvalue(i)*field%complex_f
                       coupling(1,1) =-field%gridvalue(i)*field%complex_f*b(2,1)/b(1,1)
                       !
                     else 
                       !
                       coupling(1,2) = field%gridvalue(i)*field%complex_f
                       coupling(1,1) =-field%gridvalue(i)*field%complex_f*b(2,2)/b(1,2)
                       !
                     endif 
                     ! 
                   case('L+','ABINITIO-LX')
                     !
                     ! eigen-vector 1 is for Lambda
                     !
                     coupling(1,1) = field%gridvalue(i)*field%complex_f*cmplx(0.0_rk, 1.0_rk,kind=rk)  
                     coupling(1,2) = field%gridvalue(i)*field%complex_f*cmplx(0.0_rk,-1.0_rk,kind=rk)*b(1,2)/b(2,2)
                     ! 
                   case ('DIPOLE')
                     !
                     ! eigen-vector 1 is for Lambda
                     !
                     coupling(1,1) = field%gridvalue(i)*field%complex_f*(-sqrt(0.5_rk))
                     coupling(1,2) = field%gridvalue(i)*field%complex_f*(b(1,2)/b(2,2)*sqrt(0.5_rk))
                     !
                   case default
                     !
                     write(out,"(/'molpro_duo (lambdai=0): for class ',a,' is not implemented ')") field%class
                     stop 'molpro_duo: not for this object'
                     !
                   end select
                   !
                 elseif (abs(ix_lz_y)/=abs(jx_lz_y)) then
                   !
                   select case(trim(field%class))
                     !
                   case('SPINORBIT','ABINITIO-SPINORBIT')
                     !
                     if (abs(nint(field%sigmaj-field%sigmai))/=abs(field%lambda-field%lambdaj)) then
                       !
                       write(out,"(/'molpro_duo: SO ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") &
                             field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                       stop 'molpro_duo: illegal selection rules '
                       !
                     endif
                     !
                     if ( field%sigmai<0 ) then
                       !
                       !write(out,"('molpro_duo: SO ',2i4,'; illegal reference sigmai <0 ?',2f8.1)") &
                       !      field%iref,field%jref,field%sigmai,field%sigmaj
                       !stop 'molpro_duo: illegal reference sigmai'
                       !
                     endif
                     !
                     c = field%gridvalue(i)*field%complex_f
                     !
                     coupling(1,1) =  0
                     coupling(1,2) =  c
                     coupling(2,1) =  -c*conjg(a(1,1))*b(2,2)/(conjg(a(2,1))*b(1,2))
                     coupling(2,2) =  0
                     !
                   case('DIPOLE')
                     !
                     if (abs(field%lambda-field%lambdaj)/=1) then
                       !
                       write(out,"('molpro_duo: DIPOLE ',2i4,'; illegal selection rules for lambda = ',2i4,' not +/-1')") & 
                             field%iref,field%jref,field%lambda,field%lambdaj
                       stop 'molpro_duo: illegal selection rules for transition dipole'
                       !
                     endif
                     !
                     c = -field%gridvalue(i)*field%complex_f*sqrt(0.5_rk)
                     !
                     ! maple for Lx
                     !c = -conjugate(A[1,2])*a/conjugate(A[2,2]), b = -a*B[1,2]/B[2,2], 
                     !d = B[1,2]*conjugate(A[1,2])*a/(B[2,2]*conjugate(A[2,2]))
                     ! m+ = -1/sqrt(2)l+
                     !
                     coupling(1,1) =  c
                     coupling(1,2) = -c*b(1,2)/b(2,2)
                     coupling(2,1) = -c*conjg(a(1,2))/conjg(a(2,2))
                     coupling(2,2) =  c*b(1,2)*conjg(a(1,2))/(conjg(a(2,2))*b(2,2))
                     !
                   case('L+','ABINITIO-LX')
                     !
                     if (abs(field%lambda-field%lambdaj)/=1) then
                       !
                       write(out,"('molpro_duo: L+ ',2i4,'; illegal selection rules for lambda = ',2i4,' not +/-1')") & 
                             field%iref,field%jref,field%lambda,field%lambdaj
                       stop 'molpro_duo: L+ illegal selection rules for transition dipole'
                       !
                     endif
                     !
                     ! The <x|Lx+iLy|y> element is <x|Lx|y>
                     !
                     c = field%gridvalue(i)*field%complex_f
                     !
                     ! maple:
                     !
                     !d = -b*conjugate(A[1,2])/conjugate(A[2,2]), c = B[2,2]*b*conjugate(A[1,2])/(B[1,2]*conjugate(A[2,2])), 
                     !a = -B[2,2]*b/B[1,2]
                     !
                     coupling(1,2) =  c
                     coupling(1,1) = -c*b(2,2)/b(1,2)
                     coupling(2,1) =  c*b(2,2)*conjg(a(1,2))/(conjg(a(2,2))*b(1,2))
                     coupling(2,2) = -c*conjg(a(1,2))/conjg(a(2,2))
                     !
                   case default
                     !
                     write(out,"('molpro_duo (lambdaj<>lambdai): for class ',a,' is not implemented ')") field%class
                     stop 'molpro_duo: not for this object'
                     !
                   end select
                   !
                 else
                   !
                   select case(trim(field%class))
                     !
                   case('SPINORBIT','ABINITIO-SPINORBIT')
                     !
                     if (nint(field%sigmaj-field%sigmai)/=0) then
                       !
                       write(out,"('molpro_duo: SOZ ',2i4,'; illegal selection rules, sigma = ',2f8.1,' lambda = ',2i4)") & 
                             field%iref,field%jref,field%sigmai,field%sigmaj,field%lambda,field%lambdaj
                       stop 'molpro_duo: SOZ illegal selection rules '
                       !
                     endif
                     !
                     if ( field%sigmai<0 ) then
                       !
                       !write(out,"('molpro_duo: SO ',2i4,'; illegal reference sigmai <0 ',2f8.1)") &
                       !      field%iref,field%jref,field%sigmai,field%sigmaj
                       !stop 'molpro_duo: illegal reference sigmai'
                       !
                     endif
                     !
                     c = field%gridvalue(i)*field%complex_f
                     !
                     !coupling(1,1) =  c
                     !coupling(1,2) =  0
                     !coupling(2,1) =  0
                     !coupling(2,2) = -c*b(1,2)/b(2,2)*conjg(a(1,1))/conjg(a(2,1))
                     !
                     coupling(1,1) =  0
                     coupling(1,2) =  c
                     coupling(2,1) = -conjg(a(1,1))*c*b(2,2)/(conjg(a(2,1))*b(1,2))
                     coupling(2,2) =  0
                     !
                   case('DIPOLE')
                     !
                     if (abs(field%lambda-field%lambdaj)/=0) then
                       !
                       write(out,"('molpro_duo: DMZ ',2i4,'; illegal selection rules for lambda = ',2i4,' not 0')") &
                             field%iref,field%jref,field%lambda,field%lambdaj
                       stop 'molpro_duo: illegal selection rules for DMZ'
                       !
                     endif
                     !
                     c = field%gridvalue(i)*field%complex_f
                     !
                     !coupling(1,1) =  0
                     !coupling(1,2) =  c
                     !coupling(2,1) = -conjg(a(1,1))*c*b(2,2)/(conjg(a(2,1))*b(1,2))
                     !coupling(2,2) =  0
                     !
                     coupling(1,1) =  c
                     coupling(1,2) =  0
                     coupling(2,1) =  0
                     coupling(2,2) =  c
                     !
                   case default
                     !
                     write(out,"('molpro_duo: for class ',a,' is not implemented ')") field%class
                     stop 'molpro_duo: not for this object'
                     !
                   end select
                   !
                 endif
                 !
                 f_t = matmul( conjg(transpose(a)),matmul(coupling,(b)) )
                 !
                 field%gridvalue(i) = real(f_t(1,1))
                 !
                 if (any( abs( aimag( f_t ) )>small_ ) ) then
                   !
                   write(out,"(/'molpro_duo: ',a,' ',2i3,'; duo-complex values ',8f8.1)") trim(field%class),field%iref,&
                                                                                          field%jref,f_t(:,:)
                   stop 'molpro_duo: duo-complex values ?'
                   !
                 endif
                 !
                 if (abs( real( f_t(1,1) ) )<=sqrt(small_) .and.abs(field%gridvalue(i))>sqrt(small_)) then
                   !
                   write(out,"('molpro_duo: ',a,' ',2i3,'; duo-zero values ',8f8.1)") &
                          trim(field%class),field%iref,field%jref,f_t(:,:)
                   stop 'molpro_duo: duo-zero values ?'
                   !
                 endif
                 !
               enddo
               !omp end parallel do
   
        end subroutine molpro_duo_old_2018
   
   
        !
        subroutine check_and_print_coupling(N,iverbose,fl,name)
          !
          type(fieldT),intent(in) :: fl(:)
          integer(ik),intent(in)  :: N,iverbose
          character(len=*),intent(in) :: name
          integer(ik)             :: i,istate,ngrid
          character(len=200) :: filename
          integer :: u1
            !
            if (N<=0.or.iverbose<5) return
            !
            write(out,'(/a)') trim(name)
            !
            ngrid = grid%npoints
            !
            ! double check
            !
            do i=1,N
              !
              istate = fl(i)%istate
              jstate = fl(i)%jstate
              !
              if (trim(fl(i)%class)=='SPINORBIT'.and.(abs(fl(i)%sigmai)>fl(i)%spini.or.abs(fl(i)%sigmaj)>fl(i)%spinj)) then
                 write(out,'("For N =",i4," one of sigmas (",2f8.1,") large than  spins or undefined (",2f8.1,")")') & 
                           i,fl(i)%sigmai,fl(i)%sigmaj,fl(i)%spini,fl(i)%spinj
                 stop 'illegal sigma or spin'
              endif
              if (nint(2.0_rk*fl(i)%spini)+1/=poten(istate)%multi.or. &
                  nint(2.0_rk*fl(i)%spinj)+1/=poten(jstate)%multi ) then
                 write(out,'("For N =",i3," multi (",2i3,") dont agree with either of multi (",2i3,") of states ",i2," and ",i2)') &
                           i,nint(2.0*fl(i)%spini)+1,nint(2.0*fl(i)%spinj)+1,poten(istate)%multi, &
                           poten(jstate)%multi,istate,jstate
                 stop 'illegal multi in map_fields_onto_grid'
              endif
              !
              if (  fl(i)%lambda<bad_value+1.or.fl(i)%lambdaj<bad_value+1 ) then
                 write(out,'("For N =",i3," lambdas  are undefined of states ",i2," and ",i2," ",a)') &
                           i,istate,jstate,trim(name)
                 stop 'lambdas are undefined: map_fields_onto_grid'
              endif
              !
              if (  abs(fl(i)%lambda)/=poten(istate)%lambda.or.abs(fl(i)%lambdaj)/=poten(jstate)%lambda ) then
                 write(out,'("For N =",i3," lambdas (",2i3,") dont agree with either of lambdas (",2i3,") ' // &
                                                                        'of states ",i2," and ",i2)') &
                           i,fl(i)%lambda,fl(i)%lambdaj,poten(istate)%lambda, &
                           poten(jstate)%lambda,istate,jstate
                 stop 'illegal lambdas in map_fields_onto_grid'
              endif
              !
              if (  trim(name)=="Spin-Orbit:".and.(abs(fl(i)%sigmai)>fl(i)%spini.or.abs(fl(i)%sigmaj)>fl(i)%spinj) ) then
                 write(out,&
                 '("For N =",i3," sigmas (",2f9.2,") dont agree with their spins (",2f9.2,") of states ",i2," and ",i2)') &
                           i,fl(i)%sigmai,fl(i)%sigmai,fl(i)%spini, &
                           fl(i)%spinj,istate,jstate
                 stop 'illegal sigmas in map_fields_onto_grid'
              endif
              !
              if ( trim(name)=="Spin-Orbit:".and.(fl(i)%sigmai<bad_value+1.0.or.fl(i)%sigmaj<bad_value+1) ) then
                 write(out,'("For N =",i3," sigmas are undefined for states ",i2," and ",i2," ",a)') &
                           i,istate,jstate,trim(name)
                 stop 'sigmas are undefined: map_fields_onto_grid'
              endif
              !
              if ( trim(name)=="Spin-Orbit:".and.(fl(i)%sigmai<bad_value+1.0.or.fl(i)%sigmaj<bad_value+1.0) ) then
                 write(out,'("For N =",i3," sigmas are undefined for states ",i2," and ",i2," ",a)') &
                           i,istate,jstate,trim(name)
                 stop 'sigmas are undefined: map_fields_onto_grid'
              endif
              !
              if ( trim(name)=="<L+> functions:".and.istate==jstate ) then
                 write(out,'("For N =",i3," Lx/L+ are defined for the same state ",i2," and ",i2," ",a)') &
                           i,istate,jstate,trim(name)
                 stop 'illegal  - diagonal  - L+/Lx coupling: map_fields_onto_grid'
              endif
              !
              if ( trim(name(1:6))=="Lambda".and.istate/=jstate ) then
                 write(out,'("For N =",i3," Lambda-doubling must be defined for the same state, not ",i2," and ",i2," ",a)') &
                           i,istate,jstate,trim(name)
                 stop 'illegal  - non-diagonal  - Lambda doubling coupling: map_fields_onto_grid'
              endif
              !
            enddo
            !
            do istate=1,N
              write(out,'(i4,2x,a)') istate,trim(fl(istate)%name)
            enddo
            !
            !
            ! write to file if required
            if(job%print_pecs_and_couplings_to_file .eqv. .true.) then
              ! set up filename
              write(filename, '(a)') trim(name)
              i=len_trim(filename)
              if( filename(i:i) == ':') filename(i:i)=' ' !remove trailing colons
              do i=1, len_trim(filename) ! remove/change some character in the file name
                if( filename(i:i) == ' ') filename(i:i)='_' !spaces to underscores
                if( filename(i:i) == '*') filename(i:i)='_' !stars to underscores
                if( filename(i:i) == '<') filename(i:i)='_' !< to underscores
                if( filename(i:i) == '>') filename(i:i)='_' !> to underscores
                if( filename(i:i) == '(') filename(i:i)='_' !( to underscores
                if( filename(i:i) == ')') filename(i:i)='_' !) to underscores
                if( filename(i:i) == '+') filename(i:i)='p' !+ to `p'
              enddo
              filename=trim(filename) // '.dat'
   
              call IOstart("check and_print coupling",u1)
              open(unit=u1, file=trim(filename), status='unknown',action='write')
              write(my_fmt,'(A,I0,A)') '(f18.10,', N, '(f18.9))'
              do i=1,ngrid
               write(u1 ,my_fmt) r(i),(fl(istate)%gridvalue(i),istate=1,N)
              enddo
              close(u1)
              call IOstop("check and_print coupling")
            endif
            !
            write(my_fmt,'(A,I0,A)') '("        r(Ang)  ",2x,', N, '(i9,10x))'
            write(out,my_fmt) (istate,istate=1,N)
            !
            write(my_fmt,'(A,I0,A)') '(f18.8,', N, '(1x,f18.8))'
            do i=1,ngrid
               write(out,my_fmt) r(i),(fl(istate)%gridvalue(i),istate=1,N)
            enddo
            !
        end subroutine check_and_print_coupling
        !
   
        !
        subroutine check_and_print_field(N,iverbose,fl,name)
          !
          type(fieldT),intent(in) :: fl(:)
          integer(ik),intent(in)  :: N,iverbose
          character(len=*),intent(in) :: name
          integer(ik)             :: i,istate,ngrid
          character(len=100)   :: filename ! file containing tabulated field (for plotting etc)
          integer, parameter  :: u1=1001
              !
              if (N<=0.or.iverbose<5) return
              !
              write(out,'(/a)') trim(name)
              !
              ngrid = grid%npoints
              !
              do istate=1,N
                 !
                 write(out,'(i6,2x,a30,2x,a7,a20)') istate,trim(fl(istate)%name), "type = " , trim(fl(istate)%type)
                 !
                 ! double check
                 if (fl(istate)%lambda == 0 .and.fl(istate)%parity%pm==0 ) then
                    write(out,'("Please define the +/- symmetry of the Lambda=0 state (",i4,")")') istate
                    write(out,'("It is important for the SO component")')
                    stop 'Parity (+/-) for Lambda=0 undefined'
                 endif
                 !
              enddo
              ! 
              ! write to file if required
              if(job%print_pecs_and_couplings_to_file .eqv. .true.) then
                ! set up filename for output
                write(filename, '(a)') trim(name)
                i=len_trim(filename)
                if( filename(i:i) == ':') filename(i:i)=' ' !remove trailing colons
                do i=1, len_trim(filename) ! remove/change some character in the file name
                  if( filename(i:i) == ' ') filename(i:i)='_' !spaces to underscores
                  if( filename(i:i) == '*') filename(i:i)='_' !stars to underscores
                  if( filename(i:i) == '<') filename(i:i)='_' !< to underscores
                  if( filename(i:i) == '>') filename(i:i)='_' !> to underscores
                  if( filename(i:i) == '(') filename(i:i)='_' !( to underscores
                  if( filename(i:i) == ')') filename(i:i)='_' !) to underscores
                  if( filename(i:i) == '+') filename(i:i)='p' !+ to `p'
                enddo
                filename=trim(filename) // '.dat'
                open(unit=u1, file=trim(filename), status='unknown',action='write')
                write(my_fmt, '(A,I0,A)') '(f18.10,', N, '(es22.14))'
                do i=1,ngrid
                  write(u1, my_fmt) r(i),(fl(istate)%gridvalue(i),istate=1,N)
                enddo
                close(u1)
              endif
              !
              write(my_fmt, '(A,I0,A)') '("            r(Ang)",', N, '(i22))'
              write(out,my_fmt) (istate,istate=1,N)
              !
              write(my_fmt, '(A,I0,A)') '(f18.8,', N, '(es22.8))'
              do i=1,ngrid
                 write(out,my_fmt) r(i),(fl(istate)%gridvalue(i),istate=1,N)
              enddo
              !
              ! write equilibrium properties
              write(out,*)
              write(out,'(A)') 'Equilibrium properties'
              !
              write(out, '(A18)',advance='no') 'imin = '
              do istate=1,N; write(out,'(I22)',advance='no') fl(istate)%imin ; enddo ! grid index for which poten is minimum
              write(out,*)
              !
              write(out, '(A18)',advance='no') 'r_imin / ang = '
              do istate=1,N; write(out,'(F22.10)',advance='no') fl(istate)%rimin ; enddo ! minimum of potential on the grid
              write(out,*)
              !
              write(out, '(A18)',advance='no') 'V(r_imin) / cm-1 = '
              do istate=1,N; write(out,'(F22.8)',advance='no') fl(istate)%Vimin ; enddo ! minimum of potential on the grid
              write(out,*)
              !
              write(out, '(A18)',advance='no') 'Has a single min?'
              do istate=1,N;
                 if(fl(istate)%zHasMinimum) then
                      write(out,'(A22)',advance='no') 'Yes'
                 else
                      write(out,'(A22)',advance='no') ' No'
                 endif
              enddo ! minimum of potential on the grid
              write(out,*)
              !
              write(out, '(A18)',advance='no') 'True r_e  / ang ='
              do istate=1,N;
                 if(fl(istate)%zHasMinimum) then
                   write(out,'(F22.10)',advance='no') fl(istate)%re
                 else
                   write(out, '(a22)', advance='no') 'N.A.'
                 endif
              enddo
              write(out,*)
              !
              write(out,'(/,A)') 'Derivatives at true re'
              !
              write(out, '(A18)',advance='no') 'der0, cm-1       ='
              do istate=1,N
                 if(fl(istate)%zHasMinimum) then
                   write(out,'(F22.10)',advance='no') fl(istate)%V0
                 else
                   write(out, '(a22)', advance='no') 'N.A.'
                 endif
              enddo
              write(out,*)
              !
              write(out, '(A18)',advance='no') "der1, cm-1/ang   ="
              do istate=1,N
                 if(fl(istate)%zHasMinimum) then
                   write(out,'(F22.10)',advance='no') fl(istate)%V1
                 else
                   write(out, '(a22)', advance='no') 'N.A.'
                 endif
              enddo
              write(out,*)
              !
              write(out, '(A18)',advance='no') "der2, cm-1/ang^2 ="
              do istate=1,N
                 if(fl(istate)%zHasMinimum) then
                   write(out,'(F22.10)',advance='no') fl(istate)%V2
                 else
                   write(out, '(a22)', advance='no') 'N.A.'
                 endif
              enddo
              write(out,*)
              !
              write(out, '(A18)',advance='no') "der3, cm-1/ang^3 ="
              do istate=1,N
                 if(fl(istate)%zHasMinimum) then
                   write(out,'(F22.10)',advance='no') fl(istate)%V3
                 else
                   write(out, '(a22)', advance='no') 'N.A.'
                 endif
              enddo
              write(out,*)
              !
              write(out, '(A18)',advance='no') "der4, cm-1/ang^4 ="
              do istate=1,N
                 if(fl(istate)%zHasMinimum) then
                   write(out,'(F22.10)',advance='no') fl(istate)%V4
                 else
                   write(out, '(a22)', advance='no') 'N.A.'
                 endif
              enddo
              write(out,*)
              !
              write(out,*)
              !
              write(out, '(A18)',advance='no') "Harmonic we, cm-1="
              do istate=1,N
                 if(fl(istate)%zHasMinimum) then
                   write(out,'(F22.10)',advance='no') fl(istate)%we
                 else
                   write(out, '(a22)', advance='no') 'N.A.'
                 endif
              enddo
              write(out,*)
              !
              write(out, '(A18)',advance='no') "Rotat. B0, cm-1="
              do istate=1,N
                 if(fl(istate)%zHasMinimum) then
                   write(out,'(F22.10)',advance='no') fl(istate)%B0
                 else
                   write(out, '(a22)', advance='no') 'N.A.'
                 endif
              enddo
              write(out,*)
              !
              write(out, '(A18)',advance='no') "Anharm. const. xe="
              do istate=1,N
                 if(fl(istate)%zHasMinimum) then
                   write(out,'(F22.10)',advance='no') fl(istate)%xe
                 else
                   write(out, '(a22)', advance='no') 'N.A.'
                 endif
              enddo
              write(out,*)
              !
              write(out, '(A18)',advance='no') "Coriol. ae, cm-1="
              do istate=1,N
                 if(fl(istate)%zHasMinimum) then
                   write(out,'(F22.10)',advance='no') fl(istate)%alphae
                 else
                   write(out, '(a22)', advance='no') 'N.A.'
                 endif
              enddo
              write(out,*)
              !
              write(out, '(A18)',advance='no') "Centr. De, cm-1="
              do istate=1,N
                 if(fl(istate)%zHasMinimum) then
                   write(out,'(F22.10)',advance='no') fl(istate)%Debar
                 else
                   write(out, '(a22)', advance='no') 'N.A.'
                 endif
              enddo
              write(out,*)
              !
              write(out, '(A18)',advance='no') "Y00, cm-1="
              do istate=1,N
                 if(fl(istate)%zHasMinimum) then
                   write(out,'(F22.10)',advance='no') fl(istate)%Y00
                 else
                   write(out, '(a22)', advance='no') 'N.A.'
                 endif
              enddo
              write(out,*)
              !
              write(out,*)
              write(out,'(A)') 'Approximate J=0 vibrational energy levels (no couplings) '
              write(out,'(A)') ' given by E(v, J=0) = V(re) + Y00 + we*(v+0.5) - we*xe*(v+0.5)^2'
              write(out,'(A18)') 'v'
              do i=0, 3
              write(out,'(I18)',advance='no') i
              do istate=1,N
                 if(fl(istate)%zHasMinimum) then
                   write(out,'(F22.10)',advance='no') fl(istate)%Y00+fl(istate)%V1+fl(istate)%we*(real(i,rk)+0.5_rk) - &
                                                      fl(istate)%we*fl(istate)%xe*(real(i,rk)+0.5_rk)**2
                 else
                   write(out, '(a22)', advance='no') 'N.A.'
                 endif
              enddo
              write(out,*)
              enddo
              write(out,*)
              write(out,'(a18)', advance='no') 'Spin ='
              do istate=1,N
                   write(out,'(F22.1)',advance='no') fl(istate)%spini
              enddo
              write(out,*)
              write(out,'(a18)', advance='no') '|Lambda| ='
              do istate=1,N
                   write(out,'(F22.1)',advance='no') real( fl(istate)%Lambda, rk)
              enddo
              write(out,*)
              write(out,'(a18)', advance='no') 'Physical J_min ='
              do istate=1,N
                   write(out,'(F22.1)',advance='no') fl(istate)%Omega_min
              enddo
              write(out,*)
              write(out,*)
              write(out,'(A)') 'Approximate J=J_min vibrational-rotational energy levels (no couplings)'
              write(out,'(A)') ' given by E(v, J) = E(v, J=0) B0*J*(J+1)- ae*(v+0.5)*J*(J+1) - De*[J(J+1)]^2'
              write(out,'(A18)') 'v'
              do i=0, 3
              write(out,'(I18)',advance='no') i
              do istate=1,N
                 if(fl(istate)%zHasMinimum) then
                   write(out,'(F22.10)',advance='no')  fl(istate)%Y00+fl(istate)%V1+fl(istate)%we*(real(i,rk)+0.5_rk) &
                                                      -fl(istate)%we*fl(istate)%xe*(real(i,rk)+0.5_rk)**2             &
                                                      +fl(istate)%B0*fl(istate)%Omega_min*(fl(istate)%Omega_min+1) &
                                                      -fl(istate)%alphae*(real(i,rk)+0.5_rk)*fl(istate)%Omega_min* &
                                                                                          (fl(istate)%Omega_min+1) &
                                                      -fl(istate)%Debar*(fl(istate)%Omega_min*(fl(istate)%Omega_min+1))**2
                 else
                   write(out, '(a22)', advance='no') 'N.A.'
                 endif
              enddo
              write(out,*)
              enddo
              !
              write(out,*)
              !
        end subroutine check_and_print_field
        !
    end subroutine map_F1_hyperfine_fields_onto_grid

    subroutine construct_F1_hyperfine_constant_field_matrix
        ! <state, v| hfcc(R) |state', v'>
        implicit none
        TYPE(fieldT), POINTER :: field
        INTEGER(ik) :: ilevel, jlevel, index_object, index_field

        do index_object = 1, GLOBAL_NUM_HFCC_OBJECT
            do index_field = 1, hfcc(index_object)%num_field
                field => hfcc(index_object)%field(index_field)
                allocate(field%matelem(vibrational_totalroots,vibrational_totalroots),stat=alloc)

                do ilevel = 1, vibrational_totalroots
                    do jlevel = 1, ilevel
            
                    field%matelem(ilevel,jlevel)  = &
                        sum(vibrational_contrfunc(:,ilevel) &
                            *(field%gridvalue(:)) &
                            *vibrational_contrfunc(:,jlevel))
        
                    field%matelem(jlevel,ilevel) = field%matelem(ilevel,jlevel)
        
                    end do
                end do
            end do
        end do
    end subroutine construct_F1_hyperfine_constant_field_matrix

    subroutine construct_primitive_F1_basis
        ! |state, v; Lambda; S, Sigma; J, Omega; I1; F1>
        implicit none


        REAL(rk) :: F1, J
        INTEGER(ik) :: index_F1, Ndimen_F1
        INTEGER(ik) :: index_J, index_J_min, index_J_max, start_contr_F1, end_contr_F1
        
        ! something for J_global_max vs basis(:)

        num_primitive_F1_basis = 0
        do index_F1 = 1, num_F1
            F1 = F1_global_min + REAL(index_F1 - 1, rk)
            F1_list(index_F1) = F1
            J_min = abs(F1 - I1)
            J_max = F1 + I1

            index_J_min = nint(J_min - J_global_min) + 1
            index_J_max = nint(J_max - J_global_min) + 1

            Ndimen_F1 = 0
            do index_J = index_J_min, index_J_max
                Ndimen_F1 = Ndimen_F1 + basis(index_J)%Ndimen
            end do

            primitive_F1_basis(index_F1)%Ndimen = Ndimen_F1
            ALLOCATE (primitive_F1_basis(index_F1)%icontr(Ndimen_F1))

            num_primitive_F1_basis = &
                num_primitive_F1_basis + Ndimen_F1

            end_contr_F1 = 0
            do index_J = index_J_min, index_J_max
                J = J_global_min + REAL(index_J - 1, rk)
                start_contr_F1 = end_contr_F1
                end_contr_F1 = start_contr_F1 + basis(index_J)%Ndimen
                primitive_F1_basis(index_F1)%icontr(start_contr_F1 + 1:end_contr_F1) = &
                    basis(index_J)%icontr(:)
                primitive_F1_basis(index_F1)%icontr(start_contr_F1 + 1:end_contr_F1)%Jrot = J
            end do

            primitive_F1_basis(index_F1)%icontr(:)%index_F1 = index_F1
            primitive_F1_basis(index_F1)%icontr(:)%I1 = I1
            primitive_F1_basis(index_F1)%icontr(:)%F1 = F1
        end do
    end subroutine construct_primitive_F1_basis

    subroutine construct_parity_conserved_F1_matrix( &
        index_F1, index_represCs, &
        primitive_F1_hyperfine_matrix, parity_conserved_F1_matrix, transformation_matrix)
        ! < phi_m^{tau, J}; I; F| H_hfs + H0 |phi_n^{tau, J'}; I; F>
        implicit none

        INTEGER(ik), INTENT(IN) :: index_F1, index_represCs
        REAL(rk), ALLOCATABLE, INTENT(IN) :: primitive_F1_hyperfine_matrix(:,:)
        REAL(rk), ALLOCATABLE, INTENT(INOUT) :: parity_conserved_F1_matrix(:,:)
        REAL(rk), ALLOCATABLE, INTENT(INOUT) :: transformation_matrix(:,:)

        INTEGER(ik) :: index, Ndimen_F1, Nlevels_F1
        REAL(rk), ALLOCATABLE :: eigen_value_J_in_F1(:)
        REAL(rk), ALLOCATABLE :: intermediate_matrix(:,:)

        Ndimen_F1 = primitive_F1_basis(index_F1)%Ndimen
        Nlevels_F1 = Ndimen_F1 / sym%NrepresCs
            
        ALLOCATE(eigen_value_J_in_F1(Nlevels_F1))
      
        CALL construct_F1_transformation_matrix( &
            index_F1, index_represCs, &
            transformation_matrix, eigen_value_J_in_F1)

        ! ! H_hfs^tau = Phi^\dagger * H_hfs * Phi
        ! parity_conserved_F1_matrix = &
        !     matmul(matmul(transpose(transformation_matrix), &
        !                   primitive_F1_hyperfine_matrix), &
        !            transformation_matrix)

        ALLOCATE(intermediate_matrix(Nlevels_F1, Ndimen_F1))

        ! M_intermediate = Phi^\dagger * H_hfs
        CALL dgemm('T', 'N', Nlevels_F1, Ndimen_F1, Ndimen_F1, 1.0_rk, &
                    transformation_matrix, Ndimen_F1, &
                    primitive_F1_hyperfine_matrix, Ndimen_F1, &
                    0.0_rk, intermediate_matrix, Nlevels_F1)

        ! H_hfs^tau = H_indermediate * Phi           
        CALL dgemm('N', 'N', Nlevels_F1, Nlevels_F1, Ndimen_F1, 1.0_rk, &
                    intermediate_matrix, Nlevels_F1, &
                    transformation_matrix, Ndimen_F1, &
                    0.0_rk, parity_conserved_F1_matrix, Nlevels_F1)

        DEALLOCATE(intermediate_matrix)

        ! H^tau = H0^tau + H_hfs^tau 
        do index = 1, Nlevels_F1
            parity_conserved_F1_matrix(index, index) = &
                parity_conserved_F1_matrix(index, index) &
                + eigen_value_J_in_F1(index)
        end do    
        ! wavenumber_to_MHz = wavenumber_to_MHz
        DEALLOCATE(eigen_value_J_in_F1)
    end subroutine construct_parity_conserved_F1_matrix

    subroutine construct_F1_transformation_matrix( &
        index_F1, index_represCs, &
        transformation_matrix, &
        eigen_value_J_in_F1)
        ! Construct the basis tranformation matrix
        ! (technically transition matrix), Phi, 
        ! between the primitive F1 basis set
        ! and the parity conserved F1 basis set.
        ! We don't term Phi 'transition matrix' here
        ! to avoid misunderstanding although it seems unnecessary.
        ! The eigen values of H0 are obtained in the meanwhile.

        implicit none

        INTEGER(ik), INTENT(IN) :: index_F1, index_represCs
        REAL(rk), ALLOCATABLE, INTENT(INOUT) :: transformation_matrix(:,:), &
                                                eigen_value_J_in_F1(:)

        REAL(rk) :: F1, J
        INTEGER(ik) ::  index_J, index_J_min, index_J_max, &
                        start_dimen_F1, end_dimen_F1, &
                        start_level_F1, end_level_F1
        
        transformation_matrix = 0
        eigen_value_J_in_F1 = 0

        F1 = F1_global_min + REAL(index_F1 - 1, rk)
        J_min = abs(F1 - I1)
        J_max = F1 + I1

        index_J_min = nint(J_min - J_global_min) + 1
        index_J_max = nint(J_max - J_global_min) + 1
        
        end_dimen_F1 = 0
        end_level_F1 = 0
        do index_J = index_J_min, index_J_max
            
            J = J_global_min + REAL(index_J - 1, rk)   ! test
            
            start_dimen_F1 = end_dimen_F1
            start_level_F1 = end_level_F1
            end_dimen_F1 = start_dimen_F1 &
                + eigen(index_J, index_represCs)%Ndimen
            end_level_F1 = start_level_F1 &
                + eigen(index_J, index_represCs)%Nlevels

            ! Eigenvalues of H0
            eigen_value_J_in_F1(start_level_F1 + 1 : end_level_F1) &
                = eigen(index_J, index_represCs)%val

            ! Phi
            transformation_matrix(start_dimen_F1 + 1 : end_dimen_F1, &
                                  start_level_F1 + 1 : end_level_F1) &
                = eigen(index_J, index_represCs)%vect

            ! eigen_all_F1 is the table of eigen states
            ! corresponding to all F1 values.
            eigen_all_F1(index_F1, index_represCs)% &
                    quanta(start_level_F1 + 1 : end_level_F1) &
                = eigen(index_J, index_represCs)%quanta
            eigen_all_F1(index_F1, index_represCs)% &
                    quanta(start_level_F1 + 1 : end_level_F1)%F1 &
                = F1
            eigen_all_F1(index_F1, index_represCs)% &
                    quanta(start_level_F1 + 1 : end_level_F1)%I1 &
                = I1
            eigen_all_F1(index_F1, index_represCs)% &
                    quanta(start_level_F1 + 1 : end_level_F1)%index_F1 &
                = index_F1
            eigen_all_F1(index_F1, index_represCs)% &
                    quanta(start_level_F1 + 1 : end_level_F1)%Jrot &
                = J
        end do
    end subroutine construct_F1_transformation_matrix

    subroutine construct_primitive_F1_hyperfine_matrix( &
        index_F1, primitive_F1_hyperfine_matrix)
        ! <b, J; I; F| H_hfs |k, J'; I; F>
        implicit none

        INTEGER(ik), INTENT(IN) :: index_F1
        REAL(rk), ALLOCATABLE, INTENT(INOUT) :: primitive_F1_hyperfine_matrix(:, :)

        REAL(rk) :: F1, &
                    S_bra, S_ket, Sigma_bra, Sigma_ket, &
                    J_bra, J_ket, Omega_bra, Omega_ket
        INTEGER(ik) :: state_bra, state_ket, bra, ket, v_bra, v_ket, &
                       index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        INTEGER(ik) :: Ndimen_F1, index_hfcc, index_field

        primitive_F1_hyperfine_matrix = 0

        Ndimen_F1 = primitive_F1_basis(index_F1)%Ndimen

        do ket = 1, Ndimen_F1
            call get_quanta(index_F1, ket, &
                            F1, state_ket, index_v_ket, v_ket, Lambda_ket, &
                            S_ket, Sigma_ket, J_ket, Omega_ket)

            do bra = 1, Ndimen_F1
                call get_quanta(index_F1, bra, &
                                F1, state_bra, index_v_bra, v_bra, Lambda_bra, &
                                S_bra, Sigma_bra, J_bra, Omega_bra)

                if (state_bra == state_ket) then

                    do index_hfcc = 1, GLOBAL_NUM_HFCC_OBJECT
                        do index_field = 1, hfcc(index_hfcc)%num_field
                            if ( (state_bra == hfcc(index_hfcc)%field(index_field)%istate) &
                                .and. (state_ket == hfcc(index_hfcc)%field(index_field)%jstate)) then

                                primitive_F1_hyperfine_matrix(bra, ket) = &
                                    primitive_F1_hyperfine_matrix(bra, ket) &
                                    + primitive_F1_hyperfine_matrix_element( &
                                        F1, index_hfcc, index_field, & 
                                        index_v_bra, index_v_ket, &
                                        Lambda_bra, Lambda_ket, &
                                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                        J_bra, J_ket, Omega_bra, Omega_ket)
                            end if
                        end do                  
                    end do

                end if

                ! if (ket < bra) then
                !     primitive_F1_hyperfine_matrix(ket, bra) = &
                !         primitive_F1_hyperfine_matrix(bra, ket)
                ! end if
            end do
        end do
    end subroutine construct_primitive_F1_hyperfine_matrix

    subroutine get_quanta(index_F1, index_contr_F1, &
        F1, state, index_v, v, Lambda, S, Sigma, J, Omega)
        ! Get the quanta of index_F1 && index_contr_F1
        implicit none

        INTEGER(ik), INTENT(IN) :: index_F1, index_contr_F1
        INTEGER(ik), INTENT(OUT) :: state, Lambda, index_v, v
        REAL(rk), INTENT(OUT) :: F1, S, Sigma, J, Omega
        TYPE(quantaT), POINTER :: icontr
        
        icontr => primitive_F1_basis(index_F1)%icontr(index_contr_F1)

        F1 = icontr%F1

        state = icontr%istate
        Lambda = icontr%ilambda
        index_v = icontr%ivib
        v = icontr%v
        
        S = icontr%spin
        Sigma = icontr%sigma
        J = icontr%Jrot
        Omega = icontr%omega 

    end subroutine get_quanta
    
    function primitive_F1_hyperfine_matrix_element( &
        F1, index_hfcc, index_field, & 
        index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_hfs)

        implicit none

        REAL(rk), INTENT(IN) :: F1, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        INTEGER(ik), INTENT(IN) :: index_hfcc, index_field, &
                                   index_v_bra, index_v_ket, &
                                   Lambda_bra, Lambda_ket
        REAL(rk) :: H_hfs

        ! The switch flow could be refactored if 
        ! an abstract interface is defined.
        ! This function could be removed then.
        if (index_hfcc == 1) then
            H_hfs = H_FC_element( &
                        F1, index_field, index_v_bra, index_v_ket,&
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket) 
        else if (index_hfcc == 2) then
            H_hfs =  H_IL_element( &
                        F1, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket) 
        else if (index_hfcc == 3) then
            H_hfs = H_DIP_C_element( &
                        F1, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket) 
        else if (index_hfcc == 4) then 
            H_hfs = H_DIP_D_element( &
                        F1, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket) 
        else if (index_hfcc == 5) then
            H_hfs = H_IJ_element( &
                        F1, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket) 
        else if (index_hfcc == 6) then
            H_hfs = H_EQQ0_element( &        
                        F1, index_field, index_v_bra, index_v_ket,  &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket)
        else if (index_hfcc == 7) then
            H_hfs = H_EQQ2_element( &        
                        F1, index_field, index_v_bra, index_v_ket,  &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket)
        end if

    end function primitive_F1_hyperfine_matrix_element

    function get_sign(exp) result(sign_)
        ! Evaluate (-1)^exp and return a REAL(rk) value.
        implicit none

        REAL(rk), INTENT(IN) :: exp
        REAL(rk) sign_

        if (mod(nint(exp), 2) == 0) then
            sign_ = 1.0_rk
        else
            sign_ = -1.0_rk
        end if
    end function get_sign

    function H_FC_element( &
        F1, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_FC)

        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket,&
                                   Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: F1, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket 
        
        REAL(rk) :: H_FC, rank1_reduced_matrix_element_M1
        INTEGER(ik) :: q
        
        ! Fermi contact interaction matrix elemets
        ! See Kato 1993 Eq.(102)
        ! https://doi.org/10.1246/bcsj.66.3203

        rank1_reduced_matrix_element_M1 = 0.0_rk
        do q = -1, 1
            if ((Lambda_ket == Lambda_bra) &
                .and. (nint(Sigma_ket - (Sigma_bra - q)) == 0)) then
                rank1_reduced_matrix_element_M1 = &
                    get_sign(J_bra + S_bra - Lambda_bra - 2.0_rk*Sigma_bra) &
                    * SQRT((2.0_rk*S_bra + 1.0_rk)*(2.0_rk*J_ket + 1.0_rk)) &
                    * Wigner3j(S_ket, S_bra, 1.0_rk, &
                               Sigma_ket, -Sigma_bra, REAL(q, rk)) &
                    * Wigner3j(J_ket, J_bra, 1.0_rk, &
                               Omega_ket, -Omega_bra, REAL(q, rk)) &
                    * sqrt(S_bra * (S_bra + 1.0_rk)) &
                    * hfcc(1)%field(index_field)%matelem(index_v_bra, index_v_ket)
                    
                exit
            end if
        end do 

        H_FC = magnetic_dipole_hyperfine_matrix_element( &
            F1, J_bra, J_ket, rank1_reduced_matrix_element_M1)

        
        ! H_FC = 0.0_rk
        ! do q = -1, 1
        !     if ( Lambda_ket == Lambda_bra ) then
        !         H_FC = H_FC + &
        !                 FC_bF * get_sign(J_ket+F1+I1) * &
        !                 Wigner6j(I1, J_ket, F1, J_bra, I1, 1.0_rk) * &
        !                 sqrt(I1*(I1+1.0_rk)*(2.0_rk*I1+1.0_rk)) * &
        !                 get_sign(J_bra - Omega_bra) * &
        !                 Wigner3j(J_bra, 1.0_rk, J_ket, -Omega_bra, REAL(q, rk), Omega_ket) * &
        !                 sqrt((2.0_rk*J_bra+1.0_rk)*(2.0_rk*J_ket+1.0_rk)) *&
        !                 get_sign(S_bra-Sigma_bra) *&
        !                 Wigner3j(S_bra, 1.0_rk, S_bra, -Sigma_bra, REAL(q, rk), Sigma_ket) *&
        !                 sqrt(S_bra*(S_bra+1.0_rk)*(2.0_rk*S_bra+1.0_rk))
        !     end if
        ! end do
        ! wavenumber_to_MHz = wavenumber_to_MHz
    end function H_FC_element

    function H_IL_element( &
        F1, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_IL)

        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, &
                                   Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: F1, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket 

        REAL(rk) :: H_IL, rank1_reduced_matrix_element_M1
        INTEGER(ik) :: p

        ! Nuclear spin - electron orbital angular momentum interaction
        ! See Kato 1993 Eq.(104)
        ! https://doi.org/10.1246/bcsj.66.3203  

        rank1_reduced_matrix_element_M1 = 0.0_rk
        do p = -1, 1
            if ((Lambda_ket == (Lambda_bra - p)) &
            .and. (nint(S_bra - S_ket) == 0) &
            .and. (nint(Sigma_bra - Sigma_ket) == 0)) then
                rank1_reduced_matrix_element_M1 = &
                    get_sign(J_bra - Omega_bra) &
                    * sqrt(2.0_rk * J_ket + 1.0_rk) &
                    * Wigner3j(J_ket, J_bra, 1.0_rk, &
                                Omega_ket, -Omega_bra, REAL(p, rk)) &
                    * Lambda_bra &
                    * hfcc(2)%field(index_field)%matelem(index_v_bra, index_v_ket)
                exit              
            end if
        end do
        H_IL = magnetic_dipole_hyperfine_matrix_element( &
                    F1, J_bra, J_ket, rank1_reduced_matrix_element_M1)
    end function H_IL_element

    function H_DIP_C_element( &
        F1, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_IS)

        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: F1, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket 

        REAL(rk) :: H_IS, rank1_reduced_matrix_element_M1
        INTEGER(ik) :: p, q, q1, q2

        ! Nuclear spin - electron spin dipole interaction
        ! See Kato 1993 Eq.(106)
        ! See Brown and Carrington Eq.(8.237)
        ! https://doi.org/10.1246/bcsj.66.3203

        rank1_reduced_matrix_element_M1 = 0.0_rk
        p = 0
        do q = -1, 1, 1
            if ( (Lambda_ket == (Lambda_bra - p)) &
            .and. (nint(Sigma_ket - (Sigma_bra - q)) == 0) ) then
                rank1_reduced_matrix_element_M1 = &
                    get_sign(J_bra + S_bra - Lambda_bra + p + q) &
                    * sqrt((2.0_rk * S_bra + 1.0_rk) * (2.0_rk * J_ket + 1.0_rk)) &
                    * Wigner3j(1.0_rk, 1.0_rk, 2.0_rk, &
                                REAL(-q, rk), REAL(p+q, rk), REAL(-p, rk)) &
                    * Wigner3j(S_bra, 1.0_rk, S_ket, &
                                -Sigma_bra, REAL(q, rk), Sigma_ket) &
                    * Wigner3j(J_ket, J_bra, 1.0_rk, &
                                Omega_ket, -Omega_bra, REAL(p+q, rk)) &
                    * sqrt(30.0_rk) * sqrt(S_bra * (S_bra + 1.0_rk)) &
                    * hfcc(3)%field(index_field)%matelem(index_v_bra, index_v_ket) / (-3.0_rk)
                exit 
            end if
        end do        

        H_IS = magnetic_dipole_hyperfine_matrix_element(&
            F1, J_bra, J_ket, rank1_reduced_matrix_element_M1)

        ! H_IS = 0
        ! do q = -1,1
        !     do q1 = -1, 1, 1
        !         do q2 = 0,0
        !             H_IS = H_IS + H_IS_D() * get_sign(J_ket+F1+I1) *&
        !                 Wigner6j(I1, J_ket, F1, J_bra, I1, 1.0_rk) *&
        !                 sqrt(I1*(I1+1.0_rk)*(2.0_rk*I1+1.0_rk)) *&
        !                 get_sign(J_bra-Omega_bra) *&
        !                 Wigner3j(J_bra, 1.0_rk, J_ket, -Omega_bra, REAL(q,rk), Omega_ket) *&
        !                 sqrt(30.0_rk*(2.0_rk*J_bra+1.0_rk)*(2.0_rk*J_ket+1.0_rk)) *&
        !                 get_sign(REAL(q,rk)) * &
        !                 Wigner3j(1.0_rk,2.0_rk,1.0_rk,REAL(q1,rk),REAL(q2,rk),REAL(-q,rk)) *&
        !                 get_sign(S_bra-Sigma_bra) *&
        !                 Wigner3j(S_bra,1.0_rk,S_ket,-Sigma_bra,REAL(q1,rk),Sigma_ket) *&
        !                 sqrt(S_bra*(S_bra+1.0_rk)*(2.0_rk*S_bra+1.0_rk))
        !         end do
        !     end do
        ! end do

    end function H_DIP_C_element

    function H_DIP_D_element( &
        F1, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_IS)

        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: F1, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket 

        REAL(rk) :: H_IS, rank1_reduced_matrix_element_M1
        INTEGER(ik) :: p, q, q1, q2

        ! Nuclear spin - electron spin dipole interaction
        ! See Kato 1993 Eq.(106)
        ! See Brown and Carrington Eq.(8.237)
        ! https://doi.org/10.1246/bcsj.66.3203

        rank1_reduced_matrix_element_M1 = 0.0_rk
        outer: do p = -2, 2, 4
            do q = -1, 1, 1
                if ( (Lambda_ket == (Lambda_bra - p)) &
                .and. (nint(Sigma_ket - (Sigma_bra - q)) == 0) ) then
                    rank1_reduced_matrix_element_M1 = &
                        get_sign(J_bra + S_bra - Lambda_bra + p + q) &
                        * sqrt((2.0_rk * S_bra + 1.0_rk) * (2.0_rk * J_ket + 1.0_rk)) &
                        * Wigner3j(1.0_rk, 1.0_rk, 2.0_rk, &
                                    REAL(-q, rk), REAL(p+q, rk), REAL(-p, rk)) &
                        * Wigner3j(S_bra, 1.0_rk, S_ket, &
                                    -Sigma_bra, REAL(q, rk), Sigma_ket) &
                        * Wigner3j(J_ket, J_bra, 1.0_rk, &
                                    Omega_ket, -Omega_bra, REAL(p+q, rk)) &
                        * sqrt(30.0_rk) * sqrt(S_bra * (S_bra + 1.0_rk)) &
                        * hfcc(4)%field(index_field)%matelem(index_v_bra, index_v_ket) / sqrt(6.0_rk)
                    exit outer
                end if
            end do            
        end do outer

        H_IS = magnetic_dipole_hyperfine_matrix_element(&
            F1, J_bra, J_ket, rank1_reduced_matrix_element_M1)

        ! H_IS = 0
        ! do q = -1,1
        !     do q1 = -1, 1, 1
        !         do q2 = 0,0
        !             H_IS = H_IS + H_IS_D() * get_sign(J_ket+F1+I1) *&
        !                 Wigner6j(I1, J_ket, F1, J_bra, I1, 1.0_rk) *&
        !                 sqrt(I1*(I1+1.0_rk)*(2.0_rk*I1+1.0_rk)) *&
        !                 get_sign(J_bra-Omega_bra) *&
        !                 Wigner3j(J_bra, 1.0_rk, J_ket, -Omega_bra, REAL(q,rk), Omega_ket) *&
        !                 sqrt(30.0_rk*(2.0_rk*J_bra+1.0_rk)*(2.0_rk*J_ket+1.0_rk)) *&
        !                 get_sign(REAL(q,rk)) * &
        !                 Wigner3j(1.0_rk,2.0_rk,1.0_rk,REAL(q1,rk),REAL(q2,rk),REAL(-q,rk)) *&
        !                 get_sign(S_bra-Sigma_bra) *&
        !                 Wigner3j(S_bra,1.0_rk,S_ket,-Sigma_bra,REAL(q1,rk),Sigma_ket) *&
        !                 sqrt(S_bra*(S_bra+1.0_rk)*(2.0_rk*S_bra+1.0_rk))
        !         end do
        !     end do
        ! end do


    end function H_DIP_D_element

    function H_IJ_element( &
        F1, index_field, index_v_bra, index_v_ket,  &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_IJ)

        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: F1, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket 
        
        REAL(rk) :: H_IJ, rank1_reduced_matrix_element_M1

        ! Nuclear spin-rotation interaction
        ! See Brown and Carrington (2003) Eq.(8.20)
        ! https://doi.org/10.1017/cbo9780511814808

        H_IJ = 0.0_rk

        if ( (Lambda_bra == Lambda_ket) &
            .and. (nint(S_bra - S_ket) == 0) &
            .and. (nint(Sigma_bra - Sigma_ket) == 0) &
            .and. (nint(J_bra - J_ket) == 0) ) then
            H_IJ = get_sign(J_bra + F1 + I1) &
                * Wigner6j(I1, J_bra, F1, J_bra, I1, 1.0_rk) &
                * sqrt(J_bra * (J_bra + 1.0_rk) * (2.0_rk * J_bra + 1.0_rk) &
                        * I1 * (I1 + 1.0_rk) * (2.0_rk * I1 + 1.0_rk)) &
                * hfcc(5)%field(index_field)%matelem(index_v_bra, index_v_ket)
        end if

    end function H_IJ_element
    
    function H_EQQ0_element( &        
        F1, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_EQ)
        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: F1, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket 
        
        REAL(rk) :: H_EQ, rank2_reduced_matrix_element_V1
        INTEGER(ik) :: p

        ! Nuclear electric quadrupole interaction
        ! See Kato 1993 Eq.(109)
        ! https://doi.org/10.1246/bcsj.66.3203
        ! The rule here is stronger than that of 
        ! Brown and Carrington (2003), Eq.(8.382)
        ! https://doi.org/10.1017/cbo9780511814808
        ! rank2_reduced_matrix_element_V1 = 0        
        ! do p = -2, 2, 2
        !     if ( (Lambda_ket == Lambda_bra - p) & 
        !     .and. (nint(S_ket - S_bra) == 0) &
        !     .and. (nint(Sigma_ket - Sigma_bra) == 0) ) then
        !         rank2_reduced_matrix_element_V1 = &
        !             get_sign(J_bra - Omega_bra) &
        !             * sqrt(2.0_rk * J_ket + 1.0_rk) &
        !             * Wigner3j(J_ket, J_bra, 2.0_rk, &
        !                         Omega_ket, -Omega_bra, REAL(p,rk)) &
        !             * H_EQ_E()
        !         exit
        !     end if
        ! end do

        ! H_EQ = electric_quadrupole_hyperfine_matrix_element(F1, &
        !     J_bra, J_ket, rank2_reduced_matrix_element_V1)
        ! See Brown and Carrington (2003) Eq.(9.14) Eq.(8.382)
        ! https://doi.org/10.1017/cbo9780511814808
        H_EQ = 0.0_rk
        if ( I1 < 1.0_rk ) return
        p = 0
        if ( (Lambda_ket == Lambda_bra - p) & 
        .and. (nint(S_ket - S_bra) == 0) &
        .and. (nint(Sigma_ket - Sigma_bra) == 0) ) then
            H_EQ = get_sign(J_ket + I1 + F1) &
                * Wigner6j(I1, J_bra, F1, J_ket, I1, 2.0_rk) &
                * (-0.5_rk) & ! eQ is absorbed by eQq
                / Wigner3j(I1, 2.0_rk, I1, -I1, 0.0_rk, I1)&
                * get_sign(J_bra - Omega_bra) &
                * Wigner3j(J_bra, 2.0_rk, J_ket, -Omega_bra, REAL(p, rk), Omega_ket) &
                * sqrt((2.0_rk * J_bra + 1.0_rk) * (2.0_rk * J_ket + 1.0_rk)) &
                * hfcc(6)%field(index_field)%matelem(index_v_bra, index_v_ket) / (-2.0_rk)
        end if
    end function H_EQQ0_element

    function H_EQQ2_element( &        
        F1, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_EQ)
        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: F1, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket 
        
        REAL(rk) :: H_EQ, rank2_reduced_matrix_element_V1
        INTEGER(ik) :: p

        ! Nuclear electric quadrupole interaction
        ! See Kato 1993 Eq.(109)
        ! https://doi.org/10.1246/bcsj.66.3203
        ! The rule here is stronger than that of 
        ! Brown and Carrington (2003), Eq.(8.382)
        ! https://doi.org/10.1017/cbo9780511814808
        ! rank2_reduced_matrix_element_V1 = 0        
        ! do p = -2, 2, 2
        !     if ( (Lambda_ket == Lambda_bra - p) & 
        !     .and. (nint(S_ket - S_bra) == 0) &
        !     .and. (nint(Sigma_ket - Sigma_bra) == 0) ) then
        !         rank2_reduced_matrix_element_V1 = &
        !             get_sign(J_bra - Omega_bra) &
        !             * sqrt(2.0_rk * J_ket + 1.0_rk) &
        !             * Wigner3j(J_ket, J_bra, 2.0_rk, &
        !                         Omega_ket, -Omega_bra, REAL(p,rk)) &
        !             * H_EQ_E()
        !         exit
        !     end if
        ! end do

        ! H_EQ = electric_quadrupole_hyperfine_matrix_element(F1, &
        !     J_bra, J_ket, rank2_reduced_matrix_element_V1)
        ! See Brown and Carrington (2003) Eq.(9.14) Eq.(8.382)
        ! https://doi.org/10.1017/cbo9780511814808
        H_EQ = 0.0_rk
        if ( I1 < 1.0_rk ) return
        do p = -2, 2, 4
            if ( (Lambda_ket == Lambda_bra - p) & 
            .and. (nint(S_ket - S_bra) == 0) &
            .and. (nint(Sigma_ket - Sigma_bra) == 0) ) then
                H_EQ = get_sign(J_ket + I1 + F1) &
                    * Wigner6j(I1, J_bra, F1, J_ket, I1, 2.0_rk) &
                    * (-0.5_rk) & ! eQ is absorbed by eQq
                    / Wigner3j(I1, 2.0_rk, I1, -I1, 0.0_rk, I1)&
                    * get_sign(J_bra - Omega_bra) &
                    * Wigner3j(J_bra, 2.0_rk, J_ket, -Omega_bra, REAL(p, rk), Omega_ket) &
                    * sqrt((2.0_rk * J_bra + 1.0_rk) * (2.0_rk * J_ket + 1.0_rk)) &
                    * hfcc(7)%field(index_field)%matelem(index_v_bra, index_v_ket) / (-sqrt(24.0_rk))
                exit
            end if
        end do
    end function H_EQQ2_element

    function Wigner6j(j1, j2, j3, J_1, J_2, J_3) result(w6j)
        implicit none
        REAL(rk), INTENT(IN) :: j1, j2, j3, J_1, J_2, J_3
        REAL(rk) :: w6j, threshold, tri1, tri2, tri3, tri4
        REAL(rk) :: a1, a2, a3, a4, k1, k2, k3, upper_limit, lower_limit, t

        tri1 = triangle_coefficient(j1, j2, j3) 
        tri2 = triangle_coefficient(j1, J_2, J_3) 
        tri3 = triangle_coefficient(J_1, j2, J_3) 
        tri4 = triangle_coefficient(J_1, J_2, j3)
        threshold = safe_min
        if ( (abs(tri1) < threshold) &
            .or. (abs(tri2) < threshold) &
            .or. (abs(tri3) < threshold) &
            .or. (abs(tri4) < threshold) ) then
            w6j = 0.0_rk
            return
        end if

        a1 = j1 + j2 + j3;
        a2 = j1 + J_2 + J_3;
        a3 = J_1 + j2 + J_3;
        a4 = J_1 + J_2 + j3;
        
        k1 = j1 + j2 + J_1 + J_2;
        k2 = j2 + j3 + J_2 + J_3;
        k3 = j3 + j1 + J_3 + J_1;

        lower_limit = max(a1, a2, a3, a4)
        upper_limit = min(k1, k3, k3)

        t = lower_limit
        threshold = 1E-6
        
        w6j = 0.0_rk
        do while (t < upper_limit + threshold)
            w6j = w6j + get_sign(t) * exp(log_())
            t = t + 1.0_rk
        end do

        w6j = sqrt(tri1 * tri2 * tri3 * tri4) * w6j

    contains
        REAL(rk) function log_()   
            log_ = log_gamma(t + 2.0_rk) &
                - log_gamma(t - j1 - j2 - j3 + 1.0_rk) &
                - log_gamma(t - j1 - J_2 - J_3 + 1.0_rk) & 
                - log_gamma(t - J_1 - j2 - J_3 + 1.0_rk) &
                - log_gamma(t - J_1 - J_2 - j3 + 1.0_rk) &
                - log_gamma(j1 + j2 + J_1 + J_2 - t + 1.0_rk) &
                - log_gamma(j2 + j3 + J_2 + J_3 - t + 1.0_rk) &
                - log_gamma(j3 + j1 + J_3 + J_1 - t + 1.0_rk)        
        end function log_

        function triangle_coefficient(a, b, c) result(tri)
            implicit none
            REAL(rk), INTENT(IN) :: a, b, c
            REAL(rk) :: tri
            integer(ik) :: xa
            
            tri = 0
            if ( a < 0 .or. b < 0 .or. c < 0 ) return

            do xa = nint(abs(a-b) * 2.0_rk), nint((a+b) * 2.0_rk), 2
                if ( nint(c * 2.0_rk) == xa) then
                    tri = factorial(a + b - c) &
                        * factorial(a + c - b) &
                        * factorial(b + c - a) &
                        / factorial(a + b + c + 1.0_rk)
                end if
            end do
        end function triangle_coefficient

    end function Wigner6j
    
    function factorial(n) result(f)
        implicit none
        REAL(rk) :: n
        REAL(rk) :: f

        f = GAMMA(n+1.0_rk)
    end function factorial

    function Wigner3j(j1, j2, j3, m_1, m_2, m_3) result(w3j)
        implicit none
        REAL(RK), INTENT(IN) :: j1, j2, j3, m_1, m_2, m_3
        REAL(rk) :: w3j

        w3j = three_j(j1, j2, j3, m_1, m_2, m_3)
    end function Wigner3j

    function magnetic_dipole_hyperfine_matrix_element( &
        F1, &
        J_bra, J_ket, rank1_reduced_matrix_element_M1) &
        result(H_MD1)   
        implicit none

        REAL(rk), INTENT(IN) :: F1, J_bra, J_ket, &
                                rank1_reduced_matrix_element_M1
        REAL(rk) :: H_MD1, rank1_reduced_matrix_element_I1

        ! See Kato 1993 Eq.(108)
        ! https://doi.org/10.1246/bcsj.66.3203
        rank1_reduced_matrix_element_I1 = sqrt(I1 * (I1 + 1.0_rk))

        ! See Kato 1993 Eq.(88)
        H_MD1 = get_sign(J_ket + I1 + F1) &
               * Wigner6j(J_bra, J_ket, 1.0_rk, &
                            I1, I1, F1) &
               * SQRT((2.0_rk*J_bra + 1.0_rk)*(2.0_rk*I1 + 1.0_rk)) &
               * rank1_reduced_matrix_element_M1 &
               * rank1_reduced_matrix_element_I1   
    end function magnetic_dipole_hyperfine_matrix_element

    function electric_quadrupole_hyperfine_matrix_element( &
        F1, &
        J_bra, J_ket, rank2_reduced_matrix_element_V1) &
        result(H_EQ1)

        REAL(rk), INTENT(IN) :: J_bra, J_ket, F1, rank2_reduced_matrix_element_V1
        REAL(rk) :: rank2_reduced_matrix_element_Q1, H_EQ1

        ! See Kato 1993 Eq.(111)
        ! https://doi.org/10.1246/bcsj.66.3203
        ! The equation should be further verified 
        ! eQ is absorbed by eQq
        rank2_reduced_matrix_element_Q1 = &
            0.5_rk * sqrt( (2.0_rk * I1 + 3.0_rk)*(I1 + 1.0_rk) &
                                / I1 * (2.0_rk*I1 - 1.0_rk) )
        
        ! See Kato 1993 Eq.(89)
        ! https://doi.org/10.1246/bcsj.66.3203
        H_EQ1 = get_sign(J_bra + I1 + F1) &
            * Wigner6j(J_bra, J_ket, 2.0_rk, &
                        I1, I1, F1) &
            * SQRT((2.0_rk*J_bra + 1.0_rk)*(2.0_rk*I1 + 1.0_rk)) &
            * rank2_reduced_matrix_element_V1 &
            * rank2_reduced_matrix_element_Q1
    end function electric_quadrupole_hyperfine_matrix_element

end module F1_hyperfine
