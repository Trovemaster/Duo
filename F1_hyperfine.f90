module F1_hyperfine

    use accuracy
    use lapack
    use symmetry, only : sym, correlate_to_Cs 
    use diatom_module, only: basis, three_j, faclog, eigen, jmax, job,&
                            vibrational_totalroots, vibrational_contrfunc, vibrational_quantum_number, &
                            hfcc1, hfcc2, hfcc12, hyperfine_setup, GLOBAL_NUM_HFCC_OBJECT, GLOBAL_NUM_HFCC12_OBJECT, poten, &
                            hyperfine_case, &
                            eigenT, basisT, quantaT, fieldT, FieldListT, linkT

    implicit none
    
    
    TYPE(basisT), ALLOCATABLE :: primitive_hfs_basis(:)
    TYPE(FieldListT), POINTER :: active_hfcc(:) => hfcc1
    
    REAL(rk), ALLOCATABLE :: F1_list(:)
    REAL(rk), ALLOCATABLE :: Itot_list(:)
    REAL(rk), ALLOCATABLE :: Ftot_list(:)
    
    INTEGER(ik) :: num_F1, num_Ftot, num_primitive_hfs_basis, num_represCs
    REAL(rk) :: F1_global_min, F1_global_max, Ftot_global_min, Ftot_global_max, I1, I2, Itot_max
    INTEGER(ik) :: alloc
    REAL(rk) :: J_global_min = 0.5_rk, J_global_max
    REAL(rk) :: J_min, J_max, MHz_to_wavenumber = 1.0E6_rk / vellgt, &
                wavenumber_to_MHz = vellgt / 1.0E6_rk

    TYPE(eigenT), allocatable :: eigen_all_hfs(:,:)
    REAL(rk) :: min_eigen_value_hfs
    ! LOGICAL :: F1_hyperfine_gridvalue_allocated  = .false.
    INTEGER(ik), PARAMETER :: unit_hyperfine_states = 65
    
contains

    subroutine hyperfine_structure(iverbose)   ! control different nuclear-spin cases
        implicit none
        integer(ik), INTENT(IN) :: iverbose
        
        write(out, '(/A)') "Start: hyperfine energies calculation"
        open(unit=unit_hyperfine_states, file="hyperfine.states")

        active_hfcc => hfcc1
        select case(hyperfine_case)
        case(1)
            write(out, '(/A)') " Hyperfine case: single non-zero-spin nucleus."
            call F1_hyperfine_structure(iverbose)
        case(2)
            write(out, '(/A)') " Hyperfine case: two non-zero-spin nuclei; homonuclear."
            call IF_hyperfine_structure(iverbose)
        case(3)
            write(out, '(/A)') " Hyperfine case: heteronuclear two-spin coupling."
            call F1F_hyperfine_structure(iverbose)
        case default
            stop 'ERROR: hyperfine nuclear-spin case has not been initialized.'
        end select
    end subroutine hyperfine_structure

    subroutine IF_hyperfine_structure(iverbose)
        implicit none
        integer(ik), INTENT(IN) :: iverbose
        INTEGER(ik) :: index_J, index_J_min, index_J_max, index_Itot  
        INTEGER(ik) :: ilevel, istate 
        INTEGER(ik) :: index_Ftot, index_represCs, index_represC2v, ilabel_gu, iroot, iF_ID
        INTEGER(ik) :: Ndimen_Ftot, Nlevels_Ftot, index_levels_Ftot, maxlocation, maxlocation2
        REAL(rk) :: Ftot
        REAL(rk), ALLOCATABLE :: eigen_value_Ftot(:), &
                                    parity_conserved_IF_matrix(:, :), &
                                    primitive_IF_hyperfine_matrix(:,:), &
                                    transformation_matrix(:,:), & 
                                    eigen_vector_Ftot(:, :), &
                                    eigen_value_rov_J(:) 
        TYPE(quantaT), POINTER :: quanta

        I1 = hyperfine_setup%I1
        
        Itot_max = 2.0_rk*I1   ! = 2*I1, integer

        allocate(Itot_list(int(Itot_max) + 1)) 
        do index_Itot = 1, int(Itot_max) + 1 
            Itot_list(index_Itot) = real(index_Itot-1, rk)   ! 0, 1, 2, ..., Itot_max 
        end do 

        J_global_max = jmax
        if (J_global_max < Itot_max) then
            write (out, *) "J_global_max < Itot_max"
            stop
        end if

        Ftot_global_max = J_global_max - Itot_max

        if (mod(nint(Ftot_global_max*2.0_rk), 2) == 0) then
            Ftot_global_min = 0.0_rk
        else
            Ftot_global_min = 0.5_rk
        end if

        ! The basis of J always starts from the minimum of J
        if (mod(nint(J_global_max*2.0_rk), 2) == 0) then
            J_global_min = 0.0_rk
        else
            J_global_min = 0.5_rk
        end if

        num_Ftot = nint(Ftot_global_max - Ftot_global_min) + 1
        
        CALL construct_hyperfine_constant_field_matrix(iverbose)
        

        ALLOCATE(Ftot_list(num_Ftot))
        ALLOCATE(primitive_hfs_basis(num_Ftot))

        CALL construct_primitive_IF_basis(iverbose)

        if (iverbose>=4) then
            write(out, '(/A2, A)') '', 'Start: hyperfine states'
            write(out, '(/A2, A, F22.12)') '', 'Zero point energy =', job%ZPE
        endif

        ALLOCATE(eigen_all_hfs(num_Ftot, sym%Nrepresen)) ! ! sym%Nrepresen=4, see symmetry.f90
        eigen_all_hfs(:,:)%Nlevels = 0 ! This is used in the intensity module to skip empty blocks

        iroot = 0
        min_eigen_value_hfs = safe_max
        do index_Ftot = 1, num_Ftot
            Ftot = Ftot_list(index_Ftot)
            if (iverbose>=4) then
                write(out, '(/A4, A, F6.1)') '', 'Eigen states of F =', Ftot
                write(out, '(/A6, A)') '', 'Construct primitive hyperfine matrix'
            endif

            Ndimen_Ftot = primitive_hfs_basis(index_Ftot)%Ndimen
            ALLOCATE(primitive_IF_hyperfine_matrix(Ndimen_Ftot, Ndimen_Ftot))
            CALL construct_primitive_IF_hyperfine_matrix( &
                index_Ftot, primitive_IF_hyperfine_matrix)

            if (iverbose>=4) write(out, '(A6, A/)') '', '... done'

            do index_represCs = 1, sym%NrepresCs 

                do ilabel_gu = -1, 1, 2 ! -1 for u and 1 for g 

                    index_represC2v = correlate_to_Cs(index_represCs, ilabel_gu)

                    if (iverbose>=4) then 
                        write(out, '(/A6, A, A3)') '', 'Eigen states of symmetry block ', trim(parity_sign(index_represC2v))
                    endif

                    Nlevels_Ftot = 0
                    do index_Itot = 1, size(Itot_list) 

                        if (index_represC2v <= 2) then  ! A1/A2 
                            if (mod(nint(Itot_list(index_Itot)), 2) /= 0) cycle   ! Itot = even 
                        else ! B1/B2
                            if (mod(nint(Itot_list(index_Itot)), 2) /= 1) cycle   ! Itot = odd
                        end if

                        J_min = abs(Ftot - Itot_list(index_Itot))
                        J_max = Ftot + Itot_list(index_Itot)
                        index_J_min = nint(J_min - J_global_min) + 1
                        index_J_max = nint(J_max - J_global_min) + 1

                        do index_J = index_J_min, index_J_max
                            if (eigen(index_J, index_represCs)%Nlevels > 0) then
                                do ilevel = 1, eigen(index_J, index_represCs)%Nlevels
                                    istate = eigen(index_J, index_represCs)%quanta(ilevel)%istate
                                    if (poten(istate)%parity%gu /= ilabel_gu) cycle
                                    Nlevels_Ftot = Nlevels_Ftot + 1
                                end do
                            else
                                write(out,'(A8, A, F5.1, A, A, A)') '', 'The current block with (J = ', &
                                    J_global_min + REAL(index_J - 1, rk), &
                                    ', symmetry = ', trim(parity_sign(index_represC2v)), ') is empty.'
                            end if
                        end do
                    end do

                    if (Nlevels_Ftot == 0) then
                        write(out,'(/A8, A, F5.2, A, A, A)') '', 'The current block with (F = ', Ftot, &
                            ', symmetry = ', trim(parity_sign(index_represC2v)), ') is empty.'
                        cycle     ! Skip empty F-symmetry block for all Itot
                    end if

                    ALLOCATE(eigen_value_Ftot(Nlevels_Ftot))
                    ALLOCATE(parity_conserved_IF_matrix(Nlevels_Ftot, Nlevels_Ftot))
                    ALLOCATE(transformation_matrix(Ndimen_Ftot, Nlevels_Ftot))
                    ALLOCATE(eigen_vector_Ftot(Ndimen_Ftot, Nlevels_Ftot))
                    ALLOCATE(eigen_value_rov_J(Nlevels_Ftot))

                    ALLOCATE(eigen_all_hfs(index_Ftot, index_represC2v)%quanta(Nlevels_Ftot))
                    ALLOCATE(eigen_all_hfs(index_Ftot, index_represC2v)%val(Nlevels_Ftot))
                    ALLOCATE(eigen_all_hfs(index_Ftot, index_represC2v)%vect(Ndimen_Ftot, Nlevels_Ftot))

                    if (iverbose>=4) then
                        write(out, '(/A8, A)') &
                            '', 'Construct parity conserved hyperfine matrix'
                    endif

                    CALL construct_parity_conserved_IF_matrix( &
                        index_Ftot, index_represC2v, index_represCs, ilabel_gu, &
                        Ndimen_Ftot, Nlevels_Ftot, primitive_IF_hyperfine_matrix, &
                        parity_conserved_IF_matrix, transformation_matrix, eigen_value_rov_J) 
                        ! Now, parity_conserved_IF_matrix(:, :) is the sum of parity-conserved hyperfine matrix 
                        ! under the rovibronic basis set and the eigenvalues of rovibronic matrix.
                        ! transformation_matrix(:, :) holds the rovibronic eigenvectors of Ftot, viz Phi
                    if (iverbose>=4) write(out, '(A8, A/)') '', '... done'

                    if (iverbose>=4) then
                        write(out, '(/A8, A)') &
                            '', 'Diagonalize parity conserved hyperfine matrix'
                    endif
                    CALL lapack_syev(parity_conserved_IF_matrix, eigen_value_Ftot)
                    ! Now, parity_conserved_IF_matrix(:, :) is the eigenvector matrix, viz U, under the rovibronic wavefunction.
                    if (iverbose>=4) write(out, '(A8, A/)') '', '... done'
                    
                    if (iverbose>=4) then
                        write(out, '(/A8, A)') &
                            '', "Transform wavefunction to Hund's case (a_beta) basis set"
                    endif
                    ! Psi = Phi * U
                    ! eigen_vector_Ftot = matmul(transformation_matrix, &
                    !                          parity_conserved_IF_matrix)
                    CALL dgemm('N', 'N', Ndimen_Ftot, Nlevels_Ftot, Nlevels_Ftot, 1.0_rk, &
                                transformation_matrix, Ndimen_Ftot, &
                                parity_conserved_IF_matrix, Nlevels_Ftot, &
                                0.0_rk, eigen_vector_Ftot, Ndimen_Ftot)
                    ! Now, eigen_vector_Ftot(:, :) is the eigenvector matrix under the Hund's case(a_beta) basis set.
                    if (iverbose>=4) write(out, '(A8, A/)') '', '... done'

                    eigen_all_hfs(index_Ftot, index_represC2v)%Nlevels = Nlevels_Ftot
                    eigen_all_hfs(index_Ftot, index_represC2v)%Ndimen = Ndimen_Ftot
                    eigen_all_hfs(index_Ftot, index_represC2v)%val = eigen_value_Ftot
                    eigen_all_hfs(index_Ftot, index_represC2v)%vect = eigen_vector_Ftot
                    
                    if (iverbose>=4) then
                        write(out, '(/A8, A)') '', 'Calculate energy levels'
                    endif
                    eigen_value_Ftot(:) = eigen_value_Ftot(:) - job%ZPE

                    iF_ID = 0
                    do index_levels_Ftot = 1, Nlevels_Ftot
                        maxlocation = MAXLOC(ABS(eigen_vector_Ftot(:, index_levels_Ftot)), DIM=1)
                        maxlocation2  = maxloc(abs(parity_conserved_IF_matrix(:, index_levels_Ftot)), dim=1) 
                        quanta => primitive_hfs_basis(index_Ftot)%icontr(maxlocation)

                        iroot = iroot + 1
                        iF_ID = iF_ID + 1
                        eigen_all_hfs(index_Ftot, index_represC2v)%quanta(index_levels_Ftot)%iroot = iroot
                        eigen_all_hfs(index_Ftot, index_represC2v)%quanta(index_levels_Ftot)%iF_ID = iF_ID

                        min_eigen_value_hfs = min(&
                            eigen_all_hfs(index_Ftot, index_represC2v)%val(index_levels_Ftot), &
                            min_eigen_value_hfs)

                        if (iverbose>=4) then
                            write(unit_hyperfine_states, '(I7, 3F22.12, I7, A7, A4, A5, A7, A12, I7, I4, 2A5)') & 
                                iroot, &
                                eigen_value_Ftot(index_levels_Ftot), &
                                eigen_value_rov_J(maxlocation2)-job%ZPE, & 
                                eigen_value_Ftot(index_levels_Ftot)-(eigen_value_rov_J(maxlocation2)-job%ZPE), & 
                                int(quanta%Ftot * 2.0_rk + 1.0_rk), &
                                qnum2a(quanta%Ftot, 7), &
                                trim(parity_sign(index_represC2v)), & 
                                qnum2a(quanta%Itot, 5), &  
                                qnum2a(quanta%Jrot, 7), &
                                trim(poten(quanta%istate)%name), &
                                eigen_all_hfs(index_Ftot, index_represC2v)%quanta(maxlocation2)%v, &
                                quanta%ilambda, & 
                                qnum2a(quanta%sigma, 5), &
                                qnum2a(quanta%omega, 5) 
                        endif
                    enddo 
                    if (iverbose>=4) write(out, '(A8, A/)') '', '... done'

                    DEALLOCATE(eigen_vector_Ftot)
                    DEALLOCATE(transformation_matrix)
                    DEALLOCATE(parity_conserved_IF_matrix)
                    DEALLOCATE(eigen_value_Ftot)
                    DEALLOCATE(eigen_value_rov_J) 
                end do
                
            end do

            DEALLOCATE(primitive_IF_hyperfine_matrix)
        end do

        if (iverbose>=4) write(out, '(A2, A/)') '', 'End: hyperfine states'    
        write(out, '(A/)') "End: hyperfine energy calculation"
        
        close(unit=unit_hyperfine_states)
    contains
        function parity_sign(index_represC2v)
            INTEGER(ik), INTENT(IN) :: index_represC2v
            CHARACTER(4) :: parity_sign
            if ( index_represC2v == 1 ) then
                parity_sign = '+g'
            elseif ( index_represC2v == 2 ) then
                parity_sign = '-u'
            elseif ( index_represC2v == 3 ) then
                parity_sign = '-g'
            elseif ( index_represC2v == 4 ) then
                parity_sign = '+u' 
            endif
        end function parity_sign 

    end subroutine IF_hyperfine_structure

    subroutine construct_primitive_IF_basis(iverbose)
        ! |state, v; Lambda; S, Sigma; J, Omega; I1, Itot, F>
        implicit none

        integer(ik), INTENT(IN) :: iverbose
        REAL(rk) :: Ftot, J, Itot
        INTEGER(ik) :: index_Ftot, Ndimen_Ftot
        INTEGER(ik) :: index_Itot 
        INTEGER(ik) :: index_J, index_J_min, index_J_max, start_contr_Ftot, end_contr_Ftot, index_contr_Ftot
        TYPE(quantaT), POINTER :: contr_Ftot
        
        write(out, '(/A2, A)') '', "Start: contracted Hund's case (a_beta) basis"

        num_primitive_hfs_basis = 0
        do index_Ftot = 1, num_Ftot
            Ftot = Ftot_global_min + REAL(index_Ftot - 1, rk)

            if (iverbose>=4) write(out, '(A4, A, F5.1)') '', "F = ", Ftot

            Ftot_list(index_Ftot) = Ftot

            ! the 1st loop: count Ndimen_Ftot (sum over all Itot and J)
            Ndimen_Ftot = 0
            do index_Itot = 1, int(Itot_max) + 1
                Itot = Itot_list(index_Itot) ! 
                J_min = abs(Ftot - Itot)
                J_max = Ftot + Itot ! replace I1 as Itot

                index_J_min = nint(J_min - J_global_min) + 1
                index_J_max = nint(J_max - J_global_min) + 1

                do index_J = index_J_min, index_J_max
                        ! These J values must form a continuous (contiguous) sequence with the step of 1.0: 
                        if (basis(index_J)%Ndimen < 0 .or. basis(index_J)%Ndimen > 1000000) then
                            stop 'basis(index_J)%Ndimen corrupted. Please check if the input J values are contiguous.'
                        end if
                    Ndimen_Ftot = Ndimen_Ftot + basis(index_J)%Ndimen
                end do
            end do

            if (iverbose>=4) then
                write(out, '(A4, A, I20)') &
                    '', 'Number of the contracted basis', Ndimen_Ftot
            endif 

            primitive_hfs_basis(index_Ftot)%Ndimen = Ndimen_Ftot
            ALLOCATE (primitive_hfs_basis(index_Ftot)%icontr(Ndimen_Ftot))

            num_primitive_hfs_basis = &
                num_primitive_hfs_basis + Ndimen_Ftot
            
            ! the 2nd loop: fill icontr (append blocks in fixed order)
            end_contr_Ftot = 0
            do index_Itot = 1, int(Itot_max) + 1
                Itot = Itot_list(index_Itot)
                J_min = abs(Ftot - Itot)
                J_max = Ftot + Itot
                index_J_min = nint(J_min - J_global_min) + 1
                index_J_max = nint(J_max - J_global_min) + 1

                do index_J = index_J_min, index_J_max
                    J = J_global_min + REAL(index_J - 1, rk)
                    start_contr_Ftot = end_contr_Ftot
                    end_contr_Ftot = start_contr_Ftot + basis(index_J)%Ndimen
                    primitive_hfs_basis(index_Ftot)%icontr(start_contr_Ftot + 1:end_contr_Ftot) = basis(index_J)%icontr(:)
                    primitive_hfs_basis(index_Ftot)%icontr(start_contr_Ftot + 1:end_contr_Ftot)%Jrot = J
                    primitive_hfs_basis(index_Ftot)%icontr(start_contr_Ftot + 1:end_contr_Ftot)%Itot = Itot
                end do
            end do

            primitive_hfs_basis(index_Ftot)%icontr(:)%index_Ftot = index_Ftot
            primitive_hfs_basis(index_Ftot)%icontr(:)%I1 = I1
            primitive_hfs_basis(index_Ftot)%icontr(:)%Ftot = Ftot
            if (iverbose>=4) then
                write(out, '(A6, A7, A7, A7, A7, A7, A7, A7, A7, A7, A7)') &
                    '    ', 'No.', 'F', 'Itot', 'state', 'v', 'Lambda', 'J', 'Omega', 'S', 'Sigma'
            endif
            do index_contr_Ftot = 1, Ndimen_Ftot
                contr_Ftot => primitive_hfs_basis(index_Ftot)%icontr(index_contr_Ftot)
                if (iverbose>=4) then
                    write(out, '(A6, I7, F7.1, F7.1, I7, I7, I7, F7.1, F7.1, F7.1, F7.1)') &
                        '    ', index_contr_Ftot, contr_Ftot%Ftot, contr_Ftot%Itot, &
                        contr_Ftot%istate, contr_Ftot%v, contr_Ftot%ilambda, &
                        contr_Ftot%Jrot, contr_Ftot%omega, &
                        contr_Ftot%spin, contr_Ftot%sigma
                endif
            end do

            if (iverbose>=4) write(out, '(A4, A/)') '', '... done'
        end do
       
        write(out, '(A2, A/)') '', "End: contracted Hund's case (a_beta) basis"
    end subroutine construct_primitive_IF_basis

    subroutine construct_parity_conserved_IF_matrix( &
        index_Ftot, index_represC2v, index_represCs, ilabel_gu, Ndimen_Ftot, Nlevels_Ftot,&
        primitive_IF_hyperfine_matrix, parity_conserved_IF_matrix, transformation_matrix, eigen_value_rov_J)
        ! < phi_m^{tau, J}; I; F| H_hfs + H0 |phi_n^{tau, J'}; I; F>
        implicit none

        INTEGER(ik), INTENT(IN) :: index_Ftot, index_represC2v, index_represCs, ilabel_gu, &
                                   Ndimen_Ftot, Nlevels_Ftot
        REAL(rk), ALLOCATABLE, INTENT(IN) :: primitive_IF_hyperfine_matrix(:,:)
        REAL(rk), ALLOCATABLE, INTENT(INOUT) :: parity_conserved_IF_matrix(:,:)
        REAL(rk), ALLOCATABLE, INTENT(INOUT) :: transformation_matrix(:,:)
        REAL(rk), ALLOCATABLE, INTENT(INOUT), OPTIONAL :: eigen_value_rov_J(:)  
        ! eigen_value_rov_J(:) is an optional output to print the rovibronic eigenvalues,
        ! which can be combined to save memory in future 

        INTEGER(ik) :: index    
        REAL(rk), ALLOCATABLE :: eigen_value_J_in_Ftot(:)
        REAL(rk), ALLOCATABLE :: intermediate_matrix(:,:)
            
        ALLOCATE(eigen_value_J_in_Ftot(Nlevels_Ftot))
      
        CALL construct_IF_transformation_matrix( &
            index_Ftot, index_represC2v, index_represCs, ilabel_gu, &
            transformation_matrix, eigen_value_J_in_Ftot)

        ! ! H_hfs^tau = Phi^\dagger * H_hfs * Phi
        ! parity_conserved_IF_matrix = &
        !     matmul(matmul(transpose(transformation_matrix), &
        !                   primitive_IF_hyperfine_matrix), &
        !            transformation_matrix)

        ALLOCATE(intermediate_matrix(Nlevels_Ftot, Ndimen_Ftot))

        ! M_intermediate = Phi^\dagger * H_hfs
        CALL dgemm('T', 'N', Nlevels_Ftot, Ndimen_Ftot, Ndimen_Ftot, 1.0_rk, &
                    transformation_matrix, Ndimen_Ftot, &
                    primitive_IF_hyperfine_matrix, Ndimen_Ftot, &
                    0.0_rk, intermediate_matrix, Nlevels_Ftot)

        ! H_hfs^tau = M_intermediate * Phi           
        CALL dgemm('N', 'N', Nlevels_Ftot, Nlevels_Ftot, Ndimen_Ftot, 1.0_rk, &
                    intermediate_matrix, Nlevels_Ftot, &
                    transformation_matrix, Ndimen_Ftot, &
                    0.0_rk, parity_conserved_IF_matrix, Nlevels_Ftot)

        DEALLOCATE(intermediate_matrix)

        ! H^tau = H0^tau + H_hfs^tau 
        do index = 1, Nlevels_Ftot
            parity_conserved_IF_matrix(index, index) = &
                parity_conserved_IF_matrix(index, index) &
                + eigen_value_J_in_Ftot(index)
        end do  

        if (present(eigen_value_rov_J)) then
            eigen_value_rov_J = eigen_value_J_in_Ftot
        end if 

        DEALLOCATE(eigen_value_J_in_Ftot)
    end subroutine construct_parity_conserved_IF_matrix

    subroutine construct_IF_transformation_matrix( &
        index_Ftot, index_represC2v, index_represCs, ilabel_gu, &
        transformation_matrix, &
        eigen_value_J_in_Ftot)
        ! Construct the basis transformation matrix Phi between 
        ! the primitive IF basis set and the selected rovibronic symmetry block (+g/-u/-g/+u).

        implicit none

        INTEGER(ik), INTENT(IN) :: index_Ftot, index_represC2v, index_represCs, ilabel_gu
        REAL(rk), ALLOCATABLE, INTENT(INOUT) :: transformation_matrix(:,:), &
                                                eigen_value_J_in_Ftot(:)
        REAL(rk) :: Ftot, J, Itot
        integer(ik) :: index_Itot, ilevel 
        INTEGER(ik) ::  index_J, index_J_min, index_J_max, &
                        start_dimen_Ftot, end_dimen_Ftot, end_level_Ftot
        
        transformation_matrix = 0.0_rk
        eigen_value_J_in_Ftot = 0.0_rk

        Ftot = Ftot_list(index_Ftot)

        end_dimen_Ftot = 0
        end_level_Ftot = 0

        do index_Itot = 1, size(Itot_list)

            Itot = Itot_list(index_Itot)
            J_min = abs(Ftot - Itot)
            J_max = Ftot + Itot

            index_J_min = nint(J_min - J_global_min) + 1
            index_J_max = nint(J_max - J_global_min) + 1

            do index_J = index_J_min, index_J_max

                J = J_global_min + REAL(index_J - 1, rk)

                start_dimen_Ftot = end_dimen_Ftot

                ! rows always follow the primitive IF basis
                end_dimen_Ftot = start_dimen_Ftot + basis(index_J)%Ndimen

                ! physical Itot filter
                if (index_represC2v <= 2) then  ! A1/A2
                    if (mod(nint(Itot), 2) /= 0) cycle   ! Itot = even 
                else ! B1/B2
                    if (mod(nint(Itot), 2) /= 1) cycle   ! Itot = odd
                end if

                if (eigen(index_J, index_represCs)%Nlevels == 0) cycle
                if (eigen(index_J, index_represCs)%Ndimen == 0) then
                    write(out,*) 'ERROR: Nlevels > 0 but Ndimen = 0 in construct_IF_transformation_matrix'
                    stop
                end if

                do ilevel = 1, eigen(index_J, index_represCs)%Nlevels

                    if (poten(eigen(index_J, index_represCs)%quanta(ilevel)%istate)%parity%gu /= ilabel_gu) cycle

                    end_level_Ftot = end_level_Ftot + 1

                    ! Eigenvalues of H0
                    eigen_value_J_in_Ftot(end_level_Ftot) = eigen(index_J, index_represCs)%val(ilevel)

                    ! Phi
                    transformation_matrix(start_dimen_Ftot + 1 : end_dimen_Ftot, end_level_Ftot) = &
                        eigen(index_J, index_represCs)%vect(:, ilevel)
                    
                    ! eigen_all_hfs is the table of IF eigen states for all F values.
                    eigen_all_hfs(index_Ftot, index_represC2v)%quanta(end_level_Ftot) = &
                        eigen(index_J, index_represCs)%quanta(ilevel)
                    eigen_all_hfs(index_Ftot, index_represC2v)%quanta(end_level_Ftot)%Ftot       = Ftot
                    eigen_all_hfs(index_Ftot, index_represC2v)%quanta(end_level_Ftot)%I1         = I1
                    eigen_all_hfs(index_Ftot, index_represC2v)%quanta(end_level_Ftot)%index_Ftot = index_Ftot
                    eigen_all_hfs(index_Ftot, index_represC2v)%quanta(end_level_Ftot)%Jrot       = J
                    eigen_all_hfs(index_Ftot, index_represC2v)%quanta(end_level_Ftot)%Itot       = Itot

                end do
            end do
        end do

        ! Check that the Phi transformation matrix has been filled consistently: 
        if (end_dimen_Ftot /= size(transformation_matrix,1)) then
            stop 'Phi row dimension mismatch.'
        end if

        if (end_level_Ftot /= size(transformation_matrix,2)) then
            stop 'Phi column dimension mismatch.'
        end if

    end subroutine construct_IF_transformation_matrix

    subroutine construct_primitive_IF_hyperfine_matrix( &
        index_Ftot, primitive_IF_hyperfine_matrix)
        ! <b, J; I; F| H_hfs |k, J'; I; F>
        implicit none

        INTEGER(ik), INTENT(IN) :: index_Ftot
        REAL(rk), ALLOCATABLE, INTENT(INOUT) :: primitive_IF_hyperfine_matrix(:, :)

        REAL(rk) :: Ftot, Itot_bra, Itot_ket, &
                    S_bra, S_ket, Sigma_bra, Sigma_ket, &
                    J_bra, J_ket, Omega_bra, Omega_ket
                    
        INTEGER(ik) :: state_bra, state_ket, bra, ket, v_bra, v_ket, &
                       index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        INTEGER(ik) :: Ndimen_Ftot, index_hfcc1, index_hfcc12, index_field

        primitive_IF_hyperfine_matrix = 0

        Ndimen_Ftot = primitive_hfs_basis(index_Ftot)%Ndimen
        Ftot = Ftot_list(index_Ftot)

        do ket = 1, Ndimen_Ftot
            call get_quanta(index_Ftot, ket, &
                            Ftot, state_ket, index_v_ket, v_ket, Lambda_ket, &
                            S_ket, Sigma_ket, J_ket, Omega_ket, Itot=Itot_ket)

            do bra = 1, Ndimen_Ftot
                call get_quanta(index_Ftot, bra, &
                                Ftot, state_bra, index_v_bra, v_bra, Lambda_bra, &
                                S_bra, Sigma_bra, J_bra, Omega_bra, Itot=Itot_bra)

                if (state_bra == state_ket) then

                    do index_hfcc1 = 1, GLOBAL_NUM_HFCC_OBJECT
                        do index_field = 1, hfcc1(index_hfcc1)%num_field
                            if ( (state_bra == hfcc1(index_hfcc1)%field(index_field)%istate) &
                                .and. (state_ket == hfcc1(index_hfcc1)%field(index_field)%jstate)) then

                                primitive_IF_hyperfine_matrix(bra, ket) = &
                                    primitive_IF_hyperfine_matrix(bra, ket) &
                                    + primitive_IF_hyperfine_matrix_element( &
                                        Ftot, Itot_bra, Itot_ket, index_hfcc1, index_field, &
                                        index_v_bra, index_v_ket, &
                                        Lambda_bra, Lambda_ket, &
                                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                        J_bra, J_ket, Omega_bra, Omega_ket)
                            end if
                        end do                  
                    end do

                    do index_hfcc12 = 1, GLOBAL_NUM_HFCC12_OBJECT
                        do index_field = 1, hfcc12(index_hfcc12)%num_field
                            if ( (state_bra == hfcc12(index_hfcc12)%field(index_field)%istate) &
                                .and. (state_ket == hfcc12(index_hfcc12)%field(index_field)%jstate)) then

                                select case (index_hfcc12)
                                case (1)
                                    primitive_IF_hyperfine_matrix(bra, ket) = &
                                        primitive_IF_hyperfine_matrix(bra, ket) &
                                        + Heq_II_element( &
                                            Ftot, Itot_bra, Itot_ket, index_field, &
                                            index_v_bra, index_v_ket, &
                                            Lambda_bra, Lambda_ket, &
                                            S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                            J_bra, J_ket, Omega_bra, Omega_ket)
                                case (2)
                                    primitive_IF_hyperfine_matrix(bra, ket) = &
                                        primitive_IF_hyperfine_matrix(bra, ket) &
                                        + Heq_IIS_element( &
                                            Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket, &
                                            Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                            J_bra, J_ket, Omega_bra, Omega_ket)
                                end select
                            end if
                        end do
                    end do

                end if
            end do
        end do
    end subroutine construct_primitive_IF_hyperfine_matrix

    function primitive_IF_hyperfine_matrix_element( &
        Ftot, Itot_bra, Itot_ket, index_hfcc1, index_field, &
        index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_hfs)

        implicit none

        REAL(rk), INTENT(IN) :: Ftot, Itot_bra, Itot_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        INTEGER(ik), INTENT(IN) :: index_hfcc1, index_field, &
                                   index_v_bra, index_v_ket, &
                                   Lambda_bra, Lambda_ket
        REAL(rk) :: H_hfs

        ! The switch flow could be refactored if 
        ! an abstract interface is defined.
        ! This function could be removed then.
        H_hfs = 0.0_rk
        if (index_hfcc1 == 1) then 
            H_hfs = Heq_FC_element( &
                        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket) 
        else if (index_hfcc1 == 2) then
            H_hfs = Heq_IL_element( &
                        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket)
        else if (index_hfcc1 == 3) then 
            H_hfs = Heq_DIP_C_element( &
                        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket)  
        else if (index_hfcc1 == 4) then 
            H_hfs = Heq_DIP_D_element( &
                        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket)                 
        else if (index_hfcc1 == 5) then 
            H_hfs = Heq_IJ_element( &
                        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket) 
        else if (index_hfcc1 == 6) then
            H_hfs = Heq_CIp_element( &
                        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket)
        else if (index_hfcc1 == 7) then
            H_hfs = Heq_EQQ0_element( &        
                        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket,  &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket)
        else if (index_hfcc1 == 8) then
            H_hfs = Heq_EQQ2_element( &        
                        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket,  &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket)
        end if
    end function primitive_IF_hyperfine_matrix_element

    subroutine F1F_hyperfine_structure(iverbose)
        ! Top level for the heteronuclear two-spin coupling scheme:
        ! |state, v; Lambda; S, Sigma; J, Omega; I1, F1; I2, Ftot>
        implicit none
        integer(ik), INTENT(IN) :: iverbose
        INTEGER(ik) :: index_Ftot, index_represCs, iroot, iF_ID
        INTEGER(ik) :: index_F1, index_J, index_J_min, index_J_max
        INTEGER(ik) :: Ndimen_Ftot, Nlevels_Ftot, index_levels_Ftot, maxlocation, maxlocation2
        REAL(rk) :: Ftot, F1, F1_min, F1_max
        REAL(rk), ALLOCATABLE :: eigen_value_Ftot(:), &
                                    parity_conserved_F1F_matrix(:, :), &
                                    primitive_F1F_hyperfine_matrix(:,:), &
                                    transformation_matrix(:,:), &
                                    eigen_vector_Ftot(:, :), &
                                    eigen_value_rov_J(:)
        TYPE(quantaT), POINTER :: quanta

        I1 = hyperfine_setup%I1
        I2 = hyperfine_setup%I2

        J_global_max = jmax
        if (J_global_max < I1 + I2) then
            write (out, *) "J_global_max < I1 + I2"
            stop
        end if

        Ftot_global_max = J_global_max - I1 - I2

        if (mod(nint(Ftot_global_max*2.0_rk), 2) == 0) then
            Ftot_global_min = 0.0_rk
        else
            Ftot_global_min = 0.5_rk
        end if

        ! The basis of J always starts from the minimum of J
        if (mod(nint(J_global_max*2.0_rk), 2) == 0) then
            J_global_min = 0.0_rk
        else
            J_global_min = 0.5_rk
        end if

        num_Ftot = nint(Ftot_global_max - Ftot_global_min) + 1

        CALL construct_hyperfine_constant_field_matrix(iverbose)

        ALLOCATE(Ftot_list(num_Ftot))
        ALLOCATE(primitive_hfs_basis(num_Ftot))

        CALL construct_primitive_F1F_basis(iverbose)

        if (iverbose>=4) then
            write(out, '(/A2, A)') '', 'Start: hyperfine states'
            write(out, '(/A2, A, F22.12)') '', 'Zero point energy =', job%ZPE
        endif

        ALLOCATE(eigen_all_hfs(num_Ftot, sym%NrepresCs))
        eigen_all_hfs(:,:)%Nlevels = 0

        iroot = 0
        min_eigen_value_hfs = safe_max
        do index_Ftot = 1, num_Ftot
            if (iverbose>=4) then
                write(out, '(/A4, A, F6.1)') '', 'Eigen states of Ftot =', Ftot_list(index_Ftot)
                write(out, '(/A6, A)') '', 'Construct primitive hyperfine matrix'
            endif

            Ndimen_Ftot = primitive_hfs_basis(index_Ftot)%Ndimen
            ALLOCATE(primitive_F1F_hyperfine_matrix(Ndimen_Ftot, Ndimen_Ftot))
            CALL construct_primitive_F1F_hyperfine_matrix(index_Ftot, primitive_F1F_hyperfine_matrix)

            if (iverbose>=4) write(out, '(A6, A/)') '', '... done'

            do index_represCs = 1, sym%NrepresCs
                if (iverbose>=4) then
                    write(out, '(/A6, A, A3)') '', 'Eigen states of parity', trim(parity_sign(index_represCs))
                endif

                Nlevels_Ftot = 0
                Ftot = Ftot_list(index_Ftot)
                F1_min = abs(Ftot - I2)
                F1_max = Ftot + I2

                do index_F1 = 1, nint(F1_max - F1_min) + 1
                    F1 = F1_min + real(index_F1 - 1, rk)
                    J_min = abs(F1 - I1)
                    J_max = F1 + I1
                    index_J_min = nint(J_min - J_global_min) + 1
                    index_J_max = nint(J_max - J_global_min) + 1

                    do index_J = index_J_min, index_J_max
                        if (eigen(index_J, index_represCs)%Ndimen > 0) then
                            Nlevels_Ftot = Nlevels_Ftot + eigen(index_J,index_represCs)%Nlevels
                        else
                            write(out,'(A8, A, F5.1, A, A, A)') '', 'The current block with (J = ', &
                                J_global_min + REAL(index_J - 1, rk), &
                                ', parity = ', trim(parity_sign(index_represCs)), ') is empty.'
                        end if
                    end do
                end do

                if (Nlevels_Ftot == 0) then
                    write(out,'(/A8, A, F5.2, A, A, A)') '', 'The current block with (Ftot = ', &
                        Ftot_list(index_Ftot), &
                        ', parity = ', trim(parity_sign(index_represCs)), ') is empty.'
                    cycle
                end if

                ALLOCATE(eigen_value_Ftot(Nlevels_Ftot))
                ALLOCATE(eigen_value_rov_J(Nlevels_Ftot))
                ALLOCATE(parity_conserved_F1F_matrix(Nlevels_Ftot, Nlevels_Ftot))
                ALLOCATE(transformation_matrix(Ndimen_Ftot, Nlevels_Ftot))
                ALLOCATE(eigen_vector_Ftot(Ndimen_Ftot, Nlevels_Ftot))

                ALLOCATE(eigen_all_hfs(index_Ftot, index_represCs)%quanta(Nlevels_Ftot))
                ALLOCATE(eigen_all_hfs(index_Ftot, index_represCs)%val(Nlevels_Ftot))
                ALLOCATE(eigen_all_hfs(index_Ftot, index_represCs)%vect(Ndimen_Ftot, Nlevels_Ftot))

                if (iverbose>=4) then
                    write(out, '(/A8, A)') &
                        '', 'Construct parity conserved hyperfine matrix'
                endif

                CALL construct_parity_conserved_F1F_matrix( &
                    index_Ftot, index_represCs, Ndimen_Ftot, Nlevels_Ftot, primitive_F1F_hyperfine_matrix, &
                    parity_conserved_F1F_matrix, transformation_matrix, eigen_value_rov_J)
                if (iverbose>=4) write(out, '(A8, A/)') '', '... done'

                if (iverbose>=4) then
                    write(out, '(/A8, A)') &
                        '', 'Diagonalize parity conserved hyperfine matrix'
                endif
                CALL lapack_syev(parity_conserved_F1F_matrix, eigen_value_Ftot)
                if (iverbose>=4) write(out, '(A8, A/)') '', '... done'

                if (iverbose>=4) then
                    write(out, '(/A8, A)') &
                         '', "Transform wavefunction to Hund's case (a_beta) basis set"
                endif
                CALL dgemm('N', 'N', Ndimen_Ftot, Nlevels_Ftot, Nlevels_Ftot, 1.0_rk, &
                            transformation_matrix, Ndimen_Ftot, &
                            parity_conserved_F1F_matrix, Nlevels_Ftot, &
                            0.0_rk, eigen_vector_Ftot, Ndimen_Ftot)
                if (iverbose>=4) write(out, '(A8, A/)') '', '... done'

                eigen_all_hfs(index_Ftot, index_represCs)%Nlevels = Nlevels_Ftot
                eigen_all_hfs(index_Ftot, index_represCs)%Ndimen = Ndimen_Ftot
                eigen_all_hfs(index_Ftot, index_represCs)%val = eigen_value_Ftot
                eigen_all_hfs(index_Ftot, index_represCs)%vect = eigen_vector_Ftot

                if (iverbose>=4) then
                    write(out, '(/A8, A)') '', 'Calculate energy levels'
                endif
                eigen_value_Ftot(:) = eigen_value_Ftot(:) - job%ZPE

                iF_ID = 0
                do index_levels_Ftot = 1, Nlevels_Ftot
                    maxlocation = MAXLOC(ABS(eigen_vector_Ftot(:, index_levels_Ftot)), DIM=1)
                    maxlocation2 = MAXLOC(ABS(parity_conserved_F1F_matrix(:, index_levels_Ftot)), DIM=1)
                    quanta => primitive_hfs_basis(index_Ftot)%icontr(maxlocation)

                    iroot = iroot + 1
                    iF_ID = iF_ID + 1
                    eigen_all_hfs(index_Ftot, index_represCs)%quanta(index_levels_Ftot)%iroot = iroot
                    eigen_all_hfs(index_Ftot, index_represCs)%quanta(index_levels_Ftot)%iF_ID = iF_ID

                    min_eigen_value_hfs = min(&
                        eigen_all_hfs(index_Ftot, index_represCs)%val(index_levels_Ftot), &
                        min_eigen_value_hfs)

                    if (iverbose>=4) then
                        write(unit_hyperfine_states, '(I7, 3F22.12, I7, A7, A3, A7, A7, A12, I7, I4, 2A5)') &
                            iroot, &
                            eigen_value_Ftot(index_levels_Ftot), &
                            eigen_value_rov_J(maxlocation2) - job%ZPE, &
                            eigen_value_Ftot(index_levels_Ftot) &
                                - (eigen_value_rov_J(maxlocation2) - job%ZPE), &
                            int(quanta%Ftot * 2.0_rk + 1.0_rk), &
                            qnum2a(quanta%Ftot, 7), &
                            trim(parity_sign(index_represCs)), &
                            qnum2a(quanta%F1, 7), &
                            qnum2a(quanta%Jrot, 7), &
                            trim(poten(quanta%istate)%name), &
                            eigen_all_hfs(index_Ftot, index_represCs)%quanta(maxlocation2)%v, &
                            quanta%ilambda, &
                            qnum2a(quanta%sigma, 5), &
                            qnum2a(quanta%omega, 5)
                    endif
                enddo
                if (iverbose>=4) write(out, '(A8, A/)') '', '... done'

                DEALLOCATE(eigen_vector_Ftot)
                DEALLOCATE(transformation_matrix)
                DEALLOCATE(parity_conserved_F1F_matrix)
                DEALLOCATE(eigen_value_Ftot)
                DEALLOCATE(eigen_value_rov_J)
            end do

            DEALLOCATE(primitive_F1F_hyperfine_matrix)
        end do

        if (iverbose>=4) write(out, '(A2, A/)') '', 'End: hyperfine states'
        write(out, '(A/)') "End: hyperfine energy calculation"

        close(unit=unit_hyperfine_states)
    contains
        function parity_sign(index_represCs)
            INTEGER(ik), INTENT(IN) :: index_represCs
            CHARACTER(4) :: parity_sign
            if ( index_represCs == 1 ) then
                parity_sign = '+'
            elseif ( index_represCs == 2 ) then
                parity_sign = '-'
            endif
        end function parity_sign
    end subroutine F1F_hyperfine_structure

    subroutine construct_primitive_F1F_basis(iverbose)
        ! |state, v; Lambda; S, Sigma; J, Omega; I1, F1; I2, Ftot>
        implicit none

        integer(ik), INTENT(IN) :: iverbose
        REAL(rk) :: Ftot, F1, F1_min, F1_max, J
        INTEGER(ik) :: index_Ftot, index_F1, Ndimen_Ftot
        INTEGER(ik) :: index_J, index_J_min, index_J_max, start_contr_Ftot, end_contr_Ftot, index_contr_Ftot
        TYPE(quantaT), POINTER :: contr_Ftot

        write(out, '(/A2, A)') '', "Start: contracted Hund's case (a_beta) basis"

        num_primitive_hfs_basis = 0
        do index_Ftot = 1, num_Ftot
            Ftot = Ftot_global_min + REAL(index_Ftot - 1, rk)

            if (iverbose>=4) write(out, '(A4, A, F5.1)') '', "Ftot = ", Ftot

            Ftot_list(index_Ftot) = Ftot
            F1_min = abs(Ftot - I2)
            F1_max = Ftot + I2

            Ndimen_Ftot = 0
            do index_F1 = 1, nint(F1_max - F1_min) + 1
                F1 = F1_min + real(index_F1 - 1, rk)
                J_min = abs(F1 - I1)
                J_max = F1 + I1

                index_J_min = nint(J_min - J_global_min) + 1
                index_J_max = nint(J_max - J_global_min) + 1

                do index_J = index_J_min, index_J_max
                    Ndimen_Ftot = Ndimen_Ftot + basis(index_J)%Ndimen
                end do
            end do

            if (iverbose>=4) then
                write(out, '(A4, A, I5)') &
                    '', 'Number of the contracted basis', Ndimen_Ftot
            endif

            primitive_hfs_basis(index_Ftot)%Ndimen = Ndimen_Ftot
            ALLOCATE (primitive_hfs_basis(index_Ftot)%icontr(Ndimen_Ftot))

            num_primitive_hfs_basis = &
                num_primitive_hfs_basis + Ndimen_Ftot

            end_contr_Ftot = 0
            do index_F1 = 1, nint(F1_max - F1_min) + 1
                F1 = F1_min + real(index_F1 - 1, rk)
                J_min = abs(F1 - I1)
                J_max = F1 + I1

                index_J_min = nint(J_min - J_global_min) + 1
                index_J_max = nint(J_max - J_global_min) + 1

                do index_J = index_J_min, index_J_max
                    J = J_global_min + REAL(index_J - 1, rk)
                    start_contr_Ftot = end_contr_Ftot
                    end_contr_Ftot = start_contr_Ftot + basis(index_J)%Ndimen
                    primitive_hfs_basis(index_Ftot)%icontr(start_contr_Ftot + 1:end_contr_Ftot) = &
                        basis(index_J)%icontr(:)
                    primitive_hfs_basis(index_Ftot)%icontr(start_contr_Ftot + 1:end_contr_Ftot)%Jrot = J
                    primitive_hfs_basis(index_Ftot)%icontr(start_contr_Ftot + 1:end_contr_Ftot)%F1 = F1
                end do
            end do

            primitive_hfs_basis(index_Ftot)%icontr(:)%index_Ftot = index_Ftot
            primitive_hfs_basis(index_Ftot)%icontr(:)%I1 = I1
            primitive_hfs_basis(index_Ftot)%icontr(:)%I2 = I2
            primitive_hfs_basis(index_Ftot)%icontr(:)%Ftot = Ftot
            if (iverbose>=4) then
                write(out, '(A6, A7, A7, A7, A7, A7, A7, A7, A7, A7, A7, A7)') &
                    '    ', 'No.', 'Ftot', 'F1', 'state', 'v', 'Lambda', 'J', 'Omega', 'S', 'Sigma', 'I2'
            endif
            do index_contr_Ftot = 1, Ndimen_Ftot
                contr_Ftot => primitive_hfs_basis(index_Ftot)%icontr(index_contr_Ftot)
                if (iverbose>=4) then
                    write(out, '(A6, I7, F7.1, F7.1, I7, I7, I7, F7.1, F7.1, F7.1, F7.1, F7.1)') &
                        '    ', index_contr_Ftot, contr_Ftot%Ftot, contr_Ftot%F1, &
                        contr_Ftot%istate, contr_Ftot%v, contr_Ftot%ilambda, &
                        contr_Ftot%Jrot, contr_Ftot%omega, &
                        contr_Ftot%spin, contr_Ftot%sigma, contr_Ftot%I2
                endif
            end do

            if (iverbose>=4) write(out, '(A4, A/)') '', '... done'
        end do

        write(out, '(A2, A/)') '', "End: contracted Hund's case (a_beta) basis"
    end subroutine construct_primitive_F1F_basis

    subroutine construct_parity_conserved_F1F_matrix( &
        index_Ftot, index_represCs, Ndimen_Ftot, Nlevels_Ftot,&
        primitive_F1F_hyperfine_matrix, parity_conserved_F1F_matrix, transformation_matrix, eigen_value_rov_J)
        ! < phi_m^{tau, J}; I1, F1; I2, Ftot | H_hfs + H0 | phi_n^{tau, J'}; I1, F1'; I2, Ftot>
        implicit none

        INTEGER(ik), INTENT(IN) :: index_Ftot, index_represCs
        REAL(rk), ALLOCATABLE, INTENT(IN) :: primitive_F1F_hyperfine_matrix(:,:)
        REAL(rk), ALLOCATABLE, INTENT(INOUT) :: parity_conserved_F1F_matrix(:,:)
        REAL(rk), ALLOCATABLE, INTENT(INOUT) :: transformation_matrix(:,:)
        REAL(rk), ALLOCATABLE, INTENT(INOUT), OPTIONAL :: eigen_value_rov_J(:)

        INTEGER(ik) :: index, Ndimen_Ftot, Nlevels_Ftot
        REAL(rk), ALLOCATABLE :: eigen_value_J_in_Ftot(:)
        REAL(rk), ALLOCATABLE :: intermediate_matrix(:,:)

        ALLOCATE(eigen_value_J_in_Ftot(Nlevels_Ftot))

        CALL construct_F1F_transformation_matrix( &
            index_Ftot, index_represCs, &
            transformation_matrix, eigen_value_J_in_Ftot)

        ALLOCATE(intermediate_matrix(Nlevels_Ftot, Ndimen_Ftot))

        CALL dgemm('T', 'N', Nlevels_Ftot, Ndimen_Ftot, Ndimen_Ftot, 1.0_rk, &
                    transformation_matrix, Ndimen_Ftot, &
                    primitive_F1F_hyperfine_matrix, Ndimen_Ftot, &
                    0.0_rk, intermediate_matrix, Nlevels_Ftot)

        CALL dgemm('N', 'N', Nlevels_Ftot, Nlevels_Ftot, Ndimen_Ftot, 1.0_rk, &
                    intermediate_matrix, Nlevels_Ftot, &
                    transformation_matrix, Ndimen_Ftot, &
                    0.0_rk, parity_conserved_F1F_matrix, Nlevels_Ftot)

        DEALLOCATE(intermediate_matrix)

        do index = 1, Nlevels_Ftot
            parity_conserved_F1F_matrix(index, index) = &
                parity_conserved_F1F_matrix(index, index) &
                + eigen_value_J_in_Ftot(index)
        end do

        if (present(eigen_value_rov_J)) then
            eigen_value_rov_J = eigen_value_J_in_Ftot
        end if

        DEALLOCATE(eigen_value_J_in_Ftot)
    end subroutine construct_parity_conserved_F1F_matrix

    subroutine construct_F1F_transformation_matrix( &
        index_Ftot, index_represCs, &
        transformation_matrix, &
        eigen_value_J_in_Ftot)
        ! Construct Phi between the primitive and parity-conserved basis sets.
        implicit none

        INTEGER(ik), INTENT(IN) :: index_Ftot, index_represCs
        REAL(rk), ALLOCATABLE, INTENT(INOUT) :: transformation_matrix(:,:), &
                                                eigen_value_J_in_Ftot(:)

        REAL(rk) :: Ftot, F1, F1_min, F1_max, J
        INTEGER(ik) :: index_F1, index_J, index_J_min, index_J_max, &
                        start_dimen_Ftot, end_dimen_Ftot, &
                        start_level_Ftot, end_level_Ftot

        transformation_matrix = 0.0_rk
        eigen_value_J_in_Ftot = 0.0_rk

        Ftot = Ftot_global_min + REAL(index_Ftot - 1, rk)
        F1_min = abs(Ftot - I2)
        F1_max = Ftot + I2

        end_dimen_Ftot = 0
        end_level_Ftot = 0
        do index_F1 = 1, nint(F1_max - F1_min) + 1
            F1 = F1_min + real(index_F1 - 1, rk)
            J_min = abs(F1 - I1)
            J_max = F1 + I1

            index_J_min = nint(J_min - J_global_min) + 1
            index_J_max = nint(J_max - J_global_min) + 1

            do index_J = index_J_min, index_J_max

                J = J_global_min + REAL(index_J - 1, rk)

                start_dimen_Ftot = end_dimen_Ftot
                start_level_Ftot = end_level_Ftot
                end_dimen_Ftot = start_dimen_Ftot + basis(index_J)%Ndimen

                if (eigen(index_J, index_represCs)%Nlevels == 0) cycle

                end_level_Ftot = start_level_Ftot + eigen(index_J, index_represCs)%Nlevels

                eigen_value_J_in_Ftot(start_level_Ftot + 1 : end_level_Ftot) &
                    = eigen(index_J, index_represCs)%val

                transformation_matrix(start_dimen_Ftot + 1 : end_dimen_Ftot, &
                                      start_level_Ftot + 1 : end_level_Ftot) &
                    = eigen(index_J, index_represCs)%vect

                eigen_all_hfs(index_Ftot, index_represCs)% &
                        quanta(start_level_Ftot + 1 : end_level_Ftot) &
                    = eigen(index_J, index_represCs)%quanta
                eigen_all_hfs(index_Ftot, index_represCs)% &
                        quanta(start_level_Ftot + 1 : end_level_Ftot)%Ftot &
                    = Ftot
                eigen_all_hfs(index_Ftot, index_represCs)% &
                        quanta(start_level_Ftot + 1 : end_level_Ftot)%F1 &
                    = F1
                eigen_all_hfs(index_Ftot, index_represCs)% &
                        quanta(start_level_Ftot + 1 : end_level_Ftot)%I1 &
                    = I1
                eigen_all_hfs(index_Ftot, index_represCs)% &
                        quanta(start_level_Ftot + 1 : end_level_Ftot)%I2 &
                    = I2
                eigen_all_hfs(index_Ftot, index_represCs)% &
                        quanta(start_level_Ftot + 1 : end_level_Ftot)%index_Ftot &
                    = index_Ftot
                eigen_all_hfs(index_Ftot, index_represCs)% &
                        quanta(start_level_Ftot + 1 : end_level_Ftot)%Jrot &
                    = J
            end do
        end do
    end subroutine construct_F1F_transformation_matrix

    subroutine construct_primitive_F1F_hyperfine_matrix( &
        index_Ftot, primitive_F1F_hyperfine_matrix)
        ! <b, J, F1; I2, Ftot | H_hfs | k, J', F1'; I2, Ftot>
        implicit none

        INTEGER(ik), INTENT(IN) :: index_Ftot
        REAL(rk), ALLOCATABLE, INTENT(INOUT) :: primitive_F1F_hyperfine_matrix(:, :)

        REAL(rk) :: Ftot, F1_bra, F1_ket, &
                    S_bra, S_ket, Sigma_bra, Sigma_ket, &
                    J_bra, J_ket, Omega_bra, Omega_ket
        INTEGER(ik) :: state_bra, state_ket, bra, ket, v_bra, v_ket, &
                       index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        INTEGER(ik) :: Ndimen_Ftot, index_hfcc, index_hfcc12, index_field

        primitive_F1F_hyperfine_matrix = 0.0_rk

        Ndimen_Ftot = primitive_hfs_basis(index_Ftot)%Ndimen

        do ket = 1, Ndimen_Ftot
            call get_quanta(index_Ftot, ket, &
                            Ftot, state_ket, index_v_ket, v_ket, Lambda_ket, &
                            S_ket, Sigma_ket, J_ket, Omega_ket, F1=F1_ket)

            do bra = 1, Ndimen_Ftot
                call get_quanta(index_Ftot, bra, &
                                Ftot, state_bra, index_v_bra, v_bra, Lambda_bra, &
                                S_bra, Sigma_bra, J_bra, Omega_bra, F1=F1_bra)

                if (state_bra == state_ket) then

                    do index_hfcc = 1, GLOBAL_NUM_HFCC_OBJECT
                        do index_field = 1, hfcc1(index_hfcc)%num_field
                            if ( (state_bra == hfcc1(index_hfcc)%field(index_field)%istate) &
                                .and. (state_ket == hfcc1(index_hfcc)%field(index_field)%jstate)) then

                                if (nint(F1_bra - F1_ket) == 0) then
                                    primitive_F1F_hyperfine_matrix(bra, ket) = &
                                        primitive_F1F_hyperfine_matrix(bra, ket) &
                                        + primitive_F1_hyperfine_matrix_element( &
                                            F1_bra, index_hfcc, index_field, &
                                            index_v_bra, index_v_ket, &
                                            Lambda_bra, Lambda_ket, &
                                            S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                            J_bra, J_ket, Omega_bra, Omega_ket)
                                end if
                            end if
                        end do

                        do index_field = 1, hfcc2(index_hfcc)%num_field
                            if ( (state_bra == hfcc2(index_hfcc)%field(index_field)%istate) &
                                .and. (state_ket == hfcc2(index_hfcc)%field(index_field)%jstate)) then

                                primitive_F1F_hyperfine_matrix(bra, ket) = &
                                    primitive_F1F_hyperfine_matrix(bra, ket) &
                                    + primitive_F1F_hyperfine_matrix_element_h2( &
                                        Ftot, index_hfcc, index_field, index_v_bra, index_v_ket, &
                                        Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                        J_bra, J_ket, Omega_bra, Omega_ket, F1_bra, F1_ket)
                            end if
                        end do
                    end do

                    do index_hfcc12 = 1, GLOBAL_NUM_HFCC12_OBJECT
                        do index_field = 1, hfcc12(index_hfcc12)%num_field
                            if ( (state_bra == hfcc12(index_hfcc12)%field(index_field)%istate) &
                                .and. (state_ket == hfcc12(index_hfcc12)%field(index_field)%jstate)) then

                                select case (index_hfcc12)
                                case (1)
                                    primitive_F1F_hyperfine_matrix(bra, ket) = &
                                        primitive_F1F_hyperfine_matrix(bra, ket) &
                                        + H12_II_element( &
                                            Ftot, F1_bra, F1_ket, index_field, &
                                            index_v_bra, index_v_ket, &
                                            Lambda_bra, Lambda_ket, &
                                            S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                            J_bra, J_ket, Omega_bra, Omega_ket)
                                case (2)
                                    primitive_F1F_hyperfine_matrix(bra, ket) = &
                                        primitive_F1F_hyperfine_matrix(bra, ket) &
                                        + H12_IIS_element( &
                                            Ftot, F1_bra, F1_ket, index_field, index_v_bra, index_v_ket, &
                                            Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                            J_bra, J_ket, Omega_bra, Omega_ket)
                                end select
                            end if
                        end do
                    end do

                end if
            end do
        end do
    end subroutine construct_primitive_F1F_hyperfine_matrix

    function primitive_F1F_hyperfine_matrix_element_h2( &
        Ftot, index_hfcc2, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket, F1_bra, F1_ket) &
        result(H_hfs)

        implicit none

        REAL(rk), INTENT(IN) :: Ftot, F1_bra, F1_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        INTEGER(ik), INTENT(IN) :: index_hfcc2, index_field, &
                                   index_v_bra, index_v_ket, &
                                   Lambda_bra, Lambda_ket
        REAL(rk) :: H_hfs

        H_hfs = 0.0_rk
        if (index_hfcc2 >= 1 .and. index_hfcc2 <= 4) then
            H_hfs = H2_MD_element( &
                        Ftot, index_hfcc2, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket, F1_bra, F1_ket)
        else if (index_hfcc2 == 5) then
            H_hfs = H2_IJ_element( &
                        Ftot, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket, F1_bra, F1_ket)
        else if (index_hfcc2 == 6) then
            H_hfs = H2_CIp_element( &
                        Ftot, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket, F1_bra, F1_ket)
        else if (index_hfcc2 >= 7 .and. index_hfcc2 <= 8) then
            H_hfs = H2_EQQ_element( &
                        Ftot, index_hfcc2, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket, F1_bra, F1_ket)
        end if
    end function primitive_F1F_hyperfine_matrix_element_h2

    subroutine F1_hyperfine_structure(iverbose)
        ! Top level 
        implicit none
        integer(ik), INTENT(IN) :: iverbose
        INTEGER(ik) :: index_J, index_J_min, index_J_max 
        INTEGER(ik) :: index_F1, index_represCs, iroot, iF_ID
        INTEGER(ik) :: Ndimen_F1, Nlevels_F1, index_levels_F1, maxlocation, maxlocation2
        REAL(rk), ALLOCATABLE :: eigen_value_F1(:), &
                                    parity_conserved_F1_matrix(:, :), &
                                    primitive_F1_hyperfine_matrix(:,:), &
                                    transformation_matrix(:,:), &
                                    eigen_vector_F1(:, :)
        TYPE(quantaT), POINTER :: quanta

        I1 = max(hyperfine_setup%I1,hyperfine_setup%I2)
        if (hyperfine_setup%I1>0.0_rk) then
            active_hfcc => hfcc1
        else
            active_hfcc => hfcc2
        endif
 
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

        ! The basis of J always starts from the minimum of J
        if (mod(nint(J_global_max*2.0_rk), 2) == 0) then
            J_global_min = 0.0_rk
        else
            J_global_min = 0.5_rk
        end if

        num_F1 = nint(F1_global_max - F1_global_min) + 1
        
        CALL construct_hyperfine_constant_field_matrix(iverbose)
        

        ALLOCATE(F1_list(num_F1))
        ALLOCATE(primitive_hfs_basis(num_F1))

        CALL construct_primitive_F1_basis(iverbose)

        if (iverbose>=4) then
            write(out, '(/A2, A)') '', 'Start: hyperfine states'
            write(out, '(/A2, A, F22.12)') '', 'Zero point energy =', job%ZPE
        endif

        ALLOCATE(eigen_all_hfs(num_F1, sym%NrepresCs))
        eigen_all_hfs(:,:)%Nlevels = 0 ! This is used in F1_intensity.f90 to skip the empty block

        ! write(unit_hyperfine_states, '(A7, A22, A7, A7, A7, A7, A7, A12, A7, A7, A7, A7)') &
        !     'Number', 'Energy [cm-1]', 'g', 'F', 'I', 'parity', 'J', 'state', 'v', 'Lambda', 'Sigma', 'Omega'

        iroot = 0
        min_eigen_value_hfs = safe_max
        do index_F1 = 1, num_F1
            if (iverbose>=4) then
                write(out, '(/A4, A, F6.1)') '', 'Eigen states of F =', F1_list(index_F1)
                write(out, '(/A6, A)') '', 'Construct primitive hyperfine matrix'
            endif

            Ndimen_F1 = primitive_hfs_basis(index_F1)%Ndimen
            ALLOCATE(primitive_F1_hyperfine_matrix(Ndimen_F1, Ndimen_F1))
            CALL construct_primitive_F1_hyperfine_matrix( &
                index_F1, primitive_F1_hyperfine_matrix)

            if (iverbose>=4) write(out, '(A6, A/)') '', '... done'

            do index_represCs = 1, sym%NrepresCs
                if (iverbose>=4) then
                    write(out, '(/A6, A, A3)') '', 'Eigen states of parity', trim(parity_sign(index_represCs))
                endif

                ! Previously, the primitive hyperfine matrix was split into two equal-sized parity-conserved blocks. 
                ! Now, the number of levels is counted explicitly for each J-parity block: 
                Nlevels_F1 = 0
                J_min = abs(F1_list(index_F1) - I1)
                J_max = F1_list(index_F1) + I1
                index_J_min = nint(J_min - J_global_min) + 1
                index_J_max = nint(J_max - J_global_min) + 1

                do index_J = index_J_min, index_J_max
                    if (eigen(index_J, index_represCs)%Ndimen > 0) then
                        Nlevels_F1 = Nlevels_F1 + eigen(index_J,index_represCs)%Nlevels
                    else
                        write(out,'(A8, A, F5.1, A, A, A)') '', 'The current block with (J = ', &
                            J_global_min + REAL(index_J - 1, rk), &
                            ', parity = ', trim(parity_sign(index_represCs)), ') is empty.'
                    end if
                end do

                if (Nlevels_F1 == 0) then
                    write(out,'(/A8, A, F5.2, A, A, A)') '', 'The current block with (F1 = ', &
                        F1_list(index_F1), &
                        ', parity = ', trim(parity_sign(index_represCs)), ') is empty.'
                    cycle     ! Skip empty F-parity block
                end if

                ALLOCATE(eigen_value_F1(Nlevels_F1))
                ALLOCATE(parity_conserved_F1_matrix(Nlevels_F1, Nlevels_F1))
                ALLOCATE(transformation_matrix(Ndimen_F1, Nlevels_F1))
                ALLOCATE(eigen_vector_F1(Ndimen_F1, Nlevels_F1))

                ALLOCATE(eigen_all_hfs(index_F1, index_represCs)%quanta(Nlevels_F1))
                ALLOCATE(eigen_all_hfs(index_F1, index_represCs)%val(Nlevels_F1))
                ALLOCATE(eigen_all_hfs(index_F1, index_represCs)%vect(Ndimen_F1, Nlevels_F1))

                if (iverbose>=4) then
                    write(out, '(/A8, A)') &
                        '', 'Construct parity conserved hyperfine matrix'
                endif

                CALL construct_parity_conserved_F1_matrix( &
                    index_F1, index_represCs, Ndimen_F1, Nlevels_F1, primitive_F1_hyperfine_matrix, & 
                    parity_conserved_F1_matrix, transformation_matrix)
                ! Now, parity_conserved_F1_matrix(:, :) is the sum of
                ! parity-conserved hyperfine matrix under the rovibronic basis set
                ! and the eigenvalues of rovibronic matrix.
                ! transformation_matrix(:, :) holds the rovibronic eigenvectors of F1
                ! viz Phi                
                if (iverbose>=4) write(out, '(A8, A/)') '', '... done'

                if (iverbose>=4) then
                    write(out, '(/A8, A)') &
                        '', 'Diagonalize parity conserved hyperfine matrix'
                endif
                CALL lapack_syev(parity_conserved_F1_matrix, eigen_value_F1)
                ! Now, parity_conserved_F1_matrix(:, :) is the eigenvector matrix, viz U,
                ! under the rovibronic wavefunction.
                if (iverbose>=4) write(out, '(A8, A/)') '', '... done'
                
                if (iverbose>=4) then
                    write(out, '(/A8, A)') &
                         '', "Transform wavefunction to Hund's case (a_beta) basis set"
                endif
                ! Psi = Phi * U
                ! eigen_vector_F1 = matmul(transformation_matrix, &
                !                          parity_conserved_F1_matrix)
                CALL dgemm('N', 'N', Ndimen_F1, Nlevels_F1, Nlevels_F1, 1.0_rk, &
                            transformation_matrix, Ndimen_F1, &
                            parity_conserved_F1_matrix, Nlevels_F1, &
                            0.0_rk, eigen_vector_F1, Ndimen_F1)
                ! Now, eigen_vector_F1(:, :) is the eigenvector matrix
                ! under the Hund's case(a_beta) basis set.
                if (iverbose>=4) write(out, '(A8, A/)') '', '... done'

                eigen_all_hfs(index_F1, index_represCs)%Nlevels = Nlevels_F1
                eigen_all_hfs(index_F1, index_represCs)%Ndimen = Ndimen_F1
                eigen_all_hfs(index_F1, index_represCs)%val = eigen_value_F1
                eigen_all_hfs(index_F1, index_represCs)%vect = eigen_vector_F1
                
                if (iverbose>=4) then
                    write(out, '(/A8, A)') '', 'Calculate energy levels'
                endif
                eigen_value_F1(:) = eigen_value_F1(:) - job%ZPE

                iF_ID = 0
                do index_levels_F1 = 1, Nlevels_F1
                    maxlocation = MAXLOC(ABS(eigen_vector_F1(:, index_levels_F1)), DIM=1)
                    maxlocation2 = MAXLOC(ABS(parity_conserved_F1_matrix(:, index_levels_F1)), DIM=1)
                    quanta => primitive_hfs_basis(index_F1)%icontr(maxlocation)

                    iroot = iroot + 1
                    iF_ID = iF_ID + 1
                    eigen_all_hfs(index_F1, index_represCs)%quanta(index_levels_F1)%iroot = iroot
                    eigen_all_hfs(index_F1, index_represCs)%quanta(index_levels_F1)%iF_ID = iF_ID

                    min_eigen_value_hfs = min(&
                        eigen_all_hfs(index_F1, index_represCs)%val(index_levels_F1), &
                        min_eigen_value_hfs)

                    if (iverbose>=4) then
                        write(unit_hyperfine_states, '(I7, F22.12, I7, A7, A3, A5, A7, A12, I7, I4, 2A5)') & 
                            iroot, &
                            eigen_value_F1(index_levels_F1), &
                            int(quanta%F1 * 2.0_rk + 1.0_rk), &
                            qnum2a(quanta%F1, 7), &
                            trim(parity_sign(index_represCs)), &
                            qnum2a(quanta%I1, 5), &       
                            qnum2a(quanta%Jrot, 7), &
                            trim(poten(quanta%istate)%name), &
                            eigen_all_hfs(index_F1, index_represCs)%quanta(maxlocation2)%v, &
                            quanta%ilambda, &
                            qnum2a(quanta%sigma, 5), &
                            qnum2a(quanta%omega, 5)
                    endif
                enddo 
                if (iverbose>=4) write(out, '(A8, A/)') '', '... done'

                DEALLOCATE(eigen_vector_F1)
                DEALLOCATE(transformation_matrix)
                DEALLOCATE(parity_conserved_F1_matrix)
                DEALLOCATE(eigen_value_F1)
            end do

            DEALLOCATE(primitive_F1_hyperfine_matrix)
        end do

        if (iverbose>=4) write(out, '(A2, A/)') '', 'End: hyperfine states'    
        write(out, '(A/)') "End: hyperfine energy calculation"
        
        close(unit=unit_hyperfine_states)
    contains
        function parity_sign(index_represCs)
            INTEGER(ik), INTENT(IN) :: index_represCs
            CHARACTER(4) :: parity_sign
            if ( index_represCs == 1 ) then
                parity_sign = '+'
            elseif ( index_represCs == 2 ) then
                parity_sign = '-'
            endif
        end function parity_sign
    end subroutine F1_hyperfine_structure

    subroutine construct_hyperfine_constant_field_matrix(iverbose)
        ! <state, v| hfcc(R) |state', v'>
        implicit none
        integer(ik), INTENT(IN) :: iverbose

        TYPE(fieldT), POINTER :: field
        INTEGER(ik) :: ilevel, jlevel, index_object, index_field

        write(out, '(/A2, A)') '', "Start: hyperfine constants under the vibrational bases"

        do index_object = 1, GLOBAL_NUM_HFCC_OBJECT
            do index_field = 1, hfcc1(index_object)%num_field
                field => hfcc1(index_object)%field(index_field)
                allocate(field%matelem(vibrational_totalroots,vibrational_totalroots),stat=alloc)

                if (iverbose>=4) write(out, '(A4, A)') '', trim(field%name)

                do ilevel = 1, vibrational_totalroots
                    do jlevel = 1, ilevel
            
                    field%matelem(ilevel,jlevel)  = &
                        sum(vibrational_contrfunc(:,ilevel) &
                            *(field%gridvalue(:)) &
                            *vibrational_contrfunc(:,jlevel))
        
                    field%matelem(jlevel,ilevel) = field%matelem(ilevel,jlevel)
        
                    end do
                end do

                if (iverbose>=4) write(out, '(A4, A/)') '', '... done'
            end do

            do index_field = 1, hfcc2(index_object)%num_field
                field => hfcc2(index_object)%field(index_field)
                allocate(field%matelem(vibrational_totalroots,vibrational_totalroots),stat=alloc)
                if (iverbose>=4) write(out, '(A4, A)') '', trim(field%name)
                do ilevel = 1, vibrational_totalroots
                    do jlevel = 1, ilevel

                    field%matelem(ilevel,jlevel)  = &
                        sum(vibrational_contrfunc(:,ilevel) &
                            *(field%gridvalue(:)) &
                            *vibrational_contrfunc(:,jlevel))

                    field%matelem(jlevel,ilevel) = field%matelem(ilevel,jlevel)
                    end do
                end do
                if (iverbose>=4) write(out, '(A4, A/)') '', '... done'
            end do
        end do

        do index_object = 1, GLOBAL_NUM_HFCC12_OBJECT
            do index_field = 1, hfcc12(index_object)%num_field
                field => hfcc12(index_object)%field(index_field)
                allocate(field%matelem(vibrational_totalroots,vibrational_totalroots),stat=alloc)
                if (iverbose>=4) write(out, '(A4, A)') '', trim(field%name)
                do ilevel = 1, vibrational_totalroots
                    do jlevel = 1, ilevel

                    field%matelem(ilevel,jlevel)  = &
                        sum(vibrational_contrfunc(:,ilevel) &
                            *(field%gridvalue(:)) &
                            *vibrational_contrfunc(:,jlevel))

                    field%matelem(jlevel,ilevel) = field%matelem(ilevel,jlevel)
                    end do
                end do
                if (iverbose>=4) write(out, '(A4, A/)') '', '... done'
            end do
        end do

        write(out, '(A2, A/)') '', "End: hyperfine constants under the vibrational bases"
    end subroutine construct_hyperfine_constant_field_matrix

    subroutine construct_primitive_F1_basis(iverbose)
        ! |state, v; Lambda; S, Sigma; J, Omega; I1; F1>
        implicit none

        integer(ik), INTENT(IN) :: iverbose
        REAL(rk) :: F1, J
        INTEGER(ik) :: index_F1, Ndimen_F1
        INTEGER(ik) :: index_J, index_J_min, index_J_max, start_contr_F1, end_contr_F1, index_contr_F1
        TYPE(quantaT), POINTER :: contr_F1
        
        write(out, '(/A2, A)') '', "Start: contracted Hund's case (a_beta) basis"

        num_primitive_hfs_basis = 0
        do index_F1 = 1, num_F1
            F1 = F1_global_min + REAL(index_F1 - 1, rk)

            if (iverbose>=4) write(out, '(A4, A, F5.1)') '', "F = ", F1

            F1_list(index_F1) = F1
            J_min = abs(F1 - I1)
            J_max = F1 + I1

            index_J_min = nint(J_min - J_global_min) + 1
            index_J_max = nint(J_max - J_global_min) + 1

            Ndimen_F1 = 0
            do index_J = index_J_min, index_J_max
                Ndimen_F1 = Ndimen_F1 + basis(index_J)%Ndimen
            end do

            if (iverbose>=4) then
                write(out, '(A4, A, I5)') &
                    '', 'Number of the contracted basis', Ndimen_F1
            endif 

            primitive_hfs_basis(index_F1)%Ndimen = Ndimen_F1
            ALLOCATE (primitive_hfs_basis(index_F1)%icontr(Ndimen_F1))

            num_primitive_hfs_basis = &
                num_primitive_hfs_basis + Ndimen_F1

            end_contr_F1 = 0
            do index_J = index_J_min, index_J_max
                J = J_global_min + REAL(index_J - 1, rk)
                start_contr_F1 = end_contr_F1
                end_contr_F1 = start_contr_F1 + basis(index_J)%Ndimen
                primitive_hfs_basis(index_F1)%icontr(start_contr_F1 + 1:end_contr_F1) = &
                    basis(index_J)%icontr(:)
                primitive_hfs_basis(index_F1)%icontr(start_contr_F1 + 1:end_contr_F1)%Jrot = J
            end do

            primitive_hfs_basis(index_F1)%icontr(:)%index_F1 = index_F1
            primitive_hfs_basis(index_F1)%icontr(:)%I1 = I1
            primitive_hfs_basis(index_F1)%icontr(:)%F1 = F1
            if (iverbose>=4) then
                write(out, '(A6, A7, A7, A7, A7, A7, A7, A7, A7, A7, A7)') &
                    '    ', 'No.', 'F', 'I', 'state', 'v', 'Lambda', 'J', 'Omega', 'S', 'Sigma'
            endif
            do index_contr_F1 = 1, Ndimen_F1
                contr_F1 => primitive_hfs_basis(index_F1)%icontr(index_contr_F1)
                if (iverbose>=4) then
                    write(out, '(A6, I7, F7.1, F7.1, I7, I7, I7, F7.1, F7.1, F7.1, F7.1)') &
                        '    ', index_contr_F1, contr_F1%F1, contr_F1%I1, &
                        contr_F1%istate, contr_F1%v, contr_F1%ilambda, &
                        contr_F1%Jrot, contr_F1%omega, &
                        contr_F1%spin, contr_F1%sigma 
                endif
            end do

            if (iverbose>=4) write(out, '(A4, A/)') '', '... done'
        end do
       
        write(out, '(A2, A/)') '', "End: contracted Hund's case (a_beta) basis"
    end subroutine construct_primitive_F1_basis

    subroutine construct_parity_conserved_F1_matrix( &
        index_F1, index_represCs, Ndimen_F1, Nlevels_F1,&
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

        ! H_hfs^tau = M_intermediate * Phi           
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

        DEALLOCATE(eigen_value_J_in_F1)
    end subroutine construct_parity_conserved_F1_matrix

    subroutine construct_F1_transformation_matrix( &
        index_F1, index_represCs, &
        transformation_matrix, &
        eigen_value_J_in_F1)
        ! Construct the basis tranformation matrix, Phi, 
        ! between the primitive F1 basis set
        ! and the parity conserved F1 basis set.

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
                + basis(index_J)%Ndimen
            if (eigen(index_J, index_represCs)%Nlevels == 0) cycle   
            ! if (eigen(index_J, index_represCs)%Ndimen == 0) then
            !     write(out,*) 'ERROR: Nlevels > 0 but Ndimen = 0 in construct_F1_transformation_matrix'
            !     stop
            ! end if             
            end_level_F1 = start_level_F1 &
                + eigen(index_J, index_represCs)%Nlevels

            ! Eigenvalues of H0
            eigen_value_J_in_F1(start_level_F1 + 1 : end_level_F1) &
                = eigen(index_J, index_represCs)%val

            ! Phi
            transformation_matrix(start_dimen_F1 + 1 : end_dimen_F1, &
                                  start_level_F1 + 1 : end_level_F1) &
                = eigen(index_J, index_represCs)%vect

            ! eigen_all_hfs is the table of eigen states
            ! corresponding to all F1 values.
            eigen_all_hfs(index_F1, index_represCs)% &
                    quanta(start_level_F1 + 1 : end_level_F1) &
                = eigen(index_J, index_represCs)%quanta
            eigen_all_hfs(index_F1, index_represCs)% &
                    quanta(start_level_F1 + 1 : end_level_F1)%F1 &
                = F1
            eigen_all_hfs(index_F1, index_represCs)% &
                    quanta(start_level_F1 + 1 : end_level_F1)%I1 &
                = I1
            eigen_all_hfs(index_F1, index_represCs)% &
                    quanta(start_level_F1 + 1 : end_level_F1)%index_F1 &
                = index_F1
            eigen_all_hfs(index_F1, index_represCs)% &
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
        INTEGER(ik) :: Ndimen_F1, index_hfcc1, index_field

        primitive_F1_hyperfine_matrix = 0

        Ndimen_F1 = primitive_hfs_basis(index_F1)%Ndimen

        do ket = 1, Ndimen_F1
            call get_quanta(index_F1, ket, &
                            F1, state_ket, index_v_ket, v_ket, Lambda_ket, &
                            S_ket, Sigma_ket, J_ket, Omega_ket)

            do bra = 1, Ndimen_F1
                call get_quanta(index_F1, bra, &
                                F1, state_bra, index_v_bra, v_bra, Lambda_bra, &
                                S_bra, Sigma_bra, J_bra, Omega_bra)

                if (state_bra == state_ket) then

                    do index_hfcc1 = 1, GLOBAL_NUM_HFCC_OBJECT
                        do index_field = 1, active_hfcc(index_hfcc1)%num_field
                            if ( (state_bra == active_hfcc(index_hfcc1)%field(index_field)%istate) &
                                .and. (state_ket == active_hfcc(index_hfcc1)%field(index_field)%jstate)) then

                                primitive_F1_hyperfine_matrix(bra, ket) = &
                                    primitive_F1_hyperfine_matrix(bra, ket) &
                                    + primitive_F1_hyperfine_matrix_element( &
                                        F1, index_hfcc1, index_field, & 
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

    subroutine get_quanta(index_F_block, index_contr, &
        Ftot, state, index_v, v, Lambda, S, Sigma, J, Omega, F1, Itot)
        ! Get quanta from the active hyperfine basis.
        implicit none

        INTEGER(ik), INTENT(IN) :: index_F_block, index_contr
        INTEGER(ik), INTENT(OUT) :: state, Lambda, index_v, v
        REAL(rk), INTENT(OUT) :: Ftot, S, Sigma, J, Omega
        REAL(rk), INTENT(OUT), OPTIONAL :: F1   ! optional output for heteronuclear cases
        REAL(rk), INTENT(OUT), OPTIONAL :: Itot ! optional output for equal coupling case
        TYPE(quantaT), POINTER :: icontr
        
        icontr => primitive_hfs_basis(index_F_block)%icontr(index_contr)

        if (allocated(Ftot_list)) then
            Ftot = icontr%Ftot
        else
            Ftot = icontr%F1
        end if

        state = icontr%istate
        Lambda = icontr%ilambda
        index_v = icontr%ivib
        v = icontr%v
        
        S = icontr%spin
        Sigma = icontr%sigma
        J = icontr%Jrot
        Omega = icontr%omega 
        if (present(F1)) F1 = icontr%F1
        if (present(Itot)) Itot = icontr%Itot

    end subroutine get_quanta
    
    function primitive_F1_hyperfine_matrix_element( &
        F1, index_hfcc1, index_field, & 
        index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_hfs)

        implicit none

        REAL(rk), INTENT(IN) :: F1, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        INTEGER(ik), INTENT(IN) :: index_hfcc1, index_field, &
                                   index_v_bra, index_v_ket, &
                                   Lambda_bra, Lambda_ket
        REAL(rk) :: H_hfs

        ! The switch flow could be refactored if 
        ! an abstract interface is defined.
        ! This function could be removed then.
        H_hfs = 0.0_rk
        if (index_hfcc1 == 1) then
            H_hfs = H_FC_element( &
                        F1, index_field, index_v_bra, index_v_ket,&
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket) 
        else if (index_hfcc1 == 2) then
            H_hfs =  H_IL_element( &
                        F1, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket) 
        else if (index_hfcc1 == 3) then
            H_hfs = H_DIP_C_element( &
                        F1, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket) 
        else if (index_hfcc1 == 4) then 
            H_hfs = H_DIP_D_element( &
                        F1, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket) 
        else if (index_hfcc1 == 5) then
            H_hfs = H_IJ_element( &
                        F1, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket) 
        else if (index_hfcc1 == 6) then
            H_hfs = H_CIp_element( &
                        F1, index_field, index_v_bra, index_v_ket, &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket)
        else if (index_hfcc1 == 7) then
            H_hfs = H_EQQ0_element( &        
                        F1, index_field, index_v_bra, index_v_ket,  &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket)
        else if (index_hfcc1 == 8) then
            H_hfs = H_EQQ2_element( &        
                        F1, index_field, index_v_bra, index_v_ket,  &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket)
        end if

    end function primitive_F1_hyperfine_matrix_element

    function qnum2a(x, w) result(s)
        ! Convert a quantum number to a fixed-width string of width w.
        ! Integer -> I0; otherwise -> F0.1 (half-integers).
        ! Right-aligned within width w. If it does not fit, fill with '*'.
        ! Also enforce leading zero: .5 -> 0.5, -.5 -> -0.5.
        real(rk), intent(in) :: x
        integer,  intent(in) :: w
        character(len=:), allocatable :: s
        real(rk) :: tol
        integer  :: xi
        character(len=64) :: tmp
        character(len=:), allocatable :: t

        tol = 1000.0_rk * epsilon(1.0_rk)

        ! Build minimal text in a clean buffer
        tmp = ' '
        if (abs(x - real(nint(x), rk)) < tol) then
            xi = nint(x)
            write(tmp,'(I0)') xi        ! minimal-width integer
        else
            write(tmp,'(F0.1)') x       ! minimal-width real with 1 decimal
        end if
        t = trim(adjustl(tmp))         ! minimal text (no leading/trailing blanks)

        ! Add leading zero for fractional numbers: .5 -> 0.5, -.5 -> -0.5, +.5 -> +0.5
        if (len(t) >= 2) then
            if (t(1:2) == '-.') t = '-0' // t(2:)
            if (t(1:2) == '+.') t = '+0' // t(2:)
        end if
        if (len(t) >= 1) then
            if (t(1:1) == '.')  t = '0' // t
        end if

        ! Allocate fixed-width result and right-align
        allocate(character(len=w) :: s)
        if (len(t) > w) then
            s = repeat('*', w)         ! overflow marker
        else
            s = repeat(' ', w-len(t)) // t
        end if
    end function qnum2a

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
                    * active_hfcc(1)%field(index_field)%matelem(index_v_bra, index_v_ket)
                    
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
                    * active_hfcc(2)%field(index_field)%matelem(index_v_bra, index_v_ket)
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
        INTEGER(ik) :: p, q

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
                    get_sign(J_bra - Omega_bra + S_bra - Sigma_bra + p + q) & 
                    * sqrt((2.0_rk * S_bra + 1.0_rk) * (2.0_rk * J_ket + 1.0_rk)) &
                    * Wigner3j(1.0_rk, 1.0_rk, 2.0_rk, &
                                REAL(-q, rk), REAL(p+q, rk), REAL(-p, rk)) &
                    * Wigner3j(S_bra, 1.0_rk, S_ket, &
                                -Sigma_bra, REAL(q, rk), Sigma_ket) &
                    * Wigner3j(J_ket, J_bra, 1.0_rk, &
                                Omega_ket, -Omega_bra, REAL(p+q, rk)) &
                    * sqrt(30.0_rk) * sqrt(S_bra * (S_bra + 1.0_rk)) &
                    * active_hfcc(3)%field(index_field)%matelem(index_v_bra, index_v_ket) / (3.0_rk)
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
        INTEGER(ik) :: p, q, q1

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
                        get_sign(J_bra - Omega_bra + S_bra - Sigma_bra + p + q) & 
                        * sqrt((2.0_rk * S_bra + 1.0_rk) * (2.0_rk * J_ket + 1.0_rk)) &
                        * Wigner3j(1.0_rk, 1.0_rk, 2.0_rk, &
                                    REAL(-q, rk), REAL(p+q, rk), REAL(-p, rk)) &
                        * Wigner3j(S_bra, 1.0_rk, S_ket, &
                                    -Sigma_bra, REAL(q, rk), Sigma_ket) &
                        * Wigner3j(J_ket, J_bra, 1.0_rk, &
                                    Omega_ket, -Omega_bra, REAL(p+q, rk)) &
                        * sqrt(30.0_rk) * sqrt(S_bra * (S_bra + 1.0_rk)) &
                        * active_hfcc(4)%field(index_field)%matelem(index_v_bra, index_v_ket) /(-sqrt(6.0_rk))
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
                * active_hfcc(5)%field(index_field)%matelem(index_v_bra, index_v_ket)
        end if

    end function H_IJ_element

    function H_CIp_element( &
        F1, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_CIp)

        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: F1, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        REAL(rk) :: H_CIp, reduced_A1

        ! Off-diagonal nuclear spin-rotation for nucleus 1
        reduced_A1 = rank1_reduced_matrix_element_A1( &
            Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
            J_bra, J_ket, Omega_bra, Omega_ket)

        H_CIp = get_sign(J_ket + I1 + F1) &
            * Wigner6j(I1, J_bra, F1, J_ket, I1, 1.0_rk) &
            * sqrt(I1*(I1 + 1.0_rk)*(2.0_rk*I1 + 1.0_rk)) &
            * reduced_A1 &
            * active_hfcc(6)%field(index_field)%matelem(index_v_bra, index_v_ket)

    end function H_CIp_element
    
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
        
        REAL(rk) :: H_EQ
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
                * active_hfcc(7)%field(index_field)%matelem(index_v_bra, index_v_ket) / (-2.0_rk)
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
        
        REAL(rk) :: H_EQ
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
                    * active_hfcc(8)%field(index_field)%matelem(index_v_bra, index_v_ket) / (-sqrt(24.0_rk))
                exit
            end if
        end do
    end function H_EQQ2_element

    function Heq_FC_element( &
        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket,  &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_FC)

        implicit none

        integer(ik), intent(in) :: index_field, index_v_bra, index_v_ket
        integer(ik), intent(in) :: Lambda_bra, Lambda_ket
        real(rk), intent(in) :: Ftot, Itot_bra, Itot_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket

        real(rk) :: H_FC
        integer(ik) :: q

        ! Uses b_F(R) curve stored in hfcc1(1)%field(index_field)%matelem(v_bra,v_ket)
        
        H_FC = 0.0_rk
        do q = -1, 1, 1
            if ( (Lambda_bra == Lambda_ket) &
                .and. (nint(S_bra - S_ket) == 0) &
                .and. (nint(Sigma_ket - (Sigma_bra - real(q, rk))) == 0) &
                .and. (nint(Itot_bra - Itot_ket) == 0) ) then

                H_FC = H_FC + &
                    get_sign(J_ket + Itot_bra + Ftot) &   ! (-1)^(J' + I + F)
                    * Wigner6j(Itot_bra, J_bra, Ftot, J_ket, Itot_ket, 1.0_rk) &  ! { I  J  F ; J' I' 1 }
                    * sqrt(Itot_bra*(Itot_bra + 1.0_rk)*(2.0_rk*Itot_bra + 1.0_rk)) &
                    * get_sign(J_bra - Omega_bra) &
                    * Wigner3j(J_bra, 1.0_rk, J_ket, -Omega_bra, real(q, rk), Omega_ket) &
                    * sqrt((2.0_rk*J_bra + 1.0_rk)*(2.0_rk*J_ket + 1.0_rk)) &
                    * get_sign(S_bra - Sigma_bra) &
                    * Wigner3j(S_bra, 1.0_rk, S_ket, -Sigma_bra, real(q, rk), Sigma_ket) &
                    * sqrt(S_bra*(S_bra + 1.0_rk)*(2.0_rk*S_bra + 1.0_rk)) &
                    * hfcc1(1)%field(index_field)%matelem(index_v_bra, index_v_ket)
            end if
        end do

    end function Heq_FC_element

    function Heq_IL_element( &
        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket,  &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_IL)

        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, &
                                   Lambda_bra, Lambda_ket
        REAL(rk), INTENT(in) :: Ftot, Itot_bra, Itot_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket 

        REAL(rk) :: H_IL

        H_IL = 0.0_rk 
        
        ! Selection rules from the diagonal-Λ form (q=0):
        
        if ( (Lambda_bra == Lambda_ket) &
            .and. (nint(S_bra - S_ket) == 0 ) &
            .and. (nint(Sigma_bra - Sigma_ket) == 0 ) &
            .and. (nint(Itot_bra - Itot_ket) == 0)  &
            .and. ( nint(Omega_bra - Omega_ket) == 0 )) then

        ! (-1)^{J'+F+I} { I  J  F ; J' I 1 } sqrt(I(I+1)(2I+1))
        ! × (-1)^{J-Ω} ( J 1 J' ; -Ω 0 Ω' ) sqrt((2J+1)(2J'+1))
        ! × Λ <v|a(R)|v'>
            H_IL = get_sign(J_ket + Itot_bra + Ftot) &
                * Wigner6j(Itot_bra, J_bra, Ftot, J_ket, Itot_ket, 1.0_rk) &
                * sqrt(Itot_bra * (Itot_bra + 1.0_rk) * (2.0_rk*Itot_bra + 1.0_rk)) &
                * get_sign(J_bra - Omega_bra) &
                * Wigner3j(J_bra, 1.0_rk, J_ket, -Omega_bra, 0.0_rk, Omega_ket) &
                * sqrt((2.0_rk*J_bra + 1.0_rk) * (2.0_rk*J_ket + 1.0_rk)) &
                * real(Lambda_bra, rk) &
                * hfcc1(2)%field(index_field)%matelem(index_v_bra, index_v_ket)
        end if 

    end function Heq_IL_element

    function Heq_DIP_C_element( &
        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket,  &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_IS)
        
        implicit none

        integer(ik), intent(in) :: index_field, index_v_bra, index_v_ket
        integer(ik), intent(in) :: Lambda_bra, Lambda_ket
        real(rk), intent(in) :: Ftot, Itot_bra, Itot_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        real(rk) :: H_IS
        integer(ik) :: q

        H_IS = 0.0_rk

        do q = -1, 1, 1
            if ( (Lambda_bra == Lambda_ket) &
                .and. (nint(S_bra - S_ket) == 0) &
                .and. (nint(Itot_bra - Itot_ket) == 0) &
                .and. (nint(Omega_ket - (Omega_bra - real(q, rk))) == 0) &
                .and. (nint(Sigma_ket - (Sigma_bra - real(q, rk))) == 0) ) then

                H_IS = H_IS + &
                    sqrt(30.0_rk) &
                    * get_sign(J_ket + Ftot + Itot_bra) &  ! (-1)^(J' + F + I)
                    * Wigner6j(Itot_bra, J_bra, Ftot, J_ket, Itot_ket, 1.0_rk) &  ! { I  J  F ; J' I' 1 }
                    * sqrt(Itot_bra*(Itot_bra + 1.0_rk)*(2.0_rk*Itot_bra + 1.0_rk)) &
                    * get_sign(J_bra - Omega_bra) &
                    * Wigner3j(J_bra, 1.0_rk, J_ket, -Omega_bra, real(q, rk), Omega_ket) &
                    * sqrt((2.0_rk*J_bra + 1.0_rk)*(2.0_rk*J_ket + 1.0_rk)) &
                    * get_sign(real(q, rk)) &  ! (-1)^q
                    * Wigner3j(1.0_rk, 2.0_rk, 1.0_rk, real(q, rk), 0.0_rk, -real(q, rk)) &
                    * get_sign(S_bra - Sigma_bra) &
                    * Wigner3j(S_bra, 1.0_rk, S_ket, -Sigma_bra, real(q, rk), Sigma_ket) &
                    * sqrt(S_bra*(S_bra + 1.0_rk)*(2.0_rk*S_bra + 1.0_rk)) &
                    * hfcc1(3)%field(index_field)%matelem(index_v_bra, index_v_ket)/(3.0_rk)

            end if
        end do
    end function Heq_DIP_C_element

    function Heq_DIP_D_element( &
        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket,  &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_IS)
        
        implicit none

        integer(ik), intent(in) :: index_field, index_v_bra, index_v_ket
        integer(ik), intent(in) :: Lambda_bra, Lambda_ket
        real(rk), intent(in) :: Ftot, Itot_bra, Itot_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        real(rk) :: H_IS
        integer(ik) :: q1, q2, q

        H_IS = 0.0_rk


        do q2 = -2, 2, 4
            do q1 = -1, 1, 2
                q = q1 + q2
                if (abs(q) > 1) cycle
                if ( (Lambda_ket == (Lambda_bra - q2)) &
                    .and. (nint(S_bra - S_ket) == 0) &
                    .and. (nint(Sigma_bra - Sigma_ket) == q1) &
                    .and. (nint(Omega_bra - Omega_ket) == q) &
                    .and. (nint(Itot_bra - Itot_ket) == 0) ) then
                    H_IS = H_IS + &
                        sqrt(30.0_rk) &
                        * get_sign(J_ket + Ftot + Itot_bra) &  ! (-1)^(J' + F + I)
                        * Wigner6j(Itot_bra, J_bra, Ftot, J_ket, Itot_ket, 1.0_rk) &  ! { I  J  F ; J' I' 1 }
                        * sqrt(Itot_bra*(Itot_bra + 1.0_rk)*(2.0_rk*Itot_bra + 1.0_rk)) &
                        * get_sign(J_bra - Omega_bra) &
                        * Wigner3j(J_bra, 1.0_rk, J_ket, -Omega_bra, real(q, rk), Omega_ket) &
                        * sqrt((2.0_rk*J_bra + 1.0_rk)*(2.0_rk*J_ket + 1.0_rk)) &
                        * get_sign(real(q, rk)) &  ! (-1)^q
                        * Wigner3j(1.0_rk, 2.0_rk, 1.0_rk, real(q1, rk), real(q2, rk), -real(q, rk)) &
                        * get_sign(S_bra - Sigma_bra) &
                        * Wigner3j(S_bra, 1.0_rk, S_ket, -Sigma_bra, real(q1, rk), Sigma_ket) &
                        * sqrt(S_bra*(S_bra + 1.0_rk)*(2.0_rk*S_bra + 1.0_rk)) &
                        * hfcc1(4)%field(index_field)%matelem(index_v_bra, index_v_ket)/(-sqrt(6.0_rk))
                end if
            end do            
        end do
         
    end function Heq_DIP_D_element
    
    function Heq_IJ_element( &
        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket,  &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_IJ)

        implicit none

        integer(ik), intent(in) :: index_field, index_v_bra, index_v_ket
        integer(ik), intent(in) :: Lambda_bra, Lambda_ket
        real(rk), intent(in) :: Ftot, Itot_bra, Itot_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        
        REAL(rk) :: H_IJ

        ! Nuclear spin-rotation interaction
        ! See Brown and Carrington (2003) Eq.(8.20)
        ! https://doi.org/10.1017/cbo9780511814808

        H_IJ = 0.0_rk

        ! selection rules: Lambda,S,Sigma,J,Omega conserved; and I=I'
        if ( (Lambda_bra == Lambda_ket) &
            .and. (nint(S_bra - S_ket) == 0) &
            .and. (nint(Sigma_bra - Sigma_ket) == 0) &
            .and. (nint(J_bra - J_ket) == 0) &
            .and. (nint(Omega_bra - Omega_ket) == 0) &
            .and. (nint(Itot_bra - Itot_ket) == 0) ) then

            H_IJ = get_sign(J_ket + Itot_bra + Ftot) &   ! (-1)^{J'+I+F}
                * Wigner6j(Itot_bra, J_bra, Ftot, J_ket, Itot_ket, 1.0_rk) &
                * sqrt( J_bra*(J_bra + 1.0_rk)*(2.0_rk*J_bra + 1.0_rk) &
                    * Itot_bra*(Itot_bra + 1.0_rk)*(2.0_rk*Itot_bra + 1.0_rk) ) &
                * hfcc1(5)%field(index_field)%matelem(index_v_bra, index_v_ket)
        end if
    end function Heq_IJ_element

    function Heq_CIp_element( &
        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_CIp)

        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: Ftot, Itot_bra, Itot_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        REAL(rk) :: H_CIp, reduced_A1

        H_CIp = 0.0_rk
        if (nint(Itot_bra - Itot_ket) /= 0) return

        reduced_A1 = rank1_reduced_matrix_element_A1( &
            Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
            J_bra, J_ket, Omega_bra, Omega_ket)

        H_CIp = get_sign(J_ket + Itot_bra + Ftot) &
            * Wigner6j(Itot_bra, J_bra, Ftot, J_ket, Itot_ket, 1.0_rk) &
            * sqrt(Itot_bra*(Itot_bra + 1.0_rk)*(2.0_rk*Itot_bra + 1.0_rk)) &
            * reduced_A1 &
            * hfcc1(6)%field(index_field)%matelem(index_v_bra, index_v_ket)

    end function Heq_CIp_element

    function Heq_EQQ0_element( &
        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_EQ)
        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: Ftot, Itot_bra, Itot_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket 
        REAL(rk) :: H_EQ
        INTEGER(ik) :: q

        H_EQ = 0.0_rk
        ! Nuclear quadrupole exists only if single-nucleus spin I1 >= 1
        if ( I1 < 1.0_rk ) return 
        q = 0
        if ( (Lambda_ket == Lambda_bra - q) & 
            .and. (nint(S_ket - S_bra) == 0) &
            .and. (nint(Sigma_ket - Sigma_bra) == 0) ) then
            H_EQ = get_sign(J_ket + 2.0_rk*I1 + Itot_bra + Ftot) &     ! (-1)^{J' + 2I1 + I + F}
                * get_sign(J_bra - Omega_bra) &                     ! (-1)^{J - Omega}
                * (get_sign(Itot_ket) + get_sign(Itot_bra))& 
                * sqrt( (2.0_rk*J_bra + 1.0_rk) * (2.0_rk*J_ket + 1.0_rk) ) &
                * Wigner3j(J_bra, 2.0_rk, J_ket, -Omega_bra, 0.0_rk, Omega_ket) &
                * Wigner6j(Itot_bra, J_bra, Ftot, J_ket, Itot_ket, 2.0_rk) &      ! { I  J  F ; J' I' 2 }
                * Wigner6j(I1, Itot_bra, I1, Itot_ket, I1, 2.0_rk) &            ! { I1 I  I1; I' I1 2 }
                * sqrt( (2.0_rk*Itot_bra + 1.0_rk) * (2.0_rk*Itot_ket + 1.0_rk) ) &
                / Wigner3j(I1, 2.0_rk, I1, -I1, 0.0_rk, I1) &                  ! ( I1 2 I1; -I1 0 I1 )^{-1}
                * hfcc1(7)%field(index_field)%matelem(index_v_bra, index_v_ket)/4.0_rk
        end if
    end function Heq_EQQ0_element

    function Heq_EQQ2_element( &
        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_EQ)
        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: Ftot, Itot_bra, Itot_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket 
        REAL(rk) :: H_EQ
        INTEGER(ik) :: q

        H_EQ = 0.0_rk
        ! Nuclear quadrupole exists only if single-nucleus spin I1 >= 1
        if ( I1 < 1.0_rk ) return 
        do q = -2, 2, 4
            if ( (Lambda_ket == Lambda_bra - q) & 
                .and. (nint(S_ket - S_bra) == 0) &
                .and. (nint(Sigma_ket - Sigma_bra) == 0) &
                .and. (nint(Omega_bra - Omega_ket) == q)) then
                H_EQ = get_sign(J_ket + 2.0_rk*I1 + Itot_bra + Ftot) &     ! (-1)^{J' + 2I1 + I + F}
                    * get_sign(J_bra - Omega_bra) &                     ! (-1)^{J - Omega}
                    * (get_sign(Itot_ket) + get_sign(Itot_bra))& 
                    * sqrt( (2.0_rk*J_bra + 1.0_rk) * (2.0_rk*J_ket + 1.0_rk) ) &
                    * Wigner3j(J_bra, 2.0_rk, J_ket, -Omega_bra, real(q, rk), Omega_ket) &
                    * Wigner6j(Itot_bra, J_bra, Ftot, J_ket, Itot_ket, 2.0_rk) &      ! { I  J  F ; J' I' 2 }
                    * Wigner6j(I1, Itot_bra, I1, Itot_ket, I1, 2.0_rk) &            ! { I1 I  I1; I' I1 2 }
                    * sqrt( (2.0_rk*Itot_bra + 1.0_rk) * (2.0_rk*Itot_ket + 1.0_rk) ) &
                    / Wigner3j(I1, 2.0_rk, I1, -I1, 0.0_rk, I1) &                  ! ( I1 2 I1; -I1 0 I1 )^{-1}
                    * hfcc1(8)%field(index_field)%matelem(index_v_bra, index_v_ket)/(4.0_rk*sqrt(6.0_rk))
            end if
        end do
    end function Heq_EQQ2_element


    function Heq_II_element( &
        Ftot, Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_II)
        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: Ftot, Itot_bra, Itot_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket 
        
        REAL(rk) :: H_II

        H_II = 0.0_rk
        
        if ( (Lambda_bra == Lambda_ket) &
            .and. (nint(S_bra - S_ket) == 0) &
            .and. (nint(Sigma_bra - Sigma_ket) == 0) ) then
            H_II =  -get_sign(J_ket + Itot_bra + Ftot) &               ! (-1)^{J' + I + F}
                * get_sign(J_bra - Omega_bra) &                     ! (-1)^{J - Omega}
                * sqrt( (2.0_rk*J_bra + 1.0_rk) * (2.0_rk*J_ket + 1.0_rk) &
                        * (2.0_rk*Itot_bra + 1.0_rk) * (2.0_rk*Itot_ket + 1.0_rk) ) &
                * Wigner3j(J_bra, 2.0_rk, J_ket, -Omega_bra, 0.0_rk, Omega_ket) &
                * Wigner6j(Itot_bra, J_bra, Ftot, J_ket, Itot_ket, 2.0_rk) &       ! { I  J  F ; J' I' 2 }
                * Wigner9j(I1, I1, 1.0_rk, I1, I1, 1.0_rk, Itot_ket, Itot_bra, 2.0_rk) &  ! {I1 I1 1; I1 I1 1; I' I 2}
                * (I1*(I1 + 1.0_rk)*(2.0_rk*I1 + 1.0_rk)) &
                * sqrt(30.0_rk) &
                * hfcc12(1)%field(index_field)%matelem(index_v_bra, index_v_ket)
        end if 
    end function Heq_II_element

    function Heq_IIS_element( &
        Itot_bra, Itot_ket, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) result(H_IIS)
        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: Itot_bra, Itot_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        REAL(rk) :: H_IIS

        H_IIS = 0.0_rk
        ! c4(R) I1.I2 in the (I1 I2)I, J -> F basis; here I2 = I1.
        ! Angular factor from Brown & Carrington, Eq. (8.177)
        if ( (Lambda_bra == Lambda_ket) &
            .and. (nint(S_bra - S_ket) == 0) &
            .and. (nint(Sigma_bra - Sigma_ket) == 0) &
            .and. (nint(J_bra - J_ket) == 0) &
            .and. (nint(Omega_bra - Omega_ket) == 0) &
            .and. (nint(Itot_bra - Itot_ket) == 0) ) then
            H_IIS = 0.5_rk * (Itot_bra*(Itot_bra + 1.0_rk) - 2.0_rk*I1*(I1 + 1.0_rk)) &
                * hfcc12(2)%field(index_field)%matelem(index_v_bra, index_v_ket)
        end if
    end function Heq_IIS_element

    function H2_MD_element( &
        Ftot, index_hfcc2, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket, F1_bra, F1_ket) &
        result(H_MD)

        implicit none

        INTEGER(ik), INTENT(IN) :: index_hfcc2, index_field, index_v_bra, index_v_ket, &
                                   Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: Ftot, F1_bra, F1_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        REAL(rk) :: H_MD, reduced_M2, rotational_factor, spin_factor, q_real
        INTEGER(ik) :: q

        ! nucleus-2 magnetic terms in unequal coupling scheme F1 = J + I1, F = F1 + I2.
        ! reduced_M2 is <k,J||T^1(M2)||k',J'> including sqrt(2J_ket+1)(2J_bra+1) factors
        
        H_MD = 0.0_rk
        if (I2 <= 0.0_rk .or. nint(S_bra - S_ket) /= 0) return
        if (abs(J_bra - J_ket) > 1.0_rk .or. abs(F1_bra - F1_ket) > 1.0_rk) return

        reduced_M2 = 0.0_rk
        rotational_factor = get_sign(J_bra - Omega_bra) &
            * sqrt((2.0_rk*J_bra + 1.0_rk)*(2.0_rk*J_ket + 1.0_rk))
        spin_factor = get_sign(S_bra - Sigma_bra) &
            * sqrt(S_bra*(S_bra + 1.0_rk)*(2.0_rk*S_bra + 1.0_rk))

        select case (index_hfcc2)
        case (1) ! b_F: Fermi contact
            if (Lambda_bra /= Lambda_ket) return
            do q = -1, 1
                q_real = real(q, rk)
                if (nint(Sigma_bra - Sigma_ket) /= q) cycle
                reduced_M2 = reduced_M2 + rotational_factor * spin_factor &
                    * Wigner3j(J_bra, 1.0_rk, J_ket, -Omega_bra, q_real, Omega_ket) &
                    * Wigner3j(S_bra, 1.0_rk, S_ket, -Sigma_bra, q_real, Sigma_ket)
            end do
        case (2) ! a: orbital term diagonal in Lambda
            if (Lambda_bra /= Lambda_ket .or. nint(Sigma_bra - Sigma_ket) /= 0) return
            reduced_M2 = rotational_factor * real(Lambda_bra, rk) &
                * Wigner3j(J_bra, 1.0_rk, J_ket, -Omega_bra, 0.0_rk, Omega_ket)
        case (3) ! c: electron-nuclear dipolar term diagonal in Lambda
            if (Lambda_bra /= Lambda_ket) return
            do q = -1, 1
                q_real = real(q, rk)
                if (nint(Sigma_bra - Sigma_ket) /= q) cycle
                reduced_M2 = reduced_M2 + sqrt(10.0_rk/3.0_rk) &
                    * rotational_factor * spin_factor * get_sign(q_real) &
                    * Wigner3j(J_bra, 1.0_rk, J_ket, -Omega_bra, q_real, Omega_ket) &
                    * Wigner3j(1.0_rk, 2.0_rk, 1.0_rk, q_real, 0.0_rk, -q_real) &
                    * Wigner3j(S_bra, 1.0_rk, S_ket, -Sigma_bra, q_real, Sigma_ket)
            end do
        case (4) ! d: Lambda = +/-1 components of a Pi state
            if (abs(Lambda_bra) /= 1 .or. Lambda_ket /= -Lambda_bra) return
            do q = -1, 1, 2
                q_real = real(q, rk)
                if (Lambda_bra - Lambda_ket /= 2*q) cycle
                if (nint(Sigma_bra - Sigma_ket) /= -q) cycle
                reduced_M2 = reduced_M2 - sqrt(5.0_rk) &
                    * rotational_factor * spin_factor * get_sign(q_real) &
                    * Wigner3j(J_bra, 1.0_rk, J_ket, -Omega_bra, q_real, Omega_ket) &
                    * Wigner3j(1.0_rk, 2.0_rk, 1.0_rk, -q_real, 2.0_rk*q_real, -q_real) &
                    * Wigner3j(S_bra, 1.0_rk, S_ket, -Sigma_bra, -q_real, Sigma_ket)
            end do
        end select

        ! recouple the rovibronic rank-one tensor through F1 and I2.
        ! In particular, the phase contains J_bra and twice F1_ket.
        H_MD = get_sign(1.0_rk + J_bra + I1 + I2 + 2.0_rk*F1_ket + Ftot) &
            * Wigner6j(J_ket, F1_ket, I1, F1_bra, J_bra, 1.0_rk) &
            * Wigner6j(I2, F1_ket, Ftot, F1_bra, I2, 1.0_rk) &
            * sqrt((2.0_rk*F1_bra + 1.0_rk)*(2.0_rk*F1_ket + 1.0_rk)) &
            * sqrt(I2*(I2 + 1.0_rk)*(2.0_rk*I2 + 1.0_rk)) &
            * reduced_M2 * hfcc2(index_hfcc2)%field(index_field)%matelem(index_v_bra, index_v_ket)
    end function H2_MD_element

    function H2_IJ_element( &
        Ftot, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket, F1_bra, F1_ket) &
        result(H_IJ)

        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: Ftot, F1_bra, F1_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket

        REAL(rk) :: H_IJ, rank1_reduced_matrix_element_J, &
                    rank1_reduced_matrix_element_F1, rank1_reduced_matrix_element_I2

        H_IJ = 0.0_rk

        if ( (Lambda_bra == Lambda_ket) &
            .and. (nint(S_bra - S_ket) == 0) &
            .and. (nint(Sigma_bra - Sigma_ket) == 0) &
            .and. (nint(J_bra - J_ket) == 0) &
            .and. (nint(Omega_bra - Omega_ket) == 0) ) then
            rank1_reduced_matrix_element_J = &
                sqrt(J_bra * (J_bra + 1.0_rk) * (2.0_rk * J_bra + 1.0_rk)) &
                * hfcc2(5)%field(index_field)%matelem(index_v_bra, index_v_ket)

            rank1_reduced_matrix_element_F1 = &
                get_sign(J_bra + I1 + F1_ket + 1.0_rk) &
                * sqrt((2.0_rk*F1_bra + 1.0_rk) * (2.0_rk*F1_ket + 1.0_rk)) &
                * Wigner6j(J_bra, F1_bra, I1, F1_ket, J_ket, 1.0_rk) &
                * rank1_reduced_matrix_element_J

            rank1_reduced_matrix_element_I2 = sqrt(I2 * (I2 + 1.0_rk) * (2.0_rk*I2 + 1.0_rk))

            H_IJ = get_sign(F1_ket + I2 + Ftot) &
                * Wigner6j(F1_bra, I2, Ftot, I2, F1_ket, 1.0_rk) &
                * rank1_reduced_matrix_element_F1 &
                * rank1_reduced_matrix_element_I2
        end if
    end function H2_IJ_element

    function H2_CIp_element( &
        Ftot, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket, F1_bra, F1_ket) &
        result(H_CIp)

        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: Ftot, F1_bra, F1_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        REAL(rk) :: H_CIp, reduced_A1

        ! Off-diagonal nuclear spin-rotation for nucleus 2

        reduced_A1 = rank1_reduced_matrix_element_A1( &
            Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
            J_bra, J_ket, Omega_bra, Omega_ket)

        H_CIp = get_sign(1.0_rk + J_bra + I1 + I2 + 2.0_rk*F1_ket + Ftot) &
            * Wigner6j(J_ket, F1_ket, I1, F1_bra, J_bra, 1.0_rk) &
            * Wigner6j(I2, F1_ket, Ftot, F1_bra, I2, 1.0_rk) &
            * sqrt((2.0_rk*F1_bra + 1.0_rk)*(2.0_rk*F1_ket + 1.0_rk)) &
            * sqrt(I2*(I2 + 1.0_rk)*(2.0_rk*I2 + 1.0_rk)) &
            * reduced_A1 &
            * hfcc2(6)%field(index_field)%matelem(index_v_bra, index_v_ket)

    end function H2_CIp_element

    function H2_EQQ_element( &
        Ftot, index_hfcc2, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket, F1_bra, F1_ket) &
        result(H_EQ)

        implicit none

        INTEGER(ik), INTENT(IN) :: index_hfcc2, index_field, index_v_bra, index_v_ket, &
                                   Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: Ftot, F1_bra, F1_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        REAL(rk) :: H_EQ, normalization
        INTEGER(ik) :: q, q_min, q_max, q_step

        ! Nucleus-2 electric quadrupole interaction in the coupling scheme
        ! F1 = J + I1, F = F1 + I2. 
        H_EQ = 0.0_rk
        if (I2 < 1.0_rk) return
        if (abs(J_bra - J_ket) > 2.0_rk .or. abs(F1_bra - F1_ket) > 2.0_rk) return

        if (index_hfcc2 == 7) then
            q_min = 0
            q_max = 0
            q_step = 1
            normalization = 4.0_rk
        else if (index_hfcc2 == 8) then
            q_min = -2
            q_max = 2
            q_step = 4
            normalization = 4.0_rk*sqrt(6.0_rk) ! Brown--Carrington convention for eQq2.
        end if

        do q = q_min, q_max, q_step
            if ( (Lambda_ket == Lambda_bra - q) &
                .and. (nint(S_bra - S_ket) == 0) &
                .and. (nint(Sigma_bra - Sigma_ket) == 0) &
                .and. (nint(Omega_bra - Omega_ket) == q) ) then
                H_EQ = H_EQ &
                    + get_sign(J_bra + I1 + I2 + 2.0_rk*F1_ket + Ftot + 2.0_rk) &
                    * get_sign(J_bra - Omega_bra) &
                    * Wigner6j(J_ket, F1_ket, I1, F1_bra, J_bra, 2.0_rk) &
                    * Wigner6j(I2, F1_ket, Ftot, F1_bra, I2, 2.0_rk) &
                    * sqrt((2.0_rk*F1_bra + 1.0_rk)*(2.0_rk*F1_ket + 1.0_rk)) &
                    * sqrt((2.0_rk*J_bra + 1.0_rk)*(2.0_rk*J_ket + 1.0_rk)) &
                    * Wigner3j(J_bra, 2.0_rk, J_ket, -Omega_bra, real(q, rk), Omega_ket) &
                    / Wigner3j(I2, 2.0_rk, I2, -I2, 0.0_rk, I2) &
                    * hfcc2(index_hfcc2)%field(index_field)%matelem(index_v_bra, index_v_ket) &
                    / normalization
            end if
        end do
    end function H2_EQQ_element

    function H12_II_element( &
        Ftot, F1_bra, F1_ket, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(H_II)
        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: Ftot, F1_bra, F1_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket

        REAL(rk) :: H_II

        H_II = 0.0_rk

        if ( (Lambda_bra == Lambda_ket) &
            .and. (nint(S_bra - S_ket) == 0) &
            .and. (nint(Sigma_bra - Sigma_ket) == 0) ) then
            H_II = get_sign(F1_ket + I2 + Ftot + J_bra - Omega_bra) &
                * sqrt(30.0_rk &
                    * (2.0_rk*F1_bra + 1.0_rk) * (2.0_rk*F1_ket + 1.0_rk) &
                    * (2.0_rk*J_bra + 1.0_rk) * (2.0_rk*J_ket + 1.0_rk)) &
                * sqrt(I1*(I1 + 1.0_rk)*(2.0_rk*I1 + 1.0_rk) &
                    * I2*(I2 + 1.0_rk)*(2.0_rk*I2 + 1.0_rk)) &
                * Wigner3j(J_bra, 2.0_rk, J_ket, -Omega_bra, 0.0_rk, Omega_ket) &
                * Wigner6j(I2, F1_ket, Ftot, F1_bra, I2, 1.0_rk) &
                * Wigner9j(J_bra, J_ket, 2.0_rk, I1, I1, 1.0_rk, F1_bra, F1_ket, 1.0_rk) &
                * hfcc12(1)%field(index_field)%matelem(index_v_bra, index_v_ket)
        end if
    end function H12_II_element

    function H12_IIS_element( &
        Ftot, F1_bra, F1_ket, index_field, index_v_bra, index_v_ket, &
        Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) result(H_IIS)
        implicit none

        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket, Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: Ftot, F1_bra, F1_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        REAL(rk) :: H_IIS

        H_IIS = 0.0_rk
        ! c4(R) I1.I2 in the (J I1)F1, I2 -> F basis; F1 need not be conserved.
        ! Angular factor from Brown & Carrington, Eq. (8.295)
        if ( (Lambda_bra == Lambda_ket) &
            .and. (nint(S_bra - S_ket) == 0) &
            .and. (nint(Sigma_bra - Sigma_ket) == 0) &
            .and. (nint(J_bra - J_ket) == 0) &
            .and. (nint(Omega_bra - Omega_ket) == 0) ) then
            H_IIS = get_sign(F1_ket + Ftot + I2) &
                * get_sign(F1_bra + J_ket + 1.0_rk + I1) &
                * sqrt((2.0_rk*F1_bra + 1.0_rk)*(2.0_rk*F1_ket + 1.0_rk)) &
                * Wigner6j(I2, F1_ket, Ftot, F1_bra, I2, 1.0_rk) &
                * Wigner6j(I1, F1_ket, J_bra, F1_bra, I1, 1.0_rk) &
                * sqrt(I1*(I1 + 1.0_rk)*(2.0_rk*I1 + 1.0_rk) &
                    * I2*(I2 + 1.0_rk)*(2.0_rk*I2 + 1.0_rk)) &
                * hfcc12(2)%field(index_field)%matelem(index_v_bra, index_v_ket)
        end if
    end function H12_IIS_element

    function rank1_reduced_matrix_element_A1( &
        Lambda_bra, Lambda_ket, S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(reduced_A1)

        implicit none

        INTEGER(ik), INTENT(IN) :: Lambda_bra, Lambda_ket
        REAL(rk), INTENT(IN) :: S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        INTEGER(ik) :: q, double_omega, double_omega_min, double_omega_max
        REAL(rk) :: reduced_A1, q_real, Omega_intermediate, &
                    term_J_bra, term_J_ket

        ! For off-diagonal nuclear spin-rotation; excludes the C_I' radial integral and nuclear-spin coupling factors
        ! see Eq. (A10) from Brown, Kaise, Kerr and Milton, Molecular Physics 36, 553-582 (1978)
        ! q retains the molecular-fixed component label of the original operator 
        
        reduced_A1 = 0.0_rk
        if (nint(S_bra - S_ket) /= 0) return
        if (nint(Sigma_bra - Sigma_ket) /= 0) return
        if (abs(J_bra - J_ket) > 1.0_rk + small_) return

        do q = -1, 1, 2
            if (Lambda_bra /= Lambda_ket - 2*q) cycle
            if (nint(Omega_bra - Omega_ket) /= -2*q) cycle
            q_real = real(q, rk)

            term_J_ket = 0.0_rk
            double_omega_min = -nint(2.0_rk*J_ket)
            double_omega_max =  nint(2.0_rk*J_ket)
            do double_omega = double_omega_min, double_omega_max, 2
                Omega_intermediate = 0.5_rk*real(double_omega, rk)
                term_J_ket = term_J_ket &
                    + sqrt(J_ket*(J_ket + 1.0_rk)*(2.0_rk*J_ket + 1.0_rk)) &
                    * get_sign(J_ket - Omega_intermediate) &
                    * Wigner3j(J_bra, 1.0_rk, J_ket, &
                               -Omega_bra, -q_real, Omega_intermediate) &
                    * Wigner3j(J_ket, 1.0_rk, J_ket, &
                               -Omega_intermediate, -q_real, Omega_ket)
            end do

            term_J_bra = 0.0_rk
            double_omega_min = -nint(2.0_rk*J_bra)
            double_omega_max =  nint(2.0_rk*J_bra)
            do double_omega = double_omega_min, double_omega_max, 2
                Omega_intermediate = 0.5_rk*real(double_omega, rk)
                term_J_bra = term_J_bra &
                    + sqrt(J_bra*(J_bra + 1.0_rk)*(2.0_rk*J_bra + 1.0_rk)) &
                    * get_sign(J_bra - Omega_intermediate) &
                    * Wigner3j(J_bra, 1.0_rk, J_bra, &
                               -Omega_bra, -q_real, Omega_intermediate) &
                    * Wigner3j(J_bra, 1.0_rk, J_ket, &
                               -Omega_intermediate, -q_real, Omega_ket)
            end do

            reduced_A1 = reduced_A1 - 0.5_rk &
                * get_sign(J_bra - Omega_bra) &
                * sqrt((2.0_rk*J_bra + 1.0_rk)*(2.0_rk*J_ket + 1.0_rk)) &
                * (term_J_ket + term_J_bra)
        end do

    end function rank1_reduced_matrix_element_A1

    function Wigner9j(j11, j12, j13, j21, j22, j23, j31, j32, j33) result(w9j)
        implicit none
        REAL(rk), INTENT(IN) :: j11, j12, j13, j21, j22, j23, j31, j32, j33
        REAL(rk) :: w9j
        REAL(rk) :: jmin, jmax, j_, term
        INTEGER(ik) :: double_j_min, double_j_max, double_j

        ! Intersection of triangle ranges that involve the summation variable j:
        ! (j11, j33, j), (j21, j32, j), (j12, j23, j)
        ! see (A17) of Kato's paper 
        jmin = max( abs(j11 - j33), abs(j21 - j32), abs(j12 - j23) )
        jmax = min( j11 + j33, j21 + j32, j12 + j23)

        if (jmax < jmin - 1.0e-10_rk) then
            w9j = 0.0_rk
            return
        end if

        ! Loop over 2*j as an integer so half-integers are included
        double_j_min = ceiling(2.0_rk*jmin - 100.0_rk*small_)
        double_j_max = floor(2.0_rk*jmax + 100.0_rk*small_)

        w9j = 0.0_rk
        do double_j = double_j_min, double_j_max
            j_ = 0.5_rk * real(double_j, rk)
            term = get_sign( real(double_j, rk) ) * (real(double_j, rk) + 1.0_rk) &
                * Wigner6j(j11, j21, j31,  j32, j33, j_) &
                * Wigner6j(j12, j22, j32,  j21, j_, j23) &
                * Wigner6j(j13, j23, j33,  j_, j11, j12)
            w9j = w9j + term
        end do
    end function Wigner9j

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
        upper_limit = min(k1, k2, k3) 

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
        real(rk) log_gamma
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
                    tri = log_gamma(a + b - c + 1.0_rk) &
                        + log_gamma(a + c - b + 1.0_rk) &
                        + log_gamma(b + c - a + 1.0_rk) &
                        - log_gamma(a + b + c + 2.0_rk)
                    tri = exp(tri)
                end if
            end do
        end function triangle_coefficient

    end function Wigner6j

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
