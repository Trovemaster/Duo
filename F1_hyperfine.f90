module F1_hyperfine

    use accuracy
    use lapack
    use symmetry, only : sym
    use diatom_module, only: basis, three_j, faclog, eigen, jmax, job,&
                            vibrational_totalroots, vibrational_contrfunc, vibrational_quantum_number, &
                            hfcc1, F1_hyperfine_setup, GLOBAL_NUM_HFCC_OBJECT, poten, &
                            eigenT, basisT, quantaT, fieldT, linkT

    implicit none
    
    
    TYPE(basisT), ALLOCATABLE :: primitive_F1_basis(:)
    
    REAL(rk), ALLOCATABLE :: F1_list(:)
    
    INTEGER(ik) :: num_F1, num_primitive_F1_basis, num_represCs
    REAL(rk) :: F1_global_min, F1_global_max, I1
    INTEGER(ik) :: alloc
    REAL(rk) :: J_global_min = 0.5_rk, J_global_max
    REAL(rk) :: J_min, J_max, MHz_to_wavenumber = 1.0E6_rk / vellgt, &
                wavenumber_to_MHz = vellgt / 1.0E6_rk

    TYPE(eigenT), allocatable :: eigen_all_F1(:,:)
    REAL(rk) :: min_eigen_value_F1
    LOGICAL :: F1_hyperfine_gridvalue_allocated  = .false.
    INTEGER(ik), PARAMETER :: unit_hyperfine_states = 65
    
contains

    subroutine F1_hyperfine_structrure(iverbose)
        ! Top level 
        implicit none
        integer(ik), INTENT(IN) :: iverbose
        INTEGER(ik) :: index_F1, index_represCs, iroot, iF1_ID
        INTEGER(ik) :: Ndimen_F1, Nlevels_F1, index_levels_F1, maxlocation
        REAL(rk), ALLOCATABLE :: eigen_value_F1(:), &
                                    parity_conserved_F1_matrix(:, :), &
                                    primitive_F1_hyperfine_matrix(:,:), &
                                    transformation_matrix(:,:), &
                                    eigen_vector_F1(:, :)
        TYPE(quantaT), POINTER :: quanta
        
        open(unit=unit_hyperfine_states, file="hyperfine.states")
        
        write(out, '(/A)') "Start: hyperfine energies calculation"

        I1 = F1_hyperfine_setup%I1
 
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
        
        CALL construct_F1_hyperfine_constant_field_matrix(iverbose)      
        

        ALLOCATE(F1_list(num_F1))
        ALLOCATE(primitive_F1_basis(num_F1))

        CALL construct_primitive_F1_basis(iverbose)

        if (iverbose>=4) then
            write(out, '(/A2, A)') '', 'Start: hyperfine states'
            write(out, '(/A2, A, F22.12)') '', 'Zero point energy =', job%ZPE
        endif

        ALLOCATE(eigen_all_F1(num_F1, sym%NrepresCs))

        ! write(unit_hyperfine_states, '(A7, A22, A7, A7, A7, A7, A7, A12, A7, A7, A7, A7)') &
        !     'Number', 'Energy [cm-1]', 'g', 'F', 'I', 'parity', 'J', 'state', 'v', 'Lambda', 'Sigma', 'Omega'

        iroot = 0
        min_eigen_value_F1 = safe_max
        do index_F1 = 1, num_F1
            if (iverbose>=4) then
                write(out, '(/A4, A, F6.1)') '', 'Eigen states of F =', F1_list(index_F1)
                write(out, '(/A6, A)') '', 'Construct primitive hyperfine matrix'
            endif

            Ndimen_F1 = primitive_F1_basis(index_F1)%Ndimen
            ALLOCATE(primitive_F1_hyperfine_matrix(Ndimen_F1, Ndimen_F1))
            CALL construct_primitive_F1_hyperfine_matrix( &
                index_F1, primitive_F1_hyperfine_matrix)

            if (iverbose>=4) write(out, '(A6, A/)') '', '... done'

            do index_represCs = 1, sym%NrepresCs
                if (iverbose>=4) then
                    write(out, '(/A6, A, A3)') '', 'Eigen states of parity', trim(parity_sign(index_represCs))
                endif

                Nlevels_F1 = Ndimen_F1 / sym%NrepresCs

                ALLOCATE(eigen_value_F1(Nlevels_F1))
                ALLOCATE(parity_conserved_F1_matrix(Nlevels_F1, Nlevels_F1))
                ALLOCATE(transformation_matrix(Ndimen_F1, Nlevels_F1))
                ALLOCATE(eigen_vector_F1(Ndimen_F1, Nlevels_F1))

                ALLOCATE(eigen_all_F1(index_F1, index_represCs)%quanta(Nlevels_F1))
                ALLOCATE(eigen_all_F1(index_F1, index_represCs)%val(Nlevels_F1))
                ALLOCATE(eigen_all_F1(index_F1, index_represCs)%vect(Ndimen_F1, Nlevels_F1))

                if (iverbose>=4) then
                    write(out, '(/A8, A)') &
                        '', 'Construct parity conserved hyperfine matrix'
                endif

                CALL construct_parity_conserved_F1_matrix( &
                    index_F1, index_represCs, primitive_F1_hyperfine_matrix, &
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

                eigen_all_F1(index_F1, index_represCs)%Nlevels = Nlevels_F1
                eigen_all_F1(index_F1, index_represCs)%Ndimen = Ndimen_F1
                eigen_all_F1(index_F1, index_represCs)%val = eigen_value_F1
                eigen_all_F1(index_F1, index_represCs)%vect = eigen_vector_F1
                
                if (iverbose>=4) then
                    write(out, '(/A8, A)') '', 'Calculate energy levels'
                endif
                eigen_value_F1(:) = eigen_value_F1(:) - job%ZPE

                iF1_ID = 0
                do index_levels_F1 = 1, Nlevels_F1
                    maxlocation = MAXLOC(ABS(eigen_vector_F1(:, index_levels_F1)), DIM=1)
                    quanta => primitive_F1_basis(index_F1)%icontr(maxlocation) 

                    iroot = iroot + 1
                    iF1_ID = iF1_ID + 1
                    eigen_all_F1(index_F1, index_represCs)%quanta(index_levels_F1)%iroot = iroot
                    eigen_all_F1(index_F1, index_represCs)%quanta(index_levels_F1)%iF1_ID = iF1_ID

                    min_eigen_value_F1 = min(&
                        eigen_all_F1(index_F1, index_represCs)%val(index_levels_F1), &
                        min_eigen_value_F1) 

                    if (iverbose>=4) then
                        write(unit_hyperfine_states, '(I7, F22.12, I7, F7.1, F7.1, A7, F7.1, A12, I7, I7, F7.1, F7.1)') &
                            iroot, &
                            eigen_value_F1(index_levels_F1), &
                            int(quanta%F1 * 2.0_rk + 1.0_rk), &
                            quanta%F1, &
                            quanta%I1, &
                            trim(parity_sign(index_represCs)), &
                            quanta%Jrot, &
                            trim(poten(quanta%istate)%name), &
                            quanta%v, &
                            quanta%ilambda, &
                            quanta%sigma, &
                            quanta%omega
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
    end subroutine F1_hyperfine_structrure

    subroutine construct_F1_hyperfine_constant_field_matrix(iverbose)
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
        end do

        write(out, '(A2, A/)') '', "End: hyperfine constants under the vibrational bases"
    end subroutine construct_F1_hyperfine_constant_field_matrix

    subroutine construct_primitive_F1_basis(iverbose)
        ! |state, v; Lambda; S, Sigma; J, Omega; I1; F1>
        implicit none

        integer(ik), INTENT(IN) :: iverbose
        REAL(rk) :: F1, J
        INTEGER(ik) :: index_F1, Ndimen_F1
        INTEGER(ik) :: index_J, index_J_min, index_J_max, start_contr_F1, end_contr_F1, index_contr_F1
        TYPE(quantaT), POINTER :: contr_F1
        
        write(out, '(/A2, A)') '', "Start: contracted Hund's case (a_beta) basis"

        num_primitive_F1_basis = 0
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
            if (iverbose>=4) then
                write(out, '(A6, A7, A7, A7, A7, A7, A7, A7, A7, A7, A7)') &
                    '    ', 'No.', 'F', 'I', 'state', 'v', 'Lambda', 'J', 'Omega', 'S', 'Sigma'
            endif
            do index_contr_F1 = 1, Ndimen_F1
                contr_F1 => primitive_F1_basis(index_F1)%icontr(index_contr_F1)
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

        ! H_hfs^tau = M_indermediate * Phi           
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
        INTEGER(ik) :: Ndimen_F1, index_hfcc1, index_field

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

                    do index_hfcc1 = 1, GLOBAL_NUM_HFCC_OBJECT
                        do index_field = 1, hfcc1(index_hfcc1)%num_field
                            if ( (state_bra == hfcc1(index_hfcc1)%field(index_field)%istate) &
                                .and. (state_ket == hfcc1(index_hfcc1)%field(index_field)%jstate)) then

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
            H_hfs = H_EQQ0_element( &        
                        F1, index_field, index_v_bra, index_v_ket,  &
                        Lambda_bra, Lambda_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket)
        else if (index_hfcc1 == 7) then
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
                    * hfcc1(1)%field(index_field)%matelem(index_v_bra, index_v_ket)
                    
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
                    * hfcc1(2)%field(index_field)%matelem(index_v_bra, index_v_ket)
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
                    * hfcc1(3)%field(index_field)%matelem(index_v_bra, index_v_ket) / (-3.0_rk)
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
                        * hfcc1(4)%field(index_field)%matelem(index_v_bra, index_v_ket) / sqrt(6.0_rk)
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
                * hfcc1(5)%field(index_field)%matelem(index_v_bra, index_v_ket)
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
                * hfcc1(6)%field(index_field)%matelem(index_v_bra, index_v_ket) / (-2.0_rk)
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
                    * hfcc1(7)%field(index_field)%matelem(index_v_bra, index_v_ket) / (-sqrt(24.0_rk))
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
