module F1_intensity

    use F1_hyperfine
    use accuracy
    use diatom_module, only : Ndipoles, dipoletm

    implicit none

    REAL(rk), ALLOCATABLE :: primitive_F1_reduced_TDM_matrix(:,:)
    INTEGER(ik), PARAMETER :: unit_hyperfine_transitions = 66
contains

    subroutine F1_hyperfine_intensity
        implicit none
        INTEGER(ik) :: index_F1_bra, index_F1_ket, &
                       index_represCs_bra, index_represCs_ket
        INTEGER(ik) :: Ndimen_F1_bra, Ndimen_F1_ket, &
                       Nlevels_F1_bra, Nlevels_F1_ket
        REAL(rk), ALLOCATABLE :: parity_conserved_F1_reduced_TDM_matrix(:,:), &
                                 parity_conserved_F1_transitions_matrix(:,:)
        CHARACTER(LEN=100) :: FMT1 =  "(A22,&
                                        A22,&
                                        A22,&
                                        A8,&
                                        A22,&
                                        A6,&
                                        A6,&
                                        A6,&
                                        A7,&
                                        A6,&
                                        A6,&
                                        A4,&
                                        A7,&
                                        A6,&
                                        A6,&
                                        A8,&
                                        A6,&
                                        A6,&
                                        A7,&
                                        A6,&
                                        A6,&
                                        A4,&
                                        A7,&
                                        A6,&
                                        A6)"


        CALL F1_hyperfine_structrure

        open(unit=unit_hyperfine_transitions, file="hyperfine.trans")

        ! write(unit_hyperfine_transitions, "(A12, A12, A22, A22, A22)") &
        !     "N_upper", "N_lower", "Einstein-A [s-1]", "nu [cm-1]", "S [Debye^2]"

        ! write(unit_hyperfine_transitions, FMT1) &
        !     "nu", &
        !     "Strength", &
        !     "E", &
        !     " --> ", &
        !     "E", &
        !     " :: ", &
        !     "F", &
        !     "I", &
        !     "parity", &
        !     "J", &
        !     "state", &
        !     "v", &
        !     "Lambda", &
        !     "Sigma", &
        !     "Omega", &
        !     " --> ", &
        !     "F", &
        !     "I", &
        !     "parity", &
        !     "J", &
        !     "state", &
        !     "v", &
        !     "Lambda", &
        !     "Sigma", &
        !     "Omega" 

        do index_F1_bra = 1, num_F1
            Ndimen_F1_bra = primitive_F1_basis(index_F1_bra)%Ndimen  

            do index_F1_ket = max(1, index_F1_bra - 1), min(num_F1, index_F1_bra + 1) ! F1 selection rule
                if ((index_F1_bra == index_F1_ket) .and. F1_list(index_F1_bra) == 0) cycle
                Ndimen_F1_ket = primitive_F1_basis(index_F1_ket)%Ndimen

                ALLOCATE(primitive_F1_reduced_TDM_matrix(Ndimen_F1_bra, Ndimen_F1_ket))

                CALL construct_primitive_F1_reduced_TDM_matrix( &
                        index_F1_bra, index_F1_ket)   

                do index_represCs_bra = 1, sym%NrepresCs
                    ! parity selection rule
                    index_represCs_ket = mod(index_represCs_bra, sym%NrepresCs) + 1 

                    Nlevels_F1_bra = eigen_all_F1(index_F1_bra, index_represCs_bra)%Nlevels
                    Nlevels_F1_ket = eigen_all_F1(index_F1_ket, index_represCs_ket)%Nlevels

                    ALLOCATE(parity_conserved_F1_reduced_TDM_matrix(Nlevels_F1_bra, Nlevels_F1_ket))
                    
                    CALL construct_parity_conserved_F1_reduced_TDM_matrix ! TDM: Transition dipole moment

                    ALLOCATE(parity_conserved_F1_transitions_matrix(Nlevels_F1_bra, Nlevels_F1_ket))
                    ! wavenumber_to_MHz = wavenumber_to_MHz
                    CALL construct_parity_conserved_F1_transitions_matrix
                    
                    DEALLOCATE(parity_conserved_F1_transitions_matrix)
                    DEALLOCATE(parity_conserved_F1_reduced_TDM_matrix)
                end do           
                
                DEALLOCATE(primitive_F1_reduced_TDM_matrix)
            end do            
        end do

        close(unit=unit_hyperfine_transitions)

    contains
        subroutine construct_parity_conserved_F1_reduced_TDM_matrix
        ! | <psi_m^{tau,F}| T^1(mu) |psi_n^{tau',F'}> | ** 2
        implicit none

        REAL(rk), ALLOCATABLE :: intermediate_matrix(:, :)
            ! ! ^{tau,F}M^{tau', F'} = Psi^{tau, F} * ^{F}M^{F'} * Psi^{tau', F'} 
            ! parity_conserved_F1_reduced_TDM_matrix = &
            !     matmul(matmul(transpose(eigen_all_F1(index_F1_bra, index_represCs_bra)%vect), &
            !                 primitive_F1_reduced_TDM_matrix), &
            !         eigen_all_F1(index_F1_ket, index_represCs_ket)%vect)

            ALLOCATE(intermediate_matrix(Nlevels_F1_bra, Ndimen_F1_ket))

            ! M_indermediate} = Psi^{tau, F} * ^{F}M^{F'}
            CALL dgemm('T', 'N', Nlevels_F1_bra, Ndimen_F1_ket, Ndimen_F1_bra, 1.0_rk, &
                        eigen_all_F1(index_F1_bra, index_represCs_bra)%vect, Ndimen_F1_bra, &
                        primitive_F1_reduced_TDM_matrix, Ndimen_F1_bra, &
                        0.0_rk, intermediate_matrix, Nlevels_F1_bra)

            ! ^{tau, F}M^{tau', F'} = M_indermediate * Psi^{tau', F'}           
            CALL dgemm('N', 'N', Nlevels_F1_bra, Nlevels_F1_ket, Ndimen_F1_ket, 1.0_rk, &
                        intermediate_matrix, Nlevels_F1_bra, &
                        eigen_all_F1(index_F1_ket, index_represCs_ket)%vect, Ndimen_F1_ket, &
                        0.0_rk, parity_conserved_F1_reduced_TDM_matrix, Nlevels_F1_bra)

            DEALLOCATE(intermediate_matrix)

            parity_conserved_F1_reduced_TDM_matrix = &
                parity_conserved_F1_reduced_TDM_matrix ** 2
        end subroutine construct_parity_conserved_F1_reduced_TDM_matrix
        
        subroutine construct_parity_conserved_F1_transitions_matrix
            implicit none
            INTEGER(ik) :: bra, ket, maxloc_ket, maxloc_bra
            TYPE(quantaT), POINTER :: quanta_bra, quanta_ket   
            REAL(rk)  :: EinsteinA
            
            CHARACTER(LEN=120) :: FMT2 =  "(F22.12,&
                                            E22.12,&
                                            F22.12,&
                                            A8,&
                                            F22.12,&
                                            A6,&
                                            F6.1,&
                                            F6.1,&
                                            A7,&
                                            F6.1,&
                                            I6,& 
                                            I4,&
                                            I7,&
                                            F6.1,&
                                            F6.1,&
                                            A8,&
                                            F6.1,&
                                            F6.1,&
                                            A7,&
                                            F6.1,&
                                            I6,& 
                                            I4,&
                                            I7,&
                                            F6.1,&
                                            F6.1)"

            do ket = 1, Nlevels_F1_ket                
                ! maxloc_ket = MAXLOC(ABS(eigen_all_F1(index_F1_ket, index_represCs_ket)%vect(:, ket)), DIM=1) 
                ! quanta_ket => primitive_F1_basis(index_F1_ket)%icontr(maxloc_ket)
                do bra = 1, Nlevels_F1_bra
                    ! maxloc_bra = MAXLOC(ABS(eigen_all_F1(index_F1_bra, index_represCs_bra)%vect(:, bra)), DIM=1) 
                    ! quanta_bra => primitive_F1_basis(index_F1_bra)%icontr(maxloc_bra)

                    parity_conserved_F1_transitions_matrix(bra, ket) = &
                        eigen_all_F1(index_F1_bra, index_represCs_bra)%val(bra) &
                            - eigen_all_F1(index_F1_ket, index_represCs_ket)%val(ket)
                    
                    if (parity_conserved_F1_transitions_matrix(bra, ket) < 0) cycle

                    ! Einstein A coefficients: see Eq.(23) of Western 
                    ! doi:10.1016/j.jqsrt.2016.04.010
                    ! g_u = 2 * F'+1
                    EinsteinA = 3.13618872E-7_rk &
                        * parity_conserved_F1_transitions_matrix(bra, ket) &
                        * parity_conserved_F1_transitions_matrix(bra, ket) &
                        * parity_conserved_F1_transitions_matrix(bra, ket) &
                        * parity_conserved_F1_reduced_TDM_matrix(bra, ket) &
                        / (2.0_rk * primitive_F1_basis(index_F1_bra)%icontr(1)%F1 + 1.0_rk)

                    WRITE(unit_hyperfine_transitions, "(I12, I12, E22.12, F22.12, E22.12)") &
                        eigen_all_F1(index_F1_bra, index_represCs_bra)%quanta(bra)%iroot, &
                        eigen_all_F1(index_F1_ket, index_represCs_ket)%quanta(ket)%iroot, &
                        EinsteinA, &
                        parity_conserved_F1_transitions_matrix(bra, ket), &
                        parity_conserved_F1_reduced_TDM_matrix(bra, ket)

                    ! WRITE(unit_hyperfine_transitions, FMT2) &
                    !     parity_conserved_F1_transitions_matrix(bra, ket), &
                    !     parity_conserved_F1_reduced_TDM_matrix(bra, ket), &
                    !     eigen_all_F1(index_F1_bra, index_represCs_bra)%val(bra), &
                    !     " --> ", &
                    !     eigen_all_F1(index_F1_ket, index_represCs_ket)%val(ket), &
                    !     " :: ", &
                    !     quanta_bra%F1, &
                    !     quanta_bra%I1, &
                    !     trim(parity_sign(index_represCs_bra)), &
                    !     quanta_bra%Jrot, &
                    !     quanta_bra%istate, &
                    !     quanta_bra%v, &
                    !     quanta_bra%ilambda, &
                    !     quanta_bra%sigma, &
                    !     quanta_bra%omega, &
                    !     " --> ", &
                    !     quanta_ket%F1, &
                    !     quanta_ket%I1, &
                    !     trim(parity_sign(index_represCs_ket)), &
                    !     quanta_ket%Jrot, &
                    !     quanta_ket%istate, &
                    !     quanta_ket%v, &
                    !     quanta_ket%ilambda, &
                    !     quanta_ket%sigma, &
                    !     quanta_ket%omega
                end do    
            end do 
        end subroutine construct_parity_conserved_F1_transitions_matrix

        function parity_sign(index_represCs)
            INTEGER(ik), INTENT(IN) :: index_represCs
            CHARACTER(4) :: parity_sign
            if ( index_represCs == 1 ) then
                parity_sign = '+'
            elseif ( index_represCs == 2 ) then
                parity_sign = '-'
            endif
        end function parity_sign
    end subroutine F1_hyperfine_intensity


    subroutine construct_primitive_F1_reduced_TDM_matrix( &
        index_F1_bra, index_F1_ket)

        implicit none
        
        INTEGER(ik), INTENT(IN) :: index_F1_bra, index_F1_ket
        INTEGER(ik) :: Ndimen_F1_bra, Ndimen_F1_ket, bra, ket
        REAL(rk) :: F1_bra, F1_ket, &
            S_bra, S_ket, Sigma_bra, Sigma_ket, &
            J_bra, J_ket, Omega_bra, Omega_ket
        INTEGER(ik) :: state_bra, state_ket, index_v_bra, index_v_ket, &
                        v_bra, v_ket, Lambda_bra, Lambda_ket, index_field
        

        primitive_F1_reduced_TDM_matrix = 0.0_rk

        Ndimen_F1_bra = primitive_F1_basis(index_F1_bra)%Ndimen
        Ndimen_F1_ket = primitive_F1_basis(index_F1_ket)%Ndimen

        do ket = 1, Ndimen_F1_ket
            call get_quanta(index_F1_ket, ket, &
                            F1_ket, state_ket, index_v_ket, v_ket, Lambda_ket, &
                            S_ket, Sigma_ket, J_ket, Omega_ket)

            do bra = 1, Ndimen_F1_bra
                call get_quanta(index_F1_bra, bra, &
                                F1_bra, state_bra, index_v_bra, v_bra, Lambda_bra, &
                                S_bra, Sigma_bra, J_bra, Omega_bra)

                do index_field = 1, Ndipoles
                    if ( (state_bra == dipoletm(index_field)%istate) &
                        .and. (state_ket == dipoletm(index_field)%jstate)) then
                            
                        primitive_F1_reduced_TDM_matrix(bra, ket) = &
                            primitive_F1_reduced_TDM_matrix(bra, ket)  &
                            + primitive_F1_reduced_TDM_matrix_element( &
                                F1_bra, F1_ket, &
                                index_field, index_v_bra, index_v_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket)
                    end if
                end do
            end do
        end do
    end subroutine construct_primitive_F1_reduced_TDM_matrix

    function primitive_F1_reduced_TDM_matrix_element( &
        F1_bra, F1_ket, &
        index_field, index_v_bra, index_v_ket, &
        S_bra, S_ket, Sigma_bra, Sigma_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(tdm)

        implicit none

        REAL(rk), INTENT(IN) :: F1_bra, F1_ket, &
                                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        INTEGER(ik), INTENT(IN) :: index_field, index_v_bra, index_v_ket
        INTEGER(ik) :: t, q
        REAL(rk) :: tdm

        ! Kato
        ! if ( (nint(S_ket - S_bra) == 0 ) &
        !     .and.  (nint(Sigma_bra - Sigma_ket) == 0) ) then
        !     do t = -1, 1
        !         if ( nint(Omega_bra - Omega_ket - t) == 0 ) then
        !             tdm = get_sign(-F1_bra+I1-Omega_ket) &
        !                 * sqrt((2.0_rk * F1_bra + 1.0_rk) &
        !                         * (2.0_rk * F1_ket + 1.0_rk) &
        !                         * (2.0_rk * J_bra + 1.0_rk) &
        !                         * (2.0_rk * J_ket + 1.0_rk) &
        !                         ) &
        !                 * Wigner6j(F1_bra, F1_ket, 1.0_rk, J_ket, J_bra, I1) &
        !                 * get_sign(-real(t,rk)) &
        !                 * mu_v() &
        !                 * Wigner3j(J_bra, J_ket, 1.0_rk, &
        !                             Omega_bra, -Omega_ket, -real(t,rk))
        !             exit
        !         end if
        !     end do
        ! else
        !     tdm = 0
        ! end if

        tdm = 0
        if ( (nint(S_ket - S_bra) == 0 ) &
            .and.  (nint(Sigma_bra - Sigma_ket) == 0) ) then
            do q = -1, 1
                tdm = tdm + &
                    get_sign(J_bra+I1+F1_ket+1.0_rk) &
                    * sqrt((2.0_rk * F1_bra + 1.0_rk) &
                        * (2.0_rk * F1_ket + 1.0_rk) &
                        * (2.0_rk * J_bra + 1.0_rk) &
                        * (2.0_rk * J_ket + 1.0_rk) &
                        ) &
                    * Wigner6j(J_bra, F1_bra, I1, F1_ket, J_ket, 1.0_rk)&
                    * get_sign(J_bra-Omega_bra)&
                    * Wigner3j(J_bra, 1.0_rk, J_ket, -Omega_bra, real(q,rk), Omega_ket)&
                    * dipoletm(index_field)%matelem(index_v_bra, index_v_ket)
            end do
        end if
    end function primitive_F1_reduced_TDM_matrix_element
    
end module F1_intensity