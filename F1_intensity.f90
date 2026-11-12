module F1_intensity

    use accuracy
    use symmetry, only : sym, correlate_to_Cs
    use diatom_module, only : Ndipoles, dipoletm, intensity, job, fieldT, poten, &
                            hyperfine_case, nQuadrupoles, quadrupoletm, action
    use F1_hyperfine, only : primitive_hfs_basis, &
                            F1_list, Ftot_list, eigen_all_hfs, min_eigen_value_hfs, &
                            get_sign, I1, I2, Wigner3j, Wigner6j, num_F1, num_Ftot, get_quanta

    implicit none

    INTEGER(ik), PARAMETER :: unit_hyperfine_transitions = 66
    real(rk) :: A_coef_s_1, A_coef_s_2 ! related to Einstein A coefficient of E1, E2 transitions, respectively
contains

    subroutine hyperfine_transition ! control different nuclear-spin cases
        implicit none

        write(out, '(/A)') "Start: hyperfine transitions calculation"
        open(unit=unit_hyperfine_transitions, file="hyperfine.trans")

        select case(hyperfine_case)
        case(1)
            write(out, '(/A)') " Hyperfine case: single non-zero-spin nucleus."
            if (action%quadrupole) then
                stop 'E2 transition intensity is not supported for single non-zero-spin case.'
            elseif (action%magdipole) then
                stop 'M1 transition intensity is not supported for single non-zero-spin case.'
            else
                if (Ndipoles <= 0) stop 'Error: No dipole curve is given for the intensity calculation.'
                call F1_hyperfine_intensity
            endif
        case(2)
            write(out, '(/A)') " Hyperfine case: two non-zero-spin nuclei; homonuclear."
            if (action%quadrupole) then
                if (nQuadrupoles <= 0) stop 'Error: No quadrupole curve is given for E2 intensity calculation.'
                call IF_hyperfine_E2_intensity
            elseif (action%magdipole) then
                stop 'M1 transition intensity is not supported for homonuclear case.'
            else
                if (Ndipoles <= 0) stop 'Error: No dipole curve is given for the intensity calculation.'
                call IF_hyperfine_E1_intensity
            endif
        case(3)
            write(out, '(/A)') " Hyperfine case: heteronuclear two-spin coupling."
            if (action%quadrupole) then
                if (nQuadrupoles <= 0) stop 'Error: No quadrupole curve is given for E2 intensity calculation.'
                call F1F_hyperfine_E2_intensity
            elseif (action%magdipole) then
                stop 'M1 transition intensity is not supported for the heteronuclear two-spin case.'
            else
                if (Ndipoles <= 0) stop 'Error: No dipole curve is given for the intensity calculation.'
                call F1F_hyperfine_E1_intensity
            endif
        case default
            stop 'ERROR: hyperfine nuclear-spin case has not been initialized.'
        end select
    end subroutine hyperfine_transition 

    subroutine IF_hyperfine_E1_intensity
        implicit none 
        INTEGER(ik) :: index_Ftot_bra, index_Ftot_ket, &
                        index_represC2v_bra, index_represC2v_ket
        INTEGER(ik) :: Ndimen_Ftot_bra, Ndimen_Ftot_ket, &
                        Nlevels_Ftot_bra, Nlevels_Ftot_ket
        REAL(rk), ALLOCATABLE :: primitive_IF_reduced_TDM_matrix(:,:), &
                                parity_conserved_IF_reduced_TDM_matrix(:,:), &
                                parity_conserved_IF_transitions_matrix(:,:)

        A_coef_s_1 = 64.0E-36 * pi**4 / (3.0_rk * planck)

        do index_Ftot_bra = 1, num_Ftot
            write(out, '(/A4, A, F6.1)') '', 'Start: F_bra =', Ftot_list(index_Ftot_bra)

            Ndimen_Ftot_bra = primitive_hfs_basis(index_Ftot_bra)%Ndimen

            do index_Ftot_ket = max(1, index_Ftot_bra - 1), min(num_Ftot, index_Ftot_bra + 1)
                if ((index_Ftot_bra == index_Ftot_ket) .and. Ftot_list(index_Ftot_bra) == 0) cycle

                write(out, '(/A6, A, F6.1)') '', 'Start: F_ket =', Ftot_list(index_Ftot_ket)
                Ndimen_Ftot_ket = primitive_hfs_basis(index_Ftot_ket)%Ndimen

                write(out, '(/A8, A)') '', &
                    'Construct primitive reduced transition dipole moment matrix'
                ALLOCATE(primitive_IF_reduced_TDM_matrix(Ndimen_Ftot_bra, Ndimen_Ftot_ket))

                CALL construct_primitive_IF_reduced_TDM_matrix(index_Ftot_bra, index_Ftot_ket)

                write(out, '(A8, A/)') '', '... done'

                do index_represC2v_bra = 1, sym%Nrepresen

                    ! symmetry selection rule
                    index_represC2v_ket = E1_symmetry_selection(index_represC2v_bra)

                    Nlevels_Ftot_bra = eigen_all_hfs(index_Ftot_bra, index_represC2v_bra)%Nlevels
                    Nlevels_Ftot_ket = eigen_all_hfs(index_Ftot_ket, index_represC2v_ket)%Nlevels

                    write(out, '(/A10, A)') '',&
                        'Construct parity conserved reduced transition dipole moment matrix'

                    if (Nlevels_Ftot_bra <= 0 .or. Nlevels_Ftot_ket <= 0) then
                        write(out, '(A10, A/)') '', '... skipped: empty F-symmetry block'
                        cycle
                    end if

                    ALLOCATE(parity_conserved_IF_reduced_TDM_matrix(Nlevels_Ftot_bra, Nlevels_Ftot_ket))
                    CALL construct_parity_conserved_IF_reduced_TDM_matrix
                    write(out, '(A10, A/)') '', '... done'

                    write(out, '(/A10, A)') '', 'Calculate and print transitions'
                    ALLOCATE(parity_conserved_IF_transitions_matrix(Nlevels_Ftot_bra, Nlevels_Ftot_ket))
                    CALL calc_IF_EinsteinA_coefficient
                    write(out, '(A10, A/)') '', '... done'

                    DEALLOCATE(parity_conserved_IF_transitions_matrix)
                    DEALLOCATE(parity_conserved_IF_reduced_TDM_matrix)
                end do

                DEALLOCATE(primitive_IF_reduced_TDM_matrix)

                write(out, '(/A6, A, F6.1)') '', 'End: F_ket =', Ftot_list(index_Ftot_ket)
            end do

            write(out, '(/A4, A, F6.1)') '', 'End: F_bra =', Ftot_list(index_Ftot_bra)
        end do

        close(unit=unit_hyperfine_transitions)
        write(out, '(/A)') "End: hyperfine transitions calculation"

    contains
        subroutine construct_parity_conserved_IF_reduced_TDM_matrix
            implicit none
            REAL(rk), ALLOCATABLE :: intermediate_matrix(:, :)

            ALLOCATE(intermediate_matrix(Nlevels_Ftot_bra, Ndimen_Ftot_ket))

            CALL dgemm('T', 'N', Nlevels_Ftot_bra, Ndimen_Ftot_ket, Ndimen_Ftot_bra, 1.0_rk, &
                        eigen_all_hfs(index_Ftot_bra, index_represC2v_bra)%vect, Ndimen_Ftot_bra, &
                        primitive_IF_reduced_TDM_matrix, Ndimen_Ftot_bra, &
                        0.0_rk, intermediate_matrix, Nlevels_Ftot_bra)

            CALL dgemm('N', 'N', Nlevels_Ftot_bra, Nlevels_Ftot_ket, Ndimen_Ftot_ket, 1.0_rk, &
                        intermediate_matrix, Nlevels_Ftot_bra, &
                        eigen_all_hfs(index_Ftot_ket, index_represC2v_ket)%vect, Ndimen_Ftot_ket, &
                        0.0_rk, parity_conserved_IF_reduced_TDM_matrix, Nlevels_Ftot_bra)

            DEALLOCATE(intermediate_matrix)

            parity_conserved_IF_reduced_TDM_matrix = &
                parity_conserved_IF_reduced_TDM_matrix ** 2
        end subroutine construct_parity_conserved_IF_reduced_TDM_matrix

        subroutine calc_IF_EinsteinA_coefficient
            implicit none
            INTEGER(ik) :: bra, ket
            REAL(rk) :: EinsteinA, lower_energy_lower_bound, lower_energy_upper_bound, &
                        upper_energy_lower_bound, upper_energy_upper_bound, &
                        frequency_lower_bound, frequency_upper_bound

            lower_energy_lower_bound = intensity%erange_low(1) + min_eigen_value_hfs
            lower_energy_upper_bound = intensity%erange_low(2) + min_eigen_value_hfs
            upper_energy_lower_bound = intensity%erange_upp(1) + min_eigen_value_hfs
            upper_energy_upper_bound = intensity%erange_upp(2) + min_eigen_value_hfs
            frequency_lower_bound = max(0.0_rk, intensity%freq_window(1))
            frequency_upper_bound = intensity%freq_window(2)

            do ket = 1, Nlevels_Ftot_ket
                do bra = 1, Nlevels_Ftot_bra
                    if (eigen_all_hfs(index_Ftot_ket, index_represC2v_ket)%val(ket) &
                        < lower_energy_lower_bound) cycle
                    if (eigen_all_hfs(index_Ftot_ket, index_represC2v_ket)%val(ket) &
                        > lower_energy_upper_bound) cycle
                    if (eigen_all_hfs(index_Ftot_bra, index_represC2v_bra)%val(bra) &
                        < upper_energy_lower_bound) cycle
                    if (eigen_all_hfs(index_Ftot_bra, index_represC2v_bra)%val(bra) &
                        > upper_energy_upper_bound) cycle

                    parity_conserved_IF_transitions_matrix(bra, ket) = &
                        eigen_all_hfs(index_Ftot_bra, index_represC2v_bra)%val(bra) &
                        - eigen_all_hfs(index_Ftot_ket, index_represC2v_ket)%val(ket)

                    if (parity_conserved_IF_transitions_matrix(bra, ket) < 0.0_rk) cycle
                    if (parity_conserved_IF_transitions_matrix(bra, ket) < frequency_lower_bound) cycle
                    if (parity_conserved_IF_transitions_matrix(bra, ket) > frequency_upper_bound) cycle

                    EinsteinA = A_coef_s_1 &
                        * parity_conserved_IF_transitions_matrix(bra, ket) &
                        * parity_conserved_IF_transitions_matrix(bra, ket) &
                        * parity_conserved_IF_transitions_matrix(bra, ket) &
                        * parity_conserved_IF_reduced_TDM_matrix(bra, ket) &
                        / (2.0_rk * primitive_hfs_basis(index_Ftot_bra)%icontr(1)%Ftot + 1.0_rk)

                    WRITE(unit_hyperfine_transitions, "(I12, I12, E22.12, F22.12)") &
                        eigen_all_hfs(index_Ftot_bra, index_represC2v_bra)%quanta(bra)%iroot, &
                        eigen_all_hfs(index_Ftot_ket, index_represC2v_ket)%quanta(ket)%iroot, &
                        EinsteinA, &
                        parity_conserved_IF_transitions_matrix(bra, ket)
                end do
            end do
        end subroutine calc_IF_EinsteinA_coefficient

        subroutine construct_primitive_IF_reduced_TDM_matrix( &
            index_Ftot_bra, index_Ftot_ket)

            implicit none

            INTEGER(ik), INTENT(IN) :: index_Ftot_bra, index_Ftot_ket
            INTEGER(ik) :: Ndimen_Ftot_bra, Ndimen_Ftot_ket, bra, ket
            REAL(rk) :: Ftot_bra, Ftot_ket, &
                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                J_bra, J_ket, Omega_bra, Omega_ket, Itot_bra, Itot_ket
            INTEGER(ik) :: state_bra, state_ket, index_v_bra, index_v_ket, &
                            v_bra, v_ket, Lambda_bra, Lambda_ket, index_field
            type(fieldT), pointer :: field

            primitive_IF_reduced_TDM_matrix = 0.0_rk

            Ndimen_Ftot_bra = primitive_hfs_basis(index_Ftot_bra)%Ndimen
            Ndimen_Ftot_ket = primitive_hfs_basis(index_Ftot_ket)%Ndimen

            do ket = 1, Ndimen_Ftot_ket
                call get_quanta(index_Ftot_ket, ket, &
                                Ftot_ket, state_ket, index_v_ket, v_ket, Lambda_ket, &
                                S_ket, Sigma_ket, J_ket, Omega_ket, Itot=Itot_ket)

                loop_bra: do bra = 1, Ndimen_Ftot_bra
                    call get_quanta(index_Ftot_bra, bra, &
                                    Ftot_bra, state_bra, index_v_bra, v_bra, Lambda_bra, &
                                    S_bra, Sigma_bra, J_bra, Omega_bra, Itot=Itot_bra)

                    ! if (abs(J_ket - J_bra) > 1.0_rk + 10.0_rk*small_) cycle loop_bra
                    ! if (J_ket + J_bra < 1.0_rk - 10.0_rk*small_) cycle loop_bra
                    
                    if ((nint(S_ket - S_bra) /= 0) .or. (nint(Sigma_bra - Sigma_ket) /= 0)) cycle loop_bra
                    if (abs(nint(Omega_ket - Omega_bra)) > 1) cycle loop_bra
                    if (abs(nint(Omega_ket - Omega_bra)) == 0 .and. Lambda_bra /= Lambda_ket) cycle loop_bra
                    if (abs(nint(Omega_ket - Omega_bra)) == 1 .and. abs(Lambda_bra - Lambda_ket) /= 1) cycle loop_bra

                    do index_field = 1, Ndipoles
                        field => dipoletm(index_field)

                        primitive_IF_reduced_TDM_matrix(bra, ket) = &
                            primitive_IF_reduced_TDM_matrix(bra, ket) &
                            + primitive_IF_reduced_TDM_matrix_element( &
                                Ftot_bra, Ftot_ket, &
                                field, index_v_bra, index_v_ket, &
                                state_bra, state_ket, Lambda_bra, Lambda_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket, Itot_bra, Itot_ket)
                    end do
                end do loop_bra
            end do
        end subroutine construct_primitive_IF_reduced_TDM_matrix

        function primitive_IF_reduced_TDM_matrix_element( &
            Ftot_bra, Ftot_ket, &
            field, index_v_bra, index_v_ket, &
            state_bra, state_ket, Lambda_bra, Lambda_ket, &
            J_bra, J_ket, Omega_bra, Omega_ket, Itot_bra, Itot_ket) &
            result(tdm)

            implicit none

            REAL(rk), INTENT(IN) :: Ftot_bra, Ftot_ket, Itot_bra, Itot_ket, &
                                    J_bra, J_ket, Omega_bra, Omega_ket
            INTEGER(ik), INTENT(IN) :: index_v_bra, index_v_ket, &
                                        state_bra, state_ket, Lambda_bra, Lambda_ket
            type(fieldT), pointer, intent(in) :: field

            INTEGER(ik) :: ipermute, istateI_, istateF_, ilambdaI_, Lambda_ket_, isigmav, itau
            REAL(rk) :: tdm, spinI_, spinF_, f_t

            tdm = 0.0_rk

            if (nint(Itot_bra - Itot_ket) /= 0) return

            do ipermute = 0, 1
                if (ipermute == 0) then
                    istateI_ = field%istate
                    ilambdaI_ = field%lambda
                    spinI_ = field%spini

                    istateF_ = field%jstate
                    Lambda_ket_ = field%lambdaj
                    spinF_ = field%spinj
                else
                    istateF_ = field%istate
                    Lambda_ket_ = field%lambda
                    spinF_ = field%spini

                    istateI_ = field%jstate
                    ilambdaI_ = field%lambdaj
                    spinI_ = field%spinj
                endif

                if (ipermute == 1 .and. istateI_ == istateF_ &
                    .and. ilambdaI_ == Lambda_ket_ &
                    .and. nint(spinI_ - spinF_) == 0) cycle

                if (state_bra /= istateI_ .or. state_ket /= istateF_) cycle

                do isigmav = 0, 1
                    if (isigmav == 1 .and. abs(field%lambda) + abs(field%lambdaj) == 0) cycle

                    ilambdaI_ = ilambdaI_ * (-1)**isigmav
                    Lambda_ket_ = Lambda_ket_ * (-1)**isigmav

                    if (ilambdaI_ /= Lambda_bra .or. Lambda_ket_ /= Lambda_ket) cycle

                    if (abs(Lambda_bra - Lambda_ket) > 1) cycle

                    f_t = field%matelem(index_v_bra, index_v_ket)

                    if (isigmav == 1) then
                        itau = 0
                        if (ilambdaI_ == 0 .and. poten(state_bra)%parity%pm == -1) itau = itau + 1
                        if (Lambda_ket_ == 0 .and. poten(state_ket)%parity%pm == -1) itau = itau + 1
                        f_t = f_t * (-1.0_rk)**itau
                    endif

                    tdm = get_sign(J_bra + Itot_bra + Ftot_ket + 1.0_rk) &
                        * sqrt((2.0_rk * Ftot_bra + 1.0_rk) &
                            * (2.0_rk * Ftot_ket + 1.0_rk) &
                            * (2.0_rk * J_bra + 1.0_rk) &
                            * (2.0_rk * J_ket + 1.0_rk)) &
                        * Wigner6j(J_ket, Ftot_ket, Itot_bra, Ftot_bra, J_bra, 1.0_rk) &
                        * get_sign(J_bra - Omega_bra) &
                        * Wigner3j(J_bra, 1.0_rk, J_ket, &
                            -Omega_bra, Omega_bra - Omega_ket, Omega_ket) &
                        * f_t

                enddo
            enddo
        end function primitive_IF_reduced_TDM_matrix_element

        function E1_symmetry_selection(index_represC2v) result(index_partner)
            implicit none
            INTEGER(ik), INTENT(IN) :: index_represC2v
            INTEGER(ik) :: index_partner

            select case(index_represC2v)
            case(1)
                index_partner = 2
            case(2)
                index_partner = 1
            case(3)
                index_partner = 4
            case(4)
                index_partner = 3
            case default
                stop 'ERROR: homonuclear E1 intensity requires C2v symmetry.'
            end select
        end function E1_symmetry_selection

    end subroutine IF_hyperfine_E1_intensity
    
    subroutine IF_hyperfine_E2_intensity
        implicit none
        INTEGER(ik) :: index_Ftot_bra, index_Ftot_ket, &
                        index_represC2v_bra, index_represC2v_ket
        INTEGER(ik) :: Ndimen_Ftot_bra, Ndimen_Ftot_ket, &
                        Nlevels_Ftot_bra, Nlevels_Ftot_ket
        REAL(rk) :: unitConv, vacPerm
        REAL(rk), ALLOCATABLE :: primitive_IF_reduced_TQM_matrix(:,:), &
                                 parity_conserved_IF_reduced_TQM_matrix(:,:), &
                                 parity_conserved_IF_transitions_matrix(:,:)

        ! Common prefactor for the Einstein A coefficient of E2 transitions.
        ! Includes the conversion Q [a.u.] -> Q [SI],
        ! h [erg s] -> h [J s], and nu [cm^-1] -> nu [m^-1].
        unitConv = 2.012914458d-62
        ! vacuum permittivity (NIST 2018)
        vacPerm = 8.8541878128d-12 
        A_coef_s_2 = unitConv*(8.0_rk * pi**5)/(5.0_rk * vacPerm * planck)

        do index_Ftot_bra = 1, num_Ftot
            write(out, '(/A4, A, F6.1)') '', 'Start: F_bra =', Ftot_list(index_Ftot_bra)

            Ndimen_Ftot_bra = primitive_hfs_basis(index_Ftot_bra)%Ndimen

            do index_Ftot_ket = max(1, index_Ftot_bra - 2), min(num_Ftot, index_Ftot_bra + 2)
                if (Ftot_list(index_Ftot_bra) + Ftot_list(index_Ftot_ket) < 2.0_rk - 10*small_) cycle ! E2 selection rules, F_bra+F_ket >= 2

                write(out, '(/A6, A, F6.1)') '', 'Start: F_ket =', Ftot_list(index_Ftot_ket)
                Ndimen_Ftot_ket = primitive_hfs_basis(index_Ftot_ket)%Ndimen

                write(out, '(/A8, A)') '',&
                    'Construct primitive reduced transition quadrupole moment matrix'
                ALLOCATE(primitive_IF_reduced_TQM_matrix(Ndimen_Ftot_bra, Ndimen_Ftot_ket))
                CALL construct_primitive_IF_reduced_TQM_matrix( &
                        index_Ftot_bra, index_Ftot_ket)
                write(out, '(A8, A/)') '', '... done'

                do index_represC2v_bra = 1, sym%Nrepresen

                    ! symmetry selection rule
                    index_represC2v_ket = index_represC2v_bra

                    Nlevels_Ftot_bra = eigen_all_hfs(index_Ftot_bra, index_represC2v_bra)%Nlevels
                    Nlevels_Ftot_ket = eigen_all_hfs(index_Ftot_ket, index_represC2v_ket)%Nlevels

                    write(out, '(/A10, A)') '',&
                        'Construct parity conserved reduced transition quadrupole moment matrix'
                    
                    if (Nlevels_Ftot_bra <= 0 .or. Nlevels_Ftot_ket <= 0) then
                        write(out, '(A10, A/)') '', '... skipped: empty F-symmetry block'
                        cycle
                    end if

                    ALLOCATE(parity_conserved_IF_reduced_TQM_matrix(Nlevels_Ftot_bra, Nlevels_Ftot_ket))
                    CALL construct_parity_conserved_IF_reduced_TQM_matrix
                    write(out, '(A10, A/)') '', '... done'

                    write(out, '(/A10, A)') '', 'Calculate and print transitions'
                    ALLOCATE(parity_conserved_IF_transitions_matrix(Nlevels_Ftot_bra, Nlevels_Ftot_ket))
                    CALL calc_IF_EinsteinA_coefficient
                    write(out, '(A10, A/)') '', '... done'
                    
                    DEALLOCATE(parity_conserved_IF_transitions_matrix)
                    DEALLOCATE(parity_conserved_IF_reduced_TQM_matrix)
                end do           
                
                DEALLOCATE(primitive_IF_reduced_TQM_matrix)

                write(out, '(/A6, A, F6.1)') '', 'End: F_ket =', Ftot_list(index_Ftot_ket)
            end do 
            
            write(out, '(/A4, A, F6.1)') '', 'End: F_bra =', Ftot_list(index_Ftot_bra)
        end do

        close(unit=unit_hyperfine_transitions)
        write(out, '(/A)') "End: hyperfine transitions calculation"
    contains
        subroutine construct_parity_conserved_IF_reduced_TQM_matrix
            implicit none

            REAL(rk), ALLOCATABLE :: intermediate_matrix(:, :)
            ! See Eq.(58) of DOI: 10.1021/acs.jctc.1c01244
            ! ^{tau,F}M^{tau', F'} = {Psi^{tau, F}}^\dagger * ^{F}M^{F'} * Psi^{tau', F'} 

            ALLOCATE(intermediate_matrix(Nlevels_Ftot_bra, Ndimen_Ftot_ket))

            ! M_intermediate = {Psi^{tau, F}}^\dagger * ^{F}M^{F'}
            CALL dgemm('T', 'N', Nlevels_Ftot_bra, Ndimen_Ftot_ket, Ndimen_Ftot_bra, 1.0_rk, &
                        eigen_all_hfs(index_Ftot_bra, index_represC2v_bra)%vect, Ndimen_Ftot_bra, &
                        primitive_IF_reduced_TQM_matrix, Ndimen_Ftot_bra, &
                        0.0_rk, intermediate_matrix, Nlevels_Ftot_bra)

            ! ^{tau, F}M^{tau', F'} = M_intermediate * Psi^{tau', F'}           
            CALL dgemm('N', 'N', Nlevels_Ftot_bra, Nlevels_Ftot_ket, Ndimen_Ftot_ket, 1.0_rk, &
                        intermediate_matrix, Nlevels_Ftot_bra, &
                        eigen_all_hfs(index_Ftot_ket, index_represC2v_ket)%vect, Ndimen_Ftot_ket, &
                        0.0_rk, parity_conserved_IF_reduced_TQM_matrix, Nlevels_Ftot_bra)

            DEALLOCATE(intermediate_matrix)

            ! See Eq.(57) of DOI: 10.1021/acs.jctc.1c01244
            parity_conserved_IF_reduced_TQM_matrix = &
                parity_conserved_IF_reduced_TQM_matrix ** 2
        end subroutine construct_parity_conserved_IF_reduced_TQM_matrix
        
        subroutine calc_IF_EinsteinA_coefficient
            implicit none
            INTEGER(ik) :: bra, ket  
            REAL(rk)  :: EinsteinA, lower_energy_lower_bound, lower_energy_upper_bound, &
                         upper_energy_lower_bound, upper_energy_upper_bound, &
                         frequency_lower_bound, frequency_upper_bound

            lower_energy_lower_bound = intensity%erange_low(1) + min_eigen_value_hfs
            lower_energy_upper_bound = intensity%erange_low(2) + min_eigen_value_hfs
            upper_energy_lower_bound = intensity%erange_upp(1) + min_eigen_value_hfs
            upper_energy_upper_bound = intensity%erange_upp(2) + min_eigen_value_hfs
            frequency_lower_bound = max(0.0_rk, intensity%freq_window(1))
            frequency_upper_bound = intensity%freq_window(2) 
            
            do ket = 1, Nlevels_Ftot_ket
                do bra = 1, Nlevels_Ftot_bra

                    if (eigen_all_hfs(index_Ftot_ket, index_represC2v_ket)%val(ket) &
                        < lower_energy_lower_bound) cycle
                    if (eigen_all_hfs(index_Ftot_ket, index_represC2v_ket)%val(ket) &
                        > lower_energy_upper_bound) cycle
                    if (eigen_all_hfs(index_Ftot_bra, index_represC2v_bra)%val(bra) &
                        < upper_energy_lower_bound) cycle
                    if (eigen_all_hfs(index_Ftot_bra, index_represC2v_bra)%val(bra) &
                        > upper_energy_upper_bound) cycle

                    ! The bra state is the upper state
                    parity_conserved_IF_transitions_matrix(bra, ket) = &
                        eigen_all_hfs(index_Ftot_bra, index_represC2v_bra)%val(bra) &
                            - eigen_all_hfs(index_Ftot_ket, index_represC2v_ket)%val(ket)
                    
                    if (parity_conserved_IF_transitions_matrix(bra, ket) < 0.0_rk) cycle
                    if (parity_conserved_IF_transitions_matrix(bra, ket) < frequency_lower_bound) cycle
                    if (parity_conserved_IF_transitions_matrix(bra, ket) > frequency_upper_bound) cycle

                    ! Einstein A coefficients for hyperfine quadrupole: 
                    ! see Eq.(34) of M. Šimečková 
                    ! doi:10.1016/j.jqsrt.2005.07.003
                    ! degeneracy of the upper state: g_bra = 2 * F_bra + 1
                    EinsteinA = A_coef_s_2 &
                        * parity_conserved_IF_transitions_matrix(bra, ket) &
                        * parity_conserved_IF_transitions_matrix(bra, ket) &
                        * parity_conserved_IF_transitions_matrix(bra, ket) &
                        * parity_conserved_IF_transitions_matrix(bra, ket) &
                        * parity_conserved_IF_transitions_matrix(bra, ket) &
                        * parity_conserved_IF_reduced_TQM_matrix(bra, ket) & 
                        / (2.0_rk * primitive_hfs_basis(index_Ftot_bra)%icontr(1)%Ftot + 1.0_rk)

                    WRITE(unit_hyperfine_transitions, "(I12, I12, E22.12, F22.12)") &
                        eigen_all_hfs(index_Ftot_bra, index_represC2v_bra)%quanta(bra)%iroot, &
                        eigen_all_hfs(index_Ftot_ket, index_represC2v_ket)%quanta(ket)%iroot, &
                        EinsteinA, &
                        parity_conserved_IF_transitions_matrix(bra, ket)
                end do    
            end do 
        end subroutine calc_IF_EinsteinA_coefficient

        subroutine construct_primitive_IF_reduced_TQM_matrix( &
            index_Ftot_bra, index_Ftot_ket)

            implicit none
            
            INTEGER(ik), INTENT(IN) :: index_Ftot_bra, index_Ftot_ket
            INTEGER(ik) :: Ndimen_Ftot_bra, Ndimen_Ftot_ket, bra, ket
            REAL(rk) :: Ftot_bra, Ftot_ket, &
                S_bra, S_ket, Sigma_bra, Sigma_ket, &
                J_bra, J_ket, Omega_bra, Omega_ket, Itot_ket, Itot_bra
            INTEGER(ik) :: state_bra, state_ket, index_v_bra, index_v_ket, &
                            v_bra, v_ket, Lambda_bra, Lambda_ket, index_field
            type(fieldT), pointer :: field

            primitive_IF_reduced_TQM_matrix = 0.0_rk

            Ndimen_Ftot_bra = primitive_hfs_basis(index_Ftot_bra)%Ndimen
            Ndimen_Ftot_ket = primitive_hfs_basis(index_Ftot_ket)%Ndimen

            do ket = 1, Ndimen_Ftot_ket
                call get_quanta(index_Ftot_ket, ket, &
                                Ftot_ket, state_ket, index_v_ket, v_ket, Lambda_ket, &
                                S_ket, Sigma_ket, J_ket, Omega_ket, Itot=Itot_ket)

                loop_bra: do bra = 1, Ndimen_Ftot_bra
                    call get_quanta(index_Ftot_bra, bra, &
                                    Ftot_bra, state_bra, index_v_bra, v_bra, Lambda_bra, &
                                    S_bra, Sigma_bra, J_bra, Omega_bra, Itot=Itot_bra)

                    ! Selection rules due to the Wigner-3j symbol: 
                    ! Wigner3j(J, 2, J', -Omega, q, Omega')
                    ! triangle rule for rank-2:
                    ! if (abs(J_ket - J_bra) > 2) cycle loop_bra
                    ! if (J_ket + J_bra < 2 - 10.0_rk*small_) cycle loop_bra   ! excludes J=0 <-> 0,1 etc.
                    
                    if ((nint(S_ket - S_bra) /= 0 ) .or. (nint(Sigma_bra - Sigma_ket) /= 0)) cycle loop_bra
                    ! where q = Omega - Omega' and abs(q) <= 2: 
                    if (abs(nint(Omega_ket - Omega_bra)) > 2) cycle loop_bra
                    ! Lambda-selection (assume ΔΣ=0 so ΔΩ = ΔΛ): 
                    if (abs(nint(Omega_ket - Omega_bra)) == 0 .and. Lambda_bra /= Lambda_ket) cycle loop_bra
                    if (abs(nint(Omega_ket - Omega_bra)) == 1 .and. abs(Lambda_bra - Lambda_ket) /= 1) cycle loop_bra
                    if (abs(nint(Omega_ket - Omega_bra)) == 2 .and. abs(Lambda_bra - Lambda_ket) /= 2) cycle loop_bra

                    do index_field = 1, nQuadrupoles
                        field => quadrupoletm(index_field)
        
                        primitive_IF_reduced_TQM_matrix(bra, ket) = &
                            primitive_IF_reduced_TQM_matrix(bra, ket)  &
                            + primitive_IF_reduced_TQM_matrix_element( &
                                Ftot_bra, Ftot_ket, &
                                field, index_v_bra, index_v_ket, &
                                state_bra, state_ket, Lambda_bra, Lambda_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket, Itot_bra, Itot_ket)
                    end do
                end do loop_bra
            end do
        end subroutine construct_primitive_IF_reduced_TQM_matrix

        function primitive_IF_reduced_TQM_matrix_element( &
            Ftot_bra, Ftot_ket, &
            field, index_v_bra, index_v_ket, &
            state_bra, state_ket, Lambda_bra, Lambda_ket, &
            J_bra, J_ket, Omega_bra, Omega_ket, Itot_bra, Itot_ket) & 
            result(tqm)

            implicit none

            REAL(rk), INTENT(IN) :: Ftot_bra, Ftot_ket, Itot_bra, Itot_ket, &
                                    J_bra, J_ket, Omega_bra, Omega_ket
            INTEGER(ik), INTENT(IN) :: index_v_bra, index_v_ket, &
                                        state_bra, state_ket, Lambda_bra, Lambda_ket 
            type(fieldT), pointer, intent(in) :: field
            INTEGER(ik) :: ipermute, istateI_, istateF_, ilambdaI_, Lambda_ket_, isigmav, itau
            REAL(rk) :: tqm, spinI_, spinF_, vibME

            ! Adapted from do_1st_half_linestrength() in dipole.f90 and quadrupole.f90.
            ! The half-line-strength accumulation is replaced by the hyperfine
            ! reduced quadrupole matrix element tqm.

            tqm = 0.0_rk
            if (nint(Itot_bra - Itot_ket) /= 0) return
            
            do ipermute = 0, 1
                if (ipermute==0) then 
                    istateI_ = field%istate; ilambdaI_ = field%lambda; spinI_ = field%spini
                    istateF_ = field%jstate; Lambda_ket_ = field%lambdaj; spinF_ = field%spinj
                else  ! permute
                    istateF_ = field%istate; Lambda_ket_ = field%lambda; spinF_ = field%spini
                    istateI_ = field%jstate; ilambdaI_ = field%lambdaj; spinI_ = field%spinj
                endif

                ! The permutation is needed only for non-diagonal
                ! <State,Lambda,Spin|F|State',Lambda',Spin'> terms;
                ! otherwise it would lead to double counting.
                if (ipermute==1.and.istateI_==istateF_ &
                    .and.ilambdaI_==Lambda_ket_ &
                    .and.nint(spinI_-spinF_) == 0) cycle

                ! check if we at the right electronic states
                if( state_bra/=istateI_.or.state_ket/=istateF_ ) cycle
                
                ! We should also take into account that Lambda can change sign (only Lambda>0 is given in input)
                ! In order to recover other combinations we apply the symmetry transformation
                ! laboratory fixed inversion which is equivalent to the sigmav operation 
                !                    (sigmav= 0 correspond to the unitary transformation)
                do isigmav = 0, 1
                
                    ! the permutation is only needed if at least some of the quanta is not zero. 
                    ! otherwise it should be skipped to avoid the double counting.
                    if( isigmav==1.and. abs( field%lambda ) + abs( field%lambdaj )==0 ) cycle
            
                    ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                    ilambdaI_ = ilambdaI_*(-1)**isigmav
                    Lambda_ket_ = Lambda_ket_*(-1)**isigmav
                    
                    ! proceed only if the quantum numbers of the field equal to 
                    ! the corresponding <i| and |j> quantum numbers:
                    if (ilambdaI_/=Lambda_bra.or.Lambda_ket_/=Lambda_ket) cycle
                    
                    ! E2: reinforce lambda selection rule (|ΔΛ| <= 2)
                    if (abs(Lambda_bra-Lambda_ket)>2) cycle

                    vibME = field%matelem(index_v_bra, index_v_ket)
                    
                    ! the result of the symmetry transformation:
                    if (isigmav==1) then
                        itau = 0
                        if (ilambdaI_==0.and.poten(state_bra)%parity%pm==-1) itau = itau+1
                        if (Lambda_ket_==0.and.poten(state_ket)%parity%pm==-1) itau = itau+1
                        vibME = vibME*(-1.0_rk)**(itau)
                    endif

                    tqm = get_sign(J_bra + Itot_bra + Ftot_ket + 2.0_rk) &
                        * sqrt((2.0_rk * Ftot_bra + 1.0_rk) &
                                * (2.0_rk * Ftot_ket + 1.0_rk) &
                                * (2.0_rk * J_bra + 1.0_rk) &
                                * (2.0_rk * J_ket + 1.0_rk)) &
                        * Wigner6j(J_ket, Ftot_ket, Itot_bra, Ftot_bra, J_bra, 2.0_rk)&
                        * get_sign(J_bra-Omega_bra)&
                        * Wigner3j(J_bra, 2.0_rk, J_ket, &
                            -Omega_bra, Omega_bra - Omega_ket, Omega_ket)&
                        * vibME
                enddo
            enddo
        end function primitive_IF_reduced_TQM_matrix_element

    end subroutine IF_hyperfine_E2_intensity

    subroutine F1F_hyperfine_E1_intensity
        ! E1 intensities for heteronuclear two-spin molecules.
        ! This follows the single-spin F1_E1 workflow, but the primitive basis is
        ! coupled as J+I1=F1 and then F1+I2=Ftot.  Therefore the outer blocks are
        ! labelled by Ftot, while F1 remains an internal primitive-basis quantum
        ! number and appears in the recoupling factors below.
        implicit none
        INTEGER(ik) :: index_Ftot_bra, index_Ftot_ket, &
                       index_represCs_bra, index_represCs_ket
        INTEGER(ik) :: Ndimen_Ftot_bra, Ndimen_Ftot_ket, &
                       Nlevels_Ftot_bra, Nlevels_Ftot_ket
        REAL(rk), ALLOCATABLE :: primitive_F1F_reduced_TDM_matrix(:,:), &
                                 parity_conserved_F1F_reduced_TDM_matrix(:,:), &
                                 parity_conserved_F1F_transitions_matrix(:,:)

        A_coef_s_1 = 64.0E-36 * pi**4 / (3.0_rk * planck)

        do index_Ftot_bra = 1, num_Ftot
            write(out, '(/A4, A, F6.1)') '', 'Start: F_bra =', Ftot_list(index_Ftot_bra)

            Ndimen_Ftot_bra = primitive_hfs_basis(index_Ftot_bra)%Ndimen

            ! E1 tensor rank is 1, so total F follows Delta F = 0,+/-1,
            ! excluding the 0 -> 0 case.
            do index_Ftot_ket = max(1, index_Ftot_bra - 1), min(num_Ftot, index_Ftot_bra + 1)
                if ((index_Ftot_bra == index_Ftot_ket) .and. Ftot_list(index_Ftot_bra) == 0) cycle

                write(out, '(/A6, A, F6.1)') '', 'Start: F_ket =', Ftot_list(index_Ftot_ket)
                Ndimen_Ftot_ket = primitive_hfs_basis(index_Ftot_ket)%Ndimen

                write(out, '(/A8, A)') '', &
                    'Construct primitive reduced transition dipole moment matrix'
                ALLOCATE(primitive_F1F_reduced_TDM_matrix(Ndimen_Ftot_bra, Ndimen_Ftot_ket))
                CALL construct_primitive_F1F_reduced_TDM_matrix(index_Ftot_bra, index_Ftot_ket)
                write(out, '(A8, A/)') '', '... done'

                do index_represCs_bra = 1, sym%NrepresCs
                    ! Same Cs parity selection as F1_E1: electric dipole changes
                    ! total parity, so the ket block is the opposite Cs irrep.
                    index_represCs_ket = mod(index_represCs_bra, sym%NrepresCs) + 1

                    Nlevels_Ftot_bra = eigen_all_hfs(index_Ftot_bra, index_represCs_bra)%Nlevels
                    Nlevels_Ftot_ket = eigen_all_hfs(index_Ftot_ket, index_represCs_ket)%Nlevels

                    write(out, '(/A10, A)') '', &
                        'Construct parity conserved reduced transition dipole moment matrix'
                    if (Nlevels_Ftot_bra <= 0 .or. Nlevels_Ftot_ket <= 0) then
                        write(out, '(A10, A/)') '', '... skipped: empty F-parity block'
                        cycle
                    end if

                    ALLOCATE(parity_conserved_F1F_reduced_TDM_matrix( &
                        Nlevels_Ftot_bra, Nlevels_Ftot_ket))
                    CALL construct_parity_conserved_F1F_reduced_TDM_matrix
                    write(out, '(A10, A/)') '', '... done'

                    write(out, '(/A10, A)') '', 'Calculate and print transitions'
                    ALLOCATE(parity_conserved_F1F_transitions_matrix( &
                        Nlevels_Ftot_bra, Nlevels_Ftot_ket))
                    CALL calc_F1F_EinsteinA_coefficient
                    write(out, '(A10, A/)') '', '... done'

                    DEALLOCATE(parity_conserved_F1F_transitions_matrix)
                    DEALLOCATE(parity_conserved_F1F_reduced_TDM_matrix)
                end do

                DEALLOCATE(primitive_F1F_reduced_TDM_matrix)

                write(out, '(/A6, A, F6.1)') '', 'End: F_ket =', Ftot_list(index_Ftot_ket)
            end do

            write(out, '(/A4, A, F6.1)') '', 'End: F_bra =', Ftot_list(index_Ftot_bra)
        end do

        close(unit=unit_hyperfine_transitions)
        write(out, '(/A)') "End: hyperfine transitions calculation"

    contains
        subroutine construct_parity_conserved_F1F_reduced_TDM_matrix
            implicit none
            REAL(rk), ALLOCATABLE :: intermediate_matrix(:, :)

            ALLOCATE(intermediate_matrix(Nlevels_Ftot_bra, Ndimen_Ftot_ket))

            ! Same coherent transformation as in F1_E1:
            ! M_eigen = Psi_bra^T * M_primitive * Psi_ket.
            ! The difference is that the F1F eigenvectors can mix primitive
            ! components with different k, J, and F1 inside a fixed Ftot block.
            CALL dgemm('T', 'N', Nlevels_Ftot_bra, Ndimen_Ftot_ket, Ndimen_Ftot_bra, 1.0_rk, &
                        eigen_all_hfs(index_Ftot_bra, index_represCs_bra)%vect, Ndimen_Ftot_bra, &
                        primitive_F1F_reduced_TDM_matrix, Ndimen_Ftot_bra, &
                        0.0_rk, intermediate_matrix, Nlevels_Ftot_bra)

            CALL dgemm('N', 'N', Nlevels_Ftot_bra, Nlevels_Ftot_ket, Ndimen_Ftot_ket, 1.0_rk, &
                        intermediate_matrix, Nlevels_Ftot_bra, &
                        eigen_all_hfs(index_Ftot_ket, index_represCs_ket)%vect, Ndimen_Ftot_ket, &
                        0.0_rk, parity_conserved_F1F_reduced_TDM_matrix, Nlevels_Ftot_bra)

            DEALLOCATE(intermediate_matrix)

            ! Square only after the eigenvector transformation, so all primitive
            ! amplitudes have been summed coherently.
            parity_conserved_F1F_reduced_TDM_matrix = &
                parity_conserved_F1F_reduced_TDM_matrix ** 2
        end subroutine construct_parity_conserved_F1F_reduced_TDM_matrix

        subroutine calc_F1F_EinsteinA_coefficient
            implicit none
            INTEGER(ik) :: bra, ket
            REAL(rk) :: EinsteinA, lower_energy_lower_bound, lower_energy_upper_bound, &
                        upper_energy_lower_bound, upper_energy_upper_bound, &
                        frequency_lower_bound, frequency_upper_bound

            lower_energy_lower_bound = intensity%erange_low(1) + min_eigen_value_hfs
            lower_energy_upper_bound = intensity%erange_low(2) + min_eigen_value_hfs
            upper_energy_lower_bound = intensity%erange_upp(1) + min_eigen_value_hfs
            upper_energy_upper_bound = intensity%erange_upp(2) + min_eigen_value_hfs
            frequency_lower_bound = max(0.0_rk, intensity%freq_window(1))
            frequency_upper_bound = intensity%freq_window(2)

            do ket = 1, Nlevels_Ftot_ket
                do bra = 1, Nlevels_Ftot_bra
                    if (eigen_all_hfs(index_Ftot_ket, index_represCs_ket)%val(ket) &
                        < lower_energy_lower_bound) cycle
                    if (eigen_all_hfs(index_Ftot_ket, index_represCs_ket)%val(ket) &
                        > lower_energy_upper_bound) cycle
                    if (eigen_all_hfs(index_Ftot_bra, index_represCs_bra)%val(bra) &
                        < upper_energy_lower_bound) cycle
                    if (eigen_all_hfs(index_Ftot_bra, index_represCs_bra)%val(bra) &
                        > upper_energy_upper_bound) cycle

                    parity_conserved_F1F_transitions_matrix(bra, ket) = &
                        eigen_all_hfs(index_Ftot_bra, index_represCs_bra)%val(bra) &
                        - eigen_all_hfs(index_Ftot_ket, index_represCs_ket)%val(ket)

                    if (parity_conserved_F1F_transitions_matrix(bra, ket) < 0.0_rk) cycle
                    if (parity_conserved_F1F_transitions_matrix(bra, ket) &
                        < frequency_lower_bound) cycle
                    if (parity_conserved_F1F_transitions_matrix(bra, ket) &
                        > frequency_upper_bound) cycle

                    ! Same Einstein-A expression as F1_E1, but the upper-state
                    ! degeneracy is now the total-F degeneracy 2*Ftot_bra+1.
                    EinsteinA = A_coef_s_1 &
                        * parity_conserved_F1F_transitions_matrix(bra, ket) &
                        * parity_conserved_F1F_transitions_matrix(bra, ket) &
                        * parity_conserved_F1F_transitions_matrix(bra, ket) &
                        * parity_conserved_F1F_reduced_TDM_matrix(bra, ket) &
                        / (2.0_rk * Ftot_list(index_Ftot_bra) + 1.0_rk)

                    WRITE(unit_hyperfine_transitions, "(I12, I12, E22.12, F22.12)") &
                        eigen_all_hfs(index_Ftot_bra, index_represCs_bra)%quanta(bra)%iroot, &
                        eigen_all_hfs(index_Ftot_ket, index_represCs_ket)%quanta(ket)%iroot, &
                        EinsteinA, parity_conserved_F1F_transitions_matrix(bra, ket)
                end do
            end do
        end subroutine calc_F1F_EinsteinA_coefficient

        subroutine construct_primitive_F1F_reduced_TDM_matrix( &
            index_Ftot_bra, index_Ftot_ket)

            implicit none
            INTEGER(ik), INTENT(IN) :: index_Ftot_bra, index_Ftot_ket
            INTEGER(ik) :: Ndimen_Ftot_bra, Ndimen_Ftot_ket, bra, ket
            REAL(rk) :: Ftot_bra, Ftot_ket, F1_bra, F1_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket
            INTEGER(ik) :: state_bra, state_ket, index_v_bra, index_v_ket, &
                           v_bra, v_ket, Lambda_bra, Lambda_ket, index_field
            type(fieldT), pointer :: field

            primitive_F1F_reduced_TDM_matrix = 0.0_rk

            Ndimen_Ftot_bra = primitive_hfs_basis(index_Ftot_bra)%Ndimen
            Ndimen_Ftot_ket = primitive_hfs_basis(index_Ftot_ket)%Ndimen

            do ket = 1, Ndimen_Ftot_ket
                call get_quanta(index_Ftot_ket, ket, &
                                Ftot_ket, state_ket, index_v_ket, v_ket, Lambda_ket, &
                                S_ket, Sigma_ket, J_ket, Omega_ket, F1=F1_ket)

                loop_bra: do bra = 1, Ndimen_Ftot_bra
                    call get_quanta(index_Ftot_bra, bra, &
                                    Ftot_bra, state_bra, index_v_bra, v_bra, Lambda_bra, &
                                    S_bra, Sigma_bra, J_bra, Omega_bra, F1=F1_bra)

                    ! Electronic-spin and body-fixed projection checks are copied
                    ! from F1_E1.  The extra F1F checks below enforce the two
                    ! rank-1 recoupling layers before calling the Wigner symbols.
                    if ((nint(S_ket - S_bra) /= 0) &
                        .or. (nint(Sigma_bra - Sigma_ket) /= 0)) cycle loop_bra
                    if (abs(nint(Omega_ket - Omega_bra)) > 1) cycle loop_bra
                    if (abs(nint(Omega_ket - Omega_bra)) == 0 &
                        .and. Lambda_bra /= Lambda_ket) cycle loop_bra
                    if (abs(nint(Omega_ket - Omega_bra)) == 1 &
                        .and. abs(Lambda_bra - Lambda_ket) /= 1) cycle loop_bra
                    ! The rotational 3j imposes Delta J = 0,+/-1, excluding
                    ! J=0 -> J=0.  These explicit checks avoid needless 3j calls.
                    if (abs(J_bra - J_ket) > 1.0_rk + small_) cycle loop_bra
                    if (J_bra + J_ket < 1.0_rk - small_) cycle loop_bra
                    ! New relative to single-spin F1_E1: F1 is now an intermediate
                    ! angular momentum, so the second 6j also imposes
                    ! Delta F1 = 0,+/-1, excluding F1=0 -> F1=0.
                    if (abs(F1_bra - F1_ket) > 1.0_rk + small_) cycle loop_bra
                    if (F1_bra + F1_ket < 1.0_rk - small_) cycle loop_bra

                    do index_field = 1, Ndipoles
                        field => dipoletm(index_field)

                        primitive_F1F_reduced_TDM_matrix(bra, ket) = &
                            primitive_F1F_reduced_TDM_matrix(bra, ket) &
                            + primitive_F1F_reduced_TDM_matrix_element( &
                                Ftot_bra, Ftot_ket, F1_bra, F1_ket, &
                                field, index_v_bra, index_v_ket, &
                                state_bra, state_ket, Lambda_bra, Lambda_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket)
                    end do
                end do loop_bra
            end do
        end subroutine construct_primitive_F1F_reduced_TDM_matrix

        function primitive_F1F_reduced_TDM_matrix_element( &
            Ftot_bra, Ftot_ket, F1_bra, F1_ket, &
            field, index_v_bra, index_v_ket, &
            state_bra, state_ket, Lambda_bra, Lambda_ket, &
            J_bra, J_ket, Omega_bra, Omega_ket) result(tdm)

            implicit none
            REAL(rk), INTENT(IN) :: Ftot_bra, Ftot_ket, F1_bra, F1_ket, &
                                    J_bra, J_ket, Omega_bra, Omega_ket
            INTEGER(ik), INTENT(IN) :: index_v_bra, index_v_ket, &
                                       state_bra, state_ket, Lambda_bra, Lambda_ket
            type(fieldT), pointer, intent(in) :: field
            INTEGER(ik) :: ipermute, istateI_, istateF_, ilambdaI_, &
                           Lambda_ket_, isigmav, itau
            REAL(rk) :: tdm, spinI_, spinF_, f_t

            tdm = 0.0_rk
            ! The dipole field bookkeeping below is copied from F1_E1/dipole.f90.
            do ipermute = 0, 1
                if (ipermute == 0) then
                    istateI_ = field%istate
                    ilambdaI_ = field%lambda
                    spinI_ = field%spini
                    istateF_ = field%jstate
                    Lambda_ket_ = field%lambdaj
                    spinF_ = field%spinj
                else
                    istateF_ = field%istate
                    Lambda_ket_ = field%lambda
                    spinF_ = field%spini
                    istateI_ = field%jstate
                    ilambdaI_ = field%lambdaj
                    spinI_ = field%spinj
                endif

                if (ipermute == 1 .and. istateI_ == istateF_ &
                    .and. ilambdaI_ == Lambda_ket_ &
                    .and. nint(spinI_ - spinF_) == 0) cycle

                if (state_bra /= istateI_ .or. state_ket /= istateF_) cycle

                do isigmav = 0, 1
                    if (isigmav == 1 &
                        .and. abs(field%lambda) + abs(field%lambdaj) == 0) cycle

                    ilambdaI_ = ilambdaI_ * (-1)**isigmav
                    Lambda_ket_ = Lambda_ket_ * (-1)**isigmav

                    if (ilambdaI_ /= Lambda_bra .or. Lambda_ket_ /= Lambda_ket) cycle
                    if (abs(Lambda_bra - Lambda_ket) > 1) cycle

                    f_t = field%matelem(index_v_bra, index_v_ket)

                    if (isigmav == 1) then
                        itau = 0
                        if (ilambdaI_ == 0 .and. poten(state_bra)%parity%pm == -1) itau = itau + 1
                        if (Lambda_ket_ == 0 .and. poten(state_ket)%parity%pm == -1) itau = itau + 1
                        f_t = f_t * (-1.0_rk)**itau
                    endif

                    ! Different from F1_E1: two successive recouplings are needed.
                    ! First the rank-1 dipole is recoupled from J to F1 through I1;
                    ! then from F1 to total Ftot through I2.
                    tdm = get_sign(Ftot_ket + F1_bra + I2 + 1.0_rk) &
                        * sqrt((2.0_rk * Ftot_bra + 1.0_rk) &
                            * (2.0_rk * Ftot_ket + 1.0_rk)) &
                        * Wigner6j(F1_ket, Ftot_ket, I2, &
                            Ftot_bra, F1_bra, 1.0_rk) &
                        * get_sign(F1_ket + J_bra + I1 + 1.0_rk) &
                        * sqrt((2.0_rk * F1_bra + 1.0_rk) &
                            * (2.0_rk * F1_ket + 1.0_rk) &
                            * (2.0_rk * J_bra + 1.0_rk) &
                            * (2.0_rk * J_ket + 1.0_rk)) &
                        * Wigner6j(J_ket, F1_ket, I1, &
                            F1_bra, J_bra, 1.0_rk) &
                        * get_sign(J_bra - Omega_bra) &
                        * Wigner3j(J_bra, 1.0_rk, J_ket, &
                            -Omega_bra, Omega_bra - Omega_ket, Omega_ket) &
                        * f_t
                enddo
            enddo
        end function primitive_F1F_reduced_TDM_matrix_element

    end subroutine F1F_hyperfine_E1_intensity

    subroutine F1F_hyperfine_E2_intensity
        ! E2 intensities for heteronuclear molecules with two non-zero nuclear spins.
        ! The primitive basis is coupled as J+I1=F1 and then F1+I2=Ftot.  
        ! Therefore the outer blocks are labelled by Ftot, 
        ! while F1 remains an internal quantum number in the recoupling factors.
        implicit none
        INTEGER(ik) :: index_Ftot_bra, index_Ftot_ket, &
                       index_represCs_bra, index_represCs_ket
        INTEGER(ik) :: Ndimen_Ftot_bra, Ndimen_Ftot_ket, &
                       Nlevels_Ftot_bra, Nlevels_Ftot_ket
        REAL(rk) :: unitConv, vacPerm
        REAL(rk), ALLOCATABLE :: primitive_F1F_reduced_TQM_matrix(:,:), &
                                 parity_conserved_F1F_reduced_TQM_matrix(:,:), &
                                 parity_conserved_F1F_transitions_matrix(:,:)

        ! Common prefactor for the Einstein A coefficient of E2 transitions.
        ! Includes the conversion Q [a.u.] -> Q [SI],
        ! h [erg s] -> h [J s], and nu [cm^-1] -> nu [m^-1].
        unitConv = 2.012914458d-62
        vacPerm = 8.8541878128d-12
        A_coef_s_2 = unitConv * (8.0_rk * pi**5) / (5.0_rk * vacPerm * planck)

        do index_Ftot_bra = 1, num_Ftot
            write(out, '(/A4, A, F6.1)') '', 'Start: F_bra =', Ftot_list(index_Ftot_bra)

            Ndimen_Ftot_bra = primitive_hfs_basis(index_Ftot_bra)%Ndimen

            ! E2 tensor rank is 2, so |Delta F| <= 2 and F_bra+F_ket >= 2.
            do index_Ftot_ket = max(1, index_Ftot_bra - 2), min(num_Ftot, index_Ftot_bra + 2)
                if (Ftot_list(index_Ftot_bra) + Ftot_list(index_Ftot_ket) &
                    < 2.0_rk - 10.0_rk*small_) cycle

                write(out, '(/A6, A, F6.1)') '', 'Start: F_ket =', Ftot_list(index_Ftot_ket)
                Ndimen_Ftot_ket = primitive_hfs_basis(index_Ftot_ket)%Ndimen

                write(out, '(/A8, A)') '', &
                    'Construct primitive reduced transition quadrupole moment matrix'
                ALLOCATE(primitive_F1F_reduced_TQM_matrix(Ndimen_Ftot_bra, Ndimen_Ftot_ket))
                CALL construct_primitive_F1F_reduced_TQM_matrix(index_Ftot_bra, index_Ftot_ket)
                write(out, '(A8, A/)') '', '... done'

                do index_represCs_bra = 1, sym%NrepresCs
                    ! Electric quadrupole transitions conserve total parity.
                    index_represCs_ket = index_represCs_bra

                    Nlevels_Ftot_bra = eigen_all_hfs(index_Ftot_bra, index_represCs_bra)%Nlevels
                    Nlevels_Ftot_ket = eigen_all_hfs(index_Ftot_ket, index_represCs_ket)%Nlevels

                    write(out, '(/A10, A)') '', &
                        'Construct parity conserved reduced transition quadrupole moment matrix'
                    if (Nlevels_Ftot_bra <= 0 .or. Nlevels_Ftot_ket <= 0) then
                        write(out, '(A10, A/)') '', '... skipped: empty F-parity block'
                        cycle
                    end if

                    ALLOCATE(parity_conserved_F1F_reduced_TQM_matrix( &
                        Nlevels_Ftot_bra, Nlevels_Ftot_ket))
                    CALL construct_parity_conserved_F1F_reduced_TQM_matrix
                    write(out, '(A10, A/)') '', '... done'

                    write(out, '(/A10, A)') '', 'Calculate and print transitions'
                    ALLOCATE(parity_conserved_F1F_transitions_matrix( &
                        Nlevels_Ftot_bra, Nlevels_Ftot_ket))
                    CALL calc_F1F_EinsteinA_coefficient
                    write(out, '(A10, A/)') '', '... done'

                    DEALLOCATE(parity_conserved_F1F_transitions_matrix)
                    DEALLOCATE(parity_conserved_F1F_reduced_TQM_matrix)
                end do

                DEALLOCATE(primitive_F1F_reduced_TQM_matrix)

                write(out, '(/A6, A, F6.1)') '', 'End: F_ket =', Ftot_list(index_Ftot_ket)
            end do

            write(out, '(/A4, A, F6.1)') '', 'End: F_bra =', Ftot_list(index_Ftot_bra)
        end do

        close(unit=unit_hyperfine_transitions)
        write(out, '(/A)') "End: hyperfine transitions calculation"

    contains
        subroutine construct_parity_conserved_F1F_reduced_TQM_matrix
            implicit none
            REAL(rk), ALLOCATABLE :: intermediate_matrix(:, :)

            ALLOCATE(intermediate_matrix(Nlevels_Ftot_bra, Ndimen_Ftot_ket))

            ! Coherent transformation to the hyperfine eigenstate basis:
            ! M_eigen = Psi_bra^T * M_primitive * Psi_ket.
            ! The difference is that the F1F eigenvectors can mix primitive
            ! components with different k, J, and F1 inside a fixed Ftot block.
            CALL dgemm('T', 'N', Nlevels_Ftot_bra, Ndimen_Ftot_ket, Ndimen_Ftot_bra, 1.0_rk, &
                        eigen_all_hfs(index_Ftot_bra, index_represCs_bra)%vect, Ndimen_Ftot_bra, &
                        primitive_F1F_reduced_TQM_matrix, Ndimen_Ftot_bra, &
                        0.0_rk, intermediate_matrix, Nlevels_Ftot_bra)

            CALL dgemm('N', 'N', Nlevels_Ftot_bra, Nlevels_Ftot_ket, Ndimen_Ftot_ket, 1.0_rk, &
                        intermediate_matrix, Nlevels_Ftot_bra, &
                        eigen_all_hfs(index_Ftot_ket, index_represCs_ket)%vect, Ndimen_Ftot_ket, &
                        0.0_rk, parity_conserved_F1F_reduced_TQM_matrix, Nlevels_Ftot_bra)

            DEALLOCATE(intermediate_matrix)

            ! Square only after the eigenvector transformation, so all primitive
            ! amplitudes have been summed coherently.
            parity_conserved_F1F_reduced_TQM_matrix = &
                parity_conserved_F1F_reduced_TQM_matrix ** 2
        end subroutine construct_parity_conserved_F1F_reduced_TQM_matrix

        subroutine calc_F1F_EinsteinA_coefficient
            implicit none
            INTEGER(ik) :: bra, ket
            REAL(rk) :: EinsteinA, lower_energy_lower_bound, lower_energy_upper_bound, &
                        upper_energy_lower_bound, upper_energy_upper_bound, &
                        frequency_lower_bound, frequency_upper_bound

            lower_energy_lower_bound = intensity%erange_low(1) + min_eigen_value_hfs
            lower_energy_upper_bound = intensity%erange_low(2) + min_eigen_value_hfs
            upper_energy_lower_bound = intensity%erange_upp(1) + min_eigen_value_hfs
            upper_energy_upper_bound = intensity%erange_upp(2) + min_eigen_value_hfs
            frequency_lower_bound = max(0.0_rk, intensity%freq_window(1))
            frequency_upper_bound = intensity%freq_window(2)

            do ket = 1, Nlevels_Ftot_ket
                do bra = 1, Nlevels_Ftot_bra
                    if (eigen_all_hfs(index_Ftot_ket, index_represCs_ket)%val(ket) &
                        < lower_energy_lower_bound) cycle
                    if (eigen_all_hfs(index_Ftot_ket, index_represCs_ket)%val(ket) &
                        > lower_energy_upper_bound) cycle
                    if (eigen_all_hfs(index_Ftot_bra, index_represCs_bra)%val(bra) &
                        < upper_energy_lower_bound) cycle
                    if (eigen_all_hfs(index_Ftot_bra, index_represCs_bra)%val(bra) &
                        > upper_energy_upper_bound) cycle

                    parity_conserved_F1F_transitions_matrix(bra, ket) = &
                        eigen_all_hfs(index_Ftot_bra, index_represCs_bra)%val(bra) &
                        - eigen_all_hfs(index_Ftot_ket, index_represCs_ket)%val(ket)

                    if (parity_conserved_F1F_transitions_matrix(bra, ket) < 0.0_rk) cycle
                    if (parity_conserved_F1F_transitions_matrix(bra, ket) &
                        < frequency_lower_bound) cycle
                    if (parity_conserved_F1F_transitions_matrix(bra, ket) &
                        > frequency_upper_bound) cycle

                    ! E2 Einstein-A coefficient with upper-state degeneracy
                    ! g_bra = 2*Ftot_bra+1.
                    EinsteinA = A_coef_s_2 &
                        * parity_conserved_F1F_transitions_matrix(bra, ket) &
                        * parity_conserved_F1F_transitions_matrix(bra, ket) &
                        * parity_conserved_F1F_transitions_matrix(bra, ket) &
                        * parity_conserved_F1F_transitions_matrix(bra, ket) &
                        * parity_conserved_F1F_transitions_matrix(bra, ket) &
                        * parity_conserved_F1F_reduced_TQM_matrix(bra, ket) &
                        / (2.0_rk * Ftot_list(index_Ftot_bra) + 1.0_rk)

                    WRITE(unit_hyperfine_transitions, "(I12, I12, E22.12, F22.12)") &
                        eigen_all_hfs(index_Ftot_bra, index_represCs_bra)%quanta(bra)%iroot, &
                        eigen_all_hfs(index_Ftot_ket, index_represCs_ket)%quanta(ket)%iroot, &
                        EinsteinA, parity_conserved_F1F_transitions_matrix(bra, ket)
                end do
            end do
        end subroutine calc_F1F_EinsteinA_coefficient

        subroutine construct_primitive_F1F_reduced_TQM_matrix( &
            index_Ftot_bra, index_Ftot_ket)

            implicit none
            INTEGER(ik), INTENT(IN) :: index_Ftot_bra, index_Ftot_ket
            INTEGER(ik) :: Ndimen_Ftot_bra, Ndimen_Ftot_ket, bra, ket
            REAL(rk) :: Ftot_bra, Ftot_ket, F1_bra, F1_ket, &
                        S_bra, S_ket, Sigma_bra, Sigma_ket, &
                        J_bra, J_ket, Omega_bra, Omega_ket
            INTEGER(ik) :: state_bra, state_ket, index_v_bra, index_v_ket, &
                           v_bra, v_ket, Lambda_bra, Lambda_ket, index_field
            type(fieldT), pointer :: field

            primitive_F1F_reduced_TQM_matrix = 0.0_rk

            Ndimen_Ftot_bra = primitive_hfs_basis(index_Ftot_bra)%Ndimen
            Ndimen_Ftot_ket = primitive_hfs_basis(index_Ftot_ket)%Ndimen

            do ket = 1, Ndimen_Ftot_ket
                call get_quanta(index_Ftot_ket, ket, &
                                Ftot_ket, state_ket, index_v_ket, v_ket, Lambda_ket, &
                                S_ket, Sigma_ket, J_ket, Omega_ket, F1=F1_ket)

                loop_bra: do bra = 1, Ndimen_Ftot_bra
                    call get_quanta(index_Ftot_bra, bra, &
                                    Ftot_bra, state_bra, index_v_bra, v_bra, Lambda_bra, &
                                    S_bra, Sigma_bra, J_bra, Omega_bra, F1=F1_bra)

                    ! Electronic spin is unchanged by the electric quadrupole operator.
                    if ((nint(S_ket - S_bra) /= 0) &
                        .or. (nint(Sigma_bra - Sigma_ket) /= 0)) cycle loop_bra
                    if (abs(nint(Omega_ket - Omega_bra)) > 2) cycle loop_bra
                    if (abs(J_bra - J_ket) > 2.0_rk + small_) cycle loop_bra
                    if (J_bra + J_ket < 2.0_rk - small_) cycle loop_bra
                    if (abs(F1_bra - F1_ket) > 2.0_rk + small_) cycle loop_bra
                    if (F1_bra + F1_ket < 2.0_rk - small_) cycle loop_bra

                    do index_field = 1, Nquadrupoles
                        field => quadrupoletm(index_field)

                        primitive_F1F_reduced_TQM_matrix(bra, ket) = &
                            primitive_F1F_reduced_TQM_matrix(bra, ket) &
                            + primitive_F1F_reduced_TQM_matrix_element( &
                                Ftot_bra, Ftot_ket, F1_bra, F1_ket, &
                                field, index_v_bra, index_v_ket, &
                                state_bra, state_ket, Lambda_bra, Lambda_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket)
                    end do
                end do loop_bra
            end do
        end subroutine construct_primitive_F1F_reduced_TQM_matrix

        function primitive_F1F_reduced_TQM_matrix_element( &
            Ftot_bra, Ftot_ket, F1_bra, F1_ket, &
            field, index_v_bra, index_v_ket, &
            state_bra, state_ket, Lambda_bra, Lambda_ket, &
            J_bra, J_ket, Omega_bra, Omega_ket) result(tqm)

            implicit none
            REAL(rk), INTENT(IN) :: Ftot_bra, Ftot_ket, F1_bra, F1_ket, &
                                    J_bra, J_ket, Omega_bra, Omega_ket
            INTEGER(ik), INTENT(IN) :: index_v_bra, index_v_ket, &
                                       state_bra, state_ket, Lambda_bra, Lambda_ket
            type(fieldT), pointer, intent(in) :: field
            INTEGER(ik) :: ipermute, istateI_, istateF_, ilambdaI_, &
                           Lambda_ket_, isigmav, itau
            REAL(rk) :: tqm, spinI_, spinF_, vibME

            tqm = 0.0_rk
            ! The field bookkeeping follows quadrupole.f90.
            do ipermute = 0, 1
                if (ipermute == 0) then
                    istateI_ = field%istate
                    ilambdaI_ = field%lambda
                    spinI_ = field%spini
                    istateF_ = field%jstate
                    Lambda_ket_ = field%lambdaj
                    spinF_ = field%spinj
                else
                    istateF_ = field%istate
                    Lambda_ket_ = field%lambda
                    spinF_ = field%spini
                    istateI_ = field%jstate
                    ilambdaI_ = field%lambdaj
                    spinI_ = field%spinj
                endif

                if (ipermute == 1 .and. istateI_ == istateF_ &
                    .and. ilambdaI_ == Lambda_ket_ &
                    .and. nint(spinI_ - spinF_) == 0) cycle

                if (state_bra /= istateI_ .or. state_ket /= istateF_) cycle

                do isigmav = 0, 1
                    if (isigmav == 1 &
                        .and. abs(field%lambda) + abs(field%lambdaj) == 0) cycle

                    ilambdaI_ = ilambdaI_ * (-1)**isigmav
                    Lambda_ket_ = Lambda_ket_ * (-1)**isigmav

                    if (ilambdaI_ /= Lambda_bra .or. Lambda_ket_ /= Lambda_ket) cycle
                    if (abs(Lambda_bra - Lambda_ket) > 2) cycle

                    vibME = field%matelem(index_v_bra, index_v_ket)

                    if (isigmav == 1) then
                        itau = 0
                        if (ilambdaI_ == 0 .and. poten(state_bra)%parity%pm == -1) itau = itau + 1
                        if (Lambda_ket_ == 0 .and. poten(state_ket)%parity%pm == -1) itau = itau + 1
                        vibME = vibME * (-1.0_rk)**itau
                    endif

                    ! Recouple the rank-2 quadrupole tensor successively from
                    ! J to F1 through I1 and from F1 to Ftot through I2.
                    tqm = get_sign(Ftot_ket + F1_bra + I2 + 2.0_rk) &
                        * sqrt((2.0_rk * Ftot_bra + 1.0_rk) &
                            * (2.0_rk * Ftot_ket + 1.0_rk)) &
                        * Wigner6j(F1_ket, Ftot_ket, I2, &
                            Ftot_bra, F1_bra, 2.0_rk) &
                        * get_sign(F1_ket + J_bra + I1 + 2.0_rk) &
                        * sqrt((2.0_rk * F1_bra + 1.0_rk) &
                            * (2.0_rk * F1_ket + 1.0_rk) &
                            * (2.0_rk * J_bra + 1.0_rk) &
                            * (2.0_rk * J_ket + 1.0_rk)) &
                        * Wigner6j(J_ket, F1_ket, I1, &
                            F1_bra, J_bra, 2.0_rk) &
                        * get_sign(J_bra - Omega_bra) &
                        * Wigner3j(J_bra, 2.0_rk, J_ket, &
                            -Omega_bra, Omega_bra - Omega_ket, Omega_ket) &
                        * vibME
                enddo
            enddo
        end function primitive_F1F_reduced_TQM_matrix_element

    end subroutine F1F_hyperfine_E2_intensity

    subroutine F1_hyperfine_intensity
        implicit none
        INTEGER(ik) :: index_F1_bra, index_F1_ket, &
                       index_represCs_bra, index_represCs_ket
        INTEGER(ik) :: Ndimen_F1_bra, Ndimen_F1_ket, &
                       Nlevels_F1_bra, Nlevels_F1_ket
        REAL(rk), ALLOCATABLE :: primitive_F1_reduced_TDM_matrix(:,:), &
                                 parity_conserved_F1_reduced_TDM_matrix(:,:), &
                                 parity_conserved_F1_transitions_matrix(:,:)


        A_coef_s_1 = 64.0E-36 * pi**4  / (3.0_rk * planck)

        ! write(unit_hyperfine_transitions, "(A12, A12, A22, A22, A22)") &
        !     "N_upper", "N_lower", "Einstein-A [s-1]", "nu [cm-1]", "S [Debye^2]"

        do index_F1_bra = 1, num_F1
            write(out, '(/A4, A, F6.1)') '', 'Start: F_bra =', F1_list(index_F1_bra)

            Ndimen_F1_bra = primitive_hfs_basis(index_F1_bra)%Ndimen

            do index_F1_ket = max(1, index_F1_bra - 1), min(num_F1, index_F1_bra + 1) ! F1 selection rules
                if ((index_F1_bra == index_F1_ket) .and. F1_list(index_F1_bra) == 0) cycle

                write(out, '(/A6, A, F6.1)') '', 'Start: F_ket =', F1_list(index_F1_ket)
                Ndimen_F1_ket = primitive_hfs_basis(index_F1_ket)%Ndimen

                write(out, '(/A8, A)') '',&
                    'Construct primitive reduced transition dipole moment matrix'
                ALLOCATE(primitive_F1_reduced_TDM_matrix(Ndimen_F1_bra, Ndimen_F1_ket))
                CALL construct_primitive_F1_reduced_TDM_matrix( &
                        index_F1_bra, index_F1_ket)
                write(out, '(A8, A/)') '', '... done'

                do index_represCs_bra = 1, sym%NrepresCs
                    ! parity selection rule
                    index_represCs_ket = mod(index_represCs_bra, sym%NrepresCs) + 1 

                    Nlevels_F1_bra = eigen_all_hfs(index_F1_bra, index_represCs_bra)%Nlevels
                    Nlevels_F1_ket = eigen_all_hfs(index_F1_ket, index_represCs_ket)%Nlevels

                    write(out, '(/A10, A)') '',&
                        'Construct parity conserved reduced transition dipole moment matrix'
                    if (Nlevels_F1_bra <= 0 .or. Nlevels_F1_ket <= 0) then
                        write(out, '(A10, A/)') '', '... skipped: empty F-parity block'
                        cycle
                    end if
                    ALLOCATE(parity_conserved_F1_reduced_TDM_matrix(Nlevels_F1_bra, Nlevels_F1_ket))                    
                    CALL construct_parity_conserved_F1_reduced_TDM_matrix
                    write(out, '(A10, A/)') '', '... done'

                    write(out, '(/A10, A)') '',&
                        'Calculate and print transitions'
                    ALLOCATE(parity_conserved_F1_transitions_matrix(Nlevels_F1_bra, Nlevels_F1_ket))                    
                    CALL construct_parity_conserved_F1_transitions_matrix
                    write(out, '(A10, A/)') '', '... done'
                    
                    DEALLOCATE(parity_conserved_F1_transitions_matrix)
                    DEALLOCATE(parity_conserved_F1_reduced_TDM_matrix)
                end do           
                
                DEALLOCATE(primitive_F1_reduced_TDM_matrix)

                write(out, '(/A6, A, F6.1)') '', 'End: F_ket =', F1_list(index_F1_ket)
            end do 
            
            write(out, '(/A4, A, F6.1)') '', 'End: F_bra =', F1_list(index_F1_bra)
        end do

        close(unit=unit_hyperfine_transitions)
        write(out, '(/A)') "End: hyperfine transitions calculation"
    contains
        subroutine construct_parity_conserved_F1_reduced_TDM_matrix
            ! | <psi_m^{tau,F}| T^1(mu) |psi_n^{tau',F'}> | ** 2
            implicit none

            REAL(rk), ALLOCATABLE :: intermediate_matrix(:, :)
            ! See Eq.(58) of DOI: 10.1021/acs.jctc.1c01244
            ! ! ^{tau,F}M^{tau', F'} = {Psi^{tau, F}}^\dagger * ^{F}M^{F'} * Psi^{tau', F'} 
            ! parity_conserved_F1_reduced_TDM_matrix = &
            !     matmul(matmul(transpose(eigen_all_hfs(index_F1_bra, index_represCs_bra)%vect), &
            !                 primitive_F1_reduced_TDM_matrix), &
            !         eigen_all_hfs(index_F1_ket, index_represCs_ket)%vect)

            ALLOCATE(intermediate_matrix(Nlevels_F1_bra, Ndimen_F1_ket))

            ! M_intermediate = {Psi^{tau, F}}^\dagger * ^{F}M^{F'}
            CALL dgemm('T', 'N', Nlevels_F1_bra, Ndimen_F1_ket, Ndimen_F1_bra, 1.0_rk, &
                        eigen_all_hfs(index_F1_bra, index_represCs_bra)%vect, Ndimen_F1_bra, &
                        primitive_F1_reduced_TDM_matrix, Ndimen_F1_bra, &
                        0.0_rk, intermediate_matrix, Nlevels_F1_bra)

            ! ^{tau, F}M^{tau', F'} = M_intermediate * Psi^{tau', F'}           
            CALL dgemm('N', 'N', Nlevels_F1_bra, Nlevels_F1_ket, Ndimen_F1_ket, 1.0_rk, &
                        intermediate_matrix, Nlevels_F1_bra, &
                        eigen_all_hfs(index_F1_ket, index_represCs_ket)%vect, Ndimen_F1_ket, &
                        0.0_rk, parity_conserved_F1_reduced_TDM_matrix, Nlevels_F1_bra)

            DEALLOCATE(intermediate_matrix)

            ! See Eq.(57) of DOI: 10.1021/acs.jctc.1c01244
            parity_conserved_F1_reduced_TDM_matrix = &
                parity_conserved_F1_reduced_TDM_matrix ** 2
        end subroutine construct_parity_conserved_F1_reduced_TDM_matrix
        
        subroutine construct_parity_conserved_F1_transitions_matrix
            implicit none
            INTEGER(ik) :: bra, ket  
            REAL(rk)  :: EinsteinA, lower_energy_lower_bound, lower_energy_upper_bound, &
                         upper_energy_lower_bound, upper_energy_upper_bound, &
                         frequency_lower_bound, frequency_upper_bound

            lower_energy_lower_bound = intensity%erange_low(1) + min_eigen_value_hfs
            lower_energy_upper_bound = intensity%erange_low(2) + min_eigen_value_hfs
            upper_energy_lower_bound = intensity%erange_upp(1) + min_eigen_value_hfs
            upper_energy_upper_bound = intensity%erange_upp(2) + min_eigen_value_hfs
            frequency_lower_bound = max(0.0_rk, intensity%freq_window(1))
            frequency_upper_bound = intensity%freq_window(2)
            
            do ket = 1, Nlevels_F1_ket
                do bra = 1, Nlevels_F1_bra

                    if (eigen_all_hfs(index_F1_ket, index_represCs_ket)%val(ket) &
                        < lower_energy_lower_bound) cycle
                    if (eigen_all_hfs(index_F1_ket, index_represCs_ket)%val(ket) &
                        > lower_energy_upper_bound) cycle

                    if (eigen_all_hfs(index_F1_bra, index_represCs_bra)%val(bra) &
                        < upper_energy_lower_bound) cycle
                    if (eigen_all_hfs(index_F1_bra, index_represCs_bra)%val(bra) &
                        > upper_energy_upper_bound) cycle

                    ! The bra state is the final state
                    parity_conserved_F1_transitions_matrix(bra, ket) = &
                        eigen_all_hfs(index_F1_bra, index_represCs_bra)%val(bra) &
                            - eigen_all_hfs(index_F1_ket, index_represCs_ket)%val(ket)
                    
                    if (parity_conserved_F1_transitions_matrix(bra, ket) < 0.0_rk) cycle

                    if (parity_conserved_F1_transitions_matrix(bra, ket) &
                        < frequency_lower_bound) cycle
                    if (parity_conserved_F1_transitions_matrix(bra, ket) &
                        > frequency_upper_bound) cycle

                    ! Einstein A coefficients: see Eq.(23) of Western 
                    ! doi:10.1016/j.jqsrt.2016.04.010
                    ! degeneracy of the final state: g_bra = 2 * F_bra + 1
                    EinsteinA = A_coef_s_1 &
                        * parity_conserved_F1_transitions_matrix(bra, ket) &
                        * parity_conserved_F1_transitions_matrix(bra, ket) &
                        * parity_conserved_F1_transitions_matrix(bra, ket) &
                        * parity_conserved_F1_reduced_TDM_matrix(bra, ket) &
                        / (2.0_rk * primitive_hfs_basis(index_F1_bra)%icontr(1)%F1 + 1.0_rk)

                    WRITE(unit_hyperfine_transitions, "(I12, I12, E22.12, F22.12)") &
                        eigen_all_hfs(index_F1_bra, index_represCs_bra)%quanta(bra)%iroot, &
                        eigen_all_hfs(index_F1_ket, index_represCs_ket)%quanta(ket)%iroot, &
                        EinsteinA, &
                        parity_conserved_F1_transitions_matrix(bra, ket)
                end do    
            end do 
        end subroutine construct_parity_conserved_F1_transitions_matrix

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
        type(fieldT), pointer :: field

        primitive_F1_reduced_TDM_matrix = 0.0_rk

            Ndimen_F1_bra = primitive_hfs_basis(index_F1_bra)%Ndimen
            Ndimen_F1_ket = primitive_hfs_basis(index_F1_ket)%Ndimen

        do ket = 1, Ndimen_F1_ket
            call get_quanta(index_F1_ket, ket, &
                            F1_ket, state_ket, index_v_ket, v_ket, Lambda_ket, &
                            S_ket, Sigma_ket, J_ket, Omega_ket)

            loop_bra: do bra = 1, Ndimen_F1_bra
                call get_quanta(index_F1_bra, bra, &
                                F1_bra, state_bra, index_v_bra, v_bra, Lambda_bra, &
                                S_bra, Sigma_bra, J_bra, Omega_bra)

                if ((nint(S_ket - S_bra) /= 0 ) &
                    .or. (nint(Sigma_bra - Sigma_ket) /= 0)) cycle loop_bra

                ! Selection rules due to the Wigner-3j symbol: 
                ! Wigner3j(J, 1, J', -Omega, q, Omega')
                ! where q = Omega - Omega' and abs(q) <= 1
                if (abs(nint(Omega_ket - Omega_bra)) > 1) cycle loop_bra 
                if (abs(nint(Omega_ket - Omega_bra)) == 0 &
                    .and. Lambda_bra /= Lambda_ket) cycle loop_bra
                if (abs(nint(Omega_ket - Omega_bra)) == 1 &
                    .and. abs(Lambda_bra - Lambda_ket) /= 1) cycle loop_bra 

                do index_field = 1, Ndipoles
                    field => dipoletm(index_field)
    
                    primitive_F1_reduced_TDM_matrix(bra, ket) = &
                        primitive_F1_reduced_TDM_matrix(bra, ket)  &
                        + primitive_F1_reduced_TDM_matrix_element( &
                            F1_bra, F1_ket, &
                            field, index_v_bra, index_v_ket, &
                            state_bra, state_ket, Lambda_bra, Lambda_ket, &
                            J_bra, J_ket, Omega_bra, Omega_ket)
                end do
            end do loop_bra
        end do
    end subroutine construct_primitive_F1_reduced_TDM_matrix

    function primitive_F1_reduced_TDM_matrix_element( &
        F1_bra, F1_ket, &
        field, index_v_bra, index_v_ket, &
        state_bra, state_ket, Lambda_bra, Lambda_ket, &
        J_bra, J_ket, Omega_bra, Omega_ket) &
        result(tdm)

        implicit none

        REAL(rk), INTENT(IN) :: F1_bra, F1_ket, &
                                J_bra, J_ket, Omega_bra, Omega_ket
        INTEGER(ik), INTENT(IN) :: index_v_bra, index_v_ket, &
                                    state_bra, state_ket, Lambda_bra, Lambda_ket 
        type(fieldT), pointer, intent(in) :: field
        INTEGER(ik) :: ipermute, istateI_, istateF_, ilambdaI_, Lambda_ket_, isigmav, itau
        REAL(rk) :: tdm, spinI_, spinF_, f_t

        ! The codes are copied from
        ! the do_1st_half_linestrength() subroutine in dipole.f90
        ! except the lines for 'tdm'. 

        tdm = 0.0_rk
        do ipermute = 0, 1
            if (ipermute==0) then
                istateI_ = field%istate; ilambdaI_ = field%lambda; spinI_ = field%spini
                istateF_ = field%jstate; Lambda_ket_ = field%lambdaj; spinF_ = field%spinj
            else  ! permute
                istateF_ = field%istate; Lambda_ket_ = field%lambda; spinF_ = field%spini
                istateI_ = field%jstate; ilambdaI_ = field%lambdaj; spinI_ = field%spinj
            endif

            ! however the permutation makes sense only when for non diagonal 
            ! <State,Lambda,Spin|F|State',Lambda',Spin'>
            ! otherwise it will cause a double counting:
            if (ipermute==1.and.istateI_==istateF_ &
                .and.ilambdaI_==Lambda_ket_ &
                .and.nint(spinI_-spinF_) == 0) cycle

            ! check if we at the right electronic states
            if( state_bra/=istateI_.or.state_ket/=istateF_ ) cycle
            
            ! We should also take into account that Lambda can change sign (only Lambda>0 is given in input)
            ! In order to recover other combinations we apply the symmetry transformation
            ! laboratory fixed inversion which is equivalent to the sigmav operation 
            !                    (sigmav= 0 correspond to the unitary transformation)
            do isigmav = 0, 1
            
                ! the permutation is only needed if at least some of the quanta is not zero. 
                ! otherwise it should be skipped to avoid the double counting.
                if( isigmav==1.and. abs( field%lambda ) + abs( field%lambdaj )==0 ) cycle
        
                ! do the sigmav transformations (it simply changes the sign of lambda and sigma simultaneously)
                ilambdaI_ = ilambdaI_*(-1)**isigmav
                Lambda_ket_ = Lambda_ket_*(-1)**isigmav
                
                ! proceed only if the quantum numbers of the field equal to 
                ! the corresponding <i| and |j> quantum numbers:
                if (ilambdaI_/=Lambda_bra.or.Lambda_ket_/=Lambda_ket) cycle
                
                ! check the selection rule Delta Lambda = +/1
                if (abs(Lambda_bra-Lambda_ket)>1) cycle

                f_t = field%matelem(index_v_bra, index_v_ket)
                
                ! the result of the symmetry transformation:
                if (isigmav==1) then
                    itau = 0
                    if (ilambdaI_==0.and.poten(state_bra)%parity%pm==-1) itau = itau+1
                    if (Lambda_ket_==0.and.poten(state_ket)%parity%pm==-1) itau = itau+1
                    f_t = f_t*(-1.0_rk)**(itau)
                endif

                ! See Eqs. (59) to (61) of DOI: 10.1021/acs.jctc.1c01244
                tdm = get_sign(J_bra + I1 + F1_ket + 1.0_rk) &
                    * sqrt((2.0_rk * F1_bra + 1.0_rk) &
                            * (2.0_rk * F1_ket + 1.0_rk) &
                            * (2.0_rk * J_bra + 1.0_rk) &
                            * (2.0_rk * J_ket + 1.0_rk)) &
                    * Wigner6j(J_ket, F1_ket, I1, F1_bra, J_bra, 1.0_rk)&
                    * get_sign(J_bra-Omega_bra)&
                    * Wigner3j(J_bra, 1.0_rk, J_ket, &
                        -Omega_bra, Omega_bra - Omega_ket, Omega_ket)&
                    * f_t
            enddo
        enddo
    end function primitive_F1_reduced_TDM_matrix_element
    
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


end module F1_intensity
