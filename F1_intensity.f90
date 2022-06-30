module F1_intensity

    use accuracy
    use symmetry, only : sym
    use diatom_module, only : Ndipoles, dipoletm, intensity, job, fieldT, poten
    use F1_hyperfine, only : primitive_F1_basis, F1_list, eigen_all_F1, min_eigen_value_F1, &
                            get_sign, I1, Wigner3j, Wigner6j, num_F1, get_quanta

    implicit none

    REAL(rk), ALLOCATABLE :: primitive_F1_reduced_TDM_matrix(:,:)
    INTEGER(ik), PARAMETER :: unit_hyperfine_transitions = 66
    real(rk) :: A_coef_s_1
contains

    subroutine F1_hyperfine_intensity
        implicit none
        INTEGER(ik) :: index_F1_bra, index_F1_ket, &
                       index_represCs_bra, index_represCs_ket
        INTEGER(ik) :: Ndimen_F1_bra, Ndimen_F1_ket, &
                       Nlevels_F1_bra, Nlevels_F1_ket
        REAL(rk), ALLOCATABLE :: parity_conserved_F1_reduced_TDM_matrix(:,:), &
                                 parity_conserved_F1_transitions_matrix(:,:)

        write(out, '(/A)') "Start: hyperfine transitions calculation"
        open(unit=unit_hyperfine_transitions, file="hyperfine.trans")

        A_coef_s_1 = 64.0E-36 * pi**4  / (3.0_rk * planck)

        ! write(unit_hyperfine_transitions, "(A12, A12, A22, A22, A22)") &
        !     "N_upper", "N_lower", "Einstein-A [s-1]", "nu [cm-1]", "S [Debye^2]"

        do index_F1_bra = 1, num_F1
            write(out, '(/A4, A, F6.1)') '', 'Start: F_bra =', F1_list(index_F1_bra)

            Ndimen_F1_bra = primitive_F1_basis(index_F1_bra)%Ndimen  

            do index_F1_ket = max(1, index_F1_bra - 1), min(num_F1, index_F1_bra + 1) ! F1 selection rules
                if ((index_F1_bra == index_F1_ket) .and. F1_list(index_F1_bra) == 0) cycle

                write(out, '(/A6, A, F6.1)') '', 'Start: F_ket =', F1_list(index_F1_ket)
                Ndimen_F1_ket = primitive_F1_basis(index_F1_ket)%Ndimen

                write(out, '(/A8, A)') '',&
                    'Construct primitive reduced transition dipole moment matrix'
                ALLOCATE(primitive_F1_reduced_TDM_matrix(Ndimen_F1_bra, Ndimen_F1_ket))
                CALL construct_primitive_F1_reduced_TDM_matrix( &
                        index_F1_bra, index_F1_ket)
                write(out, '(A8, A/)') '', '... done'

                do index_represCs_bra = 1, sym%NrepresCs
                    ! parity selection rule
                    index_represCs_ket = mod(index_represCs_bra, sym%NrepresCs) + 1 

                    Nlevels_F1_bra = eigen_all_F1(index_F1_bra, index_represCs_bra)%Nlevels
                    Nlevels_F1_ket = eigen_all_F1(index_F1_ket, index_represCs_ket)%Nlevels

                    write(out, '(/A10, A)') '',&
                        'Construct parity conserved reduced transition dipole moment matrix'
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
            !     matmul(matmul(transpose(eigen_all_F1(index_F1_bra, index_represCs_bra)%vect), &
            !                 primitive_F1_reduced_TDM_matrix), &
            !         eigen_all_F1(index_F1_ket, index_represCs_ket)%vect)

            ALLOCATE(intermediate_matrix(Nlevels_F1_bra, Ndimen_F1_ket))

            ! M_indermediate = {Psi^{tau, F}}^\dagger * ^{F}M^{F'}
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

            lower_energy_lower_bound = intensity%erange_low(1) + min_eigen_value_F1
            lower_energy_upper_bound = intensity%erange_low(2) + min_eigen_value_F1
            upper_energy_lower_bound = intensity%erange_upp(1) + min_eigen_value_F1
            upper_energy_upper_bound = intensity%erange_upp(2) + min_eigen_value_F1
            frequency_lower_bound = max(0.0_rk, intensity%freq_window(1))
            frequency_upper_bound = intensity%freq_window(2)
            
            do ket = 1, Nlevels_F1_ket
                do bra = 1, Nlevels_F1_bra

                    if (eigen_all_F1(index_F1_ket, index_represCs_ket)%val(ket) &
                        < lower_energy_lower_bound) cycle
                    if (eigen_all_F1(index_F1_ket, index_represCs_ket)%val(ket) &
                        > lower_energy_upper_bound) cycle

                    if (eigen_all_F1(index_F1_bra, index_represCs_bra)%val(bra) &
                        < upper_energy_lower_bound) cycle
                    if (eigen_all_F1(index_F1_bra, index_represCs_bra)%val(bra) &
                        > upper_energy_upper_bound) cycle

                    ! The bra state is the final state
                    parity_conserved_F1_transitions_matrix(bra, ket) = &
                        eigen_all_F1(index_F1_bra, index_represCs_bra)%val(bra) &
                            - eigen_all_F1(index_F1_ket, index_represCs_ket)%val(ket)
                    
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
                        / (2.0_rk * primitive_F1_basis(index_F1_bra)%icontr(1)%F1 + 1.0_rk)

                    WRITE(unit_hyperfine_transitions, "(I12, I12, E22.12, F22.12)") &
                        eigen_all_F1(index_F1_bra, index_represCs_bra)%quanta(bra)%iroot, &
                        eigen_all_F1(index_F1_ket, index_represCs_ket)%quanta(ket)%iroot, &
                        EinsteinA, &
                        parity_conserved_F1_transitions_matrix(bra, ket)
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
        type(fieldT), pointer :: field

        primitive_F1_reduced_TDM_matrix = 0.0_rk

        Ndimen_F1_bra = primitive_F1_basis(index_F1_bra)%Ndimen
        Ndimen_F1_ket = primitive_F1_basis(index_F1_ket)%Ndimen

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
    
end module F1_intensity