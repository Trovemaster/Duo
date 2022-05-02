module F1_fitting
    use accuracy
    use F1_hyperfine, only : eigen_all_F1, F1_global_min, &
                             F1_hyperfine_structrure, min_eigen_value_F1, &
                             F1_list, primitive_F1_basis
    use diatom_module, only : action, job, duo_j0, fieldmap, fieldT, &
        grid, fitting, linkT, &
        poten, spinorbit, l2, lxly, spinspin, spinspino, spinrot, &
        bobrot, diabatic, lambdaopq, lambdap2q, lambdaq, &
        l_omega_obj, s_omega_obj, quadrupoletm, hfcc1, &
        Nobjects, abinitio, eigen, basis
    ! use refinement, only: fit_indexT, object_containerT
    
    ! use nlopt_wrap, only : nlopt_opt, nlopt_func, create, destroy
    ! use nlopt_enum, only : NLOPT_SUCCESS, algorithm_from_string

    implicit none
        
    include "C:\Users\Leopard\nlopt\x64-Release\include\nlopt.f"

    type fit_indexT
        integer(ik)  :: i
        integer(ik)  :: iobject
        integer(ik)  :: iterm
        integer(ik)  :: ifield
    end type fit_indexT
  
    type object_containerT
        type(fieldT), pointer :: field
    end type object_containerT

    type states_containerT
        real(rk) :: term_value
        real(rk) :: F1 
        integer(ik) :: index_represCs
        integer(ik) :: iF1_ID
    endtype states_containerT

    integer(ik), parameter :: unit_hyperfine_states_fit = 67
    integer(ik) :: num_F1_hyperfine_states_fit, evalation_step = 0_ik

    type(states_containerT), ALLOCATABLE ::  F1_hyperfine_states_fit(:)
    real(rk), allocatable :: measured_F1_hyperfine_energies(:), &
                             calculated_F1_hyperfine_energies(:)
    integer(ik) :: num_parameters, Nobjectmax

    type(fit_indexT), allocatable :: fit_index(:)
    type(object_containerT), allocatable :: objects(:, :)

    ! type :: constraint_data
    !     real(rk) :: d(2)
    ! end type
    
contains

    subroutine F1_refine()
        implicit none

        real(rk) :: lb(num_parameters), ub(num_parameters)
        integer*8 opt
        ! double precision d1(2), d2(2)
        real(rk) :: x(num_parameters), minf, xtol
        integer(ik) :: ires, i

        write(out, '(/A)') 'Start hyperfine refinement'
        
        ! call nlo_create(opt, NLOPT_LD_MMA, num_parameters)
        call nlo_create(opt, NLOPT_LN_COBYLA, num_parameters)


        call nlo_set_min_objective(ires, opt, objective_function_F1_hyperfine_fit, 0)
        
        ! d1(1) = 2.
        ! d1(2) = 0.
        ! call nlo_add_inequality_constraint(ires, opt, myconstraint, d1, 1.D-8)
        ! d2(1) = -1.
        ! d2(2) = 1.
        ! call nlo_add_inequality_constraint(ires, opt, myconstraint, d2, 1.D-8)
        
        call nlo_set_xtol_rel(ires, opt, 1.0E-4)
        call nlo_set_ftol_rel(ires, opt, 1.0E-4);
        ! call nlo_set_ftol_abs(ires, opt, double tol)
        call nlo_set_stopval(ires, opt, 100);
        call nlo_set_maxeval(ires, opt, 10);
        ! call nlo_set_initial_step(ires, opt, dx)
        
        x = get_initial_guess()

        ! write(out, '(/A)') 'Hyperfine refinement initial guess:'
        ! do i = 1, num_parameters
        !     write(out, '(A2, A, I3, A, F24.14)') '', 'x(', i, ') = ', x(i) 
        ! enddo

        xtol = 0.01
        lb = x * (1 - xtol)
        ub = x * (1 + xtol)
        call nlo_set_lower_bounds(ires, opt, lb)
        call nlo_set_upper_bounds(ires, opt, ub)

        call nlo_optimize(ires, opt, x, minf)

        if (ires < 0) then
            write(out, '(/A2, A)') '', trim(nlopt_flag(ires))
        else
            write(out, '(/A)') 'NLOPT found minimum at'
            do i = 1, num_parameters
                write(out, '(A2, A, I3, A, F24.14)') '', 'x(', i, ') = ', x(i) 
            enddo
            write(out, '(A, F24.14)') 'The minimum is ', minf
            write(out, '(A, A)') 'The NLOPT flag is ', trim(nlopt_flag(ires))
        endif
        
        call nlo_destroy(opt)

    contains
        function nlopt_flag(ires) result(flag)
            implicit none
            integer(ik) :: ires
            character(100):: flag

            select case(ires)
            case (NLOPT_SUCCESS)
                flag = "NLOPT_SUCCESS"
            case (NLOPT_STOPVAL_REACHED)
                flag = "NLOPT_STOPVAL_REACHED"
            case (NLOPT_FTOL_REACHED)
                flag = "NLOPT_FTOL_REACHED"
            case (NLOPT_XTOL_REACHED)
                flag = "NLOPT_XTOL_REACHED"
            case (NLOPT_MAXEVAL_REACHED)
                flag = "NLOPT_MAXEVAL_REACHED"
            case (NLOPT_MAXTIME_REACHED)
                flag = "NLOPT_MAXTIME_REACHED ="
            case (NLOPT_FAILURE)
                flag = "NLOPT_FAILURE"
            case (NLOPT_INVALID_ARGS)
                flag = "NLOPT_INVALID_ARGS"
            case (NLOPT_OUT_OF_MEMORY )
                flag = "NLOPT_OUT_OF_MEMORY "
            case (NLOPT_ROUNDOFF_LIMITED)
                flag = "NLOPT_ROUNDOFF_LIMITED"
            case (NLOPT_FORCED_STOP)
                flag = "NLOPT_FORCED_STOP"

            case default
                flag = "Unknown NLOPT flag"
            endselect
            
        end function nlopt_flag
    
    end subroutine F1_refine

    function get_initial_guess() result(x)
        implicit none
        real(rk) :: x(num_parameters)

        INTEGER(ik) :: ifitpar 
        ! integer(ik)  :: i
        integer(ik)  :: iobject
        integer(ik)  :: iterm
        integer(ik)  :: ifield

        do ifitpar = 1, num_parameters
            ! i = fit_index(ifitpar)%i
            iobject = fit_index(ifitpar)%iobject
            ifield = fit_index(ifitpar)%ifield
            iterm = fit_index(ifitpar)%iterm

            x(ifitpar) = objects(iobject,ifield)%field%value(iterm) 
        enddo
        
    end function get_initial_guess

    subroutine F1_refinement_init()
        implicit none

        INTEGER(ik) :: iobject, ifield, iterm, j, Nfields
        INTEGER(ik) :: Nterms, ifield_, ipotmin, maxNfields, info, parmax
        real(rk) :: req
       
        action%save_eigen_J = .true.
        job%basis_set = 'KEEP'
        
        call read_F1_hyperfine_states_fit

        Nobjectmax = Nobjects - 3

        ipotmin = minloc(poten(1)%gridvalue,dim=1) 
        req = grid%r(ipotmin)

        maxNfields = 1
        do iobject =1,Nobjectmax
            maxNfields = max(fieldmap(iobject)%Nfields, maxNfields)
        enddo

        allocate(objects(Nobjectmax, maxNfields), stat=info)

        ifield_ = 0

        do iobject =1, Nobjectmax
            Nfields = fieldmap(iobject)%Nfields
            do ifield =1,Nfields
                select case (iobject)
                    case (1)
                        objects(iobject,ifield)%field => poten(ifield)
                    case (2)
                        objects(iobject,ifield)%field => spinorbit(ifield)
                    case (3)
                        objects(iobject,ifield)%field => l2(ifield)
                    case (4)
                        objects(iobject,ifield)%field => lxly(ifield)
                    case (5)
                        objects(iobject,ifield)%field => spinspin(ifield)
                    case (6)
                        objects(iobject,ifield)%field => spinspino(ifield)
                    case (7)
                        objects(iobject,ifield)%field => bobrot(ifield)
                    case (8)
                        objects(iobject,ifield)%field => spinrot(ifield)
                    case (9)
                        objects(iobject,ifield)%field => diabatic(ifield)
                    case (10)
                        objects(iobject,ifield)%field => lambdaopq(ifield)
                    case (11)
                        objects(iobject,ifield)%field => lambdap2q(ifield)
                    case (12)
                        objects(iobject,ifield)%field => lambdaq(ifield)
                    case (21, 22, 23, 24, 25, 26, 27)
                        objects(iobject,ifield)%field => hfcc1(iobject - 20)%field(ifield)
                    case (13)
                        objects(iobject,ifield)%field => quadrupoletm(ifield)
                end select
                
                Nterms = objects(iobject,ifield)%field%Nterms
                !
                ifield_ = ifield_ + 1
                !
                ! change to the equilibrium
                !
                if (abinitio(ifield_)%Nterms==1) then
                    !
                    abinitio(ifield_)%grid(1) = req
                    !
                    !f = analytical_field(req,objects(iobject,ifield)%field%type,objects(iobject,ifield)%field%value)
                    !
                    abinitio(ifield_)%gridvalue(1) = objects(iobject,ifield)%field%gridvalue(ipotmin)
                    
                endif 
            enddo
        enddo

        parmax = fitting%parmax
        ALLOCATE(fit_index(parmax))

        num_parameters = 0
        ifield_ = 0
        j = 0     

        do iobject =1, Nobjectmax
            Nfields = fieldmap(iobject)%Nfields
            do ifield =1, Nfields
                do iterm=1, objects(iobject,ifield)%field%Nterms
                    j = j + 1
                    if (nint(objects(iobject,ifield)%field%weight(iterm)) > 0 .and. & 
                        trim(objects(iobject,ifield)%field%type)/='GRID') then 
                        
                        num_parameters=num_parameters+1
                        
                        fit_index(num_parameters)%i = j
                        fit_index(num_parameters)%iobject = iobject
                        fit_index(num_parameters)%ifield = ifield
                        fit_index(num_parameters)%iterm = iterm
                        
                    endif
                enddo
            enddo
        enddo
    
    end subroutine F1_refinement_init

    subroutine read_F1_hyperfine_states_fit
        use input

        implicit none

        integer(ik) :: parity, index_represCs, iF1_ID
        integer(ik) :: index_F1_hyperfine_state, open_state, read_state
        real(rk) :: term_value, F1
        type(states_containerT), allocatable :: temp_states_array(:)

        num_F1_hyperfine_states_fit = 1024


        allocate(F1_hyperfine_states_fit(num_F1_hyperfine_states_fit))
        
        open(unit=unit_hyperfine_states_fit, &
            file="hyperfine_states_fit.txt", &
            action="read", &
            iostat=open_state)

            index_F1_hyperfine_state = 1
            do
                read(unit=unit_hyperfine_states_fit, &
                    fmt = "(F24.14, F10.1, I10, I10)", &
                    iostat=read_state) &
                    term_value, F1, parity, iF1_ID
                if (read_state /= 0) exit
                if (term_value < 0) cycle

                F1_hyperfine_states_fit(index_F1_hyperfine_state)%term_value = term_value
                F1_hyperfine_states_fit(index_F1_hyperfine_state)%F1 = F1
                if (parity == 1) then
                    index_represCs = 0_ik
                else
                    index_represCs = 1_ik
                end if
                F1_hyperfine_states_fit(index_F1_hyperfine_state)%index_represCs = index_represCs
                F1_hyperfine_states_fit(index_F1_hyperfine_state)%iF1_ID = iF1_ID

                index_F1_hyperfine_state = index_F1_hyperfine_state + 1

                if (index_F1_hyperfine_state > num_F1_hyperfine_states_fit) then
                    allocate(temp_states_array(num_F1_hyperfine_states_fit))
                    temp_states_array = F1_hyperfine_states_fit
                    deallocate(F1_hyperfine_states_fit)                    
                    allocate(F1_hyperfine_states_fit(num_F1_hyperfine_states_fit * 2))
                    F1_hyperfine_states_fit(1:num_F1_hyperfine_states_fit) &
                        = temp_states_array 
                    deallocate(temp_states_array)
                    num_F1_hyperfine_states_fit = num_F1_hyperfine_states_fit * 2
                end if
            enddo
        close(unit=unit_hyperfine_states_fit)

        index_F1_hyperfine_state = index_F1_hyperfine_state - 1

        if (index_F1_hyperfine_state < num_F1_hyperfine_states_fit) then
            allocate(temp_states_array(index_F1_hyperfine_state))
            temp_states_array = &
                F1_hyperfine_states_fit(1:index_F1_hyperfine_state)
            deallocate(F1_hyperfine_states_fit)                    
            allocate(F1_hyperfine_states_fit(index_F1_hyperfine_state))
            F1_hyperfine_states_fit = temp_states_array 
            deallocate(temp_states_array)
            num_F1_hyperfine_states_fit = index_F1_hyperfine_state
        end if

        allocate(measured_F1_hyperfine_energies(num_F1_hyperfine_states_fit))

        do index_F1_hyperfine_state = 1, num_F1_hyperfine_states_fit
            term_value = F1_hyperfine_states_fit(index_F1_hyperfine_state)%term_value
            F1 = F1_hyperfine_states_fit(index_F1_hyperfine_state)%F1
            index_represCs = F1_hyperfine_states_fit(index_F1_hyperfine_state)%index_represCs
            iF1_ID = F1_hyperfine_states_fit(index_F1_hyperfine_state)%iF1_ID
            
            measured_F1_hyperfine_energies(index_F1_hyperfine_state) = term_value
        enddo
    
    end subroutine read_F1_hyperfine_states_fit

    subroutine forward(x)
        implicit none

        real(rk), INTENT(IN) :: x(:)

        call update_refined_parameters(x)
        call update_linked_parameters
        call duo_j0(0_ik)
        call F1_hyperfine_structrure(0_ik)
    
    end subroutine forward

    subroutine update_refined_parameters(x)
        implicit none

        real(rk), INTENT(IN) :: x(:)
        INTEGER(ik) :: ifitpar 
        ! integer(ik)  :: i
        integer(ik)  :: iobject
        integer(ik)  :: iterm
        integer(ik)  :: ifield

        do ifitpar = 1, num_parameters
            ! i = fit_index(ifitpar)d%i
            iobject = fit_index(ifitpar)%iobject
            ifield = fit_index(ifitpar)%ifield
            iterm = fit_index(ifitpar)%iterm

            objects(iobject,ifield)%field%value(iterm) = x(ifitpar)
        enddo

    end subroutine update_refined_parameters

    subroutine update_linked_parameters
        implicit none

        integer(ik) :: ifield,iterm,iobject,Nfields
        type(linkT),pointer :: flink

        do iobject =1, Nobjectmax
            Nfields = fieldmap(iobject)%Nfields
            do ifield =1, Nfields
                ! ifield_ = ifield_ + 1
                do iterm = 1, objects(iobject,ifield)%field%Nterms
                    flink => objects(iobject,ifield)%field%link(iterm)
                    if (flink%iobject/=0) then
                        objects(iobject,ifield)%field%value(iterm) = &
                            objects(flink%iobject,flink%ifield)%field%value(flink%iparam)
                    endif 
                enddo
            enddo
        enddo
    end subroutine update_linked_parameters

    subroutine objective_function_F1_hyperfine_fit( &
        objF_val, n, x, gradient, need_gradient, func_data)
        implicit none

        integer, intent(in) :: n, need_gradient
        real(rk), intent(inout) :: x(n), gradient(n)
        real(rk), intent(out) :: objF_val
        class(*), intent(in), optional :: func_data

        real(rk) :: err(num_F1_hyperfine_states_fit)
        integer(ik) :: i

        if (need_gradient .ne. 0) then
            gradient = 0.0_rk
        endif

        evalation_step = evalation_step + 1
        write(out, '(/A, I5)') 'Step ', evalation_step

        call forward(x)
        call extract_calculated_F1_hyperfine_energies

        err = calculated_F1_hyperfine_energies - measured_F1_hyperfine_energies
        objF_val = sum(err * err)
        
        write(out, '(A, F24.14)') 'The objective function value is', objF_val
        do i = 1, num_parameters
            write(out, '(A2, A, I3, A, F24.14)') '', 'x(', i, ') = ', x(i) 
        enddo
    
        call deallocate_intermediate_variables

    end subroutine objective_function_F1_hyperfine_fit

    subroutine deallocate_intermediate_variables
        deallocate(basis)
        deallocate(eigen)
        deallocate(eigen_all_F1)
        deallocate(F1_list)
        deallocate(primitive_F1_basis)   
        deallocate(calculated_F1_hyperfine_energies)
    end subroutine deallocate_intermediate_variables

    subroutine myconstraint(val, n, x, grad, need_gradient, d)
        integer need_gradient, n
        double precision val, x(n), grad(n), d(2), a, b
        a = d(1)
        b = d(2)
        if (need_gradient.ne.0) then
           grad = 0
        endif
        val = (a*x(1) + b)**3 - x(2)
    end

    subroutine extract_calculated_F1_hyperfine_energies
        implicit none

        integer(ik) :: iF1_ID, index_represCs, index_F1, index_F1_hyperfine_state
        real(rk) :: term_value, F1

        allocate(calculated_F1_hyperfine_energies(num_F1_hyperfine_states_fit))   

        do index_F1_hyperfine_state = 1, num_F1_hyperfine_states_fit
            term_value = F1_hyperfine_states_fit(index_F1_hyperfine_state)%term_value
            F1 = F1_hyperfine_states_fit(index_F1_hyperfine_state)%F1
            index_represCs = F1_hyperfine_states_fit(index_F1_hyperfine_state)%index_represCs
            iF1_ID = F1_hyperfine_states_fit(index_F1_hyperfine_state)%iF1_ID

            index_F1 = 1 + int(F1 - F1_global_min)
            
            calculated_F1_hyperfine_energies(index_F1_hyperfine_state) &
                = eigen_all_F1(index_F1, index_represCs + 1)%val(iF1_ID)
        end do

        calculated_F1_hyperfine_energies = calculated_F1_hyperfine_energies - min_eigen_value_F1
    
    end subroutine extract_calculated_F1_hyperfine_energies

end module F1_fitting