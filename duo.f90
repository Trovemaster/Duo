   program duo

    use timer
    use accuracy
    use refinement
    use dipole
    use quadrupole
    use RWF
    !use polarizability
!    use compilation_details, only: write_compilation_details
    use header_info, only: write_logo
    use diatom_module,only : duo_j0,verbose,job,readinput,map_fields_onto_grid,action !, check_and_set_atomic_data

    use F1_hyperfine, only: F1_hyperfine_structrure
    use F1_intensity, only: F1_hyperfine_intensity
    !use F1_fitting, only: F1_refinement_init, F1_refine

    interface ! used for isatty
      function isatty(fd) bind(c)
        use iso_c_binding
        integer(c_int), value :: fd
        integer(c_int) :: isatty
      end function
    end interface

    !
     if (verbose>=3) call write_logo
 !    if (verbose>=3) call write_compilation_details
    !
    if (isatty(0) /= 0) then
      write(out, *) "DUO was called without an input file. Exiting..."
      write(out, *) "The input file should be provided after a less-than sign."
      write(out, *) "Example: duo.exe < my_input.inp"
      stop
    end if
     !
     if (verbose>=3) call print_physical_constants
     !
     if (verbose>=3) call TimerStart('Duo')
     !
     ! --------------------------------------------------------
     !
     ! initializing some constants (superfluous?) 
     !
     !call accuracyInitialize
     !
     ! read input data
     !
     call ReadInput
     !
     !call check_and_set_atomic_data(verbose)
     !

     if(action%hyperfine) then
        action%save_eigen_J = .true.
        job%basis_set = 'KEEP'
        if (action%intensity) then            
            call duo_j0
            call F1_hyperfine_structrure(verbose)
            call F1_hyperfine_intensity
            stop
        else if (action%fitting) then
            write(out, '(a)') 'Fitting and hyperfine should not be used at the same time'
            stop 'please switch-off/remove either fitting or hyperfine'
        ! else if (action%fitting) then
        !     call map_fields_onto_grid(verbose)    
        !     call F1_refinement_init
        !     CALL F1_refine
        !     stop
        else
            call duo_j0
            call F1_hyperfine_structrure(verbose)
        endif
        write(out, '(a)') '--End--'
        stop
     endif   

     
     if (action%fitting) then 
        !
        if (action%intensity) then 
           write(out, '(a)') 'Fitting and intensity should not be used at the same time'
           stop 'please switch-off/remove either fitting or intensity'
        endif 
        !
        ! Here we map all fields onto the same grid 
        call map_fields_onto_grid(verbose)
        !
        call define_jlist
        call sf_fitting
      endif 
      ! 

     if (action%intensity) then 
        action%save_eigen_J = .true.
       !
       call map_fields_onto_grid(verbose)
       !
       if (action%raman) then 
         !call pol_tranint
         !
         write(out, '(a)') '--End--'
         stop
         !
       endif
       ! 
       !call define_jlist
       if (action%quadrupole) then
         call qm_tranint
         write(out, '(a)') '--End--'
         stop
       elseif(action%RWF) then
         call Raman_wavefunction
         write(out, '(a)') '--End--'
         stop
       else
         call dm_tranint
         write(out, '(a)') '--End--'
         stop
       endif
       !
       write(out, '(a)') '--End--'
       stop
       !
     endif
     !
     select case (job%contraction)  
       !
     case ("VIB","OMEGA") 
       !
       call duo_j0
       !
     case ("SIGMA-LAMBDA","ROT") 
       !
       write(out,'("This contraction scheme is not properly tested ",a)') trim(job%contraction)
       stop 'Contraction scheme is not ready'
       !
       !call duo_lambdasigma
       !
     case default 
       !
       write(out,'("Illegal contraction scheme ",a)') trim(job%contraction)
       stop 'Illegal contraction scheme'
       !
     end select 
     !
     if (verbose>=4) call MemoryReport
     !
     if (verbose>=3) call TimerStop('Duo')
     !
     if (verbose>=4) call TimerReport
     !
    end program duo
 
