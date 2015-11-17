   program duo

    use timer
    use accuracy
    use refinement
    use dipole
!    use compilation_details, only: write_compilation_details
    use header_info, only: write_logo
    use diatom_module,only : duo_j0,verbose,job,readinput,map_fields_onto_grid,action !, check_and_set_atomic_data

    interface ! used for isatty
      function isatty(fd) bind(c)
        use iso_c_binding
        integer(c_int), value :: fd
        integer(c_int) :: isatty
      end function
    end interface


    if (isatty(0) /= 0) then
      write(out, *) "DUO was called without an input file. Exiting..."
      write(out, *) "The input file should be provided after a less-than sign."
      write(out, *) "Example: duo.exe < my_input.inp"
      stop
    end if


     !
     if (verbose>=3) call TimerStart('Duo')
     !
     ! --------------------------------------------------------
     !
     ! initializing some constants 
     !
     call accuracyInitialize
     !
     ! read input data
     !
     call ReadInput
     !
     if (verbose>=3) call write_logo
 !    if (verbose>=3) call write_compilation_details
     if (verbose>=3) call print_physical_constants
     !
     !call check_and_set_atomic_data(verbose)
     !

     if (action%fitting) then 
       !
       ! Here we map all fields onto the same grid 
       !call map_fields_onto_grid
       !
       ! Here we map all fields onto the same grid 
       call map_fields_onto_grid(verbose)
       !
       call define_jlist
       call sf_fitting
     endif 
     !
     if (action%intensity) then 
       !
       call map_fields_onto_grid(verbose)
       !
       !call define_jlist
       call dm_tranint
       !
       write(out, '(a)') '--End--'
       stop
       !
     endif 
     !
     select case (job%contraction)  
       !
     case ("VIB") 
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
 
