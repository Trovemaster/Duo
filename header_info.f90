module header_info
use accuracy
implicit none

public :: write_logo

 contains

subroutine write_logo

    ! prints compiler logo, name, options,  etc.
    write(out,'(a)') " ____                                        ____  _   _  ___"
    write(out,'(a)') "|  _ \ _ __ ___   __ _ _ __ __ _ _ __ ___   |  _ \| | | |/ _ \"
    write(out,'(a)') "| |_) | '__/ _ \ / _` | '__/ _` | '_ ` _ \  | | | | | | | | | |"
    write(out,'(a)') "|  __/| | | (_) | (_| | | | (_| | | | | | | | |_| | |_| | |_| |"
    write(out,'(a)') "|_|   |_|  \___/ \__, |_|  \__,_|_| |_| |_| |____/ \___/ \___/"
    write(out,'(a)') "                 |___/"
    write(out,'(a)') 
  !  write(out,'(a)') 'Version: 20 December 2015' 
  !  write(out,'(a)') 'Compiled with: Intel(R) Fortran Intel(R) 64 Compiler XE for applications ' & 
  !                   // 'running on Intel(R) 64, Version 12.1.2.273 Build 20111128'
  !  write(out,'(a)') 'Compilation flags: -O0 -ip -openmp -m32 -mkl=parallel -static'
    write(out,'(a)')
    write(out,'(a)') "Please refer to:"
    write(out,'(a)') " Sergei N. Yurchenko, Lorenzo Lodi, Jonathan Tennyson and Andrey V. Stolyarov, "
    write(out,'(a)') " `DUO: A general program for calculating spectra of diatomic molecules',"
    write(out,'(a)') "  Computer Physics Communication, Volume 202, May 2016, pages 262-275"
    write(out,'(a)') "  Contacts: s.yurchenko@ucl.ac.uk; l.lodi@ucl.ac.uk; j.tennyson@ucl.ac.uk; "
    write(out,'(a)') "            avstol@phys.chem.msu.ru"
!     write(out,'(a)') "University College London, Gower Street, London WC1 6BT, United Kingdom"
    write(out,'(a)') "  http://arxiv.org/abs/1601.06531"
    write(out,'(a)') "  Check https://github.com/Trovemaster for the latest version!"
    write(out,'(a)')

end subroutine write_logo

end module header_info
