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
    write(out,'(a)') "Please refer to:"
    write(out,'(a)') " Sergei N. Yurchenko, Lorenzo Lodi, Jonathan Tennyson and Andrey V. Stolyarov, "
    write(out,'(a)') " `DUO: a general program for calculating spectra of diatomic molecules',"
    write(out,'(a)') "  Computer Physics Communication, (to be submitted), 2015."
    write(out,'(a)') "  Contacts: s.yurchenko@ucl.ac.uk; l.lodi@ucl.ac.uk; j.tennyson@ucl.ac.uk; "
    write(out,'(a)') "            avstol@phys.chem.msu.ru"
!     write(out,'(a)') "University College London, Gower Street, London WC1 6BT, United Kingdom"
    write(out,'(a)')

end subroutine write_logo

end module header_info
