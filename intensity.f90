module intensity_module

  ! A single entry point for all radiative-transition calculations.
  !
  ! Electric- and magnetic-dipole calculations share the rank-1 pipeline in
  ! dipole.f90. Electric quadrupole remains a separate rank-2 pipeline.

  use accuracy,        only : out
  use diatom_module,   only : action
  use dipole,          only : rank1_tranint
  use quadrupole,      only : qm_tranint

  implicit none
  private

  public :: intensity_tranint

contains

  subroutine intensity_tranint

    implicit none

    ! INTENSITY is the umbrella action.  QUADRUPOLE and MAGDIPOLE select an
    ! operator; if neither is present, use the established electric-dipole
    ! default.  Two operator selectors in one run would be ambiguous because
    ! each pipeline owns its output files and normalisation.
    if (.not. action%intensity) then
      write(out,'(/a)') 'Error: intensity_tranint called without an INTENSITY action'
      stop 'Error: INTENSITY action is not active'
    endif

    if (action%quadrupole .and. action%magdipole) then
      write(out,'(/a)') 'Error: QUADRUPOLE and MAGDIPOLE cannot be selected together'
      stop 'Error: ambiguous intensity operator'
    endif

    if (action%quadrupole) then
      call qm_tranint
    else
      ! rank1_tranint configures the
      ! appropriate symmetry map, Einstein-A prefactor and linestrength
      ! procedure from action%magdipole before entering the transition loops.
      call rank1_tranint
    endif

  end subroutine intensity_tranint

end module intensity_module
