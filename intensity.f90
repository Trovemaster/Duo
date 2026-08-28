module intensity_module

  ! A single entry point for all radiative-transition calculations.
  !
  ! E1, M1 and E2 calculations share the rank-generic pipeline in dipole.f90.

  use accuracy,        only : out
  use diatom_module,   only : action
  use dipole,          only : intensity_pipeline_tranint

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

    ! intensity_pipeline_tranint configures E1, M1 or E2 tensor rank, symmetry map,
    ! Einstein-A normalisation, RichMol adapter and numerical linestrength
    ! procedure before entering the common state and transition loops.
    call intensity_pipeline_tranint

  end subroutine intensity_tranint

end module intensity_module
