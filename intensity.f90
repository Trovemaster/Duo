module intensity_module

  ! A single entry point for all radiative-transition calculations.
  !
  ! The three operators also differ in their rotational tensors, selection
  ! rules and Einstein-A prefactors.  Keep those numerical pipelines in
  ! their owning modules and centralise the action dispatch here.  This
  ! preserves existing results and provides one boundary from which common
  ! pipeline code can subsequently be extracted safely.

  use accuracy,        only : out
  use diatom_module,   only : action
  use dipole,          only : dm_tranint
  use quadrupole,      only : qm_tranint
  use magnetic_dipole, only : md_tranint

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
    elseif (action%magdipole) then
      call md_tranint
    else
      call dm_tranint
    endif

  end subroutine intensity_tranint

end module intensity_module
