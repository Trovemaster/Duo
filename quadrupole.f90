module quadrupole

  use accuracy,      only : ik, rk, small_
  use diatom_module, only : Intensity, basis, poten, three_j, &
                            nQuadrupoles, quadrupoletm
  use timer,         only : TimerStart, TimerStop

  implicit none
  private

  public :: qm_do_1st_half_linestrength, qm_do_LF_matrix_elements

contains

  subroutine qm_do_LF_matrix_elements(iLFa, iLFb, iunit, jI, jF, icount)
    ! calculates the components of the matrix element due to the
    ! molecule to laboratory frame transformation (i.e rotational part)
    ! when called by the shared intensity pipeline

    implicit none

    ! I/O variables
    real(rk), intent(in)                :: jI, jF
    integer(ik), intent(in)             :: iLFa, iLFb, iunit
    integer(ik), intent(out), optional  :: icount

    ! quantum numbers and indexes
    integer(ik)                         :: indJf, indJi, indMf, indMi
    integer(ik)                         :: Mf_, Mi_
    real(rk)                            :: Mf, Mi

    ! calculation variables
    real(rk)                            :: T(3, 3, 0:2, -2:2)
    integer(ik)                         :: m
    real(rk)                            :: matElem, f3j

    ! assorted
    integer(ik)                         :: icount_

    ! The matrix T gives the transformation from the spherical,
    ! irreducible representation to the Cartesian representation.
    ! Entries corresponding to rank zero and one components are zero
    ! as the quadrupole moment is zero for these in diatomic molecules.
    !
    ! Imaginary values are made real here as this does not change the
    ! value of the line strength, but we must indicate when there is an
    ! additional phase factor in the RichMol file.

    T = 0.0_rk

    ! real
    T(1, 1, 2, 2) = 0.5_rk
    T(1, 1, 2, 0) = -1.0_rk/sqrt(6.0_rk)
    T(1, 1, 2, -2) = 0.5_rk

    ! real
    T(2, 2, 2, 2) = -0.5_rk
    T(2, 2, 2, 0) = -1.0_rk/sqrt(6.0_rk)
    T(2, 2, 2, -2) = -0.5_rk

    ! real
    T(3, 3, 2, 0) = 2.0_rk/sqrt(6.0_rk)

    ! imaginary
    T(1, 2, 2, 2) = -0.5_rk
    T(1, 2, 2, -2) = 0.5_rk

    ! imaginary
    T(2, 1, 2, 2) = -0.5_rk
    T(2, 1, 2, -2) = 0.5_rk

    ! real
    T(1, 3, 2, 1) = -0.5_rk
    T(1, 3, 2, -1) = 0.5_rk

    ! real
    T(3, 1, 2, 1) = -0.5_rk
    T(3, 1, 2, -1) = 0.5_rk

    ! imaginary
    T(2, 3, 2, 1) = 0.5_rk
    T(2, 3, 2, -1) = 0.5_rk

    ! imaginary
    T(3, 2, 2, 1) = 0.5_rk
    T(3, 2, 2, -1) = 0.5_rk

    ! initialise variables
    icount_ = 0
    matElem = 0
    indJf = nint(2*jF) ! for indexing loops (jF/i may be 1/2 integer)
    indJi = nint(2*jI)

    call TimerStart('do_LF_matrix_elements')

    loop_F :  do indMf = -indJf, indJf, 2

      Mf = 0.5_rk * indMf

      ! remove imaginary factor due to (-)^Mf when jF is half-integer

      ! this method requires only +i phase factor later
      if ( mod(indJf + 1, 2) == 0 ) then
        Mf_ = nint(0.5_rk*indMf - 0.5_rk)
      else
        Mf_ = nint(0.5_rk*indMf)
      endif

      loop_I : do indMi = -indJi, indJi, 2

        Mi = 0.5_rk * indMi

        ! remove imaginary phase due to (-)^Mi when jI is half-integer

        ! this method requires only +i phase factor later
        if ( mod(indJi + 1, 2) == 0 ) then
          Mi_ = nint(0.5_rk*indMi - 0.5_rk)
        else
          Mi_ = nint(0.5_rk*indMi)
        endif

        ! M quantum number selection rule (from 3-j symbol)
        m = nint(Mi - Mf)
        if  ( abs(m) > 2 )  cycle

        ! calculate value of 3j symbol and then LF matrix elem
        f3j = three_j(jI, 2.0_rk, jF, Mi, real(m, rk), -Mf)
        matElem = T(iLFa,iLFb, 2, m) * (-1.0_rk)**Mf * f3j *&
                  sqrt( (2.0_rk*jI + 1.0_rk)*(2.0_rk*jF +  1.0_rk) )

        if ( abs(matElem) > small_ ) then
          icount_ = icount_ +  1
          if (.not. present(icount)) then
            ! write integerised Mi/f and LF matrix elem to file.=
            write (iunit, "(2(i6),3x,e24.14)") Mi_, Mf_, matElem
          endif
        endif
      enddo loop_I
    enddo loop_F

    if ( present(icount) ) icount = icount_

    call TimerStop('do_LF_matrix_elements')

  end subroutine qm_do_LF_matrix_elements

  subroutine qm_do_1st_half_linestrength(jI, jF, indI, indF, dimenI,&
                                         dimenF, vector, half_ls)
    ! calculates the first half of the linestrength when called by
    ! the shared intensity pipeline

    implicit none

    ! I/O  variables
    real(rk), intent(in)    :: jI, jF
    integer(ik), intent(in) :: indI, indF, dimenI, dimenF
    real(rk), intent(in)    :: vector(:)
    real(rk), intent(out)   :: half_ls(:)

    ! quantum numbers
    integer(ik)             :: istateI, istateF, ivibI, ivibF, &
                                ilambdaI, ilambdaF
    real(rk)                :: omegaI, omegaF, spinI, spinF, &
                                sigmaI, sigmaF

    ! psuedo quantum numbers and indexes
    integer(ik)             :: iomegaI_, iomegaF_, istateI_, istateF_, &
                                ilambdaI_, ilambdaF_
    real(rk)                :: spinI_, spinF_
    integer(ik)             :: icontrI, icontrF, indQuad, indPermute, &
                                indSigmaV

    ! calculation variables
    real(rk)                :: f3j, vibME, ls
    integer(ik)             :: itau
    !
    ! This kernel is called from the master pipeline's OpenMP region.  The
    ! global timer stack is not thread-safe, so timing belongs to the caller.
    half_ls = 0
    !
    !$omp parallel do private(&
    !$omp & icontrF,istateF,ivibF,ilambdaF,omegaF,spinF,sigmaF,iomegaF_,&
    !$omp & icontrI,istateI,ivibI,ilambdaI,omegaI,spinI,sigmaI,iomegaI_,f3j,ls,indQuad,indPermute,&
    !$omp & istateI_,ilambdaI_,spinI_,istateF_,ilambdaF_,spinF_,indSigmaV,vibME,itau) shared(half_ls) schedule(guided)
    loop_F : do icontrF = 1, dimenF

      istateF  = basis(indF)%icontr(icontrF)%istate
      ivibF    = basis(indF)%icontr(icontrF)%ivib
      ilambdaF = basis(indF)%icontr(icontrF)%ilambda
      omegaF  = basis(indF)%icontr(icontrF)%omega
      spinF   = basis(indF)%icontr(icontrF)%spin
      sigmaF  = basis(indF)%icontr(icontrF)%sigma

      ! remove imaginary factor (-1)^Omega when J is half-int
      iomegaF_ = nint(omegaF)
      if (mod(nint(2.0_rk*omegaF+1.0_rk),2)==0 ) iomegaF_ = nint((2.0_rk*omegaF-1.0_rk)*0.5_rk)         

      loop_I : do icontrI = 1, dimenI

        istateI  = basis(indI)%icontr(icontrI)%istate
        ivibI    = basis(indI)%icontr(icontrI)%ivib
        ilambdaI = basis(indI)%icontr(icontrI)%ilambda
        omegaI  = basis(indI)%icontr(icontrI)%omega
        spinI   = basis(indI)%icontr(icontrI)%spin
        sigmaI  = basis(indI)%icontr(icontrI)%sigma

        if (abs(nint(omegaF - omegaI))>2.or.nint(spinI-spinF)/=0.or.nint(sigmaI-sigmaF)/=0) cycle loop_I
        if (abs(nint(omegaF - omegaI))/=abs(ilambdaI-ilambdaF)) cycle loop_I
        
        iomegaI_ = int(omegaI)
        if (mod(nint(2.0_rk*omegaI+1.0_rk),2)==0 ) iomegaI_ = nint((2.0_rk*omegaI-1.0_rk)*0.5_rk)

        !remove imaginary factor (-1)^Omega when J is half-int
        iomegaI_ = int(omegaI)
        if (mod(nint(2.0_rk*omegaI+1.0_rk),2)==0 ) iomegaI_ = nint((2.0_rk*omegaI-1.0_rk)*0.5_rk)

        f3j = three_j(jI, 2.0_rk, jF, omegaI, omegaF - omegaI, -omegaF)
        if ( abs(f3j) < Intensity%threshold%coeff ) cycle loop_I

        ls = 0

        loop_quadpole :  do indQuad = 1, nQuadrupoles

          !field  => quadrupoletm(indQuad)

          ! we can calculate opposite matrix elements at same time
          do indPermute =  0, 1

            if  ( indPermute == 0 ) then
              istateI_  = quadrupoletm(indQuad)%istate  ; istateF_  = quadrupoletm(indQuad)%jstate
              ilambdaI_ = quadrupoletm(indQuad)%lambda  ; ilambdaF_ = quadrupoletm(indQuad)%lambdaj
              spinI_   = quadrupoletm(indQuad)%spini   ; spinF_   = quadrupoletm(indQuad)%spinj
            else
              istateI_  = quadrupoletm(indQuad)%jstate  ; istateF_  = quadrupoletm(indQuad)%istate
              ilambdaI_ = quadrupoletm(indQuad)%lambdaj ; ilambdaF_ = quadrupoletm(indQuad)%lambda
              spinI_   = quadrupoletm(indQuad)%spinj   ; spinF_   = quadrupoletm(indQuad)%spini
            endif

            ! do not double count diagonal elements
            if (      indPermute == 1 &
                .and. istateI_  == istateF_ &
                .and. ilambdaI_ == ilambdaF_ &
                .and. nint(spinI_ - spinF_) == 0 &
                ) cycle

            ! check the field I and F states match the basis loop
            if (     istateI_ /= istateI &
                .or. istateF_ /= istateF &
                ) cycle

            ! since only one of Lambda > or < 0 is given as input, we
            ! must also account for negative values. This is done using
            ! a lab-fixed inversion, equivalent to sigmav operation.
            do indSigmaV = 0, 1

              ! skip if both lambda = 0, else double counting
              if ( indSigmaV == 1 &
                  .and. abs(ilambdaF_) + abs(ilambdaI_) == 0 &
                  ) cycle

              ! apply sigmav transformation to build Lambda < 0 states
              ilambdaI_ = ilambdaI_ * (-1)**indSigmaV
              ilambdaF_ = ilambdaF_ * (-1)**indSigmaV

              ! check field states match basis states
              if (     ilambdaI_ /= ilambdaI &
                  .or. ilambdaF_ /= ilambdaF &
                  ) cycle

              ! reinforce lambda selection rule
              if  ( abs(ilambdaF - ilambdaI) > 2 ) cycle

              ! vibrational matrix element
              vibME = quadrupoletm(indQuad)%matelem(ivibI, ivibF)

              ! the Sigma electronic states (with lambda = 0) have
              ! definite +/- parities - thus, if one of the mat.
              ! elem. states is the Sigma^-, lambda = 0 state, the
              ! sign of the matrix element will change under the sigmaV
              ! inversion. This cancels out if both states are ungerade.
              if ( indSigmaV == 1 ) then

                itau = 0

                if (      ilambdaI_ == 0 &
                    .and. poten(istateI)%parity%pm == -1 &
                    ) itau = itau + 1

                if (      ilambdaF_ == 0 &
                    .and. poten(istateF)%parity%pm == -1 &
                    ) itau = itau + 1

                vibME = vibME * (-1.0_rk)**itau
              endif

              ! calculate linestrength and add to sum
              ls =  vibME * f3j * vector(icontrI)
              half_ls(icontrF) = half_ls(icontrF) + (-1.0_rk)**(iomegaI_) * ls

            enddo
          enddo
        enddo loop_quadpole
      enddo  loop_I
    enddo loop_F
    !$omp end parallel do

  end subroutine qm_do_1st_half_linestrength

end module quadrupole
