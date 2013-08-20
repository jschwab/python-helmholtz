! here is the tabular helmholtz free energy eos:
      subroutine call_helmeos(nrow, den, temp, abar, zbar)
      include 'implno.dek'
      include 'vector_eos.dek'

! tests the eos routine
!
! ionmax  = number of isotopes in the network
! xmass   = mass fraction of isotope i
! aion    = number of nucleons in isotope i
! zion    = number of protons in isotope i

      integer, intent(in) :: nrow
      double precision, intent(in), dimension(nrow) :: den, temp, abar, zbar
      !f2py INTEGER, INTENT(hide) :: nrow
      !f2py DOUBLE PRECISION, DIMENSION(nrow), INTENT(in) :: den, temp, abar, zbar

      integer :: i

! don't try and overfill the array
!      if (ninput.gt.nrowmax) then stop

! set the input vector. pipeline is only 1 element long

      jlo_eos = 1 ; jhi_eos = nrow
      do i = 1, nrow
         temp_row(i) = temp(i)
         den_row(i)  = den(i)
         abar_row(i) = abar(i)
         zbar_row(i) = zbar(i)
      end do

! read the data table and call the eos
      call read_helm_table
      call helmeos

      end


      subroutine call_helmeos_DP(nrow, den, pres, abar, zbar, tguess)
      include 'implno.dek'
      include 'vector_eos.dek'

! tests the eos routine
!
! ionmax  = number of isotopes in the network
! xmass   = mass fraction of isotope i
! aion    = number of nucleons in isotope i
! zion    = number of protons in isotope i

      integer, intent(in) :: nrow
      double precision, intent(in), dimension(nrow) :: den, pres, abar, zbar, tguess
      !f2py INTEGER, INTENT(hide) :: nrow
      !f2py DOUBLE PRECISION, DIMENSION(nrow), INTENT(in) :: den, pres, abar, zbar, tguess
      double precision, dimension(nrow) :: rerr_P, rerr_T
      double precision, dimension(nrow) :: delta_P, delta_T
      double precision, dimension(nrow) :: T_lower, T_upper
      double precision, dimension(nrow) :: Pgoal
      logical, dimension(nrow) :: NR_converged

      double precision, parameter :: temp_floor = 1e4
      double precision, parameter :: rtol = 1e-6

      integer :: i, iter
      integer, parameter :: max_iter  = 100

! don't try and overfill the array
!      if (ninput.gt.nrowmax) then stop

! read the data table
      call read_helm_table

! set the input vector. pipeline is only 1 element long

      jlo_eos = 1 ; jhi_eos = nrow
      do i = 1, nrow
         Pgoal(i) = pres(i)
         temp_row(i) = tguess(i)
         den_row(i)  = den(i)
         abar_row(i) = abar(i)
         zbar_row(i) = zbar(i)
      end do

      T_lower = temp_floor
      T_upper = 1e12
      NR_converged = .FALSE.

      ! do the NR iteration

      do iter = 1, max_iter

         call helmeos

         do i = 1, nrow

            ! if this point is converged, go to the next one
            if (NR_converged(i)) cycle

            ! energy difference
            delta_P(i) = Pgoal(i) - ptot_row(i)

            ! keep things safe with bisect-limits
            if (delta_P(i).gt.0) then
               t_lower(i) = temp_row(i)
            else
               t_upper(i) = temp_row(i)
            end if

            ! update temperature
            delta_T(i) = delta_P(i) / dpt_row(i)
            temp_row(i) = temp_row(i) + delta_T(i)

            ! if this took us out of bounds, don't let it happen
            ! choose a new point inside the interval [t_lower, t_upper]
            ! the point is in the middle of the interval (logarthmically)

            if ((temp_row(i).gt.t_upper(i)).OR.(temp_row(i).lt.t_lower(i))) then
               temp_row(i) = sqrt(t_lower(i) * t_upper(i))
            end if

            ! calculate relative errors
            rerr_P(i) = delta_P(i) / Pgoal(i)
            rerr_T(i) = delta_T(i) / temp_row(i)

            ! if we're at tolerances, end this
            if ((abs(rerr_P(i)).LE.rtol).AND.(abs(rerr_T(i)).LE.rtol)) then
               NR_converged(i) = .TRUE.
            endif

            !allow points at the temperature floor to "converge"
            if (t_upper(i).le.temp_floor * (1d0 + rtol)) then
               NR_converged(i) = .TRUE.
               temp_row(i) = temp_floor
            end if

         end do

         if (ALL(NR_converged)) exit

      end do

      ! once more, with feeling
      NR_converged = .FALSE.

      call helmeos

      end subroutine call_helmeos_DP

      subroutine call_helmeos_DS(nrow, den, entr, abar, zbar, tguess)
      include 'implno.dek'
      include 'vector_eos.dek'

! tests the eos routine
!
! ionmax  = number of isotopes in the network
! xmass   = mass fraction of isotope i
! aion    = number of nucleons in isotope i
! zion    = number of protons in isotope i

      integer, intent(in) :: nrow
      double precision, intent(in), dimension(nrow) :: den, entr, abar, zbar, tguess
      !f2py INTEGER, INTENT(hide) :: nrow
      !f2py DOUBLE PRECISION, DIMENSION(nrow), INTENT(in) :: den, entr, abar, zbar, tguess
      double precision, dimension(nrow) :: rerr_S, rerr_T
      double precision, dimension(nrow) :: delta_S, delta_T
      double precision, dimension(nrow) :: T_lower, T_upper
      double precision, dimension(nrow) :: Sgoal
      logical, dimension(nrow) :: NR_converged

      double precision, parameter :: temp_floor = 1e4
      double precision, parameter :: rtol = 1e-6

      integer :: i, iter
      integer, parameter :: max_iter  = 100

! don't try and overfill the array
!      if (ninput.gt.nrowmax) then stop

! read the data table
      call read_helm_table

! set the input vector. pipeline is only 1 element long

      jlo_eos = 1 ; jhi_eos = nrow
      do i = 1, nrow
         Sgoal(i) = entr(i)
         temp_row(i) = tguess(i)
         den_row(i)  = den(i)
         abar_row(i) = abar(i)
         zbar_row(i) = zbar(i)
      end do

      T_lower = temp_floor
      T_upper = 1e12
      NR_converged = .FALSE.

      ! do the NR iteration

      do iter = 1, max_iter

         call helmeos

         do i = 1, nrow

            ! if this point is converged, go to the next one
            if (NR_converged(i)) cycle

            ! energy difference
            delta_S(i) = Sgoal(i) - stot_row(i)

            ! keep things safe with bisect-limits
            if (delta_S(i).gt.0) then
               t_lower(i) = temp_row(i)
            else
               t_upper(i) = temp_row(i)
            end if

            ! update temperature
            delta_T(i) = delta_S(i) / dst_row(i)
            temp_row(i) = temp_row(i) + delta_T(i)

            ! if this took us out of bounds, don't let it happen
            ! choose a new point inside the interval [t_lower, t_upper]
            ! the point is in the middle of the interval (logarthmically)

            if ((temp_row(i).gt.t_upper(i)).OR.(temp_row(i).lt.t_lower(i))) then
               temp_row(i) = sqrt(t_lower(i) * t_upper(i))
            end if

            ! calculate relative errors
            rerr_S(i) = delta_S(i) / Sgoal(i)
            rerr_T(i) = delta_T(i) / temp_row(i)

            ! if we're at tolerances, end this
            if ((abs(rerr_S(i)).LE.rtol).AND.(abs(rerr_T(i)).LE.rtol)) then
               NR_converged(i) = .TRUE.
            endif

            !allow points at the temperature floor to "converge"
            if (t_upper(i).le.temp_floor * (1d0 + rtol)) then
               NR_converged(i) = .TRUE.
               temp_row(i) = temp_floor
            end if

         end do

         if (ALL(NR_converged)) exit

      end do

      ! once more, with feeling
      NR_converged = .FALSE.

      call helmeos

      end subroutine call_helmeos_DS


      subroutine call_helmeos_DE(nrow, den, ener, abar, zbar, tguess)
      include 'implno.dek'
      include 'vector_eos.dek'

! tests the eos routine
!
! ionmax  = number of isotopes in the network
! xmass   = mass fraction of isotope i
! aion    = number of nucleons in isotope i
! zion    = number of protons in isotope i

      integer, intent(in) :: nrow
      double precision, intent(in), dimension(nrow) :: den, ener, abar, zbar, tguess
      !f2py INTEGER, INTENT(hide) :: nrow
      !f2py DOUBLE PRECISION, DIMENSION(nrow), INTENT(in) :: den, ener, abar, zbar, tguess
      double precision, dimension(nrow) :: rerr_e, rerr_T
      double precision, dimension(nrow) :: delta_e, delta_T
      double precision, dimension(nrow) :: T_lower, T_upper
      double precision, dimension(nrow) :: egoal
      logical, dimension(nrow) :: NR_converged

      double precision, parameter :: temp_floor = 1e4
      double precision, parameter :: rtol = 1e-6

      integer :: i, iter
      integer, parameter :: max_iter  = 100

! don't try and overfill the array
!      if (ninput.gt.nrowmax) then stop

! read the data table
      call read_helm_table

! set the input vector. pipeline is only 1 element long

      jlo_eos = 1 ; jhi_eos = nrow
      do i = 1, nrow
         egoal(i) = ener(i) / den(i) ! eos works on specific internal energy
         temp_row(i) = tguess(i)
         den_row(i)  = den(i)
         abar_row(i) = abar(i)
         zbar_row(i) = zbar(i)
      end do

      T_lower = temp_floor
      T_upper = 1e12
      NR_converged = .FALSE.

      ! do the NR iteration

      do iter = 1, max_iter

         call helmeos

         do i = 1, nrow

            ! if this point is converged, go to the next one
            if (NR_converged(i)) cycle

            ! energy difference
            delta_E(i) = egoal(i) - etot_row(i)

            ! keep things safe with bisect-limits
            if (delta_E(i).gt.0) then
               t_lower(i) = temp_row(i)
            else
               t_upper(i) = temp_row(i)
            end if

            ! update temperature
            delta_T(i) = delta_E(i) / det_row(i)
            temp_row(i) = temp_row(i) + delta_T(i)

            ! if this took us out of bounds, don't let it happen
            ! choose a new point inside the interval [t_lower, t_upper]
            ! the point is in the middle of the interval (logarthmically)

            if ((temp_row(i).gt.t_upper(i)).OR.(temp_row(i).lt.t_lower(i))) then
               temp_row(i) = sqrt(t_lower(i) * t_upper(i))
            end if

            ! calculate relative errors
            rerr_e(i) = delta_e(i) / egoal(i)
            rerr_T(i) = delta_T(i) / temp_row(i)

            ! if we're at tolerances, end this
            if ((abs(rerr_e(i)).LE.rtol).AND.(abs(rerr_T(i)).LE.rtol)) then
               NR_converged(i) = .TRUE.
            endif

            !allow points at the temperature floor to "converge"
            if (t_upper(i).le.temp_floor * (1d0 + rtol)) then
               NR_converged(i) = .TRUE.
               temp_row(i) = temp_floor
            end if

         end do

         if (ALL(NR_converged)) exit

      end do

      ! once more, with feeling
      NR_converged = .FALSE.

      call helmeos

    end subroutine call_helmeos_DE
