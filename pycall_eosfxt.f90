    subroutine call_eosfxt(nrow, den, temp, abar, zbar)
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
      call eosfxt

      end
