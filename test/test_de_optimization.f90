module test_de_optimization

   use de_mod, only: de, objfunc, ackley, beale
   use mtmod

   implicit none

contains

   subroutine optimize_with_de(ifunc, best)

      real(8), parameter :: mut = 0.7, cplo = 0.5, cphi = 1.0
      real(8), parameter :: bmin = -5, bmax = 5
      integer, parameter :: popsize = 10, its = 1000
      integer, intent(in) :: ifunc
      real(8), intent(inout) :: best(:)
      real(8) :: vect(2)

      best = -999.d0

      vect(1) = grnd()
      vect(2) = grnd()

      call de(ifunc=ifunc)

   end subroutine optimize_with_de

   subroutine minimize_ackley(best)
      integer, parameter :: ifunc = 1
      real(8), intent(inout) :: best(2)
      call optimize_with_de(ifunc, best)
   end subroutine minimize_ackley

   subroutine minimize_beale(best)
      integer, parameter :: ifunc = 2
      real(8), intent(inout) :: best(2)
      call optimize_with_de(ifunc, best)
   end subroutine minimize_beale

end module test_de_optimization
