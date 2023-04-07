module test_de_optimization

   use de_mod, only: de, objfunc
   use mtmod

   implicit none

contains

   subroutine minimize_ackley(result)

      real(8), intent(out) :: result(2)
      real(8) :: vect(2)
      real(8), parameter :: mut = 0.7, cplo = 0.5, cphi = 1.0!, crossp = 0.8
      real(8), parameter :: bmin = -5, bmax = 5
      integer, parameter :: popsize = 10, its = 1000

      vect(1) = grnd()
      vect(2) = grnd()

      call de(1)

   end subroutine minimize_ackley

end module test_de_optimization
