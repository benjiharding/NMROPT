module test_de_optimization

   use de_mod, only: de, objfunc, ackley, beale

   implicit none

contains

   subroutine optimize_with_de(ifunc, best)

      real(8), parameter :: mut = 0.7, cplo = 0.5, cphi = 1.0
      real(8), parameter :: bmin = -5, bmax = 5
      integer, parameter :: popsize = 50, its = 10000, dims = 2
      integer, intent(in) :: ifunc
      real(8), allocatable, intent(inout) :: best(:)

      call de(dims, popsize, its, mut, cplo, cphi, bmin, bmax, best, ifunc)

   end subroutine optimize_with_de

   subroutine minimize_ackley(best)

      integer, parameter :: ifunc = 1
      real(8), allocatable, intent(inout) :: best(:)

      call optimize_with_de(ifunc, best)

   end subroutine minimize_ackley

   subroutine minimize_beale(best)

      integer, parameter :: ifunc = 2
      real(8), allocatable, intent(inout) :: best(:)
      call optimize_with_de(ifunc, best)

   end subroutine minimize_beale

end module test_de_optimization
