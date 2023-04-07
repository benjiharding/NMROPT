module test_suite1

   use testdrive, only: new_unittest, unittest_type, error_type, check
   use test_variograms, only: get_variogram_pairs, get_variogram_values
   use test_de_optimization, only: minimize_ackley, minimize_beale

   implicit none

   private

   public :: collect_suite1

   real, parameter :: EPSLON = 1e-6

contains

!> Collect all exported unit tests
   subroutine collect_suite1(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("test_vario_pairs", test_vario_pairs), &
                  new_unittest("test_vario_values", test_vario_values), &
                  new_unittest("test_ackley_function", test_ackley_function), &
                  new_unittest("test_beale_function", test_beale_function) &
                  ]

   end subroutine collect_suite1

   subroutine test_vario_pairs(error)
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: pairs(:, :)
      real(8), allocatable :: lagbins(:)
      real(8) :: diff
      integer :: true(2, 2)
      integer :: i

      true(1, :) = [1, 5]
      true(2, :) = [4, 5]

      call get_variogram_pairs(pairs, lagbins)

      diff = abs(sqrt(2.0)*0.5 - lagbins(2))

      do i = 1, size(true, dim=1)
         call check(error, pairs(i, 1), true(i, 1))
         call check(error, pairs(i, 2), true(i, 2))
         call check(error, pairs(i, 3), 2)
         call check(error, (diff .lt. EPSLON), .true.)
         if (allocated(error)) return
      end do

   end subroutine test_vario_pairs

   subroutine test_vario_values(error)
      type(error_type), allocatable, intent(out) :: error
      real(8), allocatable :: expvario(:)
      real(8) :: true, diff

      true = ((1.0 - 1.0)**2 + (0.5 - 1.0)**2)/2*0.5

      call get_variogram_values(expvario)

      diff = abs(true - expvario(1))

      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

   end subroutine test_vario_values

   subroutine test_ackley_function(error)
      type(error_type), allocatable, intent(out) :: error
      real(8) :: true(2), best(2), diff(2)
      integer :: i

      true = [0.d0, 0.d0]
      call minimize_ackley(best)
      diff = true - best

      do i = 1, size(true)
         call check(error, (diff(i) .lt. EPSLON), .true.)
         if (allocated(error)) return
      end do

   end subroutine test_ackley_function

   subroutine test_beale_function(error)
      type(error_type), allocatable, intent(out) :: error
      real(8) :: true(2), best(2), diff(2)
      integer :: i

      true = [3.d0, 0.d50]
      call minimize_beale(best)
      diff = true - best

      do i = 1, size(true)
         call check(error, (diff(i) .lt. EPSLON), .true.)
         if (allocated(error)) return
      end do

   end subroutine test_beale_function

end module test_suite1
