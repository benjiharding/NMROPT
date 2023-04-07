module test_suite2

   use testdrive, only: new_unittest, unittest_type, error_type, check
   use test_de_optimization, only: minimize_ackley, minimize_beale

   implicit none

   private

   public :: collect_suite2

   real, parameter :: EPSLON = 1e-6

contains

!> Collect all exported unit tests
   subroutine collect_suite2(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("test_ackley_function", test_ackley_function), &
                  new_unittest("test_beale_function", test_beale_function) &
                  ]

   end subroutine collect_suite2

   subroutine test_ackley_function(error)
      type(error_type), allocatable, intent(out) :: error
      real(8) :: true(2), diff(2)
      real(8), allocatable :: best(:)
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
      real(8) :: true(2), diff(2)
      real(8), allocatable :: best(:)
      integer :: i

      true = [3.d0, 0.d50]
      call minimize_beale(best)
      diff = true - best

      do i = 1, size(true)
         call check(error, (diff(i) .lt. EPSLON), .true.)
         if (allocated(error)) return
      end do

   end subroutine test_beale_function

end module test_suite2
