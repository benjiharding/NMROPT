module test_suite2

   use testdrive, only: new_unittest, unittest_type, error_type, check
   use test_de_optimization, only: minimize_ackley, minimize_beale
   use test_objective, only: calc_arbitrary_objective
   use mtmod

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
                  !new_unittest("test_objective_calculation", test_objective_calculation) &
                  ]

   end subroutine collect_suite2

   subroutine test_ackley_function(error)
      type(error_type), allocatable, intent(out) :: error
      real(8) :: true(2), diff(2)
      real(8), allocatable :: best(:)
      integer :: i

      true = [0.d0, 0.d0]
      call minimize_ackley(best)
      diff = abs(true - best)

      do i = 1, size(true)
         call check(error, (diff(i) .lt. 1e-3), .true.)
         if (allocated(error)) return
      end do

   end subroutine test_ackley_function

   subroutine test_beale_function(error)
      type(error_type), allocatable, intent(out) :: error
      real(8) :: true(2), diff(2)
      real(8), allocatable :: best(:)
      integer :: i

      true = [3.d0, 0.5d0]
      call minimize_beale(best)
      diff = abs(true - best)

      do i = 1, size(true)
         call check(error, (diff(i) .lt. 1e-3), .true.)
         if (allocated(error)) return
      end do

   end subroutine test_beale_function

   subroutine test_objective_calculation(error)
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: dims = 13
      real(8) :: opt_vector(dims), rand_vector(dims)
      real(8) :: opt_obj, rand_obj
      integer :: i

      opt_vector = [0.14647875, 0.14907116, 0.54336461, 0.07923219, 0.52928929, 0.44513836, &
                    0.989115, 0.49769444, 0.6068023, 0.97930614, 0.0, 0.0, 0.0]
      rand_vector = 0.0
      do i = 1, 10
         rand_vector(i) = grnd()
      end do

      call calc_arbitrary_objective(opt_vector, rand_vector, opt_obj, rand_obj)

      print *, "optimal (?) objective value: ", opt_obj
      print *, "random objective value: ", rand_obj

   end subroutine test_objective_calculation

end module test_suite2
