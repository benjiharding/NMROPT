module test_suite3

   use testdrive, only: new_unittest, unittest_type, error_type, check
   use test_kriging_matrices, only: solve_kriging_matrices
   use test_sgs, only: test_sequential_simulation

   implicit none

   private

   public :: collect_suite3

   real, parameter :: EPSLON = 1e-6

contains

!> Collect all exported unit tests
   subroutine collect_suite3(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("test_righthand_side_vector", test_righthand_side_vector), &
                  new_unittest("test_kriging_weight_vector", test_kriging_weight_vector), &
                  ! new_unittest("test_lefthand_side_matrix", test_lefthand_side_matrix), &
                  new_unittest("test_conditional_mean", test_conditional_mean), &
                  new_unittest("test_conditional_stdev", test_conditional_stdev), &
                  new_unittest("test_sequential", test_sequential) &
                  ]

   end subroutine collect_suite3

   subroutine test_righthand_side_vector(error)

      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ndata = 5
      real(8) :: true(ndata), diff(ndata)
      real(8) :: rhs(ndata), kwts(ndata), lhs(ndata, ndata)
      real(8) :: cmean, cstdev
      integer :: i

      true = [0.71927831, 0.80977724, 0.37184356, 0.51305866, 0.76507662]

      call solve_kriging_matrices(rhs, lhs, kwts, cmean, cstdev)
      diff = abs(true - rhs)

      do i = 1, ndata
         call check(error, (diff(i) .lt. EPSLON), .true.)
         if (allocated(error)) return
      end do

   end subroutine test_righthand_side_vector

   subroutine test_kriging_weight_vector(error)

      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ndata = 5
      real(8) :: true(ndata), diff(ndata)
      real(8) :: rhs(ndata), kwts(ndata), lhs(ndata, ndata)
      real(8) :: cmean, cstdev
      integer :: i

      true = [0.23523521, 0.4931524, -0.00926474, 0.08906618, 0.23397578]

      call solve_kriging_matrices(rhs, lhs, kwts, cmean, cstdev)
      diff = abs(true - kwts)

      do i = 1, ndata
         call check(error, (diff(i) .lt. EPSLON), .true.)
         if (allocated(error)) return
      end do

   end subroutine test_kriging_weight_vector

   subroutine test_lefthand_side_matrix(error)

      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ndata = 5
      real(8) :: true(ndata, ndata), diff
      real(8) :: rhs(ndata), kwts(ndata), lhs(ndata, ndata)
      real(8) :: cmean, cstdev
      integer :: i, j

      true(1, :) = [1., 0.54173791, 0.51131047, 0.66032732, 0.69583306]
      true(2, :) = [0.54173791, 1., 0.24976248, 0.36994101, 0.67765035]
      true(3, :) = [0.51131047, 0.24976248, 1., 0.83085248, 0.27207199]
      true(4, :) = [0.66032732, 0.36994101, 0.83085248, 1., 0.40141153]
      true(5, :) = [0.69583306, 0.67765035, 0.27207199, 0.40141153, 1.]

      call solve_kriging_matrices(rhs, lhs, kwts, cmean, cstdev)

      do i = 1, ndata
         do j = 1, ndata
            diff = abs(true(i, j) - lhs(i, j))
            call check(error, (diff .lt. EPSLON), .true.)
            if (allocated(error)) return
         end do
      end do

   end subroutine test_lefthand_side_matrix

   subroutine test_conditional_mean(error)

      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ndata = 5
      real(8) :: true, diff
      real(8) :: rhs(ndata), kwts(ndata), lhs(ndata, ndata)
      real(8) :: cmean, cstdev

      true = 2.59692242

      call solve_kriging_matrices(rhs, lhs, kwts, cmean, cstdev)
      diff = abs(true - cmean)
      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

   end subroutine test_conditional_mean

   subroutine test_conditional_stdev(error)

      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ndata = 5
      real(8) :: true, diff
      real(8) :: rhs(ndata), kwts(ndata), lhs(ndata, ndata)
      real(8) :: cmean, cstdev

      true = 0.4584716875634445

      call solve_kriging_matrices(rhs, lhs, kwts, cmean, cstdev)
      diff = abs(true - cstdev)
      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

   end subroutine test_conditional_stdev

   subroutine test_sequential(error)

      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: ndata = 3
      real(8) :: true(ndata), test(ndata), diff
      integer :: i

      true = [-0.53980719, 0.05544304, -0.89299043]

      call test_sequential_simulation(test)

      do i = 1, ndata
         diff = abs(true(i) - test(i))
         call check(error, (diff .lt. 1e-3), .true.)
         if (allocated(error)) return
      end do

   end subroutine test_sequential

end module test_suite3
