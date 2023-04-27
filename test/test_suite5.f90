module test_suite5

   use testdrive, only: new_unittest, unittest_type, error_type, check
   use test_subroutines, only: test_nscore

   implicit none

   private

   public :: collect_suite5

   real, parameter :: EPSLON = 1e-6

contains

!> Collect all exported unit tests
   subroutine collect_suite5(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("test_unweighted_normal_score_transform", test_unweighted_normal_score_transform), &
                  new_unittest("test_weighted_normal_score_transform", test_weighted_normal_score_transform) &
                  ]

   end subroutine collect_suite5

   subroutine test_unweighted_normal_score_transform(error)

      type(error_type), allocatable, intent(out) :: error
      real(8) :: mu, tmu, sig, tsig, diff(2), gmin, gmax
      integer :: i

      tmu = 0.d0
      tsig = 1.d0

      call test_nscore(mu, sig, gmin, gmax, 0, 123456)

      diff(1) = abs(tmu - mu)
      diff(2) = abs(tsig - sig)

      do i = 1, 2
         call check(error, (diff(i) .lt. 1e-3), .true.)
         if (allocated(error)) return
      end do

      call check(error, (gmin .gt. -5), .true.)
      if (allocated(error)) return

      call check(error, (gmax .lt. 5), .true.)
      if (allocated(error)) return

   end subroutine test_unweighted_normal_score_transform

   subroutine test_weighted_normal_score_transform(error)

      type(error_type), allocatable, intent(out) :: error
      real(8) :: mu, tmu, sig, tsig, diff(2), gmin, gmax
      integer :: i

      tmu = 0.d0
      tsig = 1.d0

      call test_nscore(mu, sig, gmin, gmax, 1, 654321)

      diff(1) = abs(tmu - mu)
      diff(2) = abs(tsig - sig)

      do i = 1, 2
         call check(error, (diff(i) .lt. 1e-1), .true.)
         if (allocated(error)) return
      end do

      call check(error, (gmin .gt. -5), .true.)
      if (allocated(error)) return

      call check(error, (gmax .lt. 5), .true.)
      if (allocated(error)) return

   end subroutine test_weighted_normal_score_transform

end module test_suite5
