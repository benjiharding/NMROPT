module test_suite4

   use testdrive, only: new_unittest, unittest_type, error_type, check
   use test_network_forward, only: reshape_trial_vector_to_matrices
   use types_mod, only: network
   use network_mod, only: init_network

   implicit none

   private

   public :: collect_suite4

   real, parameter :: EPSLON = 1e-6

contains

!> Collect all exported unit tests
   subroutine collect_suite4(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  new_unittest("test_weight_matrix_dimensions", test_weight_matrix_dimensions), &
                  new_unittest("test_bias_matrix_dimensions", test_bias_matrix_dimensions), &
                  new_unittest("test_weight_vector_to_matrices", test_weight_vector_to_matrices), &
                  new_unittest("test_bias_vector_to_matrices", test_bias_vector_to_matrices) &
                  ]

   end subroutine collect_suite4

   subroutine test_weight_matrix_dimensions(error)

      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nl = 3, dims = 11, af = 4
      integer :: true(nl - 1, 2)
      real(8) :: diff
      type(network) :: net
      integer :: i, j

      true(1, :) = [2, 3]
      true(2, :) = [1, 2]

      call reshape_trial_vector_to_matrices(net)

      do i = 1, nl - 1
         do j = 1, 2
            diff = abs(true(i, j) - net%layer(i)%sw(j))
            call check(error, (diff .lt. EPSLON), .true.)
            if (allocated(error)) return
         end do
      end do

   end subroutine test_weight_matrix_dimensions

   subroutine test_bias_matrix_dimensions(error)

      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nl = 3, dims = 11, af = 4
      integer :: true(nl - 1, 2)
      real(8) :: diff
      type(network) :: net
      integer :: i, j

      true(1, :) = [2, 1]
      true(2, :) = [1, 1]

      call reshape_trial_vector_to_matrices(net)

      do i = 1, nl - 1
         do j = 1, 2
            diff = abs(true(i, j) - net%layer(i)%sb(j))
            call check(error, (diff .lt. EPSLON), .true.)
            if (allocated(error)) return
         end do
      end do

   end subroutine test_bias_matrix_dimensions

   subroutine test_weight_vector_to_matrices(error)

      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nl = 3, dims = 11, af = 4
      real(8) :: diff
      type(network) :: net, true
      integer :: i, j, k

      true%nl = nl
      true%dims = dims
      true%af = af
      true%ld = [3, 2, 1]
      call init_network(true)

      true%layer(1)%nnwts(1, :) = [11., 12., 13.]
      true%layer(1)%nnwts(2, :) = [21., 22., 23.]
      true%layer(1)%nnbias(:, 1) = [91., 92.]

      true%layer(2)%nnwts(1, :) = [34., 35.]
      true%layer(2)%nnbias(:, 1) = [93.]

      call reshape_trial_vector_to_matrices(net)

      ! check weight values
      do i = 1, nl - 1
         do j = 1, true%layer(i)%sw(1)
            do k = 1, true%layer(i)%sw(2)
               diff = abs(true%layer(i)%nnwts(j, k) - net%layer(i)%nnwts(j, k))
               call check(error, (diff .lt. EPSLON), .true.)
               if (allocated(error)) return
            end do
         end do
      end do

   end subroutine test_weight_vector_to_matrices

   subroutine test_bias_vector_to_matrices(error)

      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nl = 3, dims = 11, af = 4
      real(8) :: diff
      type(network) :: net, true
      integer :: i, j, k

      true%nl = nl
      true%dims = dims
      true%af = af
      true%ld = [3, 2, 1]
      call init_network(true)

      true%layer(1)%nnwts(1, :) = [11., 12., 13.]
      true%layer(1)%nnwts(2, :) = [21., 22., 23.]
      true%layer(1)%nnbias(:, 1) = [91., 92.]

      true%layer(2)%nnwts(1, :) = [34., 35.]
      true%layer(2)%nnbias(:, 1) = [93.]

      call reshape_trial_vector_to_matrices(net)

      ! check weight values
      do i = 1, nl - 1
         do j = 1, true%layer(i)%sb(1)
            do k = 1, true%layer(i)%sb(2)
               diff = abs(true%layer(i)%nnbias(j, k) - net%layer(i)%nnbias(j, k))
               call check(error, (diff .lt. EPSLON), .true.)
               if (allocated(error)) return
            end do
         end do
      end do

   end subroutine test_bias_vector_to_matrices

end module test_suite4
