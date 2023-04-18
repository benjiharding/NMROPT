module test_suite4

   use testdrive, only: new_unittest, unittest_type, error_type, check
   use test_network_forward, only: reshape_trial_vector_to_matrices, &
                                   calculate_forward_pass
   use types_mod, only: network
   use network_mod, only: init_network, sigmoid

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
                  new_unittest("test_bias_vector_to_matrices", test_bias_vector_to_matrices), &
                  new_unittest("test_network_forward_pass_linear", test_network_forward_pass_linear), &
                  new_unittest("test_network_forward_pass_sigmoid", test_network_forward_pass_sigmoid) &
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

   subroutine test_network_forward_pass_linear(error)

      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: dims = 15, af = 4, ndata = 5, nfact = 2
      real(8) :: diff, true(ndata), test(ndata)
      real(8) :: vector(dims), X(ndata, nfact)
      type(network) :: net
      integer :: i

      true = [28., 36., 44., 52., 60.]
      vector = [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 0., 0., 0., 0., 0.]
      X(:, 1) = [1., 2., 3., 4., 5.]
      X(:, 2) = [6., 7., 8., 9., 10.]

      call calculate_forward_pass(vector, X, af, net, test)

      do i = 1, ndata
         diff = abs(test(i) - true(i))
         call check(error, (diff .lt. EPSLON), .true.)
         if (allocated(error)) return
      end do

   end subroutine test_network_forward_pass_linear

   subroutine test_network_forward_pass_sigmoid(error)

      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: nl = 2, dims = 9, af = 1, ndata = 5, nfact = 2
      real(8) :: diff, true(ndata), test(ndata)
      real(8) :: vector(dims), X(ndata, nfact)
      type(network) :: net
      integer :: i

      true = [1., 1., 1., 1., 1.]
      vector = [1., 1., 1., 1., 1., 1., 0., 0., 0.]
      X(:, 1) = [0.1, 0.2, 0.3, 0.4, 0.5]
      X(:, 2) = [-0.1, -0.2, -0.3, -0.4, -0.5]
      call calculate_forward_pass(vector, X, af, net, test)

      do i = 1, ndata
         diff = abs(test(i) - true(i))
         call check(error, (diff .lt. EPSLON), .true.)
         if (allocated(error)) return
      end do

   end subroutine test_network_forward_pass_sigmoid

end module test_suite4
