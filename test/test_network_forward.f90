module test_network_forward

   use types_mod
   use network_mod
   use subs

   implicit none

contains

   subroutine reshape_trial_vector_to_matrices(net)

      integer, parameter :: nl = 3, dims = 11, af = 4
      type(network), intent(inout) :: net
      real(8) :: vector(dims)

      net%nl = nl
      net%dims = dims
      net%af = af
      net%ld = [3, 2, 1]

      call init_network(net)

      vector = [11., 12., 13., 21., 22., 23., 34., 35., 91., 92., 93.]

      call vector_to_matrices(vector, net)

   end subroutine reshape_trial_vector_to_matrices

   subroutine calculate_forward_pass(vector, X, af, net, AL)

      integer, parameter :: nl = 4, dims = 9
      type(network), intent(inout) :: net
      real(8), intent(in) :: vector(:), X(:, :)
      integer, intent(in) :: af
      real(8), intent(out) :: AL(:)

      net%nl = nl
      net%dims = dims
      net%af = af
      net%ld = [2, 2, 2, 1]

      call init_network(net)
      call vector_to_matrices(vector, net)
      ! call build_refcdf(nsamp, yref, nnet, ttable)
      call network_forward(net, X, AL, .false.)

   end subroutine calculate_forward_pass

end module test_network_forward
