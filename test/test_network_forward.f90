module test_network_forward

   use network_mod
   use subs

   implicit none

contains

   subroutine network_forward_pass(af, nstrans)

      ! parameters
      logical, intent(in) :: nstrans
      integer, intent(in) :: af

      ! network architecture
      integer, allocatable :: layer_dims(:), iwts(:), ibias(:)
      real(8), allocatable :: Ymat(:, :), v(:)
      real(8), allocatable :: AL(:)
      integer :: dims ! dimension of vect

      ! test example
      layer_dims = [3, 2, 1]
      !   v =

      call init_network

   end subroutine network_forward_pass

end module test_network_forward
