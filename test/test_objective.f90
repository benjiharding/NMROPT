module test_objective

   use readpar_mod
   use sim_mod, only: simulate
   use network_mod, only: init_network
   use objective_mod, only: init_objective, obj_nmr
   use output_mod, only: write_files

   implicit none

contains

   subroutine calc_arbitrary_objective(vector1, vector2, obj1, obj2)

      real(8), intent(in) :: vector1(:), vector2(:)
      real(8), intent(out) :: obj1, obj2

      ! read the parfile
      call readpar

      ! initilize random generator
      call sgrnd(rseed)

      ! simulate at the data locations
      call simulate

      ! initialize network parameters
      call init_network(nnet)

      ! initialize objective function parameters
      call init_objective

      ! calc the objective value given intialization and input vector
      call obj_nmr(vector2, obj2)
      call obj_nmr(vector1, obj1)

      ! write output files
      allocate (best(size(vector1)))
      best = vector1
      call write_files

   end subroutine calc_arbitrary_objective

end module test_objective
