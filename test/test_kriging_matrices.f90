module test_kriging_matrices

   use types_mod
   use kdtree2_module
   use covasubs
   use vario_mod, only: set_sill, set_rotmatrix
   use sim_mod, only: krige

   implicit none

contains

   subroutine solve_kriging_matrices(rhs, lhs, kwts, cmean, cstdev)

      ! data realted
      integer, parameter :: ndata = 5
      real(8) :: xyz(3, ndata + 1)
      real(8) :: var(ndata + 1)

      ! grid location
      integer, parameter :: nloc = 1

      ! variogram
      integer, parameter :: nst = 1, nvario = 1
      type(variogram) :: vm(nvario)

      ! kdtree required variables
      type(kdtree2), pointer :: tree
      type(kdtree2_result), allocatable :: results(:)
      integer :: nfound

      ! normal equations
      real(8), intent(inout) :: rhs(ndata), lhs(ndata, ndata), kwts(ndata)
      real(8), intent(inout) :: cmean, cstdev
      integer :: nuse(ndata), useidx(ndata, nloc)

      ! simulation
      integer :: isim(nloc), simidx

      ! indexes
      integer :: i, j

      real(8), parameter :: EPSLON = 1e-6, MINCOV = 1e-3

      !
      ! assign some random values
      !
      ! simloc is the first location
      xyz(1, :) = [0.50, 0.58933001, 0.38304157, 1.57397092, 1.27504179, 0.07812583]
      xyz(2, :) = [0.50, 1.06117351, 0.13580072, 1.31266704, 1.15120579, 0.71562721]
      xyz(3, :) = [0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0]
      var(:) = [0.0, 1.47332501, 2.65293378, 0.95760393, 0.33950179, 3.9349273]

      ! spherical variogram model with a isotropic range of 3
      vm(1)%nst = nst
      vm(1)%c0 = 0.0
      allocate (vm(1)%it(vm(1)%nst), vm(1)%cc(vm(1)%nst), &
                vm(1)%ang1(vm(1)%nst), vm(1)%ang2(vm(1)%nst), vm(1)%ang3(vm(1)%nst), &
                vm(1)%aa(vm(1)%nst), vm(1)%anis1(vm(1)%nst), vm(1)%anis2(vm(1)%nst), &
                vm(1)%ahmin(vm(1)%nst), vm(1)%avert(vm(1)%nst))
      vm(1)%cc(1) = 1.0
      vm(1)%ang1(1) = 0.0
      vm(1)%ang2(1) = 0.0
      vm(1)%ang3(1) = 0.0
      vm(1)%aa(1) = 3.0
      vm(1)%ahmin(1) = 3.0
      vm(1)%avert(1) = 3.0
      vm(1)%it(1) = 1
      vm(1)%anis1(1) = vm(1)%ahmin(1)/max(vm(1)%aa(1), EPSLON)
      vm(1)%anis2(1) = vm(1)%avert(1)/max(vm(1)%aa(1), EPSLON)
      call set_sill(vm)
      call set_rotmatrix(vm)

      ! define the search tree
      allocate (results(ndata + 1))
      tree => kdtree2_create(input_data=xyz, dim=3, sort=.false., rearrange=.true.)

      ! all location are initially unsimulated
      isim = 0

      ! loop over locations to simulate
      do i = 1, nloc

         ! what location are we at?
         simidx = 1

         ! query tree for this location
         call kdtree2_r_nearest(tp=tree, qv=xyz(:, simidx), r2=vm(1)%aa(1)*vm(1)%aa(1), &
                                nfound=nfound, nalloc=ndata + 1, results=results)

         nuse(i) = 0
         do j = 1, nfound

            ! check if this data index is the simulation index
            if (results(j)%idx .eq. simidx) cycle

            ! ! check if this data index is simulated or not
            ! if (isim(results(j)%idx) .eq. 0) cycle ! no value here

            ! meet minimum covariance? (not collocated)
            rhs(1) = get_cov(vm(1), xyz(:, simidx), xyz(:, results(j)%idx))
            if (rhs(1) .lt. MINCOV) cycle

            ! if we got this far increment number found at location i
            nuse(i) = nuse(i) + 1
            ! track data indices found at location i
            useidx(nuse(i), i) = results(j)%idx
         end do

         ! calc matrices for normal equations
         if (nuse(i) .gt. 0) then
            call krige(i, vm(1), xyz, rhs, lhs, kwts, nuse, useidx, var, &
                       simidx, cmean, cstdev)

         else
            ! if no data distribution is N(0,1)
            cmean = 0.d0
            cstdev = 1.d0

         end if

      end do

   end subroutine solve_kriging_matrices

end module test_kriging_matrices
