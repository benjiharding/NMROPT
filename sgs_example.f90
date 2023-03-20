program main

   use mtmod
   use kdtree2_module
   use types_mod
   use vario_mod, only: set_sill, set_rotmatrix
   use covasubs
   use subs
   use constants

   implicit none

   ! data realted
   integer, parameter :: ndata = 5
   real(8) :: xyz(3, ndata)
   real(8) :: var(ndata), mu = 0.d0

   ! variogram
   integer, parameter :: nst = 1, nvario = 1
   type(variogram) :: vm(nvario)

   ! kdtree required variables
   type(kdtree2), pointer :: tree
   type(kdtree2_result), allocatable :: results(:)
   integer :: nfound

   ! normal equations
   integer, parameter :: nsearch = 2
   real(8) :: rhs(nsearch), lhs(nsearch, nsearch), kwts(nsearch)
   real(8) :: cmean, cstdev
   integer :: nuse(ndata), useidx(ndata, ndata)

   ! simulation
   real(8) :: sim(ndata)
   integer :: isim(ndata), randpath(ndata)
   real(8) :: p, xp
   integer :: simidx, ierr

   ! indexes
   integer :: i, j, k, test

   !
   ! assign some random values
   !
   xyz(1, :) = [0.58933001, 0.38304157, 1.57397092, 1.27504179, 0.07812583]
   xyz(2, :) = [1.06117351, 0.13580072, 1.31266704, 1.15120579, 0.71562721]
   xyz(3, :) = [0.d0, 0.d0, 0.d0, 0.d0, 0.d0]
   var(:) = [1.47332501, 2.65293378, 0.95760393, 0.33950179, 3.9349273]

   do i = 1, ndata
      mu = mu + var(i)
   end do
   mu = mu/ndata
   var = var - mu

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

   ! build the search tree
   allocate (results(ndata))
   tree => kdtree2_create(input_data=xyz, dim=3, sort=.true., rearrange=.true.)

   ! define random path
   randpath = [(i, i=1, ndata)]
   call shuffle(randpath)

   ! all location are initially unsimulated
   isim = 0
   useidx = 0
   nuse = 0

   ! loop over locations to simulate
   do i = 1, ndata

      ! what location are we at?
      simidx = randpath(i)

      ! query tree for this location
      call kdtree2_r_nearest(tp=tree, qv=xyz(:, simidx), r2=vm(1)%aa(1)*vm(1)%aa(1), &
                             nfound=nfound, nalloc=ndata, results=results)

      ! loop over samples found in search
      nuse(i) = 0
      do j = 1, nfound

         ! check if this data index is the simulation index
         if (results(j)%idx .eq. simidx) cycle

         ! check if this data index is simulated or not
         if (isim(results(j)%idx) .eq. 0) cycle ! no conditioning value here

         ! meet minimum covariance? (not collocated)
         rhs(1) = get_cov(vm(1), xyz(:, simidx), xyz(:, results(j)%idx))
         if (rhs(1) .lt. MINCOV) cycle

         ! if we got this far increment number found at ith location
         nuse(i) = nuse(i) + 1
         ! track data indices found at ith location
         useidx(nuse(i), i) = results(j)%idx

         ! have we met the max search?
         if (nuse(i) .ge. nsearch) exit
      end do

      ! calculate matrices for normal equations
      if (nuse(i) .gt. 0) then
         do j = 1, nuse(i)
            ! build rhs vector
            rhs(j) = get_cov(vm(1), xyz(:, simidx), xyz(:, useidx(j, i)))
            do k = j, nuse(i)
               ! diagonal
               if (j .eq. k) then
                  lhs(j, j) = 1.d0
               else
                  ! build lhs matrix
                  lhs(j, k) = get_cov(vm(1), xyz(:, useidx(k, i)), xyz(:, useidx(j, i)))
                  if (lhs(j, k) .lt. 0) stop 'ERROR: Negative covariance.'
                  lhs(k, j) = lhs(j, k)
               end if
            end do
         end do

         ! solve the kriging system - external call to LAPACK
         call solve(lhs(1:nuse(i), 1:nuse(i)), kwts(1:nuse(i)), rhs(1:nuse(i)), &
                    nuse(i), 1, test)

         ! calcualte conditional mean and stdev
         cmean = 0.d0
         cstdev = get_cov(vm(1), [0.d0, 0.d0, 0.d0], [0.d0, 0.d0, 0.d0]) ! variance
         do j = 1, nuse(i)
            cmean = cmean + kwts(j)*sim(useidx(j, i))
            cstdev = cstdev - kwts(j)*rhs(j)
         end do

      else
         ! if no data the distribution is N(0,1)
         cmean = 0.d0
         cstdev = 1.d0

      end if

      print *, cmean, cstdev

      ! draw a random number and simulate
      p = grnd()
      call gauinv(p, xp, ierr)
      sim(simidx) = xp*cstdev + cmean

      ! update this location as simulated
      isim(simidx) = 1

   end do

end program main
