module sequences_mod

   use geostat
   use types_mod, only: lags, variogram
   use covasubs, only: get_cov
   use subs

   implicit none

contains

   ! TODO: n-point connectivity

   subroutine binary_runs(arr, maxrun, cumbins)

      ! Calculate runs in 1D binary sequence

      ! parameters
      integer, intent(in) :: arr(:)
      integer, intent(in) :: maxrun

      ! result
      integer, allocatable, intent(out) :: cumbins(:)

      ! local variables
      integer, allocatable :: start_idxs(:), runs(:), diff(:), cumruns(:)
      integer, allocatable :: indices(:), indices_pad(:), tmpruns(:)
      integer :: i, j, k, idx, min_, pad, maxcr, i0

      min_ = 0
      pad = 2
      i0 = arr(1) ! first either 1 or 0

      ! test neighbours for similarity
      diff = arr(2:) - arr(1:size(arr))
      allocate (start_idxs(size(diff)))
      start_idxs(:) = 0

      ! get the start index of each unique run
      do i = 1, size(diff)
         if (diff(i) .ne. 0) then
            start_idxs(i) = i
         end if
      end do

      ! count the number of non-zero elements in the array
      idx = count(start_idxs > min_)
      allocate (indices(idx))

      ! get start indices of different runs
      j = 1
      do i = 1, size(start_idxs)
         if (start_idxs(i) > min_) then
            indices(j) = i
            j = j + 1
         end if
      end do

      ! pad array to prevent index errors
      allocate (indices_pad(idx + pad))
      indices_pad(1) = 0
      indices_pad(idx + pad) = size(arr)
      indices_pad(2:idx + pad - 1) = indices

      ! calculate runs
      tmpruns = indices_pad(2:) - indices_pad(1:idx + pad)

      ! are we considering runs below the thresholds?
      if (iruns .eq. 0) then

         if (i0 .eq. 0) then ! first run is above
            runs = tmpruns(2 :: 2)
         else ! first run is below
            runs = tmpruns(1 :: 2)
         end if

         ! are we considering runs above the thresholds?
      else if (iruns .eq. 1) then

         if (i0 .eq. 0) then ! first run is above
            runs = tmpruns(1 :: 2)
         else ! first run is below
            runs = tmpruns(2 :: 2)
         end if

         ! are we considering runs above and below?
      else if (iruns .eq. 2) then
         runs = tmpruns

      end if

      ! get max number of cumulative runs
      maxcr = 0
      do i = 1, size(runs)
         do j = 1, runs(i)
            maxcr = maxcr + j
         end do
      end do

      ! allocate the output array
      allocate (cumruns(maxcr))
      cumruns = 0

      ! calculate cumulative runs
      k = 1
      do i = 1, size(runs)
         do j = 1, runs(i)
            cumruns(k:k + j - 1) = runs(i) - j + 1
            k = k + j
         end do
      end do

      ! get the bincount of each run length
      allocate (cumbins(maxrun))
      cumbins = 0
      do i = 1, size(cumruns)
         j = cumruns(i)
         if (j > maxrun) cycle
         cumbins(j) = cumbins(j) + 1
      end do

   end subroutine binary_runs

   subroutine npoint_connect(AL_i, nstep, ndh, udhidx, phi)

      ! global n-point connectivity function

      ! parameters
      integer, intent(in) :: AL_i(:), udhidx(:)
      integer, intent(in) :: nstep, ndh

      ! result
      real(8), allocatable, intent(out) :: phi(:)

      ! local variables
      integer, allocatable :: temp_idxs(:), idxs(:)
      integer, allocatable :: iarr(:), arr(:), lcount(:)
      real(8), allocatable :: prod(:)
      integer :: i, j, k, n, nx
      real(8) :: p

      allocate (phi(nstep), prod(nstep), lcount(nstep))
      prod = 0.d0
      lcount = 0

      ! loop over connected steps
      do n = 1, nstep

         allocate (temp_idxs(n + 1))

         temp_idxs(1) = 1
         do i = 1, n
            temp_idxs(i + 1) = i
         end do

         ! number of drillholes
         do k = 1, ndh

            iarr = AL_i(udhidx(k) + 1:udhidx(k + 1))
            nx = size(iarr)

            ! are we considering connectivity above the thresholds?
            if (inpoint .gt. 0) then
               iarr = 1 - iarr
            end if

            ! number of actual lags in vector given n and nx
            do i = 1, nx - n + 1

               idxs = temp_idxs + (i - 1)
               arr = iarr(idxs)

               ! product of indicators
               p = 1.d0
               do j = 1, size(arr)
                  p = p*arr(j)
               end do

               ! update products
               prod(n) = prod(n) + p

               ! increment lag counter
               lcount(n) = lcount(n) + 1

            end do
         end do

         deallocate (temp_idxs)

      end do

      ! phi = prod/dble(lcount)
      phi = prod/size(AL_i)

   end subroutine npoint_connect

   subroutine binary_runs3d()
   end subroutine binary_runs3d

   subroutine npoint_connect3d(head, tail, izval, nstep, vm, phi)

      ! N-point connectivity function in the tertiary direction
      ! specified by the experimental variogram search parameters.
      ! The lag vector h is assumed to be the mean composite length
      ! and tolerance parameters are taken from the search parameters.

      ! parameters
      type(lags), intent(in) :: head, tail
      type(variogram), intent(in) :: vm
      integer, intent(in) :: nstep
      integer, intent(in) :: izval(:)

      ! result
      real(8), allocatable, intent(out) :: phi(:)

      ! local variables
      integer :: iz(size(izval))
      real(8) :: exp, cov, a(3), b(3)
      real(8) :: std1, std2, mu
      integer :: i, j, n, np

      allocate (phi(nstep))

      ! indicators above threshold?
      if (inpoint .gt. 0) then
         iz = 1 - izval
      else
         iz = izval
      end if

      ! phi(1) is mean of all indicators
      exp = 0.d0
      j = 0
      do i = 1, size(iz)
         exp = exp + iz(i)
         j = j + 1
      end do
      phi(1) = exp/j

      do n = 2, nstep

         np = size(head%lags(n - 1)%idxs)
         a = [xyz(:, tail%lags(n - 1)%idxs(1))]
         b = [xyz(:, head%lags(n - 1)%idxs(1))]
         std1 = stdev_int(iz(tail%lags(n - 1)%idxs))
         std2 = stdev_int(iz(head%lags(n - 1)%idxs))
         cov = get_cov(vm, a, b)
         cov = cov*std1*std2

         mu = 0.d0
         do j = 1, np
            mu = mu + iz(tail%lags(n - 1)%idxs(j)) + iz(head%lags(n - 1)%idxs(j))
         end do
         mu = mu/(2.d0*dble(np))

         phi(n) = phi(n - 1)*mu + cov

      end do

   end subroutine npoint_connect3d

end module sequences_mod
