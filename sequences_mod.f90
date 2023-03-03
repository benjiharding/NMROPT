module sequences_mod

   use constants

   implicit none

   integer :: iruns, inpoint

contains

   ! TODO: n-point connectivity

   subroutine calculate_runs(arr, maxrun, cumbins)

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

   end subroutine calculate_runs

end module sequences_mod
