program calc_runs
   implicit none

   integer, allocatable :: arr(:), diff(:)
   integer, allocatable :: start_idxs(:), runs(:), cumruns(:), cumbins(:)
   integer :: min, maxrun
   ! local variables
   integer, allocatable :: indices(:), indices_pad(:)
   integer :: i, j, ii, idx, pad, max_cum_runs

   min = 0
   pad = 2
   arr = [1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, &
          0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, &
          0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, &
          1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, &
          0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0]

   maxrun = 30

   ! test neighbours for similarity
   diff = arr(2:) - arr(1:size(arr))
   allocate (start_idxs(size(diff)))
   start_idxs(:) = 0

   do i = 1, size(diff)
      if (diff(i) .ne. 0) then
         start_idxs(i) = i
      end if
   end do

   ! count counts the number of non-zero elements in the array
   idx = count(start_idxs > min)
   allocate (indices(idx))

   ! get start indices of different runs
   j = 1
   do i = 1, size(start_idxs)
      if (start_idxs(i) > min) then
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
   runs = indices_pad(2:) - indices_pad(1:idx + pad)

   ! get max number of cumulative runs
   max_cum_runs = 0
   do i = 1, size(runs)
      do j = 1, runs(i)
         max_cum_runs = max_cum_runs + j
      end do
   end do

   allocate (cumruns(max_cum_runs))

   ! calculate cumulative runs
   ii = 1
   do i = 1, size(runs)
      do j = 1, runs(i)
         cumruns(ii:ii + j - 1) = runs(i) - j + 1
         ii = ii + j
      end do
   end do

   ! get the bincount of each run
   allocate (cumbins(maxrun))
   cumbins = 0
   do i = 1, size(cumruns)
      j = cumruns(i)
      cumbins(j) = cumbins(j) + 1
   end do

   print *, cumbins

end program calc_runs
