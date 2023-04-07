program main

   implicit none

   integer, parameter :: nstep = 5, xmax = 5, ndh = 1
   integer :: iarr(xmax)
   real(8) :: phi(nstep)
   integer, allocatable :: temp_idxs(:), idxs(:), arr(:), lcount(:)
   real(8), allocatable :: prod(:)
   integer :: i, j, k, n, nx
   real(8) :: p

   allocate (lcount(nstep), prod(nstep))
   iarr = [1, 1, 1, 0, 1]
   prod = 0.d0
   lcount = 0

   ! iarr = 1 - iarr

   ! connected steps
   do n = 1, nstep

      allocate (temp_idxs(n + 1))

      temp_idxs(1) = 1
      do i = 1, n
         temp_idxs(i + 1) = i
      end do

      ! number of drillholes
      do k = 1, ndh

         nx = size(iarr)

         ! number of actual lags given n and lenght of hole
         do i = 1, nx - n + 1

            idxs = temp_idxs + (i - 1)
            arr = iarr(idxs)

            ! product of indicators
            p = 1.d0
            do j = 1, size(arr)
               p = p*arr(j)
            end do

            print *, arr

            ! update products
            prod(n) = prod(n) + p

            ! increment lag counter
            lcount(n) = lcount(n) + 1

         end do

      end do

      deallocate (temp_idxs)

   end do

   phi = prod/dble(lcount)

   print *, phi

end program main
