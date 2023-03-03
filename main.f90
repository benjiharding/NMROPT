program main

   use readpar_mod
   use lusim_mod
   use network_mod
   use objective_mod
   use de_mod

   implicit none

   real(8), allocatable :: opt_AL(:)
   integer :: i, j

   ! read the parfile
   call readpar

   ! initilize random generator
   call sgrnd(rseed)

   ! simulate at the data locations
   call simulate

   ! write out a realization if debugging
   if (idbg .gt. 0) then
      do i = 1, ndata
         write (ldbg, "(*(g14.8,1x))") (ysimd(i, j, dbgireal), j=1, ngvarg + 1)
      end do
   end if

   ! initialize network parameters
   call init_network

   ! initialize objective function parameters
   call init_objective

   ! minimize the objective
   call optimize

   ! write out the optimized network mixture
   allocate (opt_AL(ndata))
   do i = 1, nreals
      call network_forward(ysimd(:, :, i), best, opt_AL)
      do j = 1, ndata
         write (lout, "(*(g14.8,1x))") opt_AL(j)
      end do
   end do

   close (lout)
   close (ldbg)

end program main
