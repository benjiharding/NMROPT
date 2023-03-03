module randiter_mod

   use random

   implicit none

contains

   subroutine randiter

      real(8) :: rand
      integer, parameter :: int = 10
      integer :: i, j

      do i = 1, 5
         write (*, *) "iteration: ", i
         do j = 1, 3
            rand = grnd()
            write (*, *) "rand ", j, rand, floor(rand*int + 1.d0)
         end do
      end do

   end subroutine randiter

end module randiter_mod
