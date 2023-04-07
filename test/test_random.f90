program test_random

   use randiter_mod

   implicit none

   integer, parameter :: rseed = 175686

   call init_genrand(rseed)
   call randiter

contains

end program test_random
