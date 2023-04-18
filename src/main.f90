program main

   use readpar_mod
   use sim_mod, only: simulate
   use network_mod, only: init_network
   use objective_mod, only: init_objective
   use de_mod, only: optimize
   use output_mod, only: write_files

   implicit none

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

   ! minimize the objective
   call optimize

   ! write output files
   call write_files

   ! close all output files
   close (lout)
   close (ldbg)
   close (lwts)
   close (lobj)

contains

end program main
