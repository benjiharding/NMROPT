program main

   ! Network model of regionalization optimzation. Weight
   ! and bias values are optimized using differential evolution
   ! for the user defined network architecture and specified
   ! objectives.
   !
   ! Author: Ben Harding
   ! Date: May 2023
   ! Location: Centre for Computational Geostatistics,
   ! University of Alberta, Edmonton, Canada
   ! Contact: bharding@ualberta.ca

   use readpar_mod
   use mtmod
   use sim_mod, only: simulate
   use network_mod, only: init_network
   use objective_mod, only: init_objective
   use de_mod, only: optimize
   use output_mod, only: write_files

   implicit none

   real, parameter :: VERSION = 1.000

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
   call optimize(ipara)

   ! write output files
   call write_files

   ! close all output files
   close (lout)
   close (ldbg)
   close (lwts)
   close (lobj)

   ! finished
   write (*, *) " "
   write (*, "(A,f5.3,A)") "NMROPT version ", VERSION, " finished"

end program main
