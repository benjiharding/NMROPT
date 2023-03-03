module constants

   implicit none

   real(8), parameter :: EPSLON = 1e-5
   real(8), parameter :: PI = 4*atan(1.0d0)
   real(8), parameter :: SMALLDBLE = 1d-6
   real(8), parameter :: DEG2RAD = PI/180d0
   integer, parameter :: MAXNST = 4
   integer, parameter :: MAXGNST = 1 ! max nst for Gaussian variograms

contains

end module constants
