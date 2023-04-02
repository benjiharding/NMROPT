module types_mod

   implicit none

   ! variogram model
   type variogram
      integer :: nst !number of nested structure
      real(8) :: c0 !nugget
      real(8) :: sill
      integer, allocatable :: it(:) !structure types
      real(8), allocatable :: cc(:) !structure contributions
      real(8), allocatable :: ang1(:), ang2(:), ang3(:) !structure angles
      real(8), allocatable :: aa(:), anis1(:), anis2(:) !structure anisotropy
      real(8), allocatable :: ahmin(:), avert(:) ! minor and vert ranges
      real(8), allocatable :: rm(:, :, :) ! rotation matrix
   end type variogram

   ! experimental search parameters
   type experimental
      real(8) :: azm, atol
      real(8) :: dip, dtol
      real(8) :: bandh, bandv
      real(8) :: tilt
      real(8) :: lagdis, lagtol
      integer :: nlags
   end type experimental

   ! objective values
   type objective
      real(8), allocatable :: vario(:)
      real(8), allocatable :: ivario(:)
      real(8), allocatable :: runs(:)
      real(8), allocatable :: npoint(:)
   end type objective

   ! ragged arrays for lag data
   type :: indices
      integer, allocatable :: idxs(:)
   end type indices
   type :: lags
      type(indices), allocatable :: lags(:)
   end type lags
   type :: lag_array
      type(lags), allocatable :: dirs(:)
   end type lag_array

   ! ragged arrays for variogram data
   type :: vlags
      real(8), allocatable :: vlags(:)
   end type vlags
   type :: vario_array
      type(vlags), allocatable :: dirs(:)
   end type vario_array
   type :: ivario_array
      type(vario_array), allocatable :: cuts(:)
   end type ivario_array

end module types_mod
