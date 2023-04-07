module test_variograms

   use vario_mod
   use types_mod

   implicit none

contains

   subroutine get_variogram_pairs(pairs, lagbins)

      integer, parameter :: ndata = 5
      integer, allocatable, intent(out) :: pairs(:, :)
      real(8), allocatable, intent(out) :: lagbins(:)

      ! local
      real(8) :: locs(3, ndata) ! xyz(3, :)
      real(8) :: azm, azmtol, dip, diptol, tilt, bandhorz, bandvert
      integer :: nlags
      real(8) :: lagdist, lagtol

      locs(1, :) = [0.0, 1.0, 0.0, 1.0, 0.5]
      locs(2, :) = [0.0, 0.0, 1.0, 1.0, 0.5]
      locs(3, :) = [0.0, 0.0, 0.0, 0.0, 0.0]

      azm = 45.0
      azmtol = 5.0
      dip = 0.0
      diptol = 0.0
      tilt = 0.0
      bandhorz = 1.0
      bandvert = 1.0
      nlags = 1
      lagdist = sqrt(2.0)*0.5
      lagtol = 0.1

      call vario_pairs(locs, azm, azmtol, bandhorz, dip, diptol, &
                       bandvert, tilt, nlags, lagdist, lagtol, &
                       ndata, pairs, lagbins)

   end subroutine get_variogram_pairs

   subroutine get_variogram_values(expvario)

      integer, parameter :: ndata = 5
      real(8), allocatable, intent(inout) :: expvario(:)
      integer, allocatable :: pairs(:, :)
      type(lag_array) :: heads, tails
      real(8), allocatable :: lagbins(:)
      real(8) :: var(ndata)

      call get_variogram_pairs(pairs, lagbins)

      var = [1.0, 0.0, 0.0, 0.5, 1.0]
      allocate (heads%dirs(1), tails%dirs(1))
      allocate (heads%dirs(1)%lags(1), tails%dirs(1)%lags(1))

      heads%dirs(1)%lags(1)%idxs = pairs(:, 2)
      tails%dirs(1)%lags(1)%idxs = pairs(:, 1)

      call update_vario(heads%dirs(1), tails%dirs(1), var, expvario, 1.d0)

   end subroutine get_variogram_values

end module test_variograms
