module vario_mod

   use covasubs
   use constants

   implicit none

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

contains

   subroutine indicator_transform(zval, zc, nd, ncut, iz, ivars)

      ! indicator transform of zval based on cutoffs zc

      ! parameters
      real(8), intent(in) :: zval(:), zc(:)
      integer, intent(in) :: nd, ncut

      ! result
      integer, intent(inout) :: iz(:, :)
      real(8), intent(inout) :: ivars(:)

      ! local variables
      integer :: i, j
      real(8) :: mu

      ! initalize indicators
      iz = 0

      ! iterate over cutoffs and set indicators
      do i = 1, nd
         do j = 1, ncut
            if (zval(i) .le. zc(j)) then
               iz(i, j) = 1
            end if
         end do
      end do

      ! indicator variance to scale sill
      do j = 1, ncut
         mu = 0.d0
         do i = 1, nd
            mu = mu + iz(i, j)
         end do
         mu = mu/nd
         ivars(j) = mu*(1 - mu)
      end do

   end subroutine indicator_transform

   ! subroutine update_vario(head, tail, lag, zval, nlags, expvario)

   !    ! recalculate experimental variogram using updated zval

   !    ! parameters
   !    integer, intent(in) :: head(:), tail(:), lag(:)
   !    real(8), intent(in) :: zval(:)
   !    integer, intent(in) :: nlags

   !    ! result
   !    real(8), allocatable, intent(out) :: expvario(:)

   !    ! local variables
   !    real(8) :: expv
   !    integer :: i, j, ap, np

   !    allocate (expvario(nlags))
   !    ap = size(head, dim=1)

   !    do i = 1, nlags
   !       expv = 0.d0
   !       np = 0
   !       do j = 1, ap ! dont need to iterate over all pairs on each lag...
   !          if (lag(j) .eq. i) then
   !             expv = expv + (zval(tail(j)) - zval(head(j)))**2
   !             np = np + 1
   !          else
   !             cycle
   !          end if
   !       end do
   !       expvario(i) = 1.d0/(2.d0*dble(np))*expv
   !    end do

   ! end subroutine update_vario

   subroutine update_vario(head, tail, zval, stride, expvario)

      ! recalculate experimental variogram using updated zval

      ! parameters
      type(lags), intent(in) :: head, tail
      real(8), intent(in) :: zval(:)
      integer, intent(in) :: stride(:)

      ! result
      real(8), allocatable, intent(out) :: expvario(:)

      ! local variables
      real(8) :: expv
      integer :: i, np, nl

      nl = size(head%lags)
      allocate (expvario(nl))

      do i = 1, nl
         expv = 0.d0
         np = size(head%lags(i)%idxs)
         if (np .eq. 0) cycle
         expv = sum((zval(head%lags(i)%idxs) - zval(tail%lags(i)%idxs))**2)
         expvario(i) = 1.d0/(2.d0*dble(np))*expv
      end do

   end subroutine update_vario

   subroutine vario_mse(expvario, varmodelvals, varlagdist, idwpow, mse)

      ! MSE between experimental points and variogram model

      ! parameters
      real(8), intent(inout) :: expvario(:)
      real(8), intent(in) :: varmodelvals(:), varlagdist(:)
      real(8), intent(in) :: idwpow

      ! result
      real(8), intent(out) :: mse

      ! local variables
      real(8), allocatable :: wts(:)
      integer :: nlags

      ! initilaize a few variables
      nlags = size(varlagdist)
      mse = 0.d0

      ! get weights by lag distance
      wts = inv_dist(varlagdist, idwpow, nlags)

      ! get the weighted MSE
      mse = sum(wts*(varmodelvals - expvario)**2)
      mse = 1.d0/dble(nlags)*mse/sum(wts)

   end subroutine vario_mse

   subroutine varmodelpts(nst, c0, it, cc, azm, dip, tilt, ahmax, ahmin, avert, &
                          nvargs, varlagdist, varazm, vardip, varmodelvals)
      !-----------------------------------------------------------------------
      ! This subroutine calculates variogram model points at input lags with
      ! provided lag distances, azimuths and dips using the GSLIB cova3
      ! subroutine.
      !
      ! Experimental variogram values:
      !  nvargs - number of experimental variogram points
      !  varlagdist(nvargs) - lag distances
      !  varvalue(nvargs) - variogram value
      !  varazm(nvargs) - variogram azimuth
      !  vardip(nvargs) - variogram dips
      !
      ! (C) Jared Deutsch 2014
      !-----------------------------------------------------------------------

      ! Parameters
      integer, intent(in) :: nst, nvargs
      real(kind=8), intent(in) :: c0
      integer, dimension(nst), intent(in) :: it
      real(kind=8), dimension(nst), intent(in) :: cc, azm, dip, tilt, ahmax, &
                                                  ahmin, avert
      real(kind=8), dimension(nvargs), intent(in) :: varlagdist, varazm, vardip

      ! result
      real(kind=8), dimension(nvargs), intent(out) :: varmodelvals

      ! Internal variables
      integer :: ist, ivarg

      ! GSLIB interface variables
      integer :: MAXROT
      real(kind=8), allocatable :: rotmat(:, :, :)
      real(kind=8) :: cmax, maxcov, cova, x, y, z, c0gslib
      real(kind=8), dimension(nst) :: ccgslib, aagslib, anis1, anis2

      ! GSLIB interface allocation
      MAXROT = nst
      allocate (rotmat(MAXROT, 3, 3))
      c0gslib = real(c0)

      ! GSLIB anisotropy ratios
      do ist = 1, nst
         anis1(ist) = real(ahmin(ist)/max(ahmax(ist), SMALLDBLE))
         anis2(ist) = real(avert(ist)/max(ahmax(ist), SMALLDBLE))
         ccgslib(ist) = real(cc(ist))
         aagslib(ist) = real(ahmax(ist))
      end do

      ! GSLIB rotation matrix
      rotmat = 0
      do ist = 1, nst
         call setrot(azm(ist), dip(ist), tilt(ist), &
                     anis1(ist), anis2(ist), ist, MAXROT, rotmat)
      end do

      ! Maximum covariance for calculating variogram values from cova3
      call cova3(0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 1, nst, nst, c0gslib, &
                 it, ccgslib(1:nst), aagslib, 1, MAXROT, rotmat, cmax, maxcov)

      ! Calculate variogram values
      do ivarg = 1, nvargs
         x = real(sin(DEG2RAD*varazm(ivarg))*cos(DEG2RAD*vardip(ivarg))*varlagdist(ivarg))
         y = real(cos(DEG2RAD*varazm(ivarg))*cos(DEG2RAD*vardip(ivarg))*varlagdist(ivarg))
         z = real(sin(DEG2RAD*vardip(ivarg))*varlagdist(ivarg))
         call cova3(0.d0, 0.d0, 0.d0, x, y, z, 1, nst, nst, c0gslib, &
                    it, ccgslib(1:nst), aagslib, 1, MAXROT, rotmat, cmax, cova)
         varmodelvals(ivarg) = dble(maxcov - cova)
      end do

   end subroutine varmodelpts

   subroutine vario_pairs(x, y, z, azm, atol, bandh, dip, dtol, &
                          bandv, nlags, lagdis, lagtol, nd, &
                          pairs, lagbins)

      ! data indices and lag index for valid variogram pairs

      ! parameters
      real(8), intent(in) :: x(:), y(:), z(:)
      real(8), intent(in) :: azm, atol, bandh, dip, dtol, &
                             bandv, lagdis, lagtol
      integer, intent(in) :: nlags, nd

      ! result
      integer, allocatable, intent(out) :: pairs(:, :)
      real(8), allocatable, intent(out) :: lagbins(:)

      ! local variables
      real(8), allocatable :: lagh(:)
      real(8) :: azmuth, uvxazm, uvyazm, csatol, declin, &
                 uvzdec, uvhdec, csdtol, band
      real(8) :: maxdis, dx, dy, dz, dxs, dys, dzs, dxy, &
                 h, hs, dcazm, dcdec
      integer :: i, j, k, n, np, maxpairs

      maxdis = (nlags*lagdis)**2

      ! specify angles and distances
      azmuth = (90.d0 - azm)*PI/180.d0
      uvxazm = cos(azmuth)
      uvyazm = sin(azmuth)

      if (atol .le. 0.d0) then
         csatol = cos(45.d0*PI/180.d0)
      else
         csatol = cos(atol*PI/180.d0)
      end if

      declin = (90.d0 - dip)*PI/180.d0
      uvzdec = cos(declin)
      uvhdec = sin(declin)

      if (dtol .le. 0.d0) then
         csdtol = cos(45.d0*PI/180.d0)
      else
         csdtol = cos(dtol*PI/180.d0)
      end if

      ! initialize max possible pairs and counter
      k = 0
      maxpairs = 0
      do i = 1, nd
         do j = i, nd
            maxpairs = maxpairs + 1
         end do
      end do

      ! allocate the output arrays
      allocate (pairs(maxpairs, 3)) ! head idx, tail idx, lag idx
      allocate (lagh(maxpairs), lagbins(nlags))

      ! iterate over all pairs
      do i = 1, nd
         do j = i, nd

            dx = x(j) - x(i)
            dy = y(j) - y(i)
            dz = z(j) - z(i)
            dxs = dx*dx
            dys = dy*dy
            dzs = dz*dz
            hs = dxs + dys + dzs

            if (hs .lt. EPSLON) cycle
            if (hs .gt. maxdis) cycle

            h = sqrt(max(hs, 0.d0))

            ! check the azimuth is acceptable
            dxy = sqrt(max((dxs + dys), 0.d0))
            if (dxy .lt. EPSLON) then
               dcazm = 1.0
            else
               dcazm = (dx*uvxazm + dy*uvyazm)/dxy
            end if
            if (abs(dcazm) .lt. csatol) cycle

            ! check horizontal bandwidth
            band = uvxazm*dy - uvyazm*dx
            if (abs(band) .gt. bandh) cycle

            ! check dip angle
            if (dcazm .lt. 0.d0) then
               dxy = -dxy
            end if
            if (dxy .lt. EPSLON) then
               dcdec = 0.d0
            else
               dcdec = (dxy*uvhdec + dz*uvzdec)/h
               if (abs(dcdec) .lt. csdtol) cycle
            end if

            ! check vertical bandwidth
            band = uvhdec*dz - uvzdec*dxy
            if (abs(band) .gt. bandv) cycle

            ! finaly check the lag tolerance
            ! iterate over all bins as they may overlap
            do n = 0, nlags - 1

               if (h .ge. (n*lagdis - lagtol)) then
                  if (h .le. (n*lagdis + lagtol)) then

                     ! pair is acceptable if we got this far
                     k = k + 1
                     pairs(k, 1) = i
                     pairs(k, 2) = j
                     pairs(k, 3) = n + 1
                     lagh(k) = h

                  end if
               end if
            end do

            ! end main loop
         end do
      end do

      ! only retain the populated pairs
      pairs = pairs(1:k, :)
      np = size(pairs, dim=1)

      ! get average distance in each lag bin
      do n = 1, nlags
         h = 0.d0
         k = 0
         do i = 1, np
            if (pairs(i, 3) .eq. n) then
               h = h + lagh(i)
               k = k + 1
            end if
            lagbins(n) = h/k
         end do
      end do

   end subroutine vario_pairs

   function inv_dist(lagdis, power, nlags) result(wts)

      ! inverse distance power weighting

      real(8), intent(in) :: lagdis(nlags)
      real(8), intent(in) :: power
      integer, intent(in) :: nlags
      real(8) :: wts(nlags)

      wts = 1.d0/lagdis**power
      ! wts = wts/maxval(wts)
      wts = wts/sum(wts)

   end function inv_dist

end module vario_mod
