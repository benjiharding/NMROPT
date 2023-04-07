module vario_mod

   ! use geostat
   use rotationmatrix, only: set_rmat
   use covasubs
   use subs
   use constants

   implicit none

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

   subroutine update_vario(head, tail, zval, expvario, sill)

      ! recalculate experimental variogram using updated zval

      ! parameters
      type(lags), intent(in) :: head, tail
      real(8), intent(in) :: zval(:), sill

      ! result
      real(8), allocatable, intent(out) :: expvario(:)

      ! local variables
      real(8) :: expv
      integer :: i, np, nl

      nl = size(head%lags)
      allocate (expvario(nl)) ! move this allocation outside subroutine?
      expvario = -999.d0

      do i = 1, nl
         expv = 0.d0
         np = size(head%lags(i)%idxs)
         expv = sum((zval(tail%lags(i)%idxs) - zval(head%lags(i)%idxs))**2)
         expvario(i) = expv/dble(np)
      end do

      expvario = 0.5d0*expvario
      expvario = expvario/sill

   end subroutine update_vario

   subroutine calc_expsill(var, sill)

      ! sill of traditional variogram (variance)

      real(8), intent(in) :: var(:)
      real(8), intent(out) :: sill
      real(8) :: mean, sumsqs
      integer :: i, nd

      nd = size(var)
      mean = 0.d0
      sumsqs = 0.d0

      do i = 1, nd
         mean = mean + var(i)
         sumsqs = sumsqs + var(i)*var(i)
      end do

      mean = mean/nd
      sumsqs = sumsqs/nd
      sill = sumsqs - mean*mean

   end subroutine calc_expsill

   subroutine vario_mse(expvario, varmodelvals, varlagdist, idwpow, mse)

      ! Weighted MSE between experimental points and variogram model

      ! parameters
      real(8), intent(inout) :: expvario(:)
      real(8), intent(in) :: varmodelvals(:), varlagdist(:)
      real(8), intent(in) :: idwpow

      ! result
      real(8), intent(out) :: mse

      ! local variables
      real(8), allocatable :: idw_wts(:)
      integer :: nlags

      ! initilaize a few variables
      nlags = size(varlagdist)
      mse = 0.d0

      ! get weights by lag distance
      idw_wts = inv_dist(varlagdist, idwpow, nlags)

      ! get the weighted MSE
      mse = sum(idw_wts*(varmodelvals - expvario)**2)/dble(nlags)
      ! mse = 1.d0/dble(nlags)*mse/sum(idw_wts)

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

   ! subroutine vario_pairs(x, y, z, azm, atol, bandh, dip, dtol, &
   !                        bandv, nlags, lagdis, lagtol, nd, &
   !                        pairs, lagbins)

   !    ! data indices and lag index for valid variogram pairs

   !    ! parameters
   !    real(8), intent(in) :: x(:), y(:), z(:)
   !    real(8), intent(in) :: azm, atol, bandh, dip, dtol, &
   !                           bandv, lagdis, lagtol
   !    integer, intent(in) :: nlags, nd

   !    ! result
   !    integer, allocatable, intent(out) :: pairs(:, :)
   !    real(8), allocatable, intent(out) :: lagbins(:)

   !    ! local variables
   !    real(8), allocatable :: lagh(:)
   !    real(8) :: azmuth, uvxazm, uvyazm, csatol, declin, &
   !               uvzdec, uvhdec, csdtol, band
   !    real(8) :: maxdis, dx, dy, dz, dxs, dys, dzs, dxy, &
   !               h, hs, dcazm, dcdec
   !    integer :: i, j, k, n, np, maxpairs

   !    maxdis = (nlags*lagdis)**2

   !    ! specify angles and distances
   !    azmuth = (90.d0 - azm)*PI/180.d0
   !    uvxazm = cos(azmuth)
   !    uvyazm = sin(azmuth)

   !    if (atol .le. 0.d0) then
   !       csatol = cos(45.d0*PI/180.d0)
   !    else
   !       csatol = cos(atol*PI/180.d0)
   !    end if

   !    declin = (90.d0 - dip)*PI/180.d0
   !    uvzdec = cos(declin)
   !    uvhdec = sin(declin)

   !    if (dtol .le. 0.d0) then
   !       csdtol = cos(45.d0*PI/180.d0)
   !    else
   !       csdtol = cos(dtol*PI/180.d0)
   !    end if

   !    ! initialize max possible pairs and counter
   !    k = 0
   !    maxpairs = 0
   !    do i = 1, nd
   !       do j = i, nd
   !          maxpairs = maxpairs + 1
   !       end do
   !    end do

   !    ! allocate the output arrays
   !    allocate (pairs(maxpairs, 3)) ! head idx, tail idx, lag idx
   !    allocate (lagh(maxpairs), lagbins(nlags))

   !    ! iterate over all pairs
   !    do i = 1, nd
   !       do j = i, nd

   !          dx = x(j) - x(i)
   !          dy = y(j) - y(i)
   !          dz = z(j) - z(i)
   !          dxs = dx*dx
   !          dys = dy*dy
   !          dzs = dz*dz
   !          hs = dxs + dys + dzs

   !          if (hs .lt. EPSLON) cycle
   !          if (hs .gt. maxdis) cycle

   !          h = sqrt(max(hs, 0.d0))

   !          ! check the azimuth is acceptable
   !          dxy = sqrt(max((dxs + dys), 0.d0))
   !          if (dxy .lt. EPSLON) then
   !             dcazm = 1.0
   !          else
   !             dcazm = (dx*uvxazm + dy*uvyazm)/dxy
   !          end if
   !          if (abs(dcazm) .lt. csatol) cycle

   !          ! check horizontal bandwidth
   !          band = uvxazm*dy - uvyazm*dx
   !          if (abs(band) .gt. bandh) cycle

   !          ! check dip angle
   !          if (dcazm .lt. 0.d0) then
   !             dxy = -dxy
   !          end if
   !          if (dxy .lt. EPSLON) then
   !             dcdec = 0.d0
   !          else
   !             dcdec = (dxy*uvhdec + dz*uvzdec)/h
   !             if (abs(dcdec) .lt. csdtol) cycle
   !          end if

   !          ! check vertical bandwidth
   !          band = uvhdec*dz - uvzdec*dxy
   !          if (abs(band) .gt. bandv) cycle

   !          ! finaly check the lag tolerance
   !          ! iterate over all bins as they may overlap
   !          do n = 0, nlags - 1

   !             if (h .ge. (n*lagdis - lagtol)) then
   !                if (h .le. (n*lagdis + lagtol)) then

   !                   ! pair is acceptable if we got this far
   !                   k = k + 1
   !                   pairs(k, 1) = i
   !                   pairs(k, 2) = j
   !                   pairs(k, 3) = n + 1
   !                   lagh(k) = h

   !                end if
   !             end if
   !          end do

   !          ! end main loop
   !       end do
   !    end do

   !    ! only retain the populated pairs
   !    pairs = pairs(1:k, :)
   !    np = size(pairs, dim=1)

   !    ! get average distance in each lag bin
   !    lagbins = HUGE(h)
   !    do n = 1, nlags
   !       h = 0.d0
   !       k = 0
   !       do i = 1, np
   !          if (pairs(i, 3) .eq. n) then
   !             h = h + lagh(i)
   !             k = k + 1
   !          end if
   !       end do
   !       lagbins(n) = h/k
   !    end do

   ! end subroutine vario_pairs

   subroutine vario_pairs(xyz, azm, azmtol, bandhorz, dip, diptol, &
                          bandvert, tilt, nlags, lagdist, lagtol, ndata, &
                          pairs, lagbins)

      ! This subroutine assembles variogram pair indices following the
      ! logic of varcalc.f90

      ! parameters
      real(8), intent(in) :: xyz(3, ndata)
      real(8), intent(in) :: azm, azmtol, bandhorz, dip, diptol, &
                             bandvert, lagdist, lagtol, tilt
      integer, intent(in) :: nlags, ndata

      ! result
      integer, allocatable, intent(out) :: pairs(:, :)
      real(8), allocatable, intent(out) :: lagbins(:)

      ! local variables
      real(8), allocatable :: lagh(:)
      real(8) :: xyzrot(3, ndata), dblerotindex(ndata)
      integer :: introtindex(ndata)
      real(8) :: angle, azmtolrad, diptolrad, diptoldist, azmtoldist, farthest, &
                 ymin, ymax, zmin, zmax, xdist, h
      real(8), dimension(3, 3) :: forwardrotmat, reverserotmat
      integer :: i, j, k, n, np, maxpairs, startidx
      logical :: omni

      ! Initialize index tracking for values
      do i = 1, ndata
         dblerotindex(i) = i
      end do

      ! Correct the angle tolerances to lie between 0 and 90 degrees
      omni = .false.
      if (azmtol .ge. 90d0) then
         azmtolrad = 90d0*DEG2RAD
         omni = .true.
      elseif (azmtol .lt. 0d0) then
         azmtolrad = 0d0
      else
         azmtolrad = azmtol*DEG2RAD
      end if

      if (diptol .ge. 90d0) then
         diptolrad = 90d0*DEG2RAD
         omni = .true.
      elseif (diptol .lt. 0d0) then
         diptolrad = 0d0
      else
         diptolrad = diptol*DEG2RAD
      end if

      ! Precalculate the point where the angle tolerances no longer come in to play
      if (azmtol .gt. 90d0 - SMALLDBLE) then
         azmtoldist = 0d0
      elseif (azmtol .lt. SMALLDBLE) then
         azmtoldist = BIGDBLE
      else
         azmtoldist = bandhorz/tan(azmtolrad)
      end if

      if (diptol .gt. 90d0 - SMALLDBLE) then
         diptoldist = 0d0
      elseif (diptol .lt. SMALLDBLE) then
         diptoldist = BIGDBLE
      else
         diptoldist = bandvert/tan(diptolrad)
      end if

      ! Set up the rotation matrix using the GSLIB angles
      call setrotmat(azm, dip, tilt, forwardrotmat, reverserotmat)

      ! Rotate to align the new x' axis with the variogram direction
      do i = 1, ndata
         xyzrot(:, i) = matmul(xyz(:, i), forwardrotmat)
      end do

      ! ! Sort along the direction vector
      ! call dblemodsortem(xyzrot(1, :), ndata, 3, xyzrot(2, :), xyzrot(3, :), &
      !                    dblerotindex)

      ! Track the sorted indices for computing values
      do i = 1, ndata
         introtindex(i) = int(dblerotindex(i))
      end do

      ! initialize max possible pairs and counter
      k = 0
      maxpairs = 0
      do i = 1, ndata
         do j = i, ndata
            if (i .eq. j) cycle
            maxpairs = maxpairs + 1
         end do
      end do

      ! allocate the output arrays
      allocate (pairs(maxpairs, 3)) ! head idx, tail idx, lag idx
      allocate (lagh(maxpairs), lagbins(nlags + 1))

      do i = 1, ndata - 1   ! This point, i, is the "tail" variable in GSLIB notation

         ! Farthest point to consider from the current location
         farthest = xyzrot(1, i) + nlags*lagdist + lagtol
         ! Pre-calculated bandwidth min/max y's and z's
         ymax = xyzrot(2, i) + bandhorz
         ymin = xyzrot(2, i) - bandhorz
         zmax = xyzrot(3, i) + bandvert
         zmin = xyzrot(3, i) - bandvert
         ! Loop over all possible pairs forward of the current point
         if (omni) then
            startidx = 1
         else
            startidx = i + 1
         end if

         do j = startidx, ndata   ! This point, j, is the "head" variable in GSLIB notation

            ! If we are too far away, skip checking other pairs
            if (xyzrot(1, j) .gt. farthest) cycle !exit
            ! Check horizontal bandwidths
            if (xyzrot(2, j) .gt. ymax) cycle
            if (xyzrot(2, j) .lt. ymin) cycle
            ! Check vertical bandwidths
            if (xyzrot(3, j) .gt. zmax) cycle
            if (xyzrot(3, j) .lt. zmin) cycle
            ! Check if the points are at the same point on the x' axis
            xdist = abs(xyzrot(1, j) - xyzrot(1, i))
            if (i .eq. j) then
               cycle
            else
               ! Check azimuth tolerance if it matters at this distance
               if (xdist .lt. azmtoldist) then
                  ! Check absolute azimuth tolerance
                  angle = abs(atan((xyzrot(2, j) - xyzrot(2, i))/(max(xdist, SMALLDBLE))))
                  if (angle .gt. azmtolrad) cycle
               end if
               ! Check dip tolerance if it matters at this distance
               if (xdist .lt. diptoldist) then
                  ! Check absolute dip tolerance
                  angle = abs(atan((xyzrot(3, j) - xyzrot(3, i))/(max(xdist, SMALLDBLE))))
                  if (angle .gt. diptolrad) cycle
               end if
            end if

            ! At this point all tests have been passed so place the pair in the
            ! appropriate bin - due to possible overlapping bins each bin
            ! should be checked to see if pairs fall in it
            ! Find the distance between the values
            h = distance(xyzrot(:, i), xyzrot(:, j))

            do n = 0, nlags

               ! Are we inside this bin
               if ((h .ge. (n*lagdist - lagtol)) .and. &
                   (h .le. (n*lagdist + lagtol))) then

                  ! pair is acceptable if we got this far
                  k = k + 1
                  pairs(k, 1) = i
                  pairs(k, 2) = j
                  pairs(k, 3) = n + 1
                  lagh(k) = h

               end if

            end do
         end do
      end do

      ! only retain the populated pairs
      pairs = pairs(1:k, :)
      np = size(pairs, dim=1)

      ! get average distance in each lag bin
      lagbins = HUGE(h)
      do n = 1, nlags + 1
         h = 0.d0
         k = 0
         do i = 1, np
            if (pairs(i, 3) .eq. n) then
               h = h + lagh(i)
               k = k + 1
            end if
         end do
         lagbins(n) = h/k
      end do

   end subroutine vario_pairs

   function inv_dist(lagdis, power, nlags) result(idw_wts)

      ! inverse distance power weighting

      real(8), intent(in) :: lagdis(nlags)
      real(8), intent(in) :: power
      integer, intent(in) :: nlags
      real(8) :: idw_wts(nlags)

      idw_wts = 1.d0/lagdis**power
      ! idw_wts = idw_wts/maxval(idw_wts)
      idw_wts = idw_wts/sum(idw_wts)

   end function inv_dist

   subroutine set_sill(vm)

      implicit none

      type(variogram), intent(inout) :: vm(:)
      integer :: n, i, j

      n = ubound(vm(:), 1)

      do i = 1, n
         vm(i)%sill = vm(i)%c0
         do j = 1, vm(i)%nst
            if (vm(i)%it(j) .eq. 4) then
               vm(i)%sill = vm(i)%sill + 999.D+00
            else
               vm(i)%sill = vm(i)%sill + vm(i)%cc(j)
            end if
         end do
      end do

   end subroutine set_sill

   subroutine set_rotmatrix(vm)

      implicit none

      type(variogram), intent(inout) :: vm(:)

      integer :: n, i, j, test

      n = ubound(vm(:), 1)

      do i = 1, n
         if (.not. (allocated(vm(i)%rm))) then
            allocate (vm(i)%rm(3, 3, vm(i)%nst), stat=test)
            if (test .ne. 0) stop 'ERROR: Allocation failed due to insufficient memory.'
         end if
         do j = 1, vm(i)%nst
            vm(i)%rm(:, :, j) = set_rmat([vm(i)%ang1(j), vm(i)%ang2(j), vm(i)%ang3(j)], &
                                         [1.0D+00, vm(i)%anis1(j), vm(i)%anis2(j)])
         end do
      end do

   end subroutine set_rotmatrix

end module vario_mod
