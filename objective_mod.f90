module objective_mod

   use geostat
   use sequences_mod, only: binary_runs, npoint_connect
   use vario_mod, only: update_vario, vario_mse, indicator_transform, &
                        vario_pairs, varmodelpts, set_sill
   use network_mod, only: network_forward
   use mtmod
   use subs
   use constants

   implicit none

contains

   subroutine init_objective()

      ! initialize objective targets, values and scaling components

      integer, allocatable :: tmppairs(:, :), tmpruns(:), tmparr(:)
      integer, allocatable :: headt(:), tailt(:), lagt(:)
      integer, allocatable :: lag_idxs(:), sub_idxs(:), valid_idxs(:)
      real(8), allocatable :: tmpbins(:), tmpnpt(:)

      integer :: maxlags, minlag, maxlag
      integer :: i, ii, j, k, nl

      ! data coordinates
      x = xyz(1, :)
      y = xyz(2, :)
      z = xyz(3, :)

      ! allocate data arrays
      allocate (iz(ndata, ncut), ivars(ncut))
      allocate (AL(ndata), AL_i(ndata, ncut))

      ! indicator transform of 'var'
      call indicator_transform(var, thresholds, ndata, ncut, iz, ivars)

      ! get the start/end dh indices in the mixture array
      allocate (udhidx(ndh + 1))
      udhidx = cumsum(dhlens)

      ! get variogram lag and pair data
      allocate (heads%dirs(ndir), tails%dirs(ndir))
      allocate (varazm%dirs(ndir), vardip%dirs(ndir), varlagdist%dirs(ndir))

      do i = 1, ndir
         call vario_pairs(x, y, z, expvar(i)%azm, expvar(i)%atol, expvar(i)%bandh, &
                          expvar(i)%dip, expvar(i)%dtol, expvar(i)%bandv, &
                          expvar(i)%nlags, expvar(i)%lagdis, expvar(i)%lagtol, &
                          ndata, tmppairs, tmpbins)
         headt = tmppairs(:, 1)
         tailt = tmppairs(:, 2)
         lagt = tmppairs(:, 3)

         ! ! this assumes defined lags are consecutive...
         ! ! should probably check for all valid lag indices and put them in a separate  array
         ! minlag = minval(lagt)
         ! maxlag = maxval(lagt)
         ! nl = maxlag - minlag + 1

         ! check for valid lag indices (defined lag bins)
         valid_idxs = pack([(ii, ii=1, expvar(i)%nlags)], tmpbins(:) .lt. HUGE(tmpbins))
         nl = size(valid_idxs)

         ! only allocate the number of defined lags
         allocate (varazm%dirs(i)%vlags(nl))
         allocate (vardip%dirs(i)%vlags(nl))
         allocate (varlagdist%dirs(i)%vlags(nl))
         allocate (heads%dirs(i)%lags(nl))
         allocate (tails%dirs(i)%lags(nl))

         do j = 1, nl

            ! experiemntal parameters
            varazm%dirs(i)%vlags(j) = expvar(i)%azm
            vardip%dirs(i)%vlags(j) = expvar(i)%dip
            varlagdist%dirs(i)%vlags(j) = tmpbins(valid_idxs(j))

            ! lag indices
            lag_idxs = pack([(k, k=1, size(tmppairs, dim=1))], &
                            tmppairs(:, 3) .eq. valid_idxs(j))
            call get_subsample(lag_idxs, max_pairs, sub_idxs)
            allocate (heads%dirs(i)%lags(j)%idxs(size(sub_idxs)))
            allocate (tails%dirs(i)%lags(j)%idxs(size(sub_idxs)))
            do k = 1, size(sub_idxs)
               heads%dirs(i)%lags(j)%idxs(k) = headt(sub_idxs(k))
               tails%dirs(i)%lags(j)%idxs(k) = tailt(sub_idxs(k))
            end do

         end do
      end do

      ! write out pairs if debugging
      if (idbg .gt. 0) then
         open (lprs, file="vario_pairs.out", status="UNKNOWN")
         write (lprs, "(A)") "Experimental Variogram Pairs"
         write (lprs, "(i1)") 4
         write (lprs, "(A)") "head idx"
         write (lprs, "(A)") "tail idx"
         write (lprs, "(A)") "lag idx"
         write (lprs, "(A)") "dir idx"
         do i = 1, ndir
            do j = 1, size(heads%dirs(i)%lags)
               do k = 1, size(heads%dirs(i)%lags(j)%idxs)
                  write (lprs, "(*(i7, 1x))") heads%dirs(i)%lags(j)%idxs(k), &
                     tails%dirs(i)%lags(j)%idxs(k), j, i
               end do
            end do
         end do
         close (lprs)
      end if

      ! variogram target
      allocate (target_vario%dirs(ndir))
      do i = 1, ndir
         maxlags = size(varlagdist%dirs(i)%vlags)
         allocate (target_vario%dirs(i)%vlags(maxlags))
         call varmodelpts(vmod(1)%nst, vmod(1)%c0, vmod(1)%it, vmod(1)%cc, &
                          vmod(1)%ang1, vmod(1)%ang2, vmod(1)%ang3, vmod(1)%aa, &
                          vmod(1)%ahmin, vmod(1)%avert, maxlags, &
                          varlagdist%dirs(i)%vlags, varazm%dirs(i)%vlags, &
                          vardip%dirs(i)%vlags, target_vario%dirs(i)%vlags)
      end do

      ! indicator variogram target
      allocate (target_ivario%cuts(ncut))

      do j = 1, ncut

         allocate (target_ivario%cuts(j)%dirs(ndir))

         ! rescale sill parameters
         ivmod(j)%c0 = ivmod(j)%c0*ivars(j)
         do i = 1, ivmod(j)%nst
            ivmod(j)%cc(i) = ivmod(j)%cc(i)*ivars(j)
         end do

         ! calculate the model points
         do k = 1, ndir
            maxlags = size(varlagdist%dirs(k)%vlags)
            allocate (target_ivario%cuts(j)%dirs(k)%vlags(maxlags))
            call varmodelpts(ivmod(j)%nst, ivmod(j)%c0, ivmod(j)%it, &
                             ivmod(j)%cc, ivmod(j)%ang1, ivmod(j)%ang2, &
                             ivmod(j)%ang3, ivmod(j)%aa, ivmod(j)%ahmin, &
                             ivmod(j)%avert, maxlags, varlagdist%dirs(k)%vlags, &
                             varazm%dirs(k)%vlags, vardip%dirs(k)%vlags, &
                             target_ivario%cuts(j)%dirs(k)%vlags)
         end do
      end do
      call set_sill(ivmod)

      ! cumulative runs target
      allocate (target_runs(maxrun, ncut))
      target_runs(:, :) = 0
      do i = 1, ncut
         do j = 1, ndh
            tmparr = iz(udhidx(j) + 1:udhidx(j + 1), i)
            call binary_runs(tmparr, maxrun, tmpruns)
            target_runs(:, i) = target_runs(:, i) + tmpruns
         end do
      end do

      ! npoint connectivity target
      allocate (target_npoint(nstep, ncut))
      target_npoint(:, :) = 0.d0
      do i = 1, ncut
         call npoint_connect(iz(:, i), nstep, ndh, udhidx, tmpnpt)
         target_npoint(:, i) = tmpnpt
      end do

      ! scale the components
      write (*, *) " "
      write (*, *) "Scaling objective components..."

      call obj_scale

      ! write (*, "(*(i2,a2,g14.8,1x))") (i, ":", objscale(i), i=1, 4)

      write (*, "(A*(g14.8,1x))") "variogram: ", objscale%vario(1)
      write (*, "(A*(g14.8,1x))") "indicator variogram: ", objscale%ivario(1)
      write (*, "(A*(g14.8,1x))") "runs: ", (objscale%runs(i), i=1, ncut)
      write (*, "(A*(g14.8,1x))") "npoint: ", (objscale%npoint(i), i=1, ncut)

   end subroutine init_objective

   subroutine obj_scale()

      ! scale objective components by looping over MAXPERT iterations and
      ! tracking the how much each component changes

      integer, parameter :: MAXPERT = 1000
      real(8), allocatable :: vect_denorm(:), min_b(:), max_b(:), diff(:)
      real(8), allocatable :: trial(:), trial_denorm(:)
      ! real(8) :: objinit(4), objdelta(4), rescale
      type(objective) :: objinit, objdelta
      integer, allocatable :: iarr(:), expruns(:)
      integer :: cumruns(maxrun)
      real(8), allocatable :: phi_n(:)
      integer :: i, j, ic
      real(8) :: mse

      ! objscale = 1.d0
      ! objinit = 0.d0
      ! objdelta = 0.d0
      allocate (objinit%vario(1), objinit%ivario(1), &
                objinit%runs(ncut), objinit%npoint(ncut))
      allocate (objdelta%vario(1), objdelta%ivario(1), &
                objdelta%runs(ncut), objdelta%npoint(ncut))
      allocate (objscale%vario(1), objscale%ivario(1), &
                objscale%runs(ncut), objscale%npoint(ncut))

      objscale%vario = 1.d0
      objscale%ivario = 1.d0
      objscale%runs = 1.d0
      objscale%npoint = 1.d0

      objinit%vario = 0.d0
      objinit%ivario = 0.d0
      objinit%runs = 0.d0
      objinit%npoint = 0.d0

      objdelta%vario = 0.d0
      objdelta%ivario = 0.d0
      objdelta%runs = 0.d0
      objdelta%npoint = 0.d0

      ! initial trial vector
      allocate (min_b(size(vect)), max_b(size(vect)), trial(size(vect)))
      min_b = bmin
      max_b = bmax
      diff = abs(min_b - max_b)
      vect_denorm = min_b + vect*diff

      ! the choice of the first realization here is arbitrary
      call network_forward(ysimd(:, :, 1), vect_denorm, AL)
      call indicator_transform(AL, thresholds, ndata, ncut, AL_i, ivars)

      ! initilaize a starting value for each component
      if (vario .gt. 0) then
         call obj_vario(objt_vario)
         objinit%vario(1) = objt_vario
      end if
      if (ivario .gt. 0) then
         call obj_ivario(objt_ivario)
         objinit%ivario(1) = objt_ivario
      end if
      if (runs .gt. 0) then
         do i = 1, ncut
            mse = 0.d0
            cumruns = 0
            do j = 1, ndh
               iarr = AL_i(udhidx(j) + 1:udhidx(j + 1), i)
               call binary_runs(iarr, maxrun, expruns)
               cumruns = cumruns + expruns
            end do
            objinit%runs(i) = sum((dble(target_runs(:, i)) - dble(cumruns))**2)/dble(maxrun)
         end do
      end if
      if (npoint .gt. 0) then
         do i = 1, ncut
            mse = 0.d0
            call npoint_connect(AL_i(:, i), nstep, ndh, udhidx, phi_n)
            objinit%npoint(i) = sum((target_npoint(:, i) - phi_n)**2)/dble(nstep)
         end do
      end if

      ! iterate over the pertubations
      do i = 1, MAXPERT

         ! generate a random vector within the bounds
         do j = 1, size(vect)
            trial(j) = grnd()
         end do
         trial_denorm = min_b + trial*diff

         ! evalute the random vector
         call network_forward(ysimd(:, :, 1), trial_denorm, AL)
         call indicator_transform(AL, thresholds, ndata, ncut, AL_i, ivars)

         if (vario .gt. 0) then
            call obj_vario(objt_vario)
            if (objt_vario .lt. 0.0) objt_vario = objinit%vario(1)
            objdelta%vario(1) = objdelta%vario(1) + &
                                abs(objinit%vario(1) - objt_vario)
         end if

         if (ivario .gt. 0) then
            call obj_ivario(objt_ivario)
            if (objt_ivario .lt. 0.0) objt_ivario = objinit%ivario(1)
            objdelta%ivario(1) = objdelta%ivario(1) + &
                                 abs(objinit%ivario(1) - objt_ivario)
         end if

         if (runs .gt. 0) then

            do ic = 1, ncut
               mse = 0.d0
               cumruns = 0
               do j = 1, ndh
                  iarr = AL_i(udhidx(j) + 1:udhidx(j + 1), ic)
                  call binary_runs(iarr, maxrun, expruns)
                  cumruns = cumruns + expruns
               end do
               mse = sum((dble(target_runs(:, ic)) - dble(cumruns))**2)/dble(maxrun)
               if (mse .lt. 0.0) mse = objinit%runs(ic)
               objdelta%runs(ic) = objdelta%runs(ic) + abs(objinit%runs(ic) - mse)
            end do

         end if

         if (npoint .gt. 0) then

            do ic = 1, ncut
               mse = 0.d0
               call npoint_connect(AL_i(:, ic), nstep, ndh, udhidx, phi_n)
               mse = sum((target_npoint(:, ic) - phi_n)**2)/dble(nstep)
               if (mse .lt. 0.0) mse = objinit%npoint(ic)
               objdelta%npoint(ic) = objdelta%npoint(ic) + &
                                     abs(objinit%npoint(ic) - mse)
            end do

         end if

      end do

      ! scale objective components
      if (vario .gt. 0) objscale%vario(1) = MAXPERT/objdelta%vario(1)
      if (ivario .gt. 0) objscale%ivario(1) = MAXPERT/objdelta%ivario(1)
      do ic = 1, ncut
         if (runs .gt. 0) then
            objscale%runs(ic) = MAXPERT/objdelta%runs(ic)
         end if
         if (npoint .gt. 0) then
            objscale%npoint(ic) = MAXPERT/objdelta%npoint(ic)
         end if
      end do

      ! ! user defined scaling if required
      ! objscale(1) = userfac(1)*objscale(1)
      ! objscale(2) = userfac(2)*objscale(2)
      ! objscale(3) = userfac(3)*objscale(3)
      ! objscale(4) = userfac(4)*objscale(4)

   end subroutine obj_scale

   subroutine obj_nmr(v, gobjt)

      ! network model of regionalization objective
      ! returns scalar expected objective value

      real(8), intent(in) :: v(:) ! trial vector
      real(8), intent(out) :: gobjt ! global temp obj value
      integer :: ireal

      gobjt = 0.d0
      objt_vario = 0.d0
      objt_ivario = 0.d0
      objt_runs = 0.d0
      objt_npt = 0.d0

      do ireal = 1, nreals

         call network_forward(ysimd(:, :, ireal), v, AL)
         call indicator_transform(AL, thresholds, ndata, ncut, AL_i, ivars)

         if (vario .gt. 0) call obj_vario(objt_vario)
         if (ivario .gt. 0) call obj_ivario(objt_ivario)
         if (runs .gt. 0) call obj_runs(objt_runs)
         if (npoint .gt. 0) call obj_npoint(objt_npt)

         gobjt = gobjt + objt_vario + objt_ivario + objt_runs + objt_npt

      end do

      gobjt = gobjt/nreals

   end subroutine obj_nmr

   subroutine obj_vario(objt)

      ! continuous variogram component

      real(8), intent(inout) :: objt

      real(8), allocatable :: expvario(:)
      real(8) :: mse
      integer :: i

      objt = 0.d0

      do i = 1, ndir
         call update_vario(heads%dirs(i), tails%dirs(i), AL, expvario)
         call vario_mse(expvario, target_vario%dirs(i)%vlags, &
                        varlagdist%dirs(i)%vlags, dble(idwpow), mse)
         objt = objt + mse
      end do

      objt = objt*objscale%vario(1)

   end subroutine obj_vario

   subroutine obj_ivario(objt)

      ! indicator variogram component

      real(8), intent(inout) :: objt

      real(8), allocatable :: expivario(:)
      real(8) :: mse
      integer :: i, j

      objt = 0.d0

      do j = 1, ncut
         do i = 1, ndir
            call update_vario(heads%dirs(i), tails%dirs(i), dble(AL_i(:, j)), &
                              expivario)
            call vario_mse(expivario, target_ivario%cuts(j)%dirs(i)%vlags, &
                           varlagdist%dirs(i)%vlags, dble(idwpow), mse)
            objt = objt + mse
         end do
      end do

      objt = objt*objscale%ivario(1)

   end subroutine obj_ivario

   subroutine obj_runs(objt)

      ! cumulative run frequency component

      real(8), intent(inout) :: objt

      integer, allocatable :: iarr(:), expruns(:)
      integer :: cumruns(maxrun)
      real(8) :: mse
      integer :: i, j

      objt = 0.d0

      do i = 1, ncut
         mse = 0.d0
         cumruns = 0
         do j = 1, ndh
            iarr = AL_i(udhidx(j) + 1:udhidx(j + 1), i)
            call binary_runs(iarr, maxrun, expruns)
            cumruns = cumruns + expruns
         end do
         mse = sum((dble(target_runs(:, i)) - dble(cumruns))**2)/dble(maxrun)
         objt = objt + mse*objscale%runs(i)
      end do

      ! objt = objt/ncut
      ! objt = objt*objscale(3)

   end subroutine obj_runs

   subroutine obj_npoint(objt)

      ! npoint connectivity component

      real(8), intent(inout) :: objt

      real(8), allocatable :: phi_n(:)
      real(8) :: mse
      integer :: i

      objt = 0.d0

      do i = 1, ncut
         mse = 0.d0
         call npoint_connect(AL_i(:, i), nstep, ndh, udhidx, phi_n)
         mse = sum((target_npoint(:, i) - phi_n)**2)/dble(nstep)
         objt = objt + mse*objscale%npoint(i)
      end do

      ! objt = objt/ncut
      ! objt = objt*objscale(4)

   end subroutine obj_npoint

end module objective_mod
