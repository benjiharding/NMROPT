module objective_mod

   use geostat
   use sequences_mod, only: binary_runs, npoint_connect
   use vario_mod, only: update_vario, vario_mse, indicator_transform, &
                        vario_pairs, varmodelpts, set_sill, calc_expsill
   use network_mod, only: network_forward, vector_to_matrices, &
                          calc_regularization
   use mtmod, only: grnd
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
      integer :: maxlags
      integer :: i, ii, j, k, nl

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
         call vario_pairs(xyz, expvar(i)%azm, expvar(i)%atol, expvar(i)%bandh, &
                          expvar(i)%dip, expvar(i)%dtol, expvar(i)%bandv, expvar(i)%tilt, &
                          expvar(i)%nlags, expvar(i)%lagdis, expvar(i)%lagtol, &
                          ndata, tmppairs, tmpbins)
         tailt = tmppairs(:, 1)
         headt = tmppairs(:, 2)
         lagt = tmppairs(:, 3)

         ! check for valid lag indices (defined lag bins)
         valid_idxs = pack([(ii, ii=1, expvar(i)%nlags + 1)], tmpbins(:) .lt. HUGE(tmpbins))
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

      ! variogram target
      if (vario .gt. 0) then
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
      end if

      ! indicator variogram target
      if (ivario .gt. 0) then
         allocate (target_ivario%cuts(ncut))
         do j = 1, ncut

            allocate (target_ivario%cuts(j)%dirs(ndir))

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
      end if

      ! cumulative runs target
      if (runs .gt. 0) then
         if (t_iruns .eq. 0) then ! calculate from data, not file
            allocate (target_runs(maxrun, ncut))
            target_runs(:, :) = 0
            do i = 1, ncut
               do j = 1, ndh
                  tmparr = iz(udhidx(j) + 1:udhidx(j + 1), i)
                  call binary_runs(tmparr, maxrun, tmpruns)
                  target_runs(:, i) = target_runs(:, i) + tmpruns
               end do
            end do
         end if
      end if

      ! npoint connectivity target
      if (npoint .gt. 0) then
         if (t_inpoint .eq. 0) then ! calculate from data, not file
            allocate (target_npoint(nstep, ncut))
            target_npoint(:, :) = 0.d0
            do i = 1, ncut
               call npoint_connect(iz(:, i), nstep, ndh, udhidx, tmpnpt)
               target_npoint(:, i) = tmpnpt
            end do
         end if
      end if

      ! scale the components
      write (*, *) " "
      write (*, *) "Scaling objective components..."
      write (*, *) " "

      objscale = 1.d0
      call obj_scale()

      write (*, "(A*(g14.8,1x))") "Variogram component: ", objscale(1)
      write (*, "(A*(g14.8,1x))") "Indicator variogram component: ", objscale(2)
      write (*, "(A*(g14.8,1x))") "Runs component: ", objscale(3)
      write (*, "(A*(g14.8,1x))") "N-point component: ", objscale(4)

      ! write (*, "(A*(g14.8,1x))") "Variogram component: ", objscale%vario(1)
      ! write (*, "(A*(g14.8,1x))") "Indicator variogram component: ", objscale%ivario(1)
      ! write (*, "(A*(g14.8,1x))") "Runs component: ", (objscale%runs(i), i=1, ncut)
      ! write (*, "(A*(g14.8,1x))") "N-point component: ", (objscale%npoint(i), i=1, ncut)

   end subroutine init_objective

   subroutine obj_scale()

      ! scale objective components by looping over MAXPERT iterations and
      ! tracking the how much each component changes

      integer, parameter :: MAXPERT = 1000
      real(8), allocatable :: vect(:), vect_denorm(:), min_b(:), max_b(:), diff(:)
      real(8), allocatable :: trial(:), trial_denorm(:)
      real(8) :: objinit(4), objdelta(4)
      ! type(objective) :: objinit, objdelta
      integer :: i, j, ic

      real(8) :: mse
      integer, allocatable :: iarr(:), expruns(:)
      integer :: cumruns(maxrun)
      real(8), allocatable :: phi_n(:)

      objinit = 0.d0
      objdelta = 0.d0

      ! allocate (objinit%vario(1), objinit%ivario(ncut), &
      !           objinit%runs(ncut), objinit%npoint(ncut))
      ! allocate (objdelta%vario(1), objdelta%ivario(ncut), &
      !           objdelta%runs(ncut), objdelta%npoint(ncut))
      ! allocate (objscale%vario(1), objscale%ivario(ncut), &
      !           objscale%runs(ncut), objscale%npoint(ncut))

      ! objscale%vario = 1.d0
      ! objscale%ivario = 1.d0
      ! objscale%runs = 1.d0
      ! objscale%npoint = 1.d0

      ! objinit%vario = 0.d0
      ! objinit%ivario = 0.d0
      ! objinit%runs = 0.d0
      ! objinit%npoint = 0.d0

      ! objdelta%vario = 0.d0
      ! objdelta%ivario = 0.d0
      ! objdelta%runs = 0.d0
      ! objdelta%npoint = 0.d0

      ! initial trial vector
      allocate (vect(nnet%dims))
      allocate (min_b(nnet%dims), max_b(nnet%dims), trial(nnet%dims))
      min_b = bmin
      max_b = bmax
      diff = abs(min_b - max_b)
      do i = 1, nnet%dims
         vect(i) = grnd()
      end do
      vect_denorm = min_b + vect*diff

      ! get matrices for this trial vector
      call vector_to_matrices(vect_denorm, nnet)

      ! the choice of the first realization here is arbitrary
      call network_forward(nnet, ysimd(:, :, 1), AL, .true., nnet%norm)
      call calc_expsill(AL, sill)
      call indicator_transform(AL, thresholds, ndata, ncut, AL_i, ivars)

      ! initilaize a starting value for each component
      if (vario .gt. 0) then
         call obj_vario(AL, sill, objt_vario)
         objinit(1) = objt_vario
      end if
      if (ivario .gt. 0) then
         call obj_ivario(AL_i, ivars, objt_ivario)
         objinit(2) = objt_ivario
      end if
      if (runs .gt. 0) then
         call obj_runs(AL_i, objt_runs)
         objinit(3) = objt_runs
      end if
      if (npoint .gt. 0) then
         call obj_npoint(AL_i, objt_npt)
         objinit(4) = objt_npt
      end if

      ! if (vario .gt. 0) then
      !    call obj_vario(AL, sill, objt_vario)
      !    objinit%vario(1) = objt_vario
      ! end if

      ! if (ivario .gt. 0) then
      !    call obj_ivario(AL_i, ivars, objt_ivario)
      !    objinit%ivario(1) = objt_ivario
      ! end if

      ! if (runs .gt. 0) then
      !    do i = 1, ncut
      !       mse = 0.d0
      !       cumruns = 0
      !       do j = 1, ndh
      !          iarr = AL_i(udhidx(j) + 1:udhidx(j + 1), i)
      !          call binary_runs(iarr, maxrun, expruns)
      !          cumruns = cumruns + expruns
      !       end do
      !       objinit%runs(i) = sum((dble(target_runs(:, i)) - dble(cumruns))**2)
      !    end do
      ! end if

      ! if (npoint .gt. 0) then
      !    do i = 1, ncut
      !       mse = 0.d0
      !       call npoint_connect(AL_i(:, i), nstep, ndh, udhidx, phi_n)
      !       objinit%npoint(i) = sum((target_npoint(:, i) - phi_n)**2)
      !    end do
      ! end if

      ! iterate over the pertubations
      do i = 1, MAXPERT

         ! generate a random vector within the bounds
         do j = 1, size(vect)
            trial(j) = grnd()
         end do
         trial_denorm = min_b + trial*diff

         ! get matrices for this trial vector
         call vector_to_matrices(trial_denorm, nnet)

         ! evalute the random vector
         call network_forward(nnet, ysimd(:, :, 1), AL, .true., nnet%norm)
         call calc_expsill(AL, sill)
         call indicator_transform(AL, thresholds, ndata, ncut, AL_i, ivars)

         if (vario .gt. 0) then
            call obj_vario(AL, sill, objt_vario)
            if (objt_vario .lt. 0.0) objt_vario = objinit(1)
            objdelta(1) = objdelta(1) + abs(objinit(1) - objt_vario)
         end if

         if (ivario .gt. 0) then
            call obj_ivario(AL_i, ivars, objt_ivario)
            if (objt_ivario .lt. 0.0) objt_ivario = objinit(2)
            objdelta(2) = objdelta(2) + abs(objinit(2) - objt_ivario)
         end if

         if (runs .gt. 0) then
            call obj_runs(AL_i, objt_runs)
            if (objt_runs .lt. 0.0) objt_runs = objinit(3)
            objdelta(3) = objdelta(3) + abs(objinit(3) - objt_runs)
         end if

         if (npoint .gt. 0) then
            call obj_npoint(AL_i, objt_npt)
            if (objt_npt .lt. 0.0) objt_npt = objinit(4)
            objdelta(4) = objdelta(4) + abs(objinit(4) - objt_npt)
         end if

         ! if (vario .gt. 0) then
         !    call obj_vario(AL, sill, objt_vario)
         !    if (objt_vario .lt. 0.0) objt_vario = objinit%vario(1)
         !    objdelta%vario(1) = objdelta%vario(1) + &
         !                        abs(objinit%vario(1) - objt_vario)
         ! end if

         ! if (ivario .gt. 0) then
         !    call obj_ivario(AL_i, ivars, objt_ivario)
         !    if (objt_ivario .lt. 0.0) objt_ivario = objinit%ivario(1)
         !    objdelta%ivario(1) = objdelta%ivario(1) + &
         !                         abs(objinit%ivario(1) - objt_ivario)
         ! end if

         ! if (runs .gt. 0) then
         !    do ic = 1, ncut
         !       mse = 0.d0
         !       cumruns = 0
         !       do j = 1, ndh
         !          iarr = AL_i(udhidx(j) + 1:udhidx(j + 1), ic)
         !          call binary_runs(iarr, maxrun, expruns)
         !          cumruns = cumruns + expruns
         !       end do
         !       mse = sum((dble(target_runs(:, ic)) - dble(cumruns))**2)
         !       if (mse .lt. 0.0) mse = objinit%runs(ic)
         !       objdelta%runs(ic) = objdelta%runs(ic) + abs(objinit%runs(ic) - mse)
         !    end do
         ! end if

         ! if (npoint .gt. 0) then
         !    do ic = 1, ncut
         !       mse = 0.d0
         !       call npoint_connect(AL_i(:, ic), nstep, ndh, udhidx, phi_n)
         !       mse = sum((target_npoint(:, ic) - phi_n)**2)
         !       if (mse .lt. 0.0) mse = objinit%npoint(ic)
         !       objdelta%npoint(ic) = objdelta%npoint(ic) + &
         !                             abs(objinit%npoint(ic) - mse)
         !    end do
         ! end if

      end do

      ! scale objective components
      if (vario .gt. 0) objscale(1) = MAXPERT/objdelta(1)
      if (ivario .gt. 0) objscale(2) = MAXPERT/objdelta(2)
      if (runs .gt. 0) objscale(3) = MAXPERT/objdelta(3)
      if (npoint .gt. 0) objscale(4) = MAXPERT/objdelta(4)

      ! if (vario .gt. 0) objscale%vario(1) = MAXPERT/objdelta%vario(1)
      ! if (ivario .gt. 0) objscale%ivario(1) = MAXPERT/objdelta%ivario(1)
      ! do ic = 1, ncut
      !    if (runs .gt. 0) then
      !       objscale%runs(ic) = MAXPERT/objdelta%runs(ic)
      !    end if
      !    if (npoint .gt. 0) then
      !       objscale%npoint(ic) = MAXPERT/objdelta%npoint(ic)
      !    end if
      ! end do

      ! user defined scaling if required
      objscale(1) = userfac(1)*objscale(1)
      objscale(2) = userfac(2)*objscale(2)
      objscale(3) = userfac(3)*objscale(3)
      objscale(4) = userfac(4)*objscale(4)

      ! objscale%vario(1) = userfac(1)*objscale%vario(1)
      ! objscale%ivario(2) = userfac(2)*objscale%ivario(2)
      ! do ic = 1, ncut
      !    if (runs .gt. 0) then
      !       objscale%runs(ic) = userfac(3)*objscale%runs(ic)
      !    end if
      !    if (npoint .gt. 0) then
      !       objscale%npoint(ic) = userfac(3)*objscale%npoint(ic)
      !    end if
      ! end do

   end subroutine obj_scale

   subroutine obj_nmr(v, gobjt)

      ! network model of regionalization objective
      ! returns scalar expected objective value

      real(8), intent(in) :: v(:) ! trial vector
      real(8), intent(out) :: gobjt ! global temp obj value
      real(8) :: reg
      integer :: ireal

      gobjt = 0.d0
      objt_vario = 0.d0
      objt_ivario = 0.d0
      objt_runs = 0.d0
      objt_npt = 0.d0

      ! get matrices for this trial vector
      call vector_to_matrices(v, nnet)

      ! calculate regularization values if required
      call calc_regularization(nnet, reg)

      do ireal = 1, nreals

         call network_forward(nnet, ysimd(:, :, ireal), AL, .true., nnet%norm)
         call calc_expsill(AL, sill)
         call indicator_transform(AL, thresholds, ndata, ncut, AL_i, ivars)

         if (vario .gt. 0) call obj_vario(AL, sill, objt_vario)
         if (ivario .gt. 0) call obj_ivario(AL_i, ivars, objt_ivario)
         if (runs .gt. 0) call obj_runs(AL_i, objt_runs)
         if (npoint .gt. 0) call obj_npoint(AL_i, objt_npt)

         gobjt = gobjt + objt_vario + objt_ivario + objt_runs + objt_npt

      end do

      gobjt = (gobjt + reg)/nreals

   end subroutine obj_nmr

   subroutine obj_nmr_vect(v, gobjt)

      ! network model of regionalization objective
      ! returns scalar expected objective value

      ! this version has no components, rather it is
      ! data reconstruction error

      real(8), intent(in) :: v(:) ! trial vector
      real(8), intent(out) :: gobjt ! global temp obj value
      real(8) :: sqerr, reg
      integer :: ireal

      gobjt = 0.d0

      ! get matrices for this trial vector
      call vector_to_matrices(v, nnet)

      ! calculate regularization values if required
      call calc_regularization(nnet, reg)

      do ireal = 1, nreals

         call network_forward(nnet, ysimd(:, :, ireal), AL, .true., nnet%norm)
         sqerr = sum((AL - var)**2)
         gobjt = gobjt + sqerr

      end do

      gobjt = (gobjt + reg)/nreals

   end subroutine obj_nmr_vect

   subroutine pobj_nmr(v, net, simd, gobjt)

      ! parallel network model of regionalization objective
      ! returns scalar expected objective value

      real(8), intent(in) :: v(:) ! trial vector
      type(network), intent(inout) :: net
      real(8), intent(in) :: simd(:, :, :)
      real(8), intent(out) :: gobjt ! global temp obj value
      real(8) :: tobj_vario, tobj_ivario, tobj_runs, tobj_npt
      real(8) :: mix(ndata)
      integer :: imix(ndata, ncut)
      real(8) :: reg
      real(8) :: expsill, iexpsills(ncut)
      integer :: ireal

      gobjt = 0.d0
      tobj_vario = 0.d0
      tobj_ivario = 0.d0
      tobj_runs = 0.d0
      tobj_npt = 0.d0

      ! get matrices for this trial vector
      call vector_to_matrices(v, net)

      ! calculate regularization values if required
      call calc_regularization(net, reg)

      do ireal = 1, nreals

         call network_forward(net, simd(:, :, ireal), mix, .true., nnet%norm)
         call calc_expsill(mix, expsill)
         call indicator_transform(mix, thresholds, ndata, ncut, imix, iexpsills)

         if (vario .gt. 0) call obj_vario(mix, expsill, tobj_vario)
         if (ivario .gt. 0) call obj_ivario(imix, iexpsills, tobj_ivario)
         if (runs .gt. 0) call obj_runs(imix, tobj_runs)
         if (npoint .gt. 0) call obj_npoint(imix, tobj_npt)

         gobjt = gobjt + tobj_vario + tobj_ivario + tobj_runs + tobj_npt

      end do

      gobjt = (gobjt + reg)/nreals

   end subroutine pobj_nmr

   subroutine pobj_nmr_vect(v, net, simd, gobjt)

      ! parallel network model of regionalization objective
      ! returns scalar expected objective value

      ! this version has no components, rather it is
      ! data reconstruction error

      real(8), intent(in) :: v(:) ! trial vector
      type(network), intent(inout) :: net
      real(8), intent(in) :: simd(:, :, :)
      real(8), intent(out) :: gobjt ! global temp obj value
      real(8) :: mix(ndata)
      real(8) :: sqerr, reg
      integer :: ireal

      gobjt = 0.d0

      ! get matrices for this trial vector
      call vector_to_matrices(v, net)

      ! calculate regularization values if required
      call calc_regularization(net, reg)

      do ireal = 1, nreals

         call network_forward(net, simd(:, :, ireal), mix, .true., nnet%norm)
         sqerr = sum((mix - var)**2)
         gobjt = gobjt + sqerr

      end do

      gobjt = (gobjt + reg)/nreals

   end subroutine pobj_nmr_vect

   subroutine obj_vario(mix, expsill, objt)

      ! continuous variogram component
      real(8), intent(in) :: mix(:), expsill
      real(8), intent(inout) :: objt

      real(8), allocatable :: expvario(:)
      real(8) :: mse
      integer :: i

      objt = 0.d0

      do i = 1, ndir
         call update_vario(heads%dirs(i), tails%dirs(i), mix, expvario, expsill)
         call vario_mse(expvario, target_vario%dirs(i)%vlags, &
                        varlagdist%dirs(i)%vlags, dble(idwpow), mse)
         objt = objt + mse
      end do

      objt = objt*objscale(1)

   end subroutine obj_vario

   subroutine obj_ivario(imix, iexpsills, objt)

      ! indicator variogram component
      integer, intent(in) :: imix(:, :)
      real(8), intent(in) :: iexpsills(:)
      real(8), intent(inout) :: objt

      real(8), allocatable :: expivario(:)
      real(8) :: mse
      integer :: i, j

      objt = 0.d0

      do j = 1, ncut
         do i = 1, ndir
            call update_vario(heads%dirs(i), tails%dirs(i), dble(imix(:, j)), &
                              expivario, iexpsills(j))
            call vario_mse(expivario, target_ivario%cuts(j)%dirs(i)%vlags, &
                           varlagdist%dirs(i)%vlags, dble(idwpow), mse)
            objt = objt + mse
         end do
      end do

      objt = objt*objscale(2)

   end subroutine obj_ivario

   subroutine obj_runs(imix, objt)

      ! cumulative run frequency component
      integer, intent(in) :: imix(:, :)
      real(8), intent(inout) :: objt

      integer, allocatable :: iarr(:), expruns(:)
      integer :: cumruns(maxrun)
      real(8) :: mse
      integer :: i, j

      objt = 0.d0

      do i = 1, ncut
         cumruns = 0
         do j = 1, ndh
            iarr = imix(udhidx(j) + 1:udhidx(j + 1), i)
            call binary_runs(iarr, maxrun, expruns)
            cumruns = cumruns + expruns
         end do
         mse = sum((dble(target_runs(:, i)) - dble(cumruns))**2)
         objt = objt + mse
      end do

      objt = objt*objscale(3)

   end subroutine obj_runs

   subroutine obj_npoint(imix, objt)

      ! npoint connectivity component
      integer, intent(in) :: imix(:, :)
      real(8), intent(inout) :: objt

      real(8), allocatable :: phi_n(:)
      real(8) :: mse
      integer :: i

      objt = 0.d0

      do i = 1, ncut
         call npoint_connect(imix(:, i), nstep, ndh, udhidx, phi_n)
         mse = sum((target_npoint(:, i) - phi_n)**2)
         objt = objt + mse
      end do

      objt = objt*objscale(4)

   end subroutine obj_npoint

   subroutine feature_importance(net, ysim, fimp)

      ! permutation feature importance for trained network

      type(network), intent(inout) :: net
      real(8), intent(in) :: ysim(:, :, :)
      real(8), allocatable, intent(inout) :: fimp(:)
      real(8), allocatable :: yperm(:, :), yp(:)
      real(8) :: ep, eo, reg
      integer, allocatable :: idxs(:)
      integer :: i, j, nfact

      ndata = size(ysim, dim=1)
      nfact = net%ld(1)

      allocate (yperm(ndata, nfact), yp(ndata))
      allocate (fimp(nfact))

      fimp = 0.d0

      idxs = [(i, i=1, ndata)]
      call vector_to_matrices(best, net)
      call calc_regularization(net, reg)

      ! main loop over input features
      do i = 1, net%ld(1)

         ! permute each realization of the ith feature
         do j = 1, nreals

            ! get the permuted error
            call shuffle(idxs)
            yperm = ysim(:, :, j)
            yp = yperm(idxs, i)
            yperm(:, i) = yp

            call network_forward(net, yperm, AL, .true., nnet%norm)
            call calc_expsill(AL, sill)
            call indicator_transform(AL, thresholds, ndata, ncut, AL_i, ivars)

            if (vario .gt. 0) call obj_vario(AL, sill, objt_vario)
            if (ivario .gt. 0) call obj_ivario(AL_i, ivars, objt_ivario)
            if (runs .gt. 0) call obj_runs(AL_i, objt_runs)
            if (npoint .gt. 0) call obj_npoint(AL_i, objt_npt)

            ep = objt_vario + objt_ivario + objt_runs + objt_npt

            ! get the original error
            call network_forward(net, ysim(:, :, j), AL, .true., nnet%norm)
            call calc_expsill(AL, sill)
            call indicator_transform(AL, thresholds, ndata, ncut, AL_i, ivars)

            if (vario .gt. 0) call obj_vario(AL, sill, objt_vario)
            if (ivario .gt. 0) call obj_ivario(AL_i, ivars, objt_ivario)
            if (runs .gt. 0) call obj_runs(AL_i, objt_runs)
            if (npoint .gt. 0) call obj_npoint(AL_i, objt_npt)

            eo = objt_vario + objt_ivario + objt_runs + objt_npt

            ! sum the feature importance
            fimp(i) = fimp(i) + (ep - eo)

         end do
      end do

      ! get the expected value
      fimp = fimp/nreals

   end subroutine feature_importance

end module objective_mod
