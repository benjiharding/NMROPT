module objective_mod

   use geostat
   use sequences_mod, only: binary_runs, npoint_connect
   use vario_mod, only: update_vario, vario_mse, indicator_transform, &
                        vario_pairs, varmodelpts, set_sill, calc_expsill
   use network_mod, only: network_forward, vector_to_matrices, &
                          calc_regularization, build_refcdf, &
                          transform_to_refcdf
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
      real(8) :: q
      integer :: maxlags
      integer :: i, ii, j, k, nl

      ! allocate data arrays
      allocate (iz(ndata, ncut), isills(ncut))
      allocate (AL(ndata), AL_i(ndata, ncut))

      ! indicator transform of input data
      call indicator_transform(var, thresholds, ndata, ncut, iz)

      ! get vario sills for standardization

      sill = 1.d0
      ! call calc_expsill(var, sill, vtype=1)

      do i = 1, ncut
         q = gcum(thresholds(i))
         isills(i) = q*(1.d0 - q)
         ! call calc_expsill(var, isills(i), vtype=2, cut=thresholds(i))
      end do

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

         if (size(valid_idxs, dim=1) .lt. 1) stop "No defined experimental lags"

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
         allocate (avg_vario%dirs(ndir))
         do i = 1, ndir
            maxlags = size(varlagdist%dirs(i)%vlags)
            allocate (target_vario%dirs(i)%vlags(maxlags))
            allocate (avg_vario%dirs(i)%vlags(maxlags))
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
         allocate (avg_ivario%cuts(ncut))
         do j = 1, ncut
            allocate (target_ivario%cuts(j)%dirs(ndir))
            allocate (avg_ivario%cuts(j)%dirs(ndir))
            ! calculate the model points
            do k = 1, ndir
               maxlags = size(varlagdist%dirs(k)%vlags)
               allocate (target_ivario%cuts(j)%dirs(k)%vlags(maxlags))
               allocate (avg_ivario%cuts(j)%dirs(k)%vlags(maxlags))
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
         allocate (avg_runs(maxrun, ncut))
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
         allocate (avg_npoint(nstep, ncut))
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
      write (*, "(A*(g14.8,1x))") "Data component: ", objscale(5)

   end subroutine init_objective

   subroutine obj_scale()

      ! scale objective components by looping over MAXPERT iterations and
      ! tracking the how much each component changes

      integer, parameter :: MAXPERT = 1000
      real(8), allocatable :: vect(:), vect_denorm(:), min_b(:), max_b(:), diff(:)
      real(8), allocatable :: trial(:), trial_denorm(:)
      real(8) :: objinit(5), objdelta(5)
      integer :: i, j

      objinit = 0.d0
      objdelta = 0.d0

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

      ! build transform for this trial vector and
      ! use the same ttable for all realizations
      call build_refcdf(nsamp, yref, nnet, ttable)

      ! the choice of the first realization here is arbitrary
      call network_forward(nnet, ysimd(:, :, 1), AL, .true., nnet%norm, ttable)
      call indicator_transform(AL, thresholds, ndata, ncut, AL_i)

      ! initilaize a starting value for each component
      if (vario .gt. 0) then
         call obj_vario(AL, sill, objt_vario)
         objinit(1) = objt_vario
      end if
      if (ivario .gt. 0) then
         call obj_ivario(AL_i, isills, objt_ivario)
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
      call obj_data(AL, objt_data)
      objinit(5) = objt_data

      ! iterate over the pertubations
      do i = 1, MAXPERT

         ! generate a random vector within the bounds
         do j = 1, size(vect)
            trial(j) = grnd()
         end do
         trial_denorm = min_b + trial*diff

         ! get matrices for this trial vector
         call vector_to_matrices(trial_denorm, nnet)

         ! build transform for this trial vector and
         ! use the same ttable for all realizations
         call build_refcdf(nsamp, yref, nnet, ttable)

         ! evalute the random vector
         call network_forward(nnet, ysimd(:, :, 1), AL, .true., nnet%norm, ttable)
         call indicator_transform(AL, thresholds, ndata, ncut, AL_i)

         if (vario .gt. 0) then
            call obj_vario(AL, sill, objt_vario)
            if (objt_vario .lt. 0.0) objt_vario = objinit(1)
            objdelta(1) = objdelta(1) + abs(objinit(1) - objt_vario)
         end if

         if (ivario .gt. 0) then
            call obj_ivario(AL_i, isills, objt_ivario)
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

         call obj_data(AL, objt_data)
         if (objt_data .lt. 0.0) objt_data = objinit(5)
         objdelta(5) = objdelta(5) + abs(objinit(5) - objt_data)

      end do

      ! scale objective components
      if (vario .gt. 0) objscale(1) = MAXPERT/objdelta(1)
      if (ivario .gt. 0) objscale(2) = MAXPERT/objdelta(2)
      if (runs .gt. 0) objscale(3) = MAXPERT/objdelta(3)
      if (npoint .gt. 0) objscale(4) = MAXPERT/objdelta(4)
      objscale(5) = MAXPERT/objdelta(5)

      ! user defined scaling if required
      objscale(1) = userfac(1)*objscale(1)
      objscale(2) = userfac(2)*objscale(2)
      objscale(3) = userfac(3)*objscale(3)
      objscale(4) = userfac(4)*objscale(4)

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
      objt_data = 0.d0

      ! get matrices for this trial vector
      call vector_to_matrices(v, nnet)

      ! build transform for this trial vector and
      ! use the same ttable for all realizations
      call build_refcdf(nsamp, yref, nnet, ttable)

      ! calculate regularization values if required
      call calc_regularization(nnet, reg)

      do ireal = 1, nreals

         call network_forward(nnet, ysimd(:, :, ireal), AL, .true., nnet%norm, ttable)
         call indicator_transform(AL, thresholds, ndata, ncut, AL_i)

         if (vario .gt. 0) call obj_vario(AL, sill, objt_vario)
         if (ivario .gt. 0) call obj_ivario(AL_i, isills, objt_ivario, threshwt)
         if (runs .gt. 0) call obj_runs(AL_i, objt_runs, threshwt)
         if (npoint .gt. 0) call obj_npoint(AL_i, objt_npt, threshwt)
         call obj_data(AL, objt_data)

         gobjt = gobjt + objt_vario + objt_ivario + objt_runs + objt_npt + objt_data

      end do

      gobjt = (gobjt + reg)/nreals

   end subroutine obj_nmr

   subroutine obj_nmr_avg(v, gobjt)

      ! network model of regionalization objective
      ! returns scalar expected objective value

      real(8), intent(in) :: v(:) ! trial vector
      real(8), intent(out) :: gobjt ! global temp obj value
      real(8) :: reg

      gobjt = 0.d0
      objt_vario = 0.d0
      objt_ivario = 0.d0
      objt_runs = 0.d0
      objt_npt = 0.d0
      objt_data = 0.d0

      ! get matrices for this trial vector
      call vector_to_matrices(v, nnet)

      ! build transform for this trial vector and
      ! use the same ttable for all realizations
      call build_refcdf(nsamp, yref, nnet, ttable)

      ! calculate regularization values if required
      call calc_regularization(nnet, reg)

      ! calculate component averages
      call avg_obj(ysimd, nreals, vario, ivario, runs, npoint)

      if (vario .gt. 0) call obj_vario_avg(avg_vario, objt_vario)
      if (ivario .gt. 0) call obj_ivario_avg(avg_ivario, objt_ivario, threshwt)
      if (runs .gt. 0) call obj_runs_avg(avg_runs, objt_runs, threshwt)
      if (npoint .gt. 0) call obj_npoint_avg(avg_npoint, objt_npt, threshwt)
      ! call obj_data(mix, objt_data)

      gobjt = objt_vario + objt_ivario + objt_runs + objt_npt + objt_data
      gobjt = gobjt + reg

   end subroutine obj_nmr_avg

   subroutine pobj_nmr(v, net, simd, gobjt)

      ! parallel network model of regionalization objective
      ! returns scalar expected objective value

      real(8), intent(in) :: v(:) ! trial vector
      type(network), intent(inout) :: net
      real(8), intent(in) :: simd(:, :, :)
      real(8), intent(out) :: gobjt ! global temp obj value
      real(8) :: tobj_vario, tobj_ivario, tobj_runs, tobj_npt, tobj_data
      real(8) :: mix(ndata)
      integer :: imix(ndata, ncut)
      real(8) :: reg
      ! real(8) :: expsill, iexpsills(ncut)
      integer :: ireal

      gobjt = 0.d0
      tobj_vario = 0.d0
      tobj_ivario = 0.d0
      tobj_runs = 0.d0
      tobj_npt = 0.d0
      tobj_data = 0.d0

      ! get matrices for this trial vector (updates net object)
      call vector_to_matrices(v, net)

      ! build transform for this trial vector and
      ! use the same ttable for all realizations
      call build_refcdf(nsamp, yref, net, ttable)

      ! calculate regularization values if required
      call calc_regularization(net, reg)

      do ireal = 1, nreals

         call network_forward(net, simd(:, :, ireal), mix, .true., net%norm, ttable)
         call indicator_transform(mix, thresholds, ndata, ncut, imix)

         if (vario .gt. 0) call obj_vario(mix, sill, tobj_vario)
         if (ivario .gt. 0) call obj_ivario(imix, isills, tobj_ivario, threshwt)
         if (runs .gt. 0) call obj_runs(imix, tobj_runs, threshwt)
         if (npoint .gt. 0) call obj_npoint(imix, tobj_npt, threshwt)
         ! call obj_data(mix, tobj_data)

         gobjt = gobjt + tobj_vario + tobj_ivario + tobj_runs + tobj_npt + tobj_data

      end do

      gobjt = (gobjt + reg)/nreals

   end subroutine pobj_nmr

   subroutine pobj_nmr_avg(v, net, simd, gobjt)

      ! parallel network model of regionalization objective
      ! returns scalar expected objective value

      real(8), intent(in) :: v(:) ! trial vector
      type(network), intent(inout) :: net
      real(8), intent(in) :: simd(:, :, :)
      real(8), intent(out) :: gobjt ! global temp obj value
      real(8) :: tobj_vario, tobj_ivario, tobj_runs, tobj_npt, tobj_data
      real(8) :: reg

      gobjt = 0.d0
      tobj_vario = 0.d0
      tobj_ivario = 0.d0
      tobj_runs = 0.d0
      tobj_npt = 0.d0
      tobj_data = 0.d0

      ! get matrices for this trial vector (updates net object)
      call vector_to_matrices(v, net)

      ! build transform for this trial vector and
      ! use the same ttable for all realizations
      call build_refcdf(nsamp, yref, net, ttable)

      ! calculate regularization values if required
      call calc_regularization(net, reg)

      ! calculate component averages
      call avg_obj(simd, nreals, vario, ivario, runs, npoint)

      if (vario .gt. 0) call obj_vario_avg(avg_vario, tobj_vario)
      if (ivario .gt. 0) call obj_ivario_avg(avg_ivario, tobj_ivario, threshwt)
      if (runs .gt. 0) call obj_runs_avg(avg_runs, tobj_runs, threshwt)
      if (npoint .gt. 0) call obj_npoint_avg(avg_npoint, tobj_npt, threshwt)
      ! call obj_data(mix, tobj_data)

      gobjt = tobj_vario + tobj_ivario + tobj_runs + tobj_npt + tobj_data
      gobjt = gobjt + reg

   end subroutine pobj_nmr_avg

   subroutine avg_obj(reals, nreals, avg1, avg2, avg3, avg4)

      ! calculate the average objective for each component
      ! this is done in a single subroutine call to prevent calling
      ! network_forward 4 times

      ! parameters
      real(8), intent(in) :: reals(:, :, :)
      integer, intent(in) :: nreals
      integer, intent(in) :: avg1, avg2, avg3, avg4

      ! locals
      real(8) :: mix(ndata)
      integer :: imix(ndata, ncut)
      real(8), allocatable :: expvario(:), expivario(:)
      integer, allocatable :: iarr(:), expruns(:)
      integer :: cumruns(maxrun)
      real(8), allocatable :: phi_n(:)
      integer :: ireal, i, j

      ! intialize all arrays to zero
      do i = 1, ndir
         if (avg1 .gt. 0) avg_vario%dirs(i)%vlags = 0.d0
         do j = 1, ncut
            if (avg2 .gt. 0) avg_ivario%cuts(j)%dirs(i)%vlags = 0.d0
            if (avg3 .gt. 0) avg_runs(:, j) = 0.d0
            if (avg4 .gt. 0) avg_npoint(:, j) = 0.d0
         end do
      end do

      ! iterate over realizations
      do ireal = 1, nreals

         call network_forward(nnet, reals(:, :, ireal), mix, .true., nnet%norm, ttable)
         call indicator_transform(mix, thresholds, ndata, ncut, imix)

         ! vario
         if (avg1 .gt. 0) then
            do i = 1, ndir
               call update_vario(heads%dirs(i), tails%dirs(i), mix, expvario, sill)
               avg_vario%dirs(i)%vlags = avg_vario%dirs(i)%vlags + expvario
            end do
         end if

         ! ivario
         if (avg2 .gt. 0) then
            do j = 1, ncut
               do i = 1, ndir
                  call update_vario(heads%dirs(i), tails%dirs(i), dble(imix(:, j)), &
                                    expivario, isills(j))
                  avg_ivario%cuts(j)%dirs(i)%vlags = avg_ivario%cuts(j)%dirs(i)%vlags &
                                                     + expivario
               end do
            end do
         end if

         ! runs
         if (avg3 .gt. 0) then
            do i = 1, ncut
               cumruns = 0
               do j = 1, ndh
                  iarr = imix(udhidx(j) + 1:udhidx(j + 1), i)
                  call binary_runs(iarr, maxrun, expruns)
                  cumruns = cumruns + expruns
               end do
               avg_runs(:, i) = avg_runs(:, i) + cumruns
            end do
         end if

         ! npoint
         if (avg4 .gt. 0) then
            do i = 1, ncut
               call npoint_connect(imix(:, i), nstep, ndh, udhidx, phi_n)
               avg_npoint(:, i) = avg_npoint(:, i) + phi_n
            end do
         end if
      end do

      ! take the averages
      do i = 1, ndir
         if (avg1 .gt. 0) then
            avg_vario%dirs(i)%vlags = avg_vario%dirs(i)%vlags/nreals
         end if
         do j = 1, ncut
            if (avg2 .gt. 0) then
               avg_ivario%cuts(j)%dirs(i)%vlags = avg_ivario%cuts(j)%dirs(i)%vlags/nreals
            end if
         end do
      end do

      do j = 1, ncut
         if (avg3 .gt. 0) avg_runs(:, j) = avg_runs(:, j)/nreals
         if (avg4 .gt. 0) avg_npoint(:, j) = avg_npoint(:, j)/nreals
      end do

   end subroutine avg_obj

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

   subroutine obj_vario_avg(avg_vario, objt)

      ! continuous variogram component
      type(vario_array), intent(in) :: avg_vario
      real(8), intent(inout) :: objt

      real(8) :: mse
      integer :: i

      objt = 0.d0

      do i = 1, ndir
         call vario_mse(avg_vario%dirs(i)%vlags, target_vario%dirs(i)%vlags, &
                        varlagdist%dirs(i)%vlags, dble(idwpow), mse)
         objt = objt + mse
      end do

      objt = objt*objscale(1)

   end subroutine obj_vario_avg

   subroutine obj_ivario(imix, iexpsills, objt, tw)

      ! indicator variogram component
      integer, intent(in) :: imix(:, :)
      real(8), intent(in) :: iexpsills(:)
      real(8), intent(inout) :: objt
      real(8), optional :: tw(:) ! threshold weights

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
            if (present(tw)) then
               objt = objt + mse*threshwt(j)
            else
               objt = objt + mse
            end if
         end do
      end do

      objt = objt*objscale(2)

   end subroutine obj_ivario

   subroutine obj_ivario_avg(avg_ivario, objt, tw)

      ! indicator variogram component
      type(ivario_array), intent(in) :: avg_ivario
      real(8), intent(inout) :: objt
      real(8), optional :: tw(:) ! threshold weights

      real(8) :: mse
      integer :: i, j

      objt = 0.d0

      do j = 1, ncut
         do i = 1, ndir
            call vario_mse(avg_ivario%cuts(j)%dirs(i)%vlags, target_ivario%cuts(j)%dirs(i)%vlags, &
                           varlagdist%dirs(i)%vlags, dble(idwpow), mse)
            if (present(tw)) then
               objt = objt + mse*threshwt(j)
            else
               objt = objt + mse
            end if
         end do
      end do

      objt = objt*objscale(2)

   end subroutine obj_ivario_avg

   subroutine obj_runs(imix, objt, tw)

      ! cumulative run frequency component
      integer, intent(in) :: imix(:, :)
      real(8), intent(inout) :: objt
      real(8), optional :: tw(:) ! threshold weights

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
         if (present(tw)) then
            objt = objt + mse*threshwt(i)
         else
            objt = objt + mse
         end if
      end do

      objt = objt*objscale(3)

   end subroutine obj_runs

   subroutine obj_runs_avg(avg_runs, objt, tw)

      ! cumulative run frequency component
      real(8), intent(in) :: avg_runs(:, :)
      real(8), intent(inout) :: objt
      real(8), optional :: tw(:) ! threshold weights
      real(8) :: mse
      integer :: i, j

      objt = 0.d0

      do i = 1, ncut
         mse = sum((dble(target_runs(:, i)) - avg_runs(:, i))**2)
         if (present(tw)) then
            objt = objt + mse*threshwt(i)
         else
            objt = objt + mse
         end if
      end do

      objt = objt*objscale(3)

   end subroutine obj_runs_avg

   subroutine obj_npoint(imix, objt, tw)

      ! npoint connectivity component
      integer, intent(in) :: imix(:, :)
      real(8), intent(inout) :: objt
      real(8), optional :: tw(:) ! threshold weights

      real(8), allocatable :: phi_n(:)
      real(8) :: mse
      integer :: i

      objt = 0.d0

      do i = 1, ncut
         call npoint_connect(imix(:, i), nstep, ndh, udhidx, phi_n)
         mse = sum((target_npoint(:, i) - phi_n)**2)
         if (present(tw)) then
            objt = objt + mse*threshwt(i)
         else
            objt = objt + mse
         end if
      end do

      objt = objt*objscale(4)

   end subroutine obj_npoint

   subroutine obj_npoint_avg(avg_npoint, objt, tw)

      ! npoint connectivity component
      real(8), intent(in) :: avg_npoint(:, :)
      real(8), intent(inout) :: objt
      real(8), optional :: tw(:) ! threshold weights
      real(8) :: mse
      integer :: i

      objt = 0.d0

      do i = 1, ncut
         mse = sum((target_npoint(:, i) - avg_npoint(:, i))**2)
         if (present(tw)) then
            objt = objt + mse*threshwt(i)
         else
            objt = objt + mse
         end if
      end do

      objt = objt*objscale(4)

   end subroutine obj_npoint_avg

   subroutine obj_data(mix, objt)

      ! data component to enforece network output
      ! to be positively correlated with actual data
      real(8), intent(in) :: mix(:)
      real(8), intent(inout) :: objt
      real(8) :: rho

      objt = 0.d0

      ! penalize sum of squares
      objt = sum((mix - var)**2)
      objt = objt*objscale(5)

      ! ! penalize correlation coefficient
      ! call correlation_coefficient(mix, var, ndata, rho)
      ! objt = rho_error(rho)

      objt = objt*objscale(5)

   end subroutine obj_data

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
      call build_refcdf(nsamp, yref, net, ttable)
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

            call network_forward(net, yperm, AL, .true., net%norm, ttable)
            call indicator_transform(AL, thresholds, ndata, ncut, AL_i)

            if (vario .gt. 0) call obj_vario(AL, sill, objt_vario)
            if (ivario .gt. 0) call obj_ivario(AL_i, isills, objt_ivario)
            if (runs .gt. 0) call obj_runs(AL_i, objt_runs)
            if (npoint .gt. 0) call obj_npoint(AL_i, objt_npt)

            ep = objt_vario + objt_ivario + objt_runs + objt_npt

            ! get the original error
            call network_forward(net, ysim(:, :, j), AL, .true., net%norm, ttable)
            call indicator_transform(AL, thresholds, ndata, ncut, AL_i)

            if (vario .gt. 0) call obj_vario(AL, sill, objt_vario)
            if (ivario .gt. 0) call obj_ivario(AL_i, isills, objt_ivario)
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

   function rho_error(rho) result(err)

      real(8), intent(in) :: rho
      real(8) :: err

      if (rho .gt. 1.d0 .or. rho .lt. -1.d0) then
         write (*, *) rho
         stop
      end if

      if (rho .le. 0.d0) err = 1.d0
      if (rho .gt. 0.d0) err = 1.d0 - rho

   end function rho_error

end module objective_mod
