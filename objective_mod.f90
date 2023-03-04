module objective_mod

   use vario_mod
   use sequences_mod
   use network_mod
   use lusim_mod
   use quicksort
   use subs
   use constants

   implicit none

   ! candidate vectors
   real(8), allocatable :: AL(:) ! Gaussian mixture
   integer, allocatable :: AL_i(:, :) ! indicator transform of AL (nd, ncut)

   ! objective function targets
   type(vario_array) :: target_vario ! (ndir*nlags)
   type(ivario_array) :: target_ivario ! (ncut, ndir*nlags)

   integer, allocatable :: target_runs(:, :) ! (ncut, maxrun)
   integer, allocatable :: target_npoint(:, :) ! (ncut, nstep)

   ! objective function scaling
   real(8) :: objscale(4) ! fobj scaling parameters
   real(8) :: objt_vario, objt_ivario, objt_runs, objt_npt ! temp obj values

   ! drillhole parameters
   integer, allocatable :: dhids(:), dhlens(:) ! dhids and length
   integer, allocatable :: udhids(:), udhidx(:) ! unique dhids
   integer :: ndh ! number of drillholes
   real(8), allocatable :: x(:), y(:), z(:) ! coordinates
   integer, allocatable :: iz(:, :) ! indicator transform of 'var'
   integer :: itrans ! nscore flag

   ! variogram parameters
   integer, allocatable :: pairs(:, :) ! (npairs*ndir, 3)
   integer, allocatable :: head(:), tail(:), lag(:), dir(:)
   type(lag_array) :: heads, tails
   type(vario_array) :: varlagdist

   real(8), allocatable :: thresholds(:) ! for indicator transform
   real(8), allocatable :: ivars(:) ! indicator sills (ncut)
   integer, allocatable :: nlags(:) ! multiple directions
   integer, allocatable :: stride(:, :) ! (nlags, ndir)
   real(8), allocatable :: azm(:), atol(:), bandh(:), dip(:), dtol(:), &
                           bandv(:), tilt(:), lagdis(:), lagtol(:)
   real(8), allocatable :: c0(:), cc(:), ang1(:), ang2(:), ang3(:), aa(:), &
                           anis1(:), anis2(:), ahmin(:), avert(:)
   integer, allocatable :: nst(:), it(:)
   real(8), allocatable :: ic0(:), icc(:, :), iang1(:, :), iang2(:, :), &
                           iang3(:, :), iaa(:, :), ianis1(:, :), ianis2(:, :), &
                           iahmin(:, :), iavert(:, :)
   integer, allocatable :: inst(:), iit(:, :)
   real(8) :: aa1, aa2, idwpow
   integer :: nvarg, nivarg, isill, ndir, ncut
   integer, allocatable :: udiridx(:), ulagidx(:) ! array indices
   integer :: e, f, g, h, nc

   ! run parameters
   integer :: maxrun, nstep ! maximum runs and connected steps

   ! optimization initialization
   real(8) :: bmin, bmax
   real(8) :: userfac(4)

   ! boolean flags
   integer :: vario
   integer :: ivario
   integer :: runs
   integer :: npoint
   integer :: runs_above
   integer :: conn_above

   ! output file
   integer :: lprs = 6

contains

   subroutine init_objective()

      ! initialize objective targets, values and scaling components

      integer :: numpairs(ndir)
      integer, allocatable :: tmppairs(:, :), tmpruns(:), tmparr(:)
      integer, allocatable :: headt(:), tailt(:), lagt(:), lag_idxs(:)
      real(8), allocatable :: tmpbins(:)
      type(vario_array) :: varazm, vardip

      integer :: maxlags
      integer :: i, j, k, nl

      ! data coordinates
      x = xyz(1, :)
      y = xyz(2, :)
      z = xyz(3, :)

      ! some allocations
      allocate (iz(ndata, ncut), ivars(ncut))
      allocate (AL(ndata), AL_i(ndata, ncut))
      allocate (stride(maxval(nlags), ndir))

      ! indicator transform of 'var'
      call indicator_transform(var, thresholds, ndata, ncut, iz, ivars)

      ! get the start/end dh indices in the mixture array
      allocate (udhidx(ndh + 1))
      udhidx = cumsum(dhlens)

      ! get variogram lag and pair data
      allocate (heads%dirs(ndir), tails%dirs(ndir))
      allocate (varazm%dirs(ndir), vardip%dirs(ndir), varlagdist%dirs(ndir))

      do i = 1, ndir
         call vario_pairs(x, y, z, azm(i), atol(i), bandh(i), dip(i), &
                          dtol(i), bandv(i), nlags(i), lagdis(i), lagtol(i), &
                          ndata, tmppairs, tmpbins)
         headt = tmppairs(:, 1)
         tailt = tmppairs(:, 2)
         lagt = tmppairs(:, 3)
         numpairs(i) = size(lagt)
         nl = maxval(lagt) - minval(lagt) + 1

         ! only allocate the number of defined lags
         allocate (varazm%dirs(i)%vlags(nl))
         allocate (vardip%dirs(i)%vlags(nl))
         allocate (varlagdist%dirs(i)%vlags(nl))
         allocate (heads%dirs(i)%lags(nl))
         allocate (tails%dirs(i)%lags(nl))

         do j = 1, nl

            ! experiemntal parameters
            varazm%dirs(i)%vlags(j) = azm(i)
            vardip%dirs(i)%vlags(j) = dip(i)
            varlagdist%dirs(i)%vlags(j) = tmpbins(j)

            ! lag indices
            lag_idxs = pack([(k, k=1, size(tmppairs, dim=1))], &
                            tmppairs(:, 3) .eq. j)
            allocate (heads%dirs(i)%lags(j)%idxs(size(lag_idxs)))
            allocate (tails%dirs(i)%lags(j)%idxs(size(lag_idxs)))
            do k = 1, size(lag_idxs)
               heads%dirs(i)%lags(j)%idxs(k) = headt(lag_idxs(k))
               tails%dirs(i)%lags(j)%idxs(k) = tailt(lag_idxs(k))
            end do

         end do
      end do

      ! establish stride to reduce total number of exp variogram pairs
      stride = 1

      ! ! write out pairs if debugging
      ! open (lprs, file="vario_pairs.out", status="UNKNOWN")
      ! write (lprs, "(a28)") "Experimental Variogram Pairs"
      ! write (lprs, "(i1)") 4
      ! write (lprs, "(a8)") "Head IDX"
      ! write (lprs, "(a8)") "Tail IDX"
      ! write (lprs, "(a7)") "Lag IDX"
      ! write (lprs, "(a7)") "Dir IDX"
      ! do i = 1, ndir
      !    do j = 1, size(heads%dirs(i)%lags)
      !       do k = 1, size(heads%dirs(i)%lags(j)%idxs)
      !          write (lprs, "(*(i5))") heads%dirs(i)%lags(j)%idxs(k), &
      !             tails%dirs(i)%lags(j)%idxs(k), j, i
      !       end do
      !    end do
      ! end do
      ! close (lprs)

      ! variogram target
      allocate (target_vario%dirs(ndir))
      ahmin = anis1*aa
      avert = anis2*aa
      do i = 1, ndir
         maxlags = size(varlagdist%dirs(i)%vlags)
         allocate (target_vario%dirs(i)%vlags(maxlags))
         call varmodelpts(nst(1), c0(1), it, cc, ang1, ang2, ang3, aa, ahmin, &
                          avert, maxlags, varlagdist%dirs(i)%vlags, &
                          varazm%dirs(i)%vlags, vardip%dirs(i)%vlags, &
                          target_vario%dirs(i)%vlags)
      end do

      ! indicator variogram target
      allocate (target_ivario%cuts(ncut))
      allocate (iahmin(ncut, MAXNST), iavert(ncut, MAXNST))

      do j = 1, ncut

         allocate (target_ivario%cuts(j)%dirs(ndir))

         ! rescale sill parameters
         ic0(j) = ic0(j)*ivars(j)
         do i = 1, inst(j)
            icc(j, i) = icc(j, i)*ivars(j)
         end do

         ! calculate the model points
         iahmin(j, :) = ianis1(j, :)*iaa(j, :)
         iavert(j, :) = ianis2(j, :)*iaa(j, :)

         do k = 1, ndir
            maxlags = size(varlagdist%dirs(k)%vlags)
            allocate (target_ivario%cuts(j)%dirs(k)%vlags(maxlags))
            call varmodelpts(inst(j), ic0(j), iit(j, :), icc(j, :), iang1(j, :), &
                             iang2(j, :), iang3(j, :), iaa(j, :), iahmin(j, :), &
                             iavert(j, :), maxlags, varlagdist%dirs(k)%vlags, &
                             varazm%dirs(k)%vlags, vardip%dirs(k)%vlags, &
                             target_ivario%cuts(j)%dirs(k)%vlags)
         end do
      end do

      ! runs target
      allocate (target_runs(maxrun, ncut))
      target_runs(:, :) = 0
      do i = 1, ncut
         do j = 2, ndh
            tmparr = iz(udhidx(j - 1) + 1:udhidx(j), i)
            call calculate_runs(tmparr, maxrun, tmpruns)
            target_runs(:, i) = target_runs(:, i) + tmpruns
         end do
      end do

      ! scale the components
      write (*, *) " "
      write (*, *) "Scaling objective components..."

      call obj_scale

      write (*, "(*(i2,a2,g14.8,1x))") (i, ":", objscale(i), i=1, 4)

   end subroutine init_objective

   subroutine obj_scale()

      ! scale objective components by looping over MAXPERT iterations and
      ! tracking the how much each component changes

      integer, parameter :: MAXPERT = 1000
      real(8), allocatable :: vect_denorm(:), min_b(:), max_b(:)
      real(8), allocatable :: trial(:), trial_denorm(:)
      real(8) :: objinit(4), objdelta(4), diff, rescale
      integer :: i, j

      objscale = 1.d0
      objinit = 0.d0
      objdelta = 0.d0

      ! initial trial vector
      allocate (min_b(size(vect)), max_b(size(vect)), trial(size(vect)))
      min_b = bmin
      max_b = bmax
      vect_denorm = min_b + vect*diff

      ! the choice of the first realization here is arbitrary
      call network_forward(ysimd(:, :, 1), vect_denorm, AL)
      call indicator_transform(AL, thresholds, ndata, ncut, AL_i, ivars)

      ! initilaize a starting value for each component
      if (vario .gt. 0) then
         call obj_vario(objt_vario)
         objinit(1) = objt_vario
      end if
      if (ivario .gt. 0) then
         call obj_ivario(objt_ivario)
         objinit(2) = objt_ivario
      end if
      if (runs .gt. 0) then
         call obj_runs(objt_runs)
         objinit(3) = objt_runs
      end if
      if (npoint .gt. 0) then
         call obj_npoint(objt_npt)
         objinit(4) = objt_npt
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
            if (objt_vario .lt. 0.0) objt_vario = objinit(1)
            objdelta(1) = objdelta(1) + abs(objinit(1) - objt_vario)
         end if

         if (ivario .gt. 0) then
            call obj_ivario(objt_ivario)
            if (objt_ivario .lt. 0.0) objt_ivario = objinit(2)
            objdelta(2) = objdelta(2) + abs(objinit(2) - objt_ivario)
         end if

         if (runs .gt. 0) then
            call obj_runs(objt_runs)
            if (objt_runs .lt. 0.0) objt_runs = objinit(3)
            objdelta(3) = objdelta(3) + abs(objinit(3) - objt_runs)
         end if

         if (npoint .gt. 0) then
            call obj_npoint(objt_npt)
            if (objt_npt .lt. 0.0) objt_npt = objinit(4)
            objdelta(4) = objdelta(4) + abs(objinit(4) - objt_npt)
         end if

      end do

      ! scale objective components
      if (vario .gt. 0) objscale(1) = MAXPERT/objdelta(1)
      if (ivario .gt. 0) objscale(2) = MAXPERT/objdelta(2)
      if (runs .gt. 0) objscale(3) = MAXPERT/objdelta(3)
      if (npoint .gt. 0) objscale(4) = MAXPERT/objdelta(4)

      ! rescale factor
      rescale = 0.0
      rescale = rescale + objscale(1)*objinit(1) + objscale(2)*objinit(2) &
                + objscale(3)*objinit(3) + objscale(4)*objinit(4)
      rescale = 1.d0*max(rescale, EPSLON)

      ! user defined scaling if required
      objscale(1) = userfac(1)*objscale(1)
      objscale(2) = userfac(2)*objscale(2)
      objscale(3) = userfac(3)*objscale(3)
      objscale(4) = userfac(4)*objscale(4)

      ! final scaling
      if (vario .gt. 0) objscale(1) = objscale(1)/rescale
      if (ivario .gt. 0) objscale(2) = objscale(2)/rescale
      if (runs .gt. 0) objscale(3) = objscale(3)/rescale
      if (npoint .gt. 0) objscale(4) = objscale(4)/rescale

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
         call update_vario(heads%dirs(i), tails%dirs(i), AL, stride(:, i), expvario)
         call vario_mse(expvario, target_vario%dirs(i)%vlags, &
                        varlagdist%dirs(i)%vlags, dble(idwpow), mse)
         objt = objt + mse
      end do

      objt = objt*objscale(1)

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
                              stride(:, i), expivario)
            call vario_mse(expivario, target_ivario%cuts(j)%dirs(i)%vlags, &
                           varlagdist%dirs(i)%vlags, dble(idwpow), mse)
            objt = objt + mse
         end do
      end do

      objt = objt*objscale(2)

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
         cumruns = 0
         do j = 2, ndh
            iarr = AL_i(udhidx(j - 1) + 1:udhidx(j), i)
            call calculate_runs(iarr, maxrun, expruns)
            cumruns = cumruns + expruns
         end do
         mse = sum((dble(target_runs(:, i)) - dble(cumruns))**2)/dble(maxrun)
         objt = objt + mse
      end do

      objt = objt/ncut
      objt = objt*objscale(3)

   end subroutine obj_runs

   subroutine obj_npoint(objt)

      ! npoint connectivity component

      real(8), intent(inout) :: objt

      objt = 0.d0

   end subroutine obj_npoint

end module objective_mod
