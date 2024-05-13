module output_mod

   use readpar_mod
   use network_mod
   use vario_mod
   use sequences_mod
   use objective_mod

   implicit none

contains

   subroutine write_files()

      ! conveience subroutine to write out final files

      character(256) :: tempfl
      real(8), allocatable :: opt_AL(:, :)
      integer, allocatable :: opt_AL_i(:, :, :)
      real(8) :: best_ivars(ncut, nreals)
      real(8), allocatable :: expvario(:), expnpoint(:), fimp(:)
      integer, allocatable :: expruns(:)
      integer :: cumruns(maxrun)
      integer :: i, ic, j, jj, k, kk, hidx, tidx

      allocate (opt_AL(ndata, nreals), opt_AL_i(ndata, ncut, nreals))

      ! write out the optimized network mixture
      call vector_to_matrices(best, nnet)
      call build_refcdf(nsamp, yref, nnet, ttable)
      do i = 1, nreals
         if (ifp) then
            call network_forward2(nnet, ysimd(:, :, i), opt_AL(:, i), .true., &
                                  fprec, sigwt, ttable)
         else
            call network_forward(nnet, ysimd(:, :, i), opt_AL(:, i), .true., ttable)
         end if
         call indicator_transform(opt_AL(:, i), thresholds, ndata, ncut, &
                                  opt_AL_i(:, :, i))
         do j = 1, ndata
            write (lout, "(*(g14.8,1x))") dhids(j), xyz(1, j), xyz(2, j), &
               xyz(3, j), opt_AL(j, i), (opt_AL_i(j, ic, i), ic=1, ncut)
         end do
      end do

      ! write out the optimized network weights
      do i = 1, size(best)
         write (lwts, "(*(g14.8,1x))") best(i)
      end do

      ! append factor precedence
      do i = 1, nnet%ld(1)
         write (lwts, "(*(g14.8,1x))") fprec(i)
      end do

      ! append sigmoid weighting
      do i = 1, nnet%ld(1)
         write (lwts, "(*(g14.8,1x))") sigwt(i)
      end do

      ! calculate and write out feature importance
      tempfl = trim(prefix)//"feature_importance.out"
      open (ltrg, file=tempfl, status="UNKNOWN")
      write (ltrg, "(A)") "Expected Feature Importance"
      write (ltrg, "(i1)") 2
      write (ltrg, "(A)") "Feature Index"
      write (ltrg, "(A)") "Feature Importance"
      call feature_importance(nnet, ysimd, fimp)
      do i = 1, nnet%ld(1)
         write (ltrg, "((i2,1x,g14.8,1x))") i, fimp(i)
      end do

      ! write out targets and experimental values

      if (vario .gt. 0) then

         ! write target variogram model and optimized values
         tempfl = trim(prefix)//"variogram.out"
         open (ltrg, file=tempfl, status="UNKNOWN")
         write (ltrg, "(A)") "Experimental Variogram Values"
         write (ltrg, "(i1)") 8
         write (ltrg, "(A)") "Realization"
         write (ltrg, "(A)") "Variogram Index"
         write (ltrg, "(A)") "Lag Distance"
         write (ltrg, "(A)") "Number of Pairs"
         write (ltrg, "(A)") "Variogram Value"
         write (ltrg, "(A)") "Target Value"
         write (ltrg, "(A)") "Azimuth"
         write (ltrg, "(A)") "Dip"
         do i = 1, nreals

            ! ! TEST STANDARDIZING BY REALIZATION
            ! call calc_expsill(opt_AL(:, i), sill, vtype=1)
            ! ! TEST STANDARDIZING BY REALIZATION

            do j = 1, ndir
               call update_vario(heads%dirs(j), tails%dirs(j), opt_AL(:, i), expvario, sill)
               ! do k = 1, size(heads%dirs(j)%lags)
               do k = 1, expvar(j)%nlags + 1
                  write (ltrg, "(*(g14.8,1x))") i, j, varlagdist%dirs(j)%vlags(k), &
                     expvar(j)%npairs(k), expvario(k), target_vario%dirs(j)%vlags(k), &
                     varazm%dirs(j)%vlags(k), vardip%dirs(j)%vlags(k)
               end do
            end do
         end do
         close (ltrg)
      end if

      ! write target indicator variogram models and optimized values
      if (ivario .gt. 0) then
         tempfl = trim(prefix)//"ivariogram.out"
         open (ltrg, file=tempfl, status="UNKNOWN")
         write (ltrg, "(A)") "Experimental Indicator Variogram Values"
         write (ltrg, "(i1)") 9
         write (ltrg, "(A)") "Realization"
         write (ltrg, "(A)") "Threshold"
         write (ltrg, "(A)") "Variogram Index"
         write (ltrg, "(A)") "Lag Distance"
         write (ltrg, "(A)") "Number of Pairs"
         write (ltrg, "(A)") "Variogram Value"
         write (ltrg, "(A)") "Target Value"
         write (ltrg, "(A)") "Azimuth"
         write (ltrg, "(A)") "Dip"
         do i = 1, nreals
            do ic = 1, ncut

               ! ! TEST STANDARDIZING BY REALIZATION
               ! call calc_expsill(opt_AL(:, i), isills(ic), vtype=2, cut=thresholds(ic))
               ! ! TEST STANDARDIZING BY REALIZATION

               do j = 1, ndir
                  call update_vario(heads%dirs(j), tails%dirs(j), dble(opt_AL_i(:, ic, i)), &
                                    expvario, isills(ic))
                  do k = 1, expvar(j)%nlags + 1
                     write (ltrg, "(*(g14.8,1x))") i, ic, j, varlagdist%dirs(j)%vlags(k), &
                        expvar(j)%npairs(k), expvario(k), target_ivario%cuts(ic)%dirs(j)%vlags(k), &
                        varazm%dirs(j)%vlags(k), vardip%dirs(j)%vlags(k)
                  end do
               end do
            end do
         end do
         close (ltrg)
      end if

      ! write out final runs values
      if (runs .gt. 0) then
         tempfl = trim(prefix)//"runs.out"
         open (ltrg, file=tempfl, status="UNKNOWN")
         write (ltrg, "(A)") "Experimental Runs Values"
         write (ltrg, "(i1)") 5
         write (ltrg, "(A)") "Realization"
         write (ltrg, "(A)") "Threshold"
         write (ltrg, "(A)") "Step"
         write (ltrg, "(A)") "Cumulative Runs"
         write (ltrg, "(A)") "Target Value"
         do i = 1, nreals
            do ic = 1, ncut
               cumruns = 0
               do j = 1, ndh
                  call binary_runs(opt_AL_i(udhidx(j) + 1:udhidx(j + 1), ic, i), &
                                   maxrun, expruns)
                  cumruns = cumruns + expruns
               end do
               do k = 1, maxrun
                  write (ltrg, "(*(g14.8,1x))") i, ic, k, cumruns(k), target_runs(k, ic)
               end do
            end do
         end do
         close (ltrg)
      end if

      ! write out final n-point connectivity values
      if (npoint .gt. 0) then
         tempfl = trim(prefix)//"npoint.out"
         open (ltrg, file=tempfl, status="UNKNOWN")
         write (ltrg, "(A)") "Experimental n-Point Connectivity Values"
         write (ltrg, "(i1)") 5
         write (ltrg, "(A)") "Realization"
         write (ltrg, "(A)") "Threshold"
         write (ltrg, "(A)") "Step"
         write (ltrg, "(A)") "Prob of Connection"
         write (ltrg, "(A)") "Target Value"
         do i = 1, nreals
            do ic = 1, ncut
               call npoint_connect(opt_AL_i(:, ic, i), nstep, ndh, udhidx, expnpoint)
               do k = 1, nstep
                  write (ltrg, "(*(g14.8,1x))") i, ic, k, expnpoint(k), target_npoint(k, ic)
               end do
            end do
         end do
         close (ltrg)
      end if

      if (idbg .gt. 0) then

         ! write out variogram pairs
         tempfl = trim(prefix)//"vario_pairs.out"
         open (lprs, file=tempfl, status="UNKNOWN")
         write (lprs, "(A)") "Experimental Variogram Pairs"
         write (lprs, "(i2)") 10
         write (lprs, "(A)") "Tail Index"
         write (lprs, "(A)") "Head Index"
         write (lprs, "(A)") "Lag Index"
         write (lprs, "(A)") "Direction Index"
         write (lprs, "(A)") "Tail x"
         write (lprs, "(A)") "Tail y"
         write (lprs, "(A)") "Tail z"
         write (lprs, "(A)") "Head x"
         write (lprs, "(A)") "Head y"
         write (lprs, "(A)") "Head z"
         do i = 1, ndir
            jj = size(heads%dirs(i)%lags)
            do j = 1, jj
               kk = size(heads%dirs(i)%lags(j)%idxs)
               if (kk .le. 1) cycle
               do k = 1, kk
                  tidx = tails%dirs(i)%lags(j)%idxs(k)
                  hidx = heads%dirs(i)%lags(j)%idxs(k)
                  write (lprs, "(*(g14.8,1x))") tidx, hidx, j, i, xyz(1, tidx), xyz(2, tidx), &
                     xyz(3, tidx), xyz(1, hidx), xyz(2, hidx), xyz(3, hidx)
               end do
            end do
         end do
         close (lprs)

         ! write out a realization if debugging, -1 for all
         if (dbgireal .gt. 0) then
            do i = 1, ndata
               write (ldbg, "(*(g14.8,1x))") xyz(1, j), xyz(2, j), xyz(3, j), &
                  (ysimd(i, j, dbgireal), j=1, ngvarg + 1)
            end do
         else if (dbgireal .lt. 0) then
            do i = 1, nreals
               do j = 1, ndata
                  write (ldbg, "(*(g14.8,1x))") xyz(1, j), xyz(2, j), xyz(3, j), &
                     (ysimd(j, k, i), k=1, ngvarg + 1)
               end do
            end do
         end if
      end if

   end subroutine write_files

end module output_mod
