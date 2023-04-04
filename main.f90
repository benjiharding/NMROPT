program main

   use readpar_mod
   use sim_mod
   use network_mod
   use objective_mod
   use de_mod

   implicit none

   ! read the parfile
   call readpar

   ! initilize random generator
   call sgrnd(rseed)

   ! simulate at the data locations
   call simulate

   ! initialize network parameters
   call init_network

   ! initialize objective function parameters
   call init_objective

   ! minimize the objective
   call optimize

   ! write files
   call write_files

   close (lout)
   close (ldbg)
   close (lwts)
   close (lobj)

contains

   subroutine write_files()

      ! conveience subroutine to write out final files

      real(8), allocatable :: opt_AL(:)
      integer, allocatable :: opt_AL_i(:, :)
      real(8), allocatable :: expvario(:), expnpoint(:)
      integer, allocatable :: expruns(:)
      integer :: cumruns(maxrun)
      integer :: i, ic, j, k

      allocate (opt_AL(ndata), opt_AL_i(ndata, ncut))

      ! write out the optimized network mixture
      do i = 1, nreals
         call network_forward(ysimd(:, :, i), best, opt_AL)
         do j = 1, ndata
            write (lout, "(*(g14.8,1x))") opt_AL(j)
         end do
      end do

      ! write out the optimized network weights
      do i = 1, size(best)
         write (lwts, "(*(g14.8,1x))") best(i)
      end do

      if (idbg .gt. 0) then

         if (vario .gt. 0) then
            ! write target variogram model and optimized values
            open (ltrg, file="variogram.out", status="UNKNOWN")
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
               call network_forward(ysimd(:, :, i), best, opt_AL)
               call calc_expsill(opt_AL, sill)
               do j = 1, ndir
                  call update_vario(heads%dirs(j), tails%dirs(j), opt_AL, expvario, sill)
                  do k = 1, size(heads%dirs(j)%lags)
                     write (ltrg, "(*(g14.8,1x))") i, j, varlagdist%dirs(j)%vlags(k), &
                        size(heads%dirs(j)%lags(k)%idxs), expvario(k), target_vario%dirs(j)%vlags(k), &
                        varazm%dirs(j)%vlags(k), vardip%dirs(j)%vlags(k)
                  end do
               end do
            end do
            close (ltrg)
         end if

         if (ivario .gt. 0) then
            ! write target indicator variogram models and optimized values
            open (ltrg, file="ivariogram.out", status="UNKNOWN")
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
               call network_forward(ysimd(:, :, i), best, opt_AL)
               call indicator_transform(opt_AL, thresholds, ndata, ncut, opt_AL_i, ivars)
               do ic = 1, ncut
                  do j = 1, ndir
                     call update_vario(heads%dirs(j), tails%dirs(j), dble(opt_AL_i(:, ic)), expvario, 1.d0)
                     do k = 1, size(heads%dirs(j)%lags)
                        write (ltrg, "(*(g14.8,1x))") i, ic, j, varlagdist%dirs(j)%vlags(k), &
                           size(heads%dirs(j)%lags(k)%idxs), expvario(k)/ivmod(ic)%sill, &
                           target_ivario%cuts(ic)%dirs(j)%vlags(k)/ivmod(ic)%sill, varazm%dirs(j)%vlags(k), &
                           vardip%dirs(j)%vlags(k)
                     end do
                  end do
               end do
            end do
            close (ltrg)
         end if

         if (runs .gt. 0) then
            ! write out final runs values
            open (ltrg, file="runs.out", status="UNKNOWN")
            write (ltrg, "(A)") "Experimental Runs Values"
            write (ltrg, "(i1)") 5
            write (ltrg, "(A)") "Realization"
            write (ltrg, "(A)") "Threshold"
            write (ltrg, "(A)") "Step"
            write (ltrg, "(A)") "Cumulative Runs"
            write (ltrg, "(A)") "Target Value"
            do i = 1, nreals
               call network_forward(ysimd(:, :, i), best, opt_AL)
               call indicator_transform(opt_AL, thresholds, ndata, ncut, opt_AL_i, ivars)
               do ic = 1, ncut
                  cumruns = 0
                  do j = 1, ndh
                     call binary_runs(opt_AL_i(udhidx(j) + 1:udhidx(j + 1), ic), &
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

         if (npoint .gt. 0) then
            ! write out final n-point connectivity values
            open (ltrg, file="npoint.out", status="UNKNOWN")
            write (ltrg, "(A)") "Experimental n-Point Connectivity Values"
            write (ltrg, "(i1)") 5
            write (ltrg, "(A)") "Realization"
            write (ltrg, "(A)") "Threshold"
            write (ltrg, "(A)") "Step"
            write (ltrg, "(A)") "Prob of Connection"
            write (ltrg, "(A)") "Target Value"
            do i = 1, nreals
               call network_forward(ysimd(:, :, i), best, opt_AL)
               call indicator_transform(opt_AL, thresholds, ndata, ncut, opt_AL_i, ivars)
               do ic = 1, ncut
                  call npoint_connect(opt_AL_i(:, ic), nstep, ndh, udhidx, expnpoint)
                  do k = 1, nstep
                     write (ltrg, "(*(g14.8,1x))") i, ic, k, expnpoint(k), target_npoint(k, ic)
                  end do
               end do
            end do
            close (ltrg)
         end if

         ! write out a realization if debugging, -1 for all
         if (dbgireal .gt. 0) then
            do i = 1, ndata
               write (ldbg, "(*(g14.8,1x))") (ysimd(i, j, dbgireal), j=1, ngvarg + 1)
            end do
         else if (dbgireal .lt. 0) then
            do i = 1, nreals
               do j = 1, ndata
                  write (ldbg, "(*(g14.8,1x))") (ysimd(j, k, i), k=1, ngvarg + 1)
               end do
            end do
         end if

         ! end debug IF
      end if

   end subroutine write_files

end program main
