module readpar_mod

   use geostat
   use makepar_mod
   use vario_mod, only: set_sill, set_rotmatrix
   use subs
   use constants

   implicit none

   integer :: lin, lout, ldbg, dbgireal, lwts, ltrg
   character(256) :: outfile
   character(256) :: dbgfile
   character(256) :: wtsfile
   character(256) :: objfile

contains

   subroutine readpar()

      character(256), parameter :: parfile = 'nmropt.par'
      character(256) :: datafile
      character(256) :: poolfile
      character(256) :: str
      logical :: testfl
      integer :: test, i, j, k, ic, iv
      integer :: dhcol, xyzcols(3), varcol, wtcol, ncols
      real(8) :: tmin, tmax
      integer :: nx, ny, nz
      real(8) :: xmn, ymn, zmn
      real(8) :: xsiz, ysiz, zsiz
      integer :: nxyz
      integer :: nnl
      real(8), allocatable :: tmpvar(:)

      ! unit numbers
      lin = 1
      lout = 2
      ldbg = 3
      lwts = 4
      ltrg = 8

      ! attempt to open the parfile
      do i = 1, 256
         str(i:i) = ' '
      end do
      call getarg(1, str)

      ! if the input argument is empty request a parameter file name
      if (str .eq. "") then
         write (*, "('Which parameter file do you want to use?')")
         read (*, '(a)') str
      end if
      str = trim(adjustl(str))
      if (str .eq. "") str = parfile

      ! verify if the parameter file exists
      inquire (file=str, exist=testfl)

      ! create a blank parameter file if required
      if (.not. testfl) then
         print *, "ERROR - the parameter file does not exist"
         print *
         if (str .eq. parfile) then
            print *, "        creating a blank parameter file"
            call makepar(parfile)
            print *
         end if
         stop
      end if

      ! open the parfile
      open (lin, file=str, status="old")
      read (lin, '(a4)', iostat=test) str(1:4)
      if (test .ne. 0) stop "ERROR in parameter file"

      ! find the start
      do while (str(1:4) .ne. 'STAR')
         read (lin, '(a4)', iostat=test) str(1:4)
         if (test .ne. 0) stop "ERROR in parameter file"
      end do

      ! read input data file
      write (*, *) " "
      write (*, *) " reading parameter file..."
      write (*, *) " "

      read (lin, '(a256)', iostat=test) datafile
      if (test .ne. 0) stop "ERROR in parameter file"
      call chknam(datafile, 256)
      inquire (file=datafile, exist=testfl)
      if (.not. testfl) stop "ERROR - the data file does not exist"
      write (*, "(2a)") '  data file: ', trim(adjustl(datafile))

      ! read columns for coordinates and values
      read (lin, *, iostat=test) dhcol, xyzcols, varcol, wtcol
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, "(a,6(i0,x))") '  column for dhid, x, y, z, var and wt: ', dhcol, xyzcols, varcol, wtcol
      if (dhcol .le. 0) then
         write (*, *) "ERROR: Column for dhid must be > 0. Drill hole IDs are required for computations."
         stop
      end if

      ! ! nscore flag
      ! read (lin, *, iostat=test) itrans
      ! if (test .ne. 0) stop "ERROR in parameter file"
      ! write (*, *) ' normal score transform flag: ', itrans

      ! trimming limits
      read (lin, *, iostat=test) tmin, tmax
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' trimming limits: ', tmin, tmax

      ! number of uncond. reals
      read (lin, *, iostat=test) nreals
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' number of unconditional realizations: ', nreals

      ! simulation type: 0 = LU, 1 = SGS
      read (lin, *, iostat=test) stype
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' simulation type: ', stype

      ! grid definition
      read (lin, *, iostat=test) nx, xmn, xsiz
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, '(a,i0,x,2(g13.6,x))') ' nx, xmn, xsiz: ', nx, xmn, xsiz

      read (lin, *, iostat=test) ny, ymn, ysiz
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, '(a,i0,x,2(g13.6,x))') ' ny, ymn, ysiz: ', ny, ymn, ysiz

      read (lin, *, iostat=test) nz, zmn, zsiz
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, '(a,i0,x,2(g13.6,x))') ' nz, zmn, zsiz: ', nz, zmn, zsiz

      nxyz = nx*ny*nz

      ! random number seed
      read (lin, *, iostat=test) rseed
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' random number seed: ', rseed

      ! debugging information and output
      read (lin, *, iostat=test) idbg, dbgireal
      if (test .ne. 0) stop "ERROR in parameter file"
      if (dbgireal .gt. nreals) stop "Debugging realization index cannot &
&       be greater than the number of realizations"

      read (lin, '(a256)', iostat=test) dbgfile
      if (test .ne. 0) stop "ERROR in parameter file"
      call chknam(dbgfile, 256)
      write (*, "(2a)") '  debugging file: ', trim(adjustl(dbgfile))

      if (idbg .gt. 0) then
         open (ldbg, file=dbgfile, status="UNKNOWN")
      end if

      ! network output
      read (lin, '(a256)', iostat=test) outfile
      if (test .ne. 0) stop "ERROR in parameter file"
      call chknam(outfile, 256)
      write (*, "(2a)") '  output file: ', trim(adjustl(outfile))

      ! open the output file and write headers
      open (lout, file=outfile, status="UNKNOWN")
      write (lout, "(a15)") "Network Mixture"
      write (lout, "(i1)") 1
      write (lout, "(a9)") "nmr value"

      ! network weights
      read (lin, '(a256)', iostat=test) wtsfile
      if (test .ne. 0) stop "ERROR in parameter file"
      call chknam(wtsfile, 256)
      write (*, "(2a)") '  output file: ', trim(adjustl(wtsfile))

      ! open the weight file and write headers
      open (lwts, file=wtsfile, status="UNKNOWN")
      write (lwts, "(a15)") "Network Weights"
      write (lwts, "(i1)") 1
      write (lwts, "(a6)") "weight"

      ! obj function vs iterations
      read (lin, '(a256)', iostat=test) objfile
      if (test .ne. 0) stop "ERROR in parameter file"
      call chknam(objfile, 256)
      write (*, "(2a)") '  output file: ', trim(adjustl(objfile))

      ! open the objective file and write headers
      open (lobj, file=objfile, status="UNKNOWN")
      write (lobj, "(a22)") "NMR Objective Function"
      write (lobj, "(i1)") 2
      write (lobj, "(a9)") "iteration"
      write (lobj, "(a15)") "objective value"

      ! network architecture
      read (lin, *, iostat=test) nnet%nl !nnl
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' number of network layers: ', nnet%nl !nnl

      ! allocate (layer_dims(nnl), stat=test)
      allocate (nnet%ld(nnet%nl), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      read (lin, *, iostat=test) nnet%ld !layer_dims
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, "(a,10(i0,x))") '  network layer dimensions: ', nnet%ld !layer_dims

      ! activation function
      read (lin, *, iostat=test) nnet%af !af
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, "(a,10(i0,x))") '  activation function: ', nnet%af !af

      ! Gaussian pool file
      read (lin, '(a256)', iostat=test) poolfile
      if (test .ne. 0) stop "ERROR in parameter file"
      call chknam(poolfile, 256)
      inquire (file=poolfile, exist=testfl)
      if (.not. testfl) stop "ERROR - the data file does not exist"
      write (*, "(2a)") '  file with Gaussian pool cov. struct.: ', &
         trim(adjustl(poolfile))

      ! objective components
      read (lin, *, iostat=test) vario, ivario, runs, npoint
      if (test .ne. 0) stop "ERROR in parameter file"

      read (lin, *, iostat=test) userfac
      if (test .ne. 0) stop "ERROR in parameter file"

      if (vario .gt. 0) write (*, *) ' considering a variogram'
      if (vario .gt. 0) write (*, *) '   user scaling factor = ', userfac(1)
      if (ivario .gt. 0) write (*, *) ' considering indicator variogram(s)'
      if (ivario .gt. 0) write (*, *) '   user scaling factor = ', userfac(2)
      if (runs .gt. 0) write (*, *) ' considering runs'
      if (runs .gt. 0) write (*, *) '   user scaling factor = ', userfac(3)
      if (npoint .gt. 0) write (*, *) ' considering npoint connectivity'
      if (npoint .gt. 0) write (*, *) '   user scaling factor = ', userfac(4)

      ! runs
      read (lin, *, iostat=test) iruns
      if (test .ne. 0) stop "ERROR in parameter file"
      if (runs .gt. 0) write (*, *) ' considering runs above thresholds: ', iruns

      read (lin, *, iostat=test) maxrun
      if (test .ne. 0) stop "ERROR in parameter file"
      if (runs .gt. 0) write (*, *) ' max number of runs: ', maxrun

      ! npoint connectivity
      read (lin, *, iostat=test) inpoint
      if (test .ne. 0) stop "ERROR in parameter file"
      if (npoint .gt. 0) write (*, *) ' considering npoint connectivity &
&      above thresholds: ', inpoint

      read (lin, *, iostat=test) nstep
      if (test .ne. 0) stop "ERROR in parameter file"
      if (npoint .gt. 0) write (*, *) ' max number of connected steps: ', nstep

      ! differential evolution
      read (lin, *, iostat=test) mut, cplo, cphi, popsize, its
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, '(a, 3(f5.2, x), 2(i0, x))') ' mutation factor, crossover prob. low, crossover prob. &
  &    high, pop. size, iterations: ', mut, cplo, cphi, popsize, its
      if (popsize < 4) stop 'for rand1 DE popsize must be > 3'
      read (lin, *, iostat=test) bmin, bmax
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) 'DE lower and upper bounds: ', bmin, bmax

      ! exp variogram parameters
      read (lin, *, iostat=test) ndir
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) 'number of experimental directions: ', ndir

      allocate (expvar(ndir), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      do i = 1, ndir
         read (lin, *, iostat=test) expvar(i)%azm, expvar(i)%atol, expvar(i)%bandh, &
            expvar(i)%dip, expvar(i)%dtol, expvar(i)%bandv, expvar(i)%tilt
         if (test .ne. 0) stop "ERROR in parameter file"
         write (*, *) '  azm, azmtol, bandhorz', expvar(i)%azm, expvar(i)%atol, expvar(i)%bandh
         write (*, *) '  dip, diptol, bandvert', expvar(i)%dip, expvar(i)%dtol, expvar(i)%bandv
         write (*, *) '  tilt', expvar(i)%tilt
         read (lin, *, iostat=test) expvar(i)%nlags, expvar(i)%lagdis, expvar(i)%lagtol
         if (test .ne. 0) stop "ERROR in parameter file"
         write (*, *) '  nlags, lagdist, lagtol', expvar(i)%nlags, expvar(i)%lagdis, expvar(i)%lagtol
      end do

      ! number of variograms
      read (lin, *, iostat=test) nvarg
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) '  number of variogram models: ', nvarg

      allocate (vmod(nvarg), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      ! number of indicator variograms
      read (lin, *, iostat=test) nivarg
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) '  number of indicator variogram models: ', nivarg
      ncut = nivarg

      allocate (ivmod(ncut), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      ! thresholds
      allocate (thresholds(ncut), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"
      read (lin, *, iostat=test) thresholds
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) '  thresholds: ', thresholds

      ! IDW power for weighting
      read (lin, *, iostat=test) idwpow
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) '  IDW power: ', idwpow

      ! max number of experimental pairs
      read (lin, *, iostat=test) max_pairs
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) '  Max number of experimental pairs: ', max_pairs

      ! standardize sills?
      read (lin, *, iostat=test) isill
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) '  flag to standardize sills: ', isill

      ! parse continous variogram model(s)
      do i = 1, nvarg
         read (lin, *, iostat=test) vmod(i)%nst, vmod(i)%c0
         if (test .ne. 0) stop "ERROR in parameter file"
         write (*, *) '  nst, c0: ', vmod(i)%nst, vmod(i)%c0
         if (vmod(i)%nst .gt. MAXNST) then
            write (*, *) 'nst must be less than or equal to ', MAXNST
            stop
         end if

         allocate (vmod(i)%it(vmod(i)%nst), vmod(i)%cc(vmod(i)%nst), &
                   vmod(i)%ang1(vmod(i)%nst), vmod(i)%ang2(vmod(i)%nst), vmod(i)%ang3(vmod(i)%nst), &
                   vmod(i)%aa(vmod(i)%nst), vmod(i)%anis1(vmod(i)%nst), vmod(i)%anis2(vmod(i)%nst), &
                   vmod(i)%ahmin(vmod(i)%nst), vmod(i)%avert(vmod(i)%nst), stat=test)
         if (test .ne. 0) stop "allocation failed due to insufficient memory!"

         do j = 1, vmod(i)%nst
            read (lin, *, iostat=test) vmod(i)%it(j), vmod(i)%cc(j), vmod(i)%ang1(j), vmod(i)%ang2(j), vmod(i)%ang3(j)
            if (test .ne. 0) stop "ERROR in parameter file"
            read (lin, *, iostat=test) vmod(i)%aa(j), vmod(i)%ahmin(j), vmod(i)%avert(j)
            if (test .ne. 0) stop "ERROR in parameter file"
            vmod(i)%anis1(j) = vmod(i)%ahmin(j)/max(vmod(i)%aa(j), EPSLON)
            vmod(i)%anis2(j) = vmod(i)%avert(j)/max(vmod(i)%aa(j), EPSLON)
            write (*, *) ' it, cc, ang[1,2,3]; ', vmod(i)%it(j), vmod(i)%cc(j), vmod(i)%ang1(j), &
               vmod(i)%ang2(j), vmod(i)%ang3(j)
            write (*, *) ' a1 a2 a3: ', vmod(i)%aa(j), vmod(i)%ahmin(j), vmod(i)%avert(j)
         end do
      end do

      call set_sill(vmod)
      call set_rotmatrix(vmod)

      ! parse indicator variogram models
      if (ivario .gt. 0) then
         do ic = 1, ncut
            read (lin, *, iostat=test) ivmod(ic)%nst, ivmod(ic)%c0
            if (test .ne. 0) stop "ERROR in parameter file"
            write (*, *) '  inst, ic0: ', ivmod(ic)%nst, ivmod(ic)%c0
            if (ivmod(ic)%nst .gt. MAXNST) then
               write (*, *) 'inst must be less than or equal to ', MAXNST
               stop
            end if

            allocate (ivmod(ic)%it(ivmod(ic)%nst), ivmod(ic)%cc(ivmod(ic)%nst), &
                      ivmod(ic)%ang1(ivmod(ic)%nst), ivmod(ic)%ang2(ivmod(ic)%nst), &
                      ivmod(ic)%ang3(ivmod(ic)%nst), ivmod(ic)%aa(ivmod(ic)%nst), &
                      ivmod(ic)%anis1(ivmod(ic)%nst), ivmod(ic)%anis2(ivmod(ic)%nst), &
                      ivmod(ic)%ahmin(ivmod(ic)%nst), ivmod(ic)%avert(ivmod(ic)%nst), stat=test)
            if (test .ne. 0) stop "allocation failed due to insufficient memory!"

            do j = 1, ivmod(ic)%nst
               read (lin, *, iostat=test) ivmod(ic)%it(j), ivmod(ic)%cc(j), ivmod(ic)%ang1(j), &
                  ivmod(ic)%ang2(j), ivmod(ic)%ang3(j)
               if (test .ne. 0) stop "ERROR in parameter file"
               read (lin, *, iostat=test) ivmod(ic)%aa(j), ivmod(ic)%ahmin(j), ivmod(ic)%avert(j)
               if (test .ne. 0) stop "ERROR in parameter file"
               ivmod(ic)%anis1(j) = ivmod(ic)%ahmin(j)/max(ivmod(ic)%aa(j), EPSLON)
               ivmod(ic)%anis2(j) = ivmod(ic)%avert(j)/max(ivmod(ic)%aa(j), EPSLON)
               write (*, *) ' iit, icc, iang[1,2,3]; ', ivmod(ic)%it(j), ivmod(ic)%cc(j), &
                  ivmod(ic)%ang1(j), ivmod(ic)%ang2(j), ivmod(ic)%ang3(j)
               write (*, *) ' a1 a2 a3: ', ivmod(ic)%aa(j), ivmod(ic)%ahmin(j), ivmod(ic)%avert(j)
            end do
         end do

         call set_sill(ivmod)
         call set_rotmatrix(ivmod)

      end if

      ! finished reading parameters
      close (lin)

      ! start reading the data file
      write (*, *) " "
      write (*, *) " reading data file..."
      write (*, *) " "

      open (lin, file=datafile, status='OLD')
      read (lin, *)
      read (lin, *, iostat=test) ncols
      if (test .ne. 0) stop "ERROR in data file"

      ! check column numbers
      if (dhcol .gt. ncols .or. &
          any(xyzcols .gt. ncols) .or. &
          (varcol .gt. ncols) .or. &
          (wtcol .gt. ncols)) then
         write (*, *) 'there are only ', ncols, ' columns in the data file'
         write (*, *) '  your specification is out of range'
         stop
      end if
      allocate (tmpvar(ncols), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      ! jump headers and get names
      do i = 1, ncols
         read (lin, *, iostat=test) str
         if (test .ne. 0) stop "ERROR in data file"
      end do

      ! get the number of data
      ndata = 0
      do
         read (lin, *, iostat=test) tmpvar(:)
         if (test > 0) stop "ERROR in data file"
         if (test < 0) exit
         ndata = ndata + 1
      end do

      ! allocate arrays for input data
      allocate (dhids(ndata), xyz(3, ndata), var(ndata), wts(ndata), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      ! restart data file
      rewind (lin)
      do i = 1, ncols + 2
         read (lin, *, iostat=test)
      end do

      ! read file again, but storing variable and weight
      ndata = 0
      do
         read (lin, *, iostat=test) tmpvar(:)
         if (test > 0) stop "ERROR in data file"
         if (test < 0) exit
         ndata = ndata + 1
         dhids(ndata) = tmpvar(dhcol)
         xyz(:, ndata) = tmpvar(xyzcols)
         var(ndata) = tmpvar(varcol)
         wts(ndata) = tmpvar(wtcol)
      end do

      ! assume equal weighting if not specified
      if (wtcol .eq. 0) then
         wts = 1.d0
      end if

      ! get unique drillhole ids - this array is sorted
      udhids = unique(dhids)

      ! get drillhole lengths
      ndh = size(udhids)
      allocate (dhlens(ndh))
      do i = 1, ndh
         j = 0
         do k = 1, ndata
            if (dhids(k) .eq. udhids(i)) then
               j = j + 1
            end if
            dhlens(i) = j
         end do
      end do

      ! start reading Gaussian pool file
      write (*, *) " "
      write (*, *) " reading covariance structure of Gaussian pool..."

      ! allocate arrays for the pool
      ngvarg = nnet%ld(1) - 1 !layer_dims(1) - 1
      allocate (pool(ngvarg), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      ! write out headers if debugging
      if (idbg .gt. 0) then
         write (ldbg, "(a22)") "Debugging realizations"
         write (ldbg, "(i2)") ngvarg + 1
         do iv = 1, ngvarg + 1
            write (ldbg, "(a6, i3)") "Factor", iv
         end do
      end if

      ! open the pool file
      open (lin, file=poolfile, status='OLD')

      ! parse the Gaussian variogram models
      do iv = 1, ngvarg
         read (lin, *, iostat=test) pool(iv)%nst, pool(iv)%c0
         if (test .ne. 0) stop "ERROR in parameter file"
         write (*, *) '  gnst, gc0: ', pool(iv)%nst, pool(iv)%c0
         if (pool(iv)%nst .gt. MAXGNST) then
            write (*, *) 'gnst must be equal to ', MAXGNST
            stop
         end if

         allocate (pool(iv)%it(pool(iv)%nst), pool(iv)%cc(pool(iv)%nst), &
                   pool(iv)%ang1(pool(iv)%nst), pool(iv)%ang2(pool(iv)%nst), pool(iv)%ang3(pool(iv)%nst), &
                   pool(iv)%aa(pool(iv)%nst), pool(iv)%anis1(pool(iv)%nst), pool(iv)%anis2(pool(iv)%nst), &
                   pool(iv)%ahmin(pool(iv)%nst), pool(iv)%avert(pool(iv)%nst), stat=test)
         if (test .ne. 0) stop "allocation failed due to insufficient memory!"

         do j = 1, pool(iv)%nst
            read (lin, *, iostat=test) pool(iv)%it(j), pool(iv)%cc(j), pool(iv)%ang1(j), &
               pool(iv)%ang2(j), pool(iv)%ang3(j)
            if (test .ne. 0) stop "ERROR in parameter file"
            read (lin, *, iostat=test) pool(iv)%aa(j), pool(iv)%ahmin(j), pool(iv)%avert(j)
            if (test .ne. 0) stop "ERROR in parameter file"
            pool(iv)%anis1(j) = pool(iv)%ahmin(j)/max(pool(iv)%aa(j), EPSLON)
            pool(iv)%anis2(j) = pool(iv)%avert(j)/max(pool(iv)%aa(j), EPSLON)
            write (*, *) ' git, gcc, gang[1,2,3]; ', pool(iv)%it(j), pool(iv)%cc(j), &
               pool(iv)%ang1(j), pool(iv)%ang2(j), pool(iv)%ang3(j)
            write (*, *) ' a1 a2 a3: ', pool(iv)%aa(j), pool(iv)%ahmin(j), pool(iv)%avert(j)
         end do
      end do

      call set_sill(pool)
      call set_rotmatrix(pool)

      ! finished reading data
      close (lin)

   end subroutine readpar

end module readpar_mod
