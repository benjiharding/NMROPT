module readpar_mod

   use makepar_mod
   use network_mod
   use lusim_mod
   use objective_mod
   use sequences_mod
   use de_mod
   use subs
   use constants

   implicit none

   integer :: lin, lout, ldbg, dbgireal, lwts
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
      read (lin, *, iostat=test) nnl
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' number of network layers: ', nnl

      allocate (layer_dims(nnl), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      read (lin, *, iostat=test) layer_dims
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, "(a,10(i0,x))") '  network layer dimensions: ', layer_dims

      ! activation function
      read (lin, *, iostat=test) af
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, "(a,10(i0,x))") '  activation function: ', af

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

      allocate (azm(ndir), atol(ndir), bandh(ndir), dip(ndir), dtol(ndir), &
                bandv(ndir), nlags(ndir), lagdis(ndir), lagtol(ndir), &
                tilt(ndir), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      do i = 1, ndir
         read (lin, *, iostat=test) azm(i), atol(i), bandh(i), dip(i), dtol(i), &
            bandv(i), tilt(i)
         if (test .ne. 0) stop "ERROR in parameter file"
         write (*, *) '  azm, azmtol, bandhorz', azm(i), atol(i), bandh(i)
         write (*, *) '  dip, diptol, bandvert', dip(i), dtol(i), bandv(i)
         write (*, *) '  tilt', tilt(i)
         read (lin, *, iostat=test) nlags(i), lagdis(i), lagtol(i)
         if (test .ne. 0) stop "ERROR in parameter file"
         write (*, *) '  nlags, lagdist, lagtol', nlags(i), lagdis(i), lagtol(i)
      end do

      ! number of variograms
      read (lin, *, iostat=test) nvarg
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) '  number of variogram models: ', nvarg

      allocate (nst(MAXNST), it(MAXNST), c0(MAXNST), cc(MAXNST), ang1(MAXNST), &
                ang2(MAXNST), ang3(MAXNST), aa(MAXNST), anis1(MAXNST), &
                anis2(MAXNST), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      ! number of indicator variograms
      read (lin, *, iostat=test) nivarg
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) '  number of indicator variogram models: ', nivarg
      ncut = nivarg

      allocate (inst(ncut), iit(ncut, MAXNST), ic0(ncut), &
                icc(ncut, MAXNST), iang1(ncut, MAXNST), iang2(ncut, MAXNST), &
                iang3(ncut, MAXNST), iaa(ncut, MAXNST), ianis1(ncut, MAXNST), &
                ianis2(ncut, MAXNST), stat=test)
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

      ! standardize sills?
      read (lin, *, iostat=test) isill
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) '  flag to standardize sills: ', isill

      ! parse continous variogram model(s)
      do i = 1, nvarg
         read (lin, *, iostat=test) nst(i), c0(i)
         if (test .ne. 0) stop "ERROR in parameter file"
         write (*, *) '  nst, c0: ', nst(i), c0(i)
         if (nst(i) .gt. MAXNST) then
            write (*, *) 'nst must be less than or equal to ', MAXNST
            stop
         end if
         do j = 1, nst(i)
            read (lin, *, iostat=test) it(j), cc(j), ang1(j), ang2(j), ang3(j)
            if (test .ne. 0) stop "ERROR in parameter file"
            read (lin, *, iostat=test) aa(j), aa1, aa2
            if (test .ne. 0) stop "ERROR in parameter file"
            anis1(j) = aa1/max(aa(j), EPSLON)
            anis2(j) = aa2/max(aa(j), EPSLON)
            write (*, *) ' it, cc, ang[1,2,3]; ', it(j), cc(j), ang1(j), &
               ang2(j), ang3(j)
            write (*, *) ' a1 a2 a3: ', aa(j), aa1, aa2
         end do
      end do

      ! parse indicator variogram models
      if (ivario .gt. 0) then
         do ic = 1, ncut
            read (lin, *, iostat=test) inst(ic), ic0(ic)
            if (test .ne. 0) stop "ERROR in parameter file"
            write (*, *) '  inst, ic0: ', inst(ic), ic0(ic)
            if (inst(ic) .gt. MAXNST) then
               write (*, *) 'inst must be less than or equal to ', MAXNST
               stop
            end if
            do j = 1, inst(ic)
               read (lin, *, iostat=test) iit(ic, j), icc(ic, j), iang1(ic, j), &
                  iang2(ic, j), iang3(ic, j)
               if (test .ne. 0) stop "ERROR in parameter file"
               read (lin, *, iostat=test) iaa(ic, j), aa1, aa2
               if (test .ne. 0) stop "ERROR in parameter file"
               ianis1(ic, j) = aa1/max(iaa(ic, j), EPSLON)
               ianis2(ic, j) = aa2/max(iaa(ic, j), EPSLON)
               write (*, *) ' iit, icc, iang[1,2,3]; ', iit(ic, j), icc(ic, j), &
                  iang1(ic, j), iang2(ic, j), iang3(ic, j)
               write (*, *) ' a1 a2 a3: ', iaa(ic, j), aa1, aa2
            end do
         end do
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
      ngvarg = layer_dims(1) - 1
      allocate (gnst(ngvarg), git(ngvarg, MAXGNST), gc0(ngvarg), &
                gcc(ngvarg, MAXGNST), gang1(ngvarg, MAXGNST), &
                gang2(ngvarg, MAXGNST), gang3(ngvarg, MAXGNST), &
                gaa(ngvarg, MAXGNST), ganis1(ngvarg, MAXGNST), &
                ganis2(ngvarg, MAXGNST), stat=test)
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
         read (lin, *, iostat=test) gnst(iv), gc0(iv)
         if (test .ne. 0) stop "ERROR in parameter file"
         write (*, *) '  gnst, gc0: ', gnst(iv), gc0(iv)
         if (gnst(iv) .gt. MAXGNST) then
            write (*, *) 'gnst must be equal to ', MAXGNST
            stop
         end if
         do j = 1, gnst(iv)
            read (lin, *, iostat=test) git(iv, j), gcc(iv, j), gang1(iv, j), &
               gang2(iv, j), gang3(iv, j)
            if (test .ne. 0) stop "ERROR in parameter file"
            read (lin, *, iostat=test) gaa(iv, j), aa1, aa2
            if (test .ne. 0) stop "ERROR in parameter file"
            ganis1(iv, j) = aa1/max(gaa(iv, j), EPSLON)
            ganis2(iv, j) = aa2/max(gaa(iv, j), EPSLON)
            write (*, *) ' git, gcc, gang[1,2,3]; ', git(iv, j), gcc(iv, j), &
               gang1(iv, j), gang2(iv, j), gang3(iv, j)
            write (*, *) ' a1 a2 a3: ', gaa(iv, j), aa1, aa2
         end do
      end do

      ! finished reading data
      close (lin)

   end subroutine readpar

end module readpar_mod
