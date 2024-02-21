module readpar_mod

   use geostat
   use makepar_mod
   use vario_mod, only: set_sill, set_rotmatrix
   use subs
   use constants
   use omp_lib, only: omp_get_max_threads

   implicit none

   integer :: lin, lout, ldbg, dbgireal, lwts, ltrg
   character(256) :: outfile
   character(256) :: dbgfile
   character(256) :: wtsfile
   character(256) :: objfile
   character(256) :: prefix

contains

   subroutine readpar()

      character(256), parameter :: parfile = 'nmropt.par'
      character(256) :: datafile
      character(256) :: poolfile
      character(256) :: runsfile
      character(256) :: npointfile
      character(256) :: omegafile
      character(256) :: str
      logical :: testfl
      integer :: test, i, j, k, ic, iv, tmp, L
      integer :: dhcol, xyzcols(3), varcol, wtcol, ncols
      real(8) :: tmin, tmax, minmax(4)
      real(8), allocatable :: tmpvar(:), tmpnsvar(:)

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
         write (*, *) "ERROR: Column for dhid must be > 0. Drill hole IDs are required for sequences."
         stop
      end if

      ! nscore flag
      read (lin, *, iostat=test) itrans
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' normal score transform flag: ', itrans

      ! trimming limits
      read (lin, *, iostat=test) tmin, tmax
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' trimming limits: ', tmin, tmax

      ! number of uncond. reals
      read (lin, *, iostat=test) nreals
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' number of unconditional realizations: ', nreals

      ! ! try ranking if less than 25 reals are specified?
      ! irank = 0
      ! nrank = nreals
      ! if (nreals .le. 25) then
      !    irank = 1
      !    nrank = 100 ! arbitrary, but should be enough
      ! end if

      ! simulation type: 0 = LU, 1 = SGS
      read (lin, *, iostat=test) stype
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' simulation type: ', stype

      ! ! grid definition
      ! read (lin, *, iostat=test) nx, xmn, xsiz
      ! if (test .ne. 0) stop "ERROR in parameter file"
      ! write (*, '(a,i0,x,2(g13.6,x))') ' nx, xmn, xsiz: ', nx, xmn, xsiz

      ! read (lin, *, iostat=test) ny, ymn, ysiz
      ! if (test .ne. 0) stop "ERROR in parameter file"
      ! write (*, '(a,i0,x,2(g13.6,x))') ' ny, ymn, ysiz: ', ny, ymn, ysiz

      ! read (lin, *, iostat=test) nz, zmn, zsiz
      ! if (test .ne. 0) stop "ERROR in parameter file"
      ! write (*, '(a,i0,x,2(g13.6,x))') ' nz, zmn, zsiz: ', nz, zmn, zsiz

      ! nxyz = nx*ny*nz

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
      write (lobj, "(A)") "NMR Objective Function"
      write (lobj, "(i1)") 8
      write (lobj, "(A)") "Iteration"
      write (lobj, "(A)") "Objective value"
      write (lobj, "(A)") "Objective value (unormalized)"
      write (lobj, "(A)") "x1"
      write (lobj, "(A)") "x2"
      write (lobj, "(A)") "Fitness history"
      write (lobj, "(A)") "Crossover prob."
      write (lobj, "(A)") "Mutation factor"

      ! prefix for target/exp values
      read (lin, '(a256)', iostat=test) prefix
      if (test .ne. 0) stop "ERROR in parameter file"
      call chknam(prefix, 256)
      write (*, "(2a)") '  output prefix: ', trim(adjustl(prefix))

      ! network architecture
      read (lin, *, iostat=test) nnet%nl
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) ' number of network layers: ', nnet%nl

      allocate (nnet%ld(nnet%nl), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      read (lin, *, iostat=test) nnet%ld
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, "(a,10(i0,x))") '  network layer dimensions: ', nnet%ld

      ! activation function
      read (lin, *, iostat=test) nnet%af
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, "(a,10(i0,x))") '  activation function: ', nnet%af

      ! regularization
      read (lin, *, iostat=test) nnet%ireg, nnet%regconst
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) '  regularization and constant: ', nnet%ireg, nnet%regconst

      ! batch normalization
      nnet%norm = .false.
      read (lin, *, iostat=test) tmp
      if (test .ne. 0) stop "ERROR in parameter file"
      if (tmp .gt. 0) nnet%norm = .true.
      write (*, *) '  normalize layer inputs?: ', nnet%norm

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

      ! number of cutoffs
      read (lin, *, iostat=test) ncut
      if (test .ne. 0) stop "ERROR in parameter file"

      ! open the output file and write headers
      open (lout, file=outfile, status="UNKNOWN")
      write (lout, "(A)") "Network Mixture"
      write (lout, "(i1)") 5 + ncut
      write (lout, "(A)") "dhid"
      write (lout, "(A)") "x"
      write (lout, "(A)") "y"
      write (lout, "(A)") "z"
      write (lout, "(A)") "NMR value"

      ! write indicator columns to output file
      do ic = 1, ncut
         write (lout, "(A, i1)") "Threshold ", ic
      end do

      ! thresholds
      allocate (thresholds(ncut), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"
      read (lin, *, iostat=test) thresholds
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) '  thresholds: ', thresholds

      ! threshold weighting
      allocate (threshwt(ncut), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"
      read (lin, *, iostat=test) threshwt
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) '  threshold wights: ', threshwt

      ! runs
      read (lin, *, iostat=test) iruns
      if (test .ne. 0) stop "ERROR in parameter file"
      if (runs .gt. 0) write (*, *) ' considering runs above thresholds: ', iruns

      read (lin, *, iostat=test) maxrun
      if (test .ne. 0) stop "ERROR in parameter file"
      if (runs .gt. 0) write (*, *) ' max number of runs: ', maxrun

      ! target runs from file?
      read (lin, *, iostat=test) t_iruns
      if (test .ne. 0) stop "ERROR in parameter file"
      read (lin, '(a256)', iostat=test) runsfile
      if (test .ne. 0) stop "ERROR in parameter file"
      call chknam(runsfile, 256)
      if (t_iruns .gt. 0) then
         inquire (file=runsfile, exist=testfl)
         if (.not. testfl) stop "ERROR - the data file does not exist"
         ! load targets if the file is correct
         open (ltrg, file=runsfile, status='OLD')
         read (ltrg, *)
         read (ltrg, *, iostat=test) ncols
         if (test .ne. 0) stop "ERROR in data file"
         if (ncols .ne. ncut) stop "Number of columns in target file must match &
&         number of thresholds!"
         ! restart target file
         rewind (ltrg)
         do i = 1, ncols + 2
            read (ltrg, *, iostat=test)
         end do
         allocate (target_runs(maxrun, ncut), stat=test)
         if (test .ne. 0) stop "allocation failed due to insufficient memory!"
         do i = 1, maxrun
            read (ltrg, *, iostat=test) target_runs(i, :)
         end do
         close (ltrg)
      end if

      ! npoint connectivity
      read (lin, *, iostat=test) inpoint
      if (test .ne. 0) stop "ERROR in parameter file"
      if (npoint .gt. 0) write (*, *) ' considering npoint connectivity &
&      above thresholds: ', inpoint

      read (lin, *, iostat=test) nstep
      if (test .ne. 0) stop "ERROR in parameter file"
      if (npoint .gt. 0) write (*, *) ' max number of connected steps: ', nstep

      ! target npoint from file?
      read (lin, *, iostat=test) t_inpoint
      if (test .ne. 0) stop "ERROR in parameter file"
      read (lin, '(a256)', iostat=test) npointfile
      if (test .ne. 0) stop "ERROR in parameter file"
      call chknam(npointfile, 256)
      if (t_inpoint .gt. 0) then
         inquire (file=npointfile, exist=testfl)
         if (.not. testfl) stop "ERROR - the data file does not exist"
         ! load targets if the file is correct
         open (ltrg, file=npointfile, status='OLD')
         read (ltrg, *)
         read (ltrg, *, iostat=test) ncols
         if (test .ne. 0) stop "ERROR in data file"
         if (ncols .ne. ncut) stop "Number of columns in target file must match &
&         number of thresholds!"
         ! restart target file
         rewind (ltrg)
         do i = 1, ncols + 2
            read (ltrg, *, iostat=test)
         end do
         allocate (target_npoint(nstep, ncut), stat=test)
         if (test .ne. 0) stop "allocation failed due to insufficient memory!"
         do i = 1, nstep
            read (ltrg, *, iostat=test) target_npoint(i, :)
         end do
         close (ltrg)
      end if

      ! differential evolution
      read (lin, *, iostat=test) mut, cplo, cphi, popsize, its
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, '(a, 3(f5.2, x), 2(i0, x))') ' mutation factor, crossover prob. low, crossover prob. &
  &    high, pop. size, iterations: ', mut, cplo, cphi, popsize, its
      if (popsize < 4) stop 'for rand1 DE popsize must be > 3'
      read (lin, *, iostat=test) bmin, bmax
      if (test .ne. 0) stop "ERROR in parameter file"
      write (*, *) 'DE lower and upper bounds: ', bmin, bmax
      ! allocate bounds arrays
      L = nnet%ld(1)
      allocate (min_b(2*L, 1), max_b(2*L, 1), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"
      min_b(1:L - 1, 1) = bmin
      max_b(1:L - 1, 1) = bmax
      ! constrain nugget lower/upper limit
      min_b(L, 1) = 0.d0
      max_b(L, 1) = 0.25d0

      ! read omega bounds
      read (lin, '(a256)', iostat=test) omegafile
      if (test .ne. 0) stop "ERROR in parameter file"
      call chknam(omegafile, 256)
      inquire (file=omegafile, exist=testfl)
      if (.not. testfl) stop "ERROR - the data file does not exist"
      ! load if the file is correct
      open (ltrg, file=omegafile, status='OLD')
      read (ltrg, *)
      read (ltrg, *, iostat=test) ncols
      if (test .ne. 0) stop "ERROR in data file"
      if (ncols .ne. 4) stop "Number of columns in omega file must be 2!"
      allocate (fprec(L), sigwt(L)) ! factor precedence
      ! restart file
      rewind (ltrg)
      do i = 1, ncols + 2
         read (ltrg, *, iostat=test)
      end do
      do i = 1, L - 1
         read (ltrg, *, iostat=test) minmax
         min_b(L + i, 1) = minmax(1)
         max_b(L + i, 1) = minmax(2)
         fprec(i) = minmax(3)
         sigwt(i) = minmax(4)
      end do
      ! constrain nugget contribution to be linear
      min_b(2*L, 1) = 1.d0
      max_b(2*L, 1) = 1.d0
      ! nugget precedence always last
      fprec(L) = L
      sigwt(L) = 1.d0
      close (ltrg)

      ! parallel processing
      ipara = 0
      read (lin, *, iostat=test) num_threads
      if (test .ne. 0) stop "ERROR in parameter file"
      if (num_threads .eq. 0) then
         num_threads = 1
         ipara = 0
      end if
      if (num_threads .eq. 1) write (*, *) "performing serial DE"
      if (num_threads .gt. 1) then
         ipara = 1
         write (*, *) "performing parallel DE with", num_threads, "threads"
      end if
      if (num_threads .lt. 0) then
         ipara = 1
         num_threads = omp_get_max_threads()
         write (*, *) "performing parallel DE with", num_threads, "threads"
      end if

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

      allocate (ivmod(ncut), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

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

      ! basic warning for potential errors
      do i = 1, nvarg
         if ((vmod(i)%sill - 1.d0 .gt. EPSLON) .and. (isill .eq. 1)) then
            write (*, *) "WARNING: standardizing exp variogram sill but &
&            variogram model sill is not 1.0"
         end if
      end do

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

      ! basic warning for potential errors
      do ic = 1, ncut
         if ((ivmod(ic)%sill - 1.d0 .gt. EPSLON) .and. (isill .eq. 1)) then
            write (*, *) "WARNING: standardizing exp variogram sill but &
&            variogram model sill is not 1.0"
         end if
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

      ! random despike
      do i = 1, ndata
         var(i) = var(i) + grnd()*SMALLDBLE
      end do

      ! nscore input var if required
      if (itrans .eq. 1) then
         call nscore(ndata, var, tmin, tmax, 1, wts, tmpnsvar, nsvar, test)
         if (test .gt. 0) stop "Error in normal score transform"
         var = nsvar
      end if

      ! build declustered cdf for quantile lookup
      vsort = var
      wsort = wts/sum(wts)
      vord = [(i, i=1, ndata)]
      call sortem(1, ndata, vsort, 2, wsort, vord, vsort, vsort, vsort, vsort, &
                  vsort, vsort)
      vtmp = dblecumsum(wsort)
      vcdf = vtmp(2:)
      vcdf = vcdf - vcdf(1)/2.0 ! midpoint

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
      write (*, *) " "

      ! allocate arrays for the pool
      ngvarg = nnet%ld(1) - 1
      allocate (pool(ngvarg), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      ! write out headers if debugging
      if (idbg .gt. 0) then
         write (ldbg, "(A)") "Debugging realizations"
         write (ldbg, "(i2)") ngvarg + 1 + 3 ! + nugget + coords
         write (ldbg, "(A)") "x"
         write (ldbg, "(A)") "y"
         write (ldbg, "(A)") "z"
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
