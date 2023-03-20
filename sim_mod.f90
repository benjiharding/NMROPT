module sim_mod

   use geostat
   use mtmod
   use kdtree2_module
   use types_mod
   use vario_mod, only: set_sill, set_rotmatrix
   use covasubs
   use constants
   use subs

   implicit none

contains

   subroutine simulate()

      real(8) :: c0(ndata)
      integer :: i, ireal, igv, test, ierr
      real(8) :: p, xp
      real(8) :: start, finish

      ! allocte array for realizations
      allocate (ysimd(ndata, ngvarg + 1, nreals), stat=test)
      if (test .ne. 0) stop "allocation failed due to insufficient memory!"

      write (*, *) " "
      write (*, *) "Simulating unconditional factors..."

      if (stype .eq. 0) then
         write (*, *) "using LU simulation"
      else
         write (*, *) "using sequential Gaussian simulation"
      end if

      call cpu_time(start)

      ! simulate factors
      do igv = 1, ngvarg
         write (*, *) "working on factor", igv
         if (stype .eq. 0) then
            call lusim(igv)
         else
            call sequential_sim(igv)
         end if
      end do

      ! nugget effect
      do ireal = 1, nreals

         do i = 1, ndata
            c0(i) = grnd()
         end do

         do i = 1, ndata
            p = dble(c0(i))
            call gauinv(p, xp, ierr)
            c0(i) = xp
         end do

         ! last factor is nugget
         ysimd(:, ngvarg + 1, ireal) = c0

      end do

      call cpu_time(finish)

      write (*, *) " "
      print '("simulation took ", f5.2, " minutes")', (finish - start)/60

   end subroutine simulate

   subroutine sequential_sim(igv)

      ! unconditional sequential Guassian simulation at data locations
      ! using variogram index igv

      ! kdtree required variables
      type(kdtree2), pointer :: tree
      type(kdtree2_result), allocatable :: results(:)
      integer :: nfound

      ! normal equations
      real(8), allocatable :: rhs(:), lhs(:, :), kwts(:)
      real(8) :: cmean, cstdev
      integer, allocatable :: nuse(:), useidx(:, :)

      ! simulation
      real(8), allocatable :: sim(:)
      integer, allocatable :: isim(:), randpath(:)
      real(8) :: p, xp
      integer :: simidx, ierr
      integer, parameter :: nsearch = 40

      ! indexes
      integer :: i, j, k, igv, ireal, test

      ! allocate arrays based on number of data and max search
      allocate (rhs(nsearch), lhs(nsearch, nsearch), kwts(nsearch))
      allocate (nuse(ndata), useidx(ndata, ndata))
      allocate (sim(ndata), isim(ndata), randpath(ndata))
      allocate (results(ndata))

      ! build the search tree
      tree => kdtree2_create(input_data=xyz, dim=3, sort=.true., rearrange=.true.)

      ! main loop over realizations
      do ireal = 1, nreals

         ! define random path through nodes
         randpath = [(i, i=1, ndata)]
         call shuffle(randpath)

         ! all location are initially unsimulated
         isim = 0
         useidx = 0
         nuse = 0

         ! loop over locations to simulate
         do i = 1, ndata

            ! what location are we at?
            simidx = randpath(i)

            ! query tree for this location
            call kdtree2_r_nearest(tp=tree, qv=xyz(:, simidx), r2=pool(igv)%aa(1)**2, &
                                   nfound=nfound, nalloc=ndata, results=results)

            ! loop over samples found in search
            nuse(i) = 0
            do j = 1, nfound

               ! check if this data index is the simulation index
               if (results(j)%idx .eq. simidx) cycle

               ! check if this data index is simulated or not
               if (isim(results(j)%idx) .eq. 0) cycle ! no conditioning value here

               ! meet minimum covariance? (not collocated)
               rhs(1) = get_cov(pool(igv), xyz(:, simidx), xyz(:, results(j)%idx))
               if (rhs(1) .lt. MINCOV) cycle

               ! if we got this far increment number found at ith location
               nuse(i) = nuse(i) + 1
               ! track data indices found at ith location
               useidx(nuse(i), i) = results(j)%idx

               ! have we met the max search?
               if (nuse(i) .ge. nsearch) exit

            end do

            ! calculate matrices for normal equations
            if (nuse(i) .gt. 0) then
               do j = 1, nuse(i)
                  ! build rhs vector
                  rhs(j) = get_cov(pool(igv), xyz(:, simidx), xyz(:, useidx(j, i)))
                  do k = j, nuse(i)
                     ! diagonal
                     if (j .eq. k) then
                        lhs(j, j) = 1.d0
                     else
                        ! build lhs matrix
                        lhs(j, k) = get_cov(pool(igv), xyz(:, useidx(k, i)), xyz(:, useidx(j, i)))
                        if (lhs(j, k) .lt. 0) stop 'ERROR: Negative covariance.'
                        lhs(k, j) = lhs(j, k)
                     end if
                  end do
               end do

               ! solve the kriging system - external call to LAPACK
               call solve(lhs(1:nuse(i), 1:nuse(i)), kwts(1:nuse(i)), rhs(1:nuse(i)), &
                          nuse(i), 1, test)

               ! calcualte conditional mean and stdev
               cmean = 0.d0
               cstdev = get_cov(pool(igv), [0.d0, 0.d0, 0.d0], [0.d0, 0.d0, 0.d0]) ! variance
               do j = 1, nuse(i)
                  ! cmean considers previusly simulated values
                  cmean = cmean + kwts(j)*sim(useidx(j, i))
                  cstdev = cstdev - kwts(j)*rhs(j)
               end do

            else
               ! if no data the distribution is N(0,1)
               cmean = 0.d0
               cstdev = 1.d0

            end if

            ! draw a random number and simulate
            p = grnd()
            call gauinv(p, xp, ierr)
            sim(simidx) = xp*cstdev + cmean

            ! update this location as simulated
            isim(simidx) = 1

         end do

         ysimd(:, igv, ireal) = sim

      end do

   end subroutine sequential_sim

   subroutine lusim(igv)

      ! unconditional LU simulation using Gaussian variogram index igv

      integer :: i, j, igv, ierr, ireal
      real(8) :: cova
      real(8) :: C11(ndata, ndata), L11(ndata, ndata)
      real(8) :: w2(ndata, 1), y2(ndata, 1)
      real(8) :: p, xp

      ! initialize matrices
      C11 = 0.0
      L11 = 0.0

      ! compute C11
      do i = 1, ndata
         C11(i, i) = pool(igv)%sill
         do j = i + 1, ndata
            cova = get_cov(pool(igv), xyz(:, i), xyz(:, j))
            C11(i, j) = cova
            C11(j, i) = cova
         end do
      end do

      ! compute L11
      call chol(C11, L11, ndata, ndata, ierr)

      ! simulate realizations for this variogram model
      do ireal = 1, nreals

         do i = 1, ndata
            w2(i, 1) = grnd()
         end do

         do i = 1, ndata
            p = dble(w2(i, 1))
            call gauinv(p, xp, ierr)
            w2(i, 1) = xp
         end do

         y2 = matmul(L11, w2)
         ysimd(:, igv, ireal) = y2(:, 1)

      end do

   end subroutine lusim

   subroutine chol(a, t, n, ndim, ierr)
      !-----------------------------------------------------------------------
      !
      !                      Cholesky Decomposition
      !                      **********************
      !
      ! This subroutine calculates the lower triangular matrix T which, when
      ! multiplied by its own transpose, gives the symmetric matrix A. (from
      ! "Numerical Analysis of Symmetric Matrices,"  H.R. Schwarz et al.,
      ! p. 254)
      !
      !
      !
      ! INPUT VARIABLES:
      !
      !   a(n,n)           Symmetric positive definite matrix to be
      !                      decomposed (destroyed in the calculation of t)
      !   t(n,n)           Lower triangular matrix solution
      !   n                Dimension of the system you're decomposing
      !   ndim             Dimension of matrix a (Note: In the main program,
      !                      matrix a may have been dimensioned larger than
      !                      necessary, i.e. n, the size of the system you're
      !                      decomposing, may be smaller than ndim.)
      !   ierr             Error code:  ierr=0 - no errors; ierr=1 - matrix a
      !                      is not positive definite
      !
      !
      !
      ! NO EXTERNAL REFERENCES:
      !-----------------------------------------------------------------------
      real(8) :: a(ndim, ndim), t(ndim, ndim)
      integer :: n, ndim, ierr
      integer :: i, k, ip

      ierr = 0
      !
      ! Check for positive definiteness:
      !
      do ip = 1, n
         if (a(ip, ip) .le. 0.0) then
            write (*, '(a)') 'WARNING: chol - not positive definite'
            ierr = 1
            go to 1
         end if
         t(ip, ip) = sqrt(a(ip, ip))
         if (ip .ge. n) return
         do k = ip + 1, n
            t(k, ip) = a(ip, k)/t(ip, ip)
         end do
         do i = ip + 1, n
            do k = i, n
               a(i, k) = a(i, k) - t(i, ip)*t(k, ip)
            end do
         end do
1        continue
      end do

   end subroutine chol

   subroutine linv(a, b, n, ndim)
      !-----------------------------------------------------------------------
      !
      !                Inverse of a Lower Triangular Matrix
      !                ************************************
      !
      ! This subroutine finds the inverse of a lower triangular matrix A and
      ! stores the answer in B. (from "Numerical Analysis of Symmetric
      ! Matrices,"  H.R. Schwarz et al.,)
      !
      !
      !
      ! INPUT VARIABLES:
      !
      !   a(n,n)           Lower triangular matrix to be inverted
      !   b(n,n)           the inverse
      !   n                Dimension of the matrix you're inverting
      !   ndim             Dimension of matrix a (Note: In the main program,
      !                      matrix a may have been dimensioned larger than
      !                      necessary, i.e. n, the size of the system you're
      !                      decomposing, may be smaller than ndim.)
      !
      !-----------------------------------------------------------------------
      real(8) :: a(ndim, ndim), b(ndim, ndim), sum_
      integer :: n, ndim
      integer :: i, k, j

      do i = 1, n
         if (i .gt. 1) then
            do k = 1, i - 1
               sum_ = 0.
               do j = k, i - 1
                  sum_ = sum_ + a(i, j)*b(j, k)
               end do
               b(i, k) = -sum_/a(i, i)
            end do
         end if
         b(i, i) = 1./a(i, i)
      end do

   end subroutine linv

end module sim_mod
