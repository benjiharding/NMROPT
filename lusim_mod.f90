module lusim_mod

   use mtmod
   use covasubs
   use constants
   use subs

   implicit none

   ! data parameters
   real(8), allocatable :: xyz(:, :) ! data coordinates
   integer :: ndata

   ! Guassian pool covariance structure
   integer :: ngvarg ! number of Gaussian variograms
   real(8), allocatable :: gc0(:), gcc(:, :), gang1(:, :), gang2(:, :), &
                           gang3(:, :), gaa(:, :), ganis1(:, :), ganis2(:, :)
   integer, allocatable :: gnst(:), git(:, :)

   ! output realizations at data locations
   integer :: rseed ! random seed
   integer :: nreals ! number of realizations
   real(8), allocatable :: ysimd(:, :, :) !(ndata, ngvarg + 1, nreals)

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

      call cpu_time(start)

      ! simulate factors
      do igv = 1, ngvarg
         write (*, *) "working on factor", igv
         call lusim(igv)
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

   subroutine lusim(igv)

      ! unconditional LU simulation using variogram index igv

      integer :: MAXROT
      integer :: i, j, igv, is, ierr, ireal!, n
      real(8) :: x(ndata), y(ndata), z(ndata)
      real(8), allocatable :: rotmat(:, :, :)
      real(8) :: cmax, cova, sill
      real(8) :: C11(ndata, ndata), L11(ndata, ndata)
      real(8) :: w2(ndata, 1), y2(ndata, 1)
      real(8) :: p, xp

      ! data coordinates
      x = xyz(1, :)
      y = xyz(2, :)
      z = xyz(3, :)

      ! initialize matrices
      C11 = 0.0
      L11 = 0.0

      ! set up rotation matrix for given variogram
      MAXROT = gnst(igv)
      allocate (rotmat(MAXROT, 3, 3))
      do is = 1, gnst(igv)
         call setrot(gang1(igv, is), gang2(igv, is), gang3(igv, is), &
                     ganis1(igv, is), ganis2(igv, is), is, MAXROT, rotmat)
      end do

      ! check collocated points?

      ! get cmax and the sill
      call cova3(x(1), y(1), z(1), x(1), y(1), z(1), 1, gnst(igv), MAXGNST, &
                 gc0(igv), git(igv, :), gcc(igv, :), gaa(igv, :), 1, MAXROT, &
                 rotmat, cmax, sill)

      ! compute C11
      do i = 1, ndata
         C11(i, i) = sill
         do j = i + 1, ndata
            call cova3(x(i), y(i), z(i), x(j), y(j), z(j), 1, gnst(igv), MAXGNST, &
                       gc0(igv), git(igv, :), gcc(igv, :), gaa(igv, :), 1, MAXROT, &
                       rotmat, cmax, cova)
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

end module lusim_mod
