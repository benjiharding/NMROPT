module test_subroutines

   use subs
   use mtmod

   implicit none

contains

   subroutine test_nscore(mu, sigma, gmin, gmax, iwt, seed)

      integer, parameter :: nd = 10000
      real(8), parameter :: tmin = -1e21, tmax = 1e21
      integer, intent(in) :: iwt, seed
      real(8), intent(inout) :: mu, sigma, gmin, gmax
      real(8) :: vr(nd), wt(nd), sumsqs
      real(8), allocatable :: vrg(:), tmp(:)
      integer :: i, ierror

      ! initialize a uniform random vector
      call sgrnd(seed)
      do i = 1, nd
         vr(i) = grnd()
         if (iwt .eq. 1) then
            wt(i) = grnd()
         else if (iwt .eq. 0) then
            wt(i) = 1.d0
         end if
      end do

      ! shift the mean
      vr = vr + 5

      ! add some extremes
      vr(ubound(vr)) = 50
      vr(lbound(vr)) = -10

      call nscore(nd, vr, tmin, tmax, iwt, wt, tmp, vrg, ierror)
      gmin = minval(vrg)
      gmax = maxval(vrg)

      mu = 0.d0
      sumsqs = 0.d0

      do i = 1, nd
         mu = mu + vrg(i)
         sumsqs = sumsqs + vrg(i)*vrg(i)
      end do

      mu = mu/nd
      sumsqs = sumsqs/nd
      sigma = sqrt(sumsqs - mu*mu)

   end subroutine test_nscore

end module test_subroutines
