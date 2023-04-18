module test_suite1

   use testdrive, only: new_unittest, unittest_type, error_type, check
   use test_variograms, only: get_variogram_pairs, get_variogram_values
   use types_mod
   use readpar_mod, only: readpar, expvar, vmod, lin

   implicit none

   private

   public :: collect_suite1

   real(8), parameter :: EPSLON = 1e-6

contains

!> Collect all exported unit tests
   subroutine collect_suite1(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)

      testsuite = [ &
                  ! new_unittest("test_parse_experimental_variogram", test_parse_experimental_variogram), &
                  ! new_unittest("test_parse_variogram_model", test_parse_variogram_model), &
                  new_unittest("test_vario_pairs", test_vario_pairs), &
                  new_unittest("test_vario_values", test_vario_values) &
                  ]

   end subroutine collect_suite1

   subroutine test_parse_experimental_variogram(error)
      type(error_type), allocatable, intent(out) :: error
      type(experimental) :: true(1)
      real(8) :: diff(10)
      integer :: i

      true(1)%azm = -40
      true(1)%atol = 30
      true(1)%dip = 0
      true(1)%dtol = 30
      true(1)%bandh = 1000
      true(1)%bandv = 1000
      true(1)%tilt = 0
      true(1)%lagdis = 30
      true(1)%lagtol = 17
      true(1)%nlags = 8

      call readpar()

      diff(1) = abs(true(1)%azm - expvar(1)%azm)
      diff(2) = abs(true(1)%atol - expvar(1)%atol)
      diff(3) = abs(true(1)%dip - expvar(1)%dip)
      diff(4) = abs(true(1)%dtol - expvar(1)%dtol)
      diff(5) = abs(true(1)%bandh - expvar(1)%bandh)
      diff(6) = abs(true(1)%bandv - expvar(1)%bandv)
      diff(7) = abs(true(1)%tilt - expvar(1)%tilt)
      diff(8) = abs(true(1)%lagdis - expvar(1)%lagdis)
      diff(9) = abs(true(1)%lagtol - expvar(1)%lagtol)
      diff(10) = abs(true(1)%nlags - expvar(1)%nlags)

      do i = 1, 10
         call check(error, (diff(i) .lt. EPSLON), .true.)
         if (allocated(error)) return
      end do

   end subroutine test_parse_experimental_variogram

   subroutine test_parse_variogram_model(error)

      type(error_type), allocatable, intent(out) :: error
      type(variogram) :: true(1)
      real(8) :: diff

      true(1)%nst = 3
      true(1)%c0 = 0.1
      true(1)%sill = 1.0
      true(1)%it = [2, 2, 2]
      true(1)%cc = [0.2, 0.3, 0.4]
      true(1)%ang1 = [-40, -40, -40]
      true(1)%ang2 = [0, 0, 0]
      true(1)%ang3 = [0, 0, 0]
      true(1)%aa = [50, 175, 750]
      true(1)%ahmin = [100, 125, 125]
      true(1)%avert = [50, 175, 450]
      true(1)%anis1 = [true(1)%ahmin(1)/true(1)%aa(1), true(1)%ahmin(2)/true(1)%aa(2), &
                       true(1)%ahmin(3)/true(1)%aa(3)]
      true(1)%anis2 = [true(1)%avert(1)/true(1)%aa(1), true(1)%avert(2)/true(1)%aa(2), &
                       true(1)%avert(3)/true(1)%aa(3)]

      call readpar()

      diff = abs(true(1)%nst - vmod(1)%nst)
      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

      diff = abs(true(1)%c0 - vmod(1)%c0)
      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

      diff = abs(true(1)%sill - vmod(1)%sill)
      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

      diff = sum(abs(true(1)%it - vmod(1)%it))
      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

      diff = sum(abs(true(1)%cc - vmod(1)%cc))
      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

      diff = sum(abs(true(1)%ang1 - vmod(1)%ang1))
      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

      diff = sum(abs(true(1)%ang2 - vmod(1)%ang2))
      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

      diff = sum(abs(true(1)%ang3 - vmod(1)%ang3))
      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

      diff = sum(abs(true(1)%aa - vmod(1)%aa))
      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

      diff = sum(abs(true(1)%anis1 - vmod(1)%anis1))
      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

      diff = sum(abs(true(1)%anis2 - vmod(1)%anis2))
      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

      diff = sum(abs(true(1)%ahmin - vmod(1)%ahmin))
      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

      diff = sum(abs(true(1)%avert - vmod(1)%avert))
      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

   end subroutine test_parse_variogram_model

   subroutine test_vario_pairs(error)
      type(error_type), allocatable, intent(out) :: error
      integer, allocatable :: pairs(:, :)
      real(8), allocatable :: lagbins(:)
      real(8) :: diff
      integer :: true(2, 2)
      integer :: i

      true(1, :) = [1, 5]
      true(2, :) = [4, 5]

      call get_variogram_pairs(pairs, lagbins)

      diff = abs(sqrt(2.0)*0.5 - lagbins(2))

      do i = 1, size(true, dim=1)
         call check(error, pairs(i, 1), true(i, 1))
         call check(error, pairs(i, 2), true(i, 2))
         call check(error, pairs(i, 3), 2)
         call check(error, (diff .lt. EPSLON), .true.)
         if (allocated(error)) return
      end do

   end subroutine test_vario_pairs

   subroutine test_vario_values(error)
      type(error_type), allocatable, intent(out) :: error
      real(8), allocatable :: expvario(:)
      real(8) :: true, diff

      true = ((1.0 - 1.0)**2 + (0.5 - 1.0)**2)/2*0.5

      call get_variogram_values(expvario)

      diff = abs(true - expvario(1))

      call check(error, (diff .lt. EPSLON), .true.)
      if (allocated(error)) return

   end subroutine test_vario_values

end module test_suite1
