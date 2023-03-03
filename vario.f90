program vario_test
   use vario_mod
   implicit none

   real(8), allocatable :: lags(:), wts(:)
   real(8), allocatable :: zval(:), zc(:), ivars(:)
   real(8), allocatable :: x(:), y(:), z(:)
   real(8), allocatable :: lagbins(:), expvario(:)
   integer, allocatable :: iz(:, :), pairs(:, :)
   integer :: np, nlags, ncut, ndata, nvargs
   integer :: nst
   real(8) :: c0
   integer, allocatable :: it(:)
   real(8), dimension(:), allocatable :: cc, azm, dip, tilt, ahmax, &
                                         ahmin, avert, varlagdist, &
                                         varazm, vardip, varmodelvals
   real(8) :: power, exp_azm, atol, bandh, exp_dip, dtol, &
              bandv, lagdis, lagtol, idwpow, sill, mse

   lags = [2.5, 5.0]
   nlags = size(lags)
   power = 2.5
   wts = inv_dist(lags, power, nlags)

   zval = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
   zc = [1.0, 2.0, 3.0]
   ndata = size(zval)
   ncut = size(zc)

   ! returns iz, ivars
   call indicator_transform(zval, zc, iz, ivars)

   x = [0.5, 1.5, 2.5, 0.5, 1.5, 2.5, 0.5, 1.5, 2.5]
   y = [0.5, 0.5, 0.5, 1.5, 1.5, 1.5, 2.5, 2.5, 2.5]
   z = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
   exp_azm = 45.0
   atol = 5.0
   bandh = 5.0
   exp_dip = 0.0
   dtol = 15.0
   bandv = 5.0
   lagdis = 2.5
   lagtol = 1.5

   ! returns pairs, lagbins
   call vario_pairs(x, y, z, exp_azm, atol, bandh, exp_dip, dtol, bandv, &
                    nlags, lagdis, lagtol, ndata, pairs, lagbins)

   np = size(pairs, 1)
!    do i = 1, np
!       print "(3(i5))", pairs(i, :)
!    end do
!    print *, "lag bins: ", lagbins

   ! returns expvario
   call update_vario(pairs, zval, nlags, expvario)
!    print *, "exp vario points: ", expvario

   ! vario model parameters to calc values
   nst = 3
   c0 = 0.1
   it = [2, 2, 2]
   cc = [0.2, 0.3, 0.4]
   azm = [-40.0, -40.0, -40.0]
   dip = [0.0, 0.0, 0.0]
   tilt = [0.0, 0.0, 0.0]
   ahmax = [50.0, 175.0, 750.0]
   ahmin = [100.0, 125.0, 125.0]
   avert = [50.0, 175.0, 450.0]
   varlagdist = [11.25987346006488, 34.416810340114736, 62.56047473550351, &
                 90.99785761540689, 121.11408939721576, 150.45973333661067, &
                 179.75717048749755, 209.85258862661385]
   nvargs = size(varlagdist)
   varazm = [-40.0, -40.0, -40.0, -40.0, -40.0, -40.0, -40.0, -40.0]
   vardip = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
   allocate (varmodelvals(nvargs))

   idwpow = 1.0
   sill = 1.0
   expvario = [0.299994, 0.4340728, 0.50034643, 0.5869547, 0.63732723, &
               0.75148828, 0.81216737, 0.78014661]

   ! returns varmodelvals
   call varmodelpts(nst, c0, it, cc, azm, dip, tilt, ahmax, ahmin, &
                    avert, nvargs, varlagdist, varazm, vardip, &
                    varmodelvals)
   ! returns mse
   call vario_mse(expvario, varmodelvals, varlagdist, idwpow, sill, mse)

   print *, varmodelvals, mse

end program vario_test

