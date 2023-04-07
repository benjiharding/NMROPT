module network_mod

   use geostat, only: layer_dims, iwts, ibias, vect, af, wts
   use mtmod
   use subs
   use constants

   implicit none

   ! interface for procedural pointer
   abstract interface
      function afunc(y) result(a)
         real(8), intent(in) :: y(:, :)
         real(8) :: a(size(y, 1), size(y, 2))
      end function afunc
   end interface

contains

   subroutine init_network()

      ! initialize network layer parameters and weight vector

      ! local variables
      integer, allocatable :: nwts(:), nbias(:)
      integer :: i, j, ld

      ld = size(layer_dims)
      allocate (nwts(ld - 1), nbias(ld - 1))
      allocate (iwts(ld), ibias(ld))

      ! number of connections and bias terms
      do i = 1, ld - 1
         nwts(i) = layer_dims(i)*layer_dims(i + 1)
      end do

      do j = 2, ld
         nbias(j - 1) = layer_dims(j)
      end do

      ! matrix indices from cumulative sums
      iwts(1) = 0
      iwts(2:) = nwts
      do i = 2, ld
         iwts(i) = iwts(i - 1) + nwts(i - 1)
      end do

      ibias(1) = sum(nwts)
      ibias(2:) = nbias
      do j = 2, ld
         ibias(j) = ibias(j - 1) + nbias(j - 1)
      end do

      ! allocate output array
      allocate (vect(sum(nwts) + sum(nbias)))

      do i = 1, size(vect)
         vect(i) = grnd()
      end do

   end subroutine init_network

   subroutine network_forward(Ymat, v, AL)
      ! forward pass through network

      ! iterate over layers(2:)
      ! get number of connection in each layer n^l * n^l-1
      ! get cumulative sum of connecitons to get 1D ids, ie for each matrix
      ! do the same for the bias terms
      ! on each iteration reshape the array slice to (n^l * n^l-1)
      ! do matrix math and apply activations
      ! linear activation on final layer followed by nscore

      ! parameters
      real(8), intent(in) :: Ymat(:, :), v(:)

      ! return
      real(8), intent(inout) :: AL(:)

      ! internal variables
      procedure(afunc), pointer :: f_ptr => null()
      real(8), allocatable :: Amat(:, :), A_prev(:, :)
      real(8), allocatable :: W(:, :), WL(:, :), b(:, :), bL(:, :), &
                              Zmat(:, :), ZL(:, :)
      real(8), allocatable :: vrg(:), tmp(:)
      integer :: sw(2), sb(2)
      integer :: i, ld, ierr

      Amat = Ymat
      ld = size(layer_dims)

      ! pointer to activation function
      if (af .eq. 1) then
         f_ptr => sigmoid
      else if (af .eq. 2) then
         f_ptr => hyptan
      else if (af .eq. 3) then
         f_ptr => relu
      else if (af .eq. 4) then
         f_ptr => linear
      end if

      ! hidden layers
      do i = 2, ld - 1

         ! matrix shapes
         sw = [layer_dims(i), layer_dims(i - 1)] ! wts
         sb = [layer_dims(i), 1] ! bias vector

         ! reshape weights and biases
         W = reshape(v(iwts(i - 1) + 1:iwts(i)), shape=(sw), order=[2, 1])
         b = reshape(v(ibias(i - 1) + 1:ibias(i)), shape=(sb), order=[2, 1])

         ! transpose prior to forward pass
         W = transpose(W)
         b = transpose(b)

         ! forward pass
         A_prev = Amat
         b = spread(b(1, :), 1, size(A_prev, dim=1))
         Zmat = matmul(A_prev, W) + b
         Amat = f_ptr(Zmat)

      end do

      ! output layer
      sw = [layer_dims(ld), layer_dims(ld - 1)] ! wt matrix
      sb = [layer_dims(ld), 1] ! bias vector

      WL = reshape(v(iwts(ld - 1) + 1:iwts(ld)), shape=(sw), order=[2, 1])
      WL = transpose(WL)

      bL = reshape(v(ibias(ld - 1) + 1:ibias(ld)), shape=(sb), order=[2, 1])
      bL = transpose(bL)
      bL = spread(bL(1, :), 1, size(Amat, dim=1))

      ZL = matmul(Amat, WL) + bL

      ! linear activation and reduce dims
      AL = ZL(:, 1)

      ! normal score transform
      ! call nscore(AL, wts, vrg)
      call nscore(size(AL), AL, dble(-1.0e21), dble(1.0e21), 1, wts, tmp, vrg, ierr)
      AL = vrg

   end subroutine network_forward

   ! subroutine nscore(dvar, dwts, vrg)

   !    ! parameters
   !    real(8), intent(in) :: dvar(:), dwts(:)

   !    ! result
   !    real(8), allocatable, intent(out) :: vrg(:)

   !    ! internal variables
   !    real(8), allocatable :: vr(:)  ! temp copy of 'dvar'
   !    real(8), allocatable :: wt_ns(:)
   !    real(8) :: vrrg, vrr
   !    real(8) :: u, twt, wtfac, oldcp, cp, w
   !    real(8), dimension(1) :: c, d, e, f, g, h, aa
   !    integer :: i, j, nd, ierr

   !    nd = size(dvar, dim=1)
   !    wt_ns = dwts

   !    ! allocate output array
   !    allocate (vr(nd), vrg(nd))

   !    ! add small random value for sorting
   !    do i = 1, nd
   !       u = grnd()
   !       vr(i) = dvar(i) + u*EPSLON
   !       twt = twt + wt_ns(i)
   !    end do

   !    ! sort by data value
   !    call sortem(1, nd, vr, 1, wt_ns, c, d, e, f, g, h, aa)

   !    ! compute the cumulative probabilities
   !    wtfac = 1.d0/twt
   !    oldcp = 0.d0
   !    cp = 0.d0
   !    do j = 1, nd
   !       w = wtfac*wt_ns(j)
   !       cp = cp + w
   !       wt_ns(j) = (cp + oldcp)/2.d0
   !       call gauinv(wt_ns(j), vrrg, ierr)
   !       oldcp = cp
   !       ! now, reset the weight to the normal scores value
   !       wt_ns(j) = vrrg
   !    end do

   !    ! normal scores transform
   !    do i = 1, nd
   !       u = grnd()
   !       vrr = dvar(i) + u*EPSLON
   !       call locate(vr, nd, 1, nd, vrr, j)
   !       j = min(max(1, j), (nd - 1))
   !       vrg(i) = powint(vr(j), vr(j + 1), wt_ns(j), wt_ns(j + 1), vrr, 1.d0)
   !    end do

   ! end subroutine nscore

   subroutine print_matrix(A)
      real(8), intent(in) :: A(:, :)  ! An assumed-shape dummy argument

      integer :: i

      do i = 1, size(A, 1)
         print *, A(i, :)
      end do

   end subroutine print_matrix

   function relu(yval) result(a)

      ! rectified lienar unit activation

      real(8), intent(in) :: yval(:, :)
      real(8) :: a(size(yval, 1), size(yval, 2))

      a = max(0.d0, yval)

   end function relu

   function sigmoid(yval) result(a)

      ! sigmoid activation

      real(8), intent(in) :: yval(:, :)
      real(8) :: a(size(yval, 1), size(yval, 2))

      a = 1.d0/(1.d0 + exp(-yval))

   end function sigmoid

   function hyptan(yval) result(a)

      ! hyperbolic tangent activation

      real(8), intent(in) :: yval(:, :)
      real(8) :: a(size(yval, 1), size(yval, 2))

      a = dble(tanh(yval))

   end function hyptan

   function linear(yval) result(a)

      ! linear activation

      real(8), intent(in) :: yval(:, :)
      real(8) :: a(size(yval, 1), size(yval, 2))

      a = yval

   end function linear

end module network_mod
