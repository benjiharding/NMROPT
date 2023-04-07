program test_network

   implicit none

   real(8), allocatable :: x(:), Y(:, :), yy(:)
   integer, allocatable :: connections(:), biases(:), iwts(:), ibias(:)
   integer, allocatable :: layer_dims(:)
   real(8), allocatable :: Amat(:, :), A_prev(:, :), AL(:)
   real(8), allocatable :: W(:, :), WL(:, :), b(:, :), bL(:, :), &
                           Zmat(:, :), ZL(:, :)
   integer :: sw(2), sb(2)
   integer :: i, nwts, nbias, ld

   x = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1]
   yy = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, &
         16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29]
   Y = reshape(yy, [10, 3])
   ! call print_matrix(Y)

   connections = [6, 2]
   biases = [2, 1]
   layer_dims = [3, 2, 1]

   nwts = 8
   nbias = 3
   ld = size(layer_dims)

   allocate (ibias(nbias))

   iwts = cumsum(connections)
   ibias(1) = nwts
   ibias(2:) = biases
   ibias = cumsum(ibias)
   ibias = ibias(2:)

   print *, iwts
   print *, ibias

   Amat = Y

   do i = 2, ld - 1

      sw = [layer_dims(i), layer_dims(i - 1)] ! wts
      sb = [layer_dims(i), 1] ! bias vector

      print *, sw
      print *, sb

      ! reshape weights and biases
      W = reshape(x(iwts(i - 1) + 1:iwts(i)), shape=(sw), order=[2, 1])
      b = reshape(x(ibias(i - 1) + 1:ibias(i)), shape=(sb), order=[2, 1])

      print *, ""
      call print_matrix(W)
      print *, ""
      call print_matrix(b)

      ! transpose prior to forward pass
      W = transpose(W)
      b = transpose(b)

      print *, ""
      call print_matrix(W)
      print *, ""
      call print_matrix(b)

      A_prev = Amat
      b = spread(b(1, :), 1, size(A_prev, dim=1))
      Zmat = matmul(A_prev, W) + b
      print *, ""
      call print_matrix(Zmat)

      Amat = Zmat

   end do

   ! output layer
   sw = [layer_dims(ld), layer_dims(ld - 1)] ! wt matrix
   sb = [layer_dims(ld), 1] ! bias vector

   print *, sw
   print *, sb

   WL = reshape(x(iwts(ld - 1) + 1:iwts(ld)), shape=(sw), order=[2, 1])
   WL = transpose(WL)

   bL = reshape(x(ibias(ld - 1) + 1:ibias(ld)), shape=(sb), order=[2, 1])
   bL = transpose(bL)

   print *, ""
   call print_matrix(WL)
   print *, ""
   call print_matrix(bL)

   bL = spread(bL(1, :), 1, size(Amat, dim=1))
   ZL = matmul(Amat, WL) + bL

   print *, ""
   call print_matrix(ZL)

   ! linear activation and reduce dims
   AL = ZL(:, 1)

contains

   subroutine print_matrix(A)
      real(8), intent(in) :: A(:, :)  ! An assumed-shape dummy argument

      integer :: i

      do i = 1, size(A, 1)
         print *, A(i, :)
      end do

   end subroutine print_matrix

   function cumsum(a) result(b)
      integer, intent(in) :: a(:)
      integer :: b(size(a) + 1)
      integer :: i
      b(1) = 0
      b(2:) = a
      do i = 2, size(a) + 1
         b(i) = b(i - 1) + a(i - 1)
      end do
   end function cumsum

end program test_network
