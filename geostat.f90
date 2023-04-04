module geostat

   use types_mod

   implicit none

   ! GLOBAL VARIABLE DECLARATIONS

   !
   ! objective module
   !

   ! candidate vectors
   real(8), allocatable :: AL(:) ! Gaussian mixture
   integer, allocatable :: AL_i(:, :) ! indicator transform of AL (nd, ncut)

   ! objective function targets
   type(vario_array) :: target_vario ! (ndir*nlags)
   type(ivario_array) :: target_ivario ! (ncut, ndir*nlags)

   integer, allocatable :: target_runs(:, :) ! (ncut, maxrun)
   real(8), allocatable :: target_npoint(:, :) ! (ncut, nstep)

   ! objective function scaling
   real(8) :: objscale(4) ! fobj scaling parameters
   real(8) :: objt_vario, objt_ivario, objt_runs, objt_npt ! temp obj values

   ! drillhole parameters
   integer, allocatable :: dhids(:), dhlens(:) ! dhids and length
   integer, allocatable :: udhids(:), udhidx(:) ! unique dhids
   integer :: ndh ! number of drillholes
   real(8), allocatable :: x(:), y(:), z(:) ! coordinates
   integer, allocatable :: iz(:, :) ! indicator transform of 'var'

   ! variogram parameters
   integer, allocatable :: pairs(:, :) ! (npairs*ndir, 3)
   integer, allocatable :: head(:), tail(:), lag(:), dir(:)
   integer :: max_pairs
   type(lag_array) :: heads, tails
   type(vario_array) :: varlagdist
   type(vario_array) :: varazm, vardip

   real(8), allocatable :: thresholds(:) ! for indicator transform
   real(8), allocatable :: ivars(:) ! indicator sills (ncut)
   real(8) :: sill
   integer, allocatable :: nlags(:) ! multiple directions
   type(experimental), allocatable :: expvar(:)
   type(variogram), allocatable :: vmod(:)
   type(variogram), allocatable :: ivmod(:)
   real(8) :: idwpow
   integer :: nvarg, nivarg, isill, ndir, ncut
   integer, allocatable :: udiridx(:), ulagidx(:) ! array indices
   integer :: e, f, g, h, nc

   ! run parameters
   integer :: maxrun, nstep ! maximum runs and connected steps

   ! optimization initialization
   real(8) :: bmin, bmax
   real(8) :: userfac(4)

   ! boolean flags
   integer :: vario
   integer :: ivario
   integer :: runs
   integer :: npoint
   integer :: runs_above
   integer :: conn_above

   ! output file
   integer :: idbg, lprs = 7 ! cant be 6

   !
   ! network module
   !

   ! network architecture
   integer, allocatable :: layer_dims(:), iwts(:), ibias(:)
   real(8), allocatable :: vect(:) ! candidate vector
   integer :: af ! activation function

   ! data parameters
   real(8), allocatable :: var(:), wts(:) ! variable and weights

   !
   ! simulation module
   !

   ! data parameters
   real(8), allocatable :: xyz(:, :) ! data coordinates
   integer :: ndata

   ! Guassian pool covariance structure
   integer :: ngvarg ! number of Gaussian variograms
   type(variogram), allocatable :: pool(:)

   ! output realizations at data locations
   integer :: rseed ! random seed
   integer :: nreals ! number of realizations
   integer :: stype ! simulation type
   real(8), allocatable :: ysimd(:, :, :) !(ndata, ngvarg + 1, nreals)

   !
   ! DE module
   !
   real(8) :: mut, cplo, cphi, crossp
   integer :: popsize, its
   real(8), allocatable :: best(:)

   ! objective iterations output file
   integer :: lobj = 5

   !
   ! sequences module
   !
   integer :: iruns, inpoint

end module geostat
