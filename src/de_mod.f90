module de_mod

   use geostat, only: mut, cplo, cphi, crossp, popsize, its, best, &
                      bmin, bmax, lobj, nnet
   use objective_mod, only: obj_nmr
   use types_mod
   use mtmod
   use subs
   use constants

   implicit none

contains

   subroutine optimize()

      real(8) :: start, finish

      write (*, *) " "
      write (*, *) "Starting differential evolution..."
      write (*, *) " "

      call cpu_time(start)

      call de(nnet%dims, popsize, its, mut, cplo, cphi, bmin, bmax, best)

      call cpu_time(finish)

      write (*, *) " "
      print '("Optimization took ", f5.2, " minutes")', (finish - start)/60

   end subroutine optimize

   subroutine de(dims, popsize, its, mut, cplo, cphi, bmin, bmax, best, ifunc)

      ! differential evolution

      integer, parameter :: maxnochange = 200

      ! inputs
      integer, intent(in) :: dims, popsize, its
      real(8), intent(in) :: mut, cplo, cphi, bmin, bmax
      real(8), allocatable, intent(inout) :: best(:)
      integer, optional, intent(in) :: ifunc

      ! local variables
      integer :: i, j, k, best_idx, idxs(3), func, nochange
      real(8), allocatable :: pop(:, :), pop_denorm(:, :)
      real(8), allocatable :: min_b(:, :), max_b(:, :), diff(:, :)
      real(8), allocatable :: mutant(:)
      real(8), allocatable :: trial(:), trial_denorm(:), cr(:)
      real(8), allocatable :: fitness(:), fobj(:)
      real(8) :: fit, ftry

      ! interface for passing arbitrary objective functions
      if (.not. present(ifunc)) then
         func = 0
      else
         func = ifunc
      end if

      ! allocate the population and bounds
      allocate (best(dims))
      allocate (pop(dims, popsize))
      allocate (min_b(dims, 1), max_b(dims, 1))
      min_b = bmin
      max_b = bmax

      ! allocate the trial vectors
      allocate (trial(dims))
      allocate (cr(dims))

      ! allocate objective counter
      allocate (fobj(its))
      fobj = huge(fit)

      ! initialize random population
      do i = 1, dims
         do j = 1, popsize
            pop(i, j) = grnd()
         end do
      end do

      ! calculate the denormalized population
      diff = abs(min_b - max_b)
      pop_denorm = spread(min_b(:, 1), 2, popsize) + &
                   pop*spread(diff(:, 1), 2, popsize)

      ! evaluate the fitness of each member
      allocate (fitness(popsize))
      fitness = 0.d0
      do i = 1, popsize
         ! call obj_nmr(pop_denorm(:, i), fit)
         ! fitness(i) = fit
         fitness(i) = objfunc(pop_denorm(:, i), func)
      end do

      ! get the current best
      best_idx = minloc(fitness, dim=1)
      best = pop_denorm(:, best_idx)

      nochange = 0

      ! main loop
      do i = 1, its

         nochange = nochange + 1

         ! write out the current value?
         if (modulo(i, 100) .eq. 0) then
            write (*, *) " working on DE iteration", i
            if (i .gt. 1) then
               write (*, *) " current objective value", fobj(i - 1)/fobj(1)
            end if
         end if

         ! crossover prob b/w cpho and cphi
         ! crossp = grnd()*(cphi - cplo) + cplo
         crossp = cphi + (cplo - cphi)*(1 - i/its)**4

         ! loop over population and evolve
         do j = 1, popsize

            ! sample without replacement
            idxs = random_sample(popsize, j)

            ! mutuate the selected vectors
            mutant = rand1_mutation(idxs, dims, pop, mut)

            ! crossover and get the trial vector
            trial = crossover(mutant, crossp, pop(:, j), dims)

            ! denormalize the trial and evaluate
            ftry = 0.d0
            trial_denorm = min_b(:, 1) + trial*diff(:, 1)
            ! call obj_nmr(trial_denorm, ftry)
            ftry = objfunc(trial_denorm, func)

            ! did we improve?
            if (ftry .lt. fitness(j)) then
               fitness(j) = ftry
               pop(:, j) = trial

               ! do we need to update the global best?
               if (ftry .lt. fitness(best_idx)) then
                  best_idx = j
                  best = trial_denorm
                  nochange = 0
               end if
            end if

         end do

         ! store this iterations best value
         fobj(i) = fitness(best_idx)

         ! have we reached the nochange threshold?
         if (nochange .gt. maxnochange) then

            write (*, *) " random restart on iteration", i

            ! randomly mutate the population
            do k = 1, popsize

               ! exclude the current best from mutation
               if (k .eq. best_idx) cycle

               idxs = random_sample(popsize, k)
               mutant = rand1_mutation(idxs, dims, pop, mut)
               pop(:, k) = crossover(mutant, crossp, pop(:, k), dims)

            end do

            nochange = 0

         end if

         if (.not. present(ifunc)) then
            ! write out the values
            write (lobj, "(1(i0,1x),1(g14.8,1x))") i, fobj(i)/fobj(1)
         end if

         ! are we done early?
         if (fobj(i) .le. 1e-10) return

      end do

   end subroutine de

   function random_sample(popsz, target) result(idxs)

      integer, intent(in) :: popsz, target
      integer :: idxs(3)
      real(8) :: rand(3)
      integer :: i

      rand = 0.d0
      do i = 1, 3
         rand(i) = grnd()
      end do
      idxs = floor(rand*popsz + 1.d0)

      do while (idxs(1) .eq. target)
         idxs(1) = floor(grnd()*popsz + 1.d0)
      end do

      do while ((idxs(2) .eq. target) .or. (idxs(2) .eq. idxs(1)))
         idxs(2) = floor(grnd()*popsz + 1.d0)
      end do

      do while ((idxs(3) .eq. target) .or. (idxs(3) .eq. idxs(1)) &
                .or. (idxs(3) .eq. idxs(2)))
         idxs(3) = floor(grnd()*popsz + 1.d0)
      end do

   end function random_sample

   function rand1_mutation(idxs, dims, population, f) result(xmut)

      integer, intent(in) :: idxs(:), dims
      real(8), intent(in) :: population(:, :)
      real(8) :: a(dims), b(dims), c(dims), xmut(dims), f

      a = population(:, idxs(1))
      b = population(:, idxs(2))
      c = population(:, idxs(3))

      xmut = a + f*(b - c)

      ! clip it between [0, 1]
      xmut = min(max(xmut, 0.d0), 1.d0)

   end function rand1_mutation

   function crossover(xmut, crsp, target_vect, dims) result(xtrial)

      real(8), intent(in) :: xmut(:), target_vect(:), crsp
      integer, intent(in) :: dims
      real(8) :: xtrial(dims), cr(dims)
      integer :: i, cp

      ! crossover probability
      do i = 1, dims
         cr(i) = grnd()
      end do

      ! get the trial vector
      cp = 0
      do i = 1, dims
         if (cr(i) .lt. crsp) then
            xtrial(i) = xmut(i)
            cp = cp + 1
         else
            xtrial(i) = target_vect(i)
         end if
      end do

      ! if nothing is less than CR, take a random sample
      if (cp .eq. 0) then
         do i = 1, dims
            xtrial(i) = target_vect(1 + int(popsize*grnd()))
         end do
      end if

   end function crossover

   !
   ! unit testing purposes
   !
   function objfunc(z, ifunc) result(v)

      real(8), intent(in) :: z(:)
      real(8) :: v
      integer :: ifunc

      ! default behaviour
      if (ifunc .eq. 0) call obj_nmr(z, v)

      ! unit testing
      if (ifunc .eq. 1) call ackley(z, v)
      if (ifunc .eq. 2) call beale(z, v)

   end function objfunc

   subroutine ackley(z, v)

      real(8), intent(in) :: z(2)
      real(8), intent(out) :: v

      real(8), parameter :: PI = 4*atan(1.d0)
      real(8) :: a, b, c, x, y

      a = 20.0d0
      b = 0.20d0
      c = 2*PI
      x = z(1)
      y = z(2)

      v = -a*exp(-b*sqrt(0.5d0*(x**2 + y**2))) - &
          exp(0.5d0*(cos(c*x) + cos(c*y))) + a + exp(1.d0)

   end subroutine ackley

   subroutine beale(z, v)

      real(8), intent(in) :: z(2)
      real(8), intent(out) :: v

      real(8), parameter :: PI = 4*atan(1.d0)
      real(8) :: a, b, c, x, y

      a = 1.5d0
      b = 2.25d0
      c = 2.625d0
      x = z(1)
      y = z(2)

      v = (a - x + x*y)**2 + (b - x + x*(y*y))**2 + (c - x + x*(y*y*y))**2

   end subroutine beale

end module de_mod
