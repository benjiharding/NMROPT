module de_mod

   use geostat
   use objective_mod, only: obj_nmr, obj_nmr_avg, pobj_nmr, pobj_nmr_avg
   use types_mod
   use mtmod, only: grnd
   use constants
   use omp_lib

   implicit none

contains

   subroutine optimize(ipara)

      integer, intent(in) :: ipara
      real(8) :: start, finish

      write (*, *) " "
      write (*, *) "Starting differential evolution..."
      write (*, *) " "

      if (ipara .eq. 0) then

         call cpu_time(start)

         call de(nnet%dims, popsize, its, mut, cplo, cphi, bmin, bmax, best)

         call cpu_time(finish)

      else if (ipara .eq. 1) then

         call omp_set_num_threads(num_threads)

         start = omp_get_wtime()

         call pde(nnet%dims, popsize, its, mut, cplo, cphi, bmin, bmax, best)

         finish = omp_get_wtime()

      end if

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
      real(8), allocatable :: trial(:), trial_denorm(:)
      real(8), allocatable :: fitness(:), fobj(:)
      real(8) :: fit, ftry, crossp

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
         crossp = cphi + (cplo - cphi)*(1 - dble(i)/dble(its))**4

         ! loop over population and evolve
         do j = 1, popsize

            ! select candidates for mutation
            idxs = selection(popsize, j)

            ! mutuate the selected vector
            mutant = c2b1_mutation(idxs, best_idx, j, dims, pop, mut, mut)

            ! crossover and get the trial vector
            trial = crossover(mutant, crossp, pop(:, j), dims)

            ! denormalize the trial and evaluate
            ftry = 0.d0
            trial_denorm = min_b(:, 1) + trial*diff(:, 1)
            ftry = objfunc(trial_denorm, func)

            ! did we improve?
            if (ftry .lt. fitness(j)) then
               fitness(j) = ftry
               pop(:, j) = trial
               pop_denorm(:, j) = trial_denorm

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
               idxs = selection(popsize, k)
               mutant = c2b1_mutation(idxs, best_idx, k, dims, pop, mut, mut)
               pop(:, k) = crossover(mutant, crossp, pop(:, k), dims)

            end do

            nochange = 0

         end if

         if (.not. present(ifunc)) then
            ! write out the values
            ! write (lobj, "(1(i0,1x),1(g14.8,1x))") i, fobj(i)/fobj(1)
            write (lobj, "(1(i0,1x),*(g14.8,1x))") i, fobj(i)/fobj(1), fobj(i), &
               pop_denorm(1, 1), pop_denorm(2, 1), fitness(1), crossp, mut
         end if

         ! are we done early?
         if (fobj(i) .le. 1e-10) return

      end do

   end subroutine de

   subroutine pde(dims, popsize, its, mut, cplo, cphi, bmin, bmax, best, ifunc)

      ! parallel differential evolution

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
      real(8), allocatable :: trial(:)
      real(8), allocatable :: fitness(:), fobj(:)
      real(8) :: crossp

      ! parallel local variables
      real(8), allocatable :: mutant_loc(:), pfit(:)
      real(8), allocatable :: trial_loc(:), trial_denorm_loc(:)
      real(8), allocatable :: trials(:, :), trials_denorm(:, :)
      integer :: idx_loc, idxs_loc(3)
      integer :: id, first, last
      real(8) :: pf, pfmin
      integer :: pidx

      ! interface for passing arbitrary objective functions
      if (.not. present(ifunc)) then
         func = 0
      else
         func = ifunc
      end if

      ! allocate the population and bounds
      allocate (best(dims))
      allocate (pop(dims, popsize))
      allocate (trials(dims, popsize), trials_denorm(dims, popsize))
      allocate (min_b(dims, 1), max_b(dims, 1))
      min_b = bmin
      max_b = bmax

      ! allocate the trial vectors
      allocate (trial(dims))

      ! allocate objective counter
      allocate (fobj(its))
      fobj = huge(pfit)

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
      allocate (fitness(popsize), pfit(popsize))
      fitness = 0.d0
      do i = 1, popsize
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
         crossp = cphi + (cplo - cphi)*(1 - dble(i)/dble(its))**4

         !
         ! begin parallel region
         !
         !$omp PARALLEL DEFAULT(NONE) &
         !$omp FIRSTPRIVATE(nnet, ysimd) &
         !$omp PRIVATE(mutant_loc, trial_loc, trial_denorm_loc, idxs_loc,  &
         !$omp idx_loc, pf, id, first, last, avg_vario, avg_ivario, avg_runs, &
         !$omp avg_npoint) &
         !$omp SHARED(dims, pop, popsize, mut, crossp, min_b, diff, func,  &
         !$omp num_threads, pfit, best_idx, trials, trials_denorm, &
         !$omp yref, ttable, sill, isills, thresholds)

         id = omp_get_thread_num()
         first = (id*popsize)/num_threads + 1
         last = ((id + 1)*popsize)/num_threads

         ! evaluate the fitness for this thread's indices
         do idx_loc = first, last

            idxs_loc = selection(popsize, idx_loc)
            mutant_loc = c2b1_mutation(idxs_loc, best_idx, idx_loc, dims, pop, mut, mut)
            trial_loc = crossover(mutant_loc, crossp, pop(:, idx_loc), dims)
            trial_denorm_loc = min_b(:, 1) + trial_loc*diff(:, 1)
            ! call pobj_nmr(trial_denorm_loc, nnet, ysimd, pf)
            call pobj_nmr_avg(trial_denorm_loc, nnet, ysimd, pf)
            pfit(idx_loc) = pf
            trials(:, idx_loc) = trial_loc
            trials_denorm(:, idx_loc) = trial_denorm_loc

         end do
         !$omp END PARALLEL

         ! update population for next iteration
         do j = 1, popsize
            if (pfit(j) .lt. fitness(j)) then
               fitness(j) = pfit(j)
               pop(:, j) = trials(:, j)
            end if
         end do

         ! update the best
         pfmin = minval(pfit, dim=1)
         pidx = minloc(pfit, dim=1)
         if (pfmin .lt. fitness(best_idx)) then
            best_idx = pidx
            best = trials_denorm(:, pidx)
            nochange = 0
         end if

         ! store this iterations best value
         fobj(i) = fitness(best_idx)

         ! have we reached the nochange threshold?
         if (nochange .gt. maxnochange) then

            write (*, *) " random restart on iteration", i

            ! randomly mutate the population
            do k = 1, popsize

               ! exclude the current best from mutation
               if (k .eq. best_idx) cycle
               idxs = selection(popsize, k)
               mutant = c2b1_mutation(idxs, best_idx, k, dims, pop, mut, mut)
               pop(:, k) = crossover(mutant, crossp, pop(:, k), dims)

            end do

            nochange = 0

         end if

         if (.not. present(ifunc)) then
            ! write out the values
            ! write (lobj, "(1(i0,1x),1(g14.8,1x))") i, fobj(i)/fobj(1)
            write (lobj, "(1(i0,1x),*(g14.8,1x))") i, fobj(i)/fobj(1), fobj(i), &
               pop_denorm(1, 1), pop_denorm(2, 1), fitness(1), crossp, mut
         end if

         ! are we done early?
         if (fobj(i) .le. 1e-10) return

      end do

   end subroutine pde

   function selection(popsz, target) result(idxs)

      integer, intent(in) :: popsz, target
      integer :: idxs(3)
      real(8) :: rand(3)
      integer :: i

      ! three random idxs that != target
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

   end function selection

   function rand1_mutation(idxs, dims, population, f) result(xmut)

      ! DE/rand/1

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

   function c2b1_mutation(idxs, best_idx, curr_idx, dims, population, f1, f2) &
      result(xmut)

      ! DE/current-to-best/1

      integer, intent(in) :: idxs(:), best_idx, curr_idx, dims
      real(8), intent(in) :: population(:, :)
      real(8) :: a(dims), b(dims)
      real(8) :: xmut(dims), xcurr(dims), xbest(dims), f1, f2

      a = population(:, idxs(1))
      b = population(:, idxs(2))
      xcurr = population(:, curr_idx)
      xbest = population(:, best_idx)

      xmut = xcurr + f1*(xbest - xcurr) + f2*(a - b)

      ! clip it between [0, 1]
      xmut = min(max(xmut, 0.d0), 1.d0)

   end function c2b1_mutation

   function crossover(xmut, cp, xtarget, dims) result(xtrial)

      real(8), intent(in) :: xmut(:), xtarget(:), cp
      integer, intent(in) :: dims
      real(8) :: xtrial(dims), cr(dims)
      integer :: i, irand

      ! crossover probability
      do i = 1, dims
         cr(i) = grnd()
      end do

      ! random index to prevent copying target completely
      irand = floor(grnd()*dims + 1.d0)

      ! get the trial vector
      do i = 1, dims
         if ((cr(i) .lt. cp) .or. (irand .eq. i)) then
            xtrial(i) = xmut(i)
         else
            xtrial(i) = xtarget(i)
         end if
      end do

   end function crossover

   !
   ! unit testing purposes
   !
   function objfunc(z, ifunc) result(v)

      real(8), intent(in) :: z(:)
      real(8) :: v
      integer :: ifunc

      ! default behaviour
      ! if (ifunc .eq. 0) call obj_nmr(z, v)
      if (ifunc .eq. 0) call obj_nmr_avg(z, v)

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
