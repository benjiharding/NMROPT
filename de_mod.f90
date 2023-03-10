module de_mod

   use objective_mod
   use network_mod
   use mtmod
   use subs
   use constants

   implicit none

   ! DE parameters
   real(8) :: mut, cplo, cphi, crossp
   integer :: popsize, its
   real(8), allocatable :: best(:)

   ! objective iterations output file
   integer :: lobj = 5

contains

   subroutine optimize()

      real(8) :: start, finish

      write (*, *) " "
      write (*, *) "Starting differential evolution..."
      write (*, *) " "

      call cpu_time(start)

      call de

      call cpu_time(finish)

      write (*, *) " "
      print '("optimization took ", f5.2, " minutes")', (finish - start)/60

   end subroutine optimize

   subroutine de()

      ! differential evolution

      integer :: i, j, k, dims, cp, best_idx, idxs(3)
      real(8), allocatable :: pop(:, :), pop_denorm(:, :)
      real(8), allocatable :: min_b(:, :), max_b(:, :), diff(:, :)
      real(8), allocatable :: a(:), b(:), c(:), mutant(:)
      real(8), allocatable :: trial(:), trial_denorm(:), cr(:)
      real(8), allocatable :: fitness(:), fobj(:)
      real(8) :: fit, ftry, rand(3)

      ! allocate the population and bounds
      dims = size(vect)
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
         call obj_nmr(pop_denorm(:, i), fit)
         fitness(i) = fit
      end do

      ! get the current best
      best_idx = minloc(fitness, dim=1)
      best = pop_denorm(:, best_idx)

      ! main loop
      do i = 1, its

         ! write out the current value?
         if (modulo(i, 100) .eq. 0) then
            write (*, *) " working on DE iteration", i
            if (i .gt. 1) write (*, *) " current objective value", fobj(i - 1)
         end if

         ! crossover prob b/w cpho and cphi
         crossp = grnd()*(cphi - cplo) + cplo

         ! loop over population and evolve
         do j = 1, popsize

            ! sample without replacement
            rand = 0.d0
            do k = 1, 3
               rand(k) = grnd()
            end do
            idxs = floor(rand*popsize + 1.d0)

            do while (idxs(1) .eq. j)
               idxs(1) = floor(grnd()*popsize + 1.d0)
            end do

            do while ((idxs(2) .eq. j) .or. (idxs(2) .eq. idxs(1)))
               idxs(2) = floor(grnd()*popsize + 1.d0)
            end do

            do while ((idxs(3) .eq. j) .or. (idxs(3) .eq. idxs(1)) &
                      .or. (idxs(3) .eq. idxs(2)))
               idxs(3) = floor(grnd()*popsize + 1.d0)
            end do

            ! mutuate the selected vectors
            a = pop(:, idxs(1))
            b = pop(:, idxs(2))
            c = pop(:, idxs(3))
            mutant = a + mut*(b - c)

            ! clip it between [0, 1]
            mutant = min(max(mutant, 0.d0), 1.d0)

            ! crossover probability
            do k = 1, dims
               cr(k) = grnd()
            end do

            ! get the trial vector
            cp = 0
            do k = 1, dims
               if (cr(k) .lt. crossp) then
                  trial(k) = mutant(k)
                  cp = cp + 1
               else
                  trial(k) = pop(k, j)
               end if
            end do

            ! if nothing is less than CR, take a random sample
            if (cp .eq. 0) then
               do k = 1, dims
                  trial(k) = pop(1 + int(popsize*grnd()), j)
               end do
            end if

            ! denormalize the trial and evaluate
            ftry = 0.d0
            trial_denorm = min_b(:, 1) + trial*diff(:, 1)
            call obj_nmr(trial_denorm, ftry)

            ! did we improve?
            if (ftry .lt. fitness(j)) then
               fitness(j) = ftry
               pop(:, j) = trial

               ! do we need to update the global best?
               if (ftry .lt. fitness(best_idx)) then
                  best_idx = j
                  best = trial_denorm
               end if
            end if

         end do

         ! store this iterations best value
         fobj(i) = fitness(best_idx)

         ! write out the values
         write (lobj, "(1(i0,1x),1(g14.8,1x))") i, fobj(i)

         ! are we done early?
         if (fobj(i) .le. 0.0000000000000001D0) exit

      end do

   end subroutine de

end module de_mod
