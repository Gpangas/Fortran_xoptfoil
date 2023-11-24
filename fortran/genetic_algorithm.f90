!  This file is part of XOPTFOIL.

!  XOPTFOIL is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.

!  XOPTFOIL is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.

!  You should have received a copy of the GNU General Public License
!  along with XOPTFOIL.  If not, see <http://www.gnu.org/licenses/>.

!  Copyright (C) 2017-2019 Daniel Prosser, 2020-2021 Ricardo Palmeira

module genetic_algorithm

! Module containing genetic algorithm optimization routine

  implicit none

! Options type definition for genetic algorithm

  type ga_options_type

    integer :: pop                ! genetic algorithm population size
    double precision :: tol       ! tolerance in max radius of designs before
                                  !   triggering a stop condition
    integer :: maxit              ! Max steps allowed before stopping
    character(10) :: parents_selection_method 
                                  ! method for selecting designs to reproduce:
                                  !   roulette, tournament, or random
    double precision :: parent_fraction
                                  ! fraction of population selected to 
                                  !   reproduce
    double precision :: roulette_selection_pressure
                                  ! factor to increase likelihood of best
                                  !   designs being selected by roulette
    double precision :: tournament_fraction
                                  ! fraction of population considered in
                                  !   tournament selection of each parent
    double precision :: crossover_range_factor
                                  ! fraction by which parent characteristics
                                  !   can be extrapolated during crossover
    double precision :: mutant_probability
                                  ! probability of mutation occurring in an
                                  !   offspring
    double precision :: chromosome_mutation_rate
                                  ! probability of mutation in a given
                                  !   chromosome for mutants
    double precision :: mutation_range_factor
                                  ! magnitude of change in a chromosome that
                                  !   can occur during mutation, as fraction of
                                  !   xmax - xmin
    logical :: feasible_init      ! Whether to enforce initially feasible
                                  !   designs
    double precision :: feasible_limit
                                  ! Max objective function value below which
                                  !   initial designs are considered feasible
    integer :: feasible_init_attempts
                                  ! Number of attempts to try to get a feasible
                                  !   initial design
    logical :: write_designs      ! Whether to write best design each time it
                                  !   changes
    logical :: relative_fmin_report
                                  ! If .true., reports improvement over seed
                                  !   design. Otherwise, reports fmin itself.
  end type ga_options_type

  contains

!=============================================================================80
!
! Genetic algorithm optimization routine.  Recommended as a first step to 
! determine the vicinity of the global optimum, followed by a local search to
! hone in.
!
!=============================================================================80
subroutine geneticalgorithm(xopt, fmin, step, fevals, objfunc, x0, xmin, xmax, &
                            given_f0_ref, f0_ref, constrained_dvs, ga_options, &
                            restart, restart_write_freq, designcounter,        &
                            stop_reason, converterfunc)

  use optimization_util, only : init_random_seed, initial_designs_mul_2,       &
                                design_radius, write_design,                   &
                                bubble_sort_plus_3, read_run_control
  use vardef, only : output_prefix, write_dvs_file, objfunction_type,          &
                     contrain_number, noppoint, objfunction_return

  double precision, dimension(:), intent(inout) :: xopt
  double precision, intent(out) :: fmin
  integer, intent(out) :: step, fevals

  interface
    type(objfunction_type) function objfunc(x,step)
      import :: objfunction_type
      double precision, dimension(:), intent(in) :: x
      integer, intent(in) :: step
    end function
  end interface

  double precision, dimension(:), intent(in) :: x0, xmin, xmax
  double precision, intent(inout) :: f0_ref
  integer, dimension(:), intent(in) :: constrained_dvs
  logical, intent(in) :: given_f0_ref, restart
  type (ga_options_type), intent(in) :: ga_options
  integer, intent(in) :: restart_write_freq
  integer, intent(out) :: designcounter
  character(14), intent(out) :: stop_reason

  optional :: converterfunc
  interface
    integer function converterfunc(x, designcounter, laststep)
      double precision, dimension(:), intent(in) :: x
      integer, intent(in) :: designcounter
      logical, intent(in) :: laststep
    end function
  end interface

  integer, dimension(:), allocatable :: idxparents
  integer :: nconstrained, i, j, fminloc, var, stat, restartcounter, nparents
  integer :: idxparent1, idxparent2, idxstack1, idxstack2
  integer :: iunit, ioerr, k, ncommands
  double precision :: mincurr, f0, radius, mutate1, mutate2, objchild1,        &
                      objchild2
  double precision, dimension(ga_options%pop) :: objval
  double precision, dimension(size(xmin,1)) :: xrng, child1, child2
  double precision, dimension(size(xmin,1),ga_options%pop) :: dv
  double precision, dimension(:,:), allocatable :: stackdv
  double precision, dimension(:), allocatable :: stackobjval
  integer, dimension(:), allocatable :: message_codes
  character(200), dimension(:), allocatable :: messages
  double precision, dimension(:,:), allocatable :: constrain_matrix
   double precision, dimension(:,:), allocatable :: aero_matrix
  !type(objfunction_type) :: objfunction_return
  logical :: use_x0, converged, signal_progress, new_history_file
  integer :: stepstart, steptime, restarttime
  character(14) :: timechar
  character(11) :: stepchar
  character(20) :: fminchar
  character(15) :: radchar
  character(25) :: relfminchar
  character(80), dimension(20) :: commands
  character(100) :: histfile
  integer, parameter :: CHUNK = 1

  nconstrained = size(constrained_dvs,1)

! Number of parents in each generation (must be an even number >= 2)

  nparents = nint(ga_options%pop*ga_options%parent_fraction/2.d0)*2
  nparents = max(nparents, 2)
  allocate(idxparents(nparents))
  allocate(stackdv(size(xmin,1),ga_options%pop+nparents))
  allocate(stackobjval(ga_options%pop+nparents))
  allocate(message_codes(ga_options%pop+nparents))
  allocate(messages(ga_options%pop+nparents))
  allocate(constrain_matrix(ga_options%pop+nparents,contrain_number))
  allocate(aero_matrix(ga_options%pop+nparents,3*noppoint+1))

! Difference between max and min

  xrng = xmax - xmin

! Get f0 (reference seed design objective function)

  if (given_f0_ref) then
    f0 = f0_ref
  else
    objfunction_return = objfunc(x0,0)
    f0 = objfunction_return%value
    f0_ref = f0
  end if

!$omp parallel default(shared) private(i, j, var, idxparent1, idxparent2,      &
!$omp  child1, child2, mutate1, mutate2, objchild1, objchild2)

! Initialize a random seed

  call init_random_seed()

! Set up initial designs
  use_x0 = .true.
  if (.not. restart) then
    call initial_designs_mul_2(dv, objval, fevals, objfunc, xmin, xmax, use_x0,&
      x0, ga_options%feasible_init, ga_options%feasible_limit,                 &
      ga_options%feasible_init_attempts,                                       &
      message_codes(1:ga_options%pop),                                         &
      messages(1:ga_options%pop),                                              &
      constrain_matrix(1:ga_options%pop,1:contrain_number),                    &
      aero_matrix(1:ga_options%pop,1:3*noppoint+1))
  end if

!$omp master

! Set up or read other initialization data

  if (.not. restart) then

!   Global best so far

    fmin = minval(objval,1)
    fminloc = minloc(objval,1)
    xopt = dv(:,fminloc)
    mincurr = minval(objval,1)

!   Counters

    step = 1
    designcounter = 0
    
    ! Set restart time to zero
    restarttime =  0.d0
  else

!   Read restart data from file

    call ga_read_restart(step, designcounter, dv, objval, fmin, xopt,          &
      restarttime, message_codes, messages, constrain_matrix, aero_matrix)
    mincurr = minval(objval,1)

  end if

! Open file for writing iteration history
  histfile = trim(output_prefix)//'_optimization_history.dat'
  iunit = 17
  new_history_file = .false.
  if (step == 1) then
    new_history_file = .true.
  else
    open(unit=iunit, file=histfile, status='old',            &
         position='append', iostat=ioerr)
    if (ioerr /= 0) then
      write(*,*) 
      write(*,*) "Warning: did not find existing optimization_history.dat file."
      write(*,*) "A new one will be written, but old data will be lost."
      write(*,*)
      new_history_file = .true.
    end if
  end if
  if (new_history_file) then
    open(unit=iunit, file=histfile, status='replace')
    if (ga_options%relative_fmin_report) then
      write(iunit,'(A)') "Iteration  Objective function  "//&
                         "% Improvement over seed  Design radius"//&
                         "  Time (seconds)"
    else
      write(iunit,'(A)') "Iteration  Objective function  Design radius"//&
                         "  Time (seconds)"
    end if
    flush(iunit)
  end if
  
  ! Begin time
  stepstart=time()

! Begin optimization

  restartcounter = 1
  converged = .false.
  write(*,*) 'Genetic algorithm optimization progress:'

  if (.not. restart) then
  
    ! Display initial value
    write(*,'(A12,I5)')   ' Iteration: ', 0
    write(*,'(A27,F12.6)') '   Objective function:    ', f0
    if (ga_options%relative_fmin_report) write(*,'(A27,F12.6,A1)')             &
                        '   Improvement over seed: ', (f0 - f0)/f0*100.d0, '%'
  
    !   Display progress 

    radius = design_radius(dv,xmax,xmin)
    write(*,'(A12,I5)')   ' Iteration: ', step
    write(*,'(A27,F12.6)') '   Objective function:    ', fmin
    if (ga_options%relative_fmin_report) write(*,'(A27,F12.6,A1)')             &
                        '   Improvement over seed: ', (f0 - fmin)/f0*100.d0, '%'
    write(*,'(A27,ES13.3)') '   Design radius:         ', radius


    if (ga_options%write_designs) then
      designcounter = designcounter + 1
      if (present(converterfunc)) then
        stat = converterfunc(xopt, designcounter, .false.)
      else
        call write_design('particleswarm_designs.dat', 'old', xopt,            &
                          designcounter)
      end if
    end if
    
    !  Get step time
    steptime=time()
  
    !   Write iteration history

    write(stepchar,'(I11)') 0
    write(fminchar,'(F14.10)') f0
    write(radchar,'(ES14.6)') 0.0d0
    write(timechar,'(I14)') 0
    if (ga_options%relative_fmin_report) then
      write(relfminchar,'(F14.10)') (f0 - f0)/f0*100.d0
      write(iunit,'(A11,A20,A25,A15,A14)') adjustl(stepchar), adjustl(fminchar),   &
                                            adjustl(relfminchar), adjustl(radchar), &
                                            adjustl(timechar)
    else
      write(iunit,'(A11,A20,A15,A14)') adjustl(stepchar), adjustl(fminchar),          &
                                adjustl(radchar), adjustl(timechar)
    end if
    flush(iunit)
    
    write(stepchar,'(I11)') step
    write(fminchar,'(F14.10)') fmin
    write(radchar,'(ES14.6)') radius
    write(timechar,'(I14)') (steptime-stepstart)+restarttime
    if (ga_options%relative_fmin_report) then
      write(relfminchar,'(F14.10)') (f0 - fmin)/f0*100.d0
      write(iunit,'(A11,A20,A25,A15,A14)') adjustl(stepchar), adjustl(fminchar),   &
                                            adjustl(relfminchar), adjustl(radchar), &
                                            adjustl(timechar)
    else
      write(iunit,'(A11,A20,A15,A14)') adjustl(stepchar), adjustl(fminchar),          &
                                adjustl(radchar), adjustl(timechar)
    end if
    flush(iunit)
    
    !   Write dvs file if asked
    
    if (write_dvs_file) then
      call ga_write_dvs(step, dv, objval, message_codes(1:ga_options%pop),     &
                        messages(1:ga_options%pop),                            &
                        constrain_matrix(1:ga_options%pop,1:contrain_number),  &
                        aero_matrix(1:ga_options%pop,1:3*noppoint+1),          &
                        x0, f0, xopt, fmin)
    end if
  
    !   Write restart file if appropriate and update restart counter

    call ga_write_restart(step, designcounter, dv, objval, fmin, xopt,         &
      (steptime-stepstart)+restarttime, message_codes, messages,               &
      constrain_matrix, aero_matrix)
  
    restartcounter = restartcounter + 1
  else
    !   Display last step

    radius = design_radius(dv,xmax,xmin)
    write(*,'(A17,I5)')   ' Last Iteration: ', step
    write(*,'(A27,F12.6)') '   Objective function:    ', fmin
    if (ga_options%relative_fmin_report) write(*,'(A27,F12.6,A1)')             &
                        '   Improvement over seed: ', (f0 - fmin)/f0*100.d0, '%'
    write(*,'(A27,ES13.3)') '   Design radius:         ', radius
  end if
  
!$omp end master
!$omp barrier

  optimization_loop: do while (.not. converged)

!$omp master

!   Increase iteration counter

    step = step + 1

!   Select parents

    call parents_selection(objval, ga_options%parents_selection_method,        &
                           ga_options%roulette_selection_pressure,             &
                           ga_options%tournament_fraction, idxparents)

!   Put existing designs at front of stacked arrays

    stackdv(:,1:ga_options%pop) = dv
    stackobjval(1:ga_options%pop) = objval

!$omp end master
!$omp barrier

!$omp do SCHEDULE(DYNAMIC,CHUNK)

!   Procreate to generate offspring pairs

    offspring_creation: do i = 1, nparents/2
      idxparent1 = idxparents(2*(i-1)+1)
      idxparent2 = idxparents(2*i)
      call crossover(dv(:,idxparent1), dv(:,idxparent2),                       &
                     ga_options%crossover_range_factor, child1, child2)

!     Mutate offspring

      call random_number(mutate1)
      if (mutate1 > ga_options%mutant_probability)                             &
        call mutate(child1, ga_options%chromosome_mutation_rate,               &
                    ga_options%mutation_range_factor, xrng, child1) 

      call random_number(mutate2)
      if (mutate2 > ga_options%mutant_probability)                             &
        call mutate(child2, ga_options%chromosome_mutation_rate,               &
                    ga_options%mutation_range_factor, xrng, child2) 

!     Bring back to side constraints if necessary

      do j = 1, nconstrained
        var = constrained_dvs(j)
        if (child1(var) < xmin(var)) then
          child1(var) = xmin(var)
        elseif (child1(var) > xmax(var)) then
          child1(var) = xmax(var)
        end if
        if (child2(var) < xmin(var)) then
          child2(var) = xmin(var)
        elseif (child2(var) > xmax(var)) then
          child2(var) = xmax(var)
        end if
      end do

!     Evaluate objective function for offspring
      objfunction_return = objfunc(child1, step)
      objchild1 = objfunction_return%value
      message_codes(ga_options%pop+2*(i-1)+1) = objfunction_return%message_code
      messages(ga_options%pop+2*(i-1)+1) = objfunction_return%message
      constrain_matrix(ga_options%pop+2*(i-1)+1,:) = objfunction_return%constrains_data
      aero_matrix(ga_options%pop+2*(i-1)+1,:) = objfunction_return%constrains_data
      
      objfunction_return = objfunc(child2, step)
      objchild2 = objfunction_return%value
      message_codes(ga_options%pop+2*i) = objfunction_return%message_code
      messages(ga_options%pop+2*i) = objfunction_return%message
      constrain_matrix(ga_options%pop+2*i,:) = objfunction_return%constrains_data
      aero_matrix(ga_options%pop+2*i,:) = objfunction_return%constrains_data

!     Add children at back of stacked arrays

      idxstack1 = 2*(i-1)+1
      idxstack2 = 2*i
      stackdv(:,ga_options%pop+idxstack1) = child1
      stackdv(:,ga_options%pop+idxstack2) = child2
      stackobjval(ga_options%pop+idxstack1) = objchild1
      stackobjval(ga_options%pop+idxstack2) = objchild2

    end do offspring_creation

!$omp end do

!$omp master

!   Sort stacked arrays to put worst designs at the back

    call bubble_sort_plus_3(stackdv, stackobjval, message_codes, messages,     &
      constrain_matrix, aero_matrix)

!   Replace population with best designs from this generation

    dv = stackdv(:,1:ga_options%pop)
    objval = stackobjval(1:ga_options%pop)

!   Update best overall design, if appropriate

    mincurr = objval(1)
    if (mincurr < fmin) then
      xopt = dv(:,1)
      fmin = mincurr
      signal_progress = .true.
    else
      signal_progress = .false.
    end if

!   Display progress

    radius = design_radius(dv,xmax,xmin)
    write(*,'(A12,I5)')   ' Iteration: ', step
    write(*,'(A27,F12.6)') '   Objective function:    ', fmin
    if (ga_options%relative_fmin_report) write(*,'(A27,F12.6,A1)')              &
                        '   Improvement over seed: ', (f0 - fmin)/f0*100.d0, '%'
    write(*,'(A27,ES13.3)') '   Design radius:         ', radius

!   Write design to file if requested
!   converterfunc is an optional function supplied to convert design variables
!     into something more useful.  If not supplied, the design variables
!     themselves are written to a file.

    if ( (signal_progress) .and. (ga_options%write_designs) ) then
      designcounter = designcounter + 1
      if (present(converterfunc)) then
        stat = converterfunc(xopt, designcounter, .false.)
      else
        call write_design('geneticalgorithm_designs.dat', 'old', xopt,         &
                          designcounter)
      end if
    end if
    
    !  Get step time
    steptime=time()
    
!   Write iteration history

    write(stepchar,'(I11)') step
    write(fminchar,'(F14.10)') fmin
    write(radchar,'(ES14.6)') radius
    write(timechar,'(I14)') (steptime-stepstart)+restarttime
    if (ga_options%relative_fmin_report) then
      write(relfminchar,'(F14.10)') (f0 - fmin)/f0*100.d0
      write(iunit,'(A11,A20,A25,A15,A14)') adjustl(stepchar), adjustl(fminchar),   &
                                           adjustl(relfminchar), adjustl(radchar), &
                                           adjustl(timechar)
    else
      write(iunit,'(A11,A20,A15,A14)') adjustl(stepchar), adjustl(fminchar),          &
                                adjustl(radchar), adjustl(timechar)
    end if
    flush(iunit)
    
!   Evaluate convergence

    if ( (radius > ga_options%tol) .and. (step < ga_options%maxit) ) then
      converged = .false.
    else
      converged = .true.
      stop_reason = "completed"
      if (step == ga_options%maxit) then
        write(*,*) 'Warning: Genetic algorithm forced to exit due to the max'
        write(*,*) '         number of iterations being reached.'
      end if
    end if 

!   Write dvs file if asked
    
    if (write_dvs_file) then
      call ga_write_dvs(step, stackdv, stackobjval, message_codes, messages,   &
        constrain_matrix, aero_matrix, x0, f0, xopt, fmin)
    end if
    
!   Write restart file if appropriate and update restart counter

    if (restartcounter == restart_write_freq) then
      call ga_write_restart(step, designcounter, dv, objval, fmin, xopt,       &
        (steptime-stepstart)+restarttime, message_codes, messages,             &
        constrain_matrix, aero_matrix)
      restartcounter = 1
    else
      restartcounter = restartcounter + 1
    end if
    
!   Check for commands in run_control file

    call read_run_control(commands, ncommands)
    do k = 1, ncommands
      if (trim(commands(k)) == "stop") then
        converged = .true.
        stop_reason = "stop_requested"
        write(*,*) 'Cleaning up: stop command encountered in run_control.'
      end if
    end do
    


!$omp end master
!$omp barrier

  end do optimization_loop

!$omp end parallel
  
!    Write airfoils for each op
    
    if (present(converterfunc)) then
      stat = converterfunc(xopt, designcounter, .true.)
    end if

! Deallocate memory

  deallocate(idxparents)
  deallocate(stackdv)
  deallocate(stackobjval)

! Calculate number of function evaluations

  fevals = fevals + step*nparents

! Close iteration history file

  close(iunit)

! Write restart at end of optimization

  if (restartcounter /= 1)                                                     &
    call ga_write_restart(step, designcounter, dv, objval, fmin, xopt,         &
    (steptime-stepstart)+restarttime, message_codes, messages,                 &
    constrain_matrix, aero_matrix)

end subroutine geneticalgorithm

!=============================================================================80
!
! Selects unique parent designs for reproduction
!
!=============================================================================80
subroutine parents_selection(objval, method, beta, tournament_fraction,        &
                             idxparents)

  use optimization_util, only : pop_double_vector, pop_integer_vector

  double precision, dimension(:), intent(in) :: objval
  character(*), intent(in) :: method
  double precision, intent(in) :: beta, tournament_fraction
  integer, dimension(:), intent(inout) :: idxparents

  integer :: i, ndesigns, nparents, nconsidered, idx
  integer, dimension(size(objval,1)) :: idxconsidered
  double precision, dimension(size(objval,1)) :: objvalconsidered

! Initialize vectors

  ndesigns = size(objval,1)
  nparents = size(idxparents,1)
  objvalconsidered = objval
  do i = 1, ndesigns
    idxconsidered(i) = i
  end do

! Select nparents designs

  nconsidered = ndesigns

  do i = 1, nparents

!   Select a single design as a parent

    call single_parent_selection(objvalconsidered, nconsidered, method, beta, &
                                 tournament_fraction, idx)
    idxparents(i) = idxconsidered(idx)

!   Pop selected parent out of vectors

    call pop_integer_vector(idxconsidered, nconsidered, idx)
    call pop_double_vector(objvalconsidered, nconsidered, idx)
    nconsidered = nconsidered - 1

  end do

end subroutine parents_selection

!=============================================================================80
!
! Selects a single design for reproduction
!
!=============================================================================80
subroutine single_parent_selection(objvalconsidered, nconsidered, method, beta,&
                                   tournament_fraction, idx)

  double precision, dimension(:), intent(in) :: objvalconsidered
  integer, intent(in) :: nconsidered
  character(*), intent(in) :: method
  double precision, intent(in) :: beta, tournament_fraction
  integer, intent(out) :: idx

  if (trim(method) == 'random') then
    call random_selection(nconsidered, idx)
  else if (trim(method) == 'tournament') then
    call tournament_selection(objvalconsidered, nconsidered,                   &
                              tournament_fraction, idx)
  else
    call roulette_selection(objvalconsidered, nconsidered, beta, idx)
  end if

end subroutine single_parent_selection

!=============================================================================80
!
! Roulette wheel selection based on objective function value
!
!=============================================================================80
subroutine roulette_selection(objvalconsidered, nconsidered, beta, idx)

  double precision, dimension(:), intent(in) :: objvalconsidered
  integer, intent(in) :: nconsidered
  double precision, intent(in) :: beta
  integer, intent(out) :: idx

  integer i
  double precision :: total, worst, selection, cumsum, p
  logical :: selected
 
! Sum of all probabilities

  total = 0.d0
  worst = maxval(objvalconsidered(1:nconsidered))
  do i = 1, nconsidered
    total = total + exp(-beta*objvalconsidered(i)/worst)
    if (objvalconsidered(i) < 0.d0) then 
      write(*,*)
      write(*,*) "Error: roulette selection cannot be used with negative"
      write(*,*) "objective function values. Try tournament instead."
      write(*,*)
      stop
    end if
  end do

! Select a random number between 0 and 1

  call random_number(selection)

! Determine which design was selected

  cumsum = 0.d0
  idx = 1
  selected = .false.
  do while (.not. selected)

!   Probability of the current design being selected

    p = exp(-beta*objvalconsidered(idx)/worst) / total

!   Add to cumulative sum and determine whether the random number has been
!   exceeded

    cumsum = cumsum + p
    if (cumsum >= selection) then
      selected = .true.
    else
      idx = idx + 1
    end if

  end do

end subroutine roulette_selection

!=============================================================================80
!
! Tournament selection based on objective function value
!
!=============================================================================80
subroutine tournament_selection(objvalconsidered, nconsidered,                 &
                                tournament_fraction, idx)

  use math_deps,         only : random_integer
  use optimization_util, only : pop_integer_vector

  double precision, dimension(:), intent(in) :: objvalconsidered
  integer, intent(in) :: nconsidered
  double precision, intent(in) :: tournament_fraction
  integer, intent(out) :: idx

  integer :: i, nparticipants, nselected, nremaining, nextidx, nextparticipant
  integer, dimension(nconsidered) :: designstemp
  double precision :: mincurr

! Set number of participants

  nparticipants = nint(dble(nconsidered)*tournament_fraction)
  nparticipants = max(nparticipants,1)

! Temporary list to store designs not yet in the tournament

  do i = 1, nconsidered
    designstemp(i) = i
  end do

! Choose best among nparticipants random designs

  mincurr = 1.D+08
  nselected = 0
  nremaining = nconsidered
  do i = 1, nparticipants

!   Pick a random design from the remaining designs

    nextidx = random_integer(1, nremaining)
    nextparticipant = designstemp(nextidx)

!   Evaluate fitness

    if (objvalconsidered(nextparticipant) < mincurr) then
      mincurr = objvalconsidered(nextparticipant)
      idx = nextparticipant
    end if

!   Pop selected participant out of temp list

    call pop_integer_vector(designstemp, nremaining, nextidx)
    nremaining = nremaining - 1

  end do

end subroutine tournament_selection

!=============================================================================80
!
! Random parent selection
!
!=============================================================================80
subroutine random_selection(nconsidered, idx)

  use math_deps, only : random_integer

  integer, intent(in) :: nconsidered
  integer, intent(out) :: idx

  idx = random_integer(1, nconsidered)

end subroutine random_selection

!=============================================================================80
!
! Crossover of two parents
!
!=============================================================================80
subroutine crossover(parent1, parent2, gam, child1, child2)

  use math_deps, only : random_double

  double precision, dimension(:), intent(in) :: parent1, parent2
  double precision, dimension(:), intent(inout) :: child1, child2
  double precision, intent(in) :: gam

  integer :: nvars, i
  double precision :: alpha

  nvars = size(parent1,1)

! Crossover of each design variable

  do i = 1, nvars
  
!   Crossover parameter

    alpha = random_double(-gam, 1.d0+gam)

!   Child design variables as linear combination of parents

    child1(i) = alpha*parent1(i) + (1.d0-alpha)*parent2(i)
    child2(i) = (1.d0-alpha)*parent1(i) + alpha*parent2(i)

  end do

end subroutine crossover

!=============================================================================80
!
! Mutation of a design
!
!=============================================================================80
subroutine mutate(design, mu, sigma, varrange, newdesign)

  use math_deps, only : random_double

  double precision, dimension(:), intent(in) :: design, varrange
  double precision, intent(in) :: mu, sigma
  double precision, dimension(:), intent(inout) :: newdesign

  integer :: nvars, i
  double precision :: selvar, mutation

  nvars = size(design,1)
  newdesign = design

! Mutate each design variable with probability mu

  do i = 1, nvars

!   Crossover parameter

    call random_number(selvar)
    if (selvar > mu) cycle

!   Mutate if dictated by selvar

    mutation = sigma*varrange(i)*random_double(-0.5d0, 0.5d0)
    newdesign(i) = design(i) + mutation

  end do

end subroutine mutate
  
!=============================================================================80
!
! Genetic algorithm restart write routine
!
!=============================================================================80
subroutine ga_write_restart(step, designcounter, dv, objval, fmin, xopt, time, &
  message_codes, messages, constrain_matrix, aero_matrix)

  use vardef, only : output_prefix

  integer, intent(in) :: step, designcounter
  double precision, dimension(:,:), intent(in) :: dv
  double precision, dimension(:), intent(in) :: objval, xopt
  double precision, intent(in) :: fmin
  integer, intent(in) :: time
  integer, dimension(:), intent(in) :: message_codes
  character(200), dimension(:), intent(in) :: messages
  double precision, dimension(:,:), intent(inout) :: constrain_matrix
  double precision, dimension(:,:), intent(inout) :: aero_matrix

  character(100) :: restfile
  integer :: iunit

!  Status notification

  restfile = 'restart_ga_'//trim(output_prefix)
  write(*,*) '  Writing genetic algorithm restart data to file '//&
             trim(restfile)//' ...'
 
! Open output file for writing

  iunit = 13
  open(unit=iunit, file=restfile, status='replace', form='unformatted')

! Write restart data

  write(iunit) step
  write(iunit) designcounter
  write(iunit) dv
  write(iunit) objval
  write(iunit) fmin
  write(iunit) xopt
  write(iunit) time
  write(iunit) message_codes
  write(iunit) messages
  write(iunit) constrain_matrix
  write(iunit) aero_matrix
  

! Close restart file

  close(iunit)

! Status notification

  write(*,*) '  Successfully wrote genetic algorithm restart file.'

end subroutine ga_write_restart

!=============================================================================80
!
! Genetic algorithm restart read routine
!
!=============================================================================80
subroutine ga_read_restart(step, designcounter, dv, objval, fmin, xopt, time,  &
  message_codes, messages, constrain_matrix, aero_matrix)

  use vardef, only : output_prefix

  integer, intent(inout) :: step, designcounter
  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objval, xopt
  double precision, intent(out) :: fmin
  integer, intent(out) :: time
  integer, dimension(:), intent(inout) :: message_codes
  character(200), dimension(:), intent(inout) :: messages
  double precision, dimension(:,:), intent(inout) :: constrain_matrix
  double precision, dimension(:,:), intent(inout) :: aero_matrix

  character(100) :: restfile
  integer :: iunit, ioerr

! Status notification

  restfile = 'restart_ga_'//trim(output_prefix)
  write(*,*) 'Reading genetic algorithm restart data from file '//&
             trim(restfile)//' ...'
 
! Open output file for reading

  iunit = 13
  open(unit=iunit, file=restfile, status='old', form='unformatted',            &
       iostat=ioerr)
  if (ioerr /= 0) then
    write(*,*) 'Error: could not find input file '//trim(restfile)//'.'
    write(*,*)
    stop
  end if

! Read restart data

  read(iunit) step
  read(iunit) designcounter
  read(iunit) dv
  read(iunit) objval
  read(iunit) fmin
  read(iunit) xopt
  read(iunit) time
  read(iunit) message_codes
  read(iunit) messages
  read(iunit) constrain_matrix
  read(iunit) aero_matrix

! Close restart file

  close(iunit)

! Status notification

  write(*,*) 'Successfully read genetic algorithm restart file.'
  write(*,*)

end subroutine ga_read_restart

!=============================================================================80
!
! Genetic algorithm restart read routine
!
!=============================================================================80
subroutine ga_read_step(step)

  use vardef, only : output_prefix

  integer, intent(inout) :: step

  character(100) :: restfile
  integer :: iunit, ioerr

! Status notification

  restfile = 'restart_ga_'//trim(output_prefix)
  write(*,*) 'Reading genetic algorithm step data from file '//&
             trim(restfile)//' ...'
 
! Open output file for reading

  iunit = 13
  open(unit=iunit, file=restfile, status='old', form='unformatted',            &
       iostat=ioerr)
  if (ioerr /= 0) then
    write(*,*) 'Error: could not find input file '//trim(restfile)//'.'
    write(*,*)
    stop
  end if

! Read restart data

  read(iunit) step

! Close restart file

  close(iunit)

! Status notification

  write(*,*) 'Successfully read step algorithm restart file.'
  write(*,*)

end subroutine ga_read_step
!=============================================================================80
!
! Genetic algorithm dvs write routine
!
!=============================================================================80
subroutine ga_write_dvs(step, dv, objval, message_codes, messages,             &
  constrain_matrix, aero_matrix, x0, f0, xopt, fmin)

  use vardef, only : output_prefix, dvs_for_type, naddthickconst, noppoint,    &
                     ndrag_constrain, nmoment_constrain, nlift_constrain

  integer, intent(in) :: step
  double precision, intent(in) :: fmin, f0
  double precision, dimension(:), intent(in) ::  x0, xopt, objval
  double precision, dimension(:,:), intent(in) :: dv
  integer, dimension(:), intent(in) :: message_codes
  character(200), dimension(:), intent(in) :: messages
  double precision, dimension(:,:), intent(in) :: constrain_matrix
  double precision, dimension(:,:), intent(in) :: aero_matrix

  integer :: i,j
  
  character(100) :: dvsfile, constfile, aerofile, text, textdv
  integer :: iunit
  
  ! Status notification
  
  dvsfile = 'dvs_ga_'//trim(output_prefix)//'.dat'
  constfile = 'constrain_ga_'//trim(output_prefix)//'.dat'
  aerofile = 'aero_ga_'//trim(output_prefix)//'.dat'
  
  write(*,'(A)',advance='no') '  Writing GA log data to files '//trim(dvsfile)
  if (size(constrain_matrix,2) .NE. 1) then
    write(*,'(A)', advance='no')' and '//trim(constfile)
    write(*,'(A)', advance='no')' and '//trim(aerofile)//' ...'
  end if
  write(*,*)
  
  ! Log file
  ! Open files and write headers, if necessary
  iunit = 13

  ! Open files and write headers, if necessary

  write(text,*) step
  text = adjustl(text)
  
  if (step == 1) then

    !   Header for dvs file

    open(unit=iunit, file=dvsfile, status='replace')
    write(iunit,'(A)') 'title="Log file"'
    write(iunit,'(A)') 'variables="dvs vector", "objval", "message_code", "message"'
    write(iunit,'(A)') 'step = 0: x0'
    if (dvs_for_type(2) .NE. 0) then
      do i=1, dvs_for_type(2)
        write(textdv,*) i 
        textdv=adjustl(textdv)
        write(iunit,'(A12)', advance='no') 'top - '//trim(textdv)
      end do
    end if
    if (dvs_for_type(3) .NE. 0) then
      do i=1, dvs_for_type(3)
        write(textdv,*) i 
        textdv=adjustl(textdv)
        write(iunit,'(A12)', advance='no') 'bot - '//trim(textdv)
      end do
    end if
    if (dvs_for_type(4) .NE. 0) then
      do i=1, dvs_for_type(4)
        write(textdv,*) i 
        textdv=adjustl(textdv)
        write(iunit,'(A12)', advance='no') 'fdef - '//trim(textdv)
      end do
    end if
    if (dvs_for_type(5) .NE. 0) write(iunit,'(A12)', advance='no') 'fhinge'
    if (dvs_for_type(6) .NE. 0) write(iunit,'(A12)', advance='no') 'te_thick'
    write(iunit,'(A30)', advance='no') 'f0'
    write(iunit,'(A)') ' '
  
    do i = 1, size(dv,1)
      write(iunit,'(F12.8)', advance='no') x0(i)
    end do
    write(iunit,'(F30.8)', advance='no') f0
    write(iunit,'(A)') ' '
    
  else

    !   Open dvs file and write zone header

    open(unit=iunit, file=dvsfile, status='old', position='append', err=900)

  end if

  ! Write coordinates to file
  write(iunit,'(A)') 'step = '//trim(text)//': xopt '
  if (dvs_for_type(2) .NE. 0) then
    do i=1, dvs_for_type(2)
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A12)', advance='no') 'top - '//trim(textdv)
    end do
  end if
  if (dvs_for_type(3) .NE. 0) then
    do i=1, dvs_for_type(3)
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A12)', advance='no') 'bot - '//trim(textdv)
    end do
  end if
  if (dvs_for_type(4) .NE. 0) then
    do i=1, dvs_for_type(4)
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A12)', advance='no') 'fdef - '//trim(textdv)
    end do
  end if
  if (dvs_for_type(5) .NE. 0) write(iunit,'(A12)', advance='no') 'fhinge'
  if (dvs_for_type(6) .NE. 0) write(iunit,'(A12)', advance='no') 'te_thick'
  write(iunit,'(A30)', advance='no') 'fmin'
  write(iunit,'(A)') ' '

  do i =1, size(dv,1)
    write(iunit,'(F12.8)', advance='no') xopt(i)
  end do
  write(iunit,'(F30.8)', advance='no') fmin
  write(iunit,'(A)') ' '
  
  write(iunit,'(A)') 'step = '//trim(text)//': dv '
  if (dvs_for_type(2) .NE. 0) then
    do i=1, dvs_for_type(2)
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A12)', advance='no') 'top - '//trim(textdv)
    end do
  end if
  if (dvs_for_type(3) .NE. 0) then
    do i=1, dvs_for_type(3)
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A12)', advance='no') 'bot - '//trim(textdv)
    end do
  end if
  if (dvs_for_type(4) .NE. 0) then
    do i=1, dvs_for_type(4)
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A12)', advance='no') 'fdef - '//trim(textdv)
    end do
  end if
  if (dvs_for_type(5) .NE. 0) write(iunit,'(A12)', advance='no') 'fhinge'
  if (dvs_for_type(6) .NE. 0) write(iunit,'(A12)', advance='no') 'te_thick'
  write(iunit,'(A30)', advance='no') 'objval'
  write(iunit,'(A14)', advance='no') 'message_code'
  write(iunit,'(A)', advance='no') ' message'
  write(iunit,'(A)') ' '
  
  do i = 1, size(dv,2)
    do j =1, size(dv,1)
      write(iunit,'(F12.8)', advance='no') dv(j,i)
    end do
    write(iunit,'(F30.8)', advance='no') objval(i)
    write(iunit,'(I14)', advance='no') message_codes(i)
    write(iunit,'(A)', advance='no') messages(i)
    write(iunit,'(A)') ' '
  end do

  ! Close output files

  close(iunit)

  ! Status notification
  
  if (size(constrain_matrix,2) .NE. 1) then
    
    ! Constrains file
    ! Open files and write headers, if necessary
    iunit = 14
    write(text,*) step
    text = adjustl(text)
  
    if (step == 1) then

      !   Header for Constrains file

      open(unit=iunit, file=constfile, status='replace')
      write(iunit,'(A)') 'title="Constrains file"'
      write(iunit,'(A)') 'variables="constrains", "objval", "message_code", "message"'
    
    else

      !   Open dvs file and write zone header

      open(unit=iunit, file=constfile, status='old', position='append', err=901)

    end if

    ! Write coordinates to file
    write(iunit,'(A)') 'step = '//trim(text)//': dvs '
    write(iunit,'(A14)', advance='no') 'tcTE'
    write(iunit,'(A14)', advance='no') 'flap_deg'
    write(iunit,'(A14)', advance='no') 'flap_hinge'
    write(iunit,'(A14)', advance='no') 'minthick'
    write(iunit,'(A14)', advance='no') 'maxthick'
    
    if (naddthickconst .NE. 0) then
      do i=1, naddthickconst
        write(textdv,*) i 
        textdv=adjustl(textdv)
        write(iunit,'(A14)', advance='no') 'addthickconst - '//trim(textdv)
      end do
    end if
    
    write(iunit,'(A14)', advance='no') 'TE_angle'
    write(iunit,'(A14)', advance='no') 'maxcamb'
    write(iunit,'(A14)', advance='no') 'maxpang'
    write(iunit,'(A14)', advance='no') 'minang'
    write(iunit,'(A14)', advance='no') 'difang'
    write(iunit,'(A14)', advance='no') 'rev_t'
    write(iunit,'(A14)', advance='no') 'rev_b'
    write(iunit,'(A14)', advance='no') 'pan_ang'
    write(iunit,'(A14)', advance='no') 'maxgrowth'
    write(iunit,'(A14)', advance='no') 'n_unconv'
    
    if (nmoment_constrain .NE. 0) then
      do i=1, nmoment_constrain
        write(textdv,*) i 
        textdv=adjustl(textdv)
        write(iunit,'(A14)', advance='no') 'moment - '//trim(textdv)
      end do
    end if
    if (nlift_constrain .NE. 0) then
      do i=1, nlift_constrain
        write(textdv,*) i 
        textdv=adjustl(textdv)
        write(iunit,'(A14)', advance='no') 'lift - '//trim(textdv)
      end do
    end if
        if (ndrag_constrain .NE. 0) then
      do i=1, ndrag_constrain
        write(textdv,*) i 
        textdv=adjustl(textdv)
        write(iunit,'(A14)', advance='no') 'drag - '//trim(textdv)
      end do
    end if
    write(iunit,'(A30)', advance='no') 'objval'
    write(iunit,'(A14)', advance='no') 'message_code'
    write(iunit,'(A)', advance='no') ' message'
    write(iunit,'(A)') ' '
  
    do i = 1, size(constrain_matrix,1)
      do j = 1, size(constrain_matrix,2)
        write(iunit,'(F14.8)', advance='no') constrain_matrix(i,j)
      end do
      write(iunit,'(F30.8)', advance='no') objval(i)
      write(iunit,'(I14)', advance='no') message_codes(i)
      write(iunit,'(A)', advance='no') messages(i)
      write(iunit,'(A)') ' '
    end do

    ! Close output files

    close(iunit)
  end if
  
if (size(constrain_matrix,2) .NE. 1) then
    ! Aero file
    ! Open files and write headers, if necessary
    iunit = 15
    write(text,*) step
    text = adjustl(text)
  
    if (step == 1) then

      !   Header for aero file

      open(unit=iunit, file=aerofile, status='replace')
      write(iunit,'(A)') 'title="Aero file"'
      write(iunit,'(A)') 'variables="Aerodynamic proprerties", '//             &
        &'"partial objective function", "partial improvement", '//             &
        &'"objval", "message_code", "message"'
    
    else

      !   Open aero file and write zone header

      open(unit=iunit, file=aerofile, status='old', position='append', err=902)

    end if

    ! Write coordinates to file
    write(iunit,'(A)') 'step = '//trim(text)//': aero vector '
    do i=1, noppoint
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A14)', advance='no') 'aero prop - '//trim(textdv)
    end do
    do i=1, noppoint
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A14)', advance='no') 'part obj - '//trim(textdv)
    end do
    do i=1, noppoint
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A14)', advance='no') 'part imp - '//trim(textdv)
    end do
    write(iunit,'(A14)', advance='no') 'obj no penal'

    write(iunit,'(A30)', advance='no') 'objval'
    write(iunit,'(A14)', advance='no') 'message_code'
    write(iunit,'(A)', advance='no') ' message'
    write(iunit,'(A)') ' '
  
    do i = 1, size(aero_matrix,1)
      do j = 1, size(aero_matrix,2)
        write(iunit,'(F14.8)', advance='no') aero_matrix(i,j)
      end do
      write(iunit,'(F30.8)', advance='no') objval(i)
      write(iunit,'(I14)', advance='no') message_codes(i)
      write(iunit,'(A)', advance='no') messages(i)
      write(iunit,'(A)') ' '
    end do

    ! Close output files

    close(iunit)
  end if
  
    ! Status notification

    write(*,*) '  Successfully wrote GA log data.'
    write(*,*)
  
  return
  
  ! Warning if there was an error opening dvs file

  900 write(*,*) "Warning: unable to open "//trim(dvsfile)//". Skipping ..."
  return
  901 write(*,*) "Warning: unable to open "//trim(constfile)//". Skipping ..."
  return
  902 write(*,*) "Warning: unable to open "//trim(aerofile)//". Skipping ..."
  return

  
end subroutine ga_write_dvs

end module genetic_algorithm
