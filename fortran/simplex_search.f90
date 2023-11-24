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

module simplex_search

! Module containing simplex search optimization routine

  implicit none

! Options type for direct searches

  type ds_options_type

    double precision :: tol       ! tolerance in simplex radius before
                                  !   triggering a stop condition
    integer :: maxit              ! Max steps allowed before stopping
    logical :: write_designs      ! Whether to write best design each time it
                                  !   changes
    logical :: relative_fmin_report
                                  ! If .true., reports improvement over seed
                                  !   design. Otherwise, reports fmin itself.
  end type ds_options_type

  contains

!=============================================================================80
!
! Nelder-Mead simplex search algorithm
!
!=============================================================================80
subroutine simplexsearch(xopt, fmin, step, fevals, objfunc, x0, given_f0_ref,  &
                         f0_ref, ds_options, restart, restart_write_freq,      &
                         indesigncounter, instep, converterfunc)

  use optimization_util, only : bubble_sort, design_radius_simplex, write_design,      &
                                read_run_control, bubble_sort_plus

  use vardef, only : output_prefix, objfunction_type, write_dvs_file,          &
                     contrain_number

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

  double precision, dimension(:), intent(in) :: x0
  double precision, intent(inout) :: f0_ref
  logical, intent(in) :: given_f0_ref, restart
  type (ds_options_type), intent(in) :: ds_options
  integer, intent(in) :: restart_write_freq
  integer, intent(in), optional :: indesigncounter, instep

  optional :: converterfunc
  interface
    integer function converterfunc(x, designcounter, laststep)
      double precision, dimension(:), intent(in) :: x
      integer, intent(in) :: designcounter
      logical, intent(in) :: laststep
    end function
  end interface

  double precision, dimension(size(x0,1),size(x0,1)+1) :: dv
  double precision, dimension(size(x0,1)+1) :: objvals
  double precision, dimension(size(x0,1)) :: xcen, xr, xe, xc
  
  character(200), dimension(size(x0,1)+1) :: messages
  integer, dimension(size(x0,1)+1) :: message_codes

  double precision :: rho, xi, gam, sigma, fr, fe, fc, f0, mincurr, radius
  integer :: i, j, nvars, stat, designcounter, restartcounter, iunit, ioerr,   &
             prevsteps, k, ncommands
  type(objfunction_type) :: objfunction_return
  logical :: converged, needshrink, signal_progress, new_history_file
  integer :: stepstart, steptime, restarttime
  character(14) :: timechar
  character(3) :: filestat
  character(11) :: stepchar
  character(20) :: fminchar
  character(15) :: radchar
  character(25) :: relfminchar
  character(80), dimension(20) :: commands
  character(100) :: histfile

! Standard Nelder-Mead constants

  rho = 1.d0
  xi = 2.d0
  gam = 0.5d0
  sigma = 0.5d0

! Set up or read initialzation data

  nvars = size(x0,1)

  if (.not. restart) then

!   Get f0 (reference seed design objective function)

    if (given_f0_ref) then
      f0 = f0_ref
    else
      objfunction_return = objfunc(x0,0)
      f0 = objfunction_return%value
      f0_ref = f0
    end if

!   Set up initial simplex

    fevals = 0
    do j = 1, nvars
      do i = 1, nvars
        if (i == j) then
          if (x0(i) == 0.d0) then
            dv(i,j) = 0.00025d0
          else
            dv(i,j) = 1.05d0*x0(i)
          end if
        else
          dv(i,j) = x0(i)
        end if
      end do
      objfunction_return = objfunc(dv(:,j),1)
      objvals(j) = objfunction_return%value
      messages(j) = objfunction_return%message
      message_codes(j) = objfunction_return%message_code
      fevals = fevals + 1
    end do

    dv(:,nvars+1) = x0
    objfunction_return = objfunc(x0,0)
    objvals(nvars+1) = objfunction_return%value
    messages(nvars+1) = objfunction_return%message
    message_codes(nvars+1) = objfunction_return%message_code
    fevals = fevals + 1

!   Counters

    step = 0
    if (.not. present(indesigncounter)) then
      designcounter = 0
    else
      designcounter = indesigncounter
    end if
    if (.not. present(instep)) then
      prevsteps = 0
    else
      prevsteps = instep
    end if

    ! Set restart time to zero
    restarttime =  0.d0
  else

!   Get initial simplex and counters from restart file
    if (.not. present(instep)) then
      prevsteps = 0
    else
      prevsteps = instep
    end if
    
    call simplex_read_restart(step, designcounter, dv, objvals, f0, fevals,    &
                              restarttime, message_codes, messages)

  end if

! Initial minimum value

  fmin = minval(objvals)
  mincurr = fmin

! Open file for writing iteration history
  histfile = trim(output_prefix)//'_optimization_history.dat'
  iunit = 17
  new_history_file = .false.
  if ( (prevsteps == 0) .and. (step == 0) ) then
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
    if (ds_options%relative_fmin_report) then
      write(iunit,'(A)') "Iteration  Objective function  "//&
                         "% Improvement over seed  Design radius"//&
                         "  Time (seconds)"
    else
      write(iunit,'(A)') "Iteration  Objective function  Design radius"//&
                         "  Time (seconds)"
    end if
    flush(iunit)
  end if

  if (step + prevsteps .EQ. 0 ) then
    ! Display initial value
    write(*,'(A12,I5)')   ' Iteration: ', 0
    write(*,'(A27,F12.6)') '   Objective function:    ', f0
    if (ds_options%relative_fmin_report) write(*,'(A27,F12.6,A1)')             &
                        '   Improvement over seed: ', (f0 - f0)/f0*100.d0, '%'
  end if
  
  ! Begin time
  stepstart=time()

  ! Iterative procedure for optimization

  restartcounter = 1
  needshrink = .false.
  converged = .false.
  write(*,*) 'Simplex optimization progress:'

  main_loop: do while (.not. converged)

    step = step + 1
    ! '>=' instead of '==', to make sure it doesn't run after restart when max
    ! iterations has already been reached
    if (step >= ds_options%maxit) converged = .true.

!   Sort according to ascending objective function value

    call bubble_sort_plus(dv, objvals, message_codes, messages)
    mincurr = objvals(1)

!   Update fmin if appropriate

    if (mincurr < fmin) then
      fmin = mincurr
      signal_progress = .true.
    else
      signal_progress = .false.
    end if

!   Check for convergence

    radius = design_radius_simplex(dv)
    if (radius < ds_options%tol) converged = .true.

!   Display progress

    write(*,'(A12,I5)')   ' Iteration: ', step + prevsteps
    write(*,'(A27,F12.6)') '   Objective function:    ', fmin
    if (ds_options%relative_fmin_report) write(*,'(A27,F12.6,A1)')              &
                        '   Improvement over seed: ', (f0 - fmin)/f0*100.d0, '%'
    write(*,'(A27,ES13.3)') '   Design radius:         ', radius

!   Write design to file if requested
!   converterfunc is an optional function supplied to convert design variables
!     into something more useful.  If not supplied, the design variables
!     themselves are written to a file.
    
    if (ds_options%write_designs .and. designcounter == 1) then
      filestat = 'new'
    else
      filestat = 'old'
    end if

    if ( (signal_progress) .and. (ds_options%write_designs) ) then
      designcounter = designcounter + 1
      if (present(converterfunc)) then
        stat = converterfunc(dv(:,1), designcounter, .false.)
      else
        call write_design('simplex_designs.dat', filestat, dv(:,1),            &
                          designcounter)
      end if
    end if
    
    !  Get step time
    steptime=time()
    
!   Write iteration history

    flush(iunit)
    write(stepchar,'(I11)') step + prevsteps
    write(fminchar,'(F14.10)') fmin
    write(radchar,'(ES14.6)') radius
    write(timechar,'(I14)') (steptime-stepstart)+restarttime
    if (ds_options%relative_fmin_report) then
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
      call simplex_write_dvs(step, dv, objvals, message_codes,        &
                             messages, x0, f0, dv(:,1), fmin)
    end if
    
!   Write restart file if appropriate and update restart counter

    if (restartcounter == restart_write_freq) then
      ! 'step' to correct the number on simplex restart
      call simplex_write_restart(step, designcounter, dv, objvals,   &
        f0,fevals, (steptime-stepstart)+restarttime, message_codes,  &
        messages)
      restartcounter = 1
    else
      restartcounter = restartcounter + 1
    end if

!   Check for commands in run_control file

    call read_run_control(commands, ncommands)
    do k = 1, ncommands
      if (trim(commands(k)) == "stop") then
        converged = .true.
        write(*,*) 'Cleaning up: stop command encountered in run_control.'
      end if
    end do

!   Compute the centroid of the best nvals designs

    xcen(:) = 0.d0
    do i = 1, nvars
      xcen = xcen + dv(:,i)
    end do
    xcen = xcen/dble(nvars)

!   Compute the reflection point and evaluate its objective function value

    xr = (1.d0 + rho)*xcen - rho*dv(:,nvars+1)
    objfunction_return = objfunc(xr, step)
    fr = objfunction_return%value
    fevals = fevals + 1

    expand_or_contract: if (objvals(1) <= fr .and. fr < objvals(nvars)) then

!      Accept reflection point

       dv(:,nvars+1) = xr
       objvals(nvars+1) = fr
       messages(nvars+1) = objfunction_return%message
       message_codes(nvars+1) = objfunction_return%message_code
       cycle

    elseif (fr < objvals(1)) then

!     Expand

      xe = (1.d0 + rho*xi)*xcen - rho*xi*dv(:,nvars+1)
      objfunction_return = objfunc(xe, step)
      fe = objfunction_return%value
      fevals = fevals + 1
      if (fe < fr) then
        dv(:,nvars+1) = xe
        objvals(nvars+1) = fe
        messages(nvars+1) = objfunction_return%message
        message_codes(nvars+1) = objfunction_return%message_code
      else
        dv(:,nvars+1) = xr
        objvals(nvars+1) = fr
      end if
      cycle

    elseif (fr >= objvals(nvars)) then

!     Outside contraction

        contraction: if (fr < objvals(nvars+1)) then

          xc = (1.d0 + rho*gam)*xcen - rho*gam*dv(:,nvars+1)
          objfunction_return = objfunc(xc, step)
          fc = objfunction_return%value
          fevals = fevals + 1

          if (fc < fr) then
            dv(:,nvars+1) = xc
            objvals(nvars+1) = fc
            messages(nvars+1) = objfunction_return%message
            message_codes(nvars+1) = objfunction_return%message_code
            needshrink = .false.
          else
            needshrink = .true.
          end if

!       Inside contraction

        else

          xc = (1.d0 - gam)*xcen + gam*dv(:,nvars+1)
          objfunction_return = objfunc(xc, step)
          fc = objfunction_return%value
          fevals = fevals + 1

          if (fc < objvals(nvars+1) ) then
            dv(:,nvars+1) = xc
            objvals(nvars+1) = fc
            messages(nvars+1) = objfunction_return%message
            message_codes(nvars+1) = objfunction_return%message_code
            needshrink = .false.
          else
            needshrink = .true.
          end if

        end if contraction

!       Shrink

        shrink: if (needshrink) then

          do i = 2, nvars + 1
            dv(:,i) = dv(:,1) + sigma*(dv(:,i) - dv(:,1))
            objfunction_return = objfunc(dv(:,i), step)
            objvals(i) = objfunction_return%value
            messages(i) = objfunction_return%message
            message_codes(i) = objfunction_return%message_code
            fevals = fevals + 1
          end do
          cycle

        else

          cycle

        end if shrink

    end if expand_or_contract

  end do main_loop

! Sort one more time according to ascending objective function value

  call bubble_sort(dv, objvals)
  xopt = dv(:,1)
  fmin = objvals(1)
  
!    Write airfoils for each op
    
  if (present(converterfunc)) then
    stat = converterfunc(xopt, designcounter, .true.)
  end if

! Check for convergence one more time

  radius = design_radius_simplex(dv)

! Display warning if max iterations are reached

  if (step == ds_options%maxit .and. (radius >= ds_options%tol)) then
    write(*,*) 'Warning: Simplex optimizer forced to exit due to the max number'
    write(*,*) '         of iterations being reached.'
  end if

! Close iteration history file

  close(iunit)

! Write restart at end of optimization

  if (restartcounter /= 1)                                                     &
    call simplex_write_restart(step, designcounter, dv, objvals, f0,           &
                               fevals, (steptime-stepstart)+restarttime,       &
                               message_codes, messages)

end subroutine simplexsearch

!=============================================================================80
!
! Simplex restart write routine
!
!=============================================================================80
subroutine simplex_write_restart(step, designcounter, dv, objvals, f0, fevals, &
                                 time, message_codes, messages)

  use vardef, only : output_prefix

  integer, intent(in) :: step, designcounter, fevals
  double precision, dimension(:,:), intent(in) :: dv
  double precision, dimension(:), intent(in) :: objvals
  double precision, intent(in) :: f0
  integer, intent(in) :: time
  integer, dimension(:), intent(in) :: message_codes
  character(200), dimension(:), intent(in) :: messages

  character(100) :: restfile
  integer :: iunit

! Status notification

  restfile = 'restart_simplex_'//trim(output_prefix)
  write(*,*) '  Writing simplex restart data to file '//trim(restfile)//' ...'

! Open output file for writing

  iunit = 13
  open(unit=iunit, file=restfile, status='replace', form='unformatted')

! Write restart data

  write(iunit) step
  write(iunit) designcounter
  write(iunit) dv
  write(iunit) objvals
  write(iunit) f0
  write(iunit) fevals
  write(iunit) time
  write(iunit) message_codes
  write(iunit) messages

! Close restart file

  close(iunit)

! Status notification

  write(*,*) '  Successfully wrote simplex restart file.'

end subroutine simplex_write_restart

!=============================================================================80
!
! Particle swarm restart read routine
!
!=============================================================================80
subroutine simplex_read_restart(step, designcounter, dv, objvals, f0, fevals,  &
                                 time, message_codes, messages)
  use vardef, only : output_prefix

  integer, intent(out) :: step, designcounter, fevals
  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objvals
  double precision, intent(out) :: f0
  integer, intent(out) :: time
  integer, dimension(:), intent(inout) :: message_codes
  character(200), dimension(:), intent(inout) :: messages

  character(100) :: restfile
  integer :: iunit, ioerr

! Status notification

  restfile = 'restart_simplex_'//trim(output_prefix)
  write(*,*) 'Reading simplex restart data from file '//trim(restfile)//' ...'

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
  read(iunit) objvals
  read(iunit) f0
  read(iunit) fevals
  read(iunit) time
  read(iunit) message_codes
  read(iunit) messages

! Close restart file

  close(iunit)

! Status notification

  write(*,*) 'Successfully read simplex restart data.'
  write(*,*)

end subroutine simplex_read_restart

!=============================================================================80
!
! Genetic algorithm dvs write routine
!
!=============================================================================80
subroutine simplex_write_dvs(step, dv, objval, message_codes, messages, x0, f0,&
                         xopt, fmin)

  use vardef, only : output_prefix, dvs_for_type

  integer, intent(in) :: step
  double precision, intent(in) :: fmin, f0
  double precision, dimension(:), intent(in) ::  x0, xopt, objval
  double precision, dimension(:,:), intent(in) :: dv
  integer, dimension(:), intent(in) :: message_codes
  character(200), dimension(:), intent(in) :: messages

  integer :: i,j
  
  character(100) :: dvsfile, text, textdv
  integer :: iunit
  
  ! Status notification

  dvsfile = 'dvs_simplex_'//trim(output_prefix)//'.dat'
  write(*,*) '  Writing Simplex dvs data to file '//trim(dvsfile)//' ...'
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
  
  return
  
  ! Warning if there was an error opening dvs file

  900 write(*,*) "Warning: unable to open "//trim(dvsfile)//". Skipping ..."
  return
  
end subroutine simplex_write_dvs  
  
end module simplex_search
