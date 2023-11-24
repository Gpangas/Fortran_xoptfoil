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

module optimization_driver

! Contains subroutines to set options and conditions for optimization and to
! issue optimizer calls

  implicit none

  contains

!=============================================================================80
!
! Preprocessing for non-aerodynamic optimization
!
!=============================================================================80
subroutine matchfoils_preprocessing(matchfoil_file)

  use vardef,             only : airfoil_type, xmatcht, xmatchb, zmatcht,      &
                                 zmatchb, xseedt, xseedb, zseedt, zseedb,      &
                                 symmetrical, tcTE, TE_spec, xltTE, match_foil
  use memory_util,        only : deallocate_airfoil
  use airfoil_operations, only : get_seed_airfoil, get_split_points,           &
                                 split_airfoil, my_stop
  use math_deps,          only : interp_vector, spline_interp
  use naca,               only : naca_options_type

  character(*), intent(in) :: matchfoil_file

  type(naca_options_type) :: dummy_naca_options
  integer :: pointst, pointsb, i
  double precision, dimension(:), allocatable :: zttmp, zbtmp
  double precision :: xoffmatch, zoffmatch, scale_match, angle_match, tcTE_match
  double precision ::TE_seed, TE_match, TE_ratio, TE_dif
  
  write(*,*) 'Note: using the optimizer to match the seed airfoil to the '
  write(*,*) 'airfoil about to be loaded.'
  write(*,*)

! Check if symmetrical airfoil was requested (not allowed for this type)

  if (symmetrical)                                                             &
    call my_stop("Symmetrical airfoil constraint not permitted for non-"//&
                 "aerodynamic optimizations.")

! Load airfoil to match

  call get_seed_airfoil('from_file', matchfoil_file, dummy_naca_options,       &
                        match_foil, xoffmatch, zoffmatch, scale_match, angle_match, tcTE_match)

! Split match_foil into upper and lower halves

  call get_split_points(match_foil, pointst, pointsb, .false.)
  allocate(xmatcht(pointst))
  allocate(zmatcht(pointst))
  allocate(xmatchb(pointsb))
  allocate(zmatchb(pointsb))
  call split_airfoil(match_foil, xmatcht, xmatchb, zmatcht, zmatchb, .false.)

! Interpolate x-vals of foil to match to seed airfoil points to x-vals

  pointst = size(xseedt,1)
  pointsb = size(xseedb,1)
  allocate(zttmp(pointst))
  allocate(zbtmp(pointsb))
  zttmp(pointst) = zmatcht(size(zmatcht,1))
  zbtmp(pointsb) = zmatchb(size(zmatchb,1))
    
  ! Set Trailing edge
  TE_match=zmatcht(size(zmatcht,1))-zmatchb(size(zmatchb,1))
  write(*,*) '    Trailing Edge of match foil = ', TE_match
  
  if (trim(TE_spec) == 'specify') then
    do i = 1, size(xmatcht,1)
      if (xmatcht(i) .GT. xltTE) then
        zmatcht(i) = zmatcht(i) + (xmatcht(i)-xltTE) * (tcTE-TE_match) / 2.0d0
      end if
    end do
    do i = 1, size(xmatchb,1)
      if (xmatchb(i) .GT. xltTE) then
        zmatchb(i) = zmatchb(i) - (xmatchb(i)-xltTE) * (tcTE-TE_match) / 2.0d0
      end if
    end do
  
  end if
  
  ! compare TE of match with seed
  TE_seed=zseedt(pointst)-zseedb(pointsb)
  TE_match=zmatcht(size(zmatcht,1))-zmatchb(size(zmatchb,1))
  TE_dif=TE_seed-TE_match
  TE_ratio=TE_seed/TE_match
  
  write(*,*) 'New Trailing Edge of match foil = ', TE_match
  
  ! interpolate points
  
  call spline_interp(size(xmatcht,1), xmatcht, zmatcht,                        &
                     size(xmatchb,1), xmatchb, zmatchb,                        &
                     pointst-1, xseedt(1:pointst-1), zttmp(1:pointst-1),       &
                     pointsb-1, xseedb(1:pointsb-1), zbtmp(1:pointsb-1))

! Re-set coordinates of foil to match from interpolated points
    
  deallocate(xmatcht)
  deallocate(zmatcht)
  deallocate(xmatchb)
  deallocate(zmatchb)
  allocate(xmatcht(pointst))
  allocate(zmatcht(pointst))
  allocate(xmatchb(pointsb))
  allocate(zmatchb(pointsb))
  xmatcht = xseedt
  xmatchb = xseedb
  zmatcht = zttmp
  zmatchb = zbtmp

  deallocate(match_foil%x,match_foil%z)
  allocate(match_foil%x(pointst+pointsb+1),match_foil%z(pointst+pointsb+1))
  
  match_foil%npoint=pointst+pointsb+1
  
  do i = 1, pointst
    match_foil%x(i) = xmatcht(pointst-i+1)
    match_foil%z(i) = zmatcht(pointst-i+1)
  end do
  do i = 1, pointsb-1
    match_foil%x(i+pointst) = xmatchb(i+1)
    match_foil%z(i+pointst) = zmatchb(i+1)
  end do
  
  
! Deallocate temporary arrays

  deallocate(zttmp)
  deallocate(zbtmp)

end subroutine matchfoils_preprocessing

!=============================================================================80
!
! Subroutine to drive the optimization
!
!=============================================================================80
subroutine optimize(search_type, global_search, local_search, constrained_dvs, &
                    pso_options, ga_options, ds_options, restart,              &
                    restart_write_freq, optdesign, f0_ref, fmin, steps, fevals)

  use vardef,             only : output_prefix, global_search_stat,            &
                                 local_search_stat, restart_stat, objfunction_type,&
                                 contrain_number, objfunction_return, noppoint
  use particle_swarm,     only : pso_options_type, particleswarm, pso_read_step
  use genetic_algorithm,  only : ga_options_type, geneticalgorithm, ga_read_step
  use simplex_search,     only : ds_options_type, simplexsearch
  use airfoil_evaluation, only : objective_function,                           &
                                 objective_function_nopenalty, write_function, &
                                 write_function_restart_cleanup
  use parametrization,    only : parametrization_init, parametrization_maxmin
  
  character(*), intent(in) :: search_type, global_search, local_search
  type(pso_options_type), intent(in) :: pso_options
  type(ga_options_type), intent(in) :: ga_options
  type(ds_options_type), intent(in) :: ds_options
  double precision, dimension(:), intent(inout) :: optdesign
  double precision, intent(out) :: f0_ref, fmin
  integer, intent(in) :: restart_write_freq
  integer, dimension(:), intent(in) :: constrained_dvs
  integer, intent(out) :: steps, fevals
  logical, intent(in) :: restart

  double precision, dimension(size(optdesign,1)) :: xmin, xmax, x0
  logical :: restart_temp, write_designs
  integer :: stepsg, fevalsg, stepsl, fevalsl, stat, iunit, ioerr, designcounter
  character(100) :: restart_status_file
  character(19) :: restart_status
  character(14) :: stop_reason
  double precision :: f0_penalty
  character(200) :: message
  double precision, dimension(contrain_number) :: constrain_vector
  double precision, dimension(3*noppoint+1) :: aero_vector

  ! Set search types
  
  global_search_stat=global_search
  local_search_stat=local_search
  
! Delete existing run_control file and rewrite it

  iunit = 23
  open(unit=iunit, file='run_control_'//trim(output_prefix), status='replace')
  close(iunit)

! Restart status file setup

  iunit = 15
  restart_status_file = 'restart_status_'//trim(output_prefix)

! Perform optimization: global, local, or global + local

  stepsg = 0
  fevalsg = 0
  stepsl = 0
  fevalsl = 0
  designcounter = 0

! Set initial design

  call parametrization_init(optdesign, x0)

  write(*,*)
  write(*,*) 'x0: '
  write(*,*) x0
  
! Set up mins and maxes
  
  call parametrization_maxmin(optdesign, xmin, xmax)
  
  write(*,*)
  write(*,*) 'xmin: '
  write(*,*) xmin
  write(*,*) 'xmax: '
  write(*,*) xmax
  
! Compute f0_ref, ignoring penalties for violated constraints
  objfunction_return = objective_function_nopenalty(x0) 
  f0_ref = objfunction_return%value
  message = objfunction_return%message
  constrain_vector = objfunction_return%constrains_data
  aero_vector = objfunction_return%aero_data
 
  write(*,*)
  write(*,*) 'f0_ref:     ', f0_ref
  write(*,*) 'message: ', message
  
  !write(*,*) 'constrain_vector', constrain_vector
  !write(*,*) 'aero_vector', aero_vector
  
! Compute with penalty
  objfunction_return = objective_function(x0,0) 
  f0_penalty = objfunction_return%value
  message = objfunction_return%message
  constrain_vector = objfunction_return%constrains_data
  aero_vector = objfunction_return%aero_data
  
  write(*,*)
  write(*,*) 'f0_penalty: ', f0_penalty
  write(*,*) 'message: ', message
  write(*,*) 'constrain_vector', constrain_vector
  write(*,*) 'aero_vector', aero_vector
  
! Wait for user to check information
  write(*,*)
  write(*,*) 'Press Enter to continue'
  read(*,*) 

! Set default restart status (global or local optimization) from user input

  if (trim(search_type) == 'global_and_local' .or. trim(search_type) ==        &
      'global') then
    restart_status = 'global_optimization'
  else
    restart_status = 'local_optimization'
  end if

  restart_stat=restart_status
  
! Read restart status from file for restart case

  if (restart) then

!   Open status file

    open(unit=iunit, file=restart_status_file, status='old', iostat=ioerr)

!   Read or issue warning if file is not found

    if (ioerr /= 0) then
      write(*,*) 'Warning: could not find restart status file '//              &
                 trim(restart_status_file)//'.'
      write(*,*) 'Restarting with '//trim(restart_status)//'.'
    else
      read(iunit,*) restart_status
    end if

!   Close status file

    close(iunit)

  end if

! Design coordinates/polars output handling

  write_designs = .false.
  if ( (trim(search_type) == 'global_and_local') .or.                          &
       (trim(search_type) == 'global') ) then
    if ( (pso_options%write_designs) .or. (ga_options%write_designs) )         &
      write_designs = .true.
  else
    if (ds_options%write_designs) write_designs = .true.
  end if
  
  ! Set up stepsg
  
  if ((restart_status == 'local_optimization') .and.                           &
    (trim(search_type) == 'global_and_local')) then 
    if (trim(global_search) == 'genetic_algorithm') then
      call ga_read_step(stepsg)
    elseif (trim(global_search) == 'particle_swarm') then
      call pso_read_step(stepsg)
    end if
  elseif(restart_status == 'local_optimization') then
    stepsg=0
  end if
  
! Write seed airfoil coordinates and polars to file

  if (write_designs) then
    if (.not. restart) then

!     Analyze and write seed airfoil
  
      stat = write_function(x0, 0, .false.) 

    else

!     Remove unused entries in design polars and coordinates from previous run
 
      stat = write_function_restart_cleanup(restart_status, global_search,     &
                                            local_search, stepsg)

    end if
  end if
  
  !stop
  
! Set temporary restart variable

  restart_temp = restart

! Global optimization

  if (trim(restart_status) == 'global_optimization') then

!   Write restart status to file

    open(unit=iunit, file=restart_status_file, status='replace')
    write(iunit,'(A)') trim(restart_status)
    close(iunit)

    if (trim(global_search) == 'particle_swarm') then

!     Particle swarm optimization

      call particleswarm(optdesign, fmin, stepsg, fevalsg, objective_function, &
                         x0, xmin, xmax, .true., f0_ref, constrained_dvs,      &
                         pso_options, restart_temp, restart_write_freq,        &
                         designcounter, stop_reason, write_function)

    else if (trim(global_search) == 'genetic_algorithm') then

!     Genetic algorithm optimization

      call geneticalgorithm(optdesign, fmin, stepsg, fevalsg,                  &
                            objective_function, x0, xmin, xmax, .true.,        &
                            f0_ref, constrained_dvs, ga_options, restart_temp, &
                            restart_write_freq, designcounter, stop_reason,    &
                            write_function)

    end if

!   Update restart status and turn off restarting for local search

    if ( (stop_reason == "completed") .and.                                    &
         (trim(search_type) == 'global_and_local') )                           &
        restart_status = 'local_optimization'
    restart_temp = .false.

  end if

! Local optimization
  
  if (restart_status == 'local_optimization') then

    if (trim(local_search) == 'simplex') then

!     Write optimization status to file

      open(unit=iunit, file=restart_status_file, status='replace')
      write(iunit,'(A)') trim(restart_status)
      close(iunit)

!     Simplex optimization

      if (trim(search_type) == 'global_and_local') then
        x0 = optdesign  ! Copy x0 from global search result
      end if

      call simplexsearch(optdesign, fmin, stepsl, fevalsl, objective_function, &
                         x0, .true., f0_ref, ds_options, restart_temp,         &
                         restart_write_freq, designcounter, stepsg,            &
                         write_function)

    end if

  end if

! Total number of steps and function evaluations

  steps = stepsg + stepsl
  fevals = fevalsg + fevalsl

! Write stop_monitoring command to run_control file

  iunit = 23
  open(unit=iunit, file='run_control_'//trim(output_prefix), status='old', position='append',        &
       iostat=ioerr)
  if (ioerr /= 0) open(unit=iunit, file='run_control', status='new')
  write(iunit,'(A)') "stop_monitoring"
  close(iunit)

end subroutine optimize

!=============================================================================80
!
! Writes final airfoil design to a file 
!
!=============================================================================80
subroutine write_final_design(optdesign, f0, fmin)

  use vardef
  use memory_util,        only : allocate_airfoil, deallocate_airfoil
  use airfoil_operations, only : airfoil_write
  use parametrization,    only : parametrization_dvs
  use airfoil_evaluation, only : get_last_design_parameters, get_last_airfoil
  use xfoil_driver,       only : run_xfoil

  double precision, dimension(:), intent(in) :: optdesign
  double precision, intent(in) :: f0, fmin

  double precision, dimension(noppoint) :: alpha, lift, drag, moment, &
                                           xtrt, xtrb
  double precision, dimension(noppoint) :: actual_flap_degrees
  double precision :: ffact, fxfact, actual_x_flap, tefact, actual_tcTE
  integer :: dvtbnd1, dvtbnd2, dvbbnd1, dvbbnd2, nmodest, nmodesb, nptt, nptb, i
  integer :: flap_idx, flap_idi, dvcounter, iunit
  type(airfoil_type) :: final_airfoil
  character(80) :: output_file, aero_file
  character(80) :: text, text2
  character(80) :: flapnote

  nmodest = nshapedvtop
  nmodesb = nshapedvbot
  nptt = size(xseedt,1)
  nptb = size(xseedb,1)

! Set modes for top and bottom surfaces
  
  dvtbnd1 = 1
  dvtbnd2 = nshapedvtop
  dvbbnd1 = dvtbnd2 + 1
  dvbbnd2 = nshapedvtop + nshapedvbot
! Overwrite lower DVs for symmetrical airfoils (they are not used)

  if (symmetrical) then
    dvbbnd1 = 1
    dvbbnd2 = dvtbnd2
  end if

  if (flap_optimization_only) then
    dvtbnd1 = 0
    dvtbnd2 = 0
    dvbbnd1 = 0
    dvbbnd2 = 0
  end if
    
  ! Get actual trailing edge based on design variable
  if (int_tcTE_spec == 1) then
    tefact = 1.d0/(max_tcTE - min_tcTE)
    actual_tcTE = optdesign(dvbbnd2+nflap_optimize+int_x_flap_spec+         &
      int_tcTE_spec)/tefact + min_tcTE
  else
    actual_tcTE=tcTE
  end if
  
  final_airfoil%npoint = nptt + nptb - 1
  
  call get_last_airfoil(restart_stat, global_search_stat, local_search_stat,   &
    final_airfoil)
  
! Use Xfoil to analyze final design

  if (.not. match_foils) then
    
!   Get actual flap angles based on design variables

    ffact = 1.d0/(max_flap_degrees - min_flap_degrees)
    actual_flap_degrees(1:noppoint) = flap_degrees(1:noppoint)
    dvcounter = dvbbnd2 + 1
    do i = 1, nflap_optimize
      flap_idx = flap_optimize_points(i)
      actual_flap_degrees(flap_idx) = optdesign(dvcounter)/ffact+min_flap_degrees
      dvcounter = dvcounter + 1
    end do
!   Set identical flap angles
    do i = 1, nflap_identical
      flap_idi = flap_identical_points(i)
      actual_flap_degrees(flap_idi) = actual_flap_degrees(flap_identical_op(flap_idi))
    end do
!   Get flap chord based on design variables
    fxfact = 1.d0/(max_flap_x - min_flap_x)
    actual_x_flap = optdesign(dvcounter)/fxfact + min_flap_x
    
!   Run xfoil for requested operating points
    !   Get lift, drag, moment, alpha, xtrt, xtrb from file

    call get_last_design_parameters(restart_stat, global_search_stat,          &
      local_search_stat, lift, drag, moment, alpha, xtrt, xtrb)

!   Write summary to screen and file

    aero_file = trim(output_prefix)//'_performance_summary.dat'
    iunit = 13
    open(unit=iunit, file=aero_file, status='replace')

    write(*,*)
    write(*,'(A)') " Optimal airfoil performance summary"
    write(iunit,'(A)') " Optimal airfoil performance summary"
    write(*,'(A)') " ----------------------------------------------------------"
    write(iunit,'(A)')                                                         &
                   " ----------------------------------------------------------"
    if (int_x_flap_spec == 1) then
      write(*,'(A,F9.5)') " Flap x hinge position: ", actual_x_flap
      write(iunit,'(A,F9.5)') " Flap x hinge position: ", actual_x_flap
      write(*,'(A,F9.5)') " Flap 1-x hinge position: ", 1.-actual_x_flap
      write(iunit,'(A,F9.5)') " Flap 1-x hinge position: ", 1.-actual_x_flap
      write(*,*)
      write(iunit,*)
    end if
    
    if (int_tcTE_spec == 1) then
      write(*,'(A,F9.5)') " Trailing Edge thickness: ", actual_tcTE
      write(iunit,'(A,F9.5)') " Trailing Edge thickness: ", actual_tcTE
      write(*,*)
      write(iunit,*)
    end if
    
    do i = 1, noppoint
      write(text,*) i
      text = adjustl(text)
      if (flap_selection(i) == "specify") then
        flapnote = " (specified)"
        flapnote = " (specified)"
      elseif (flap_selection(i) == "identical") then
        write(text2,*) flap_identical_op(i)
        text2 = adjustl(text2)
        flapnote = " (identical to "//trim(text2)//")"
      else
        flapnote = " (optimized)"
      end if
      write(*,'(A)') " Operating point "//trim(text)
      write(iunit,'(A)') " Operating point "//trim(text)
      write(*,'(A18,ES9.3)') " Reynolds number: ", reynolds(i)
      write(iunit,'(A18,ES9.3)') " Reynolds number: ", reynolds(i)
      write(*,'(A14,F9.5)') " Mach number: ", mach(i)
      write(iunit,'(A14,F9.5)') " Mach number: ", mach(i)
      write(*,'(A8,F9.5)') " ncrit: ", ncrit_pt(i)
      write(iunit,'(A8,F9.5)') " ncrit: ", ncrit_pt(i)
      if (use_flap) then
        write(*,'(A25,F9.5,A)') " Flap setting (degrees): ",                 &
                                  actual_flap_degrees(i), trim(flapnote)
        write(iunit,'(A25,F9.5,A)') " Flap setting (degrees): ",             &
                                  actual_flap_degrees(i), trim(flapnote)
      endif
      write(*,'(A18,F9.5)') " Angle of attack: ", alpha(i) 
      write(iunit,'(A18,F9.5)') " Angle of attack: ", alpha(i) 
      write(*,'(A19,F9.5)') " Lift coefficient: ", lift(i)
      write(iunit,'(A19,F6.4)') " Lift coefficient: ", lift(i)
      write(*,'(A19,F9.5)') " Drag coefficient: ", drag(i)
      write(iunit,'(A19,F9.5)') " Drag coefficient: ", drag(i)
      write(*,'(A21,F9.5)') " Moment coefficient: ", moment(i)
      write(iunit,'(A21,F9.5)') " Moment coefficient: ", moment(i)
      write(*,'(A21,F9.5)') " Top transition x/c: ", xtrt(i)
      write(iunit,'(A21,F9.5)') " Top transition x/c: ", xtrt(i)
      write(*,'(A24,F9.5)') " Bottom transition x/c: ", xtrb(i)
      write(iunit,'(A24,F9.5)') " Bottom transition x/c: ", xtrb(i)
      if (i /= noppoint) then
        write(*,*)
        write(iunit,*)
      end if
    end do

    write(*,*)
    write(*,'(A43F10.4A1)') " Objective function improvement over seed: ",      &
                           (f0 - fmin)/f0*100.d0, "%" 
    write(iunit,*)
    write(iunit,'(A43F10.4A1)') " Objective function improvement over seed: ",  &
                           (f0 - fmin)/f0*100.d0, "%" 

    close(iunit)

    write(*,*)
    write(*,*) "Optimal airfoil performance summary written to "               &
               //trim(aero_file)//"."

  end if

! Write airfoil to file

  output_file = trim(output_prefix)//'.dat'
  call airfoil_write(output_file, output_prefix, final_airfoil)

! Deallocate final airfoil

  call deallocate_airfoil(final_airfoil)

end subroutine write_final_design

end module optimization_driver
