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
!  2023-2024 Guilherme Pangas
    
module input_output

! Module with subroutines for reading and writing of files

  implicit none

#ifndef PACKAGE_VERSION
  #define PACKAGE_VERSION ""
#endif

  contains

!=============================================================================80
!
! Subroutine to read inputs from namelist file
!
!=============================================================================80
subroutine read_inputs(input_file, search_type, global_search, local_search,   &
                       seed_airfoil, nparameters_top,             &
                       nparameters_bot, restart, restart_write_freq,            &
                       naca_options, pso_options, ga_options, &
                       ds_options, matchfoil_file)

  use vardef
  use particle_swarm,     only : pso_options_type
  use genetic_algorithm,  only : ga_options_type
  use simplex_search,     only : ds_options_type
  use airfoil_operations, only : my_stop
  use airfoil_evaluation, only : xfoil_options, xfoil_geom_options
  use naca,               only : naca_options_type
  use math_deps,          only : sort_vector
  use parametrization,    only : parametrization_constrained_dvs
 
  character(*), intent(in) :: input_file
  character(80), intent(out) :: search_type, global_search, local_search,      &
                                seed_airfoil, matchfoil_file
  integer, intent(out) :: nparameters_top, nparameters_bot
  integer, dimension(max_addthickconst) :: sort_idxs
  double precision, dimension(max_addthickconst) :: temp_thickmin, temp_thickmax
  double precision, dimension(max_op_points) :: op_point_start, op_point_end,  &
                                                op_point_step
  type(naca_options_type), intent(out) :: naca_options
  type(pso_options_type), intent(out) :: pso_options
  type(ga_options_type), intent(out) :: ga_options
  type(ds_options_type), intent(out) :: ds_options

  logical :: viscous_mode, silent_mode, feasible_init,                         &
             reinitialize, restart, write_designs, reflexed
  integer :: init_number_points
  double precision :: init_al0, init_cl0, init_initial_position
  character(20) :: init_type, init_dist
  
  integer :: restart_write_freq, pso_pop, bl_maxit, npan, feasible_init_attempts
  integer :: ga_pop, simplex_maxit, ga_maxit, pso_maxit
  double precision :: maxt, xmaxt, maxc, xmaxc, design_cl, a, leidx
  double precision :: pso_tol, simplex_tol, ncrit, xtript, xtripb, vaccel
  double precision :: cvpar, cterat, ctrrat, xsref1, xsref2, xpref1, xpref2
  double precision :: feasible_limit, pso_speed_limit
  double precision :: ga_tol, parent_fraction, roulette_selection_pressure,    &
                      tournament_fraction, crossover_range_factor,             &
                      mutant_probability, chromosome_mutation_rate,            &
                      mutation_range_factor
  integer :: nmoment_constraint, nxtr_opt
  integer :: i, j, iunit, ioerr, iostat1
  character(30) :: text
  character(3) :: family
  character(10) :: pso_convergence_profile, parents_selection_method
  character :: choice
 
  namelist /optimization_options/ search_type, global_search, local_search,    &
            seed_airfoil, airfoil_file, shape_functions, nparameters_top,      &
            nparameters_bot, flap_optimization_only, abs_initial_perturb,      &
            rel_initial_perturb, penalty_limit_initial, penalty_limit_end,     &
            penalty_factor, allow_seed_penalties,  min_bump_width,             &
            kulfan_bussoletti_LEM,                                             &
            b_spline_degree, b_spline_xtype, b_spline_distribution, restart,   &
            restart_write_freq, write_designs, write_cp_file, write_bl_file,   &
            write_dvs_file, progress_per_eval, number_threads
  namelist /operating_conditions/ noppoint, op_mode, op_point, op_point_start, &
            op_point_end, op_point_step, reynolds, mach, use_flap,  x_flap,    &
            flap_connection, connection_apply, connection_radius, x_flap_spec, &
            y_flap, y_flap_spec, TE_spec, tcTE, xltTE, flap_selection,         &
            flap_identical_op, flap_degrees, weighting, optimization_type,     &
            use_previous_op, target_value, ncrit_pt, optimization_correlation
  namelist /constraints/ min_thickness, max_thickness, moment_constraint_type, &
                         min_moment, lift_constraint_type,                     &
                         min_lift, drag_constraint_type,                       &
                         max_drag, max_growth_seed_mult,                       &
                         min_leading_edge_angle, max_leading_edge_angle,       &
                         dif_leading_edge_angle, min_te_angle,                 &
                         te_angle_x_apply, lift_check_tol, drag_check_tol,     &
                         max_panel_angle, check_curvature,                     &
                         max_curv_reverse_top, max_curv_reverse_bot,           &
                         curv_threshold, symmetrical, min_flap_degrees,        &
                         max_flap_degrees, min_flap_x, max_flap_x, max_tcTE,   &
                         min_tcTE, max_camber, min_camber, naddthickconst,     &
                         addthick_x, addthick_min, addthick_max
  namelist /naca_airfoil/ family, maxt, xmaxt, maxc, xmaxc, design_cl, a,      &
                          leidx, reflexed
  namelist /initialization/ feasible_init, feasible_limit,                     &
                            feasible_init_attempts
  namelist /particle_swarm_options/ pso_pop, pso_tol, pso_maxit,               &
                                    pso_speed_limit, pso_convergence_profile
  namelist /genetic_algorithm_options/ ga_pop, ga_tol, ga_maxit,               &
            parents_selection_method, parent_fraction,                         &
            roulette_selection_pressure, tournament_fraction,                  &
            crossover_range_factor, mutant_probability,                        &
            chromosome_mutation_rate, mutation_range_factor
  namelist /simplex_options/ simplex_tol, simplex_maxit
  namelist /xfoil_run_options/ ncrit, xtript, xtripb, viscous_mode,            &
            silent_mode, bl_maxit, vaccel, reinitialize, init_type,            &
            init_number_points, init_al0, init_cl0, init_initial_position,     &
            init_dist
  namelist /xfoil_paneling_options/ npan, cvpar, cterat, ctrrat, xsref1,       &
            xsref2, xpref1, xpref2
  namelist /matchfoil_options/ match_foils, matchfoil_file
  namelist /aircraft_data/ weight_i, S_w, A_w, e_w, S_expose, thrust_coeff,      &
                           height, width, length, f_skin_roughness,            &
                           wetted_area, interference_factor, Cd_ld,            &
                           tail_config, tail_chord, t_c_ratio, max_t_x,        &
                           t_skin_roughness, tail_surface_area, add_drag
  namelist /take_off_data/ h_take_off, A_1, miu, S_g, weight_empty,            &
                           weight_payload_ref
  namelist /climb_data/ acel_to_climb, time_climb, h_climb, dh, V_0_climb,     &
                        points_coeff
  namelist /dash_data/ acel_to_dash, time_dash, h_dash, V_0_dash, time_extra,  &
                       dist_ref_dash
  namelist /turn_data/ activation_turn, h_turn, n_turn, dash_leg

  choice = 'n'

! Open input file

  iunit = 12
  open(unit=iunit, file=input_file, status='old', iostat=ioerr)
  if (ioerr /= 0)                                                              &
    call my_stop('Could not find input file '//trim(input_file)//'.')

! Set defaults for main namelist options

  search_type = 'global_and_local'
  global_search = 'particle_swarm'
  local_search = 'simplex'
  seed_airfoil = 'naca'
  shape_functions = 'hicks-henne'
  min_bump_width = 0.1d0
  kulfan_bussoletti_LEM = .false.
  b_spline_degree = 3      
  b_spline_xtype = 1       
  b_spline_distribution = 3
  nparameters_top = 4
  nparameters_bot = 4
  flap_optimization_only = .false.
  abs_initial_perturb = 0.025d0
  rel_initial_perturb = 0.0d0
  penalty_limit_initial = 1.0D-4
  penalty_limit_end = 1.0D-4
  penalty_factor = 1.0D0
  allow_seed_penalties = .false.
  restart = .false.
  restart_write_freq = 20
  write_designs = .true.
  write_cp_file = .false.
  write_bl_file = .false.
  write_dvs_file = .false.
  progress_per_eval = 'full'
  number_threads = 0

! Read main namelist options

  rewind(iunit)
  read(iunit, iostat=iostat1, nml=optimization_options)
  call namelist_check('optimization_options', iostat1, 'warn')

! Error checking and setting search algorithm options

  if (trim(search_type) /= 'global_and_local' .and. trim(search_type) /=       &
      'global' .and. trim(search_type) /= 'local')                             &
    call my_stop("search_type must be 'global_and_local', 'global', "//   &
                 "or 'local'.")

! Set defaults for operating conditions and constraints

  noppoint = 1
  use_flap = .false.
  flap_connection = 'sharp'
  connection_apply = 'none' 
  connection_radius = 0.02
  x_flap = 0.75d0
  x_flap_spec = 'specify' ! specify, optimize 
  y_flap = 0.d0
  y_flap_spec = 'y/c'
  TE_spec = 'use_seed'    !'specify', 'use_seed' or 'optimize' 
  tcTE = 1.0E-4           
  xltTE = 0.8             
  
  op_mode(:) = 'spec-cl'
  optimization_type(:) = 'min-drag'
  optimization_correlation(:) = 'none'
  use_previous_op(:) = .false.
  op_point(:) = 1.0d0
  op_point_start(:) = 1.0d0
  op_point_end(:) = 1.0d0
  op_point_step(:) = 1.0d0
  target_value(:) = 1.0d0
  reynolds(:) = 1.0D+05
  mach(:) = 0.d0
  flap_selection(:) = 'specify'
  flap_identical_op(:) = 0  
  flap_degrees(:) = 0.d0
  weighting(:) = 1.d0
  ncrit_pt(:) = -1.d0

  min_thickness = 0.06d0
  max_thickness = 1000.d0
  min_camber = -0.1d0
  max_camber = 0.1d0
  moment_constraint_type(:) = 'none'
  min_moment(:) = -1.d0
  lift_constraint_type(:) = 'none'
  min_lift(:) = -1.d0
  drag_constraint_type(:) = 'none'
  max_drag(:) = 1.d0
  max_growth_seed_mult = 2.0 
  min_leading_edge_angle = 60.00
  max_leading_edge_angle = 89.99
  dif_leading_edge_angle = 20.00
  min_te_angle = 5.00
  te_angle_x_apply = 0.5
  max_panel_angle = 25.0
  lift_check_tol = 0.2
  drag_check_tol = 0.2
  check_curvature = .false.
  max_curv_reverse_top = 1
  max_curv_reverse_bot = 1
  curv_threshold = 0.30d0
  symmetrical = .false.
  min_flap_degrees = -5.d0
  max_flap_degrees = 15.d0
  min_flap_x = 0.7d0
  max_flap_x = 0.9d0
  max_tcTE = 0.1
  min_tcTE = 0.0
  naddthickconst = 0
  addthick_x(:) = 0.01d0
  addthick_min(:) = -1000.d0
  addthick_max(:) = 1000.d0

! Read operating conditions and constraints

  rewind(iunit)
  read(iunit, iostat=iostat1, nml=operating_conditions)
  call namelist_check('operating_conditions', iostat1, 'stop')
  rewind(iunit)
  read(iunit, iostat=iostat1, nml=constraints)
  call namelist_check('constraints', iostat1, 'stop')
  
! Set defaults for aircraft_data namelist options
  weight_i = 30
  S_w = 2
  A_w = 6
  e_w = 0.9
  S_expose = 1.8
  thrust_coeff(1) = 1
  thrust_coeff(2) = 1
  thrust_coeff(3) = 1
  fuselage%height = 0.2
  fuselage%width = 0.2
  fuselage%length = 0.8
  fuselage%skin_roughness = 0.1
  fuselage%wetted_area = 1
  fuselage%interference_factor = 1.2
  Cd_ld = 0.03
  add_drag = 0
  tail%config = 1
  tail%chord(1) = 0.15
  tail%t_c_ratio(1) = 0.5
  tail%max_t_x(1) = 0.12
  tail%skin_roughness(1) = 0.3
  tail%surface_area(1) = 0.1
  
! Set defaults for take_off_data namelist options
  h_take_off = 200
  A_1 = 1.2
  miu = 0.03
  S_g = 60
  weight_empty = 10
  weight_payload_ref = 50
  
! Set defaults for climb_data namelist options
  acel_to_climb = .false.
  time_climb = 60
  h_climb = 300
  dh = 100
  V_0_climb = 0 
  points_coeff(1) = -3.92E-5
  points_coeff(2) = 1.08E-2
  points_coeff(3) = -1.156
  points_coeff(4) = 64.2
  points_coeff(5) = -537
 
! Set defaults for dash_data namelist options
  acel_to_dash = .true.         
  time_dash = 60
  h_dash = 300
  V_0_dash = 0
  time_extra = 0
  
! Set defaults for turn_data namelist options
  activation_turn = .false.
  h_turn = 300
  n_turn = 2
  dash_leg = 300
  
!read aircraft_data, take_off_data, climb_data and dash_data    
  
  rewind(iunit)
  read(iunit, iostat=iostat1, nml=aircraft_data)
  call namelist_check('aircraft_data', iostat1, 'warn')   
  rewind(iunit)
  read(iunit, iostat=iostat1, nml=take_off_data)
  call namelist_check('take_off_data', iostat1, 'warn')
  rewind(iunit)
  read(iunit, iostat=iostat1, nml=climb_data)
  call namelist_check('climb_data', iostat1, 'warn')
  rewind(iunit)
  read(iunit, iostat=iostat1, nml=dash_data)
  call namelist_check('dash_data', iostat1, 'warn')
  rewind(iunit)
  read(iunit, iostat=iostat1, nml=turn_data)
  call namelist_check('turn_data', iostat1, 'warn')
  
!Aircraft_data in derived type
  weight = weight_i
  fuselage%height = height
  fuselage%width = width
  fuselage%length = length
  fuselage%skin_roughness = f_skin_roughness
  fuselage%wetted_area = wetted_area
  fuselage%interference_factor = interference_factor
    
  tail%config = tail_config
  do i = 1, tail%config
      tail%chord(i) = tail_chord(i)
      tail%t_c_ratio(i) = t_c_ratio(i)
      tail%max_t_x(i) = max_t_x(i)
      tail%skin_roughness(i) =  t_skin_roughness(i)
      tail%surface_area(i) = tail_surface_area(i)
  end do
    
!take_off_data in derived type
  take_off%h = h_climb
  take_off%A_1 = A_1
  take_off%miu = miu
  take_off%S_g = S_g
  take_off%weight_empty = weight_empty
  take_off%weight_payload_ref = weight_payload_ref
  take_off%CL_max = 0
  take_off%CL_run = 0
  take_off%CD_run = 0
  take_off%V_to = 0
  take_off%V_run = 0
  take_off%points = 0
    
!climb_data in derived type
  climb%acel = acel_to_climb
  climb%time = time_climb
  climb%h = h_climb
  climb%dh = dh
  climb%V_0 = V_0_climb
  do i=1, 5
    climb%points_coeff(i) = points_coeff(i)
  end do
  climb%RC_max = 0
  climb%V = 0
  climb%Cl = 0
  climb%D = 0
  climb%T = 0
  climb%t_acel = 0
  climb%points = 0
  
!dash_data in derived type
  dash%acel = acel_to_dash
  dash%time = time_dash
  dash%h = h_dash
  dash%dist_ref = dist_ref_dash
  dash%V_0 = V_0_dash
  dash%t_ex = time_extra
  dash%V_max = 0
  dash%t_acel_d = 0 
  dash%dist_acel_d = 0 
  dash%dist = 0
  dash%points = 0
  
!turn_data in derived type 
  turn%activation = activation_turn
  turn%h = h_turn
  turn%n = n_turn
  turn%dash_leg = dash_leg
  turn%V= 0
  turn%radius = 0
  
! use_previous_op set all op options to previous if true
  
  do i=2, noppoint
    if(use_previous_op(i)) then
      !optimization_type stays the same
      op_mode(i) = op_mode(i-1)
      op_point(i) = op_point(i-1)
      op_point_start(i) = op_point_start(i-1)
      op_point_end(i) = op_point_end(i-1)
      op_point_step(i) = op_point_step(i-1)
      target_value(i) = target_value(i-1)
      reynolds(i) = reynolds(i-1)
      mach(i) = mach(i-1)
      if(flap_selection(i-1) .EQ. 'identical') then
        flap_selection(i) = flap_selection(i-1)
        flap_identical_op(i) = flap_identical_op(i-1)
      else
        flap_selection(i) = 'identical'
        flap_identical_op(i) = i-1
      end if
      flap_degrees(i) = flap_degrees(i-1)
      !weighting stays the same
      ncrit_pt(i) = ncrit_pt(i-1)
    end if
  end do
  
! Set Search variable   
  
  op_search%noppoint = 0
  do i=1, noppoint
    if ((optimization_type(i) .eq. 'max-lift-search') .or.                     &
        ((optimization_type(i) .eq. 'take-off') .and.                          &
        (optimization_correlation(i) .eq. 'init'))) then
      op_search%noppoint = op_search%noppoint + 1
    end if
  end do
  
  if (op_search%noppoint .ne. 0) then
    
    allocate(op_search%oppoints(op_search%noppoint))
    allocate(op_search%op_start(op_search%noppoint))
    allocate(op_search%op_end(op_search%noppoint))
    allocate(op_search%op_step(op_search%noppoint))
  
    j=1
    do i=1, noppoint
      if ((optimization_type(i) .eq. 'max-lift-search') .or.                     &
        ((optimization_type(i) .eq. 'take-off') .and.                          &
        (optimization_correlation(i) .eq. 'init'))) then
        op_search%oppoints(j) = i
        op_search%op_start(j) = op_point_start(i)
        op_search%op_end(j) = op_point_end(i)
        op_search%op_step(j) = op_point_step(i)
        j=j+1
      end if
    end do
    
  end if
    
  
! Normalize weightings for operating points

  weighting = weighting/sum(weighting(1:noppoint))

! Ask about removing pitching moment constraints for symmetrical optimization

  if (symmetrical) then
    nmoment_constraint = 0
    do i = 1, noppoint
      if (trim(moment_constraint_type(i)) /= 'none')                           &
        nmoment_constraint = nmoment_constraint + 1
    end do
    if (nmoment_constraint > 0) choice = ask_moment_constraints()
    if (choice == 'y') moment_constraint_type(:) = 'none'
  end if

! Sort thickness constraints in ascending x/c order

  if (naddthickconst > 0) then
    call sort_vector(addthick_x(1:naddthickconst), sort_idxs(1:naddthickconst))
    temp_thickmin = addthick_min
    temp_thickmax = addthick_max
    do i = 1, naddthickconst
      addthick_min(i) = temp_thickmin(sort_idxs(i))
      addthick_max(i) = temp_thickmax(sort_idxs(i))
    end do
  end if

! Set defaults for naca airfoil options
 
  family = '4'
  maxt = 0.1d0
  xmaxt = 0.3d0
  maxc = 0.d0
  xmaxc = 0.3d0
  design_cl = 0.3d0
  a = 1.d0
  leidx = 6.d0
  reflexed = .false.

! Read naca airfoil options and put them into derived type

  if ( (seed_airfoil == 'naca') .or. (seed_airfoil == 'NACA') .or.             &
       (seed_airfoil == 'Naca') ) then
    rewind(iunit)
    read(iunit, iostat=iostat1, nml=naca_airfoil)
    call namelist_check('naca_airfoil', iostat1, 'warn')
    
    airfoil_file = 'naca_family.dat'
    naca_options%family = family
    naca_options%maxt = maxt
    naca_options%xmaxt = xmaxt
    naca_options%maxc = maxc
    naca_options%xmaxc = xmaxc
    naca_options%design_cl = design_cl
    naca_options%a = a
    naca_options%leidx = leidx
    naca_options%reflexed = reflexed
  end if

! Set default initialization options

  feasible_init = .true.
  feasible_limit = 5.0D+04
  feasible_init_attempts = 1000

! Read initialization parameters

  rewind(iunit)
  read(iunit, iostat=iostat1, nml=initialization)
  call namelist_check('initialization', iostat1, 'warn')

! Set default particle swarm options

  pso_pop = 40
  pso_tol = 1.D-04
  pso_maxit = 700
  pso_speed_limit = 0.025d0
  pso_convergence_profile = 'exhaustive'

! Set default genetic algorithm options

  ga_pop = 80
  ga_tol = 1.D-04
  ga_maxit = 700
  parents_selection_method = 'tournament'
  parent_fraction = 0.5d0
  roulette_selection_pressure = 8.d0
  tournament_fraction = 0.025d0
  crossover_range_factor = 0.5d0
  mutant_probability = 0.4d0
  chromosome_mutation_rate = 0.01d0
  mutation_range_factor = 0.2d0

! Set default simplex search options

  simplex_tol = 1.0D-05
  simplex_maxit = 1000
  
  epsexit_linear = .true. 
  
  ! Read only the required namelist
  if (trim(search_type) == 'global_and_local' .or. trim(search_type) ==        &
      'global') then

    if (trim(global_search) == 'particle_swarm') then

!     Read PSO options and put them into derived type

      rewind(iunit)
      read(iunit, iostat=iostat1, nml=particle_swarm_options)
      call namelist_check('particle_swarm_options', iostat1, 'warn')
      pso_options%pop = pso_pop
      pso_options%tol = pso_tol
      pso_options%maxspeed = pso_speed_limit
      pso_options%maxit = pso_maxit
      pso_options%convergence_profile = pso_convergence_profile
      pso_options%feasible_init = feasible_init
      pso_options%feasible_limit = feasible_limit
      pso_options%feasible_init_attempts = feasible_init_attempts
      pso_options%write_designs = write_designs
      if (.not. match_foils) then
        pso_options%relative_fmin_report = .true.
      else
        pso_options%relative_fmin_report = .false.
      end if
      maxit = pso_maxit
    else if (trim(global_search) == 'genetic_algorithm') then

!     Read genetic algorithm options and put them into derived type

      rewind(iunit)
      read(iunit, iostat=iostat1, nml=genetic_algorithm_options)
      call namelist_check('genetic_algorithm_options', iostat1, 'warn')
      ga_options%pop = ga_pop
      ga_options%tol = ga_tol
      ga_options%maxit = ga_maxit
      ga_options%parents_selection_method = parents_selection_method
      ga_options%parent_fraction = parent_fraction
      ga_options%roulette_selection_pressure = roulette_selection_pressure
      ga_options%tournament_fraction = tournament_fraction
      ga_options%crossover_range_factor = crossover_range_factor
      ga_options%mutant_probability = mutant_probability
      ga_options%chromosome_mutation_rate = chromosome_mutation_rate
      ga_options%mutation_range_factor = mutation_range_factor
      ga_options%feasible_init = feasible_init
      ga_options%feasible_limit = feasible_limit
      ga_options%feasible_init_attempts = feasible_init_attempts
      ga_options%write_designs = write_designs
      if (.not. match_foils) then
        ga_options%relative_fmin_report = .true.
      else
        ga_options%relative_fmin_report = .false.
      end if
      maxit = ga_options%maxit
    else
      call my_stop("Global search type '"//trim(global_search)//               &
                   "' is not available.")
    end if
  end if

  if (trim(search_type) == 'global_and_local' .or. trim(search_type) ==        &
      'local') then

    if (trim(local_search) == 'simplex') then

!     Read simplex search options and put them into derived type

      rewind(iunit)
      read(iunit, iostat=iostat1, nml=simplex_options)
      call namelist_check('simplex_options', iostat1, 'warn')
      ds_options%tol = simplex_tol
      ds_options%maxit = simplex_maxit
      ds_options%write_designs = write_designs
      if (.not. match_foils) then
        ds_options%relative_fmin_report = .true.
      else
        ds_options%relative_fmin_report = .false.
      end if
      if (trim(search_type) == 'local') maxit = ds_options%maxit
    else
      call my_stop("Local search type '"//trim(local_search)//   &
                   "' is not available.")
    end if

  end if 

! Set default xfoil aerodynamics and paneling options

  ncrit = 9.d0
  xtript = 1.d0
  xtripb = 1.d0
  viscous_mode = .true.
  silent_mode = .true.
  bl_maxit = 100
  vaccel = 0.01d0
  reinitialize = .true.
  init_type = 'unconverged'
  init_number_points = 5
  init_al0 = 0.0
  init_cl0 = 0.5
  init_initial_position = 0.5
  init_dist = 'linear'
  
  npan = 160
  cvpar = 1.d0
  cterat = 0.15d0
  ctrrat = 0.2d0
  xsref1 = 1.d0
  xsref2 = 1.d0
  xpref1 = 1.d0
  xpref2 = 1.d0

! Read xfoil options

  rewind(iunit)
  read(iunit, iostat=iostat1, nml=xfoil_run_options)
  call namelist_check('xfoil_run_options', iostat1, 'warn')
  rewind(iunit)
  read(iunit, iostat=iostat1, nml=xfoil_paneling_options)
  call namelist_check('xfoil_paneling_options', iostat1, 'warn')

  ! Set trasition points according to flap_connection 

  if (use_flap .and. (trim(flap_connection) .NE. 'smooth') .and.     &
      ( (trim(connection_apply) .EQ. 'trip_wire') .OR.               &
        (trim(connection_apply) .EQ. 'both'     ) )      ) then
    if (trim(flap_connection) .EQ. 'smooth-top') then
      xtripb=x_flap
    elseif (trim(flap_connection) .EQ. 'smooth-bot') then
      xtript=x_flap
    else
      xtript=x_flap
      xtripb=x_flap
    end if
  end if
 
  ! Ask about removing turbulent trips for max-xtr optimization
  
  nxtr_opt = 0
  if ( (xtript < 1.d0) .or. (xtripb < 1.d0) )then
    do i = 1, noppoint
      if (trim(optimization_type(i)) == "max-xtr") nxtr_opt = nxtr_opt + 1
    end do
 
    if (nxtr_opt > 0) choice = ask_forced_transition()
    if (choice == 'y') then
      xtript = 1.d0
      xtripb = 1.d0
    end if
  end if

! Put xfoil options into derived types

  xfoil_options%ncrit = ncrit
  xfoil_options%xtript = xtript
  xfoil_options%xtripb = xtripb
  xfoil_options%viscous_mode = viscous_mode
  xfoil_options%silent_mode = silent_mode
  xfoil_options%maxit = bl_maxit
  xfoil_options%vaccel = vaccel
  xfoil_options%reinitialize = reinitialize
  xfoil_options%init_type = init_type
  xfoil_options%init_number_points = init_number_points
  xfoil_options%init_al0 = init_al0
  xfoil_options%init_cl0 = init_cl0
  xfoil_options%init_initial_position = init_initial_position
  xfoil_options%init_dist = init_dist
  

  xfoil_geom_options%npan = npan
  xfoil_geom_options%cvpar = cvpar
  xfoil_geom_options%cterat = cterat
  xfoil_geom_options%ctrrat = ctrrat
  xfoil_geom_options%xsref1 = xsref1
  xfoil_geom_options%xsref2 = xsref2
  xfoil_geom_options%xpref1 = xpref1
  xfoil_geom_options%xpref2 = xpref2

! Set per-point ncrit if not specified in namelist

  do i = 1, noppoint
    if (ncrit_pt(i) == -1.d0) ncrit_pt(i) = ncrit
  end do

! Option to match seed airfoil to another instead of aerodynamic optimization

  match_foils = .false.
  matchfoil_file = 'none'
  read(iunit, iostat=iostat1, nml=matchfoil_options)
  call namelist_check('matchfoil_options', iostat1, 'warn')

! Close the input file

  close(iunit)
  
  ! Avaliate int_kulfan_bussoletti_LEM  
  int_kulfan_bussoletti_LEM = 0
  if (kulfan_bussoletti_LEM) int_kulfan_bussoletti_LEM = 1
  
! Store number and identification of operating points restricted by moment 
  nmoment_constrain = 0
  if (.not. match_foils) then
    do i = 1, noppoint
      if (trim(moment_constraint_type(i)) /= 'none') then
        nmoment_constrain = nmoment_constrain + 1
        moment_constrain_points(nmoment_constrain) = i
      end if
    end do
  end if
  
! Store number and identification of operating points restricted by drag 
  ndrag_constrain = 0
  if (.not. match_foils) then
    do i = 1, noppoint
      if (trim(drag_constraint_type(i)) /= 'none') then
        ndrag_constrain = ndrag_constrain + 1
        drag_constrain_points(ndrag_constrain) = i
      end if
    end do
  end if
  
! Store number and identification of operating points restricted by lift 
  nlift_constrain = 0
  if (.not. match_foils) then
    do i = 1, noppoint
      if (trim(lift_constraint_type(i)) /= 'none') then
        nlift_constrain = nlift_constrain + 1
        lift_constrain_points(nlift_constrain) = i
      end if
    end do
  end if
  
! Store operating points where flap setting will be optimized
  nflap_optimize = 0
  if ((use_flap) .and. (.not. match_foils)) then
    do i = 1, noppoint
      if (flap_selection(i) == 'optimize') then
        nflap_optimize = nflap_optimize + 1
        flap_optimize_points(nflap_optimize) = i
      end if
    end do
  end if

! Store operating points where flap setting will be identical
  nflap_identical = 0
  if (use_flap .and. (.not. match_foils)) then
    do i = 1, noppoint
      if (flap_selection(i) == 'identical') then
        nflap_identical = nflap_identical + 1
        flap_identical_points(nflap_identical) = i
      end if
    end do
  end if
  
  ! Avaliate type of flap chord
  if ((trim(x_flap_spec) == 'optimize') .and. (use_flap) .and.                 &
    (.not. match_foils) .and. (nflap_optimize/=0)) then
    int_x_flap_spec = 1
  else
    int_x_flap_spec = 0
  end if  
  
! Avaliate type of trailing edge
  if (trim(TE_spec) == 'optimize') then
    int_tcTE_spec = 1
  else
    int_tcTE_spec = 0
  end if  
  
  if (.not. match_foils) then
    contrain_number = 1+nflap_optimize+3+naddthickconst+10+                    &
                               nmoment_constrain+nlift_constrain+ndrag_constrain
  else
    contrain_number = 1
  end if
  
! Echo namelist options for checking purposes

  write(*,*)
  write(*,*) 'Echoing program options:'
  write(*,*)

! Optimization options namelist

  write(*,'(A)') " &optimization_options"
  write(*,*) " search_type = '"//trim(search_type)//"'"
  write(*,*) " global_search = '"//trim(global_search)//"'"
  write(*,*) " local_search = '"//trim(local_search)//"'"
  write(*,*) " seed_airfoil = '"//trim(seed_airfoil)//"'"
  write(*,*) " airfoil_file = '"//trim(airfoil_file)//"'"
  write(*,*) " shape_functions = '"//trim(shape_functions)//"'"
  write(*,*) " min_bump_width = ", min_bump_width
  write(*,*) " kulfan_bussoletti_LEM = ", kulfan_bussoletti_LEM
  write(*,*) " b_spline_degree = ", b_spline_degree
  write(*,*) " b_spline_xtype = ", b_spline_xtype
  write(*,*) " b_spline_distribution = ", b_spline_distribution
  write(*,*) " nparameters_top = ", nparameters_top
  write(*,*) " nparameters_bot = ", nparameters_bot
  write(*,*) " flap_optimization_only = ", flap_optimization_only  
  write(*,*) " abs_initial_perturb = ", abs_initial_perturb
  write(*,*) " rel_initial_perturb = ", rel_initial_perturb
  write(*,*) " penalty_limit_initial = ", penalty_limit_initial
  write(*,*) " penalty_limit_end = ", penalty_limit_end
  write(*,*) " penalty_factor = ", penalty_factor
  write(*,*) " allow_seed_penalties = ", allow_seed_penalties
  write(*,*) " restart = ", restart
  write(*,*) " restart_write_freq = ", restart_write_freq
  write(*,*) " write_designs = ", write_designs
  write(*,*) " write_cp_file = ", write_cp_file
  write(*,*) " write_bl_file = ", write_bl_file
  write(*,*) " write_dvs_file = ", write_dvs_file
  write(*,*) " progress_per_eval = ", "'"//trim(progress_per_eval)//"'"
  write(*,*) " number_threads = ", number_threads
  write(*,'(A)') " /"
  write(*,*)

! Operating conditions namelist

  write(*,'(A)') " &operating_conditions"
  write(*,*) " noppoint = ", noppoint
  write(*,*) " use_flap = ", use_flap
  write(*,*) " flap_connection = ", flap_connection
  write(*,*) " connection_apply = ", connection_apply
  write(*,*) " connection_radius = ", connection_radius
  write(*,*) " x_flap = ", x_flap
  write(*,*) " x_flap_spec = "//trim(x_flap_spec)
  write(*,*) " y_flap = ", y_flap
  write(*,*) " y_flap_spec = "//trim(y_flap_spec)
  write(*,*) " TE_spec = "//trim(TE_spec)
  write(*,*) " tcTE = ", tcTE
  write(*,*) " xltTE = ", xltTE
  write(*,*)
  do i = 1, noppoint
    write(text,*) i
    text = adjustl(text)
    write(*,*) " optimization_type("//trim(text)//") = '"//                    &
               trim(optimization_type(i))//"'"
    write(*,*) " optimization_correlation("//trim(text)//") = '"//             &
               trim(optimization_correlation(i))//"'"
    write(*,*) " op_mode("//trim(text)//") = '"//trim(op_mode(i))//"'"
    write(*,*) " op_point("//trim(text)//") = ", op_point(i)
    write(*,*) " use_previous_op("//trim(text)//") = ", use_previous_op(i)
    write(*,*) " op_point_start("//trim(text)//") = ", op_point_start(i)
    write(*,*) " op_point_end("//trim(text)//") = ", op_point_end(i)
    write(*,*) " op_point_step("//trim(text)//") = ", op_point_step(i)
    write(*,*) " target_value("//trim(text)//") = ", target_value(i)
    write(*,'(A,es17.8)') "  reynolds("//trim(text)//") = ", reynolds(i)
    write(*,*) " mach("//trim(text)//") = ", mach(i)
    write(*,*) " flap_selection("//trim(text)//") = '"//                       &
               trim(flap_selection(i))//"'"
    write(*,*) " flap_identical_op("//trim(text)//") = ", flap_identical_op(i)
    write(*,*) " flap_degrees("//trim(text)//") = ", flap_degrees(i)
    write(*,*) " weighting("//trim(text)//") = ", weighting(i)
    if (ncrit_pt(i) /= -1.d0)                                                  &
      write(*,*) " ncrit_pt("//trim(text)//") = ", ncrit_pt(i)
    if (i < noppoint) write(*,*)
  end do
  write(*,'(A)') " /"
  write(*,*)

! Constraints namelist

  write(*,'(A)') " &constraints"
  write(*,*) " min_thickness = ", min_thickness
  write(*,*) " max_thickness = ", max_thickness
  write(*,*) " max_growth_seed_mult = ", max_growth_seed_mult
  write(*,*) " min_leading_edge_angle = ", min_leading_edge_angle
  write(*,*) " max_leading_edge_angle = ", max_leading_edge_angle
  write(*,*) " dif_leading_edge_angle = ", dif_leading_edge_angle
  write(*,*) " min_te_angle = ", min_te_angle
  write(*,*) " te_angle_x_apply = ", te_angle_x_apply
  write(*,*) " max_panel_angle = ", max_panel_angle
  write(*,*) " lift_check_tol = ", lift_check_tol
  write(*,*) " drag_check_tol = ", drag_check_tol
  write(*,*) " check_curvature = ", check_curvature
  write(*,*) " max_curv_reverse_top = ", max_curv_reverse_top
  write(*,*) " max_curv_reverse_bot = ", max_curv_reverse_bot
  write(*,*) " curv_threshold = ", curv_threshold
  write(*,*) " symmetrical = ", symmetrical
  write(*,*) " min_flap_degrees = ", min_flap_degrees
  write(*,*) " max_flap_degrees = ", max_flap_degrees
  write(*,*) " min_flap_x = ", min_flap_x
  write(*,*) " max_flap_x = ", max_flap_x
  write(*,*) " min_tcTE = ", min_tcTE
  write(*,*) " max_tcTE = ", max_tcTE
  write(*,*) " min_camber = ", min_camber
  write(*,*) " max_camber = ", max_camber
  write(*,*)
  do i = 1, noppoint
    write(text,*) i
    text = adjustl(text)
    write(*,*) " moment_constraint_type("//trim(text)//") = "//                &
               trim(moment_constraint_type(i))
    write(*,*) " min_moment("//trim(text)//") = ", min_moment(i)
  write(*,*) " lift_constraint_type("//trim(text)//") = "//                &
               trim(lift_constraint_type(i))
    write(*,*) " min_lift("//trim(text)//") = ", min_lift(i)
  write(*,*) " drag_constraint_type("//trim(text)//") = "//                &
               trim(drag_constraint_type(i))
    write(*,*) " max_drag("//trim(text)//") = ", max_drag(i)
  end do
  write(*,*)
  write(*,*) " naddthickconst = ", naddthickconst
  do i = 1, naddthickconst
    write(text,*) i
    text = adjustl(text)
    write(*,*) " addthick_x("//trim(text)//") = ", addthick_x(i)
    write(*,*) " addthick_min("//trim(text)//") = ", addthick_min(i)
    write(*,*) " addthick_max("//trim(text)//") = ", addthick_max(i)
  end do
  write(*,'(A)') " /"
  write(*,*)
  
! Aircraft_data namelist
  
  write(*,'(A)') "&aircraft_data"
  write(*,*) " weight = ", weight_i
  write(*,*) " S_w = ", S_w
  write(*,*) " A_w = ", A_w
  write(*,*) " e_w = ", e_w
  write(*,*) " S_expose = ", S_expose
  write(*,*) " thrust_coeff(1) = ", thrust_coeff(1)
  write(*,*) " thrust_coeff(2) = ", thrust_coeff(2)
  write(*,*) " thrust_coeff(3) = ", thrust_coeff(3)
  write(*,*) " fuselage%height = ", fuselage%height
  write(*,*) " fuselage%width = ", fuselage%width
  write(*,*) " fuselage%length = ", fuselage%length
  write(*,*) " fuselage%skin_roughness = ", fuselage%skin_roughness
  write(*,*) " fuselage%wetted_area = ", fuselage%wetted_area
  write(*,*) " fuselage%interderence_factor = ", fuselage%interference_factor
  write(*,*) " Cd_ld = ", Cd_ld
  write(*,*) " tail%tail_config = ", tail%config
  write(*,*) " add_drag = ", add_drag
  write(*,*)
  do i = 1, tail%config
      write(text,*) i
      text = adjustl(text)
      write(*,*) " tail%chord("//trim(text)//") = ", tail%chord(i)
      write(*,*) " tail%t_c_ratio("//trim(text)//") = ", tail%t_c_ratio(i)
      write(*,*) " tail%max_t_x("//trim(text)//") = ", tail%max_t_x(i)
      write(*,*) " tail%skin_roughness("//trim(text)//") = ", tail%skin_roughness(i)
      write(*,*) " tail%surface_area("//trim(text)//") = ", tail%surface_area(i)
      if (i < tail%config) write(*,*)
  end do
  write(*,'(A)') "/"
  write(*,*)
  
! Take_off_data namelist
    
  write(*,'(A)') "&take_off_data"
  write(*,*) " h = ", take_off%h
  write(*,*) " A_1 = ", take_off%A_1
  write(*,*) " miu = ", take_off%miu
  write(*,*) " S_g = ", take_off%S_g
  write(*,*) " weight_empty = ", take_off%weight_empty
  write(*,*) " weight_payload_ref = ", take_off%weight_payload_ref
  write(*,'(A)') "/"
  write(*,*)
  
! Climb_data namelist
   
  write(*,'(A)') "&climb_data"
  write(*,*) " acel_to_climb = ", climb%acel
  write(*,*) " time = ", climb%time
  write(*,*) " h = ", climb%h
  write(*,*) " dh = ", climb%dh
  write(*,*) " V_0 = ", climb%V_0
  do i=1, 5
      write(text,*) i
      text = adjustl(text)
      write(*,*) " points_coeff("//trim(text)//") = ", climb%points_coeff(i)
  end do
  write(*,'(A)') "/"
  write(*,*)
  
! Dash_data namelist
   
  write(*,'(A)') "&dash_data"
  write(*,*) " acel = ", dash%acel
  write(*,*) " time = ", dash%time
  write(*,*) " h = ", dash%h
  write(*,*) " dist_ref = ", dash%dist_ref
  write(*,*) " V_0 = ", dash%V_0
  write(*,*) " t_ex = ", dash%t_ex
  write(*,'(A)') "/"
  write(*,*)
  
  ! Turn_data namelist
   
  write(*,'(A)') "&Turn_data"
  write(*,*) " activation = ", turn%activation
  write(*,*) " h = ", turn%h
  write(*,*) " n = ", turn%n
  write(*,*) " dash_leg = ", turn%dash_leg
  write(*,'(A)') "/"
  write(*,*)

! naca_airfoil namelist

  write(*,'(A)') " &naca_airfoil"
  write(*,*) " family = "//trim(adjustl(family))
  write(*,*) " maxt = ", maxt
  write(*,*) " xmaxt = ", xmaxt
  write(*,*) " maxc = ", maxc
  write(*,*) " xmaxc = ", xmaxc
  write(*,*) " design_cl = ", design_cl
  write(*,*) " a = ", a
  write(*,*) " leidx = ", leidx
  write(*,*) " reflexed = ", reflexed
  write(*,'(A)') " /"
  write(*,*)

! Initialization namelist

  write(*,'(A)') " &initialization"
  write(*,*) " feasible_init = ", feasible_init
  write(*,*) " feasible_limit = ", feasible_limit
  write(*,*) " feasible_init_attempts = ", feasible_init_attempts
  write(*,'(A)') " /"
  write(*,*)

! Optimizer namelists
! Only the used namelists are written

  if (trim(search_type) == 'global_and_local' .or. trim(search_type) ==        &
      'global') then

    if (trim(global_search) == 'particle_swarm') then

!     Particle swarm namelist

      write(*,'(A)') " &particle_swarm_options"
      write(*,*) " pso_pop = ", pso_options%pop
      write(*,*) " pso_tol = ", pso_options%tol
      write(*,*) " pso_maxit = ", pso_options%maxit
      write(*,*) " pso_speed_limit = ", pso_options%maxspeed
      write(*,*) " pso_convergence_profile = ", pso_options%convergence_profile
      write(*,'(A)') " /"
      write(*,*)

    else if (trim(global_search) == 'genetic_algorithm') then

!     Genetic algorithm options

      write(*,'(A)') " &genetic_algorithm_options"
      write(*,*) " ga_pop = ", ga_options%pop
      write(*,*) " ga_tol = ", ga_options%tol
      write(*,*) " ga_maxit = ", ga_options%maxit
      write(*,*) " parents_selection_method = ",                               &
                 parents_selection_method
      write(*,*) " parent_fraction = ", ga_options%parent_fraction 
      write(*,*) " roulette_selection_pressure = ",                            &
                 ga_options%roulette_selection_pressure
      write(*,*) " tournament_fraction = " , ga_options%tournament_fraction
      write(*,*) " crossover_range_factor = ", ga_options%crossover_range_factor
      write(*,*) " mutant_probability = ", ga_options%mutant_probability
      write(*,*) " chromosome_mutation_rate = ",                               &
                 ga_options%chromosome_mutation_rate
      write(*,*) " mutation_range_factor = ", ga_options%mutation_range_factor
      write(*,'(A)') " /"
      write(*,*)

    end if

  end if

  if (trim(search_type) == 'global_and_local' .or. trim(search_type) ==        &
      'local') then

    if(trim(local_search) == 'simplex') then

!     Simplex search namelist

      write(*,'(A)') " &simplex_options"
      write(*,*) " simplex_tol = ", ds_options%tol
      write(*,*) " simplex_maxit = ", ds_options%maxit
      write(*,'(A)') " /"
      write(*,*)

    end if

  end if

! Xfoil run options namelist

  write(*,'(A)') " &xfoil_run_options"
  write(*,*) " ncrit = ", xfoil_options%ncrit
  write(*,*) " xtript = ", xfoil_options%xtript
  write(*,*) " xtripb = ", xfoil_options%xtripb
  write(*,*) " viscous_mode = ", xfoil_options%viscous_mode
  write(*,*) " silent_mode = ", xfoil_options%silent_mode
  write(*,*) " bl_maxit = ", xfoil_options%maxit
  write(*,*) " vaccel = ", xfoil_options%vaccel
  write(*,*) " reinitialize = ", xfoil_options%reinitialize
  write(*,*) " init_type = ", xfoil_options%init_type
  write(*,*) " init_number_points = ", xfoil_options%init_number_points
  write(*,*) " init_al0 = ", xfoil_options%init_al0
  write(*,*) " init_cl0 = ", xfoil_options%init_cl0
  write(*,*) " init_initial_position = ", xfoil_options%init_initial_position
  write(*,*) " init_dist = ", xfoil_options%init_dist
  write(*,'(A)') " /"
  write(*,*)

! Xfoil paneling options namelist

  write(*,'(A)') " &xfoil_paneling_options"
  write(*,*) " npan = ", xfoil_geom_options%npan
  write(*,*) " cvpar = ", xfoil_geom_options%cvpar
  write(*,*) " cterat = ", xfoil_geom_options%cterat
  write(*,*) " ctrrat = ", xfoil_geom_options%ctrrat
  write(*,*) " xsref1 = ", xfoil_geom_options%xsref1
  write(*,*) " xsref2 = ", xfoil_geom_options%xsref2
  write(*,*) " xpref1 = ", xfoil_geom_options%xpref1
  write(*,*) " xpref2 = ", xfoil_geom_options%xpref2
  write(*,'(A)') " /"
  write(*,*)

! Matchfoil options

  write(*,'(A)') " &matchfoil_options"
  write(*,*) " match_foils = ", match_foils
  write(*,*) " matchfoil_file = '"//trim(matchfoil_file)//"'"
  write(*,'(A)') " /"
  write(*,*)

! Echo namelist options for checking purposes
  open(unit=100, file='echo_'//trim(output_prefix)//".txt", status='replace') 
  
  write(100,*) '! Echoing program options:'
  write(100,*)

! Optimization options namelist

  write(100,'(A)') "&optimization_options"
  write(100,*) " search_type = '"//trim(search_type)//"'"
  write(100,*) " global_search = '"//trim(global_search)//"'"
  write(100,*) " local_search = '"//trim(local_search)//"'"
  write(100,*) " seed_airfoil = '"//trim(seed_airfoil)//"'"
  write(100,*) " airfoil_file = '"//trim(airfoil_file)//"'"
  write(100,*) " shape_functions = '"//trim(shape_functions)//"'"
  write(100,*) " min_bump_width = ", min_bump_width
  write(100,*) " kulfan_bussoletti_LEM = ", kulfan_bussoletti_LEM
  write(100,*) " b_spline_degree = ", b_spline_degree
  write(100,*) " b_spline_xtype = ", b_spline_xtype
  write(100,*) " b_spline_distribution = ", b_spline_distribution
  write(100,*) " nparameters_top = ", nparameters_top
  write(100,*) " nparameters_bot = ", nparameters_bot
  write(100,*) " flap_optimization_only = ", flap_optimization_only  
  write(100,*) " abs_initial_perturb = ", abs_initial_perturb
  write(100,*) " rel_initial_perturb = ", rel_initial_perturb
  write(100,*) " penalty_limit_initial = ", penalty_limit_initial
  write(100,*) " penalty_limit_end = ", penalty_limit_end
  write(100,*) " penalty_factor = ", penalty_factor
  write(100,*) " allow_seed_penalties = ", allow_seed_penalties
  write(100,*) " restart = ", restart
  write(100,*) " restart_write_freq = ", restart_write_freq
  write(100,*) " write_designs = ", write_designs
  write(100,*) " write_cp_file = ", write_cp_file
  write(100,*) " write_bl_file = ", write_bl_file
  write(100,*) " write_dvs_file = ", write_dvs_file
  write(100,*) " progress_per_eval = ", "'"//trim(progress_per_eval)//"'"
  write(100,*) " number_threads = ", number_threads
  write(100,'(A)') "/"
  write(100,*)

! Operating conditions namelist

  write(100,'(A)') "&operating_conditions"
  write(100,*) " noppoint = ", noppoint
  write(100,*) " use_flap = ", use_flap
  write(100,*) " flap_connection = '"//trim(flap_connection)//"'"
  write(100,*) " connection_apply = '"//trim(connection_apply)//"'"
  write(100,*) " connection_radius = ", connection_radius
  write(100,*) " x_flap = ", x_flap
  write(100,*) " x_flap_spec = '"//trim(x_flap_spec)//"'"
  write(100,*) " y_flap = ", y_flap
  write(100,*) " y_flap_spec = '"//trim(y_flap_spec)//"'"
  write(100,*) " TE_spec = '"//trim(TE_spec)//"'"
  write(100,*) " tcTE = ", tcTE
  write(100,*) " xltTE = ", xltTE
  write(100,*)
  do i = 1, noppoint
    write(text,*) i
    text = adjustl(text)
    write(100,*) " optimization_type("//trim(text)//") = '"//                    &
               trim(optimization_type(i))//"'"
    write(100,*) " optimization_correlation("//trim(text)//") = '"//             &
               trim(optimization_correlation(i))//"'"
    write(100,*) " op_mode("//trim(text)//") = '"//trim(op_mode(i))//"'"
    write(100,*) " op_point("//trim(text)//") = ", op_point(i)
    write(100,*) " use_previous_op("//trim(text)//") = ", use_previous_op(i)
    write(100,*) " op_point_start("//trim(text)//") = ", op_point_start(i)
    write(100,*) " op_point_end("//trim(text)//") = ", op_point_end(i)
    write(100,*) " op_point_step("//trim(text)//") = ", op_point_step(i)
    write(100,*) " target_value("//trim(text)//") = ", target_value(i)
    write(100,'(A,es17.8)') "  reynolds("//trim(text)//") = ", reynolds(i)
    write(100,*) " mach("//trim(text)//") = ", mach(i)
    write(100,*) " flap_selection("//trim(text)//") = '"//                       &
               trim(flap_selection(i))//"'"
    write(100,*) " flap_identical_op("//trim(text)//") = ", flap_identical_op(i)
    write(100,*) " flap_degrees("//trim(text)//") = ", flap_degrees(i)
    write(100,*) " weighting("//trim(text)//") = ", weighting(i)
    if (ncrit_pt(i) /= -1.d0)                                                  &
      write(100,*) " ncrit_pt("//trim(text)//") = ", ncrit_pt(i)
    if (i < noppoint) write(100,*)
  end do
  write(100,'(A)') "/"
  write(100,*)

! Constraints namelist

  write(100,'(A)') "&constraints"
  write(100,*) " min_thickness = ", min_thickness
  write(100,*) " max_thickness = ", max_thickness
  write(100,*) " max_growth_seed_mult = ", max_growth_seed_mult
  write(100,*) " min_leading_edge_angle = ", min_leading_edge_angle
  write(100,*) " max_leading_edge_angle = ", max_leading_edge_angle
  write(100,*) " dif_leading_edge_angle = ", dif_leading_edge_angle
  write(100,*) " min_te_angle = ", min_te_angle
  write(100,*) " te_angle_x_apply = ", te_angle_x_apply
  write(100,*) " max_panel_angle = ", max_panel_angle
  write(100,*) " lift_check_tol = ", lift_check_tol
  write(100,*) " drag_check_tol = ", drag_check_tol
  write(100,*) " check_curvature = ", check_curvature
  write(100,*) " max_curv_reverse_top = ", max_curv_reverse_top
  write(100,*) " max_curv_reverse_bot = ", max_curv_reverse_bot
  write(100,*) " curv_threshold = ", curv_threshold
  write(100,*) " symmetrical = ", symmetrical
  write(100,*) " min_flap_degrees = ", min_flap_degrees
  write(100,*) " max_flap_degrees = ", max_flap_degrees
  write(100,*) " min_flap_x = ", min_flap_x
  write(100,*) " max_flap_x = ", max_flap_x
  write(100,*) " min_tcTE = ", min_tcTE
  write(100,*) " max_tcTE = ", max_tcTE
  write(100,*) " min_camber = ", min_camber
  write(100,*) " max_camber = ", max_camber
  write(100,*)
  do i = 1, noppoint
    write(text,*) i
    text = adjustl(text)
    write(100,*) " moment_constraint_type("//trim(text)//") = '"//                &
               trim(moment_constraint_type(i))//"'"
    write(100,*) " min_moment("//trim(text)//") = ", min_moment(i)
  write(100,*) " lift_constraint_type("//trim(text)//") = '"//                &
               trim(lift_constraint_type(i))//"'"
    write(100,*) " min_lift("//trim(text)//") = ", min_lift(i)
  write(100,*) " drag_constraint_type("//trim(text)//") = '"//                &
               trim(drag_constraint_type(i))//"'"
    write(100,*) " max_drag("//trim(text)//") = ", max_drag(i)
  end do
  write(100,*)
  write(100,*) " naddthickconst = ", naddthickconst
  do i = 1, naddthickconst
    write(text,*) i
    text = adjustl(text)
    write(100,*) " addthick_x("//trim(text)//") = ", addthick_x(i)
    write(100,*) " addthick_min("//trim(text)//") = ", addthick_min(i)
    write(100,*) " addthick_max("//trim(text)//") = ", addthick_max(i)
  end do
  write(100,'(A)') "/"
  write(100,*)
  
  ! Aircraft_data namelist
  
  write(100,'(A)') "&aircraft_data"
  write(100,*) " weight = ", weight_i
  write(100,*) " S_w = ", S_w
  write(100,*) " A_w = ", A_w
  write(100,*) " e_w = ", e_w
  write(100,*) " S_expose = ", S_expose
  write(100,*) " thrust_coeff(1) = ", thrust_coeff(1)
  write(100,*) " thrust_coeff(2) = ", thrust_coeff(2)
  write(100,*) " thrust_coeff(3) = ", thrust_coeff(3)
  write(100,*) " fuselage%height = ", fuselage%height
  write(100,*) " fuselage%width = ", fuselage%width
  write(100,*) " fuselage%length = ", fuselage%length
  write(100,*) " fuselage%skin_roughness = ", fuselage%skin_roughness
  write(100,*) " fuselage%wetted_area = ", fuselage%wetted_area
  write(100,*) " fuselage%interderence_factor = ", fuselage%interference_factor
  write(100,*) " Cd_ld = ", Cd_ld
  write(100,*) " tail%tail_config = ", tail%config
  write(100,*) " add_drag = ", add_drag
  write(100,*)
  do i = 1, tail%config
    write(text,*) i
    text = adjustl(text)
    write(100,*) " tail%chord("//trim(text)//") = ", tail%chord(i)
    write(100,*) " tail%t_c_ratio("//trim(text)//") = ", tail%t_c_ratio(i)
    write(100,*) " tail%max_t_x("//trim(text)//") = ", tail%max_t_x(i)
    write(100,*) " tail%skin_roughness("//trim(text)//") = ", tail%skin_roughness(i)
    write(100,*) " tail%surface_area("//trim(text)//") = ", tail%surface_area(i)
    if (i < tail%config) write(100,*)
  end do
  write(100,'(A)') "/"
  write(100,*)
  
! Take_off_data namelist
    
  write(100,'(A)') "&take_off_data"
  write(100,*) " h = ", take_off%h
  write(100,*) " A_1 = ", take_off%A_1
  write(100,*) " miu = ", take_off%miu
  write(100,*) " S_g = ", take_off%S_g
  write(100,*) " weight_empty = ", take_off%weight_empty
  write(100,*) " weight_payload_ref = ", take_off%weight_payload_ref
  write(100,'(A)') "/"
  write(100,*)
  
! Climb_data namelist
   
  write(100,'(A)') "&climb_data"
  write(100,*) " acel_to_climb = ", climb%acel
  write(100,*) " time = ", climb%time
  write(100,*) " h = ", climb%h
  write(100,*) " dh = ", climb%dh
  write(100,*) " V_0 = ", climb%V_0
  do i=1, 5
    write(text,*) i
    text = adjustl(text)
    write(100,*) " points_coeff("//trim(text)//") = ", climb%points_coeff(i)
  end do
  write(100,'(A)') "/"
  write(100,*)
  
! Dash_data namelist
   
  write(100,'(A)') "&dash_data"
  write(100,*) " acel = ", dash%acel
  write(100,*) " time = ", dash%time
  write(100,*) " h = ", dash%h
  write(100,*) " dist_ref = ", dash%dist_ref
  write(100,*) " V_0 = ", dash%V_0
  write(100,*) " t_ex = ", dash%t_ex
  write(100,'(A)') "/"
  write(100,*)
  
! Turn_data namelist
   
  write(100,'(A)') "&Turn_data"
  write(100,*) " activation = ", turn%activation
  write(100,*) " h = ", turn%h
  write(100,*) " n = ", turn%n
  write(100,*) " dash_leg = ", turn%dash_leg
  write(100,'(A)') "/"
  write(100,*)

! naca_airfoil namelist

  write(100,'(A)') "&naca_airfoil"
  write(100,*) " family = "//trim(adjustl(family))
  write(100,*) " maxt = ", maxt
  write(100,*) " xmaxt = ", xmaxt
  write(100,*) " maxc = ", maxc
  write(100,*) " xmaxc = ", xmaxc
  write(100,*) " design_cl = ", design_cl
  write(100,*) " a = ", a
  write(100,*) " leidx = ", leidx
  write(100,*) " reflexed = ", reflexed
  write(100,'(A)') "/"
  write(100,*)

! Initialization namelist

  write(100,'(A)') "&initialization"
  write(100,*) " feasible_init = ", feasible_init
  write(100,*) " feasible_limit = ", feasible_limit
  write(100,*) " feasible_init_attempts = ", feasible_init_attempts
  write(100,'(A)') "/"
  write(100,*)

! Optimizer namelists
! Attention not all namelists are read but all will be writen, must use the 
! inputs

! Particle swarm namelist

  write(100,'(A)') "&particle_swarm_options"
  write(100,*) " pso_pop = ", pso_pop
  write(100,*) " pso_tol = ", pso_tol
  write(100,*) " pso_maxit = ", pso_maxit
  write(100,*) " pso_speed_limit = ", pso_speed_limit
  write(100,*) " pso_convergence_profile = '"//                                &
    trim(pso_convergence_profile)//"'"
  write(100,'(A)') "/"
  write(100,*)


! Genetic algorithm options

  write(100,'(A)') "&genetic_algorithm_options"
  write(100,*) " ga_pop = ", ga_pop
  write(100,*) " ga_tol = ", ga_tol
  write(100,*) " ga_maxit = ", ga_maxit
  write(100,*) " parents_selection_method = '"//                               &
    trim(parents_selection_method)//"'"
  write(100,*) " parent_fraction = ", parent_fraction 
  write(100,*) " roulette_selection_pressure = ",                              &
    roulette_selection_pressure
  write(100,*) " tournament_fraction = " , tournament_fraction
  write(100,*) " crossover_range_factor = ", crossover_range_factor
  write(100,*) " mutant_probability = ", mutant_probability
  write(100,*) " chromosome_mutation_rate = ",                                 &
    chromosome_mutation_rate
  write(100,*) " mutation_range_factor = ", mutation_range_factor
  write(100,'(A)') "/"
  write(100,*)

! Simplex search namelist

  write(100,'(A)') "&simplex_options"
  write(100,*) " simplex_tol = ", simplex_tol
  write(100,*) " simplex_maxit = ", simplex_maxit
  write(100,'(A)') "/"
  write(100,*)


! Xfoil run options namelist

  write(100,'(A)') "&xfoil_run_options"
  write(100,*) " ncrit = ", xfoil_options%ncrit
  write(100,*) " xtript = ", xfoil_options%xtript
  write(100,*) " xtripb = ", xfoil_options%xtripb
  write(100,*) " viscous_mode = ", xfoil_options%viscous_mode
  write(100,*) " silent_mode = ", xfoil_options%silent_mode
  write(100,*) " bl_maxit = ", xfoil_options%maxit
  write(100,*) " vaccel = ", xfoil_options%vaccel
  write(100,*) " reinitialize = ", xfoil_options%reinitialize
  write(100,*) " init_type = '"//trim(xfoil_options%init_type)//"'"
  write(100,*) " init_number_points = ", xfoil_options%init_number_points
  write(100,*) " init_al0 = ", xfoil_options%init_al0
  write(100,*) " init_cl0 = ", xfoil_options%init_cl0
  write(100,*) " init_initial_position = ", xfoil_options%init_initial_position
  write(100,*) " init_dist = '"//trim(xfoil_options%init_dist)//"'"
  write(100,'(A)') "/"
  write(100,*)

! Xfoil paneling options namelist

  write(100,'(A)') "&xfoil_paneling_options"
  write(100,*) " npan = ", xfoil_geom_options%npan
  write(100,*) " cvpar = ", xfoil_geom_options%cvpar
  write(100,*) " cterat = ", xfoil_geom_options%cterat
  write(100,*) " ctrrat = ", xfoil_geom_options%ctrrat
  write(100,*) " xsref1 = ", xfoil_geom_options%xsref1
  write(100,*) " xsref2 = ", xfoil_geom_options%xsref2
  write(100,*) " xpref1 = ", xfoil_geom_options%xpref1
  write(100,*) " xpref2 = ", xfoil_geom_options%xpref2
  write(100,'(A)') "/"
  write(100,*)

! Matchfoil options

  write(100,'(A)') "&matchfoil_options"
  write(100,*) " match_foils = ", match_foils
  write(100,*) " matchfoil_file = '"//trim(matchfoil_file)//"'"
  write(100,'(A)') "/"
  write(100,*)
  
  close(100)
  
! Check that inputs are reasonable
  write(*,'(A)') "| Checking if inputs are reasonable"
  write(*,'(A)') "|"
  
! Optimization settings

  if (trim(seed_airfoil) /= 'from_file' .and.                                  &
      trim(seed_airfoil) /= 'naca' .and.                                       &
      trim(seed_airfoil) /= 'from_variables')                                  &
    call my_stop("seed_airfoil must be 'from_file' or 'naca' or "//            &
      &"'from_variables'.")
  if (trim(seed_airfoil) .EQ. 'from_variables' .AND.                           &
     (trim(shape_functions) /= 'kulfan-bussoletti' .AND.                       &
      trim(shape_functions) /= 'bezier-parsec'))                               &
    call my_stop("seed_airfoil must not be 'from_variables' when using "//     &
      "'naca' or 'hicks-henne' or 'b-spline' shape_functions")
  if (trim(shape_functions) /= 'hicks-henne' .and.                             &
      trim(shape_functions) /= 'naca' .and.                                    &
      trim(shape_functions) /= 'kulfan-bussoletti' .and.                       &
      trim(shape_functions) /= 'b-spline' .and.                                &
      trim(shape_functions) /= 'bezier-parsec')                                &
    call my_stop("shape_functions must be 'hicks-henne' or 'naca' or "//       &
      "'kulfan-bussoletti' or 'b-spline' or 'bezier-parsec'.")
  if (nparameters_top < 0)                                                     &
    call my_stop("nparameters_top must be >= 0.")
  if (nparameters_bot < 0)                                                     &
    call my_stop("nparameters_bot must be >= 0.")
  if (abs_initial_perturb <= 0.d0)                                             &
    call my_stop("abs_initial_perturb must be > 0.")
  if (rel_initial_perturb < 0.d0)                                              &
    call my_stop("rel_initial_perturb must be >= 0.")
  if (penalty_limit_initial < 0.d0)                                            &
    call my_stop("penalty_limit_initial must be >= 0.")
  if (penalty_limit_end < 0.d0)                                                &
    call my_stop("penalty_limit_end must be >= 0.")
  if (penalty_factor < 0.d0)                                                   &
    call my_stop("penalty_factor must be >= 0.")
  if (min_bump_width <= 0.d0)                                                  &
    call my_stop("min_bump_width must be > 0.")
  if (b_spline_degree < 2)                                                     &
    call my_stop("b_spline_degree must be >= 2")
  if (b_spline_xtype /= 1 .and. b_spline_xtype /= 2)                           &
    call my_stop("b_spline_xtype must be 1 or 2")
  if (b_spline_distribution /= 1 .and. b_spline_distribution /= 2 .and.        &
      b_spline_distribution /= 3)                                              &
    call my_stop("b_spline_distribution must be 1 or 2 or 3")
  if (restart_write_freq < 2)                                                  &
    call my_stop("restart_write_freq must be >= 2.")
  if ( (flap_optimization_only) .and. ((.not. use_flap) .or. (match_foils)) )  &
    call my_stop("flap_optimization_only can only be used with use_flap=.true. &
      &and match_foils=.false.")
  if ( (flap_optimization_only) .and. (nflap_optimize==0))                     &
    call my_stop("flap_optimization_only can only be used with with a number   &
      &of flaps to optimize diferent from 0")
  if (trim(progress_per_eval) /= 'full' .and.                                  &
      trim(progress_per_eval) /= 'none' .and.                                  &
      trim(progress_per_eval) /= 'part')                                    &
    call my_stop("progress_per_eval must be 'full' or 'none' or 'part'.")    
  
! Operating points

  if (noppoint < 1) call my_stop("noppoint must be > 0.")
  if (noppoint > max_op_points) then
     write(text,*) max_op_points
     text = adjustl(text)
     call my_stop("noppoints must be <= "//trim(text)//".")
  end if
  if (use_flap .and.                                                           &
      trim(flap_connection) /= 'sharp' .and.                                   &
      trim(flap_connection) /= 'smooth-top' .and.                              &
      trim(flap_connection) /= 'smooth-bot' .and.                              &
      trim(flap_connection) /= 'smooth')                                       &
    call my_stop("flap_connection must be 'sharp', 'smooth-top', "//           &
    &            "'smooth-bot' or 'smooth'.")
  if (use_flap .and.                                                           &
      trim(connection_apply) /= 'both' .and.                                   &
      trim(connection_apply) /= 'geometric' .and.                              &
      trim(connection_apply) /= 'trip_wire' .and.                              &
      trim(connection_apply) /= 'none')                                        &
    call my_stop("connection_apply must be 'none', 'geometric', 'trip_wire'"// &
    &            " or 'both'")
  if (use_flap .and.                                                           &
      trim(connection_apply) /= 'trip_wire' .and.                              &
      trim(connection_apply) /= 'none'  .and.                                  &
      connection_radius .LE. 0.0d0)                                            &
    call my_stop("connection_radius must be > 0.")
  if ((nxtr_opt .NE. 0) .AND. (trim(connection_apply) /= 'none' .OR.           &
                               trim(connection_apply) /= 'geometric'))         &
    call my_stop("connection_apply must be 'none' or 'geometric' if using "//  &
    &            "max-xtr optimization type.")
  if ((use_flap) .and. (x_flap <= 0.0)) call my_stop("x_flap must be > 0.")
  if ((use_flap) .and. (x_flap >= 1.0)) call my_stop("x_flap must be < 1.")
  if ((use_flap) .and. (x_flap_spec /= 'specify')                              &
    .and. (x_flap_spec /= 'optimize'))                                         &
    call my_stop("x_flap_spec must be 'specify' or 'optimize'.")
  if ((trim(x_flap_spec) /= 'specify') .and. (x_flap < min_flap_x))            &
    call my_stop("x_flap must be >= min_flap_x.")
  if ((trim(x_flap_spec) /= 'specify') .and. (x_flap > max_flap_x))            &
        call my_stop("x_flap must be <= max_flap_x.")
  if ((use_flap) .and. (y_flap_spec /= 'y/c') .and. (y_flap_spec /= 'y/t'))    &
    call my_stop("y_flap_spec must be 'y/c' or 'y/t'.")
  if (trim(TE_spec) /= 'use_seed' .and.                    &
      trim(TE_spec) /= 'specify' .and.                       &
      trim(TE_spec) /= 'optimize')                               &
      call my_stop("TE_spec must be 'use_seed', 'specify', or 'optimize'.")
  if (tcTE < 0.0) call my_stop("tcTE must be > 0.")
  if (xltTE < 0.0) call my_stop("xltTE must be > 0.")
  if (xltTE >= 1.0) call my_stop("xltTE must be < 1")
  
  if (use_previous_op(1)) call my_stop("op 1 can not have a previous op.")
  do i = 1, noppoint
    if (trim(op_mode(i)) /= 'spec-cl' .and. trim(op_mode(i)) /= 'spec-al')     &
      call my_stop("op_mode must be 'spec-al' or 'spec-cl'.")
    if (reynolds(i) <= 0.d0) call my_stop("reynolds must be > 0.")
    if (mach(i) < 0.d0) call my_stop("mach must be >= 0.")
    if (trim(flap_selection(i)) /= 'specify' .and.                             &
        trim(flap_selection(i)) /= 'optimize' .and.                            &
        trim(flap_selection(i)) /= 'identical')                                &
      call my_stop("flap_selection must be 'specify' or 'optimize' or "//      &
                   "'identical'.")
    
    if (flap_degrees(i) < -90.d0) call my_stop("flap_degrees must be > -90.")
    if (flap_degrees(i) > 90.d0) call my_stop("flap_degrees must be < 90.")
    if (weighting(i) <= 0.d0) call my_stop("weighting must be > 0.")
    if (trim(optimization_type(i)) /= 'min-drag' .and.                         &
      trim(optimization_type(i)) /= 'max-glide' .and.                          &
      trim(optimization_type(i)) /= 'min-sink' .and.                           &
      trim(optimization_type(i)) /= 'max-lift' .and.                           &
      trim(optimization_type(i)) /= 'max-xtr' .and.                            &
      trim(optimization_type(i)) /= 'target-lift' .and.                        &
      trim(optimization_type(i)) /= 'target-drag' .and.                        &
      trim(optimization_type(i)) /= 'target-moment' .and.                      &
      trim(optimization_type(i)) /= 'target-xtrt' .and.                        &
      trim(optimization_type(i)) /= 'target-xtrb' .and.                        &
      trim(optimization_type(i)) /= 'target-glide' .and.                       &
      trim(optimization_type(i)) /= 'target-sink'.and.                         &
      trim(optimization_type(i)) /= 'max-lift-slope'.and.                      &
      trim(optimization_type(i)) /= 'max-lift-search'.and.                     &
      trim(optimization_type(i)) /= 'take-off'.and.                            &
      trim(optimization_type(i)) /= 'climb'.and.                               &
      trim(optimization_type(i)) /= 'dash'.and.                                &
      trim(optimization_type(i)) /= 'turn')                                    &
      call my_stop("optimization_type must be 'min-drag', 'max-glide', "//     &
                   "min-sink', 'max-lift', 'max-xtr', "//                      &
                   "'target-lift', 'target-drag', 'target-moment', "//         &
                   "'target-xtrt', 'target-xtrb', 'target-glide', "//          &
                   "'target-sink', 'max-lift-slope', 'max-lift-search', "//    &
                   "'take-off', 'climb' or 'dash'.")
    if (trim(optimization_type(i)) == 'max-lift-slope' .and. (noppoint == 1))  &
      call my_stop("at least two operating points are required for to "//      &
                   "maximize lift curve slope.")
    if (((trim(optimization_type(i)) == 'take-off') .or.                       &
        (trim(optimization_type(i)) == 'climb') .or.                           &
        (trim(optimization_type(i)) == 'dash') .or.                            &
        (trim(optimization_type(i)) == 'turn')) .and. (noppoint == 1))         &
      call my_stop("at least two operating points are required for to "//      &
                   "to evaluate performance ('take-off', 'climb', 'dash', "//  &
                   "'turn').")
    if (((trim(optimization_type(i)) == 'take-off') .or.                       &
        (trim(optimization_type(i)) == 'climb') .or.                           &
        (trim(optimization_type(i)) == 'dash') .or.                            &
        (trim(optimization_type(i)) == 'turn')) .and.                          &
        (trim(optimization_correlation(i)) == 'none'))                         &
      call my_stop("optimization_correlationmust be different from 'none' "//  &
                   "to evaluate performance ('take-off', 'climb', 'dash', "//  &
                   "'turn').")
    if ((trim(optimization_correlation(i)) /= 'init') .and.                    &
        (trim(optimization_correlation(i)) /= 'mid') .and.                     &
        (trim(optimization_correlation(i)) /= 'end') .and.                     &
        (trim(optimization_correlation(i)) /= 'none'))                         &
      call my_stop("optimization_correlation must be 'none', 'init', "//       &
                   "'mid' or 'end'.")
    if (ncrit_pt(i) <= 0.d0) call my_stop("ncrit_pt must be > 0 or -1.")
  end do

  ! Check if flap_identical_op does not refer to a point in flap_identical_points

  do i = 1, nflap_identical
    do j = 1, nflap_identical
      if (flap_identical_op(flap_identical_points(j)) .EQ.                     &
        flap_identical_points(i)) then
        write(text,*) flap_identical_points(j)
        text = adjustl(text)
        call my_stop("identical flap angle must refer to a 'specify' or "//    &
                     "'optimize' flap type at oppoint("//trim(text)//")")
      end if
    end do
  end do
  
! Constraints

  if (min_thickness <= 0.d0) call my_stop("min_thickness must be > 0.")
  if (max_thickness <= 0.d0) call my_stop("max_thickness must be > 0.")
  if (min_thickness >= max_thickness)                                          &
    call my_stop("min_thickness must be < max_thickness.")
  do i = 1, noppoint
    if (trim(moment_constraint_type(i)) /= 'use_seed' .and.                    &
      trim(moment_constraint_type(i)) /= 'specify' .and.                       &
      trim(moment_constraint_type(i)) /= 'none')                               &
      call my_stop("moment_constraint_type must be 'use_seed', 'specify', "//  &
                 "or 'none'.")
    if (trim(lift_constraint_type(i)) /= 'use_seed' .and.                    &
      trim(lift_constraint_type(i)) /= 'specify' .and.                       &
      trim(lift_constraint_type(i)) /= 'none')                               &
      call my_stop("lift_constraint_type must be 'use_seed', 'specify', "//  &
                 "or 'none'.")
    if (trim(drag_constraint_type(i)) /= 'use_seed' .and.                    &
      trim(drag_constraint_type(i)) /= 'specify' .and.                       &
      trim(drag_constraint_type(i)) /= 'none')                               &
      call my_stop("drag_constraint_type must be 'use_seed', 'specify', "//  &
                 "or 'none'.")
  end do
  
  if (max_growth_seed_mult < 1.d0) call my_stop("max_growth_seed_mult"//       &
    &" should be >= 1.", "warm")
  if (min_leading_edge_angle < 0.d0) call my_stop("min_leading_edge_angle"//   &
    &" must be >= 0.")
  if (max_leading_edge_angle > 90.00d0) call my_stop("max_leading_edge_angle"//&
    &" must be <= 90.00")
  if (dif_leading_edge_angle < 10.d0) call my_stop("dif_leading_edge_angle"//  &
    &" should be >= 10, 20 is recommended.", "warn")
  if (min_te_angle < 0.d0) call my_stop("min_te_angle must be >= 0.")
  if (te_angle_x_apply < 0.5d0) call my_stop("te_angle_x_apply"//             &
    &" should be after max thickness", "warn")
  if (max_panel_angle /= 25.d0) call my_stop(" recommended value for"//        &
    &" max_panel_angle is 25. degrees", "warn")
  if (check_curvature .and. (curv_threshold <= 0.d0))                          &
    call my_stop("curv_threshold must be > 0.")
  if (check_curvature .and. (max_curv_reverse_top < 0))                        &
    call my_stop("max_curv_reverse_top must be >= 0.")
  if (check_curvature .and. (max_curv_reverse_bot < 0))                        &
    call my_stop("max_curv_reverse_bot must be >= 0.")
  if (symmetrical)                                                             &
    write(*,*) "Mirroring top half of seed airfoil for symmetrical constraint."
  if (min_flap_degrees >= max_flap_degrees)                                    &
    call my_stop("min_flap_degrees must be < max_flap_degrees.")
  if (min_flap_degrees <= -90.d0)                                              &
    call my_stop("min_flap_degrees must be > -90.")
  if (max_flap_degrees >= 90.d0)                                               &
    call my_stop("max_flap_degrees must be < 90.")
  if (min_flap_x >= max_flap_x)                                    &
    call my_stop("min_flap_x must be < max_flap_x.")
  if (min_flap_x <= 0.d0)                                              &
    call my_stop("min_flap_x must be > 0.")
  if (max_flap_x >= 1.d0)                                               &
    call my_stop("max_flap_x must be < 1.")
  if (min_tcTE >= max_tcTE)                                    &
    call my_stop("min_tcTE must be < max_tcTE.")
  if (min_tcTE < 0.d0)                                              &
    call my_stop("min_tcTE must be >= 0.")
  if (max_tcTE >= 1.d0)                                               &
    call my_stop("max_tcTE must be < 1.")
  if (min_camber >= max_camber)                                                &
    call my_stop("min_camber must be < max_camber.")
  
  if (naddthickconst > max_addthickconst) then
     write(text,*) max_addthickconst
     text = adjustl(text)
     call my_stop("naddthickconst must be <= "//trim(text)//".")
  end if
  do i = 1, naddthickconst
    if (addthick_x(i) <= 0.d0) call my_stop("addthick_x must be > 0.")
    if (addthick_x(i) >= 1.d0) call my_stop("addthick_x must be < 1.")
    if (addthick_min(i) <= 0.d0) call my_stop("addthick_min must be > 0.")
    if (addthick_max(i) <= 0.d0) call my_stop("addthick_max must be > 0.")
    if (addthick_min(i) >= addthick_max(i))                                    &
      call my_stop("addthick_min must be < addthick_max.")
  end do

! Naca airfoil options

  select case (adjustl(family))
    case ('4', '4M', '5', '63', '64', '65', '66', '67', '63A', '64A', '65A')
      continue
    case default
      call my_stop("Unrecognized NACA airfoil family.")
  end select
  if (maxt <= 0.d0) call my_stop("maxt must be > 0.")
  if ( (xmaxt < 0.d0) .or. (xmaxt > 1.d0) )                                    &
    call my_stop("xmaxt must be >= 0 and <= 1.")
  if ( (xmaxc < 0.d0) .or. (xmaxc > 1.d0) )                                    &
    call my_stop("xmaxc must be >= 0 and <= 1.")
  if ( (a < 0.d0) .or. (a > 1.d0) )                                            &
    call my_stop("a must be >= 0 and <= 1.")
  if (leidx <= 0.d0) call my_stop("leidx must be > 0.")

! Initialization options
    
  if ((feasible_limit <= 0.d0) .and. feasible_init)                            &
    call my_stop("feasible_limit must be > 0.")
  if ((feasible_init_attempts < 1) .and. feasible_init)                        &
    call my_stop("feasible_init_attempts must be > 0.")

! Optimizer options

  if (trim(search_type) == 'global' .or.                                       &
       trim(search_type) == 'global_and_local') then

    if (trim(global_search) == 'particle_swarm') then

!     Particle swarm options

      if (pso_pop < 1) call my_stop("pso_pop must be > 0.")
      if (pso_tol <= 0.d0) call my_stop("pso_tol must be > 0.")
      if (pso_maxit < 1) call my_stop("pso_maxit must be > 0.")
      if (pso_speed_limit <= 0.d0)                                             &
        call my_stop("pso_speed_limit must be > 0.")
      if ( (trim(pso_convergence_profile) /= "quick") .and.                    &
           (trim(pso_convergence_profile) /= "exhaustive") .and.               &
           (trim(pso_convergence_profile) /= "standard") )                     &
        call my_stop("pso_convergence_profile must be 'exhaustive' "//&
                     "or 'quick' or 'standard'.")

    else if (trim(global_search) == 'genetic_algorithm') then

!     Genetic algorithm options

      if (ga_pop < 1) call my_stop("ga_pop must be > 0.")
      if (ga_tol <= 0.d0) call my_stop("ga_tol must be > 0.")
      if (ga_maxit < 1) call my_stop("ga_maxit must be > 0.")
      if ( (trim(parents_selection_method) /= "roulette") .and.                &
           (trim(parents_selection_method) /= "tournament") .and.              &
           (trim(parents_selection_method) /= "random") )                      &
        call my_stop("parents_selection_method must be 'roulette', "//&
                     "'tournament', or 'random'.")
      if ( (parent_fraction <= 0.d0) .or. (parent_fraction > 1.d0) )           &
        call my_stop("parent_fraction must be > 0 and <= 1.")
      if (roulette_selection_pressure <= 0.d0)                                 &
        call my_stop("roulette_selection_pressure must be > 0.")
      if ( (tournament_fraction <= 0.d0) .or. (tournament_fraction > 1.d0) )   &
        call my_stop("tournament_fraction must be > 0 and <= 1.")
      if (crossover_range_factor < 0.d0)                                       &
        call my_stop("crossover_range_factor must be >= 0.")
      if ( (mutant_probability < 0.d0) .or. (mutant_probability > 1.d0) )      &
        call my_stop("mutant_probability must be >= 0 and <= 1.") 
      if (chromosome_mutation_rate < 0.d0)                                     &
        call my_stop("chromosome_mutation_rate must be >= 0.")
      if (mutation_range_factor < 0.d0)                                        &
        call my_stop("mutation_range_factor must be >= 0.")

    end if

  end if

  if (trim(search_type) == 'local' .or.                                        &
       trim(search_type) == 'global_and_local') then

!   Simplex options

    if (simplex_tol <= 0.d0) call my_stop("simplex_tol must be > 0.")
    if (simplex_maxit < 1) call my_stop("simplex_maxit must be > 0.")  
  
  end if

! XFoil run options

  if (ncrit < 0.d0) call my_stop("ncrit must be >= 0.")
  if (xtript < 0.d0 .or. xtript > 1.d0)                                        &
    call my_stop("xtript must be >= 0. and <= 1.")
  if (xtripb < 0.d0 .or. xtripb > 1.d0)                                        &
    call my_stop("xtripb must be >= 0. and <= 1.")
  if (bl_maxit < 1) call my_stop("bl_maxit must be > 0.")
  if (vaccel < 0.d0) call my_stop("vaccel must be >= 0.")
  if (trim(init_type) /= 'always' .and. trim(init_type) /= 'never' .and.       &
    trim(init_type) /= 'unconverged')                                          &
    call my_stop("init_type must be 'always', 'unconverged' or 'never'")
  if (init_number_points < 1) call my_stop("init_number_points must be > 0.")
  if (init_initial_position < 0.d0 .or. init_initial_position > 1.d0)                                        &
    call my_stop("init_initial_position must be >= 0. and <= 1.")
  if (trim(init_dist) /= 'linear' .and. trim(init_dist) /= 'sine')             &
      call my_stop("init_dist must be 'linear' or 'sine'")
  
! XFoil paneling options

  if (npan < 20) call my_stop("npan must be >= 20.")
  if (cvpar <= 0.d0) call my_stop("cvpar must be > 0.")
  if (cterat <= 0.d0) call my_stop("cterat must be > 0.")
  if (ctrrat <= 0.d0) call my_stop("ctrrat must be > 0.")
  if (xsref1 < 0.d0) call my_stop("xsref1 must be >= 0.")
  if (xsref2 < xsref1) call my_stop("xsref2 must be >= xsref1")
  if (xsref2 > 1.d0) call my_stop("xsref2 must be <= 1.")
  if (xpref1 < 0.d0) call my_stop("xpref1 must be >= 0.")
  if (xpref2 < xpref1) call my_stop("xpref2 must be >= xpref1")
  if (xpref2 > 1.d0) call my_stop("xpref2 must be <= 1.")

  write(*,'(A)') "|"
  write(*,'(A)') "| End of Checking "
end subroutine read_inputs

!=============================================================================80
!
! Prints error and stops or warns for bad namelist read
!
!=============================================================================80
subroutine namelist_check(nmlname, errcode, action_missing_nml)

  character(*), intent(in) :: nmlname
  integer, intent(in) :: errcode
  character(*), intent(in) :: action_missing_nml

  if (errcode < 0) then
    write(*,*)
    if (trim(action_missing_nml) == 'warn') then
      write(*,'(A)') 'Warning: namelist '//trim(nmlname)//&
                     ' not found in input file.'
      write(*,'(A)') 'Using default values.'
      write(*,*)
    else
      write(*,'(A)') 'Warning: namelist '//trim(nmlname)//&
                     ' is required and was not found in input file.'
      write(*,*)
      stop
    end if
  else if (errcode > 0) then
    write(*,*)
    write(*,'(A)') 'Error: unrecognized variable in namelist '//trim(nmlname)//&
                   '.'
    write(*,'(A)') 'See User Guide for correct variable names.'
    write(*,*)
    stop
  else
    continue
  end if

end subroutine namelist_check

!=============================================================================80
!
! Reads command line arguments for input file name and output file prefix
!
!=============================================================================80
subroutine read_clo(input_file, output_prefix, exename)

  use airfoil_operations, only : my_stop

  character(*), intent(inout) :: input_file, output_prefix
  character(*), intent(in), optional :: exename

  character(80) :: arg, exeprint
  integer i, nargs
  logical getting_args

  if (present(exename)) then
    exeprint = exename
  else
    exeprint = "xoptfoil"
  end if 

  nargs = iargc()
  if (nargs > 0) then
    getting_args = .true.
  else
    getting_args = .false.
  end if

  i = 1
  do while (getting_args)
    call getarg(i, arg) 

    if (trim(arg) == "-i") then
      if (i == nargs) then
        call my_stop("Must specify an input file with -i option.")
      else
        call getarg(i+1, input_file)
        i = i+2
      end if
    else if (trim(arg) == "-o") then
      if (i == nargs) then
        call my_stop("Must specify an output prefix with -o option.")
      else
        call getarg(i+1, output_prefix)
        i = i+2
      end if
    else if ( (trim(arg) == "-v") .or. (trim(arg) == "--version") ) then
      call print_version()
      stop
    else if ( (trim(arg) == "-h") .or. (trim(arg) == "--help") ) then
      call print_usage(exeprint)
      stop
    else
      write(*,'(A)') "Unrecognized option "//trim(arg)//"."
      call print_usage(exeprint)
      stop 1
    end if

    if (i > nargs) getting_args = .false.
  end do

end subroutine read_clo

!=============================================================================80
!
! Prints version information
!
!=============================================================================80
subroutine print_version()

  write(*,'(A)') "Xoptfoil "//trim(PACKAGE_VERSION)
  write(*,'(A)') "Copyright (C) 2017-2019 Daniel Prosser, "//                  &
                 "2020-2021 Ricardo Palmeira"
  write(*,'(A)') "License GPLv3+: GNU GPL version 3 or later "//               &
                 "<http://gnu.org/licenses/gpl.html>"
  write(*,'(A)') "This is free software: you are free to change and "//        &
                 "redistribute it."
  write(*,'(A)') "There is NO WARRANTY, to the extent permitted by law."

end subroutine print_version

!=============================================================================80
!
! Prints usage information
!
!=============================================================================80
subroutine print_usage(exeprint)

  character(*), intent(in) :: exeprint

  write(*,'(A)') "Usage: "//trim(exeprint)//" [OPTION]"
  write(*,'(A)')
  write(*,'(A)') "Options:"
  write(*,'(A)') "  -i input_file     Specify a non-default input file"
  write(*,'(A)') "  -o output_prefix  Specify a non-default output prefix"
  write(*,'(A)') "  -h, --help        Display usage information and exit"
  write(*,'(A)') "  -v, --version     Display Xoptfoil version and exit"
  write(*,'(A)')
  write(*,'(A)') "Refer to the PDF user guide for complete input help."
  write(*,'(A)')
  write(*,'(A)') "Home page: https://github.com/r-palmeira-pt/Xoptfoil"
  write(*,'(A)') "Development page: https://github.com/r-palmeira-pt/Xoptfoil"
  write(*,'(A)') "Report bugs using the issue reporting system at either "//   &
                 "of the above sites."

end subroutine print_usage

!=============================================================================80
!
! Asks user to turn off pitching moment constraints
!
!=============================================================================80
function ask_moment_constraints()

  character :: ask_moment_constraints
  logical :: valid_choice

! Get user input

  valid_choice = .false.
  do while (.not. valid_choice)
  
    write(*,*)
    write(*,'(A)') 'Warning: pitching moment constraints not recommended for '
    write(*,'(A)', advance='no') 'symmetrical airfoil optimization. '//&
                                 'Turn them off now? (y/n): '
    read(*,'(A)') ask_moment_constraints

    if ( (ask_moment_constraints == 'y') .or.                                  &
         (ask_moment_constraints == 'Y') ) then
      valid_choice = .true.
      ask_moment_constraints = 'y'
      write(*,*)
      write(*,*) "Setting moment_constraint_type(:) = 'none'."
    else if ( (ask_moment_constraints == 'n') .or.                             &
         (ask_moment_constraints == 'N') ) then
      valid_choice = .true.
      ask_moment_constraints = 'n'
    else
      write(*,'(A)') 'Please enter y or n.'
      valid_choice = .false.
    end if

  end do

end function ask_moment_constraints

!=============================================================================80
!
! Asks user to turn off forced transition
!
!=============================================================================80
function ask_forced_transition()

  character :: ask_forced_transition
  logical :: valid_choice

! Get user input

  valid_choice = .false.
  do while (.not. valid_choice)
  
    write(*,*)
    write(*,'(A)') 'Warning: using max-xtr optimization but xtript or xtripb'
    write(*,'(A)', advance='no') 'is less than 1. Set them to 1 now? (y/n): '
    read(*,'(A)') ask_forced_transition

    if ( (ask_forced_transition == 'y') .or.                                  &
         (ask_forced_transition == 'Y') ) then
      valid_choice = .true.
      ask_forced_transition = 'y'
      write(*,*)
      write(*,*) "Setting xtript and xtripb to 1."
    else if ( (ask_forced_transition == 'n') .or.                             &
         (ask_forced_transition == 'N') ) then
      valid_choice = .true.
      ask_forced_transition = 'n'
    else
      write(*,'(A)') 'Please enter y or n.'
      valid_choice = .false.
    end if

  end do

end function ask_forced_transition

end module input_output
