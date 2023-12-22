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

module vardef

! Data structures for airfoil optimization code

  implicit none

  type airfoil_type

    integer :: npoint
    double precision, dimension(:), allocatable :: x, z ! Airfoil coordinates
    double precision :: xle, zle                        ! Leading edge coords
    integer :: leclose                                  ! Index closest to LE
    integer :: addpoint_loc                             ! Whether to add point 
                                                        !  for LE before or
                                                        !  after leclose

  end type airfoil_type

  type op_search_type

    integer :: noppoint = 0
    integer, dimension(:), allocatable :: oppoints
    double precision, dimension(:), allocatable :: op_start
    double precision, dimension(:), allocatable :: op_end
    double precision, dimension(:), allocatable :: op_step
    
  end type op_search_type
  
  type objfunction_type

    double precision :: value
    integer :: message_code
    character(200) :: message
    double precision, dimension(:), allocatable :: constrains_data
    double precision, dimension(:), allocatable :: aero_data
    
  end type objfunction_type
  
  type test_type
    double precision :: value
    type(objfunction_type) :: obj
  end type test_type
  
  type fuselage_type

    double precision :: height, width
    double precision :: length, skin_roughness
    double precision :: wetted_area, interference_factor
    
  end type fuselage_type
    
  type tail_type

    integer :: config
    double precision, dimension(2) :: chord, t_c_ratio, max_t_x,               &     
			     skin_roughness, surface_area
  end type tail_type 
  
  type take_off_type
      
    double precision :: h, A_1, miu, S_g, weight_empty, weight_payload_ref 
    double precision :: CL_max, CL_run, CD_run, V_to, V_run, points
    
  end type take_off_type
  
  type climb_type
      
    logical :: acel
    double precision :: time, h, dh, V_0
    double precision, dimension(5) :: points_coeff
    double precision :: RC_max, V, Cl, D, T, t_acel, points
    
  end type climb_type
  
  type dash_type
    
    logical :: acel
    double precision :: time, h, dist_ref, V_0, t_ex
    double precision :: V_max, t_acel_d, dist_acel_d, dist, points
    
  end type dash_type
  
  type turn_type
    
    logical :: activation
    double precision :: h, n, dash_leg
    double precision :: V, radius
    
  end type turn_type
  
! Main return variable for all optimizations
  type(objfunction_type) :: objfunction_return
  !$omp threadprivate(objfunction_return)
  
! Global variables (mainly needed to preserve generality of optimization
! routines)
  
  integer :: noppoint
  integer, parameter :: max_op_points = 100
  double precision, dimension(:), allocatable :: xseedt, xseedb, zseedt, zseedb
  character(7), dimension(max_op_points) :: op_mode
  character(9), dimension(max_op_points) :: flap_selection
  double precision, dimension(max_op_points) :: op_point, reynolds, mach,      &
                                 flap_degrees, weighting, scale_factor,        &
                                 ncrit_pt, target_value
  type(op_search_type) :: op_search
  integer, dimension(max_op_points) :: flap_identical_op
  double precision :: x_flap, y_flap, tcTE, xltTE, connection_radius 
  character(3) :: y_flap_spec
  character(8) :: x_flap_spec, TE_spec ! added x_flap_spec and TE_spec
  integer :: int_x_flap_spec, int_tcTE_spec           ! added int_x_flap_spec
  logical :: use_flap, flap_optimization_only
  character(20) :: flap_connection, connection_apply 
  character(15), dimension(max_op_points) :: optimization_type, optimization_correlation
  logical, dimension(max_op_points) :: use_previous_op
  integer, dimension(6) :: dvs_for_type
  integer :: nflap_optimize, nflap_identical, nmoment_constrain,               &
             ndrag_constrain, nlift_constrain, contrain_number
                                     ! Number of operating points where flap 
                                     !   setting will be optimized or identical
  integer, dimension(max_op_points) :: flap_optimize_points,                   &
                                       flap_identical_points,                  &
           moment_constrain_points, drag_constrain_points, lift_constrain_points
  double precision :: xoffset, zoffset, foilscale, foilangle

  type(airfoil_type) :: curr_foil, match_foil
  ! added min_flap_x, max_flap_x, max_tcTE, min_tcTE
  double precision :: min_thickness, max_thickness, min_te_angle,              &
                      growth_allowed, min_flap_degrees, max_flap_degrees,      &
                      min_flap_x, max_flap_x, max_tcTE, min_tcTE, min_camber,  &
                      max_camber, max_growth_seed_mult, min_leading_edge_angle,&
                      max_leading_edge_angle, dif_leading_edge_angle,          &
                      te_angle_x_apply, lift_check_tol, drag_check_tol,        &
                      max_panel_angle
  double precision :: curv_threshold
  integer :: max_curv_reverse_top, max_curv_reverse_bot
  ! added lift and drag constraints
  character(8), dimension(max_op_points) :: moment_constraint_type,            &
                                          lift_constraint_type,                &
                                          drag_constraint_type
  character(4) :: progress_per_eval
  character(80) :: airfoil_file
  double precision, dimension(max_op_points) :: min_moment, min_lift, max_drag

  character(17) :: shape_functions
  double precision, dimension(:), allocatable :: xmatcht, xmatchb, zmatcht,    &
                                                 zmatchb
  logical :: match_foils
  logical :: check_curvature
  logical :: symmetrical
  logical :: write_dvs_file, write_cp_file, write_bl_file
  integer :: number_threads

  integer :: nparams_top, nparams_bot, nshapedvtop, nshapedvbot
  double precision, dimension(:), allocatable :: modest_seed
  double precision, dimension(:), allocatable :: modesb_seed
  
  double precision :: abs_initial_perturb, rel_initial_perturb, tcTE_seed
  double precision :: penalty_limit_initial, penalty_limit_end, penalty_factor
  double precision :: min_bump_width
  logical :: kulfan_bussoletti_LEM, epsexit_linear
  integer :: int_kulfan_bussoletti_LEM
  integer :: b_spline_degree, b_spline_xtype, b_spline_distribution
  double precision, dimension(:), allocatable :: upointst, upointsb, xcontrolt,&
                                                 xcontrolb
  logical ::  first_run_xfoil, allow_seed_penalties
  character(80) :: output_prefix

  integer :: naddthickconst
  integer :: maxit
  integer, parameter :: max_addthickconst = 10 ! remove limit?
  double precision, dimension(max_addthickconst) :: addthick_x, addthick_min,  &
                                                    addthick_max
  character(19) :: restart_stat
  character(80) :: global_search_stat, local_search_stat
  !add aircraft data
  type(fuselage_type) :: fuselage
  type(tail_type) :: tail
  type(take_off_type) :: take_off
  type(climb_type) :: climb
  type(dash_type) :: dash
  type(turn_type) :: turn
  
  double precision :: weight, weight_i
  double precision :: A_w, e_w, S_w, S_expose
  double precision, dimension(3) :: thrust_coeff
  double precision :: height, width
  double precision :: length, f_skin_roughness
  double precision :: wetted_area, interference_factor
  double precision :: tail_config
  double precision, dimension(2) :: tail_chord, t_c_ratio, max_t_x
  double precision, dimension(2) :: t_skin_roughness, tail_surface_area
  double precision :: Cd_ld, add_drag
  double precision :: h_take_off, A_1, miu, S_g
  double precision :: weight_empty, weight_payload_ref
  logical :: acel_to_climb, acel_to_dash
  double precision :: time_climb, h_climb, dh, V_0_climb
  double precision, dimension(5) :: points_coeff
  double precision :: time_dash, h_dash, V_0_dash, time_extra, dist_ref_dash
  logical :: activation_turn
  double precision :: h_turn, n_turn, dash_leg
    
!$omp threadprivate(curr_foil)

end module vardef
