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

module xfoil_driver

! Contains subroutines to use XFoil to analyze an airfoil

  implicit none

  type xfoil_options_type

    double precision :: ncrit          !Critical ampl. ratio
    double precision :: xtript, xtripb !Trip locations
    logical :: viscous_mode                       
    logical :: silent_mode             !Toggle xfoil screen write
    integer :: maxit                   !Iterations for BL calcs
    double precision :: vaccel         !Xfoil BL convergence accelerator
    logical :: reinitialize            !Reinitialize BLs at every operating
                                       !  point (recommended for optimization)
    character(20) :: init_type         !Interval initialization types:
                      !  'always' - initialize to get op point
                      !  'unconverged' - initialize if op point unconverged
                      !  'never' - do not initialize
    integer :: init_number_points      !Number of points to use in initialization
    double precision :: init_al0       !Alpha reference value
    double precision :: init_cl0       !Cl reference value
    double precision :: init_initial_position !Adimensional inicial interval position
                      !Interval is between reference value (0) and op point (1) 
    character(20) :: init_dist         !Distribution of init points. Either 'linear' or 'sine'

  end type xfoil_options_type

  type xfoil_geom_options_type

    integer :: npan
    double precision :: cvpar, cterat, ctrrat, xsref1, xsref2, xpref1, xpref2

  end type xfoil_geom_options_type

  type xfoil_file_options_type
    
    integer :: design_number = 0      !Design number
    logical :: polar = .false.        !Save polar file ?
    logical :: cp = .false.           !Save cp file ?
    logical :: bl = .false.           !Save bl file ?
    logical :: write_all_airfoils = .false. !Save all airfoils with flaps
  
  end type xfoil_file_options_type
  
  
  contains

!=============================================================================80
!
! Subroutine to smooth an airfoil using Xfoil's PANGEN subroutine
!
!=============================================================================80
subroutine smooth_paneling(foilin, geom_options, foilout)

  use xfoil_inc
  use vardef, only : airfoil_type

  type(airfoil_type), intent(in) :: foilin
  type(xfoil_geom_options_type), intent(in) :: geom_options
  type(airfoil_type), intent(out) :: foilout
  
  integer :: i
  logical :: needs_cleanup

! Some things that need to be allocated for XFoil PANGEN

  needs_cleanup = .false.
  if (.not. allocated(W1)) then
    allocate(W1(6*IQX))
    allocate(W2(6*IQX))
    allocate(W3(6*IQX))
    allocate(W4(6*IQX))
    allocate(W5(6*IQX))
    allocate(W6(6*IQX))
    needs_cleanup = .true.
  end if

! Set some things that Xfoil may need to do paneling

  PI = 4.d0*atan(1.d0)
  HOPI = 0.5d0/PI
  QOPI = 0.25d0/PI
  SIG(:) = 0.d0
  NW = 0
  AWAKE = 0.d0
  LWDIJ = .false.
  LIPAN = .false.
  LBLINI = .false.
  WAKLEN = 1.d0
  GAM(:) = 0.d0
  SIGTE = 0.d0
  GAMTE = 0.d0
  SIGTE_A = 0.d0
  GAMTE_A = 0.d0
  SILENT_MODE = .TRUE.

! Set xfoil airfoil and paneling options

  call xfoil_set_airfoil(foilin)
  call xfoil_set_paneling(geom_options)

! Smooth paneling with PANGEN

  call PANGEN(.NOT. SILENT_MODE)

! Put smoothed airfoil coordinates into derived type

  foilout%npoint = geom_options%npan
  allocate(foilout%x(geom_options%npan))
  allocate(foilout%z(geom_options%npan))
  do i = 1, geom_options%npan
    foilout%x(i) = X(i)
    foilout%z(i) = Y(i)
  end do

! Deallocate memory that is not needed anymore

  if (needs_cleanup) then
    deallocate(W1)
    deallocate(W2)
    deallocate(W3)
    deallocate(W4)
    deallocate(W5)
    deallocate(W6)
  end if
  
end subroutine smooth_paneling

!=============================================================================80
!
! Subroutine to apply a flap deflection to the buffer airfoil and set it as the
! current airfoil.  For best results, this should be called after PANGEN.
!
!=============================================================================80
subroutine xfoil_apply_flap_deflection(xflap, yflap, y_flap_spec, degrees)

  use xfoil_inc
  use vardef, only : flap_connection, connection_apply, connection_radius
  
  double precision, intent(in) :: xflap, yflap, degrees
  character(3), intent(in) :: y_flap_spec
  
  double precision :: xradius
  integer :: int_option
  integer y_flap_spec_int

  if (y_flap_spec == 'y/c') then
    y_flap_spec_int = 0
  else
    y_flap_spec_int = 1
  end if

! Apply flap deflection
  
  !  Set radius 
  xradius = connection_radius
  
  !  Find the function to use
  if (trim(flap_connection) .EQ. 'sharp') int_option = 1
  if (trim(flap_connection) .EQ. 'smooth-top') int_option = 3
  if (trim(flap_connection) .EQ. 'smooth-bot') int_option = 2
  if (trim(flap_connection) .EQ. 'smooth') int_option = 4
  
  if (trim(connection_apply) .EQ. 'none') int_option = 1
  if (trim(connection_apply) .EQ. 'trip_wire') int_option = 1
  
  !  Apply the function
  select case (int_option)
    case (1)
      call FLAP(xflap, yflap, y_flap_spec_int, degrees)
    case (2)
      call FLAP_L(xflap, yflap, y_flap_spec_int, degrees, xradius)
    case (3)
      call FLAP_U(xflap, yflap, y_flap_spec_int, degrees, xradius)
    case (4)
      call FLAP_LU(xflap, yflap, y_flap_spec_int, degrees, xradius)
    case default
      write(*,*) 'Error in flap selection.'
      stop
  end select

end subroutine xfoil_apply_flap_deflection

!=============================================================================80
!
! Subroutine to get the max panel angle in an airfoil using
! Xfoil's CANG subroutine
!
!=============================================================================80
subroutine get_max_panel_angle(foilin, panel_angle)

  use vardef, only : airfoil_type

  type(airfoil_type), intent(in) :: foilin
  double precision, intent(out) :: panel_angle
  
  integer :: NB, INTMAX
  double precision :: ANGMAX
  double precision, dimension(foilin%npoint) :: XB, YB 
  
  NB = foilin%npoint
  XB(1:NB) = foilin%x
  YB(1:NB) = foilin%z
  
  CALL CANG(XB,YB,NB,0,INTMAX,ANGMAX)
  
  panel_angle = ANGMAX
  
end subroutine get_max_panel_angle

subroutine operating_point_analysis(op_mode, op_point,                         &
  xfoil_options_viscous_mode, xfoil_options_maxit,                             &
  op_lift, op_drag, op_moment, op_viscrms, op_alpha, op_xtrt, op_xtrb)
  
  use xfoil_inc
  
  character(7), intent(in) :: op_mode
  double precision, intent(in) :: op_point
  logical, intent(in) :: xfoil_options_viscous_mode
  integer, intent(in) :: xfoil_options_maxit
  double precision, intent(out) :: op_lift, op_drag, op_moment, op_viscrms
  double precision, intent(out), optional :: op_alpha, op_xtrt, op_xtrb


  if (op_mode == 'spec-al') then

    call xfoil_specal(op_point, xfoil_options_viscous_mode,                    &
                      xfoil_options_maxit, op_lift, op_drag, op_moment)

  elseif (op_mode == 'spec-cl') then

    call xfoil_speccl(op_point, xfoil_options_viscous_mode,                    &
                      xfoil_options_maxit, op_lift, op_drag, op_moment)

  else

    write(*,*)
    write(*,*) "Error in xfoil_driver: op_mode must be 'spec-al' or "//    &
      &        "'spec-cl'"
    write(*,*)
    stop

  end if

  ! Convergence check

  op_viscrms = RMSBL
  
  ! Get optional outputs

  if (present(op_alpha)) op_alpha = ALFA/DTOR
  if (present(op_xtrt)) op_xtrt = XOCTR(1)
  if (present(op_xtrb)) op_xtrb = XOCTR(2)

end subroutine operating_point_analysis
  
subroutine operating_point_analysis_with_init(op_mode, op_point,               &
  xfoil_options_viscous_mode, xfoil_options_maxit,                             &
  xfoil_options_init_number_points, xfoil_options_init_dist,                   &
  xfoil_options_init_initial_position,                                         &
  xfoil_options_init_al0, xfoil_options_init_cl0, xfoil_options_silent_mode,   &
  first_run_xfoil,                                                             &
  op_lift, op_drag, op_moment, op_viscrms, op_alpha, op_xtrt, op_xtrb)

  use xfoil_inc
  
  character(7), intent(in) :: op_mode
  double precision, intent(in) :: op_point
  logical, intent(in) :: xfoil_options_viscous_mode
  integer, intent(in) :: xfoil_options_maxit, xfoil_options_init_number_points
  character(20), intent(in) :: xfoil_options_init_dist
  double precision, intent(in) :: xfoil_options_init_initial_position
  double precision, intent(in) :: xfoil_options_init_al0, xfoil_options_init_cl0
  logical, intent(in) :: xfoil_options_silent_mode, first_run_xfoil

  double precision, intent(out) :: op_lift, op_drag, op_moment, op_viscrms
  double precision, intent(out), optional :: op_alpha, op_xtrt, op_xtrb

  integer :: naddpoints
  double precision, dimension(:), allocatable :: addpoints , addlift, adddrag, &
    addmoment, addalpha, addxtrt, addxtrb, addviscrms, addweight
  logical, dimension(:), allocatable :: addconverged
  double precision :: lower, upper
  character(30) :: text, text1, text2
  character(150) :: message

  integer :: j
  
  naddpoints = xfoil_options_init_number_points 
      
  allocate(addpoints(naddpoints), addlift(naddpoints),                   &
            adddrag(naddpoints), addmoment(naddpoints),                   &
            addalpha(naddpoints), addxtrt(naddpoints),                    &
            addxtrb(naddpoints), addviscrms(naddpoints),                  &
            addconverged(naddpoints), addweight(naddpoints))
        
  ! Loop all points in sequence init
  do j = 1, naddpoints
          
    ! Get limits
    lower = xfoil_options_init_initial_position
    upper = 1.0d0
          
    ! Linear sequence between 0 and 1
    addweight(j) = real(j-1,8)*(1.0d0)/real(naddpoints-1,8)
          
    ! Set distribution and scale to lower and upper limits
    if (trim(xfoil_options_init_dist) .EQ. 'sine') then
      addweight(j) = sin(addweight(j) * pi/2.0d0)
      addweight(j) = addweight(j)*(upper-lower)+lower
    elseif (trim(xfoil_options_init_dist) .EQ. 'linear') then
      addweight(j) = addweight(j)*(upper-lower)+lower
    else
      write(*,*)
      write(*,*) "Error in xfoil_driver: init_dist must be 'linear' or "//  &
      &          "'sine'"
      write(*,*)
      stop
    end if
          
    ! Run XFOIL
    if (op_mode == 'spec-al') then
      addpoints(j) = op_point - (1.0d0-addweight(j)) *        &
        abs(op_point-xfoil_options_init_al0) *                &
        sign(1.d0, op_point-xfoil_options_init_al0)
            
      call xfoil_specal(addpoints(j), xfoil_options_viscous_mode,          &
                  xfoil_options_maxit, addlift(j), adddrag(j), addmoment(j))

    elseif (op_mode == 'spec-cl') then
      addpoints(j) = op_point - (1.0d0-addweight(j)) *        &
        abs(op_point-xfoil_options_init_cl0) *                &
        sign(1.d0, op_point-xfoil_options_init_cl0)
            
      call xfoil_speccl(addpoints(j), xfoil_options_viscous_mode,          &
                  xfoil_options_maxit, addlift(j), adddrag(j), addmoment(j))

    else

      write(*,*)
      write(*,*) "Error in xfoil_driver: op_mode must be 'spec-al' or "//  &
      &          "'spec-cl'"
      write(*,*)
      stop

    end if

    ! Get optional outputs
        
    if (LVCONV) addconverged(j) = .true.
        
    if (present(op_alpha)) addalpha(j) = ALFA/DTOR
    if (present(op_xtrt)) addxtrt(j) = XOCTR(1)
    if (present(op_xtrb)) addxtrb(j) = XOCTR(2)
        
    addviscrms(j) = RMSBL
    
    if (RMSBL .GT. 1.0E-4) then
      LIPAN = .false.
      LBLINI = .false.
    end if
  
  end do
      
  ! Get op point values
  op_lift=addlift(naddpoints)
  op_drag=adddrag(naddpoints)
  op_moment=addmoment(naddpoints)
  op_alpha=addalpha(naddpoints)
  op_xtrt=addxtrt(naddpoints)
  op_xtrb=addxtrb(naddpoints)
  op_viscrms = addviscrms(naddpoints)
      
  if (.not. xfoil_options_silent_mode .or. first_run_xfoil) then
    do j = 1, naddpoints
      write(text,*) j
      write(text1,'(F8.4)') addalpha(j)
      write(text2,'(F8.6)') addlift(j)
      text = adjustl(text)
      if (addconverged(j) .AND. addviscrms(j) < 1.0D-04) then
        message = '   Initialization point '//trim(text)//' converged at '  &
          &//trim(text1)//' degrees with Cl of '//trim(text2)
      elseif (.not. addconverged(j) .OR. addviscrms(j) > 1.0D-04) then
        message = '   Initialization point '//trim(text)//' did not '//     &
          &'converge at '//trim(text1)//' degrees with Cl of '//trim(text2)
      end if
    write(*,*) trim(message)
    end do

  end if
      
  deallocate(addpoints, addlift, adddrag, addmoment, addalpha, addxtrt,    &
              addxtrb, addviscrms, addconverged, addweight)
  
end subroutine operating_point_analysis_with_init

subroutine operating_point_analysis_sequence(op_mode,                &
  op_search_op_start, op_search_op_end, op_search_op_step,                     &
  xfoil_options_viscous_mode, xfoil_options_maxit,                             &
  xfoil_options_silent_mode, first_run_xfoil,                                  &
  op_lift, op_drag, op_moment, op_viscrms, op_alpha, op_xtrt, op_xtrb)

  use xfoil_inc
  
  character(7), intent(in) :: op_mode
  logical, intent(in) :: xfoil_options_viscous_mode
  integer, intent(in) :: xfoil_options_maxit
  double precision, intent(in) :: op_search_op_start, op_search_op_end,        &
    op_search_op_step
  logical, intent(in) :: xfoil_options_silent_mode, first_run_xfoil
  
  double precision, intent(out) :: op_lift, op_drag, op_moment, op_viscrms
  double precision, intent(out), optional :: op_alpha, op_xtrt, op_xtrb

  integer :: naddpoints
  double precision, dimension(:), allocatable :: addpoints , addlift, adddrag, &
    addmoment, addalpha, addxtrt, addxtrb, addviscrms
  logical, dimension(:), allocatable :: addconverged
  character(30) :: text, text1, text2
  character(150) :: message
  
  integer :: j,k
  
  naddpoints = int((op_search_op_end-op_search_op_start)/op_search_op_step+0.5)+1
  if (naddpoints .LT. 1) then
    write(*,*)
    write(*,*) "Error in xfoil_driver op_point: start, end and step must "&
      &"result in a valid interval"
    write(*,*)
    stop
  end if
      
  allocate(addpoints(naddpoints), addlift(naddpoints), adddrag(naddpoints),&
            addmoment(naddpoints), addalpha(naddpoints),                    &
            addxtrt(naddpoints), addxtrb(naddpoints),                       &
            addviscrms(naddpoints), addconverged(naddpoints))
      
  do j = 1, naddpoints
        
    addpoints(j) = op_search_op_start + op_search_op_step*real(j-1,8)
      
    if (op_mode == 'spec-al') then

      call xfoil_specal(addpoints(j), xfoil_options_viscous_mode,          &
                  xfoil_options_maxit, addlift(j), adddrag(j), addmoment(j))

    elseif (op_mode == 'spec-cl') then

      call xfoil_speccl(addpoints(j), xfoil_options_viscous_mode,          &
                  xfoil_options_maxit, addlift(j), adddrag(j), addmoment(j))

    else

      write(*,*)
      write(*,*) "Error in xfoil_driver: op_mode must be 'spec-al' or "//  &
      &          "'spec-cl'"
      write(*,*)
      stop

    end if

    ! Get optional outputs
        
    if (LVCONV) addconverged(j) = .true.
        
    addalpha(j) = ALFA/DTOR
    addxtrt(j) = XOCTR(1)
    addxtrb(j) = XOCTR(2)
        
    addviscrms(j) = RMSBL
    
    if (RMSBL .GT. 1.0E-4) then
      LIPAN = .false.
      LBLINI = .false.
    end if
  end do
      
  k=0
  op_lift=0.0d0
  do j = 1, naddpoints
    if (.not. isnan(addlift(j)) .and. addviscrms(j) .LT. 1.0D-04) then
      if (addlift(j) .GT. op_lift) then
        op_lift=addlift(j)
        op_drag=adddrag(j)
        op_moment=addmoment(j)
        op_alpha=addalpha(j)
        op_xtrt=addxtrt(j)
        op_xtrb=addxtrb(j)
        op_viscrms = addviscrms(j)
        k=k+1
      end if
    end if
  end do
      
  if (k .eq. 0) then
    op_lift = -1.D+08
    op_drag = 1.D+08
    op_moment = -1.D+08
    op_viscrms = 1.D+08
  end if
      
  if (.not. xfoil_options_silent_mode .or. first_run_xfoil) then
    do j = 1, naddpoints
      write(text,*) j
      write(text1,'(F8.4)') addalpha(j)
      write(text2,'(F8.6)') addlift(j)
      text = adjustl(text)
      if (addconverged(j)) then
        message = '   Search point '//trim(text)//' converged at '            &
          &//trim(text1)//' degrees with Cl of '//trim(text2)
      elseif (.not. addconverged(j)) then
        message = '   Search point '//trim(text)//' did not converge at '     &
          &//trim(text1)//' degrees.'
      end if
      write(*,*) trim(message)
    end do
    write(text1,'(F8.4)') op_lift
    write(text2,'(F8.4)') op_alpha
    message = '   Best lift of '//trim(text1)//' at '//trim(text2)//' degrees'
    write(*,*) trim(message)

  end if
      
  deallocate(addpoints, addlift, adddrag, addmoment, addalpha, addxtrt,    &
              addxtrb, addviscrms, addconverged)
  
end subroutine operating_point_analysis_sequence

!=============================================================================80
!
! Subroutine to get Cl, Cd, Cm for an airfoil from Xfoil at given operating
! conditions.  Reynolds numbers and mach numbers should be specified for each
! operating point.  Additionally, op_mode determines whether each point is run
! at a constant alpha or cl - use 'spec-al' for specified alpha and 'spec-cl'
! for specified cl.  
! 
! Outputs:
!   alpha, Cl, Cd, Cm each operating point
!   viscrms: rms for viscous calculations (check for convergence)
!
!=============================================================================80
subroutine run_xfoil(foil, geom_options, operating_points, op_modes, op_search,&
                     use_previous_op,                                          &
                     reynolds_numbers, mach_numbers, use_flap, x_flap, y_flap, &
                     y_flap_spec, flap_degrees, xfoil_options, file_options,   &
                     lift, drag, moment, viscrms, alpha, xtrt, xtrb,           &
                     ncrit_per_point)

  use xfoil_inc
  use vardef,    only : airfoil_type, first_run_xfoil, op_search_type, output_prefix

  type(airfoil_type), intent(in) :: foil
  type(xfoil_geom_options_type), intent(in) :: geom_options
  double precision, dimension(:), intent(in) :: operating_points,              &
                                                reynolds_numbers, mach_numbers,&
                                                flap_degrees
  logical, dimension(:), intent(in) :: use_previous_op
  double precision, intent(in) :: x_flap, y_flap
  character(3), intent(in) :: y_flap_spec
  logical, intent(in) :: use_flap
  character(7), dimension(:), intent(in) :: op_modes
  type(op_search_type), intent(in) :: op_search
  type(xfoil_options_type), intent(in) :: xfoil_options
  type(xfoil_file_options_type), intent(in) :: file_options
  double precision, dimension(size(operating_points,1)), intent(out) :: lift,  &
                                                           drag, moment, viscrms
  double precision, dimension(size(operating_points,1)), intent(out),          &
                                                   optional :: alpha, xtrt, xtrb
  double precision, dimension(:), intent(in), optional :: ncrit_per_point
  
  integer :: start_cp_bl        !Replace file if point number is this

  integer :: i, j, noppoint
  integer :: current_search_point
 
  logical :: use_search
  logical, dimension(size(operating_points,1)) :: point_converged, point_fixed 
  character(30) :: text
  character(150) :: message

  start_cp_bl = 1
  current_search_point = 0
  
  if (.not. xfoil_options%silent_mode .or. first_run_xfoil) then
    write(*,*) 
    write(*,*) ' Analyzing aerodynamics using the XFOIL engine ...'
  end if 

  ! Check to make sure xfoil is initialized

  if (.not. allocated(AIJ)) then
    write(*,*) "Error: xfoil is not initialized!  Call xfoil_init() first."
    stop
  end if
  
  ! Set default Xfoil parameters

  call xfoil_defaults(xfoil_options)

  point_converged(:) = .false.
  point_fixed(:) = .false.

  noppoint = size(operating_points,1)

  ! Set paneling options

  call xfoil_set_paneling(geom_options)

  ! Set airfoil and smooth paneling

  if (.not. use_flap) then
    call xfoil_set_airfoil(foil)
    call PANGEN(.not. SILENT_MODE)
  end if

  ! Run xfoil for requested operating points

  lift(:) = 0.d0
  drag(:) = 0.d0
  moment(:) = 0.d0
  viscrms(:) = 0.d0

  ! Run xfoil for requested operating points

  if (op_search%noppoint .NE. 0) then
    current_search_point = 1 ! for search optimization start
  else
    current_search_point = 0
  end if
  
  run_oppoints: do i = 1, noppoint

    if (use_previous_op(i)) then
      j = i - 1
      lift(i)=lift(j)
      drag(i)=drag(j)
      moment(i)=moment(j)
      alpha(i)=alpha(j)
      xtrt(i)=xtrt(j)
      xtrb(i)=xtrb(j)
      viscrms(i) = viscrms(j)
      point_converged(i) = point_converged(j)
      point_fixed(i) = point_fixed(j)
      
      go to 800
    end if
    
    ! Reset airfoil, smooth paneling, and apply flap deflection
    
    if (use_flap) then
      call xfoil_set_airfoil(foil)
      call PANGEN(.not. SILENT_MODE)
      if (file_options%write_all_airfoils .and. i .eq. 1) then
        open(unit=1300, file=trim(output_prefix)//'_base.dat', status='replace')
        write(1300,'(A)') trim(output_prefix)//'_base'
        do j = 1, N
          write(1300,'(2F12.6)') X(j), Y(j)
        end do
        close(1300)
      end if
      call xfoil_apply_flap_deflection(x_flap, y_flap, y_flap_spec,            &
                                       flap_degrees(i))
      if (file_options%write_all_airfoils) then
        write(text,'(I0)') i
        open(unit=1300, file=trim(output_prefix)//'_at_op_'//trim(text)//'.dat', status='replace')
        write(1300,'(A)') trim(output_prefix)//'_at_op_'//trim(text)
        do j = 1, N
          write(1300,'(2F12.6)') X(j), Y(j)
        end do
        close(1300)
      end if
    end if

    REINF1 = reynolds_numbers(i)
    call MINFSET(mach_numbers(i))

    if (xfoil_options%reinitialize) then
      LIPAN = .false.
      LBLINI = .false.
    end if

    ! Set compressibility parameters from MINF

    CALL COMSET

    ! Set ncrit per point

    if (present(ncrit_per_point)) ACRIT = ncrit_per_point(i)
    
    ! Check if search optimization type in point is to be used
    if (op_search%noppoint .EQ. 0 .OR.                                         &
      current_search_point .GT. op_search%noppoint) then
      use_search = .false.
    else
      if (i .ne. op_search%oppoints(current_search_point)) then
        use_search = .false.
      else
        use_search = .true.
      end if
    end if
    
    ! Analyse operating point
    if (.not. use_search) then
      
      ! Check if op point is to be analysed without initialization
      if (trim(xfoil_options%init_type) .EQ. 'never' .or.                      &
        trim(xfoil_options%init_type) .EQ. 'unconverged') then
      
        call operating_point_analysis(op_modes(i), operating_points(i),        &
          xfoil_options%viscous_mode, xfoil_options%maxit, lift(i), drag(i),   &
          moment(i), viscrms(i), alpha(i), xtrt(i), xtrb(i))
            
        if (.not. isnan(viscrms(i))) then
          if (viscrms(i).LT.1.0E-4) then
            point_converged(i) = .true.
          else
            point_converged(i) = .false.
          end if
        end if
        
      end if
      
      ! Check if op point is to be reanalysed with initialization 
      if (xfoil_options%viscous_mode .and. .not. point_converged(i) .and.      &
          trim(xfoil_options%init_type) .EQ. 'unconverged') then
        
        call operating_point_analysis_with_init(op_modes(i),                   &
          operating_points(i), xfoil_options%viscous_mode, xfoil_options%maxit,&
          xfoil_options%init_number_points, xfoil_options%init_dist,           &
          xfoil_options%init_initial_position,                                 &
          xfoil_options%init_al0, xfoil_options%init_cl0,                      &
          xfoil_options%silent_mode, first_run_xfoil,                          &
          lift(i), drag(i), moment(i), viscrms(i), alpha(i), xtrt(i), xtrb(i))
        
        if (.not. isnan(viscrms(i))) then
          if (viscrms(i).LT.1.0E-4) then
            point_fixed(i) = .true.
          else
            point_fixed(i) = .false.
          end if
        end if

      end if
          
      ! Check if op point is to be analysed with initialization   
      if (trim(xfoil_options%init_type) .EQ. 'always') then
        
        call operating_point_analysis_with_init(op_modes(i),                   &
          operating_points(i), xfoil_options%viscous_mode, xfoil_options%maxit,&
          xfoil_options%init_number_points, xfoil_options%init_dist,           &
          xfoil_options%init_initial_position,                                 &
          xfoil_options%init_al0, xfoil_options%init_cl0,                      &
          xfoil_options%silent_mode, first_run_xfoil,                          &
          lift(i), drag(i), moment(i), viscrms(i), alpha(i), xtrt(i), xtrb(i))
        
        if (.not. isnan(viscrms(i))) then
          if (viscrms(i).LT.1.0E-4) then
            point_converged(i) = .true.
          else
            point_converged(i) = .false.
          end if
        end if
        
      end if
    
    else

      call operating_point_analysis_sequence(op_modes(i), &
        op_search%op_start(current_search_point),                              &
        op_search%op_end(current_search_point),                                &
        op_search%op_step(current_search_point),                               &
        xfoil_options%viscous_mode, xfoil_options%maxit,                       &
        xfoil_options%silent_mode, first_run_xfoil,                            &
        lift(i), drag(i), moment(i), viscrms(i), alpha(i), xtrt(i), xtrb(i))
      
      current_search_point = current_search_point + 1
      
      if (.not. isnan(viscrms(i))) then
        if (viscrms(i).LT.1.0E-4) then
          point_converged(i) = .true.
        else
          point_converged(i) = .false.
        end if
      end if
      
    end if

800 continue
  
    if (.not. use_search) then
      if (file_options%polar) call write_polar(i, file_options%design_number,    &
        x_flap, flap_degrees(i))
      if (file_options%cp) call write_polar_cp(i, file_options%design_number,    &
        x_flap, flap_degrees(i), start_cp_bl)
      if (file_options%bl) call write_polar_bl(i, file_options%design_number,    &
        x_flap, flap_degrees(i), start_cp_bl)
    else
      if (file_options%polar) call write_polar_search(i, file_options%design_number,    &
        x_flap, flap_degrees(i), alpha(i), lift(i), moment(i),drag(i), drag(i), xtrt(i), xtrb(i))
      if (i .EQ. 1) then
        start_cp_bl = 2 ! start op for cp and bl can not be search type
      end if
    end if
  end do run_oppoints
    
    ! Final check for NaNs

  do i = 1, noppoint
    if (isnan(lift(i))) then
      lift(i) = -1.D+08
      viscrms(i) = 1.D+08
    end if
    if (isnan(drag(i))) then
      drag(i) = 1.D+08
      viscrms(i) = 1.D+08
    end if
    if (isnan(moment(i))) then
      moment(i) = -1.D+08
      viscrms(i) = 1.D+08
    end if
    if (isnan(viscrms(i))) then
      viscrms(i) = 1.D+08
    end if
  end do

  ! Print warnings about unconverged points

  if (.not. xfoil_options%silent_mode .or. first_run_xfoil) then

    write(*,*)

    do i = 1, noppoint
  
      write(text,*) i
      text = adjustl(text)
  
      if (point_converged(i) .AND. viscrms(i) < 1.0D-04) then
        message = '   Operating point '//trim(text)//' converged.'
      elseif (.not. point_converged(i) .and. point_fixed(i) .AND. viscrms(i) < 1.0D-04) then
        message = '   Operating point '//trim(text)//' initially did not '//      &
                  'converge but was fixed.'
      elseif ((.not. point_converged(i) .and. .not. point_fixed(i)) .OR. viscrms(i) > 1.0D-04) then
        message = '   Operating point '//trim(text)//' initially did not '//      &
                  'converge and was not fixed.'
      end if
  
      write(*,*) trim(message)
    end do
  end if

end subroutine run_xfoil

!=============================================================================80
!
! Runs Xfoil at a specified angle of attack
! Assumes airfoil geometry, reynolds number, and mach number have already been 
! set in Xfoil.
!
!=============================================================================80
subroutine xfoil_specal(angle_of_attack, viscous_mode, maxit, lift, drag,      &
                        moment)

  use xfoil_inc

  double precision, intent(in) :: angle_of_attack
  logical, intent(in) :: viscous_mode
  integer, intent(in) :: maxit
  double precision, intent(out) :: lift, drag, moment

! Inviscid calculations for specified angle of attack

  LALFA = .TRUE.
  ALFA = angle_of_attack*DTOR
  call SPECAL
  if (abs(ALFA-AWAKE) .GT. 1.0E-5) LWAKE  = .false.
  if (abs(ALFA-AVISC) .GT. 1.0E-5) LVCONV = .false.
  if (abs(MINF-MVISC) .GT. 1.0E-5) LVCONV = .false.

! Viscous calculations (if requested)

  if (viscous_mode) call VISCAL(maxit)

! Outputs

  lift = CL
  moment = CM
  if (viscous_mode) then
    drag = CD
  else
    drag = CDP
  end if

end subroutine xfoil_specal
                        
!=============================================================================80
!
! Writes XFOIL polars using xfoil format
!
!=============================================================================80
subroutine write_polar_xfoil(point_number)

  use xfoil_inc
  use vardef, only : output_prefix, airfoil_file

  integer, intent(inout) :: point_number
  
  double precision :: alpha, lift, moment,drag, drag_pressure, xtrt, xtrb
  integer :: i
  character(100) :: text1

! Get values to create polar
  alpha = ALFA/DTOR
  lift = CL
  moment = CM
  drag = CD
  drag_pressure = CDP
  xtrt = XOCTR(1)
  xtrb = XOCTR(2)
  
! Replace file if point_number is 1

  if (point_number .EQ. 1) then
    open(unit=321, file=trim(output_prefix)//'_polar_file.log',                &
      status='replace')

    write(321,'(A)') '  '
    write(321,'(A)') '       XFOIL_ONLY as part of XOPTFOIL'
    write(321,'(A)') '  '
    write(text1,*) airfoil_file
    text1 = adjustl(text1)
    i=len_trim(text1)
    write(321,'(A)') ' Calculated polar for: '//text1(1:i-4)
    write(321,'(A)') '  '
    write(321,'(A)') ' 1 1 Reynolds number fixed          Mach number fixed         '
    write(321,'(A)') '  '
    write(321,'(A9,F7.3,A13,F7.3,A11)') ' xtrf =  ',XSTRIP(1),' (top)       ', &
      XSTRIP(2),' (bottom)  '
    write(321,'(A9,F7.3,A13,F7.3,A17,F7.3)') ' Mach =  ',MINF1,'     Re =    ',&
      REINF1/10.0**6,' e 6     Ncrit = ',ACRIT
    write(321,'(A)') '  '
    write(321,'(A)') '   alpha    CL        CD       CDp       CM    '//       &
      & ' Top_Xtr  Bot_Xtr'
    write(321,'(A)') '  ------ -------- --------- --------- -------- '//       &
      & '-------- --------'

    close(321)
  end if
  
! Write values to polar file
  open(unit=300, file=trim(output_prefix)//'_polar_file.log', status='old',    &
    position='append')
  write(300,1000) alpha, lift, drag, drag_pressure, moment, xtrt, xtrb
1000 format(F8.3,F9.4,F10.5,F10.5,F9.4,F9.4,F9.4)
  close(300)  
  
end subroutine write_polar_xfoil


!=============================================================================80
!
! Writes XFOIL polars
!
!=============================================================================80
subroutine write_polar(point_number, design_number, x_flap, flap_degrees)

  use xfoil_inc
  use vardef, only : output_prefix

  integer, intent(in) :: point_number
  integer, intent(in) :: design_number
  double precision, intent(in) :: x_flap, flap_degrees
  
  double precision :: alpha, lift, moment,drag, drag_pressure, xtrt, xtrb

  character(100) :: text

! Get values to create polar
  alpha = ALFA/DTOR
  lift = CL
  moment = CM
  drag = CD
  drag_pressure = CDP
  xtrt = XOCTR(1)
  xtrb = XOCTR(2)
  
  write(text,*) design_number
  text = adjustl(text)
  
  if (point_number .EQ. 1) then
    write(*,*) "  Writing polars for design number "//trim(text)//             &
               " to file "//trim(output_prefix)//"_design_polars.dat"//" ..."
  end if

! Replace file if point_number is 1

  if (design_number .EQ. 0 .AND. point_number .EQ. 1) then
    open(unit=321, file=trim(output_prefix)//'_design_polars.dat',             &
      status='replace')
    write(321,'(A)') ' Polar file '
    write(321,'(A)') '   alpha    CL        CD       CDp       CM    '//       &
      & ' Top_Xtr  Bot_Xtr Number     Re/10^6    Mach   Top_Xtrip Bot_Xtrip'// &
      & ' Ncrit FxHinge  Fdefl'
    write(321,'(A)') 'zone t="Seed airfoil polar"'

    close(321)
  elseif (point_number .EQ. 1) then
    open(unit=310, file=trim(output_prefix)//'_design_polars.dat',             &
      status='old', position='append')
    write(310,'(A)') 'zone t="Polars", SOLUTIONTIME='//trim(text)
    close(310)
  end if
  
! Write values to polar file
  open(unit=300, file=trim(output_prefix)//'_design_polars.dat', status='old', &
    position='append')
  write(300,1000) alpha, lift, drag, drag_pressure, moment, xtrt, xtrb,        &
    point_number, REINF1/10**6, MINF1, XSTRIP(1), XSTRIP(2), ACRIT, x_flap,    &
    flap_degrees
    1000 format(F8.3,F9.4,F10.5,F10.5,F9.4,F9.4,F9.4,I7,F13.6,F9.4,F10.4,F10.4,&
         &      F6.2, F7.3, F8.2)
  close(300)  
  
end subroutine write_polar          

!=============================================================================80
!
! Writes XFOIL polars for search type op
!
!=============================================================================80
subroutine write_polar_search(point_number, design_number, x_flap, flap_degrees,&
  alpha, lift, moment,drag, drag_pressure, xtrt, xtrb)

  use xfoil_inc
  use vardef, only : output_prefix

  integer, intent(in) :: point_number
  integer, intent(in) :: design_number
  double precision, intent(in) :: x_flap, flap_degrees
  
  double precision, intent(in) :: alpha, lift, moment,drag, drag_pressure, xtrt, xtrb

  character(100) :: text

! Get values to create polar
  !alpha = ALFA/DTOR
  !lift = CL
  !moment = CM
  !drag = CD
  !drag_pressure = CDP
  !xtrt = XOCTR(1)
  !xtrb = XOCTR(2)
  
  write(text,*) design_number
  text = adjustl(text)
  
  if (point_number .EQ. 1) then
    write(*,*) "  Writing polars for design number "//trim(text)//             &
               " to file "//trim(output_prefix)//"_design_polars.dat"//" ..."
  end if

! Replace file if point_number is 1

  if (design_number .EQ. 0 .AND. point_number .EQ. 1) then
    open(unit=321, file=trim(output_prefix)//'_design_polars.dat',             &
      status='replace')
    write(321,'(A)') ' Polar file '
    write(321,'(A)') '   alpha    CL        CD       CDp       CM    '//       &
      & ' Top_Xtr  Bot_Xtr Number     Re/10^6    Mach   Top_Xtrip Bot_Xtrip'// &
      & ' Ncrit FxHinge  Fdefl'
    write(321,'(A)') 'zone t="Seed airfoil polar"'

    close(321)
  elseif (point_number .EQ. 1) then
    open(unit=310, file=trim(output_prefix)//'_design_polars.dat',             &
      status='old', position='append')
    write(310,'(A)') 'zone t="Polars", SOLUTIONTIME='//trim(text)
    close(310)
  end if
  
! Write values to polar file
  open(unit=300, file=trim(output_prefix)//'_design_polars.dat', status='old', &
    position='append')
  write(300,1000) alpha, lift, drag, drag_pressure, moment, xtrt, xtrb,        &
    point_number, REINF1/10**6, MINF1, XSTRIP(1), XSTRIP(2), ACRIT, x_flap,    &
    flap_degrees
    1000 format(F8.3,F9.4,F10.5,F10.5,F9.4,F9.4,F9.4,I7,F13.6,F9.4,F10.4,F10.4,&
         &      F6.2, F7.3, F8.2)
  close(300)  
  
end subroutine write_polar_search   

!=============================================================================80
!
! Writes XFOIL polars and cp distribution
! Based on CPDUMP subroutine from XFOIL V6.99
!
!=============================================================================80
subroutine write_polar_cp(point_number, design_number, x_flap, flap_degrees,   &
  start_number)

  use xfoil_inc
  use vardef, only : output_prefix

  integer, intent(in) :: point_number
  integer, intent(in) :: design_number
  double precision, intent(in) :: x_flap, flap_degrees
  integer, intent(in) :: start_number
  
  double precision :: alpha, lift, moment,drag, drag_pressure, xtrt, xtrb
  !double precision :: BETA, BFAC, CPINC, DEN, CPCOM
  integer :: i
  character(100) :: text

! Get values to create polar
  alpha = ALFA/DTOR
  lift = CL
  moment = CM
  drag = CD
  drag_pressure = CDP
  xtrt = XOCTR(1)
  xtrb = XOCTR(2)
  
! Replace file if point_number is start_number

if (design_number .EQ. 0 .AND. point_number .EQ. start_number) then
    open(unit=321, file=trim(output_prefix)//'_design_cp.dat',                 &
      status='replace')
    write(321,'(A)') ' Polar file and Cp distributions'
    write(321,'(A)') 'zone t="Seed airfoil polar"'

    close(321)
  elseif (point_number .EQ. start_number) then
    open(unit=310, file=trim(output_prefix)//'_design_cp.dat', status='old',   &
      position='append')
    write(text,*) design_number
    text = adjustl(text)
    write(310,'(A)') 'zone t="Polars", SOLUTIONTIME='//trim(text)
    close(310)
  end if
  
  
! Write values to polar file
  open(unit=300, file=trim(output_prefix)//'_design_cp.dat', status='old',     &
    position='append')
  
  WRITE(300,'(A)')
  write(300,'(A)') '   alpha    CL        CD       CDp       CM    '//         &
    & ' Top_Xtr  Bot_Xtr Number     Re/10^6    Mach   Top_Xtrip Bot_Xtrip'//   &
    & ' Ncrit FxHinge  Fdefl NPanel'
    
  write(300,1000) alpha, lift, drag, drag_pressure, moment, xtrt, xtrb,        &
    point_number, REINF1/10**6, MINF1, XSTRIP(1), XSTRIP(2), ACRIT, x_flap,    &
    flap_degrees, N
    1000 format(F8.3,F9.4,F10.5,F10.5,F9.4,F9.4,F9.4,I7,F13.6,F9.4,F10.4,F10.4,&
         &      F6.2, F7.3, F8.2, I7)
  
  ! Write Cp distribution
  WRITE(300,'(A)')
  WRITE(300,'(A)') '     x        y        Cpi      Cpv  '

  !CALL COMSET
  !
  !BETA = SQRT(1.0 - MINF**2)
  !BFAC = 0.5*MINF**2 / (1.0 + BETA)

  DO I=1, N
    !CPINC = 1.0 - (GAM(I)/QINF)**2
    !DEN = BETA + BFAC*CPINC
    !CPCOM = CPINC / DEN
        
    WRITE(300,8500) X(I), Y(I), CPI(I), CPV(I)
8500   FORMAT(1X,4F9.5)
  end do
  WRITE(300,'(A)')

  close(300)  
  
end subroutine write_polar_cp          

!=============================================================================80
!
! Writes XFOIL polars and bl distributions
! Based on BLDUMP subroutine from XFOIL V6.99
!
!=============================================================================80
subroutine write_polar_bl(point_number, design_number, x_flap, flap_degrees,   &
  start_number)

  use xfoil_inc
  use xbl_inc
  
  use vardef, only : output_prefix

  integer, intent(in) :: point_number
  integer, intent(in) :: design_number
  double precision, intent(in) :: x_flap, flap_degrees
  integer, intent(in) :: start_number
  
  double precision :: alpha, lift, moment,drag, drag_pressure, xtrt, xtrb
  double precision :: AMSQ, CDIS, CF, CT, DS, dummy, H, HK, TH, UE, UI
  integer :: i, IBL, IS
  character(100) :: text

! Get values to create polar
  alpha = ALFA/DTOR
  lift = CL
  moment = CM
  drag = CD
  drag_pressure = CDP
  xtrt = XOCTR(1)
  xtrb = XOCTR(2)
  
! Replace file if point_number is 1

if (design_number .EQ. 0 .AND. point_number .EQ. start_number) then
    open(unit=321, file=trim(output_prefix)//'_design_bl.dat',                 &
      status='replace')
    write(321,'(A)') ' Polar file and BL distributions'
    write(321,'(A)') 'zone t="Seed airfoil polar"'

    close(321)
  elseif (point_number .EQ. start_number) then
    open(unit=310, file=trim(output_prefix)//'_design_bl.dat', status='old',   &
      position='append')
    write(text,*) design_number
    text = adjustl(text)
    write(310,'(A)') 'zone t="Polars", SOLUTIONTIME='//trim(text)
    close(310)
  end if
  
! Write values to polar file
  open(unit=300, file=trim(output_prefix)//'_design_bl.dat', status='old',     &
    position='append')
  
  WRITE(300,'(A)')
  write(300,'(A)') '   alpha    CL        CD       CDp       CM    '//         &
    & ' Top_Xtr  Bot_Xtr Number     Re/10^6    Mach   Top_Xtrip Bot_Xtrip'//   &
    & ' Ncrit FxHinge  Fdefl  NP+NW'
    
  write(300,1000) alpha, lift, drag, drag_pressure, moment, xtrt, xtrb,        &
    point_number, REINF1/10**6, MINF1, XSTRIP(1), XSTRIP(2), ACRIT, x_flap,    &
    flap_degrees, N+NW
    1000 format(F8.3,F9.4,F10.5,F10.5,F9.4,F9.4,F9.4,I7,F13.6,F9.4,F10.4,F10.4,&
         &      F6.2, F7.3, F8.2, I7)
  
  ! Write BL distribution
  write(300,'(A)')
  WRITE(300,'(A)') '     s        x        y        Nx       Ny     Ue/Vinf'// &
  & '    Dstar     Theta      Cf         H        CD       CT'

  CALL COMSET
  DO I=1, N
    IS = 1
    IF(GAM(I) .LT. 0.0) IS = 2

    IF(LIPAN .AND. LVISC) THEN
      IF(IS.EQ.1) THEN
        IBL = IBLTE(IS) - I + 1
      ELSE
        IBL = IBLTE(IS) + I - N
      ENDIF
      DS = DSTR(IBL,IS)
      TH = THET(IBL,IS)
      CF =  TAU(IBL,IS)/(0.5*QINF**2)
      IF(TH.EQ.0.0) THEN
        H = 1.0
      ELSE
        H = DS/TH
      ENDIF
    ELSE
      DS = 0.
      TH = 0.
      CF = 0.
      H = 1.0
    ENDIF
    UE = (GAM(I)/QINF)*(1.0-TKLAM) / (1.0 - TKLAM*(GAM(I)/QINF)**2)
    AMSQ = UE*UE*HSTINV / (GAMM1*(1.0 - 0.5*UE*UE*HSTINV))
    CALL HKIN( H, AMSQ, HK, DUMMY, DUMMY)

    CDIS = DIS(IBL,IS)/QINF**3
    CT = CTAU(IBL,IS)
    WRITE(300,8500) S(I), X(I), Y(I), NX(I), NY(I), UE, DS, TH, CF, HK, CDIS, CT
  end do


8500   FORMAT(1X,6F9.5,3F10.6,F10.3,2F10.6)

  IF(LWAKE) THEN
    IS = 2
    DO I=N+1, N+NW
      IBL = IBLTE(IS) + I - N
      DS = DSTR(IBL,IS)
      TH = THET(IBL,IS)
      H = DS/TH
      CF = 0.
      UI = UEDG(IBL,IS)
      UE = (UI/QINF)*(1.0-TKLAM) / (1.0 - TKLAM*(UI/QINF)**2)
      AMSQ = UE*UE*HSTINV / (GAMM1*(1.0 - 0.5*UE*UE*HSTINV))
      CALL HKIN( H, AMSQ, HK, DUMMY, DUMMY)

      CDIS = DIS(IBL,IS)/QINF**3
      CT = CTAU(IBL,IS)
      WRITE(300,8500) S(I), X(I), Y(I), UE, DS, TH, CF, HK, CDIS, CT
    end do
  ENDIF
  write(300,'(A)')
  close(300)  
  
end subroutine write_polar_bl          

!=============================================================================80
!
! Runs Xfoil at a specified lift coefficient
! Assumes airfoil geometry, reynolds number, and mach number have already been 
! set in Xfoil.
!
!=============================================================================80
subroutine xfoil_speccl(cl_spec, viscous_mode, maxit, lift, drag, moment)

  use xfoil_inc

  double precision, intent(in) :: cl_spec
  logical, intent(in) :: viscous_mode
  integer, intent(in) :: maxit
  double precision, intent(out) :: lift, drag, moment

! Inviscid calculations for specified lift coefficient

  LALFA = .FALSE.
  ALFA = 0.d0
  CLSPEC = cl_spec
  call SPECCL
  if (abs(ALFA-AWAKE) .GT. 1.0E-5) LWAKE  = .false.
  if (abs(ALFA-AVISC) .GT. 1.0E-5) LVCONV = .false.
  if (abs(MINF-MVISC) .GT. 1.0E-5) LVCONV = .false.

! Viscous calculations (if requested)

  if (viscous_mode) call VISCAL(maxit)

! Outputs

  lift = CL
  moment = CM
  if (viscous_mode) then
    drag = CD
  else
    drag = CDP
  end if

end subroutine xfoil_speccl

!=============================================================================80
!
! Allocates xfoil variables that may be too big for the stack in OpenMP
!
!=============================================================================80
subroutine xfoil_init()

  use xfoil_inc

! Allocate variables that may be too big for the stack in OpenMP

  allocate(AIJ(IQX,IQX))
  allocate(BIJ(IQX,IZX))
  allocate(DIJ(IZX,IZX))
  allocate(CIJ(IWX,IQX))
  allocate(IPAN(IVX,ISX))
  allocate(ISYS(IVX,ISX))
  allocate(W1(6*IQX))
  allocate(W2(6*IQX))
  allocate(W3(6*IQX))
  allocate(W4(6*IQX))
  allocate(W5(6*IQX))
  allocate(W6(6*IQX))
  allocate(VTI(IVX,ISX))
  allocate(XSSI(IVX,ISX))
  allocate(UINV(IVX,ISX))
  allocate(UINV_A(IVX,ISX))
  allocate(UEDG(IVX,ISX))
  allocate(THET(IVX,ISX))
  allocate(DSTR(IVX,ISX))
  allocate(CTAU(IVX,ISX))
  allocate(MASS(IVX,ISX))
  allocate(TAU(IVX,ISX))
  allocate(DIS(IVX,ISX))
  allocate(CTQ(IVX,ISX))
  allocate(DELT(IVX,ISX))
  allocate(TSTR(IVX,ISX))
  allocate(USLP(IVX,ISX))
  allocate(VM(3,IZX,IZX))
  allocate(VA(3,2,IZX))
  allocate(VB(3,2,IZX))
  allocate(VDEL(3,2,IZX))

end subroutine xfoil_init

!=============================================================================80
!
! Initializes xfoil variables
!
!=============================================================================80
subroutine xfoil_defaults(xfoil_options)

  use xfoil_inc

  type(xfoil_options_type), intent(in) :: xfoil_options

  N = 0
  SILENT_MODE = xfoil_options%silent_mode
  PI = 4.d0*atan(1.d0)
  HOPI = 0.5d0/PI
  QOPI = 0.25d0/PI
  DTOR = PI/180.d0
  QINF = 1.d0
  SIG(:) = 0.d0
  QF0(:) = 0.d0
  QF1(:) = 0.d0
  QF2(:) = 0.d0
  QF3(:) = 0.d0
  NW = 0
  RETYP = 1
  MATYP = 1
  GAMMA = 1.4d0
  GAMM1 = GAMMA - 1.d0
  XCMREF = 0.25d0
  YCMREF = 0.d0
  LVISC = xfoil_options%viscous_mode
  AWAKE = 0.d0
  AVISC = 0.d0
  ITMAX = xfoil_options%maxit
  LWDIJ = .false.
  LIPAN = .false.
  LBLINI = .false.
  ACRIT = xfoil_options%ncrit
  IDAMP = 0
  XSTRIP(1) = xfoil_options%xtript
  XSTRIP(2) = xfoil_options%xtripb
  VACCEL = xfoil_options%vaccel
  WAKLEN = 1.d0
  PSIO = 0.d0
  GAMU(:,:) = 0.d0
  GAM(:) = 0.d0
  SIGTE = 0.d0
  GAMTE = 0.d0
  SIGTE_A = 0.d0
  GAMTE_A = 0.d0
  APANEL(:) = 0.d0

! Set boundary layer calibration parameters

  call BLPINI

end subroutine xfoil_defaults

!=============================================================================80
!
! Sets airfoil for xfoil
!
!=============================================================================80
subroutine xfoil_set_airfoil(foil)

  use xfoil_inc, only : XB, YB, NB
  use vardef,    only : airfoil_type

  type(airfoil_type), intent(in) :: foil

  NB = foil%npoint
  XB(1:NB) = foil%x
  YB(1:NB) = foil%z

end subroutine xfoil_set_airfoil

!=============================================================================80
!
! Sets xfoil paneling options
!
!=============================================================================80
subroutine xfoil_set_paneling(geom_options)

  use xfoil_inc, only : NPAN, CVPAR, CTERAT, CTRRAT, XSREF1, XSREF2, XPREF1,   &
                        XPREF2

  type(xfoil_geom_options_type), intent(in) :: geom_options

  NPAN = geom_options%npan
  CVPAR = geom_options%cvpar
  CTERAT = geom_options%cterat
  CTRRAT = geom_options%ctrrat
  XSREF1 = geom_options%xsref1
  XSREF2 = geom_options%xsref2
  XPREF1 = geom_options%xpref1
  XPREF2 = geom_options%xpref2
  
end subroutine xfoil_set_paneling

!=============================================================================80
!
! Deallocates memory in xfoil
!
!=============================================================================80
subroutine xfoil_cleanup()

  use xfoil_inc

! Deallocate variables

  deallocate(AIJ)
  deallocate(BIJ)
  deallocate(DIJ)
  deallocate(CIJ)
  deallocate(IPAN)
  deallocate(ISYS)
  deallocate(W1)
  deallocate(W2)
  deallocate(W3)
  deallocate(W4)
  deallocate(W5)
  deallocate(W6)
  deallocate(VTI)
  deallocate(XSSI)
  deallocate(UINV)
  deallocate(UINV_A)
  deallocate(UEDG)
  deallocate(THET)
  deallocate(DSTR)
  deallocate(CTAU)
  deallocate(MASS)
  deallocate(TAU)
  deallocate(DIS)
  deallocate(CTQ)
  deallocate(DELT)
  deallocate(TSTR)
  deallocate(USLP)
  deallocate(VM)
  deallocate(VA)
  deallocate(VB)
  deallocate(VDEL)

end subroutine xfoil_cleanup

!=============================================================================80
!
! Gets thickness and camber information for the current airfoil
!
!=============================================================================80
subroutine xfoil_geometry_info(maxt, xmaxt, maxc, xmaxc)

  use xfoil_inc, only : THICKB, XTHICKB, CAMBR, XCAMBR

  double precision, intent(out) :: maxt, xmaxt, maxc, xmaxc

  maxt = THICKB
  xmaxt = XTHICKB
  maxc = CAMBR
  xmaxc = XCAMBR 

end subroutine xfoil_geometry_info

end module xfoil_driver
