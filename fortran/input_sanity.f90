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
    
module input_sanity

  implicit none

  contains

!=============================================================================80
!
! Checks that the seed airfoil passes all constraints, sets scale factors for
! objective functions at each operating point, and optionally sets the
! minimum allowable pitching moment, lift and maximum drag.
!
!=============================================================================80
subroutine check_seed()

  use vardef
  use math_deps,          only : interp_vector, nu_curvature, derv1f1, derv1b1,&
                                 spline_interp_z, spline_interp_t
  use xfoil_driver,       only : run_xfoil, get_max_panel_angle
  use airfoil_evaluation, only : xfoil_options, xfoil_geom_options, file_options
  use performance_evaluation

  double precision, dimension(noppoint) :: lift, drag, moment, viscrms, alpha, &
                                           xtrt, xtrb
  double precision, dimension(:), allocatable :: x_interp, thickness, camber
  double precision, dimension(:), allocatable :: zt_interp, zb_interp
  double precision, dimension(size(xseedt,1)) :: curvt
  double precision, dimension(size(xseedb,1)) :: curvb
  double precision, dimension(naddthickconst) :: add_thickvec
  double precision :: penaltyval, tegap, gapallow, maxthick, heightfactor, maxcamb, panel_angle
  double precision :: panang1, panang2, maxpanang, minpanang, curv1, curv2
  double precision :: checkval, len1, len2, growth1, growth2, xtrans, ztrans
  double precision, dimension(noppoint) :: actual_flap_degrees
  double precision :: pi, te_angle, actual_min_te_angle, max_growth_top,       &
                      max_growth_bot
  integer :: i, nptt, nptb, nreversalst, nreversalsb, nptint, flap_idi
  integer :: ninit, nend, ninit1, nend1
  character(30) :: text, text2
  character(15) :: opt_type
  logical :: addthick_violation, side

  penaltyval = 0.d0
  pi = acos(-1.d0)
  nptt = size(xseedt,1)
  nptb = size(xseedb,1)

  ! Format coordinates in a single loop in derived type. Also remove translation
  ! and scaling to ensure Cm_x=0.25 doesn't change.

  do i = 1, nptt
    curr_foil%x(i) = xseedt(nptt-i+1)!/foilscale - xoffset
    curr_foil%z(i) = zseedt(nptt-i+1)!/foilscale - zoffset
  end do
  do i = 1, nptb-1
    curr_foil%x(i+nptt) = xseedb(i+1)!/foilscale - xoffset
    curr_foil%z(i+nptt) = zseedb(i+1)!/foilscale - zoffset
  end do
  
  write(*,*)
  write(*,*) ' Checking to make sure seed airfoil passes all constraints ...'
  write(*,*)

  ! Interpolate either bottom surface to top surface x locations or vice versa
  ! to determine thickness

  if (xseedt(nptt) <= xseedb(nptb)) then
    allocate(x_interp(nptt))
    allocate(zt_interp(nptt))
    allocate(zb_interp(nptt))
    allocate(thickness(nptt))
    allocate(camber(nptt))
    nptint = nptt
    side=.false.
    call spline_interp_z(curr_foil%npoint,curr_foil%x,curr_foil%z,xseedt,zb_interp(1:nptt),side)
    !call interp_vector(xseedb, zseedb, xseedt, zb_interp)
    x_interp = xseedt
    zt_interp = zseedt
  else
    allocate(x_interp(nptb))
    allocate(zt_interp(nptb))
    allocate(zb_interp(nptb))
    allocate(thickness(nptb))
    allocate(camber(nptb))
    nptint = nptb
    side=.true.
    call spline_interp_z(curr_foil%npoint,curr_foil%x,curr_foil%z,xseedb,zt_interp(1:nptb),side)
    !call interp_vector(xseedt, zseedt, xseedb, zt_interp)
    x_interp = xseedb
    zb_interp = zseedb
  end if

! Compute thickness parameters

  tegap = zseedt(nptt) - zseedb(nptb)
  maxthick = 0.d0
  maxcamb = 0.d0
  heightfactor = tan(min_te_angle*acos(-1.d0)/180.d0)
  actual_min_te_angle = 2.d0

  do i = 2, nptint - 1

!   Thickness array and max thickness
    
    thickness(i) = zt_interp(i) - zb_interp(i)
    if (thickness(i) > maxthick) maxthick = thickness(i)
    
    camber(i) = ( zt_interp(i) + zb_interp(i) ) / 2.0d0
    if (camber(i) > maxcamb) maxcamb = camber(i)

!   Check if thinner than specified wedge angle on back half of airfoil
    if (x_interp(i) > te_angle_x_apply) then
      gapallow = tegap + heightfactor * (x_interp(nptint) - x_interp(i))
      TE_angle = atan((thickness(i) - tegap) / (x_interp(nptint) - x_interp(i)))
      if (TE_angle .LT. actual_min_TE_angle) actual_min_TE_angle = TE_angle
      if (thickness(i) < gapallow) then
        xtrans = x_interp(i)!/foilscale - xoffset
        write(text,'(F8.4)') xtrans
        text = adjustl(text)
        write(*,*) "   Detected too thin at x = "//trim(text)
        penaltyval = penaltyval + (gapallow - thickness(i))/0.1d0
      end if
    end if

  end do
  
  write(text,'(F8.4)') actual_min_te_angle/pi*180.d0
  text = adjustl(text)
  write(text2,'(F8.4)') te_angle_x_apply
  text2 = adjustl(text2)
  write(*,*) "   Min TE angle is "//trim(text)//" from "//trim(text2)//" to 1.0"

! Too thin on specified back

  if (penaltyval > 0.d0)                                                       &
     call ask_stop("Seed airfoil is thinner than min_te_angle near the "//&
                   "trailing edge.", allow_seed_penalties)
  penaltyval = 0.d0

! Check additional thickness constraints

  if (naddthickconst > 0) then
    call spline_interp_t(curr_foil%npoint,curr_foil%x,curr_foil%z,             &
                         addthick_x(1:naddthickconst),add_thickvec)
    !call interp_vector(x_interp, thickness,                                    &
    !                   addthick_x(1:naddthickconst), add_thickvec)

    addthick_violation = .false.
    do i = 1, naddthickconst
      write(text,'(F8.4)') addthick_x(i)
      text = adjustl(text)
      write(text2,'(F8.4)') add_thickvec(i)
      text2 = adjustl(text2)
      write(*,*) "   Thickness at x = "//trim(text)//": "//trim(text2)
      if ( (add_thickvec(i) < addthick_min(i)) .or.                            &
           (add_thickvec(i) > addthick_max(i)) ) then
        addthick_violation = .true.
      end if
    end do

    if (addthick_violation)                                                    &
      call ask_stop("Seed airfoil violates one or more thickness constraints.",&
        allow_seed_penalties)
  end if
  
! Free memory

  deallocate(x_interp)
  deallocate(zt_interp)
  deallocate(zb_interp)
  deallocate(thickness)
  deallocate(camber)

! Max thickness too low or too high

  write(text,'(F8.4)') maxthick
  text = adjustl(text)
  write(*,*) "   Max Thickness: "//trim(text)
  
  if (maxthick < min_thickness) then
    call ask_stop("Seed airfoil violates min_thickness constraint.", allow_seed_penalties)
  end if

  if (maxthick > max_thickness) then
    call ask_stop("Seed airfoil violates max_thickness constraint.", allow_seed_penalties)
  end if
  
! Camber too high or too low

  write(text,'(F8.4)') maxcamb
  text = adjustl(text)
  write(*,*) "   Max Camber: "//trim(text)
  
  if (maxcamb > max_camber) then
    call ask_stop("Seed airfoil violates max_camber constraint.", allow_seed_penalties)
  end if

  if (maxcamb < min_camber) then
    call ask_stop("Seed airfoil violates min_camber constraint.", allow_seed_penalties)
  end if
  
  ! Too blunt, sharp or disparate leading edge

  panang1 = atan((zseedt(2)-zseedt(1))/(xseedt(2)-xseedt(1))) *                &
            180.d0/acos(-1.d0)
  panang2 = atan((zseedb(1)-zseedb(2))/(xseedb(2)-xseedb(1))) *                &
            180.d0/acos(-1.d0)
  maxpanang = max(panang2,panang1)
  minpanang = min(panang2,panang1)
  write(text,'(F8.4)') maxpanang
  text = adjustl(text)
  write(*,*) "   Max LE panel angle: "//trim(text)//" degrees"
  if (maxpanang > max_leading_edge_angle) then
    call ask_stop("Seed airfoil's leading edge is too blunt.", allow_seed_penalties)
  end if
  write(text,'(F8.4)') minpanang
  text = adjustl(text)
  write(*,*) "   Min LE panel angle: "//trim(text)//" degrees"
  if (minpanang < min_leading_edge_angle) then
    call ask_stop("Seed airfoil's leading edge is too sharp.", allow_seed_penalties)
  end if
  write(text,'(F8.4)') abs(panang1 - panang2)
  text = adjustl(text)
  write(*,*) "   Dif LE panel angle: "//trim(text)//" degrees"
  if (abs(panang1 - panang2) > dif_leading_edge_angle) then
    call ask_stop("Seed airfoil's leading edge is too disparate.", allow_seed_penalties)
  end if
  
! Check for curvature reversals

  if (check_curvature) then

!   Compute curvature on top and bottom surfaces

    curvt = nu_curvature(nptt, xseedt, zseedt)
    curvb = nu_curvature(nptb, xseedb, zseedb)

!   Check number of reversals that exceed the threshold

    nreversalst = 0
    curv1 = 0.d0
    do i = 2, nptt - 1
      if (abs(curvt(i)) >= curv_threshold) then
        curv2 = curvt(i)
        if (curv2*curv1 < 0.d0) then
          xtrans = xseedt(i)!/foilscale - xoffset
          write(text,'(F8.4)') xtrans
          text = adjustl(text)
          ztrans = zseedt(i)!/foilscale - zoffset
          write(text2,'(F8.4)') ztrans
          text2 = adjustl(text2)
          write(*,*) "   Curvature reversal on top surface near (x, z) = ("//&
                         trim(text)//", "//trim(text2)//")"
          write(text,'(F8.4)') curvt(i)
          text = adjustl(text)
          write(*,*) "   Curvature: "//trim(text)
          nreversalst = nreversalst + 1
        end if
        curv1 = curv2
      end if
    end do

    write(text,'(I6)') nreversalst
    write(*,*) "   Number of curvature reversals on top surface: "// trim(text)
    if (nreversalst > max_curv_reverse_top)                                    &
      call ask_stop("Seed airfoil violates max_curv_reverse_top constraint.", allow_seed_penalties)

!   Bottom surface

    nreversalsb = 0
    curv1 = 0.d0
    do i = 2, nptb - 1
      if (abs(curvb(i)) >= curv_threshold) then
        curv2 = curvb(i)
        if (curv2*curv1 < 0.d0) then
          xtrans = xseedb(i)!/foilscale - xoffset
          write(text,'(F8.4)') xtrans
          text = adjustl(text)
          ztrans = zseedb(i)!/foilscale - zoffset
          write(text2,'(F8.4)') ztrans
          text2 = adjustl(text2)
          write(*,*) "   Curvature reversal on bot surface near (x, z) = ("//&
                         trim(text)//", "//trim(text2)//")"
          write(text,'(F8.4)') curvb(i)
          text = adjustl(text)
          write(*,*) "   Curvature: "//trim(text)
          nreversalsb = nreversalsb + 1
        end if
        curv1 = curv2
      end if
    end do
    write(text,'(I6)') nreversalsb
    write(*,*) "   Number of curvature reversals on bot surface: "// trim(text)
    if (nreversalsb > max_curv_reverse_bot)                                    &
      call ask_stop("Seed airfoil violates max_curv_reverse_bot constraint.", allow_seed_penalties)

  end if
  
  ! Penalty for too large panel angles
  
  call get_max_panel_angle(curr_foil, panel_angle)
  
  write(text,'(F8.4)') panel_angle
  text = adjustl(text)
  write(*,*) "   Max panel angle: "//trim(text)
  
  if (panel_angle > max_panel_angle) then
    call ask_stop("Seed airfoil panel angles are too large. Try adjusting "//&
                  "xfoil_paneling_options.", allow_seed_penalties)
  end if
  
! Get allowable panel growth rate

  growth_allowed = 0.d0

! Top surface growth rates
  max_growth_top = 0.0d0

  len1 = sqrt((xseedt(2)-xseedt(1))**2.d0 + (zseedt(2)-zseedt(1))**2.d0)
  do i = 2, nptt - 1
    len2 = sqrt((xseedt(i+1)-xseedt(i))**2.d0 + (zseedt(i+1)-zseedt(i))**2.d0)
    growth1 = len2/len1
    growth2 = len1/len2
    if (max(growth1,growth2) > max_growth_top)          &
        max_growth_top = max(growth1,growth2)
    if (max_growth_seed_mult * max(growth1,growth2) > growth_allowed)          &
        growth_allowed = max_growth_seed_mult*max(growth1,growth2)
    len1 = len2
  end do
  
  write(text,'(F8.4)') max_growth_top
  text = adjustl(text)
  write(*,*) "   Max growth rate at top: "//trim(text)

! Bottom surface growth rates

  len1 = sqrt((xseedb(2)-xseedb(1))**2.d0 + (zseedb(2)-zseedb(1))**2.d0)
  do i = 2, nptb - 1
    len2 = sqrt((xseedb(i+1)-xseedb(i))**2.d0 + (zseedb(i+1)-zseedb(i))**2.d0)
    growth1 = len2/len1
    growth2 = len1/len2
    if (max(growth1,growth2) > max_growth_bot)          &
        max_growth_bot = max(growth1,growth2)
    if (max_growth_seed_mult * max(growth1,growth2) > growth_allowed)          &
        growth_allowed = max_growth_seed_mult*max(growth1,growth2)
    len1 = len2
  end do
  
  write(text,'(F8.4)') max_growth_bot
  text = adjustl(text)
  write(*,*) "   Max growth rate at bot: "//trim(text)
  
  write(text,'(F8.4)') growth_allowed
  text = adjustl(text)
  write(*,*) "   Max growth rate allowed: "//trim(text)
  
! Check for bad combinations of operating conditions and optimization types

  do i = 1, noppoint
    write(text,*) i
    text = adjustl(text)

    opt_type = optimization_type(i)
    if ((op_point(i) <= 0.d0) .and. (op_mode(i) == 'spec-cl')) then
      if ( (trim(opt_type) /= 'min-drag') .and.                                &
           (trim(opt_type) /= 'max-xtr') .and.                                 &
           (trim(opt_type) /= 'max-lift-slope') ) then
        write(*,*) "Error: operating point "//trim(text)//" is at Cl = 0. "//  &
                 "Cannot use '"//trim(opt_type)//"' optimization in this case."
        write(*,*) 
        stop
      end if
    elseif ((op_mode(i) == 'spec-cl') .and.                                    &
            ((optimization_type(i) .eq. 'max-lift-search') .or.                &
            ((optimization_type(i) .eq. 'take_off') .and.                      &
            (optimization_correlation(i) .eq. 'init')))) then
      write(*,*) "Error: Cl is specified for operating point "//trim(text)//   &
                 ". Cannot use 'max-lift' optimization type in this case."
      write(*,*) 
      stop
    end if

  end do

  ! Get actual flap angles based on input variables

  actual_flap_degrees(1:noppoint) = flap_degrees(1:noppoint)
  
  ! Set identical flap angles
  do i = 1, nflap_identical
    flap_idi = flap_identical_points(i)
    actual_flap_degrees(flap_idi) = actual_flap_degrees(flap_identical_op(flap_idi))
  end do

  ! Set trasition points according to flap_connection 

  if (use_flap .and. (trim(flap_connection) .NE. 'smooth') .and.     &
      ( (trim(connection_apply) .EQ. 'trip_wire') .OR.               &
        (trim(connection_apply) .EQ. 'both'     ) )      ) then
    if (trim(flap_connection) .EQ. 'smooth-top') then
      xfoil_options%xtripb=x_flap
    elseif (trim(flap_connection) .EQ. 'smooth-bot') then
      xfoil_options%xtript=x_flap
    else
      xfoil_options%xtript=x_flap
      xfoil_options%xtripb=x_flap
    end if
  end if
  
! Analyze airfoil at requested operating conditions with Xfoil
  
  first_run_xfoil=.true.
  call run_xfoil(curr_foil, xfoil_geom_options, op_point(1:noppoint),          &
                 op_mode(1:noppoint), op_search, use_previous_op(1:noppoint),  &
                 reynolds(1:noppoint),                                         &
                 mach(1:noppoint), use_flap, x_flap, y_flap, y_flap_spec,      &
                 actual_flap_degrees(1:noppoint), xfoil_options, file_options, &
                 lift, drag, moment, viscrms, alpha, xtrt, xtrb, ncrit_pt)
  first_run_xfoil=.false.
  write(*,*)
    
! Check for unconverged points

  do i = 1, noppoint
    if (viscrms(i) > 1.0D-04) then
      write(text,*) i
      text = adjustl(text)
      call ask_stop("Xfoil calculations did not converge for operating "//&
                    "point "//trim(text)//".", .false.)
    end if
  end do

! Set aerodinamic constraint or check for violation of specified constraint

  do i = 1, noppoint
	! moment contraint
    if (trim(moment_constraint_type(i)) == 'use_seed') then
      write(text,'(F8.4)') moment(i)
      text = adjustl(text)
      write(text2,*) i
      text2 = adjustl(text2)
      write(*,*) "   Moment at op "//trim(text2)//" : "//trim(text)
      min_moment(i) = moment(i)
    elseif (trim(moment_constraint_type(i)) == 'specify') then
      write(text,'(F8.4)') moment(i)
      text = adjustl(text)
      write(text2,*) i
      text2 = adjustl(text2)
      write(*,*) "   Moment at op "//trim(text2)//": "//trim(text)
      if (moment(i) < min_moment(i)) then
        write(text,*) i
        text = adjustl(text)
        call ask_stop("Seed airfoil violates min_moment constraint for "//&
                      "operating point "//trim(text)//".", allow_seed_penalties)
      end if
    end if
	! lift contraint
    if (trim(lift_constraint_type(i)) == 'use_seed') then
      write(text,'(F8.4)') lift(i)
      text = adjustl(text)
      write(text2,*) i
      text2 = adjustl(text2)
      write(*,*) "   Lift at op "//trim(text2)//" : "//trim(text)
      min_lift(i) = lift(i)
    elseif (trim(lift_constraint_type(i)) == 'specify') then
      write(text,'(F8.4)') lift(i)
      text = adjustl(text)
      write(text2,*) i
      text2 = adjustl(text2)
      write(*,*) "   Lift at op "//trim(text2)//" : "//trim(text)
      if (lift(i) < min_lift(i)) then
        write(text,*) i
        text = adjustl(text)
        call ask_stop("Seed airfoil violates min_lift constraint for "//&
                      "operating point "//trim(text)//".", allow_seed_penalties)
      end if
    end if
	! drag contraint
    if (trim(drag_constraint_type(i)) == 'use_seed') then
      write(text,'(F8.4)') drag(i)
      text = adjustl(text)
      write(text2,*) i
      text2 = adjustl(text2)
      write(*,*) "   Drag at op "//trim(text2)//" : "//trim(text)
      max_drag(i) = drag(i)
    elseif (trim(drag_constraint_type(i)) == 'specify') then
      write(text,'(F8.4)') drag(i)
      text = adjustl(text)
      write(text2,*) i
      text2 = adjustl(text2)
      write(*,*) "   Drag at op "//trim(text2)//" : "//trim(text)
      if (drag(i) > max_drag(i)) then
        write(text,*) i
        text = adjustl(text)
        call ask_stop("Seed airfoil violates max_drag constraint for "//&
                      "operating point "//trim(text)//".", allow_seed_penalties)
      end if
    end if
  end do

! Evaluate objectives to establish scale factors for each point

  do i = 1, noppoint
    write(text,*) i
    text = adjustl(text)

    if (lift(i) <= 0.d0 .and. (trim(optimization_type(i)) == 'min-sink' .or.   &
        trim(optimization_type(i)) == 'max-glide') ) then
      write(*,*) "Error: operating point "//trim(text)//" has Cl <= 0. "//     &
                 "Cannot use "//trim(optimization_type(i))//" optimization "// &
                 "in this case."
      write(*,*)
      stop
    end if

    if (trim(optimization_type(i)) == 'min-sink') then
      checkval = drag(i)/lift(i)**1.5d0
    elseif (trim(optimization_type(i)) == 'max-glide') then
      checkval = drag(i)/lift(i)
    elseif (trim(optimization_type(i)) == 'min-drag') then
      checkval = drag(i)
    elseif (trim(optimization_type(i)) == 'max-lift') then
      checkval = 1.d0/lift(i)
    elseif (trim(optimization_type(i)) == 'max-xtr') then
      checkval = 1.d0/(0.5d0*(xtrt(i)+xtrb(i))+0.1d0)  ! Ensure no division by 0
    elseif (trim(optimization_type(i)) == 'max-lift-slope') then

!     Maximize dCl/dalpha (0.1 factor to ensure no division by 0)

      checkval = 0.d0
      if (i < noppoint) then
        if (alpha(i+1) > alpha(i)) then
          checkval = derv1f1(lift(i+1), lift(i),                               &
                             (alpha(i+1)-alpha(i)+0.1d0)*pi/180.d0)
        else
          checkval = derv1b1(lift(i+1), lift(i),                               &
                             (alpha(i)-alpha(i+1)+0.1d0)*pi/180.d0)
        end if
      end if

      if (i > 1) then
        if (alpha(i) > alpha(i-1)) then
          checkval = checkval + derv1b1(lift(i-1), lift(i),                    &
                                        (alpha(i)-alpha(i-1)+0.1d0)*pi/180.d0) 
        else
          checkval = checkval + derv1f1(lift(i-1), lift(i),                    &
                                        (alpha(i-1)-alpha(i)+0.1d0)*pi/180.d0) 
        end if
      end if
      if ( (i < noppoint) .and. (i > 1) ) checkval = checkval/2.d0 

!     4*pi factor is to move singularity location. Without it, negative
!     objective function occurs for Cl_a < 0. With this factor, it occurs for
!     Cl_a < -4*pi (hopefully never encountered).

      checkval = 1.d0/(checkval + 4.d0*pi)
    elseif (trim(optimization_type(i)) == 'target-lift') then
      
      checkval = (1.D-9 +               &
        (lift(i)-target_value(i))**2.d0 / (lift(i) + 1.D-9) )
      write(*,*) checkval
    elseif (trim(optimization_type(i)) == 'target-drag') then
      
      checkval = (1.D-9 +               &
        (drag(i)-target_value(i))**2.d0 / (drag(i) + 1.D-9) )
    elseif (trim(optimization_type(i)) == 'target-moment') then
      
      checkval = (1.D-9 +               &
        (moment(i)-target_value(i))**2.d0 /                  &
        (moment(i) + 1.D-9) )
    elseif (trim(optimization_type(i)) == 'target-xtrt') then
      
      checkval = (1.D-9 +               &
        (xtrt(i)-target_value(i))**2.d0 / (xtrt(i) + 1.D-9) )
    elseif (trim(optimization_type(i)) == 'target-xtrb') then
      
      checkval = (1.D-9 +               &
        (xtrb(i)-target_value(i))**2.d0 / (xtrb(i) + 1.D-9) )
    elseif (trim(optimization_type(i)) == 'target-glide') then
      
      checkval = (1.D-9 +               &
        (lift(i)/drag(i)-target_value(i))**2.d0 /            &
        (lift(i)/drag(i) + 1.D-9) )
    elseif (trim(optimization_type(i)) == 'target-sink') then
      
      checkval = (1.D-9 +               &
        (lift(i)**1.5d0/drag(i)-target_value(i))**2.d0 /     &
        (lift(i)**1.5d0/drag(i) + 1.D-9) )
    elseif (trim(optimization_type(i)) == 'max-lift-search') then

      checkval = 1.d0/lift(i)
    elseif (trim(optimization_type(i)) == 'take-off') then
       
      if(trim(optimization_correlation(i)) == 'init')then
        ninit = i
        checkval = 1.D-9
      end if
      if(trim(optimization_correlation(i)) == 'end') then
        nend = i
        checkval = 1.d0/1000 
      end if
    elseif (trim(optimization_type(i)) == 'climb') then
        
      checkval = 1.D-9 
      if(trim(optimization_correlation(i)) == 'init')then
        ninit = i
      end if
      if(trim(optimization_correlation(i)) == 'end') then
        nend = i
        checkval = 1.d0/1000 
      end if
    elseif (trim(optimization_type(i)) == 'dash') then
        
      checkval = 1.D-9 
      if(trim(optimization_correlation(i)) == 'init')then
        ninit = i
      end if
      if(trim(optimization_correlation(i)) == 'end') then
        nend = i
        checkval = 1.d0/1000 
      end if
    elseif (trim(optimization_type(i)) == 'turn') then
        
      checkval = 1.D-9 
      if(trim(optimization_correlation(i)) == 'init')then
        ninit1 = i
      end if
      if(trim(optimization_correlation(i)) == 'end') then
        nend1 = i 
      end if
    else
      write(*,*)
      write(*,*) "Error: requested optimization_type for operating point "//   &
                 trim(text)//" not recognized."
      stop
    end if
    scale_factor(i) = 1.d0/checkval
  end do
  
  write(*,*)
  write(*,*) ' Airfoil passed all constraints'
  write(*,*)
  
end subroutine check_seed

!=============================================================================80
!
! Asks user to stop or continue
!
!=============================================================================80
subroutine ask_stop(message, allow)

  character(*), intent(in) :: message
  logical, intent(in) :: allow

  character :: choice
  logical :: valid_choice

  if (allow) then
    write(*,'(A)') '    Warning: '//trim(message)
  else
    ! Get user input

    valid_choice = .false.
    do while (.not. valid_choice)
  
      write(*,'(A)') '    Warning: '//trim(message)
      write(*,'(A)', advance='no') '   Continue anyway? (y/n): '
      read(*,'(A)') choice

      if ( (choice == 'y') .or. (choice == 'Y') ) then
        valid_choice = .true.
        choice = 'y'
      else if ( ( choice == 'n') .or. (choice == 'N') ) then
        valid_choice = .true.
        choice = 'n'
      else
        write(*,'(A)') '    Please enter y or n.'
        valid_choice = .false.
      end if

    end do

  ! Stop or continue

    write(*,*)
    if (choice == 'n') stop
  end if
  write(*,*)


end subroutine ask_stop
  
end module input_sanity
