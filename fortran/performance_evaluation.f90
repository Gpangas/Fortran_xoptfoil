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

!  Copyright (C) 2017-2019 Daniel Prosser, 2020-2021 Ricardo Palmeira,
!  2023-2024 Guilherme Pangas
    
module performance_evaluation

! Module containing Aircraft Polar routine
    
  implicit none
    
    contains
    
!=============================================================================80
!
! Calculate Performance
!
!=============================================================================80         
subroutine calculate_performance(moment, drag, lift, alpha, viscrms, points)

  use vardef, only: noppoint, optimization_type, weight, climb, dash, turn,    &
                    weight_min, RC_min, dash_V_min, turn_V_min
  double precision, dimension(:), intent(in) :: moment, drag, lift, alpha,     &
                                                viscrms
  integer :: i, ninit_1, ninit_2, nend_1, nend_2
  double precision, dimension(:), intent(out) :: points
  
  do i = 1, noppoint
  
  !   Extra checks for really bad designs
    if (viscrms(i) >= 1.d0) then
      points(:) = 0.d0
      return
    end if  
  
    if (trim(optimization_type(i)) == 'take-off') then
        
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'take-off')) ninit_1 = i
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'take-off')) then
        nend_1 = i
        call take_off_evaluation(lift(ninit_1), lift(nend_1), drag(nend_1),    &
          points(1))
        if(weight_min > weight) return
      end if
    elseif (trim(optimization_type(i)) == 'climb') then
       
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'climb')) ninit_1 = i               
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'climb')) then
        nend_1 = i
        call climb_evaluation(ninit_1, nend_1, drag, lift, points(2)) 
        if(RC_min > climb%RC_max) return
        if(0 > climb%t_accel) return
      end if
    elseif (trim(optimization_type(i)) == 'dash') then
       
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'dash')) ninit_1 = i
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'dash')) then
        nend_1 = i
        call dash_evaluation(ninit_1, nend_1, ninit_2, nend_2, drag, lift,     &
          points(3))
        if(dash_V_min > dash%V_max) return
        if(0 > dash%t_accel) return
        if(0 > dash%dist_accel) return
        if(0 > dash%dist) return
        if(turn_V_min > turn%V) return
      end if
    elseif (trim(optimization_type(i)) == 'turn') then
      
      if((i == 1) .or. (trim(optimization_type(i-1)) /= 'turn')) ninit_2 = i
      if((i == noppoint) .or. (trim(optimization_type(i+1)) /= 'turn'))        &
        nend_2 = i
    end if
  end do
  
end subroutine calculate_performance
!=============================================================================80
!
! Take-off evaluation
!
!=============================================================================80  
subroutine take_off_evaluation(Cl_int, Cl_end, Cd_end, points)

  use aircraft_polar
  use vardef, only : weight, weight_i, S_w, take_off, climb

  double precision, intent(in) :: Cl_int, Cl_end, Cd_end
  double precision, intent(out) :: points
  
  double precision, parameter :: n=1, g = 9.80665
  double precision :: converge
  double precision :: rho
  double precision :: Cl_max, V_s
  double precision :: Cl_run, Cd_run, thurst
  double precision :: S_g_i
  double precision :: weight_payload
  
  Cl_max = 0.9 * Cl_int
  call rho_calculation(take_off%h, rho)
  weight = weight_i
  V_s = sqrt(abs(2*weight/(rho*Cl_max*S_w)))  
  take_off%V_to = take_off%A_1 * V_s
  take_off%V_run = take_off%V_to/sqrt(2.0_8)
  call aircraft_data_take_off(Cl_end, Cd_end, n, take_off%h, take_off%V_run,   &
      Cl_run, Cd_run, thurst)
  converge = 1
  do while (converge > 1E-3)
    V_s = sqrt(abs(2*weight/(rho*Cl_max*S_w)))  
    take_off%V_to = take_off%A_1 * V_s
    take_off%V_run = take_off%V_to/sqrt(2.0_8)

    S_g_i = ((take_off%A_1**2)*(weight/S_w))/(rho*g*((thurst/weight-           & 
      take_off%miu)*Cl_max+((take_off%A_1**2)/2)*(take_off%miu*Cl_run-Cd_run)))
    
    converge = abs((take_off%S_g-S_g_i) / take_off%S_g)
    weight = weight*(1 + 0.5*(take_off%S_g-S_g_i) / take_off%S_g)
  end do 
 
  climb%V_0 = take_off%V_to
  weight_payload = weight - take_off%weight_empty
  points = weight_payload/take_off%weight_payload_ref * 1000 

end subroutine take_off_evaluation
!=============================================================================80
!
! Climb evaluation
!
!=============================================================================80
subroutine climb_evaluation(oppoint_init, oppoint_end, drag, lift, points)

  use aircraft_polar
  use vardef, only : weight, climb, dash
  
  integer, intent(in) :: oppoint_init, oppoint_end
  double precision, dimension(:), intent(in) :: drag, lift
  double precision, intent(out) :: points
  
  integer :: i, n_oppoint_climb, i_RC_max
  double precision, parameter :: n = 1, g = 9.80665
  double precision, allocatable :: T(:), D(:), V(:), RC(:)
  double precision :: h, t_climb
  n_oppoint_climb = oppoint_end - oppoint_init + 1
  
  allocate(T(n_oppoint_climb), D(n_oppoint_climb), V(n_oppoint_climb),           &
           RC(n_oppoint_climb))

  do i=1, n_oppoint_climb
    call aircraft_data_eq(lift(oppoint_init+i-1), drag(oppoint_init+i-1), n,     &
                           climb%h, V(i), D(i), T(i))
     RC(i) = (T(i)*V(i)-D(i)*V(i))/weight
     if(i .EQ. 1)then
       climb%RC_max = RC(1)
       i_RC_max = 1
     elseif(climb%RC_max < RC(i))then
       climb%RC_max = RC(i)
       i_RC_max = i
     end if   
  end do
  
  climb%V = V(i_RC_max)
  climb%Cl = lift(i_RC_max)
  
  climb%t_accel = 0.d0
  if(climb%accel)then
    climb%t_accel = (V(1) - climb%V_0)/((T(1)-D(1))*(g/weight))
    if(i_RC_max .NE. 1)then
      accel_to_climb: do i=2, n_oppoint_climb
        climb%t_accel = climb%t_accel + (V(i)-V(i-1)) / (((T(i)-D(i))+(T(i-1)-D(i-1)))/2*    &
                  (g/weight))
        if(climb%t_accel .GT. climb%time)then
          points = 0
          return
        end if
        if(i_RC_max .EQ. i) exit accel_to_climb
      end do accel_to_climb
    end if
  end if
    
  t_climb = climb%t_accel + climb%dh/climb%RC_max
  
  if(t_climb .gt. climb%time)then
    h = climb%dh - (t_climb - climb%time)/climb%RC_max
    points = climb%points_coeff(1)*h**4 + climb%points_coeff(2)*h**3 +         &  
              climb%points_coeff(3)*h**2 + climb%points_coeff(4)*h +           &
              climb%points_coeff(5)
    points = points/1.203
  else
    points = 1000
  end if
  
  if(t_climb .lt. climb%time) dash%t_ex = climb%time- t_climb
  
  dash%V_0 = V(i_RC_max)
  
  deallocate(T, D, V, RC)
  
end subroutine climb_evaluation
!=============================================================================80
!
! Dash  evaluation
!
!=============================================================================80
subroutine dash_evaluation(oppoint_init_d, oppoint_end_d, oppoint_init_t,      &
                           oppoint_end_t, drag, lift, points)

  use aircraft_polar
  use vardef, only : weight, dash, turn
  
  integer, intent(in) :: oppoint_init_d, oppoint_end_d
  integer, intent(in) :: oppoint_init_t, oppoint_end_t 
  double precision, dimension(:), intent(in) :: drag, lift
  double precision, intent(out) :: points
  
  integer :: i, j, n_oppoint_dash, n_oppoint_turn
  double precision, parameter :: n = 1, g = 9.80665, pi = 3.1415926
  double precision, allocatable :: T_d(:), D_d(:), V_d(:)
  double precision, allocatable :: T_t(:), D_t(:), V_t(:)
  double precision :: t_accel_i, turn_perimeter
  
  n_oppoint_dash = oppoint_end_d - oppoint_init_d + 1
  n_oppoint_turn = oppoint_end_t - oppoint_init_t + 1
  
  allocate(T_d(n_oppoint_dash), D_d(n_oppoint_dash), V_d(n_oppoint_dash))
  allocate(T_t(n_oppoint_turn), D_t(n_oppoint_turn), V_t(n_oppoint_turn))
  
  !Calculate Velocity, Drag and Trusth for every dash operating point
  do i=1, n_oppoint_dash
    call aircraft_data_eq(lift(oppoint_init_d+i-1), drag(oppoint_init_d+i-1),  &
                            n, dash%h, V_d(i), D_d(i), T_d(i))
  end do
  
  !Calculate Max velocity and acelleration time and distance
  dash%t_accel = 0.d0
  dash%dist_accel = 0.d0
  if((T_d(1)-D_d(1)) .le. 0)then
    points = 0.d0
    return
  !criar msg de erro para Cl definido demasiado baixo
  else
    dash%V_max = V_d(1)
    dash%t_accel = (V_d(1) - dash%V_0)/((T_d(1)-D_d(1))*(g/weight))
    if(dash%t_accel .gt. dash%t_ex) dash%dist_accel = (V_d(1) + dash%V_0)/2 *  &
      dash%t_accel
  end if
  
  accel_to_dash: do i = 2, n_oppoint_dash
    if((T_d(i)-D_d(i)) .gt. 0)then
      dash%V_max = V_d(i)
      t_accel_i = (V_d(i) - V_d(i-1))/(((T_d(i)-D_d(i))+(T_d(i-1)-D_d(i-1)))/2*&
        (g/weight))
      dash%t_accel = dash%t_accel + t_accel_i
      if(dash%t_accel .gt. dash%t_ex)then
        if(dash%t_ex .gt. (dash%t_accel-t_accel_i))then
          dash%dist_accel = dash%dist_accel + (V_d(i) + V_d(i-1))/2 *          & 
            (dash%t_accel-dash%t_ex)
        else
          dash%dist_accel = dash%dist_accel + (V_d(i) + V_d(i-1))/2 * t_accel_i  
        end if
      end if
    else
      dash%V_max = V_d(i)-(T_d(i)-D_d(i))/(((T_d(i)-D_d(i))-(T_d(i-1)-D_d(i-1))&
        )/(V_d(i)-V_d(i-1)))
      t_accel_i = (dash%V_max-V_d(i-1))/((T_d(i-1)-D_d(i-1))/2*(g/weight))
      dash%t_accel = dash%t_accel + t_accel_i
      if(dash%t_accel .gt. dash%t_ex)then
        if(dash%t_ex .gt. (dash%t_accel-t_accel_i))then
          dash%dist_accel = dash%dist_accel + (dash%V_max + V_d(i-1))/2 *      & 
            (dash%t_accel-dash%t_ex)
        else
          dash%dist_accel = dash%dist_accel + (dash%V_max + V_d(i-1))/2 *      &
            t_accel_i  
        end if
      end if
      exit accel_to_dash
    end if
  end do accel_to_dash
  
  !Turn ponderation
  if(turn%activation)then
    do i=1, n_oppoint_turn
      call aircraft_data_eq(lift(oppoint_init_t+i-1), drag(oppoint_init_t+i-1),&
                            turn%n, turn%h, V_t(i), D_t(i), T_t(i))  
    end do
    
    !Check if Turn velocity is in between Cl range 
    if((T_t(1)-D_t(1)) .lt. 0.d0)then
      points = 0.d0
      return 
    end if 
    
    !Calculate Turn Velocity
    V_turn_calc: do i = 1, n_oppoint_turn
      if((T_t(i)-D_t(i)) .gt. 0.d0)then
        turn%V = V_t(i)
      else
        turn%V = V_t(i)-(T_t(i)-D_t(i))/(((T_t(i)-D_t(i))-(T_t(i+1)-D_t(i+1)))/&
                   (V_t(i)-V_t(i+1)))
        exit V_turn_calc
      end if
    end do V_turn_calc

    turn%radius = turn%V**2/(g*SQRT(turn%n**2-1))
    turn_perimeter = 2*pi*turn%radius
    if(dash%accel)then
      if(dash%t_accel .le. dash%t_ex)then
        dash%dist = (turn%dash_leg + turn_perimeter)/(turn%dash_leg/dash%V_max+&
               turn_perimeter/turn%V) * (dash%time) 
      else
        dash%dist = (turn%dash_leg + turn_perimeter)/(turn%dash_leg/dash%V_max+&
              turn_perimeter/turn%V) * (dash%time - (dash%t_accel-dash%t_ex))  &
              + dash%dist_accel  
      end if
    else
      dash%dist = (turn%dash_leg + turn_perimeter)/(turn%dash_leg/dash%V_max + &
              turn_perimeter/turn%V) * dash%time
    end if
  else
    if(dash%accel)then
      if(dash%t_accel .le. dash%t_ex)then
        dash%dist = dash%time*dash%V_max 
      else
        dash%dist = dash%dist_accel + dash%V_max*(dash%time-(dash%t_accel-     &
          dash%t_ex))
      end if  
    else
      dash%dist = dash%time*dash%V_max
    end if
  end if
  
  points = dash%dist/dash%dist_ref*1000
  
  if(points .lt. 0.d0)then
     points = 0.d0  
  end if
  
  deallocate(T_d, D_d, V_d, T_t, D_t, V_t)

  
end subroutine dash_evaluation
                             
                        
!=============================================================================80
!
! Caculate air density for a certain altitude, (ISA)
!
!=============================================================================80
subroutine rho_calculation(h, rho)
  
  double precision, intent(in) :: h
  double precision, intent(out) :: rho
  
  double precision :: T, p
  double precision, parameter :: T_0 = 228.15, lambda = -6.5E-3
  double precision, parameter :: p_0 = 101325
  double precision, parameter :: rho_0 = 1.225
  double precision, parameter :: g_0 = 9.80665, R = 287.05307
  
  T = T_0 + lambda*h
  p = p_0 * (T/T_0)**(-g_0/(lambda*R))
  rho = rho_0 * (T/T_0)**(-g_0/(lambda*R)-1)

end subroutine rho_calculation

end module performance_evaluation   