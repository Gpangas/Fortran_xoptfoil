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
  double precision :: V_to, V_run
  double precision :: Cl_run, Cd_run, thurst
  double precision :: S_g_i
  double precision :: weight_payload
  
  Cl_max = 0.9 * Cl_int
  call rho_calculation(take_off%h, rho)
  weight = weight_i
  V_s = sqrt(abs(2*weight/(rho*Cl_max*S_w)))  
  V_to = take_off%A_1 * V_s
  V_run = V_to/sqrt(2.0_8)
  call aircraft_data_take_off(Cl_end, Cd_end, n, take_off%h, V_run, Cl_run,  &
                                 Cd_run, thurst)
  converge = 1
  do while (converge > 0.1)
    V_s = sqrt(abs(2*weight/(rho*Cl_max*S_w)))  
    V_to = take_off%A_1 * V_s
    V_run = V_to/sqrt(2.0_8)

    S_g_i = ((take_off%A_1**2)*(weight/S_w))/(rho*g*((thurst/weight-take_off%miu)* &
             Cl_max+((take_off%A_1**2)/2)*(take_off%miu*Cl_run-Cd_run)))
    
    converge = abs((take_off%S_g-S_g_i) / take_off%S_g)
    weight = weight*(1 + 0.5*(take_off%S_g-S_g_i) / take_off%S_g)
  end do 
 
  climb%V_0 = V_to
  weight_payload = weight - take_off%weight_empty
  points = weight_payload/take_off%weight_payload_ref * 1000
  
  take_off%CL_max = Cl_max
  take_off%CL_run = Cl_run
  take_off%CD_run = Cd_run
  take_off%V_to = V_to
  take_off%V_run = V_run
  take_off%points = points  

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
  double precision :: RC_max, h
  double precision :: t_acel, t_climb
  n_oppoint_climb = oppoint_end - oppoint_init + 1
  
  allocate(T(n_oppoint_climb), D(n_oppoint_climb), V(n_oppoint_climb),           &
           RC(n_oppoint_climb))

  do i=1, n_oppoint_climb
    call aircraft_data_eq(lift(oppoint_init+i-1), drag(oppoint_init+i-1), n,     &
                           climb%h, V(i), D(i), T(i))
     RC(i) = (T(i)*V(i)-D(i)*V(i))/weight
     if(i .EQ. 1)then
       RC_max = RC(1)
       i_RC_max = 1
     elseif(RC_max < RC(i))then
       RC_max = RC(i)
       i_RC_max = i
     end if   
  end do
  
  t_acel = 0.d0
  if(climb%acel)then
    t_acel = (V(1) - climb%V_0)/((T(1)-D(1))*(g/weight))
    if(i_RC_max .NE. 1)then
      acel_to_climb: do i=2, n_oppoint_climb
        t_acel = t_acel + (V(i)-V(i-1)) / (((T(i)-D(i))+(T(i-1)-D(i-1)))/2*    &
                  (g/weight))
        if(t_acel .GT. climb%time)then
          points = 0
          return
        end if
        if(i_RC_max .EQ. i) exit acel_to_climb
      end do acel_to_climb
    end if
  end if
    
  t_climb = t_acel + climb%dh/RC_max
  
  if(t_climb .gt. climb%time)then
    h = climb%dh - (t_climb - climb%time)/RC_max
    points = climb%points_coeff(1)*h**4 + climb%points_coeff(2)*h**3 +         &  
              climb%points_coeff(3)*h**2 + climb%points_coeff(4)*h +           &
              climb%points_coeff(5)
    points = points/1.203
  else
    points = 1000
  end if
  
  if(t_climb .lt. climb%time) dash%t_ex = climb%time- t_climb
  
  dash%V_0 = V(i_RC_max)
  
  climb%RC_max = RC_max
  climb%V = V(i_RC_max)
  climb%Cl = lift(i_RC_max)
  climb%D = D(i_RC_max)
  climb%T = T(i_RC_max)
  climb%t_acel = t_acel
  climb%points = points
  
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
  double precision :: V_max, t_acel_d, t_acel_i, dist, dist_acel_d
  double precision :: V_turn, turn_radius, turn_perimeter
  
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
  t_acel_d = 0.d0
  dist_acel_d = 0.d0
  if((T_d(1)-D_d(1)) .le. 0)then
    points = 0.d0
    return
  !criar msg de erro para Cl definido demasiado baixo
  else
    V_max = V_d(1)
    t_acel_d = (V_d(1) - dash%V_0)/((T_d(1)-D_d(1))*(g/weight))
    if(t_acel_d .gt. dash%t_ex) dist_acel_d = (V_d(1) + dash%V_0)/2 * t_acel_d
  end if
  
  acel_to_dash: do i = 2, n_oppoint_dash
    if((T_d(i)-D_d(i)) .gt. 0)then
      V_max = V_d(i)
      t_acel_i = (V_d(i) - V_d(i-1))/(((T_d(i)-D_d(i))+(T_d(i-1)-D_d(i-1)))/2*(g/weight))
      t_acel_d = t_acel_d + t_acel_i
      if(t_acel_d .gt. dash%t_ex)then
        if(dash%t_ex .gt. (t_acel_d-t_acel_i))then
          dist_acel_d = dist_acel_d + (V_d(i) + V_d(i-1))/2 * (t_acel_d-       &
                        dash%t_ex)
        else
          dist_acel_d = dist_acel_d + (V_d(i) + V_d(i-1))/2 * t_acel_i  
        end if
      end if
    else
      V_max = V_d(i)-(T_d(i)-D_d(i))/(((T_d(i)-D_d(i))-(T_d(i-1)-D_d(i-1)))/   &
              (V_d(i)-V_d(i-1)))
      t_acel_i = (V_max-V_d(i-1))/((T_d(i-1)-D_d(i-1))/2*(g/weight))
      t_acel_d = t_acel_d + t_acel_i
      if(t_acel_d .gt. dash%t_ex)then
        if(dash%t_ex .gt. (t_acel_d-t_acel_i))then
          dist_acel_d = dist_acel_d + (V_max + V_d(i-1))/2 *(t_acel_d-dash%t_ex)
        else
          dist_acel_d = dist_acel_d + (V_max + V_d(i-1))/2 *t_acel_i  
        end if
      end if
      exit acel_to_dash
    end if
  end do acel_to_dash
  
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
        V_turn = V_t(i)
      else
        V_turn = V_t(i)-(T_t(i)-D_t(i))/(((T_t(i)-D_t(i))-(T_t(i+1)-D_t(i+1)))/&
                   (V_t(i)-V_t(i+1)))
        exit V_turn_calc
      end if
    end do V_turn_calc

    turn_radius = V_turn**2/(g*SQRT(turn%n**2-1))
    turn_perimeter = 2*pi*turn_radius
    if(dash%acel)then
      if(t_acel_d .le. dash%t_ex)then
        dist = (turn%dash_leg + turn_perimeter)/(turn%dash_leg/V_max +         &
               turn_perimeter/V_turn) * (dash%time) 
      else
        dist = (turn%dash_leg + turn_perimeter)/(turn%dash_leg/V_max +         &
              turn_perimeter/V_turn) * (dash%time - (t_acel_d-dash%t_ex))      &
              + dist_acel_d  
      end if
    else
      dist = (turn%dash_leg + turn_perimeter)/(turn%dash_leg/V_max +           &
              turn_perimeter/V_turn) * dash%time
    end if
  else
    if(dash%acel)then
      if(t_acel_d .le. dash%t_ex)then
        dist = dash%time*V_max 
      else
        dist = dist_acel_d + V_max*(dash%time - (t_acel_d-dash%t_ex))
      end if  
    else
      dist = dash%time*V_max
    end if
  end if
  
  points = dist/dash%dist_ref*1000
  
  if(points .lt. 0.d0)then
     points = 0.d0  
  end if
  
  deallocate(T_d, D_d, V_d, T_t, D_t, V_t)
  
  dash%V_max = V_max
  dash%t_acel_d = t_acel_d
  dash%dist_acel_d = dist_acel_d
  turn%V = V_turn
  turn%radius = turn_radius
  dash%dist = dist
  dash%points = points
  
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