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
    
module aircraft_polar

! Module containing Aircraft Polar routine
    
  implicit none
    
  contains
    
!=============================================================================80
!
! Caculate Aircraft's Drag, Velocity and Thurst
!
!=============================================================================80
subroutine aircraft_data_eq(Cl, Cd, n, h, V, drag_t, thurst)
    
  use vardef, only : weight, A_w, e_w, S_w, S_expose,fuselage, Cd_ld,          &
                        add_drag, thrust_coeff
  
  double precision, intent(in) :: Cl, Cd, n, h
  double precision, intent(out) :: V, drag_t, thurst
  
  double precision, parameter :: pi = 3.1415926536
  double precision :: Cl_a, Cl_i, alpha_i, converge
  double precision :: Cd_t, Cd_i, Cd_fus, Cd_tail
  double precision :: fuselage_d
  double precision :: rho, miu
  
  Cl_a = Cl
  fuselage_d = 1.3*(fuselage%height*fuselage%width)**0.625                      &
                /(fuselage%height+fuselage%width)**0.25
  
  converge = 1 
  do while (converge > 1E-5)
    alpha_i =  Cl_a/(pi*A_w*e_w)
    CL_i = Cl * cos(alpha_i)*(S_expose/S_w)
    converge = abs(Cl_a - Cl_i)
    Cl_a = Cl_a - 0.5*(Cl_a - Cl_i)
  end do
  
  Cd_i = (Cl_a**2)/(pi*A_w*e_w)
  
  call air_data(h, rho, miu)
  
  V = Sqrt((2*n*weight)/(rho*S_w*Cl_a))
  
  call fuselage_drag(V, rho, miu, Cd_fus)
  
  call tail_drag(V, rho, miu, Cd_tail)
  
  Cd_t = Cd + Cd_i + Cd_fus + Cd_tail + Cd_ld + add_drag
  
  drag_t = 0.5*rho*(V**2)*S_w*Cd_t
  
  thurst = thrust_coeff(1)*V**2 + thrust_coeff(2)*V + thrust_coeff(3)
  
end subroutine aircraft_data_eq
!=============================================================================80
!
! Caculate Aircraft's Lift, Drag and Thurst from define velocity
!
!=============================================================================80
subroutine aircraft_data_take_off(Cl, Cd, n, h, V, Cl_a, Cd_t, thurst)
    
  use vardef, only : A_w, e_w, S_w, S_expose,fuselage, Cd_ld, add_drag,        &
                        thrust_coeff
  
  double precision, intent(in) :: Cl, Cd, n, h, V
  double precision, intent(out) :: Cl_a, Cd_t, thurst
  
  double precision, parameter :: pi = 3.1415926536
  double precision :: Cl_i, alpha_i, converge
  double precision :: Cd_i, Cd_fus, Cd_tail
  double precision :: fuselage_d
  double precision :: rho, miu
  
  Cl_a = Cl
  fuselage_d = 1.3*(fuselage%height*fuselage%width)**0.625                     &
                /(fuselage%height+fuselage%width)**0.25
  
  converge = 1 
  do while (converge > 1E-5)
    alpha_i =  Cl_a/(pi*A_w*e_w)
    CL_i = Cl * cos(alpha_i)*(S_expose/S_w)
    converge = abs(Cl_a - Cl_i)
    Cl_a = Cl_a - 0.5*(Cl_a - Cl_i)
  end do
  
  Cd_i = (Cl_a**2)/(pi*A_w*e_w)
  
  call air_data(h, rho, miu)
  
  call fuselage_drag(V, rho, miu, Cd_fus)
  
  call tail_drag(V, rho, miu, Cd_tail)
  
  Cd_t = Cd + Cd_i + Cd_fus + Cd_tail + Cd_ld + add_drag
  
  
  thurst = thrust_coeff(1)*V**2 + thrust_coeff(2)*V + thrust_coeff(3)
  
end subroutine aircraft_data_take_off
!=============================================================================80
!
! Caculate Fuselage's Drag
!
!=============================================================================80
subroutine fuselage_drag(V, rho, miu, Cd_fus)

  use vardef, only : weight, A_w, e_w, S_w, fuselage 

  double precision, intent(in) :: V, rho, miu
  double precision, intent(out) :: Cd_fus

  double precision, parameter :: pi = 3.1415926536
  double precision :: re_i, re_cutoff, re_fus
  double precision :: fuselage_d, f, ff, C_f

  f = fuselage%length/sqrt(abs((4/pi)*fuselage%height*fuselage%width))
  
  ff = 1 + 60/f**3 + f/100
  
  re_i = V*rho*fuselage%length/miu
  
  re_cutoff = 38.21*(fuselage%length/fuselage%skin_roughness)**1.053
  
  if(re_i .GT. re_cutoff) then
    re_fus = re_cutoff
  else
    re_fus = re_i
  end if
  
  C_f = 0.455/log10(re_fus)**2.58
  
  Cd_fus =  C_f*ff*fuselage%wetted_area*fuselage%interference_factor/S_w

end subroutine fuselage_drag
!=============================================================================80
!
! Caculate Tail's Drag
!
!=============================================================================80
subroutine tail_drag(V, rho, miu, Cd_tail)

  use vardef, only : weight, A_w, e_w, S_w, tail 

  double precision, intent(in) :: V, rho, miu
  double precision, intent(out) :: Cd_tail
  
  integer :: i
  double precision, parameter :: pi = 3.1415926536
  double precision :: re_i, re_cutoff, re_tail
  double precision :: Q, ff, C_f, S_wetted
  
  if(tail%config .EQ. 1) then
    Q = 1.03
  else
    Q = 1.04
  end if
  
  Cd_tail = 0.d0
  
  do i = 1, tail%config
      
    re_i = V*rho*tail%chord(i)/miu
    re_cutoff = 38.21*(tail%chord(i)/tail%skin_roughness(i))**1.053
    
    if(re_i .GT. re_cutoff) then
      re_tail = re_cutoff
    else
      re_tail = re_i
    end if
    
    ff = 1 + 0.6/tail%max_t_x(i)*tail%t_c_ratio(i) + 100*tail%t_c_ratio(i)**4
    
    C_f = 0.5*(0.455/log10(re_tail)**2.58 + 1.328/sqrt(abs(re_tail)))
    
    S_wetted = (1.977 + 0.52*tail%t_c_ratio(i)) * tail%surface_area(i)
    
    if(tail%config .EQ. 1) then
      S_wetted = 2*S_wetted
    end if
  
    Cd_tail = Cd_tail + C_f*ff*S_wetted*Q/S_w
    
  end do

end subroutine tail_drag
!=============================================================================80
!
! Caculate air density and viscosity for a certain altitude, (ISA)
!
!=============================================================================80
subroutine air_data(h, rho, miu)
  
  double precision, intent(in) :: h
  double precision, intent(out) :: rho, miu
  
  double precision :: T, p
  double precision, parameter :: T_0 = 288.15, lambda = -6.5E-3
  double precision, parameter :: p_0 = 101325
  double precision, parameter :: rho_0 = 1.225
  double precision, parameter :: miu_ref = 1.78E-5
  double precision, parameter :: g_0 = 9.80665, R = 287.05307
  
  T = T_0 + lambda*h
  p = p_0 * (T/T_0)**(-g_0/(lambda*R))
  rho = rho_0 * (T/T_0)**(-g_0/(lambda*R)-1)
  miu = miu_ref * (T/T_0)**(3/4)

end subroutine air_data

end module aircraft_polar