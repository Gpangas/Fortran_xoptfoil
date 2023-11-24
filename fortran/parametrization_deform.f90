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
  
module parametrization_deform

! Contains subroutines to create an airfoil shape from design variables

  implicit none

! Shape functions for creating airfoil shapes (top and bottom)

  double precision, dimension(:,:), pointer :: top_shape_function
  double precision, dimension(:,:), pointer :: bot_shape_function

!$omp threadprivate(top_shape_function)
!$omp threadprivate(bot_shape_function)

  contains

!=============================================================================80
!
! Allocates memory for shape functions
!
!=============================================================================80
subroutine allocate_shape_functions(nmodest, nmodesb, npointst, npointsb)

  integer, intent(in) :: nmodest, nmodesb, npointst, npointsb

  allocate(top_shape_function(nmodest,npointst))
  allocate(bot_shape_function(nmodesb,npointsb))

  !   Initialize shape functions

  top_shape_function(:,:) = 0.d0
  bot_shape_function(:,:) = 0.d0
end subroutine allocate_shape_functions

!=============================================================================80
!
! Deallocates memory for shape functions
!
!=============================================================================80
subroutine deallocate_shape_functions

  deallocate(top_shape_function)
  deallocate(bot_shape_function)

end subroutine deallocate_shape_functions

!=============================================================================80
!
! Populates shape function arrays
! For Hicks-Hene shape functions, number of elements in modes must be a 
! multiple of 3.
!
!=============================================================================80
subroutine create_shape(x, modes, shapetype, shape_function)

  double precision, dimension(:), intent(in) :: x, modes
  character(*), intent(in) :: shapetype
  double precision, dimension(:,:), intent(inout) :: shape_function

  shape_switch: if (trim(shapetype) == 'naca') then
    
    call NACA_shape(x, modes, shape_function)
    
  elseif (trim(shapetype) == 'hicks-henne') then
    
    call HH_shape(x, modes, shape_function)
  
  else

    write(*,*)
    write(*,*) 'Shape function '//trim(shapetype)//' not recognized.'
    write(*,*)
    stop

  end if shape_switch

end subroutine create_shape

!=============================================================================80
!
! Populates shape function arrays for Hicks-Hene shape functions,
! number of elements in modes must be a multiple of 3.
!
!=============================================================================80
subroutine HH_shape(x, modes, shape_function)

  use vardef, only : abs_initial_perturb, min_bump_width

  double precision, dimension(:), intent(in) :: x, modes
  double precision, dimension(:,:), intent(inout) :: shape_function

  integer :: npt, nmodes, i, j, counter1
  double precision :: power1, st, t1, t2, t1fact, t2fact, pi
  double precision :: chord, xle, xs

  npt = size(x,1)
  chord = x(npt) - x(1)
  xle = x(1)
  
  nmodes = size(modes,1)/3
  t1fact = abs_initial_perturb/(1.d0 - 0.001d0)
  t2fact = abs_initial_perturb/(10.d0 - min_bump_width)
  pi = acos(-1.d0)

  do i = 1, nmodes

    !     Extract strength, bump location, and width

    counter1 = 3*(i-1)
    st = modes(counter1+1)

    t1 = modes(counter1+2)/t1fact
    t2 = modes(counter1+3)/t2fact

    !     Check for problems with bump location and width parameters

    if (t1 <= 0.d0) t1 = 0.001d0
    if (t1 >= 1.d0) t1 = 0.999d0
    if (t2 <= 0.d0) t2 = 0.001d0

    !     Create shape function

    power1 = log10(0.5d0)/log10(t1)
    do j = 2, npt-1
      xs = (x(j)-xle)/chord
      shape_function(i,j) = st*sin(pi*xs**power1)**t2
    end do

  end do
end subroutine HH_shape
!=============================================================================80
!
! Populates shape function arrays for NACA shape functions
!
!=============================================================================80
subroutine NACA_shape(x, modes, shape_function)

  double precision, dimension(:), intent(in) :: x, modes
  double precision, dimension(:,:), intent(inout) :: shape_function

  integer :: npt, nmodes, i, j, counter1, counter2
  double precision :: power1, power2, dvscale
  double precision :: chord, xle, xs

  npt = size(x,1)
  chord = x(npt) - x(1)
  xle = x(1)

  nmodes = size(modes,1)

!   Create naca shape functions

  do j = 1, npt
    xs = (x(j)-xle)/chord
    shape_function(1,j) = sqrt(xs) - xs
  end do

  counter1 = 1
  counter2 = 1

  do i = 2, nmodes

!     Whole-powered shapes

    if (counter2 == 1) then

      power1 = dble(counter1)
      do j = 1, npt
        xs = (x(j)-xle)/chord
        shape_function(i,j) = xs**(power1)*(1.d0 - xs)
      end do
      counter2 = 2

!     Fractional-powered shapes

    else

      power1 = 1.d0/dble(counter1 + 2)
      power2 = 1.d0/dble(counter1 + 1)
      do j = 1, npt
        xs = (x(j)-xle)/chord
        shape_function(i,j) = xs**power1 - xs**power2
      end do
      counter2 = 1
      counter1 = counter1 + 1
       
    end if

  end do

!   Normalize shape functions

  do i = 1, nmodes
    dvscale = 1.d0/abs(maxval(shape_function(i,:)))
    shape_function(i,:) = shape_function(i,:)*dvscale*modes(i)
  end do

end subroutine NACA_shape
!=============================================================================80
!
! Creates an airfoil surface by perturbing an input "seed" airfoil
! Using Hicks-Henne bump functions
!
!=============================================================================80
subroutine HH_airfoil(xt_seed, zt_seed, xb_seed, zb_seed, modest, modesb,      &
                          zt_new, zb_new, symmetrical)

  double precision, dimension(:), intent(in) :: xt_seed, zt_seed, xb_seed,     &
                                                zb_seed
  double precision, dimension(:), intent(in) :: modest, modesb
  double precision, dimension(:), intent(inout) :: zt_new, zb_new
  logical, intent(in) :: symmetrical

  integer :: i, nmodest, nmodesb, npointst, npointsb
  
  nmodest = size(modest,1)/3
  nmodesb = size(modesb,1)/3
  
  npointst = size(zt_seed,1)
  npointsb = size(zb_seed,1)

! Create shape functions for Hicks-Henne

 !   Create shape functions for top

  call create_shape(xt_seed, modest, 'hicks-henne', top_shape_function)

  !   Create shape functions for bottom

  call create_shape(xb_seed, modesb, 'hicks-henne', bot_shape_function)

! Top surface

  zt_new = zt_seed
  do i = 1, nmodest
    zt_new = zt_new + top_shape_function(i,1:npointst)
  end do

! Bottom surface

  if (.not. symmetrical) then
    zb_new = zb_seed
    do i = 1, nmodesb
      zb_new = zb_new + bot_shape_function(i,1:npointsb)
    end do

! For symmetrical airfoils, just mirror the top surface

  else
    zb_new = -zt_new
  end if

end subroutine HH_airfoil

!=============================================================================80
!
! Creates an airfoil surface by perturbing an input "seed" airfoil
! Using NACA bump functions
!
!=============================================================================80
  subroutine NACA_airfoil(xt_seed, zt_seed, xb_seed, zb_seed, modest, modesb,  &
                          zt_new, zb_new, symmetrical)

  double precision, dimension(:), intent(in) :: xt_seed, zt_seed, xb_seed,     &
                                                zb_seed
  double precision, dimension(:), intent(in) :: modest, modesb
  double precision, dimension(:), intent(inout) :: zt_new, zb_new
  logical, intent(in) :: symmetrical

  integer :: i, nmodest, nmodesb, npointst, npointsb
  double precision :: strength

  nmodest = size(modest,1)
  nmodesb = size(modesb,1)

  npointst = size(zt_seed,1)
  npointsb = size(zb_seed,1)

   !   Create shape functions for top

  call create_shape(xt_seed, modest, 'naca', top_shape_function)

  !   Create shape functions for bottom

  call create_shape(xb_seed, modesb, 'naca', bot_shape_function)
  
! Top surface

  zt_new = zt_seed
  do i = 1, nmodest
    strength = modest(i)
    zt_new = zt_new + strength*top_shape_function(i,1:npointst)
  end do

! Bottom surface

  if (.not. symmetrical) then
    zb_new = zb_seed
    do i = 1, nmodesb
      strength = modesb(i)
      zb_new = zb_new + strength*bot_shape_function(i,1:npointsb)
    end do

! For symmetrical airfoils, just mirror the top surface

  else
    zb_new = -zt_new
  end if
  
end subroutine NACA_airfoil
end module parametrization_deform