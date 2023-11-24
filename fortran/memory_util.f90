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

module memory_util

! Module containing utilities for memory management

  implicit none

  contains

!=============================================================================80
!
! Allocates memory for buffer airfoil
!
!=============================================================================80
subroutine allocate_airfoil(foil)

  use vardef, only : airfoil_type

  type(airfoil_type), intent(inout) :: foil

  integer :: npoint
 
  npoint = foil%npoint
  allocate(foil%x(npoint))
  allocate(foil%z(npoint))

end subroutine allocate_airfoil

!=============================================================================80
!
! Deallocates memory for buffer airfoil
!
!=============================================================================80
subroutine deallocate_airfoil(foil)

  use vardef, only : airfoil_type

  type(airfoil_type), intent(inout) :: foil

  deallocate(foil%x)
  deallocate(foil%z)

end subroutine deallocate_airfoil

!=============================================================================80
!
! Allocates memory for airfoil optimization
!
!=============================================================================80
subroutine allocate_airfoil_data()

  use xfoil_driver,       only : xfoil_init
  use vardef,             only : nparams_top, nparams_bot, shape_functions,    &
                                 xseedt, xseedb, curr_foil, objfunction_return,&
                                 contrain_number, noppoint
  use parametrization,    only : create_shape_functions, parametrization_dvs

  double precision, dimension(:), allocatable :: modest, modesb
  integer :: ndvs_top, ndvs_bot

! Allocate shape function setup arrays
  
  call parametrization_dvs(nparams_top, nparams_bot, shape_functions, ndvs_top, ndvs_bot)
    
  allocate(modest(ndvs_top))
  allocate(modesb(ndvs_bot))
  
  modest(:) = 0.d0
  modesb(:) = 0.d0

! Allocate private memory for airfoil optimization on each thread

!$omp parallel default(shared)

! For NACA, this will create the shape functions.  For Hicks-Henne,
! it will just allocate them.

  call create_shape_functions(xseedt, xseedb, modest, modesb, shape_functions)

! Allocate memory for working airfoil on each thread

  curr_foil%npoint = size(xseedt,1) + size(xseedb,1) - 1
  call allocate_airfoil(curr_foil)
  
! allocate objfunction_return
  allocate(objfunction_return%constrains_data(contrain_number))
  allocate(objfunction_return%aero_data(3*noppoint+1))

! Allocate memory for xfoil

  call xfoil_init()

!$omp end parallel

! Deallocate shape function setup arrays

  deallocate(modest)
  deallocate(modesb)

end subroutine allocate_airfoil_data

!=============================================================================80
!
! Frees memory used during airfoil optimization
!
!=============================================================================80
subroutine deallocate_airfoil_data()

  use parametrization,    only : deallocate_parametrization
  use vardef,             only : curr_foil, objfunction_return
  use xfoil_driver,       only : xfoil_cleanup

!$omp parallel default(shared)

  call deallocate_parametrization()
  call deallocate_airfoil(curr_foil)
  ! deallocate objfunction_return
  deallocate(objfunction_return%constrains_data)
  deallocate(objfunction_return%aero_data)
  call xfoil_cleanup()

!$omp end parallel

end subroutine deallocate_airfoil_data

end module memory_util
