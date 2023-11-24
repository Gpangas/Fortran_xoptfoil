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

program xfoil_only

! Runs xfoil by for an airfoil, but still uses the same input file as xoptfoil
! Doesn't do any transformations or scaling on the input airfoil

  use vardef
  use input_output,       only : read_inputs, read_clo
  use particle_swarm,     only : pso_options_type
  use genetic_algorithm,  only : ga_options_type
  use simplex_search,     only : ds_options_type
  use airfoil_evaluation, only : xfoil_options, xfoil_geom_options, file_options
  use memory_util,        only : deallocate_airfoil
  use airfoil_operations, only : load_airfoil
  use naca,               only : naca_options_type, naca_456
  use xfoil_driver,       only : run_xfoil, xfoil_geometry_info, xfoil_init,   &
                                 xfoil_cleanup

  implicit none

  type(airfoil_type) :: foil
  character(80) :: search_type, global_search, local_search, seed_airfoil,     &
                   matchfoil_file
  character(80) :: input_file
  type(naca_options_type) :: naca_options
  type(pso_options_type) :: pso_options
  type(ga_options_type) :: ga_options
  type(ds_options_type) :: ds_options
  integer :: restart_write_freq, i
  logical :: restart
  double precision, dimension(:), allocatable :: alpha, lift, drag, moment,    &
                                                 viscrms, xtrt, xtrb
  double precision :: maxt, xmaxt, maxc, xmaxc
  character(100) :: text, text1

! Set default names and read command line arguments

  input_file = 'inputs.txt'
  output_prefix = 'optfoil'
  call read_clo(input_file, output_prefix, "xfoil_only")

! Read inputs from namelist file

  call read_inputs(input_file, search_type, global_search, local_search,       &
                   seed_airfoil, nparams_top, nparams_bot,       &
                   restart, restart_write_freq, naca_options, &
                   pso_options, ga_options, ds_options, matchfoil_file)

! Allocate some things

  allocate(alpha(noppoint))
  allocate(lift(noppoint))
  allocate(drag(noppoint))
  allocate(moment(noppoint))
  allocate(viscrms(noppoint))
  allocate(xtrt(noppoint))
  allocate(xtrb(noppoint))

! Get airfoil to analyze, but don't do any transformations

  if (trim(seed_airfoil) == "from_file") then
    call load_airfoil(airfoil_file, foil)
  else if (trim(seed_airfoil) == "naca") then
    call naca_456(naca_options, 200, foil)
  end if

! Allocate xfoil variables

  call xfoil_init()

! Run xfoil

  file_options%design_number = 0
  file_options%polar = .true.
  file_options%cp = write_cp_file
  file_options%bl = write_bl_file
  file_options%write_all_airfoils = .true.
  
  call run_xfoil(foil, xfoil_geom_options, op_point(1:noppoint),               &
                 op_mode(1:noppoint), op_search, use_previous_op(1:noppoint),  &
                 reynolds(1:noppoint),                                         &
                 mach(1:noppoint), use_flap, x_flap, y_flap, y_flap_spec,      &
                 flap_degrees(1:noppoint), xfoil_options, file_options, lift,  &
                 drag, moment, viscrms, alpha, xtrt, xtrb, ncrit_pt(1:noppoint))

! Get geometry info from xfoil

  call xfoil_geometry_info(maxt, xmaxt, maxc, xmaxc)

! Display a summary of geometry results

  open(unit=122, file=trim(output_prefix)//'_xfoil_only.dat', status='replace')
  open(unit=321, file=trim(output_prefix)//'_xfoil_only.log', status='replace')
  
  write(*,*)
  write(*,*) 'Airfoil geometry information from Xfoil: '
  write(*,*)
  write(*,*) 'Max thickness: ', maxt
  write(*,*) 'Location of max thickness: ', xmaxt
  write(*,*) 'Max camber: ', maxc
  write(*,*) 'Location of max camber: ', xmaxc
  
  write(122,*)
  write(122,*) 'Airfoil geometry information from Xfoil: '
  write(122,*)
  write(122,*) 'Max thickness: ', maxt
  write(122,*) 'Location of max thickness: ', xmaxt
  write(122,*) 'Max camber: ', maxc
  write(122,*) 'Location of max camber: ', xmaxc
  
  write(321,'(A)') '  '
  write(321,'(A)') '       XFOIL         Version 6.99'
  write(321,'(A)') '  '
  write(text1,*) airfoil_file
  text1 = adjustl(text1)
  i=len_trim(text1)
  write(321,'(A)') ' Calculated polar for: '//text1(1:i-4)
  write(321,'(A)') '  '
  write(321,'(A)') ' 1 1 Reynolds number fixed          Mach number fixed         '
  write(321,'(A)') '  '
  write(321,'(A9,F7.3,A13,F7.3,A11)') ' xtrf =  ',xfoil_options%xtript,' (top)       ',xfoil_options%xtripb,' (bottom)  '
  write(321,'(A9,F7.3,A13,F7.3,A17,F7.3)') ' Mach =  ',mach(1),'     Re =    ',reynolds(1)/10.0**6,' e 6     Ncrit = ',ncrit_pt(1)
  write(321,'(A)') '  '
  write(321,'(A)') '   alpha    CL        CD       CDp       CM     Top_Xtr  Bot_Xtr'
  write(321,'(A)') '  ------ -------- --------- --------- -------- -------- --------'


! Display a summary of aero results

  write(*,*) 
  write(*,*) 'Aerodynamic information from Xfoil: '
  
  write(122,*) 
  write(122,*) 'Aerodynamic information from Xfoil: '
  
  write(*,*)
  write(*,'(6A14)') 'alpha(i)', 'lift(i)', 'drag(i)', 'moment(i)', 'xtrt(i)', 'xtrb(i)'
  
  write(122,*)
  write(122,'(6A14)') 'alpha(i)', 'lift(i)', 'drag(i)', 'moment(i)', 'xtrt(i)', 'xtrb(i)'
  do i = 1, noppoint

    write(text,*) i
    text = adjustl(text)

    if (viscrms(i) > 1.D-04) write(*,'(A)')                                  &
      ' Warning: operating point '//trim(text)//' did not converge.'
    write(*,'(6F14.8)')  alpha(i), lift(i), drag(i), moment(i), xtrt(i), xtrb(i)
    
    if (viscrms(i) > 1.D-04) write(122,'(A)')                                  &
      ' Warning: operating point '//trim(text)//' did not converge.'
    write(122,'(6F14.8)')  alpha(i), lift(i), drag(i), moment(i), xtrt(i), xtrb(i)
    
    if (viscrms(i) < 1.D-04) write(321,'(F8.3,F9.4,F10.5,F10.5,F9.4,F9.4,F9.4)') &
      alpha(i), lift(i), drag(i), drag(i), moment(i), xtrt(i), xtrb(i)

  end do

  close(122)
  close(321)
  
! Deallocate xfoil variables

  call xfoil_cleanup()

! Deallocate some things

  call deallocate_airfoil(foil)
  deallocate(alpha)
  deallocate(lift)
  deallocate(drag)
  deallocate(moment)
  deallocate(viscrms)
  deallocate(xtrt)
  deallocate(xtrb)

end program xfoil_only
