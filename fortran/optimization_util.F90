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

module optimization_util

! Module containing optimization routines

  implicit none

  contains

!=============================================================================80
!
! Initializes a random seed (subroutine from gcc.gnu.org)
!
!=============================================================================80
subroutine init_random_seed()

! For ifort compatibility
#ifdef intel_compilers
  use ifport, only : getpid  
#endif

  integer, dimension(:), allocatable :: myseed
  integer :: i, n, un, istat, dt(8), pid, t(2), s
  integer(8) :: count, tms
  
  call random_seed(size = n)
  allocate(myseed(n))

  ! First try if the OS provides a random number generator

  un = 18
  open(newunit=un, file="/dev/urandom", access="stream",                       &
       form="unformatted", action="read", status="old", iostat=istat)

  if (istat == 0) then

     read(un) myseed
     close(un)

  else

     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.

     call system_clock(count)
     if (count /= 0) then
        t = transfer(count, t)
     else
        call date_and_time(values=dt)
        tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24 * 60 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
        t = transfer(tms, t)
     end if

     s = ieor(t(1), t(2))
     pid = getpid() + 1099279 ! Add a prime
     s = ieor(s, pid)
     if (n >= 3) then
        myseed(1) = t(1) + 36269
        myseed(2) = t(2) + 72551
        myseed(3) = pid
        if (n > 3) then
           myseed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
        end if
     else
        myseed = s + 37 * (/ (i, i = 0, n - 1 ) /)
     end if

  end if

  call random_seed(put=myseed)
  deallocate(myseed)

end subroutine init_random_seed

!=============================================================================80
!
! Creates initial designs and tries to make them feasible, if desired
!
!=============================================================================80
subroutine initial_designs_mul(dv, objval, fevals, objfunc, xmin, xmax, use_x0,&
                           x0, feasible_init, feasible_limit, attempts,        &
                           message_codes, messages, constrain_matrix)
  use vardef, only : objfunction_type
    
  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objval
  integer, intent(out) :: fevals
  double precision, dimension(:), intent(in) :: xmin, xmax, x0
  logical, intent(in) :: use_x0, feasible_init
  double precision, intent(in) :: feasible_limit
  integer, intent(in) :: attempts
  integer, dimension(:), intent(inout) :: message_codes
  character(200), dimension(:), intent(inout) :: messages
  double precision, dimension(:,:), intent(inout) :: constrain_matrix
  integer, parameter :: CHUNK = 1

  interface
    type(objfunction_type) function objfunc(x,step)
      import :: objfunction_type
      double precision, dimension(:), intent(in) :: x
      integer, intent(in) :: step
    end function
  end interface

  integer :: i, pop, nvars, initcount
  character(30) :: text1, text2, text3
  double precision, dimension(:), allocatable :: randvec1, designstore
  double precision :: minstore
  type(objfunction_type) :: objfunction_return

! Initial settings and memory allocation

  nvars = size(dv,1)
  pop = size(dv,2)
  allocate(randvec1(nvars))
  allocate(designstore(nvars))

!$omp master

  fevals = pop

! Set up matrix of random numbers for initial designs
  
  call random_number(dv)

! Initial population of designs set between xmin and xmax
  
  write(*,*)
  write(*,*) 'Generating and evaluating initial designs ...'
  write(*,*)
  
!$omp end master
!$omp barrier

  if (use_x0) then
    dv(:,1) = x0
    objfunction_return = objfunc(x0,1)
    objval(1) = objfunction_return%value
    message_codes(1) = objfunction_return%message_code
    messages(1) = objfunction_return%message
    constrain_matrix(1,:) = objfunction_return%constrains_data
!$omp do SCHEDULE(DYNAMIC,CHUNK)
    do i = 2, pop
      dv(:,i) = (xmax - xmin)*dv(:,i) + xmin
      objfunction_return = objfunc(dv(:,i),1)
      objval(i) = objfunction_return%value
      message_codes(i) = objfunction_return%message_code
      messages(i) = objfunction_return%message
      constrain_matrix(i,:) = objfunction_return%constrains_data
    end do
!$omp end do
  else
!$omp do SCHEDULE(DYNAMIC,CHUNK)
    do i = 1, pop
      dv(:,i) = (xmax - xmin)*dv(:,i) + xmin
      objfunction_return = objfunc(dv(:,i),1)
      objval(i) = objfunction_return%value
      message_codes(i) = objfunction_return%message_code
      messages(i) = objfunction_return%message
      constrain_matrix(i,:) = objfunction_return%constrains_data
    end do
!$omp end do
  end if

! Enforce initially feasible designs

  if (feasible_init) then

    write(text1,*) attempts
    text1 = adjustl(text1)
!$omp master
    write(*,*) 'Checking initial designs for feasibility ...'
    write(*,*) '  (using a max of '//trim(text1)//' initialization attempts)'
    write(*,*)
!$omp end master

!$omp do SCHEDULE(DYNAMIC,CHUNK)
    do i = 1, pop 

      write(text2,*) i
      text2 = adjustl(text2)
      initcount = 0
      minstore = objval(i)
      designstore = dv(:,i)

!     Take a number of tries to fix infeasible designs

      do while ((initcount <= attempts) .and.                                  &
               (objval(i) >= feasible_limit))
        call random_number(randvec1)
        dv(:,i) = (xmax - xmin)*randvec1 + xmin
        objfunction_return = objfunc(dv(:,i),1)
        objval(i) = objfunction_return%value
        message_codes(i) = objfunction_return%message_code
        messages(i) = objfunction_return%message
        constrain_matrix(i,:) = objfunction_return%constrains_data
        if (objval(i) < minstore) then
          minstore = objval(i)
          designstore = dv(:,i)
        end if
        initcount = initcount + 1
!$omp critical
        fevals = fevals + 1
!$omp end critical
      end do

!     Pick the best design tested if a feasible design was not found

      if ((initcount > attempts) .and. (objval(i) >= feasible_limit)) then
        dv(:,i) = designstore
        objval(i) = minstore
      end if

!     Write a message about the feasibility of initial designs

      if ((initcount > attempts) .and. (objval(i) >= feasible_limit)) then
        write(*,*) ' Design '//trim(text2)//' is infeasible and was not'//     &
                   ' fixed within '//trim(text1)//' reinitialization attempts.'
      elseif ((initcount <= attempts) .and. (initcount > 0) .and.              &
              (objval(i) < feasible_limit)) then
        write(text3,*) initcount
        text3 = adjustl(text3)
        write(*,*) ' Design '//trim(text2)//' was initially infeasible but'//  &
                   ' was fixed after '//trim(text3)//                          &
                   ' reinitialization attempts.'
      else
        write(*,*) ' Design '//trim(text2)//' is feasible.' 
      end if

    end do

!$omp end do

!$omp master
    write(*,*)
!$omp end master

  end if

! Memory deallocation

  deallocate(randvec1)
  deallocate(designstore)

end subroutine initial_designs_mul

!=============================================================================80
!
! Creates initial designs and tries to make them feasible, if desired
!
!=============================================================================80
subroutine initial_designs_mul_2(dv, objval, fevals, objfunc, xmin, xmax, use_x0,&
                           x0, feasible_init, feasible_limit, attempts,        &
                           message_codes, messages, constrain_matrix, aero_matrix)
  use vardef, only : objfunction_type
    
  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objval
  integer, intent(out) :: fevals
  double precision, dimension(:), intent(in) :: xmin, xmax, x0
  logical, intent(in) :: use_x0, feasible_init
  double precision, intent(in) :: feasible_limit
  integer, intent(in) :: attempts
  integer, dimension(:), intent(inout) :: message_codes
  character(200), dimension(:), intent(inout) :: messages
  double precision, dimension(:,:), intent(inout) :: constrain_matrix
  double precision, dimension(:,:), intent(inout) :: aero_matrix
  integer, parameter :: CHUNK = 1

  interface
    type(objfunction_type) function objfunc(x,step)
      import :: objfunction_type
      double precision, dimension(:), intent(in) :: x
      integer, intent(in) :: step
    end function
  end interface

  integer :: i, pop, nvars, initcount
  character(30) :: text1, text2, text3
  double precision, dimension(:), allocatable :: randvec1, designstore
  double precision :: minstore
  type(objfunction_type) :: objfunction_return

! Initial settings and memory allocation

  nvars = size(dv,1)
  pop = size(dv,2)
  allocate(randvec1(nvars))
  allocate(designstore(nvars))

!$omp master

  fevals = pop

! Set up matrix of random numbers for initial designs
  
  call random_number(dv)

! Initial population of designs set between xmin and xmax
  
  write(*,*)
  write(*,*) 'Generating and evaluating initial designs ...'
  write(*,*)
  
!$omp end master
!$omp barrier

  if (use_x0) then
    dv(:,1) = x0
    objfunction_return = objfunc(x0,1)
    objval(1) = objfunction_return%value
    message_codes(1) = objfunction_return%message_code
    messages(1) = objfunction_return%message
    constrain_matrix(1,:) = objfunction_return%constrains_data
    aero_matrix(1,:) = objfunction_return%aero_data
!$omp do SCHEDULE(DYNAMIC,CHUNK)
    do i = 2, pop
      dv(:,i) = (xmax - xmin)*dv(:,i) + xmin
      objfunction_return = objfunc(dv(:,i),1)
      objval(i) = objfunction_return%value
      message_codes(i) = objfunction_return%message_code
      messages(i) = objfunction_return%message
      constrain_matrix(i,:) = objfunction_return%constrains_data
      aero_matrix(i,:) = objfunction_return%aero_data
    end do
!$omp end do
  else
!$omp do SCHEDULE(DYNAMIC,CHUNK)
    do i = 1, pop
      dv(:,i) = (xmax - xmin)*dv(:,i) + xmin
      objfunction_return = objfunc(dv(:,i),1)
      objval(i) = objfunction_return%value
      message_codes(i) = objfunction_return%message_code
      messages(i) = objfunction_return%message
      constrain_matrix(i,:) = objfunction_return%constrains_data
      aero_matrix(i,:) = objfunction_return%aero_data
    end do
!$omp end do
  end if

! Enforce initially feasible designs

  if (feasible_init) then

    write(text1,*) attempts
    text1 = adjustl(text1)
!$omp master
    write(*,*) 'Checking initial designs for feasibility ...'
    write(*,*) '  (using a max of '//trim(text1)//' initialization attempts)'
    write(*,*)
!$omp end master

!$omp do SCHEDULE(DYNAMIC,CHUNK)
    do i = 1, pop 

      write(text2,*) i
      text2 = adjustl(text2)
      initcount = 0
      minstore = objval(i)
      designstore = dv(:,i)

!     Take a number of tries to fix infeasible designs

      do while ((initcount <= attempts) .and.                                  &
               (objval(i) >= feasible_limit))
        call random_number(randvec1)
        dv(:,i) = (xmax - xmin)*randvec1 + xmin
        objfunction_return = objfunc(dv(:,i),1)
        objval(i) = objfunction_return%value
        message_codes(i) = objfunction_return%message_code
        messages(i) = objfunction_return%message
        constrain_matrix(i,:) = objfunction_return%constrains_data
        aero_matrix(i,:) = objfunction_return%aero_data
        if (objval(i) < minstore) then
          minstore = objval(i)
          designstore = dv(:,i)
        end if
        initcount = initcount + 1
!$omp critical
        fevals = fevals + 1
!$omp end critical
      end do

!     Pick the best design tested if a feasible design was not found

      if ((initcount > attempts) .and. (objval(i) >= feasible_limit)) then
        dv(:,i) = designstore
        objval(i) = minstore
      end if

!     Write a message about the feasibility of initial designs

      if ((initcount > attempts) .and. (objval(i) >= feasible_limit)) then
        write(*,*) ' Design '//trim(text2)//' is infeasible and was not'//     &
                   ' fixed within '//trim(text1)//' reinitialization attempts.'
      elseif ((initcount <= attempts) .and. (initcount > 0) .and.              &
              (objval(i) < feasible_limit)) then
        write(text3,*) initcount
        text3 = adjustl(text3)
        write(*,*) ' Design '//trim(text2)//' was initially infeasible but'//  &
                   ' was fixed after '//trim(text3)//                          &
                   ' reinitialization attempts.'
      else
        write(*,*) ' Design '//trim(text2)//' is feasible.' 
      end if

    end do

!$omp end do

!$omp master
    write(*,*)
!$omp end master

  end if

! Memory deallocation

  deallocate(randvec1)
  deallocate(designstore)

end subroutine initial_designs_mul_2
  
  
!=============================================================================80
!
! Creates initial designs and tries to make them feasible, if desired
!
!=============================================================================80
subroutine initial_designs(dv, objval, fevals, objfunc, xmin, xmax, use_x0,    &
                           x0, feasible_init, feasible_limit, attempts,        &
                           message_codes, messages)
  use vardef, only : objfunction_type
    
  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objval
  integer, intent(out) :: fevals
  double precision, dimension(:), intent(in) :: xmin, xmax, x0
  logical, intent(in) :: use_x0, feasible_init
  double precision, intent(in) :: feasible_limit
  integer, intent(in) :: attempts
  integer, dimension(:), intent(inout) :: message_codes
  character(200), dimension(:), intent(inout) :: messages
  integer, parameter :: CHUNK = 1

  interface
    type(objfunction_type) function objfunc(x,step)
      import :: objfunction_type
      double precision, dimension(:), intent(in) :: x
      integer, intent(in) :: step
    end function
  end interface

  integer :: i, pop, nvars, initcount
  character(30) :: text1, text2, text3
  double precision, dimension(:), allocatable :: randvec1, designstore
  double precision :: minstore
  type(objfunction_type) :: objfunction_return

! Initial settings and memory allocation

  nvars = size(dv,1)
  pop = size(dv,2)
  allocate(randvec1(nvars))
  allocate(designstore(nvars))

!$omp master

  fevals = pop

! Set up matrix of random numbers for initial designs
  
  call random_number(dv)

! Initial population of designs set between xmin and xmax
  
  write(*,*)
  write(*,*) 'Generating and evaluating initial designs ...'
  write(*,*)
  
!$omp end master
!$omp barrier

  if (use_x0) then
    dv(:,1) = x0
    objfunction_return = objfunc(x0,1)
    objval(1) = objfunction_return%value
    message_codes(1) = objfunction_return%message_code
    messages(1) = objfunction_return%message
!$omp do SCHEDULE(DYNAMIC,CHUNK)
    do i = 2, pop
      dv(:,i) = (xmax - xmin)*dv(:,i) + xmin
      objfunction_return = objfunc(dv(:,i),1)
      objval(i) = objfunction_return%value
      message_codes(i) = objfunction_return%message_code
      messages(i) = objfunction_return%message
    end do
!$omp end do
  else
!$omp do SCHEDULE(DYNAMIC,CHUNK)
    do i = 1, pop
      dv(:,i) = (xmax - xmin)*dv(:,i) + xmin
      objfunction_return = objfunc(dv(:,i),1)
      objval(i) = objfunction_return%value
      message_codes(i) = objfunction_return%message_code
      messages(i) = objfunction_return%message
    end do
!$omp end do
  end if

! Enforce initially feasible designs

  if (feasible_init) then

    write(text1,*) attempts
    text1 = adjustl(text1)
!$omp master
    write(*,*) 'Checking initial designs for feasibility ...'
    write(*,*) '  (using a max of '//trim(text1)//' initialization attempts)'
    write(*,*)
!$omp end master

!$omp do SCHEDULE(DYNAMIC,CHUNK)
    do i = 1, pop 

      write(text2,*) i
      text2 = adjustl(text2)
      initcount = 0
      minstore = objval(i)
      designstore = dv(:,i)

!     Take a number of tries to fix infeasible designs

      do while ((initcount <= attempts) .and.                                  &
               (objval(i) >= feasible_limit))
        call random_number(randvec1)
        dv(:,i) = (xmax - xmin)*randvec1 + xmin
        objfunction_return = objfunc(dv(:,i),1)
        objval(i) = objfunction_return%value
        message_codes(i) = objfunction_return%message_code
        messages(i) = objfunction_return%message
        if (objval(i) < minstore) then
          minstore = objval(i)
          designstore = dv(:,i)
        end if
        initcount = initcount + 1
!$omp critical
        fevals = fevals + 1
!$omp end critical
      end do

!     Pick the best design tested if a feasible design was not found

      if ((initcount > attempts) .and. (objval(i) >= feasible_limit)) then
        dv(:,i) = designstore
        objval(i) = minstore
      end if

!     Write a message about the feasibility of initial designs

      if ((initcount > attempts) .and. (objval(i) >= feasible_limit)) then
        write(*,*) ' Design '//trim(text2)//' is infeasible and was not'//     &
                   ' fixed within '//trim(text1)//' reinitialization attempts.'
      elseif ((initcount <= attempts) .and. (initcount > 0) .and.              &
              (objval(i) < feasible_limit)) then
        write(text3,*) initcount
        text3 = adjustl(text3)
        write(*,*) ' Design '//trim(text2)//' was initially infeasible but'//  &
                   ' was fixed after '//trim(text3)//                          &
                   ' reinitialization attempts.'
      else
        write(*,*) ' Design '//trim(text2)//' is feasible.' 
      end if

    end do

!$omp end do

!$omp master
    write(*,*)
!$omp end master

  end if

! Memory deallocation

  deallocate(randvec1)
  deallocate(designstore)

end subroutine initial_designs

!=============================================================================80
!
! Creates initial designs and tries to make them feasible, if desired
!
!=============================================================================80
subroutine initial_designs_original(dv, objval, fevals, objfunc, xmin, xmax, use_x0,    &
                           x0, feasible_init, feasible_limit, attempts)

  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objval
  integer, intent(out) :: fevals
  double precision, dimension(:), intent(in) :: xmin, xmax, x0
  logical, intent(in) :: use_x0, feasible_init
  double precision, intent(in) :: feasible_limit
  integer, intent(in) :: attempts

  interface
    double precision function objfunc(x)
      double precision, dimension(:), intent(in) :: x
    end function
  end interface

  integer :: i, pop, nvars, initcount
  character(30) :: text1, text2, text3
  double precision, dimension(:), allocatable :: randvec1, designstore
  double precision :: minstore

! Initial settings and memory allocation

  nvars = size(dv,1)
  pop = size(dv,2)
  allocate(randvec1(nvars))
  allocate(designstore(nvars))

!$omp master

  fevals = pop

! Set up matrix of random numbers for initial designs
  
  call random_number(dv)

! Initial population of designs set between xmin and xmax

  write(*,*) 'Generating and evaluating initial designs ...'
  write(*,*)
  
!$omp end master
!$omp barrier

  if (use_x0) then
    dv(:,1) = x0
    objval(1) = objfunc(x0)
!$omp do
    do i = 2, pop
      dv(:,i) = (xmax - xmin)*dv(:,i) + xmin
      objval(i) = objfunc(dv(:,i))
    end do
!$omp end do
  else
!$omp do
    do i = 1, pop
      dv(:,i) = (xmax - xmin)*dv(:,i) + xmin
      objval(i) = objfunc(dv(:,i))
    end do
!$omp end do
  end if

! Enforce initially feasible designs

  if (feasible_init) then

    write(text1,*) attempts
    text1 = adjustl(text1)
!$omp master
    write(*,*) 'Checking initial designs for feasibility ...'
    write(*,*) '  (using a max of '//trim(text1)//' initialization attempts)'
    write(*,*)
!$omp end master


!$omp do
    do i = 1, pop

      write(text2,*) i
      text2 = adjustl(text2)
      initcount = 0
      minstore = objval(i)
      designstore = dv(:,i)

!     Take a number of tries to fix infeasible designs

      do while ((initcount <= attempts) .and.                                  &
               (objval(i) >= feasible_limit))
        call random_number(randvec1)
        dv(:,i) = (xmax - xmin)*randvec1 + xmin
        objval(i) = objfunc(dv(:,i))
        if (objval(i) < minstore) then
          minstore = objval(i)
          designstore = dv(:,i)
        end if
        initcount = initcount + 1
!$omp critical
        fevals = fevals + 1
!$omp end critical
      end do

!     Pick the best design tested if a feasible design was not found

      if ((initcount > attempts) .and. (objval(i) >= feasible_limit)) then
        dv(:,i) = designstore
        objval(i) = minstore
      end if

!     Write a message about the feasibility of initial designs

      if ((initcount > attempts) .and. (objval(i) >= feasible_limit)) then
        write(*,*) ' Design '//trim(text2)//' is infeasible and was not'//     &
                   ' fixed within '//trim(text1)//' reinitialization attempts.'
      elseif ((initcount <= attempts) .and. (initcount > 0) .and.              &
              (objval(i) < feasible_limit)) then
        write(text3,*) initcount
        text3 = adjustl(text3)
        write(*,*) ' Design '//trim(text2)//' was initially infeasible but'//  &
                   ' was fixed after '//trim(text3)//                          &
                   ' reinitialization attempts.'
      else
        write(*,*) ' Design '//trim(text2)//' is feasible.' 
      end if

    end do

!$omp end do

!$omp master
    write(*,*)
!$omp end master

  end if

! Memory deallocation

  deallocate(randvec1)
  deallocate(designstore)

end subroutine initial_designs_original
                           
!=============================================================================80
!
! Computes max radius of designs (used for evaluating convergence
!
!=============================================================================80
function design_radius_simplex(dv)

  use math_deps, only : norm_2

  double precision, dimension(:,:), intent(in) :: dv
  double precision design_radius_simplex

  integer :: i, ndesigns
  double precision, dimension(size(dv,1)) :: design_centroid
  double precision :: radius
 
  ! Compute centroid of designs
  ndesigns = size(dv,2)
  design_centroid(:) = 0.d0
  do i = 1, ndesigns
    design_centroid = design_centroid + dv(:,i)
  end do
  design_centroid = design_centroid / dble(ndesigns)

  ! Compute max design radius

  design_radius_simplex = 0.d0
  do i = 1, ndesigns
    radius = norm_2(dv(:,i) - design_centroid)
    if (radius > design_radius_simplex) design_radius_simplex = radius
  end do
  
end function                           
                           
!=============================================================================80
!
! Computes max radius of designs (used for evaluating convergence
!
!=============================================================================80
function design_radius(dv, xmax, xmin)

  use math_deps, only : norm_2

  double precision, dimension(:,:), intent(in) :: dv
  double precision, dimension(:), intent(in) :: xmax, xmin
  double precision design_radius

  integer :: i, ndesigns
  double precision, dimension(size(dv,1)) :: design_centroid
  double precision, dimension(size(dv,1),size(dv,2)) :: dv_scaled
  double precision :: radius
 
  ndesigns = size(dv,2)
  
  ! Scale x to [-1,1]
  do i = 1, ndesigns
    dv_scaled(:,i)=-1+2*(dv(:,i)-xmin)/(xmax-xmin)
  end do
    
  ! Compute centroid of designs
  design_centroid(:) = 0.d0
  do i = 1, ndesigns
    design_centroid = design_centroid + dv_scaled(:,i)
  end do
  design_centroid = design_centroid / dble(ndesigns)

  !! Compute max design radius
  !
  !design_radius = 0.d0
  !do i = 1, ndesigns
  !  radius = norm_2(dv_scaled(:,i) - design_centroid)
  !  if (radius > design_radius) design_radius = radius
  !end do
  
  ! Compute mean design radius
  
  design_radius = 0.d0
  do i = 1, ndesigns
    radius = norm_2(dv_scaled(:,i) - design_centroid) / dble(size(dv,1))
    design_radius = design_radius + radius
  end do
  design_radius = design_radius / dble(ndesigns)
end function

!=============================================================================80
!
! Sorts a set of designs according to their objective function value
!
!=============================================================================80
subroutine bubble_sort(dv, objvals)

  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objvals

  double precision, dimension(size(dv,1),size(dv,2)) :: tempdv
  double precision, dimension(size(dv,2)) :: tempvals
  integer, dimension(size(dv,2)) :: finalorder, temporder
  integer :: nvars, ndesigns, i, sortcounter
  logical :: sorted

  nvars = size(dv,1)
  ndesigns = size(dv,2)

! Set up indexing array

  do i = 1, ndesigns
    finalorder(i) = i
  end do
  temporder = finalorder

! Bubble sorting algorithm

  sorted = .false.
  tempvals = objvals
  do while (.not. sorted)

    sortcounter = 0
    do i = 1, ndesigns - 1
      if (objvals(i+1) < objvals(i)) then

!       Flip the order of these elements. temp arrays are to preserve values.

        tempvals(i) = objvals(i+1)
        tempvals(i+1) = objvals(i)
        temporder(i) = finalorder(i+1)
        temporder(i+1) = finalorder(i)
        finalorder(i) = temporder(i)
        finalorder(i+1) = temporder(i+1)
        objvals(i) = tempvals(i)
        objvals(i+1) = tempvals(i+1)
        sortcounter = sortcounter + 1

      end if
    end do
    if (sortcounter == 0) sorted = .true.
    
  end do

! Use indexing array to rearrange order of designs

  do i = 1, ndesigns
    tempdv(:,i) = dv(:,finalorder(i))
  end do
  dv = tempdv

  end subroutine bubble_sort
  
!=============================================================================80
!
! Sorts a set of designs according to their objective function value
! while keeping the message and code
!
!=============================================================================80
subroutine bubble_sort_plus(dv, objvals, message_codes, messages)

  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objvals
  integer, dimension(:), intent(inout) :: message_codes
  character(200), dimension(:), intent(inout) :: messages

  double precision, dimension(size(dv,1),size(dv,2)) :: tempdv
  double precision, dimension(size(dv,2)) :: tempvals
  integer, dimension(size(dv,2)) :: tempmessage_codes
  character(200), dimension(size(dv,2)) :: tempmessages
  integer, dimension(size(dv,2)) :: finalorder, temporder
  integer :: nvars, ndesigns, i, sortcounter
  logical :: sorted

  nvars = size(dv,1)
  ndesigns = size(dv,2)

! Set up indexing array

  do i = 1, ndesigns
    finalorder(i) = i
  end do
  temporder = finalorder

! Bubble sorting algorithm

  sorted = .false.
  tempvals = objvals
  do while (.not. sorted)

    sortcounter = 0
    do i = 1, ndesigns - 1
      if (objvals(i+1) < objvals(i)) then

!       Flip the order of these elements. temp arrays are to preserve values.

        tempvals(i) = objvals(i+1)
        tempvals(i+1) = objvals(i)
        temporder(i) = finalorder(i+1)
        temporder(i+1) = finalorder(i)
        finalorder(i) = temporder(i)
        finalorder(i+1) = temporder(i+1)
        objvals(i) = tempvals(i)
        objvals(i+1) = tempvals(i+1)
        sortcounter = sortcounter + 1

      end if
    end do
    if (sortcounter == 0) sorted = .true.
    
  end do

! Use indexing array to rearrange order of designs

  do i = 1, ndesigns
    tempdv(:,i) = dv(:,finalorder(i))
    tempmessage_codes(i) = message_codes(finalorder(i))
    tempmessages(i) = messages(finalorder(i))
  end do
  dv = tempdv
  message_codes = tempmessage_codes
  messages = tempmessages

end subroutine bubble_sort_plus  
  
!=============================================================================80
!
! Sorts a set of designs according to their objective function value
! while keeping the message and code
!
!=============================================================================80
subroutine bubble_sort_plus_2(dv, objvals, message_codes, messages, constrain_matrix)

  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objvals
  integer, dimension(:), intent(inout) :: message_codes
  character(200), dimension(:), intent(inout) :: messages
  double precision, dimension(:,:), intent(inout) :: constrain_matrix

  double precision, dimension(size(dv,1),size(dv,2)) :: tempdv
  double precision, dimension(size(dv,2)) :: tempvals
  integer, dimension(size(dv,2)) :: tempmessage_codes
  character(200), dimension(size(dv,2)) :: tempmessages
  double precision, dimension(size(dv,2),size(constrain_matrix,2)) :: tempconstrain_matrix
  integer, dimension(size(dv,2)) :: finalorder, temporder
  integer :: nvars, ndesigns, i, sortcounter
  logical :: sorted

  nvars = size(dv,1)
  ndesigns = size(dv,2)

! Set up indexing array

  do i = 1, ndesigns
    finalorder(i) = i
  end do
  temporder = finalorder

! Bubble sorting algorithm

  sorted = .false.
  tempvals = objvals
  do while (.not. sorted)

    sortcounter = 0
    do i = 1, ndesigns - 1
      if (objvals(i+1) < objvals(i)) then

!       Flip the order of these elements. temp arrays are to preserve values.

        tempvals(i) = objvals(i+1)
        tempvals(i+1) = objvals(i)
        temporder(i) = finalorder(i+1)
        temporder(i+1) = finalorder(i)
        finalorder(i) = temporder(i)
        finalorder(i+1) = temporder(i+1)
        objvals(i) = tempvals(i)
        objvals(i+1) = tempvals(i+1)
        sortcounter = sortcounter + 1

      end if
    end do
    if (sortcounter == 0) sorted = .true.
    
  end do

! Use indexing array to rearrange order of designs

  do i = 1, ndesigns
    tempdv(:,i) = dv(:,finalorder(i))
    tempmessage_codes(i) = message_codes(finalorder(i))
    tempmessages(i) = messages(finalorder(i))
    tempconstrain_matrix(i,:) = constrain_matrix(finalorder(i),:)
  end do
  dv = tempdv
  message_codes = tempmessage_codes
  messages = tempmessages
  constrain_matrix = tempconstrain_matrix

end subroutine bubble_sort_plus_2

!=============================================================================80
!
! Sorts a set of designs according to their objective function value
! while keeping the message and code
!
!=============================================================================80
subroutine bubble_sort_plus_3(dv, objvals, message_codes, messages,            &
  constrain_matrix, aero_matrix)

  double precision, dimension(:,:), intent(inout) :: dv
  double precision, dimension(:), intent(inout) :: objvals
  integer, dimension(:), intent(inout) :: message_codes
  character(200), dimension(:), intent(inout) :: messages
  double precision, dimension(:,:), intent(inout) :: constrain_matrix
  double precision, dimension(:,:), intent(inout) :: aero_matrix

  double precision, dimension(size(dv,1),size(dv,2)) :: tempdv
  double precision, dimension(size(dv,2)) :: tempvals
  integer, dimension(size(dv,2)) :: tempmessage_codes
  character(200), dimension(size(dv,2)) :: tempmessages
  double precision, dimension(size(dv,2),size(constrain_matrix,2)) :: tempconstrain_matrix
  double precision, dimension(size(dv,2),size(aero_matrix,2)) :: tempaero_matrix
  integer, dimension(size(dv,2)) :: finalorder, temporder
  integer :: nvars, ndesigns, i, sortcounter
  logical :: sorted

  nvars = size(dv,1)
  ndesigns = size(dv,2)

! Set up indexing array

  do i = 1, ndesigns
    finalorder(i) = i
  end do
  temporder = finalorder

! Bubble sorting algorithm

  sorted = .false.
  tempvals = objvals
  do while (.not. sorted)

    sortcounter = 0
    do i = 1, ndesigns - 1
      if (objvals(i+1) < objvals(i)) then

!       Flip the order of these elements. temp arrays are to preserve values.

        tempvals(i) = objvals(i+1)
        tempvals(i+1) = objvals(i)
        temporder(i) = finalorder(i+1)
        temporder(i+1) = finalorder(i)
        finalorder(i) = temporder(i)
        finalorder(i+1) = temporder(i+1)
        objvals(i) = tempvals(i)
        objvals(i+1) = tempvals(i+1)
        sortcounter = sortcounter + 1

      end if
    end do
    if (sortcounter == 0) sorted = .true.
    
  end do

! Use indexing array to rearrange order of designs

  do i = 1, ndesigns
    tempdv(:,i) = dv(:,finalorder(i))
    tempmessage_codes(i) = message_codes(finalorder(i))
    tempmessages(i) = messages(finalorder(i))
    tempconstrain_matrix(i,:) = constrain_matrix(finalorder(i),:)
    tempaero_matrix(i,:) = aero_matrix(finalorder(i),:)
  end do
  dv = tempdv
  message_codes = tempmessage_codes
  messages = tempmessages
  constrain_matrix = tempconstrain_matrix
  aero_matrix = tempaero_matrix

end subroutine bubble_sort_plus_3
  
  

!=============================================================================80
!
! Pops item out of a vetor.  Note: doesn't actually change size of vector, just
! shuffles data so that the first nitems-1 entries represent the new vector.
!
!=============================================================================80
subroutine pop_double_vector(vector, nitems, popidx)

  double precision, dimension(:), intent(inout) :: vector
  integer, intent(in) :: nitems, popidx

  integer :: i
  double precision, dimension(size(vector,1)) :: tempvector

! Copy input vector

  tempvector = vector

! Populate output vector

  do i = 1, nitems-1
    if (i < popidx) then
      vector(i) = tempvector(i)
    else
      vector(i) = tempvector(i+1)
    end if
  end do

end subroutine pop_double_vector

!=============================================================================80
!
! Pops item out of a vetor.  Note: doesn't actually change size of vector, just
! shuffles data so that the first nitems-1 entries represent the new vector.
!
!=============================================================================80
subroutine pop_integer_vector(vector, nitems, popidx)

  integer, dimension(:), intent(inout) :: vector
  integer, intent(in) :: nitems, popidx

  integer :: i
  integer, dimension(size(vector,1)) :: tempvector

! Copy input vector

  tempvector = vector

! Populate output vector

  do i = 1, nitems-1
    if (i < popidx) then
      vector(i) = tempvector(i)
    else
      vector(i) = tempvector(i+1)
    end if
  end do

end subroutine pop_integer_vector

!=============================================================================80
!
! Writes design variables to file
!
!=============================================================================80
subroutine write_design(filename, filestat, variables, counter)

  character(*), intent(in) :: filename, filestat
  double precision, dimension(:), intent(in) :: variables
  integer, intent(in) :: counter

  integer, save :: iunit
  integer :: nvars, i
  character(30) :: text

  nvars = size(variables,1)
  iunit = 17

! Open the file and write to it if requested

  if (trim(filestat) == 'new') then
    open(unit=iunit, file=filename, status='replace')
    write(text,*) nvars
    text = adjustl(text)
    write(iunit,'(A)') 'Number of variables: '//trim(text)
  else
    open(unit=iunit, file=filename, status='old', position='append')
  end if

! Write iteration number and the design variables to file

  write(text,*) counter
  text = adjustl(text)
  write(iunit,'(A)') 'Design number '//trim(text)
  do i = 1, nvars
    write(iunit,'(es25.16)') variables(i)
  end do

! Close the file 

  close(iunit)

end subroutine write_design

!=============================================================================80
!
! Reads commands from run_control file and clears any unrecognized commands
!
!=============================================================================80
subroutine read_run_control(commands, ncommands)
  
  use vardef, only : output_prefix  

  character(80), dimension(:), intent(inout) :: commands
  integer, intent(out) :: ncommands

  character(80) :: buffer
  integer :: rcunit, ioerr, i
  
  commands(:) = ""
  ncommands = 0

  rcunit = 18
  open(unit=rcunit, file='run_control_'//trim(output_prefix), status='old', iostat=ioerr, err=501)
  if (ioerr /= 0) then
    return
  else
    do while (1 .eq. 1)
      read(rcunit,'(A)',end=501) buffer
      if ( (trim(buffer) == "stop") .or.                                       &
           (trim(buffer) == "stop_monitoring") ) then
        commands(ncommands+1) = buffer
        ncommands = ncommands + 1
      else
        write(*,*) 'Warning: unrecognized command "'//trim(buffer)//           &
                   '" in run_control.' 
      end if
    end do
  end if
   
  write(*,*) "Warning: error encountered while reading run_control. Skipping."
  return
501 close(rcunit)

  open(unit=rcunit, file='run_control_'//trim(output_prefix), status='replace', err=502)
  do i = 1, ncommands
    write(rcunit,'(A)') commands(ncommands)
  end do
  close(rcunit) 

  return

502 write(*,*) "Warning: error encountered while reading run_control. Skipping."

end subroutine read_run_control
  
subroutine write_progress_cmd(step, f0, fmin, radius, relative_fmin_report)
  
  integer, intent(in) :: step
  double precision, intent(in) :: f0, fmin, radius
  logical, intent(in) :: relative_fmin_report
  
  write(*,*)
  write(*,'(A12,I5)')   ' Iteration: ', step
  write(*,'(A27,F12.6)') '   Objective function:    ', fmin
  if (relative_fmin_report) write(*,'(A27,F12.6,A1)')             &
                      '   Improvement over seed: ', (f0 - fmin)/f0*100.d0, '%'
  write(*,'(A27,ES13.3)') '   Design radius:         ', radius

end subroutine write_progress_cmd
  
subroutine write_progress_per_eval(i, control)

  use vardef, only : progress_per_eval
  use, intrinsic :: iso_fortran_env, only: OUTPUT_UNIT
  
  integer, intent(in) :: i, control 
    
  if (progress_per_eval .NE. 'none') then
    if (progress_per_eval .EQ. 'full' .AND. control .EQ. 1) then
      write(*,'(I5)', advance='no') i
      write(*,'(A22)', advance='no') 'Start function, '
    end if
    if (progress_per_eval .EQ. 'full' .AND. control .EQ. 2) then
      write(*,'(I5)', advance='no') i
      write(*,'(A22)', advance='no') 'end the function, '
    end if
    if (progress_per_eval .EQ. 'full' .AND. control .EQ. 3) then
      write(*,'(I5)', advance='no') i
      write(*,'(A22)') 'assign values.'
    end if
    if (progress_per_eval .EQ. 'part' .AND. control .EQ. 3) then
      write(*,'(A1)', advance='no') '.'
    end if
  end if

end subroutine write_progress_per_eval
  
subroutine write_history(step, f0, fmin, radius, relative_fmin_report,         &
  stepstart, steptime, restarttime, control)
  
  use vardef, only : output_prefix

  integer, intent(in) :: step
  double precision, intent(in) :: f0, fmin, radius
  logical, intent(in) :: relative_fmin_report
  integer, intent(in) :: stepstart, steptime, restarttime
  integer, intent(in) :: control
  
  integer :: iunit, ioerr
  logical :: new_history_file
  character(100) :: histfile
  character(14) :: timechar
  character(11) :: stepchar
  character(20) :: fminchar
  character(15) :: radchar
  character(25) :: relfminchar
  
  histfile = trim(output_prefix)//'_optimization_history.dat'
  iunit = 17
  
  ! Only header
  if (control .EQ. 1) then
  ! Open file for writing iteration history
    new_history_file = .false.
    if (step == 1) then
      new_history_file = .true.
    else
      open(unit=iunit, file=histfile, status='old',            &
           position='append', iostat=ioerr)
      if (ioerr /= 0) then
        write(*,*) 
        write(*,*) "Warning: did not find existing optimization_history.dat file."
        write(*,*) "A new one will be written, but old data will be lost."
        write(*,*)
        new_history_file = .true.
      end if
    end if
    if (new_history_file) then
      open(unit=iunit, file=histfile, status='replace')
      if (relative_fmin_report) then
        write(iunit,'(A)') "Iteration  Objective function  "//&
                           "% Improvement over seed  Design radius"//&
                           "  Time (seconds)"
      else
        write(iunit,'(A)') "Iteration  Objective function  Design radius"//&
                           "  Time (seconds)"
      end if
      flush(iunit)
    end if
    close(iunit)
  end if
  
  !Only values
  if (control .EQ. 2) then
    open(unit=iunit, file=histfile, status='old', position='append',           &
      iostat=ioerr)
    if (ioerr /= 0) then
      write(*,*) 
      write(*,*) "Warning: did not find existing optimization_history.dat file."
      write(*,*) "A new one will be written, but old data will be lost."
      write(*,*)
      new_history_file = .true.
    end if
    
    write(stepchar,'(I11)') step
    write(fminchar,'(F14.10)') fmin
    write(radchar,'(ES14.6)') radius
    write(timechar,'(I14)') (steptime-stepstart)+restarttime
    if (relative_fmin_report) then
      write(relfminchar,'(F14.10)') (f0 - fmin)/f0*100.d0
      write(iunit,'(A11,A20,A25,A15,A14)') adjustl(stepchar),                  &
        adjustl(fminchar), adjustl(relfminchar), adjustl(radchar),             &
        adjustl(timechar)
    else
      write(iunit,'(A11,A20,A15,A14)') adjustl(stepchar), adjustl(fminchar),   &
        adjustl(radchar), adjustl(timechar)
    end if
    flush(iunit)
    
    close(iunit)
  end if
  
end subroutine write_history


!=============================================================================80
!
! Dvs write routine
!
!=============================================================================80
subroutine write_dvs(name_len, name, step, dv, objval, message_codes, messages,&
                         constrain_matrix, aero_matrix, x0, f0, xopt, fmin)

  use vardef, only : output_prefix, dvs_for_type, naddthickconst, noppoint,    &
                     ndrag_constrain, nmoment_constrain, nlift_constrain

  integer, intent(in) :: name_len
  character(name_len), intent(in)  :: name
  integer, intent(in) :: step
  double precision, intent(in) :: fmin, f0
  double precision, dimension(:), intent(in) ::  x0, xopt, objval
  double precision, dimension(:,:), intent(in) :: dv
  integer, dimension(:), intent(in) :: message_codes
  character(200), dimension(:), intent(in) :: messages
  double precision, dimension(:,:), intent(in) :: constrain_matrix
  double precision, dimension(:,:), intent(in) :: aero_matrix

  integer :: i,j
  
  character(100) :: dvsfile, constfile, aerofile, text, textdv
  integer :: iunit
  
  ! Status notification

  dvsfile = 'dvs_'//trim(name)//'_'//trim(output_prefix)//'.dat'
  constfile = 'constrain_'//trim(name)//'_'//trim(output_prefix)//'.dat'
  aerofile = 'aero_'//trim(name)//'_'//trim(output_prefix)//'.dat'
  
  write(*,'(A)',advance='no') '  Writing '//trim(name)//' log data to files '//trim(dvsfile)
  if (size(constrain_matrix,2) .NE. 1) then
    write(*,'(A)', advance='no')' and '//trim(constfile)
    write(*,'(A)', advance='no')' and '//trim(aerofile)//' ...'
  end if
  write(*,*)
  
  ! Log file
  ! Open files and write headers, if necessary
  iunit = 13
  write(text,*) step
  text = adjustl(text)
  
  if (step == 1) then

    !   Header for dvs file

    open(unit=iunit, file=dvsfile, status='replace')
    write(iunit,'(A)') 'title="Design variables file"'
    write(iunit,'(A)') 'variables="dvs vector", "objval", "message_code", "message"'
    write(iunit,'(A)') 'step = 0: x0'
    if (dvs_for_type(2) .NE. 0) then
      do i=1, dvs_for_type(2)
        write(textdv,*) i 
        textdv=adjustl(textdv)
        write(iunit,'(A12)', advance='no') 'top - '//trim(textdv)
      end do
    end if
    if (dvs_for_type(3) .NE. 0) then
      do i=1, dvs_for_type(3)
        write(textdv,*) i 
        textdv=adjustl(textdv)
        write(iunit,'(A12)', advance='no') 'bot - '//trim(textdv)
      end do
    end if
    if (dvs_for_type(4) .NE. 0) then
      do i=1, dvs_for_type(4)
        write(textdv,*) i 
        textdv=adjustl(textdv)
        write(iunit,'(A12)', advance='no') 'fdef - '//trim(textdv)
      end do
    end if
    if (dvs_for_type(5) .NE. 0) write(iunit,'(A12)', advance='no') 'fhinge'
    if (dvs_for_type(6) .NE. 0) write(iunit,'(A12)', advance='no') 'te_thick'
    write(iunit,'(A30)', advance='no') 'f0'
    write(iunit,'(A)') ' '
  
    do i = 1, size(dv,1)
      write(iunit,'(F12.8)', advance='no') x0(i)
    end do
    write(iunit,'(F30.8)', advance='no') f0
    write(iunit,'(A)') ' '
    
  else

    !   Open dvs file and write zone header

    open(unit=iunit, file=dvsfile, status='old', position='append', err=900)

  end if

  ! Write coordinates to file
  write(iunit,'(A)') 'step = '//trim(text)//': xopt '
  if (dvs_for_type(2) .NE. 0) then
    do i=1, dvs_for_type(2)
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A12)', advance='no') 'top - '//trim(textdv)
    end do
  end if
  if (dvs_for_type(3) .NE. 0) then
    do i=1, dvs_for_type(3)
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A12)', advance='no') 'bot - '//trim(textdv)
    end do
  end if
  if (dvs_for_type(4) .NE. 0) then
    do i=1, dvs_for_type(4)
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A12)', advance='no') 'fdef - '//trim(textdv)
    end do
  end if
  if (dvs_for_type(5) .NE. 0) write(iunit,'(A12)', advance='no') 'fhinge'
  if (dvs_for_type(6) .NE. 0) write(iunit,'(A12)', advance='no') 'te_thick'
  write(iunit,'(A30)', advance='no') 'fmin'
  write(iunit,'(A)') ' '

  do i =1, size(dv,1)
    write(iunit,'(F12.8)', advance='no') xopt(i)
  end do
  write(iunit,'(F30.8)', advance='no') fmin
  write(iunit,'(A)') ' '
  
  write(iunit,'(A)') 'step = '//trim(text)//': dv '
  if (dvs_for_type(2) .NE. 0) then
    do i=1, dvs_for_type(2)
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A12)', advance='no') 'top - '//trim(textdv)
    end do
  end if
  if (dvs_for_type(3) .NE. 0) then
    do i=1, dvs_for_type(3)
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A12)', advance='no') 'bot - '//trim(textdv)
    end do
  end if
  if (dvs_for_type(4) .NE. 0) then
    do i=1, dvs_for_type(4)
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A12)', advance='no') 'fdef - '//trim(textdv)
    end do
  end if
  if (dvs_for_type(5) .NE. 0) write(iunit,'(A12)', advance='no') 'fhinge'
  if (dvs_for_type(6) .NE. 0) write(iunit,'(A12)', advance='no') 'te_thick'
  write(iunit,'(A30)', advance='no') 'objval'
  write(iunit,'(A14)', advance='no') 'message_code'
  write(iunit,'(A)', advance='no') ' message'
  write(iunit,'(A)') ' '
  
  do i = 1, size(dv,2)
    do j =1, size(dv,1)
      write(iunit,'(F12.8)', advance='no') dv(j,i)
    end do
    write(iunit,'(F30.8)', advance='no') objval(i)
    write(iunit,'(I14)', advance='no') message_codes(i)
    write(iunit,'(A)', advance='no') messages(i)
    write(iunit,'(A)') ' '
  end do

  ! Close output files

  close(iunit)

  ! Status notification
  
  if (size(constrain_matrix,2) .NE. 1) then
    ! Constrains file
    ! Open files and write headers, if necessary
    iunit = 14
    write(text,*) step
    text = adjustl(text)
  
    if (step == 1) then

      !   Header for Constrains file

      open(unit=iunit, file=constfile, status='replace')
      write(iunit,'(A)') 'title="Constrains file"'
      write(iunit,'(A)') 'variables="constrains", "objval", "message_code", "message"'
    
    else

      !   Open dvs file and write zone header

      open(unit=iunit, file=constfile, status='old', position='append', err=901)

    end if

    ! Write coordinates to file
    write(iunit,'(A)') 'step = '//trim(text)//': constrains vector '
    write(iunit,'(A14)', advance='no') 'tcTE'
    if (dvs_for_type(4) .NE. 0) then
      do i=1, dvs_for_type(4)
        write(textdv,*) i 
        textdv=adjustl(textdv)
        write(iunit,'(A14)', advance='no') 'flap_deg - '//trim(textdv)
      end do
    end if
    write(iunit,'(A14)', advance='no') 'flap_hinge'
    write(iunit,'(A14)', advance='no') 'minthick'
    write(iunit,'(A14)', advance='no') 'maxthick'
    
    if (naddthickconst .NE. 0) then
      do i=1, naddthickconst
        write(textdv,*) i 
        textdv=adjustl(textdv)
        write(iunit,'(A14)', advance='no') 'addthickconst - '//trim(textdv)
      end do
    end if
    
    write(iunit,'(A14)', advance='no') 'TE_angle'
    write(iunit,'(A14)', advance='no') 'maxcamb'
    write(iunit,'(A14)', advance='no') 'maxpang'
    write(iunit,'(A14)', advance='no') 'minang'
    write(iunit,'(A14)', advance='no') 'difang'
    write(iunit,'(A14)', advance='no') 'rev_t'
    write(iunit,'(A14)', advance='no') 'rev_b'
    write(iunit,'(A14)', advance='no') 'pan_ang'
    write(iunit,'(A14)', advance='no') 'maxgrowth'
    write(iunit,'(A14)', advance='no') 'n_unconv'
    
    if (nmoment_constrain .NE. 0) then
      do i=1, nmoment_constrain
        write(textdv,*) i 
        textdv=adjustl(textdv)
        write(iunit,'(A14)', advance='no') 'moment - '//trim(textdv)
      end do
    end if
    if (nlift_constrain .NE. 0) then
      do i=1, nlift_constrain
        write(textdv,*) i 
        textdv=adjustl(textdv)
        write(iunit,'(A14)', advance='no') 'lift - '//trim(textdv)
      end do
    end if
        if (ndrag_constrain .NE. 0) then
      do i=1, ndrag_constrain
        write(textdv,*) i 
        textdv=adjustl(textdv)
        write(iunit,'(A14)', advance='no') 'drag - '//trim(textdv)
      end do
    end if
    write(iunit,'(A30)', advance='no') 'objval'
    write(iunit,'(A14)', advance='no') 'message_code'
    write(iunit,'(A)', advance='no') ' message'
    write(iunit,'(A)') ' '
  
    do i = 1, size(constrain_matrix,1)
      do j = 1, size(constrain_matrix,2)
        write(iunit,'(F14.8)', advance='no') constrain_matrix(i,j)
      end do
      write(iunit,'(F30.8)', advance='no') objval(i)
      write(iunit,'(I14)', advance='no') message_codes(i)
      write(iunit,'(A)', advance='no') messages(i)
      write(iunit,'(A)') ' '
    end do

    ! Close output files

    close(iunit)
  end if

  
  if (size(constrain_matrix,2) .NE. 1) then
    ! Aero file
    ! Open files and write headers, if necessary
    iunit = 15
    write(text,*) step
    text = adjustl(text)
  
    if (step == 1) then

      !   Header for aero file

      open(unit=iunit, file=aerofile, status='replace')
      write(iunit,'(A)') 'title="Aero file"'
      write(iunit,'(A)') 'variables="Aerodynamic proprerties", '//             &
        &'"partial objective function", "partial improvement", '//             &
        &'"objval", "message_code", "message"'
    
    else

      !   Open aero file and write zone header

      open(unit=iunit, file=aerofile, status='old', position='append', err=902)

    end if

    ! Write coordinates to file
    write(iunit,'(A)') 'step = '//trim(text)//': aero vector '
    do i=1, noppoint
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A14)', advance='no') 'aero prop - '//trim(textdv)
    end do
    do i=1, noppoint
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A14)', advance='no') 'part obj - '//trim(textdv)
    end do
    do i=1, noppoint
      write(textdv,*) i 
      textdv=adjustl(textdv)
      write(iunit,'(A14)', advance='no') 'part imp - '//trim(textdv)
    end do
    write(iunit,'(A14)', advance='no') 'obj no penal'

    write(iunit,'(A30)', advance='no') 'objval'
    write(iunit,'(A14)', advance='no') 'message_code'
    write(iunit,'(A)', advance='no') ' message'
    write(iunit,'(A)') ' '
  
    do i = 1, size(aero_matrix,1)
      do j = 1, size(aero_matrix,2)
        write(iunit,'(F14.8)', advance='no') aero_matrix(i,j)
      end do
      write(iunit,'(F30.8)', advance='no') objval(i)
      write(iunit,'(I14)', advance='no') message_codes(i)
      write(iunit,'(A)', advance='no') messages(i)
      write(iunit,'(A)') ' '
    end do

    ! Close output files

    close(iunit)
  end if
  
    ! Status notification

    write(*,*) '  Successfully wrote '//trim(name)//' log data.'
  
  return
  
  ! Warning if there was an error opening dvs file

  900 write(*,*) "Warning: unable to open "//trim(dvsfile)//". Skipping ..."
  return
  901 write(*,*) "Warning: unable to open "//trim(constfile)//". Skipping ..."
  return
  902 write(*,*) "Warning: unable to open "//trim(aerofile)//". Skipping ..."
  return
        
  
end subroutine write_dvs
  
end module optimization_util
