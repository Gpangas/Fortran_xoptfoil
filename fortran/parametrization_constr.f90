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
  
module parametrization_constr

! Contains subroutines to create an airfoil shape from design variables

  implicit none
! Shape functions for creating airfoil shapes (top and bottom)

  double precision, dimension(:,:), allocatable :: top_b_matrix
  double precision, dimension(:,:), allocatable :: bot_b_matrix
  
  contains

!/////////////////////////////////////////////////////////////////////////////80
!
! KBP subroutines
!
!/////////////////////////////////////////////////////////////////////////////80

! ----------------------------------------------------------------------------80
!
! Subroutine that implements the Kulfan-Bussoletti Parameterization,     
! weights to coordinates.
!
! ----------------------------------------------------------------------------80
  
  subroutine KBP_airfoil(Xu_seed, Zu_seed, Xl_seed, Zl_seed, Wu, Wl, Zu, Zl,   &
    symmetrical, tTE)
  use vardef, only : int_kulfan_bussoletti_LEM
  
  implicit none

  real*8, intent (in) :: tTE
  real*8, dimension(:), intent (in) :: Wu ! weights
  real*8, dimension(:), intent (in) :: Wl ! weights
  logical, intent(in) :: symmetrical
  real*8, dimension(:), intent (in) :: Xu_seed, Zu_seed, Xl_seed, Zl_seed
  real*8, dimension(size(Xu_seed,1)), intent (out) :: Zu
  real*8, dimension(size(Xl_seed,1)), intent (out) :: Zl
  integer :: i, nPointst, nPointsb
  integer :: nu, nl ! degree of polynomial

  nu=size(Wu,1)-1-int_kulfan_bussoletti_LEM
  nl=size(Wu,1)-1-int_kulfan_bussoletti_LEM
  nPointst=size(Xu_seed,1)
  nPointsb=size(Xl_seed,1)
  
  ! Calculates the Z coordinate of the airfoil
  do i=1,nPointst
    call KBP_Point(Wu, tTE/2.0d0, nu, Xu_seed(i), Zu(i),                       &
                                                      int_kulfan_bussoletti_LEM)
  end do
  if (.not. symmetrical) then
    do i=1,nPointsb
      call KBP_Point(Wl, -tTE/2.0d0, nl, Xl_seed(i), Zl(i),                    &
                                                      int_kulfan_bussoletti_LEM)
    end do
    
  ! For symmetrical airfoils, just mirror the top surface
  else
    do i = 1, npointsb
      Zl(i) = -Zu(i)
    end do
  end if

  
end subroutine  KBP_airfoil
  
! ----------------------------------------------------------------------------80
!
! Subroutine that computes one point in the the Kulfan-Bussoletti 
! Parameterization.
!
! ----------------------------------------------------------------------------80

subroutine KBP_Point(weight,tTE,BPO,x,z, int_LEM)
  use math_deps, only : surface_function
  
  implicit none

  integer, intent (in) ::       BPO            ! BPO->Bernstein Polynomial Order
  real*8, intent (in) ::        x, tTE         ! non-dimenstional x and z_TE
  real*8, intent (out) ::       z              ! non-dimenstional z
  integer, intent (in) ::       int_LEM        ! use LEM? yes (1) or no (0)
  real*8, dimension(BPO+1+int_LEM) ::   weight ! weights

  z         = (x**0.5)*(1.0d0-x)*surface_function(x,weight,BPO, int_LEM)+x*tTE

end subroutine KBP_Point

! ----------------------------------------------------------------------------80
!
! Subroutine that implements the Kulfan-Bussoletti Parameterization, 
! coordinates to weights.
!
! ----------------------------------------------------------------------------80
subroutine KBP_init(Xu, Zu, Xl, Zl, Wu, Wl)
  
  use vardef, only : int_kulfan_bussoletti_LEM

  implicit none

  real*8, dimension(:), intent(in) ::    Xu, Xl, Zu, Zl    ! airfoil coordinates
  real*8, dimension(:), intent(inout) :: Wu ! weights of Bernstein polinomials
  real*8, dimension(:), intent(inout) :: Wl ! weights of Bernstein polinomials

  integer :: nu, nl ! n is the order of the polynomial
  integer :: mu, ml
  real*8 ::  tTE

  ! Set the value of zcTE, dimensionless trailing edge thickness
  mu       = size(Zu,1)
  ml       = size(Zl,1)
  tTE    = Zu(mu)-Zl(ml)

  nu = size(Wu)-1-int_kulfan_bussoletti_LEM
  nl = size(Wl)-1-int_kulfan_bussoletti_LEM
  
  ! Set the weights of the Bernstein Polynomials.
  call KBParameterization_fitting(mu,Xu,Zu,nu,tTE/2.0d0,                       &
                                                   int_kulfan_bussoletti_LEM,Wu)
  call KBParameterization_fitting(ml,Xl,Zl,nl,-tTE/2.0d0,                      &
                                                   int_kulfan_bussoletti_LEM,Wl)

end subroutine KBP_init

! ----------------------------------------------------------------------------80
!
! Fits a KB Parameterization of order n to the XZ coordinates of 
! a given airfoil (top or bottom half) using least squares to 
! calculate the weights.
!
! (input) ndata     -> no. of data points
! (input) X, Z      -> airfoil coordinates
! (input) zcTE      -> dimensionless trailing edge coordinate
! (input) n         -> degree of polynomial
! (input) int_LEM   -> whether to use LEM: yes (1) or no (0)
! (input/output) W  -> weight vector 
!
! ----------------------------------------------------------------------------80

subroutine KBParameterization_fitting(ndata,X,Z,n,zcTE,int_LEM,W)

  use math_deps, only : fx, BINCOEF

  implicit none
  
  integer, intent(in) ::                                  ndata
  double precision, intent (in), dimension(ndata) ::      X, Z
  double precision, intent (in) ::                        zcTE
  integer, intent (in) ::                                 n
  integer, intent (in) ::                                 int_LEM
  double precision, intent (inout), dimension(:) ::       W
  
  double precision, dimension(ndata,n+1+int_LEM) ::       G !G.W=F
  double precision, dimension(ndata) ::                   F
  
  double precision, dimension(n+1+int_LEM,ndata) ::       G_trans
  
  double precision, dimension(n+1+int_LEM,n+1+int_LEM) :: A
  double precision, dimension(n+1+int_LEM) ::             B
  
  double precision, dimension(n+1+int_LEM,1) ::           B_send
  double precision, dimension(n+1+int_LEM) ::             S
  double precision ::                                     RCOND
  integer ::                                              RANK, LWORK, INFO
  double precision, dimension(5*(n+1+int_LEM)) ::         WORK
  
  integer ::                                              i,j
    
  F  = 0.0d0
  G  = 0.0d0  

  do i=1,ndata
    F(i)=Z(i)-X(i)*zcTE
    do j=0,n
      if ( (X(i) .EQ. 0.d0) .or. (X(i) .EQ. 1.d0)) then
        G(i,j+1) = 0.0d0
      else
        G(i,j+1) = fx(X(i)) * BINCOEF(n,j) * (X(i)**j)*((1.0d0-X(i))**(n-j))
      end if
    end do
    if (int_LEM .EQ. 1) then
      G(i,n+1+int_LEM) = fx(X(i)) * (X(i)**0.5d0)*((1.0d0-X(i))**(n-0.5d0))
    end if
  end do
    
  A  = 0.0d0
  B  = 0.0d0

  G_trans=transpose(G)

  A=matmul(G_trans, G)
  B=matmul(G_trans, F)
   
  B_send(:,1) = B(:)
  RCOND = 1.0e-6
  LWORK = 5*(n+1+int_LEM)
  call DGELSS(n+1+int_LEM, n+1+int_LEM, 1, A, n+1+int_LEM, B_send, n+1+int_LEM,&
    S, RCOND, RANK, WORK, LWORK, INFO )
  
  W = B_send(:,1)
   
end subroutine KBParameterization_fitting

!/////////////////////////////////////////////////////////////////////////////80
!
! BSP subroutines
!
!/////////////////////////////////////////////////////////////////////////////80

! ----------------------------------------------------------------------------
! Subroutine that implements the b-spline Parameterization, weights to coordinates.
subroutine BSP_airfoil(ut_seed, xt_seed, zt_seed, ub_seed, xb_seed, zb_seed,   &
xmodest, xmodesb, modest, modesb, zt_new, zb_new,            &
symmetrical , tTE)
  
  use vardef, only : b_spline_xtype, b_spline_degree
  use math_deps, only : sort_vector

  implicit none

  real*8, dimension(:), intent (in) :: xmodest, xmodesb, modest, modesb
  logical, intent(in) :: symmetrical
  real*8, intent(in) :: tTE
  real*8, dimension(:), intent (in) ::ut_seed, xt_seed, zt_seed
  real*8, dimension(:), intent (in) ::ub_seed, xb_seed, zb_seed
  real*8, dimension(size(xt_seed,1)), intent (out) :: zt_new
  real*8, dimension(size(xb_seed,1)), intent (out) :: zb_new
  
  real*8, dimension(size(xt_seed,1)) :: xt_new
  real*8, dimension(size(xb_seed,1)) :: xb_new
  integer :: nPointst, nPointsb, nmodest, nmodesb
  real*8, dimension(:), allocatable :: zmodest_use
  real*8, dimension(:), allocatable :: zmodesb_use
  real*8, dimension(:), allocatable :: xmodest_use
  real*8, dimension(:), allocatable :: xmodesb_use
  real*8, dimension(size(xt_seed,1)) :: ut_use
  real*8, dimension(size(xb_seed,1)) :: ub_use

  nmodest=size(xmodest,1)
  nmodesb=size(xmodesb,1)
  
  nPointst=size(xt_seed,1)
  nPointsb=size(xb_seed,1)
  
  if (b_spline_xtype .EQ. 1) then
    allocate(zmodest_use(nmodest),zmodesb_use(nmodesb))
    
    ! Construst full control points vector
    zmodest_use(1)=0.0d0
    zmodest_use(nmodest)=+tTE/2.0d0
    zmodest_use(2:nmodest-1)=modest
  
    zmodesb_use(1)=0.0d0
    zmodesb_use(nmodesb)=-tTE/2.0d0
    zmodesb_use(2:nmodesb-1)=modesb
  
    ! Calculates the X, Z coordinates of the airfoil
    
    !call BSpline_C2A(b_spline_degree,nmodest,xmodest,zmodest_use,nPointst,ut_seed,xt_new,zt_new)
    call BSpline_C2A_b_matrix(nmodest,xmodest,zmodest_use,nPointst,top_b_matrix,xt_new,zt_new)
    
    if (.not. symmetrical) then
      !call BSpline_C2A(b_spline_degree,nmodesb,xmodesb,zmodesb_use,nPointsb,ub_seed,xb_new,zb_new)
      call BSpline_C2A_b_matrix(nmodesb,xmodesb,zmodesb_use,nPointsb,bot_b_matrix,xb_new,zb_new)
    else
      ! For symmetrical airfoils, just mirror the top surface
      xb_new = xt_new
      zb_new = -zt_new
    end if
    
  else    
    allocate(zmodest_use(nmodest),zmodesb_use(nmodesb))
    allocate(xmodest_use(nmodest),xmodesb_use(nmodesb))
    
    ! Construst full control points vectors
    xmodest_use(1)=0.0d0
    xmodest_use(2)=0.0d0
    xmodest_use(nmodest)=xt_seed(nPointst)
    xmodest_use(3:nmodest-1)=modest(1:nmodest-3)
  
    call sort_vector(xmodest_use)
    
    xmodesb_use(1)=0.0d0
    xmodesb_use(2)=0.0d0
    xmodesb_use(nmodesb)=xb_seed(nPointsb)
    xmodesb_use(3:nmodesb-1)=modesb(1:nmodesb-3)
    
    call sort_vector(xmodesb_use)
    
    zmodest_use(1)=0.0d0
    zmodest_use(nmodest)=+tTE/2.0d0
    zmodest_use(2:nmodest-1)=modest(nmodest-2:2*nmodest-5)
  
    zmodesb_use(1)=0.0d0
    zmodesb_use(nmodesb)=-tTE/2.0d0
    zmodesb_use(2:nmodesb-1)=modesb(nmodest-2:2*nmodest-5)
    
    ! get U
    
    call GET_U_NEWTON(b_spline_degree,nmodest,xmodest_use,nPointst,xt_seed,ut_use)
    call GET_U_NEWTON(b_spline_degree,nmodesb,xmodesb_use,nPointsb,xb_seed,ub_use)
    
    ! Calculates the X, Z coordinates of the airfoil
    
    call BSpline_C2A(b_spline_degree,nmodest,xmodest_use,zmodest_use,nPointst,ut_use,xt_new,zt_new)
    !call BSpline_C2A_b_matrix(nmodest,xmodest_use,zmodest_use,nPointst,top_b_matrix,xt_new,zt_new)
    
    if (.not. symmetrical) then
      call BSpline_C2A(b_spline_degree,nmodesb,xmodesb_use,zmodesb_use,nPointsb,ub_use,xb_new,zb_new)
      !call BSpline_C2A_b_matrix(nmodesb,xmodesb_use,zmodesb_use,nPointsb,bot_b_matrix,xb_new,zb_new)
    else
      ! For symmetrical airfoils, just mirror the top surface
      xb_new = xt_new
      zb_new = -zt_new
    end if
    
    ! Interpolate to x seed
  end if
  
end subroutine  BSP_airfoil

!***************************************************************
! Given n control points of a b-spline and parametric point vector u, it
! obtains the x, y coordinates of points distributed along the b-spline.
!***************************************************************
subroutine BSpline_C2A(n,n_cp,xin,yin,n_p,uspline,xspline,yspline)
  
  ! Modules.

  ! Variable declaration.
  implicit none

  ! Parameters.

  ! Input variables.
  integer, intent(in) ::                 n !degree of the spline
  integer, intent(in) ::                 n_cp          ! no. of control points
  integer, intent(in) ::                 n_p          ! no. of desired spline points
  real*8, dimension(n_cp), intent(in) :: xin        ! x-coordinate of control points
  real*8, dimension(n_cp), intent(in) :: yin        ! y-coordinate of control points
  real*8, dimension(n_p), intent(in)  :: uspline    ! u coordinate of spline

  ! Output variables.
  real*8, dimension(n_p), intent(out)  :: xspline    ! x coordinate of spline
  real*8, dimension(n_p), intent(out)  :: yspline    ! y coordinate of spline

  ! Local variables.
  integer ::i,j,k        ! counters
  real*8, dimension(n_cp) :: x! x-coordinate of control points
  real*8, dimension(n_cp) :: y! y-coordinate of control points
  real*8, dimension(n_p,2) ::spline! values of x,y on the spline
  
  real*8, dimension(n+1):: b! local spline values
  real*8, dimension(n_cp+n+1):: t! knot vector
  
  !Knot vector
  do i =1,n_cp+n+1
    if (i .LE. n+1) then
      t(i)=0.d0
    elseif (i .GT. n_cp) then 
      t(i)=1.d0
    else
      !uniform knot vector
      t(i)=t(i-1)+1/(dble(n_cp-n))
    end if
  end do
  
  ! Determine control points vector
  x(:)        = xin
  y(:)        = yin


  ! Compute the b-spline.
  spline(:,1)=0.d0
  spline(:,2)=0.d0
  spline(1,1)=x(1)
  spline(1,2)=y(1)
  spline(n_p,1)=x(n_cp)
  spline(n_p,2)=y(n_cp)
  
  j             = 2

  !for each interval
  do i=n+1,n_cp
    !for u in t(i)<u<t(i+1)
    do while(uspline(j).LT.t(i+1))
      call d_BSpline(t,n,i,uspline(j),b)
      do k=1,n+1
        spline(j,1)=spline(j,1)+b(k)*x(i-n-1+k)
        spline(j,2)=spline(j,2)+b(k)*y(i-n-1+k)
      end do

      j         = j+1
    end do
  end do
  xspline(:)    = spline(:,1)
  yspline(:)    = spline(:,2)
  
end subroutine BSpline_C2A

!***************************************************************
! Given n control points of a b-spline and the b matrix, it
! obtains the x, y coordinates of points distributed along the b-spline.
!***************************************************************
subroutine BSpline_C2A_b_matrix(n_cp,xin,yin,n_p,BMATRIX,xspline,yspline)
  
  ! Modules.

  ! Variable declaration.
  implicit none

  ! Parameters.

  ! Input variables.
  integer, intent(in) ::                 n_cp          ! no. of control points
  integer, intent(in) ::                 n_p          ! no. of desired spline points
  real*8, dimension(n_cp), intent(in) :: xin        ! x-coordinate of control points
  real*8, dimension(n_cp), intent(in) :: yin        ! y-coordinate of control points
  real*8, dimension(n_p,n_cp), intent(in)  :: BMATRIX    ! b_matrix

  ! Output variables.
  real*8, dimension(n_p), intent(out)  :: xspline    ! x coordinate of spline
  real*8, dimension(n_p), intent(out)  :: yspline    ! y coordinate of spline

  ! Local variables.
  real*8, dimension(n_cp) :: x! x-coordinate of control points
  real*8, dimension(n_cp) :: y! y-coordinate of control points
  
  ! Determine control points vector
  x(:)        = xin
  y(:)        = yin

  ! Compute the b-spline.
  
  xspline   = matmul(BMATRIX,x)
  yspline   = matmul(BMATRIX,y)
  
end subroutine BSpline_C2A_b_matrix

!******************************************************************
! Routine that implements De Boor's algorithm, present in:
! A Practical Guide to Splines, by Carl de Boor
! Equivalent to the BSPLVB subroutne 
!******************************************************************
subroutine d_BSpline(t,k,i,x,b)

  ! Modules.

  ! Variable declaration.
  implicit none

  ! Parameters.

  ! Input variables.
  real*8, dimension(:), intent(in) :: t !knot vector
  integer, intent(in) :: k !degree of spline
  integer, intent(in) :: i !knot position, t(i)<x<t(i+1)
  real*8, intent(in) :: x !position to evaluate spline

  ! Output variables.
  real*8, dimension(k+1),intent(out) :: b !vector with the splines value

  ! Local variables.
  real*8, dimension(k) ::delta_R, delta_L
  integer :: j, r
  real*8 :: saved, term
  
  ! Begining of the code.
  b(1)=1.d0
  do j=1,k
    delta_R(j)=t(i+j)-x
    delta_L(j)=x-t(i+1-j)
    saved=0.d0
    do r=1,j
      term=b(r)/(delta_R(r)+delta_L(j+1-r))
      b(r)=saved+delta_R(r)*term
      saved=delta_L(j+1-r)*term
    end do
    b(j+1)=saved
  end do
  
end subroutine d_BSpline  

! ----------------------------------------------------------------------------
! Subroutine that implements the B-Spline Parameterization, coordinates to weights.
subroutine BSP_init(xseedt, xseedb, zseedt, zseedb, modest, modesb)

  use vardef, only : nparams_top, nparams_bot,b_spline_degree, b_spline_xtype, &
                    b_spline_distribution, upointst, upointsb, xcontrolt,      &
                    xcontrolb
  
  implicit none
  
  real*8, dimension(:), intent (in) :: xseedt, xseedb, zseedt, zseedb
  real*8, dimension(:), allocatable, intent (out) :: modest
  real*8, dimension(:), allocatable, intent (out) :: modesb
  
  real*8, dimension(:), allocatable:: modest_full
  real*8, dimension(:), allocatable:: modesb_full

  integer:: npointst, npointsb
  
  npointst=size(xseedt,1)
  npointsb=size(xseedb,1)
  
  allocate(upointst(npointst), upointsb(npointsb))
  allocate(xcontrolt(nparams_top), xcontrolb(nparams_bot))
  
  xcontrolt=0.d0
  xcontrolb=0.d0
  
  if (b_spline_xtype .EQ. 1) then
    call allocate_b_matrix(nparams_top, nparams_bot, npointst, npointsb)
    allocate(modest(nparams_top-2), modesb(nparams_bot-2))
    allocate(modest_full(nparams_top), modesb_full(nparams_bot))
    
    call BSpline_A2C_fixedx(b_spline_degree,b_spline_distribution,npointst,    &
    xseedt,zseedt,upointst,nparams_top,xcontrolt,modest_full,top_b_matrix)
  
    modest=modest_full(2:nparams_top-1)
    
    call BSpline_A2C_fixedx(b_spline_degree,b_spline_distribution,npointsb,    &
    xseedb,zseedb,upointsb,nparams_bot,xcontrolb,modesb_full,bot_b_matrix)
    
    modesb=modesb_full(2:nparams_bot-1)
  else
    call allocate_b_matrix(nparams_top, nparams_bot, npointst, npointsb)
    allocate(modest(2*nparams_top-5), modesb(2*nparams_bot-5))
    allocate(modest_full(nparams_top*2), modesb_full(nparams_bot*2))

    call BSpline_A2C_fixedx(b_spline_degree,b_spline_distribution,npointst,    &
    xseedt,zseedt,upointst,nparams_top,modest_full(1:nparams_top),              &
    modest_full(1+nparams_top:2*nparams_top),top_b_matrix)

   !call BSpline_A2C_freex(b_spline_degree,b_spline_distribution,npointst,     &
   !xseedt,zseedt,upointst,nparams_top,modest_full(1:nparams_top),              &
   !modest_full(1+nparams_top:2*nparams_top))
    
    !x
    modest(1:nparams_top-3)=modest_full(3:nparams_top-1)
    !z
    modest(nparams_top-2:2*nparams_top-5)=modest_full(nparams_top+2:2*nparams_top-1)
    
    call BSpline_A2C_fixedx(b_spline_degree,b_spline_distribution,npointsb,    &
    xseedb,zseedb,upointsb,nparams_bot,modesb_full(1:nparams_bot),              &
    modesb_full(1+nparams_bot:2*nparams_bot),bot_b_matrix)

   !call BSpline_A2C_freex(b_spline_degree,b_spline_distribution,npointsb,     &
   !xseedb,zseedb,upointsb,nparams_bot,modesb_full(1:nparams_bot),              &
   !modesb_full(1+nparams_bot:2*nparams_bot))
    
    !x
    modesb(1:nparams_bot-3)=modesb_full(3:nparams_bot-1)
    !z
    modesb(nparams_bot-2:2*nparams_bot-5)=modesb_full(nparams_bot+2:2*nparams_bot-1)
  end if

end subroutine BSP_init

! -------------------------------------------------------------------
! Subroutine to determine control points from its  airfoil b-spline points.
! X position is known and fixed
! U is computed from X and xspline, newton method
! Z is computed from U and zspline, min squares
subroutine BSpline_A2C_fixedx(n,option,nspline,xspline,zspline,uspline,n_cp,x, &
z,BMATRIX)
  
  ! Variables declaration.
  implicit none

  ! Input variables.
  integer, intent(in) :: n ! degree of spline
  integer, intent(in) :: option !type of x distribution
                             ! 1 = use cosine distribution	
									           ! 2 = use non- uniform cosine distribution
  integer, intent(in) :: n_cp ! number of control points
  integer, intent(in) :: nspline ! number of spline points
  real*8, intent(in) :: xspline(nspline) ! longitudinal spline point position
  real*8, intent(in) :: zspline(nspline) ! vertical spline point position

  ! Output variables.
  real*8, intent(out) :: uspline(nspline) !parametric value
  real*8, intent(out) :: x(n_cp)! x-coordinate of control points
  real*8, intent(out) :: z(n_cp)! z-coordinate of control points
  real*8, intent(out):: BMATRIX(nspline,n_cp)
  
  ! Get x control points
  x(1)=0.0d0
  call SetDistribution(option,n_cp-1,x(2:n_cp))

  !call SetDistribution(option,n_cp,x)

  ! Get u from x control points and x airfoil coordinates

  call GET_U_NEWTON(n,n_cp,x,nspline,xspline,uspline)

  call GET_Z_MINSQUARES(nspline,uspline,zspline,n,n_cp,z,BMATRIX)
  
end subroutine BSpline_A2C_fixedx  

! -------------------------------------------------------------------
! Subroutine to determine control points from its  airfoil b-spline points.
! U position is known and fixed
! X is computed from U and xspline, min squares
! Z is computed from U and zspline, min squares
subroutine BSpline_A2C_freex(n,option,nspline,xspline,zspline,uspline,n_cp,x,z)
  
  ! Variables declaration.
  implicit none

  ! Input variables.
  integer, intent(in) :: n ! degree of spline
  integer, intent(in) :: option !type of x distribution
                             ! 1 = use cosine distribution	
									           ! 2 = use non- uniform cosine distribution
  integer, intent(in) :: n_cp ! number of control points
  integer, intent(in) :: nspline ! number of spline points
  real*8, intent(in) :: xspline(nspline) ! longitudinal spline point position
  real*8, intent(in) :: zspline(nspline) ! vertical spline point position

  ! Output variables.
  real*8, intent(out) :: uspline(nspline) !parametric value
  real*8, intent(out) :: x(n_cp)! x-coordinate of control points
  real*8, intent(out) :: z(n_cp)! z-coordinate of control points
  
  !Local variables.
  real*8:: BSPLINE_FULL(nspline,n_cp)
  ! Set U position
  call SetDistribution(option,nspline,uspline)

  ! Get X, computed from U and xspline, min squares
  call GET_X_MINSQUARES(nspline,uspline,xspline,n,n_cp,x,BSPLINE_FULL)

  ! Get Z, computed from U and zspline, min squares
  call GET_Z_MINSQUARES(nspline,uspline,zspline,n,n_cp,z,BSPLINE_FULL)
  
end subroutine BSpline_A2C_freex  

!--------------------------------------------------------------------
! Get u from x control points and x airfoil coordinates
subroutine GET_U_NEWTON(n,n_cp,x,nspline,xspline,u)
  ! Variables declaration.
  implicit none

  ! Input variables.
  integer, intent(in) :: n !degree of spline
  integer, intent(in) :: n_cp! number of control points
  real*8, intent(in) :: x(n_cp)! x-coordinate of control points
  integer, intent(in) :: nspline! number of spline points
  real*8, intent(in) :: xspline(nspline)! longitudinal spline point position

  ! Output variables.
  real*8, intent(out):: u(nspline)!parametric value

  ! Local variables
  integer:: i,j
  real*8:: x_new,dx_new,u_dummy,beta
  
  ! Code

  do i=1,nspline
    j=1
    if (i .EQ. 1) then
      u(i)=0.0D0
    elseif (i .EQ. nspline) then
      u(i)=1.0D0
    else
      !Newton's method to get u
      u(i)=xspline(i)
      x_new=2.0D0 !just to fail condition
      ! limited to 1000 iterations, if convergence is impossible
      do while( ((abs(xspline(i)-x_new)/abs(xspline(i))) .GT. 1.0e-10) .and. (j .LT. 1000) )
        !get new x
        call B_Spline_Single(n,n_cp,x,u(i),x_new)
        !get derivative of new x
        call B_Spline_Single_d(n,n_cp,x,u(i),dx_new)
        !update u
        u_dummy=u(i)-(x_new-xspline(i))/dx_new !new value
        !check limits
        if(u_dummy .LT. 0.d0) u_dummy=0.d0!+1.0e-12
        if(u_dummy .GT. 1.d0) u_dummy=1.d0!-1.0e-12
        !mann sheme (Iterative Algorithms)
        beta=(1.d0/dble(j))**0.5 !just a divergent series
        u(i)=(dble(1)-beta)*u(i)+(beta)*u_dummy
        j=j+1
      end do
    end if
  end do

end subroutine GET_U_NEWTON  

!--------------------------------------------------------------------
! Get p control points from u parametric vector and p airfoil coordinates
subroutine GET_Z_MINSQUARES(nspline,uspline,pspline,n,n_cp,pcontrol,BSPLINE_FULL)
  
  ! Variables declaration.
  implicit none

  ! Input variables.
  integer, intent(in) :: n_cp! number of control points
  integer, intent(in) :: n !degree of spline
  integer, intent(in) :: nspline! number of spline points
  real*8, intent(in) :: pspline(nspline)! p spline point position
  real*8, intent(in):: uspline(nspline)!parametric spline vector

  ! Output variables.
  real*8, intent(out) :: pcontrol(n_cp)! p-coordinate of control points
  real*8, intent(out):: BSPLINE_FULL(nspline,n_cp)

  ! Local variables
  integer:: i,j
  real*8:: P_FULL(nspline)!, BSPLINE_FULL(nspline,n_cp)
  real*8:: BSPLINE_T(n_cp-2,nspline-2), BSPLINE(nspline-2,n_cp-2)
  real*8:: B_DUMMY(n+1)
  real*8:: p(nspline-2)
  real*8:: A(n_cp-2,n_cp-2), B(n_cp-2)
  real*8, dimension(n_cp+n+1):: t 
  
  double precision, dimension(n_cp-2,1) ::           B_send
  double precision, dimension(n_cp-2) ::             S
  double precision ::                                RCOND
  integer ::                                         RANK, LWORK, INFO
  double precision, dimension(5*(n_cp-2)) ::         WORK
  
  !Knot vector
  do i =1,n_cp+n+1
    if (i .LE. n+1) then
      t(i)=0.d0
    elseif (i .GT. n_cp) then
      t(i)=1.d0
    else
      !uniform knot vector
      t(i)=t(i-1)+1/(dble(n_cp-n))
    end if
  end do
  
  ! Assign known points
  pcontrol=0.d0
  pcontrol(1)=pspline(1)
  pcontrol(n_cp)=pspline(nspline)
  
  ! Create BSPLINE
  BSPLINE_FULL=0.d0
  j=1
  do i=1,nspline
    if(.NOT. uspline(i).LE.t(j+n+1)) j=j+1
    call d_BSpline(t,n,j+n,uspline(i),B_DUMMY)
    BSPLINE_FULL(i,j:j+n)=B_DUMMY
  end do
  
  !Create other matrixes
  BSPLINE=BSPLINE_FULL(2:nspline-1,2:n_cp-1)
  BSPLINE_T=transpose(BSPLINE)
  A=matmul(BSPLINE_T,BSPLINE)

  P_FULL=pspline
  p=P_FULL(2:nspline-1)-pcontrol(1)*BSPLINE_FULL(2:nspline-1,1) &
                       -pcontrol(n_cp)*BSPLINE_FULL(2:nspline-1,n_cp)
  B=matmul(BSPLINE_T,p)
  
  !Solve A.X=B
  B_send(:,1) = B(:)
  RCOND = 1.0e-6
  LWORK = 5*(n_cp-2)
  call DGELSS(n_cp-2, n_cp-2, 1, A, n_cp-2, B_send, n_cp-2,          &
    S, RCOND, RANK, WORK, LWORK, INFO )
  
  pcontrol(2:n_cp-1)  = B_send(:,1)
  
end subroutine GET_Z_MINSQUARES

!--------------------------------------------------------------------
! Get p control points from u parametric vector and p airfoil coordinates
subroutine GET_X_MINSQUARES(nspline,uspline,pspline,n,n_cp,pcontrol,BSPLINE_FULL)
  
  ! Variables declaration.
  implicit none

  ! Input variables.
  integer, intent(in) :: n_cp! number of control points
  integer, intent(in) :: n !degree of spline
  integer, intent(in) :: nspline! number of spline points
  real*8, intent(in) :: pspline(nspline)! p spline point position
  real*8, intent(in):: uspline(nspline)!parametric spline vector

  ! Output variables.
  real*8, intent(out) :: pcontrol(n_cp)! p-coordinate of control points
  real*8, intent(out):: BSPLINE_FULL(nspline,n_cp)

  ! Local variables
  integer:: i,j
  real*8:: P_FULL(nspline)!, BSPLINE_FULL(nspline,n_cp)
  real*8:: BSPLINE_T(n_cp-3,nspline-3), BSPLINE(nspline-3,n_cp-3)
  real*8:: B_DUMMY(n+1)
  real*8:: p(nspline-3)
  real*8:: A(n_cp-3,n_cp-3), B(n_cp-3)
  real*8, dimension(n_cp+n+1):: t 
  
  double precision, dimension(n_cp-3,1) ::           B_send
  double precision, dimension(n_cp-3) ::             S
  double precision ::                                RCOND
  integer ::                                         RANK, LWORK, INFO
  double precision, dimension(5*(n_cp-3)) ::         WORK
  
  !Knot vector
  do i =1,n_cp+n+1
    if (i .LE. n+1) then
      t(i)=0.d0
    elseif (i .GT. n_cp) then
      t(i)=1.d0
    else
      !uniform knot vector
      t(i)=t(i-1)+1/(dble(n_cp-n))
    end if
  end do
  
  ! Assign known points
  pcontrol=0.d0
  pcontrol(1)=0.d0
  pcontrol(2)=0.d0
  pcontrol(n_cp)=pspline(nspline)
  
  write(*,*) pcontrol(1), pcontrol(n_cp)
  
  ! Create BSPLINE
  BSPLINE_FULL=0.d0
  j=1
  do i=1,nspline
    if(.NOT. uspline(i).LE.t(j+n+1)) j=j+1
    call d_BSpline(t,n,j+n,uspline(i),B_DUMMY)
    BSPLINE_FULL(i,j:j+n)=B_DUMMY
  end do
  
  !Create other matrixes
  BSPLINE=BSPLINE_FULL(3:nspline-1,3:n_cp-1)
  BSPLINE_T=transpose(BSPLINE)
  A=matmul(BSPLINE_T,BSPLINE)

  P_FULL=pspline
  p=P_FULL(3:nspline-1)-pcontrol(1)*BSPLINE_FULL(3:nspline-1,1) &
                       -pcontrol(2)*BSPLINE_FULL(3:nspline-1,2) &
                       -pcontrol(n_cp)*BSPLINE_FULL(3:nspline-1,n_cp)
  B=matmul(BSPLINE_T,p)

  !Solve A.X=B
  B_send(:,1) = B(:)
  RCOND = 1.0e-6
  LWORK = 5*(n_cp-3)
  call DGELSS(n_cp-3, n_cp-3, 1, A, n_cp-3, B_send, n_cp-3,          &
    S, RCOND, RANK, WORK, LWORK, INFO )
  
  pcontrol(3:n_cp-1)  = B_send(:,1)
  
end subroutine GET_X_MINSQUARES

!***************************************************************
! Given n control points of a b-spline and the parametric point u, it
! obtains the x coordinate of a single point distributed along the b-spline.
!***************************************************************
subroutine B_Spline_Single(n,n_cp,xin,uspline,xspline)
  
  ! Modules.

  ! Variable declaration.
  implicit none

  ! Parameters.

  ! Input variables.
  integer, intent(in) ::                 n !degree of the spline
  integer, intent(in) ::                 n_cp          ! no. of control points
  real*8, dimension(n_cp), intent(in) :: xin        ! x-coordinate of control points
  real*8, intent(in)  :: uspline    ! u coordinate of spline

  ! Output variables.
  real*8, intent(out)::xspline! x coordinate of spline

  ! Local variables.
  integer ::i,k! counters
  real*8, dimension(n_cp) :: x! x-coordinate of control points
  real*8 ::spline! value of x on the spline
  
  real*8, dimension(n+1):: b
  real*8, dimension(n_cp+n+1):: t 
  
  !Knot vector
  do i =1,n_cp+n+1
    if (i .LE. n+1) then
      t(i)=0.d0
    elseif (i .GT. n_cp) then
      t(i)=1.d0
    else
      !uniform knot vector
      t(i)=t(i-1)+1/(dble(n_cp-n))
    end if
  end do

  ! Determine control points vector
  x(:)        = xin

  ! Compute the b-spline
  spline=0.d0
  mainloop: do i=n+1,n_cp
    if((t(i).LE.uspline).AND.(uspline.LE.t(i+1))) then
      call d_BSpline(t,n,i,uspline,b)
      do k=1,n+1
        spline=spline+b(k)*x(i-n-1+k)
      end do
      exit mainloop
    end if
  end do mainloop
  
  xspline    = spline
  
  end subroutine B_Spline_Single
  
!***************************************************************
! Given n control points of a b-spline and the parametric point u, it obtains
! the x coordinate of a single point distributed along the b-spline derivative.
!***************************************************************
subroutine B_Spline_Single_d(n,n_cp,xin,uspline,xspline)
  
  ! Modules.

  ! Variable declaration.
  implicit none

  ! Parameters.

  ! Input variables.
  integer, intent(in) ::n!degree of the spline
  integer, intent(in) ::n_cp! no. of control points
  real*8, dimension(n_cp), intent(in) :: xin! x-coordinate of control points
  real*8, intent(in)  :: uspline    ! u coordinate of spline

  ! Output variables.
  real*8, intent(out)::xspline! x coordinate of spline

  ! Local variables.
  integer ::i,k! counters
  real*8, dimension(n_cp) :: x! x-coordinate of control points
  real*8 ::spline! values of x on the spline
  
  real*8, dimension(n+1):: b
  real*8, dimension(n_cp+n+1-2):: t 
  
  !Knot vector
  do i =1,n_cp+n+1-2
    if (i .LE. n) then
      t(i)=0.d0
    elseif (i .GT. n_cp-1) then
      t(i)=1.d0
    else
      !uniform knot vector
      t(i)=t(i-1)+1/(dble(n_cp-n))
    end if
  end do

  ! Determine control points vector
  x(:)        = xin

  ! Compute the b-spline.
  spline=0.d0
  mainloop: do i=n,n_cp-1
    !for u in t(i)<u<t(i+1)
    if((t(i).LE.uspline).AND.(uspline.LE.t(i+1))) then
      call d_BSpline(t,n-1,i,uspline,b)
      do k=1,n
        spline=spline+b(k)*(x(i-n+1+k)-x(i-n+k))*(n-1)
      end do
      exit mainloop
    end if
  end do mainloop
  xspline    = spline
end subroutine B_Spline_Single_d  

! -------------------------------------------------------------------
! Determines new x-coordinate distribution on airfoil chord.
subroutine SetDistribution(option,n,x)

  ! Modules.
  
  ! Variables declaration.
  implicit none

  ! Parameters.
!  real*8, parameter ::                    small = 0.00001     ! small real*8

  ! Input variables.
  integer, intent(in) ::option! 1 = use cosine distribution	
! 2 = use non- uniform cosine distribution
! 3 = use cosine - LE & constant TE distribution
  integer, intent(in) ::n! number of points

  ! Input/output variables.
  real*8, intent(out) ::x(n)! chordwise position of points

  ! Output variables.

  ! Local variables.
  integer ::i! index
  real*8 ::beta! angle interval to determine x
  real*8 ::Pi
  real*8 ::a! progression constant
  real*8 ::teta(n)! angle at x
  real*8 ::sum! sum of powers of a

  ! Chordwise position.
  Pi              = 4.0D0*datan(1.0D0)

  select case(option)
  case(1)
	  ! Cosine distribution.
      beta            = Pi/(n-1)
      do i=1,n
        x(i)          = 0.5D0*(1.0D0-dcos((i-1)*beta))
      end do
  case(2)
    ! Non-uniform cosine distribution.
    a             = 1.0d0+1.0d0/dble(n-1)
    sum           = 0.0d0
    do i=2,n
      sum         = sum+a**(i-2)
    end do
    beta          = Pi/sum
    teta(1)       = 0.0d0
    x(1)          = 0.0d0
    do i=2,n
      teta(i)     = teta(i-1)+a**(i-2)*beta
      x(i)        = 0.5D0*(1.0D0-dcos(teta(i)))
    end do
    x(n)          = 1.0D0
  case(3)
    ! Linear
    do i=1,n
      x(i)          = (dble(i-1))/dble(n-1)
    end do
  end select
!      do i=1,n
!        write(*,'(I6,2F12.6)') i,x(i),x(i)-x(max(1,i-1))
!      end do
!      read(*,*)

  return

end subroutine SetDistribution

!=============================================================================80
!
! Allocates memory for b_matrix
!
!=============================================================================80
subroutine allocate_b_matrix(nmodest, nmodesb, npointst, npointsb)

  integer, intent(in) :: nmodest, nmodesb, npointst, npointsb
  
  if (.not. allocated(top_b_matrix)) then
    allocate(top_b_matrix(npointst,nmodest+2))
  end if
  
  if (.not. allocated(bot_b_matrix)) then
    allocate(bot_b_matrix(npointsb,nmodesb+2))
  end if

  !   Initialize shape functions

  top_b_matrix(:,:) = 0.d0
  bot_b_matrix(:,:) = 0.d0
end subroutine allocate_b_matrix

!=============================================================================80
!
! Deallocates memory for b_matrix
!
!=============================================================================80
subroutine deallocate_b_matrix


  if (allocated(top_b_matrix)) then
    deallocate(top_b_matrix)
  end if
  
  if (allocated(bot_b_matrix)) then
    deallocate(bot_b_matrix)
  end if
  
end subroutine deallocate_b_matrix

!/////////////////////////////////////////////////////////////////////////////80
!
! BPP subroutines
!
!/////////////////////////////////////////////////////////////////////////////80

! ----------------------------------------------------------------------------
! Subroutine that implements the Bezier-PARSEC Parameterization, weights to coordinates.
subroutine BPP_airfoil( xPt, xPb, modest, modesb, zPt_new, zPb_new, tTE)

  implicit none

  real*8, dimension(:), intent (in) :: modest, modesb

  real*8, intent(in) ::    xPt(:),xPb(:)                     ! x values of data points
  real*8, intent(in) ::    tTE
  real*8, intent(out) ::   zPt_new(size(xPt,1)), zPb_new(size(xPb,1)) ! z values of data points
  
  integer ::               nPt, nPb                          ! no. of data points
  
  real*8 ::                Xt_max, Yt_max, Kt_max            ! Position and magnitude of the maximum thickness
  real*8 ::                Xc_max, Yc_max, Kc_max            ! Position and magnitude of the maximum camber
  real*8 ::                Rle                               ! Leading edge radius
  real*8 ::                Gamma_le, Alpha_te                ! Leading and trailing edge angle
  real*8 ::                Beta_te                           ! Trailing edge wedge angle
  real*8 ::                Zte, dZte                         ! Trailing edge vertical coordinate and thickness

  real*8, dimension(4) ::  Xt_le, Xt_te, Yt_le, Yt_te          ! control points for the thickness bezier curves
  real*8, dimension(4) ::  Xc_le, Xc_te, Yc_le, Yc_te          ! control points for the camber bezier curves
  
  integer :: error_code_t, error_code_c
  
  nPt=size(xPt,1)
  nPb=size(xPb,1)
  
  Xt_max=modest(1)
  Yt_max=modest(2)
  Kt_max=modest(3)
  Xc_max=modesb(1)
  Yc_max=modesb(2)
  Kc_max=modesb(3)            
  Rle=modest(4)                              
  Gamma_le=modesb(4)
  Alpha_te=modesb(5)                
  Beta_te=modest(5)                           
  Zte = 0.0d0
  dZte = tTE
  
  ! Calculating the control points of the bezier curves.
  call SetThicknessControlPoints(Xt_le, Xt_te, Yt_le, Yt_te,                   &
                                 Xt_max, Yt_max, Kt_max, Rle, Beta_te, dZte,   &
                                 error_code_t)

  call SetCamberControlPoints(Xc_le, Xc_te, Yc_le, Yc_te,                      &
                              Xc_max, Yc_max, Kc_max, Gamma_le, Alpha_te, Zte, &
                              error_code_c)

  ! Calculate the Z values top
  call BPP_Get_z(Xt_le,Xt_te,Yt_le,Yt_te,Xc_le,Xc_te,Yc_le,Yc_te,xPt,zPt_new,nPt,0)
  
  ! Calculate the Z values bot
  call BPP_Get_z(Xt_le,Xt_te,Yt_le,Yt_te,Xc_le,Xc_te,Yc_le,Yc_te,xPb,zPb_new,nPb,1)
  
  if (error_code_t /= 0 .OR. error_code_c /= 0) then
    zPt_new=0.d0
    zPb_new=0.d0
  end if
    
end subroutine BPP_airfoil

! -------------------------------------------------------------------
subroutine SetThicknessControlPoints(Xt_le, Xt_te, Yt_le, Yt_te,               &
                                    Xt_max, Yt_max, Kt_max, Rle, Beta_te, dZte,&
                                    error_code)

  implicit none

  real*8, intent(in) :: Xt_max, Yt_max, Kt_max ! Position and magnitude of the 
                                               ! maximum thickness
  real*8, intent(in) :: Rle                    ! Leading edge radius
  real*8, intent(in) :: Beta_te                ! Trailing edge wedge angle
  real*8, intent(in) :: dZte                   ! Trailing edge thickness
  real*8 :: Rt                                 ! Thickness auxiliar parameter

  real*8, dimension(4), intent(inout) :: Xt_le, Xt_te, Yt_le, Yt_te  ! Control 
                                        ! points for the thickness bezier curves
  integer, intent(out) :: error_code
  integer :: error
  
  error = 0
  
  ! Computes the Rt auxiliar parameter.
  call Rt_calc(Rt, Kt_max, Xt_max, Yt_max, Rle, error_code)

  ! Set the values of the leading edge Bezier curve control points.
  Xt_le(1)=0.0d0
  Xt_le(2)=0.0d0
  Xt_le(3)=Rt
  Xt_le(4)=Xt_max

  !! Check the validity of the control points' values.
  call BPPCheckControlPoints(Xt_le,0.0d0,Xt_max,error)

  if (error /= 0) error_code=error
  
  Yt_le(1)=0.0d0
  Yt_le(2)=3.0d0*Kt_max*(Xt_max-Rt)**2/2.0d0+Yt_max
  Yt_le(3)=Yt_max
  Yt_le(4)=Yt_max

  call BPPCheckControlPoints(Yt_le,0.0d0,Yt_max,error)

  if (error /= 0) error_code=error
  
  ! Set the values of the trailing edge Bezier curve control points.
  Xt_te(1)=Xt_max
  Xt_te(2)=2.0d0*Xt_max-Rt
  Xt_te(3)=1.0d0+(dZte-(3.0d0*Kt_max*(Xt_max-Rt)**2./2.0d0+Yt_max)) *          &
                                                              1.0d0/tan(Beta_te)
  Xt_te(4)=1.0d0

  call BPPCheckControlPoints(Xt_te,Xt_max,1.0d0,error)

  if (error /= 0) error_code=error
  
  Yt_te(1)=Yt_max
  Yt_te(2)=Yt_max
  Yt_te(3)=3.0d0*Kt_max*(Xt_max-Rt)**2/2.0d0+Yt_max
  Yt_te(4)=dZte

  call BPPCheckControlPoints(Yt_te,Yt_max,0.0d0,error)

  if (error /= 0) error_code=error
  
end subroutine SetThicknessControlPoints

! -------------------------------------------------------------------
subroutine SetCamberControlPoints(Xc_le, Xc_te, Yc_le, Yc_te,                  &
                              Xc_max, Yc_max, Kc_max, Gamma_le, Alpha_te, Zte, &
                              error_code)
  implicit none

  real*8, intent(in) :: Xc_max, Yc_max, Kc_max ! Position and magnitude of the 
                                             ! maximum camber
  real*8, intent(in) :: Gamma_le, Alpha_te   ! Leading and trailing edge angle
  real*8, intent(in) :: Zte                  ! Trailing edge vertical coordinate
  real*8 :: Rc                               ! Camber auxiliar parameter

  real*8, dimension(4), intent(inout) :: Xc_le, Xc_te, Yc_le, Yc_te  ! Control
                                          ! points for the camber bezier curves
  integer, intent(out) :: error_code
  integer :: error
  
  if (Yc_max.LE.0.0d0) then
      Rc  = 0.0d0
    else
      call Rc_calc(Rc,Kc_max,Gamma_le,Alpha_te,Zte,Yc_max,error_code)
  end if

  ! Set the values of the leading edge Bezier curve control points.
  Xc_le(1)=0.0d0
  Xc_le(2)=Rc/tan(gamma_le)
  Xc_le(3)=Xc_max-(2.0d0*(Rc-Yc_max)/(3.0d0*Kc_max))**0.5
  Xc_le(4)=Xc_max

  call BPPCheckControlPoints(Xc_le,0.0d0,Xc_max,error)

  if (error /= 0) error_code=error
  
  Yc_le(1)=0.0d0
  Yc_le(2)=Rc
  Yc_le(3)=Yc_max
  Yc_le(4)=Yc_max

  call BPPCheckControlPoints(Yc_le,0.0d0,Yc_max,error)

  if (error /= 0) error_code=error
  
  ! Set the values of the trailing edge Bezier curve control points.
  Xc_te(1)=Xc_max
  Xc_te(2)=Xc_max+(2.0d0*(Rc-Yc_max)/(3.0d0*Kc_max))**0.5
  Xc_te(3)=1.0d0+(Zte-Rc)/(tan(Alpha_te))
  Xc_te(4)=1.0d0

  call BPPCheckControlPoints(Xc_te,Xc_max,1.0d0,error)

  if (error /= 0) error_code=error
  
  Yc_te(1)=Yc_max
  Yc_te(2)=Yc_max
  Yc_te(3)=Rc
  Yc_te(4)=Zte

  call BPPCheckControlPoints(Yc_te,Yc_max,0.0d0,error)

  if (error /= 0) error_code=error
  
end subroutine SetCamberControlPoints

! -------------------------------------------------------------------
subroutine Rt_calc(Rt, Kt_max, Xt_max, Yt_max, Rle, error_code)

  use PolynomialRoots

  implicit none

  real*8, intent(in) :: Kt_max, Xt_max, Yt_max, Rle
  real*8, intent(out) :: Rt
  integer, intent(out) :: error_code
  real*8, dimension(5) :: a                   ! Quartic polynomial coefficients in the form a(5)*x^4+a(4)*x^3+a(3)*x^2+a(2)*x+a(1)
  complex*16, dimension(4) :: X_root          ! Roots of the polynomial

  integer :: i

  error_code=0
  
  ! Set the parameters for the solution of the quartic equation
  a(5)=27.0d0*(Kt_max**2)/4.0d0
  a(4)=-27.0d0*(Kt_max**2)*Xt_max
  a(3)=9.0d0*Kt_max*Yt_max+(81.0d0*(Kt_max**2)*(Xt_max**2)/2.0d0)
  a(2)=2.0d0*Rle-18.0d0*Kt_max*Xt_max*Yt_max-27.0d0*(Kt_max**2)*(Xt_max**3)
  a(1)=3.0d0*(Yt_max**2)+9.0d0*Kt_max*(Xt_max**2)*Yt_max +                     &
                                          (27.0d0*(Kt_max**2)*(Xt_max**4)/4.0d0)

  ! Get the roots of the polynomial
  call QuarticRoots(a,X_root)

  ! Select the correct root
  Rt=2.0d0*Xt_max  ! Presets the solution

  do i=1,4
    if(AImag(X_root(i)).NE.0) cycle                                                    ! Discard the complex roots
    if(Real(X_root(i)).GT.max(0.d0,Xt_max-(-2.0d0*Yt_max/(3.0d0*Kt_max))**0.5)) then   ! Discard real roots smaller than the lower boundary
      if(Real(X_root(i)).LT.Rt) then                                                   ! Make sure the used root is the smallest within the boundaries
        Rt=Real(X_root(i))
      end if
    end if
  end do
  
  
  if(Rt.GE.Xt_max) then
    error_code=1
    !Rt=huge(1)
    Rt=Xt_max*Rand()
  end if

end subroutine Rt_calc

! -------------------------------------------------------------------
subroutine Rc_calc(Rc,Kc_max,Gamma_le,Alpha_te,Zte,Yc_max,error_code)

  implicit none

  real*8, intent(in) :: Kc_max, Gamma_le, Alpha_te, Zte, Yc_max
  integer, intent(out) :: error_code
  real*8, intent(out) :: Rc
  real*8 :: Aux1, Aux2, Aux3, Aux4                  !Auxiliar parameters

  error_code=0
  
  Aux4=(16.0d0+6.0d0*Kc_max*((1.0d0/tan(Gamma_le))+(1.0d0/tan(alpha_te))) *    &
    (1.0d0 - Yc_max*((1.0d0/tan(Gamma_le))+(1.0d0/tan(alpha_te))) +            &
                                                        Zte/tan(alpha_te)))**0.5
  
  Aux1=4*(16.0d0+6.0d0*Kc_max*((1.0d0/tan(Gamma_le))+(1.0d0/tan(alpha_te))) *  &
    (1.0d0 - Yc_max*((1.0d0/tan(Gamma_le))+(1.0d0/tan(alpha_te))) +            &
                                                        Zte/tan(alpha_te)))**0.5
  
  Aux2=(16.0d0+3.0d0*Kc_max*((1.0d0/tan(Gamma_le))+(1.0d0/tan(alpha_te))) *    &
                                                      (1.0d0+Zte/tan(alpha_te)))
  Aux3=3.0d0*Kc_max*((1.0d0/tan(Gamma_le))+(1.0d0/tan(alpha_te)))**2.0d0

  if ((Aux2+Aux1)/Aux3 .GT. 0 .AND. ((Aux2+Aux1)/Aux3).LT.Yc_max) then
      Rc=(Aux2+Aux1)/Aux3
    elseif (((Aux2-Aux1)/Aux3).GT.0.AND.((Aux2-Aux1)/Aux3).LT.Yc_max) then
      Rc=(Aux2-Aux1)/Aux3
    elseif (Yc_max.EQ.0) then
        Rc=0
    else
      Rc=Yc_max*Rand()
      !Rc=-1
      error_code=2
  end if

end subroutine Rc_calc

! -------------------------------------------------------------------
subroutine BPP_Get_z(Xt_le,Xt_te,Yt_le,Yt_te,Xc_le,Xc_te,Yc_le,Yc_te,          &
                     X,Z,nPoints,UpperLower_identifier)

  implicit none

  Real*8, dimension(4), intent(in) :: Xt_le, Xt_te, Yt_le, Yt_te,              &
                                      Xc_le, Xc_te, Yc_le, Yc_te
  integer, intent(in) :: nPoints
  integer, intent(in) :: UpperLower_identifier   ! Parameter that identifies the
                                      ! use of the upper, 0, or lower surface, 1
  Real*8, dimension(npoints), intent(in) :: X    ! Airfoil coordinates
  Real*8, dimension(npoints), intent(inout) :: Z ! Airfoil coordinates

  integer, dimension(2) :: LEte_identifier       ! Parameter that identifies the
                 ! use of the leading edge or trailing edge bezier curve for the 
                 ! (thickness,camber)

  integer :: i
  real*8, dimension(2) :: t                      ! Parameter for the calculation
                                ! of the bezir curves for the (thickness,camber)

  do i=1,nPoints      
    call BPP_Get_t(Xt_le,Xt_te,Xc_le,Xc_te,LEte_identifier,t,X(i))
    call BPP_Calc_z(Yt_le,Yt_te,Yc_le,Yc_te,LEte_identifier,t,Z(i),       &
      UpperLower_identifier)
  end do

end subroutine BPP_Get_z

! -------------------------------------------------------------------
! Returns the value of the parameter of the bezier functions
subroutine BPP_Get_t(Xt_le,Xt_te,Xc_le,Xc_te,LEte_identifier,t,X)

  use PolynomialRoots

  implicit none

  Real*8, dimension(4), intent(in) :: Xt_le, Xt_te, Xc_le, Xc_te
  integer, dimension(2), intent(inout) :: LEte_identifier
  real*8, dimension(2), intent(inout) :: t
  real*8, intent(in) :: X
  Complex*16, dimension(3) :: X_root ! cubic roots
  real*8, dimension(4) :: a          ! coeffitients to solve the cubic equation

  integer :: i

  !preset t to induce error
  t(1)=2.0d0
  t(2)=2.0d0

  ! thickness bezier curve

  a(1)=Xt_le(1)-X
  a(2)=3.0d0*Xt_le(2)-3.0d0*Xt_le(1)
  a(3)=3.0d0*Xt_le(3)-6.0d0*Xt_le(2)+3.0d0*Xt_le(1)
  a(4)=Xt_le(4)-3.0d0*Xt_le(3)+3.0d0*Xt_le(2)-Xt_le(1)

  call CubicRoots(a,X_root)

  do i=1,3
    if (aimag(X_root(i)).NE.0.0d0) cycle
    if (real(X_root(i)).GT.1.00001d0.OR.real(X_root(i)).LT.0.0d0) cycle
    if (real(X_root(i)).NE.real(X_root(i))) cycle
    t(1)=Real(X_root(i))
    LEte_identifier(1)=0       ! Identifies the leading edge to be used in the thickness section
  end do

  if (t(1).GT.1.0d0.OR.t(1).LT.0.0d0) then
    a(1)=Xt_te(1)-X
    a(2)=3.0d0*Xt_te(2)-3.0d0*Xt_te(1)
    a(3)=3.0d0*Xt_te(3)-6.0d0*Xt_te(2)+3.0d0*Xt_te(1)
    a(4)=Xt_te(4)-3.0d0*Xt_te(3)+3.0d0*Xt_te(2)-Xt_te(1)

    call CubicRoots(a,X_root)

    do i=1,3
        if (aimag(X_root(i)).NE.0.0d0) cycle
        if (real(X_root(i)).GT.1.00001d0.OR.real(X_root(i)).LT.0.0d0) cycle
        if (real(X_root(i)).NE.real(X_root(i))) cycle
        t(1)=Real(X_root(i))
        LEte_identifier(1)=1       ! Identifies the trailing edge to be used in the thickness section
    end do

  end if

  if (t(1).LT.0.0d0-1.E-6.OR.t(1).GT.1.0d0+1.E-6) LEte_identifier(1)=-1       ! Identifies error in the thickness section

  ! Camber Bezier curve

  a(1)=Xc_le(1)-X
  a(2)=3.0d0*Xc_le(2)-3.0d0*Xc_le(1)
  a(3)=3.0d0*Xc_le(3)-6.0d0*Xc_le(2)+3.0d0*Xc_le(1)
  a(4)=Xc_le(4)-3.0d0*Xc_le(3)+3.0d0*Xc_le(2)-Xc_le(1)

  call CubicRoots(a,X_root)

    do i=1,3
      if (aimag(X_root(i)).NE.0.0d0) cycle
      if (real(X_root(i)).GT.1.00001d0.OR.real(X_root(i)).LT.0.0d0) cycle
      if (real(X_root(i)).NE.real(X_root(i))) cycle
      t(2)=Real(X_root(i))
      LEte_identifier(2)=0       ! Identifies the leading edge to be used in the thickness section
    end do

  if (t(2).GT.1.0d0.OR.t(2).LT.0.0d0) then
      a(1)=Xc_te(1)-X
      a(2)=3.0d0*Xc_te(2)-3.0d0*Xc_te(1)
      a(3)=3.0d0*Xc_te(3)-6*Xc_te(2)+3.0d0*Xc_te(1)
      a(4)=Xc_te(4)-3.0d0*Xc_te(3)+3.0d0*Xc_te(2)-Xc_te(1)

      call CubicRoots(a,X_root)

      do i=1,3
          if (aimag(X_root(i)).NE.0.0d0) cycle
          if (real(X_root(i)).GT.1.00001d0.OR.real(X_root(i)).LT.0.0d0) cycle
          if (real(X_root(i)).NE.real(X_root(i))) cycle
          t(2)=Real(X_root(i))
          LEte_identifier(2)=1       ! Identifies the trailing edge to be used in the thickness section
      end do

  end if

  if ((t(2).LT.0.0d0-1.E-6).OR.(t(2).GT.1.0d0+1.E-6)) LEte_identifier(2)=-1       ! Identifies error in the thickness section

end subroutine BPP_Get_t

! -------------------------------------------------------------------
subroutine BPP_calc_z(Yt_le,Yt_te,Yc_le,Yc_te,                                 &
                      LEte_identifier,t,Z,UpperLower_identifier)

  implicit none

  Real*8, dimension(4), intent(in) :: Yt_le, Yt_te, Yc_le, Yc_te
  integer, dimension(2), intent(in) :: LEte_identifier
  real*8, dimension(2), intent(in) :: t
  real*8, intent(inout) :: Z
  integer, intent(in) :: UpperLower_identifier
  real*8 :: Z_thickness, Z_camber                     !Partial values of Z

  if (LEte_identifier(1).EQ.0) then
      Z_thickness=((1-t(1))**3)*Yt_le(1)+                                      &
                  3*((1-t(1))**2)*t(1)*Yt_le(2)+                               &
                  3*(1-t(1))*(t(1)**2)*Yt_le(3)+                               &
                  (t(1)**3)*Yt_le(4)
  elseif (LEte_identifier(1).EQ.1) then
      Z_thickness=((1-t(1))**3)*Yt_te(1)+                                      &
                  3*((1-t(1))**2)*t(1)*Yt_te(2)+                               &
                  3*(1-t(1))*(t(1)**2)*Yt_te(3)+                               &
                  (t(1)**3)*Yt_te(4)
  end if

  if (LEte_identifier(2).EQ.0) then
    Z_camber=((1-t(2))**3)*Yc_le(1)+                                           &
             3*((1-t(2))**2)*t(2)*Yc_le(2)+                                    &
             3*(1-t(2))*(t(2)**2)*Yc_le(3)+                                    &
             (t(2)**3)*Yc_le(4)
  elseif (LEte_identifier(2).EQ.1) then
    Z_camber=((1-t(2))**3)*Yc_te(1)+                                           &
             3*((1-t(2))**2)*t(2)*Yc_te(2)+                                    &
             3*(1-t(2))*(t(2)**2)*Yc_te(3)+                                    &
             (t(2)**3)*Yc_te(4)
  end if

  if (UpperLower_identifier.EQ.0) then
    Z=Z_camber+Z_thickness/2.0d0
  else
    Z=Z_camber-Z_thickness/2.0d0
  end if

end subroutine BPP_calc_z

! -------------------------------------------------------------------
subroutine BPPCheckControlPoints(Vec,Min_ini,Max_ini,error)

  implicit none

  Real*8, dimension(4), intent(inout) :: Vec      ! Input vector with the control points variables
  real*8, intent(in) :: Min_ini, Max_ini          ! Minimum and maximum of the interval
  integer, intent(out) :: error
  
  Real*8 :: Min, Max
  integer :: i
  
  error=0
  
  if (Min_ini.NE.Min_ini) then
    Min=0.5d0
    error=1
  else
    Min=Min_ini
  end if

  if (Max_ini.NE.Max_ini) then
    error=1  
    Max=0.5d0
  else
    Max=Max_ini
  end if

  do i=1,4
    if (abs(Vec(i)).GT.1.0d0.OR.Vec(i).NE.Vec(i)) then      ! Vec(i).NE.Vec(i) checks for NaN
      error=1
      Vec(i)=Min+(Max-Min)*(i-1)/3
    end if
  end do

end subroutine BPPCheckControlPoints

! ----------------------------------------------------------------------------
! Subroutine that implements the Bezier-PARSEC Parameterization, coordinates to weights.
subroutine BPP_init(xseedt, xseedb, zseedt, zseedb, modest, modesb)

  use math_deps, only : nu_first_derivative, nu_curvature

  real*8, intent(in) :: xseedt(:), xseedb(:), zseedt(:), zseedb(:)
  
  real*8, intent(out) :: modest(5) ! Xt_max, Yt_max, Kt_max, Rle     , Beta_te
  real*8, intent(out) :: modesb(5) ! Xc_max, Yc_max, Kc_max, Gamma_le, Alpha_te

  real*8:: xBPP(12)               ! Xt_max,Yt_max,Kt_max, Xc_max,Yc_max,Kc_max,&
                                  ! Rle, Gamma_le,Alpha_te,Beta_te, Zte,dZte
  real*8, dimension(size(xseedt,1)) :: zseedb_interpolated
  integer :: nPt, nPb
  
  real*8, dimension(size(xseedt,1)) :: zthick, zcamber, curvaturethick,        &
    curvaturecamber, first_derivative_thick, first_derivative_camber 
  real*8, dimension(size(xseedt,1)+size(xseedb,1)-1) :: xspline, zspline
  integer :: i
  
  nPt=size(xseedt,1)
  nPb=size(xseedb,1)
  
  ! interpolate lower surface points to upper surface points  
  call d_SplineAirfoilInterpolation(nPt,xseedt,zseedt,nPb,xseedb,zseedb,zseedb_interpolated)
  
  ! calculate thickness and camber distributions
  
  zthick=zseedt-zseedb_interpolated
  zcamber=(zseedt+zseedb_interpolated)/2.0d0
  
  do i=1,nPt
    xspline(i)=xseedt(nPt+1-i)
    zspline(i)=zseedt(nPt+1-i)
  end do
  
  do i=2,nPb
    xspline(i-1+nPt)=xseedb(i)
    zspline(i-1+nPt)=zseedb(i)
  end do
    
  curvaturethick=nu_curvature(nPt, xseedt, zthick)
  curvaturecamber=nu_curvature(nPt, xseedt, zcamber)
  
  !curvaturethick=curv_spline(nPt, xseedt, zthick)
  !curvaturecamber=curv_spline(nPt, xseedt, zcamber)
  
  first_derivative_thick=nu_first_derivative(nPt, xseedt, zthick)
  first_derivative_camber=nu_first_derivative(nPt, xseedt, zcamber)
  
  ! initialize xBPP
  xBPP=0.d0
  
  ! get Xt_max,Yt_max,Kt_max

  do i=1,nPt
    if (xBPP(2) .LT. zthick(i)) then
      xBPP(1)=xseedt(i)
      xBPP(2)=zthick(i)
      xBPP(3)=curvaturethick(i)
    end if
  end do
  
  write(*,*) 'Xt_max', xBPP(1)
  write(*,*) 'Yt_max', xBPP(2)
  write(*,*) 'Kt_max', xBPP(3)
  
  ! get Xc_max,Yc_max,Kc_max

  do i=1,nPt
    if (xBPP(5) .LT. zcamber(i)) then
      xBPP(4)=xseedt(i)
      xBPP(5)=zcamber(i)
      xBPP(6)=curvaturecamber(i)
    end if
  end do
  
  write(*,*) 'Xc_max', xBPP(4)
  write(*,*) 'Yc_max', xBPP(5)
  write(*,*) 'Kc_max', xBPP(6)
  
  ! get Rle
  
  xBPP(7) = 1.d0/curvaturethick(1)
  write(*,*) 'Rle', xBPP(7)
  
  ! get Gamma_le,Alpha_te,Beta_te
  
  do i = 1, nPt
    if (xseedt(i) .LT. 0.001d0) then
      xBPP(8) = atan(first_derivative_camber(i))
    end if
  end do
  
  do i = 1, nPt
    if (xseedt(nPt-i+1) .GT. 0.98d0) then
      xBPP(9)=-atan(first_derivative_camber(nPt-i+1))
      xBPP(10)=-atan(first_derivative_thick(nPt-i+1))
    end if
  end do
  
  write(*,*) 'Gamma_le', xBPP(8)/acos(-1.0d0)*180.d0
  write(*,*) 'Alpha_te', xBPP(9)/acos(-1.0d0)*180.d0
  write(*,*) 'Beta_te',xBPP(10)/acos(-1.0d0)*180.d0
  
  ! get Zte,dZte
  
  xBPP(11)= zcamber(nPt)
  xBPP(12)= zthick(nPt)
  write(*,*) 'Zte', xBPP(11)
  write(*,*) 'dZte',xBPP(12)
  
  ! set output
  modest(1:3)=xBPP(1:3)
  modest(4)=xBPP(7)
  modest(5)=xBPP(10)
  !modest(6)=xBPP(12)
  
  modesb(1:3)=xBPP(4:6)
  modesb(4:5)=xBPP(8:9)

end subroutine BPP_init

 !**************************************************************************
! Subroutine to interpolate values in a set of x,y data, given a value x
!   and using spline interpolation. Note: Spline coefficients are computed.
!**************************************************************************
subroutine d_SplineAirfoilInterpolation(ntop,xtop,ytop,nbot,xbot,ybot,ybot_at_xtop)

  ! Modules.
!  use IMSL

  ! Variable declaration.
  implicit none

  ! Parameters.

  ! Input variables.
  integer, intent(in) ::                        ntop, nbot           ! no. of data points
  real*8, dimension(ntop), intent(in) ::       xtop,ytop     ! x,y coordinates
  real*8, dimension(nbot), intent(in) ::       xbot,ybot     ! x,y coordinates

  ! Output variables.
  real*8, dimension(ntop), intent(out) ::          ybot_at_xtop               ! y coordinate of interpolated point

  ! Local variables.
  integer ::                                    i
  integer :: N
  real*8, dimension(ntop+nbot-1) :: S, X, XP, Y, YP
  real*8 :: SLE, SOPP, SI
  logical :: SILENT_MODE

  interface
    double precision function SEVAL(SS,X,XS,S,N)
      integer, intent(in) :: N
      double precision, intent(in) :: SS
      double precision, dimension(N), intent(in) :: X, XS, S
    end function SEVAL
  end interface 
  
  SILENT_MODE=.true.
  N = ntop+nbot-1
  do i=1,ntop
  X(i)=xtop(ntop+1-i)
  Y(i)=ytop(ntop+1-i)
  end do
  
  do i=2,nbot
  X(i-1+ntop)=xbot(i)
  Y(i-1+ntop)=ybot(i)
  end do

  CALL SCALC(X,Y,S,N)
  CALL SEGSPL(X,XP,S,N)
  CALL SEGSPL(Y,YP,S,N)
  CALL LEFIND(SLE,X,XP,Y,YP,S,N,SILENT_MODE)
  
  i=1
  SI=S(i)
  do while(S(i) .LE. SLE)
    call SOPPS(SOPP, SI, X,XP,Y,YP,S,N, SLE, SILENT_MODE)
    ybot_at_xtop(ntop+1-i) = SEVAL(SOPP,Y,YP,S,N)
    i=i+1
    SI=S(i)
  end do
  return
  
end subroutine d_SplineAirfoilInterpolation

function curv_spline(N, X, Y)

  integer, intent(in) :: N
  double precision, dimension(N), intent(in) :: X, Y
  double precision, dimension(N) :: curv_spline

 ! Local variables.
  integer :: I
  real*8, dimension(N) :: S, XP, YP
  real*8 :: SBREF

  interface
    double precision function CURV(SS,X,XS,Y,YS,S,N)
      integer, intent(in) :: N
      double precision, intent(in) :: SS
      double precision, dimension(N), intent(in) :: X, XS, Y,YS, S
    end function CURV
  end interface 
  
  CALL SCALC(X,Y,S,N)
  CALL SEGSPL(X,XP,S,N)
  CALL SEGSPL(Y,YP,S,N)

!---- normalizing length (~ chord)
  SBREF = (S(N)-S(1))

!---- set up curvature array
  DO I = 1, N
    curv_spline(I) = CURV(S(I),X,XP,Y,YP,S,N) * SBREF
  ENDDO
  
end function curv_spline

function der1_spline(N, X, Y)

  integer, intent(in) :: N
  double precision, dimension(N), intent(in) :: X, Y
  double precision, dimension(N) :: der1_spline

 ! Local variables.
  integer :: I
  real*8, dimension(N) :: S, XP, YP

  interface
    double precision function DEVAL(SS,X,XS,S,N)
      integer, intent(in) :: N
      double precision, intent(in) :: SS
      double precision, dimension(N), intent(in) :: X, XS, S
    end function DEVAL
  end interface 
  
  CALL SCALC(X,Y,S,N)
  CALL SEGSPL(X,XP,S,N)
  CALL SEGSPL(Y,YP,S,N)

!---- set up curvature array
  DO I = 1, N
    der1_spline(I) = DEVAL(S(I),Y,YP,S,N) / DEVAL(S(I),X,XP,S,N)
  ENDDO
  
end function der1_spline

end module parametrization_constr
  