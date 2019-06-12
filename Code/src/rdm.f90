module rdm

contains

! -----------------------------------------------------------------------------------
! initialize a random seed from the system clock at every run (fortran 95 code)
! -----------------------------------------------------------------------------------
subroutine init_random_seed
use glbl_prmtrs
integer :: i, n
double precision :: clock
integer, dimension(:), allocatable :: seed
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call random_seed(size = n) ! to get n, the size needed by the proc to initiate the seed
allocate(seed(n))
call cpu_time(clock)
seed = int(1.d9*clock) + 37 * (/ (i - 1, i = 1, n) /)
call random_seed(put = seed) ! now, set the seed to this random value
deallocate(seed)

end subroutine init_random_seed
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine get_flat(xmin,xmax,x)
use glbl_prmtrs
double precision, intent(in) :: xmin, xmax
double precision, intent(out) :: x
double precision :: zeta
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call random_number(zeta)
x=xmin+(xmax-xmin)*zeta

end subroutine get_flat
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine get_squared(x)
use glbl_prmtrs
double precision, intent(out) :: x
double precision :: zeta
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call random_number(zeta)
x=dsqrt(zeta)

end subroutine get_squared
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine get_exp(x)
use glbl_prmtrs
double precision, intent(out) :: x
double precision :: zeta
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call random_number(zeta)
x=-dlog(zeta)

end subroutine get_exp
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Von Neumann acceptance-rejectance method (aka rejection sampling) to pick up
! initial radial positions of clumps such as nvr2=cst w/ n # density of clumps
! and v a prescribed velocity law.
! N.B. : not sure that we need to normalize f here since we are only interested
! in knowing whether it is > or < ...
! -----------------------------------------------------------------------------------
subroutine accptnce_rjctnce(Ncl,r,R_cl,dens_cl)
use glbl_prmtrs
use mod_func
integer, intent(in) :: Ncl
double precision, intent(out) :: r(1:Ncl), R_cl(1:Ncl), dens_cl(1:Ncl)
integer :: i
double precision :: rmin, rmax, fmax, v, fnorm
! r and f picked up @ random
double precision :: rr, ff
! f evaluated @ rr picked up @ random
double precision :: fr
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! We want nvr2=cst ie (dN/4pir2dr) x vr2 = cst ie (dN/dr) as 1/v
! => the distribution function here is f(r) = 1 / ( v(r) x fnorm )

rmin = rini_
rmax = dist_max_cl_ ! we pick up clumps up to dist_max_cl

! Compute the integral of f over r=1 to rmax to normalize
call num_int_steps(100000,'one_over_v',rmin,rmax,fnorm,'linear')

fmax = one_over_v(rmin) / fnorm

! Pick up r between 1. and rmax and f between 0. and fmax
do i=1,Ncl
  fr=0.d0
  ff=1.d0 ! just to make sure we enter in the first do loop
  do while (ff>fr)
    call get_flat(rmin,rmax,rr)
    call get_flat(0.d0,fmax,ff)
    fr = one_over_v(rr) / fnorm
    ! call get_v(rr,v)
    ! fr = 1.d0 / ( v * rr**2.d0 * fnorm )
  enddo
  r(i)=rr
  ! Deduce clump radius & density
  if (rad_evol_=='lorenzo') then
    R_cl(i)=lorenzo_rad(rr)
    dens_cl(i)=lorenzo_rho(rr)
  elseif (rad_evol_=='jon') then
    R_cl(i)=jon_rad(rr)
    dens_cl(i)=jon_rho(rr)
  endif
enddo

end subroutine accptnce_rjctnce
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Numerical integration by step, with or without logarithmic stretching
! The user specifies a string 'which_func' which gives the name of a function defined
! in mod_func.f90.
! -----------------------------------------------------------------------------------
subroutine num_int_steps(Nstep,which_func,xmin,xmax,intgrl,logsteps)
use glbl_prmtrs
use mod_func
use miscellaneous

abstract interface
  function func (z)
     double precision :: func
     double precision, intent (in) :: z
  end function func
end interface

procedure (func), pointer :: f_ptr => null ()

integer, intent(in) :: Nstep
double precision, intent(in) :: xmin, xmax
character(LEN=*), intent(in) :: which_func
character(LEN=*), intent(in), optional :: logsteps
double precision, intent(out) :: intgrl

integer :: i
double precision :: dx, x, q
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if (xmax<xmin) then
  call crash('Beware, reverted bounds in integration.')
endif

if (logsteps=='log' .and. (xmin<0.d0 .or. xmax<0.d0)) &
  call crash('No log spacing for numerical integration possible if negative values in range.')

if (which_func=='one_over_vr2') then
  f_ptr => one_over_vr2
elseif (which_func=='one_over_v') then
  f_ptr => one_over_v
elseif (which_func=='r2_over_v') then
  f_ptr => r2_over_v
elseif (which_func=='r2_over_v_arcsin') then
  f_ptr => r2_over_v_arcsin
elseif (which_func=='r2v_23_over_v') then
  f_ptr => r2v_23_over_v
  elseif (which_func=='r2v_23_over_v_arcsin') then
    f_ptr => r2v_23_over_v_arcsin
else
  call crash('Unknown function. Define it in mod_func.f90 first.')
endif

intgrl = 0.d0

! Default is integration over linear steps
if ( (.not. present(logsteps)) .or. &
     (logsteps=='linear') ) then
  dx = (xmax-xmin)/dble(Nstep)
  x  = xmin+dx/2.d0
  do i=1,Nstep
    intgrl = intgrl + f_ptr(x) * dx
    x = x + dx
  enddo
elseif (logsteps=='log') then
  ! Stretching factor
  q   = (xmax/xmin)**(1.d0/dble(Nstep))
  dx  = xmin * (q-1.d0)
  x   = xmin + dx/2.d0 ! centers still half-way in-between the edges (reverse not true)
  do i=1,Nstep
    intgrl = intgrl + f_ptr(x) * dx
    x  =  x * q
    dx = dx * q
  enddo
else
  call crash('Unknown stepping in numerical integration (should be linear or log)')
endif

end subroutine num_int_steps
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine get_one_over_r2v(x)
use glbl_prmtrs
double precision, intent(out) :: x
double precision :: zeta, r, integral, x90, integral0
integer :: Nstep=1000 ! # of log steps before x90
double precision :: f0, f2, dx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

x90=1.d0/(1.d0-0.9d0**(1.d0/beta_)) ! radius to reach 90% of terminal speed

! Temporary fix, not satisfying, to compute "a" normalization x x x x x x
x=rini_
dx=2.d0*(rini_-1.d0)
call compute_deriv_nth_order(x,f0,f2)
integral0=f0*dx+(f2/2.d0)*(2.d0*(dx/2.d0)**3.d0)/3.d0
do while(x<100.d0)
  if (x<x90) then
    x=x*(x90/1.d0)**(1.d0/dble(Nstep))
    dx=dx*(x90/1.d0)**(1.d0/dble(Nstep))
  else
    x=x+dx
  endif
  call compute_deriv_nth_order(x,f0,f2)
  integral0=integral0+f0*dx+(f2/2.d0)*(2.d0*(dx/2.d0)**3.d0)/3.d0
enddo
! x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x x

call random_number(zeta)
! Since the function to integrate, 1/(r2v), is not analytically solvable,
! we rely on the numerical evaluation of the function using
! a local 2nd order expansion of the function
x=rini_ ! integral cell center
dx=2.d0*(rini_-1.d0) ! integral cell width
call compute_deriv_nth_order(x,f0,f2)
! The contribution of odd orders to the integral is zero
! since x is in-between its left and right edges
integral=f0*dx+(f2/2.d0)*(2.d0*(dx/2.d0)**3.d0)/3.d0
do while(integral/integral0<zeta)
  ! print*, x, zeta, integral, x<x90
  ! compute the integral
  if (x<x90) then
    x=x*(x90/1.d0)**(1.d0/dble(Nstep))
    dx=dx*(x90/1.d0)**(1.d0/dble(Nstep))
  else
    x=x+dx
  endif
  call compute_deriv_nth_order(x,f0,f2)
  integral=integral+f0*dx+(f2/2.d0)*(2.d0*(dx/2.d0)**3.d0)/3.d0
enddo

end subroutine get_one_over_r2v
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Compute 2nd order derivative of function f0=1/(r2v)
! -----------------------------------------------------------------------------------
subroutine compute_deriv_nth_order(x,f0,f2)
use glbl_prmtrs
use mod_wind
double precision, intent(in) :: x
double precision, intent(out) :: f0, f2
double precision :: v
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
call get_v(x,v)
f0=1.d0/(x**2.d0*v)
f2= 6.d0              /(x**4.d0*v)&
   -6.d0*beta_        /(x**5.d0*v**((beta_+1.d0)/beta_))&
   +beta_*(beta_+1.d0)/(x**6.d0*v**((beta_+2.d0)/beta_))

end subroutine compute_deriv_nth_order
! -----------------------------------------------------------------------------------

end module rdm
