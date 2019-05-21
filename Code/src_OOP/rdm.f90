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
