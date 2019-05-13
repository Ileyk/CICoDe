module mod_integration
contains

! -----------------------------------------------------------------------------------
subroutine rk4(r,dt_dyn)
use glbl_prmtrs
use mod_wind
double precision, intent(in) :: dt_dyn
double precision, intent(inout) :: r
double precision :: v
double precision:: rt, vt, vm
double precision:: dt2, dt6

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

dt2=half*dt_dyn
dt6=dt_dyn/6.d0

call get_v(r,v)
rt=r+dt2*v
call get_v(rt,vt)
rt=r+dt2*vt
call get_v(rt,vm)
rt=r+dt_dyn*vm
vm=vm+vt
call get_v(rt,vt)
r=r+dt6*(v+vt+2.d0*vm)

end subroutine rk4
! -----------------------------------------------------------------------------------


end module mod_integration
