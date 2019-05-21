module mod_wind
contains

! -----------------------------------------------------------------------------------
subroutine get_v(r,v)
use glbl_prmtrs
double precision, intent(in)  :: r
double precision, intent(out) :: v
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

v=(1.d0-1.d0/r)**beta_

end subroutine get_v
! -----------------------------------------------------------------------------------

end module mod_wind
