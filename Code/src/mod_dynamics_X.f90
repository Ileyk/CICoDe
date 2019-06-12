module mod_dynamics_X
contains

! -----------------------------------------------------------------------------------
! Gives the Cartesian coordinates (x,y,z) of the compact object w.r.t. the star.
! Phase 0 is taken as the moment when the object is
! the closest from us (inferior conjunction).
! z is along the LOS, with z increasing towards us.
! CHEAT
! Works only for a circular edge-on orbit for now
! Phase is the angular position (in radians)
! -----------------------------------------------------------------------------------
subroutine get_pos_X(phase,pos_X)
use glbl_prmtrs
double precision, intent(in)  :: phase
double precision, intent(out) :: pos_X(3)
double precision :: x, y, z
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

x=a_*dsin(phase)
y=0.d0
z=a_*dcos(phase)

pos_X(1)= x
pos_X(2)= y*dsin(inclnsn_)+z*dcos(inclnsn_)
pos_X(3)=-y*dcos(inclnsn_)+z*dsin(inclnsn_)

end subroutine get_pos_X
! -----------------------------------------------------------------------------------


end module mod_dynamics_X
