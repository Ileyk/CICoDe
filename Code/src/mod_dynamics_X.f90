module mod_dynamics_X
contains

! -----------------------------------------------------------------------------------
! Gives the positions (x,y,z) of the compact object with respect to the star.
! Phase 0 is taken as the moment when the object is
! the closest from us (inferior conjunction).
! z is along the LOS, with z increasing towards us.
! CHEAT
! Works only for a circular edge-on orbit for now
! -----------------------------------------------------------------------------------
subroutine get_pos_X(phase,pos_X)
use glbl_prmtrs
double precision, intent(in)  :: phase
double precision, intent(out) :: pos_X(3)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

pos_X(1)=a_*dsin(phase)
pos_X(2)=0.d0
pos_X(3)=a_*dcos(phase)

end subroutine get_pos_X
! -----------------------------------------------------------------------------------


end module mod_dynamics_X
