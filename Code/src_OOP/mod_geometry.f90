module mod_geometry
contains

! -----------------------------------------------------------------------------------
! To determine whether a given clump is on a ray which intercepts the surface
! built by the LOS integrated over an orbit. If not, reject and pick up another.
! Stricly speaking, one should account for the size of the clump. Here,
! I just take a fiducial "twice the size of a clump @ 2 stellar radii".
! I also set that beyond a fiducial distance to the stellar center,
! taken to be the distance beyond which clumps are deleted, the interception
! does not matter anymore.
! -----------------------------------------------------------------------------------
subroutine is_relevant(th,rlvnt)
use glbl_prmtrs
double precision, intent(in)  :: th
logical, intent(out) :: rlvnt
! double precision :: b ! impact parameter of the ray w.r.t. LOS surface
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! th
! inclnsn
! omega_i

! if (b>2.d0*clump_rad_) rlvnt = .false.

end subroutine is_relevant
! -----------------------------------------------------------------------------------


end module mod_geometry
