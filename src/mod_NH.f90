module mod_NH
contains

! -----------------------------------------------------------------------------------
! Compute NH
! -----------------------------------------------------------------------------------
subroutine compute_NH(Ncl,pos_cl,R_cl,dens_cl,NH)
use glbl_prmtrs
use mod_dynamics_X
integer, intent(in) :: Ncl
double precision, intent(in)  :: pos_cl(Ncl,3), R_cl(Ncl), dens_cl(Ncl)
double precision, intent(out) :: NH(Nphases_,Nphases_)
double precision :: pos_X(3), phase
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

dphase=(2.d0*dpi/dble(Nphases_))
do i=1,Nphases_
  phase=dphase/2.d0+dble(i-1)*dphase
  call get_pos_X(phase,pos_X)
  do while (?)
    ! Work uniform segment by uniform segment (ie wind or clump)
    ! 1st, where is the X-ray source? In the wind or in a clump?

    ! Then, what is the distance along the LOS before it changes?

    ! Account for this section, and subdivide it if too long

    NH=NH+local_dens*dl
  enddo

enddo

end subroutine compute_NH
! -----------------------------------------------------------------------------------


end module mod_NH
