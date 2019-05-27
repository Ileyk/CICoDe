module mod_NH
contains

! -----------------------------------------------------------------------------------
! Compute and save NH.
! Beware, pos_cl is in spherical!
! -----------------------------------------------------------------------------------
subroutine compute_NH(Ncl,pos_cl,R_cl,dens_cl)
use glbl_prmtrs
use mod_dynamics_X
use IO
use mod_cart_sph
integer, intent(in) :: Ncl
double precision, intent(in)  :: pos_cl(Ncl,3), R_cl(Ncl), dens_cl(Ncl)
double precision :: pos_cl_cart(Ncl,3)
double precision :: NH(Nphases_)
double precision :: pos_X(3), phase, dphase, h
double precision :: b2 ! impact parameter squared of clumps
double precision :: bToStr ! impact parameter (squared) of the star
integer :: i, k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Get the Cartesian positions of the clumps
call sph_2_cart(Ncl,pos_cl,pos_cl_cart)

NH=0.d0

dphase=(2.d0*dpi/dble(Nphases_))
phase=0.d0 ! beware here, phases go from 0 to 2pi-dpi (no half-step)

do k=1,Nphases_

  ! Where is the X-ray source @ this phase?
  call get_pos_X(phase,pos_X)

  ! Is it eclipsed?
  if (0.d0>pos_X(3)) then
    bToStr = pos_X(1)**2.d0 + pos_X(2)**2.d0
    if (bToStr<Rstar_**2.d0) then
      NH(k)=bigdble
      continue
    endif
  endif

  do i=1,Ncl
    ! 1st condition : the clump is between the X-ray source and the observer
    if (pos_cl_cart(i,3)>pos_X(3)) then
      b2 = (pos_cl_cart(i,1)-pos_X(1))**2.d0 + (pos_cl_cart(i,2)-pos_X(2))**2.d0
      ! 2nd condition : the impact parameter is smaller than the clump radius
      if (b2<R_cl(i)**2.d0) then
        h=dsqrt(R_cl(i)**2.d0-b2)
        NH(k)=NH(k)+h*dens_cl(i)
      endif
    endif
  enddo

  ! Save the position w/ phase origin phi=0 @ t=0
  ! (ie the diagonale of the NH(Nphases_,Nphases_) matrix)
  ! == 1 when called just after initialization
  ! +=1 until == Nphases_ that we do not save since redundant
  ! w/ initial state (both have phase = 0).
  if (k==int(t_/dt_+smalldble)+1) then
    call save_pos_X(pos_X)
  endif

  phase=phase+dphase

enddo

call save_NH(NH)

end subroutine compute_NH
! -----------------------------------------------------------------------------------


end module mod_NH
