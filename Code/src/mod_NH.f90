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
integer :: i, k, j_bin, N_bin
double precision :: rrr(Ncl) !, NH_bin(int(dist_max_cl_*10.d0))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

N_bin = int((dist_max_cl_/Rstar_)*10.d0)

rrr(:)=pos_cl(:,1)

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
        ! print*, k, rrr(i), dsqrt(b2), R_cl(i), h, minval(rrr), maxval(rrr)

        ! radial shells logarithmically subdivided in 10*dist_max_cl_/R_star segments
        ! (eg 100 for dist_max_cl_=10*R_star) from 1 to dist_max_cl_
        j_bin = 1 + int( dlog(rrr(i)/1.d0) / ((dlog(dist_max_cl_/1.d0))/dble(dist_max_cl_*10.d0+0.00000001d0)) )
        ! print*, j_bin, pos_cl(i,1), rrr(i), dist_max_cl_, dble(dist_max_cl_*10.d0)
        ! if (j_bin<1 .or. j_bin>int(dist_max_cl_*10.d0)) then
        !   print*, 'wtf', j_bin
        !   stop
        ! endif
        ! NH_bin(j_bin)=NH_bin(j_bin)+h*dens_cl(i)
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
    ! Compute the histogram of NH along the LOS (subdivided in log intervals)
    ! to evaluate where we should stop the integration of NH
    ! call save_NH_hist_LOS()
  endif

  phase=phase+dphase

enddo

call save_NH(NH)

end subroutine compute_NH
! -----------------------------------------------------------------------------------


end module mod_NH
