module mod_NH
contains

! -----------------------------------------------------------------------------------
! Compute and save NH.
! Beware, pos_cl is in spherical!
! -----------------------------------------------------------------------------------
subroutine compute_NH(Ncl,pos_cl,R_cl,dens_cl,is_init)
use glbl_prmtrs
use mod_dynamics_X
use IO
use mod_cart_sph
use miscellaneous
integer, intent(in) :: Ncl
double precision, intent(in)  :: pos_cl(Ncl,3), R_cl(Ncl), dens_cl(Ncl)
logical, intent(in), optional :: is_init
double precision :: pos_cl_cart(Ncl,3)
double precision :: NH(Nphases_)
double precision :: pos_X(3), phase, dphase, h
double precision :: b2 ! impact parameter squared of clumps
double precision :: bToStr ! impact parameter (squared) of the star
integer :: i, k, j_bin, N_bin
double precision :: q ! common ratio of geometric serie
double precision, allocatable :: NH_bin(:)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! radial shells logarithmically subdivided in 10 segments from Rstar_ to 10Rstar_
! and then, we scale up to dist_max_cl_ such as the common ratio remains the same
! (ie same subdivision for the initial range)
N_bin = int( (10.d0*Rstar_/Rstar_)*10.d0 * ( dlog(dist_max_cl_/Rstar_) / dlog(10.d0*Rstar_/Rstar_) ) )
q = (dist_max_cl_/Rstar_)**(1.d0/dble(N_bin))
allocate(NH_bin(N_bin))
NH_bin=0.d0

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
        ! @ initialization, compute the relative contribution of each radial bin
        ! to the final NH for phase = 0 (ie k==1).
        if (present(is_init) .and. k==1) then ! by default, if not present, not initialization
          ! which radial shell does the clump belong to?
          ! j_bin = 1 + int( dlog(pos_cl(i,1)/Rstar_) / ((dlog(dist_max_cl_/Rstar_))/dble(dist_max_cl_*10.d0+0.00000001d0)) )
          j_bin = 1 + int( dlog(pos_cl(i,1)/Rstar_) / dlog(q) ) ! ((dlog(dist_max_cl_/Rstar_))/dble(dist_max_cl_*10.d0+0.00000001d0)) )
          ! print*, k, j_bin, pos_cl(i,1), dist_max_cl_, dble(dist_max_cl_*10.d0), maxval(pos_cl(:,1))
          ! If the clump is @ a radius r < Rstar_, we have a serious problem
          if (j_bin<1) then
            print*, j_bin, '<', 1
            call crash("Clump under the lower bound of histogram (which is r=Rstar_)")
          ! If the clump is @ a radius r > dist_max_cl_, it is because it got into this
          ! region in the immediately previous time step and will be destroyed immediately
          ! after "call compute_NH" => do not account for it in the histogram.
        elseif (j_bin>N_bin) then ! int(dist_max_cl_*10.d0)) then
            continue
          endif
          NH_bin(j_bin)=NH_bin(j_bin)+h*dens_cl(i)
        endif
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

if (present(is_init)) call save_NH_shells(N_bin,NH_bin)

deallocate(NH_bin)

end subroutine compute_NH
! -----------------------------------------------------------------------------------


end module mod_NH
