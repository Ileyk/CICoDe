module mod_cicode
contains

! -----------------------------------------------------------------------------------
subroutine cicode
use glbl_prmtrs
use mod_init
use IO
integer :: Ncl
double precision, allocatable :: pos_cl(:,:), R_cl(:), dens_cl(:)
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

t0_=0.d0
t_=t0_
Ncl=Ncl0_

allocate(pos_cl(Ncl,3),R_cl(Ncl),dens_cl(Ncl))

call set_ini_clumps(Ncl,pos_cl,R_cl,dens_cl)

! Save the initial spherical coordinates of the clumps, along with their radii
call save_histograms(Ncl,pos_cl,R_cl)

call save_pos(Ncl,pos_cl,R_cl)

do while (t_<time_max_)

  call advance(Ncl,pos_cl,R_cl,dens_cl)

  ! CHEAT to debug
  if (t_==t0_+dt_) then
    call save_pos(Ncl,pos_cl,R_cl)
    ! call compute_NH
    t0_=t0_+dt_
  endif

enddo

close(1)

end subroutine cicode
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine advance(Ncl,pos_cl,R_cl,dens_cl)
use glbl_prmtrs
use mod_clumps
use mod_integration
integer, intent(inout) :: Ncl
double precision, intent(inout) :: pos_cl(Ncl,3), R_cl(Ncl), dens_cl(Ncl)
double precision :: pos_cl_old(Ncl,3)
integer :: i
double precision :: dt_dyn, r
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call get_dt(Ncl,pos_cl,R_cl,dt_dyn)

! call add_delete_clumps(dt_dyn,Ncl,pos_cl,R_cl,dens_cl)

! Save previous positions, to be used in expand_clumps
pos_cl_old=pos_cl

do i=1,Ncl
  r=pos_cl(i,1)
  call rk4(r,dt_dyn)
  pos_cl(i,1)=r
enddo

call expand_clumps(Ncl,pos_cl_old,pos_cl,R_cl,dens_cl)

! call merge_clumps(Ncl,pos_cl,R_cl)

t_ = t_ + dt_dyn

end subroutine advance
! -----------------------------------------------------------------------------------

end module mod_cicode
