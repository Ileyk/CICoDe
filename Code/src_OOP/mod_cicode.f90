module mod_cicode
contains

! -----------------------------------------------------------------------------------
subroutine cicode
use glbl_prmtrs
use mod_init
use IO
use miscellaneous
integer :: Ncl
double precision, allocatable :: pos_cl(:,:), R_cl(:), dens_cl(:)
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call cpu_time(chrono_2)

t0_=0.d0
t_=t0_
Ncl=Ncl0_

allocate(pos_cl(Ncl,3),R_cl(Ncl),dens_cl(Ncl))

call set_ini_clumps(Ncl,pos_cl,R_cl,dens_cl)

! Save the initial spherical coordinates of the clumps, along with their radii
call save_histograms(Ncl,pos_cl,R_cl)

! Save the initial Cartesian coordinates of the clumps, along with their radii
call save_pos(Ncl,pos_cl,R_cl)

call chrono(chrono_2,chrono_2_mess)

call cpu_time(chrono_3)

do while (t_<time_max_)

  call advance(Ncl,pos_cl,R_cl,dens_cl)

  ! CHEAT to debug
  if (t_==t0_+dt_) then
    call save_pos(Ncl,pos_cl,R_cl)
    ! call compute_NH
    t0_=t0_+dt_
  endif

enddo

call chrono(chrono_3,chrono_3_mess)

close(1)

end subroutine cicode
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine advance(Ncl,pos_cl,R_cl,dens_cl)
use glbl_prmtrs
use mod_clumps
use mod_integration
integer, intent(inout) :: Ncl
double precision, allocatable, intent(inout) :: pos_cl(:,:), R_cl(:), dens_cl(:)
double precision, allocatable :: pos_cl_old(:,:)
integer :: i
double precision :: dt_dyn, r
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call get_dt(Ncl,pos_cl,R_cl,dt_dyn)

call add_delete_clumps(dt_dyn,Ncl,pos_cl,R_cl,dens_cl)

! Save previous positions, to be used in expand_clumps
allocate(pos_cl_old(Ncl,3)) ! now that clumps have been deleted / added
pos_cl_old=pos_cl

do i=1,Ncl
  r=pos_cl(i,1)
  call rk4(r,dt_dyn)
  pos_cl(i,1)=r
enddo

call expand_clumps(Ncl,pos_cl_old,pos_cl,R_cl,dens_cl)

deallocate(pos_cl_old)

if (do_merge_) call merge_clumps(Ncl,pos_cl,R_cl,dens_cl)

t_ = t_ + dt_dyn

end subroutine advance
! -----------------------------------------------------------------------------------

end module mod_cicode
