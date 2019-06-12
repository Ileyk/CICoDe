module mod_cicode
contains

! -----------------------------------------------------------------------------------
subroutine cicode
use glbl_prmtrs
use mod_init
use IO
use miscellaneous
use mod_NH
use mod_clumps
integer :: Ncl
double precision, allocatable :: pos_cl(:,:), R_cl(:), dens_cl(:)
integer :: i
character(LEN=200) :: string
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call cpu_time(chrono_2)

if (restart_indx_<0) then ! if no restart

  t0_=0.d0
  t_=t0_
  Ncl=Ncl0_

  allocate(pos_cl(Ncl,3),R_cl(Ncl),dens_cl(Ncl))

  call set_ini_clumps(Ncl,pos_cl,R_cl,dens_cl)

  ! delete immediately the unappropriate (ie never on LOS) clumps
  ! (but since is_init=.true., do not add any clump)
  call add_delete_clumps(0.d0,Ncl,pos_cl,R_cl,dens_cl,.true.)

  ! Save the initial spherical coordinates of the clumps, along with their radii
  call save_histograms(Ncl,pos_cl,R_cl)

  ! Save the initial Cartesian coordinates of the clumps, along with their radii
  call save_pos(Ncl,pos_cl,R_cl,dens_cl)

  ! Compute the first NH, in initial position.
  ! Since is_init==.true., we will also compute the relative contribution
  ! of each radial layer to NH, to justify where we set dist_max_cl_
  call compute_NH(Ncl,pos_cl,R_cl,dens_cl,.true.)

else

  ! Read time & # of clumps from positions file "restart_indx_"
  call load_clumps(Ncl,pos_cl,R_cl,dens_cl)

  t_=t0_

endif

if (init_only_) call crash('Success, initialization finished.')

call chrono(chrono_2,chrono_2_mess)

print*, 'Initialization over - - - - - - - - - -'

call cpu_time(chrono_3)

do while (t_<time_max_)

  call advance(Ncl,pos_cl,R_cl,dens_cl)

  if (t_==t0_+dt_) then
    if (.not. plot_cl_only_) call compute_NH(Ncl,pos_cl,R_cl,dens_cl)
    ! Is it time to save positions and instantaneous NH?
    if (int((t_+smalldble)/dt_)/(Nphases_/Nsave_)/=int((t_-smalldble)/dt_)/(Nphases_/Nsave_)) then
      write(string,'((I6.6),(a),(I6.6))') int((t_+smalldble)/dt_)/(Nphases_/Nsave_), '/', Nsave_
      call followup(string)
      call save_pos(Ncl,pos_cl,R_cl,dens_cl)
    endif
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

! print*, Ncl

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
