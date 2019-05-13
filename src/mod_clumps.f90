module mod_clumps
contains

! -----------------------------------------------------------------------------------
! Compute dynamical dt, dt_dyn, as the minimum self-crossing time
! pos_cl : coordinates in spherical units
! -----------------------------------------------------------------------------------
subroutine get_dt(Ncl,pos_cl,R_cl,dt_dyn)
use glbl_prmtrs
use mod_wind
integer, intent(in) :: Ncl
double precision, intent(in)  :: pos_cl(Ncl,3), R_cl(Ncl)
double precision, intent(out) :: dt_dyn
double precision :: r1_cl, v_cl, dt_tmp
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

dt_dyn=bigdble

do i=1,Ncl
  r1_cl=pos_cl(i,1)
  call get_v(r1_cl,v_cl)
  dt_tmp=R_cl(i)/v_cl
  dt_dyn=min(dt_dyn,dt_tmp)
enddo

! If needed, shorten dt_dyn so as we match one of the phase-time positions,
! separated by dt_
! t0_ is the previous phase-time position
if ( (t_+dt_dyn-(t0_+dt_)) * (t_-(t0_+dt_)) < 0.d0 ) then
  dt_dyn=t0_+dt_-t_
endif

end subroutine get_dt
! -----------------------------------------------------------------------------------

! ! -----------------------------------------------------------------------------------
! ! Here, for the moment, all the clumps are initialized :
! !   - with the same size and density
! !   - @ the same position rini_ such as v @ rini_ = 1% of terminal wind speed
! ! For now, the criteria to delete clumps is (too) simple : r1 > 10 orb. sep.
! ! Notice that when clumps are added, we first give a chance to the flow to advance
! ! before computing who merged with who.
! ! -----------------------------------------------------------------------------------
! subroutine add_delete_clumps(dt_dyn,Ncl,pos_cl,R_cl,dens_cl)
! use glbl_prmtrs
! use rdm
! use miscellaneous
! use mod_wind
! double precision, intent(in) :: dt_dyn
! integer, intent(inout) :: Ncl
! double precision, intent(inout) :: pos_cl(Ncl,3), R_cl(Ncl), dens_cl(Ncl)
! ! # of clumps added and deleted respectively
! integer :: dNcl_p, dNcl_m
! integer :: i
! double precision :: pos_cl_tmp(Ncl,3), R_cl_tmp(Ncl), dens_cl_tmp(Ncl)
! double precision :: th, ph
! logical :: deleted(Ncl)
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! deleted=.false.
!
! ! Who should we delete?
! do i=1,Ncl
!   if (pos_cl(i,1)>10.d0*a_) then
!     deleted(i)=.true.
!   endif
! enddo
!
! ! # of clumps to delete
! dNcl_m=count(deleted .eqv. .true.)
!
! ! # of clumps to add
! dNcl_p = int ( dt_dyn * (mass_fraction_*Mdot_) / clump_mass_ )
! if (dNcl_p==0) call crash("Beware, # of clumps created over dt_dyn < 1...")
!
! ! Store previous arrays...
! ! pos_cl_tmp=pos_cl
! ! R_cl_tmp=R_cl
! ! dens_cl_tmp=dens_cl
! ! ! ... before deallocating them...
! ! deallocate(pos_cl,R_cl,dens_cl)
! ! ! ... and reallocating them w/ the new size
! ! allocate(pos_cl(Ncl,3),R_cl(Ncl),dens_cl(Ncl))
!
! call delete_clumps(deleted,Ncl,pos_cl,R_cl,dens_cl)
!
! ! New total # of clumps
! Ncl=Ncl+dNcl_p
!
! call add_clumps(dNcl_p,Ncl,pos_cl,R_cl,dens_cl)
!
! end subroutine add_delete_clumps
! ! -----------------------------------------------------------------------------------
!
! ! -----------------------------------------------------------------------------------
! ! Delete clumps ie paste only those which are not to be deleted,
! ! starting @ the beginning of arrays of size Ncl
! ! (which does not include the added clumps yet)
! ! - deleted : contains the information on which clump should be deleted
! ! - Ncl_tmp : size of the tmp arrays in which have been stored the values
! !   of the initial arrays
! ! -----------------------------------------------------------------------------------
! subroutine delete_clumps(deleted,Ncl,pos_cl,R_cl,dens_cl)
! use glbl_prmtrs
! use mod_wind
! use miscellaneous
! integer, intent(inout) :: Ncl
! logical, intent(in) :: deleted(Ncl)
! double precision, intent(inout) :: pos_cl(Ncl,3), R_cl(Ncl), dens_cl(Ncl)
! double precision :: pos_cl_tmp(Ncl,3), R_cl_tmp(Ncl), dens_cl_tmp(Ncl)
! integer :: i, j, dNcl
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! dNcl=count(deleted .eqv. .true.)
!
! ! Store previous arrays...
! pos_cl_tmp=pos_cl
! R_cl_tmp=R_cl
! dens_cl_tmp=dens_cl
! ! ... before deallocating them...
! deallocate(pos_cl,R_cl,dens_cl)
! Ncl=Ncl-dNcl
! ! ... and reallocating them w/ the new size
! allocate(pos_cl(Ncl,3),R_cl(Ncl),dens_cl(Ncl))
!
! j=1
! do i=1,Ncl+dNcl
!   if (.not. deleted(i)) then
!     pos_cl(j,:)=pos_cl_tmp(j,:)
!     R_cl(j)=R_cl_tmp(j)
!     dens_cl(j)=dens_cl_tmp(j)
!     j=j+1
!   endif
! enddo
!
! if (j/=Ncl+1) call crash("Problem w/ the # of clumps being deleted")
!
! end subroutine delete_clumps
! ! -----------------------------------------------------------------------------------
!
! ! -----------------------------------------------------------------------------------
! ! Add dNcl new clumps @ the end of the arrays of size Ncl
! ! (which already includes the added clumps)
! ! -----------------------------------------------------------------------------------
! subroutine add_clumps(dNcl,Ncl,pos_cl,R_cl,dens_cl)
! use glbl_prmtrs
! use rdm
! use mod_wind
! integer, intent(in) :: dNcl
! integer, intent(inout) :: Ncl
! double precision, intent(inout) :: pos_cl(Ncl,3), R_cl(Ncl), dens_cl(Ncl)
! integer :: i
! double precision :: th, ph
! logical :: rlvnt=.true.
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! do i=Ncl-dNcl+1,Ncl
!   pos_cl(i,1)=rini_
!   ! Pick up @ random angular position (over dphi and d(cos(th)) )
!   call get_flat(-1.d0,1.d0,th)
!   if (rlvnt_clumps_) then
!     call is_relevant(th,rlvnt)
!     if (.not. rlvnt) PICK UP NEW TH AND PH
!   endif
!   call get_flat(0.d0,2.d0*dpi,ph)
!   ! CHEAT to debug
!   pos_cl(i,2)=dpi/2.d0 ! dacos(th)
!   pos_cl(i,3)=ph
!   R_cl(i)=clump_rad_
!   dens_cl(i)=clump_dens_
! enddo
!
! end subroutine add_clumps
! ! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Account for the new position and expansion such as
!   - each clump mass density dilutes as 1/r2v
!   - the mass of each clump is conserved so as dens_cl x R_cl**3 = cst
! to compute the new radii and densities
! -----------------------------------------------------------------------------------
subroutine expand_clumps(Ncl,pos_cl_old,pos_cl,R_cl,dens_cl)
use glbl_prmtrs
use mod_wind
use miscellaneous
integer, intent(in) :: Ncl
double precision, intent(in) :: pos_cl(Ncl,3), pos_cl_old(Ncl,3)
double precision, intent(inout) :: R_cl(Ncl), dens_cl(Ncl)
double precision :: r, r_old, v, v_old
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

do i=1,Ncl
  r=pos_cl(i,1)
  r_old=pos_cl_old(i,1)
  call get_v(r,v)
  call get_v(r_old,v_old)
  R_cl(i)=R_cl(i)    * (r/r_old)**(2.d0/3.d0)*(v/v_old)**(1.d0/3.d0)
  dens_cl(i)=dens_cl(i)/((r/r_old)**2.d0*(v/v_old))
enddo

end subroutine expand_clumps
! -----------------------------------------------------------------------------------

! ! -----------------------------------------------------------------------------------
! ! In the way we compute the mergers, we assume "serial" mergers possible ie
! ! "A merges w/ B and B merges w/ C but A does not merge w/ C?
! ! Then we say that A merges w/ C too" to avoid order dependence
! ! (although there is still some because of d<R_cl(i)+R_cl(j) and change of
! ! position conditions)
! ! Although it is moved in position, a merged clump is always assumed to be "relevant"
! ! -----------------------------------------------------------------------------------
! subroutine merge_clumps(Ncl,pos_cl,R_cl,dens_cl)
! use glbl_prmtrs
! use mod_wind
! use miscellaneous
! integer, intent(inout) :: Ncl
! double precision, intent(inout) :: pos_cl(Ncl,3), R_cl(Ncl), dens_cl(Ncl)
! integer :: pairs_of_mergers(Ncl)
! integer :: dNcl=0, i, j, ii, jj
! double precision :: d
! logical :: deleted(Ncl)
! ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! deleted=.false.
! pairs_of_mergers=0
! do i=1,Ncl-1
!   do j=i+1,Ncl ! Check all the clumps "ahead" only since merging is commutative
!     ! Compute distance w/ any other clump
!     ! CHEAT
!     d=dsqrt((pos_cl(i,1)*dsin(pos_cl(i,2))*dcos(pos_cl(i,3))&
!             -pos_cl(j,1)*dsin(pos_cl(j,2))*dcos(pos_cl(j,3)))**2.d0+&
!             (pos_cl(i,1)*dsin(pos_cl(i,2))*dsin(pos_cl(i,3))&
!             -pos_cl(j,1)*dsin(pos_cl(j,2))*dsin(pos_cl(j,3)))**2.d0+&
!             (pos_cl(i,1)*dcos(pos_cl(i,2))&
!             -pos_cl(j,1)*dcos(pos_cl(j,2)))**2.d0)
!     ! If smaller than sum of radii, merge
!     if (d<R_cl(i)+R_cl(j)) then
!       ! By convention, the merged clump is stored into # i
!       ! while clump # j is deleted
!       if (deleted(j)) then ! beware, this clump has already been merged! => > 2 clumps merger
!         ! w/ a clump k (=pairs_of_mergers(j)) < i => put all in k
!         jj=pairs_of_mergers(j)
!         if (deleted(jj)) call crash("Pb with deleted(j) in merge clumps")
!         deleted(i)=.true.
!         pairs_of_mergers(i)=jj
!         R_cl(jj)=R_cl(jj)*2.d0**(1.d0/3.d0)
!         pos_cl(jj,:)=0.5d0*(pos_cl(i,:)+pos_cl(jj,:))
!       elseif (deleted(i)) then ! beware, this clump has already been merged! => > 2 clumps merger
!         ii=pairs_of_mergers(i)
!         if (deleted(ii)) call crash("Pb with deleted(i) in merge clumps")
!         deleted(j)=.true.
!         pairs_of_mergers(j)=ii
!         R_cl(ii)=R_cl(ii)*2.d0**(1.d0/3.d0)
!         pos_cl(ii,:)=0.5d0*(pos_cl(ii,:)+pos_cl(j,:))
!       else
!         pairs_of_mergers(i)=j ! j merges w/ i
!         pairs_of_mergers(j)=i ! i merges w/ j
!         deleted(j)=.true.
!       endif
!       dNcl=dNcl+1
!     endif
!   enddo
! enddo
!
! call delete_clumps(deleted,Ncl,pos_cl,R_cl,dens_cl)
!
! end subroutine merge_clumps
! ! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
! Beware, the initialization assumes no merger and so should be different
! from the stochastic steady-state which is reached later on.
! -----------------------------------------------------------------------------------
subroutine set_ini_radii(Ncl,pos_cl,R_cl,dens_cl)
use glbl_prmtrs
use rdm
use mod_wind
integer, intent(in) :: Ncl
double precision, intent(out) :: pos_cl(Ncl,3), R_cl(Ncl), dens_cl(Ncl)
double precision :: pos_cl_old(Ncl,3), R_cl_old(Ncl), dens_cl_old(Ncl)
double precision :: r1
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

do i=1,Ncl
  call get_one_over_r2v(r1)
  pos_cl(i,1)=r1
  print*, i, r1
  write(7,*), r1
enddo

stop

pos_cl_old(:,2)=pos_cl(:,2)
pos_cl_old(:,3)=pos_cl(:,3)
! since the clump properties set in the par file are given @ 2 stellar radii
pos_cl_old(:,1)=2.d0

R_cl=clump_rad_
dens_cl=clump_dens_

call expand_clumps(Ncl,pos_cl_old,pos_cl,R_cl,dens_cl)

end subroutine set_ini_radii
! -----------------------------------------------------------------------------------

end module mod_clumps
