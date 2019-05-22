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

! -----------------------------------------------------------------------------------
! Here, for the moment, all the clumps are initialized :
!   - with the same size and density
!   - @ the same position rini_ such as v @ rini_ = 1% of terminal wind speed
! For now, the criteria to delete clumps is (too) simple : r1 > dist_max_cl_
! Notice that when clumps are added, we first give a chance to the flow to advance
! before computing who merged with who.
! -----------------------------------------------------------------------------------
subroutine add_delete_clumps(dt_dyn,Ncl,pos_cl,R_cl,dens_cl)
use glbl_prmtrs
use rdm
use miscellaneous
use mod_wind
double precision, intent(in) :: dt_dyn
integer, intent(inout) :: Ncl
double precision, allocatable, intent(inout) :: pos_cl(:,:), R_cl(:), dens_cl(:)
! # of clumps added and deleted respectively
integer :: dNcl_p, dNcl_m
integer :: i
double precision :: pos_cl_tmp(Ncl,3), R_cl_tmp(Ncl), dens_cl_tmp(Ncl)
double precision :: th, ph
logical :: deleted(Ncl)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

deleted=.false.

! Who should we delete?
do i=1,Ncl
  if (pos_cl(i,1)>dist_max_cl_) then
    deleted(i)=.true.
  endif
enddo

! # of clumps to delete
! dNcl_m=count(deleted .eqv. .true.)

! ! # of clumps to add
cmltd_clump_ = cmltd_clump_ + dt_dyn * Ndot_
dNcl_p = int ( cmltd_clump_ )
cmltd_clump_ = cmltd_clump_ - dble(int(cmltd_clump_))
if (cmltd_clump_>1.d0 .or. cmltd_clump_<0.d0) call crash('cmltd_clump_ can not be > 1 or < 0')

call delete_clumps(deleted,Ncl,pos_cl,R_cl,dens_cl)

! New total # of clumps
! Ncl=Ncl+dNcl_p

call add_clumps(dNcl_p,Ncl,pos_cl,R_cl,dens_cl)

end subroutine add_delete_clumps
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Delete clumps ie paste only those which are not to be deleted,
! starting @ the beginning of arrays of size Ncl
! (which does not include the added clumps yet)
! - deleted : contains the information on which clump should be deleted
! - Ncl_tmp : size of the tmp arrays in which have been stored the values
!   of the initial arrays
! -----------------------------------------------------------------------------------
subroutine delete_clumps(deleted,Ncl,pos_cl,R_cl,dens_cl)
use glbl_prmtrs
use mod_wind
use miscellaneous
integer, intent(inout) :: Ncl
logical, intent(in) :: deleted(Ncl)
double precision, allocatable, intent(inout) :: pos_cl(:,:), R_cl(:), dens_cl(:)
double precision :: pos_cl_tmp(Ncl,3), R_cl_tmp(Ncl), dens_cl_tmp(Ncl)
integer :: i, j, dNcl
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

dNcl=count(deleted .eqv. .true.)
if (dNcl==0) return ! no clump to be deleted

! Store previous arrays...
pos_cl_tmp=pos_cl
R_cl_tmp=R_cl
dens_cl_tmp=dens_cl
! ... before deallocating them...
deallocate(pos_cl,R_cl,dens_cl)
Ncl=Ncl-dNcl
! ... and reallocating them w/ the new size
allocate(pos_cl(Ncl,3),R_cl(Ncl),dens_cl(Ncl))

j=1
do i=1,Ncl+dNcl
  if (.not. deleted(i)) then
    pos_cl(j,:)=pos_cl_tmp(i,:)
    R_cl(j)=R_cl_tmp(i)
    dens_cl(j)=dens_cl_tmp(i)
    j=j+1
  endif
enddo

if (j/=Ncl+1) call crash("Problem w/ the # of clumps being deleted")

end subroutine delete_clumps
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Add dNcl new clumps @ the end of the arrays of size Ncl
! (which already includes the added clumps)
! -----------------------------------------------------------------------------------
subroutine add_clumps(dNcl,Ncl,pos_cl,R_cl,dens_cl)
use glbl_prmtrs
use rdm
use mod_wind
integer, intent(in) :: dNcl
integer, intent(inout) :: Ncl
double precision, allocatable, intent(inout) :: pos_cl(:,:), R_cl(:), dens_cl(:)
double precision :: pos_cl_tmp(Ncl,3), R_cl_tmp(Ncl), dens_cl_tmp(Ncl)
integer :: i
double precision :: th, ph
logical :: rlvnt=.true.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Nothing to do
if (dNcl==0) return

! Store previous arrays...
pos_cl_tmp=pos_cl
R_cl_tmp=R_cl
dens_cl_tmp=dens_cl
! ... before deallocating them...
deallocate(pos_cl,R_cl,dens_cl)
Ncl=Ncl+dNcl
! ... and reallocating them w/ the new size
allocate(pos_cl(Ncl,3),R_cl(Ncl),dens_cl(Ncl))

! Put again old values...
do i=1,Ncl-dNcl
  pos_cl(i,:)=pos_cl_tmp(i,:)
  R_cl(i)=R_cl_tmp(i)
  dens_cl(i)=dens_cl_tmp(i)
enddo

! ... and create new clumps
do i=Ncl-dNcl+1,Ncl
  pos_cl(i,1)=rini_
  ! Pick up @ random angular position (over dphi and d(cos(th)) )
  call get_flat(-1.d0,1.d0,th)
  ! if (rlvnt_clumps_) then
  !   call is_relevant(th,rlvnt)
  !   if (.not. rlvnt) PICK UP NEW TH AND PH
  ! endif
  call get_flat(0.d0,2.d0*dpi,ph)
  ! CHEAT to debug
  pos_cl(i,2)=dpi/2.d0 ! dacos(th)
  pos_cl(i,3)=ph
  R_cl(i)=clump_rad_ini_
  dens_cl(i)=clump_dens_ini_
enddo

end subroutine add_clumps
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Account for the new position and expansion such as
!   - each clump mass density dilutes as 1/r2v
!   - the mass of each clump is conserved so as dens_cl x R_cl**3 = cst
! to compute the new radii and densities.
! See notes : Lorenzo Ducci + 09 & Jon Sundqvist + 12 disagree on the evolution of the
! clump radius with the distance to the star (=> consequences for mass density too).

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
  select case(rad_evol_)
  case('lorenzo')
    R_cl(i)=R_cl(i)    * (r/r_old)**(2.d0/3.d0)*(v/v_old)**(1.d0/3.d0)
    dens_cl(i)=dens_cl(i)/((r/r_old)**2.d0*(v/v_old))
  case('jon') ! density evolution such as constant clump mass enforced
    R_cl(i)=R_cl(i)    * (r/r_old)
    dens_cl(i)=dens_cl(i)/(r/r_old)**3.d0
  case default
    call crash('Which evolution of clump radii with distance to star?')
  end select
enddo

end subroutine expand_clumps
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! In the way we compute the mergers, we assume "serial" mergers possible ie
! "A merges w/ B and B merges w/ C but A does not merge w/ C?
! Then we say that A merges w/ C too" to avoid order dependence
! (although there is still some because of d<R_cl(i)+R_cl(j) and change of
! position conditions)
! Although it is moved in position, a merged clump is always assumed to be "relevant"
! The merged radius is obtained based on this argument : the density of all clumps
! is the same so (M1+M2)/Rnew^3 = M1/R1^3 = M2/R2^3 and M2/M1 = (R2/R1)^3.
! BEWARE : this subroutine does not fully work yet, as some suspicious snapshots
! of tools/type_merge/mega_merge.gif show (there are false negative). Check associated
! par file.
! -----------------------------------------------------------------------------------
subroutine merge_clumps(Ncl,pos_cl,R_cl,dens_cl)
use glbl_prmtrs
use mod_wind
use miscellaneous
integer, intent(inout) :: Ncl
double precision, allocatable, intent(inout) :: pos_cl(:,:), R_cl(:), dens_cl(:)
integer :: pairs_of_mergers(Ncl)
integer :: dNcl=0, i, j, ii, jj, ii_tmp
double precision :: d
logical :: deleted(Ncl)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

deleted=.false.
pairs_of_mergers=0
do i=1,Ncl-1
  do j=i+1,Ncl ! Check all the clumps "ahead" only since merging is commutative
    ! Compute distance w/ any other clump
    ! CHEAT
    d=dsqrt((pos_cl(i,1)*dsin(pos_cl(i,2))*dcos(pos_cl(i,3))&
            -pos_cl(j,1)*dsin(pos_cl(j,2))*dcos(pos_cl(j,3)))**2.d0+&
            (pos_cl(i,1)*dsin(pos_cl(i,2))*dsin(pos_cl(i,3))&
            -pos_cl(j,1)*dsin(pos_cl(j,2))*dsin(pos_cl(j,3)))**2.d0+&
            (pos_cl(i,1)*dcos(pos_cl(i,2))&
            -pos_cl(j,1)*dcos(pos_cl(j,2)))**2.d0)
    ! If smaller than sum of radii, merge
    if (d<R_cl(i)+R_cl(j)) then
      print*, i, j
      ! By convention, the merged clump is stored into # i
      ! while clump # j is deleted
      if (deleted(i) .and. deleted(j)) then ! beware, this clump has already been merged! => > 2 clumps merger
        cycle
      elseif (deleted(j)) then ! beware, this clump has already been merged! => > 2 clumps merger
        ! w/ a clump k (=pairs_of_mergers(j)) < i => put all in k
        jj=pairs_of_mergers(j)
        ! To find the root one in which they are all stored
        do while(deleted(jj))
          jj=pairs_of_mergers(jj)
        enddo
        ! print*, 'j', jj, deleted(jj)
        ! if (deleted(jj)) call crash("Pb with deleted(j) in merge clumps")
        deleted(i)=.true.
        pairs_of_mergers(i)=jj
        ii=ii_tmp
        ii=jj
        jj=ii_tmp
      elseif (deleted(i)) then ! beware, this clump has already been merged! => > 2 clumps merger
        ii=pairs_of_mergers(i)
        ! print*, 'i', ii, deleted(ii)
        ! To find the root one in which they are all stored
        do while(deleted(ii))
          ii=pairs_of_mergers(ii)
        enddo
        ! if (deleted(ii)) call crash("Pb with deleted(i) in merge clumps")
        deleted(j)=.true.
        pairs_of_mergers(j)=ii
        jj=j
      else
        pairs_of_mergers(i)=j ! j merges w/ i
        pairs_of_mergers(j)=i ! i merges w/ j
        deleted(j)=.true.
        ii=i
        jj=j
      endif
      select case(type_merge_)
      case('volume_x_2')
        R_cl(ii)=R_cl(ii)*(1.d0+(R_cl(jj)/R_cl(ii))**3.d0)**(1.d0/3.d0) ! R_cl(ii)*2.d0**(1.d0/3.d0)
      case('density_x_2')
        dens_cl(ii)=dens_cl(ii)*2.d0
      case default
        call crash('unknown merging type')
      end select
      pos_cl(ii,:)=0.5d0*(pos_cl(ii,:)+pos_cl(jj,:))
      dNcl=dNcl+1
    endif
  enddo
enddo

print*, '- - - '

call delete_clumps(deleted,Ncl,pos_cl,R_cl,dens_cl)
! if (Ncl/=size(pairs_of_mergers)) then
!    print*, pairs_of_mergers
!   print*, 'merger'
! endif

end subroutine merge_clumps
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
! Beware, the initialization assumes no merger and so should be different
! from the stochastic steady-state which is reached later on.
! -----------------------------------------------------------------------------------
subroutine set_ini_radii(Ncl,pos_cl,R_cl,dens_cl)
use glbl_prmtrs
use rdm
use mod_wind
use miscellaneous
integer, intent(in) :: Ncl
double precision, intent(out) :: pos_cl(Ncl,3), R_cl(Ncl), dens_cl(Ncl)
double precision :: pos_cl_old(Ncl,3), R_cl_old(Ncl), dens_cl_old(Ncl)
double precision :: r1, v1
integer :: i
! speed @ 2 stellar radii, the reference point @ which
! clump_rad and clump_dens are given
double precision :: v2strrad
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

call get_v(2.d0,v2strrad)

do i=1,Ncl
  call get_one_over_r2v(r1)
  pos_cl(i,1)=r1
  ! print*, i, r1
  ! write(7,*), r1
  call get_v(r1,v1)
  select case(rad_evol_)
  case('lorenzo')
    R_cl(i)=clump_rad_*(r1/2.d0)**(2.d0/3.d0)*(v1/v2strrad)**(1.d0/3.d0)
    dens_cl(i)=clump_dens_ / ((r1/2.d0)**2.d0*(v1/v2strrad))
  case('jon') ! density evolution such as constant clump mass enforced
    R_cl(i)=clump_rad_*(r1/2.d0)
    dens_cl(i)=clump_dens_ / (r1/2.d0)**3.d0
  case default
    call crash('Which evolution of clump radii with distance to star?')
  end select

enddo

pos_cl_old(:,2)=pos_cl(:,2)
pos_cl_old(:,3)=pos_cl(:,3)
! since the clump properties set in the par file are given @ 2 stellar radii
pos_cl_old(:,1)=2.d0

call expand_clumps(Ncl,pos_cl_old,pos_cl,R_cl,dens_cl)

end subroutine set_ini_radii
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Compute an underestimate of the # of clumps in a spherical shell between
! the stellar radius (r=1) and the distance beyond which the clumps are
! deleted (r=dist_max_cl_). Underestimate since the speed is < vinf (=1)
! but I do not manage to integrate 1/(1-1/r)**beta @ the r=1 edge...
! so I artificially compute Ncl_max_ as 10 times this estimate.
! -----------------------------------------------------------------------------------
subroutine get_clump_number
use miscellaneous
use glbl_prmtrs
use rdm
character(LEN=11) :: string
double precision :: Dt ! dx, ddx, ddt, x, q, v
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Numerical integration to evaluate Dt, the time required to cross
! the distance between the stellar surface and dist_max_cl_
call num_int_steps(100000,'one_over_v',1.d0,dist_max_cl_,Dt,'log')

! Beware, int of a very large double produces a negative integer...
! hence the if condition on the product and not on Ncl_max_
Ncl_max_ = int ( Ndot_ * Dt )

if (Ndot_ * Dt < 1.d9) then
  write (string, "(I10)") Ncl_max_
else
  call crash('More than a billion clumps... are you sure you want to proceed?')
endif

call followup("The # of clumps estimated in the simulation space is ~ "//trim(string))

Ncl0_ = Ncl_max_

! Ncl_max_ = 10 * Ncl_max_

! if (Ncl_max_>Ncl_pointer_) call crash("Pointer on object 'clumps' too small to contain all the clumps")

end subroutine get_clump_number
! -----------------------------------------------------------------------------------

end module mod_clumps
