module mod_init
  implicit none

contains

! -----------------------------------------------------------------------------------
!> Initialize global variables
subroutine initialization()
use glbl_prmtrs
use rdm
use miscellaneous
use mod_clumps
use mod_wind
! speed @ 2 stellar radii, the reference point @ which
! clump_rad and clump_dens are given
double precision :: v2strrad
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Dimension
Rstar_ = Rstar_ * rsolar ! cm
vinf_  = vinf_ * 1.d5 ! cm/s
Mdot_  = Mdot_ * msun / secinyear
clump_mass_ = clump_mass_ / (Mdot_*Rstar_/vinf_)

Per_=Per_*secinday/(Rstar_/vinf_) ! orbital period in code units
a_  =a_*Rstar_    /Rstar_

! if (time_max_/=1.d0) call crash('Pas certain encore que ça puisse se faire ça...')
time_max_=time_max_*Per_

! CHEAT
! dt_=Per_/dble(Nphases_)
dt_ = time_max_/dble(Nphases_)

Rstar_ = Rstar_ / Rstar_
vinf_  = vinf_  / vinf_
Mdot_  = Mdot_  / Mdot_

! # of clumps flowing through any sphere per unit time
Ndot_ = mass_fraction_ * Mdot_ / clump_mass_

! Assuming spherical clumps
clump_dens_=3.d0*clump_mass_/(4.d0*dpi*clump_rad_**3.d0)

! Compute initialization radial position rini of clumps,
! set as where v_beta = 1% of terminal wind speed
vini_=0.01d0
rini_=1.d0/(1.d0-vini_**(1.d0/beta_))

if (dabs(rini_-1.d0)<1.d-12) call crash("Velocity profile suspiciously steep, beware")

! Speed @ 2 stellar radii, the reference position for clump_rad and clump_dens
call get_v(2.d0,v2strrad)

! Deduce clump radius and density @ rini
select case (rad_evol_)
  case('lorenzo')
    clump_rad_ini_  = clump_rad_  *  (rini_/2.d0)**(2.d0/3.d0)*(vini_/v2strrad)**(1.d0/3.d0) ! clump_rad is @ 2 stellar radii by convention
    clump_dens_ini_ = clump_dens_ / ((rini_/2.d0)**2.d0*(vini_/v2strrad)) ! idem
  case('jon')
    clump_rad_ini_  = clump_rad_  * (rini_/2.d0) ! clump_rad is @ 2 stellar radii by convention
    clump_dens_ini_ = clump_dens_ / (rini_/2.d0)**3.d0 ! idem
  case default
    call crash('No evolution of clump radius specified')
end select

if (.not. deterministic_) call init_random_seed

! Preliminary estimate of the # of clumps in the simulation space
call get_clump_number

! It is where we store the fraction of clumps, to make sure we do not have
! to add a # of clumps > 1 @ each dynamical time scale
cmltd_clump_=0.d0

end subroutine initialization
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
! Pick up clumps @ random in angles (th,ph) and pick up r
! such as conservation of # of clumps ie nvr2=cst w/ n # density of clumps.
! -----------------------------------------------------------------------------------
subroutine set_ini_clumps(Ncl,pos_cl,R_cl,dens_cl)
use glbl_prmtrs
use mod_clumps
use rdm
use mod_geometry
use mod_wind
use mod_func
integer, intent(in) :: Ncl
double precision, intent(inout) :: pos_cl(Ncl,3), R_cl(Ncl), dens_cl(Ncl)
logical :: rlvnt=.true.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Add all the Ncl clumps
! CHEAT
! call add_clumps(Ncl,Ncl,pos_cl,R_cl,dens_cl)
integer :: i
double precision :: th, ph
double precision :: r1, v1

double precision :: r(Ncl)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

do i=1,Ncl
  ! Pick up @ random angular position (over dphi and d(cos(th)) )
  call get_flat(-1.d0,1.d0,th)
  if (rlvnt_clumps_) then
    call is_relevant(th,rlvnt)
    ! CHEAT
    ! if (.not. rlvnt) PICK UP NEW TH AND PH
  endif
  call get_flat(0.d0,2.d0*dpi,ph)
  ! write(7,*), th, ph
  ! CHEAT to debug
  pos_cl(i,2)=dpi/2.d0 ! dacos(th)
  pos_cl(i,3)=ph
enddo

r=pos_cl(:,1)
call accptnce_rjctnce(Ncl,r,R_cl,dens_cl)
pos_cl(:,1)=r

end subroutine set_ini_clumps
! -----------------------------------------------------------------------------------



end module mod_init
