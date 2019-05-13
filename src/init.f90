module mod_init
  implicit none

contains

! -----------------------------------------------------------------------------------
!> Initialize global variables
subroutine initialization()
use glbl_prmtrs
use rdm
use miscellaneous
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Dimension
Rstar_ = Rstar_ * rsolar ! cm
vinf_  = vinf_ * 1.d5 ! cm/s
Mdot_  = Mdot_ * msun / secinyear
clump_mass_ = clump_mass_ / (vinf_/(Mdot_*Rstar_))

Per_=Per_*secinday/(Rstar_/vinf_) ! orbital period in code units
dt_=Per_/dble(Nphases_)

if (time_max_/=1.d0) call crash('Pas certain encore que ça puisse se faire ça...')
time_max_=time_max_*Per_

Rstar_ = Rstar_ / Rstar_
vinf_  = vinf_  / vinf_
Mdot_  = Mdot_  / Mdot_

! Compute initialization radius of clumps,
! set as where v_beta = 1% of terminal wind speed
rini_=1.d0/(1.d0-0.01d0**(1.d0/beta_))

! Assuming spherical clumps
clump_dens_=3.d0*clump_mass_/(4.d0*dpi*clump_rad_**3.d0)

if (.not. deterministic_) call init_random_seed

end subroutine initialization
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
! Pick up clumps @ random in angles (th,ph) but also in r, such as conservation of
! # of clumps ie nvr2=cst w/ n # density of clumps.
! -----------------------------------------------------------------------------------
subroutine set_ini_clumps(Ncl,pos_cl,R_cl,dens_cl)
use glbl_prmtrs
use mod_clumps
use rdm
use mod_geometry
integer, intent(in) :: Ncl
double precision, intent(inout) :: pos_cl(Ncl,3), R_cl(Ncl), dens_cl(Ncl)
logical :: rlvnt=.true.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Add all the Ncl clumps
! CHEAT
! call add_clumps(Ncl,Ncl,pos_cl,R_cl,dens_cl)
integer :: i
double precision :: r, th, ph
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
  R_cl(i)=clump_rad_
  dens_cl(i)=clump_dens_
enddo

! CHEAT For now, we simply pick up r @ random between rini & 10. ...
! But overwrites the radial distances to stellar center
! call set_ini_radii(Ncl,pos_cl,R_cl,dens_cl)
do i=1,Ncl
  call get_flat(rini_,10.d0,r)
  pos_cl(i,1)=r
  write(7,*), r
enddo

end subroutine set_ini_clumps
! -----------------------------------------------------------------------------------


end module mod_init
