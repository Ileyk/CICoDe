module mod_func
  use glbl_prmtrs
contains

! -----------------------------------------------------------------------------------
function v_wind (x)
double precision :: v_wind
double precision, intent (in) :: x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
v_wind = (1.d0-1.d0/x)**beta_
return
end function v_wind
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Same v_wind as above but w/ Rstar_ divided by rmax_plot_
! as input parameter.
! -----------------------------------------------------------------------------------
function v_wind_prm (x)
double precision :: v_wind_prm
double precision, intent (in) :: x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
v_wind_prm = (1.d0-(Rstar_/rmax_plot_)/x)**beta_
return
end function v_wind_prm
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
function one_over_vr2 (x)
double precision :: one_over_vr2
double precision, intent (in) :: x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
one_over_vr2 = 1.d0 / (x**2.d0*v_wind(x))
return
end function one_over_vr2
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
function one_over_v (x)
double precision :: one_over_v
double precision, intent (in) :: x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
one_over_v = 1.d0 / v_wind(x)
return
end function one_over_v
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Gives the radius of a clump at radial position r, w/ clump_rad_ the reference set
! at 2 stellar radii, according to Lorenzo's formula.
! -----------------------------------------------------------------------------------
function lorenzo_rad (x)
double precision :: lorenzo_rad
double precision, intent (in) :: x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
lorenzo_rad = clump_rad_ * (x/2.d0)**(2.d0/3.d0) * (v_wind(x)/v_wind(2.d0))**(1.d0/3.d0)
return
end function lorenzo_rad
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Gives the radius of a clump at radial position r, w/ clump_rad_ the reference set
! at 2 stellar radii, according to Jon's formula.
! -----------------------------------------------------------------------------------
function jon_rad (x)
double precision :: jon_rad
double precision, intent (in) :: x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
jon_rad = clump_rad_ * (x/2.d0)
return
end function jon_rad
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Gives the density of a clump at radial position r, w/ clump_rad_ the reference set
! at 2 stellar radii, according to Lorenzo's formula.
! -----------------------------------------------------------------------------------
function lorenzo_rho (x)
double precision :: lorenzo_rho
double precision, intent (in) :: x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
lorenzo_rho = clump_dens_ / ((x/2.d0)**2.d0*(v_wind(x)/v_wind(2.d0)))
return
end function lorenzo_rho
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Gives the density of a clump at radial position r, w/ clump_rad_ the reference set
! at 2 stellar radii, according to Jon's formula.
! -----------------------------------------------------------------------------------
function jon_rho (x)
double precision :: jon_rho
double precision, intent (in) :: x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
jon_rho = clump_dens_ / (x/2.d0)**3.d0
return
end function jon_rho
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Gives the radial profile of the number density of the clumps
! -----------------------------------------------------------------------------------
function clump_nbr_dens (x)
double precision :: clump_nbr_dens
double precision, intent (in) :: x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
clump_nbr_dens = Ndot_ / (4.d0*dpi*x**2.d0*v_wind(x))
return
end function clump_nbr_dens
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Gives the radial profile of the mean distance between the clumps centers
! -----------------------------------------------------------------------------------
function mean_dist (x)
double precision :: mean_dist
double precision, intent (in) :: x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
mean_dist = clump_nbr_dens(x)**(-1.d0/3.d0)
return
end function mean_dist
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
function  r2_over_v(x)
double precision :: r2_over_v, x_star
double precision, intent (in) :: x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x_star=Rstar_/rmax_plot_
r2_over_v = x**2.d0 / v_wind_prm(x)
return
end function r2_over_v
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
function  r2_over_v_arcsin(x)
use miscellaneous
double precision :: r2_over_v_arcsin
double precision, intent (in) :: x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (x<1.d0) call crash('Function not defined on this range.')
r2_over_v_arcsin = r2_over_v(x) * dsin(0.5d0*dasin(1.d0/x))**2.d0
return
end function r2_over_v_arcsin
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
function  r2v_23_over_v(x)
double precision :: r2v_23_over_v, x_star
double precision, intent (in) :: x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
x_star=Rstar_/rmax_plot_
r2v_23_over_v = (v_wind_prm(x)/v_wind(2.d0))**(2.d0/3.d0) * (x**2.d0)**(2.d0/3.d0) / v_wind_prm (x)
return
end function r2v_23_over_v
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
function  r2v_23_over_v_arcsin(x)
use miscellaneous
double precision :: r2v_23_over_v_arcsin
double precision, intent (in) :: x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (x<1.d0) call crash('Function not defined on this range.')
r2v_23_over_v_arcsin = r2v_23_over_v(x) * dsin(0.5d0*dasin(1.d0/x))**2.d0
return
end function r2v_23_over_v_arcsin
! -----------------------------------------------------------------------------------

end module mod_func
