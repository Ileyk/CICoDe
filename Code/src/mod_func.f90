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
function one_over_vr2 (x)
double precision :: one_over_vr2
double precision, intent (in) :: x
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
one_over_vr2 = 1.d0 / (x**2.d0*v_wind(x)) ! x**3.d0-x**2.d0-5.d0*x+3.d0 ! 8.0 * x
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

end module mod_func
