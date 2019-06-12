module IO
  implicit none
  public

contains

! -----------------------------------------------------------------------------------
!> Read the command line arguments passed to amrvac
subroutine read_arguments()
  use mod_kracken
  use glbl_prmtrs

  integer                          :: len, ier
  integer                          :: ibegin
  integer                          :: iterm

  ! Specify the options and their default values
  call kracken('cmd','-i amrvac.par -if ' // undefined // &
       ' -slice 0 -collapse 0 --help .false. -convert .false.')

  ! Get the par file(s)
  call retrev('cmd_i', prmtr_fl, len, ier)

  ! Split the input files, in case multiple were given
  call get_fields_string(prmtr_fl, " ,'"""//char(9))

end subroutine read_arguments
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
!> Read in the user-supplied parameter-file
subroutine read_par_files()
  use miscellaneous
  ! use draw ! to determine position of L1
  use glbl_prmtrs
  integer :: Nsave, restart_indx
  double precision :: beta, a, Mdot, clump_mass, clump_rad
  double precision :: mass_fraction, Rstar, vinf, time_max, Per, dist_max_cl, NSspin
  double precision :: rmax_plot, inclnsn, omega_i
  logical :: deterministic, do_merge, init_only, plot_cl_only
  character(LEN=400) :: type_merge, rad_evol

  integer :: unitpar=9
  logical            :: file_exists

  namelist /savelist/ Nsave, restart_indx

  namelist /endlist/ time_max

  namelist /star/ Rstar

  namelist /wind/ beta, vinf, Mdot

  namelist /clumps/ clump_mass, clump_rad, mass_fraction, dist_max_cl, &
    do_merge, type_merge, rad_evol

  namelist /orbit/ a, Per, NSspin

  namelist /LOS/ inclnsn, omega_i

  namelist /numerics/ deterministic, init_only, plot_cl_only

  namelist /plot/ rmax_plot

  open(3,file=trim(err_fl))

  Nsave = 32
  restart_indx = -1 ! by default, no restart

  time_max = 1.d0

  Rstar = 20.d0

  beta = 0.5d0
  vinf = 1.d2
  Mdot = 1.d-6

  ! Ncl0          = 1000 ! DEPRECIATED
  clump_mass    = 1.d17
  clump_rad     = 0.01d0
  mass_fraction = 0.1d0
  dist_max_cl   = 5.d0
  do_merge      = .false.
  type_merge    = undefined
  rad_evol      = 'lorenzo'

  a = 1.8d0
  Per = 9.d0
  NSspin = 300.d0

  ! default is edge-on
  inclnsn = 90.d0
  omega_i =  0.0

  rmax_plot = -1.d0

  deterministic = .true.
  init_only     = .false.
  ! Parameter to bypass the computation of NH and get only the
  ! positions of "all" the clumps (=> do not delete those useless for NH).
  ! To produce fancy animated GIF w/ positions_3D.py.
  plot_cl_only  = .false.

  print *, "Reading " // trim(prmtr_fl)

  ! Check whether the file exists
  inquire(file=trim(prmtr_fl), exist=file_exists)

  if (.not. file_exists) &
    call crash("The parameter file " // trim(prmtr_fl) // " does not exist")

  open(unitpar, file=trim(prmtr_fl), status='old')
  ! Try to read in the namelists. They can be absent or in a different
  ! order, since we rewind before each read.
  rewind(unitpar)
       rewind(unitpar)
       read(unitpar, savelist, end=101)

101    rewind(unitpar)
       read(unitpar, endlist, end=102)

102    rewind(unitpar)
       read(unitpar, star, end=103)

103    rewind(unitpar)
       read(unitpar, wind, end=104)

104    rewind(unitpar)
       read(unitpar, clumps, end=105)

105    rewind(unitpar)
       read(unitpar, orbit, end=106)

106    rewind(unitpar)
       read(unitpar, LOS, end=107)

107    rewind(unitpar)
       read(unitpar, numerics, end=108)

108    rewind(unitpar)
       read(unitpar, plot, end=109)

109  close(unitpar)

  ! DEPRECIATED : now, Ncl0 is computed based on Ndot_ and the time to cross the distance
  ! between rini_ and dist_max_cl_.
  ! Ncl0=0 => dt_dyn ini too big => spurious "shell of clumps" being sent @ the beginning
  ! (OK once they leave the simulation space)
  ! if (Ncl0==0) call crash('Possible but beware, the first dt_dyn is meaningless.')

  ! To do also : compute a preliminary likelihood of "runaway" merger
  ! if type_merge='volume_x_2'
  ! CHEAT
  if (do_merge) call crash('Not fully functional yet : are you sure?')

  if (do_merge .and. (type_merge==undefined)) &
    call crash('Which type of merger?')

  if ((.not. do_merge) .and. (type_merge/=undefined)) &
    call crash('Why did you specify a type of merger?')

  if (Nsave>9999) call crash('Too many saved files')

  if (init_only .and. restart_indx>-1) &
    call crash('No restart possible if you want just the initialization.')

  if (inclnsn<0.d0) &
    call crash('Beware, inclination < 0 => initial time @ superior instead of inferior conjunction. Risky?')

  ! default value
  if (rmax_plot<0.d0) &
    rmax_plot     = 1.2*a

  if (rmax_plot<1.d0 .or. rmax_plot>dist_max_cl) &
    call crash('Distance up to which plot w/ positions_3D.py too small or too large.')


  ! Convert from degrees to radians
  inclnsn = inclnsn * dpi/180.d0
  omega_i = omega_i * dpi/180.d0

  Nsave_=Nsave
  restart_indx_=restart_indx

  time_max_=time_max

  Rstar_=Rstar

  beta_=beta
  vinf_=vinf
  Mdot_=Mdot

  ! Ncl0_          = Ncl0
  clump_mass_    = clump_mass
  clump_rad_     = clump_rad
  mass_fraction_ = mass_fraction
  dist_max_cl_   = dist_max_cl
  do_merge_      = do_merge
  type_merge_    = type_merge
  rad_evol_      = rad_evol

  a_      = a
  Per_    = Per
  NSspin_ = NSspin

  inclnsn_ = inclnsn
  omega_i_ = omega_i

  deterministic_=deterministic
  init_only_    =init_only
  plot_cl_only_ =plot_cl_only

  rmax_plot_    = rmax_plot

end subroutine read_par_files
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine give_filenames
use glbl_prmtrs

fldr='./output/'

! adim_file=trim(fldr)//'parameters'
! param_file=trim(fldr)//'param' ! where you store the 4 degrees of freedom
!
! ! fldrs
pos_fl=trim(fldr)//'positions'
dis_fl=trim(fldr)//'initial_distributions'
prsty_fl=trim(fldr)//'porosity'
NH_fl=trim(fldr)//'NH'
NH_shells_fl=trim(fldr)//'NH_shells'
posX_fl=trim(fldr)//'posX'
orb_fl=trim(fldr)//'orbit'

log_fl=trim(fldr)//'log'
! scale_file=trim(fldr)//'scale'
! ! par_file=trim(fldr)//'par'
! big_picture=trim(fldr)//'rk4_traj'
! surf_pts=trim(fldr)//'pts'
! arrival_map=trim(fldr)//'map'
out_fl=trim(fldr)//'out'
err_fl=trim(fldr)//'err'
! hist_file=trim(fldr)//'vel_at_ROCHE_surf_hist'
! mass_wind_1D=trim(fldr)//'mass_wind_1D'
! theta_logangle=trim(fldr)//'theta_mesh'
!
! draw_all_in_plane=trim(fldr)//'draw_all_in_plane'
! draw_good_in_plane=trim(fldr)//'draw_good_in_plane'

! resultats=trim(fldr)//'results'

! the updated version of out, directuly used in dimensionizing.f90 in vualatv
! dim_file =trim(fldr)//'dim'  ! where you store the dimensioned computed quantities (but still in (GM2,roche1,L) units)

end subroutine give_filenames
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine clean_outputs
use glbl_prmtrs

! Since lines are appended in log_file each time (see in miscellaneous.f90),
! we need to erase the previous one first
call system("rm -f "//log_fl)
! et tant qu a faire...
call system("rm -f "//err_fl)
call system("rm -f "//pos_fl)
call system("rm -f "//dis_fl)
call system("rm -f "//prsty_fl)
call system("rm -f "//NH_fl)
call system("rm -f "//NH_shells_fl)
call system("rm -f "//posX_fl)
call system("rm -f "//orb_fl)

call system("rm -f "//pos_fl//"_*.dat")

end subroutine clean_outputs
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! Save the radial profile of clump radius and mean distance between centers.
! Deduce the profile of the distance between clump shells compared to clump radii.
! -----------------------------------------------------------------------------------
subroutine save_prsty_prfle
use glbl_prmtrs
use mod_func
integer :: i, Nstep=1000
double precision :: x, dx, q, xmin, xmax
logical :: file_exists
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

INQUIRE(FILE=prsty_fl, EXIST=file_exists)
if (.not. file_exists) then
  open(unit=1,file=prsty_fl)
  write(1,'(a)') 'r | d (mean dist.) | Rcl | d/Rcl - 2'
  close(1)
endif

open(unit=1,file=prsty_fl,access='append')

xmin=rini_
xmax=dist_max_cl_
! Stretching factor
q   = (xmax/xmin)**(1.d0/dble(Nstep))
dx  = xmin * (q-1.d0)
x   = xmin + dx/2.d0 ! centers still half-way in-between the edges (reverse not true)
do i=1,Nstep
  select case(rad_evol_)
  case('lorenzo')
    write(1,'(200(1pe12.4))') x, mean_dist(x), lorenzo_rad(x), mean_dist(x)/lorenzo_rad(x) - 2.d0
  case('jon')
    write(1,'(200(1pe12.4))') x, mean_dist(x), jon_rad(x), mean_dist(x)/jon_rad(x) - 2.d0
  end select
  x  =  x * q
  dx = dx * q
enddo

close(1)

end subroutine save_prsty_prfle
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine save_pos(Ncl,pos_cl,R_cl,dens_cl)
use glbl_prmtrs
use mod_dynamics_X
integer, intent(in) :: Ncl
double precision, intent(in) :: pos_cl(Ncl,3), R_cl(Ncl), dens_cl(Ncl)
integer :: i
logical :: file_exists
character(LEN=80) :: filename
double precision :: phase, pos_X(1:3)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write(filename,'(I4.4)') save_index_
filename = trim(pos_fl) // '_' // trim(filename) // '.dat'

! Write header
! Save time & # of clumps for restart
open(unit=1,file=trim(filename))
write(1,'((d20.14),(I10))') t_, Ncl
write(1,'(a)') '(x,y,z) position of the X-ray source'
phase=2.d0*dpi*(t_/Per_)
call get_pos_X(phase,pos_X)
write(1,'(200(e12.4))') (pos_X(i), i=1,3)
write(1,'(a)') 'x | y | z | R | rho'
close(1)

! INQUIRE(FILE=pos_fl, EXIST=file_exists)
! if (.not. file_exists) then
!   open(unit=1,file=pos_fl)
!   write(1,'(a)') 'x | y | z | R'
!   write(1,'((a),(2I12))') 'Nsave ', Nsave_
!   close(1)
! endif

! Write data
! Save Cartesian coordinates of clumps
open(unit=1,file=trim(filename),access='append')
do i=1, Ncl
  write(1,'(200(1pe12.4))') pos_cl(i,1)*dsin(pos_cl(i,2))*dcos(pos_cl(i,3)), &
                            pos_cl(i,1)*dsin(pos_cl(i,2))*dsin(pos_cl(i,3)), &
                            pos_cl(i,1)*dcos(pos_cl(i,2)), R_cl(i), dens_cl(i)
enddo
close(1)

save_index_ = save_index_ + 1

end subroutine save_pos
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
! If it is a restart, we load the positions, radii and densities of the clumps
! from the data file. 3 reasons why there is a drift of the result which end up
! being quite different with and without restart :
!   - float in saved and loaded file /= double precision
!   - even if random seed is the same (ie deterministic=T), not the same pick up
!   - cmltd_clump_ is not saved, we start again from cmltd_clump_ = 0
! -----------------------------------------------------------------------------------
subroutine load_clumps(Ncl,pos_cl,R_cl,dens_cl)
use glbl_prmtrs
use mod_cart_sph
use miscellaneous
integer, intent(out) :: Ncl
double precision, allocatable, intent(out) :: pos_cl(:,:), R_cl(:), dens_cl(:)
integer :: i
logical :: file_exists
double precision :: x(1:3), r(1:3)
character(LEN=80) :: filename
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

write(filename,'(I4.4)') restart_indx_
filename = trim(pos_fl) // '_' // trim(filename) // '.dat'

INQUIRE(FILE=filename, EXIST=file_exists)
if (.not. file_exists) then
  call crash('File '//trim(filename)//' does not exist, no restart from it possible.')
endif

! Read time & # of clumps from header
! Save time & # of clumps for restart
open(unit=1,file=trim(filename))
read(1,'((d20.14),(I10))') t0_, Ncl
! Skip header
read(1,*)
read(1,*)
read(1,*)

allocate(pos_cl(Ncl,3), R_cl(Ncl), dens_cl(Ncl))

! Read data
! Save Cartesian coordinates of clumps
do i=1, Ncl
  read(1,'(200(1pe12.4))') x(1), x(2), x(3), R_cl(i), dens_cl(i)
  call cart_to_sph(x)
  pos_cl(i,1)=x(1)
  pos_cl(i,2)=x(2)
  pos_cl(i,3)=x(3)
enddo
close(1)

end subroutine load_clumps
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine save_histograms(Ncl,pos_cl,R_cl)
use glbl_prmtrs
integer, intent(in) :: Ncl
double precision, intent(in) :: pos_cl(Ncl,3), R_cl(Ncl)
integer :: i
logical :: file_exists
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

INQUIRE(FILE=dis_fl, EXIST=file_exists)
if (.not. file_exists) then
  open(unit=2,file=dis_fl)
  write(2,'(a)') 'r | th | ph | R'
  write(2,'((a),(2I12))') 'Nclumps0 ', Ncl0_
  close(2)
endif
open(unit=2,file=dis_fl,access='append')
do i=1, Ncl
  write(2,'(200(1pe12.4))') pos_cl(i,1), pos_cl(i,2), pos_cl(i,3), R_cl(i)
enddo
close(2)

end subroutine save_histograms
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine save_NH(NH)
use glbl_prmtrs
double precision, intent(in) :: NH(Nphases_)
integer :: i, j, line_restart
logical :: file_exists
character(LEN=480) :: string, command
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! INQUIRE(FILE=NH_fl, EXIST=file_exists)
! if (.not. file_exists) then
!   open(unit=1,file=NH_fl)
!   write(1,'(a)') 't | NH(ph1) | NH(ph2) | ... '
!   close(1)
! endif

! In case of restart, delete all the lines beyond the appropriate index
! If restart AND first time save_NH is called
if (restart_indx_>0 .and. dabs( t_ - dt_*(restart_indx_*(Nphases_/Nsave_)+1) ) < smalldble) then
  ! 1. Look for the line of restart
  ! a. ... to do so, grep the index just after restart and save the result in a tmp file
  write(string,'(I8.8)') restart_indx_*(Nphases_/Nsave_)+1
  ! write(string,'(I1.1)') 7
  print*, trim(string)
  command = "grep -rni " // trim(string) // " output/NH > output/tmp"
  ! command = "grep -rni " // trim(string) // " test > output/tmp"
  call system(command)
  print*, trim(command)
  ! b. read the result in the tmp file
  open(1,file='output/tmp')
  read(1,'(a)') string
  close(1)
  ! c. extract the appropriate line number, line_restart
  print*, trim(string)
  i=index(string,':') ! separators of grep to delimit the line number
  j=index(string,':',back=.true.)
  read(string(i+1:j-1),'(I8.8)') line_restart
  print*, i, j
  print*, line_restart, string(i+1:j-1)
  ! d. remove all lines after (and including) line_restart
  write(string,'(I8.8)') line_restart
  command="sed -i '' '"//trim(string)//",$d' "//trim(NH_fl)
  ! command="sed -i '' '"//trim(string)//",$d' test"
  print*, trim(command)
  call system(command)
  ! write(1,*) 'pif'
  ! remove the tmp file
  call system("rm output/tmp")
  ! print*,
  ! print*, dt_*restart_indx_*(Nphases_/Nsave_)
endif

open(unit=1,file=NH_fl,access='append')

write(1,'(a,1pe12.4,a,I8.8,a,I8.8)') 'time ', t_, ' | phase bin ', int((t_+smalldble)/dt_), ' in ', Nphases_ !, (NH(i),i=1,Nphases_)
do i=1, Nphases_
  write(1,'(200(1pe12.4))') NH(i)
enddo
write(1,*) ' '

close(1)

end subroutine save_NH
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine save_NH_shells(N_bin,NH_bin)
use glbl_prmtrs
integer, intent(in) :: N_bin
double precision, intent(in) :: NH_bin(N_bin)
double precision :: NH_cmltd
integer :: i
logical :: file_exists
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

INQUIRE(FILE=NH_shells_fl, EXIST=file_exists)
if (.not. file_exists) then
  open(unit=2,file=NH_shells_fl)
  write(2,'(a)') 'r_center | r_left | r_right | NH | cumulated NH'
  close(2)
endif
open(unit=2,file=NH_shells_fl,access='append')
NH_cmltd=0.d0
do i=1, N_bin
  NH_cmltd=NH_cmltd+NH_bin(i)
  write(2,'(200(1pe12.4))') Rstar_*(dist_max_cl_/Rstar_)**((dble(i-1)+0.5d0)/dble(N_bin)), &
                            Rstar_*(dist_max_cl_/Rstar_)**(dble(i-1)        /dble(N_bin)), &
                            Rstar_*(dist_max_cl_/Rstar_)**(dble(i  )        /dble(N_bin)), &
                            NH_bin(i), NH_cmltd
enddo
close(2)

end subroutine save_NH_shells
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
subroutine save_pos_X(pos_X)
use glbl_prmtrs
double precision, intent(in) :: pos_X(3)
integer :: i
logical :: file_exists
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

INQUIRE(FILE=posX_fl, EXIST=file_exists)
if (.not. file_exists) then
  open(unit=1,file=posX_fl)
  write(1,'(a)') 't | x | y | z '
  close(1)
endif
open(unit=1,file=posX_fl,access='append')

write(1,'(200(1pe12.4))') t_, (pos_X(i),i=1,3)

close(1)

end subroutine save_pos_X
! -----------------------------------------------------------------------------------

! -----------------------------------------------------------------------------------
subroutine save_orbit
use glbl_prmtrs
use mod_dynamics_X
double precision :: phase, pos_X(1:3)
integer :: i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! Header
open(unit=1,file=orb_fl)
write(1,'(a)') 'phase | x | y | z '
close(1)

open(unit=1,file=orb_fl)

phase=0.
do i=1,1000
  call get_pos_X(phase,pos_X)
  write(1,'(4(1pe12.4))') phase, pos_X(1), pos_X(2), pos_X(3)
  phase=phase+2.d0*dpi/1000.d0
enddo

close(1)

end subroutine save_orbit
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
!> Routine to find entries in a string
subroutine get_fields_string(line, delims)
  !> The line from which we want to read
  character(len=*), intent(inout) :: line
  !> A string with delimiters. For example delims = " ,'"""//char(9)
  character(len=*), intent(in)    :: delims

  integer :: ixs_start
  integer :: ixs_end
  integer :: ix, ix_prev

  ix_prev = 0

  ! Find the starting point of the next entry (a non-delimiter value)
  ix = verify(line(ix_prev+1:), delims)

  ixs_start = ix_prev + ix ! This is the absolute position in 'line'

  ! Get the end point of the current entry (next delimiter index minus one)
  ix = scan(line(ixs_start+1:), delims) - 1

  if (ix == -1) then              ! If there is no last delimiter,
    ixs_end = len(line) ! the end of the line is the endpoint
  else
    ixs_end = ixs_start + ix
  end if

  line = line(ixs_start:ixs_end)

end subroutine get_fields_string
! -----------------------------------------------------------------------------------

end module IO
