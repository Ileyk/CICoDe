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
  integer :: Nphases, Ncl0
  double precision :: beta, a, Mdot, clump_mass, clump_rad
  double precision :: mass_fraction, Rstar, vinf, time_max, Per, dist_max_cl
  logical :: deterministic, rlvnt_clumps, do_merge
  character(LEN=400) :: type_merge, rad_evol

  integer :: unitpar=9
  logical            :: file_exists

  namelist /savelist/ Nphases

  namelist /endlist/ time_max

  namelist /star/ Rstar

  namelist /wind/ beta, vinf, Mdot

  namelist /clumps/ Ncl0, clump_mass, clump_rad, mass_fraction, dist_max_cl, &
    do_merge, type_merge, rad_evol

  namelist /orbit/ a, Per

  namelist /numerics/ deterministic, rlvnt_clumps

  open(3,file=trim(err_fl))

  Nphases = 32

  time_max = 1.d0

  Rstar = 20.d0

  beta = 0.5d0
  vinf = 1.d2
  Mdot = 1.d-6

  Ncl0          = 1000
  clump_mass    = 1.d17
  clump_rad     = 0.01d0
  mass_fraction = 0.1d0
  dist_max_cl   = 5.d0
  do_merge      = .false.
  type_merge    = undefined
  rad_evol      = 'lorenzo'

  a = 1.8d0
  Per = 9.d0

  deterministic = .true.
  rlvnt_clumps  = .false.

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
       read(unitpar, numerics, end=107)

107  close(unitpar)

  ! To do also : compute a preliminary likelihood of "runaway" merger
  ! if type_merge='volume_x_2'
  ! CHEAT
  if (do_merge) call crash('Not fully functional yet : are you sure?')

  if (do_merge .and. (type_merge==undefined)) &
    call crash('Which type of merger?')

  if ((.not. do_merge) .and. (type_merge/=undefined)) &
    call crash('Why did you specify a type of merger?')

  Nphases_=Nphases

  time_max_=time_max

  Rstar_=Rstar

  beta_=beta
  vinf_=vinf
  Mdot_=Mdot

  Ncl0_          = Ncl0
  clump_mass_    = clump_mass
  clump_rad_     = clump_rad
  mass_fraction_ = mass_fraction
  dist_max_cl_   = dist_max_cl
  do_merge_      = do_merge
  type_merge_    = type_merge
  rad_evol_      = rad_evol

  a_   = a
  Per_ = Per

  deterministic_=deterministic
  rlvnt_clumps_ =rlvnt_clumps

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

! Since lines are appended in log_file each time (see in miscellaneous.f90),
! we need to erase the previous one first
call system("rm -f "//log_fl)
! et tant qu a faire...
call system("rm -f "//err_fl)
call system("rm -f "//pos_fl)
call system("rm -f "//dis_fl)

end subroutine give_filenames
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
subroutine save_pos(Ncl,pos_cl,R_cl)
use glbl_prmtrs
integer, intent(in) :: Ncl
double precision, intent(in) :: pos_cl(Ncl,3), R_cl(Ncl)
integer :: i
logical :: file_exists
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

INQUIRE(FILE=pos_fl, EXIST=file_exists)
if (.not. file_exists) then
  open(unit=1,file=pos_fl)
  write(1,'(a)') 'x | y | R'
  write(1,'((a),(2I12))') 'Nphases Nclumps0 ', Nphases_, Ncl0_
  close(1)
endif
open(unit=1,file=pos_fl,access='append')
do i=1, Ncl
  write(1,'(200(1pe12.4))') pos_cl(i,1)*dcos(pos_cl(i,3)), pos_cl(i,1)*dsin(pos_cl(i,3)), R_cl(i)
enddo
write(1,'(a)') 'xxx'
close(1)

end subroutine save_pos
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
