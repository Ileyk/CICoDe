! Compilation - - - - - - -
! see makefile

! Normalization
! - length normalized to the stellar radius
! - speed to the speed @ infinity
! - mass rate to the stellar mass loss rate

program main
use IO
use glbl_prmtrs
use miscellaneous
use mod_cicode
use mod_init

call cpu_time(chrono_0)

call cpu_time(chrono_1)

!> to index output files (eg log file where the followup subroutine prints)
call give_filenames

call read_arguments()
call read_par_files()

call initialization

call chrono(chrono_1,chrono_1_mess)

print*, 'Core mod starts - - - - - - - - - -'

call cicode

print*, 'Terminated - - - - - - - - - - - -'

call chrono(chrono_0,chrono_0_mess)

end program main
