program main
use mod_dynamics_X
use mod_physicaldata

integer :: max_clumps=1000, icl

! Just after having read the par file
! Equivalent to initialize_vars in initialize_amrvac in mod_initialize.t
allocate(pcl(max_clumps))

! Equivalent to initlevelone in amrini.t
do icl=1,Ncl_up
  call alloc_node(icl)
  call something(pcl(icl)%x)
  if (icl==7) call dealloc_node(icl)
enddo


end program main
