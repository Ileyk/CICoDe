module mod_physicaldata
   implicit none
   save

   integer :: Ncl_up=10

   type clump
      !> ID of a clump
      integer :: icl=-1
      !> 3D position of each clump
      double precision, dimension(:), allocatable :: x
      !> Radius of each clump
      double precision :: rad
      !> Density of each clump
      double precision :: rho
   end type clump

   !> array of physical states for all blocks on my processor
   type(clump), dimension(:), allocatable, target :: pcl

end module mod_physicaldata
