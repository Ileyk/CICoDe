!=============================================================================
subroutine alloc_node(icl)
use mod_physicaldata
integer, intent(in) :: icl
!-----------------------------------------------------------------------------
if(.not. allocated(pcl(icl)%x)) then
  ! allocate arrays for clumps
  call alloc_state(pcl(icl))
endif

pcl(icl)%x=0.d0
pcl(icl)%icl=icl

end subroutine alloc_node
!=============================================================================

!=============================================================================
subroutine alloc_state(s)
use mod_physicaldata
type(clump) :: s
!-----------------------------------------------------------------------------

allocate(s%x(3))
! allocate(s%rad(1:Ncl))
! allocate(s%rho(1:Ncl))

end subroutine alloc_state
!=============================================================================

!=============================================================================
subroutine dealloc_node(icl)
use mod_physicaldata
integer, intent(in) :: icl
!-----------------------------------------------------------------------------
if (icl==0) then
   ! call mpistop("trying to delete a non-existing clump in dealloc_node")
end if
call dealloc_state(pcl(icl))

end subroutine dealloc_node
!=============================================================================

!=============================================================================
subroutine dealloc_state(s)
use mod_physicaldata
type(clump) :: s
!-----------------------------------------------------------------------------

deallocate(s%x)
! deallocate(s%rad)
! deallocate(s%rho)

end subroutine dealloc_state
!=============================================================================
