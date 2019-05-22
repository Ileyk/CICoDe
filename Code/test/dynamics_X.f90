module mod_dynamics_X
contains

! -----------------------------------------------------------------------------------
subroutine something(x)
double precision, intent(inout) :: x(1:3)
integer :: i, j
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

! print*, 'shape at the beginning of something :', shape(pos)

do i=1,3
  x(i)=dble(i)
  ! print*, i, x(i)
enddo

! print*, 'shape at the end of something :', shape(pos)

end subroutine something
! -----------------------------------------------------------------------------------

end module mod_dynamics_X
