module mod_cart_sph
use glbl_prmtrs
contains

! -----------------------------------------------------------------------------------
! Convert from spherical to cartesian the positions of the clumps
! -----------------------------------------------------------------------------------
subroutine sph_2_cart(Ncl,pos_cl,pos_cl_cart)
integer, intent(in) :: Ncl
double precision, intent(in) :: pos_cl(1:Ncl,1:3)
double precision, intent(out) :: pos_cl_cart(1:Ncl,1:3)
double precision :: r, t, p, x, y, z
integer :: i

do i=1,Ncl
  r=pos_cl(i,1)
  t=pos_cl(i,2)
  p=pos_cl(i,3)
  x=r*dsin(t)*dcos(p)
  y=r*dsin(t)*dsin(p)
  z=r*dcos(t)
  pos_cl_cart(i,1)=x
  pos_cl_cart(i,2)=y
  pos_cl_cart(i,3)=z
enddo

end subroutine sph_2_cart
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
! Transform Cartesian coordinates (x) into spherical ones (x too)
! -----------------------------------------------------------------------------------
subroutine cart_to_sph(x)
use miscellaneous
double precision, intent(inout) :: x(1:3)
double precision :: r(1:3) ! where we store the Cartesian coordinates

r=x
x(1)=dsqrt(r(1)**two+r(2)**two+r(3)**two)
x(2)=datan(dsqrt(r(1)**two+r(2)**two)/r(3))
x(3)=datan(r(2)/r(1))
if (r(3)<zero) x(2)=x(2)+dpi
if ((r(1) < zero .and. r(2) > zero) .or. (r(1) < zero .and. r(2) < zero)) then
   x(3)=x(3)+dpi
endif
if (r(1) > zero .and. r(2) < zero) then
   x(3)=x(3)+dpi*two
endif

if ( x(1)**2.d0 - (r(1)**2.d0+r(2)**2.d0+r(3)**2.d0) &
        > smalldble ) then
  call crash('Beware, basic conversion from Cartesian to spherical coordinates criteria not matched')
endif

end subroutine cart_to_sph
! -----------------------------------------------------------------------------------


! ! -----------------------------------------------------------------------------------
! subroutine cart_to_sph_1pt(x,w)
! double precision, intent(inout) :: x(1:3), w(1:nw) ! We only use the 3 velocities ie w(1:3)
! double precision :: r(1:3), z(1:nw) ! We only use the 3 velocities ie z(1:3)
! r=x
! z=w
! x(1)=dsqrt(r(1)**two+r(2)**two+r(3)**two)
! x(2)=datan(dsqrt(r(1)**two+r(2)**two)/r(3))
! x(3)=datan(r(2)/r(1))
! if (r(3)<zero) x(2)=x(2)+dpi
! if ((r(1) < zero .and. r(2) > zero) .or. (r(1) < zero .and. r(2) < zero)) then
!    x(3)=x(3)+dpi
! endif
! if (r(1) > zero .and. r(2) < zero) then
!    x(3)=x(3)+dpi*two
! endif
! w(1)=z(1)*(r(1)/x(1))+z(2)*(r(2)/x(1))+z(3)*(r(3)/x(1))
! w(2)=z(1)*(r(1)/dsqrt(r(1)**two+r(2)**two))*(r(3)/x(1))+&
!         z(2)*(r(2)/dsqrt(r(1)**two+r(2)**two))*(r(3)/x(1))-&
!         z(3)*(dsqrt(r(1)**two+r(2)**two)/x(1))
! w(3)=-z(1)*(r(2)/dsqrt(r(1)**two+r(2)**two))+z(2)*(r(1)/dsqrt(r(1)**two+r(2)**two))
! ! print*, dsqrt(r(1)**two+r(2)**two+r(3)**two), dsqrt(z(1)**two+z(2)**two+z(3)**two)
! if ( dabs(w(1)**two+w(2)**two+w(3)**two-(z(1)**two+z(2)**two+z(3)**two))>1.d-8 ) then
!    write(3,*), 'Error#12 : basic sanity check not passed',&
!                dabs(w(1)**two+w(2)**two+w(3)**two-(z(1)**two+z(2)**two+z(3)**two)),&
!                dabs(w(1)**two+w(2)**two+w(3)**two-(z(1)**two+z(2)**two+z(3)**two))
!    stop
! endif
! end subroutine cart_to_sph_1pt
! ! -----------------------------------------------------------------------------------


end module mod_cart_sph
