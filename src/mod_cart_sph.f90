module mod_cart_sph
use miscellaneous
contains

! -----------------------------------------------------------------------------------
subroutine sph_to_cart(x,w)
include 'def.f90'
double precision, intent(inout) :: x(1:3), w(1:nw) ! We only use the 3 velocities ie w(1:3)
double precision :: r(1:3), z(1:nw) ! We only use the 3 velocities ie z(1:3)
r=x
z=w
x(1)=r(1)*dsin(r(2))*dcos(r(3))
x(2)=r(1)*dsin(r(2))*dsin(r(3))
x(3)=r(1)*dcos(r(2))
! For a radial launch (TO CHANGE IF WE ACCOUNT FOR A NORTH/SOUTH THETA VELOCITY Z(2))
w(1)=z(1)*dsin(r(2))*dcos(r(3))-z(3)*dsin(r(3))+z(2)*dcos(r(2))*dcos(r(3))
w(2)=z(1)*dsin(r(2))*dsin(r(3))+z(3)*dcos(r(3))+z(2)*dcos(r(2))*dsin(r(3))
w(3)=z(1)*dcos(r(2))-z(2)*dsin(r(2))
end subroutine sph_to_cart
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
subroutine cart_to_sph(x,w,ind,rotated)
include 'def.f90'
integer, intent(in) :: ind
double precision, intent(inout) :: x(ind,1:3), w(ind,1:nw) ! We only use the 3 velocities ie w(1:3)
double precision :: r(ind,1:3), z(ind,1:nw) ! We only use the 3 velocities ie z(1:3)
logical, intent(in), optional :: rotated
r=x
z=w
x(:,1)=dsqrt(r(:,1)**two+r(:,2)**two+r(:,3)**two)
x(:,2)=datan(dsqrt(r(:,1)**two+r(:,2)**two)/r(:,3))
x(:,3)=datan(r(:,2)/r(:,1))
where (r(:,3)<zero) x(:,2)=x(:,2)+dpi
! if (Lorb .EQ. 'par') then ! for 'prp', we want to remain between - and + pi/2
where ((r(:,1) < zero .and. r(:,2) > zero) .or. (r(:,1) < zero .and. r(:,2) < zero))
   x(:,3)=x(:,3)+dpi
endwhere
where (r(:,1) > zero .and. r(:,2) < zero)
   x(:,3)=x(:,3)+dpi*two
endwhere
! where (x(:,3)<-dpi/two)
!   x(:,3)=x(:,3)+two*dpi
! endwhere
! endif
! if (present(rotated) .and. rotated) then
!
! else
w(:,1)=z(:,1)*(r(:,1)/x(:,1))+z(:,2)*(r(:,2)/x(:,1))+z(:,3)*(r(:,3)/x(:,1))
w(:,2)=z(:,1)*(r(:,1)/dsqrt(r(:,1)**two+r(:,2)**two))*(r(:,3)/x(:,1))+&
       z(:,2)*(r(:,2)/dsqrt(r(:,1)**two+r(:,2)**two))*(r(:,3)/x(:,1))-&
       z(:,3)*(dsqrt(r(:,1)**two+r(:,2)**two)/x(:,1))
w(:,3)=-z(:,1)*(r(:,2)/dsqrt(r(:,1)**two+r(:,2)**two))+z(:,2)*(r(:,1)/dsqrt(r(:,1)**two+r(:,2)**two))
! endif
if ( any(dabs(w(:,1)**two+w(:,2)**two+w(:,3)**two-(z(:,1)**two+z(:,2)**two+z(:,3)**two))>1.d-8 )) then
   write(3,*), 'Error#12 : basic sanity check not passed',&
           maxval(dabs(w(:,1)**two+w(:,2)**two+w(:,3)**two-(z(:,1)**two+z(:,2)**two+z(:,3)**two))),&
           maxloc(dabs(w(:,1)**two+w(:,2)**two+w(:,3)**two-(z(:,1)**two+z(:,2)**two+z(:,3)**two)))
   stop
endif

end subroutine cart_to_sph
! -----------------------------------------------------------------------------------


! -----------------------------------------------------------------------------------
subroutine cart_to_sph_1pt(x,w)
include 'def.f90'
double precision, intent(inout) :: x(1:3), w(1:nw) ! We only use the 3 velocities ie w(1:3)
double precision :: r(1:3), z(1:nw) ! We only use the 3 velocities ie z(1:3)
r=x
z=w
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
w(1)=z(1)*(r(1)/x(1))+z(2)*(r(2)/x(1))+z(3)*(r(3)/x(1))
w(2)=z(1)*(r(1)/dsqrt(r(1)**two+r(2)**two))*(r(3)/x(1))+&
        z(2)*(r(2)/dsqrt(r(1)**two+r(2)**two))*(r(3)/x(1))-&
        z(3)*(dsqrt(r(1)**two+r(2)**two)/x(1))
w(3)=-z(1)*(r(2)/dsqrt(r(1)**two+r(2)**two))+z(2)*(r(1)/dsqrt(r(1)**two+r(2)**two))
! print*, dsqrt(r(1)**two+r(2)**two+r(3)**two), dsqrt(z(1)**two+z(2)**two+z(3)**two)
if ( dabs(w(1)**two+w(2)**two+w(3)**two-(z(1)**two+z(2)**two+z(3)**two))>1.d-8 ) then
   write(3,*), 'Error#12 : basic sanity check not passed',&
               dabs(w(1)**two+w(2)**two+w(3)**two-(z(1)**two+z(2)**two+z(3)**two)),&
               dabs(w(1)**two+w(2)**two+w(3)**two-(z(1)**two+z(2)**two+z(3)**two))
   stop
endif
end subroutine cart_to_sph_1pt
! -----------------------------------------------------------------------------------


end module mod_cart_sph
