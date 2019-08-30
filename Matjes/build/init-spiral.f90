module m_init_spiral
use m_derived_types
implicit none

private
public :: init_spiral

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the starting configuration as a spin spiral
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_spiral(io,fname,my_lattice,my_motif,mode_name,start,end)
use m_get_position
use m_vector
use m_io_utils
use m_convert
type (lattice), intent(inout) :: my_lattice
type (cell), intent(in) :: my_motif
integer, intent(in) :: io,start,end
character(len=*), intent(in) :: fname,mode_name
! internal variables
real(kind=8) :: qvec(3),Rq(3),Iq(3),dumy_vec(3),kvec(3,3),r(3,3),ss
integer :: i_z,i_y,i_x,i_m,Nx,Ny,Nz,Nmag,size_mag
real(kind=8), allocatable :: position(:,:,:,:,:)
character(len=30) :: variable_name

kvec=my_lattice%astar
r=my_lattice%areal
Rq=(/0.0,0.0,1.0/)
Iq=(/1.0,0.0,0.0/)

variable_name=convert('qvec_',mode_name)
call get_parameter(io,fname,variable_name,3,qvec)
dumy_vec=qvec(1)*kvec(1,:)+qvec(2)*kvec(2,:)+qvec(3)*kvec(3,:)
qvec=dumy_vec

variable_name=convert('Rq_',mode_name)
call get_parameter(io,fname,variable_name,3,Rq,1.0d0)
dumy_vec=Rq(1)*r(1,:)+Rq(2)*r(2,:)+Rq(3)*r(3,:)
ss=norm(dumy_vec)
Rq=dumy_vec/ss

variable_name=convert('Iq_',mode_name)
call get_parameter(io,fname,variable_name,3,Iq,1.0d0)
dumy_vec=Iq(1)*r(1,:)+Iq(2)*r(2,:)+Iq(3)*r(3,:)
ss=norm(dumy_vec)
Iq=dumy_vec/ss

! get the position of the sites on the lattice
Nx=my_lattice%dim_lat(1)
Ny=my_lattice%dim_lat(2)
Nz=my_lattice%dim_lat(3)
nmag=count(my_motif%atomic(:)%moment.gt.0.0d0)
size_mag=size(my_motif%atomic(:))

allocate(position(3,Nx,Ny,Nz,Nmag))
call get_position(position,my_lattice%dim_lat,my_lattice%areal,my_motif)

do i_m=1,size_mag
   if (my_motif%atomic(i_m)%moment.lt.1.0d-8) cycle
   do i_z=1,Nz
      do i_y=1,Ny
         do i_x=1,Nx
!normal spin spiral
            my_lattice%l_modes(i_x,i_y,i_z,i_m)%w(start:end)=( cos( dot_product(qvec,position(:,i_x,i_y,i_z,i_m)) )*Rq+ &
         sin( dot_product(qvec,position(:,i_x,i_y,i_z,i_m)) )*Iq)
! inomegenous spin spiral in mulitlayer
!        if (l.eq.1) Spin(4:6,i,j,k,l)=(/0.0d0,0.0d0,1.0d0/)
! 3x3 cubic
! define by theta and phi
!        coefftheta=dot_product(qvec,Spin(1:3,i,j,k,l))
!        coeffphi=sqrt(Spin(1,i,j,k,l)**2+Spin(2,i,j,k,l)**2)
!        if (abs(coeffphi).gt.1.0d-5) coeffphi=Spin(1,i,j,k,l)/sqrt(Spin(1,i,j,k,l)**2+Spin(2,i,j,k,l)**2)
!        Spin(4:6,i,j,k,l)=(/coeffphi*sin(coefftheta), &
!         sqrt(1-coeffphi**2)*sin(coefftheta), &
!         cos(coefftheta) /)
         enddo
      enddo
   enddo
enddo

deallocate(position)

end subroutine init_spiral

end module m_init_spiral
