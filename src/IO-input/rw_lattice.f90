subroutine rw_lattice(my_lattice)
use m_vector, only : norm,cross
use m_constants, only : pi
use m_derived_types
use m_io_files_utils
use m_io_utils
use, intrinsic  ::  ISO_FORTRAN_ENV, only: output_unit
implicit none
! out
type(lattice), intent(out) :: my_lattice
! dummy variable
integer :: j,dimension(3),io
integer :: N_site
logical :: Periodic_log(3),fexist
real (kind=8) :: a(3),ss,net(3,3)
real(8)         :: rot_angles(3)
real(8)         :: rot_mat(3,3)
integer         :: i

write(output_unit,'(a)') 'reading the lattice parameters in the input file'
inquire(file='input',exist=fexist)
if(.not. fexist)then
    write(*,*) "Did not find input file, you should at least supply some lattices"
    STOP "Supply input file"
endif
io=open_file_read('input')

call get_parameter(io,'input','Periodic_log',3,Periodic_log)
call get_parameter(io,'input','alat',3,a)

call get_parameter(io,'input','lattice',3,3,net)
do j=1,3
   ss=norm(net(j,:))
   if (ss.lt.1.0d-8) stop 'lattice vectors are crap'
   net(j,:)=net(j,:)/ss
   a(j)=a(j)*ss
enddo

rot_angles=0.0d0
call get_parameter(io,'input','rotate_lattice',3,rot_angles)
if(any(rot_angles/=0.d00))then
    write(output_unit,'(A)') 'Found "lattice_rotangle" input which does a Euler angle rotation around ZXZ:'
    write(output_unit,'(3X,A,F16.8,A)') '1. rotation (z):', rot_angles(1),' pi'
    write(output_unit,'(3X,A,F16.8,A)') '2. rotation (x):', rot_angles(2),' pi'
    write(output_unit,'(3X,A,F16.8,A)') '3. rotation (z):', rot_angles(3),' pi'
    rot_angles=rot_angles*pi
    associate( s1 => sin(rot_angles(1)), s2 => sin(rot_angles(2)), s3 => sin(rot_angles(3)),&
             & c1 => cos(rot_angles(1)), c2 => cos(rot_angles(2)), c3 => cos(rot_angles(3)))
        rot_mat(1,1)=  c1 * c3      - c2 * s1 * s3 
        rot_mat(2,1)=  c3 * s1      + c1 * c2 * s3 
        rot_mat(3,1)=  s2 * s3
        rot_mat(1,2)= -c1 * s3      - c2 * c3 * s1 
        rot_mat(2,2)=  c1 * c2 * c3 - s1 * s3
        rot_mat(3,2)=  c3 * s2
        rot_mat(1,3)=  s1 * s2
        rot_mat(2,3)= -c1 * s2
        rot_mat(3,3)=  c2
    end associate
    write(*,*) 'rot_mat=',rot_mat(:,:)
    write(output_unit,'(A,3(/3X,3F16.8)/)') 'Initial lattice vectors:',net
    net=matmul(net,rot_mat)
    write(output_unit,'(A,3(/3X,3F16.8)/)') 'Rotated lattice vectors:',net
endif


call get_parameter(io,'input','Nsize',3,dimension)
if (dimension(1)*dimension(2)*dimension(3).eq.0) then
   write(output_unit,*) 'dimension of lattice must be >0 along all directions'
   stop
elseif((dimension(1).ne.dimension(2)).and.(dimension(2).ne.1)) then
   write(output_unit,*) 'asymetric cell choosen'
elseif((dimension(1).ne.dimension(2)).and.(dimension(2).eq.1)) then
   write(output_unit,*) 'chain geometry choosen'
endif
N_site=dimension(1)*dimension(2)*dimension(3)

call close_file('input',io)


!! put the magnetic lattice in order
Call my_lattice%init_geo(net,a,dimension,Periodic_log)

end subroutine rw_lattice
