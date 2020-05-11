module m_grp_sym
use m_basic_types, only : symop
use m_derived_types, only : cell
use m_vector, only : norm
use m_constants, only : pi
use m_user_info
use m_io_files_utils

! private variable
!type(symop),allocatable :: sym_operations(:)

private
public :: get_group,read_symmetries,get_num_sym_file
contains
!
! find all the symmetry operations that apply to the lattice and the uni cell
! it requires several subroutine and functions
! subroutine: get_latt_sym   gets all the symmetry operations that leave the lattice invariant
!             get_pt_sym     gets all the symmetry operations that leave the lattice and the motif invariant
!             write_symmetries     write all the symmetry operations into 'symmetries.out'
! function:   get_symop      get all possible symmetry operations
!             look_translation    find the equivalent position in the unit cell
!
subroutine get_group(areal,my_motif,periodic,dim_lat)
implicit none
real(kind=8), intent(in) :: areal(3,3)
type (cell), intent(in) :: my_motif
logical, intent(in) :: periodic(:)
integer, intent(in) :: dim_lat(:)
! internal variables
integer :: number_sym,sym_index(64),io_sym
real(kind=8) :: time
! first step determine the lattice symmetries

time=0.0d0
call user_info(6,time,'calculating the symmetry operations',.false.)

number_sym=0
sym_index=0

call get_latt_sym(areal,number_sym,sym_index)

write(6,'(/,a,I2,a,/)') 'The lattice has  ',number_sym,'  symmetrie operations'

call get_pt_sym(areal,number_sym,sym_index,my_motif,periodic,dim_lat)

io_sym=open_file_write('symmetries.out')

call write_symmetries(io_sym,number_sym,sym_index(1:number_sym))

call close_file('symmetries.out',io_sym)

call user_info(6,time,'done',.true.)

call write_symmetries(6,number_sym,sym_index)

end subroutine get_group

! find the symmetry of the unit cell
! Take the lattice vectors R(3,3) with basis vector ordered in lines
! apply a symmetry operation P on it P.R with vectors in column
! Calculate P.R-R
! find the linear combination so that (P.R-R) must be 0
! Continue until there are no symmetry operation left
!
subroutine get_latt_sym(areal,number_sym,sym_index)
implicit none
real(kind=8), intent(in) :: areal(3,3)
integer, intent(inout) :: number_sym,sym_index(:)
! internal
type(symop) :: all_sym_op(64)
integer :: i,j
real(kind=8) :: rtest(3,3),areal_rot(3,3)
logical :: found(3)

call get_symop(all_sym_op)

number_sym=0
! loop over all symmetry
do i=1,64

!
! I multiply the rotation matrice written in line with the lattice vector written in lign too
! I have to take the transpose of alat to have it written in column
! the matrix rtest is then written in column in (column,line)
!
      rtest=matmul(all_sym_op(i)%mat,transpose(areal))


!      if ((.not.periodic(3)).and.(rtest(3,3).ne.areal(3,3))) cycle

!
! I have to take the transpose of rtest to have it written in lign again
!
      areal_rot=transpose(rtest)

      do j=1,3
         found(j)=look_translation(areal_rot(j,:),areal,(/.true.,.true.,.true./),(/2,2,2/))
      enddo

      if (all(found)) then
         number_sym=number_sym+1
         sym_index(number_sym)=i
      endif

enddo

write(6,'(a)') ''
write(6,'(a,I4,a)') 'number of lattice symmetries found  ', number_sym, '/64'
do j=1,number_sym
   write(6,'(a)') all_sym_op(sym_index(j))%name
   do i=1,3
      write(6,'(3f8.3)') all_sym_op(sym_index(j))%mat(i,:)
   enddo
enddo
write(6,'(a)') ''

end subroutine get_latt_sym





! find the symmetry of each atomic position in the unit cell compatible with the lattice symmetry
! Take the lattice vectors R(3,3) with basis vector ordered in lines
! Take to atom position rho(3)
! Calculate R.rho
! apply a symmetry operation P on it P.(R.rho)
! Calculate P.(R.rho)-R
! find the linear combination so that (P.(R.rho)-R) must be 0
! Continue until there are no symmetry operation left
!
subroutine get_pt_sym(areal,number_sym,sym_index,my_motif,periodic,dim_lat)
implicit none
real(kind=8), intent(in) :: areal(3,3)
integer, intent(in) :: dim_lat(:)
logical, intent(in) :: periodic(:)
integer, intent(inout) :: number_sym,sym_index(:)
type(cell), intent(in) :: my_motif
!internal
integer :: natom,i,j,i_sim,new_index(64),k
type(symop) :: all_sym_op(64)
real(kind=8) :: test_vec(3),pos(3)
logical :: found

natom=size(my_motif%atomic)
new_index=0

call get_symop(all_sym_op)

write(6,'(I4,a)') natom, ' atoms found in the unit cell'

do i=1,natom
   if (norm(my_motif%atomic(i)%position).lt.1.0d-8) then
      write(6,'(a,I4,a)') 'atom ',i,' has the position (0,0,0)'
      write(6,'(a)') 'no non-symorphic translation'
      cycle
   endif

   do j=1,number_sym
      i_sim=sym_index(j)

! position in cartesian coordinate
! the basis vectors have the fomat (line, column)
! take the transpose before making the product

      pos=matmul(transpose(areal),my_motif%atomic(i)%position)
      test_vec=matmul(all_sym_op(i_sim)%mat,pos)

      found=look_translation(test_vec,areal,periodic,dim_lat,pos)

      if (.not.found) then
         sym_index(j)=0
      endif

   enddo

!
! reorder the symmetries
!
   j=0
   do k=1,number_sym
     if (sym_index(k).eq.0) cycle
     j=j+1
     new_index(j)=sym_index(k)
   enddo
   number_sym=j
   sym_index=new_index

enddo

j=0
do i=1,number_sym

   if (sym_index(i).eq.0) cycle
   j=j+1
   new_index(j)=sym_index(i)

enddo

number_sym=j
sym_index=new_index
write(6,'(a,I3,a)') 'space group has ',number_sym,'  symmetry operations'

end subroutine get_pt_sym

!
! find the R(3) lattice vectors so that P.R-R=0
!
function look_translation(areal_rot,areal,periodic,dim_lat,translation)
implicit none
real(kind=8), intent(in) :: areal_rot(3),areal(3,3)
real(kind=8), intent(in), optional :: translation(3)
integer, intent(in) :: dim_lat(:)
logical :: look_translation,periodic(:)
!internal
integer :: u,v,w
real(kind=8) :: test_vec(3),eps(3)

if (present(translation)) then
   eps=translation
else
   eps=0.0d0
endif

look_translation=.false.
do u=-2,2,1
   do v=-2,2,1
      do w=-2,2,1

         test_vec=areal_rot-eps
         if ((periodic(1)).and.(dim_lat(1).ne.1)) test_vec=test_vec+real(u)*areal(1,:)

         if ((periodic(2)).and.(dim_lat(2).ne.1)) test_vec=test_vec+real(v)*areal(2,:)

         if ((periodic(3)).and.(dim_lat(3).ne.1)) test_vec=test_vec+real(w)*areal(3,:)

!
! be very carefull here! If the positions are not given with enough precision, the symetries will not be found
!

         if (norm(test_vec).lt.1.0d-6) then
            look_translation=.true.
            return
         endif

      enddo
   enddo
enddo

end function look_translation




!
! subroutine that write down the symmetry operations that were found
!
subroutine write_symmetries(io_sym,number_sym,sym_index)
implicit none
integer, intent(in) :: io_sym,number_sym,sym_index(:)
! internal variables
integer :: i,i_sym,j
type(symop) :: all_sym_op(64)

call get_symop(all_sym_op)

write(io_sym,'(I4)') number_sym
do i=1,number_sym
   i_sym=sym_index(i)
   write(io_sym,'(a)') '!'
   write(io_sym,'(a)') all_sym_op(i_sym)%name
   do j=1,3
      write(io_sym,'(3(f14.10,3x))') all_sym_op(i_sym)%mat(j,:)
   enddo
enddo

end subroutine write_symmetries






!
! read the number of symetry operations
!
integer function get_num_sym_file()
implicit none
! internal
integer :: io_sym,n_sym

get_num_sym_file=0
io_sym=open_file_read('symmetries.out')

read(io_sym,*) n_sym

get_num_sym_file=n_sym

if (get_num_sym_file.eq.0) then
   write(6,'(a)') 'problem reading symmetries.out'
   stop
endif

call close_file('symmetries.out',io_sym)

end function get_num_sym_file

subroutine read_symmetries(n_sym,symmetries)
implicit none
integer, intent(in) :: n_sym
type(symop), intent(out) :: symmetries(n_sym)
! internal
integer :: io_sym,i,j


io_sym=open_file_read('symmetries.out')

read(io_sym,*)
do i=1,n_sym
   read(io_sym,*)
   read(io_sym,*) symmetries(i)%name
   do j=1,3
      read(io_sym,*) symmetries(i)%mat(j,:)
   enddo
enddo

call close_file('symmetries.out',io_sym)

end subroutine







subroutine get_symop(rotmat)
! **********************************************
! This subroutine defines the rotation matrices for
! all the 32 point groups and names them after
! J.F. Cornwell (Group Theory??) second edition
! Appendix D, p 324-325
!
!  rotmat , real*8(3,3,64) : All 64 3x3 rotation matriced of the point groups
!  rotname , char*10(64)   : Names of the 64 rotation matrixes after J.F. Cornwell
!
! the basis vectors have the fomat (line,column)
! *********************************************
implicit none
integer i1,is
type(symop),intent(out) :: ROTMAT(64)
real(kind=8) :: RTHREE,HALF

RTHREE = SQRT(3.d0)/2.d0
HALF = 0.5d0
! set to zero
do i1=1,64
     ROTMAT(i1)%mat = 0.d0
end do
!
      ROTMAT(1)%mat(1,1) =  1.d0
      ROTMAT(1)%mat(2,2) =  1.d0
      ROTMAT(1)%mat(3,3) =  1.d0
      ROTMAT(1)%name = 'E'
!
      ROTMAT(2)%mat(1,2) =  1.d0
      ROTMAT(2)%mat(2,3) = -1.d0
      ROTMAT(2)%mat(3,1) = -1.d0
      ROTMAT(2)%name = 'C3alfa'
!
      ROTMAT(3)%mat(1,2) = -1.d0
      ROTMAT(3)%mat(2,3) = -1.d0
      ROTMAT(3)%mat(3,1) =  1.d0
      ROTMAT(3)%name = 'C3beta '
!
      ROTMAT(4)%mat(1,2) = -1.d0
      ROTMAT(4)%mat(2,3) =  1.d0
      ROTMAT(4)%mat(3,1) = -1.d0
      ROTMAT(4)%name = 'C3gamma'
!
      ROTMAT(5)%mat(1,2) = 1.d0
      ROTMAT(5)%mat(2,3) = 1.d0
      ROTMAT(5)%mat(3,1) = 1.d0
      ROTMAT(5)%name = 'C3delta '
!
      ROTMAT(6)%mat(1,3) = -1.d0
      ROTMAT(6)%mat(2,1) =  1.d0
      ROTMAT(6)%mat(3,2) = -1.d0
      ROTMAT(6)%name = 'C3alfa-1'
!
      ROTMAT(7)%mat(1,3) =  1.d0
      ROTMAT(7)%mat(2,1) = -1.d0
      ROTMAT(7)%mat(3,2) = -1.d0
      ROTMAT(7)%name = 'C3beta-1 '
!
      ROTMAT(8)%mat(1,3) = -1.d0
      ROTMAT(8)%mat(2,1) = -1.d0
      ROTMAT(8)%mat(3,2) =  1.d0
      ROTMAT(8)%name = 'C3gamma-1'
!
      ROTMAT(9)%mat(1,3) =  1.d0
      ROTMAT(9)%mat(2,1) =  1.d0
      ROTMAT(9)%mat(3,2) =  1.d0
      ROTMAT(9)%name = 'C3delta-1'
!
      ROTMAT(10)%mat(1,1) =  1.d0
      ROTMAT(10)%mat(2,2) = -1.d0
      ROTMAT(10)%mat(3,3) = -1.d0
      ROTMAT(10)%name = 'C2x'
!
      ROTMAT(11)%mat(1,1) = -1.d0
      ROTMAT(11)%mat(2,2) =  1.d0
      ROTMAT(11)%mat(3,3) = -1.d0
      ROTMAT(11)%name = 'C2y'
!
      ROTMAT(12)%mat(1,1) = -1.d0
      ROTMAT(12)%mat(2,2) = -1.d0
      ROTMAT(12)%mat(3,3) =  1.d0
      ROTMAT(12)%name = 'C2z'
!
      ROTMAT(13)%mat(1,1) =  1.d0
      ROTMAT(13)%mat(2,3) =  1.d0
      ROTMAT(13)%mat(3,2) = -1.d0
      ROTMAT(13)%name = 'C4x'
!
      ROTMAT(14)%mat(1,3) = -1.d0
      ROTMAT(14)%mat(2,2) =  1.d0
      ROTMAT(14)%mat(3,1) =  1.d0
      ROTMAT(14)%name = 'C4y '
!
      ROTMAT(15)%mat(1,2) =  1.d0
      ROTMAT(15)%mat(2,1) = -1.d0
      ROTMAT(15)%mat(3,3) =  1.d0
      ROTMAT(15)%name = 'C4z'
!
      ROTMAT(16)%mat(1,1) =  1.d0
      ROTMAT(16)%mat(2,3) = -1.d0
      ROTMAT(16)%mat(3,2) =  1.d0
      ROTMAT(16)%name = 'C4x-1 '
!
      ROTMAT(17)%mat(1,3) =  1.d0
      ROTMAT(17)%mat(2,2) =  1.d0
      ROTMAT(17)%mat(3,1) = -1.d0
      ROTMAT(17)%name = 'C4y-1'
!
      ROTMAT(18)%mat(1,2) = -1.d0
      ROTMAT(18)%mat(2,1) =  1.d0
      ROTMAT(18)%mat(3,3) =  1.d0
      ROTMAT(18)%name = 'C4z-1'
!
      ROTMAT(19)%mat(1,2) =  1.d0
      ROTMAT(19)%mat(2,1) =  1.d0
      ROTMAT(19)%mat(3,3) = -1.d0
      ROTMAT(19)%name = 'C2a'
!
      ROTMAT(20)%mat(1,2) = -1.d0
      ROTMAT(20)%mat(2,1) = -1.d0
      ROTMAT(20)%mat(3,3) = -1.d0
      ROTMAT(20)%name = 'C2b'
!
      ROTMAT(21)%mat(1,3) =  1.d0
      ROTMAT(21)%mat(2,2) = -1.d0
      ROTMAT(21)%mat(3,1) =  1.d0
      ROTMAT(21)%name = 'C2c'
!
      ROTMAT(22)%mat(1,3) = -1.d0
      ROTMAT(22)%mat(2,2) = -1.d0
      ROTMAT(22)%mat(3,1) = -1.d0
      ROTMAT(22)%name = 'C2d'
!
      ROTMAT(23)%mat(1,1) = -1.d0
      ROTMAT(23)%mat(2,3) =  1.d0
      ROTMAT(23)%mat(3,2) =  1.d0
      ROTMAT(23)%name = 'C2e'
!
      ROTMAT(24)%mat(1,1) = -1.d0
      ROTMAT(24)%mat(2,3) = -1.d0
      ROTMAT(24)%mat(3,2) = -1.d0
      ROTMAT(24)%name = 'C2f'

do i1=1,24
   ROTMAT(i1+24)%mat = -ROTMAT(i1)%mat
   ROTMAT(i1+24)%name = 'I'//ROTMAT(i1)%name(1:9)
end do
!
!
!*********************************************
! Trigonal and hexagonal groups
!*********************************************
!
      ROTMAT(49)%mat(1,1) = -HALF
      ROTMAT(49)%mat(1,2) =  RTHREE
      ROTMAT(49)%mat(2,1) = -RTHREE
      ROTMAT(49)%mat(2,2) = -HALF
      ROTMAT(49)%mat(3,3) =  1.d0
      ROTMAT(49)%name = 'C3z'
!
      ROTMAT(50)%mat(1,1) = -HALF
      ROTMAT(50)%mat(1,2) = -RTHREE
      ROTMAT(50)%mat(2,1) =  RTHREE
      ROTMAT(50)%mat(2,2) = -HALF
      ROTMAT(50)%mat(3,3) =  1.d0
      ROTMAT(50)%name = 'C3z-1'
!
      ROTMAT(51)%mat(1,1) =  HALF
      ROTMAT(51)%mat(1,2) =  RTHREE
      ROTMAT(51)%mat(2,1) = -RTHREE
      ROTMAT(51)%mat(2,2) =  HALF
      ROTMAT(51)%mat(3,3) =  1.d0
      ROTMAT(51)%name = 'C6z'
!
      ROTMAT(52)%mat(1,1) =  HALF
      ROTMAT(52)%mat(1,2) = -RTHREE
      ROTMAT(52)%mat(2,1) =  RTHREE
      ROTMAT(52)%mat(2,2) =  HALF
      ROTMAT(52)%mat(3,3) =  1.d0
      ROTMAT(52)%name = 'C6z-1'
!
      ROTMAT(53)%mat(1,1) = -HALF
      ROTMAT(53)%mat(1,2) =  RTHREE
      ROTMAT(53)%mat(2,1) =  RTHREE
      ROTMAT(53)%mat(2,2) =  HALF
      ROTMAT(53)%mat(3,3) = -1.d0
      ROTMAT(53)%name = 'C2A'
!
      ROTMAT(54)%mat(1,1) = -HALF
      ROTMAT(54)%mat(1,2) = -RTHREE
      ROTMAT(54)%mat(2,1) = -RTHREE
      ROTMAT(54)%mat(2,2) =  HALF
      ROTMAT(54)%mat(3,3) = -1.d0
      ROTMAT(54)%name = 'C2B'
!
      ROTMAT(55)%mat(1,1) =  HALF
      ROTMAT(55)%mat(1,2) = -RTHREE
      ROTMAT(55)%mat(2,1) = -RTHREE
      ROTMAT(55)%mat(2,2) = -HALF
      ROTMAT(55)%mat(3,3) = -1.d0
      ROTMAT(55)%name = 'C2C'
!
      ROTMAT(56)%mat(1,1) =  HALF
      ROTMAT(56)%mat(1,2) =  RTHREE
      ROTMAT(56)%mat(2,1) =  RTHREE
      ROTMAT(56)%mat(2,2) = -HALF
      ROTMAT(56)%mat(3,3) = -1.d0
      ROTMAT(56)%name = 'C2D'

do is=1,8
   ROTMAT(56+is)%mat = -ROTMAT(48+is)%mat
   ROTMAT(56+is)%name = 'I'//ROTMAT(48+is)%name(1:9)
end do
      
!-----------------------------------------
end subroutine

end module m_grp_sym
