module m_sym_utils
use m_rotation_matrix, only : rotate_matrix,check_rotate_matrix,rotation_matrix_real
implicit none

private
public :: look_translation,rotate_force,rotate_exchange

interface look_translation
   module procedure look_translation_vector,look_translation_lattice
end interface look_translation

contains

subroutine rotate_exchange(mat_out,mat_in,rotmat)
   real(8), intent(out)   :: mat_out(:,:)
   real(8), intent(in)    :: mat_in(:,:)
   real(8), intent(in)    :: rotmat(:,:)

   real(8)                :: mat(3,3),sym_part(3,3),antisym_part(3,3)

! decompose in symmetric and antisymmetric part

   sym_part=(mat_in+transpose(mat_in))/2.0d0
   antisym_part=(mat_in-transpose(mat_in))/2.0d0

! rotate only the antisymmetric part

   call rotate_matrix(mat,antisym_part,rotmat)

   mat_out=sym_part+mat

end subroutine

subroutine rotate_force(mat_out,mat_in,rotmat)
   real(8), intent(out)   :: mat_out(:,:)
   real(8), intent(in)    :: mat_in(:,:)
   real(8), intent(in)    :: rotmat(:,:)

! rotate only the full tensor

   call rotate_matrix(mat_out,mat_in,rotmat)

end subroutine

!
! find the R(3) lattice vectors so that P.R-R=0
!
function look_translation_vector(areal_rot,areal,periodic,dim_lat,translation)
use m_vector, only : norm
implicit none
real(kind=8), intent(in) :: areal_rot(3),areal(3,3)
real(kind=8), intent(in), optional :: translation(3)
integer, intent(in) :: dim_lat(:)
logical :: look_translation_vector,periodic(:)
!internal
integer :: u,v,w
real(kind=8) :: test_vec(3),eps(3)

if (present(translation)) then
   eps=translation
else
   eps=0.0d0
endif

look_translation_vector=.false.
do u=-2,2,1
   do v=-2,2,1
      do w=-2,2,1

         test_vec=areal_rot-eps
         if ((periodic(1)).or.(dim_lat(1).ne.1)) then
            test_vec=test_vec+real(u)*areal(1,:)
         else
            if (abs(dot_product(test_vec,areal(1,:))).gt.1.0d-6) test_vec=test_vec-areal(1,:)
         endif

         if ((periodic(2)).or.(dim_lat(2).ne.1)) then
            test_vec=test_vec+real(v)*areal(2,:)
         else
            if (abs(dot_product(test_vec,areal(2,:))).gt.1.0d-6) test_vec=test_vec-areal(2,:)
         endif

         if ((periodic(3)).or.(dim_lat(3).ne.1)) then
            test_vec=test_vec+real(w)*areal(3,:)
         else
            if (abs(dot_product(test_vec,areal(3,:))).gt.1.0d-6) test_vec=test_vec-areal(3,:)
         endif

!
! be very carefull here! If the positions are not given with enough precision, the symetries will not be found
!

         if (norm(test_vec).lt.1.0d-6) then
            look_translation_vector=.true.
            return
         endif

      enddo
   enddo
enddo
end function look_translation_vector




function look_translation_lattice(areal_rot,areal,periodic,dim_lat,translation)
use m_vector, only : norm
implicit none
real(kind=8), intent(in) :: areal_rot(3,3),areal(3,3)
real(kind=8), intent(in), optional :: translation(3)
integer, intent(in) :: dim_lat(:)
logical, intent(in) :: periodic(:)
logical :: look_translation_lattice
!internal
integer :: u,v,w,i,j
real(kind=8) :: test_vec(3),eps(3)
logical :: found(3)

if (present(translation)) then
   eps=translation
else
   eps=0.0d0
endif

look_translation_lattice=.false.
found=.false.

! when no periodicity is present, the axis along the non-periodic direction must be invariant
do i=1,3
  if (.not.(periodic(i)).and.(norm(areal_rot(:,i)-areal(:,i)).gt.1.0d-6)) then
    look_translation_lattice=.false.
    return
  endif
enddo

! this will be checked only if the axis along the open boundary are invariant by the symetry operations

do u=-2,2,1
   do v=-2,2,1
      do w=-2,2,1

         test_vec=real((/u,v,w/),8)
         test_vec=matmul(areal,test_vec)-eps

         do i=1,3
            if (norm(areal_rot(:,i)-test_vec).lt.1.0d-6) found(i)=.true.
         enddo

      enddo
   enddo
enddo

look_translation_lattice=all(found)

end function look_translation_lattice


end module
