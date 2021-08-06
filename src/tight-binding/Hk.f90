module m_Hk
!slow and general implementation of k-space Hamiltonian
use m_derived_types, only: lattice
use m_H_tb_public
use m_tb_types ,only: parameters_TB_IO_H 
use m_ham_init_type ,only: parameters_ham_init 
use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none

private
public Hk_inp_t, Hk_evec, HK_eval

type Hk_inp_t
    type(H_tb_coo),allocatable  :: H(:)
    real(8),allocatable         :: diffR(:,:)
contains
    procedure :: combine => Hk_inp_combine  !combine 2 Hk_inp_t by copying
    procedure :: destroy => Hk_inp_destroy
end type

contains

subroutine Hk_eval(Hk_inp,k,h_io,eval)
    type(Hk_inp_t),intent(in)               :: Hk_inp
    real(8),intent(in)                      :: k(3)
    type(parameters_TB_IO_H),intent(in)     :: h_io
    real(8),intent(out),allocatable         :: eval(:)
    class(H_tb),allocatable                 :: H

    Call get_Hk(Hk_inp,k,h_io,H)
    Call H%get_eval(eval)
    Call H%destroy
end subroutine

subroutine Hk_evec(Hk_inp,k,h_io,eval,evec)
    type(Hk_inp_t),intent(in)               :: Hk_inp
    real(8),intent(in)                      :: k(3)
    type(parameters_TB_IO_H),intent(in)     :: h_io
    real(8),allocatable,intent(out)         :: eval(:)
    complex(8),allocatable,intent(out)      :: evec(:,:)
    class(H_tb),allocatable                 :: H

    Call get_Hk(Hk_inp,k,h_io,H)
    Call H%get_evec(eval,evec)
    Call H%destroy
end subroutine

subroutine get_Hk(Hk_inp,k,h_io,H)
    !unfolds the Hk_inp data to and Hamiltonian (H) at a k-point (k)
    type(Hk_inp_t),intent(in)   :: Hk_inp
    real(8),intent(in)          :: k(3)
    type(parameters_TB_IO_H),intent(in)     :: h_io
    class(H_tb),allocatable,intent(out)     :: H
    class(H_tb),allocatable                 :: Htmp

    type(parameters_ham_init)   :: hinit
    integer     :: i_ham,N_ham
    complex(8),allocatable  ::  val(:)
    integer,allocatable     ::  row(:),col(:)

    real(8)  :: phase

    Call set_H(H,h_io)
    allocate(Htmp,mold=H)
    N_ham=size(HK_inp%H)
    do i_ham=1,N_ham
        phase=dot_product(HK_inp%diffR(:,i_ham),k)
        Call Hk_inp%H(i_ham)%get_hinit(hinit)
        Call Hk_inp%H(i_ham)%get_par(val,row,col)
        val=val*cmplx(cos(phase),sin(phase),8)
        Call Htmp%init_coo(val,row,col,hinit)
        Call H%add(Htmp)
        Call Htmp%destroy()
        deallocate(val,row,col)
    enddo
end subroutine

subroutine Hk_inp_combine(Hout,H1,H2)
    class(Hk_inp_t),intent(inout)   ::  Hout
    type(Hk_inp_t),intent(in)       ::  H1, H2

    integer     :: size_in(2), size_out, i

    size_in=[size(H1%H),size(H2%H)]
    size_out=sum(size_in)
    if(allocated(Hout%H)) deallocate(Hout%H)
    if(allocated(Hout%diffR)) deallocate(Hout%diffR)
    allocate(Hout%H(size_out),Hout%diffR(3,size_out))
    do i=1,size_in(1)
        Call H1%H(i)%copy(Hout%H(i))
    enddo
    Hout%diffR(:,1:size_in(1))=H1%diffR
    do i=1,size_in(2)
        Call H2%H(i)%copy(Hout%H(i+size_in(1)))
    enddo
    Hout%diffR(:,size_in(1)+1:size_out)=H2%diffR
end subroutine

subroutine Hk_inp_destroy(this)
    class(Hk_inp_t),intent(inout)   ::  this
    integer ::  i

    do i=1,size(this%H)
        Call this%H(i)%destroy()
    enddo
    deallocate(this%H)
    deallocate(this%diffR)
end subroutine

end module


