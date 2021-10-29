module m_construction_Hk
use m_FT_Ham_base
use m_FT_Ham_dense
use m_FT_Ham_coo_rtok_base
use m_setH_util
use m_FT_Ham_coo
implicit none

private
public get_Hk

contains

subroutine get_Hk(Hk_inp,k,H_out_k)
    !unfolds the real Hamiltonian into a 2 real Hamiltonians (one for the real and for the complex part)
    type(H_inp_real_to_k),intent(in)     :: Hk_inp(:)
    real(8),intent(in)                   :: k(3)
    class(FT_Ham_base),intent(inout)     :: H_out_k
    class(FT_Ham_base),allocatable       :: Htmp

    integer     :: i_ham,N_ham,i_shell
    complex(8)  ::  val
    integer,allocatable     ::  row(:),col(:)

    real(8)  :: phase


    N_ham=size(HK_inp)
    allocate(Htmp,mold=H_out_k)
!
    do i_ham=1,N_ham

        do i_shell=1,size(Hk_inp(i_ham)%H)

            phase=dot_product(Hk_inp(i_ham)%diffR(:,i_shell),k)
!            call Hk_inp(i_ham)%H(i_shell)%eval_single(1)
!            Call Hk_inp%H(i_ham)%get_hinit(hinit)
!            Call Hk_inp%H(i_ham)%get_par(val,row,col)

            val=val*cmplx(cos(phase),sin(phase),8)    ! multiply the value of the energy matrix with the exp(i*phase)

!            Call H_out_k%init(val,row,col)    ! put the value in the Hamiltonian into a temporary Hamiltonian Htmp

!            Call H_out_k%add(Htmp_real) ! add all the contributions of the Hamiltonian Htmp to the final Hamiltonian H_out_k

        enddo

!        Call Htmp%destroy()     ! destroy the temporary Hamiltonian
        deallocate(val,row,col)

    enddo
end subroutine

end module m_construction_Hk
