module m_construction_Hk
use m_FT_Ham_base
use m_FT_Ham_coo
implicit none

private
public :: get_Hk

contains

subroutine get_Hk(Hk_inp,k)
    !unfolds the real Hamiltonian into a 2 real Hamiltonians (one for the real and for the complex part)
    type(H_inp_real_to_k),intent(in)   :: Hk_inp(:)
    real(8),intent(in)                 :: k(3)
!    type(H_inp_k_coo),intent(inout)      :: H_out_k
!    type(H_inp_k_coo),allocatable      :: Htmp_real,Htmp_compl

    integer     :: i_ham,N_ham,i_shell
    complex(8),allocatable  ::  val(:)
    integer,allocatable     ::  row(:),col(:)

    real(8)  :: phase

    N_ham=size(HK_inp)
!    allocate(H_out_k(N_ham),mold=Hk_inp)

    do i_ham=1,N_ham

!       Call set_H(H_out_k(i_ham),h_io)
!       allocate(H_out_k,mold=H)

        do i_shell=1,size(Hk_inp(i_ham)%H)

            phase=dot_product(Hk_inp(i_ham)%diffR(:,i_shell),k)
            write(*,*) i_shell, Hk_inp(i_ham)%diffR(:,i_shell)
            write(*,*) phase
            pause
!            Call Hk_inp%H(i_ham)%get_hinit(hinit)
!            Call Hk_inp%H(i_ham)%get_par(val,row,col)

            val=val*cmplx(cos(phase),sin(phase),8)

!            Call H_out_k%init_coo(val,row,col,hinit)

!            Call H_out_k%add(Htmp_real)
!            Call H_compl%add(Htmp_compl)

        enddo

!        Call Htmp_real%destroy()
!        Call Htmp_compl%destroy()
        deallocate(val,row,col)

    enddo
end subroutine

end module m_construction_Hk
