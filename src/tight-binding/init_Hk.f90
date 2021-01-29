module m_init_Hk
use m_derived_types, only: lattice,k_grid_t
use m_H_tb_public
use m_tb_types ,only: parameters_TB_Hsolve,parameters_TB_IO_H, parameters_ham_init 
use m_types_tb_h_inp 
use m_neighbor_type, only: neighbors
use m_init_H
use m_delta_onsite
use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none

private
public get_Hk_inp, Hk_inp_t, Hk_evec, HK_eval

type Hk_inp_t
    type(H_tb_coo),allocatable  :: H(:)
    real(8),allocatable         :: diffR(:,:)
contains
    procedure :: combine => Hk_inp_combine  !combine 2 Hk_inp_t by copying
    procedure :: destroy => Hk_inp_destroy
end type

contains

subroutine get_Hk_inp(lat,h_io,Hk_inp)
    use m_init_HR, only: get_delta
    !get Hk_inp, which contains the Hamiltonian data to construct the Hamiltonian at arbitrary k
    type(lattice),intent(in)                :: lat
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(Hk_inp_t),intent(out)              :: Hk_inp

    type(parameters_TB_IO_H)                :: h_io_scf

    if(h_io%use_scf)then    
        H_io_scf=h_io
        if(.true.)then 
            !use reciprocal-space self-consistent delta
            Call get_Hk_inp_conv(lat,h_io,Hk_inp,h_io_scf%del)
        else
            !use real-space self-consistent delta
            Call get_delta(lat,h_io,h_io_scf%del)
            Call get_H(lat,h_io_scf,Hk_inp%H,diffR=Hk_inp%diffR)
        endif
    else
        Call get_H(lat,h_io,Hk_inp%H,diffR=Hk_inp%diffR)
    endif
end subroutine

subroutine get_Hk_inp_conv(lat,h_io,H_inp,del)
    type(lattice),intent(in)                :: lat
    type(parameters_TB_IO_H),intent(in)     :: h_io 
    type(Hk_inp_t),intent(out)              :: H_inp
    type(Hdelta),intent(inout)              :: del

    type(Hk_inp_t)                          :: H_inp_nc
    type(Hk_inp_t)                          :: H_inp_sc

    type(k_grid_t)                          :: k_grid

    character(len=3)    :: i_str
    complex(8)          :: delta_sum(2)
    real(8)             :: diff
    integer             :: i

    !get constant non-superconducting terms
    Call get_H_TB(lat,h_io,H_inp_nc%H,.true.,diffR=H_inp_nc%diffR)
    !get initial superconducting term
    Call get_H_sc(lat,h_io,H_inp_sc%H,del,diffR=H_inp_sc%diffR)
    !set self-consistent k-grid
    Call k_grid%set(lat%a_sc_inv,h_io%scf_kgrid)

    !self-consistent loop
    delta_sum=(0.0d0,0.0d0)
    diff=h_io%scf_diffconv*2.0d0
    do i=1,h_io%scf_loopmax
        if(diff<h_io%scf_diffconv)then
            write(output_unit,'(A,E16.8)') "Reached set convergence criterium for delta-convergence:",h_io%scf_diffconv
            exit
        endif
        Call H_inp%combine(H_inp_nc,H_inp_sc)
        Call H_inp_sc%destroy()
        Call get_delta_scf(H_inp,k_grid,h_io,del)
        Call get_H_sc(lat,h_io,H_inp_sc%H,del,diffR=H_inp_sc%diffR)

        delta_sum(2)=sum(del%delta)/size(del%delta)
        diff=abs(delta_sum(2)-delta_sum(1))
        delta_sum(1)=delta_sum(2)
        write(output_unit,'(2(A,I6)3XA,2E16.8,3XA,E16.8)') 'Finished delta scf loop',i,' of',h_io%scf_loopmax,'  delta-av sum :',delta_sum(1), '  delta-av diff:',diff
        if(h_io%scf_print)then
            write(i_str,'(I0.3)') i
            Call del%print_file("delta_onsite_"//i_str//".dat",lat)
        endif
    enddo
    if(diff>=h_io%scf_diffconv)then
        write(error_unit,'(2(AE14.8))') "WARNING, self-consistent delta in k-space did not reach convergence criteria with ",diff,' instead of ',h_io%scf_diffconv
    endif
    Call H_inp%combine(H_inp_nc,H_inp_sc)
    Call H_inp_nc%destroy()
    Call H_inp_sc%destroy()

    Call del%print_file("delta_onsite_final.dat",lat)
    write(output_unit,'(A,2E16.8)') "Final real delta extremal values:", minval(del%delta%re),maxval(del%delta%re)
end subroutine 

subroutine get_delta_scf(H_inp,k_grid,h_io,del)
    class(Hk_inp_T)             :: H_inp
    type(k_grid_t),intent(in)   :: k_grid
    type(parameters_TB_IO_H),intent(in)     :: h_io 
    type(Hdelta),intent(inout)  :: del

    integer         ::  ik,Nk
    real(8)         ::  k(3)

    del%delta=(0.0d0,0.0d0)
    Nk=k_grid%get_Nk()
    do ik=1,Nk  !MPI parallelization over this loop?
        k=k_grid%get_k(ik)
        Call add_delta_k_scf(H_inp,k,h_io,del)
    enddo
    del%delta=del%delta/real(Nk,8)
end subroutine


subroutine add_delta_k_scf(H_inp,k,h_io,del)
! add delta contribution for a single k as called from get_delta_scf
! use equation 3.28 (Jian-Xin Zhu book, DOI: 10.1007/978-3-319-31314-6 ), using the symmetries only using positive energies (below eq.2.26) analogous to Eq.2.43 as used in real space 
    class(Hk_inp_T)             :: H_inp
    real(8),intent(in)          :: k(3)
    type(parameters_TB_IO_H),intent(in)     :: h_io 
    type(Hdelta),intent(inout)  :: del

    integer                     :: ndim,nspin, nBdG, ncell
    integer                     :: orb
    integer                     :: i_cell
    real(8),allocatable         :: eval(:)
    complex(8),allocatable      :: evec(:,:)
    integer                     :: j,iE
    real(8)                     :: temp !temperature
    logical                     :: vanish
#ifdef CPP_BLAS
    complex(8)          :: arr(size(del%delta,1),2)
#endif

    nspin=h_io%nspin
    nBdG =h_io%norb*h_io%nspin*h_io%ncell
    ndim =h_io%norb*h_io%nspin
    ncell= size(del%delta,1)

    Call Hk_evec(H_inp,k,h_io,eval,evec)
    Call restict_solution_positive(eval,evec,h_io%scf_Ecut,vanish)
    if(vanish) return !no energy value in considered energy-range

    !modify eval array to contain tanh(E_i/(2*k_b*T)) instead of E_i
    temp=0.0d0
    if(temp<1.0d-5)then
        eval=1.0d0
    else
        eval=tanh(eval*0.5d0/temp)
    endif
    do j=1,size(del%delta,2)
        orb=(del%orb(j)-1)*nspin
#if 0
        do i_cell=1,ncell
            do iE=1,size(eval)
                Associate (u_up=>evec((i_cell-1)*ndim+orb+1     ,iE)&
                        & ,u_dn=>evec((i_cell-1)*ndim+orb+2     ,iE)&
                        & ,v_up=>evec((i_cell-1)*ndim+orb+1+nBdG,iE)&
                        & ,v_dn=>evec((i_cell-1)*ndim+orb+2+nBdG,iE))
                    !make this more efficient with vector operations, blas
                    del%delta(i_cell,j)=del%delta(i_cell,j)+(u_up*conjg(v_dn)+u_dn*conjg(v_up))*eval(ie)*del%V(j)*0.5d0 
                end associate
            enddo
        enddo
#else   
        !might make sense to employ a library like mkl vzmulbyconji
        do iE=1,size(eval)
            del%delta(:,j)=del%delta(:,j)+del%V(j)*0.5d0*eval(ie)*( &
 &               evec(orb+1:nBdG:ndim,iE)*conjg(evec(orb+2+nBdG:2*nBdG:ndim,iE)) &
 &            +  evec(orb+2:nBdG:ndim,iE)*conjg(evec(orb+1+nBdG:2*nBdG:ndim,iE)))
        enddo
#endif
    enddo
end subroutine

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


