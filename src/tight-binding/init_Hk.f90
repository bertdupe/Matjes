module m_init_Hk
!module which contains the initialization for the tight-binding Hamiltonian in k-space (Hk_inp_t, slower type with all data in separate arrays)
!most routines in the module are only needed for getting the superconducting delta-parameters self-consistently, otherwise only the normal get_H of init_H.f90 is sufficient
use m_derived_types, only: lattice,k_grid_t
use m_H_tb_public
use m_tb_types ,only: parameters_TB_IO_H 
use m_ham_init_type ,only: parameters_ham_init 
use m_types_tb_h_inp 
use m_neighbor_type, only: neighbors
use m_init_H
use m_delta_onsite
use m_Hk
use m_kgrid_int, only: get_kmesh, kmesh_t 
use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none

private
public get_Hk_inp

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

    class(kmesh_t),allocatable              :: k_grid

    character(len=3)    :: i_str
    complex(8)          :: delta_sum(2)
    real(8)             :: diff
    integer             :: i

    !get constant non-superconducting terms
    Call get_H_TB(lat,h_io,H_inp_nc%H,.true.,diffR=H_inp_nc%diffR)
    !get initial superconducting term
    Call get_H_sc(lat,h_io,H_inp_sc%H,del,diffR=H_inp_sc%diffR)
    !set self-consistent k-grid
    Call get_kmesh(k_grid,lat,h_io%scf_kgrid,h_io%fname_kmesh)

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
        write(error_unit,'(2(A,E14.8))') "WARNING, self-consistent delta in k-space did not reach convergence criteria with ",diff,' instead of ',h_io%scf_diffconv
    endif
    Call H_inp%combine(H_inp_nc,H_inp_sc)
    Call H_inp_nc%destroy()
    Call H_inp_sc%destroy()

    Call del%print_file("delta_onsite_final.dat",lat)
    write(output_unit,'(A,2E16.8)') "Final real delta extremal values:", minval(del%delta%re),maxval(del%delta%re)
end subroutine 

subroutine get_delta_scf(H_inp,k_grid,h_io,del)
    class(Hk_inp_T)             :: H_inp
    class(kmesh_t),intent(in)   :: k_grid
    type(parameters_TB_IO_H),intent(in)     :: h_io 
    type(Hdelta),intent(inout)  :: del

    integer         ::  ik,Nk
    real(8)         ::  k(3)
    real(8)         ::  normalize

    del%delta=(0.0d0,0.0d0)
    Nk=k_grid%get_Nk()
    do ik=1,Nk  !MPI parallelization over this loop?
        k=k_grid%get_k(ik)
        Call add_delta_k_scf(H_inp,k,h_io,del)
    enddo
    del%delta=del%delta*k_grid%get_normalize()
end subroutine


subroutine add_delta_k_scf(H_inp,k,h_io,del)
! add delta contribution for a single k as called from get_delta_scf
! use equation 3.28 (Jian-Xin Zhu book, DOI: 10.1007/978-3-319-31314-6 ), using the symmetries only using positive energies (below eq.2.26) analogous to Eq.2.43 as used in real space 
    class(Hk_inp_T)             :: H_inp
    real(8),intent(in)          :: k(3)
    type(parameters_TB_IO_H),intent(in)     :: h_io 
    type(Hdelta),intent(inout)  :: del

    integer                     :: ndim,nspin, nBdG, ncell
    integer                     :: orb(2)
    integer                     :: i_cell
    real(8),allocatable         :: eval(:)
    complex(8),allocatable      :: evec(:,:)
    integer                     :: j,iE
    real(8)                     :: temp !temperature
    logical                     :: vanish

    nspin=h_io%nspin
    nBdG =h_io%norb*h_io%nspin*h_io%ncell
    ndim =h_io%norb*h_io%nspin
    ncell= size(del%delta,1)

    Call Hk_evec(H_inp,k,h_io,eval,evec)
    Call restict_solution_positive(eval,evec,h_io%scf_Ecut,vanish)
    if(vanish) return !no energy value in considered energy-range

    !modify eval array to contain tanh(E_i/(2*k_b*T)) instead of E_i
    !temp=0.0d0
    !if(temp<1.0d-5)then
        eval=1.0d0
    !else
    !    eval=tanh(eval*0.5d0/temp)
    !endif

    do j=1,size(del%delta,2)
        orb=(del%orb(:,j)-1)*nspin
#if 1
        do i_cell=1,ncell
            do iE=1,size(eval)
                Associate (u1_up=>evec((i_cell-1)*ndim+orb(1)+1     ,iE)&
                        & ,u1_dn=>evec((i_cell-1)*ndim+orb(1)+2     ,iE)&
                        & ,v1_up=>evec((i_cell-1)*ndim+orb(1)+1+nBdG,iE)&
                        & ,v1_dn=>evec((i_cell-1)*ndim+orb(1)+2+nBdG,iE)&
                        & ,u2_up=>evec((i_cell-1)*ndim+orb(2)+1     ,iE)&
                        & ,u2_dn=>evec((i_cell-1)*ndim+orb(2)+2     ,iE)&
                        & ,v2_up=>evec((i_cell-1)*ndim+orb(2)+1+nBdG,iE)&
                        & ,v2_dn=>evec((i_cell-1)*ndim+orb(2)+2+nBdG,iE))
                    !make this more efficient with vector operations, blas
                    del%delta(i_cell,j)=del%delta(i_cell,j)+(u1_up*conjg(v2_dn)+u2_dn*conjg(v1_up)+u1_dn*conjg(v2_up)+u2_up*conjg(v1_dn))*eval(ie)*del%V(j)*0.25d0  !(Eq.2.26)
                end associate
                !test using different symmetries
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
end module


