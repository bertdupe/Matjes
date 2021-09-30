module m_init_Hr
!module which contains the initialization for the tight-binding Hamiltonian in r-space (H_tb_coo)
!most routines in the module are only needed for getting the superconducting delta-parameters self-consistently, otherwise only the normal get_H of init_H.f90 is sufficient
use m_derived_types, only: lattice
use m_H_tb_public
use m_tb_types ,only: parameters_TB_IO_H
use m_types_tb_h_inp 
use m_neighbor_type, only: neighbors
use m_init_H
use m_delta_onsite
use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none
private
public get_Hr, get_delta

contains


subroutine get_Hr(lat,h_io,H)
    !main routine to get the real-space Hamiltonian (H)
    type(lattice),intent(in)                :: lat
    type(parameters_TB_IO_H),intent(in)     :: h_io
    class(H_tb),allocatable,intent(out)     :: H

    type(Hdelta)                            :: del
    type(H_tb_coo),allocatable              :: H_list(:)
    integer         ::  i

    if(h_io%use_scf)then 
        !use self-conistent superconducting delta
        del=h_io%del
        Call get_Hr_conv(lat,h_io,H,del)
    else
        !get Hamiltonian without self-consistent delta (this should be the normal case)
        Call set_H(H,h_io)  !allocates H to the wanted implementation of H_tb

        Call get_H(lat,h_io,H_list) !gets the different interactions in separate Hamiltonians (init_H.f90)
        do i=1,size(H_list)
            Call H%add(H_list(i))   !add all the Hamiltonians to the output H
        enddo
    endif
end subroutine

subroutine get_delta(lat,h_io,del)
    !subroutine to only get self-consistent delta (necessary only for k-space calculation which uses the self-consistent delta obtained from real-space)
    type(lattice),intent(in)                :: lat
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(Hdelta),intent(out)                :: del

    class(H_tb),allocatable :: H
    if(h_io%use_scf)then !use self-conistent version
        del=h_io%del
        Call get_Hr_conv(lat,h_io,H,del)
    else
        STOP "It does not make sense to call this routine if no self-consistent delta is considered"
    endif
end subroutine

subroutine get_Hr_conv(lat,h_io,H,del)
    type(lattice),intent(in)                :: lat
    type(parameters_TB_IO_H),intent(in)     :: h_io 
    class(H_tb),allocatable,intent(out)     :: H
    type(Hdelta),intent(inout)              :: del

    type(H_tb_coo),allocatable              :: H_nc_tmp(:)
    class(H_tb),allocatable                 :: H_nc
    type(H_tb_coo),allocatable              :: H_sc_tmp(:)
    class(H_tb),allocatable                 :: H_sc
    class(H_tb),allocatable                 :: H_sum
    integer         ::  i,j
    complex(8)          :: delta_sum(2)
    real(8)             :: diff
    character(len=3)    :: i_str

    !get constant non-superconducting terms
    Call get_H_TB(lat,h_io,H_nc_tmp,.true.)
    Call set_H(H_nc,h_io)
    do i=1,size(H_nc_tmp)
        Call H_nc%add(H_nc_tmp(i))
        Call H_nc_tmp(i)%destroy()
    enddo
    deallocate(H_nc_tmp)

    !get initial superconducting term
    Call get_H_sc(lat,h_io,H_sc_tmp,del)
    Call set_H(H_sc,h_io)
    do i=1,size(H_sc_tmp)
        Call H_sc%add(H_sc_tmp(i))
        Call H_sc_tmp(i)%destroy()
    enddo

    !self-consistent loop
    Call set_H(H_sum,h_io)
    Call H_sum%add(H_nc)
    Call H_sum%add(H_sc)
    delta_sum=(0.0d0,0.0d0)
    diff=h_io%scf_diffconv*2.0d0
    do i=1,h_io%scf_loopmax
        if(diff<h_io%scf_diffconv)then
            write(output_unit,'(A,E16.8)') "Reached set convergence criterium for delta-convergence:",h_io%scf_diffconv
            exit
        endif
        Call get_delta_scf(H_sum,h_io%scf_Ecut,del)
        delta_sum(2)=sum(del%delta)/size(del%delta)
        
        deallocate(H_sc_tmp)
        Call get_H_sc(lat,h_io,H_sc_tmp,del)
        if(h_io%scf_print)then
            write(i_str,'(I0.3)') i
            Call del%print_file("delta_onsite_"//i_str//".dat",lat)
        endif
        Call H_sc%destroy()
        do j=1,size(H_sc_tmp)
            Call H_sc%add(H_sc_tmp(j))
            Call H_sc_tmp(j)%destroy()
        enddo
        Call H_nc%copy(H_sum)
        Call H_sum%add(H_sc)

        diff=abs(delta_sum(2)-delta_sum(1))
        delta_sum(1)=delta_sum(2)
        write(output_unit,'(2(A,I6)3XA,2E16.8,3XA,E16.8)') 'Finished delta scf loop',i,' of',h_io%scf_loopmax,'  delta-av sum :',delta_sum(1), '  delta-av diff:',diff
    enddo
    Call del%print_file("delta_onsite_final.dat",lat)
    write(output_unit,'(A,2E16.8)') "Final real delta extremal values:", minval(del%delta%re),maxval(del%delta%re)

    allocate(H,mold=H_sum)
    Call H_sum%mv(H)
end subroutine

subroutine get_delta_scf(H_in,emax,del)
    !so far only working for onsite Hdelta
    class(H_tb),allocatable     :: H_in
    type(Hdelta),intent(inout)  :: del
    real(8),intent(in)          :: emax

    complex(8),allocatable      :: evec(:,:)
    real(8),allocatable         :: eval(:)
    integer                     :: ndim,nspin, nBdG
    integer                     :: orb(2)
    integer                     :: i_cell
    integer :: j,iE
    real(8)                     :: temp !temperature


    Call H_in%get_evec(eval,evec)
    Call restict_solution_positive(eval,evec,emax)
    ndim =H_in%ndim
    nspin=H_in%nspin
    nBdG=H_in%norb*H_in%nspin*H_in%ncell

    temp=0.0d0
    if(temp<1.0d-5)then
        eval=1.0d0
    else
        eval=tanh(eval*0.5d0/temp)
    endif

    ERROR STOP "UPDATE TO DIFFERENT ORBITALS"
    !del%delta=(0.0d0,0.0d0)
    !do j=1,size(del%delta,2)
    !    orb=(del%orb(:,j)-1)*nspin
    !    do i_cell=1,size(del%delta,1)
    !        do iE=1,size(eval)
    !            Associate (u_up=>evec((i_cell-1)*ndim+orb+1     ,iE)&
    !                    & ,u_dn=>evec((i_cell-1)*ndim+orb+2     ,iE)&
    !                    & ,v_up=>evec((i_cell-1)*ndim+orb+1+nBdG,iE)&
    !                    & ,v_dn=>evec((i_cell-1)*ndim+orb+2+nBdG,iE))
    !                !make this more efficient with vector operations, blas
    !                del%delta(i_cell,j)=del%delta(i_cell,j)+(u_up*conjg(v_dn)+u_dn*conjg(v_up))*eval(ie) !(eq. 2.43 Jian-Xin Zhu book, DOI: 10.1007/978-3-319-31314-6 )
    !            end associate
    !        enddo
    !    enddo
    !    del%delta(:,j)=del%delta(:,j)*del%V(j)*0.5d0
    !enddo
end subroutine

end module

