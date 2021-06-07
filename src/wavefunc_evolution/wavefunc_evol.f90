module m_wavefunc_evol
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
use m_derived_types, only : lattice,number_different_order_parameters
use m_derived_types, only : io_parameter,simulation_parameters
use m_hamiltonian_collection, only: hamiltonian
use m_dyna_io, only: dyna_input
use mpi_basic, only: mpi_type
use m_H_tb_public
use m_constants, only:k_b
implicit none
private
public wavefunc_evolution

contains
subroutine wavefunc_evolution(lat,ext_param,H,comm)
    type(lattice), intent(inout)                :: lat
    type(simulation_parameters), intent(inout)  :: ext_param
    type(hamiltonian),intent(inout)             :: H
    type(mpi_type),intent(in)                   :: comm

   if(comm%ismas)then
       write(output_unit,'(a)') 'entering into the wavefunction evoluation routine'
       if(comm%Np>1)then
          write(error_unit,'(//A)') "WARNING, running wavefunction_evolution with MPI parallelization which is not implemented."
          write(error_unit,'(A)')   "         Only one thead will to anything."
       endif
       Call wavefunc_evolution_run(lat,ext_param,H,comm)
   endif
end subroutine


subroutine wavefunc_evolution_run(lat,ext_param,H,comm)
    !wrapper to first initialize all spin-dynamics parameters and distribute on the different threads
    use m_TB_types, only: parameters_TB
    use m_tightbinding_r, only: tightbinding_r
    use m_rw_TB, only:  rw_TB
    use m_init_Hr
    type(lattice), intent(inout)                :: lat
    type(simulation_parameters), intent(inout)  :: ext_param
    type(hamiltonian),intent(inout)             :: H
    type(mpi_type),intent(in)                   :: comm

    type(parameters_TB)         :: tb_par
    class(H_tb),allocatable     :: tbH


    !read tight-binding io parameter from input and set TB_params(m_tb_params)
    call rw_TB(tb_par,'input')
    Call tb_par%init(lat)
    Call get_Hr(lat,tb_par%io_H,tbH)
    Call init_wavefunc_EF(lat,tbH,300.0d0*k_b)




    ERROR STOP "CONTINUE WAVEFUNC_EVOLUTION"
end subroutine

subroutine init_wavefunc_EF(lat,tbH,temp)
    type(lattice),intent(inout) :: lat
    class(H_tb),intent(in)      :: tbH
    real(8),intent(in)          :: temp

    real(8),allocatable         ::  eigval(:)
    complex(8),allocatable      ::  eigvec(:,:)
    real(8),parameter           ::  cutoff=30.0d0
    integer                     ::  i_min,i_max,i,j 

    lat%w%all_modes_c=cmplx(0.0d0,0.0d0,8)
    Call tbH%get_evec(eigval,eigvec)
    write(*,*) 'eigval'
    write(*,'(8F16.8)') eigval
    eigval=(eigval)/temp    !E_F=0.0

    i_min=1
    do i=1,size(eigval)
        if(eigval(i)>-cutoff)then
            i_min=i
            exit
        endif
    enddo

    i_max=size(eigval)
    do i=i_min,size(eigval)
        if(eigval(i)>cutoff)then
            i_max=max(1,i-1)
            exit
        endif
    enddo
    eigval(i_min:i_max)=exp(eigval(i_min:i_max))
    eigval(i_min:i_max)=1.0d0/(1.0d0+eigval(i_min:i_max))

    !add all deep occupied states
    do i=1,i_min-1
        lat%w%all_modes_c=lat%w%all_modes_c+eigvec(:,i)
    enddo

    !add all states near fermi-surface with fermi-dirac weight
    write(*,*) i_min,i_max
    do i=i_min,i_max
        lat%w%all_modes_c=lat%w%all_modes_c+eigval(i)*eigvec(:,i)
        write(*,*) i, eigval(i)
    enddo
    do i=1,size(eigvec,2)
        write(*,'(/I6,I6)') i
        write(*,'(2F16.8)') eigvec(:,i)
    enddo
    write(*,*)
    write(*,'(8F16.8)') sum(eigvec(:,1:8),2)
    write(*,'(8F16.8)') dot_product(sum(eigvec(:,1:8),2),sum(eigvec(:,1:8),2))
end subroutine
end module
