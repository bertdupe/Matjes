module m_wavefunc_evol
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
use m_derived_types, only : lattice,number_different_order_parameters
use m_derived_types, only : io_parameter,simulation_parameters
use m_hamiltonian_collection, only: hamiltonian
use m_dyna_io, only: dyna_input
use mpi_basic, only: mpi_type
use m_H_tb_public
use m_constants, only:k_b
use m_pos_op
use m_ham_init_type, only: parameters_ham_init 
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
    use m_init_Hr
    use m_constants, only: hbar
    type(lattice), intent(inout)                :: lat
    type(simulation_parameters), intent(inout)  :: ext_param
    type(hamiltonian),intent(inout)             :: H
    type(mpi_type),intent(in)                   :: comm

    type(parameters_TB)         :: tb_par
    class(H_tb),allocatable     :: tbH
    type(parameters_ham_init)   :: hinit

    integer         :: nspin
    integer         :: i_outer
    integer         :: i_inner
    integer         :: l
    integer         :: dimH

    integer,parameter   ::  Norder=4
    real(8)             ::  c(Norder)=[0.0d0,0.5d0,0.5d0,1.0d0]
    real(8)             ::  b(Norder)=[1.0d0,2.0d0,2.0d0,1.0d0]/6.0d0
    real(8)             ::  a(Norder,Norder)= reshape([0.0d0, 0.0d0, 0.0d0, 0.0d0 &
                                                  ,0.5d0, 0.0d0, 0.0d0, 0.0d0 &
                                                  ,0.0d0, 0.5d0, 0.0d0, 0.0d0 &
                                                  ,0.0d0, 0.0d0, 1.0d0, 0.0d0]&
                                                  ,[Norder,Norder])
    complex(8),allocatable  :: wavefunc(:,:)
    complex(8),allocatable  :: wavefunc_tmp(:)

    real(8)                 :: time, time_inner !local time after full update and within inner loop
    real(8)                 :: timestep !timestep in fs
    complex(8)              :: prefH    !prefactor for Hamiltonian

    real(8)                 :: energy
    integer                 :: io_dat

    complex(8),allocatable  :: H_loc(:,:)
    complex(8),allocatable  :: pos(:,:)
    complex(8),allocatable  :: current(:,:,:)
    complex(8),allocatable  :: pos_spn(:,:,:)   !pos only for spin up/down

    real(8),allocatable     :: current_av(:)
    real(8),allocatable     :: pos_av(:)
    integer     ::  i


    !read tight-binding io parameter from input and set TB_params(m_tb_params)
    if(comm%ismas)then
        call tb_par%read_file('input')
        Call tb_par%init(lat)
    endif
    Call tb_par%bcast(comm)


    Call get_Hr(lat,tb_par%io_H,tbH)
    Call tbH%get_hinit(hinit)
    nspin=hinit%nspin

    !Call init_wavefunc_EF(lat,tbH,300.0d0*k_b)

    lat%w%all_modes_c=(0.0d0,0.0d0)
    !lat%w%all_modes_c(101)=(1.0d0,0.0d0)
    !lat%w%all_modes_c(102)=(1.0d0,0.0d0)
    lat%w%all_modes_c(1)=(1.0d0,0.0d0)
!    lat%w%all_modes_c(2)=(1.0d0,0.0d0)


    timestep=1.0d-3
    dimH=size(lat%w%all_modes_c)
    prefH=cmplx(0.0d0,-1.0d0/hbar,8)


    select type(tbH)
    class is(H_tb_dense)
        allocate(H_loc(dimH,dimH))
        Call tbH%get_H(H_loc)
!        do i=1,dimH
!            H_loc(i,i)=H_loc(i,i)+i*(-1.0d-3,0.0d0)
!        enddo
!        Call tbH%set_H(H_loc)
    class default
        STOP "IMPLEMENT GET_H for all tight-binding Hamiltonian or implement some sparse solution"
    end select

    allocate(pos,mold=H_loc) 
    allocate(current(dimH,dimH,nspin),source=(0.0d0,0.0d0))
    allocate(pos_spn(dimH,dimH,nspin),source=(0.0d0,0.0d0))
    allocate(current_av(nspin),source=0.0d0)
    allocate(pos_av(nspin),source=0.0d0)
    pos=(0.0d0,0.0d0)
    if(nspin==1)then
        do i=1,dimH
            pos(i,i)=(0.0d0,1.0d0)*i
        enddo
    else
        do i=1,dimH
            pos(i,i)=(0.0d0,1.0d0)*((i+1)/2)
        enddo
    endif
!    do i=1,dimH-1
!        pos(i,i+1)=(0.0d0,0.1d0)
!        pos(i+1,i)=(0.0d0,0.1d0)
!    enddo
    !pos(dimH,1)=(0.00d0,0.1d0)
    !pos(1,dimH)=(0.00d0,0.1d0)
    if(nspin==1)then
        current(:,:,1)=matmul(H_loc,pos)-matmul(pos,H_loc)
        pos_spn(:,:,1)=pos*(0.0d0,-1.0d0)
    else
        current(:,:,1)=matmul(H_loc,pos)-matmul(pos,H_loc)
        current(2::2,:,1)=(0.0d0,0.0d0)
        current(:,2::2,1)=(0.0d0,0.0d0)
        pos_spn(:,:,1)=pos*(0.0d0,-1.0d0)
        pos_spn(2::2,:,1)=(0.0d0,0.0d0)
        pos_spn(:,2::2,1)=(0.0d0,0.0d0)
        current(:,:,2)=matmul(H_loc,pos)-matmul(pos,H_loc)
        current(1::2,:,2)=(0.0d0,0.0d0)
        current(:,1::2,2)=(0.0d0,0.0d0)
        pos_spn(:,:,2)=pos*(0.0d0,-1.0d0)
        pos_spn(1::2,:,2)=(0.0d0,0.0d0)
        pos_spn(:,1::2,2)=(0.0d0,0.0d0)
    endif

    time=0.d0
    open(newunit=io_dat,file='wavefunc_evol.dat')
    allocate(wavefunc(dimH,Norder))
    allocate(wavefunc_tmp(dimH))
    do i_outer=1,1000000
        do i_inner=1,Norder
            time_inner=time+timestep*c(i_inner)
            !update external t-parameters if neccessary

            wavefunc_tmp=lat%w%all_modes_c
            do l=1,i_inner-1
                wavefunc_tmp=wavefunc_tmp+timestep*a(l,i_inner)*wavefunc(:,l)
            enddo
            Call tbH%mult_r(wavefunc_tmp,wavefunc(:,i_inner),alpha=prefH)
        enddo
        do l=1,Norder
            lat%w%all_modes_c=lat%w%all_modes_c+timestep*b(l)*wavefunc(:,l)
!            lat%w%all_modes_c(size(lat%w%all_modes_c))=(0d0,0d0) !lose everything at end
        enddo
        time=time+timestep

        if(mod(i_outer,100)==1)then 
            wavefunc_tmp=lat%w%all_modes_c
            Call tbH%mult_r(lat%w%all_modes_c,wavefunc_tmp)
            energy=DOT_PRODUCT(lat%w%all_modes_c,wavefunc_tmp)
            write(io_dat,*) time,energy,norm2(abs(lat%w%all_modes_c))
            do i=1,nspin
                wavefunc_tmp=matmul(current(:,:,i),lat%w%all_modes_c)
                current_av(i)=real(DOT_PRODUCT(wavefunc_tmp,lat%w%all_modes_c))
                wavefunc_tmp=matmul(pos_spn(:,:,i),lat%w%all_modes_c)
                pos_av(i)=real(DOT_PRODUCT(wavefunc_tmp,lat%w%all_modes_c))/DOT_PRODUCT(lat%w%all_modes_c,lat%w%all_modes_c)
            enddo
            write(432,'(1000E18.8E3)') time,norm2(abs(lat%w%all_modes_c)), energy, pos_av,current_av,abs(lat%w%all_modes_c)**2
        endif


    enddo
    close(io_dat)



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
    do i=i_min,i_max
        lat%w%all_modes_c=lat%w%all_modes_c+eigval(i)*eigvec(:,i)
    enddo
end subroutine
end module