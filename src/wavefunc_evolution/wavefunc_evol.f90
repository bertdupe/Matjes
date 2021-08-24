module m_wavefunc_evol
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
use m_io_wavefunc_evol, only: io_wavefunc_evol
use m_derived_types, only : lattice,number_different_order_parameters
use m_derived_types, only : io_parameter,simulation_parameters
use m_hamiltonian_collection, only: hamiltonian
use mpi_basic, only: mpi_type
use m_H_tb_public
use m_constants, only:k_b
use m_ham_init_type, only: parameters_ham_init 
use m_TB_types, only: parameters_TB
implicit none
private
public wavefunc_evolution

contains
subroutine wavefunc_evolution(lat,comm)
    type(lattice), intent(inout)                :: lat
    type(mpi_type),intent(in)                   :: comm

    type(parameters_TB)         :: io_tb    !tight-binding Hamiltonian parameters
    type(io_wavefunc_evol)      :: io_e     !wave-function evolution io-parameters
    integer                     :: io       !io-unit

   if(comm%ismas)then
       write(output_unit,'(a)') 'entering into the wavefunction evoluation routine'
       if(comm%Np>1)then
          write(error_unit,'(//A)') "WARNING, running wavefunction_evolution with MPI parallelization which is not implemented."
          write(error_unit,'(A)')   "         Only one thead will to anything."
       endif
       !get io-input
       open(newunit=io,file='input')
       Call io_e%read_file('input',io)
       close(io)
       call io_tb%read_file('input')
       Call wavefunc_evolution_run(lat,io_e,io_tb,comm)
   endif
end subroutine


subroutine wavefunc_evolution_run(lat,io_e,io_TB,comm)
    !evolve wave-function by mutiplying the Hamiltonian to it, not sure if this is actually correct
    !so far everything is dense and thus rather slow
    use m_tightbinding_r, only: tightbinding_r
    use m_init_Hr
    use m_constants, only: hbar
    type(lattice), intent(inout)                :: lat
    type(io_wavefunc_evol),intent(inout)        :: io_e     !wave-function evolution io-parameters
    type(parameters_TB),intent(inout)           :: io_tb    !tight-binding Hamiltonian parameters
    type(mpi_type),intent(in)                   :: comm     !mpi-communicator for parallelization that is not implemented on this level

    class(H_tb),allocatable     :: tbH      !real-space Hamiltonian
    type(parameters_ham_init)   :: hinit    !type to get shape parameters from Hamiltonian input

    integer         :: nspin    !nspin=1 -> no spin, nspin=2-> spin polarized calculation
    integer         :: i_outer
    integer         :: i_inner
    integer         :: l
    integer         :: dimH     !dimension of Hamiltonian

    !4th order Runke-Kutta parameters
    integer,parameter   ::  Norder=4
    real(8),parameter   ::  c(Norder)=[0.0d0,0.5d0,0.5d0,1.0d0]
    real(8),parameter   ::  b(Norder)=[1.0d0,2.0d0,2.0d0,1.0d0]/6.0d0
    real(8),parameter   ::  a(Norder,Norder)= reshape([0.0d0, 0.0d0, 0.0d0, 0.0d0 &
                                                  ,0.5d0, 0.0d0, 0.0d0, 0.0d0 &
                                                  ,0.0d0, 0.5d0, 0.0d0, 0.0d0 &
                                                  ,0.0d0, 0.0d0, 1.0d0, 0.0d0]&
                                                  ,[Norder,Norder])
    complex(8),allocatable  :: wavefunc(:,:)    !temporary wave-function at each step of the Runke-Kutta iteration
    complex(8),allocatable  :: wavefunc_tmp(:)  !additional temporary wave-function array for iteration

    real(8)                 :: time, time_inner !local time after full update and within inner loop
    complex(8)              :: prefH    !prefactor for Hamiltonian

    real(8)                 :: energy   !energy
    integer                 :: io_dat, io_check!IO-units for output

    complex(8),allocatable  :: H_loc(:,:)       !Hamiltonian-dense (dimH,dimH) (needed for current operator)
    complex(8),allocatable  :: pos(:,:)         !position (dimH,dimH) 
    complex(8),allocatable  :: current(:,:,:)   !spin-dependent current operator
    complex(8),allocatable  :: pos_spn(:,:,:)   !spin-dependent position operator

    real(8),allocatable     :: current_av(:)    !current average  (nspin)
    real(8),allocatable     :: pos_av(:)        !position average (nspin)
    integer     ::  i


    !set TB_params(m_tb_params) and get real-space Hamiltonian tbH
    Call io_tb%init(lat)
    Call get_Hr(lat,io_tb%io_H,tbH)
    Call tbH%get_hinit(hinit)
    nspin=hinit%nspin

    !some weird initialization
    !Call init_wavefunc_EF(lat,tbH,300.0d0*k_b)
    !lat%w%all_modes_c=(0.0d0,0.0d0)
    !lat%w%all_modes_c(1)=(1.0d0,0.0d0)
    !if(nspin==2) lat%w%all_modes_c(2)=(1.0d0,0.0d0)


    !set some parameters
    dimH=size(lat%w%all_modes_c)
    prefH=cmplx(0.0d0,-1.0d0/hbar,8)

    select type(tbH)
    class is(H_tb_dense)
        allocate(H_loc(dimH,dimH))
        Call tbH%get_H(H_loc)
!        !add increasing potential (i.e. E-field)
!        do i=1,dimH
!            H_loc(i,i)=H_loc(i,i)+i*(-1.0d-3,0.0d0)
!        enddo
!        Call tbH%set_H(H_loc)
    class default
        STOP "IMPLEMENT GET_H for all tight-binding Hamiltonian or implement some sparse solution"
    end select

    !prepare operators and array results(1D chain)
    allocate(current(dimH,dimH,nspin),source=(0.0d0,0.0d0))
    allocate(pos_spn(dimH,dimH,nspin),source=(0.0d0,0.0d0))
    allocate(current_av(nspin),source=0.0d0)
    allocate(pos_av(nspin),source=0.0d0)

    !get position operator of chain, with imaginary part as in operator
    allocate(pos(dimH,dimH),source=(0.0d0,0.0d0)) 
    if(nspin==1)then
        do i=1,dimH
            pos(i,i)=(0.0d0,1.0d0)*i
        enddo
    else
        do i=1,dimH
            pos(i,i)=(0.0d0,1.0d0)*((i+1)/2)
        enddo
    endif

    !get spin-dependent position and current operator
    if(nspin==1)then
        pos_spn(:,:,1)=pos*(0.0d0,-1.0d0)
        current(:,:,1)=matmul(H_loc,pos)-matmul(pos,H_loc)
    else
        !split current and position operator into spin parts
        pos_spn(:,:,1)=pos*(0.0d0,-1.0d0)
        pos_spn(2::2,:,1)=(0.0d0,0.0d0)
        pos_spn(:,2::2,1)=(0.0d0,0.0d0)
        current(:,:,1)=matmul(H_loc,pos)-matmul(pos,H_loc)
        current(2::2,:,1)=(0.0d0,0.0d0)
        current(:,2::2,1)=(0.0d0,0.0d0)
        pos_spn(:,:,2)=pos*(0.0d0,-1.0d0)
        pos_spn(1::2,:,2)=(0.0d0,0.0d0)
        pos_spn(:,1::2,2)=(0.0d0,0.0d0)
        current(:,:,2)=matmul(H_loc,pos)-matmul(pos,H_loc)
        current(1::2,:,2)=(0.0d0,0.0d0)
        current(:,1::2,2)=(0.0d0,0.0d0)
    endif
    deallocate(H_loc, pos)

    time=0.d0
    open(newunit=io_dat  ,file='wavefunc_evol.dat')
    open(newunit=io_check,file='wavefunc_check.dat')
    allocate(wavefunc(dimH,Norder))
    allocate(wavefunc_tmp(dimH))
    do i_outer=1,io_e%duration

        !write info of current state
        if(mod(i_outer,io_e%Efreq)==1)then 
            wavefunc_tmp=lat%w%all_modes_c
            Call tbH%mult_r(lat%w%all_modes_c,wavefunc_tmp)
            energy=DOT_PRODUCT(lat%w%all_modes_c,wavefunc_tmp)
            write(io_check,*) time,energy,norm2(abs(lat%w%all_modes_c))
            do i=1,nspin
                wavefunc_tmp=matmul(current(:,:,i),lat%w%all_modes_c)
                current_av(i)=real(DOT_PRODUCT(wavefunc_tmp,lat%w%all_modes_c))
                wavefunc_tmp=matmul(pos_spn(:,:,i),lat%w%all_modes_c)
                pos_av(i)=real(DOT_PRODUCT(wavefunc_tmp,lat%w%all_modes_c))/DOT_PRODUCT(lat%w%all_modes_c,lat%w%all_modes_c)
            enddo
            write(io_dat,'(1000E18.8E3)') time,norm2(abs(lat%w%all_modes_c)), energy, pos_av,current_av,abs(lat%w%all_modes_c)**2
        endif

        !Runge-Kutta iteration step
        do i_inner=1,Norder
            time_inner=time+io_e%timestep*c(i_inner)
            wavefunc_tmp=lat%w%all_modes_c
            do l=1,i_inner-1
                wavefunc_tmp=wavefunc_tmp+io_e%timestep*a(l,i_inner)*wavefunc(:,l)
            enddo
            Call tbH%mult_r(wavefunc_tmp,wavefunc(:,i_inner),alpha=prefH)
        enddo
        do l=1,Norder
            lat%w%all_modes_c=lat%w%all_modes_c+io_e%timestep*b(l)*wavefunc(:,l) !update state
        enddo
        time=time+io_e%timestep  !update time

    enddo
    close(io_dat)
    close(io_check)
    write(output_unit,'(A)') "Finished wave-function evolution routine"
end subroutine

subroutine init_wavefunc_EF(lat,tbH,temp)
    !initialize wave-function based on Fermi energy, not really sensible
    type(lattice),intent(inout) :: lat
    class(H_tb),intent(in)      :: tbH
    real(8),intent(in)          :: temp !temperature

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
