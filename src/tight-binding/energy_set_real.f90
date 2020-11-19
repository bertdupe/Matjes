module m_energy_set_real
use m_derived_types, only: lattice
use m_tb_types
use m_get_table_nn,only :get_table_nn
implicit none
private
public set_Hr_dense_nc

contains

    subroutine set_Hr_dense_nc(lat,h_par,h_io,mode_mag,Hr)
        type(lattice),intent(in)                :: lat
        type(parameters_TB_Hsolve),intent(in)   :: h_par
        type(parameters_TB_IO_H),intent(in)     :: h_io
        real(8),intent(in)                      :: mode_mag(:,:)
        complex(8),allocatable,intent(inout)    :: Hr(:,:)

        if(.not. allocated(Hr))then
           allocate(Hr(h_par%dimH,h_par%dimH))
        endif
        if(size(Hr,1)/=h_par%dimH.or.size(Hr,2)/=h_par%dimH) STOP "Hr has wrong size"  !could easily reallocate, but this should never happen, I guess
        Hr=cmplx(0.0d0,0.0d0, 8)
        Call set_Hr_ee(lat,h_io,h_par,Hr)
        if(any(h_io%Jsd /= 0.0d0))then
            Call set_Jsd(h_par,mode_mag,h_io%Jsd,Hr)
        endif
    end subroutine 

    subroutine set_Hr_ee(lat,h_io,h_par,Hr)
        !extract the real space Hamiltonian Hr from the electronic part in energy
        type(lattice),intent(in)                :: lat
        type(parameters_TB_Hsolve),intent(in)   :: h_par
        type(parameters_TB_IO_H),intent(in)     :: h_io
        complex(8),allocatable,intent(inout)    :: Hr(:,:)

        integer                 :: dimH
        integer                 :: N_neighbours,Ncell,dim_mode
        integer, allocatable    :: indexNN(:),tableNN(:,:,:,:,:,:)
        integer                 :: i_cell,i_orb,i_s,i_neigh,i_vois
        integer                 :: ind_orb,i_other
        integer                 :: offset
        integer                 :: ilat_1(3),ilat_2(3),i1,i2
        complex(8)              :: H_entry

        N_neighbours = size(h_io%hopping,3)
        Call get_table_nn(lat,N_neighbours,indexNN,tableNN)
        dimH=h_par%dimH
        Ncell = lat%Ncell
        dim_mode=h_par%norb*h_par%nspin
        if(.not. allocated(Hr))then
           allocate(Hr(dimH,dimH),source=cmplx(0.0d0,0.0d0,8))
        endif
        !set onsite terms
        do i_orb=1,h_io%nb_orbitals
            ind_orb=(i_orb-1)*h_io%nb_spin
            do i_s=1,h_io%nb_spin
                if(h_io%onsite(i_s,i_orb)==0.0d0) cycle
                H_entry=cmplx(h_io%onsite(i_s,i_orb),0.0d0,8)
                do i_cell=1,Ncell
                    i1=i_s+ind_orb+(i_cell-1)*dim_mode
                    Hr(i1,i1)=Hr(i1,i1)+H_entry
                enddo
            enddo
        enddo
        !super ugly with tables, translate i_x,i_y,i_z to one index?, also loop order is questionable
        !hopping terms
        do i_neigh=1,N_neighbours
            offset=sum(indexNN(1:i_neigh-1))
            do i_orb=1,h_io%nb_orbitals
                ind_orb=(i_orb-1)*h_io%nb_spin
                do i_s=1,h_io%nb_spin
                    if(h_io%hopping(i_s,i_orb,i_neigh)==0.0d0) cycle
                    H_entry=cmplx(h_io%hopping(i_s,i_orb,i_neigh),0.0d0,8)
                    do i_cell=1,Ncell
                        ilat_1=lat%index_1_3(i_cell)
                        i1=i_s+ind_orb+(i_cell-1)*dim_mode
                        do i_vois=1,indexNN(i_neigh)
                            if(tableNN(5,i_vois,ilat_1(1),ilat_1(2),ilat_1(3),1)/=1) cycle
                            ilat_2=tableNN(1:3,i_vois+offset,ilat_1(1),ilat_1(2),ilat_1(3),1)
                            i_other=lat%index_m_1(ilat_2)
                            i2=i_s+ind_orb+(i_other-1)*dim_mode
                            Hr(i1,i2)=Hr(i1,i2)+H_entry
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine 

    subroutine set_Jsd(h_par,mode_mag,Jsd,Hr)
        !adds the Jsd coupling to a local real-space Hamiltonian Hr
        type(parameters_TB_Hsolve),intent(in)        ::  h_par
        real(8),intent(in)                          ::  Jsd(:)
        real(8),intent(in)                       ::  mode_mag(:,:)
        complex(8),allocatable,intent(inout)        ::  Hr(:,:)

        integer                 ::  i,j
        integer                 ::  i1,i2
        integer                 ::  a,b
        
        complex(8),allocatable  ::  add_Jsd(:,:)

        !initial checks
        if(h_par%nspin/=2) STOP "nspin has to be 2 to use Jsd-coupling"
        if(h_par%nsc/=1) STOP "nsc/=1 in set_jsd"
        if(size(Hr,1)/=h_par%dimH.or.size(Hr,2)/=h_par%dimH) STOP "Hr has wrong size" 
        if(size(Jsd) /= h_par%norb) STOP "JSD has wrong size"
        if(size(mode_mag,1)>3) ERROR STOP "will not work with nmag>1"

        allocate(add_Jsd(h_par%nsite,h_par%nsite))
        do i=1,h_par%ncell
            add_Jsd=cmplx(0.0d0,0.0d0,8)
            do j=1,h_par%norb
                i1=j*2-1
                i2=j*2
                !could be done be elegant with pauli matrices...
                add_Jsd(i1,i1)=Jsd(j)*cmplx( mode_mag(3,i), 0.0d0        ,8)
                add_Jsd(i2,i1)=Jsd(j)*cmplx( mode_mag(1,i), mode_mag(2,i),8)
                add_Jsd(i1,i2)=Jsd(j)*cmplx( mode_mag(1,i),-mode_mag(2,i),8)
                add_Jsd(i2,i2)=Jsd(j)*cmplx(-mode_mag(3,i), 0.0d0        ,8)
            enddo
            a=(i-1)*h_par%nsite+1
            b=i*h_par%nsite
            Hr(a:b,a:b) = Hr(a:b,a:b)+add_Jsd
        enddo 
    end subroutine 

end module
