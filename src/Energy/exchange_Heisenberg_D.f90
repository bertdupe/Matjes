module m_exchange_heisenberg_D
use m_input_H_types, only: io_H_D
implicit none
private
public read_D_input, get_exchange_D
contains
subroutine read_D_input(io_param,fname,io)
    use m_io_utils
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_D),intent(out)        :: io
    !internal
    integer         :: max_entry
    real(8)         :: const_mult

    max_entry=max_ind_variable(io_param,'DMI_',fname)
    if(max_entry==0)then
        io%is_set=.false.
        return
    endif
    allocate(io%val(max_entry),source=0.0d0)
    call get_coeff(io_param,fname,'DMI_',io%val)
    const_mult=-1.0d0   !-1 for consistency with older version
    call get_parameter(io_param,fname,'c_DMI',const_mult)
    io%val=io%val*const_mult
    io%is_set=maxval(abs(io%val))/=0.0d0
end subroutine


subroutine get_exchange_D(Ham,io,tableNN,indexNN,lat,DM_vector)
    !get coupling in t_H Hamiltonian format
    use m_H_public
    use m_derived_types
    use m_setH_util,only: get_coo

    class(t_H),intent(inout)    :: Ham
    type(io_H_D),intent(in)     :: io
    integer,intent(in)          :: tableNN(:,:,:,:,:,:) !!tableNN(5,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
    integer,intent(in)          :: indexNN(:)
    type(lattice),intent(in)    :: lat
    real(8), intent(in)         :: DM_vector(:,:,:)

    !ME_parameters
    real(8), allocatable :: D(:)

    !local Hamiltonian
    real(8),allocatable  :: Hloc(:,:),Hsh(:,:)
    !local Hamiltonian in coo format
    real(8),allocatable  :: val_tmp(:)
    integer,allocatable  :: ind_tmp(:,:)

    class(t_H),allocatable    :: Ham_tmp
    integer     :: Nshell,Ncell
    integer     :: shape_tableNN(6)
    integer     :: i_sh,i_vois
    integer     :: Nvois,offset

    integer,allocatable :: line(:,:)
    integer :: ilat_1(3),ilat_2(3)
    integer :: i_x,i_y,i_z

    if(io%is_set)then
        Call get_Htype(Ham_tmp)
        D=io%val
        Nshell=size(D)
        Ncell=lat%Ncell
        shape_tableNN=shape(tableNN)
        if(shape_tableNN(6)/=1) STOP "implement several mag atoms for exchange_J"
        if(lat%M%dim_mode/=3) STOP "lat%M%dim_mode!=0, implement several mag atoms for exchange_J"

        allocate(Hloc(lat%M%dim_mode,lat%M%dim_mode),source=0.0d0)!local Hamiltonian modified for each neighbor
        allocate(Hsh(lat%M%dim_mode,lat%M%dim_mode),source=0.0d0)!local Hamiltonian modified for each shell

        do i_sh=1,Nshell
            if(D(i_sh)==0.0d0) cycle
            Hsh=0.0d0
            Hsh(2,1)= D(i_sh)
            Hsh(3,1)=-D(i_sh)
            Hsh(1,2)=-D(i_sh)
            Hsh(3,2)= D(i_sh)
            Hsh(1,3)= D(i_sh)
            Hsh(2,3)=-D(i_sh)
            
            Nvois=indexNN(i_sh)
            if(allocated(line)) deallocate(line)
            allocate(line(1,Ncell),source=0) 
            offset=sum(indexNN(1:i_sh-1))
            do i_vois=1,Nvois
                if(norm2(DM_vector(i_vois+offset,:,1))<1.0d-8)then 
                    write(6,'(A)') "WARNING, skipping DMI-neighbors since the DM-vector is zero"
                    cycle
                endif
                Call get_loc_Ham_D(DM_vector(i_vois+offset,:,1),Hsh,Hloc)
                Call get_coo(Hloc,val_tmp,ind_tmp)
                do i_z=1,shape_tableNN(5)
                  do i_y=1,shape_tableNN(4)
                    do i_x=1,shape_tableNN(3)
                        if(tableNN(5,i_vois,i_x,i_y,i_z,1)/=1) cycle
                        ilat_1=[i_x,i_y,i_z]
                        ilat_2=tableNN(1:3,i_vois+offset,i_x,i_y,i_z,1)
                        line(1,lat%index_m_1(ilat_1))=lat%index_m_1(ilat_2)
                    enddo
                  enddo
                enddo
                !add hamiltonian to output Hamiltonian
                Call Ham_tmp%init_1(line,val_tmp,ind_tmp,[1,1],lat,2)
                deallocate(val_tmp,ind_tmp)
                Call Ham%add(Ham_tmp)
                Call Ham_tmp%destroy()
            enddo
        enddo
        Ham%desc="antisymmetric magnetic exchange"
    endif

end subroutine 

subroutine get_loc_Ham_D(D,Op_DMI,H)
    !subroutine that applied the local DMI-vector D to the neighbor-independent(within shell) Hamiltonion(OP_DMI)
    !to obtain the DMI operator for the neighbor(H)
    real(8), intent(in) :: D(:),Op_DMI(:,:)
    real(8), intent(inout) :: H(:,:)
    H(1,2)=Op_DMI(1,2)*D(3)
    H(1,3)=Op_DMI(1,3)*D(2)
    H(2,3)=Op_DMI(2,3)*D(1)
    
    H(2,1)=Op_DMI(2,1)*D(3)
    H(3,1)=Op_DMI(3,1)*D(2)
    H(3,2)=Op_DMI(3,2)*D(1)
end subroutine 


end module
