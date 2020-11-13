module m_exchange_heisenberg_J
use m_input_H_types, only: io_H_J
implicit none
private
public read_J_input, get_exchange_J
contains
subroutine read_J_input(io_param,fname,io)
    use m_io_utils
    use m_input_H_types, only: io_H_aniso
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_J),intent(out)        :: io
    !internal
    integer         :: max_entry
    real(8)         :: const_mult

    max_entry=max_ind_variable(io_param,'J_',fname)
    if(max_entry==0)then
        io%is_set=.false.
        return
    endif
    io%is_set=.true.
    allocate(io%val(max_entry),source=0.0d0)
    call get_coeff(io_param,fname,'J_',io%val)
    const_mult=1.0d0
    call get_parameter(io_param,fname,'c_Jij',const_mult)
    io%val=io%val*const_mult
end subroutine


subroutine get_exchange_J(Ham,io,tableNN,indexNN,lat)
    !get coupling in t_H Hamiltonian format
    use m_H_public
    use m_derived_types
    use m_setH_util,only: get_coo

    class(t_H),intent(inout)    :: Ham
    type(io_H_J),intent(in)     :: io
    integer,intent(in)          :: tableNN(:,:,:,:,:,:) !!tableNN(5,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
    integer,intent(in)          :: indexNN(:)
    type(lattice),intent(in)    :: lat

    !ME_parameters
    real(8), allocatable :: J(:)

    !local Hamiltonian
    real(8),allocatable  :: Htmp(:,:)
    !local Hamiltonian in coo format
    real(8),allocatable  :: val_tmp(:)
    integer,allocatable  :: ind_tmp(:,:)

    class(t_H),allocatable    :: Ham_tmp
    integer     :: Nshell,Ncell
    integer     :: shape_tableNN(6)
    integer     :: i_sh,i_vois
    logical     :: l_sym,l_asym 
    integer     :: Nvois,offset

    integer,allocatable :: line(:,:)
    integer :: ilat_1(3),ilat_2(3)
    integer :: i_x,i_y,i_z

    if(io%is_set)then
        Call get_Htype(Ham_tmp)
        J=io%val
        Nshell=size(J)
        Ncell=lat%Ncell
        shape_tableNN=shape(tableNN)
        if(shape_tableNN(6)/=1) STOP "implement several mag atoms for exchange_J"
        if(lat%M%dim_mode/=3) STOP "lat%M%dim_mode!=0, implement several mag atoms for exchange_J"

        allocate(Htmp(lat%M%dim_mode,lat%M%dim_mode))!local Hamiltonian modified for each shell/neighbor

        do i_sh=1,Nshell
            if(J(i_sh)==0.0d0) cycle
            Htmp=0.0d0
            Htmp(1,1)=J(i_sh)
            Htmp(2,2)=J(i_sh)
            Htmp(3,3)=J(i_sh)
            Call get_coo(Htmp,val_tmp,ind_tmp)
            
            Nvois=indexNN(i_sh)
            if(allocated(line)) deallocate(line)
            allocate(line(Nvois,Ncell),source=0) 
            offset=sum(indexNN(1:i_sh-1))
            do i_vois=1,Nvois
                do i_z=1,shape_tableNN(5)
                  do i_y=1,shape_tableNN(4)
                    do i_x=1,shape_tableNN(3)
                        if(tableNN(5,i_vois,i_x,i_y,i_z,1)/=1) cycle
                        ilat_1=[i_x,i_y,i_z]
                        ilat_2=tableNN(1:3,i_vois+offset,i_x,i_y,i_z,1)
                        line(i_vois,lat%index_m_1(ilat_1))=lat%index_m_1(ilat_2)
                    enddo
                  enddo
                enddo
            enddo
            !add hamiltonian to output Hamiltonian
            Call Ham_tmp%init_1(line,val_tmp,ind_tmp,[1,1],lat)
            deallocate(val_tmp,ind_tmp)
            Call Ham%add(Ham_tmp)
            Call Ham_tmp%destroy()
        enddo
    endif

end subroutine 

end module
