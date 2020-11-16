module m_coupling_ME_J
use m_input_H_types, only: io_H_ME_J
use m_io_utils, only: get_parameter,get_coeff,number_nonzero_coeff,max_ind_variable
implicit none

private
public :: get_coupling_ME_J, read_ME_J_input

contains

subroutine read_ME_J_input(io_param,fname,io)
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_ME_J),intent(out)       :: io
    !internal
    integer :: max_entry
    real(8) :: const_mult

    max_entry=max_ind_variable(io_param,'ME_sym_',fname)
    if (max_entry==0) then
        io%is_set=.false.
        return
    endif
    allocate(io%val(max_entry),source=0.0d0)
    call get_coeff(io_param,fname,'ME_sym_',io%val)
    const_mult=-1.0d0 !minus one to be consistent with older version
    call get_parameter(io_param,fname,'c_ME',const_mult)
    io%val=const_mult*io%val
    io%is_set=maxval(abs(io%val))/=0.0d0
end subroutine

subroutine get_coupling_ME_J(Ham,io,tableNN,indexNN,lat)
    !get coupling in t_H Hamiltonian format
    use m_H_public
    use m_derived_types
    use m_setH_util,only: get_coo

    class(t_H),intent(inout)    :: Ham
    type(io_H_ME_J),intent(in)  :: io
    integer,intent(in)          :: tableNN(:,:,:,:,:,:) !!tableNN(5,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
    integer,intent(in)          :: indexNN(:)
    type(lattice),intent(in)    :: lat

    !local parameters
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
    integer     :: Nvois,offset

    real(8)     :: diff_pos(3)

    integer,allocatable :: connect(:,:)
    integer :: N_line
    integer :: ilat_1(3),ilat_2(3)
    integer :: i_x,i_y,i_z

    if(io%is_set)then
        if(lat%E%dim_mode==0) STOP "E-field has to be set when using coupling_ME"
        Call get_Htype(Ham_tmp)
        J=io%val
        Nshell=size(J)
        Ncell=lat%Ncell
        shape_tableNN=shape(tableNN)
        if(shape_tableNN(6)/=1) STOP "implement several mag atoms for ME-coupling"
        if(lat%M%dim_mode/=3) STOP "lat%M%dim_mode!=0, implement several mag atoms for ME-coupling"
        allocate(Htmp(lat%M%dim_mode,lat%E%dim_mode*lat%M%dim_mode))!local Hamiltonian modified for each shell/neighbor

        do i_sh=1,Nshell
            if(J(i_sh)==0.0d0) cycle
            Htmp=0.0d0
            Htmp(1,1)=J(i_sh)     ! Mi_x Mj_x E_x
            Htmp(2,5)=J(i_sh)     ! Mi_y Mj_y E_y
            Htmp(3,9)=J(i_sh)     ! Mi_z Mj_z E_z
            Call get_coo(Htmp,val_tmp,ind_tmp)
            Nvois=indexNN(i_sh)
            if(allocated(connect)) deallocate(connect)
            allocate(connect(2,Ncell*Nvois),source=0) ! at most Ncell*Nvois connections for each neighbor
            offset=sum(indexNN(1:i_sh-1))
            N_line=0
            do i_vois=1,Nvois
            !get lattice sites that have to be connected
                do i_z=1,shape_tableNN(5)
                  do i_y=1,shape_tableNN(4)
                    do i_x=1,shape_tableNN(3)
                        if(tableNN(5,i_vois,i_x,i_y,i_z,1)/=1) cycle
                        N_line=N_line+1
                        ilat_1=[i_x,i_y,i_z]
                        connect(1,N_line)=lat%index_m_1(ilat_1)
                        ilat_2=tableNN(1:3,i_vois+offset,i_x,i_y,i_z,1)
                        connect(2,N_line)=lat%index_m_1(ilat_2)
                    enddo
                  enddo
                enddo
            enddo
            !add hamiltonian to output Hamiltonian
            Call Ham_tmp%init_mult_2(connect(:,:N_line),val_tmp,ind_tmp,[1],[1,2],lat)
            deallocate(val_tmp,ind_tmp)
            Call Ham%add(Ham_tmp)
            Call Ham_tmp%destroy()
        enddo
        Ham%desc="symmetric magnetoelectric coupling"
    endif

end subroutine 


end module 
