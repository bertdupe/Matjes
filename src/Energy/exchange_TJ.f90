module m_exchange_TJ
use m_input_H_types, only: io_H_TJ
implicit none
private
public read_TJ_input, get_exchange_TJ
contains
subroutine read_TJ_input(io_param,fname,io)
    use m_io_utils
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_TJ),intent(out)        :: io
    !internal
    integer         :: max_entry
    real(8)         :: const_mult

    max_entry=max_ind_variable(io_param,'TJ_',fname)
    if(max_entry==0)then
        io%is_set=.false.
        return
    endif
    allocate(io%val(max_entry),source=0.0d0)
    call get_coeff(io_param,fname,'TJ_',io%val)
    const_mult=-1.0d0
    call get_parameter(io_param,fname,'c_TJ',const_mult)
    io%val=io%val*const_mult
    io%is_set=maxval(abs(io%val))/=0.0d0
end subroutine


subroutine get_exchange_TJ(Ham,io,tableNN,indexNN,lat)
    !get coupling in t_H Hamiltonian format
    use m_H_public
    use m_derived_types
    use m_setH_util,only: get_coo

    class(t_H),intent(inout)    :: Ham
    type(io_H_TJ),intent(in)    :: io
    integer,intent(in)          :: tableNN(:,:,:,:,:,:) !!tableNN(5,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
    integer,intent(in)          :: indexNN(:)
    type(lattice),intent(in)    :: lat

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

    integer,allocatable :: connect(:,:)
    integer :: ilat_1(3),ilat_2(3)
    integer :: i_x,i_y,i_z
    integer :: N_line

    if(io%is_set)then
        if(lat%T%dim_mode==0) STOP "temperature-field has to be set when using exchange_TJ"
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
            if(allocated(connect)) deallocate(connect)
            allocate(connect(2,Ncell*Nvois),source=0) ! at most Ncell*Nvois connections for each neighbor
            offset=sum(indexNN(1:i_sh-1))
            N_line=0
            do i_vois=1,Nvois
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
            Call Ham_tmp%init_mult_2(connect(:,:N_line),val_tmp,ind_tmp,[1],[1,4,4],lat)
            deallocate(val_tmp,ind_tmp)
            Call Ham%add(Ham_tmp)
            Call Ham_tmp%destroy()
        enddo
        Ham%desc="T^2 M^2 exchange"
    endif

end subroutine 

end module
