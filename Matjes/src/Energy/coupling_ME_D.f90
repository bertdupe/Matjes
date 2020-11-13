module m_coupling_ME_D
use m_input_H_types, only: io_H_ME_D
use m_io_utils, only: get_parameter,get_coeff,number_nonzero_coeff,max_ind_variable
implicit none

private
public :: get_coupling_ME_D, read_ME_D_input

contains

subroutine read_ME_D_input(io_param,fname,io)
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_ME_D),intent(out)     :: io
    !internal
    integer :: max_entry
    real(8) :: const_mult

    max_entry=max_ind_variable(io_param,'ME_antisym_',fname)
    if (max_entry==0) then
        io%is_set=.false.
        return
    endif
    allocate(io%val(max_entry),source=0.0d0)
    call get_coeff(io_param,fname,'ME_antisym_',io%val)
    const_mult=-1.0d0 !minus one to be consistent with older version
    call get_parameter(io_param,fname,'c_ME',const_mult)
    io%val=const_mult*io%val
    io%is_set=maxval(abs(io%val))/=0.0d0
end subroutine

subroutine get_coupling_ME_D(Ham,io,tableNN,indexNN,lat)
    !get coupling in t_H Hamiltonian format
    use m_H_public
    use m_derived_types
    use m_setH_util,only: get_coo

    class(t_H),intent(inout)    :: Ham
    type(io_H_ME_D),intent(in)    :: io
    integer,intent(in)          :: tableNN(:,:,:,:,:,:) !!tableNN(5,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
    integer,intent(in)          :: indexNN(:)
    type(lattice),intent(in)    :: lat

    !local parameters
    !ME_parameters
    real(8), allocatable :: D(:)

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

    integer     :: ind1(4),ind2(4)
    real(8)     :: diff_pos(3)

    integer,allocatable :: connect(:,:)
    integer :: N_line
    integer :: ilat_1(3),ilat_2(3)
    integer :: i_x,i_y,i_z

    if(io%is_set)then
        if(lat%E%dim_mode==0) STOP "E-field has to be set when using coupling_ME"
        Call get_Htype(Ham_tmp)
        D=io%val
        Nshell=size(D)
        Ncell=product(lat%dim_lat)
        shape_tableNN=shape(tableNN)
        if(shape_tableNN(6)/=1) STOP "implement several mag atoms for ME-coupling"
        if(lat%M%dim_mode/=3) STOP "lat%M%dim_mode!=0, implement several mag atoms for ME-coupling"

        allocate(Htmp(lat%M%dim_mode,lat%E%dim_mode*lat%M%dim_mode))!local Hamiltonian modified for each shell/neighbor
        allocate(connect(2,Ncell),source=0) ! at most Ncell connections for each neighbor

        do i_sh=1,Nshell
            if(D(i_sh)==0.0d0) cycle
            
            Nvois=indexNN(i_sh)
            offset=sum(indexNN(1:i_sh-1))
            do i_vois=1,Nvois
            !get local Hamiltonian for given neighbor
                ind1=[1,1,1,1]
                ind2=tableNN(1:4,i_vois+offset,1,1,1,1)
                diff_pos=lat%pos_diff_ind(ind1,ind2) 
                where(abs(diff_pos)<norm2(diff_pos)*1.0d-8) diff_pos=0.0d0
                diff_pos=diff_pos/norm2(diff_pos)
                !explicitly assuming M only has dim_mode 3
                !fast index is M for second index of Htmp
                Htmp=0.0d0
                !use lagrange identity to write as
                !(mi.Ei)(mj.r)-(mj.Ei)(mi.r)

                !better with function resolving Mi*E
                !drop index of E in comment
                Htmp(2,1)= D(i_sh)*diff_pos(2)    ! Mi_x Mj_y E_x r_y
                Htmp(3,1)= D(i_sh)*diff_pos(3)    ! Mi_x Mj_z E_x r_z
                Htmp(1,5)= D(i_sh)*diff_pos(1)    ! Mi_y Mj_x E_y r_x
                Htmp(3,5)= D(i_sh)*diff_pos(3)    ! Mi_y Mj_z E_y r_z
                Htmp(1,9)= D(i_sh)*diff_pos(1)    ! Mi_z Mj_x E_z r_x
                Htmp(2,9)= D(i_sh)*diff_pos(2)    ! Mi_z Mj_y E_z r_y

                Htmp(1,2)=-D(i_sh)*diff_pos(2)    ! Mi_y Mj_x E_x r_y
                Htmp(1,3)=-D(i_sh)*diff_pos(3)    ! Mi_z Mj_x E_x r_z
                Htmp(2,4)=-D(i_sh)*diff_pos(1)    ! Mi_x Mj_y E_y r_x
                Htmp(2,6)=-D(i_sh)*diff_pos(3)    ! Mi_z Mj_y E_y r_z
                Htmp(3,7)=-D(i_sh)*diff_pos(1)    ! Mi_x Mj_z E_z r_x
                Htmp(3,8)=-D(i_sh)*diff_pos(2)    ! Mi_y Mj_z E_z r_y
                Call get_coo(Htmp,val_tmp,ind_tmp)

            !get lattice sites that have to be connected
                N_line=0
                connect=0
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

            !add hamiltonian to output Hamiltonian
                Call Ham_tmp%init_mult_2(connect(:,:N_line),val_tmp,ind_tmp,[1],[1,2],lat)
                deallocate(val_tmp,ind_tmp)
                Call Ham%add(Ham_tmp)
                Call Ham_tmp%destroy()
            enddo
        enddo
    endif

end subroutine 
end module 
