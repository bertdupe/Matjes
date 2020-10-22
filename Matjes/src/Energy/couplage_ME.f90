module m_couplage_ME
use m_symmetry_operators
use m_lattice, only : my_order_parameters
use m_Hamiltonian_variables, only : coeff_ham_inter_spec
use m_convert
implicit none
type(coeff_ham_inter_spec), target, public, protected :: ME

private
public :: get_ham_ME,get_number_EM_DMI,get_coupling_ME

contains


integer function get_number_EM_DMI(fname)
use m_io_files_utils
use m_io_utils
character(len=*), intent(in) ::fname
! internal variables
integer :: io

io=open_file_read(fname)
get_number_EM_DMI=count_variables(io,'ME_antisym_',fname)
call close_file(fname,io)

end function get_number_EM_DMI


subroutine get_ME_input(fname,c_ham,sym,asym,N_me)
    use m_io_files_utils
    use m_io_utils
    character(len=*), intent(in)        ::  fname
    real(8), allocatable,intent(out)    ::  sym(:),asym(:)
    real(8),intent(out)                 ::  c_ham
    integer,intent(out)                 ::  N_me

    integer :: neighbor_sym,io_param,neighbor_asym

    neighbor_sym=0
    neighbor_asym=0
    io_param=open_file_read(fname)
    c_ham=1.0d0
    call get_parameter(io_param,fname,'c_ME',c_ham)
    !
    ! count the ME coefficients if present
    !
    neighbor_sym=max_ind_variable(io_param,'ME_sym_',fname)
    if (neighbor_sym.ne.0) then
       allocate(sym(neighbor_sym),source=0.0d0)
       call get_coeff(io_param,fname,'ME_sym_',sym)
       neighbor_sym=number_nonzero_coeff(sym,'ME_sym')
    endif
    
    neighbor_asym=max_ind_variable(io_param,'ME_antisym_',fname)
    if (neighbor_asym.ne.0) then
       allocate(asym(neighbor_asym),source=0.0d0)
       call get_coeff(io_param,fname,'ME_antisym_',asym)
       neighbor_asym=number_nonzero_coeff(asym,'ME_antisym')
    endif
    
    call close_file(fname,io_param)
    
    N_ME=max(neighbor_sym,neighbor_asym)
    if (N_ME==0) then
       return
    else
       ME%i_exist=.true.
       write(6,'(a)') 'WARNING!!!! ME coupling found. You need E_ext different from 0'
    endif

end subroutine


subroutine get_ham_ME(fname,dim_ham)
    integer, intent(in) :: dim_ham
    character(len=*), intent(in) ::fname
    ! internal
    integer :: N_ME
    real(kind=8), allocatable :: ME_local_sym(:),ME_local_antisym(:)
    real(8) :: c_ham
    ! magnetization
    integer :: x_start,x_end
    ! electric field
    integer :: y_start,y_end
    ! slope
    integer :: i,j
    character(len=50) :: form
    
    ME%name='magnetoelectric'
    ME%order=3
    
    Call get_ME_input(fname,c_ham,ME_local_sym,ME_local_antisym,N_ME)
    ME%c_ham=c_ham
    
    allocate(ME%ham(N_ME))
    do i=1,N_ME
       allocate(ME%ham(i)%H(dim_ham,dim_ham**2))
       ME%ham(i)%H=0.0d0
    enddo
    
    call get_borders('magnetic',x_start,x_end,'Efield',y_start,y_end,my_order_parameters)
#if 0
    ! get the symmetric ME effect
    do i=1,n_ME
      ! get diagonal terms
      ME%ham(i)%H(y_start,x_start)=ME_local_sym(i)*ME%c_ham
      ME%ham(i)%H(y_start+1,x_start+dim_ham+1)=ME_local_sym(i)*ME%c_ham
      ME%ham(i)%H(y_start+2,x_start+2*(dim_ham+1))=ME_local_sym(i)*ME%c_ham
    
      ! get the off-diagonal terms
    !  ME%ham(i)%H(y_start,x_start+1)=ME_local_antisym(i)*ME%c_ham
      ME%ham(i)%H(y_start+1,x_start+1)=ME_local_antisym(i)*ME%c_ham
    
    !  ME%ham(i)%H(y_start,x_start+2)=-ME_local_antisym(i)*ME%c_ham
      ME%ham(i)%H(y_start+2,x_start+2)=-ME_local_antisym(i)*ME%c_ham
    
    !  ME%ham(i)%H(y_start,x_start+dim_ham)=-ME_local_antisym(i)*ME%c_ham
      ME%ham(i)%H(y_start+1,x_start+dim_ham)=-ME_local_antisym(i)*ME%c_ham
    
    !  ME%ham(i)%H(y_start+1,x_start+dim_ham+2)=ME_local_antisym(i)*ME%c_ham
    !  ME%ham(i)%H(y_start+2,x_start+dim_ham+2)=ME_local_antisym(i)*ME%c_ham
    
    !  ME%ham(i)%H(y_start,x_start+2*dim_ham)=ME_local_antisym(i)*ME%c_ham
      ME%ham(i)%H(y_start+2,x_start+2*dim_ham)=ME_local_antisym(i)*ME%c_ham
    
    !  ME%ham(i)%H(y_start+1,x_start+2*dim_ham+1)=-ME_local_antisym(i)*ME%c_ham
    !  ME%ham(i)%H(y_start+2,x_start+2*dim_ham+1)=-ME_local_antisym(i)*ME%c_ham
    enddo
    
    form=convert('(',dim_ham,'(f12.8,2x))')
    
    write(6,'(a)') ''
    write(6,'(a)') 'Magnetoelectric Hamiltonian of order 3'
    do i=1,N_ME
      write(6,'(a,I3)') 'Shell  ',i
      do j=1,dim_ham**2
        write(6,form) ME%ham(i)%H(:,j)
      enddo
      write(6,'(a)') ''
    enddo
#endif
    
end subroutine get_ham_ME

subroutine get_coupling_ME(Ham,tableNN,indexNN,lat)
    !get coupling  in t_H Hamiltonian format
    !so far ME has to be set before
    use m_Htype_gen
    use m_derived_types
    use m_setH_util,only: get_coo

    class(t_H),intent(inout)    :: Ham
    integer, intent(in)         :: tableNN(:,:,:,:,:,:) !!tableNN(5,N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(my_motif%i_mom)
    integer, intent(in)         :: indexNN(:)
    type(lattice),intent(in)    :: lat

    !local parameters
    !ME_parameters
    real(8)              :: c_ham
    real(8), allocatable :: sym(:),asym(:)
    integer              :: N_ME

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

    if(ME%i_exist)then
        Call get_ME_input('input',c_ham,sym,asym,N_ME)
        sym=sym*c_ham;asym=asym*c_ham

        Call get_Htype(Ham_tmp)
        Nshell=max(size(sym),size(asym))
        Ncell=product(lat%dim_lat)
        shape_tableNN=shape(tableNN)
        if(shape_tableNN(6)/=1) STOP "implement several mag atoms for ME-coupling"
        if(lat%M%dim_mode/=3) STOP "lat%M%dim_mode!=0, implement several mag atoms for ME-coupling"

        !size so far only terms like E \cdot M -> same size as MxM Hamiltonian 
        allocate(Htmp(lat%M%dim_mode,lat%E%dim_mode*lat%M%dim_mode))!local Hamiltonian modified for each shell/neighbor
        allocate(connect(2,Ncell),source=0) ! at most Ncell connections for each neighbor

        do i_sh=1,Nshell
            l_sym=.False.;l_asym=.False.
            if(i_sh<=size( sym)) l_sym = sym(i_sh)/=0.0d0
            if(i_sh<=size(asym)) l_asym=asym(i_sh)/=0.0d0
            if(.not.(l_sym.or.l_asym)) cycle
            
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
                if(l_sym)then 
                    Htmp(1,1)=sym(i_sh)     ! Mi_x Mj_x E_x
                    Htmp(2,5)=sym(i_sh)     ! Mi_y Mj_y E_y
                    Htmp(3,9)=sym(i_sh)     ! Mi_z Mj_z E_z
                endif
                if(l_asym)then
                    !use lagrange identity to write as
                    !(mi.Ei)(mj.r)-(mj.Ei)(mi.r)

                    !better with function resolving Mi*E
                    !drop index of E in comment
                    Htmp(2,1)= asym(i_sh)*diff_pos(2)    ! Mi_x Mj_y E_x r_y
                    Htmp(3,1)= asym(i_sh)*diff_pos(3)    ! Mi_x Mj_z E_x r_z
                    Htmp(1,5)= asym(i_sh)*diff_pos(1)    ! Mi_y Mj_x E_y r_x
                    Htmp(3,5)= asym(i_sh)*diff_pos(3)    ! Mi_y Mj_z E_y r_z
                    Htmp(1,9)= asym(i_sh)*diff_pos(1)    ! Mi_z Mj_x E_z r_x
                    Htmp(2,9)= asym(i_sh)*diff_pos(2)    ! Mi_z Mj_y E_z r_y

                    Htmp(1,2)=-asym(i_sh)*diff_pos(2)    ! Mi_y Mj_x E_x r_y
                    Htmp(1,3)=-asym(i_sh)*diff_pos(3)    ! Mi_z Mj_x E_x r_z
                    Htmp(2,4)=-asym(i_sh)*diff_pos(1)    ! Mi_x Mj_y E_y r_x
                    Htmp(2,6)=-asym(i_sh)*diff_pos(3)    ! Mi_z Mj_y E_y r_z
                    Htmp(3,7)=-asym(i_sh)*diff_pos(1)    ! Mi_x Mj_z E_z r_x
                    Htmp(3,8)=-asym(i_sh)*diff_pos(2)    ! Mi_y Mj_z E_z r_y

                endif
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


end module m_couplage_ME
