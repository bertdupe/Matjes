module m_energy_set_real
use m_energy_commons, only : energy
use m_basic_types, only : vec_point

private
public set_Hr,H_add_Jsd


contains

    subroutine set_Hr(dimH,Hr,Tb_ext)
        !extract the real space Hamiltonian Hr from the electronic part in energy
        integer,intent(in)           ::  dimH
        complex(8),intent(inout)     ::  Hr(dimH,dimH)
        integer,intent(in)           ::  TB_ext(2)

        integer                 ::  N_neighbours,N_cells,dim_mode
        integer                 ::  i,j,k

        N_neighbours = size( energy%line, 1 )
        N_cells = size( energy%line, 2 )
        dim_mode=Tb_ext(2)-Tb_ext(1)+1
        Hr = cmplx(0.0d0,0.0d0, kind=8)
        do i=1, N_cells
            do j=1, N_neighbours
                k = energy%line(j, i) 
                Hr((i-1)*dim_mode+1:i*dim_mode, (k-1)*dim_mode+1:k*dim_mode) = energy%value(j,k)%order_op(1)%Op_loc(TB_ext(1):TB_ext(2),TB_ext(1):TB_ext(2))
            enddo
        enddo
    end subroutine 

    subroutine H_add_Jsd(dimH,Hr,Tb_ext,mode_mag,Jsd)
        !adds the Jsd coupling to a local real-space Hamiltonian Hr
        integer,intent(in)           ::  dimH
        complex(8),intent(inout)     ::  Hr(dimH,dimH)
        integer,intent(in)           ::  TB_ext(2)
        real(8),intent(in)           ::  Jsd(:)
        type(vec_point),intent(in)   ::  mode_mag(:)

        integer                 ::  i,j
        integer                 ::  i1,i2
        integer                 ::  N_neighbours,N_cells,dim_mode
        integer                 ::  a,b
        
        complex(8),allocatable  ::  add_Jsd(:,:)
        integer                 ::  dim_mode_red

        N_neighbours = size(energy%line,1)
        N_cells = size(energy%line,2)
        dim_mode=Tb_ext(2)-Tb_ext(1)+1
        dim_mode_red=dim_mode/2
        if(size(Jsd) /= dim_mode_red) STOP "JSD has wrong size"
        allocate(add_Jsd(dim_mode,dim_mode))

        do i=1,N_cells
            add_Jsd=cmplx(0.0d0,0.0d0,8)
            do j=1,dim_mode_red
                i1=j
                i2=j+dim_mode_red
                !could be done be elegant with pauli matrices...
                add_Jsd(i1,i1)=Jsd(j)*cmplx( mode_mag(i)%w(3), 0.0d0           ,8)
                add_Jsd(i2,i1)=Jsd(j)*cmplx( mode_mag(i)%w(1), mode_mag(i)%w(2),8)
                add_Jsd(i1,i2)=Jsd(j)*cmplx( mode_mag(i)%w(1),-mode_mag(i)%w(2),8)
                add_Jsd(i2,i2)=Jsd(j)*cmplx(-mode_mag(i)%w(3), 0.0d0           ,8)
            enddo
            a=(i-1)*dim_mode+1
            b=i*dim_mode
            Hr(a:b,a:b) = Hr(a:b,a:b)+add_Jsd
        enddo 
    end subroutine 


end module
