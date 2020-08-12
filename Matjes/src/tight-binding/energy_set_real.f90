module m_energy_set_real
use m_energy_commons, only : energy
use m_basic_types, only : vec_point
use m_tb_types
implicit none
private
public set_Hr_dense_nc


contains

    subroutine set_Hr_dense_nc(h_par,mode_mag,Hr)
        use m_tb_params, only : TB_params
        type(parameters_TB_Hsolve),intent(in)     ::  h_par
        type(vec_point),intent(in)   ::  mode_mag(:)
        complex(8),allocatable,intent(inout)    ::  Hr(:,:)

        if(.not. allocated(Hr))then
           allocate(Hr(h_par%dimH,h_par%dimH))
        endif
        if(size(Hr,1)/=h_par%dimH.or.size(Hr,2)/=h_par%dimH) STOP "Hr has wrong size"  !could easily reallocate, but this should never happen, I guess
        Hr=cmplx(0.0d0,0.0d0, kind=8)
        Call set_Hr_ee(h_par,Hr)
        if(any(TB_params%io_H%Jsd /= 0.0d0))then
            Call set_Jsd(h_par,mode_mag,TB_params%io_H%Jsd,Hr)
        endif
    end subroutine 

    subroutine set_Hr_ee(h_par,Hr)
        !extract the real space Hamiltonian Hr from the electronic part in energy
        type(parameters_TB_Hsolve),intent(in)     ::  h_par
        complex(8),allocatable,intent(inout)    ::  Hr(:,:)

        integer                 ::  N_neighbours,N_cells,dim_mode
        integer                 ::  i,j,k

        N_neighbours = size( energy%line, 1 )
        N_cells = size( energy%line, 2 )
        dim_mode=h_par%pos_ext(2)-h_par%pos_ext(1)+1
        if(.not. allocated(Hr))then
           allocate(Hr(h_par%dimH,h_par%dimH))
        endif
        if(size(Hr,1)/=h_par%dimH.or.size(Hr,2)/=h_par%dimH) STOP "Hr has wrong size"  !could easily reallocate, but this should never happen, I guess
        Hr=cmplx(0.0d0,0.0d0, kind=8)
        do i=1, N_cells
            do j=1, N_neighbours
                k = energy%line(j, i) 
                Hr((i-1)*dim_mode+1:i*dim_mode, (k-1)*dim_mode+1:k*dim_mode) = energy%value(j,k)%order_op(1)%Op_loc(h_par%pos_ext(1):h_par%pos_ext(2),h_par%pos_ext(1):h_par%pos_ext(2))
            enddo
        enddo
    end subroutine 

    subroutine set_Jsd(h_par,mode_mag,Jsd,Hr)
        !adds the Jsd coupling to a local real-space Hamiltonian Hr
        type(parameters_TB_Hsolve),intent(in)        ::  h_par
        real(8),intent(in)                          ::  Jsd(:)
        type(vec_point),intent(in)                  ::  mode_mag(:)
        complex(8),allocatable,intent(inout)        ::  Hr(:,:)

        integer                 ::  i,j
        integer                 ::  i1,i2
        integer                 ::  N_neighbours,N_cells,dim_mode
        integer                 ::  a,b
        
        complex(8),allocatable  ::  add_Jsd(:,:)
        integer                 ::  dim_mode_red

        N_neighbours = size(energy%line,1)
        N_cells = size(energy%line,2)
        dim_mode=h_par%pos_ext(2)-h_par%pos_ext(1)+1
        dim_mode_red=dim_mode/2

        if(size(Hr,1)/=h_par%dimH.or.size(Hr,2)/=h_par%dimH) STOP "Hr has wrong size" 
        if(size(Jsd) /= dim_mode_red) STOP "JSD has wrong size"

        allocate(add_Jsd(dim_mode,dim_mode))
        do i=1,N_cells
            add_Jsd=cmplx(0.0d0,0.0d0,8)
            do j=1,dim_mode_red
                i1=j*2-1
                i2=j*2
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
