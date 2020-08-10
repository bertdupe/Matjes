module m_energy_set_real
use m_energy_commons, only : energy
use m_basic_types, only : vec_point

private
public set_Hr,Hr_eigval,Hr_eigvec,get_Hr

!large electronic Hamiltonian
complex(8),allocatable  ::  Hr(:,:)

!setup is now a bit stupid in that Hr has to be saved and copied for the lapack routines, which doubles their memory at some points...

contains

    subroutine get_Hr(dimH,Hr_out)
        integer,intent(in)          ::  dimH
        complex(8),intent(out)      ::  Hr_out(dimH,dimH)
       
        if(.not.allocated(Hr)) STOP "Hr is not allocated but get_Hr is called"
        if(dimH/=size(Hr,1)) STOP "dimensions of Hr seems to be wrong evaluating the eigenvalues"
        Hr_out=Hr
    end subroutine 

    subroutine set_Hr(dimH,Tb_ext,mode_mag)
        use m_tb_params, only : TB_params
        integer,intent(in)           ::  dimH
        integer,intent(in)           ::  TB_ext(2)
        type(vec_point),intent(in)   ::  mode_mag(:)

        if(.not. allocated(Hr))then
           allocate(Hr(dimH,dimH))
        endif
        if(size(Hr,1)/=dimH.or.size(Hr,2)/=dimH) STOP "Hr has wrong size"  !could easily reallocate, but this should never happen, I guess
        Hr=cmplx(0.0d0,0.0d0, kind=8)
        Call set_Hr_ee(dimH,Tb_ext)
        if(any(TB_params%io_H%Jsd /= 0.0d0))then
            Call set_Jsd(dimH,tb_ext,mode_mag,TB_params%io_H%Jsd)
        endif
    end subroutine 

    subroutine set_Hr_ee(dimH,Tb_ext)
        !extract the real space Hamiltonian Hr from the electronic part in energy
        integer,intent(in)           ::  dimH
        integer,intent(in)           ::  TB_ext(2)

        integer                 ::  N_neighbours,N_cells,dim_mode
        integer                 ::  i,j,k

        N_neighbours = size( energy%line, 1 )
        N_cells = size( energy%line, 2 )
        dim_mode=Tb_ext(2)-Tb_ext(1)+1
        if(.not. allocated(Hr))then
           allocate(Hr(dimH,dimH))
        endif
        if(size(Hr,1)/=dimH.or.size(Hr,2)/=dimH) STOP "Hr has wrong size"  !could easily reallocate, but this should never happen, I guess
        Hr=cmplx(0.0d0,0.0d0, kind=8)
        do i=1, N_cells
            do j=1, N_neighbours
                k = energy%line(j, i) 
                Hr((i-1)*dim_mode+1:i*dim_mode, (k-1)*dim_mode+1:k*dim_mode) = energy%value(j,k)%order_op(1)%Op_loc(TB_ext(1):TB_ext(2),TB_ext(1):TB_ext(2))
            enddo
        enddo
    end subroutine 

    subroutine set_Jsd(dimH,Tb_ext,mode_mag,Jsd)
        !adds the Jsd coupling to a local real-space Hamiltonian Hr
        integer,intent(in)           ::  dimH
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

        if(size(Hr,1)/=dimH.or.size(Hr,2)/=dimH) STOP "Hr has wrong size" 
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

    subroutine Hr_eigval(dimH,eigval)
        integer,intent(in)          ::  dimH
        real(8),intent(out)         ::  eigval(dimH)

        complex(8)                  :: H_loc(dimH,dimH)
        complex(kind=8)             :: WORK(2*dimH)
        real(kind=8)                :: RWORK(3*dimH-2)
        integer                     :: info
        
        if(.not.allocated(Hr)) STOP "Hr is not allocated but Hr_eigval is called"
        if(dimH/=size(Hr,1)) STOP "dimensions of Hr seems to be wrong evaluating the eigenvalues"
        H_loc=Hr
        call ZHEEV( 'N', 'U', dimH, H_loc, dimH, eigval, WORK, size(Work), RWORK, INFO )

    end subroutine

    subroutine Hr_eigvec(dimH,eigvec,eigval)
        integer,intent(in)          ::  dimH
        complex(8),intent(out)      ::  eigvec(dimH,dimH)
        real(8),intent(out)         ::  eigval(dimH)

        complex(kind=8)             :: WORK(2*dimH)
        real(kind=8)                :: RWORK(3*dimH-2)
        integer                     :: info

        if(.not.allocated(Hr)) STOP "Hr is not allocated but Hr_eigvec is called"
        if(dimH/=size(Hr,1)) STOP "dimensions of Hr seems to be wrong evaluating the eigenvectors"
        eigvec=Hr
        call ZHEEV( 'V', 'U', dimH, eigvec, dimH, eigval, WORK, size(Work), RWORK, INFO )

    end subroutine 

end module
