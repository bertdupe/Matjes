module m_energy_set_real_sparse
use m_energy_commons, only : energy
use m_basic_types, only : vec_point

!
! THIS DOES NOT WORK SINCE I DIDN'T LINK A COMPLEX SPARSE EIGENVALUE SOLVER
!

private
public set_Hr,set_Jsd,Hr_eigval,Hr_eigvec



contains

    subroutine set_Hr(dimH,Tb_ext)
        use m_eigen_interface, only : eigen_set_H_e
        !extract the real space Hamiltonian Hr from the electronic part in energy
        integer,intent(in)      :: dimH
        integer,intent(in)      :: TB_ext(2)

        integer                 :: nnz
        integer                 :: N_persite

        integer                 :: N_neighbours,N_cells,dim_mode
        real(8)                 :: E
        integer                 :: i,j,k,ii,i2,k2

!sparse Hamiltonian
        complex(8),allocatable  :: val(:)
        integer,allocatable     :: rowind(:),colind(:)

        N_neighbours = size( energy%line, 1 )
        N_cells = size( energy%line, 2 )
        dim_mode=Tb_ext(2)-Tb_ext(1)+1

        !check nnz only from one site ( might not work with inhomogeneous Hamiltonians)
        N_persite=0
        do i=1,N_neighbours
            k = energy%line(i, 1) 
            N_persite=N_persite+count(energy%value(i,k)%order_op(1)%Op_loc(TB_ext(1):TB_ext(2),TB_ext(1):TB_ext(2)) /= 0.0d0)
        enddo
        nnz=N_persite*N_cells

        !set the hamiltonian in coordinate matrix storage
        allocate(rowind(nnz),colind(nnz),source=0)
        allocate(val(nnz),source=cmplx(0.0d0,0.0d0,8))
        ii=1
        do i=1, N_cells
            do j=1, N_neighbours
                k = energy%line(j, i)
                do i2=TB_ext(1),TB_ext(2)
                    do k2=TB_ext(1),TB_ext(2)
                        E=energy%value(j,k)%order_op(1)%Op_loc(i2,k2)
                        if(E/=0.0d0)then
                            rowind(ii)=dim_mode*(i-1)+i2-TB_ext(1)+1
                            colind(ii)=dim_mode*(k-1)+k2-TB_ext(1)+1
                            val(ii)=cmplx(E,0.0d0,8)
                            ii=ii+1
                        endif
                    enddo
                enddo
            enddo
        enddo

        !c++ format
        colind=colind-1
        rowind=rowind-1
        Call eigen_set_H_e(size(rowind),dimH,rowind,colind,val) 
        deallocate(rowind,colind,val)
   end subroutine 

    subroutine set_Jsd(dimH,Tb_ext,mode_mag,Jsd)
        use m_eigen_interface, only : eigen_set_H_e_jsd
        !adds the Jsd coupling to a local real-space Hamiltonian Hr
        integer,intent(in)           ::  dimH
        integer,intent(in)           ::  TB_ext(2)
        real(8),intent(in)           ::  Jsd(:)
        type(vec_point),intent(in)   ::  mode_mag(:)

        integer                 ::  N_neighbours,N_cells,dim_mode,dim_mode_red
        
        integer                 :: nnz
        complex(8),allocatable  :: val(:)
        integer,allocatable     :: rowind(:),colind(:)

        complex(8)              :: jsdmat(2,2)
        integer                 :: matind(2)

        integer                 ::  i,j
        integer                 ::  i1,i2,ii

        N_neighbours = size(energy%line,1)
        N_cells = size(energy%line,2)
        dim_mode=Tb_ext(2)-Tb_ext(1)+1
        dim_mode_red=dim_mode/2

        nnz=N_cells*dim_mode_red*4
        allocate(rowind(nnz),colind(nnz),source=0)
        allocate(val(nnz),source=cmplx(0.0d0,0.0d0,8))

        ii=1
        do i=1,N_cells
            do j=1,dim_mode_red
                matind=[j,j+dim_mode_red] + ((i-1)*dim_mode)
                jsdmat(1,1)=Jsd(j)*cmplx( mode_mag(i)%w(3), 0.0d0           ,8)
                jsdmat(2,1)=Jsd(j)*cmplx( mode_mag(i)%w(1), mode_mag(i)%w(2),8)
                jsdmat(1,2)=Jsd(j)*cmplx( mode_mag(i)%w(1),-mode_mag(i)%w(2),8)
                jsdmat(2,2)=Jsd(j)*cmplx(-mode_mag(i)%w(3), 0.0d0           ,8)
                do i1=1,2
                    do i2=1,2
                        rowind(ii)=matind(i1)
                        colind(ii)=matind(i2)
                        val(ii)=jsdmat(i1,i2)
                        ii=ii+1
                    enddo
                enddo
            enddo
        enddo 


        colind=colind-1
        rowind=rowind-1
        Call eigen_set_H_e_jsd(size(rowind),dimH,rowind,colind,val) 
        deallocate(rowind,colind,val)
        STOP "JSD SET"
    end subroutine 

    subroutine Hr_eigval(dimH,eigval)
        integer,intent(in)          ::  dimH
        real(8),intent(out)         ::  eigval(dimH)

        complex(8)                  :: H_loc(dimH,dimH)
        complex(kind=8)             :: WORK(2*dimH)
        real(kind=8)                :: RWORK(3*dimH-2)
        integer                     :: info
        
      !  if(.not.allocated(Hr)) STOP "Hr is not allocated but Hr_eigval is called"
      !  if(dimH/=size(Hr,1)) STOP "dimensions of Hr seems to be wrong evaluating the eigenvalues"
      !  H_loc=Hr
      !  call ZHEEV( 'N', 'U', dimH, H_loc, dimH, eigval, WORK, size(Work), RWORK, INFO )

    end subroutine

    subroutine Hr_eigvec(dimH,eigvec,eigval)
        integer,intent(in)          ::  dimH
        complex(8),intent(out)      ::  eigvec(dimH,dimH)
        real(8),intent(out)         ::  eigval(dimH)

        complex(kind=8)             :: WORK(2*dimH)
        real(kind=8)                :: RWORK(3*dimH-2)
        integer                     :: info

      !  if(.not.allocated(Hr)) STOP "Hr is not allocated but Hr_eigvec is called"
      !  if(dimH/=size(Hr,1)) STOP "dimensions of Hr seems to be wrong evaluating the eigenvectors"
      !  eigvec=Hr
      !  call ZHEEV( 'V', 'U', dimH, eigvec, dimH, eigval, WORK, size(Work), RWORK, INFO )

    end subroutine 
end module
