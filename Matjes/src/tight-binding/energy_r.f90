module m_energy_r
    use m_energy_commons, only : energy
    implicit none

    private
    public :: get_eigenval_r
    contains

    subroutine get_eigenval_r(dimH,tb_ext,eigval)
        integer,intent(in)  :: dimH
        integer,intent(in)  :: tb_ext(2)
        real(8),intent(out) :: eigval(dimH)

        complex(8)          :: Hr(dimH,dimH)

        Call set_Hr(dimH,Hr,tb_ext)
        Call get_eigval(dimH,Hr,eigval)
        Call write_eigval(eigval)

    end subroutine


    subroutine write_eigval(eigval)
        use m_io_files_utils, only: close_file,open_file_write
        real(8),intent(in)     :: eigval(:)
        integer                :: i,io

        io=open_file_write('eigval.dat')
        do i=1,size(eigval)
           write(io,'(E16.8)') eigval(i)
        enddo
        call close_file('eigval.dat',io)

    end subroutine 

    subroutine set_Hr(dimH,Hr,Tb_ext)
        integer,intent(in)      ::  dimH
        complex(8)              ::  Hr(dimH,dimH)
        integer,intent(in)      ::  TB_ext(2)

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

    end subroutine set_Hr

    subroutine get_eigval(dimH,Hr,eigval)
        integer,intent(in)          ::  dimH
        complex(8),intent(inout)    ::  Hr(dimH)
        real(8),intent(out)         ::  eigval(dimH)

        complex(kind=8)             :: WORK(2*dimH)
        real(kind=8)                :: RWORK(3*dimH-2)
        integer                     :: info

        call ZHEEV( 'N', 'U', dimH, Hr, dimH, eigval, WORK, size(Work), RWORK, INFO )

    end subroutine get_eigval

end module m_energy_r
