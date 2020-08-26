module m_save_state_r
implicit none
private
public :: TB_write_states_r,TB_read_states_r
integer,parameter     :: save_version=1
integer,parameter     :: io_unit=700    !arbitrary value hopefully not openen 
character(len=*),parameter      :: filen='TB_solution_r.dat'
contains

subroutine TB_write_states_r(eigval,eigvec)
    !subroutine writing previously calculated eigenvalues and eigenvectors to file
    real(8),intent(in),allocatable  :: eigval(:)
    complex(8),intent(in),allocatable  :: eigvec(:,:)
    integer             ::  dimH,numE
    logical             ::  leigval,leigvec

    logical             :: is_open
    integer             :: data_size !estimate data_size

    dimH=0; numE=0;data_size=0
    leigval=allocated(eigval)
    leigvec=allocated(eigvec)
    if(leigval)then
        numE=size(eigval)
        data_size=data_size+STORAGE_SIZE(eigval)*size(eigval)
    endif
    if(leigvec)then
        data_size=data_size+STORAGE_SIZE(eigvec)*size(eigvec)
        dimH=size(eigvec,1)
    endif
    write(*,*) 'Estimating size for '//filen//'=',data_size/8/1024/1024,'Mb'
    if(.not.(leigval.or.leigvec)) return

    inquire(io_unit,OPENED=is_open)
    if(is_open)then
        write(*,*) 'Could not write '//filen//' as io-unit is already opened'
        write(*,*) 'Change the code'
        return
    endif
    open(io_unit,file=filen,form='unformatted')
    write(io_unit) save_version
    write(io_unit) leigval,leigvec
    write(io_unit) dimH,numE
    if(leigval) write(io_unit) eigval
    if(leigvec) write(io_unit) eigvec
    close(io_unit)
end subroutine


subroutine TB_read_states_r(eigval,eigvec,success)
    !subroutine reading previously calculated eigenvalues and eigenvectors from file
    real(8),intent(out),allocatable     :: eigval(:)
    complex(8),intent(out),allocatable  :: eigvec(:,:)
    logical,intent(out) :: success
    !additional io-values
    integer             :: dimH,numE
    integer             :: version
    logical             :: leigval,leigvec
    !internal variables
    logical             :: is_open

    inquire(file=filen,exist=success)
    if(success)then
        inquire(io_unit,OPENED=is_open)
        if(is_open)then
            write(*,*) 'Could not read '//filen//' as io-unit is already opened'
            write(*,*) 'Change the code'
            return
        endif
        open(io_unit,file=filen,form='unformatted')
        read(io_unit) version
        if(version==1)then
            read(io_unit) leigval,leigvec
            read(io_unit) dimH,numE
            if(leigval)then
                allocate(eigval(numE),source=0.0d0)
                read(io_unit) eigval
            endif
            if(leigvec)then
                allocate(eigvec(dimH,numE),source=cmplx(0.0d0,0.0d0,8))
                read(io_unit) eigvec
            endif
        else
            write(*,*) "Unsupported version of "//filen//" input"
            write(*,*) "version",version
            success=.false.
        endif
        close(io_unit)
    else
        write(*,*) 'Did not find file '//filen//' to read eigenvalues and eigenvectors'
    endif
end subroutine

end module
