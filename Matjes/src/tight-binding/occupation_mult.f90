module m_occupation_mult
use m_tb_types
use m_distribution, only: int_distrib,fermi_distrib,dE_fermi_distrib
use m_occupation, only: calc_occupation
public occupation_mult
private
integer,parameter           ::  io_unit=701
character(len=*),parameter  ::  subdir='data_state_dist/'

contains
subroutine occupation_mult(h_par,io_par,eigval,eigvec)
    type(parameters_TB_Hsolve),intent(in)       :: h_par
    type(parameters_TB_IO_OCC_MULT),intent(in)  :: io_par
    real(8),intent(in)                          :: eigval(:)
    complex(8),intent(in)                       :: eigvec(:,:)
    
    !input parameters
    real(8)                     ::  E_ext(2),dE
    real(8)                     ::  kt

    !local parameters
    integer                     ::  NE
    real(8),allocatable         ::  E(:)
    integer                     ::  i
    character(len=3)            ::  i_char
    procedure(int_distrib),pointer  :: dist_ptr => null()

    !check input parameters
    if(io_par%dE<0.or.io_par%dE>io_par%E_ext(2)-io_par%E_ext(1))then
        write(*,'(/A)') 'Energy io parametes for mult_occ are inconsistent'
        write(*,*) 'Ef_ext=',io_par%E_ext,'dE=',io_par%dE
        write(*,'(A/)') 'Aborting occupation_mult'
        return
    endif
    if(io_par%kt<=0.0d0)then
        write(*,'(/A)') 'kt parametes for mult_occ is smaller 0 (default)'
        write(*,'(A)') 'this parameter has to be set'
        write(*,'(A/)') 'Aborting occupation_mult'
        return
    endif
    dE=io_par%dE
    E_ext=io_par%E_ext
    kt=io_par%kt
!get fermi energies to consider
    NE=(E_ext(2)-E_ext(1))/dE+1
    !more would require larger i_char, but this is probably 
    !a reasonable check anyways to prevent IO-madness
    if(NE>=999)then   
        write(*,'(/A)') 'Too many energies chosen to plot occupation'
        write(*,*) 'Wrong E_ext or dE?'
        write(*,*) 'E_ext=',E_ext,'dE=',dE,'NE=',NE
        write(*,'(A/)') 'Aborting occupation_mult'
        return
    endif
    allocate(E(NE),source=0.0d0)
    do i=1,NE
        E(i)=E_ext(1)+(i-1)*dE
    enddo

!prepare and check the folder for the output
    call system('mkdir -p '//subdir)
    open(io_unit,file=subdir//'/test',iostat=i)
    if(i/=i)then
       write(*,*) 'failed to access subdirectory'//subdir
       write(*,*) 'Cannot write state_occupations'
       return
    endif
    close(io_unit)
    call system('rm  '//subdir//'/test')
!Calculate all
    do i=1,NE
        write(i_char,'(I0.3)') i
        dist_ptr=>fermi_distrib
        Call calc_occupation(h_par,eigvec,eigval,E(i),kt,subdir//'occ_fermi_'//i_char//'.dat',dist_ptr)
        dist_ptr=>dE_fermi_distrib
        Call calc_occupation(h_par,eigvec,eigval,E(i),kt,subdir//'occ_dfermi_'//i_char//'.dat',dist_ptr)
    enddo
 end subroutine
 end module
