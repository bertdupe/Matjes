module m_dos_sc
use m_TB_types
implicit none
private
public calc_dos_sc
contains


subroutine calc_dos_sc(eigval,eigvec,io_dos,fname)
    !subroutine to calculate the density of states
    use m_io_files_utils, only: close_file,open_file_write
    real(8),intent(in)                    ::  eigval(:)
    complex(8),intent(in)                 ::  eigvec(:,:)
    type(parameters_TB_IO_DOS),intent(in) ::  io_dos
    character(len=*)    ::  fname

    integer             ::  NE,iE
    real(8),allocatable ::  dos(:),Eval(:)
    integer             ::  dimH

    integer             ::  i,io

    dimH=size(eigval)

    Ne=int((io_dos%E_ext(2)-io_dos%E_ext(1))/io_dos%dE)+1
    allocate(dos(Ne),Eval(Ne),source=0.0d0)
    do iE=1,Ne
        Eval(iE)=(iE-1)*io_dos%dE+io_dos%E_ext(1)
    enddo

    do i=dimH/2+1,dimH
        Call add_dos(eigval(i),eigvec(:,i),Ne,Eval,io_dos,dos)
    enddo
    dos=dos/real(size(eigval,1)/2)

    io=open_file_write(fname)
    do i=1,size(dos)
       write(io,'(2E16.8)') Eval(i),dos(i)
    enddo
    call close_file(fname,io)
end subroutine

subroutine add_dos(val,eigvec,Ne,Eval,io_dos,dos)
    real(8),intent(in)      ::  val
    integer,intent(in)      ::  Ne
    real(8),intent(in)      ::  Eval(Ne)
    type(parameters_TB_IO_DOS),intent(in) :: io_dos
    complex(8),intent(in)   ::  eigvec(:)
    real(8),intent(inout)   ::  dos(Ne)
    
    real(8)                 ::  pref
    real(8),parameter   ::  dist_inc=5.0d0
    integer             ::  i_min,i_max
    integer             ::  i
    integer             ::  dimH

    !need to add correct prefactors from eigvec
    pref=1.0d0
    dimH=size(eigvec)

    !u-part of BdG
    i_min=int(((val-io_dos%sigma*dist_inc)-io_dos%E_ext(1))/io_dos%dE)+1
    i_max=int(((val+io_dos%sigma*dist_inc)-io_dos%E_ext(1))/io_dos%dE)
    pref=dot_product(eigvec(1:dimH/2),eigvec(1:dimH/2))
    Call add_gauss(val,pref,Eval(i_min:i_max),dos(i_min:i_max),io_dos%sigma)

    
    !v-part of BdG
    i_min=int(((-val-io_dos%sigma*dist_inc)-io_dos%E_ext(1))/io_dos%dE)+1
    i_max=int(((-val+io_dos%sigma*dist_inc)-io_dos%E_ext(1))/io_dos%dE)
    pref=dot_product(eigvec(dimH/2+1:dimH),eigvec(dimH/2+1:dimH))
    Call add_gauss(-val,pref,Eval(i_min:i_max),dos(i_min:i_max),io_dos%sigma)

end subroutine

subroutine add_gauss(val,pref,E,dos,sigma)
    !add to dos:
    !gauss distribution from single energy point val into the spacing supplied by E with std. sigma and an additional prefactor pref 
    real(8),intent(in)     ::  val,pref,sigma
    real(8),intent(in)     ::  E(:)
    real(8),intent(inout)  ::  dos(:)
    real(8)                ::  tmp(size(dos))
    tmp=(val-E)**2
    tmp=-tmp*0.5d0/sigma/sigma
    tmp=exp(tmp)
    tmp=tmp/sqrt(2.0d0*3.14159265359d0)/sigma
    dos=dos+pref*tmp
end subroutine

end module
