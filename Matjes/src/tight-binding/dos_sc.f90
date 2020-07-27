module m_dos_sc
implicit none
private
public write_dos_sc
contains


subroutine write_dos_sc(eigval,eigvec,fname)
    use m_io_files_utils, only: open_file_read,close_file
    use m_io_utils, only: get_parameter
    real(8),intent(in)      ::  eigval(:)
    complex(8),intent(in)   ::  eigvec(:,:)
    character(len=*)        ::  fname

    integer                     :: io_input
    logical                     :: do_dos
    real(8)                     :: E_ext(2),sigma,dE

    io_input=open_file_read('input')
    do_dos=.False.
    call get_parameter(io_input,'input','do_dos_r',do_dos)
    if(do_dos)then
        E_ext=[-1.0d0,1.0d0]
        dE=0.01d0
        sigma=0.01d0
        call get_parameter(io_input,'input','dos_sigma',sigma)
        call get_parameter(io_input,'input','dos_E_ext',2,E_ext)
        call get_parameter(io_input,'input','dos_dE',dE)
        call close_file('input',io_input)
        Call calc_dos_sc(eigval,eigvec,E_ext,dE,sigma,fname)
    endif
end subroutine


subroutine calc_dos_sc(eigval,eigvec,E_ext,dE,sigma,fname)
    !subroutine to calculate the density of states
    use m_io_files_utils, only: close_file,open_file_write
    real(8),intent(in)     ::  eigval(:)
    complex(8),intent(in)  ::  eigvec(:,:)
    real(8),intent(in)  ::  E_ext(2),sigma,dE
    character(len=*)    ::  fname

    integer             ::  NE,iE
    real(8),allocatable ::  dos(:),Eval(:)
    integer             ::  dimH

    integer                     :: i,io

    dimH=size(eigval)

    Ne=int((E_ext(2)-E_ext(1))/dE)+1
    allocate(dos(Ne),Eval(Ne),source=0.0d0)
    do iE=1,Ne
        Eval(iE)=(iE-1)*dE+E_ext(1)
    enddo

    do i=dimH/2+1,dimH
        Call add_dos(eigval(i),eigvec(:,i),Ne,Eval,E_ext,dE,sigma,dos)
    enddo
    dos=dos/real(size(eigval,1)/2)

    io=open_file_write(fname)
    do i=1,size(dos)
       write(io,'(2E16.8)') Eval(i),dos(i)
    enddo
    call close_file(fname,io)
end subroutine

subroutine add_dos(val,eigvec,Ne,Eval,E_ext,dE,sigma,dos)
    real(8),intent(in)      ::  val,sigma
    integer,intent(in)      ::  Ne
    real(8),intent(in)      ::  Eval(Ne),E_ext(2),dE
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
    i_min=int(((val-sigma*dist_inc)-E_ext(1))/dE)+1
    i_max=int(((val+sigma*dist_inc)-E_ext(1))/dE)
    pref=dot_product(eigvec(1:dimH/2),eigvec(1:dimH/2))
    Call add_gauss(val,pref,Eval(i_min:i_max),dos(i_min:i_max),sigma)

    
    !v-part of BdG
    i_min=int(((-val-sigma*dist_inc)-E_ext(1))/dE)+1
    i_max=int(((-val+sigma*dist_inc)-E_ext(1))/dE)
    pref=dot_product(eigvec(dimH/2+1:dimH),eigvec(dimH/2+1:dimH))
    Call add_gauss(-val,pref,Eval(i_min:i_max),dos(i_min:i_max),sigma)

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
