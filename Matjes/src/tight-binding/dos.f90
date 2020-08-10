module m_dos
!this module contains a badly designed dos calculation using a gauss smearing, but so far it was sufficiently fast
use m_TB_types, only: parameters_TB_IO_DOS
implicit none
private
public calc_dos
contains

subroutine calc_dos(eigval,io_dos,fname)
    !subroutine to calculate the density of states
    !eigval input has to be sorted
    use m_io_files_utils, only: close_file,open_file_write
    real(8),intent(in)  ::  eigval(:)
    type(parameters_TB_IO_DOS) :: io_dos
    character(len=*)    ::  fname


    integer             ::  NE,iE
    real(8),allocatable ::  dos(:),Eval(:)

    integer                     :: i,io

    Ne=int((io_dos%E_ext(2)-io_dos%E_ext(1))/io_dos%dE)+1
    allocate(dos(Ne),Eval(Ne),source=0.0d0)
    do iE=1,Ne
        Eval(iE)=(iE-1)*io_dos%dE+io_dos%E_ext(1)
    enddo

    do iE=1,Ne
        Call get_dos(eigval,Eval(iE),dos(iE),io_dos%sigma)
    enddo
    dos=dos/real(size(eigval,1))

    io=open_file_write(fname)
    do i=1,size(dos)
       write(io,'(2E16.8)') Eval(i),dos(i)
    enddo
    call close_file(fname,io)
end subroutine


subroutine get_dos(val,E,res,sigma)
    real(8),intent(in)  ::  val(:),E,sigma
    real(8),intent(out) ::  res

    real(8)             ::  tmp(size(val,1))

    Call gauss_dist(val,E,sigma,tmp)
    res=sum(tmp)

end subroutine

subroutine gauss_dist(val,mu,sigma,res)
    !gauss distribution using 1/\sqrt{2*\pi*sigma^2}*e^{-(val-mu)^2/(2*sigma^2)} for all values in val respective to one given mu and sigma
    real(8),intent(in)  ::  val(:),mu,sigma
    real(8),intent(out) ::  res(size(val,1))
    integer             ::  i

    !only consider vals which are within [mu-dist_inc*sigma,mu+dist_inc*sigma], because all 
    !other results in small res  anyways and a too small exponent becomes numerically problematic
    real(8),parameter   ::  dist_inc=5.0d0

    integer             ::  i_min,i_max,Ne

    
    Ne=size(val,1)
    res=0.0d0
    if(val(Ne) <=mu-dist_inc*sigma.or.val(1) >=mu+dist_inc*sigma) return
    i_min=1
    do i=1,Ne
        if(val(i) >=mu-dist_inc*sigma)then
            i_min=i
            exit
        endif
    enddo
    i_max=i_min
    do i=Ne,i_min,-1
        if(val(i) <=mu+dist_inc*sigma)then
            i_max=i
            exit
        endif
    enddo
    if(i_min.eq.i_max) return
    res(i_min:i_max)=(val(i_min:i_max)-mu)*(val(i_min:i_max)-mu)
    res(i_min:i_max)=-res(i_min:i_max)*0.5d0/sigma/sigma
    res(i_min:i_max)=exp(res(i_min:i_max))
    res=res/sqrt(2.0d0*3.14159265359d0)/sigma

end subroutine
end module
