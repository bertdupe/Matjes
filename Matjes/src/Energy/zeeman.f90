module m_zeeman
use m_symmetry_operators
use m_lattice, only : my_order_parameters
use m_Hamiltonian_variables, only : coeff_ham_inter_spec
type(coeff_ham_inter_spec), target, public, protected :: zeeman

private
public :: get_ham_zeeman, get_zeeman_H

contains

subroutine get_ham_zeeman(fname,dim_ham,Ms)
use m_io_files_utils
use m_vector, only : norm
use m_io_utils
use m_convert
use m_constants, only : mu_B,mu_0
implicit none
real(kind=8), intent(in) :: Ms
integer, intent(in) :: dim_ham
character(len=*), intent(in) ::fname
! internal
integer :: io_param
real(kind=8) :: h_ext(3)
! anisotropy
integer :: x_start,x_end
! electric field
integer :: y_start,y_end
integer :: i
character(len=50) :: form

h_ext=0.0d0
zeeman%name='zeeman'
zeeman%order=2

io_param=open_file_read(fname)
call get_parameter(io_param,fname,'c_zeeman',zeeman%c_ham)
! count the magnetic field if present
call get_parameter(io_param,fname,'H_ext',3,h_ext)

call close_file(fname,io_param)

if (norm(H_ext).ge.1.0d-8) then
   zeeman%i_exist=.true.
  else
   return
endif

allocate(zeeman%ham(1))
allocate(zeeman%ham(1)%H(dim_ham,dim_ham))

zeeman%ham(1)%H=0.0d0


call get_borders('magnetic',x_start,x_end,'Bfield',y_start,y_end,my_order_parameters)

! get the zeeman Hamiltonian (get the symmetric tensor)
call get_symmetric_Op(zeeman%ham(1)%H,mu_0*mu_B*Ms*zeeman%c_ham,x_start,y_start,x_end,y_end)
! make it asymmetric
zeeman%ham(1)%H(:,y_start:y_end)=0.0d0

form=convert('(',dim_ham,'(f12.8,2x))')
write(6,'(a)') ''
write(6,'(a)') 'Zeeman tensor'
do i=1,dim_ham
    write(6,form) zeeman%ham(1)%H(:,i)
enddo
write(6,'(a)') ''

end subroutine get_ham_zeeman


subroutine get_zeeman_H(Ham,lat)
    !get zeeman in t_H Hamiltonian format
    !so far zeeman has to be set before
    !unsymmetric in that it only includes the B M basis and not the revers ( similar to previous implementation)
    use m_Htype_gen
    use m_derived_types
    use m_setH_util,only: get_coo

    class(t_H),intent(inout)    :: Ham
    type(lattice),intent(in)    :: lat

    !describe local Hamiltonian
    real(8),allocatable         :: Htmp(:,:)
    !input for setting t_H
    real(8),allocatable         :: val_tmp(:)
    integer,allocatable         :: ind_tmp(:,:)
    integer,allocatable         :: line(:,:)

    integer :: i
    integer :: x_start,x_end
    integer :: y_start,y_end
    integer :: Ncell


    if(zeeman%i_exist)then
        Ncell=product(lat%dim_lat)
        allocate(line(1,Ncell),source=0) !1, since always onsite
        call get_borders('Bfield',x_start,x_end,'magnetic',y_start,y_end,my_order_parameters)
        allocate(Htmp,source=zeeman%ham(1)%H(x_start:x_end,y_start:y_end))
    
        !get local Hamiltonian in coo format
        Call get_coo(Htmp,val_tmp,ind_tmp)
       
        !simple on-site term
        do i=1,Ncell
            line(1,i)=i
        enddo
    
        Call Ham%init_1(line,val_tmp,ind_tmp,[3,1],lat)
    endif
end subroutine

end module m_zeeman
