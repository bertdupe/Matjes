module m_H_type_manual
!Hamiltonian type specifications to use Bertrands manual matrix multiplication
use m_H_type

type,extends(t_H) :: t_H_manual
    private
    integer, allocatable :: line(:,:)
	real(8), allocatable, dimension(:,:) :: all_E
contains
	procedure :: eval_single=>eval_single
	procedure :: eval_all=>eval_all
	procedure :: set_H=>set_H
end type
private
public t_H,t_H_manual
contains 


subroutine set_H(this,energy_in,lat)
    use m_derived_types, only: operator_real_order_N
	use m_derived_types, only: lattice
	!need to get rid of dim_mode input
	class(t_H_manual),intent(inout)  :: this
	type(operator_real_order_N)     :: energy_in
	type(lattice),intent(in)    	:: lat

	integer	:: N,i

	allocate(this%line,source=energy_in%line)
	N=size(this%line(:,1))

	allocate(this%all_E(lat%dim_mode,lat%dim_mode*N),source=0.0d0)
	do i=1,N
	   this%all_E(:,(i-1)*lat%dim_mode+1:i*lat%dim_mode)=transpose(energy_in%value(i,1)%order_op(1)%Op_loc)
	enddo
end subroutine 

subroutine eval_single(this,E,i_m,lat)
	use m_derived_types, only: lattice
	! input
	class(t_H_manual),intent(in)    :: this
	type(lattice), intent(in) 		:: lat
	integer, intent(in) 			:: i_m
	! output
	real(kind=8), intent(out) 		:: E
	! internal
	real(8) :: all_vectors(lat%dim_mode*size(this%line(:,1)))
	integer :: i,N,j
	
	N=size(this%line(:,i_m))
	do i=1,N
	   j=this%line(i,i_m)
	   all_vectors((i-1)*lat%dim_mode+1:i*lat%dim_mode)=lat%ordpar%all_l_modes(j)%w
	enddo
	E=dot_product(lat%ordpar%all_l_modes(i_m)%w , matmul( this%all_E , all_vectors ) )
end subroutine 


subroutine eval_all(this,E,lat)
    use m_derived_types, only: lattice
	class(t_H_manual),intent(in)  :: this
    type(lattice), intent(in) :: lat
    real(8), intent(out) :: E
    ! internal
    real(8) :: E_i
    integer :: i,N,j
    
    E=0.0d0
    do i=1,product(lat%dim_lat)
        call this%eval_single(E_i,i,lat)
        E=E+E_i !might be dangerous for very large arrays (accuracy)
    enddo
    
end subroutine 

end module
