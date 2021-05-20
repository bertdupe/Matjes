module m_H_eigen_mem
#ifdef CPP_EIGEN
!Hamiltonian type specifications extending the Eigen version to also save the transpose 
!in order to make some operations faster
use m_derived_types, only: lattice, number_different_order_parameters
use m_eigen_H_interface
use m_H_coo_based
use m_H_eigen, only: t_h_eigen
use m_mode_construction_rankN_sparse_col

type,extends(t_H_eigen) :: t_H_eigen_mem
    type(C_PTR)     ::  H_T=c_null_ptr  !transpose Hamiltonian in eigen
contains
    !necessary t_H routines
    procedure :: set_from_Hcoo

    procedure :: add_child 
    procedure :: destroy_child    
    procedure :: copy_child 

    procedure :: optimize
    procedure :: mult_l
    procedure :: mult_l_cont
    procedure :: mult_l_disc
!    procedure :: mult_r_single,mult_l_single   !can probably be made more efficient, but I am not sure if it is still necessary
    procedure :: mult_r_ind 
    procedure :: mult_r_disc_disc
    procedure :: get_ind_mult_l

    !MPI
    procedure :: send
    procedure :: recv
    procedure :: distribute
    procedure :: bcast
end type

interface t_H_mkl_csr
    procedure :: dummy_constructor
end interface 

private
public t_H,t_H_eigen_mem
contains 

type(t_H_eigen_mem) function dummy_constructor()
    !might want some initialization for H and descr, but should work without
    !continue 
end function 

subroutine mult_r_disc_disc(this,ind_r,vec_r,ind_l,vec_out)
    class(t_H_eigen_mem),intent(in) :: this
    integer,intent(in)              :: ind_r(:)
    real(8),intent(in)              :: vec_r(:)
    integer,intent(in)              :: ind_l(:)
    real(8),intent(inout)           :: vec_out(:)

    Call eigen_mult_l_disc_disc(this%H_T,size(ind_r),ind_r,vec_r,size(ind_l),ind_l,vec_out)
end subroutine

subroutine get_ind_mult_l(this,ind_in,N_out,ind_out)
    !get the indicies in ind_out(:N_out) of the right vector which are necessary to 
    !get all components that 
    class(t_H_eigen_mem),intent(in) :: this
    integer,intent(in)              :: ind_in(:)
    integer,intent(out)             :: N_out
    integer,intent(inout)           :: ind_out(:)

    N_out=size(ind_out)
    Call eigen_H_get_ind_mult_r(this%H_T,size(ind_in),ind_in,N_out,ind_out)
    ind_out(1:N_out)=ind_out(1:N_out)+1
end subroutine 

subroutine mult_l_cont(this,bnd,vec,res)
    !multiply to the right with a continuous section of the right vector
    class(t_H_eigen_mem),intent(in) :: this
    integer,intent(in)              :: bnd(2)
    real(8),intent(in)              :: vec(bnd(2)-bnd(1)+1)
    real(8),intent(inout)           :: res(:)   !result matrix-vector product

    Call eigen_H_mult_mat_vec_cont(this%H_T,bnd(1),bnd(2),vec,res)
end subroutine 

subroutine mult_l_disc(this,N,ind,vec,res)
    !multiply to the right with a discontinuous section of the right vector
    class(t_H_eigen_mem),intent(in) :: this
    integer,intent(in)              :: N
    integer,intent(in)              :: ind(N)
    real(8),intent(in)              :: vec(N)
    real(8),intent(inout)           :: res(:)   !result matrix-vector product

    Call eigen_H_mult_mat_vec_disc(this%H_T,N,ind,vec,res)
end subroutine 

subroutine mult_r_ind(this,lat,N,ind_out,vec_out)
    class(t_H_eigen_mem),intent(in) :: this
    type(lattice),intent(in)        :: lat
    integer,intent(in)              :: N
    integer,intent(in)              :: ind_out(N)
    real(8),intent(out)             :: vec_out(N)

    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)

    Call this%mode_r%get_mode(lat,modes,vec)
    Call eigen_H_mult_l_ind(this%H_T,modes,N,ind_out,vec_out)
end subroutine

subroutine mult_l(this,lat,res)
    use m_derived_types, only: lattice
    class(t_H_eigen_mem),intent(in) :: this
    type(lattice), intent(in)       :: lat
    real(8), intent(inout)          :: res(:)
    ! internal
    real(8),pointer            :: modes(:)
    real(8),allocatable,target :: vec(:)

    Call this%mode_l%get_mode(lat,modes,vec)
    Call eigen_H_mult_mat_vec(this%H_T,modes,res)
    if(allocated(vec)) deallocate(vec)
end subroutine 

subroutine optimize(this)
    class(t_H_eigen_mem),intent(inout)   :: this

    Call this%t_H_eigen%optimize()
    !THERE MIGHT BE SOMETHING TO DO HERE
    continue 
end subroutine

subroutine copy_child(this,Hout)
    class(t_H_eigen_mem),intent(in)     :: this
    class(t_H_base),intent(inout)        :: Hout
    
    Call this%t_H_eigen%copy_child(Hout)
    select type(Hout)
    class is(t_H_eigen_mem)
        Call eigen_H_copy(this%H_T,Hout%H_T) 
    class default
        STOP "Cannot copy t_H_eigen_mem type with Hamiltonian that is not a class of t_H_eigen_mem"
    end select
end subroutine

subroutine add_child(this,H_in)
    class(t_H_eigen_mem),intent(inout)    :: this
    class(t_H_base),intent(in)             :: H_in

    Call this%t_H_eigen%add_child(H_in)

    select type(H_in)
    class is(t_H_eigen_mem)
        Call eigen_H_add(this%H_T,H_in%H_T)
    class default
        STOP "Cannot add t_H_eigen_mem type with Hamiltonian that is not a class of t_H_eigen_mem"
    end select
end subroutine 

subroutine destroy_child(this)
    class(t_H_eigen_mem),intent(inout)    :: this

    Call this%t_H_eigen%destroy_child()
    if(this%is_set()) Call eigen_H_destroy(this%H_T)
end subroutine

subroutine set_from_Hcoo(this,H_coo)
    class(t_H_eigen_mem),intent(inout)  :: this
    type(t_H_coo),intent(inout)     :: H_coo

    Call this%t_H_eigen%set_from_Hcoo(H_coo)
    Call eigen_get_transpose(this%H,this%H_T)
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!            MPI ROUTINES           !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine send(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_eigen_mem),intent(in) :: this
    integer,intent(in)              :: ithread
    integer,intent(in)              :: tag
    integer,intent(in)              :: com

#ifdef CPP_MPI
    Call this%t_H_eigen%send(ithread,tag,com)
    Call eigen_H_send(ithread,tag,this%H_T,com) 
#else
    continue
#endif
end subroutine

subroutine recv(this,ithread,tag,com)
    use mpi_basic                
    class(t_H_eigen_mem),intent(inout)  :: this
    integer,intent(in)                  :: ithread
    integer,intent(in)                  :: tag
    integer,intent(in)                  :: com

#ifdef CPP_MPI
    Call this%t_H_eigen%recv(ithread,tag,com)
    Call eigen_H_recv(ithread,tag,this%H_T,com) 
#else
    continue
#endif
end subroutine

subroutine bcast(this,comm)
    use mpi_basic                
    class(t_H_eigen_mem),intent(inout)  ::  this
    type(mpi_type),intent(in)           ::  comm
#ifdef CPP_MPI
    Call this%t_H_eigen%bcast(comm)
    Call eigen_H_bcast(comm%id,comm%mas,comm%ismas,this%H_T,comm%com) 
#else
    continue
#endif
end subroutine 

subroutine distribute(this,comm)
    use mpi_basic                
    class(t_H_eigen_mem),intent(inout)  ::  this
    type(mpi_type),intent(in)           ::  comm
#ifdef CPP_MPI
    Call this%t_H_eigen%distribute(comm)
    Call eigen_H_distribute(comm%id,comm%mas,comm%ismas,this%H_T,comm%com)
#else
    continue
#endif
end subroutine 

#endif 
end module
