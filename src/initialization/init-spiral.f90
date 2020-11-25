module m_init_spiral
use m_derived_types
implicit none

private
public :: init_spiral

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the starting configuration as a spin spiral
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine init_spiral(io,fname,lat,ordname,dim_mode,state)
    use m_io_utils, only: get_parameter
    use m_util_init, only: get_pos_vec
    integer,intent(in)              :: io       !init-file io-unit
    character(*),intent(in)         :: fname    !init-file name 
    type(lattice), intent(in)       :: lat      !entire lattice containing geometric information
    character(*),intent(in)         :: ordname  !name of the order parameter
    integer,intent(in)              :: dim_mode !dimension of the order parameter in each cell
    real(8),pointer,intent(inout)   :: state(:) !pointer the the order parameter

    real(8)         :: qvec(3),Rq(3),Iq(3)
    real(8),allocatable,target :: pos(:)
!    real(8),allocatable ::  position(:)
    real(8),pointer :: pos_3(:,:),state_3(:,:)
    integer         :: i
    integer         :: nmag
   
    qvec=0.0d0
    Rq=[0.0d0,0.0d0,1.0d0]
    Iq=[1.0d0,0.0d0,0.0d0]
    
    call get_parameter(io,fname,'qvec_'//ordname,3,qvec)
    qvec=matmul(qvec,lat%astar)

    call get_parameter(io,fname,'Rq_'//ordname,3,Rq,1.0d0)
    Rq=matmul(Rq,lat%areal)
    Rq=Rq/norm2(Rq)
    
    call get_parameter(io,fname,'Iq_'//ordname,3,Iq,1.0d0)
    Iq=matmul(Iq,lat%areal)
    Iq=Iq/norm2(Iq)

    Call get_pos_vec(lat,dim_mode,ordname,pos)

    pos_3(1:3,1:size(pos)/3)=>pos
    state_3(1:3,1:size(pos)/3)=>state
    do i=1,size(state_3,2)
        state_3(:,i)=(cos(dot_product(qvec,pos_3(:,i)))*Rq+ &
                      sin(dot_product(qvec,pos_3(:,i)))*Iq)
    enddo
    nullify(pos_3,state_3)
    deallocate(pos)
end subroutine 
end module m_init_spiral
