module m_init_DW
use m_derived_types
implicit none
private
public :: init_DW
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the starting configuration as a domain wall along the x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_DW(io,fname,lat,ordname,dim_mode,state,init_conf)
    use m_io_utils, only: get_parameter
    use m_util_init, only: get_pos_vec
    use m_constants, only : pi
    integer,intent(in)              :: io       !init-file io-unit
    character(*),intent(in)         :: fname    !init-file name 
    type(lattice), intent(in)       :: lat      !entire lattice containing geometric information
    character(*),intent(in)         :: ordname  !name of the order parameter
    integer,intent(in)              :: dim_mode !dimension of the order parameter in each cell
    real(8),pointer,intent(inout)   :: state(:) !pointer the the order parameter
    real(8), intent(in)             :: init_conf(:)
    ! internal variables
    real(8),allocatable,target :: pos(:)
    real(8),pointer :: pos_3(:,:),state_3(:,:)
    real(8)         :: dw_pos(3),normal(3)  !position on domain wall, normal to domain wall
    real(8),allocatable :: dist(:)
    integer             :: i,j,size_unit_cell
    real(8)             :: length       !length of domain wall
    
    dw_pos=lat%a_sc(1,:)*0.5d0
    normal=[lat%areal(2,2),-lat%areal(2,1),0.0d0]
    normal=normal/norm2(normal)
    length=10*norm2(lat%areal(1,:))
    size_unit_cell=size(init_conf)

    Call get_pos_vec(lat,dim_mode,ordname,pos)
    pos_3(1:3,1:size(pos)/3)=>pos
    state_3(1:3,1:size(pos)/3)=>state

    allocate(dist(size(pos_3,2)),source=0.0d0)
    do i=1,size(pos_3,2)
        pos_3(:,i)=pos_3(:,i)-dw_pos
        dist(i)=dot_product(pos_3(:,i),normal)
        state_3(:,i)=dist(i)
    enddo

    dist=dist*pi/length
    dist=dist+0.5d0*pi

    state=0.d0
    do i=1,size(dist)
        if(dist(i)>=pi)then
            state_3(3,i)=1.0d0
        elseif(dist(i)<=0)then
            state_3(3,i)=-1.0d0
        else
           state_3(1,i)=sin(dist(i))
           state_3(2,i)=0.0d0
           state_3(3,i)=-1.0d0*cos(dist(i))
        endif
        j=mod(i-1,size_unit_cell)+1
        state_3(:,i)=state_3(:,i)*init_conf(j)/abs(init_conf(j))
    enddo
    
    nullify(pos_3,state_3)
    deallocate(pos)
end subroutine
end module m_init_DW
