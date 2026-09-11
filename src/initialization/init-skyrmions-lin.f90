module m_init_sky_lin
use m_derived_types
implicit none
private
public :: init_sky_lin
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! sets the magnetization according the skyrmion structure of:
!   Garnier, M., Mesaros, A., & Simon, P. (2019). 
!   Topological superconductivity with deformable magnetic skyrmions.
!   Communications Physics, 2(1), 126.
!   https://doi.org/10.1038/s42005-019-0226-5
! i.e.
! relative to some center we sets
! n(r)=[sin(f(r))*cos(q*phi),sin(f(r))*sin(q*phi),cos(f(r))]
! with f(r) being a linear function with a given slope(multplied with 2*pi from the input),
! q being an integer giving the azimuthal winding number,
! and phi being the polar angle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_sky_lin(io,fname,lat,ordname,dim_mode,state,init_conf)
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

    real(8),allocatable,target :: pos(:)
    real(8),pointer :: pos_3(:,:),state_3(:,:)
    real(8)     ::  f_r,qphi
    integer     :: i,j,size_unit_cell

    real(8)     ::  slope,center(3)
    integer     ::  q

    !get parameters for punch
    center=0.0d0
    q=1
    slope=10.0
    size_unit_cell=size(init_conf)
    !older input without ordername? (legacy?)
    call get_parameter(io,fname,'skylin_pos',3,center)
    call get_parameter(io,fname,'skylin_slope',slope)
    call get_parameter(io,fname,'skylin_q',q)
    call get_parameter(io,fname,'skylin_pos_'//ordname,3,center)
    call get_parameter(io,fname,'skylin_slope_'//ordname,slope)
    call get_parameter(io,fname,'skylin_q_'//ordname,q)
    slope=2.0d0*pi/slope
    
    Call get_pos_vec(lat,dim_mode,ordname,pos)
    pos_3(1:3,1:size(pos)/3)=>pos
    state_3(1:3,1:size(pos)/3)=>state

    do i=1,size(pos_3,2)
        pos_3(:,i)=pos_3(:,i)-center
    enddo
    
    do i=1,size(pos_3,2)
        f_r=norm2(pos_3(:,i))*slope
        qphi=atan2(pos_3(2,i),pos_3(1,i))*real(q,8)
        state_3(1,i)= sin(f_r)*cos(qphi)
        state_3(2,i)= sin(f_r)*sin(qphi)
        state_3(3,i)= cos(f_r)
        j=mod(i-1,size_unit_cell)+1
        state_3(:,i)=state_3(:,i)*init_conf(j)/abs(init_conf(j))
    enddo

    nullify(pos_3,state_3)
    deallocate(pos)
end subroutine 
end module 
