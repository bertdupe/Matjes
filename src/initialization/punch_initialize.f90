module m_punch_init
use m_derived_types
implicit none
private
public :: punch_lattice
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sets magnetization to zero outside if inserted region
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine punch_lattice(lat,fname_in)
    !loops through all order parameters and checks if there is some punching to be done
    use m_type_lattice,only: number_different_order_parameters,order_parameter_name
    use m_util_init,only: get_pos
    use m_io_files_utils, only: open_file_read,close_file
    type(lattice), intent(inout)    :: lat
    character(len=*), intent(in),optional :: fname_in
    integer :: io

    character(*),parameter              :: fname_default='init.config'
    character(:), allocatable           :: fname
    logical :: exists
    integer :: i
    integer :: dim_mode
    real(8),pointer,contiguous ::  state(:)

    if(present(fname_in))then
        fname=fname_in
    else
        fname=fname_default
    endif
    inquire(file=fname,exist=exists)
    if(.not.exists)then
        write(*,'(/AAA/A/)') 'FILE: "',trim(fname),'" NOT FOUND','SKIPPING init_punch'
        return
    endif
    io=open_file_read(fname)
    do i=1,number_different_order_parameters
        dim_mode=lat%get_order_dim(i,ignore=.true.)
        if(dim_mode<1) cycle
        Call lat%set_order_point(i,state)
        Call punch_order(io,fname,lat,trim(order_parameter_name(i)),dim_mode,state)
        nullify(state)
    enddo
    call close_file(fname,io)
end subroutine


subroutine punch_order(io,fname,lat,ordname,dim_mode,state)
!punches out an area of an order parameter (i.e. sets magnetization to 0 outside of region)
    use m_io_utils,only: get_parameter
    use m_util_init,only: get_pos
    integer,intent(in)              :: io       !init-file io-unit
    character(*),intent(in)         :: fname    !init-file name 
    type(lattice), intent(in)       :: lat      !entire lattice containing geometric information
    character(*),intent(in)         :: ordname  !name of the order parameter
    integer,intent(in)              :: dim_mode !dimension of the order parameter in each cell
    real(8),pointer,intent(inout)   :: state(:) !pointer the the order parameter

    real(8)             :: pos(3),normal(3),radius,width,height,angle
    character(len=30)   :: configuration
    logical             :: do_punch

    real(8),allocatable,target :: pos_arr(:)
    real(8),pointer :: pos_3(:,:),state_pt(:,:)
    integer         :: dim_state

    do_punch=.false.
    call get_parameter(io,fname,'punch_'//ordname,do_punch)
    if(.not.do_punch) return
    !get parameters for punch
    pos=0.0d0
    normal=[1.0d0,0.0d0,0.0d0]
    radius=Huge(1.0d0)
    width=Huge(1.0d0)
    height=Huge(1.0d0)
    angle=0.0d0
    configuration=""
    !general input parameter for all punches? Could also be removed but might break older input
    call get_parameter(io,fname,'punch_config',configuration)
    call get_parameter(io,fname,'punch_pos'   ,3,pos)
    call get_parameter(io,fname,'punch_radius',radius)
    call get_parameter(io,fname,'punch_width' ,width)
    call get_parameter(io,fname,'punch_height',height)
    call get_parameter(io,fname,'punch_angle' ,angle)
    call get_parameter(io,fname,'punch_normal' ,normal)

    !specialized input parameter for this order parameter?
    call get_parameter(io,fname,'punch_config_'//ordname,configuration)
    call get_parameter(io,fname,'punch_pos_'   //ordname,3,pos)
    call get_parameter(io,fname,'punch_radius_'//ordname,radius)
    call get_parameter(io,fname,'punch_width_' //ordname,width)
    call get_parameter(io,fname,'punch_height_'//ordname,height)
    call get_parameter(io,fname,'punch_angle_' //ordname,angle)
    call get_parameter(io,fname,'punch_normal_' //ordname,normal)

    Call get_pos(lat,dim_mode,ordname,pos_arr,dim_state)
    state_pt(1:dim_state,1:size(state)/dim_state)=>state
    pos_3(1:3,1:size(pos_arr)/3)=>pos_arr
    if(size(pos_3,2)/=size(state_pt,2)) ERROR STOP "unexpected shapes of position and state array, programming mistake?"

!warning, the punch routines can change the positions(pos_3,pos_arr)
    select case (adjustl(configuration))
        case('circle')
            call punch_texture_circle(pos,radius,pos_3,state_pt)
        case('hexagon')
            Call punch_texture_hexagon(pos,width,angle,pos_3,state_pt)
        case('rectangle')
            Call punch_texture_rectangle(pos,width,height,angle,pos_3,state_pt)
        case('heaviside')
            Call punch_texture_heaviside(pos,normal,pos_3,state_pt)
        case("")
            write(*,'(//3A)') 'Failing to do punch_order for: ',ordname
            write(*,'(A)') 'Punching is requested (punch_), but no punch_config is set'
            STOP "Could not punch order parameter"
        case default
            write(*,'(//2A)') "Failing to do punch_order for: ",ordname
            write(*,'(3A)') 'Provided punch_config: "',trim(configuration),'" is not implemented'
            STOP "Could not punch order parameter"
    end select
end subroutine 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ROUTINES IMPLEMENTING PARTICULAR SHAPES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine punch_texture_circle(center,radius,pos,state)
!sets all order parameters to 0 that are outside of a cirle with radius radius around then center
!(only counting the distance in xy space)
    real(8), intent(in)         :: center(3)    !center of circle region
    real(8),intent(in)          :: radius !radius of circular region
    real(8),intent(inout)       :: pos(:,:)
    real(8),intent(inout)       :: state(:,:)
    ! Internal variables
    logical :: remove(size(pos,2))
    integer :: i
   
    Call center_pos(center,pos)
    remove=norm2(pos(1:2,:),1)>radius
    forall(i=1:size(pos,2), remove(i) ) state(:,i)=0.0d0
end subroutine

subroutine punch_texture_rectangle(center,width,height,angle_in,pos,state)
!sets all order parameters to 0 that are outside of an rectangle with width, height around center and rotated by angle_in(along z)
!(only counting the distance in xy space)
    !shape parameters:
    real(8),intent(in)          :: center(3)
    real(8),intent(in)          :: width
    real(8),intent(in)          :: height
    real(8),intent(in)          :: angle_in !angle in degree
    !lattice parameters:
    real(8),intent(inout)       :: pos(:,:)
    real(8),intent(inout)       :: state(:,:)
    ! Internal variables
    integer                 :: i
    logical                 :: remove(size(pos,2))

    Call center_pos(center,pos)
    if(angle_in/=0.0d0) Call rotate_pos(angle_in,pos)

    remove=abs(pos(1,:))>width*0.5d0.or.abs(pos(2,:))>height*0.5d0
    forall(i=1:size(pos,2), remove(i)) state(:,i)=0.0d0
end subroutine

subroutine punch_texture_heaviside(center,normal,pos,state)
!sets all order parameters to 0 for which (position-center) projected along the normal vector is larger than 0
    !shape parameters:
    real(8),intent(in)          :: center(3)    !center position  of heaviside function
    real(8),intent(in)          :: normal(3)    !normal direction of heaviside function
    !lattice parameters:
    real(8),intent(inout)       :: pos(:,:)
    real(8),intent(inout)       :: state(:,:)
    ! Internal variables
    integer                 :: i
    logical                 :: remove(size(pos,2))

    Call center_pos(center,pos)
    remove=matmul(normal,pos)>0.0d0
    forall(i=1:size(pos,2), remove(i)) state(:,i)=0.0d0
end subroutine


subroutine punch_texture_hexagon(center,width,angle_in,pos,state)
    !shape parameters:
    real(8), intent(in)         :: center(3)
    real(8),intent(in)          :: width,angle_in !minor width 
    !lattice parameters:
    real(8),intent(inout)       :: pos(:,:)
    real(8),intent(inout)       :: state(:,:)

    logical     :: remove(size(pos,2))
    integer     :: i
    
    Call center_pos(center,pos)
    if(angle_in/=0.0d0) Call rotate_pos(angle_in,pos)

    !only check upper quadrant (absolute values)
    pos=abs(pos)/width 
    !check if outside maximal rectangle
    remove=pos(1,:)>=sqrt(3.0d0)*0.5d0.or.pos(2,:)>=1.0d0
    !check is below 
    remove=remove.or.pos(2,:)-1.0d0+pos(1,:)/sqrt(3.0d0)>0.0d0
    
    forall(i=1:size(pos,2), remove(i)) state(:,i)=0.0d0
end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!SMALL HELPER FUNCTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure subroutine center_pos(center,pos)
    real(8),intent(in)          :: center(3)
    real(8),intent(inout)       :: pos(:,:)
    integer     ::  i

    do i=1,size(pos,2)
        pos(:,i)=pos(:,i)-center
    enddo
end subroutine

pure subroutine rotate_pos(angle_in,pos)
    use m_constants, only : pi
    real(8),intent(in)          :: angle_in !angle in degree
    real(8),intent(inout)       :: pos(:,:)
    real(8)     :: angle,rot(3,3)

    !rotate position
    angle=angle_in*2.0d0*pi/360.0
    rot(:,1)=[ cos(angle),-sin(angle),0.0d0]
    rot(:,2)=[ sin(angle), cos(angle),0.0d0]
    rot(:,3)=[      0.0d0,      0.0d0,1.0d0]
    pos=matmul(rot,pos)
end subroutine
end module 
