module m_init_bobber
implicit none
private
public :: init_bobber
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the starting configuration as an isolated skyrmion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_bobber(io,fname,lat,ordname,dim_mode,state,init_conf)
    use m_io_utils,only: get_parameter
    use m_util_init, only: get_pos_vec, get_skyrmion
    use m_derived_types,only: lattice
    integer,intent(in)              :: io       !init-file io-unit
    character(*),intent(in)         :: fname    !init-file name
    type(lattice), intent(in)       :: lat      !entire lattice containing geometric information
    character(*),intent(in)         :: ordname  !name of the order parameter
    integer,intent(in)              :: dim_mode !dimension of the order parameter in each cell
    real(8),pointer,intent(inout)   :: state(:) !pointer the the order parameter
    real(8), intent(in)             :: init_conf(:)
    ! internal variables
    real(8),allocatable,target          :: pos(:)
    real(8),pointer                     :: pos_3(:,:),state_3(:,:)
    real(8)                             :: pos_sky(3)
    real(8), dimension(:), allocatable  :: tab_XSky,tab_YSky
    real(8)         :: R0,coeffx,coeffy,starx,stary,chi,qSklattice
    integer         :: NSkyAdd
    integer         :: Nsk(2),i

    qSklattice=0.02d0
    call get_parameter(io,fname,'qSklattice',qSklattice) !legecy without ordname?
    call get_parameter(io,fname,'qSklattice_'//ordname,qSklattice) !legecy without ordname?

    Nsk=nint(qSklattice*lat%dim_lat(1:2))
    NSkyAdd=product(Nsk)

    allocate(tab_XSky(NSkyAdd),tab_YSky(NSkyAdd),source=0.d0)

    R0=5.0d0
    chi=-1.0d0
    coeffx=1.0d0; coeffy=1.0d0
    starx=1.0d0 ; stary=1.0d0

    call find_XYsky(tab_XSky,tab_YSky,NSkyAdd,qSklattice,lat%dim_lat,lat%areal)

    call get_parameter(io,fname,'RSky_'//ordname,R0)
    call get_parameter(io,fname,'coeffx_'//ordname,coeffx)
    call get_parameter(io,fname,'coeffy_'//ordname,coeffy)
    call get_parameter(io,fname,'starx_'//ordname,starx)
    call get_parameter(io,fname,'stary_'//ordname,stary)
    call get_parameter(io,fname,'chirality_'//ordname,chi)

    Call get_pos_vec(lat,dim_mode,ordname,pos)
    pos_3(1:3,1:size(pos)/3)=>pos
    state_3(1:3,1:size(pos)/3)=>state

    do i=1,NSkyAdd
        pos_sky=[tab_XSky(i),tab_YSky(i),0.0d0]
        call get_skyrmion(pos_sky,R0,coeffx,coeffy,starx,stary,chi,pos_3,state_3,lat,init_conf)

    enddo

    deallocate(tab_XSky,tab_YSky)

end subroutine

!!! find the positions of all the skyrmions in the lattice

subroutine find_XYsky(XSky,YSky,Nadd,qskx,dim_lat,net)
use m_vector, only : norm
implicit none
integer, intent(in) :: Nadd,dim_lat(:)
real(kind=8), intent(in) :: qskx,net(:,:)
real(kind=8), intent(inout) :: XSky(Nadd),YSky(Nadd)
! internal
integer :: i_x,i_y,N_x,N_y
integer :: k

N_x=nint(qskx*dim_lat(1))
N_y=nint(qskx*dim_lat(2))

k=0
do i_x=1,N_x
   do i_y=1,N_y
      k=k+1
      XSky(k)=dble(dim_lat(1)*(i_x-1))/dble(N_x)*net(1,1)+dble(dim_lat(2)*(i_y-1))/dble(N_y)*net(2,1)
      YSky(k)=dble(dim_lat(1)*(i_x-1))/dble(N_x)*net(1,2)+dble(dim_lat(2)*(i_y-1))/dble(N_y)*net(2,2)
   enddo
enddo

if (k.ne.Nadd) then
   write(6,'(a)') 'problem in the lattice of skyrmion'
   write(6,'(a)') 'check SkX_utils routine'
   stop
endif

end subroutine find_XYsky

end module m_init_bobber
