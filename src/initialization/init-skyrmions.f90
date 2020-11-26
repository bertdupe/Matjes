module m_init_Sk
use m_derived_types
implicit none
private
public :: init_Sk
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the starting configuration as an isolated skyrmion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_Sk(io,fname,lat,ordname,dim_mode,state)
    use m_io_utils, only: get_parameter
    use m_util_init, only: get_pos_vec,get_skyrmion
    integer,intent(in)              :: io       !init-file io-unit
    character(*),intent(in)         :: fname    !init-file name 
    type(lattice), intent(in)       :: lat      !entire lattice containing geometric information
    character(*),intent(in)         :: ordname  !name of the order parameter
    integer,intent(in)              :: dim_mode !dimension of the order parameter in each cell
    real(8),pointer,intent(inout)   :: state(:) !pointer the the order parameter

    integer :: NSkyAdd, N_Pinning
    Logical, dimension(:), allocatable :: TabPinning
    real(8), dimension(:), allocatable :: tab_XSky,tab_YSky,tab_RSky
    real(8), dimension(:), allocatable :: tab_coeffx,tab_coeffy,tab_starx,tab_stary
    real(8), dimension(:), allocatable :: chirality
    integer :: i
    real(8) :: R0,coeffx,coeffy,starx,stary,chi
    character(len=30) :: variable_name

    real(8),allocatable,target :: pos(:)
    real(8),pointer :: pos_3(:,:),state_3(:,:)
    integer         :: nmag
    real(8)         :: pos_sky(3)
    
    !!!! default variables !!!!
    NSkyAdd=1
    N_pinning=0
    
    call get_parameter(io,fname,'NSkyAdd_'//ordname,NSkyAdd)
    
    allocate(tab_XSky(NSkyAdd),tab_YSky(NSkyAdd),tab_RSky(NSkyAdd),tab_coeffx(NSkyAdd),tab_coeffy(NSkyAdd),tab_starx(NSkyAdd),tab_stary(NSkyAdd),chirality(NSkyAdd),source=1.0d0)
    tab_XSky=50.0d0
    tab_YSky=50.0d0
    tab_RSky=5.0d0
    chirality=-1.0d0
    !!!!
    
    call get_parameter(io,fname,'XSky_'//ordname,NSkyAdd,tab_XSky)
    call get_parameter(io,fname,'YSky_'//ordname,NSkyAdd,tab_YSky)
    call get_parameter(io,fname,'RSky_'//ordname,NSkyAdd,tab_RSky)
    call get_parameter(io,fname,'coeffx_'//ordname,NSkyAdd,tab_coeffx)
    call get_parameter(io,fname,'coeffy_'//ordname,NSkyAdd,tab_coeffy)
    call get_parameter(io,fname,'starx_'//ordname,NSkyAdd,tab_starx)
    call get_parameter(io,fname,'stary_'//ordname,NSkyAdd,tab_stary)
    call get_parameter(io,fname,'chirality_'//ordname,NSkyAdd,chirality)

    Call get_pos_vec(lat,dim_mode,ordname,pos)
    pos_3(1:3,1:size(pos)/3)=>pos
    state_3(1:3,1:size(pos)/3)=>state
    
    do i=1,NSkyAdd
       R0=tab_RSky(i)
       coeffx=tab_coeffx(i)
       coeffy=tab_coeffy(i)
       starx=tab_starx(i)
       stary=tab_stary(i)
       chi=chirality(i)
       pos_sky=[tab_XSky(i),tab_YSky(i),0.0d0]
       call get_skyrmion(pos_sky,R0,coeffx,coeffy,starx,stary,chi,pos_3,state_3,lat)
    enddo
    
    call get_parameter(io,fname,'N_pinning_'//ordname,N_pinning)
    if (N_pinning.ne.0) then
        ERROR STOP "pinning has to be implemented in new version" ! I don't really understand what it does
       !allocate(TabPinning(N_pinning))
       !call get_parameter(io,fname,'pinning',NSkyAdd,tabpinning)
       !do i=1,NSkyAdd
       !   if (tabpinning(i)) call pin_skyrmion(x0,y0,my_lattice,my_motif)
       !enddo
    endif
!    
!    if (N_pinning.ne.0) deallocate(TabPinning)
    deallocate(tab_XSky,tab_YSky,tab_RSky,tab_coeffx,tab_coeffy,tab_starx,tab_stary,chirality)

end subroutine 


! =====================================================================
! SUBROUTINE to Pin a spin in the lattice
! =====================================================================

subroutine pin_skyrmion(X0,Y0,my_lattice,my_motif)
use m_derived_types
use m_get_position
implicit none
type (lattice), intent(in) :: my_lattice
type(t_cell), intent(in) :: my_motif
real(kind=8), intent(in) ::X0,Y0
! internal variables
Integer::i_x,i_y,i_z,i_m
real(kind=8) :: dist, xpoint, ypoint, xpboc1, xpboc2, ypboc1, ypboc2
real(kind=8) :: DBottom,DTop,DLeft,DRight,DC1,DC2,DC3,DC4
real(kind=8) :: diml1, diml2
! Variables to detect the closest point on the lattice to the Skyrmion
! Cemter, in order to apply the mask and force its magnetization to be
! down.
real(kind=8) :: distprevious, MinDist, net(3,3)
Integer :: MinIndex(2), dim_lat(3),nmag,Nx,Ny,Nz
! get the position of the sites
real(kind=8), allocatable :: position(:,:,:,:,:)

ERROR STOP "pin_skyrmions not updated to Patrick-version"

dim_lat=my_lattice%dim_lat
net=my_lattice%areal
diml1 = real(Dim_lat(1),8)
diml2 = real(Dim_lat(2),8)

Nx=dim_lat(1)
Ny=dim_lat(2)
Nz=dim_lat(3)
nmag=count(my_motif%atomic(:)%moment.gt.0.0d0)

allocate(position(3,Nx,Ny,Nz,nmag))
position=0.0d0

call get_position(position,dim_lat,net,my_motif)

minIndex = 0
distprevious = 1.0d150
! Check that R0 is positive

do i_m=1,nmag
   do i_z=1,dim_lat(3)
      do i_y=1,dim_lat(2)
         do i_x=1,dim_lat(1)

       xpoint = position(1,i_x,i_y,i_z,i_m)-x0
       ypoint = position(2,i_x,i_y,i_z,i_m)-y0
       dist = SQRT( xpoint**2 + ypoint**2 )
! Debug =======================================================
!      WRITE(*,*) '==========================================='
!       WRITE(*,*)  ' x = ', Spin(i,1), &
!                   ' y = ', Spin(i, 2),' dist = ', dist
! =============================================================
! To test for periodic boundary conditions
       xpboc1 = xpoint-diml1
       ypboc1 = ypoint-diml2
       xpboc2 = xpoint+diml1
       ypboc2 = ypoint+diml2
       DBottom = SQRT(xpoint**2 +ypboc1**2 )
       DTop    = SQRT(xpoint**2 +ypboc2**2 )
       DLeft   = SQRT(xpboc1**2 +ypoint**2 )
       DRight  = SQRT(xpboc2**2 +ypoint**2 )
       DC1     = SQRT(xpboc1**2 +ypboc1**2 )
       DC2     = SQRT(xpboc2**2 +ypboc2**2 )
       DC3     = SQRT(xpboc1**2 +ypboc2**2 )
       DC4     = SQRT(xpboc2**2 +ypboc1**2 )

             MinDist = MIN(dist, DBottom, Dtop, DLeft, DRight, &
                           DC1, DC2, DC3, DC4)
             if (distprevious.gt.mindist) then
                 minIndex = (/i_x,i_y/)
                 distprevious = MinDist
             endif

         enddo
      enddo
   enddo
enddo


deallocate(position)

end subroutine pin_skyrmion

end module m_init_Sk
