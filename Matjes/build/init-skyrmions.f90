module m_init_Sk
use m_derived_types
implicit none
private
public :: init_Sk,get_skyrmion
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the starting configuration as an isolated skyrmion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_Sk(io,fname,my_lattice,my_motif,mode_name,start,end)
use m_vector
use m_io_utils
use m_convert
type (lattice), intent(inout) :: my_lattice
type (cell), intent(in) :: my_motif
integer, intent(in) :: io,start,end
character(len=*), intent(in) :: fname,mode_name
! internal variables
integer :: NSkyAdd, N_Pinning
Logical, dimension(:), allocatable :: TabPinning
real(kind=8), dimension(:), allocatable :: tab_XSky,tab_YSky,tab_RSky
real(kind=8), dimension(:), allocatable :: tab_coeffx,tab_coeffy,tab_starx,tab_stary
real(kind=8), dimension(:), allocatable :: chirality
integer :: i
real(kind=8) :: x0,y0,R0,coeffx,coeffy,starx,stary,chi
character(len=30) :: variable_name

variable_name=convert('NSkyAdd_',mode_name)
call get_parameter(io,fname,variable_name,NSkyAdd)
variable_name=convert('N_pinning_',mode_name)
call get_parameter(io,fname,variable_name,N_pinning)

allocate(tab_XSky(NSkyAdd),tab_YSky(NSkyAdd),tab_RSky(NSkyAdd),tab_coeffx(NSkyAdd),tab_coeffy(NSkyAdd),tab_starx(NSkyAdd),tab_stary(NSkyAdd),chirality(NSkyAdd))

variable_name=convert('XSky_',mode_name)
call get_parameter(io,fname,variable_name,NSkyAdd,tab_XSky)
variable_name=convert('YSky_',mode_name)
call get_parameter(io,fname,variable_name,NSkyAdd,tab_YSky)
variable_name=convert('RSky_',mode_name)
call get_parameter(io,fname,variable_name,NSkyAdd,tab_RSky)
variable_name=convert('coeffx_',mode_name)
call get_parameter(io,fname,variable_name,NSkyAdd,tab_coeffx)
variable_name=convert('coeffy_',mode_name)
call get_parameter(io,fname,variable_name,NSkyAdd,tab_coeffy)
variable_name=convert('starx_',mode_name)
call get_parameter(io,fname,variable_name,NSkyAdd,tab_starx)
variable_name=convert('stary_',mode_name)
call get_parameter(io,fname,variable_name,NSkyAdd,tab_stary)
variable_name=convert('chirality_',mode_name)
call get_parameter(io,fname,variable_name,NSkyAdd,chirality)

if (N_pinning.ne.0) allocate(TabPinning(N_pinning))

do i=1,NSkyAdd
   X0=tab_XSky(i)
   Y0=tab_YSky(i)
   R0=tab_RSky(i)
   coeffx=tab_coeffx(i)
   coeffy=tab_coeffy(i)
   starx=tab_starx(i)
   stary=tab_stary(i)
   chi=chirality(i)

   call get_skyrmion(x0,y0,R0,coeffx,coeffy,starx,stary,chi,my_lattice,my_motif,start,end)

enddo

if (N_pinning.ne.0) then
   call get_parameter(io,fname,'pinning',NSkyAdd,tabpinning)
   do i=1,NSkyAdd
      if (tabpinning(i)) call pin_skyrmion(x0,y0,my_lattice,my_motif)
   enddo
endif

if (N_pinning.ne.0) deallocate(TabPinning)
deallocate(tab_XSky,tab_YSky,tab_RSky,tab_coeffx,tab_coeffy,tab_starx,tab_stary,chirality)

end subroutine init_Sk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the set the magnetization for a particular isolated skyrmion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_skyrmion(x0,y0,R0,coeffx,coeffy,starx,stary,chirality,my_lattice,my_motif,start,end)
use m_constants, only : pi
use m_vector, only : norm,distance,phi
use m_get_position
Implicit None
type (lattice), intent(inout) :: my_lattice
type (cell), intent(in) :: my_motif
! (x0, y0) are the coordinates of the center of the skyrmion that we
! want to create0
real(kind=8), intent(in) :: x0, y0
! R0 is the Skyrmion radius
real(kind=8),intent(in) :: R0
! parameter to define if it is a skyrmions, antiskyrmions, Q=...
real(kind=8),intent(in) :: coeffx,coeffy,starx,stary
! define the right or left hand chirality
real(kind=8),intent(in) :: chirality
integer, intent(in) :: start,end
! Internal variables
real(kind=8) :: Theta
real(kind=8) :: Psi
Integer::i_x,i_y,i_z,i_m,dim_lat(3),nmag
real(kind=8) :: dist, net(3,3)
real(kind=8) ::epsnull
real(kind=8) :: diml1,diml2
! Parameter to create a skyrmion woth theta profile according the paper
! of Kubetzka PRB 67, 020401(R) (2003)
real(kind=8) :: widt, cen, dummyvar
! Variables to adjust the variable widt and cen according to Lilley
! criterion for Bloch walls dimension
real(kind=8) :: Inflexparam, ThetaInflex, DThetaInflex
! get the position of the sites
real(kind=8), allocatable :: position(:,:,:,:,:)
integer :: Nx,Ny,Nz
! Variables to detect the closest point on the lattice to the Skyrmion
! Cemter, in order to apply the mask and force its magnetization to be
! down.
real(kind=8) :: MinDist
Integer :: MinIndex(2)

write(6,'(a)') "User defined skyrmion selected"

dim_lat=my_lattice%dim_lat
net=my_lattice%areal
nmag=count(my_motif%i_mom)
Nx=dim_lat(1)
Ny=dim_lat(2)
Nz=dim_lat(3)

allocate(position(3,Nx,Ny,Nz,nmag))
position=0.0d0

if (chirality.gt.0.0d0) then
   write(6,'(a)') "right hand chirality selected"
else
   write(6,'(a)') "left hand chirality selected"
endif

write(6,*) "v_phi(x,y)",starx,stary
write(6,*) "xPi(x,y)",coeffx,coeffy

epsnull = 1.0d-8
! Debug =========================================
!       WRITE(*,*) 'x0 = ', x0, ' y0 = ', y0
! ===============================================
! fixes the shape of the skyrmion
dummyvar=1.0d0

!       widt = 2.0d0*R0/(pi*COSH(0.5d0*dummyvar))
inflexparam = 2.0d0 *ACOSH( 0.5d0*SQRT(2.0d0+SQRT(6.0d0+2.0d0*COSH( 2.0d0*dummyvar ) ) ) )
ThetaInflex = pi(1.0d0)-ASIN( TANH (Inflexparam-dummyvar))-ASIN( TANH (Inflexparam+dummyvar))
DThetaInflex =  -(1.0d0/COSH(Inflexparam-dummyvar) + 1.0d0/COSH(Inflexparam+dummyvar))
widt = 2.d0*R0/(Inflexparam-ThetaInflex/DThetaInflex)

cen = 0.5d0*widt*dummyvar

Theta = 0.0d0
Psi   = 0.0d0
diml1 = DBLE(Dim_lat(1))
diml2 = DBLE(Dim_lat(2))

minIndex = 0
MinDist = 1.0d10
! Check that R0 is positive
If (R0.gt. 0.0d0) then

   call get_position(position,dim_lat,net,my_motif)

   do i_m=1,nmag
      Do i_z=1,dim_lat(3)
         Do i_y=1,dim_lat(2)
            Do i_x=1,dim_lat(1)

        dist=distance(position(1:2,i_x,i_y,i_z,i_m),x0,y0,dim_lat,net)

          If (dist.lt.(1.4d0*R0)) then

             Theta = pi(1.0d0)-Asin( tanh( 2*(dist-cen)/widt ))- &
                               Asin( tanh( 2*(dist+cen)/widt ))
!             Phi   = ACOS( xpoint/dist )

! We parameterize the magnetization vector in spherical coordinates
! although the spatial coordinate is parameterized by cylindrical
! coordinates (See the thesis of Leonov for further details, page 32)

! Be careful. The Psi coordinate, as stated in Leonov's thesis, depends
! on the crystal symmetries.
! to write a skyrmion, x->-x i.e. a symmetry around the y axis
             Psi = Phi(position(1:2,i_x,i_y,i_z,i_m),x0,y0,dim_lat,net)+pi(1.0d0)

             my_lattice%l_modes(i_x,i_y,i_z,i_m)%w(start) = chirality * Sin(Theta) * Cos( starx*Psi + coeffx*pi(1.0d0))
             my_lattice%l_modes(i_x,i_y,i_z,i_m)%w(start+1) = chirality * Sin(Theta) * Sin( stary*Psi + coeffy*pi(1.0d0))
             my_lattice%l_modes(i_x,i_y,i_z,i_m)%w(end) = Cos(Theta)

           else
           cycle
          endif


             if (MinDist.gt.dist) then
                 minIndex = (/i_x,i_y/)
                 MinDist = dist
             endif

            enddo
         enddo
      enddo
   enddo

Endif
deallocate(position)

end subroutine get_skyrmion

! =====================================================================
! SUBROUTINE to Pin a spin in the lattice
! =====================================================================

subroutine pin_skyrmion(X0,Y0,my_lattice,my_motif)
use m_derived_types
use m_get_position
implicit none
type (lattice), intent(in) :: my_lattice
type (cell), intent(in) :: my_motif
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

dim_lat=my_lattice%dim_lat
net=my_lattice%areal
diml1 = DBLE(Dim_lat(1))
diml2 = DBLE(Dim_lat(2))

Nx=dim_lat(1)
Ny=dim_lat(2)
Nz=dim_lat(3)
nmag=count(my_motif%i_mom)

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
