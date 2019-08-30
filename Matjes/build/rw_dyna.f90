module m_dynamic
real(kind=8) :: timestep,damping,temps,hfin(3),max_error
integer :: duration,times,tstart,htimes,hstart,Efreq
real(kind=8) :: storque,hstep(3)
logical :: stmtorque,marche,hsweep,Ffield,i_Efield,stmtemp
real(kind=8), dimension(:,:,:), allocatable :: htor
! integration type
integer :: integtype
end module m_dynamic

subroutine rw_dyna(N,my_lattice)
use m_dynamic
use m_constants
use m_derived_types
use m_io_files_utils
use m_io_utils
implicit none
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: N(3)
! internal
integer :: io
integer :: fin,i_x,i_y,i_z,i
character(len=30) :: str,dummy
real(kind=8) :: xmed,ymed,ri,h,net(3,3)

net=my_lattice%areal

allocate(htor(N(1),N(2),N(3)))

stmtemp=.False.
storque=0.0d0
hstep=0.0d0
damping=0.0d0
hfin=0.0d0
htimes=0
hstart=0
hfin=0.0d0
timestep=0.0d0
temps=0.0d0
stmtorque=.False.
marche=.False.
hsweep=.False.
Ffield=.False.
i_Efield=.False.
Efreq=1

io=open_file_read('input')

call get_parameter(io,'input','integration',integtype)
call get_parameter(io,'input','timestep',timestep)
call get_parameter(io,'input','Efreq',Efreq)
call get_parameter(io,'input','duration',duration)
call get_parameter(io,'input','STMtemp',stmtemp)
call get_parameter(io,'input','damping',damping)

call get_parameter(io,'input','Ffield',Ffield)
call get_parameter(io,'input','Efield',i_Efield)

rewind(io)
do
 read (io,'(a)',iostat=fin) str
 if (fin /= 0) exit
 str= trim(adjustl(str))
 if (len_trim(str)==0) cycle

 if ( str(1:6) == 'Hsweep') then
     backspace(io)
     read(io,*) dummy, hsweep,(hstep(i),i=1,3),htimes,hstart,(hfin(i),i=1,3)
 endif
 if ( str(1:9) == 'stmtorque') then
     backspace(io)
     read(io,*) dummy, stmtorque, storque, ri, h
 endif
enddo

call close_file('input',io)

htor=0.0d0

if (stmtorque) then

   write(*,*) 'SPSTM torque set up'
   xmed=dble(N(1))/2.0d0*net(1,1)+dble(N(2))/2.0d0*net(2,1)
   ymed=dble(N(1))/2.0d0*net(1,2)+dble(N(2))/2.0d0*net(2,2)

   do i_x=1,N(1)
     do i_y=1,N(2)
       do i_z=1,N(3)

       htor(i_x,i_y,i_z)=dexp(-dsqrt((xmed-my_lattice%l_modes(i_x,i_y,i_z,1)%w(1))**2+ &
           (ymed-my_lattice%l_modes(i_x,i_y,i_z,1)%w(2))**2+h**2)/ri)

       enddo
     enddo
   enddo
endif
end subroutine rw_dyna
