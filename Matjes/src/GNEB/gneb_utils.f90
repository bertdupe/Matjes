module m_gneb_utils
use m_vector, only : calc_ang,norm
use m_path
use m_write_spin
use m_local_energy
use m_createspinfile
use m_basic_types, only : vec_point, vec
use m_derived_types, only : io_parameter,lattice
use m_gneb_parameters, only : do_norm_rx,en_zero
use m_eval_Beff
use m_rotation, only : rotation_axis,rotate
use m_operator_pointer_utils
use m_energyfield
use m_spline
use m_io_gneb
use m_projection
use m_tangent
use m_local_energy

private
public :: find_path,find_path_ci,find_SP,find_SP_conf
contains

subroutine find_path(nim,N_cell,dt,mass,kappa,ftol,itrmax,every,rx,ene,dene,path,my_lattice,io_simu)
implicit none
type(io_parameter), intent(in) :: io_simu
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: nim
real(kind=8), intent(in) :: mass             !> mass of the point
real(kind=8), intent(in) :: dt               !> timestep
real(kind=8), intent(in) :: kappa            !> spring constant
real(kind=8), intent(in) :: ftol            !> spring constant
integer, intent(in) :: itrmax,N_cell
real(kind=8), target, intent(inout) :: path(:,:,:)
integer, intent(in) :: every !< Save path every 'every' step
real(kind=8), intent(out) :: rx(nim) !< Reaction coordinate
real(kind=8), intent(out) :: ene(nim) !< Energy of the images
real(kind=8), intent(out) :: dene(nim) !< Derivative of the energy with respect to reaction coordinate
! internal
integer :: iomp,i_nim,itr
real(kind=8), allocatable, dimension(:,:,:) :: vel,fxyz1,fxyz2,ax,coo
real(kind=8), allocatable, dimension(:,:) :: tau_i,tau,ang
real(kind=8), allocatable, dimension(:) :: ftmp,veltmp
real(kind=8) :: fchk,u0,u(nim),fpp(nim)
real(kind=8) :: pathlen(nim)
real(kind=8) :: fv,fd,fp,E_int,norm_local
integer :: i,imax,dim_mode,dim_mode_mag,maximum(3)
logical :: found
!!!!!!!!!!! allocate the pointers to find the path
type(vec_point),allocatable,dimension(:,:) :: all_mode_path,magnetic_mode_path
!!!!!!!!!!!!!!!

! initialize all pointers
allocate(all_mode_path(N_cell,nim),magnetic_mode_path(N_cell,nim))
found=.false.
do i_nim=1,nim
   call associate_pointer(magnetic_mode_path(:,i_nim),path(:,:,i_nim),'magnetic',found)
enddo
call get_B_matrix(my_lattice%dim_mode)
call get_E_matrix(my_lattice%dim_mode)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! allocate the pointers for the B-field and the energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i_nim=1,nim
  do i=1,N_cell
    all_mode_path(i,i_nim)%w=>path(:,i,i_nim)
  enddo
enddo

dim_mode_mag=size(magnetic_mode_path(1,1)%w)
allocate(vel(dim_mode_mag,N_cell,nim),fxyz1(dim_mode_mag,N_cell,nim),fxyz2(dim_mode_mag,N_cell,nim),ax(dim_mode_mag,N_cell,nim),coo(dim_mode_mag,N_cell,nim),ang(N_cell,nim))
allocate(tau_i(dim_mode_mag,N_cell),tau(dim_mode_mag,N_cell))
vel=0.0d0
fxyz1=0.0d0
fxyz2=0.0d0
ax=0.0d0
coo=0.0d0
ang=0.0d0
tau_i=0.0d0
tau=0.0d0

fchk=1d0+ftol
itr=1
dim_mode=my_lattice%dim_mode
allocate(ftmp(dim_mode),veltmp(dim_mode_mag))
ftmp=0.0d0
veltmp=0.0d0

call the_path(nim,magnetic_mode_path,pathlen)
   
do i_nim=1,nim

   u(i_nim) = 0d0
   do iomp=1,N_cell

      coo(:,iomp,i_nim) = magnetic_mode_path(iomp,i_nim)%w

      call calculate_Beff(ftmp,iomp,all_mode_path(:,i_nim))

      call project_force(ftmp(1:3),magnetic_mode_path(iomp,i_nim)%w,fxyz1(:,iomp,i_nim))

      call local_energy(E_int,iomp,all_mode_path(:,i_nim))

      u(i_nim) = u(i_nim) + E_int

   enddo
enddo
   

                        
if (en_zero=='I') then
      u0 = u(1)
elseif (en_zero=='F') then
      u0 = u(nim)
else
      u0 = 0d0
end if
   
!   !fname_work = 'triven_'//trim(adjustl(simid))
!      !open(777,file = fname_work, access = 'sequential', action = 'write',status = 'replace')

! Calculation of the tangent trajectory of constant energy
call tang(nim,1,coo,tau_i)
fpp = 0d0
fpp(1)=sum( fxyz1(:,:,1) * tau_i )

! calculation of the tangent trajectory for all the images
do i_nim=2,nim-1
   call tang(i_nim,coo,u,tau)

   call tang(nim,i_nim,coo,tau_i)
         
   fp = 0d0
   fpp(i_nim) = 0d0

   fpp(i_nim)=sum( fxyz1(:,:,i_nim) *tau_i )
   fp = sum( fxyz1(:,:,i_nim) * tau )

         
   do iomp=1,N_cell
     fxyz1(:,iomp,i_nim) = fxyz1(:,iomp,i_nim) - tau(:,iomp)*fp + kappa*tau(:,iomp)*(pathlen(i_nim+1)+pathlen(i_nim-1)-2d0*pathlen(i_nim))
   end do

end do

call tang(nim,nim,coo,tau_i)
fpp(nim) = 0d0

fpp(nim) = sum( fxyz1(:,:,nim)*tau_i )

      
fchk=0d0
imax = 1
fchk=maxval( abs(fxyz1) )
maximum=maxloc(fxyz1)
imax=maximum(3)


itr=1
      
open(99,file = 'force_mep.txt',access = 'sequential',action='write',status='replace')
close(99)
      
open(99,file = 'force_mep.txt', access = 'sequential', action = 'write',status = 'old',position = 'append')
write(99,'(i12,a,es16.8E3,a,i3)',advance = 'no') itr,'   ',fchk,'   ',imax
      
close(99)
call write_en(nim,pathlen,u-u0,-fpp,pathlen(nim),'en_path.in',do_norm_rx)
      
call write_path(path)
      
      
!======================MAIN LOOP============================
!
! This part is the integrator
!
!
write(6,'(/a)') 'Main loop'
write(6,'(2(a,I8))') 'maximum iteration  ', itrmax, ' iteration number', itr
write(6,'(2(a,f16.12)/)') 'distance to tolerance ', fchk, ' tolerance', ftol

do while ((fchk.gt.ftol).and.(itr.le.itrmax))
   do i_nim=2,nim-1
      u(i_nim) = 0d0
      do iomp=1,N_cell
            
         !print *,'i)nim:',i_nim
         ax(:,iomp,i_nim) = rotation_axis(coo(:,iomp,i_nim),fxyz1(:,iomp,i_nim))
                        
         coo(:,iomp,i_nim) = magnetic_mode_path(iomp,i_nim)%w+vel(:,iomp,i_nim)*dt+0.5d0*fxyz1(:,iomp,i_nim)/mass*dt*dt

         norm_local=norm(coo(:,iomp,i_nim))
         coo(:,iomp,i_nim)=coo(:,iomp,i_nim)/norm_local
         ang(iomp,i_nim) = calc_ang(coo(:,iomp,i_nim),magnetic_mode_path(iomp,i_nim)%w)
         magnetic_mode_path(iomp,i_nim)%w = coo(:,iomp,i_nim)
         !print *,'ang:',ang(iomp,i_nim)
      end do
   end do
            
         
   do i_nim=2,nim-1
      u(i_nim) = 0d0

      do iomp=1,N_cell

         call calculate_Beff(ftmp,iomp,all_mode_path(:,i_nim))

         call project_force(ftmp(1:3),magnetic_mode_path(iomp,i_nim)%w,fxyz1(:,iomp,i_nim))

         call local_energy(E_int,iomp,all_mode_path(:,i_nim))

         u(i_nim) = u(i_nim) + E_int

      enddo

   end do
         
         
   call the_path(nim,magnetic_mode_path,pathlen)
         
   !print *,'pathlen:',pathlen


   do i_nim=2,nim-1
     call tang(i_nim,coo,u,tau)
     call tang(nim,i_nim,coo,tau_i)
         
     fp = 0d0
     fpp(i_nim) = 0d0

     fpp(i_nim)=sum( fxyz2(:,:,i_nim) * tau_i )
     fp = sum( fxyz2(:,:,i_nim) * tau )

     do iomp=1,N_cell
       fxyz2(:,iomp,i_nim) = fxyz2(:,iomp,i_nim) - tau(:,iomp)*fp + kappa*tau(:,iomp)*(pathlen(i_nim+1)+pathlen(i_nim-1)-2d0*pathlen(i_nim))
     enddo
   enddo
   !stop

   do i_nim=2,nim-1
      do iomp=1,N_cell
         call rotate(vel(:,iomp,i_nim),ax(:,iomp,i_nim),ang(iomp,i_nim),veltmp)
         call rotate(fxyz1(:,iomp,i_nim),ax(:,iomp,i_nim),ang(iomp,i_nim),ftmp(1:3))

         vel(:,iomp,i_nim) = veltmp + 0.5d0*(ftmp(1:3)+fxyz2(:,iomp,i_nim))/mass*dt

      enddo
   enddo
         
         
   fv = 0d0
   fd = 0d0

   fv= sum( vel * fxyz2 )
   fd= sum( fxyz2**2 )
         
   if (fv<0d0) then
     vel = 0d0
   else
     vel=fxyz2*fv/fd
   end if
         
         
         
         
   imax=1
   fchk=0d0
   fxyz1 = fxyz2

   fchk=maxval( abs(fxyz1) )
   maximum=maxloc(fxyz1)
   imax=maximum(3)
         
         
         
   call tang(nim,1,coo,tau_i)
   fpp(1) = sum( fxyz1(:,:,1)*tau_i )
         
         
   call tang(nim,nim,coo,tau_i)
   fpp(nim) = sum( fxyz1(:,:,nim)*tau_i )
         
         
   itr=itr+1

   if (mod(itr,every).eq.0) then
      call prn_gneb_progress(itr, itrmax, fchk, imax,'N',0)
      
      open(99,file = 'force_mep.txt', access = 'sequential', action = 'write',status = 'old',position = 'append')
      write(99,'(i12,a,es16.8E3,a,i3)',advance = 'no') itr,'   ',fchk,'   ',imax
      close(99)
            
      call write_en(nim,pathlen,u-u0,-fpp,pathlen(nim),'en_path.out',do_norm_rx)
      call write_path(path)
            
         
   end if
end do

if (itr>itrmax) then
   write(6,'(a)') 'WARNING: exceeded maximum iterations in GNEB'
end if
   
      

ene = u-u0
dene = -fpp
rx = pathlen
   

end subroutine find_path





!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! climbing image

subroutine find_path_ci(nim,N_cell,dt,mass,kappa,ftol,itrmax,every,rx,ene,dene,ci,path,my_lattice,io_simu)
implicit none
type(io_parameter), intent(in) :: io_simu
type(lattice), intent(in) :: my_lattice
integer, intent(in) :: nim
real(kind=8), intent(in) :: mass             !> mass of the point
real(kind=8), intent(in) :: dt               !> timestep
real(kind=8), intent(in) :: kappa            !> spring constant
real(kind=8), intent(in) :: ftol            !> spring constant
integer, intent(in) :: itrmax,N_cell
real(kind=8), target, intent(inout) :: path(:,:,:)
integer, intent(in) :: every !< Save path every 'every' step
real(kind=8), intent(out) :: rx(nim) !< Reaction coordinate
real(kind=8), intent(out) :: ene(nim) !< Energy of the images
real(kind=8), intent(out) :: dene(nim) !< Derivative of the energy with respect to reaction coordinate
 integer, intent(out) :: ci
! internal
integer :: iomp,i_nim,itr
real(kind=8), allocatable, dimension(:,:,:) :: vel,fxyz1,fxyz2,ax,coo
real(kind=8), allocatable, dimension(:,:) :: tau_i,tau,ang
real(kind=8), allocatable, dimension(:) :: ftmp,veltmp
real(kind=8) :: fchk,u0,u(nim),fpp(nim)
real(kind=8) :: pathlen(nim)
real(kind=8) :: fv,fd,fp,E_int,norm_local
integer :: i,imax,dim_mode,dim_mode_mag,maximum(3)
logical :: found
!!!!!!!!!!! allocate the pointers to find the path
type(vec_point),allocatable,dimension(:,:) :: all_mode_path,magnetic_mode_path
!!!!!!!!!!!!!!!

! initialize all pointers
allocate(all_mode_path(N_cell,nim),magnetic_mode_path(N_cell,nim))
found=.false.
do i_nim=1,nim
   call associate_pointer(magnetic_mode_path(:,i_nim),path(:,:,i_nim),'magnetic',found)
enddo
call get_B_matrix(my_lattice%dim_mode)
call get_E_matrix(my_lattice%dim_mode)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! allocate the pointers for the B-field and the energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i_nim=1,nim
  do i=1,N_cell
    all_mode_path(i,i_nim)%w=>path(:,i,i_nim)
  enddo
enddo

dim_mode_mag=size(magnetic_mode_path(1,1)%w)
allocate(vel(dim_mode_mag,N_cell,nim),fxyz1(dim_mode_mag,N_cell,nim),fxyz2(dim_mode_mag,N_cell,nim),ax(dim_mode_mag,N_cell,nim),coo(dim_mode_mag,N_cell,nim),ang(N_cell,nim))
allocate(tau_i(dim_mode_mag,N_cell),tau(dim_mode_mag,N_cell))
vel=0.0d0
fxyz1=0.0d0
fxyz2=0.0d0
ax=0.0d0
coo=0.0d0
ang=0.0d0
tau_i=0.0d0
tau=0.0d0
ci=1

fchk=1d0+ftol
itr=1
dim_mode=my_lattice%dim_mode
allocate(ftmp(dim_mode),veltmp(dim_mode_mag))
ftmp=0.0d0
veltmp=0.0d0

call the_path(nim,magnetic_mode_path,pathlen)

do i_nim=1,nim

   u(i_nim) = 0d0
   do iomp=1,N_cell

      coo(:,iomp,i_nim) = magnetic_mode_path(iomp,i_nim)%w

      call calculate_Beff(ftmp,iomp,all_mode_path(:,i_nim))

      call update_B(magnetic_mode_path(iomp,i_nim)%w,0.3d0,ftmp(1:3))

      call project_force(ftmp(1:3),magnetic_mode_path(iomp,i_nim)%w,fxyz1(:,iomp,i_nim))

      call local_energy(E_int,iomp,all_mode_path(:,i_nim))

      u(i_nim) = u(i_nim) + E_int

   enddo
   if (u(i_nim).gt.u(ci)) then
      ci = i_nim
   end if
enddo

write(6,'(a,2x,I4)') 'first guess image saddle point', ci


if (en_zero=='I') then
      u0 = u(1)
elseif (en_zero=='F') then
      u0 = u(nim)
else
      u0 = 0d0
end if

!   !fname_work = 'triven_'//trim(adjustl(simid))
!      !open(777,file = fname_work, access = 'sequential', action = 'write',status = 'replace')

! Calculation of the tangent trajectory of constant energy
call tang(nim,1,coo,tau_i)
fpp = 0d0
fpp(1)=sum( fxyz1(:,:,1) * tau_i )

! calculation of the tangent trajectory for all the images
do i_nim=2,nim-1
   call tang(i_nim,coo,u,tau)

   call tang(nim,i_nim,coo,tau_i)

   fp = 0d0
   fpp(i_nim) = 0d0

   fpp(i_nim)=sum( fxyz1(:,:,i_nim) *tau_i )
   fp = sum( fxyz1(:,:,i_nim) * tau )


   if (i_nim.eq.ci) then
     do iomp=1,N_cell
       fxyz1(:,iomp,i_nim) = fxyz1(:,iomp,i_nim) - 2.0d0*tau(:,iomp)*fp
     end do
   else
     do iomp=1,N_cell
       fxyz1(:,iomp,i_nim) = fxyz1(:,iomp,i_nim) - tau(:,iomp)*fp + kappa*tau(:,iomp)*(pathlen(i_nim+1)+pathlen(i_nim-1)-2d0*pathlen(i_nim))
     end do
   endif

end do

call tang(nim,nim,coo,tau_i)
fpp(nim) = 0d0

fpp(nim) = sum( fxyz1(:,:,nim)*tau_i )


fchk=0d0
imax = 1
fchk=maxval( abs(fxyz1) )
maximum=maxloc(fxyz1)
imax=maximum(3)


itr=1

open(99,file = 'force_mep.txt',access = 'sequential',action='write',status='replace')
close(99)

open(99,file = 'force_mep.txt', access = 'sequential', action = 'write',status = 'old',position = 'append')
write(99,'(i12,a,es16.8E3,2(a,i3))',advance = 'no') itr,'   ',fchk,'   ',imax,'   ',ci

close(99)
call write_en(nim,pathlen,u-u0,-fpp,pathlen(nim),'en_path.in',do_norm_rx)

call write_path(path)


!======================MAIN LOOP============================
!
! This part is the integrator
!
!
write(6,'(/a)') 'Main loop'
write(6,'(2(a,I8))') 'maximum iteration  ', itrmax, ' iteration number', itr
write(6,'(2(a,f16.12)/)') 'distance to tolerance ', fchk, ' tolerance', ftol

do while ((fchk.gt.ftol).and.(itr.le.itrmax))
   ci=1
   do i_nim=2,nim-1
      u(i_nim) = 0d0
      do iomp=1,N_cell

         !print *,'i)nim:',i_nim
         ax(:,iomp,i_nim) = rotation_axis(coo(:,iomp,i_nim),fxyz1(:,iomp,i_nim))

         coo(:,iomp,i_nim) = magnetic_mode_path(iomp,i_nim)%w+vel(:,iomp,i_nim)*dt+0.5d0*fxyz1(:,iomp,i_nim)/mass*dt**2

         norm_local=norm(coo(:,iomp,i_nim))
         coo(:,iomp,i_nim)=coo(:,iomp,i_nim)/norm_local
         ang(iomp,i_nim) = calc_ang(coo(:,iomp,i_nim),magnetic_mode_path(iomp,i_nim)%w)
         magnetic_mode_path(iomp,i_nim)%w = coo(:,iomp,i_nim)
         !print *,'ang:',ang(iomp,i_nim)
      end do
   end do


   do i_nim=2,nim-1
      u(i_nim) = 0d0

      do iomp=1,N_cell

         call calculate_Beff(ftmp,iomp,all_mode_path(:,i_nim))

         call update_B(magnetic_mode_path(iomp,i_nim)%w,0.3d0,ftmp(1:3))

         call project_force(ftmp(1:3),magnetic_mode_path(iomp,i_nim)%w,fxyz1(:,iomp,i_nim))

         call local_energy(E_int,iomp,all_mode_path(:,i_nim))

         u(i_nim) = u(i_nim) + E_int

      enddo

      if (u(i_nim).gt.u(ci)) then
         ci = i_nim
      end if
   end do


   call the_path(nim,magnetic_mode_path,pathlen)

   !print *,'pathlen:',pathlen


   do i_nim=2,nim-1
     call tang(i_nim,coo,u,tau)
     call tang(nim,i_nim,coo,tau_i)

     fp = 0d0
     fpp(i_nim) = 0d0

     fpp(i_nim)=sum( fxyz2(:,:,i_nim) * tau_i )
     fp = sum( fxyz2(:,:,i_nim) * tau )

     if (i_nim == ci) then
       do iomp=1,N_cell
         fxyz2(:,iomp,i_nim) = fxyz2(:,iomp,i_nim) - 2.0d0*tau(:,iomp)*fp
       enddo
     else
       do iomp=1,N_cell
         fxyz2(:,iomp,i_nim) = fxyz2(:,iomp,i_nim) - tau(:,iomp)*fp + kappa*tau(:,iomp)*(pathlen(i_nim+1)+pathlen(i_nim-1)-2d0*pathlen(i_nim))
       enddo
     endif
   enddo
   !stop

   do i_nim=2,nim-1
      do iomp=1,N_cell
         call rotate(vel(:,iomp,i_nim),ax(:,iomp,i_nim),ang(iomp,i_nim),veltmp)
         call rotate(fxyz1(:,iomp,i_nim),ax(:,iomp,i_nim),ang(iomp,i_nim),ftmp(1:3))

         vel(:,iomp,i_nim) = veltmp + 0.5d0*(ftmp(1:3)+fxyz2(:,iomp,i_nim))/mass*dt

      enddo
   enddo


   fv = 0d0
   fd = 0d0

   fv= sum( vel * fxyz2 )
   fd= sum( fxyz2**2 )

   if (fv<0d0) then
     vel = 0d0
   else
     vel=fxyz2*fv/fd
   end if




   imax=1
   fchk=0d0
   fxyz1 = fxyz2

   fchk=maxval( abs(fxyz1) )
   maximum=maxloc(fxyz1)
   imax=maximum(3)



   call tang(nim,1,coo,tau_i)
   fpp(1) = sum( fxyz1(:,:,1)*tau_i )


   call tang(nim,nim,coo,tau_i)
   fpp(nim) = sum( fxyz1(:,:,nim)*tau_i )


   itr=itr+1

   if (mod(itr,every).eq.0) then
      call prn_gneb_progress(itr, itrmax, fchk, imax,'Y',ci)

      open(99,file = 'force_mep.txt', access = 'sequential', action = 'write',status = 'old',position = 'append')
      write(99,'(i12,a,es16.8E3,2(a,i3))',advance = 'no') itr,'   ',fchk,'   ',imax, '   ',ci
      close(99)

      call write_en(nim,pathlen,u-u0,-fpp,pathlen(nim),'en_path.out',do_norm_rx)
      call write_path(path)


   end if
end do

if (itr>itrmax) then
   write(6,'(a)') 'WARNING: exceeded maximum iterations in GNEB'
end if



ene = u-u0
dene = -fpp
rx = pathlen


end subroutine find_path_ci





!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Small routines to do the dirty work
!
!
!
   
!> Calculate axes of rotation based on magnetic configuration and effective fields
subroutine calc_axis_all(m_in,f_in,ax)
use m_rotation, only : rotation_axis
implicit none
real(kind=8), intent(in) :: m_in(:,:),f_in(:,:)
real(kind=8), intent(out) :: ax(:,:)
! internal variables
real(kind=8) :: x(3),y(3)
integer :: i,N_cell

N_cell=size(m_in)
do i=1,N_cell
   x = m_in(:,i)
   y = f_in(:,i)
            
   ax(:,i) = rotation_axis(x,y)

end do
   
end subroutine calc_axis_all
   



subroutine find_SP(nim,rx,c,imax,drx0)
implicit none
  
integer, intent(in) :: nim                   !< Number of images
real(kind=8), intent(in) :: rx(nim)         !< Reaction coordinate
real(kind=8), intent(in) :: c(4,nim)        !< Coefficients of the piecewise Hermite polynomials
integer, intent(out) :: imax
real(kind=8),intent(out) :: drx0
! internal variables
real(kind=8) :: eps=epsilon(drx0),d,dl,x1,x2,q
integer :: i, ci
      
ci = 1
do i=1,nim
   if (c(1,i)>c(1,ci)) then
      ci = i
   end if
end do
      
!print *, 'ci:',ci
      
if ((abs(c(2,ci))<eps).or.((ci==1).and.(c(2,ci)<0d0)).or.((ci==nim).and.(c(2,ci)>0d0))) then
   drx0 = 0d0
   imax = ci
else
   if (c(2,ci)>0d0) then
      imax = ci
   elseif (c(2,ci)<0d0) then
      imax = ci-1
   end if
   dl = rx(imax+1)-rx(imax)
   d = c(3,imax)*c(3,imax)-3d0*c(4,imax)*c(2,imax)
   if (d<0d0) then
      write(6,*) 'Energy maximum has not been found!'
   else
      q = -(c(3,imax)+dsign(1d0,c(3,imax))*sqrt(d))
      x1 = q/(c(4,imax)*3d0)
      x2 = c(2,imax)/q
      print *,'x1,x2:',x1,x2
      if ((x1>0d0).and.(x1<dl)) then
         drx0 = x1
      elseif ((x2>0d0).and.(x2<dl)) then
         drx0 = x2
      else
         write(6,*) 'Energy maximum has not been found!'
      end if
   end if
end if
  
end subroutine find_SP
   
   
subroutine find_SP_conf(spini,spinf,l1,l2,lsp,spinsp)
implicit none
real(kind=8), intent(in) :: spini(:,:)
real(kind=8), intent(in) :: spinf(:,:)
real(kind=8), intent(in) :: l1,l2,lsp
real(kind=8), intent(inout) :: spinsp(:,:)
! internal variables
integer :: i,shape_spin(2)
real(kind=8) :: ax(3)
      
ax = 0d0
shape_spin=shape(spini)

do i=1,shape_spin(2)
   spinsp(:,i) = spini(:,i)
   call find_SP_conf_one(spini(:,i),spinf(:,i),l1,l2,lsp,ax,spinsp(:,i))
end do

end subroutine find_SP_conf


subroutine find_SP_conf_one(ni,nf,l1,l2,lsp,ax,nsp)
use m_get_random
use m_vector, only : calc_ang,norm
use m_constants, only : pi
use m_rotation, only : rotation_axis
implicit none
real(kind=8), intent(in) :: ni(3),nf(3),l1,l2,lsp
real(kind=8), intent(inout) :: ax(3)
real(kind=8), intent(out) :: nsp(3)
! internal variables
real(kind=8) :: dl,theta,angle,tmp,vec(3),eps=epsilon(angle),pr,norm_local
integer :: i,j
      
      
      
dl = lsp-l1
      
angle = calc_ang(ni,nf)
if (angle<eps) then
   tmp=0d0
   do i = 1,3
      vec(i) = (nf(i) - ni(i))*dl/(l2-l1)
      nsp(i) = ni(i) + vec(i)
      tmp = tmp + nsp(i)*nsp(i)
   end do
   tmp = dsqrt(tmp)
   nsp(:) = nsp(:)/tmp
               
elseif (dabs(angle-pi(1.0d0))<eps) then
!write(*,*) 'i am here!'
         
   pr = 0.0d0
   do j=1,3
      pr = pr + ax(j)*ni(j)
   end do
   ax(:) = ax(:) - pr*ni(:)
   tmp = norm(ax)
         
   do while (tmp<eps)
      pr = 0d0
      do j=1,3
         ax(j) = dsign(1d0,2d0*get_rand_classic()-1d0)*(get_rand_classic()+1d0)
         pr = pr + ax(j)*ni(j)
      end do
      ax(:) = ax(:) - pr*ni(:)
      tmp = norm(ax)
   end do
   norm_local=norm(ax)
   ax=ax/norm_local
         
   theta = pi(1.0d0)*dl/(l2-l1)
         
   nsp(1) = ni(1)*dcos(theta) + dsin(theta)*(ax(2)*ni(3)-ax(3)*ni(2))
   nsp(2) = ni(2)*dcos(theta) - dsin(theta)*(ax(1)*ni(3)-ax(3)*ni(1))
   nsp(3) = ni(3)*dcos(theta) + dsin(theta)*(ax(1)*ni(2)-ax(2)*ni(1))
   norm_local=norm(nsp)
   nsp=nsp/norm_local
         
   
else
   ax = rotation_axis(ni,nf)

   theta = angle*dl/(l2-l1)
    
   nsp(1) = ni(1)*dcos(theta) + dsin(theta)*(ax(2)*ni(3)-ax(3)*ni(2))
   nsp(2) = ni(2)*dcos(theta) - dsin(theta)*(ax(1)*ni(3)-ax(3)*ni(1))
   nsp(3) = ni(3)*dcos(theta) + dsin(theta)*(ax(1)*ni(2)-ax(2)*ni(1))
   norm_local=norm(nsp)
   nsp=nsp/norm_local
        
end if
   
end subroutine find_SP_conf_one
   

end module m_gneb_utils
