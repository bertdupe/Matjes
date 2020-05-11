module m_gneb_utils
use m_vector, only : calc_ang,norm
use m_rotation, only : rotate
use m_path
use m_write_spin
use m_local_energy
use m_createspinfile
use m_basic_types, only : vec_point, vec
use m_derived_types, only : io_parameter,lattice
use m_gneb_parameters, only : do_norm_rx,en_zero
use m_eval_Beff
use m_rotation, only : rotation_axis
use m_operator_pointer_utils
use m_energyfield
use m_spline
use m_io_gneb
use m_projection
use m_tangent
use m_local_energy

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
type(vec), allocatable :: vel(:,:),fxyz1(:,:),fxyz2(:,:),ax(:,:),tau_i(:),tau(:),coo(:,:)
real(kind=8), allocatable :: ang(:,:)
real(kind=8) :: fchk,veltmp(3),u0,u(nim),fpp(nim)
real(kind=8), allocatable :: ftmp(:)
real(kind=8) :: pathlen(nim)
real(kind=8) :: fv,fd,fp,E_int,norm_local
integer :: i,imax,dim_mode
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


allocate(vel(N_cell,nim),fxyz1(N_cell,nim),fxyz2(N_cell,nim),ax(N_cell,nim),coo(N_cell,nim),ang(N_cell,nim))
allocate(tau_i(N_cell),tau(N_cell))

! NOT DONE:
! one has to carefully attribute the different pointers

fchk=1d0+ftol
itr=1
dim_mode=my_lattice%dim_mode
allocate(ftmp(dim_mode))
ftmp=0.0d0

call the_path(nim,magnetic_mode_path,pathlen)
   
do i_nim=1,nim

   u(i_nim) = 0d0
   do iomp=1,N_cell

      coo(iomp,i_nim)%w = magnetic_mode_path(iomp,i_nim)%w

      call calculate_Beff(ftmp,iomp,all_mode_path(:,i_nim))

      call project_force(ftmp(1:3),magnetic_mode_path(iomp,i_nim)%w,fxyz1(iomp,i_nim)%w)

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
fpp(1) = 0d0
do iomp=1,N_cell
   fpp(1) = fpp(1) + dot_product(fxyz1(iomp,1)%w,tau_i(iomp)%w)
end do
            
      
      

! calculation of the tangent trajectory for all the images
do i_nim=2,nim-1
   call tang(i_nim,coo,u,tau)

   call tang(nim,i_nim,coo,tau_i)
         
   fp = 0d0
   fpp(i_nim) = 0d0

   do iomp=1,N_cell
      fpp(i_nim) = fpp(i_nim) + dot_product(fxyz1(iomp,i_nim)%w,tau_i(iomp)%w)
      fp = fp + dot_product(fxyz1(iomp,i_nim)%w,tau(iomp)%w)
   end do

         
   do iomp=1,N_cell
     fxyz1(iomp,i_nim)%w = fxyz1(iomp,i_nim)%w - tau(iomp)%w*fp + kappa*tau(iomp)%w*(pathlen(i_nim+1)+pathlen(i_nim-1)-2d0*pathlen(i_nim))
   end do
end do
      
call tang(nim,nim,coo,tau_i)
fpp(nim) = 0d0

do iomp=1,N_cell
   fpp(nim) = fpp(nim) + dot_product(fxyz1(iomp,nim)%w,tau_i(iomp)%w)
end do

      

      
fchk=0d0
imax = 1
do i_nim=1,nim
   do iomp=1,N_cell
      do i=1,3
         if (dabs(fxyz1(iomp,i_nim)%w(i))>fchk) then
            fchk = dabs(fxyz1(iomp,i_nim)%w(i))
            imax = i_nim
         end if
      end do
   end do
end do
                  
      
itr=1
      
open(99,file = 'force_mep.txt',access = 'sequential',action='write',status='replace')
close(99)
      
open(99,file = 'force_mep.txt', access = 'sequential', action = 'write',status = 'old',position = 'append')
write(99,'(i12,a,es16.8E3,a,i3)',advance = 'no') itr,'   ',fchk,'   ',imax
      
close(99)
do i=1,nim
   ene(i) = (u(i)-u0)
   dene(i) = -fpp(i)
   rx(i) = pathlen(i)
end do
call write_en(nim,rx,ene,dene,rx(nim),'en_path.in',do_norm_rx)
      
call write_path(path)
      
      
!======================MAIN LOOP============================
!
! This part is the integrator
!
!

do while ((fchk>ftol).and.(itr<=itrmax))
   do i_nim=2,nim-1
      u(i_nim) = 0d0
      do iomp=1,N_cell
            
         !print *,'i)nim:',i_nim
         ax(iomp,i_nim)%w = rotation_axis(coo(iomp,i_nim)%w,fxyz1(iomp,i_nim)%w)
                        
         coo(iomp,i_nim)%w = magnetic_mode_path(iomp,i_nim)%w+vel(iomp,i_nim)%w*dt+0.5d0*fxyz1(iomp,i_nim)%w/mass*dt*dt

         norm_local=norm(coo(iomp,i_nim)%w)
         coo(iomp,i_nim)%w=coo(iomp,i_nim)%w/norm_local
         ang(iomp,i_nim) = calc_ang(coo(iomp,i_nim)%w,magnetic_mode_path(iomp,i_nim)%w)
         magnetic_mode_path(iomp,i_nim)%w = coo(iomp,i_nim)%w
         !print *,'ang:',ang(i_x,i_y,i_z,i_m)
      end do
   end do
            
         
   do i_nim=2,nim-1
      u(i_nim) = 0d0

      do iomp=1,N_cell

         call calculate_Beff(ftmp,iomp,all_mode_path(:,i_nim))

         call project_force(ftmp(1:3),magnetic_mode_path(iomp,i_nim)%w,fxyz1(iomp,i_nim)%w)

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
         
      do iomp=1,N_cell
         fpp(i_nim) = fpp(i_nim) + dot_product(fxyz2(iomp,i_nim)%w,tau_i(iomp)%w)
         fp = fp + dot_product(fxyz2(iomp,i_nim)%w,tau(iomp)%w)
      end do

      do iomp=1,N_cell
         fxyz2(iomp,i_nim)%w = fxyz2(iomp,i_nim)%w - tau(iomp)%w*fp + kappa*tau(iomp)%w*(pathlen(i_nim+1)+pathlen(i_nim-1)-2d0*pathlen(i_nim))
      enddo
   enddo
   !stop

   do i_nim=2,nim-1
      do iomp=1,N_cell
         call rotate(vel(iomp,i_nim)%w,ax(iomp,i_nim)%w,ang(iomp,i_nim),veltmp)
         call rotate(fxyz1(iomp,i_nim)%w,ax(iomp,i_nim)%w,ang(iomp,i_nim),ftmp)

         vel(iomp,i_nim)%w = veltmp + 0.5d0*(ftmp(1:3)+fxyz2(iomp,i_nim)%w)/mass*dt

      enddo
   enddo
         
         
   fv = 0d0
   fd = 0d0
   do i_nim=1,nim
      do iomp=1,N_cell
                        
         fv = fv + dot_product(vel(iomp,i_nim)%w,fxyz2(iomp,i_nim)%w)
         fd = fd + dot_product(fxyz2(iomp,i_nim)%w,fxyz2(iomp,i_nim)%w)

      end do
   end do

         
   if (fv<0d0) then
      do i_nim=1,nim
         do iomp=1,N_cell
            vel(iomp,i_nim)%w = 0d0
         end do
      end do
   else
      do i_nim=1,nim
         do iomp=1,N_cell
            vel(iomp,i_nim)%w = fxyz2(iomp,i_nim)%w*fv/fd
         end do
      end do
   end if
         
         
         
         
   imax=1
   fchk=0d0

   do i_nim=2,nim-1
      do iomp=1,N_cell
         do i=1,3
            fxyz1(iomp,i_nim)%w(i) = fxyz2(iomp,i_nim)%w(i)
            if (dabs(fxyz1(iomp,i_nim)%w(i))>fchk) then
               fchk = dabs(fxyz1(iomp,i_nim)%w(i))
               imax = i_nim
            end if
         end do
      end do
   end do
         
         
         
         
   call tang(nim,1,coo,tau_i)
   fpp(1) = 0d0
   do iomp=1,N_cell
      fpp(1) = fpp(1) + dot_product(fxyz1(iomp,1)%w,tau_i(iomp)%w)
   end do
         
         
   call tang(nim,nim,coo,tau_i)
   fpp(nim) = 0d0
   do iomp=1,N_cell
      fpp(nim) = fpp(nim) + dot_product(fxyz1(iomp,nim)%w,tau_i(iomp)%w)
   end do
         
         
         
   itr=itr+1

   if (mod(itr,every).eq.0) then
      call prn_gneb_progress(itr, itrmax, fchk, imax,'N',0)
      
      open(99,file = 'force_mep.txt', access = 'sequential', action = 'write',status = 'old',position = 'append')
      write(99,'(i12,a,es16.8E3,a,i3)',advance = 'no') itr,'   ',fchk,'   ',imax
      close(99)
            
      do i=1,nim
         ene(i) = (u(i)-u0)
         dene(i) = -fpp(i)
         rx(i) = pathlen(i)
      end do
      call write_en(nim,rx,ene,dene,rx(nim),'en_path.out',do_norm_rx)
      call write_path(path)
            
         
   end if
end do

if (itr>itrmax) then
   write(*,*) 'WARNING: exceeded maximum iterations in GNEB'
end if
   
      
      
      
do i=1,nim
   ene(i) = (u(i)-u0)
   dene(i) = -fpp(i)
   rx(i) = pathlen(i)
end do
   

end subroutine find_path







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
type(vec), intent(in) :: m_in(:),f_in(:)
type(vec), intent(out) :: ax(:)
! internal variables
real(kind=8) :: x(3),y(3)
integer :: i,N_cell

N_cell=size(m_in)
do i=1,N_cell
   x = m_in(i)%w
   y = f_in(i)%w
            
   ax(i)%w = rotation_axis(x,y)

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
use mtprng
use m_vector, only : calc_ang,norm
use m_constants, only : pi
use m_rotation, only : rotation_axis
implicit none
real(kind=8), intent(in) :: ni(3),nf(3),l1,l2,lsp
real(kind=8), intent(inout) :: ax(3)
real(kind=8), intent(out) :: nsp(3)
! internal variables
type(mtprng_state) :: state
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
         ax(j) = dsign(1d0,2d0*mtprng_rand_real1(state)-1d0)*(mtprng_rand_real1(state)+1d0)
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
