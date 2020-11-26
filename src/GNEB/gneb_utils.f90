module m_gneb_utils
use m_vector, only : calc_ang,norm,normalize,project
use m_path, only: the_path
use m_basic_types, only : vec_point, vec
use m_derived_types, only : io_parameter,lattice
use m_gneb_parameters, only : do_norm_rx,en_zero
use m_rotation, only : rotation_axis,rotate
use m_io_gneb, only: write_path,prn_gneb_progress,write_en
use m_tangent, only: tang
use m_input_types,only: GNEB_input
use m_type_lattice,only: lattice
use m_H_public, only: t_H,energy_all
implicit none
private
public :: find_path,find_SP,find_SP_conf

contains

subroutine find_path(nim,N_cell,ftol,rx,ene,dene,images,io_simu,io_gneb,Hams,ci_out)
    use m_precision, only: truncate
    use m_Beff_H,only: get_B
    type(io_parameter), intent(in)  :: io_simu
    real(8),intent(in)              :: ftol
    integer, intent(in)             :: nim
    type(GNEB_input)                :: io_gneb
    integer, intent(in)             :: N_cell
    type(lattice), intent(inout)    :: images(:)
    real(8), intent(out)            :: rx(nim) !< Reaction coordinate
    real(8), intent(out)            :: ene(nim) !< Energy of the images
    real(8), intent(out)            :: dene(nim) !< Derivative of the energy with respect to reaction coordinate
    class(t_H),intent(in)           :: Hams(:)
    integer,intent(out),optional    :: ci_out
    ! internal
    integer     :: N_mag
    integer(8)  :: itr
    real(8), allocatable, dimension(:,:,:) :: vel,ax
    real(8), allocatable, dimension(:,:) :: tau_i,tau,ang
    real(8), allocatable, dimension(:) :: ftmp,veltmp
    real(8)             :: fchk,energy_ref,energy(nim),fpp(nim)
    real(8)             :: pathlen(nim)
    real(8)             :: fv,fd,fp
    integer             :: imax,dim_mode_mag,maximum(2)
    real(8),target,allocatable      :: force1(:,:)
    real(8),target,allocatable      :: force2(:,:)
    integer             :: im
    real(8),pointer,contiguous     :: force1_3(:,:),force2_3(:,:)
    real(8),pointer,contiguous     :: force1_mode(:,:),force2_mode(:,:)
    real(8),pointer,contiguous     :: force_all_3(:,:,:)
    real(8),pointer,contiguous     :: M3(:,:)
    real(8),allocatable,target :: magarr_tmp(:,:)
    integer             :: ci
    real(8), allocatable, dimension(:,:) :: vel_tmp_arr
    logical             :: converged

    ci=0
    converged=.false.
    if(present(ci_out)) ci=1
    energy_ref = 0.0d0

    dim_mode_mag=images(1)%M%dim_mode
    allocate(ftmp(dim_mode_mag),veltmp(dim_mode_mag),source=0.0d0)
    N_mag=dim_mode_mag/3
    allocate(vel(dim_mode_mag,N_cell,nim),ax(dim_mode_mag,N_cell,nim),ang(N_cell,nim),source=0.0d0)
    allocate(vel_tmp_arr(dim_mode_mag,N_cell),source=0.0d0)
    allocate(tau_i(dim_mode_mag,N_cell),tau(dim_mode_mag,N_cell),source=0.0d0)
    allocate(force1(dim_mode_mag*N_cell,nim),force2(dim_mode_mag*N_cell,nim),source=0.d0)
    allocate(magarr_tmp(3,N_mag*N_cell),source=0.0d0)
    
    fchk=1d0+ftol
   
    if(dim_mode_mag/=3) ERROR STOP "THE_PATH AND PROBABLY MORE PARTS OF THIS ROUTINE MIGHT NOT MAKE SENSE FOR MORE THAN ONE NMAG"
    call the_path(images,pathlen)

    energy = 0d0
    do im=1,nim
        force1_3(1:3,1:N_mag*N_cell)=>force1(:,im)
        M3(1:3,1:N_mag*N_cell)=>images(im)%M%modes

        energy(im)=energy_all(Hams,images(im))
        Call get_B(Hams,images(im),force1(:,im))
        !Call normalize(force1_3) !temporary taken out, think about normalizations
        Call project(force1_3,M3)
    enddo
    if(present(ci_out))then
        ci=maxloc(energy,dim=1)  !change ci if an image is higher in energy
        write(6,'(a,2x,I4)') 'first guess image saddle point', ci
    endif
    
    if (io_gneb%en_zero=='I') energy_ref = energy(1)
    if (io_gneb%en_zero=='F') energy_ref = energy(nim)
    
    fpp = 0d0
    ! Calculation of the tangent trajectory of constant energy
    call tang(1,images,tau_i)
    force1_mode(1:dim_mode_mag,1:N_cell)=>force1(:,1)
    fpp(1)=sum(force1_mode*tau_i)

    call tang(nim,images,tau_i)
    force1_mode(1:dim_mode_mag,1:N_cell)=>force1(:,nim)
    fpp(nim)=sum(force1_mode*tau_i)
    
    ! calculation of the tangent trajectory for all the images
    do im=2,nim-1
        force1_mode(1:dim_mode_mag,1:N_cell)=>force1(:,im)
        call tang(im,images,energy,tau)
        call tang(im,images  ,tau_i)
        fpp(im)=sum(force1_mode*tau_i)
        fp=sum(force1_mode*tau)
        if(im/=ci)then
            force1_mode = force1_mode - tau*fp + io_gneb%spring*tau*(pathlen(im+1)+pathlen(im-1)-2d0*pathlen(im))
        else
            force1_mode = force1_mode - 2.0d0*tau*fp 
        endif
    end do
    
    fchk=maxval(abs(force1))
    maximum=maxloc(force1)
    imax=maximum(2)
    itr=1
          
    open(99,file = 'force_mep.txt', access = 'sequential', action = 'write')
    if(present(ci_out))then
        write(99,'(i12,a,E20.12E3,2(a,i3))',advance = 'no') itr,'   ',fchk,'   ',imax,'   ',ci
    else
        write(99,'(i12,a,E20.12E3,a,i3)',advance = 'no') itr,'   ',fchk,'   ',imax
    endif
    close(99)
    call write_en(nim,pathlen,energy-energy_ref,-fpp,pathlen(nim),'en_path.in',io_gneb%do_norm_rx)
    call write_path(images)
          
    !======================MAIN LOOP============================
    !
    ! This part is the integrator
    !
    !
    write(6,'(/a)') 'Main loop'
    write(6,'(2(a,I22))') 'maximum iteration  ', io_gneb%itrmax, ' iteration number', itr
    write(6,'(2(a,E20.12E3)/)') 'distance to tolerance ', fchk, ' tolerance', ftol

    do itr=1,io_gneb%itrmax
        do im=2,nim-1
            M3(1:3,1:N_mag*N_cell)=>images(im)%M%modes
            force1_3(1:3,1:N_mag*N_cell)=>force1(:,im)
            ax(:,:,im)= rotation_axis(M3,force1_3)
            magarr_tmp=M3+vel(:,:,im)*io_gneb%dt+0.5d0*force1_3/io_gneb%mass*io_gneb%dt**2
            Call normalize(magarr_tmp,1.0d-30)
            ang(:,im)=calc_ang(magarr_tmp,M3)
            M3=magarr_tmp
        enddo

        do im=2,nim-1
            force2_3(1:3,1:N_mag*N_cell)=>force2(:,im)
            M3(1:3,1:N_mag*N_cell)=>images(im)%M%modes

            energy(im)=energy_all(Hams,images(im)) !total energy at each image
            Call get_B(Hams,images(im),force2(:,im)) !get effective field (eV)
            !Call normalize(force2_3) !temporary taken out, think about normalizations
            Call project(force2_3,M3)
        enddo
        if(present(ci_out)) ci=maxloc(energy,dim=1)  !change ci if an image is higher in energy

        call the_path(images,pathlen) !polygeodesic length of the path in computed and stored in pathlen (rad?)

        do im=2,nim-1
            force2_mode(1:dim_mode_mag,1:N_cell)=>force2(:,im)
            call tang(im,images,energy,tau) ! gets tau: normalized tangent to the path at image i_nim, projected onto tangent space
            call tang(im,images  ,tau_i)    ! outputs tau_i at image i_nim?
            fpp(im)=sum(force2_mode*tau_i)  !force is Beff projected on tangent space
            fp=sum(force2_mode*tau)
            if(im/=ci)then
                force2_mode = force2_mode - tau*fp + io_gneb%spring*tau*(pathlen(im+1)+pathlen(im-1)-2d0*pathlen(im))
            else
                force2_mode = force2_mode - 2.0d0*tau*fp 
            endif
        end do
    
        do im=2,nim-1
            force1_3(1:3,1:N_mag*N_cell)=>force1(:,im)
            force2_3(1:3,1:N_mag*N_cell)=>force2(:,im)

            call rotate(vel(:,:,im),ax(:,:,im),ang(:,im),vel_tmp_arr)
            vel(:,:,im)=vel_tmp_arr
            call rotate(force1_3,ax(:,:,im),ang(:,im),vel_tmp_arr)
            vel(:,:,im)=vel(:,:,im)+0.5d0*(vel_tmp_arr+force2_3)/io_gneb%mass*io_gneb%dt
        enddo

        force_all_3(1:dim_mode_mag,1:N_cell,1:nim)=>force2
        fv= sum(vel * force_all_3 )/real(nim,8)
        fd= sum(force2**2)/real(nim,8)
          
        if (fv<0d0) then
          vel = 0d0
        else
          vel=force_all_3*fv/fd
        end if

        !truncate to prevent underflows
        Call truncate(vel)
        Call truncate(force2)

        force1=force2
        fchk=maxval(abs(force1))
        if (mod(itr,io_gneb%every).eq.0) then
            maximum=maxloc(force1)
            imax=maximum(2)
              
            call tang(1,images,tau_i)
            force1_mode(1:dim_mode_mag,1:N_cell)=>force1(:,1)
            fpp(1)=sum(force1_mode*tau_i)

            call tang(nim,images,tau_i)
            force1_mode(1:dim_mode_mag,1:N_cell)=>force1(:,nim)
            fpp(nim)=sum(force1_mode*tau_i)

            if(present(ci_out))then
                call prn_gneb_progress(itr, io_gneb%itrmax, fchk, imax,'Y',ci)
            else
                call prn_gneb_progress(itr, io_gneb%itrmax, fchk, imax,'N',0)
            endif
            
            open(99,file = 'force_mep.txt', access = 'sequential', action = 'write',status = 'old',position = 'append')
            if(present(ci_out))then
                write(99,'(i12,a,E20.12E3,2(a,i3))',advance = 'no') itr,'   ',fchk,'   ',imax,'   ',ci
            else
                write(99,'(i12,a,E20.12E3,a,i3)',advance = 'no') itr,'   ',fchk,'   ',imax
            endif
            close(99)
                  
            call write_en(nim,pathlen,energy-energy_ref,-fpp,pathlen(nim),'en_path.out',io_gneb%do_norm_rx)
            call write_path(images)
        end if

        converged=fchk.lt.ftol
        if(converged)then
            write(6,'(//A//)') "GNEB convergence threshold reached"
            exit
        endif
    end do
    
    if (.not.converged) then
       write(6,'(///a///)') 'WARNING: exceeded maximum iterations in GNEB'
    end if

    call tang(1,images,tau_i)
    force1_mode(1:dim_mode_mag,1:N_cell)=>force1(:,1)
    fpp(1)=sum(force1_mode*tau_i)

    call tang(nim,images,tau_i)
    force1_mode(1:dim_mode_mag,1:N_cell)=>force1(:,nim)
    fpp(nim)=sum(force1_mode*tau_i)
    
    ene = energy-energy_ref
    dene = -fpp
    rx = pathlen
    if(present(ci_out)) ci_out=ci
end subroutine find_path

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
   
   
subroutine find_SP_conf(lat1,lat2,l1,l2,lsp,spinsp)
    implicit none
    type(lattice),intent(in)    :: lat1,lat2
    real(8), intent(in)         :: l1,l2,lsp
    real(8), intent(inout)      :: spinsp(:,:)
    ! internal variables
    integer :: i
    real(kind=8) :: ax(3)

    ax = 0d0
    do i=1,lat1%Ncell
       spinsp(:,i) = lat1%M%modes_v(:,i)
       call find_SP_conf_one(lat1%M%modes_v(:,i),lat2%M%modes_v(:,i),l1,l2,lsp,ax,spinsp(:,i))
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
                   
    elseif (dabs(angle-pi)<eps) then
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
             
       theta = pi*dl/(l2-l1)
             
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
