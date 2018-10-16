module m_gneb_utils
   use m_write_spin
   use m_eval_Beff
   use m_local_energy
   use m_createspinfile
   use m_derived_types

   character(len=35) :: momfile_i                             !< File name for initial magn. configuration
   character(len=35) :: momfile_f                             !< File name for final magn. configuration
   character(len=35) :: restartfile_if                        !< Name of the file with initial and final states needed to restart energy minimization calculations
   character(len=35) :: restartfile_path                      !< Name of the file with the path needed to restart GNEB calculations
   real(kind=8) :: ind_tol
   real(kind=8) :: amp_rnd                                   !< Amplitude of random perturbation of the components of magnetic moments
   real(kind=8) :: amp_rnd_path
   integer :: minalgo                                           !< Minimization algorithm (1=VPO,...)
   integer :: minitrmax                                         !< Maximum number of iterations in energy minimization procedure
   real(kind=8) :: minftol                                     !< Convergence criterion in mRy
   integer :: mintraj_step                                      !< Save configuration every 'mintraj_step' step during energy minimization
   real(kind=8) :: vpodt                                       !< VPO timestep
   real(kind=8) :: vpomass                                     !< VPO mass
   integer :: initpath                                          !< Path initialization method (1=geodesic, 2=read from file,...)
   real(kind=8) :: spring                                      !< Magnitude of the spring constant in GNEB calculations
   integer :: mepitrmax                                         !< Maximum number of iterations in MEP finding procedure
   real(kind=8) :: mepftol                                     !< Convergence criterion in the GNEB method, in mRy
   real(kind=8) :: mepftol_ci                                  !< Convergence criterion in the CI-GNEB method, in mRy
   integer :: meptraj_step                                      !< Save configuration every 'meptraj_step' step during MEP finding procedure
   character(len=1) :: do_gneb                                  !< Do GNEB calculations (Y/N)
   character(len=1) :: do_gneb_ci                               !< Do CI-GNEB calculations (Y/N)
   character(len=1) :: do_norm_rx                               !< Normalize reaction coordinate (Y/N)
   character(len=1) :: en_zero                                  !< Level of zero energy. 'I' - initial state; 'F' - final state; 'N' - 0.0
   integer :: sample_num                                         !< Number of samples in the interpolated curve
   integer :: nim
   
contains

subroutine find_path(i_DM,i_four,i_biq,i_dip,EA,nim,dt,mass,kappa,path, &
   & shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,ftol,itrmax,every,h_ext,my_lattice, &
   & rx,ene,dene)

implicit none
   type(lattice), intent(in) :: my_lattice
   logical, intent(in) :: i_DM,i_four,i_biq,i_dip
   real(kind=8), intent(in) :: EA(3)    ! easy axis
   integer, intent(in) :: nim,shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
   real(kind=8), intent(in) :: mass             !> mass of the point
   real(kind=8), intent(in) :: dt               !> timestep
   real(kind=8), intent(in) :: kappa            !> spring constant
   real(kind=8), intent(in) :: ftol            !> spring constant
   integer, intent(in) :: itrmax
   real(kind=8), intent(inout) :: path(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
   integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4)),indexNN(shape_index(1),shape_index(2))
   integer, intent(in) :: every !< Save path every 'every' step
   real(kind=8), intent(in) :: h_ext(3)   ! external magnetic field
   real(kind=8), intent(out) :: rx(nim) !< Reaction coordinate
   real(kind=8), intent(out) :: ene(nim) !< Energy of the images
   real(kind=8), intent(out) :: dene(nim) !< Derivative of the energy with respect to reaction coordinate
   integer :: i_m,i_z,i_x,i_y,i_nim,itr
   real(kind=8) :: vel(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   real(kind=8) :: fxyz1(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   real(kind=8) :: fxyz2(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   real(kind=8) :: fchk,ftmp(3),veltmp(3),u0,u(nim),fpp(nim),tau_i(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
   real(kind=8) :: tau(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
   real(kind=8) :: pathlen(nim),ax(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim),ang(shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   real(kind=8) :: coo(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   real(kind=8) :: fv,fd,fp,E_int
   integer :: i,imax
   
   fchk=1d0+ftol
   itr=1
   fxyz1 = 0d0
   fxyz2 = 0d0
   call the_path(nim,shape_spin,path,pathlen)
   
   do i_nim=1,nim

      u(i_nim) = 0d0
      do i_m=1,shape_spin(5)
         do i_z=1,shape_spin(4)
            do i_y=1,shape_spin(3)
               do i_x=1,shape_spin(2)
                  vel(:,i_x,i_y,i_z,i_m,i_nim) = 0d0
                  do i = 1,3
                     coo(i,i_x,i_y,i_z,i_m,i_nim) = path(3+i,i_x,i_y,i_z,i_m,i_nim)
                  end do
                  call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,ftmp, &
                    &  path(:,:,:,:,:,i_nim),shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext,my_lattice)
          
                  call project_force(ftmp,path(4:6,i_x,i_y,i_z,i_m,i_nim),fxyz1(:,i_x,i_y,i_z,i_m,i_nim))

                  call local_energy(E_int,i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,path(:,:,:,:,:,i_nim), &
                       & shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,my_lattice)

                  u(i_nim) = u(i_nim) + E_int
               end do
            end do
         end do
      end do
   end do

   !stop
   
   
   
             

                        
   if (en_zero=='I') then
      u0 = u(1)
   elseif (en_zero=='F') then
      u0 = u(nim)
   else
      u0 = 0d0
   end if
   
!   !fname_work = 'triven_'//trim(adjustl(simid))
!      !open(777,file = fname_work, access = 'sequential', action = 'write',status = 'replace')
      
      call tang(nim,shape_spin,1,coo,tau_i)
      fpp(1) = 0d0
      do i_m=1,shape_spin(5)
         do i_z=1,shape_spin(4)
            do i_y=1,shape_spin(3)
               do i_x=1,shape_spin(2)
                  do i=1,3
                     fpp(1) = fpp(1) + fxyz1(i,i_x,i_y,i_z,i_m,1)*tau_i(i,i_x,i_y,i_z,i_m)
                  end do
               end do
            end do
         end do
      end do
            
      
      

      
      do i_nim=2,nim-1
         call tang_spec(nim,shape_spin,i_nim,coo,u,tau)

         call tang(nim,shape_spin,i_nim,coo,tau_i)
         
         fp = 0d0
         fpp(i_nim) = 0d0
         
         do i_m=1,shape_spin(5)
            do i_z=1,shape_spin(4)
               do i_y=1,shape_spin(3)
                  do i_x=1,shape_spin(2)
                     do i=1,3
                        fpp(i_nim) = fpp(i_nim) + fxyz1(i,i_x,i_y,i_z,i_m,i_nim)*tau_i(i,i_x,i_y,i_z,i_m)
                        fp = fp + fxyz1(i,i_x,i_y,i_z,i_m,i_nim)*tau(i,i_x,i_y,i_z,i_m)
                     end do
                  end do
               end do
            end do
         end do
         
         
         do i_m=1,shape_spin(5)
            do i_z=1,shape_spin(4)
               do i_y=1,shape_spin(3)
                  do i_x=1,shape_spin(2)
                     do i=1,3
                        !print *,fxyz1(i,i_x,i_y,i_z,i_m,i_nim)
                        fxyz1(i,i_x,i_y,i_z,i_m,i_nim) = fxyz1(i,i_x,i_y,i_z,i_m,i_nim) - tau(i,i_x,i_y,i_z,i_m)*fp + kappa*tau(i,i_x,i_y,i_z,i_m)*(pathlen(i_nim+1)+pathlen(i_nim-1)-2d0*pathlen(i_nim))
                     end do
                  end do
               end do
            end do
         end do
         

      end do
      
      call tang(nim,shape_spin,nim,coo,tau_i)
      fpp(nim) = 0d0
      do i_m=1,shape_spin(5)
         do i_z=1,shape_spin(4)
            do i_y=1,shape_spin(3)
               do i_x=1,shape_spin(2)
                  do i=1,3
                     fpp(nim) = fpp(nim) + fxyz1(i,i_x,i_y,i_z,i_m,nim)*tau_i(i,i_x,i_y,i_z,i_m)
                  end do
               end do
            end do
         end do
      end do
      

      
      fchk=0d0
      imax = 1
      do i_nim=1,nim
         do i_m=1,shape_spin(5)
            do i_z=1,shape_spin(4)
               do i_y=1,shape_spin(3)
                  do i_x=1,shape_spin(2)
                     do i=1,3
                        if (dabs(fxyz1(i,i_x,i_y,i_z,i_m,i_nim))>fchk) then
                           fchk = dabs(fxyz1(i,i_x,i_y,i_z,i_m,i_nim))
                           imax = i_nim
                        end if
                     end do
                  end do
               end do
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
      
      call write_path(nim,shape_spin,path)
      
      
!======================MAIN LOOP============================
      do while ((fchk>ftol).and.(itr<=itrmax))
         do i_nim=2,nim-1
            u(i_nim) = 0d0
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
            
                        !print *,'i)nim:',i_nim
                        ax(:,i_x,i_y,i_z,i_m,i_nim) = calc_axis(coo(:,i_x,i_y,i_z,i_m,i_nim),fxyz1(:,i_x,i_y,i_z,i_m,i_nim))
                        
                        coo(:,i_x,i_y,i_z,i_m,i_nim) = path(4:6,i_x,i_y,i_z,i_m,i_nim)+vel(:,i_x,i_y,i_z,i_m,i_nim)*dt+0.5d0*fxyz1(:,i_x,i_y,i_z,i_m,i_nim)/mass*dt*dt

                        call normalize_vec(3,coo(:,i_x,i_y,i_z,i_m,i_nim))
                        ang(i_x,i_y,i_z,i_m,i_nim) = calc_ang(coo(:,i_x,i_y,i_z,i_m,i_nim),path(4:6,i_x,i_y,i_z,i_m,i_nim))
                        path(4:6,i_x,i_y,i_z,i_m,i_nim) = coo(:,i_x,i_y,i_z,i_m,i_nim)
                        !print *,'ang:',ang(i_x,i_y,i_z,i_m)
                     end do
                  end do
               end do
            end do
         end do
            
         
         do i_nim=2,nim-1
            u(i_nim) = 0d0
            
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
                        
                        call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,ftmp, &
                              &  path(:,:,:,:,:,i_nim),shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext,my_lattice)

                        call project_force(ftmp,path(4:6,i_x,i_y,i_z,i_m,i_nim),fxyz2(:,i_x,i_y,i_z,i_m,i_nim))

                        call local_energy(E_int,i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,path(:,:,:,:,:,i_nim), &
                             &shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,my_lattice)

                        u(i_nim) = u(i_nim) + E_int

                        !print *,'ene:',u(i_nim)
                     end do
                  end do
               end do
            end do
            
            
         end do
         
         
         
         call the_path(nim,shape_spin,path,pathlen)
         
         !print *,'pathlen:',pathlen
         
         do i_nim=2,nim-1
            
            
            
            
            
            call tang_spec(nim,shape_spin,i_nim,coo,u,tau)
            
            call tang(nim,shape_spin,i_nim,coo,tau_i)
         
            fp = 0d0
            fpp(i_nim) = 0d0
         
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
                        do i=1,3
                           fpp(i_nim) = fpp(i_nim) + fxyz2(i,i_x,i_y,i_z,i_m,i_nim)*tau_i(i,i_x,i_y,i_z,i_m)
                           fp = fp + fxyz2(i,i_x,i_y,i_z,i_m,i_nim)*tau(i,i_x,i_y,i_z,i_m)
                        end do
                     end do
                  end do
               end do
            end do
         
         
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
                        do i=1,3
                           fxyz2(i,i_x,i_y,i_z,i_m,i_nim) = fxyz2(i,i_x,i_y,i_z,i_m,i_nim) - tau(i,i_x,i_y,i_z,i_m)*fp + kappa*tau(i,i_x,i_y,i_z,i_m)*(pathlen(i_nim+1)+pathlen(i_nim-1)-2d0*pathlen(i_nim))
                        end do
                     end do
                  end do
               end do
            end do
         end do
         !stop
         fv = 0d0
         fd = 0d0
         do i_nim=2,nim-1
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
                        call rotate_vec(vel(:,i_x,i_y,i_z,i_m,i_nim),ax(:,i_x,i_y,i_z,i_m,i_nim),ang(i_x,i_y,i_z,i_m,i_nim),veltmp)
                        call rotate_vec(fxyz1(:,i_x,i_y,i_z,i_m,i_nim),ax(:,i_x,i_y,i_z,i_m,i_nim),ang(i_x,i_y,i_z,i_m,i_nim),ftmp)
                        !vel(:,i_x,i_y,i_z,i_m,i_nim) = vel(:,i_x,i_y,i_z,i_m,i_nim) + 0.5d0*(fxyz1(:,i_x,i_y,i_z,i_m,i_nim)+fxyz2(:,i_x,i_y,i_z,i_m,i_nim))/mass*dt
                        vel(:,i_x,i_y,i_z,i_m,i_nim) = veltmp + 0.5d0*(ftmp+fxyz2(:,i_x,i_y,i_z,i_m,i_nim))/mass*dt
                        do i=1,3
                           fv = fv + vel(i,i_x,i_y,i_z,i_m,i_nim)*fxyz2(i,i_x,i_y,i_z,i_m,i_nim)
                           fd = fd + fxyz2(i,i_x,i_y,i_z,i_m,i_nim)*fxyz2(i,i_x,i_y,i_z,i_m,i_nim)
                        end do
                     end do
                  end do
               end do
            end do
         end do
         
         
         fv = 0d0
         fd = 0d0
         do i_nim=1,nim
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
                        
                        do i=1,3
                           fv = fv + vel(i,i_x,i_y,i_z,i_m,i_nim)*fxyz2(i,i_x,i_y,i_z,i_m,i_nim)
                           fd = fd + fxyz2(i,i_x,i_y,i_z,i_m,i_nim)*fxyz2(i,i_x,i_y,i_z,i_m,i_nim)
                        end do
                     end do
                  end do
               end do
            end do
         end do
         
         if (fv<0d0) then
            do i_nim=1,nim
               do i_m=1,shape_spin(5)
                  do i_z=1,shape_spin(4)
                     do i_y=1,shape_spin(3)
                        do i_x=1,shape_spin(2)
                           do i=1,3
                              vel(i,i_x,i_y,i_z,i_m,i_nim) = 0d0
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         else
            do i_nim=1,nim
               do i_m=1,shape_spin(5)
                  do i_z=1,shape_spin(4)
                     do i_y=1,shape_spin(3)
                        do i_x=1,shape_spin(2)
                           do i=1,3
                              vel(i,i_x,i_y,i_z,i_m,i_nim) = fxyz2(i,i_x,i_y,i_z,i_m,i_nim)*fv/fd
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end if
         
         
         
         
         imax=1
         fchk=0d0

         do i_nim=2,nim-1
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
                        do i=1,3
                           fxyz1(i,i_x,i_y,i_z,i_m,i_nim) = fxyz2(i,i_x,i_y,i_z,i_m,i_nim)
                           if (dabs(fxyz1(i,i_x,i_y,i_z,i_m,i_nim))>fchk) then
                              fchk = dabs(fxyz1(i,i_x,i_y,i_z,i_m,i_nim))
                              imax = i_nim
                           end if
                        end do
                     end do
                  end do
               end do
            end do
         end do
         
         
         
         
         
         call tang(nim,shape_spin,1,coo,tau_i)
         fpp(1) = 0d0
         do i_m=1,shape_spin(5)
            do i_z=1,shape_spin(4)
               do i_y=1,shape_spin(3)
                  do i_x=1,shape_spin(2)
                     do i=1,3
                        fpp(1) = fpp(1) + fxyz1(i,i_x,i_y,i_z,i_m,1)*tau_i(i,i_x,i_y,i_z,i_m)
                     end do
                  end do
               end do
            end do
         end do
         
         
         call tang(nim,shape_spin,nim,coo,tau_i)
         fpp(nim) = 0d0
         do i_m=1,shape_spin(5)
            do i_z=1,shape_spin(4)
               do i_y=1,shape_spin(3)
                  do i_x=1,shape_spin(2)
                     do i=1,3
                        fpp(nim) = fpp(nim) + fxyz1(i,i_x,i_y,i_z,i_m,nim)*tau_i(i,i_x,i_y,i_z,i_m)
                     end do
                  end do
               end do
            end do
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
            call write_path(nim,shape_spin,path)
            
         
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



subroutine find_path_ci(i_DM,i_four,i_biq,i_dip,EA,nim,dt,mass,kappa,path, &
   & shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,ftol,itrmax,every,h_ext,my_lattice, &
   & rx,ene,dene,ci)

implicit none
   type(lattice), intent(in) :: my_lattice
   logical, intent(in) :: i_DM,i_four,i_biq,i_dip
   real(kind=8), intent(in) :: EA(3)    ! easy axis
   integer, intent(in) :: nim,shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
   real(kind=8), intent(in) :: mass             !> mass of the point
   real(kind=8), intent(in) :: dt               !> timestep
   real(kind=8), intent(in) :: kappa            !> spring constant
   real(kind=8), intent(in) :: ftol            !> spring constant
   integer, intent(in) :: itrmax
   real(kind=8), intent(inout) :: path(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
   integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4)),indexNN(shape_index(1),shape_index(2))
   integer, intent(in) :: every !< Save path every 'every' step
   real(kind=8), intent(in) :: h_ext(3)   ! external magnetic field
   real(kind=8), intent(out) :: rx(nim) !< Reaction coordinate
   real(kind=8), intent(out) :: ene(nim) !< Energy of the images
   real(kind=8), intent(out) :: dene(nim) !< Derivative of the energy with respect to reaction coordinate
   integer, intent(out) :: ci
   integer :: i_m,i_z,i_x,i_y,i_nim,itr
   real(kind=8) :: vel(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   real(kind=8) :: fxyz1(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   real(kind=8) :: fxyz2(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   real(kind=8) :: fchk,ftmp(3),veltmp(3),u0,u(nim),fpp(nim),tau_i(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
   real(kind=8) :: tau(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
   real(kind=8) :: pathlen(nim),ax(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim),ang(shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   real(kind=8) :: coo(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   real(kind=8) :: fv,fd,fp,E_int
   integer :: i,imax
   
   fchk=1d0+ftol
   itr=1
   
   call the_path(nim,shape_spin,path,pathlen)
   
   ci = 1
   fxyz1 = 0d0
   fxyz2 = 0d0
   do i_nim=1,nim

      u(i_nim) = 0d0
      do i_m=1,shape_spin(5)
         do i_z=1,shape_spin(4)
            do i_y=1,shape_spin(3)
               do i_x=1,shape_spin(2)
                  vel(:,i_x,i_y,i_z,i_m,i_nim) = 0d0
                  do i = 1,3
                     coo(i,i_x,i_y,i_z,i_m,i_nim) = path(3+i,i_x,i_y,i_z,i_m,i_nim)
                  end do
                  
                  
                  
          
                  call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,ftmp, &
                    &  path(:,:,:,:,:,i_nim),shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext,my_lattice)
                  
                  call project_force(ftmp,path(4:6,i_x,i_y,i_z,i_m,i_nim),fxyz1(:,i_x,i_y,i_z,i_m,i_nim))

                  call local_energy(E_int,i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,path(:,:,:,:,:,i_nim), &
                       & shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,my_lattice)

                  u(i_nim) = u(i_nim) + E_int
               end do
            end do
         end do
      end do
      if (u(i_nim)>u(ci)) then
         ci = i_nim
      end if
   end do
   
   print *,'ci:',ci
   
   !stop
   
   
   
             

                        
   if (en_zero=='I') then
      u0 = u(1)
   elseif (en_zero=='F') then
      u0 = u(nim)
   else
      u0 = 0d0
   end if
   
!   !fname_work = 'triven_'//trim(adjustl(simid))
!      !open(777,file = fname_work, access = 'sequential', action = 'write',status = 'replace')
      
      call tang(nim,shape_spin,1,coo,tau_i)
      fpp(1) = 0d0
      do i_m=1,shape_spin(5)
         do i_z=1,shape_spin(4)
            do i_y=1,shape_spin(3)
               do i_x=1,shape_spin(2)
                  do i=1,3
                     fpp(1) = fpp(1) + fxyz1(i,i_x,i_y,i_z,i_m,1)*tau_i(i,i_x,i_y,i_z,i_m)
                  end do
               end do
            end do
         end do
      end do
            
      
      

      
      do i_nim=2,nim-1
         call tang_spec(nim,shape_spin,i_nim,coo,u,tau)

         call tang(nim,shape_spin,i_nim,coo,tau_i)
         
         fp = 0d0
         fpp(i_nim) = 0d0
         
         do i_m=1,shape_spin(5)
            do i_z=1,shape_spin(4)
               do i_y=1,shape_spin(3)
                  do i_x=1,shape_spin(2)
                     do i=1,3
                        fpp(i_nim) = fpp(i_nim) + fxyz1(i,i_x,i_y,i_z,i_m,i_nim)*tau_i(i,i_x,i_y,i_z,i_m)
                        fp = fp + fxyz1(i,i_x,i_y,i_z,i_m,i_nim)*tau(i,i_x,i_y,i_z,i_m)
                     end do
                  end do
               end do
            end do
         end do
         
         if (i_nim==ci) then
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
                        do i=1,3
                           !print *, fxyz1(i,i_x,i_y,i_z,i_m,i_nim)
                           fxyz1(i,i_x,i_y,i_z,i_m,i_nim) = fxyz1(i,i_x,i_y,i_z,i_m,i_nim) - 2d0*tau(i,i_x,i_y,i_z,i_m)*fp 
                           !print *, fxyz1(i,i_x,i_y,i_z,i_m,i_nim)
                        end do
                     end do
                  end do
               end do
            end do
         else
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
                        do i=1,3
                           !print *, fxyz1(i,i_x,i_y,i_z,i_m,i_nim)
                           fxyz1(i,i_x,i_y,i_z,i_m,i_nim) = fxyz1(i,i_x,i_y,i_z,i_m,i_nim) - tau(i,i_x,i_y,i_z,i_m)*fp + kappa*tau(i,i_x,i_y,i_z,i_m)*(pathlen(i_nim+1)+pathlen(i_nim-1)-2d0*pathlen(i_nim))
                        end do
                     end do
                  end do
               end do
            end do
         end if
         

      end do
      
      call tang(nim,shape_spin,nim,coo,tau_i)
      fpp(nim) = 0d0
      do i_m=1,shape_spin(5)
         do i_z=1,shape_spin(4)
            do i_y=1,shape_spin(3)
               do i_x=1,shape_spin(2)
                  do i=1,3
                     fpp(nim) = fpp(nim) + fxyz1(i,i_x,i_y,i_z,i_m,nim)*tau_i(i,i_x,i_y,i_z,i_m)
                  end do
               end do
            end do
         end do
      end do
      

      
      fchk=0d0
      imax = 1
      do i_nim=1,nim
         do i_m=1,shape_spin(5)
            do i_z=1,shape_spin(4)
               do i_y=1,shape_spin(3)
                  do i_x=1,shape_spin(2)
                     do i=1,3
                        if (dabs(fxyz1(i,i_x,i_y,i_z,i_m,i_nim))>fchk) then
                           fchk = dabs(fxyz1(i,i_x,i_y,i_z,i_m,i_nim))
                           imax = i_nim
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end do
                  
      
      itr=1
      
      open(99,file = 'force_mep.txt',access = 'sequential',action='write',status='replace')
      close(99)
      
      open(99,file = 'force_mep.txt', access = 'sequential', action = 'write',status = 'old',position = 'append')
      write(99,'(i12,a,es16.8E3,a,i3,a,i3)',advance = 'no') itr,'   ',fchk,'   ',imax,'   ',ci
      
      close(99)
      do i=1,nim
               ene(i) = (u(i)-u0)
               dene(i) = -fpp(i)
               rx(i) = pathlen(i)
      end do
      call write_en(nim,rx,ene,dene,rx(nim),'en_path.in',do_norm_rx)
      
      call write_path(nim,shape_spin,path)
      
      
!======================MAIN LOOP============================
      do while ((fchk>ftol).and.(itr<=itrmax))
         ci = 1
         do i_nim=2,nim-1
            u(i_nim) = 0d0
            !print *,'i_nim:',i_nim
            !print *,' '
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
            
                        
                        ax(:,i_x,i_y,i_z,i_m,i_nim) = calc_axis(coo(:,i_x,i_y,i_z,i_m,i_nim),fxyz1(:,i_x,i_y,i_z,i_m,i_nim))
                        
                        coo(:,i_x,i_y,i_z,i_m,i_nim) = path(4:6,i_x,i_y,i_z,i_m,i_nim)+vel(:,i_x,i_y,i_z,i_m,i_nim)*dt+0.5d0*fxyz1(:,i_x,i_y,i_z,i_m,i_nim)/mass*dt*dt

                        call normalize_vec(3,coo(:,i_x,i_y,i_z,i_m,i_nim))
                        ang(i_x,i_y,i_z,i_m,i_nim) = calc_ang(coo(:,i_x,i_y,i_z,i_m,i_nim),path(4:6,i_x,i_y,i_z,i_m,i_nim))
                        path(4:6,i_x,i_y,i_z,i_m,i_nim) = coo(:,i_x,i_y,i_z,i_m,i_nim)
                        !print *,'ang:',ang(i_x,i_y,i_z,i_m)
                     end do
                  end do
               end do
            end do
            
            
         end do
            
         
         do i_nim=2,nim-1
            u(i_nim) = 0d0
            
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
                        
                        call calculate_Beff(i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,ftmp, &
                              &  path(:,:,:,:,:,i_nim),shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext,my_lattice)

                        call project_force(ftmp,path(4:6,i_x,i_y,i_z,i_m,i_nim),fxyz2(:,i_x,i_y,i_z,i_m,i_nim))

                        call local_energy(E_int,i_DM,i_four,i_biq,i_dip,EA,i_x,i_y,i_z,i_m,path(:,:,:,:,:,i_nim), &
                             &shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,h_ext,my_lattice)

                        u(i_nim) = u(i_nim) + E_int

                        !print *,'ene:',u(i_nim)
                     end do
                  end do
               end do
            end do
            
            if (u(i_nim)>u(ci)) then
               ci = i_nim
            end if
            
         end do
         
         
         
         call the_path(nim,shape_spin,path,pathlen)
         
         !print *,'pathlen:',pathlen
         
         do i_nim=2,nim-1
            
            
            
            
            
            call tang_spec(nim,shape_spin,i_nim,coo,u,tau)
            
            call tang(nim,shape_spin,i_nim,coo,tau_i)
         
            fp = 0d0
            fpp(i_nim) = 0d0
         
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
                        do i=1,3
                           fpp(i_nim) = fpp(i_nim) + fxyz2(i,i_x,i_y,i_z,i_m,i_nim)*tau_i(i,i_x,i_y,i_z,i_m)
                           fp = fp + fxyz2(i,i_x,i_y,i_z,i_m,i_nim)*tau(i,i_x,i_y,i_z,i_m)
                        end do
                     end do
                  end do
               end do
            end do
         
            if (i_nim == ci) then
               do i_m=1,shape_spin(5)
                  do i_z=1,shape_spin(4)
                     do i_y=1,shape_spin(3)
                        do i_x=1,shape_spin(2)
                           do i=1,3
                              fxyz2(i,i_x,i_y,i_z,i_m,i_nim) = fxyz2(i,i_x,i_y,i_z,i_m,i_nim) - 2d0*tau(i,i_x,i_y,i_z,i_m)*fp
                           end do
                        end do
                     end do
                  end do
               end do
            else
               do i_m=1,shape_spin(5)
                  do i_z=1,shape_spin(4)
                     do i_y=1,shape_spin(3)
                        do i_x=1,shape_spin(2)
                           do i=1,3
                              fxyz2(i,i_x,i_y,i_z,i_m,i_nim) = fxyz2(i,i_x,i_y,i_z,i_m,i_nim) - tau(i,i_x,i_y,i_z,i_m)*fp + kappa*tau(i,i_x,i_y,i_z,i_m)*(pathlen(i_nim+1)+pathlen(i_nim-1)-2d0*pathlen(i_nim))
                           end do
                        end do
                     end do
                  end do
               end do
            end if
         end do
         !stop
         fv = 0d0
         fd = 0d0
         do i_nim=2,nim-1
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
                        call rotate_vec(vel(:,i_x,i_y,i_z,i_m,i_nim),ax(:,i_x,i_y,i_z,i_m,i_nim),ang(i_x,i_y,i_z,i_m,i_nim),veltmp)
                        call rotate_vec(fxyz1(:,i_x,i_y,i_z,i_m,i_nim),ax(:,i_x,i_y,i_z,i_m,i_nim),ang(i_x,i_y,i_z,i_m,i_nim),ftmp)
                        !vel(:,i_x,i_y,i_z,i_m,i_nim) = vel(:,i_x,i_y,i_z,i_m,i_nim) + 0.5d0*(fxyz1(:,i_x,i_y,i_z,i_m,i_nim)+fxyz2(:,i_x,i_y,i_z,i_m,i_nim))/mass*dt
                        vel(:,i_x,i_y,i_z,i_m,i_nim) = veltmp + 0.5d0*(ftmp+fxyz2(:,i_x,i_y,i_z,i_m,i_nim))/mass*dt
                        do i=1,3
                           fv = fv + vel(i,i_x,i_y,i_z,i_m,i_nim)*fxyz2(i,i_x,i_y,i_z,i_m,i_nim)
                           fd = fd + fxyz2(i,i_x,i_y,i_z,i_m,i_nim)*fxyz2(i,i_x,i_y,i_z,i_m,i_nim)
                        end do
                     end do
                  end do
               end do
            end do
         end do
         
         
         fv = 0d0
         fd = 0d0
         do i_nim=1,nim
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
                        
                        do i=1,3
                           fv = fv + vel(i,i_x,i_y,i_z,i_m,i_nim)*fxyz2(i,i_x,i_y,i_z,i_m,i_nim)
                           fd = fd + fxyz2(i,i_x,i_y,i_z,i_m,i_nim)*fxyz2(i,i_x,i_y,i_z,i_m,i_nim)
                        end do
                     end do
                  end do
               end do
            end do
         end do
         
         if (fv<0d0) then
            do i_nim=1,nim
               do i_m=1,shape_spin(5)
                  do i_z=1,shape_spin(4)
                     do i_y=1,shape_spin(3)
                        do i_x=1,shape_spin(2)
                           do i=1,3
                              vel(i,i_x,i_y,i_z,i_m,i_nim) = 0d0
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         else
            do i_nim=1,nim
               do i_m=1,shape_spin(5)
                  do i_z=1,shape_spin(4)
                     do i_y=1,shape_spin(3)
                        do i_x=1,shape_spin(2)
                           do i=1,3
                              vel(i,i_x,i_y,i_z,i_m,i_nim) = fxyz2(i,i_x,i_y,i_z,i_m,i_nim)*fv/fd
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end if
         
         
         
         
         imax=1
         fchk=0d0

         do i_nim=2,nim-1
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
                        do i=1,3
                           fxyz1(i,i_x,i_y,i_z,i_m,i_nim) = fxyz2(i,i_x,i_y,i_z,i_m,i_nim)
                           if (dabs(fxyz1(i,i_x,i_y,i_z,i_m,i_nim))>fchk) then
                              fchk = dabs(fxyz1(i,i_x,i_y,i_z,i_m,i_nim))
                              imax = i_nim
                           end if
                        end do
                     end do
                  end do
               end do
            end do
         end do
         
         
         
         
         
         call tang(nim,shape_spin,1,coo,tau_i)
         fpp(1) = 0d0
         do i_m=1,shape_spin(5)
            do i_z=1,shape_spin(4)
               do i_y=1,shape_spin(3)
                  do i_x=1,shape_spin(2)
                     do i=1,3
                        fpp(1) = fpp(1) + fxyz1(i,i_x,i_y,i_z,i_m,1)*tau_i(i,i_x,i_y,i_z,i_m)
                     end do
                  end do
               end do
            end do
         end do
         
         
         call tang(nim,shape_spin,nim,coo,tau_i)
         fpp(nim) = 0d0
         do i_m=1,shape_spin(5)
            do i_z=1,shape_spin(4)
               do i_y=1,shape_spin(3)
                  do i_x=1,shape_spin(2)
                     do i=1,3
                        fpp(nim) = fpp(nim) + fxyz1(i,i_x,i_y,i_z,i_m,nim)*tau_i(i,i_x,i_y,i_z,i_m)
                     end do
                  end do
               end do
            end do
         end do
         
         
         
         itr=itr+1

         if (mod(itr,every).eq.0) then
            call prn_gneb_progress(itr, itrmax, fchk, imax,'Y',ci)
      
            open(99,file = 'force_mep.txt', access = 'sequential', action = 'write',status = 'old',position = 'append')
            write(99,'(i12,a,es16.8E3,a,i3,a,i3)',advance = 'no') itr,'   ',fchk,'   ',imax,'   ',ci
            close(99)
            
            do i=1,nim
               ene(i) = (u(i)-u0)
               dene(i) = -fpp(i)
               rx(i) = pathlen(i)
            end do
            call write_en(nim,rx,ene,dene,rx(nim),'en_path.out',do_norm_rx)
            call write_path(nim,shape_spin,path)
            
         
         end if
      end do
   
      if (itr>itrmax) then
         write(*,*) 'WARNING: exceeded maximum iterations in CI-GNEB'
      end if
   
      
      
      
      do i=1,nim
         ene(i) = (u(i)-u0)
         dene(i) = -fpp(i)
         rx(i) = pathlen(i)
      end do
   

end subroutine find_path_ci


subroutine prn_gneb_progress(itr, itrmax, fchk, imax,do_ci,ci)
implicit none
   integer, intent(in) :: itr, itrmax
   real(kind=8), intent(in) :: fchk
   integer, intent(in) :: imax,ci
   character(len=1) :: do_ci
   character(35) :: num,num1,num2,num3
      
   write(num,'(I8)') idnint(real(itr,kind=8)/real(itrmax,kind=8)*100d0)
   write(num1,'(es16.8E3)') fchk
   write(num2,'(i8)') imax
   write(num3,'(i8)') ci
      
   if (do_ci.eq.'Y') then
      write (*,'(2x,8a)') 'MP  ',trim(adjustl(num)),'% of itrmax.   fchk: ',trim(adjustl(num1)),'   imax: ',trim(adjustl(num2)),'   ci: ',trim(adjustl(num3))
   else
      write (*,'(2x,6a)') 'MP  ',trim(adjustl(num)),'% of itrmax.   fchk: ',trim(adjustl(num1)),'   imax: ',trim(adjustl(num2))
   end if
  
  
end subroutine prn_gneb_progress


!> Estimate tangent to the path according to the special definition: doi:10.1016/j.cpc.2015.07.001
   subroutine tang_spec(nim,shape_spin,im,coo,u,tau)
   implicit none
   integer, intent(in) :: im,nim,shape_spin(5)
   real(kind=8), intent(in) :: coo(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim), u(nim)
   real(kind=8), intent(out) :: tau(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
   real(kind=8) :: u1, u2, u3,taup(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),taum(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),dumin,dumax,tmp,tau_tmp(3)
   integer :: i_m,i_x,i_y,i_z,i
   u1=u(im-1)
   u2=u(im)
   u3=u(im+1)
   !print *,'ene:',u1,u2,u3
   if (u3>u2.and.u2>u1) then
      !print *,'I am here!!!'
      do i_m=1,shape_spin(5)
         do i_z=1,shape_spin(4)
            do i_y=1,shape_spin(3)
               do i_x=1,shape_spin(2)
                  do i=1,3
                     tau(i,i_x,i_y,i_z,i_m) = coo(i,i_x,i_y,i_z,i_m,im+1)-coo(i,i_x,i_y,i_z,i_m,im)
                  end do
               end do
            end do
         end do
      end do
   elseif (u1>u2.and.u2>u3) then
      !print *,'I am here'
      do i_m=1,shape_spin(5)
         do i_z=1,shape_spin(4)
            do i_y=1,shape_spin(3)
               do i_x=1,shape_spin(2)
                  do i=1,3
                     tau(i,i_x,i_y,i_z,i_m) = coo(i,i_x,i_y,i_z,i_m,im)-coo(i,i_x,i_y,i_z,i_m,im-1)
                  end do
               end do
            end do
         end do
      end do
   else 
      !print *,'I am here!'
      do i_m=1,shape_spin(5)
         do i_z=1,shape_spin(4)
            do i_y=1,shape_spin(3)
               do i_x=1,shape_spin(2)
                  do i=1,3
                     taup(i,i_x,i_y,i_z,i_m) = coo(i,i_x,i_y,i_z,i_m,im+1)-coo(i,i_x,i_y,i_z,i_m,im)
                     taum(i,i_x,i_y,i_z,i_m) = coo(i,i_x,i_y,i_z,i_m,im)-coo(i,i_x,i_y,i_z,i_m,im-1)
                  end do
               end do
            end do
         end do
      end do
      dumax=dabs(u3-u2)
      dumin=dabs(u1-u2)
      if (dumax<dumin) then
         tmp=dumax
         dumax=dumin
         dumin=tmp
      end if
      if (u3>u1) then
        
         tau=dumax*taup+dumin*taum
      else
        !print *,'imhere', 'u1 = ', u1, 'u2 = ', u2,'u3 = ', u3
         tau=dumin*taup+dumax*taum
      end if
   end if

   tmp = 0d0
      do i_m=1,shape_spin(5)
         do i_z=1,shape_spin(4)
            do i_y=1,shape_spin(3)
               do i_x=1,shape_spin(2)
                  tau_tmp(:) = tau(:,i_x,i_y,i_z,i_m)
                  !print *,'tau_tmp:',tau_tmp
                  !print *,'norm_vec:',norm_vec(3,tau(:,i_x,i_y,i_z,i_m))
                  call project_force(tau_tmp,coo(:,i_x,i_y,i_z,i_m,im),tau(:,i_x,i_y,i_z,i_m))
                  !print *,'norm_vec:',norm_vec(3,coo(:,i_x,i_y,i_z,i_m,im))
                  !print *,' '
                  do i=1,3
                     tmp = tmp+tau(i,i_x,i_y,i_z,i_m)*tau(i,i_x,i_y,i_z,i_m)
                  end do
               end do
            end do
         end do
      end do
      
      tmp = dsqrt(tmp)
      !print *,'tmp:',tmp
      tau = tau/tmp
   
end subroutine tang_spec
   
   subroutine tang(nim,shape_spin,im,coo,tau)
   implicit none
   integer, intent(in) :: im,nim,shape_spin(5)
   real(kind=8), intent(in) :: coo(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   real(kind=8), intent(out) :: tau(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
   real(kind=8) :: tmp,tau_tmp(3)
   integer :: i_m,i_x,i_y,i_z,i
   
   if (im==1) then
      tau = coo(:,:,:,:,:,im+1)-coo(:,:,:,:,:,im)
   elseif (im==nim) then
      tau = coo(:,:,:,:,:,im)-coo(:,:,:,:,:,im-1)
   else
      tau = coo(:,:,:,:,:,im+1)-coo(:,:,:,:,:,im-1)
   end if

   
   tmp = 0d0
      do i_m=1,shape_spin(5)
         do i_z=1,shape_spin(4)
            do i_y=1,shape_spin(3)
               do i_x=1,shape_spin(2)
                  tau_tmp(:) = tau(:,i_x,i_y,i_z,i_m)
                  call project_force(tau_tmp,coo(:,i_x,i_y,i_z,i_m,im),tau(:,i_x,i_y,i_z,i_m))
                  do i=1,3
                     tmp = tmp+tau(i,i_x,i_y,i_z,i_m)*tau(i,i_x,i_y,i_z,i_m)
                  end do
               end do
            end do
         end do
      end do
      
      tmp = dsqrt(tmp)
      tau = tau/tmp
   
end subroutine tang



!> Convert effective fields to the format used in VPO
subroutine project_force(beff,coo,fxyz)
      implicit none 
      real(kind=8), intent(in) :: beff(3),coo(3)
      real(kind=8), intent(out) :: fxyz(3)
      integer :: i
      real(kind=8) :: tmp
      
      tmp = 0d0
      do i=1,3

         tmp = tmp+beff(i)*coo(i)
      end do
      
      do i=1,3
         fxyz(i) = beff(i) - tmp*coo(i)
      end do
   
end subroutine project_force

!> Calculate poly-geodesic length of the path
subroutine the_path(nim,shape_spin,path,pathlen)
   implicit none
   integer, intent(in) :: nim,shape_spin(5)
   real(kind=8), intent(in) :: path(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   real(kind=8), intent(out) :: pathlen(nim)
   real(kind=8) :: tmp,l1
   integer :: i_m,i_x,i_y,i_z,i_nim

      pathlen(1) = 0d0
      
      do i_nim=2,nim
         tmp = 0d0
         do i_m=1,shape_spin(5)
            do i_z=1,shape_spin(4)
               do i_y=1,shape_spin(3)
                  do i_x=1,shape_spin(2)
                     l1 = calc_ang(path(4:6,i_x,i_y,i_z,i_m,i_nim-1),path(4:6,i_x,i_y,i_z,i_m,i_nim))
                     
                     tmp = tmp+l1*l1
                  end do
               end do
            end do
         end do
         pathlen(i_nim) = pathlen(i_nim-1) + dsqrt(tmp)
      end do
      

end subroutine the_path

subroutine set_gneb_defaults()
   
implicit none
   momfile_i     = 'momfile_i'
   momfile_f     = 'momfile_f'
   restartfile_if = 'nofile'
   restartfile_path = 'nofile'
   amp_rnd = 0.0
   amp_rnd_path = 0.0

   minftol = 0.000000001d0
   mintraj_step = 100
   vpodt = 0.01d0
   vpomass = 1d0
   minitrmax = 10000000
    
      !Parameters for GNEB calculations
   initpath = 1
   spring = 0.5d0
   mepftol = 0.001d0
   mepftol_ci = 0.00001d0
   mepitrmax = 10000000
   meptraj_step = 100
   do_gneb = 'Y'
   do_gneb_ci = 'N'
   do_norm_rx = 'N'
   en_zero = 'N'
   nim = 10
    
    !Parameters for energy interpolation along the MEP
   sample_num = 500
end subroutine set_gneb_defaults


subroutine read_gneb_parameters()

    implicit none
    
    integer, parameter :: io=357
    character(len=50) :: keyword,cache
    integer :: key_len,i_err,i_errb
    logical :: comment,exi

!    character(len=1) :: OPT_flag_str, ip_adapt_flag_str, OPT_printcores_flag_str
    
    inquire (file='GNEB.in',exist=exi)
    if (.not. exi) then
      write(6,*) 'no input file for the GNEB method!'
      STOP
    endif
    
    
    
    
    open (io,file='GNEB.in',form='formatted',status='old',action='read')
    
    do
10     continue
       ! Read file character for character until first whitespace
       keyword=""
       call read_keyword(io,keyword,key_len,i_errb)

       ! check for comment markers (currently % and #)
       comment=(scan(trim(keyword),'%')==1).or.(scan(trim(keyword),'#')==1).or.&
            (scan(trim(keyword),'*')==1).or.(scan(trim(keyword),'=')==1.or.&
            (scan(trim(keyword),'!')==1))

       if (comment) then
          read(io,*)
       else
       ! Parse keyword
       keyword=trim(keyword)
       select case(keyword)


          
       case('momfile_i')
          read(io,'(a)',iostat=i_err) cache
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          momfile_i=trim(adjustl(cache))
          
       case('momfile_f')
          read(io,'(a)',iostat=i_err) cache
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          momfile_f=trim(adjustl(cache))
          
       case('restartfile_if')
          read(io,'(a)',iostat=i_err) cache
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          restartfile_if=trim(adjustl(cache))
          
       case('restartfile_path')
          read(io,'(a)',iostat=i_err) cache
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          restartfile_path=trim(adjustl(cache))
       
       case('spring')
          read(io,*,iostat=i_err) spring
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          
       case('initpath')
          read(io,*,iostat=i_err) initpath
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
       
       case('amp_rnd')
          read(io,*,iostat=i_err) amp_rnd
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          
       case('amp_rnd_path')
          read(io,*,iostat=i_err) amp_rnd_path
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
           
       case('min_itrmax')
          read(io,*,iostat=i_err) minitrmax
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          
       case('mintraj_step')
          read(io,*,iostat=i_err) mintraj_step
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          
       case('min_ftol')
          read(io,*,iostat=i_err) minftol
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          
       case('mep_itrmax')
          read(io,*,iostat=i_err) mepitrmax
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          
       case('meptraj_step')
          read(io,*,iostat=i_err) meptraj_step
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          
       case('mep_ftol')
          read(io,*,iostat=i_err) mepftol
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          
       case('mep_ftol_ci')
          read(io,*,iostat=i_err) mepftol_ci
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          
       case('do_gneb')
          read(io,'(a)',iostat=i_err) do_gneb
          if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword), ' data',i_err
          
       case('do_gneb_ci')
          read(io,'(a)',iostat=i_err) do_gneb_ci
          if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword), ' data',i_err
          
       case('do_norm_rx')
          read(io,'(a)',iostat=i_err) do_norm_rx
          if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword), ' data',i_err
          
       case('en_zero')
          read(io,'(a)',iostat=i_err) en_zero
          if(i_err/=0) write(*,*) 'ERROR: Reading ', trim(keyword), ' data',i_err
          
       case('vpo_dt')
          read(io,*,iostat=i_err) vpodt
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          
       case('vpo_mass')
          read(io,*,iostat=i_err) vpomass
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err
          
       case('sample_num')
          read(io,*,iostat=i_err) sample_num
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err

       case('nim')
          read(io,*,iostat=i_err) nim
          if(i_err/=0) write(*,*) 'ERROR: Reading ',trim(keyword),' data',i_err


      case default
          if(len(trim(keyword))>0) then
             print *,"Keyword '",trim(keyword),"' is not recognized"
             !call ErrorHandling_ERROR("Keyword '"//trim(keyword)//"' is not recognized")
             read(io,*)
          end if

      end select
      end if

      ! End of file
      if (i_errb==20) goto 20
      ! End of row
      if (i_errb==10) goto 10
    end do

20  continue
      close(io)
    return
  end subroutine read_gneb_parameters

subroutine read_keyword(io,keyword,keyword_len,i_err)
    ! 
implicit none
    !
    integer, intent(in) :: io
    character(len=*), intent(out) :: keyword  
    integer, intent(out) :: keyword_len   
    integer, intent(out) :: i_err  
    !
    logical :: rd_done,rd_start
    !
    rd_done=.false.
    rd_start=.false.
    keyword_len=0
    do while(.not.rd_done.and.keyword_len<len(keyword))
       keyword_len=keyword_len+1
       read(io,'(a1)',advance='no',end=20,eor=10) keyword(keyword_len:keyword_len)
       rd_start=rd_start.or.(keyword(keyword_len:keyword_len)/=" ")
       rd_done=rd_start.and.(keyword(keyword_len:keyword_len)==" ")
    end do
    ! happy ending
    i_err=0
    keyword=adjustl(keyword(1:keyword_len)//'')
    return
    ! final word
10  continue
    i_err=10
    keyword=adjustl(keyword(1:keyword_len)//'')
    !        print *,'<<<<',trim(keyword),'>>>>>'
    return
    ! end of file
20  continue
    i_err=20
    keyword=adjustl(keyword(1:keyword_len)//'')
    return
    !
  end subroutine read_keyword 
  
subroutine geodesic_path_one(state,nim,ni,nf,ax,path)
use mtprng
implicit none
      integer, intent(in) :: nim
      real(kind=8), intent(in) :: ni(3),nf(3)
      real(kind=8), intent(inout) :: ax(3)
      real(kind=8), intent(inout) :: path(3,nim)
      type(mtprng_state), intent(inout) :: state
      real(kind=8), parameter :: pi = 3.14159265358979323d0
      real(kind=8) :: dtheta,theta,angle,tmp,vec(3),eps=epsilon(angle),pr
      integer :: i,j
      
      angle = calc_ang(ni,nf)
      if (angle<eps) then
         do i = 1,3
            vec(i) = (nf(i) - ni(i))/(nim-1)
         end do
         path(:,1) = ni(:)
         path(:,nim) = nf(:)
         do i=2,nim-1
            do j=1,3
               path(j,i) = ni(j) + real((i-1),kind=8)*vec(j)
            end do
            call normalize_vec(3,path(:,i))
         end do
               
      elseif (dabs(angle-pi)<eps) then
      !   write(*,*) 'i am here!'
         pr = 0d0
         do j=1,3
            pr = pr + ax(j)*ni(j)
         end do
         ax(:) = ax(:) - pr*ni(:)
         tmp = norm_vec(3,ax)
         do while (tmp<eps)
            pr = 0d0
            do j=1,3
               ax(j) = dsign(1d0,2d0*mtprng_rand_real1(state)-1d0)*(mtprng_rand_real1(state)+1d0)
               pr = pr + ax(j)*ni(j)
            end do
            ax(:) = ax(:) - pr*ni(:)
            tmp = norm_vec(3,ax)
         end do
         call normalize_vec(3,ax)
         dtheta = pi/(nim-1)
         path(:,1) = ni(:)
         path(:,nim) = nf(:)
      
         do i=2,nim-1
            theta = (i-1)*dtheta
            path(1,i) = ni(1)*dcos(theta) + dsin(theta)*(ax(2)*ni(3)-ax(3)*ni(2))
            path(2,i) = ni(2)*dcos(theta) - dsin(theta)*(ax(1)*ni(3)-ax(3)*ni(1))
            path(3,i) = ni(3)*dcos(theta) + dsin(theta)*(ax(1)*ni(2)-ax(2)*ni(1))
            call normalize_vec(3,path(:,i))
         end do
   
      else
         ax = calc_axis(ni,nf)
      

         dtheta = angle/(nim-1)

         path(:,1) = ni(:)
         path(:,nim) = nf(:)

   

         do i=2,nim-1
            theta = (i-1)*dtheta
            path(1,i) = ni(1)*dcos(theta) + dsin(theta)*(ax(2)*ni(3)-ax(3)*ni(2))
            path(2,i) = ni(2)*dcos(theta) - dsin(theta)*(ax(1)*ni(3)-ax(3)*ni(1))
            path(3,i) = ni(3)*dcos(theta) + dsin(theta)*(ax(1)*ni(2)-ax(2)*ni(1))
            call normalize_vec(3,path(:,i))
         end do   
   
      end if
   
end subroutine geodesic_path_one


   subroutine geodesic_path(state,nim,amp_rnd,spini,spinf,shape_spin,path)
   use mtprng
   implicit none
   type(mtprng_state), intent(inout) :: state
   integer, intent(in) :: nim,shape_spin(5)
   real(kind=8), intent(in) :: spini(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),spinf(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),amp_rnd
   real(kind=8), intent(inout) :: path(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   integer :: i,i_x,i_y,i_z,i_m,i_nim
   real(kind=8) :: u(3),path_one(3,nim),ax(3)
      
      ax = 0d0
      
      do i_m=1,shape_spin(5)
       do i_z=1,shape_spin(4)
        do i_y=1,shape_spin(3)
         do i_x=1,shape_spin(2)
            do i_nim = 1,nim
               do i=1,shape_spin(1)
                  path(i,i_x,i_y,i_z,i_m,i_nim) = spini(i,i_x,i_y,i_z,i_m)
               end do
            end do
            call geodesic_path_one(state,nim,spini(4:6,i_x,i_y,i_z,i_m),spinf(4:6,i_x,i_y,i_z,i_m),ax,path_one)
            do i_nim=1,nim
               path(4:6,i_x,i_y,i_z,i_m,i_nim) = path_one(1:3,i_nim)
            end do
            do i_nim = 2,nim-1
               u(1)=2d0*(mtprng_rand_real1(state)-1d0)
               u(2)=2d0*(mtprng_rand_real1(state)-1d0)
               u(3)=2d0*(mtprng_rand_real1(state)-1d0)
               path(4:6,i_x,i_y,i_z,i_m,i_nim) = path(4:6,i_x,i_y,i_z,i_m,i_nim)+amp_rnd*u
               call normalize_vec(3,path(4:6,i_x,i_y,i_z,i_m,i_nim))
            end do
            
         end do
        end do
       end do
      end do
      
      
   end subroutine geodesic_path

   !> Calculate angle between two 3-vectors
   real(kind=8) function calc_ang(n1,n2)
   implicit none
   real(kind=8), intent(in) :: n1(3), n2(3) !n1 and n2 have to be normalized
   !real(kind=8) :: calc_ang
   real(kind=8) :: n(3),prod,tmp
   
      n(1) = n1(2)*n2(3)-n1(3)*n2(2)
      n(2) =-n1(1)*n2(3)+n1(3)*n2(1)
      n(3) = n1(1)*n2(2)-n1(2)*n2(1)
      
      tmp = dsqrt(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))
      prod = n1(1)*n2(1)+n1(2)*n2(2)+n1(3)*n2(3)
      
      calc_ang = datan2(tmp,prod)
      return
   end function calc_ang
   
      !> Find norm of a vector
   real(kind=8) function norm_vec(N,vec)
   implicit none
   integer, intent(in) :: N
   real(kind=8), intent(in) :: vec(N)
   real(kind=8) :: tmp
   integer :: i
      tmp = 0d0
      do i=1,N
         !print *,vec(i)
         tmp = tmp + vec(i)*vec(i)
      end do
      !print *,' '
      norm_vec = dsqrt(tmp)
   
   end function norm_vec
   
      !> Normalize vector
   subroutine normalize_vec(N,vec)
   implicit none
   integer, intent(in) :: N
   real(kind=8), intent(inout) :: vec(N)
   real(kind=8) :: tmp
   integer :: i
   
      tmp = norm_vec(N,vec)
   
      do i=1,N
         vec(i) = vec(i)/tmp
      end do
      
   end subroutine normalize_vec
   
!> Rotate 3-vector v_in by angle ang around axis ax
subroutine rotate_vec(v_in,ax,ang,v_out)
   implicit none
   real(kind=8), intent(in) :: v_in(3),ax(3),ang
   real(kind=8), intent(out) :: v_out(3)
   integer :: i
   real(kind=8) :: tmp,sinang,cosang
   
   tmp = 0d0
   do i = 1,3
      tmp = tmp + ax(i)*v_in(i)
   end do
   
   sinang = dsin(ang)
   cosang = dcos(ang)
   
   v_out(1) = v_in(1)*cosang + sinang*(ax(2)*v_in(3)-ax(3)*v_in(2))+(1d0-cosang)*tmp*ax(1)
   v_out(2) = v_in(2)*cosang - sinang*(ax(1)*v_in(3)-ax(3)*v_in(1))+(1d0-cosang)*tmp*ax(2)
   v_out(3) = v_in(3)*cosang + sinang*(ax(1)*v_in(2)-ax(2)*v_in(1))+(1d0-cosang)*tmp*ax(3)
end subroutine rotate_vec
   
!> Calculate axes of rotation based on magnetic configuration and effective fields
   subroutine calc_axis_all(shape_spin,m_in,f_in,ax)
   implicit none
   integer, intent(in) :: shape_spin(5)
   real(kind=8), intent(in) :: m_in(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),f_in(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
   real(kind=8), intent(out) :: ax(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
   real(kind=8) :: x(3),y(3),axi(3)
   integer :: i_x,i_y,i_z,i_m
   
      do i_m=1,shape_spin(5)
       do i_z=1,shape_spin(4)
        do i_y=1,shape_spin(3)
         do i_x=1,shape_spin(2)
            x = m_in(:,i_x,i_y,i_z,i_m)
      
            y = f_in(:,i_x,i_y,i_z,i_m)
            
      
            axi = calc_axis(x,y)
      
            ax(:,i_x,i_y,i_z,i_m) = axi

         end do
        end do
       end do
      end do
   
   end subroutine calc_axis_all

      !> Calculate axis of rotation based on two 3-vectors
   function calc_axis(m_in,f_in)
   implicit none
   real(kind=8), parameter :: fact=10000000d0
   real(kind=8), intent(in) :: m_in(3),f_in(3)
   real(kind=8), dimension(3) :: calc_axis
   real(kind=8) :: x(3),tmp,y(3),a,b,c,eps=epsilon(a)
   
     
      
      tmp = norm_vec(3,f_in)
      if (tmp<eps) then
         a = f_in(1)*fact
         b = f_in(2)*fact
         c = f_in(3)*fact
      else
         a = f_in(1)
         b = f_in(2)
         c = f_in(3)
      end if
    
      
      y(1) = dsign(a,f_in(1))
      y(2) = dsign(b,f_in(2))
      y(3) = dsign(c,f_in(3))
      
      
      
      x(1) = m_in(2)*y(3)-m_in(3)*y(2)
      x(2) =-m_in(1)*y(3)+m_in(3)*y(1)
      x(3) = m_in(1)*y(2)-m_in(2)*y(1)

      tmp = dsqrt(x(1)*x(1)+x(2)*x(2)+x(3)*x(3))

      calc_axis(1) = x(1)/tmp
      calc_axis(2) = x(2)/tmp
      calc_axis(3) = x(3)/tmp

   end function calc_axis
   
subroutine find_SP(nim,rx,c,imax,drx0)
      implicit none
  
      integer, intent(in) :: nim                   !< Number of images
      real(kind=8), intent(in) :: rx(nim)         !< Reaction coordinate
      real(kind=8), intent(in) :: c(4,nim)        !< Coefficients of the piecewise Hermite polynomials
      integer, intent(out) :: imax
      real(kind=8),intent(out) :: drx0
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
   
   
subroutine find_SP_conf(state,shape_spin,spini,spinf,l1,l2,lsp,spinsp)
use mtprng
implicit none
   type(mtprng_state), intent(inout) :: state
   integer, intent(in) :: shape_spin(5)
   real(kind=8), intent(in) :: spini(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
   real(kind=8), intent(in) :: spinf(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
   real(kind=8), intent(in) :: l1,l2,lsp
   real(kind=8), intent(inout) :: spinsp(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
   integer :: i,i_x,i_y,i_z,i_m
   real(kind=8) :: ax(3)
      
      ax = 0d0
      
      do i_m=1,shape_spin(5)
       do i_z=1,shape_spin(4)
        do i_y=1,shape_spin(3)
         do i_x=1,shape_spin(2)
            do i=1,3
               spinsp(i,i_x,i_y,i_z,i_m) = spini(i,i_x,i_y,i_z,i_m)
            end do
            call find_SP_conf_one(state,spini(4:6,i_x,i_y,i_z,i_m),spinf(4:6,i_x,i_y,i_z,i_m),l1,l2,lsp,ax,spinsp(4:6,i_x,i_y,i_z,i_m))
         end do
        end do
       end do
      end do
end subroutine find_SP_conf


subroutine find_SP_conf_one(state,ni,nf,l1,l2,lsp,ax,nsp)
use mtprng
      implicit none
      type(mtprng_state), intent(inout) :: state
      real(kind=8), intent(in) :: ni(3),nf(3),l1,l2,lsp
      real(kind=8), intent(inout) :: ax(3)
      real(kind=8), intent(out) :: nsp(3)
      real(kind=8), parameter :: pi = 3.14159265358979323d0
      real(kind=8) :: dl,theta,angle,tmp,vec(3),eps=epsilon(angle),pr
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
         !write(*,*) 'i am here!'
         
         pr = 0.0d0
         do j=1,3
            pr = pr + ax(j)*ni(j)
         end do
         ax(:) = ax(:) - pr*ni(:)
         tmp = norm_vec(3,ax)
         
         do while (tmp<eps)
            pr = 0d0
            do j=1,3
               ax(j) = dsign(1d0,2d0*mtprng_rand_real1(state)-1d0)*(mtprng_rand_real1(state)+1d0)
               pr = pr + ax(j)*ni(j)
            end do
            ax(:) = ax(:) - pr*ni(:)
            tmp = norm_vec(3,ax)
         end do
         call normalize_vec(3,ax)
         
         theta = pi*dl/(l2-l1)
         
         nsp(1) = ni(1)*dcos(theta) + dsin(theta)*(ax(2)*ni(3)-ax(3)*ni(2))
         nsp(2) = ni(2)*dcos(theta) - dsin(theta)*(ax(1)*ni(3)-ax(3)*ni(1))
         nsp(3) = ni(3)*dcos(theta) + dsin(theta)*(ax(1)*ni(2)-ax(2)*ni(1))
         call normalize_vec(3,nsp)
         
   
      else
         ax = calc_axis(ni,nf)
      

         theta = angle*dl/(l2-l1)

    
         nsp(1) = ni(1)*dcos(theta) + dsin(theta)*(ax(2)*ni(3)-ax(3)*ni(2))
         nsp(2) = ni(2)*dcos(theta) - dsin(theta)*(ax(1)*ni(3)-ax(3)*ni(1))
         nsp(3) = ni(3)*dcos(theta) + dsin(theta)*(ax(1)*ni(2)-ax(2)*ni(1))
         call normalize_vec(3,nsp)
        
      end if
   
end subroutine find_SP_conf_one
   
subroutine read_inifin(state,file_ini,file_fin,amp_rnd,shape_spin,spini,spinf)
   use mtprng
   implicit none
   integer, parameter :: io=358
   type(mtprng_state), intent(inout) :: state
   character(len=*), intent(in) :: file_ini,file_fin
   integer, intent(in) :: shape_spin(5)
   real(kind=8), intent(in):: amp_rnd
   real(kind=8), intent(inout) :: spini(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),spinf(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
   logical :: exi
   integer :: i_x,i_y,i_z,i_m
   real(kind=8) :: u(3)
   
   inquire (file=file_ini,exist=exi)
    if (.not. exi) then
      write(6,*) 'no input file for the initial state!'
      STOP
    endif
    
    inquire (file=file_fin,exist=exi)
    if (.not. exi) then
      write(6,*) 'no input file for the final state!'
      STOP
    endif
    
    open (io,file=file_ini,form='formatted',status='old',action='read')
       do i_x=1,shape_spin(2)
        do i_y=1,shape_spin(3)
         do i_z=1,shape_spin(4)
          do i_m=1,shape_spin(5)
            read(io,*) spini(1:shape_spin(1),i_x,i_y,i_z,i_m)
            u(1)=2d0*(mtprng_rand_real1(state)-1d0)
            u(2)=2d0*(mtprng_rand_real1(state)-1d0)
            u(3)=2d0*(mtprng_rand_real1(state)-1d0)
            spini(4:6,i_x,i_y,i_z,i_m) = spini(4:6,i_x,i_y,i_z,i_m)+amp_rnd*u
            call normalize_vec(3,spini(4:6,i_x,i_y,i_z,i_m))
          end do  
         end do
        end do
      end do
    close(io)
    
    open (io,file=file_fin,form='formatted',status='old',action='read')
       do i_x=1,shape_spin(2)
        do i_y=1,shape_spin(3)
         do i_z=1,shape_spin(4)
          do i_m=1,shape_spin(5)
            read(io,*) spinf(1:shape_spin(1),i_x,i_y,i_z,i_m)
            u(1)=2d0*(mtprng_rand_real1(state)-1d0)
            u(2)=2d0*(mtprng_rand_real1(state)-1d0)
            u(3)=2d0*(mtprng_rand_real1(state)-1d0)
            spinf(4:6,i_x,i_y,i_z,i_m) = spinf(4:6,i_x,i_y,i_z,i_m)+amp_rnd*u
            call normalize_vec(3,spinf(4:6,i_x,i_y,i_z,i_m))
          end do
         end do
        end do
      end do
    close(io)
     
end subroutine read_inifin
   
   
subroutine write_path(nim,shape_spin,path)
   implicit none
   
   integer, intent(in) :: nim,shape_spin(5)
   real(kind=8), intent(in) :: path(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   integer :: i_nim
   
      do i_nim=1,nim
         call WriteSpinAndCorrFile('image-GNEB',i_nim,path(:,:,:,:,:,i_nim),shape_spin)
         call CreateSpinFile('povray-GNEB',i_nim,path(:,:,:,:,:,i_nim),shape_spin)
      end do
end subroutine write_path


subroutine read_path(state,fname_part,nim,amp_rnd,shape_spin,path,exists)
use mtprng
   implicit none
   integer, parameter :: io=358
   type(mtprng_state), intent(inout) :: state
   character(len=*), intent(in) :: fname_part
   integer, intent(in) :: nim,shape_spin(5)
   real(kind=8), intent(in) :: amp_rnd
   real(kind=8), intent(inout) :: path(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim)
   logical, intent(out) :: exists
   integer :: i_x,i_y,i_z,i_m,i_nim
   real(kind=8) :: u(3)
   character(len=8) :: num
   character(len=50) :: fname
   
   do i_nim = 1,nim
      write(num,'(I8)') i_nim
      fname = trim(adjustl(fname_part))//'_'//trim(adjustl(num))//'.dat'
      inquire(file=fname,exist=exists)
      if (exists) then
         open (io,file=trim(adjustl(fname)),form='formatted',status='old',action='read')
         do i_x=1,shape_spin(2)
            do i_y=1,shape_spin(3)
               do i_z=1,shape_spin(4)
                  do i_m=1,shape_spin(5)
                     read(io,*) path(1:shape_spin(1),i_x,i_y,i_z,i_m,i_nim)
                     if ((i_nim.ne.1).and.(i_nim.ne.nim)) then
                        u(1)=2d0*(mtprng_rand_real1(state)-1d0)
                        u(2)=2d0*(mtprng_rand_real1(state)-1d0)
                        u(3)=2d0*(mtprng_rand_real1(state)-1d0)
                        path(4:6,i_x,i_y,i_z,i_m,i_nim) = path(4:6,i_x,i_y,i_z,i_m,i_nim)+amp_rnd*u
                        call normalize_vec(3,path(4:6,i_x,i_y,i_z,i_m,i_nim))
                     end if
                  end do
               end do
            end do
         end do
         close(io)
      else
         write(*,*) 'ERROR: File ',trim(adjustl(fname)), ' does not exist. Path not loaded.'
         path = 0d0
         return
      end if
   end do

   write(*,*) 'Path loaded.'
  
end subroutine read_path


      !> Print the path to file. Ensembles correspond to images in GNEB method
subroutine write_en(n,x,y,dy,x0,filn,do_norm_rx)

      !.. Implicit declarations
      implicit none

      integer, intent(in) :: n     !< Number of samples
      real(kind=8), intent(in) :: x(n)        !< Reaction coordinate
      real(kind=8), intent(in) :: y(n)        !< Energy
      real(kind=8), intent(in) :: dy(n)       !< Derivative of the energy wrt x
      real(kind=8), intent(in) :: x0        !< Normalization for the reaction coordinate
      character(*), intent(in) :: filn             !< filename 
      character(len=1), intent(in) :: do_norm_rx   !< normalize reaction coordinate (Y/N)
      integer :: i
      real(kind=8) :: norm
      
         if (do_norm_rx=='Y') then
            norm = x0
         elseif (do_norm_rx=='N') then
            norm = 1d0
         else
            write(6,*) 'Invalid value for do_norm_rx!'
            STOP
            
         end if
         
         open(4, file=filn, access = 'sequential',action = 'write', status = 'replace')
    do i=1,n
          write (4,10002) x(i)/norm, y(i), dy(i)
    end do
    close(4)
    return
10002 format (es16.8E3,2x,es16.8E3,2x,es16.8E3)

  end subroutine write_en
   
   
subroutine spline_hermite_set ( ndata, tdata, c )
!

  implicit none
!
  integer, intent(in) :: ndata
!
  real(kind=8), intent(inout) :: c(4,ndata)
  real(kind=8) :: divdif1, divdif3,dt 
  integer i
  real(kind=8), intent(in) :: tdata(ndata)
!
  do i = 1, ndata-1
    dt = tdata(i+1) - tdata(i)
    divdif1 = ( c(1,i+1) - c(1,i) ) / dt
    divdif3 = c(2,i) + c(2,i+1) - 2d0* divdif1
    c(3,i) = ( divdif1 - c(2,i) - divdif3 ) / dt
    c(4,i) = divdif3 / (dt*dt)
  end do

  c(3,ndata) = 0d0
  c(4,ndata) = 0d0


end subroutine spline_hermite_set 

subroutine spline_hermite_val ( ndata, tdata, c, tval, sval, dsval )
!

!
  implicit none
!
  integer, intent(in) :: ndata
!
  real(kind=8), intent(in) :: c(4,ndata),tdata(ndata),tval
  real(kind=8) :: dt
  integer :: left
  real(kind=8), intent(out) :: sval, dsval

!
!  Find the interval [ TDATA(LEFT), TDATA(RIGHT) ] that contains
!  or is nearest to TVAL.
!
  call rvec_bracket ( ndata, tdata, tval, left)
!
!  Evaluate the cubic polynomial.
!
  dt = tval - tdata(left)

  sval = c(1,left) + dt * ( c(2,left) + dt * ( c(3,left) + dt * c(4,left) ) )
  dsval = c(2,left) + dt*(2d0*c(3,left) + dt*3d0*c(4,left))

  
end subroutine spline_hermite_val



subroutine rvec_bracket ( n, x, xval, left)

!
  implicit none
!
  integer, intent(in) :: n
   real(kind=8), intent(in) :: x(n),xval
	
  integer i
  integer, intent(out) :: left
  
!
  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      
      return
    end if

   end do

  left = n - 1
  


end subroutine rvec_bracket

   subroutine hermite_fit(n,nn,x,y,dydx,xx,yy,dyy,c)
      implicit none
      integer, intent(in) :: n,nn
      real(kind=8), intent(in) :: x(n),y(n),dydx(n)
      real(kind=8), intent(out) :: xx(nn),yy(nn),dyy(nn),c(4,n)
      
      real(kind=8) :: dx
      integer :: i
      
      do i =1,n
         c(1,i) = y(i)
         c(2,i) = dydx(i)
      end do
      
      dx = (x(n)-x(1))/(nn-1)
      
      call spline_hermite_set(n,x,c)
      
      do i=1,nn
         xx(i) = x(1) + (i-1)*dx
         call spline_hermite_val(n,x,c,xx(i),yy(i),dyy(i))
      end do
      
   
   end subroutine hermite_fit
   

end module m_gneb_utils
