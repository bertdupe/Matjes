subroutine GNEB(state,i_biq,i_dm,i_four,i_dip,EA,h_ext, &
    & indexNN,shape_index,shape_spin,tableNN,shape_tableNN,masque,shape_masque,N_cell)!,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,state)
      use m_fieldeff
      use m_energy
      use mtprng
      use m_gneb_utils
      use m_write_spin
      use m_createspinfile
      implicit none
      logical, intent(in) :: i_biq,i_dm,i_four,i_dip
      real(kind=8), intent(in) :: EA(3),h_ext(3)
      integer, intent(in) :: shape_spin(5),shape_index(2),shape_tableNN(6),shape_masque(4),N_cell
      type(mtprng_state), intent(inout) :: state
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4)),indexNN(shape_index(1),shape_index(2))
      integer :: iseed,clock_int,i
      character(9) :: clock
      real(kind=8) :: spini(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8) :: spinf(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8) :: spinsp(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
      real(kind=8), allocatable :: path(:,:,:,:,:,:),rx(:),ene(:),dene(:),xx(:),yy(:),dyy(:),c(:,:)
      real(kind=8) :: rx0(1),drx0,ene0(1),dene0(1)
      integer :: imax,i_x,i_y,i_z,i_m,ci
      logical :: exists
      character(len=8) :: num

      
      spini=0.0d0
      spinf=0.0d0
      spinsp=0.0d0

!      call SYSTEM_CLOCK(COUNT=clock_int)
!      write(clock,'(I9)') clock_int
!      read(clock(6:9),*)  iseed
!      call mtprng_init(iseed,state)

      call set_gneb_defaults()
      call read_gneb_parameters()
      
      allocate(path(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5),nim),rx(nim),ene(nim),dene(nim))
      
      path=0.0d0

      print *,nim
      
      
      
      
      if (initpath==2) then
         write (*,'(a)') "Initial guess for the path from file "
         call read_path(state,trim(adjustl(restartfile_path)),nim,amp_rnd_path,shape_spin,path,exists)
         if (exists) then
            write(num,'(I8)') nim
            momfile_i = trim(adjustl(restartfile_path))//'_1.dat'
            momfile_f = trim(adjustl(restartfile_path))//'_'//trim(adjustl(num))//'.dat'
            call read_inifin(state,trim(adjustl(momfile_i)),trim(adjustl(momfile_f)),amp_rnd,shape_spin,spini,spinf)
            call WriteSpinAndCorrFile('SpinSTM_GNEB_ini.dat',spini,shape_spin)
            call CreateSpinFile('povray_GNEB_ini.dat',spini,shape_spin)
            call WriteSpinAndCorrFile('SpinSTM_GNEB_fin.dat',spinf,shape_spin)
            call CreateSpinFile('povray_GNEB_fin.dat',spinf,shape_spin)
            
            write (*,'(a)') "Relaxing the first image..."
     
            call minimize(i_biq,i_dm,i_four,i_dip,.true.,mintraj_step,EA, &
              & spini,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,N_cell,h_ext)
            write (*,'(a)') "Done!"
            
            write (*,'(a)') "Relaxing the last image..."
            call minimize(i_biq,i_dm,i_four,i_dip,.true.,mintraj_step,EA, &
              & spinf,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,N_cell,h_ext)
            write (*,'(a)') "Done!"
         
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
                        do i=1,shape_spin(1)
                           path(i,i_x,i_y,i_z,i_m,1) = spini(i,i_x,i_y,i_z,i_m)
                           path(i,i_x,i_y,i_z,i_m,nim) = spinf(i,i_x,i_y,i_z,i_m)
                        end do
                     end do
                  end do
               end do
            end do
         
         else
            call read_inifin(state,trim(adjustl(momfile_i)),trim(adjustl(momfile_f)),amp_rnd,shape_spin,spini,spinf)
            call WriteSpinAndCorrFile('SpinSTM_GNEB_ini.dat',spini,shape_spin)
            call CreateSpinFile('povray_GNEB_ini.dat',spini,shape_spin)
            call WriteSpinAndCorrFile('SpinSTM_GNEB_fin.dat',spinf,shape_spin)
            call CreateSpinFile('povray_GNEB_fin.dat',spinf,shape_spin)
            
            write (*,'(a)') "Relaxing the first image..."
            call minimize(i_biq,i_dm,i_four,i_dip,.true.,mintraj_step,EA, &
              & spini,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,N_cell,h_ext)
            write (*,'(a)') "Done!"
            
            write (*,'(a)') "Relaxing the last image..."
            call minimize(i_biq,i_dm,i_four,i_dip,.true.,mintraj_step,EA, &
              & spinf,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,N_cell,h_ext)
            write (*,'(a)') "Done!"
            
            
            call geodesic_path(state,nim,amp_rnd_path,spini,spinf,shape_spin,path)
            write (*,'(a)') "Geodesic path generated"
         end if
      else
         call read_inifin(state,trim(adjustl(momfile_i)),trim(adjustl(momfile_f)),amp_rnd,shape_spin,spini,spinf)
         call WriteSpinAndCorrFile('SpinSTM_GNEB_ini.dat',spini,shape_spin)
         call CreateSpinFile('povray_GNEB_ini.dat',spini,shape_spin)
         call WriteSpinAndCorrFile('SpinSTM_GNEB_fin.dat',spinf,shape_spin)
         call CreateSpinFile('povray_GNEB_fin.dat',spinf,shape_spin)
         
         write (*,'(a)') "Relaxing the first image..."

         call minimize(i_biq,i_dm,i_four,i_dip,.true.,mintraj_step,EA, &
              & spini,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,N_cell,h_ext)
         write (*,'(a)') "Done!"
              
         write (*,'(a)') "Relaxing the last image..."
         call minimize(i_biq,i_dm,i_four,i_dip,.true.,mintraj_step,EA, &
              & spinf,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,N_cell,h_ext)
         write (*,'(a)') "Done!"
      
         call geodesic_path(state,nim,amp_rnd_path,spini,spinf,shape_spin,path)
         write (*,'(a)') "Geodesic path generated"
      end if
      
      call write_path(nim,shape_spin,path)
      
      if (do_gneb=='Y') then
         write (*,'(a)') "GNEB calculation in progress..."
         call find_path(i_DM,i_four,i_biq,i_dip,EA,nim,vpodt,vpomass,spring,path, &
                  & shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,mepftol,mepitrmax,meptraj_step,h_ext, &
                  & rx,ene,dene)
         
         write (*,'(a)') "Done!"
      end if
      
      if (do_gneb_ci=='Y') then
         write (*,'(a)') "CI-GNEB calculation in progress..."
         call find_path_ci(i_DM,i_four,i_biq,i_dip,EA,nim,vpodt,vpomass,spring,path, &
                  & shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,mepftol_ci,mepitrmax,meptraj_step,h_ext, &
                  & rx,ene,dene,ci)
         print *,'ci:',ci
         write (*,'(a)') "Done!"
      end if
      
      call write_en(nim,rx,ene,dene,rx(nim),'en_path.out',do_norm_rx)
      call write_path(nim,shape_spin,path)
      
      allocate(xx(sample_num),yy(sample_num),dyy(sample_num),c(4,nim))
      
      call hermite_fit(nim,sample_num,rx,ene,dene,xx,yy,dyy,c)
      

      
      call write_en(sample_num,xx,yy,dyy,rx(nim),'enfit_path.out',do_norm_rx)
      
      
      if (do_gneb_ci=='Y') then 
         rx0(1) = rx(ci)
         do i_m=1,shape_spin(5)
            do i_z=1,shape_spin(4)
               do i_y=1,shape_spin(3)
                  do i_x=1,shape_spin(2)
                     do i=1,shape_spin(1)
                        spinsp(i,i_x,i_y,i_z,i_m) = path(i,i_x,i_y,i_z,i_m,ci)
                     end do
                  end do
               end do
            end do
         end do
      else
      
         call find_SP(nim,rx,c,imax,drx0)
         rx0(1) = rx(imax)+drx0
         if (imax==nim) then
            do i_m=1,shape_spin(5)
               do i_z=1,shape_spin(4)
                  do i_y=1,shape_spin(3)
                     do i_x=1,shape_spin(2)
                        do i=1,shape_spin(1)
                           spinsp(i,i_x,i_y,i_z,i_m) = path(i,i_x,i_y,i_z,i_m,imax)
                        end do
                     end do
                  end do
               end do
            end do
         else
            call find_SP_conf(state,shape_spin,path(:,:,:,:,:,imax),path(:,:,:,:,:,imax+1),rx(imax),rx(imax+1),rx(imax)+drx0,spinsp)
         end if
      end if
      
      
      
         
      call spline_hermite_val(nim,rx,c,rx0(1),ene0(1),dene0(1))
      
      call write_en(1,rx0,ene0,dene0,rx(nim),'ensp_path.out',do_norm_rx)
      
      
      
      call WriteSpinAndCorrFile('image-GNEB-saddlepoint.dat',spinsp,shape_spin)
      call CreateSpinFile('povray-GNEB-saddlepoint.dat',spinsp,shape_spin)
      
      

      deallocate(path,rx,ene,dene,xx,yy,dyy,c)
end subroutine GNEB 
