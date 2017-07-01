      subroutine field_sd(j,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,h_ext)
      use m_fieldeff
      use m_parameters, only : i_DM,i_biq,i_four,i_dip,EA,Periodic_log
      use m_constants, only : pi
      use m_vector, only : norm,cross
      use m_dynamic, only : Ipol,torque_FL,torque_AFL
      implicit none

! Routine that calculates the Torque, the Energy and the force ensity fields
!
!
!
!
!

      integer, intent(in) :: j
! internals
      integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
      integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
      integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
      integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
      real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),h_ext(3)
      real(kind=8), dimension(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)) :: map_field,map_Tfield
      integer :: iomp,i_x,i_y,i_z,i_m,k,i
      character(len=30) :: fname,toto
!povray stuff
      real(kind=8) :: anglx,anglz,phi_color,widthc,Delta,Rc,Gc,Bc,dumy,max_force,max_Tforce
      real(kind=8) :: steptor(3),force(3),dmdr(3,3),f_in(3),f_ext(3),Torque_in(3),Torque_ext(3)
      integer :: ipu,ipv

      max_force=0.0d0
      widthc=5.0d0
      Delta =PI(2.0d0/3.0d0)
      map_field=0.0d0
      map_Tfield=0.0d0
      max_Tforce=0.0d0
      steptor=0.0d0

#ifndef CPP_OPENMP
      write(fname,'(I12)') j
      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'all_Torque_field_',(toto(i:i),i=1, &
         len_trim(toto)),'.dat'
      OPEN(72,FILE=fname,action='write',status='unknown',form='formatted')
      write(72,'(a)') '# x y z Sx Sy Sz Exch Zeeman Ani DM'
#endif

#ifdef CPP_OPENMP
!$OMP parallel DO default(shared) private(i_x,i_y,i_z,i_m)
#endif
       do i_m=1,shape_spin(5)
        do i_z=1,shape_spin(4)
         do i_y=1,shape_spin(3)
          do i_x=1,shape_spin(2)
        steptor=cross(spin(4:6,i_x,i_y,i_z,i_m),Ipol)

        map_field(1:3,i_x,i_y,i_z,i_m)=Bexch(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN) &
   &     +BZ(i_x,i_y,i_z,i_m,spin,shape_spin,h_ext) &
   &     +Bani(i_x,i_y,i_z,i_m,EA,spin,shape_spin,masque,shape_masque) &
   &     +efield(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
        if (i_DM) map_field(1:3,i_x,i_y,i_z,i_m)=map_field(1:3,i_x,i_y,i_z,i_m)+Bdm(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
        if (i_biq) map_field(1:3,i_x,i_y,i_z,i_m)=map_field(1:3,i_x,i_y,i_z,i_m)+Bbiqd(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)
        if (i_four) map_field(1:3,i_x,i_y,i_z,i_m)=map_field(1:3,i_x,i_y,i_z,i_m)+Bfour(i_x,i_y,i_z,i_m,spin,shape_spin,masque,shape_masque)
        if (i_dip) map_field(1:3,i_x,i_y,i_z,i_m)=map_field(1:3,i_x,i_y,i_z,i_m)+Bdip(i_x,i_y,i_z,i_m,spin,shape_spin)
        if (norm(map_field(1:3,i_x,i_y,i_z,i_m)).gt.max_force) max_force=norm(map_field(1:3,i_x,i_y,i_z,i_m))

           map_Tfield(1:3,i_x,i_y,i_z,i_m)=cross(map_field(1:3,i_x,i_y,i_z,i_m),spin(4:6,i_x,i_y,i_z,i_m))
        if (norm(map_Tfield(1:3,i_x,i_y,i_z,i_m)).gt.max_Tforce) max_force=norm(map_Tfield(1:3,i_x,i_y,i_z,i_m))

        write(72,'(26(E20.10,2x))') Spin(1:6,i_x,i_y,i_z,i_m),Bexch(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN) &
    &     ,BZ(i_x,i_y,i_z,i_m,spin,shape_spin,h_ext),Bani(i_x,i_y,i_z,i_m,EA,spin,shape_spin,masque,shape_masque) &
    &     ,Bdm(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN), &
    &      dot_product(Spin(4:6,i_x,i_y,i_z,i_m),Bexch(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)), &
    &      dot_product(Spin(4:6,i_x,i_y,i_z,i_m),BZ(i_x,i_y,i_z,i_m,spin,shape_spin,h_ext)),dot_product(Spin(4:6,i_x,i_y,i_z,i_m),Bani(i_x,i_y,i_z,i_m,EA,spin,shape_spin,masque,shape_masque)), &
    &      dot_product(Spin(4:6,i_x,i_y,i_z,i_m),Bdm(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN)),torque_FL*steptor, &
    &      torque_FL*torque_AFL*cross(steptor,spin(4:6,i_x,i_y,i_z,i_m))

          enddo
         enddo
        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif

#ifndef CPP_OPENMP
       close(72)
#endif

#ifndef CPP_OPENMP
      write(fname,'(I12)') j
      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'all_Force_field_',(toto(i:i),i=1, &
         len_trim(toto)),'.dat'
      OPEN(73,FILE=fname,action='write',status='unknown',form='formatted')
      write(73,'(a)') '# x y z Sx Sy Sz F_in_u F_in_v F_in_uv F_ext_u F_ext_v F_ext_uv F_u F_v F_uv'
#endif

#ifdef CPP_OPENMP
!$OMP parallel DO default(shared) private(i_x,i_y,i_z,i_m)
#endif
       do i_m=1,shape_spin(5)
        do i_z=1,shape_spin(4)
         do i_y=1,shape_spin(3)
          do i_x=1,shape_spin(2)

          force=0.0d0
          dmdr=0.0d0
          f_in=0.0d0
          f_ext=0.0d0
          Torque_in=0.0d0
          Torque_ext=0.0d0

        if (Periodic_log(1).and.Periodic_log(2)) then
         ipu=mod(i_x+i-1+shape_spin(2),shape_spin(2))+1
         ipv=mod(i_y+j-1+shape_spin(3),shape_spin(3))+1
        else
         if (i_x+i.eq.shape_spin(2)) then
          ipu=i_x
          else
          ipu=i_x+i
         endif
         if (i_y+j.eq.shape_spin(3)) then
          ipv=i_y
          else
          ipv=i_y+1
         endif
        endif

! calculations of the Wronskian

        dmdr(:,1)=spin(4:6,ipu,i_y,i_z,i_m)-(dot_product(spin(4:6,ipu,i_y,i_z,i_m),spin(4:6,i_x,i_y,i_z,i_m))*spin(4:6,i_x,i_y,i_z,i_m))
        dmdr(:,2)=spin(4:6,i_x,ipv,i_z,i_m)-(dot_product(spin(4:6,i_x,ipv,i_z,i_m),spin(4:6,i_x,i_y,i_z,i_m))*spin(4:6,i_x,i_y,i_z,i_m))
        dmdr(:,3)=spin(4:6,ipu,ipv,i_z,i_m)-(dot_product(spin(4:6,ipu,ipv,i_z,i_m),spin(4:6,i_x,i_y,i_z,i_m))*spin(4:6,i_x,i_y,i_z,i_m))

! calculations on the internal forces (Exch, ani,Dm...)

        Torque_in=map_field(1:3,i_x,i_y,i_z,i_m)-BZ(i_x,i_y,i_z,i_m,spin,shape_spin,h_ext)


        f_in(1)=dot_product(Torque_in,dmdr(:,1))
        f_in(2)=dot_product(Torque_in,dmdr(:,2))
        f_in(3)=dot_product(Torque_in,dmdr(:,3))

! calculations on the external forces (Zeeman, Torques...)

        Torque_ext=BZ(i_x,i_y,i_z,i_m,spin,shape_spin,h_ext)
        Torque_ext=Torque_ext+cross(spin(4:6,i_x,i_y,i_z,i_m),Ipol)

        f_ext(1)=dot_product(Torque_ext,dmdr(:,1))
        f_ext(2)=dot_product(Torque_ext,dmdr(:,2))
        f_ext(3)=dot_product(Torque_ext,dmdr(:,3))

! total forces

        force=f_in+f_ext

        write(73,'(15(E20.10,2x))') Spin(1:6,i_x,i_y,i_z,i_m),-f_in,-f_ext,-force

          enddo
         enddo
        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif

#ifndef CPP_OPENMP
       close(73)
#endif

! find the max of the force

      write(fname,'(I12)') j
      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'Ffield_',(toto(i:i),i=1, &
         len_trim(toto)),'.dat'
      OPEN(70,FILE=fname,action='write',status='unknown',form='formatted')

      write(fname,'(I12)') j
      toto=trim(adjustl(fname))
      write(fname,'(a,18a,a)')'Tfield_',(toto(i:i),i=1, &
         len_trim(toto)),'.dat'
      OPEN(71,FILE=fname,action='write',status='unknown',form='formatted')

        do i_m=1,shape_spin(5)
         do i_z=1,shape_spin(4)
          do i_y=1,shape_spin(3)
           do i_x=1,shape_spin(2)


        if (norm(map_field(:,i_x,i_y,i_z,i_m)).gt.(1.0d-8*max_force)) then
         dumy=map_field(3,i_x,i_y,i_z,i_m)/norm(map_field(:,i_x,i_y,i_z,i_m))
         anglz=acos(dumy)*180.0d0/pi(1.0d0)
        else
         anglz=0.0d0
        endif

        anglx=atan2(map_field(2,i_x,i_y,i_z,i_m),map_field(1,i_x,i_y,i_z,i_m))
        anglx=anglx*180.0d0/pi(1.0d0)

        phi_color=pi(anglx/300.0d0*2.0d0)
        Rc = widthc*(cos(phi_color+0*Delta))
        if (Rc.lt.0.000001d0)  Rc=0.0d0
        Gc = widthc*(cos(phi_color+1*Delta))
        if (Gc.lt.0.000001d0)  Gc=0.0d0
        Bc = widthc*(cos(phi_color+2*Delta))
        if (Bc.lt.0.000001d0)  Bc=0.0d0

        write(70,'(9(a,f16.8),a)') 'Vector(',norm(map_field(1:3,i_x,i_y,i_z,i_m))/max_force,',', &
     & anglz,',',anglx,',',Spin(1,i_x,i_y,i_z,i_m),',',Spin(2,i_x,i_y,i_z,i_m),',',Spin(3,i_x,i_y,i_z,i_m),',', &
     & Rc,',',Bc,',',Gc,')'

!         Write(70,'(9(f14.8,2x))') (Spin(k,i_x,i_y,i_z,i_m),k=1,3),(map_field(k,i_x,i_y,i_z,i_m),k=1,3)

        if (norm(map_Tfield(:,i_x,i_y,i_z,i_m)).gt.1.0d-8) then
         dumy=map_Tfield(3,i_x,i_y,i_z,i_m)/norm(map_Tfield(:,i_x,i_y,i_z,i_m))
         anglz=acos(dumy)*180.0d0/pi(1.0d0)
        else
         anglz=0.0d0
        endif

        anglx=atan2(map_Tfield(2,i_x,i_y,i_z,i_m),map_Tfield(1,i_x,i_y,i_z,i_m))
        anglx=anglx*180.0d0/pi(1.0d0)

        phi_color=pi(anglx/300.0d0*2.0d0)
        Rc = widthc*(cos(phi_color+0*Delta))
        if (Rc.lt.0.000001d0)  Rc=0.0d0
        Gc = widthc*(cos(phi_color+1*Delta))
        if (Gc.lt.0.000001d0)  Gc=0.0d0
        Bc = widthc*(cos(phi_color+2*Delta))
        if (Bc.lt.0.000001d0)  Bc=0.0d0

        write(71,'(9(a,f16.8),a)') 'Vector(',norm(map_Tfield(1:3,i_x,i_y,i_z,i_m))/max_force,',', &
     & anglz,',',anglx,',',Spin(1,i_x,i_y,i_z,i_m),',',Spin(2,i_x,i_y,i_z,i_m),',',Spin(3,i_x,i_y,i_z,i_m),',', &
     & Rc,',',Bc,',',Gc,')'

!         Write(71,'(9(f14.8,2x))') (Spin(k,i_x,i_y,i_z,i_m),k=1,3),(map_Tfield(k,i_x,i_y,i_z,i_m),k=1,3)

           enddo
          enddo
         enddo
        enddo
       close(71)
       close(70)

      end subroutine field_sd
