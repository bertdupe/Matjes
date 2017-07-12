	  subroutine print_line()
	  write(6,*) '#####################################################'
	  end subroutine



! ###################################################################
! ################ Beginn des Programms: ############################
! ###################################################################

	  subroutine spstm
      use m_lattice, only : spin
      use m_rw_lattice, only : dim_lat,motif,net
      use m_parameters, only : spstmonly,Periodic_log
      use m_vector, only : norm
      use m_constants, only : pi
      use m_sym_utils, only : order_zaxis,pos_nei
      implicit none
      integer :: i,j,k,h,xmax,ymax,n_atoms,stat,x_index,y_index
      integer :: i_x,i_y,i_z,i_m
      integer :: n_raster,file_index,spin_control,y_number,lauf
      real(kind=8) :: xstart,xende,ystart,yende,x,y,maximum,minimum
      real(kind=8), dimension(:,:), allocatable :: daten
      real(kind=8), dimension(:,:), allocatable :: erg,field
      real(kind=8), dimension(:,:), allocatable :: image
      real(kind=8), dimension(:,:), allocatable :: inbox,fft
      real(kind=8), dimension(28) :: config
      real(kind=8) :: height,tipx,tipy,r,x_atom,y_atom,p_eff,tipz
      real(kind=8) :: kappa,cosinus,cos_theta,vec_x1,vec_x2,vec_y1,vec_y2
	  integer :: x_units,y_units,n_spins,test,maximal,cut
	  real(kind=8), dimension(:,:), allocatable :: spins
	  real(kind=8) :: b3,a3,tamr,gamm,x_length,y_length,verh
	  character(len=100) :: datei,ratio,si1,si2,si,xlength,ylength
	  character(len=100) :: bound,over1,over,max_array,blaaa
	  character(len=3) :: nummer
	  character :: komma
      character(len=10**5) :: string,punkte,pov
	  real(kind=8) :: max_z,min_z,abs_z,vec_len,latt,lx,ly,lz,dir1(3),dir2(3)
	  real(kind=8) :: orbital,ys,ye,xs,xe,xl,yl,betr,x_scan,y_scan
	  real(kind=8), dimension(:), allocatable :: corrugation,dat,speicher1
	  real(kind=8) :: speicher,cos_phi,he
	  integer :: x_ind,y_ind,orb_choice
	  integer, dimension(7) :: elemente
	  integer, dimension(2) :: nn
	  character(len=8) :: spx,spy,spz,lax,lay,xposi,yposi,ro,gr,bl,zpo
	  character(len=8) :: laz,xrange,yrange,hx,hy
	  real(kind=8) :: spinx,spiny,spinz,zposi,povscale,yscan,xscan
      real(kind=8) :: xmed,ymed,dx,dy,vec_test(2),vec_norm,cutdown
      logical :: i_dir,i_old,i_user_spin
!     colors
      real(kind=8) :: phi_color, Delta, widthc
      real(kind=8) :: Rc,Gc,Bc
!     calculating the angles
      real(kind=8) :: angly,anglz,orderz
! scal line part
      real(kind=8),allocatable :: radius(:,:,:),sample(:,:)
      real(kind=8),allocatable :: vec_test_dist(:)
      integer :: erreur,n_count
! taking care of the periodic boundary conditions
      real(kind=8) :: x_n,y_n
! coponent of the fields
      real(kind=8) :: fieldx,fieldy,fieldz

	  i_dir=.False.
	  i_old=.False.
	  i_user_spin=.False.
! ###################################################################
! ################ Vorgabe der Rahmenbedingungen: ###################
! ###################################################################


      orderz=abs(dble(order_zaxis(net)))
      if (orderz.eq.4.0d0) orderz=2.0d0
      if (orderz.eq.6.0d0) orderz=3.0d0

! Einlesen der Konfigurationsdatei 'config.dat'
	  open(38,file='config.dat')
	  do i=1,28
	  read(38,*,iostat=stat) config(i)
	  if (stat .ne. 0) then
	  call print_line()
	  write(6,'(a)') 'Error reading "config.dat"'
	  call print_line()
	  stop
	  endif
	  enddo
	  close(38)

	  inquire (file='direction',exist=i_dir)
	  if (i_dir) then
	   open(38,file='direction')
	   read(38,*) (dir1(i),i=1,3)
	   read(38,*) (dir2(i),i=1,3)
	  endif

      inquire (file='user_def_input.dat',exist=i_user_spin)

! Definition der Rahmenbedingungen aus 'config.dat'
	  tipx=config(1)
	  tipy=config(2)
	  tipz=config(3)
	  vec_norm=sqrt(sum(config**2))
	  if (vec_norm.gt.1.0d-5) then
	   tipx=tipx/vec_norm
	   tipy=tipy/vec_norm
	   tipz=tipz/vec_norm
	  endif

	  kappa=config(4)! Kappa in 1/Angstrom
	  height=config(5)  ! Abstand zwischen Tip und Probe in Angstrom

	  xmax=2**config(6)
	  ymax=2**config(7)! Auflösung des betrachteten Bereichs ergibt sich zu xmax*ymax


       xstart=config(8)
       xende=config(9)
       ystart=config(10)
       yende=config(11)

        xmed=(xstart+xende)/2.0d0
        ymed=(ystart+yende)/2.0d0
        dx=abs(xstart-xende)/2.0d0
        dy=abs(ystart-yende)/2.0d0

       write(6,*) "you are using ",xmax," , ",ymax," points in the x,y direction" 
       write(6,*) "density of points in the x direction: ", xmax/dx
       write(6,*) "density of points in the y direction: ", ymax/dy

      x_units=config(12) ! Anzahl der auszugebenden Bildausschnitte in x-Richtung
	  y_units=config(13) ! """""  y-Richtung

 	  P_eff=config(14)
	
      gamm=config(15)
      lx=config(16)
      ly=config(17)
      lz=config(18)
      latt=sqrt(lx**2+ly**2+lz**2)

      maximal=int(config(19))  ! Maximale Höhe/Breite der .png Bilder

      vec_len=config(20) ! Skalierungsfaktor für den Spinvektorplot
      orb_choice=config(21) ! Auswahl des Tip-Orbitals

      xs=config(22) !Start- und Endpunkt der Scanline
      ys=config(23)
      xe=config(24)
      ye=config(25)
      cut=config(27)
      cutdown=config(28)
	
      xscan=xe-xs
      yscan=ye-ys

	
      n_raster=0  ! Anzahl der Atome innerhalb des abgerasterten Bereiches




! ###################################################################
! ################ Einlesen der Oberfläche ##########################
! ###################################################################
!------------------------------------------------------------------------
! debut du if portant l utilisation du programme pour tracer uniquement l image STM
! donnee dans input.dat
! Bestimmung der Anzahl gegebener Atome:
      if (spstmonly) then
      open(42,file='input.dat')
      n_atoms=0
       do
        read(42, '(A)', iostat=stat)
        if (stat /= 0) exit
        n_atoms=n_atoms+1
        end do 
	  close(42)

      if(n_atoms .eq. 0) then
      call print_line()
      write(6,'(a)') 'File "input.dat" is empty.'
      call print_line()
      endif

! Ausgabe der Anzahl eingelesener Atome
      call print_line
        write(6,'(a,I10)') 'number of given atoms:',n_atoms
      call print_line

! Initialisierung der Datenfelder:

      allocate(daten(n_atoms,6))
!	allocate(spins(n_atoms,3))

	spin_control=0 ! Überprüfung, ob die Spinkonfiguration
			   ! in 2D oder 3D geplottet werden soll.
			   ! Falls es verschiedene z-Komponenten gibt, setze spin_control 1

! Auslesen von 'input.dat' in daten(:,:) und Erstellung der Liste spins(:,:) mit
! allen unterschiedlichen Spin-Konfigurationen
! if the input.dat contains only 5 values then it is the old format
      daten=1.0d0
      open(41,file='input.dat')
      inquire(file='old-format',exist=i_old)
      if (i_old) then
       do i=1,n_atoms
        read(41,*) (daten(i,j),j=1,2),(daten(i,j),j=4,6)
       enddo
      else
       do i=1,n_atoms
        read(41,*) (daten(i,j),j=1,6)
       enddo
      endif
      close(41)
! prend les spins non pas dans input.dat mais comme sortie du code MC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
      else
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      allocate(daten(product(dim_lat),6))
      daten=0.0d0

      n_atoms=product(dim_lat)
      i_m=size(spin,5)
      do i_x=1,dim_lat(1)
       do i_y=1,dim_lat(2)
        do i_z=1,dim_lat(3)
        i=i_x+(i_y-1)*dim_lat(1)+(i_z-1)*dim_lat(1)*dim_lat(2)
        daten(i,:)=spin(1:6,i_x,i_y,i_z,i_m)
        enddo
       enddo
      enddo
      endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc

	   if (daten(1,6).ne.0.0) spin_control=1

! count the number of atoms that should be included in the domain
! check that spins are in the box
      if (all(periodic_log)) then
       n_raster=n_atoms
       allocate(spins(n_atoms,3))
       allocate(inbox(n_atoms,2))
       inbox(1:n_atoms,1)=daten(:,1)
       inbox(1:n_atoms,2)=daten(:,2)
       spins(1:n_atoms,:)=daten(:,4:6)
      else
       n_raster=0
       do i=1,n_atoms
        if ((daten(i,1).lt.xende).and.(daten(i,1).gt.xstart).and.(daten(i,2).gt.ystart).and.(daten(i,2).lt.yende)) n_raster=n_raster+1
       enddo

       allocate(spins(n_raster,3))
       allocate(inbox(n_raster,2))
       n_raster=0
       do i=1,n_atoms
        if ((daten(i,1).lt.xende).and.(daten(i,1).gt.xstart).and.(daten(i,2).gt.ystart).and.(daten(i,2).lt.yende)) then
         n_raster=n_raster+1
         inbox(n_raster,:)=daten(n_raster,1:2)
         spins(n_raster,:)=daten(n_raster,4:6)
        endif
       enddo
      n_atoms=n_raster
      endif

      n_spins=1






! ################################################################################
! ################# Voreinstellungen zur Bildausgabe #############################
! ################################################################################


! Festlegen des zu plottenden Bildausschnitts:

	if (dabs((xende-xstart)*x_units) .lt. 10) then
	write(xlength,'(I2)') int((xende-xstart)*x_units)
	elseif (dabs((xende-xstart)*x_units) .lt. 100) then
	write(xlength,'(I3)') int((xende-xstart)*x_units)
	elseif (dabs((xende-xstart)*x_units) .lt. 1000) then
	write(xlength,'(I4)') int((xende-xstart)*x_units)
	elseif (dabs((xende-xstart)*x_units) .lt. 10000) then
	write(xlength,'(I5)') int((xende-xstart)*x_units)
	else
	write(xlength,'(I6)') int((xende-xstart)*x_units)
	endif


       
	if (dabs((yende-ystart)*y_units) .lt. 10) then
	write(ylength,'(I2)') int((yende-ystart)*y_units)
	elseif (dabs((yende-ystart)*y_units) .lt. 100) then
	write(ylength,'(I3)') int((yende-ystart)*y_units)
	elseif (dabs((yende-ystart)*y_units) .lt. 1000) then
	write(ylength,'(I4)') int((yende-ystart)*y_units)
	elseif (dabs((yende-ystart)*y_units) .lt. 10000) then
	write(ylength,'(I5)') int((yende-ystart)*y_units)
	else
	write(ylength,'(I6)') int((yende-ystart)*y_units)
	endif


       verh=dabs((yende-ystart)*y_units/((xende-xstart)*x_units))
	write(ratio,*) verh

! Berechne Größe der *.png Bilder:
       if(verh .lt. 1.0) then
       x_length=maximal
       y_length=x_length*verh
       else
       y_length=maximal
       x_length=y_length/verh
       endif

       if (dabs(x_length) .lt. 10) then
       write(si1,'(I2)') int(x_length)
       elseif(x_length .lt. 100) then
       write(si1,'(I3)') int(x_length)
       elseif(x_length .lt. 1000) then
       write(si1,'(I4)') int(x_length)
       elseif(x_length .lt. 10000) then
       write(si1,'(I5)') int(x_length)
       else
       write(si1,'(I6)') int(x_length)
       endif

       if (dabs(y_length) .lt. 10) then
       write(si2,'(I2)') int(y_length)
       elseif(y_length .lt. 100) then
       write(si2,'(I3)') int(y_length)
       elseif(y_length .lt. 1000) then
       write(si2,'(I4)') int(y_length)
       elseif(y_length .lt. 10000) then
       write(si2,'(I5)') int(y_length)
       else
       write(si2,'(I6)') int(y_length)
       endif

       bound='[0:'//trim(xlength)//'] [0:'//trim(ylength)//']'




! ############################################################################
! ########## POV-Ray Ausgabe der eingegebenen Spinkonfiguration ##############
! ############################################################################
! #                                                                          #
! #                         http://www.povray.org/                           #
! #                                                                          #
! ############################################################################


	 povscale=0.8*sqrt(dabs((yende-ystart)*(xende-xstart))/n_raster)

	 zposi=max((xende-xstart)*dble(x_units),(yende-ystart)*dble(y_units))*1.5
!        zposi=max(xmed*x_units,ymed*y_units)*1.5
	 write(zpo,'(I8)') int(zposi)

	 write(xposi,'(I8)') int((xende+xstart)/2.0*x_units)
	 write(yposi,'(I8)') int((yende+ystart)/2.0*y_units)

	 open(78,file='3d_spins.pov')
	 write(78,'(a)') '#include "colors.inc"'
	 write(78,'(a)') '#include "stones.inc"'
	 write(78,'(a)') 'background{White}'

	 write(78,*) 'camera {'
	 write(78,*) 'location <'//trim(xposi)//','//trim(yposi)//','//trim(zpo)//'>'
	 write(78,*) 'look_at <'//trim(xposi)//','//trim(yposi)//', 0>'
	 write(78,*) 'up <1,0,0>'
	 write(78,*) 'right <1,0,0>'
	 write(78,*) '}'


	 do k=1,n_raster

	  write(78,*) 'cone {'
	  write(78,*) '<0.0,0.0,1.0>, 0.0'
	  write(78,*) '<0.0,0.0,-1.0>, 0.18'

        if (dabs(Spins(k,3)).lt.1.0d0) then
          angly=acos(Spins(k,3))*180.0d0/pi(1.0d0)
        else
          angly=90.0d0-dsign(90.0d0,Spins(k,3))
        endif

        anglz=atan2(Spins(k,2),Spins(k,1))
        anglz=anglz*180.0d0/pi(1.0d0)

!       Calcualting the color as a function of the angle in or
!       out of the plane
        phi_color=pi(angly/300.0d0*2.0d0)
        Rc = 5.0d0*(cos(phi_color+0*PI(2.0d0/3.0d0)))
        if (Rc.lt.0.000001d0)  Rc=0.0d0
        Gc = 5.0d0*(cos(phi_color+1*PI(2.0d0/3.0d0)))
        if (Gc.lt.0.000001d0)  Gc=0.0d0
        Bc = 5.0d0*(cos(phi_color+2*PI(2.0d0/3.0d0)))
        if (Bc.lt.0.000001d0)  Bc=0.0d0

      write(78,*) 'rotate   < 0.0 ,' ,angly, ' , 0.0 >'
      write(78,*) 'rotate   < 0.0 , 0.0 ,',anglz,' >'
      write(78,*) 'translate   < ',inbox(k,1),' , ',inbox(k,2),' , 0.0 >'
      write(78,*) 'texture{pigment  {color rgb <',Rc,',',Bc,',',Gc,'>}}}'

	enddo

	write(78,*) 'light_source { <'//trim(xposi)//','//trim(yposi)//', 50> color White}'

    close(78)

! #########################################################################
! ########### Ende der Pov-Ray Ausgabe ####################################
! #########################################################################











! Ausgabe von Atomen/ Spins innerhalb des betrachteten Bereichs:

	call print_line()
	write(6,*) 'number of spins inside the box:',n_raster
	call print_line()


	allocate(erg(xmax,ymax))
	allocate(field(xmax,ymax))

	erg=0.0d0
	field=0.0d0

! ###################################################################
! ############# Beginn der Rechnung #################################
! ###################################################################
!debut du calcul de la densite d etats dans le vide

! Abrastern des ausgewählten Bereiches:

    allocate(radius(3,xmax,ymax))
    allocate(sample(xmax,ymax))
    allocate(vec_test_dist(n_raster))
    if (.not.i_dir) then
     dir1=net(1,:)
     dir2=net(2,:)
    endif
    sample=0.0d0

    if (abs(net(1,2)*net(2,1)-net(1,1)*net(2,2)).lt.1.0d-8) then
     write(*,*) 'error in the input lattice'
     stop
    endif

	do i=1,xmax
	 do j=1,ymax

      radius(:,i,j)=(ystart+(yende-ystart)*dble(j-1)/dble(ymax))*dir2+ &
       (xstart+(xende-xstart)*dble(i-1)/dble(xmax))*dir1

      if (all(periodic_log)) then
! if I am periodic, replace the tip in the sample
!we are looking for x_n, y_n, position of the tip in the simulation box

       y_n=(net(1,2)*radius(1,i,j)-net(1,1)*radius(2,i,j))/(net(1,2)*net(2,1)-net(1,1)*net(2,2))
       x_n=(radius(1,i,j)-y_n*net(2,1))/net(1,1)

       if (i_user_spin) then
        if ((x_n.ge.maxval(inbox(:,1),1)).or.(x_n.lt.0.0d0)) x_n=x_n-floor(x_n/maxval(inbox(:,1),1))*maxval(inbox(:,1),1)
        if ((y_n.ge.maxval(inbox(:,2),1)).or.(y_n.lt.0.0d0)) y_n=y_n-floor(y_n/maxval(inbox(:,2),1))*maxval(inbox(:,2),1)
       else
        if ((x_n.ge.dble(dim_lat(1))).or.(x_n.lt.0.0d0)) x_n=x_n-floor(x_n/dble(dim_lat(1)))*dble(dim_lat(1))
        if ((y_n.ge.dble(dim_lat(2))).or.(y_n.lt.0.0d0)) y_n=y_n-floor(y_n/dble(dim_lat(2)))*dble(dim_lat(2))
       endif

       radius(:,i,j)=x_n*net(1,:)+y_n*net(2,:)

      endif

     enddo
    enddo

    do i=1,xmax
     do j=1,ymax

     vec_test_dist=1000.0d0
         vec_test=0.0d0
         n_count=0
! F�ur jeden der xmax*ymax Punkte wird nun der Tunnelstrombeitrag
! aller n_atoms eingegeben Atome aufsummiert:

	  do k=1,n_raster

!	inbox(k,1) ! x-Koord des k-ten Atoms
!	inbox(k,2) ! y_koord """"""""""""

       if (all(periodic_log)) then
        if (i_user_spin) then
! case of an input given by a colleague and the code can not handle the lattice well
         vec_test=pos_nei(radius(1:2,i,j),inbox(k,:),maxval(inbox(:,1),1),maxval(inbox(:,2),1))
         r=sqrt(vec_test(1)**2+vec_test(2)**2)
        else
         vec_test=pos_nei(radius(1:2,i,j),inbox(k,:),net,dim_lat)
         r=sqrt(vec_test(1)**2+vec_test(2)**2) ! Abstand zw. Tip und k-tem Atom
        endif
       else
        vec_test=radius(1:2,i,j)-inbox(k,1:2)
	    r=sqrt((radius(1,i,j)-inbox(k,1))**2+(radius(2,i,j)-inbox(k,2))**2) ! Abstand zw. Tip und k-tem Atom
       endif

	  if (r .lt. cut) then

        n_count=n_count+1
        vec_test_dist(n_count)=sqrt(r**2+height**2)

	    cos_theta=cosinus(tipx,tipy,tipz,spins(k,1),spins(k,2),spins(k,3)) ! Winkel zwischen ...
      ! ... Tip-Magnetisierung und Magnetisierung des k-ten Atoms

! #########################################################################
! Auswahl des Spitzen-Orbitals:
      select case(orb_choice)
	  case(1)
	   orbital=1.0 ! s Orbital
	  case(2)
       orbital=(kappa/vec_test_dist(n_count))**2*(vec_test(1))**2                                !p_x Orbital
	  case(3)
	   orbital=(kappa/vec_test_dist(n_count))**2*(vec_test(2))**2                                !p_y Orbital
	  case(4)
	   orbital=(kappa/vec_test_dist(n_count))**2*(height)**2                                  !p_z Orbital
	  case(5)
	   orbital=kappa**2*(vec_test(1)**2*(1.0/vec_test_dist(n_count)**3+kappa/vec_test_dist(n_count)**2)-1.0/vec_test_dist(n_count))**2   !d_xx Orbital
	  case(6)
	   orbital=kappa**2*((vec_test(2))**2*(1.0/vec_test_dist(n_count)**3+kappa/vec_test_dist(n_count)**2)-1.0/vec_test_dist(n_count))**2   !d_yy Orbital
	  case(7)
	   orbital=kappa**2*(height**2*(1.0/vec_test_dist(n_count)**3+kappa/vec_test_dist(n_count)**2)-1.0/vec_test_dist(n_count))**2        !d_zz Orbital
	  case(8)
	   orbital=(kappa**2/vec_test_dist(n_count)**4)*(vec_test(1))**2*(vec_test(2))**2*(kappa+1.0/vec_test_dist(n_count))**2              !d_xy Orbital
	  case(9)
	   orbital=(kappa**2/vec_test_dist(n_count)**4)*(vec_test(1))**2*(height)**2*(kappa+1.0/vec_test_dist(n_count))**2               !d_xz Orbital
	  case(10)
	   orbital=(kappa**2/vec_test_dist(n_count)**4)*(height)**2*(vec_test(2))**2*(kappa+1.0/vec_test_dist(n_count))**2                !d_yz Orbital
	  case(11)
	   orbital=(kappa/vec_test_dist(n_count))**2*(0.3*(vec_test(1))-0.7*(vec_test(2)))**2           !sonstiges

	  case default
	   call print_line
	   write(6,'(a)') 'Ungultige Tip-Orbital Option:',orb_choice
	   call print_line
	   stop
	  end select
	
! ##########################################################################
! principal partie pour l'extrapolation de la densite



!          cos_phi=(lx*spins(k,1)+ly*spins(k,2)+lz*spins(k,3))/(latt*sqrt(spins(k,1)**2+spins(k,2)**2+spins(k,3)**2))
       cos_phi=spins(k,3)/sqrt(spins(k,1)**2+spins(k,2)**2+spins(k,3)**2)
	   tamr=gamm*cos_phi**2

	   erg(i,j)=erg(i,j)+exp(-2.0*kappa*vec_test_dist(n_count))*(1.0+tamr+p_eff*cos_theta)*orbital
	   fieldx=3.0d0*vec_test(1)*(vec_test(1)*spins(k,1)+vec_test(2)*spins(k,2)+height*spins(k,3))-(r**2+height**2)*spins(k,1)
	   fieldy=3.0d0*vec_test(2)*(vec_test(1)*spins(k,1)+vec_test(2)*spins(k,2)+height*spins(k,3))-(r**2+height**2)*spins(k,2)
	   fieldz=3.0d0*height*(vec_test(1)*spins(k,1)+vec_test(2)*spins(k,2)+height*spins(k,3))-(r**2+height**2)*spins(k,3)
	   field(i,j)=field(i,j)+sqrt(fieldx**2+fieldy**2+fieldz**2)/pi(4.0d0)/vec_test_dist(n_count)**5
           
	  endif
     end do ! end loop n_raster

     if (minval(vec_test_dist).lt.sqrt(cutdown**2+height**2)) then
      sample(i,j)=1.0d0
     endif

	end do

! Ausgabe einer Statusmeldung alle 10 Spalten:
	 if(mod(i,10).eq.0) then
	  write(6,'(2(a,2x,i10))') 'column',i,' of',xmax
	 endif

	end do

        if (all(periodic_log)) then
         minimum=minval(erg)
        else
         minimum=minval(erg,mask=sample.eq.1.0d0)
        endif

        maximum=maxval(erg)

        write(6,'(a,I10)') "Atoms taken into account for the image",int(sum(sample))
        write(6,'(a,2E20.12E3)') "extrema of the contrast in the domain ",maximum,minimum
        write(6,'(a,2E20.12E3)') "extrema of the contrast ", minval(erg)

        do i=1,xmax
         do j=1,ymax
          if (sample(i,j).eq.0.0d0) erg(i,j)=minimum 
         enddo
        enddo
! ###################################################################
! ############# Auswertung, Ausgabe der Intensität ##################
! ###################################################################

! Initialisierung des Bild-Arrays:
	allocate(image(xmax*x_units,ymax*y_units))

! delta-I
         erg=erg-minimum
        
! Erzeuge Array 'image' mit x_units*y_units berechneten Bildausschnitten
	do i=1,x_units*xmax
	 do j=1,y_units*ymax
          x_index=i
          y_index=j
!          do k=1,x_units
!           if (i .gt. (k-1)*xmax) x_index=i-(k-1)*xmax
!          enddo
!          do k=1,y_units
!           if (j .gt. (k-1)*ymax) y_index=j-(k-1)*ymax
!          enddo
!          if (erg(x_index,y_index).gt.0.0d0) then
	   image(i,j)=erg(x_index,y_index)/(2.0d0*kappa*(minimum+0.5d0*(maximum-minimum)))
!      image(i,j)=erg(x_index,y_index)/((minimum+0.5d0*(maximum-minimum)))
!          else
!           image(i,j)=0.0d0
!          endif
	 enddo
	enddo


! Augabe der Ergebnismatrix:
	open(23,file='image.dat')
	do k=1,y_units*ymax
#ifdef CPP_GFORTRAN
         write(23,'(E20.10E3'//repeat(',E20.10E3',x_units*xmax-1)//')') image(:,k)
#else
         write(23,'(E20.10E3,'//repeat(',E20.10E3',x_units*xmax-1)//')') image(:,k)
#endif
    enddo
	close(23)

        open(23,file='sample.dat')
        do k=1,ymax
#ifdef CPP_GFORTRAN
         write(23,'(I2'//repeat(',I2',xmax-1)//')') int(sample(:,k))
#else
         write(23,'(I2,'//repeat(',I2',xmax-1)//')') int(sample(:,k))
#endif
        enddo
        close(23)


    open(23,file='mathematica.dat')
    do i=1,xmax
     do j=1,ymax
      write(23,'(3(E20.10E3))') (radius(k,i,j),k=1,2),erg(i,j)/(2.0d0*kappa*(minimum+0.5d0*(maximum-minimum)))
     enddo
    enddo
    close(23)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! part of the dipolar field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if (all(periodic_log)) then
         minimum=minval(field)
        else
         minimum=minval(field,mask=sample.eq.1.0d0)
       endif

        maximum=maxval(field)

        write(6,'(a,I10)') "Atoms taken into account for the image",int(sum(sample))
        write(6,'(a,2E20.12E3)') "extrema of the contrast of the emmited field in the domain ",maximum,minimum
        write(6,'(a,2E20.12E3)') "extrema of the contrast ", minval(field)

        do i=1,xmax
         do j=1,ymax
          if (sample(i,j).eq.0.0d0) field(i,j)=minimum
         enddo
        enddo
! ###################################################################
! ############# Auswertung, Ausgabe der Intensität ##################
! ###################################################################

         field=field-minimum

! Erzeuge Array 'image' mit x_units*y_units berechneten Bildausschnitten
    do i=1,x_units*xmax
     do j=1,y_units*ymax
          x_index=i
          y_index=j
       image(i,j)=field(x_index,y_index)/(minimum+0.5d0*(maximum-minimum))
     enddo
    enddo


! Augabe der Ergebnismatrix:
    open(23,file='dipolarfield.dat')
    do k=1,y_units*ymax
#ifdef CPP_GFORTRAN
         write(23,'(E20.10E3'//repeat(',E20.10E3',x_units*xmax-1)//')') image(:,k)
#else
         write(23,'(E20.10E3,'//repeat(',E20.10E3',x_units*xmax-1)//')') image(:,k)
#endif
    enddo
    close(23)


    open(23,file='dip-mathematica.dat')
    do i=1,xmax
     do j=1,ymax
      write(23,'(3(E20.10E3))') (radius(k,i,j),k=1,2),field(i,j)/(minimum+0.5d0*(maximum-minimum))
     enddo
    enddo
    close(23)

    deallocate(radius)








	







! ###################################################################
! ######### Erzeugung des Plots der Spinkonfiguration  in 2D ########
! ###################################################################
        si=trim(si1)//trim(',')//trim(si2)
! uberprufung, ob es sich um 2D oder 3D Plot handelt:
	if(spin_control .eq. 0) then
! Erzeuge Bild mit Spinausrichtungen des betrachteten Bereichs zum Vergleich:

! Lege nun fur jede Spinausrichtung eine Datei mit Spinvektoren an:
	do i=1,n_raster

	 if (i .lt. 10) then
	  write(nummer,'(I1)') i
	 elseif (i .lt. 100) then
	  write(nummer,'(I2)') i
	 elseif (i .lt. 1000) then
	  write(nummer,'(I3)') i
	 elseif (i .lt. 10000) then
	  write(nummer,'(I4)') i
	 else
	  write(nummer,'(I5)') i
	 end if
	 datei='spin'//trim(nummer)//'.dat'
	 open(10,file=datei)

! uberprufe jedes Atom im Raster auf i-te Spinkonfiguration:
	  do k=1,n_raster

! Position des Atoms ist Mittelpunkt des Vektors:
	   vec_x1=inbox(k,1)-xstart-0.5*spins(k,1)*vec_len
	   vec_x2=spins(k,1)*vec_len
!
	   vec_y1=inbox(k,2)-ystart-0.5*spins(k,2)*vec_len
	   vec_y2=spins(k,2)*vec_len

! Dubliziere jeweils den Raster-Bereich x_units*y_units mal, um den identischen
! Bereich wie in image.png auszugeben:
	   do j=1,x_units
	    do h=1,y_units
	     write(10,*) vec_x1+(j-1)*(xende-xstart),vec_y1+(h-1)*(yende-ystart),vec_x2,vec_y2
!       write(10,*),vec_x1,vec_y1,vec_x2,vec_y2
	    enddo
	   enddo
	  end do
! Schliesse die i-te Spinkonfigurationsdatei :
	close(10)
	enddo





! Erzeugung der .plot Datei zur Ausgabe der Spinkonfiguration:
	komma=','
	open(19,file='spins.plot')
!	write(19,*),'unset border'
	write(19,*) 'unset xtics'
	write(19,*) 'unset ytics' 
	write(19,*) 'set tmargin 0' 
	write(19,*) 'set lmargin 0' 
	write(19,*) 'set rmargin 0' 
	write(19,*) 'set bmargin 0' 
	write(19,*) 'unset key'
	write(19,*) 'set size ratio '//trim(ratio)
	write(19,*) 'set terminal png size '//trim(si)
	write(19,*) 'set output "spins.png" '
	string='plot '//trim(bound)
	do k=1,n_raster
	
	if (k .lt. 10) then
	 write(nummer,'(I1)') k
	elseif (k .lt. 100) then
	 write(nummer,'(I2)') k
	elseif (k .lt. 1000) then
	 write(nummer,'(I3)') k
	elseif (k .lt. 10000) then
	 write(nummer,'(I4)') k
	else
	 write(nummer,'(I5)') k
	end if

	if (k .eq. n_raster) komma=' '
	string=trim(string)//' "spin'//trim(nummer)//'.dat" w vector lw 1'//komma

	enddo
	write(19,*) string
	close(19)



! Erzeugung der .png Datei der Spinkonfiguration
!	call system ('gnuplot "spins.plot" ')





	endif





! ###################################################################
! ########### Erzeugung des Plots der Spinkonfiguration in 3D #######
! ###################################################################

	
	if(spin_control .eq. 1) then

! Ermittle zunächst maximalen Betrag der z-Komponente der Magnetisierung:

	max_z=maxval(spins(:,3))
	min_z=minval(spins(:,3))
	abs_z=max(abs(min_z),abs(max_z))

	elemente(:)=0
	do k=1,n_raster

! Zuordnung der *.dat Datei, je nach Betrag der z-Komponente
	 do j=0,4
	  if((abs(spins(k,3)) .ge. j*abs_z/5.0).and.(abs(spins(k,3)) .le. (j+1.0)*abs_z/5.0)) then
	   file_index=j
	   if(spins(k,3) .lt. 0.0) file_index=-file_index
	  endif
	 enddo

	 betr=sqrt(spins(k,1)**2+spins(k,2)**2)

	 if(betr .ne. 0.0) then

	if((file_index .eq. -1) .or. (file_index .eq. 1)) then
	datei='spin3.dat'
	elemente(3)=elemente(3)+1
	elseif((file_index .eq. -2) .or. (file_index .eq. -3)) then
	datei='spin2.dat'
	elemente(2)=elemente(2)+1
	elseif((file_index .eq. -4) .or. (file_index .eq. -5)) then
	datei='spin1.dat'
	elemente(1)=elemente(1)+1
	elseif((file_index .eq. 2) .or. (file_index .eq. 3)) then
	datei='spin4.dat'
	elemente(4)=elemente(4)+1
	else
	datei='spin5.dat'
	elemente(5)=elemente(5)+1
	endif

	open(80,file=datei)
! Normierung der Spin-Vektoren in der x-y-Ebene
	spins(k,1)=spins(k,1)/betr
	spins(k,2)=spins(k,2)/betr

	vec_x1=inbox(k,1)-xstart-1.0d0-0.5*spins(k,1)*vec_len
!      vec_x1=inbox(k,1)
	vec_x2=spins(k,1)*vec_len

	vec_y1=inbox(k,2)-ystart-1.0d0-0.5*spins(k,2)*vec_len
!       vec_y1=inbox(k,2)
	vec_y2=spins(k,2)*vec_len
!       write(*,*) 'toto'
!       write(80,*),vec_x1,vec_y1,vec_x2,vec_y2

	do i=1,x_units
	do j=1,y_units
	write(80,*) vec_x1+(i-1.0)*(xende-xstart),vec_y1+(j-1.0)*(yende-ystart),vec_x2,vec_y2
	enddo
	enddo

	close(80)

	else
	if(spins(k,3) .gt. 0.0) then
	open(80,file='plus.dat')
	do i=1,x_units
	do j=1,y_units
	write(80,*) inbox(k,1)-xstart+(i-1.0)*(xende-xstart),inbox(k,2)+(j-1.0)*(yende-ystart)-ystart
	enddo
	enddo
	close(80)
	elemente(6)=elemente(6)+1
	else
	open(80,file='minus.dat')
	do i=1,x_units
	do j=1,y_units
	write(80,*) inbox(k,1)-xstart+(i-1.0)*(xende-xstart),inbox(k,2)+(j-1.0)*(yende-ystart)-ystart
	enddo
	enddo
	close(80)
	elemente(7)=elemente(7)+1
	endif
	endif

	enddo


! Erzeugung der .plot Datei zur Ausgabe der Spinkonfiguration:
	komma=' '
	open(19,file='spins.plot')

	write(19,*) 'set xrange [0:'//trim(xlength)//']'
	write(19,*) 'set yrange [0:'//trim(ylength)//']'
	write(19,*) 'unset xtics'
	write(19,*) 'unset ytics' 
	write(19,*) 'set tmargin 0' 
	write(19,*) 'set lmargin 0' 
	write(19,*) 'set rmargin 0' 
	write(19,*) 'set bmargin 0' 
	write(19,*) 'unset key'
	write(19,*) 'set terminal png size '//trim(si)
	write(19,*) 'set output "spins.png" '
	
	punkte='plot '//trim(bound)

	if(elemente(1) .ne. 0) then
	punkte=trim(punkte)//trim(komma)//' "./spin1.dat" w vector lw 1 lc rgb "#ff0000"'
	komma=','
	endif

	if(elemente(2) .ne. 0) then
	punkte=trim(punkte)//trim(komma)//' "./spin2.dat" w vector lw 1 lc rgb "#ff9900"'
	komma=','
	endif

	if(elemente(3) .ne. 0) then
	punkte=trim(punkte)//trim(komma)//' "./spin3.dat" w vector lw 1 lc rgb "#0000ff"'
	komma=','
	endif

	if(elemente(4) .ne. 0) then
	punkte=trim(punkte)//trim(komma)//' "./spin4.dat" w vector lw 1 lc rgb "#00ffff"'
	komma=','
	endif

	if(elemente(5) .ne. 0) then
	punkte=trim(punkte)//trim(komma)//' "./spin5.dat" w vector lw 1 lc rgb "#00ff00"'
	komma=','
	endif

	if(elemente(6) .ne. 0) then
	punkte=trim(punkte)//trim(komma)//' "./plus.dat" lc rgb "#00ff00" pt 19'
	komma=','
	endif

	if(elemente(7) .ne. 0) then
	punkte=trim(punkte)//trim(komma)//' "./minus.dat" lc rgb "#ff0000" pt 2'
	komma=','
	endif

	write(19,*) trim(punkte)
	close(19)
	call system ('gnuplot "./spins.plot" ')

	endif









	




! Scan-Linien:

!	lauf=int(sqrt((xe-xs)**2+(ye-ys)**2))*2  !Anzahl Messpunkte
!	speicher=0.0
!
!	open(97,file='scan.dat')
!
!	k=0
!	do i=1,10*lauf
!
!	x_scan=xs+dble(i-1)/10.0d0/dble(lauf)*abs(xe-xs)
!	y_scan=ys+dble(i-1)/10.0d0/dble(lauf)*abs(ye-ys)
!        x_ind=minloc(abs(inbox(:,1)-x_scan),1)
!        y_ind=minloc(abs(inbox(:,2)-y_scan),1)
!
!        speicher=speicher+image(x_ind,y_ind)
!
!	if(mod(i,10) .eq. 0) then
!
!	write(97,*) sqrt((x_scan-xs)**2+(y_scan-ys)**2)/sqrt((xe-xs)**2+(ye-ys)**2), &
!     & (speicher/10.0)/(2.0*kappa* &
!     & (minimum+0.5*(maximum-minimum)))
!
!	speicher=0.0
!        k=k+1
!	endif
!
!	enddo


	
! Ausgabe der Corrugation Amplitude, I_0 in der Mitte von maximum und minimum:
	call print_line()
	write(6,*) 'corrugation amplitude:',(maximum-minimum)/(2.0*kappa*(minimum+0.5*(maximum-minimum)))


! Ausgabe der Corrugation Amplitude, I_0 bei minimum, als Vergleich für Genauigkeit ...
! ... der linearen Näherung
	call print_line()
	if (dabs(minimum).gt.1.0d-8) write(6,*) 'corrugation amplitude (vergleich):',(maximum-minimum)/(2.0*kappa*minimum)
	call print_line()



! Aufruf der diskreten FFT:

	allocate(fft(xmax,ymax))
	allocate(dat(2*xmax*ymax))
	allocate(speicher1(xmax))

	dat(:)=0.0d0


	do j=1,ymax
	do i=1,xmax
	if((i .le. xmax) .and. (j .le. ymax)) then
	dat(2*i-1+(j-1)*2*xmax)=erg(i,j)
	endif

	enddo
	enddo

	nn(1)=xmax
	nn(2)=ymax
	call fourn(dat,nn,2,1)


	open(77,file='fft.dat')
	do i=1,ymax
	do k=1,xmax
	fft(k,i)=sqrt(dat((i-1)*2*xmax+2*k-1)**2+1.0*dat((i-1)*2*xmax+2*k)**2)
	enddo
	do j=1,xmax/2
	he=fft(j,i)
	fft(j,i)=fft(xmax/2+j,i)
	fft(xmax/2+j,i)=he
	enddo
	enddo

	do j=1,ymax/2
	speicher1(:)=fft(:,j)
	fft(:,j)=fft(:,ymax/2+j)
	fft(:,ymax/2+j)=speicher1(:)
	enddo

	do j=1,ymax
        write(blaaa,'(a,I5,12a)') '(',xmax,'(2X,f20.15))'
	write(77,trim(blaaa)) fft(:,j)
	enddo
	close(77)





	open(77,file='make_image.plot')
	write(77,*) 'set autoscale fix'
	write(77,*) 'unset colorbox'
	write(77,*) 'unset key'
	write(77,*) 'set terminal png size '//trim(si)
	write(77,*) 'set output "image_tics.png"'
	write(77,*) 'plot "image.dat" matrix w image'
!	write(77,*),'unset border' 
	write(77,*) 'set size ratio '//trim(ratio)
	write(77,*) 'unset xtics'
	write(77,*) 'unset ytics '
	write(77,*) 'set tmargin 0'
	write(77,*) 'set lmargin 0'
	write(77,*) 'set rmargin 0'
	write(77,*) 'set bmargin 0'
	write(77,*) 'set output "image.png"'

	write(77,*) 'replot'
	




	write(77,*) 'set output "fft.png"'
	write(77,*) 'plot "fft.dat" matrix w image'

! Neue sw Ausgabe:

	write(77,*) 'set output "bw.png"'
	write(77,*) 'set palette grey'
	WRITE(77,*) 'plot "image.dat" matrix w image'

	write(77,*) 'set output "bw_scale.png" '

	write(77,*) 'set colorbox'

	write(77,*) 'set colorbox'

	write(77,*) 'set rmargin 6'

	
	write(77,*) 'set autoscale fix'

	write(77,*) 'replot'




	close(77)



	call system ('gnuplot ./make_image.plot ')
! Aufraumen des Ordners ./data/ nach Ausführung des Programms:
!	call system('rm -rf ./data/*.dat')
!	call system('rm -rf ./data/*.plot')
!	call system('rm -rf ./*.plot')
	end subroutine




	subroutine fourn(dat,nn,ndim,isign)
! Subroutine fur diskrete FFT
! entnommen aus:
!
! Numerical Recipes in Fortran
! Cambridge University Press, 1992
! W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery
	implicit none
	integer :: isign,ndim,nn(ndim)
	real(kind=8) :: dat(*)
	integer :: i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2, &
     & ip3,k1,k2,n,nprev,nrem,ntot
	real(kind=8) :: tempi,tempr
	real(kind=8) :: theta,wi,wpi,wpr,wr,wtemp
	ntot=1
	do idim=1,ndim
	ntot=ntot*nn(idim)
	enddo
	nprev=1
	do idim=1,ndim
	n=nn(idim)
	nrem=ntot/(n*nprev)
	ip1=2*nprev
	ip2=ip1*n
	ip3=ip2*nrem
	i2rev=1
	do i2=1,ip2,ip1
	if(i2 .lt. i2rev) then
	do i1=i2,i2+ip1-2,2
	do i3=i1,ip3,ip2
	i3rev=i2rev+i3-i2
	tempr=dat(i3)
	tempi=dat(i3+1)
	dat(i3)=dat(i3rev)
	dat(i3+1)=dat(i3rev+1)
	dat(i3rev)=tempr
	dat(i3rev+1)=tempi
	enddo
	enddo
	endif
	ibit=ip2/2
1	if ((ibit .ge. ip1) .and. (i2rev .gt. ibit)) then
	i2rev=i2rev-ibit
	ibit=ibit/2
	goto 1
	endif
	i2rev=i2rev+ibit
	enddo
	ifp1=ip1
2	if(ifp1 .lt. ip2) then
	ifp2=2*ifp1
	theta=isign*6.28318530717959d0/(ifp2/ip1)
	wpr=-2.d0*sin(0.5d0*theta)**2
	wpi=sin(theta)
	wr=1.d0
	wi=0.d0
	do i3=1,ifp1,ip1
	do i1=i3,i3+ip1-2,2
	do i2=i1,ip3,ifp2
	k1=i2
	k2=k1+ifp1
	tempr=sngl(wr)*dat(k2)-sngl(wi)*dat(k2+1)
	tempi=sngl(wr)*dat(k2+1)+sngl(wi)*dat(k2)
	dat(k2)=dat(k1)-tempr
	dat(k2+1)=dat(k1+1)-tempi
	dat(k1)=dat(k1)+tempr
	dat(k1+1)=dat(k1+1)+tempi     
	enddo
	enddo
	wtemp=wr
	wr=wr*wpr-wi*wpi+wr
	wi=wi*wpr+wtemp*wpi+wi
	enddo
	ifp1=ifp2
	goto 2
	endif
	nprev=n*nprev
	enddo
	return
	end

	real(kind=8) function cosinus(a1,a2,a3,b1,b2,b3)
	implicit none
	real(kind=8) :: a1,a2,a3,b1,b2,b3
	cosinus=(a1*b1+a2*b2+a3*b3)/(sqrt(a1**2+a2**2+a3**2)*sqrt(b1**2+b2**2+b3**2))
	return
	end function



