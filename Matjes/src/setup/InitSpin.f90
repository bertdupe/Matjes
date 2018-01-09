! ===============================================================
      SUBROUTINE InitSpin(r,motif,world,state)
      USE m_lattice, only : spin
      use m_vector, only : cross,norm
      use m_rw_lattice, only : dim_lat
      use m_constants
      use m_derived_types
      use mtprng
#ifdef CPP_MPI
      use m_mpi
      use m_mpi_prop, only : irank,MPI_COMM
#endif
      Implicit none
! variables that come in
      real (kind=8), intent(in) :: r(3,3)
      integer, intent(in) :: world(:)
      type (cell), intent(in) :: motif
      type(mtprng_state), intent(inout) :: state
!     Slope Indexes for three dim spins
      INTEGER :: i_spin,i_coord,i_lat,j_lat,fin,i,j,k,l,ierr
      LOGICAL :: i_exi,tag,heavy,reorder,i_alat,domainwall
!     Absolute value of a spin
      real (kind=8) :: SpinAbs,ss,dumy(3),choice,mu_S,alat(3),coefftheta,coeffphi, alpha
      real (kind=8), dimension(3) :: qvec,Rq,Iq,rpol,kx,ky,kz
      character(len=10) :: str,dummy
      integer, parameter  :: io=9
      integer :: dw_position,  o

#ifdef CPP_MPI
      include 'mpif.h'
      if (irank.eq.0) then 
#endif
      mu_S=motif%mom(1)
      Spin(7,:,:,:,:)=mu_S

      tag=.False.
      i_exi=.False.
      i_alat=.False.
      domainwall=.False.
!     Check if spin lattice input is present
      inquire(file='SpinSTMi.dat',exist=i_exi)
!     Reset of Spinlattice
      Spin(1:6,:,:,:,:)=0.0d0
      if (i_exi) then
       write(6,'(a)') 'Read spin lattice from SpinSTMi.dat'
       open(io,file='SpinSTMi.dat',form='formatted',status='old')
!cc check for old format
       open (1,file='inp',form='formatted',status='old',action='read')
       rewind(1)
       do
        read (1,'(a)',iostat=fin) str
        if (fin /= 0) exit
        str= trim(adjustl(str))
        if (len_trim(str)==0) cycle

        if (str(1:1) == '#' ) cycle 
         if ( str(1:9) == 'oldformat') then
          backspace(1)
           read(1,*) dummy, tag
         endif
       enddo
       close(1)
!cccccccccccccccccc end reading for old format
       inquire(file='alat.in',exist=i_alat)
       if (i_alat) then
        open(142,file='alat.in',form='formatted',status='old')
        read(142,*) (alat(l),l=1,3)
        close(142)
       endif

       if (.not.tag) then

         do k=1,dim_lat(3)
          do j=1,dim_lat(2)
           do i=1,dim_lat(1)
         read(io,*) ((Spin(j_lat,i,j,k,l),j_lat=1,6),l=1,count(motif%i_m))
           enddo
          enddo
         enddo

        if (i_alat) then
         do l=1,count(motif%i_m)
          do k=1,dim_lat(3)
           do j=1,dim_lat(2)
            do i=1,dim_lat(1)
           Spin(1:3,i,j,k,l)=Spin(1:3,i,j,k,l)*alat
            enddo
           enddo
          enddo
         enddo
        endif
!truc bizare
!        do i=1,dim_lat(1)
!        do j=1,dim_lat(2)
!         Spin(4:6,i,j,1,1)=Spin(4:6,i,72,1,1)
!         Spin(4:6,i,j,1,2)=Spin(4:6,i,72,1,2)
!
!        enddo
!        enddo
! Check for reordering the lattice
      reorder=.False.
      inquire(file='reorder.lattice',exist=reorder)
      if (reorder) call order_lattice(spin,motif,r,dim_lat)
       else
        do j=1,dim_lat(2)
         do i=1,dim_lat(1)
          do k=1,dim_lat(3)
         read(io,*) dummy,dummy,(Spin(j_lat,i,j,k,1), j_lat=4,6)
         spin(1:3,i,j,k,1)=dble(i-1)*r(1,:)+dble(j-1)*r(2,:)+dble(k-1)*r(3,:)
          enddo
         enddo
        enddo
       endif
      else
       do k=1,dim_lat(3)
        do j=1,dim_lat(2)
         do i=1,dim_lat(1)
          do l=1,size(motif%i_m)
           if (.not.(motif%i_m(l))) cycle

! fix the spin direction as random
          Do i_spin=4,6
#ifdef CPP_MRG
      Choice=mtprng_rand_real1(state)
#else
      CALL RANDOM_NUMBER(Choice)
#endif
            Spin(i_spin,i,j,k,l)=2.0d0*(choice-0.5d0)

          End do

!fix the coordinates. This part takes care of the motif of the structure
         Spin(1:3,i,j,k,l)=r(1,:)*(dble(i-1)+motif%pos(l,1))+r(2,:)*(dble(j-1)+motif%pos(l,2))+ &
          r(3,:)*(dble(k-1)+motif%pos(l,3))

       enddo
       enddo
       enddo
       enddo
      endif

!     Check if spin lattice input is present
      i_exi=.False.
      inquire(file='spiral.start',exist=i_exi)
!     Reset of Spinlattice
      if (i_exi) then
      Spin(4:6,:,:,:,:)=0.0d0
      write(6,'(/,a,/)') 'Starting configuration as spin spiral'
      open(io,file='spiral.start',form='formatted',status='old')
      rewind(io)
      do
      read (io,'(a)',iostat=fin) str
        if (fin /= 0) exit
        str= trim(adjustl(str))
        if (len_trim(str)==0) cycle

        if (str(1:1) == '#' ) cycle
!cccc We start to read the input
        if ( str(1:4) == 'qvec') then
           backspace(io)
           read(io,*) dummy,(qvec(i),i=1,3)
           kx(:) = qvec(1)*pi(2.0d0)*cross(r(2,:),r(3,:))/dot_product(r(1,:),cross(r(2,:),r(3,:)))
           ky(:) = qvec(2)*pi(2.0d0)*cross(r(3,:),r(1,:))/dot_product(r(1,:),cross(r(2,:),r(3,:)))
           kz(:) = qvec(3)*pi(2.0d0)*cross(r(1,:),r(2,:))/dot_product(r(1,:),cross(r(2,:),r(3,:)))
           qvec=kx+ky+kz
         endif
        if ( str(1:2) == 'Rq') then
           backspace(io)
           read(io,*) dummy,(Rq(i),i=1,3)
           dumy=Rq(1)*r(1,:)+Rq(2)*r(2,:)+Rq(3)*r(3,:)
           ss=norm(dumy)
           Rq=dumy/ss
         endif
        if ( str(1:2) == 'Iq') then
           backspace(io)
           read(io,*) dummy,(Iq(i),i=1,3)
           dumy=Iq(1)*r(1,:)+Iq(2)*r(2,:)+Iq(3)*r(3,:)
           ss=norm(dumy)
           Iq=dumy/ss
         endif
        if ( str(1:9) == 'heavyside') then
           backspace(io)
           read(io,*) dummy,heavy
         endif
        if ( str(1:10) == 'domainwall') then
           backspace(io)
           read(io,*) dummy,domainwall
         endif
      enddo
      close(io)

! To initialize the spin, the Phd of Ph Kurz uses the convention Si=cos()Rq-sin()Iq. If Iq is pointing along
! the q vector, it then rotates in the sense of -q. I take the convention Si=cos()Rq+sin()Iq

      if ( heavy ) then
       do k=1,dim_lat(3)
        do j=1,dim_lat(2)
         do i=1,dim_lat(1)
          do l=1,size(motif%i_m)
          if (.not.(motif%i_m(l))) cycle

        Spin(4:5,i,j,k,l)=0.0d0
        Spin(6,i,j,k,l)=1-2*(2*i/dim_lat(1))
          enddo
         enddo
        enddo
       enddo

      else
      do l=1,size(motif%i_m)
      if (.not.(motif%i_m(l))) cycle
       do k=1,dim_lat(3)
        do j=1,dim_lat(2)
         do i=1,dim_lat(1)
!normal spin spiral
         Spin(4:6,i,j,k,l)=(dcos(dot_product(qvec,Spin(1:3,i,j,k,l)))*Rq+ &
         dsin(dot_product(qvec,Spin(1:3,i,j,k,l)))*Iq)
! inomegenous spin spiral in mulitlayer
!        if (l.eq.1) Spin(4:6,i,j,k,l)=(/0.0d0,0.0d0,1.0d0/)
! 3x3 cubic
! define by theta and phi
!        coefftheta=dot_product(qvec,Spin(1:3,i,j,k,l))
!        coeffphi=sqrt(Spin(1,i,j,k,l)**2+Spin(2,i,j,k,l)**2)
!        if (abs(coeffphi).gt.1.0d-5) coeffphi=Spin(1,i,j,k,l)/sqrt(Spin(1,i,j,k,l)**2+Spin(2,i,j,k,l)**2)
!        Spin(4:6,i,j,k,l)=(/coeffphi*sin(coefftheta), &
!         sqrt(1-coeffphi**2)*sin(coefftheta), &
!         cos(coefftheta) /)
          enddo
         enddo
        enddo
       enddo

      endif
      endif

      if ( domainwall ) then
       do i=1,dim_lat(1)
        if ( (2*i/dim_lat(1)) == 1) then
         dw_position=i
         exit
        endif
       enddo

       do i=1,dim_lat(1)
        do o=-4,4
        alpha=real(5+o)/10.0d0*acos(-1.0d0)
        write(*,*) o, alpha
         if ( i+o == dw_position ) then
         do j=1,dim_lat(2)
          do k=1,dim_lat(3)
           do l=1,size(motif%i_m)
           if (.not.(motif%i_m(l))) cycle

        Spin(4,i,j,k,l)=sin(alpha)
        Spin(6,i,j,k,l)=-1.0d0*cos(alpha)
           enddo
          enddo
         enddo
         endif
        enddo
       enddo

      endif


!alpha, o


! make sure that the spins have a norm 1

      do l=1,count(motif%i_m)
       do k=1,dim_lat(3)
        do j=1,dim_lat(2)
         do i=1,dim_lat(1)

         ss=sqrt(Spin(4,i,j,k,l)**2+Spin(5,i,j,k,l)**2+Spin(6,i,j,k,l)**2)

         Spin(4:6,i,j,k,l)=Spin(4:6,i,j,k,l)/ss

         enddo
        enddo
       enddo
      enddo

#ifdef CPP_MPI
      endif

      spin=bcast(Spin,7,dim_lat(1),dim_lat(2),dim_lat(3),count(motif%i_m),MPI_REAL8,0,MPI_COMM)

#endif

#ifdef CPP_DEBUG
      do l=1,count(motif%i_m)
       do k=1,dim_lat(3)
        do j=1,dim_lat(2)
         do i=1,dim_lat(1)
          write(*,*) spin(6,i,j,k,l)
         enddo
        enddo
       enddo
      enddo
      stop
#endif

      END SUBROUTINE InitSpin
! ===============================================================


