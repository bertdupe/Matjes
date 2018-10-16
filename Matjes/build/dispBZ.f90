      integer function fullBZ(d,N_ip,N_op,phase,my_lattice,motif)
      use m_vector, only : cross,norm
      use m_sym_utils, only : mesh_bz,order_zaxis
      use m_constants, only : pi
      use m_parameters, only :c_Ji,c_jB,DM,c_dm,J_ij,J_B,i_DM,i_biq,J_il,i_sliptun
      use m_lattice, only : indexNN
      use m_derived_types
      implicit none
! internal variable
      type(lattice), intent(in) :: my_lattice
      type(cell), intent(in) :: motif
      integer, intent(in) :: N_ip,N_op,phase
      real(kind=8), intent(in) :: d(max(N_ip,N_op),phase)
! mesh for the BZ
      real(kind=8), allocatable :: mesh(:,:)
      integer :: nkpt,fin,N_sampling
      real(kind=8) :: E_q0
! dummy
      integer :: i1,k,i,l,n_sym,i_n_sym,worked,dim_lat(3)
      real(kind=8) :: kv0(3,3),E_int,Rq(3),Iq(3),qvec(3),rotation(3,3)
      real(kind=8) :: rotangle,ss,net(3,3)
      real(kind=8), allocatable :: location(:,:)
      logical :: exists
      character(len=5) :: str

      net=my_lattice%areal
      dim_lat=my_lattice%dim_lat

      fullBZ=0
      Rq=(/0.0d0,0.0d0,1.0d0/)
      E_q0=c_Ji*sum(J_ij(1:N_ip,1)*dble(indexNN(1:N_ip,1)))
      if (i_sliptun) then
       E_q0=E_q0+c_Ji*sum(J_ij(1:N_ip,2)*dble(indexNN(1:N_ip,1)))
       else
       E_q0=E_q0+c_Ji*sum(J_ij(1:N_ip,1)*dble(indexNN(1:N_ip,1)))
      endif
      if (phase.ne.1) E_q0=E_q0+c_Ji*sum(J_il(1:N_op)*dble(indexNN(1:N_op,2)))*dble(phase)
      if (i_biq) E_q0=E_q0+c_JB*J_B*dble(indexNN(1,1))

!! reciprocal lattice vector
      kv0(1,:)=pi(2.0d0)*cross(net(2,:),net(3,:))/dot_product(net(1,:),cross(net(2,:),net(3,:)))
      kv0(2,:)=pi(2.0d0)*cross(net(3,:),net(1,:))/dot_product(net(1,:),cross(net(2,:),net(3,:)))
      kv0(3,:)=pi(2.0d0)*cross(net(1,:),net(2,:))/dot_product(net(1,:),cross(net(2,:),net(3,:)))

      n_sym=abs(order_zaxis(kv0))

      OPEN(2,FILE='inp',action='read',status='old',form='formatted')
      do
       read(2,'(a)',iostat=fin) str
       if (fin /= 0) exit
       if ( str(1:5) == 'dispe') then
        backspace(2)
        read(2,*) str, str, N_sampling
       endif
      enddo
      close(2)

      inquire (file='kpoints',exist=exists)
      if (.not.exists) worked = mesh_bz(n_sym,kv0,N_sampling)

      OPEN(2,FILE='kpoints',action='read',status='old',form='formatted')

      nkpt=0
      do
       read(2,'(a)',iostat=fin) str
       if (fin /= 0) exit
       nkpt=nkpt+1
      end do
      rewind(2)
      allocate(mesh(3,nkpt))
      do l=1,nkpt
       read (2,*) (mesh(i,l),i=1,3)
      end do
      close(2)

      OPEN(1,FILE='dispersion.dat',action='write',status='unknown',form='formatted')

      write(1,'(4f14.8)') 0.0d0,0.0d0,0.0d0

      i_n_sym=0
      do while (i_n_sym.lt.n_sym)
!! E=S^2 sum_j J_j exp(-iqR_j)
       rotangle = dble(i_n_sym)*pi(2.0d0)/dble(n_sym)

       rotation=reshape((/dcos(rotangle),-dsin(rotangle),0.0d0, &
        dsin(rotangle),dcos(rotangle),0.0d0,&
        0.0d0,0.0d0,1.0d0/),(/3,3/))

       do i1=2,nkpt

        qvec=matmul(transpose(rotation),mesh(:,i1))

        Iq=qvec
        ss=norm(Iq)
        if (dabs(ss).lt.1.0d-7) ss=1.0d0
        Iq=Iq/ss

        E_int=0.0d0
!---------------------
! normal case
        do k=1,N_ip

         allocate(location(3,indexNN(k,1)))
         location=pos_nei(d(k,1),k,net,indexNN(k,1))

         l=1
          do while(l.lt.(indexNN(k,1)+1))

           E_int=E_int+c_Ji*J_ij(k,1)*dcos(dot_product(qvec,location(:,l)))
           if (i_sliptun) then
            E_int=E_int+c_Ji*J_ij(k,2)*dcos(dot_product(qvec,location(:,l)))
            else
            E_int=E_int+c_Ji*J_ij(k,1)*dcos(dot_product(qvec,location(:,l)))
           endif
           if ((k.eq.1).and.(i_biq)) E_int=E_int+c_JB*J_B*dcos(2.0d0*dot_product(qvec,location(:,l)))
           if ((k.le.count(abs(DM(:,1))>1.0d-7)).and.(i_DM)) E_int=E_int+C_DM*DM(k,1)*0.5d0* &
             dsin(dot_product(qvec,location(:,l)))*dot_product(qvec,location(:,l)) &
             /norm(qvec)/norm(location(:,l))

           l=l+1
          enddo
          deallocate(location)
         enddo
!---------------------
!---------------------
! case of the superlattice
        if (phase.eq.2) then
         do k=1,N_op

          allocate(location(3,indexNN(k,2)))
          location=pos_nei_SL(d(k,2),k,net,indexNN(k,2),motif)


          l=1
           do while(l.lt.(indexNN(k,2)+1))

            E_int=E_int+c_Ji*J_il(k)*dcos(dot_product(qvec,location(:,l)))*dble(phase)

            l=l+1
           enddo
          deallocate(location)
         enddo
        endif
!---------------------
         write(1,'(4f14.8)') (qvec(i),i=1,2),(E_int-E_q0)*1000.0d0

       enddo

      i_n_sym=i_n_sym+1
      enddo
      close(1)
      call system ('rm "kpoints" ')
      stop

      contains

      function pos_nei(d,i_nei,r,Nb_nei)
      use m_vector, only : norm
      implicit none
      real(kind=8) :: pos_nei(3,Nb_nei)
      integer, intent(in) :: i_nei,Nb_nei
      real(kind=8), intent(in) :: d,r(3,3)
!!! dummy
      real(kind=8) :: vec
      integer :: i,j,i_Nb_nei

      i_Nb_nei=1

      i=-i_nei
      do while (i.lt.(i_nei+1))
      j=-i_nei
       do while (j.lt.(i_nei+1))

        vec=norm(dble(i)*r(1,:)+dble(j)*r(2,:))

        if ((dabs(vec-d).lt.1.0d-8).and.(i_Nb_nei.lt.Nb_nei+1)) then
         pos_nei(:,i_Nb_nei)=dble(i)*r(1,:)+dble(j)*r(2,:)

         i_Nb_nei=1+i_Nb_nei
        endif
         
       j=j+1
       enddo
      i=i+1
      enddo

      end function pos_nei

      function pos_nei_SL(d,i_nei,r,Nb_nei,motif)
      use m_vector, only : norm
      use m_derived_types, only : cell
      implicit none
      real(kind=8) :: pos_nei_SL(3,Nb_nei)
      real(kind=8), intent(in) :: d,r(3,3)
      integer, intent(in) :: i_nei,Nb_nei
      type(cell), intent(in) :: motif
!!! dummy
      real(kind=8) :: vec
      integer :: i,j,i_Nb_nei,i_m

      i_Nb_nei=1

      do i_m=1,size(motif%i_mom)
      if ((.not.motif%i_mom(i_m)).or.(dabs(sum(motif%pos(i_m,:))).lt.1.0d-8)) cycle
       i=-i_nei
       do while (i.lt.(i_nei+1))
        j=-i_nei
        do while (j.lt.(i_nei+1))

         vec=norm((dble(i)+motif%pos(i_m,1))*r(1,:)+(dble(j)+motif%pos(i_m,2))*r(2,:)+ &
          motif%pos(i_m,3)*r(3,:))

         if ((dabs(vec-d).lt.1.0d-8).and.(i_Nb_nei.lt.Nb_nei+1)) then
          pos_nei_SL(:,i_Nb_nei)=(dble(i)+motif%pos(i_m,1))*r(1,:)+(dble(j)+motif%pos(i_m,2))*r(2,:)+ &
           motif%pos(i_m,3)*r(3,:)

          i_Nb_nei=1+i_Nb_nei
         endif

        j=j+1
        enddo
       i=i+1
       enddo
      enddo

      end function pos_nei_SL

      end function fullBZ
