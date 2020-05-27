      module m_topocharge_all
      use m_derived_types
      interface topo
       module procedure initialize_topo
      end interface topo

      contains

      subroutine initialize_topo(spin,shape_spin,masque,qeulerp,qeulerm,my_lattice)
      use m_topocharge_local, only : local_topo
#ifdef CPP_MPI
      use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
#endif
      implicit none
      integer, intent(in) :: shape_spin(:),masque(:,:,:,:)
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      type(lattice), intent(in) :: my_lattice
      real(kind=8), intent(out) :: qeulerp,qeulerm
! dumy
      integer :: i_x,i_y,i_z,i_m,n_system
      real(kind=8) :: qp,qm

#ifndef CPP_MPI
      integer :: Xstart,Xstop,Ystart,Ystop,Zstart,Zstop

      Xstart=1
      Xstop=shape_spin(2)
      Ystart=1
      Ystop=shape_spin(3)
      Zstart=1
      Zstop=shape_spin(4)
#endif

      n_system=my_lattice%n_system
      qp=0.0d0
      qm=0.0d0
      qeulerp=0.0d0
      qeulerm=0.0d0

! calculate the topological charge and the vortivity for the different types of system.

      do i_m=1,shape_spin(5)
#ifdef CPP_OPENMP
!$OMP parallel DO REDUCTION(+:qeulerp,qeulerm) private(i_x,i_y,i_z,qm,qp) default(shared)
#endif

       do i_z=Zstart,Zstop
        do i_y=Ystart,Ystop
         do i_x=Xstart,Xstop

         if (masque(1,i_x,i_y,i_z).eq.0) cycle

         select case(n_system)
          case(2)
           call local_topo(i_x,i_y,qm,qp,spin,shape_spin,my_lattice)
          case(3)
           call local_topo(i_x,i_y,i_z,qm,qp,spin,shape_spin,my_lattice)
          case default
           write(6,'(a)') 'case not taken into account into topocharge_all.f90'
           stop
          end select
         qeulerp=qeulerp+qp
         qeulerm=qeulerm+qm

         enddo
        enddo
       enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif
      enddo

      end subroutine initialize_topo

      end module m_topocharge_all
