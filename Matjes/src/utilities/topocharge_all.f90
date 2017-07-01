      module m_topocharge_all
      interface topo
       module procedure initialize_topo
      end interface topo

      contains

      subroutine initialize_topo(spin,shape_spin,masque,shape_masque,qeulerp,qeulerm)
      use m_rw_lattice, only : n_system
      use m_topocharge_local, only : local_topo
#ifdef CPP_MPI
      use m_make_box, only : Xstart,Xstop,Ystart,Ystop,Zstart,Zstop
#endif
      implicit none
      integer, intent(in) :: shape_spin(:),shape_masque(:),masque(:,:,:,:)
      real(kind=8), intent(in) :: spin(:,:,:,:,:)
      real(kind=8), intent(out) :: qeulerp,qeulerm
! dumy
      integer :: i_x,i_y,i_z,i_m
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

         select case(n_system)
          case(2)
           call local_topo(i_x,i_y,qm,qp,spin,shape_spin,masque,shape_masque)
          case(22)
           call local_topo(i_x,i_y,i_m,qm,qp,spin,shape_spin,masque,shape_masque)
          case(32)
           call local_topo(i_x,i_y,i_z,i_m,qm,qp,spin,shape_spin,masque,shape_masque)
          case default
           call local_topo(i_x,i_y,qm,qp,spin,shape_spin,masque,shape_masque)
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
