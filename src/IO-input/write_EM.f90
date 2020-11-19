       module m_write_EM
       interface write_EM
        module procedure write_EM_paratemp
        module procedure write_EM_serial
!        module procedure write_EM_init
       end interface write_EM
       contains

!       subroutine write_EM_init(fname)
!       implicit none
!       integer, intent(in) :: N_cell,n_thousand,n_ghost
!       real(kind=8), intent(in) :: vortex_mpi(:,:),qeulerp_mpi(:),qeulerm_mpi(:),kt_mpi(:), &
!      & E_mpi(:),E_sq_mpi(:),M_sq_mpi(:),M(:,:)
!! internal variable
!       integer :: i
!       real(kind=8) :: cell,C,chi,E_av,steps,qp,qm,vortex(3),E_err,M_err,M_local(3)
!
!       write(6,'(a)') "----------------"
!       write(6,'(a)') "write in EM.dat before changing the temperature set (if applicable)"
!       write(6,'(a)') "----------------"
!
!       cell=dble(N_cell)
!       steps=dble(n_thousand)
!
!       do i=1,isize/n_ghost
!
!        M_local=M(:,i)/steps/cell
!        C=(E_sq_mpi(i)/steps-(E_mpi(i)/steps)**2)/kT_mpi(i)**2/cell
!        chi=(M_sq_mpi(i)/steps-(norm(M(:,i))/steps)**2)/kT_mpi(i)/cell
!        E_av=E_mpi(i)/cell/steps
!        E_err=sqrt(abs(E_sq_mpi(i)-E_mpi(i)**2/steps)/(steps-1.0d0))/dble(N_cell)
!        M_err=sqrt(abs(M_sq_mpi(i)-norm(M(:,i))**2/steps)/(steps-1.0d0))/dble(N_cell)
!        qp=qeulerp_mpi(i)/pi/4.0d0/steps
!        qm=qeulerm_mpi(i)/pi/4.0d0/steps
!        vortex=vortex_mpi(:,i)/3.0d0/sqrt(3.0d0)/steps
!
!        Write(7,'(15(E20.10E3,2x),E20.10E3)') kt_mpi(i)/k_B,E_av,E_err,C, &
!     & norm(M_local(:)), M_local(:),M_err,chi,vortex,qp+qm,qp,qm
!       enddo
!
!        Write(7,'(a)') ' '
!
!       end subroutine write_EM_init

       subroutine write_EM_paratemp(N_cell,n_thousand,vortex_mpi,qeulerp_mpi,qeulerm_mpi,kt_mpi,E_mpi,E_sq_mpi,M_sq_mpi,M,n_ghost,isize)
       use m_vector, only : norm
       use m_constants, only : k_B,pi
       implicit none
       integer, intent(in) :: N_cell,n_thousand,n_ghost,isize
       real(kind=8), intent(in) :: vortex_mpi(:,:),qeulerp_mpi(:),qeulerm_mpi(:),kt_mpi(:), &
      & E_mpi(:),E_sq_mpi(:),M_sq_mpi(:),M(:,:)
! internal variable
       integer :: i
       real(kind=8) :: cell,C,chi,E_av,steps,qp,qm,vortex(3),E_err,M_err,M_local(3)

       write(6,'(a)') "----------------"
       write(6,'(a)') "write in EM.dat before changing the temperature set (if applicable)"
       write(6,'(a)') "----------------"

       cell=dble(N_cell)
       steps=dble(n_thousand)

       do i=1,isize/n_ghost

        M_local=M(:,i)/steps/cell
        C=(E_sq_mpi(i)/steps-(E_mpi(i)/steps)**2)/kT_mpi(i)**2/cell
        chi=(M_sq_mpi(i)/steps-(norm(M(:,i))/steps)**2)/kT_mpi(i)/cell
        E_av=E_mpi(i)/cell/steps
        E_err=sqrt(abs(E_sq_mpi(i)-E_mpi(i)**2/steps)/(steps-1.0d0))/dble(N_cell)
        M_err=sqrt(abs(M_sq_mpi(i)-norm(M(:,i))**2/steps)/(steps-1.0d0))/dble(N_cell)
        qp=qeulerp_mpi(i)/pi/4.0d0/steps
        qm=qeulerm_mpi(i)/pi/4.0d0/steps
        vortex=vortex_mpi(:,i)/3.0d0/sqrt(3.0d0)/steps

        Write(7,'(15(E20.10E3,2x),E20.10E3)') kt_mpi(i)/k_B,E_av,E_err,C, &
     & norm(M_local(:)), M_local(:),M_err,chi,vortex,qp+qm,qp,qm
       enddo

        Write(7,'(a)') ' '

       end subroutine write_EM_paratemp

       subroutine write_EM_serial(N_cell,n_thousand,vortex_mpi,qeulerp_mpi,qeulerm_mpi,kt_mpi,E_mpi,E_sq_mpi,M_sq_mpi,M,n_ghost)
       use m_vector, only : norm
       use m_constants, only : k_B,pi
       use m_mpi_prop, only : isize
       implicit none
       integer, intent(in) :: N_cell,n_thousand,n_ghost
       real(kind=8), intent(in) :: vortex_mpi(:,:),qeulerp_mpi(:),qeulerm_mpi(:),kt_mpi(:), &
      & E_mpi(:),E_sq_mpi(:),M_sq_mpi(:),M(:,:)
! internal variable
       integer :: i
       real(kind=8) :: cell,C,chi,E_av,steps,qp,qm,vortex(3),E_err,M_err,M_local(3)

       write(6,'(a)') "----------------"
       write(6,'(a)') "write in EM.dat before changing the temperature set (if applicable)"
       write(6,'(a)') "----------------"

       cell=dble(N_cell)
       steps=dble(n_thousand)

       do i=1,isize/n_ghost

        M_local=M(:,i)/steps/cell
        C=(E_sq_mpi(i)/steps-(E_mpi(i)/steps)**2)/kT_mpi(i)**2/cell
        chi=(M_sq_mpi(i)/steps-(norm(M(:,i))/steps)**2)/kT_mpi(i)/cell
        E_av=E_mpi(i)/cell/steps
        E_err=sqrt(abs(E_sq_mpi(i)-E_mpi(i)**2/steps)/(steps-1.0d0))/dble(N_cell)
        M_err=sqrt(abs(M_sq_mpi(i)-norm(M(:,i))**2/steps)/(steps-1.0d0))/dble(N_cell)
        qp=qeulerp_mpi(i)/pi/4.0d0/steps
        qm=qeulerm_mpi(i)/pi/4.0d0/steps
        vortex=vortex_mpi(:,i)/3.0d0/sqrt(3.0d0)/steps

        Write(7,'(15(E20.10E3,2x),E20.10E3)') kt_mpi(i)/k_B,E_av,E_err,C, &
     & norm(M_local(:)), M_local(:),M_err,chi,vortex,qp+qm,qp,qm
       enddo

        Write(7,'(a)') ' '

       end subroutine write_EM_serial

       end module
