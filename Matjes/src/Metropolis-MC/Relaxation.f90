module m_relaxation
implicit none
contains
!
! ===============================================================
SUBROUTINE Relaxation(lat,N_cell,n_sizerelax,n_relaxation,T_relax,E_total,E,Magnetization,qeulerp,qeulerm,kt,acc,rate,nb,cone,ising,equi,overrel,sphere,underrel,print_relax,Hams)
use mtprng
use m_Corre
use m_constants, only : k_b
use m_topocharge_all
use m_store_relaxation
use m_derived_types, only : point_shell_Operator,lattice
use m_basic_types, only : vec_point
use m_modes_variables, only : point_shell_mode
use m_H_public
use m_topo_commons
use m_io_utils
use m_io_files_utils
use m_convert
use m_MCstep
#ifdef CPP_MPI
      use m_mpi_prop, only : isize,irank_working,MPI_COMM,MPI_COMM_BOX,irank_box
#endif
Implicit none
! input
type(lattice),intent(inout)    :: lat
!type(vec_point), intent(inout) :: mode(N_cell)
real(kind=8), intent(inout) :: qeulerp,qeulerm,cone,acc,rate,E_total,magnetization(3),E(8)
real(kind=8), intent(inout) :: nb
real(kind=8), intent(in) :: kT
integer, intent(in) :: n_relaxation,T_relax,N_cell,n_sizerelax
logical, intent(in) :: ising,equi,overrel,sphere,underrel,print_relax
class(t_H), intent(in) :: Hams(:)
! a big table
real(kind=8) :: Relax(18,n_sizerelax),dumy(5)
!     Slope Index
integer :: i_relaxation,i_MC,io_relax
! dummy
integer :: i,j,n_w_step
character(len=50) :: fname

#ifdef CPP_MPI

      include 'mpif.h'

      if (irank_working.eq.0) write(6,'(/,a,/)') "starting relaxation"
#else
      write(6,'(/,a,/)') "starting relaxation"
#endif

Relax=0.0d0
n_w_step=n_relaxation/n_sizerelax

!     Monte Carlo steps for thermal equilibrium
!     The first time it will take longer
!     -----------------------------------------------------------------

do i_relaxation=1,n_relaxation
!         T_relax_1 is probably larger then T_relax, this is because
!         the last step might take more time to relax than an
!         the step to an unordered structure
            !Relaxation of the System
   Do i_MC=1,T_relax*N_cell
      Call MCStep(lat,N_cell,E_total,E,Magnetization,kt,acc,rate,nb,cone,ising,equi,overrel,sphere,underrel,Hams)
   enddo

            !In case T_relax set to zero at least one MCstep is done
   Call MCStep(lat,N_cell,E_total,E,Magnetization,kt,acc,rate,nb,cone,ising,equi,overrel,sphere,underrel,Hams)

! calculate the topocharge
   dumy=get_charge()
   qeulerp=dumy(1)
   qeulerm=dumy(2)


! Write the Equilibrium files

   if (print_relax) then
#ifdef CPP_MPI
      if (i_ghost) then

         E_int=E
         M_int=Magnetization
         qp_int=qeulerp
         qm_int=qeulerm
         qeulerm=0.0d0
         qeulerp=0.0d0
         Magnetization=0.0d0
         E=0.0d0

         call mpi_reduce(E,E_int,8,MPI_REAL8,MPI_SUM,0,MPI_COMM_BOX,ierr)
         call mpi_reduce(Magnetization,M_int,3,MPI_REAL8,MPI_SUM,0,MPI_COMM_BOX,ierr)
         call mpi_reduce(qeulerp,qp_int,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_BOX,ierr)
         call mpi_reduce(qeulerm,qm_int,1,MPI_REAL8,MPI_SUM,0,MPI_COMM_BOX,ierr)

      endif
#endif

      if (mod(i_relaxation,n_w_step).eq.0) call store_relaxation(Relax,i_relaxation,dble(i_relaxation), &
             &  sum(E)/dble(N_cell),E,dble(N_cell),kt,Magnetization,rate,cone,qeulerp,qeulerm)

   endif

enddo   ! enddo over the relaxation loop
!!!!!!!!!!!!!***************************!!!!!!
!!!!!!!!!!!!!***************************!!!!!!

write(6,'(/,a,f8.4,2x,a,/)') 'System is relaxed for T= ',kT/k_B,'Kelvin'

!print the Equilibrium files
if (print_relax) then

! numbering of the files
    fname=convert('Equilibriumphi-',kT/k_B,'.dat')
    io_relax=open_file_write(fname)

    Write(io_relax,'(a)') "#   1:n_MC  2:Etot  3:T  4:M  5:rate  6:cone  7:Q  8:Exch  9:Zeeman  10:Ani &
   &         11:4S  12:DM  13:biq  14:dip  15:stoner  16:Chi(E)(t,T)  17:Chi(M)(t,T)  18:Chi(T)(t)"

    Relax(16,:)=correlation(Relax(2,:),n_sizerelax)
    Relax(17,:)=correlation(Relax(4,:),n_sizerelax)
    Relax(18,:)=correlation(Relax(7,:),n_sizerelax)

    Do i=1,n_sizerelax
         Write(io_relax,'(i10,18(2x,E20.10E3))') int(Relax(1,i)),(Relax(j,i),j=2,18)
    enddo

    call close_file(fname,io_relax)

    write(6,'(/,a,f8.4,2x,a,/)') 'Equilibrium files are written for T= ',kT/k_B,'Kelvin'

endif


END SUBROUTINE Relaxation
end module
