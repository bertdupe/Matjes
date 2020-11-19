      module m_gather_reduce
      interface mpi_average
       module procedure mpi_average_ghost
       module procedure mpi_average_noghost
       module procedure mpi_average_gra
      end interface mpi_average
      interface end_reduce
       module procedure end_reduce_noghost
       module procedure end_reduce_ghost
      end interface end_reduce
      interface end_gather
       module procedure gather
       module procedure gather_separate
      end interface end_gather
      contains

! function that reduces all the quantities in the code
!---------------------------------
      subroutine end_reduce_noghost(vortex_mpi,qeulerp_mpi,qeulerm_mpi,Im_mpi,Qim_mpi,mpi_l, &
     & E_av,C,M_err,E_err,chi,M,isize,MPI_COMM)
      use m_rw_lattice, only : dim_lat
      use m_parameters, only : n_Tsteps
      use m_mpi
      implicit none
!intent(in)
      integer, intent(in) :: isize,MPI_COMM
!intent(inout)
      real(kind=8), intent(inout) :: vortex_mpi(:,:),qeulerp_mpi(:)
      real(kind=8), intent(inout) :: qeulerm_mpi(:),Im_mpi(:,:)
      real(kind=8), intent(inout) :: Qim_mpi(:,:),mpi_l(:,:)
      real(kind=8), intent(inout) :: E_av(:),C(:),M_err(:),E_err(:),chi(:),M(:,:)

! reduce the scalar
       E_av=reduce(E_av,n_Tsteps,MPI_COMM)/dble(isize)
       qeulerp_mpi=reduce(qeulerp_mpi,n_Tsteps,MPI_COMM)/dble(isize)
       qeulerm_mpi=reduce(qeulerm_mpi,n_Tsteps,MPI_COMM)/dble(isize)
       chi=reduce(chi,n_Tsteps,MPI_COMM)/dble(isize)
       C=reduce(C,n_Tsteps,MPI_COMM)/dble(isize)
       M_err=reduce(M_err,n_Tsteps,MPI_COMM)/dble(isize)
       E_err=reduce(E_err,n_Tsteps,MPI_COMM)/dble(isize)

! reduce the vector
      vortex_mpi=reduce(vortex_mpi,4,n_Tsteps,MPI_COMM)/dble(isize)
      M=reduce(M,3,n_Tsteps,MPI_COMM)/dble(isize)
      Im_mpi(:,1)=reduce(Im_mpi(:,1),n_Tsteps,MPI_COMM)/dble(isize)
      Qim_mpi(:,1)=reduce(Qim_mpi(:,1),n_Tsteps,MPI_COMM)/dble(isize)
      mpi_l=reduce(mpi_l,3,n_Tsteps,MPI_COMM)/dble(isize)

      end subroutine end_reduce_noghost

! function that reduces all the quantities in the code
!---------------------------------
      subroutine end_reduce_ghost(vortex_mpi,qeulerp_mpi,qeulerm_mpi,Im_mpi,Qim_mpi,mpi_l, &
     & E_av,M_err,E_err,M,isize,MPI_COMM)
      use m_rw_lattice, only : dim_lat
      use m_parameters, only : n_Tsteps
      use m_mpi
      implicit none
!intent(in)
      integer, intent(in) :: isize,MPI_COMM
!intent(inout)
      real(kind=8), intent(inout) :: vortex_mpi(:),qeulerp_mpi
      real(kind=8), intent(inout) :: qeulerm_mpi,Im_mpi
      real(kind=8), intent(inout) :: Qim_mpi,mpi_l(:)
      real(kind=8), intent(inout) :: E_av,M_err,E_err,M(:)

! reduce the scalar
       E_av=reduce(E_av,MPI_COMM)
       qeulerp_mpi=reduce(qeulerp_mpi,MPI_COMM)
       qeulerm_mpi=reduce(qeulerm_mpi,MPI_COMM)
       M_err=reduce(M_err,MPI_COMM)
       E_err=reduce(E_err,MPI_COMM)

! reduce the vector
      vortex_mpi=reduce(vortex_mpi,4,MPI_COMM)
      M=reduce(M,3,MPI_COMM)
      Im_mpi=reduce(Im_mpi,MPI_COMM)
      Qim_mpi=reduce(Qim_mpi,MPI_COMM)
      mpi_l=reduce(mpi_l,3,MPI_COMM)

      end subroutine end_reduce_ghost

! function that gather all the quantities in the code
! in this version, each proc calculates the quantities for one temperature
!---------------------------------
      subroutine gather(shape_masque,vortex_mpi,qeulerp_mpi,qeulerm_mpi,Im_mpi,Qim_mpi,kt_mpi,mpi_l, &
     & E_av,C,M_err,E_err,chi,M,isize,MPI_COMM)
      use m_rw_lattice, only : dim_lat
      use m_parameters, only : n_Tsteps,kt
      use m_mpi
      implicit none
! intent(in)
      integer, intent(in) :: isize,MPI_COMM
! intent(inout)
!      real(kind=8), intent(inout) :: vortex_mpi(4,n_Tsteps),qeulerp_mpi(n_Tsteps)
!      real(kind=8), intent(inout) :: qeulerm_mpi(n_Tsteps),Im_mpi(n_Tsteps,dim_lat(3))
!      real(kind=8), intent(inout) :: Qim_mpi(n_Tsteps,dim_lat(3)),kt_mpi(n_Tsteps)
!      real(kind=8), intent(inout) :: mpi_l(3,n_Tsteps)
!      real(kind=8), intent(inout) :: E_av(n_Tsteps),C(n_Tsteps),M_err(n_Tsteps),E_err(n_Tsteps),chi(n_Tsteps),M(3,n_Tsteps)
!      integer, intent(in) :: shape_tableNN(6)
      real(kind=8), intent(inout) :: vortex_mpi(:,:),qeulerp_mpi(:)
      real(kind=8), intent(inout) :: qeulerm_mpi(:),Im_mpi(:,:)
      real(kind=8), intent(inout) :: Qim_mpi(:,:),kt_mpi(:)
      real(kind=8), intent(inout) :: mpi_l(:,:)
      real(kind=8), intent(inout) :: E_av(:),C(:),M_err(:),E_err(:),chi(:),M(:,:)
      integer, intent(in) :: shape_masque(:)
!!dumy
      integer :: k,ierr,taille
      real(kind=8) :: toto,toto_4(4),toto_IQ(dim_lat(3)),Nsite

      include 'mpif.h'

       Nsite=dble(shape_masque(2)*shape_masque(3)*shape_masque(4))

! reduce the scalar
       toto=E_av(1)/Nsite
       call mpi_gather(toto,1,MPI_REAL8,E_av,1,MPI_REAL8,0,MPI_COMM,ierr)
       toto=qeulerp_mpi(1)
       call mpi_gather(toto,1,MPI_REAL8,qeulerp_mpi,1,MPI_REAL8,0,MPI_COMM,ierr)
       toto=qeulerm_mpi(1)
       call mpi_gather(toto,1,MPI_REAL8,qeulerm_mpi,1,MPI_REAL8,0,MPI_COMM,ierr)
       toto=M_err(1)
       call mpi_gather(toto,1,MPI_REAL8,M_err,1,MPI_REAL8,0,MPI_COMM,ierr)
       toto=E_err(1)
       call mpi_gather(toto,1,MPI_REAL8,E_err,1,MPI_REAL8,0,MPI_COMM,ierr)
       toto=chi(1)
       call mpi_gather(toto,1,MPI_REAL8,chi,1,MPI_REAL8,0,MPI_COMM,ierr)
       toto=C(1)
       call mpi_gather(toto,1,MPI_REAL8,C,1,MPI_REAL8,0,MPI_COMM,ierr)
       call mpi_gather(kt,1,MPI_REAL8,kt_mpi,1,MPI_REAL8,0,MPI_COMM,ierr)
! reduce vector
       taille=size(vortex_mpi,1)
       toto_4(1:taille)=vortex_mpi(:,1)
       call mpi_gather(toto_4,taille,MPI_REAL8,vortex_mpi(:,:),taille,MPI_REAL8,0,MPI_COMM,ierr)
       toto_4(1:3)=M(:,1)/Nsite
       call mpi_gather(toto_4,3,MPI_REAL8,M(:,:),3,MPI_REAL8,0,MPI_COMM,ierr)
       taille=size(Qim_mpi,2)
       toto_IQ=Qim_mpi(1,:)
       call mpi_gather(toto_IQ,taille,MPI_REAL8,Qim_mpi(:,:),taille,MPI_REAL8,0,MPI_COMM,ierr)
       taille=size(Im_mpi,2)
       toto_IQ=Im_mpi(1,:)
       call mpi_gather(toto_IQ,taille,MPI_REAL8,Im_mpi(:,:),taille,MPI_REAL8,0,MPI_COMM,ierr)
       taille=size(mpi_l,1)
       toto_4(1:3)=mpi_l(:,1)
       call mpi_gather(toto_4,taille,MPI_REAL8,mpi_l(:,:),taille,MPI_REAL8,0,MPI_COMM,ierr)

      end subroutine gather

!---------------------------------

! function that gather all the quantities in the code
! in this version, each proc calculates the quantities for one temperature
!---------------------------------
      subroutine gather_separate(kT_all,E_av,E_err_av,C_av,M_sum_av,M_err_av,chi_M,vortex_av,chi_Q,qeulerp_av,qeulerm_av,chi_l, &
                          & size_table,irank,isize,MPI_COMM)
      use m_mpi
      implicit none
! intent(in)
      integer, intent(in) :: isize,MPI_COMM,size_table,irank
! intent(inout)
      real(kind=8), intent(inout) :: M_sum_av(:,:),M_err_av(:,:),chi_M(:,:),vortex_av(:,:),chi_l(:,:),chi_Q(:,:)
      real(kind=8), intent(inout) :: kT_all(:),E_av(:),E_err_av(:),C_av(:),qeulerp_av(:),qeulerm_av(:)
!!dumy
      integer :: ierr,istart,istop
      real(kind=8) :: transfer_1D(size_table),transfer_2D(3,size_table)

      include 'mpif.h'

      transfer_1D=0.0d0
      transfer_2D=0.0d0

      istart=1+(irank*isize)
      istop=((irank+1)*isize)

! reduce the 1 dimensional vector
!
      transfer_1D(istart:istop)=kT_all(1:isize)
      kT_all=reduce(transfer_1D,size_table,MPI_COMM)
      transfer_1D=0.0d0

      transfer_1D(istart:istop)=E_av(1:isize)
      E_av=reduce(transfer_1D,size_table,MPI_COMM)
      transfer_1D=0.0d0

      transfer_1D(istart:istop)=E_err_av(1:isize)
      E_err_av=reduce(transfer_1D,size_table,MPI_COMM)
      transfer_1D=0.0d0

      transfer_1D(istart:istop)=C_av(1:isize)
      C_av=reduce(transfer_1D,size_table,MPI_COMM)
      transfer_1D=0.0d0

      transfer_1D(istart:istop)=qeulerp_av(1:isize)
      qeulerp_av=reduce(transfer_1D,size_table,MPI_COMM)
      transfer_1D=0.0d0

      transfer_1D(istart:istop)=qeulerm_av(1:isize)
      qeulerm_av=reduce(transfer_1D,size_table,MPI_COMM)


! reduce the 2 dimensional vector
!
      transfer_2D(:,istart:istop)=M_sum_av(:,1:isize)
      M_sum_av=reduce(transfer_2D,3,size_table,MPI_COMM)
      transfer_2D=0.0d0

      transfer_2D(:,istart:istop)=M_err_av(:,1:isize)
      M_err_av=reduce(transfer_2D,3,size_table,MPI_COMM)
      transfer_2D=0.0d0

      transfer_2D(:,istart:istop)=chi_M(:,1:isize)
      chi_M=reduce(transfer_2D,3,size_table,MPI_COMM)
      transfer_2D=0.0d0

      transfer_2D(:,istart:istop)=vortex_av(:,1:isize)
      vortex_av=reduce(transfer_2D,3,size_table,MPI_COMM)
      transfer_2D=0.0d0

      transfer_2D(:,istart:istop)=chi_l(:,1:isize)
      chi_l=reduce(transfer_2D,3,size_table,MPI_COMM)
      transfer_2D=0.0d0

      transfer_2D(:,istart:istop)=chi_Q(:,1:isize)
      chi_Q=reduce(transfer_2D,3,size_table,MPI_COMM)
      transfer_2D=0.0d0

      end subroutine gather_separate

!---------------------------------
      subroutine mpi_average_gra(n_kt,vortex_mpi,qeulerp_mpi,qeulerm_mpi,map_vort,mapvort_mpi &
     &  ,mapeul_mpi,map_toto,Im_mpi,Qim_mpi,Im,Qim,qeulerm,qeulerp,vortex, &
     &  E_sum,E_sq_sum,M_sum,M_sq_sum,C,chi,E_err,M,E_av,M_err,N_site)
      use m_constants, only : pi
#ifdef CPP_MPI
      use m_parameters, only : n_Tsteps,kt,Total_MC_Steps,i_average
      use m_mpi_prop, only : irank,irank_box,isize,MPI_COMM
      use m_mpi
#else
      use m_parameters, only : n_Tsteps,kt,Total_MC_Steps
#endif
      use m_rw_lattice, only : dim_lat
      use m_vector, only : norm
      implicit none
!intent(inout)
      real(kind=8), intent(inout) :: vortex_mpi(:,:),qeulerp_mpi(:)
      real(kind=8), intent(inout) :: qeulerm_mpi(:),mapvort_mpi(:,:,:,:,:)
      real(kind=8), intent(inout) :: mapeul_mpi(:,:,:,:),Im_mpi(:,:)
      real(kind=8), intent(inout) :: Qim_mpi(:,:)
      real(kind=8), intent(inout) :: C(:),chi(:),E_err(:),M(:,:),E_av(:),M_err(:)
!intent(in)
      integer, intent(in) :: n_kt
      real(kind=8), intent(in) :: map_vort(:,:,:,:),map_toto(:,:,:)
      real(kind=8), intent(in) :: qeulerp,qeulerm,Im(:),Qim(:),N_site
      real(kind=8), intent(in) :: vortex(:)
      real(kind=8), intent(in) :: E_sum,E_sq_sum,M_sum(:),M_sq_sum
! dummy
      real(kind=8) :: MC_Steps
      integer :: i_x,i_y,i_z

      MC_Steps=dble(Total_MC_Steps)

#ifdef CPP_MPI
      if (i_average) then
       MC_Steps=dble(Total_MC_Steps*isize)
       C(n_kT)=(allreduce(E_sq_sum,MPI_COMM)/MC_Steps-(allreduce(E_sum,MPI_COMM)/MC_Steps)**2)/kT**2/N_site
       chi(n_kT)=(allreduce(M_sq_sum,MPI_COMM)/MC_Steps-norm(allreduce(M_sum,3,MPI_COMM))**2/MC_Steps**2)/kT/N_site
       E_err(n_kT)=dsqrt(dabs(allreduce(E_sq_sum,MPI_COMM)-allreduce(E_sum**2,MPI_COMM)/MC_Steps)/(MC_Steps-1.0d0))/N_site
       M_err(n_kT)=dsqrt(dabs(allreduce(M_sq_sum,MPI_COMM)-norm(allreduce(M_sum,3,MPI_COMM))**2/MC_Steps)/(MC_Steps-1.0d0))/N_site
       E_av(n_kT)=allreduce(E_sum,MPI_COMM)/N_site/MC_Steps
       M(:,n_kT)=allreduce(M_sum,3,MPI_COMM)/N_site/MC_Steps
       qeulerp_mpi(n_kT)=allreduce(qeulerp,MPI_COMM)/MC_Steps/pi/4.0d0
       qeulerm_mpi(n_kT)=allreduce(qeulerm,MPI_COMM)/MC_Steps/pi/4.0d0
       vortex_mpi(:,n_kT)=allreduce(vortex,4,MPI_COMM)/MC_Steps/3.0d0/dsqrt(3.0d0)
       Im_mpi(n_kT,:)=allreduce(Im,dim_lat(3),MPI_COMM)/MC_Steps
       Qim_mpi(n_kT,:)=allreduce(Qim,dim_lat(3),MPI_COMM)/MC_Steps
      else
       C(n_kT)=(E_sq_sum/MC_Steps-(E_sum/MC_Steps)**2)/kT**2/N_site
       chi(n_kT)=(M_sq_sum/MC_Steps-(norm(M_sum)/MC_Steps)**2)/kT/N_site
       E_err(n_kT)=sqrt(abs(E_sq_sum-E_sum**2/MC_Steps)/(MC_Steps-1.0d0))/N_site
       M_err(n_kT)=sqrt(abs(M_sq_sum-norm(M_sum)**2/MC_Steps)/(MC_Steps-1.0d0))/N_site
       E_av(n_kT)=E_sum/N_site/MC_Steps
       M(:,n_kT)=M_sum/N_site/MC_Steps
       qeulerp_mpi(n_kT)=qeulerp/MC_Steps/pi/4.0d0
       qeulerm_mpi(n_kT)=qeulerm/MC_Steps/pi/4.0d0
       vortex_mpi(:,n_kT)=vortex/MC_Steps/3.0d0/dsqrt(3.0d0)
       Im_mpi(n_kT,:)=Im/MC_Steps
       Qim_mpi(n_kT,:)=Qim/MC_Steps
      endif
#else
       C(n_kT)=(E_sq_sum/MC_Steps-(E_sum/MC_Steps)**2)/kT**2/N_site
       chi(n_kT)=(M_sq_sum/MC_Steps-(norm(M_sum)/MC_Steps)**2)/kT/N_site
       E_err(n_kT)=sqrt(abs(E_sq_sum-E_sum**2/MC_Steps)/(MC_Steps-1.0d0))/N_site
       M_err(n_kT)=sqrt(abs(M_sq_sum-norm(M_sum)**2/MC_Steps)/(MC_Steps-1.0d0))/N_site
       E_av(n_kT)=E_sum/N_site/MC_Steps
       M(:,n_kT)=M_sum/N_site/MC_Steps
       qeulerp_mpi(n_kT)=qeulerp/MC_Steps/pi/4.0d0
       qeulerm_mpi(n_kT)=qeulerm/MC_Steps/pi/4.0d0
       vortex_mpi(:,n_kT)=vortex/MC_Steps/3.0d0/dsqrt(3.0d0)
       Im_mpi(n_kT,:)=Im/MC_Steps
       Qim_mpi(n_kT,:)=Qim/MC_Steps
#endif

#ifdef CPP_MPI
      if (irank.eq.0) then
#endif
      do i_z=1,dim_lat(3)
       do i_y=1,dim_lat(2)
        do i_x=1,dim_lat(1)
        mapvort_mpi(i_x,i_y,i_z,:,n_kT)=map_vort(i_x,i_y,i_z,:)/MC_Steps/3.0d0/dsqrt(3.0d0)
        mapeul_mpi(i_x,i_y,i_z,n_kT)=map_toto(i_x,i_y,i_z)/MC_Steps/pi/4.0d0
        enddo
       enddo
      enddo
#ifdef CPP_MPI
      endif
#endif

      end subroutine mpi_average_gra

!---------------------------------
      subroutine mpi_average_ghost(vortex_mpi,qeulerp_mpi,qeulerm_mpi &
     &  ,Im_mpi,Qim_mpi,Im,Qim,qeulerm,qeulerp,vortex, &
     &  E_sum,E_sq_sum,M_sum,M_sq_sum,C,chi,E_err,M,E_av,M_err,isize,MPI_COMM)
      use m_constants, only : pi
      use m_parameters, only : n_Tsteps,kt,Total_MC_Steps
      use m_vector, only : norm
      use m_rw_lattice, only : dim_lat
      use m_mpi
      implicit none
!intent(inout)
      real(kind=8), intent(inout) :: vortex_mpi(:),qeulerp_mpi
      real(kind=8), intent(inout) :: qeulerm_mpi
      real(kind=8), intent(inout) :: Im_mpi(:),Qim_mpi(:)
      real(kind=8), intent(inout) :: C,chi,E_err,M(:),E_av,M_err
!intent(in)
      integer, intent(in) :: isize,MPI_COMM
      real(kind=8), intent(in) :: qeulerp,qeulerm,Im(:),Qim(:)
      real(kind=8), intent(in) :: vortex(:)
      real(kind=8), intent(in) :: E_sum,E_sq_sum,M_sum(:),M_sq_sum
! dummy
      real(kind=8) ::  MC_Steps

      MC_Steps=dble(Total_MC_Steps)

      C=(allreduce(E_sq_sum,MPI_COMM)/MC_Steps-(allreduce(E_sum,MPI_COMM)/MC_Steps)**2)/kT**2
      chi=(allreduce(M_sq_sum,MPI_COMM)/MC_Steps-norm(allreduce(M_sum,3,MPI_COMM))**2/MC_Steps**2)/kT
      E_err=sqrt(dabs(allreduce(E_sq_sum,MPI_COMM)-allreduce(E_sum**2,MPI_COMM)/MC_Steps)/(MC_Steps-1.0d0))
      M_err=sqrt(dabs(allreduce(M_sq_sum,MPI_COMM)-norm(allreduce(M_sum,3,MPI_COMM))**2/MC_Steps)/(MC_Steps-1.0d0))
      E_av=allreduce(E_sum,MPI_COMM)/MC_Steps
      M(:)=allreduce(M_sum,3,MPI_COMM)/MC_Steps
      qeulerp_mpi=allreduce(qeulerp,MPI_COMM)/MC_Steps/pi/4.0d0
      qeulerm_mpi=allreduce(qeulerm,MPI_COMM)/MC_Steps/pi/4.0d0
      vortex_mpi(:)=allreduce(vortex,4,MPI_COMM)/MC_Steps/3.0d0/dsqrt(3.0d0)
      Im_mpi(:)=allreduce(Im,dim_lat(3),MPI_COMM)/MC_Steps
      Qim_mpi(:)=allreduce(Qim,dim_lat(3),MPI_COMM)/MC_Steps

      end subroutine mpi_average_ghost

      !---------------------------------
      subroutine mpi_average_noghost(n_kt,vortex_mpi,qeulerp_mpi,qeulerm_mpi &
     &  ,Im_mpi,Qim_mpi,Im,Qim,qeulerm,qeulerp,vortex, &
     &  E_sum,E_sq_sum,M_sum,M_sq_sum,C,chi,E_err,M,E_av,M_err,N_site)
      use m_constants, only : pi
#ifdef CPP_MPI
      use m_parameters, only : n_Tsteps,kt,Total_MC_Steps,i_average
      use m_mpi_prop, only : isize,MPI_COMM
#else
      use m_parameters, only : n_Tsteps,kt,Total_MC_Steps
#endif
      use m_rw_lattice, only : dim_lat
      use m_lattice, only : masque
      use m_vector, only : norm
      use m_mpi
      implicit none
!intent(inout)
      real(kind=8), intent(inout) :: vortex_mpi(:,:),qeulerp_mpi(:)
      real(kind=8), intent(inout) :: qeulerm_mpi(:)
      real(kind=8), intent(inout) :: Im_mpi(:,:),Qim_mpi(:,:)
      real(kind=8), intent(inout) :: C(:),chi(:),E_err(:),M(:,:),E_av(:),M_err(:)
!intent(in)
      integer, intent(in) :: n_kt
      real(kind=8), intent(in) :: qeulerp,qeulerm,Im(:),Qim(:),N_site
      real(kind=8), intent(in) :: vortex(:)
      real(kind=8), intent(in) :: E_sum,E_sq_sum,M_sum(:),M_sq_sum
! dummy
      real(kind=8) ::  MC_Steps

      MC_Steps=dble(Total_MC_Steps)

#ifdef CPP_MPI
      if (i_average) then
       MC_Steps=dble(Total_MC_Steps)
       C(n_kT)=(allreduce(E_sq_sum,MPI_COMM)/MC_Steps-(allreduce(E_sum,MPI_COMM)/MC_Steps)**2)/kT**2/N_site
       chi(n_kT)=(allreduce(M_sq_sum,MPI_COMM)/MC_Steps-norm(allreduce(M_sum,3,MPI_COMM))**2/MC_Steps**2)/kT/N_site
       E_err(n_kT)=sqrt(dabs(allreduce(E_sq_sum,MPI_COMM)-allreduce(E_sum**2,MPI_COMM)/MC_Steps)/(MC_Steps-1.0d0))/N_site
       M_err(n_kT)=sqrt(dabs(allreduce(M_sq_sum,MPI_COMM)-norm(allreduce(M_sum,3,MPI_COMM))**2/MC_Steps)/(MC_Steps-1.0d0))/N_site
       E_av(n_kT)=allreduce(E_sum,MPI_COMM)/MC_Steps
       M(:,n_kT)=allreduce(M_sum,3,MPI_COMM)/MC_Steps
       qeulerp_mpi(n_kT)=allreduce(qeulerp,MPI_COMM)/MC_Steps/pi/4.0d0
       qeulerm_mpi(n_kT)=allreduce(qeulerm,MPI_COMM)/MC_Steps/pi/4.0d0
       vortex_mpi(:,n_kT)=allreduce(vortex,4,MPI_COMM)/MC_Steps/3.0d0/dsqrt(3.0d0)
       Im_mpi(n_kT,:)=allreduce(Im,dim_lat(3),MPI_COMM)/MC_Steps
       Qim_mpi(n_kT,:)=allreduce(Qim,dim_lat(3),MPI_COMM)/MC_Steps
      else
       C(n_kT)=(E_sq_sum/MC_Steps-(E_sum/MC_Steps)**2)/kT**2/N_site
       chi(n_kT)=(M_sq_sum/MC_Steps-(norm(M_sum)/MC_Steps)**2)/kT/N_site
       E_err(n_kT)=sqrt(dabs(E_sq_sum-E_sum**2/MC_Steps)/(MC_Steps-1.0d0))/N_site
       M_err(n_kT)=sqrt(dabs(M_sq_sum-norm(M_sum)**2/MC_Steps)/(MC_Steps-1.0d0))/N_site
       E_av(n_kT)=E_sum/MC_Steps
       M(:,n_kT)=M_sum/MC_Steps
       qeulerp_mpi(n_kT)=qeulerp/MC_Steps/pi/4.0d0
       qeulerm_mpi(n_kT)=qeulerm/MC_Steps/pi/4.0d0
       vortex_mpi(:,n_kT)=vortex/MC_Steps/3.0d0/dsqrt(3.0d0)
       Im_mpi(n_kT,:)=Im/MC_Steps
       Qim_mpi(n_kT,:)=Qim/MC_Steps
      endif
#else
       C(n_kT)=(E_sq_sum/MC_Steps-(E_sum/MC_Steps)**2)/kT**2/N_site
       chi(n_kT)=(M_sq_sum/MC_Steps-(norm(M_sum)/MC_Steps)**2)/kT/N_site
       E_err(n_kT)=sqrt(dabs(E_sq_sum-E_sum**2/MC_Steps)/(MC_Steps-1.0d0))/N_site
       M_err(n_kT)=sqrt(dabs(M_sq_sum-norm(M_sum)**2/MC_Steps)/(MC_Steps-1.0d0))/N_site
       E_av(n_kT)=E_sum/MC_Steps
       M(:,n_kT)=M_sum/MC_Steps
       qeulerp_mpi(n_kT)=qeulerp/MC_Steps/pi/4.0d0
       qeulerm_mpi(n_kT)=qeulerm/MC_Steps/pi/4.0d0
       vortex_mpi(:,n_kT)=vortex/MC_Steps/3.0d0/dsqrt(3.0d0)
       Im_mpi(n_kT,:)=Im/MC_Steps
       Qim_mpi(n_kT,:)=Qim/MC_Steps
#endif

      end subroutine mpi_average_noghost

      end module m_gather_reduce
!------------------------------
