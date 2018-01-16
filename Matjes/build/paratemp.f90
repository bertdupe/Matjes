     module m_paratemp
       interface paratemp
          module procedure serial_paratemp
       end interface paratemp

       interface calculate_diffusion
          module procedure calculate_diffusion_serial
       end interface calculate_diffusion
#ifdef CPP_MPI
       interface paratemp_gather
          module procedure paratemp_gather_mpi
       end interface paratemp_gather

       interface paratemp_scatter
          module procedure paratemp_scatter_mpi
       end interface paratemp_scatter
#endif
       contains

!--------------------------------------------
! parallel tempering as described in Katzgraber et al PRE and cond-mat 30 mars 2006 arXiv:cond-mat/0602085v3

!--------------------------------------------
!--------------------------------------------
       subroutine calculate_diffusion_serial(kT_all,kt_updated,nup,ndown,Nsuccess,n_thousand,isize,i_optTset)
       use m_fit
       use m_constants, only : k_B
#ifdef CPP_MPI
       use m_mpi_prop, only : MPI_COMM
#endif
       implicit none
       integer, intent(in) :: n_thousand,isize
       real(kind=8), intent(in) :: nup(:),ndown(:),Nsuccess(:),kT_all(:)
       logical, intent(in) :: i_optTset
       real(kind=8), intent(inout) :: kT_updated(:)
!      internals
! save the fraction
       real(kind=8) :: frac(isize),order_frac(2,isize)
! dT for the integration
       real(kind=8) :: deltaT(isize)
! derivative if the fraction
       real(kind=8) :: Dfrac(isize-1)
! diffusivity of the images
       real(kind=8) :: diffusivity(isize-1)
! integral and all the necessary stuff
       real(kind=8) :: this_dt,integral,Kprime,tau
       real(kind=8) :: etaTprim(isize-1)
! fit of the fraction
       real(kind=8) :: fit_frac(2)
! number of temperatures actually usefull
       integer :: Ntemp_used
! counters
       integer :: i,j,l

       frac=0.0d0
       order_frac=-1.0d0
       etaTprim=0.0d0

#ifdef CPP_DEBUG
       write(6,*) 'nup and ndown'
       write(6,*) nup
       write(6,*) '-------------'
       write(6,*) ndown
#endif


!calculation of the fraction of images going up
       do i=1,isize
           if (nup(i)+ndown(i).ge.1.0d0) then
               frac(i)=nup(i)/(nup(i)+ndown(i))
           else
               write(6,'(a)') 'some up+down fractions are 0'
               write(6,'(a)') 'try to reduce the temperature range to favor overlapping between replicas'
               write(6,*) (frac(j),j=1,i)
               write(6,*) (nup(j),j=1,i)
               write(6,*) (ndown(j),j=1,i)
#ifdef CPP_MPI
               if (i_optTset) call mpi_abort(MPI_COMM,errcode,ierr)
#else
               if (i_optTset) stop
#endif
           endif
       enddo

       fit_frac=fit(frac,isize)

! starting a lots of checks

       if ((frac(1)-1.0d0).gt.1.0d-8) then
           write(6,'(a)') 'the fraction is not 1 at lower temperature'
           write(6,'(a)') 'you can not optimise the temperature set in these conditions'
           write(6,'(a)') 'try to increase N_thousand'
#ifdef CPP_MPI
           if (i_optTset) call mpi_abort(MPI_COMM,errcode,ierr)
#else
           if (i_optTset) stop
#endif
       endif

       if (frac(isize).gt.1.0d-8) then
           write(6,'(a)') 'the fraction is not 0 at higher temperature'
           write(6,'(a)') 'you can not optimise the temperature set in these conditions'
           write(6,'(a)') 'try to increase N_thousand'
#ifdef CPP_MPI
           if (i_optTset) call mpi_abort(MPI_COMM,errcode,ierr)
#else
           if (i_optTset) stop
#endif
       endif

! if the temperature is not optimized, write the fraction and go out
       if (.not.i_optTset) then
           do i=1,isize-1
               write(21,'(3(E20.12E3,2x))') kT_all(i)/k_B,frac(i),(Nsuccess(i)/n_thousand)
           enddo
           write(21,'(2(E20.12E3,2x))') kT_all(isize)/k_B,frac(isize)
           write(21,'(a)') ''
           kT_updated=kT_all
           return
       endif
!!!!!!!!!!!!!!!!!!!!! i_optTset part !!!!!!!!!!!!!!!

! sort out the fraction which are equal
       order_frac(1,1)=kT_all(1)
       order_frac(2,1)=frac(1)
       l=1
       do i=2,isize
           j=1
           do while (j.le.l)
              if ((order_frac(2,j).le.frac(i)).and.(order_frac(2,j).ne.frac(1))) then
                  order_frac(2,j)=frac(i)
                  order_frac(1,j)=kT_all(i)
                  l=j
                  write(6,'(a)') "WARNING: the fraction is increasing with temperature"
                  write(6,'(a)') "It should always decrease: increase n_thousand"
                  exit
              else
                  if ((j.eq.l).and.(order_frac(2,l).ne.frac(i))) then
                      l=l+1
                      order_frac(2,l)=frac(i)
                      order_frac(1,l)=kT_all(i)
                      exit
                  endif
              endif
              j=j+1
              if (order_frac(2,l).lt.1.0d-8) exit
           enddo
           if (order_frac(2,l).lt.1.0d-8) then
               order_frac(1,l)=kT_all(isize)
               exit
           endif
       enddo
! number of relevant temperatures
       Ntemp_used=l

#ifdef CPP_DEBUG
       do i=1,isize
           write(*,*) i,order_frac(2,i),frac(i)
       enddo
       write(*,*) '-------'
       do i=1,isize
           write(*,*) i,order_frac(1,i)/k_b,kT_all(i)/k_b
       enddo
       write(*,*) '-------'
#endif

! Calculation of DeltaT
       deltaT=0.0d0
       do i=1,Ntemp_used-1
          deltaT(i)=order_frac(1,i+1)-order_frac(1,i)
       enddo

! derivative of the fraction with a 3 point linear regression
! the idea is to get Dfrac<0 all the time. We are sure that the algo is correct because by construction
!F(T1)=1 and F(Tn)=0 and T1<Tn. The only problem is to chose enough point in the middle to allow the correct linear regression
       Dfrac=0.0d0
       do i=1,Ntemp_used-1
          Dfrac(i)=(order_frac(2,i+1)-order_frac(2,i))/deltaT(i)
       enddo

!calculation of the diffusivity
       diffusivity=0.0d0
       do i=1,Ntemp_used-1
          diffusivity(i)=-deltaT(i)/Dfrac(i)
       enddo

! calculate the renormalization constant. The number of temperature is not constant anymore
       Kprime=0.0d0
       do i=1,Ntemp_used-1
          Kprime=Kprime+deltaT(i)/sqrt(diffusivity(i))
       enddo
       Kprime=dble(isize-1)/Kprime

! calculate DTprime
       etaTprim=0.0d0
       do i=1,Ntemp_used-1
          etaTprim(i)=Kprime/sqrt(diffusivity(i))
       enddo

       kt_updated(1)=kT_all(1)
       kt_updated(isize)=kT_all(isize)

       l=1
       j=1
       integral=0.0d0
       this_dt=0.0d0
       do i=1,Ntemp_used-1
          do while ((integral.le.dble(l)).and.(l.lt.(isize-1)))
              integral=integral+deltaT(i)*etaTprim(i)
              if (integral.gt.dble(l)) then
                 integral=integral-deltaT(i)*etaTprim(i)
                 tau=(dble(l)-integral)/etaTprim(i)
                 if (tau.gt.(deltaT(i)+this_dt)) then
                    write(6,'(a)') 'inconsistent new delta T'
                    write(6,'(2(a,f8.4,2x))') 'delta T=',deltaT(i)/k_b,'tau=',tau/k_b
#ifdef CPP_MPI
                    call mpi_abort(MPI_COMM,errcode,ierr)
#else
                    stop
#endif
                 endif
              deltaT(i)=deltaT(i)-tau
              this_dt=this_dt+tau
              kt_updated(l+1)=this_dt+kt_updated(l)
              integral=integral+tau*etaTprim(i)
              l=l+1
              this_dt=0.0d0
              else
                 this_dt=this_dt+deltaT(i)
                 exit
              endif
          enddo
       enddo

       if (integral.ne.dble(isize-2)) then
          write(6,'(a)') "inconsistent new temperature step"
#ifdef CPP_MPI
          call mpi_abort(MPI_COMM,errcode,ierr)
#else
          stop
#endif
       endif

       do i=1,isize-1
            write(21,'(7(E20.12E3,2x))') kt_updated(i)/k_B,kT_all(i)/k_B,frac(i),diffusivity(i),Dfrac(i),Nsuccess(i)/dble(n_thousand),etaTprim(i)
       enddo

       write(21,'(3(E20.12E3,2x),a,/)') kt_updated(isize)/k_B,kT_all(isize)/k_B,frac(isize),'   #  ""  ""  ""  "" '

       end subroutine calculate_diffusion_serial





!--------------------------------------------
!--------------------------------------------
       subroutine serial_paratemp(state,label,nup,ndown,nstart,Nsuccess,kt_saved,image_temp,E_temp,isize,kt_all)
       use m_constants, only : k_B
       use mtprng
       implicit none
       integer, intent(in) :: nstart,isize
       real(kind=8), intent(in) :: kt_saved(:),kt_all(:)
       integer, intent(inout) :: image_temp(:)
       real(kind=8), intent(inout) :: nup(:),ndown(:)
       real(kind=8), intent(inout) :: label(:),Nsuccess(:)
       real(kind=8), intent(inout) :: E_temp(:)
       type(mtprng_state), intent(inout) :: state
!      internals
       integer :: iii
       real(kind=8) :: choice,delta
       real(kind=8) :: kTn,kti,Ei,En
       logical :: accept
       integer :: i,i_image,i_start


       accept=.False.

! count is running over the number of temperature

       do i=1,isize

#ifdef CPP_MPI
! since all the temperatures can be mixed and on different processors in MPI, one has to order the current depending on the
! temperature
! CAUTIOUS: the currents must match the temperatures order. it means that the loop is done on the replicas in MPI and on the
! temperature in serial
          kt=kt_all(i)
          i_image=minloc(abs(kt-kt_saved),1)
#else
! replica nb i_image with the temperature i
          i_image=image_temp(i)
#endif

          if (label(i_image).gt.0.5d0) nup(i)=nup(i)+1.0d0
          if (label(i_image).lt.-0.5d0) ndown(i)=ndown(i)+1.0d0
       enddo

! in case there are an odd number of replicas, change the start from 1 to 2
       i_start=mod(nstart,2)+1
! count is running over the number of temperatures
       do i=i_start,isize-1,2

#ifdef CPP_MPI
! cautious, in the case of MPI, the temperatures in kt_saved can be in a different order
          kti=kt_all(i)
          i_image=minloc(abs(kti-kt_saved),1)
          Ei=E_temp(i_image)

          ktn=kt_all(i+1)
          i_image=minloc(abs(ktn-kt_saved),1)
          En=E_temp(i_image)
#else
          kti=kt_saved(i)
          Ei=E_temp(i)

          ktn=kt_saved(i+1)
          En=E_temp(i+1)
#endif
          delta=(1/ktn-1/kti)*(En-Ei)

! exp(-100) is the maximum number that can be caculated it seams.
          if (delta.lt.-100.0d0) delta=-100.0d0
!                delta=1.0d0

          if (exp(delta).ge.1.0d0) then
             accept=.True.
          else
#ifdef CPP_MRG
             Choice=mtprng_rand_real1(state)
#else
             CALL RANDOM_NUMBER(Choice)
#endif
             if (exp(delta).gt.Choice) accept=.True.
          endif

          if (accept) then

#ifdef CPP_MPI
            iii=minloc(abs(kti-kt_saved),1)
            iin=minloc(abs(ktn-kt_saved),1)
            image_temp(iii)=iin
            image_temp(iin)=iii
#else
! position of the replica with the temperature i+1
            iii=image_temp(i+1)
! swap the temperature of replica i and i+1
            image_temp(i+1)=image_temp(i)
! position of the replica with the temperature i
            image_temp(i)=iii
#endif
! update the number of success of temperature i
            Nsuccess(i)=Nsuccess(i)+1.0d0

#ifdef CPP_MPI
! cautious, this loop is on the replicas in MPI
            i_image=iin
#else
! update the label of the replicas that have been swap
! image i has the temperature kti
            i_image=image_temp(i)
#endif
            if ((abs(kti-kt_all(isize))/k_B).lt.1.0d-8) label(i_image)=-1.0d0
            if ((abs(kti-kt_all(1))/k_B).lt.1.0d-8) label(i_image)=1.0d0

#ifdef CPP_MPI
! cautious, this loop is on the replicas in MPI
            i_image=iii
#else
! image i+1 has the temperature ktn
            i_image=image_temp(i+1)
#endif
            if ((abs(ktn-kt_all(isize))/k_B).lt.1.0d-8) label(i_image)=-1.0d0
            if ((abs(ktn-kt_all(1))/k_B).lt.1.0d-8) label(i_image)=1.0d0

! end if accept
          endif

! reinitialize accept to False
          accept=.False.

       enddo

       end subroutine serial_paratemp

#ifdef CPP_MPI
!--------------------------------------------
! utils for the mpi part

!--------------------------------------------
!--------------------------------------------
       subroutine paratemp_gather_mpi(istart,istop,irank,nRepProc,size_table,MPI_COMM,kt_updated,E_temp,image_temp, &
       &   qeulerp_av,qeulerm_av,Q_sq_sum_av,Qp_sq_sum_av,Qm_sq_sum_av,vortex_av,E_sum_av,E_sq_sum_av,M_sum_av,M_sq_sum_av, &
       &   label_world)
       use m_mpi, only : gather
       implicit none
       integer, intent(in) :: istart,istop,nRepProc,size_table,MPI_COMM,irank
       integer, intent(inout) :: image_temp(:)
       real(kind=8), intent(inout) :: kt_updated(:),E_temp(:)
       real(kind=8), intent(inout) :: qeulerp_av(:),qeulerm_av(:),Q_sq_sum_av(:),Qp_sq_sum_av(:),Qm_sq_sum_av(:),E_sum_av(:),E_sq_sum_av(:)
       real(kind=8), intent(inout) :: vortex_av(:,:),M_sum_av(:,:),M_sq_sum_av(:,:)
       real(kind=8), intent(inout) :: label_world(:)
! dummy variables
       integer :: transfer_I_1D(size_table)
       integer :: ierr

       include 'mpif.h'

       transfer_I_1D=0

! part important for the parallel tempering
       kt_updated=gather(kt_updated(1:nRepProc),nRepProc,size_table,0,MPI_COMM)
       E_temp=gather(E_temp(1:nRepProc),nRepProc,size_table,0,MPI_COMM)

       transfer_I_1D(1:nRepProc)=image_temp(1:nRepProc)+irank*nRepProc
       call mpi_gather(transfer_I_1D(1:nRepProc),nRepProc,MPI_INT,image_temp,nRepProc,MPI_INT,0,MPI_COMM,ierr)

! label and current
       label_world=gather(label_world(1:nRepProc),nRepProc,size_table,0,MPI_COMM)

! variables to reorder after words
       qeulerp_av=gather(qeulerp_av(1:nRepProc),nRepProc,size_table,0,MPI_COMM)
       qeulerm_av=gather(qeulerm_av(1:nRepProc),nRepProc,size_table,0,MPI_COMM)
       Q_sq_sum_av=gather(Q_sq_sum_av(1:nRepProc),nRepProc,size_table,0,MPI_COMM)
       Qp_sq_sum_av=gather(Qp_sq_sum_av(1:nRepProc),nRepProc,size_table,0,MPI_COMM)
       Qm_sq_sum_av=gather(Qm_sq_sum_av(1:nRepProc),nRepProc,size_table,0,MPI_COMM)
       E_sum_av=gather(E_sum_av(1:nRepProc),nRepProc,size_table,0,MPI_COMM)
       E_sq_sum_av=gather(E_sq_sum_av(1:nRepProc),nRepProc,size_table,0,MPI_COMM)

! 2D variables
       vortex_av=gather(vortex_av(:,1:nRepProc),3,nRepProc,size_table,0,MPI_COMM)
       M_sum_av=gather(M_sum_av(:,1:nRepProc),3,nRepProc,size_table,0,MPI_COMM)
       M_sq_sum_av=gather(M_sq_sum_av(:,1:nRepProc),3,nRepProc,size_table,0,MPI_COMM)



       end subroutine paratemp_gather_mpi

!--------------------------------------------
!--------------------------------------------
       subroutine paratemp_scatter_mpi(isize,irank,nRepProc,kt_updated,image_temp,MPI_COMM, &
       &   qeulerp_av,qeulerm_av,Q_sq_sum_av,Qp_sq_sum_av,Qm_sq_sum_av,vortex_av,E_sum_av,E_sq_sum_av,M_sum_av,M_sq_sum_av, &
       &   label_world)
       use m_mpi, only : scatter
       implicit none
       integer, intent(in) :: nRepProc,isize,irank,MPI_COMM
       integer, intent(inout) :: image_temp(:)
       real(kind=8), intent(inout) :: kt_updated(:)
       real(kind=8), intent(inout) :: qeulerp_av(:),qeulerm_av(:),Q_sq_sum_av(:),Qp_sq_sum_av(:),Qm_sq_sum_av(:),E_sum_av(:),E_sq_sum_av(:)
       real(kind=8), intent(inout) :: vortex_av(:,:),M_sum_av(:,:),M_sq_sum_av(:,:)
       real(kind=8), intent(inout) :: label_world(:)
! dummy variables
       real(kind=8) :: kt_local(isize*nRepProc)
       real(kind=8) :: qeulerp_av_local(isize*nRepProc),qeulerm_av_local(isize*nRepProc),Q_sq_sum_av_local(isize*nRepProc), &
        &   Qp_sq_sum_av_local(isize*nRepProc),Qm_sq_sum_av_local(isize*nRepProc),E_sum_av_local(isize*nRepProc),E_sq_sum_av_local(isize*nRepProc)
       real(kind=8) :: vortex_av_local(3,isize*nRepProc),M_sum_av_local(3,isize*nRepProc),M_sq_sum_av_local(3,isize*nRepProc)
       real(kind=8) :: label_world_local(isize*nRepProc)
       integer :: image_temp_local(isize*nRepProc)
       integer :: i_image,i,j,k,ierr

       include 'mpif.h'

       image_temp_local=0
       kt_local=0.0d0
       qeulerp_av_local=0.0d0
       qeulerm_av_local=0.0d0
       Q_sq_sum_av_local=0.0d0
       Qp_sq_sum_av_local=0.0d0
       Qm_sq_sum_av_local=0.0d0
       E_sum_av_local=0.0d0
       E_sq_sum_av_local=0.0d0
       vortex_av_local=0.0d0
       M_sum_av_local=0.0d0
       M_sq_sum_av_local=0.0d0
       label_world_local=0.0d0

       if (irank.eq.0) then

! the labels are attached to the replicas not the temperatures and should not be rearrange
       label_world_local=label_world

! all quantitites attached to the temperature should be rearranged

          do j=1,isize
              do i=1,nRepProc
                 k=i+(j-1)*nRepProc
                 i_image=image_temp(k)
                 kt_local(k)=kt_updated(i_image)
                 image_temp_local(k)=mod(k-1,nRepProc)+1
! reorder the data to match the temperatures
                 qeulerp_av_local(k)=qeulerp_av(i_image)
                 qeulerm_av_local(k)=qeulerm_av(i_image)
                 Q_sq_sum_av_local(k)=Q_sq_sum_av(i_image)
                 Qp_sq_sum_av_local(k)=Qp_sq_sum_av(i_image)
                 Qm_sq_sum_av_local(k)=Qm_sq_sum_av(i_image)
                 E_sum_av_local(k)=E_sum_av(i_image)
                 E_sq_sum_av_local(k)=E_sq_sum_av(i_image)

              enddo
          enddo

         endif

        image_temp(1:nRepProc)=scatter(image_temp_local,nRepProc,isize,0,MPI_COMM)
        kt_updated(1:nRepProc)=scatter(kt_local,nRepProc,isize,0,MPI_COMM)
        qeulerp_av(1:nRepProc)=scatter(qeulerp_av_local,nRepProc,isize,0,MPI_COMM)
        qeulerm_av(1:nRepProc)=scatter(qeulerm_av_local,nRepProc,isize,0,MPI_COMM)
        Q_sq_sum_av(1:nRepProc)=scatter(Q_sq_sum_av_local,nRepProc,isize,0,MPI_COMM)
        Qm_sq_sum_av(1:nRepProc)=scatter(Qm_sq_sum_av_local,nRepProc,isize,0,MPI_COMM)
        Qp_sq_sum_av(1:nRepProc)=scatter(Qp_sq_sum_av_local,nRepProc,isize,0,MPI_COMM)
        E_sum_av(1:nRepProc)=scatter(E_sum_av_local,nRepProc,isize,0,MPI_COMM)
        E_sq_sum_av(1:nRepProc)=scatter(E_sq_sum_av_local,nRepProc,isize,0,MPI_COMM)
        label_world(1:nRepProc)=scatter(label_world_local,nRepProc,isize,0,MPI_COMM)


         kt_updated(nRepProc+1:)=0.0d0
         image_temp(nRepProc+1:)=0


       end subroutine paratemp_scatter_mpi

#endif

       end module
