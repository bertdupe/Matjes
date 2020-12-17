module m_paratemp
use mpi_distrib_v
!--------------------------------------------
! parallel tempering as described in Katzgraber et al PRE and cond-mat 30 mars 2006 arXiv:cond-mat/0602085v3

!--------------------------------------------
!--------------------------------------------
type paratemp_track
    !variable to keep track of what is where on the different parallel tempering steps
    ! parallel tempering label. can take 3 values
    
    !in space of images
    integer,allocatable :: label(:)    ! notifies if this image has last been at lowest(1) or highest(-1) temperature, or has never reached the extrema(0) (0: no label; +1: up; -1: down)
    !in space of temperatures
    real(8),allocatable :: kT_all(:)        ! contains all temperatures T/k_b ordered from smallest to largest 
    real(8),allocatable :: kT_updated(:)    ! used to save updated temperature set 
    integer,allocatable :: nup(:),ndown(:)  ! number of times the temperature has had the up or down label
    integer,allocatable :: Nsuccess(:)      ! keeps track how often a given temperature has moved up (only used for output)
    integer,allocatable :: image_temp(:)    ! gives position of ordered temperature in full image array ( image_temp(1)=2 -> lattice with 1st temperature is at lattice site 2)
end type

private
public :: calculate_diffusion,paratemp_track,alloc_paratemp_track, reorder_temperature

interface calculate_diffusion
    module procedure calculate_diffusion_paratemp_track
    module procedure calculate_diffusion_serial
end interface calculate_diffusion

contains


subroutine alloc_paratemp_track(this,size_table)
    !allocate paratemp_track to the correct size (should only be called on master)
    type(paratemp_track),intent(inout)      :: this
    integer,intent(in)                      :: size_table
    integer                                 :: i

    allocate(this%image_temp(size_table))
    allocate(this%nup(size_table),this%ndown(size_table),this%label(size_table),this%Nsuccess(size_table),source=0) 
    allocate(this%kt_all(size_table),this%kT_updated(size_table),source=0.0d0)
    ! temperature
     ! at the beginning, the first image has the first temperature, the second image, the second ...
    this%image_temp=[(i,i=1,size_table)]
end subroutine

subroutine reorder_temperature(energy,measure,com,nstart,track)
    !subroutine which swaps neighboring temperatures as necessary for the parallel tempering
    use mpi_util 
    use mpi_basic
    use m_mc_exp_val
    use m_constants, only : k_B
    use m_sort
    real(8),intent(inout)                   :: energy(:)    !energy of each image in space of images
    type(exp_val),intent(inout)             :: measure(:)   !Temperature & thermodyn. vals (swap order to swap images)
    class(mpi_distv),intent(in)             :: com          !mpi-communicator which contains one entry for each temperature
    integer,intent(in)                      :: nstart       !outer loop index to alternate between comparing (1<>2,3<>4,..) with (2<>3,4<>5,...)
    type(paratemp_track),intent(inout)      :: track        !various parameters to keep track of parallel tempering evolution
!internals
    integer :: NT   !number of temperatures
    real(8) :: kt_loc(2)        !energies of the compared pairs of lattices
    real(8) :: E_loc(2)         !temperatures of the compared pairs of lattices
    integer :: i_loc(2)         !position in lattices-array of compared pair
    integer         :: i,i_start,i_temp
    type(exp_val)   :: tmp_measure

    !gather all necessary parameters (energies, temperatures(in measure), and expectation values (measure)
    Call gatherv(energy ,com) 
    Call gatherv(measure,com)

    !Swap temperatures if necessary
    if(com%ismas)then
        NT=size(measure)
        !increment nup, ndown for diffusion 
        do i=1,NT
            if (track%label(track%image_temp(i))== 1) track%nup(i)  =track%nup(i)  +1
            if (track%label(track%image_temp(i))==-1) track%ndown(i)=track%ndown(i)+1
        enddo
        
        i_start=mod(nstart,2)+1! in case there are an odd number of replicas, change the start from 1 to 2
        do i_temp=i_start,NT-1,2
            !check if lattices with temperature i_temp and i_temp+1 should be exchanges
            kt_loc=track%kt_all(i_temp:i_temp+1)
            i_loc =track%image_temp(i_temp:i_temp+1)
            E_loc =energy(i_loc)
            if (accept(kt_loc,E_loc)) then
                !Swap the measurements and temperatures
                tmp_measure=measure(i_loc(2))
                measure(i_loc(2))=measure(i_loc(1))
                measure(i_loc(1))=tmp_measure
                track%image_temp(i_temp)=i_loc(2)
                track%image_temp(i_temp+1)=i_loc(1)

                track%Nsuccess(i_temp)=track%Nsuccess(i_temp)+1
            endif
        enddo
        !update label
        track%label(track%image_temp( 1))= 1
        track%label(track%image_temp(Nt))=-1
    endif
    !distribute the rearanged expectation values again
    Call scatterv(measure,com)
end subroutine

function accept(kt,E)
    use m_get_random, only: get_rand_classic
    real(8),intent(in)  ::  kt(2)
    real(8),intent(in)  ::  E(2)
    logical             ::  accept
    real(8)     :: delta
    real(8)     :: choice

    delta=(1.0d0/kt(2)-1.0d0/kt(1))*(E(2)-E(1))
    accept=delta>0.0d0
    if(.not.accept)then
        delta=max(delta,-200.0d0) !prevent underflow in exp
        choice=get_rand_classic()
        accept=exp(delta).gt.Choice
    endif
end function


subroutine calculate_diffusion_paratemp_track(v,relaxation_steps,i_optTset)
    !this simply unwraps the input of the paratemp_track parameter to the old subroutine I did not touch
    type(paratemp_track),intent(inout)      :: v
    integer, intent(in) :: relaxation_steps
    logical, intent(in) :: i_optTset

    call calculate_diffusion(v%kT_all,v%kt_updated,v%nup,v%ndown,v%Nsuccess,relaxation_steps,size(v%kt_all),i_optTset)
end subroutine

subroutine calculate_diffusion_serial(kT_all,kt_updated,nup,ndown,Nsuccess,n_swapT,isize,i_optTset)
       use m_fit
       use m_constants, only : k_B
       implicit none
       integer, intent(in) :: n_swapT,isize
       real(kind=8), intent(in) :: kT_all(:)
       integer, intent(in) :: Nsuccess(:),nup(:),ndown(:)
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
           if (nup(i)+ndown(i)/=0) then
               frac(i)=real(nup(i),8)/real(nup(i)+ndown(i),8)
           else
               write(6,'(a)') 'some up+down fractions are 0'
               write(6,'(a)') 'try to reduce the temperature range to favor overlapping between replicas'
               write(6,*) (frac(j),j=1,i)
               write(6,*) (nup(j),j=1,i)
               write(6,*) (ndown(j),j=1,i)
               if (i_optTset) stop
           endif
       enddo

       fit_frac=fit(frac,isize)

! starting a lots of checks

       if ((frac(1)-1.0d0).gt.1.0d-8) then
           write(6,'(a)') 'the fraction is not 1 at lower temperature'
           write(6,'(a)') 'you can not optimise the temperature set in these conditions'
           write(6,'(a)') 'try to increase N_thousand'
           if (i_optTset) stop
       endif

       if (frac(isize).gt.1.0d-8) then
           write(6,'(a)') 'the fraction is not 0 at higher temperature'
           write(6,'(a)') 'you can not optimise the temperature set in these conditions'
           write(6,'(a)') 'try to increase N_thousand'
!#ifdef CPP_MPI
!           if (i_optTset) call mpi_abort(MPI_COMM,errcode,ierr)
!#else
           if (i_optTset) stop
!#endif
       endif

! if the temperature is not optimized, write the fraction and go out
       if (.not.i_optTset) then
           do i=1,isize-1
               write(21,'(3(E20.12E3,2x))') kT_all(i)/k_B,frac(i),(real(Nsuccess(i),8)/n_swapT)
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
!#ifdef CPP_MPI
!                    call mpi_abort(MPI_COMM,errcode,ierr)
!#else
                    stop
!#endif
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
!#ifdef CPP_MPI
!          call mpi_abort(MPI_COMM,errcode,ierr)
!#else
          stop
!#endif
       endif

       do i=1,isize-1
            write(21,'(7(E20.12E3,2x))') kt_updated(i)/k_B,kT_all(i)/k_B,frac(i),diffusivity(i),Dfrac(i),real(Nsuccess(i),8)/dble(n_swapT),etaTprim(i)
       enddo

       write(21,'(3(E20.12E3,2x),a,/)') kt_updated(isize)/k_B,kT_all(isize)/k_B,frac(isize),'   #  ""  ""  ""  "" '

end subroutine calculate_diffusion_serial
end module
