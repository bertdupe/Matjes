module m_hamiltonian_collection
use m_fft_H_base,only: fft_H
use m_H_public, only: t_H, get_Htype_N
use m_fft_H_public
use m_work_ham_single, only: work_ham_single, work_mode

use m_H_type,only : len_desc
use m_derived_types, only: lattice
use mpi_basic
use mpi_util
use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit
private
public  ::  hamiltonian

type    ::  hamiltonian
    logical                         :: is_set=.false.
    class(t_H),allocatable          :: H(:)
    class(fft_H),allocatable        :: H_fft(:)
    real(8),allocatable             :: H_fft_tmparr(:,:)   !temporary array for effective field from fft
    type(work_mode)                 :: work
    
    logical                         :: is_para(2)=.false. !signifies if any of the internal mpi-parallelizations have been initialized
    type(mpi_type)                  :: com_global
    type(mpi_distv)                 :: com_outer
    type(mpi_type)                  :: com_inner
    integer                         :: NH_total=0  !global size of H
    integer                         :: NH_local=0  !local size of H (relevant with parallelization)
    integer                         :: NHF_total=0 !global size of H_fft
    integer                         :: NHF_local=0 !local size of H_fft (relevant with parallelization, which is not there yet so ==NHF_total)
    integer                         :: N_total=0   !total number of Hamiltonian entries NH_total+NHF_total

    character(len=len_desc),allocatable ::  desc_master(:)
contains
    !creation routines
    procedure   ::  init_H_mv
    procedure   ::  init_H_cp

    !energy getting routines
    procedure   :: energy
    procedure   :: energy_resolved
    procedure   :: energy_single
    procedure   :: energy_distrib
    
    !derivative getting routines
    procedure   :: get_eff_field
    procedure   :: get_eff_field_single

    !parallelization routines
    procedure   :: is_master
    procedure   :: bcast => bcast_hamil        !allows to bcast the Hamiltonian along one communicator before internal parallelization is done
    procedure   :: distribute   !distributes the different Hamiltonian entries to separate threads

    !utility routines
    procedure   :: get_single_work
!    procedure   :: set_work_mode

    !small access routines
    procedure   :: size_H
    procedure   :: get_desc
end type

contains


subroutine distribute(this,com_in)
    use mpi_distrib_v
    !subroutine which distributes the H-Hamiltonians from the master of comm to have as few as possible H per thread
    !expects the hamiltonian of the comm-master to be initialized correctly
    class(hamiltonian),intent(inout)    :: this
    type(mpi_type),intent(in)           :: com_in
#ifdef CPP_MPI

    !integer         :: div, color
    integer         :: ierr

    integer         :: iH,ithread,i

    type(mpi_distv)             :: com_outer
    type(mpi_type)              :: com_inner

    class(t_H),allocatable      :: H_tmp(:)
    class(fft_H),allocatable    :: H_fft_tmp

    if(com_in%ismas)then
        if(.not.this%is_set) ERROR STOP "Cannot distribute hamiltonian that as not been initialized"
        if(allocated(this%H_fft))then
            write(error_unit,*) "WARNING, MPI-distributions for H_fft works only with FFTW_MPI"
        endif
    endif

    Call bcast(this%N_total,com_in)
    Call bcast(this%NH_total,com_in)
    Call bcast(this%NHF_total,com_in)
    this%NH_local=this%NH_total
    this%NHF_local=this%NHF_total

    if(com_in%Np.eq.1) return
    !decide how to parallelize Hamiltonian
    ! parallelized non-FFT Hamiltonian
    if (this%NH_total.eq.0) then
         Call get_two_level_comm(com_in,this%NHF_total,com_outer,com_inner)
    else
         Call get_two_level_comm(com_in,this%NH_total,com_outer,com_inner)
    endif

    if(com_outer%Np>1.and.com_inner%ismas) this%is_para(1)=.true.
    Call reduce_lor(this%is_para(1),com_in)

    !distribute the Hamiltonian to the masters of the inner parallelization
    ! the Hamiltonian are sent as a whole on each of the proc
    if (this%NH_total.ne.0) then
      if (com_inner%ismas)then
         if(com_outer%ismas)then
!            !save descriptions for later easier io-access
!            allocate(this%desc_master(this%N_total))
!            this%desc_master(1:this%NH_total)=this%H(:)%desc
!            do i=1,this%NHF_total
!                Call this%H_fft(i)%get_desc(this%desc_master(this%NH_total+i))
!            enddo

            !send Hamiltonians normal Hamiltonians 
            do ithread=2,com_outer%Np
                do i=1,com_outer%cnt(ithread)
                    iH=com_outer%displ(ithread)+i
                    Call this%H(iH)%send(ithread-1,i,com_outer%com)
                enddo
            enddo
            !keep own ham entry
            this%NH_local=com_outer%cnt(1)
            Call get_Htype_N(H_tmp,this%NH_local)
            do i=1,this%NH_local
                Call this%H(i)%copy(H_tmp(i))    !moving would be much smarter, but has to be implemented
            enddo
            do i=1,size(this%H)
                Call this%H(i)%destroy()
            enddo
            deallocate(this%H)
            Call move_alloc(H_tmp,this%H)
         else
            this%NH_local=com_outer%cnt(com_outer%id+1)
            Call get_Htype_N(this%H,this%NH_local)
            do i=1,this%NH_local
                Call this%H(i)%recv(0,i,com_outer%com)
            enddo
         endif
         this%is_set=.true.
      endif

      if(com_inner%Np>1)then
        this%is_para(2)=.true.
        Call bcast(this%NH_local,com_inner)
        if(.not.com_inner%ismas) Call get_Htype_N(this%H,this%NH_local)
        do i=1,this%NH_local
            Call this%H(i)%distribute(com_inner)
        enddo
        this%is_set=.true.
      endif
    endif   ! end if (this%NH_total.ge.1)

    ! parallalize the FFT Hamiltonian: distribute the FFT_Ham on comm_inner and parallalize then on comm_outer
    ! only FFT Hamiltonians are present
    if (this%NHF_total.ne.0) then

      ! send the FFT Ham on all proc of the inner comms = contains all the masters of the outer comms
      ! first try, parallelize on the inner comm
      if (com_inner%Np.gt.1) then

         this%is_para(2)=.true.

         if (.not.com_inner%ismas) Call get_fft_H_N(this%H_fft,this%NHF_local)
         do i=1,this%NHF_local
             call this%H_fft(i)%bcast(com_inner)
         enddo

      endif

      this%is_set=.true.

    endif     ! end if (this%NHF_total.ge.1)


! get the global communicator for the Hamiltonians without FFT
    this%com_outer=com_outer
    this%com_inner=com_inner
    this%com_global=com_in

    Call set_work_mode(this)
#else
    continue
#endif
end subroutine

subroutine bcast_hamil(this,comm)
    !bcast assuming Hamiltonian has not been scattered yet 
    class(hamiltonian),intent(inout)    :: this
    type(mpi_type),intent(in)           :: comm
#ifdef CPP_MPI
    integer     ::  i
    integer     ::  N_tmp(5)

    if(comm%ismas)then
        if(any(this%is_para))then
            write(error_unit,'(3/A)') "Cannot broadcast Hamiltonian, since it appearst the Hamiltonian already has been scattered"  !world master only contains a part of the full Hamiltonian
            Error STOP
        endif
    endif

    N_tmp=[this%N_total,this%NH_total,this%NHF_total,this%NH_local,this%NHF_local]
    Call bcast(N_tmp,comm)
    this%N_total =N_tmp(1); this%NH_total=N_tmp(2); this%NHF_total=N_tmp(3); this%NH_local=N_tmp(4); this%NHF_local=N_tmp(5)

    if(this%NH_total>1)then
        if(.not.comm%ismas) Call get_Htype_N(this%H,this%NH_total)
        do i=1,this%NH_total
            Call this%H(i)%bcast(comm)
        enddo
    endif
    if(this%NHF_total>1)then
        if(.not.comm%ismas) Call get_fft_H_N(this%H_fft,this%NHF_total)
        do i=1,this%NHF_total
            Call this%H_fft(i)%bcast(comm)
        enddo
    endif
    Call bcast_alloc(this%H_fft_tmparr,comm) 
    if(.not.comm%ismas) Call set_work_mode(this)
#else
    continue
#endif
end subroutine

function size_H(this) result(N)
    class(Hamiltonian),intent(in)   :: this
    integer ::  N

    N=this%N_total
end function

function get_desc(this,i) result(desc)
    use m_H_type, only: len_desc
    class(Hamiltonian),intent(in)       :: this
    integer,intent(in)                  :: i
    character(len=len_desc)             :: desc

    if(i<1.or.i>this%N_total)then
        write(error_unit,'(3/A)') "Cannot get desciption of Hamiltonian, as the index is not with the bounds of the H-array"
        write(error_unit,'(A,I6)')  "Wanted index: ", i
        write(error_unit,'(A,I6)')  "H-arr size  : ", this%N_total
        ERROR STOP
    endif
!    if(i==this%NH_total.and.this%H_fft%is_set())then
!        desc="dipolar"
!        return
!    endif
    if(this%is_para(1))then
        if(this%com_outer%ismas)then
            desc=this%desc_master(i)
        else
            ERROR STOP "CAN ONLY USE GET_DESC FROM MASTER THREAD"
        endif
    else
        if(i<=this%NH_total)then
            desc=this%H(i)%desc
        else
            Call this%H_fft(i-this%NH_total)%get_desc(desc)
        endif
    endif
end function

subroutine init_H_mv(this,lat,Harr,H_fft)
    !initializes the Hamiltonian by moving the H array (thus destroying Harr)
    class(hamiltonian),intent(inout)                :: this
    type(lattice),intent(in)                        :: lat
    class(t_H),allocatable,intent(inout)            :: Harr(:)
    class(fft_H),allocatable,intent(inout),optional :: H_fft(:)
  
    if(.not.allocated(Harr))then
        write(error_unit,'(3/A)') "Cannot initialize Hamiltonian, since Harr-input is not allocated"
        ERROR STOP
    endif
    Call move_alloc(Harr,this%H)
    this%NH_total=size(this%H)
    this%NH_local=this%NH_total
    if(present(H_fft))then
        if(allocated(H_fft))then
            Call move_alloc(H_fft,this%H_fft)
! Bertrand made a change here
! check the note in init_H_cp - this part may not work
            allocate(this%H_fft_tmparr(3*max(lat%nmag,lat%nph),lat%ncell))   !has to be adjusted as well, probably work on same arrays anyways... 
            this%NHF_total=size(this%H_fft)
            this%NHF_local=this%NHF_total
        endif
    endif
    this%N_total=this%NH_total+this%NHF_total
    Call set_work_mode(this)
    this%is_set=.true.
end subroutine

subroutine init_H_cp(this,lat,Harr,H_fft)
    use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit
    !initializes the Hamiltonian by moving the H array (thus destroying Harr)
    class(hamiltonian),intent(inout)                :: this
    type(lattice),intent(in)                        :: lat
    class(t_H),allocatable,intent(in)               :: Harr(:)
    class(fft_H),allocatable,intent(in),optional    :: H_fft(:)
    integer     ::  i
   
    if(.not.allocated(Harr))then
        write(error_unit,'(3/A)') "Cannot initialize Hamiltonian, since Harr-input is not allocated"
        ERROR STOP
    endif
    allocate(this%H,mold=Harr)
    do i=1,size(Harr)
        Call Harr(i)%copy(this%H(i))
    enddo
    this%NH_total=size(this%H)
    this%NH_local=this%NH_total
    if(present(H_fft))then
        if(allocated(H_fft))then
            allocate(this%H_fft,mold=H_fft)
            do i=1,size(H_fft)
                Call H_fft(i)%copy(this%H_fft(i))
            enddo
! Bertrand made a change here
! this part may not work when nph is different from nmag and are both different from 0
! this is a temporary array but it is used on all phonon and magnetic Hamiltonian which might be a problem
            allocate(this%H_fft_tmparr(3*max(lat%nmag,lat%nph),lat%ncell))   !has to be adjusted as well, probably work on same arrays anyways... 
            this%NHF_total=size(this%H_fft)
            this%NHF_local=this%NHF_total
        endif
    endif
    this%N_total=this%NH_total+this%NHF_total
    Call set_work_mode(this)
    this%is_set=.true.
end subroutine


subroutine energy_distrib(this,lat,order,Edist)
    !gets the energy at the sites, so far very wastefull with the memory
    use mpi_distrib_v
    class(hamiltonian),intent(inout)    :: this
    type(lattice), intent(in)           :: lat
    integer,intent(in)                  :: order
    real(8),allocatable,intent(inout)   :: Edist(:,:)

    integer     :: iH
    integer     :: mult(this%com_outer%Np)
    
    if(allocated(Edist))then
        if(size(Edist,1)/=lat%Ncell*lat%site_per_cell(order).and.size(Edist,2)==this%N_total) deallocate(Edist)
    endif
    if(.not.allocated(Edist)) allocate(Edist(lat%Ncell*lat%site_per_cell(order),this%N_total),source=0.0d0)
    do iH=1,this%NH_local
        Call this%H(iH)%energy_dist(lat,order,this%work,Edist(:,iH))
    enddo
    if(this%is_para(2)) Call reduce_sum(Edist,this%com_inner)
    if(this%is_para(1))then
        mult=0
        mult(this%com_outer%id+1)=size(Edist,1)*this%NH_local
        Call reduce_sum(mult,this%com_outer)
        Call gatherv(Edist,mult,this%com_outer)
    endif
    if(order==1)then
        do iH=1,this%NHF_local
            Call this%H_fft(iH)%get_E_distrib(lat,this%H_fft_tmparr,Edist(:,this%NH_total+iH))
        enddo
    endif

end subroutine

subroutine energy_resolved(this,lat,E)
    use mpi_distrib_v, only: gatherv
    !get contribution-resolved energies
    !only correct on outer master thread
    class(hamiltonian),intent(inout):: this
    type (lattice),intent(in)       :: lat
    real(8),intent(out)             :: E(this%N_total)

    integer     ::  iH

    E=0.0d0
    do iH=1,this%NH_local
        Call this%H(iH)%eval_all(E(iH),lat,this%work)
    enddo

    do iH=1,this%NHF_total
        Call this%H_fft(ih)%get_E(lat,this%H_fft_tmparr,E(this%NH_total+iH))
    enddo

    if(this%is_para(2)) Call reduce_sum(E,this%com_inner)
    if(this%is_para(1).and.this%com_inner%ismas) Call gatherv(E,this%com_outer)

!    stop 'toto'
end subroutine

function energy(this,lat)result(E)
    !returns the total energy of the Hamiltonian array
    class(hamiltonian),intent(inout):: this
    type(lattice),intent(in)        :: lat
    real(8)                         :: E

    real(8)     ::  tmp_E(this%N_total)

    Call this%energy_resolved(lat,tmp_E)
    E=sum(tmp_E)
end function

subroutine energy_single(this,i_m,order,lat,work,E)
    use m_derived_types, only: number_different_order_parameters
    use mpi_distrib_v
    !returns the total energy caused by a single entry !needs some updating 
    class(hamiltonian),intent(inout)    :: this
    integer,intent(in)                  :: i_m
    type (lattice),intent(in)           :: lat
    integer,intent(in)                  :: order
    type(work_ham_single),intent(inout) ::  work
    real(8),intent(out)                 :: E

    real(8)     ::  tmp_E(this%N_total)
    integer     ::  iH
#ifdef CPP_DEBUG 
    if(any(this%is_para)) ERROR STOP "IMPLEMENT" !I is rather unlikely that a parallelization here can make this faster
#endif
    E=0.0d0
    do iH=1,this%NH_local
        Call this%H(iH)%eval_single(order)%calc(this%H(iH),tmp_E(iH),i_m,lat,work)
    enddo
    if (this%NH_local.ne.0) tmp_E(:this%NH_local)=tmp_E(:this%NH_local)*real(this%H(:)%mult_M_single,8)

    ! take care of the dipole dipole Hamiltonian if necessary
    do iH=1,this%NHF_total
       Call this%H_fft(ih)%get_E_single(lat,i_m,tmp_E(this%NH_total+iH))
    enddo

    E=sum(tmp_E)

end subroutine 

subroutine get_eff_field(this,lat,B,Ham_type)
    !calculates the effective internal magnetic field acting on the magnetization for the dynamics
    class(hamiltonian),intent(inout)    :: this
    type (lattice),intent(in)           :: lat    !lattice containing current order-parameters 
    real(8),intent(out)                 :: B(:)
    integer,intent(in)                  :: Ham_type   !integer that decides with respect to which mode the Hamiltonians derivative shall be obtained [1,number_different_order_parameters]

    integer     :: iH

    B=0.0d0
    do iH=1,this%NH_local
        Call this%H(iH)%deriv(Ham_type)%get(this%H(iH),lat,B,this%work)
    enddo

    if(any(this%is_para)) Call reduce_sum(B,this%com_global)
    B=-B    !field is negative derivative

    do iH=1,this%NHF_local
        Call this%H_fft(iH)%get_H(lat,this%H_fft_tmparr)
        B=B-2.0d0*reshape(this%H_fft_tmparr,shape(B))
    enddo
end subroutine

subroutine get_eff_field_single(this,lat,site,B,work,Ham_type,tmp)
    !calculates the effective internal field acting on a single site
    class(hamiltonian),intent(inout)    :: this
    type (lattice),intent(in)           :: lat    !lattice containing current order-parameters 
    integer,intent(in)                  :: site
    real(8),intent(out)                 :: B(:)
    type(work_ham_single),intent(inout) :: work     !work arrays
    integer,intent(in)                  :: Ham_type   !integer that decides with respect to which mode the Hamiltonians derivative shall be obtained [1,number_different_order_parameters]
    real(8),intent(inout)               :: tmp(size(B))

    integer     :: iH
    B=0.0d0
    do iH=1,this%NH_local
        Call this%H(iH)%deriv(Ham_type)%get_single(this%H(iH),lat,site,work,B,tmp)
    enddo

    if(any(this%is_para)) Call reduce_sum(B,this%com_global)
    B=-B    !field is negative derivative

    do iH=1,this%NHF_local
        Call this%H_fft(iH)%get_H_single(lat,site,tmp)
        B=B-2.0d0*tmp
    enddo
end subroutine


function is_master(this)result(master)
    class(hamiltonian),intent(in)       :: this
    logical                             :: master

    if(this%is_para(1))then
        master=this%com_outer%ismas
    else
        if(this%is_para(2))then
            master=this%com_inner%ismas
        else
            master=.true.
        endif
    endif
end function

subroutine get_single_work(this,order,work)
    class(hamiltonian),intent(inout)    :: this
    integer,intent(in)                  :: order
    type(work_ham_single),intent(inout) :: work

    type(work_ham_single)               :: work_arr(this%NH_local)
    integer     :: iH

    do iH=1,this%NH_local
        Call this%H(iH)%set_work_single(work_arr(iH),order)
    enddo
    Call work%set_max(work_arr)
    do iH=1,this%NH_local
        Call work_arr(iH)%destroy()
    enddo
end subroutine

subroutine set_work_mode(this)
    class(hamiltonian),intent(inout)    :: this

    type(work_mode)     :: work_arr(this%NH_local)
    integer     :: iH

    do iH=1,this%NH_local
        Call this%H(iH)%set_work_mode(work_arr(iH))
    enddo
    if (this%NH_local.ne.0) Call this%work%set_max(work_arr)
    do iH=1,this%NH_local
        Call work_arr(iH)%destroy()
    enddo
end subroutine

subroutine allocate_H_fft_tmparr(H_fft_tmparr,H_fft,lat)
    ! Helper-routine to allocate the h_fft_tmparr with the correct maximal size 
    use m_type_lattice, only: number_different_order_parameters
    real(8),intent(inout),allocatable   :: H_fft_tmparr(:,:)   !temporary array for effective field from fft
    type(lattice),intent(in)            :: lat
    class(fft_H),intent(in)             :: H_fft(:)

    logical                     :: order_mask(size(lat%dim_modes))
    integer                     :: size1,i

    do i=1,size(H_fft)
        order_mask(H_fft(i)%order) = .true.
    enddo

    size1 = maxval(lat%dim_modes,order_mask)
    allocate(H_fft_tmparr(size1,lat%ncell))
end subroutine

end module
