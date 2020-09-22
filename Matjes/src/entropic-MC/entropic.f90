!
! Entropic sampling Monte Carlo
!
!
subroutine entropic(my_lattice,my_motif,io_simu,ext_param)
use mtprng
use m_choose_spin
use m_constants, only : pi,k_B
use m_relaxtyp
use m_sampling
use m_derived_types
use m_lattice, only : my_order_parameters
use m_local_energy
use m_total_energy
use m_operator_pointer_utils
implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variable part
type(lattice), intent(inout) :: my_lattice
type(cell), intent(in) :: my_motif
type(io_parameter), intent(in) :: io_simu
type(simulation_parameters), intent(in) :: ext_param
!------------------------------------------
! internal variable
type(mtprng_state) :: state
integer :: Ilat(4),n_steps,i_step,n_sweep,i,i_loop,n_loop
real(kind=8) :: mu_s,cone,choice
logical :: accept,i_exit
! coordinate of the spin in internal coordinate
integer :: iomp
real(kind=8) :: S_new(3),S_old(3)
! energy variable
real(kind=8) :: E_min,E_max,E_delta,E_test
real(kind=8) :: E_new,E_old,E_total,kTfin,kTini,Et
integer,allocatable :: count_E(:)
! histogram variables
integer :: n_histo,n_old,n_new,min_histo,i_histo
real(kind=8),allocatable :: lng(:)
real(kind=8) :: lnf,lng_new,lng_old,delta_lng,flatness,mean_histo
! thermodynamical quantitites
real(kind=8),allocatable :: energy(:),Z(:),temperatures(:),free_E(:),entropy(:),energy_sq(:),heat_cap(:)
real(kind=8) :: kT,DProb,DE
integer :: i_T,N_cell,n_Tsteps,dim_mode
logical :: sphere,equi,ising,i_magnetic,underrel,overrel

type(vec_point),pointer :: all_mode(:)
type(vec_point),allocatable,dimension(:) :: mode_magnetic

!call mtprng_init(37, state)

sphere=.true.
equi=.false.
ising=.false.
underrel=.false.
overrel=.false.
n_Tsteps=10
! initialize variables
allocate(energy(n_Tsteps),Z(n_Tsteps),temperatures(n_Tsteps),free_E(n_Tsteps),entropy(n_Tsteps),energy_sq(n_Tsteps),heat_cap(n_Tsteps))
N_cell=product(my_lattice%dim_lat)
mu_s=0.0d0
E_min=0.0d0
E_max=0.0d0
E_new=0.0d0
E_old=0.0d0
E_delta=0.0d0
n_histo=90
allocate(lng(n_histo))
energy=0.0d0
heat_cap=0.0d0
energy_sq=0.0d0
entropy=0.0d0
free_E=0.0d0
Z=0.0d0
temperatures=0.0d0
lng=0.0d0
cone=pi(1.0d0)
S_new=0.0d0
S_old=0.0d0
n_steps=100000000
n_loop=20
lnf=0.50d0
lng_new=0.0d0
lng_old=0.0d0
delta_lng=0.0d0
n_old=0
n_new=0
accept=.False.
allocate(count_E(n_histo))
count_E=0
E_test=0.0d0
flatness=0.9d0
i_exit=.False.
kTfin=1.0d0
kTini=1.0d0
dim_mode=my_lattice%dim_mode

all_mode(1:N_cell)=>my_lattice%ordpar%all_l_modes

! magnetization
do i=1,size(my_order_parameters)
  if ('magnetic'.eq.trim(my_order_parameters(i)%name)) then
   allocate(mode_magnetic(N_cell))
   call dissociate(mode_magnetic,N_cell)
   call associate_pointer(mode_magnetic,all_mode,'magnetic',i_magnetic)
  endif
enddo

call set_E_matrix(my_lattice%dim_mode)

! initializing the variables above
Call sum_energy(E_total,my_lattice)
write(6,'(a,2x,E20.12E3)') 'Initial Total Energy (eV)',E_total/real(N_cell)


write(6,'(a)') "in entropic sampling subroutine"
! find the energy range

E_min=-8.0d-3*real(N_cell)
E_max=8.0000001d-3*real(N_cell)
#ifdef CPP_DEBUG
if ((E_total-E_min).lt.0.0d0) then
    write(6,'(a,f12.7)') 'the minimum of energy ', E_min
    write(6,'(a,f12.7)') 'is smaller than the energy ', E_total
    stop
endif
#endif

E_delta=(E_max-E_min)/real(n_histo)
write(6,'(a,E20.10E3)') 'The size of the energy bin is ', E_delta

! open the results file

OPEN(7,FILE='histogram.dat',action='write',status='unknown',position='append',form='formatted')

do i_loop=1,n_loop
   do i_step=1,n_steps
! pic a spin
      call choose_spin(iomp,N_cell)

      S_old=mode_magnetic(iomp)%w
!---------------------------------------
! here are the different sampling
! first the sphere sampling
!---------------------------------------
      if (sphere) then
         S_new=sphereft(mode_magnetic(iomp)%w,cone)
!---------------------------------------
! Algorithm square sampling
!---------------------------------------
      elseif (equi) then
         S_new=equirep(mode_magnetic(iomp)%w,cone)
      endif
!---------------------------------
! different relaxation process
!---------------------------------
      if (underrel) then
         S_new=underrelax(iomp,my_lattice)
      elseif (overrel) then
         S_new=overrelax(iomp,my_lattice)
      endif

      if (ising) S_new=-mode_magnetic(iomp)%w

! calculate the energy before and after

      call local_energy(E_old,iomp,my_lattice)

      mode_magnetic(iomp)%w=S_new

      call local_energy(E_new,iomp,my_lattice)

      mode_magnetic(iomp)%w=S_old
      E_test=E_total+E_new-E_old

! position of the energy difference in energy space
   ! check if the energy is smaller that the minimum of energy
!             if ((E_test.lt.E_min).or.(E_test.gt.E_max)) cycle

   ! locate the position of E
      n_old=INT((E_total-E_min)/E_delta)+1
      n_new=INT((E_total+E_new-E_old-E_min)/E_delta)+1

   ! density of states
      lng_old=lng(n_old)
      lng_new=lng(n_new)
      delta_lng=lng_old-lng_new

   ! Do I accept the move or not???
      if (delta_lng.ge.1.0d2) then
         accept=.True.
      else

#ifdef CPP_MRG
         choice=mtprng_rand_real1(state)
#else
         CALL RANDOM_NUMBER(choice)
#endif

         if (choice.lt.1.0d-9) then
            accept=.True.
         else
            if (log(choice).lt.delta_lng) accept=.True.
         endif

      endif

      if (accept) then
         mode_magnetic(iomp)%w=S_new
         lng(n_new)=lng(n_new)+lnf
         count_E(n_new)=count_E(n_new)+1
         E_total=E_total+E_new-E_old
         accept=.False.
!                length=length-int(S_old(3))+int(S_new(3))
      else
         count_E(n_old)=count_E(n_old)+1
         lng(n_old)=lng(n_old)+lnf
      endif

!             vector((N_cell+length)/2+1)=vector((N_cell+length)/2+1)+1
! check the flatness of the histogramm
! The minimum of energy and the max must be taken out
      n_sweep=sum(count_E)
      mean_histo=real(n_sweep)/real(N_cell)
      min_histo=minval(count_E)

#ifdef CPP_DEBUG
      if (mod(i_step,100).eq.0) then
         write(*,*) 'Check'
         write(*,*) lng
         write(*,*) count_E
         write(*,*) min_histo,flatness*mean_histo,N_cell,n_sweep,mean_histo
      endif
#endif

      if ((mod(i_step,1000).eq.0).and.(real(min_histo).gt.flatness*mean_histo)) then
         write(6,'(a,2x,I10,2x,a)') 'histogramm flat enough after',i_step,'steps'
         write(6,'(a/)') 'reset the histogramm and the parameter f'
         write(6,'(/a,I10/)') 'end of loop number', i_loop
         exit
      endif

      if (i_step.eq.n_steps) then
         write(6,'(a)') 'The simulation did not converge - Histogram not flat'
         i_exit=.true.
         exit
      endif
   enddo   ! end over the loop of the energy histogram

   write(6,'(/a,I10/)') "writing the results for loop", i_loop

   do i=1,n_histo
        write(7,'(3(E20.10E3,2x),I10)') (E_min+real(i-1)*E_delta)/real(N_cell),lng(i),real(count_E(i))/real(N_cell),count_E(i)
   enddo
   write(7,'(a)') ''

   if (i_exit) stop 'stoping the simulation'

! update of f and reset of the the energy histogram
   lnf=lnf/2.0d0
   count_E=0

enddo

close(7)

#ifdef CPP_DEBUG
      OPEN(666,FILE='test.dat',action='write',status='unknown',form='formatted')
      do i=1,N_cell+1
        write(666,'(I10)') vector(i)
      enddo
      close(666)
#endif

write(6,'(/a)') "-------------------------"
write(6,'(a/)') "writing the final results"

! renormalizing the density g(E)
!      lng=lng-log(real(N_cell))
delta_lng=minval(lng)

OPEN(7,FILE='lngE.dat',action='write',status='unknown',form='formatted')
do i=1,n_histo
   write(7,'(2(E20.10E3,2x))') (E_min+real(i-1)*E_delta)/real(N_cell),lng(i)-delta_lng
enddo
close(7)

! calculate thermodynamical quantitites

do i_T=1,n_Tsteps
   kT=(kTfin-kTini)/real(n_Tsteps-1)*real(i_T-1)+kTini
   temperatures(i_T)=kT/k_B
enddo

do i_T=1,n_Tsteps
   kT=temperatures(i_T)*k_B
   do i_histo=1,n_histo

      DE = (real(i_histo-1)*E_delta+E_min)
      DProb=lng(i_histo)-delta_lng-(DE+E_max)/kt

      Z(i_T)=Z(i_T)+exp(DProb)
      energy(i_T)=energy(i_T)+DE*exp(DProb)
      energy_sq(i_T)=energy_sq(i_T)+(DE)**2*exp(DProb)

   enddo
! renormalize the energy
   energy(i_T)=energy(i_T)/Z(i_T)
   energy_sq(i_T)=energy_sq(i_T)/Z(i_T)
   free_E(i_T)=-kT*log(Z(i_T))
   entropy(i_T)=k_B*(log(Z(i_T))+(energy(i_T)+E_max)/kT)
   heat_cap(i_T)=(energy_sq(i_T)-energy(i_T)**2)/kT**2*k_B
enddo

OPEN(7,FILE='EM.dat',action='write',status='unknown',form='formatted')
write(7,'(a)') '#1:T   2:Z(T)   3:E(T)   4:F(B,T)   5:S(B,T)   6:C_V(B,T)'
do i_T=1,n_Tsteps
   write(7,'(6(E20.10E3,2x))') temperatures(i_T),Z(i_T),energy(i_T),free_E(i_T),entropy(i_T),heat_cap(i_T)
enddo
close(7)

nullify(all_mode)
write(6,'(a)') "done with entropic"

end subroutine entropic
