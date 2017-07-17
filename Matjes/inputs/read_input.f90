module m_parameters
implicit none
!    #  *** general ***
!    #  H_ext: external magnetic field in x y z
      real(kind=8), Dimension(1:3) :: H_ext
!    #  number of temperature steps
      integer :: n_Tsteps
!    #  T_ini: initial temperature
!    #  T_fin: final temperature
      real(kind=8) :: T_ini, T_fin
!    #  i_ising: Ising model (default: Heisenberg-model)
      logical :: i_ising
!#######################################
!    # *** Spinstructure ***
!    #  type_lattice: lattice type. choose of: square, hexagonal, userDefined (here you need an additional input file, the inputs here are ignored)
!    #  type_system: monolayer, multilayer, superlattice
!    #  n_unitCells: number of unit cells along the lattice vectors
!    #  periodic boundary conditions in x,y,z direction
!    #  n_atomsUnitCell: Atoms in unit cell, magnetic and nonmagnetic for DMI
!    #  pos_atom1: position in unitcell and magnetic moment of atom 1. Respectively for further atoms in same unit cell.
      logical :: i_hexagonalLattice, i_squareLattice, i_userDefinedLattice
      logical :: i_monolayer, i_multilayer, i_superlattice
      real(kind=8), Dimension(1:3) :: n_unitCells
      logical :: i_periodic_boundary(3)
      integer :: n_atomsUnitCell
      real(kind=8), Dimension(1:4) :: pos_atom1, pos_atom2, pos_atom3
!#######################################
!    #  *** setup of the starting config ***
!    #  *!If a SpinSTMi.dat file, the inputs here will be ignored!*
!    #  qvector: propagation vector of the spin spiral
!    #  Orientation of magnetization m at position r will be created as:
!    #  m(r) -- Rq*cos(qvector*r)+Iq*sin(qvector*r)
!    #  heavyside: create a domainwall
      real(kind=8), Dimension(1:3) :: qvector, Rq, Iq
      logical :: heavyside
!#######################################
!    #  *** Coefficients for the Hamiltonian ***
!    #  J_i: Exchange interaction of i-th nearest neighbour (all togehther saved in J_ij)
!    #  J_ml1: interaction in multilayer (intra-bilayer exchange interaction)
!    #  J_sl: superlattice: interaction between multilayers (inter-bilayer exchange interaction)
      real(kind=8) :: J_1, J_2, J_3, J_4, J_5, J_6, J_7, J_8, J_9, J_10, J_11, J_12, J_ml1, J_ml2, J_ml3, J_ml4, J_sl
! Dzyaloshinsky-Moriya interaction
      real(kind=8) :: DM(3)
! easy axis
      real(kind=8) :: EA(3)
! Anisotropie Vector
      real(kind=8), Dimension(1:3) :: D_ani
! dipole-dipole interaction
      logical :: i_dipdip
! stoner parameter
      logical :: i_stoner
      real(kind=8) :: Ist
! biquadratic spin coupling constant
      logical :: i_biq
      real(kind=8) :: J_biq

! four spin interaction
      real(kind=8) :: K_4spin
! constants in the sum
      real(kind=8) :: c_Ji,c_DM,c_JB,c_Ki,c_ani

!#############################################
!    #  *** type of simulation ***
!    #  please choose one of the following types for your simulation (default: metropolis):
!    #  metropolis, paratemp, spindynamic, gneb
!    #  i_warnings: optional, print warnings during the run
      logical :: i_metropolis, i_parallelTempering, i_spinDynamic, i_gneb
      logical :: i_warnings
!########################
!    ###  Monte Carlo ###
!    #  n_measurement: measurement steps on which you make the average
!    #  t_autocorrelation: defines the autocorrelation steps between the measurement steps
!    #  t_equilibration is multiplied by the number of spins
!    #  i_writeEquilibration: whether equilibrationfile is written: default: .F.
!    #  n_writeEquilibration: how often, should be .lessEqual. t_equilibration
!    #  type_equilibration: type of equilibration (normal, ??? )
!    #  type_sampling: type of sampling (sphere,???)
!    #  i_optimizationTset: optimization of the temperature set
!    #  n_optimizationSteps: number of optimization steps
!    #  coneAngle: starting angle in which the new spin is choosen
    integer :: n_measurement,t_autocorrelation,t_equilibration,n_writeEquilibration,n_optimizationSteps
    logical :: i_writeEquilibration, i_equiNormal, i_samplingSphere,i_optimizationTset
    real(kind=8) :: coneAngle

!########################
!    ### spin dynamic ###
!    # integration: which type of solver? (1: Euler; 2: SIA; 3: Heun; 4: SIB; 5: SIB with error correction, 6: SIB without temperature)
!    # timestep in fs
!    # duration in units of dt
!    # Efreq
!    # damping without units
!    # resolution for povray script
!    # torque
!    # damping torque
!    # Ipol unit vector of polarisation of spin current
!    # adia(batic) and nonadia(batic) term of SST
!    # stmtorque
!    # step
!    # rampe
!    # Hsweep
!    # Ffield
!    # Efield
!    # STMtemp
      integer :: integration, timestep, duration, Efreq
      real(kind=8) :: damping, torque, damping_torque
      integer :: resolution(2)
      real(kind=8) :: Ipol(3)
      real(kind=8) :: adia, nonadia
      logical :: i_stmtorque, i_step, i_rampe, i_Hsweep
      real(kind=8) :: stmtorque(3),step(2), rampe(3), Hsweep (8)
      logical :: Ffield, Efield, STMtemp

!#######################################
!    #  *** MPI parameters ***
!    #  i_ghost: logical decompose the lattice into for sublattice to pick up the spins in parallel
!    #  n_ghost: number of sublattices for ghost. Please choose carefully.
!    #  mpi_algo_separate: the temperatures were simulated simultaneously
!    #  n_TProc: How many temperatures should be simulated on one Proc?
!    #  *! Make sure: n_Tsteps/n_recProc= number of processors !*
      logical :: i_ghost, mpi_algo_separate
      integer :: n_ghost,n_TProc
!#######################################
!    #  *** What do you want to print additionally after a metropolis simulation? ****
!    #  i_CalTheta: prints the angles of the spins
!    #  i_CalEnergy: prints the different energy contributions spin resolved
!    #  i_printSpinse: print the Spinstructure for povray
!    #  n_printSpinse: frequency of sprinted spinstructures (important for spin dynamics)
!    #  i_topoCharDensity: print the topological charge density
!    #  Cor_log: ???
!    #  qorien: ???
!    #  gra_fft: fourier transform of what ???
!    #  dispersion: how is it plotted
!    #  i_Tsweep: Is is still needed? ???
      logical :: i_CalTheta, i_CalEnergy, i_printSpinse,i_topoCharDensity, i_gra_fft, Cor_log,qorien, i_dispersion, i_Tsweep
      integer :: n_printSpinse, n_dispersion
      integer :: n_gra_fft(2)

end module m_parameters




! Subroutine to read the input file of Matjes
program read_input
 use m_parameters
 implicit none

! what could and should be read
character(len=25) :: type_simulation
character(len=25) :: type_equilibration
character(len=25) :: type_sampling
character(len=25) :: type_lattice
character(len=25) :: type_system

logical :: i_sphere, i_normal


namelist /general/ type_simulation, H_ext,n_Tsteps,T_ini,T_fin, i_ising,i_warnings
namelist /MC/ type_sampling, type_equilibration, n_measurement,t_autocorrelation,t_equilibration,n_writeEquilibration, &
 &n_optimizationSteps, i_writeEquilibration, i_equiNormal, i_samplingSphere,i_optimizationTset, coneAngle
namelist /SD/ integration, timestep, duration, Efreq, damping, torque, damping_torque, resolution, Ipol, adia, nonadia, &
              & i_stmtorque, i_step, i_rampe, i_Hsweep, stmtorque, step, rampe, Hsweep, Ffield, Efield, STMtemp
!namelist /GNEB/
namelist /MPI/ i_ghost, mpi_algo_separate, n_ghost, n_TProc
namelist /LAT/ type_lattice, type_system, n_unitCells, i_periodic_boundary, n_atomsUnitCell, pos_atom1, pos_atom2, &
      &pos_atom3, qvector,Rq,Iq, heavyside

namelist /HAM/ J_1, J_2, J_3, J_4, J_5, J_6, J_7, J_8, J_9, J_10, J_11, J_12, J_ml1, J_ml2, J_ml3, J_ml4, J_sl,&
 & EA, D_ani, DM, i_dipdip, Ist, J_biq, K_4spin, c_Ji,c_DM,c_JB,c_Ki,c_ani
namelist /OPT/ i_CalTheta, i_CalEnergy, i_printSpinse,i_topoCharDensity, i_gra_fft, Cor_log, qorien, i_dispersion, &
       & i_Tsweep, n_printSpinse, n_dispersion, n_gra_fft

! needed for the reading process
integer :: fin,i,k
integer, parameter :: io=1
character(len=100) :: str

! all the parameters are initialized now
! general
type_simulation = 'metropolis'
H_ext = (/0.0, 0.0, 0.0/)
n_Tsteps = 1
T_ini = 1.0d0
T_fin = 10.0d0
i_ising = .false.
i_warnings = .true.
! MC
type_sampling = 'sphere'
type_equilibration = 'normal'
n_measurement = 1000
t_autocorrelation=1
t_equilibration = 1000
n_writeEquilibration = 10
n_optimizationSteps= 5
i_writeEquilibration = .true.
i_equiNormal = .true.
i_samplingSphere = .true.
i_optimizationTset= .false.
coneAngle = 1.0d0
! SD
integration = 1
timestep = 1
duration = 1
Efreq = 1
damping = 0.0d0
torque = 0.0d0
damping_torque = 0.0d0
resolution = (/100,100/)
Ipol= (/ 0.0, 0.0, 1.0/)
adia= 0.0
nonadia= 0.0
i_stmtorque = .false.
stmtorque= (/1.0, 1.0, 0.0/)
i_step= .false.
step= (/1000, 1000/)
i_rampe= .false.
rampe= (/-10.0, 10.0, 1000.0/)
i_Hsweep= .false.
Hsweep= (/0.0, 0.0, -10.0, 1000.0, 1500.0, 0.0, 0.0, 0.0 /)
Ffield= .false.
Efield= .false.
STMtemp= .false.
! GNEB

! MPI
i_ghost = .false.
mpi_algo_separate = .false.
n_ghost = 1
n_TProc = 1

! Lattice
type_lattice = 'hexagonal'
type_system = 'monolayer'
n_unitCells = (/1, 1, 1 /)
i_periodic_boundary = (/.true., .true., .true. /)
n_atomsUnitCell = 1
pos_atom1 = (/0.0, 0.0, 0.0, 0.0/)
pos_atom2 = (/0.0, 0.0, 0.0, 0.0/)
pos_atom3 = (/0.0, 0.0, 0.0, 0.0/)
qvector = (/0.0,0.0,0.0/)
Iq = (/0.0,0.0,0.0/)
Rq=(/0.0,0.0,1.0/)
heavyside = .false.
! Hamiltonian
J_1= 0.5d-3
J_2= 0.0d-3
J_3= 0.0d-3
J_4= 0.0d-3
J_5= 0.0d-3
J_6= 0.3d-3
J_7= 0.0d-3
J_8= 0.0d-3
J_9= 0.0d-3
J_10= 0.0d-3
J_11= 0.0d-3
J_12= 0.0d-3
J_ml1= 0.0d0
J_ml2= 0.0d0
J_ml3= 0.0d0
J_ml4 = 0.0d0
J_sl= 0.0d0
J_biq= 0.04d-0
EA= (/0.0d0, 0.0d0, 1.0d0/)
D_ani=(/ 0.000d-0, 0.0d0, 0.0d0/)
K_4spin= 0.0d-3
DM= (/ 0.0d-3, 0.0d-3, 0.0d-3 /)
i_dipdip= .false.
Ist = 0.0
c_Ji= -1.0d0
c_DM= -1.0d0
c_JB= -1.0d0
c_Ki= -1.0d0
c_ani= 1.0d0

! Options
i_CalTheta= .false.
i_CalEnergy= .false.
i_printSpinse= .true.
n_printSpinse= 1
i_topoCharDensity= .false.
Cor_log= .false.
qorien= .false.
i_gra_fft= .false.
n_gra_fft= (/100, 100 /)
i_dispersion= .false.
n_dispersion= 50
i_Tsweep= .false.
k=0

open (io,file='input',form='formatted', status= 'old', action='read')

rewind(io)
 do
  k=k+1
  read(io,'(a)',iostat=fin) str
  if (fin /= 0) exit           ! exit at the last row??
  str=trim(adjustl(str))
  if (len_trim(str)==0) cycle  ! ignore empty lines
  if (str(1:1) == '#') cycle   ! ignore Comments

  read(str,nml=general)
  read(str,nml=LAT)
  if (trim(type_lattice)=='hexagonal') then
   i_squareLattice=.false.
   i_hexagonalLattice=.true.
   i_userDefinedLattice=.false.
  elseif (trim(type_lattice)=='square') then
   i_squareLattice=.true.
   i_hexagonalLattice=.false.
   i_userDefinedLattice=.false.
  elseif ((trim(type_lattice)=='userDefined').or.(trim(type_lattice)=='userdefined')) then
   i_squareLattice=.false.
   i_hexagonalLattice=.false.
   i_userDefinedLattice=.true.
!  else
!   write(6,*) 'Please check type of lattice'
  endif
  if (trim(type_system)=='monolayer') then
   i_monolayer=.true.
   i_multilayer=.false.
   i_superlattice=.false.
  elseif (trim(type_system)=='multilayer') then
   i_monolayer=.false.
   i_multilayer=.true.
   i_superlattice=.false.
  elseif (trim(type_system)=='superlattice') then
   i_monolayer=.false.
   i_multilayer=.false.
   i_superlattice=.true.
!  else
!   write(6,*) 'Please check type of system'
  endif
  read(str,nml=HAM)
  read(str,nml=MPI)
 enddo

 if (trim(type_simulation)=='metropolis') then
   read(str,nml=MC)
   i_metropolis=.true.
   if (trim(type_sampling)=='sphere') i_samplingSphere=.true.
   if (trim(type_equilibration)=='normal') i_equiNormal=.true.

 elseif  (trim(type_simulation)=='paratemp') then
   read(str,nml=MC)
   i_parallelTempering=.true.
   if (trim(type_sampling)=='sphere') i_sphere=.true.
   if (trim(type_equilibration)=='normal') i_normal=.true.

 elseif  (trim(type_simulation)=='spindynamic') then
   read(str,nml=SD)
   i_spinDynamic=.true.
! elseif  (trim(type_simulation)=='GNEB') then
!   read(str,nml=GNEB)
!   i_GNEB = .true.
! else
!  write(6,*) 'Please check the type of Simulation!'
!  call abort
 endif

close(io)


write(*,nml=general)
write(*,HAm)
write(*,MC)
write(*,SD)
write(*,LAT)

end program read_input
