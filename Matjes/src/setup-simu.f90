      subroutine setup_simu(my_simu,io_simu,state)
      use m_lattice
      use m_rw_lattice
      use m_voisins
      use m_setup_DM
      use m_parameters
      use m_sym_utils
      use m_derived_types
      use m_table_dist
      use m_mapping
      use m_indexation
      use m_arrange_neigh
      use m_topoplot
      use mtprng
#ifdef CPP_MPI
      use m_make_box
      use m_split_work
      use m_mpi_prop, only : isize,irank,irank_working,N,start
#endif
#ifndef CPP_BRUTDIP
      use m_setup_dipole
#endif
      implicit none
! this subroutine is used only to setup the simulation box
! it reads first the parameters of the simulation i.e. inp file
! then it reads the lattice
! this order aims at not taking care of too many neigbours if the Jij are set to 0
      type(io_parameter), intent(out) :: io_simu
      type(type_simu), intent(out) :: my_simu
      type(mtprng_state), intent(inout) :: state
! variable of the system
      real (kind=8), allocatable :: tabledist(:,:)
      real (kind=8), allocatable :: map_vort(:,:,:,:),map_toto(:,:,:)
      integer :: N_Nneigh,phase,tot_N_Nneigh
!nb of neighbors in the z direction
      integer :: Nei_z
!nb of neighbors in the superlattice
      integer :: Nei_il
!     indeNN defined in maindata
! dummy variable
      integer :: i,j,k,l
!checking various files
      logical :: exists,topoonly,i_exi,i_usestruct
! check the allocation of memory
      integer :: alloc_check

! innitialisation of dummy
      phase=1
      Nei_z=0
      i_usestruct=.False.
      alloc_check=0

#ifdef CPP_MPI
if (irank.eq.0) write(6,'(/,a)') 'reading the setup routine'
#else
write(6,'(/,a)') 'reading the setup routine'
#endif
! read the important inputs
#ifdef CPP_MPI
if (irank.eq.0) write(6,'(a)') 'reading the lattice.in file'
#else
write(6,'(a)') 'reading lattice.in file'
#endif
      call rw_lattice()

#ifdef CPP_MPI
if (irank.eq.0) write(6,'(a)') 'reading the inp file'
#else
write(6,'(a,/)') 'reading the inp file'
#endif
      call inp_rw(my_simu,io_simu,N_Nneigh,phase,Nei_z,Nei_il)

#ifdef CPP_MPI
if (irank.eq.0) write(6,'(a)') 'allocating the spin, table of neighbors and the index...'
#else
write(6,'(a)') 'allocating the spin, table of neighbors and the index...'
#endif

      call create_lattice(dim_lat,motif)

! of ghost is there, then make the box
#ifdef CPP_MPI
      if (i_ghost) then

       if (irank.eq.0) then
        write(6,'(a)') "domain decomposition parallelisation selected"
        write(6,'(I5,a)') n_ghost," domains split"
       endif

       if (size(world).eq.1) then
        call box_1D(dim_lat(1),n_ghost,Periodic_log,N_Nneigh/2)
        N_site=(Xstop-Xstart+1)*count(motif%i_m)
       elseif (size(world).eq.2) then
        call box_2D(dim_lat(1),dim_lat(2),n_ghost,Periodic_log,N_Nneigh/2)
        N_site=(Xstop-Xstart+1)*count(motif%i_m)*(Ystop-Ystart+1)
       else
        call box_3D(dim_lat(1),dim_lat(2),dim_lat(3),n_ghost,Periodic_log)
        N_site=(Xstop-Xstart+1)*count(motif%i_m)*(Ystop-Ystart+1)*(Zstop-Zstart+1)
       endif

       if (irank.eq.0) write(6,'(a)') "--------------------------"
       if (irank/n_ghost.eq.0) write(6,'(a,I5,3(a,2I8))') "on process ",irank," domain: X",Xstart,Xstop,", Y: ",Ystart,Ystop, ", Z: ",Zstart,Zstop

      else
      Xstart=1
      Xstop=dim_lat(1)
      Ystart=1
      Ystop=dim_lat(2)
      Zstart=1
      Zstop=dim_lat(3)
      N=(/Xstop,Ystop,Zstop/)
      start=0
      irank_working=irank
      endif
#endif

      allocate(indexNN(max(N_Nneigh,Nei_z),phase),stat=alloc_check)
      if (alloc_check.ne.0) write(6,'(a)') 'out of memory cannot allocate indexNN'
      allocate(tabledist(max(N_Nneigh,Nei_z),phase),stat=alloc_check)
      if (alloc_check.ne.0) write(6,'(a)') 'out of memory cannot allocate table of distance'
      indexNN=0
      tabledist=0.0d0

! calculate distance of NN and table of nearest neighbors
#ifdef CPP_MPI
      if (irank.eq.0) write(6,'(a)') 'Calculating the table of distances'
#else
      write(6,'(a)') 'Calculating the table of distances'
#endif
      if ((count(motif%i_m).eq.1).and.(phase.eq.1)) then
       call Tdist(net,N_Nneigh,world,phase,motif,tabledist(:,1))
      else
       call Tdist(net,N_Nneigh,Nei_z,Nei_il,world,phase,motif,tabledist(:,:))
      endif

#ifdef CPP_MPI
      if (irank.eq.0) write(6,'(a)') 'Calculating the number of neighbors for each shells'
#else
      write(6,'(a)') 'Calculating the number of neighbors for each shells'
#endif
      if (phase.eq.1) then
       call numneigh(N_Nneigh,tabledist(:,1),net,Nei_z,world,motif,phase,indexNN(:,1))
       else
       call numneigh(N_Nneigh,tabledist(:,:),net,Nei_z,Nei_il,world,motif,phase,indexNN(:,:))
      endif

! prepare the lattice
#ifdef CPP_MPI
      if (irank.eq.0) write(6,'(a)') 'initializing the spin structure'
#else
      write(6,'(a,/)') 'initializing the spin structure'
#endif
       call InitSpin(net,motif,world,state)

! allocate table of neighbors and masque
       if (phase.eq.1) then
        tot_N_Nneigh=sum(indexNN(1:N_Nneigh,1),1)
        else
        tot_N_Nneigh=sum(indexNN(1:N_Nneigh,1:phase))
       endif

       allocate(masque(tot_N_Nneigh+1,dim_lat(1),dim_lat(2),dim_lat(3)),stat=alloc_check)
       if (alloc_check.ne.0) write(6,'(a)') 'out of memory cannot allocate masque'

#ifdef CPP_MPI
       allocate(tableNN(4,tot_N_Nneigh,Xstart:Xstop,Ystart:Ystop,Zstart:Zstop,count(motif%i_m)),stat=alloc_check)
       if (alloc_check.ne.0) write(6,'(a)') 'out of memory cannot allocate tableNN'
#else
       allocate(tableNN(4,tot_N_Nneigh,dim_lat(1),dim_lat(2),dim_lat(3),count(motif%i_m)),stat=alloc_check)
       if (alloc_check.ne.0) write(6,'(a)') 'out of memory cannot allocate tableNN'
#endif
       tableNN=0

!-------------------------------------------------
!mapping of the neighbours
#ifdef CPP_MPI
if (irank.eq.0) write(6,'(a)') 'Calculating the table of neighbors (this can take time)'
#else
write(6,'(a)') 'Calculating the table of neighbors (this can take time)'
#endif
      if (size(world).eq.1) then
       call mapping(net,tot_N_Nneigh,tabledist(:,1),N_Nneigh,Nei_il,phase,motif,indexNN(:,1), &
     & tableNN(:,:,:,1,1,1))
      elseif (size(world).eq.2) then
       if ((phase.eq.1).and.(count(motif%i_m).eq.1)) then
        call mapping(net,tot_N_Nneigh,tabledist(:,1),N_Nneigh,Nei_il,phase,motif,indexNN(:,1), &
     & tableNN(:,:,:,:,1,1))
       elseif (phase.eq.1) then
        call mapping(net,tot_N_Nneigh,tabledist(:,1),N_Nneigh,Nei_il,phase,motif,indexNN(:,1), &
     & tableNN(:,:,:,:,1,:))
       endif
       if (phase.eq.2) call mapping(net,tot_N_Nneigh,tabledist(:,:),N_Nneigh,Nei_il,phase,motif,indexNN(:,:), &
     & tableNN(:,:,:,:,1,:))
      else
       if (phase.eq.1) then
        call mapping(net,tot_N_Nneigh,tabledist(:,1),N_Nneigh,Nei_il,phase,motif,indexNN(:,1), &
     & tableNN(:,:,:,:,:,:))
       else
        call mapping(net,tot_N_Nneigh,tabledist(:,:),N_Nneigh,Nei_il,Nei_z,phase,motif,indexNN(:,:), &
     & tableNN(:,:,:,:,:,:))
       endif

      endif
#ifdef CPP_MPI
if (irank.eq.0) write(6,'(a)') 'done'
#else
write(6,'(a)') 'done'
#endif

!-------------------------------------------------
! take care of the periodic boundary conditions
      if (all(Periodic_log)) then
       masque=1
      else
#ifdef CPP_MPI
if (irank.eq.0) write(6,'(a)') 'seting up the masque for non-periodic systems'
#else
write(6,'(a)') 'seting up the masque for non-periodic systems'
#endif
       call periodic(indexNN(:,1),sum(indexNN(:,1))+1,tabledist(:,1),N_Nneigh,Periodic_log(:),masque,spin,tableNN)
      endif

!check for the structure.xyz file for shape modification
      inquire (file='structure.xyz',exist=i_usestruct)
      if (i_usestruct) call user_def_struct(masque,spin, &
     & dim_lat,count(motif%i_m),tot_N_Nneigh+1,net)

! setup the different parameters for the interaction DM, 4-spin... 
! hard part: setup DM interaction
! I send an example neighbor table to setup the DM

      if (i_DM) then
#ifdef CPP_MPI
if (irank.eq.0) write(6,'(a)') 'Calculating the DM vectors for each shells'
#else
write(6,'(a)') 'Calculating the DM vectors for each shells'
#endif
       i=count(abs(DM(:,1))>1.0d-8)
       allocate(DM_vector(sum(indexNN(1:i,1)),3,phase))
       if (size(world).eq.1) then
        DM_vector=setup_DM(sum(indexNN(1:i,1)),i,indexNN(:,1),net,motif,dim_lat,world,phase, &
         1)
       elseif (size(world).eq.2) then
        DM_vector=setup_DM(sum(indexNN(1:i,1)),i,indexNN(:,1),net,motif,dim_lat,world,phase, &
         1,1)
       else
        if (phase.eq.1) then
         DM_vector=setup_DM(sum(indexNN(1:i,1)),i,indexNN(:,1),net,motif,dim_lat,world,phase, &
         1,1,1)
        else
         DM_vector=setup_DM(sum(indexNN(1:i,1)),i,indexNN(:,1),net,motif,dim_lat,world,phase, &
         1,1)
        endif
       endif

! the DM and the neighbors have to turn in the same direction. therefore the matrix of the neighbors
! have to be rearranged
#ifdef CPP_MPI
if (irank.eq.0) write(6,'(a)') 'Re-aranging the position of the DM vectors'
#else
write(6,'(a)') 'Re-aranging the position of the DM vectors'
#endif
       if (size(world).eq.1) then
        call arrange_neigh(DM_vector(:,:,1),tableNN(:,:,:,1,1,:),indexNN(:,1),dim_lat,net,phase,tot_N_Nneigh,world,tabledist(:,1),masque)
       elseif (size(world).eq.2) then
        if (phase.eq.1) then
         call arrange_neigh(DM_vector(:,:,1),tableNN(:,:,:,:,1,:),indexNN(:,1),dim_lat,net,phase,tot_N_Nneigh,world,tabledist(:,1),masque(:,:,:,1))
         else
         call arrange_neigh(DM_vector,tableNN(:,:,:,:,1,:),indexNN(:,:),dim_lat,net,phase,tot_N_Nneigh,world,tabledist(:,:),masque)
        endif
       else
        if (phase.eq.1) then
         call arrange_neigh(DM_vector(:,:,1),tableNN(:,:,:,:,:,:),indexNN(:,1),dim_lat,net,phase,tot_N_Nneigh,world,tabledist(:,1),masque)
        else
         call arrange_neigh(DM_vector(:,:,1),tableNN(:,:,:,:,:,:),indexNN(:,:),dim_lat,net,phase,tot_N_Nneigh,world,tabledist(:,:),masque)
        endif
       endif
      else
       allocate(DM_vector(1,1,1))
       DM_vector=0.0d0
      endif

      if (i_four) then
#ifdef CPP_MPI
if (irank.eq.0) write(6,'(a)') 'Calculating the position of the corners of the unit cell for the 4-spin'
#else
write(6,'(a)') 'Calculating the position of the corners of the unit cell for the 4-spin'
#endif
       call givemecorners(net,abs(order_zaxis(net)))
      else
       allocate(corners(1,1))
      endif

#ifdef CPP_MPI
if (irank.eq.0) write(6,'(a)') 'dealing with the z-direction'
#else
write(6,'(a)') 'dealing with the z-direction'
#endif
!! check the z direction structure
!      call setup_zdir(phase,tot_N_Nneigh,motif)

!!--------------------------------------------------------------------

!c Calculation of the dispertion in the full BZ
! this is here because we have the table of distances
      if (dispersion) then
       call fullBZ(tabledist,N_Nneigh,Nei_z,phase)
      endif

      deallocate(tabledist)

!!!! part of the FFT dipole dipole interaction
#ifndef CPP_BRUTDIP
! prepare the demag tensor
      if (i_dip) then
       write(6,'(a)') "preparing spatial contribution to dipole convolution "
       call setup_dipole(dim_lat,Periodic_log,net,count(motif%i_m),size(world))
       do k=1,dim_lat(3)*count(motif%i_m)
        do j=1,dim_lat(2)
         do i=1,dim_lat(1)
         mmatrix(1:3,i,j,k)=spin(4:6,i,j,k,mod(k-1,count(motif%i_m))+1)*spin(7,i,j,k,1)
         enddo
        enddo
       enddo
       write(6,'(a)') 'FFT dipole dipole is set up'
       else
       allocate(ntensor(1,1,1,1),mmatrix(1,1,1,1),ctrans(1,1,1),rtrans(1,1,1), &
     & hcomplex(1,1,1,1),mcomplex(1,1,1,1),hreal(1,1,1,1))

       mmatrix=0.0d0
       rtrans=0.0d0
       hreal=0.0d0
       ctrans=dcmplx(0.d0,0.d0)
       mcomplex=dcmplx(0.d0,0.d0)
       hcomplex=dcmplx(0.d0,0.d0)
       ntensor=dcmplx(0.d0,0.d0)
      endif
!#endif
#endif

#ifdef CPP_MPI
if (irank.eq.0) write(6,'(/,a/)') 'the setup of the simulation is over'
#else
write(6,'(/,a,/)') 'the setup of the simulation is over'
#endif

!!!!!!!!!!!!!! end of the setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Simulation fo the STM stuff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Does the input file of Tobias is there
      inquire (file='config.dat',exist=spstmL)
      if (spstmL) then
      write(6,'(a)') 'The program of spmstm will be used to draw spins'
!Does spmstm will be used without MC
      inquire (file='input.dat',exist=spstmonly)
      if (spstmonly) then
      write(6,'(a)') 'The program of spmstm will be used to draw spins with & 
     &  the configuration specified in input.dat'
      call spstm
      write(6,'(a)') 'STM images were done'
      stop
      endif
      endif

!calculate topocharge only
      inquire (file='topoonly.dat',exist=topoonly)
      inquire(file='SpinSTMi.dat',exist=i_exi)
      if ((topoonly).and.(i_exi)) then
       write(6,*) 'calculate the topological charge only with &
     &  the configuration specified in SpinSTMi.dat'

      allocate(map_vort(dim_lat(1),dim_lat(2),dim_lat(3),3))
      allocate(map_toto(dim_lat(1),dim_lat(2),dim_lat(3)))
      map_toto=0.0d0
      map_vort=0.0d0
! check for a thrid non periodic dimension inside table hexa
      call InitSpin
      if (size(world).eq.2) then
        if (count(motif%i_m).eq.1) then
!         call topo_map(spin(4:6,:,:,1,1),map_vort(:,:,1,:),map_toto(:,:,1))
        else
!         call topo_map(spin(4:6,:,:,1,:),map_vort(:,:,1,:),map_toto(:,:,1))
        endif
      endif
      write(6,'(3(a,f10.6,3x))') 'Q=',sum(map_toto), &
       'Q+= ',sum(map_toto,mask=map_toto>0.0d0), &
       'Q-= ',sum(map_toto,mask=map_toto<0.0d0)
      deallocate(map_vort,map_toto)
      stop
      endif

      end subroutine setup_simu