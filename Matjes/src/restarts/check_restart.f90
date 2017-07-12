    module m_check_restart
    interface check_restart_read
     module procedure restart_read_1D
     module procedure restart_file
     module procedure restart_read_fraction
     module procedure restart_read_data
     module procedure restart_read_replicas
     module procedure restart_read_spin
    end interface
    interface check_restart_write
     module procedure restart_write_fraction
     module procedure restart_write_data
     module procedure restart_write_1D
     module procedure restart_write_replicas
     module procedure restart_write_spin
    end interface
    contains
!   modules that restart from the different previous configurations
!
!
!-------------------------------------------------------------------
! check if a file is present or absent
    subroutine restart_file(fname,i_test)
    implicit none
    character(len=*), intent(in) :: fname
    logical, intent(inout) :: i_test
! internal
    character(len=50) :: name_int
    integer :: i

    name_int=trim(fname)

    inquire(file=name_int,exist=i_test)

    end subroutine restart_file


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!=================  restart with a given spin structure   ==========
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine restart_read_spin(fname,spin,shape_spin)
    implicit none
    integer, intent(in) :: shape_spin(:)
    character(len=*), intent(in) :: fname
    real(kind=8), intent(inout) :: spin(:,:,:,:,:)
! internal
    character(len=50) :: name_int
    integer :: j,i_x,i_y,i_z,i_m
    integer, parameter :: io=156

    name_int=trim(fname)

    write(6,'(/,a,/)') 'reading spin structure from spin.in'

    open(io,file=name_int,action='read',form='formatted',status='old')

    do i_m=1,shape_spin(5)
     do i_z=1,shape_spin(4)
      do i_y=1,shape_spin(3)
       do i_x=1,shape_spin(2)

         read(io,*) (spin(j,i_x,i_y,i_z,i_m),j=1,7)

       enddo
      enddo
     enddo
    enddo

    close(io)

    end subroutine restart_read_spin

    subroutine restart_write_spin(fname,spin,shape_spin)
    implicit none
    integer, intent(in) :: shape_spin(:)
    character(len=*), intent(in) :: fname
    real(kind=8), intent(in) :: spin(:,:,:,:,:)
! internal
    character(len=50) :: name_int
    integer :: j,i_x,i_y,i_z,i_m
    integer, parameter :: io=156

    name_int=trim(fname)

    write(6,'(/,a,/)') 'write the spin sctructure in spin.out'

    open(io,file=name_int,action='write',form='formatted',status='old',position='append')

    do i_m=1,shape_spin(5)
     do i_z=1,shape_spin(4)
      do i_y=1,shape_spin(3)
       do i_x=1,shape_spin(2)

         write(io,'(7E20.10E3)') (spin(j,i_x,i_y,i_z,i_m),j=1,7)

       enddo
      enddo
     enddo
    enddo

    close(io)

    end subroutine restart_write_spin


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!=================    restart with replicas              ==========
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine restart_read_replicas(fname,replicas,shape_spin,isize)
    implicit none
    integer, intent(in) :: isize,shape_spin(:)
    character(len=*), intent(in) :: fname
    real(kind=8), intent(inout) :: replicas(:,:,:,:,:,:)
! internal
    character(len=50) :: name_int
    integer :: i,j,i_x,i_y,i_z,i_m
    integer, parameter :: io=156

    name_int=trim(fname)

    write(6,'(/,a,/)') 'initialize the replicas from replicas.in'

    open(io,file=name_int,action='read',form='formatted',status='old')

    do i=1,isize
     read(io,*)
     do i_m=1,shape_spin(5)
      do i_z=1,shape_spin(4)
       do i_y=1,shape_spin(3)
        do i_x=1,shape_spin(2)

         read(io,*) (replicas(j,i_x,i_y,i_z,i_m,i),j=1,7)

        enddo
       enddo
      enddo
     enddo
    enddo

    close(io)

    end subroutine restart_read_replicas

    subroutine restart_write_replicas(fname,replicas,shape_spin,isize)
    implicit none
    integer, intent(in) :: isize,shape_spin(:)
    character(len=*), intent(in) :: fname
    real(kind=8), intent(in) :: replicas(:,:,:,:,:,:)
! internal
    character(len=50) :: name_int
    integer :: i,j,i_x,i_y,i_z,i_m
    integer, parameter :: io=156

    name_int=trim(fname)

    write(6,'(/,a,/)') 'save the replicas in replicas.out'

    open(io,file=name_int,action='write',form='formatted',status='old',position='append')

    do i=1,isize
     write(io,'(a,2x,I6)') 'replicas number',i
     do i_m=1,shape_spin(5)
      do i_z=1,shape_spin(4)
       do i_y=1,shape_spin(3)
        do i_x=1,shape_spin(2)

         write(io,'(7E20.10E3)') (replicas(j,i_x,i_y,i_z,i_m,i),j=1,7)

        enddo
       enddo
      enddo
     enddo
    enddo

    close(io)

    end subroutine restart_write_replicas



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!=================    restart with the measured data      ==========
!  M_sum, E_sum...
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine restart_read_data(fname,E_sum_av,E_sq_sum_av,E_err_av, &
      &   M_sum_av,M_sq_sum_av,M_err_av, qeulerp_av,qeulerm_av,Q_sq_sum_av,vortex_av, &
      &   isize)
    implicit none
    integer, intent(in) :: isize
    character(len=*), intent(in) :: fname
    real(kind=8), intent(inout) :: E_sum_av(:),E_sq_sum_av(:),E_err_av(:),qeulerp_av(:),qeulerm_av(:),Q_sq_sum_av(:)
    real(kind=8), intent(inout) :: M_sum_av(:,:),M_sq_sum_av(:,:),M_err_av(:,:),vortex_av(:,:)
! internal
    character(len=50) :: name_int
    integer :: i,j
    integer, parameter :: io=156

    name_int=trim(fname)

    write(6,'(/,a,/)') 'initialize the measured data (M_sum, E, E**2) from file'

    open(io,file=name_int,action='read',form='formatted',status='old')

! energies
    do i=1,isize
        read(io,*) E_sum_av(i),E_sq_sum_av(i),E_err_av(i)
    enddo
! magnetizations
    do i=1,isize
        read(io,*) (M_sum_av(j,i),j=1,3),(M_sq_sum_av(j,i),j=1,3),(M_err_av(j,i),j=1,3)
    enddo
! topological charges
    do i=1,isize
        read(io,*) qeulerp_av(i),qeulerm_av(i),Q_sq_sum_av(i)
    enddo
! vortices
    do i=1,isize
        read(io,*) (vortex_av(j,i),j=1,3)
    enddo

    close(io)

    end subroutine restart_read_data

    subroutine restart_write_data(fname,E_sum_av,E_sq_sum_av,E_err_av, &
      &   M_sum_av,M_sq_sum_av,M_err_av, qeulerp_av,qeulerm_av,Q_sq_sum_av,vortex_av, &
      &   isize)
    implicit none
    integer, intent(in) :: isize
    character(len=*), intent(in) :: fname
    real(kind=8), intent(in) :: E_sum_av(:),E_sq_sum_av(:),E_err_av(:),qeulerp_av(:),qeulerm_av(:),Q_sq_sum_av(:)
    real(kind=8), intent(in) :: M_sum_av(:,:),M_sq_sum_av(:,:),M_err_av(:,:),vortex_av(:,:)
! internal
    character(len=50) :: name_int
    integer :: i,j
    integer, parameter :: io=156

    name_int=trim(fname)

    write(6,'(/,a,/)') 'save the measured data (M_sum, E, E**2) in measured_data.out'

    open(io,file=name_int,action='write',form='formatted',status='old',position='append')

! energies
    do i=1,isize
        write(io,'(3E20.10E3)') E_sum_av(i),E_sq_sum_av(i),E_err_av(i)
    enddo
! magnetizations
    do i=1,isize
        write(io,'(9E20.10E3)') (M_sum_av(j,i),j=1,3),(M_sq_sum_av(j,i),j=1,3),(M_err_av(j,i),j=1,3)
    enddo
! topological charges
    do i=1,isize
        write(io,'(3E20.10E3)') qeulerp_av(i),qeulerm_av(i),Q_sq_sum_av(i)
    enddo
! vortices
    do i=1,isize
        write(io,'(3E20.10E3)') (vortex_av(j,i),j=1,3)
    enddo

    close(io)

    end subroutine restart_write_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!=================    restart of the parallel tempering    ==========
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine restart_read_fraction(fname,nup,ndown,success_PT,label_world,image_temp,isize)
    implicit none
    integer, intent(in) :: isize
    character(len=*), intent(in) :: fname
    integer, intent(inout) :: image_temp(:)
    real(kind=8), intent(inout) :: nup(:),ndown(:),success_PT(:),label_world(:)
! internal
    character(len=50) :: name_int
    integer :: i
    integer, parameter :: io=156

    name_int=trim(fname)

    write(6,'(/,a,/)') 'initialize the fraction from file'

    open(io,file=name_int,action='read',form='formatted',status='old')
    do i=1,isize
        read(io,*) nup(i),ndown(i),success_PT(i),label_world(i),image_temp(i)
     enddo
    close(io)

    end subroutine restart_read_fraction

    subroutine restart_write_fraction(fname,nup,ndown,success_PT,label_world,image_temp,isize)
    implicit none
    integer, intent(in) :: isize
    character(len=*), intent(in) :: fname
    integer, intent(in) :: image_temp(:)
    real(kind=8), intent(in) :: nup(:),ndown(:),success_PT(:),label_world(:)
! internal
    character(len=50) :: name_int
    integer :: i
    integer, parameter :: io=156

    name_int=trim(fname)

    write(6,'(/,a,/)') 'save the fraction data in fraction.out'

    open(io,file=name_int,action='write',form='formatted',status='old',position='append')
    do i=1,isize
        write(io,'(4E20.10E3,I6)') nup(i),ndown(i),success_PT(i),label_world(i),image_temp(i)
     enddo
    close(io)

    end subroutine restart_write_fraction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!=================    restart of the temperature    ================
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! reread of the temperature
    subroutine restart_read_1D(fname,kt_all,isize)
    use m_constants, only : k_b
    implicit none
    integer, intent(in) :: isize
    character(len=*), intent(in) :: fname
    real(kind=8), intent(inout) :: kt_all(isize)
! internal
    integer :: i
    integer, parameter :: io=156
    logical :: i_exist
    character(len=50) :: name_int

    name_int=trim(fname)
    i_exist=.False.

    inquire(file=name_int,exist=i_exist)
    if (i_exist) then

     write(6,'(/,a,/)') 'restarting the calculations with the temperatures from temperature.in'
     open(io,file=name_int,form='formatted',status='old')
     do i=1,isize
      read(io,*) kt_all(i)
      kt_all(i)=kt_all(i)*k_b
     enddo
     close(io)
    endif

    end subroutine restart_read_1D


! write the temperature
    subroutine restart_write_1D(fname,kt_all,isize)
    use m_constants, only : k_b
    implicit none
    integer, intent(in) :: isize
    character(len=*), intent(in) :: fname
    real(kind=8), intent(in) :: kt_all(isize)
! internal
    integer :: i
    integer, parameter :: io=156
    character(len=50) :: name_int

    name_int=trim(fname)

    write(6,'(/,a,/)') 'writting the final temperature set in temperature.out'

    open(io,FILE=name_int,action='write',status='unknown',form='formatted')

    do i=1,isize
         write(io,'(E20.10E3)') kT_all(i)/k_b
    enddo

    close(io)


    end subroutine restart_write_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!===================================================================
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end module m_check_restart
