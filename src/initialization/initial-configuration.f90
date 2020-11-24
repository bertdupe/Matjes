module m_init_config
use m_derived_types
use m_init_spiral
use m_init_DW
use m_init_Sk
use m_init_sky_lin, only: init_sky_lin
use m_init_Sklattice
use m_init_constant
use m_init_punch, only: init_punch
use m_io_utils
use m_init_heavyside
use m_io_files_utils
use m_init_random_config
use m_lattice, only : my_order_parameters
use m_convert

private
public :: init_config_new, init_config
contains


subroutine init_config_new(lat,fname_in)
    !loops through all order parameters and checks if there is some initial configuration supplied from the file input
    use m_type_lattice,only: number_different_order_parameters,order_parameter_name
    type(lattice), intent(inout) :: lat
    character(len=*), intent(in),optional :: fname_in
    integer :: io

    character(*),parameter              :: fname_default='init.config'
    character(:), allocatable           :: fname
    logical :: exists
    integer :: i
    integer :: dim_mode
    real(8),pointer ::  state(:)

    if(present(fname_in))then
        fname=fname_in
    else
        fname=fname_default
    endif
    inquire(file=fname,exist=exists)
    if(.not.exists)then
        write(*,'(///A)') 'FILE:',fname," NOT FOUND"
        write(*,'(A///)') 'SKIPPING init_config'
        return
    endif
    io=open_file_read(fname)
    do i=1,number_different_order_parameters
        dim_mode=lat%get_order_dim(i,ignore=.true.)
        if(dim_mode<1) cycle
        Call lat%set_order_point(i,state)
        Call set_init(io,fname,lat,trim(order_parameter_name(i)),dim_mode,state)
        nullify(state)
    enddo
    call close_file(fname,io)
end subroutine

subroutine set_init(io,fname,lat,ordname,dim_mode,state)
    integer,intent(in)              :: io       !init-file io-unit
    character(*),intent(in)         :: fname    !init-file name 
    type(lattice), intent(in)       :: lat      !entire lattice containing geometric information
    character(*),intent(in)         :: ordname  !name of the order parameter
    integer,intent(in)              :: dim_mode !dimension of the order parameter in each cell
    real(8),pointer,intent(inout)   :: state(:) !pointer the the order parameter

    character(30)                   :: configuration

    configuration=""
    call get_parameter(io,fname,'configuration_'//ordname,configuration)

    select case (adjustl(configuration))
        case('spiral')
            call init_spiral(io,fname,lat,ordname,dim_mode,state)
        case('domainwall')
            call init_DW(io,fname,lat,ordname,dim_mode,state)
        case('heavyside')
            call init_heavyside(lat,dim_mode,state)
        case('skyrmion')
            call init_spiral(io,fname,lat,ordname,dim_mode,state)
            call init_Sk(io,fname,lat,ordname,dim_mode,state)
        case('skyrmion_lin')
            Call init_Sky_lin(io,fname,lat,ordname,dim_mode,state)
        case('skyrmionla')
            call init_spiral(io,fname,lat,ordname,dim_mode,state)
            call init_Sk_lattice(io,fname,lat,ordname,dim_mode,state)
        case('random')
            write(6,'(3a)') 'random configuration for ',ordname,' was chosen'
            call init_random_config(dim_mode,state)
        case default
            write(6,'(3a)') 'No initial configuration for ',ordname,' was found'
    end select
end subroutine 

!!!!!!!!!!!
! find the different starting configurations for the order parameters
subroutine init_config(fname,my_lattice,my_motif,ext_param)
implicit none
character(len=*), intent(in) :: fname
type (lattice), intent(inout) :: my_lattice
type(t_cell), intent(in) :: my_motif
type (simulation_parameters), intent(in) :: ext_param
! internal variables
integer :: io,N_mode,index_mode,i
character(len=30) :: seed_name,mode_name,configuration
integer :: nconfig

nconfig=0
seed_name='configuration_'
io=open_file_read(fname)
N_mode=size(my_order_parameters)

do i=1,N_mode

  configuration='constant'
  mode_name=convert(seed_name,my_order_parameters(i)%name)
  call get_parameter(io,fname,mode_name,configuration)

  select case (adjustl(configuration))
    case('spiral')
      call init_spiral(io,fname,my_lattice,my_motif,my_order_parameters(i)%name,my_order_parameters(i)%start,my_order_parameters(i)%end)
    case('domainwall')
      call init_DW(my_lattice,my_motif,my_order_parameters(i)%start,my_order_parameters(i)%end)
    case('heavyside')
      call init_heavyside(my_lattice,my_motif,my_order_parameters(i)%start,my_order_parameters(i)%end)
    case('skyrmion')
      call init_spiral(io,fname,my_lattice,my_motif,my_order_parameters(i)%name,my_order_parameters(i)%start,my_order_parameters(i)%end)
      call init_Sk(io,fname,my_lattice,my_motif,my_order_parameters(i)%name,my_order_parameters(i)%start,my_order_parameters(i)%end)
    case('skyrmion_lin')
      Call init_Sky_lin(io,fname,my_lattice,my_motif,my_order_parameters(i)%start,my_order_parameters(i)%end)
    case('skyrmionla')
      call init_spiral(io,fname,my_lattice,my_motif,my_order_parameters(i)%name,my_order_parameters(i)%start,my_order_parameters(i)%end)
      call init_Sk_lattice(io,fname,my_lattice,my_motif,my_order_parameters(i)%name,my_order_parameters(i)%start,my_order_parameters(i)%end)
    case('random')
      write(6,'(a)') 'random configuration was choosen'
      call init_random_config(my_lattice,my_motif,my_order_parameters(i)%start,my_order_parameters(i)%end)
    case default
      call init_constant_config(my_lattice,my_order_parameters(i)%name,my_order_parameters(i)%start,my_order_parameters(i)%end,ext_param)
  end select
  Call init_punch(io,fname,my_lattice,my_motif,my_order_parameters(i)%name,my_order_parameters(i)%start,my_order_parameters(i)%end)

  !Call copy_init(my_lattice,my_order_parameters(i)%name,my_order_parameters(i)%start,my_order_parameters(i)%end)
enddo


call close_file(fname,io)

end subroutine init_config

subroutine copy_init(my_lattice,mode_name,start,end)
 !copy initialization to specific orderparameter in lattice (has to be moved further upwards, i.e. directly initialize the new array)
    type(lattice),intent(inout)     ::  my_lattice
    character(len=30),intent(in)    ::  mode_name
    integer,intent(in)              :: start,end

    integer                         :: ix,iy,iz,im
    if(adjustl(mode_name)=='magnetic')then
        do im=1,size(my_lattice%ordpar%l_modes,4)
            do iz=1,size(my_lattice%ordpar%l_modes,3)
                do iy=1,size(my_lattice%ordpar%l_modes,2)
                    do ix=1,size(my_lattice%ordpar%l_modes,1)
                        my_lattice%M%l_modes(ix,iy,iz,im)%w=my_lattice%ordpar%l_modes(ix,iy,iz,im)%w(start:end)
                    enddo
                enddo
            enddo
        enddo
    elseif(adjustl(mode_name)=='Efield')then
        do im=1,size(my_lattice%ordpar%l_modes,4)
            do iz=1,size(my_lattice%ordpar%l_modes,3)
                do iy=1,size(my_lattice%ordpar%l_modes,2)
                    do ix=1,size(my_lattice%ordpar%l_modes,1)
                        my_lattice%E%l_modes(ix,iy,iz,im)%w=my_lattice%ordpar%l_modes(ix,iy,iz,im)%w(start:end)
                    enddo
                enddo
            enddo
        enddo
    elseif(adjustl(mode_name)=='Bfield')then
        do im=1,size(my_lattice%ordpar%l_modes,4)
            do iz=1,size(my_lattice%ordpar%l_modes,3)
                do iy=1,size(my_lattice%ordpar%l_modes,2)
                    do ix=1,size(my_lattice%ordpar%l_modes,1)
                        my_lattice%B%l_modes(ix,iy,iz,im)%w=my_lattice%ordpar%l_modes(ix,iy,iz,im)%w(start:end)
                    enddo
                enddo
            enddo
        enddo
    elseif(adjustl(mode_name)=='temperature')then
        write(*,*) "NOT READING TEMPERATURE FROM INIT, PUT BACK IN"
    !    do im=1,size(my_lattice%ordpar%l_modes,4)
    !        do iz=1,size(my_lattice%ordpar%l_modes,3)
    !            do iy=1,size(my_lattice%ordpar%l_modes,2)
    !                do ix=1,size(my_lattice%ordpar%l_modes,1)
    !                    my_lattice%T%l_modes(ix,iy,iz,im)%w=my_lattice%ordpar%l_modes(ix,iy,iz,im)%w(start:end)
    !                enddo
    !            enddo
    !        enddo
    !    enddo
    else
        write(*,*) 'mode_name: ',mode_name
        STOP 'unexpected mode_name in initial-configuration'
    endif

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!
! function that finds the index of the mode of interest
!!!!!!!!!!!!!!!!!!!!!!!!!
integer function find_mode(mode_name)
implicit none
character(len=*), intent(in) :: mode_name
!internal
integer :: i,n_mode_total

find_mode=-1
n_mode_total=size(my_order_parameters)

do i=1,n_mode_total
   if (trim(my_order_parameters(i)%name).eq.trim(mode_name)) find_mode=i
enddo

if (find_mode.eq.-1) stop 'ERROR: Mode can not be found in initial-config'

end function find_mode

end module m_init_config
