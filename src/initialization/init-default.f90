module m_init_default
use m_input_types, only: extpar_input
implicit none
private
public init_default
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Default initialization for all order parameters
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_default(ordname,dim_mode,state,extpar_io)
    character(*),intent(in)         :: ordname  !name of the order parameter
    integer,intent(in)              :: dim_mode !dimension of the order parameter in each cell
    real(8),pointer,intent(inout)   :: state(:) !pointer the the order parameter
    type(extpar_input),intent(in)   :: extpar_io

    real(8),pointer     ::  state_v(:,:)
    integer             ::  i
    real(8)             ::  def_val(dim_mode)

    state_v(1:dim_mode,1:size(state)/dim_mode)=> state
    call get_default_value(ordname,dim_mode,extpar_io,def_val)
    do i=1,size(state_v,2)
        state_v(:,i)=def_val
    enddo
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine that finds the default value,based on ext_param and ordname
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_default_value(ordname,dim_mode,extpar_io,def_val)
    character(len=*), intent(in)        :: ordname
    integer,intent(in)                  :: dim_mode
    type(extpar_input),intent(in)       :: extpar_io
    real(8),intent(out)                 :: def_val(dim_mode)

    integer ::  i
    
    select case (ordname)
        case('magnetic')
            do i=1,dim_mode/3
                def_val(1+(i-1)*3:i*3)=[0.0d0,0.0d0,1.0d0]
            enddo
        case('Efield')
            def_val=extpar_io%E
        case('Bfield')
            def_val=extpar_io%H
        case('temperature')
            def_val=extpar_io%T(1)
        case('phonon')
            def_val=0.0
        case('wavefunction')
            def_val=0.0
        case default
            write(*,'(//2A)') "COULD NOT OBTAIN DEFAULT VALUE FOR ORDER: ",ordname
            write(*,'(A)') "IMPLEMENT A DEFAULT VALUE"
            stop 'error in get_default_value'
    end select

end subroutine

end module 
