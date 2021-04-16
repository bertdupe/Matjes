module m_mc_therm_val
use m_mc_exp_val, only: exp_val
use m_mc_track_val, only: track_val
private
public :: therm_val,therm_set
public :: thermo_gatherv
public :: thermo_print, thermo_print_init
!public :: exp_val, measure_print_thermo, measure_add,measure_eval,print_av
!public :: measure_print_thermo_init
!!public MPI_routines
!public :: measure_scatterv, measure_gatherv, measure_bcast, measure_reduce
integer,protected       :: MPI_custom_type          ! mpi variable for direct mpi operations with this type 
logical,protected       :: MPI_custom_set=.false.   !check is MPI_custom_type is set
!
!!THESE ENTRIES HAVE TO BE UPDATED EVERYTIME EXP_VAL IS MODIFIED, OTHERWISE MPI-STUFF WILL FAIL
integer,parameter       :: blocks(*)=[1, 1,3,4, 1,3,1,3, 1,1,3,3, 1,1,1,1] !size of each element of therm_val
integer,parameter       :: bnd_real(2) =[ 1, 12] !initial and final entry of real numbers
integer,parameter       :: bnd_cmplx(2)=[13, 16] !initial and final entry of complex numbers
integer,parameter       :: bnd_int(2)=[-1,-1] !initial and final entry of complex numbers
integer,parameter       :: N_entry=size(blocks)

type therm_val                                                                                             
    !! contribution of the different energy parts
    real(8) :: kt=0.0d0 !temperature should always be at first position so it does not get added measure_reduce
    ! thermodynamical quantities
    real(8) :: C=0.0d0
    real(8) :: chi_M(3)=0.0d0
    real(8) :: chi_Q(4)=0.0d0
    ! Energy/Magnetization and errors
    real(8) :: E=0.0d0
    real(8) :: M(3)=0.0d0
    real(8) :: E_err=0.0d0
    real(8) :: M_err(3)=0.0d0
    ! energy and so on
    real(8) :: qeulerp=0.0d0
    real(8) :: qeulerm=0.0d0
    real(8) :: vortex(3)=0.0d0
    real(8) :: chi_l(3)=0.0d0   !probably quite obsolete
    
    complex(8)  :: MjpMim =cmplx(0.0d0,0.0d0,8) !<Mj+Mi->
    complex(8)  :: MipMip =cmplx(0.0d0,0.0d0,8) !<Mi+Mi+>
    complex(8)  :: MipMim =cmplx(0.0d0,0.0d0,8) !<Mi+Mi-> (is real)
    complex(8)  :: MipMjp =cmplx(0.0d0,0.0d0,8) !<Mi+Mj+>
end type

contains 

subroutine therm_set(this,measure,Cor_log,N_cell_in)
    use m_constants, only : pi
    !routine which calculates all thermodynamic quantities in this from the 
    !tracked expectation values
    class(therm_val),intent(out)    :: this
    class(exp_val),intent(in)       :: measure
    logical, intent(in)             :: Cor_log
    integer,intent(in)              :: N_cell_in
    !internal
    real(8) :: av_site,av_Nadd,av_kt
    
    !set normalization variables
    av_site=1.0d0/real(N_cell_in,8)
    av_Nadd=1.0d0/real(measure%N_add,8)
    av_kt=1.0d0/measure%kt

    this%kt=measure%kt

    this%C    =(measure%E_sq*av_Nadd-(measure%E*av_Nadd)**2)*av_kt**2*av_site
    this%chi_M=(measure%M_sq*av_Nadd-(measure%M*av_Nadd)**2)*av_kt   *av_site
    if (measure%N_add.gt.1)then
        this%E_err=sqrt(abs(measure%E_sq-(measure%E)**2*av_Nadd)/real(measure%N_add-1,8))*av_site
        this%M_err=sqrt(abs(measure%M_sq-(measure%M)**2*av_Nadd)/real(measure%N_add-1,8))*av_site
    endif
    this%E=measure%E*av_Nadd*av_site
    this%M=measure%M*av_Nadd*av_site

    this%qeulerp=measure%Qp*av_Nadd/pi/4.0d0
    this%qeulerm=measure%Qp*av_Nadd/pi/4.0d0
    this%vortex=measure%vortex*av_Nadd/3.0d0/sqrt(3.0d0)

    this%chi_Q(1)=-((measure%Qp-measure%Qm)**2 - measure%Q_sq*av_Nadd/16.0d0/pi**2)*av_kT!/(measure%qeulerp_av-measure%qeulerm_av)! do in the postprocessing 
    this%chi_Q(2)=-(measure%Qp**2-measure%Qp_sq*av_Nadd/16.0d0/pi**2)*av_kT!/measure%qeulerp_av
    this%chi_Q(3)=-(measure%Qm**2-measure%Qm_sq*av_Nadd/16.0d0/pi**2)*av_kT!/measure%qeulerm_av
    this%chi_Q(4)=((measure%Qm_sq*av_Nadd/16.0d0/pi**2-measure%Qm**2)* &
              &    (measure%Qp_sq*av_Nadd/16.0d0/pi**2-measure%Qp**2))*av_kT!/(measure%qeulerm_av*measure%qeulerp_av)

    !fluctuation parameters
    if(allocated(measure%MjpMim_ij)) this%MjpMim=sum(measure%MjpMim_ij)*av_Nadd*av_site  !not entirely sure about this term
    this%MipMip=measure%MipMip*av_Nadd*av_site
    this%MipMim=measure%MipMim*av_Nadd*av_site
    this%MipMjp=measure%MipMjp*av_Nadd*av_site

    if (Cor_log) this%chi_l=measure%N_add !what is this supposed to do?
    Call print_av(this)
end subroutine

subroutine print_av(this)
    use, intrinsic :: iso_fortran_env, only : output_unit
    class(therm_val),intent(inout)    :: this

    write(output_unit,'(5(a,f18.9,2x))') 'M= ',norm2(this%M), &
      & 'E= ',this%E,'Q+= ',this%qeulerp,'Q-= ',-this%qeulerm,'Q= ',this%qeulerp-this%qeulerm
end subroutine

subroutine thermo_print_init(io_unit)
    use m_io_files_utils, only: open_file_write
    integer,intent(inout)   ::  io_unit

    io_unit=open_file_write("EM.dat")
    Write(io_unit,'(34(a,15x))') '# 1:T','2:E_av','3:E_err','4:C','5:M','6:Mx','7:My','8:Mz', &
          & '9:M_err_x','10:M_err_y','11:M_err_z','12:chi_x','13:chi_y','14:chi_z','15:vx', &
          & '16:vy','17:vz','18:qeuler','19:Chi_q','20:Q+','21:Chi_qp','22:Q-','23:Chi_qm', &
          & '24:l_x','25:l_y','26:l_z','27:Chi_QpQm', &
          & '28: Re(Mj+Mi-)','29: Im(Mj+Mi-)','30: Re(Mi+Mi+)','31: Im(Mi+Mi+)','32: Re(Mi+Mi-)', &
          & '33: Re(Mi+Mj+)','34: Im(Mi+Mj+)'
end subroutine

subroutine thermo_print(this,io_unit_in)
    !print thermodynamic quantities 
    use m_constants, only : k_b
    class(therm_val),intent(inout)  :: this(:)
    integer,optional                :: io_unit_in
    integer     ::  io_unit,i
    real(8) ::  Q

    if(present(io_unit_in))then
        io_unit=io_unit_in
    else
        Call thermo_print_init(io_unit)
    endif
    ! write the data into a file
    do i=1,size(this)
        Q=-this(i)%qeulerm+this(i)%qeulerp
        Write(io_unit,'(34(E20.10E3,2x),E20.10E3)') &
         &             this(i)%kT/k_B  ,this(i)%E  , this(i)%E_err, this(i)%C,&
         &             norm2(this(i)%M),this(i)%M, this(i)%M_err,&
         &             this(i)%chi_M   ,this(i)%vortex, Q, this(i)%chi_Q(1),&
         &             this(i)%qeulerp ,this(i)%chi_Q(2), -this(i)%qeulerm, this(i)%chi_Q(3), this(i)%chi_l(1:3), this(i)%chi_Q(4),&
         &             this(i)%MjpMim  ,this(i)%MipMip, real(this(i)%MipMim,8), this(i)%MipMjp
    enddo
    if(.not.present(io_unit_in)) close(io_unit) 
end subroutine

subroutine thermo_gatherv(this,com)
    use mpi_basic
    type(therm_val),intent(inout)   :: this(:)
    class(mpi_distv),intent(in)     :: com
#ifdef CPP_MPI    
    integer                         :: ierr

    if(.not.MPI_custom_set)then
        Call set_MPI_type(blocks,bnd_real,bnd_cmplx,bnd_int,MPI_custom_type)
        MPI_custom_set=.true.
    endif
    Call MPI_Gatherv(this(1),com%cnt(com%id+1),MPI_custom_type,this(1),com%cnt,com%displ,MPI_custom_type,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

end module
