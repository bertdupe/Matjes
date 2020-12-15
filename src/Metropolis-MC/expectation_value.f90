module m_mc_exp_val
use m_constants, only : pi
use m_mc_track_val, only: track_val
use m_fluct,only: fluct_parameters, eval_fluct
implicit none
private
public :: exp_val, measure_print_thermo, measure_add,measure_eval,print_av
public :: measure_print_thermo_init
!public MPI_routines
public :: measure_scatterv,measure_gatherv
integer,parameter       :: N_entry=28
integer,protected       :: MPI_exp_val ! mpi variable for direct mpi exp_val operations                 
logical,protected       :: MPI_set=.false.
type exp_val                                                                                             
    !! contribution of the different energy parts
    !real(8) :: E_decompose(8) !not really used here
    real(8) :: kt=0.0d0
    ! thermodynamical quantities
    real(8) :: C_av=0.0d0
    real(8) :: chi_M(3)=0.0d0
    real(8) :: chi_Q(4)=0.0d0
    ! errors on the different quantities
    real(8) :: E_err_av=0.0d0
    real(8) :: M_err_av(3)=0.0d0
    ! sums
    real(8) :: M_sq_sum_av(3)=0.0d0
    real(8) :: E_sum_av=0.0d0
    real(8) :: E_sq_sum_av=0.0d0
    real(8) :: Q_sq_sum_av=0.0d0
    real(8) :: Qp_sq_sum_av=0.0d0
    real(8) :: Qm_sq_sum_av=0.0d0
    ! energy and so on
    real(8) :: E_av=0.0d0
    real(8) :: qeulerp_av=0.0d0
    real(8) :: qeulerm_av=0.0d0
    real(8) :: M_av(3)=0.0d0
    real(8) :: M_sum_av(3)=0.0d0
    real(8) :: vortex_av(3)=0.0d0
    real(8) :: chi_l(3)=0.0d0

	!<Mj+Mi->
    complex(8)  :: MjpMim_sum=cmplx(0.0d0,0.0d0,8)
    complex(8)  :: MjpMim_av=cmplx(0.0d0,0.0d0,8)

	!<Mi+Mi+>
    complex(8)  :: MipMip_sum=cmplx(0.0d0,0.0d0,8)
    complex(8)  :: MipMip_av=cmplx(0.0d0,0.0d0,8)

	!<Mi+Mi-> (is real)
    complex(8)  :: MipMim_sum=cmplx(0.0d0,0.0d0,8)
    complex(8)  :: MipMim_av=cmplx(0.0d0,0.0d0,8)

	!<Mi+Mj+>
    complex(8)  :: MipMjp_sum=cmplx(0.0d0,0.0d0,8)
    complex(8)  :: MipMjp_av=cmplx(0.0d0,0.0d0,8)

    integer :: N_add=0 !counts how often values have been added

end type

contains 

subroutine measure_set_MPI_type()
#ifdef CPP_MPI    
    use mpi_basic
    integer                                     ::  blocks(N_entry)
    integer(kind=MPI_ADDRESS_KIND)              ::  displ(N_entry)
    integer                                     ::  types(N_entry)
    integer(kind=MPI_ADDRESS_KIND)              ::  LB,extend_real,extend_int,extend_cmplx
    integer(kind=MPI_ADDRESS_KIND)              ::  extend
    integer                                     ::  ierr,i

    blocks=[1, 1,3,4, 1,3, 3,1,1,1,1,1, 1,1,1,3,3,3,3, 1,1,1,1,1,1,1,1, 1]
    types=MPI_DOUBLE_PRECISION
    types(20:)=MPI_DOUBLE_COMPLEX
    types(N_entry)=MPI_INT
    CALL MPI_TYPE_GET_EXTENT(MPI_INT,lb,extend_int,ierr)
    CALL MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION,lb,extend_real,ierr)
    CALL MPI_TYPE_GET_EXTENT(MPI_DOUBLE_COMPLEX,lb,extend_cmplx,ierr)
    displ=0
    do i=1,N_entry-1
        if(types(i)==MPI_DOUBLE_PRECISION)then
            extend=extend_real
        elseif(types(i)==MPI_INT)then
            extend=extend_int
        elseif(types(i)==MPI_DOUBLE_COMPLEX)then
            extend=extend_cmplx
        else
            STOP "UNEXPECTED types"
        endif
        displ(i+1)=displ(i)+blocks(i)*extend
    enddo

    CALL MPI_TYPE_CREATE_STRUCT(N_entry,blocks,displ,types,MPI_exp_val,ierr)
    CALL MPI_TYPE_COMMIT(MPI_exp_val,ierr)
#else
    continue
#endif
end subroutine

subroutine measure_gatherv(this,com)
    use mpi_basic
    type(exp_val),intent(inout)     :: this(:)
    class(mpi_distv),intent(in)     :: com
#ifdef CPP_MPI    
    integer                         :: ierr

    if(.not.MPI_set)then
        Call measure_set_MPI_type()
        MPI_set=.true.
    endif
    Call MPI_Gatherv(this(1),com%cnt(com%id+1),MPI_exp_val,this(1),com%cnt,com%displ,MPI_exp_val,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

subroutine measure_scatterv(this,com)
    use mpi_basic
    type(exp_val),intent(inout)     :: this(:)
    class(mpi_distv),intent(in)     :: com
!    integer,intent(in)              :: displ(com%Np)
!    integer,intent(in)              :: cnt(com%Np)
#ifdef CPP_MPI    
    integer                         :: ierr

    if(.not.MPI_set)then
        Call measure_set_MPI_type()
        MPI_set=.true.
    endif
    Call MPI_Scatterv(this(1),com%cnt,com%displ,MPI_exp_val,this(1),com%cnt(com%id+1),MPI_exp_val,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

subroutine measure_print_thermo_init(io_unit)
    use m_io_files_utils, only: open_file_write
    integer,intent(inout)     ::  io_unit

    io_unit=open_file_write("EM.dat")
    Write(io_unit,'(34(a,15x))') '# 1:T','2:E_av','3:E_err','4:C','5:M','6:Mx','7:My','8:Mz', &
          & '9:M_err_x','10:M_err_y','11:M_err_z','12:chi_x','13:chi_y','14:chi_z','15:vx', &
          & '16:vy','17:vz','18:qeuler','19:Chi_q','20:Q+','21:Chi_qp','22:Q-','23:Chi_qm', &
          & '24:l_x','25:l_y','26:l_z','27:Chi_QpQm', &
		  & '28: Re(Mj+Mi-)','29: Im(Mj+Mi-)','30: Re(Mi+Mi+)','31: Im(Mi+Mi+)','32: Re(Mi+Mi-)', &
		  & '33: Re(Mi+Mj+)','34: Im(Mi+Mj+)'
end subroutine


subroutine measure_print_thermo(this,com,io_unit_in)
    use m_constants, only : k_b
    use mpi_basic, only: mpi_type
    class(exp_val),intent(inout)    :: this(:)
    class(mpi_type),intent(in)      :: com
    integer,optional                :: io_unit_in
    integer     ::  io_unit,i

    if(com%ismas)then
        if(present(io_unit_in))then
            io_unit=io_unit_in
        else
            Call measure_print_thermo_init(io_unit)
        endif
        ! write the data into a file
        do i=1,size(this)
            Write(io_unit,'(34(E20.10E3,2x),E20.10E3)') this(i)%kT/k_B ,this(i)%E_av, this(i)%E_err_av, this(i)%C_av,&
             &             norm2(this(i)%M_av),this(i)%M_av, this(i)%M_err_av,&
             &             this(i)%chi_M, this(i)%vortex_av, -this(i)%qeulerm_av+this(i)%qeulerp_av, this(i)%chi_Q(1),&
             &             this(i)%qeulerp_av, this(i)%chi_Q(2), -this(i)%qeulerm_av, this(i)%chi_Q(3), this(i)%chi_l(:), this(i)%chi_Q(4),&
			 &			   this(i)%MjpMim_av, this(i)%MipMip_av, real(this(i)%MipMim_av,8), this(i)%MipMjp_av
        enddo
        if(.not.present(io_unit_in)) close(io_unit) 
    endif
end subroutine

subroutine measure_add(this,lat,state_prop,Q_neigh,fluct_val)
    use m_topo_commons
    use m_derived_types,only: lattice
    class(exp_val),intent(inout)            :: this
    type(lattice),intent(in)                :: lat
    type(track_val),intent(in)              :: state_prop
    integer,intent(in)                      :: Q_neigh(:,:)
    type(fluct_parameters),intent(in)      :: fluct_val

    !put that into state_prop as well?
    real(8)     :: dumy(5),qeulerp,qeulerm


    this%N_add=this%N_add+1

    this%E_sum_av=this%E_sum_av+state_prop%E_total
    this%E_sq_sum_av=this%E_sq_sum_av+state_prop%E_total**2

    this%M_sum_av=this%M_sum_av+state_prop%Magnetization
    this%M_sq_sum_av=this%M_sq_sum_av+state_prop%Magnetization**2
    ! calculate the topocharge
    dumy=get_charge(lat,Q_neigh)
    qeulerp=dumy(1)
    qeulerm=-dumy(2)

    this%qeulerp_av=this%qeulerp_av+qeulerp
    this%qeulerm_av=this%qeulerm_av+qeulerm
    this%Q_sq_sum_av=this%Q_sq_sum_av+(qeulerp-qeulerm)**2
    this%Qp_sq_sum_av=this%Qp_sq_sum_av+qeulerp**2
    this%Qm_sq_sum_av=this%Qm_sq_sum_av+qeulerm**2
    this%vortex_av=this%vortex_av+dumy(3:5)

    if(fluct_val%l_use) Call eval_fluct(this%MjpMim_sum, &
                                       &this%MipMip_sum, &
                                       &this%MipMim_sum, &
                                       &this%MipMjp_sum, &
                                       &lat,fluct_val)
!	write(*,*) 'in measure_add,   MjpMim_sum = ' , this%MjpMim_sum
end subroutine


subroutine measure_eval(this,Cor_log, N_cell_in)
    class(exp_val),intent(inout)    :: this
    logical, intent(in)             :: Cor_log
    integer,intent(in)              :: N_cell_in
    real(8) :: av_site,av_Nadd,av_kt

    av_site=1.0d0/real(N_cell_in,8)
    av_Nadd=1.0d0/real(this%N_add,8)
    av_kt=1.0d0/this%kt

    this%C_av=(this%E_sq_sum_av*av_Nadd-(this%E_sum_av*av_Nadd)**2)*av_kt**2*av_site
    this%chi_M=(this%M_sq_sum_av(:)*av_Nadd-(this%M_sum_av(:)*av_Nadd)**2)*av_kt*av_site
    if (this%N_add.gt.1) this%E_err_av=sqrt(abs(this%E_sq_sum_av-(this%E_sum_av)**2*av_Nadd)/real(this%N_add-1,8))*av_site
    if (this%N_add.gt.1) this%M_err_av(:)=sqrt(abs(this%M_sq_sum_av(:)-this%M_sum_av(:)**2*av_Nadd)/real(this%N_add-1,8))*av_site
    this%E_av=this%E_sum_av*av_Nadd*av_site
    this%M_av=this%M_sum_av*av_Nadd*av_site
    this%qeulerp_av=this%qeulerp_av*av_Nadd/pi/4.0d0
    this%qeulerm_av=this%qeulerm_av*av_Nadd/pi/4.0d0

    this%chi_Q(1)=-((this%qeulerp_av-this%qeulerm_av)**2 - this%Q_sq_sum_av*av_Nadd/16.0d0/pi**2)*av_kT!/(this%qeulerp_av-this%qeulerm_av)! do in the postprocessing 
    this%chi_Q(2)=-(this%qeulerp_av**2-this%Qp_sq_sum_av*av_Nadd/16.0d0/pi**2)*av_kT!/this%qeulerp_av
    this%chi_Q(3)=-(this%qeulerm_av**2-this%Qm_sq_sum_av*av_Nadd/16.0d0/pi**2)*av_kT!/this%qeulerm_av
    this%chi_Q(4)=((this%Qm_sq_sum_av*av_Nadd/16.0d0/pi**2-this%qeulerm_av**2)* &
              &    (this%Qp_sq_sum_av*av_Nadd/16.0d0/pi**2-this%qeulerp_av**2))*av_kT!/(this%qeulerm_av*this%qeulerp_av)

    this%vortex_av=this%vortex_av*av_Nadd/3.0d0/sqrt(3.0d0)

	this%MjpMim_av=this%MjpMim_sum*av_Nadd*av_site
	this%MipMip_av=this%MipMip_sum*av_Nadd*av_site
	this%MipMim_av=this%MipMim_sum*av_Nadd*av_site
	this%MipMjp_av=this%MipMjp_sum*av_Nadd*av_site
	!write(*,*) "in measure_eval,   MjpMim_av = " , this%MjpMim_av

    if (Cor_log) this%chi_l=this%N_add !what is this supposed to do?
    Call print_av(this)
END subroutine

subroutine print_av(this)
    use, intrinsic :: iso_fortran_env, only : output_unit
    class(exp_val),intent(inout)    :: this

    write(output_unit,'(5(a,f18.9,2x))') 'M= ',norm2(this%M_av), &
      & 'E= ',this%E_av,'Q+= ',this%qeulerp_av,'Q-= ',-this%qeulerm_av,'Q= ',this%qeulerp_av-this%qeulerm_av
end subroutine

end module
