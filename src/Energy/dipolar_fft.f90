module m_dipolar_fft
!module which contains the type to calculate the magnetic dipolar interaction using the discrete fourier-transformation of the FFTW3 library.

!The main procedure is implemented along: N. Hayashi et al., Jpn. J. Appl. Phys. 35 6065 (1996)
!The basic idea is that the effective field H of the dipolar interaction (DDI) is a convolution between the magnetic moments (M) and the demagnetization tensor (K).
!Hence the effective field can be calculated as the inverse fourier transform of the product between the fourier-transformed magnetic moments and the fourier-transformed demagnetization tensor.
!The FT of K has to be calculated only once.

!In case of periodic boundary conditions that is relatively straight-forwards, only the contributions of the periodic images of K have to be summed up to a certain cutoff.

!In the case with open boundaries the case is a bit more difficult and the dimensions of M and K have to be doubled to cover all interactions by the discrete convolution.
!However, I seemingly failed to understand how exactly the enlarged K has to be constructed and implemented something of which I am not sure if it is what is described in the paper, however it seems to work and makes sense to me.
!if K(i) corresponds to the interaction between the M at the 0-site and at the -i-site ((-N+1,...,-1,0,1,...N-1)-sites), then I implemented the K-array as
![K(0),K(1),...,K(N-1),0,K(-N+1),...,K(-1)], instead of what is in the paper which I am not sure I correctly understand
!the corresponding M array is [M(0),M(1),...,M(N-1),0,...,0]
!notice that K(1) is actually a (3*Nmag,3*Nmag) tensor and M(1) is [Mx,My,Mz], so for H there also has to be a matrix.vector product.
!The FT in of each component is done using the howmany-feature from FFTW3.

!The K (3*Nmag,3*Nmag) components are expressed as:  K(\alpha,\beta)= 1/2*\mu_0 /(4*\pi) \sum_{ij,i/=j} \mu_i \mu_j \frac{3 f({\alpha,\beta}) - \delta_{\alpha,\beta}}{r^3}
!where \mu_0 is the vacuum permeability, \mu_i is the magetic moment of the atom_i, r is the distance between site i and j, f{\alpha,\beta} gives the normalized position difference product of the \alpha and \beta component, and 
!\delta_{\alpha,\beta} is the kroneker delta between \alpha and \beta
!This corresponds to the conventional DDI Hamiltonian definitions, replacing the M_{\alpha} M_{\beta} with the compenents only.

use m_input_H_types, only: io_H_dipole_direct
use m_derived_types, only: lattice
use m_fftw3
use m_dipolar_util, only: dip_pref, get_supercell_vec
use m_dipolar_fft_internal, only: int_set_M, int_get_H
implicit none

private
public dipolar_fft, get_dip

type        ::  dipolar_fft
    logical         :: set=.false.  !all parameters have been initialized
    integer         :: N_rep(3)=0   !size in each dimension of each considered field
    type(c_ptr)     :: plan_mag_F   !FFTW plan M_n -> M_F
    type(c_ptr)     :: plan_H_I     !FFTW plan H_F -> H_n
    real(C_DOUBLE),allocatable              ::  M_n(:,:)    !magnetization in normal-space
    complex(C_DOUBLE_COMPLEX),allocatable   ::  M_F(:,:)    !magnetization in fourier-space

    complex(C_DOUBLE_COMPLEX),allocatable   ::  K_F(:,:,:)  !demagnetization tensor in fourier-space

    complex(C_DOUBLE_COMPLEX),allocatable   ::  H_F(:,:)    !effective field fourier-space 
    real(C_DOUBLE),allocatable              ::  H_n(:,:)    !effective field normal-space
    procedure(int_set_M), pointer,nopass    ::  M_internal => null()    !function to set internal magnetization depending on periodic/open boundaries
    procedure(int_get_H), pointer,nopass    ::  H_internal => null()    !function to get effective field depending on periodic/open boundaries
contains
    procedure           :: set_M            !set internal magnetization in normal-space from lattice
    procedure           :: init_internal    !initialize internal procedures

    procedure,public    :: get_H            !get effective field
    procedure,public    :: get_E            !get energy
    procedure,public    :: get_E_distrib    !get energy-distribution in Nmag*Ncell-space
    procedure,public    :: is_set           !returns set

end type


contains

pure function is_set(this)result(set)
    class(dipolar_fft),intent(in)       :: this
    logical     ::  set
    set=this%set
end function

subroutine init_internal(this,periodic)
    !initialize internal procedures M_internal and H_internal
    use m_dipolar_fft_internal
    class(dipolar_fft),intent(inout)    :: this
    logical,intent(in)                  :: periodic(3)  !T: periodic boundary, F: open boundary

    if(all(periodic))then
        this%M_internal=>set_M_period_TTT
        this%H_internal=>set_H_period_TTT
    elseif(all(periodic(1:2)))then
        this%M_internal=>set_M_period_TTF
        this%H_internal=>set_H_period_TTF
    elseif(all(periodic(1:1)))then
        this%M_internal=>set_M_period_TFF
        this%H_internal=>set_H_period_TFF
    else
        this%M_internal=>set_M_period_FFF
        this%H_internal=>set_H_period_FFF
    endif
end subroutine

subroutine set_M(this,lat)
    !set this%M_n from lat%M%modes according to this%M_internal
    class(dipolar_fft),intent(inout)    ::  this
    type(lattice),intent(in)            ::  lat

    Call this%M_internal(this%M_n,lat%M%modes,lat%dim_lat,this%N_rep,size(this%M_n,1))
end subroutine

subroutine get_H(this,lat,Hout)
    !get effective field H 
    class(dipolar_fft),intent(inout)    ::  this
    type(lattice),intent(in)            ::  lat
    real(8),intent(inout)               ::  Hout(:,:)
#ifdef CPP_FFTW3
    integer ::  i,j,l

    Call this%set_M(lat)
    Call fftw_execute_dft_r2c(this%plan_mag_F, this%M_n, this%M_F)

    this%H_F=cmplx(0.0d0,0.0d0,8)
    do j=1,size(this%M_F,2)
        do i=1,lat%Nmag*3
            do l=1,lat%Nmag*3
                this%H_F(i,j)=this%H_F(i,j)+this%K_F(i,l,j)*this%M_F(l,j)
            enddo
        enddo
    enddo
    Call fftw_execute_dft_c2r(this%plan_H_I, this%H_F, this%H_n)
    Call this%H_internal(this%H_n,Hout,lat%dim_lat,this%N_rep,size(Hout,1))
#else
    ERROR STOP "dipolar_fft%get_H requires CPP_FFTW"
#endif
end subroutine

subroutine get_E_distrib(this,lat,Htmp,E)
    class(dipolar_fft),intent(inout)    ::  this
    type(lattice),intent(in)            ::  lat
    real(8),intent(inout)               ::  Htmp(:,:)
    real(8),intent(out)                 ::  E(:)

    Call this%get_H(lat,Htmp)
    Htmp=Htmp*lat%M%modes_v
    E=sum(reshape(Htmp,[3,lat%Nmag*lat%Ncell]),1)*2.0d0 !not sure about *2.0d0
end subroutine


subroutine get_E(this,lat,Htmp,E)
    class(dipolar_fft),intent(inout)    ::  this
    type(lattice),intent(in)            ::  lat
    real(8),intent(inout)               ::  Htmp(:,:)
    real(8),intent(out)                 ::  E

    Call this%get_H(lat,Htmp)
    Htmp=Htmp*lat%M%modes_v
    E=sum(Htmp)
end subroutine

subroutine read_dip_input(io_unit,fname,io)
    use m_io_utils
    integer,intent(in)                      :: io_unit
    character(len=*), intent(in)            :: fname
    type(io_H_dipole_direct),intent(out)    :: io

    Call get_parameter(io_unit,fname,'mag_dip_fft',io%is_set) 
    Call get_parameter(io_unit,fname,'mag_dip_period_cut',io%period_cutoff) 
end subroutine

subroutine get_dip(dip,lat)
    type(dipolar_fft),allocatable,intent(inout) :: dip(:)
    type(lattice),intent(in)                    :: lat

    type(io_H_dipole_direct)    :: io
    integer                     :: io_unit

    open(newunit=io_unit,file='input')
    Call read_dip_input(io_unit,'input',io)
    close(io_unit)
    Call get_dipolar(dip,io,lat)
end subroutine

subroutine get_dipolar(dip,io,lat)
    !main routine to initialize the diplar_fft type
    use m_constants, only : pi
!$  use omp_lib    
    type(dipolar_fft),intent(inout),allocatable :: dip(:) 
    type(io_H_dipole_direct),intent(in)         :: io
    type(lattice),intent(in)                    :: lat

    integer         :: Ncell            !number of unit-cells (product of dim_lat)
    integer         :: dim_lat(3)       !number of unit-cells in each direction in supercell
    integer         :: Nmag             !number of magnetic atoms per unit-cell
    logical         :: period(3)        !consider as periodic or open boundary condition along each direction (T:period, F:open)
                                        ! (dim_lat(i)=1->period(i)=T, since the calculation in the periodic case is easier, but choice of supercell_vec still does not consider periodicity)
    integer         :: N_rep(3)         !number of states in each direction in the fourier transformation
    integer         :: Nk_tot           !number of state considered in FT (product of N_rep)

    real(8),allocatable :: Karr(:,:,:)  !K-operator to be FT'd (1:3*Nmag,1:3*Nmag,1:Nk_tot)

    real(8),allocatable :: supercell_vec(:,:)   !supercell lattice vectors whose periodicity are considered (3,:)
    integer             :: i_per                !index to loop over supercell_vec(1:3,*)

    real(8),allocatable :: magmom(:)    !magnetic moment per magnetic atom within unit-cell
    real(8)             :: magmom_prod(3*lat%nmag,3*lat%nmag)  !product of magmom_3 for all atom combinations 
    real(8)             :: magmom_3(3*lat%nmag) !magnetic moments each repeated thrice

    integer         :: i_k              !index for Karr(:,:,*) (encodes flattened N_rep(1:3))
    integer         :: i1,i2,i3         !indices specifying difference in unit-cell space
    integer         :: iat, jat         !intizes for atoms
    integer         :: ibnd(2), jbnd(2) !boundaries in K(*,*,:) for the respective atoms
    real(8)         :: pos_base(3)      !difference vector between 2 sites without atomic position or periodic-image
    real(8)         :: alat(3,3)        !real-space lattice
    real(8)         :: diff_base(3)     !difference vector between 2 sites without periodic image
    real(8)         :: diff(3)          !difference vector between 2 sites
    real(8)         :: diffnorm         !norm of difference vector
    real(8)         :: diffunit(3)      !normalized difference vector 
    real(8),parameter   ::  eps=1.0d-30 !small value to avoid division by zero

    integer,allocatable     :: ind_mag(:)   !cell atoms with magnetic moment indices
    real(8)                 :: pos_mag(1:3,lat%nmag)       !position vector of all magnetic atoms (3*Nmag*Ncell)
!    real(8),pointer             :: pos3(:,:)    !position vector of all magnetic atoms (3,Nmag*Ncell)

    integer         ::  Kbd(2,3)    !boundaries of the K-operator

    real(8)         :: tmp_2D(3,3)  !temporary values for a pair set of magnetic atoms
    integer         :: i            !multi-purpose loop index

    real(8)         :: iarr_real(3) !real vector version for i1,i2,i3, which gives difference vector due to supercell difference
    integer         :: i_mult(3)    !convenience multiplier to get i_k from i_add
    integer         :: i_add(3) !position in N_rep(:) array for considered entry

    !FFTW util
    integer(C_INT)  :: N_rep_rev(3)     !reversed N_rep necessary for fftw3 (col-major -> row-major)
    integer(C_int)  :: howmany  !dimension of quantitiy which is fourier-transformed (see FFTW3)
    type(c_ptr)     :: plan_K_F !plan for fourier transformation of K

    if(io%is_set)then
#ifdef CPP_FFTW3
        !set some initial parameters locally for convencience
        Nmag=lat%nmag
        Ncell=lat%Ncell
        dim_lat=lat%dim_lat
        period=lat%periodic.or.dim_lat==1
        alat=transpose(lat%areal)
        N_rep=dim_lat
        !set K-boundaries for periodic boundaries
        Kbd(1,:)=[0,0,0]
        Kbd(2,:)=dim_lat-1
        !set K-boundaries for open boundaries if open
        do i=1,3
            if(.not.period(i))then
                N_rep(i)=2*N_rep(i)
                Kbd(:,i)=[-dim_lat(i)+1,dim_lat(i)-1]
            endif
        enddo
        Nk_tot=product(N_rep)
        N_rep_rev=N_rep(size(N_rep):1:-1)
        allocate(dip(1))
        dip(1)%N_rep=N_rep
        Call dip(1)%init_internal(period)

!$      Call fftw_plan_with_nthreads(omp_get_max_threads())
        !allocate space for magnetic moment in real and fourier space and prepare plan for fourier-transformation
        allocate(dip(1)%M_N(3*Nmag,Nk_tot),source=0.0d0)
        allocate(dip(1)%M_F(3*Nmag,Nk_tot),source=cmplx(0.0d0,0.0d0,8))
        howmany=int(3*Nmag,C_int)
        dip(1)%plan_mag_F= fftw_plan_many_dft_r2c(int(3,C_INT), N_rep_rev, howmany,&
                                                 &dip(1)%M_n,   N_rep_rev,&
                                                 &howmany,      int(1,C_int), &
                                                 &dip(1)%M_F,   N_rep_rev,&
                                                 &howmany,      int(1,C_int), &
                                                 &FFTW_FORWARD+FFTW_MEASURE+FFTW_PATIENT)

        !get positions of magnetic atoms in unit_cell
        Call lat%cell%ind_mag_all(ind_mag)
        do i=1,Nmag
            pos_mag(:,i)=lat%cell%atomic(ind_mag(i))%position
        enddo
        deallocate(ind_mag)

        !get supercell difference vectors
        Call get_supercell_vec(supercell_vec,lat,io%period_cutoff)

        !prepare K-matrix in real-space
        allocate(Karr(3*Nmag,3*Nmag,Nk_tot),source=0.0d0)
        !get position dependent parts of K-tensor
        i_mult=[(product(N_rep(:i-1)),i=1,3)]
        do i3=Kbd(1,3),Kbd(2,3)
            i_add(3)=(i3-floor(real(i3,8)/N_rep(3))*N_rep(3))*i_mult(3)+1
            do i2=Kbd(1,2),Kbd(2,2)
                i_add(2)=(i2-floor(real(i2,8)/N_rep(2))*N_rep(2))*i_mult(2)
                do i1=Kbd(1,1),Kbd(2,1)
                    i_add(1)=(i1-floor(real(i1,8)/N_rep(1))*N_rep(1))*i_mult(1)
                    i_k=sum(i_add)
                    iarr_real=[real(i1,8),real(i2,8),real(i3,8)]
                    pos_base=matmul(alat,iarr_real)
                    do iat=1,Nmag
                        ibnd=[(iat-1)*3+1,iat*3]
                        do jat=1,Nmag
                            jbnd=[(jat-1)*3+1,jat*3]
                            diff_base=pos_base+pos_mag(:,iat)-pos_mag(:,jat)
                            do i_per=1,size(supercell_vec,2)    !loop over different supercells
                                diff=diff_base+supercell_vec(:,i_per)
                                diffnorm=norm2(diff)
                                diffunit=real(diff/cmplx(diffnorm,eps,8),8)
                                tmp_2D(1,1)=3.0d0*diffunit(1)*diffunit(1)-1.0d0     !xx
                                tmp_2D(2,1)=3.0d0*diffunit(2)*diffunit(1)           !yx 
                                tmp_2D(3,1)=3.0d0*diffunit(3)*diffunit(1)           !zx 
                                tmp_2D(1,2)=3.0d0*diffunit(1)*diffunit(2)           !xy 
                                tmp_2D(2,2)=3.0d0*diffunit(2)*diffunit(2)-1.0d0     !yy
                                tmp_2D(3,2)=3.0d0*diffunit(3)*diffunit(2)           !zy 
                                tmp_2D(1,3)=3.0d0*diffunit(1)*diffunit(3)           !xz 
                                tmp_2D(2,3)=3.0d0*diffunit(2)*diffunit(3)           !yz 
                                tmp_2D(3,3)=3.0d0*diffunit(3)*diffunit(3)-1.0d0     !zz
                                Karr(ibnd(1):ibnd(2),jbnd(1):jbnd(2),i_k)=Karr(ibnd(1):ibnd(2),jbnd(1):jbnd(2),i_k)+real(tmp_2D/(cmplx(diffnorm**3,eps,8)),8)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        !multiply magnetic moment magnitude product to K-tensor
        Call lat%cell%get_mag_magmom(magmom)
        magmom_3=[(magmom((i-1)/3+1),i=1,3*Nmag)]   !unfold magnetic moment times 3 for all magnetization directions
        magmom_prod=spread(magmom_3,1,3*Nmag)*spread(magmom_3,2,3*Nmag) !get all products of magnetic moments in (3*Nmag,3*Nmag)-format for direct multiplication
        do i_k=1,Nk_tot
            Karr(:,:,i_k)=Karr(:,:,i_k)*magmom_prod
        enddo
        !get correct prefractor for K
        Karr=-Karr*dip_pref*0.5d0*0.25d0/pi

        !calculate fourier transform of K and save it in dipole-type 
        howmany=int(9*Nmag**2,C_int)
        allocate(dip(1)%K_F(3*Nmag,3*Nmag,Nk_tot),source=cmplx(0.0d0,0.0d0,8))
        plan_K_F= fftw_plan_many_dft_r2c(int(3,C_INT), N_rep_rev,  howmany,&
                                        &Karr,         N_rep_rev,&
                                        &howmany,      int(1,C_int),&
                                        &dip(1)%K_F,   N_rep_rev,&
                                        &howmany,      int(1,C_int),&
                                        &FFTW_FORWARD+FFTW_ESTIMATE)
        Call fftw_execute_dft_r2c(plan_K_F, Karr, dip(1)%K_F)
        dip(1)%K_F=dip(1)%K_F/real(product(N_rep_rev),8)
        Call fftw_destroy_plan(plan_K_F)

        !allocate space for effective field and initialize plan for inverse fourier transformation there 
        allocate(dip(1)%H_F(3*Nmag,Nk_tot),source=cmplx(0.0d0,0.0d0,8))
        allocate(dip(1)%H_n(3*Nmag,Nk_tot),source=0.0d0)
        howmany=int(3*Nmag,C_int)
        dip(1)%plan_H_I= fftw_plan_many_dft_c2r(int(3,C_INT), N_rep_rev, howmany,&
                                               &dip(1)%H_F,   N_rep_rev,&
                                               &howmany,      int(1,C_int), &
                                               &dip(1)%H_n,   N_rep_rev,&
                                               &howmany,      int(1,C_int), &
                                               &FFTW_BACKWARD+FFTW_MEASURE+FFTW_PATIENT)

        dip(1)%set=.true.
#else
        ERROR STOP "CANNOT USE FFTW DIPOL (mag_dip_fft) without FFTW (CPP_FFTW3)"
#endif
    endif
end subroutine 


end module
