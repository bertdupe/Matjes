module m_dipolar_fft
use m_input_H_types, only: io_H_dipole_direct
use m_derived_types, only: lattice
use m_fftw3
use m_dipolar_util, only: dip_pref, get_supercell_vec
implicit none



private
public read_dip_input, get_dipolar,dipolar_fft, get_dip

type        ::  dipolar_fft
    logical     ::  is_set=.false.
    type(c_ptr) :: plan_mag_F
!    type(c_ptr) :: plan_mag_I   !necessary?
    type(c_ptr) :: plan_H_I
    real(C_DOUBLE),allocatable              ::  M_n(:,:)
    complex(C_DOUBLE_COMPLEX),allocatable   ::  M_F(:,:)

    complex(C_DOUBLE_COMPLEX),allocatable   ::  K_F(:,:,:)

    complex(C_DOUBLE_COMPLEX),allocatable   ::  H_F(:,:) 
    real(C_DOUBLE),allocatable              ::  H_n(:,:) 
contains
    procedure   ::  set_M
    procedure   ::  get_H
    procedure   ::  get_E
    procedure   ::  get_E_distrib
end type


contains

subroutine set_M(this,lat)
    class(dipolar_fft),intent(inout)    ::  this
    type(lattice),intent(in)            ::  lat

    if(any(.not.lat%periodic.and.lat%dim_lat>1)) ERROR STOP "IMPLEMENT"
    this%M_n=lat%M%modes_v
end subroutine

subroutine get_H(this,lat,Hout)
    class(dipolar_fft),intent(inout)    ::  this
    type(lattice),intent(in)            ::  lat
    real(8),intent(inout)               ::  Hout(:,:)
    integer ::  i,j,l

    Call this%set_M(lat)
!    write(*,*) this%M_n
    Call fftw_execute_dft_r2c(this%plan_mag_F, this%M_n, this%M_F)

!    write(*,*) this%M_n
!    write(*,*) this%M_F
!    STOP

    this%H_F=cmplx(0.0d0,0.0d0,8)
    do j=1,size(this%M_F,2)
        do i=1,lat%Nmag*3
            do l=1,lat%Nmag*3
                this%H_F(i,j)=this%H_F(i,j)+this%K_F(i,l,j)*this%M_F(l,j)
            enddo
        enddo
    enddo
    Call fftw_execute_dft_c2r(this%plan_H_I, this%H_F, this%H_n)
    Hout=this%H_n
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

subroutine read_dip_input(io_param,fname,io)
    use m_io_utils
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_dipole_direct),intent(out)        :: io

    Call get_parameter(io_param,fname,'mag_dip_fft',io%is_set) 
    Call get_parameter(io_param,fname,'mag_dip_period_cut',io%period_cutoff) 
end subroutine

subroutine get_dip(dip,lat)
    type(dipolar_fft),allocatable,intent(inout) :: dip(:)
    type(lattice),intent(in)                    :: lat

    type(io_H_dipole_direct)    ::  io
    integer     :: io_unit

    open(newunit=io_unit,file='input')
    Call read_dip_input(io_unit,'input',io)
    close(io_unit)
    Call get_dipolar(dip,io,lat)
end subroutine

subroutine get_dipolar(dip,io,lat)
    use m_constants, only : pi
!$  use omp_lib    
    type(dipolar_fft),intent(inout),allocatable :: dip(:) 
    type(io_H_dipole_direct),intent(in)         :: io
    type(lattice),intent(in)                    :: lat

    integer         :: Ncell, dim_lat(3)
    integer         :: Nmag
    logical         :: period(3)
    integer         :: N_rep(3)
    integer(C_INT)  :: N_rep_rev(3)

    real(8),allocatable :: Karr(:,:,:)
    integer         ::  Nk_tot  !number of Dipole operator entries

    real(8),allocatable         :: supercell_vec(:,:)   !supercell lattice vectors whose periodicity are considered (3,:)
    integer         :: ind_zero

    real(8),allocatable         :: magmom(:)    !magnetic moment per magnetic atom within unit-cell

    integer         :: ii,i1,i2,i3, l, i
    integer         :: iat, jat
    integer         :: ibnd(2), jbnd(2)
    real(8)         :: pos_base_arr(3,3)
    real(8)         :: pos_base(3)
    real(8)         :: alat(3,3)
    real(8)         :: diff(3), diff_base(3)
    real(8)         :: diffunit(3), diffnorm
    real(8),parameter   ::  eps=1.0d-30

    real(8),allocatable,target  :: pos(:)       !position vector of all magnetic atoms (3*Nmag*Ncell)
    real(8),pointer             :: pos3(:,:)    !position vector of all magnetic atoms (3,Nmag*Ncell)

    real(8)         :: tmp_2D(3,3)
    real(8)         :: magmom_prod(3*lat%nmag,3*lat%nmag)
    real(8)         :: magmom_3(3*lat%nmag)

    integer         :: iarr(3)
    real(8)         :: iarr_real(3)
    integer         :: div(3),modul(3)


    integer(C_int)  :: howmany

    type(c_ptr) :: plan_K_F

    real(8),allocatable ::  Hout(:,:)


    if(io%is_set)then
        Nmag=lat%nmag
        Ncell=lat%Ncell
        period=lat%periodic
        dim_lat=lat%dim_lat
        alat=transpose(lat%areal)
        N_rep=dim_lat
        Nk_tot=product(N_rep)
        N_rep_rev=N_rep(size(N_rep):1:-1)
        if(any(.not.period.and.dim_lat>1)) ERROR STOP "ONLY PERIODIC CASE SO FAR"
!$      Call fftw_plan_with_nthreads(omp_get_max_threads())

        !get positions
        Call lat%get_pos_mag(pos)
        pos3(1:3,1:size(pos)/3)=>pos

        allocate(dip(1))


        !set magnetization plans
        allocate(dip(1)%M_N(3*Nmag,Nk_tot),source=0.0d0)
        allocate(dip(1)%M_F(3*Nmag,Nk_tot),source=cmplx(0.0d0,0.0d0,8))

        howmany=int(3*Nmag,C_int)
        dip(1)%plan_mag_F= fftw_plan_many_dft_r2c(int(3,C_INT), N_rep_rev, howmany,&
                                                 &dip(1)%M_n,   N_rep_rev,&
                                                 &howmany,      int(1,C_int), &
                                                 &dip(1)%M_F,   N_rep_rev,&
                                                 &howmany,      int(1,C_int), &
                                                 &FFTW_FORWARD+FFTW_MEASURE+FFTW_PATIENT)
!        dip(1)%plan_mag_I= fftw_plan_many_dft_c2r(int(3,C_INT), N_rep_rev, howmany,&
!                                                 &dip(1)%M_F,   N_rep_rev,&
!                                                 &howmany,      int(1,C_int), &
!                                                 &dip(1)%M_n,   N_rep_rev,&
!                                                 &howmany,      int(1,C_int), &
!                                                 &FFTW_BACKWARD+FFTW_MEASURE+FFTW_PATIENT)

!         Call dip(1)%set_M(lat)
!        write(*,*) "M PRE"
!        write(*,'(3F16.8)') dip(1)%M_n
!        write(*,*)
!
!        Call fftw_execute_dft_r2c(dip(1)%plan_mag_F, dip(1)%M_n, dip(1)%M_F)
!        dip(1)%M_n=0.0d0
!        Call fftw_execute_dft_c2r(dip(1)%plan_mag_I, dip(1)%M_F, dip(1)%M_n)
!        write(*,*) "M POST"
!        write(*,'(3F16.8)') dip(1)%M_n/real(product(N_rep_rev),8)


        !get supercell difference vectors
        Call get_supercell_vec(supercell_vec,lat,io%period_cutoff)
        ind_zero=minloc(norm2(supercell_vec,1),1)


        !prepare magnetic moment magnitudes
        Call lat%cell%get_mag_magmom(magmom)
        magmom_3=[(magmom((i-1)/3+1),i=1,3*Nmag)]   !unfold magnetic moment times 3 for all magnetization directions
        magmom_prod=spread(magmom_3,1,3*Nmag)*spread(magmom_3,2,3*Nmag) !get all products of magnetic moments in (3*Nmag,3*Nmag)-format for direct multiplication

        allocate(Karr(3*Nmag,3*Nmag,Nk_tot),source=0.0d0)
        
        !modul=[(product(dim_lat(:i)),i=1,3)]
        !div=[(product(dim_lat(:i-1)),i=1,3)]
        !do ii=1,Nk_tot
        !    iarr=modulo((ii-1),modul)/div
        !    iarr_real=real(iarr,8)
        !    pos_base=matmul(alat,iarr_real)
        !    do iat=1,Nmag
        !        ibnd=[(iat-1)*3+1,iat*3]
        !        do jat=1,Nmag
        !            jbnd=[(jat-1)*3+1,jat*3]
        !            diff_base=pos_base+pos3(:,iat)-pos3(:,iat)
        !            do l=1,size(supercell_vec,2)    !loop over different supercells
        !                diff=diff_base+supercell_vec(:,l)
        !                diffnorm=norm2(diff)
        !                diffunit=real(diff/cmplx(diffnorm,eps,8),8)
        !                tmp_2D(1,1)=3.0d0*diffunit(1)*diffunit(1)-1.0d0     !xx
        !                tmp_2D(2,1)=3.0d0*diffunit(2)*diffunit(1)           !yx 
        !                tmp_2D(3,1)=3.0d0*diffunit(3)*diffunit(1)           !zx 
        !                tmp_2D(1,2)=3.0d0*diffunit(1)*diffunit(2)           !xy 
        !                tmp_2D(2,2)=3.0d0*diffunit(2)*diffunit(2)-1.0d0     !yy
        !                tmp_2D(3,2)=3.0d0*diffunit(3)*diffunit(2)           !zy 
        !                tmp_2D(1,3)=3.0d0*diffunit(1)*diffunit(3)           !xz 
        !                tmp_2D(2,3)=3.0d0*diffunit(2)*diffunit(3)           !yz 
        !                tmp_2D(3,3)=3.0d0*diffunit(3)*diffunit(3)-1.0d0     !zz
        !                Karr(ibnd(1):ibnd(2),jbnd(1):jbnd(2),ii)=Karr(ibnd(1):ibnd(2),jbnd(1):jbnd(2),ii)+real(tmp_2D/(cmplx(diffnorm**3,eps,8)),8)
        !            enddo
        !        enddo
        !    enddo
        !enddo
        !old version without Nmag>3
        do i3=1,dim_lat(3)
            pos_base_arr(:,3)=alat(:,3)*real(i3-1,8)
            do i2=1,dim_lat(2)
                pos_base_arr(:,2)=alat(:,2)*real(i2-1,8)
                do i1=1,dim_lat(1)
                    pos_base_arr(:,1)=alat(:,1)*real(i1-1,8)
                    ii=i1+(i2-1)*dim_lat(1)+(i3-1)*dim_lat(1)*dim_lat(2)
                    pos_base=sum(pos_base_arr,2)
                    do l=1,size(supercell_vec,2)    !loop over different supercells
                        diff=pos_base+supercell_vec(:,l)
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
                        Karr(:,:,ii)=Karr(:,:,ii)+real(tmp_2D/(cmplx(diffnorm**3,eps,8)),8)
                    enddo
                    Karr(:,:,ii)=Karr(:,:,ii)*magmom_prod
                enddo
            enddo
        enddo
        Karr=-Karr*dip_pref*0.5d0*0.25d0/pi
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


        allocate(dip(1)%H_F(3*Nmag,Nk_tot),source=cmplx(0.0d0,0.0d0,8))
        allocate(dip(1)%H_n(3*Nmag,Nk_tot),source=0.0d0)
        howmany=int(3*Nmag,C_int)
        dip(1)%plan_H_I= fftw_plan_many_dft_c2r(int(3,C_INT), N_rep_rev, howmany,&
                                               &dip(1)%H_F,   N_rep_rev,&
                                               &howmany,      int(1,C_int), &
                                               &dip(1)%H_n,   N_rep_rev,&
                                               &howmany,      int(1,C_int), &
                                               &FFTW_BACKWARD+FFTW_MEASURE+FFTW_PATIENT)


        allocate(Hout,mold=dip(1)%M_n)
        Call dip(1)%get_H(lat,Hout)
        dip(1)%is_set=.true.
        if(Nmag/=1) ERROR STOP "FINISH NMAG>1" 
    endif
end subroutine 


!subroutine get_supercell_vec(supercell_vec,lat,period)
!    real(8),intent(inout),allocatable   :: supercell_vec(:,:)
!    type(lattice),intent(in)            :: lat
!    integer,intent(in)                  :: period(3)
!
!    integer             ::  bnd_ext(2,3)
!    integer             ::  Nrep(3)
!    integer             ::  Ntot
!    integer             ::  i,ii, i3,i2,i1
!    real(8)             ::  vec(3,3)
!
!    if(allocated(supercell_vec)) deallocate(supercell_vec)
!    bnd_ext=0
!    do i=1,3
!        if(lat%periodic(i)) bnd_ext(:,i)=[-period(i),period(i)]
!        Nrep(i)=bnd_ext(2,i)-bnd_ext(1,i)+1
!    enddo
!    Ntot=product(Nrep)
!    allocate(supercell_vec(3,Ntot))
!    ii=0
!    do i3=1,Nrep(3)
!        vec(:,3)=lat%a_sc(3,:)*real(bnd_ext(1,3)+i3-1,8)
!        do i2=1,Nrep(2)
!            vec(:,2)=lat%a_sc(2,:)*real(bnd_ext(1,2)+i2-1,8)
!            do i1=1,Nrep(1)
!                vec(:,1)=lat%a_sc(1,:)*real(bnd_ext(1,1)+i1-1,8)
!                ii=ii+1
!                supercell_vec(:,ii)=sum(vec,2)
!            enddo
!        enddo
!    enddo
!    if(ii/=Ntot) ERROR STOP "unexpected size for supercell vectors in dipolar calculation"
!end subroutine

end module
