module m_dipolar_magnetic
!module which contains the type to calculate the magnetic dipolar interaction using the discrete fourier-transformation of the FFTW3 library (and the slow explicit routine get_dipolar not described here anymore).

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
!\delta_{\alpha,\beta} is the kronecker delta between \alpha and \beta
!This corresponds to the conventional DDI Hamiltonian definitions, replacing the M_{\alpha} M_{\beta} with the compenents only.


use, intrinsic :: iso_fortran_env, only : output_unit
use m_input_H_types, only: io_H_dipole
use m_derived_types, only: lattice
use m_lattice_position, only: get_pos_mag
implicit none
private
public read_dip_input, get_dipolar, get_dipolar_fft

character(len=*),parameter  :: ham_desc="magnetic dipolar"

!Matjes length scale is nm, magnetic moments are in bohr magneton mu_b, energy shall be in eV
!mu_0*mu_b^2/(nm^3)  * J_to_eV
!mu_0=1.25663706212d-6  kg*m/s^2/A^2
!mu_b=9.2740100783d-24  A*m^2
!nm=1.0d-9  m
!J_to_eV=1/1.602176634d-19 
real(8),parameter   ::  dip_pref=6.745817653234234066975d-7

contains

subroutine read_dip_input(io_param,fname,io)
    use m_io_utils
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_dipole),intent(out)   :: io

    Call get_parameter(io_param,fname,'mag_dip_use',io%is_set) 
    if(.not.io%is_set) return
    Call get_parameter(io_param,fname,'mag_dip_period_cut',io%period_cutoff) 
    io%fft=.true.
    Call get_parameter(io_param,fname,'mag_dip_fft',io%fft) 
end subroutine

subroutine get_dipolar_fft(dip,io,lat)
    !main routine to describe the dipolar interaction using the discrete fourier transformation
    use m_constants, only : pi
    use m_fft_H_base, only: fft_H
!$  use omp_lib    
    class(fft_H),intent(inout)              :: dip
    type(io_H_dipole),intent(in)            :: io
    type(lattice),intent(in)                :: lat

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

    integer         ::  Kbd(2,3)    !boundaries of the K-operator

    real(8)         :: tmp_2D(3,3)  !temporary values for a pair set of magnetic atoms
    integer         :: i            !multi-purpose loop index

    real(8)         :: iarr_real(3) !real vector version for i1,i2,i3, which gives difference vector due to supercell difference
    integer         :: i_mult(3)    !convenience multiplier to get i_k from i_add
    integer         :: i_add(3) !position in N_rep(:) array for considered entry

    if(io%is_set)then
        write(output_unit,'(/2A)') "Start setting fft-Hamiltonian: ", ham_desc
        !set some initial parameters locally for convencience
        Nmag=lat%nmag
        period=lat%periodic.or.lat%dim_lat==1
        alat=transpose(lat%areal)

        !set shape-dependent quantities of fft_H and get Kdb,N_rep
        Call dip%init_shape(3*lat%nmag,period,lat%dim_lat,Kbd,N_rep)
        Nk_tot=product(N_rep)

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
        Call dip%init_op(3*Nmag,Karr,ham_desc)
    endif
end subroutine 



subroutine get_dipolar(Ham,io,lat)
    !poor mans approach to get dipolar energy manually by an in principle dense matrix
    !the implementation is rather studid since it uses a sparse matrix format to save a dense matrix (also it doesn't check for zeros...)
    !this scales really terribly to large systems and should only be used for very small systems or as a check for the fft dipolar inplementation
    use m_H_public
    use m_mode_public
    use m_setH_util, only: get_coo
    use m_neighbor_type, only: neighbors
    use m_constants, only : pi

    class(t_H),intent(inout)            :: Ham  !Hamiltonian in which all contributions are added up
    type(io_H_dipole),intent(in)        :: io   !input parameters 
    type(lattice),intent(in)            :: lat

    real(8),allocatable ::  val(:)              !value array of coo-matrix
    integer,allocatable ::  rowind(:),colind(:) !indice arrays of coo-matrix
    integer             ::  nnz                 !number of entries in coo-arrays

    real(8),allocatable,target  :: pos(:)       !position vector of all magnetic atoms (3*Nmag*Ncell)
    real(8),pointer             :: pos3(:,:)    !position vector of all magnetic atoms (3,Nmag*Ncell)
    real(8) :: diff(3), diffunit(3), diffnorm    !difference vector, unit-vector, and norm

    integer                     :: Nmag         !number of magnetic moments
    real(8),allocatable         :: magmom(:)    !magnetic moment per magnetic atom within unit-cell
    real(8)                     :: mag_i,mag_j  !local magnetic moment in loop

    real(8),allocatable         :: supercell_vec(:,:)   !supercell lattice vectors whose periodicity are considered (3,:)
    integer                     :: ind_zero !index where supercell_vec contains the 0 entry (0.0,0.0,0.0)

    real(8)     :: tmp_val(9)   !temporary Hamiltonian values for a i-j pair
    integer     :: i, j, l, ii  !some loop parameters
    integer     :: N            !number of total magnetic atoms in entire supercell


    if(io%is_set)then
        write(output_unit,'(/2A)') "Start setting Hamiltonian: ", ham_desc
        Nmag=lat%nmag
        N=Nmag*lat%Ncell

        !get positions
        Call get_pos_mag(lat,pos)
        pos3(1:3,1:size(pos)/3)=>pos

        !get supercell difference vectors
        Call get_supercell_vec(supercell_vec,lat,io%period_cutoff)
        ind_zero=minloc(norm2(supercell_vec,1),1)

        !get magnetic moment magnitudes
        Call lat%cell%get_mag_magmom(magmom)

        !prepare coo-matrix
        nnz=N*N*9
        allocate(val(nnz),source=0.0d0)
        allocate(rowind(nnz),colind(nnz),source=0)


        !fill all Hamiltonian entries
        ii=1
        do i=1,N
            mag_i=magmom(modulo(i-1,Nmag)+1)
            do j=1,N
                mag_j=magmom(modulo(j-1,Nmag)+1)
                rowind(ii:ii+8)=(i-1)*3+[1,1,1,2,2,2,3,3,3] 
                colind(ii:ii+8)=(j-1)*3+[1,2,3,1,2,3,1,2,3] 
                do l=1,size(supercell_vec,2)    !loop over different supercells
                    if(j==i.and.l==ind_zero) cycle  !no entry for really same site
                    diff=pos3(:,j)-pos3(:,i)+supercell_vec(:,l)
                    diffnorm=norm2(diff)
                    diffunit=diff/diffnorm
                    tmp_val(1)=3.0d0*diffunit(1)*diffunit(1)-1.0d0     !xx
                    tmp_val(2)=3.0d0*diffunit(2)*diffunit(1)           !yx 
                    tmp_val(3)=3.0d0*diffunit(3)*diffunit(1)           !zx 
                    tmp_val(4)=3.0d0*diffunit(1)*diffunit(2)           !xy 
                    tmp_val(5)=3.0d0*diffunit(2)*diffunit(2)-1.0d0     !yy
                    tmp_val(6)=3.0d0*diffunit(3)*diffunit(2)           !zy 
                    tmp_val(7)=3.0d0*diffunit(1)*diffunit(3)           !xz 
                    tmp_val(8)=3.0d0*diffunit(2)*diffunit(3)           !yz 
                    tmp_val(9)=3.0d0*diffunit(3)*diffunit(3)-1.0d0     !zz
                    val(ii:ii+8)=val(ii:ii+8)+tmp_val/(diffnorm**3)*mag_j*mag_i
                enddo
                ii=ii+9
            enddo
        enddo
        val=-val*dip_pref*0.5d0*0.25d0/pi

        !initialize Hamiltonian array with calculated parameters
        Call Ham%init_coo(rowind,colind,val,[Nmag*3,Nmag*3],"M","M",lat,2)
        Ham%desc=ham_desc
        !set modes
        Call mode_set_rank1(Ham%mode_l,lat,"M")
        Call mode_set_rank1(Ham%mode_r,lat,"M")
    endif
end subroutine 

subroutine get_supercell_vec(supercell_vec,lat,period)
    use m_derived_types, only: lattice
    real(8),intent(inout),allocatable   :: supercell_vec(:,:)
    type(lattice),intent(in)            :: lat
    integer,intent(in)                  :: period(3)

    integer             ::  bnd_ext(2,3)
    integer             ::  Nrep(3)
    integer             ::  Ntot
    integer             ::  i,ii, i3,i2,i1
    real(8)             ::  vec(3,3)

    if(allocated(supercell_vec)) deallocate(supercell_vec)
    bnd_ext=0
    do i=1,3
        if(lat%periodic(i)) bnd_ext(:,i)=[-period(i),period(i)]
        Nrep(i)=bnd_ext(2,i)-bnd_ext(1,i)+1
    enddo
    Ntot=product(Nrep)
    allocate(supercell_vec(3,Ntot))
    ii=0
    do i3=1,Nrep(3)
        vec(:,3)=lat%a_sc(3,:)*real(bnd_ext(1,3)+i3-1,8)
        do i2=1,Nrep(2)
            vec(:,2)=lat%a_sc(2,:)*real(bnd_ext(1,2)+i2-1,8)
            do i1=1,Nrep(1)
                vec(:,1)=lat%a_sc(1,:)*real(bnd_ext(1,1)+i1-1,8)
                ii=ii+1
                supercell_vec(:,ii)=sum(vec,2)
            enddo
        enddo
    enddo
    if(ii/=Ntot) ERROR STOP "unexpected size for supercell vectors in dipolar calculation"
end subroutine


end module
