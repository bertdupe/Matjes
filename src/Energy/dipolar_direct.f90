module m_dipolar_direct
use m_input_H_types, only: io_H_dipole_direct
use m_derived_types, only: lattice
implicit none
private
public read_dip_input, get_dipolar

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
    type(io_H_dipole_direct),intent(out)        :: io

    Call get_parameter(io_param,fname,'mag_dip_direct',io%is_set) 
    Call get_parameter(io_param,fname,'mag_dip_period_cut',io%period_cutoff) 
!    Call get_parameter(io_param,fname,'mag_dip_dist_cut',io%dist_cutoff) 

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
    type(io_H_dipole_direct),intent(in) :: io   !input parameters 
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
        Nmag=lat%nmag
        N=Nmag*lat%Ncell

        !get positions
        Call lat%get_pos_mag(pos)
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
            write(*,*) 'DERP',i,N
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
        Ham%desc="dipolar direct"
        !set modes
        Call mode_set_rank1(Ham%mode_l,lat,"M")
        Call mode_set_rank1(Ham%mode_r,lat,"M")
    endif
end subroutine 

subroutine get_supercell_vec(supercell_vec,lat,period)
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
