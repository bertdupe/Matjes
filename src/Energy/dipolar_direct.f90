module m_dipolar_direct
use m_input_H_types, only: io_H_dipole_direct
use m_derived_types, only: lattice
use m_dipolar_util, only: dip_pref, get_supercell_vec
implicit none
private
public read_dip_input, get_dipolar
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
end module
