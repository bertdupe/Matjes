module m_dipolar_direct
use m_input_H_types, only: io_H_dipole_direct
use m_derived_types, only: lattice
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
    Call get_parameter(io_param,fname,'mag_dip_dist_cut',io%dist_cutoff) 

end subroutine

subroutine get_dipolar(Ham,io,lat)
    !get coupling in t_H Hamiltonian format
    use m_H_public
    use m_mode_public
    use m_setH_util, only: get_coo
    use m_neighbor_type, only: neighbors

    class(t_H),intent(inout)            :: Ham  !Hamiltonian in which all contributions are added up
    type(io_H_dipole_direct),intent(in) :: io
    type(lattice),intent(in)            :: lat


!    !local Hamiltonian
!    real(8),allocatable  :: Htmp(:,:)   !local Hamiltonian in (dimmode(1),dimmode(2))-basis
!    !local Hamiltonian in coo format
!    real(8),allocatable  :: val_tmp(:)
!    integer,allocatable  :: ind_tmp(:,:)
!
    real(8),allocatable ::  diff_vec_arr(:,:,:)
    real(8),allocatable ::  one_over_r3(:,:)

    integer     :: i,j,ii
    integer     ::  N
!
!    integer         :: i_atpair,N_atpair    !loop parameters which atom-type connection are considered (different neighbor types)
!    integer         :: i_dist,N_dist        !loop parameters which  connection are considered (different neighbor types)
!    integer         :: i_pair           !loop keeping track which unique connection between the same atom types is considered (indexes "number shells" in neighbors-type)
!    integer         :: i_shell          !counting the number of unique connection for given atom types and a distance
!    integer         :: connect_bnd(2)   !indices keeping track of which pairs are used for the particular connection
!    type(neighbors) :: neigh            !all neighbor information for a given atom-type pair
!    real(8)         :: J                !magnitude of Hamiltonian parameter
!    integer         :: atind_mag(2)     !index of considered atom in basis of magnetic atoms (1:Nmag)
    real(8),allocatable ::  val(:)
    integer,allocatable ::  rowind(:),colind(:)
!
    if(io%is_set)then
        Call get_diff_vec_arr(diff_vec_arr,lat)
        N=size(diff_vec_arr,2)
        allocate(one_over_r3(N.N))
        Call get_1_over_R3(diff_vec_arr,size(one_over_r3),one_over_r3)
        do i=1,size(one_over_r3,1)
            one_over_r3(i,i)=0.0d0
        enddo
        allocate(val((N-1)*N)
        ii=0
        do i=1,N
            do j=1,N
                if(j==i) cycle
                ii=ii+1
                rowind(ii)=j
                colind(ii)=j
                val

            enddo
        enddo
        

        write(*,*) minval(one_over_r3),maxval(one_over_r3)


        
!        Call get_Htype(Ham_tmp)
!        allocate(Htmp(lat%M%dim_mode,lat%M%dim_mode))!local Hamiltonian modified for each shell/neighbor
!        do i_atpair=1,N_atpair
!            !loop over different connected atom types
!            Call neigh%get(io%pair(i_atpair)%attype,io%pair(i_atpair)%dist,lat)
!            N_dist=size(io%pair(i_atpair)%dist)
!            i_pair=0
!            connect_bnd=1 !initialization for lower bound
!            do i_dist=1,N_dist
!                !loop over distances (nearest, next nearest,... neighbor)
!                J=io%pair(i_atpair)%val(i_dist)
!                do i_shell=1,neigh%Nshell(i_dist)
!                    !loop over all different connections with the same distance
!                    i_pair=i_pair+1
!
!                    !set local Hamiltonian in basis of magnetic orderparameter
!                    atind_mag(1)=lat%cell%ind_mag(neigh%at_pair(1,i_pair))
!                    atind_mag(2)=lat%cell%ind_mag(neigh%at_pair(2,i_pair))
!                    Htmp=0.0d0
!                    Htmp(atind_mag(1)*3-2,atind_mag(2)*3-2)=J
!                    Htmp(atind_mag(1)*3-1,atind_mag(2)*3-1)=J
!                    Htmp(atind_mag(1)*3  ,atind_mag(2)*3  )=J
!                    connect_bnd(2)=neigh%ishell(i_pair)
!                    Htmp=-Htmp !flip sign corresponding to previous implementation
!                    Call get_coo(Htmp,val_tmp,ind_tmp)
!
!                    !fill Hamiltonian type
!                    Call Ham_tmp%init_connect(neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val_tmp,ind_tmp,"MM",lat,2)
!                    deallocate(val_tmp,ind_tmp)
!                    Call Ham%add(Ham_tmp)
!                    Call Ham_tmp%destroy()
!                    connect_bnd(1)=connect_bnd(2)+1
!                enddo 
!            enddo
!        enddo
        Ham%desc="dipolar direct"
        !set modes
        Call mode_set_rank1(Ham%mode_l,lat,"M")
        Call mode_set_rank1(Ham%mode_r,lat,"M")
    endif

    ERROR STOP "FINISH IMPLEMENT DIPOLAR DIRECT"

end subroutine 

subroutine get_1_over_R3(dist,Ndist,div)
    integer,intent(in)      :: Ndist
    real(8),intent(in)      :: dist(3,Ndist)
    real(8),intent(out)     :: div(Ndist)

    div=norm2(dist,dim=1)
    div=div**3
    div=1.0d0/div
end subroutine

subroutine get_diff_vec_arr(diff_vec,lat)
    real(8),intent(inout),allocatable   :: diff_vec(:,:,:)
    type(lattice),intent(in)            :: lat
    integer     ::  N
    real(8),allocatable,target  ::  pos(:)
    real(8),pointer             :: pos3(:,:)
    integer ::  i,j

    if(allocated(diff_vec)) deallocate(diff_vec)
    N=lat%Ncell*lat%nmag
    allocate(diff_vec(3,N,N),source=0.0d0)

    Call lat%get_pos_mag(pos)
    pos3(1:3,1:N)=>pos

!    do i=1,N
!        diff_vec(:,:,i)=pos3-spread(pos3(:,i),2,N)
!    enddo
    do i=1,N
        do j=1,N
            diff_vec(:,j,i)=pos3(:,j)-pos3(:,i)
        enddo
    enddo
    
    nullify(pos3)
    if(allocated(pos)) deallocate(pos)
end subroutine
end module
