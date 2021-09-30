module m_fluct
use m_derived_types,only: lattice
!use m_get_position
use m_neighbor_type
implicit none
private
public fluct_parameters,init_fluct_parameter,eval_fluct, eval_fluct_spatial, print_fluct_spatial

type fluct_parameters
    logical                 :: l_use=.False. 
    type(neighbors)         :: neigh
    integer                 :: ind(2)  !index of positive and negative nearest neighbor in neigh%diff_vec array 
contains
    procedure get_nneigh
end type
contains

function get_nneigh(this)result(nneigh)
    class(fluct_parameters),intent(in)  :: this
    integer                             :: nneigh

    nneigh=this%neigh%Nshell(2) !2 because one is on-site, 2 next nearest neighbor as decided in init_fluct_parameter
end function

subroutine init_fluct_parameter(fluct_val,lat,l_use_in,direction)
    use, intrinsic :: iso_fortran_env, only : output_unit
    type(fluct_parameters)     :: fluct_val
    type(lattice),intent(in)    :: lat
    logical,intent(in)          :: l_use_in
    real(8),intent(in)          :: direction(3)

    integer                     :: bnd(2)
    real(8),allocatable         :: diff_vec(:,:)
    real(8),allocatable         :: proj(:) 

    fluct_val%l_use=l_use_in
    if(fluct_val%l_use)then
        Call fluct_val%neigh%get([1,1],[0,1],lat)   !so far only fluctuations implemented between next nearest neighbors of first atom type
        bnd=[fluct_val%neigh%Nshell(1)+1,sum(fluct_val%neigh%Nshell(1:2))]
        allocate(diff_vec(3,bnd(2)-bnd(1)+1))
        diff_vec=fluct_val%neigh%diff_vec(:,bnd(1):bnd(2))
        allocate(proj(size(diff_vec,2)),source=0.0d0)
        proj=matmul(direction,diff_vec)
        fluct_val%ind(1)=maxloc(proj,1)
        proj=matmul(diff_vec(:,fluct_val%ind(1)),diff_vec)
        fluct_val%ind(2)=minloc(proj,1)
        fluct_val%ind=fluct_val%ind+bnd(1)-1

        write(output_unit,'(//A)') "Fluction neighbors are in the following order:"
        write(output_unit,'(3F16.8)') diff_vec
        write(output_unit,'(/A)') "Considering the following neighbors:"
        write(output_unit,'(3F16.8)') fluct_val%neigh%diff_vec(:,fluct_val%ind)

    endif
end subroutine

subroutine eval_fluct(  MipMip_sum, MipMim_sum, MipMjp_sum, MjpMim_sum, lat,fluct_val)
    use m_vector, only : cross
    complex(8),intent(inout)           :: MipMip_sum, MipMim_sum, MipMjp_sum, MjpMim_sum
    type(lattice),intent(in)           :: lat
    type(fluct_parameters),intent(in)  :: fluct_val

    integer         :: i_sh,i_nei
    integer         :: i,j
    real(8)         :: Mi(3), Mj(3)
	complex(8)		:: MipMip, MipMim, MipMjp, MjpMim
    integer         :: neigh_start
    integer         :: ind

    if(.not.fluct_val%l_use) return

    MipMip=(0.0d0,0.0d0)
    MipMjp=(0.0d0,0.0d0)
    MipMim=(0.0d0,0.0d0)
    MjpMim=(0.0d0,0.0d0)

    !do onsite-terms
    neigh_start=1
    do i_sh=1,fluct_val%neigh%Nshell(1)
        do i_nei=neigh_start,fluct_val%neigh%ishell(i_sh) 
            i=fluct_val%neigh%pairs(1,i_nei)
            Mi=lat%M%modes_3(:,i)
            MipMip=MipMip+cmplx(Mi(1)**2-Mi(2)**2, 2.0d0*Mi(1)*Mi(2),8)
            MipMim=MipMim+cmplx(Mi(1)**2+Mi(2)**2,     0.0d0        ,8)
        enddo
        neigh_start=fluct_val%neigh%ishell(i_sh)+1
    enddo

    !do nearest neighbor terms
    do i_sh=1,size(fluct_val%ind)
        ind=fluct_val%ind(i_sh)
        do i_nei=fluct_val%neigh%ishell(ind-1)+1, fluct_val%neigh%ishell(ind) 
            i=fluct_val%neigh%pairs(1,i_nei)
            j=fluct_val%neigh%pairs(2,i_nei)
            Mi=lat%M%modes_3(:,i)
            Mj=lat%M%modes_3(:,j)
            MipMjp=MipMjp + cmplx(Mi(1)*Mj(1)-Mi(2)*Mj(2)  ,  Mi(1)*Mj(2)+Mi(2)*Mj(1),8)     !why is here a minus?
            MjpMim=MjpMim + cmplx(Mi(1)*Mj(1)+Mi(2)*Mj(2)  ,  Mi(1)*Mj(2)-Mi(2)*Mj(1),8)
        enddo
    enddo

    MipMip_sum=MipMip_sum+MipMip
    MipMim_sum=MipMim_sum+MipMim
    MipMjp_sum=MipMjp_sum+MipMjp
    MjpMim_sum=MjpMim_sum+MjpMim
end subroutine


subroutine eval_fluct_spatial( MjpMim_ij_sum,lat,fluct_val)
    use m_vector, only : cross
    complex(8),intent(inout)           :: MjpMim_ij_sum(:,:)
    type(lattice),intent(in)           :: lat
    type(fluct_parameters),intent(in)  :: fluct_val

    integer         :: i_sh,i_nei
    integer         :: i,j
    real(8)         :: Mi(3), Mj(3) !temporary arrays for the chosen magnetizations
    integer         :: ind

    if(.not.fluct_val%l_use) return

    do i_sh=1,size(fluct_val%ind)
        ind=fluct_val%ind(i_sh)
        do i_nei=fluct_val%neigh%ishell(ind-1)+1, fluct_val%neigh%ishell(ind) 
            i=fluct_val%neigh%pairs(1,i_nei)
            j=fluct_val%neigh%pairs(2,i_nei)
            Mi=lat%M%modes_3(:,i)
            Mj=lat%M%modes_3(:,j)
            MjpMim_ij_sum(i_sh,i)=MjpMim_ij_sum(i_sh,i)+cmplx(Mi(1)*Mj(1)+Mi(2)*Mj(2)  ,  Mi(1)*Mj(2)-Mi(2)*Mj(1),8)
        enddo
    enddo
end subroutine

!!!! print spatial distribution for <Mj+Mi->
subroutine print_fluct_spatial(N_add,temp,fluct,com)
    use m_constants, only : k_b
    use m_convert
    use mpi_util, only: reduce_sum, mpi_type
    integer,intent(in)              :: N_add        !on com master thread number of times the fluct array has been added (already reduced of the mpi-communicator)
    real(8),intent(in)              :: temp         !on com master thread temperature
    complex(8),intent(inout)        :: fluct(:,:)   !summed up spatial fluctuation contributions (initially each thread only contains its own contributions, on output destroyed)
    class(mpi_type),intent(in)      :: com          !inner mpi-communicator that contains all threads with the same temperature
    integer             :: io
    character(len=50)   :: frm

    Call reduce_sum(fluct,com)
    if(com%ismas)then
        fluct=fluct/real(N_add,8)  !normalization (number of times the fluct array has been added with a configuration)
        
        frm=convert('(',size(fluct,1),'(E20.12E3,2x))')

        open(newunit=io,file=convert('fluct_Re_MjpMim_per_site_',temp,'.dat'))
        write(io,frm) real(fluct)
        close(io)

        open(newunit=io,file=convert('fluct_Im_MjpMim_per_site_',temp,'.dat'))
        write(io,frm) aimag(fluct)
        close(io)
    endif
    fluct=(0.0d0,0.0d0) !set to zero since the values now are nonsensical anyways
end subroutine


end module 

