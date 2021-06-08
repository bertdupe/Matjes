module m_set_Hamiltonians
use m_H_public
implicit none
private
public :: set_Hamiltonians

contains

subroutine set_Hamiltonians(Ham_res,Ham_comb,keep_res,H_io,lat)
    use m_derived_types
    use m_input_H_types 
    use m_get_position,only :get_position_ND_to_1D 
    
    use m_symmetry_operators
    use m_anisotropy_heisenberg,only: get_anisotropy_H
    !use m_exchange_heisenberg,only: get_exchange_H,exchange
    use m_zeeman,only: get_zeeman_H
    use m_exchange_heisenberg_J, only: get_exchange_J
    use m_exchange_heisenberg_D, only: get_exchange_D
    use m_coupling_ME_J,only: get_coupling_ME_J
    use m_coupling_ME_D,only: get_coupling_ME_D
    use m_harmonic_phonon,only: get_Forces_F
    use m_exchange_TJ,only: get_exchange_TJ
    use m_ASR_phonon, only: get_ASR_Ph
    use m_Mag_Biq, only: get_Mag_Biq
    use m_4spin, only: get_4spin
    use m_dipolar_magnetic, only: get_dipolar
    class(t_H),allocatable,intent(out)  :: Ham_res(:)
    class(t_H),allocatable,intent(out)  :: Ham_comb(:)
    logical,intent(in)                  :: keep_res ! keeps the Ham_res terms allocated
    type(io_h),intent(in)               :: H_io
    type(lattice), intent(in)           :: lat

    integer :: i_H,N_ham
    logical :: use_Ham(12)


    use_ham(1)=H_io%J%is_set.and..not.H_io%J%fft
    use_ham(2)=H_io%D%is_set.and..not.H_io%D%fft
    use_ham(3)=H_io%aniso%is_set.and..not.H_io%aniso%fft
    use_ham(4)=H_io%zeeman%is_set
    use_ham(5)=H_io%ME_J%is_set
    use_ham(6)=H_io%ME_D%is_set
    use_ham(7)=H_io%F%is_set
    use_ham(8)=H_io%TJ%is_set
    use_ham(9)=H_io%ASR_Ph%is_set
    use_ham(10)=H_io%M_biq%is_set
    use_ham(11)=H_io%sp4%is_set
    use_ham(12)=H_io%dip%is_set.and..not.H_io%dip%fft

    N_ham=count(use_ham)
    Call get_Htype_N(Ham_res,N_ham)
    i_H=1 
    !exchange_J (without DMI)
    if(use_ham(1))then
        Call get_exchange_J(Ham_res(i_H),H_io%J,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !exchange_D (only DMI)
    if(use_ham(2))then
        Call get_exchange_D(Ham_res(i_H),H_io%D,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !anisotropy
    if(use_ham(3))then
        Call get_anisotropy_H(Ham_res(i_H),H_io%aniso,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !zeeman
    if(use_ham(4))then
        Call get_zeeman_H(Ham_res(i_H),H_io%zeeman,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !ME-coupling symmetric (J)
    if(use_ham(5))then
        Call get_coupling_ME_J(Ham_res(i_H),H_io%ME_J,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !ME-coupling antisymmetric (D)
    if(use_ham(6))then
        Call get_coupling_ME_D(Ham_res(i_H),H_io%ME_D,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !Harmonic phonon (F)
    if(use_ham(7))then
        Call get_Forces_F(Ham_res(i_H),H_io%F,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !TJ coupling
    if(use_ham(8))then
        Call get_exchange_TJ(Ham_res(i_H),H_io%TJ,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !ASR for the phonon
    if(use_ham(9))then
        Call get_ASR_Ph(Ham_res(i_H),H_io%ASR_Ph,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !Magnetic biquadratic
    if(use_ham(10))then
        Call get_Mag_Biq(Ham_res(i_H),H_io%M_biq,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !4-spin interaction
    if(use_ham(11))then
        Call get_4spin(Ham_res(i_H),H_io%sp4,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif
    !direct calculation of dipolar interaction
    if(use_ham(12))then
        Call get_dipolar(Ham_res(i_H),H_io%dip,lat)
        if(Ham_res(i_H)%is_set()) i_H=i_H+1
    endif

    do i_H=1,N_ham
        if(.not. Ham_res(i_h)%is_set()) STOP "not all Hamiltonians are set"
        Call Ham_res(i_h)%optimize()
        Call Ham_res(i_h)%finish_setup()
    enddo

    Call combine_Hamiltonians(keep_res,Ham_res,Ham_comb)
    do i_H=1,size(Ham_comb)
        Call Ham_comb(i_h)%optimize()
        Call Ham_comb(i_h)%finish_setup()
    enddo
end subroutine

subroutine combine_Hamiltonians(keep_in,H_in,H_comb)
    !subroutine which combines Hamiltonians from the input array, if they act on the same space
    !THIS NEEDS TO BE UPDATED IF THE SPACE IS NOT UNIQUE ANYMORE (eg. terms in one dimension that are not onsite)
    logical,intent(in)                      :: keep_in
    class(t_H),intent(inout)                :: H_in(:)
    class(t_H),allocatable,intent(inout)    :: H_comb(:)
    !internal
    integer         :: N_in,N_out
    integer         :: max_dim
    integer,allocatable ::  space(:,:,:)    !save all space indicators
    integer         ::  i_unique(size(H_in))
    integer         ::  i,j

    !find out which operators are in which space
    N_in=size(H_in)
    max_dim=0
    do i=1,N_in
        max_dim=max(max_dim,size(H_in(i)%op_l))
        max_dim=max(max_dim,size(H_in(i)%op_r))
    enddo
    allocate(space(max_dim,2,N_in),source=0)
    do i=1,N_in
        do j=1,size(H_in(i)%op_l)
            space(j,1,i)=H_in(i)%op_l(j)
        enddo
        do j=1,size(H_in(i)%op_r)
            space(j,2,i)=H_in(i)%op_r(j)
        enddo
    enddo
    N_out=0
    outer:do j=1,N_in
        do i=1,j-1
            if(all(space(:,:,j)==space(:,:,i)))then
                if(H_in(j)%same_space(H_in(i)))then
                    i_unique(j)=i_unique(i)
                    cycle outer
                endif
            endif
        enddo
        N_out=N_out+1
        i_unique(j)=N_out
    enddo outer

    Call get_Htype_N(H_comb,N_out)
    do i=1,N_in
        Call H_comb(i_unique(i))%add(H_in(i))
        if(.not.keep_in) Call H_in(i)%destroy()
    enddo

end subroutine
end module
