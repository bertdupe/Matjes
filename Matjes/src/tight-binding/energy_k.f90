module m_energy_k
    use m_tb_types
    use m_basic_types, only : vec_point
    use m_energy_commons, only : energy
    use m_fftw, only: get_FFT
    use m_J_sd_exchange
    use m_energy_r, only: set_Hr

    implicit none

    private
    public :: diagonalise_H_k_list,set_dist_neigh,get_energy_kpts !, fermi_distrib, compute_Etot
    contains

        subroutine get_energy_kpts(klist,h_par,pos,mode_mag,eigval)
            real(8),intent(in)                      ::  klist(:,:)
            type(parameters_TB_Hsolve),intent(in)     ::  h_par
            real(8),intent(in)                      ::  pos(:,:)
            type(vec_point),intent(in)              ::  mode_mag(:)
            real(8),intent(inout),allocatable       ::  eigval(:,:)

            complex(8),allocatable      ::  Hr(:,:)
            integer                     ::  dim_mode
            
            !set real space Hamiltonian
            !if used more often, set this in advance
            Call set_Hr(h_par,mode_mag,Hr)
!            Call get_Hr(h_par%dimH,Hr)
            !get actual eigenvalues
            allocate( eigval(h_par%dimH,size(klist,2)),source=0.0d0)
            dim_mode=h_par%pos_ext(2)-h_par%pos_ext(1)+1
            Call diagonalise_H_k_list(Hr,pos,klist,dim_mode,-1.0d0,eigval)
        end subroutine

        subroutine set_dist_neigh(dist_neigh,pos)
            real(8),intent(out),allocatable   ::  dist_neigh(:,:)
            real(8),intent(in)                  ::  pos(:,:)
            integer :: N_neighbours,i,j

            N_neighbours = size( energy%line, 1 )
            allocate(dist_neigh(3,N_neighbours),source=0.0d0)
            do i=1, N_neighbours
               j = energy%line(i, 1)
               dist_neigh(:,i) = pos(:,j)
            enddo
        end subroutine

        subroutine diagonalise_H_k_list(Hr,pos,kpts,dim_mode,sense,eigval)
            external :: ZHEEV
            complex(8),intent(in)   ::  Hr(:,:)
            real(8),intent(in)      ::  pos(:,:)
            real(8),intent(in)      ::  kpts(:,:)
            integer, intent(in)     ::  dim_mode
            real(kind=8), intent(in) :: sense
            real(kind=8), intent(inout),allocatable :: eigval(:,:)

            integer :: i
            integer :: N,Nkpt, INFO
            complex(kind=8), allocatable :: WORK(:), H_complex(:,:)
            real(kind=8), allocatable :: RWORK(:)

            Nkpt=size(kpts,2)
            N = size(Hr, 1)
            if(.not. allocated(eigval)) allocate(eigval(N,Nkpt),source=0.0d0)
            if(size(eigval,2)/= Nkpt) STOP "eigval array does not have correct dimension(Nkpts)"
            allocate( H_complex(N,N), RWORK(max(1,3*N-2)), WORK(2*N))
            do i=1,Nkpt
                H_complex=get_FFT(Hr, pos, kpts(:,i), dim_mode, sense)
                call ZHEEV( 'N', 'U', N, H_complex, N, eigval(:,i), WORK, size(Work), RWORK, INFO )
            enddo
        end subroutine

end module m_energy_k
