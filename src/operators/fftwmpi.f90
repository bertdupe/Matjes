module m_fftwmpi
use m_fftw3
!    ! kmesh in internal units
    integer, protected, public :: N_kpoint(3)
    real(kind=8), protected, public, allocatable :: kmesh(:,:)

    private
    contains
#if 0
#ifdef CPP_FFTWMPI
        !!!!!!!!!!!!!!!!!!!!!!
        ! FFT dipole
        !!!!!!!!!!!!!!!!!!!!!!
        subroutine fft_mpi_depol_r2c(Nx,Ny,Nz,rmat,cmat,rtrans,ctrans,plan)
            use, intrinsic :: iso_c_binding
            implicit none
            integer, intent(in) :: Nx,Ny,Nz
            real(kind=8), intent(in) :: rmat(:,:,:,:)
            complex, intent(out) :: cmat(:,:,:,:)
            type(C_PTR), intent(in) :: plan
            real(c_double), intent(inout), pointer :: rtrans(:,:,:)
            complex(c_double_complex), intent(inout), pointer :: ctrans(:,:,:)

            ! Internal variables
            integer :: i,j,k,l
            integer :: N
            integer(C_INTPTR_T) :: i, j, alloc_local, local_Nz, local_Nz_offset

            N=size(rmat,1)

            !   get local data size and allocate (note dimension reversal)
            alloc_local = fftw_mpi_local_size_3d(Nz, Ny, Nx/2+1, MPI_COMM_WORLD,local_Nz, local_Nz_offset)

            ctrans=fftw_alloc_complex(alloc_local)
            call c_f_pointer(cdata, data, [Nx/2,Ny,local_Nz])

            ! split the real space in data
            rtrans = fftw_alloc_real(2*alloc_local)
            call c_f_pointer(rdata, in, [2*(Nx/2+1),Ny,local_Nz])

            ! create the plan for the FFTW_MPI
            plan = fftw_mpi_plan_dft_r2c_3d(Nz, Ny, Nx, rtrans, ctrans, MPI_COMM_WORLD,FFTW_MEASURE)

            do i=1,N
                rtrans(1:Nx,1:Ny,1:Nz)=rmat(i,1:Nx,1:Ny,1:Nz)
                call fftw_mpi_execute_dft_r2c_3d(plan,rtrans,ctrans)
                cmat(i,1:Nx,1:Ny,1:Nz)=ctrans(1:Nx,1:Ny,1:Nz)
            enddo

            call fftw_destroy_plan(plan)
            call fftw_free(cdata)

        end subroutine fft_mpi_depol_r2c



        !!!!!!!!!!!!!!!!!!!!!!
        ! FFT dipole
        !!!!!!!!!!!!!!!!!!!!!!
        subroutine fft_mpi_depol_c2r(Nx,Ny,Nz,cmat,rmat,alpha,rtrans,ctrans,plan)
            use, intrinsic :: iso_c_binding
            implicit none
            integer, intent(in) :: Nx,Ny,Nz
            integer,intent(in) :: alpha
            real(kind=8), intent(out) :: rmat(:,:,:,:)
            complex, intent(in) :: cmat(:,:,:,:)
            type(C_PTR), intent(in) :: plan
            real(c_double), intent(inout) :: rtrans(:,:,:)
            complex(c_double_complex), intent(inout) :: ctrans(:,:,:)

            ! Internal variables
            integer :: i,j,k,l
            integer :: N

            ctrans=dcmplx(0.d0,0.d0)
            rtrans=0.0d0

            N=size(cmat,1)

            do i=1,N
                ctrans(1:Nx,1:Ny,1:Nz)=cmat(i,1:Nx,1:Ny,1:Nz)
                call fftw_mpi_execute_dft_c2r(plan,ctrans,rtrans)
                rmat(i,1:Nx,1:Ny,1:Nz)=rtrans(1:Nx,1:Ny,1:Nz)/dble(Nx*Ny*Nz)
            enddo
            !      call fftw_destroy_plan(plan)
        end subroutine fft_mpi_depol_c2r
#endif



        ! Calculate the Fourrier transform of a field of pointers
        subroutine calculate_fft_mpi_vec_point(field,pos,sense,dim_mode,FFT)
            use m_derived_types, only : vec_point
            implicit none
            complex(kind=8), intent(inout) :: FFT(:,:)
            type(vec_point), intent(in) :: field(:)
            real(kind=8), intent(in) :: sense       ! should be + or - one
            real(kind=8), intent(in) :: pos(:,:)
            integer, intent(in) :: dim_mode

            ! Internal variables
            integer :: i,Nsize,N_k

            Nsize=size(field)
            N_k=size(kmesh,2)
            FFT=0.0d0

            do i=1,N_k
                FFT(:,i)=get_FFT(kmesh(:,i),sense,pos,field,Nsize,dim_mode)
            enddo
        end subroutine calculate_fft_vec_point



        ! Calculate the Fourrier transform of a field of reals
        subroutine calculate_fft_mpi_matrix(field,pos,sense,FFT)
            use m_get_position
            implicit none
            real(kind=8), intent(in) :: field(:,:),pos(:,:)
            real(kind=8), intent(in) :: sense       ! should be + or - one
            complex(kind=8), intent(inout) :: FFT(:,:)

            ! Internal variables
            integer :: i,N_k,N(2)

            N=shape(field)
            N_k=size(kmesh,2)

            do i=1,N_k
                FFT(:,i)=get_FFT(field,sense,kmesh(:,i),pos,N(1))
            enddo
        end subroutine calculate_fft_mpi_matrix



        ! Calculate the Fourrier transform of the D(r) matrix
        subroutine calculate_FFT_mpi_Dr(pos,sense,dim_lat)
            use m_vector, only : norm
            implicit none
            real(kind=8), intent(in) :: sense       ! should be + or - one
            real(kind=8), intent(in) :: pos(:,:)
            integer, intent(in) :: dim_lat(:)

            ! Internal variables
            integer :: Nksize,Nsize,i
            real(kind=8) :: r(3),norm_int
            real(kind=8), allocatable :: pos_D(:,:,:)

            Nksize=product(N_kpoint)
            Nsize=product(dim_lat)
            allocate(FFT_pos_D(3,3,Nksize),pos_D(3,3,Nsize))
            FFT_pos_D=0.0d0

            do i=1,Nsize
                r=pos(:,i)
                norm_int=norm(r)
                if (norm_int.gt.1.0d-8) then
                    pos_D(:,1,i)=(/ r(1)**2 - sum(r**2)/3 , r(1)*r(2) , r(1)*r(3) /)/norm_int**5
                    pos_D(:,2,i)=(/ r(1)*r(2) , r(2)**2 - sum(r**2)/3 , r(2)*r(3) /)/norm_int**5
                    pos_D(:,3,i)=(/ r(1)*r(3) , r(2)*r(3) , r(3)**2 - sum(r**2)/3 /)/norm_int**5
                endif
            enddo

            do i=1,Nksize
                FFT_pos_D(:,:,i)=get_FFT(pos_D,sense,kmesh(:,i),pos,Nsize)
            enddo
        end subroutine calculate_FFT_mpi_Dr



        ! Function that calculate the FFT of a matrix of vec_point for one component
        function get_FFT_mpi_vec_point(kvec,sense,pos,field,Nsize,dim_mode)
            use m_derived_types, only : vec_point
            implicit none
            integer, intent(in) :: Nsize,dim_mode
            real(kind=8), intent(in) :: kvec(:),pos(:,:),sense
            type(vec_point), intent(in) :: field(:)
            complex(kind=8) :: get_FFT_vec_point(dim_mode)

            ! Internal variables
            real(kind=8) :: phase,r(3)
            integer :: k,j

            get_FFT_vec_point=0.0d0

            do j=1,Nsize
                r=pos(:,j)
                phase=dot_product(kvec,r)
                do k=1,dim_mode
                    get_FFT_vec_point(k)=get_FFT_vec_point(k)+cmplx(field(j)%w(k)*cos(sense*phase),field(j)%w(k)*sin(sense*phase),8)
                enddo
            enddo
            get_FFT_vec_point=get_FFT_vec_point/sqrt(real(Nsize,8))
        end function get_FFT_mpi_vec_point



        ! Function that calculate the FFT of the D(r) matrix of real for one component
        function get_FFT_mpi_Dr_matrix(pos_D,sense,kvec,real_dist,Nsize)
            implicit none
            integer, intent(in) :: Nsize
            real(kind=8), intent(in) :: pos_D(:,:,:),sense,kvec(:),real_dist(:,:)
            complex(kind=8) :: get_FFT_Dr_matrix(3,3)

            ! Internal variables
            integer :: i,l,m
            real(kind=8) :: phase,r(3)

            get_FFT_Dr_matrix=0.0d0

            do i=1,Nsize
                r=real_dist(:,i)
                phase=dot_product(kvec,r)
                do l=1,3
                    do m=1,3
                        get_FFT_Dr_matrix(m,l)=get_FFT_Dr_matrix(m,l)+cmplx(pos_D(m,l,i)*cos(sense*phase),pos_D(m,l,i)*sin(sense*phase),8)
                    enddo
                enddo
            enddo
            get_FFT_Dr_matrix=get_FFT_Dr_matrix/real(Nsize,8)
        end function get_FFT_mpi_Dr_matrix



        ! Function that calculate the FFT of a matrix of real for one component
        function get_FFT_mpi_matrix(field,sense,kvec,real_dist,dim_mode)
            implicit none
            real(kind=8), intent(in) :: field(:,:),sense,kvec(:),real_dist(:,:)
            integer, intent(in) :: dim_mode
            complex(kind=8) :: get_FFT_matrix(dim_mode)

            ! Internal variables
            integer :: shape_pos(2),i,l
            real(kind=8) :: alpha,r(3)

            shape_pos=shape(field)

            get_FFT_matrix=0.0d0

            do i=1,shape_pos(2)
                r=real_dist(:,i)
                alpha=dot_product(kvec,r)
                do l=1,shape_pos(1)
                    get_FFT_matrix(l)=get_FFT_matrix(l)+cmplx(field(l,i)*cos(sense*alpha),field(l,i)*sin(sense*alpha),8)
                enddo
            enddo
            get_FFT_matrix=get_FFT_matrix/real(shape_pos(2),8)
        end function get_FFT_mpi_matrix


        !Function computing the product X*e^{ikr}*H(r-r')*X*e^{ikr'} for a given r and kvec
        !all_mode contains the values in X
        !pos contains the positions of the sites and neighbours of the sites
        !kvec is the k-vector multiplying the position in the exponential
        !sense is the sense of the transform (-1.0 or +1.0)
        function get_E_k_local(all_energy, pos, all_vectors, X_site, pos_k, size_all_vectors, dim_mode, sense)
            implicit none
            integer :: pos_k, dim_mode
            integer :: size_all_vectors
            complex(kind=8) :: all_energy(:,:), all_vectors(:), X_site(:)
            real(kind=8) :: pos(:,:), sense !pos is the array of all r-r'
            complex(kind=8) :: get_E_k_local

            ! Internal variable
            integer :: i
            complex(kind=8) :: intermediate_sum(size_all_vectors)
            real(kind=8) :: alpha

            get_E_k_local=0.0d0

            do i=1, size_all_vectors, dim_mode
                alpha = dot_product( pos(:,(i/dim_mode)+1), kmesh(:,pos_k) )
                intermediate_sum(i:i+dim_mode-1) = all_vectors(i:i+dim_mode-1)*cmplx( cos(sense*alpha), sin(sense*alpha),8 )!gives exp[i*k*(r-r')] (all r' and one particular r)
            enddo

            get_E_k_local=dot_product( X_site, matmul(all_energy, intermediate_sum) )
        end function get_E_k_local



        ! Function computing the Fourier transform of an input Hamiltonian
        ! Input:
        !   _ all_E_k: input Hamiltonian that will be Fourier transformed
        !   _ n_lines: non-zero elements in all_E_k
        !   _ pos: vector containing the position of all non-zero elements in the Hamiltonian
        !   _ pos_k: position in the kmesh of the k-vector needed in the Fourier transform
        !   _ dim_mode: length of the order parameter
        !   _ sense: gives the sense of the transform (-1.0d0 ==> direct, +1.0d0 ==> indirect)
        ! Output:
        !   _ Fourier transform of the input all_E_k
        function Fourier_transform_H(all_E_k, pos, pos_k, dim_mode, sense)result(FT)
            use m_energy_commons, only : energy
            implicit none
            integer,intent(in)         :: pos_k, dim_mode
            real(kind=8),intent(in)    :: pos(:, :), sense !pos is the array of all r-r'
            complex(kind=8),intent(in) :: all_E_k(:, :)
            complex(kind=8)            :: FT( size(all_E_k, 1), size(all_E_k, 2) )

            FT=Fourier_H_at_k(all_E_k, pos,kmesh(:, pos_k), dim_mode, sense)

        end function Fourier_transform_H




        ! Function computing the Fourier transform of an input Hamiltonian
        ! Input:
        !   _ all_E_k: input Hamiltonian that will be Fourier transformed
        !   _ n_lines: non-zero elements in all_E_k
        !   _ pos: vector containing the position of all non-zero elements in the Hamiltonian
        !   _ kpt: k-point
        !   _ dim_mode: length of the order parameter
        !   _ sense: gives the sense of the transform (-1.0d0 ==> direct, +1.0d0 ==> indirect)
        ! Output:
        !   _ Fourier transform of the input all_E_k
        function Fourier_H_at_k(Hr, pos, kpt, dim_mode, sense)result(FT)
            !most basic implementation which does not work for interactions larger than the unit-cell size
            !and which has to loop over the entire Hr matrix
            use m_energy_commons, only : energy
            implicit none
            integer          :: dim_mode
            real(kind=8)     :: kpt(3)
            real(kind=8)     :: pos(:, :), sense !pos is positions(3,Ncell)
            complex(kind=8)  :: Hr(:, :)
            complex(kind=8)  :: FT( size(Hr, 1), size(Hr, 2) )

            ! Internal variable
            integer :: i, j, nblines_energy, nbcols_energy, tmp
            real(kind=8) :: alpha
            integer :: aa,ab,ba,bb
            nblines_energy = size(energy%line, 1)
            nbcols_energy = size(energy%line, 2)

            FT = cmplx(0.0d0,0.0d0,kind=8)
            do i=1, nbcols_energy
               do j=1, nblines_energy
                    tmp = energy%line(j,i)
                    alpha = sense*dot_product( pos(:, j), kpt )
                    aa=(i-1)*dim_mode+1
                    ab=i*dim_mode
                    ba=(tmp-1)*dim_mode+1
                    bb=tmp*dim_mode
                    FT( aa:ab, ba:bb)=FT(aa:ab,ba:bb)+Hr(aa:ab,ba:bb)*cmplx(cos(alpha), sin(alpha),kind=8)
                enddo
            enddo

        end function



#endif

end module m_fftwmpi
