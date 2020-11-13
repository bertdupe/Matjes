module m_fftw

!    ! kmesh in internal units
    integer, protected, public :: N_kpoint(3)
    real(kind=8), protected, public, allocatable :: kmesh(:,:)

!    ! FFT of the dipolar martix D(r)
!    complex(kind=8), allocatable, protected, public :: FFT_pos_D(:,:,:)
!
!    interface calculate_fft
!      module procedure calculate_fft_vec_point,calculate_FFT_matrix
!    end interface
!
!    interface get_FFT
!      module procedure get_FFT_vec_point,get_FFT_matrix,get_FFT_Dr_matrix,get_E_k_local, Fourier_transform_H,Fourier_H_at_k
!    end interface

    private
    public :: set_k_mesh
!    public :: calculate_fft,calculate_FFT_Dr,get_FFT, Fourier_H_at_k
    contains
#if 0
#ifdef CPP_FFTW
        !!!!!!!!!!!!!!!!!!!!!!
        ! FFT dipole
        !!!!!!!!!!!!!!!!!!!!!!
        subroutine fft_depol_r2c(Nx,Ny,Nz,rmat,cmat,rtrans,ctrans,plan)
            use, intrinsic :: iso_c_binding
            implicit none
            integer, intent(in) :: Nx,Ny,Nz
            real(kind=8), intent(in) :: rmat(:,:,:,:)
            complex, intent(out) :: cmat(:,:,:,:)
            type(C_PTR), intent(in) :: plan
            real(c_double), intent(inout) :: rtrans(:,:,:)
            complex(c_double_complex), intent(inout) :: ctrans(:,:,:)

            ! Internal variables
            integer :: i,j,k,l
            integer :: N

            include 'fftw3.f03'

            N=size(rmat,1)

            ctrans=dcmplx(0.d0,0.d0)
            rtrans=0.0d0

            do i=1,N
                rtrans(1:Nx,1:Ny,1:Nz)=rmat(i,1:Nx,1:Ny,1:Nz)
                call fftw_execute_dft_r2c(plan,rtrans,ctrans)
                cmat(i,1:Nx,1:Ny,1:Nz)=ctrans(1:Nx,1:Ny,1:Nz)
            enddo
        end subroutine fft_depol_r2c



        !!!!!!!!!!!!!!!!!!!!!!
        ! FFT dipole
        !!!!!!!!!!!!!!!!!!!!!!
        subroutine fft_depol_c2r(Nx,Ny,Nz,cmat,rmat,alpha,rtrans,ctrans,plan)
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

            include 'fftw3.f03'

            ctrans=dcmplx(0.d0,0.d0)
            rtrans=0.0d0

            N=size(cmat,1)

            do i=1,N
                ctrans(1:Nx,1:Ny,1:Nz)=cmat(i,1:Nx,1:Ny,1:Nz)
                call fftw_execute_dft_c2r(plan,ctrans,rtrans)
                rmat(i,1:Nx,1:Ny,1:Nz)=rtrans(1:Nx,1:Ny,1:Nz)/dble(Nx*Ny*Nz)
            enddo
            !      call fftw_destroy_plan(plan)
        end subroutine fft_depol_c2r
#endif



        ! Calculate the Fourrier transform of a field of pointers
        subroutine calculate_fft_vec_point(field,pos,sense,dim_mode,FFT)
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
        subroutine calculate_fft_matrix(field,pos,sense,FFT)
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
        end subroutine calculate_fft_matrix



        ! Calculate the Fourrier transform of the D(r) matrix
        subroutine calculate_FFT_Dr(pos,sense,dim_lat)
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
        end subroutine calculate_FFT_Dr



        ! Function that calculate the FFT of a matrix of vec_point for one component
        function get_FFT_vec_point(kvec,sense,pos,field,Nsize,dim_mode)
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
            get_FFT_vec_point=get_FFT_vec_point/sqrt(real(Nsize))
        end function get_FFT_vec_point



        ! Function that calculate the FFT of the D(r) matrix of real for one component
        function get_FFT_Dr_matrix(pos_D,sense,kvec,real_dist,Nsize)
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
            get_FFT_Dr_matrix=get_FFT_Dr_matrix/real(Nsize)
        end function get_FFT_Dr_matrix


    
        ! Function that calculate the FFT of a matrix of real for one component
        function get_FFT_matrix(field,sense,kvec,real_dist,dim_mode)
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
            get_FFT_matrix=get_FFT_matrix/real(shape_pos(2))
        end function get_FFT_matrix


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

        ! Get the mesh for the Fourrier transform
        subroutine set_k_mesh(fname,my_lattice)
            use m_kmesh
            use m_derived_types, only : lattice
            use m_io_utils
            use m_io_files_utils
            implicit none
            character(len=*), intent(in) :: fname
            type(lattice), intent(in) :: my_lattice

            ! Internal variables
            integer :: io_input,test
            integer :: Nkpoint,i,iomp
            logical :: i_plot,i_file

            N_kpoint=my_lattice%dim_lat
            i_plot=.false.
            io_input=open_file_read(fname)
            call get_parameter(io_input,fname,'kmesh',3,N_kpoint)
            call get_parameter(io_input,fname,'write_kmesh',i_plot)
            call close_file(fname,io_input)

            if(allocated(kmesh))then
                write(*,*) "WARNING: kmesh is already allocated"
                return
            endif

            i_file=.false.
            inquire(file='kpoints',exist=i_file)
            if (i_file) then
                write(*,'(A)') "Reading kpoint mesh from file: kpoints"
                io_input=open_file_read('kpoints')
                Nkpoint=0
                do 
                    read(io_input,*,iostat=test) 
                    if(test/=0) exit
                    Nkpoint=Nkpoint+1
                enddo
                rewind(io_input)
                allocate(kmesh(3,Nkpoint),stat=test)
                if (test.ne.0) return
                do iomp=1,Nkpoint
                    read(io_input,*) (kmesh(i,iomp),i=1,3)
                enddo
                call close_file('kpoints',io_input)
                !PB it might make sense to write the gridsize into the kpoints file...
                N_kpoint=[Nkpoint,1,1]
            else
                Nkpoint=product(N_kpoint)
                ! Check if the variable kmesh is allocated
                allocate(kmesh(3,Nkpoint),stat=test)
                if (test.ne.0) return
                kmesh=0.0d0
                call get_kmesh(N_kpoint,kmesh,i_plot)
                !PB: rescale Bz with lattice size !might fail with hexagonal unit cell??? <- CHECK THAT
                do iomp=1,Nkpoint
                    kmesh(:,iomp)=kmesh(:,iomp)/my_lattice%dim_lat
                enddo
            endif

            ! Convert the kmesh in internal units
            do iomp=1,Nkpoint
                kmesh(:,iomp)=matmul(kmesh(:,iomp),my_lattice%astar)
            enddo
        end subroutine set_k_mesh
end module m_fftw
