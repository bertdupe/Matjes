module m_highsym
    use m_basic_types, only : vec_point
    use m_io_utils, only: get_parameter
    use m_io_files_utils, only: close_file,open_file_write
    use m_type_lattice, only : lattice
    use m_energy_k,only : get_energy_kpts
    use m_tb_types
    implicit none
    private 
    public :: plot_highsym_kpts,set_highs_path

    integer                 ::  N_kpts
    real(8),allocatable     ::  kpts(:,:)
    real(8),allocatable     ::  kpts_dist(:)
    real(8),allocatable     ::  highs_dist(:)

    contains


    subroutine plot_highsym_kpts(h_par,pos,mode_mag,E_F)
        type(parameters_TB_Hsolve),intent(in)     ::  h_par
        real(8),intent(in)          :: pos(:,:)
        type(vec_point),intent(in)  :: mode_mag(h_par%dimH)
        real(8),intent(in)          :: E_F !fermi energy

        real(8),allocatable :: eigval(:,:)

        if(.not. allocated(kpts)) return
        write(*,'(A,I6,A)') "Calculate ",N_kpts,' kpoints on the high symmetry path'
        Call get_energy_kpts(kpts,h_par,pos,mode_mag,eigval)
        !eigval=eigval-E_F !adjust by calculated fermi energy only confuses comparison with analytic calculations
        Call print_highsym('highs_plot',eigval)
        Call highsym_clear_path()

    end subroutine


    subroutine print_highsym(fname,eigval)
        !TODO: set in some labels...
        character(len=*), intent(in) :: fname
        real(8),intent(in)      ::  eigval(:,:)

        integer                 ::  i,j,io
        io=open_file_write(fname//".dat")
        do i=1,N_kpts
           write(io,'(2E16.8)') (kpts_dist(i),eigval(j,i),j=1,size(eigval,1))
        enddo
        call close_file(fname//".dat",io)
        io=open_file_write(fname//".gnu")
        write(io,'(A)') 'set xtics( \'
        do i=1,size(highs_dist,1)-1
            write(io,"(A,E16.8,A)") '"label"',highs_dist(i),', \'
        enddo
        write(io,"(A,E16.8,A)") '"label"',highs_dist(size(highs_dist,1)), '\'
        write(io,'(A)') ')'
        write(io,'(A)') 'set grid x linetype -1'
        write(io,'(A)') 'set ylabel "E-E_F"'
        write(io,'(A)') 'set nokey'
        write(io,'(A)') "plot 'highs_plot.dat'"
        write(io,'(A)') 'pause -1 "Hit any key to continue"'
        !should somehow also add labels
        call close_file(fname//".gnu",io)

    end subroutine

    subroutine highsym_clear_path()
        deallocate(kpts)
        deallocate(kpts_dist)
        deallocate(highs_dist)
    end subroutine 

    subroutine set_highs_path(my_lattice,highs)
        type(lattice), intent(in)                   :: my_lattice
        type(parameters_TB_IO_highs),intent(in)     :: highs

        integer                 ::  N_highsym
        real(8),allocatable     ::  k_highs(:,:)
        integer                 ::  i,j
        integer,allocatable     ::  Npath_k(:)
        real(8)                 ::  distance
        integer                 ::  i_k
        real(8)                 ::  diff_vec(3)
        integer,allocatable     ::  highs_ind(:)

        if(highs%N_highsym < 1)then
            return
        endif

        N_highsym=highs%N_highsym
        allocate(k_highs,source=highs%k_highs)
        !set high symmetry kpoints in correct units
        do i=1,N_highsym
            k_highs(:,i)=k_highs(:,i)/highs%k_highs_frac(i)
        enddo
        do i=1,N_highsym
            k_highs(:,i)=matmul(k_highs(:,i),my_lattice%astar)
        enddo
        !adjust with dimension of real-space cell decreasing BZ
        !might only works for rectangular lattice?
        do i=1,N_highsym
            k_highs(:,i)=k_highs(:,i)/real(my_lattice%dim_lat,8)
        enddo

        !set kpts array with all kpoints along the path
        allocate(Npath_k(N_highsym-1),source=1)
        do i=1,N_highsym-1 
            distance=norm2(k_highs(:,i+1)-k_highs(:,i))
            Npath_k(i)=max(nint(distance/highs%aim_dist)-1,0)
        enddo
        N_kpts=N_highsym+sum(Npath_k)
        allocate(kpts(3,N_kpts),source=0.0d0)
        allocate(highs_ind(N_highsym),source=1)
        kpts(:,1)=k_highs(:,1)
        i_k=2
        do i=1,N_highsym-1
            diff_vec=(k_highs(:,i+1)-k_highs(:,i))/real(Npath_K(i)+1,kind=8)
            do j=1,Npath_K(i)
                kpts(:,i_k)=k_highs(:,i)+real(j,kind=8)*diff_vec
                i_k=i_k+1
            enddo
            kpts(:,i_k)=k_highs(:,i+1)
            highs_ind(i+1)=i_k
            i_k=i_k+1
        enddo
        allocate(kpts_dist(N_kpts),source=0.0d0)
        do i=2,N_kpts
            kpts_dist(i)=kpts_dist(i-1)+norm2(kpts(:,i)-kpts(:,i-1))
        enddo
        allocate(highs_dist(N_highsym))
        highs_dist=kpts_dist(highs_ind)

    end subroutine 


end module m_highsym
