module m_highsym
    use m_derived_types, only : lattice
    use m_Hk
    use m_tb_types
    implicit none
    private 
    public :: plot_highsym_kpts,set_highs_path,plot_highsym_kpts_sc, highsym_clear_path, plot_highsym_kpts_proj

    integer                 ::  N_kpts
    real(8),allocatable     ::  kpts(:,:)
    real(8),allocatable     ::  kpts_dist(:)
    real(8),allocatable     ::  highs_dist(:)

    contains

    subroutine plot_highsym_kpts_sc(Hk_inp,h_io)
        type(Hk_inp_t),intent(in)               :: Hk_inp
        type(parameters_TB_IO_H),intent(in)     :: h_io

        real(8),allocatable     :: eval(:)
        complex(8),allocatable  :: evec(:,:)
        character(len=*),parameter  ::  fname="highs_plot_sc"
        integer             ::  io,i,j
        integer             ::  bnd
        real(8)             ::  weight

        if(.not. allocated(kpts)) return
        write(*,'(A,I6,A)') "Calculate ",N_kpts,' kpoints on the high symmetry path'
        !get and print data
        open(newunit=io,file=fname//'.dat')
        bnd=Hk_inp%H(1)%dimH/2
        do i=1,N_kpts
            Call Hk_evec(Hk_inp,kpts(:,i),h_io,eval,evec)
            do j=1,size(eval)
                weight=real(dot_product(evec(1:bnd,j),evec(1:bnd,j)),8)
                write(io,'(3E16.8)') kpts_dist(i),eval(j),weight
            enddo
            deallocate(eval,evec) !allocate at beginning
        enddo
        close(io)

        !write help function
        open(newunit=io,file=fname//'.py')
        write(io,'(A)') 'import numpy as np'
        write(io,'(A)') 'import matplotlib.pyplot as plt'
        write(io,'(A)') 'from matplotlib import colors'
        write(io,'(A)') ''
        write(io,'(A)') 'c="red"'
        write(io,'(A)') 'cmap = colors.LinearSegmentedColormap.from_list("incr_alpha", [(0, (*colors.to_rgb(c),0)), (1, c)])'
        write(io,'(A)') 'dat=np.loadtxt("highs_plot_sc.dat")'
        write(io,'(A)') 'fig,ax =plt.subplots()'
        write(io,'(A)') 'im=ax.scatter(dat[:,0],dat[:,1],c=dat[:,2],s=1.0,edgecolor=None,cmap=cmap,vmin=0.0,vmax=1.0,zorder=10)'
        write(io,'(A)') 'ax.set_xlim(np.amin(dat[:,0]),np.amax(dat[:,0]))'
        write(io,'(A)') 'ax.set_xticks(['
        write(io,'(E16.8,",")') (highs_dist(i) ,i=1,size(highs_dist))
        write(io,'(A)') '])'
        write(io,'(A)') 'ax.set_xticklabels(['
        write(io,'(A,",")') ('"label"' ,i=1,size(highs_dist))
        write(io,'(A)') '])'
        write(io,'(A)') 'ax.set_ylabel(r"E-E$_F$ [eV]")'
        write(io,'(A)') 'ax.grid(axis="x",color="black",zorder=1)'
        write(io,'(A)') 'cbar=plt.colorbar(im,ax=ax)'
        write(io,'(A)') 'cbar.set_label("electron weight")'
        write(io,'(A)') 'plt.show()'
        write(io,'(3A)') 'fig.savefig("',fname,'.png")'
        close(io)
    end subroutine


    subroutine plot_highsym_kpts(Hk_inp,h_io)
        type(Hk_inp_t),intent(in)               :: Hk_inp
        type(parameters_TB_IO_H),intent(in)     :: h_io

        real(8),allocatable :: eval(:)
        character(len=*),parameter  ::  fname="highs_plot"
        integer             ::  io,i,j

        if(.not. allocated(kpts)) return
        write(*,'(A,I6,A)') "Calculate ",N_kpts,' kpoints on the high symmetry path'
        !get and print data
        open(newunit=io,file=fname//'.dat')
        do i=1,N_kpts
            Call Hk_eval(Hk_inp,kpts(:,i),h_io,eval) 
            write(io,'(2E16.8)') (kpts_dist(i),eval(j),j=1,size(eval))
            deallocate(eval)   !allocate at beginning
        enddo
        close(io)

        !write help function
        open(newunit=io,file=fname//'.gnu')
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
        close(io)
    end subroutine

    subroutine plot_highsym_kpts_proj(Hk_inp,h_io)
        type(Hk_inp_t),intent(in)               :: Hk_inp
        type(parameters_TB_IO_H),intent(in)     :: h_io

        real(8),allocatable     :: eval(:)
        complex(8),allocatable  :: evec(:,:)
        character(len=*),parameter  ::  fname="highs_plot_proj"
        integer             ::  io,i,j,l
        real(8),allocatable ::  weight(:)
        integer             ::  norb,ncell,dimH
        character(len=8)    ::  Norb_char

        if(.not. allocated(kpts)) return
        write(*,'(A,I6,A)') "Calculate ",N_kpts,' kpoints on the high symmetry path'
        norb=h_io%norb
        ncell=h_io%ncell
        dimH=h_io%dimH
        allocate(weight(norb))
        write(Norb_char,'(I8)') norb
        !get and print data
        open(newunit=io,file=fname//'.dat')
        do i=1,N_kpts
            Call Hk_evec(Hk_inp,kpts(:,i),h_io,eval,evec) 
            do j=1,size(eval)
                do l=1,norb !not very elegant...
!                    write(*,*) evec(l::ncell,j)
                    weight(l)=real(dot_product(evec(l::norb,j),evec(l::norb,j)),8)
!                    write(*,*) weight(l)
!                    STOP
                enddo
                write(io,'(2E16.8X1'//norb_char//'E16.8)') kpts_dist(i),eval(j),weight
            enddo
            deallocate(eval,evec)   !allocate at beginning
        enddo
        close(io)
    end subroutine


    subroutine highsym_clear_path()
        deallocate(kpts)
        deallocate(kpts_dist)
        deallocate(highs_dist)
    end subroutine 

    subroutine set_highs_path(lat,highs)
        use m_constants, only : pi
        type(lattice), intent(in)                   :: lat
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
            k_highs(:,i)=matmul(k_highs(:,i),lat%a_sc_inv)
        enddo
        k_highs=k_highs*2.0d0*pi    !no 2*pi factor in a_sc_inv

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
