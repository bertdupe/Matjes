module m_occupation
use m_fermi,only: fermi_distrib
use m_tb_types
use m_io_files_utils, only: close_file,open_file_write
use m_distribution, only: int_distrib
private
public calc_occupation,calc_occupation_sc




contains 
    
    subroutine calc_occupation(h_par,eigvec,eigval,E_f_in,kt,fname,dist_ptr)
        type(parameters_TB_Hsolve),intent(in)     ::  h_par
        real(8),intent(in)      ::  eigval(:),E_f_in,kt
        complex(8),intent(in)   ::  eigvec(:,:)
        character(len=*),intent(in) ::  fname
		procedure(int_distrib),pointer,intent(in)	:: dist_ptr

        real(8),allocatable     ::  occ(:)
        integer                 ::  n_occ 
        real(8)                 ::  E_F


        n_occ=h_par%nspin*h_par%ncell*h_par%norb
        allocate(occ(n_occ),source=0.0d0)
        E_f=E_F_in

        !calc occupations for normal conducting or sc case
        if(h_par%nsc==2)then
            !!Fermi energy doesn't really make much sense here
            !if(E_f<=1d-6.or.E_f>=1d-6)then
            !    write(*,*) "Warning, setting E_F=0 for SC"
            !endif
            !E_f=0.0d0
            Call calc_occupation_sc(eigvec,eigval,E_f,kt,occ,dist_ptr)
        elseif(h_par%nsc==1)then
            Call calc_occupation_nc(eigvec,eigval,E_f,kt,occ,dist_ptr)
        else
            write(*,*) 'h_par%nsc=',h_par%nsc
            STOP 'unexpected value for h_par%nsc'
        endif

        !write output
        Call write_output(h_par,fname,occ)

    end subroutine

    subroutine write_output(h_par,fname,occ)
        !sum up orbitals of each cell and print out each spin-channel 
        type(parameters_TB_Hsolve),intent(in)     ::  h_par
        character(len=*),intent(in) ::  fname
        real(8),intent(in)          ::  occ(:)
        
        integer                     :: i,io
        integer                     :: ist,ien !start/end
        integer                     :: per_site

        per_site=h_par%nspin*h_par%norb
        io=open_file_write(fname)
        if(h_par%nspin==2)then
            do i=1,h_par%ncell
                ist=1+(i-1)*per_site
                ien=ist+per_site-1
                write(io,'(2E16.8)') sum(occ(ist:ien:2)),sum(occ(ist+1:ien:2))
            enddo
        else
            do i=1,h_par%ncell
                ist=1+(i-1)*per_site
                ien=ist+per_site-1
                write(io,'(E16.8)') sum(occ(ist:ien))
            enddo
        endif
        call close_file(fname,io)
    end subroutine 

    subroutine calc_occupation_nc(eigvec,eigval,E_f,kt,occ,dist_ptr)
        real(8),intent(in)      ::  eigval(:),E_f,kt
        complex(8),intent(in)   ::  eigvec(:,:)
        real(8),intent(out)     ::  occ(:)
		procedure(int_distrib),pointer,intent(in)	:: dist_ptr

        integer                 ::  i
        real(8)                 ::  fd

        occ=0.0d0
        do i=1,size(eigval)
            fd=dist_ptr(E_f,eigval(i),kt)
            occ=occ+fd*real(conjg(eigvec(:,i))*eigvec(:,i),kind=8)
        enddo
    end subroutine


    subroutine calc_occupation_sc(eigvec,eigval,E_f,kt,occ,dist_ptr)
        real(8),intent(in)      ::   eigval(:),E_f,kt
        complex(8),intent(in)   ::   eigvec(:,:)
        real(8),intent(out)     ::   occ(:)
		procedure(int_distrib),pointer,intent(in)	:: dist_ptr

        integer                 ::  i
        integer                 ::  n_state
        integer                 ::  dimH,half
        real(8)                 ::  fd

        dimH=size(eigvec,1)
        half=dimH/2
        if(size(occ)/=half) STOP "sizes in calc_occupation_sc not matching"
        n_state=size(eigval)
        do i=1,size(eigval)
            if(eigval(i)>0.0d0)then
                n_state=i-1 
                exit
            endif
        enddo
        if(n_state<1) STOP "n_state<1 in calc_occupation_sc, eigval not sorted or all possitive?"
        do i=1,n_state
            fd=dist_ptr(E_f,eigval(i),kt)
            occ=occ+fd*real(conjg(eigvec(1:half,i))*eigvec(1:half,i),kind=8)
            fd=dist_ptr(E_f,-eigval(i),kt)
            occ=occ+fd*real(conjg(eigvec(half+1:dimH,i))*eigvec(half+1:dimH,i),kind=8)
        enddo

    end subroutine


end module
