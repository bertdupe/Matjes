module m_fermi_dos
!module which contains routines that allow to calculate number of states at the fermi-energy resolved in k-space
!All energies are broadened with a gaussian distribution and the distributions from all relevant energy levels  at the zero energy are added up (fermi_dos_nc)
!alternatively, additionally the projections on list of basic orbitals (fermi_dos_proj_nc) or on all orbitals (fermi_dos_projall_nc) can be calculated 
!So far this only works for the normal conducting case (_nc), because at E=0 the superconducting case doesn't exactly makes sense

use m_derived_types, only : t_cell,lattice
use m_get_position, only: calculate_distances,get_position
use m_highsym, only: plot_highsym_kpts,set_highs_path
use m_TB_types
use m_tb_k_public       !mode that contains the more efficient TB k-space type
private
public fermi_dos_nc, fermi_dos_proj_nc, fermi_dos_projall_nc
real(8),parameter       ::  dist_inc=8.0d0  !how many sigma away from my the energy entries are still considered

contains
subroutine fermi_dos_nc(HK,h_io,lat,dos_io,work)
    use m_derived_types, only: k_grid_t
    class(H_k_base),intent(inout)           :: Hk
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_dos), intent(in)  :: dos_io
    type(work_ham),intent(inout)            :: work

    type(k_grid_t)              :: k_grid

    real(8),allocatable         :: dos_weight(:)    !resulting weight at each k
    integer                     :: Nk               !number of k-points
    real(8),allocatable         :: eval(:)          !eigenvalues
    real(8)                     :: sigma            !gauss smearing
    integer                     :: ibnd(2), Nentry  !variables controlling which leveln are considered
    integer                     :: Nin,Nout         !maximal and output number of eigenvalues
    integer                     :: i,ik, io
    real(8)                     :: k(3)             !temporary k-point

    Call k_grid%set(lat%a_sc_inv,dos_io%kgrid)
    Nk=k_grid%get_Nk()

    sigma=dos_io%sigma
    Nin=Hk%get_size_eval()
    allocate(eval(Nin),source=0.0d0)
    allocate(dos_weight(Nk),source=0.d0)
    do ik=1,Nk
        k=k_grid%get_K(ik)
        Call Hk%get_eval(k,Nin,eval,Nout,work) 
        Call get_bnd(Nout,eval,sigma,ibnd,Nentry)
        if(Nentry<1) cycle
        Call get_gauss(eval(ibnd(1):ibnd(2)),0.0d0,sigma)
        dos_weight(ik)=sum(eval(ibnd(1):ibnd(2)))
    enddo
    deallocate(eval)

    open(newunit=io,file='fermidos.dat')
    do ik=1,Nk
        k=k_grid%get_K(ik)
        write(io,*) k, dos_weight(ik)
    enddo
    close(io)
end subroutine

subroutine fermi_dos_projall_nc(HK,h_io,lat,dos_io,work)
    use m_derived_types, only: k_grid_t
    class(H_k_base),intent(inout)           :: Hk
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_dos), intent(in)  :: dos_io
    type(work_ham),intent(inout)            :: work

    type(k_grid_t)                          :: k_grid

    real(8),allocatable         :: dos_weight(:,:)
    integer                     :: Nk
    real(8)   ,allocatable      :: eval(:)
    complex(8),allocatable      :: evec(:,:)
    real(8)                     :: sigma
    integer                     :: ibnd(2)
    integer                     :: Nin,Nout         !maximal and output number of eigenvalues
    integer                     :: i,ik, io
    real(8)                     :: k(3)

    integer :: dimH,ndim,Nentry
    real(8),allocatable ::  pref(:,:)

    Call k_grid%set(lat%a_sc_inv,dos_io%kgrid)
    Nk=k_grid%get_Nk()

    sigma=dos_io%sigma
    Nin=Hk%get_size_eval()
    dimH=Hk%get_dimH()
    allocate(eval(Nin),source=0.0d0)
    allocate(evec(dimH,Nin),source=(0.0d0,0.0d0))
    allocate(dos_weight(dimH,Nk),source=0.d0)
    allocate(pref(dimH,dimH))

    do ik=1,Nk
        k=k_grid%get_K(ik)
        Call Hk%get_evec(k,Nin,eval,evec,Nout,work) 
        Call get_bnd(Nout,eval,sigma,ibnd,Nentry)
        if(Nentry<1) cycle
        do i=1,dimH
            pref(1:Nentry,i)=conjg(evec(i,ibnd(1):ibnd(2)))*evec(i,ibnd(1):ibnd(2))
        enddo
        Call get_gauss(eval(ibnd(1):ibnd(2)),0.0d0,sigma)
        dos_weight(:,ik)=matmul(eval(ibnd(1):ibnd(2)),pref(1:Nentry,:))
    enddo
    deallocate(eval,evec)

    !this separation might result in too large dos_weight arrays -> possible to put the IO directly above
    open(newunit=io,file='fermidos_allproj.dat')
    do ik=1,Nk
        k=k_grid%get_K(ik)
        write(io,*) k, dos_weight(:,ik)
    enddo
    close(io)
end subroutine


subroutine fermi_dos_proj_nc(HK,h_io,lat,dos_io,work)
    use m_derived_types, only: k_grid_t
    class(H_k_base),intent(inout)           :: Hk
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_dos), intent(in)  :: dos_io
    type(work_ham),intent(inout)            :: work

    type(k_grid_t)              :: k_grid

    real(8),allocatable         :: dos_weight(:,:)
    integer                     :: Nk
    real(8)   ,allocatable      :: eval(:)
    complex(8),allocatable      :: evec(:,:)
    real(8)                     :: sigma
    integer                     :: ibnd(2)
    integer                     :: Nin,Nout         !maximal and output number of eigenvalues
    integer                     :: i,ik, io
    real(8)                     :: k(3)

    integer,allocatable         :: proj_ind(:)
    integer                     :: dimH,Nproj,Nentry
    real(8),allocatable         :: pref(:,:)

    proj_ind=dos_io%fermi_orb
    Call k_grid%set(lat%a_sc_inv,dos_io%kgrid)
    Nk=k_grid%get_Nk()

    Nproj=size(proj_ind)
    sigma=dos_io%sigma
    Nin=Hk%get_size_eval()
    dimH=Hk%get_dimH()
    allocate(eval(Nin),source=0.0d0)
    allocate(evec(dimH,Nin),source=(0.0d0,0.0d0))
    allocate(dos_weight(Nproj,Nk),source=0.d0)
    allocate(pref(dimH,Nproj))
    do ik=1,Nk
        k=k_grid%get_K(ik)
        Call Hk%get_evec(k,Nin,eval,evec,Nout,work) 
        Call get_bnd(Nout,eval,sigma,ibnd,Nentry)
        if(Nentry<1) cycle
        do i=1,Nproj
            pref(1:Nentry,i)=conjg(evec(proj_ind(i),ibnd(1):ibnd(2)))*evec(proj_ind(i),ibnd(1):ibnd(2))
        enddo
        Call get_gauss(eval(ibnd(1):ibnd(2)),0.0d0,sigma)
        dos_weight(:,ik)=matmul(eval(ibnd(1):ibnd(2)),pref(1:Nentry,:))
    enddo
    deallocate(eval,evec)

    open(newunit=io,file='fermidos_proj.dat')
    do ik=1,Nk
        k=k_grid%get_K(ik)
        write(io,*) k, dos_weight(:,ik)
    enddo
    close(io)
end subroutine

subroutine get_bnd(Nev,eval,sigma,ibnd,Nentry)
    integer,intent(in)  :: Nev
    real(8),intent(in)  :: eval(Nev)
    real(8),intent(in)  :: sigma
    integer,intent(out) :: ibnd(2),Nentry
    integer ::  i
    
    ibnd=[Nev+1,0]
    do i=1,Nev
        if(eval(i)+dist_inc*sigma>0.d0)then
            ibnd(1)=i 
            exit
        endif
    enddo
    do i=Nev,1,-1
        if(eval(i)-dist_inc*sigma<0.d0)then
            ibnd(2)=i 
            exit
        endif
    enddo
    Nentry=ibnd(2)-ibnd(1)+1
end subroutine

subroutine get_gauss(val,mu,sigma)
    use m_constants, only : pi
    real(8),intent(inout)  ::  val(:)
    real(8),intent(in)     ::  sigma
    real(8),intent(in)     ::  mu
    val=(val-mu)**2
    val=-val*0.5d0/sigma/sigma
    val=exp(val)
    val=val/sqrt(2.0d0*pi)/sigma
end subroutine

end module 
