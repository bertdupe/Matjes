module m_fermi_dos
use m_derived_types, only : t_cell,lattice
use m_get_position, only: calculate_distances,get_position
use m_highsym, only: plot_highsym_kpts,set_highs_path
use m_TB_types
use m_Hk
private
public fermi_dos_nc, fermi_dos_proj_nc, fermi_dos_projall_nc
real(8),parameter       ::  dist_inc=8.0d0  !how many sigma away from my the energy entries are still considered

contains
subroutine fermi_dos_nc(HK_inp,h_io,lat,dos_io)
    use m_derived_types, only: k_grid_t
    type(Hk_inp_t),intent(in)               :: Hk_inp
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_dos), intent(in)  :: dos_io

    type(k_grid_t)                          :: k_grid

    real(8),allocatable                     :: dos_weight(:)
    integer                                 :: Nk
    real(8)   ,allocatable                  :: eval(:)
    complex(8),allocatable                  :: evec(:,:)
    real(8)                                 :: sigma
    integer                                 :: ibnd(2), Nentry
    integer ::  i,ik, io
    real(8) :: k(3)

    Call k_grid%set(lat%a_sc_inv,dos_io%kgrid)
    Nk=k_grid%get_Nk()
    allocate(dos_weight(Nk),source=0.d0)
    sigma=dos_io%sigma
    do ik=1,Nk
        k=k_grid%get_K(ik)
        Call Hk_eval(Hk_inp,k,h_io,eval) 
        Call get_bnd(eval,sigma,ibnd,Nentry)
        if(Nentry<1) cycle
        Call get_gauss(eval(ibnd(1):ibnd(2)),0.0d0,sigma)
        dos_weight(ik)=sum(eval(ibnd(1):ibnd(2)))
        deallocate(eval)
    enddo

    open(newunit=io,file='fermidos.dat')
    do ik=1,Nk
        k=k_grid%get_K(ik)
        write(io,*) k, dos_weight(ik)
    enddo
    close(io)
end subroutine

subroutine fermi_dos_projall_nc(HK_inp,h_io,lat,dos_io)
    use m_derived_types, only: k_grid_t
    type(Hk_inp_t),intent(in)               :: Hk_inp
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_dos), intent(in)  :: dos_io

    type(k_grid_t)                          :: k_grid

    real(8),allocatable                     :: dos_weight(:,:)
    integer                                 :: Nk
    real(8)   ,allocatable                  :: eval(:)
    complex(8),allocatable                  :: evec(:,:)
    real(8)                                 :: sigma
    integer                                 :: ibnd(2)
    integer ::  i,ik, io
    real(8) :: k(3)

    integer :: dimH,ndim,Nentry
    real(8),allocatable ::  pref(:,:)

    Call k_grid%set(lat%a_sc_inv,dos_io%kgrid)
    Nk=k_grid%get_Nk()
    sigma=dos_io%sigma
    dimH=Hk_inp%H(1)%dimH
    ndim=Hk_inp%H(1)%ndim
    allocate(dos_weight(dimH,Nk),source=0.d0)
    allocate(pref(dimH,dimH))
    do ik=1,Nk
        k=k_grid%get_K(ik)
        Call Hk_evec(Hk_inp,k,h_io,eval,evec) 
        Call get_bnd(eval,sigma,ibnd,Nentry)
        if(Nentry<1) cycle
        do i=1,dimH
            pref(1:Nentry,i)=conjg(evec(i,ibnd(1):ibnd(2)))*evec(i,ibnd(1):ibnd(2))
        enddo
        Call get_gauss(eval(ibnd(1):ibnd(2)),0.0d0,sigma)
        dos_weight(:,ik)=matmul(eval(ibnd(1):ibnd(2)),pref(1:Nentry,:))
        deallocate(eval)
    enddo

    open(newunit=io,file='fermidos_allproj.dat')
    do ik=1,Nk
        k=k_grid%get_K(ik)
        write(io,*) k, dos_weight(:,ik)
    enddo
    close(io)
end subroutine


subroutine fermi_dos_proj_nc(HK_inp,h_io,lat,dos_io)
    use m_derived_types, only: k_grid_t
    type(Hk_inp_t),intent(in)               :: Hk_inp
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_dos), intent(in)  :: dos_io

    type(k_grid_t)                          :: k_grid

    real(8),allocatable                     :: dos_weight(:,:)
    integer                                 :: Nk
    real(8)   ,allocatable                  :: eval(:)
    complex(8),allocatable                  :: evec(:,:)
    real(8)                                 :: sigma
    integer                                 :: ibnd(2)
    integer ::  i,ik, io
    real(8) :: k(3)

    integer,allocatable   :: proj_ind(:)
    integer :: dimH,ndim,Nproj,Nentry
    real(8),allocatable ::  pref(:,:)

    proj_ind=dos_io%fermi_orb
    Call k_grid%set(lat%a_sc_inv,dos_io%kgrid)
    Nk=k_grid%get_Nk()
    sigma=dos_io%sigma
    dimH=Hk_inp%H(1)%dimH
    ndim=Hk_inp%H(1)%ndim
    if(dimH/=ndim) ERROR STOP "implement including orbitals from different unit-cells with another loop"
    Nproj=size(proj_ind)
    allocate(dos_weight(Nproj,Nk),source=0.d0)
    allocate(pref(dimH,Nproj))
    do ik=1,Nk
        k=k_grid%get_K(ik)
        Call Hk_evec(Hk_inp,k,h_io,eval,evec) 
        Call get_bnd(eval,sigma,ibnd,Nentry)
        if(Nentry<1) cycle
        do i=1,Nproj
            pref(1:Nentry,i)=conjg(evec(proj_ind(i),ibnd(1):ibnd(2)))*evec(proj_ind(i),ibnd(1):ibnd(2))
        enddo
        Call get_gauss(eval(ibnd(1):ibnd(2)),0.0d0,sigma)
        dos_weight(:,ik)=matmul(eval(ibnd(1):ibnd(2)),pref(1:Nentry,:))
        deallocate(eval)
    enddo

    open(newunit=io,file='fermidos_proj.dat')
    do ik=1,Nk
        k=k_grid%get_K(ik)
        write(io,*) k, dos_weight(:,ik)
    enddo
    close(io)
end subroutine

subroutine get_bnd(eval,sigma,ibnd,Nentry)
    real(8),intent(in)  :: eval(:)
    real(8),intent(in)  :: sigma
    integer,intent(out) :: ibnd(2),Nentry
    integer ::  i
    
    ibnd=[size(eval)+1,0]
    do i=1,size(eval)
        if(eval(i)+dist_inc*sigma>0.d0)then
            ibnd(1)=i 
            exit
        endif
    enddo
    do i=size(eval),1,-1
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
