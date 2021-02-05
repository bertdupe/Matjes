module m_fermi_dos
use m_derived_types, only : t_cell,lattice
use m_get_position, only: calculate_distances,get_position
use m_highsym, only: plot_highsym_kpts,set_highs_path
use m_TB_types
use m_init_Hk
private
public fermi_dos_nc
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
    real(8),allocatable                     :: eval(:)
    real(8)                                 :: sigma
    integer                                 :: ibnd(2)
    integer ::  i,ik, io
    real(8) :: k(3)

    Call k_grid%set(lat%a_sc_inv,dos_io%kgrid)
    Nk=k_grid%get_Nk()
    allocate(dos_weight(Nk),source=0.d0)
    sigma=dos_io%sigma
    do ik=1,Nk
!        if(.false.) write(output_unit,'(2(AI6))') 'start fermi dos k', ik,' of',Nk
        k=k_grid%get_K(ik)
        Call Hk_eval(Hk_inp,k,h_io,eval) 
        ibnd=[1,size(eval)]
        do i=1,size(eval)
            if(eval(i)+dist_inc*sigma>0.d0)then
                ibnd(1)=max(1,i-1) 
                exit
            endif
        enddo
        do i=ibnd(1),size(eval)
            if(eval(i)-dist_inc*sigma>0.d0)then
                ibnd(2)=max(ibnd(1),i-1) 
                exit
            endif
        enddo

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
