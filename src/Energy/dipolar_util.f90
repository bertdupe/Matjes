module m_dipolar_util
public
!Matjes length scale is nm, magnetic moments are in bohr magneton mu_b, energy shall be in eV
!mu_0*mu_b^2/(nm^3)  * J_to_eV
!mu_0=1.25663706212d-6  kg*m/s^2/A^2
!mu_b=9.2740100783d-24  A*m^2
!nm=1.0d-9  m
!J_to_eV=1/1.602176634d-19 

real(8),parameter   ::  dip_pref=6.745817653234234066975d-7

contains

subroutine get_supercell_vec(supercell_vec,lat,period)
    use m_derived_types, only: lattice
    real(8),intent(inout),allocatable   :: supercell_vec(:,:)
    type(lattice),intent(in)            :: lat
    integer,intent(in)                  :: period(3)

    integer             ::  bnd_ext(2,3)
    integer             ::  Nrep(3)
    integer             ::  Ntot
    integer             ::  i,ii, i3,i2,i1
    real(8)             ::  vec(3,3)

    if(allocated(supercell_vec)) deallocate(supercell_vec)
    bnd_ext=0
    do i=1,3
        if(lat%periodic(i)) bnd_ext(:,i)=[-period(i),period(i)]
        Nrep(i)=bnd_ext(2,i)-bnd_ext(1,i)+1
    enddo
    Ntot=product(Nrep)
    allocate(supercell_vec(3,Ntot))
    ii=0
    do i3=1,Nrep(3)
        vec(:,3)=lat%a_sc(3,:)*real(bnd_ext(1,3)+i3-1,8)
        do i2=1,Nrep(2)
            vec(:,2)=lat%a_sc(2,:)*real(bnd_ext(1,2)+i2-1,8)
            do i1=1,Nrep(1)
                vec(:,1)=lat%a_sc(1,:)*real(bnd_ext(1,1)+i1-1,8)
                ii=ii+1
                supercell_vec(:,ii)=sum(vec,2)
            enddo
        enddo
    enddo
    if(ii/=Ntot) ERROR STOP "unexpected size for supercell vectors in dipolar calculation"
end subroutine

end module
