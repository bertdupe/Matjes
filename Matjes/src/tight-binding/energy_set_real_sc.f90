module m_energy_set_real_sc
use m_energy_commons, only : energy
use m_basic_types, only : vec_point
use m_tb_types
use m_energy_set_real, only: set_Hr_dense_nc
implicit none

private
public set_Hr_dense_sc


!basis of Hamiltonian is as follows
!  inner index -> outer index
!  spin -> orbital ->lattice-site -> electron/hole
! basis is c_up, c_down,c_up^+,c_down^+ for the right basis part

contains


    subroutine set_Hr_dense_sc(h_par,mode_mag,Hr)
        !extract the real space Hamiltonian Hr from the electronic part in energy
        use m_tb_params, only : TB_params
        type(parameters_TB_Hsolve),intent(in)     ::  h_par
        type(vec_point),intent(in)   :: mode_mag(:)
        complex(8),allocatable,intent(inout)    ::  Hr(:,:)
        complex(8),allocatable                  :: Hr_nc(:,:)

        type(parameters_TB_Hsolve)     ::  h_par_nc
        integer                 ::  dimH_nc,dim_mode

        if(.not. allocated(Hr))then
           allocate(Hr(h_par%dimH,h_par%dimH))
        endif
        if(size(Hr,1)/=h_par%dimH.or.size(Hr,2)/=h_par%dimH) STOP "Hr has wrong size"  !could easily reallocate, but this should never happen, I guess

        !get Hamiltonian without superconductivity
        h_par_nc=h_par
        h_par_nc%nsc=1
        Call h_par_nc%upd()
        Call set_Hr_dense_nc(h_par_nc,mode_mag,Hr_nc)

        !fill local Hamiltonian and add SC_delta
        Hr=cmplx(0.0d0,0.0d0, kind=8)
        Hr(1:h_par_nc%dimH,1:h_par_nc%dimH)=Hr_nc
        Hr(h_par_nc%dimH+1:h_par%dimH,h_par_nc%dimH+1:h_par%dimH)=-Hr_nc
        dim_mode=h_par%pos_ext(2)-h_par%pos_ext(1)+1
        Call set_delta(TB_params%io_H%delta,dim_mode,Hr)

        if(h_par%rearrange) Call rearange_H(h_par%dimH,Hr)
    end subroutine 

    subroutine set_delta(delta,dim_mode,Hr)
        !use h_par as well
        complex(8),intent(in)       ::  delta(:)
        complex(8),intent(inout)    ::  Hr(:,:) !set checks here for dimension
        integer,intent(in)          ::  dim_mode
        integer                     ::  n_cells

        integer         ::  dimH_nc,dim_mode_red
        integer         ::  i_cell,i_orb
        integer         ::  i_up,i_dn,i_up_dg,i_dn_dg

        N_cells = size(energy%line,2)
        dim_mode_red=dim_mode/2
        dimH_nc=size(Hr,1)/2
        do i_cell=1,N_cells
            do i_orb=1,dim_mode_red
                i_up=2*(i_cell-1)*i_orb+2*i_orb-1
                i_dn=i_up+1
                i_up_dg=i_up+dimH_nc
                i_dn_dg=i_dn+dimH_nc
                Hr(i_up_dg,i_dn)=Hr(i_up_dg,i_dn)-conjg(delta(i_orb))
                hr(i_dn_dg,i_up)=Hr(i_dn_dg,i_up)+conjg(delta(i_orb))
                Hr(i_up,i_dn_dg)=Hr(i_up,i_dn_dg)+delta(i_orb)
                Hr(i_dn,i_up_dg)=Hr(i_dn,i_up_dg)-delta(i_orb)
            enddo
        enddo
    end subroutine

    subroutine rearange_H(dimH,Hr)
        !super stupid implementation to swap
        integer,intent(in)              ::  dimH
        complex(8),intent(inout)        ::  Hr(dimH,dimH)
        complex(8)                      ::  tmp(dimH,dimH)

        integer                         ::  i
        integer                         ::  ind1,ind2
    
        tmp=Hr
        do i=1,dimH/2
            Hr(:,2*i-1)=tmp(:,i)
            Hr(:,2*i)=tmp(:,i+dimH/2)
        enddo
        do i=1,dimH/2
            Hr(2*i-1,:)=tmp(i,:)
            Hr(2*i,:)=tmp(i+dimH/2,:)
        enddo
    end subroutine

end module
