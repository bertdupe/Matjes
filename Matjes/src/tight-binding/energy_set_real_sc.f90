module m_energy_set_real_sc
use m_energy_commons, only : energy
use m_basic_types, only : vec_point

private
public set_Hr,Hr_eigval,Hr_eigvec,get_Hr

!large electronic Hamiltonian
complex(8),allocatable  ::  Hr(:,:)

!setup is now a bit stupid in that Hr has to be saved and copied for the lapack routines, which doubles their memory at some points...


!basis of Hamiltonian is as follows
!  inner index -> outer index
!  spin -> orbital ->lattice-site -> electron/hole
! basis is c_up, c_down,c_up^+,c_down^+ for the right basis part

contains

    subroutine get_Hr(dimH,Hr_out)
        integer,intent(in)          ::  dimH
        complex(8),intent(out)      ::  Hr_out(dimH,dimH)
        
        if(.not.allocated(Hr)) STOP "Hr is not allocated but get_Hr is called"
        if(dimH/=size(Hr,1)) STOP "dimensions of Hr seems to be wrong for getting Hr"
        Hr_out=Hr
    end subroutine 

    subroutine set_Hr(dimH,Tb_ext,mode_mag)
        !extract the real space Hamiltonian Hr from the electronic part in energy
        use m_energy_set_real, only: get_Hr_nc=>get_Hr, set_Hr_nc=>set_Hr
        use m_tb_params, only : TB_params
        integer,intent(in)           :: dimH
        integer,intent(in)           :: TB_ext(2)
        type(vec_point),intent(in)   :: mode_mag(:)
        complex(8)                   :: Hr_nc(dimH/2,dimH/2)

        integer                 ::  dimH_nc,dim_mode

        if(.not. allocated(Hr))then
           allocate(Hr(dimH,dimH))
        endif
        if(size(Hr,1)/=dimH.or.size(Hr,2)/=dimH) STOP "Hr has wrong size"  !could easily reallocate, but this should never happen, I guess
        Hr=cmplx(0.0d0,0.0d0, kind=8)
        dimH_nc=dimH/2
        Call set_Hr_nc(dimH_nc,tb_ext,mode_mag)
        Call get_Hr_nc(dimH_nc,Hr_nc)
        Hr(1:dimH_nc,1:dimH_nc)=Hr_nc
        Hr(dimH_nc+1:dimH,dimH_nc+1:dimH)=-Hr_nc
        dim_mode=Tb_ext(2)-Tb_ext(1)+1
        Call set_delta(TB_params%io_H%delta,dim_mode)
    end subroutine 

    subroutine set_delta(delta,dim_mode)
        complex(8)          ::  delta(:)
        integer,intent(in)  ::  dim_mode
        integer             ::  n_cells

        integer         ::  dimH_nc,dimH_mag,dim_mode_red
        integer         ::  i_cell,i_orb
        integer         ::  i_up,i_dn,i_up_dg,i_dn_dg

        N_cells = size(energy%line,2)
        dim_mode_red=dim_mode/2
        dimH_mag=dim_mode_red*N_cells
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

    subroutine Hr_eigval(dimH,eigval)
        integer,intent(in)          ::  dimH
        real(8),intent(out)         ::  eigval(dimH)

        complex(8)                  :: H_loc(dimH,dimH)
        complex(kind=8)             :: init_WORK(2*dimH)
        real(kind=8)                :: RWORK(3*dimH-2)
        integer                     :: info,l_work
        complex(kind=8),allocatable :: WORK(:)
        
        if(.not.allocated(Hr)) STOP "Hr is not allocated but Hr_eigval is called"
        if(dimH/=size(Hr,1)) STOP "dimensions of Hr seems to be wrong evaluating the eigenvalues"
        H_loc=Hr
        call ZHEEV( 'V', 'U', dimH, H_loc, dimH, eigval, init_WORK, -1, RWORK, INFO )
        l_work=int(init_work(1))
        allocate(work(l_work),source=cmplx(0.0d0,0.0d0,8))
        call ZHEEV( 'N', 'U', dimH, H_loc, dimH, eigval, WORK, l_work, RWORK, INFO )

    end subroutine
#if 0
    subroutine Hr_eigvec(dimH,eigvec,eigval)
        integer,intent(in)          ::  dimH
        complex(8),intent(out)      ::  eigvec(dimH,dimH)
        real(8),intent(out)         ::  eigval(dimH)

        complex(kind=8)             :: init_WORK(2*dimH)
        real(kind=8)                :: RWORK(3*dimH-2)
        integer                     :: info,l_work
        complex(kind=8),allocatable :: WORK(:)

        if(.not.allocated(Hr)) STOP "Hr is not allocated but Hr_eigvec is called"
        if(dimH/=size(Hr,1)) STOP "dimensions of Hr seems to be wrong evaluating the eigenvectors"
        eigvec=Hr
        eigval=0.0d0
        call ZHEEV( 'V', 'U', dimH, eigvec, dimH, eigval, init_WORK, -1, RWORK, INFO )
        l_work=int(init_work(1))
        allocate(work(l_work),source=cmplx(0.0d0,0.0d0,8))
        call ZHEEV( 'V', 'U', dimH, eigvec, dimH, eigval, WORK, l_work, RWORK, INFO )
    end subroutine 


#else
    subroutine Hr_eigvec(dimH,eigvec,eigval)
        integer,intent(in)          ::  dimH
        complex(8),intent(out)      ::  eigvec(dimH,dimH)
        real(8),intent(out)         ::  eigval(dimH)

        integer                     :: info,lwork
        integer                     :: lrwork,liwork
        integer,allocatable         :: iwork(:)
        complex(kind=8),allocatable :: WORK(:)
        real(8),allocatable         :: RWORK(:)

        integer                     :: tmp_iwork(1)
        real(8)                     :: tmp_rwork(1)
        complex(8)                  :: tmp_work(1)

        if(.not.allocated(Hr)) STOP "Hr is not allocated but Hr_eigvec is called"
        if(dimH/=size(Hr,1)) STOP "dimensions of Hr seems to be wrong evaluating the eigenvectors"
        eigvec=Hr
        eigval=0.0d0
        call ZHEEVD( 'V', 'U', dimH, eigvec, dimH, eigval, tmp_WORK, -1, tmp_RWORK, -1,tmp_IWORK,-1,INFO )
        lwork=int(tmp_work(1))
        LIWORK=tmp_IWORK(1)
        LRWORK=int(tmp_rwork(1))
        allocate(work(lwork),source=cmplx(0.0d0,0.0d0,8))
        allocate(iwork(liwork),source=0)
        allocate(rwork(lrwork),source=0.0d0)
        call ZHEEVD( 'V', 'U', dimH, eigvec, dimH, eigval, WORK, lwork, RWORK, LRWORK,IWORK,LIWORK,INFO )
    end subroutine 
#endif
end module
