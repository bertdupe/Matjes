module m_init_Hr
use m_derived_types, only: lattice
use m_H_tb_public
use m_tb_types ,only: parameters_TB_Hsolve,parameters_TB_IO_H
implicit none
private
public get_all_Hr

contains

subroutine get_all_Hr(lat,h_io,H)
    type(lattice),intent(in)                :: lat
    type(parameters_TB_IO_H),intent(in)     :: h_io
    class(H_tb),allocatable,intent(out)     :: H

    Call set_H(H,h_io)
    if(allocated(h_io%hop))then
        Call get_Hr_hop(lat,h_io,H)
        if(allocated(h_io%Jsd)) Call get_Hr_Jsd(lat,h_io,H)
        if(allocated(h_io%del)) Call get_Hr_delta(lat,h_io,H)
    endif

end subroutine


subroutine get_Hr_hop(lat,h_io,H)
    !extract the real space Hamiltonian Hr from the electronic part in energy
    !TODO, make neigh%get more efficient by first sorting though h_io%hop to get all required connections for an atom type pair at once
    use m_neighbor_type, only: neighbors
    type(lattice),intent(in)                :: lat
    type(parameters_TB_IO_H),intent(in)     :: h_io
    class(H_tb),intent(inout)               :: H

    integer         :: i_hop,i_pair

    type(neighbors) :: neigh            !all neighbor information for a given atom-type pair
    integer         :: at_pair(2)       !pair of atoms locally considered
    integer         :: orb(2)           !orbital offset in basic unit-cell
    integer         :: spin(2)          !spin offset in basic unit-cell
    integer         :: ind(2)           !index in basic unit-cell 
    integer         :: connect_bnd(2)   !indices keeping track of which pairs are used for the particular connection
    complex(8)      :: val(1)

    type(H_tb_coo)    :: Htmp         !local Hamiltonian to save 

    do i_hop=1,size(h_io%hop)
        Call neigh%get(h_io%hop(i_hop)%attype,[h_io%hop(i_hop)%dist],lat)

        connect_bnd=1
        do i_pair=1,neigh%Nshell(1)
            connect_bnd(2)=neigh%ishell(i_pair)
            at_pair=neigh%at_pair(:,i_pair)
            orb=h_io%norb_at_off(at_pair)+h_io%hop(i_hop)%orbital
            spin=h_io%hop(i_hop)%spin
            ind=spin+(orb-1)*h_io%nspin
            val=cmplx(h_io%hop(i_hop)%val ,0.0d0,8)

            !electron/electron-term
            Call Htmp%init_connect(neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val,reshape(ind,[2,1]),h_io)
            Call H%add(Htmp)
            Call Htmp%destroy()
            !hole/hole-term
            if(h_io%nsc==2) Call add_BdG(H,neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val,reshape(ind,[2,1]),h_io)
            connect_bnd(1)=connect_bnd(2)+1
        enddo
    enddo
end subroutine 


subroutine get_Hr_delta(lat,h_io,H)
    !TODO, make neigh%get more efficient by first sorting though h_io%del to get all required connections for an atom type pair at once
    use m_neighbor_type, only: neighbors
    type(lattice),intent(in)                :: lat
    type(parameters_TB_IO_H),intent(in)     :: h_io
    class(H_tb),intent(inout)               :: H

    integer         :: i_del,i_pair

    type(neighbors) :: neigh            !all neighbor information for a given atom-type pair
    integer         :: at_pair(2)       !pair of atoms locally considered
    integer         :: orb(2)           !orbital offset in basic unit-cell
    integer         :: ind(2,2)         !index in basic unit-cell 
    integer         :: connect_bnd(2)   !indices keeping track of which pairs are used for the particular connection
    complex(8)      :: val(2)
    integer         :: nBdG             !length to next BdG sector

    type(H_tb_coo)    :: Htmp         !local Hamiltonian to save 

    nBdG=h_io%norb*h_io%nspin*h_io%ncell
    do i_del=1,size(h_io%del)
        Call neigh%get(h_io%del(i_del)%attype,[h_io%del(i_del)%dist],lat)
        connect_bnd=1
        do i_pair=1,neigh%Nshell(1)
            connect_bnd(2)=neigh%ishell(i_pair)
            at_pair=neigh%at_pair(:,i_pair)
            orb=h_io%norb_at_off(at_pair)+h_io%del(i_del)%orbital
            ind(:,1)=(orb-1)*h_io%nspin+[1,2] !spin_up, spin_dn
            ind(:,2)=(orb-1)*h_io%nspin+[2,1] !spin_dn, spin_up
            !electron/hole-term
            val=h_io%del(i_del)%val
            Call Htmp%init_connect(neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val,ind,h_io,[0,nBdG])
            Call H%add(Htmp)
            Call Htmp%destroy()
            !hole/electron-term
            val=conjg(h_io%del(i_del)%val)
            Call Htmp%init_connect(neigh%pairs(:,connect_bnd(1):connect_bnd(2)),val,ind,h_io,[nBdG,0])
            Call H%add(Htmp)
            Call Htmp%destroy()

            connect_bnd(1)=connect_bnd(2)+1
        enddo
    enddo
end subroutine 


subroutine add_BdG(H,connect,Hval_in,Hval_ind,io)
    !adds the lower right quadrant of the BdG-equation from the upper quadrant
    use m_derived_types, only: lattice,op_abbrev_to_int
    class(H_tb),intent(inout)       :: H
    complex(8),intent(in)           :: Hval_in(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)              :: Hval_ind(2,size(Hval_in))
    integer,intent(in)              :: connect(:,:)  !(2,Nentries) index in (1:Ncell) basis of both connected sites 
    type(parameters_TB_IO_H),intent(in) :: io

    complex(8)     :: Hval(size(Hval_in))  !all entries between 2 cell sites of considered orderparameter

    integer     :: i

    type(H_tb_coo)    :: Htmp         !local Hamiltonian to save 

    hval=conjg(hval_in)
    forall(i=1:size(Hval), modulo(Hval_ind(1,i),2)==modulo(Hval_ind(2,i),2)) Hval(i)=-Hval(i)     !minus for spin-conserving terms
    Call Htmp%init_connect(connect,hval,hval_ind,io,[io%norb*io%nspin*io%ncell,io%norb*io%nspin*io%ncell])
    Call H%add(Htmp)
    Call Htmp%destroy()
end subroutine

subroutine get_Hr_Jsd(lat,h_io,H)
    type(lattice),intent(in)                :: lat
    type(parameters_TB_IO_H),intent(in)     :: h_io
    class(H_tb),intent(inout)               :: H

    integer,allocatable     :: at_arr(:)        !atom locally considered
    integer                 :: at
    integer                 :: orb_offset       !orbital offset in basic unit-cell
    integer                 :: ndim

    integer     ::  nnz
    complex(8),allocatable  ::  val(:)
    integer,allocatable     ::  row(:),col(:)


    type(H_tb_coo)    :: Htmp         !local Hamiltonian to save 
    integer ::  i
    integer ::  i_jsd,i_at,i_nnz,i_cell
    integer :: mag_offset

    ndim=h_io%norb*h_io%nspin
    do i_jsd=1,size(h_io%Jsd)
        at_arr=pack([(i,i=1,size(lat%cell%atomic))],lat%cell%atomic%type_id==h_io%jsd(i_jsd)%attype)
        nnz=size(at_arr)*h_io%ncell*4   !4=2**2 spin entries
        allocate(val(nnz),source=(0.0d0,0.0d0))
        allocate(row(nnz),col(nnz),source=0)
        i_nnz=0
        do i_at=1,size(at_arr)
            at=at_arr(i_at)
            mag_offset=3*count(lat%cell%atomic(1:at-1)%moment/=0.0d0)
            orb_offset=(h_io%norb_at_off(at)+h_io%jsd(i_jsd)%orbital-1)*h_io%nspin
            do i_cell=1,h_io%ncell
                row(i_nnz+1:i_nnz+4)=(i_cell-1)*ndim+orb_offset+[1,2,1,2]
                col(i_nnz+1:i_nnz+4)=(i_cell-1)*ndim+orb_offset+[1,1,2,2]
                val(i_nnz+1)=cmplx( lat%M%modes_v(mag_offset+3,i_cell), 0.0d0                             ,8)
                val(i_nnz+2)=cmplx( lat%M%modes_v(mag_offset+1,i_cell), lat%M%modes_v(mag_offset+2,i_cell),8)
                val(i_nnz+3)=cmplx( lat%M%modes_v(mag_offset+1,i_cell),-lat%M%modes_v(mag_offset+2,i_cell),8)
                val(i_nnz+4)=cmplx(-lat%M%modes_v(mag_offset+3,i_cell), 0.0d0                             ,8)
                i_nnz=i_nnz+4
            enddo
        enddo
        val=val*h_io%jsd(i_jsd)%val
        !electron/electron-term
        Call Htmp%init_coo(val,row,col,h_io)
        Call H%add(Htmp)
        Call Htmp%destroy()
        !hole/hole-term
        if(h_io%nsc==2) Call add_BdG_coo(H,val,row,col,h_io)
        deallocate(at_arr)
        deallocate(val,row,col)
    enddo
end subroutine 

subroutine add_BdG_coo(H,val,col,row,io)
    !adds the lower right quadrant of the BdG-equation from the upper quadrant
    !values get destroyed
    class(H_tb),intent(inout)       :: H
    complex(8),intent(inout)        :: val(:) 
    integer,intent(in)              :: row(size(val)),col(size(val))
    type(parameters_TB_IO_H),intent(in) :: io
    integer     :: i

    type(H_tb_coo)    :: Htmp         !local Hamiltonian to save 

    val=conjg(val)
    forall(i=1:size(val),modulo(row(i),2)==modulo(col(i),2)) val(i)=-val(i)   !minus for spin-conserving terms
    Call Htmp%init_coo(val,row,col,io,[io%norb*io%nspin*io%ncell,io%norb*io%nspin*io%ncell])
    Call H%add(Htmp)
    Call Htmp%destroy()
end subroutine



end module
