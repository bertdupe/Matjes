module m_init_H
use m_derived_types, only: lattice
use m_H_tb_public
use m_types_tb_h_inp 
use m_delta_onsite
use,intrinsic :: iso_fortran_env, only : output_unit,error_unit
 
use m_ham_init_type, only: parameters_ham_init
use m_neighbor_type, only: neighbors

private
public :: get_H, get_H_TB, get_H_sc
contains

subroutine get_H(lat,h_io,H_out,diffR)
    !routine that returns arrays of all Hamiltonian-components (optionally with the real-space differences) either in the normal- or in the BdG-space
    !if a superconducting delta is used, all Hamiltonians are transformed to the BdG space
    type(lattice),intent(in)                            :: lat
    type(parameters_TB_IO_H),intent(in)                 :: h_io
    type(H_tb_coo),allocatable,intent(inout)            :: H_out(:)
    real(8),allocatable,intent(out),optional            :: diffR(:,:)

    type(H_tb_coo),allocatable  :: H(:)
    real(8),allocatable         :: diffR_loc(:,:)
    integer                     :: i
    !get all normal-conducting Hamiltonian components
    Call get_H_TB(lat,h_io,H_out,diffR=diffR)   

    !get superconducting Hamiltonian entries
    if(allocated(h_io%del_io))then 
        allocate(H(1))
        Call get_H_sc(lat,h_io,H,h_io%del,diffR=diffR_loc)
        Call H_append(H_out,H)
        if(present(diffR)) Call append_arr(diffR,diffR_loc)
    endif

    !translate all Hamiltonians to the BdG space
    if(any(H_out%nsc>1))then    
        do i=1,size(H_out)
            Call H_out(i)%to_BdG()
        enddo
    endif
end subroutine 


subroutine get_H_TB(lat,h_io,H_out,sc,diffR)
    !sets the Hamiltonian terms without superconductivity,
    !i.e. only Block-diagonal terms in BdG-space
    type(lattice),intent(in)                            :: lat
    type(parameters_TB_IO_H),intent(in)                 :: h_io         !Hamiltonian input parameters
    type(H_tb_coo),allocatable,intent(inout)            :: H_out(:)     !Hamiltonian array gathering all found interactions
    logical,intent(in),optional                         :: sc           !transform terms to BdG-space
    real(8),allocatable,intent(out),optional            :: diffR(:,:)   !real-space difference lattice vectors H_out entry

    type(H_tb_coo),allocatable  :: H(:)             !temporary Hamiltonian array 
    real(8),allocatable         :: diffR_loc(:,:)   !corresponding temporary real-space difference arrays
    integer                     :: i

    if(h_io%wann_io%nrpts/=0)then
        Call get_wann_arr(lat,h_io,H,diffR=diffr_loc)
        Call H_append(H_out,H)
        if(present(diffR)) Call append_arr(diffR,diffR_loc)
        if(allocated(diffR_loc)) deallocate(diffR_loc)
    endif

    !Hopping term  (input: TB_hopping)
    if(allocated(h_io%hop_io))then
        Call get_Hop_arr(lat,h_io,H,diffR=diffR_loc)
        Call H_append(H_out,H)
        if(present(diffR)) Call append_arr(diffR,diffR_loc)
        if(allocated(diffR_loc)) deallocate(diffR_loc)
    endif

    !Jsd-coupling term  (input: TB_Jsd)
    if(allocated(h_io%jsd))then
        allocate(H(1))
        Call get_H_Jsd(lat,h_io,H(1))
        Call H_append(H_out,H)
        if(present(diffR))then
            allocate(diffR_loc(3,1),source=0.0d0) !single onsite term, so difference is 0
            Call append_arr(diffR,diffR_loc)
        endif
    endif

    !defect term,potential at certain site  (input: TB_defect)
    if(allocated(h_io%defect))then
        allocate(H(1))
        Call get_H_defect(lat,h_io,H(1))
        Call H_append(H_out,H)
        if(present(diffR))then
            allocate(diffR_loc(3,1),source=0.0d0) !single onsite term, so difference is 0
            Call append_arr(diffR,diffR_loc)
        endif
    endif

    !diagonal term adjusting the fermi energy (input: TB_Efermi)
    if(h_io%Efermi/=0.0d0)then
        allocate(H(1))
        Call get_H_fermi(lat,h_io,H(1))
        Call H_append(H_out,H)
        if(present(diffR))then
            allocate(diffR_loc(3,1),source=0.0d0) !single onsite term, so difference is 0
            Call append_arr(diffR,diffR_loc)
        endif
    endif

    !transform all entries to BdG-space already at this level
    if(present(sc))then
        if(sc)then
            do i=1,size(H_out)
                Call H_out(i)%to_BdG()
            enddo
        endif
    endif
end subroutine 

subroutine get_H_sc(lat,h_io,H_out,del,diffR)
    type(lattice),intent(in)                            :: lat
    type(parameters_TB_IO_H),intent(in)                 :: h_io
    type(H_tb_coo),allocatable,intent(inout)            :: H_out(:)
    type(Hdelta),intent(in)                             :: del
    real(8),allocatable,intent(out),optional            :: diffR(:,:)

    if(.not.allocated(h_io%del_io)) STOP "h_io%del_io should be allocated when calculating delta"
    if(allocated(H_out))then
        if(size(H_out)/=1) STOP "get_H_sc H_out has wrong size"
    else
        allocate(H_out(1)) !increase if more than onsite-delta
    endif
    Call get_H_delta_onsite(lat,h_io,H_out(1),del)
    if(present(diffR)) allocate(diffR(3,1),source=0.0d0) !single onsite term, so difference is 0
end subroutine 


subroutine get_H_delta_onsite(lat,h_io,H,del)
    type(lattice),intent(in)                :: lat
    type(parameters_TB_IO_H),intent(in)     :: h_io
    class(H_tb),intent(inout)               :: H
    type(Hdelta),intent(in)                 :: del

    integer         :: i_del,i_cell, i_nnz

    type(parameters_ham_init)   ::  hinit   !type containing variables defining shape of Hamiltonian
    integer         :: orb(2)           !orbital offset in basic unit-cell
    integer         :: nBdG             !length to next BdG sector
    integer         :: ndim             !length to next unit-cell
    integer         :: add_row(4),add_col(4)

    type(H_tb_coo)  :: Htmp             !local Hamiltonian to save 
    integer,allocatable ::  row(:),col(:)
    complex(8),allocatable  :: val(:)

    ndim=h_io%norb*h_io%nspin
    Call hinit%init(h_io)
    nBdG=hinit%norb*hinit%nspin*hinit%ncell
    allocate(row(lat%ncell*4),col(lat%ncell*4),val(lat%ncell*4))
    !order of entries 
    add_row=[1     , 2     , nBdG+1, nBdG+2]
    add_col=[nBdG+2, nBdG+1, 2     , 1     ]
    do i_del=1,size(del%orb,2)
        orb=(del%orb(:,i_del)-1)*hinit%nspin
        i_nnz=0
        do i_cell=1,lat%ncell
            row(i_nnz+1:i_nnz+4)=(i_cell-1)*ndim+orb(1)+add_row
            col(i_nnz+1:i_nnz+4)=(i_cell-1)*ndim+orb(2)+add_col
            val(i_nnz+1:i_nnz+4)=del%delta(i_cell,i_del)
            val(i_nnz+3:i_nnz+4)=conjg(val(i_nnz+3:i_nnz+4))
            i_nnz=i_nnz+4
        enddo
        Call Htmp%init_coo(val,row,col,hinit)
        Call H%add(Htmp)
        Call Htmp%destroy()
    enddo
end subroutine 

subroutine get_H_fermi(lat,h_io,H)
    type(lattice),intent(in)                    :: lat
    type(parameters_TB_IO_H),intent(in)         :: h_io
    type(H_tb_coo),intent(inout)                :: H

    type(parameters_ham_init)   ::  hinit   !type containing variables defining shape of Hamiltonian
    integer                 ::  nnz,i_nnz
    complex(8),allocatable  ::  val(:)
    integer,allocatable     ::  row(:),col(:)

    Call hinit%init(h_io)
    Hinit%nsc=1 !set hoppings without BdG
    nnz=h_io%norb*h_io%nspin*h_io%ncell
    val=[(-h_io%Efermi,i_nnz=1,nnz)]
    row=[(i_nnz,i_nnz=1,nnz)]
    col=[(i_nnz,i_nnz=1,nnz)]
    Call H%init_coo(val,row,col,hinit)
    deallocate(val,row,col)
end subroutine 


subroutine get_H_defect(lat,h_io,H)
    type(lattice),intent(in)                    :: lat
    type(parameters_TB_IO_H),intent(in)         :: h_io
    type(H_tb_coo),intent(inout)                :: H

    integer                 :: at
    integer                 :: orb_offset       !orbital offset in basic unit-cell
    integer                 :: mag_offset
    integer                 :: ndim

    integer                 :: nnz
    complex(8),allocatable  :: val(:)
    integer,allocatable     :: row(:),col(:)
    real(8)                 :: mag(3)
    logical                 :: has_spin

    type(parameters_ham_init)   ::  hinit   !type containing variables defining shape of Hamiltonian
    type(H_tb_coo)    :: Htmp         !local Hamiltonian to save 
    integer ::  i_defect,i_nnz,i_cell

    has_spin=h_io%nspin>1
    ndim=h_io%norb*h_io%nspin
    Call hinit%init(h_io)
    Hinit%nsc=1 !set hoppings without BdG
    nnz=size(h_io%defect)
    if(has_spin) nnz=size(h_io%defect)*4

    allocate(val(nnz),source=(0.0d0,0.0d0))
    allocate(row(nnz),col(nnz),source=0)
    i_nnz=0
    if(has_spin)then
        do i_defect=1,size(h_io%defect)
            i_cell=lat%index_m_1(h_io%defect(i_defect)%site)
            at=h_io%defect(i_defect)%atom
            mag_offset=3*count(lat%cell%atomic(1:at-1)%moment/=0.0d0)
            orb_offset=(h_io%norb_at_off(at)+h_io%defect(i_defect)%orbital-1)*hinit%nspin
            row(i_nnz+1:i_nnz+4)=(i_cell-1)*ndim+orb_offset+[1,2,1,2]
            col(i_nnz+1:i_nnz+4)=(i_cell-1)*ndim+orb_offset+[1,1,2,2]
            mag=lat%M%modes_v(mag_offset+1:mag_offset+3,i_cell)

            val(i_nnz+1)=cmplx( mag(3)*h_io%defect(i_defect)%mag, 0.0d0                             ,8)
            val(i_nnz+2)=cmplx( mag(1)*h_io%defect(i_defect)%mag, mag(2)*h_io%defect(i_defect)%mag  ,8)
            val(i_nnz+3)=cmplx( mag(1)*h_io%defect(i_defect)%mag,-mag(2)*h_io%defect(i_defect)%mag  ,8)
            val(i_nnz+4)=cmplx(-mag(3)*h_io%defect(i_defect)%mag, 0.0d0                             ,8)
            val(i_nnz+1)=val(i_nnz+1)+cmplx(h_io%defect(i_defect)%nonmag,0.0d0,8)
            val(i_nnz+4)=val(i_nnz+4)+cmplx(h_io%defect(i_defect)%nonmag,0.0d0,8)
            i_nnz=i_nnz+4
        enddo
    else
        do i_defect=1,size(h_io%defect)
            i_cell=lat%index_m_1(h_io%defect(i_defect)%site)
            at=h_io%defect(i_defect)%atom
            orb_offset=(h_io%norb_at_off(at)+h_io%defect(i_defect)%orbital-1)*hinit%nspin
            row(i_nnz+1)=(i_cell-1)*ndim+orb_offset+1
            col(i_nnz+1)=(i_cell-1)*ndim+orb_offset+1
            val(i_nnz+1)=val(i_nnz+1)+cmplx(h_io%defect(i_defect)%nonmag,0.0d0,8)
            i_nnz=i_nnz+1
        enddo
    endif
    Call Htmp%init_coo(val,row,col,hinit)
    Call H%add(Htmp)
    Call Htmp%destroy()
    deallocate(val,row,col)
end subroutine 



subroutine get_H_Jsd(lat,h_io,H)
    type(lattice),intent(in)                    :: lat
    type(parameters_TB_IO_H),intent(in)         :: h_io
    type(H_tb_coo),intent(inout)                :: H

    integer,allocatable     :: at_arr(:)        !atom locally considered
    integer                 :: at
    integer                 :: orb_offset       !orbital offset in basic unit-cell
    integer                 :: ndim

    integer                 ::  nnz
    complex(8),allocatable  ::  val(:)
    integer,allocatable     ::  row(:),col(:)

    type(parameters_ham_init)   ::  hinit   !type containing variables defining shape of Hamiltonian
    type(H_tb_coo)    :: Htmp         !local Hamiltonian to save 
    integer ::  i
    integer ::  i_jsd,i_at,i_nnz,i_cell
    integer :: mag_offset

    ndim=h_io%norb*h_io%nspin
    Call hinit%init(h_io)
    Hinit%nsc=1 !set hoppings without BdG
    do i_jsd=1,size(h_io%Jsd)
        at_arr=pack([(i,i=1,size(lat%cell%atomic))],lat%cell%atomic%type_id==h_io%jsd(i_jsd)%attype)
        nnz=size(at_arr)*h_io%ncell*4   !4=2**2 spin entries
        allocate(val(nnz),source=(0.0d0,0.0d0))
        allocate(row(nnz),col(nnz),source=0)
        i_nnz=0
        do i_at=1,size(at_arr)
            at=at_arr(i_at)
            mag_offset=3*count(lat%cell%atomic(1:at-1)%moment/=0.0d0)
            orb_offset=(h_io%norb_at_off(at)+h_io%jsd(i_jsd)%orbital-1)*hinit%nspin
            do i_cell=1,hinit%ncell
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
        Call Htmp%init_coo(val,row,col,hinit)
        Call H%add(Htmp)
        Call Htmp%destroy()
        deallocate(at_arr)
        deallocate(val,row,col)
    enddo
end subroutine 


subroutine get_wann_arr(lat,h_io,Hhop,diffR)
    !gets all hopping Hamiltonians according to the input from h_io%hop
    !saves each entries with differing atom connections and connection direction separately in the Hhop-array
    !if neigh_out is specified, returns the calculated neighbor entires
    !if diffR is specified, returns the difference vector from the first atom the second for each Hamiltonian (used to calculate H in k-space)
    use m_ham_arrange
    use m_setH_util, only: get_coo
    type(lattice),intent(in)                            :: lat
    type(parameters_TB_IO_H),intent(in)                 :: h_io
    type(H_tb_coo),allocatable,intent(inout)            :: Hhop(:)
    real(8),allocatable,intent(out),optional            :: diffR(:,:)

    type(parameters_ham_init)   :: hinit        !type containing variables defining shape of Hamiltonian
    integer                     :: N_R, irpts   !number and loop variable of cell-differences
    integer,allocatable         :: ipair(:,:)   !saves which supercells are connected by the interaction (indices)
    logical,allocatable         :: H_is_set(:)  !save which Hamiltonian entry as been set
    !parameters for Hamiltonian at given cell-difference in coo sparse matrix format
    complex(8),allocatable      :: val(:)       !coo Hamiltonian values
    integer,allocatable         :: ind(:,:)     !coo Hamiltonian indices

    !prepare parameters and allocate arrays
    Call hinit%init(h_io)
    Hinit%nsc=1     !set hoppings without BdG
    N_r=h_io%wann_io%nrpts  !number of real-space differences (number of different Hamiltonian entreis)
    allocate(Hhop(N_r))
    allocate(H_is_set(N_r),source=.false.)
    if(present(diffR)) allocate(diffR(3,N_r),source=0.0d0)

    !loop over all different Hamiltonian differences
    do irpts=1,N_r
        Call get_coo(h_io%wann_io%H(:,:,irpts),val,ind)                 !transform dense wannier-input to sparse coo-matrix
        Call unfold_neighbors(h_io%wann_io%irvec(:,irpts), lat, ipair)  !get ipair which encodes distribution in the supercell on current interaction
        if(allocated(ipair))then
            Call Hhop(irpts)%init_connect(ipair,val,ind,hinit)          !initialize Hamiltonian
            H_is_set(irpts)=.true.
            deallocate(ipair)
        endif
        if(present(diffR)) diffR(:,irpts)=matmul(real(h_io%wann_io%irvec(:,irpts),8),lat%areal)
    enddo

    !remove unused entries (Hhop,diffR)
    if(.not.all(H_is_set))  Call H_pick_entries(H_is_set,Hhop,diffR)
end subroutine

subroutine H_pick_entries(H_use, Harr, diffR)
    !subroutine, which takes an Hamiltonian array and optionally a corresponding difference vector array and removes all entries
    !for whose index the H_use array is false
    logical,intent(in)                          :: H_use(:)     !true for entries that shall remain in the arrays
    type(H_tb_coo),allocatable,intent(inout)    :: Harr(:)      !Hamiltonian array of which the entries encoded by H_use are retained (resizing)
    real(8),allocatable,intent(inout),optional  :: diffR(:,:)   !difference vector array of which the entries encoded by H_use are retained (resizing)
    type(H_tb_coo),allocatable      :: Harr_tmp(:)      !temporary Hamiltonian array storage
    real(8),allocatable             :: diffR_tmp(:,:)   !temporary difference vector array storage
    integer     ::  N_keep,i_keep, i                    !loop integers

    if(size(Harr)/=size(H_use))then
        write(error_unit,'(/)')
        write(error_unit,'(A,I8)') "H_use array size:", size(H_use)
        write(error_unit,'(A,I8)') "Harr  array size:", size(Harr)
        ERROR STOP "CALLED H_pick_entries with incompatible Hamiltonian lengths"
    endif
    !move entries of Hamiltonian array
    N_keep=count(H_use)
    Call move_alloc(Harr,Harr_tmp)
    allocate(Harr(N_keep))
    i_keep=0
    do i=1,size(H_use)
        if(H_use(i))then
            i_keep=i_keep+1
            Call Harr_tmp(i)%mv(Harr(i_keep))
        endif
    enddo
    do i=1,size(Harr_tmp)
        Call Harr_tmp(i)%destroy()
    enddo
    deallocate(Harr_tmp)

    if(present(diffR))then
        if(any(shape(diffR,2)/=[3,size(H_use)]))then
            write(error_unit,'(A,2I6)') "diffR shape :", shape(diffR)
            write(error_unit,'(A,2I6)') "wanted shape:", 3,size(H_use)
            ERROR STOP "CALLED H_pick_entries with incompatible Hamiltonian lengths"
        endif
        !move entries of difference vector array
        Call move_alloc(diffR,diffR_tmp)
        allocate(diffR(3,N_keep))
        i_keep=0
        do i=1,size(H_use)
            if(H_use(i))then
                i_keep=i_keep+1
                diffR(:,i_keep)=diffR_tmp(:,i)
            endif
        enddo
        deallocate(diffR_tmp)
    endif
end subroutine

subroutine unfold_neighbors(diff,lat,ipair)
    integer,intent(in)                  ::  diff(3)
    type(lattice),intent(in)            ::  lat
    integer,intent(inout),allocatable   ::  ipair(:,:)

    integer         :: imax(2,3)
    logical,save    :: said=.false.
    integer         :: N_entry
    integer         :: i3_1(3),i3_2(3)  !lattice sites in ([1,dimlat(1)],[1,dimlat(2)],[1,dimlat(3)])-basis
    integer         :: i1,i2,i3
    integer         :: ii

    !get boundaries of connection in each direction obeying periodicity
    !for each atom1 with a site within imax the connection to atom2 exists
    imax=reshape([1,lat%dim_lat(1),1,lat%dim_lat(2),1,lat%dim_lat(3)],[2,3])
    if(.not.lat%periodic(1).and.diff(1)<0) imax(1,1)=1             -diff(1)
    if(.not.lat%periodic(1).and.diff(1)>0) imax(2,1)=lat%dim_lat(1)-diff(1)
    if(.not.lat%periodic(2).and.diff(2)<0) imax(1,2)=1             -diff(2)
    if(.not.lat%periodic(2).and.diff(2)>0) imax(2,2)=lat%dim_lat(2)-diff(2)
    if(.not.lat%periodic(3).and.diff(3)<0) imax(1,3)=1             -diff(3)
    if(.not.lat%periodic(3).and.diff(3)>0) imax(2,3)=lat%dim_lat(3)-diff(3)

    if(any(imax(1,:)>imax(2,:)))then
        if(.not.said) write(error_unit,'(/A/)') "WARNING, some Wannier Hamiltonian entries are disregarded because of missing periodicity"
        said=.true.
        return
    endif
    
    N_entry=(imax(2,1)-imax(1,1)+1)*(imax(2,2)-imax(1,2)+1)*(imax(2,3)-imax(1,3)+1)
    allocate(ipair(2,N_entry),source=0)
    ii=0
    !fill all pairs within the imax boundaries
    do i3=imax(1,3),imax(2,3)
        i3_1(3)=i3
        do i2=imax(1,2),imax(2,2)
            i3_1(2)=i2
            do i1=imax(1,1),imax(2,1)
                i3_1(1)=i1
                i3_2=modulo(i3_1+diff-1,lat%dim_lat)+1
                ii=ii+1
                ipair(1,ii)=lat%index_m_1(i3_1)
                ipair(2,ii)=lat%index_m_1(i3_2)
            enddo
        enddo
    enddo
end subroutine

subroutine get_Hop_arr(lat,h_io,Hhop,neigh_out,diffR)
    !gets all hopping Hamiltonians according to the input from h_io%hop
    !saves each entries with differing atom connections and connection direction separately in the Hhop-array
    !if neigh_out is specified, returns the calculated neighbor entires (for so far unimplemented later reuse) (TODO,->intent(inout), check if correct and then reuse)
    !if diffR is specified, returns the difference vector from the first atom the second for each Hamiltonian (used to calculate H in k-space)
    use m_ham_arrange
    use m_setH_util, only: get_coo
    type(lattice),intent(in)                            :: lat
    type(parameters_TB_IO_H),intent(in)                 :: h_io
    type(H_tb_coo),allocatable,intent(inout)            :: Hhop(:)
    type(neighbors),intent(out),allocatable,optional    :: neigh_out(:)
    real(8),allocatable,intent(out),optional            :: diffR(:,:)

    type(neighbors),allocatable     :: neigh(:)
    type(parameters_ham_init)       :: hinit   !type containing variables defining shape of Hamiltonian
    integer :: Natt_pair !number of unique atom-type pairs
    integer :: N_pair   !number of unique atom pairs
    integer :: attpair(2) !atom (type) pair ids
    integer,allocatable     :: dist(:)
    integer :: i_att, i_dist,i_hop, i_shell, i_neigh_shell

    Natt_pair=size(h_io%hop%at)
    allocate(neigh(Natt_pair))
    N_pair=0
    do i_att=1,Natt_pair
        dist=h_io%hop%at(i_att)%dist%dist
        attpair=h_io%hop%at(i_att)%atpair
        Call neigh(i_att)%get(attpair,dist,lat)
        Call neigh(i_att)%prt()
        N_pair=N_pair+sum(neigh(i_att)%Nshell)
        deallocate(dist)
    enddo

    !get Hamiltonians 
    Call hinit%init(h_io)
    Hinit%nsc=1 !set hoppings without BdG
    allocate(Hhop(N_pair))
    if(present(diffR)) allocate(diffR(3,N_pair),source=0.0d0)
    i_hop=0
    do i_att=1,Natt_pair
        do i_dist=1,size(h_io%hop%at(i_att)%dist)
            do i_shell=1,neigh(i_att)%Nshell(i_dist)
                i_hop=i_hop+1
                i_neigh_shell=i_shell+sum(neigh(i_att)%Nshell(1:i_dist-1))
                Call set_Hop(i_neigh_shell,neigh(i_att),h_io%hop%at(i_att)%dist(i_dist)%Hloc,h_io,hinit,Hhop(i_hop))
                if(present(diffR)) diffR(:,i_hop)=neigh(i_att)%diff_vec(:,i_neigh_shell)
            enddo
        enddo
    enddo

    if(present(neigh_out)) Call move_alloc(neigh,neigh_out)
end subroutine

subroutine set_Hop(i_shell,neigh,Hloc,h_io,hinit,H)
    !routine to get supplied sparse Hamiltonian values, get correct indices for the neighbor shell, and initialize Hamiltonian
    use m_setH_util, only: get_coo
    integer                     :: i_shell  !shell index of neigh considered here
    type(neighbors),intent(in)  :: neigh    !neighbor type knowing connections and atom pairs (here consider i_shell)
    complex(8),intent(in)       :: Hloc(:,:)
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(parameters_ham_init),intent(in)    :: hinit   !type containing variables defining shape of Hamiltonian
    type(H_tb_coo)                          :: H

    complex(8),allocatable  :: val(:)
    integer,allocatable     :: ind(:,:)
    integer,allocatable     :: ind_loc(:,:)
    integer                 :: i1,i2    !pair boundary indices

    i1=1
    if(i_shell>1) i1=neigh%ishell(i_shell-1)+1
    i2=neigh%ishell(i_shell)

    Call get_coo(Hloc,val,ind)
    allocate(ind_loc,mold=ind)
    ind_loc=ind+spread(h_io%norb_at_off(neigh%at_pair(:,i_shell))*hinit%nspin,2,size(ind,2)) !add offset for actually considered atoms (and not atom-types)
    Call H%init_connect(neigh%pairs(:,i1:i2),val,ind_loc,hinit)
end subroutine

subroutine append_arr(arr,arr_append)
    !takes 2 rank-2 real-arrays with the first size in the 1. dimension and appends the second array to the first (with some allocation and moving)
    real(8),intent(inout),allocatable       :: arr(:,:), arr_append(:,:)
    real(8),allocatable :: tmp(:,:)

    if(.not.allocated(arr))then
        Call move_alloc(arr_append,arr)
        return
    endif
    if(size(arr,1)/=size(arr_append,1)) STOP "CANNOT append real arrays whose first ranks differ"
    allocate(tmp(size(arr,1),size(arr,2)+size(arr_append,2)))
    tmp(:,1:size(arr,2))=arr
    tmp(:,size(arr,2)+1:size(tmp,2))=arr_append
    deallocate(arr,arr_append)
    Call move_alloc(tmp,arr)
end subroutine

end module
