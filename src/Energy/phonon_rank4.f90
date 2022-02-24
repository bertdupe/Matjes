module m_phonon_rank4
use, intrinsic :: iso_fortran_env, only : output_unit
use m_input_H_types, only: io_H_PH4
use m_derived_types, only: lattice
use m_neighbor_pair, only: pair_dat_t, get_pair_dat_U
implicit none

character(len=*),parameter  :: ham_desc="4-spin interaction"
private
public :: read_PH4_input, get_PH4

contains

subroutine read_PH4_input(io_param,fname,io)
    use m_io_read_util, only: set_pos_entry, check_further_entry
    use, intrinsic :: iso_fortran_env, only : output_unit,error_unit
    integer,intent(in)              :: io_param
    character(len=*), intent(in)    :: fname
    type(io_H_Ph4),intent(out)      :: io
    character(len=*),parameter      :: var_name="U_Phonon_r4"

    logical     :: success
    integer     :: i,N,stat

    integer     :: at_type
    real(8)     :: val
    character(len=100)  :: str

    Call set_pos_entry(io_param,fname,var_name,success)
    read(io_param,*,iostat=stat) str !advance one line
    if(.not.success) return
    N=0
    do
        read(io_param,'(A100)',iostat=stat) str
        read(str,*,iostat=stat) at_type, val
        if(stat/=0) exit
        N=N+1
    enddo
    if(N<1)then
        write(error_unit,'(2/A)') 'Found entry "'//var_name//'" in file "'//fname//'", but no corresponding data was supplied'
#ifndef CPP_SCRIPT
        ERROR STOP "INPUT PROBABLY WRONG (disable with CPP_SCRIPT preprocessor flag)"
#endif
        return
    endif
    allocate(io%at_type(N),source=-1)
    allocate(io%val(N),source=0.0d0)
    do i=1,N+1
        backspace(io_param)
    enddo
    do  i=1,N
        read(io_param,'(A100)',iostat=stat) str
        read(str,*) io%at_type(i),io%val(i)
    enddo
    Call check_further_entry(io_param,fname,var_name)

    if(all(io%val==0))then
        write(error_unit,'(/2A/A/)') "WARNING, Found no nonzero entries for: ",var_name,' although the keyword is specified'
#ifndef CPP_SCRIPT
        ERROR STOP "INPUT PROBABLY WRONG (disable with CPP_SCRIPT preprocessor flag)"
#endif
        return
    endif

    io%is_set=.true.
end subroutine

subroutine get_PH4(Ham,io,lat)
    !get rank 4 phonon interaction in t_H Hamiltonian format
    use m_H_public, only: t_H, get_Htype
    use m_mode_public
    use m_setup_4spin
    use m_coo_mat

    class(t_H),intent(inout)    :: Ham
    type(io_H_Ph4),intent(in)   :: io
    type(lattice),intent(in)    :: lat

    type(pair_dat_t),allocatable    :: pair_dat(:)

    integer         ::  N4phonon, i_4phonon
    integer         :: i
    integer         :: N_con    !number of different connections per unit-cell combined in modes
    integer         :: Nplaq_base, iplaq_base

    !infomation for about 4-spin
    integer,allocatable     :: con_unfold(:,:)     !all sites of the base-points of the connections
    integer,allocatable     :: N_entry(:)          !number of unfolded plaques per base-plaque
    integer,allocatable     :: connect_neigh(:,:)  !index of neighbor step

    !variables to construct Hamiltonian
    integer                 :: dim_mode(2)
    integer                 :: neigh(6) !local neighbor indices
    real(8)                 :: K
    integer,allocatable     :: connect(:,:)
    class(t_H),allocatable  :: Ham_tmp    !temporary Hamiltonian type used to add up Ham
    real(8)                 :: val_tmp(9)
    integer                 :: ind_tmp(2,9)
        !helper to get all neighbor indices of the (xx,yx,zx,yx,yy,yz,zx,zy,zz)-combinations
    integer,parameter   :: ii_jj_comb(2,9)=reshape([1,1, 2,1, 3,1,  1,2, 2,2, 3,2,  1,3, 2,3, 3,3],[2,9])

    !parameter to construct the modes
    type(coo_mat)       :: mat(2)       !mode construction matrices of left/right side of Hamiltonian ( first left, reused for right)
    integer             :: dim_mat(2),nnz
    integer,allocatable :: row(:,:),col(:,:)
    real(8),allocatable :: val(:,:)
    integer             :: bnd(2)
    integer,allocatable :: pairs(:,:)
    integer,allocatable :: at_Ph(:)
    integer,allocatable :: pos_offset(:,:)

    if(io%is_set)then
        write(output_unit,'(/2A)') "Start setting Hamiltonian: ", ham_desc
        Call get_Htype(Ham_tmp)
        N4phonon=size(io%val)
        Call get_pair_dat_U(lat,io%at_type,[1,2],pair_dat)    !get all pairs of same atom type for nearest and next-nearest neighbor
        N_con=sum(pair_dat%Npair)
        dim_mode=[3*N_con,3*N_con]
        do i_4phonon=1,N4phonon
            K=io%val(i_4phonon)
            Call get_4spin_dat(lat,pair_dat(i_4phonon),con_unfold,N_entry,connect_neigh)
            connect_neigh=connect_neigh+sum(pair_dat(:i_4phonon-1)%Npair) !adjust to mode considering other 4-spin interactions
            Nplaq_base=size(N_entry)
            do iplaq_base=1,Nplaq_base
                if(N_entry(iplaq_base)==0) cycle
                bnd=[sum(N_entry(:iplaq_base-1))+1,sum(N_entry(:iplaq_base))]
                allocate(connect(2,N_entry(iplaq_base)),source=0)
                neigh=(connect_neigh(:,iplaq_base)-1)*3 !+1->x +2->y +3->z
                 !12 34 connection
                 connect(:,:)=con_unfold(1:2,bnd(1):bnd(2))
                 val_tmp=-K
                 ind_tmp=spread(neigh(1:2),2,9)+ii_jj_comb
                 Call Ham_tmp%init_mult_connect_2(connect,val_tmp,ind_tmp,"UU","UU",lat,4,dim_mode_in=dim_mode)
                 Call Ham%add(Ham_tmp)
                 Call Ham_tmp%destroy()

                 !14 23 connection
                 connect(:,:)=con_unfold(3:4,bnd(1):bnd(2))
                 val_tmp=-K
                 ind_tmp=spread(neigh(3:4),2,9)+ii_jj_comb
                 Call Ham_tmp%init_mult_connect_2(connect,val_tmp,ind_tmp,"UU","UU",lat,4,dim_mode_in=dim_mode)
                 Call Ham%add(Ham_tmp)
                 Call Ham_tmp%destroy()

                 !13 24 connection
                 connect(:,:)=con_unfold(5:6,bnd(1):bnd(2))
                 val_tmp=K
                 ind_tmp=spread(neigh(5:6),2,9)+ii_jj_comb
                 Call Ham_tmp%init_mult_connect_2(connect,val_tmp,ind_tmp,"UU","UU",lat,4,dim_mode_in=dim_mode)
                 Call Ham%add(Ham_tmp)
                 Call Ham_tmp%destroy()
                deallocate(connect)
            enddo
            deallocate(con_unfold,N_entry,connect_neigh)
        enddo

        !SET STATES
        !combine the pair data from all considered 4-spin interactions
        allocate(pairs(5,N_con),source=0)
        do i_4phonon=1,N4phonon
            bnd=[sum(pair_dat(:i_4phonon-1)%Npair)+1,sum(pair_dat(:i_4phonon)%Npair)]
            pairs(:,bnd(1):bnd(2))=pair_dat(i_4phonon)%pairs
            pairs(1,bnd(1):bnd(2))=pair_dat(i_4phonon)%ind_ph(pairs(1,bnd(1):bnd(2)))
            pairs(2,bnd(1):bnd(2))=pair_dat(i_4phonon)%ind_ph(pairs(2,bnd(1):bnd(2)))
        enddo
        !set sizes, allocate coo-matrix variables, and set trivial values
        dim_mat=[3*N_con*lat%Ncell,lat%U%dim_mode*lat%Ncell]
        nnz=dim_mat(1)
        allocate(row(nnz,2),col(nnz,2),source=0)
        allocate(val(nnz,2),source=1.0d0)
        row(:,1)=[(i,i=1,nnz)]
        row(:,2)=row(:,1)
        allocate(at_Ph(N_con),source=0)
        allocate(pos_offset(3,N_con),source=0)

        at_Ph=pairs(1,:)
        pos_offset=spread([0,0,0],2,N_con)
        Call set_mode_pair(lat%Ncell,N_con,lat,at_Ph,pos_offset,col(:,1))

        at_Ph=pairs(2,:)
        pos_offset=pairs(3:5,:)
        Call set_mode_pair(lat%Ncell,N_con,lat,at_Ph,pos_offset,col(:,2))

        !both M-states have same construction procedure
        Call mat(1)%init(dim_mat,nnz,row(:,1),col(:,1),val(:,1))
        Call mat(2)%init(dim_mat,nnz,row(:,2),col(:,2),val(:,2))
        Call mode_set_rankN_sparse(Ham%mode_l,"UU",lat,mat,1)

        !recreate mat since it is destroyed in initialization
        Call mat(1)%init(dim_mat,nnz,row(:,1),col(:,1),val(:,1))
        Call mat(2)%init(dim_mat,nnz,row(:,2),col(:,2),val(:,2))
        Call mode_set_rankN_sparse(Ham%mode_r,"UU",lat,mat,1)
        !END SET STATES
        Ham%desc=ham_desc
    endif
end subroutine


subroutine set_mode_pair(Ncell,Ncon,lat,at_ph,pos_offset,col)
    !sets the column indices for a mode array by unfolding the supercell of the atom with the at_mag magnetic index
    !and assuming an position offset of pos_offset
    integer,intent(in)              :: Ncell
    integer,intent(in)              :: Ncon
    type(lattice),intent(in)        :: lat
    integer,intent(in)              :: at_ph(Ncon)
    integer,intent(in)              :: pos_offset(3,Ncon)
    integer,intent(out)             :: col(3,Ncon,Ncell)

    integer     :: ind3_lat(3,Ncell)
    integer     :: ind3_unfold(3,Ncon,Ncell)
    integer     :: ind1_unfold(Ncon,Ncell)
    integer     :: dim_lat(3)
    integer     :: i1,i2,i3, ind3(3),ii

    integer     :: dim_lat_spread(3,Ncell)

    dim_lat=lat%dim_lat
    dim_lat_spread=spread(dim_lat,2,Ncell)

    ii=0    !unfold full lattice
    do i3=1,dim_lat(3)
        ind3(3)=i3
        do i2=1,dim_lat(2)
            ind3(2)=i2
            do i1=1,dim_lat(1)
                ind3(1)=i1
                ii=ii+1
                ind3_lat(:,ii)=ind3
            enddo
        enddo
    enddo

    ind3_unfold=spread(ind3_lat,2,Ncon)
    ind3_unfold=ind3_unfold+spread(pos_offset,3,Ncell)
    Call lat%fold_3_arr(Ncell*Ncon,ind3_unfold)
    Call lat%index_3_1_arr(Ncell*Ncon,ind3_unfold,ind1_unfold)  !get 1d lattice sites
    ind1_unfold=(ind1_unfold-1)*lat%nph+spread(at_Ph,2,Ncell) !get 1d magnetic atom entry
    ind1_unfold=ind1_unfold-1
    ind1_unfold=ind1_unfold*3
    do i2=1,Ncell
        do i1=1,Ncon
            col(1,i1,i2)=ind1_unfold(i1,i2)+1
            col(2,i1,i2)=ind1_unfold(i1,i2)+2
            col(3,i1,i2)=ind1_unfold(i1,i2)+3
        enddo
    enddo
end subroutine

end module m_phonon_rank4
