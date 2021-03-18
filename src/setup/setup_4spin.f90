module m_setup_4spin
use m_derived_types, only: lattice
use m_neighbor_pair, only: pair_dat_t
implicit none
private
public :: get_4spin_dat
contains

subroutine get_4spin_dat(lat,pair_dat,con_unfold,N_entry,connect_neigh)
    !subroutine that gives the information about the connections of all plaques in the supercell 
    use m_neighbor_type
    type(lattice),intent(in)            :: lat
    type(pair_dat_t),intent(in)         :: pair_dat
    integer,intent(inout),allocatable   :: con_unfold(:,:)     !all sites of the base-points of the connections
    integer,intent(inout),allocatable   :: N_entry(:)          !number of unfolded plaques per base-plaque
    integer,intent(inout),allocatable   :: connect_neigh(:,:)  !index of neighbor step

    integer,allocatable :: plaq_base(:,:,:) !plaque information without superlattice consideration
    integer,allocatable :: used_con(:,:)     !(2,Nplaq), which connections are used 
    integer,allocatable :: connect_base(:,:,:)  !site (at,pos) of base point of connection
    integer :: Nplaq_base

    !delete possible previous data
    if(allocated(con_unfold))       deallocate(con_unfold)
    if(allocated(N_entry))          deallocate(N_entry)
    if(allocated(connect_neigh))    deallocate(connect_neigh)
    
    !find out which plaques have one edge in the basic unit-cell
    Call get_plaq_base(pair_dat,plaq_base,used_con)

    !get all connections required for each basic plaque
    Call get_plaq_base_connect(pair_dat,plaq_base,used_con,connect_base,connect_neigh)
    deallocate(used_con)

    !unfold the plaques (i.e. connections) to the entire supercell obeying the periodicity
    Nplaq_base=size(connect_neigh,2)
    allocate(N_entry(Nplaq_base),source=0)
    Call unfold_plaq_connect(lat,Nplaq_base,plaq_base,connect_base,con_unfold,N_entry)
    deallocate(plaq_base,connect_base)
end subroutine

subroutine unfold_plaq_connect(lat,Nplaq,plaq_base,base,con_unfold,N_entry)
    type(lattice),intent(in)            :: lat
    integer,intent(in)                  :: Nplaq
    integer,intent(in)                  :: plaq_base(4,4,Nplaq) !corner sites of each plaque (iat/x/y/z,edge1/2/3/4,plaque)
    integer,intent(in)                  :: base(4,6,Nplaq)      !(4,6,Nplaq) site (at,pos) of base point of connection
    integer,intent(inout),allocatable   :: con_unfold(:,:)      !(6,Nplaq_unf) unfolded base points for all plaques
    integer,intent(out)                 :: N_entry(Nplaq)       !number of unfolded plaques for each base plaque

    integer     :: Nplaq_unf
    integer     :: Ncell
    integer     :: iplaq
    integer     :: ext_diff(3,2)
    integer     :: imax(2,3)
    integer     :: i1,i2,i3
    integer     :: i_vec(3)
    integer     :: ii

    integer,allocatable ::  site3_unfold(:,:,:)
    integer,allocatable ::  site_unfold(:,:)

    Ncell=lat%ncell
    allocate(site3_unfold(3,6,Ncell),source=0)
    allocate(site_unfold(6,Ncell*Nplaq),source=0)
    N_entry=0
    do iplaq=1,Nplaq
        !get boundaries for included plaques in superlattice
        ext_diff(:,1)=minval(plaq_base(2:4,:,iplaq),2)
        ext_diff(:,2)=maxval(plaq_base(2:4,:,iplaq),2)
        imax=reshape([1,lat%dim_lat(1),1,lat%dim_lat(2),1,lat%dim_lat(3)],[2,3])
        if(.not.lat%periodic(1)) imax(:,1)=[1-ext_diff(1,1), lat%dim_lat(1)-ext_diff(1,2)]
        if(.not.lat%periodic(2)) imax(:,2)=[1-ext_diff(2,1), lat%dim_lat(2)-ext_diff(2,2)]
        if(.not.lat%periodic(3)) imax(:,3)=[1-ext_diff(3,1), lat%dim_lat(3)-ext_diff(3,2)]

        ii=0
        do i3=imax(1,3),imax(2,3)
            i_vec(3)=i3
            do i2=imax(1,2),imax(2,2)
                i_vec(2)=i2
                do i1=imax(1,1),imax(2,1)
                    i_vec(1)=i1
                    ii=ii+1
                    site3_unfold(:,:,ii)=spread(i_vec,2,6)+base(2:4,:,iplaq)
                enddo
            enddo
        enddo
        N_entry(iplaq)=ii
        Call lat%fold_3_arr(6*ii,site3_unfold)
        Call lat%index_3_1_arr(6*ii,site3_unfold,site_unfold(:,sum(N_entry(:iplaq-1))+1:sum(N_entry)))
    enddo
    deallocate(site3_unfold)

    Nplaq_unf=sum(N_entry)
    allocate(con_unfold(6,Nplaq_unf),source=site_unfold(:,1:Nplaq_unf))
end subroutine

subroutine get_plaq_base_connect(pair_dat,plaq_base,used_con,connect_base,connect_neigh)
    !subroutine which returns the 6 connections constituting each plaque encoded in connect_neigh and connect_base
    !the connections are:
    !1: site 1 -> site 2
    !2: site 4 -> site 2
    !3: site 1 -> site 4
    !4: site 2 -> site 3
    !5: site 1 -> site 3
    !6: site 2 -> site 4
    !, where site 1 -> site 2 -> site 3 -> site3 -> site 1 (or the inverse direction) goes along the edges of the plaque
    !The input contains the used connections (used_con) which contain the basic connections which give the edges of the plaq
    !site 1 -> site 2 : connection 1
    !site 2 -> site 3 : connection 2
    !site 3 -> site 4 :-connection 1  (mainly use site 4 -> site 3 : connection 1)
    !site 4 -> site 1 :-connection 2  (mainly use site 1 -> site 4 : connection 2)

    type(pair_dat_t),intent(in)         :: pair_dat             !necessary information about neighbors
    integer,intent(inout),allocatable   :: plaq_base(:,:,:)     !encodes atom sites which contintute plaques(4,4,N_plaq)
                                                                !(x,:,:) gives [atom_id,site_x,site_y,site_z] (with sites relative to origin)
                                                                !(:,x,:) gives the 4 different sites constituting the plaque
                                                                !(:,:,x) gives the different plaques that can be constructed from the origin
    integer,intent(in)                  :: used_con(:,:)        !(2,Nplaq), which connections are used 
    integer,allocatable,intent(inout)   :: connect_neigh(:,:)   !(6,Nplaq) index of the 6 different neighbor steps in space of pair_dat%pairs(:,X)
    integer,allocatable,intent(inout)   :: connect_base(:,:,:)  !(4,6,Nplaq) site (at,pos) of base point of connection

    integer             :: pos(4,4)     !temporary array saving the position of the constituting sites for a plaque
    integer             :: Nplaq, iplaq !number of plaques found
    integer             :: Ncon         !total number of nearest neighbor connections starting at an atom
    integer             :: Nbase        !half number of nearest neighbor connections (to exclude the negative pairs)

    integer             :: N1,N2        !total number of nearest/next nearest connections
    integer,allocatable :: mat_nnn(:,:) !(N1,N1) gives index of            next nearest neighbor interaction when combining 2 nearest neighbor interactions encoded by their index (otherwise 0)
    integer,allocatable :: mat_nn(:,:)  !(N1,N1) gives index of                 nearest neighbor interaction when combining 2 nearest neighbor interactions encoded by their index (otherwise 0)
    integer,allocatable :: mat(:,:)     !(N1,N1) gives index of nearest of next nearest neighbor interaction when combining 2 nearest neighbor interactions encoded by their index (otherwise 0)
    integer             :: con1,con2    !temporary values for connections

    Nplaq=size(plaq_base,3)
    Ncon=pair_dat%Npair_at(1)
    Nbase=Ncon/2

    !get matrix which contains the connection index for combining 2 nearest neighbor connections
    N1=pair_dat%Nshell(1)
    N2=pair_dat%Nshell(2)
    allocate(mat_nnn(N1,N1),mat_nn(N1,N1),source=0)
    Call get_neigh_arr_nnn(N1,N2,pair_dat%pairs(:,1:N1),pair_dat%pairs(:,N1+1:N1+N2),mat_nnn)
    Call get_neigh_arr_nn(N1,pair_dat%pairs(:,1:N1),mat_nn)
    allocate(mat,source=mat_nn)
    mat=mat+mat_nnn
    deallocate(mat_nn,mat_nnn)

    !fill required data
    allocate(connect_neigh(6,Nplaq))
    allocate(connect_base(4,6,Nplaq))
    do iplaq=1,Nplaq
        !help variable (for readability)
        pos(:,1)=plaq_base(:,1,iplaq)
        pos(:,2)=plaq_base(:,2,iplaq)
        pos(:,3)=plaq_base(:,3,iplaq)
        pos(:,4)=plaq_base(:,4,iplaq)
        !initial sites for respective connection:
        connect_base(:,1,iplaq)=pos(:,1)
        connect_base(:,2,iplaq)=pos(:,4)
        connect_base(:,3,iplaq)=pos(:,1)
        connect_base(:,4,iplaq)=pos(:,2)
        connect_base(:,5,iplaq)=pos(:,1)
        connect_base(:,6,iplaq)=pos(:,2)
        !easy nearest neighbor connections along the edge:
        connect_neigh(1,iplaq)=used_con(1,iplaq)+Ncon*(connect_base(1,1,iplaq)-1)  !connection site 1-> site 2
        connect_neigh(2,iplaq)=used_con(1,iplaq)+Ncon*(connect_base(1,2,iplaq)-1)  !connection site 4-> site 3
        connect_neigh(3,iplaq)=used_con(2,iplaq)+Ncon*(connect_base(1,3,iplaq)-1)  !connection site 1-> site 4
        connect_neigh(4,iplaq)=used_con(2,iplaq)+Ncon*(connect_base(1,4,iplaq)-1)  !connection site 2-> site 3
        !more complicated connections crossing the plaque:
        con1=used_con(1,iplaq)+Ncon*(pos(1,1)-1)    !connection site 1-> site 2
        con2=used_con(2,iplaq)+Ncon*(pos(1,2)-1)    !connection site 2-> site 3
        connect_neigh(5,iplaq)=mat(con1,con2)       !combine both connections to site 1 -> site 3
        con1=mod(used_con(1,iplaq)-1+Nbase,Ncon)+1  +Ncon*(pos(1,2)-1)   !connection site 2-> site 1  (mod-stuff: get negative connection 1 from site 2)
        con2=used_con(2,iplaq)+Ncon*(pos(1,1)-1)    !connection site 1-> site 4
        connect_neigh(6,iplaq)=mat(con1,con2)       !combine both connections to site 2 -> site 4
    enddo
    !some checks
    if(any(connect_neigh(5,:)==0)) ERROR STOP "Trying failed to get one of the connections from site 1 to site 3"
    if(any(connect_neigh(6,:)==0)) ERROR STOP "Trying failed to get one of the connections from site 2 to site 4"
    if(any(connect_neigh     ==0)) ERROR STOP "Unexpected fail fo getting neighbors"
end subroutine                                                

subroutine get_plaq_base(pair_dat,plaq_base,used_con)
    !gets all infomation about which plaques have a site in the basic unit-cell
    type(pair_dat_t),intent(in)         :: pair_dat          !necessary information about neighbors
    integer,intent(inout),allocatable   :: plaq_base(:,:,:)  !encodes atom sites which contintute plaques(4,4,N_plaq)
                                                             !(x,:,:) gives [atom_id,site_x,site_y,site_z] (with sites relative to origin)
                                                             !(:,x,:) gives the 4 different sites constituting the plaque
                                                             !(:,:,x) gives the different plaques that can be constructed from the origin
    integer,intent(inout),allocatable   :: used_con(:,:)     !(2,Nplaq), which connections are used 
    integer     :: plaq_tmp(4,4)        !temporary plaq indices
    integer     :: Nat, iat             !different atoms per unit-cell which constitute the plaques
    integer     :: Ncon, icon           !number of nearest neighbor connections
    integer     :: Nbase, ibase         !half number of nearest neighbor connections (to exclude the negative pairs)
    integer     :: Nplaq, iplaq         !number of plaques found
    integer,allocatable ::  bnd_nn(:,:) !(2:Nat) boundaries of nearest neighbor connections of each atoms in space of pair_dat%pairs(1:5,x)
    
    !initial parameters
    Nat=pair_dat%Nat
    Ncon=pair_dat%Npair_at(1)
    Nbase=Ncon/2
    Nplaq=Nat*Nbase*(Nbase-1)*2
    allocate(plaq_base(4,4,Nplaq),source=0)
    allocate(used_con(2,Nplaq),source=0)
    allocate(bnd_nn(2,Nat))
    do iat=1,Nat
        bnd_nn(1,iat)=sum(pair_dat%Npair_at(:iat-1))+1
        bnd_nn(2,iat)=sum(pair_dat%Npair_at(:iat))
    enddo

    !check that the different nearest neighbors are all correctly ordered and equivalent
    Call check_pair(pair_dat,Nat,Ncon,Nbase)    
   
    !fill the different basis plaques
    iplaq=0
    do iat=1,Nat
        plaq_tmp(:,1)=[iat,0,0,0]
        do ibase=1,Nbase
            plaq_tmp(:,2)=pair_dat%pairs(2:,bnd_nn(1,iat)-1+ibase)
            do icon=1,Ncon
                if(mod(icon-1,Nbase)+1<=ibase) cycle
                plaq_tmp(:,3)=pair_dat%pairs(2:,bnd_nn(1,plaq_tmp(1,2))-1+icon)
                plaq_tmp(2:4,3)=plaq_tmp(2:4,3)+plaq_tmp(2:4,2)
                plaq_tmp(:,4)=pair_dat%pairs(2:,bnd_nn(1,plaq_tmp(1,1))-1+icon)
                iplaq=iplaq+1
                used_con(:,iplaq)=[ibase,icon]
                plaq_base(:,:,iplaq)=plaq_tmp
            enddo
            plaq_tmp(:,2)=pair_dat%pairs(2:,bnd_nn(1,iat)-1+ibase+Nbase)
            do icon=1,Ncon
                if(mod(icon-1,Nbase)+1<=ibase) cycle
                plaq_tmp(:,3)=pair_dat%pairs(2:,bnd_nn(1,plaq_tmp(1,2))-1+icon)
                plaq_tmp(2:4,3)=plaq_tmp(2:4,3)+plaq_tmp(2:4,2)
                plaq_tmp(:,4)=pair_dat%pairs(2:,bnd_nn(1,plaq_tmp(1,1))-1+icon)
                iplaq=iplaq+1
                used_con(:,iplaq)=[ibase+Nbase,icon]
                plaq_base(:,:,iplaq)=plaq_tmp
            enddo
        enddo
    enddo
end subroutine

subroutine check_pair(pair_dat,Nat,Ncon,Nbase)
    !some checks to make sure that the neighboring infomation is supplied as it should be (later asumptions of setup)
    type(pair_dat_t),intent(in) :: pair_dat
    integer,intent(in)      ::  Nat,Ncon,Nbase
    integer ::  iat,ibase, icon,i1,i2

    if(maxval(pair_dat%Npair_at(:Nat))/=minval(pair_dat%Npair_at(:Nat))) STOP "maxval and minval Npair_at differ, probably the structure does not make sense in case of 4-spin interaction"
    do iat=1,Nat
        do ibase=1,Nbase
            if(norm2(pair_dat%diff_vec(:,(iat-1)*Ncon+ibase)+pair_dat%diff_vec(:,(iat-1)*Ncon+ibase+Nbase))>1.0d-6*norm2(pair_dat%diff_vec(:,(iat-1)*Ncon+ibase)))then
                STOP "nearest neighbor connecting vectors do not include negative element or are orderer incorrectly"
            endif
        enddo
    enddo
    do i1=1,Nat
        do i2=1,Nat
            do icon=1,Ncon
                if(norm2(pair_dat%diff_vec(:,(i1-1)*Ncon+icon)-pair_dat%diff_vec(:,(i2-1)*Ncon+icon))>1.0d-6*norm2(pair_dat%diff_vec(:,(i1-1)*Ncon+icon)))then
                    STOP "nearest neighbors of different atoms are not ordered in the same way"
                endif
            enddo
        enddo
    enddo
end subroutine

subroutine get_neigh_arr_nn(N1,pairs1,mat)
    !get matrix which contains index of nearest neighbor interaction by combining 2 nearest neighbor interactions according to their indices
    !super inefficient implementation which will scale terribly with the number of atoms in the unit-cell, but it is not called very often...
    integer,intent(in)      :: N1
    integer,intent(in)      :: pairs1(5,N1)
    integer,intent(out)     :: mat(N1,N1) 

    integer     :: i,j,l
    integer     :: ind_sum(3)

    mat=0
    do i=1,N1
        do j=1,N1
            if(pairs1(2,i)/=pairs1(1,j)) cycle
            do l=1,N1
                if(pairs1(1,i)/=pairs1(1,l)) cycle
                if(pairs1(2,j)/=pairs1(2,l)) cycle
                ind_sum=pairs1(3:,i)+pairs1(3:,j)
                if(any(ind_sum/=pairs1(3:,l))) cycle
                mat(i,j)=l
            enddo
        enddo
    enddo
end subroutine

subroutine get_neigh_arr_nnn(N1,N2,pairs1,pairs2,mat)
    !get matrix which contains index of next nearest neighbor interaction by combining 2 nearest neighbor interactions according to their indices
    !super inefficient implementation which will scale terribly with the number of atoms in the unit-cell, but it is not called very often...
    integer,intent(in)      :: N1,N2
    integer,intent(in)      :: pairs1(5,N1)
    integer,intent(in)      :: pairs2(5,N2)
    integer,intent(out)     :: mat(N1,N1) 

    integer     :: i,j,l
    integer     :: ind_sum(3)

    mat=0
    do i=1,N1
        do j=1,N1
            if(pairs1(2,i)/=pairs1(1,j)) cycle
            do l=1,N2
                if(pairs1(1,i)/=pairs2(1,l)) cycle
                if(pairs1(2,j)/=pairs2(2,l)) cycle
                ind_sum=pairs1(3:,i)+pairs1(3:,j)
                if(any(ind_sum/=pairs2(3:,l))) cycle
                mat(i,j)=l+N1
            enddo
        enddo
    enddo
end subroutine
end module
