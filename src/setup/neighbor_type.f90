module m_neighbor_type
use m_derived_types
implicit none
private
public :: neighbors, get_neigh_distances


type :: neighbors
    integer,allocatable ::  Nshell(:)   !how many different connections are there for each distance  (1:number distance)
    integer,allocatable ::  ishell(:)   !how many entries are there for each shell (1:number shells)
    integer,allocatable ::  at_pair(:,:)    !what are the at_pairs for the shell (2,1:number shells)
    integer,allocatable ::  pairs(:,:)  !what cell pairs in (1:Ncell) for the pair (2,1:number connections)

    contains
    procedure   ::  get   => get_neighbors
    procedure   ::  unset 
end type

contains
subroutine unset(neigh)
    class(neighbors),intent(inout) :: neigh

    if(allocated(neigh%Nshell )) deallocate(neigh%Nshell )
    if(allocated(neigh%ishell )) deallocate(neigh%ishell )
    if(allocated(neigh%at_pair)) deallocate(neigh%at_pair)
    if(allocated(neigh%pairs  )) deallocate(neigh%pairs  )
end subroutine

subroutine get_neighbors(neigh,atid,neighval,lat)
    class(neighbors),intent(out)    :: neigh
    integer,intent(in)              :: atid(2)      !atom type index of the atom types whose neighbors shall be found
    integer,intent(in)              :: neighval(:)  !every entry means that this neighbor is requires (2= second nearest neighbor etc.)
    type(lattice),intent(in)        :: lat          !the all-knowing lattice providing all required information

    integer,allocatable     :: id1(:),id2(:)
    integer                 :: Nat(2)
    real(8),allocatable     :: atpos1(:,:),atpos2(:,:)
    integer                 :: i
    real(8)                 :: distance(size(neighval))
    integer,allocatable     :: pairs(:,:)
    integer,allocatable     :: Nshell(:)

    Call lat%cell%ind_attype(atid(1),id1)
    Call lat%cell%ind_attype(atid(2),id2)
    Nat=[size(id1),size(id2)]
    allocate(atpos1(3,Nat(1)),atpos2(3,Nat(2)),source=0.0d0)
    do i=1,Nat(1)
        atpos1(:,i)=lat%cell%atomic(id1(i))%position
    enddo
    do i=1,Nat(2)
        atpos2(:,i)=lat%cell%atomic(id2(i))%position
    enddo
    !get all neighbor connections(pair_ind,Nshell) and the distances
    Call get_neigh_distances(atpos1,atpos2,neighval,lat,pairs,Nshell,distance)

    !set pair_ind atom indices to real cell atomic indices
    pairs(1,:)=id1(pairs(1,:))
    pairs(2,:)=id2(pairs(2,:))

    !unfold the neighbors to all lattice sites obeying periodicity and fill neighbors-type (neigh)
    Call unfold_neighbors(pairs,Nshell,lat,neigh)

end subroutine

subroutine unfold_neighbors(pairs_in,Nshell_in,lat,neigh_out)
    type(lattice),intent(in)    :: lat          !the all-knowing lattice providing all required information
    integer,intent(in)          :: pairs_in(:,:)
    integer,intent(in)          :: Nshell_in(:)
    type(neighbors),intent(out) :: neigh_out

    integer                     :: pairs(2,size(pairs_in,2)*lat%Ncell)
    integer                     :: ishell(sum(Nshell_in)),at_pair(2,sum(Nshell_in))

    integer                     :: ii
    integer                     :: i_shell,i_pair_in,i_pair_done
    integer                     :: i_1,i_2  !lattice sites in [1,Ncell] basis
    integer                     :: i3_1(3),i3_2(3)  !lattice sites in ([1,dimlat(1)],[1,dimlat(2)],[1,dimlat(3)])-basis
    integer                     :: i1,i2,i3
    integer                     :: imax(2,3)

    ii=0
    i_pair_done=0
    do i_shell=1,size(Nshell_in)
        do i_pair_in=i_pair_done+1,i_pair_done+Nshell_in(i_shell)
            at_pair(:,i_pair_in)=pairs_in(1:2,i_pair_in)
            imax=reshape([1,lat%dim_lat(1),1,lat%dim_lat(2),1,lat%dim_lat(3)],[2,3])
            if(.not.lat%periodic(1).and.pairs_in(2+1,i_pair_in)<0) imax(1,1)=1             -pairs_in(2+1,i_pair_in)
            if(.not.lat%periodic(1).and.pairs_in(2+1,i_pair_in)>0) imax(2,1)=lat%dim_lat(1)-pairs_in(2+1,i_pair_in)
            if(.not.lat%periodic(2).and.pairs_in(2+2,i_pair_in)<0) imax(1,2)=1             -pairs_in(2+2,i_pair_in)
            if(.not.lat%periodic(2).and.pairs_in(2+2,i_pair_in)>0) imax(2,2)=lat%dim_lat(2)-pairs_in(2+2,i_pair_in)
            if(.not.lat%periodic(3).and.pairs_in(2+3,i_pair_in)<0) imax(1,3)=1             -pairs_in(2+3,i_pair_in)
            if(.not.lat%periodic(3).and.pairs_in(2+3,i_pair_in)>0) imax(2,3)=lat%dim_lat(3)-pairs_in(2+3,i_pair_in)

            do i3=imax(1,3),imax(2,3)
                i3_1(3)=i3
                do i2=imax(1,2),imax(2,2)
                    i3_1(2)=i2
                    do i1=imax(1,1),imax(2,1)
                        i3_1(1)=i1
                        i3_2=modulo(i3_1+pairs_in(3:5,i_pair_in)-1,lat%dim_lat)+1
                        ii=ii+1
                        pairs(1,ii)=lat%index_m_1(i3_1)
                        pairs(2,ii)=lat%index_m_1(i3_2)
                    enddo
                enddo
            enddo
            ishell(i_pair_in)=ii
            i_pair_done=i_pair_done+1
        enddo
    enddo

    allocate(neigh_out%Nshell,source=Nshell_in)
    allocate(neigh_out%ishell,source=ishell)
    allocate(neigh_out%at_pair,source=at_pair)
    allocate(neigh_out%pairs,source=pairs(:,:ii))
end subroutine

subroutine get_neigh_distances(atpos1,atpos2,neighval,lat,pair_ind,N_shell,dist_out)
    !get the neighval's distances between the atoms in atpos1 and atpos2
    use m_sort
    !subroutine which gets the neighval's distances between the atoms in atpos1 and atpos2
    real(8),intent(in)          :: atpos1(:,:)
    real(8),intent(in)          :: atpos2(:,:)
    integer,intent(in)          :: neighval(:)  !every entry means that this neighbor is requires (2= second nearest neighbor etc.)
    type(lattice),intent(in)    :: lat          !the all-knowing lattice providing all required information
    real(8),intent(out)         :: dist_out(size(neighval))
    integer,allocatable     :: pair_ind(:,:) !integer array containing the information how which atoms are corrected by the distance([[ia1,ia2,ix,iy,iz],[:]])
                                             ! distance corresponds to atpos1(ia1)-atpos2(ia2)+lat1*ix+lat2*iy+lat3*iz
    integer,allocatable     :: N_shell(:)   !number of shells per considered neighbor

    integer                 :: Ncheck   !how many unit cells in each direction are checked for neighbors (could definitely be reduced in some cases)
    integer,allocatable     :: ind1(:),ind2(:),ind3(:)
    real(8),allocatable     :: all_dist(:)      !contains all distances (first unsorted, later sorted)
    integer,allocatable     :: all_distint(:)   !integer indices giving original index after sorting of all_dist
    integer                 :: i1,i2,i3,ia,ii,i,j
    real(8)                 :: lat_diff(3,3),diff(3)

    real(8)                 :: at_diff(3,size(atpos1,2)*size(atpos2,2)) 
    real(8)                 :: vmax,eps
    integer                 :: dist_start(maxval(neighval)+1)   !starting index in the sorted all_dist
    integer                 :: orig_ind !original index of considerer connection
    integer                 :: multi(5) !size of each index part
    integer                 :: prod_multi(5) !size of each index part


    ii=0
    do i1=1,size(atpos1,2)
        do i2=1,size(atpos2,2)
            ii=ii+1
            at_diff(:,ii)=atpos1(:,i1)-atpos2(:,i2)
        enddo
   enddo

    Ncheck=maxval(neighval) !this might actually be chosen much smaller, depending of the geometry...
    if(lat%dim_lat(1)>1.or.lat%periodic(1))then
        ind1=[(-Ncheck+i,i=0,2*Ncheck)]
    else
        ind1=[0]
    endif
    if(lat%dim_lat(2)>1.or.lat%periodic(2))then
        ind2=[(-Ncheck+i,i=0,2*Ncheck)]
    else
        ind2=[0]
    endif
    if(lat%dim_lat(3)>1.or.lat%periodic(3))then
        ind3=[(-Ncheck+i,i=0,2*Ncheck)]
    else
        ind3=[0]
    endif
    allocate(all_dist(size(atpos1,2)*size(atpos2,2)*size(ind1)*size(ind2)*size(ind3)),source=0.0d0)
    ii=0
    do i3=1,size(ind3)
        lat_diff(:,3)=lat%areal(3,:)*real(ind3(i3),8)
        do i2=1,size(ind2)
            lat_diff(:,2)=lat%areal(2,:)*real(ind2(i2),8)
            do i1=1,size(ind1)
                lat_diff(:,1)=lat%areal(1,:)*real(ind1(i1),8)
                diff=sum(lat_diff,2)
                do ia=1,size(at_diff,2)
                    ii=ii+1
                    all_dist(ii)=norm2(at_diff(:,ia)+diff)
                enddo
            enddo
       enddo
    enddo

    !sort to get distances easily 
    allocate(all_distint(size(all_dist)),source=0)
    eps=minval(norm2(lat%areal,1))*1.0d-8
    Call sort(size(all_dist),all_dist,all_distint,eps)

    !get starting indices of new distance set in sorted distances
    ii=0; i=0
    vmax=0.0d0 !do not include zero distances if considering same atom types
    do while(ii<size(dist_start))
        i=i+1
        if(all_dist(i)>vmax+eps)then
            ii=ii+1
            dist_start(ii)=i
            vmax=all_dist(i)
        endif
    enddo

    !fill the distance array and number of shells found for each distance 
    dist_out=[(all_dist(dist_start(neighval(i))),i=1,size(neighval))]
    N_shell=[(dist_start(neighval(i)+1)-dist_start(neighval(i)),i=1,size(neighval))]

    !get pair_ind by translating the relevant entries from the all_distint array to the atom indices and lattice displacements
    allocate(pair_ind(5,sum(N_shell)),source=0)
    ii=0
    multi=[size(atpos2,2), size(atpos1,2), size(ind1), size(ind2), size(ind3)]
    prod_multi=[(product(multi(1:i-1)),i=1,5)]
    do i=1,size(neighval)
        ia=dist_start(neighval(i))
        do j=0,N_shell(i)-1
            ii=ii+1
            orig_ind=all_distint(ia+j)
            pair_ind(:,ii)=modulo((orig_ind-1)/prod_multi,multi)+1
        enddo
    enddo
    do i=1,size(pair_ind,2)!translate lattice displacements to be correct with respect to origin
        pair_ind(3,i)=ind1(pair_ind(3,i))
        pair_ind(4,i)=ind2(pair_ind(4,i))
        pair_ind(5,i)=ind3(pair_ind(5,i))
    enddo

    !write(*,*) N_shell
    !do i=1,size(pair_ind,2)
    !    write(*,*) pair_ind(:,i), norm2(atpos1(:,pair_ind(1,i))-atpos2(:,pair_ind(2,i))+matmul(pair_ind(3:5,i),lat%areal))
    !enddo
end subroutine
end module
