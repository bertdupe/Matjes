module m_mapping
use m_get_position
use m_arrange_neigh, only : get_map,fill_map

interface mapping_geometry
    module procedure mapping_1D
    module procedure mapping_2D
    module procedure mapping_2D_SL
    module procedure mapping_2D_motif
    module procedure mapping_3D
    module procedure mapping_3D_motif_SL
end interface mapping_geometry

private
public :: mapping
contains

subroutine mapping(tabledist,N_Nneigh,my_motif,indexNN,tableNN,my_lattice)
use m_derived_types, only : lattice,cell
implicit none
integer, intent(inout) :: tableNN(:,:,:,:,:,:)
type(lattice), intent(in) :: my_lattice
type(cell), intent(in) :: my_motif
real(kind=8), intent(in) :: tabledist(:,:)
integer, intent(in) :: N_Nneigh,indexNN(:,:)
! internal
integer :: phase,Nei_il,Nei_z,shape_tabdist(2),i,j

shape_tabdist=shape(tabledist)
phase=1
!
! take care of the mapping
!
call get_map(indexNN)

if (shape_tabdist(2).gt.1) phase=2

if (size(my_lattice%world).eq.1) then

    call mapping_geometry(tabledist(:,1),N_Nneigh,indexNN(:,1),tableNN(:,:,:,1,1,1),my_lattice)

elseif (size(my_lattice%world).eq.2) then
    if ((phase.eq.1).and.(count(my_motif%atomic(:)%moment.gt.0.0d0).eq.1)) then
        call mapping_geometry(tabledist(:,1),N_Nneigh,indexNN(:,1),tableNN(:,:,:,:,1,1),my_lattice)

    elseif (phase.eq.1) then
        call mapping_geometry(tabledist(:,1),N_Nneigh,indexNN(:,1),tableNN(:,:,:,:,1,:),my_lattice)

    endif

    if (phase.eq.2) then
        Nei_il=count(tabledist(:,2).ne.0)
        call mapping_geometry(tabledist(:,:),N_Nneigh,Nei_il,indexNN(:,:),tableNN(:,:,:,:,1,:),my_lattice)

    endif

else
     if (phase.eq.1) then
        call mapping_geometry(tabledist(:,1),N_Nneigh,indexNN(:,1),tableNN(:,:,:,:,:,:),my_lattice)

     else
        Nei_il=count(tabledist(:,2).ne.0)
        Nei_z=count(tabledist(:,3).ne.0)
        call mapping_geometry(tabledist(:,:),N_Nneigh,Nei_il,Nei_z,indexNN(:,:),tableNN(:,:,:,:,:,:),my_lattice)

     endif

endif

call write_out(tableNN,indexNN)

end subroutine mapping

! establish the table of neighbors
! it gives the indices of the neighbors of spin i_s
! mapping is a table which first contains
! first entry i_s: key of the cell
! second n is the neighbor index written on a line
! third column gives: 1-x;2-y;3-z position
subroutine mapping_1D(d,nei,indexNN,tableNN,my_lattice)
use m_vector , only : norm
use m_derived_types, only : lattice
#ifdef CPP_MPI
use m_make_box, only : Xstart
#endif
implicit none
! variable that come in
integer, intent(in) :: nei
type(lattice), intent(in) :: my_lattice
real(kind=8), intent(in) :: d(:)
integer, intent(in) :: indexNN(:)
! value of the function
integer, intent(inout) :: tableNN(:,:,:)
! external blas
! 3D coordinate ix,iy,iz of 1d coordinate k
integer :: i_x,Xstop
! dummy variable
integer :: i,l,i_Nei,avant,dim_lat(3)
integer :: v_x,ok
real (kind=8) :: vec(3),dist,r(3,3)
logical :: found

avant=0
Xstop=size(tableNN,3)
do i=1,3
  r(:,i)=my_lattice%areal(i,:)
enddo
dim_lat=my_lattice%dim_lat

do i_nei=1,nei
    do i_x=1,Xstop

       l=1
       do i=i_nei,-i_nei,-1

          vec=0.0d0
          found=.false.
          ! suppose that the neighbour should be taken into account
          ok=1
          call test_neighbour(vec,v_x,ok,i_x,i,Xstop,r(:,1),my_lattice%boundary(1))
          dist=norm(vec)
          call associate_neighbour(tableNN(:,avant+l,i_x),v_x,1,1,1,ok,l,d(i_nei),dist,found)

          if ((i_x.eq.1).and.(found)) call fill_map(i_nei,l-1,vec)

       enddo

    enddo
    avant=avant+indexNN(i_Nei)
enddo

end subroutine mapping_1D

!subroutine to treat the neighbour with one atom in the 2D unit cell
!
subroutine mapping_2D(d,nei,indexNN,tableNN,my_lattice)
use m_vector , only : norm
use m_derived_types, only : lattice
#ifdef CPP_MPI
use m_make_box, only : Xstart,Ystart
#endif
#ifdef CPP_OPENMP
use omp_lib
#endif
implicit none
integer, intent(in) :: nei
type(lattice), intent(in) :: my_lattice
real(kind=8), intent(in) :: d(:)
integer, intent(in) :: indexNN(:)
! value of the function
integer, intent(inout) :: tableNN(:,:,:,:)
! external blas
! 3D coordinate ix,iy,iz of 1d coordinate k
integer :: i_x,i_y,Xstop,Ystop
! size of the tableNN
integer :: shape_tableNN(4)
! dummy variable
integer :: i,j,l,i_Nei,avant,dim_lat(3)
integer :: v_x,v_y,ok
real (kind=8) :: vec(3),dist,r(3,3)
logical :: found
#ifndef CPP_MPI
integer, parameter ::  Xstart=1
integer, parameter ::  Ystart=1
#endif
#ifdef CPP_OPENMP
integer :: ithread,nthreads

nthreads=omp_get_num_procs()
call omp_set_num_threads(nthreads)
#endif

avant=0
shape_tableNN=shape(tableNN)
Xstop=shape_tableNN(3)
Ystop=shape_tableNN(4)
dim_lat=my_lattice%dim_lat
do i=1,3
  r(:,i)=my_lattice%areal(i,:)
enddo

do i_nei=1,nei

   do i_y=1,Ystop
      do i_x=1,Xstop

         l=1

         do i=i_nei,-i_nei,-1
            do j=i_nei,-i_nei,-1

               vec=0.0d0
               found=.false.
               ! suppose that the neighbour should be taken into account
               ok=1
               call test_neighbour(vec,v_x,ok,i_x,i,Xstop,r(:,1),my_lattice%boundary(1))
               call test_neighbour(vec,v_y,ok,i_y,j,Ystop,r(:,2),my_lattice%boundary(2))
               dist=norm(vec)

               call associate_neighbour(tableNN(:,avant+l,i_x,i_y),v_x,v_y,1,1,ok,l,d(i_nei),dist,found)
               call check_l(l-1,indexNN(i_Nei))

               if ((i_x.eq.1).and.(i_y.eq.1).and.(found)) call fill_map(i_nei,l-1,vec)

            enddo
         enddo

      enddo
   enddo

   avant=avant+indexNN(i_Nei)

enddo

end subroutine mapping_2D

!subroutine to treat the SL case
!
subroutine mapping_2D_SL(d,nei,Nei_il,indexNN,tableNN,my_lattice)
use m_vector , only : norm
use m_derived_types, only : lattice
#ifdef CPP_MPI
use m_make_box, only : Xstart,Ystart
#endif
implicit none
integer, intent(in) :: nei,Nei_il
type(lattice), intent(in) :: my_lattice
real(kind=8), intent(in) :: d(:,:)
integer, intent(in) :: indexNN(:,:)
! value of the function
integer, intent(inout) :: tableNN(:,:,:,:,:)
! external blas
! 3D coordinate ix,iy,iz of 1d coordinate k
integer :: i_x,i_y,i_m,Xstop,Ystop,Mstop
integer :: v_x,v_y,ok
! dummy variable
integer :: i,j,l,i_Nei,avant,i_p,i_phase,dim_lat(3)
real (kind=8) :: vec(3),dist,r(3,3)
logical :: found
#ifndef CPP_MPI
integer, parameter ::  Xstart=1
integer, parameter ::  Ystart=1
#endif

avant=0
Ystop=size(tableNN,4)
Xstop=size(tableNN,3)
Mstop=size(tableNN,5)
do i=1,3
  r(:,i)=my_lattice%areal(i,:)
enddo
dim_lat=my_lattice%dim_lat

i_phase=1

do i_nei=1,nei
   do i_y=1,Ystop
      do i_x=1,Xstop
         do i_m=1,Mstop

            l=1
            do i=i_nei,-i_nei,-1
               do j=i_nei,-i_nei,-1
                  do i_p=1,Mstop

                     vec=0.0d0
                     found=.false.
                     ! suppose that the neighbour should be taken into account
                     ok=1

                     call test_neighbour(vec,v_x,ok,i_x,i,Xstop,r(:,1),my_lattice%boundary(1))
                     call test_neighbour(vec,v_y,ok,i_y,j,Ystop,r(:,2),my_lattice%boundary(2))
                     dist=norm(vec)

                     call associate_neighbour(tableNN(:,avant+l,i_x,i_y,i_m),v_x,v_y,1,i_p,ok,l,d(i_nei,i_phase),dist,found)

                     if ((i_x.eq.1).and.(i_y.eq.1).and.(i_m.eq.1).and.(found)) call fill_map(i_nei,l-1,vec)

                  enddo
               enddo
            enddo

         enddo
      enddo
   enddo
   avant=avant+indexNN(i_Nei,i_phase)
enddo

i_phase=2

do i_nei=1,Nei_il
   do i_y=1,Ystop
      do i_x=1,Xstop
         do i_m=1,Mstop

            l=1
            do i=i_nei,-i_nei,-1
               do j=i_nei,-i_nei,-1
                  do i_p=1,Mstop

                     vec=0.0d0
                     found=.false.
                     ! suppose that the neighbour should be taken into account
                     ok=1

                     call test_neighbour(vec,v_x,ok,i_x,i,Xstop,r(:,1),my_lattice%boundary(1))
                     call test_neighbour(vec,v_y,ok,i_y,j,Ystop,r(:,2),my_lattice%boundary(2))
                     dist=norm(vec)

                     call associate_neighbour(tableNN(:,avant+l,i_x,i_y,i_m),v_x,v_y,1,i_p,ok,l,d(i_nei,i_phase),dist,found)

                     if ((i_x.eq.1).and.(i_y.eq.1).and.(i_m.eq.1).and.(found)) call fill_map(i_nei,l-1,vec)

                  enddo
               enddo
            enddo

         enddo
      enddo
   enddo
   avant=avant+indexNN(i_Nei,i_phase)
enddo

end subroutine mapping_2D_SL

!subroutine to treat more than one atom in the unit cell
!
subroutine mapping_2D_motif(d,nei,indexNN,tableNN,my_lattice)
use m_vector , only : norm
use m_derived_types, only : lattice
#ifdef CPP_MPI
use m_make_box, only : Xstart,Ystart
#endif
implicit none
integer, intent(in) :: nei
type(lattice), intent(in) :: my_lattice
real(kind=8), intent(in) :: d(:)
integer, intent(in) :: indexNN(:)
! value of the function
integer, intent(inout) :: tableNN(:,:,:,:,:)
! external blas
! 3D coordinate ix,iy,iz of 1d coordinate k
integer :: i_x,i_y,i_m,Xstop,Ystop,Mstop
integer :: v_x,v_y,ok
! dummy variable
integer :: i,j,l,i_Nei,avant,i_p,dim_lat(3)
real (kind=8) :: vec(3),dist,r(3,3)
logical :: found
#ifndef CPP_MPI
integer, parameter ::  Xstart=1
integer, parameter ::  Ystart=1
#endif

avant=0
Xstop=size(tableNN,3)
Ystop=size(tableNN,4)
Mstop=size(tableNN,5)
dim_lat=my_lattice%dim_lat
do i=1,3
  r(:,i)=my_lattice%areal(i,:)
enddo

do i_nei=1,nei
   do i_y=1,Ystop
      do i_x=1,Xstop
         do i_m=1,Mstop

            l=1
            do i=i_nei,-i_nei,-1
               do j=i_nei,-i_nei,-1
                  do i_p=1,Mstop

                     vec=0.0d0
                     found=.false.
                     ! suppose that the neighbour should be taken into account
                     ok=1

                     call test_neighbour(vec,v_x,ok,i_x,i,Xstop,r(:,1),my_lattice%boundary(1))
                     call test_neighbour(vec,v_y,ok,i_y,j,Ystop,r(:,2),my_lattice%boundary(2))
                     dist=norm(vec)

                     call associate_neighbour(tableNN(:,avant+l,i_x,i_y,i_m),v_x,v_y,1,i_p,ok,l,d(i_nei),dist,found)

                     if ((i_x.eq.1).and.(i_y.eq.1).and.(i_m.eq.1).and.(found)) call fill_map(i_nei,l-1,vec)

                  enddo
               enddo
            enddo

         enddo
      enddo
   enddo
   avant=avant+indexNN(i_Nei)
enddo


end subroutine mapping_2D_motif

subroutine mapping_3D(d,nei,indexNN,tableNN,my_lattice)
use m_vector , only : norm
use m_derived_types, only : lattice
#ifdef CPP_MPI
use m_make_box, only : Xstart,Ystart,Zstart
#endif
implicit none
! variable that come in
integer, intent(in) :: nei
type(lattice), intent(in) :: my_lattice
real(kind=8), intent(in) :: d(:)
integer, intent(in) :: indexNN(:)
! value of the function
integer :: tableNN(:,:,:,:,:,:)
! external blas
! 3D coordinate ix,iy,iz of 1d coordinate k
 integer :: i_x,i_y,i_z,i_m,Xstop,Ystop,Zstop,Mstop
 integer :: v_x,v_y,v_z,ok
! dummy variable
integer :: i,j,k,l,i_Nei,avant,i_p,dim_lat(3)
real (kind=8) :: vec(3),dist,r(3,3)
logical :: found
#ifndef CPP_MPI
integer, parameter ::  Xstart=1
integer, parameter ::  Ystart=1
integer, parameter ::  Zstart=1
#endif

avant=0
Xstop=size(tableNN,3)
Ystop=size(tableNN,4)
Zstop=size(tableNN,5)
Mstop=size(tableNN,6)
dim_lat=my_lattice%dim_lat
do i=1,3
  r(:,i)=my_lattice%areal(i,:)
enddo

do i_nei=1,nei
   do i_z=1,Zstop
      do i_y=1,Ystop
         do i_x=1,Xstop
            do i_m=1,Mstop

! I selected a cell (i_x,i_y,i_z) and inside one of the atom i_m
! now we have to go along x,y and z direction and check all the distances
               l=1
               do i=i_nei,-i_nei,-1
                  do j=i_nei,-i_nei,-1
                     do k=i_nei,-i_nei,-1
                        do i_p=1,Mstop

                           vec=0.0d0
                           found=.false.
                           ! suppose that the neighbour should be taken into account
                           ok=1

                           call test_neighbour(vec,v_x,ok,i_x,i,Xstop,r(:,1),my_lattice%boundary(1))
                           call test_neighbour(vec,v_y,ok,i_y,j,Ystop,r(:,2),my_lattice%boundary(2))
                           call test_neighbour(vec,v_z,ok,i_z,k,Zstop,r(:,3),my_lattice%boundary(3))
                           dist=norm(vec)

                           call associate_neighbour(tableNN(:,avant+l,i_x,i_y,i_z,i_m),v_x,v_y,v_z,i_p,ok,l,d(i_nei),dist,found)

                           if ((i_x.eq.1).and.(i_y.eq.1).and.(i_m.eq.1).and.(i_z.eq.1).and.(found)) call fill_map(i_nei,l-1,vec)

                        enddo
                     enddo
                  enddo
               enddo

            enddo
         enddo
      enddo
   enddo
   avant=avant+indexNN(i_Nei)
enddo

end subroutine mapping_3D
!
!
!input variables are
!net,tot_N_Nneigh,tabledist(:,:),N_Nneigh,Nei_il,Nei_z,phase,motif,indexNN(:,:),tableNN(:,:,:,:,:,:)

subroutine mapping_3D_motif_SL(d,nei,Nei_il,Nei_z,indexNN,tableNN,my_lattice)
use m_vector , only : norm
use m_derived_types, only : cell,lattice
#ifdef CPP_MPI
use m_make_box, only : Xstart,Ystart,Zstart
#endif
implicit none
! variable that come in
integer, intent(in) :: nei,Nei_il,Nei_z
type(lattice), intent(in) :: my_lattice
real(kind=8), intent(in) :: d(:,:)
integer, intent(in) :: indexNN(:,:)
! value of the function
integer :: tableNN(:,:,:,:,:,:)
! external blas
! 3D coordinate ix,iy,iz of 1d coordinate k
integer :: i_x,i_y,i_z,i_m,Xstop,Ystop,Zstop,Mstop
integer :: v_x,v_y,v_z,ok
! dummy variable
integer :: i,j,k,l,i_Nei,avant,i_p,i_phase,dim_lat(3)
real (kind=8) :: vec(3),dist,r(3,3)
logical :: found
#ifndef CPP_MPI
integer, parameter ::  Xstart=1
integer, parameter ::  Ystart=1
integer, parameter ::  Zstart=1
#endif

avant=0
Xstop=size(tableNN,3)
Ystop=size(tableNN,4)
Zstop=size(tableNN,5)
Mstop=size(tableNN,6)
do i=1,3
  r(:,i)=my_lattice%areal(i,:)
enddo
dim_lat=my_lattice%dim_lat

i_phase=1

do i_nei=1,nei
   do i_y=1,Ystop
      do i_x=1,Xstop
         do i_z=1,Zstop
            do i_m=1,Mstop

               l=1
               do i=i_nei,-i_nei,-1
                  do j=i_nei,-i_nei,-1
                     do i_p=1,Mstop

                        vec=0.0d0
                        found=.false.
                        ! suppose that the neighbour should be taken into account
                        ok=1

                        call test_neighbour(vec,v_x,ok,i_x,i,Xstop,r(:,1),my_lattice%boundary(1))
                        call test_neighbour(vec,v_y,ok,i_y,j,Ystop,r(:,2),my_lattice%boundary(2))
                        dist=norm(vec)

                        call associate_neighbour(tableNN(:,avant+l,i_x,i_y,i_z,i_m),v_x,v_y,i_z,i_p,ok,l,d(i_nei,i_phase),dist,found)

                        if ((i_x.eq.1).and.(i_y.eq.1).and.(i_m.eq.1).and.(i_z.eq.1).and.(found)) call fill_map(i_nei,l-1,vec)

                     enddo
                  enddo
               enddo

            enddo
         enddo
      enddo
   enddo
   avant=avant+indexNN(i_Nei,i_phase)
enddo

i_phase=2

do i_nei=1,Nei_il
   do i_y=1,Ystop
      do i_x=1,Xstop
         do i_z=1,Zstop
            do i_m=1,Mstop

               l=1
               do i=i_nei,-i_nei,-1
                  do j=i_nei,-i_nei,-1
                     do i_p=1,Mstop

                        vec=0.0d0
                        found=.false.
                        ! suppose that the neighbour should be taken into account
                        ok=1

                        call test_neighbour(vec,v_x,ok,i_x,i,Xstop,r(:,1),my_lattice%boundary(1))
                        call test_neighbour(vec,v_y,ok,i_y,j,Ystop,r(:,2),my_lattice%boundary(2))
                        dist=norm(vec)

                        call associate_neighbour(tableNN(:,avant+l,i_x,i_y,i_z,i_m),v_x,v_y,i_z,i_p,ok,l,d(i_nei,i_phase),dist,found)

                        if ((i_x.eq.1).and.(i_y.eq.1).and.(i_m.eq.1).and.(i_z.eq.1).and.(found)) call fill_map(i_nei,l-1,vec)

                     enddo
                  enddo
               enddo

            enddo
         enddo
      enddo
   enddo
   avant=avant+indexNN(i_Nei,i_phase)
enddo

i_phase=3

do i_nei=1,Nei_z
   do i_y=1,Ystop
      do i_x=1,Xstop
         do i_z=1,Zstop
            do i_m=1,Mstop

               l=1
               do i=i_nei,-i_nei,-1
                  do j=i_nei,-i_nei,-1
                     do k=i_nei,-i_nei,-1
                        do i_p=1,Mstop

                           vec=0.0d0
                           found=.false.
                           ! suppose that the neighbour should be taken into account
                           ok=1

                           call test_neighbour(vec,v_x,ok,i_x,i,Xstop,r(:,1),my_lattice%boundary(1))
                           call test_neighbour(vec,v_y,ok,i_y,j,Ystop,r(:,2),my_lattice%boundary(2))
                           call test_neighbour(vec,v_z,ok,i_z,k,Zstop,r(:,3),my_lattice%boundary(3))
                           dist=norm(vec)

                           call associate_neighbour(tableNN(:,avant+l,i_x,i_y,i_z,i_m),v_x,v_y,v_z,i_p,ok,l,d(i_nei,i_phase),dist,found)

                           if ((i_x.eq.1).and.(i_y.eq.1).and.(i_m.eq.1).and.(i_z.eq.1).and.(found)) call fill_map(i_nei,l-1,vec)

                        enddo
                     enddo
                  enddo
               enddo

            enddo
         enddo
      enddo
   enddo
   avant=avant+indexNN(i_Nei,i_phase)
enddo

end subroutine mapping_3D_motif_SL




! subroutine that associates the neighbour in the table of neighbours
subroutine associate_neighbour(tableNN,v_x,v_y,v_z,v_m,ok,l,d_ref,d_test,found)
implicit none
integer, intent(inout) :: tableNN(5),l
logical, intent(inout) :: found
real(kind=8), intent(in) :: d_ref,d_test
integer, intent(in) :: v_x,v_y,v_z,v_m,ok
! internal variables

if (dabs(d_ref-d_test).lt.1.0d-8) then
   tableNN(1)=v_x
   tableNN(2)=v_y
   tableNN(3)=v_z
   tableNN(4)=v_m
   tableNN(5)=ok
   l=l+1
   found=.true.
endif

end subroutine associate_neighbour


! subroutine that tests the neighbour
subroutine test_neighbour(vec,v,ok,i_x,i,stop,r,period)
implicit none
real(kind=8), intent(inout) :: vec(3)
integer, intent(out) :: v
integer, intent(inout) :: ok
real(kind=8), intent(in) :: r(3)
integer, intent(in) :: i_x,i,stop
logical, intent(in) :: period
! internal
integer :: test

test=i_x+i
! make a translation of r if the test>stop
vec=vec+real(i)*r
! use the periodic boundary conditions if necessary
v=periodic(test,stop)

if (((test.lt.1).or.(test.gt.stop)).and.(.not.period)) ok=0

end subroutine test_neighbour

! return the position of the neighbor in internal units with the periodic boundary conditions
integer function periodic(i,N)
implicit none
integer, intent(in) :: i,N

periodic=mod(i-1+N,N)+1

end function periodic

subroutine check_l(l,Nmax)
implicit none
integer, intent(in) :: l,Nmax

if (l.gt.Nmax) then
   write(6,'(a)') 'Error in the Mapping routine'
   stop
endif

end subroutine check_l





!
! Subroutine that writes out the table of neighbours for selected sites
!
!
subroutine write_out(tableNN,indexNN)
implicit none
integer, intent(in) :: tableNN(:,:,:,:,:,:),indexNN(:,:)
! internal
integer :: i_nei,i_shell
integer :: shape_tableNN(6)

shape_tableNN=shape(tableNN)

write(6,'(a)') ''
write(6,'(a)') 'Table of neighbors for selected sites'
write(6,'(a)') '(cut out at shell 2)'
write(6,'(a)') '-------------------------------------'
write(6,'(a)') '(1,1,1,1)'

do i_shell=1,size(indexNN,1)
  write(6,'(a,I4)') 'Shell  ', i_shell
  do i_nei=1,indexNN(i_shell,1)
    write(6,'(5(a,I8,2X))') 'i_x=',tableNN(1,i_nei,1,1,1,1),'i_y=',tableNN(2,i_nei,1,1,1,1), 'i_z=',tableNN(3,i_nei,1,1,1,1), 'i_m=',tableNN(4,i_nei,1,1,1,1), 'present', tableNN(5,i_nei,1,1,1,1)
  enddo
  if (i_shell.eq.3) exit
enddo

write(6,'(a)') '-------------------------------------'
write(6,'(a)') '(N/2+1,N/2+1,N/2+1,N/2+1)'
do i_shell=1,size(indexNN,1)
  write(6,'(a,I4)') 'Shell  ', i_shell
  do i_nei=1,indexNN(i_shell,1)
    write(6,'(5(a,I8,2X))') 'i_x=',tableNN(1,i_nei,shape_tableNN(3)/2+1,shape_tableNN(4)/2+1,shape_tableNN(5)/2+1,shape_tableNN(6)/2+1), &
        & 'i_y=',tableNN(2,i_nei,shape_tableNN(3)/2+1,shape_tableNN(4)/2+1,shape_tableNN(5)/2+1,shape_tableNN(6)/2+1), &
        & 'i_z=',tableNN(3,i_nei,shape_tableNN(3)/2+1,shape_tableNN(4)/2+1,shape_tableNN(5)/2+1,shape_tableNN(6)/2+1), &
        & 'i_m=',tableNN(4,i_nei,shape_tableNN(3)/2+1,shape_tableNN(4)/2+1,shape_tableNN(5)/2+1,shape_tableNN(6)/2+1), &
        & 'present', tableNN(5,i_nei,shape_tableNN(3)/2+1,shape_tableNN(4)/2+1,shape_tableNN(5)/2+1,shape_tableNN(6)/2+1)
  enddo
  if (i_shell.eq.3) exit
enddo

write(6,'(a)') '-------------------------------------'
write(6,'(a)') '(N/2+1,N/2+1,N/2+1,N/2+1)'
do i_shell=1,size(indexNN,1)
  write(6,'(a,I4)') 'Shell  ', i_shell
  do i_nei=1,indexNN(i_shell,1)
    write(6,'(5(a,I8,2X))') 'i_x=',tableNN(1,i_nei,shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6)), &
        & 'i_y=',tableNN(2,i_nei,shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6)), &
        & 'i_z=',tableNN(3,i_nei,shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6)), &
        & 'i_m=',tableNN(4,i_nei,shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6)), &
        & 'present', tableNN(5,i_nei,shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
  enddo
  if (i_shell.eq.3) exit
enddo

end subroutine

end module m_mapping
