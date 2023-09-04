module m_lattice_position
!module which contains several routines to get or print position information for the lattice type
use m_type_lattice, only : lattice
use m_convert
private
public get_pos_mag, get_pos_center, get_pos_ph, get_pos_at_Id
public print_pos_ind

contains

subroutine print_pos_ind(lat,ind,fname)
    !print position of atoms referenced by ind (atom number in unit-cell) to file (fname)
    type(lattice),intent(in)    :: lat
    integer,intent(in)          :: ind(:)
    character(len=*),intent(in) :: fname
    integer             ::  io,N_int
    real(8),allocatable ::  pos(:)
    character(len=30)  :: format_out

    open(newunit=io,file=fname)
    Call get_pos_ind(lat,ind,pos)
    N_int=size(pos)/lat%Ncell/3
    format_out=convert('(',3*N_int,'E16.8,x)')   ! to use of you want to have all positions of the sites one fu on a line
!    format_out=convert('(',3,'E16.8,x)')
    write(io,format_out) pos
    close(io)
    deallocate(pos)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!   Get positions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_pos_atom(this,pos_out)
    !get the real-space position of all atom in the supercell
    class(lattice),intent(in)               :: this
    real(8),allocatable,intent(inout)       :: pos_out(:)

    integer     :: ind(size(this%cell%atomic))
    integer     :: i

    ind=[(i,i=1,size(this%cell%atomic))]
    Call get_pos_ind(this,ind,pos_out)
end subroutine

subroutine get_pos_at_Id(this,Id,pos_out)
    !get the real-space position of each magnetic atom in the supercell
    class(lattice),intent(in)               :: this
    integer       ,intent(in)               :: Id
    real(8),allocatable,intent(inout)       :: pos_out(:)

    integer, allocatable            :: ind(:)

    Call this%cell%ind_attype(Id,ind)
    Call get_pos_ind(this,ind,pos_out)
end subroutine

subroutine get_pos_mag(this,pos_out)
    !get the real-space position of each magnetic atom in the supercell
    class(lattice),intent(in)               :: this
    real(8),allocatable,intent(inout)       :: pos_out(:)

    integer, allocatable            :: ind(:)

    Call this%cell%ind_mag_all(ind)
    Call get_pos_ind(this,ind,pos_out)
end subroutine

subroutine get_pos_ph(this,pos_out)
    !get the real-space position of each magnetic atom in the supercell
    class(lattice),intent(in)               :: this
    real(8),allocatable,intent(inout)       :: pos_out(:)

    integer, allocatable            :: ind(:)

    Call this%cell%ind_M_all(ind)
    Call get_pos_ind(this,ind,pos_out)

end subroutine

subroutine get_pos_center(this,pos_out)
    !returns the position of the center of each unit-cell
    class(lattice),intent(in)               :: this
    real(8),allocatable,intent(out),target  :: pos_out(:)

    real(8),pointer,contiguous      :: pos(:,:,:,:)
    real(8),pointer,contiguous      :: pos3(:,:)
    real(8),allocatable             :: pos_lat(:,:,:,:)
    real(8)                         :: pos_center(3) 
    integer                         :: i

    call get_pos_cell(pos_lat,this%dim_lat,this%areal)
    allocate(pos_out(3*product(this%dim_lat)),source=0.0d0)
    pos(1:3,1:this%dim_lat(1),1:this%dim_lat(2),1:this%dim_lat(3))=>pos_out
    pos3(1:3,1:product(this%dim_lat))=>pos_out
    pos_center=sum(this%areal,1)*0.5d0
    pos=pos_lat
    do i=1,size(pos3,2)
        pos3(:,i)=pos3(:,i)+pos_center
    enddo
    deallocate(pos_lat)
    nullify(pos,pos3)
end subroutine

subroutine get_pos_ind(this,ind,pos_out)
    !returns the position of the atomic positions within in super-cell for all 
    !atoms whose indices given by ind (number in the unit-cell)
    use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit
    class(lattice),intent(in)               :: this
    integer,intent(in)                      :: ind(:)
    real(8),allocatable,intent(inout)       :: pos_out(:)

    real(8), allocatable            :: pos_lat(:,:,:,:)
    real(8), allocatable            :: pos_at(:,:)
    integer, allocatable            :: ind_at_mag(:)
    integer     :: shape_unfold(4)
    integer     :: Nind,Nat
    integer     :: i,i1,i2,i3,j

    if(allocated(pos_out)) deallocate(pos_out)
    Nind=size(ind)
    Nat=size(this%cell%atomic)
    if(any(ind<1).or.any(ind>Nat))then
        write(error_unit,'(//A,I5,A)') "Input indices for get_pos_ind are not within [1,",Nat,"]"
        write(error_unit,'(A)') "Indicies:"
        write(error_unit,'(3X,5I8)') ind
        ERROR STOP 
    endif
    shape_unfold(1)=Nind
    shape_unfold(2:4)=this%dim_lat
    call get_pos_cell(pos_lat,this%dim_lat,this%areal)
    allocate(pos_out(3*product(shape_unfold)),source=0.0d0)

    allocate(pos_at(3,Nind),source=0.0d0)
    do i=1,Nind
        pos_at(:,i)=this%cell%atomic(ind(i))%position
    enddo

    Call unfold_position(shape_unfold,pos_at,pos_lat,pos_out)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!   Utility functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine unfold_position(pos_shape,pos_at,pos_lat,pos_out)
    integer,intent(in)      ::  pos_shape(4)                                                        !dimensions [number atoms, dim_lat(1:3)]
    real(8),intent(in)      ::  pos_at(3,pos_shape(1))                                              !positions of considered atoms
    real(8),intent(in)      ::  pos_lat(3,pos_shape(2),pos_shape(3),pos_shape(4))                   !positions of considered lattice sites
    real(8),intent(inout)   ::  pos_out(3,pos_shape(1),pos_shape(2),pos_shape(3),pos_shape(4))      !output unfolded positions

    pos_out=spread(pos_lat,2,pos_shape(1))  !repeat the full lattice for each atomic entry (pos_shape(1))
    Call add_util(pos_shape(1)*3,pos_at,product(pos_shape(2:4)),pos_out)    !add the pos_at for each entry in pos_out(:.:,x,x,x)
end subroutine

subroutine add_util(N,arr_add,Nrep,arr_out)
    integer,intent(in)      :: N,Nrep
    real(8),intent(in)      :: arr_add(N)
    real(8),intent(inout)   :: arr_out(N,Nrep)
    integer ::  i

    do i=1,Nrep
        arr_out(:,i)=arr_out(:,i)+arr_add
    enddo
end subroutine

subroutine get_pos_cell(pos,dim_lat,lat)
    !get the positions of the 0 point of each cell within the super-cell
    real(8), intent(inout),allocatable  :: pos(:,:,:,:)
    integer, intent(in)                 :: dim_lat(:)
    real(8), intent(in)                 :: lat(:,:)
    ! internal variables
    real(8)     :: tmp(3,3)
    integer     :: i_z,i_y,i_x
   
    if(allocated(pos)) deallocate(pos)
    allocate(pos(3,dim_lat(1),dim_lat(2),dim_lat(3)),source=0.0d0)
    do i_z=1,dim_lat(3)
        tmp(:,3)=lat(3,:)*real(i_z-1,8)
        do i_y=1,dim_lat(2)
            tmp(:,2)=lat(2,:)*real(i_y-1,8)
            do i_x=1,dim_lat(1)
                tmp(:,1)=lat(1,:)*real(i_x-1,8)
                pos(1:3,i_x,i_y,i_z)=sum(tmp,dim=2)
            enddo
        enddo
    enddo
end subroutine 


end module
