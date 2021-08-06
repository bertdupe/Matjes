module m_util_init
implicit none
public

contains

subroutine get_pos(lat,dim_mode,ordname,pos,dim_state)
    use m_derived_types,only: lattice 
    type(lattice), intent(in)       :: lat      !entire lattice containing geometric information
    integer,intent(in)              :: dim_mode
    character(*),intent(in)         :: ordname  !name of the order parameter
    real(8),allocatable,intent(out) :: pos(:)
    integer,intent(out)             :: dim_state
    integer     :: nmag
    
    nmag=lat%nmag
    if(dim_mode==1.and.nmag==1)then
        !choose position of magnetic atoms for initialization
        Call lat%get_pos_center(pos)
        dim_state=1
    elseif(dim_mode==1.and.nmag>1)then
        !use center of cell as position
        Call lat%get_pos_mag(pos)
        dim_state=1
    elseif(dim_mode==nmag*3)then
        !choose position of magnetic atoms for initialization
        Call lat%get_pos_mag(pos)
        dim_state=3
    elseif(dim_mode==3)then
        !use center of cell as position
        Call lat%get_pos_center(pos)
        dim_state=3
        write(*,*) "CANNOT USE get_pos FOR",ordname,"as its dim_modes is most unexpected" !add something with positions of all atoms, and multiplied by 3? (risky mixup in natom=3*nmag
        return
    endif
end subroutine

subroutine get_pos_vec(lat,dim_mode,ordname,pos)
    use m_derived_types,only: lattice 
    type(lattice), intent(in)       :: lat      !entire lattice containing geometric information
    integer,intent(in)              :: dim_mode
    character(*),intent(in)         :: ordname  !name of the order parameter
    real(8),allocatable,intent(out) :: pos(:)
    integer     :: nmag
    
    nmag=lat%nmag
    if(dim_mode==nmag*3)then
        !choose position of magnetic atoms for initialization
        Call lat%get_pos_mag(pos)
    elseif(dim_mode==3)then
        !use center of cell as position
        Call lat%get_pos_center(pos)
    else
        write(*,*) "CANNOT USE GET_POS_VEC FOR",ordname,"as its dim_mode is neither 3 nor 3*nmag"
        return
    endif
end subroutine

subroutine get_skyrmion(pos_sky,R0,coeffx,coeffy,starx,stary,chirality,pos,state,lat)
    !subroutine to add a single skyrmion a lattice given by (pos,state,lat)
    use m_constants, only : pi
    use m_derived_types, only: lattice
    real(8),intent(in)          :: pos_sky(3)
    real(8),intent(in)          :: R0 ! R0 is the Skyrmion radius
    real(8),intent(in)          :: coeffx,coeffy,starx,stary ! parameter to define if it is a skyrmions, antiskyrmions, Q=...
    real(8),intent(in)          :: chirality  ! define the right or left hand chirality
    real(8),intent(in)          :: pos(:,:)   !position of all entries(3,(Nmag)*Ncell)
    real(8),intent(inout)       :: state(:,:) !states corresponding to the positions (3,(Nmag)*Ncell)
    type(lattice), intent(in)   :: lat      !entire lattice containing geometric information
    ! Internal variables
    real(8)     :: theta,psi
    real(8)     :: diff(3)
    integer     :: i
    real(8)     :: dist
    ! Parameter to create a skyrmion woth theta profile according the paper
    ! of Kubetzka PRB 67, 020401(R) (2003)
    real(8)     :: widt, cen
    ! Variables to adjust the variable widt and cen according to Lilley
    ! criterion for Bloch walls dimension
    real(8)     :: Inflexparam, ThetaInflex, DThetaInflex
   
    If (R0 <= 0.0d0) return
    if(any(shape(pos)/=shape(state))) ERROR STOP "incorrect input for get_skyrmion_new, programming error" 
    if(size(pos,1)/=3) ERROR STOP "position shape wrong, programming error"

    write(6,'(a)') "User defined skyrmion selected"
    if (chirality.gt.0.0d0) then
       write(6,'(a)') "right hand chirality selected"
    else
       write(6,'(a)') "left hand chirality selected"
    endif
    write(6,*) "v_phi(x,y)",starx,stary
    write(6,*) "xPi(x,y)",coeffx,coeffy
    
    inflexparam = 2.0d0 *ACOSH( 0.5d0*SQRT(2.0d0+SQRT(6.0d0+2.0d0*COSH( 2.0d0) ) ) )
    ThetaInflex = pi-ASIN( TANH (Inflexparam-1.0d0))-ASIN( TANH (Inflexparam+1.0d0))
    DThetaInflex =  -(1.0d0/COSH(Inflexparam-1.0d0) + 1.0d0/COSH(Inflexparam+1.0d0))
    widt = 2.d0*R0/(Inflexparam-ThetaInflex/DThetaInflex)
    cen = 0.5d0*widt

    do i=1,size(pos,2)
        diff=pos(:,i)-pos_sky
        Call lat%min_diffvec(diff)
        dist=norm2(diff(1:2))
        !terrible parameter clculation directly translate from old implementation
        If (dist>=(1.4d0*R0)) cycle
        Theta = pi-Asin(tanh( 2*(dist-cen)/widt ))- &
                   Asin(tanh( 2*(dist+cen)/widt ))
        if(dist>1.0d-8)then
            psi = acos(diff(1)/dist)+pi
            if(diff(2)/dist<0.0d0) psi=-psi+2.0d0*pi
        else
            psi=0.0d0
        endif

        state(1,i) = chirality * Sin(Theta) * Cos( starx*Psi + coeffx*pi)
        state(2,i) = chirality * Sin(Theta) * Sin( stary*Psi + coeffy*pi)
        state(3,i) = Cos(Theta)
    enddo
end subroutine 


end module

