module m_init_2Q
use m_derived_types
use, intrinsic :: iso_fortran_env, only : error_unit
implicit none

private
public :: init_2Q

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the starting configuration as a spin spiral
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine init_2Q(io,fname,lat,ordname,dim_mode,state,init_conf)
    use m_io_utils, only: get_parameter
    use m_util_init, only: get_pos_vec
    integer,intent(in)              :: io       !init-file io-unit
    character(*),intent(in)         :: fname    !init-file name 
    type(lattice), intent(in)       :: lat      !entire lattice containing geometric information
    character(*),intent(in)         :: ordname  !name of the order parameter
    integer,intent(in)              :: dim_mode !dimension of the order parameter in each cell
    real(8),pointer,intent(inout)   :: state(:) !pointer the the order parameter
    real(8),intent(in)              :: init_conf(:)

    real(8)         :: q1(3),q2(3),qnorm
    real(8)         :: qk(3),qm(3)
!    real(8)         :: qvec(3),Rq(3),Iq(3),norm,qnorm(3)
    real(8),allocatable,target :: pos(:)
!    real(8),allocatable ::  position(:)
    real(8),pointer :: pos_3(:,:),state_3(:,:)
    real(8)         :: theta
    real(8)         :: phi
    integer         :: Nsite
    integer         :: i,j,size_unit_cell
   
    q1=0.0d0
    q2=0.0d0
    qnorm=0.0
    size_unit_cell=size(init_conf)

    call get_parameter(io,fname,'Q1_'//ordname,3,q1)
    call get_parameter(io,fname,'Q2_'//ordname,3,q2)
    call get_parameter(io,fname,'Qnorm_'//ordname,qnorm)
    !if(norm2(q1)==0.0d0) ERROR STOP "Q1 has to differ from 0"
    !if(norm2(q2)==0.0d0) ERROR STOP "Q2 has to differ from 0"
    !write(*,*)'Warning: Q1 and Q2 should be normalized by the user for Qnorm to be their actual norm.'
    !if(qnorm>0.0d0)then
  
    !else
    !    q1(:)=0.0d0
    !    q2(:)=0.0d0
    !endif
    !write(*,*)'before a* matrix mult: Q1=',q1,' Q2=',q2
    !write(*,*)'lat%astar=',lat%astar
    q1=matmul(q1,lat%astar)
    q2=matmul(q2,lat%astar)
 	 !q1=q1*qnorm
    !q2=q2*qnorm


	
    qk=q1 !(q1+q2)*0.5d0
    qm=q2 !(q1-q2)*0.5d0
	write(*,*)'Q1=',q1,' Q2=',q2,' QM=',qm,' QK=', qk
    Call get_pos_vec(lat,dim_mode,ordname,pos)
    Nsite=size(pos)/3
    pos_3(1:3,1:Nsite)=>pos
    state_3(1:3,1:Nsite)=>state
    
	 write(*,*)'Warning: mx and my have been swapped in 2Q state formula. my=sin(phi)'
    do i=1,Nsite

        phi  =dot_product(qm,pos_3(:,i))
        theta=dot_product(qk,pos_3(:,i))

        state_3(1,i)=cos(phi)*sin(theta)
        state_3(2,i)=sin(phi)
        state_3(3,i)=cos(phi)*cos(theta)
                 ! write(*,*)'pos_3(:,i)=',pos_3(:,i)
            	!	write(*,*)'phi=',phi
            		!write(*,*)'theta=',theta
            		!write(*,*)'mx=', state_3(1,i)
            		!write(*,*)'my=', state_3(2,i)
            		!write(*,*)'mz=', state_3(3,i)
        j=mod(i-1,size_unit_cell)+1
        state_3(:,i)=state_3(:,i)*init_conf(j)/abs(init_conf(j))
       ! write(*,*)'then state_3(:,i)=',state_3(:,i)
    enddo

    nullify(pos_3,state_3)
    deallocate(pos)
end subroutine

end module 
