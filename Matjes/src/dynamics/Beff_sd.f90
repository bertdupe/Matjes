      module m_eval_Beff
      interface calculate_Beff
       module procedure normal
      end interface calculate_Beff
      contains
! subroutine that calculates the field
! dE/DM
!
!
!
!--------------------------------------------------------------
! for normal
subroutine normal(i_x,i_y,i_z,i_m,B, &
         &   spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN, &
         my_lattice,h_int,Hamiltonian)
use m_fieldeff
use m_derived_types
implicit none
! input variable
type(Coeff_Ham), intent(in) :: Hamiltonian
integer, intent(in) :: i_x,i_y,i_z,i_m
integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4)
integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
real(kind=8), intent(in) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5)),h_int(3)
type(lattice), intent(in) :: my_lattice
! output of the function
real(kind=8), intent(out) :: B(:)
logical :: i_dip
! internals

i_dip=.False.
if (masque(1,i_x,i_y,i_z).eq.0) then
   B=0.0d0
   return
endif

#ifdef CPP_SUM
B=BZ(i_x,i_y,i_z,i_m,spin,shape_spin,h_int)  &
 &   +efield(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,Hamiltonian%c_Ji,Hamiltonian%exchange)
 if (Hamiltonian%i_exch) B=B+Bexch(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,Hamiltonian%c_Ji)
 if (Hamiltonian%i_ani) B=B+Bani(i_x,i_y,i_z,i_m,spin,shape_spin,Hamiltonian%c_ani,Hamiltonian%ani)
 if (Hamiltonian%i_DM) B=B+Bdm(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,Hamiltonian%c_DM,Hamiltonian%DMI)
 if (Hamiltonian%i_four) B=B+Bfour(i_x,i_y,i_z,i_m,spin,shape_spin,masque,shape_masque,my_lattice,Hamiltonian%c_Ki,Hamiltonian%fours)
 if (Hamiltonian%i_biq) B=B+Bbiqd(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,Hamiltonian%c_JB,Hamiltonian%biq)
 if (i_dip) B=B+Bdip(i_x,i_y,i_z,i_m,spin,shape_spin)
#endif

#ifdef CPP_BLAS

#endif

!#ifdef CPP_DEBUG
      write(*,*) Bexch(i_x,i_y,i_z,i_m,spin,shape_spin,indexNN,shape_index,masque,shape_masque,tableNN,shape_tableNN,Hamiltonian%c_Ji)
      write(*,*) BZ(i_x,i_y,i_z,i_m,spin,shape_spin,h_int)
      write(*,*) Bani(i_x,i_y,i_z,i_m,spin,shape_spin,Hamiltonian%c_ani,Hamiltonian%ani)
      stop
!#endif
      end subroutine normal

      end module m_eval_Beff
