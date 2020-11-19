subroutine pimc(state,i_biq,i_dm,i_four,i_dip,EA,h_ext, &
                    & spin,shape_spin,tableNN,shape_tableNN,masque,shape_masque,indexNN,shape_index,N_cell)
use mtprng
implicit none
logical, intent(in) :: i_biq,i_dm,i_four,i_dip
integer, intent(in) :: shape_index(2),shape_spin(5),shape_tableNN(6),shape_masque(4),N_cell
real(kind=8), intent(inout) :: spin(shape_spin(1),shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
type(mtprng_state), intent(inout) :: state
integer, intent(in) :: tableNN(shape_tableNN(1),shape_tableNN(2),shape_tableNN(3),shape_tableNN(4),shape_tableNN(5),shape_tableNN(6))
integer, intent(in) :: masque(shape_masque(1),shape_masque(2),shape_masque(3),shape_masque(4))
integer, intent(in) :: indexNN(shape_index(1),shape_index(2))
real(kind=8), intent(in) :: EA(3),h_ext(3)

write(6,'(a)') 'you are in the Path integral Monte Carlo routine'


end subroutine pimc
