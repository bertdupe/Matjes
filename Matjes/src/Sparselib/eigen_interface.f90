module m_eigen_interface
use,intrinsic :: iso_c_binding
public

interface
  subroutine eigen_set_all_E(size_1,size_2,all_E ) bind( c, name="eigen_set_all_E" )
    use, intrinsic :: iso_c_binding
    integer( kind = c_int ), value :: size_1,size_2
    real( kind = c_double )        :: all_E(*)
  end subroutine


  subroutine eigen_matmul_allE(size_1,vec_1,size_2,vec_2,res) bind( c, name="eigen_matmul_allE" )
    use, intrinsic :: iso_c_binding
    integer( kind = c_int ), value,intent(in)  :: size_1,size_2
    real( kind = c_double ), intent(in)        :: vec_1(*),vec_2(*)
    real( kind = c_double ), intent(out)       :: res
  end subroutine


  subroutine eigen_set_H(Nentry,dimH,ind1,ind2,H ) bind( c, name="eigen_set_H" )
    use, intrinsic :: iso_c_binding
    integer( kind = c_int ), value :: Nentry,dimH
    integer( kind = c_int )        :: ind1(*),ind2(*)
    real( kind = c_double )        :: H(*)
  end subroutine


  subroutine eigen_set_B(Nentry,dimH,ind1,ind2,B ) bind( c, name="eigen_set_B" )
    use, intrinsic :: iso_c_binding
    integer( kind = c_int ), value :: Nentry,dimH
    integer( kind = c_int )        :: ind1(*),ind2(*)
    real( kind = c_double )        :: B(*)
  end subroutine

  subroutine eigen_eval_H(dimH,vec,res ) bind( c, name="eigen_eval_H" )
    use, intrinsic :: iso_c_binding
    integer( kind = c_int ), intent(in),value  :: dimH
    real( kind = c_double ), intent(in)        :: vec(*)
    real( kind = c_double ), intent(out)       :: res
  end subroutine

  subroutine eigen_eval_B(dimH,vec,res ) bind( c, name="eigen_eval_B" )
    use, intrinsic :: iso_c_binding
    integer( kind = c_int ), intent(in),value  :: dimH
    real( kind = c_double ), intent(in)        :: vec(*)
    real( kind = c_double ), intent(out)       :: res(dimH)
  end subroutine


  subroutine eigen_set_H_e(Nentry,dimH,ind1,ind2,H ) bind( c, name="eigen_set_H_e" )
    use, intrinsic :: iso_c_binding
    integer( kind = c_int ), value          :: Nentry,dimH
    integer( kind = c_int )                 :: ind1(*),ind2(*)
    complex( kind = c_double_complex )      :: H(*)
  end subroutine

  subroutine eigen_set_H_e_jsd(Nentry,dimH,ind1,ind2,H ) bind( c, name="eigen_set_H_e_jsd" )
    use, intrinsic :: iso_c_binding
    integer( kind = c_int ), value          :: Nentry,dimH
    integer( kind = c_int )                 :: ind1(*),ind2(*)
    complex( kind = c_double_complex )      :: H(*)
  end subroutine
end interface

end module 
