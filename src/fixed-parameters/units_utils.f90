module m_units_utils
! all constants conversion are in constant.f90
use m_constants
use, intrinsic :: iso_fortran_env, only : output_unit
implicit none

public
contains

subroutine unit_convert(var_nam,name,in,out,alpha)
real(8), intent(in)           :: in
real(8), intent(out)          :: out
character(len=*), intent(in)  :: name
character(len=*), intent(in)  :: var_nam
real(8), intent(in), optional :: alpha

select case (name)
   case ('')
     out=in
   case ('Joule')
     out=in*J_to_eV
     write(output_unit,'(5(a,2X)/)') "converting",var_nam,"from",'Joule',"to eV"
   case ('meter')
     out=in*m_to_nm
     write(output_unit,'(5(a,2X)/)') "converting",var_nam,"from",'meter',"to nm"
   case ('Angstrom')
     out=in*A_to_nm
     write(output_unit,'(5(a,2X)/)') "converting",var_nam,"from",'Angstrom',"to nm"
   case ('Bohr')
     out=in*bohr_to_nm
     write(output_unit,'(5(a,2X)/)') "converting",var_nam,"from",'Bohr',"to nm"
   case ('kg')
     out=in*kg_to_amu
     write(output_unit,'(5(a,2X)/)') "converting",var_nam,"from",'kg',"to amu"
   case ('second')
     out=in*s_to_fs
     write(output_unit,'(5(a,2X)/)') "converting",var_nam,"from",'second',"to fs"
  case ('microC/cm2')
     out=in*1.0d-20*alpha
     write(output_unit,'(5(a,2X)/)') "converting",var_nam,"from",'microC/cm2',"to C.nm"
   case default
     write(output_unit,'(a)') "units not found. Choose between:"
     write(output_unit,'(a)') "Joule meter Angstrom Bohr kg second"
     stop
end select

end subroutine


end module m_units_utils
