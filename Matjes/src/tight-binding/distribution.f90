module m_distribution
implicit none
public

    abstract interface
     	function int_distrib(E_F, energy, kt)result(res)
         real(8)            :: res
         real(8),intent(in) :: E_F
         real(8),intent(in) :: kt
         real(8),intent(in) :: energy
       end function 
    end interface

contains 

     ! Function implementing the FD distribution for a given energy
     function fermi_distrib(E_F, energy, kt)result(res)
         real(8)                  :: res
         real(kind=8),intent(in) :: E_F
         real(kind=8),intent(in) :: kt
         real(kind=8),intent(in) :: energy

         real(8)                 ::  exp_val
         real(8)                 ::  cutoff=30.0d0

         exp_val=(energy-E_F)/kt
         if(exp_val>cutoff)then
             res=0.0d0
         elseif(exp_val<-cutoff)then
             res=1.0d0
         else
             res = 1.0d0/( 1.0d0 + exp(exp_val ) )
         endif
     end function fermi_distrib


    ! Function implementing the FD distribution derivative for a given energy
    function dE_fermi_distrib(E_F, energy, kt) result(res)
        real(8)            :: res
        real(8),intent(in) :: E_F
        real(8),intent(in) :: kt
        real(8),intent(in) :: energy

        real(8)                 ::  exp_val
        real(8)                 ::  cutoff=30.0d0

        exp_val=(energy-E_F)/kt
        if(abs(exp_val)>cutoff)then
            res=0.0d0
        else
            exp_val=exp(exp_val)
            res= -exp_val/(1.0d0 + exp_val)/(1.0d0 + exp_val)/kT
        endif
    end function


end module
