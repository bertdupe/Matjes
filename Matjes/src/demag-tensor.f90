! function that computes the demag tenser for a given r vector
! the tensor has the form
! demag_tensor(1)=Nxx
! demag_tensor(2)=Nyy
! demag_tensor(3)=Nzz
! demag_tensor(4)=Nxy
! demag_tensor(5)=Nyz
! demag_tensor(6)=Nxz
! from Newell et al Journ. Geophys. Res.: Solid Earth 98, 9551

      module m_dten
      contains

      function demag_tensor(r)
      use m_constants, only : pi
      use m_vector, only : norm
      implicit none
      real(kind=8) :: r(3)
      real(kind=8), parameter :: alpha=6.74582d-7
! alpha=mu0*muB/a**3; a in nm. The result is an energy in eV
      real(kind=8), dimension(6) :: demag_tensor

      if (norm(r).eq.0) then
       demag_tensor=alpha/3.0d0
!       demag_tensor=0.0d0
       return
      endif

      demag_tensor(1)=r(1)**2+norm(r)**2/3.0d0
      demag_tensor(2)=r(2)**2+norm(r)**2/3.0d0
      demag_tensor(3)=r(3)**2+norm(r)**2/3.0d0
      demag_tensor(4)=r(1)*r(2)
      demag_tensor(5)=r(2)*r(3)
      demag_tensor(6)=r(1)*r(3)

      demag_tensor=3.0d0*alpha*demag_tensor/norm(r)**5/pi(4.0d0)

      end function  demag_tensor

      end module m_dten
