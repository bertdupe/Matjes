      module m_Hamiltonian
      contains

      subroutine setup_Hamilton()
      use m_parameters
      use m_derived_types
      implicit none
! internal
      real(kind=8) :: energy,J1
      integer :: ic,il,i
      real(kind=8), save, target :: J_loc(3,3)
      real(kind=8), save, target :: zero(3,3)
      type(vec), save, target :: S(3)
      type(H_loc) :: voisin(3,3)
      type(vec_point) :: S2(3)
      logical :: test


      J_loc=0.0d0
      zero=0.0d0
      do i=1,3
         S(i)%w=(/0.0d0,0.0d0,1.0d0/)
         S2(i)%w => S(i)%w
      enddo

      J1=2.0d0

      J_loc(1,1)=1.0d0
      J_loc(2,2)=1.0d0
      J_loc(3,3)=1.0d0

      do ic=1,3
         do il=1,3
            nullify(voisin(il,ic)%H_loc)
         enddo
      enddo

      test = associated(voisin(1,2)%H_loc)
      if (.not.test) voisin(1,2)%H_loc => J_loc

      test = associated(voisin(2,1)%H_loc)
      if (.not.test) voisin(2,1)%H_loc => J_loc

      test = associated(voisin(2,3)%H_loc)
      if (.not.test) voisin(2,3)%H_loc => J_loc

      test = associated(voisin(3,2)%H_loc)
      if (.not.test) voisin(3,2)%H_loc => J_loc


!  make a routine of matrix multiply

      energy=0.0d0
      do ic=1,3
         do il=1,3
            test = associated(voisin(il,ic)%H_loc)
            if (.not.test) cycle
            energy = energy + get_E_local(S2(il)%w,voisin(il,ic)%H_loc,S(ic)%w,J1)
         enddo
      enddo

      end subroutine setup_Hamilton

      real(kind=8) function get_E_local(U,M,V,Jij)
      implicit none
      real(kind=8), intent(in) :: Jij
      real(kind=8), intent(in) :: V(3),U(3),M(3,3)
! internal
      real(kind=8) :: temp
      integer :: i,j

      get_E_local=0.0d0

      do i=1,3
         temp=0.0d0
         if (U(i).eq.0.0d0) cycle
         do j=1,3
            temp=temp+M(j,i)*V(j)
         enddo
         get_E_local=get_E_local+U(i)*temp*Jij
      enddo

      end function get_E_local

      end module m_Hamiltonian
