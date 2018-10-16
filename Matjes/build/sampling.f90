      module m_sampling

      contains
      function sphereft(old_S,coco)
      use m_vector, only : norm
      use mtprng
      implicit none
      real(kind=8), intent(in) :: coco
      real(kind=8) , dimension(3), intent(in) :: old_s
      real(kind=8) , dimension(3) :: sphereft
! internal variable
      type(mtprng_state) :: state
      real(kind=8) :: costheta, sintheta, Choice, phi, new_s(3),ss
      real(kind=8) :: costhetaini, sinthetaini, cosphiini, sinphiini

#ifdef CPP_MRG
        choice=mtprng_rand_real1(state)
        costheta=Choice+dcos(coco)*(1-Choice)
        sintheta=dsqrt(1-costheta**2)
        choice=mtprng_rand_real1(state)
        phi=2*dacos(-1.0d0)*Choice
#else

        CALL RANDOM_NUMBER(Choice)
        costheta=Choice+dcos(coco)*(1-Choice)
        sintheta=dsqrt(1-costheta**2)
        CALL RANDOM_NUMBER(Choice)
        phi=2*dacos(-1.0d0)*Choice
#endif

! Sx,Sy,Sz new spin in the referential Sini||z
!        S_store(1)=sintheta*cos(phi)
!        S_store(2)=sintheta*sin(phi)
!C        S_store(3)=costheta

!C store cos(theta) and sin(theta) for the old spin
        costhetaini=old_S(3)
        sinthetaini=dsqrt(old_s(1)**2+old_s(2)**2)
        if ((old_s(1)**2+old_s(2)**2)<0.00001) then
          cosphiini=1.0d0
          sinphiini=0.0d0
        else
          cosphiini=old_s(1)/(dsqrt(old_s(1)**2+old_s(2)**2))
          sinphiini=old_s(2)/(dsqrt(old_s(1)**2+old_s(2)**2))
        endif

!C Rotation around z-axis on angle fi0=Heff_xy^Oy and
!C rotation around y-axis on angle theta0=Heff^Oz
!C The rotation matrix:
!C                  |cos(f) -sin(f)*cos(t)  sin(f)*sin(t)|
!C R=(Rx(t)*Rz(f))= |sin(f)  cos(f)*cos(t) -cos(f)*sin(t)|
!C                  |  0         sin(t)         cos(t)   |
!C the rotations assumed to be conter-clockwise !

        new_s(2)=sintheta*dcos(phi)*cosphiini+ &
      (sintheta*dsin(phi)*costhetaini+ &
      costheta*sinthetaini)*sinphiini
        new_s(1)=-sintheta*dcos(phi)*sinphiini+ &
      (sintheta*dsin(phi)*costhetaini+ &
      costheta*sinthetaini)*cosphiini
        new_s(3)=-sintheta*dsin(phi)*sinthetaini+ &
      costheta*costhetaini

      ss=norm(new_s)
      sphereft=new_s/ss

      end function sphereft
!!!!
! equirepqrtition of the spins
      function equirep(old_S,coco)
      use mtprng
      implicit none
      double precision, intent(in) :: coco
      double precision , dimension(3), intent(in) :: old_s
      double precision , dimension(3) :: equirep
!internal variable
      type(mtprng_state) :: state
      double precision :: Skal_store, choice, new_s(3)
      integer :: i_store

   12 Continue
   11 Continue
      new_s=0.0d0
      Skal_store=0.0d0
!C       Generate random Spin
#ifdef CPP_MRG
      Do i_store=1,3
       choice=mtprng_rand_real1(state)
       new_s(i_store)=2.0d0*(Choice-0.5d0)
       Skal_store=Skal_store+new_s(i_store)**2
      End do
#else
      Do i_store=1,3
        CALL RANDOM_NUMBER(Choice)
        new_s(i_store)=2.0d0*(Choice-0.5d0)
        Skal_store=Skal_store+new_s(i_store)**2
      End do
#endif

      Do i_store=1,3
        new_s(i_store)=new_s(i_store)/(dsqrt(Skal_store))
      End do

      if (dsqrt(Skal_store).gt.1.0d0) &
#ifdef CPP_MPI
        new_s=old_s
#else
        goto 11
#endif

!C     Calculate, if the angle between
!C     the new and the old spin is accepted
      Skal_store=dacos(dot_product(old_s,new_s))

!C     if the spin is not turned by an allowed angle, try again
      if ((Skal_store.GT.coco).or.(Skal_store.LT.-0.0001)) &
#ifdef CPP_MPI
        new_s=old_s
#else
      goto 12
#endif
       equirep=new_s
      end function equirep

      end module
