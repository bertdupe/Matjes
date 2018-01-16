! =====================================================================
! SUBROUTINE Creation
! Creates a skyrmion at point (x0,y0) with radius 6w. The way to
! construct the skyrmion is inspired from the paper of Kubetzka et al.,
! whose paper's reference is PRB67,020401(R) (published in 2003).
! Author : Charles Paillard
! Starting Writing date : May 28th, 2013
! Last modification : June 5th, 2013
! Take Periodic Boundary Conditions into account
! Version : 3.0
! =====================================================================

      SUBROUTINE Creation(x0,y0,R0,ToPinOrNotToPin,coeff,star,my_lattice)

      use m_derived_types
      use m_skyrmion
      use m_constants, only : pi
      use m_parameters, only : H_ext,EA
      use m_lattice, only : spin,masque
      use m_vector, only : norm,distance,phi
      use m_parameters, only : c_DM,DM
      Implicit None
      type(lattice), intent(in) :: my_lattice
! (x0, y0) are the coordinates of the center of the skyrmion that we
! want to create0
      real(kind=8), intent(in) :: x0, y0
! R0 is the Skyrmion radius
      real(kind=8),intent(in) :: R0
! coeff defines if it is akyrmion or an antiskyrmion
! coeff=1 -> antiskyrmion and coeff=0 -> skyrmion
      real(kind=8),intent(in) :: coeff
! star=3 -> starskyrmion and star=1 sinon
      real(kind=8),intent(in) :: star
! To know if we pin (i.e. put the mask to zero) a spin
      Logical, intent(in) ::ToPinOrNotToPin
! Intermediate variables
      real(kind=8) :: Theta
      real(kind=8) :: Psi
      Integer::i_x,i_y,i_z,i_m,dim_lat(3)
      real(kind=8) :: dist,net(3,3)
      real(kind=8) ::epsnull
      real(kind=8) :: diml1,diml2
! Parameter to create a skyrmion woth theta profile according the paper
! of Kubetzka PRB 67, 020401(R) (2003)
      real(kind=8) :: widt, cen, dummyvar
! Variables to adjust the variable widt and cen according to Lilley
! criterion for Bloch walls dimension
      real(kind=8) :: Inflexparam, ThetaInflex, DThetaInflex
! Variables to detect the closest point on the lattice to the Skyrmion
! Cemter, in order to apply the mask and force its magnetization to be
! down.
      real(kind=8) :: MinDist
      Integer :: MinIndex(2)
! define the right or left hand chirality
      real(kind=8) :: chirality


      dim_lat=my_lattice%dim_lat
      net=my_lattice%areal
      if (abs(DM(1,1)).lt.1.0d-8) then
       chirality=1.0d0
      else
       chirality=-sign(1.0d0,c_DM*DM(1,1))
      endif

      if (chirality.gt.0.0d0) then
       write(6,'(a)') "right hand chirality selected"
      else
       write(6,'(a)') "left hand chirality selected"
      endif

       epsnull = 1.0d-8
! Debug =========================================
!       WRITE(*,*) 'x0 = ', x0, ' y0 = ', y0
! ===============================================
       If (norm(H_ext).gt.epsnull) then
               
          dummyvar = norm(H_ext)
       else
          dummyvar = norm(EA)
       endif
!       widt = 2.0d0*R0/(pi*COSH(0.5d0*dummyvar))
        Inflexparam = 2.0d0 *ACOSH( 0.5d0*SQRT(2.0d0+SQRT(6.0d0+2.0d0*COSH( &
                      2.0d0*dummyvar ) ) ) )
        ThetaInflex = pi(1.0d0)-ASIN( TANH (Inflexparam-dummyvar))&
                                  -ASIN( TANH (Inflexparam+dummyvar))
        DThetaInflex =  -(1.0d0/COSH(Inflexparam-dummyvar) + &
                                1.0d0/COSH(Inflexparam+dummyvar))
        widt = 2.d0*R0/(Inflexparam-ThetaInflex/DThetaInflex)

       cen = 0.5d0*widt*dummyvar

       Theta = 0.0d0
       Psi   = 0.0d0
       diml1 = DBLE(Dim_lat(1))
       diml2 = DBLE(Dim_lat(2))
       
       minIndex = 0
       MinDist = 1.0d10
! Check that R0 is positive
       If (R0.gt. 0.0d0) then
                             
       do i_m=1,size(spin,5)
        Do i_z=1,dim_lat(3)
         Do i_y=1,dim_lat(2)
          Do i_x=1,dim_lat(1)

        dist=distance(Spin(1:2,i_x,i_y,i_z,i_m),x0,y0,dim_lat,net)

          If (dist.lt.(1.4d0*R0)) then

             Theta = pi(1.0d0)-Asin( tanh( 2*(dist-cen)/widt ))- &
                               Asin( tanh( 2*(dist+cen)/widt ))
!             Phi   = ACOS( xpoint/dist )

! We parameterize the magnetization vector in spherical coordinates
! although the spatial coordinate is parameterized by cylindrical
! coordinates (See the thesis of Leonov for further details, page 32)

! Be careful. The Psi coordinate, as stated in Leonov's thesis, depends
! on the crystal symmetries.
! to write a skyrmion, x->-x i.e. a symmetry around the y axis
             Psi = Phi(Spin(1:2,i_x,i_y,i_z,i_m),x0,y0,dim_lat,net)+pi(1.0d0)
           else
           cycle
          endif

             Spin(4,i_x,i_y,i_z,i_m) = chirality * Sin(Theta) * Cos( star*(Psi + coeff*pi(1.0d0)))
             Spin(5,i_x,i_y,i_z,i_m) = chirality * Sin(Theta) * Sin( star*Psi )
             Spin(6,i_x,i_y,i_z,i_m) = Cos(Theta)

             if (MinDist.gt.dist) then
                 minIndex = (/i_x,i_y/)
                 MinDist  = dist
             endif

          enddo
         enddo
        enddo
       EndDo

       Endif
       
       If (ToPinOrNotToPin) then 
       masque(1,minIndex(1),minIndex(2),1) = 0
       if (coeff.lt.0.0d0) then
         WRITE(*,*) 'The Center of an artificially created skyrmion is', &
                  ' pinned at index ', minIndex, ' with location ', &
                  ' x = ', Spin(1,minIndex(1),minIndex(2),1,1), &
                  ' y = ', Spin(2,minIndex(1),minIndex(2),1,1)
       else
         WRITE(*,*) 'The Center of an artificially created antiskyrmion is', &
                  ' pinned at index ', minIndex, ' with location ', &
                  ' x = ', Spin(1,minIndex(1),minIndex(2),1,1), &
                  ' y = ', Spin(2,minIndex(1),minIndex(2),1,1)
       endif
       endif
! Debug
!       WRITE(*,*) 'Masque : ', (masque(minIndex, i), i=1,6)

      END SUBROUTINE Creation

! =====================================================================
! SUBROUTINE usercreation
! Creates a user defined skyrmion at point (x0,y0) with radius 6w. The way to
! construct the skyrmion is inspired from the paper of Kubetzka et al.,
! whose paper's reference is PRB67,020401(R) (published in 2003).
! Author : Charles Paillard and Bertrand Dupe
! Starting Writing date : February 9th, 2015
! Last modification : February 12th, 2015
! Take Periodic Boundary Conditions into account
! Version : 3.0
! =====================================================================

      SUBROUTINE usercreation(x0,y0,R0,my_lattice)

      use m_derived_types
      use m_skyrmion, only : coeffx,coeffy,starx,stary
      use m_constants, only : pi
      use m_parameters, only : H_ext,EA
      use m_lattice, only : spin
      use m_vector, only : norm,distance,phi
      use m_parameters, only : c_DM,DM
      Implicit None
      type(lattice), intent(in) :: my_lattice
! (x0, y0) are the coordinates of the center of the skyrmion that we
! want to create0
      real(kind=8), intent(in) :: x0, y0
! R0 is the Skyrmion radius
      real(kind=8),intent(in) :: R0
! Intermediate variables
      real(kind=8) :: Theta
      real(kind=8) :: Psi
      Integer::i_x,i_y,i_z,i_m,dim_lat(3)
      real(kind=8) :: dist, net(3,3)
      real(kind=8) ::epsnull
      real(kind=8) :: diml1,diml2
! Parameter to create a skyrmion woth theta profile according the paper
! of Kubetzka PRB 67, 020401(R) (2003)
      real(kind=8) :: widt, cen, dummyvar
! Variables to adjust the variable widt and cen according to Lilley
! criterion for Bloch walls dimension
      real(kind=8) :: Inflexparam, ThetaInflex, DThetaInflex
! Variables to detect the closest point on the lattice to the Skyrmion
! Cemter, in order to apply the mask and force its magnetization to be
! down.
      real(kind=8) :: MinDist
      Integer :: MinIndex(2)
! define the right or left hand chirality
      real(kind=8) :: chirality

      write(6,'(a)') "User defined skyrmion selected"

      dim_lat=my_lattice%dim_lat
      net=my_lattice%areal

      if (abs(DM(1,1)).lt.1.0d-8) then
       chirality=1.0d0
      else
       chirality=-sign(1.0d0,c_DM*DM(1,1))
      endif

      if (chirality.gt.0.0d0) then
       write(6,'(a)') "right hand chirality selected"
      else
       write(6,'(a)') "left hand chirality selected"
      endif

       write(6,*) "v_phi(x,y)",starx,stary
       write(6,*) "xPi(x,y)",coeffx,coeffy

       epsnull = 1.0d-8
! Debug =========================================
!       WRITE(*,*) 'x0 = ', x0, ' y0 = ', y0
! ===============================================
       If (norm(H_ext).gt.epsnull) then

          dummyvar = norm(H_ext)
       else
          dummyvar = norm(EA)
       endif
!       widt = 2.0d0*R0/(pi*COSH(0.5d0*dummyvar))
        Inflexparam = 2.0d0 *ACOSH( 0.5d0*SQRT(2.0d0+SQRT(6.0d0+2.0d0*COSH( &
                      2.0d0*dummyvar ) ) ) )
        ThetaInflex = pi(1.0d0)-ASIN( TANH (Inflexparam-dummyvar))&
                                  -ASIN( TANH (Inflexparam+dummyvar))
        DThetaInflex =  -(1.0d0/COSH(Inflexparam-dummyvar) + &
                                1.0d0/COSH(Inflexparam+dummyvar))
        widt = 2.d0*R0/(Inflexparam-ThetaInflex/DThetaInflex)

       cen = 0.5d0*widt*dummyvar

       Theta = 0.0d0
       Psi   = 0.0d0
       diml1 = DBLE(Dim_lat(1))
       diml2 = DBLE(Dim_lat(2))

       minIndex = 0
       MinDist = 1.0d10
! Check that R0 is positive
       If (R0.gt. 0.0d0) then

       do i_m=1,size(spin,5)
        Do i_z=1,dim_lat(3)
         Do i_y=1,dim_lat(2)
          Do i_x=1,dim_lat(1)

        dist=distance(Spin(1:2,i_x,i_y,i_z,i_m),x0,y0,dim_lat,net)

          If (dist.lt.(1.4d0*R0)) then

             Theta = pi(1.0d0)-Asin( tanh( 2*(dist-cen)/widt ))- &
                               Asin( tanh( 2*(dist+cen)/widt ))
!             Phi   = ACOS( xpoint/dist )

! We parameterize the magnetization vector in spherical coordinates
! although the spatial coordinate is parameterized by cylindrical
! coordinates (See the thesis of Leonov for further details, page 32)

! Be careful. The Psi coordinate, as stated in Leonov's thesis, depends
! on the crystal symmetries.
! to write a skyrmion, x->-x i.e. a symmetry around the y axis
             Psi = Phi(Spin(1:2,i_x,i_y,i_z,i_m),x0,y0,dim_lat,net)+pi(1.0d0)
           else
           cycle
          endif

             Spin(4,i_x,i_y,i_z,i_m) = chirality * Sin(Theta) * Cos( starx*Psi + coeffx*pi(1.0d0))
             Spin(5,i_x,i_y,i_z,i_m) = chirality * Sin(Theta) * Sin( stary*Psi + coeffy*pi(1.0d0))
             Spin(6,i_x,i_y,i_z,i_m) = Cos(Theta)

             if (MinDist.gt.dist) then
                 minIndex = (/i_x,i_y/)
                 MinDist = dist
             endif

          enddo
         enddo
        enddo
       EndDo

       Endif

! Debug
!       WRITE(*,*) 'Masque : ', (masque(minIndex, i), i=1,6)

      END SUBROUTINE usercreation

! =====================================================================
! SUBROUTINE targetcreation as defined by Leonov et al EPJ web conf 2014
! Creates a target skyrmions at point (x0,y0) with radius 6w. The way to
! construct the skyrmion is inspired from the paper of Kubetzka et al.,
! whose paper's reference is PRB67,020401(R) (published in 2003).
! Author : Bertrand Dupe
! Starting Writing date : June 21th, 2017
! Last modification :
! Take Periodic Boundary Conditions into account
! Version : 3.0
! =====================================================================

      SUBROUTINE targetcreation(x0,y0,R0,ToPinOrNotToPin,coeff,star,my_lattice)

      use m_derived_types
      use m_skyrmion
      use m_constants, only : pi
      use m_parameters, only : H_ext,EA
      use m_lattice, only : spin,masque
      use m_vector, only : norm,distance,phi
      use m_parameters, only : c_DM,DM
      Implicit None
      type(lattice), intent(in) :: my_lattice
! (x0, y0) are the coordinates of the center of the skyrmion that we
! want to create0
      real(kind=8), intent(in) :: x0, y0
! R0 is the Skyrmion radius
      real(kind=8),intent(in) :: R0
! coeff defines if it is akyrmion or an antiskyrmion
! coeff=1 -> antiskyrmion and coeff=0 -> skyrmion
      real(kind=8),intent(in) :: coeff
! star=3 -> starskyrmion and star=1 sinon
      real(kind=8),intent(in) :: star
! To know if we pin (i.e. put the mask to zero) a spin
      Logical, intent(in) ::ToPinOrNotToPin
! Intermediate variables
      real(kind=8) :: Theta
      real(kind=8) :: Psi
      Integer::i_x,i_y,i_z,i_m,dim_lat(3)
      real(kind=8) :: dist, net(3,3)
      real(kind=8) ::epsnull
      real(kind=8) :: diml1,diml2
! Parameter to create a skyrmion woth theta profile according the paper
! of Kubetzka PRB 67, 020401(R) (2003)
      real(kind=8) :: widt, cen, dummyvar
! Variables to adjust the variable widt and cen according to Lilley
! criterion for Bloch walls dimension
      real(kind=8) :: Inflexparam, ThetaInflex, DThetaInflex
! Variables to detect the closest point on the lattice to the Skyrmion
! Cemter, in order to apply the mask and force its magnetization to be
! down.
      real(kind=8) :: MinDist
      Integer :: MinIndex(2)
! define the right or left hand chirality
      real(kind=8) :: chirality

      dim_lat=my_lattice%dim_lat
      net=my_lattice%areal

      if (abs(DM(1,1)).lt.1.0d-8) then
       chirality=1.0d0
      else
       chirality=-sign(1.0d0,c_DM*DM(1,1))
      endif

      if (chirality.gt.0.0d0) then
       write(6,'(a)') "right hand chirality selected"
      else
       write(6,'(a)') "left hand chirality selected"
      endif

       epsnull = 1.0d-8
! Debug =========================================
!       WRITE(*,*) 'x0 = ', x0, ' y0 = ', y0
! ===============================================
       If (norm(H_ext).gt.epsnull) then

          dummyvar = norm(H_ext)
       else
          dummyvar = norm(EA)
       endif
!       widt = 2.0d0*R0/(pi*COSH(0.5d0*dummyvar))
        Inflexparam = 2.0d0 *ACOSH( 0.5d0*SQRT(2.0d0+SQRT(6.0d0+2.0d0*COSH( &
                      2.0d0*dummyvar ) ) ) )
        ThetaInflex = pi(1.0d0)-ASIN( TANH (Inflexparam-dummyvar))&
                                  -ASIN( TANH (Inflexparam+dummyvar))
        DThetaInflex =  -(1.0d0/COSH(Inflexparam-dummyvar) + &
                                1.0d0/COSH(Inflexparam+dummyvar))
        widt = 2.d0*R0/(Inflexparam-ThetaInflex/DThetaInflex)

       cen = 0.5d0*widt*dummyvar

       Theta = 0.0d0
       Psi   = 0.0d0
       diml1 = DBLE(Dim_lat(1))
       diml2 = DBLE(Dim_lat(2))

       minIndex = 0
       MinDist = 1.0d10
! Check that R0 is positive
       If (R0.gt. 0.0d0) then

       do i_m=1,size(spin,5)
        Do i_z=1,dim_lat(3)
         Do i_y=1,dim_lat(2)
          Do i_x=1,dim_lat(1)

        dist=distance(Spin(1:2,i_x,i_y,i_z,i_m),x0,y0,dim_lat,net)

          If (dist.lt.(1.4d0*R0)) then

             Theta = pi(1.0d0)-Asin( tanh( 2*(dist-cen)/widt ))- &
                               Asin( tanh( 2*(dist+cen)/widt ))
!             Phi   = ACOS( xpoint/dist )

! We parameterize the magnetization vector in spherical coordinates
! although the spatial coordinate is parameterized by cylindrical
! coordinates (See the thesis of Leonov for further details, page 32)

! Be careful. The Psi coordinate, as stated in Leonov's thesis, depends
! on the crystal symmetries.
! to write a skyrmion, x->-x i.e. a symmetry around the y axis
             Psi = Phi(Spin(1:2,i_x,i_y,i_z,i_m),x0,y0,dim_lat,net)+pi(1.0d0)
           else
! We parameterize the magnetization vector of the boundary of the target skyrmion as teh sk case in spherical coordinates
! although the spatial coordinate is parameterized by cylindrical
! coordinates (See the thesis of Leonov for further details, page 32)

! Be careful. The Psi coordinate, as stated in Leonov's thesis, depends
! on the crystal symmetries.
! to write a skyrmion, x->-x i.e. a symmetry around the y axis
             Theta = pi(1.0d0)-Asin( tanh( 2*(-2.4d0*R0+(dist-cen))/widt ))- &
                               Asin( tanh( 2*(-2.4d0*R0+(dist-cen))/widt ))
             Psi = Phi(Spin(1:2,i_x,i_y,i_z,i_m),x0,y0,dim_lat,net)+pi(1.0d0)
          endif

             Spin(4,i_x,i_y,i_z,i_m) = chirality * Sin(Theta) * Cos( star*(Psi + coeff*pi(1.0d0)))
             Spin(5,i_x,i_y,i_z,i_m) = chirality * Sin(Theta) * Sin( star*Psi )
             Spin(6,i_x,i_y,i_z,i_m) = Cos(Theta)

             if (MinDist.gt.dist) then
                 minIndex = (/i_x,i_y/)
                 MinDist  = dist
             endif

          enddo
         enddo
        enddo
       EndDo

       Endif

       If (ToPinOrNotToPin) then
       masque(1,minIndex(1),minIndex(2),1) = 0
       if (coeff.lt.0.0d0) then
         WRITE(*,*) 'The Center of an artificially created skyrmion is', &
                  ' pinned at index ', minIndex, ' with location ', &
                  ' x = ', Spin(1,minIndex(1),minIndex(2),1,1), &
                  ' y = ', Spin(2,minIndex(1),minIndex(2),1,1)
       else
         WRITE(*,*) 'The Center of an artificially created antiskyrmion is', &
                  ' pinned at index ', minIndex, ' with location ', &
                  ' x = ', Spin(1,minIndex(1),minIndex(2),1,1), &
                  ' y = ', Spin(2,minIndex(1),minIndex(2),1,1)
       endif
       endif
! Debug
!       WRITE(*,*) 'Masque : ', (masque(minIndex, i), i=1,6)

      END SUBROUTINE targetcreation

!======================================================================
!Function that destroys the skyrmion located at site a (x,y) with a
!radius R_S. (x,y) and R_S are input parameters, given by the user.
!Author : Charles Paillard
!Starting writing date : 24.05.2013 (May 24th, 2013)
! Last modification : May 31st, 2013
! For periodic boundary conditions
!Version: 2.1
!======================================================================

      SUBROUTINE  Annihilation(x_SkyCenter,y_SkyCenter, R_S,my_lattice)

      use m_derived_types
      use m_lattice, only : spin
      use m_skyrmion
      Implicit None
      type(lattice), intent(in) :: my_lattice
      real(kind=8), intent(in) :: R_S, x_SkyCenter,y_SkyCenter
      real(kind=8) :: R_SSquare
      Integer::i_x,i_y,i_z,i_m,dim_lat(3)
      real(kind=8) :: dist, xpoint, ypoint, xpboc1, xpboc2, ypboc1, ypboc2
      real(kind=8) :: DBottom,DTop,DLeft,DRight,DC1,DC2,DC3,DC4
      real(kind=8) :: diml1, diml2
      
      dist  = 0.0d0
      dim_lat=my_lattice%dim_lat
      R_SSquare = R_S**2
      diml1 = DBLE(Dim_lat(1))
      diml2 = DBLE(Dim_lat(2))

! Check that R0 is positive      
      If (R_S.gt. 0.0d0) then
   
      do i_m=1,size(spin,5)
       do i_x=1,dim_lat(1)
        do i_y=1,dim_lat(2)
         do i_z=1,dim_lat(3)
      
       xpoint  = Spin(1,i_x,i_y,i_z,i_m)-x_SkyCenter
       ypoint  = Spin(2,i_x,i_y,i_z,i_m)-y_SkyCenter
       dist    = xpoint**2 + ypoint**2 
! To test for periodic boundary conditions
       xpboc1  = xpoint-diml1
       ypboc1  = ypoint-diml2
       xpboc2  = xpoint+diml1
       ypboc2  = ypoint+diml2
       DBottom = xpoint**2 +ypboc1**2
       DTop    = xpoint**2 +ypboc2**2
       DLeft   = xpboc1**2 +ypoint**2
       DRight  = xpboc2**2 +ypoint**2
       DC1     = xpboc1**2 +ypboc1**2
       DC2     = xpboc2**2 +ypboc2**2
       DC3     = xpboc1**2 +ypboc2**2
       DC4     = xpboc2**2 +ypboc1**2

         if ((dist.le.R_SSquare).or.(DBottom.le.R_SSquare)     &
              .or.(DTop.le.R_SSquare).or.(DLeft.le.R_SSquare)  &
              .or.(DRight.le.R_SSquare).or.(DC1.le.R_SSquare)  &
              .or.(DC2.le.R_SSquare).or.(DC3.le.R_SSquare).or. &
                  (DC4.le.R_SSquare)) then
         Spin(4,i_x,i_y,i_z,i_m) = 0.0d0
         Spin(5,i_x,i_y,i_z,i_m) = 0.0d0
         Spin(6,i_x,i_y,i_z,i_m) = 1.0d0
         endif

         enddo
        enddo
       enddo
      enddo
      endif
!      WRITE(*,*) 'OK Annihilation'

      END SUBROUTINE Annihilation

! =====================================================================      
! SUBROUTINE to Pin a spin in the lattice
! =====================================================================
      
      SUBROUTINE Pinned(X0,Y0,R0,PinOrNot,my_lattice)

      use m_derived_types
      use m_lattice, only : spin,masque
      use m_skyrmion
      implicit none
      type(lattice), intent(in) :: my_lattice
      real(kind=8), intent(in) ::X0,Y0,R0
      Logical, intent(in) ::PinOrNot
      Integer::i_x,i_y,i_z,i_m
      real(kind=8) :: dist, xpoint, ypoint, xpboc1, xpboc2, ypboc1, ypboc2
      real(kind=8) :: DBottom,DTop,DLeft,DRight,DC1,DC2,DC3,DC4
      real(kind=8) :: diml1, diml2
! Variables to detect the closest point on the lattice to the Skyrmion
! Cemter, in order to apply the mask and force its magnetization to be
! down.
      real(kind=8) :: distprevious, MinDist
      Integer :: MinIndex(2), dim_lat(3)

       dim_lat=my_lattice%dim_lat
       diml1 = DBLE(Dim_lat(1))
       diml2 = DBLE(Dim_lat(2))
       
       minIndex = 0
       distprevious = 1.0d150
! Check that R0 is positive
       If (R0.gt. 0.0d0) then
                             

       do i_m=1,size(spin,5)
        do i_x=1,dim_lat(1)
         do i_y=1,dim_lat(2)
          do i_z=1,dim_lat(3)
       
       xpoint = Spin(1,i_x,i_y,i_z,i_m)-x0
       ypoint = Spin(2,i_x,i_y,i_z,i_m)-y0
       dist = SQRT( xpoint**2 + ypoint**2 )
! Debug =======================================================
!      WRITE(*,*) '==========================================='
!       WRITE(*,*)  ' x = ', Spin(i,1), &
!                   ' y = ', Spin(i, 2),' dist = ', dist
! =============================================================
! To test for periodic boundary conditions
       xpboc1 = xpoint-diml1
       ypboc1 = ypoint-diml2
       xpboc2 = xpoint+diml1
       ypboc2 = ypoint+diml2
       DBottom = SQRT(xpoint**2 +ypboc1**2 )
       DTop    = SQRT(xpoint**2 +ypboc2**2 )
       DLeft   = SQRT(xpboc1**2 +ypoint**2 )
       DRight  = SQRT(xpboc2**2 +ypoint**2 )
       DC1     = SQRT(xpboc1**2 +ypboc1**2 )
       DC2     = SQRT(xpboc2**2 +ypboc2**2 )
       DC3     = SQRT(xpboc1**2 +ypboc2**2 )
       DC4     = SQRT(xpboc2**2 +ypboc1**2 )

             MinDist = MIN(dist, DBottom, Dtop, DLeft, DRight, &
                           DC1, DC2, DC3, DC4)
             if (distprevious.gt.mindist) then 
                 minIndex = (/i_x,i_y/)
                 distprevious = MinDist
             endif

          enddo
         enddo
        enddo
       enddo

       Endif
       
       If (PinOrNot) then 
!       masque(minIndex, :) = 0.0d0
       masque(1,minIndex(1),minIndex(2),1)=0
       WRITE(*,*) 'A point on the lattice is', &
                  ' pinned at index ', minIndex, ' with location ', &
                  ' x = ', Spin(1,minIndex(1),minIndex(2),1,1), &
                  ' y = ', Spin(2,minIndex(1),minIndex(2),1,1)
       endif
! Debug
!       WRITE(*,*) 'Masque : ', (masque(minIndex, i), i=1,6)

      END SUBROUTINE Pinned      



