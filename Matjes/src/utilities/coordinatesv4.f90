!===================================================================
! Program calculating the spherical coordinates of a vector using an 
! arbitrary vector as the polar reference
! Author : Charles Paillard
! date : 06/06/2013
! Version :4.0
! Home-made Gram-Schmidt
! Preparing integration into MC_code
!===================================================================

      SUBROUTINE SphericalCoordinates(Spins,shape_spin,angle_sum,VectPolarRef)

! The file FSpinInp will be referred to with the Unit number 86 in
! Read/Write Statements
! The file FSpinOut has unit number 89
! The format in which those file should be read or write is 5f14.8. This
! format is referred to by the number 99999
      use m_constants, only : pi
      Implicit None
      integer, intent(in) :: shape_spin(5)
! Table containing the sum of the angles and spin components
      real(kind=8), intent(in) ::Spins(3,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
! axis of references
     real(kind=8), intent(in) :: VectPolarRef(:)
! Average of the polar angles
      real(kind=8), intent(inout) ::angle_sum(2,shape_spin(2),shape_spin(3),shape_spin(4),shape_spin(5))
! Basis of orthononormal vectors with VectPolarRef as the polar vector
! reference (i.e. defining the polar angle theta in spherical
! coordinates) 
      real(kind=8), dimension(3,3)::ortho_matrix
      ! Other variables
      Logical::iexists
      Integer::counter, column, i,j,k,i_m
      real(kind=8), dimension(3)::Vectx, Vecty, colx, coly
      real(kind=8), dimension(:,:),allocatable:: ThetaTab
      real(kind=8)::testnull, modulus, Alpha, Beta, Gamm, mod2, &
          mod3, var, testtheta
      real(kind=8):: varx, vary

      testnull= 1.0d-4

      Alpha   = VectPolarRef(1)
      Beta    = VectPolarRef(2)
      Gamm    = VectPolarRef(3)
      modulus = SQRT(Alpha**2+Beta**2+ Gamm**2)
      Alpha   = Alpha/modulus
      Beta    = Beta/modulus
      Gamm    = Gamm/modulus

! Calling the Gram Schmidt process


! check that the polar vector is not along x or y
      Vectx(1) =  1.0d0
      Vectx(2) =  0.0d0
      Vectx(3) =  0.0d0
      Vecty(1) =  0.0d0
      Vecty(2) =  1.0d0
      Vecty(3) =  0.0d0
      colx     = dot_product(Vectx,VectPolarRef)
      coly     = dot_product(Vecty,VectPolarRef)
       
      if ((dABS(colx(1)-Vectx(1)).LE.testnull).AND. &
          (dABS(colx(2)-Vectx(2)).LE.testnull).AND. &
          (dABS(colx(3)-Vectx(3)).LE.testnull) ) then

       do i_m=1,shape_spin(5)
        do k=1,shape_spin(4)
         do j=1,shape_spin(3)
          do i=1,shape_spin(2)

! Computing theta and phi
           angle_Sum(1,i,j,k,i_m) = angle_Sum(1,i,j,k,i_m)+ACOS(Spins(3,i,j,k,i_m))
           if (abs(Spins(3,i,j,k,i_m)).gt.0.9999d0) then
            angle_Sum(2,i,j,k,i_m) = angle_Sum(2,i,j,k,i_m)+0.0d0
            cycle
           endif

           if (abs(Spins(1,i,j,k,i_m)).gt.testnull) then
            angle_Sum(2,i,j,k,i_m) = angle_Sum(2,i,j,k,i_m)+dATAN(Spins(2,i,j,k,i_m)/ &
                               Spins(1,i,j,k,i_m))
           else
             angle_Sum(2,i,j,k,i_m) = angle_Sum(2,i,j,k,i_m)+dATAN(Spins(2,i,j,k,i_m)/ &
                               Spins(1,i,j,k,i_m))+pi(1.0d0)
           endif

           enddo
          enddo
         enddo
        enddo

      else

!	Now, Gram Schmidt to orthonormalize VectPolarRef,x and y.

! 	We have to fill the matrix that we want to orthogonalize and we
! 	normalize the z  vector to 1.

        do i=1, 3
         ortho_matrix(i,3) = VectPolarRef(i)/modulus
        enddo

! Computing the new y vector
      mod2 = dSQRT(1.0d0-Beta**2)
      ortho_matrix(1,2) = -Alpha*Beta/mod2
      ortho_matrix(2,2) = (1.d0 - Beta**2)/mod2
      ortho_matrix(3,2) = -Gamm*Beta/mod2
      var               = ortho_matrix(1,2) 
      ortho_matrix(1,1) = 1.0d0 - var**2-Alpha**2
      ortho_matrix(2,1) = -var*ortho_matrix(2,2) - Alpha*Beta
      ortho_matrix(3,1) = -var* ortho_matrix(2,3)- Alpha*Gamm
! Normalizing the x vector to 1.
      mod3 = dSQRT(ortho_matrix(1,1)**2+ortho_matrix(1,2)**2+ &
                  ortho_matrix(1,3)**2   )
      ortho_matrix(1,1) = ortho_matrix(1,1)/mod3
      ortho_matrix(2,1) = ortho_matrix(2,1)/mod3
      ortho_matrix(3,1) = ortho_matrix(3,1)/mod3
      
        do i_m=1,shape_spin(5)
         do k=1,shape_spin(4)
          do j=1,shape_spin(3)
           do i=1,shape_spin(2)

         testtheta= ACOS(Spins(1,i,j,k,i_m)*ortho_matrix(1,3)+&
                              Spins(2,i,j,k,i_m)*ortho_matrix(2,3)+ &
                              Spins(3,i,j,k,i_m)*ortho_matrix(3,3))

         angle_Sum(1,i,j,k,i_m) =angle_Sum(1,i,j,k,i_m)+ testtheta
! Now we want to have phi between 0 and 2*pi to have usual spherical
! coordinates
! If the x component is negative, phi = pi - ATAN
! If the x component is >0 and the y component is <0, then phi =
! 2*pi+ATAN                
                
           varx = Spins(1,i,j,k,i_m)*ortho_matrix(1,1)+ &
                  Spins(2,i,j,k,i_m)*ortho_matrix(2,1)+ &
                  Spins(3,i,j,k,i_m)*ortho_matrix(3,1)
           vary = Spins(1,i,j,k,i_m)*ortho_matrix(1,2)+ &
                  Spins(2,i,j,k,i_m)*ortho_matrix(2,2)+ &
                  Spins(3,i,j,k,i_m)*ortho_matrix(3,2)

           if (testtheta.lt.testnull) then
            angle_Sum(2,i,j,k,i_m)=angle_Sum(2,i,j,k,i_m)+0.0d0
            cycle
           endif

           if (varx.gt.0.0d0) then
            angle_Sum(2,i,j,k,i_m)=angle_Sum(2,i,j,k,i_m)+ATAN(vary/varx)
           else
            angle_Sum(2,i,j,k,i_m)=angle_Sum(2,i,j,k,i_m)+ATAN(vary/varx)+pi(1.0d0)
           endif

!            if ((dSIGN(1.0d0,vary).LE.testnull).AND. &
!             (dSIGN(1.0d0,-varx).LE.testnull)) then
!              angle_Sum(2,i,j,k,i_m)=angle_Sum(2,i,j,k,i_m)+ pi(2.0d0)+ATAN(vary/varx)
!             elseif ( dSIGN(1.0d0,varx).LE.testnull) then
!              angle_Sum(2,i,j,k,i_m)=angle_Sum(2,i,j,k,i_m)+pi(1.0d0)-ATAN(vary/varx)
!             else
!              angle_Sum(2,i,j,k,i_m)=angle_Sum(2,i,j,k,i_m)+ATAN(vary/varx)
!            endif
!           endif
             

            enddo
           enddo
          enddo
         enddo

      endif

!      PRINT *, 'SUBROUTINE PolarCoordinatesSpin was run successfully '
!      PRINT *, 'You will find your results in file ', FSpinOut
      END SUBROUTINE SphericalCoordinates
