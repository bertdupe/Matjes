       module m_fit_old
       interface fit
        module procedure fit_1d
        module procedure fit_2d
        module procedure fit_Nd
       end interface fit
       interface get_a
        module procedure a_2d
        module procedure a_1d
       end interface get_a
       contains

       function fit_1d(Y,N)
       implicit none
       real(kind=8) :: Y(:)
       integer :: N
! internals
! y=ax+b
! a=fit1_2d(2) and b=fit1_2d(1)
       real(kind=8) :: fit_1d(2)
! internals
       real(kind=8) :: X(N)
       integer :: Nx,Ny,i
       real(kind=8) :: avex,avey,sumxy,sumx2,st,sr,syx,r2,sumx,sumy

       X=(/(i,i=1,N)/)

       Nx=size(X)
       Ny=size(Y)

       if (Nx.ne.Ny) stop 'problem in the data for the 2D fit'

       sumxy=0.0d0
       sumx2=0.0d0
       sumx=0.0d0
       sumy=0.0d0
       do i=1,Nx
        sumx=sumx+X(i)
        sumy=sumy+Y(i)
        sumx2=sumx2+X(i)*X(i)
        sumxy=sumxy+X(i)*Y(i)
       enddo
       avex=sumx/dble(Nx)
       avey=sumy/dble(Nx)

       fit_1d(2)=(dble(Nx)*sumxy-sumx*sumy)/(dble(Nx)*sumx2-sumx*sumx)
       fit_1d(1)=avey-fit_1d(2)*avex

       st=0.0d0
       sr=0.0d0
       do i=1,Nx
        st=st+(Y(i)-avey)**2
        sr=sr+(Y(i)-fit_1d(2)*X(i)-fit_1d(1))**2
       enddo

       syx=sqrt(sr/(dble(Nx-2)))
       r2=(st-sr)/st

       write(6,'(2(a,E20.12))') '1D (y=a*x+b) regression found with a=', fit_1d(2), '  and b=', fit_1d(1)
       write(6,'(a,E20.12,a)') 'coefficient of determination (r2)', r2, '  (perfect fit r2=1)'
       write(6,'(a,E20.12,a)') 'standard error', syx, '  (perfect fit S=0)'

       end function fit_1d

!---------------------------------------------------

       function fit_2d(X,Y)
       implicit none
       real(kind=8) :: X(:),Y(:)
! y=ax+b
! a=fit1_2d(2) and b=fit1_2d(1)
       real(kind=8) :: fit_2d(2)
! internals
       integer :: Nx,Ny,i
       real(kind=8) :: avex,avey,sumxy,sumx2,st,sr,syx,r2,sumx,sumy

       Nx=size(X)
       Ny=size(Y)

       if (Nx.ne.Ny) stop 'problem in the data for the 2D fit'

       sumxy=0.0d0
       sumx2=0.0d0
       sumx=0.0d0
       sumy=0.0d0
       do i=1,Nx
        sumx=sumx+X(i)
        sumy=sumy+Y(i)
        sumx2=sumx2+X(i)*X(i)
        sumxy=sumxy+X(i)*Y(i)
       enddo
       avex=sumx/dble(Nx)
       avey=sumy/dble(Nx)

       fit_2d(2)=(dble(Nx)*sumxy-sumx*sumy)/(dble(Nx)*sumx2-sumx*sumx)
       fit_2d(1)=avey-fit_2d(2)*avex

       st=0.0d0
       sr=0.0d0
       do i=1,Nx
        st=st+(Y(i)-avey)**2
        sr=sr+(Y(i)-fit_2d(2)*X(i)-fit_2d(1))**2
       enddo

       syx=sqrt(sr/(dble(Nx-2)))
       r2=(st-sr)/st

       write(6,'(2(a,E20.12))') '1D (y=a*x+b) regression found with a=', fit_2d(2), '  and b=', fit_2d(1)
       write(6,'(a,E20.12,a)') 'coefficient of determination (r2)', r2, '  (perfect fit r2=1)'
       write(6,'(a,E20.12,a)') 'standard error', syx, '  (perfect fit S=0)'

       end function fit_2d

!---------------------------------------------------

       function fit_nd(X,Y,N)
       implicit none
       real(kind=8) :: X(:),Y(:)
       integer :: N
! y=sum_i a_i x**(i-1)
! the coefficients are in the same order than as x**(i-1)
       real(kind=8) :: fit_nd(N)
! internals
       integer :: Nx,Ny,i,j,k
! we solve the system so that MAT(N,N)=B(N)
! at column=cte, for all lines a_i=cte
       real(kind=8) :: MAT(N,N),B(N)
       real(kind=8) :: st,sr,avey,avex,r2,syx

       fit_nd=0.0d0
       Nx=size(X)
       Ny=size(Y)


       if (Nx.ne.Ny) stop 'problem in the data for the ND fit'
       if (N.ge.Ny) stop 'more variables that data point'

       avex=sum(X)/dble(Nx)
       avey=sum(Y)/dble(Ny)

       MAT=0.0d0
       B=0.0d0
       do i=1,N   ! column
        do j=1,N    ! line
         do k=1,Nx
          MAT(i,j)=MAT(i,j)+X(k)**(i-1+j-1)
         enddo
        enddo
        do k=1,Nx
         B(i)=B(i)+Y(i)*X(i)**(i-1)
        enddo
       enddo

       st=0.0d0
       sr=0.0d0
       do i=1,Nx
        st=st+(Y(i)-avey)**2
        sr=sr+(Y(i)-sum((/(fit_nd(j)*X(i)**(j-1),j=1,N)/)))**2
       enddo

       syx=sqrt(sr/(dble(Nx-(N+1))))
       r2=(st-sr)/st

       write(6,'(2(a,E20.12))') '1D (y=a*x+b) regression found with a=', fit_nd(1), '  and b=', fit_nd(2)
       write(6,'(a,E20.12,a)') 'coefficient of determination (r2)', r2, '  (perfect fit r2=1)'
       write(6,'(a,E20.12,a)') 'standard error', syx, '  (perfect fit S=0)'

       end function fit_nd

!---------------------------------------------------

       real(kind=8) function a_2d(X,Y)
       implicit none
       real(kind=8) :: X(:),Y(:)
! y=ax+b
! internals
       integer :: Nx,Ny,i
       real(kind=8) :: avex,avey,sumxy,sumx2,sumx,sumy

       Nx=size(X)
       Ny=size(Y)

       if (Nx.ne.Ny) stop 'problem in the data for the 2D fit'

       sumxy=0.0d0
       sumx2=0.0d0
       sumx=0.0d0
       sumy=0.0d0
       do i=1,Nx
        sumx=sumx+X(i)
        sumy=sumy+Y(i)
        sumx2=sumx2+X(i)*X(i)
        sumxy=sumxy+X(i)*Y(i)
       enddo
       avex=sumx/dble(Nx)
       avey=sumy/dble(Nx)

       a_2d=(dble(Nx)*sumxy-sumx*sumy)/(dble(Nx)*sumx2-sumx*sumx)

       end function a_2d

!---------------------------------------------------

       real(kind=8) function a_1d(Y)
       implicit none
       real(kind=8) :: Y(:)
! y=ax+b
! internals
       integer :: Ny,i
       real(kind=8) :: sumxy,sumx2,sumx,sumy

       Ny=size(Y)

       sumxy=0.0d0
       sumx2=0.0d0
       sumx=0.0d0
       sumy=0.0d0
       do i=1,Ny
        sumx=sumx+dble(i)
        sumy=sumy+Y(i)
        sumx2=sumx2+dble(i**2)
        sumxy=sumxy+dble(i)*Y(i)
       enddo

       a_1d=(dble(Ny)*sumxy-sumx*sumy)/(dble(Ny)*sumx2-sumx**2)

       end function a_1d

       end module m_fit_old
