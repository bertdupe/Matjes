module m_eigen_val_vec

interface Jacobi
  module procedure :: eigenvalvec_real,eigenvalvec_cmplx
end interface

private
public :: Jacobi

contains


subroutine eigenvalvec_real(a,x,abserr,n)
!===========================================================
! Evaluate eigenvalues and eigenvectors
! of a real symmetric matrix a(n,n): a*x = lambda*x
! method: Jacoby method for symmetric matrices
! Alex G. (December 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output ...
! a(i,i) - eigenvalues
! x(i,j) - eigenvectors
! comments ...
!===========================================================
implicit none
integer, intent(in) :: n
real(kind=8), intent(in) :: abserr
real(kind=8), intent(inout) :: a(:,:),x(:,:)
! internal
integer :: i, j, k
real(kind=8) :: b2, bar
real(kind=8) :: beta, coeff, c, s, cs, sc

! initialize x(i,j)=0, x(i,i)=1
! *** the array operation x=0.0 is specific for Fortran 90/95
x = 0.0
do i=1,n
  x(i,i) = 1.0
end do

! find the sum of all off-diagonal elements (squared)
b2 = 0.0
do i=1,n
  do j=1,n
    if (i.ne.j) b2 = b2 + a(i,j)**2
  end do
end do

if (b2 <= abserr) return

! average for off-diagonal elements /2
bar = 0.5*b2/float(n*n)

do while (b2.gt.abserr)
  do i=1,n-1
    do j=i+1,n
      if (a(j,i)**2 <= bar) cycle  ! do not touch small elements
      b2 = b2 - 2.0*a(j,i)**2
      bar = 0.5*b2/float(n*n)
! calculate coefficient c and s for Givens matrix
      beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
      coeff = 0.5*beta/sqrt(1.0+beta**2)
      s = sqrt(max(0.5+coeff,0.0))
      c = sqrt(max(0.5-coeff,0.0))
! recalculate rows i and j
      do k=1,n
        cs =  c*a(i,k)+s*a(j,k)
        sc = -s*a(i,k)+c*a(j,k)
        a(i,k) = cs
        a(j,k) = sc
      end do
! new matrix a_{k+1} from a_{k}, and eigenvectors
      do k=1,n
        cs =  c*a(k,i)+s*a(k,j)
        sc = -s*a(k,i)+c*a(k,j)
        a(k,i) = cs
        a(k,j) = sc
        cs =  c*x(k,i)+s*x(k,j)
        sc = -s*x(k,i)+c*x(k,j)
        x(k,i) = cs
        x(k,j) = sc
      end do
    end do
  end do
end do

end subroutine


! CEigensystem.F
! diagonalization of a complex n-by-n matrix using the Jacobi algorithm
! code adapted from the "Handbook" routines for complex A
! (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
! this file is part of the Diag library
! last modified 24 Aug 15 th

!***********************************************************************
!* CEigensystem diagonalizes a general complex n-by-n matrix.
!* Input:   n, A = n-by-n matrix
!* Output:  d = vector of eigenvalues,
!*      U = transformation matrix,
!*      sort= sort eigenvalues based on their real parts (0 or 1)
!* these fulfill
!*  d = U A U^-1,  A = U d U^-1,  U A = d U  (UCOLS=0),
!*  d = U^-1 A U,  A = U^-1 d U,  A U = U d  (UCOLS=1).

subroutine eigenvalvec_cmplx(EPS, n, A,ldA, d, U,ldU, sort)
implicit none
integer, intent(in) :: n, ldA, ldU, sort
real(kind=8), intent(in) :: EPS
complex(kind=16), intent(inout) :: A(:,:) ! size lda,N
complex(kind=16), intent(out) :: U(:,:), d(:)  ! size ldU,N
! internal
integer :: p, q, j
real(kind=8) :: red, thresh, norm, off
complex(kind=16) :: delta, t, s, invc, sx, sy, tx, ty
complex(kind=16) :: x, y
complex(kind=16) :: ev(2,n),UL(ldU,n)
integer :: sweep

do p = 1, n
   ev(1,p) = 0.0d0
   ev(2,p) = A(p,p)
   d(p) = ev(2,p)
enddo

do p = 1, n
   do q = 1, n
      U(q,p) = 0.0d0
   enddo
   U(p,p) = 1.0d0
enddo

UL=U

red = .01D0/n**4

do sweep = 1, 50
   off = 0.0d0
   do q = 2, n
      do p = 1, q - 1
         off = off + abs(A(p,q))**2 + abs(A(q,p))**2
      enddo
   enddo

   if( off .lt. EPS ) then
     exit
   else
     if (sweep.eq.50) then
        write(6,'(a)') "Bad convergence in CEigensystem"
        exit
     endif
   endif

   thresh = 0.0d0
   if( sweep .lt. 4 ) thresh = off*red

   do q = 2, n
      do p = 1, q - 1
         off = abs(A(p,q))**2 + abs(A(q,p))**2
         if( sweep .gt. 4 .and. off .lt. EPS*(abs(ev(2,p))**2 + abs(ev(2,q))**2) ) then
            A(p,q) = 0
            A(q,p) = 0
         else if( off .gt. thresh ) then
            delta = A(p,q)*A(q,p)
            x = .5D0*(ev(2,p) - ev(2,q))
            y = sqrt(x**2 + delta)
            t = x - y
            s = x + y
            if( abs(t)**2 .lt. abs(s)**2 ) t = s

            t = 1/t
            delta = delta*t
            ev(1,p) = ev(1,p) + delta
            ev(2,p) = d(p) + ev(1,p)
            ev(1,q) = ev(1,q) - delta
            ev(2,q) = d(q) + ev(1,q)

            invc = sqrt(delta*t + 1)
            s = t/invc
            t = t/(invc + 1)
            sx = s*A(p,q)
            ty = t*A(p,q)
            sy = s*A(q,p)
            tx = t*A(q,p)

            do j = 1, n
              x = A(j,p)
              y = A(j,q)
              A(j,p) = x + sy*(y - ty*x)
              A(j,q) = y - sx*(x + tx*y)
              x = A(p,j)
              y = A(q,j)
              A(p,j) = x + sx*(y - tx*x)
              A(q,j) = y - sy*(x + ty*y)
            enddo

            A(p,q) = 0.0d0
            A(q,p) = 0.0d0

            do j = 1, n
              x = UL(p,j)
              y = UL(q,j)

              UL(p,j) = x + sx*(y - tx*x)
              UL(q,j) = y - sy*(x + ty*y)

            enddo

          else
            cycle
         endif
      enddo
   enddo

   do p = 1, n
      ev(1,p) = 0
      d(p) = ev(2,p)
   enddo
enddo


! normalize the eigenvectors
do p = 1, n
   norm = 0.0d0
   do q = 1, n
      norm = norm + abs(UL(p,q))**2
   enddo
   norm = 1/sqrt(norm)
   do q = 1, n
      UL(p,q) = UL(p,q)*norm
   enddo
enddo

if( sort .eq. 0 ) return

! sort the eigenvalues by their real part
do p = 1, n - 1
   j = p
   t = d(p)
   do q = p + 1, n
      if( sort*(Real(t) - Real(d(q))) .gt. 0 ) then
         j = q
         t = d(q)
      endif
   enddo

   if( j .ne. p ) then
      d(j) = d(p)
      d(p) = t
      do q = 1, n
         x = UL(p,q)
         UL(p,q) = UL(j,q)
         UL(j,q) = x
      enddo
   endif
enddo

end


end module m_eigen_val_vec
