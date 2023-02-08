MODULE m_random_number_library
! Inspired from: http://www.johndcook.com/julia_rng.html
! Original author in julia : John D Cook
! coded : Sukhbinder in fortran
! Date : 28th Feb 2012
!
! in order of appearance
! Uniform random Number Generators in Fortran
! Random Sample from normal (Gaussian) distribution
! Random smaple from an exponential distribution
! Return a random sample from a gamma distribution
! return a random sample from an inverse gamma random variable
! return a sample from a Weibull distribution
! return a random sample from a Cauchy distribution
! return a random sample from a Student t distribution
! return a random sample from a Laplace distribution
! return a random sample from a log-normal distribution
! return a random sample from a beta distribution
!
!

use m_constants, only : pi
use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit, output_unit
use m_random_base
implicit none



abstract interface

    function int_rnd_distrib(this) result(c)
       import ranbase
       class(ranbase)   :: this
       real(kind=8)     :: c
    end function

end interface

private
public :: select_rnd_in_lib,int_rnd_distrib

CONTAINS


subroutine select_rnd_in_lib(name,rnd_distrib)
character(len=3),intent(in)   :: name
procedure(int_rnd_distrib),pointer,intent(out)  :: rnd_distrib

select case (name)
   case('uni')
      write(output_unit,'(/a)') "uniform distribution selected"
      write(output_unit,'(a/)') "reading parameters min_rnd_val+x*(max_rnd_val-min_rnd_val)"
      rnd_distrib => rand_uniform

   case('nor')
      write(output_unit,'(/a/)') "normal (Gaussian) distribution"
      rnd_distrib => rand_normal

   case('exp')
      write(output_unit,'(/a/)') "exponential distribution selected"
      rnd_distrib => rand_exponential

   case('gam')
      write(output_unit,'(/a)') "gamma distribution selected"
      write(output_unit,'(a)') "Implementation based on -A Simple Method for Generating Gamma Variables-"
      write(output_unit,'(a)') "Vol 26, No 3, September 2000, pages 363-372"
      write(output_unit,'(a/)') "sigma=>scale and mean=>shape"
      rnd_distrib => rand_gamma_int

   case('chi')
      write(output_unit,'(/a)') "chi square distribution selected"
      write(output_unit,'(a/)') "mean=>dof"
      rnd_distrib => rand_chi_square

   case('iga')
      write(output_unit,'(/a)') "inverse gamma distribution selected"
      write(output_unit,'(a/)') "sigma=>scale and mean=>shape"
      rnd_distrib => rand_inverse_gamma

   case('wei')
      write(output_unit,'(/a)') "Weibull distribution selected"
      write(output_unit,'(a/)') "sigma=>scale and mean=>shape"
      rnd_distrib => rand_weibull

   case('cau')
      write(output_unit,'(/a)') "Cauchy distribution selected"
      write(output_unit,'(a/)') "sigma=>scale and mean=>median"
      rnd_distrib => rand_cauchy

   case('stu')
      write(output_unit,'(/a)') "Student t distribution selected"
      write(output_unit,'(a/)') "mean=>dof"
      rnd_distrib => rand_cauchy

   case('lap')
      write(output_unit,'(/a)') "laplace distribution selected"
      write(output_unit,'(a/)') "sigma=>scale and mean=>mean"
      rnd_distrib => rand_laplace

   case('log')
      write(output_unit,'(/a)') "log-normal distribution selected"
      write(output_unit,'(a/)') "sigma=>scale and mean=>mu"
      rnd_distrib => rand_log_normal

   case('bet')
      write(output_unit,'(/a)') "Beta distribution selected"
      write(output_unit,'(a/)') "sigma=>b and mean=>a"
      rnd_distrib => rand_beta

       case default
         stop 'please select a distribution for your random numbers'
    end select

end subroutine










!
! uniform random Number Generators in Fortran
!
FUNCTION rand_uniform(this) RESULT(c)
class(ranbase) :: this
real(kind=8) :: a,b,c,temp

   a=this%min_rnd_val
   b=this%max_rnd_val
   CALL this%rand_get(temp)

   c= a+temp*(b-a)
END FUNCTION

!
! Random Sample from normal (Gaussian) distribution
!
FUNCTION rand_normal(this) RESULT(c)
class(ranbase) :: this
real(kind=8) :: mean,stdev,c
mean=this%mean
stdev=this%sigma

c=rand_normal_int(mean,stdev,this)

END FUNCTION

FUNCTION rand_normal_int(mean,stdev,this) RESULT(c)
class(ranbase) :: this
real(kind=8) :: mean,stdev,c,temp(2),r,theta
c=0.0d0

IF(stdev <= 0.0d0) THEN

  WRITE(output_unit,'(a)') "Standard Deviation must be +ve"
ELSE
  CALL this%rand_get(temp(1))
  CALL this%rand_get(temp(2))
  r=(-2.0d0*log(temp(1)))**0.5
  theta = 2.0d0*PI*temp(2)
  c= mean+stdev*r*sin(theta)
END IF
END FUNCTION





!
! Random sample from an exponential distribution
!
FUNCTION rand_exponential(this) RESULT(c)
class(ranbase) :: this
real(kind=8) :: mean,c,temp
mean=this%mean
c=0.0d0
IF (mean <= 0.0d0) THEN

   WRITE(output_unit,'(a)') "mean must be positive"
ELSE
   CALL this%rand_get(temp)
   c=-mean*log(temp)
END IF
END FUNCTION

!
! Return a random sample from a gamma distribution
!
FUNCTION rand_gamma_int(this) RESULT(c)
class(ranbase) :: this
real(kind=8) :: SHAPE,scale,c
SHAPE=this%mean
scale=this%sigma
c=0.0d0

c=rand_gamma(shape, SCALE,this)

END FUNCTION

RECURSIVE FUNCTION rand_gamma(shape, SCALE,this) RESULT(ans)
class(ranbase) :: this
real(kind=8) SHAPE,scale,u,w,d,c,x,xsq,g,ans,v
IF (shape <= 0.0d0) THEN

   WRITE(output_unit,'(a)') "Shape PARAMETER must be positive"
END IF
IF (scale <= 0.0d0) THEN

   WRITE(output_unit,'(a)') "Scale PARAMETER must be positive"
END IF
!
! ## Implementation based on "A Simple Method for Generating Gamma Variables"
! ## by George Marsaglia and Wai Wan Tsang.
! ## ACM Transactions on Mathematical Software

! ## Vol 26, No 3, September 2000, pages 363-372.
!
IF (shape >= 1.0d0) THEN
   d = SHAPE - 1.0d0/3.0d0
   c = 1.0d0/(9.0d0*d)**0.5
   DO while (.true.)
      x = rand_normal_int(0.0d0, 1.0d0, this)
      v = 1.0 + c*x
      DO while (v <= 0.0d0)
         x = rand_normal_int(0.0d0, 1.0d0,this)
         v = 1.0d0 + c*x
      END DO

      v = v*v*v
      CALL this%rand_get(u)
      xsq = x*x
      IF ((u < 1.0d0 -.0331d0*xsq*xsq) .OR.  &
         (log(u) < 0.5d0*xsq + d*(1.0d0 - v + log(v))) )then
         ans=scale*d*v
         RETURN
      END IF

   END DO
ELSE
   g = rand_gamma(shape+1.0d0, 1.0d0,this)
   CALL this%rand_get(w)
   ans=scale*g*(w)**(1.0d0/shape)
   RETURN
END IF

END FUNCTION
!
! ## return a random sample from a chi square distribution
! ## with the specified degrees of freedom
!
FUNCTION rand_chi_square(this) RESULT(ans)
class(ranbase) :: this
real(kind=8) ans,dof
dof=this%mean
ans=rand_gamma(0.5d0, 2.0d0*dof,this)
END FUNCTION

!
! ## return a random sample from an inverse gamma random variable
!
FUNCTION rand_inverse_gamma(this) RESULT(ans)
class(ranbase) :: this
real(kind=8) SHAPE,scale,ans
SHAPE=this%mean
scale=this%sigma

! ## If X is gamma(shape, scale) then
! ## 1/Y is inverse gamma(shape, 1/scale)
ans= 1.0d0 / rand_gamma(shape, 1.0d0 / SCALE, this)
END FUNCTION
!
!## return a sample from a Weibull distribution
!

FUNCTION rand_weibull(this) RESULT(ans)
class(ranbase) :: this
real(kind=8) SHAPE,scale,temp,ans
SHAPE=this%mean
scale=this%sigma
IF (shape <= 0.0d0) THEN

   WRITE(output_unit,'(a)') "Shape PARAMETER must be positive"
END IF
IF (scale <= 0.0d0) THEN

   WRITE(output_unit,'(a)') "Scale PARAMETER must be positive"
END IF
CALL this%rand_get(temp)
ans= SCALE * (-log(temp))**(1.0 / SHAPE)
END FUNCTION

!
!## return a random sample from a Cauchy distribution
!
FUNCTION rand_cauchy(this) RESULT(ans)
class(ranbase) :: this
real(kind=8) ans,median,scale,p
median=this%mean
scale=this%sigma

IF (scale <= 0.0d0) THEN

   WRITE(output_unit,'(a)') "Scale PARAMETER must be positive"
END IF
CALL this%rand_get(p)
ans = median + SCALE*tan(PI*(p - 0.5))
END FUNCTION

!
!## return a random sample from a Student t distribution
!
!FUNCTION rand_student_t(dof) RESULT(ans)
FUNCTION rand_student_t(this) RESULT(ans)
class(ranbase) :: this
real(kind=8) ans,dof,y1,y2
dof=this%mean
IF (dof <= 0.d0) THEN

   WRITE(output_unit,'(a)') "Degrees of freedom must be positive"
END IF
!
! ## See Seminumerical Algorithms by Knuth
y1 = rand_normal_int(0.0d0, 1.0d0,this)
y2 = rand_chi_square(this)
ans= y1 / (y2 / DOf)**0.50d0
!

END FUNCTION

!
!## return a random sample from a Laplace distribution
!## The Laplace distribution is also known as the double exponential distribution.
!
!FUNCTION rand_laplace(mean, SCALE)  RESULT(ans)
FUNCTION rand_laplace(this)  RESULT(ans)
class(ranbase) :: this
real(kind=8) ans,mean,scale,u
mean=this%mean
scale=this%sigma
IF (scale <= 0.0d0) THEN

   WRITE(output_unit,'(a)') "Scale PARAMETER must be positive"
END IF
CALL this%rand_get(u)
IF (u < 0.5d0) THEN

   ans = mean + SCALE*log(2.0*u)
ELSE
   ans = mean - SCALE*log(2*(1-u))
END IF

END FUNCTION

!
! ## return a random sample from a log-normal distribution
!
FUNCTION rand_log_normal(this) RESULT(ans)
class(ranbase) :: this
real(kind=8) ans
ans= EXP(rand_normal(this))
END FUNCTION

!
! ## return a random sample from a beta distribution
!
FUNCTION rand_beta(this) RESULT(ans)
class(ranbase) :: this
real(kind=8) a,b,ans,u,v
a=this%mean
b=this%sigma
IF ((a <= 0.0d0) .OR. (b <= 0.0d0)) THEN

   WRITE(output_unit,'(a)') "Beta PARAMETERs must be positive"
END IF

! ## There are more efficient methods for generating beta samples.
! ## However such methods are a little more efficient and much more complicated.
! ## For an explanation of why the following method works, see
! ## http://www.johndcook.com/distribution_chart.html#gamma_beta

u = rand_gamma(a, 1.0d0,this)
v = rand_gamma(b, 1.0d0,this)
ans = u / (u + v)
END FUNCTION

END MODULE m_random_number_library
