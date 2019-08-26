MODULE deviates_module
  USE rnd_array_module
  IMPLICIT NONE
!... For Gaussian
  INTEGER :: iset=0
  REAL(KIND=8) :: gset,x1_gauss,x2_gauss,xln_gauss,xlnsphere_r
  LOGICAL :: force_iset=.false.
!... For Lorentzian
  REAL(KIND=8) :: xold=-1.d0

CONTAINS

  SUBROUTINE uniform_pm1(x)
    REAL(KIND=8) :: x
    CALL wrap_one_rnd(x)
    x=2.d0*(x-0.5d0)
    RETURN
  END SUBROUTINE uniform_pm1
  
  SUBROUTINE uniform_range(x,xmin,xmax)
    REAL(KIND=8) :: x,xmin,xmax
    CALL wrap_one_rnd(x)
    x=xmin+(xmax-xmin)*x
    RETURN
  END SUBROUTINE uniform_range

  SUBROUTINE uniform3d(length,vcen,v)
    REAL(KIND=8) :: length,vcen(3),v(3),x1,x2,x3
    CALL wrap_one_rnd(x1)
    CALL wrap_one_rnd(x2)
    CALL wrap_one_rnd(x3)
    v(1)=length*(x1-0.5d0)+vcen(1)
    v(2)=length*(x2-0.5d0)+vcen(2)
    v(3)=length*(x3-0.5d0)+vcen(3)
    RETURN
  END SUBROUTINE uniform3d

  SUBROUTINE unidisc_deviate(radius,v2cen,v2)
    REAL(KIND=8) :: radius,v2cen(2),v2(2),x1,x2,r,x,y
    REAl(KIND=8) :: tau=6.283185307179586476925d0! 2*pi
    CALL wrap_one_rnd(x1)
    CALL wrap_one_rnd(x2)
    r=sqrt(x1)
    x=r*cos(tau*x2)
    y=r*sin(tau*x2)
    v2(1)=radius*x+v2cen(1)
    v2(2)=radius*y+v2cen(2)
    RETURN
  END SUBROUTINE unidisc_deviate

  SUBROUTINE unisphere_deviate(radius,vcen,v)
    REAL(KIND=8) :: radius,vcen(3),v(3),x1,x2,x3,r,x,y,z
    REAl(KIND=8) :: tau=6.283185307179586476925d0! 2*pi
    REAl(KIND=8) :: third=1.d0/3.d0
    LOGICAL :: issrfc
    issrfc=radius<0
    CALL wrap_one_rnd(x1)
    CALL wrap_one_rnd(x2)
    CALL wrap_one_rnd(x3)
    IF (issrfc) THEN
       r=1.d0
    ELSE
       r=x1**third
    ENDIF
    z=(1.d0-2.d0*x2)*r
    x=sqrt(r*r-z*z)*cos(tau*x3)
    y=sqrt(r*r-z*z)*sin(tau*x3)
    v(1)=abs(radius)*x+vcen(1)
    v(2)=abs(radius)*y+vcen(2)
    v(3)=abs(radius)*z+vcen(3)
    RETURN
  END SUBROUTINE unisphere_deviate

  SUBROUTINE gauss_deviate(xgauss)
!    USES gfortran intrinsic PRNG
!    integer :: aseed(12),i!debug
    REAL(KIND=8) :: xgauss
    REAL(KIND=8) :: fac,rsq,v1,v2
    REAL(KIND=8) :: x1,x2
!    iset=0
    if (iset.eq.0.or.force_iset) then
1      CALL wrap_one_rnd(x1)
       CALL wrap_one_rnd(x2)
!       call random_seed(get=aseed)
!       print *,'rnd_call+inicount=',rnd_callcount,rnd_initcount
!       do i=1,12
!          print *,'GAUSS_DEVIATE: i,aseed =',i,aseed(i)
!       enddo
! export
       x1_gauss=x1
       x2_gauss=x2
!       
       v1=2.d0*x1-1.d0
       v2=2.d0*x2-1.d0
       rsq=v1**2+v2**2
       if (rsq.ge.1.d0.or.rsq.eq.0.d0) goto 1
       fac=sqrt(-2.d0*log(rsq)/rsq)
       gset=v1*fac
       xgauss=v2*fac
       iset=1
    else
       xgauss=gset
       iset=0
!       print *,'x1,x2,xgauss =',x1,x2,xgauss
    endif
    return
  END SUBROUTINE gauss_deviate

  SUBROUTINE lognormal(u,sigma,xln)
    REAL(KIND=8) :: u,sigma,xsn,xn,xln
    CALL gauss_deviate(xsn)
    xn=u+sigma*xsn
    xln=exp(xn)
    xln_gauss=xsn
    RETURN
  END SUBROUTINE lognormal
  
  SUBROUTINE lognormal_hybride(u,sigma_abs,sigma_rel,xhybride)
!---- Hybride lognormal + normal distribution
!     Note that u must be greater than 0.    
    REAL(KIND=8) :: u,lu,sigma_abs,sigma_rel,xsn,xln,xhybride
    lu=log(u)
    CALL lognormal(lu,sigma_rel,xln)
    CALL gauss_deviate(xsn)
    xhybride=xln+sigma_abs*xsn
    RETURN
  END SUBROUTINE lognormal_hybride
  
  SUBROUTINE lognormal_sphere(u,sigma,v)
    REAl(KIND=8) :: tau=6.283185307179586476925d0! 2*pi
    REAL(KIND=8) :: u,sigma,v(3),x1,x2,xln,x,y,z
    call lognormal(u,sigma,xln)
    CALL wrap_one_rnd(x1)
    CALL wrap_one_rnd(x2)
    z=(1.d0-2.d0*x1)
    x=sqrt(1.d0-z*z)*cos(tau*x2)
    y=sqrt(1.d0-z*z)*sin(tau*x2)
    v(1)=xln*x
    v(2)=xln*y
    v(3)=xln*z
    xlnsphere_r=xln
    RETURN
  END SUBROUTINE lognormal_sphere

  SUBROUTINE gauss3d(sigma,v)
!---- 3D Gaussian distribution of one random particle.
!     Becomes Maxwell-Boltzmann for sigma = sqrt(k*T/mu) = sqrt(R*T/M)
    REAL(KIND=8) :: xgauss(3),sigma,v(3)!,vsq(n)
    CALL gauss_deviate(xgauss(1))!; v(1)=xgauss*sigma
    CALL gauss_deviate(xgauss(2))!; v(2)=xgauss*sigma
    CALL gauss_deviate(xgauss(3))!; v(3)=xgauss*sigma
    v(:)=xgauss(:)*sigma
    RETURN
  END SUBROUTINE gauss3d

  SUBROUTINE gauss_truncated(maxsigma,xgauss)
    REAL(KIND=8) :: maxsigma,xgauss
    xgauss=2.d0*maxsigma
    IF (maxsigma<=0.d0) THEN!criterion turned off
       CALL gauss_deviate(xgauss)
       RETURN
    ENDIF
    DO WHILE (xgauss>maxsigma)
       CALL gauss_deviate(xgauss)
    END DO
    RETURN
  END SUBROUTINE gauss_truncated

  SUBROUTINE gauss3d_truncated(sigma,maxsigma,v)
!---- 3D Gaussian distribution of one random particle.
!     Becomes Maxwell-Boltzmann for sigma = sqrt(k*T/mu) = sqrt(R*T/M)
    REAL(KIND=8) :: xgauss(3),sigma,maxsigma,v(3)!,vsq(n)
    CALL gauss_truncated(maxsigma,xgauss(1))!; v(1)=xgauss*sigma
    CALL gauss_truncated(maxsigma,xgauss(2))!; v(2)=xgauss*sigma
    CALL gauss_truncated(maxsigma,xgauss(3))!; v(3)=xgauss*sigma
    v(:)=xgauss(:)*sigma
    RETURN
  END SUBROUTINE gauss3d_truncated

  SUBROUTINE swexp_gauss(u,sigma,xsw)
    REAL(KIND=8) :: u,sigma,xsw,xgauss,xsn
    CALL gauss_deviate(xgauss)
    xsw=1.d0/(exp(u-sigma*xgauss)+1.d0)
    RETURN
  END SUBROUTINE swexp_gauss
  
  SUBROUTINE swexp_gauss_borders(u,sigma,a,b,xsw)
    REAL(KIND=8) :: u,sigma,a,b,xsw,xgauss,xdum
    CALL gauss_deviate(xgauss)
    xdum=1.d0/(exp(u-sigma*xgauss)+1.d0)
    xsw=xdum*(b-a)+a
    RETURN
  END SUBROUTINE swexp_gauss_borders

  SUBROUTINE lorentz_deviate(xlorentz)
    REAL(KIND=8) :: xlorentz,xgauss1,xgauss2
1   CALL gauss_deviate(xgauss1)
    IF (xold<=0.d0) THEN
       CALL gauss_deviate(xgauss2)
    ELSE
       xgauss2=xold
    ENDIF
    xold=xgauss1!recycle
    IF (xgauss2==0.d0) GOTO 1
! The ratio of two independent Gaussian deviates is Lorentz-distributed
    xlorentz=xgauss1/xgauss2
    RETURN
  END SUBROUTINE lorentz_deviate

  SUBROUTINE logistic_deviate(xlogistic)
    REAL(KIND=8) :: u,sigma,xlogistic
    REAL(KIND=8) :: x
    CALL wrap_one_rnd(x)
    xlogistic=log(x/(1.-x))
    RETURN
  END SUBROUTINE logistic_deviate

  SUBROUTINE sech_deviate(xsech)
    REAL(KIND=8) :: x,xsech
    REAL(KIND=8) :: pi=3.141592653589793238D0
    CALL wrap_one_rnd(x)
    xsech=2.d0/pi*log(tan(pi/2.d0*x))
    RETURN
  END SUBROUTINE sech_deviate

!**** Distributions with additional parameters ****

  SUBROUTINE voigt_deviate(sigma,gamma,u,v,xgauss,xlorentz,xvoigt)
! A Voigt profile is a convolution of Gaussian with Lorentzian profile
    REAL(KIND=8) :: xgauss,xlorentz,xvoigt,sigma,gamma,u,v
    CALL gauss_deviate(xgauss)
    CALL lorentz_deviate(xlorentz)
    xgauss=xgauss*sigma+u
    xlorentz=xlorentz*gamma+v
    xvoigt=xgauss+xlorentz
    RETURN
  END SUBROUTINE voigt_deviate

  SUBROUTINE weibull_distribution(lambda,k,xweibull)
    real(kind=8) :: xflat,xweibull,lambda,k
    call wrap_one_rnd(xflat)
    xweibull=1./lambda*(-log(1.d0-xflat))**(1.d0/k)
    return
  END SUBROUTINE weibull_distribution

END MODULE deviates_module
