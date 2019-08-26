!==== vec3arithm
      SUBROUTINE smult3(a,k,v)
!---- Skalaren-Multiplikation
!     reflexiv moeglich
      IMPLICIT NONE
      DOUBLE PRECISION a(3),k,v(3)
      v(1) = k*a(1)
      v(2) = k*a(2)
      v(3) = k*a(3)
      RETURN
      END

      SUBROUTINE psca3(a,b,s)
!---- Skalarprodukt
      IMPLICIT NONE
      DOUBLE PRECISION a(3),b(3),s
      s = a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
      RETURN
      END

      SUBROUTINE pvec3(a,b,v)
!---- Vektorprodukt a x b = v
      IMPLICIT NONE
      DOUBLE PRECISION a(3),b(3),v(3)
      v(1) = a(2)*b(3)-a(3)*b(2)
      v(2) = a(3)*b(1)-a(1)*b(3)
      v(3) = a(1)*b(2)-a(2)*b(1)
      RETURN
      END

      SUBROUTINE vsum3(a,b,v)
!---- Vektorsumme
!     reflexiv moeglich
      IMPLICIT NONE
      DOUBLE PRECISION a(3),b(3),v(3)
      v(1) = a(1)+b(1)
      v(2) = a(2)+b(2)
      v(3) = a(3)+b(3)
      RETURN
      END

      SUBROUTINE vdiff3(a,b,v)
!---- Vektordifferenz (wie vdn3, aber einfacher)
!     reflexiv moeglich
      IMPLICIT NONE
      DOUBLE PRECISION a(3),b(3),v(3)
      v(1) = a(1)-b(1)
      v(2) = a(2)-b(2)
      v(3) = a(3)-b(3)
      RETURN
      END

      SUBROUTINE vksum3(a,b,ka,kb,v)
!---- Vektorsumme mit Skalarfaktor
!     reflexiv moeglich
      IMPLICIT NONE
      DOUBLE PRECISION a(3),b(3),ka,kb,v(3)
      v(1) = ka*a(1)+kb*b(1)
      v(2) = ka*a(2)+kb*b(2)
      v(3) = ka*a(3)+kb*b(3)
      RETURN
      END

!==== vec3norm

      SUBROUTINE vnsq3(a,absa,asq)
!---- Standard-Vektornorm
      IMPLICIT NONE
      DOUBLE PRECISION a(3),absa,asq
      asq = a(1)**2+a(2)**2+a(3)**2
      absa=sqrt(asq)
      RETURN
      END

      SUBROUTINE vdn3(a,b,v,d,dsq)
!---- Vektordistanz
!     reflexiv moeglich
      IMPLICIT NONE
      DOUBLE PRECISION a(3),b(3),v(3),d,dsq
      v(1) = a(1)-b(1)
      v(2) = a(2)-b(2)
      v(3) = a(3)-b(3)
      CALL vnsq3(v,d,dsq)
      RETURN
      END

      SUBROUTINE ev3(a,k,kea,absa,asq)
!---- Skalierter Einheitsvektor
!     kea = k*e_a
      IMPLICIT NONE
      DOUBLE PRECISION a(3),k,ka,kea(3),absa,asq
      CALL vnsq3(a,absa,asq)
!      CALL smult3(a,k/absa,kea)
      ka=k/absa
      kea(1) = ka*a(1)
      kea(2) = ka*a(2)
      kea(3) = ka*a(3)
      RETURN
      END

!==== vec3rot
      SUBROUTINE vxrot(vec0,vec1,phi)
!---- Drehung um x-Achse
      IMPLICIT NONE
      DOUBLE PRECISION vec0(1:3),vec1(1:3),phi,c,s
      c = cos(phi)
      s = sin(phi)
      vec1(1) =   vec0(1)
      vec1(2) =             c*vec0(2) - s*vec0(3)
      vec1(3) =             s*vec0(2) + c*vec0(3)
      RETURN
      END

      SUBROUTINE vyrot(vec0,vec1,phi)
!---- Drehung um y-Achse
      IMPLICIT NONE
      DOUBLE PRECISION vec0(1:3),vec1(1:3),phi,c,s
      c = cos(phi)
      s = sin(phi)
      vec1(1) = c*vec0(1)           + s*vec0(3)
      vec1(2) =             vec0(2)
      vec1(3) =-s*vec0(1)           + c*vec0(3)
      RETURN
      END

      SUBROUTINE vzrot(vec0,vec1,phi)
!---- Drehung um z-Achse
      IMPLICIT NONE
      DOUBLE PRECISION vec0(1:3),vec1(1:3),phi,c,s
      c = cos(phi)
      s = sin(phi)
      vec1(1) = c*vec0(1) - s*vec0(2)
      vec1(2) = s*vec0(1) + c*vec0(2)
      vec1(3) =                         vec0(3)
      RETURN
      END

      SUBROUTINE vnrot(v,vz,vnew,theta,phi)
!---- Rotation via neuem z-Vektor: Vektor in xyz-System wird so gedreht,
!     dass die neue x'y'-Ebene den Normalenvektor vz hat.
!     Verwendet vec3rot.f
      IMPLICIT NONE
      DOUBLE PRECISION v(3),vz(3),vnew(3),x,y,z,r,theta,phi
      DOUBLE PRECISION v1(3)
      x=vz(1)
      y=vz(2)
      z=vz(3)
      CALL ct2sph(x,y,z,r,theta,phi)
      CALL vyrot(v,v1,theta)
      CALL vzrot(v1,vnew,phi)
      RETURN
      END

      SUBROUTINE rotzxz(vec0,vec1,alpha,beta,gamma)
!---- Euler rotation extrinsic forward
!     (swap alpha and gamma for intrinsic rotation)
      IMPLICIT NONE
      DOUBLE PRECISION vec0(1:3),vec1(1:3),alpha,beta,gamma
      DOUBLE PRECISION veca(1:3),vecb(1:3)
      CALL vzrot(vec0,veca,alpha)
      CALL vxrot(veca,vecb,beta)
      CALL vzrot(vecb,vec1,gamma)
      RETURN
      END

      SUBROUTINE torzxz(vec0,vec1,alpha,beta,gamma)
!---- Euler rotation extrinsic backward
!     (swap alpha and gamma for intrinsic rotation)
      IMPLICIT NONE
      DOUBLE PRECISION vec0(1:3),vec1(1:3),alpha,beta,gamma
      DOUBLE PRECISION veca(1:3),vecb(1:3)
      CALL vzrot(vec0,veca,-gamma)
      CALL vxrot(veca,vecb,-beta)
      CALL vzrot(vecb,vec1,-alpha)
      RETURN
      END

      SUBROUTINE rotzyz(vec0,vec1,alpha,beta,gamma)
!---- Euler rotation extrinsic forward
!     (swap alpha and gamma for intrinsic rotation)
      IMPLICIT NONE
      DOUBLE PRECISION vec0(1:3),vec1(1:3),alpha,beta,gamma
      DOUBLE PRECISION veca(1:3),vecb(1:3)
      CALL vzrot(vec0,veca,alpha)
      CALL vyrot(veca,vecb,beta)
      CALL vzrot(vecb,vec1,gamma)
      RETURN
      END

      SUBROUTINE torzyz(vec0,vec1,alpha,beta,gamma)
!---- Euler rotation extrinsic backward
!     (swap alpha and gamma for intrinsic rotation)
      IMPLICIT NONE
      DOUBLE PRECISION vec0(1:3),vec1(1:3),alpha,beta,gamma
      DOUBLE PRECISION veca(1:3),vecb(1:3)
      CALL vzrot(vec0,veca,-gamma)
      CALL vyrot(veca,vecb,-beta)
      CALL vzrot(vecb,vec1,-alpha)
      RETURN
      END

!==== vec3ctf

      SUBROUTINE sph2ct(r,theta,phi,x,y,z)
!---- Spherical to Cartesian coordinates (theta from pole)
      IMPLICIT NONE
      DOUBLE PRECISION r,phi,theta,x,y,z
      x = r*cos(phi)*sin(theta)
      y = r*sin(phi)*sin(theta)
      z = r*cos(theta)
      RETURN
      END

      SUBROUTINE ct2sph(x,y,z,r,theta,phi)
!---- Cartesian to Spherical coordinates (theta from pole)
      IMPLICIT NONE
      DOUBLE PRECISION r,rho,phi,theta,x,y,z,xn
      r    = sqrt(x*x+y*y+z*z)
      rho  = sqrt(x*x+y*y)
      theta= acos(z/r)
      IF (rho.gt.0.d0) THEN
         xn   = max(min(x/(rho),1.d0),-1.d0)
      ELSE
         xn   = 0.d0
      ENDIF
      phi  = acos(xn)*sign(1.d0,y)
      RETURN
      END

!==== vec3radial (NEW)

      SUBROUTINE vrad3(v,vr,rsq,tsq)
!     v  = input vector
!     vr = radial reference vector
!     rsq= squared radial norm
!     tsq= squared tangential norm
!     Uses vec3arithm.f, vec3norm.f
      IMPLICIT NONE
      DOUBLE PRECISION v(3),vr(3),rsq,tsq, vsq,absvr,s,ddum
      CALL vnsq3(v,ddum,vsq)
      CALL vnsq3(vr,absvr,ddum)
      CALL psca3(v,vr,s)
      rsq=(s/absvr)**2
      tsq=vsq-rsq
      RETURN
      END
