!
      SUBROUTINE wrap_xysmooth(mode,bmod,x,y,x0,xn,y0,yn,nbx,nby,hx,hy,        &
     &kw,  xtab,ytab,qbin)
      IMPLICIT NONE
      INTEGER nbmax
      PARAMETER (nbmax=1200)
      INTEGER mode,bmod,nbx,nby
!      DOUBLE PRECISION :: pi4=7.8539816339744831d-1
      DOUBLE PRECISION x,y,x0,xn,y0,yn,hx,hy,kw
      DOUBLE PRECISION xtab(nbmax),ytab(nbmax),qbin(nbmax,nbmax)
      IF (hx.le.0.d0.or.hy.le.0.d0) THEN
! factor of kw to approx. match xysmooth levels
         CALL xybin(mode,bmod,x,y,x0,xn,y0,yn,nbx,nby,kw,                      &
     &        xtab,ytab,qbin)
      ELSE
         CALL xysmooth(mode,bmod,x,y,x0,xn,y0,yn,nbx,nby,hx,hy,                &
     &        kw,  xtab,ytab,qbin)
      ENDIF
      RETURN
      END

! Secondary composition file - don't include new features here but in original source files!
!
! xybin052.f
      SUBROUTINE xybin(mode,bmod,x,y,x0,xn,y0,yn,nbx,nby,kw,            &
     &xtab,ytab,qbin)
!******** 2-D BINNING ALGORITHM v 0.52***
!     Based on xbin v 0.401
!     Minor change in xybin.dim.f (2013-10-16): increased nbmax from 1000 to 2000
!     Changes since v.0.51:
!     - moved qcnt and iqb to optional COMMON block
!     - added qref as the bin midpoint closes to (x,y)
!
!     XYBIN Calculates for both x and y direction the number of events
!     within interval partition
!     [u_i,u_i+1[, where u_i = x0 + i*(xn-x0)/nb
!     for i=0...nb
!
!     Example for 4 bins and x0=-10, x1=10
!
!     bin#  | 1 | 2 | 3 | 4 |
!     x    -10 -5   0   5   10
!     u_i   0   1   2   3   4
!
!     Bins are labelled as follows:
!     Bin #j := [x_j-db/2, x_j+db/2[ where x_j = (u_{j-1}+u_j)/2
!     and db the bin width. u_j are not explicitely used in the routine.
!
!     OPTIONS
!     mode = setup option: 0 - setup bin positions ybin(1,*) (first each)
!                          1 - re-use from given ybin
!                         -1 - just get bin index (and setup)
!     bmod = boundary option: 0 - strictly excludes right interval sup
!                             1 - counts event if equal to sup of nth bin
!     INPUT:
!     x,y    = event abscissa
!     x0,xn  = overall interval (from left side of first bin to right side of last bin)
!     y0,yn  = same, for y
!     nbx,nby= number of x,y partitions
!     kw     = weighting factor (normally = 1)
!     qcnt   = see OUTPUT; qcnt(0) used here as qcnt weighting switch
!                             >=1: kb=kw*db
!                            else: kb=kw
!
!     OUTPUT:
!     xtab = Midpoint abscissa table in x
!     ytab = Midpoint abscissa table in y
!     qbin = output array; (i,j) = content of bin at position (i,j)
!     qcnt = (0) weighted (by kw AND bin width) total number outside qbin
!            (1)    "       "     "    within range
!     ibq  = index of actual bin
!            Boundary: ibq(1)=0 for x<u_0, ibq(1)=nbx+1 for x>(=) u_nbx
!                      ibq(2)=0 for y<v_0, ibq(2)=nby+1 for x>(=) v_nby
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER nbmax
      PARAMETER (nbmax=1200)
      INTEGER nbx,nby,mode,bmod,jbx,jby!,ncall
      DOUBLE PRECISION x,y,x0,xn,y0,yn,kw,qbin(nbmax,nbmax)
      DOUBLE PRECISION u,v,dbx,dby,kb,xtab(nbmax),ytab(nbmax)
!---- Auxiliary variables export
      DOUBLE PRECISION qcnt(0:1),qref(2)
      INTEGER ibq(2)
      LOGICAL isinq
      COMMON /xybaux/ qcnt,qref,ibq,isinq
!---- Debug variables import/export
      LOGICAL dotnk,verbous,dokb
      COMMON /xbntnk/ dotnk,verbous
      INTEGER i,j!for dotnk
!---- Optional call count
!      SAVE ncall
!      DATA ncall/0/
!      ncall=ncall+1
      IF (x0.ge.xn.OR.y0.ge.yn) THEN
         WRITE(*,*) 'XBIN ERROR: xy0 < xyn violated! --> Abort'
         WRITE(*,*) 'x0 = ',x0,'          xn = ',xn
         WRITE(*,*) 'y0 = ',y0,'          xn = ',yn
         STOP
      ENDIF
!---- Setup
      dbx=(xn-x0)/nbx!calculate this always (e.g. for qbin's of different size)
      dby=(yn-y0)/nby
      IF (dbx.le.0.d0.or.dby.le.0.d0) THEN
         WRITE(*,*) 'XYBIN Error: stepsize must be positive!'
         WRITE(*,*) 'x0,xn,nbx,dbx =',x0,xn,nbx,' ',dbx
         WRITE(*,*) 'y0,yn,nby,dby =',y0,yn,nby,' ',dby
         STOP
      ENDIF
      IF (mode.eq.0) THEN
         CALL setxyb(x0,xn,y0,yn,nbx,nby,dbx,dby,                       &
     &        xtab,ytab,qbin,qcnt,dokb, dotnk,verbous)
      ENDIF
      IF (dokb) THEN
         kb=kw*dbx*dby
      ELSE
         kb=kw
      ENDIF
!---- Daily use
      u=x-x0
      v=y-y0
      ibq(1)=max(min(int(u/dbx+1.d0),nbx+1),0)
      ibq(2)=max(min(int(v/dby+1.d0),nby+1),0)
      CALL getbin(x,x0,dbx,ibq(1),qref(1))
      CALL getbin(y,y0,dby,ibq(2),qref(2))
      IF (bmod.eq.1 .and. x.eq.xn) ibq(1)=nbx
      IF (bmod.eq.1 .and. y.eq.yn) ibq(2)=nby
      isinq=     ibq(1).ge.1.and.ibq(1).le.nbx                          &
     &     .AND.ibq(2).ge.1.and.ibq(2).le.nby
!      print *,'x,y,ibq =',x,y,ibq,isinq
c.... Return if index-only mode is set
      IF (mode.lt.0.or.kw.eq.0.d0) RETURN
c.... Continue else 
!      if(x0.eq.0.d0) then
      if (dotnk) then
         WRITE(*,*) 'x,y,xref,yref =',x,y,qref(1),qref(2)
         jbx=max(min(ibq(1),nbx),1)
         jby=max(min(ibq(2),nby),1)
!         write(*,1000) 'x,x0,xn,xtab,db,nb,ib =',                       &
!     &        x,x0,xn,xtab(jbx),dbx,nbx,nbx,ibq(1)
!         write(*,1000) 'y,y0,yn,ytab,db,nb,ib =',                       &
!     &        y,y0,yn,ytab(jby),dby,nby,nby,ibq(2)
!         if(ib.eq.0) pause
!         if(ib.le.1) write(*,*) 'ib,x = ',ib,x
!         write(*,*) 'ibq(1,2),nb = ',ibq(1),ibq(2)!,nbx,nby
!         write(*,*) 'ibq(2),y = ',ibq(2),y
      endif
      IF (isinq) THEN
         qcnt(1)=qcnt(1)+kb
!         qbin(ibq(1),ibq(2))=qbin(ibq(1),ibq(2))+kw
         qbin(ibq(1),ibq(2))=qbin(ibq(1),ibq(2))+kw
         IF (dotnk) WRITE(*,2000) 'Inside:  x,y,ibq(1),ibq(2) =',       &
     &    x,y,ibq(1),ibq(2)
      ELSE
         qcnt(0)=qcnt(0)+kb
!         IF (dotnk) WRITE(*,2000) 'Outside: x,y,ibq(1),ibq(2) =',       &
!     &        x,y,ibq(1),ibq(2)
      ENDIF
      RETURN
      IF (.not.dotnk) RETURN
!---- Debugging stuff
!      WRITE(*,3000)'Q: ibq(1),ibq(2),qb,qc =',                          &
!     &     ibq(1),ibq(2),qbin(ibq(1),ibq(2)),qcnt(1)
      DO 200 j=1,nby
         DO 100 i=1,nbx
            WRITE(*,'(2(I3),3(F16.8,2X))') i,j,xtab(i),ytab(j),qbin(i,j)
 100     CONTINUE
         WRITE(*,*)
 200  CONTINUE
 1000 FORMAT(A,1X,5(F9.3,2X),2(2X,I6))
 2000 FORMAT(A,1X,2(F9.3,2X),2(2X,I6),2X,F12.6)
 3000 FORMAT(A,1X,2(2X,I6),3(2X,F12.6))
      RETURN
      END

      SUBROUTINE setxyb(x0,xn,y0,yn,nbx,nby,dbx,dby,                    &
     &     xtab,ytab,qbin,qcnt,dokb, dotnk,verbous)
      IMPLICIT NONE
      INTEGER nbmax
      PARAMETER (nbmax=1200)
      INTEGER nbx,nby,i,j
      DOUBLE PRECISION x0,xn,y0,yn,dbx,dby,qbin(nbmax,nbmax),           &
     &xtab(nbmax),ytab(nbmax),qcnt(0:1)
      LOGICAL dokb,verbous,dotnk
      IF (verbous) write(*,*) 'XBIN: Initializing...'
      DO 12 j=1,nby
         ytab(j)=y0+(dble(j)-0.5d0)*dby
         if (dotnk) then
!            write(*,*) 'j,y0,dby,ytab =',j,y0,dby,ytab(j)
         endif
         DO 11 i=1,nbx
!--- xtab needs to be set up only once --> test for i=1
            IF (j.eq.1) xtab(i)=x0+(dble(i)-0.5d0)*dbx
!--- Calculation of bin midpoints
!     xtab(i)   x abscissa equal to bin midpoint
!     ytab(j)   y abscissa equal to bin midpoint
!     qbin(i,j) use small non-zero value to avoid log crashing
!            qbin(i,j)=1.d-99
!     or just use zero instead
            qbin(i,j)=0.d0
            if (dotnk) then
!               write(*,*) 'i,xi,yi =',i,xtab(i),ytab(j)
            endif
 11      CONTINUE
 12   CONTINUE
      dokb=qcnt(0).ge.1.d0
      qcnt(0)=0.d0
      qcnt(1)=0.d0
      RETURN
      END

      SUBROUTINE getbin(s,s0,db,ib,sib)
      INTEGER ib
      DOUBLE PRECISION s,s0,db,dbint,sib
      DOUBLE PRECISION t
      t=s-s0
      dbint=aint(t/db)
      IF (db.NE.0.D0) THEN
         IF (t/db.EQ.dbint.OR.t.GE.0.D0) THEN
            ib = dbint+1.D0
         ELSE
            ib = dbint
         ENDIF
      ELSE
         ib = t
      ENDIF
      sib=s0+(dble(ib)-0.5d0)*db
!      write(*,*) 'sib =',sib
      RETURN
      END

!#########################################################################

!xysmooth052.f
      SUBROUTINE xysmooth(mode,bmod,x,y,x0,xn,y0,yn,nbx,nby,hx,hy,        &
     &kw,  xtab,ytab,qbin)
      IMPLICIT NONE
      INTEGER nbmax
      PARAMETER (nbmax=1200)
      INTEGER ip,np,npmax,nwmax
      PARAMETER (npmax=1002001)
      INTEGER mode,bmod,qmod,nbx,nby,npxy(0:3)
      DOUBLE PRECISION x,y,xp,yp,x0,xn,y0,yn,h,hx,hy,kw,kp,dx,dy
      DOUBLE PRECISION xtab(nbmax),ytab(nbmax),qbin(nbmax,nbmax)
      DOUBLE PRECISION posmxy(0:2,npmax)
!---- Auxiliary variables export
      DOUBLE PRECISION qcnt(0:1),qref(2)
      INTEGER ibq(2)
      LOGICAL isinq
      COMMON /xybaux/ qcnt,qref,ibq,isinq
      LOGICAL dotnk,verbous
      COMMON /xbntnk/ dotnk,verbous
!---- 
      SAVE dx,dy
      DATA nwmax/1001/!compatibility
      IF (mode.le.0) THEN
         WRITE(*,*) 'Initializing SMOOTHXYBIN'
         dx=(xn-x0)/nbx
         dy=(yn-y0)/nby
      ENDIF
      IF (hx.le.0.d0.AND.hy.le.0.d0) THEN
         WRITE(*,*) 'XYSMOOTH ERROR: hxy must be > 0'
         WRITE(*,*) 'hx,hy =',hx,hy
         STOP
      ELSE IF (hy.le.0.d0.OR.hy.le.0.d0) THEN
         WRITE(*,*) 'XYSMOOTH WARNING: hxy must be > 0'
         WRITE(*,*) 'hx,hy =',hx,hy
         hx=max(hx,hy)
         hy=hx
         WRITE(*,*) 'Using hx = hy =',hx
      ENDIF
      CALL xybin(mode,bmod,x,y,x0,xn,y0,yn,nbx,nby,0.d0,                &
     &     xtab,ytab,qbin)
               
c      CALL spreadwkernxy(x,y,qref(1),qref(2),dx,dy,h,kw,nwmax,          &
c     &     posmxy,npxy)
      CALL smoothwkern2d(x,y,qref(1),qref(2),dx,dy,hx,hy,kw,             &
     &     posmxy,npxy)
      np=npxy(0)
      DO ip=1,np
         xp=posmxy(1,ip)
         yp=posmxy(2,ip)
         kp=posmxy(0,ip)
         CALL xybin(mode,bmod,xp,yp,x0,xn,y0,yn,nbx,nby,kp,             &
     &        xtab,ytab,qbin)
         mode=1
      ENDDO
      RETURN
      END

!     2D spread of wml85 and other kernels (SPH)
      SUBROUTINE smoothwkern2d(x,y,xref,yref,dx,dy,hx,hy,m,pmxy,npxy)
!     New variation of spreadkernxy that use individual hx, hy
!     instead of a common h.
!     INPUT
!     x,y    : x,y position of kernel (e.g. SPH)
!     x/yref : x,y reference position (usually the center of a nearby bin)
!     dx,dy  : x,y bin width
!     hx,hy  : kernel radius (=1/2 wml radius) in x and y direction
!     m      : mass or any other quantity of the kernel
!     OUTPUT
!     pmxy(0:2,1:np): list of subkernels into the kernel has been split up
!                     0: kernel value
!                     1: px
!                     2: py
!     npxy(0)       : np = total number of subkernels = nx*ny
!         (1)       : nx = number of x columns
!         (2)       : ny = number of y rows
!
!     USES wml85.f
      IMPLICIT NONE
      INTEGER i,j,p,pin,npmax,nxymax,NMAXX,nx,ny,np,npxy(0:3)
      PARAMETER (NMAXX=1002001,nxymax=1001)
!      PARAMETER (NMAXX=4004001,nxymax=2001)
      INTEGER ndim,wsel
      COMMON /wkpars/ ndim,wsel
!      DOUBLE PRECISION wml85
      DOUBLE PRECISION wkernel
      DOUBLE PRECISION x,y,xref,yref,dx,dy,hx,hy,m,pmxy(0:2,NMAXX)
      DOUBLE PRECISION xmin,xmax,ymin,ymax,xlo,xhi,ylo,yhi
      DOUBLE PRECISION rxy,rhx,rhy,hu,rhmaxx,rhmaxy
      DOUBLE PRECISION px,py,mxy,snorm           !,unorm
!.... wkernel.f common block
      DOUBLE PRECISION rsigma(-1:3,0:7),hrad(0:7)
      COMMON /waux/ rsigma,hrad
      COMMON /smkexport/ rhx,rhy,hu,xlo,ylo,xhi,yhi
      DATA ndim,wsel /2,4/! 2D as default for x-y images
!      INTEGER indx
!      COMMON /pdebug/ indx
      CALL wscale
      rhmaxx=0.5d0*DBLE(nxymax)*dx
      rhmaxy=0.5d0*DBLE(nxymax)*dy
!      rh=2.d0*h
      rhx=min(hrad(wsel)*hx,rhmaxx)!avoid too large radii
      rhy=min(hrad(wsel)*hy,rhmaxy)
!      hu=0.5d0*rh
      hu=rhx/hrad(wsel)!use trimmed x radius as reference scale
!      write(*,*) 'hu =',hu
!      stop
      xmin=x-rhx
      xmax=x+rhx
      ymin=y-rhy
      ymax=y+rhy
!      write(*,*) 'rh =',rh
!      write(*,*) 'rh,x,xmin,xmax =',rh,x,xmin,xmax
!      write(*,*) 'rh,y,ymin,ymax =',rh,y,ymin,ymax
      xlo=xref-dx*int((xref-xmin)/dx)!round towards xref by dx
      xhi=xref+dx*int((xmax-xref)/dx)!round towards xref by dx
      ylo=yref-dy*int((yref-ymin)/dy)!round towards yref by dy
      yhi=yref+dy*int((ymax-yref)/dy)!round towards yref by dy
      nx=min(nint((xhi-xlo)/dx)+1,nxymax)!another overflow trap
      ny=min(nint((yhi-ylo)/dy)+1,nxymax)
      snorm=0.d0
      p=0
      pin=0
!      IF (min(nx*dx,ny*dy).lt.2.d0*rh) THEN
!         WRITE(*,*) 'Matrix too small for kernel'
!         WRITE(*,*) 'nx*dx,ny*dy,2*rh =',nx*dx,ny*dy,2.d0*rh
!         WRITE(*,*) 'nx,ny,dx,dy =',nx,ny,dx,dy
!         WRITE(*,*) 'xhi,yhi =',xhi,yhi
!         WRITE(*,*) 'xlo,ylo =',xlo,ylo
!         WRITE(*,*) 'Dx,Dy   =',xhi-xlo,yhi-ylo
!         WRITE(*,*) 'nxymax =',nxymax
!         STOP
!      ENDIF
      DO i=1,nx
         px=xlo+(i-1)*dx
         DO j=1,ny
            py=ylo+(j-1)*dy
            p=p+1
            pmxy(1,p)=px
            pmxy(2,p)=py
!            rxy=sqrt((px-xref)**2+(py-yref)**2)!why did I use the bin center instead of kernel center here?
!            rxy=sqrt((px-x)**2+(py-y)**2)      !use true kernel center rather than closest bin center for kernel function(?)
            rxy=sqrt((px-x)**2+((py-y)*hx/hy)**2)!scale y position to hu=trimmed hx
            mxy=wkernel(rxy,hu,ndim,wsel)
!            write(60,*) px,py
            IF (rxy.le.rhx) THEN
!            IF (mxy.gt.0.d0) THEN
               pin=pin+1
!               write(66,*) px,py
            ENDIF
            pmxy(0,p)=mxy
            snorm=snorm+mxy
!            if(indx.ge.988) write(*,*) 'indx,i,snorm =',indx,i,snorm
         ENDDO
      ENDDO
      np=p
!      write(*,*) 'nx,ny,np,nx*ny =',nx,ny,np,nx*ny
      npxy(0)=np
      npxy(1)=nx
      npxy(2)=ny
      npxy(3)=pin
!      unorm=1.d0/dble(pin)
!     Normalisation
      IF (snorm.gt.0.d0) THEN
         DO i=1,np
            pmxy(0,i)=m*pmxy(0,i)/snorm
         ENDDO
      ENDIF
      RETURN
      END

! Subroutine version of wkernel.f
! Role of function wkernel should be the same, but now
! it wraps the routine ckernel (calculate kernel)

      SUBROUTINE ckernel(r,h,ndim,wsel, kio,w)
!     Kernel functions, including those by
!     Monaghan & Lattanzio (A&A 1985)
!     (set wsel=4 for the cubic spline)
!
!     INPUT:
!     -- r   : distance from kernel center
!     -- h   : scaling radius of kernel (total radius = 2h)
!     -- ndim: dim switch: dim = 1,2,3
!     -- wsel: kernel selector:
!              0: simple hard switch
!              1: linear
!              2: Gaussian (Gingold & Monaghan 1977)
!              3: M_3 spline (Monaghan & Lattanzio 1985)
!              4: M_4 spline (           "             )
!              5: M_5 spline (UNDER CONSTRUCTION)
!              6: M_6 spline (UNDER CONSTRUCTION)
!              7: Lucy (1977) just for completeness...
!
!     OUTPUT (subroutine only):
!     -- kio=k : volume scaling factor
!     -- w     : smoothing weight (kernel function = w*k)
!
!     Effective kernel volume V = 1/k
!         |n| | V_wml        | V_simple     | V_wml/V_simple
!         ----+--------------+--------------+---------------
!          1  | 3/2    * h   | 2      * h   | 3/4  = 0.75
!          2  | 7/10*pi* h^2 | pi     * h^2 | 7/10 = 0.70
!          3  | pi     * h^3 | 4/3*pi * h^3 | 3/4  = 0.75
!     
!     Relationship between wml85 and the Gaussian function:
!     h = 2/3*sqrt(2*pi) => wml85(0,h,1) = fgauss(0,u=0,sigma=1)
!     As a rule of the thumb: sigma=0.6h yields a near-fitting Gaussian.
      IMPLICIT NONE
      INTEGER ndim,wsel,ndold,wsold
      DOUBLE PRECISION r,h,hold,xi,w,k,kio
!--- Kernel switch
!--- To be improved: wmlmod=1 => M_4
!                             =2 => M_3
!                            n<0 => switch (es before)
!---- scaling constants
      DOUBLE PRECISION aconst(1:3,0:7)
      COMMON /wconst/ aconst
      SAVE ndold,wsold,hold,k
      DATA ndold,wsold,hold /-1,-1,1.d0/
      CALL wscale
      IF (ndim.lt.1.OR.ndim.gt.3) THEN
         WRITE(*,*) 'Bad dimension (ndim must be within 1 and 3).'
         STOP
      ENDIF
      IF (wsel.lt.0.OR.wsel.gt.7) THEN
         WRITE(*,*) 'Bad selector (wsel must be within 0 and 4).'
         STOP
      ENDIF
!---- Calculate k only when it is needed
      IF (ndim.ne.ndold.OR.wsel.ne.wsold.OR.h.ne.hold) THEN
         ndold=ndim
         wsold=wsel
         hold=h
         k=aconst(ndim,wsel)/h**ndim
      ENDIF
!      k=1.d0/h**ndim
!===== COMPUTE WEIGHTING FOR GIVEN TEST PARTICLE
!      xi=abs(r/h)
      xi=r/h
      IF (wsel.eq.0) THEN
!---- Hard switch!Changed radius from 1 to 1/2
         IF (xi.lt.0.5d0) THEN
            w = 1.d0
         ELSE
            w = 0.d0
         ENDIF
      ELSE IF (wsel.eq.1) THEN
!---- Linear interpolation!Changed radius from 2 to 1
         IF (xi.lt.1.d0) THEN
            w = 1.d0-xi
         ELSE
            w = 0.d0
         ENDIF
      ELSE IF (wsel.eq.2) THEN
!---- Gingold & Monaghan (1977): Gauss-like
         w = exp(-xi**2)
      ELSE IF (wsel.eq.3) THEN
!---- Monaghan & Lattanzio (1985): M_3
         IF (xi.lt.0.5d0) THEN
!            w = (1.5d0-xi)**2 - 3.d0*(0.5d0-xi)**2
!....equivalent, but more efficient notation
            w = 1.5d0-2.d0*xi**2
         ELSE IF (xi.ge.0.5d0 .AND. xi.le.1.5d0) THEN
            w = (1.5d0-xi)**2
         ELSE
            w = 0.d0
         ENDIF
      ELSE IF (wsel.eq.4) THEN
!---- Monaghan & Lattanzio (1985): M_4
         IF (xi.lt.1.d0) THEN
!            w = 0.25d0*(2.d0-xi)**3-(1.d0-xi)**3
!....equivalent, but more efficient notation
            w = 1.d0-1.5d0*xi**2+0.75d0*xi**3
!....exfactorisation of 1/4 (equivalent, but looks nicer in math book)
!            w = 4.d0-6d0*xi**2+3.d0*xi**3
         ELSE IF (xi.ge.1.d0 .AND. xi.le.2.d0) THEN
            w = 0.25d0*(2.d0-xi)**3
!....exfactorisation of 1/4
!            w = (2.d0-xi)**3
         ELSE
            w = 0.d0
         ENDIF
      ELSE IF (wsel.eq.5) THEN
!---- M_5
         IF (xi.le.0.5d0) THEN
            w = (2.5d0-xi)**4 - 5.d0*(1.5d0-xi)**4 + 1.d1*(0.5d0-xi)**4
         ELSE IF (xi.gt.0.5d0 .AND. xi.le.1.5d0) THEN
            w = (2.5d0-xi)**4 - 5.d0*(1.5d0-xi)**4
         ELSE IF (xi.gt.1.5d0 .AND. xi.le.2.5d0) THEN
            w = (2.5d0-xi)**4
         ELSE
            w = 0.d0
         ENDIF
      ELSE IF (wsel.eq.6) THEN
!---- M_6
         IF (xi.le.1.d0) THEN
            w = (3.d0-xi)**5 - 6.d0*(2.d0-xi)**5 + 15.d0*(1.d0-xi)**5
         ELSE IF (xi.gt.1.d0 .AND. xi.le.2.d0) THEN
            w = (3.d0-xi)**5 - 6.d0*(2.d0-xi)**5
         ELSE IF (xi.gt.2.d0 .AND. xi.le.3.d0) THEN
            w = (3.d0-xi)**5
         ELSE
            w = 0.d0
         ENDIF
      ELSE IF (wsel.eq.7) THEN
!---- Lucy (1977)
         IF (xi.lt.1.d0) THEN
            w = (1.d0+3.d0*xi)*(1.d0-xi)**3
         ELSE
            w = 0.d0
         ENDIF
      ENDIF
      kio=k
!      write(*,*) 'h,ndim,wsel,aconst =',h,ndim,wsel,aconst(ndim,wsel)
!      write(*,*) 'ndim,wsel,k =',ndim,wsel,k
!      write(*,*) 'r,h,w,k =',r,h,w,k
      RETURN
      END

      SUBROUTINE wscale
!==== Set all the scaling constants
!     aconst(ndim,wsel): scaling constant for kernel value where
!                        ndim is the dimension and wsel the kernel ID
!     rsigma(hmod,wsel): scaling factors for h,
!                        hmod=0: no scaling (dummy)
!                             1: max = Gauss peak
!                             2: integral over +-h matches 1 sigma area
!                                (for 1D only)
!     hrad(wsel)       : kernel radius in units of h
!                        Gaussian: 3 sigma
!     rsigma is not used by wkernel but can be called externally.
!     sigma_68.5% = rsigma(mode,wsel)!changed!!!
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION aconst(1:3,0:7)
      DOUBLE PRECISION pi,gpeak,ss1sig
!==== MATHEMATICAL CONSTANTS
      PARAMETER(pi=3.141592653589793238d0,                              &
     &       gpeak=3.989422804014326779d-1,                             &
     &      ss1sig=0.68268949213708d0)
      COMMON /wconst/ aconst
!.... rsigma is in COMMON block for export purposes
      DOUBLE PRECISION rsigma(-1:3,0:7),hrad(0:7)
      COMMON /waux/ rsigma,hrad
!.... First call: first = .true.
      LOGICAL first
      SAVE first
      DATA first /.true./
      DATA rsigma /40*1.d0/
      IF (.not.first) RETURN
!==== SET KERNEL SCALING FACTORS
!---- Hard switch
      aconst(1,0)=1.d0     !dim=1
      aconst(2,0)=4.d0/pi  !dim=2
      aconst(3,0)=6.d0/pi  !dim=3
!---- Linear
      aconst(1,1)=1.d0
      aconst(2,1)=3.d0/pi
      aconst(3,1)=3.d0/pi
!---- Gaussian (Gingold & Monaghan 1977)
      aconst(1,2)=1.d0/sqrt(pi)
      aconst(2,2)=1.d0/pi
      aconst(3,2)=1.d0/sqrt(pi)**3
!---- WML85 M_3 spline
      aconst(1,3)=1.d0/2.d0
      aconst(2,3)=16.d0/13.d0/pi!bugfixed!
      aconst(3,3)=1.d0/pi
!---- WML85 M_4 spline
      aconst(1,4)=2.d0/3.d0
      aconst(2,4)=10.d0/7.d0/pi
      aconst(3,4)=1.d0/pi
!.... for exfactorisation of 1/4 (equivalent but looks nicer in math book)
!      aconst(1,4)=1.d0/6.d0
!      aconst(2,4)=5.d0/14.d0/pi
!      aconst(3,4)=1.d0/4.d0/pi
!---- M_5 spline
      aconst(1,5)=1.d0/24.d0
      aconst(2,5)=96.d0/1199.d0/pi
      aconst(3,5)=1.d0/20.d0/pi
!---- M_6 spline
      aconst(1,6)=1.d0/120.d0
      aconst(2,6)=7.d0/478.d0/pi
      aconst(3,6)=1.d0/120.d0/pi
!---- Lucy (1977)
      aconst(1,7)=5.d0/4.d0
      aconst(2,7)=5.d0/pi
      aconst(3,7)=105.d0/16.d0/pi
!==== SET H SCALING FACTORS (external only)
      rsigma(-1,0)=gpeak
      rsigma(1,0)=0.5d0*ss1sig
      rsigma(2,0)=0.4131251299960698d0
      rsigma(3,0)=0.4402618732382132d0
      rsigma(-1,1)=2.d0*gpeak
      rsigma(1,1)=1.d0-sqrt(1.d0-ss1sig)!2.d0*(1.d0-sqrt(1.d0-ss1sig))
      rsigma(2,1)=0.6243571921829532d0!1.2487143843659063d0
      rsigma(3,1)=0.7176558542175803d0!1.4353117084351605d0
      rsigma(-1,2)=sqrt(0.5d0)
      rsigma(1,2)=sqrt(0.5d0)
      rsigma(2,2)=1.0713890350611710d0
      rsigma(3,2)=1.3279194968562078d0
!      rsigma(2,2)=sqrt(0.5d0)!unfixed scaling problem, simple circumvention
!      rsigma(3,2)=sqrt(0.5d0)
      rsigma(-1,3)=gpeak/0.75d0
      rsigma(1,3)=0.5162866400795781d0
      rsigma(2,3)=0.7569743606039291d0
      rsigma(3,3)=0.9035198713931222d0
      rsigma(-1,4)=3.d0/2.d0*gpeak
      rsigma(1,4)=0.5931516369275666d0
      rsigma(2,4)=0.8718129610194552d0
      rsigma(3,4)=1.0507263852356394d0
      rsigma(-1,5)=192.d0/115.d0*gpeak
      rsigma(1,5)=0.6591294551795218d0
      rsigma(2,5)=0.9757593665238691d0
      rsigma(3,5)=1.1816797769119092d0
      rsigma(-1,6)=gpeak/0.55d0
      rsigma(1,6)=0.7194534143822457d0
      rsigma(2,6)=1.0691541338904875d0
      rsigma(3,6)=1.3001513895068315d0
      rsigma(-1,7)=5.d0/4.d0/gpeak
      rsigma(1,7)=0.3194382856229125d0
      rsigma(2,7)=0.4682752067368575d0
      rsigma(3,7)=0.5577272142501907d0
!==== SET MAX OR RECOMMENDED KERNEL RADIUS
      hrad(0)=0.5d0
      hrad(1)=1.d0
      hrad(2)=3.d0
      hrad(3)=1.5d0
      hrad(4)=2.d0
      hrad(5)=2.5d0
      hrad(6)=3.d0
      hrad(7)=1.d0
!==== FINISH
      first=.false.
      RETURN
      END

      SUBROUTINE showhscale
      IMPLICIT NONE
      INTEGER i,j
      DOUBLE PRECISION rsigma(-1:3,0:7),hrad(0:7)
      COMMON /waux/ rsigma,hrad
      CALL wscale
      WRITE(*,*)
      WRITE(*,*) 'wsel=0 : hard switch'
      WRITE(*,*) 'wsel=1 : linear interpolation'
      WRITE(*,*) 'wsel=2 : Gauss-type (times constant)'
      WRITE(*,*) 'wsel=3 : M_3 spline ML85'
      WRITE(*,*) 'wsel=4 : M_4 spline ML85'
      WRITE(*,*) 'wsel=5 : Lucy (1977)'
      WRITE(*,*)
      WRITE(*,*) 'hmode=0 : no rescaling'
      WRITE(*,*) 'hmode=-1: same peak height as Gauss standard'
      WRITE(*,*) 'hmode=d : same 1-sigma integral in d dimensions'
      WRITE(*,*)
      WRITE(*,*) 'rsigma: =sigma/h; scale h by rsigma to get sigma.'
      WRITE(*,*) 'hrad  : =R/h; Kernel radius R in units of h'
      WRITE(*,*)
      WRITE(*,*) 'Table 1: 1/rsigma(i,j) = 1 sigma; kernel radii R/h:'
      WRITE(*,10) 'i=wsel\j=hmode',(j,j=-1,3),'R/h'
      DO j=0,7
         WRITE(*,11) j,(1.d0/rsigma(i,j),i=-1,3),hrad(j)
      ENDDO
      WRITE(*,*)
      WRITE(*,*) 'Table 2: h in units of sigma:'
      WRITE(*,10) 'i=wsel\j=hmode',(j,j=-1,3)
      DO j=0,7
         WRITE(*,11) j,(rsigma(i,j),i=-1,3)
      ENDDO
      WRITE(*,*)
      WRITE(*,*) 'Table 3: R in units of sigma:'
      WRITE(*,10) 'i=wsel\j=hmode',(j,j=-1,3)
      DO j=0,7
         WRITE(*,11) j,(hrad(j)*rsigma(i,j),i=-1,3)
      ENDDO
      WRITE(*,*)
 10   FORMAT(A14,5(4X,I10),6X,A)
 11   FORMAT(I14,5(4X,F10.8),6X,F3.1)
      RETURN
      END

      DOUBLE PRECISION FUNCTION wkernel(r,hi,ndim,wsel)
      IMPLICIT NONE
      INTEGER ndim,wsel
      DOUBLE PRECISION r,hi,w,k,h
      h=abs(hi)
      CALL ckernel(r,h,ndim,wsel, k,w)
      IF (hi.lt.0.d0) THEN
         wkernel=w              !plain weight
      ELSE
         wkernel=k*w            !scaled for integral=1
      ENDIF
      RETURN
      END
