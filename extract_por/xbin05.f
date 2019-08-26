      SUBROUTINE xbin(imod,bmod,x,x0,xn,nb,kw,ybin,iybin,ib)
!******** BINNING ALGORITHM v 0.50 ****
!     XBIN Calculates the number of events within interval partition
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
!     imod = setup option: 0 - setup bin positions ybin(1,*) (first each)
!                              (will then be set to 1 automatically)
!                         >0 - re-use from given ybin
!                          2 - same as 1 + cumulative bins
!                         <0 - just get bin index (and setup)
!     bmod = boundary option: 0 - strictly excludes right interval sup
!                             1 - counts event if equal to sup of nth bin
!     INPUT:
!     x    = event abscissa
!     x0,xn= overall interval (from left side of first bin to right side of last bin)
!     nb   = number of partitions (set negative for re-initialization)
!     kw   = weighting factor (normally equal 1)
!
!     OUTPUT:
!     ybin = output array
!            (1,i) = ith midpoint abscissa, i=1...nb
!            (2,i) = ith bin value
!            (3,i) = spare array (used for cumulative binning in imod=2)
!     iybin = unweighted integer binning array
!             bins are advanced by 1 if kw != 0, else unchanged
!     ib    = index of actual bin
!             Boundary: ib=0 for x<u_0, ib=nb+1 for x>(=) u_nb
!
!     AUXILIARY VARIABLES/SWITCHES:
!     Copy corresponding COMMON block into calling program
!     dokb  = ycnt weighting switch    (input)
!                      true: kb=kw*db
!                     false: kb=kw
!     ycnt  = (0) weighted (by kw AND bin width) total number < x0
!             (1)    "       "     "    within range
!             (2)    "       "     "    >(=) xn
!     iycnt = unweighted count of inside/outside hits (if kw!=0)
!     nkweq0= Number if calls with kw=0 (inside + outside bins)
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER nb,imod,bmod,iybin(nb),ib,jb!,ncall
      DOUBLE PRECISION x,u,x0,xn,db,kw,ybin(3,nb),kb
!      SAVE db
      LOGICAL add2y
!.... Auxiliary variables: Copy this block into calling program to use
      LOGICAL dokb!used for ycnt weighting only
      INTEGER iycnt(0:2),nkweq0
      DOUBLE PRECISION ycnt(0:2)
      COMMON /xbaux/ ycnt,iycnt,nkweq0,dokb
!.... Debugging switch
      LOGICAL dotnk
      COMMON /xbntnk/ dotnk
!      SAVE ncall
!      DATA ncall/0/
!      ncall=ncall+1
      IF (x0.ge.xn) THEN
         WRITE(*,*) 'XBIN ERROR: x0 < xn violated! --> Abort'
         WRITE(*,*) 'x0 = ',x0,'          xn = ',xn
         STOP
      ENDIF
!---- Setup
      db=(xn-x0)/nb!calculate this always (e.g. for ybin's of different size)
      IF (imod.eq.0) THEN
         CALL setxb(x0,xn,nb,db, ybin,ycnt,iybin,iycnt, dotnk)
         nkweq0=0
         imod=1!unexpected reset may be fatal
      ENDIF
!---- Switches
      add2y=kw.ne.0.d0
      IF (add2y) nkweq0=nkweq0+1
      IF (dokb) THEN
         kb=kw*db
      ELSE
         kb=kw
      ENDIF
!---- Daily use
      u=x-x0
      ib=max(min(int(u/db+1.d0),nb+1),0)
      IF (bmod.eq.1 .and. x.eq.xn) ib=nb
!.... Return if index-only imod is set
      IF (imod.lt.0) RETURN
!.... Continue else 
      jb=max(min(ib,nb),1)
!      if(x0.eq.0.d0) then
      if (dotnk) then
         write(*,1000) 'x,x0,xn,yb,db,nb,ib =',                         &
     &        x,x0,xn,ybin(1,jb),db,nb,ib
!         if(ib.eq.0) pause
!         if(ib.le.1) write(*,*) 'ib,x = ',ib,x
      endif
      IF (ib.le.0) THEN
         ycnt(0)=ycnt(0)+kb
         IF (add2y) iycnt(0)=iycnt(0)+1
      ELSE IF (ib.gt.nb) THEN
         ycnt(2)=ycnt(2)+kb
         IF (add2y) iycnt(2)=iycnt(2)+1
      ELSE
         ycnt(1)=ycnt(1)+kb
         IF (add2y) iycnt(1)=iycnt(1)+1
         ybin(2,ib)=ybin(2,ib)+kw
         IF (add2y) iybin(ib)=iybin(ib)+1
      ENDIF
!      write(*,*) 'ib,nb,iycnt =',ib,nb,iycnt(0),iycnt(1),iycnt(2),add2y
!---- Create cumulative histogram (e.g. for final call of xbin)
      IF (imod.eq.2) THEN
         CALL cumbin(nb,ybin)
      ENDIF
 1000 FORMAT(A,1X,5(F9.6,2X),2(2X,I6))
      RETURN
      END

      SUBROUTINE setxb(x0,xn,nb,db, ybin,ycnt,iybin,iycnt, dotnk)
      IMPLICIT NONE
      INTEGER nb,i
      INTEGER iybin(nb),iycnt(0:2)
      DOUBLE PRECISION x0,xn,db,ybin(3,nb),ycnt(0:2)
      LOGICAL dotnk
!      write(*,*) 'XBIN: Initializing...'
      DO 1 i=1,nb
!--- Calculation of bin midpoints
!     ybin(1,i) abscissa equal to bin midpoint
!     ybin(2,*) use small non-zero value to avoid log crashing
         ybin(1,i)=x0+(dble(i)-0.5d0)*db
         ybin(2,i)=1.d-99
         ybin(3,i)=1.d-99
         iybin(i)=0
         if (dotnk) then
!         if (xn-x0.le.1.01) then
!            write(*,*) 'i,ybin012 = ',i,ybin(1,i)-db/2,                 &
!     &                                  ybin(1,i),                      &
!     &                                  ybin(1,i)+db/2
!         endif
            write(*,*) 'i,xi =',i,ybin(1,i)
!         if(ib.eq.0) pause
!         if(ib.le.1) write(*,*) 'ib,x = ',ib,x
         endif
 1    CONTINUE
      ycnt(0)=0.d0
      ycnt(1)=0.d0
      ycnt(2)=0.d0
      iycnt(0)=0
      iycnt(1)=0
      iycnt(2)=0
      RETURN
      END

      SUBROUTINE cumbin(nb,ybin)
      IMPLICIT NONE
      INTEGER nb,i
      DOUBLE PRECISION ybin(3,nb)
      ybin(3,1)=ybin(2,1)
      DO 5 i=2,nb
         ybin(3,i)=ybin(3,i-1)+ybin(2,i)
 5    CONTINUE
      RETURN
      END
!**** Change Log v0.50
!     Comment character changed from 'c' to '!' for F90 compatibility
!     Changes since v0.41:
!     - Changed transfer list;
!       added iycnt to it and moved ycnt to COMMON block
!     - transfer dokb via COMMON block (since only needed for ycnt)
!     Changes since v0.402:
!     - imod is now set to 1 automatically after setup
!       (to prevent unexpected overwriting of arrays)
!     - mode renamed into imod ("mode" is a too common name)
!     Changes during v0.401:
!     - (comments only) add description for cumulative bind (mode=2)
!     Changes since v0.40:
!     - Extended "dotnk" debug ouputs
!     Changes since v0.32:
!     - Extended histogram array for cumulative column ybin(3,*) for mode=2
!     - Restricted setup to mode=0 (previously mode<=0)
!     - Now located in Tools folder (formerly IMF)
!     Changes since v0.31:
!     - Optional db-weighing of ynct
!     - lower initial ybin(2,*); now 10^-99 instead of 10^-30
!     - more comments...

