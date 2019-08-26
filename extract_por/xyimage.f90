module xyimage
  use commons4ramses
  implicit none
! image dimension
  INTEGER,PARAMETER :: nxybmax=1200
  integer :: nbx_xyimage,nby_xyimage
  real(kind=8) :: area_pixel,dx_pixel,dy_pixel
  real(kind=8) :: cx_xyimage=0.d0,cy_xyimage=0.d0
  real(kind=8) :: wx_xyimage=100.d0,wy_xyimage=100.d0
  real(kind=8) :: hx_xyimage=1.d0,hy_xyimage=1.d0!smoothing length xy
  real(kind=8) :: hxpix_xyimage=2.d0,hypix_xyimage=2.d0!smoothing length xypixel
  real(kind=8) :: rgb_scallevel=0.05d0,rgb_maxlevel=255.d0,rgb_xgamma=1.d0,rgb_scallog=1.d0
  real(kind=8),dimension(nxybmax) :: xtab,ytab
  real(kind=8),dimension(nxybmax,nxybmax) :: qbin,qtbin,qbin_r,qbin_g,qbin_b
contains

  subroutine get_xycoords(i,edgemode,x,y)
! Subroutine added at 2016-06-20
    integer :: i,edgemode
    real(kind=8) :: x,y
    real(kind=8),dimension(3) :: pos
    pos(1)=all_data(1,i)
    pos(2)=all_data(2,i)
    pos(3)=all_data(3,i)
! changed 2015-12-15: add option for y-z projection
       if (edgemode==1) then
!    x-z
          x=pos(1)
          y=pos(3)
       else if (edgemode>=2) then
!    y-z
          x=pos(2)
          y=pos(3)
       else
!    x-y
          x=pos(1)
          y=pos(2)
       endif
    return
  end subroutine get_xycoords
  
  subroutine get_xyimage(edgemode,massbinmode,nskip_particles)
    integer :: edgemode,massbinmode,mode,bmod=0,i,j,nskip_particles
    real(kind=8) :: x0,xn,y0,yn!xy range
    real(kind=8) :: kw=1.d0,wa=1.d0,tage,kwa,kha=1.d0
    real(kind=8),dimension(3) :: pos
    real(kind=8) :: x,y
    INTEGER ibq,isinq
    DOUBLE PRECISION qcnt(0:1),qref(2)
    COMMON /xybaux/ qcnt,qref,ibq,isinq
    qcnt=0.d0
    x0=cx_xyimage-wx_xyimage/2.d0
    xn=x0+wx_xyimage
    y0=cy_xyimage-wy_xyimage/2.d0
    yn=y0+wy_xyimage
    hx_xyimage=hxpix_xyimage*wx_xyimage/nbx_xyimage
    hy_xyimage=hypix_xyimage*wy_xyimage/nby_xyimage
    dx_pixel=wx_xyimage/nbx_xyimage
    dy_pixel=wy_xyimage/nby_xyimage
    area_pixel=dx_pixel*dy_pixel
    mode=0
    do i=1,ndata
       if (i<=nskip_particles) cycle
!       if (i==10000*int(i/10000)) print *,'i =',i,' of ',ndata
!       pos(1)=all_data(1,i)
!       pos(2)=all_data(2,i)
!       pos(3)=all_data(3,i)
!! changed 2015-12-15: add option for y-z projection
!       if (edgemode==1) then
!!    x-z
!          x=pos(1)
!          y=pos(3)
!       else if (edgemode>=2) then
!!    y-z
!          x=pos(2)
!          y=pos(3)
!       else
!!    x-y
!          x=pos(1)
!          y=pos(2)
!       endif
       call get_xycoords(i,edgemode,x,y)
       tage=tnow_myr-all_data(11,i)
       tage=tnow_myr-all_data(11,i)
       if (age_exp>0.and.age_thres/=0.d0) then
          wa=(tage/abs(age_thres))**age_exp
       else
          wa=1.d0
       endif
       if (flg_agespread>=1) then
          kha=1.d0/sqrt(wa+agespread_heps)
       else
          kha=1.d0
       endif
       if (massbinmode>=1) then
          kw=all_data(7,i)/area_pixel
       else
          kw=1.d0
       endif
       if ((age_thres>0.d0).and.(tage>age_thres.or.all_data(11,i)<=0.d0)) then
! age_mode only: skip dm particles and stars older than age_threshold
          kw=0.d0
       endif
       if (massbinmode==2.and.all_data(11,i)<0.d0) then
          kw=0.d0
       endif
!       kw=1.d0!test
!       if (i<=1) print *,'hx_xyimage,hy_xyimage =',hx_xyimage,hy_xyimage,massbinmode
       call wrap_xysmooth(mode,bmod,x,y,x0,xn,y0,yn,nbx_xyimage,nby_xyimage,&
            hx_xyimage*kha,hy_xyimage*kha,kw, xtab,ytab,qbin)
       mode=1
    enddo
    if (massbinmode==2) then
!       print *,'ndata,nbx_xyimage,nby_xyimage =',ndata,nbx_xyimage,nby_xyimage
       mode=0
       qcnt=0.d0
       do i=1,ndata
          if (i<=nskip_particles) cycle
          call get_xycoords(i,edgemode,x,y)
          tage=tnow_myr-all_data(11,i)
          if (age_exp>0.and.age_thres/=0.d0) then
             wa=(tage/abs(age_thres))**age_exp
          else
             wa=1.d0
          endif
          if (flg_agespread>=1) then
             kha=1.d0/sqrt(wa+agespread_heps)
          else
             kha=1.d0
          endif
          if (flg_agespread==1) then
             kwa=wa
          else
             kwa=1.d0
          endif
          if (all_data(11,i)>=0.d0) then
             kw=tage*kwa
          else
             kw=0.d0
          endif
          kw=kw*all_data(7,i)/area_pixel
!          if (kw>0.d0) print '(a,3(e14.4))','kw =',tnow_myr,all_data(11,i),kw
          call wrap_xysmooth(mode,bmod,x,y,x0,xn,y0,yn,nbx_xyimage,nby_xyimage,&
               hx_xyimage*kha,hy_xyimage*kha,kw, xtab,ytab,qtbin)
          mode=1
       enddo
       do j=1,nby_xyimage
          do i=1,nbx_xyimage
             if (qtbin(i,j)>0.d0) then
                qtbin(i,j)=qtbin(i,j)/max(qbin(i,j),1.d-30)
             else
                qtbin(i,j)=1.d12
             endif
          enddo
       enddo
    endif
    return
  end subroutine get_xyimage

  subroutine get_rgb_xyimage(edgemode,massbinmode,nskip_particles)
    integer :: edgemode,massbinmode,mode_r,mode_g,mode_b,bmod=0,i,nskip_particles
    real(kind=8) :: x0,xn,y0,yn!xy range
    real(kind=8) :: kw=1.d0
    real(kind=8),dimension(3) :: pos
    real(kind=8) :: x,y
    x0=cx_xyimage-wx_xyimage/2.d0
    xn=x0+wx_xyimage
    y0=cy_xyimage-wy_xyimage/2.d0
    yn=y0+wy_xyimage
    hx_xyimage=hxpix_xyimage*wx_xyimage/nbx_xyimage
    hy_xyimage=hypix_xyimage*wy_xyimage/nby_xyimage
    dx_pixel=wx_xyimage/nbx_xyimage
    dy_pixel=wy_xyimage/nby_xyimage
    area_pixel=dx_pixel*dy_pixel
    mode_r=0;mode_g=0;mode_b=0
    do i=1,ndata
       if (i<=nskip_particles) cycle
       pos(1)=all_data(1,i)
       pos(2)=all_data(2,i)
       pos(3)=all_data(3,i)
!       x=pos(1)
!       if (edgemode>=1) then
!          y=pos(3)
!       else
!          y=pos(2)
!       endif
! changed 2015-12-15: add option for y-z projection
       if (edgemode==1) then
!    x-z
          x=pos(1)
          y=pos(3)
       else if (edgemode>=2) then
!    y-z
          x=pos(2)
          y=pos(3)
       else
!    x-y
          x=pos(1)
          y=pos(2)
       endif
       if (massbinmode==1) then
          kw=all_data(7,i)/area_pixel
       else
          kw=1.d0
       endif
!       if (i<=1) print *,'hx_xyimage,hy_xyimage =',hx_xyimage,hy_xyimage,massbinmode
       if (i>=partition_commons(1,1).and.i<=partition_commons(2,1)) then
! xtab, ytab identical for all colors
          call wrap_xysmooth(mode_r,bmod,x,y,x0,xn,y0,yn,nbx_xyimage,nby_xyimage,&
               hx_xyimage,hy_xyimage,kw, xtab,ytab,qbin_r)
          mode_r=1
       else if (i>=partition_commons(1,2).and.i<=partition_commons(2,2)) then
          call wrap_xysmooth(mode_g,bmod,x,y,x0,xn,y0,yn,nbx_xyimage,nby_xyimage,&
               hx_xyimage,hy_xyimage,kw, xtab,ytab,qbin_g)
          mode_g=1
       else if (i>=partition_commons(1,3).and.i<=partition_commons(2,3)) then
          call wrap_xysmooth(mode_b,bmod,x,y,x0,xn,y0,yn,nbx_xyimage,nby_xyimage,&
               hx_xyimage,hy_xyimage,kw, xtab,ytab,qbin_b)
          mode_b=1
       endif
    enddo
    return
  end subroutine get_rgb_xyimage

  subroutine write_xyimage(imagefile)
    character(len=80) :: imagefile
    integer :: i,j
    real(kind=8) :: x,y,f,avage
    open(20,file=trim(imagefile))
    do j=1,nby_xyimage
       y=ytab(j)
       do i=1,nbx_xyimage
          x=xtab(i)
          f=qbin(i,j)
          avage=qtbin(i,j)
!            f=log10(max(qbin(i,j),1.e-4))
          write(20,2000) x,y,f,avage
!          if (f/=1.d12.and.f/=0.d0) then
!          if (abs(f)<1.d3.and.abs(f)>0.d0) then
!             print *,'f=',f
!             stop
!          endif
       enddo
!         print *,'j = ',j
!         print *,'j,ytab(j)=',j,' ',ytab(j),' ',f
       write(20,*)
    enddo
    close(20)
2000 FORMAT(2(F12.6,2X),2(2X,1PE12.4))
    return
  end subroutine write_xyimage

  subroutine write_rgb_xyimage(imagefile,satmode)
    character(len=80) :: imagefile
    integer :: i,j,satmode,rgbint(3),srgbint(3)
    real(kind=8) :: x,y,f,rgblin(3),srgb(3),rgbsat(3),srgbsat(3),rgb255(3),srgb255(3),maxrgb
    maxrgb=0.d0
    open(20,file=trim(imagefile))
    if (rgb_xgamma<=0.d0) then
       rgb_scallog=1.d1**(-rgb_xgamma)
    endif
    do j=1,nby_xyimage
       y=ytab(j)
       do i=1,nbx_xyimage
          x=xtab(i)
          if (rgb_xgamma<=0.d0) then
             rgblin(1)=rgb_scallevel*log10(rgb_scallog*max(qbin_r(i,j),1.d-99))
             rgblin(2)=rgb_scallevel*log10(rgb_scallog*max(qbin_g(i,j),1.d-99))
             rgblin(3)=rgb_scallevel*log10(rgb_scallog*max(qbin_b(i,j),1.d-99))
          else
             rgblin(1)=rgb_scallevel*qbin_r(i,j)**rgb_xgamma
             rgblin(2)=rgb_scallevel*qbin_g(i,j)**rgb_xgamma
             rgblin(3)=rgb_scallevel*qbin_b(i,j)**rgb_xgamma
          endif
          f=sum(rgblin)
          call srgb_gamma(rgblin,srgb,1)
          maxrgb=max(maxrgb,maxval(rgblin))
          if (satmode>=1) then
             call oversaturate(rgblin,rgbsat)
             call oversaturate(srgb,srgbsat)
!             call srgb_gamma(rgbsat,srgbsat,1)
             rgb255(:)=rgb_maxlevel*max(min(rgbsat(:),1.d0),0.d0)
             srgb255(:)=rgb_maxlevel*max(min(srgbsat(:),1.d0),0.d0)
          else
             rgb255(:)=rgb_maxlevel*max(min(rgblin(:),1.d0),0.d0)
             srgb255(:)=rgb_maxlevel*max(min(srgb(:),1.d0),0.d0)
          endif
          rgbint(:)=nint(rgb255(:))
          srgbint(:)=nint(srgb255(:))
!          write(20,2000) x,y,f,rgbint,srgbint
          write(20,3000) x,y,f,rgb255,srgb255
!         f=log10(max(qbin(i,j),1.e-4))
       enddo
!       print *,'j = ',j
!       print *,'j,ytab(j)=',j,' ',ytab(j),' ',f
       write(20,*)
    enddo
    close(20)
    print *,'Maximum RGB value = ',rgb_maxlevel*maxrgb
!2000 FORMAT(2(F12.6,2X),7(2X,1PE12.4))
2000 FORMAT(2(F12.6,2X),2X,E12.4,2x,3(2x,I5),2x,3(2x,I5))
!2100 FORMAT(2(F12.6,2X),2X,E12.4,2x,3(2x,I5))
3000 FORMAT(2(F12.6,2X),2X,E12.4,2x,3(2x,F8.2),2x,3(2x,F8.2))
!3100 FORMAT(2(F12.6,2X),2X,E12.4,2x,3(2x,F8.2))
    return
  end subroutine write_rgb_xyimage

  subroutine oversaturate(rgbraw,rgbsat)
    integer :: i
    real(kind=8),dimension(3) :: rgbraw,rgbsat,rgbclip
    real(kind=8) :: s,s1,headroom(3)
    s=max(min(maxval(rgbraw),2.d0),0.d0)
    s1=max(s-1.d0,0.d0)
!    print *,'using oversaturate'
    do i=1,3
       rgbclip(i)=max(min(rgbraw(i),1.d0),0.d0)
       if (s>1.d0) then
          headroom(i)=1.d0-rgbclip(i)
          rgbsat(i)=rgbclip(i)+s1*headroom(i)
!          print *,'i,headroom(i),s1 =',i,headroom(i),s1
!          if (s1>0.2d0) print *,'i,s1,headroom,rgbsat,s1 =',i,s1,headroom(i),rgbsat(i)
       else
          rgbsat(i)=max(min(rgbraw(i),1.d0),0.d0)
       endif
    enddo
    return
  end subroutine oversaturate

  SUBROUTINE srgb_gamma(Jin,Jout,direction)
    IMPLICIT NONE
    INTEGER :: direction,i
    REAL(KIND=8) :: Jin(3),Jout(3)
    REAL(KIND=8) :: Jrgblin(3),Jsrgb(3),c,clin
    REAL(KIND=8) :: a=5.5d-2
    IF (direction>=1) THEN
! RGB_lin to RGB_gamma
       DO i=1,3
          clin=Jin(i)
          IF (clin<=0.0031308d0) THEN
             c=12.92d0*clin
          ELSE
             c=(1.d0+a)*clin**(1.d0/2.4d0)-a
          ENDIF
          Jout(i)=c
       ENDDO
    ELSE IF (direction<=-1) THEN
! RGB_gamma to RGB_lin
       DO i=1,3
          c=Jin(i)
          IF (c<=0.04045d0) THEN
             clin=c/12.92d0
          ELSE
             clin=((c+a)/(1.d0+a))**2.4d0
          ENDIF
          Jout(i)=clin
       ENDDO
    ELSE
! No gamma correction
       Jout=Jin
    ENDIF
    RETURN
  END SUBROUTINE srgb_gamma

end module xyimage
