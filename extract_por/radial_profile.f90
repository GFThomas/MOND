module radial_profile
  use commons4ramses
  implicit none
  real(kind=8) :: rbinmax,mass,dbindata,bincenter(3)
  integer,parameter :: nbinmax=10000,nlagrangemax=10
  integer :: nbindata,nlagrange
  integer,dimension(2,ndatamax) :: bindex
  real(kind=8),dimension(20,nbinmax) :: rbindata
  real(kind=8),dimension(nlagrangemax) :: aquantile,arlagrange

contains

  subroutine write_data(filename)
    integer :: ib
    character(len=80) :: filename
    open(45,file=trim(filename))
    write(45,4501) '1:rbin','2:mass/tot','3:mcum/tot','4:volshell','5:rhomass','6:sigma_r','7:sigma_t','8:sigma_3D','9:vcirc'
    do ib=1,nbindata
!       write(45,4500) rbindata(1,ib),nint(rbindata(2:3,ib)),rbindata(4:8,ib)
       if (useacc) then
          write(45,4500) rbindata(1,ib),rbindata(2:3,ib)/mass,rbindata(4:8,ib),rbindata(13,ib)
       else
          write(45,4500) rbindata(1,ib),rbindata(2:3,ib)/mass,rbindata(4:8,ib)
       endif
!       write(*,'(a,1p,2(1x,e12.4))') 'rbindata(2,ib),mass =',rbindata(2,ib),mass
    enddo
9   close(45)
!4500 format(  f16.4,2x,1p,8(2x,e16.8))
4500 format(  f16.4,2x,2(2x,f10.8),1p,7(2x,e16.8e3))
4501 format('#',a15,2x,2(2x,a10),     7(2x,a16))
  end subroutine write_data

  subroutine import_data(ndata_import,pdata_import,id_import)
    integer :: ndata_import,i,ii
    integer,dimension(ndatamax) :: id_import
    real(kind=8),dimension(11,ndatamax) :: pdata_import
    ndata=ndata_import
    do i=1,ndata
       ii=ai_sorted(i)
       all_data(1:11,i)=pdata_import(1:11,ii)
       bindex(1,i)=id_import(ii)
!       print *,'pos =',pdata_import(1:3,ii)
!       print *,'i,ii =',i,ii
!       if (i>=20001.and.i<=22124) print *,'pos =',all_data(1:3,i)
!       if (ii>=22125.and.ii<=23126) print *,'pos =',all_data(1:3,i)
!       if (i>=20001.and.i<=22124) print *,'pos =',pdata_import(1:3,i)
    enddo
    return
  end subroutine import_data

  subroutine radial_binning
    integer :: i,imod=0,bmod=0,nb,iybin(nbinmax),ib
    integer :: i1mod=0,iy1bin(nbinmax),i1b
    real(kind=8) :: r,r0,rn,kw,ybin(3,nbinmax),y1bin(3,nbinmax),db,ri,ro,pos(3),ddum
    r0=0.d0
    rn=rbinmax
    nb=nbindata; db=(rn-r0)/nb
    dbindata=db
!    mass=all_data(7,1)!all particles are equal - no longer true!
    mass=0.d0
!    if (boxlength>0) then
!       boxcenter=(/1.d0,1.d0,1.d0/)*boxlength/2.d0
!    endif
    call set_boxcenter
    bincenter=boxcenter
    imod=0
    i1mod=0
    do i=1,ndata
       mass=mass+all_data(7,i)!re-define 'mass' as 'total mass'; maybe useful some day
       kw=1.d0*all_data(7,i)!all masses are equal, but some are more more equal than others
       pos(1:3)=all_data(1:3,i)-bincenter(1:3)
       call vnsq3(pos,r,ddum)
       call xbin(imod,bmod,r,r0,rn,nb,kw,ybin,iybin,ib)
!       call xbin(i1mod,bmod,r,r0,rn,nb,1.d0,y1bin,iy1bin,i1b)
       bindex(1,i)=id_data(i)
       bindex(2,i)=ib!store index of bin in which the ith particle is in
       imod=2
       i1mod=2
    enddo
    do ib=1,nb
       rbindata(1:3,ib)=ybin(1:3,ib)
       ri=rbindata(1,ib)-dbindata/2.d0
       ro=ri+dbindata
       if (voldim==1) then
          rbindata(4,ib)=(ro-ri)
       else if (voldim==2) then
          rbindata(4,ib)=pi*(ro**2-ri**2)
       else if (voldim==3) then
          rbindata(4,ib)=4.d0/3.*pi*(ro**3-ri**3)
       endif
    enddo
    return
  end subroutine radial_binning

  subroutine radial_density
    integer :: i,ib
    real(kind=8) :: r0,rn,ri,ro,dvol,rho
!    r0=rbindata(1,0)-dbindata/2.d0
!    rn=rbindata(1,nbindata)+dbbindata/2.d0
    do ib=1,nbindata
       ri=rbindata(1,ib)-dbindata/2.d0
       ro=ri+dbindata
       dvol=rbindata(4,ib)
       rbindata(5,ib)=rbindata(2,ib)/dvol!mass factor already in radial_binning routine
    enddo
    return
  end subroutine radial_density

  subroutine radial_vdispersion
    implicit none
    integer :: i,ii,ib,nshell
    real(kind=8) :: r0,rn,ri,ro,dvol,rho,pos(3),vel(3),acc(3),vtsq,vrsq,atsq,arsq
    real(kind=8) :: ssrq,sstq,rvacom(9),ssmass,rsq,rabs,vsq,vabs,asq,aabs
!    do ib=1,nbindata
!       ri=rbindata(1,ib)-dbindata/2.d0
!       ro=ri+dbindata
!    enddo
    rbindata(6,:)=0.d0! reset 6th column dedicated for sigma_r
    rbindata(7,:)=0.d0! reset 7th column dedicated for sigma_t
    rbindata(8,:)=0.d0! reset 8th column dedicated for sigma
    rbindata(9,:)=0.d0! reset 9th column dedicated for asigma_r
    rbindata(10,:)=0.d0! reset 10th column dedicated for asigma_t
    rbindata(11,:)=0.d0! reset 11th column dedicated for asigma
!    rbindata(14,:)=0.d0! dummy
    ssrq=0.d0
    sstq=0.d0
    rvacom=0.d0
    ssmass=0.d0
    do i=1,ndata
       rvacom(1:3)=rvacom(1:3)+all_data(1:3,i)*all_data(7,i)
       rvacom(4:6)=rvacom(4:6)+all_data(4:6,i)*all_data(7,i)
       rvacom(7:9)=rvacom(7:9)+all_data(8:10,i)*all_data(7,i)
       ssmass=ssmass+all_data(7,i)
    enddo
    rvacom(1:9)=rvacom(1:9)/ssmass
    if (flgcom4bin<=0) then
       rvacom=0.d0!disable COM correction for binning
    endif
    do i=1,ndata
       ii=bindex(1,i)
       ib=bindex(2,i)
       pos(1:3)=all_data(1:3,i)-rvacom(1:3)
       vel(1:3)=all_data(4:6,i)-rvacom(4:6)
       acc(1:3)=all_data(8:10,i)-rvacom(7:9)
       call vrad3(vel,pos,vrsq,vtsq)
!       if (i<=100) then
!          rsq=pos(1)**2+pos(2)**2+pos(3)**2
!          vsq=vel(1)**2+vel(2)**2+vel(3)**2
!          print *,'rabs,vabs =',sqrt(rsq),sqrt(vsq),sqrt(vrsq+vtsq)
!       endif
       if (useacc) call vrad3(acc,pos,arsq,atsq)
       rbindata(6,ib)=rbindata(6,ib)+all_data(7,i)*vrsq
       rbindata(7,ib)=rbindata(7,ib)+all_data(7,i)*vtsq/2.d0!account for vt being 2-dim
       rbindata(8,ib)=rbindata(8,ib)+all_data(7,i)*(vrsq+vtsq)
!       print '(a,1p,2(1x,e12.4))','rbindata6,vrsq =',rbindata(6,ib),vrsq
       if (useacc) then
          rbindata(9,ib) =rbindata(9,ib) +all_data(7,i)*arsq
          rbindata(10,ib)=rbindata(10,ib)+all_data(7,i)*atsq/2.d0!account for vt being 2-dim
          rbindata(11,ib)=rbindata(11,ib)+all_data(7,i)*(arsq+atsq)
       endif
       call vnsq3(pos,rabs,rsq)
       call vnsq3(vel,vabs,vsq)
       call vnsq3(acc,aabs,asq)
!       rbindata(14,ib)=rbindata(14,ib)+all_data(7,i)*vsq
       rbindata(12,ib)=rabs!store abs radius
!       if (i==1) print '(a,9(1x,1pe12.4))','rabs,rbin3,rbin9,a =',&
!            rabs,rbindata(3,ib),rbindata(9,ib),aabs
    enddo
! convert square-v sums into velocity dispersion (=standard deviation of vr,vt)
    do ib=1,nbindata
       rbindata(6,ib)=sqrt(rbindata(6,ib)/rbindata(2,ib))
       rbindata(7,ib)=sqrt(rbindata(7,ib)/rbindata(2,ib))
       rbindata(8,ib)=sqrt(rbindata(8,ib)/rbindata(2,ib))
       if (useacc) then
          rbindata(9,ib)=sqrt(rbindata(9,ib)/rbindata(2,ib))
          rbindata(10,ib)=sqrt(rbindata(10,ib)/rbindata(2,ib))
          rbindata(11,ib)=sqrt(rbindata(11,ib)/rbindata(2,ib))
       endif
!       rbindata(14,ib)=sqrt(rbindata(14,ib)/rbindata(2,ib))
!       vsq=rbindata(9,ib)*rbindata(12,ib); vabs=sqrt(vsq)
       vsq=rbindata(9,ib)*rbindata(12,ib)/scalacc/scalr
!       vsq=rbindata(9,ib)*rbindata(12,ib)/scalr
       vabs=sqrt(vsq)*scalv
       rbindata(13,ib)=vabs
!       print *,'scalv =',scalv
!       if (ib<=100) then
!          print '(a,1p,9(1x,e12.4))','rbin3,9,12,vabs=',rbindata(3,ib),rbindata(9,ib),rbindata(12,ib),vabs,rbindata(14,ib)
!          print '(a,1p,9(1x,e12.4))','rbin,14=',rbindata(1,ib),rbindata(14,ib)
!       endif
    enddo
    return
  end subroutine radial_vdispersion

!  subroutine radial_rotcurve
!    integer :: i,ii,ib,nshell
!    real(kind=8) :: r0,rn,ri,ro,dvol,rho,pos(3),vel(3),vtsq,vrsq
!    real(kind=8) :: ssrq,sstq
!  end subroutine radial_rotcurve

  subroutine check_angmom
    integer :: i
    real(kind=8) :: lpart(3),ltot(3),pos(3),vel(3),absltot,abslpart,meanlpart,ddum
    ltot=0.d0
    meanlpart=0.d0
    do i=1,ndata
       pos(1:3)=all_data(1:3,i)
       vel(1:3)=all_data(4:6,i)
       call pvec3(pos,vel,lpart)
       call vnsq3(lpart,abslpart,ddum)
       meanlpart=meanlpart+abslpart/ndata
       ltot(1:3)=ltot(1:3)+lpart(1:3)
    enddo
    call vnsq3(ltot,absltot,ddum)
    print '(a,3(1x,1pe16.8))','L_tot        =',ltot
    print '(a,1x,1pe16.8)',   '||L||        =',absltot
    print '(a,1x,1pe16.8)',   '<||l||>      =',meanlpart
    print '(a,1x,1pe16.8)',   '||L||/Poisson=',absltot/(meanlpart*sqrt(dble(ndata)))
    return
  end subroutine check_angmom

  subroutine lagrange_radius(quantile,rlagrange)
    implicit none
    integer :: ib,flgfound
    real(kind=8) :: quantile,rlagrange
    flgfound=0
    do ib=1,nbindata
       if (flgfound==0.and.ib>1.and.rbindata(3,ib)>=quantile*rbindata(3,nbindata)) then
          if (rbindata(3,ib)==quantile*rbindata(3,nbindata)) then
             rlagrange=rbindata(1,ib)
          else
             rlagrange=0.5d0*(rbindata(1,ib)+rbindata(1,ib-1))
          endif
          flgfound=1
       endif
    enddo
    return
  end subroutine lagrange_radius

  subroutine all_lagrange_radii(iu,ofile)
    implicit none
    integer :: i,iu,ofile
    character(len=80) :: fname_lagrange
!    open(46,file=trim(fname_lagrange))
    do i=1,nlagrange
       call lagrange_radius(aquantile(i),arlagrange(i))
!       write(iu,100) arlagrange(i)
    enddo
    write(iu,100) ofile,tnow_myr,arlagrange(1:nlagrange)
!    close(46)
100 format(i5.5,1x,f6.0,1x,10(1x,f16.4))
    return
  end subroutine all_lagrange_radii

end module radial_profile
