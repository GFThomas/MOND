module partial_com
  use commons4ramses
  implicit none
  integer,dimension(2,3) :: partition_n
  real(kind=8),dimension(6,3) :: all_com,all_dcen
  logical :: first_write_comdata=.true.

contains

  subroutine get_partition(n1,n2,n3)
    integer :: i,n1,n2,n3
    partition_n(1,1)=1
    partition_n(2,1)=n1
    partition_n(1,2)=n1+1
    partition_n(2,2)=n1+n2
    partition_n(1,3)=n1+n2+1
    partition_n(2,3)=n1+n2+n3
    if (partition_n(2,3)>ndata) then
       print *,'GET_PARTITION warning: partition exceeds ndata'
    endif
    call export_partition!dirty hack to export partition data to new commons module
    return
  end subroutine get_partition

  subroutine export_partition
    partition_commons=partition_n
    return
  end subroutine export_partition

  subroutine get_com(igroup)
    integer :: i,n,igroup,nstart,nend
    real(kind=8),dimension(6) :: pv,com
    com=0.d0
    nstart=min(partition_n(1,igroup),ndata)
    nend  =min(partition_n(2,igroup),ndata)
    n=nend-nstart+1
    print '(a,3(1x,i7))','n,nlo,nhi =',n,nstart,nend
    do i=nstart,nend
       pv(1:6)=all_data(1:6,i)
       com(1:6)=com(1:6)+pv(1:6)/n
!       if (igroup==3) print *,'pos =',pv(1:3)
!       if (pv(2)<1600.d0) print *,'igroup =',igroup
!       if (igroup==1.and.pv(2)<1500.d0) print *,'pos =',pv(1:3)
    enddo
    all_com(1:6,igroup)=com(1:6)
    return
  end subroutine get_com

  subroutine get_dcen(igroup,n,frac,nstat)
!EXPERIMENTAL!
    integer :: i,ii,n,nlo,nhi,igroup,nstat
    real(kind=8) :: frac,mthres=1.d99!safely take all particle masses
    real(kind=8),dimension(n) :: om
    real(kind=8),dimension(3,n) :: opos
    real(kind=8),dimension(3) :: posc
    real(kind=8) :: rdens,avdens
    nlo=partition_n(1,igroup)
    nhi=partition_n(2,igroup)
!    n=nhi-nlo+1
    do i=1,n
       ii=i+nlo-1
       om(i)=all_data(7,ii)
       opos(1:3,i)=all_data(1:3,ii)
!       if (100*(i/100)==i) print *,'opos =',opos(1:3,i)
    enddo
    print '(a,i1,a)','Calculation of density center ',igroup,'...'
    print '(a,3(1x,i7))','n,nlo,nhi =',n,nlo,nhi
!    print *,'opos =',opos(1:3,1)
    call densitycenter(3,OM,OPOS,n,FRAC,MTHRES,POSC,RDENS,AVDENS,NSTAT)
    all_dcen(1:3,igroup)=posc(1:3)
    all_dcen(4:6,igroup)=0.d0
    return
  end subroutine get_dcen
  
  subroutine write_comdata(iu,ifile,dt)
    integer :: iu,ifile
    real(kind=8) :: dt
    real(kind=8),dimension(3) :: dist_com
    dist_com(1)=sqrt((all_com(1,1)-all_com(1,2))**2+(all_com(2,1)-all_com(2,2))**2+(all_com(3,1)-all_com(3,2))**2)
    dist_com(2)=sqrt((all_com(1,2)-all_com(1,3))**2+(all_com(2,2)-all_com(2,3))**2+(all_com(3,2)-all_com(3,3))**2)
    dist_com(3)=sqrt((all_com(1,3)-all_com(1,1))**2+(all_com(2,3)-all_com(2,1))**2+(all_com(3,3)-all_com(3,1))**2)
    if (first_write_comdata) then
       write(iu,1001) '#t/Myr','cx1','cy1','cz1','cx2','cy2','cz2','cx3','cy3','cz3','d(1:2)','d(2:3)','d(3:1)'
       first_write_comdata=.false.
    endif
    write(iu,1000) nint((ifile-1)*dt),all_com(1:3,1),all_com(1:3,2),all_com(1:3,3),dist_com(1:3)
1000 format(i7,2x,3(2x,3(2x,f12.4)),4x,3(2x,f12.4))
1001 format(a7,2x,3(2x,3(2x,a12)),  4x,3(2x,a12))
    return
  end subroutine write_comdata

end module partial_com
