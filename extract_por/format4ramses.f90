module format4ramses
  use commons4ramses
  use deviates_module
  implicit none
  integer :: com_mode,sort_mode,skip_mode
  integer :: np1,np2,ncalls_splitpops
  real(kind=8) :: mp1,mp2,mstar
  real(kind=8),dimension(6) :: com_rv
  integer,parameter :: npmax=5001000
  logical :: com_correction_firstrun=.true.
  real(kind=8) :: com_cpu0,com_cpu1

contains

  subroutine read_ramses_info(infofile)
    character(len=200) :: infofile
    character(len=15) :: chdum
    real(kind=8) :: tsu
    integer :: i,iu=98
    open(iu,file=trim(infofile))
    do i=1,8
       read(iu,*)
    enddo
    read(iu,'(a15,e21.15)') chdum,tsu
    do i=1,6
       read(iu,*)
    enddo
    read(iu,'(a15,e21.15)') chdum,slu
    read(iu,'(a15,e21.15)') chdum,sdu
    read(iu,'(a15,e21.15)') chdum,stu
    close(iu)
    tnow_myr=tsu*stu/3.15576d13
    return
  end subroutine read_ramses_info

  subroutine read_ramses_part(fileid,nthreads,nptot,all_part,id_part)
    character(len=200) :: fileid,infile
    character(len=5) :: nchar
    integer :: i,j,iu,ithread,nthreads,ndata,ndim,idum,npart,ierr
    integer :: nptot
    real(kind=8) :: ddum
    real(kind=8),dimension(11,npmax) :: all_part
!    real(kind=8),dimension(3,npmax) :: all_acc
    integer,dimension(npmax) :: id_part
    real(kind=8),allocatable,dimension(:)::xdp
    integer,allocatable,dimension(:)::idp,ill
    nptot=0
    do ithread=1,nthreads
       iu=10+ithread
       write(nchar,'(i5.5)') ithread
       infile=trim(fileid)//nchar
!       print *,'ithread,nchar =',ithread,'; "',trim(infile),'"'
       open(iu,file=trim(infile),form='unformatted',status='old',iostat=ierr)
       if (ierr/=0) then
          close(iu)
          print '(a,a,a)','Error wile opening file "',trim(infile),'"'
          print '(a,i2.2,a,i2)','iu',iu,': ierr =',ierr
          call close_program
!          exit
       endif
       read(iu) idum
       read(iu) ndim
       read(iu) npart
       read(iu) idum
       read(iu) idum
       read(iu) ddum
       read(iu) ddum
       read(iu) idum
!       print *,'npart =',npart
       allocate(xdp(1:npart))
       allocate(idp(1:npart))
       allocate(ill(1:npart))
!       ndata=0
!       do i=1,npart
!          read(iu,end=7) ddum
!          print *,'i =',i
!          ndata=ndata+1
!       enddo
!7      continue
!       print *,'ndata =',ndata
!       cycle
! positions
       do j=1,ndim
          read(iu) xdp
          do i=1,npart
!             read(iu) all_part(j,i+nptot)
             all_part(j,i+nptot)=xdp(i)*scalr
!             print *,'j,i =',j,i,all_part(j,i+nptot)
          enddo
       enddo
! velocities
       do j=1,ndim
          read(iu) xdp
          do i=1,npart
!             read(iu) all_part(j+3,i+nptot)
             all_part(j+3,i+nptot)=xdp(i)*scalv
          enddo
       enddoprint '(a,a,a)','Error wile opening file "',trim(infile),'"'
!       print *,'scalv =',scalv
! accelerations !NEW!
       if (existacc) then
          do j=1,ndim
             read(iu) xdp
             do i=1,npart
                all_part(j+7,i+nptot)=xdp(i)*scalacc
             enddo
          enddo
       endif
! mass
       read(iu) xdp
       do i=1,npart
!          read(iu) all_part(7,i+nptot)
          all_part(7,i+nptot)=xdp(i)*scalm
       enddo
! ID
       read(iu) idp
       do i=1,npart
!          read(iu) id_part(i+nptot)
          id_part(i+nptot)=idp(i)
       enddo
       if (usetp) then
! timestamp (birth epoch) !UNDER CONSTRUCTION!
          read(iu,end=7,err=8) ill
          read(iu,end=7,err=8) xdp
          do i=1,npart
             if (xdp(i)>0.d0) then
                all_part(11,i+nptot)=xdp(i)*scaltp
             else
                all_part(11,i+nptot)=-1.d15!older than the universe is safe
             endif
!                if (xdp(i)>0.d0) print *,'**** tp>0. ****'
          enddo
       endif
       goto 9
7      continue
       !print *,'End of file --> close'
       goto 9
8      print *,'File read erprint '(a,a,a)','Error wile opening file "',trim(infile),'"'ror -->close'
9      close(iu)
       nptot=nptot+npart
       print *,'Read thread No.',ithread,'; particle number so far:',nptot, ' (added ',npart,')'
       deallocate(xdp)
       deallocate(idp)
       deallocate(ill)
    enddo
  end subroutine read_ramses_part

  subroutine write_ramses_part(outfile,gpauxfile,nptot,all_part,id_part)
    character(len=200) :: outfile,gpauxfile
    character(len=5) :: nchar
    integer :: i,j,ii,ndim,idum
    integer :: nptot
    integer,dimension(npmax) :: status_part
    real(kind=8),dimension(11,npmax) :: all_part
!    real(kind=8),dimension(3,npmax) :: all_acc
    real(kind=8) :: ddum,vdum(3)
    integer,dimension(npmax) :: id_part
!    logical :: write_output=.true.
    call com_correction(nptot,all_part)
    call rotate_ramses_part(nptot,all_part)
    if (sort_mode==1) then
       call rqsort(nptot,id_part,ai_sorted)
    else
       do i=1,nptot
          ai_sorted(i)=i
       enddo
    endif
    call trim_ramses_part(nptot,all_part,status_part)
    open(90,file=trim(outfile))
    do i=1,nptot
       if (skip_mode>0.and.i/=skip_mode*(i/skip_mode)) cycle
       if (sort_mode==1) then
          ii=ai_sorted(i)
       else
          ii=i
       endif
       if (status_part(ii)==0) then
          if (flgtrim==2) then
             ddum=boxlength-boxrim
             call uniform3d(ddum,boxcenter,vdum)
             all_part(1:3,ii)=vdum
          else
             cycle
          endif
       endif
       if (usetp) then
          write(90,1090) all_part(1:7,ii),id_part(ii),all_part(8:10,ii),all_part(11,ii),tnow_myr-all_part(11,ii)
       else if (useacc) then
          write(90,1090) all_part(1:7,ii),id_part(ii),all_part(8:10,ii)
       else
          !This is the really important line for normal options (write ASCII files without acceleration and birth time).
          write(90,1090) all_part(1:7,ii),id_part(ii)
       endif
!       if (i>=1199000) write(*,1090) all_part(1:7,ii),id_part(ii)
    enddo
    close(90)
    open(91,file=trim(gpauxfile))
    write(91,*) 'tnow_myr=',tnow_myr
    write(91,'(a,i5,a)') 'tnow_lbl="',nint(tnow_myr),' Myr"'
    close(91)
1090 format(7(1x,e17.9),2x,i8,3(1x,e17.9),2(1x,e17.9))
  end subroutine write_ramses_part

  subroutine split_populations(nptot,all_part,id_part)
    !requires usetp to be true!
    implicit none
    integer :: nptot,i,i1,i2
    integer,dimension(npmax) :: id_part
    real(kind=8),dimension(11,npmax) :: all_part
    real(kind=8) :: accdum
    logical :: is_star(npmax)
    ncalls_splitpops=ncalls_splitpops+1
    nstars=0
    do i=1,nptot
       is_star(i)=all_part(11,i)>0.d0!positive age
       if (is_star(i)) nstars=nstars+1
    enddo
    if (ndm<0) then
       ndm=nptot-nstars
       print *,'ndm,nstars =',ndm,nstars
    endif
    i1=0
    i2=ndm
    ! re-write id_part to remove RAMSES' index clutter
    do i=1,nptot
       if (usetp) then
          if (is_star(i)) then!advance star index (appended to dm particles block)
             i2=i2+1
             !          id_part(i)=i2
             id_part(i)=id_part(i)+ndm
          else!advance dm index
             i1=i1+1
             !          id_part(i)=i1
          endif
       endif
       !debug: warning about zero accelerations (maybe normal for newborn stars)
!       if (ncalls_splitpops>1) then
!          accdum=sqrt(all_part(8,i)**2+all_part(9,i)**2+all_part(10,i)**2)
!          if (accdum==0.d0) print *,'Zero acceleration found at i =',i
!       endif
    enddo
    return
  end subroutine split_populations

  subroutine com_correction(nptot,all_part)
    integer :: nptot,i
    real(kind=8),dimension(11,npmax) :: all_part
    if (com_mode==0) return
    if (com_correction_firstrun.or.com_mode>0) then
       com_rv=0.d0
       do i=1,nptot
          com_rv(1:6)=com_rv(1:6)+all_part(1:6,i)/nptot!assume equal masses here
       enddo
       com_correction_firstrun=.false.
    endif
    print '(a,3(1x,f16.8))','COM_RV(1:3) =',com_rv(1:3)
    print '(a,3(1x,f16.8))','COM_RV(4:6) =',com_rv(4:6)
    do i=1,nptot
       if (abs(com_mode)==1.or.abs(com_mode)==3) then
          all_part(1:3,i)=all_part(1:3,i)-com_rv(1:3)
       endif
       if (abs(com_mode)==2.or.abs(com_mode)==3) then
          all_part(4:6,i)=all_part(4:6,i)-com_rv(4:6)
       endif
    enddo
    return
  end subroutine com_correction

  subroutine rotate_ramses_part(nptot,all_part)
    integer :: nptot,i
    real(kind=8) :: ddum,eps=1.d-10
    real(kind=8),dimension(11,npmax) :: all_part
    real(kind=8),dimension(3) :: p0,v0,p1,v1
    ddum=abs(a_rot)+abs(b_rot)+abs(c_rot)
    if (ddum<eps) return!do nothing if angles are (almost) zero
    do i=1,nptot
       p0(1:3)=all_part(1:3,i)
       v0(1:3)=all_part(4:6,i)
       call rotzxz(p0,p1,a_rot,b_rot,c_rot)
       call rotzxz(v0,v1,a_rot,b_rot,c_rot)
       all_part(1:3,i)=p1(1:3)
       all_part(4:6,i)=v1(1:3)
    enddo
    return
  end subroutine rotate_ramses_part

  subroutine trim_ramses_part(nptot,all_part,status_part)
    integer :: i,nptot,nout
    integer,dimension(npmax) :: status_part
    real(kind=8),dimension(11,npmax) :: all_part
    real(kind=8) :: dquad
    logical :: isout=.false.
    status_part=0
    nout=0
    call set_boxcenter
    do i=1,nptot
       dquad=max(abs(all_part(1,i)-boxcenter(1)),&
                 abs(all_part(2,i)-boxcenter(2)),&
                 abs(all_part(3,i)-boxcenter(3)))
       isout=dquad>=boxlength/2.d0-boxrim
       if (isout.and.flgtrim>=1) then
          status_part(i)=0
          nout=nout+1
!          print *,'pos =',all_part(1:3,i)-boxcenter(1:3),dquad
!          print *,'bc =',boxcenter
       else
          status_part(i)=1
       endif
    enddo
    if (nout>0) then
       if (flgtrim==2) then
          print *,'particles outside box - redistributed: ',nout
       else
          print *,'particles outside box - skipped: ',nout
       endif
    endif
    return
  end subroutine trim_ramses_part

!**********************************************************************
!     PIKAIA-Version of Quicksort
!     made F90-ready and converted to double precision by I. Thies in 2013,
!     adopted for format4ramses (integer array) in 2015
  subroutine rqsort(n,a,p)
!======================================================================
!     Return integer array p which indexes array a in increasing order.
!     Array a is not disturbed.  The Quicksort algorithm is used.
!
!     B. G. Knapp, 86/12/23
!
!     Reference: N. Wirth, Algorithms and Data Structures,
!     Prentice-Hall, 1986
!======================================================================
    implicit none

!     Input:
    integer :: n
!    real(kind=8) :: a(n)
    integer :: a(n)

!     Output:
    integer :: p(n)

!     Constants
    integer,parameter :: LGN=32, Q=11
!        (LGN = log base 2 of maximum n;
!         Q = smallest subfile to use quicksort on)

!     Local:
    real(kind=8)      x
    integer :: stackl(LGN),stackr(LGN),s,t,l,m,r,i,j

!     Initialize the stack
    stackl(1)=1
    stackr(1)=n
    s=1

!     Initialize the pointer array
    do 1 i=1,n
       p(i)=i
1   enddo

2   if (s.gt.0) then
       l=stackl(s)
       r=stackr(s)
       s=s-1

3      if ((r-l).lt.Q) then

!           Use straight insertion
          do 6 i=l+1,r
             t = p(i)
             x = a(t)
             do 4 j=i-1,l,-1
                if (a(p(j)).le.x) goto 5
                p(j+1) = p(j)
4            enddo
             j=l-1
5            p(j+1) = t
6         enddo
       else

!           Use quicksort, with pivot as median of a(l), a(m), a(r)
          m=(l+r)/2
          t=p(m)
          if (a(t).lt.a(p(l))) then
             p(m)=p(l)
             p(l)=t
             t=p(m)
          endif
          if (a(t).gt.a(p(r))) then
             p(m)=p(r)
             p(r)=t
             t=p(m)
             if (a(t).lt.a(p(l))) then
                p(m)=p(l)
                p(l)=t
                t=p(m)
             endif
          endif

!           Partition
          x=a(t)
          i=l+1
          j=r-1
7         if (i.le.j) then
8            if (a(p(i)).lt.x) then
                i=i+1
                goto 8
             endif
9            if (x.lt.a(p(j))) then
                j=j-1
                goto 9
             endif
             if (i.le.j) then
                t=p(i)
                p(i)=p(j)
                p(j)=t
                i=i+1
                j=j-1
             endif
             goto 7
          endif

!           Stack the larger subfile
          s=s+1
          if ((j-l).gt.(r-i)) then
             stackl(s)=l
             stackr(s)=j
             l=i
          else
             stackl(s)=i
             stackr(s)=r
             r=j
          endif
          goto 3
       endif
       goto 2
    endif
    return
  end subroutine rqsort

end module format4ramses
