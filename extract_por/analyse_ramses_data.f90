  program analyse_ramses_data
    use format4ramses
    use radial_profile
    use partial_com
    use xyimage
!    use commons4ramses, only : ndata,useacc,scalr,scalv,scalacc,scalm,scaltp,age_thres
    implicit none
    integer :: i,nfile,ifile,ofile,nptot,iloop,nloop=1,nskipped,&
         analysis_mode,flg_compart,flg_image,edgemode,massbinmode,nskip_particles
    integer :: update_mode,n1,n2,n3,nthreads,iucompart=55,iutimes=56,iulagrange=57,nstat,idum
    character(len=5) :: ncharin,ncharout,ncharx
    character(len=30) :: parfile
    character(len=200) :: fmain,fpath,fnametrunc,infofile,xfile,outfile,gpauxfile,&
         radialfile,rlagrangefile,imagefile
    integer,parameter :: nmax=5001000
    real(kind=8) :: dt=1.d2,frac,adeg,bdeg,cdeg,rdum
    real(kind=8),dimension(11,nmax) :: all_part
    integer,dimension(nmax) :: id_part
    logical :: file_exists,prev_existed=.false.
    call cpu_time(com_cpu0)
    CALL get_command_argument(1, parfile)
    IF (LEN_TRIM(parfile)==0) parfile='fmtramses.par'
    open(10,file=trim(parfile))
    read(10,*) fmain
    read(10,*) nfile
    read(10,*) nthreads
    read(10,*) update_mode
    read(10,*) com_mode
    read(10,*) sort_mode
    read(10,*) skip_mode
    read(10,*) scalr
    read(10,*) scalv
    read(10,*) scalm
    read(10,*) scalacc,scaltp
    read(10,*) existacc
    read(10,*) adeg,bdeg,cdeg
    read(10,*)
    read(10,*) analysis_mode,voldim
    read(10,*) flg_splitpops,ndm
    read(10,*)
    read(10,*) rbinmax,nbindata
    read(10,*) boxcenter(1:3)
    read(10,*) boxlength,flgtrim,boxrim
    read(10,*) flghalfbox,flgcom4bin
    read(10,*)
    read(10,*) flg_compart,frac,dt
    read(10,*) n1,n2,n3
    read(10,*) nstat
    read(10,*) 
    read(10,*) flg_image,edgemode,massbinmode
    read(10,*) age_thres,age_exp,flg_agespread,agespread_heps
    read(10,*) cx_xyimage,cy_xyimage
    read(10,*) wx_xyimage,wy_xyimage
    read(10,*) nbx_xyimage,nby_xyimage
    read(10,*) hxpix_xyimage,hypix_xyimage
    read(10,*) nskip_particles
    read(10,*) rgb_scallevel,rgb_maxlevel,rgb_xgamma
    read(10,*)
    nlagrange=0
    do i=1,nlagrangemax
       read(10,*) rdum
       if (rdum<=0.d0.or.rdum>1.d0) exit
       nlagrange=i
       aquantile(i)=rdum
    enddo
    close(10)
    if (rgb_scallevel<=0.d0) rgb_scallevel=1.d0
    if (rgb_maxlevel<=0.d0) rgb_maxlevel=255.d0
    a_rot=adeg*pi/180.d0
    b_rot=bdeg*pi/180.d0
    c_rot=cdeg*pi/180.d0
    useacc=scalacc>0.d0
    usetp=scaltp>0.d0
    if (wy_xyimage<=0.d0) then
       wy_xyimage=wx_xyimage
    endif
    if (nby_xyimage<=0) then
       nby_xyimage=nbx_xyimage
    endif
    if (nfile<=0) then
       nfile=abs(nfile)
       nloop=nfile
       iloop=0; ofile=0
    else
       nloop=nfile-1
       iloop=nfile-1; ofile=nfile-1
    endif
    if (analysis_mode==1.or.analysis_mode<0) then
       open(iutimes,file='times.dat')!read from beginning
    else
       open(iutimes,file='times.dat',position='append')!append new data
    endif
    if (analysis_mode>=1) then
       open(iulagrange,file='rlagrange.dat',position='append')!append new data
       write(iulagrange,'(a5,1x,a6,1x,10(1x,f16.6))') '#iout','t_myr',aquantile(1:nlagrange)
    endif
    if (flg_compart>=1.and.analysis_mode>=1) then
       open(iucompart,file='partial-com.dat',position='append')
    endif
    nskipped=0
!    do iloop=1,nloop+nskipped
    do while (iloop<=nloop)
       iloop=iloop+1
       if (nloop==1) then
          ifile=nfile
       else
          ifile=iloop
          ofile=ofile+1
       endif
       write(ncharin,'(i5.5)') ifile
       write(ncharout,'(i5.5)') ofile
       fpath=trim(fmain)//'output_'//ncharin//'/'
       fnametrunc=trim(fpath)//'part_'//ncharin//'.out'
       infofile=trim(fpath)//'info_'//ncharin//'.txt'
       outfile='part_'//ncharout//'.asc'
       gpauxfile='gpaux_'//ncharout//'.asc'
       print *,'outfile = ',trim(outfile)
!       print *,'fpath     =',fpath
!       print *,'fnametrunc=',fnametrunc
!    stop
       inquire(file=trim(outfile),exist=file_exists)
       if (analysis_mode==1.or.analysis_mode<0) then
          if (.not.file_exists) exit
          call read_data(outfile)
          read(iutimes,*) idum,tnow_myr
       else
          if (update_mode>=0) then
             if (file_exists) then
                if (update_mode==1) then
                   nskipped=nskipped+1
                endif
                print '(2a)','Skipping existing file: ',outfile
                prev_existed=.true.
                first_write_comdata=.false.
                cycle
             else if (.not.file_exists.and.prev_existed.and.update_mode>=1) then
! this triggers the re-numberung and also skips the first continuation file,
!  avoiding dupes.
                iloop=iloop-nskipped
                ofile=ofile-1
                prev_existed=.false.
                cycle
             endif
          endif
          call read_ramses_part(fnametrunc,nthreads,nptot,all_part,id_part)
          call read_ramses_info(infofile)
!          if (flg_splitpops>=1) then
          !this re-writes id_part!
          call split_populations(nptot,all_part,id_part)
!          endif
          write(iutimes,'(a5,2x,f9.2)') ncharout,tnow_myr
          call write_ramses_part(outfile,gpauxfile,nptot,all_part,id_part)
          print *,'Processed files ',trim(fnametrunc),'*'
          print *,'Written file    ',outfile
!          print *,'*********************************************'
          if (analysis_mode>=1.OR.flg_image>=1) then
             call import_data(nptot,all_part,id_part)
          endif
       endif
!       if (analysis_mode>=1) then
!          radialfile='radprofile_'//ncharout//'.dat'
!          rlagrangefile='rlagrange_'//ncharout//'.dat'
!       endif
       if (flg_image>=1) imagefile='image_'//ncharout//'.dat'
       if (analysis_mode>=1) then
          radialfile='radprofile_'//ncharout//'.dat'
!          rlagrangefile='rlagrange_'//ncharout//'.dat'
          print *,'radialfile =',trim(radialfile)
          call radial_binning
          call radial_density
          call radial_vdispersion
          call check_angmom
          call write_data(radialfile)
          call all_lagrange_radii(iulagrange,ofile)
!---- Partial COM calculation
          if (flg_compart>=1) then
             call get_partition(n1,n2,n3)
             if (flg_compart==2) then
                call get_dcen(1,n1,frac,nstat)
                call get_dcen(2,n2,frac,nstat)
                call get_dcen(3,n3,frac,nstat)
                all_com=all_dcen
             else
                call get_com(1)
                call get_com(2)
                call get_com(3)
             endif
             call write_comdata(iucompart,ifile,dt)
          endif
       endif
!---- Image table file
       if (flg_image>=1) then
          print *,'imagefile = ',trim(imagefile)
          if (flg_image==1) then
             call get_xyimage(edgemode,massbinmode,nskip_particles)
             call write_xyimage(imagefile)
          else
             call get_rgb_xyimage(edgemode,massbinmode,nskip_particles)
             call write_rgb_xyimage(imagefile,flg_image-2)
          endif
       endif
       print *,'******************************************************************************************'
       print *,'******************************************************************************************'
       print *
       call flush
    enddo
    if (flg_compart>=1.and.analysis_mode>=1) then
       close(iucompart)
    endif
    close(iutimes)
    close(iulagrange)
    call close_program
  end program analyse_ramses_data

  subroutine close_program
    use format4ramses
    call cpu_time(com_cpu1)
    print *
    print '("CPU time = ",f10.3," seconds.")',com_cpu1-com_cpu0
    print *
    stop
  end subroutine close_program
