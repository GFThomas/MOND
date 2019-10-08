module commons4ramses
  implicit none
  real(kind=8) :: pi=3.141592653589793238d0
  integer :: flgtrim=0,flghalfbox=0,flgcom4bin=0
  integer :: ndata,ndm,nstars,voldim,flg_agespread,flg_splitpops
  integer,parameter :: ndatamax=5001000
  integer,dimension(2,3) :: partition_commons
  real(kind=8) :: tnow_myr,age_thres,age_exp,agespread_heps
  real(kind=8) :: scalr,scalv,scalacc,scalm,scaltp
  real(kind=8) :: boxlength=1000.d0,boxrim=0.001d0
  real(kind=8) :: a_rot,b_rot,c_rot
  real(kind=8) :: stu,slu,sdu
  real(kind=8),dimension(3) :: boxcenter
  integer,dimension(ndatamax) :: id_data,ai_sorted
  real(kind=8),dimension(11,ndatamax) :: all_data
  logical :: first_set_boxcenter=.true.,existacc=.true.,useacc=.false.,usetp=.false.

contains

  subroutine read_data(filename)
    integer :: i
    character(len=200) :: filename
    open(15,file=trim(filename))
    call nextline(15)
    ndata=0
    do i=1,ndatamax
       ai_sorted(i)=i
       read(15,*,end=9) all_data(1:7,i),id_data(i)
       ndata=ndata+1
    enddo
9   close(15)
  end subroutine read_data

  subroutine set_boxcenter
    if (.not.first_set_boxcenter) return
    if (flghalfbox>0) then
       boxcenter=(/1.d0,1.d0,1.d0/)*boxlength/2.d0
       print *,'Re-set boxcenter to boxlength/2'
    endif
    first_set_boxcenter=.false.
  end subroutine set_boxcenter

end module commons4ramses
