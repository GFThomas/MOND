  MODULE rnd_array_module
    IMPLICIT NONE
    INTEGER :: rnd_ioseed=0,rnd_cmode=0,rnd_kseed=37,rnd_dimseed
    INTEGER :: rnd_callcount=0,rnd_initcount=0
    INTEGER :: rnd_verbose=0

  CONTAINS
    
    SUBROUTINE wrap_get_rnd(nx,xxx)
! Version with seed hidden in common block
      IMPLICIT NONE
!      INTEGER iseed
!      COMMON /ioseed/ iseed
      INTEGER nx
      REAL(KIND=8) :: xxx(nx)
!      write(*,*) 'rnd_ioseed =',rnd_ioseed
      CALL get_rnd_array(rnd_ioseed,nx,xxx)
      RETURN
    END SUBROUTINE wrap_get_rnd

    SUBROUTINE wrap_one_rnd(x)
! Simple random number, no array
! Version with seed hidden in common block
      IMPLICIT NONE
!      INTEGER iseed
!      COMMON /ioseed/ iseed
      REAL(KIND=8) :: xxx(1),x
      CALL get_rnd_array(rnd_ioseed,1,xxx)
      x=xxx(1)
      RETURN
    END SUBROUTINE wrap_one_rnd

    SUBROUTINE wrap_ind_rnd(x,ix,irange)
! Simple integer random number, no array
! Version with seed hidden in common block
      IMPLICIT NONE
!      INTEGER iseed
!      COMMON /ioseed/ iseed
      INTEGER :: ix,irange,imax=2147483647
      REAL(KIND=8) :: xxx(1),x,dimax=4294967296.d0
      CALL get_rnd_array(rnd_ioseed,1,xxx)
      x=xxx(1)
      IF (irange==0) THEN
         ix=int(imax*x)
      ELSE IF (irange<0) THEN
         ix=floor((x-0.5d0)*dimax)
      ELSE
         ix=int(irange*x)
      ENDIF
      RETURN
    END SUBROUTINE wrap_ind_rnd

    SUBROUTINE get_rnd_array(iseed,nx,xxx)
      IMPLICIT NONE
      INTEGER :: iseed,nx
      REAL(KIND=8) :: xxx(nx),xdum(1)
      rnd_callcount=rnd_callcount+1
      IF (iseed.NE.-1) CALL init_random_seed(iseed)
      CALL RANDOM_NUMBER(xxx)
      IF (rnd_cmode>=1) THEN!experimental closed interval mode (req. one more rnd call)
         CALL RANDOM_NUMBER(xdum)
         IF (xdum(1)<0.5d0) THEN
            xxx(1:nx)=1.d0-xxx(1:nx)
         ENDIF
      ENDIF
      iseed=-1!seed is now initialized
      RETURN
    END SUBROUTINE get_rnd_array

    SUBROUTINE get_rnd_number(iseed,x)
! Simple random number, no array
      IMPLICIT NONE
      INTEGER :: iseed
      REAL(KIND=8) :: xxx(1),x
      CALL get_rnd_array(iseed,1,xxx)
      x=xxx(1)
      RETURN
    END SUBROUTINE get_rnd_number

    SUBROUTINE scale_rnd_array(nx,xsc,xlo,xhi)
! Version with seed hidden in common block
      IMPLICIT NONE
!      INTEGER iseed
!      COMMON /ioseed/ iseed
      INTEGER :: nx
      REAL(KIND=8) :: xxx(nx),xsc(nx),xlo,xhi
!      write(*,*) 'rnd_ioseed =',rnd_ioseed
      CALL get_rnd_array(rnd_ioseed,nx,xxx)
      xsc(:)=xlo+(xhi-xlo)*xxx(:)
      RETURN
    END SUBROUTINE scale_rnd_array

    SUBROUTINE init_random_seed(iseed)
      INTEGER :: i, n, clock, iseed
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
      rnd_initcount=rnd_initcount+1
      CALL RANDOM_SEED(size = n)
      rnd_dimseed=n
      ALLOCATE(seed(n))
      if (rnd_verbose>=1) write(*,*) 'n =',n
      IF (iseed<=0) THEN
         CALL SYSTEM_CLOCK(COUNT=clock)
      ELSE
         clock = iseed
      ENDIF
!      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      seed = clock + rnd_kseed * (/ (i - 1, i = 1, n) /)
      if (rnd_verbose>=1) then
         do i=1,n
            print *,'i,seed(i) =',i,seed(i)
         enddo
      endif
      CALL RANDOM_SEED(PUT = seed)
      
      DEALLOCATE(seed)
    END SUBROUTINE init_random_seed

    SUBROUTINE iter_random_seed(iter,iseed)
      implicit none
      integer :: i,j,iter,iseed,oseed,n
      real(kind=8) :: dimax=4294967296.d0
      real(kind=8),dimension(:),allocatable :: xxx
      INTEGER, DIMENSION(:), ALLOCATABLE :: aseed
      call init_random_seed(iseed)
      print *,'iter,iseed,rnd_dimseed=',iter,iseed,rnd_dimseed
      iseed=-1
      allocate(aseed(1:rnd_dimseed))
      allocate(xxx(1:rnd_dimseed))
      do i=1,iter
         call RANDOM_NUMBER(xxx)
         do j=1,rnd_dimseed
            aseed(j)=floor((xxx(j)-0.5d0)*dimax)
         enddo
      enddo
      if (rnd_verbose>=1) then
         do j=1,rnd_dimseed
            print *,'i,aseed(j) =',j,aseed(j)
         enddo
      endif
      if (iter>0) call RANDOM_SEED(put=aseed)
      return
    END SUBROUTINE iter_random_seed

  END MODULE rnd_array_module
