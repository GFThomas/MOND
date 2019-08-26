!C***********************************************************************
!C
!C
                          SUBROUTINE nextline(iufile)
!C
!C
!C***********************************************************************
!C
!C    Read on unit ufile until next line without comment-character
!C           and with some text
!C------------------------------------------------------------------
!C 
      CHARACTER*100 line
      CHARACTER*1 first
      INTEGER iufile,ufile,nskip,i
      LOGICAL verbose
      ufile=abs(iufile)
      verbose=iufile<0
      nskip=0
10    READ(ufile,'(a)',ERR=50,END=60) line
!Cold      WRITE(first,'(a1)') line(1:1)
      DO 20 i=1,100
         first=line(i:i)
         IF (first.ne. ' '.AND. first.ne.char(0)) EXIT
 20   CONTINUE
!c      first = line(1:1)
!c      WRITE(*,'(5A)') 'line = "',line,'", first = "',first,'"'
      IF(first.NE.'#' .AND. first.NE.'C' .AND. first.NE.'c' .AND.       &
     &   first.NE.';' .AND. first.NE.'%' .AND. first.NE.'!' .AND.       &
     &   line.NE.' ') THEN
         BACKSPACE ufile
         IF (verbose) THEN
            WRITE(*,*)
            WRITE(*,*) 'BACKSPACE: Skipped',nskip,' lines'
            WRITE(*,*)
         ENDIF
         R E T U R N
      ELSE
         nskip=nskip+1
         GOTO 10
      ENDIF

50    CONTINUE
      PRINT*,' nextline: error in reading line from unit ',ufile
!C          ----
      STOP
60    CONTINUE
      PRINT*,' nextline: end of file (may be regular)'
      BACKSPACE ufile!let calling routine finish reading
      RETURN
      END
