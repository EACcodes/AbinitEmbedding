subroutine checkfile(fname)
!-------------------------------------
! if file exits this subroutine will
! not exit, until the file is removed
!-------------------------------------
implicit none

character(len=300) fname
character(len=500) stmp
logical :: ex

WRITE(stmp,'(a,a)') 'checking file:  ', trim(fname)
!CALL c_wrtlog(stmp)

DO 
  INQUIRE(file=fname, exist=ex)
  IF ( ex == .true. )  THEN
    CALL SYSTEM('sleep 1')
    CYCLE
  ENDIF
  EXIT
ENDDO


END subroutine


subroutine WaitForFile(fname)
!-------------------------------------
! if file does not exit, we wait
! otherwise we return 
!-------------------------------------
implicit none

character(len=300) fname
character(len=500) stmp
logical :: ex

WRITE(stmp,'(a,a)') 'waiting for file: ', trim(fname)
CALL c_wrtlog(stmp)

DO 
  INQUIRE(file=fname, exist=ex)
  IF ( ex == .FALSE. )  THEN
    CALL SYSTEM('sleep 1')
    CYCLE
  ELSE
    EXIT
  ENDIF
ENDDO

WRITE(stmp,'(a,a)') 'waitForFile:  found file => ', fname
CALL wrtout(06,stmp,'COLL')

END subroutine WaitForFile

subroutine c_checkFileLen(fname, recNum)
!------------------------------------------------------
! fname is supposed to be a single-column data
! file. recNum is the desired number of data in fname
! Code will stop if the length of fname is not equal 
! to recNum.
!------------------------------------------------------

implicit none

character(len=*), intent(in) :: fname
character(len=500)   :: msg
integer, intent(in ) :: recNum  ! record number to be compared to
integer              :: i,fs,n
real(kind=8) :: tmp

!>>>>>>>> FUNCTION <<<<<<<<<<!
write(msg,'(a,a)') '(CheckFileLen) checking file = ',TRIM(fname)
call c_wrtlog(msg)
open (unit=111,file=TRIM(fname),status='old',action='read',iostat=fs)
if (fs/=0) then
  call c_wrtlog("(CheckFileLen) Error in Opening File. STOP!")
  call leave_new('COLL')
endif

! Start to count 
n = 0
do while (1==1)
  read(111,*,iostat=fs) tmp
  if (fs/=0) then
    exit
  else
    n = n + 1
  endif
enddo
close (111)

! test record number
if (n/=recNum) then
  write(msg,*) "Record number of ", TRIM(fname)," does not match ", recNum, " STOP!"
  call c_wrtlog(msg)
  stop
else
  write(msg,'(a,a,I15,a,I15,a)')& 
   ADJUSTL(TRIM(fname))," has record number:", n, " == ", recNum, "(target), test passed. "
  call c_wrtlog(msg)
endif

end subroutine c_checkFileLen
