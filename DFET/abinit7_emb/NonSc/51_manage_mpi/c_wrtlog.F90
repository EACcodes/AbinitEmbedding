#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! Writes output for Chen's code, into invks_log, contains info from abinit

subroutine c_wrtlog(msg)

  use interfaces_12_hide_mpi
  use interfaces_14_hidewrite, except_this_one => wrtout

  implicit none

  character(len=*), intent(in) :: msg;
  INTEGER :: me,ierr,nproc
  INTEGER,save :: master = 0
  
  me = 0

!Determine who I am
 call xmpi_me(me)
 call xmpi_nproc(nproc,ierr)

 if(me/=master) RETURN
  ! only master node write to invKS_log file
  ! open file  ...
  open(unit=2222,access='append',action='write',status='unknown', file='invKS_log',form='formatted');

  IF (msg(1:13) .eq. 'WELCOMEMSG_BF') THEN
    WRITE(2222,*)'-----------------------------------------------------------------'
    write(2222,*)'-            Embedding Potential Solver (ver. 2.3)              -'
    write(2222,*)'-                   Release  March/12/2010                      -'
    write(2222,*)'-                      (2009-2011)                              -'
    write(2222,*)'-                                                               -'
    write(2222,*)'-                        *  *  *                                -'
    write(2222,*)'-                                                               -'
    write(2222,*)'-                   Author: Chen Huang                          -'
    write(2222,*)'-                                                               -'
    write(2222,*)'-             Under GNU General Public License                  -'
    write(2222,*)'-           Based on the Version 6.0.3 of ABINIT                -'
    write(2222,*)'-           Support: NCPP (parallel and serial)                 -'
    write(2222,*)'-                    PAW  (serial)                              -'
    write(2222,*)'-----------------------------------------------------------------'
  ELSE
    write (2222, '(A)') trim(msg);
  ENDIF
  CLOSE (2222)

end subroutine
