subroutine c_wval1(W,nfft,etotal_cluster,etotal_env,extpot, &
                  ref_rhor,cellvol,usepaw)

  use defs_abitypes

  implicit none

  integer,intent(in) :: nfft,usepaw;
  real(dp), intent(in) :: & 
    extpot(nfft), &
    ref_rhor(nfft), & 
    etotal_cluster,etotal_env, &
    cellvol
  real(dp), intent(out) :: W
  real(dp) :: lagrange
  character(len=500) :: msg;

! >>> function body begins <<<!

! =============================
! Direct calculation for W
!
  if (usepaw==1) then 
    lagrange = dot_product(extpot,-ref_rhor)*cellvol/dble(nfft)
  else
    lagrange = dot_product(extpot,-ref_rhor)*cellvol/dble(nfft)
  endif
  W = -(etotal_cluster+etotal_env+lagrange)

  call c_wrtlog('--------------- Summary (c_Wval)---------------')
  WRITE(msg,'(a,ES16.8)') 'etotal_cluster : ', etotal_cluster;    call c_wrtlog(msg);
  WRITE(msg,'(a,ES16.8)') 'etotal_env     : ', etotal_env;        call c_wrtlog(msg);
  WRITE(msg,'(a,ES16.8)') 'Lagrange       : ', lagrange;          call c_wrtlog(msg);
  WRITE(msg,'(a,ES16.8)') 'W (no penalty) : ', W;                 call c_wrtlog(msg);
 
end subroutine

