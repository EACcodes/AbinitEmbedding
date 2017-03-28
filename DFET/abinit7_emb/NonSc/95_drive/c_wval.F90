subroutine c_Wval(W,nfftf,etotal_cluster,etotal_env,extpot, &
                  ref_rhor,cluster_den,env_den,cellvol,usepaw)

  use defs_abitypes

  implicit none

  integer,intent(in)   :: nfftf,usepaw;
  real(dp), intent(in) :: & 
    extpot(nfftf), &
    etotal_cluster,etotal_env, &
    ref_rhor(nfftf), & 
    cluster_den(nfftf),env_den(nfftf), &
    cellvol
  real(dp), intent(out) :: W
  real(dp) :: lagrange,emb_cluster,emb_env
  character(len=500) :: msg;

! >>> function body begins <<<!

! =============================
! Direct calculation for W
!
!! NOTE: 
!! make sure to include entropy term in the W
!!

  if (usepaw==1) then 
    lagrange = dot_product(extpot,cluster_den+env_den-ref_rhor)*cellvol/dble(nfftf)
  else
    emb_cluster = dot_product(extpot,cluster_den)*cellvol/dble(nfftf)
    emb_env     = dot_product(extpot,env_den    )*cellvol/dble(nfftf)
    lagrange    = dot_product(extpot,-ref_rhor)*cellvol/dble(nfftf)
  endif
  W = -(etotal_cluster+etotal_env+lagrange)

  call c_wrtlog('--------------- Summary (c_Wval)---------------')
  WRITE(msg,'(a,ES16.8)') ' etotal_cluster : ', etotal_cluster;    call c_wrtlog(msg);
  WRITE(msg,'(a,ES16.8)') ' etotal_env     : ', etotal_env;        call c_wrtlog(msg);
  WRITE(msg,'(a,ES16.8)') ' Lagrange       : ', lagrange;          call c_wrtlog(msg);
  WRITE(msg,'(a,ES16.8)') ' W (no penalty) : ', W;                 call c_wrtlog(msg);
  write(msg,'(a)')''; call c_wrtlog(msg)
  write(msg,'(a)')'KS-DFT energy for cluster (no embedding term)'; call c_wrtlog(msg);
  write(msg,'(a,es16.8)')'  cluster:',etotal_cluster - emb_cluster; call c_wrtlog(msg);
  write(msg,'(a,es16.8)')'  env    :',etotal_env     - emb_env    ; call c_wrtlog(msg);
  write(msg,'(a,es16.8)')'  vemb*(clu_rho+env_rho-ref_rho)  :', & 
    dot_product(extpot,cluster_den+env_den-ref_rhor)*cellvol/dble(nfftf) ; call c_wrtlog(msg);
  write(msg,'(a)')''; call c_wrtlog(msg)
  write(msg,'(a)')''; call c_wrtlog(msg)
 
end subroutine

