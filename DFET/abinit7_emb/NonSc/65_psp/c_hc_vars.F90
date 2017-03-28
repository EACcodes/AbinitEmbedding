module c_hc_vars

 use defs_basis
 use defs_datatypes

 implicit none

!=========================
! control parameters 
!=========================
 integer :: & 
   whoAmI, &    ! environment (2) or cluster (1)
   bLoadDen, &  ! load Density nref (true > 0, false < 0)
   bLoadVr, &   ! load Potential for restart  (true > 0, false < 0)
   bAllowSym , &  ! Switch on/off symmetries, you should switch this off
   bNCPP_PHI      ! unknown at the moment

!==================================
! internal arrays 
!==================================
 real(dp),allocatable :: extpot(:)       ! the embedding potential to obtain

 real(dp) :: hc_rmax            ! radial integration cutoff for PAW integration
 integer, allocatable :: in_paw(:)     ! =1 : if the point is in paw of one atom
                                       ! =-1: else
 integer :: & 
   hc_meshsz, &   ! max of mesh_size of all atoms
   hc_lsize, &    ! max of l + 1 for nonzero Gaunt coeff
   hc_lmax, &     ! Maximum value of angular momentum l+1
   hc_lmn2size    ! Maximum lmn2_size of all atoms
 
 real(dp) :: & 
   ext_energy     ! energy due to extpot, i.e. \int(extpot*hc_rho_paw(r))

 real(dp), allocatable :: & 
   hc_rho1(:,:,:), & ! n1
   hc_trho1(:,:,:), &! tilde-n1
   hc_rhor_paw(:)          ! cluster/env: trhor+rho1-trho1, no core density
 
end module
