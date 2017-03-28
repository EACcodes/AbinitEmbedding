module orbital_shifts

  use get_xc_potential
  use fourierhlp

  use defs_basis
  use defs_datatypes
  use defs_abitypes
  use defs_wvltypes
  use m_xmpi

  use m_hamiltonian
  use m_pawcprj

  !! Module to calculate the orbital shifts in Hyman OEP approach
  !!
  !! Following Chen Huang's hack, we expect the VTORHO to call 
  !! scfcv -> vtorho -> vtowf 
  !! vtowf is hacked to call this function if solve_orbital_shift is true

contains


  !{\src2tex{textfont=tt}}
  !!****f* ABINIT/cgwf
  !! NAME
  !! cgwf
  !!
  !! FUNCTION
  !! Update all wavefunction |C>, non self-consistently.
  !! also compute the corresponding H|C> and Vnl|C> (and S|C> if paw).
  !! Uses a conjugate-gradient algorithm.
  !! In case of paw, resolves a generalized eigenproblem using an
  !!  overlap matrix (not used for norm conserving psps).
  !!
  !! COPYRIGHT
  !! Copyright (C) 1998-2011 ABINIT group (DCA, XG, GMR, MT)
  !! This file is distributed under the terms of the
  !! GNU General Public License, see ~abinit/COPYING
  !! or http://www.gnu.org/copyleft/gpl.txt .
  !! For the initials of contributors, 
  !! see ~abinit/doc/developers/contributors.txt .
  !!
  !! INPUTS
  !!  berryopt == 4: electric field is on;
  !!              5: magnetic field is on;
  !!              all other values, no field is present
  !!  chkexit= if non-zero, check whether the user wishes to exit
  !!  cpus = CPU time limit
  !!  dimffnl=second dimension of ffnl (1+number of derivatives)
  !!  ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
  !!  filnam_ds1=name of input file (used for exit checking)
  !!  filstat=name of the status file
  !!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
  !!  icg=shift to be applied on the location of data in the array cg
  !!  ikpt=number of the k-point
  !!  inonsc=index of non self-consistent loop
  !!  isppol=spin polarization currently treated
  !!  kg_k(3,npw)=coordinates of planewaves in basis sphere.
  !!  kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
  !!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
  !!        =if useylm=0, max number of (l,n)   comp. over all type of psps
  !!  matblk=dimension of the array ph3d
  !!  mband =maximum number of bands
  !!  mcg=second dimension of the cg array
  !!  mgfft=maximum size of 1D FFTs
  !!  mkgq = second dimension of pwnsfacq
  !!  mpi_enreg=informations about MPI parallelization
  !!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
  !!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
  !!  mpw=maximum dimensioned size of npw
  !!  natom=number of atoms in cell.
  !!  nband=number of bands.
  !!  nbdblock=number of bands in a block
  !!  nkpt=number of k points
  !!  nline=number of line minimizations per band.
  !!  nloalg(5) data concerning nonlop application
  !!  npw=number of planewaves in basis sphere at given k.
  !!  npwarr(nkpt)=number of planewaves in basis at this k point
  !!  nspinor=number of spinorial components of the wavefunctions
  !!  ntypat=number of types of atoms in cell.
  !!  nvloc=final dimension of vlocal (usually 1, but 4 for non-collinear
  !!  n4,n5,n6 used for dimensionning of vlocal
  !!  ortalg=governs the choice of the algorithm for orthogonalisation.
  !!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
  !!  prtvol=control print volume and debugging output
  !!  pwind(pwind_alloc,2,3) = array used to compute
  !!           the overlap matrix smat between k-points (see initberry.f)
  !!  pwind_alloc = first dimension of pwind
  !!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
  !!                           (see initberry.f)
  !!  pwnsfacq(2,mkgq) = phase factors for the nearest neighbours of the
  !!                     current k-point (electric field, MPI //)
  !!  tolwfr=tolerance on largest wf residual
  !!  vlocal(n4,n5,n6,nvloc)= local potential in real space, 
  !!             on the augmented fft grid
  !! to the derivative of XC energy with respect to
  !!    kinetic energy density, in real space, 
  !! on the augmented fft grid. (optional argument)
  !!  wfoptalg=govern the choice of algorithm for wf optimisation
  !!   (0, 1, 10 and 11 : in the present routine, usual CG algorithm ;
  !!   (2 and 3 : use shifted square Hamiltonian)
  !!
  !! OUTPUT
  !!  dphase_k(3) = change in Zak phase for the current k-point 
  !! in case berryopt = 4 (electric field)
  !!  resid(nband)=wf residual for new states=|(H-e)|C>|^2 (hartree^2)
  !! matrix expressed in sthe WFs subspace
  !! Hamiltonian expressed in sthe WFs subspace
  !!
  !! SIDE EFFECTS
  !!  cg(2,mcg)
  !!    at input =wavefunction <G|C band,k> coefficients for ALL bands
  !!    at output same as input except that
  !!      the current band, with number 'band' has been updated
  !!      calculations (see initberry.f)
  !!  if(gs_hamk%usepaw==1)
  !!   gsc(2,mgsc)=<G|S|C band,k> coefficients for ALL bands
  !!               where S is the overlap matrix (used only for paw)
  !!
  !! NOTES
  !!  cg should not be filtered and normalized : it should already
  !!   be OK at input !
  !!  Not sure that that the generalized eigenproblem (when gs_hamk%usepaw=1)
  !!   is compatible with wfoptalg=2 or 3 (use of shifted square
  !!   Hamiltonian) - to be verified
  !!
  !! PARENTS
  !!      vtowfk
  !!
  !! CHILDREN
  !!      bestwfs,chkexi,dotprod_g,etheta,getghc,leave_new,linemin
  !!      mpi_bcast,mpi_recv,mpi_send,precon,projbd,smatrix,
  !!      sqnorm_g,status,timab,wrtout,xcomm_init,xsum_mpi
  !!
  !! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

  subroutine orthogonalize(A, B, istwf_k, mpi_enreg, N )

    implicit none

    real*8, intent(inout) :: A(:,:)
    real*8, intent(in) :: B(:,:)

    type(MPI_type), intent(inout) :: mpi_enreg

    integer, intent(in):: istwf_k
    integer, intent(in):: N

!!!!!!!!

    real*8 dotr, doti

    integer ipw

!!!!!!!!!
    ! project out the B part
!
    call dotprod_g(dotr,doti,istwf_k,&
              mpi_enreg,N,2,B,A)          
    !
    do ipw=1,N
       A(1,ipw) = A(1,ipw) - dotr * B(1,ipw) + doti * B(2,ipw)
       A(2,ipw) = A(2,ipw) - doti * B(1,ipw) - dotr * B(2,ipw)
    enddo

  end subroutine orthogonalize

  subroutine cgwf_orbital_shift(berryopt,cg,chkexit,cpus,&
       &                dimffnl,eigen,&
       &                ffnl,filnam_ds1,filstat,&
       &                gs_hamk,icg,ikpt,&
       &                isppol,kg_k,kinpw,lmnmax,matblk,mband,&
       &                mcg,mgfft,mkgq,mpi_enreg,mpsang,&
       &                mpssoang,mpw,natom,nband,nbdblock,nkpt,nline,&
       &                nloalg,npw,npwarr,nspinor,&
       &                ntypat,nvloc,n4,n5,n6,occ,ortalg,&
       &                paral_kgb,ph3d,prtvol,pwind,pwind_alloc,pwnsfac,&
       &                pwnsfacq,run_params, tolwfr,vhartree_new, &
       &                vhartree_old, vlocal,wfoptalg) 

#if defined HAVE_MPI2
    use mpi
#endif

    !This section has been created automatically by the script Abilint (TD).
    !Do not modify the following lines by hand.
    use interfaces_14_hidewrite
    use interfaces_16_hideleave
    use interfaces_18_timing
    use interfaces_32_util
    use interfaces_51_manage_mpi
    use interfaces_53_abiutil
    use interfaces_53_spacepar
    use interfaces_59_io_mpi
    use interfaces_65_nonlocal
    use interfaces_66_paw
    use interfaces_66_wfs
    use interfaces_67_common
    !End of the abilint section

    implicit none

#if defined HAVE_MPI1
    include 'mpif.h'
#endif

    !Arguments ------------------------------------

    integer, intent(in) :: berryopt,chkexit,dimffnl,icg
    integer, intent(in) :: ikpt,isppol,lmnmax,matblk
    integer, intent(in) :: mband,mcg,mgfft,mkgq
    integer, intent(in) :: mpsang,mpssoang,mpw,n4
    integer, intent(in) :: n5,n6,natom,nband,nbdblock,nkpt
    integer, intent(in) :: nline,npw,nspinor,ntypat
    integer, intent(in) :: nvloc,ortalg,paral_kgb,prtvol,pwind_alloc
    integer, intent(in) :: wfoptalg

    real(dp), intent(in) :: cpus,tolwfr

    character(len=fnlen), intent(in) :: filnam_ds1,filstat

    type(MPI_type), intent(inout) :: mpi_enreg
    type(gs_hamiltonian_type), intent(inout) :: gs_hamk

    integer, intent(in) :: kg_k(3,npw),nloalg(5),npwarr(nkpt)
    integer, intent(in) :: pwind(pwind_alloc,2,3)

    real(dp), intent(in) :: ffnl(npw,dimffnl,lmnmax,ntypat)
    real(dp), intent(in) :: kinpw(npw)
    real(dp), intent(in) :: pwnsfac(2,pwind_alloc),pwnsfacq(2,mkgq)

    real(dp), intent(in) :: cg(2,mcg)
    real(dp), intent(inout) :: ph3d(2,npw,matblk),vlocal(n4,n5,n6,nvloc)

    real(dp), intent(in) :: vhartree_new(gs_hamk%nfft), &
         & vhartree_old(gs_hamk%nfft)

    real(dp), intent(in) :: occ(mband*nkpt*2)
    real(dp), intent(in) :: eigen(mband*nkpt*2)

    type (run_parameter_type), intent(inout) :: run_params

    !Local variables-------------------------------

    real(dp) maxocc

#if defined HAVE_MPI
    integer :: tag
# else
    integer :: openexit
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    integer icg_shift, ipw

    !! orbital shift coefficients recip space
    REAL(dp), allocatable :: orb_shift(:,:) 

    integer iband, ibandmin, ibandmax, ibdblock ! counter for bands
    integer cptopt

    ! Working variables for CG
    integer icg_step !! loop variable for CG loop

    integer ntot_k !! number of grid points for this k-point

    integer nblock, iblock !! loop over blocks of bands

    integer sij_opt !! paw option, set to zero here\
    
    real(dp) lambda ! dummy lambda for get_ghc, not needed here

    complex(dp) zg_alpha

    real(dp) cg_beta, dotr, doti, dnorm, dnorm_old, normb

    real(dp), allocatable :: potloc(:,:)     ! real-space tmp potential

    real(dp), allocatable :: rr_b(:,:)     ! real-space CG right-hand side
    real(dp), allocatable :: ri_b(:,:)     ! real-space CG right-hand side
    real(dp), allocatable :: wf_real(:,:)  ! real-space KS orbital, real part
    real(dp), allocatable :: wf_imag(:,:)  ! real-space KS orbital, imag part

    real(dp), allocatable :: cg_p(:,:)     ! helper array in CG
    real(dp), allocatable :: cg_resid(:,:) ! helper array in CG
    real(dp), allocatable :: cg_Ap(:,:)    ! helper array in CG

    real(dp), allocatable :: cg_b(:,:) ! right hand side of orbshift equation

    real(dp), allocatable :: cwavef(:,:)     ! current c coefficient vector

    real(dp), allocatable :: gvnlc(:,:) ! intermediate storage for <g|V_nl|c>

    real(dp) gsc_dummy(1,1) ! Dummy parameter, would be needed
    ! for PAW calculations if overlap <G|S|C> is desired

    type(cprj_type) :: cprj_dum(1,1) ! dummy argument that is
    ! required for doing PAW

    real(dp) tmp

    integer indx_k, indx_bk, i, save_paral_compil_fft

    integer wfopta10

    maxocc = maxval(occ)

    wfopta10=mod(wfoptalg,10)

    write(2244,*) "Inside orbital-shifts"

    sij_opt=0  ! this must be set to zero, and is PAW only

    nblock=(nband-1)/nbdblock+1

    ! Reserve space for Orbital shifts 
    ! vxc (the KS potential)
    !  and the orb_shift vector g
    allocate(&
         gvnlc(2,npw*nspinor), &           ! intermediate storage for <g|V_nl|c>
         orb_shift(2,npw*nspinor), &       ! the orbital shifts
         cg_b(2,npw*nspinor), &            ! k-space right hand side in CG
         potloc(gs_hamk%nfft,nvloc), &     ! real-space tmp potential
         rr_b(gs_hamk%nfft,nspinor), &     ! real-space realpart RHS in CG
         ri_b(gs_hamk%nfft,nspinor), &     ! real-space imagpart RHS in CG
         cg_resid(2,npw*nspinor), &        ! residual in CG 
         cg_p(2,npw*nspinor), &            ! p-vektor CG
         cg_Ap(2,npw*nspinor),&            ! A . p vektor in CG
         cwavef(2,npw*nspinor),&           ! k-space KS orbital
         wf_real(gs_hamk%nfft,nspinor),&   ! real-space KS orbital, realpart
         wf_imag(gs_hamk%nfft,nspinor))    ! real-space KS orbital, imagpart

    write(2244,*) "Entering nblock/nband loop, nband = ", nband
    !Loop over blocks of bands.
    !In the standard band-sequential algorithm, nblock=nband.

    !obtain access index for occ, eigen:
    indx_k = 0
    if (isppol == 2) then
       do i=1,nkpt
          indx_k = indx_k + run_params%nband(i) 
       enddo
    endif
    do i=1,ikpt-1
       indx_k = indx_k + run_params%nband(i + (isppol-1)*nkpt) 
    enddo

! Need to do this to prevent abinit from hanging 
! in the dotprod_g calls: there, if paral_compil_fft
! is 1, abinit wants to exchange dotprod_g results - WHY?
!    save_paral_compil_fft = mpi_enreg%paral_compil_fft
!    mpi_enreg%paral_compil_fft = 0

    do iblock=1,nblock

       !  Loop over bands in a block
       !  This loop can be MPI-parallelized, 
       !  over processors attached to the same k point
       ibandmin=1+(iblock-1)*nbdblock
       ibandmax=min(iblock*nbdblock,nband)

       !  =====================================================================================
       ! now comes the loop over all bands
       ! for each band, solve for the orbital shift
       ! and update the dE/dV gradient
       ! we can do this in parallel,
       ! we just have to collect all the calculated dEdU gradients
       ! at the end

       do iband=ibandmin,ibandmax

          indx_bk = indx_k + iband 

          ibdblock=iband-(iblock-1)*nbdblock
          ! index for the wave-function coefficient


          !    In MPI-parallelisation, one should determine if the present
          !    band is treated by the present processor ...
#if defined HAVE_MPI

          write(2244,*) "Check whether I am responsible for",iband
          write(2244,*) "mpi_enreg%paralbd",mpi_enreg%paralbd
          write(2244,*) "mpi_enreg%proc_distrb",mpi_enreg%proc_distrb
          if (mpi_enreg%paralbd >= 1) then
             if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/= mpi_enreg%me)&
                  &         then
                if (mpi_enreg%me_kpt/=0) then
                   write(2244,*) "No"
                   cycle
                else
                   write(2244,*) "Yes"
                end if
             end if
          end if

#endif

          write(2244,*) "Occupation of band ",iband," is ",&
               occ(indx_bk)

!!          we do not need to treat empty bands

          if (occ(indx_bk) < 1d-5) then
             cycle
          endif

          write(2244,*) "Set up orbital shift equations"
          ! 1. Set up orbital shift equations:

          write(2244,*) "NPW ",npw, "gs_hamk%nfft", gs_hamk%nfft

          ! INITIALIZE CG 

          ntot_k = npw*nspinor

          orb_shift = 0.0

          write(2244,*) "original wavefunction"
          !    Extraction of the original wavefunction vector
          
          icg_shift=npw*nspinor*(iband-1)+icg

          do ipw=1,ntot_k
             cwavef(:,ipw)=cg(:,ipw+icg_shift)
          end do

          ! Obtain wf in real space
          write(2244,*) "fourier transform k->r"

          write(2244,*) "max/min cwavef",&
               minval(cwavef),maxval(cwavef)

          call recip_to_cmplx(mpi_enreg,paral_kgb,nkpt,npw, ikpt,&
               gs_hamk%istwf_k,nspinor,mgfft,npw,gs_hamk%ucvol,gs_hamk%ngfft,&
               kg_k,npwarr(ikpt),cwavef,wf_real,wf_imag)

          write(2244,*) "done"
          write(2244,*) "max/min wfreal",&
               minval(wf_real),maxval(wf_real)
          write(2244,*) "max/min wf_imag",&
               minval(wf_imag),maxval(wf_imag)

          call dotprod_g(dotr,doti,gs_hamk%istwf_k,&
               mpi_enreg,ntot_k,2,cwavef,cwavef)

          write(2244,*) "norm cwavef",dotr, doti

          ! Calculate (Vxc - Uxc)*psi in real space
          
          ! (1) get the uxc potential, which will in general
          ! depend on the k-point and band

          tmp = run_params%wtk(ikpt) * run_params%gradscal * &
               occ(indx_bk) * gs_hamk%ucvol / dble(gs_hamk%nfft)
          
          if (nspinor==1) then             
             call get_uxc(wf_real(:,1), rr_b(:,1), wf_imag(:,1), ri_b(:,1),&
                  & gs_hamk%nfft, isppol, iband, ikpt,&
                  & run_params)
             
!             potloc(:,1) = (run_params%orb_shift_vxc(:,isppol)&
!                  - potloc(:,1))*tmp

             potloc(:,1) = (run_params%orb_shift_vxc(:,isppol)&
                  - run_params%vhrtr_prefac*vhartree_old(:))*tmp
          !        right-hand side of orbital shift equation
          !        b is difference of uxc and vxc
             ri_b(:,1) = potloc(:,1) * wf_imag(:,1) - ri_b(:,1)*tmp
             rr_b(:,1) = potloc(:,1) * wf_real(:,1) - rr_b(:,1)*tmp
          else
             ! (spinor = 2)
             if (nvloc==4) then ! (full 2x2 problem)
!!     V are stored as : V^11, V^22, Re[V^12], Im[V^12] (complex, hermitian)
stop
!                call get_uxc(potloc(:,1), gs_hamk%nfft, 1, iband, &
!                     &ikpt, run_params)
!                call get_uxc(potloc(:,2), gs_hamk%nfft, 2, iband, &
!                     &ikpt, run_params)
!                call get_uxc(potloc(:,3), gs_hamk%nfft, 3, iband, &
!                     &ikpt, run_params)
!                call get_uxc(potloc(:,4), gs_hamk%nfft, 4, iband, &
!                     &ikpt, run_params)

                potloc(:,:) = (run_params%orb_shift_vxc(:,:) - potloc(:,:))*tmp

                rr_b(:,1) = &
                     potloc(:,1) * wf_real(:,1) + &
                     potloc(:,3) * wf_real(:,2) - &
                     potloc(:,4) * wf_imag(:,2) 
                ri_b(:,1) = &
                     potloc(:,1) * wf_imag(:,1) + &
                     potloc(:,3) * wf_imag(:,2) + &
                     potloc(:,4) * wf_real(:,2) 
                rr_b(:,2) = &
                     potloc(:,2) * wf_real(:,2) + &
                     potloc(:,3) * wf_real(:,1) + &
                     potloc(:,4) * wf_imag(:,1) 
                ri_b(:,2) = &
                     potloc(:,2) * wf_imag(:,2) + &
                     potloc(:,3) * wf_imag(:,1) - &
                     potloc(:,4) * wf_real(:,1) 

             else !spin-orbit
                stop
!                call get_uxc(potloc(:,1), gs_hamk%nfft, 1, iband, &
!                     &ikpt, run_params)

                potloc(:,1) = (run_params%orb_shift_vxc(:,1) - potloc(:,1))*tmp

                rr_b(:,1) = potloc(:,1) * wf_real(:,1) 
                ri_b(:,1) = potloc(:,1) * wf_imag(:,1) 
                rr_b(:,2) = potloc(:,1) * wf_real(:,2) 
                ri_b(:,2) = potloc(:,1) * wf_imag(:,2) 
             endif
          endif

          write(2244,*) "max/min rr_b",&
               minval(rr_b),maxval(rr_b)
          write(2244,*) "max/min ri_b",&
               minval(ri_b),maxval(ri_b)

          write(2244,*) "fourier transform r->k"
          ! convert rr_b to reciprocal space

!    mpi_enreg%paral_compil_fft = save_paral_compil_fft

          call cmplx_to_recip(mpi_enreg,paral_kgb,nkpt,npw,ikpt,&
             gs_hamk%istwf_k, nspinor,mgfft,gs_hamk%nfft,gs_hamk%ucvol,&
             gs_hamk%ngfft, kg_k,npwarr(ikpt),rr_b,ri_b,cg_b)  

!    save_paral_compil_fft = mpi_enreg%paral_compil_fft
!    mpi_enreg%paral_compil_fft = 0

          ! project out the cwavef part

          call orthogonalize( cg_b, cwavef, gs_hamk%istwf_k, &
               mpi_enreg, ntot_k)
          !! scale
          cg_resid = cg_b

          cg_p      = cg_resid  

          write(2244,*) "max/min cg_resid",&
               minval(cg_resid),maxval(cg_resid)

          call dotprod_g(dnorm,doti,gs_hamk%istwf_k,&
               mpi_enreg,ntot_k,1,cg_resid,cg_resid)

          normb = dnorm          

          write(2244,*) "entering CG loop for ",&
               eigen(indx_bk)
          !CG LOOP
          do icg_step = 1, 1000

             if (dnorm < run_params%CG_accuracy * normb) then 
                write(2244,*)&
                     'orbital shifts CG converged to ',dnorm/normb
                exit
             endif

             ! now compute <g|H|cg_p>
             cptopt=0
             call getghc(cptopt,cg_p,cprj_dum,dimffnl,ffnl,filstat,&
                  cg_Ap,gsc_dummy,gs_hamk,gvnlc,kg_k,&
                  kinpw,lambda,mpi_enreg,natom,1,npw,nspinor,&
                  paral_kgb,ph3d,prtvol,sij_opt,0,0,vlocal)

             ! <ei|cg_p>
             call dotprod_g(dotr,doti,gs_hamk%istwf_k,&
                  mpi_enreg,ntot_k,2,cwavef,cg_p)

             ! obtain the complete A|cg_p> = (H - ei + alpha |ei><ei|)|cg_p>
             do ipw=1,ntot_k
                cg_Ap(1,ipw) = cg_Ap(1,ipw) - eigen(indx_bk) * &
                     cg_p(1,ipw) + cwavef(1,ipw) * dotr - cwavef(2,ipw) * doti
                cg_Ap(2,ipw) = cg_Ap(2,ipw) - eigen(indx_bk) * &
                     cg_p(2,ipw) + cwavef(2,ipw) * dotr + cwavef(1,ipw) * doti
             enddo

!             call orthogonalize( cg_Ap, cwavef, gs_hamk%istwf_k, &
!                  mpi_enreg, ntot_k)

             call dotprod_g(dotr,doti,gs_hamk%istwf_k,&
                  mpi_enreg,ntot_k,2,cg_p,cg_Ap)

             zg_alpha = dnorm / dcmplx(dotr,doti)
             ! zg_alpha = < cg_resid | cg_resid > / < cg_p | cg_Ap >

             do ipw=1,ntot_k
                orb_shift(1,ipw) = orb_shift(1,ipw) +  &
                     dreal(zg_alpha) * cg_p(1,ipw) - &
                     dimag(zg_alpha) * cg_p(2,ipw)  
                orb_shift(2,ipw) = orb_shift(2,ipw) +  &
                     dimag(zg_alpha) * cg_p(1,ipw) + &
                     dreal(zg_alpha) * cg_p(2,ipw) 
                ! orb_shift = orb_shift + cg_alpha * cg_p

                cg_resid(1,ipw)  = cg_resid(1,ipw) -  &
                     dreal(zg_alpha) * cg_Ap(1,ipw) + &
                     dimag(zg_alpha) * cg_Ap(2,ipw) 
                cg_resid(2,ipw)  = cg_resid(2,ipw) -  &
                     dimag(zg_alpha) * cg_Ap(1,ipw) - &
                     dreal(zg_alpha) * cg_Ap(2,ipw) 
                ! cg_resid = cg_resid - cg_alpha * cg_Ap

             enddo

             dnorm_old = dnorm

             call dotprod_g(dnorm,doti,gs_hamk%istwf_k,&
                  mpi_enreg,ntot_k,1,cg_resid,cg_resid)
             
             cg_beta = dnorm / dnorm_old

             cg_p    = cg_resid + cg_beta * cg_p

          enddo
          !END CG LOOP

          call dotprod_g(dotr,doti,gs_hamk%istwf_k,&
              mpi_enreg,ntot_k,1,orb_shift,orb_shift)          
          write(2244,*) "orbital shift in recipspc(max/min/norm) ",&
               maxval(orb_shift), minval(orb_shift), dotr

          ! project out the part collinear to |c>
          call orthogonalize( orb_shift, cwavef, gs_hamk%istwf_k, &
               mpi_enreg, ntot_k)

          call dotprod_g(dotr,doti,gs_hamk%istwf_k,&
              mpi_enreg,ntot_k,1,orb_shift,orb_shift)          
          write(2244,*) "orbital shift after projection(max/min/norm) ",&
               maxval(orb_shift), minval(orb_shift), dotr

          ! Check solution

          call getghc(cptopt,orb_shift,cprj_dum,dimffnl,ffnl,filstat,&
               cg_Ap,gsc_dummy,gs_hamk,gvnlc,kg_k,&
               kinpw,lambda,mpi_enreg,natom,1,npw,nspinor,&
               paral_kgb,ph3d,prtvol,sij_opt,0,0,vlocal)

          do ipw=1,ntot_k
             cg_Ap(1,ipw) = cg_Ap(1,ipw) - eigen(indx_bk) * &
                  orb_shift(1,ipw) - cg_b(1,ipw) 
             cg_Ap(2,ipw) = cg_Ap(2,ipw) - eigen(indx_bk) * &
                  orb_shift(2,ipw) - cg_b(2,ipw)
          enddo

          call dotprod_g(dotr,doti,gs_hamk%istwf_k,&
              mpi_enreg,ntot_k,1,cg_Ap, cg_Ap)          

          write(2244,*) "Solution for CG ORB EQ correct by", dsqrt(dotr)
             
          ! end Check solution

          call recip_to_cmplx(mpi_enreg,paral_kgb,nkpt,npw,&
               ikpt, gs_hamk%istwf_k,nspinor,mgfft,npw,gs_hamk%ucvol,&
               gs_hamk%ngfft, kg_k,npwarr(ikpt),orb_shift,rr_b,ri_b)

          write(2244,*) "updating Gradient with rrb (max/min) ",&
               maxval(rr_b), minval(rr_b)
          write(2244,*) "updating Gradient with rib (max/min) ",&
               maxval(ri_b), minval(ri_b)

          write(2244,*) "KPT weight ", run_params%wtk(ikpt)

          ! Update total energy gradient
          ! *2 for H.c.
          run_params%dEdu(:,isppol) = run_params%dEdu(:,isppol) + & 
                (rr_b(:,1) * wf_real(:,1) + ri_b(:,1)*wf_imag(:,1)) * 2.d0


          if (nspinor==2) then
             if (nvloc==4) then
             ! V22
                run_params%dEdu(:,2) = run_params%dEdu(:,2) + & 
                     (rr_b(:,2) * wf_real(:,2) + ri_b(:,2)*wf_imag(:,2)) * 2.d0

                ! Re V12
                run_params%dEdu(:,3) = run_params%dEdu(:,3) + & 
                     (rr_b(:,1) * wf_real(:,2) + ri_b(:,1)*wf_imag(:,2) + &
                     rr_b(:,2) * wf_real(:,1) + ri_b(:,2)*wf_imag(:,1))
                ! Im V12, not sure about sign here
                run_params%dEdu(:,4) = run_params%dEdu(:,4) + & 
                     (rr_b(:,1) * wf_imag(:,2) - rr_b(:,2)*wf_imag(:,1) + &
                     ri_b(:,2) * wf_real(:,1) - ri_b(:,1)*wf_real(:,2))

             else !spin-orbit
                run_params%dEdu(:,isppol) = run_params%dEdu(:,isppol) + & 
                     (rr_b(:,2) * wf_real(:,2) + ri_b(:,2)*wf_imag(:,2)) * 2.d0                
             endif
          endif

       enddo
       !END LOOP OVER BANDS
    enddo
    !END LOOP OVER BAND BLOCKS

    write(2244,*) "Finnished with the orbital shift stuff, cleanup"

!    mpi_enreg%paral_compil_fft = save_paral_compil_fft

    deallocate(gvnlc, cg_b, orb_shift, rr_b, ri_b, cg_resid, cg_p, cg_Ap, &
         cwavef, wf_real, wf_imag, potloc)    

    write(2244,*) "Finnished with cleanup, exit orbital_shifts now"

  end subroutine cgwf_orbital_shift
  !!***
  
end module orbital_shifts
