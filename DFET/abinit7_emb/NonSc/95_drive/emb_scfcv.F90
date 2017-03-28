!{\src2tex{textfont=tt}}
!!****f* ABINIT/scfcv
!! NAME
!! scfcv
!!
!! FUNCTION
!! Self-consistent-field convergence.
!! Conducts set of passes or overall iterations of preconditioned
!! conjugate gradient algorithm to converge wavefunctions to
!! ground state and optionally to compute forces and energy.
!! This routine is called to compute forces for given atomic
!! positions or else to do non-SCF band structures.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (XG, GMR, AR, MKV, MT, FJ, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cpus= cpu time limit in seconds
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs for the "coarse" grid (see NOTES below)
!!   | mkmem =number of k points which can fit in memory; set to 0 if use disk
!!   | mpw=maximum dimensioned size of npw.
!!   | natom=number of atoms in cell.
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   |      for the "coarse" grid (see NOTES below)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  ecore=core psp energy (part of total energy) (hartree)
!!  fatvshift=factor to multiply dtset%atvshift
!!  iapp=indicates the eventual suffix to be appended to the generic output root
!!         if 0 : no suffix to be appended (called directly from gstate)
!!         if positive : append "_TIM//iapp" (called from move or brdmin)
!!         if -1 : append "_TIM0" (called from brdmin)
!!         if -2, -3, -4, -5: append "_TIMA", ... ,"_TIMD", (called from move)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)= # atoms of each type.
!!  ndtpawuj=size of dtpawuj
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points and spins
!!
!! SIDE EFFECTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=updated wavefunctions;  if mkmem>=nkpt,
!!         these are kept in a disk file.
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!  dtpawuj(ndtpawuj)= data used for the automatic determination of U (relevant only for PAW+U)
!!      calculations (see initberry.f)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  initialized= if 0 the initialization of the gstate run is not yet finished
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!       for the "fine" grid (see NOTES below)
!!  occ(mband*nkpt*nsppol)=occupation number for each band (often 2) at each k point
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!   (should be made a pure output quantity)
!!  rhog(2,nfftf)=array for Fourier transform of electron density
!!  rhor(nfftf,nspden)=array for electron density in el./bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
!!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
!!  wffnew,wffnow=struct info for wf disk files.
!!  wvl <type(wvl_data)>=all wavelets data.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)= at input, previous reduced dimensionless atomic coordinates
!!                     at output, current xred is transferred to xred_old
!!
!! NOTES
!! It is worth to explain THE USE OF FFT GRIDS:
!! ============================================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! PARENTS
!!      brdmin,delocint,diisrelax,gstate,moldyn,move,pawuj_drive,scphon
!!
!! CHILDREN
!!      abi_etsf_init,afterscfloop,berryphase_new,calc_xc_ep,chkdilatmx
!!      chkpawovlp,cprj_alloc,cprj_free,ctocprj,energies_init,energy,etotfor
!!      extraprho,fappnd,fourdp,fresid,getcut,getmpw,getng,getph,initylmg
!!      int2char4,ioarr,kpgio,leave_new,leave_test,metric,newrho,newvtr
!!      nhatgrid,odamix,out_geometry_xml,out_resultsgs_xml,outscfcv,pawdenpot
!!      pawdij,pawmknhat,prtene,rhohxc,rhotov,scprqt,setnoccmmp,setrhoijpbe0
!!      setsym,setup_positron,setvtr,sphereboundary,status,symdij,symrhg,symzat
!!      timab,vtorho,vtorhorec,vtorhotf,wrtout,wvl_mkrho,wvl_newvtr
!!      wvl_wfsinp_reformat,xcomm_world,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine emb_scfcv(atindx,atindx1,cg,cpus,dtefield,dtfil,dtpawuj,&
&  dtset,ecore,eigen,electronpositron,fatvshift,hdr,iapp,indsym,initialized,&
&  irrzon,kg,mpi_enreg,nattyp,ndtpawuj,nfftf,npwarr,occ,&
&  paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&  phnons,psps,pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
&  scf_history,symrec,taug,taur,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)

!=========== chen ====================
 use c_hc_vars
!=========== end of hack =============

 use defs_basis
 use defs_datatypes
 use m_wffile
 use defs_abitypes
 use defs_wvltypes
 use defs_parameters
 use defs_rectypes
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype
 use m_rec
 use m_io_tools, only : flush_unit
 use m_paw_dmft, only: paw_dmft_type
#if defined HAVE_ETSF_IO
 use etsf_io
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_12_hide_mpi
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_27_toolbox_oop
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_abiutil
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_56_xc
 use interfaces_57_iovars
 use interfaces_61_ionetcdf
 use interfaces_62_iowfdenpot
 use interfaces_65_nonlocal
 use interfaces_66_paw
 use interfaces_67_common
 use interfaces_68_recursion
 use interfaces_68_rsprc
 use interfaces_79_seqpar_mpi
 use interfaces_95_drive, except_this_one => scfcv
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iapp,ndtpawuj,pwind_alloc
 integer,intent(inout) :: initialized,nfftf
 real(dp),intent(in) :: cpus,ecore,fatvshift
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(efield_type),intent(inout) :: dtefield
 type(electronpositron_type),pointer :: electronpositron
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(inout) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(recursion_type),intent(inout) :: rec_set
 type(results_gs_type),intent(inout) :: results_gs
 type(scf_history_type),intent(inout) :: scf_history
 type(wffile_type),intent(inout) :: wffnew,wffnow
 type(wvl_data),intent(inout) :: wvl
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
 integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
!no_abirules
 integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  !(nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise)
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
 integer, intent(in) :: nattyp(psps%ntypat),npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
 integer, intent(in) :: symrec(3,3,dtset%nsym)
 real(dp), intent(inout) :: cg(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
 real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  !(nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise)
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
 real(dp), intent(in) :: rprimd(3,3)
 real(dp), pointer :: rhog(:,:),rhor(:,:)
 real(dp), pointer :: taug(:,:),taur(:,:)
 real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: xred(3,dtset%natom),xred_old(3,dtset%natom)
 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(macro_uj_type),intent(inout) :: dtpawuj(0:ndtpawuj)
 type(pawrhoij_type), intent(inout) :: pawrhoij(mpi_enreg%natom*psps%usepaw)
 type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(paw_dmft_type), intent(inout) :: paw_dmft

!Local variables -------------------------
!Variables for partial dos calculation
!scalars
 integer,parameter :: level=110,response=0
 integer :: accessfil,afford,choice
 integer :: computed_forces,cplex,dbl_nnsclo,dielop,dielstrt,dimdmat
 integer :: fformr,forces_needed
!integer :: dtset_iprcel
 integer :: iatom,ider,idir,ierr,iexit,ii,ikpt,impose_dmat
 integer :: initialized0,iorder_cprj,ipert,ipositron,ir,isave,iscf10,ispden
 integer :: ispmix,istep,istep_mix,itypat,izero
 integer :: lm_size,lmax_diel,lmn2_size,lpawumax
 integer :: mgfftdiel,mgfftf,moved_atm_inside,moved_rhor,n1xccc
 integer :: n3xccc,n_fftgr,n_index,nele,nfftdiel,nfftmix,nhatgrdim,nk3xc,nkxc
 integer :: npawmix,npwdiel,nstep,nzlmopt,offset,optberry,optene,optgr0
 integer :: optgr1,optgr2,option,optrad,optres,optxc,prtden,prtkden,prtfor,quit
 integer :: quit_sum,rdwr,rdwrpaw,spaceComm,stress_needed,unit_out
 integer :: usecprj,usexcnhat,v_size
 real(dp) :: boxcut,compch_fft,compch_sph,deltae,diecut,diffor,ecut
 real(dp) :: ecutf,ecutsus,edum,elast,etotal,fermie,gsqcut
 real(dp) :: maxfor,res2,residm,ucvol,val_max
 real(dp) :: val_min,vxcavg,vxcavg_dum
 !logical :: ex
 character(len=4) :: tag
 character(len=500) :: message
 character(len=fnlen) :: fildata,kgnam
 type(MPI_type) :: mpi_enreg_diel
 type(energies_type) :: energies
!arrays
 integer :: dt_ngfft(3),ngfft(18),ngfftdiel(18),ngfftf(18),ngfftmix(18),npwarr_diel(1)
 integer :: npwtot_diel(1)
 integer,allocatable :: dimcprj(:),gbound_diel(:,:),i_rhor(:),i_vresid(:)
 integer,allocatable :: i_vrespc(:),i_vtrial(:),irrzondiel(:,:,:),kg_diel(:,:)
 integer,allocatable :: lmn_size(:)
 integer,allocatable :: indsym_dum(:,:,:),symrec_dum(:,:,:)
 real(dp) :: dielar(7),dphase(3),dummy2(6),favg(3),gmet(3,3),gprimd(3,3),k0(3)
 real(dp) :: kpt_diel(3),pel(3),pel_cg(3),pelev_dum(3),pion(3),ptot(3)
 real(dp) :: rhodum(1),rmet(3,3),strsxc(6),strten(6),tollist(12)
 real(dp) :: tsec(2),vnew_mean(dtset%nspden),vres_mean(dtset%nspden)
 real(dp),allocatable :: dielinv(:,:,:,:,:),dtn_pc(:,:)
 real(dp),allocatable :: f_atm(:,:,:),f_fftgr(:,:,:),f_paw(:,:)
 real(dp),allocatable :: fcart(:,:),forold(:,:),fred(:,:),gresid(:,:)
 real(dp),allocatable :: grewtn(:,:),grhf(:,:),grnl(:),grxc(:,:)
 real(dp),allocatable :: kxc(:,:),nhat(:,:),nhatgr(:,:,:),nvresid(:,:)
 real(dp),allocatable :: ph1d(:,:),ph1ddiel(:,:),ph1df(:,:)
 real(dp),allocatable :: phnonsdiel(:,:,:),shiftvector(:)
 real(dp),allocatable :: susmat(:,:,:,:,:),synlgr(:,:)
 real(dp),allocatable :: vhartr(:),vpsp(:),vtrial(:,:)
 real(dp),allocatable :: vxc(:,:),workr(:,:),xccc3d(:),ylmdiel(:,:)
 real(dp),pointer :: elfr(:,:),grhor(:,:,:),lrhor(:,:)
 type(cprj_type),allocatable :: cprj(:,:)
 type(paw_an_type),allocatable :: paw_an(:)
 type(paw_ij_type),allocatable :: paw_ij(:)
 type(pawfgrtab_type),allocatable,save :: pawfgrtab(:)


 !============================ chen ===================
 integer:: &
   lmn, &
   c_itype
 real(dp):: &
   tmp1, tmp2, tmp3, &
   xcart(3,dtset%natom),intg
 real(dp),allocatable::  & 
   extpot_lm(:,:,:), &    ! (l,m) moments of extpot at each atom, on a radial mesh
   extpot_dij(:,:)        ! store dij due to \int[v_emb*n1]

 !=====================================================

! *********************************************************************

!DEBUG
!write(6,*)' scfcv : enter'
!ENDDEBUG

#if defined DEBUG_MODE
 write(message,'(a)')' scfcv : enter '
 call wrtout(std_out,message,'COLL')
 call flush_unit(std_out)
#endif

 call timab(20,1,tsec)
 call timab(54,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 if(dtset%prtvol==-level)then
   write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' scfcv : enter '
   call wrtout(std_out,message,'COLL')
 end if

!######################################################################
!Initializations - Memory allocations
!----------------------------------------------------------------------

 call status(0,dtfil%filstat,iexit,level,'allocate/init ')

 dielstrt=0

!Save some variables from dataset definition
 nstep=dtset%nstep
!dtset_iprcel = dtset%iprcel
 ecut=dtset%ecut
 ecutf=ecut;if (psps%usepaw==1.and.pawfgr%usefinegrid==1) ecutf=dtset%pawecutdg
 iscf10=mod(dtset%iscf,10)
 tollist(1)=dtset%tolmxf;tollist(2)=dtset%tolwfr
 tollist(3)=dtset%toldff;tollist(4)=dtset%toldfe
 tollist(6)=dtset%tolvrs;tollist(7)=dtset%tolrff

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Some variables need to be initialized/nullify at start
 quit=0 ; dbl_nnsclo=0 ;
 dielop=0 ; strsxc=zero
 deltae=zero ; elast=zero ;
 if (dtset%positron==0.or.initialized==0) then
   call energies_init(energies)
 else
   energies%e0_electronpositron =results_gs%energies%e0_electronpositron
   energies%e_electronpositron  =results_gs%energies%e_electronpositron
   energies%edc_electronpositron=results_gs%energies%edc_electronpositron
 end if
 if (dtset%nstep==0) energies%e_fermie=results_gs%energies%e_fermie
 energies%e_corepsp = ecore / ucvol
 fermie=energies%e_fermie
 isave=0 !initial index of density protection file
 optres=merge(0,1,dtset%iscf<10)
 usexcnhat=0;usecprj=0
 initialized0=initialized
 ipert=0;idir=0;cplex=1
 istep_mix=1
 ipositron=electronpositron_calctype(electronpositron)

!Stresses and forces flags
 forces_needed=0;prtfor=0
 if ((dtset%optforces==1.or.dtset%ionmov==4.or.abs(tollist(3))>tiny(0._dp))) then
   if (dtset%iscf>0.and.nstep>0) forces_needed=1
   if (nstep==0) forces_needed=2
   prtfor=1
 else if (dtset%iscf>0.and.dtset%optforces==2) then
   forces_needed=2
 end if
 stress_needed=0
 if (dtset%optstress>0.and.dtset%iscf>0.and.dtset%prtstm==0.and. &
& (nstep>0.or.dtfil%ireadwf==1)) stress_needed=1

!This is only needed for the tddft routine, and does not
!correspond to the intented use of results_gs (should be only
!for output of scfcv
 etotal  =results_gs%etotal

!Get FFT grid(s) sizes (be careful !)
!See NOTES in the comments at the beginning of this file.
 ngfft(:)=dtset%ngfft(:)
 if (psps%usepaw==1) then
   mgfftf=pawfgr%mgfft;ngfftf(:)=pawfgr%ngfft(:)
 else
   mgfftf=dtset%mgfft;ngfftf(:)=ngfft(:)
 end if

!We create Output files when required
 if (dtset%accesswff == 3) then
#if defined HAVE_ETSF_IO
!  Compute this lmn_size stuff
   allocate(lmn_size(psps%npsp))
   if(psps%usepaw==1) then
     lmn_size(:) = pawtab(1:psps%npsp)%lmn_size
!    DC: in PAW, the scalar quantities like the densities are on the fine grid.
     dt_ngfft(:)=ngfftf(1:3)
   else
     lmn_size(:) = psps%lmnmax
     dt_ngfft(:)=dtset%ngfft(1:3)
   end if
!  Create an ETSF file for each required files
   if (dtset%prtden /= 0) then
!    Case of density.
     call abi_etsf_init(dt_ngfft,dtset,dtfil%fnameabo_app_den, 1, .true., lmn_size, psps, wvl%wfs)
   end if
   if (dtset%prtelf /= 0) then
!    Case of electron localization function
     call abi_etsf_init(dt_ngfft,dtset,dtfil%fnameabo_app_elf, 1, .true., lmn_size, psps, wvl%wfs)
     if (dtset%nspden==2) then
!      Case of spin-dependent electron localization function
       call abi_etsf_init(dt_ngfft,dtset,dtfil%fnameabo_app_elf_up, 1, .true., lmn_size, psps, wvl%wfs)
       call abi_etsf_init(dt_ngfft,dtset,dtfil%fnameabo_app_elf_down, 1, .true., lmn_size, psps, wvl%wfs)
     end if
   end if
   if (dtset%prtgden /= 0) then
!    Case of gradient of electron density.
     call abi_etsf_init(dt_ngfft,dtset, dtfil%fnameabo_app_gden1, 1, .true., lmn_size, psps, wvl%wfs)
     call abi_etsf_init(dt_ngfft,dtset, dtfil%fnameabo_app_gden2, 1, .true., lmn_size, psps, wvl%wfs)
     call abi_etsf_init(dt_ngfft,dtset, dtfil%fnameabo_app_gden3, 1, .true., lmn_size, psps, wvl%wfs)
   end if
   if (dtset%prtkden /= 0) then
!    Case of kinetic energy density.
     call abi_etsf_init(dt_ngfft,dtset,dtfil%fnameabo_app_kden, 1, .true., lmn_size, psps, wvl%wfs)
   end if
   if (dtset%prtlden /= 0) then
!    Case of Laplacian of electron density.
     call abi_etsf_init(dt_ngfft,dtset,dtfil%fnameabo_app_lden, 1, .true., lmn_size, psps, wvl%wfs)
   end if
   if (dtset%prtwf == 1) then
!    Case of wavefunctions.
     call abi_etsf_init(dt_ngfft,dtset, dtfil%fnameabo_app_wfk, 2, .true., lmn_size, psps, wvl%wfs)
   end if
   if (dtset%prtvxc /= 0) then
!    Case of Exchange-correlation potential.
     call abi_etsf_init(dt_ngfft,dtset,dtfil%fnameabo_app_vxc, 24, .true., lmn_size, psps, wvl%wfs)
   end if
!  FIXME: append other possibilities.
!  * ETSFIO cases of only correlation or only exchange are not abinit
!  options
!  * the fix for VHA,VHXC,POT,STM is dirty: they are flagged as
!  exchange correlation pot files, except stm, which is flagged density.
!  START dirty treatment
   if (dtset%prtvha /= 0) then
!    Case of Hartree potential.
     call abi_etsf_init(dt_ngfft,dtset,dtfil%fnameabo_app_vha, 24, .true., lmn_size, psps, wvl%wfs)
   end if
   if (dtset%prtvhxc /= 0) then
!    Case of Hartree+XC potential.
     call abi_etsf_init(dt_ngfft,dtset,dtfil%fnameabo_app_vhxc, 24, .true., lmn_size, psps, wvl%wfs)
   end if
   if (dtset%prtpot /= 0) then
!    Case of total potential.
     call abi_etsf_init(dt_ngfft,dtset,dtfil%fnameabo_app_pot, 24, .true., lmn_size, psps, wvl%wfs)
   end if
   if (dtset%prtstm /= 0) then
!    Case of STM output.
     call abi_etsf_init(dt_ngfft,dtset,dtfil%fnameabo_app_stm, 1, .true., lmn_size, psps, wvl%wfs)
   end if
   deallocate(lmn_size)
#endif
 end if

!Entering a scfcv loop, printing data to XML file if required.
 if (mpi_enreg%me == 0 .and. dtset%prtxml == 1) then
!  scfcv() will handle a scf loop, so we output the scfcv markup.
   write(ab_xml_out, "(A)") '    <scfcvLoop>'
   write(ab_xml_out, "(A)") '      <initialConditions>'
!  We output the geometry of the dataset given in argument.
!  xred and rprimd are given independently since dtset only
!  stores original and final values.
   call out_geometry_XML(dtset, 4, dtset%natom, rprimd, xred)
   write(ab_xml_out, "(A)") '      </initialConditions>'
 end if

!Examine tolerance criteria, and eventually  print a line to the output
!file (with choice=1, the only non-dummy arguments of scprqt are
!nstep, tollist and iscf - still, diffor and res2 are here initialized to 0)
 choice=1 ; diffor=zero ; res2=zero
 allocate(fcart(3,dtset%natom),fred(3,dtset%natom))
 fred(:,:)=zero
 fcart(:,:)=results_gs%fcart(:,:) ! This is a side effect ...
!results_gs should not be used as input of scfcv
 call scprqt(choice,cpus,deltae,diffor,dtset,&
& eigen,etotal,favg,fcart,energies%e_fermie,dtfil%fnameabo_app_eig,dtfil%filnam_ds(1),&
& initialized0,dtset%iscf,istep,dtset%kptns,maxfor,moved_atm_inside,mpi_enreg,&
& dtset%nband,dtset%nkpt,nstep,occ,optres,&
& prtfor,quit,res2,resid,residm,response,tollist,psps%usepaw,&
& vxcavg,dtset%wtk,xred)

!Various allocations (potentials, gradients, ...)
 allocate(forold(3,dtset%natom),grnl(3*dtset%natom),gresid(3,dtset%natom),&
& grewtn(3,dtset%natom),grxc(3,dtset%natom),synlgr(3,dtset%natom))
 allocate(ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom),ph1df(2,3*(2*mgfftf+1)*dtset%natom))
 allocate(vhartr(nfftf),vtrial(nfftf,dtset%nspden),vpsp(nfftf),vxc(nfftf,dtset%nspden))
 forold(:,:)=zero ; gresid(:,:)=zero ; pel(:)=zero
 n1xccc=0;if (psps%n1xccc/=0) n1xccc=psps%n1xccc
 n3xccc=0;if (psps%n1xccc/=0) n3xccc=nfftf
 allocate(xccc3d(n3xccc))

!Allocations/initializations for PAW only
 lpawumax=-1
 if(psps%usepaw==1) then

!  Variables/arrays related to the fine FFT grid
   allocate(nhat(nfftf,dtset%nspden));if (nstep==0) nhat=zero
   allocate(pawfgrtab(dtset%natom))
   do iatom=1,dtset%natom
     pawfgrtab(iatom)%cplex=cplex
     pawfgrtab(iatom)%nspden=pawrhoij(iatom)%nspden
     pawfgrtab(iatom)%l_size=pawtab(dtset%typat(iatom))%lcut_size
     pawfgrtab(iatom)%nfgd=0;allocate(pawfgrtab(iatom)%ifftsph(0))
     pawfgrtab(iatom)%rfgd_allocated=0;allocate(pawfgrtab(iatom)%rfgd(0,0))
     pawfgrtab(iatom)%gylm_allocated=0;allocate(pawfgrtab(iatom)%gylm(0,0))
     pawfgrtab(iatom)%gylmgr_allocated=0;allocate(pawfgrtab(iatom)%gylmgr(0,0,0))
     pawfgrtab(iatom)%gylmgr2_allocated=0;allocate(pawfgrtab(iatom)%gylmgr2(0,0,0))
     pawfgrtab(iatom)%nhatfr_allocated=0;allocate(pawfgrtab(iatom)%nhatfr(0,0))
     pawfgrtab(iatom)%vlocgr_allocated=0;allocate(pawfgrtab(iatom)%vlocgr(0,0))
   end do
   compch_fft=-1.d5
   usexcnhat=maxval(pawtab(:)%vlocopt)
   if (usexcnhat==0.and.dtset%ionmov==4.and.dtset%iscf<10) then
     write(message, '(a,a,a,a)' ) ch10,&
&     ' scfcv :  ERROR -',ch10,&
&     '  You cannot simultaneously use ionmov=4 and such a PAW psp file !'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  Variables/arrays related to the PAW spheres
   allocate(paw_ij(dtset%natom),paw_an(dtset%natom),dimcprj(dtset%natom))
   compch_sph=-1.d5
   do iatom=1,dtset%natom
     itypat=dtset%typat(iatom)
     lmn2_size=pawtab(itypat)%lmn2_size
     lm_size=pawtab(itypat)%lcut_size**2
     paw_an(iatom)%cplex        =cplex
     paw_an(iatom)%has_vxc      =0
     paw_an(iatom)%has_kxc      =0
     paw_an(iatom)%has_vxcval   =0
     paw_an(iatom)%has_vhartree =0
     paw_an(iatom)%angl_size=pawang%angl_size
     paw_an(iatom)%mesh_size=pawtab(itypat)%mesh_size
     paw_an(iatom)%nspden   =dtset%nspden
     paw_an(iatom)%lm_size  =lm_size
     allocate(paw_an(iatom)%lmselect(lm_size))
     paw_ij(iatom)%cplex    =cplex
     paw_ij(iatom)%cplex_dij=dtset%nspinor
     paw_ij(iatom)%nspden   =dtset%nspden
     paw_ij(iatom)%nsppol   =dtset%nsppol
     paw_ij(iatom)%lmn_size =pawtab(itypat)%lmn_size
     paw_ij(iatom)%lmn2_size=lmn2_size
     paw_ij(iatom)%ndij     =max(dtset%nspinor**2,paw_ij(iatom)%nspden)
     allocate(paw_ij(iatom)%dij(paw_ij(iatom)%cplex_dij*lmn2_size,paw_ij(iatom)%ndij))
     paw_ij(iatom)%has_dijxc=0;nullify(paw_ij(iatom)%dijxc)
     paw_ij(iatom)%has_dijxc_val=0;nullify(paw_ij(iatom)%dijxc_val)
     paw_ij(iatom)%has_dijU=0;nullify(paw_ij(iatom)%dijU)
     paw_ij(iatom)%has_dijfr=0;nullify(paw_ij(iatom)%dijfr)
     dimcprj(iatom)=pawtab(itypat)%lmn_size
     if (dtset%pawspnorb>0) then
       paw_ij(iatom)%has_dijso=1
       allocate(paw_ij(iatom)%dijso(paw_ij(iatom)%cplex_dij*lmn2_size,paw_ij(iatom)%ndij))
     else
       paw_ij(iatom)%has_dijso=0;nullify(paw_ij(iatom)%dijso)
     end if
     if (dtset%iscf==22) then
       paw_ij(iatom)%has_dijhat=1
       allocate(paw_ij(iatom)%dijhat(paw_ij(iatom)%cplex_dij*lmn2_size,paw_ij(iatom)%ndij))
     else
       paw_ij(iatom)%has_dijhat=0;nullify(paw_ij(iatom)%dijhat)
     end if
     if (pawtab(itypat)%usepawu>0) then
       lpawumax=max(pawtab(itypat)%lpawu,lpawumax)
       allocate(paw_ij(iatom)%noccmmp(paw_ij(iatom)%cplex_dij,2*pawtab(itypat)%lpawu+1,&
&       2*pawtab(itypat)%lpawu+1,paw_ij(iatom)%ndij))
       allocate(paw_ij(iatom)%nocctot(paw_ij(iatom)%ndij))
     end if
     if (pawtab(itypat)%useexexch>0) then
       allocate(paw_ij(iatom)%vpawx(1,lmn2_size,paw_ij(iatom)%nspden))
     end if
   end do

 end if ! PAW

!WVL - since wavelets change the size of the box, dont
!need dilatmax.
 if (dtset%usewvl == 0) then
!  Check that the possible change of unit cell size has not lead
!  to a too large increase
   call chkdilatmx(dtset%dilatmx,rprimd,dtset%rprimd_orig(1:3,1:3,1))
 end if


!Several parameters and arrays for the SCF mixing:
!These arrays are needed only in the self-consistent case
 if (dtset%iscf>0) then
   dielar(1)=dtset%diecut;dielar(2)=dtset%dielng
   dielar(3)=dtset%diemac;dielar(4)=dtset%diemix
   dielar(5)=dtset%diegap;dielar(6)=dtset%dielam
   dielar(7)=dtset%diemix;if (dtset%iscf>=10) dielar(7)=dtset%diemixmag
   allocate(nvresid(nfftf,dtset%nspden));if (nstep==0) nvresid=zero
   allocate(dtn_pc(3,dtset%natom))
   if(iscf10==1) then
!    For iscf==1, five additional vectors are needed
!    The index 1 is attributed to the old trial potential,
!    The new residual potential, and the new
!    preconditioned residual potential receive now a temporary index
!    The indices number 4 and 5 are attributed to work vectors.
     n_fftgr=5 ; n_index=1
     allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
     i_vtrial(1)=1 ; i_vresid(1)=2 ; i_vrespc(1)=3
   else if(iscf10==2) then
!    For iscf==2, three additional vectors are needed.
!    The index number 1 is attributed to the old trial vector
!    The new residual potential, and the new preconditioned
!    residual potential, receive now a temporary index.
     n_fftgr=3 ; n_index=1
     allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
     i_vtrial(1)=1 ; i_vresid(1)=2 ; i_vrespc(1)=3
   else if(iscf10==3) then
!    For iscf==3 , four additional vectors are needed.
!    The index number 1 is attributed to the old trial vector
!    The new residual potential, and the new and old preconditioned
!    residual potential, receive now a temporary index.
     n_fftgr=4 ; n_index=2
     allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
     i_vtrial(1)=1 ; i_vresid(1)=2 ; i_vrespc(1)=3 ; i_vrespc(2)=4
   else if (iscf10==4) then
!    For iscf==4 , six additional vectors are needed.
!    The indices number 1 and 2 are attributed to two old trial vectors
!    The new residual potential, and the new and two old preconditioned
!    residual potentials, receive now a temporary index.
     n_fftgr=6 ; n_index=3
     allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
     i_vtrial(1)=1 ; i_vtrial(2)=2 ; i_vresid(1)=3
     i_vrespc(1)=4 ; i_vrespc(2)=5 ; i_vrespc(3)=6
   else if((iscf10==5).or.(iscf10==6)) then
!    For iscf==5 or 6, ten additional vectors are needed
!    The index number 1 is attributed to the old trial vector
!    The index number 6 is attributed to the search vector
!    Other indices are attributed now. Altogether ten vectors
     n_fftgr=10 ; n_index=3
     allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
     i_vtrial(1)=1 ; i_vresid(1)=2 ; i_vrespc(1)=3 ; i_vresid(2)=4 ; i_vrespc(2)=5
     i_vresid(3)=7 ; i_vrespc(3)=8 ; i_rhor(2)=9 ; i_rhor(3)=10
   else if(iscf10==7) then
!    For iscf==7, lot of additional vectors are needed
!    The index number 1 is attributed to the old trial vector
!    The index number 2 is attributed to the old residual
!    The indices number 2 and 3 are attributed to two old precond. residuals
!    Other indices are attributed now.
     n_fftgr=2+2*dtset%npulayit ; n_index=1+dtset%npulayit
     allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
     i_vrespc(dtset%npulayit+1)=2*dtset%npulayit+1; i_vresid(1)=2*dtset%npulayit+2
     do ii=1,dtset%npulayit
       i_vtrial(ii)=2*ii-1 ; i_vrespc(ii)=2*ii
     end do
   end if ! iscf cases
!  The next arrays are needed if iscf==5 and ionmov==4,
!  but for the time being, they are always allocated
   allocate(grhf(3,dtset%natom),f_atm(3,dtset%natom,n_fftgr))
!  Additional allocation for mixing within PAW
   npawmix=0
   if(psps%usepaw==1) then
     do iatom=1,dtset%natom
       itypat=dtset%typat(iatom)
       lmn2_size=pawtab(itypat)%lmn2_size
       pawrhoij(iatom)%use_rhoijres=1
       allocate(pawrhoij(iatom)%rhoijres(pawrhoij(iatom)%cplex*lmn2_size,pawrhoij(iatom)%nspden))
       do ispden=1,pawrhoij(iatom)%nspden
         pawrhoij(iatom)%rhoijres(:,ispden)=zero
       end do
       allocate(pawrhoij(iatom)%kpawmix(pawtab(itypat)%lmnmix_sz))
       pawrhoij(iatom)%lmnmix_sz=pawtab(itypat)%lmnmix_sz
       pawrhoij(iatom)%kpawmix=pawtab(itypat)%kmix
       npawmix=npawmix+pawrhoij(iatom)%nspden*pawtab(itypat)%lmnmix_sz*pawrhoij(iatom)%cplex
     end do
   end if
 end if ! iscf>0

!Here, allocate arrays for computation of susceptibility and dielectric matrix or for TDDFT
 if( (nstep>0 .and. dtset%iscf>0) .or. dtset%iscf==-1 ) then !MF
!  The dielectric stuff is performed in sequential mode.
!  Set mpi_enreg_diel accordingly
   mpi_enreg_diel%paral_compil_fft=0
   mpi_enreg_diel%paral_compil_kpt=0
   mpi_enreg_diel%mode_para='n'
   mpi_enreg_diel%me=0
   mpi_enreg_diel%me_fft=0
   mpi_enreg_diel%me_kpt=0
   mpi_enreg_diel%nproc_fft=1
   mpi_enreg_diel%fft_option_lob=0
   mpi_enreg_diel%paral_fft=0
   nullify(mpi_enreg_diel%bandfft_kpt)
   nullify(mpi_enreg_diel%tab_kpt_distrib)
   mpi_enreg_diel%flag_ind_kg_mpi_to_seq = 0
!  Here, for TDDFT, artificially set iprcel . Also set a variable to reduce
!  the memory needs.
   afford=1
   if(dtset%iscf==-1) then
!    dtset%iprcel=21
     afford=0
   end if

!  First compute dimensions
   if(dtset%iprcel>=21 .or. dtset%iscf==-1)then
!    With dielop=1, the matrices will be computed when istep=dielstrt
!    With dielop=2, the matrices will be computed when istep=dielstrt and 1
     dielop=1
     if(dtset%iprcel>=41)dielop=2
     if((dtset%iprcel >= 71).and.(dtset%iprcel<=79)) dielop=0 !RSkerker preconditioner do not need the susceptibility matrix
!    Immediate computation of dielectric matrix
     dielstrt=1
!    Or delayed computation
     if(modulo(dtset%iprcel,100)>21 .and. modulo(dtset%iprcel,100)<=29)dielstrt=modulo(dtset%iprcel,100)-20
     if(modulo(dtset%iprcel,100)>31 .and. modulo(dtset%iprcel,100)<=39)dielstrt=modulo(dtset%iprcel,100)-30
     if(modulo(dtset%iprcel,100)>41 .and. modulo(dtset%iprcel,100)<=49)dielstrt=modulo(dtset%iprcel,100)-40
     if(modulo(dtset%iprcel,100)>51 .and. modulo(dtset%iprcel,100)<=59)dielstrt=modulo(dtset%iprcel,100)-50
     if(modulo(dtset%iprcel,100)>61 .and. modulo(dtset%iprcel,100)<=69)dielstrt=modulo(dtset%iprcel,100)-60
!    Get diecut, and the fft grid to be used for the susceptibility computation
     diecut=abs(dtset%diecut)
     if( dtset%diecut<0.0_dp )then
       ecutsus=ecut
     else
       ecutsus= ( sqrt(ecut) *0.5_dp + sqrt(diecut) *0.25_dp )**2
     end if
!    Impose sequential calculation
     ngfftdiel(1:3)=0 ; ngfftdiel(7)=100 ; ngfftdiel(9)=0; ngfftdiel(8)=dtset%ngfft(8);ngfftdiel(10:18)=0
     if(dtset%iscf==-1)ngfftdiel(7)=102
     write(6,*) 'call getng diel'
     call getng(dtset%boxcutmin,ecutsus,gmet,mpi_enreg_diel%me_fft,mgfftdiel,nfftdiel,ngfftdiel,&
&     mpi_enreg_diel%nproc_fft,dtset%nsym,mpi_enreg_diel%fft_option_lob,mpi_enreg_diel%paral_fft,dtset%symrel)
!    Compute the size of the dielectric matrix
     kpt_diel(1:3)=(/ 0.0_dp, 0.0_dp, 0.0_dp /)
     call getmpw(diecut,dtset%exchn2n3d,gmet,(/1/),kpt_diel,&
&     mpi_enreg_diel,npwdiel,1)
     lmax_diel=0
     if (psps%usepaw==1) then
       do ii=1,dtset%ntypat
         lmax_diel=max(lmax_diel,pawtab(ii)%lcut_size)
       end do
     end if
   else
     npwdiel=1
     mgfftdiel=1
     nfftdiel=1
     lmax_diel=0
   end if

!  Now, performs allocation
   allocate(dielinv(2,npwdiel*afford,dtset%nspden,npwdiel,dtset%nspden))
   allocate(susmat(2,npwdiel*afford,dtset%nspden,npwdiel,dtset%nspden))
   allocate(kg_diel(3,npwdiel))
   allocate(gbound_diel(2*mgfftdiel+8,2))
   allocate(irrzondiel(nfftdiel**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
   allocate(phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
   allocate(ph1ddiel(2,3*(2*mgfftdiel+1)*dtset%natom*psps%usepaw))
   allocate(ylmdiel(npwdiel,lmax_diel**2))
!  Then, compute the values of different arrays
   if(dielop>=1)then
     call status(0,dtfil%filstat,iexit,level,'kpgio(sus)    ')
!    Note : kgnam is dummy, npwarr_diel is dummy, npwtot_diel is dummy
!    This kpgio call for going from the suscep FFT grid to the diel sphere
     npwarr_diel(1)=npwdiel

     call kpgio(diecut,dtset%exchn2n3d,gmet,(/1/),kg_diel,kgnam,&
&     kpt_diel,1,(/1/),1,'COLL',mpi_enreg_diel,npwdiel,&
&     npwarr_diel,npwtot_diel,dtset%nsppol,tmp_unit)
     call sphereboundary(gbound_diel,1,kg_diel,mgfftdiel,npwdiel)
     if (dtset%nsym>1 .and. dtset%iscf>0 ) then
!      Should replace this initialization of irrzondiel and phnonsdiel through setsym by a direct call to irrzg
       allocate(indsym_dum(4,dtset%nsym,dtset%natom),symrec_dum(3,3,dtset%nsym))
       call setsym(indsym_dum,irrzondiel,dtset%iscf,dtset%natom,&
&       nfftdiel,ngfftdiel,dtset%nspden,dtset%nsppol,dtset%nsym,phnonsdiel,&
&       dtset%symafm,symrec_dum,dtset%symrel,dtset%tnons,dtset%typat,xred)
       deallocate(indsym_dum,symrec_dum)
     end if
     if (psps%usepaw==1) then
       call getph(atindx,dtset%natom,ngfftdiel(1),ngfftdiel(2),ngfftdiel(3),ph1ddiel,xred)
       call initylmg(gprimd,kg_diel,kpt_diel,1,mpi_enreg_diel,lmax_diel,npwdiel,dtset%nband,1,&
&       npwarr_diel,dtset%nsppol,0,rprimd,tmp_unit,tmp_unit,ylmdiel,rhodum)
     end if
   end if

 else
   npwdiel=1
   mgfftdiel=1
   nfftdiel=1
 end if

 call status(0,dtfil%filstat,iexit,level,'further allocs')

!The big array containing functions defined on the fft grid is allocated now.
!Note however, that a zero value of mffmem will cause allocation with
!the third dimension set to zero ! In this case, another temporary
!will be used inside newvtr.
!This array is needed only in the self-consistent case
 if(nstep>0 .and. dtset%iscf>0) then
   if (psps%usepaw==0) npawmix=0
   if (psps%usepaw==1) allocate(f_paw(npawmix,n_fftgr*dtset%mffmem))
   if (psps%usepaw==1.and.dtset%pawmixdg==0) then
     ispmix=2;nfftmix=dtset%nfft;ngfftmix(:)=ngfft(:)
   else
     ispmix=1;nfftmix=nfftf;ngfftmix(:)=ngfftf(:)
   end if
   allocate(f_fftgr(ispmix*nfftmix,dtset%nspden,n_fftgr*dtset%mffmem))
 end if

 nkxc=0
!TDDFT - For a first coding
 if (dtset%nfreqsus>0 .and. dtset%ikhxc==0)nkxc=0 !MF no xc kernel
 if (dtset%nfreqsus>0 .and. dtset%ikhxc==1)nkxc=0 !MF no xc kern, but (later) RPA ok
 if (dtset%nfreqsus>0 .and. dtset%ikhxc==2)nkxc=1 !MF LDA xc kernel + (later) RPA
 if (dtset%iscf==-1 .and. dtset%nspden==1) nkxc=2
 if (dtset%iscf==-1 .and. dtset%nspden==2) nkxc=3
!Eventually need kxc-LDA when susceptibility matrix has to be computed
 if (dtset%iscf>0.and.modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79)) nkxc=min(2*dtset%nspden-1,3)
!Eventually need kxc-LDA for residual forces (when density mixing is selected)
 if (dtset%iscf>=10.and.dtset%usewvl==0.and.forces_needed>0 .and. &
& abs(dtset%iprcch)>=1.and.abs(dtset%iprcch)<=6.and.abs(dtset%iprcch)/=5) then
   if (dtset%xclevel==1.or.dtset%iprcch>=0) nkxc=min(2*dtset%nspden-1,3)
   if (dtset%xclevel==2.and.dtset%nspden==2.and.dtset%iprcch<0) nkxc=23
 end if
 allocate(kxc(nfftf,nkxc))

!This flag will be set to 1 just before an eventual change of atomic
!positions inside the iteration, and set to zero when the consequences
!of this change are taken into account.
 moved_atm_inside=0
!This flag will be set to 1 if the forces are computed inside the iteration.
 computed_forces=0

 if(dtset%wfoptalg==2)then
   allocate(shiftvector((dtset%mband+2)*dtset%nkpt))
   val_min=-1.0_dp
   val_max=zero
 else
   allocate(shiftvector(1))
 end if

 call status(0,dtfil%filstat,iexit,level,'berryphase    ')

!!PAW+DMFT: allocate structured datatype paw_dmft if dtset%usedmft=1
!call init_sc_dmft(dtset%dmftbandi,dtset%dmftbandf,dtset%mband,dtset%nkpt,&
!&  dtset%nsppol,dtset%usedmft,paw_dmft,dtset%usedmft)
!call print_sc_dmft(paw_dmft)

!Electric field initializations: initialize pel_cg(:) and p_ion(:)
 if (dtset%berryopt == 4) then
!  berrys phase with PAW needs cprj
   if (psps%usepaw==1) then
     usecprj=1
     allocate(cprj(dtset%natom,dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol))
     call cprj_alloc(cprj,0,dimcprj)
     iatom=0 ; iorder_cprj=1 ! cprj are not ordered
     call ctocprj(atindx,cg,1,cprj,gmet,gprimd,&
&     iatom,idir,iorder_cprj,dtset%istwfk,kg,dtset%kptns,&
&     dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,&
&     dtset%natom,nattyp,dtset%nband,dtset%natom,ngfft,dtset%nkpt,dtset%nloalg,&
&     npwarr,dtset%nspinor,dtset%nsppol,dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,&
&     ucvol,dtfil%unpaw,dtfil%unkg,dtfil%unylm,0,wffnow,xred,ylm,ylmgr)
   end if
   unit_out=0;if (dtset%prtvol >= 10) unit_out=ab_out
   optberry=1     ! compute polarization only
   pel_cg(:) = zero;pelev_dum=zero
   call berryphase_new(atindx1,cg,cprj,dtefield,dtfil,dtset,gprimd,hdr,kg,&
&   dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,npwarr,dtset%nspinor,&
&   dtset%nsppol,psps%ntypat,dtset%nkpt,optberry,pawang,pawrad,pawrhoij,pawtab,&
&   pel_cg,pelev_dum,pion,psps,pwind,&
&   pwind_alloc,pwnsfac,rprimd,dtset%typat,ucvol,&
&   unit_out,usecprj,psps%usepaw,wffnow,xred,psps%ziontypat)
!  deallocate cprj
   if (psps%usepaw==1) then
     usecprj=0
     call cprj_free(cprj)
     deallocate(cprj)
   end if
 end if

 if (dtset%iscf==22) energies%h0=zero
 call timab(54,2,tsec)

! =========================== chen hack ================
!! if (.false.) then
 if (psps%usepaw==1) then 
!   Initialize 
!    (1) extpot_lm
!    (2) extpot_dij
! 
!   arrays:
!    extpot_lm(hc_meshsz,hc_lsize**2,natom)     : (l,m) moments of extpot on a radial mesh
!    extpot_dij(hc_meshsz,hc_lmnsize,natom)    : dij components due to extpot 
!    hc_rho1(hc_meshsz,hc_lsize**2,dtset%natom) : rho1 for each atom
!    hc_trho1(hc_meshsz,hc_lsize**2,dtset%natom): \tilde rho1 for each atom
!
!   integers:
!    hc_meshsz    : max of mesh_size of all atoms
!    hc_lmax      : Maximum value of angular momentum l+1
!    hc_lsize     : Maximum value of angular momentun + 1 for non zero gaunt coeff
!    hc_lmn2size  : max of lmn2_size of all atoms
!
   hc_rmax=0.d0
   hc_meshsz=0         ! for mesh where extpot is projected
   hc_lmn2size=0       ! for extpot_dij 
   do iatom=1,dtset%natom
     c_itype=dtset%typat(iatom)
     if (hc_rmax<pawrad(c_itype)%rmax) hc_rmax=pawrad(c_itype)%rmax
     if (hc_meshsz<pawrad(c_itype)%mesh_size) hc_meshsz=pawrad(c_itype)%mesh_size
     if (hc_lmn2size<pawrhoij(c_itype)%lmn2_size) hc_lmn2size=pawrhoij(c_itype)%lmn2_size
   enddo
   hc_lmax=pawang%l_max  ! Maximum value of angular momentum l+1
   hc_lsize=pawang%l_size_max

   allocate(extpot_lm (hc_meshsz,hc_lsize**2,dtset%natom))       ! (l,m) moments of extpot
   allocate(extpot_dij(hc_lmn2size,dtset%natom))                 ! dij due to extpot
   allocate(hc_rho1   (hc_meshsz,hc_lsize**2,dtset%natom))       ! rho1
   allocate(hc_trho1  (hc_meshsz,hc_lsize**2,dtset%natom))       ! \tilde rho1

   write(message,'(a)')'(chen/scfcv) ------------------ '; call wrtout(std_out,message,"COLL")
   write(message,'(a,I6)')   '   dtset%iscf  =',dtset%iscf; call wrtout(std_out,message,"COLL")
   write(message,'(a,F8.4)') '   hc_rmax     =',hc_rmax;   call wrtout(std_out,message,"COLL")
   write(message,'(a,I6)')   '   hc_meshsz   =',hc_meshsz;   call wrtout(std_out,message,"COLL")
   write(message,'(a,I6)')   '   hc_lmax     =',hc_lmax;   call wrtout(std_out,message,"COLL")
   write(message,'(a,I6)')   '   hc_lsize    =',hc_lsize;   call wrtout(std_out,message,"COLL")
   write(message,'(a,I6)')   '   hc_lmn2size =',hc_lmn2size;   call wrtout(std_out,message,"COLL")

!  calculate (l,m) moments of the incoming extpot: extpot_lm
!
   call xredxcart(dtset%natom,1,rprimd,xcart,xred)
   k0=zero
   call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,6,k0,ngfftf)
   call c_uni2rad(mpi_enreg,dtset%paral_kgb,nfftf,ngfftf, & 
     dtset%natom,dtset%typat,psps,xcart,gprimd,pawrad,extpot_lm,ucvol,gsqcut)

!DEBUG
!  if (mpi_enreg%me==4) then 
!    open(unit=111,file='extpot_lm',action='write')
!    write(111,*) extpot_lm
!    close(111)
!  endif
!  call leave_new('COLL')
!ENDDEBUG  

!  calculate D_ij due to extpot, extpot_dij does not change during SCF (below)
!
   call c_extpot_paw_dij(dtset%natom,psps%ntypat,mpi_enreg,pawang,pawtab,paw_ij,paw_an,pawrad,dtset%typat, &
          extpot_lm,extpot_dij)

 endif ! if psps%usepaw = 1
!
!==================== End of hack ===================

!######################################################################
!PERFORM ELECTRONIC ITERATIONS
!######################################################################

!Offer option of computing total energy with existing
!wavefunctions when nstep<=0, else do nstep iterations
!Note that for non-self-consistent calculations, this loop will be exited
!after the first call to vtorho
!Pass through the first routines even when nstep==0

 do istep=1,max(1,nstep)

   if (moved_atm_inside==1 .or. istep==1) then
!    ######################################################################
!    The following steps are done once for a given set of atomic
!    coordinates or for the nstep=1 case
!    ----------------------------------------------------------------------

!    Eventually symmetrize atomic coordinates over space group elements:
     call status(istep,dtfil%filstat,iexit,level,'call symzat   ')
     call symzat(indsym,dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,xred)

     if (dtset%usewvl == 0) then
!      Get cut-off for g-vectors
       if (psps%usepaw==1) then
         write(message,'(2a)') ch10,' FFT (fine) grid used in SCF cycle:'
         call wrtout(std_out,message,'COLL')
       end if
       k0=zero
       call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,6,k0,ngfftf)

!      Compute structure factor phases and large sphere cut-off (gsqcut):
       call status(istep,dtfil%filstat,iexit,level,'call getph    ')
       call getph(atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred)
       if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
         call getph(atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)
       else
         ph1df(:,:)=ph1d(:,:)
       end if
     end if

!    Initialization of atomic data for PAW
     if (psps%usepaw==1) then
!      Check for non-overlapping spheres
       call status(istep,dtfil%filstat,iexit,level,'call chkpawovlp')
       call chkpawovlp(dtset%natom,psps%ntypat,dtset%pawovlp,pawtab,rmet,dtset%typat,xred)
!      Identify parts of the rectangular grid where the density has to be calculated
       optgr0=dtset%pawstgylm;optgr1=0;optgr2=0;optrad=1-dtset%pawstgylm
       if (forces_needed==1.or.(dtset%xclevel==2.and.dtset%pawnhatxc>0.and.usexcnhat>0)) then
         optgr1=dtset%pawstgylm;if (stress_needed==1) optrad=1
       end if
       call status(istep,dtfil%filstat,iexit,level,'call nhatgrid ')
       call nhatgrid(atindx1,gmet,mpi_enreg,dtset%natom,dtset%natom,nattyp,ngfftf,psps%ntypat,&
&       optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,rprimd,ucvol,xred)
     end if

!    If we are inside SCF cycle or inside dynamics over ions,
!    we have to translate the density of previous iteration
     moved_rhor=0

!============================ chen ==================================
! In fact, one does not need initialize rhor for an inverting-KS case
!    
!     if (initialized/=0.and.dtset%usewvl == 0.and.ipositron/=1.and. &
!&     (abs(dtset%iprcch)==2.or.abs(dtset%iprcch)==5.or.abs(dtset%iprcch)==6)) then
!       moved_rhor=1
!       if (abs(dtset%iprcch)==2) then
!         option=2;allocate(workr(nfftf,dtset%nspden))
!         call status(istep,dtfil%filstat,iexit,level,'call fresid   ')
!         call fresid(dtset,gresid,mpi_enreg,nfftf,ngfftf,&
!&         psps%ntypat,option,pawtab,rhor,rprimd,&
!&         ucvol,workr,xred,xred_old,psps%znuclpsp)
!         rhor=workr;deallocate(workr)
!       else if (abs(dtset%iprcch)==5.or.abs(dtset%iprcch)==6) then
!         call status(istep,dtfil%filstat,iexit,level,'call extraprho')
!         scf_history%icall=scf_history%icall+1
!         call extraprho(atindx,atindx1,cg,dtset,gmet,gprimd,gsqcut,scf_history%icall,kg,mgfftf,mpi_enreg,&
!&         psps%mqgrid_vl,nattyp,nfftf,ngfftf,npwarr,psps%ntypat,pawrhoij,pawtab,&
!&         ph1df,psps,psps%qgrid_vl,rhor,rprimd,scf_history,ucvol,psps%usepaw,&
!&         xred,xred_old,ylm,psps%ziontypat,psps%znuclpsp)
!       end if
!       call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
!     end if
!================================ end of hack ============

   end if ! moved_atm_inside==1 .or. istep==1

!  Initialize/update data in the electron-positron case
   if (dtset%positron<0.or.(dtset%positron>0.and.istep==1)) then
     call setup_positron(atindx,atindx1,cg,dtfil,dtset,ecore,eigen,etotal,electronpositron,&
&     energies,gmet,forces_needed,fred,grewtn,gsqcut,hdr,initialized0,indsym,istep,istep_mix,kg,&
&     kxc,mgfftf,mpi_enreg,n3xccc,nattyp,nfftf,ngfftf,nhat,nkxc,npwarr,nvresid,occ,optres,&
&     paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,pel,ph1df,ph1d,pion,psps,rhog,rhor,&
&     rprimd,stress_needed,strsxc,symrec,ucvol,usexcnhat,vhartr,vpsp,vxc,wffnow,&
&     xccc3d,xred,ylm,ylmgr)
     ipositron=electronpositron_calctype(electronpositron)
   end if
   if ((moved_atm_inside==1 .or. istep==1).or. &
&   (dtset%positron<0.and.istep_mix==1)) then

!    PAW only: we sometimes have to compute compensation density
!    and eventually add it to density from WFs
     nhatgrdim=0
     if (psps%usepaw==1.and.(dtset%positron>=0.or.ipositron/=1) &
&     .and.((usexcnhat==0) &
&     .or.(dtset%xclevel==2.and.(dtfil%ireadwf/=0.or.dtfil%ireadden/=0.or.initialized/=0)) &
&     .or.(dtfil%ireadwf/=0.and.dtfil%ireadden==0.and.initialized==0))) then
       call timab(558,1,tsec)
       nhatgrdim=0;if (dtset%xclevel==2) nhatgrdim=usexcnhat*dtset%pawnhatxc
       ider=2*nhatgrdim
       if (nhatgrdim>0) allocate(nhatgr(nfftf,dtset%nspden,3))
       izero=0;k0(:)=zero
       call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,mpi_enreg,dtset%natom,dtset%natom,&
&       nfftf,ngfftf,nhatgrdim,dtset%nspden,psps%ntypat,dtset%paral_kgb,pawang,pawfgrtab,&
&       nhatgr,nhat,pawrhoij,pawrhoij,pawtab,k0,rprimd,ucvol)
       if (dtfil%ireadwf/=0.and.dtfil%ireadden==0.and.initialized==0) then
         rhor(:,:)=rhor(:,:)+nhat(:,:)
         call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
       end if
       call timab(558,2,tsec)
     end if


!    The following steps have been gathered in the setvtr routine:
!    - get Ewald energy and Ewald forces
!    - compute local ionic pseudopotential vpsp
!    - eventually compute 3D core electron density xccc3d
!    - eventually compute vxc and vhartr
!    - set up vtrial

!    DEBUG
!     write(6,*)' scfcv : before setvtr, energies%e_hartree=',energies%e_hartree
!    ENDDEBUG

     call status(istep,dtfil%filstat,iexit,level,'call setvtr   ')
     if (dtset%usewvl == 0) then
       optene = 4 * optres
       if(dtset%iscf==-3)optene=4
     else
!      We need the Hartree energy for the wavefunctions mixing
       optene = 1
     end if
     call setvtr(atindx1,dtset,energies,gmet,gprimd,grewtn,gsqcut,&
&     istep,kxc,mgfftf,moved_atm_inside,moved_rhor,mpi_enreg,&
&     nattyp,nfftf,ngfftf,nhat,nhatgr,nhatgrdim,nkxc,psps%ntypat,n1xccc,n3xccc,&
&     optene,pawtab,ph1df,psps,rhog,rhor,rmet,rprimd,strsxc,&
&     ucvol,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,&
&     xccc3d,xred,electronpositron=electronpositron,taug=taug,taur=taur)
!    DEBUG
     write(6,*)' scfcv : after setvtr, energies%e_hartree=',energies%e_hartree
!    ENDDEBUG

!======================= chen: setvtr sets e_localpsp ============================
!    NOTICE: 
!    For total energy, we only need to add local energy related to extpot
!    manully, thanks to the use of double-counting scheme in 
!    etotfor() subroutine, we do not need to worry about on-site
!    energies due to extpot: \int(extpot(n1-tn1), since extpot is a linear operator, 
!    therefore no double-counting correction is needed.
!
     if (psps%usepaw==1) then
       energies%e_localpsp=energies%e_localpsp+ & 
          ucvol/REAL(nfftf,kind=dp)*sum((rhor(:,1)-nhat(:,1))*extpot)
     endif

     IF (whoAmI==1 .or. whoAmI==2) THEN 
       write(message,'(a)') & 
        '             ITER      ETOT             DELTA(E)   residm    res2'
       call c_wrtlog(message)
       IF (psps%usepaw ==1 ) then 
       ELSE
!         Adding contribution from embedding potential
!           
!        for norm-conserving case: 
!        (1) Add extpot to vpsp (ionic local potential)
!        (2) Add extpot to vtrial 
         call c_wrtlog('NCPP: extpot added')
         vpsp=vpsp+extpot
         vtrial(:,1)=vtrial(:,1)+extpot(:)
         if (dtset%nsppol==2) then 
           vtrial(:,2)=vtrial(:,2)+extpot
         endif
       ENDIF
     ENDIF
!==================== end of hack ===============================

     if (nhatgrdim>0.and.nstep>0) deallocate(nhatgr)

!    Recursion Initialisation
     if(dtset%userec==1 .and. istep==1)  then
       rec_set%quitrec = 0
!      --At any step calculate the metric
       call Init_MetricRec(rec_set%inf,rec_set%nl%nlpsp,rmet,ucvol,rprimd,xred,&
&       dtset%ngfft(1:3),dtset%natom,rec_set%debug)
       if(initialized==0)  call first_rec(dtset,psps,rec_set)
     end if


!    End the condition of atomic position change or istep==1
   end if

!  ######################################################################
!  The following steps are done at every iteration
!  ----------------------------------------------------------------------

!  PAW: Compute energies and potentials in the augmentation regions (spheres)
!  Compute pseudopotential strengths (Dij quantities)
   if (psps%usepaw==1)then
!    "on-site" energies, potentials, densities computation
     nzlmopt=0;if (istep_mix==2.and.dtset%pawnzlm>0) nzlmopt=-1
     if (istep_mix>2) nzlmopt=dtset%pawnzlm
     option=0;if (dtset%iscf>0.and.dtset%iscf<10.and.nstep>0) option=1
     do iatom=1,dtset%natom
       itypat=dtset%typat(iatom)
       v_size=paw_an(iatom)%lm_size;if (dtset%pawxcdev==0) v_size=paw_an(iatom)%angl_size
       paw_ij(iatom)%has_dijhartree=1
       paw_an(iatom)%has_vxc=1
       allocate(paw_ij(iatom)%dijhartree(pawtab(itypat)%lmn2_size))
       allocate(paw_an(iatom)%vxc1 (pawtab(itypat)%mesh_size,v_size,paw_an(iatom)%nspden))
       allocate(paw_an(iatom)%vxct1(pawtab(itypat)%mesh_size,v_size,paw_an(iatom)%nspden))
       if (pawtab(itypat)%useexexch>0) allocate(paw_an(iatom)%vxc_ex(pawtab(itypat)%mesh_size,v_size,paw_an(iatom)%nspden))
     end do
!    Local exact exch.: impose occ. matrix if required
     if (dtset%useexexch>0) &
&     call setrhoijpbe0(dtset,psps%indlmn,initialized0,istep,istep_mix,psps%lmnmax,&
&     dtset%natom,dtset%ntypat,pawrhoij,pawtab,dtset%typat)
     call status(istep,dtfil%filstat,iexit,level,'call pawdenpot')
!    Computation of on-site densities/potentials/energies
     call pawdenpot(compch_sph,energies%e_paw,energies%e_pawdc,ipert,dtset%ixc,mpi_enreg,dtset%natom,dtset%natom,&
&     dtset%nspden,psps%ntypat,nzlmopt,option,dtset%paral_kgb,paw_an,paw_an,paw_ij,pawang,dtset%pawprtvol,pawrad,&
&     pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%spnorbscl,dtset%xclevel,psps%znuclpsp,&
&     electronpositron=electronpositron)
!    PAW+U: impose density matrix if required
     if (dtset%usepawu>0.and.(ipositron/=1)) then
       impose_dmat=0
       if ((istep<=abs(dtset%usedmatpu)).and.(dtset%usedmatpu<0.or.initialized0==0)) impose_dmat=1
       if (impose_dmat==1.or.dtset%dmatudiag/=0) then
         dimdmat=0;if (impose_dmat==1) dimdmat=2*lpawumax+1
         call setnoccmmp(0,dimdmat,dtset%dmatpawu(1:dimdmat,1:dimdmat,1:dtset%nsppol*dtset%nspinor,1:dtset%natpawu*impose_dmat),&
&         dtset%dmatudiag,impose_dmat,indsym,dtset%natom,dtset%natpawu,&
&         dtset%nspinor,dtset%nsppol,dtset%nsym,dtset%ntypat,paw_ij,pawang,dtset%pawprtvol,&
&         pawrhoij,pawtab,dtset%spinat,dtset%symafm,0,dtset%usepawu)
!        Reinitalize mixing if PAW+U and occupation matrix now allowed to change
!        For experimental purpose...
         if ((dtset%userib==1234).and.(istep==abs(dtset%usedmatpu)+1).and.(dtset%usedmatpu<0.or.initialized0==0)) istep_mix=1
       end if
     end if
!    Dij computation
     call status(istep,dtfil%filstat,iexit,level,'call pawdij   ')
     call pawdij(cplex,dtset,dtset%enunit,fatvshift,ipert,mpi_enreg,dtset%natom,dtset%natom,nfftf,ngfftf,dtset%nspden,psps%ntypat,&
&     dtset%paral_kgb,paw_an,paw_ij,pawang,pawfgrtab,dtset%pawprtvol,pawrad,dtset%pawspnorb,pawtab,dtset%pawxcdev,&
&     dtset%typat,ucvol,vtrial,vxc,electronpositron=electronpositron)

     do iatom=1,dtset%natom
       itypat=dtset%typat(iatom)
       deallocate(paw_ij(iatom)%dijhartree,paw_an(iatom)%vxc1,paw_an(iatom)%vxct1)
       paw_an(iatom)%has_vxc=0
       paw_ij(iatom)%has_dijhartree=0
       if (pawtab(itypat)%useexexch>0) deallocate(paw_an(iatom)%vxc_ex)
     end do
     call status(istep,dtfil%filstat,iexit,level,'call symdij   ')
     call symdij(gprimd,psps%indlmn,indsym,ipert,psps%lmnmax,dtset%natom,dtset%nsym,psps%ntypat,0,paw_ij,pawang,&
&     dtset%pawprtvol,rprimd,dtset%symafm,symrec,dtset%typat)
     if (paw_ij(1)%has_dijhat>0) &
&     call symdij(gprimd,psps%indlmn,indsym,ipert,psps%lmnmax,dtset%natom,dtset%nsym,psps%ntypat,1,paw_ij,pawang,&
&     dtset%pawprtvol,rprimd,dtset%symafm,symrec,dtset%typat)

!=============== chen ==================================     
     if (mpi_enreg%nproc_atom/=1) then
       print *,'src/95_drive/scfcv.F90: nproc_atom/=1, stop'
       call leave_new("COLL")
     endif
     do iatom=1,dtset%natom
       lmn=paw_ij(iatom)%lmn2_size
       paw_ij(iatom)%dij(1:lmn,1)=paw_ij(iatom)%dij(1:lmn,1)+extpot_dij(1:lmn,iatom) 
     enddo
!========================== end of hack ================

   end if

!  Write out occupancies to dtpawuj-dataset
   if (dtset%usepawu>0.and.dtset%macro_uj>0.and.istep>1.and.ipositron/=1) then
     call pawuj_red(dtset,dtpawuj,fatvshift,dtset%natom,dtset%ntypat,paw_ij,pawtab,ndtpawuj)
   end if

!  No need to continue and call vtorho, when nstep==0
   if(nstep==0)exit

!  ######################################################################
!  The following steps are done only when nstep>0
!  ----------------------------------------------------------------------
   call timab(56,1,tsec)
   call status(istep,dtfil%filstat,iexit,level,'loop istep    ')

   if(dtset%iscf>0)then
     write(message, '(a,a,i4)' )ch10,' ITER STEP NUMBER  ',istep
     call wrtout(std_out,message,'COLL')
   end if

!  The next flag says whether the xred have to be changed in the current iteration
   moved_atm_inside=0
   if(dtset%ionmov==4 .and. mod(iapp,2)/=1 .and. dtset%iscf>0 )moved_atm_inside=1
   if(dtset%ionmov==5 .and. iapp/=1 .and. istep==1 .and. dtset%iscf>0)moved_atm_inside=1

!  The next flag says whether the forces have to be computed in the current iteration
   computed_forces=0
   if ((dtset%optforces==1 .and. dtset%usewvl == 0).or.(moved_atm_inside==1)) computed_forces=1
   if (abs(tollist(3))>tiny(0._dp)) computed_forces=1
   if (dtset%iscf<0) computed_forces=0
   if ((istep==1).and.(dtset%optforces/=1)) then
     if (moved_atm_inside==1) then
       write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&       ' scfcv : WARNING -',ch10,&
&       '  Although the computation of forces during electronic iterations',ch10,&
&       '  was not required by user, it is done (required by the',ch10,&
&       '  choice of ionmov input parameter).'
       call wrtout(std_out,message,'COLL')
     end if
     if (abs(tollist(3))+abs(tollist(7))>tiny(0._dp)) then
       write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&       ' scfcv : WARNING -',ch10,&
&       '  Although the computation of forces during electronic iterations',ch10,&
&       '  was not required by user, it is done (required by the',ch10,&
&       '  "toldff" or "tolrff" tolerance criteria).'
       call wrtout(std_out,message,'COLL')
     end if
   end if
   if ((istep==1).and.(dtset%optforces==1).and. dtset%usewvl == 1) then
     write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&     ' scfcv : WARNING -',ch10,&
&     '  Although the computation of forces during electronic iterations',ch10,&
&     '  was required by user, it has been disable since the tolerence',ch10,&
&     '  is not on forces (force computation is expensive in wavelets).'
     call wrtout(std_out,message,'COLL')
   end if

   call timab(56,2,tsec)

!  ######################################################################
!  Compute the density rho from the trial potential
!  ----------------------------------------------------------------------

!  Compute the density from the trial potential
   if (dtset%tfkinfunc==0) then

!======================= chen: add extpot to vtrial ===============================
     if (psps%usepaw==1) then 
       vtrial(:,1)=vtrial(:,1)+extpot(:)
     endif
!======================= end of hack ==============================================

     call status(istep,dtfil%filstat,iexit,level,'call vtorho   ')
     call vtorho(afford,atindx,atindx1,cg,compch_fft,cpus,dbl_nnsclo,&
&     dielop,dielstrt,dphase,dtefield,dtfil,dtset,&
&     eigen,electronpositron,energies,etotal,gbound_diel,&
&     gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
&     istep,istep_mix,kg,kg_diel,kxc,lmax_diel,mgfftdiel,mpi_enreg,&
&     psps%mpsang,dtset%natom,nattyp,nfftf,nfftdiel,ngfftdiel,nhat,nkxc,&
&     npwarr,npwdiel,res2,dtset%nspinor,psps%ntypat,nvresid,occ,computed_forces,&
&     optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
&     pwind,pwind_alloc,pwnsfac,resid,residm,&
&     rhog,rhor,rmet,rprimd,susmat,symrec,taug,taur,&
&     ucvol,wffnew,wffnow,vtrial,wvl,xred,ylm,ylmdiel)
     call status(istep,dtfil%filstat,iexit,level,'after vtorho  ')

!======================= chen: remove extpot to vtrial ============================
     if (psps%usepaw==1) then
       vtrial(:,1)=vtrial(:,1)-extpot(:)
     endif
!======================= end of hack ==============================================

   elseif (dtset%tfkinfunc==1) then
     write(6,*)'WARNING : THOMAS FERMI'
     call vtorhotf(dtfil,dtset,energies%e_kinetic,energies%e_nonlocalpsp,energies%entropy,energies%e_fermie,&
&     gprimd,grnl,irrzon,mpi_enreg,dtset%natom,nfftf,dtset%nspden,dtset%nsppol,dtset%nsym,phnons,&
&     rhog,rhor,rprimd,ucvol,vtrial)
     residm=zero
     energies%e_eigenvalues=zero
   end if

!  Recursion method
   if(dtset%userec==1)then
     call vtorhorec(dtset,&
&     energies%e_kinetic,energies%e_nonlocalpsp,energies%entropy,energies%e_eigenvalues,energies%e_fermie,&
&     grnl,initialized,irrzon,mpi_enreg,nfftf,phnons,&
&     rhog,rhor,vtrial,rec_set,istep-nstep,rprimd,gprimd)
     residm=zero
   end if

   if(dtset%wfoptalg==2)then
     do ikpt=1,dtset%nkpt
       shiftvector(1+(ikpt-1)*(dtset%mband+2))=val_min
       shiftvector(2+(ikpt-1)*(dtset%mband+2):ikpt*(dtset%mband+2)-1)=&
&       eigen((ikpt-1)*dtset%mband+1:ikpt*dtset%mband)
       shiftvector(ikpt*(dtset%mband+2))=val_max
     end do
   end if

!  ######################################################################
!  Skip out of step loop if non-SCF (completed)
!  ----------------------------------------------------------------------

!  Indeed, nstep loops have been done inside vtorho
   if (dtset%iscf<=0) exit

!  ######################################################################
!  In case of density mixing or wavelet handling, compute the total energy
!  ----------------------------------------------------------------------
   if (dtset%iscf>=10 .or. dtset%usewvl == 1) then
     if (dtset%usewvl == 0) then
       optene = 1 ! use double counting scheme
     else if (dtset%iscf/=22) then
       optene = 0 ! use direct scheme for computation of energy
     else
       optene = -1
     end if

!    if the mixing is the ODA mixing, compute energy and new density here
     if (dtset%iscf==22) then
       call odamix(deltae,dtset,dtefield%efield_dot,elast,energies,&
&       etotal,gprimd,gsqcut,kxc,mpi_enreg,&
&       nfftf,ngfftf,nhat,nkxc,psps%ntypat,nvresid,n3xccc,optres,&
&       paw_ij,paw_an,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,pel_cg,pion,psps,rhog,rhor,rprimd,strsxc,&
&       taug,taur,ucvol,psps%usepaw,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,xccc3d)
     end if

!======= chen: double-counting scheme (optene==1) is needed for etotfor() ===============
     if (psps%usepaw==1 .and. optene/=1) then
       write(message,*)'scfcv.F90:  double-counting scheme is a MUST for etotfor(), ERROR, STOP'
       call wrtout(std_out,message,'COLL')
       call leave_new('COLL')
     endif
!=================== end of hack ========================================================

!    If the density mixing is required, compute the total energy here
     call etotfor(atindx1,deltae,diffor,dtset,dtefield%efield_dot,elast,&
&     electronpositron,energies,&
&     etotal,favg,fcart,forold,fred,gresid,grewtn,grhf,grnl,&
&     grxc,gsqcut,indsym,kxc,maxfor,mgfftf,mpi_enreg,nattyp,&
&     nfftf,ngfftf,nhat,nkxc,psps%ntypat,&
&     nvresid,n1xccc,n3xccc,optene,computed_forces,optres,&
&     pawang,pawfgrtab,pawrhoij,pawtab,pel_cg,ph1df,pion,psps,rhog,rhor,rprimd,symrec,synlgr,&
&     psps%usepaw,usexcnhat,vhartr,vpsp,vxc,xccc3d,xred)

   end if  ! iscf>=10 or usewvl==1

!  ######################################################################
!  In case of density mixing, check the exit criterion
!  ----------------------------------------------------------------------
   if (dtset%iscf>=10 .and. dtset%usewvl == 0) then
!    Check exit criteria
     call timab(52,1,tsec)
     call status(istep,dtfil%filstat,iexit,level,'call scprqt   ')
     choice=2
     call scprqt(choice,cpus,deltae,diffor,dtset,&
&     eigen,etotal,favg,fcart,energies%e_fermie,dtfil%fnameabo_app_eig,dtfil%filnam_ds(1),&
&     initialized0,dtset%iscf,istep,dtset%kptns,maxfor,moved_atm_inside,mpi_enreg,&
&     dtset%nband,dtset%nkpt,nstep,occ,optres,&
&     prtfor,quit,res2,resid,residm,response,tollist,psps%usepaw,&
&     vxcavg,dtset%wtk,xred,electronpositron=electronpositron)

!========== chen: Print interation information for PAW ===========
     write(message,'(a,i3,1p,g22.14,5es10.3)') & 
       '(scfcv)      ',istep,etotal,deltae,residm,res2
     call c_wrtlog(message)
!================= end of hack ===================================

!    Exit criteria for the recursion method
     if(dtset%userec==1.and.rec_set%quitrec==2)quit=1

     if (istep==nstep) quit=1
     call timab(52,2,tsec)

!    If criteria in scprqt say to quit, then exit the loop over istep.
     if(mpi_enreg%paral_compil_kpt==1)then
       quit_sum=quit
       call xcomm_world(mpi_enreg,spaceComm)
       call xsum_mpi(quit_sum,spaceComm,ierr)
       if (quit_sum > 0) quit=1
     end if ! mpi_enreg%paral_compil_kpt==1
     if (quit==1) exit
   end if

!  ######################################################################
!  Mix the total density (if required)
!  ----------------------------------------------------------------------
   if (dtset%iscf>=10 .and.dtset%iscf/=22.and. dtset%usewvl == 0) then
!    If LDA dielectric matrix is used for preconditionning, has to update here Kxc
     if (nkxc>0.and.modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79) &
&     .and.((istep==1.or.istep==dielstrt).or.(dtset%iprcel>=100))) then
       optxc=10

!      to be adjusted for the call to rhohxc
       nk3xc=1
       call rhohxc(dtset,edum,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,ngfftf,nhat,&
&       psps%usepaw,nhatgr,0,nkxc,nk3xc,dtset%nspden,n3xccc,optxc,rhog,&
&       rhor,rprimd,dummy2,0,vhartr,vxc,vxcavg_dum,xccc3d,taug=taug,taur=taur)
     end if

     call status(istep,dtfil%filstat,iexit,level,'call newrho   ')
     call newrho(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,dtn_pc,dtset,etotal,fcart,pawfgr%fintocoa,dtfil%fnametmp_fft,&
&     f_atm,f_fftgr,f_paw,gmet,grhf,gsqcut,initialized,&
&     ispmix,istep_mix,i_rhor,i_vresid,i_vrespc,i_vtrial,kg_diel,kxc,mgfftf,pawfgr%coatofin,&
&     moved_atm_inside,mpi_enreg,nattyp,nfftf,nfftmix,ngfftf,ngfftmix,nkxc,npawmix,npwdiel,&
&     nvresid,psps%ntypat,n_fftgr,n_index,n1xccc,pawrhoij,pawtab,&
&     ph1df,psps,rhog,rhor,rprimd,susmat,psps%usepaw,vtrial,xred)
   end if   ! iscf>=10

!  ######################################################################
!  Additional computation in case of an electric field
!  ----------------------------------------------------------------------

!  In case of an electric field calculation, need polarization
!  to compute electric enthalpy instead of energy.

!  if (psps%usepaw==0.and.dtset%berryopt == 4) then
   if (dtset%berryopt == 4) then
!    When using symmetry, it is costly to update polarization from changes in Zak
!    phases. It is better to call berryphase here.
!    Update polarization by adding increment from the SCF step
!    pel_cg(1:3) = pel_cg(1:3) + dtefield%sdeg*dphase(1:3)/two_pi
!    ptot(1:3) = pel_cg(1:3) + pion(1:3)
!    write(message,'(6(a),3(e16.9,2x),a,a,3(e16.9,2x),a,a,3(e16.9,2x))')ch10,&
!    &    ' scfcv: Polarization from accumulated change in Berry phase:',ch10,&
!    &    ' (reduced coordinates, a. u., without correcting for branch cuts)',ch10,&
!    &    '     Electronic: ', (pel_cg(ii), ii = 1, 3), ch10,&
!    &    '     Ionic:      ', (pion(ii), ii = 1, 3), ch10, &
!    &    '     Total:      ', (ptot(ii), ii = 1, 3)
!    call wrtout(std_out,message,'COLL')
!    berrys phase with PAW needs cprj
     if (psps%usepaw==1) then
       usecprj=1
       allocate(cprj(dtset%natom,dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol))
       call cprj_alloc(cprj,0,dimcprj)
       iatom=0 ; iorder_cprj=1 ! cprj are not ordered
       call ctocprj(atindx,cg,1,cprj,gmet,gprimd,&
&       iatom,idir,iorder_cprj,dtset%istwfk,kg,dtset%kptns,&
&       dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,&
&       dtset%natom,nattyp,dtset%nband,dtset%natom,ngfft,dtset%nkpt,dtset%nloalg,&
&       npwarr,dtset%nspinor,dtset%nsppol,dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,&
&       ucvol,dtfil%unpaw,dtfil%unkg,dtfil%unylm,0,wffnow,xred,ylm,ylmgr)
     end if
     pelev_dum=zero
     call berryphase_new(atindx1,cg,cprj,dtefield,dtfil,dtset,&
&     gprimd,hdr,kg,&
&     dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,npwarr,dtset%nspinor,&
&     dtset%nsppol,psps%ntypat,dtset%nkpt,optberry,pawang,pawrad,pawrhoij,&
&     pawtab,pel_cg,pelev_dum,pion,psps,pwind,&
&     pwind_alloc,pwnsfac,rprimd,dtset%typat,ucvol,&
&     unit_out,usecprj,psps%usepaw,wffnow,xred,psps%ziontypat)
     ptot(:) = pel_cg(:) + pion(:)
     write(message,'(6(a),3(e16.9,2x),a,a,3(e16.9,2x),a,a,3(e16.9,2x))')ch10,&
&     ' scfcv: New value of the polarization:',ch10,&
&     ' (reduced coordinates, a. u.)',ch10,&
&     '     Electronic: ', (pel_cg(ii), ii = 1, 3), ch10,&
&     '     Ionic:      ', (pion(ii), ii = 1, 3), ch10, &
&     '     Total:      ', (ptot(ii), ii = 1, 3)
     call wrtout(std_out,message,'COLL')
!    deallocate cprj
     if (psps%usepaw==1) then
       usecprj=0
       call cprj_free(cprj)
       deallocate(cprj)
     end if
   end if       ! berryopt

!  ######################################################################
!  Compute the new potential from the trial density
!  ----------------------------------------------------------------------

!  Set XC computation flag
   optxc=1
   if (nkxc>0) then
     if (dtset%nfreqsus>0) optxc=2
     if (dtset%iscf<0) optxc=2
     if (modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79).and. &
&     dtset%iscf<10.and. &
&     (dtset%iprcel>=100.or.istep==1.or.istep==dielstrt)) optxc=2
     if (dtset%iscf>=10.and.dtset%iprcch/=0.and.abs(dtset%iprcch)/=5) optxc=2
     if (optxc==2.and.dtset%xclevel==2.and.nkxc==3-2*mod(dtset%nspden,2)) optxc=12
   end if

   if (dtset%iscf/=22) then
!    PAW: eventually recompute compensation density (and gradients)
     nhatgrdim=0
     if (psps%usepaw==1) then
       ider=-1;if (dtset%iscf>=10.and.((dtset%xclevel==2.and.dtset%pawnhatxc>0).or.usexcnhat==0)) ider=0
       if (dtset%xclevel==2.and.dtset%pawnhatxc>0.and.usexcnhat>0) ider=ider+2
       if (ipositron==1) ider=-1
       if (ider>0) then
         nhatgrdim=1;allocate(nhatgr(nfftf,dtset%nspden,3))
       end if
       if (ider>=0) then
         call timab(558,1,tsec)
         izero=0;k0(:)=zero
         call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,mpi_enreg,dtset%natom,dtset%natom,nfftf,ngfftf,&
&         nhatgrdim,dtset%nspden,psps%ntypat,dtset%paral_kgb,pawang,pawfgrtab,nhatgr,nhat,&
&         pawrhoij,pawrhoij,pawtab,k0,rprimd,ucvol)
         call timab(558,2,tsec)
       end if
     end if
!    Compute new potential from the trial density
     call status(istep,dtfil%filstat,iexit,level,'call rhotov')
     optene=2*optres;if(psps%usepaw==1) optene=2

     call rhotov(dtset,energies,gprimd,gsqcut,kxc,mpi_enreg,nfftf,ngfftf, &
&     nhat,nhatgr,nhatgrdim,nkxc,nvresid,n3xccc,&
&     optene,optres,optxc,&
&     rhog,rhor,rprimd,strsxc,ucvol,psps%usepaw,usexcnhat,&
&     vhartr,vnew_mean,vpsp,vres_mean,res2,vtrial,vxcavg,vxc,xccc3d,&
&     electronpositron=electronpositron,taug=taug,taur=taur)

!======================= chen: add local energy due to extpot ===========
!
!================= 
!    rhotov() also computes
!      e_localpsp=local psp energy (hartree)
!      e_hartree=Hartree part of total energy (hartree units)
!      e_xc=exchange-correlation energy (hartree)
!      e_xcdc=exchange-correlation double-counting energy (hartree)
!==================

!    NOTICE: 
!    For total energy, we only need to add local energy related to extpot
!    manully, thanks to the use of double-counting scheme in 
!    etotfor() subroutine, we do not need to worry about on-site
!    energies due to extpot: \int(extpot(n1-tn1), since extpot is a linear operator, 
!    therefore no double-counting correction is needed.
     if (psps%usepaw==1) then
        energies%e_localpsp=energies%e_localpsp+ucvol/nfftf*sum((rhor(:,1)-nhat(:,1))*extpot)
     endif
!========================== end of hack =================================     

   end if

!  If the xred have to be changed in the current iteration, they has to be saved
   if(dtset%iextrapwf==1) scf_history%rprimd(:,:)=rprimd(:,:)

   if(moved_atm_inside==1) xred_old(:,:)=xred(:,:)


!  ######################################################################
!  Check exit criteria in case of potential mixing
!  ----------------------------------------------------------------------
   if (dtset%iscf<10 .and. dtset%usewvl == 0) then

!    If the potential mixing is required, compute the total energy here
!    PAW: has to compute here spherical terms
     if (psps%usepaw==1) then
       nzlmopt=0;if (istep_mix==1.and.dtset%pawnzlm>0) nzlmopt=-1
       if (istep_mix>1) nzlmopt=dtset%pawnzlm
       option=2
       do iatom=1,dtset%natom
         allocate(paw_ij(iatom)%dijhartree(pawtab(dtset%typat(iatom))%lmn2_size))
         paw_ij(iatom)%has_dijhartree=1
       end do
       call pawdenpot(compch_sph,energies%e_paw,energies%e_pawdc,ipert,dtset%ixc,mpi_enreg,dtset%natom,dtset%natom,&
&       dtset%nspden,psps%ntypat,nzlmopt,option,dtset%paral_kgb,paw_an,paw_an,paw_ij,pawang,dtset%pawprtvol,pawrad,&
&       pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%spnorbscl,dtset%xclevel,&
&       psps%znuclpsp,electronpositron=electronpositron)
       do iatom=1,dtset%natom
         deallocate(paw_ij(iatom)%dijhartree)
         paw_ij(iatom)%has_dijhartree=0
       end do
     end if

     call status(istep,dtfil%filstat,iexit,level,'call etotfor  ')
     call etotfor(atindx1,deltae,diffor,dtset,&
&     dtefield%efield_dot,elast,electronpositron,energies,&
&     etotal,favg,fcart,forold,fred,gresid,grewtn,grhf,grnl,&
&     grxc,gsqcut,indsym,kxc,maxfor,mgfftf,mpi_enreg,nattyp,&
&     nfftf,ngfftf,nhat,nkxc,dtset%ntypat,nvresid,n1xccc, &
&     n3xccc,0,computed_forces,optres,&
&     pawang,pawfgrtab,pawrhoij,pawtab,pel_cg,ph1df,pion,psps,rhog,rhor,rprimd,symrec,synlgr,&
&     psps%usepaw,usexcnhat,vhartr,vpsp,vxc,xccc3d,xred)
   end if

!  ######################################################################
!  Check exit criteria in case of potential mixing or wavelet handling
!  ----------------------------------------------------------------------
   if (dtset%iscf<10 .or. dtset%usewvl == 1) then
!    Check exit criteria
     call timab(52,1,tsec)
     call status(istep,dtfil%filstat,iexit,level,'call scprqt   ')
     choice=2
     call scprqt(choice,cpus,deltae,diffor,dtset,&
&     eigen,etotal,favg,fcart,energies%e_fermie,dtfil%fnameabo_app_eig,dtfil%filnam_ds(1),&
&     initialized0,dtset%iscf,istep,dtset%kptns,maxfor,moved_atm_inside,mpi_enreg,&
&     dtset%nband,dtset%nkpt,nstep,occ,optres,&
&     prtfor,quit,res2,resid,residm,response,tollist,psps%usepaw,&
&     vxcavg,dtset%wtk,xred,electronpositron=electronpositron)
     if (istep==nstep.and.psps%usepaw==1) quit=1
     call timab(52,2,tsec)

!========== chen: print information for NCPP case  =====================
   if (psps%usepaw==0) then 
     write(message,'(a,i3,1p,g22.14,5es10.3)') & 
       '(scfcv)      ',istep,etotal,deltae,residm,res2
     call c_wrtlog(message)
   endif
!================= end of hack =========================================

!    exit criteria for the recursion
     if(dtset%userec==1 .and. rec_set%quitrec==2) quit=1

!    If criteria in scprqt say to quit, then exit the loop over istep.
     if(mpi_enreg%paral_compil_kpt==1)then
       quit_sum=quit
       call xcomm_world(mpi_enreg,spaceComm)
       call xsum_mpi(quit_sum,spaceComm,ierr)
       if (quit_sum > 0) quit=1
     end if ! mpi_enreg%paral_compil_kpt==1
     if (quit==1) then
       do ispden=1,dtset%nspden
         vtrial(:,ispden)=vtrial(:,ispden)+nvresid(:,ispden)+vres_mean(ispden)
       end do
       exit
     end if
   end if

!  ######################################################################
!  Mix the potential (if required) - Check exit criteria
!  ----------------------------------------------------------------------
   if (dtset%iscf<10 .and. dtset%usewvl /= 1) then

!    Precondition the residual and forces, then determine the new vtrial
!    (Warning: the (H)xc potential may have been subtracted from vtrial)
     call status(istep,dtfil%filstat,iexit,level,'call newvtr   ')
     call newvtr(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,&
&     dtn_pc,dtset,energies%e_fermie,etotal,fcart,pawfgr%fintocoa,&
&     f_atm,f_fftgr,f_paw,gmet,grhf,gsqcut,initialized,ispmix,&
&     istep_mix,i_rhor,i_vresid,i_vrespc,i_vtrial,&
&     kg_diel,kxc,mgfftf,pawfgr%coatofin,&
&     moved_atm_inside,mpi_enreg,nattyp,nfftf,nfftmix,&
&     nhat,nhatgr,nhatgrdim,ngfftf,ngfftmix,nkxc,npawmix,npwdiel,&
&     nstep,psps%ntypat,n_fftgr,n_index,n1xccc,optres,optxc,&
&     pawrhoij, &
&     ph1df,psps,rhor,rprimd,susmat,psps%usepaw,&
&     vhartr,vnew_mean,vpsp,nvresid,&
&     vtrial,vxc,xred,atindx1,cg,deltae,&
&     dtfil,energies%e_eigenvalues,eigen,energies%e_kinetic,&
&     energies%e_nonlocalpsp,kg,nfftf,ngfftf,npwarr,n3xccc,occ,optene,&
&     pawfgr,pawtab,resid,rhog,usexcnhat,wffnow,ylm,dtset%nspinor,xccc3d)
   end if   ! iscf<10

!  No potential mixing in wavelet, direct minimisation scheme!
   if (dtset%usewvl == 1) then
     call status(istep,dtfil%filstat,iexit,level,'call wvl_newvtr')
     call wvl_newvtr(dtset, mpi_enreg, nele, offset, vhartr, vpsp, vtrial, vxc)
   end if


!  ######################################################################
!  END MINIMIZATION ITERATIONS
!  ######################################################################

!  The initialisation of the gstate run should be done when this point is reached
   initialized=1

!  This is to save the density for restart.
   if (mpi_enreg%paral_compil_kpt==0.or.mpi_enreg%me==0) then
     prtden=dtset%prtden
     if (prtden<0) then
       if (mod(istep-1,abs(prtden))==0) then
         isave=isave+1
         call status(0,dtfil%filstat,iexit,level,'call ioarr-den')
         rdwr=2 ; fformr=52 ; rdwrpaw=0
         call int2char4(mod(isave,2),tag)
         fildata=trim(dtfil%fnametmp_app_den)//'_'//tag
         accessfil = 0
         call ioarr(accessfil,rhor, dtset, etotal,fformr,fildata,hdr, mpi_enreg, &
&         nfftf,pawrhoij,rdwr,rdwrpaw)
       end if
     end if
     prtkden=dtset%prtkden
     if (prtkden<0) then
       if (mod(istep-1,abs(prtkden))==0) then
         isave=isave+1
         call status(0,dtfil%filstat,iexit,level,'call ioarr-kden')
         rdwr=2 ; fformr=52 ; rdwrpaw=0
         call int2char4(mod(isave,2),tag)
         fildata=trim(dtfil%fnametmp_app_kden)//'_'//tag
         accessfil = 0
         call ioarr(accessfil,taur, dtset, etotal,fformr,fildata,hdr, mpi_enreg, &
&         nfftf,pawrhoij,rdwr,rdwrpaw)
       end if
     end if
   end if

   if (nhatgrdim>0) deallocate(nhatgr)

   istep_mix=istep_mix+1
   if (ipositron/=0) electronpositron%istep_scf=electronpositron%istep_scf+1
 end do ! istep

 if (quit==1.and.nstep==1) initialized=1

!######################################################################
!Case nstep==0: compute energy based on incoming wf
!----------------------------------------------------------------------

 if(nstep==0) then
   optene=2*psps%usepaw+optres
   energies%entropy=results_gs%energies%entropy  !MT20070219: entropy is not recomputed in routine energy
   call energy(atindx,atindx1,cg,compch_fft,dtfil,dtset,electronpositron,energies,eigen,&
&   etotal,gsqcut,indsym,irrzon,kg,mpi_enreg,nattyp,nfftf,ngfftf,nhat,nhatgr,nhatgrdim,&
&   npwarr,dtset%nspinor,n3xccc,occ,optene,paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,&
&   phnons,ph1d,psps,resid,rhog,rhor,rprimd,strsxc,symrec,taug,taur,usexcnhat,vhartr,vtrial,&
&   vpsp,vxc,wffnow,wvl%wfs,xccc3d,xred,ylm)
   if (nhatgrdim>0) deallocate(nhatgr)
 end if ! nstep==0

!######################################################################
!Additional steps after SC iterations, including force, stress, polarization calculation
!----------------------------------------------------------------------

 if (dtset%userec==1) then
   call status(0,dtfil%filstat,iexit,level,'call prtene   ')
   call prtene(dtset,energies,ab_out,psps%usepaw)
   call prtene(dtset,energies,6,psps%usepaw)
 end if
 call timab(60,1,tsec)
 call status(0,dtfil%filstat,iexit,level,'endloop istep ')

!PAW: in some cases, need to recompute <p_lmn|Cnk> projected WF:
!should be output from vtorho (but is "type-sorted" in vtorho, not here)...
 if (psps%usepaw==1.and. &
& (dtset%prtwant==2.or.dtset%prtwant==3.or.dtset%prtnabla>0.or.dtset%prtdos==3 &
& .or.dtset%berryopt/=0.or.dtset%kssform==3.or.dtset%pawfatbnd>0.or.dtset%pawprtwf>0)) then
   usecprj=1
   allocate(cprj(dtset%natom,dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol))
   call cprj_alloc(cprj,0,dimcprj)
   iatom=0 ; iorder_cprj=1 ! cprj are not ordered
   call ctocprj(atindx,cg,1,cprj,gmet,gprimd,&
&   iatom,idir,iorder_cprj,dtset%istwfk,kg,dtset%kptns,&
&   dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,&
&   dtset%natom,nattyp,dtset%nband,dtset%natom,ngfft,dtset%nkpt,dtset%nloalg,&
&   npwarr,dtset%nspinor,dtset%nsppol,dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,&
&   ucvol,dtfil%unpaw,dtfil%unkg,dtfil%unylm,0,wffnow,xred,ylm,ylmgr)
 end if

!SHOULD CLEAN THE ARGS OF THIS ROUTINE
 call status(0,dtfil%filstat,iexit,level,'afterscfloop  ')
 call afterscfloop(atindx,atindx1,cg,computed_forces,cprj,cpus,&
& deltae,diffor,dtefield,dtfil,dtset,eigen,electronpositron,elfr,energies,etotal,&
& favg,fcart,forold,fred,gresid,grewtn,grhf,grhor,&
& grxc,gsqcut,hdr,indsym,irrzon,&
& istep,kg,kxc,lrhor,maxfor,mgfftf,&
& moved_atm_inside,mpi_enreg,&
& n3xccc,nattyp,&
& nfftf,ngfft,ngfftf,nhat,nkxc,npwarr,nvresid,&
& occ,optres,optxc,paw_an,paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,pel,pel_cg,&
& ph1d,ph1df,phnons,pion,prtfor,psps,pwind,pwind_alloc,pwnsfac,res2,resid,residm,results_gs,&
& rhog,rhor,rprimd,stress_needed,strsxc,strten,symrec,synlgr,taug,taur,tollist,&
& usecprj,usexcnhat,vhartr,vpsp,vxc,vxcavg,wffnow,wvl,xccc3d,xred,xred_old,ylm,ylmgr)

 call timab(60,2,tsec)

!######################################################################
!All calculations in scfcv are finished. Printing section
!----------------------------------------------------------------------


 call status(istep,dtfil%filstat,iexit,level,'call outscfcv ')
 
 call outscfcv(atindx1,cg,compch_fft,compch_sph,cprj,dimcprj,dtfil,dtset,&
& ecut,eigen,electronpositron,elfr,etotal,energies%e_fermie,gmet,gprimd,grhor,hdr,istep_mix,kg,&
& lrhor,dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,&
& nattyp,nfftf,ngfftf,nhat,dtset%nkpt,npwarr,dtset%nspden,dtset%nspinor,dtset%nsppol,&
& dtset%nsym,psps%ntypat,n3xccc,occ,&
& pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,paw_an,paw_ij,dtset%prtvol,psps,rhor,rprimd,&
& taur,ucvol,usecprj,usexcnhat,wffnow,vhartr,vtrial,vxc,xccc3d,xred)

!================ chen: obtain energy due to external pot (extpot) ==================
! We need to update rho1 and tilde_rho1 with the final rho_ij from vtorho()
! 
 if (psps%usepaw==1) then
   call c_paw_onsite_densities(mpi_enreg,&
     dtset%natom,dtset%natom,dtset%nspden,psps%ntypat,dtset%paral_kgb,paw_an,&
     paw_ij,pawang,pawrad,pawrhoij,pawtab,hc_rho1,hc_trho1)
!   
!  compute ext_pot*rho1 and ext_pot*trho1 (tilde-rho1)
!
   call c_int_extpot_rho1(dtset%natom,dtset%ntypat,paw_an,pawrhoij,pawrad,hc_rho1,extpot_lm,tmp1)
   call c_int_extpot_rho1(dtset%natom,dtset%ntypat,paw_an,pawrhoij,pawrad,hc_trho1,extpot_lm,tmp2)
   ext_energy = sum(extpot*(rhor(:,1)-nhat(:,1)))*ucvol/dble(nfftf) + tmp1-tmp2
   write(message,'(a,g22.14,a)') "(chen/scfcv): ext_energy     =",ext_energy," (evaluated on radial mesh)"; call wrtout(std_out,message,'COLL')
   write(message,'(a,g22.14,a)') "              on-site energy =",tmp1-tmp2, " (evaluated on radial mesh)";  call wrtout(std_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 endif

!========== end of hack =============================================================

 if(associated(elfr))then
   deallocate(elfr)
   nullify(elfr)
 end if

 if(associated(grhor))then
   deallocate(grhor)
   nullify(grhor)
 end if

 if(associated(lrhor))then
   deallocate(lrhor)
   nullify(lrhor)
 end if

 if(mpi_enreg%paral_compil_kpt==1)then
   call timab(61,1,tsec)
   call leave_test(mpi_enreg)
   call timab(61,2,tsec)
 end if

 call timab(60,1,tsec)

!Transfer eigenvalues computed by BigDFT in afterscfloop to eigen.
 if (dtset%usewvl == 1) then
   if (dtset%nsppol == 1) then
     eigen = wvl%wfs%eval
   else
     eigen(1:wvl%wfs%nstates_up) = wvl%wfs%eval(1:wvl%wfs%nstates_up)
     eigen(dtset%mband + 1:dtset%mband + wvl%wfs%nstates_dn) = &
&     wvl%wfs%eval(wvl%wfs%nstates_up + 1:wvl%wfs%nstates_up + wvl%wfs%nstates_dn)
   end if
 end if

!Debugging : print the different parts of rhor, as well as vxc
!MPIWF Warning : this should not be parallelized over space, leave this debugging feature as such.
 if(dtset%prtvol==-level)then
   write(message,'(a)') '   ir     vxc(ir)     rhor(ir)     '
   call wrtout(std_out,message,'COLL')
   do ir=1,nfftf
     if(ir<=11 .or. mod(ir,301)==0 )then
       write(message,'(i5,a,2es13.6)')ir,' ',vxc(ir,1),rhor(ir,1)
       call wrtout(std_out,message,'COLL')
       if(dtset%nspden==2)then
         write(message,'(a,2es13.6)')'      ',vxc(ir,2),rhor(ir,2)
         call wrtout(std_out,message,'COLL')
       end if
     end if
   end do
 end if

!Structured debugging : if prtvol=-level, stop here.
 if(dtset%prtvol==-level)then
   write(message,'(a1,a,a1,a,i1,a)') ch10,' scfcv : exit ',&
&   ch10,'  prtvol=-',level,', debugging mode => stop '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!######################################################################
!Deallocate memory and save results
!----------------------------------------------------------------------

 call status(0,dtfil%filstat,iexit,level,'deallocate    ')
 if (psps%usepaw==1.and. &
& (dtset%prtwant==2.or.dtset%prtnabla>0.or.dtset%prtdos==3.or.dtset%berryopt/=0.or.dtset%kssform==3)) then
   usecprj=0
   call cprj_free(cprj)
   deallocate(cprj)
 end if
 deallocate(fcart,fred,forold)
 deallocate(grnl,gresid,grewtn,grxc,synlgr)
 deallocate(ph1d,ph1df)
 deallocate(vhartr,vtrial,vpsp,vxc,xccc3d)
 deallocate(kxc,shiftvector)
 if (dtset%iscf>0) then
   deallocate(dtn_pc,f_atm,grhf,nvresid)
   if(nstep>0) deallocate(f_fftgr)
   if(nstep>0.and.psps%usepaw==1) deallocate(f_paw)
   deallocate(i_rhor,i_vtrial,i_vresid,i_vrespc)
 end if
 if((nstep>0.and.dtset%iscf>0).or.dtset%iscf==-1) then
   deallocate(dielinv,gbound_diel)
   deallocate(irrzondiel,kg_diel)
   deallocate(phnonsdiel,susmat)
   deallocate(ph1ddiel,ylmdiel)
 end if
 if (psps%usepaw==1) then
   deallocate(nhat)
   do iatom=1,dtset%natom
     itypat=dtset%typat(iatom)
!    Free pawfgrtab
     if (associated(pawfgrtab(iatom)%ifftsph))deallocate(pawfgrtab(iatom)%ifftsph)
     if (associated(pawfgrtab(iatom)%rfgd))    deallocate(pawfgrtab(iatom)%rfgd)
     if (associated(pawfgrtab(iatom)%gylm))    deallocate(pawfgrtab(iatom)%gylm)
     if (associated(pawfgrtab(iatom)%gylmgr))  deallocate(pawfgrtab(iatom)%gylmgr)
     if (associated(pawfgrtab(iatom)%gylmgr2))deallocate(pawfgrtab(iatom)%gylmgr2)
     if (associated(pawfgrtab(iatom)%nhatfr )) deallocate(pawfgrtab(iatom)%nhatfr)
     if (associated(pawfgrtab(iatom)%vlocgr )) deallocate(pawfgrtab(iatom)%vlocgr)
     pawfgrtab(iatom)%nfgd=0;pawfgrtab(iatom)%rfgd_allocated=0
     pawfgrtab(iatom)%gylm_allocated=0;pawfgrtab(iatom)%gylmgr_allocated=0
     pawfgrtab(iatom)%gylmgr2_allocated=0;pawfgrtab(iatom)%nhatfr_allocated=0
     pawfgrtab(iatom)%vlocgr_allocated=0
!    Free paw_an
     deallocate(paw_an(iatom)%lmselect,paw_ij(iatom)%dij)
     if (paw_an(iatom)%has_vxc>0) deallocate(paw_an(iatom)%vxc1,paw_an(iatom)%vxct1)
     if (paw_an(iatom)%has_vxcval>0) deallocate(paw_an(iatom)%vxc1_val,paw_an(iatom)%vxct1_val)
     if (paw_an(iatom)%has_vhartree>0) deallocate(paw_an(iatom)%vh1,paw_an(iatom)%vht1)
     if (paw_ij(iatom)%has_dijfr>0) deallocate(paw_ij(iatom)%dijfr)
     if (paw_ij(iatom)%has_dijhat>0) deallocate(paw_ij(iatom)%dijhat)
     if (paw_ij(iatom)%has_dijxc>0) deallocate(paw_ij(iatom)%dijxc)
     if (paw_ij(iatom)%has_dijxc_val>0) deallocate(paw_ij(iatom)%dijxc_val)
     if (paw_ij(iatom)%has_dijso>0) deallocate(paw_ij(iatom)%dijso)
     if (paw_ij(iatom)%has_dijU>0) deallocate(paw_ij(iatom)%dijU)
     if (pawtab(itypat)%usepawu>0) deallocate(paw_ij(iatom)%noccmmp,paw_ij(iatom)%nocctot)
     if (pawtab(itypat)%useexexch>0) deallocate(paw_ij(iatom)%vpawx)
     paw_an(iatom)%has_vxc=0;paw_an(iatom)%has_vxcval=0;paw_an(iatom)%has_vhartree=0
     paw_ij(iatom)%has_dijfr=0;paw_ij(iatom)%has_dijhat=0
     paw_ij(iatom)%has_dijxc=0;paw_ij(iatom)%has_dijxc_val=0
     paw_ij(iatom)%has_dijso=0;paw_ij(iatom)%has_dijU=0
     if (dtset%iscf>0) then
       pawrhoij(iatom)%lmnmix_sz=0
       pawrhoij(iatom)%use_rhoijres=0
       deallocate(pawrhoij(iatom)%kpawmix)
       deallocate(pawrhoij(iatom)%rhoijres)
     end if
   end do
   deallocate(pawfgrtab,paw_an,paw_ij)
   deallocate(dimcprj)
 end if

!Restore some variables in the dtset
!Here, for TDDFT, iprcel was artificially set.
!if(dtset%iscf==-1) then
!dtset%iprcel = dtset_iprcel
!end if

!===============================
!    chen: deallocate arrays ...
 if (psps%usepaw==1) then
   deallocate(extpot_lm)
   deallocate(extpot_dij)
   deallocate(hc_rho1)
   deallocate(hc_trho1)
 endif


 if (mpi_enreg%me == 0 .and. dtset%prtxml == 1) then
   write(ab_xml_out, "(A)") '      <finalConditions>'
!  We output the final result given in results_gs
   call out_resultsgs_XML(dtset, 4, results_gs, psps%usepaw)
   write(ab_xml_out, "(A)") '      </finalConditions>'
   write(ab_xml_out, "(A)") '    </scfcvLoop>'
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(60,2,tsec)
 call timab(20,2,tsec)

#if defined DEBUG_MODE
 write(message,'(a)')' scfcv : exit '
 call wrtout(std_out,message,'COLL')
 call flush_unit(std_out)
#endif

end subroutine emb_scfcv
!!***
