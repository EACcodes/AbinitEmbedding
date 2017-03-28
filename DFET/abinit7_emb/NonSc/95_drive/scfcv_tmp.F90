!{\src2tex{textfont=tt}}
!!****f* ABINIT/scfcv_tmp
!! NAME
!! scfcv_tmp
!!
!! FUNCTION
!! WARNING : Temporary wrapper to scfcv 
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
!!  nspinor=number of spinorial components of the wavefunctions
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
!!  scf_in <type(scf_in_type)>=intent(in) to scfcv, additional to the ones in dtset, dtfil, dtgs .
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points and spins
!!
!! SIDE EFFECTS
!!  acell(3)=length scales of primitive translations (bohr)
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

subroutine scfcv_tmp(acell,atindx,atindx1,cg,cpus,dtefield,dtfil,dtpawuj,&
&  dtset,ecore,eigen,electronpositron,hdr,iapp,indsym,initialized,&
&  irrzon,kg,mpi_enreg,nattyp,ndtpawuj,nfftf,npwarr,nspinor,occ,&
&  paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&  phnons,psps,pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
&  scf_history,scf_in,symrec,taug,taur,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)


 use defs_basis
 use defs_datatypes
 use m_wffile
 use defs_abitypes
 use defs_scftypes
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
 use interfaces_67_common
 use interfaces_79_seqpar_mpi
 use interfaces_95_drive, except_this_one => scfcv_tmp
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iapp,ndtpawuj,pwind_alloc
 integer,intent(inout) :: initialized,nfftf
 integer,intent(in) :: nspinor
 real(dp),intent(in) :: cpus,ecore
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(inout) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(efield_type),intent(inout) :: dtefield
 type(electronpositron_type),pointer :: electronpositron
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(inout) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(recursion_type),intent(inout) :: rec_set
 type(results_gs_type),intent(inout) :: results_gs
 type(scf_history_type),intent(inout) :: scf_history
 type(scf_in_type),intent(in) :: scf_in
 type(wffile_type),intent(inout) :: wffnew,wffnow
 type(wvl_data),intent(inout) :: wvl
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
 integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
!no_abirules
 integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  !(nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise)
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
 integer, intent(in) :: nattyp(psps%ntypat),npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
 integer, intent(inout) :: symrec(3,3,dtset%nsym)
 real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
 real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  !(nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise)
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
 real(dp), intent(inout) :: acell(3),rprimd(3,3)
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
!integer,parameter :: level=110,response=0
!integer :: accessfil,afford,choice
!integer :: computed_forces,cplex,dbl_nnsclo,dielop,dielstrt,dimdmat
!integer :: dtset_iprcel,fformr,forces_needed
!integer :: iatom,ider,idir,ierr,iexit,ii,ikpt,impose_dmat
!integer :: initialized0,iorder_cprj,ipert,ipositron,ir,isave,iscf10,ispden
!integer :: ispmix,istep,istep_mix,itypat,izero
!integer :: lm_size,lmax_diel,lmn2_size,lpawumax
!integer :: mgfftdiel,mgfftf,moved_atm_inside,moved_rhor,n1xccc
!integer :: n3xccc,n_fftgr,n_index,nele,nfftdiel,nfftmix,nhatgrdim,nk3xc,nkxc
!integer :: npawmix,npwdiel,nstep,nzlmopt,offset,optberry,optene,optgr0
!integer :: optgr1,optgr2,option,optrad,optres,optxc,prtden,prtkden,prtfor,quit
!integer :: quit_sum,rdwr,rdwrpaw,spaceComm,stress_needed,unit_out
!integer :: usecprj,usexcnhat,v_size
!real(dp) :: boxcut,compch_fft,compch_sph,deltae,diecut,diffor,ecut
!real(dp) :: ecutf,ecutsus,edum,elast,etotal,fermie,gsqcut
!real(dp) :: maxfor,res2,residm,ucvol,val_max
!real(dp) :: val_min,vxcavg,vxcavg_dum
!!logical :: ex
!character(len=4) :: tag
!character(len=500) :: message
!character(len=fnlen) :: fildata,kgnam
!type(MPI_type) :: mpi_enreg_diel
!type(energies_type) :: energies
!arrays
!integer :: ngfft(18),ngfftdiel(18),ngfftf(18),ngfftmix(18),npwarr_diel(1)
!integer :: npwtot_diel(1)
!integer,allocatable :: dimcprj(:),gbound_diel(:,:),i_rhor(:),i_vresid(:)
!integer,allocatable :: i_vrespc(:),i_vtrial(:),irrzondiel(:,:,:),kg_diel(:,:)
!integer,allocatable :: lmn_size(:)
!real(dp) :: dielar(7),dphase(3),dummy2(6),favg(3),gmet(3,3),gprimd(3,3),k0(3)
!real(dp) :: kpt_diel(3),pel(3),pel_cg(3),pelev_dum(3),pion(3),ptot(3)
!real(dp) :: rhodum(1),rmet(3,3),strsxc(6),strten(6),tollist(12)
!real(dp) :: tsec(2),vnew_mean(dtset%nspden),vres_mean(dtset%nspden)
!real(dp),allocatable :: dielinv(:,:,:,:,:),dtn_pc(:,:)
!real(dp),allocatable :: f_atm(:,:,:),f_fftgr(:,:,:),f_paw(:,:)
!real(dp),allocatable :: fcart(:,:),forold(:,:),fred(:,:),gresid(:,:)
!real(dp),allocatable :: grewtn(:,:),grhf(:,:),grnl(:),grxc(:,:)
!real(dp),allocatable :: kxc(:,:),nhat(:,:),nhatgr(:,:,:),nvresid(:,:)
!real(dp),allocatable :: ph1d(:,:),ph1ddiel(:,:),ph1df(:,:)
!real(dp),allocatable :: phnonsdiel(:,:,:),shiftvector(:)
!real(dp),allocatable :: susmat(:,:,:,:,:),synlgr(:,:)
!real(dp),allocatable :: vhartr(:),vpsp(:),vtrial(:,:)
!real(dp),allocatable :: vxc(:,:),workr(:,:),xccc3d(:),ylmdiel(:,:)
!real(dp),pointer :: elfr(:,:),grhor(:,:,:),lrhor(:,:)
!type(cprj_type),allocatable :: cprj(:,:)
!type(paw_an_type),allocatable :: paw_an(:)
!type(paw_ij_type),allocatable :: paw_ij(:)
!type(pawfgrtab_type),allocatable,save :: pawfgrtab(:)
! *********************************************************************

 if(.false.)write(6,*)nspinor

!WVL - reformat the wavefunctions in the case of xred != xred_old
 if (dtset%usewvl == 1 .and. maxval(xred_old - xred) > zero) then
!  WVL - Before running scfcv, on non-first geometry step iterations, we need
!  to reformat the wavefunctions, taking into acount the new
!  coordinates.
!  We prepare to change rhog (to be removed) and rhor.
   deallocate(rhog)
   deallocate(rhor)

   call wvl_wfsinp_reformat(acell, dtset, mpi_enreg, psps, &
&   rprimd, wvl, xred, xred_old)
   nfftf = dtset%nfft

   allocate(rhog(2, dtset%nfft))
   allocate(rhor(2, dtset%nfft))
   call wvl_mkrho(dtset, mpi_enreg, occ, rhor, wvl%wfs)
 end if

!Prepare the names of the auxiliary files whose name depend on the itimimage, iimage and itime loops.
 call dtfil_init2(dtfil,iapp,mpi_enreg)

 call scfcv(atindx,atindx1,cg,cpus,dtefield,dtfil,dtpawuj,&
& dtset,ecore,eigen,electronpositron,scf_in%fatvshift,hdr,iapp,indsym,initialized,&
& irrzon,kg,mpi_enreg,nattyp,ndtpawuj,nfftf,npwarr,occ,&
& paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,&
& phnons,psps,pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
& scf_history,symrec,taug,taur,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)

end subroutine scfcv_tmp
!!***
