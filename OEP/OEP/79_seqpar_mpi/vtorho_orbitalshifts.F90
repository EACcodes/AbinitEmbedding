module vtoorbitalshifts

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_profiling
 use m_xmpi
 use m_errors
 use m_wffile
 use m_efield
 use m_bandfft_kpt
#if defined HAVE_MPI2
 use mpi
#endif

 use m_header,   only : hdr_update, hdr_skip, hdr_io
 use m_pawrhoij, only : pawrhoij_type, rhoij_alloc, rhoij_free

 use m_energies,           only : energies_type
 use m_hamiltonian,        only : init_hamiltonian,destroy_hamiltonian,&
&                                 load_paw_hamiltonian,finalize_hamiltonian,gs_hamiltonian_type
 use m_electronpositron,   only : electronpositron_type,electronpositron_calctype
 use m_paw_dmft,           only : paw_dmft_type,init_dmft,destroy_dmft,print_dmft,saveocc_dmft
 use m_crystal,            only : init_crystal, destroy_crystal, crystal_structure
 use m_oper,               only : oper_type,init_oper,destroy_oper
 use m_io_tools,           only : flush_unit

 use get_xc_potential

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'vtorho'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_53_spacepar
 use interfaces_56_recipspace
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_62_occeig
 use interfaces_62_wvl_wfs
 use interfaces_65_nonlocal
 use interfaces_66_paw
 use interfaces_66_wfs
 use interfaces_67_common
 use interfaces_68_dmft
 use interfaces_77_suscep
 use interfaces_79_seqpar_mpi
!End of the abilint section

 use m_linalg_interfaces
 use m_cgtools
  
 use orbital_shifts 

contains

  !{\src2tex{textfont=tt}}
  !!****f* ABINIT/vtowfk
  !! NAME
  !! vtowfk
  !!
  !! FUNCTION
  !! This routine compute the partial density at a given k-point,
  !! for a given spin-polarization, from a fixed Hamiltonian
  !! but might also simply compute eigenvectors and eigenvalues at this k point
  !!
  !! COPYRIGHT
  !! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MT)
  !! This file is distributed under the terms of the
  !! GNU General Public License, see ~abinit/COPYING
  !! or http://www.gnu.org/copyleft/gpl.txt .
  !!
  !! INPUTS
  !!  cpus                  = cpu time limit in seconds
  !!  dimffnl=second dimension of ffnl (1+number of derivatives)
  !!  dtfil <type(datafiles_type)>=variables related to files
  !!  dtset <type(dataset_type)>=all input variables for this dataset
  !!  fixed_occ=true if electronic occupations are fixed (occopt<3)
  !!  ffnl(npw_k,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
  !!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
  !!  icg=shift to be applied on the location of data in the array cg
  !!  ikpt=number of the k-point
  !!  iscf=(<= 0  =>non-SCF), >0 => SCF
  !!  isppol isppol=1 for unpolarized, 2 for spin-polarized
  !!  kg_k(3,npw_k)=reduced planewave coordinates.
  !!  kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
  !!  kpg_k(npw,nkpg)= (k+G) components (only if useylm=1)
  !!  mcg=second dimension of the cg array
  !!  mkgq = second dimension of pwnsfacq
  !!  mpi_enreg=informations about MPI parallelization
  !!  mpw=maximum dimensioned size of npw
  !!  natom=number of atoms in cell.
  !!  nband_k=number of bands at this k point for that spin polarization
  !!  nkpg=second dimension of kpg_k (0 if useylm=0)
  !!  nkpt=number of k points.
  !!  nnsclo_now=number of non-self-consistent loops for the current vtrial
  !!             (often 1 for SCF calculation, =nstep for non-SCF calculations)
  !!  npw_k=number of plane waves at this k point
  !!  npwarr(nkpt)=number of planewaves in basis at this k point
  !!  ntypat=number of types of atoms in unit cell.
  !!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
  !!  optforces=option for the computation of forces
  !!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
  !!  prtvol=control print volume and debugging output
  !!  pwind(pwind_alloc,2,3)= array used to compute
  !!           the overlap matrix smat between k-points (see initberry.f)
  !!  pwind_alloc= first dimension of pwind
  !!  pwnsfac(2,pwind_alloc)= phase factors for non-symmorphic translations
  !!                          (see initberry.f)
  !!  pwnsfacq(2,mkgq)= phase factors for the nearest neighbours of the
  !!                    current k-point (electric field, MPI //)
  !!  usebanfft=flag for band-fft parallelism
  !!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
  !!  vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
  !!
  !! OUTPUT
  !!  dphase_k(3)=change in Zak phase for the current k-point
  !!  eig_k(nband_k)=array for holding eigenvalues (hartree)
  !!  ek_k(nband_k)=contribution from each band to kinetic energy, at this k-point
  !!  ek_k_nd(2,nband_k,nband_k*use_dmft)=contribution to kinetic energy, including non-diagonal terms, at this k-point (usefull if use_dmft)
  !!  resid_k(nband_k)=residuals for each band over all k points,
  !!                   BEFORE the band rotation.
  !!  ==== if optforces>0 ====
  !!    grnl_k(3*natom,nband_k)=nonlocal gradients, at this k-point
  !!  ==== if (gs_hamk%usepaw==0) ====
  !!    enl_k(nband_k)=contribution from each band to nonlocal pseudopotential part of total energy, at this k-point
  !!  ==== if (gs_hamk%usepaw==1) ====
  !!
  !! SIDE EFFECTS
  !!  cg(2,mcg)=updated wavefunctions
  !!  rhoaug(n4,n5,n6,nvloc)= density in electrons/bohr**3, on the augmented fft grid.
  !!                    (cumulative, so input as well as output). Update only
  !!                    for occopt<3 (fixed occupation numbers)
  !!
  !! PARENTS
  !!      vtorho
  !!
  !! CHILDREN
  !!      cgwf,dsymm,fourwf,fxphas,lobpcgcciiwf
  !!      lobpcgiiwf,lobpcgwf,meanvalue_g,nonlop,prep_fourwf,prep_nonlop
  !!      pw_orthon,subdiago,timab,wrtout,xsum_mpi,zhemm
  !!
  !! NOTES
  !!  Only the mod((iband-1)/mpi_enreg%bandpp,mpi_enreg%nproc_band) projectors
  !!  are stored on each proc.
  !!
  !! SOURCE
  
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
  
#include "abi_common.h"
  
  !{\src2tex{textfont=tt}}
  !!****f* ABINIT/vtorho
  !! NAME
  !! vtorho
  !!
  !! FUNCTION
  !! This routine compute the new density from a fixed potential (vtrial)
  !! but might also simply compute eigenvectors and eigenvalues.
  !! The main part of it is a wf update over all k points.
  !!
  !! COPYRIGHT
  !! Copyright (C) 1998-2012 ABINIT group (DCA, XG, GMR, MF, AR, MM, MT, FJ, MB, MT)
  !! This file is distributed under the terms of the
  !! GNU General Public License, see ~abinit/COPYING
  !! or http://www.gnu.org/copyleft/gpl.txt .
  !! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
  !!
  !! INPUTS
  !!  afford=used to dimension susmat
  !!  atindx(natom)=index table for atoms (see gstate.f)
  !!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
  !!  cpus= cpu time limit in seconds
  !!  dbl_nnsclo=if 1, will double the value of dtset%nnsclo
  !!  dielop= if positive, the dielectric matrix must be computed.
  !!  dielstrt=number of the step at which the dielectric preconditioning begins.
  !!  dtfil <type(datafiles_type)>=variables related to files
  !!  dtset <type(dataset_type)>=all input variables for this dataset
  !!   | mband=maximum number of bands
  !!   | mgfft=maximum size of 1D FFTs
  !!   | mkmem =number of k points which can fit in memory; set to 0 if use disk
  !!   | mpw=maximum dimensioned size of npw
  !!   | nfft=(effective) number of FFT grid points (for this processor)
  !!   | nkpt=number of k points.
  !!   | nspden=number of spin-density components
  !!   | nsppol=1 for unpolarized, 2 for spin-polarized
  !!   | nsym=number of symmetry elements in space group
  !!   | typat= array of types of the natoms
  !!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
  !!  etotal=total energy (Ha) - only needed for tddft
  !!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for the dielectric matrix
  !!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
  !!  gprimd(3,3)=dimensional reciprocal space primitive translations
  !!   (3x3 tensor) and grads wrt atomic coordinates (3*natom)
  !!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
  !!  hdr <type(hdr_type)>=the header of wf, den and pot files
  !!  indsym(4,nsym,natom)=indirect indexing array for atom labels
  !!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
  !!  irrzondiel(nfftdiel**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data for diel matrix
  !!                                     nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
  !!  istep=index of the number of steps in the routine scfcv
  !!  istep_mix=index of the number of steps for the SCF mixing (can be <istep)
  !!  kg(3,mpw*mkmem)=reduced planewave coordinates.
  !!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
  !!  kxc(nfftf,nkxc)=exchange-correlation kernel, needed only if nkxc/=0 .
  !!  lmax_diel=1+max. value of l angular momentum used for dielectric matrix
  !!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
  !!  mgfftdiel=maximum size of 1D FFTs, for the computation of the dielectric matrix
  !!  mpi_enreg=informations about MPI parallelization
  !!  my_natom=number of atoms treated by current processor
  !!  natom=number of atoms in cell.
  !!  nattyp(ntypat)= # atoms of each type.
  !!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
  !!         (nfftf=nfft for norm-conserving potential runs)
  !!  nfftdiel=number of fft grid points for the computation of the diel matrix
  !!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
  !!                see ~abinit/doc/input_variables/vargs.htm#ngfft
  !!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
  !!  npwarr(nkpt)=number of planewaves in basis at this k point
  !!  ntypat=number of types of atoms in unit cell.
  !!  optforces=option for the computation of forces (0: no force;1: forces)
  !!  optres=0: the new value of the density is computed in place of the input value
  !!         1: only the density residual is computed ; the input density is kept
  !!  paw_dmft  <type(paw_dmft_type)>= paw+dmft related data
  !!  paw_ij(my_natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
  !!  pawang <type(pawang)>=paw angular mesh and related data
  !!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
  !!  pawfgrtab(my_natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
  !!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
  !!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
  !!                                    nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
  !!  phnonsdiel(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases,
  !!   for diel matr
  !!                                     nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
  !!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
  !!  ph1ddiel(2,3*(2*mgfftdiel+1)*natom*usepaw)=one-dimensional structure factor information
  !!                                             for the dielectric matrix
  !!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
  !!  pwind(pwind_alloc,2,3) = array used to compute
  !!           the overlap matrix smat between k-points (see initberry.f)
  !!  pwind_alloc = first dimension of pwind
  !!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
  !!                           (see initberry.f)
  !!  rmet(3,3)=real space metric (bohr**2)
  !!  rprimd(3,3)=dimensional primitive vectors
  !!  symrec(3,3,nsym)=symmetry operations in reciprocal space
  !!  ucvol=unit cell volume in bohr**3.
  !!  wffnew,wffnow=unit numbers for wf disk files.
  !!  vtrial(nfftf,nspden)=INPUT potential Vtrial(r).
  !!  xred(3,natom)=reduced dimensionless atomic coordinates
  !!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
  !!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
  !!  ylmdiel(npwdiel,lmax_diel**2)= real spherical harmonics for each G and k point
  !!                                 for the dielectric matrix
  !!
  !! OUTPUT
  !!  dphase(3) : dphase(idir) = accumulated change in the string-averaged
  !!     Zak phase along the idir-th direction caused by the update of all
  !!     the occupied Bloch states at all the k-points (only if finite electric field)
  !!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
  !!  resid(mband*nkpt*nsppol)=residuals for each band over all k points.
  !!  residm=maximum value from resid array (except for nbdbuf highest bands)
  !!   the susceptibility (or density-density response) matrix in reciprocal space
  !!  === if optforces>0 ===
  !!    grnl(3*natom)=stores grads of nonlocal energy wrt length scales
  !!  ==== if optres==1
  !!    nres2=square of the norm of the residual
  !!    nvresid(nfftf,nspden)=density residual
  !!  ==== if psps%usepaw==1
  !!
  !! SIDE EFFECTS
  !!  cg(2,mpw*dtset%nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
  !!   At output contains updated wavefunctions coefficients;
  !!    if nkpt>1, these are kept in a disk file.
  !!  energies <type(energies_type)>=storage for energies computed here :
  !!   | e_eigenvalues=Sum of the eigenvalues - Band energy (Hartree)
  !!   | e_kinetic=kinetic energy part of total energy
  !!   | e_nonlocalpsp=nonlocal pseudopotential part of total energy
  !!   | e_fermie=fermi energy (Hartree)
  !!  occ(mband*nkpt*nsppol)=occupation number for each band for each k.
  !!      (input if insulator - occopt<3 - ; output if metallic)
  !!  pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
  !!  rhog(2,nfftf)=Fourier transform of total electron density
  !!  rhor(nfftf,nspden)=total electron density (el/bohr**3)
  !!  taug(2,nfftf*dtset%usekden)=array for Fourier transform of kinetic energy density
  !!  taur(nfftf,nspden*dtset%usekden)=array for kinetic energy density
  !!  wvl <type(wvl_data)>=wavelets structures in case of wavelets basis.
  !!
  !! PARENTS
  !!      scfcv
  !!
  !! CHILDREN
  !!      calcdensph,clsopn,datafordmft,destroy_crystal
  !!      destroy_dmft,destroy_hamiltonian,destroy_oper,dmft_solve
  !!      eigensystem_info,evaltoocc,fftpac,finalize_hamiltonian,flush_unit
  !!      gpu_finalize_ffnl_ph3d,gpu_finalize_ham_data,gpu_update_ffnl_ph3d
  !!      hdr_io,hdr_skip,hdr_update,init_crystal,init_dmft,init_hamiltonian
  !!      init_oper,last_orthon,leave_test,load_paw_hamiltonian,mag_loc_k,magcart
  !!      mkffnl,mkkin,mkkpg,mkrho,newocc,pawmkrho,pawmkrhoij,ph1d3d
  !!      prep_bandfft_tabs,print_dmft,prteigrs,prtrhomxmn,rdnpw,rhoij_alloc
  !!      rhoij_free,rwwf,saveocc_dmft,setnoccmmp,sphereboundary,sqnorm_v,status
  !!      suscep_stat,symrhg,tddft,testsusmat,timab,transgrid,update_mmat,vtowfk
  !!      wffkg,write_energies,wrtout,wvl_hpsitopsi,wvl_nl_gradient,wvl_psitohpsi
  !!      xallgather_mpi,xbarrier_mpi,xdefineoff,xmax_mpi,xrecv_mpi,xred2xcart
  !!      xsend_mpi,xsum_mpi
  !!
  !! NOTES
  !!  Be careful to the meaning of nfft (size of FFT grids):
  !!   - In case of norm-conserving calculations the FFT grid is the usual FFT grid.
  !!   - In case of PAW calculations:
  !!     Two FFT grids are used; one with nfft points (coarse grid) for
  !!     the computation of wave functions ; one with nfftf points
  !!     (fine grid) for the computation of total density.
  !!
  !!  The total electronic density (rhor,rhog) is divided into two terms:
  !!   - The density related to WFs =Sum[Psi**2]
  !!
  !!  The parallelisation needed for the electric field should be
  !!  made an independent subroutine, so that this routine could be put
  !!  back in the 95_drive directory.
  !!
  !! SOURCE
  
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
  
#include "abi_common.h"
  
  subroutine vtorho_orbitalshifts(afford,atindx,atindx1,cg,&
  &           cpus,dbl_nnsclo, &
  &           dielop,dielstrt,dphase,dtfil,dtset,&
  &           eigen,electronpositron,energies,etotal,gbound_diel,&
  &           gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
  &           istep,istep_mix,kg,kg_diel,kxc,lmax_diel,&
  &           mcg,mgfftdiel,mpi_enreg,&
  &           my_natom,natom,nattyp,nfftf,nfftdiel,ngfftdiel,nkxc,&
  &           npwarr,npwdiel,nres2,ntypat,nvresid,occ,optforces,&
  &           optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,&
  &           phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
  &           pwind,pwind_alloc,pwnsfac,resid,residm,rhog,rhor,&
  &           rmet,rprimd,run_params,symrec,taug,taur,&
  &           ucvol,wffnew,wffnow,vhartr_new, vhartr_old,vtrial,&
  &           wvl,xred,ylm,ylmgr,ylmdiel)
  
#if defined HAVE_DFT_BIGDFT
  use BigDFT_API, only : last_orthon, SMEARING_DIST_ERF, evaltoocc
#endif

   implicit none
  
#if defined HAVE_MPI1
   include 'mpif.h'
#endif
  
  !Arguments -------------------------------
   integer, intent(in) :: afford,dbl_nnsclo,dielop,dielstrt,istep,istep_mix,lmax_diel,mcg,mgfftdiel
   integer, intent(in) :: my_natom,natom,nfftf,nfftdiel,nkxc,npwdiel
   integer, intent(in) :: ntypat,optforces,optres,pwind_alloc
   real(dp), intent(in) :: cpus,etotal,gsqcut,ucvol
   real(dp), intent(out) :: nres2,residm
   type(MPI_type), intent(inout) :: mpi_enreg
   type(datafiles_type), intent(in) :: dtfil
   type(dataset_type), intent(in) :: dtset
   type(electronpositron_type),pointer :: electronpositron
   type(energies_type), intent(inout) :: energies
   type(hdr_type), intent(inout) :: hdr
   type(paw_dmft_type), intent(inout)  :: paw_dmft
   type(pawang_type), intent(in) :: pawang
   type(pawfgr_type), intent(in) :: pawfgr
   type(pseudopotential_type), intent(in) :: psps
   type(wffile_type), intent(inout) :: wffnew,wffnow
   type(wvl_data), intent(inout) :: wvl

   type (run_parameter_type) run_params

   real(dp) vhartr_new(nfftf),vhartr_old(nfftf)

   integer, intent(in) :: atindx(natom),atindx1(natom),gbound_diel(2*mgfftdiel+8,2)
   integer, intent(in) :: indsym(4,dtset%nsym,natom)
   integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
   integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
   integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem),kg_diel(3,npwdiel),nattyp(ntypat),ngfftdiel(18),npwarr(dtset%nkpt)
   integer, intent(in) :: pwind(pwind_alloc,2,3),symrec(3,3,dtset%nsym)
   real(dp), intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*dtset%mgfft+1)*natom)
   real(dp), intent(in) :: ph1ddiel(2,(3*(2*mgfftdiel+1)*natom)*psps%usepaw)
   real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
   real(dp), intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
   real(dp), intent(in) :: pwnsfac(2,pwind_alloc),rmet(3,3),rprimd(3,3)
   real(dp), intent(inout) :: vtrial(nfftf,dtset%nspden)
   real(dp), intent(inout) :: xred(3,natom)
   real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
   real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
   real(dp), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
   real(dp), intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
   real(dp), intent(out) :: dphase(3),grnl(3*natom)
   real(dp), intent(out) :: nvresid(nfftf,dtset%nspden),resid(dtset%mband*dtset%nkpt*dtset%nsppol)
   real(dp), intent(in) :: cg(2,mcg)
   real(dp), intent(inout) :: kxc(nfftf,nkxc)
   real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
   real(dp), intent(inout) :: rhog(2,nfftf),rhor(nfftf,dtset%nspden)
   real(dp), intent(inout) :: taug(2,nfftf*dtset%usekden),taur(nfftf,dtset%nspden*dtset%usekden)
   type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
   type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
   type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*psps%usepaw)
   type(pawtab_type),intent(in)  :: pawtab(ntypat*psps%usepaw)
  
  !Local variables-------------------------------
   integer,parameter :: level=111,tim_mkrho=2
   integer,save :: nwarning=0
   integer :: bantot,bdtot_index,count,count1,counter
   integer :: cplex,dest,dimffnl,enunit
  ! integer :: dimenl1 ! used in magnetization
   integer :: fform,formeig,i1,i2,i3,ia,iatom,iband,iband1,ibdkpt
   integer :: ibg,icg,icg1,icg2,icp1,icp2,ider,idir
   integer :: ierr,iexit,ifft,ifor,ifor1,ii,ikg,ikg1,ikg2,ikpt
   integer :: ikptf,ikpt1f,ikpt1i
   integer :: ikpt_loc,ikpt1,ikpt_this_proc,ikxc,ilm,index1,imagn,ipert,iplex
   integer :: iproc,ir,iscf,ispden,isppol,istwf_k,itypat,jkpt
   integer :: jsppol,mbdkpsp,mb2dkpsp
   integer :: mcg_disk,me_distrb,mkgq,muig
   integer :: mwarning,my_nspinor,n1,n2,n3,n4,n5,n6,nband_eff
   integer :: nband_k,nbuf,ndatarecv,neglect_pawhat,nfftot,nkpg,nkpt1,nnn,nnsclo_now
   integer :: nproc_distrb,npw_k,npw_k1,nsp,nspden_rhoij,option,prtvol
   integer :: rdwr,spaceComm_distrb,tag,tim_rwwf,usetimerev
   integer :: my_source, his_source, jkptf, jkpt1f, jkpt1i
   logical :: berryflag,computesusmat,fixed_occ
   logical :: locc_test,remove_inv
   real(dp) :: dmft_ldaocc,dummy
   real(dp) :: edmft,ebandlda,ebanddmft,ebandldatot,ekindmft,ekindmft2,ekinlda
   real(dp) :: emax,min_occ,vxcavg_dum,strsxc(6)
   character(len=500) :: message
   type(gs_hamiltonian_type) :: gs_hamk
   type(wffile_type) :: wfftmp
   integer,allocatable :: ikpt_recv(:),kg_dum(:,:),kg_k(:,:)
   integer,allocatable :: flag_send(:,:), flag_receive(:)
   real(dp) :: dielar(7),dphase_k(3),kpoint(3),qpt(3),rhodum(1),tsec(2),ylmgr_dum(1)
   real(dp),allocatable :: EigMin(:,:),buffer(:,:),buffer1(:),buffer2(:)
   real(dp),allocatable :: cg_disk(:,:),cgrkxc(:,:),cgrvtrial(:,:),doccde(:)
   real(dp),allocatable :: dphasek(:,:),eig_dum(:),ek_k(:),ek_k_nd(:,:,:),eknk(:),eknk_nd(:,:,:,:,:)
   real(dp),allocatable :: enl_k(:),enlnk(:),ffnl(:,:,:,:),grnl_k(:,:), xcart(:,:)
   real(dp),allocatable :: grnlnk(:,:),kinpw(:),kpg_k(:,:),occ_dum(:),occ_k(:),ph3d(:,:,:)
   real(dp),allocatable :: pwnsfacq(:,:),resid_k(:),rhoaug(:,:,:,:),rhowfg(:,:),rhowfr(:,:)
   real(dp),allocatable :: vlocal(:,:,:,:),vlocal_tmp(:,:,:)
   real(dp),allocatable :: zshift(:),ylm_k(:,:)
   type(oper_type) :: lda_occup
   type(pawrhoij_type),allocatable :: pawrhoij_unsym(:)
   type(crystal_structure) :: cryst_struc
   integer,allocatable :: idum1(:),idum3(:,:,:)
   real(dp),allocatable :: rdum2(:,:),rdum4(:,:,:,:)


   integer wfoptalg, wfopta10
  
  ! *********************************************************************
  
   DBG_ENTER("COLL")
  
  !Keep track of total time spent in vtorho
   call timab(980,1,tsec)
   call timab(981,1,tsec)
  
  !Structured debugging if prtvol==-level
   prtvol=dtset%prtvol
   if(prtvol==-level)then
     write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' vtorho : enter '
     call wrtout(std_out,message,'COLL')
   end if
  
   call status(0,dtfil%filstat,iexit,level,'allocate/init ')
  
  !Init MPI for kpt communicator (.i.e. communicator for that k-point)
   spaceComm_distrb=mpi_enreg%comm_cell
   me_distrb=xcomm_rank(spaceComm_distrb)
   nproc_distrb=xcomm_size(spaceComm_distrb)
   if ((xmpi_paral==1).and.(mpi_enreg%paral_kgb==1)) me_distrb=mpi_enreg%me_kpt
   if (mpi_enreg%me_img/=0) nwarning=mwarning+1
  
  !Test size of FFT grids (1 grid in norm-conserving, 2 grids in PAW)
   if ((psps%usepaw==1.and.pawfgr%nfft/=nfftf).or.(psps%usepaw==0.and.dtset%nfft/=nfftf)) then
     MSG_BUG('wrong values for nfft, nfftf!')
   end if
  
  !Test optforces (to prevent memory overflow)
   if (optforces/=0.and.optforces/=1) then
     write(message,'(a,i0)')' wrong value for optforces = ',optforces
     MSG_BUG(message)
   end if
  
  !Debugging : print vtrial and rhor
  !MPIWF Warning : this should not be parallelized over space, leave this debugging feature as such.
   if(prtvol==-level)then
     if (psps%usepaw==0) then
       n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
     else
       n1=pawfgr%ngfft(1) ; n2=pawfgr%ngfft(2) ; n3=pawfgr%ngfft(3)
     end if
     write(message,'(a)') '   ir              vtrial(ir)     rhor(ir) '
     call wrtout(std_out,message,'COLL')
     do ir=1,nfftf
       i3=(ir-1)/n1/n2
       i2=(ir-1-i3*n1*n2)/n1
       i1=ir-1-i3*n1*n2-i2*n1
       write(message,'(i5,3i3,a,2es13.6)')ir,i1,i2,i3,' ',vtrial(ir,1),rhor(ir,1)
       call wrtout(std_out,message,'COLL')
       if(dtset%nspden>=2)then
         write(message,'(a,2es13.6)')'               ',vtrial(ir,2),rhor(ir,2)
         call wrtout(std_out,message,'COLL')
       end if
     end do
   end if
  
  !WVL - Branching with a separate vtorho procedure in wavelet. Should be merge in vtorho later on.
   if (dtset%usewvl == 1) then
     ABI_ALLOCATE(xcart,(3, dtset%natom))
     call xred2xcart(dtset%natom, rprimd, xcart, xred)
  !  This loop is for diaganolisation scheme.
     do ii = 1, dtset%nnsclo - 1, 1
       call wvl_hpsitopsi(dtset, energies, ii, mpi_enreg, residm, wvl)
       call wvl_psitohpsi(dtset%diemix, energies%e_exactX, energies%e_xc, &
  &     energies%e_hartree, energies%e_kinetic, energies%e_localpsp, &
  &     energies%e_nonlocalpsp, energies%e_sicdc, istep, ii, dtset%iscf, &
  &     dtset%ixc, me_distrb, dtset%natom, dtset%nfft, nproc_distrb, dtset%nspden, &
  &     nres2, .false., energies%e_vxc, wvl, xcart, strsxc)
     end do
     if (dtset%nnsclo > 0) then
#if defined HAVE_DFT_BIGDFT
       call write_energies(ii,0,wvl%energs,0.d0,0.d0,"final")
       call last_orthon(me_distrb, nproc_distrb, ii, wvl%wfs%ks, wvl%energs%evsum, .true.)
       call evaltoocc(me_distrb, nproc_distrb, .false., dtset%tsmear, wvl%wfs%ks%orbs, &
  &     SMEARING_DIST_ERF)
       energies%e_fermie = wvl%wfs%ks%orbs%efermi
       energies%e_eigenvalues = energies%e_kinetic + energies%e_localpsp + &
  &     energies%e_nonlocalpsp
       energies%e_xcdc = zero
       call eigensystem_info(me_distrb, nproc_distrb,0.d0,&
       wvl%wfs%ks%Lzd%Glr%wfd%nvctr_c+7*wvl%wfs%ks%Lzd%Glr%wfd%nvctr_f,&
       wvl%wfs%ks%orbs,wvl%wfs%ks%psi)
#else
       write(message, '(a,a,a,a)' ) ch10,&
  &     ' rhotov: BigDFT library is not compiled.', ch10, &
  &     '   Action, used the flag --enable-bigdft when configuring.'
       MSG_ERROR(message)
#endif
     else
       call wvl_hpsitopsi(dtset, energies, istep, mpi_enreg, residm, wvl)
     end if
     if (optforces == 1) then
       call wvl_nl_gradient(grnl, mpi_enreg, dtset%natom, rprimd, wvl, xcart)
     end if
     ABI_DEALLOCATE(xcart)
     call timab(980,2,tsec)
     return
   end if
  
  !WVL - Following is done in plane waves.
   iscf=dtset%iscf
   n1=dtset%ngfft(1); n2=dtset%ngfft(2); n3=dtset%ngfft(3)
   fixed_occ=(dtset%occopt<3.or.electronpositron_calctype(electronpositron)==1)
  
   energies%e_eigenvalues = zero
   energies%e_kinetic     = zero
   energies%e_nonlocalpsp = zero
   grnl(:)=zero
   resid(:) = zero ! JWZ 13 May 2010. resid and eigen need to be fully zeroed each time before use
   bdtot_index=0
   ibg=0;icg=0
   mbdkpsp=dtset%mband*dtset%nkpt*dtset%nsppol
   if(paw_dmft%use_dmft==1) mb2dkpsp=2*dtset%mband*dtset%mband*dtset%nkpt*dtset%nsppol
  
   ABI_ALLOCATE(eknk,(mbdkpsp))
   ABI_ALLOCATE(kg_k,(3,dtset%mpw))
   ABI_ALLOCATE(eknk_nd,(dtset%nsppol,dtset%nkpt,2,dtset%mband,dtset%mband*paw_dmft%use_dmft))
   ABI_ALLOCATE(EigMin,(2,dtset%mband))
   ABI_ALLOCATE(grnlnk,(3*natom,mbdkpsp*optforces))
   if (psps%usepaw==0)  then
     ABI_ALLOCATE(enlnk,(mbdkpsp))
   end if
  
   if(paw_dmft%use_dmft==1) eknk_nd=zero
   eknk(:)=zero;if (optforces>0) grnlnk(:,:)=zero
   if (psps%usepaw==0) enlnk(:)=zero
  
  !Initialize rhor if needed; store old rhor
   if(iscf>0 .or. iscf==-3) then
     if (optres==1) nvresid=rhor
     if (psps%usepaw==0) then
       rhor=zero
     else
       ABI_ALLOCATE(rhowfr,(dtset%nfft,dtset%nspden))
       ABI_ALLOCATE(rhowfg,(2,dtset%nfft))
       rhowfr(:,:)=zero
     end if
   end if
  
  !Set max number of non-self-consistent loops nnsclo_now for use in vtowfk
   if(iscf<=0)then
     nnsclo_now=dtset%nstep
   else if(iscf>0)then
     if(dtset%nnsclo>0) then
       nnsclo_now=dtset%nnsclo
     else if(dtset%nnsclo<=0)then
       nnsclo_now=1
       if(istep<=2)nnsclo_now=2
     end if
     if(dbl_nnsclo==1)then
  !    DEBUG
  !    write(std_out,*)' vtorho : use doubled nnsclo '
  !    ENDDEBUG
       nnsclo_now=nnsclo_now*2
     end if
   end if
   if(dtset%wfoptalg==2)nnsclo_now=40  ! UNDER DEVELOPMENT
  
   write(message, '(a,i3,a,i3,i2,i3)' ) ' vtorho : nnsclo_now=',nnsclo_now,&
  & ', note that nnsclo,dbl_nnsclo,istep=',dtset%nnsclo,dbl_nnsclo,istep
   call wrtout(std_out,message,'COLL')
  
  
   n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)
   my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
  
  
  
  !Prepare wf files for reading if dtset%mkmem==0
   if (dtset%mkmem==0) then
  
  !  Close files, and then reopen them
  !  (this is supposedly helpful for use of networked workstations
  !  and also sets up for later addition of a checkpoint facility
  !  for restarting crashed jobs)
  !  clsopn automatically checks to see whether file is scratch
  !  file and if so, does not close and open it.
     call clsopn(wffnow)
  
  
  
  !  Read wffnow header
     call hdr_skip(wffnow,ierr)
  
  
  !  Define offsets, in case of MPI I/O
     formeig=0
     call xdefineOff(formeig,wffnow,mpi_enreg,dtset%nband,npwarr,dtset%nspinor,dtset%nsppol,dtset%nkpt)
  
  
     call clsopn(wffnew)
  
  
  !  Update the content of the header (evolving variables)
     bantot=hdr%bantot ; dummy=1.0d20
  !  WARNING! fermie is used before set.
     call hdr_update(bantot,dummy,energies%e_fermie,hdr,dtset%natom,&
  &   residm,rprimd,occ,pawrhoij,psps%usepaw,xred,&
  &   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
  
  
  !  Write the content of hdr to the new wf file
     rdwr=2 ; fform=2
     if (wffnew%accesswff /= IO_MODE_NETCDF) then
       call hdr_io(fform,hdr,rdwr,wffnew)
  !    XG120513 Badly resolved conflict ?
  !    #if defined HAVE_TRIO_NETCDF
  !    else if (wffnew%accesswff == IO_MODE_NETCDF) then
  !    call hdr_io_netcdf(fform,hdr,rdwr,wffnew)
  !    
  !    call ini_wf_netcdf(dtset%mpw,wffnew%unwff,0)
  !    #endif
     else
       write(message,'(a,i0,a)')"accesswff = ",wffnew%accesswff," not coded"
       MSG_ERROR(message)
     end if
  
  
  !  Define offsets, in case of MPI I/O
     formeig=0
     call WffKg(wffnew,1)
     call xdefineOff(formeig,wffnew,mpi_enreg,dtset%nband,npwarr,dtset%nspinor,dtset%nsppol,dtset%nkpt)
  
     mcg_disk=dtset%mpw*my_nspinor*dtset%mband
     ABI_ALLOCATE(cg_disk,(2,mcg_disk))
   end if
  
  
  
  !============================================
  !==== Initialize most of the Hamiltonian ====
  !============================================
  !1) Allocate all arrays and initialize quantities that do not depend on k and spin.
  !2) Perform the setup needed for the non-local factors:
  !* Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
  !* PAW: Initialize the overlap coefficients and allocate the Dij coefficients.
  
   call init_hamiltonian(gs_hamk,psps,pawtab,dtset%nspinor,dtset%nspden,natom,ntypat,&
  & dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,&
  & ph1d=ph1d,electronpositron=electronpositron,use_gpu_cuda=dtset%use_gpu_cuda)
  
  !Electric and magnetic fields: set flag to turn on various behaviors
   berryflag = .FALSE.
   if (dtset%berryopt == 4 .or. (abs(dtset%berryopt) == 5) .or. &
  & dtset%berryopt == 6 .or. dtset%berryopt == 7 .or.  &
  & dtset%berryopt == 14 .or. dtset%berryopt == 16 .or. dtset%berryopt == 17) berryflag = .TRUE.    !!HONG
  
   nkpt1 = dtset%nkpt
  
   ikpt_loc = 0
  
   ABI_ALLOCATE(rhoaug,(n4,n5,n6,gs_hamk%nvloc))
   ABI_ALLOCATE(vlocal,(n4,n5,n6,gs_hamk%nvloc))
   rhoaug=zero
  
   call timab(981,2,tsec)
  
  !LOOP OVER SPINS
   do isppol=1,dtset%nsppol
  
     call timab(982,1,tsec)
  
     if (dtset%nsppol==2) then
       write(message,*)' ****  In vtorho for isppol=',isppol
       call wrtout(std_out,message,'COLL')
     end if
  
  !  original line
  !  if ((mpi_enreg%paral_compil_kpt == 0).or.(.not.berryflag)) ikpt_loc = 0
     ikpt_loc = 0
  
  !  DEBUG feature
  !  write(std_out,'(a,i4,L8,i4)')' debug--paral_compil_kpt, berryflag, ikpt_loc : ',mpi_enreg%paral_compil_kpt,&
  !  &        berryflag,ikpt_loc
  !  END DEBUG feature
  
  !  Rewind kpgsph data file if needed:
     if (dtset%mkmem==0) rewind dtfil%unkg
     if (dtset%mkmem==0.and.psps%useylm==1) rewind dtfil%unylm
     ikg=0
  
  !  Set up local potential vlocal with proper dimensioning, from vtrial
  !  Also take into account the spin.
     if(dtset%nspden/=4)then
       if (psps%usepaw==0.or.pawfgr%usefinegrid==0) then
         call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,vtrial,vlocal,2)
       else
         ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
         call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
         call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal,2)
         ABI_DEALLOCATE(cgrvtrial)
       end if
     else
       ABI_ALLOCATE(vlocal_tmp,(n4,n5,n6))
       if (psps%usepaw==0.or.pawfgr%usefinegrid==0) then
         do ispden=1,dtset%nspden
           call fftpac(ispden,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,vtrial,vlocal_tmp,2)
           vlocal(:,:,:,ispden)=vlocal_tmp(:,:,:)
         end do
       else
         ABI_ALLOCATE(cgrvtrial,(dtset%nfft,dtset%nspden))
         call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
         do ispden=1,dtset%nspden
           call fftpac(ispden,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal_tmp,2)
           vlocal(:,:,:,ispden)=vlocal_tmp(:,:,:)
         end do
         ABI_DEALLOCATE(cgrvtrial)
       end if
       ABI_DEALLOCATE(vlocal_tmp)
     end if
     rhoaug(:,:,:,:)=zero
  
  !  Continue to initialize the Hamiltonian
     call load_paw_hamiltonian(gs_hamk,isppol,paw_ij,&
  &   mpi_comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
  
     call timab(982,2,tsec)
  
  !  BIG FAT k POINT LOOP
  !  MVeithen: I had to modify the structure of this loop in order to implement MPI // of the electric field
  
  !  note that the loop here differs from the similar one in berryphase_new.F90.
  !  here, ikpt_loc numbers the kpts treated by the current processor.
  !  in berryphase_new.F90, ikpt_loc ALSO includes info about value of isppol.
  
     ikpt = 0
  
  
  
     do while (ikpt_loc < nkpt1)
  
       call timab(997,1,tsec)
  
       if (.not.berryflag) then
         ikpt_loc = ikpt_loc + 1
         ikpt = ikpt_loc
       else
         if (ikpt_loc < dtset%mkmem) ikpt = ikpt + 1
         if ((ikpt > dtset%nkpt).and.(ikpt_loc < dtset%mkmem)) exit
       end if
  
       dphase_k(:) = zero
       counter=100*ikpt+isppol
       call status(counter,dtfil%filstat,iexit,level,'loop ikpt     ')
       nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
       istwf_k=dtset%istwfk(ikpt)
       npw_k=npwarr(ikpt)
       mkgq = 1
  
       if(xmpi_paral==1)then
  
          if(proc_distrb_cycle(mpi_enreg%proc_distrb,&
               &ikpt,1,nband_k,isppol,me_distrb)) then

             resid(1+bdtot_index : nband_k+bdtot_index) = zero
             bdtot_index=bdtot_index+nband_k
             cycle  ! Skip the rest of the k-point loop
          end if
    
       end if ! parallel kpt
  
       ABI_ALLOCATE(pwnsfacq,(2,mkgq))
  
       call timab(984,1,tsec)
  
  !    Complete the initialization of the Hamiltonian
       call finalize_hamiltonian(gs_hamk,npw_k,istwf_k,dtset%kptns(:,ikpt))
  
       ABI_ALLOCATE(ek_k,(nband_k))
       ABI_ALLOCATE(ek_k_nd,(2,nband_k,nband_k*paw_dmft%use_dmft))
       ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
       ABI_ALLOCATE(zshift,(nband_k))
       ABI_ALLOCATE(grnl_k,(3*natom,nband_k*optforces))
  
       if (psps%usepaw==0)  then
         ABI_ALLOCATE(enl_k,(nband_k))
       end if
  
       ek_k(:)=zero
       if(paw_dmft%use_dmft==1) ek_k_nd(:,:,:)=zero
       if (optforces>0) grnl_k(:,:)=zero
       if (psps%usepaw==0) enl_k(:)=zero
       kpoint(:)=dtset%kptns(:,ikpt)
       zshift(:)=dtset%eshift
  
       if (dtset%mkmem==0) then
  !      Read (k+G) basis sphere data (same for each spin)
         nsp=dtset%nspinor
  
         call rdnpw(ikpt,isppol,nband_k,npw_k,nsp,0,dtfil%unkg)
  !      Read k+g data
         read (dtfil%unkg) kg_k(1:3,1:npw_k)
         call sphereboundary(gs_hamk%gbound,istwf_k,kg_k,dtset%mgfft,npw_k)
  !      Eventually read spherical harmonics
  
         if (psps%useylm==1) then
           read(dtfil%unylm)
           read(dtfil%unylm) ((ylm_k(muig,ilm),muig=1,npw_k),ilm=1,psps%mpsang*psps%mpsang)
         end if
  
  !      Read the wavefunction block for ikpt,isppol
         call status(counter,dtfil%filstat,iexit,level,'read wfs      ')
         tim_rwwf=1
         ABI_ALLOCATE(eig_dum,(dtset%mband))
         ABI_ALLOCATE(kg_dum,(3,0))
         ABI_ALLOCATE(occ_dum,(dtset%mband))
  
         call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,dtset%mband,mcg_disk,mpi_enreg,nband_k,nband_k,&
  &       npw_k,my_nspinor,occ_dum,-2,0,tim_rwwf,wffnow)
  
         ABI_DEALLOCATE(eig_dum)
         ABI_DEALLOCATE(kg_dum)
         ABI_DEALLOCATE(occ_dum)
  
       else
         kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
         if (mpi_enreg%paral_kgb/=1.or.istep<=1) then
           call sphereboundary(gs_hamk%gbound,istwf_k,kg_k,dtset%mgfft,npw_k)
         end if
  
  
         if (psps%useylm==1) then
           do ilm=1,psps%mpsang*psps%mpsang
             ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
  
           end do
         end if
  
       end if  ! End if for choice governed by dtset%mkmem
  
  !    Set up remaining of the Hamiltonian
  
  !    Compute (1/2) (2 Pi)**2 (k+G)**2:
       ABI_ALLOCATE(kinpw,(npw_k))
       call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,kg_k,kinpw,kpoint,npw_k)
  
       ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))
       if(dtset%nloalg(1)>0) then ! Here, precomputation of ph3d
         if (mpi_enreg%paral_kgb/=1.or.istep<=1) then
           call ph1d3d(1,natom,kg_k,gs_hamk%matblk,natom,npw_k,n1,n2,n3,gs_hamk%phkxred,ph1d,ph3d)
         end if
       end if
  
  !    Compute (k+G) vectors (only if useylm=1)
       nkpg=3*optforces*dtset%nloalg(5)
       ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
       if ((mpi_enreg%paral_kgb/=1.or.istep<=1).and.nkpg>0) then
         call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
       end if
  
  !    Compute nonlocal form factors ffnl at all (k+G):
  
       ider=0;idir=0;dimffnl=1
       ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,ntypat))
       if (mpi_enreg%paral_kgb/=1.or.istep<=1) then
         call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
  &       gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
  &       psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
  &       npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,&
  &       psps%usepaw,psps%useylm,ylm_k,ylmgr)
       end if
  
  !    Transpose the ffnl, kinpw, kpg and ph3d arrays.
       if (mpi_enreg%paral_kgb==1.and.istep<=1) then
         call prep_bandfft_tabs(dimffnl,ffnl,gs_hamk%gbound,ikpt,kinpw,&
  &       kpoint,psps%lmnmax,gs_hamk%matblk,dtset%mgfft,dtset%mkmem,mpi_enreg,nkpg,npw_k,ntypat,1,ph3d)
       end if
  
       call status(counter,dtfil%filstat,iexit,level,'call vtowfk   ')
  
  
  
       if(dtset%use_gpu_cuda==1) then
         if(mpi_enreg%paral_kgb==1) then
           ikpt_this_proc = mpi_enreg%my_kpttab(ikpt)
           ndatarecv      = bandfft_kpt(ikpt_this_proc)%ndatarecv
#if defined HAVE_GPU_CUDA
           call gpu_update_ffnl_ph3d(bandfft_kpt(ikpt_this_proc)%ph3d_gather, &
  &         bandfft_kpt(ikpt_this_proc)%ffnl_gather, &
  &         ndatarecv,dimffnl,psps%lmnmax,ntypat,natom)
#endif
         else
#if defined HAVE_GPU_CUDA
           call gpu_update_ffnl_ph3d(ph3d,ffnl,npw_k,dimffnl,psps%lmnmax,ntypat,natom)
#endif
         end if
       end if
  
       call timab(984,2,tsec)
  
       wfoptalg=dtset%wfoptalg; wfopta10=mod(wfoptalg,10)
       istwf_k=gs_hamk%istwf_k

  !    Compute the eigenvalues, wavefunction, residuals,
  !    contributions to kinetic energy, nonlocal energy, forces,
  !    and update of rhor to this k-point and this spin polarization.
       if(dtset%mkmem/=0)then

          call cgwf_orbital_shift(dtset%berryopt,cg,&
               & dtset%chkexit,cpus,dimffnl,eigen,&
               & ffnl,dtfil%filnam_ds(1),dtfil%filstat,&
               & gs_hamk,icg,ikpt,&
               & isppol,kg_k,kinpw,gs_hamk%lmnmax,gs_hamk%matblk,dtset%mband,&
               & mcg,dtset%mgfft,mkgq,mpi_enreg,gs_hamk%mpsang,&
               & gs_hamk%mpssoang,dtset%mpw,natom,nband_k,dtset%nbdblock,&
               & dtset%nkpt,dtset%nline,dtset%nloalg,&
               & npw_k,npwarr,dtset%nspinor,&
               & ntypat,gs_hamk%nvloc,dtset%ngfft(4),&
               & dtset%ngfft(5),dtset%ngfft(6),occ,dtset%ortalg,&
               & dtset%paral_kgb,ph3d,prtvol,pwind,pwind_alloc,pwnsfac,&
               & pwnsfacq,run_params,dtset%tolwfr,vhartr_new, vhartr_old,&
               & vlocal,wfoptalg)

       else if(dtset%mkmem==0)then

          write(0,*) "out-of-core orbital shifts not implemented yet"

          stop

       end if
  
       call status(counter,dtfil%filstat,iexit,level,'after vtowfk  ')
  
       call timab(985,1,tsec)
  
#if defined HAVE_GPU_CUDA
       if(dtset%use_gpu_cuda==1) then
         call gpu_finalize_ffnl_ph3d()
       end if
#endif
       ABI_DEALLOCATE(ffnl)
       ABI_DEALLOCATE(kinpw)
       ABI_DEALLOCATE(kpg_k)
       ABI_DEALLOCATE(ph3d)
       ABI_DEALLOCATE(pwnsfacq)
       call timab(985,2,tsec)
  
       ABI_DEALLOCATE(ek_k)
       ABI_DEALLOCATE(ek_k_nd)
       ABI_DEALLOCATE(grnl_k)
       ABI_DEALLOCATE(ylm_k)
       ABI_DEALLOCATE(zshift)
       if (psps%usepaw==0)  then
         ABI_DEALLOCATE(enl_k)
       end if
  
  
  !    Keep track of total number of bands (all k points so far, even for
  !    k points not treated by me)
       bdtot_index=bdtot_index+nband_k
  
  !    Also shift array memory if dtset%mkmem/=0
       if (dtset%mkmem/=0) then
         ibg=ibg+my_nspinor*nband_k
         icg=icg+npw_k*my_nspinor*nband_k
         ikg=ikg+npw_k
       end if
  
     end do ! End big k point loop
  
     call status(counter,dtfil%filstat,iexit,level,'after k loop  ')
  
     call timab(986,1,tsec)
  
   end do ! End loop over spins  
  
   call status(counter,dtfil%filstat,iexit,level,'after spinloop')
  
   if(xmpi_paral==1)then
     call timab(987,1,tsec)
     call leave_test()
     call wrtout(std_out,' vtorho: loop on k-points and spins done in parallel','COLL')
     call timab(987,2,tsec)
   end if
  
  
   call timab(988,1,tsec)
  !electric field: compute string-averaged change in Zak phase
  !along each direction, store it in dphase(idir)
  
  !ji: it is not convenient to do this anymore. Remove. Set dphase(idir)=0.0_dp.
  !eventually, dphase(idir) will have to go...
  
   if (dtset%berryopt == 4 .or. dtset%berryopt == 6 .or. dtset%berryopt == 7  .or.  &
  & dtset%berryopt ==14 .or. dtset%berryopt ==16 .or. dtset%berryopt ==17 ) dphase(:) = zero  !!HONG
  
  !In case of MPI // of a finite field calculation, send dphasek to all cpus
   if ((xmpi_paral == 1).and.(dtset%berryopt == 4 .or. dtset%berryopt == 6 .or. dtset%berryopt == 7  .or.  &
  & dtset%berryopt ==14 .or. dtset%berryopt ==16 .or. dtset%berryopt ==17)) then   !!HONG
     call xsum_mpi(dphasek,spaceComm_distrb,ierr)
   end if
  
   if (dtset%berryopt == 4 .or. dtset%berryopt == 6 .or. dtset%berryopt == 7  .or.  &
  & dtset%berryopt ==14 .or. dtset%berryopt ==16 .or. dtset%berryopt ==17 ) then  !!HONG
     ABI_DEALLOCATE(dphasek)
   end if
  
  
   call destroy_hamiltonian(gs_hamk)
  
#if defined HAVE_GPU_CUDA
   if(dtset%use_gpu_cuda==1) then
     call gpu_finalize_ham_data()
   end if
#endif
  
   ABI_DEALLOCATE(rhoaug)
   ABI_DEALLOCATE(vlocal)
   if(dtset%mkmem==0) then
     ABI_DEALLOCATE(cg_disk)
   end if
  
   ABI_ALLOCATE(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
   doccde(:)=zero !MF initialize
  
   call timab(994,1,tsec)
  
   ABI_DEALLOCATE(eknk)
   ABI_DEALLOCATE(eknk_nd)
   ABI_DEALLOCATE(kg_k)
   ABI_DEALLOCATE(grnlnk)
   if (psps%usepaw==0)  then
     ABI_DEALLOCATE(enlnk)
   end if
  
   if(psps%usepaw==1.and.(iscf>0.or.iscf==-3))  then
     ABI_DEALLOCATE(rhowfr)
     ABI_DEALLOCATE(rhowfg)
   end if
  
   call timab(994,2,tsec)
  
  
   call status(0,dtfil%filstat,iexit,level,'deallocate    ')
  
   ABI_DEALLOCATE(doccde)
   ABI_DEALLOCATE(EigMin)
    
  
  
  !Rotate labels of disk files when wf i/o is used
   if (dtset%mkmem==0) then
     wfftmp=wffnow ; wffnow=wffnew ; wffnew=wfftmp
   end if
  
   call status(0,dtfil%filstat,iexit,level,'exit          ')
  
  !Structured debugging : if prtvol=-level, stop here.
   if (prtvol==-level) then
     write(message,'(a,i0,a)')' vtorho: exit. prtvol=-',level,', debugging mode => stop '
     MSG_ERROR(message)
   end if
  
   call status(0,dtfil%filstat,iexit,level,'exit          ')
  
   call timab(980,2,tsec)
  
   DBG_EXIT("COLL")
  
 end subroutine vtorho_orbitalshifts
!!***
end module vtoorbitalshifts
