!!****m* ABINIT/interfaces_95_drive
!! NAME
!! interfaces_95_drive
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/95_drive
!!
!! COPYRIGHT
!! Copyright (C) 2010 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!! 
!!
!! SOURCE

module interfaces_95_drive

 implicit none

interface
 subroutine afterscfloop(atindx,atindx1,cg,computed_forces,cprj,cpus,&  
  &  deltae,diffor,dtefield,dtfil,dtset,eigen,&  
  &  electronpositron,elfr,energies,etotal,&  
  &  favg,fcart,forold,fred,gresid,grewtn,grhf,grhor,&  
  &  grxc,gsqcut,hdr,indsym,irrzon,&  
  &  istep,kg,kxc,lrhor,maxfor,mgfftf,&  
  &  moved_atm_inside,mpi_enreg,&  
  &  n3xccc,nattyp,&  
  &  nfftf,ngfft,ngfftf,nhat,nkxc,npwarr,nvresid,&  
  &  occ,optres,optxc,paw_an,paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,pel,pel_cg,&  
  &  ph1d,ph1df,phnons,pion,prtfor,psps,pwind,pwind_alloc,pwnsfac,res2,resid,residm,results_gs,&  
  &  rhog,rhor,rprimd,stress_needed,strsxc,strten,symrec,synlgr,taug,taur,tollist,&  
  &  usecprj,usexcnhat,vhartr,vpsp,vxc,vxcavg,wffnow,wvl,xccc3d,xred,xred_old,ylm,ylmgr)
  use m_electronpositron
  use defs_abitypes
  use m_wffile
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(inout) :: computed_forces
  integer,intent(in) :: istep
  integer,intent(in) :: mgfftf
  integer,intent(in) :: moved_atm_inside
  integer,intent(in) :: n3xccc
  integer,intent(in) :: nfftf
  integer,intent(in) :: nkxc
  integer,intent(in) :: optres
  integer,intent(in) :: optxc
  integer,intent(in) :: prtfor
  integer,intent(in) :: pwind_alloc
  integer,intent(in) :: stress_needed
  integer,intent(in) :: usecprj
  integer,intent(in) :: usexcnhat
  real(dp),intent(in) :: cpus
  real(dp),intent(in) :: deltae
  real(dp),intent(inout) :: diffor
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(electronpositron_type),pointer :: electronpositron
  type(energies_type),intent(inout) :: energies
  real(dp),intent(inout) :: etotal
  real(dp),intent(in) :: gsqcut
  type(hdr_type),intent(inout) :: hdr
  real(dp),intent(inout) :: maxfor
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: res2
  real(dp),intent(in) :: residm
  type(results_gs_type),intent(inout) :: results_gs
  real(dp),intent(inout) :: vxcavg
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data),intent(inout) :: wvl
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: ngfftf(18)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(inout) :: cg(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  type(cprj_type),intent(in) :: cprj(dtset%natom,dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol*usecprj)
  real(dp),intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),pointer :: elfr(:,:)
  real(dp),intent(out) :: favg(3)
  real(dp),intent(out) :: fcart(3,dtset%natom)
  real(dp),intent(inout) :: forold(3,dtset%natom)
  real(dp),intent(out) :: fred(3,dtset%natom)
  real(dp),intent(out) :: gresid(3,dtset%natom)
  real(dp),intent(in) :: grewtn(3,dtset%natom)
  real(dp),intent(out) :: grhf(3,dtset%natom)
  real(dp),pointer :: grhor(:,:,:)
  real(dp),intent(out) :: grxc(3,dtset%natom)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  integer,intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  real(dp),intent(out) :: kxc(nfftf,nkxc)
  real(dp),pointer :: lrhor(:,:)
  integer,intent(in) :: nattyp(dtset%ntypat)
  real(dp),intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(inout) :: nvresid(nfftf,dtset%nspden)
  real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(paw_an_type),intent(inout) :: paw_an(dtset%natom*psps%usepaw)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%ntypat*psps%usepaw)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)
  real(dp),intent(inout) :: pel(3)
  real(dp),intent(inout) :: pel_cg(3)
  real(dp),intent(inout) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(inout) :: ph1df(2,3*(2*dtset%mgfft+1)*dtset%natom)
  real(dp),intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  real(dp),intent(inout) :: pion(3)
  integer,intent(in) :: pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(in) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(inout) :: rhog(2,nfftf)
  real(dp),intent(inout) :: rhor(nfftf,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: strsxc(6)
  real(dp),intent(out) :: strten(6)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  real(dp),intent(out) :: synlgr(3,dtset%natom)
  real(dp),pointer :: taug(:,:)
  real(dp),pointer :: taur(:,:)
  real(dp),intent(in) :: tollist(12)
  real(dp),intent(inout) :: vhartr(nfftf)
  real(dp),intent(in) :: vpsp(nfftf)
  real(dp),intent(inout) :: vxc(nfftf,dtset%nspden)
  real(dp),intent(inout) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(out) :: xred_old(3,dtset%natom)
  real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine afterscfloop
end interface

interface
 subroutine bethe_salpeter(acell,codvsn,Dtfil,Dtset,iexit,MPI_enreg,Pawang,Pawrad,Pawtab,Psps,rprim,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(inout) :: iexit
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(inout) :: Dtset
  type(mpi_type),intent(inout) :: MPI_enreg
  type(pawang_type),intent(inout) :: Pawang
  type(pseudopotential_type),intent(inout) :: Psps
  character(len=6),intent(in) :: codvsn
  type(pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Psps%usepaw)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: xred(3,Dtset%natom)
 end subroutine bethe_salpeter
end interface

interface
 subroutine brdmin(acell,atindx,atindx1,cg,cpus,dtefield,dtfil,&  
  &  dtset,ecore,eigen,electronpositron,hdr,indsym,initialized,irrzon,&  
  &  kg,mpi_enreg,mxfh,&  
  &  nattyp,nfftf,npwarr,nspinor,nxfh,occ,&  
  &  paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,pwind,pwind_alloc&  
  &  ,pwnsfac,rec_set,resid,results_gs,&  
  &  rhog,rhor,rprim,scf_history,scf_in,symrec,taug,taur,wffnew,wffnow,vel,wvl,&  
  &  xfhist,xred,xred_old,ylm,ylmgr)
  use defs_wvltypes
  use m_paw_dmft
  use defs_abitypes
  use defs_scftypes
  use defs_basis
  use defs_rectypes
  use defs_datatypes
  use m_electronpositron
  use m_wffile
  implicit none
  integer,intent(inout) :: initialized
  integer,intent(in) :: mxfh
  integer,intent(inout) :: nfftf
  integer,intent(inout) :: nspinor
  integer,intent(inout) :: nxfh
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(inout) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecore
  type(electronpositron_type),pointer :: electronpositron
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(recursion_type),intent(inout) :: rec_set
  type(results_gs_type),intent(out) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(scf_in_type),intent(inout) :: scf_in
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data),intent(inout) :: wvl
  real(dp), intent(inout) :: acell(3)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp), intent(out) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
  integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym), &
  &         2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(inout) :: pawrhoij(mpi_enreg%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(inout) :: rprim(3,3)
  integer, intent(inout) :: symrec(3,3,dtset%nsym)
  real(dp), pointer :: taug(:,:)
  real(dp), pointer :: taur(:,:)
  real(dp), intent(in) :: vel(3,dtset%natom)
  real(dp), intent(inout) :: xfhist(3,dtset%natom+4,2,mxfh)
  real(dp), intent(inout) :: xred(3,dtset%natom)
  real(dp), intent(inout) :: xred_old(3,dtset%natom)
  real(dp), intent(inout) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(inout) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine brdmin
end interface

interface
 subroutine delocint(acell,atindx,atindx1,cg,cpus,dtefield,dtfil,&  
  &  dtset,ecore,eigen,electronpositron,hdr,indsym,initialized,irrzon,&  
  &  kg,mpi_enreg,mxfh,&  
  &  nattyp,nfftf,npwarr,nspinor,nxfh,occ,&  
  &  paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,pwind&  
  &  ,pwind_alloc,pwnsfac,rec_set,resid,results_gs,&  
  &  rhog,rhor,rprim,scf_history,scf_in,symrec,taug,taur,wffnew,wffnow,vel,wvl,&  
  &  xfhist,xred,xred_old,ylm,ylmgr)
  use defs_wvltypes
  use m_paw_dmft
  use defs_abitypes
  use defs_scftypes
  use defs_basis
  use defs_rectypes
  use defs_datatypes
  use m_electronpositron
  use m_wffile
  implicit none
  integer,intent(inout) :: initialized
  integer,intent(in) :: mxfh
  integer,intent(inout) :: nfftf
  integer,intent(inout) :: nspinor
  integer,intent(inout) :: nxfh
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(inout) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecore
  type(electronpositron_type),pointer :: electronpositron
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(recursion_type),intent(inout) :: rec_set
  type(results_gs_type),intent(out) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(scf_in_type),intent(inout) :: scf_in
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data),intent(inout) :: wvl
  real(dp), intent(inout) :: acell(3)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp), intent(out) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
  integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym), &
  &         2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(inout) :: pawrhoij(mpi_enreg%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(inout) :: rprim(3,3)
  integer, intent(inout) :: symrec(3,3,dtset%nsym)
  real(dp), pointer :: taug(:,:)
  real(dp), pointer :: taur(:,:)
  real(dp), intent(in) :: vel(3,dtset%natom)
  real(dp), intent(inout) :: xfhist(3,dtset%natom+4,2,mxfh)
  real(dp), intent(inout) :: xred(3,dtset%natom)
  real(dp), intent(inout) :: xred_old(3,dtset%natom)
  real(dp), intent(inout) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(inout) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine delocint
end interface

interface
 subroutine diisRelax(acell, atindx, atindx1, cg, cpus, dtefield,&  
  &  dtfil, dtset, ecore, eigen, electronpositron, hdr, iapp, indsym, initialized,&  
  &  irrzon, kg, mpi_enreg, nattyp, nfftf, npwarr, nspinor, occ, paw_dmft, pawang,&  
  &  pawfgr, pawrad, pawrhoij, pawtab, phnons, psps, pwind, pwind_alloc, pwnsfac,&  
  &  rec_set, resid, results_gs, rhog, rhor, rprimd, scf_history, scf_in, symrec, taug, taur,&  
  &  wffnew, wffnow, wvl, xred, xred_old, ylm, ylmgr)
  use defs_wvltypes
  use m_paw_dmft
  use defs_abitypes
  use defs_scftypes
  use defs_basis
  use defs_rectypes
  use defs_datatypes
  use m_electronpositron
  use m_wffile
  implicit none
  integer,intent(in) :: iapp
  integer,intent(inout) :: initialized
  integer, intent(inout) :: nfftf
  integer,intent(inout) :: nspinor
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(inout) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecore
  type(electronpositron_type),pointer :: electronpositron
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(recursion_type),intent(inout) :: rec_set
  type(results_gs_type),intent(inout) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(scf_in_type),intent(inout) :: scf_in
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data), intent(inout) :: wvl
  real(dp), intent(inout) :: acell(3)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
  integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym), &
  &         2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(inout) :: pawrhoij(mpi_enreg%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(inout) :: rprimd(3,3)
  integer, intent(inout) :: symrec(3,3,dtset%nsym)
  real(dp), pointer :: taug(:,:)
  real(dp), pointer :: taur(:,:)
  real(dp), intent(inout) :: xred(3,dtset%natom)
  real(dp), intent(inout) :: xred_old(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine diisRelax
end interface

interface
 subroutine driver(codvsn,cpui,dtsets,filnam,filstat,&  
  &  mpi_enreg,ndtset,ndtset_alloc,npsp,pspheads,results_out)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: ndtset
  integer,intent(in) :: ndtset_alloc
  integer,intent(in) :: npsp
  character(len=6),intent(in) :: codvsn
  real(dp),intent(in) :: cpui
  character(len=fnlen),intent(in) :: filstat
  type(mpi_type),intent(inout) :: mpi_enreg
  type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
  character(len=fnlen),intent(in) :: filnam(5)
  type(pspheader_type),intent(in) :: pspheads(npsp)
  type(results_out_type),intent(inout) :: results_out(0:ndtset_alloc)
 end subroutine driver
end interface

interface
 subroutine dtfil_init1(dtfil,dtset,filnam,filstat,idtset,jdtset_,mpi_enreg,ndtset)
  use defs_basis
  use defs_abitypes
  implicit none
  integer, intent(in) :: idtset
  integer, intent(in) :: ndtset
  type(datafiles_type),intent(out) :: dtfil
  type(dataset_type),intent(in) :: dtset
  character(len=fnlen),intent(in) :: filstat
  type(mpi_type),intent(in) :: mpi_enreg
  character(len=fnlen),intent(in) :: filnam(5)
  integer :: jdtset_(0:ndtset)
 end subroutine dtfil_init1
end interface

interface
 subroutine dtfil_init2(dtfil,iapp,mpi_enreg)
  use defs_abitypes
  implicit none
  integer, intent(in) :: iapp
  type(datafiles_type),intent(inout) :: dtfil
  type(mpi_type),intent(in) :: mpi_enreg
 end subroutine dtfil_init2
end interface

interface
 subroutine elpolariz(atindx1,cg,cprj,dtefield,dtfil,dtset,etotal,enefield,gprimd,hdr,&  
  &  kg,mband,mkmem,mpi_enreg,mpw,natom,nattyp,nkpt,&  
  &  npwarr,nspinor,nsppol,ntypat,pawang,pawrad,pawrhoij,pawtab,&  
  &  pel,pel_cg,pelev,pion,psps,pwind,pwind_alloc,&  
  &  pwnsfac,rprimd,ucvol,usecprj,wffnow,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: ntypat
  integer,intent(in) :: pwind_alloc
  integer,intent(in) :: usecprj
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(inout) :: enefield
  real(dp),intent(inout) :: etotal
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  type(cprj_type),intent(in) :: cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)
  real(dp),intent(in) :: gprimd(3,3)
  integer,intent(in) :: kg(3,mpw*mkmem)
  integer,intent(in) :: nattyp(ntypat)
  integer,intent(in) :: npwarr(nkpt)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(in) :: pawrhoij(natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)
  real(dp),intent(inout) :: pel(3)
  real(dp),intent(in) :: pel_cg(3)
  real(dp),intent(inout) :: pelev(3)
  real(dp),intent(inout) :: pion(3)
  integer,intent(in) :: pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine elpolariz
end interface

interface
 subroutine gstate(acell,codvsn,cpui,dtfil,dtset,iexit,&  
  &  mpi_enreg,npwtot,nspinor,occ,pawang,pawrad,pawtab,psps,results_gs,rprim,vel,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(inout) :: iexit
  integer,intent(inout) :: nspinor
  character(len=6),intent(in) :: codvsn
  real(dp),intent(in) :: cpui
  type(datafiles_type),intent(inout) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(inout) :: pawang
  type(pseudopotential_type),intent(inout) :: psps
  type(results_gs_type),intent(inout) :: results_gs
  real(dp),intent(inout) :: acell(3)
  integer,intent(out) :: npwtot(dtset%nkpt)
  real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(inout) :: rprim(3,3)
  real(dp),intent(inout) :: vel(3,dtset%natom)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine gstate
end interface

interface
 subroutine gstateimg(acell_img,codvsn,cpui,dtfil,dtset,etotal_img,fcart_img,fred_img,iexit,&  
  &  mpi_enreg,npwtot,nspinor,occ_img,pawang,pawrad,pawtab,psps,&  
  &  rprim_img,strten_img,vel_img,xred_img)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(inout) :: iexit
  integer,intent(inout) :: nspinor
  character(len=6),intent(in) :: codvsn
  real(dp),intent(in) :: cpui
  type(datafiles_type),intent(inout) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(inout) :: pawang
  type(pseudopotential_type),intent(inout) :: psps
  real(dp),intent(inout) :: acell_img(3,dtset%nimage)
  real(dp), intent(out) :: etotal_img(dtset%nimage)
  real(dp), intent(out) :: fcart_img(3,dtset%natom,dtset%nimage)
  real(dp), intent(out) :: fred_img(3,dtset%natom,dtset%nimage)
  integer,intent(out) :: npwtot(dtset%nkpt)
  real(dp),intent(inout) :: occ_img(dtset%mband*dtset%nkpt*dtset%nsppol,dtset%nimage)
  type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(inout) :: rprim_img(3,3,dtset%nimage)
  real(dp), intent(out) :: strten_img(6,dtset%nimage)
  real(dp),intent(inout) :: vel_img(3,dtset%natom,dtset%nimage)
  real(dp),intent(inout) :: xred_img(3,dtset%natom,dtset%nimage)
 end subroutine gstateimg
end interface

interface
 subroutine gw_driver(idtset,jdtset_,ndtset,acell,codvsn,filnam,Dtfil,Dtset,iexit,MPI_enreg,&  
  &  Pawang,Pawrad,Pawtab,Psps,rprim,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: idtset
  integer,intent(inout) :: iexit
  integer,intent(in) :: ndtset
  type(datafiles_type),intent(inout) :: Dtfil
  type(dataset_type),intent(inout) :: Dtset
  type(mpi_type),intent(inout) :: MPI_enreg
  type(pawang_type),intent(inout) :: Pawang
  type(pseudopotential_type),intent(inout) :: Psps
  character(len=6),intent(in) :: codvsn
  type(pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Psps%usepaw)
  real(dp),intent(in) :: acell(3)
  character(len=fnlen),intent(in) :: filnam(5)
  integer,intent(in) :: jdtset_(0:ndtset)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: xred(3,Dtset%natom)
 end subroutine gw_driver
end interface

interface
 subroutine iofn1(filnam,filstat,mpi_enreg)
  use defs_basis
  use defs_abitypes
  implicit none
  character(len=fnlen), intent(out) :: filstat
  type(mpi_type), intent(in) :: mpi_enreg
  character(len=fnlen), intent(out) :: filnam(5)
 end subroutine iofn1
end interface

interface
 subroutine isotemp(amass,dtion,dtset,ekin,ktemp,mttk_vars,vel)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  real(dp),intent(in) :: dtion
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(out) :: ekin
  real(dp),intent(in) :: ktemp
  type(mttk_type) :: mttk_vars
  real(dp),intent(in) :: amass(dtset%natom)
  real(dp),intent(inout) :: vel(3,dtset%natom)
 end subroutine isotemp
end interface

interface
 subroutine isopress(amass,dtion,dtset,ekin,ktemp,strten,strtarget,ucvol,mttk_vars,vel,vlogv)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  real(dp),intent(in) :: dtion
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(out) :: ekin
  real(dp),intent(in) :: ktemp
  type(mttk_type) :: mttk_vars
  real(dp),intent(in) :: ucvol
  real(dp),intent(inout) :: vlogv
  real(dp),intent(in) :: amass(dtset%natom)
  real(dp),intent(inout) :: strtarget(6)
  real(dp),intent(inout) :: strten(6)
  real(dp),intent(inout) :: vel(3,dtset%natom)
 end subroutine isopress
end interface

interface
 subroutine isostress(amass,dtion,dtset,ekin,ktemp,strten,strtarget,ucvol,vel,mttk_vars)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  real(dp),intent(in) :: dtion
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(out) :: ekin
  real(dp),intent(in) :: ktemp
  type(mttk_type) :: mttk_vars
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: amass(dtset%natom)
  real(dp),intent(in) :: strtarget(6)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(inout) :: vel(3,dtset%natom)
 end subroutine isostress
end interface

interface
 subroutine loop3dte(blkflg,cg,cgindex,dtfil,dtset,d3lo,&  
  &  etotal,gmet,gprimd,gsqcut,&  !gsqcut_eff
  &  hdr,kg,kneigh,kg_neigh,kptindex,kpt3,kxc,k3xc,mband,mgfft,mkmem,mkmem_max,mk1mem,&  
  &  mpert,mpi_enreg,mpw,mvwtk,natom,nfft,nkpt,nkpt3,nkxc,nk3xc,nneigh,nspinor,nsppol,&  
  &  npwarr,occ,psps,pwind,&  
  &  rfpert,rmet,rprimd,ucvol,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mk1mem
  integer,intent(in) :: mkmem
  integer,intent(in) :: mkmem_max
  integer,intent(in) :: mpert
  integer,intent(in) :: mpw
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nk3xc
  integer,intent(in) :: nkpt
  integer,intent(in) :: nkpt3
  integer,intent(in) :: nkxc
  integer,intent(in) :: nneigh
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(inout) :: etotal
  real(dp),intent(in) :: gsqcut
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  integer,intent(out) :: blkflg(3,mpert,3,mpert,3,mpert)
  real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  integer,intent(in) :: cgindex(nkpt,nsppol)
  real(dp),intent(out) :: d3lo(2,3,mpert,3,mpert,3,mpert)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),intent(in) :: k3xc(nfft,nk3xc)
  integer,intent(in) :: kg(3,mk1mem*mpw)
  integer,intent(in) :: kg_neigh(30,nkpt,3)
  integer,intent(in) :: kneigh(30,nkpt)
  real(dp),intent(in) :: kpt3(3,nkpt3)
  integer,intent(in) :: kptindex(2,nkpt3)
  real(dp),intent(in) :: kxc(nfft,nkxc)
  real(dp),intent(in) :: mvwtk(30,nkpt)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
  integer,intent(in) :: pwind(mpw,nneigh,mkmem)
  integer,intent(in) :: rfpert(3,mpert,3,mpert,3,mpert)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine loop3dte
end interface

interface
 subroutine loper3(atindx,atindx1,blkflg,codvsn,cpus,dimcprj,doccde,&  
  &  ddkfil,dtfil,dtset,dyew,dyfrlo,dyfrnl,dyfrx1,dyfrx2,d2bbb,d2lo,d2nl,&  
  &  eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,&  
  &  etotal,fermie,gsqcut_eff,iexit,indsym,kxc,&  
  &  mkmem,mkqmem,mk1mem,mpert,mpi_enreg,mpsang,nattyp,&  
  &  nfftf,nkpt,nkxc,nspden,nspinor,nsym,occ,&  
  &  paw_an,paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,&  
  &  pertsy,prtbbb,psps,rfpert,rhog,rhor,symq,symrec,timrev,&  
  &  usecprj,vtrial,vxcavg,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer, intent(out) :: iexit
  integer, intent(in) :: mk1mem
  integer, intent(in) :: mkmem
  integer, intent(in) :: mkqmem
  integer, intent(in) :: mpert
  integer, intent(in) :: mpsang
  integer, intent(in) :: nfftf
  integer, intent(in) :: nkpt
  integer, intent(in) :: nkxc
  integer, intent(in) :: nspden
  integer, intent(inout) :: nspinor
  integer, intent(in) :: nsym
  integer, intent(in) :: prtbbb
  integer, intent(in) :: timrev
  integer, intent(in) :: usecprj
  character(len=6), intent(in) :: codvsn
  real(dp), intent(in) :: cpus
  type(datafiles_type), intent(in) :: dtfil
  type(dataset_type), intent(in) :: dtset
  real(dp), intent(out) :: etotal
  real(dp), intent(inout) :: fermie
  real(dp), intent(in) :: gsqcut_eff
  type(mpi_type), intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type), intent(inout) :: psps
  real(dp), intent(in) :: vxcavg
  integer, intent(out) :: ddkfil(3)
  integer, intent(in) :: atindx(dtset%natom)
  integer, intent(in) :: atindx1(dtset%natom)
  integer, intent(out) :: blkflg(3,mpert,3,mpert)
  real(dp), intent(out) :: d2bbb(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb)
  real(dp), intent(out) :: d2lo(2,3,mpert,3,mpert)
  real(dp), intent(out) :: d2nl(2,3,mpert,3,mpert)
  integer, intent(in) :: dimcprj(dtset%natom*psps%usepaw)
  real(dp), intent(in) :: doccde(dtset%mband*nkpt*dtset%nsppol)
  real(dp), intent(in) :: dyew(2,3,dtset%natom,3,dtset%natom)
  real(dp), intent(in) :: dyfrlo(3,3,dtset%natom)
  real(dp), intent(in) :: dyfrnl(3,3,dtset%natom)
  real(dp), intent(in) :: dyfrx1(2,3,dtset%natom,3,dtset%natom)
  real(dp), intent(in) :: dyfrx2(3,3,dtset%natom)
  real(dp), intent(in) :: eltcore(6,6)
  real(dp), intent(in) :: elteew(6+3*dtset%natom,6)
  real(dp), intent(in) :: eltfrhar(6,6)
  real(dp), intent(in) :: eltfrkin(6,6)
  real(dp), intent(in) :: eltfrloc(6+3*dtset%natom,6)
  real(dp), intent(in) :: eltfrnl(6+3*dtset%natom,6)
  real(dp), intent(in) :: eltfrxc(6+3*dtset%natom,6)
  integer, intent(in) :: indsym(4,nsym,dtset%natom)
  real(dp), intent(in) :: kxc(nfftf,nkxc)
  integer, intent(in) :: nattyp(dtset%ntypat)
  real(dp), intent(in) :: occ(dtset%mband*nkpt*dtset%nsppol)
  type(paw_an_type),intent(inout) :: paw_an(dtset%natom*psps%usepaw)
  type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  integer, intent(in) :: pertsy(3,mpert)
  integer, intent(in) :: rfpert(mpert)
  real(dp), intent(in) :: rhog(2,nfftf)
  real(dp), intent(in) :: rhor(nfftf,nspden)
  integer, intent(in) :: symq(4,2,nsym)
  integer, intent(in) :: symrec(3,3,nsym)
  real(dp), intent(inout) :: vtrial(nfftf,nspden)
  real(dp), intent(inout) :: xred(3,dtset%natom)
 end subroutine loper3
end interface

interface
 subroutine moldyn(acell,amass,atindx,atindx1,cg,cpus,&  
  &  dtefield,dtfil,dtset,ecore,eigen,electronpositron,hdr,indsym,initialized,&  
  &  irrzon,kg,mpi_enreg,mxfh,nattyp,nfftf,npwarr,nspinor,nxfh,occ,&  
  &  paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,&  
  &  phnons,psps,pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprim,&  
  &  scf_history,scf_in,symrec,taug,taur,wffnew,wffnow,vel,wvl,xfhist,xred,xred_old,ylm,ylmgr)
  use defs_wvltypes
  use m_paw_dmft
  use defs_abitypes
  use defs_scftypes
  use defs_basis
  use defs_rectypes
  use defs_datatypes
  use m_electronpositron
  use m_wffile
  implicit none
  integer,intent(inout) :: initialized
  integer,intent(in) :: mxfh
  integer,intent(inout) :: nfftf
  integer,intent(inout) :: nspinor
  integer,intent(inout) :: nxfh
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(inout) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecore
  type(electronpositron_type),pointer :: electronpositron
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type) :: paw_dmft
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(recursion_type),intent(inout) :: rec_set
  type(results_gs_type),intent(out) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(scf_in_type),intent(inout) :: scf_in
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data),intent(inout) :: wvl
  real(dp),intent(inout) :: acell(3)
  real(dp),intent(in) :: amass(dtset%natom)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp),intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp),intent(out) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
  integer,intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym), &
  &         2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer,intent(in) :: nattyp(psps%ntypat)
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type),intent(inout) :: pawrhoij(mpi_enreg%natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer,intent(in) :: pwind(pwind_alloc,2,3)
  real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp),intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),pointer :: rhog(:,:)
  real(dp),pointer :: rhor(:,:)
  real(dp),intent(inout) :: rprim(3,3)
  integer,intent(inout) :: symrec(3,3,dtset%nsym)
  real(dp),pointer :: taug(:,:)
  real(dp),pointer :: taur(:,:)
  real(dp),intent(inout) :: vel(3,dtset%natom)
  real(dp),intent(inout) :: xfhist(3,dtset%natom+4,2,mxfh)
  real(dp),intent(inout) :: xred(3,dtset%natom)
  real(dp),intent(inout) :: xred_old(3,dtset%natom)
  real(dp),intent(inout) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp),intent(inout) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine moldyn
end interface

interface
 subroutine move(acell,amass,atindx,atindx1,cg,cpus,dtefield,dtfil,dtset,&  
  &  ecore,eigen,electronpositron,hdr,indsym,initialized,irrzon,&  
  &  kg,mpi_enreg,&  
  &  nattyp,nfftf,npwarr,nspinor,occ,&  
  &  paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,&  
  &  phnons,psps,pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&  
  &  scf_history,scf_in,symrec,taug,taur,wffnew,wffnow,vel,wvl,xred,xred_old,ylm,ylmgr)
  use defs_wvltypes
  use m_paw_dmft
  use defs_abitypes
  use defs_scftypes
  use defs_basis
  use defs_rectypes
  use defs_datatypes
  use m_electronpositron
  use m_wffile
  implicit none
  integer,intent(inout) :: initialized
  integer,intent(inout) :: nfftf
  integer,intent(inout) :: nspinor
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(inout) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecore
  type(electronpositron_type),pointer :: electronpositron
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type) :: paw_dmft
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(recursion_type),intent(inout) :: rec_set
  type(results_gs_type),intent(out) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(scf_in_type),intent(inout) :: scf_in
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data),intent(inout) :: wvl
  real(dp),intent(inout) :: acell(3)
  real(dp), intent(in) :: amass(dtset%natom)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp), intent(out) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
  integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym), &
  &         2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(inout) :: pawrhoij(mpi_enreg%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(inout) :: rprimd(3,3)
  integer, intent(inout) :: symrec(3,3,dtset%nsym)
  real(dp), pointer :: taug(:,:)
  real(dp), pointer :: taur(:,:)
  real(dp), intent(inout) :: vel(3,dtset%natom)
  real(dp), intent(inout) :: xred(3,dtset%natom)
  real(dp), intent(inout) :: xred_old(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine move
end interface

interface
 subroutine nonlinear(codvsn,dtfil,dtset,etotal,iexit,&  
  &  mband,mgfft,mkmem,mpi_enreg,mpw,natom,nfft,nkpt,npwtot,nspden,&  
  &  nspinor,nsppol,nsym,occ,pawrad,pawtab,psps,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(in) :: iexit
  integer,intent(in) :: mband
  integer,intent(in) :: mgfft
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpw
  integer,intent(inout) :: natom
  integer,intent(in) :: nfft
  integer,intent(inout) :: nkpt
  integer,intent(inout) :: nspden
  integer,intent(inout) :: nspinor
  integer,intent(inout) :: nsppol
  integer,intent(inout) :: nsym
  character(len=6),intent(in) :: codvsn
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(inout) :: etotal
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(inout) :: psps
  integer,intent(out) :: npwtot(nkpt)
  real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
  type(pawrad_type),intent(inout) :: pawrad(psps%ntypat,psps%usepaw)
  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat,psps%usepaw)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine nonlinear
end interface

interface
 subroutine outscfcv(atindx1,cg,compch_fft,compch_sph,cprj,dimcprj,dtfil,dtset,&  
  &  ecut,eigen,electronpositron,elfr,etotal,fermie,gmet,gprimd,grhor,hdr,istep_mix,kg,&  
  &  lrhor,mband,mgfftc,mkmem,mpi_enreg,mpsang,mpw,natom,&  
  &  nattyp,nfft,ngfft,nhat,nkpt,npwarr,nspden,nspinor,nsppol,nsym,ntypat,n3xccc,occ,&  
  &  pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,paw_an,paw_ij,prtvol,psps,rhor,rprimd,&  
  &  taur,ucvol,usecprj,usexcnhat,wffnow,vhartr,vtrial,vxc,xccc3d,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  use m_wffile
  use m_electronpositron
  implicit none
  integer,intent(in) :: istep_mix
  integer,intent(in) :: mband
  integer,intent(in) :: mgfftc
  integer,intent(in) :: mkmem
  integer,intent(in) :: mpsang
  integer,intent(in) :: mpw
  integer,intent(in) :: n3xccc
  integer,intent(in) :: natom
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nspden
  integer,intent(in) :: nspinor
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: prtvol
  integer,intent(in) :: usecprj
  integer,intent(in) :: usexcnhat
  real(dp),intent(in) :: compch_fft
  real(dp),intent(in) :: compch_sph
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: ecut
  type(electronpositron_type),pointer :: electronpositron
  real(dp),intent(inout) :: etotal
  real(dp),intent(in) :: fermie
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: ucvol
  type(wffile_type),intent(inout) :: wffnow
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: atindx1(natom)
  real(dp),intent(inout) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
  type(cprj_type),intent(inout) :: cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)
  integer,intent(in) :: dimcprj(natom*usecprj)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
  real(dp),pointer :: elfr(:,:)
  real(dp),intent(in) :: gmet(3,3)
  real(dp),intent(in) :: gprimd(3,3)
  real(dp),pointer :: grhor(:,:,:)
  integer,intent(in) :: kg(3,mpw*mkmem)
  real(dp),pointer :: lrhor(:,:)
  integer,intent(in) :: nattyp(ntypat)
  real(dp),intent(inout) :: nhat(nfft,nspden*psps%usepaw)
  integer,intent(in) :: npwarr(nkpt)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  type(paw_an_type),intent(inout) :: paw_an(natom*psps%usepaw)
  type(paw_ij_type),intent(inout) :: paw_ij(natom)
  type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
  type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type),intent(in) :: pawrhoij(natom*psps%usepaw)
  type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
  real(dp),intent(inout) :: rhor(nfft,nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),pointer :: taur(:,:)
  real(dp),intent(in) :: vhartr(nfft)
  real(dp),intent(inout) :: vtrial(nfft,nspden)
  real(dp),intent(inout) :: vxc(nfft,nspden)
  real(dp),intent(in) :: xccc3d(n3xccc)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine outscfcv
end interface

interface
 subroutine papi_init()
  implicit none
 end subroutine papi_init
end interface

interface
 subroutine pawuj_drive(acell,atindx,atindx1,cg,cpus,dtefield,dtfil,&  
  &  dtset,ecore,eigen,electronpositron,hdr,indsym,initialized,&  
  &  irrzon,kg,mpi_enreg,nattyp,nfftf,npwarr,nspinor,occ,&  
  &  paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,&  
  &  phnons,psps,pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&  
  &  scf_history,scf_in,symrec,taug,taur,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)
  use defs_wvltypes
  use m_paw_dmft
  use defs_abitypes
  use defs_scftypes
  use defs_basis
  use defs_rectypes
  use defs_datatypes
  use m_electronpositron
  use m_wffile
  implicit none
  integer,intent(inout) :: initialized
  integer,intent(inout) :: nfftf
  integer,intent(inout) :: nspinor
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(inout) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecore
  type(electronpositron_type),pointer :: electronpositron
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type) :: paw_dmft
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(inout) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(recursion_type),intent(inout) :: rec_set
  type(results_gs_type),intent(inout) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(scf_in_type),intent(inout) :: scf_in
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data),intent(inout) :: wvl
  real(dp), intent(inout) :: acell(3)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
  integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym), &
  &         2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(inout) :: pawrhoij(mpi_enreg%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(inout) :: rprimd(3,3)
  integer, intent(inout) :: symrec(3,3,dtset%nsym)
  real(dp), pointer :: taug(:,:)
  real(dp), pointer :: taur(:,:)
  real(dp), intent(inout) :: xred(3,dtset%natom)
  real(dp), intent(inout) :: xred_old(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine pawuj_drive
end interface

interface
 subroutine respfn(codvsn,cpui,dtfil,dtset,etotal,iexit,&  
  &  mkmems,mpi_enreg,npwtot,&  
  &  nspinor,occ,pawang,pawrad,pawtab,psps,xred)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(inout) :: iexit
  integer,intent(inout) :: nspinor
  character(len=6),intent(in) :: codvsn
  real(dp),intent(in) :: cpui
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: etotal
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pawang_type),intent(inout) :: pawang
  type(pseudopotential_type),intent(inout) :: psps
  integer,intent(in) :: mkmems(3)
  integer,intent(inout) :: npwtot(dtset%nkpt)
  real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine respfn
end interface

interface
 subroutine scfcv(atindx,atindx1,cg,cpus,dtefield,dtfil,dtpawuj,&  
  &  dtset,ecore,eigen,electronpositron,fatvshift,hdr,iapp,indsym,initialized,&  
  &  irrzon,kg,mpi_enreg,nattyp,ndtpawuj,nfftf,npwarr,occ,&  
  &  paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,&  
  &  phnons,psps,pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&  
  &  scf_history,symrec,taug,taur,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)
  use defs_wvltypes
  use m_paw_dmft
  use defs_abitypes
  use m_wffile
  use defs_basis
  use defs_rectypes
  use defs_datatypes
  use m_electronpositron
  implicit none
  integer,intent(in) :: iapp
  integer,intent(inout) :: initialized
  integer,intent(in) :: ndtpawuj
  integer,intent(inout) :: nfftf
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: ecore
  type(electronpositron_type),pointer :: electronpositron
  real(dp),intent(in) :: fatvshift
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(inout) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(recursion_type),intent(inout) :: rec_set
  type(results_gs_type),intent(inout) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data),intent(inout) :: wvl
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  type(macro_uj_type),intent(inout) :: dtpawuj(0:ndtpawuj)
  real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
  integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2, &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(inout) :: pawrhoij(mpi_enreg%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(in) :: rprimd(3,3)
  integer, intent(in) :: symrec(3,3,dtset%nsym)
  real(dp), pointer :: taug(:,:)
  real(dp), pointer :: taur(:,:)
  real(dp), intent(inout) :: xred(3,dtset%natom)
  real(dp), intent(inout) :: xred_old(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine scfcv
end interface

interface
 subroutine scfcv_tmp(acell,atindx,atindx1,cg,cpus,dtefield,dtfil,dtpawuj,&  
  &  dtset,ecore,eigen,electronpositron,hdr,iapp,indsym,initialized,&  
  &  irrzon,kg,mpi_enreg,nattyp,ndtpawuj,nfftf,npwarr,nspinor,occ,&  
  &  paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,&  
  &  phnons,psps,pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&  
  &  scf_history,scf_in,symrec,taug,taur,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)
  use defs_wvltypes
  use m_paw_dmft
  use defs_abitypes
  use defs_scftypes
  use defs_basis
  use defs_rectypes
  use defs_datatypes
  use m_electronpositron
  use m_wffile
  implicit none
  integer,intent(in) :: iapp
  integer,intent(inout) :: initialized
  integer,intent(in) :: ndtpawuj
  integer,intent(inout) :: nfftf
  integer,intent(in) :: nspinor
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(inout) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecore
  type(electronpositron_type),pointer :: electronpositron
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type), intent(inout) :: paw_dmft
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type),intent(inout) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(recursion_type),intent(inout) :: rec_set
  type(results_gs_type),intent(inout) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(scf_in_type),intent(in) :: scf_in
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data),intent(inout) :: wvl
  real(dp), intent(inout) :: acell(3)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  type(macro_uj_type),intent(inout) :: dtpawuj(0:ndtpawuj)
  real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
  integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym), &
  &         2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(inout) :: pawrhoij(mpi_enreg%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(inout) :: rprimd(3,3)
  integer, intent(inout) :: symrec(3,3,dtset%nsym)
  real(dp), pointer :: taug(:,:)
  real(dp), pointer :: taur(:,:)
  real(dp), intent(inout) :: xred(3,dtset%natom)
  real(dp), intent(inout) :: xred_old(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine scfcv_tmp
end interface

interface
 subroutine scphon(acell,amass, atindx, atindx1, cg, cpus, dtefield,&  
  &  dtfil, dtset, ecore, eigen, electronpositron, hdr, indsym, initialized,&  
  &  irrzon, kg, mpi_enreg, nattyp, nfftf, npwarr, nspinor, occ, pawang,&  
  &  paw_dmft,pawfgr, pawrad, pawrhoij, pawtab, phnons, psps, pwind, pwind_alloc, pwnsfac,&  
  &  rec_set, resid, results_gs, rhog, rhor,rprimd, scf_history, scf_in,symrec, taug,taur,&  
  &  wffnew, wffnow, wvl, xred, xred_old, ylm, ylmgr)
  use defs_wvltypes
  use m_paw_dmft
  use defs_abitypes
  use defs_scftypes
  use defs_basis
  use defs_rectypes
  use defs_datatypes
  use m_electronpositron
  use m_wffile
  implicit none
  integer,intent(inout) :: initialized
  integer,intent(inout) :: nfftf
  integer,intent(inout) :: nspinor
  integer,intent(in) :: pwind_alloc
  real(dp),intent(in) :: cpus
  type(efield_type),intent(inout) :: dtefield
  type(datafiles_type),intent(inout) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ecore
  type(electronpositron_type),pointer :: electronpositron
  type(hdr_type),intent(inout) :: hdr
  type(mpi_type),intent(inout) :: mpi_enreg
  type(paw_dmft_type) :: paw_dmft
  type(pawang_type),intent(in) :: pawang
  type(pawfgr_type) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(recursion_type),intent(inout) :: rec_set
  type(results_gs_type),intent(inout) :: results_gs
  type(scf_history_type),intent(inout) :: scf_history
  type(scf_in_type),intent(inout) :: scf_in
  type(wffile_type),intent(inout) :: wffnew
  type(wffile_type),intent(inout) :: wffnow
  type(wvl_data),intent(inout) :: wvl
  real(dp), intent(inout) :: acell(3)
  real(dp), intent(in) :: amass(dtset%natom)
  integer,intent(in) :: atindx(dtset%natom)
  integer,intent(in) :: atindx1(dtset%natom)
  real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
  real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
  integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
  integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym), &
  &         2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
  integer, intent(in) :: nattyp(psps%ntypat)
  integer, intent(in) :: npwarr(dtset%nkpt)
  real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawrhoij_type), intent(inout) :: pawrhoij(mpi_enreg%natom*psps%usepaw)
  type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym), &
  &         (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
  integer, intent(in) :: pwind(pwind_alloc,2,3)
  real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
  real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp), pointer :: rhog(:,:)
  real(dp), pointer :: rhor(:,:)
  real(dp), intent(inout) :: rprimd(3,3)
  integer, intent(inout) :: symrec(3,3,dtset%nsym)
  real(dp), pointer :: taug(:,:)
  real(dp), pointer :: taur(:,:)
  real(dp), intent(inout) :: xred(3,dtset%natom)
  real(dp), intent(inout) :: xred_old(3,dtset%natom)
  real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
  real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 end subroutine scphon
end interface

interface
 subroutine scphon_phonon_init (dtfil,natom_primitive_cell,&  
  &  nphononq,phonon_eigvec_ref,phonon_eigval_ref)
  use defs_basis
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom_primitive_cell
  integer,intent(in) :: nphononq
  type(datafiles_type),intent(in) :: dtfil
  real(dp),intent(out) :: phonon_eigval_ref(3*natom_primitive_cell,nphononq)
  real(dp),intent(out) :: phonon_eigvec_ref(2,3*natom_primitive_cell,3*natom_primitive_cell,nphononq)
 end subroutine scphon_phonon_init
end interface

interface
 subroutine scphon_qpoint_init (nphononq,phononq,supercell_multiplicity)
  use defs_basis
  implicit none
  integer,intent(in) :: nphononq
  integer,intent(in) :: supercell_multiplicity(3)
  real(dp),intent(out) :: phononq(3,nphononq)
 end subroutine scphon_qpoint_init
end interface

interface
 subroutine scphon_ft_fcart(sqrt_amass_pcell,fcart,natom,natom_primitive_cell,nphononq,phononq,&  
  &  forces_on_atoms_ft,pcell_atom_in_supercell,supercell_vectors)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: natom_primitive_cell
  integer,intent(in) :: nphononq
  real(dp),intent(in) :: fcart(3,natom)
  real(dp),intent(out) :: forces_on_atoms_ft(2,3*natom_primitive_cell,nphononq)
  integer,intent(in) :: pcell_atom_in_supercell(natom)
  real(dp),intent(in) :: phononq(3,nphononq)
  real(dp),intent(in) :: sqrt_amass_pcell(natom_primitive_cell)
  real(dp),intent(in) :: supercell_vectors(3,natom)
 end subroutine scphon_ft_fcart
end interface

interface
 subroutine scphon_new_frequencies(forces_on_atoms_ft,istep,natom_primitive_cell,&  
  &  normal_mode_displacements,nphononq,nsym_primitive_cell,pcell,phonon_eigvec_ref,&  
  &  phonon_eigval2_averaged,phonon_eigval,phononq,qsym_map)
  use defs_basis
  use m_primcell_ddb_info
  implicit none
  integer,intent(in) :: istep
  integer,intent(in) :: natom_primitive_cell
  integer,intent(in) :: nphononq
  integer,intent(in) :: nsym_primitive_cell
  type(primcell_ddb_info),intent(inout) :: pcell
  real(dp),intent(in) :: forces_on_atoms_ft(2,3*natom_primitive_cell,nphononq)
  real(dp),intent(in) :: normal_mode_displacements(3*natom_primitive_cell,nphononq)
  real(dp),intent(out) :: phonon_eigval(3*natom_primitive_cell,nphononq)
  real(dp),intent(inout) :: phonon_eigval2_averaged(3*natom_primitive_cell,nphononq)
  real(dp),intent(in) :: phonon_eigvec_ref(2,3*natom_primitive_cell,3*natom_primitive_cell,nphononq)
  real(dp),intent(in) :: phononq(3,nphononq)
  integer,intent(in) :: qsym_map(nphononq,nsym_primitive_cell,2)
 end subroutine scphon_new_frequencies
end interface

interface
 subroutine scphon_freq_to_normmode (minusq_map,natom_primitive_cell,normal_mode_displacements,&  
  &  nphononq,phonon_eigval,scphon_temp)
  use defs_basis
  implicit none
  integer,intent(in) :: natom_primitive_cell
  integer,intent(in) :: nphononq
  real(dp),intent(in) :: scphon_temp
  integer,intent(in) :: minusq_map(nphononq)
  real(dp),intent(out) :: normal_mode_displacements(3*natom_primitive_cell,nphononq)
  real(dp),intent(in) :: phonon_eigval(3*natom_primitive_cell,nphononq)
 end subroutine scphon_freq_to_normmode
end interface

interface
 subroutine scphon_update_xcart (sqrt_amass_pcell,cartesian_displacements,natom,&  
  &  natom_primitive_cell,normal_mode_displacements,&  
  &  nphononq,pcell_atom_in_supercell,phonon_eigvec_ref,phononq,supercell_vectors,xcart,xcart0)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: natom_primitive_cell
  integer,intent(in) :: nphononq
  real(dp),intent(inout) :: cartesian_displacements(3,natom)
  real(dp),intent(in) :: normal_mode_displacements(3*natom_primitive_cell,nphononq)
  integer,intent(in) :: pcell_atom_in_supercell(natom)
  real(dp),intent(in) :: phonon_eigvec_ref(2,3*natom_primitive_cell,3*natom_primitive_cell,nphononq)
  real(dp),intent(in) :: phononq(3,nphononq)
  real(dp),intent(in) :: sqrt_amass_pcell(natom_primitive_cell)
  real(dp),intent(in) :: supercell_vectors(3,natom)
  real(dp),intent(inout) :: xcart(3,natom)
  real(dp),intent(in) :: xcart0(3,natom)
 end subroutine scphon_update_xcart
end interface

interface
 subroutine scphon_build_qsym_map(nphononq,nsym_primitive_cell,phononq,&  
  &  qsym_map,symrec_primitive_cell)
  use defs_basis
  implicit none
  integer,intent(in) :: nphononq
  integer,intent(in) :: nsym_primitive_cell
  real(dp),intent(in) :: phononq(3,nphononq)
  integer,intent(out) :: qsym_map(nphononq,nsym_primitive_cell,2)
  integer,intent(in) :: symrec_primitive_cell(3,3,nsym_primitive_cell)
 end subroutine scphon_build_qsym_map
end interface

interface
 subroutine scphon_supercell_vectors_init(natom,natom_primitive_cell,&  
  &  pcell,pcell_atom_in_supercell,&  
  &  supercell_multiplicity,supercell_vectors,xred)
  use defs_basis
  use m_primcell_ddb_info
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: natom_primitive_cell
  type(primcell_ddb_info),intent(inout) :: pcell
  integer,intent(in) :: supercell_multiplicity(3)
  integer,intent(out) :: pcell_atom_in_supercell(natom)
  real(dp),intent(out) :: supercell_vectors(3,natom)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine scphon_supercell_vectors_init
end interface

interface
 subroutine scphon_check_fcart(cartesian_displacements,fcart,natom)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  real(dp),intent(in) :: cartesian_displacements(3,natom)
  real(dp),intent(in) :: fcart(3,natom)
 end subroutine scphon_check_fcart
end interface

interface
 subroutine scphon_free_energy(free_energy,istep,t_phonon_dos,scphon_temp)
  use defs_basis
  use m_phdos
  implicit none
  integer,intent(in) :: istep
  real(dp),intent(out) :: free_energy
  real(dp),intent(in) :: scphon_temp
  type(phonon_dos_type),intent(inout) :: t_phonon_dos
 end subroutine scphon_free_energy
end interface

interface
 subroutine scphon_make_phonon_dos (dos_smearing,natom_primitive_cell,&  
  &  nfreq_int,nphononq,maxfreq,minfreq,phonon_dos,phonon_eigval)
  use defs_basis
  implicit none
  integer,intent(in) :: natom_primitive_cell
  integer,intent(in) :: nfreq_int
  integer,intent(in) :: nphononq
  real(dp),intent(in) :: dos_smearing
  real(dp),intent(out) :: maxfreq
  real(dp),intent(out) :: minfreq
  real(dp),intent(out) :: phonon_dos(nfreq_int)
  real(dp),intent(in) :: phonon_eigval(3*natom_primitive_cell,nphononq)
 end subroutine scphon_make_phonon_dos
end interface

interface
 subroutine scphon_interpolate_phonon_and_dos (natom_primitive_cell,&  
  &  nphononq,pcell,t_phonon_dos,phonon_eigval,phonon_eigvec_ref,&  
  &  phononq,supercell_multiplicity)
  use defs_basis
  use m_primcell_ddb_info
  use m_phdos
  implicit none
  integer,intent(in) :: natom_primitive_cell
  integer,intent(in) :: nphononq
  type(primcell_ddb_info),intent(inout) :: pcell
  type(phonon_dos_type),intent(inout) :: t_phonon_dos
  integer,intent(in) :: supercell_multiplicity(3)
  real(dp),intent(in) :: phonon_eigval(3*natom_primitive_cell,nphononq)
  real(dp),intent(in) :: phonon_eigvec_ref(2,3*natom_primitive_cell,3*natom_primitive_cell,nphononq)
  real(dp),intent(in) :: phononq(3,nphononq)
 end subroutine scphon_interpolate_phonon_and_dos
end interface

interface
 subroutine scphon_freq_to_dynmat(dynmat,natom_primitive_cell,&  
  &  nphononq,phonon_eigval,phonon_eigvec_ref)
  use defs_basis
  implicit none
  integer,intent(in) :: natom_primitive_cell
  integer,intent(in) :: nphononq
  real(dp),intent(out) :: dynmat(2,3,natom_primitive_cell,3,natom_primitive_cell,nphononq)
  real(dp),intent(in) :: phonon_eigval(3*natom_primitive_cell,nphononq)
  real(dp),intent(in) :: phonon_eigvec_ref(2,3*natom_primitive_cell,3*natom_primitive_cell,nphononq)
 end subroutine scphon_freq_to_dynmat
end interface

interface
 subroutine scphon_dynmat_to_freq2(dynmat,natom_primitive_cell,&  
  &  nphononq,phonon_eigval,phonon_eigvec_ref)
  use defs_basis
  implicit none
  integer, intent(in) :: natom_primitive_cell
  integer, intent(in) :: nphononq
  real(dp),intent(in) :: dynmat(2,3,natom_primitive_cell,3,natom_primitive_cell,nphononq)
  real(dp),intent(out) :: phonon_eigval(3*natom_primitive_cell,nphononq)
  real(dp),intent(in) :: phonon_eigvec_ref(2,3*natom_primitive_cell,3*natom_primitive_cell,nphononq)
 end subroutine scphon_dynmat_to_freq2
end interface

interface
 subroutine print_phonfreq(istep,natom_primitive_cell,nphononq,phonon_eigval)
  use defs_basis
  implicit none
  integer,intent(in) :: istep
  integer,intent(in) :: natom_primitive_cell
  integer,intent(in) :: nphononq
  real(dp),intent(in) :: phonon_eigval(3*natom_primitive_cell,nphononq)
 end subroutine print_phonfreq
end interface

interface
 subroutine screening(acell,codvsn,Dtfil,Dtset,iexit,MPI_enreg,Pawang,Pawrad,Pawtab,Psps,rprim)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(inout) :: iexit
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(inout) :: Dtset
  type(mpi_type),intent(inout) :: MPI_enreg
  type(pawang_type),intent(inout) :: Pawang
  type(pseudopotential_type),intent(inout) :: Psps
  character(len=6),intent(in) :: codvsn
  type(pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Dtset%usepaw)
  type(pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Dtset%usepaw)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
 end subroutine screening
end interface

interface
 subroutine sigma(acell,codvsn,Dtfil,Dtset,iexit,MPI_enreg,Pawang,Pawrad,Pawtab,Psps,rprim,xred,converged)
  use defs_basis
  use defs_abitypes
  use defs_datatypes
  implicit none
  integer,intent(inout) :: iexit
  type(datafiles_type),intent(in) :: Dtfil
  type(dataset_type),intent(inout) :: Dtset
  type(mpi_type),intent(inout) :: MPI_enreg
  type(pawang_type),intent(inout) :: Pawang
  type(pseudopotential_type),intent(inout) :: Psps
  character(len=6),intent(in) :: codvsn
  logical,intent(out) :: converged
  type(pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Psps%usepaw)
  type(pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Psps%usepaw)
  real(dp),intent(in) :: acell(3)
  real(dp),intent(in) :: rprim(3,3)
  real(dp),intent(in) :: xred(3,Dtset%natom)
 end subroutine sigma
end interface

interface
 subroutine testfi(etotal,filnam,filstat,fred,natom,strten,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  real(dp),intent(in) :: etotal
  character(len=fnlen),intent(in) :: filstat
  character(len=fnlen),intent(in) :: filnam(5)
  real(dp),intent(in) :: fred(3,natom)
  real(dp),intent(in) :: strten(6)
  real(dp),intent(in) :: xred(3,natom)
 end subroutine testfi
end interface

interface
 subroutine timana(mpi_enreg,natom,nband,ndtset,nfft,nkpt,npwtot,nsppol,timopt, papiopt)
  use defs_abitypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ndtset
  integer,intent(in) :: nfft
  integer,intent(in) :: nkpt
  integer,intent(in) :: nsppol
  integer,intent(in) :: papiopt
  integer,intent(in) :: timopt
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: nband(nkpt*nsppol)
  integer,intent(in) :: npwtot(nkpt)
 end subroutine timana
end interface

end module interfaces_95_drive
!!***
