!=====================================================================
!  chen: compute on-site densities which is on a radial mesh.
!        and then stored on hc_rho1, and hc_trho1 arrays.
!        based on ABINIT's subroutine pawdenpot.F90
!
!
!=====================================================================
!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawdenpot
!! NAME
!! pawdenpot
!!
!! FUNCTION
!! Compute different (PAW) energies, densities and potentials (or potential-like quantities)
!! inside PAW spheres
!! Can also compute first-order densities potentials and second-order energies (RF calculations).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  ipert=index of perturbation (used only for RF calculation ; set ipert<=0 for GS calculations.
!!  ixc= choice of exchange-correlation scheme (see above, and below)
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms on current process, size of PAW arrays
!!  natom_tot=total number of atoms in cell
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms in unit cell.
!!  nzlmopt= if -1, compute all LM-moments of densities
!!                  initialize "lmselect" (index of non-zero LM-moments of densities)
!!           if  0, compute all LM-moments of densities
!!                  force "lmselect" to .true. (index of non-zero LM-moments of densities)
!!           if  1, compute only non-zero LM-moments of densities (stored before)
!!  option=0: compute both energies and potentials
!!         1: compute only potentials
!!         2: compute only energies
!!  paral_kgb=Flag related to the kpoint-band-fft parallelism
!!  paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  paw_an0(natom) <type(paw_an_type)>=paw arrays given on angular mesh for Ground-State
!                                      used only if ipert>0; must be set equal to paw_an for GS calc.
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawspnorb=flag: 1 if spin-orbit coupling is activated
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  ucvol=unit cell volume (bohr^3)
!!  xclevel= XC functional level
!!  znucl(ntypat)=gives the nuclear charge for all types of atoms
!!
!! OUTPUT
!!  paw_ij(natom)%dijhartree(cplex*lmn2_size)=Hartree contribution to dij;
!!                                      Enters into calculation of hartree energy
!!  ==== if option=0 or 2
!!    epaw=contribution to total energy from the PAW "on-site" part
!!    epawdc=contribution to total double-counting energy from the PAW "on-site" part
!!  ==== if option=0 or 2 and ipert<=0
!!    compch_sph=compensation charge integral inside spheres computed over spherical meshes
!!  ==== if (option=0 or 1) and paw_an(:)%has_vxc=1
!!    paw_an(natom)%vxc1[m](cplex*mesh_size,:,nspden)=XC potential calculated from "on-site" density
!!    paw_an(natom)%vxct1[m](cplex*mesh_size,:,nspden)=XC potential calculated from "on-site" pseudo density
!!    ==== if paw_an(iatom_tot)%has_vxcval==1 compute also XC potentials neglecting core charge
!!      paw_an(natom)%vxc1_val[m](cplex*mesh_size,:nspden)=XC potential calculated from spherical valence density
!!      paw_an(natom)%vxct1_val[m](cplex*mesh_size,:nspden)=XC potential calculated from spherical valence pseudo density
!!  ==== if nzlmopt==-1,
!!    paw_an(iatom_tot)%lnmselect(lm_size,nspden)=select the non-zero LM-moments of rho1 and trho1
!!  ==== if paw_an(:)%has_vhartree=1
!!    paw_an(natom)%vh1(cplex*mesh_size,1,1)=Hartree total potential calculated from "on-site" density
!!  ==== if pawspnorb>0
!!    paw_ij(natom)%dijso(cplex_dij*lmn2_size,nspden)=spin-orbit contribution to dij
!!
!! NOTES
!!  Response function calculations:
!!    In order to compute first- or second-order qunatities, paw_an (resp. paw_ij) datastructures
!!    must contain first-order quantities, namely paw_an1 (resp. paw_ij1).
!!
!! PARENTS
!!      bethe_salpeter,odamix,paw_qpscgw,respfn,scfcv,scfcv3,screening,sigma
!!
!! CHILDREN
!!      deducer0,leave_new,pawdensities,pawdijhartree,pawdijso,pawuenergy,pawxc
!!      pawxc3,pawxcm,pawxcm3,pawxcmpositron,pawxcpositron,pawxenergy,pawxpot
!!      poisson,setnoccmmp,timab,wrtout,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine c_paw_onsite_densities(mpi_enreg,&
& natom,natom_tot,nspden,ntypat,paral_kgb,paw_an,&
& paw_ij,pawang,pawrad,pawrhoij,pawtab,hc_rho1,hc_trho1)

!========= chen =========
 use c_hc_vars, only : & 
   hc_meshsz, & 
   hc_lsize

 use defs_basis
 use defs_datatypes
 use defs_abitypes

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in)    :: natom,natom_tot,nspden,ntypat,paral_kgb
 type(MPI_type),intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 type(paw_an_type),intent(in) :: paw_an(natom)
 type(paw_ij_type),intent(in) :: paw_ij(natom)
 type(pawrad_type),intent(in)    :: pawrad(ntypat)
 type(pawrhoij_type),intent(in)  :: pawrhoij(natom)
 type(pawtab_type),intent(in)    :: pawtab(ntypat)
! out vars
 real(dp),intent(out) :: hc_rho1(hc_meshsz,hc_lsize**2,natom)
 real(dp),intent(out) :: hc_trho1(hc_meshsz,hc_lsize**2,natom)

!Local variables ---------------------------------------
!scalars
 integer :: cplex,iatom,itypat,itypat0,nzlmopt
 integer :: lm_size,lmn2_size,mesh_size
 integer :: opt_compch
 character(len=500) :: message
 real(dp) :: compch_sph
!arrays
 logical,allocatable  :: lmselect_cur(:)
 real(dp),allocatable :: one_over_rad2(:),nhat1(:,:,:)
 real(dp),allocatable :: rho1(:,:,:),trho1(:,:,:)

! *************************************************************************

 write(message,'(a)') '(chen/c_paw_onsite_densities) enter '
 call wrtout(std_out,message,'COLL')

!==========================================================
!================ Big loop on atoms =======================
!==========================================================

 do iatom=1,natom

   itypat=pawrhoij(iatom)%itypat
   lmn2_size=paw_ij(iatom)%lmn2_size
   lm_size=paw_an(iatom)%lm_size
   mesh_size=pawrad(itypat)%mesh_size
   cplex=1

!  Allocations of "on-site" densities
   allocate(rho1 (cplex*mesh_size,lm_size,nspden))
   allocate(trho1(cplex*mesh_size,lm_size,nspden))
   allocate(nhat1(cplex*mesh_size,lm_size,nspden))
   allocate(lmselect_cur(lm_size));lmselect_cur(:)=.true.

!  Store some usefull quantities
   itypat0=0;if (iatom>1) itypat0=pawrhoij(iatom-1)%itypat

!  ! TO BE REMOVED
!  if (itypat/=itypat0) then
   allocate(one_over_rad2(mesh_size))
   one_over_rad2(2:mesh_size)=one/pawrad(itypat)%rad(2:mesh_size)**2

!  ===== Compute "on-site" densities (n1, ntild1, nhat1) =====

   nzlmopt=1
   opt_compch = 0 
   call pawdensities(compch_sph,cplex,iatom,lmselect_cur,paw_an(iatom)%lmselect,lm_size,&
&   nhat1,nspden,nzlmopt,one_over_rad2,opt_compch,0,-1,1,pawang,-1,&
&   pawrad(itypat),pawrhoij(iatom),pawtab(itypat),rho1,trho1)

!=================== chen: store rho1 and trho1 ===================
   if (hc_lsize**2<lm_size) then 
     print *,"pawdenpot.F90: hc_lsize**2<lm_size, our hc_lsize**2 should be big enough. ERROR! STOP"
     stop
   endif
   if (hc_meshsz<mesh_size) then 
     print *,"pawdenpot.F90: hc_meshsz<mesh_size, ERROR! STOP"
     stop
   endif
   hc_rho1 (:,:,iatom)=0.d0
   hc_trho1(:,:,iatom)=0.d0
   hc_rho1 (1:mesh_size,1:lm_size,iatom)=rho1 (1:mesh_size,1:lm_size,1)
   hc_trho1(1:mesh_size,1:lm_size,iatom)=trho1(1:mesh_size,1:lm_size,1)
   write(message,'(a,I3)')'(chen/c_paw_onsite_densities.F90) stored rho1 and trho1 for atom:',iatom
   call wrtout(std_out,message,'COLL')
!=================== end of hack ==================================
   
   deallocate(one_over_rad2)
   deallocate(lmselect_cur)
   deallocate(nhat1)
   deallocate(trho1)
   deallocate(rho1)

 enddo ! iatom

 write(message,'(a)')'(chen/c_paw_onsite_densities) leave.'
 call wrtout(std_out,message,'COLL')

end subroutine c_paw_onsite_densities
!!***
