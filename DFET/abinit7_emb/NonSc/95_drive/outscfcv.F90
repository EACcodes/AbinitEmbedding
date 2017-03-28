! {\src2tex{textfont=tt}}
!!****f* ABINIT/outscfcv
!! NAME
!! outscfcv
!!
!! FUNCTION
!! Output routine for the scfcv.F90 routine
!!
!! COPYRIGHT
!! Copyright (C) 2005-2010 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions (see also side effects)
!!  compch_fft=compensation charge, from FFT grid
!!  compch_sph=compensation charge, from sphere
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecrpj)=<p_lmn|Cnk> coefficients for each WF |Cnk>
!!   and each |p_lmn> non-local projector. See also side effects
!!  dimcprj(natom*usecprj)=array of dimensions of array cprj (not ordered)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  elfr(nfft,nspden(+1))=electron localization function, real space.
!!   (+1) if spin-polarized in order to get total, spin up and spin down elf
!!  etotal=total energy
!!  fermie= Fermi energy
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  grhor(nfft,nspden,3)= gradient of electron density in electrons/bohr**4, real space
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  istep_mix=index of the number of steps for the SCF mixing (can be <istep)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  ------Removed in beautification because unused MS ------
!!  kssform=govern the Kohn-Sham Structure file format
!!  --------------------------------------------------------
!!  lrhor(nfft,nspden)= Laplacian of electron density in electrons/bohr**5, real space
!!  mband=maximum number of bands
!!  mgfftc=maximum size of 1D FFTs for the PAW coarse grid
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw.
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nfft=(effective) number of FFT grid points (for this processor) (see NOTES at beginning of scfcv)
!!  ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!  nhat(nfft,nspden*usepaw)= compensation charge density  (PAW)
!!  nkpt=number of k points.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetries in space group
!!  ntypat=number of types of atoms in unit cell.
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  occ(mband*nkpt*nsppol)=occupation number for each band (usually 2) for each k.
!!  paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr(natom) <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!            note:structure factors are given on the coarse grid for PAW
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  ------Removed in beautification because unused MS ------
!!  rhog(nfft,nspden)=total electron density in electrons/bohr**3, reciprocal space.
!!  --------------------------------------------------------
!!  rhor(nfft,nspden)=total electron density in electrons/bohr**3, real space.
!!  rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!  taur(nfft,nspden)=total kinetic energy density in bohr**(-5), real space.
!!  ucvol=unit cell volume (bohr**3)
!!  usecprj=1 if cprj datastructure has been allocated
!!  usexcnhat= flag controling use of compensation density in the computation of Vxc
!!  vhartr(nfft)=Hartree potential
!!  vxc(nfft,nspden)=xc potential
!!  ------Removed in beautification because unused MS ------
!!  vxcavg=vxc average
!!  --------------------------------------------------------
!!  wffnow=information about wf disk file
!!  vtrial(nfft,nspden)=the trial potential
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  If prtwant==3 the following quantitities are updated using the unitary transformation
!!  defining the QP amplitudes in terms of the KS basis set:
!!   cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
!!   cprj(natom,nspinor*mband*mkmem*nsppol*usecrpj)=<p_lmn|Cnk> coefficients for each WF |Cnk>
!!   and each |p_lmn> non-local projector
!!
!! NOTES
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      abi_etsf_electrons_put,abi_etsf_geo_put,bonds_lgth_angles,calc_cs
!!      calc_efg,calc_fc,calcdensph,denfgr,dos_degeneratewfs,ioarr,leave_new
!!      mati3inv,mlwfovlp,mlwfovlp_qp,optics_paw,optics_paw_core,out1dm,outkss
!!      outwant,partial_dos_fractions,partial_dos_fractions_paw,pawmkaewf
!!      pawprt,poslifetime,printbxsf,prt_cml2,prtfatbands,read_atomden
!!      tetrahedron,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine outscfcv(atindx1,cg,compch_fft,compch_sph,cprj,dimcprj,dtfil,dtset,&
& ecut,eigen,electronpositron,elfr,etotal,fermie,gmet,gprimd,grhor,hdr,istep_mix,kg,&
& lrhor,mband,mgfftc,mkmem,mpi_enreg,mpsang,mpw,natom,&
& nattyp,nfft,ngfft,nhat,nkpt,npwarr,nspden,nspinor,nsppol,nsym,ntypat,n3xccc,occ,&
& pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,paw_an,paw_ij,prtvol,psps,rhor,rprimd,&
& taur,ucvol,usecprj,usexcnhat,wffnow,vhartr,vtrial,vxc,xccc3d,xred)

!===================== chen ===============
 use c_hc_vars
 use interfaces_12_hide_mpi
!===================== chen ===============

 use defs_basis
 use defs_datatypes
 use m_wffile 
 use defs_abitypes
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype
 use defs_parameters

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_47_xml
 use interfaces_61_ionetcdf
 use interfaces_62_iowfdenpot
 use interfaces_62_occeig
 use interfaces_66_paw
 use interfaces_67_common
 use interfaces_68_gw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep_mix,mband,mgfftc,mkmem,mpsang,mpw,n3xccc,natom,nfft !,kssform
 integer,intent(in) :: nkpt,nspden,nsppol,nsym,ntypat,prtvol,usecprj,usexcnhat
 integer,intent(in) :: nspinor
 real(dp),intent(in) :: compch_fft,compch_sph,ecut,fermie,ucvol
 real(dp),intent(inout) :: etotal!,vxcavg
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow
!arrays
 integer,intent(in) :: atindx1(natom),dimcprj(natom*usecprj)
 integer,intent(in) :: kg(3,mpw*mkmem),nattyp(ntypat),ngfft(18),npwarr(nkpt)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol)
 real(dp),intent(in) :: rprimd(3,3),vhartr(nfft),xccc3d(n3xccc)
 real(dp),intent(inout) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(inout) :: nhat(nfft,nspden*psps%usepaw)!,rhog(nfft,nspden)
 real(dp),intent(inout) :: rhor(nfft,nspden),vtrial(nfft,nspden)
 real(dp),intent(inout) :: vxc(nfft,nspden),xred(3,natom)
 real(dp),pointer :: elfr(:,:),grhor(:,:,:),lrhor(:,:),taur(:,:)
 type(cprj_type),intent(inout) :: cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)
 type(paw_an_type),intent(inout) :: paw_an(natom*psps%usepaw)
 type(paw_ij_type),intent(inout) :: paw_ij(natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: accessfil,coordn,fatbands_flag,fformr,fformv
 integer :: ierr,ifft,ii,ikpt,ispden,isppol,isym
 integer :: m_dos_flag,mbesslang,ndosfraction,nzlmopt
 integer :: occopt,partial_dos_flag,paw_dos_flag,pawfatbnd,prt1dm
 integer :: prtcml,prtcs,prtden,prtdos,prtefg,prtelf,prtfc,prtgden,prtgeo,prtkden,prtlden,prtnabla
 integer :: prtpot,prtstm,prtvha,prtvhxc,prtvxc,rdwr,rdwrpaw,timrev
 logical :: use_afm
 real(dp) :: norm
 character(len=500) :: message
!arrays
 integer,allocatable :: symrec(:,:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: dos_fractions(:,:,:,:),dos_fractions_m(:,:,:,:),dos_fractions_average_m(:,:,:,:)
 real(dp),allocatable :: dos_fractions_paw1(:,:,:,:)
 real(dp),allocatable :: dos_fractions_pawt1(:,:,:,:),eigen2(:)
 real(dp),allocatable :: eigen2bxsf(:,:,:),elfr_down(:,:),elfr_up(:,:)
 real(dp),allocatable :: rhor_paw(:,:),rhor_paw_core(:,:),rhor_paw_val(:,:),vwork(:,:)
 type(pawrhoij_type),allocatable :: pawrhoij_dum(:)

! *************************************************************************

!DEBUG
!write(6,*)' outscfcv : enter '
!ENDDEBUG

 if (usecprj==0.and.psps%usepaw==1.and. &
& (dtset%prtwant==2.or.dtset%prtwant==3.or.dtset%prtnabla>0.or.dtset%prtdos==3 &
& .or.dtset%kssform==3.or.dtset%pawfatbnd>0.or.dtset%pawprtwf>0)) then
   write (message,'(8a)')ch10,&
&   ' outscfcv : ERROR- ',ch10,&
&   ' cprj datastructure must be allocated',ch10,&
&   ' with options prtwant=2,3, prtnabla>0, prtdos>3, kssform==3, pawfatbnd>0, pawprtwf>0',ch10,&
&   ' Action: change pawusecp input keyword.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

!wannier interface
 if (dtset%prtwant==2) then

   call mlwfovlp(atindx1,cg,cprj,dtset,dtfil,eigen,gprimd,hdr,kg,&
&   mband,mgfftc,mkmem,mpi_enreg,mpw,natom,&
&   nattyp,nfft,ngfft,nkpt,npwarr,nspinor,nsppol,ntypat,&
&   pawang,pawrad,pawtab,prtvol,psps,rprimd,ucvol,xred)

 else if (dtset%prtwant==3) then

!  Convert cg and eigen to GW quasiparticle wave functions and eigenvalues in mlwfovlp_qp
   allocate(eigen2(mband*nkpt*nsppol))
   eigen2=eigen

   call mlwfovlp_qp(cg,cprj,dtset,dtfil,eigen2,mband,mkmem,mpw,natom,&
&   nkpt,npwarr,nspden,nspinor,nsppol,ntypat,Hdr,pawtab,rprimd,MPI_enreg)

!  Call Wannier90
   call mlwfovlp(atindx1,cg,cprj,dtset,dtfil,eigen2,gprimd,hdr,kg,&
&   mband,mgfftc,mkmem,mpi_enreg,mpw,natom,&
&   nattyp,nfft,ngfft,nkpt,npwarr,nspinor,nsppol,ntypat,&
&   pawang,pawrad,pawtab,prtvol,psps,rprimd,ucvol,xred)

!  this is the old implementation, risky due to unpredictable size effects
!  now eigen is not overwritten, one should use other ways to print the GW corrections
!  eigen=eigen2
   deallocate(eigen2)
 end if !prtwant

!
!if accesswff == 2 then set all outputs to netcdf format
!if accesswff == 3 then set all outputs to ETSF format
!
 accessfil = 0
 if (dtset%accesswff == 2) accessfil = 1
 if (dtset%accesswff == 3) accessfil = 3
 if (dtset%accesswff == 1) accessfil = 4

 occopt=dtset%occopt;

 pawfatbnd=dtset%pawfatbnd
 prtden=dtset%prtden ; prtpot=dtset%prtpot ; prtgeo=dtset%prtgeo
 prtcml=dtset%prtcml ; prtdos=dtset%prtdos ; prtstm=dtset%prtstm
 prt1dm=dtset%prt1dm ; prtvha=dtset%prtvha ; prtvhxc=dtset%prtvhxc
 prtvxc=dtset%prtvxc ; prtnabla=dtset%prtnabla; prtefg=dtset%prtefg
 prtcs=dtset%prtcs   ; prtfc=dtset%prtfc ; prtkden=dtset%prtkden
 prtelf=dtset%prtelf ; prtgden=dtset%prtgden; prtlden=dtset%prtlden

!Warnings :
!- core charge is excluded from the charge density;
!- the potential is the INPUT vtrial.
 if(  mpi_enreg%paral_compil_kpt==0                         .or. &
& (mpi_enreg%me==0 .and. mpi_enreg%paral_compil_fft==0 ) .or. &
& (mpi_enreg%paral_compil_fft==1 .and. mpi_enreg%me_band==0 .and. mpi_enreg%me_kpt==0)) then

!  We output the density.
   if (prtden/=0) then
     rdwr=2 ; fformr=52 ; rdwrpaw=0
!    We create the file name.
     call ioarr(accessfil, rhor, dtset, etotal, fformr, dtfil%fnameabo_app_den, hdr, mpi_enreg, &
&     nfft, pawrhoij_dum, rdwr, rdwrpaw)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_den, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_den)
     end if
   end if

 end if ! if master

!============================================================================
!============================ chen: obtain rhor_paw, and exit ===============
!============================================================================
   if (psps%usepaw==1) then
     allocate(rhor_paw(pawfgr%nfft,nspden))
     if (bNCPP_PHI>0) then 
       call c_smooth_denfgr(mpi_enreg,natom,nspden,nhat,ntypat,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,&
&                  rhor,rhor_paw,prtvol,psps,dtset%typat)
       call c_wrtlog("(outscfcv) === using c_smooth_denfgr, replace AE with NCPP ===")
     else
       call denfgr(mpi_enreg,natom,nspden,nhat,ntypat,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,&
&                rhor,rhor_paw,prtvol,psps,dtset%typat)
       call c_wrtlog("(outscfcv) === using original denfgr ===")
     endif
     !
     ! rhor_paw is only MPI-safe for the master node (due to denfgr())
     ! we need to broadcast master's rhor_paw to all the other nodes
     !
     hc_rhor_paw(:)=rhor_paw(:,1)
     if (mpi_enreg%nproc>1) then
       call xbarrier_mpi(mpi_enreg%spaceComm,ierr)
       call xcast_mpi(hc_rhor_paw,0,mpi_enreg%spaceComm,ierr)
       call xbarrier_mpi(mpi_enreg%spaceComm,ierr)
     endif
     norm = SUM(hc_rhor_paw(:))*ucvol/PRODUCT(pawfgr%ngfft(1:3))
     write(message,'(a,g22.14)') '(chen/outscfcv) hc_rhor_paw   - NORM OF DENSITY: ',norm; call wrtout(6,message,'COLL')
     write(message,'(a,g22.14)') '                rhor-hat      - NORM OF DENSITY: ',sum(rhor(:,1)-nhat(:,1))*ucvol/PRODUCT(pawfgr%ngfft(1:3)); call wrtout(6,message,'COLL')
     write(message,'(a)')'(chen/outscfcv) done hc_rhor_paw. leave outscfcv.F90.';  call wrtout(std_out,message,'COLL')

     deallocate(rhor_paw)
     write(message,*)'chen: return from outscfcv.F90()'
     return 
   endif
!============================ end of hack ===================================   

!! MS - Printing of PAWDEN parallellised and several possible options
!!      included
!We output the total electron density in the PAW case
!this requires removing nhat from rhor and making PAW on-site corrections
 if (dtset%pawprtden>0 .and. psps%usepaw==1) then
!  pawprtden 1 --> output PAW valence density
!  "     2 --> output PAW valence+core density
!  "     3 --> output core, valence and full atomic protodensity
!  "     4 --> options 1+3
!  "     5 --> options 2+3
   if (dtset%pawprtden/=3) then ! calc PAW valence density
     allocate(rhor_paw(pawfgr%nfft,nspden))
     call denfgr(mpi_enreg,natom,nspden,nhat,ntypat,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,&
&     rhor,rhor_paw,prtvol,psps,dtset%typat)
!    Check normalisation
     if (prtvol>9) then
       norm = SUM(rhor_paw(:,1))*ucvol/PRODUCT(pawfgr%ngfft(1:3))
       write(message,'(a,F8.4)') '  PAWDEN - NORM OF DENSITY: ',norm
       call wrtout(6,message,'COLL')
     end if
   end if
   if (dtset%pawprtden>1) then ! We will need the core density
     allocate(rhor_paw_core(pawfgr%nfft,nspden))
     call read_atomden(mpi_enreg,natom,nspden,ntypat,pawfgr,rhor_paw_core,&
&     dtset%typat,rprimd,xred,prtvol,file_prefix='core   ')
!    Check normalisation
     if (prtvol>9) then
       norm = SUM(rhor_paw_core(:,1))*ucvol/PRODUCT(pawfgr%ngfft(1:3))
       write(message,'(a,F8.4)') '  ATMDEN - NORM OF CORE DENSITY: ', norm
       call wrtout(6,message,'COLL')
     end if
   end if
   if (dtset%pawprtden>2) then ! We will need the valence protodensity
     allocate(rhor_paw_val(pawfgr%nfft,nspden))
     call read_atomden(mpi_enreg,natom,nspden,ntypat,pawfgr,rhor_paw_val,&
&     dtset%typat,rprimd,xred,prtvol,file_prefix='valence')
!    Check normalisation
     if (prtvol>9) then
       norm = SUM(rhor_paw_val(:,1))*ucvol/PRODUCT(pawfgr%ngfft(1:3))
       write(message,'(a,F8.4)') '  ATMDEN - NORM OF VALENCE PROTODENSITY: ', norm
       call wrtout(6,message,'COLL')
     end if
   end if
   if(  mpi_enreg%paral_compil_kpt==0                         .or. &
&   (mpi_enreg%me==0 .and. mpi_enreg%paral_compil_fft==0 ) .or. &
&   (mpi_enreg%paral_compil_fft==1 .and. mpi_enreg%me_band==0 .and. mpi_enreg%me_kpt==0)) then ! if master
     if (dtset%pawprtden/=3) then
       if (dtset%pawprtden==2.or.dtset%pawprtden==5) rhor_paw = rhor_paw + rhor_paw_core
!      PAWDEN
       rdwr=2 ; fformr=52 ; rdwrpaw=0
       call ioarr(accessfil, rhor_paw, dtset, etotal, fformr, dtfil%fnameabo_app_pawden, hdr, mpi_enreg, &
&       pawfgr%nfft, pawrhoij_dum, rdwr, rdwrpaw)
       if ( accessfil == 3 ) then
!        Complete the geometry informations with missing values from hdr_io().
         call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_pawden, psps)
!        Complete the electrons definition with missing values from hdr_io().
         call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_pawden)
       end if
     end if
     if (dtset%pawprtden>2) then
!      ATMDEN_CORE
       rdwr=2 ; fformr=52 ; rdwrpaw=0
       call ioarr(accessfil, rhor_paw_core, dtset, etotal, fformr, dtfil%fnameabo_app_atmden_core, hdr, mpi_enreg, &
&       pawfgr%nfft, pawrhoij_dum, rdwr, rdwrpaw)
       if ( accessfil == 3 ) then
         call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_atmden_core, psps)
         call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_atmden_core)
       end if
!      We create the file name for valence protodensity. ATMDEN_VAL
       call ioarr(accessfil, rhor_paw_val, dtset, etotal, fformr, dtfil%fnameabo_app_atmden_val, hdr, mpi_enreg, &
&       pawfgr%nfft, pawrhoij_dum, rdwr, rdwrpaw)
       if ( accessfil == 3 ) then
         call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_atmden_val, psps)
         call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_atmden_val)
       end if
!      We create the file name for full protodensity. ATMDEN_FULL
       rhor_paw_val = rhor_paw_val + rhor_paw_core
       call ioarr(accessfil, rhor_paw_val, dtset, etotal, fformr, dtfil%fnameabo_app_atmden_full, hdr, mpi_enreg, &
&       pawfgr%nfft, pawrhoij_dum, rdwr, rdwrpaw)
       if ( accessfil == 3 ) then
         call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_atmden_full, psps)
         call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_atmden_full)
       end if
     end if
   end if ! if master
   if (allocated(rhor_paw)) deallocate(rhor_paw)
   if (allocated(rhor_paw_core)) deallocate(rhor_paw_core)
   if (allocated(rhor_paw_val)) deallocate(rhor_paw_val)
 end if ! if paw+pawprtden

 if(  mpi_enreg%paral_compil_kpt==0                         .or. &
& (mpi_enreg%me==0 .and. mpi_enreg%paral_compil_fft==0 ) .or. &
& (mpi_enreg%paral_compil_fft==1 .and. mpi_enreg%me_band==0 .and. mpi_enreg%me_kpt==0)) then ! if master

!  We output the electron localization function ELF
   if (prtelf/=0) then
     rdwr=2 ; fformr=52 ; rdwrpaw=0
     call ioarr(accessfil,elfr, dtset, etotal,fformr,dtfil%fnameabo_app_elf,hdr, mpi_enreg, &
&     nfft,pawrhoij_dum,rdwr,rdwrpaw)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_elf, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_elf)
     end if
     if (nspden==2)then
       allocate(elfr_up(nfft,nspden))
       elfr_up(:,:) = zero
       do ifft=1,nfft
         elfr_up(ifft,1) = elfr(ifft,2)
       end do
!      ELF_UP
       call ioarr(accessfil,elfr_up, dtset, etotal,fformr,dtfil%fnameabo_app_elf_up,hdr, mpi_enreg, &
&       nfft,pawrhoij_dum,rdwr,rdwrpaw)
       if ( accessfil == 3 ) then
!        Complete the geometry informations with missing values from hdr_io().
         call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_elf_up, psps)
!        Complete the electrons definition with missing values from hdr_io().
         call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_elf_up)
       end if
       allocate(elfr_down(nfft,nspden))
       elfr_down(:,:) = zero
       do ifft=1,nfft
         elfr_down(ifft,1) = elfr(ifft,3)
       end do
!      ELF_DOWN'
       call ioarr(accessfil,elfr_down, dtset, etotal,fformr,dtfil%fnameabo_app_elf_down,hdr, mpi_enreg, &
&       nfft,pawrhoij_dum,rdwr,rdwrpaw)
       if ( accessfil == 3 ) then
!        Complete the geometry informations with missing values from hdr_io().
         call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_elf_down, psps)
!        Complete the electrons definition with missing values from hdr_io().
         call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_elf_down)
       end if
       deallocate(elfr_up,elfr_down)
     end if
   end if

!  We output the gradient of density
   if (prtgden/=0) then
     rdwr=2 ; fformr=52 ; rdwrpaw=0
     call ioarr(accessfil,grhor(:,:,1), dtset, etotal,fformr,dtfil%fnameabo_app_gden1,hdr, mpi_enreg, &
&     nfft,pawrhoij_dum,rdwr,rdwrpaw)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_gden1, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_gden1)
     end if
     call ioarr(accessfil,grhor(:,:,2), dtset, etotal,fformr,dtfil%fnameabo_app_gden2,hdr, mpi_enreg, &
&     nfft,pawrhoij_dum,rdwr,rdwrpaw)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_gden2, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_gden2)
     end if
     call ioarr(accessfil,grhor(:,:,3), dtset, etotal,fformr,dtfil%fnameabo_app_gden3,hdr, mpi_enreg, &
&     nfft,pawrhoij_dum,rdwr,rdwrpaw)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_gden3, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_gden3)
     end if
   end if

!  We output the total kinetic energy density KDEN
   if (prtkden/=0) then
     rdwr=2 ; fformr=52 ; rdwrpaw=0
     call ioarr(accessfil,taur, dtset, etotal,fformr,dtfil%fnameabo_app_kden,hdr, mpi_enreg, &
&     nfft,pawrhoij_dum,rdwr,rdwrpaw)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_kden, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_kden)
     end if
   end if

!  We output the Laplacian of density
   if (prtlden/=0) then
     rdwr=2 ; fformr=52 ; rdwrpaw=0
     call ioarr(accessfil,lrhor, dtset, etotal,fformr,dtfil%fnameabo_app_lden,hdr, mpi_enreg, &
&     nfft,pawrhoij_dum,rdwr,rdwrpaw)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_lden, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_lden)
     end if
   end if

!  We handle the output of wavefunctions. WFK
   if (dtset%prtwf == 1) then
!    In ETSF, some geometric informations are required for wave functions files.
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_wfk, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_wfk)
     end if
   end if

!  POT
   if (prtpot>0) then
     rdwr=2 ; fformv=102 ; rdwrpaw=0
!    MJV note: why is accessfil forced to 0???? This disables the writing of ETSF
!    format potentials!
!    
!    set to 1 for netcdf output
     accessfil = 0
     call ioarr(accessfil,vtrial, dtset, etotal,fformv,dtfil%fnameabo_app_pot,hdr, mpi_enreg, &
&     nfft,pawrhoij_dum,rdwr,rdwrpaw)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_pot, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_pot)
     end if
   end if

   if (prtgeo>0) then
     coordn=prtgeo
     if ( accessfil == 3 ) then
       call abi_etsf_geo_put(dtset,dtfil%fnameabo_app, psps)
     else
       call bonds_lgth_angles(coordn,dtfil%fnameabo_app_geo,natom,psps%ntypat,&
&       rprimd,dtset%typat,xred,dtset%znucl)
     end if
   end if

   if (prtcml>0) then
     call prt_cml2(dtfil%fnameabo_app_cml_xml,natom,dtset%nsym,psps%ntypat,&
&     rprimd,dtset%spgroup,dtset%symrel,dtset%tnons,dtset%typat,xred,dtset%znucl)
   end if

!  STM
   if (prtstm>0) then
     rdwr=2 ; fformr=52 ; rdwrpaw=0
!    set to 1 for netcdf output
     call ioarr(accessfil,rhor, dtset, etotal,fformr,dtfil%fnameabo_app_stm,hdr, mpi_enreg, &
&     nfft,pawrhoij_dum,rdwr,rdwrpaw)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_stm, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_stm)
     end if
   end if

   if (prt1dm>0) then
     call out1dm(dtfil%fnameabo_app_1dm,natom,nfft,ngfft,nspden,psps%ntypat,&
&     rhor,rprimd,dtset%typat,ucvol,vtrial,xred,dtset%znucl)
   end if

!  VHA
   if (prtvha>0) then
     rdwr=2 ; fformv=102 ; rdwrpaw=0
!    set to 1 for netcdf output
     allocate(vwork(nfft,nspden))
     do ispden=1,nspden
       vwork(:,ispden)=vhartr(:)
     end do
     call ioarr(accessfil,vwork, dtset, etotal,fformv,dtfil%fnameabo_app_vha,hdr, mpi_enreg, &
&     nfft,pawrhoij_dum,rdwr,rdwrpaw)
     deallocate(vwork)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_vha, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_vha)
     end if
   end if

!  VHXC
   if (prtvhxc>0) then
     rdwr=2 ; fformv=102 ; rdwrpaw=0
!    set to 1 for netcdf output
     allocate(vwork(nfft,nspden))
     do ispden=1,nspden
       vwork(:,ispden)=vhartr(:)+vxc(:,ispden)
     end do
     call ioarr(accessfil,vwork, dtset, etotal,fformv,dtfil%fnameabo_app_vhxc,hdr, mpi_enreg, &
&     nfft,pawrhoij_dum,rdwr,rdwrpaw)
     deallocate(vwork)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_vhxc, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_vhxc)
     end if
   end if

!  VXC
   if (prtvxc>0) then
     rdwr=2 ; fformv=102 ; rdwrpaw=0
!    set to 1 for netcdf output
     call ioarr(accessfil,vxc, dtset, etotal,fformv,dtfil%fnameabo_app_vxc,hdr, mpi_enreg, &
&     nfft,pawrhoij_dum,rdwr,rdwrpaw)
     if ( accessfil == 3 ) then
!      Complete the geometry informations with missing values from hdr_io().
       call abi_etsf_geo_put(dtset, dtfil%fnameabo_app_vxc, psps)
!      Complete the electrons definition with missing values from hdr_io().
       call abi_etsf_electrons_put(dtset, dtfil%fnameabo_app_vxc)
     end if
   end if

 end if ! if master

!Generate DOS using the tetrahedron method
 partial_dos_flag = 0
 if (prtdos>=2.or.pawfatbnd>0) then

   if(prtdos==2)partial_dos_flag = 0
   if(prtdos==3)partial_dos_flag = 1
   m_dos_flag=0
   if (partial_dos_flag==1) m_dos_flag=dtset%prtdosm
   paw_dos_flag=0
   if (psps%usepaw==1.and.partial_dos_flag==1.and.dtset%pawprtdos>=1) paw_dos_flag=1
   fatbands_flag=0
   if(pawfatbnd>0.and.m_dos_flag==0) fatbands_flag=1
   if(m_dos_flag==1.and.pawfatbnd>0)then
     write(message,'(4a)') ch10,&
&     ' chkinp: WARNING -',ch10,&
&     ' pawfatbnd>0  and prtdosm>0 are not compatible '
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
!    to remove this, one should compute everything in the same basis (cubic)
   end if


!  mjv : initialization is needed as mbesslang is used for allocation below
   mbesslang = 1
   if(partial_dos_flag==1.or.fatbands_flag==1)then
     mbesslang = 5
     ndosfraction=dtset%natsph*mbesslang
   else
     ndosfraction = 1
     mbesslang = 0
   end if

!  For other types of partial DOSs, should use a pointer or something
!  to be able to allocate dos_fractions inside partial_dos_fractions. XG20030506 : Mmmm... not sure !
   allocate(dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction))
   if (m_dos_flag==1.or.fatbands_flag==1) then
     allocate(dos_fractions_m(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*mbesslang))
     allocate(dos_fractions_average_m(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*mbesslang))
   end if
   if (psps%usepaw==1.and.(partial_dos_flag==1)) &
&   allocate(dos_fractions_paw1(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction),&
&   dos_fractions_pawt1(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction))
   if( partial_dos_flag==1.or.fatbands_flag==1)then
!    Generate fractions for partial DOSs if needed
!    partial_dos 1,2,3,4  give different decompositions
     if ((psps%usepaw==0.or.dtset%pawprtdos/=2).and.partial_dos_flag==1) then
       call partial_dos_fractions(cg,dos_fractions,dos_fractions_m,dtfil,dtset,hdr,mbesslang,mpi_enreg, &
&       m_dos_flag,ndosfraction,partial_dos_flag,wffnow)
     else
       dos_fractions=zero;if (m_dos_flag==1.or.fatbands_flag==1) dos_fractions_m=zero
     end if
     if (psps%usepaw==1) then
       call partial_dos_fractions_paw(atindx1,cprj,dimcprj,dos_fractions,dos_fractions_m,&
&       dos_fractions_paw1,dos_fractions_pawt1,dtfil,dtset,fatbands_flag,psps%indlmn,&
&       psps%lmnmax,mbesslang,mkmem,mpi_enreg,m_dos_flag,ndosfraction,&
&       paw_dos_flag,pawrad,pawtab)
     end if
     if(m_dos_flag==1)then
       call dos_degeneratewfs(dos_fractions_m,dos_fractions_average_m,&
&       eigen,mband,dtset%nband,ndosfraction*mbesslang,dtset%nkpt,dtset%nsppol)
     end if
   else
     dos_fractions(:,:,:,1)=one
   end if

!  Here, computation of fatbands for the k-point given. _FATBANDS
   if(pawfatbnd>0.and.fatbands_flag==1) then
     call prtfatbands (dos_fractions_m,dtset,dtfil%fnameabo_app_fatbands,fermie,eigen,&
&     mbesslang,m_dos_flag,ndosfraction,pawfatbnd,pawtab)
   end if

!  Here, computation and output of DOS and partial DOS  _DOS
   if(fatbands_flag==0) then
     call tetrahedron (dos_fractions,dos_fractions_average_m,dos_fractions_paw1,dos_fractions_pawt1,&
&     dtset,fermie,eigen,dtfil%fnameabo_app_dos,mbesslang,m_dos_flag,ndosfraction,paw_dos_flag,rprimd)
   end if

   deallocate(dos_fractions)
   if (m_dos_flag==1.or.fatbands_flag==1) deallocate(dos_fractions_m,dos_fractions_average_m)
   if (psps%usepaw==1.and.(partial_dos_flag==1))&
&   deallocate(dos_fractions_paw1,dos_fractions_pawt1)

 end if ! prtdos > 1

!Output of integrated density inside atomic spheres
 if (dtset%prtdensph==1) call calcdensph(gmet,mpi_enreg,natom,nfft,ngfft,nspden,&
& ntypat,dtset%ratsph,rhor,rprimd,dtset%typat,ucvol,xred)

!If PAW, provide additional outputs
 if (psps%usepaw==1) then
!  Output of compensation charge
   if (usexcnhat>0) then    !if (dtset%nstep>0.or.dtfil%ireadwf/=0)
     write(message, '(4a)' )ch10,' PAW TEST:',ch10,&
&     ' ==== Compensation charge inside spheres ============'
     if (compch_sph>-1.d4.and.compch_fft>-1.d4) &
&     write(message, '(3a)' ) trim(message),ch10,&
     ' The following values must be close to each other ...'
     if (compch_sph>-1.d4) write(message, '(3a,f22.15)' ) trim(message),ch10,&
&     ' Compensation charge over spherical meshes = ',compch_sph
     if (compch_fft>-1.d4) then
       if (pawfgr%usefinegrid==1) then
         write(message, '(3a,f22.15)' ) trim(message),ch10,&
&         ' Compensation charge over fine fft grid    = ',compch_fft
       else
         write(message, '(3a,f22.15)' ) trim(message),ch10,&
&         ' Compensation charge over fft grid         = ',compch_fft
       end if
     end if
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if
!  Output of pseudopotential strength Dij and augmentation occupancies Rhoij
   call pawprt(dtset,psps%indlmn,psps%lmnmax,paw_ij,pawrhoij,pawtab,electronpositron=electronpositron)
 end if

!PAW + output for optical conductivity   _OPT and _OPT2
 if (psps%usepaw==1.and.prtnabla>0) then
   call optics_paw(atindx1,cg,cprj,dimcprj,dtfil,dtset,gprimd,psps%indlmn,kg,psps%lmnmax,&
&   mband,mkmem,mpi_enreg,mpsang,mpw,natom,nkpt,npwarr,nspinor,nsppol,pawrad,pawtab,wffnow)
   if (prtnabla>1) then
     call optics_paw_core(atindx1,cprj,dimcprj,dtfil,dtset,psps%indlmn,psps%lmnmax,&
&     mband,mkmem,mpi_enreg,mpsang,natom,nkpt,nspinor,nsppol,pawrad,pawtab)
   end if
 end if
 if (prtnabla<0) then
  call optics_vloc(cg,dtfil,dtset,gprimd,kg,&
&      mband,mkmem,mpi_enreg,mpw,nkpt,npwarr,nspinor,nsppol,wffnow)
 end if

!Optionally provide output for AE wavefunctions (only for PAW)
 if (psps%usepaw==1 .and. dtset%pawprtwf==1) then
   call pawmkaewf(Dtset,natom,mpw,nspinor,mband,nkpt,mkmem,nsppol,ntypat,Dtset%nband,Dtset%istwfk,npwarr,Dtset%kptns,&
&   Dtset%paral_kgb,Dtset%ngfftdg,kg,dimcprj,Pawfgrtab,Pawrad,Pawtab,&
&   Psps,Hdr,Dtfil,Dtset%typat,eigen,occ,cg,Cprj,Wffnow,MPI_enreg,ierr)
 end if

!Optionally provide output for the GW part of ABINIT
 if (dtset%nbandkss/=0) then
   call timab(233,1,tsec) ! outkss(Total)
   call outkss(dtfil,dtset,ecut,gmet,gprimd,hdr,&
&   dtset%kssform,mband,mgfftc,mkmem,mpi_enreg,mpsang,mpw,natom,&
&   nfft,nkpt,npwarr,nspinor,nspden,nsppol,nsym,psps%ntypat,occ,pawtab,pawfgr,paw_ij,&
&   prtvol,psps,rprimd,wffnow,vtrial,xred,cg,usecprj,cprj,eigen,ierr)
   call timab(233,2,tsec) ! outkss(Total)
   if (ierr/=0) then
     message=' outscfcv: outkss returned a non zero status error, check log'
     call wrtout(ab_out,message,'COLL')
   end if
 end if

!Optionally provide output for  positron life time calculation
 if (electronpositron_calctype(electronpositron)/=0) then
   nzlmopt=dtset%pawnzlm;if (istep_mix<=2) nzlmopt=0
   call poslifetime(dtset,electronpositron,gprimd,mpi_enreg,n3xccc,nfft,ngfft,nzlmopt,&
&   paw_an,pawang,pawrad,pawrhoij,pawtab,rhor,ucvol,xccc3d)
 end if

!Optionally provide output for WanT
 if (dtset%prtwant==1) then
   call outwant(dtfil,dtset,eigen,cg,kg,npwarr,mband,mpi_enreg,nkpt,nsppol,&
&   nspinor,mkmem,mpw,wffnow,dtset%prtwant)
 end if

!Optionally provide output for chemical shielding calculation
 if (prtcs > 0) then
   call calc_cs(dtset%corecs,natom,ntypat,occopt,pawang,pawrad,pawrhoij,pawtab,&
&   prtcs,dtset%typat,psps%usepaw)
 end if

!Optionally provide output for electric field gradient calculation
 if (prtefg > 0) then
   call calc_efg(gprimd,natom,nfft,ngfft,nhat,nspden,ntypat,dtset%paral_kgb,pawang,pawrad,pawrhoij,pawtab,&
&   dtset%ptcharge,prtefg,dtset%quadmom,rhor,rprimd,dtset%typat,ucvol,psps%usepaw,xred,psps%zionpsp)
 end if

!Optionally provide output for Fermi-contact term at nuclear positions
 if (prtfc > 0) then
   call calc_fc(natom,ntypat,pawrad,pawrhoij,pawtab,psps,dtset%typat)
 end if

!Optionally provide Xcrysden output for the Fermi surface
!* Only master enters this part.
 if (dtset%prtfsurf==1.and.MPI_enreg%me==0) then

!  warning; it wont work if nband is not constant!
   allocate(eigen2bxsf(mband,nkpt,nsppol))
   ii=0
   do isppol=1,nsppol
     do ikpt=1,nkpt
       eigen2bxsf(:,ikpt,isppol) = eigen(1+ii:mband+ii)
       ii=ii+mband
     end do
   end do

!  Invert symrel => symrec for k-points
   allocate(symrec(3,3,nsym))
   do isym=1,nsym
     call mati3inv(dtset%symrel(:,:,isym),symrec(:,:,isym))
   end do
!  define whether time reversal can be used or not.
!  FIXME: Here we might have a problem if kptopt==0, we really need to introduce a
!  logical flag defining whether time reversal can be used or not.
   timrev=1 ; if (Dtset%kptopt>=3) timrev=0
   use_afm=(Dtset%nsppol==1.and.Dtset%nspden==2)

!  _BXSF
   call printbxsf(eigen2bxsf,zero,fermie,gprimd,dtset%kptrlatt,mband,dtset%nkpt,hdr%kptns,&
&   nsym,use_afm,symrec,Dtset%symafm,timrev,hdr%nsppol,dtset%shiftk,dtset%nshiftk,dtfil%fnameabo_app_bxsf)
   deallocate(symrec,eigen2bxsf)
 end if ! prtfsurf==1

!DEBUG
!write(std_out,*)' outscfcv : exit'
!stop
!ENDDEBUG

end subroutine outscfcv
!!***
