!{\src2tex{textfont=tt}}
!!****f* ABINIT/gstate
!! NAME
!! gstate
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations by CG minimization.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2010 ABINIT group (DCA, XG, GMR, JYR, MKV, MT, FJ, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  codvsn=code version
!!  cpui=initial CPU time
!!  nspinor=number of spinorial components of the wavefunctions
!!
!! OUTPUT
!!  npwtot(nkpt) = total number of plane waves at each k point
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!
!! SIDE EFFECTS
!!  acell(3)=unit cell length scales (bohr)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | mband =maximum number of bands (IN)
!!   | mgfft =maximum single fft dimension (IN)
!!   | mkmem =maximum number of k points which can fit in core memory (IN)
!!   | mpw   =maximum number of planewaves in basis sphere (large number) (IN)
!!   | natom =number of atoms in unit cell (IN)
!!   | nfft  =(effective) number of FFT grid points (for this processor) (IN)
!!   | nkpt  =number of k points (IN)
!!   | nspden=number of spin-density components (IN)
!!   | nsppol=number of channels for spin-polarization (1 or 2) (IN)
!!   | nsym  =number of symmetry elements in space group
!!  iexit= exit flag
!!  mpi_enreg=MPI-parallelisation information (some already initialized,
!!   some others to be initialized here)
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   Before entering the first time in gstate, a significant part of
!!   psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,
!!     ntypat,n1xccc,usepaw,useylm, and the arrays dimensioned to npsp
!!   All the remaining components of psps are to be initialized in the call
!!   to pspini .
!!   The next time the code enters gstate, psps might be identical to the
!!   one of the previous dtset, in which case, no reinitialisation is scheduled
!!   in pspini.f .
!!  rprim(3,3)=dimensionless real space primitive translations
!!  vel(3,natom)=value of velocity
!!  xred(3,natom) = reduced atomic coordinates
!!
!! NOTES
!! USE OF FFT GRIDS:
!! =================
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
!! In case of wavelets:
!! --------------------
!!    - Only the usual FFT grid (defined by wvl_crmult) is used.
!!      It is defined by nfft, ngfft, mgfft, ... This is strictly not
!!      an FFT grid since its dimensions are not suited for FFTs. They are
!!      defined by wvl_setngfft().
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! TODO
!! Not yet possible to use restartxf in parallel when localrdwf==0
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      blok8,brdmin,bstruct_clean,bstruct_init,chkexi,CleanRec,clnmpi_fft,clnmpi_gs
!!      clnup1,clnup2,delocint,destroy_electronpositron,diisrelax,energies_init,fconv,fixsym,fourdp
!!      getph,handle_ncerr,hdr_clean,hdr_init,hdr_update,indgrid,initberry
!!      initmpi_fft,initmpi_gs,initrhoij,initro,InitRec,initylmg,init_electronpositron,int2char4,inwffil
!!      ioarr,ioddb8,jellium,kpgio,leave_new,mkrho,moldyn,move,newocc,outqmc
!!      outwf,outxfhist,pawinit,pawpuxinit,prtene,psddb8,psolver_kernel,pspini
!!      rhoij_copy,scfcv,scphon,setsym,setsymrhoij,setup1,setup2,status,timab
!!      transgrid,wffclose,wffdelete,wffopen,wffreadskiprec,wrtout
!!      wvl_atoms_data_set_psp,wvl_free_type_proj,wvl_free_type_wfs
!!      wvl_init_type_proj,wvl_init_type_wfs,wvl_mkrho,wvl_setboxgeometry
!!      wvl_setngfft,xcomm_world,xme_init,xmpi_nproc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine gstate(acell,codvsn,cpui,dtfil,dtset,iexit,&
& mpi_enreg,npwtot,nspinor,occ,pawang,pawrad,pawtab,psps,results_gs,rprim,vel,xred)

!=========================== chen =========
 use c_hc_vars
 use m_io_tools, only : flush_unit

 ! for 3D spline
 use EZspline_obj
 use EZspline
!========================================== 


 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_scftypes
 use defs_wvltypes
 use defs_parameters
 use defs_rectypes
 use m_errors
 use m_wffile
 use m_rec
#if defined HAVE_NETCDF
 use netcdf
#endif

 use m_paw_dmft,         only : init_sc_dmft,destroy_sc_dmft,print_sc_dmft,paw_dmft_type
 use m_electronpositron, only : electronpositron_type,init_electronpositron,destroy_electronpositron, &
&                               electronpositron_calctype
 use m_header,           only : hdr_init, hdr_clean
 use m_ebands,           only : bstruct_init, bstruct_clean

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_12_hide_mpi
 use interfaces_14_hidewrite
 use interfaces_16_hideleave
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_53_abiutil
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_57_iovars
 use interfaces_59_io_mpi
 use interfaces_62_iowfdenpot
 use interfaces_62_occeig
 use interfaces_62_poisson
 use interfaces_62_wvl_wfs
 use interfaces_65_psp
 use interfaces_66_paw
 use interfaces_66_wfs
 use interfaces_67_common
 use interfaces_72_response
 use interfaces_79_seqpar_mpi
 use interfaces_95_drive, except_this_one => gstate
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: iexit,nspinor
 real(dp),intent(in) :: cpui
 character(len=6),intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(inout) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
 type(results_gs_type),intent(inout) :: results_gs
!arrays
 integer,intent(out) :: npwtot(dtset%nkpt)
 real(dp),intent(inout) :: acell(3),occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(inout) :: rprim(3,3),vel(3,dtset%natom),xred(3,dtset%natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!Define file format for different type of files. Presently,
!only one file format is supported for each type of files, but this might
!change soon ...
!2   for wavefunction file, new format (version 2.0 and after)    (fform)   NOT USED
!52  for density rho(r)       (fformr)
!102 for potential V(r) file. (fformv)  NOT USED
!scalars
 integer,parameter :: formeig=0,level=101,ndtpawuj=0,response=0
 integer,save :: nsym_old=-1
 integer :: accessfil,ask_accurate,bantot,blktyp,choice,fformr=52,fullinit
 integer :: gscase,iapp,iatom,idir,ierr,ii,indx
 integer :: initialized,ionmov,ios,ir,iscf,itime,itimexit,itypat
 integer :: ixfh,ixx,izero,master,me,mgfftf,mpert,msize,mu,mxfh,nblok
 integer :: ncerr,ncid_hdr,nfftf,nfftftot,nfftot,nproc,nspden_rhoij
 integer :: ntime,ntypat,nxfh,nxfhr,openexit,option,optorth,prtvol,psp_gencond
 integer :: pwind_alloc,rdwr,rdwrpaw,restartxf,spaceworld,tim_mkrho
 integer :: vrsddb
 real(dp) :: cpus,diecut_eff,ecore,ecut_eff,ecutdg_eff,etot,fermie
 real(dp) :: gsqcut_eff,gsqcutc_eff,residm,tolwfr,ucvol
 logical :: read_wf_or_den,use_scf_history
 character(len=500) :: message
 character(len=fnlen) :: ddbnm,dscrpt,filnam
 type(bandstructure_type) :: bstruct
 type(efield_type) :: dtefield
 type(electronpositron_type),pointer :: electronpositron
 type(hdr_type) :: hdr
 type(macro_uj_type) :: dtpawuj(ndtpawuj)
 type(paw_dmft_type) :: paw_dmft
 type(pawfgr_type) :: pawfgr
 type(recursion_type) ::rec_set
 type(scf_in_type) :: scf_in
 type(scf_history_type) :: scf_history
 type(wffile_type) :: wff1,wffnew,wffnow
 type(wvl_data) :: wvl

!arrays
 integer,save :: paw_gencond(6)=(/-1,-1,-1,-1,-1,-1/)
 integer :: ngfft(18),ngfftf(18)
 integer,allocatable :: atindx(:),atindx1(:),blkflg(:),indsym(:,:,:)
 integer,allocatable :: irrzon(:,:,:),kg(:,:),nattyp(:),symrec(:,:,:)
 integer,allocatable,target :: npwarr(:)
 integer,pointer :: npwarr_(:),pwind(:,:,:)
 real(dp) :: blknrm(3),blkqpt(9),gmet(3,3),gprimd(3,3)
 real(dp) :: rmet(3,3),rprimd(3,3),tsec(2)
 real(dp),allocatable :: amass(:),blkval(:,:),cg(:,:),doccde(:)
 real(dp),allocatable :: eigen(:),ph1df(:,:),phnons(:,:,:),resid(:),rhowfg(:,:)
 real(dp),allocatable :: rhowfr(:,:),spinat_dum(:,:),start(:,:),work(:)
 real(dp),allocatable :: xfhist(:,:,:,:),xred_old(:,:)
 real(dp),allocatable :: ylm(:,:),ylmgr(:,:,:)
 real(dp),pointer :: kernel_dummy(:),pwnsfac(:,:),rhog(:,:),rhor(:,:),taug(:,:),taur(:,:)
 type(pawrhoij_type),allocatable :: pawrhoij(:)


!================= chen ==============
 character(len=500)   :: msg,structure,path
 integer,dimension(8) :: dateTime
 integer :: invt_quit,ix,iy,iz,idx,i,j,k,myicount,out_step,target_step, nline=0
 integer :: n1,n2,n3,maxcount,tmp_maxcount,tmp_out_step
 integer :: filestatus , & 
            abnormal_exit,& 
            env_status, & 
            cluster_status

 !------------------------------
 integer :: optmethod = 1

 ! BFGS ----------------------
 integer,parameter :: mmax=5
 integer           :: BFGS_print=-1
 character(len=500):: task,BFGS_csave
 integer           :: BFGS_isave(44)
 integer,allocatable :: BFGS_nbd(:)
 real(dp),allocatable :: BFGS_l(:),BFGS_u(:)
 real(dp),allocatable :: BFGS_wa(:)
 real(dp),allocatable :: BFGS_iwa(:)
 real(dp) :: BFGS_dsave(29)
 logical  :: BFGS_lsave(4)

 ! conjugate Gradient --------------
 integer  :: isave(2),cg_stat
 real(dp) :: norm_dk,norm_gk,norm_y,dTy,cg_beta,grad_stp,stp,dsave(13)
 real(dp) :: xtol,ftol,gtol,stpmax,stpmin

 character(len=300) :: & 
   filename

 !Penalty function related ...
 real(dp) :: & 
   penLambda = -1.d0,        & 
   penfun    = 0.d0

 real(dp) :: & 
   cgW, & 
   stopTol, &
   etotal_cluster, & 
   etotal_env, &
   ek_cluster, ek_env, & 
   enl_cluster, enl_env, &
   eeig_cluster, & 
   eeig_env

 real(dp),allocatable :: & 
   extpot_c(:), &       ! extpot on coarse grid
   yk(:),  &            ! for conjugate gradient
   cgD(:), &            ! for conjugate gradient
   cg_extpot0(:), &     ! for conjugate gradient
   cgG_old(:), &        ! for conjugate gradient
   cgG(:), & 
   cgG_fg(:), &         ! gradient on fine grid
   cgG_fg3d(:,:,:),&        ! gradient on fine grid
   lapvfg(:,:), &
   ref_rhor(:), &       ! for densities
   cluster_den(:), & 
   env_den(:), & 
   tmpfg(:,:), & 
   tmpg(:,:)

 !======= for mpi ======
 logical :: i_am_master
 integer :: mpi_ier

 !======= back up vars before scfcv ====
 type(scf_history_type) :: bk_scf_history
 integer :: bk_initialized,nfft

 !======= for EZspline =====
 type(EZspline3_r8) :: f_spl ! 3-d object/real*8
 integer spn1, spn2, spn3, sp_ier, bcs1(2), bcs2(2), bcs3(2)
 real*8 :: x_int, y_int, z_int, f_int

! ***********************************************************************
!============== chen ========
 nfft = dtset%nfft

 ! for 3D-spline
 spn1=dtset%ngfftdg(1)+1
 spn2=dtset%ngfftdg(2)+1
 spn3=dtset%ngfftdg(3)+1
 bcs1=-1 ! boundary conditions, -1: periodic cond.
 bcs2=-1 ! boundary conditions
 bcs3=-1 ! boundary conditions
!============================ 

 DBG_ENTER("COLL")

 call timab(32,1,tsec)
 call timab(33,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'enter         ')

 if (mpi_enreg%me == 0 .and. dtset%prtxml == 1) then
!  gstate() will handle a dataset, so we output the dataSet markup.
   write(ab_xml_out, "(A)") '  <dataSet>'
!  We output the variables of the dataset given in argument.
!  call outvarsXML()
 end if

!Set up mpi informations from the dataset
 if (dtset%usewvl == 0) then
   mpi_enreg%paral_level=2
   call initmpi_gs(dtset,mpi_enreg)
   call initmpi_fft(dtset,mpi_enreg)
   call initmpi_atom(dtset,mpi_enreg)
 end if

!Define FFT grid(s) sizes (be careful !)
!See NOTES in the comments at the beginning of this file.
 ngfft(:)=dtset%ngfft(:)
 if (psps%usepaw==1) then
   if (dtset%pawecutdg >= 1.0000001_dp*dtset%ecut) then
     pawfgr%usefinegrid=1
     nfftf=dtset%nfftdg;mgfftf=dtset%mgfftdg;ngfftf(:)=dtset%ngfftdg(:)
     nfftot=ngfft(1)*ngfft(2)*ngfft(3)
     nfftftot=ngfftf(1)*ngfftf(2)*ngfftf(3)
     allocate(pawfgr%coatofin(nfftot),pawfgr%fintocoa(nfftftot))
     call indgrid(pawfgr%coatofin,pawfgr%fintocoa,nfftot,nfftftot,ngfft,ngfftf)
   else !this is a simple transfer, this can be done in parallel with only local info
     pawfgr%usefinegrid=0
     nfftf=dtset%nfft;mgfftf=dtset%mgfft;ngfftf(:)=dtset%ngfft(:)
     allocate(pawfgr%coatofin(dtset%nfft),pawfgr%fintocoa(dtset%nfft))
     do ii=1,dtset%nfft
       pawfgr%coatofin(ii)=ii;pawfgr%fintocoa(ii)=ii
     end do
   end if
   pawfgr%nfftc=dtset%nfft;pawfgr%mgfftc=dtset%mgfft;pawfgr%ngfftc(:)=dtset%ngfft(:)
   pawfgr%nfft =nfftf     ;pawfgr%mgfft=mgfftf      ;pawfgr%ngfft(:)=ngfftf(:)
   ecutdg_eff = dtset%pawecutdg * (dtset%dilatmx)**2
   ecut_eff   = dtset%ecut      * (dtset%dilatmx)**2
 else
   pawfgr%usefinegrid=0
   nfftf=dtset%nfft;mgfftf=dtset%mgfft;ngfftf(:)=dtset%ngfft(:)
   allocate(pawfgr%coatofin(0),pawfgr%fintocoa(0))
   ecut_eff= dtset%ecut * (dtset%dilatmx)**2
   ecutdg_eff=ecut_eff
 end if

!
!If dtset%accesswff == 2 set all array outputs to netcdf format
!
 accessfil = 0
 if (dtset%accesswff == 1) then
   accessfil = 4
 end if
 if (dtset%accesswff == 2) then
   accessfil = 1
 end if
 if (dtset%accesswff == 3) then
   accessfil = 3
 end if

!Init spaceworld
 call xcomm_world(mpi_enreg,spaceworld)
 master =0
!Define nproc
 call xmpi_nproc(nproc,ierr)

!Structured debugging if prtvol==-level
 prtvol=dtset%prtvol
 if(prtvol==-level)then
   write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,&
&   ' gstate : enter , debug mode '
   call wrtout(std_out,message,'COLL')
 end if

 ntime=dtset%ntime

!Option input variables
 ionmov   =dtset%ionmov
 iscf     =dtset%iscf
 restartxf=dtset%restartxf

 call status(0,dtfil%filstat,iexit,level,'call setup1   ')

 initialized=0
 ecore=zero

 results_gs%grewtn(:,:)=zero
 call energies_init(results_gs%energies)
 results_gs%pel(1:3)   =zero

!Set up for iterations
 allocate(amass(dtset%natom))
 call setup1(acell,amass,dtset%amu,bantot,&
& ecutdg_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,dtset%iboxcut,dtset%intxc,&
& dtset%natom,dtset%nband,ngfftf,ngfft,dtset%nkpt,dtset%nqpt,dtset%nsppol,psps%ntypat,&
& dtset%qptn,response,rmet,rprim,rprimd,dtset%typat,ucvol,psps%usepaw)

 allocate(npwarr(dtset%nkpt))
 if (dtset%usewvl == 0 .and. dtset%tfkinfunc /= 2) then
   call status(0,dtfil%filstat,iexit,level,'call kpgio    ')
!  Set up the basis sphere of planewaves
   allocate(kg(3,dtset%mpw*dtset%mkmem))
   if (mpi_enreg%mode_para=='b') then
     call prep_kpgio(dtset%accesswff,ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg,dtset%kptns,&
&     dtfil%fnametmp_kg,dtset%mgfft,dtset%mkmem,'PERS',mpi_enreg,dtset%mpw,dtset%nband,dtset%nkpt,&
&     npwarr,npwtot,dtset%nsppol,dtfil%unkg)
   else
     call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg,dtfil%fnametmp_kg, &
&     dtset%kptns,dtset%mkmem,dtset%nband,dtset%nkpt,'PERS',mpi_enreg,&
&     dtset%mpw,npwarr,npwtot,dtset%nsppol,dtfil%unkg)
   end if
 else
   npwarr(:) = 0
   npwtot(:) = 0
 end if

!Set up the Ylm for each k point
 if ( dtset%tfkinfunc /= 2) then
   allocate(ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm))
   allocate(ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm))
   if (psps%useylm==1) then
     if(dtset%mkmem==0) open(dtfil%unylm,file=dtfil%fnametmp_ylm,form='unformatted',status='unknown')
     call status(0,dtfil%filstat,iexit,level,'call initylmg ')
     option=0
     if (dtset%prtstm==0.and.iscf>0.and.dtset%positron/=1) option=1
     call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,&
&     psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
&     npwarr,dtset%nsppol,option,rprimd,dtfil%unkg,dtfil%unylm,ylm,ylmgr)
   end if
 end if
 call timab(33,2,tsec)
 call timab(701,1,tsec)

!Open and read pseudopotential files
 call status(0,dtfil%filstat,iexit,level,'call pspini   ')

 call pspini(dtset,ecore,psp_gencond,gsqcutc_eff,gsqcut_eff,level,pawrad,pawtab,psps,rprimd)

 call timab(701,2,tsec)
 call timab(33,1,tsec)

!In case of isolated computations, ecore must set to zero
!because its contribution is counted in the ewald energy
!as the ion-ion interaction.
 if (dtset%icoulomb == 1) ecore = zero

!WVL - Now that psp data are available, we compute rprimd, acell...
!from the atomic positions.
 if (dtset%usewvl == 1) then
   call wvl_atoms_data_set_psp(dtset, psps)
   call wvl_setBoxGeometry(acell, dtset, mpi_enreg%me, psps%gth_params%radii_cf, &
&   rprimd, xred)
   rprim(:, :)    = reshape((/ &
&   real(1., dp), real(0., dp), real(0., dp), &
&   real(0., dp), real(1., dp), real(0., dp), &
&   real(0., dp), real(0., dp), real(1., dp) /), (/ 3, 3 /))
   call wvl_setngfft(dtset, mpi_enreg)
   nfftf          = dtset%nfft
   mgfftf         = dtset%mgfft
   ngfftf(:)      = dtset%ngfft(:)
 end if

!Initialize band structure datatype
 allocate(doccde(bantot),eigen(bantot))
 doccde(:)=zero ; eigen(:)=zero
 if (dtset%paral_kgb/=0) then     !  We decide to store total npw in bstruct,
   allocate(npwarr_(dtset%nkpt))   !  instead of npw on current proc
   npwarr_(:)=npwarr(:)
   call xsum_mpi(npwarr_,mpi_enreg%commcart,ierr)
 else
   npwarr_ => npwarr
 end if
 call bstruct_init(bantot,bstruct,dtset%nelect,doccde,eigen,dtset%istwfk,dtset%kptns,&
& dtset%nband,dtset%nkpt,npwarr_,dtset%nsppol,dtset%nspinor,dtset%tphysel,&
& dtset%tsmear,dtset%occopt,occ,dtset%wtk)
 deallocate(doccde,eigen)
 if (dtset%paral_kgb/=0) deallocate(npwarr_)
 nullify(npwarr_)

!Initialize PAW atomic occupancies
 if (psps%usepaw==1) then
   allocate(pawrhoij(mpi_enreg%natom))
   nspden_rhoij=dtset%nspden;if (dtset%pawspnorb>0.and.nspinor==2) nspden_rhoij=4
   call initrhoij(dtset%pawcpxocc,psps%indlmn,dtset%lexexch,psps%lmnmax,&
&   dtset%lpawu,mpi_enreg,dtset%natom,mpi_enreg%natom,nspden_rhoij,&
&   dtset%nsppol,dtset%ntypat,pawrhoij,pawtab,dtset%spinat,dtset%typat)
 else
!  it should be sufficient to add the following line (+ deallocation) to avoid an uninitialized jump in hdr_update
   allocate(pawrhoij(0))
 end if

!Initialize header
 gscase=0
 call hdr_init(bstruct,codvsn,dtset,hdr,pawtab,gscase,psps)

!Update header, with evolving variables, when available
!Here, rprimd, xred and occ are available
 etot=hdr%etot ; fermie=hdr%fermie ; residm=hdr%residm
 call hdr_update(bantot,etot,fermie,hdr,dtset%natom,&
& residm,rprimd,occ,mpi_enreg,pawrhoij,psps%usepaw,xred)

!Clean band structure datatype (should use it more in the future !)
 call bstruct_clean(bstruct)

 call status(0,dtfil%filstat,iexit,level,'call inwffil  ')

 allocate(cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol))
 allocate(eigen(dtset%mband*dtset%nkpt*dtset%nsppol))
 allocate(resid(dtset%mband*dtset%nkpt*dtset%nsppol))
 eigen(:)=0.0_dp ; resid(:)=0.0_dp
!mpi_enreg%paralbd=0 ; ask_accurate=0
 ask_accurate=0
!WVL - Branching, allocating wavefunctions as wavelets.
 if (dtset%usewvl == 1) then
!  Create access arrays for wavefunctions and allocate wvl%wfs%psi (other arrays
!  are left unallocated).
   call wvl_init_type_wfs(dtset, mpi_enreg, psps, rprimd, wvl%wfs, xred)
!  We transfer wavelets informations to the hdr structure.
#if defined HAVE_BIGDFT
   hdr%nwvlarr(1) = wvl%wfs%keys%nvctr_c
   hdr%nwvlarr(2) = 7 * wvl%wfs%keys%nvctr_f
#endif
!  Create access arrays for projectors and allocate them.
!  Compute projectors from each atom.
   call wvl_init_type_proj(dtset, mpi_enreg, wvl%projectors, psps, rprimd, xred)
 end if

 read_wf_or_den=(iscf<=0.or.dtfil%ireadwf/=0.or.(dtfil%ireadden/=0.and.dtset%positron<=0))

!RECURSION -  initialization
 if(initialized==0 .and. dtset%userec==1) then
   call InitRec(dtset,mpi_enreg,nfftf,rec_set,maxval(psps%indlmn(3,:,:)))
 end if

!Initialize wavefunctions.
!Warning : ideally, results_gs%fermie and results_gs%residm
!should not be initialized here. One might make them separate variables.
 if(dtset%tfkinfunc /=2) then
   wff1%unwff=dtfil%unwff1
   optorth=1   !if (psps%usepaw==1) optorth=0
   if(psps%usepaw==1 .and. dtfil%ireadwf==1)optorth=0

   call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen,dtset%exchn2n3d,&
&   formeig,gmet,hdr,dtfil%ireadwf,dtset%istwfk,kg,dtset%kptns,&
&   dtset%localrdwf,dtset%mband,dtset%mkmem,mpi_enreg,&
&   dtset%mpw,dtset%nband,ngfft,dtset%nkpt,npwarr,&
&   nspinor,dtset%nsppol,dtset%nsym,occ,&
&   optorth,psps,prtvol,rprimd,dtset%symafm,dtset%symrel,dtset%tnons,&
&   dtfil%unkg,wff1,wffnow,dtfil%unwff1,dtfil%unwft1,&
&   dtfil%fnamewffk,dtfil%fnametmp_wf1,wvl)

 end if
 if (psps%usepaw==1.and.dtfil%ireadwf==1)then
   call rhoij_copy(hdr%pawrhoij,pawrhoij,mpi_enreg=mpi_enreg)
!  Has to update header again (because pawrhoij has changed)
!  MT 2007-10-22: Why ?
   call hdr_update(bantot,etot,fermie,hdr,dtset%natom,&
&   residm,rprimd,occ,mpi_enreg,pawrhoij,psps%usepaw,xred)
 end if

!DEBUG
!write(6,*)' gstate : stop after inwffil for test memory leak '
!call hdr_clean(hdr)
!return
!ENDDEBUG

!Initialize xf history (should be put in inwffil)
 nxfh=0
 if(restartxf>=1 .and. dtfil%ireadwf==1)then

!  Should exchange the data about history in parallel localrdwf==0
   if(mpi_enreg%paral_compil_kpt==1 .and. dtset%localrdwf==0)then
     write(message, '(a,a,a,a,a,a)' )ch10,&
&     ' gstate : BUG -',ch10,&
&     '  It is not yet possible to use non-zero restartxf,',ch10,&
&     '  in parallel, when localrdwf=0. Sorry for this ...'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   end if

!  DEBUG
!  write(std_out,*)'gstate before outxfhist'
!  END DEBUG

   allocate(xfhist(3,dtset%natom+4,2,0))
   call outxfhist(nxfh,dtset%natom,mxfh,xfhist,2,wff1,ios)
   deallocate(xfhist)

!  DEBUG
!  write(std_out,*)'gstate after outxfhist'
!  END DEBUG

   if(ios>0)then
     write(message, '(a,a,a,a,a,a)' )ch10,&
&     ' gstate : BUG -',ch10,&
&     '  An error occurred reading the input wavefunction file,',ch10,&
&     '  with restartxf=1.'
     call wrtout(std_out,message,'COLL')
     call leave_new('COLL')
   else if(ios==0)then
     write(message, '(a,a,i4,a)' )ch10,&
&     ' gstate : reading',nxfh,' (x,f) history pairs from input wf file.'
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
   end if
!  WARNING : should check that restartxf is not negative
!  WARNING : should check that restartxf /= only when dtfil%ireadwf is activated
 end if

!Allocate the xf history array : takes into account the existing
!pairs, minus those that will be discarded, then those that will
!be computed, governed by ntime, and some additional pairs
!(needed when it will be possible to use xfhist for move.f)
 mxfh=(nxfh-restartxf+1)+ntime+5
 allocate(xfhist(3,dtset%natom+4,2,mxfh))
!WARNING : should check that the number of atoms in the wf file and natom are the same

!Initialize the xf history array
 if(nxfh>=restartxf .and. nxfh>0)then
!  Eventually skip some of the previous history
   if(restartxf>=2)then
     do ixfh=1,restartxf-1
       call WffReadSkipRec(ios,1,wff1)
     end do
   end if

!  Read and store the relevant history
   nxfhr=nxfh-restartxf+1
   call outxfhist(nxfhr,dtset%natom,mxfh,xfhist,3,wff1,ios)
 end if

!Determine whether SCF history has to be used
 use_scf_history=(ionmov>0.and.dtset%usewvl==0.and. &
& (abs(dtset%iprcch)==5.or.abs(dtset%iprcch)==6))

!Close wff1, if it was ever opened (in inwffil)
 if (dtfil%ireadwf==1) then
   call WffClose(wff1,ierr)
 end if

!Initialize second wavefunction file if needed
 if(dtset%mkmem==0 .and. dtset%nstep/=0) then
   write(message, '(a,i4,a,a)' )&
&   ' gstate about to open unit',dtfil%unwft2,' for file=',trim(dtfil%fnametmp_wf2)
   call wrtout(std_out,message,'PERS')

#if defined HAVE_NETCDF
   if(dtset%accesswff==2) then
!    Create empty netCDF file
     ncerr = nf90_create(path=trim(dtfil%fnametmp_wf2), cmode=NF90_CLOBBER, ncid=ncid_hdr)
     call handle_ncerr(ncerr," create netcdf wavefunction file")
     ncerr = nf90_close(ncid_hdr)
     call handle_ncerr(ncerr," close netcdf wavefunction file")
   else if(dtset%accesswff==3) then
     write (std_out,*) "FIXME: ETSF I/O support in gstate"
   end if
#endif
   call xme_init(mpi_enreg,me)
   call WffOpen(dtset%accesswff,spaceworld,dtfil%fnametmp_wf2,ierr,wffnew,master,me,dtfil%unwft2)
 end if

 call status(0,dtfil%filstat,iexit,level,'call setup2   ')

!DEBUG
!write(std_out,*)'gstate before setup2'
!END DEBUG

!Further setup
 allocate(start(3,dtset%natom))
 call setup2(dtset,iscf,npwtot,start,wvl%wfs,xred)

!Allocation for forces and atomic positions
 allocate(xred_old(3,dtset%natom))
 xred_old = xred

!Do symmetry stuff only for nsym>1
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 allocate(irrzon(nfftot**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 allocate(phnons(2,nfftot**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 irrzon(:,:,:)=0
 allocate(indsym(4,dtset%nsym,dtset%natom),symrec(3,3,dtset%nsym))

 if (dtset%nsym>1) then

   call status(0,dtfil%filstat,iexit,level,'call setsym   ')

   call setsym(indsym,irrzon,iscf,dtset%natom,&
&   nfftot,ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
&   phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)

!  Make sure dtset%iatfix does not break symmetry
   call status(0,dtfil%filstat,iexit,level,'call fixsym   ')

   call fixsym(dtset%iatfix,indsym,dtset%natom,dtset%nsym)

 else

!  The symrec array is used by initberry even in case nsym = 1
   symrec(:,:,1) = 0
   symrec(1,1,1) = 1 ; symrec(2,2,1) = 1 ; symrec(3,3,1) = 1

 end if

!Timing for initialisation period
 call timab(33,2,tsec)
 call timab(34,1,tsec)

!Compute new occupation numbers, in case wavefunctions and eigenenergies
!were read from disk, occupation scheme is metallic (this excludes iscf=-1),
!and occupation numbers are required by iscf
 if( dtfil%ireadwf==1 .and. &
& (dtset%occopt>=3.and.dtset%occopt<=7) .and. &
& (iscf>0 .or. iscf==-3) .and. dtset%positron/=1 ) then

   call status(0,dtfil%filstat,iexit,level,'call newocc   ')
   allocate(doccde(dtset%mband*dtset%nkpt*dtset%nsppol))
!  Warning : ideally, results_gs%entropy should not be set up here XG 20011007
!  Warning : ideally, results_gs%fermie should not be set up here XG 20011007
!  Do not take into account the possible STM bias
   call newocc(doccde,eigen,results_gs%energies%entropy,&
&   results_gs%energies%e_fermie,&
&   dtset%fixmom,dtset%mband,dtset%nband,&
&   dtset%nelect,dtset%nkpt,nspinor,dtset%nsppol,occ,&
&   dtset%occopt,prtvol,zero,dtset%tphysel,dtset%tsmear,dtset%wtk)
   deallocate(doccde)

 else
!  Warning : ideally, results_gs%entropy should not be set up here XG 20011007
   results_gs%energies%entropy=zero
 end if

!Generate an index table of atoms, in order for them to be used
!type after type.
 ntypat=psps%ntypat
 allocate(atindx(dtset%natom),atindx1(dtset%natom),nattyp(ntypat))
 indx=1
 do itypat=1,ntypat
   nattyp(itypat)=0
   do iatom=1,dtset%natom
     if(dtset%typat(iatom)==itypat)then
       atindx(iatom)=indx
       atindx1(indx)=iatom
       indx=indx+1
       nattyp(itypat)=nattyp(itypat)+1
     end if
   end do
 end do

!Compute structure factor phases for current atomic pos:
 if ((.not.read_wf_or_den).or.use_scf_history) then
   allocate(ph1df(2,3*(2*mgfftf+1)*dtset%natom))
   call status(0,dtfil%filstat,iexit,level,'call getph    ')
   call getph(atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)
 end if

!Initialize paw_dmft, even if neither dmft not paw are used
!write(6,*) "dtset%usedmft",dtset%usedmft
 call init_sc_dmft(dtset%dmftbandi,dtset%dmftbandf,dtset%mband,dtset%nband,dtset%nkpt,nspinor,&
& dtset%nsppol,occ,dtset%usedmft,paw_dmft,dtset%usedmft)
!write(6,*) "paw_dmft%use_dmft",paw_dmft%use_dmft

!PAW: 1- Initialize values for several arrays unchanged during iterations
!2- Initialize data for LDA+U
!3- Eventually open temporary storage file
 if(psps%usepaw==1) then
!  1-
   if (psp_gencond==1.or.&
&   paw_gencond(1)/=dtset%pawlcutd .or.paw_gencond(2)/=dtset%pawlmix  .or.&
&   paw_gencond(3)/=dtset%pawnphi  .or.paw_gencond(4)/=dtset%pawntheta.or.&
&   paw_gencond(5)/=dtset%pawspnorb.or.paw_gencond(6)/=dtset%pawxcdev) then
     call timab(553,1,tsec)
     diecut_eff=abs(dtset%diecut)*dtset%dilatmx**2
     call pawinit(diecut_eff,psps%indlmn,dtset%pawlcutd,dtset%pawlmix,psps%lmnmax,psps%mpsang,&
&     dtset%pawnphi,dtset%nsym,dtset%pawntheta,psps%ntypat,&
&     pawang,pawrad,dtset%pawspnorb,pawtab,dtset%pawxcdev)
     paw_gencond(1)=dtset%pawlcutd ; paw_gencond(2)=dtset%pawlmix
     paw_gencond(3)=dtset%pawnphi  ; paw_gencond(4)=dtset%pawntheta
     paw_gencond(5)=dtset%pawspnorb; paw_gencond(6)=dtset%pawxcdev
     call timab(553,2,tsec)
   else
     if (pawtab(1)%has_kij  ==1) pawtab(1:psps%ntypat)%has_kij  =2
     if (pawtab(1)%has_nabla==1) pawtab(1:psps%ntypat)%has_nabla=2
   end if
   if (psp_gencond==1.or.nsym_old/=dtset%nsym) then
     call setsymrhoij(gprimd,pawang%l_max-1,dtset%nsym,dtset%pawprtvol,&
&     rprimd,symrec,pawang%zarot)
     nsym_old=dtset%nsym
   end if
   psps%n1xccc=maxval(pawtab(1:psps%ntypat)%usetcore)
!  2-Initialize and compute data for LDA+U, EXX, or LDA+DMFT
   pawtab(:)%usepawu=0
   pawtab(:)%useexexch=0
   pawtab(:)%exchmix=zero
   if(paw_dmft%use_dmft==1) call print_sc_dmft(paw_dmft,dtset%pawprtvol)
   if (dtset%usepawu>0.or.dtset%useexexch>0.or.paw_dmft%use_dmft>0) then
     call pawpuxinit(dtset%dmatpuopt,dtset%exchmix,dtset%jpawu,dtset%lexexch,dtset%lpawu,&
&     psps%indlmn,psps%lmnmax,dtset%normpawu,ntypat,pawang,dtset%pawprtvol,pawrad,pawtab,dtset%upawu,&
&     dtset%usedmft,dtset%useexexch,dtset%usepawu)
   end if
!  3-Eventually open temporary storage file
   if(dtset%mkmem==0) then
     open(dtfil%unpaw,file=dtfil%fnametmp_paw,form='unformatted',status='unknown')
     rewind(unit=dtfil%unpaw)
   end if
 end if

!DEBUG
!write(std_out,*)'gstate before call of initberry'
!END DEBUG

!Initialize (eventually) electron-positron data
 if (dtset%positron/=0) then
   allocate(electronpositron)
   electronpositron%calctype=0
   call init_electronpositron(dtfil%ireadwf,dtset,electronpositron,mpi_enreg,nfftf,pawrhoij,pawtab)
 else
   nullify (electronpositron)
 end if

!Electric field: initialization
 dtefield%has_qijb = 0
 if ((dtset%berryopt < 0).or.(dtset%berryopt == 4)) then
   nullify(pwind,pwnsfac)
   call initberry(dtefield,dtset,gmet,gprimd,kg,&
&   dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,&
&   dtset%nkpt,npwarr,dtset%nsppol,&
&   dtset%nsym,dtset%ntypat,occ,pawang,pawrad,pawtab,&
&   pwind,pwind_alloc,pwnsfac,rprimd,symrec,psps%usepaw)
 else
   pwind_alloc = 1
   allocate(pwind(pwind_alloc,2,3),pwnsfac(2,pwind_alloc))
 end if


!Get starting charge density : rhor as well as rhog
 allocate(rhog(2,nfftf),rhor(nfftf,dtset%nspden))
!Also initialize the kinetic energy density
 allocate(taug(2,nfftf*dtset%usekden),taur(nfftf,dtset%nspden*dtset%usekden))
 if (iscf>0) then
   if(dtfil%ireadden/=0.and.dtset%positron<=0)then

     rdwr=1;rdwrpaw=psps%usepaw;if(dtfil%ireadwf/=0) rdwrpaw=0
     call ioarr(accessfil,rhor,dtset,results_gs%etotal,fformr,dtfil%fildensin,hdr,&
&     mpi_enreg, nfftf,pawrhoij,rdwr,rdwrpaw)

     if(dtfil%ireadkden/=0 .and. dtset%usekden==1 )then
!      write(6,*)' gstate : before ioarr, 1 '
       call ioarr(accessfil,taur,dtset,results_gs%etotal,fformr,dtfil%filkdensin,hdr,&
&       mpi_enreg, nfftf,pawrhoij,rdwr,rdwrpaw)
     end if
     if (rdwrpaw/=0) then
       call hdr_update(bantot,etot,fermie,hdr,dtset%natom,&
&       residm,rprimd,occ,mpi_enreg,pawrhoij,psps%usepaw,xred)
     end if
!    Compute up+down rho(G) by fft
     allocate(work(nfftf))
     work(:)=rhor(:,1)
     call fourdp(1,rhog,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
     if(dtset%usekden==1)then
       work(:)=taur(:,1)
       call fourdp(1,taug,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
     end if
     deallocate(work)

   else if(dtfil%ireadwf/=0)then
     izero=0
!    Obtain the charge density from wfs that were read previously
!    Be careful: in PAW, rho does not include the compensation
!    density (to be added in scfcv.F90) !
     call status(0,dtfil%filstat,iexit,level,'call mkrho    ')
!    tim_mkrho=1 ; mpi_enreg%paralbd=0
     tim_mkrho=1
     if (psps%usepaw==1) then
       allocate(rhowfg(2,dtset%nfft),rhowfr(dtset%nfft,dtset%nspden))
!      write(6,*) "mkrhogstate"
       call mkrho(cg,dtset,gprimd,irrzon,kg,&
&       mpi_enreg,npwarr,nspinor,occ,paw_dmft,phnons,rhowfg,rhowfr,rprimd,tim_mkrho,ucvol,&
&       dtfil%unkg,wffnow,wvl%wfs)
       call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,pawfgr,rhowfg,rhog,rhowfr,rhor)
       deallocate(rhowfg,rhowfr)
     else
       call mkrho(cg,dtset,gprimd,irrzon,kg,&
&       mpi_enreg,npwarr,nspinor,occ,paw_dmft,phnons,rhog,rhor,rprimd,tim_mkrho,ucvol,&
&       dtfil%unkg,wffnow,wvl%wfs)
       if(dtset%usekden==1)then
         call mkrho(cg,dtset,gprimd,irrzon,kg,&
&         mpi_enreg,npwarr,nspinor,occ,paw_dmft,phnons,taug,taur,rprimd,tim_mkrho,ucvol,&
&         dtfil%unkg,wffnow,wvl%wfs,option=1)
       end if

     end if

   else if(dtfil%ireadwf==0.and.dtset%positron/=1)then

!    Crude, but realistic initialisation of the density
!    There is not point to compute it from random wavefunctions
!    except with wavelets.
     call status(0,dtfil%filstat,iexit,level,'call initro   ')
     if (dtset%usewvl == 0) then
       call initro(atindx,dtset%densty,gmet,gsqcut_eff,psps%usepaw,mgfftf,mpi_enreg,psps%mqgrid_vl,&
&       dtset%natom,nattyp,nfftf,ngfftf,dtset%nspden,ntypat,dtset%paral_kgb,&
&       pawtab,ph1df,psps%qgrid_vl,rhog,rhor,&
&       dtset%spinat,ucvol,psps%usepaw,dtset%ziontypat,dtset%znucl)
!      Update initialized density taking into account jellium slab
       if(dtset%jellslab/=0) then
         option=2; allocate(work(nfftf))
         call jellium(gmet,gsqcut_eff,mpi_enreg,nfftf,ngfftf,dtset%nspden,&
&         option,dtset%paral_kgb,dtset%slabwsrad,rhog,rhor,rprimd,work,dtset%slabzbeg,dtset%slabzend)
         deallocate(work)
       end if ! of usejell
!      Kinetic energy density initialized to zero (used only in metaGGAs ... )
       if(dtset%usekden==1)then
         taur=zero ; taug=zero
       end if
     else if (dtset%usewvl/=0) then
       call wvl_mkrho(dtset, mpi_enreg, occ, rhor, wvl%wfs)
     end if

   end if

 else if ((iscf==-1.or.iscf==-2.or.iscf==-3).and.dtset%positron<=0) then

   call status(0,dtfil%filstat,iexit,level,'call ioarr    ')
!  Read rho(r) from a disk file
   rdwr=1;rdwrpaw=psps%usepaw
!  Note : results_gs%etotal is read here,
!  and might serve in the tddft routine, but it is contrary to the
!  intended use of results_gs ...
!  Warning : should check the use of results_gs%fermie
!  Warning : should check the use of results_gs%residm
!  One might make them separate variables.

   call ioarr(accessfil,rhor,dtset, results_gs%etotal,fformr,dtfil%fildensin,hdr,&
&   mpi_enreg,nfftf,pawrhoij,rdwr,rdwrpaw)
   if(dtfil%ireadkden/=0 .and. dtset%usekden==1)then
     call ioarr(accessfil,taur,dtset, results_gs%etotal,fformr,dtfil%filkdensin,hdr,&
&     mpi_enreg,nfftf,pawrhoij,rdwr,rdwrpaw)
   end if

!  Compute up+down rho(G) by fft
   call status(0,dtfil%filstat,iexit,level,'call fourdp   ')
   allocate(work(nfftf))
   work(:)=rhor(:,1)
   call fourdp(1,rhog,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
   if(dtset%usekden==1)then
     work(:)=taur(:,1)
     call fourdp(1,taug,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
   end if
   deallocate(work)

 else ! Disallowed value for iscf
   write(message,'(a,i12,a)')'  iscf has the disallowed value=',iscf,'.'
   MSG_BUG(message)
 end if

!Debugging : print the different parts of rhor
!MPIWF Warning : this should not be parallelized over space, leave this debugging feature as such.
 if(prtvol==-level)then
   write(message,'(a)') '   ir     rhor(ir)     '
   call wrtout(std_out,message,'COLL')
   do ir=1,nfftf
     if(ir<=11 .or. mod(ir,301)==0 )then
       write(message,'(i5,a,es13.6)')ir,' ',rhor(ir,1)
       call wrtout(std_out,message,'COLL')
       if(dtset%nsppol==2)then
         write(message,'(a,es13.6)')'      ',rhor(ir,2)
         call wrtout(std_out,message,'COLL')
       end if
     end if
   end do
 end if

!If needed, allocate and initialize SCF history variables
 if (use_scf_history) then
!  Allocations
   scf_history%history_size=2
   scf_history%natom=dtset%natom
   scf_history%nfft=nfftf
   scf_history%nspden=dtset%nspden
   scf_history%alpha=zero
   scf_history%beta=zero
   allocate(scf_history%hindex(scf_history%history_size));scf_history%hindex(:)=0
   allocate(scf_history%deltarhor(nfftf,dtset%nspden,scf_history%history_size))
   allocate(scf_history%xreddiff(3,dtset%natom,scf_history%history_size))
   allocate(scf_history%atmrho_last(nfftf))
   if (psps%usepaw==1) allocate(scf_history%pawrhoij(mpi_enreg%natom,scf_history%history_size))
!  If rhor is an atomic density, just store it in history
   if (.not.read_wf_or_den) then
     scf_history%atmrho_last(:)=rhor(:,1)
   else
!    If rhor is not an atomic density, has to compute rho_at(r)
     allocate(rhowfg(2,nfftf),rhowfr(nfftf,1))
     allocate(spinat_dum(3,dtset%natom));spinat_dum=zero
     call initro(atindx,dtset%densty,gmet,gsqcut_eff,psps%usepaw,mgfftf,mpi_enreg,&
&     psps%mqgrid_vl,dtset%natom,nattyp,nfftf,ngfftf,1,ntypat,dtset%paral_kgb,pawtab,&
&     ph1df,psps%qgrid_vl,rhowfg,rhowfr,spinat_dum,ucvol,&
&     psps%usepaw,dtset%ziontypat,dtset%znucl)
     scf_history%atmrho_last(:)=rhowfr(:,1)
     deallocate(rhowfg,rhowfr,spinat_dum)
   end if
   scf_history%usecg=0
   if(dtset%iextrapwf==1) then
     scf_history%usecg=1
     scf_history%icall=0
     allocate(scf_history%cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol,scf_history%history_size))
     allocate(scf_history%cprj(dtset%natom,nspinor*dtset%mband*dtset%mkmem*dtset%nsppol,scf_history%history_size))
   end if
 else
   scf_history%history_size=0
 end if

 if ((.not.read_wf_or_den).or.use_scf_history) deallocate(ph1df)

 scf_in%fatvshift=one

 call status(0,dtfil%filstat,iexit,level,'end gstate(1) ')

 if(prtvol==-level)then
   write(message,'(a1,a,a1,a,i1,a)') ch10,&
&   ' gstate : before scfcv, move or brdmin ',&
&   ch10,'  prtvol=-',level,', debugging mode => stop '
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
 end if

 call timab(34,2,tsec)
!Check whether exiting was required by the user.
!If found then do not start minimization steps
!At this first call to chkexi, initialize cpus, if it
!is non-zero (which would mean that no action has to be taken)
!Should do this in driver ...
 cpus=dtset%cpus
 if(abs(cpus)>1.0d-5)cpus=cpus+cpui
 openexit=1 ; if(dtset%chkexit==0) openexit=0
 call chkexi(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)
!If immediate exit, and wavefunctions were not read, must zero eigenvalues
 if (iexit/=0) then
   eigen(:)=zero
 end if
 if (iexit==0) then

!  #######################################################################

!  If atoms are not being moved and U should not be determined, use scfcv directly; else
!  call move, pawuj_drive or brdmin which in turn calls scfcv.

   call timab(35,1,tsec)

   write(message,'(a,80a)')ch10,('=',mu=1,80)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   if (ionmov==0) then

!    Should merge this call with the call for ionmov==4 and 5
     iapp=0
!    mpi_enreg%paralbd=0

     call status(0,dtfil%filstat,iexit,level,'call scfcv    ')

     if (dtset%macro_uj==0) then







 !====================================================
 !============ CHEN TRAP THE CODE: START =============
 !====================================================

 i_am_master = .true.
 if (nproc>1) then
   if (mpi_enreg%me/=0) then
     i_am_master = .false.
   endif
 endif
 
!!DEBUG 
 call xbarrier_mpi(spaceworld,mpi_ier)
 write(*,'(a,I3,a,g22.14,a,L)') & 
  " (chen/gstate) MPI information >>> mpi_enreg%me:",mpi_enreg%me,&
  "  sum(den)=",sum(rhor)*ucvol/dble(nfftf),& 
  "   master? => ",i_am_master
 call flush_unit(std_out)
 call xbarrier_mpi(spaceworld,mpi_ier)
 if (i_am_master) print *,'(chen) ======= above should only have one TRUE ========'
!!END DEBUG

 allocate(BFGS_nbd(nfft))
 allocate(BFGS_l(nfft),BFGS_u(nfft))
 allocate(BFGS_wa(2*mmax*nfft+4*nfft+12*mmax*mmax+12*mmax))
 allocate(BFGS_iwa(3*nfft))
 allocate(lapvfg(nfftf,1),cgG(nfft),cgG_fg(nfftf),cgG_fg3d(spn1,spn2,spn3))
 allocate(ref_rhor(nfftf),cluster_den(nfftf),env_den(nfftf))
 allocate(extpot_c(nfft),extpot(nfftf))
 allocate(hc_rhor_paw(nfftf))                    ! rhor + rho1-trho1, total true density
 if (psps%usepaw==1) then 
   allocate(in_paw(nfftf))
 endif
 if (optmethod==1 .or. optmethod==2) then ! do conjugate gradient
   allocate(yk(nfft),cgD(nfft),cgG_old(nfft),cg_extpot0(nfft))
 endif
 allocate(tmpfg(2,dtset%nfft),tmpg(2,nfftf))
 
! initialize parameters:
 bAllowsym = 0;
 bLoadVr = 1;

! BFGS ....
 task='START'
 BFGS_nbd = 0   ! unbounded

! CG ----------------
 cg_stat = 1
 ftol = 1e-4  ! new energy should be smaller than the older one.
 gtol = .2d0  ! 0.1 is for CG which needs a more exact line-search,
              ! gtol MUST > ftol to have wolfe condition always satisified
 xtol = 1.E-12
 stpmin = 0.d0
 stpmax = 1e5

 maxcount = 200;
 n1=ngfftf(1); 
 n2=ngfftf(2); 
 n3=ngfftf(3);

! Welcome Messages
 call c_wrtlog(' ')
 call c_wrtlog('WELCOMEMSG_BF')
 call c_wrtlog(' ')
 
!==========================================
!print date time info.
 call date_and_time(VALUES=dateTime)    !today
 write(msg, '(A,i2,A,i2,A,i4,A,i2.2,A,i2.2,A,i2.2,A,i3.3)') "START AT: ", &
   dateTime(3),"/",  dateTime(2),"/",  dateTime(1)," : ",             &
   dateTime(5), ":", dateTime(6), ":", dateTime(7), ".", dateTime(8)
 call c_wrtlog(msg);
 call c_wrtlog('ABINIT parameters:');
 write (message, '(A,I10,A,I5,A,I5,A,I5)')  & 
   'nfftf=',nfftf,' n1=',n1,' n2=',n2,' n3=',n3
 call c_wrtlog(message);
 write (message, '(A,E15.7,A,9I3)') & 
   'ecut =',dtset%ecut,' nkpt=',dtset%kptrlatt
 call c_wrtlog(message);
 write (message, '(a,I2,a,E15.7,a,E15.7)') & 
   'occ  =',dtset%occopt,' tsmear=',dtset%tsmear,' ucvol=', ucvol
 call c_wrtlog(message);
 call c_wrtlog("pseudopotenitals:")
 do ii=1,psps%ntypat
   write(message,'(a,a)')'psp file: ',psps%filpsp(ii)
   call c_wrtlog(message)
 enddo

 ! Get current directory
 CALL c_wrtlog("PATH:")
 CALL GetCWD(path)
 CALL c_wrtlog(path)

 IF ( path(LEN_TRIM(path)-6:LEN_TRIM(path)) .EQ. 'cluster' ) THEN 
   whoAmI = 1
   CALL c_wrtlog(" ======= I am the cluster ======")
 ELSEIF ( path(LEN_TRIM(path)-2:LEN_TRIM(path)) .EQ. 'env' ) THEN 
   whoAmI = 2
   CALL c_wrtlog(" ======= I am the env ======")
 ELSE
   CALL c_wrtlog("ERROR! you must name the folder to be 'cluster' or 'env'. STOP")
   call leave_new('COLL')
 ENDIF

 ! Load parameters 
 CALL LoadPara(whoAmI,bLoadDen,bLoadVr,maxcount, & 
   bNCPP_PHI,penLambda,bAllowsym,& 
   optmethod,structure,stopTol,out_step)

 target_step = out_step

 ! ===== ECHO for input parameters =====!
 call c_wrtlog("")
 CALL C_WRTLOG("User parameters:")
 WRITE(message,'(a,I5)')    "whoAmI        :", whoAmI         ;  CALL C_WRTLOG(message)
 WRITE(message,'(a,I5)')    "bLoadDen      :", bLoadDen       ;  CALL C_WRTLOG(message)
 WRITE(message,'(a,I5)')    "bLoadVr       :", bLoadVr        ;  CALL C_WRTLOG(message)
 WRITE(message,'(a,I5)')    "maxcount      :", maxcount       ;  CALL C_WRTLOG(message)
 WRITE(message,'(a,I5)')    "bNCPP_PHI     :", bNCPP_PHI      ;  CALL C_WRTLOG(message)
 WRITE(message,'(a,ES12.4)')"penLambda     :", penLambda      ;  CALL C_WRTLOG(message)
 WRITE(message,'(a,I5)')    "bAllowsym     :", bAllowsym      ;  CALL C_WRTLOG(message)
 WRITE(message,'(a,I5)')    "optmethod     :", optmethod      ;  CALL C_WRTLOG(message)
 WRITE(message,'(a,ES12.4)')"stopTol       :", stopTol        ;  CALL C_WRTLOG(message)
 WRITE(message,'(a,I5)')    "out_step      :", out_step       ;  CALL C_WRTLOG(message)

 !=== PRINT penalty function infomation ===
 if (optmethod==1) then
   write(message,'(a)') 'use conjugate gradient (Polak and Ribire,i.e. PRP+) beta=max(beta,0), more robust.'
   call c_wrtlog(message)
 else if (optmethod==2) then
   write(message,'(a)') 'use conjugate gradient (Hager and Zhang), more fast.'
   call c_wrtlog(message)
 else if (optmethod==11) then 
   write(message,'(a)') 'use BFGS.'
   call c_wrtlog(message)
 else 
   write(message,*)'optmethod is not defined, optmethod=',optmethod
   call c_wrtlog(message)
   call leave_new('COLL')
 endif
    

 !====== PRINT ALL MPI-INFORMATIONS ======
 call c_wrtlog("")
 call c_wrtlog("MPI informations:")
 write(message,'(a,I4)')"dtset%para_kgb         :",dtset%paral_kgb;  call c_wrtlog(message)
 write(message,'(a,I4)')"mpi_enreg%nproc        :",mpi_enreg%nproc;  call c_wrtlog(message)
 write(message,'(a,I4)')"mpi_enreg%nproc_atom   :",mpi_enreg%nproc_atom;  call c_wrtlog(message)
 call c_wrtlog("")

 !==== some checks ========!
 if (mpi_enreg%nproc_atom>1) then
   call c_wrtlog('mpi_enreg%nproc_atom>1, this is not implemented yet for PAW, stop')
   call leave_new('COLL')
 endif
 
 !==== Check the 'nsym' flag in the ABINIT input file ===
 if (bAllowSym == 0) then
   if (dtset%nsym/=1) then
     call c_wrtlog("Do the symmtry or not, bAllowsym=0 in param.in, but nsym!=1, confused. stop.")
     call leave_new('COLL')
   end if
 else
   call c_wrtlog('You set the bAllowSym = 1')
   if (dtset%nsym /= 1) then
     call c_wrtlog('nsym > 1, symmetry is considered.')
   else
     call c_wrtlog('The dtset%nsym=1 but you set bAllowSym/=0 Confusing. STOP')
     call leave_new('COLL')
   end if
 end if

 !==== Load refden.paw only if i am cluster ====
 IF ( whoAmI == 1) THEN 
   if (bLoadDen == 1) then
     if (psps%usepaw==1) then
       filename = 'refden.paw'   ! PAW's *o_PAWDEN file, with pawprtden=1
     else
       filename = 'refden.ncpp'  ! Norm conserving density file with prtden=1
     endif
     if (i_am_master) call c_checkFileLen(filename,nfftf)
     call xbarrier_mpi(spaceworld,mpi_ier)
     open(unit=1111, access="sequential", action="read", &
       file=filename, form="formatted", status="old");
     read(1111,*)ref_rhor;
     close(1111);
     if (psps%usepaw==1) then
       write(msg,'(a,ES18.10)') "loaded refden.paw, electron number:",& 
       sum(ref_rhor)*ucvol/dble(nfftf)
     else
       write(msg,'(a,ES18.10)') "loaded refden.ncpp, electron number:",& 
       sum(ref_rhor)*ucvol/dble(nfftf)
     endif
     call c_wrtlog(msg)
     write(msg,'(2(a,ES12.4))') & 
       "max(refden)=",maxval(ref_rhor), & 
       " min(refden)=",minval(ref_rhor)
     call c_wrtlog(msg)
   else
     call c_wrtlog('not load the refden.paw, be cautious!')
   endif
 ENDIF

 !=== Initialize vtrial ===
 if (bLoadVr==1) then
   if (whoAmI==1) then 
     filename = 'pot.in'
   else
     filename = '../cluster/pot.in'
   endif
   call c_wrtlog("(gstate) to load pot.in file")
   call c_checkFileLen(filename,nfftf)
   call xbarrier_mpi(spaceworld,mpi_ier)
   open(unit=1111, access="sequential", action="read", & 
     file=filename, form="formatted", status="old");
   do i=1,nfftf
     read(1111,*)extpot(i)
   enddo
   close (1111)
   write(message,*)  & 
     "set pot.in to extpot. max/min(extpot):", & 
     maxval(extpot),minval(extpot)
   extpot_c = extpot
   call c_wrtlog(message)
   call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,tmpfg,tmpg,extpot_c,extpot)
   call xbarrier_mpi(spaceworld,mpi_ier)
 else
   extpot   = 0.d0 !set to zero as the initial guess
   extpot_c = 0.d0 !set to zero as the initial guess
   call c_wrtlog('as an initial guess, extpot and extpot_c are set to zero.')
 endif

! DEBUG
!   open(unit=111,action='read',file='vxc.pot')
!   PRINT *,"DEBUG: chen read in vxc.pot"
!   open(unit=111,action='read',file='vtrial_debug'); 
!   print *, "read in vtrial_debug"
!   read(111,*)extpot
!   close(111)
!   extpot=extpot+10.d0
!   call c_wrtlog('DEBUG ============== extpot=extpot+10.d0 =========')
! END DEBUG

 !--------------------------------------------------------------------------------------
 !----------------------- Embedding Potential Loop -------------------------------------
 !--------------------------------------------------------------------------------------
 myicount = -1
 DO 
1997 continue   

   myicount=myicount+1
   call c_wrtlog('')
   call c_wrtlog('')
   write (message, '(a,i4,a)')  & 
    '============ tried embedding potential for ', &
     myicount,' times ============ ';
   call c_wrtlog(message);       
   tmp_maxcount = maxcount
   tmp_out_step = out_step
   CALL LoadPara(whoAmI,bLoadDen,bLoadVr,maxcount, & 
          bNCPP_PHI,penLambda,bAllowsym,  &
          optmethod,structure,stopTol,out_step)
   if ( tmp_maxcount /= maxcount ) then
     write(message,'(a,I4)') "new maxcount:", maxcount
     CALL c_wrtlog(message)
   endif
   if (tmp_out_step/=out_step) then
     target_step=target_step-tmp_out_step+out_step
     write(message,'(a,I4,a,I4)') "new out_step:",out_step,' next print step:',target_step
     CALL c_wrtlog(message)
   endif


   !back up vars
   bk_initialized = initialized
   bk_scf_history = scf_history   
   
   write(message,'(a,2es12.4,a,es12.4)')&
     "(gstate) to do SCF: max/min(extpot):", maxval(extpot),minval(extpot), & 
     ' max-min:',maxval(extpot)-minval(extpot)
   call c_wrtlog(message)
   open(unit=2222,access='append',action='write',status='unknown', file='invKS_log',form='formatted');
!   call c_print_time(2222)
   close(2222)
   call xbarrier_mpi(spaceworld,mpi_ier)

!=============================
! Chen: ABINIT's SCF 

       call scfcv_tmp(acell,atindx,atindx1,cg,cpus,dtefield,dtfil,dtpawuj,dtset,&
&       ecore,eigen,electronpositron,hdr,iapp,indsym,initialized,&
&       irrzon,kg,mpi_enreg,nattyp,ndtpawuj,nfftf,npwarr,nspinor,occ,&
&       paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
&       pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
&       scf_history,scf_in,symrec,taug,taur,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)

! chen: end of ABINIT's SCF
!=============================

   call xbarrier_mpi(spaceworld,mpi_ier)
   write(message,'(a)') "(gstate) finished SCF ----- "
   call c_wrtlog(message)
   open(unit=2222,access='append',action='write',status='unknown', file='invKS_log',form='formatted');
!   call c_print_time(2222)
   close(2222)

   initialized = bk_initialized
   scf_history = bk_scf_history

   !------------!
   ! Cluster    !
   !------------!
   IF (whoAmI == 1) THEN 
     if (psps%usepaw==1) then
       cluster_den = hc_rhor_paw(:)
       etotal_cluster = results_gs%etotal-ext_energy
     else
       cluster_den = rhor(:,1)
       etotal_cluster = results_gs%etotal
     endif
! tmp codes
     write(message,*) ' cluster_total_energy = ', etotal_cluster
     call c_wrtlog(message)
     write(message,*) ' cluster_den * embpot = ', dot_product(cluster_den,extpot)*ucvol/dble(nfftf)
     call c_wrtlog(message)
     write(message,*) ' E_cluster            = ', etotal_cluster - dot_product(cluster_den,extpot)*ucvol/dble(nfftf), ' (w/o embpot)'
     call c_wrtlog(message)
! end tmp codes

     eeig_cluster = results_gs%energies%e_eigenvalues
     ek_cluster   = results_gs%energies%e_kinetic
     enl_cluster  = results_gs%energies%e_nonlocalpsp

     !=========================================
     ! Read env.status file, sync cluster and env (every node reads)
     filename = 'env.status'
     CALL WaitForFile(filename)
     OPEN(unit=199,file=filename,action='read',iostat=filestatus)
     IF (fileStatus/=0) THEN
       CALL c_wrtlog('(gstate): Error in opening env.status, code stop')
       CALL leave_new('COLL')
     ENDIF
     DO
       read(199,*,iostat=filestatus)env_status, message
       close (199)
       if (filestatus/=0 .or. message(1:3)/='end' ) then
         write(message,'(a,I3,a)')'(gstate): error in reading env.status, iostat=',filestatus,' sleep for 2 sec.'
         call c_wrtlog(message)
         call system('sleep 2')
       else
         exit
       endif
       open(unit=199,file=filename,action='read',iostat=filestatus)
       if (filestatus/=0) then
         call c_wrtlog('(gstate): error in opening env.status (site A), code stop')
         call leave_new('coll')
       endif
     ENDDO ! loop to read env.status file
     IF (env_status==myicount) THEN
       CALL c_wrtlog("cluster: test of synchronization is passed.")
     ELSE
       write(message,*)'cluster iter=',myicount, & 
       '/    env iter=',env_status, & 
        ' sync failed! stop'
       call c_wrtlog(message)
       call leave_new("COLL")
     ENDIF
     call xbarrier_mpi(spaceworld,mpi_ier)
     if (i_am_master) then 
       OPEN(unit=199,file=filename)
       CLOSE(199,status='DELETE')
       call c_wrtlog("env.status is removed")
     endif
     call xbarrier_mpi(spaceworld,mpi_ier)

     
     !======================== 
     ! Read in env.density ( every node reads )
     filename = 'env.density'
     CALL WaitForFile(filename)
     OPEN(unit=199,file=filename,action='read',iostat=filestatus)
     IF (filestatus/=0)THEN
       CALL c_wrtlog('ERROR in reading env.density, stop')
       CALL leave_new('COLL')
      ENDIF
     call c_wrtlog("reading env.density")
     READ(199,*,iostat=filestatus) env_den
     if (filestatus /= 0) then 
       call c_wrtlog('readinng env.density, filestatus/=0, stop')
       CALL leave_new('COLL')
     endif
     CLOSE(199)
     WRITE(message,'(a,2ES16.4)') "(gstate) max/min(cluster_den) :",maxval(cluster_den),minval(cluster_den);CALL c_wrtlog(message)
     WRITE(message,'(a,2ES16.4)') "(gstate) max/min(env_den)     :",maxval(env_den),minval(env_den); CALL c_wrtlog(message)
     call xbarrier_mpi(spaceworld,mpi_ier)
     if (i_am_master) then
       OPEN(unit=199,file=filename)
       CLOSE(199,status='DELETE')
       call c_wrtlog("env.density is removed.")
     endif
     call xbarrier_mpi(spaceworld,mpi_ier)
     write(message,'(a,g22.14,a)')'(gstate) norm of density=',sum(cluster_den)/dble(nfftf)*ucvol,' (cluster)'
     call c_wrtlog(message)
     write(message,'(a,g22.14,a)')'(gstate) norm of density=',sum(env_den)/dble(nfftf)*ucvol, ' (env)'
     call c_wrtlog(message)
     
     !=======================
     ! Read env.eeig file (every node reads)
     filename = 'env.eeig'
     CALL WaitForFile(filename)
     OPEN(unit=199, file=filename,action="read",iostat=filestatus)
     IF (filestatus /=0 )THEN
       CALL c_wrtlog("Error in opening  env.eeig file. error stop!")
       CALL leave_new('COLL')
     ENDIF
     READ(199,*) etotal_env
     READ(199,*) eeig_env
     READ(199,*) ek_env
     READ(199,*) enl_env
     CLOSE(199)
     call xbarrier_mpi(spaceworld,mpi_ier)
     if (i_am_master) then 
       OPEN(unit=199,file=filename)
       CLOSE(199,status='DELETE')
       call c_wrtlog("env.eeig is removed")
     endif
     call xbarrier_mpi(spaceworld,mpi_ier)
     CALL c_wrtLog("loaded and removed env.eeig file")
     WRITE(message,'(a,g22.14)') '   -- etotal_cluster   = ', etotal_cluster;   CALL c_wrtlog(message)
     WRITE(message,'(a,g22.14)') '   -- etotal_env       = ', etotal_env;       CALL c_wrtlog(message)
     WRITE(message,'(a,g22.14)') '   -- eeig_cluster     = ', eeig_cluster;     CALL c_wrtlog(message)
     WRITE(message,'(a,g22.14)') '   -- eeig_env         = ', eeig_env;         CALL c_wrtlog(message)
     WRITE(message,'(a,g22.14)') '   -- ek_cluster       = ', ek_cluster;       CALL c_wrtlog(message)
     WRITE(message,'(a,g22.14)') '   -- ek_env           = ', ek_env;           CALL c_wrtlog(message)
     if (psps%usepaw==0) then
       WRITE(message,'(a,g22.14)') '   -- enl_cluster      = ', enl_cluster;        CALL c_wrtlog(message)
       WRITE(message,'(a,g22.14)') '   -- enl_env          = ', enl_env;            CALL c_wrtlog(message)
     endif
     
   ENDIF  !whoAmI==1
 
   !-----------------!
   ! Enviroment      !
   !-----------------!
   IF ( whoAmI == 2) THEN
     if (psps%usepaw==1) then
       env_den = hc_rhor_paw(:)  
       etotal_env = results_gs%etotal-ext_energy
     else
       env_den = rhor(:,1)
       etotal_env = results_gs%etotal
     endif
     eeig_env = results_gs%energies%e_eigenvalues
     ek_env   = results_gs%energies%e_kinetic
     enl_env  = results_gs%energies%e_nonlocalpsp

     write(message,*) ' env_den * embpot = ', dot_product(env_den,extpot)*ucvol/dble(nfftf)
     call c_wrtlog(message)
     write(message,*) ' E_env            = ', etotal_env - dot_product(env_den,extpot)*ucvol/dble(nfftf), ' (w/o embpot)'
     call c_wrtlog(message)
 
    !=====================
    ! Write density file
    if (i_am_master) then   ! master
      filename = "../cluster/env.density"
      CALL checkfile(filename)  ! if the file exists, wait until it is removed
      print *,"(gstate) ../cluster/env.density is not there, it is OK to create a new one"
      OPEN(unit=199, file='../cluster/env.density',action="write",iostat=filestatus)
      IF (filestatus/=0)THEN
        CALL c_wrtlog(" error in writing env.density file by env. error stop!")
        CALL leave_new('COLL')
      ENDIF
      WRITE(199,*) env_den
      CLOSE(199)
      WRITE(message,'(a,2ES16.4)') "max/min(env_den):",maxval(env_den),minval(env_den)
      CALL c_wrtlog(message)
      write(message,'(a,es16.4)') "norm of env_den:",sum(env_den)/dble(nfftf)*ucvol
      call c_wrtlog(message)
      CALL c_wrtlog("dumped env.density")
    endif
    write(*,'(a,g22.14,a,i3)') "(chen/gstate) norm of env_den:",sum(env_den)/dble(nfftf)*ucvol," me:",mpi_enreg%me
    call xbarrier_mpi(spaceworld,mpi_ier)
    
    !============================
    ! write eeig/ek/enl to file
    if (i_am_master) then
      filename = "../cluster/env.eeig"
      CALL checkfile(filename)  ! if the file exists, wait until it is removed by cluster
      print *,"(gstate) ../cluster/env.eeig is not there, it is OK to create a new one"
      OPEN(unit=199, file='../cluster/env.eeig',action="write",iostat=filestatus)
      IF (filestatus/=0)THEN
        CALL c_wrtlog(" error in writing env.eeig file by env. error stop!")
        CALL leave_new('COLL')
      ENDIF
      WRITE(199,*)etotal_env
      WRITE(199,*)eeig_env
      WRITE(199,*)ek_env
      WRITE(199,*)enl_env 
!      call c_print_time(199)
      CLOSE(199)
      WRITE(message,'(a,g22.14)') '(env) etotal_env = ', etotal_env;    CALL c_wrtlog(message)
      WRITE(message,'(a,g22.14)') '(env) eeig_env   = ', eeig_env;      CALL c_wrtlog(message)
      WRITE(message,'(a,g22.14)') '(env) ek_env     = ', ek_env;        CALL c_wrtlog(message)
      WRITE(message,'(a,g22.14)') '(env) enl_env    = ', enl_env;       CALL c_wrtlog(message)
      CALL c_wrtlog("dumped env.eeig")
    endif  ! master
    call xbarrier_mpi(spaceworld,mpi_ier)
    
    !=========================
    ! write env.status file
    if (i_am_master) then
      call write_status( myicount, 'ENV' )
!      call c_print_time(199)
      CALL c_wrtlog("dumped env.status")
    endif ! master
    call xbarrier_mpi(spaceworld,mpi_ier)

    !========================
    ! Read ../cluster/cluster.status  (every node read)
    filename = '../cluster/cluster.status'
    CALL WaitForFile(filename)
    OPEN(unit=199,file=filename,action='read',iostat=filestatus)
    IF ( filestatus/=0 ) THEN
      CALL c_wrtlog("error in open ../cluster/cluster.status, stop")
      CALL leave_new('COLL')
    ENDIF
    do
      read(199,*,iostat=filestatus)cluster_status , message
      close(199)
      if (filestatus/=0 .or. message(1:3)/='end' ) then
        write(message,'(a,I3,a)')'(gstate): error in reading cluster.status, iostat=',filestatus,' sleep for 2 sec.'
        call c_wrtlog(message)
        call system('sleep 2')
      else
        exit
      endif
      open(unit=199,file=filename,action='read',iostat=filestatus)
      if (filestatus/=0) then
        call c_wrtlog('(gstate): error in opening cluster.status, code stop')
        call leave_new('coll')
      endif
    enddo ! loop to read cluster.status file
    call xbarrier_mpi(spaceworld,mpi_ier)
    if ( cluster_status==1 ) then 
      if (i_am_master) then 
        open(unit=199,file=filename)
        close(199,status='delete')
      endif
      call c_wrtlog('env: ../cluster/cluster.status is removed')
      call c_wrtlog("env: ok message from cluster.status (removed)")
    elseif ( cluster_status==-1 )then 
      if (i_am_master) then
        open(unit=199,file=filename)
        close(199,status='delete')
      endif
      call c_wrtlog('env: ../cluster/cluster.status is removed')
      call c_wrtlog("env: get exit message from cluster.status, leaving...")
      call xbarrier_mpi(spaceworld,mpi_ier)
      goto 99999
    endif
    call xbarrier_mpi(spaceworld,mpi_ier)
    
    !=======================
    ! Read ../cluster/new_extpot.dat (every node reads)
    filename = '../cluster/new_extpot.tmp'
    call waitforfile(filename)
    open(unit=199,file=filename,action='read',iostat=filestatus,form="unformatted")
    if ( filestatus/=0 ) then
      call c_wrtlog('error in opening new_extpot.tmp. code stop.')
      call leave_new('coll')
    endif
    call c_wrtlog("reading new_extpot.tmp ...")
    read(199,iostat=filestatus)extpot(:)
    if (filestatus/=0) then 
      call c_wrtlog("reading new_extpot.tmp, filestatus/=0, error, code stop.")
      close (199)
      call leave_new('coll')
    endif
    close(199)
    write(message,'(a,2es12.4)') & 
      "max/min_extpot:",maxval(extpot),minval(extpot)
    call c_wrtlog(message)
    call c_wrtlog('env: loaded new_extpot.tmp file.')

    call xbarrier_mpi(spaceworld,mpi_ier)
    cycle  ! environment will cycle from here

  endif ! whoAmI == 2
 
  !==============================================================
  ! Code below is for cluster ONLY
  !==============================================================
  
  if (psps%usepaw==1) then
    call c_output_info_test(nfftf,env_den,cluster_den,ref_rhor,ucvol,stopTol,psps%usepaw,myicount,invt_quit,in_paw)
  else
    call c_output_info_test(nfftf,env_den,cluster_den,ref_rhor,ucvol,stopTol,psps%usepaw,myicount,invt_quit)
  endif
  if (invt_quit>0) then 
    if (i_am_master) then 
      ! write Status file to signal env to exit
      filename = "cluster.status"
      CALL checkfile(filename)  ! if the file exists, wait until it is removed
      open(unit=199,file='cluster.status',action='write',iostat=filestatus)
      if ( filestatus/=0 ) then
        CALL c_wrtlog("Error in creating cluster.status file, STOP")
        CALL leave_new('COLL')
      endif
      write(199,'(i)')-1
      close(199)
    endif
    call xbarrier_mpi(spaceworld,mpi_ier)
    exit
  endif
  
  !===============
  ! get G
  cgG_fg=ref_rhor-(cluster_den+env_den)                 !assume spin number = 1
  if (penlambda>0) then
    call laplacian(gprimd,mpi_enreg,nfftf,1,ngfftf,dtset%paral_kgb,rdfuncr=extpot,laplacerdfuncr=lapvfg)
    cgG_fg=cgG_fg-2._DP*penLambda*lapvfg(:,1)
    write(message,'(a,ES12.4,a,ES12.4)') & 
     '(with penalty) max/min cgG_fg :',maxval(cgG_fg),"  ",minval(cgG_fg)
    call c_wrtlog(message)
  endif
  ! transform from fine --> coarse grid
  if (psps%usepaw==1) then
    call coarse_to_fine(dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),nfft,cgG, & 
                  ngfftf(1),ngfftf(2),ngfftf(3),nfftf,cgG_fg,-1)
  else
    cgG=cgG_fg
  endif
  write (message,'(a,ES12.4,a,ES12.4,a)') & 
     '(coarse grid)  max/min cgG    :',maxval(cgG),"  ",minval(cgG),' (splined)'
  call c_wrtlog(message)

  !===============
  ! get W 
  call c_wval(cgW,nfftf,etotal_cluster,etotal_env,extpot, & 
              ref_rhor,cluster_den,env_den,ucvol,psps%usepaw);
  penfun = 0.d0
  if (penlambda>0) Then
    call laplacian(gprimd,mpi_enreg,nfftf,1,ngfftf,dtset%paral_kgb,rdfuncr=extpot,laplacerdfuncr=lapvfg)
    penfun=-penLambda*dot_product(extpot,lapvfg(:,1))*ucvol/dble(nfftf)
    call c_wrtlog('W: added penalty function.')
  endif
  cgW = cgW + penfun
  write(message,'(a,I4,a,g22.14,a,ES12.4)') & 
    '(gstate) iter:',myiCount,' cgW(with penalty)=',cgW,' penalty=', penfun
  call c_wrtlog(message)
 
100  continue  
  !----------------------------------------
  ! Optimizer
  !----------------------------------------
  if (optmethod==1 .or. optmethod==2) then 
    !========================================
    ! conjugate gradient 
    !========================================

    ! on the start of cg
    if (cg_stat==1) then
      cg_stat    = 0
      cg_extpot0 = extpot_c ! backup current extpot
      cgG_old    = cgG      ! backup old gradient
      cgD        = -cgG     ! d0=-g0
      stp        = 1.0d0
      task       = 'START'
    endif

    grad_stp = dot_product(cgG,cgD)  ! dW/d(stp)
    !---------------
    ! line search 
    if (task.eq.'START') then 
      write(message,'(3(a,ES16.8),a)')'(gstate) dcsrch W=',cgW,' stp=',0.0d0,' grad_stp=',grad_stp, ' <--- NEW LINE SEARCH'
    else 
      write(message,'(3(a,ES16.8))')'(gstate) dcsrch W=',cgW,' stp=',stp,' grad_stp=',grad_stp
    endif
    call c_wrtlog(message)
    call dcsrch(cgW,grad_stp,stp,ftol,gtol,xtol,stpmin,stpmax,task,isave,dsave )

    if (task(1:4) .eq. 'CONV') then
      call c_wrtlog('(gstate) dcsrch converged, prepare next cg direction.')

      ! Polak-Ribiere
      if (optmethod==1) then
        yk = cgG-cgG_old
        cg_beta = dot_product(cgG,yk)/(dot_product(cgG_old,cgG_old)) 
        cg_beta = max(cg_beta,0.d0)
      endif

      ! Hager and Zhang
      if (optmethod==2) then
        yk     = cgG-cgG_old
        dTy    = dot_product(cgD,yk)        ! use old D_k here
        norm_y = sqrt(dot_product(yk,yk))
        cg_beta= 1.d0/dTy * dot_product(yk-2.d0*cgD*norm_y**2/dTy,cgG)
        
        norm_dk= sqrt(dot_product(cgD,cgD))
        norm_gk= sqrt(dot_product(cgG_old,cgG_old))
        cg_beta= max(cg_beta,-1.d0/(norm_dk*min(0.d0,norm_gk)))
      endif

      cg_extpot0 = extpot_c            ! backup current extpot
      cgG_old    = cgG                 ! backup old gradient
      cgD        = -cgG + cg_beta*cgD  ! obtain new cg direction: D_k+1
      stp        = stp*2.d0            ! use this stp as initial guess for next line search
      task       = 'START'

      ! exceed maximum loops?
      if (myicount+1 >= maxcount) then
        call c_wrtlog('max cg iter reached. write out density and potential.')
        call write_status( -1, 'CLUSTER' )
        exit
      endif

      ! dump informations every out_step bfgs steps
      if ( myicount>=target_step .and. i_am_master) then  ! only master writes
        target_step = out_step + target_step
        if (psps%usepaw==1) then
          ! coarse to fine grid transform
          call coarse_to_fine(dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),nfft,extpot_c, & 
                         ngfftf(1),ngfftf(2),ngfftf(3),nfftf,extpot,+1)
        endif
        call dump_extpot(nfftf,extpot,myicount)
      endif 

      !goto dcsrch() again to initiate a new line search in new CG direction
      call xbarrier_mpi(spaceworld,mpi_ier)
      goto 100

    else if (task(1:4) .eq. 'WARN') then
      call c_wrtlog('(gstate) dcsrch => task=WARN, exit')
      write(message,*)'(gstate) detailed message=',task
      call c_wrtlog(message)
      exit

    else if (task(1:5) .eq. 'ERROR') then
      call c_wrtlog('(gstate) dcsrch ==> task=ERROR, exit')
      exit

    else if (task(1:2) .eq. 'FG') then
      extpot_c = cg_extpot0+stp*cgD
      write(message,'(a,ES12.4)')'(gstate) more line search. to try stp:',stp
      call c_wrtlog(message)
    endif

  else if (optmethod==11) then
    !=========================
    ! BFGS code
    !=========================
    ! for spin-polarized case, 
    ! if extpot=0 and mag_field=0, we might have degenerary (due to high symmetry)
    ! this would not happen for non-spin-polarzied case
    ! To avoid it, in the first iteration, we just take one step
    ! to remove this degeneracy
    if (dtset%nsppol==2 .and. myicount==1 .and. bloadVr==0) then 
      ! we use a step 0.05 here
      extpot = extpot + cgG*0.05d0 
      call c_wrtlog('')
      call c_wrtlog('gstate: nsp==2 and myicount==1, skip BFGS in the 1st iteration to avoid degeneracy')
      call c_wrtlog('')
      goto 88888  ! skip the following BFGS
    endif

    call setulb(nfft,mmax,extpot_c,BFGS_l,BFGS_u,BFGS_nbd,cgW,cgG,0._DP,0._DP, & 
              BFGS_wa,BFGS_iwa,task,BFGS_print,BFGS_csave,BFGS_lsave,BFGS_isave,BFGS_dsave)
 
    if (task(1:5) .EQ. 'NEW_X') then
      !-------------------
      ! task = NEW_X
      CALL c_wrtlog('BFGS: TASK=NEW_X')
 
      ! exceed maximum loops?
      if (myicount+1 > maxcount) then
        call c_wrtlog('Max BFGS iter reached. Write out density and potential.')
        call write_status( -1, 'CLUSTER' )
        exit
      endif

      !dump informations every out_step bfgs steps
      if ( myicount>=target_step .and. i_am_master) then  ! only master writes
        target_step = out_step + target_step
        if (psps%usepaw==1) then
          call coarse_to_fine(dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),dtset%nfft,extpot_c, & 
                              ngfftf(1),ngfftf(2),ngfftf(3),nfftf,extpot,+1)
        endif
        call dump_extpot(nfftf,extpot,myicount)
      endif 

      !Goto setulb() again
      call xbarrier_mpi(spaceworld,mpi_ier)
      goto 100 
  
    ELSE IF (task(1:2) .EQ. 'FG') then
      !task = FG
      call c_wrtlog('BFGS: TASK=FG'); 
      ! Jin this is to avoid local minimum....
      if (nline>=10) then
         task = 'START'
         nline    = 0
         goto 100
      endif
      nline = nline + 1
    ELSE 
      call c_wrtlog('What is Task?, I should not get to this point, EXIT anyway')
      write(message,*)'task ==>',task(1:10)
      call c_wrtlog(message)
      call write_status( -1, 'CLUSTER' )
      call c_wrtlog(task)
      call xbarrier_mpi(spaceworld,mpi_ier)
      EXIT
    ENDIF  ! BGFS if(task)

  else
    print *,'optmethod is not defined, optmethod=',optmethod
    call leave_new('COLL')
  endif
  
  !====================================
  ! update extpot, interpolate extpot 
  ! from extpot_c: transgrid()
!! DEBUG
!!  extpot_c = extpot_c - sum(extpot_c)/dble(nfft)
!! ENDDEBUG
  if (psps%usepaw==1) then
    write(message,'(a,es12.4,a,es12.4)')'(before interp) max/min extpot_c:',maxval(extpot_c),"  ",minval(extpot_c)
    call c_wrtlog(message)
    call coarse_to_fine(dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),dtset%nfft,extpot_c, & 
                         ngfftf(1),ngfftf(2),ngfftf(3),nfftf,extpot,+1)
    write(message,'(a,es12.4,a,es12.4)')'(after  interp) max/min extpot  :',maxval(extpot),"  ",minval(extpot)
    call c_wrtlog(message)
  else 
    extpot = extpot_c
  endif

88888  continue   

  !========= write extpot to disk for environment =========== 
  call xbarrier_mpi(spaceworld,mpi_ier)
  IF (whoAmI == 1) THEN 
    if (i_am_master) then   ! only master node writes
      filename = "new_extpot.tmp"
      open(unit=199,file=filename)
      close(199,status='delete')
      OPEN(unit=199,file=filename,action='write',iostat=filestatus,form='unformatted')
      IF ( filestatus/=0 ) THEN
        call c_wrtlog('error in writing new_extpot.tmp. stop')
        CALL leave_new('COLL')
      ENDIF
      WRITE(199)extpot(:)
      CLOSE(199)
      WRITE(message,'(a,2ES12.4)')  & 
        "dumpped new_extpot.tmp, max/min:", & 
        maxval(extpot),minval(extpot)
      CALL c_wrtlog(message)

      ! write cluster.status file
      call write_status( 1, 'CLUSTER' )
!      call c_print_time(199)
    endif ! master writes only
  ENDIF
  call xbarrier_mpi(spaceworld,mpi_ier)

END DO  
!---------------------------------------------------------
! chen: end of embedding potential loop 
!---------------------------------------------------------

99999 continue 

 IF (whoAmI==1 .and. i_am_master) THEN   ! master writes
   call c_wrtlog('Writing final extpot.')
!   open(unit=100, access="sequential", action="write", & 
!    file='extpot.final', form="formatted", status="replace");   
   open(unit=100, file='extpot.final')
!   do i=1,nfftf
!     write(100,*)extpot(i)
!   enddo 
   rewind(100)
   write(100,'(e22.14)') extpot(:)
   close(100)
   call c_wrtlog("extpot.final is dumped.")
 ENDIF
 call xbarrier_mpi(spaceworld,mpi_ier)

 call c_wrtlog('')
 call c_wrtlog('-------------------------------------------------------')
 write (message, *) 'Kinetic Energy = ', results_gs%energies%e_kinetic; call c_wrtlog(message)
 write (message, *) 'Hartree Energy = ', results_gs%energies%e_hartree; call c_wrtlog(message)
 write (message, *) 'XC      Energy = ', results_gs%energies%e_xc; call c_wrtlog(message)
 write (message, *) 'Total   Energy = ', results_gs%etotal; call c_wrtlog(message)
 call c_wrtlog('-------------------------------------------------------')

 call c_wrtlog('CODE END AT:')
 call date_and_time(VALUES=dateTime)    !today
 write(msg, '(A,i2,A,i2,A,i4,A,i2.2,A,i2.2,A,i2.2,A,i3.3)') "", &
             dateTime(3),"/",  dateTime(2),"/",  dateTime(1)," : ",             &
             dateTime(5), ":", dateTime(6), ":", dateTime(7), ".", dateTime(8)
 call c_wrtlog(msg);

 !================================
 ! deallocate arrays
 deallocate(BFGS_nbd,BFGS_l,BFGS_u,BFGS_wa,BFGS_iwa)
 deallocate(lapvfg,cgG,cgG_fg,cgG_fg3d)
 deallocate(ref_rhor,extpot,cluster_den,env_den)
 deallocate(hc_rhor_paw)
 if(psps%usepaw) then 
   deallocate(in_paw)
 endif
 if (optmethod==1 .or. optmethod==2) then
   deallocate(yk,cgD,cgG_old,cg_extpot0)
 endif


!===================================================================================
!============================== Chen: END OF MY CODE ===============================
!===================================================================================








     else
!      Conduct determination of U
       call pawuj_drive(acell,atindx,atindx1,cg,cpus,dtefield,dtfil,dtset,&
&       ecore,eigen,electronpositron,hdr,indsym,initialized,&
&       irrzon,kg,mpi_enreg,nattyp,nfftf,npwarr,nspinor,occ,&
&       paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
&       pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
&       scf_history,scf_in,symrec,taug,taur,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)
     end if

   else if (ionmov==1) then
!    Conduct molecular dynamics, with or without viscous damping

     call status(0,dtfil%filstat,iexit,level,'call move     ')
!    mpi_enreg%paralbd=0
     call move(acell,amass,atindx,atindx1,cg,cpus,dtefield,dtfil,dtset,&
&     ecore,eigen,electronpositron,hdr,indsym,initialized,irrzon,&
&     kg,mpi_enreg,&
&     nattyp,nfftf,npwarr,nspinor,occ,&
&     paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&     phnons,psps,pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
&     scf_history,scf_in,symrec,taug,taur,wffnew,wffnow,vel,wvl,xred,xred_old,ylm,ylmgr)

   else if (ionmov==2 .or. ionmov==3) then

!    Apply Broyden method for structural optimization, as
!    implemented by Jean-Christophe Charlier (May 1992)

     call status(0,dtfil%filstat,iexit,level,'call brdmin   ')
!    mpi_enreg%paralbd=0

     call brdmin(acell,atindx,atindx1,cg,cpus,dtefield,dtfil,dtset,&
&     ecore,eigen,electronpositron,hdr,indsym,initialized,irrzon,&
&     kg,mpi_enreg,mxfh,&
&     nattyp,nfftf,npwarr,nspinor,nxfh,occ,&
&     paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab, phnons,psps,&
&     pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprim,&
&     scf_history,scf_in,symrec,taug,taur,wffnew,wffnow,vel,wvl,xfhist,xred,xred_old,ylm,ylmgr)
!    call mkrdim(acell,rprim,rprimd)

   else if (ionmov==4 .or. ionmov==5) then

     do itime=1,ntime

       call status(itime,dtfil%filstat,iexit,level,'call scfcv(mv)')

       if(ionmov==4)then
         if(mod(itime,2)==1)then
           write(message, '(a,a,i3,a)' ) ch10,' STEP NUMBER ',itime,&
&           ' : OPTIMIZE ELECTRONS ------------------------------------'
         else
           write(message, '(a,a,i3,a)' ) ch10,' STEP NUMBER ',itime,&
&           ' : OPTIMIZE ELECTRONS AND IONS ---------------------------'
         end if
       else
         write(message, '(a,a,i3,a)' ) ch10,' STEP NUMBER ',itime,&
&         ' : SIMPLE RELAXATION -------------------------------------'
       end if
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')

!      In this case, iapp is simply itime
       iapp=itime
!      mpi_enreg%paralbd=0
       call scfcv_tmp(acell,atindx,atindx1,cg,cpus,dtefield,dtfil,dtpawuj,dtset,ecore,&
&       eigen,electronpositron,hdr,iapp,indsym,initialized,irrzon,kg,mpi_enreg,&
&       nattyp,ndtpawuj,nfftf,npwarr,nspinor,occ,paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&       phnons,psps,pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprimd,&
&       scf_history,scf_in,symrec,taug,taur,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)

       if(mod(itime,2)==1)then
!        When the SCF cycle dealt with electrons only,
!        check whether forces are below tolerance; if so, exit
!        from the itime loop
         itimexit=0 ; if(itime==ntime)itimexit=1
         call fconv(results_gs%fcart,dtset%iatfix,itimexit,itime,dtset%natom,&
&         ntime,0,1.0_dp,dtset%strtarget,results_gs%strten,dtset%tolmxf)
       end if
       if (itimexit/=0) exit

!      Check whether exiting was required by the user.
!      If found then beat a hasty exit from time steps
       if(dtset%chkexit==0) then
         openexit=0
       else
         openexit=1
       end if
       call chkexi(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)
       if (iexit/=0) then
         iexit=0   ! In order not to exit of dataset loop automatically
         exit
       end if

     end do

   else if ( (ionmov>=6 .and. ionmov<=9) .or. (ionmov>=12 .and. ionmov<=14) ) then

!    Molecular dynamics, using Verlet algorithm (ionmov=6)
!    or fake molecular dynamics for minimisation (ionmov=7)
!    or true molecular dynamics with Nose thermostat (ionmov=8)
!    or Langevin dynamics (ionmov=9) or Fei Zhang algorithm (ionmov=12)
!    or Molecular dynamics using the Runge-Kutta-Nyström methods SRKNa14
!    parametrized by Blanes and Moans (ionmov14)

     call status(0,dtfil%filstat,iexit,level,'call moldyn   ')

     call moldyn(acell,amass,atindx,atindx1,cg,cpus,dtefield,dtfil,&
&     dtset,ecore,eigen,electronpositron,hdr,indsym,initialized,&
&     irrzon,kg,mpi_enreg,mxfh,&
&     nattyp,nfftf,npwarr,nspinor,nxfh,occ,&
&     paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&     phnons,psps,pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprim,&
&     scf_history,scf_in,symrec,taug,taur,wffnew,wffnow,vel,wvl,xfhist,xred,xred_old,ylm,ylmgr)
!    call mkrdim(acell,rprim,rprimd)

   else if (ionmov == 10) then

     call delocint(acell,atindx,atindx1,cg,cpus,dtefield,dtfil,&
&     dtset,ecore,eigen,electronpositron,hdr,indsym,initialized,irrzon,&
&     kg,mpi_enreg,mxfh,&
&     nattyp,nfftf,npwarr,nspinor,nxfh,occ,&
&     paw_dmft,pawang,pawfgr,pawrad,pawrhoij,pawtab,phnons,psps,&
&     pwind,pwind_alloc,pwnsfac,rec_set,resid,results_gs,rhog,rhor,rprim,&
&     scf_history,scf_in,symrec,taug,taur,wffnew,wffnow,vel,wvl,xfhist,xred,xred_old,ylm,ylmgr)
!    call mkrdim(acell,rprim,rprimd)

   else if (ionmov == 20) then
!    Ground state call.
     iapp = 0
!    Ionic positions relaxation using DIIS. This algorithm is fast
!    and converge to the nearest singular point (where gradient vanishes).
!    This is a good algorithm to precisely tune saddle-points.
     call diisRelax(acell, atindx, atindx1, cg, cpus, dtefield, &
&     dtfil, dtset, ecore, eigen, electronpositron, hdr, iapp, indsym, initialized, &
&     irrzon, kg, mpi_enreg, nattyp, nfftf, npwarr, nspinor, occ, paw_dmft, pawang, &
&     pawfgr, pawrad, pawrhoij, pawtab, phnons, psps, pwind, pwind_alloc, pwnsfac, &
&     rec_set,resid, results_gs, rhog, rhor, rprimd, scf_history, scf_in, symrec, taug,taur,&
&     wffnew, wffnow, wvl, xred, xred_old, ylm, ylmgr)

   else if (ionmov == 30) then

     call scphon(acell, amass, atindx, atindx1, cg, cpus, dtefield, &
&     dtfil, dtset, ecore, eigen, electronpositron, hdr, indsym, initialized, &
&     irrzon, kg, mpi_enreg, nattyp, nfftf, npwarr, nspinor, occ, pawang, &
&     paw_dmft,pawfgr, pawrad, pawrhoij, pawtab, phnons, psps, pwind, pwind_alloc,pwnsfac, &
&     rec_set, resid, results_gs, rhog, rhor, rprimd, scf_history, scf_in, symrec, taug,taur, &
&     wffnew, wffnow, wvl, xred, xred_old, ylm, ylmgr)

   else ! Not an allowed option
     write(message, '(a,i12,a,a)' )&
&     ' Disallowed value for ionmov=',ionmov,ch10,&
&     ' Allowed values are 0 to 5.'
     MSG_BUG(message)
   end if

   call timab(35,2,tsec)

!  #####################################################################

!  End of the check of hasty exit
 end if

 call timab(36,1,tsec)

 write(message, '(80a,a,a,a,a)' ) ('=',mu=1,80),ch10,ch10,&
& ' ----iterations are completed or convergence reached----',ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

!Close the unneeded temporary data files, if any.
!Other files are closed in clnup1.
 if (dtset%mkmem==0) then
   close (unit=dtfil%unkg,status='delete')
   if (psps%useylm==1) close (unit=dtfil%unylm,status='delete')
   if (psps%usepaw==1) close (unit=dtfil%unpaw,status='delete')
   call WffDelete(wffnew,ierr)
 end if

!Will be put here later.
!!$ ! WVL - maybe compute the tail corrections to energy
!!$ if (dtset%tl_radius > real(0, dp)) then
!!$    ! Store xcart for each atom
!!$    allocate(xcart(3, dtset%natom))
!!$    call xredxcart(dtset%natom, 1, rprimd, xcart, xred)
!!$    ! Use the tails to improve energy precision.
!!$    call wvl_tail_corrections(dtset, results_gs%energies, results_gs%etotal, &
!!$         & mpi_enreg, occ, psps, vtrial, wvl, xcart)
!!$    deallocate(xcart)
!!$ end if

!Update the header, before using it
 call hdr_update(bantot,results_gs%etotal,results_gs%energies%e_fermie,hdr,dtset%natom,&
& results_gs%residm,rprimd,occ,mpi_enreg,pawrhoij,psps%usepaw,xred)

!if(dtset%tfkinfunc/=2) then
 call status(0,dtfil%filstat,iexit,level,'call outwf    ')

 if(dtset%nqpt==0)filnam=dtfil%fnameabo_wfk
 if(dtset%nqpt==1)filnam=dtfil%fnameabo_wfq
 call outwf(cg,dtset,eigen,filnam,hdr,kg,dtset%kptns,&
& dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,mxfh,dtset%natom,dtset%nband,&
& dtset%nkpt,npwarr,nspinor,dtset%nsppol,dtset%nstep,nxfh,&
& occ,resid,response,dtfil%unwff2,wffnow,wvl%wfs,xfhist)
!end if

 if(dtset%prtwf==2)then
   call outqmc(cg,dtset,eigen,gprimd,hdr,kg,dtset%mband,&
&   dtset%mkmem,mpi_enreg,dtset%mpw,dtset%nkpt,npwarr,&
&   nspinor,dtset%nsppol,occ,psps,results_gs)
 end if

 call status(0,dtfil%filstat,iexit,level,'call clnup1   ')

 call clnup1(acell,dtset%dosdeltae,dtset,eigen,dtset%enunit,&
& results_gs%energies%e_fermie,dtfil%fnameabo_dos,dtfil%fnameabo_eig,&
& results_gs%fred,dtset%iatfix,iscf,dtset%kptns,dtset%kptopt,dtset%mband,&
& mpi_enreg,dtset%natom,dtset%nband,nfftf,ngfftf,&
& dtset%nkpt,dtset%nspden,nspinor,dtset%nsppol,dtset%nstep,occ,dtset%occopt,&
& dtset%prtdos,dtset%prteig,dtset%optforces,dtset%prtstm,prtvol,resid,rhor,&
& rprimd,dtset%tphysel,dtset%tsmear,results_gs%vxcavg,dtset%wtk,xred)

 if ( (iscf>0 .or. iscf==-3) .and. dtset%prtstm==0) then
   call status(0,dtfil%filstat,iexit,level,'call prtene   ')
   call prtene(dtset,results_gs%energies,ab_out,psps%usepaw)
 end if

!Open the formatted derivative database file, and write the
!preliminary information
!In the // case, only one processor writes the energy and
!the gradients to the DDB

 if ((mpi_enreg%me==0).and.((iscf > 0).or.&
& (dtset%berryopt == -1).or.(dtset%berryopt) == -3)) then

   call status(0,dtfil%filstat,iexit,level,'call ioddb8_ou')
   vrsddb=010929
   dscrpt=' Note : temporary (transfer) database '
   ddbnm=trim(dtfil%filnam_ds(4))//'_DDB'
!  tolwfr must be initialized here, but it is a dummy value
   tolwfr=1.0_dp
   call ioddb8_out (dscrpt,ddbnm,dtset%natom,dtset%mband,&
&   dtset%nkpt,dtset%nsym,psps%ntypat,dtfil%unddb,vrsddb,&
&   acell,dtset%amu,dtset%dilatmx,dtset%ecut,dtset%ecutsm,&
&   dtset%intxc,dtset%iscf,dtset%ixc,dtset%kpt,dtset%kptnrm,&
&   dtset%natom,dtset%nband,ngfft,dtset%nkpt,dtset%nspden,nspinor,&
&   dtset%nsppol,dtset%nsym,psps%ntypat,occ,dtset%occopt,&
&   rprim,dtset%sciss,dtset%spinat,dtset%symafm,dtset%symrel,&
&   dtset%tnons,tolwfr,dtset%tphysel,dtset%tsmear,&
&   dtset%typat,dtset%wtk,xred,psps%ziontypat,dtset%znucl)

   if (iscf > 0) then
     nblok = 2          ! 1st blok = energy, 2nd blok = gradients
   else
     nblok = 1
   end if
   fullinit = 0 ; choice=2
   call psddb8 (choice,psps%dimekb,psps%ekb,fullinit,psps%indlmn,&
&   psps%lmnmax,psps%lnmax,nblok,&
&   psps%ntypat,dtfil%unddb,psps%pspso,psps%usepaw,psps%useylm,vrsddb)

   mpert = dtset%natom + 6   ; msize = 3*mpert
   allocate(blkflg(msize),blkval(2,msize))

   blkflg(:) = 0       ; blkval(:,:) = zero
   blkqpt(:) = zero    ; blknrm(:) = one

!  Write total energy to the DDB
   if (iscf > 0) then
     blktyp = 0
     blkval(1,1) = results_gs%etotal
     blkflg(1) = 1
     call blok8(blkflg,blknrm,blkqpt,blktyp,blkval,choice,dtset%mband,&
&     mpert,msize,dtset%nkpt,dtfil%unddb)
   end if

!  Write gradients to the DDB
   blktyp = 4
   blkflg(:) = 0       ; blkval(:,:) = zero
   indx = 0
   if (iscf > 0) then
     do iatom = 1, dtset%natom
       do idir = 1, 3
         indx = indx + 1
         blkflg(indx) = 1
         blkval(1,indx) = results_gs%fred(idir,iatom)
       end do
     end do
   end if

   indx = 3*dtset%natom + 3
   if ((abs(dtset%berryopt) == 1).or.(abs(dtset%berryopt) == 3)) then
     do idir = 1, 3
       indx = indx + 1
       if (dtset%rfdir(idir) == 1) then
         blkflg(indx) = 1
         blkval(1,indx) = results_gs%pel(idir)
       end if
     end do
   end if

   indx = 3*dtset%natom + 6
   if (iscf > 0) then
     blkflg(indx+1:indx+6) = 1
     blkval(1,indx+1:indx+6) = results_gs%strten(1:6)
   end if

   call blok8(blkflg,blknrm,blkqpt,blktyp,blkval,choice,dtset%mband,&
&   mpert,msize,dtset%nkpt,dtfil%unddb)

   deallocate(blkflg,blkval)

!  Close DDB
   close(dtfil%unddb)

 end if

 if (dtset%nstep>0 .and. dtset%prtstm==0 .and. dtset%positron/=1) then
   call status(0,dtfil%filstat,iexit,level,'call clnup2   ')
   call clnup2(psps%n1xccc,results_gs%fred,results_gs%gresid,&
&   results_gs%grewtn,&
&   results_gs%grxc,iscf,dtset%natom,dtset%optforces,dtset%optstress,prtvol,start,&
&   results_gs%strten,results_gs%synlgr,xred)
 end if

 call status(0,dtfil%filstat,iexit,level,'deallocate    ')

!Deallocate arrays
 deallocate(amass,atindx,atindx1,cg,eigen,indsym)
 deallocate(irrzon,npwarr,nattyp,phnons,resid)
 deallocate(rhog,rhor,start,symrec,taug,taur,xfhist,xred_old)
 deallocate(pawfgr%fintocoa,pawfgr%coatofin)
 if(dtset%tfkinfunc/=2) deallocate(ylm,ylmgr)

 if (psps%usepaw==1) then
   do iatom=1,dtset%natom
     deallocate(pawrhoij(iatom)%rhoijp,pawrhoij(iatom)%rhoijselect)
   end do
 end if
 deallocate(pawrhoij) ! this is allocated to 0 length in non PAW case

 if (use_scf_history) then
   if (psps%usepaw==1) then
     do ii=1,scf_history%history_size
       if (scf_history%hindex(ii)>0) then
         ixx=scf_history%hindex(ii)
         do iatom=1,dtset%natom
           deallocate(scf_history%pawrhoij(iatom,ixx)%rhoijselect,&
&           scf_history%pawrhoij(iatom,ixx)%rhoijp)
         end do
       end if
     end do
     deallocate(scf_history%pawrhoij)
   end if
   deallocate(scf_history%hindex,scf_history%deltarhor,scf_history%xreddiff,&
&   scf_history%atmrho_last)
   if (dtset%iextrapwf==1) then
     deallocate(scf_history%cg)
     if (psps%usepaw==1) then
       do ii=1,scf_history%history_size
         call cprj_free(scf_history%cprj(:,:,ii))
       end do
       deallocate(scf_history%cprj)
     end if
   end if
 end if
!RShaltaf: Changed to include SBC
 if (dtset%icoulomb > 0) then
!  Ask to deallocate the kernel part of Poisson's solver
!  Arguments are dummy ones since iaction == 0.
   call psolver_kernel(dtset, 0, kernel_dummy, mpi_enreg, rprimd)
 end if

!PAW+U
 if (dtset%usepawu>0.or.dtset%useexexch>0) then
   do itypat=1,ntypat
     if((dtset%lpawu(itypat)/=-1).or.(dtset%lexexch(itypat)/=-1)) deallocate(pawtab(itypat)%lnproju)
     if((dtset%lpawu(itypat)/=-1).or.(dtset%lexexch(itypat)/=-1)) deallocate(pawtab(itypat)%phiphjint)
     if((dtset%lpawu(itypat)/=-1).or.(dtset%lexexch(itypat)/=-1)) deallocate(pawtab(itypat)%ph0phiint)
   end do
 end if

!PAW+DMFT
!write(6,*) "before destroy_dmft", paw_dmft%use_dmft
 call destroy_sc_dmft(paw_dmft)

!Destroy electronpositron datastructure
 call destroy_electronpositron(electronpositron)

!Deallocating the basis set.
 if (dtset%usewvl == 1) then
   call wvl_free_type_wfs(wvl%wfs)
   call wvl_free_type_proj(wvl%projectors)
 else
   if(dtset%tfkinfunc /=2 )then
     deallocate(kg)
   end if
 end if

 if ((dtset%berryopt < 0).or.(dtset%berryopt == 4)) then
   deallocate(pwind,pwnsfac)
   deallocate(dtefield%ikpt_dk,dtefield%idxkstr)
   deallocate(dtefield%sflag,dtefield%cgindex,dtefield%kgindex)
   deallocate(dtefield%fkptns,dtefield%indkk_f2ibz,dtefield%i2fbz)
   if (mpi_enreg%paral_compil_kpt == 1) then
     deallocate(mpi_enreg%kptdstrb)
     if (dtset%berryopt == 4) then
       deallocate(mpi_enreg%kptdstrbi,dtefield%cgqindex,dtefield%nneigh)
     end if
   end if
 else
   deallocate(pwind,pwnsfac)
 end if


 if (dtset%berryopt == 4) deallocate(dtefield%smat)

!deallocate Recursion
 if (dtset%userec == 1)  call CleanRec(rec_set)
!Clean the header
 call hdr_clean(hdr)

 if (mpi_enreg%me == 0 .and. dtset%prtxml == 1) then
!  The dataset given in argument has been treated, then
!  we output its variables.
!  call outvarsXML()
!  gstate() will handle a dataset, so we output the dataSet markup.
   write(ab_xml_out, "(A)") '  </dataSet>'
 end if

 if (dtset%usewvl == 0) then
!  Clean the MPI informations
   call clnmpi_bandfft(mpi_enreg)
   call clnmpi_gs(mpi_enreg)
   call clnmpi_fft(mpi_enreg)
   call clnmpi_atom(mpi_enreg)
 else
!  Clean the wavelet part.
   if (associated(mpi_enreg%nscatterarr)) then
     deallocate(mpi_enreg%nscatterarr)
   end if
   if (associated(mpi_enreg%ngatherarr)) then
     deallocate(mpi_enreg%ngatherarr)
   end if
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(36,2,tsec)
 call timab(32,2,tsec)

 DBG_EXIT("COLL")

end subroutine gstate


!=========================== chen =============
! dump extpot to file 
!
subroutine dump_extpot(nfftf,extpot,iter)
 implicit none
 real(8) :: extpot(nfftf)
 integer :: nfftf,iter,i
 character (len=200) :: message

 write(message,'(I3)') iter 
 open (unit=100, access="sequential", action="write", & 
!    file='extpot.'//adjustl(trim(message)), form="formatted", status="replace");   
     file='extpot.tmp');
 do i=1,nfftf
   write (100,*)extpot(i)
 end do 
 close (100)
 call c_wrtlog("extpot.{#} is dumped")
end subroutine  dump_extpot

!====================================
! output status of inverting KS eqn.
! and test if quit criterion
!
subroutine c_output_info_test(nfftf,env_den,cluster_den,ref_rhor,ucvol,stopTol,usepaw,myicount,invt_quit, &
  in_paw)  ! optional

  implicit none
  integer, intent(in)        :: nfftf,usepaw,myicount
  integer, intent(out)       :: invt_quit
  real(kind=8), intent(in)   :: env_den(nfftf),cluster_den(nfftf),ref_rhor(nfftf),ucvol,stopTol
  integer,optional,intent(in):: in_paw(nfftf)

  !
  ! local vars
  !
  character(len=500)::message
  integer :: ii
  real(kind=8) ::  & 
   c_ratio,c_ratio_max,c_ratio_min, & 
   c_rho_diff_paw,c_rho_diff_paw_max,c_rho_diff_paw_min,&
   c_rho_diff_outpaw,c_rho_diff_outpaw_max,c_rho_diff_outpaw_min
 
  !
  ! output information for current step
  !
  if (usepaw==0) then
    WRITE(message,'(a,i3,3(a,ES12.4))')"step:", myicount, & 
    " rho_diff_integral = ",SUM(ABS(env_den+cluster_den-ref_rhor))*ucvol/dble(nfftf), & 
    " max_diff=", maxval(ABS(env_den+cluster_den-ref_rhor)),  & 
    " stopTol=",stopTol
    CALL c_wrtlog(message)
  else if (usepaw==1) then
    c_ratio_max=-1e4
    c_ratio_min=1e4
    c_rho_diff_paw_max = -1e4
    c_rho_diff_paw_min = 1e4
    c_rho_diff_outpaw_max = -1e4
    c_rho_diff_outpaw_min = 1e4
    do ii=1,nfftf
      if (in_paw(ii)>0) then
        ! inside PAW sphere
        c_ratio=abs( (env_den(ii)+cluster_den(ii)-ref_rhor(ii))/ref_rhor(ii) )
        if (c_ratio_max<c_ratio) c_ratio_max=c_ratio
        if (c_ratio_min>c_ratio) c_ratio_min=c_ratio
        c_rho_diff_paw = abs(env_den(ii)+cluster_den(ii)-ref_rhor(ii))
        if (c_rho_diff_paw_max<c_rho_diff_paw) c_rho_diff_paw_max=c_rho_diff_paw
        if (c_rho_diff_paw_min>c_rho_diff_paw) c_rho_diff_paw_min=c_rho_diff_paw
      else
        ! outside PAW sphere
        c_rho_diff_outpaw = abs(env_den(ii)+cluster_den(ii)-ref_rhor(ii))
        if (c_rho_diff_outpaw_max<c_rho_diff_outpaw) c_rho_diff_outpaw_max=c_rho_diff_outpaw
        if (c_rho_diff_outpaw_min>c_rho_diff_outpaw) c_rho_diff_outpaw_min=c_rho_diff_outpaw
      endif
    enddo
    write(message,'(a,I4,a,es12.4,a,es12.4,a,a,es12.4,a)') & 
      '(gstate) step:', myicount,& 
      ' max_den_diff:',c_rho_diff_paw_max,' (in PAW)', & 
      c_rho_diff_outpaw_max, ' (out of PAW)', & 
      ' max_diff_ratio:',c_ratio_max,' (in PAW)'
    call c_wrtlog(message)
  else
    print *,'(c_output_info_test) undefined psps%usepaw:',usepaw,' code stop'
    call leave_new('COLL')
  endif

  !-------------------------
  ! Test density convergence
  !
  invt_quit = -1
  if (maxval(ABS(env_den+cluster_den-ref_rhor))<stopTol) then
    write(message,'(a,ES12.4,a)')"max density difference < ", stopTol, ". Done!"
    CALL c_wrtLog(message)
    invt_quit = 1
  endif

  return
end subroutine c_output_info_test

!-------------------------------------------------------------------
! do coarse grid ==> fine grid interpolation in real space
! using interpolation , using periodic boundary condition
!
!  dc: array on coarse grid ==> dc(nfft)
!  df: array on fine gird   ==> df(nfftf)
!  nfft  = nc1*nc2*nc3
!  nfftf = nf1*nf2*nf3
!  
!  flag = 1 :   coarse -> fine 
!       = -1:   fine   -> coarse
!
subroutine coarse_to_fine(nc1,nc2,nc3,nfft,dc, & 
                          nf1,nf2,nf3,nfftf,df,flag)
  ! for 3D spline
  use EZspline_obj
  use EZspline

  implicit none 
  
  integer,intent(in) :: flag
  integer,intent(in) :: nfft,nfftf,nc1,nc2,nc3,nf1,nf2,nf3
  real(kind=8)  :: dc(nfft)     ! data on coarse grid
  real(kind=8)  :: df(nfftf)    ! data on fine grid
  
  !======= local variables =========
  integer :: n_old1, n_old2, n_old3
  integer :: n_new1, n_new2, n_new3
  integer :: i,j,k,idx,ix,iy,iz
  real*8,parameter :: two_pi = 6.28318531d0
  real*8,allocatable :: work3d(:,:,:)
  character (len=500) :: message

  !======= for EZspline =====
  type(EZspline3_r8) :: f_spl ! 3-d object/real*8
  integer spn1, spn2, spn3, sp_ier, bcs1(2), bcs2(2), bcs3(2)
  real*8 :: x_int, y_int, z_int, f_int
  
  !========= function begins =======
  
  if (nc1*nc2*nc3/=nfft) then
    print *,'nc1*nc2*nc3/=nfft, code stop'
    call leave_new('COLL')
    return
  endif
  if (nf1*nf2*nf3/=nfftf) then
    print *,'nf1*nf2*nf3/=nfftf, code stop'
    call leave_new('COLL')
    return
  endif

  write(message,'(a,I3)')'get in coarse_to_fine(),flag=',flag
  call wrtout(6,message,'COLL')

  if (flag==1) then
    ! from coarse grid ===> fine grid
    n_old1=nc1; n_old2=nc2;  n_old3=nc3;
    n_new1=nf1; n_new2=nf2;  n_new3=nf3;
    write(message,'(a,2es14.4)')'   min/max data_coarse =',minval(dc),maxval(dc)
    call wrtout(6,message,'COLL')
  elseif(flag==-1) then
    ! from fine grid ===> coarse grid
    n_old1=nf1; n_old2=nf2;  n_old3=nf3;
    n_new1=nc1; n_new2=nc2;  n_new3=nc3;
    write(message,'(a,2es14.4)')'   min/max data_fine   =',minval(df),maxval(df)
    call wrtout(6,message,'COLL')
  else
    print *,'flag is not defined, flag=',flag,& 
    ' in coarse_to_fine subroutine, code stop'
    call leave_new('COLL')
    return
  endif

  spn1=n_old1+1
  spn2=n_old2+1
  spn3=n_old3+1
  bcs1=-1  ! boundary condition, -1 ==> periodic 
  bcs2=-1
  bcs3=-1

  ! prepare 3d data (coarse grid)
  allocate(work3d(spn1,spn2,spn3))

  do k=1,n_old3+1
    do j=1,n_old2+1
      do i=1,n_old1+1
        ix=i;iy=j;iz=k
        if (i==n_old1+1) ix=1
        if (j==n_old2+1) iy=1
        if (k==n_old3+1) iz=1
        idx=ix+(iy-1)*n_old1+(iz-1)*n_old2*n_old1
        if (flag==1)  work3d(i,j,k)=dc(idx)  ! coarse grid
        if (flag==-1) work3d(i,j,k)=df(idx)  ! fine grid
      enddo
    enddo
  enddo


  ! coarse grid ===> fine grid , by 3d interpolation
  call EZspline_init(f_spl,spn1,spn2,spn3,bcs1,bcs2,bcs3,sp_ier)
  call EZspline_error(sp_ier)
  f_spl%isHermite = 1  ! use Akima Hermite is required  
  call EZspline_setup(f_spl,work3d,sp_ier) ! set up coefficients
  call EZspline_error(sp_ier)

  do k=1,n_new3
    do j=1,n_new2
      do i=1,n_new1
        x_int = two_pi*dble(i-1)/dble(n_new1)  
        y_int = two_pi*dble(j-1)/dble(n_new2)  
        z_int = two_pi*dble(k-1)/dble(n_new3)  
        call EZspline_interp(f_spl,x_int,y_int,z_int,f_int,sp_ier)
        call EZspline_error(sp_ier)
        idx=i+(j-1)*n_new1+(k-1)*n_new2*n_new1
        if (flag==1) df(idx)=f_int
        if (flag==-1)dc(idx)=f_int
      enddo
    enddo
  enddo

  call EZspline_free(f_spl,sp_ier)
  call EZspline_error(sp_ier)

  deallocate(work3d)

  if (flag==1) then 
    write(message,'(a,2es14.4,a)')'   min/max data_fine   =',minval(df),maxval(df), ' (interpolated)'
    call wrtout(6,message,'COLL')
  else
    write(message,'(a,2es14.4,a)')'   min/max data_coarse =',minval(dc),maxval(dc), ' (interpolated)'
    call wrtout(6,message,'COLL')
  endif

  write(message,'(a,I3)')'leave coarse_to_fine(),flag=',flag
  call wrtout(6,message,'COLL')

  return 
end subroutine coarse_to_fine

!===========================================
! write status file of cluster or env.
! these status files will communicate cluster or env
!
!   flag           : the code to send
!   cluster_or_env : 'COL' or 'ENV'
!
! by Chen Huang  May/30/2010
!
subroutine write_status( flag, cluster_or_env )

  implicit none
  
  character(len=200),intent(in) :: cluster_or_env
  character(len=200) :: filename
  integer :: flag, filestatus

  if (cluster_or_env(1:3)=='CLU') then 
    filename = "cluster.status"
  elseif (cluster_or_env(1:3) =='ENV') then 
    filename = "../cluster/env.status"
  else
    print *,'filename = ',filename,' not defined ,stop'
    call leave_new('COLL')
    return
  endif

  CALL checkfile(filename) 
  open(unit=199,file=filename,action='write',iostat=filestatus)
  if ( filestatus/=0 ) then
    print *, 'filename=',filename,' (write_status) error on open, <=src/95_drive/gstate.F90'
    call leave_new('COLL')
    return
  endif
  write(199,'(i,a)') flag, '  end'
  close(199)
  return
end subroutine write_status
!!***
