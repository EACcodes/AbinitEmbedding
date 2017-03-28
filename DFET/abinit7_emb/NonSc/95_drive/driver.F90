!{\src2tex{textfont=tt}}
!!****f* ABINIT/driver
!! NAME
!! driver
!!
!! FUNCTION
!! Driver for ground state, response function, susceptibility, screening
!! and sigma calculations. The present routine drives the following operations.
!! An outer loop allows computation related to different data sets.
!! For each data set, either a GS calculation, a RF calculation,
!! a SUS calculation, a SCR calculation or a SIGMA calculation is made.
!! In both cases, the input variables are transferred in the proper variables,
!! selected big arrays are allocated, then the gstate, respfn or suscep subroutines are called.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2010 ABINIT group (XG,MKV,MM,MT,FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! codvsn= code version
!! cpui=initial CPU time
!! dtsets(0:ndtset_alloc)=<type datasets_type>contains all input variables
!! filnam(5)=character strings giving file names
!! filstat=character strings giving name of status file
!! mpi_enreg=informations about MPI parallelization
!! ndtset=number of datasets
!! ndtset_alloc=number of datasets, corrected for allocation of at
!!               least one data set.
!! npsp=number of pseudopotentials
!! pspheads(npsp)=<type pspheader_type>all the important information from the
!!   pseudopotential file header, as well as the psp file name
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  Input/Output
!! results_out(0:ndtset_alloc)=<type results_out_type>contains the results
!!   needed for outvars, including evolving variables
!!   Default values are set up in the calling routine
!!
!! NOTES
!! The array filnam is used for the name of input and output files,
!! and roots for generic input, output or temporary files.
!! Pseudopotential file names are set in pspini and pspatm,
!! using another name. The name filstat will be needed beyond gstate to check
!! the appearance of the "exit" flag, to make a hasty exit, as well as
!! in order to output the status of the computation.
!!
!! TODO
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      chkdilatmx,chkexi,echo_xc_name,getdim_nloc,gstate,leave_new
!!      mkfilename,mkrdim,nonlinear,respfn,screening,sigma
!!      status,suscep,timab,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine driver(codvsn,cpui,dtsets,filnam,filstat,&
&                 mpi_enreg,ndtset,ndtset_alloc,npsp,pspheads,results_out)

 use defs_basis
 use defs_parameters
 use defs_datatypes
 use defs_abitypes
 use m_errors
#ifdef HAVE_ABI_TIMER
 use m_timer
#endif
#if defined HAVE_LIBXC
 use libxc_functionals
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_42_geometry
 use interfaces_53_abiutil
 use interfaces_56_xc
 use interfaces_57_iovars
 use interfaces_59_io_mpi
 use interfaces_65_psp
 use interfaces_66_paw
 use interfaces_77_suscep
 use interfaces_93_rdm
 use interfaces_95_drive, except_this_one => driver
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndtset,ndtset_alloc,npsp
 real(dp),intent(in) :: cpui
 character(len=6),intent(in) :: codvsn
 character(len=fnlen),intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 character(len=fnlen),intent(in) :: filnam(5)
 type(dataset_type),intent(in) :: dtsets(0:ndtset_alloc)
 type(pspheader_type),intent(in) :: pspheads(npsp)
 type(results_out_type),intent(inout) :: results_out(0:ndtset_alloc)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=2
 integer,save :: dimekb_old=-1,lmnmax_old=-1,lnmax_old=-1,mqgridff_old=0
 integer,save :: mqgridvl_old=0,ntypat_old=-1,paw_size_old=-1,usepaw_old=-1
 integer :: idtset,iexit,iget,ii,iimage,iimage_get
 integer :: ipsp,ireadwf,jdtset
 integer :: jdtset_status,lmnmax,lmnmaxso,lnmax,lnmaxso!,mband
 integer :: mdtset,mtypalch,mu,mpsang,mxnimage !,mdtsetmpsang,mgfft,mk1mem,mkmem,mkqmem,mpssoang,mpw
 integer :: nimage,n1xccc,openexit,paw_size,prtvol !,natom,nfft,nkpt,nspden,nsppol,nsym
 integer :: usepaw
 real(dp) :: etotal
 character(len=500) :: message
 logical :: converged
 type(dataset_type) :: dtset
 type(datafiles_type) :: dtfil
 type(pawang_type) :: pawang
 type(pseudopotential_type) :: psps
 type(vardims_type) :: abidims
!arrays
 integer :: mkmems(3)
 integer,allocatable :: jdtset_(:),npwtot(:)
 real(dp) :: acell(3),rprim(3,3),rprimd(3,3),rprimdget(3,3),strten(6),tsec(2)
 real(dp),allocatable :: acell_img(:,:),rprim_img(:,:,:)
 real(dp),allocatable :: fcart(:,:),fred(:,:)
 real(dp),allocatable :: fcart_img(:,:,:),fred_img(:,:,:)
 real(dp),allocatable :: etotal_img(:),strten_img(:,:)
 real(dp),allocatable :: miximage(:,:)
 real(dp),allocatable :: occ(:),vel(:,:),xcart(:,:),xred(:,:),xredget(:,:)
 real(dp),allocatable :: occ_img(:,:),vel_img(:,:,:),xred_img(:,:,:)
 type(pawrad_type),allocatable :: pawrad(:)
 type(pawtab_type),allocatable :: pawtab(:)

!******************************************************************

 DBG_ENTER("COLL")

 ireadwf=0

 call timab(100,1,tsec)
 call status(0,filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 prtvol=dtsets(1)%prtvol
 if(prtvol==-level)then
   write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,&
&   ' driver : enter , debug mode '
   call wrtout(std_out,message,'COLL')
 end if

 mdtset=99

 if(ndtset>mdtset)then
   write(message,'(a,i2,a,i5,a)')&
&   '  The maximal allowed ndtset is ',mdtset,' while the input value is ',ndtset,'.'
   MSG_BUG(message)
 end if

 mtypalch=dtsets(1)%ntypalch
 do ii=1,ndtset_alloc
   mtypalch=max(dtsets(ii)%ntypalch,mtypalch)
 end do

!Allocation of some arrays independent of the dataset
 allocate(psps%filpsp(npsp))
 allocate(psps%pspcod(npsp))
 allocate(psps%pspdat(npsp))
 allocate(psps%pspso(npsp))
 allocate(psps%pspxc(npsp))
 allocate(psps%title(npsp))
 allocate(psps%zionpsp(npsp))
 allocate(psps%znuclpsp(npsp))
 call psp2params_init(psps%gth_params, npsp)

 psps%filpsp(1:npsp)=pspheads(1:npsp)%filpsp
 psps%pspcod(1:npsp)=pspheads(1:npsp)%pspcod
 psps%pspdat(1:npsp)=pspheads(1:npsp)%pspdat
 psps%pspso(1:npsp)=pspheads(1:npsp)%pspso
 psps%pspxc(1:npsp)=pspheads(1:npsp)%pspxc
 psps%title(1:npsp)=pspheads(1:npsp)%title
 psps%zionpsp(1:npsp)=pspheads(1:npsp)%zionpsp
 psps%znuclpsp(1:npsp)=pspheads(1:npsp)%znuclpsp

 allocate(jdtset_(0:ndtset))
 if(ndtset/=0)then
   jdtset_(:)=dtsets(0:ndtset)%jdtset
 else
   jdtset_(0)=0
 end if

!*********************************************************************
!Big loop on datasets

!Do loop on idtset (allocate statements are present)
 do idtset=1,ndtset_alloc

#ifdef HAVE_ABI_TIMER
   call init_timer()
#endif

   jdtset=dtsets(idtset)%jdtset ; if(ndtset==0)jdtset=1

   if(ndtset>=2)then
     jdtset_status=jdtset
   else
     jdtset_status=0
   end if

   call status(jdtset_status,filstat,iexit,level,'loop jdtset   ')

   write(message,'(a,80a,a,a,i2,a,66a,a)') ch10,&
&   ('=',mu=1,80),ch10,&
&   '== DATASET ',jdtset,' ',('=',mu=1,66),ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'PERS')     ! PERS is choosen to make debugging easier

!  Determine here whether the calculation is PAW
   usepaw  =0
   if (pspheads(1)%pspcod==7) usepaw=1  ! If paw, all pspcod necessarily are 7 (see iofn2)

!  Copy input variables into a local dtset.
   call dtsetCopy(dtset, dtsets(idtset))

!  Set other values
   dtset%jdtset = jdtset
   dtset%ndtset = ndtset

!  Copy input values
   mkmems(1)        = dtset%mkmem
   mkmems(2)        = dtset%mkqmem
   mkmems(3)        = dtset%mk1mem
   psps%optnlxccc   = dtset%optnlxccc
   psps%mqgrid_ff   = dtset%mqgrid
   if (usepaw == 1) then
     psps%mqgrid_vl = dtset%mqgriddg
   else
     psps%mqgrid_vl = dtset%mqgrid
   end if

   mpi_enreg%paral_compil_respfn=dtset%paral_rf
   mpi_enreg%ngroup_respfn=dtset%ngroup_rf

   
   allocate(occ(dtset%mband*dtset%nkpt*dtset%nsppol))
   allocate(vel(3,dtset%natom) )
   allocate(xred(3,dtset%natom))

   nimage=dtset%nimage
   allocate(acell_img(3,nimage))
   allocate(occ_img(dtset%mband*dtset%nkpt*dtset%nsppol,nimage))
   allocate(rprim_img(3,3,nimage))
   allocate(vel_img(3,dtset%natom,nimage))
   allocate(xred_img(3,dtset%natom,nimage))

   acell_img(:,:)   = dtset%acell_orig(:,1:nimage)
   rprim_img(:,:,:) = dtset%rprim_orig(:,:,1:nimage)
   vel_img  (:,:,:) = dtset%vel_orig(:,1:dtset%natom,1:nimage)
   xred_img (:,:,:) = dtset%xred_orig(:,1:dtset%natom,1:nimage)

!  Note that occ is not supposed to depend on the image in the input file.
!  However, the results will depend on the image. And occ can be used to reinitialize the
!  next dataset, or the next timimage. 
   do iimage=1,nimage
     occ_img  (:,iimage)   = dtset%occ_orig(1:dtset%mband*dtset%nkpt*dtset%nsppol)
   end do

!  ****************************************************************************
!  Treat the file names (get variables)

   call status(jdtset_status,filstat,iexit,level,'filenames     ')

   call dtfil_init1(dtfil,dtset,filnam,filstat,idtset,jdtset_,mpi_enreg,ndtset)

!  ****************************************************************************
!  Treat other get variables

!  If multi dataset mode, and already the second dataset,
!  treatment of other get variables.
   if( ndtset>1 .and. idtset>1 )then

     call status(jdtset_status,filstat,iexit,level,'get variables ')

     mxnimage=maxval(dtsets(1:ndtset_alloc)%nimage)
     allocate(miximage(mxnimage,mxnimage))

     call find_getdtset(dtsets,dtset%getocc,'getocc',idtset,iget,miximage,mxnimage,ndtset_alloc)
     if(iget/=0)then
       do iimage=1,nimage
         occ_img(:,iimage)=zero
         do iimage_get=1,dtsets(iget)%nimage
           occ_img(:,iimage)=occ_img(:,iimage)+&
&           miximage(iimage,iimage_get)*results_out(iget)%occ(1:dtset%mband*dtset%nkpt*dtset%nsppol,iimage_get)
         end do
       end do
     end if

!    Getcell has to be treated BEFORE getxcart
!    since acell and rprim will be used
     call find_getdtset(dtsets,dtset%getcell,'getcell',idtset,iget,miximage,mxnimage,ndtset_alloc)
     if(iget/=0)then
       do iimage=1,nimage
         acell_img(:,iimage)=zero
         rprim_img(:,:,iimage)=zero
         do iimage_get=1,dtsets(iget)%nimage
           acell_img(:,iimage)=acell_img(:,iimage)+&
&           miximage(iimage,iimage_get)*results_out(iget)%acell(:,iimage_get)
           rprim_img(:,:,iimage)=rprim_img(:,:,iimage)+&
&           miximage(iimage,iimage_get)*results_out(iget)%rprim(:,:,iimage_get)
!          Check that the new acell and rprim are consistent with the input dilatmx
           call mkrdim(acell_img(:,iimage),rprim_img(:,:,iimage),rprimd)
           call chkdilatmx(dtset%dilatmx,rprimd,dtset%rprimd_orig(1:3,1:3,iimage))
         end do
       end do
     end if

     call find_getdtset(dtsets,dtset%getxred,'getxred',idtset,iget,miximage,mxnimage,ndtset_alloc)
     if(iget/=0)then
       do iimage=1,nimage
         xred_img(:,:,iimage)=zero
         do iimage_get=1,dtsets(iget)%nimage
           xred_img(:,:,iimage)=xred_img(:,:,iimage)+&
&           miximage(iimage,iimage_get)*results_out(iget)%xred(:,1:dtset%natom,iimage_get)
         end do
       end do
     end if

     call find_getdtset(dtsets,dtset%getxcart,'getxcart',idtset,iget,miximage,mxnimage,ndtset_alloc)
     if(iget/=0)then
       allocate(xcart(3,dtset%natom),xredget(3,dtset%natom))
       do iimage=1,nimage
         xred_img(:,:,iimage)=zero
         do iimage_get=1,dtsets(iget)%nimage
!          Compute xcart of the previous dataset
           rprimdget(:,:)=results_out(iget)%rprimd(:,:,iimage_get)
           xredget (:,:)=results_out(iget)%xred(:,1:dtset%natom,iimage_get)
           call xredxcart(dtset%natom,1,rprimdget,xcart,xredget)
!          xcart from previous dataset is computed. Now, produce xred for the new dataset,
!          with the new acell and rprim ...
           call mkrdim(acell_img(:,iimage),rprim_img(:,:,iimage),rprimd)
           call xredxcart(dtset%natom,-1,rprimd,xcart,xredget(:,:))
           xred_img(:,:,iimage)=xred_img(:,:,iimage)+&
&           miximage(iimage,iimage_get)*xredget(:,:)
         end do
       end do
       deallocate(xcart,xredget)
     end if

     call find_getdtset(dtsets,dtset%getvel,'getvel',idtset,iget,miximage,mxnimage,ndtset_alloc)
     if(iget/=0)then
       do iimage=1,nimage
         vel_img(:,:,iimage)=zero
         do iimage_get=1,dtsets(iget)%nimage
           vel_img(:,:,iimage)=vel_img(:,:,iimage)+&
&           miximage(iimage,iimage_get)*results_out(iget)%vel(:,1:dtset%natom,iimage_get)
         end do
       end do
     end if

     deallocate(miximage)

   end if

!  ****************************************************************************
!  Treat the pseudopotentials : initialize the psps variable

   call status(jdtset_status,filstat,iexit,level,'init psps     ')

!  Note that mpsang is the max of 1+lmax, with minimal value 1 (even for local psps, at present)
!  mpsang=max(maxval(pspheads(1:npsp)%lmax)+1,1) ! might not work with HP compiler
!  n1xccc=maxval(pspheads(1:npsp)%xccc)
   mpsang=1
   n1xccc=pspheads(1)%xccc
   do ii=1,npsp
     mpsang=max(pspheads(ii)%lmax+1,mpsang)
     n1xccc=max(pspheads(ii)%xccc,n1xccc)
   end do

!  Determine the maximum number of projectors, for the set of pseudo atom
   call getdim_nloc(lmnmax,lmnmaxso,lnmax,lnmaxso,dtset%mixalch,npsp,dtset%npspalch,&
&   dtset%ntypat,dtset%ntypalch,pspheads)

   psps%mpsang   = mpsang
   psps%mtypalch = mtypalch
   psps%npsp     = npsp
   psps%npspalch = dtset%npspalch
   psps%ntypat   = dtset%ntypat
   psps%ntypalch = dtset%ntypalch
   psps%ntyppure = dtset%ntyppure
   psps%n1xccc   = n1xccc

!  Set the flag for reciprocal space or real space calculations
   psps%vlspl_recipSpace = (dtset%icoulomb /= 1)
!  changed by RShaltaf
   psps%positron = dtset%positron
   psps%usepaw  =usepaw
   psps%useylm  =dtset%useylm

   allocate(psps%algalch(psps%ntypalch))
   allocate(psps%mixalch(psps%npspalch,psps%ntypalch))
   psps%algalch(1:psps%ntypalch)=dtset%algalch(1:psps%ntypalch)
   psps%mixalch(1:psps%npspalch,1:psps%ntypalch)=dtset%mixalch(1:psps%npspalch,1:psps%ntypalch)

!  Set mpspso and psps%pspso
!  Warning : mpspso might be different for each dataset.
   psps%mpspso=1
   do ipsp=1,dtset%npsp
     if(dtset%nspinor==1)then
       psps%pspso(ipsp)=0
     else
       if(dtset%so_psp(ipsp)/=1)then
         psps%pspso(ipsp)=dtset%so_psp(ipsp)
       else
         psps%pspso(ipsp)=pspheads(ipsp)%pspso
       end if
       if(psps%pspso(ipsp)/=0)psps%mpspso=2
     end if
!    Ideally the following line should not exist, but at present, the space has to be booked
     if(pspheads(ipsp)%pspso/=0)psps%mpspso=2
   end do

!  Set mpssoang, lmnmax, lnmax
   if(psps%mpspso==1)then
     psps%mpssoang=psps%mpsang
     psps%lmnmax  =lmnmax
     psps%lnmax   =lnmax
   else
     psps%mpssoang=2*psps%mpsang-1
     psps%lmnmax=lmnmaxso
     psps%lnmax=lnmaxso
   end if
   if (psps%useylm==0) then
     psps%lmnmax=psps%lnmax
   end if

!  Set dimekb
   if (psps%usepaw==0) then
     psps%dimekb=psps%lnmax
   else
     psps%dimekb=psps%lmnmax*(psps%lmnmax+1)/2
   end if

!  The following arrays are often not deallocated before the end of the dtset loop
!  and might keep their content from one dataset to the other,
!  if the conditions are fulfilled
   if(dimekb_old/=psps%dimekb .or. ntypat_old/=dtset%ntypat .or. &
&   usepaw_old/=psps%usepaw) then
     if(idtset/=1)deallocate(psps%ekb)
     allocate(psps%ekb(psps%dimekb,dtset%ntypat*(1-psps%usepaw)))
     dimekb_old=psps%dimekb
   end if
   if(lmnmax_old/=psps%lmnmax .or. ntypat_old/=dtset%ntypat)then
     if(idtset/=1)deallocate(psps%indlmn)
     allocate(psps%indlmn(6,psps%lmnmax,dtset%ntypat))
     lmnmax_old=psps%lmnmax
   end if
   if(mqgridff_old/=psps%mqgrid_ff .or. lnmax_old/=psps%lnmax .or. ntypat_old/=dtset%ntypat)then
     if(idtset/=1)deallocate(psps%ffspl,psps%qgrid_ff)
     allocate(psps%ffspl(psps%mqgrid_ff,2,psps%lnmax,dtset%ntypat))
     allocate(psps%qgrid_ff(psps%mqgrid_ff))
     mqgridff_old=psps%mqgrid_ff
     lnmax_old=psps%lnmax
   end if
   if(mqgridvl_old/=psps%mqgrid_vl .or. ntypat_old/=dtset%ntypat)then
     if(idtset/=1)deallocate(psps%vlspl,psps%qgrid_vl)
     if (idtset/=1 .and. .not.psps%vlspl_recipSpace) then
       deallocate(psps%dvlspl)
     end if
     allocate(psps%vlspl(psps%mqgrid_vl,2,dtset%ntypat),psps%qgrid_vl(psps%mqgrid_vl))
     if (.not.psps%vlspl_recipSpace) then
       allocate(psps%dvlspl(psps%mqgrid_vl,2,dtset%ntypat))
     end if
     mqgridvl_old=psps%mqgrid_vl
   end if
   if(ntypat_old/=dtset%ntypat.or. usepaw_old/=psps%usepaw)then
     if(idtset/=1)deallocate(psps%xccc1d)
     allocate(psps%xccc1d(n1xccc*(1-psps%usepaw),6,dtset%ntypat))
     usepaw_old=psps%usepaw
   end if
   if(ntypat_old/=dtset%ntypat)then
     if(idtset/=1)deallocate(psps%xcccrc,psps%ziontypat,psps%znucltypat)
     allocate(psps%xcccrc(dtset%ntypat))
     allocate(psps%znucltypat(dtset%ntypat))
     allocate(psps%ziontypat(dtset%ntypat))
     ntypat_old=dtset%ntypat
!    The correct dimension of pawrad/tab is ntypat.
!    In case of paw, no alchemical psp is allowed, so npsp=ntypat
!    However, in case of alchemical psps, pawrad/tab(ipsp) is invoked in
!    pspini. So, in order to avoid any problem, declare pawrad/tab
!    at paw_size=max(ntypat,npsp).
     paw_size=0;if (psps%usepaw==1) paw_size=max(dtset%ntypat,npsp)
   end if
   psps%ziontypat(:)=dtset%ziontypat(:)

!  ****************************************************************************
!  PAW allocations.

   call status(jdtset_status,filstat,iexit,level,'PAW allocs    ')

   if (paw_size/=paw_size_old) then
     if(idtset/=1) then
       call pawalloc(dtset,idtset,mpsang,psps%mqgrid_vl,npsp,2,paw_size,paw_size_old,&
&       pawang,pawrad,pawtab,pspheads)
       deallocate(pawrad,pawtab)
     end if
     allocate(pawrad(paw_size),pawtab(paw_size))
   end if
   call pawalloc(dtset,idtset,mpsang,psps%mqgrid_vl,npsp,1,paw_size,paw_size_old,&
&   pawang,pawrad,pawtab,pspheads)
   paw_size_old=paw_size

!  ****************************************************************************
!  WVL allocations.

   call status(jdtset_status,filstat,iexit,level,'WVL allocs    ')

!  Set up mpi informations from the dataset
   if (dtset%usewvl == 0) then
     nullify(mpi_enreg%nscatterarr, mpi_enreg%ngatherarr)
   else
!    WVL - data distribution
     allocate(mpi_enreg%nscatterarr(0:mpi_enreg%nproc - 1, 4))
     allocate(mpi_enreg%ngatherarr(0:mpi_enreg%nproc - 1, 2))
   end if

!  ****************************************************************************
!  ETSF I/O

   call status(jdtset_status,filstat,iexit,level,'ETSF IO       ')

!  Fill in abidims structure
   abidims%mband    = dtset%mband
   abidims%mproj    = 1                 ! FIXME
   abidims%mpsang   = mpsang
   abidims%mpw      = dtset%mpw
   abidims%natom    = dtset%natom
   abidims%natsph   = dtset%natsph
   abidims%nberry   = dtset%nberry
   abidims%nconeq   = dtset%nconeq
   abidims%nfft     = dtset%nfft
   abidims%nfreqsus = dtset%nfreqsus
   abidims%ngrid1   = dtset%ngfft(1)
   abidims%ngrid2   = dtset%ngfft(2)
   abidims%ngrid3   = dtset%ngfft(3)
   abidims%nkpt     = dtset%nkpt
   abidims%nkptgw   = dtset%nkptgw
   abidims%npsp     = npsp
   abidims%npspalch = dtset%npspalch
   abidims%nqptdm   = dtset%nqptdm
   abidims%nshiftk  = dtset%nshiftk
   abidims%nspden   = dtset%nspden
   abidims%nspinor  = dtset%nspinor
   abidims%nsppol   = dtset%nsppol
   abidims%nsym     = dtset%nsym
   abidims%ntypat   = dtset%ntypat
   abidims%ntypalch = dtset%ntypalch
   abidims%wfs_dim1 = -1
   abidims%wfs_dim2 = -1
   abidims%npw_tiny = -1


!  ****************************************************************************
!  At this stage, all the data needed for the treatment of one dataset
!  have been transferred from multi-dataset arrays.

   iexit=0

!  Smaller integer arrays :
   allocate(npwtot(dtset%nkpt))
   allocate(fcart(3,dtset%natom))
   allocate(fred(3,dtset%natom))

!  Initialize these to zero (needed when hasty exit)
   etotal=zero
   strten(:)=zero
   fcart(:,:)=zero ; fred(:,:)=zero

   if(dtset%optdriver/=0)then
     acell(:)=acell_img(:,1)
     occ(:)=occ_img(:,1)
     rprim(:,:)=rprim_img(:,:,1)
     vel(:,:)=vel_img(:,:,1)
     xred(:,:)=xred_img(:,:,1)
   end if

   call echo_xc_name (dtset%ixc)

#if defined HAVE_LIBXC
   if (dtset%ixc<0) then
     call libxc_functionals_init(dtset%ixc,dtset%nspden)
     if (libxc_functionals_isgga() .or. libxc_functionals_ismgga()) dtset%xclevel=2
   end if
#endif

   select case(dtset%optdriver)

     case(0)

       call status(jdtset_status,filstat,iexit,level,'call gstateimg')

       allocate(fcart_img(3,dtset%natom,nimage))
       allocate(fred_img(3,dtset%natom,nimage))
       allocate(etotal_img(nimage))
       allocate(strten_img(6,nimage))

       call gstateimg(acell_img,codvsn,cpui,dtfil,dtset,etotal_img,fcart_img,fred_img,&
&       iexit,mpi_enreg,npwtot,dtset%nspinor,&
&       occ_img,pawang,pawrad,pawtab,psps,rprim_img,strten_img,vel_img,xred_img)

       call status(jdtset_status,filstat,iexit,level,'after gstateimg')

     case(1)

       call status(jdtset_status,filstat,iexit,level,'call respfn   ')

       mpi_enreg%paral_level=2
       call respfn(codvsn,cpui,dtfil,dtset,etotal,iexit,mkmems,mpi_enreg,&
&       npwtot,dtset%nspinor,occ,pawang,pawrad,pawtab,psps,xred)

       call status(jdtset_status,filstat,iexit,level,'after respfn  ')

     case(2)

       call status(jdtset_status,filstat,iexit,level,'call suscep   ')

       mpi_enreg%paral_level=2
       call suscep(dtfil,dtset,iexit,&
&       dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,dtset%nfft,dtset%nkpt,&
&       dtset%nspden,dtset%nspinor,dtset%nsppol,dtset%nsym,occ,xred)

       call status(jdtset_status,filstat,iexit,level,'after suscep  ')

     case(RUNL_SCREENING)

       call status(jdtset_status,filstat,iexit,level,'call screening')

       call screening(acell,codvsn,dtfil,dtset,iexit,mpi_enreg,pawang,pawrad,pawtab,psps,rprim)

       call status(jdtset_status,filstat,iexit,level,'after screenin')

     case(RUNL_SIGMA)

       call status(jdtset_status,filstat,iexit,level,'call sigma    ')

       call sigma(acell,codvsn,dtfil,dtset,iexit,mpi_enreg,pawang,pawrad,pawtab,psps,rprim,xred,converged)

       call status(jdtset_status,filstat,iexit,level,'after sigma   ')

     case (RUNL_SCGW)

       call status(jdtset_status,filstat,iexit,level,'call gw_driver')

       call gw_driver(idtset,jdtset_,ndtset,acell,codvsn,filnam,dtfil,dtset,iexit,mpi_enreg,&
&       pawang,pawrad,pawtab,psps,rprim,xred)

       call status(jdtset_status,filstat,iexit,level,'after gw_driver')

     case(5)

       call status(jdtset_status,filstat,iexit,level,'call nonlinear   ')

       mpi_enreg%paral_level=2
       call nonlinear(codvsn,dtfil,dtset,etotal,iexit,&
&       dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,dtset%nfft,dtset%nkpt,npwtot,dtset%nspden,&
&       dtset%nspinor,dtset%nsppol,dtset%nsym,occ,pawrad,pawtab,psps,xred)

       call status(jdtset_status,filstat,iexit,level,'after nonlinear  ')

     case(6)

       write(message, '(a,i12,a,a,a,a)' )&
&       '  The optdriver value 6 has been disabled since ABINITv6.0.',ch10,&
&       '  Action : modify optdriver in the input file.'
       MSG_ERROR(message)

     case (RUNL_BSE) ! Bethe-Salpeter
       call status(jdtset_status,filstat,iexit,level,'call bethe_salpeter')
       
       call bethe_salpeter(acell,codvsn,dtfil,dtset,iexit,mpi_enreg,pawang,pawrad,pawtab,psps,rprim,xred)
       
       call status(jdtset_status,filstat,iexit,level,'after bethe_salpeter')

     case(RUNL_RDM)

       call rdm(acell,dtfil,dtset,pawtab,mpi_enreg,rprim)

       case default

!      Error the choice is either 0 -> gstate, 1 -> respfn, 2 -> suscep,
!      3 -> screening, 4 -> sigma,  5 -> nonlinear, 6 -> wannier
       write(message,'(a,i12,4a)')&
&       '  Unknown value for the variable optdriver: ',dtset%optdriver,ch10,&
&       '  This is not allowed. ',ch10,&
&       '  Action : modify optdriver in the input file.'
       MSG_ERROR(message)
   end select

!  ****************************************************************************

!  Transfer of multi dataset outputs from temporaries :
!  acell, xred, occ rprim, and vel might be modified from their
!  input values
!  etotal, fcart, fred, and strten have been computed
!  npwtot was already computed before, but is stored only now

   if(dtset%optdriver==0)then
     do iimage=1,nimage
       results_out(idtset)%acell(:,iimage)                =acell_img(:,iimage)
       results_out(idtset)%etotal(iimage)                 =etotal_img(iimage)
       results_out(idtset)%rprim(:,:,iimage)              =rprim_img(:,:,iimage)
       call mkrdim(acell_img(:,iimage),rprim_img(:,:,iimage),rprimd)
       results_out(idtset)%rprimd(:,:,iimage)             =rprimd(:,:)
       results_out(idtset)%strten(:,iimage)                =strten_img(:,iimage)
       results_out(idtset)%fcart(1:3,1:dtset%natom,iimage)=fcart_img(:,:,iimage)
       results_out(idtset)%fred(1:3,1:dtset%natom,iimage) =fred_img(:,:,iimage)
       results_out(idtset)%npwtot(1:dtset%nkpt,iimage)    =npwtot(1:dtset%nkpt)
       results_out(idtset)%occ(1:dtset%mband*dtset%nkpt*dtset%nsppol,iimage)=&
&       occ_img(1:dtset%mband*dtset%nkpt*dtset%nsppol,iimage)
       results_out(idtset)%vel(:,1:dtset%natom,iimage)    =vel_img(:,:,iimage)
       results_out(idtset)%xred(:,1:dtset%natom,iimage)   =xred_img(:,:,iimage)
     end do
     deallocate(etotal_img,fcart_img,fred_img,strten_img)
   else
     results_out(idtset)%acell(:,1)                =acell(:)
     results_out(idtset)%etotal(1)                 =etotal
     results_out(idtset)%rprim(:,:,1)              =rprim(:,:)
     call mkrdim(acell,rprim,rprimd)
     results_out(idtset)%rprimd(:,:,1)             =rprimd(:,:)
     results_out(idtset)%npwtot(1:dtset%nkpt,1)    =npwtot(1:dtset%nkpt)
     results_out(idtset)%occ(1:dtset%mband*dtset%nkpt*dtset%nsppol,1)=&
&     occ(1:dtset%mband*dtset%nkpt*dtset%nsppol)
     results_out(idtset)%xred(:,1:dtset%natom,1)   =xred(:,:)
   end if
   deallocate(fcart,fred)
   deallocate(acell_img,occ_img,rprim_img,vel_img,xred_img)

#if defined HAVE_LIBXC
   if (dtset%ixc<0) then
     call libxc_functionals_end()
   end if
#endif
   call dtsetFree(dtset)

   deallocate(psps%algalch,psps%mixalch)

   deallocate(occ)
   deallocate(vel,xred)
   deallocate(npwtot)

#ifdef HAVE_ABI_TIMER
   call print_timer()
   call destroy_timer()
#endif

   if(iexit/=0)exit

!  Check whether exiting was required by the user.
!  If found then beat a hasty exit from time steps
   openexit=1 ; if(dtset%chkexit==0) openexit=0

   call chkexi(zero,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)

   if (iexit/=0)exit

!  End do loop on idtset (allocate statements are present -
!  an exit statement is present)
 end do

!*********************************************************************

 deallocate(psps%qgrid_ff)
 deallocate(psps%qgrid_vl)
 deallocate(psps%xcccrc)
 deallocate(psps%xccc1d)
 deallocate(psps%vlspl)
 if (.not.psps%vlspl_recipSpace) then
   deallocate(psps%dvlspl)
 end if
 deallocate(psps%ekb)
 deallocate(psps%indlmn)
 deallocate(psps%filpsp)
 deallocate(psps%pspcod)
 deallocate(psps%pspdat)
 deallocate(psps%pspso)
 deallocate(psps%pspxc)
 deallocate(psps%title)
 deallocate(psps%znuclpsp)
 deallocate(psps%znucltypat)
 deallocate(psps%zionpsp)
 deallocate(psps%ziontypat)
 deallocate(psps%ffspl)
 call psp2params_free(psps%gth_params)

!PAW deallocation
 call pawalloc(dtset,idtset,mpsang,psps%mqgrid_vl,npsp,3,paw_size,paw_size_old,&
& pawang,pawrad,pawtab,pspheads)
 deallocate(pawrad,pawtab)

 deallocate(jdtset_)

 call status(0,filstat,iexit,level,'exit          ')
 call timab(100,2,tsec)

 DBG_EXIT("COLL")

end subroutine driver
!!***
