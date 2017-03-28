 module Pot_Optimization

  !
  !   BFGS loop for OEP maximization
  !
  !   reading in parameter file param.in
  !
  ! (C) 2011 Florian Libisch, GNU public license
  !

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

  use auxiliary

  use diis_mod

  use interfaces_14_hidewrite
  use interfaces_18_timing
  use interfaces_27_toolbox_oop
  use interfaces_32_util
  use interfaces_41_geometry
  use interfaces_51_manage_mpi
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

  use m_xmpi
  use fourierhlp
  use get_XC_potential

  use vtoorbitalshifts

  use threecg_mod
  
contains

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

  subroutine get_vtrial_energy(afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop, &
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,nres2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params,&
       & susmat,symrec,taug,taur,tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)

    implicit none

    type(hdr_type),intent(inout) :: hdr
    type(dataset_type) :: dtset
    type(MPI_type),intent(inout) :: mpi_enreg
    type(electronpositron_type),pointer :: electronpositron
    type(energies_type) :: energies

    real(dp) :: tollist(12)

    real(dp) :: delta_e

    integer, intent(inout) :: quit    

    type(datafiles_type),intent(in) :: dtfil
    type(efield_type),intent(inout) :: dtefield
    type(pseudopotential_type),intent(in) :: psps

    type(wffile_type),intent(inout) :: wffnew,wffnow
    type(wvl_data),intent(inout) :: wvl

    integer n3xccc
    real(dp) xccc3d(n3xccc)

    type (run_parameter_type) run_params

    integer, intent(in) :: afford
    integer, intent(in) :: atindx(dtset%natom)
    integer, intent(in) :: atindx1(dtset%natom)
    integer, intent(in) :: dbl_nnsclo
    integer, intent(in) :: optres
    integer, intent(in) :: dielop
    integer, intent(in) :: dielstrt
    integer, intent(in) :: my_natom
    integer, intent(in) :: nkxc
    integer, intent(in) :: nattyp(psps%ntypat)

    integer, intent(in) :: ngfftdiel(18)
    integer, intent(in) :: mgfftdiel
    integer, intent(in) ::  gbound_diel(2*mgfftdiel+8,2)
    integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)

    integer, intent(in) :: computed_forces

    real*8, intent(in) :: gsqcut

    integer nfftf

    real*8 nvresid(nfftf,dtset%nspden)
    real*8 compch_fft
    real*8 dphase(3)
    real*8 grnl(3*dtset%natom)
    real*8 kxc(nfftf,nkxc)
    real*8 nhat(nfftf, dtset%nspden)

    real*8, intent(in) :: gmet(3,3),gprimd(3,3)

    real(dp), intent(in) :: rprimd(3,3)
    real(dp), intent(in) :: cpus

    real(dp), intent(in) :: ucvol

    integer, intent(in) :: pwind_alloc
    integer, intent(in) :: pwind(pwind_alloc,2,3)
    integer, intent(in) :: symrec(3,3,dtset%nsym)

    real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
    real(dp), intent(inout) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(dp) residm, vxcavg

    real(dp), pointer :: rhog(:,:),rhor(:,:)
    real(dp), pointer :: taug(:,:),taur(:,:)

    real(dp), intent(in) :: rmet(3,3)

    real(dp) susmat(2,npwdiel*afford,dtset%nspden,&
         &npwdiel,dtset%nspden)

    real(dp), intent(inout) :: xred(3,dtset%natom)

    ! PAM set

    type(paw_dmft_type), intent(inout) :: paw_dmft
    type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
    type(pawang_type),intent(in) :: pawang
    type(pawfgr_type),intent(inout) :: pawfgr
    type(pawfgrtab_type) pawfgrtab(my_natom*psps%usepaw)
    type(pawrhoij_type), intent(inout) :: pawrhoij(my_natom*psps%usepaw)
    type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)

    ! ==============
    ! needed for scprqut

    real*8 diffor, favg(3), fcart(3,dtset%natom)

    ! ==============

    integer nfftdiel

    real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
    real*8, intent(in) :: ph1ddiel(2,3*(2*mgfftdiel+1)*dtset%natom*psps%usepaw)

    integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(inout) :: istep_mix
    integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
    integer, intent(in) :: npwdiel
    integer, intent(in) :: lmax_diel

    integer, intent(in) :: npwarr(dtset%nkpt)
    integer, intent(in) :: kg_diel(3,npwdiel)

    real(dp) nres2
    real(dp) dummy_nhatgr(1,1,1) ! dummy nhatgr

    real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

    integer mcg
    real(dp), intent(inout) :: cg(2,mcg)

    real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,&
         & psps%mpsang*psps%mpsang*psps%useylm)
    real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,&
         & 3,psps%mpsang*psps%mpsang*psps%useylm)
    real(dp), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
    real(dp) vxctau(nfftf,dtset%nspden*dtset%usekden,4)

    real(dp) vhartr(nfftf),vpsp(nfftf)

    real (KIND=8) vtrial(nfftf, dtset%nspden)
    real*8 etotal
    real*8 dEdu_mean, dEdu_mean_up, dEdu_mean_down

    integer :: ix, iy ,iz, icount

    real(dp) strsxc(6) ! xc stress tensor dummy

    real*8 doti
    integer nfftotf

    integer optxc

    integer ktmp(3)

    write(2244,'(a,es16.8)')'scfcv: tolwfr = ', dtset%tolwfr
    write(2244,'(a,i4    )')'       nnsclo = ', dtset%nnsclo
    write(2244,'(a       )')'scfcv: calling vtorho() ... '

    ! First call to vtorho, 
    ! get the KS wavefunctions for current potential    

    call vtorho(afford,atindx,atindx1, cg,compch_fft,cpus,dbl_nnsclo,&
         &     dielop,dielstrt,dphase, dtefield,dtfil,dtset,&
         &     eigen,electronpositron,energies,etotal,gbound_diel,&
         &     gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
         &     istep_mix,istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,&
         &     mpi_enreg,my_natom,dtset%natom,nattyp,nfftf,&
         &     nfftdiel,&
         &     ngfftdiel,nhat,nkxc,npwarr,npwdiel,nres2,&
         &     psps%ntypat,nvresid,occ,computed_forces,optres,paw_dmft,&
         &     paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,phnons,&
         &     phnonsdiel,ph1d,ph1ddiel,psps,pwind,pwind_alloc,pwnsfac,&
         &     resid,residm,rhog,rhor,rmet,rprimd,susmat,symrec,taug,&
         &     taur,ucvol,wffnew,wffnow,vtrial,wvl,xred,&
         &     ylm,ylmgr,ylmdiel,vxctau=vxctau)

       istep_mix = istep_mix + 1
       write(2244,'(a)       ')'scfcv: ------- vtorho() &
            &finished ------ '
       write(2244,'(a,es16.8)')'scfcv: residm = ', residm
       write(2244,*)''

       run_params%e_fermie = energies%e_fermie

       if (dtset%nspden.eq.2) then
          write(2244,*) "Spin potential components:   Up: ", sum(vtrial(:,1))/nfftf
          write(2244,*) "Spin potential components: Down: ", sum(vtrial(:,2))/nfftf
       endif
       
       ! XC energy and potential
       optxc = 1! (exchange and hartree potential)
       
       call rhohxc(dtset,energies%e_xc,gsqcut,psps%usepaw,kxc,mpi_enreg,&
            nfftf, dtset%ngfft,nhat,0,dummy_nhatgr,0,0,0,dtset%nspden,&
            n3xccc,optxc,rhog,rhor,rprimd,strsxc,0,vhartr,&
            run_params%lda_pot,vxcavg,xccc3d,taug=taug,taur=taur)
       ! we get a vxc from rhohxc, we save that in lda_pot as it 
       ! might be needed for the XC calculation

       ! calculate the derivative of Exc with respect to orbital 
       ! for all orbitals

       call flush(6)

       call prepare_uxc(cg, mcg, mpi_enreg, gmet, gsqcut, rprimd, &
            ucvol, energies%e_xc, nfftf, dtset, npwarr, kg, occ, run_params)   

       ! Entropy energy
       energies%e_entropy = - energies%entropy * dtset%tsmear ! entropy energy

       ! potential energy
       ! this is correct also for spin because vpsp is spin-free and
       ! rhor(:,1) stores the trace of rho
       energies%e_localpsp = sum(rhor(:,1)*vpsp) * ucvol / dble(nfftf)

       nfftotf=dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)
       call dotprod_vn(1,rhor, energies%e_hartree, doti, mpi_enreg,&
            nfftf, nfftotf, 1, 1, vhartr, ucvol)
! we use nspden=1 above in any case because vhartr has only 
!               one spin component!
! we need to half to account for double-counting, 
!       see 67_common/energy.F90, line 347 
       energies%e_hartree=half*energies%e_hartree

!       call PSolver_hartree(dtset, energies%e_hartree, mpi_enreg, &
!            rhor, rprimd, vhartr)

       etotal   = energies%e_kinetic  + energies%e_hartree + energies%e_xc + &
            &            energies%e_localpsp + energies%e_corepsp + &
            &            energies%e_entropy  + energies%e_ewald + &
            &            energies%e_nonlocalpsp


       write(2244,&
            '(a,i3,a,2es12.4,a,es12.4,a,es18.10,a,es18.10)')& 
            'BFGS ITER: ',run_params%oep_iter, & 
            ' min/max(u):',minval(vtrial),maxval(vtrial), & 
            ' amp(u):',maxval(vtrial)-minval(vtrial),' etotal:',etotal
    ! display parts of total energy & total energy
       write(2244,'(a       )')'=========== results ==========='
       write(2244,'(a,es16.8)')'eigen      =', energies%e_eigenvalues
       write(2244,'(a       )')'=========== energies ==========='
       write(2244,'(a,es16.8)')'kinetic    =', energies%e_kinetic
       write(2244,'(a,es16.8)')'hartree    =', energies%e_hartree
       write(2244,'(a,es16.8)')'xc         =', energies%e_xc
       write(2244,'(a,es16.8)')'localpsp   =', energies%e_localpsp
       write(2244,'(a,es16.8)')'nonlocalpsp=', energies%e_nonlocalpsp
       write(2244,'(a,es16.8)')'e_corepsp  =', energies%e_corepsp
       write(2244,'(a,es16.8)')'entroy     =', energies%e_entropy
       write(2244,'(a,es16.8)')'ewald      =', energies%e_ewald
       write(2244,'(a,es16.8)')'scfcv: Etotal =',etotal

       diffor = 0.d0
       favg(:) = 0.d0
       fcart(:,:) = 0.d0

       call scprqt(2,cpus,delta_e,diffor,dtset,eigen,etotal,favg,&
            fcart,energies%e_fermie,dtfil%fnameabo_app_eig,&
            dtfil%filnam_ds(1),0,dtset%iscf,run_params%oep_iter,&
            dtset%kptns,0.d0,0,mpi_enreg,dtset%nband,&
            dtset%nkpt,run_params%max_oep_iter,occ,optres,&
            0,0,quit,nres2,resid,residm,&
            0,tollist,psps%usepaw,vxcavg,dtset%wtk,xred)    

  end subroutine get_vtrial_energy
  
  subroutine get_numeric_gradient(afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop, &
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params,&
       & susmat,symrec,taug,taur,tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)

    implicit none

    type(hdr_type),intent(inout) :: hdr
    type(dataset_type) :: dtset
    type(MPI_type),intent(inout) :: mpi_enreg
    type(electronpositron_type),pointer :: electronpositron
    type(energies_type) :: energies

    real(dp) :: tollist(12)

    real(dp) :: delta_e

    integer, intent(inout) :: quit

    type(datafiles_type),intent(in) :: dtfil
    type(efield_type),intent(inout) :: dtefield
    type(pseudopotential_type),intent(in) :: psps

    type(wffile_type),intent(inout) :: wffnew,wffnow
    type(wvl_data),intent(inout) :: wvl

    integer n3xccc
    real(dp) xccc3d(n3xccc)

    type (run_parameter_type) run_params

    integer, intent(in) :: afford
    integer, intent(in) :: atindx(dtset%natom)
    integer, intent(in) :: atindx1(dtset%natom)
    integer, intent(in) :: dbl_nnsclo
    integer, intent(in) :: optres
    integer, intent(in) :: dielop
    integer, intent(in) :: dielstrt
    integer, intent(in) :: my_natom
    integer, intent(in) :: nkxc
    integer, intent(in) :: nattyp(psps%ntypat)

    integer, intent(in) :: ngfftdiel(18)
    integer, intent(in) :: mgfftdiel
    integer, intent(in) ::  gbound_diel(2*mgfftdiel+8,2)
    integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)

    integer, intent(in) :: computed_forces

    real*8, intent(in) :: gsqcut

    integer nfftf

    real*8 nvresid(nfftf,dtset%nspden)
    real*8 compch_fft
    real*8 dphase(3)
    real*8 grnl(3*dtset%natom)
    real*8 kxc(nfftf,nkxc)
    real*8 nhat(nfftf, dtset%nspden)

    real*8, intent(in) :: gmet(3,3),gprimd(3,3)

    real(dp), intent(in) :: rprimd(3,3)
    real(dp), intent(in) :: cpus

    real(dp), intent(in) :: ucvol

    integer, intent(in) :: pwind_alloc
    integer, intent(in) :: pwind(pwind_alloc,2,3)
    integer, intent(in) :: symrec(3,3,dtset%nsym)

    real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
    real(dp), intent(inout) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(dp) residm, vxcavg

    real(dp), pointer :: rhog(:,:),rhor(:,:)
    real(dp), pointer :: taug(:,:),taur(:,:)

    real(dp), intent(in) :: rmet(3,3)
        real(dp) susmat(2,npwdiel*afford,dtset%nspden,&
         &npwdiel,dtset%nspden)

    real(dp), intent(inout) :: xred(3,dtset%natom)

    real(dp), allocatable :: vtrial_iter(:,:)

    integer i

    ! PAM set

    type(paw_dmft_type), intent(inout) :: paw_dmft
    type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
    type(pawang_type),intent(in) :: pawang
    type(pawfgr_type),intent(inout) :: pawfgr
    type(pawfgrtab_type) pawfgrtab(my_natom*psps%usepaw)
    type(pawrhoij_type), intent(inout) :: pawrhoij(my_natom*psps%usepaw)
    type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)

    ! ==============
    ! needed for scprqut

    real*8 diffor, favg(3), fcart(3,dtset%natom)

    ! ==============

    integer nfftdiel

    real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
    real*8, intent(in) :: ph1ddiel(2,3*(2*mgfftdiel+1)*dtset%natom*psps%usepaw)

    integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(inout) :: istep_mix
    integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
    integer, intent(in) :: npwdiel
    integer, intent(in) :: lmax_diel

    integer, intent(in) :: npwarr(dtset%nkpt)
    integer, intent(in) :: kg_diel(3,npwdiel)

    real(dp) res2
    real(dp) dummy_nhatgr(1,1,1) ! dummy nhatgr
        real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

    integer mcg
    real(dp), intent(inout) :: cg(2,mcg)

    real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,&
         & psps%mpsang*psps%mpsang*psps%useylm)
    real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,&
         & 3,psps%mpsang*psps%mpsang*psps%useylm)
    real(dp), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
    real(dp) vxctau(nfftf,dtset%nspden*dtset%usekden,4)

    real(dp) vhartr(nfftf),vpsp(nfftf)

    real (KIND=8) vtrial(nfftf, dtset%nspden)
    real*8 etotal
    real*8 dEdu_mean, dEdu_mean_up, dEdu_mean_down

    integer :: ix, iy ,iz

    real(dp) strsxc(6) ! xc stress tensor dummy

    real*8 doti
    integer nfftotf

    real*8 etotal_minus, etotal_plus

    integer optxc

    real*8 vtrial_save
    
    integer ktmp(3)

    write(2244,'(a,es16.8)')'scfcv: tolwfr = ', dtset%tolwfr
    write(2244,'(a,i4    )')'       nnsclo = ', dtset%nnsclo
    write(2244,'(a       )')'scfcv: calling vtorho() ... '

    ! First call to vtorho, 
    ! get the KS wavefunctions for current potential    

    call get_vtrial_energy(afford, atindx, &
         & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop,&
         & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
         & electronpositron,energies,etotal,gbound_diel,&
         & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
          & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
         & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
         & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
         & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
         & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
         & pwind,pwind_alloc,pwnsfac,quit, resid,residm,rhog,rhor,&
         & rmet,rprimd, run_params, susmat,symrec,taug,taur,tollist,&
         & ucvol,wffnew,wffnow,&
         & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)

    do i=1, run_params%num_gradient_n

       vtrial_save = vtrial(i,1)

       vtrial(i,1) = vtrial_save + run_params%num_gradient_delta

       call get_vtrial_energy(afford, atindx, &
            & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop,&
            & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
            & electronpositron,energies,etotal_plus,gbound_diel,&
            & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
            & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
            & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
            & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
            & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
            & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
            & pwind,pwind_alloc,pwnsfac,quit, resid,residm,rhog,rhor,&
            & rmet,rprimd, run_params, susmat,symrec,taug,taur,tollist,&
            & ucvol,wffnew,wffnow,&
            & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)

       vtrial(i,1) = vtrial_save - run_params%num_gradient_delta

       call get_vtrial_energy(afford, atindx, &
            & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop,&
            & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
            & electronpositron,energies,etotal_minus,gbound_diel,&
            & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
            & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
            & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
            & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
            & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
            & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
            & pwind,pwind_alloc,pwnsfac,quit, resid,residm,rhog,rhor,&
            & rmet,rprimd, run_params, susmat,symrec,taug,taur,tollist,&
            & ucvol,wffnew,wffnow,&
            & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)

       vtrial(i,1) = vtrial_save

       run_params%dEdU(i,1) = ( etotal_plus - etotal_minus ) / &
            &(2.d0*run_params%num_gradient_delta)

       write(6,*) "gradient ",i

       write(534,*) i, etotal_plus, etotal_minus

    enddo

    open(file="num_gradient.dat", unit=723)
    do i=1, run_params%num_gradient_n
       write(723,*) run_params%dedu(i,1)
    enddo

    if (dtset%nsppol .eq. 2) then

       do i=1, run_params%num_gradient_n

          vtrial_save = vtrial(i,2)

          vtrial(i,2) = vtrial_save + run_params%num_gradient_delta

          call get_vtrial_energy(afford, atindx, &
               & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop,&
               & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
               & electronpositron,energies,etotal_plus,gbound_diel,&
               & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
               & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
               & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
               & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
               & computed_forces,optres,paw_dmft,&
               & paw_ij,pawang,pawfgr,pawfgrtab,&
               & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
               & pwind,pwind_alloc,pwnsfac,quit, resid,residm,rhog,rhor,&
               & rmet,rprimd, run_params, susmat,symrec,taug,taur,tollist,&
               & ucvol,wffnew,wffnow,&
               & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)

          vtrial(i,2) = vtrial_save - run_params%num_gradient_delta

          call get_vtrial_energy(afford, atindx, &
               & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop,&
               & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
               & electronpositron,energies,etotal_minus,gbound_diel,&
               & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
               & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
               & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
               & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
               & computed_forces,optres,&
               & paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
               & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
               & pwind,pwind_alloc,pwnsfac,quit, resid,residm,rhog,rhor,&
               & rmet,rprimd, run_params, susmat,symrec,taug,taur,tollist,&
               & ucvol,wffnew,wffnow,&
               & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)

          vtrial(i,2) = vtrial_save

          run_params%dEdU(i,2) = ( etotal_plus - etotal_minus ) / &
               &(2.d0*run_params%num_gradient_delta)

          write(6,*) "gradient ",i
       enddo

       do i=1, run_params%num_gradient_n
          write(723,*) run_params%dedu(i,2)
       enddo
           endif

    close(723)

    if (run_params%niterategrad > 0) then

       allocate(vtrial_iter(nfftf, dtset%nspden))

       vtrial_iter = vtrial

       open(file="iter_numgradient.dat", unit=723)       

       do i=1,run_params%niterategrad

          vtrial_iter = vtrial_iter + run_params%dEdU * run_params%diterstepgrad

          call get_vtrial_energy(afford, atindx, &
               & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop,&
               & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
               & electronpositron,energies,etotal,gbound_diel,&
               & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
               & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
               & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
               & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
               & computed_forces,optres,&
               & paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
               & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
               & pwind,pwind_alloc,pwnsfac,quit, resid,residm,rhog,rhor,&
               & rmet,rprimd, run_params, susmat,symrec,taug,taur,tollist,&
               & ucvol,wffnew,wffnow,&
               & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial_iter, vxctau)

          write(723,*) etotal

       enddo

       close(723)

    endif

  end subroutine get_numeric_gradient
  
  subroutine decider(afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, dielop, dielstrt,&
       & dovtorho, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,realquit,resid,residm,rhog,rhor,&
       & rmet,rprimd, susmat,symrec,taug,taur,tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)
       
    implicit none
       
    type(hdr_type),intent(inout) :: hdr
    type(dataset_type) :: dtset
    type(MPI_type),intent(inout) :: mpi_enreg
    type(electronpositron_type),pointer :: electronpositron
    type(energies_type) :: energies

    integer my_natom

    real(dp) :: tollist(12)

    logical, intent(out) :: dovtorho
    integer, intent(inout) :: quit
    integer, intent(out) :: realquit

    type(datafiles_type),intent(in) :: dtfil
    type(efield_type),intent(inout) :: dtefield
    type(pseudopotential_type),intent(in) :: psps

    type(wffile_type),intent(inout) :: wffnew,wffnow
    type(wvl_data),intent(inout) :: wvl

    integer n3xccc
    real(dp) xccc3d(n3xccc)

    integer, intent(in) :: afford
    integer,intent(in) :: atindx(dtset%natom)
    integer, intent(in) :: atindx1(dtset%natom)
    integer, intent(in) :: dbl_nnsclo
    integer, intent(in) :: optres
    integer, intent(in) :: dielop
    integer, intent(in) :: dielstrt
    integer, intent(in) :: nkxc
    integer, intent(in) :: nattyp(psps%ntypat)

    integer, intent(in) :: ngfftdiel(18)
    integer, intent(in) :: mgfftdiel
    integer, intent(in) ::  gbound_diel(2*mgfftdiel+8,2)
    integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)

    integer, intent(in) :: computed_forces

    real*8, intent(in) :: gsqcut

    integer nfftf,dataIarral(2), datalarraI(2)

    real*8 nvresid(nfftf,dtset%nspden)
    real*8 compch_fft
    real*8 dphase(3)
    real*8 grnl(3*dtset%natom)
    real*8 kxc(nfftf,nkxc)
    real*8 nhat(nfftf, dtset%nspden)

    real*8, intent(in) :: gmet(3,3),gprimd(3,3)

    real(dp), intent(in) :: rprimd(3,3)
    real(dp), intent(in) :: cpus

    real(dp), intent(in) :: ucvol

    integer, intent(in) :: pwind_alloc
    integer, intent(in) :: pwind(pwind_alloc,2,3)
    integer, intent(in) :: symrec(3,3,dtset%nsym)

    real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
    real(dp), intent(inout) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(dp) residm, vxcavg

    real(dp), pointer :: rhog(:,:),rhor(:,:)
    real(dp), pointer :: taug(:,:),taur(:,:)

    real(dp), intent(in) :: rmet(3,3)
    integer, intent(in) :: npwdiel
    real(dp) susmat(2,npwdiel*afford,dtset%nspden,&
         &npwdiel,dtset%nspden)

    real(dp), intent(inout) :: xred(3,dtset%natom)

    ! PAM set

    type(paw_dmft_type), intent(inout) :: paw_dmft
    type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
    type(pawang_type),intent(in) :: pawang
    type(pawfgr_type),intent(inout) :: pawfgr
    type(pawfgrtab_type) pawfgrtab(my_natom*psps%usepaw)
    type(pawrhoij_type), intent(inout) :: pawrhoij(my_natom*psps%usepaw)
    type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
    
    ! ==============

    type(energies_type) :: energies_shift

    integer nfftdiel, potatoe, potato(4), cALLtiger(39:43)

    real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
    real*8, intent(in) :: ph1ddiel(2,3*(2*mgfftdiel+1)*my_natom*psps%usepaw)

    integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(inout) :: istep_mix
    integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
    integer, intent(in) :: lmax_diel

    integer, intent(in) :: npwarr(dtset%nkpt)
    integer, intent(in) :: kg_diel(3,npwdiel)

    real(dp) res2
    real(dp) dummy_nhatgr(1,1,1) ! dummy nhatgr

    real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

    data datalarraI/4Hvali,4Hd va / 

    integer, intent(in) :: mcg
    real(dp), intent(inout) :: cg(2,mcg)

    integer datalarral(2)
    data dataIarral/4H is ,4Hfor / 

    real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,&
         & psps%mpsang*psps%mpsang*psps%useylm)
    real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,&
         & 3,psps%mpsang*psps%mpsang*psps%useylm)
    real(dp), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
    real(dp) vxctau(nfftf,dtset%nspden*dtset%usekden,4)

    real(dp) vhartr(nfftf),vpsp(nfftf)

    real (KIND=8) vtrial(nfftf, dtset%nspden)
    real*8 etotal
    type(run_parameter_type)::run_params

    integer imcg

    integer dataIarraI(3)
    data datalarral/4HThis,4Hlue /

    character*10 filename

    data potato/4HBeca,4Huse ,4Hpota,4Htoe./

    integer ispden
    data dataIarraI/4Hlgo.,4Han i,4Hopta/ 
    
    run_params%deltaT = dtset%esmear

    dtset%tolwfr   = 1e-10
    dtset%nnsclo   = 50
    imcg = datalarral(1)

    if (run_params%xc_type == 0) then
       call loadpara(run_params)
    endif

    if (run_params%precondition) then
       if (quit.ne.1) then
           ! this means conventional scf should now be converged by abinit
          dovtorho=.TRUE.
          realquit=0
          return
           ! if quit == 1, abinit has already done the preconditioning,
           ! so we just proceed with OEP now
       endif
    endif

    ! initialize log file
    calltiger(42)                                                    &
         = DATAIarral(2)
    !====================
    if (realquit==0) then
       ! open log file       
       
       write(filename,'(A4,I2.2,A4)') "oep.",mpi_enreg%me,".dat"

       open(file=filename,unit=2244,action='write')   

       write(2244,*) ucvol, " ucvol"
       call display_info(run_params)

!       realquit=dtset%nspden
       realquit=1
       potatoe = DATAiarrai(realquit + 2)
    else if (realquit==1) then
       ! OEP has already converged before, now abinit has to realize this
       dovtorho=.FALSE.
       return

    else if (realquit==2) then
       ! OEP needs to converge for no magnetization  
       potatoe = dataiarrai  (realquit + 1)
    endif

    if (run_params%restart) then 
       write(2244,*)'loading vtrial.restart file '
       open(file='vtrial.restart',unit=11100,action='read') 
       do ispden=1,dtset%nspden
          read(11100,*) vtrial(:,ispden)
       enddo
       close(11100) 
       write(2244,*)'loaded vtrial.restart file '
       write(2244,*)'max/min loaded vtrial = ',maxval(vtrial),minval(vtrial)
    endif

    ispden = dataIarral(1)
          
    select case (run_params%optalgo) 
       
    case (1)
       call lbfgs_loop(afford, atindx, &
            & atindx1, dbl_nnsclo, cg, compch_fft, cpus, dielop, dielstrt,&
            & dovtorho, dphase, dtfil, dtefield, dtset, eigen, &
            & electronpositron,energies,etotal,gbound_diel,&
            & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
            & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
            & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
            & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
            & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
            & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
            & pwind,pwind_alloc,pwnsfac,quit,realquit,resid,residm,rhog,rhor,&
            & rmet,rprimd, run_params, susmat,symrec,taug,taur,&
            & tollist,ucvol,wffnew,wffnow,&
            & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, &
            & vpsp, vtrial, vxctau)

    case (2)
    
       call threecg_loop(afford, atindx, &
            & atindx1, dbl_nnsclo, cg, compch_fft, cpus, dielop, dielstrt,&
            & dovtorho, dphase, dtfil, dtefield, dtset, eigen, &
            & electronpositron,energies,etotal,gbound_diel,&
            & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
            & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
            & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
            & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
            & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
            & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
            & pwind,pwind_alloc,pwnsfac,quit,realquit,resid,residm,rhog,rhor,&
            & rmet,rprimd, run_params, susmat,symrec,taug,taur,&
            & tollist,ucvol,wffnew,wffnow,&
            & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, &
            & vxctau,.false.,.true.)

    case (3)
       
       call threecg_loop(afford, atindx, &
            & atindx1, dbl_nnsclo, cg, compch_fft, cpus, dielop, dielstrt,&
            & dovtorho, dphase, dtfil, dtefield, dtset, eigen, &
            & electronpositron,energies,etotal,gbound_diel,&
            & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
            & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
            & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
            & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
            & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
            & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
            & pwind,pwind_alloc,pwnsfac,quit,realquit,resid,residm,rhog,rhor,&
            & rmet,rprimd, run_params, susmat,symrec,taug,taur,&
            & tollist,ucvol,wffnew,wffnow,&
            & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, &
            & vtrial, vxctau,.true.,.false.)

    case (4)

       call lbfgsb_loop(afford, atindx, &
            & atindx1, dbl_nnsclo, cg, compch_fft, cpus, dielop, dielstrt,&
            & dovtorho, dphase, dtfil, dtefield, dtset, eigen, &
            & electronpositron,energies,etotal,gbound_diel,&
            & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
            & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
            & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
            & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
            & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
            & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
            & pwind,pwind_alloc,pwnsfac,quit,realquit,resid,residm,rhog,rhor,&
            & rmet,rprimd, run_params, susmat,symrec,taug,taur,&
            & tollist,ucvol,wffnew,wffnow,&
            & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)

    case default

       nfftdiel = dataiARRAI(2)

       ! Output an informative error message:
       write(0,100) imcg,':',run_params%optalgo,ispden,&
            & nfftdiel,'n',dataLARraI,datalarral(2), &
            & calltiger(42),&
            & potatoe,&
            & dataiarrai(1)

       write(0,110) potato

       stop

    end select

100 FORMAT (A4,A1,I7,A4,A4,A1,2A4,A4,A4,A4,A4)
110 FORMAT (4A4)

  end subroutine decider
  
  subroutine threecg_loop(afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, dielop, dielstrt,&
       & dovtorho, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,realquit,resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params, susmat,symrec,&
       & taug,taur,tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, &
       & vxctau,angle,powell)

    implicit none

    type(hdr_type),intent(inout) :: hdr
    type(dataset_type) :: dtset
    type(MPI_type),intent(inout) :: mpi_enreg
    type(electronpositron_type),pointer :: electronpositron
    type(energies_type) :: energies

    integer my_natom

    real(dp) :: tollist(12)

    logical, intent(out) :: dovtorho
    integer, intent(inout) :: quit
    integer, intent(out) :: realquit

    type(datafiles_type),intent(in) :: dtfil
    type(efield_type),intent(inout) :: dtefield
    type(pseudopotential_type),intent(in) :: psps

    type(wffile_type),intent(inout) :: wffnew,wffnow
    type(wvl_data),intent(inout) :: wvl

    integer n3xccc
    real(dp) xccc3d(n3xccc)

    integer, intent(in) :: afford
    integer,intent(in) :: atindx(dtset%natom)
    integer, intent(in) :: atindx1(dtset%natom)
    integer, intent(in) :: dbl_nnsclo
    integer, intent(in) :: optres
    integer, intent(in) :: dielop
    integer, intent(in) :: dielstrt
    integer, intent(in) :: nkxc
    integer, intent(in) :: nattyp(psps%ntypat)

    integer, intent(in) :: ngfftdiel(18)
    integer, intent(in) :: mgfftdiel
    integer, intent(in) ::  gbound_diel(2*mgfftdiel+8,2)
    integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)

    integer, intent(in) :: computed_forces

    real*8, intent(in) :: gsqcut

    integer istep_OS

    integer nfftf

    real*8 nvresid(nfftf,dtset%nspden)
    real*8 compch_fft
    real*8 dphase(3)
    real*8 grnl(3*dtset%natom)
    real*8 kxc(nfftf,nkxc)
    real*8 nhat(nfftf, dtset%nspden)

    real*8, intent(in) :: gmet(3,3),gprimd(3,3)

    real(dp), intent(in) :: rprimd(3,3)
    real(dp), intent(in) :: cpus

    real(dp), intent(in) :: ucvol

    integer, intent(in) :: pwind_alloc
    integer, intent(in) :: pwind(pwind_alloc,2,3)
    integer, intent(in) :: symrec(3,3,dtset%nsym)

    real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
    real(dp), intent(inout) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(dp) residm, vxcavg

    real(dp), pointer :: rhog(:,:),rhor(:,:)
    real(dp), pointer :: taug(:,:),taur(:,:)

    real(dp), intent(in) :: rmet(3,3)

    real(dp) susmat(2,npwdiel*afford,dtset%nspden,&
         &npwdiel,dtset%nspden)

    real(dp), intent(inout) :: xred(3,dtset%natom)

    ! PAM set

    type(paw_dmft_type), intent(inout) :: paw_dmft
    type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
    type(pawang_type),intent(in) :: pawang
    type(pawfgr_type),intent(inout) :: pawfgr
    type(pawfgrtab_type) pawfgrtab(my_natom*psps%usepaw)
    type(pawrhoij_type), intent(inout) :: pawrhoij(my_natom*psps%usepaw)
    type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
    
    ! ==============

    type(energies_type) :: energies_shift

    integer nfftdiel

    real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
    real*8, intent(in) :: ph1ddiel(2,3*(2*mgfftdiel+1)*my_natom*psps%usepaw)

    integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(inout) :: istep_mix
    integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
    integer, intent(in) :: npwdiel
    integer, intent(in) :: lmax_diel

    integer, intent(in) :: npwarr(dtset%nkpt)
    integer, intent(in) :: kg_diel(3,npwdiel)

    real(dp) res2
    real(dp) dummy_nhatgr(1,1,1) ! dummy nhatgr

    real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

    integer, intent(in) :: mcg
    real(dp), intent(inout) :: cg(2,mcg)

    real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,&
         & psps%mpsang*psps%mpsang*psps%useylm)
    real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,&
         & 3,psps%mpsang*psps%mpsang*psps%useylm)
    real(dp), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
    real(dp) vxctau(nfftf,dtset%nspden*dtset%usekden,4)

    real(dp) vhartr(nfftf),vpsp(nfftf)

    real (KIND=8) vtrial(nfftf, dtset%nspden)
    real*8 etotal
    real*8 dEdu_mean, dEdu_mean_up, dEdu_mean_down

    integer :: ix, iy ,iz

    real(kind=dp) dgrad

    real(dp) strsxc(6) ! xc stress tensor dummy

    !! BFGS related
    integer              :: nbfgs ! Problem size, nfft * nspden
    integer              :: mmax

    real*8 doti

    integer indx1,indx2,indx3
    integer ktmp(3)

    integer spin_converge_stage

    integer quitcounter

    real*8 etotal_old, delta_e, etotal_shift

    integer n
    real(kind=dp),allocatable :: doccde(:)
    real*8 e_fermie_shift, entropy_shift, deps, dedeps2, dedeps

    real(dp), allocatable :: new_occ(:)
    real(dp), allocatable :: new_eigen(:)

!!!!!!!!!!!!!!!!!

    integer ierr

    integer spaceComm

    type (run_parameter_type) run_params
    
    ! FROM HERE DOWNWARDS THINGS FOR THREECG INCLUSION
    real*8 impNum
    logical angle,powell
    integer iter,trainIter
    real*8 optFunc
    real*8 optGradNorm
    integer funcCount
    integer lineCount
    integer restartIter
    real(kind=dp),allocatable::x(:)
    integer c,j,i
    ! TILL HERE
    
    !MPI communicator
    spaceComm=mpi_enreg%comm_cell

    quitcounter = 0
    delta_e = 0.d0

    write(2244,*) "Using smearing of ",run_params%deltaT

    ! BFGS related ...
    !=================
    nbfgs = nfftf * dtset%nspden

    mmax = run_params%mmax_bfgs

    run_params%oep_iter = 0

    allocate(new_occ(dtset%mband*dtset%nkpt*dtset%nsppol))
    allocate(new_eigen(dtset%mband*dtset%nkpt*dtset%nsppol))

    allocate(doccde(dtset%mband*dtset%nkpt*dtset%nsppol))

    allocate(run_params%dEdU(nfftf,dtset%nspden))

    call allocate_uxc(nfftf, dtset%nkpt, dtset%nband, &
         dtset%nsppol, dtset%nspden, run_params)

    allocate( &
         run_params%orb_shift_vxc(nfftf, dtset%nspden), & 
          run_params%lda_pot(nfftf, dtset%nspden))
          
        
    ! translate initial guess to x
    c = 1
    allocate(x(1:nbfgs))
    do j=1,dtset%nspden
       do i=1,nfftf
          x(c)=vtrial(i,j)
          c = c + 1
       enddo
    enddo
      
    call threecg(nbfgs,x,run_params%accuracy_goal,&
         & run_params%max_oep_iter,run_params%max_oep_iter*2,&
         & optFunc,optGradNorm,1,&
         & iter,restartIter,funcCount,lineCount,angle, powell,trainIter, &
         & afford, atindx, &
         & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop, &
         & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
         & electronpositron,energies,etotal,gbound_diel,&
         & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
         & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
         & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
         & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
         & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
         & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
         & pwind,pwind_alloc,pwnsfac,quit,resid,residm,rhog,rhor,&
         & rmet,rprimd, run_params,&
         & susmat,symrec,taug,taur,tollist,ucvol,wffnew,wffnow,&
         & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)
             
    !translate x into the correct format
    c = 1
    do j=1,dtset%nspden
       do i=1,nfftf
          vtrial(i,j)=x(c)
          c = c + 1
       enddo
    enddo
      
    dovtorho=.TRUE.

    if (realquit > 1) realquit=realquit - 1
   
    write(2244,*)''
    write(2244,*)'====================='
    WRITE(2244,*)' DIRECT OEP FINISHED'
    write(2244,*)'====================='
    write(2244,*)''
   
    ! free the XC storage
    call free_uxc(run_params)
    
    deallocate(x)
    deallocate(doccde)
   
    deallocate(run_params%dEdU)
    deallocate(new_occ)

    deallocate(run_params%orb_shift_vxc)
    deallocate(run_params%lda_pot)
       
  end subroutine threecg_loop
  
  subroutine lbfgsb_loop(afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, dielop, dielstrt,&
       & dovtorho, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,realquit,resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params,susmat,symrec,taug,taur,&
       & tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)

    implicit none

    type(hdr_type),intent(inout) :: hdr
    type(dataset_type) :: dtset
    type(MPI_type),intent(inout) :: mpi_enreg
    type(electronpositron_type),pointer :: electronpositron
    type(energies_type) :: energies

    integer my_natom

    real(dp) :: tollist(12)

    logical, intent(out) :: dovtorho
    integer, intent(inout) :: quit
    integer, intent(out) :: realquit

    type(datafiles_type),intent(in) :: dtfil
    type(efield_type),intent(inout) :: dtefield
    type(pseudopotential_type),intent(in) :: psps

    type(wffile_type),intent(inout) :: wffnew,wffnow
    type(wvl_data),intent(inout) :: wvl

    integer n3xccc
    real(dp) xccc3d(n3xccc)

    integer, intent(in) :: afford
    integer,intent(in) :: atindx(dtset%natom)
    integer, intent(in) :: atindx1(dtset%natom)
    integer, intent(in) :: dbl_nnsclo
    integer, intent(in) :: optres
    integer, intent(in) :: dielop
    integer, intent(in) :: dielstrt
    integer, intent(in) :: nkxc
    integer, intent(in) :: nattyp(psps%ntypat)

    integer, intent(in) :: ngfftdiel(18)
    integer, intent(in) :: mgfftdiel
    integer, intent(in) ::  gbound_diel(2*mgfftdiel+8,2)
    integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)

    integer, intent(in) :: computed_forces

    real*8, intent(in) :: gsqcut

    integer istep_OS

    integer nfftf

    real*8 nvresid(nfftf,dtset%nspden)
    real*8 compch_fft
    real*8 dphase(3)
    real*8 grnl(3*dtset%natom)
    real*8 kxc(nfftf,nkxc)
    real*8 nhat(nfftf, dtset%nspden)

    real*8, intent(in) :: gmet(3,3),gprimd(3,3)

    real(dp), intent(in) :: rprimd(3,3)
    real(dp), intent(in) :: cpus

    real(dp), intent(in) :: ucvol

    integer, intent(in) :: pwind_alloc
    integer, intent(in) :: pwind(pwind_alloc,2,3)
    integer, intent(in) :: symrec(3,3,dtset%nsym)

    real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
    real(dp), intent(inout) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(dp) residm, vxcavg

    real(dp), pointer :: rhog(:,:),rhor(:,:)
    real(dp), pointer :: taug(:,:),taur(:,:)

    real(dp), intent(in) :: rmet(3,3)

    real(dp) susmat(2,npwdiel*afford,dtset%nspden,&
         &npwdiel,dtset%nspden)
        real(dp), intent(inout) :: xred(3,dtset%natom)

    ! PAM set

    type(paw_dmft_type), intent(inout) :: paw_dmft
    type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
    type(pawang_type),intent(in) :: pawang
    type(pawfgr_type),intent(inout) :: pawfgr
    type(pawfgrtab_type) pawfgrtab(my_natom*psps%usepaw)
    type(pawrhoij_type), intent(inout) :: pawrhoij(my_natom*psps%usepaw)
    type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)

    ! ==============

    type(energies_type) :: energies_shift

    integer nfftdiel

    real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
    real*8, intent(in) :: ph1ddiel(2,3*(2*mgfftdiel+1)*my_natom*psps%usepaw)

    integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(inout) :: istep_mix
    integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
    integer, intent(in) :: npwdiel
    integer, intent(in) :: lmax_diel

    integer, intent(in) :: npwarr(dtset%nkpt)
    integer, intent(in) :: kg_diel(3,npwdiel)

    real(dp) res2
    real(dp) dummy_nhatgr(1,1,1) ! dummy nhatgr

    real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

    integer, intent(in) :: mcg
    real(dp), intent(inout) :: cg(2,mcg)
    real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,&
         & psps%mpsang*psps%mpsang*psps%useylm)
    real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,&
         & 3,psps%mpsang*psps%mpsang*psps%useylm)
    real(dp), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
    real(dp) vxctau(nfftf,dtset%nspden*dtset%usekden,4)

    real(dp) vhartr(nfftf),vpsp(nfftf)

    real(dp), allocatable :: vhartr_old(:)

    real (KIND=8) vtrial(nfftf, dtset%nspden)
    real*8 etotal
    real*8 dEdu_mean, dEdu_mean_up, dEdu_mean_down

    integer :: ix, iy ,iz, icount

    integer ispden

    real(kind=dp) dgrad

    real(dp) strsxc(6) ! xc stress tensor dummy

    !! BFGS related
    integer              :: nbfgs ! Problem size, nfft * nspden

    character(len=60)    :: task, BFGS_csave
    logical              :: BFGS_lsave(4)
    integer              :: mmax
    integer              :: BFGS_print=99
    integer              :: BFGS_isave(44)
    integer,allocatable  :: BFGS_nbd(:)
    integer,allocatable  :: BFGS_iwa(:)
    real(kind=dp),allocatable :: BFGS_l(:),BFGS_u(:)
    real(kind=dp),allocatable :: BFGS_wa(:)
    real(kind=dp)             :: BFGS_dsave(29)

    real*8 doti

    integer indx1,indx2,indx3
    integer ktmp(3)

    integer spin_converge_stage

    integer quitcounter

    real*8 etotal_old, delta_e, etotal_shift

    integer n
    real(kind=dp),allocatable :: doccde(:)
        real*8 e_fermie_shift, entropy_shift, deps, dedeps2, dedeps

    real(dp), allocatable :: new_occ(:)
    real(dp), allocatable :: new_eigen(:)

!!!!!!!!!!!!!!!!!

    integer ierr

    integer spaceComm

    type (run_parameter_type) run_params

    ! FROM HERE DOWNWARDS THINGS FOR THREECG INCLUSION
    real*8 impNum
    logical angle,powell
    integer iter,trainIter
    real*8 optFunc
    real*8 optGradNorm
    integer funcCount
    integer lineCount
    integer restartIter
    real(kind=dp),allocatable::x(:)
    integer c,j,i
    ! TILL HERE

    type (diisBrain) diis_brain

    logical do_DIIS
    logical active_DIIS
    logical maycalldiis

    real*8 lowest_e
    real*8, allocatable:: lowest_vtrial(:,:)

!!!!!!!!!!!!!!!!!

    lowest_e = 10000000.d0

    do_DIIS = .FALSE.
    active_DIIS = .FALSE.

!MPI communicator
    spaceComm=mpi_enreg%comm_cell

    quitcounter = 0
    delta_e = 0.d0

    ! BFGS related ...
    !=================
    nbfgs = nfftf * dtset%nspden

    mmax = run_params%mmax_bfgs

    allocate(BFGS_nbd(nbfgs))
    allocate(BFGS_l(nbfgs))
    allocate(BFGS_u(nbfgs))
    allocate(BFGS_wa(2*mmax*nbfgs+5*nbfgs+12*mmax*mmax+12*mmax))
    allocate(BFGS_iwa(3*nbfgs))

    allocate(lowest_vtrial(nfftf, dtset%nspden))

    allocate(vhartr_old(nfftf))

    BFGS_nbd(:) = 0

    run_params%oep_iter = 0

    allocate(new_occ(dtset%mband*dtset%nkpt*dtset%nsppol))
    allocate(new_eigen(dtset%mband*dtset%nkpt*dtset%nsppol))

    allocate(doccde(dtset%mband*dtset%nkpt*dtset%nsppol))

    allocate(run_params%dEdU(nfftf,dtset%nspden))

    call allocate_uxc(nfftf, dtset%nkpt, dtset%nband, &
         dtset%nsppol, dtset%nspden, run_params)

    allocate( &
         run_params%orb_shift_vxc(nfftf, dtset%nspden), &
          run_params%lda_pot(nfftf, dtset%nspden))

    if (run_params%num_gradient_n .eq. -1) &
         & run_params%num_gradient_n = dtset%nfft

    dtset%tolwfr   = 1e-10
    dtset%nnsclo   = 50

    !=========================================================
    ! BFGS optimization loop
    !=========================================================

    ! initialize the LBFGSB library
    task = 'START'
    call setulb(nbfgs,mmax,vtrial,BFGS_l,BFGS_u,BFGS_nbd,&
         & etotal,run_params%dEdu,0._DP,0._DP, &
         & BFGS_wa,BFGS_iwa,task,BFGS_print,BFGS_csave, &
         & BFGS_lsave,BFGS_isave,BFGS_dsave)

    if (run_params%diis_threshold > 0.d0) then
       diis_brain%problemDim = nbfgs
!       diis_brain%problemDim = 10

       diis_brain%maxPercent = run_params%diis_maxPercent
       diis_brain%maxExtraPoints = run_params%diis_maxExtraPoints
       diis_brain%extraPoints = run_params%diis_ExtraPoints

       do_DIIS = .TRUE.
      
       write(2244,*) "Initializing DIIS, size ", diis_brain%problemDim

       call initializeDIIS(diis_brain)
    endif


    if (run_params%num_gradient_n > 0) then
       call get_numeric_gradient(afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop,&
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit, resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params, susmat,symrec,taug,taur,tollist,&
       & ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)
    endif

    !=========================================================
    ! BFGS optimization loop
    !=========================================================

    do

       vhartr_old = vhartr

       call get_vtrial_energy(afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop,&
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit, resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params, susmat,symrec,taug,taur,tollist,&
       & ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)
              ! 1. Set up orbital shift equations:

       do i=1,nfftf
          write(843,*) vhartr_old(i), vhartr(i)-vhartr_old(i)
       enddo

       write(843,*) ' '
       write(843,*) ' '

       ! Reserve space for Orbital shifts 
       ! vxc (the KS potential)
       ! uxc (the real energy derivative dEXC/dphir)
       !  and the orb_shift vector g

       do ispden=1,min(dtset%nspden,2)
          run_params%orb_shift_vxc(:,ispden) = vtrial(:,ispden) -  &
                   vhartr(:) - vpsp(:)
          ! just the current trial XC potential 
       enddo
       if (dtset%nspden==4) run_params%orb_shift_vxc(:,3:4)=vtrial(:,3:4)

90    continue

       !! second call to vtorho, trap code in orbital_shifts
       !! ( vtorho -> vtowfk -> cgwf_orbital_shift )
       !! calculate the orbital shifts there
       
       !! reset the gradient
       
       run_params%dEdU = 0.d0

       write(2244,*) "Calling vtorho second time"
       
       istep_os = 1
       
       call vtorho_orbitalshifts(afford,atindx,atindx1, &
            &     cg,cpus,dbl_nnsclo,&
            &     dielop,dielstrt,dphase, dtfil,dtset,&
            &     eigen,electronpositron,energies,etotal,gbound_diel,&
            &     gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
            &     istep_os,istep_os,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,&
            &     mpi_enreg,my_natom,dtset%natom,nattyp,nfftf,nfftdiel,&
            &     ngfftdiel,nkxc,npwarr,npwdiel,res2,&
            &     psps%ntypat,nvresid,occ,computed_forces,optres,paw_dmft,&
            &     paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,phnons,&
            &     phnonsdiel,ph1d,ph1ddiel,psps,pwind,pwind_alloc,pwnsfac,&
            &     resid,residm,rhog,rhor,rmet,rprimd,run_params,&
            &     symrec,taug,taur,ucvol,wffnew,wffnow,&
            &     vhartr, vhartr_old, vtrial,wvl,xred,&
            &     ylm,ylmgr,ylmdiel)


#if defined HAVE_MPI
       if (mpi_enreg%nproc > 1) then
          call xsum_mpi(run_params%dEdu,spaceComm,ierr)
       end if
#endif

       do iy = 1,10
          write(231,*) vtrial(iy,1)
       enddo
       write(231,*) ' '
       write(231,*) ' '

       maycalldiis = .FALSE.

       if (etotal < lowest_e) then
          lowest_e = etotal
          lowest_vtrial = vtrial
          if (do_DIIS) then
             maycalldiis = .TRUE.
             call notifyDIIS(diis_brain, vtrial(:,1))
          endif
       endif
       
       ! correct for spin density

       if (dtset%nspden == 2) then

          deps = run_params%spindeps

!******
          n = dtset%nkpt * dtset%mband

          new_eigen(1:n) = eigen(1:n) + deps
          new_eigen(n+1:2*n) = eigen(n+1:2*n) - deps

          ! use new_occ as placeholder for the changed occupations
          call newocc(doccde,new_eigen,entropy_shift,e_fermie_shift,&
               dtset%spinmagntarget,dtset%mband,dtset%nband,dtset%nelect,&
               dtset%nkpt,dtset%nspinor,dtset%nsppol,new_occ,&
               dtset%occopt,dtset%prtvol,dtset%stmbias,dtset%tphysel,&
               dtset%tsmear,dtset%wtk)

          icount = 0
          dEdeps2 = 0.d0
          do iy=1,dtset%nsppol
             do iz=1,dtset%nkpt
                do ix=1,dtset%mband
                   dEdeps2 = dEdeps2 + (new_occ(icount)*&
                        new_eigen(icount)-&
                        occ(icount)*eigen(icount))*dtset%wtk(iz)
                   icount = icount+1
                enddo
             enddo
          enddo
          
          dEdeps2 = -dEdeps2 / (2.d0 * deps)

!******

          dEdeps = dEdeps2
!******
!
!             vtrial(:,1) = vtrial(:,1) + deps
!             vtrial(:,2) = vtrial(:,2) - deps
!
!             energies_shift%e_corepsp = energies%e_corepsp
!             energies_shift%e_ewald = energies%e_ewald
!
!             call get_vtrial_energy(afford, atindx, atindx1, &
!                  & dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop,&
!                  & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
!                  & electronpositron,energies_shift,etotal_shift,gbound_diel,&
!                  & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
!                  & istep_mix,kg,kg_diel,kxc,lmax_diel,mgfftdiel,mpi_enreg,&
!                  & nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
!                  & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
!                  & computed_forces,optres,paw_dmft,paw_ij,pawang,&
!                  & pawfgr,pawfgrtab,&
!                  & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
!                  & pwind,pwind_alloc,pwnsfac,quit, resid,residm,rhog,rhor,&
!                  & rmet,rprimd, susmat,symrec,taug,taur,tollist,ucvol,&
!                  & wffnew,wffnow,&
!                  & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial)
!             
!             vtrial(:,1) = vtrial(:,1) - deps
!             vtrial(:,2) = vtrial(:,2) + deps
!
!             dEdeps = (etotal_shift - etotal)/deps             
!
!             write(785,*) dedeps,dedeps2

!*******

          run_params%dEdu(:,1) = run_params%dEdu(:,1) + &
               & dEdeps * run_params%spingradmult / dble(dtset%nfft)
          run_params%dEdu(:,2) = run_params%dEdu(:,2) - &
               & dEdeps * run_params%spingradmult / dble(dtset%nfft)
          
          dEdu_mean_up = sum(run_params%dEdu(:,1))/dtset%nfft
          dEdu_mean_down = sum(run_params%dEdu(:,2))/dtset%nfft
          
       endif
       
       dEdu_mean = sum(run_params%dEdu(:,:))
       run_params%dEdu(:,:) = run_params%dEdu(:,:) - dEdu_mean
       write(2244,*) "Gradient max/min is ",&
            maxval(run_params%dEdu), minval(run_params%dEdu)

       if (run_params%num_gradient_n > 0) then

          open(file="analytic_gradient.dat", unit=723)

          do i=1, dtset%nfft   !run_params%num_gradient_n
             write(723,*) run_params%dedu(i,1)
          enddo

          if (dtset%nsppol .eq. 2) then

             do i=1, dtset%nfft   !run_params%num_gradient_n
                write(723,*) run_params%dedu(i,2)
             enddo
          endif

          close(723)

          if (run_params%niterategrad > 0) then

             open(file="iter_analyticgradient.dat", unit=723)       

             write(723,*) etotal

             do i=1,run_params%niterategrad
                    
                vtrial = vtrial + run_params%dEdU * run_params%diterstepgrad

                call get_vtrial_energy(afford, atindx, &
                  & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop,&
                  & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
                  & electronpositron,energies,etotal,gbound_diel,&
                  & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
                  & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
                  & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
                  & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
                  & computed_forces,optres,&
                  & paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
                  & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
                  & pwind,pwind_alloc,pwnsfac,quit, resid,residm,rhog,rhor,&
                  & rmet,rprimd, run_params, susmat,symrec,taug,taur,tollist,&
                  & ucvol,wffnew,wffnow,&
                  & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)

                write(723,*) etotal

             enddo

             close(723)

          endif

          stop

       endif


100    continue

!       dgrad = dsqrt(dot_product(run_params%dEdU,run_params%dEdU)*ucvol/nfftf)

       call setulb(nbfgs,mmax,vtrial,BFGS_l,BFGS_u,BFGS_nbd,&
            & etotal,run_params%dEdu,0._DP,0._DP, &
            & BFGS_wa,BFGS_iwa,task,BFGS_print,BFGS_csave, &
            & BFGS_lsave,BFGS_isave,BFGS_dsave)

       write(2244,*) "setulb returns task ", task(1:5)

       ! check BFGS flag
       if (task(1:5) .EQ. 'NEW_X') then
          write(2244,'(a)')'BFGS: TASK=NEW_X. New ITER ------- '
          ! dump files to disk
          ! call out_put_pots(nfftf,vtrial(:,1),vpsp,vhartr)
          goto 100           ! goto BFGS subroutine

       elseif (task(1:2) .EQ. 'FG') then
          write(2244,'(a)')'BFGS: TASK=FG, try new vtrial(r)'

          delta_e = etotal_old - etotal

          write(2244,*) "Delta e: ", delta_e, " comp ", run_params%accuracy_goal

          if (active_DIIS .and. abs(delta_e) > run_params%diis_alarmthreshold) then
             write(2244,*) "Delta e suspiciously increased with DIIS"
             write(2244,*) "reset DIIS state"
             call resetDIIS(diis_brain)
             active_DIIS = .FALSE.
          endif

          if (maycalldiis .and. abs(delta_e) < run_params%diis_threshold) then
             active_DIIS = .TRUE.
             write(2244,*) "Calling DIIS"

                ! we do not need to check whether DIIS has enough points,
                ! it will do that for us
             call getDIISPrediction(diis_brain, vtrial(:,1))

          else
             active_DIIS = .FALSE.
          endif
          maycalldiis = .FALSE.

          if (abs(delta_e) < run_params%accuracy_goal) then
             quitcounter = quitcounter + 1
             if (quitcounter > 2) then
                if (lowest_e < etotal) then
                   vtrial = lowest_vtrial
                   quitcounter = 0
                else
                   write(2244,*) "Convergence reached"
                   exit
                endif
             endif
          else
             quitcounter = 0
          endif

          etotal_old = etotal

       else
          ! unknown task, exit the program
          write(2244,'(a)')'what is task?'
          write(2244,'(a)')'task from bfgs => ',task(1:10)
          write(2244,'(a)')'abnormal task from BFGS, exit code'
          exit
       endif  ! BFGS test: task


       run_params%oep_iter = run_params%oep_iter + 1

       if (run_params%oep_iter>=run_params%max_oep_iter) exit

    enddo ! end of BFGS loop

    ! Compute the final energy with the converged potential

    dovtorho=.TRUE.

    if (realquit > 1) realquit=realquit - 1

    write(2244,*)''
    write(2244,*)'====================='
    WRITE(2244,*)' DIRECT OEP FINISHED'
    write(2244,*)'====================='
    write(2244,*)''

    ! free the XC storage
    call free_uxc(run_params)

    deallocate(doccde)
    deallocate(vhartr_old)

    deallocate(run_params%dEdU)
    deallocate(new_occ)
    deallocate(new_eigen)

    deallocate(run_params%orb_shift_vxc)
    deallocate(run_params%lda_pot)

    deallocate(BFGS_nbd)
    deallocate(BFGS_l)
    deallocate(BFGS_u)
    deallocate(BFGS_wa)
    deallocate(BFGS_iwa)

    deallocate(lowest_vtrial)

    if (do_diis) call cleanBrain(diis_brain)

  end subroutine lbfgsb_loop

  subroutine lbfgs_loop(afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, dielop, dielstrt,&
       & dovtorho, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,realquit,resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params,susmat,symrec,&
       & taug,taur,tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)

    implicit none

    type(hdr_type),intent(inout) :: hdr
    type(dataset_type) :: dtset
    type(MPI_type),intent(inout) :: mpi_enreg
    type(electronpositron_type),pointer :: electronpositron
    type(energies_type) :: energies

    integer my_natom

    real(dp) :: tollist(12)

    logical, intent(out) :: dovtorho
    integer, intent(inout) :: quit
    integer, intent(out) :: realquit

    type(datafiles_type),intent(in) :: dtfil
    type(efield_type),intent(inout) :: dtefield
    type(pseudopotential_type),intent(in) :: psps

    type(wffile_type),intent(inout) :: wffnew,wffnow
    type(wvl_data),intent(inout) :: wvl

    integer n3xccc
    real(dp) xccc3d(n3xccc)

    integer, intent(in) :: afford
    integer,intent(in) :: atindx(dtset%natom)
    integer, intent(in) :: atindx1(dtset%natom)
    integer, intent(in) :: dbl_nnsclo
    integer, intent(in) :: optres
    integer, intent(in) :: dielop
    integer, intent(in) :: dielstrt
    integer, intent(in) :: nkxc
    integer, intent(in) :: nattyp(psps%ntypat)

    integer, intent(in) :: ngfftdiel(18)
    integer, intent(in) :: mgfftdiel
    integer, intent(in) ::  gbound_diel(2*mgfftdiel+8,2)
    integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)

    integer, intent(in) :: computed_forces

    real*8, intent(in) :: gsqcut

    integer istep_OS

    integer nfftf

    real*8 nvresid(nfftf,dtset%nspden)
    real*8 compch_fft
    real*8 dphase(3)
    real*8 grnl(3*dtset%natom)
    real*8 kxc(nfftf,nkxc)
    real*8 nhat(nfftf, dtset%nspden)

    real*8, intent(in) :: gmet(3,3),gprimd(3,3)

    real(dp), intent(in) :: rprimd(3,3)
    real(dp), intent(in) :: cpus

    real(dp), intent(in) :: ucvol

    integer, intent(in) :: pwind_alloc
    integer, intent(in) :: pwind(pwind_alloc,2,3)
    integer, intent(in) :: symrec(3,3,dtset%nsym)

    real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
    real(dp), intent(inout) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(dp) residm, vxcavg

    real(dp), pointer :: rhog(:,:),rhor(:,:)
    real(dp), pointer :: taug(:,:),taur(:,:)

    real(dp), intent(in) :: rmet(3,3)

    real(dp) susmat(2,npwdiel*afford,dtset%nspden,&
         &npwdiel,dtset%nspden)

    real(dp), intent(inout) :: xred(3,dtset%natom)

    ! PAM set

    type(paw_dmft_type), intent(inout) :: paw_dmft
    type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
    type(pawang_type),intent(in) :: pawang
    type(pawfgr_type),intent(inout) :: pawfgr
    type(pawfgrtab_type) pawfgrtab(my_natom*psps%usepaw)
    type(pawrhoij_type), intent(inout) :: pawrhoij(my_natom*psps%usepaw)
    type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
    
    ! ==============

    type(energies_type) :: energies_shift

    integer nfftdiel

    real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
    real*8, intent(in) :: ph1ddiel(2,3*(2*mgfftdiel+1)*my_natom*psps%usepaw)

    integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(inout) :: istep_mix
    integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
    integer, intent(in) :: npwdiel
    integer, intent(in) :: lmax_diel

    integer, intent(in) :: npwarr(dtset%nkpt)
    integer, intent(in) :: kg_diel(3,npwdiel)

    real(dp) res2
    real(dp) dummy_nhatgr(1,1,1) ! dummy nhatgr

    real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

    integer, intent(in) :: mcg
    real(dp), intent(inout) :: cg(2,mcg)

    real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,&
         & psps%mpsang*psps%mpsang*psps%useylm)
    real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,&
         & 3,psps%mpsang*psps%mpsang*psps%useylm)
    real(dp), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
    real(dp) vxctau(nfftf,dtset%nspden*dtset%usekden,4)

    real(dp) vhartr(nfftf),vpsp(nfftf)

    real(dp), allocatable :: vhartr_old(:)

    real (KIND=8) vtrial(nfftf, dtset%nspden)
    real*8 etotal
    real*8 dEdu_mean, dEdu_mean_up, dEdu_mean_down

    integer :: ix, iy ,iz, icount

    integer ispden

    real(kind=dp) dgrad

    real(dp) strsxc(6) ! xc stress tensor dummy

    !! BFGS related
    integer              :: nbfgs ! Problem size, nfft * nspden
    integer              :: mmax

    real*8 doti

    integer indx1,indx2,indx3
    integer ktmp(3)

    integer spin_converge_stage

    integer quitcounter

    real*8 etotal_old, delta_e, etotal_shift

    integer n
    real(kind=dp),allocatable :: doccde(:)
    real*8 e_fermie_shift, entropy_shift, deps, dedeps2, dedeps

    real(dp), allocatable :: new_occ(:)
    real(dp), allocatable :: new_eigen(:)

!!!!!!!!!!!!!!!!!

    integer ierr

    integer spaceComm

    character*10 filename

    type (run_parameter_type) run_params
        
    ! FROM HERE ON DOWNWARDS THINGS FOR THE L-BFGS
    real(kind=dp),allocatable::diag(:)
    real(kind=dp),allocatable::work(:)
    integer::lbfgs_flag
    ! TILL HERE

!!!!!!!!!!!!!!!!!

!MPI communicator
    spaceComm=mpi_enreg%comm_cell

    dtset%tolwfr   = 1e-10
    dtset%nnsclo   = 50
    quitcounter = 0
    delta_e = 0.d0

    write(2244,*) "Using smearing of ",run_params%deltaT

    ! BFGS related ...
    !=================
    nbfgs = nfftf * dtset%nspden

    mmax = run_params%mmax_bfgs

    run_params%oep_iter = 0

    allocate(vhartr_old(nfftf))

    allocate(new_occ(dtset%mband*dtset%nkpt*dtset%nsppol))
    allocate(new_eigen(dtset%mband*dtset%nkpt*dtset%nsppol))

    allocate(doccde(dtset%mband*dtset%nkpt*dtset%nsppol))

    allocate(run_params%dEdU(nfftf,dtset%nspden))

    call allocate_uxc(nfftf, dtset%nkpt, dtset%nband, &
         dtset%nsppol, dtset%nspden, run_params)

    allocate( &
         run_params%orb_shift_vxc(nfftf, dtset%nspden), & 
          run_params%lda_pot(nfftf, dtset%nspden))
          

    !=========================================================
    ! BFGS optimization loop
    !=========================================================

    dtset%tolwfr   = 1e-10
    dtset%nnsclo   = 50

    bfgsloop: do 

       vhartr_old = vhartr
      
       call get_vtrial_energy(afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop,&
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit, resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params, susmat,symrec,taug,taur,tollist,&
       & ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)

       ! 1. Set up orbital shift equations:

       ! Reserve space for Orbital shifts 
       ! vxc (the KS potential)
       ! uxc (the real energy derivative dEXC/dphir)
       !  and the orb_shift vector g

       do ispden=1,min(dtset%nspden,2)
          run_params%orb_shift_vxc(:,ispden) = vtrial(:,ispden) -  &
                   vhartr(:) - vpsp(:)
          ! just the current trial XC potential 
       enddo
       if (dtset%nspden==4) run_params%orb_shift_vxc(:,3:4)=vtrial(:,3:4)      

90    continue           

       !! second call to vtorho, trap code in orbital_shifts
       !! ( vtorho -> vtowfk -> cgwf_orbital_shift )
       !! calculate the orbital shifts there
       
       !! reset the gradient
       
       run_params%dEdU = 0.d0

       write(2244,*) "Calling vtorho second time"

       istep_os = 1

       call vtorho_orbitalshifts(afford,atindx,atindx1, &
            &     cg,cpus,dbl_nnsclo,&
            &     dielop,dielstrt,dphase, dtfil,dtset,&
            &     eigen,electronpositron,energies,etotal,gbound_diel,&
            &     gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
            &     istep_os,istep_os,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,&
            &     mpi_enreg,my_natom,dtset%natom,nattyp,nfftf,nfftdiel,&
            &     ngfftdiel,nkxc,npwarr,npwdiel,res2,&
            &     psps%ntypat,nvresid,occ,computed_forces,optres,paw_dmft,&
            &     paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,phnons,&
            &     phnonsdiel,ph1d,ph1ddiel,psps,pwind,pwind_alloc,pwnsfac,&
            &     resid,residm,rhog,rhor,rmet,rprimd,run_params,&
            &     symrec,taug,taur,ucvol,wffnew,wffnow,&
            &     vhartr, vhartr_old, vtrial,wvl,xred,&
            &     ylm,ylmgr,ylmdiel)


#if defined HAVE_MPI
       if (mpi_enreg%nproc > 1) then
          call xsum_mpi(run_params%dEdu,spaceComm,ierr)
       end if
#endif
       ! correct for spin density

       if (dtset%nspden == 2) then             

          deps = 0.001d0

          !******
          n = dtset%nkpt * dtset%mband

          new_eigen(1:n) = eigen(1:n) + deps
          new_eigen(n+1:2*n) = eigen(n+1:2*n) - deps

          ! use new_occ as placeholder for the changed occupations
          call newocc(doccde,new_eigen,entropy_shift,e_fermie_shift,&
               dtset%spinmagntarget,dtset%mband,dtset%nband,dtset%nelect,&
               dtset%nkpt,dtset%nspinor,dtset%nsppol,new_occ,&
               dtset%occopt,dtset%prtvol,dtset%stmbias,dtset%tphysel,&
               dtset%tsmear,dtset%wtk)

          icount = 0 
          dEdeps2 = 0.d0
          do iy=1,dtset%nsppol
             do iz=1,dtset%nkpt
                do ix=1,dtset%mband
                   dEdeps2 = dEdeps2 + (new_occ(icount)*&
                        new_eigen(icount)-&
                        occ(icount)*eigen(icount))*dtset%wtk(iz)
                   icount = icount+1
                enddo
             enddo
          enddo

          dEdeps2 = -dEdeps2 / (2.d0 * deps)

          !******

          dEdeps = dEdeps2

          !******
          !
          !             vtrial(:,1) = vtrial(:,1) + deps
          !             vtrial(:,2) = vtrial(:,2) - deps
          !
          !             energies_shift%e_corepsp = energies%e_corepsp
          !             energies_shift%e_ewald = energies%e_ewald
          !
          !             call get_vtrial_energy(afford, atindx, atindx1, &
          !                  & dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop,&
          !                  & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
          !                  & electronpositron,energies_shift,etotal_shift,gbound_diel,&
          !                  & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
          !                  & istep_mix,kg,kg_diel,kxc,lmax_diel,mgfftdiel,mpi_enreg,&
          !                  & nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
          !                  & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&


          !              WHATS WITH RES2 / NRES2 ??


          !                  & computed_forces,optres,paw_dmft,paw_ij,pawang,&
          !                  & pawfgr,pawfgrtab,&
          !                  & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
          !                  & pwind,pwind_alloc,pwnsfac,quit, resid,residm,rhog,rhor,&
          !                  & rmet,rprimd, susmat,symrec,taug,taur,tollist,ucvol,&
          !                  & wffnew,wffnow,&
          !                  & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial)
          !             
          !             vtrial(:,1) = vtrial(:,1) - deps
          !             vtrial(:,2) = vtrial(:,2) + deps
          !
          !             dEdeps = (etotal_shift - etotal)/deps             
          !
          !             write(785,*) dedeps,dedeps2

          !*******

          run_params%dEdu(:,1) = run_params%dEdu(:,1) + dEdeps * 0.1
          run_params%dEdu(:,2) = run_params%dEdu(:,2) - dEdeps * 0.1

          dEdu_mean_up = sum(run_params%dEdu(:,1))/dtset%nfft
          dEdu_mean_down = sum(run_params%dEdu(:,2))/dtset%nfft

       endif

       dEdu_mean = sum(run_params%dEdu(:,:))
       run_params%dEdu(:,:) = run_params%dEdu(:,:) - dEdu_mean
       write(2244,*) "Gradient max/min is ",&
            maxval(run_params%dEdu), minval(run_params%dEdu)


!100    continue           

!       dgrad = dsqrt(dot_product(run_params%dEdU,run_params%dEdU)*ucvol/nfftf)

       call lbfgs(nbfgs,mmax,vtrial,etotal,run_params%dEdu,&
            .false.,diag,-1,run_params%accuracy_goal,1E-20,work,lbfgs_flag)
            

       write(2244,*) "bfgs returns task ", lbfgs_flag

       ! check BFGS flag
       if(lbfgs_flag .eq. -1) then
          write(2244,*) "LBFGS FAILS",lbfgs_flag
          stop
       else if (lbfgs_flag .eq. -2) then
          write(2244,*) "LBFGS: DIAG IS NOT POSITIVE",lbfgs_flag
          stop
       else if (lbfgs_flag .eq. -3) then
          write(2244,*) "LBFGS IMPROPER INPUT",lbfgs_flag
          stop
       else if (lbfgs_flag .eq. 0) then
           write(2244,*) "Convergence reached"
           exit bfgsloop
       else if (lbfgs_flag .eq. 1) then
          ! function and gradient evaluated
          etotal_old = etotal
          write(2244,'(a)')'BFGS: TASK=FG, try new vtrial(r)'

          delta_e = etotal_old - etotal

          write(2244,*) "Delta e: ", delta_e, " comp ", run_params%accuracy_goal


          !if (abs(delta_e) < run_params%accuracy_goal) then
          !   quitcounter = quitcounter + 1
          !   if (quitcounter > 2) then
          !      write(2244,*) "Convergence reached"
          !      exit bfgsloop
          !   endif
          !else
          !   quitcounter = 0
          !endif
       else if (lbfgs_flag .eq. 2) then
          write(2244,*) "LBFGS ASKS TO PROVIDE DIAG WHICH WE CAN'T ", lbfgs_flag
          stop   
       else
          write(2244,*) "LBFGS RETURNS UNKNOWN FLAG ", lbfgs_flag
          stop
       endif
              
       run_params%oep_iter = run_params%oep_iter + 1

       if (run_params%oep_iter>=run_params%max_oep_iter) exit bfgsloop

    enddo bfgsloop ! end of BFGS loop

    ! Compute the final energy with the converged potential

    dovtorho=.TRUE.

    if (realquit > 1) realquit=realquit - 1
    
    write(2244,*)''
    write(2244,*)'====================='
    WRITE(2244,*)' DIRECT OEP FINISHED'
    write(2244,*)'====================='
    write(2244,*)''
    
    ! free the XC storage
    call free_uxc(run_params)

    deallocate(doccde)
    deallocate(vhartr_old)
    
    deallocate(run_params%dEdU)
    deallocate(new_occ)
    deallocate(new_eigen)

    deallocate(run_params%orb_shift_vxc)
    deallocate(run_params%lda_pot)

    deallocate(diag)
    deallocate(work)

  end subroutine lbfgs_loop


  !=============================================
  ! Chen: compute dExc/dphir for LDA (PW92)
  !
  !  by Chen Huang  (May/21/2011)
  !=============================================
  subroutine xcdev_ldapw92(n1,n2,n3,nfft,rhor,ucvol,xcpot,xcenergy)
    implicit none 

    integer,intent(in)       :: nfft,n1,n2,n3
    real(kind=8),intent(in)  :: rhor(nfft),ucvol
    real(kind=8),intent(out) :: xcpot(nfft),xcenergy
    real(kind=8)             :: rhor3d(n1,n2,n3),xcpot3d(n1,n2,n3,2)

    ! test
    if (n1*n2*n3/=nfft) then 
       write(2244,*)'xcdev_ldapw92: n1*n2*n3/=nfft'
       stop
    endif

    rhor3d = reshape(rhor,(/n1,n2,n3/))

    call LSDAPW92(n1,n2,n3,rhor3d,xcpot3d,xcenergy,ucvol/dble(nfft))

    xcpot = reshape(xcpot3d(:,:,:,1),(/nfft/))

  end subroutine xcdev_ldapw92

  !================================================================
  ! PW92 code from PROFESS, compute the xc of Perdew and Wang 1992
  !
  !  by Chen Huang (May/21/2011)
  !================================================================
  SUBROUTINE LSDAPW92(n1,n2,n3,rhor,LDAPotential,LDAEnergy,dVol)
    IMPLICIT NONE

    !>> EXTERNAL VARIABLES <<!
    integer,parameter        :: dp=8
    integer,intent(in)       :: n1,n2,n3 
    real(kind=dp),intent(in) :: rhor(n1,n2,n3),dVol
    real(kind=dp),intent(out):: LDAPotential(n1,n2,n3,2)
    real(kind=dp),intent(out):: LDAEnergy

    !>> INTERNAL VARIABLES <<!
    real(dp) :: rho(n1,n2,n3,2) ! consider spin

    real(kind=dp), parameter :: &
         pi = 3.1415926d0, &
         p75vpi = 0.75d0/pi, &
         ax = -0.7385588d0, &
         fzz = 1.709921d0, &        
         gamma = 0.5198421d0, &
         one =1.d0, &
         two = 2.d0, &
         three = 3.d0, &
         four = 4.d0

    REAL(kind=DP) :: &
         rs, &                              ! (3/(4.pi.rho))^1/3
         zet, &                             ! Spin-polarization
         ex(2), vx(2),exc, &                           ! exchange energy
         ec, vc(2), &                              ! correlation energy
         eu,f,z4,ecrs,eurs,eczet,comm,fz,d, &  ! work variables
         ep,eprs,alfrsm,alfm,ac2,third         ! work variables

    INTEGER :: &
         ix, iy, iz                           ! Dummy counters

    !! we ignore spin now
    !! 

    rho(:,:,:,1) = rhor/2.d0
    rho(:,:,:,2) = rhor/2.d0

    third = one/three
    ac2 = four/three
    exc= 0.d0
    DO iz=1, n3
       DO iy=1, n2
          DO ix=1, n1
             ! exchange
             d=two*rho(ix,iy,iz,1)
             ex(1)= ax*d**third
             vx(1)=ex(1)*ac2
             d=two*rho(ix,iy,iz,2)
             ex(2)= ax*d**third
             vx(2)=ex(2)*ac2
             ! local correlation
             d=sum(rho(ix,iy,iz,:))
             rs=(p75vpi/d)**third
             call spn_gcor(0.0310907d0,0.21370d0,7.5957d0,3.5876d0, &
                  1.6382d0,0.49294d0,1.00d0,rs,eu,eurs)
             call spn_gcor(0.01554535d0,0.20548d0,14.1189d0,6.1977d0, &
                  3.3662d0,0.62517d0,1.00D0,rs,ep,eprs)
             call spn_gcor(0.0168869d0,0.11125d0,10.357d0,3.6231d0, &
                  0.88026d0,0.49671d0,1.00d0,rs,alfm,alfrsm)
             zet=(rho(ix,iy,iz,1)-rho(ix,iy,iz,2))/d
             f = ((one+zet)**ac2+(one-zet)**ac2-two)/gamma
             z4 = zet**4
             ec = eu*(one-f*z4)+ep*f*z4-alfm*f*(one-z4)/fzz
             ecrs = eurs*(one-f*z4)+eprs*f*z4-alfrsm*f*(one-z4)/fzz
             fz = ac2*((one+zet)**third-(one-zet)**third)/gamma
             eczet = four*(zet**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu &
                  -(one-z4)*alfm/fzz)
             comm = ec -rs*ecrs/three-zet*eczet
             vc(1) = comm + eczet
             vc(2) = comm - eczet

             !           vxc(i,1)=two*(vx(1)+vc(1))
             !           vxc(i,2)=two*(vx(2)+vc(2))
             ec=(rho(ix,iy,iz,1)*ex(1)+rho(ix,iy,iz,2)*ex(2)+ec*d)
             exc = exc + ec
             LDAPotential(ix,iy,iz,1)=(vx(1)+vc(1)) !- 0.5d0 * emf
             LDAPotential(ix,iy,iz,2)=(vx(2)+vc(2)) !+ 0.5d0 * emf
          ENDDO
       ENDDO
    ENDDO
    LDAEnergy = dVol*exc

    write(6,'(a,f)')'leave LSDAPW92(), e_lda=',LDAEnergy
  END SUBROUTINE LSDAPW92

  !---------------------------------------------------------------
  !
  ! This subroutine computes the local correlation energy and
  ! potential for the Perdew-Wang exchange-correlation scheme.
  !
  !---------------------------------------------------------------
  SUBROUTINE spn_gcor(a,a1,b1,b2,b3,b4,p,rs,gg,ggrs)

    implicit none
    !
    ! Input/Output variables:
    !
    integer,parameter :: dp=8
    real(dp), intent(in) :: a,a1,b1,b2,b3,b4,p,rs
    real(dp), intent(out) :: gg,ggrs
    !
    ! Work variables:
    !
    real(dp) :: p1,q0,rsp,q1,q2,q3,rs12,rs32,two, one,three
    !---------------------------------------------------------------
    one = 1.d0
    two = 2.d0
    three = 3.d0
    p1 = p + 1.d0
    q0 = -two*a*(one+a1*rs)
    rs12 = sqrt(rs)
    rs32 = rs12**3
    rsp = rs**p
    q1 = two*a*(b1*rs12+b2*rs+b3*rs32+b4*rs*rsp)
    q2 = log(one+one/q1)
    gg = q0*q2
    q3 = a*(b1/rs12+two*b2+three*b3*rs12+two*b4*p1*rsp)
    ggrs = -two*a*a1*q2-q0*q3/(q1**2+q1)

  END SUBROUTINE spn_gcor


  !==============================================
  ! chen: code
  ! display the spectrum 
  !==============================================
  subroutine disp_eigenvalues(nband,ei,occ,oep_unit)

    implicit none

    integer,intent(in) :: nband,oep_unit
    real(kind=8),intent(in) :: ei(nband)
    real(kind=8),intent(in) :: occ(nband)

    integer :: ii

    write(oep_unit,*)' '
    write(oep_unit,*)'---------------------------------------------------- '
    write(oep_unit,*)' band     eigenvalue        occ '
    do ii=1,nband
       if ( abs(occ(ii))>1e-10 ) then 
          write(oep_unit,'(a,i3,a,es16.8,a,es12.4)')'  ',ii,'  ',ei(ii),' ',occ(ii)
       endif
    enddo
    write(oep_unit,*)'---------------------------------------------------- '
    write(oep_unit,*)' '

  end subroutine disp_eigenvalues

  !!
  !! output various potentials to disk for analysis
  !!  
  !!  Chen Huang , June/1/2011
  !!
  subroutine out_put_pots(nfft,vtrial,vpsp,vhart)
    implicit none

    integer,intent(in)      :: nfft
    real(kind=8),intent(in) :: vtrial(nfft),vpsp(nfft),vhart(nfft)

    integer :: ii

    open(file='dump_vtrial',unit=199,action='write')
    open(file='dump_vpsp'  ,unit=200,action='write')
    open(file='dump_vxc'   ,unit=201,action='write')

    ! write(2244,'(a,i)')'(out_put_pots) dumping vtrial, vpsp, and vxc. nfft=',nfft

    do ii=1,nfft
       write(199,*)vtrial(ii)
       write(200,*)vpsp(ii)
       write(201,*)vtrial(ii)-vpsp(ii)-vhart(ii)
    enddo

    close(199) 
    close(200)
    close(201)
    write(2244,'(a)')'(out_put_pots) dump_vxc and dump_vtrial files dumped.'

  end subroutine out_put_pots


  !!
  !! load parameters from param.in file
  !!
  subroutine  loadpara( run_params )

    implicit none 

    type (run_parameter_type) run_params

    type (T_Input_data) id

    type (diisBrain) db

    integer indx
    
    call open_input("param.in", "OEP.log",id)
    
    run_params%precondition = check_bool(id, "precond", indx)
    ! 1: run conventional abinit scf before OEP
    ! 0: do not, directly start for OEP from initial guess

    run_params%restart      = check_bool(id, "restart", indx)
    ! 1: real space EXX code will be performed every iter. for virial theorem
    ! 0: every iter, EXX is evaluated with recip code (much faster)
    run_params%check_virial = check_bool(id, "chkvirl", indx)

    ! xc_type:
    ! 1  : LDA(PW92)
    ! 2  : XX nonperiodic
    ! 3  : XX periodic
    indx=-1
    call demand_i_param(id, "xc_type", indx, run_params%xc_type)

    ! realspc:
    ! 0: do chen_EXX in recip space
    ! 1: do chen_EXX in realspace
    ! 2: do the EXX energy in recip space, but potential in real space
    indx=-1
    call demand_i_param(id, "realspc", indx, run_params%do_realspace)

    !Scaling of analytic gradient
    indx=-1
    if (.NOT.check_r_param(id, "gradscl", indx, run_params%gradscal)) then
       run_params%gradscal=1.d0
    endif

    ! maximum number of oep iterations
    indx=-1
    call demand_i_param(id, "maxNoep", indx, run_params%max_oep_iter)

    ! number of Hessian corrections in LBFGS
    indx=-1
    if (.NOT.check_i_param(id, "mmaxBFG", indx, run_params%mmax_bfgs)) then
       run_params%mmax_bfgs=7
    endif
    
    ! which optimization algorithm to use
    ! 1: L-BFGS
    ! 2: ThreeCG (Powell restart)
    ! 3: ThreeCG (angle restart)
    ! 4: L-BFGS-B
    indx=-1!    
    if (.NOT.check_i_param(id, "optalgo", indx, run_params%optalgo)) then
       run_params%optalgo=4
    endif

    indx=-1
    call demand_r_param(id, "E-mstep", indx, run_params%finite_E)

    if (.not.check_i_param(id, "NumGr_N", indx, run_params%num_gradient_n)) then
       run_params%num_gradient_n = 0
       run_params%num_gradient_delta = 0.d0
    else
       indx=-1
       call demand_r_param(id, "NumGr_D", indx, run_params%num_gradient_delta)
    endif

    if (.NOT.check_r_param(id, "vhrtrpf", indx, run_params%vhrtr_prefac)) then
       run_params%vhrtr_prefac=0.d0
    endif
    
    if (.NOT.check_r_param(id, "ds_thrs", indx, run_params%diis_threshold)) then
       run_params%diis_threshold=0.d0
    endif

    if (.NOT.check_r_param(id, "ds_alrm", indx, run_params%diis_alarmthreshold)) then
       run_params%diis_alarmthreshold=run_params%diis_threshold * 2.d0
    endif

    if (.NOT.check_r_param(id, "ds_maxP", indx, run_params%diis_maxPercent)) then
       run_params%diis_maxPercent=db%maxPercent
    endif

    if (.NOT.check_i_param(id, "ds_extP", indx, run_params%diis_extraPoints)) then
       run_params%diis_extraPoints=db%extraPoints
    endif

    if (.NOT.check_i_param(id, "ds_mxEP", indx, run_params%diis_maxExtraPoints)) then
       run_params%diis_maxExtraPoints=db%maxExtraPoints
    endif

    if (.not.check_i_param(id, "GrIterN", indx, run_params%niterategrad)) then
       run_params%niterategrad = 0
       run_params%diterstepgrad = 0.d0
    else
       indx=-1
       call demand_r_param(id, "GrIterD", indx, run_params%diterstepgrad)
    endif

    if (.not.check_r_param(id, "CGaccur", indx, run_params%CG_accuracy)) &
         run_params%CG_accuracy = 1d-6

    if (.not.check_r_param(id, "spndeps", indx, run_params%spindeps)) &
         run_params%spindeps = 0.001d0

    if (.not.check_r_param(id, "spngrad", indx, run_params%spingradmult)) &
         run_params%spingradmult = 0.2d0

    indx=-1
    call demand_r_param(id, "acuracy", indx, run_params%accuracy_goal)

    close(1001)

    if (run_params%max_oep_iter<0) then 
       write(2244,*)'max_oep_iter <=0 stop'
       stop
    endif
    

    call free_ID(id)

  end subroutine loadpara

  subroutine display_info(run_params)

    implicit none

    type (run_parameter_type), intent(in):: run_params

    write(2244,'(a)')'========================================================'
    WRITE(2244,'(a)')'                       fast OEP          '
    WRITE(2244,'(a)')'                                       '          
    WRITE(2244,'(a)')'                     (VERSION  4.1)       '          
    WRITE(2244,'(a)')'                                       '          
    WRITE(2244,'(a)')'    REVISION HISTORY                          '
    WRITE(2244,'(a)')'     * incorporated parallelism *'
    WRITE(2244,'(a)')'     * transferred to newest abinit version *'
    WRITE(2244,'(a)')'     * periodic exact exchange *'
    WRITE(2244,'(a)')'     * includes DIIS *'
    WRITE(2244,'(a)')'                                 '
    WRITE(2244,'(a)')'     Florian Libisch (April/2013)          '
    WRITE(2244,'(a)')'                                 '
    WRITE(2244,'(a)')'     DIIS-solver by Johannes M. Dieterich (April/2013)  '
    WRITE(2244,'(a)')'                                 '
    write(2244,'(a)')'========================================================'

    write(2244,*) "optimization alg   : optalgo : ", run_params%optalgo
    write(2244,*) "restart            : restart : ", run_params%restart
    write(2244,*) "check_virial       : chkvirl : ", run_params%check_virial
    write(2244,*) "gradscal           : gradscl : ", run_params%gradscal
    write(2244,*) "vhrtr_prefac       : vhrtrpf : ", run_params%vhrtr_prefac
    write(2244,*) "xc_type            : xc_type : ", run_params%xc_type
    write(2244,*) "do_realspace       : realspc : ", run_params%do_realspace
    write(2244,*) "max_oep_iter       : maxNoep : ", run_params%max_oep_iter
    write(2244,*) "num_gradient_delta : NumGr_D : ", &
                                 & run_params%num_gradient_delta
    write(2244,*) "num_gradient_n     : NumGr_N : ", run_params%num_gradient_n
    write(2244,*) "mmax_bfgs          : mmaxBFG : ", run_params%mmax_bfgs
    write(2244,*) "niterategrad       : GrIterN : ", run_params%niterategrad
    write(2244,*) "diterstepgrad      : GrIterD : ", run_params%diterstepgrad
    write(2244,*) "spingradmult       : spngrad : ", run_params%spingradmult
    write(2244,*) "spindeps           : spndeps : ", run_params%spindeps
    write(2244,*) "CG_accuracy        : CGaccur : ", run_params%CG_accuracy
    write(2244,*) "finite_E           : E-mstep : ", run_params%finite_E
    write(2244,*) "accuracy_goal      : acuracy : ", run_params%accuracy_goal
    write(2244,*) "precondition       : precond : ", run_params%precondition

    write(2244,*) "diis_threshold     : ds_thrs : ", run_params%diis_threshold
    write(2244,*) "diis_alarmthreshold: ds_alrm : ", &
                                 & run_params%diis_alarmthreshold
    write(2244,*) "diis_maxPercent    : ds_maxP : ", run_params%diis_maxPercent
    write(2244,*) "diis_extraPoints   : ds_extP : ", run_params%diis_extraPoints
    write(2244,*) "diis_maxExtraPoints: ds_mxEP : ", &
                                 & run_params%diis_maxExtraPoints


  end subroutine display_info
  
end module Pot_Optimization
