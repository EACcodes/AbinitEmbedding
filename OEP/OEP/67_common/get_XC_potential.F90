module get_xc_potential

  use defs_basis
  use defs_datatypes
  use defs_abitypes
  use m_xmpi
  use fourierhlp
  use spencerxx

  type run_parameter_type

     logical restart
     logical check_virial     
     logical precondition

     real*8 diis_threshold
     real*8 diis_alarmthreshold
     real*8 diis_maxPercent
     integer diis_extraPoints
     integer diis_maxExtraPoints

     integer xc_type
     integer do_realspace
     integer oep_iter
     integer max_oep_iter

     integer mmax_bfgs
     
     integer optalgo

     integer num_gradient_n

     integer niterategrad

     real*8 gradscal

     real*8 diterstepgrad

     real*8 vhrtr_prefac

     real*8 num_gradient_delta

     real*8 spingradmult
     real*8 spindeps

     real*8 CG_accuracy

     real*8 e_fermie
     real*8 deltaT
     
     real*8 finite_E

     real*8 accuracy_goal 
     ! When to stop the BFGS iteration

     REAL(dp), allocatable :: orb_shift_vxc(:,:) ! KS XC potential 

     !! Energy gradient with respect to potential
     REAL(dp), allocatable :: dEdu(:,:)

     real*8, allocatable :: uxc_data(:,:,:,:)       ! uxc potential
     real*8, allocatable :: uxc_im_data(:,:,:,:)       ! uxc potential, imaginary part

     real*8, allocatable :: lda_pot(:,:)       ! xc potential used by abinit

     integer, allocatable :: nband(:)

     real(dp), allocatable :: wtk(:)

  end type run_parameter_type

contains

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


  subroutine chen_exx(mpi_enreg,npw,nfft,norb,gmet,gsqcut,rprimd, & 
       cellvol,occ,ngfftf,dtset,cg,mcg,kg,ex,dex_dphir,do_realspace)

    use defs_basis
    use defs_datatypes
    use defs_abitypes
    use defs_scftypes
    use defs_parameters

    implicit none

    !!  type(dataset_type),intent(in) :: dtset
    type(dataset_type) :: dtset

    type(MPI_type),intent(in) :: mpi_enreg

    integer , intent(in)  ::  & 
         do_realspace, &  ! flag for real space or recip space
                                ! 0: everything in recip space
                                ! 1: everything in real space
                                ! 2: only xc energy in real space, 
                                !    xc_pot in recip space

         mcg, &           ! max size of cg
         npw, &           ! number of plane waves
         nfft, &          ! number of point in real space
         norb, &          ! number of orbitals to be consider
         ngfftf(18)    

    integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)

    real(dp), intent(in) :: cg(2,mcg)

    real(dp),intent(in) ::  & 
         rprimd(3,3), & 
         gmet(3,3),  &                  ! metric of g-space
         gsqcut,     &                  ! cutoff of g^2
         cellvol, &                     ! cell volumn
         occ(norb)

    real(dp),intent(out) :: & 
         ex, &                          ! E^HF_x energy
         dex_dphir(nfft,norb)           ! output: dE^HF_x / d( phir_n )

    ! local vars 
    !===============
    real(dp), allocatable  :: tmpr(:)      ! temp vars in real space
    real(dp), allocatable  :: tmpr_out(:)  ! temp vars for output of PSolver

    real(dp) :: tmp_ehart,   & 
         tmp_qphon(3)=0.d0, &   ! no use
         mu = 0.d0              ! no use

    real (dp), allocatable :: phir(:,:)     ! orbital in real space

    integer  :: usepaw = 0, & 
         m,n,old_icoulomb, i

    integer nproc, me, comm, icg_shift, iband

    real*8 hg(3)

    !! function begins 
    !!==================

    !!
    !! only for closed shell system
    !! cannot be used for H atom or other open-shell systems
    !! in those cases, the EXX needs to be written carefully
    !!

    icg_shift=0
 
    allocate(phir(nfft, norb))

    do iband=1,norb
             
       call recip_to_real(mpi_enreg,dtset%paral_kgb,dtset%nkpt,&
            npw,1, dtset%istwfk(1),dtset%nspinor,&
            dtset%mgfft,npw,cellvol, dtset%ngfft,kg(:,1:npw), npw,&
            cg(:,icg_shift+1:icg_shift+npw*dtset%nspinor),&
            phir(:,iband))
          
       icg_shift=icg_shift + npw*dtset%nspinor
             
    enddo

    if(mpi_enreg%paral_kgb==1) then
       comm=mpi_enreg%comm_fft
       me=mpi_enreg%me_fft
       nproc=mpi_enreg%nproc_fft
    else
       comm=mpi_enreg%comm_cell
       nproc=xcomm_size(comm)
       me=xcomm_rank(comm)
    end if

    do i=1,3
       hg(i) = dsqrt(&
            &  rprimd(1,i)**2 + rprimd(2,i)**2 + rprimd(3,i)**2 &
            & )/dtset%ngfft(i)
    enddo
    write(2244,'(a)')'get in chen_exx()...'

    allocate(tmpr(nfft))      ! temp vars in real space
    allocate(tmpr_out(nfft))  ! temp vars for output of PSolver

    old_icoulomb   = dtset%icoulomb  ! back up icoulomb

    if (do_realspace==1) then 
       write(2244,'(a)')'(chen_exx) NOTICE: real space : xc_energy and xc_pot'
       dtset%icoulomb = 1
    else if (do_realspace==2) then 
       write(2244,'(a)')'(chen_exx) NOTICE: real space : xc_energy, recip space: xc_pot'
       dtset%icoulomb = 1
    else if (do_realspace==0) then 
       write(2244,'(a)')'(chen_exx) NOTICE: recip space : xc_energy and xc_pot '
    endif

    ! E^HF_x, real space integral, eq.1 ref[2]
    !------------------------------------------
    ex = 0.d0
    do m = 1,norb
       do n = 1,norb 
          tmpr = phir(:,m)*phir(:,n)
          ! PSolver will do recip or realspace 
          ! depending on dtset%icouloum automatically
!          call PSolver_hartree(dtset,tmp_ehart,mpi_enreg,tmpr,rprimd,tmpr_out)

          call PSolver_hartree(tmp_ehart,hg,dtset%icoulomb,&
          me, dtset%nfft, dtset%ngfft(1:3), nproc, dtset%nscforder,&
          dtset%nspden, tmpr, tmpr_out, dtset%usewvl)
 
          tmpr = tmpr_out
          ex = ex - occ(n)*occ(m)/4.d0* & 
               sum(phir(:,m)*phir(:,n)*tmpr)*cellvol/dble(nfft)
       enddo  ! orbital
    enddo  ! orbital

    if (do_realspace==2) then 
       dtset%icoulomb = old_icoulomb
    endif

    write(2244,'(a,es16.8)')'(chen_exx): final_ex=', ex

    ! dEx/d phir, eq.11 ref[2]
    !-------------------------
    dex_dphir = 0.d0
    do n = 1,norb
       do m = 1,norb
          tmpr = phir(:,n)*phir(:,m)
          ! PSolver will do recip or realspace depending
          ! on dtset%icouloum automatically

          call PSolver_hartree(tmp_ehart,hg,dtset%icoulomb,&
          me, dtset%nfft, dtset%ngfft(1:3), nproc, dtset%nscforder,&
          dtset%nspden, tmpr, tmpr_out, dtset%usewvl)
!          call PSolver_hartree(dtset,tmp_ehart,mpi_enreg,tmpr,rprimd,tmpr_out)

          tmpr = tmpr_out
          dex_dphir(:,n) = dex_dphir(:,n) - occ(n)*occ(m)*phir(:,m)*tmpr
       enddo
       write(2244,'(a,2es12.4,a,i3)') '(chen_exx): min/max(d(EXX)/d(phir))=', & 
            minval(dex_dphir(:,n)),maxval(dex_dphir(:,n)),' <- orbital:',n
    enddo

    dtset%icoulomb = old_icoulomb

    do iband=1,norb
       dex_dphir(:,iband) = dex_dphir(:,iband) & 
            / ( 2.d0 * phir(:,iband) * occ(iband) )
    enddo


    deallocate(tmpr)
    deallocate(tmpr_out)
    deallocate(phir)

    write(2244,'(a,es16.8)')'leave chen_exx()'

  end subroutine chen_exx

  subroutine flo_exx(mpi_enreg,nfft,norb1,norb2,gmet,gsqcut,rprimd,wtk, & 
       &cellvol,phir1,phii1,phir2,phii2,occ1,occ2,&
       &ngfftf,dtset,ex,dex_dphir,dex_dphir_im,do_realspace)

    use defs_basis
    use defs_datatypes
    use defs_abitypes
    use defs_scftypes
    use defs_parameters

    implicit none

    !!  type(dataset_type),intent(in) :: dtset
    type(dataset_type) :: dtset

    type(MPI_type),intent(in) :: mpi_enreg

    integer , intent(in)  ::  & 
         do_realspace, &  ! flag for real space or recip space
                                ! 0: everything in recip space
                                ! 1: everything in real space
                                ! 2: only xc energy in real space, 
                                !    xc_pot in recip space
         
         nfft, &          ! number of point in real space
         norb1, &          ! number of orbitals to be consider
         norb2, &          ! number of orbitals to be consider
         ngfftf(18)    

    real(dp),intent(in) ::  & 
         wtk, &                         ! weight factor for k-points
         rprimd(3,3), & 
         gmet(3,3),  &                  ! metric of g-space
         gsqcut,     &                  ! cutoff of g^2
         cellvol, &                     ! cell volumn
         occ1(norb1), &                   ! occupation numbers
         occ2(norb2), &                   ! occupation numbers
         phir1(nfft,norb1), &             ! orbital in real space on k-point 1
         phii1(nfft,norb1), &             ! orbital in real space on k-point 1
         phir2(nfft,norb2), &             ! orbital in real space on k-point 1
         phii2(nfft,norb2)                ! orbital in real space on k-point 2

    real(dp),intent(out) :: & 
         ex, &                          ! E^HF_x energy
         dex_dphir(nfft,norb1), &       ! output: dE^HF_x / d( phir_n )
         dex_dphir_im(nfft,norb1)       ! output: dE^HF_x / d( phir_n )

    ! local vars 
    !===============
    real(dp), allocatable  :: tmpr(:)      ! temp vars in real space
    real(dp), allocatable  :: tmpr_out(:)  ! temp vars for output of PSolver
    real(dp), allocatable  :: tmpi(:)      ! temp vars in real space
    real(dp), allocatable  :: tmpi_out(:)  ! temp vars for output of PSolver

    real(dp) :: tmp_ehart,   & 
         tmp_qphon(3)=0.d0, &   ! no use
         mu = 0.d0              ! no use

    integer  :: usepaw = 0, & 
         n1,n2,old_icoulomb

    integer nproc, me, comm, i

    real*8 hg(3)


    !! function begins 
    !!==================

    !!
    !! only for closed shell system
    !! cannot be used for H atom or other open-shell systems
    !! in thoes cases, the EXX should be write carefully
    !!

    if(mpi_enreg%paral_kgb==1) then
       comm=mpi_enreg%comm_fft
       me=mpi_enreg%me_fft
       nproc=mpi_enreg%nproc_fft
    else
       comm=mpi_enreg%comm_cell
       nproc=xcomm_size(comm)
       me=xcomm_rank(comm)
    end if

    do i=1,3
       hg(i) = dsqrt(&
            &  rprimd(1,i)**2 + rprimd(2,i)**2 + rprimd(3,i)**2 &
            & )/dtset%ngfft(i)
    enddo

    write(2244,'(a)')'get in flo_exx()...'

    allocate(tmpr(nfft))      ! temp vars in real space
    allocate(tmpr_out(nfft))  ! temp vars for output of PSolver
    allocate(tmpi(nfft))      ! temp vars in real space
    allocate(tmpi_out(nfft))  ! temp vars for output of PSolver

    old_icoulomb   = dtset%icoulomb  ! back up icoulomb

    if (do_realspace==1) then 
       write(2244,'(a)')'(flo_exx) NOTICE: real space : xc_energy and xc_pot'
       dtset%icoulomb = 1
    else if (do_realspace==2) then 
       write(2244,'(a)')'(flo_exx) NOTICE: real space : xc_energy, recip space: xc_pot'
       dtset%icoulomb = 1
    else if (do_realspace==0) then 
       write(2244,'(a)')'(flo_exx) NOTICE: recip space : xc_energy and xc_pot '
    endif

    ! E^HF_x, real space integral, eq.1 ref[2]
    !------------------------------------------
    ex = 0.d0
    dex_dphir = 0.d0
    dex_dphir_im = 0.d0

    do n1 = 1,norb1
       do n2 = 1,norb2
          tmpr =   phir1(:,n1)*phir2(:,n2) + phii1(:,n1)*phii2(:,n2)
          tmpi =   phii1(:,n1)*phir2(:,n2) - phir1(:,n1)*phii2(:,n2)
          ! PSolver will do recip or realspace 
          ! depending on dtset%icouloum automatically

          call PSolver_hartree(tmp_ehart,hg,dtset%icoulomb,&
          me, dtset%nfft, dtset%ngfft(1:3), nproc, dtset%nscforder,&
          dtset%nspden, tmpr, tmpr_out, dtset%usewvl)

          call PSolver_hartree(tmp_ehart,hg,dtset%icoulomb,&
          me, dtset%nfft, dtset%ngfft(1:3), nproc, dtset%nscforder,&
          dtset%nspden, tmpi, tmpi_out, dtset%usewvl)

          ex = ex - occ1(n1)*occ2(n2)/4.d0* wtk *& 
               (sum(tmpr(:)*tmpr_out)+&
               &sum(tmpi(:)*tmpi_out))*cellvol/dble(nfft)

          dex_dphir(:,n1) = dex_dphir(:,n1) &
               & - occ1(n1)*occ2(n2)*phir2(:,n2)*tmpr_out &
               & + occ1(n1)*occ2(n2)*phii2(:,n2)*tmpi_out
          dex_dphir_im(:,n1) = dex_dphir_im(:,n1) &
               & - occ1(n1)*occ2(n2)*phir2(:,n2)*tmpi_out &
               & - occ1(n1)*occ2(n2)*phii2(:,n2)*tmpr_out

       enddo  ! orbital
    enddo  ! orbital

    if (do_realspace==2) then 
       dtset%icoulomb = old_icoulomb
    endif

    write(2244,'(a,es16.8)')'(flo_exx): final_ex=', ex
!
!    ! dEx/d phir, eq.11 ref[2]
!    !-------------------------
!    do n1 = 1,norb1
!       do n2 = 1,norb2
!          tmpr =   phir1(:,n1)*phir2(:,n2) + phii1(:,n1)*phii2(:,n2)
!          tmpi =   phii1(:,n1)*phir2(:,n2) - phir1(:,n1)*phii2(:,n2)
!          ! PSolver will do recip or realspace depending
!          ! on dtset%icouloum automatically
!
!          call PSolver_hartree(tmp_ehart,hg,dtset%icoulomb,&
!          me, dtset%nfft, dtset%ngfft(1:3), nproc, dtset%nscforder,&
!          dtset%nspden, tmpr, tmpr_out, dtset%usewvl)
!          call PSolver_hartree(tmp_ehart,hg,dtset%icoulomb,&
!          me, dtset%nfft, dtset%ngfft(1:3), nproc, dtset%nscforder,&
!          dtset%nspden, tmpi, tmpi_out, dtset%usewvl)
!
!          dex_dphir(:,n1) = dex_dphir(:,n1) &
!               & - occ1(n1)*occ2(n2)*phir2(:,n2)*tmpr_out &
!               & + occ1(n1)*occ2(n2)*phii2(:,n2)*tmpi_out
!          dex_dphir_im(:,n1) = dex_dphir_im(:,n1) &
!               & - occ1(n1)*occ2(n2)*phir2(:,n2)*tmpi_out &
!               & - occ1(n1)*occ2(n2)*phii2(:,n2)*tmpr_out
!
!       enddo
!       write(2244,'(a,2es12.4,a,i3)') '(flo_exx): min/max(d(EXX)/d(phir))=', & 
!            minval(dex_dphir(:,n1)),maxval(dex_dphir(:,n1)),' <- orbital:',n1
!    enddo

    dtset%icoulomb = old_icoulomb

    deallocate(tmpr)
    deallocate(tmpr_out)
    deallocate(tmpi)
    deallocate(tmpi_out)

    write(2244,'(a,es16.8)')'leave flo_exx()'

  end subroutine flo_exx

  subroutine flosum_exx(cg, mcg, mpi_enreg, gmet, gsqcut, rprimd, &
       ucvol, exc, nfft, dtset, npwarr, kg, occ, run_params)

    implicit none

    type(MPI_type),intent(inout) :: mpi_enreg
    real*8, intent(in) :: gmet(3,3)
    real*8, intent(in) :: gsqcut
    real(dp), intent(in) :: rprimd(3,3)
    real(dp), intent(in) :: ucvol
    real(dp), intent(inout) :: exc

    integer mcg
    real(dp), intent(in) :: cg(2,mcg)

    integer, intent(in) :: nfft

    type (run_parameter_type) run_params

    type(dataset_type) :: dtset
   
    real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

    integer, intent(in) :: npwarr(dtset%nkpt)

    integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)

    real*8, allocatable :: phir(:,:,:)
    real*8, allocatable :: phii(:,:,:)

    real*8 wtk

    integer iband, icg_shift, npw, ikpt, ikg, norb

    integer ikpt1, ikpt2
    integer ikg1, ikg2, nb1, nb2

    icg_shift=0

    allocate(phir(nfft, dtset%mband, dtset%nkpt))
    allocate(phii(nfft, dtset%mband, dtset%nkpt))

    ikg = 0

    icg_shift = 0

    do ikpt=1, dtset%nkpt

       npw = npwarr(ikpt)

       do iband=1,dtset%nband(ikpt)
             
          call recip_to_cmplx(mpi_enreg,dtset%paral_kgb,dtset%nkpt,&
               npw,ikpt, dtset%istwfk(ikpt),dtset%nspinor,&
               dtset%mgfft,npw,ucvol,&
               dtset%ngfft,kg(:,ikg+1:ikg+npw), npwarr(ikpt),&
               cg(:,icg_shift+1:icg_shift+&
               npw*dtset%nspinor),phir(:,iband,ikpt),phii(:,iband,ikpt))

          icg_shift=icg_shift + npw*dtset%nspinor
       enddo
          
       ikg = ikg + npw
             
    enddo

    ikg1 = 0

    do ikpt1=1, dtset%nkpt
       ikg2 = 0
       do ikpt2=1, dtset%nkpt

          nb1 = dtset%nband(ikpt1)
          nb2 = dtset%nband(ikpt2)
          wtk = dtset%wtk(ikpt1)*dtset%wtk(ikpt2)

          call flo_exx(mpi_enreg,nfft,nb1,nb2,gmet,gsqcut,rprimd, & 
               wtk,ucvol,phir(:,:,ikpt1),phii(:,:,ikpt1),&
               phir(:,:,ikpt2),phii(:,:,ikpt2),&
               occ(ikg1+1),occ(ikg2+1),&
               dtset%ngfft,dtset,exc,run_params%uxc_data(:,:,:,ikpt1),&
               run_params%uxc_im_data(:,:,:,ikpt1),run_params%do_realspace)

          ikg2 = ikg2 + dtset%nband(ikpt2)
       enddo
       do iband=1,nb1
          run_params%uxc_data(:,1,iband,ikpt1) = &
               run_params%uxc_data(:,1,iband,ikpt1) & 
               / ( 2.d0 * occ(iband+ikg1) )
          
          run_params%uxc_im_data(:,1,iband,ikpt1) = &
               run_params%uxc_im_data(:,1,iband,ikpt1) & 
               / ( 2.d0 * occ(iband+ikg1) )
          
          !                  / ( 2.d0 * phir(:,iband,ikpt1) * &
       enddo
       ikg1 = ikg1 + dtset%nband(ikpt1)
    enddo
    
    deallocate(phir)
    deallocate(phii)
  end subroutine flosum_exx

  subroutine allocate_uxc(nfft, nkpt, nband, nsppol, nspden, run_params)

    implicit none

    integer, intent(in) :: nfft
    integer, intent(in) :: nkpt
    integer, intent(in) :: nband(nkpt)
    integer, intent(in) :: nsppol
    integer, intent(in) :: nspden

    type (run_parameter_type) run_params

    integer ikpt, n_tot_band

    allocate(run_params%nband(nkpt*nsppol))
    allocate(run_params%wtk(nkpt))
 
    if (run_params%xc_type == 1) then
       ! do not allocate anything, lda_pot is already there
       ! and exc already has the correct value
    elseif (run_params%xc_type == 2) then ! exact exchange, nonperiodic
       n_tot_band = 0
       do ikpt=1,nkpt
          n_tot_band = n_tot_band + nband(ikpt)
       enddo
       allocate(run_params%uxc_data(nfft, nspden, n_tot_band,nkpt))
    elseif (run_params%xc_type == 3) then ! exact exchange, periodic
       n_tot_band = 0
       do ikpt=1,nkpt
          n_tot_band = n_tot_band + nband(ikpt)
       enddo
       allocate(run_params%uxc_data(nfft, nspden, n_tot_band,nkpt))
       allocate(run_params%uxc_im_data(nfft, nspden, n_tot_band,nkpt))
    elseif (run_params%xc_type == 4) then ! exact exchange, Spencer et al.
       n_tot_band = 0
       do ikpt=1,nkpt
          n_tot_band = n_tot_band + nband(ikpt)
       enddo
       allocate(run_params%uxc_data(nfft, nspden, n_tot_band,nkpt))
       allocate(run_params%uxc_im_data(nfft, nspden, n_tot_band,nkpt))
    else
       write(0,*) "unknown XC specifier, run_params%xc_type"
       write(0,*) "see 67_common/get_XC_potential.F90 for a list"
       stop
    endif

  end subroutine allocate_uxc

  subroutine prepare_uxc(cg, mcg, mpi_enreg, gmet, gsqcut, rprimd, &
       ucvol, exc, nfft, dtset, npwarr, kg, occ, run_params)

    implicit none

    type(MPI_type),intent(inout) :: mpi_enreg
    real*8, intent(in) :: gmet(3,3)
    real*8, intent(in) :: gsqcut
    real(dp), intent(in) :: rprimd(3,3)
    real(dp), intent(in) :: ucvol
    real(dp), intent(inout) :: exc

    integer mcg
    real(dp), intent(in) :: cg(2,mcg)

    integer, intent(in) :: nfft

    type (run_parameter_type) run_params

    type(dataset_type) :: dtset
   
    real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

    integer, intent(in) :: npwarr(dtset%nkpt)

    integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)

    real*8 wtk

    integer iband, icg_shift, npw, ikpt, ikg, norb

    integer ikpt1, ikpt2
    integer ikg1, ikg2, nb1, nb2


    run_params%nband = dtset%nband
    run_params%wtk = dtset%wtk

    if (run_params%xc_type == 1) then
       ! do not do anything, lda_pot is already there
       ! and exc already has the correct value
    elseif (run_params%xc_type == 2) then ! exact exchange, nonperiodic

       npw = npwarr(1)

       call chen_exx(mpi_enreg,npw,nfft,dtset%nband(1),&
            gmet,gsqcut,rprimd, & 
            ucvol,occ,dtset%ngfft,dtset,cg,mcg,kg,&
            exc,run_params%uxc_data,run_params%do_realspace)

    elseif (run_params%xc_type == 3) then ! exact exchange, periodic

       call flosum_exx(cg, mcg, mpi_enreg, gmet, gsqcut, rprimd, &
       ucvol, exc, nfft, dtset, npwarr, kg, occ, run_params)       

    elseif (run_params%xc_type == 4) then ! exact exchange, 
       !periodic, Spencer et al.

       call Spencer_exx(cg, dtset, gmet, exc, run_params%uxc_data,&
            & kg, mcg, mpi_enreg, nfft, npwarr, occ, ucvol)

    else
       write(0,*) "unknown XC specifier, run_params%xc_type"
       write(0,*) "see 67_common/get_XC_potential.F90 for a list"
       stop
    endif

  end subroutine prepare_uxc

  subroutine get_uxc(wf_real,rr_b, wf_imag,ri_b,&
       & nfft, ispden, iband, ikpt, run_params)

    implicit none

    real*8, intent(in) :: wf_real(nfft)
    real*8, intent(in) :: wf_imag(nfft)
    real*8, intent(out) :: rr_b(nfft)
    real*8, intent(out) :: ri_b(nfft)

    integer, intent(in) :: nfft
    integer, intent(in) :: ispden

    integer, intent(in) :: iband ! current band
    integer, intent(in) :: ikpt  ! current k-point

    type (run_parameter_type), intent(in) :: run_params

    if (run_params%xc_type == 1) then
       ri_b(:) = run_params%lda_pot(:,ispden) * wf_imag(:)
       rr_b(:) = run_params%lda_pot(:,ispden) * wf_real(:) ! LDA potential does not depend
       ! on kpoint or band
    elseif (run_params%xc_type == 2) then ! exact exchange, non-periodic
       ri_b(:) = run_params%uxc_data(:,ispden,iband,ikpt) * wf_imag(:)
       rr_b(:) = run_params%uxc_data(:,ispden,iband,ikpt) * wf_real(:)
    elseif (run_params%xc_type == 3) then ! exact exchange, periodic
       ri_b(:) = run_params%uxc_im_data(:,ispden,iband,ikpt)
       rr_b(:) = run_params%uxc_data(:,ispden,iband,ikpt)
    elseif (run_params%xc_type == 4) then ! exact exchange, periodic
       ri_b(:) = run_params%uxc_im_data(:,ispden,iband,ikpt)
       rr_b(:) = run_params%uxc_data(:,ispden,iband,ikpt)
    else
       write(0,*) "unknown XC specifier, run_params%xc_type"
       write(0,*) "see 67_common/get_XC_potential.F90 for a list"
       stop
    endif

  end subroutine get_uxc

  subroutine free_uxc( run_params )

    implicit none

    type (run_parameter_type) run_params

    deallocate(run_params%nband)

    if (run_params%xc_type == 1) then
       ! did not allocate anything in the first place
    elseif (run_params%xc_type == 2) then ! exact exchange, non-periodic
       deallocate(run_params%uxc_data)
    elseif (run_params%xc_type == 3) then ! exact exchange, periodic
       deallocate(run_params%uxc_data)
       deallocate(run_params%uxc_im_data)
    elseif (run_params%xc_type == 4) then ! exact exchange, periodic
       deallocate(run_params%uxc_data)
       deallocate(run_params%uxc_im_data)
    else
       write(0,*) "unknown XC specifier, run_params%xc_type"
       write(0,*) "see 67_common/get_XC_potential.F90 for a list"
       stop
    endif

  end subroutine free_uxc

end module get_xc_potential
