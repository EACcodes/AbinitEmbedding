module spencerXX

  use defs_basis
  use defs_datatypes
  use defs_abitypes
  use m_xmpi
  use fourierhlp
  use m_geometry,      only : normv

  contains

  subroutine Spencer_exx(cg, dtset, gmet, ex, dex_dphir, kg, mcg, mpi_enreg,&
       & nfft, npwarr, occ, ucvol)

    use defs_basis
    use defs_datatypes
    use defs_abitypes
    use defs_scftypes
    use defs_parameters

    !! use a truncated Coulomb potential to effectively
    !! calculated EXX in reciprocal space for periodic systems,
    !! see J. Spencer and A. Alavi, PRB 77, 193110 (2008)

    implicit none

    type(dataset_type),intent(in) :: dtset

    real*8, intent(in) :: gmet(3,3)
    
    type(MPI_type),intent(in) :: mpi_enreg

    integer, intent(in) :: mcg
    real*8, intent(in) :: cg(2,mcg) 

    integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
    
    integer, intent(in) :: npwarr(dtset%nkpt)

    real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

    real*8, intent(in) :: ucvol

    real(dp),intent(out) :: & 
         ex, &                          ! E^HF_x energy
         dex_dphir(dtset%nfft,dtset%mband)   ! output: dE^HF_x / d( phir_n )

!    =================
    ! Local variables

!    Ypsilon(G) = sum_G'  cg(G') cg(G+G')
!               = F[ F [cg] . F [cg]]

    integer ikg, icg_shift, iband
    integer ikpt, ikpt1, ikpt2
    integer nspinor, nfft, norb1, norb2, npw

    real*8, allocatable:: wf_real(:,:,:,:), wf_imag(:,:,:,:)

    real*8 v_GKK

    real*8, allocatable:: Ymnkk(:,:)

    integer kpt1(3), kpt2(3)

    real*8, allocatable:: denominator(:)

    real*8 rcut

    real*8, allocatable :: tmpr(:)
    real*8, allocatable :: tmpi(:)

    integer ipw, m,n, maxnpw

    integer indx_k1, indx_k2

 ! === Dimensional primitive translations rprimd (from input), gprimd, metrics and unit cell volume ===

    nfft = dtset%nfft
    nspinor = dtset%nspinor

    allocate(tmpr(nfft))
    allocate(tmpi(nfft))

    allocate(wf_real(nfft, nspinor, dtset%mband, dtset%nkpt))
    allocate(wf_imag(nfft, nspinor, dtset%mband, dtset%nkpt))
    allocate(denominator(dtset%mpw))
    allocate(Ymnkk(2,dtset%mpw))

    ikg = 0
    
    icg_shift = 0

    do ikpt=1, dtset%nkpt
       
       npw = npwarr(ikpt)

       do iband=1,dtset%nband(ikpt)
          
          call recip_to_cmplx(mpi_enreg,dtset%paral_kgb,dtset%nkpt,npw, &
               ikpt,dtset%istwfk(ikpt),&              
               dtset%nspinor,dtset%mgfft,nfft,ucvol,&
               dtset%ngfft,kg(:,ikg+1:ikg+npw),npwarr(ikpt),&
               cg(:,icg_shift+1:icg_shift+npw*dtset%nspinor),&
               wf_real(:,1,iband,ikpt),wf_imag(:,1,iband,ikpt))

          icg_shift=icg_shift + npw*dtset%nspinor
       enddo
       
       ikg = ikg + npw
             
    enddo

    rcut = (ucvol*npwarr(1)*3.d0/four_pi)**third

    Ex = 0

    indx_k1 = 0
    do ikpt1=1,dtset%nkpt

       norb1 = dtset%nband(ikpt1)

       indx_k2 = 0

       do ikpt2=1,dtset%nkpt

          norb2 = dtset%nband(ikpt2)

          kpt1(:) = dtset%kptns(:,ikpt1)
          kpt2(:) = dtset%kptns(:,ikpt2)    

          do ipw=1,npwarr(1)

!             v_GKK=(kg(1,ipw) - kpt1(1) + kpt2(1))**2 + &
!                  &(kg(2,ipw) - kpt1(2) + kpt2(2))**2 + &
!                  &(kg(3,ipw) - kpt1(3) + kpt2(3))**2  
             


             v_GKK=normv(kpt2(:) - kpt1(:) + kg(:,ipw),gmet,'G')
           
             if (v_GKK.eq.0) then
                denominator(ipw) = two_pi*rcut**2
             else
                denominator(ipw) = four_pi / (dble(v_GKK))**2*&
                     &(1.d0 - dcos(dble(v_GKK)*rcut))
             endif
          enddo
    
          do m=1,norb1
             do n=1,norb2
                !Calculate co-density real part
                tmpr(:) = wf_real(:,1,m,ikpt1)*wf_real(:,1,n,ikpt2) + &
                     &wf_imag(:,1,m,ikpt1)*wf_imag(:,1,n,ikpt2) 
                !imaginary part
                tmpi(:) = wf_real(:,1,m,ikpt1)*wf_imag(:,1,n,ikpt2) - &
                     &wf_imag(:,1,m,ikpt1)*wf_real(:,1,n,ikpt2)
             
                call cmplx_to_recip(mpi_enreg,dtset%paral_kgb,&
                     & dtset%nkpt,npw,1,&
                     & dtset%istwfk(1), nspinor,dtset%mgfft,dtset%nfft,&
                     & ucvol,dtset%ngfft, kg,npwarr(1),tmpr,tmpi,Ymnkk)

                do ipw=1,npwarr(1)
                   Ex = Ex + (Ymnkk(1,ipw)**2 + Ymnkk(2,ipw)**2) &
                        &  *denominator(ipw ) *occ(indx_k1+m) *occ(indx_k2+n)
                enddo
          
             enddo
          enddo
          indx_k2 = indx_k2 + norb2
       enddo
       indx_k1 = indx_k1 + norb1
    enddo

    deallocate(wf_real)
    deallocate(wf_imag)
    deallocate(denominator)
    deallocate(Ymnkk)

    deallocate(tmpi)
    deallocate(tmpr)

  end subroutine Spencer_exx

end module spencerXX
