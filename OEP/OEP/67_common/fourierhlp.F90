module fourierhlp

  use defs_basis
  use defs_datatypes
  use defs_abitypes
  use m_xmpi
#if defined HAVE_MPI2
  use mpi
#endif

contains

  !!=====================================================
  !! complex convert recip_input from recip space to real space
  !!=====================================================
  subroutine recip_to_cmplx(mpi_enreg,paral_kgb,nkpt,npw,ikpt,istwf_k, & 
       nspinor,mgfft,nfft,ucvol,ngfft,&
       kg_k,npw_k,recip_input,res_real,res_imag)  
    use defs_datatypes
    use defs_abitypes  
    use interfaces_53_ffts  
    ! In FORTRAN, if optinal vars are used, we must use interface!!
    ! otherwise the optional vars will be screwed up
    ! here, we must include the interface for fourwf() subroutine
    ! this cost me four hours!

    implicit none 

    integer,intent(in)  :: npw,nspinor,nfft,nkpt,ngfft(18), & 
         istwf_k,mgfft,ikpt,kg_k(3,npw), & 
         npw_k,paral_kgb

    type(MPI_type):: mpi_enreg                 

    real(kind=8) :: recip_input(2,npw,nspinor)
    real(kind=8), intent(out)  :: res_real(nfft,nspinor)
    real(kind=8), intent(out)  :: res_imag(nfft,nspinor)

    real(kind=8)  :: ucvol
    real(kind=8)  :: fftg_out(2,npw)

    ! local vars 
    integer       :: gbound(2*mgfft+8,2)
    real(kind=8)  :: tmp_real(1,1,1), & 
         weight,dummy(2,1)

    real(kind=8), allocatable :: wfraug(:,:,:,:)

    integer i1,i2,i3,n1,n2,n3,nd1,nd2,nd3,me_fft,nproc_fft

    me_fft=ngfft(11) 
    nproc_fft=ngfft(10)

    n1  = ngfft(1)
    n2  = ngfft(2)
    n3  = ngfft(3)
    nd1 = ngfft(4)
    nd2 = ngfft(5)
    nd3 = ngfft(6)

    allocate(wfraug(2,ngfft(4),ngfft(5),ngfft(6)))

    ! get gbound
    call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)       
    weight = 1.d0/ucvol

    ! recip -> real space
    
    call fourwf(1,tmp_real,recip_input(:,:,1),dummy,wfraug,&
         gbound,gbound,istwf_k,&
         kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,1,ngfft(4),&
         ngfft(5),ngfft(6),0, paral_kgb,0,weight,weight)

    wfraug = wfraug / sqrt(ucvol)

    ! transfer wfraug(n4,n5,n6) array ==> tmp_phir(nfft,1) array

! this is from fftpac. Have not found how to call
! fftpac correctly so it does this for full complex
! realspace wf needed with k-point sampling
    do i3=1,n3
       if (((i3-1)/(n3/nproc_fft))==me_fft) then
          do i2=1,n2
             do i1=1,n1
                res_real(i1+n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1)),1)=&
                     wfraug(1,i1,i2,i3)
                res_imag(i1+n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1)),1)=&
                     wfraug(2,i1,i2,i3)
             end do
          end do
       endif
    end do

    if (nspinor==2) then
       call fourwf(1,tmp_real,recip_input(:,:,2),&
            dummy,wfraug,gbound,gbound,istwf_k,&
            kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,1,ngfft(4),&
            ngfft(5),ngfft(6),0, paral_kgb,0,weight,weight)

       do i3=1,n3
          if (((i3-1)/(n3/nproc_fft))==me_fft) then
             do i2=1,n2
                do i1=1,n1
                   res_real(i1+n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1)),2)=&
                        wfraug(1,i1,i2,i3)
                   res_imag(i1+n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1)),2)=&
                        wfraug(2,i1,i2,i3)
                end do
             end do
          endif
       end do
    endif


    deallocate(wfraug)

    return     
  end subroutine recip_to_cmplx

  !!=====================================================
  !! convert recip_input from recip space to real space
  !!=====================================================
  subroutine recip_to_real(mpi_enreg,paral_kgb,nkpt,npw,ikpt,istwf_k, & 
       nspinor,mgfft,nfft,ucvol,ngfft,&
       kg_k,npw_k,recip_input,res_real)  
    use defs_datatypes
    use defs_abitypes  
    use interfaces_53_ffts  
    ! In FORTRAN, if optinal vars are used, we must use interface!!
    ! otherwise the optional vars will be screwed up
    ! here, we must include the interface for fourwf() subroutine
    ! this cost me four hours!

    implicit none 

    integer,intent(in)  :: npw,nspinor,nfft,nkpt,ngfft(18), & 
         istwf_k,mgfft,ikpt,kg_k(3,npw), & 
         npw_k,paral_kgb

    type(MPI_type):: mpi_enreg                 

    real(kind=8) :: recip_input(2,npw*nspinor)
    real(kind=8), intent(out)  :: res_real(nfft)

    real(kind=8)  :: ucvol
    real(kind=8)  :: fftg_out(2,npw)

    ! local vars 
    integer       :: gbound(2*mgfft+8,2)
    real(kind=8)  :: tmp_real(1,1,1), & 
         weight,dummy(2,1)

    real(kind=8), allocatable :: wfraug(:,:,:,:)

    allocate(wfraug(2,ngfft(4),ngfft(5),ngfft(6)))

    ! get gbound
    call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)       
    weight = 1.d0/ucvol

    ! recip -> real space
    
    call fourwf(1,tmp_real,recip_input,dummy,wfraug,gbound,gbound,istwf_k,&
         kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,1,ngfft(4),&
         ngfft(5),ngfft(6),0, paral_kgb,0,weight,weight)

    wfraug = wfraug / sqrt(ucvol)

    ! transfer wfraug(n4,n5,n6) array ==> tmp_phir(nfft,1) array

    call fftpac(1,1,ngfft(1),ngfft(2),ngfft(3),&
         ngfft(4),ngfft(5),ngfft(6),ngfft,res_real,wfraug(1,:,:,:),1)

    deallocate(wfraug)

    return     
  end subroutine recip_to_real

  !!=====================================================
  !! convert input from real space to recip space
  !!=====================================================
  subroutine real_to_recip(mpi_enreg,paral_kgb,nkpt,npw,ikpt,istwf_k, & 
       nspinor,mgfft,nfft,ucvol,ngfft,&
       kg_k,npw_k,input_real,res_recip)  
    use defs_datatypes
    use defs_abitypes  
    use interfaces_53_ffts  ! In FORTRAN, if optinal vars are used, we must use interface!!
    ! otherwise the optional vars will be screwed up
    ! here, we must include the interface for fourwf() subroutine
    ! this cost me four hours!

    implicit none 

    integer,intent(in)  :: npw,nspinor,nfft,nkpt,ngfft(18), & 
         istwf_k,mgfft,ikpt,kg_k(3,npw), & 
         npw_k,paral_kgb

    type(MPI_type):: mpi_enreg                 

    real(kind=8) :: input_real(nfft)
    real(kind=8), intent(out)  :: res_recip(2,npw*nspinor)

    real(kind=8)  :: ucvol

    ! local vars 
    integer       :: gbound(2*mgfft+8,2)
    real(kind=8)  :: tmp_real(1,1,1),&
         weight,dummy(2,1)
    
    real(kind=8), allocatable  :: wfraug(:,:,:,:)

    allocate(wfraug(2,ngfft(4),ngfft(5),ngfft(6)))

    ! get gbound
    call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)       
    weight = 1.d0/ucvol

    ! real space -> recip space
    wfraug = 0.d0
    call fftpac(1,1,ngfft(1),ngfft(2),ngfft(3)&
         ,ngfft(4),ngfft(5),ngfft(6),ngfft,input_real,wfraug(1,:,:,:),2)
    wfraug = wfraug * sqrt(ucvol)
    call fourwf(1,tmp_real,res_recip,res_recip,wfraug,&
         gbound,gbound,istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,npw,&
         ngfft(4),ngfft(5),ngfft(6),3,paral_kgb,0,weight,weight)       

    deallocate(wfraug)

    return     
  end subroutine real_to_recip

  !!=====================================================
  !! convert input from real space to recip space
  !!=====================================================
  subroutine cmplx_to_recip(mpi_enreg,paral_kgb,nkpt,npw,ikpt,istwf_k, & 
       nspinor,mgfft,nfft,ucvol,ngfft,&
       kg_k,npw_k,input_real,input_imag,res_recip)  
    use defs_datatypes
    use defs_abitypes  
    use interfaces_53_ffts  ! In FORTRAN, if optinal vars are used, we must use interface!!
    ! otherwise the optional vars will be screwed up
    ! here, we must include the interface for fourwf() subroutine
    ! this cost me four hours!

    implicit none 

    integer,intent(in)  :: npw,nspinor,nfft,nkpt,ngfft(18), & 
         istwf_k,mgfft,ikpt,kg_k(3,npw), & 
         npw_k,paral_kgb

    type(MPI_type):: mpi_enreg                 

    real(kind=8) :: input_real(nfft,nspinor)
    real(kind=8) :: input_imag(nfft,nspinor)
    real(kind=8), intent(out)  :: res_recip(2,npw,nspinor)

    real(kind=8)  :: ucvol

    ! local vars 
    integer       :: gbound(2*mgfft+8,2)
    real(kind=8)  :: tmp_real(1,1,1),&
         weight,dummy(2,1)
    
    real(kind=8), allocatable  :: wfraug(:,:,:,:)

    integer i1,i2,i3,n1,n2,n3,nd1,nd2,nd3,me_fft,nproc_fft

    me_fft=ngfft(11) 
    nproc_fft=ngfft(10)

    n1  = ngfft(1)
    n2  = ngfft(2)
    n3  = ngfft(3)
    nd1 = ngfft(4)
    nd2 = ngfft(5)
    nd3 = ngfft(6)

    allocate(wfraug(2,ngfft(4),ngfft(5),ngfft(6)))

    gbound = 0
    ! get gbound
    call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)       
    weight = 1.d0/ucvol

    ! real space -> recip space
    wfraug = 0.d0

! this is from fftpac. Have not found how to call
! fftpac correctly so it does this for full complex
! realspace wf needed with k-point sampling
    do i3=1,n3
       if (((i3-1)/(n3/nproc_fft))==me_fft) then
          do i2=1,n2
             do i1=1,n1
                wfraug(1,i1,i2,i3)=&
                     input_real(i1+n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1)),1)
                wfraug(2,i1,i2,i3)=&
                     input_imag(i1+n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1)),1)
             end do
          end do
       endif
    end do

    wfraug = wfraug * sqrt(ucvol)
    call fourwf(1,tmp_real,res_recip,res_recip,wfraug,gbound,gbound,istwf_k,&
         kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,npw,&
         ngfft(4),ngfft(5),ngfft(6),3,paral_kgb,0,weight,weight)       

    if (nspinor == 2) then
       wfraug = 0.d0

       ! this is from fftpac. Have not found how to call
       ! fftpac correctly so it does this for full complex
       ! realspace wf needed with k-point sampling
       do i3=1,n3
          if (((i3-1)/(n3/nproc_fft))==me_fft) then
             do i2=1,n2
                do i1=1,n1
                   wfraug(1,i1,i2,i3)=&
                        input_real(i1+n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1)),2)
                   wfraug(2,i1,i2,i3)=&
                        input_imag(i1+n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1)),2)
                end do
             end do
          endif
       end do

       wfraug = wfraug * sqrt(ucvol)
       call fourwf(1,tmp_real,res_recip,res_recip,wfraug,gbound,gbound,istwf_k,&
            kg_k,kg_k,mgfft,mpi_enreg,1,ngfft,npw_k,npw,&
            ngfft(4),ngfft(5),ngfft(6),3,paral_kgb,0,weight,weight)       
    endif

    deallocate(wfraug)

    return     
  end subroutine cmplx_to_recip

end module fourierhlp
