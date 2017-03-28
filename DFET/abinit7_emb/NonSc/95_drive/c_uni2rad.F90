!--------------------------------------------------------------
! project a 3D arrasy, defined on uniform grid in space, onto
! spherical harmonics cetenerd at each atom.
!
!  Chen Huang
!--------------------------------------------------------------
subroutine c_uni2rad(mpi_enreg,paral_kgb,nfftf,ngfftf, & 
   natom,typat,psps,xcart,gprimd,pawrad,extpot_lm,ucvol,gsqcut)

  use c_hc_vars

  use defs_basis
  use defs_datatypes
  use defs_abitypes

  use interfaces_12_hide_mpi
  use m_io_tools, only : flush_unit
 
  implicit none
  
  type(MPI_type),intent(in) :: mpi_enreg;
  integer, intent(in):: & 
    paral_kgb, &
    ngfftf(18), & 
    nfftf, &                      ! number of grid points in real space
    natom                         ! number of atoms

  real(dp),intent(in) :: & 
    gsqcut, &
    ucvol, &                        ! volume of cell
    gprimd(3,3), &                  ! 
    xcart(3,natom)                  ! coordiantes of atoms

  integer,intent(in) :: typat(natom)       ! array for atom types
  type(pseudopotential_type),intent(in) :: psps
  type(pawrad_type),intent(in):: pawrad(psps%ntypat)
  real(dp),intent(out):: extpot_lm(hc_meshsz,hc_lsize**2,natom) 

  !>>>>>>>> local vars <<<<<<<<<!

  character(len=500) :: message
  integer :: ia,ilm,ir,totgpt
  integer :: iq,qx,qy,qz,mm,ll,idx,klm
  integer :: normchoice,option
  real(dp):: c_norm                       ! function for norm calculation
  complex(dp) :: cplx,imag
  real(dp)::  & 
    potg(2,ngfftf(1)*ngfftf(2)*ngfftf(3)), &  
    phase, &
    extpot_lm_im(hc_meshsz,hc_lsize**2,natom), &  ! store imaginary part of extpot_lm for debug purpose only
    gcart(ngfftf(1),ngfftf(2),ngfftf(3),3), &
    qvec(3,ngfftf(1)*ngfftf(2)*ngfftf(3)), &
    qnorm(ngfftf(1)*ngfftf(2)*ngfftf(3))
  real(dp) :: & 
    besp,bespp, &
    cos_phase,sin_phase, &
    rmesh_int, &
    coord(3), & 
    tmp1,tmp2
  integer::& 
    ii,atype,mesh_size, &
    rmesh_pts
  real(dp),allocatable:: & 
    slmrq(:), &
    bes(:,:), & 
    bes_slmrq(:,:), &
    qr(:),& 
    re(:,:),  & 
    im(:,:),  & 
    ylmr(:,:),  &
    fit_result(:), & 
    vtmp(:,:), &
    vtmp_im(:,:), &
    alpha(:), &
    rmesh(:), & 
    ire(:,:),iim(:,:)

! bessel interpolation
  integer  :: nqr,lmsize
  integer,allocatable :: index1(:)
  real(dp) :: qr_max,qr_spa,xx
  real(dp),allocatable :: y_bes(:,:)

! parallel over atoms
  logical :: do_atom
  integer :: & 
    iatom , &  ! loop index for atoms
    my_atom, & 
    atom_ptr_parallel, & 
    natom_proc, &   ! how many atoms per proc
    ierr

!!DEBUG
  character(len=100) :: chtmp

!********************************************************
!        function begins 
  write(message,'(a)') '(chen/c_uni2rad) enter '
  call wrtout(std_out,message,'COLL')

! Initialize local variable.
! extpot_lm is calculated on a radial mesh (even spacing) first
! and then interpolted onto log-mesh later
!
  rmesh_int=min(0.1d0,(ucvol/dble(nfftf))**0.33333_dp) ! bohr
  rmesh_pts=ceiling(hc_rmax/rmesh_int)+3    ! a uniform grid
  allocate(rmesh(rmesh_pts))
  do ii=1,rmesh_pts
    rmesh(ii)=rmesh_int*real(ii-1,kind=dp)
  enddo

  lmsize=hc_lsize**2    ! hc_lsize: maximum value of angular momentum + 1 leading to nonzero Gaunt coeff.
  extpot_lm_im = 0.d0
  extpot_lm    = 0.d0

  allocate(slmrq(lmsize))
  allocate(vtmp(rmesh_pts,lmsize),vtmp_im(rmesh_pts,lmsize),qr(rmesh_pts))
  allocate(re(rmesh_pts,lmsize),im(rmesh_pts,lmsize))
  allocate(bes_slmrq(rmesh_pts,lmsize))
  allocate(bes(rmesh_pts,lmsize),index1(rmesh_pts),alpha(rmesh_pts))

! information from gridgcart.F90 --------
!
! In FFT algorithms, the transformed data is ordered such that zero comes first, then
! increasing positive frequencies, then *decreasing* negative frequencies. The cases of
! even and odd data must be distinguished, as follows.
! If ngfft(i) is even, then components ngfft(i)/2 and -ngfft(i)/2 are identical and
! both held in data point ngfft(i)/2 + 1 in the output vector. If ngfft(i) is odd, they
! are offset. Example:
!
  call fourdp(1,potg,extpot,-1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)

  write(message,'(a,F12.6,a)')"   hc_rmax   =",hc_rmax," bohr";  call wrtout(std_out,message,"COLL")
  write(message,'(a,F12.6,a)')"   rmesh_int =",rmesh_int," bohr";  call wrtout(std_out,message,"COLL")
  write(message,'(a,F12.6,a)')"   rmesh_max =",rmesh(rmesh_pts)," bohr";  call wrtout(std_out,message,"COLL")
  write(message,'(a,I5)')     "   rmesh_pts =",rmesh_pts;call wrtout(std_out,message,"COLL")
 
  totgpt=ngfftf(1)*ngfftf(2)*ngfftf(3)
  call gridgcart(gcart,gprimd,ngfftf)
  
  do qz=1,ngfftf(3)
    do qy=1,ngfftf(2)
      do qx=1,ngfftf(1)
        idx=qx+(qy-1)*ngfftf(1)+(qz-1)*ngfftf(1)*ngfftf(2)
!       ------------------------------------------------------
!       NOTICE: abinit's G-vector is used as
!       (2*pi)g*r in the exponent in Fourier Transform
!       so here we add two_pi to gcart
!
        qvec(:,idx)=two_pi*gcart(qx,qy,qz,:)
        qnorm(idx)=c_norm(qvec(:,idx))
      enddo
    enddo
  enddo

! ---------------------------------------------------
! calculate real spherical harmonics => Y_lmr
! initylmr.F90 uses appendix equation (A.5)
! in M. Torrent, Computational Material Science, (2007)
! to compuate slmr, in fact stored in ylmr
!
  normchoice=1  ! the input vector is NOT normalized
  option    =1  ! compute only Ylm_r
  allocate(ylmr(lmsize,ngfftf(1)*ngfftf(2)*ngfftf(3)))
  call initylmr(hc_lsize,normchoice,totgpt,qnorm,option,qvec,ylmr)

! ==========================================
! prepare spline for bessel(qr) to make 
! the evaluation of bessel function faster later
!
  qr_spa = 0.01d0
  qr_max = maxval(qnorm)*maxval(rmesh)
  nqr    = ceiling(qr_max/qr_spa)+10
  allocate(y_bes(nqr,lmsize))
  do ll=0,hc_lsize-1
    do mm=-ll,ll
      do ii=1,nqr
        xx=qr_spa*dble(ii-1)
        call jbessel(y_bes(ii,ll**2+ll+1+mm),besp,bespp,ll,0,xx)
      enddo
    enddo
  enddo
  write(message,'(a,F12.6,a)')'   qr_max    = ',qr_max," "; call wrtout(std_out,message,'COLL')
  write(message,'(a,I10)')    '   nqr       = ',nqr;    call wrtout(std_out,message,'COLL')

! ========================================
! prepare imaginary_unit^L array
!
  allocate(ire(rmesh_pts,lmsize),iim(rmesh_pts,lmsize))
  do ll=0,hc_lsize-1
    do mm=-ll,ll
      call imag_l(ll,tmp1,tmp2)
      ire(:,ll**2+ll+1+mm)=tmp1
      iim(:,ll**2+ll+1+mm)=tmp2
    enddo
  enddo

  if (mpi_enreg%nproc>1) then 
    write(message,'(a)')'(chen/c_uni2rad) start to loop over atoms ---------'
    call wrtout(std_out,message,"COLL")
    atom_ptr_parallel = 1
    call c_print_time(std_out)
    call flush_unit(std_out)
    call xbarrier_mpi(mpi_enreg%spaceComm,ierr)
  endif

! =============================
! loop over atoms
!
  do ia=1,natom
    
    !------------------------------------
    ! parallel: dispatch atoms to nodes
    if (mpi_enreg%nproc>1) then
      if (atom_ptr_parallel>natom) then 
        call c_print_time(std_out)
        write(message,'(a)')'(chen/c_uni2rad) finished loop over atoms ---------'
        call wrtout(std_out,message,"COLL")
        call flush_unit(std_out)
        exit
      endif
      my_atom = atom_ptr_parallel+mpi_enreg%me
      if (my_atom<=natom) then
        do_atom = .true.
        iatom = my_atom
        write(*,'(a,I3,a,I3,a,I3)')    & 
          " group:",atom_ptr_parallel, & 
          " atom_index ",iatom,' ==> node ',mpi_enreg%me
      else
        do_atom = .false.
      endif
      atom_ptr_parallel=atom_ptr_parallel+mpi_enreg%nproc
      call xbarrier_mpi(mpi_enreg%spaceComm,ierr)
      call flush_unit(std_out)
      call xbarrier_mpi(mpi_enreg%spaceComm,ierr)
      if (do_atom==.false.) exit
    else
      iatom = ia
    endif
    
    !-----------------------------------------------
    !------------ main code begin ------------------
    !-----------------------------------------------

    coord=xcart(:,iatom)
    
    if(mpi_enreg%nproc==1) then 
      write(message,'(a,I4)')'   doing extpot_lm for atom:',iatom
      call wrtout(std_out,message,'COLL')
      call c_print_time(std_out)
    endif
    vtmp    = 0.d0
    vtmp_im = 0.d0
    
    !=======================
    !loop over q points
    do iq=1,ngfftf(1)*ngfftf(2)*ngfftf(3)
 
      !skip q points outside g-sphere
      if (qnorm(iq)<sqrt(gsqcut*two_pi**2)) then
        
        slmrq(:)=ylmr(:,iq)                   ! S_lm(q)      
        qr(:)=qnorm(iq)*rmesh(:)              ! |q|.|r-R|
        index1(:)=floor(qr(:)/qr_spa)+1       ! index for searching in y_bes array
        alpha(:) =qr(:)/qr_spa-dble(index1(:)-1)  
        do ir=1,rmesh_pts
          bes(ir,:)=y_bes(index1(ir),:)*(1.d0-alpha(ir))+y_bes(index1(ir)+1,:)*alpha(ir) ! linear approximation
        enddo
        do klm=1,lmsize
          bes_slmrq(:,klm)=bes(:,klm)*slmrq(klm)
        enddo
        phase=dot_product(qvec(1:3,iq),coord(1:3))  ! q \dot R
        cos_phase=cos(phase)
        sin_phase=sin(phase)
        re=(cos_phase*potg(1,iq)-sin_phase*potg(2,iq))*bes_slmrq
        im=(cos_phase*potg(2,iq)+sin_phase*potg(1,iq))*bes_slmrq
        vtmp    = vtmp    + (re*ire-im*iim)
        vtmp_im = vtmp_im + (re*iim+ire*im)

      endif

    enddo 
    ! q-space
    !---------------------------

    vtmp    = four_pi*vtmp
    vtmp_im = four_pi*vtmp_im
    if (mpi_enreg%nproc==1) call c_print_time(std_out)
    
    !===========================
    ! interpolate our even spacing mesh to obtain extpot_lm array
    atype=typat(iatom)
    mesh_size=pawrad(atype)%mesh_size
    allocate(fit_result(mesh_size))
    do klm=1,lmsize
      call c_driver_spline(rmesh_pts,rmesh,vtmp(:,klm),mesh_size,pawrad(atype)%rad(1:mesh_size),fit_result)
       extpot_lm(1:mesh_size,klm,iatom)=fit_result(1:mesh_size)
      call c_driver_spline(rmesh_pts,rmesh,vtmp_im(:,klm),mesh_size,pawrad(atype)%rad(1:mesh_size),fit_result)
       extpot_lm_im(1:mesh_size,klm,iatom)=fit_result(1:mesh_size)
    enddo
    deallocate(fit_result)

  enddo !iatom
  
  !============================================
  ! parallel: gether all the results from nodes 
  if (mpi_enreg%nproc>1) then
    call xbarrier_mpi(mpi_enreg%spaceComm,ierr)
    call xsum_mpi(extpot_lm,mpi_enreg%spaceComm,ierr)
    call xbarrier_mpi(mpi_enreg%spaceComm,ierr)
    call xsum_mpi(extpot_lm_im,mpi_enreg%spaceComm,ierr)
    call xbarrier_mpi(mpi_enreg%spaceComm,ierr)
    if(mpi_enreg%me==0) then 
      print *,"*** (chen) extpot_lm and extpot_lm_im are reduced over nodes ***"
      call flush_unit(std_out)
    endif
  endif


!! DEBUG
!  if (.true. .and. mpi_enreg%me==0) then
!    do ii=1,natom
!      write(chtmp,*)ii
!      open(file='extpot_lm_l_0_atom_'//ADJUSTL(trim(chtmp)),unit=100,action='write')
!      open(file='extpot_lm_l_1_atom_'//ADJUSTL(trim(chtmp)),unit=101,action='write')
!      open(file='extpot_lm_l_2_atom_'//ADJUSTL(trim(chtmp)),unit=102,action='write')
!      do ir=1,hc_meshsz
!        write(100,*)extpot_lm(ir,1,1)  ! L=0,m=0
!        write(101,*)extpot_lm(ir,3,1)  ! L=1,m=0
!        write(102,*)extpot_lm(ir,7,1)  ! L=2,m=0
!      enddo ! rmesh
!      close(100)
!      close(101)
!      close(102)
!    enddo ! atoms
!  endif
! stop
!
!!ENDDEBUG  
  
  deallocate(rmesh,slmrq,qr,re,im)
  deallocate(bes_slmrq,bes,index1,ire,iim)
  deallocate(y_bes)
  deallocate(vtmp)
  deallocate(vtmp_im)
  deallocate(ylmr,alpha)

  write(message,'(2(a,ES12.4))')"   max_extpot    =",maxval(extpot_lm),   " min_extpot    =",minval(extpot_lm)
  call wrtout(std_out,message,'COLL')
  write(message,'(2(a,ES12.4))')"   max_extpot_im =",maxval(extpot_lm_im)," min_extpot_im =",minval(extpot_lm_im)
  call wrtout(std_out,message,'COLL')
  write(message,'(a)') '(chen): leave c_uni2rad()'
  call wrtout(std_out,message,'COLL')
  call xbarrier_mpi(mpi_enreg%spaceComm,ierr)

end subroutine c_uni2rad


!=============================================
! Calcualte Norm of a incoming vector
!=============================================
function c_norm(vec)
 use defs_basis
 use defs_datatypes

 implicit none

 real(dp),intent(in) :: vec(3)
 real(dp) :: c_norm
 c_norm = sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
 return
end function c_norm

!=============================================
! Calcualte (imaginary_unit)^L
!=============================================
subroutine imag_l(l,ire,iim)
  use defs_basis
  use defs_datatypes

  implicit none
  ! real and imag part of the L power 
  ! of the imaginary unit i

  integer :: l
  real(dp):: ire, iim

  if (mod(l,4)==0) then  !i^0
    ire=1.d0
    iim=0.d0
  endif
  if (mod(l,4)==1) then  !i^1
    ire=0.d0
    iim=1.0d0
  endif
  if (mod(l,4)==2) then !i^2
    ire=-1.d0
    iim=0.d0
  endif
  if (mod(l,4)==3) then !i^3
    ire=0.d0
    iim=-1.d0
  endif
      
end subroutine imag_l

