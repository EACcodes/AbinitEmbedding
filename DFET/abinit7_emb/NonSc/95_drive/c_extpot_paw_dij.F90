!===================================================================
! Compute the d_ij due to the external potential (extpot):
!
! Then these d_ij can then be combined with PAW's own d_ij during SCF
!
! Reference: pawdij.F90
!
! Chen Huang
!===================================================================
subroutine c_extpot_paw_dij(natom,ntypat,mpi_enreg,pawang,pawtab,paw_ij,paw_an,pawrad,typat, &
       extpot_lm,extpot_dij)

  use c_hc_vars

  use defs_basis
  use defs_datatypes
  use defs_abitypes

  implicit none
  integer,intent(in):: & 
    ntypat, &
    natom          

  integer,intent(in)  :: typat(natom)
  real(dp),intent(in) :: extpot_lm(hc_meshsz,hc_lsize**2,natom)
  type(MPI_type),intent(in)   :: mpi_enreg
  type(paw_an_type),intent(in):: paw_an(natom)
  type(paw_ij_type),intent(in):: paw_ij(natom)
  type(pawtab_type),intent(in):: pawtab(ntypat)
  type(pawang_type),intent(in):: pawang
  type(pawrad_type),intent(in):: pawrad(ntypat)
  real(dp),intent(out) :: extpot_dij(hc_lmn2size,natom)

! =====================
! local variables
  
  character(len=500) :: message
  integer :: ii,iatom,klm,klm1,klmn,kln,itypat,isel
  integer :: mesh_size,lmn2_size,ij_size
  real(dp):: tmp
  real(dp),allocatable :: vij1(:),vijtot(:),ff(:)
  integer,allocatable :: idum(:),indklmn(:,:)


!================================================
! based on pawdij.F90
! where vxc is expanded over (l,m) moments
! 

  write(message,'(a)')"(chen): enter c_extpot_paw_dij()"
  call wrtout(std_out,message,'COLL')


!  write(*,'(a,i2,a,i5,a,i5)')& 
!   '(my_dij with mpi_enreg) me=>',mpi_enreg%me,' nproc_atom=',mpi_enreg%nproc_atom, " natom=",mpi_enreg%natom
!  print *,'(my_lm) extpot_lm=',extpot_lm
  
  extpot_dij = zero

!!=========================================
! some checks, need to be removed later
  if (mpi_enreg%nproc_atom/=1) then
    print *, "src/95_drive/c_extpot_paw_dij.F90: nproc_atom/=1, stop"
    call leave_new('COLL')
  endif
!!=========================================  

! =======================================
! Loop over atoms
! =======================================
  do iatom=1,natom

    itypat   =typat(iatom)
    lmn2_size=paw_ij(iatom)%lmn2_size    ! (l1m1n1,l2m2n2) pairs
    ij_size  =pawtab(itypat)%ij_size     ! (l1n1,l2n2) pairs
    mesh_size=pawrad(itypat)%mesh_size

    allocate(vij1(ij_size))
    allocate(vijtot(lmn2_size))
    allocate(indklmn(6,lmn2_size))
    allocate(ff(1:mesh_size))

    indklmn(:,:)=pawtab(itypat)%indklmn(:,:)
    vijtot= zero
    vij1  = zero

!   =======================================================
!   ===== sum over (l,m) moments of extpot_LM
!   =======================================================
    do klm=1,hc_lsize**2

      vij1=zero
!     ==========================================
!     loop over (i,j) pair, because radial part!
      do kln=1,ij_size 
        ff(1:mesh_size)= extpot_lm(1:mesh_size,klm,iatom) * &
           (pawtab(itypat)%phiphj(1:mesh_size,kln)- & 
            pawtab(itypat)%tphitphj(1:mesh_size,kln))
        call simp_gen(vij1(kln),ff,pawrad(itypat))
      end do
      
!     ============================================
!     loop over the (n1l1,n2l2) pair.
!     just all the PAW basis pairs inside a sphere
!     remember: lmn2_size=lmn_size*(lmn_size+1)/2
!
      do klmn=1,lmn2_size
        klm1=indklmn(1,klmn)              ! select (l1m1,l2m2) for Phi_i and Phi_j
        isel=pawang%gntselect(klm,klm1)   ! select (LM,l1m1,l2m2) pair, gntselect>0 if Gaunt coeff is non-zero
        kln=indklmn(2,klmn)               ! select (i,j) pair for vij1(i,j)
        tmp=0.d0
        if (isel>0) tmp=vij1(kln)*pawang%realgnt(isel)
        vijtot(klmn)=vijtot(klmn)+tmp
      end do 
      ! Loop klmn (i,j) pair

    end do
!   =====================================
!   End of Loop klm (l,m) moments of extpot
!   =====================================

    extpot_dij(1:lmn2_size,iatom)=vijtot

    write(message,'(a,I3)')    ' --- for atom:',iatom;  call wrtout(std_out,message,'COLL')
    write(message,'(a,ES12.4)')'   max(extpot_dij) = ', maxval(extpot_dij(1:lmn2_size,iatom));    call wrtout(std_out,message,'COLL')
    write(message,'(a,ES12.4)')'   min(extpot_dij) = ', minval(extpot_dij(1:lmn2_size,iatom));    call wrtout(std_out,message,'COLL')

    deallocate(vij1)
    deallocate(vijtot)
    deallocate(indklmn)
    deallocate(ff)

  enddo ! Loop iatom

  write(message,'(a)')"(chen/c_extpot_paw_dij) leave"
  call wrtout(std_out,message,'COLL')

endsubroutine c_extpot_paw_dij
