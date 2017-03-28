!!-------------------------------------------------------
!!  Chen : this subroutine is based on denfgr()
!!         it uses smoothed PHI(r) to make the projected
!!         densities smoother, not the AE one
!! 
!!-------------------------------------------------------
!{\src2tex{textfont=tt}}
!!****f* ABINIT/denfgr
!! NAME
!! denfgr
!!
!! FUNCTION
!!  Construct complete electron density on fine grid, by removing nhat
!!  and adding PAW corrections
!!
!! COPYRIGHT
!!   Copyright (C) 2005-2010 ABINIT group (JWZ)
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~ABINIT/COPYING
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!   natom= number of atoms in cell
!!   nspden= number of spin densities
!!   ntypat= number of types of atoms in the cell
!!   nhat(pawfgr%nfft,nspden)= compensation charge density used in PAW
!!   pawfgr <type(pawfgr_type)>= data about the fine grid
!!   pawfgrtab(natom) <type(pawfgrtab_type)>= data about the fine grid around each atom
!!   pawrad(ntypat) <type(pawrad_type)>= radial mesh data for each type of atom
!!   pawrhoij(natom) <type(pawrhoij_type)>= rho_ij data for each atom
!!   pawtab(ntypat) <type(pawtab_type)>= PAW functions around each type of atom
!!   psps <type(pseudopotential_type)>= basic pseudopotential data
!!   rhor(pawfgr%nfft,nspden)= input density ($\tilde{n}+\hat{n}$ in PAW case)
!!   typat(natom)= list of atom types
!!
!! OUTPUT
!! rhor_paw(pawfgr%nfft,nspden)= full electron density on the fine grid
!!
!! NOTES
!!   In PAW calculations, the valence density present in rhor includes the
!!   compensation charge density $\hat{n}$, and also doesn't include the on-site
!!   PAW contributions. For post-processing and proper visualization it is necessary
!!   to use the full electronic density, which is what this subroutine constructs.
!!   Specifically, it removes $\hat{n}$ from rhor, and also computes the on-site PAW
!!   terms. This is nothing other than the proper PAW treatment of the density
!!   operator $|\mathbf{r}\rangle\langle\mathbf{r}|$, and yields the formula
!!   $$\tilde{n}+\sum_{ij}\rho_ij\left[\varphi_i(\mathbf{r})\varphi_j(\mathbf{r})-
!!   \tilde{\varphi}_i(\mathbf{r})\tilde{\varphi}_j(\mathbf{r})\right]$$
!!   Notice that this formula is expressed on the fine grid, and requires
!!   interpolating the PAW radial functions onto this grid, as well as calling
!!   initylmr in order to get the angular functions on the grid points.
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      initylmr,sort_dp,spline,splint,wrtout,xbarrier_mpi,xcomm_init
!!      xmaster_init,xme_init,xmpi_nproc,xsum_master
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine c_smooth_denfgr(MPI_enreg,natom,nspden,nhat,ntypat,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,rhor,rhor_paw,prtvol,psps,typat)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
 use interfaces_12_hide_mpi
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
 use interfaces_32_util
 use interfaces_51_manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nspden,ntypat,prtvol
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
!arrays
 type(MPI_type),intent(in) :: MPI_enreg
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: nhat(pawfgr%nfft,nspden)
 real(dp),intent(in) :: rhor(pawfgr%nfft,nspden)
 real(dp),intent(out) :: rhor_paw(pawfgr%nfft,nspden)
 type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 character(len=500) :: message
 integer :: delta,iatom,ierr,ifgd,ifftsph,inl,inrm,ipsang,irhoij
 integer :: ispden,itypat,il,im,ilm,iln,ilmn
 integer :: jl,jlm,jln,jm,j0lmn,jlmn
 integer :: klmn,master,my_start_indx,my_end_indx
 integer :: nfgd,nnl,normchoice,nprocs,option,rank,remainder,spaceComm
 real(dp) :: phj,phi,rR,tphj,tphi,ybcbeg,ybcend
!arrays
 integer,allocatable :: nrm_ifftsph(:)
 real(dp) :: ylmgr(3,3,0)
 real(dp),allocatable :: diag(:),nrm(:),phigrd(:,:),tphigrd(:,:),ylm(:,:),ypp(:)
 real(dp),allocatable :: rhor_paw_in_loop(:,:)

!!============== chen ========= 
 integer :: iost,ibasis,ii,jj,msz,max_mesh_size,max_basis_size,nrow
 character(len=200) :: ctmp,file_name
 real(dp) :: ftmp,aa,bb
 integer  :: ntmp
 real(dp),allocatable :: sm_phi(:,:,:)

! ************************************************************************

 DBG_ENTER("COLL")

 allocate(rhor_paw_in_loop(pawfgr%nfft,nspden))

!Initialise rhor_paw
 rhor_paw = zero

!Initialise and check parallell execution
 call xcomm_init  (MPI_enreg,spaceComm)
 call xmpi_nproc  (nprocs,ierr)
 call xmaster_init(MPI_enreg,master)
 call xme_init    (MPI_enreg,rank)

!=============== chen: read in smooth_phi files for each atom ======================
! Set up smooth (norm-conserving) partials
!
  write(message,*)"chen/c_smooth_denfgr: ================================================================"; call wrtout(std_out,message,'COLL')
  write(message,*)"chen/c_smooth_denfgr: =============== In src/66_paw/c_smooth_denfgr.F90 =============="; call wrtout(std_out,message,'COLL')
  write(message,*)"chen/c_smooth_denfgr: ================================================================"; call wrtout(std_out,message,'COLL')
  max_mesh_size =0
  max_basis_size=0
  do itypat=1,psps%ntypat
    if (max_mesh_size<pawrad(itypat)%mesh_size)& 
      max_mesh_size=pawrad(itypat)%mesh_size
    if (max_basis_size<pawtab(itypat)%basis_size)& 
      max_basis_size=pawtab(itypat)%basis_size
  enddo
  write(message,'(a,i)')"chen/c_smooth_denfgr: max_mesh_size  = ",max_mesh_size;  call wrtout(std_out,message,'COLL')
  write(message,'(a,i)')"chen/c_smooth_denfgr: max_basis_size = ",max_basis_size; call wrtout(std_out,message,'COLL')
  allocate(sm_phi(max_mesh_size,max_basis_size,psps%ntypat))

  sm_phi=0.d0
  write(message,'(a)')'chen/c_smooth_denfgr: reading smooth_phi ------------'
  call wrtout(std_out,message,'COLL')
  
  !--------------------------
  ! loop over atom types
  !
  do itypat=1,psps%ntypat
    
    ! prepare file names for this atom type
    !
    if (psps%znucltypat(itypat)<10) then
      write(file_name,'(I1,a)') int(psps%znucltypat(itypat)),'_smooth_phi'
    elseif(psps%znucltypat(itypat)>=10  .and. psps%znucltypat(itypat)<=99) then
      write(file_name,'(I2,a)') int(psps%znucltypat(itypat)),'_smooth_phi'
    elseif(psps%znucltypat(itypat)>=100) then
      write(file_name,'(I3,a)') int(psps%znucltypat(itypat)),'_smooth_phi'
    endif
    write(message,'(a,I3,a,a)') 'chen/c_smooth_denfgr: loading atom_type:',itypat,' => file: ',file_name
    call wrtout(std_out,message,'COLL')

    open(unit=111,file=adjustl(trim(file_name)),action='read',iostat=iost)
    if (iost/=0) then
      write(message,*)'src/66_paw/c_smooth_denfgr.F90: error in open file',adjustl(trim(file_name)),' stop'
      call wrtout(std_out,message,'COLL')
      call leave_new('COLL')
    endif

    read(111,*) ctmp  ! skip header
    msz=pawrad(itypat)%mesh_size
    
    !--------------------
    ! loop over basis
    do ibasis=1,pawtab(itypat)%basis_size
      read(111,*) ctmp,ntmp,ctmp,nrow
      write(message,'(a,I3,a,I5)')'  => reading basis:',ibasis,' nrow=',nrow
      call wrtout(std_out,message,'COLL')
      if (msz/=nrow) then 
        print *,'mesh_size/=nrow, error, in src/66_paw/c_smooth_denfgr.F90, must equal.'
        call leave_new('COLL')
      endif
      do jj=1,msz
        read(111,*)ftmp,ftmp,sm_phi(jj,ibasis,itypat) ! read in smooth-phi
      enddo
    enddo 
    ! end of loop over basises
    !----------------------------

    close(111)

  enddo 
  ! end of loop over atom_types
  !------------------------

!!DEBUG
!  allocate(load_data(350),load_rr(350))
!  open(file='out_vncpp_h068',unit=111,action='read')
!  read(111,*)(load_data(ii),ii=1,350)
!  close(111)
!  load_rr(1)=0.d0
!  aa=0.38106E-05
!  bb=0.35000E-01
!  do jj=2,350
!    load_rr(jj)=AA*exp(BB*dble(jj-2))
!  enddo
!  ! interpolate to sm_phi array
!  msz=pawrad(1)%mesh_size
!  call c_driver_spline(350,load_rr,load_data,msz,pawrad(1)%rad,sm_phi(1:msz,1,1))
!!!ENDDEBUG  
!=============== END OF HACK ===============

!------------------------------------
!loop over atoms in cell
 do iatom = 1, natom

   itypat = typat(iatom)
   nfgd = pawfgrtab(iatom)%nfgd ! number of points in the fine grid for this PAW sphere
   nnl = pawtab(itypat)%basis_size ! number of nl elements in PAW basis

   rhor_paw_in_loop = zero ! Zero the lopp dummy array

!  Division of fine grid points among processors
   if (nprocs==1) then ! Make sure everything runs with one proc
     write(message,'(a)') '  In denfgr - number of processors:     1'
     call wrtout(06,message,'COLL')
     write(message,'(a)') '  Calculation of PAW density done in serial'
     call wrtout(06,message,'COLL')
     write(message,'(a,I6)') '  Number of fine grid points:',nfgd
     call wrtout(06,message,'COLL')
     my_start_indx = 1
     my_end_indx = nfgd
   else ! Divide up the fine grid points among the processors
     if (rank==master) then
       write(message,'(a,I4)') '  In denfgr - number of processors: ',nprocs
       call wrtout(06,message,'COLL')
       write(message,'(a)') '  Calculation of PAW density done in parallel'
       call wrtout(06,message,'COLL')
       write(message,'(a,I6)') '  Number of fine grid points:',nfgd
       call wrtout(06,message,'COLL')
     end if
!    Divide the fine grid points among the processors
     delta = int(floor(real(nfgd)/real(nprocs)))
     remainder = nfgd-nprocs*delta
     my_start_indx = 1+rank*delta
     my_end_indx = (rank+1)*delta
!    Divide the remainder points among the processors
!    by shuffling indices
     if ((rank+1)>remainder) then
       my_start_indx = my_start_indx + remainder
       my_end_indx = my_end_indx + remainder
     else
       my_start_indx = my_start_indx + rank
       my_end_indx = my_end_indx + rank + 1
     end if
     if (prtvol>9) then
       write(message,'(a,I6)') '  My index Starts at: ',my_start_indx
       call wrtout(06,message,'PERS')
       write(message,'(a,I6)') '             Ends at: ',my_end_indx
       call wrtout(06,message,'PERS')
       write(message,'(a,I6)') '               # pts: ',my_end_indx+1-my_start_indx
       call wrtout(06,message,'PERS')
     end if
   end if

   write(message,'(a,I3,a,I3)') '  Entered loop for atom: ',iatom,' of:',natom
   call wrtout(06,message,'PERS')

!  obtain |r-R| values on fine grid
   allocate(nrm(nfgd))
   do ifgd=1, nfgd
     nrm(ifgd) = sqrt(dot_product(pawfgrtab(iatom)%rfgd(:,ifgd),pawfgrtab(iatom)%rfgd(:,ifgd)))
   end do ! these are the |r-R| values

!  compute Ylm for each r-R vector.
!  ----
   ipsang = 1 + (pawtab(itypat)%l_size - 1)/2 ! recall l_size=2*l_max+1
   allocate(ylm(ipsang*ipsang,nfgd))
   option = 1 ! compute Ylm(r-R) for vectors
   normchoice = 1 ! use computed norms of input vectors
   call initylmr(ipsang,normchoice,nfgd,nrm,option,pawfgrtab(iatom)%rfgd,ylm,ylmgr)

!  in order to do spline fits, the |r-R| data must be sorted
!  ----
   allocate(nrm_ifftsph(nfgd))
   nrm_ifftsph(:) = pawfgrtab(iatom)%ifftsph(:) ! copy of indices of points, to be rearranged by sort_dp
   call sort_dp(nfgd,nrm,nrm_ifftsph,tol8) ! sort the nrm points, keeping track of which goes where

!  now make spline fits of phi and tphi  onto the fine grid around the atom
!  ----
   allocate(phigrd(nfgd,nnl),tphigrd(nfgd,nnl))
   allocate(ypp(pawrad(itypat)%mesh_size),diag(pawrad(itypat)%mesh_size))


   do inl = 1, nnl

!!================================= chen hack in ===========================
!!    spline phi onto points
!!     ypp(:) = zero; diag(:) = zero; ybcbeg = zero; ybcend = zero;
!!     call spline(pawrad(itypat)%rad,pawtab(itypat)%phi(:,inl),pawrad(itypat)%mesh_size,ybcbeg,ybcend,ypp)
!!     call splint(pawrad(itypat)%mesh_size,pawrad(itypat)%rad,pawtab(itypat)%phi(:,inl),ypp,nfgd,nrm,phigrd(:,inl))
!!---------------       
!! ABINIT info:  real(dp), pointer :: phi(:,:)
!!      ! phi(mesh_size, basis_size)
!!      ! Gives, on the radial grid, the paw all electron wavefunctions

     write(message,'(a,I2,a,I2)')'chen/c_smooth_denfgr: sm_phi->phigrd, atom type:',itypat,' basis:',inl
     call wrtout(std_out,message,'COLL')

     ! inl: the basis index
     ypp(:) = zero; diag(:) = zero; ybcbeg = zero; ybcend = zero;
     msz=pawrad(itypat)%mesh_size
     call spline(pawrad(itypat)%rad,sm_phi(1:msz,inl,itypat),pawrad(itypat)%mesh_size,ybcbeg,ybcend,ypp)
     call splint(pawrad(itypat)%mesh_size,pawrad(itypat)%rad,sm_phi(1:msz,inl,itypat),ypp,nfgd,nrm,phigrd(:,inl))

!    next splint tphi onto points
     ypp(:) = zero; diag(:) = zero; ybcbeg = zero; ybcend = zero;
     call spline(pawrad(itypat)%rad,pawtab(itypat)%tphi(:,inl),pawrad(itypat)%mesh_size,ybcbeg,ybcend,ypp)
     call splint(pawrad(itypat)%mesh_size,pawrad(itypat)%rad,pawtab(itypat)%tphi(:,inl),ypp,nfgd,nrm,tphigrd(:,inl))

   end do ! end loop over nnl basis functions
   deallocate(ypp,diag)

!  loop over basis elements for this atom
!  because we have to store things like <phi|r'><r'|phi>-<tphi|r'><r'|tphi> at each point of the
!  fine grid, there is no integration, and hence no simplifications of the Y_lm's. That's why
!  we have to loop through the basis elements in exhaustive detail, rather than just a loop over
!  lmn2_size or something comparable.
!  ----
   if (prtvol>9) then
     write(message,'(a,I3)') '  Entering j-loop over basis elements for atom:',iatom
     call wrtout(06,message,'PERS')
   end if

   do jlmn=1,pawtab(itypat)%lmn_size

     if (prtvol>9) then
       write(message,'(2(a,I3))') '  Element:',jlmn,' of:',pawtab(itypat)%lmn_size
       call wrtout(06,message,'PERS')
     end if

     jl= psps%indlmn(1,jlmn,itypat)
     jm=psps%indlmn(2,jlmn,itypat)
     jlm = psps%indlmn(4,jlmn,itypat)
     jln=psps%indlmn(5,jlmn,itypat)
     j0lmn=jlmn*(jlmn-1)/2

     if (prtvol>9) then
       write(message,'(a,I3)') '  Entering i-loop for j:',jlmn
       call wrtout(06,message,'PERS')
     end if

     do ilmn=1,jlmn

       if (prtvol>9) then
         write(message,'(2(a,I3))') '    Element:',ilmn,' of:',jlmn
         call wrtout(06,message,'PERS')
       end if

       il= psps%indlmn(1,ilmn,itypat)
       im=psps%indlmn(2,ilmn,itypat)
       iln=psps%indlmn(5,ilmn,itypat)
       ilm = psps%indlmn(4,ilmn,itypat)
       klmn=j0lmn+ilmn

       if (prtvol>9) then
         write(message,'(a)') '    Entering loop over nonzero elems of rhoij'
         call wrtout(06,message,'PERS')
       end if

!      Loop over non-zero elements of rhoij
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         if (klmn==pawrhoij(iatom)%rhoijselect(irhoij)) then ! rho_ij /= 0 for this klmn

           do ifgd=my_start_indx, my_end_indx ! loop over fine grid points in current PAW sphere
             ifftsph = pawfgrtab(iatom)%ifftsph(ifgd) ! index of the point on the grid

!            have to retrieve the spline point to use since these were sorted
             do inrm=1, nfgd
               if(nrm_ifftsph(inrm) == ifftsph) exit ! have found nrm point corresponding to nfgd point
             end do ! now inrm is the index of the sorted nrm vector to use

!            avoid division by zero
             if(nrm(inrm) > zero) then
               rR = nrm(inrm) ! value of |r-R| in the following

!              recall that <r|phi>=u(r)*Slm(r^)/r
               phj  = phigrd(inrm,jln)*ylm(jlm,ifgd)/rR
               phi  = phigrd(inrm,iln)*ylm(ilm,ifgd)/rR
               tphj = tphigrd(inrm,jln)*ylm(jlm,ifgd)/rR
               tphi = tphigrd(inrm,iln)*ylm(ilm,ifgd)/rR

               do ispden=1,nspden
                 if (pawrhoij(iatom)%cplex == 1) then
                   rhor_paw_in_loop(ifftsph,ispden) = rhor_paw_in_loop(ifftsph,ispden) + &
&                   pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(irhoij,ispden)*&
&                   (phj*phi - tphj*tphi)

                 else
                   rhor_paw_in_loop(ifftsph,ispden) = rhor_paw_in_loop(ifftsph,ispden) + &
&                   pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(2*irhoij-1,ispden)*&
&                   (phj*phi - tphj*tphi)

                 end if ! end check on cplex rhoij
               end do ! end loop over nsdpen
             end if ! avoid |r-R| = 0

           end do ! end loop over nfgd
         end if ! end selection on rhoij /= 0
       end do ! end loop over non-zero rhoij
     end do ! end loop over ilmn atomic basis states
   end do ! end loop over jlmn atomic basis states

   deallocate(nrm,nrm_ifftsph,phigrd,tphigrd,ylm)

!  If we are working in parallel, collect all
!  contributions to the loop dummy array
   if (nprocs>1) then
     call xbarrier_mpi(spaceComm,ierr)
     call xsum_master(rhor_paw_in_loop,master,spaceComm,ierr)
     if (rank==master) rhor_paw = rhor_paw + rhor_paw_in_loop
     call xbarrier_mpi(spaceComm,ierr)
     write(message,'(a)') '  *** contributions to rhor_paw summed ***'
     call wrtout(std_out,message,'PERS')
   else
     rhor_paw = rhor_paw + rhor_paw_in_loop
   end if

 end do     ! Loop on atoms

!Add the plane-wave contribution \tilde{n} and remove \hat{n}
 rhor_paw = rhor_paw + rhor - nhat

 if (allocated(rhor_paw_in_loop)) deallocate(rhor_paw_in_loop)

!!=====================================
! chen: 
 deallocate(sm_phi)
!!===================================== 

 DBG_EXIT("COLL")

 end subroutine c_smooth_denfgr
!!***
