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

 subroutine denfgr(MPI_enreg,natom,nspden,nhat,ntypat,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,rhor,rhor_paw,prtvol,psps,typat)

 use c_hc_vars, only: & 
   in_paw 

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

! ************************************************************************

 DBG_ENTER("COLL")

 allocate(rhor_paw_in_loop(pawfgr%nfft,nspden))
 in_paw = -1

!Initialise rhor_paw
 rhor_paw = zero

!Initialise and check parallell execution
 call xcomm_init  (MPI_enreg,spaceComm)
 call xmpi_nproc  (nprocs,ierr)
 call xmaster_init(MPI_enreg,master)
 call xme_init    (MPI_enreg,rank)

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

!    spline phi onto points
     ypp(:) = zero; diag(:) = zero; ybcbeg = zero; ybcend = zero;
     call spline(pawrad(itypat)%rad,pawtab(itypat)%phi(:,inl),pawrad(itypat)%mesh_size,ybcbeg,ybcend,ypp)
     call splint(pawrad(itypat)%mesh_size,pawrad(itypat)%rad,pawtab(itypat)%phi(:,inl),ypp,nfgd,nrm,phigrd(:,inl))

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

!!================== chen =================
             in_paw(ifftsph)=1
!!=========== end of modification =========

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

 DBG_EXIT("COLL")

 end subroutine denfgr
!!***
