subroutine c_int_extpot_rho1(natom,ntypat,paw_an,pawrhoij,pawrad,rho1,extpot_lm,intg)
!===============================================================
! Do the 3D integral: \int v(r) n1(r) dr^3
! v and n1 are three-dimensional variable in the space.
! however n1 is in fact inside the PAW sphere, so we do 
! integral on radial mesh instead, we need to expand 
! v(r) and n1(r) over (l,m) moments, then the final 
! expression is :
!   \sum_(l,m) \int_0^rc [r**2 * v(|r|) n1(|r|)]
!
!  Chen Huang
!
!===============================================================
 
 use c_hc_vars

 use defs_basis
 use defs_datatypes
 use defs_abitypes

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat
 real(dp),intent(out) :: intg

!arrays
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(paw_an_type),intent(in) :: paw_an(natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 real(dp),intent(in) :: extpot_lm(hc_meshsz,hc_lsize**2,natom)
 real(dp),intent(in) :: rho1(hc_meshsz,hc_lsize**2,natom)

!Local variables ---------------------------------------
  character(len=500) :: message
  integer :: ia,itypat,ll,mm,mesh_size,idxl,l_max,lm_size,klm
  real(dp):: tmp,res(natom)
  real(dp),allocatable :: func(:)

! *************************************************************************

  write(message,'(a)')'(chen) enter c_int_extpot_rho1 '
  call wrtout(std_out,message,'COLL')

!==========================================================
!================ Big loop on atoms =======================
!==========================================================
 res = 0.d0
 do ia=1,natom
   
   lm_size=paw_an(ia)%lm_size
   itypat=pawrhoij(ia)%itypat
   mesh_size=pawrad(itypat)%mesh_size
   allocate(func(1:mesh_size))
    
   do klm=1,lm_size
     func=pawrad(itypat)%rad(1:mesh_size)**2 & 
         *rho1(1:mesh_size,klm,ia)*extpot_lm(1:mesh_size,klm,ia)
     call simp_gen(tmp,func,pawrad(itypat))
     res(ia)=res(ia)+tmp
   enddo ! loop over klm

   write(message,'(a,I4,a,ES16.8)')& 
     '(chen/c_int_extpot_rho1) atom:',ia,' ==> result:',res(ia)
   call wrtout(std_out,message,'COLL')

   deallocate(func)
 end do
!==========================================================
!=========== End loop on atoms ============================
!==========================================================
  intg=sum(res)

  write(message,'(a)')'(chen) exit c_int_extpot_rho1'
  call wrtout(std_out,message,'COLL')

end subroutine c_int_extpot_rho1
!!***

function  c_round_num(x)

  use defs_datatypes

  implicit none
  real(dp) :: x,y
  integer :: c_round_num

  c_round_num = floor(x)
  y = real(c_round_num,kind=dp)
  if (x-y>0.50d0) then
    c_round_num = c_round_num + 1
  endif
  return 

end function c_round_num

