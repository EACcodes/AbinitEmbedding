!======================================================
! This is a very simple interface to do spline
! 
!  Input:
!   old_num:  nodes number
!   old_x  :  nodes
!   old_y  :  values to be splined
!   num_new:  new number of points to calculate
!   new_x  :  x-coordinates to splined
!   new_y  :  splined values returned
!
! Chen Huang
!======================================================
subroutine c_driver_spline(old_num,old_x,old_y,num_new,new_x,new_y)
!
  use defs_basis
  use defs_datatypes
  use defs_abitypes

  implicit none
  integer :: old_num,num_new
  real(dp) :: old_x(old_num),old_y(old_num)
  real(dp) :: new_x(num_new),new_y(num_new)
  real(dp) :: ypp(old_num)

!!DEBUG
!  print *, "old_num=",old_num
!  print *, "num_new=",num_new
!  print *, "old_x=",old_x
!  print *, "old_y=",old_y
!  print *, "new_x=",new_x
!!END DEBUG

  call spline(old_x,old_y,old_num,0.d0,0.d0,ypp)

!!DEBUG
!  print *, "ypp="
!  print *, ypp
!!ENDDEBUG
  
  call splint(old_num,old_x,old_y,ypp,num_new,new_x,new_y)

!  print *, "new_y=",new_y
  
  return
end subroutine c_driver_spline

