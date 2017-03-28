SUBROUTINE LoadPara(whoAmI,bLoadDen,bLoadVr,maxcount, & 
   bNCPP_PHI,penLambda,bAllowsym,optmethod,& 
   structure,stopTol,out_step)
!------------------------------------------
! Read in param.in file , load parameters
!------------------------------------------
   implicit none

   Integer :: & 
     whoAmI, &
     bLoadDen,     &
     bLoadVr,  &         
     bAllowSym, &
     optmethod, &
     bNCPP_PHI, &
     maxcount, & 
     out_step

   real(kind=8) ::  & 
     penLambda, & 
     stopTol

   character(len=500) :: & 
     msg, &
     structure
   
   !>>>>>>>>> FUNCTION BEGINS <<<<<<<<<!

   !Load parameters for job ----------------------------------
   open (unit=1111, access="sequential", action="read", file='../cluster/param.in', form="formatted", status="old");
   
   read (1111, *) bLoadDen,msg     ! 1: code will load density file, 0: will not load
   read (1111, *) bLoadVr,msg      ! 1: code will load potential file, 0: code will not load potential file
   read (1111, *) maxcount,msg     ! Max loop for our CG optimizer to find V_{bulk}, usually set to 200
   read (1111, *) bNCPP_PHI,msg ! 1: We do nonlocal inverting on the density
   read (1111, *) penLambda,msg    ! number: coeffcient for penalty function
   read (1111, *) bAllowsym,msg    ! 1: code allow symmetry (for bVatom=1 case), 0: code does not allow symmetry (for bVatom=0 case)
   read (1111, *) optmethod,msg    ! 1: cg, 2: BFGS
   read (1111, *) structure,msg    ! FCC or BCC or SC or DIA 
   read (1111, *) stopTol,msg      ! if ek change is smaller than stopTol twice successively, our job is done. If your system is metal just set this value to be 1.d-10 which is very small and will do nothing to the job, otherwise set to be 5.d-5 for one atom in a unit cell
   read (1111, *) out_step,msg     ! if ek change is smaller than stopTol twice successively, our job is done. If your system is metal just set this value to be 1.d-10 which is very small and will do nothing to the job, otherwise set to be 5.d-5 for one atom in a unit cell

   close (1111)

END SUBROUTINE LoadPara
