module threecg_mod

  !
  !   BFGS loop for OEP maximization
  !
  !   reading in parameter file param.in
  !
  ! (C) 2011 Florian Libisch, GNU public license
  !

  use defs_basis
  use defs_datatypes
  use m_wffile
  use defs_abitypes
  use defs_wvltypes
  use defs_parameters
  use defs_rectypes
  use m_electronpositron, only : electronpositron_type,electronpositron_calctype
  use m_rec
  use m_io_tools, only : flush_unit
  use m_paw_dmft, only: paw_dmft_type

  use auxiliary

  use interfaces_14_hidewrite
  use interfaces_18_timing
  use interfaces_27_toolbox_oop
  use interfaces_32_util
  use interfaces_41_geometry
  use interfaces_51_manage_mpi
  use interfaces_53_ffts
  use interfaces_56_recipspace
  use interfaces_56_xc
  use interfaces_57_iovars
  use interfaces_61_ionetcdf
  use interfaces_62_iowfdenpot
  use interfaces_65_nonlocal
  use interfaces_66_paw
  use interfaces_67_common
  use interfaces_68_recursion
  use interfaces_68_rsprc
  use interfaces_79_seqpar_mpi
  use interfaces_95_drive

  use m_xmpi
  use fourierhlp
  use get_XC_potential

  use vtoorbitalshifts

contains

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!********************************************************** Line 1
!
!                   December 2, 2011
!      ==========================================  
!                         
!
!                  Main Program THREECG
!              ============================    
!
!       Optimal Design with Composite Materials 
!     ===========================================
!
!   *** Application 4.4 from MINPACK-2
!   *** B.M. Averick, R.G. Carter, J.J. More, G.L. Xue
!   *** The MINPACK-2 test problem collection.
!   *** Argonne National Laboratory, Preprint MCS-P153-0692
!   *** June 1992. pp.36-39.
!
!  
!   THREECG is a subroutine dedicated to compute the minimizer of 
!   a differentiable function with a large number of variables. 
!   It is supposed that we have the algebraic expression of the
!   function and its gradient.    
!
!   THREECG implements an accelerated conjugate gradient algorithm 
!   with three terms, that at each iteration both the descent and 
!   the conjugacy conditions are guaranteed.  
!
!   THREECG implements an algorithm which is a modificatio of the
!   Hestenes and Stiefel conjugate gradient algorithm, or of that
!   of Hager and Zhang in such a way that the search direction is
!   descent and it satisfies the conjugacy condition. 
!   These properties are independent of the line search.
!
!
!   The algorithm is described in:
!   N. Andrei, A simple three-term conjugate gradient algorithm for
!   unconstrained optimization.
!   Journal of Computational and Applied Mathematics, vol.241 (2013)
!   pp. 19-29.
!     
!   Remark:
!   This paper is dedicated to the memory of Neculai Gheorghiu 
!   (1930-2009), my professor of Mathematical Analysis, Faculty of
!   Mathematics, "Alexandru Ioan Cuza" University, Iasi, Romania.
! 
!-------------------------------------------------------------------------
!
!   The algorithm is defined as:
!
!   x(k+1) = x(k) + alpha*d(k)
!
!   where alpha is computed by the Wolfe line search conditions.
!
!  
!   The search direction of this algorithm is computed as:
!
!   d(k+1) = -g(k+1) - delta(k)*s(k) - etha(k)*y(k)
!
!   where:
!                   
!                   ||y(k)||^2    s(k)T*g(k+1)   y(k)T*g(k+1)
!   delta(k) = (1 + ----------) * ------------ - ------------ ,
!                   y(k)T*s(k)     y(k)T*s(k)     y(k)T*s(k)   
!
!                  
!             s(k)T*g(k+1)
!   etha(k) = ------------.
!              y(k)T*s(k)                                 
!------------------------------------------------------------------------
!                                                          Neculai Andrei 
!                                      Research Institute for Informatics
!                           Center for Advanced Modeling and Optimization  
!                                                     and
!                                          Academy of Romanian Scientists
!*************************************************************************
!*
!       write(4,15)                                                      
! 15    format(4x,'***************************************************',/,&
!             4x,'* THREECG                                       ***',/,&                                                                   
!             4x,'* A simple three-term conjugate gradient        ***',/,&
!             4x,'* algorithm with Guaranteed Descent and         ***',/,& 
!             4x,'* Conjugacy conditions                          ***',/,&  
!             4x,'* --------------------------------------------- ***',/,&
!             4x,'* Project: FCGA                                 ***',/,&
!             4x,'* The Fastest Conjugate Gradient Algorithm      ***',/,&
!             4x,'*                                               ***',/,&
!             4x,'* Dr. Neculai Andrei                            ***',/,&
!             4x,'* Research Institute for Informatics            ***',/,&    
!             4x,'* Bucharest - Romania                           ***',/,&
!             4x,'***************************************************',/) 


!*            *** Conjugate Gradient Algorithms Project *** 
!*           ===============================================
!*
!*                                         
!*
!*
!*-------------------------------------------------------------------
!*                           Subroutine THREECG                                
!*                         ======================
!*                                   
!*
!*                             Neculai Andrei
!*                  Research Institute for Informatics
!*             Center for Advanced Modeling and Optimization
!*              8-10, Averescu Avenue, Bucharest 1, Romania
!*                        E-mail: nandrei@ici.ro
!*                           voice: 666.58.70
!*                                and
!*                     Academy of Romanian Scientists 
!*              Science and Information Technology Section              
!*            54, Splaiul Independentei, Bucharest 5, Romania
!*
!*
!*
!* /-----------------------------------------------------------------\
!* | THREECG is a subroutine dedicated to compute the minimizer of   |
!* | a differentiable function with a large number of variables.     |
!* |                                                                 |
!* | This subroutine is accompanied by subroutine "LineSearch" which |
!* | implements the Wolfe line search. Both these subroutines belong |
!* | to THREECG package.                                             |
!* |                                                                 |
!* | The user must provide a subroutine to evaluate the function     |
!* | value and its gradient in a point. The name of this subroutine  |
!* | is EVALFG.                                                      |
!* | The algebraic expression of the functions considered in this    |
!* | program, as well as their Fortran code are presented into the   |
!* | paper:                                                          |
!* | N. Andrei, "An unconstrained optimization test functions        |
!* | collection", Advanced Modeling and Optimization, vol.10, No.1,  |
!* | (2008) pp.147-161.                                              | 
!* | http://www.ici.ro/camo/journal/v10n1.htm                        |
!* |                                                                 |
!* | There are some facilities for the user to specify:              |
!* | 1) The termination criterion.                                   |
!* | 2) Convergence tolerance for gradient.                          |
!* | 3) Convergence tolerance for function value.                    |
!* | 4) Maximum number of iterations in LineSearch subroutine.       |
!* |                                                                 |
!* |-----------------------------------------------------------------|      
!* |                                                                 |
!* | The calling sequence of THREECG is:                             |
!* |                                                                 |
!* |   subroutine threecg(n,x,epsg,maxiter,maxfg,f,gnorm,stoptest,   |
!* |                      iter,irs,fgcnt,lscnt,angle, powell,        |
!* |                      nexp)                                      |
!* |                                                                 |
!* |Input parameters:                                                |
!* |=================                                                |
!* |n          (integer) number of variables.                        |
!* |x          (double)  starting guess, length n. On output         |
!* |                     contains the solution.                      |
!* |epsg       (double)  convergence tolerance for gradient.         |
!* |maxiter    (integer) maximum number of iterations.               |
!* |maxfg      (integer) maxumum number of function and its gradient |
!* |                     evaluations.                                |             
!* |stoptest = option parameter for selection of                     | 
!* |            stopping criterion:                                  |
!* |            if stoptest = 1 then consider the following test:    | 
!* |               if(ginf .le. epsg)                                |
!* |            if stoptest = 2 then consider the following test:    |
!* |               if(gnorm .le. epsg)                               |
!* |               where:                                            |
!* |               ginf  = infinite norm of gradient g(xk),          | 
!* |               gnorm = norm-2 of gradient g(xk).                 | 
!* |angle      (logical) parameter specifying the angle criterion of |
!* |                     restart.                                    |
!* |powell     (logical) parameter specifying the Powell criterion of|
!* |                     restart.                                    |   
!* |nexp       (integer) parameter specifying the number of the      |
!* |                     problem considered in a train of experiments|
!* |                                                                 |
!* |                                                                 |
!* |Output parameters:                                               |
!* |==================                                               |
!* |f          (double)  function value in final (optimal) point.    |
!* |gnorm      (double)  norm-2 of gradient in final point.          |
!* |iter       (integer) number of iterations to get the final point.|
!* |irs        (integer) number of restart iterations.               |
!* |fgcnt      (integer) number of function evaluations.             |
!* |lscnt      (integer) number of line searches.                    |
!* |-----------------------------------------------------------------|
!* |                                                                 |
!* |                                                                 |
!* |Calling subroutines:                                             |
!* |====================                                             | 
!* |Subroutine THREECG is calling two subroutines:                   |
!* |EVALFG     an user subroutine (function and gradient),           |
!* |LINESEARCH a package subroutine.                                 |
!* |                                                                 |   
!* |The user must supply a subroutine with the function and its      |
!* |gradient:                                                        |
!* |  call evalfg(n,x,fx,grad,nexp)                                  |
!* |where:                                                           |
!* |  n    (integer)  number of variables.                           |
!* |  x    (double)   the current iteration.                         |
!* |  fx   (double)   function value in point x.                     |
!* |  grad (double)   array with gradient of function in point x.    |
!* |  nexp (integer)  parameter specifying the number of the         |
!* |                  problem considered in a train of experiments.  |
!* |                                                                 |
!* |                                                                 |
!* |The Wolfe line search is implemented in the subroutine           |
!* |LINESEARCH, which belongs to the package.                        |
!* |The calling sequence of this subroutine is as follows:           |            
!* | call LineSearch (n,x,f,d,gtd,dnorm,alpha,xnew,fnew,gnew,sigma,  |
!* |                  fgcnt,lscnt,lsflag, nexp)                      |
!* |where:                                                           |
!* |  n      (integer)  number of variables.                         |
!* |  x      (double)   the current iteration.                       |
!* |  f      (double)   function value in current point.             |
!* |  d      (double)   array with search direction.                 |
!* |  gtd    (double)   scalar: grad'*d.                             |   
!* |  dnorm  (double)   2 norm of d.                                 |
!* |  alpha  (double)   step length (given by the LineSearch).       |
!* |  xnew   (double)   array with the new estimation of variables.  |
!* |  fnew   (double)   function value in xnew.                      |
!* |  gnew   (double)   array with gradient in xnew.                 |
!* |  sigma  (double)   parameter sigma in the second Wolfe line     |
!* |                    search condition. (input parameter)          |
!* |  fgcnt  (integer)  number of function evaluations.              |
!* |  lscnt  (integer)  number of line searches.                     |
!* |  lsflag (integer)  parameter for abnormal Line Search           |
!* |                    Termination. If the # of iterations in       |
!* |                    LineSearch is greater than 20 then lsflag=1  |
!* |  nexp   (integer)  parameter specifying the number of the       |
!* |                    problem considered in a train of experiments |
!* |                                                                 |
!* |-----------------------------------------------------------------|
!* |                                           Neculai Andrei, 2011  |
!* \-----------------------------------------------------------------/
!*           
!*
!*
!*********************************************************************

      subroutine threecg(n,x,epsg,maxiter,maxfg,f,gnorm,stoptest,&
                       iter,irs,fgcnt,lscnt,angle, powell,nexp, &
      afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop, &
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,nres2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params,&
       & susmat,symrec,taug,taur,tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)

      parameter(ia=1000000)

!*     SCALAR ARGUMENTS
      integer n,iter,irs,fgcnt,lscnt,maxiter,maxfg
      integer stoptest, nexp
      real(kind=8) epsg, f,gnorm     
      logical angle, powell 

      common /acca/epsm
      
!*     ARRAY ARGUMENTS
      real(kind=8) x(n)

!*     LOCAL SCALARS
      integer i,lsflag
      real(kind=8) fnew,alpha,gtg,& 
            dnorm,dnormprev, ginf,&
            gtd, gtgp, stg, yty, epsm,&
            yts, ytg, etha, delta, bn,sigma, lambda
      integer nx, ny
      
!*     LOCAL ARRAYS
      real(kind=8) xnew(ia),g(ia),gnew(ia),d(ia),&
             y(ia), s(ia)
      type(hdr_type),intent(inout) :: hdr
    type(dataset_type) :: dtset
    type(MPI_type),intent(inout) :: mpi_enreg
    type(electronpositron_type),pointer :: electronpositron
    type(energies_type) :: energies

    real(kind=8) :: tollist(12)

    real(kind=8) :: delta_e

    integer, intent(inout) :: quit    

    type(datafiles_type),intent(in) :: dtfil
    type(efield_type),intent(inout) :: dtefield
    type(pseudopotential_type),intent(in) :: psps

    type(wffile_type),intent(inout) :: wffnew,wffnow
    type(wvl_data),intent(inout) :: wvl

    integer n3xccc
    real(kind=8) xccc3d(n3xccc)

    type (run_parameter_type) run_params

    integer, intent(in) :: afford
    integer, intent(in) :: atindx(dtset%natom)
    integer, intent(in) :: atindx1(dtset%natom)
    integer, intent(in) :: dbl_nnsclo
    integer, intent(in) :: optres
    integer, intent(in) :: dielop
    integer, intent(in) :: dielstrt
    integer, intent(in) :: my_natom
    integer, intent(in) :: nkxc
    integer, intent(in) :: nattyp(psps%ntypat)

    integer, intent(in) :: ngfftdiel(18)
    integer, intent(in) :: mgfftdiel
    integer, intent(in) ::  gbound_diel(2*mgfftdiel+8,2)
    integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)

    integer, intent(in) :: computed_forces

    real*8, intent(in) :: gsqcut

    integer nfftf

    logical needgradient

    real*8 nvresid(nfftf,dtset%nspden)
    real*8 compch_fft
    real*8 dphase(3)
    real*8 grnl(3*dtset%natom)
    real*8 kxc(nfftf,nkxc)
    real*8 nhat(nfftf, dtset%nspden)

    real*8, intent(in) :: gmet(3,3),gprimd(3,3)

    real(kind=8), intent(in) :: rprimd(3,3)
    real(kind=8), intent(in) :: cpus

    real(kind=8), intent(in) :: ucvol

    integer, intent(in) :: pwind_alloc
    integer, intent(in) :: pwind(pwind_alloc,2,3)
    integer, intent(in) :: symrec(3,3,dtset%nsym)

    real(kind=8), intent(in) :: pwnsfac(2,pwind_alloc)
    real(kind=8), intent(inout) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(kind=8) residm, vxcavg

    real(kind=8), pointer :: rhog(:,:),rhor(:,:)
    real(kind=8), pointer :: taug(:,:),taur(:,:)

    real(kind=8), intent(in) :: rmet(3,3)

    real(kind=8) susmat(2,npwdiel*afford,dtset%nspden,&
         &npwdiel,dtset%nspden)

    real(kind=8), intent(inout) :: xred(3,dtset%natom)

    ! PAM set

    type(paw_dmft_type), intent(inout) :: paw_dmft
    type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
    type(pawang_type),intent(in) :: pawang
    type(pawfgr_type),intent(inout) :: pawfgr
    type(pawfgrtab_type) pawfgrtab(my_natom*psps%usepaw)
    type(pawrhoij_type), intent(inout) :: pawrhoij(my_natom*psps%usepaw)
    type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)

    real(dp), allocatable :: vh_o(:)

    ! ==============
    ! needed for scprqut

    real*8 diffor, favg(3), fcart(3,dtset%natom)

    ! ==============

    integer nfftdiel

    real(kind=8), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
    real*8, intent(in) :: ph1ddiel(2,3*(2*mgfftdiel+1)*dtset%natom*psps%usepaw)

    integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(inout) :: istep_mix
    integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
    integer, intent(in) :: npwdiel
    integer, intent(in) :: lmax_diel

    integer, intent(in) :: npwarr(dtset%nkpt)
    integer, intent(in) :: kg_diel(3,npwdiel)

    real(kind=8) nres2
    real(kind=8) dummy_nhatgr(1,1,1) ! dummy nhatgr

    real(kind=8), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(kind=8), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

    integer mcg
    real(kind=8), intent(inout) :: cg(2,mcg)

    real(kind=8), intent(in) :: ylm(dtset%mpw*dtset%mkmem,&
         & psps%mpsang*psps%mpsang*psps%useylm)
    real(kind=8), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,&
         & 3,psps%mpsang*psps%mpsang*psps%useylm)
    real(kind=8), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
    real(kind=8) vxctau(nfftf,dtset%nspden*dtset%usekden,4)

    real(kind=8) vhartr(nfftf),vpsp(nfftf)

    real (KIND=8) vtrial(nfftf, dtset%nspden)
    real*8 etotal
    real*8 dEdu_mean, dEdu_mean_up, dEdu_mean_down

    integer :: ix, iy ,iz, icount

    real(kind=8) strsxc(6) ! xc stress tensor dummy

    real*8 doti
    integer nfftotf

    integer optxc

    integer ktmp(3)
            
!* Initialization
    allocate(vh_o(nfftf))

!c                            Optimal Design with Composite Materials 
!c                            =======================================         
        lambda=0.008d0
        nx=1000
        ny=1000
      
      n5 = mod(n,5)
      n6 = n5 + 1
      
      iter   = 0
      irs    = 0
      fgcnt  = 0
      lscnt  = 0         
      ibeta  = 0
      itheta = 0
      
      call evalfg(n,x,f,g, nexp, &
      afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop, &
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,nres2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params,&
       & susmat,symrec,taug,taur,tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vh_o, vpsp, vtrial, vxctau)  
!      call dodcfg(nx,ny,x,f,g,'FG',lambda)
      fgcnt = fgcnt + 1

        gtg = 0.0d0
        do i = 1,n5
          d(i) = - g(i)
          gtg  = gtg + g(i) ** 2
        end do    
        do i = n6,n,5
          d(i)   = -g(i)
          d(i+1) = -g(i+1)
          d(i+2) = -g(i+2)
          d(i+3) = -g(i+3)
          d(i+4) = -g(i+4)
          gtg  = gtg + g(i)**2 + g(i+1)**2 + g(i+2)**2 + &
                               g(i+3)**2 + g(i+4)**2   
        end do
      gnorm = sqrt( gtg )

      gtd   = -gtg
      dnorm = gnorm                   
    

      if ( gnorm .gt. 0.0d0 ) then
          alpha = 1.0d0 / dnorm
      end if 

!* Initial value of parameter sigma. 

      sigma = 0.8d0       
    

    
!* --------------------------------   Main loop   --------------------    
!*====================================================================

110    continue
    
      
!*------------------------------------  STOP test section
!*                                      =================

      if(iter .eq. 0) go to 91

        if(stoptest .eq. 1) then
          ginf=dabs(g(1))
          do i=2,n5
            if(dabs(g(i))   .gt. ginf) ginf = dabs(g(i))
          end do   
          do i=n6,n,5
            if(dabs(g(i))   .gt. ginf) ginf = dabs(g(i))
            if(dabs(g(i+1)) .gt. ginf) ginf = dabs(g(i+1))
            if(dabs(g(i+2)) .gt. ginf) ginf = dabs(g(i+2))
            if(dabs(g(i+3)) .gt. ginf) ginf = dabs(g(i+3))
            if(dabs(g(i+4)) .gt. ginf) ginf = dabs(g(i+4))
          end do  
          if(ginf .le. epsg) go to 999    
        end if
!*
        if(stoptest .eq. 2) then
          if(gnorm .le. epsg) go to 999
        end if 
!*      
91    continue

      

!*---------------------------------- Increment iteration section
!*                                   ===========================

          iter = iter + 1  
          if(iter .gt. maxiter) go to 999         
          



!*---------------------------------- Line search section
!*                                   ===================                                           
!*
!* Determine the step length ALPHA and the new point XNEW, as well as
!* the function value in xnew, FNEW, and the gradient in xnew, GNEW.


          call LineSearch(n,x,f,d,gtd,dnorm,alpha,xnew,fnew,gnew,&
                         sigma,fgcnt,lscnt,lsflag, n5,n6,nexp, &
      afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop, &
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,nres2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params,&
       & susmat,symrec,taug,taur,tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vh_o, vpsp, vtrial, vxctau)

          if(fgcnt .gt. maxfg) go to 999

      
      
!*
!*---------------------------------- Acceleration section
!*                                   ====================
!*
!* (We use an/bn only if b is different from zero)
!* (an=gtd)
      
      bn=0.d0
      do i = 1,n5       
        bn = bn + (g(i)-gnew(i))*d(i)
      end do   
      do i = n6,n,5
        bn= bn + (g(i)  -gnew(i))*d(i)+&
                (g(i+1)-gnew(i+1))*d(i+1)+&
                (g(i+2)-gnew(i+2))*d(i+2)+&
                (g(i+3)-gnew(i+3))*d(i+3)+&
                (g(i+4)-gnew(i+4))*d(i+4)   
      end do
!*
      if(dabs(bn) .gt. epsm) then
          do i=1,n5
            xnew(i)   = x(i)   + (gtd/bn)*alpha*d(i)
          end do   
          do i=n6,n,5
            xnew(i)   = x(i)   + (gtd/bn)*alpha*d(i)
            xnew(i+1) = x(i+1) + (gtd/bn)*alpha*d(i+1)
            xnew(i+2) = x(i+2) + (gtd/bn)*alpha*d(i+2)
            xnew(i+3) = x(i+3) + (gtd/bn)*alpha*d(i+3)
            xnew(i+4) = x(i+4) + (gtd/bn)*alpha*d(i+4)
          end do
          call evalfg(n,xnew,fnew,gnew, nexp, &
      afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop, &
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,nres2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params,&
       & susmat,symrec,taug,taur,tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vh_o, vpsp, vtrial, vxctau)
!          call dodcfg(nx,ny,xnew,fnew,gnew,'FG',lambda)
          fgcnt = fgcnt + 1
          if(fgcnt .gt. maxfg) go to 999
      end if


!*
!*
!*---------------------------------- Prepare some scalar products
!*                                   ============================
!*              
         gtg = 0.d0
         gtgp= 0.d0          
         ytg = 0.d0 
         stg = 0.d0
         yts = 0.d0  
         yty = 0.d0 
       do i = 1,n5    
         s(i) = xnew(i) - x(i)
         y(i) = gnew(i) - g(i)
         ytg = ytg + gnew(i)*y(i)
         stg = stg + gnew(i)*s(i)
         yts = yts + y(i)*s(i)
         yty = yty + y(i)*y(i)
         gtgp = gtgp + gnew(i)*g(i)
           x(i) = xnew(i)
           g(i) = gnew(i)
         gtg  = gtg + g(i) * g(i)             
       end do  

       do i = n6,n,5
         s(i)   = xnew(i)   - x(i) 
         s(i+1) = xnew(i+1) - x(i+1)
         s(i+2) = xnew(i+2) - x(i+2)
         s(i+3) = xnew(i+3) - x(i+3)
         s(i+4) = xnew(i+4) - x(i+4)  
         y(i)   = gnew(i)   - g(i)
         y(i+1) = gnew(i+1) - g(i+1)
         y(i+2) = gnew(i+2) - g(i+2)
         y(i+3) = gnew(i+3) - g(i+3)
         y(i+4) = gnew(i+4) - g(i+4)
         ytg = ytg + gnew(i)*y(i)+gnew(i+1)*y(i+1)+gnew(i+2)*y(i+2)+&
                    gnew(i+3)*y(i+3)+gnew(i+4)*y(i+4)    
         stg = stg + gnew(i)*s(i)+gnew(i+1)*s(i+1)+gnew(i+2)*s(i+2)+&
                    gnew(i+3)*s(i+3)+gnew(i+4)*s(i+4)    
         yts = yts + y(i)*s(i)+y(i+1)*s(i+1)+y(i+2)*s(i+2)+&
                              y(i+3)*s(i+3)+y(i+4)*s(i+4)
         gtgp = gtgp + gnew(i)*g(i)+gnew(i+1)*g(i+1)+gnew(i+2)*g(i+2)+&
                     gnew(i+3)*g(i+3)+gnew(i+4)*g(i+4)
         yty = yty + y(i)*y(i)+y(i+1)*y(i+1)+y(i+2)*y(i+2)+&
                              y(i+3)*y(i+3)+y(i+4)*y(i+4)

           x(i)   = xnew(i)
           x(i+1) = xnew(i+1)
           x(i+2) = xnew(i+2)
           x(i+3) = xnew(i+3)
           x(i+4) = xnew(i+4)  
           g(i)   = gnew(i)
           g(i+1) = gnew(i+1)
           g(i+2) = gnew(i+2)
           g(i+3) = gnew(i+3)        
           g(i+4) = gnew(i+4)                                   
           gtg  = gtg + g(i)*g(i)+g(i+1)*g(i+1)+g(i+2)*g(i+2)+&
                       g(i+3)*g(i+3)+g(i+4)*g(i+4)      
       end do
!*    
          gnorm= sqrt( gtg )                  
!*
          f = fnew
          dnormprev = dnorm  
                
!*
!* -------------------------------- Delta and etha computation  
!*                                  ==========================
!*
!*
!*   Now delta computation:

         if(dabs(yts) .gt. epsm) then 
          delta = (1.d0+yty/yts)*stg/yts-ytg/yts
         else
          delta = 1.d0 
         end if  
!*
!*   Now etha computation
!*         
         if(yts .gt. epsm) then
           etha = stg/yts
         else
           etha  = 0.d0              
         end if                 
      
!*                        
!* --------------------------------- Direction computation
!*                                   =====================
!*
         dnorm = 0.0d0
         gtd   = 0.0d0
         do i = 1,n5
           d(i)  = -g(i) - delta*s(i) - etha*y(i) 
           dnorm =   dnorm + d(i) ** 2
           gtd   =   gtd + g(i) * d(i)  
         end do
         
         do i = n6,n,5  
           d(i)    = -g(i)   - delta * s(i)   - etha * y(i)     
           d(i+1)  = -g(i+1) - delta * s(i+1) - etha * y(i+1)
           d(i+2)  = -g(i+2) - delta * s(i+2) - etha * y(i+2)  
           d(i+3)  = -g(i+3) - delta * s(i+3) - etha * y(i+3)  
           d(i+4)  = -g(i+4) - delta * s(i+4) - etha * y(i+4)  

           dnorm   =   dnorm + d(i)**2+d(i+1)**2+d(i+2)**2+&
                                      d(i+3)**2+d(i+4)**2
           gtd     =   gtd + g(i)*d(i)+g(i+1)*d(i+1)+g(i+2)*d(i+2)+&
                                      g(i+3)*d(i+3)+g(i+4)*d(i+4) 
         end do

         dnorm = sqrt( dnorm )   
    
!*       
!*------------------------------------ END direction computation
!*             
!*   RESTART CRITERIA
!*   ================
!*
!* Angle Restart Test
!*
        if(angle) then
          if ( gtd .gt. -1.0d-03 * gnorm * dnorm ) then
              irs = irs + 1
              do i = 1,n
                d(i) = -g(i)
              end do
              dnorm =  gnorm
              gtd   = -gtg  
          end if        
        end if

!*
!* Beale-Powell restart test
!*             
        if(powell) then
          if(dabs(gtgp) .gt. 0.2d0*dabs(gtg)) then
              irs = irs + 1
              do i = 1,n5
                d(i) = -g(i) 
              end do     
              do i = n6,n,5
                d(i)   = -g(i)     
                d(i+1) = -g(i+1)   
                d(i+2) = -g(i+2)   
                d(i+3) = -g(i+3)   
                d(i+4) = -g(i+4)            
              end do              
              dnorm =   gnorm
              gtd   =  -gtg  
          end if            
        end if  
!*------------------------------------------ Prepare first trial
!*                                           of steplength          
!*                                           ===================
          
        if(dnorm .ne. 0.d0) then
          alpha = alpha * dnormprev / dnorm
        else
          alpha = 1.d0
        end if    
!*
      go to 110


!*------------------------------------------ End of main loop
!*                                           ================

999   continue
      
!*     
      deallocate(vh_o)
      
      return 

      end subroutine

!*-------------------------------------------- END THREECG subroutine





!c******************************************************************

      subroutine LineSearch (n,x,f,d,gtd,dnorm,alpha,xnew,fnew,gnew,&
                            sigma,fgcnt,lscnt,lsflag, n5,n6,nexp, &
      afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop, &
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,nres2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params,&
       & susmat,symrec,taug,taur,tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vh_o, vpsp, vtrial, vxctau)
       implicit none

!C     This is the one-dimensional line search used in CONMIN

!C     SCALAR ARGUMENTS
      integer n,fgcnt,lscnt,lsflag, nexp,n5,n6
      real(kind=8) f,gtd,dnorm,alpha,fnew

!C     ARRAY ARGUMENTS
      real(kind=8) x(n),d(n),xnew(n),gnew(n)

!C     LOCAL SCALARS
      integer i,lsiter
      real(kind=8) alphap,alphatemp,fp,dp,gtdnew,a,b,sigma
      real*8 lambda
      integer nx, ny
      type(hdr_type),intent(inout) :: hdr
    type(dataset_type) :: dtset
    type(MPI_type),intent(inout) :: mpi_enreg
    type(electronpositron_type),pointer :: electronpositron
    type(energies_type) :: energies

    real(kind=8) :: tollist(12)

    real(kind=8) :: delta_e

    integer, intent(inout) :: quit    

    type(datafiles_type),intent(in) :: dtfil
    type(efield_type),intent(inout) :: dtefield
    type(pseudopotential_type),intent(in) :: psps

    type(wffile_type),intent(inout) :: wffnew,wffnow
    type(wvl_data),intent(inout) :: wvl

    integer n3xccc
    real(kind=8) xccc3d(n3xccc)

    type (run_parameter_type) run_params

    integer, intent(in) :: afford
    integer, intent(in) :: atindx(dtset%natom)
    integer, intent(in) :: atindx1(dtset%natom)
    integer, intent(in) :: dbl_nnsclo
    integer, intent(in) :: optres
    integer, intent(in) :: dielop
    integer, intent(in) :: dielstrt
    integer, intent(in) :: my_natom
    integer, intent(in) :: nkxc
    integer, intent(in) :: nattyp(psps%ntypat)

    integer, intent(in) :: ngfftdiel(18)
    integer, intent(in) :: mgfftdiel
    integer, intent(in) ::  gbound_diel(2*mgfftdiel+8,2)
    integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)

    integer, intent(in) :: computed_forces

    real*8, intent(in) :: gsqcut

    integer nfftf

    logical needgradient

    real*8 nvresid(nfftf,dtset%nspden)
    real*8 compch_fft
    real*8 dphase(3)
    real*8 grnl(3*dtset%natom)
    real*8 kxc(nfftf,nkxc)
    real*8 nhat(nfftf, dtset%nspden)

    real*8, intent(in) :: gmet(3,3),gprimd(3,3)

    real(kind=8), intent(in) :: rprimd(3,3)
    real(kind=8), intent(in) :: cpus

    real(kind=8), intent(in) :: ucvol

    integer, intent(in) :: pwind_alloc
    integer, intent(in) :: pwind(pwind_alloc,2,3)
    integer, intent(in) :: symrec(3,3,dtset%nsym)

    real(kind=8), intent(in) :: pwnsfac(2,pwind_alloc)
    real(kind=8), intent(inout) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(kind=8) residm, vxcavg

    real(kind=8), pointer :: rhog(:,:),rhor(:,:)
    real(kind=8), pointer :: taug(:,:),taur(:,:)

    real(kind=8), intent(in) :: rmet(3,3)

    real(kind=8) susmat(2,npwdiel*afford,dtset%nspden,&
         &npwdiel,dtset%nspden)

    real(kind=8), intent(inout) :: xred(3,dtset%natom)

    ! PAM set

    type(paw_dmft_type), intent(inout) :: paw_dmft
    type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
    type(pawang_type),intent(in) :: pawang
    type(pawfgr_type),intent(inout) :: pawfgr
    type(pawfgrtab_type) pawfgrtab(my_natom*psps%usepaw)
    type(pawrhoij_type), intent(inout) :: pawrhoij(my_natom*psps%usepaw)
    type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)

    ! ==============
    ! needed for scprqut

    real*8 diffor, favg(3), fcart(3,dtset%natom)

    ! ==============

    real*8 vh_o(nfftf)

    integer nfftdiel

    real(kind=8), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
    real*8, intent(in) :: ph1ddiel(2,3*(2*mgfftdiel+1)*dtset%natom*psps%usepaw)

    integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(inout) :: istep_mix
    integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
    integer, intent(in) :: npwdiel
    integer, intent(in) :: lmax_diel

    integer, intent(in) :: npwarr(dtset%nkpt)
    integer, intent(in) :: kg_diel(3,npwdiel)

    real(kind=8) nres2
    real(kind=8) dummy_nhatgr(1,1,1) ! dummy nhatgr

    real(kind=8), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(kind=8), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

    integer mcg
    real(kind=8), intent(inout) :: cg(2,mcg)

    real(kind=8), intent(in) :: ylm(dtset%mpw*dtset%mkmem,&
         & psps%mpsang*psps%mpsang*psps%useylm)
    real(kind=8), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,&
         & 3,psps%mpsang*psps%mpsang*psps%useylm)
    real(kind=8), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
    real(kind=8) vxctau(nfftf,dtset%nspden*dtset%usekden,4)

    real(kind=8) vhartr(nfftf),vpsp(nfftf)

    real (KIND=8) vtrial(nfftf, dtset%nspden)
    real*8 etotal
    real*8 dEdu_mean, dEdu_mean_up, dEdu_mean_down

    integer :: ix, iy ,iz, icount

    real(kind=8) strsxc(6) ! xc stress tensor dummy

    real*8 doti
    integer nfftotf

    integer optxc

    integer ktmp(3)
    integer max$ls
      
      nexp=1
                 
      nx=1000
      ny=1000
      lambda=0.008d0


      lsflag = 0                                              
      
!* Maximum number of LineSearch is max$ls (now=6)     

      max$ls=8          
      
      alphap = 0.0d0
      fp     = f
      dp     = gtd

      do i = 1,n5
        xnew(i) = x(i) + alpha * d(i)
      end do     
      do i = n6,n,5
        xnew(i)   = x(i)   + alpha * d(i)
        xnew(i+1) = x(i+1) + alpha * d(i+1)
        xnew(i+2) = x(i+2) + alpha * d(i+2)
        xnew(i+3) = x(i+3) + alpha * d(i+3)
        xnew(i+4) = x(i+4) + alpha * d(i+4)
      end do    
!c1
         call evalfg(n,xnew,fnew,gnew, nexp, &
      afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop, &
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,nres2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params,&
       & susmat,symrec,taug,taur,tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vh_o, vpsp, vtrial, vxctau)
!         call dodcfg(nx,ny,xnew,fnew,gnew,'FG',lambda)
         fgcnt = fgcnt + 1

      gtdnew = 0.0d0
      do i = 1,n5
       gtdnew = gtdnew + gnew(i) * d(i)
      end do     
      do i = n6,n,5
       gtdnew = gtdnew + gnew(i)*d(i)+gnew(i+1)*d(i+1)+gnew(i+2)*d(i+2)+&
                        gnew(i+3)*d(i+3)+gnew(i+4)*d(i+4)
      end do
          
      lsiter = 0                                          

 10   if ( alpha * dnorm .gt. 1.0d-30 .and. lsiter .lt. max$ls .and. &
        .not. ( gtdnew .eq. 0.0d0 .and. fnew .lt. f ) .and. &
        ( ( fnew .gt. f + 1.0d-04 * alpha * gtd .or. &
        dabs( gtdnew / gtd ) .gt. sigma ) .or. ( lsiter .eq. 0 .and. &
        dabs( gtdnew / gtd ) .gt. 0.50d0 ) ) ) then

 20       if ( alpha * dnorm .gt. 1.0d-30 .and. fnew .gt. f .and. &
              gtdnew .lt. 0.0d0 ) then

              alpha = alpha / 3.0d0

              do i = 1,n5
                xnew(i) = x(i) + alpha * d(i)
              end do     
              do i = n6,n,5                
                xnew(i)   = x(i)   + alpha * d(i)
                xnew(i+1) = x(i+1) + alpha * d(i+1)
                xnew(i+2) = x(i+2) + alpha * d(i+2)
                xnew(i+3) = x(i+3) + alpha * d(i+3)
                xnew(i+4) = x(i+4) + alpha * d(i+4)
              end do  
!c2
                call evalfg(n,xnew,fnew,gnew, nexp, &
      afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop, &
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,nres2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params,&
       & susmat,symrec,taug,taur,tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vh_o, vpsp, vtrial, vxctau)
!                call dodcfg(nx,ny,xnew,fnew,gnew,'FG',lambda)
                fgcnt = fgcnt + 1

              gtdnew = 0.0d0
              do i = 1,n5
                gtdnew = gtdnew + gnew(i) * d(i)
              end do     
              do i = n6,n,5              
                gtdnew = gtdnew + gnew(i)*d(i)+gnew(i+1)*d(i+1)+ &
                                 gnew(i+2)*d(i+2)+gnew(i+3)*d(i+3)+ &
                                 gnew(i+4)*d(i+4)     
              end do

              alphap = 0.0d0
              fp     = f
              dp     = gtd

              goto 20

          end if
                 
          a = dp + gtdnew - 3.0d0 * ( fp - fnew ) / ( alphap - alpha )
          b = a ** 2 - dp * gtdnew

          if ( b .gt. 0.0d0 ) then
              b = sqrt( b )
          else
              b = 0.0d0
          end if

          alphatemp = alpha - ( alpha - alphap ) * ( gtdnew + b - a ) / &
                     ( gtdnew - dp + 2.0d0 * b )

          if ( gtdnew / dp .le. 0.0d0 ) then

              if ( 0.99d0 * max( alpha, alphap ) .lt. alphatemp .or. &
                 alphatemp .lt. 1.01d0 * min( alpha, alphap ) ) then
                  alphatemp = ( alpha + alphap ) / 2.0d0
              end if

          else

              if ( gtdnew .lt. 0.0d0 .and. &
                 alphatemp .lt. 1.01d0 * max( alpha, alphap ) ) then
                  alphatemp = 2.0d0 * max( alpha, alphap )
              end if

              if ( ( gtdnew .gt. 0.0d0 .and. &
                 alphatemp .gt. 0.99d0 * min( alpha, alphap ) ) .or. &
                 alphatemp .lt. 0.0d0 ) then
                  alphatemp = min( alpha, alphap ) / 2.0d0
              end if

          end if

          alphap = alpha
          fp     = fnew
          dp     = gtdnew

          alpha = alphatemp

          do i = 1,n5
            xnew(i) = x(i) + alpha * d(i)
          end do                           
          do i = n6,n,5
            xnew(i)   = x(i)   + alpha * d(i)
            xnew(i+1) = x(i+1) + alpha * d(i+1)
            xnew(i+2) = x(i+2) + alpha * d(i+2)
            xnew(i+3) = x(i+3) + alpha * d(i+3)
            xnew(i+4) = x(i+4) + alpha * d(i+4)
          end do   
!c3
            call evalfg(n,xnew,fnew,gnew, nexp, &
      afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop, &
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,nres2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params,&
       & susmat,symrec,taug,taur,tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vh_o, vpsp, vtrial, vxctau)
!            call dodcfg(nx,ny,xnew,fnew,gnew,'FG',lambda)
            fgcnt = fgcnt + 1

          gtdnew = 0.0d0
          do i = 1,n5
            gtdnew = gtdnew + gnew(i) * d(i)
          end do                              
          do i = n6,n,5
            gtdnew = gtdnew + gnew(i)*d(i)+gnew(i+1)*d(i+1)+ &
                             gnew(i+2)*d(i+2)+gnew(i+3)*d(i+3)+ &     
                             gnew(i+4)*d(i+4)
          end do

          lsiter = lsiter + 1

        goto 10

      end if

      if ( lsiter .ge. max$ls ) then
          lsflag = 1
      end if

      if ( lsiter .ne. 0 ) then
          lscnt = lscnt + 1
      end if

      return

      end subroutine
!*---------------------------------- End LineSearch subroutine

                                           
                                           

                                
                                
!*-----------------------------------------------------------
!*  Date created       : May 30, 1995
!*  Date last modified : May 30, 1995
!*
!*  Subroutine for execution time computation.
!*
!*-----------------------------------------------------------
!*
	subroutine exetim(tih,tim,tis,tic, tfh,tfm,tfs,tfc)
!*
	  integer*4 tih,tim,tis,tic
	  integer*4 tfh,tfm,tfs,tfc
!*
	  integer*4 ti,tf
	  integer*4 ch,cm,cs
	  data ch,cm,cs/360000,6000,100/
!*
	  ti=tih*ch+tim*cm+tis*cs+tic
	  tf=tfh*ch+tfm*cm+tfs*cs+tfc
	  tf=tf-ti
	  tfh=tf/ch
	  tf=tf-tfh*ch
	  tfm=tf/cm
	  tf=tf-tfm*cm
	  tfs=tf/cs
	  tfc=tf-tfs*cs
!*
	  return
	end subroutine
!*---------------------------------------------- End of EXETIM

      subroutine evalfg(n,xnew,fnew,gnew, nexp, &
      afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop, &
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,nres2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params,&
       & susmat,symrec,taug,taur,tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vh_o, vpsp, vtrial, vxctau)
      implicit none
      integer::n
      real(kind=8),dimension(:)::xnew,gnew
      real(kind=8)::fnew
      integer::nexp
    type(hdr_type),intent(inout) :: hdr
    type(dataset_type) :: dtset
    type(MPI_type),intent(inout) :: mpi_enreg
    type(electronpositron_type),pointer :: electronpositron
    type(energies_type) :: energies

    real(dp) vh_o(nfftf)

    real(kind=8) :: tollist(12)

    real(kind=8) :: delta_e

    integer, intent(inout) :: quit    

    type(datafiles_type),intent(in) :: dtfil
    type(efield_type),intent(inout) :: dtefield
    type(pseudopotential_type),intent(in) :: psps

    type(wffile_type),intent(inout) :: wffnew,wffnow
    type(wvl_data),intent(inout) :: wvl

    integer n3xccc
    real(kind=8) xccc3d(n3xccc)

    type (run_parameter_type) run_params

    integer, intent(in) :: afford
    integer, intent(in) :: atindx(dtset%natom)
    integer, intent(in) :: atindx1(dtset%natom)
    integer, intent(in) :: dbl_nnsclo
    integer, intent(in) :: optres
    integer, intent(in) :: dielop
    integer, intent(in) :: dielstrt
    integer, intent(in) :: my_natom
    integer, intent(in) :: nkxc
    integer, intent(in) :: nattyp(psps%ntypat)

    integer, intent(in) :: ngfftdiel(18)
    integer, intent(in) :: mgfftdiel
    integer, intent(in) ::  gbound_diel(2*mgfftdiel+8,2)
    integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)

    integer, intent(in) :: computed_forces

    real*8, intent(in) :: gsqcut

    integer nfftf

    real*8 nvresid(nfftf,dtset%nspden)
    real*8 compch_fft
    real*8 dphase(3)
    real*8 grnl(3*dtset%natom)
    real*8 kxc(nfftf,nkxc)
    real*8 nhat(nfftf, dtset%nspden)

    real*8, intent(in) :: gmet(3,3),gprimd(3,3)

    real(kind=8), intent(in) :: rprimd(3,3)
    real(kind=8), intent(in) :: cpus

    real(kind=8), intent(in) :: ucvol

    integer, intent(in) :: pwind_alloc
    integer, intent(in) :: pwind(pwind_alloc,2,3)
    integer, intent(in) :: symrec(3,3,dtset%nsym)

    real(kind=8), intent(in) :: pwnsfac(2,pwind_alloc)
    real(kind=8), intent(inout) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(kind=8) residm, vxcavg

    real(kind=8), pointer :: rhog(:,:),rhor(:,:)
    real(kind=8), pointer :: taug(:,:),taur(:,:)

    real(kind=8), intent(in) :: rmet(3,3)

    real(kind=8) susmat(2,npwdiel*afford,dtset%nspden,&
         &npwdiel,dtset%nspden)

    real(kind=8), intent(inout) :: xred(3,dtset%natom)

    ! PAM set

    type(paw_dmft_type), intent(inout) :: paw_dmft
    type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
    type(pawang_type),intent(in) :: pawang
    type(pawfgr_type),intent(inout) :: pawfgr
    type(pawfgrtab_type) pawfgrtab(my_natom*psps%usepaw)
    type(pawrhoij_type), intent(inout) :: pawrhoij(my_natom*psps%usepaw)
    type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)

    ! ==============
    ! needed for scprqut

    real*8 diffor, favg(3), fcart(3,dtset%natom)

    ! ==============

    integer nfftdiel

    real(kind=8), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
    real*8, intent(in) :: ph1ddiel(2,3*(2*mgfftdiel+1)*dtset%natom*psps%usepaw)

    integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(inout) :: istep_mix
    integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
    integer, intent(in) :: npwdiel
    integer, intent(in) :: lmax_diel

    integer, intent(in) :: npwarr(dtset%nkpt)
    integer, intent(in) :: kg_diel(3,npwdiel)

    real(kind=8) nres2
    real(kind=8) dummy_nhatgr(1,1,1) ! dummy nhatgr

    real(kind=8), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(kind=8), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

    integer mcg
    real(kind=8), intent(inout) :: cg(2,mcg)

    real(kind=8), intent(in) :: ylm(dtset%mpw*dtset%mkmem,&
         & psps%mpsang*psps%mpsang*psps%useylm)
    real(kind=8), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,&
         & 3,psps%mpsang*psps%mpsang*psps%useylm)
    real(kind=8), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
    real(kind=8) vxctau(nfftf,dtset%nspden*dtset%usekden,4)

    real(kind=8) vhartr(nfftf),vpsp(nfftf)

    real (KIND=8) vtrial(nfftf, dtset%nspden)
    real*8 etotal
    real*8 dEdu_mean, dEdu_mean_up, dEdu_mean_down

    integer :: ix, iy ,iz, icount

    real(kind=8) strsxc(6) ! xc stress tensor dummy

    real*8 doti
    integer nfftotf

    integer spaceComm

    integer optxc

    integer ktmp(3)
    
    logical,parameter::needgradient=.true.
    integer c,i,j,ispden,istep_os
    real*8 deps,dedeps,dedeps2,entropy_shift,e_fermie_shift
    real(kind=dp),allocatable :: doccde(:)
    real(dp), allocatable :: new_occ(:)
   
    integer ierr
 
    allocate(new_occ(dtset%mband*dtset%nkpt*dtset%nsppol))
    allocate(doccde(dtset%mband*dtset%nkpt*dtset%nsppol))
   
    spaceComm = mpi_enreg%comm_cell
 
    c = 1
    do j=1,dtset%nspden
       do i=1,nfftf
          vtrial(i,j)=xnew(c)
          c = c + 1
       enddo
    enddo
    
    call get_vtrial_energy_threecg(afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop,&
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,nres2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit, resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params, susmat,symrec,taug,taur,tollist,&
       & ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)

    vh_o = vhartr

       ! 1. Set up orbital shift equations:

       ! Reserve space for Orbital shifts 
       ! vxc (the KS potential)
       ! uxc (the real energy derivative dEXC/dphir)
       !  and the orb_shift vector g

       do ispden=1,min(dtset%nspden,2)
          run_params%orb_shift_vxc(:,ispden) = vtrial(:,ispden) -  &
                   vhartr(:) - vpsp(:)
          ! just the current trial XC potential 
       enddo
       if (dtset%nspden==4) run_params%orb_shift_vxc(:,3:4)=vtrial(:,3:4)      

90    continue           

       if (needgradient) then

          !! second call to vtorho, trap code in orbital_shifts
          !! ( vtorho -> vtowfk -> cgwf_orbital_shift )
          !! calculate the orbital shifts there

          !! reset the gradient
              
          run_params%dEdU = 0.d0

          write(2244,*) "Calling vtorho second time"

          istep_os = 1

          call vtorho_orbitalshifts(afford,atindx,atindx1, &
            &     cg,cpus,dbl_nnsclo,&
            &     dielop,dielstrt,dphase, dtfil,dtset,&
            &     eigen,electronpositron,energies,etotal,gbound_diel,&
            &     gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
            &     istep_os,istep_os,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,&
            &     mpi_enreg,my_natom,dtset%natom,nattyp,nfftf,nfftdiel,&
            &     ngfftdiel,nkxc,npwarr,npwdiel,nres2,&
            &     psps%ntypat,nvresid,occ,computed_forces,optres,paw_dmft,&
            &     paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,phnons,&
            &     phnonsdiel,ph1d,ph1ddiel,psps,pwind,pwind_alloc,pwnsfac,&
            &     resid,residm,rhog,rhor,rmet,rprimd,run_params,&
            &     symrec,taug,taur,ucvol,wffnew,wffnow,vtrial,vhartr, vh_o, wvl,xred,&
            &     ylm,ylmgr,ylmdiel)


#if defined HAVE_MPI
          if (mpi_enreg%nproc > 1) then
             call xsum_mpi(run_params%dEdu,spaceComm,ierr)
          end if
#endif
          ! correct for spin density

          if (dtset%nspden == 2) then             

             deps = 0.001d0

!******
             n = dtset%nkpt * dtset%mband

             eigen(1:n) = eigen(1:n) + deps
             eigen(n+1:2*n) = eigen(n+1:2*n) - deps

             ! use new_occ as placeholder for the changed occupations
             call newocc(doccde,eigen,entropy_shift,e_fermie_shift,&
                  dtset%spinmagntarget,dtset%mband,dtset%nband,dtset%nelect,&
                  dtset%nkpt,dtset%nspinor,dtset%nsppol,new_occ,&
                  dtset%occopt,dtset%prtvol,dtset%stmbias,dtset%tphysel,&
                  dtset%tsmear,dtset%wtk)

             icount = 0 
             dEdeps2 = 0.d0
             do iy=1,dtset%nsppol
                do iz=1,dtset%nkpt
                   do ix=1,dtset%mband
                      dEdeps2 = dEdeps2 + (new_occ(icount)*&
                           eigen(icount)-&
                           occ(icount)*eigen(icount))*dtset%wtk(iz)
                      icount = icount+1
                   enddo
                enddo
             enddo
          
             dEdeps2 = -dEdeps2 / (2.d0 * deps)

!******

             dEdeps = dEdeps2

!******
!
!             vtrial(:,1) = vtrial(:,1) + deps
!             vtrial(:,2) = vtrial(:,2) - deps
!
!             energies_shift%e_corepsp = energies%e_corepsp
!             energies_shift%e_ewald = energies%e_ewald
!
!             call get_vtrial_energy(afford, atindx, atindx1, &
!                  & dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop,&
!                  & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
!                  & electronpositron,energies_shift,etotal_shift,gbound_diel,&
!                  & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
!                  & istep_mix,kg,kg_diel,kxc,lmax_diel,mgfftdiel,mpi_enreg,&
!                  & nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
!                  & npwarr,npwdiel,res2,nvresid,occ,n3xccc, xccc3d,&


!              WHATS WITH RES2 / NRES2 ??


!                  & computed_forces,optres,paw_dmft,paw_ij,pawang,&
!                  & pawfgr,pawfgrtab,&
!                  & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
!                  & pwind,pwind_alloc,pwnsfac,quit, resid,residm,rhog,rhor,&
!                  & rmet,rprimd, susmat,symrec,taug,taur,tollist,ucvol,&
!                  & wffnew,wffnow,&
!                  & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial)
!             
!             vtrial(:,1) = vtrial(:,1) - deps
!             vtrial(:,2) = vtrial(:,2) + deps
!
!             dEdeps = (etotal_shift - etotal)/deps             
!
!             write(785,*) dedeps,dedeps2

!*******

             run_params%dEdu(:,1) = run_params%dEdu(:,1) + dEdeps * 0.2
             run_params%dEdu(:,2) = run_params%dEdu(:,2) - dEdeps * 0.2

             dEdu_mean_up = sum(run_params%dEdu(:,1))/dtset%nfft
             dEdu_mean_down = sum(run_params%dEdu(:,2))/dtset%nfft

          endif

          dEdu_mean = sum(run_params%dEdu(:,:))
          run_params%dEdu(:,:) = run_params%dEdu(:,:) - dEdu_mean
          write(2244,*) "Gradient max/min is ",&
               maxval(run_params%dEdu), minval(run_params%dEdu)
       
       endif

       c = 1
       do j=1,dtset%nspden
          do i=1,nfftf
             gnew(c) = run_params%dEdu(i,j)
             c = c + 1
          enddo
       enddo
       
       fnew = etotal
       
       deallocate(new_occ)
       deallocate(doccde)
       
    end subroutine evalfg
      
    subroutine get_vtrial_energy_threecg(afford, atindx, &
       & atindx1, dbl_nnsclo, cg, compch_fft, cpus, delta_e, dielop, &
       & dielstrt, dphase, dtfil, dtefield, dtset, eigen, &
       & electronpositron,energies,etotal,gbound_diel,&
       & gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
       & istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,mpi_enreg,&
       & my_natom,nattyp,nfftdiel,ngfftdiel,nfftf,nhat,nkxc,&
       & npwarr,npwdiel,nres2,nvresid,occ,n3xccc, xccc3d,&
       & computed_forces,optres,paw_dmft,paw_ij,pawang,pawfgr,pawfgrtab,&
       & pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
       & pwind,pwind_alloc,pwnsfac,quit,resid,residm,rhog,rhor,&
       & rmet,rprimd, run_params,&
       & susmat,symrec,taug,taur,tollist,ucvol,wffnew,wffnow,&
       & wvl,xred,ylm,ylmgr,ylmdiel, vhartr, vpsp, vtrial, vxctau)

    implicit none

    type(hdr_type),intent(inout) :: hdr
    type(dataset_type) :: dtset
    type(MPI_type),intent(inout) :: mpi_enreg
    type(electronpositron_type),pointer :: electronpositron
    type(energies_type) :: energies

    real(dp) :: tollist(12)

    real(dp) :: delta_e

    integer, intent(inout) :: quit    

    type(datafiles_type),intent(in) :: dtfil
    type(efield_type),intent(inout) :: dtefield
    type(pseudopotential_type),intent(in) :: psps

    type(wffile_type),intent(inout) :: wffnew,wffnow
    type(wvl_data),intent(inout) :: wvl

    integer n3xccc
    real(dp) xccc3d(n3xccc)

    type (run_parameter_type) run_params

    integer, intent(in) :: afford
    integer, intent(in) :: atindx(dtset%natom)
    integer, intent(in) :: atindx1(dtset%natom)
    integer, intent(in) :: dbl_nnsclo
    integer, intent(in) :: optres
    integer, intent(in) :: dielop
    integer, intent(in) :: dielstrt
    integer, intent(in) :: my_natom
    integer, intent(in) :: nkxc
    integer, intent(in) :: nattyp(psps%ntypat)

    integer, intent(in) :: ngfftdiel(18)
    integer, intent(in) :: mgfftdiel
    integer, intent(in) ::  gbound_diel(2*mgfftdiel+8,2)
    integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)

    integer, intent(in) :: computed_forces

    real*8, intent(in) :: gsqcut

    integer nfftf

    logical needgradient

    real*8 nvresid(nfftf,dtset%nspden)
    real*8 compch_fft
    real*8 dphase(3)
    real*8 grnl(3*dtset%natom)
    real*8 kxc(nfftf,nkxc)
    real*8 nhat(nfftf, dtset%nspden)

    real*8, intent(in) :: gmet(3,3),gprimd(3,3)

    real(dp), intent(in) :: rprimd(3,3)
    real(dp), intent(in) :: cpus

    real(dp), intent(in) :: ucvol

    integer, intent(in) :: pwind_alloc
    integer, intent(in) :: pwind(pwind_alloc,2,3)
    integer, intent(in) :: symrec(3,3,dtset%nsym)

    real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
    real(dp), intent(inout) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(dp) residm, vxcavg

    real(dp), pointer :: rhog(:,:),rhor(:,:)
    real(dp), pointer :: taug(:,:),taur(:,:)

    real(dp), intent(in) :: rmet(3,3)

    real(dp) susmat(2,npwdiel*afford,dtset%nspden,&
         &npwdiel,dtset%nspden)

    real(dp), intent(inout) :: xred(3,dtset%natom)

    ! PAM set

    type(paw_dmft_type), intent(inout) :: paw_dmft
    type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
    type(pawang_type),intent(in) :: pawang
    type(pawfgr_type),intent(inout) :: pawfgr
    type(pawfgrtab_type) pawfgrtab(my_natom*psps%usepaw)
    type(pawrhoij_type), intent(inout) :: pawrhoij(my_natom*psps%usepaw)
    type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)

    ! ==============
    ! needed for scprqut

    real*8 diffor, favg(3), fcart(3,dtset%natom)

    ! ==============

    integer nfftdiel

    real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),&
         &(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    real*8, intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
    real*8, intent(in) :: ph1ddiel(2,3*(2*mgfftdiel+1)*dtset%natom*psps%usepaw)

    integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym),2,&
         &    (dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
    integer, intent(inout) :: istep_mix
    integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
    integer, intent(in) :: npwdiel
    integer, intent(in) :: lmax_diel

    integer, intent(in) :: npwarr(dtset%nkpt)
    integer, intent(in) :: kg_diel(3,npwdiel)

    real(dp) nres2
    real(dp) dummy_nhatgr(1,1,1) ! dummy nhatgr

    real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)

    real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

    integer mcg
    real(dp), intent(inout) :: cg(2,mcg)

    real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,&
         & psps%mpsang*psps%mpsang*psps%useylm)
    real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,&
         & 3,psps%mpsang*psps%mpsang*psps%useylm)
    real(dp), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
    real(dp) vxctau(nfftf,dtset%nspden*dtset%usekden,4)

    real(dp) vhartr(nfftf),vh_o(nfftf),vpsp(nfftf)

    real (KIND=8) vtrial(nfftf, dtset%nspden)
    real*8 etotal
    real*8 dEdu_mean, dEdu_mean_up, dEdu_mean_down

    integer :: ix, iy ,iz, icount

    real(dp) strsxc(6) ! xc stress tensor dummy

    real*8 doti
    integer nfftotf

    integer optxc

    integer ktmp(3)

    write(2244,'(a,es16.8)')'scfcv: tolwfr = ', dtset%tolwfr
    write(2244,'(a,i4    )')'       nnsclo = ', dtset%nnsclo
    write(2244,'(a       )')'scfcv: calling vtorho() ... '

    ! First call to vtorho, 
    ! get the KS wavefunctions for current potential    

    call vtorho(afford,atindx,atindx1, cg,compch_fft,cpus,dbl_nnsclo,&
         &     dielop,dielstrt,dphase, dtefield,dtfil,dtset,&
         &     eigen,electronpositron,energies,etotal,gbound_diel,&
         &     gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
         &     istep_mix,istep_mix,kg,kg_diel,kxc,lmax_diel,mcg,mgfftdiel,&
         &     mpi_enreg,my_natom,dtset%natom,nattyp,nfftf,&
         &     nfftdiel,&
         &     ngfftdiel,nhat,nkxc,npwarr,npwdiel,nres2,&
         &     psps%ntypat,nvresid,occ,computed_forces,optres,paw_dmft,&
         &     paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,phnons,&
         &     phnonsdiel,ph1d,ph1ddiel,psps,pwind,pwind_alloc,pwnsfac,&
         &     resid,residm,rhog,rhor,rmet,rprimd,susmat,symrec,taug,&
         &     taur,ucvol,wffnew,wffnow,vtrial,wvl,xred,&
         &     ylm,ylmgr,ylmdiel,vxctau=vxctau)

       istep_mix = istep_mix + 1
       write(2244,'(a)       ')'scfcv: ------- vtorho() &
            &finished ------ '
       write(2244,'(a,es16.8)')'scfcv: residm = ', residm
       write(2244,*)''

       ! back up the obtained wavefunction and occupations

       run_params%e_fermie = energies%e_fermie

       if (dtset%nspden.eq.2) then
          write(2244,*) "Spin potential components:   Up: ", sum(vtrial(:,1))/nfftf
          write(2244,*) "Spin potential components: Down: ", sum(vtrial(:,2))/nfftf
       endif
       
       ! XC energy and potential
       optxc = 1! (exchange and hartree potential)
       
       call rhohxc(dtset,energies%e_xc,gsqcut,psps%usepaw,kxc,mpi_enreg,&
            nfftf, dtset%ngfft,nhat,0,dummy_nhatgr,0,0,0,dtset%nspden,&
            n3xccc,optxc,rhog,rhor,rprimd,strsxc,0,vhartr,&
            run_params%lda_pot,vxcavg,xccc3d,taug=taug,taur=taur)
       ! we get a vxc from rhohxc, we save that in lda_pot as it 
       ! might be needed for the XC calculation

       ! calculate the derivative of Exc with respect to orbital 
       ! for all orbitals

       call flush(6)

       call prepare_uxc(cg, mcg, mpi_enreg, gmet, gsqcut, rprimd, &
            ucvol, energies%e_xc, nfftf, dtset, npwarr, kg, occ, run_params)   

       ! Entropy energy
       energies%e_entropy = - energies%entropy * dtset%tsmear ! entropy energy

       ! potential energy
       ! this is correct also for spin because vpsp is spin-free and
       ! rhor(:,1) stores the trace of rho
       energies%e_localpsp = sum(rhor(:,1)*vpsp) * ucvol / dble(nfftf)

       nfftotf=dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)
       call dotprod_vn(1,rhor, energies%e_hartree, doti, mpi_enreg,&
            nfftf, nfftotf, 1, 1, vhartr, ucvol)
! we use nspden=1 above in any case because vhartr has only 
!               one spin component!
! we need to half to account for double-counting, 
!       see 67_common/energy.F90, line 347 
       energies%e_hartree=half*energies%e_hartree

!       call PSolver_hartree(dtset, energies%e_hartree, mpi_enreg, &
!            rhor, rprimd, vhartr)

       etotal   = energies%e_kinetic  + energies%e_hartree + energies%e_xc + &
            &            energies%e_localpsp + energies%e_corepsp + &
            &            energies%e_entropy  + energies%e_ewald + &
            &            energies%e_nonlocalpsp


       write(2244,&
            '(a,i3,a,2es12.4,a,es12.4,a,es18.10,a,es18.10)')& 
            'BFGS ITER: ',run_params%oep_iter, & 
            ' min/max(u):',minval(vtrial),maxval(vtrial), & 
            ' amp(u):',maxval(vtrial)-minval(vtrial),' etotal:',etotal
    ! display parts of total energy & total energy
       write(2244,'(a       )')'=========== results ==========='
       write(2244,'(a,es16.8)')'eigen      =', energies%e_eigenvalues
       write(2244,'(a       )')'=========== energies ==========='
       write(2244,'(a,es16.8)')'kinetic    =', energies%e_kinetic
       write(2244,'(a,es16.8)')'hartree    =', energies%e_hartree
       write(2244,'(a,es16.8)')'xc         =', energies%e_xc
       write(2244,'(a,es16.8)')'localpsp   =', energies%e_localpsp
       write(2244,'(a,es16.8)')'nonlocalpsp=', energies%e_nonlocalpsp
       write(2244,'(a,es16.8)')'e_corepsp  =', energies%e_corepsp
       write(2244,'(a,es16.8)')'entroy     =', energies%e_entropy
       write(2244,'(a,es16.8)')'ewald      =', energies%e_ewald
       write(2244,'(a,es16.8)')'scfcv: Etotal =',etotal

       diffor = 0.d0
       favg(:) = 0.d0
       fcart(:,:) = 0.d0

       call scprqt(2,cpus,delta_e,diffor,dtset,eigen,etotal,favg,&
            fcart,energies%e_fermie,dtfil%fnameabo_app_eig,&
            dtfil%filnam_ds(1),0,dtset%iscf,run_params%oep_iter,&
            dtset%kptns,0.d0,0,mpi_enreg,dtset%nband,&
            dtset%nkpt,run_params%max_oep_iter,occ,optres,&
            0,0,quit,nres2,resid,residm,&
            0,tollist,psps%usepaw,vxcavg,dtset%wtk,xred)    

  end subroutine get_vtrial_energy_threecg
!
!C                                                  Neculai Andrei
!C                                       Last line THREECG package
!C================================================================
!* Last Line
end module threecg_mod
