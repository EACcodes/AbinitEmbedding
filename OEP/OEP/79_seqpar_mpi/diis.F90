! x#define DEBUG_DIIS
!Copyright (c) 2013, J. M. Dieterich
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
!
!    * Redistributions of source code must retain the above copyright
!      notice, this list of conditions and the following disclaimer.
!
!    * Redistributions in binary form must reproduce the above copyright
!      notice, this list of conditions and the following disclaimer in the
!      documentation and/or other materials provided with the distribution.
!
!    * All advertising materials mentioning features or use of this software
!      must display the following acknowledgement:
!
!      This product includes software developed at Princeton University (USA)
!      by J. M. Dieterich.
!
!    * Neither the name of Princeton University nor the names of its contributors
!      may be used to endorse or promote products derived from this software
!      without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE AUTHOR(S) ''AS IS'' AND ANY
!EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE AUTHOR(S) BE LIABLE FOR ANY
!DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

! This module should provide for a rather standalone version of the DIIS algorithm
! by P. Pulay. In difference to Pulay it uses the simpler error vector formulation
! as introduced by E. A. Carter (based on the delta between two successive solution 
! candidates p).
! PLEASE NOTE: WE REQUIRE BLAS/LAPACK FUNCTIONALITY FOR THIS!
! NAMELY, WE NEED DDOT, DPOTRF, DPOTRS TO COMPILE
! author: Johannes M. Dieterich
! version: 2013-04-08
! Please report bugs/problems/suggestions to jmd@ogolem.org
module diis_mod
    
    type diisBrain
        real(kind=8)::maxPercent = 0.1  ! TUNABLE: maximum allowed difference between the DIIS prediction and original point, default 50%
        integer::extraPoints = 4        ! TUNABLE: how many DIIS updates we must have before starting prediction (default: 4)
        integer::maxExtraPoints = 4     ! TUNABLE: how many DIIS updates we maximally use and keep in memory (default: 4)
        integer::problemDim = -1        ! TUNABLE: the problem dimension
        logical::autoNorm = .false.      ! TUNABLE: wether or not we norm the remaining parts of the B mat (-1 in the original paper) (default: true)
        integer::updateCounter          ! internal counter for updates to the DIIS memory (only to be touched by this module!)
        real(kind=8),dimension(:,:),allocatable::errorVecs  ! storage of error vectors (only to be touched by this module!)
        real(kind=8),dimension(:,:),allocatable::vectors    ! storage of vectors (only to be touched by this module!)
    end type diisBrain
    
    contains
    
    subroutine initializeDIIS(brain)
        implicit none
        
        type(diisBrain),intent(inout)::brain
        integer::allocstat
        
#ifdef DEBUG_DIIS
        write(*,*) "DIIS: Problem dimensionality for DIIS: ",brain%problemDim,brain%maxExtraPoints
#endif
        
        brain%updateCounter = 0
        allocate(brain%errorVecs(1:brain%problemDim,1:brain%maxExtraPoints),stat=allocstat)
        if(allocstat /= 0) then
            write(*,*) "ERROR: Failure to allocate error vector matrix in DIIS module. Stopping. ",allocstat
            stop
        endif
        allocate(brain%vectors(1:brain%problemDim,1:brain%maxExtraPoints+1),stat=allocstat)
        if(allocstat /= 0) then
            write(*,*) "ERROR: Failure to allocate vector matrix in DIIS module. Stopping. ",allocstat
            stop
        endif
    end subroutine initializeDIIS
    
    subroutine notifyDIIS(brain,p)
        implicit none
        
        type(diisBrain),intent(inout)::brain
        real(kind=8),dimension(brain%problemDim),intent(in)::p
        integer::i
        
        brain%updateCounter = brain%updateCounter +1
        
#ifdef DEBUG_DIIS
        write(*,*) "DIIS: Carrying out update ", brain%updateCounter
#endif
        
        if(brain%updateCounter > brain%maxExtraPoints+1) then
#ifdef DEBUG_DIIS
            write(*,*) "DIIS: Moving things around."
#endif
            ! we need to move a couple of things around
            do i=2,brain%maxExtraPoints
                brain%errorVecs(:,i-1) = brain%errorVecs(:,i)
            enddo
            do i=2,brain%maxExtraPoints+1
                brain%vectors(:,i-1) = brain%vectors(:,i)
            enddo
            brain%updateCounter = brain%maxExtraPoints+1
        endif
        
        brain%vectors(:,brain%updateCounter) = p(:)
        if(brain%updateCounter > 1) then
            brain%errorVecs(:,brain%updateCounter-1) = brain%vectors(:,brain%updateCounter)-brain%vectors(:,brain%updateCounter-1)
        endif
        
#ifdef DEBUG_DIIS
        write(*,*) "DIIS: Update counter is ",brain%updateCounter
        write(*,*) "DIIS: After update error vectors are "
        write(*,*) brain%errorVecs
        write(*,*) "DIIS: and vectors are "
        write(*,*) brain%vectors
#endif
        
    end subroutine notifyDIIS
    
    subroutine getDIISPrediction(brain,p)
        implicit none
        
        type(diisBrain),intent(inout)::brain
        real(kind=8),dimension(brain%problemDim),intent(inout)::p
        integer::i,j,allocstat,deallocstat,info
        real(kind=8)::tmp,tmp2,norm
        real(kind=8),dimension(:,:),allocatable::bmat
        real(kind=8),dimension(:),allocatable::xvec
        integer,dimension(:),allocatable::ipiv
        real(kind=8),external::ddot
        
#ifdef DEBUG_DIIS
        write(*,*) "DIIS: updateCounter and minimum extrapolation points ",brain%updateCounter,brain%extraPoints        
#endif
        if(brain%updateCounter < brain%extraPoints+1) return ! not yet enough knowledge in memory
        
        ! allocate some scratch space
        allocate(bmat(brain%updateCounter,brain%updateCounter),xvec(brain%updateCounter),ipiv(brain%updateCounter),stat=allocstat)
        if(allocstat /= 0) then
            write(*,*) "ERROR: Failure to allocate scratch space in DIIS module. Stopping. ",allocstat
            stop
        endif
        
        ! DIIS magic: build system of linear equations and solve
        bmat = 1.d200 ! (makes detection of errors easier)
        tmp2 = 0.0d0
        do i=1,brain%updateCounter-1
            do j=i,brain%updateCounter-1
                tmp = ddot(brain%problemDim,brain%errorVecs(:,i),1,brain%errorVecs(:,j),1)
                bmat(i,j) = tmp
                bmat(j,i) = tmp
                ! since the dot product is on real vectors, the matrix is symmetric
            enddo
            tmp2 = tmp2 + abs(bmat(i,i))
            !bmat(brain%updateCounter,i) = -1.0d0
            bmat(i,brain%updateCounter) = -1.0d0
        enddo
        
        norm = -1.0d0
        if(brain%autoNorm) then
            norm = -tmp2/(brain%updateCounter-1)
        endif
        
        do i=1,brain%updateCounter-1
            bmat(brain%updateCounter,i) = norm
            bmat(i,brain%updateCounter) = norm
        enddo
        bmat(brain%updateCounter,brain%updateCounter) = 0.0d0

#ifdef DEBUG_DIIS
        write(*,*) "DIIS: Finished building system of linear equations. Matrix is (only lower half relevant)"
        write(*,*) bmat
#endif
        
        xvec(:) = 0.0d0
        xvec(brain%updateCounter) = norm
        
        !factorize with LAPACK (bmat is symmetric)
        call dgetrf(brain%updateCounter,brain%updateCounter,bmat,brain%updateCounter,ipiv,info)
        if(info < 0) then
            write(*,*) "ERROR: DGETRF failed with exit code ",info
            stop
        else if(info > 0) then
#ifdef DEBUG_DIIS
        write(*,*) "DIIS: Linear dependency. Forgetting DIIS wisdom. ",info        
#endif
            ! linear dependency detected, forget some wisdom
            if(info < brain%updateCounter) then
                call forgetWisdomDIIS(brain,info)
                deallocate(bmat,xvec,ipiv,stat=deallocstat)
                if(deallocstat /= 0) then
                    write(*,*) "ERROR: Problem deallocating scratch arrays."
                    stop
                endif
                return
            else
                write(*,*) "ERROR: The last vector of bmat should never be a problem (as we cannot correct for that)."
                stop
            endif
        endif
        
#ifdef DEBUG_DIIS
        write(*,*) "DIIS: Finished factorization. Entering solution with LAPACK."        
#endif
        
        ! solve with lapack
        call dgetrs('N',brain%updateCounter,1,bmat,brain%updateCounter,ipiv,xvec,brain%updateCounter,info)
        if(info /= 0) then
            write(*,*) "ERROR: DGETRS failed with exit code ",info
            stop
        endif
        
#ifdef DEBUG_DIIS
        write(*,*) "DIIS: Finished solution. Solution vector:"
        write(*,*) xvec
#endif
        
        ! DIIS magic: use the xvec information to get a new p
        ! since we want to check on the fly, this loop ordering is a bit inefficient (to be changed!)
        do i=1,brain%problemDim
            tmp = 0.0d0
            do j=1,brain%updateCounter-1
                tmp = tmp + xvec(j)*brain%vectors(i,j+1)
            enddo
            if((p(i).eq.0.0d0) .and.(abs(tmp-p(i)) < brain%maxPercent)) then
                ! edge case that happens (criterion not percent anymore!)
#ifdef DEBUG_DIIS
                write(*,*) "DIIS: Accepting prediction ",i,p(i),tmp        
#endif
                p(i) = tmp
            else if(abs((tmp-p(i))/p(i)) < brain%maxPercent) then
#ifdef DEBUG_DIIS
                write(*,*) "DIIS: Accepting prediction ",i,p(i),tmp        
#endif                
                p(i) = tmp
            !else: new predicted step is too big. ignoring.
#ifdef DEBUG_DIIS
            else
                write(*,*) "DIIS: Discarding prediction ",i,p(i),tmp
#endif
            endif
        enddo

        deallocate(bmat, xvec, ipiv, stat = deallocstat)
        if (deallocstat /= 0) then
            write(*, *) "ERROR: Problem deallocating scratch arrays."
            stop
        endif

    end subroutine getDIISPrediction
    
    subroutine resetDIIS(brain)
        implicit none
        
        type(diisBrain),intent(inout)::brain
        
        brain%updateCounter = 0
    end subroutine resetDIIS
    
    subroutine forgetLastWisdomDIIS(brain)
        implicit none
        
        type(diisBrain),intent(inout)::brain
        
        brain%updateCounter = brain%updateCounter-1
    end subroutine forgetLastWisdomDIIS
    
    subroutine forgetWisdomDIIS(brain,which)
        implicit none
        
        type(diisBrain),intent(inout)::brain
        integer,intent(in)::which
        
        integer::i
        
        if(which < 1 .or. which > brain%maxExtraPoints) then
            write(*,*) "DIIS: Trying to remove an vector out of range: ",which
            stop
        endif
        
        ! but also copy
        do i=which+1,brain%updateCounter
            brain%errorVecs(:,i-1) = brain%errorVecs(:,i)
        enddo
        do i=which+2,brain%updateCounter
            brain%vectors(:,i-1) = brain%vectors(:,i)
        enddo
        
        brain%updateCounter = brain%updateCounter-1
        
    end subroutine forgetWisdomDIIS
    
    subroutine cleanBrain(brain)
        implicit none
        
        type(diisBrain),intent(inout)::brain
        integer::deallocstat
        
        brain%updateCounter = 0
        deallocate(brain%errorVecs,stat=deallocstat)
        if(deallocstat /= 0) then
            write(*,*) "ERROR: Failure to deallocate error vector matrix in DIIS module. Stopping. ",deallocstat
            stop
        endif
        deallocate(brain%vectors,stat=deallocstat)
        if(deallocstat /= 0) then
            write(*,*) "ERROR: Failure to deallocate vector matrix in DIIS module. Stopping. ",deallocstat
            stop
        endif
    end subroutine cleanBrain
end module diis_mod
