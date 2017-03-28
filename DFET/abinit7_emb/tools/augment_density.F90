module wf_files

  !This module reads in abinit wavefunctiondata files

  !According to the Abinit help file, the WF file contains ( see
  !http://www.abinit.org/documentation/helpfiles/for-v6.8/users/abinit_help.html#wavefctfile)
  
  !bantot=0                                    <-- counts over all bands
  !       index=0                                     <-- index for the wavefunction location
  !       do isppol=1,nsppol
  !        do ikpt=1,nkpt
  !         write(unit) npw,nspinor,nband                    <-- for each k point
  !         write(unit) kg(1:3,1:npw)                        <-- plane wave reduced coordinates
  !         write(unit) eigen(1+bantot:nband+bantot),        <-- eigenvalues for this k point
  !                     occ(1+bantot:nband+bantot)           <-- occupation numbers for this k point
  !         do iband=1,nband
  !          write(unit) (cg(ii+index),ii=1,2*npw*nspinor)   <-- wavefunction coefficients
  !         enddo                                            for a single band and k point
  !         bantot=bantot+nband
  !         index=index+2*npw*nspinor*nband
  !        enddo
  !       enddo

   type psp_type
     character*132 :: title
     double precision :: znuclpsp,zionpsp
     integer :: pspso,pspdat,pspcod,pspxc
     integer lmn_size
   END type psp_type

   type Header_type

     character*6 :: codvsn
     integer headform,fform
     integer bantot,date,intxc,ixc,natom,ngfft(3)
     integer nkpt,npsp,nspden,nspinor,nsppol,nsym,ntypat
     integer occopt,pertcase,usepaw
     integer usewvl, cplex
     integer lmax,lloc,mmax

     double precision :: acell(3),ecut,ecutdg,ecutsm,ecut_eff
     double precision :: qptn(3),rprimd(3,3),stmbias,tphysel,tsmear
     
     integer, allocatable :: istwfk(:) ! nkpt

     integer, allocatable :: nband(:)  ! nkpt x nsppol
         ! access as      ikpt+(isppol-1)*nkpt

     integer, allocatable :: npwarr(:) ! nkpt
     integer, allocatable :: so_psp(:) ! npsp
     integer, allocatable :: symafm(:) ! nsym
     integer, allocatable :: symrel(:,:,:) ! 3 x 3 x nsym
     integer, allocatable :: typat(:)  ! natom
     integer, allocatable :: nrhoijsel(:) ! nspden
     integer, allocatable :: rhoijselect(:,:) ! * x nspden
     real*8, allocatable :: kpt(:,:)          ! 3 x nkpt
     real*8, allocatable :: occ(:)            ! bantot
     real*8, allocatable :: tnons(:,:)        ! 3 x nsym
     real*8, allocatable :: znucltypat(:)     ! ntypat
     real*8, allocatable :: wtk(:)            ! nkpt

     real*8 :: residm,etotal,fermie
     real*8, allocatable :: xred(:,:) ! 3 x natom
     real*8, allocatable :: rhoij(:,:) ! * x nspden

     type (psp_type), allocatable :: psp(:)

   END type Header_type

   type wavefunction_type
      integer npw
      integer nspinor
      integer nband

      integer, pointer :: kg(:,:)

      complex*16, pointer :: cg(:)
      real*8, pointer :: eigen(:)
      real*8, pointer :: occ(:)
   END type wavefunction_type

   type state_data_type      
      
      integer nsppol
      integer nkpt

      integer, pointer :: kg(:,:)

      complex*16, pointer :: cg(:)
      real*8, pointer :: eigen(:)
      real*8, pointer :: occ(:)

      type (wavefunction_type), allocatable :: wft(:,:)
   end type state_data_type

contains

  subroutine read_header( file_id, header)
  
  implicit none
  
  type (Header_type) header

  integer, intent(IN) :: file_id

  integer ipsp
  
   write(6,*) "Reading in header data"
  
   read(file_id) header%codvsn,header%headform,header%fform 
  
   write(6,*) "forms:",header%headform, header%fform 
  
   read(file_id) header%bantot,header%date,header%intxc,header%ixc,header%natom, &
   & header%ngfft(1:3), header%nkpt, header%nspden, header%nspinor, &
   & header%nsppol,header%nsym,header%npsp,header%ntypat,header%occopt, &
   & header%pertcase,header%usepaw, header%ecut,header%ecutdg, &
   & header%ecutsm, header%ecut_eff, header%qptn(1:3), &
   & header%rprimd(1:3,1:3), header%stmbias, header%tphysel, &
   & header%tsmear, header%usewvl
  
   write(6,*) "Sizes read in ....", header%npsp, "additional records"
  
   allocate(header%istwfk(header%nkpt))
   allocate(header%nband(header%nkpt*header%nsppol))
   allocate(header%npwarr(header%nkpt))
   allocate(header%so_psp(header%npsp))
   allocate(header%symafm(header%nsym))
   allocate(header%symrel(3,3,header%nsym))
   allocate(header%typat(header%natom))
   allocate(header%kpt(3,header%nkpt))
   allocate(header%occ(header%bantot))
   allocate(header%tnons(3,header%nsym))
   allocate(header%znucltypat(header%ntypat))
   allocate(header%wtk(header%nkpt))
  
   read(file_id) header%istwfk(1:header%nkpt), &
   & header%nband(1:header%nkpt*header%nsppol), &
   & header%npwarr(1:header%nkpt),header%so_psp(1:header%npsp), &
   & header%symafm(1:header%nsym),header%symrel(1:3,1:3,1:header%nsym), &
   & header%typat(1:header%natom), header%kpt(1:3,1:header%nkpt), &
   & header%occ(1:header%bantot),header%tnons(1:3,1:header%nsym), &
   & header%znucltypat(1:header%ntypat),header%wtk(1:header%nkpt)
  
   allocate(header%psp(header%npsp))
  
   do ipsp=1,header%npsp
      ! (npsp lines, 1 for each pseudopotential ; npsp=ntypat, 
      ! except if alchemical pseudo-atoms)
   write(6,*) "next", ipsp
      read(file_id) header%psp(ipsp)%title, header%psp(ipsp)%znuclpsp, &
           & header%psp(ipsp)%zionpsp, header%psp(ipsp)%pspso, &
           & header%psp(ipsp)%pspdat, header%psp(ipsp)%pspcod, &
           & header%psp(ipsp)%pspxc, header%psp(ipsp)%lmn_size
   write(6,*) "crash", ipsp
   enddo
  
   allocate(header%xred(3,header%natom))
  
   !(in case of usepaw==0, final record: residm, coordinates, total energy, Fermi energy)
   read(file_id) header%residm, header%xred(1:3,1:header%natom),&
        &   header%etotal,header%fermie
  
  !(in case of usepaw==1, there are some additional records)
   if (header%usepaw==1)then
  
      write (6,*) "USEPAW NOT IMPLEMENTED"
      stop
  
  !  read(file_id)( pawrhoij(iatom)%nrhoijsel(1:nspden),iatom=1,natom), cplex, nspden
  !  read(file_id)((pawrhoij(iatom)%rhoijselect(1:      nrhoijsel(ispden),ispden),ispden=1,nspden),iatom=1,natom),&
  !&                 ((pawrhoij(iatom)%rhoijp     (1:cplex*nrhoijsel(ispden),ispden),ispden=1,nspden),iatom=1,natom)
   endif
  
  END subroutine read_header

  subroutine write_header( file_id, header)
  
  implicit none
  
  type (Header_type) header

  integer, intent(IN) :: file_id

  integer ipsp
  
   write(6,*) "Writing in header data"
  
   write(file_id) header%codvsn,header%headform,header%fform 
  
   write(file_id) header%bantot,header%date,header%intxc,header%ixc,header%natom, &
   & header%ngfft(1:3), header%nkpt, header%nspden, header%nspinor, &
   & header%nsppol,header%nsym,header%npsp,header%ntypat,header%occopt, &
   & header%pertcase,header%usepaw, header%ecut,header%ecutdg, &
   & header%ecutsm, header%ecut_eff, header%qptn(1:3), &
   & header%rprimd(1:3,1:3), header%stmbias, header%tphysel, &
   & header%tsmear, header%usewvl
  
   write(file_id) header%istwfk(1:header%nkpt), &
   & header%nband(1:header%nkpt*header%nsppol), &
   & header%npwarr(1:header%nkpt),header%so_psp(1:header%npsp), &
   & header%symafm(1:header%nsym),header%symrel(1:3,1:3,1:header%nsym), &
   & header%typat(1:header%natom), header%kpt(1:3,1:header%nkpt), &
   & header%occ(1:header%bantot),header%tnons(1:3,1:header%nsym), &
   & header%znucltypat(1:header%ntypat),header%wtk(1:header%nkpt)
  
   do ipsp=1,header%npsp
      ! (npsp lines, 1 for each pseudopotential ; npsp=ntypat, 
      ! except if alchemical pseudo-atoms)
      write(file_id) header%psp(ipsp)%title, header%psp(ipsp)%znuclpsp, &
           & header%psp(ipsp)%zionpsp, header%psp(ipsp)%pspso, &
           & header%psp(ipsp)%pspdat, header%psp(ipsp)%pspcod, &
           & header%psp(ipsp)%pspxc, header%psp(ipsp)%lmn_size
   enddo
  
   allocate(header%xred(3,header%natom))
  
   !(in case of usepaw==0, final record: residm, coordinates, total energy, Fermi energy)
   write(file_id) header%residm, header%xred(1:3,1:header%natom),&
        &   header%etotal,header%fermie
  
  !(in case of usepaw==1, there are some additional records)
   if (header%usepaw==1)then
  
      write (6,*) "USEPAW NOT IMPLEMENTED"
      stop
  
  !  write(file_id)( pawrhoij(iatom)%nrhoijsel(1:nspden),iatom=1,natom), cplex, nspden
  !  write(file_id)((pawrhoij(iatom)%rhoijselect(1:      nrhoijsel(ispden),ispden),ispden=1,nspden),iatom=1,natom),&
  !&                 ((pawrhoij(iatom)%rhoijp     (1:cplex*nrhoijsel(ispden),ispden),ispden=1,nspden),iatom=1,natom)
   endif
  
  END subroutine write_header

  subroutine free_header(header)
  
  implicit none
  
  type (Header_type) header
  
  integer ipsp
  
   deallocate(header%istwfk)
   deallocate(header%nband)
   deallocate(header%npwarr)
   deallocate(header%so_psp)
   deallocate(header%symafm)
   deallocate(header%symrel)
   deallocate(header%typat)
   deallocate(header%kpt)
   deallocate(header%occ)
   deallocate(header%tnons)
   deallocate(header%znucltypat)
   deallocate(header%wtk)
   deallocate(header%xred)
  
  END subroutine free_header
  
  subroutine read_unformatted( file_id, A, header )
    ! read either a density or an potential file, and write into A
    ! A has to be of size A(header%cplex*headery%ngfft(1)*headery%ngfft(2)*headery%ngfft(3),header%nspden)
    implicit none

    integer, intent(IN) :: file_id
  
    type (Header_type), intent(IN) :: header

    real*8, intent(OUT) :: A(header%cplex*header%ngfft(1)*header%ngfft(2)*header%ngfft(3),header%nspden)

    integer ispden, ir

    do ispden=1,header%nspden 
       read(file_id)(A(ir,ispden),ir=1,header%cplex*header%ngfft(1)*header%ngfft(2)*header%ngfft(3))       
    enddo

  end subroutine read_unformatted
  
  subroutine write_unformatted( file_id, A, header )
    ! write either a density or an potential file, and write into A
    ! A has to be of size A(header%cplex*header%ngfft(1)*header%ngfft(2)*header%ngfft(3),header%nspden)
    implicit none

    integer, intent(IN) :: file_id
  
    type (Header_type), intent(IN) :: header

    real*8, intent(IN) :: A(header%cplex*header%ngfft(1)*header%ngfft(2)*header%ngfft(3),header%nspden)

    integer ispden, ir

    do ispden=1,header%nspden 
       write(file_id)(A(ir,ispden),ir=1,header%cplex*header%ngfft(1)*header%ngfft(2)*header%ngfft(3))       
    enddo

  end subroutine write_unformatted
  
  subroutine write_density_file( lun, rho, header )
    ! write a density file in abinit format
    ! rho has to be of size A(header%cplex*header%ngfft(1)*header%ngfft(2)*header%ngfft(3),header%nspden)
    ! header has to be filled

    implicit none

    integer lun
  
    type (Header_type), intent(IN) :: header

    real*8, intent(IN) :: rho(header%cplex*header%ngfft(1)*header%ngfft(2)*header%ngfft(3),header%nspden)

    call write_header(31, header)
    call write_unformatted(31, rho, header)

  end subroutine write_density_file

END module wf_files

program tester

  use wf_files

  implicit none

  type (Header_type) :: header
  type (Header_type) :: header_big
  real*8 , allocatable :: rho(:,:)
  real*8 , allocatable :: rho_big(:,:)

  integer i,j,k
  integer i_big, j_big, k_big

  integer index, index_big

  integer ispden, nspden

  open(31,file="abinit.low.DEN", status="OLD", form="UNFORMATTED")
  call read_header(31, header)
  allocate(rho(header%cplex*&
       header%ngfft(1)*&
       header%ngfft(2)*&
       header%ngfft(3),header%nspden))

  call read_unformatted(31, rho, header)
  close(31)

  open(32,file="abinit.high.DEN", status="OLD", form="UNFORMATTED")
  call read_header(32, header_big)
  close(32)

  allocate(rho_big(header_big%cplex*&
       header_big%ngfft(1)*&
       header_big%ngfft(2)*&
       header_big%ngfft(3),header_big%nspden))

  if ( header_big%nspden .ne. header%nspden) then
     write (0,*) "Spdens not equal - ERROR "
     stop
  endif

  write(6,*) "Success!"

  ! now we can change the density 

  index_big = 0

  do ispden = 1, header_big%nspden
     do i_big = 1, header_big%ngfft(3)
        
        i = MOD(i_big-1, header%ngfft(3))

        do j_big = 1, header_big%ngfft(2)
           
           j = MOD(j_big-1, header%ngfft(2))
           index = header%ngfft(1) * ( j + i * header%ngfft(2) )
           
           do k_big = 1, header_big%ngfft(1)
           
              k = MOD(k_big-1, header%ngfft(1))+1
              index_big = index_big + 1

              rho_big(index_big,ispden) = rho(k + index,ispden)

           enddo
        enddo
     enddo
  enddo
  
  open(31,file="abinit.result.DEN", status="NEW", form="UNFORMATTED")
  
  call write_density_file(31,rho_big,header_big);

  close(31)

  deallocate(rho)

  call free_header(header)

end program tester
