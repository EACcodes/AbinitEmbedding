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
  
  subroutine write_density_file( filename, rho, header )
    ! write a density file in abinit format
    ! rho has to be of size A(header%cplex*header%ngfft(1)*header%ngfft(2)*header%ngfft(3),header%nspden)
    ! header has to be filled

    implicit none

    character, intent(IN) :: filename(:)
  
    type (Header_type), intent(IN) :: header

    real*8, intent(IN) :: rho(header%cplex*header%ngfft(1)*header%ngfft(2)*header%ngfft(3),header%nspden)

    open(31,file=filename, status="NEW", form="UNFORMATTED")
  
    call write_header(31, header)
    call write_unformatted(31, rho, header)

    close(31)

  end subroutine write_density_file

  subroutine read_wavefunction( file_id, states, header )
  
  implicit none
  
  integer, intent(IN) :: file_id
  
  type (state_data_type), intent(OUT), target :: states
  type (Header_type), intent(IN) :: header

  type (wavefunction_type), pointer:: cwf
  
  integer isppol, ikpt, iband
  
  integer ii

  integer index_cg, index_kg, bantot, cg_size, kg_size
     
  states%nsppol = header%nsppol
  states%nkpt = header%nkpt

  allocate(states%wft(states%nsppol, states%nkpt))
  
  index_kg = 0
  index_cg = 0
  bantot = 0

  write(6,*) "reading in WF for nsppol ", header%nsppol
  write(6,*) "and nkpt ", header%nkpt
  
  allocate(states%eigen(header%bantot))
  allocate(states%occ(header%bantot))
  
  cg_size = 0
  kg_size = 0
  do isppol=1,header%nsppol
     do ikpt=1,header%nkpt
        cg_size = cg_size + &
             & header%npwarr(ikpt)*header%nband(ikpt+(isppol-1)*header%nkpt)
        kg_size = kg_size + header%npwarr(ikpt)
     enddo
  enddo
  allocate(states%cg(cg_size*header%nspinor))
  allocate(states%kg(3,kg_size))
  
  do isppol=1,header%nsppol
     do ikpt=1,header%nkpt
        cwf => states%wft(isppol, ikpt)
  
        read(file_id) cwf%npw, cwf%nspinor, cwf%nband
  
        if (header%nband(ikpt+(isppol-1)*header%nkpt) .ne. cwf%nband) then
           write(6,*) "The band number is terribly wrong"
           stop         
        endif
        if (header%nspinor .ne. cwf%nspinor) then
           write(6,*) "Something spins terribly wrong"
  
           stop         
        endif
        if (header%npwarr(ikpt) .ne. cwf%npw) then
           write(6,*) "The plane wave number is terribly wrong"
           stop         
        endif
        write(6,*) "read cwf sizes",cwf%npw,cwf%nspinor,cwf%nband
  
        cwf%eigen => states%eigen(bantot+1:)
        cwf%occ => states%occ(bantot+1:)
        cwf%cg => states%cg(index_cg+1:)
        cwf%kg => states%kg(1:3,index_kg+1:)

        read(file_id) cwf%kg(1:3,1:cwf%npw)
        write(6,*) "read one kg"
        read(file_id) cwf%eigen(1:cwf%nband), cwf%occ(1:cwf%nband)
        write(6,*) "read eigen and occ"
  
        do iband=1,cwf%nband
           write(6,*) "read iband",iband, index_cg
           read(file_id) cwf%cg(1:cwf%npw*cwf%nspinor)
           index_cg = index_cg + cwf%npw*cwf%nspinor
        enddo
        index_kg = index_kg + cwf%npw
        bantot = bantot + cwf%nband
     enddo
  enddo
        
  END subroutine read_wavefunction
  
  subroutine free_wavefunction(states)

    implicit none

    type (state_data_type), target :: states

    deallocate( states%cg )
    deallocate( states%kg )
    deallocate( states%occ )
    deallocate( states%eigen )

  end subroutine free_wavefunction

END module wf_files

program tester

  use wf_files

  type (state_data_type) :: states
  type (Header_type) :: header

  open(21,file="t11o_WFK", status="OLD", form="UNFORMATTED")

  call read_header(21, header)
  call read_wavefunction(21, states, header)

  write(6,*) "Success!"

  write(6,*) "read in ", states%nsppol, 'x', states%nkpt , "state groups"

  write(6,*) "Eigenvalues", states%eigen(1)
  write(6,*) "Eigenvalues", states%eigen(2)

  call free_wavefunction(states)
  call free_header(header)

end program tester
