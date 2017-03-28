
MODULE auxiliary

  type T_Input_data
     integer np

     character*7, pointer::   names(:)
     character*22, pointer::  texts(:)
     logical, pointer::      parsed(:)

  END type T_Input_data

contains

  subroutine open_input(filename,logfile,id)

    implicit none

    character(LEN=*) filename
    character(LEN=*) logfile

    type (T_Input_data) id

    integer i

    integer NMax, indx

    NMax = 1000

    i = 1

    allocate(id%names(NMax))
    allocate(id%texts(NMax))
    allocate(id%parsed(NMax))

    open(20,file=filename,status='OLD')
    open(24,file=logfile)

10  continue

    read(20,160,END=20) id%names(i),id%texts(i)

    write(6,*) id%names(i),id%texts(i)
    id%parsed(i) = .FALSE.

    i = i + 1

    if (i.gt.NMax) then
       write(0,*) "Too long input file!"
       stop
    endif

    goto 10

20  continue

    close(20)

    id%np = i - 1

    return 



 55   format(A7,X,D22.15)

 65   format(25X,I5)
 66   format(8X,D22.15)
 67   format(26X,A4)
 68   format(20X,A10)
 69   format(10X,A20)


 75   format(A8,'                 ',I5)
 76   format(A8,D22.15)
 77   format(A8,'                  ',A4)
 78   format(A8,'            ',A10)
 79   format(A8,'  ',A20)

 80   format('rectp180tmn.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 81   format('rectp180tmn.',I2.2,'.',I2.2,'.',A4,'.absq')
 90   format('rectp180rmn.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 91   format('rectp180rmn.',I2.2,'.',I2.2,'.',A4,'.absq')
 92   format('rectp180trtot.',       I3.3,'.',A4,'.absq')

 180  format('recpt.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 181  format('recpt.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.absq')
 190  format('recpr.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 191  format('recpr.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.absq')
 192  format('recptot.mp.'       ,I2.2,'.',I3.3,'.',A4,'.absq')

 193  format('Bilddata.',A4,'.',I8.8,'.dat')
 194  format('pic.',A4,'.',I5.5,'.',I2.2,'.jpg') 
 195  format('Bildproj.',A4,'.',I8.8,'.dat')
 196  FORMAT('Eigenstate.',A4,'.',I4.4,'.dat')

 93   format('IC= ',I3,' ,T= ',E11.4,' ,R= ',E11.4,' Unitarity: ',E13.6)

 160  format(A7,'=',A22)

 165  format(18X,I5)
 166  format(D22.15)
 167  format(18X,A4)
 168  format(12X,A10)
 169  format(2X,A20)


 175  format(A7,'=            '     ,I10)
 176  format(A7,'=',              D22.15)
 177  format(A7,'=                  ',A4)
 178  format(A7,'=            '     ,A10)
 179  format(A7,'=  '               ,A20)

  END subroutine open_input

  subroutine free_ID(id)

    implicit none

    type (T_Input_data) id

    integer i

    do i=1,id%np
       if (.not.id%parsed(i)) then
          write(24,160) id%names(i),id%texts(i)
       endif
    enddo


    do i=1,id%np
       if (.not.id%parsed(i)) then
          write(0,*) 'Warning!'
          write(0,*) 'Unknown parameter detected! ', id%names(i)
          write(24,*) 'Warning!'
          write(24,*) 'Unknown parameter detected! ', id%names(i)
       endif
    enddo

    deallocate(id%names)
    deallocate(id%texts)
    deallocate(id%parsed)

 55   format(A7,X,D22.15)

 65   format(25X,I5)
 66   format(8X,D22.15)
 67   format(26X,A4)
 68   format(20X,A10)
 69   format(10X,A20)


 75   format(A8,'                 ',I5)
 76   format(A8,D22.15)
 77   format(A8,'                  ',A4)
 78   format(A8,'            ',A10)
 79   format(A8,'  ',A20)

 80   format('rectp180tmn.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 81   format('rectp180tmn.',I2.2,'.',I2.2,'.',A4,'.absq')
 90   format('rectp180rmn.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 91   format('rectp180rmn.',I2.2,'.',I2.2,'.',A4,'.absq')
 92   format('rectp180trtot.',       I3.3,'.',A4,'.absq')

 180  format('recpt.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 181  format('recpt.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.absq')
 190  format('recpr.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 191  format('recpr.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.absq')
 192  format('recptot.mp.'       ,I2.2,'.',I3.3,'.',A4,'.absq')

 193  format('Bilddata.',A4,'.',I8.8,'.dat')
 194  format('pic.',A4,'.',I5.5,'.',I2.2,'.jpg') 
 195  format('Bildproj.',A4,'.',I8.8,'.dat')
 196  FORMAT('Eigenstate.',A4,'.',I4.4,'.dat')

 93   format('IC= ',I3,' ,T= ',E11.4,' ,R= ',E11.4,' Unitarity: ',E13.6)

 160  format(A7,'=',A22)

 165  format(18X,I5)
 166  format(D22.15)
 167  format(18X,A4)
 168  format(12X,A10)
 169  format(2X,A20)


 175  format(A7,'=            '     ,I10)
 176  format(A7,'=',              D22.15)
 177  format(A7,'=                  ',A4)
 178  format(A7,'=            '     ,A10)
 179  format(A7,'=  '               ,A20)

  end subroutine free_ID

  function contains(id,name,indx)

    implicit none

    LOGICAL contains

    type (T_Input_data) id

    integer indx

    character*7 name

    integer i

    contains = .FALSE.

    indx = 0

    do i=1,id%np
       if (id%names(i).eq.name) then
          if (contains) then
             write(6,*) 'Warning: Parameter ',name, &
                  ' found twice in inputfile, last occurence used '
          endif
          contains = .TRUE.
          indx = i
       endif
    enddo

    RETURN

  END function contains

  subroutine demand_i_param(id,name,indx,ival)

    implicit none

    type (T_Input_data) id

    character*7 name

    integer indx

    integer ival, bla

    if (indx.eq.-1) then
       if (contains(id,name,bla)) then
          read(id%texts(bla),*) ival
          if (.NOT.id%parsed(bla)) write(24,175) id%names(bla),ival
          id%parsed(bla) = .TRUE.
       else
          write(0,*) "Parameter ",name," not found "
          stop
       endif
    else
       if (id%names(indx).eq.name) then
          read(id%texts(indx),*) ival
          if (.NOT.id%parsed(indx)) write(24,175) id%names(indx),ival
          id%parsed(indx) = .TRUE.
       else
          write(0,*) "Parameter ",name," expected at position ",indx
          stop
       endif
    endif



 55   format(A7,X,D22.15)

 65   format(25X,I5)
 66   format(8X,D22.15)
 67   format(26X,A4)
 68   format(20X,A10)
 69   format(10X,A20)


 75   format(A8,'                 ',I5)
 76   format(A8,D22.15)
 77   format(A8,'                  ',A4)
 78   format(A8,'            ',A10)
 79   format(A8,'  ',A20)

 80   format('rectp180tmn.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 81   format('rectp180tmn.',I2.2,'.',I2.2,'.',A4,'.absq')
 90   format('rectp180rmn.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 91   format('rectp180rmn.',I2.2,'.',I2.2,'.',A4,'.absq')
 92   format('rectp180trtot.',       I3.3,'.',A4,'.absq')

 180  format('recpt.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 181  format('recpt.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.absq')
 190  format('recpr.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 191  format('recpr.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.absq')
 192  format('recptot.mp.'       ,I2.2,'.',I3.3,'.',A4,'.absq')

 193  format('Bilddata.',A4,'.',I8.8,'.dat')
 194  format('pic.',A4,'.',I5.5,'.',I2.2,'.jpg') 
 195  format('Bildproj.',A4,'.',I8.8,'.dat')
 196  FORMAT('Eigenstate.',A4,'.',I4.4,'.dat')

 93   format('IC= ',I3,' ,T= ',E11.4,' ,R= ',E11.4,' Unitarity: ',E13.6)

 160  format(A7,'=',A22)

 165  format(18X,I5)
 166  format(D22.15)
 167  format(18X,A4)
 168  format(12X,A10)
 169  format(2X,A20)


 175  format(A7,'=            '     ,I10)
 176  format(A7,'=',              D22.15)
 177  format(A7,'=                  ',A4)
 178  format(A7,'=            '     ,A10)
 179  format(A7,'=  '               ,A20)

  END subroutine demand_i_param

  function check_i_param(id,name,indx,ival)

    implicit none

    type (T_Input_data) id

    character*7 name

    integer indx

    integer ival

    logical check_i_param

    if (contains(id,name,indx)) then
       read(id%texts(indx),*) ival
       if (.NOT.id%parsed(indx)) write(24,175) id%names(indx),ival
       id%parsed(indx) = .TRUE.
       check_i_param = .TRUE.
    else
       check_i_param = .FALSE.
    endif



 55   format(A7,X,D22.15)

 65   format(25X,I5)
 66   format(8X,D22.15)
 67   format(26X,A4)
 68   format(20X,A10)
 69   format(10X,A20)


 75   format(A8,'                 ',I5)
 76   format(A8,D22.15)
 77   format(A8,'                  ',A4)
 78   format(A8,'            ',A10)
 79   format(A8,'  ',A20)

 80   format('rectp180tmn.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 81   format('rectp180tmn.',I2.2,'.',I2.2,'.',A4,'.absq')
 90   format('rectp180rmn.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 91   format('rectp180rmn.',I2.2,'.',I2.2,'.',A4,'.absq')
 92   format('rectp180trtot.',       I3.3,'.',A4,'.absq')

 180  format('recpt.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 181  format('recpt.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.absq')
 190  format('recpr.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 191  format('recpr.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.absq')
 192  format('recptot.mp.'       ,I2.2,'.',I3.3,'.',A4,'.absq')

 193  format('Bilddata.',A4,'.',I8.8,'.dat')
 194  format('pic.',A4,'.',I5.5,'.',I2.2,'.jpg') 
 195  format('Bildproj.',A4,'.',I8.8,'.dat')
 196  FORMAT('Eigenstate.',A4,'.',I4.4,'.dat')

 93   format('IC= ',I3,' ,T= ',E11.4,' ,R= ',E11.4,' Unitarity: ',E13.6)

 160  format(A7,'=',A22)

 165  format(18X,I5)
 166  format(D22.15)
 167  format(18X,A4)
 168  format(12X,A10)
 169  format(2X,A20)


 175  format(A7,'=            '     ,I10)
 176  format(A7,'=',              D22.15)
 177  format(A7,'=                  ',A4)
 178  format(A7,'=            '     ,A10)
 179  format(A7,'=  '               ,A20)

  END function check_i_param

  function check_r_param(id,name,indx,rval)

    implicit none

    type (T_Input_data) id

    character*7 name

    integer indx

    real*8 rval

    logical check_r_param

    if (contains(id,name,indx)) then
       read(id%texts(indx),*) rval
       if (.NOT.id%parsed(indx)) write(24,176) id%names(indx),rval
       id%parsed(indx) = .TRUE.
       check_r_param = .TRUE.
    else
       check_r_param = .FALSE.
    endif



 55   format(A7,X,D22.15)

 65   format(25X,I5)
 66   format(8X,D22.15)
 67   format(26X,A4)
 68   format(20X,A10)
 69   format(10X,A20)


 75   format(A8,'                 ',I5)
 76   format(A8,D22.15)
 77   format(A8,'                  ',A4)
 78   format(A8,'            ',A10)
 79   format(A8,'  ',A20)

 80   format('rectp180tmn.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 81   format('rectp180tmn.',I2.2,'.',I2.2,'.',A4,'.absq')
 90   format('rectp180rmn.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 91   format('rectp180rmn.',I2.2,'.',I2.2,'.',A4,'.absq')
 92   format('rectp180trtot.',       I3.3,'.',A4,'.absq')

 180  format('recpt.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 181  format('recpt.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.absq')
 190  format('recpr.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 191  format('recpr.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.absq')
 192  format('recptot.mp.'       ,I2.2,'.',I3.3,'.',A4,'.absq')

 193  format('Bilddata.',A4,'.',I8.8,'.dat')
 194  format('pic.',A4,'.',I5.5,'.',I2.2,'.jpg') 
 195  format('Bildproj.',A4,'.',I8.8,'.dat')
 196  FORMAT('Eigenstate.',A4,'.',I4.4,'.dat')

 93   format('IC= ',I3,' ,T= ',E11.4,' ,R= ',E11.4,' Unitarity: ',E13.6)

 160  format(A7,'=',A22)

 165  format(18X,I5)
 166  format(D22.15)
 167  format(18X,A4)
 168  format(12X,A10)
 169  format(2X,A20)


 175  format(A7,'=            '     ,I10)
 176  format(A7,'=',              D22.15)
 177  format(A7,'=                  ',A4)
 178  format(A7,'=            '     ,A10)
 179  format(A7,'=  '               ,A20)

  END function check_r_param

  subroutine demand_r_param(id,name,indx,rval)

    implicit none

    type (T_Input_data) id

    character*7 name

    integer indx,bla

    real*8 rval

    if (indx.eq.-1) then
       if (contains(id,name,bla)) then
          read(id%texts(bla),*) rval
          if (.NOT.id%parsed(bla)) write(24,176) id%names(bla),rval
          id%parsed(bla) = .TRUE.
       else
          write(0,*) "Parameter ",name," not found "
          stop
       endif
    else
       if (id%names(indx).eq.name) then
          read(id%texts(indx),*) rval
          if (.NOT.id%parsed(indx)) write(24,176) id%names(indx),rval
          id%parsed(indx) = .TRUE.
       else
          write(0,*) "Parameter ",name," expected at position ",indx
          stop
       endif
    endif



 55   format(A7,X,D22.15)

 65   format(25X,I5)
 66   format(8X,D22.15)
 67   format(26X,A4)
 68   format(20X,A10)
 69   format(10X,A20)


 75   format(A8,'                 ',I5)
 76   format(A8,D22.15)
 77   format(A8,'                  ',A4)
 78   format(A8,'            ',A10)
 79   format(A8,'  ',A20)

 80   format('rectp180tmn.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 81   format('rectp180tmn.',I2.2,'.',I2.2,'.',A4,'.absq')
 90   format('rectp180rmn.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 91   format('rectp180rmn.',I2.2,'.',I2.2,'.',A4,'.absq')
 92   format('rectp180trtot.',       I3.3,'.',A4,'.absq')

 180  format('recpt.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 181  format('recpt.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.absq')
 190  format('recpr.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 191  format('recpr.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.absq')
 192  format('recptot.mp.'       ,I2.2,'.',I3.3,'.',A4,'.absq')

 193  format('Bilddata.',A4,'.',I8.8,'.dat')
 194  format('pic.',A4,'.',I5.5,'.',I2.2,'.jpg') 
 195  format('Bildproj.',A4,'.',I8.8,'.dat')
 196  FORMAT('Eigenstate.',A4,'.',I4.4,'.dat')

 93   format('IC= ',I3,' ,T= ',E11.4,' ,R= ',E11.4,' Unitarity: ',E13.6)

 160  format(A7,'=',A22)

 165  format(18X,I5)
 166  format(D22.15)
 167  format(18X,A4)
 168  format(12X,A10)
 169  format(2X,A20)


 175  format(A7,'=            '     ,I10)
 176  format(A7,'=',              D22.15)
 177  format(A7,'=                  ',A4)
 178  format(A7,'=            '     ,A10)
 179  format(A7,'=  '               ,A20)

  END subroutine demand_r_param

  function check_bool(id,name,indx)

    implicit none

    logical check_bool

    type (T_Input_data) id

    character*7 name

    integer indx,bla

    if (contains(id,name,indx)) then
       write(24,177) id%names(indx),'true'
       id%parsed(indx) = .TRUE.
       check_bool = .TRUE.
    else
       check_bool = .FALSE.
    endif



 55   format(A7,X,D22.15)

 65   format(25X,I5)
 66   format(8X,D22.15)
 67   format(26X,A4)
 68   format(20X,A10)
 69   format(10X,A20)


 75   format(A8,'                 ',I5)
 76   format(A8,D22.15)
 77   format(A8,'                  ',A4)
 78   format(A8,'            ',A10)
 79   format(A8,'  ',A20)

 80   format('rectp180tmn.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 81   format('rectp180tmn.',I2.2,'.',I2.2,'.',A4,'.absq')
 90   format('rectp180rmn.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 91   format('rectp180rmn.',I2.2,'.',I2.2,'.',A4,'.absq')
 92   format('rectp180trtot.',       I3.3,'.',A4,'.absq')

 180  format('recpt.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 181  format('recpt.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.absq')
 190  format('recpr.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.cmpl')
 191  format('recpr.mp.',I2.2,'.',I2.2,'.',I2.2,'.',A4,'.absq')
 192  format('recptot.mp.'       ,I2.2,'.',I3.3,'.',A4,'.absq')

 193  format('Bilddata.',A4,'.',I8.8,'.dat')
 194  format('pic.',A4,'.',I5.5,'.',I2.2,'.jpg') 
 195  format('Bildproj.',A4,'.',I8.8,'.dat')
 196  FORMAT('Eigenstate.',A4,'.',I4.4,'.dat')

 93   format('IC= ',I3,' ,T= ',E11.4,' ,R= ',E11.4,' Unitarity: ',E13.6)

 160  format(A7,'=',A22)

 165  format(18X,I5)
 166  format(D22.15)
 167  format(18X,A4)
 168  format(12X,A10)
 169  format(2X,A20)


 175  format(A7,'=            '     ,I10)
 176  format(A7,'=',              D22.15)
 177  format(A7,'=                  ',A4)
 178  format(A7,'=            '     ,A10)
 179  format(A7,'=  '               ,A20)

  END function check_bool

END MODULE AUXILIARY
