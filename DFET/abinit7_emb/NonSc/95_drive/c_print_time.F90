subroutine c_print_time(io)
  implicit none
  integer,intent(in) :: io
  integer,dimension(8) :: dateTime
  character(len=500) :: message

  call date_and_time(VALUES=dateTime)    !today
  write(message, '(i2.2,A,i2.2,A,i2.2,A,i3.3)')dateTime(5), ":", dateTime(6), ":", dateTime(7), ".", dateTime(8)
  call wrtout(io,message,"COLL")
end subroutine c_print_time

