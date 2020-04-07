program temp

  implicit none

  integer*8::bit
  integer*8::i,j
  integer*8::LUT

  open(unit=10,file='LUT.txt',action='write',status='replace')
  do i=0,2**20-1
    LUT=0
    do j=0,20
      LUT=LUT+iand(i,2**j)**3
    enddo
!    write(*,*)LUT
    write(10,*)LUT
  enddo
  close(10)

end program temp
