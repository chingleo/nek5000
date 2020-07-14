program main
    implicit none

    integer :: i,k
    character(len=40) :: filename
    character(len=40) :: filenum
    character(len=40) :: num

    do k = 1,9,1
    write(num,*) k
    open(20, file = "a"//trim(adjustl(num))//"_file.list", status = 'unknown')
    write(20,'(a)') '2'
    do i = 2,3,1
    write(filenum,'(I5.5)') i
    filename = '../avg2/a'//trim(adjustl(num))//'_disk3d0.f'//trim(adjustl(filenum))
    write(20,'(a)') filename
    enddo
    close(20)
    enddo

    end program main
