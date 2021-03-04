module read_gem_file
  !Reads in a gem file
  !Outputs x and phi
use precision_set
  implicit none
  contains
    subroutine read_gem(file_name,x,phi)
      character(len=*),intent(in) :: file_name
      real(dp), intent(out), allocatable, dimension(:) :: phi
      real(dp), intent(out), allocatable, dimension(:) :: x
      integer :: status ! Checks if end of file.
      character(len=30) :: line ! Used to test the line.
      integer :: x_length = 0 ! Records file length.
      integer :: x_pos,phi_pos = 0 ! Current position of x and phi
      open(1,file=file_name,iostat=status)
      ! First test for length of file.
      do
        read(1,*,iostat=status) line
        ! print*, 'status = (if not zero then it failed)',line
        if (status /= 0) exit ! exit if end of file (or fail).
        if (line == 'Position') then
          x_length = x_length + 1
        end if
      end do
      ! Go back to start of file.
      rewind (unit = 1)
      ! Now set allocate the data from the file.
      ! First set the arrays to correct size.
      allocate(phi(1:x_length)) ! x and phi length the same
      allocate(x(1:x_length))
      ! Now loop over file.
      open(1,file=file_name,iostat=status)
      do
        read(1,*,iostat=status) line
        if (status /= 0) exit ! exit if end of file (or fail).
        ! Look for the x values
        if (line == 'Position') then
          x_pos = x_pos + 1
          print*, 'x_pos=',x_pos
          backspace (unit = 1) ! goes back to start of line.
          read(1,*,iostat=status) line,line,x(x_pos) ! x value in the 3rd column
        ! Look for the phi values
        else if (line == 'Total') then
          phi_pos = phi_pos + 1
          backspace (unit = 1)
          read(1,*,iostat=status) line,phi(phi_pos) ! x value in the 2nd column
        end if
      end do
      close(1)
    end subroutine read_gem
end module read_gem_file
