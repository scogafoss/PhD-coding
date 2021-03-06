MODULE line_class
  use precision_set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Filename: line_class.f90                                                 !!
!!                                                                           !!
!!  Dependant files: precision_set.f90                                       !!
!!                                                                           !!
!!  Author: Oliver Conway                             Start date: 12/01/2021 !!
!!                                                                           !!
!!  Purpose: Class to read in and store dimensional information from file.   !!
!!                                                                           !!
!!  Revisions:                                                               !!
!!    05/02/2021: Updated how line was read in. Added line ID                !!
!!    30/04/2021: Added delta to variables                                   !!
!!    07/05/2021: Added dimension to store what direction line points        !!  
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE
! Type definition
TYPE,PUBLIC :: line ! This will be the name we instantiate
! Instance variables.
PRIVATE
character(80) :: id = 'default' ! Line id
real(dp) :: start ! Start of region
real(dp) :: length ! Length of region
integer :: steps ! Number of nodes
integer :: geomtype ! Type of geometry (0:slab, 1:cylindrical, 2:spherical)
real(dp) :: delta ! Stores the step size
integer :: dimension ! 1 = x/r, 2 = y, 3 = z

CONTAINS
! Bound procedures
PROCEDURE,PUBLIC :: set_variables => set_variables_sub ! Allows user to input desired variables
PROCEDURE,PUBLIC :: read_variables => read_variables_sub ! Reads in variables from file
PROCEDURE,PUBLIC :: get_start => get_start_fn ! Returns start of region
PROCEDURE,PUBLIC :: get_length => get_length_fn ! Returns length of region
PROCEDURE,PUBLIC :: get_steps => get_steps_fn ! Returns number of nodes
PROCEDURE,PUBLIC :: get_geomtype => get_geomtype_fn ! Returns geomtype
procedure,public ::set_id => set_id_sub ! Allows name to be set
procedure,public ::get_id => get_id_fn ! Returns ID
procedure,public :: get_delta => get_delta_fn ! Returns delta
procedure,public :: set_steps => set_steps_sub ! Directly sets the steps - used for periodic BC
procedure,public :: get_dimension => get_dimension_fn ! Returns dimension of the line
END TYPE line
! Restrict access to the actual procedure names
PRIVATE :: set_variables_sub
private :: read_variables_sub
private :: get_start_fn
private :: get_length_fn
private :: get_steps_fn
private :: get_geomtype_fn
private :: set_id_sub
private :: get_id_fn
private :: get_delta_fn
private :: set_steps_sub
private :: get_dimension_fn
! Now add methods
CONTAINS

SUBROUTINE set_variables_sub(this, start, length, steps, geomtype, dimension)
  !
  ! Subroutine to set the variables
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(line) :: this ! Line object
  real(dp),INTENT(IN) :: start
  real(dp),INTENT(IN) :: length
  integer,INTENT(IN) :: steps
  integer,INTENT(IN) :: geomtype
  integer, INTENT(IN) :: dimension
  ! Save data
  this%start = start
  this%length = length
  this%steps = steps
  this%geomtype = geomtype
  this%delta = this%length / this%steps
  this%dimension = dimension
END SUBROUTINE set_variables_sub

SUBROUTINE read_variables_sub(this, filename)
  !
  ! Subroutine to read in line data from input deck
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(line) :: this ! Matrix object
  character(len=*),INTENT(IN) :: filename
  integer :: status ! Checks if end of file.
  integer :: lineskip ! skips desired number of lines.
  character(80) :: line ! Used to test the line.
  real(dp) :: start, length
  integer :: steps, geomtype, dimension
  open(10,file=filename,iostat=status)
  ! Read through file.
  do
  	read(10,'(A)',iostat=status) line ! (10,'(A)') <-- the 10 indicates write to file 10, the '(A)' indicates read the full line as a string
  	if (status /= 0) exit ! exit if end of file (or fail).
    ! Check if the line has a defined name
    if (this%id == 'default') then ! If no defined name
    	if (line == 'start, length, steps, ref (lxx), dim (x/r,y,z)') then
    		do lineskip = 1,3 ! Skips three lines.
    			read(10,*)
    	  end do
    		read(10,*,iostat=status) start,length,steps,line,line
        this%start = start
        this%length = length
        this%steps = steps
      else if (line == 'Geometry - 0: slab, 1: cylindrical, 2: spherical') then
        read(10,*,iostat=status) line,geomtype
        this%geomtype = geomtype
    	end if
    else ! If there is a defined ID
      if (line == 'start, length, steps, ref (lxx), dim (x/r,y,z)') then
        do ! Loop until the name of the line is found, or the end of the section is reached.
          read(10,'(A)',iostat=status) line
          if (index(line,trim(this%id))/=0) then ! If the line contains the name we are looking for
            backspace (unit = 10) ! Goes back to start of line
            read(10,*,iostat=status) start,length,steps,line,dimension
            this%start = start
            this%length = length
            this%steps = steps
            this%delta = length / steps
            this%dimension = dimension
            exit
          end if
          if (index(line,'---')/=0) stop 'No match for line name in input deck'
        end do
      else if (line == 'Geometry - 0: slab, 1: cylindrical, 2: spherical') then
        read(10,*,iostat=status) line,geomtype
        this%geomtype = geomtype
      end if
    end if
  end do
  close(10)
  ! this%delta = this%length / this%steps
END SUBROUTINE read_variables_sub

real(dp) FUNCTION get_start_fn(this)
  !
  ! Function to return start of region
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(line),INTENT(IN) :: this ! Line object
  get_start_fn = this%start
END FUNCTION get_start_fn

real(dp) FUNCTION get_length_fn(this)
  !
  ! Function to return start of region
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(line),INTENT(IN) :: this ! Line object
  get_length_fn = this%length
END FUNCTION get_length_fn

integer FUNCTION get_steps_fn(this)
  !
  ! Function to return start of region
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(line),INTENT(IN) :: this ! Line object
  get_steps_fn = this%steps
END FUNCTION get_steps_fn

integer FUNCTION get_geomtype_fn(this)
  !
  ! Function to return start of region
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(line),INTENT(IN) :: this ! Line object
  get_geomtype_fn = this%geomtype
END FUNCTION get_geomtype_fn

subroutine set_id_sub (this, id)
  !
  ! Subroutine to define the name
  !
  class(line) :: this
  character(len=*),INTENT(IN) :: id
  this%id = id
end subroutine set_id_sub

character(80) FUNCTION get_id_fn(this)
  !
  ! Function to return right BC
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(line),INTENT(IN) :: this ! Line object
  get_id_fn = this%id
END FUNCTION get_id_fn

subroutine set_steps_sub(this,steps)
  !
  ! Subroutine to set the steps
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(line) :: this ! Line object
  integer,INTENT(IN) :: steps
  ! Save datum
  this%steps = steps !Note this preserves the original delta as this is just for periodic
end subroutine set_steps_sub

real(dp) FUNCTION get_delta_fn(this)
  !
  ! Function to return right BC
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(line),INTENT(IN) :: this ! Line object
  get_delta_fn = this%delta
END FUNCTION get_delta_fn

integer FUNCTION get_dimension_fn(this)
  !
  ! Function to return right BC
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(line),INTENT(IN) :: this ! Line object
  get_dimension_fn = this%dimension
END FUNCTION get_dimension_fn

END MODULE line_class
