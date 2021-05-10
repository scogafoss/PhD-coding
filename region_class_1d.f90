MODULE region_class_1d
  use line_class
  use region_class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Filename: region_class_1d.f90                                            !!
!!                                                                           !!
!!  Dependant files: precision_set.f90, line_class.f90, material_class.f90,  !!
!!    region_class.f90                                                       !!
!!                                                                           !!
!!  Author: Oliver Conway                             Start date: 11/02/2021 !!
!!                                                                           !!
!!  Purpose: Class to hold pointers to lines and materials in a region       !!
!!                                                                           !!
!!  Revisions:                                                               !!
!!    16/02/2021: Added function to return nu-sigma_f                        !!
!!    10/03/2021: Added functions to return prob, removal, delta and scatter !!
!!    10/05/2021: Reshuffled functions which belonged in region_class.f90    !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE
! Type definition
TYPE,PUBLIC,extends(region) :: region_1d ! This will be the name we instantiate
! Instance variables.
character(80) :: line_id ! line label
type(line),pointer :: lines ! Points to line for region
CONTAINS
! Bound procedures
procedure,public :: set_line_id => set_line_id_sub ! Sets line ID
procedure,public :: associate_line => associate_line_sub ! Points to line
procedure,public :: get_start => get_start_fn ! Returns the start x
procedure,public :: get_length => get_length_fn ! Returns length
procedure,public :: get_steps => get_steps_fn ! Returns steps
procedure,public :: get_geomtype => get_geomtype_fn ! Returns geomtype
procedure,public :: get_line_id => get_line_id_fn ! Returns line ID
procedure,public :: get_delta => get_delta_fn ! Returns delta
procedure,public :: set_steps => set_steps_sub ! Sets the number of steps in line class - used for periodic BC
END TYPE region_1d
! Restrict access to the actual procedure names
private :: set_line_id_sub
private :: associate_line_sub
private :: get_start_fn
private :: get_length_fn
private :: get_steps_fn
private :: get_geomtype_fn
private :: get_line_id_fn
private :: get_delta_fn
private :: set_steps_sub
! Now add methods
CONTAINS

  subroutine set_line_id_sub(this,id)
    !
    ! Subroutine to set material id
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region_1d) :: this ! Region object
    character(len=*),INTENT(IN) :: id
    ! Save data
    this%line_id = id
  end subroutine set_line_id_sub

  subroutine associate_line_sub(this,lines)
    !
    ! Subroutine to point to line
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region_1d) :: this ! Region object
    type(line), target, intent(in) :: lines ! line object
    ! Save data
    this%lines => lines
  end subroutine associate_line_sub

  real(dp) FUNCTION get_start_fn(this)
    !
    ! Function to get start
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region_1d),INTENT(IN) :: this ! region object
    if(.not.associated(this%lines)) stop 'Error no line associated with region'
    get_start_fn = this%lines%get_start()
  END FUNCTION get_start_fn

  real(dp) FUNCTION get_length_fn(this)
    !
    ! Function to get length of region
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region_1d),INTENT(IN) :: this ! region object
    if(.not.associated(this%lines)) stop 'Error no line associated with region'
    get_length_fn = this%lines%get_length()
  END FUNCTION get_length_fn

  real(dp) FUNCTION get_delta_fn(this)
    !
    ! Function to get delta of region
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region_1d),INTENT(IN) :: this ! region object
    if(.not.associated(this%lines)) stop 'Error no line associated with region'
    get_delta_fn = this%lines%get_delta()
  END FUNCTION get_delta_fn

  integer FUNCTION get_steps_fn(this)
    !
    ! Function to get number of steps
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region_1d),INTENT(IN) :: this ! region object
    if(.not.associated(this%lines)) stop 'Error no line associated with region'
    get_steps_fn = this%lines%get_steps()
  END FUNCTION get_steps_fn

  integer FUNCTION get_geomtype_fn(this)
    !
    ! Function to return geomtype
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region_1d),INTENT(IN) :: this ! region object
    if(.not.associated(this%lines)) stop 'Error no line associated with region'
    get_geomtype_fn = this%lines%get_geomtype()
  END FUNCTION get_geomtype_fn

  character(80) FUNCTION get_line_id_fn(this)
    !
    ! Function to return id
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region_1d),INTENT(IN) :: this ! region object
    get_line_id_fn = this%line_id
  END FUNCTION get_line_id_fn

  subroutine set_steps_sub(this,steps)
    !
    ! Function to return fission probability
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region_1d),INTENT(IN) :: this ! region object
    INTEGER,INTENT(IN) :: steps
    if(.not.associated(this%lines)) stop 'Error no line associated with region (set_steps)'
    call this%lines%set_steps(steps)
  end subroutine set_steps_sub
  
end module region_class_1d
