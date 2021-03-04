MODULE region_class_1d
  use precision_set
  use line_class
  use material_class
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
PROCEDURE,PUBLIC :: get_left_boundary => get_left_boundary_fn ! Returns left BC
PROCEDURE,PUBLIC :: get_right_boundary => get_right_boundary_fn ! Returns right BC
PROCEDURE,PUBLIC :: get_left_albedo => get_left_albedo_fn ! Returns left albedo
PROCEDURE,PUBLIC :: get_right_albedo => get_right_albedo_fn ! Returns right albedo
PROCEDURE,PUBLIC :: get_surface_source => get_surface_source_fn ! Returns surface flux
procedure,public :: get_fission => get_fission_fn ! Returns nu-sigma_f
END TYPE region_1d
! Restrict access to the actual procedure names
private :: set_line_id_sub
private :: associate_line_sub
private :: get_start_fn
private :: get_length_fn
private :: get_steps_fn
private :: get_geomtype_fn
private :: get_line_id_fn
private :: get_left_boundary_fn
private :: get_right_boundary_fn
private :: get_left_albedo_fn
private :: get_right_albedo_fn
private :: get_surface_source_fn
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
  character(80) FUNCTION get_left_boundary_fn(this)
    !
    ! Function to return left boundary
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region_1d),INTENT(IN) :: this ! region object
    if(.not.associated(this%materials)) stop 'Error no material associated with region'
    get_left_boundary_fn = this%materials%get_left_boundary()
  END FUNCTION get_left_boundary_fn
  character(80) FUNCTION get_right_boundary_fn(this)
    !
    ! Function to return right boundary
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region_1d),INTENT(IN) :: this ! region object
    if(.not.associated(this%materials)) stop 'Error no material associated with region'
    get_right_boundary_fn = this%materials%get_right_boundary()
  END FUNCTION get_right_boundary_fn
  real(dp) FUNCTION get_left_albedo_fn(this)
    !
    ! Function to return left albedo
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region_1d),INTENT(IN) :: this ! region object
    if(.not.associated(this%materials)) stop 'Error no material associated with region'
    get_left_albedo_fn = this%materials%get_left_albedo()
  END FUNCTION get_left_albedo_fn
  real(dp) FUNCTION get_right_albedo_fn(this)
    !
    ! Function to return right albedo
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region_1d),INTENT(IN) :: this ! region object
    if(.not.associated(this%materials)) stop 'Error no material associated with region'
    get_right_albedo_fn = this%materials%get_right_albedo()
  END FUNCTION get_right_albedo_fn
  real(dp) FUNCTION get_surface_source_fn(this)
    !
    ! Function to return surface flux
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region_1d),INTENT(IN) :: this ! region object
    if(.not.associated(this%materials)) stop 'Error no material associated with region'
    get_surface_source_fn = this%materials%get_surface_source()
  END FUNCTION get_surface_source_fn
  real(dp) FUNCTION get_fission_fn(this)
    !
    ! Function to return surface flux
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region_1d),INTENT(IN) :: this ! region object
    if(.not.associated(this%materials)) stop 'Error no material associated with region'
    get_fission_fn = this%materials%get_fission()
  END FUNCTION get_fission_fn
end module region_class_1d
