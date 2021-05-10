MODULE region_class
  use precision_set
  use material_class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Filename: region_class.f90                                             !!
!!                                                                           !!
!!  Dependant files: precision_set.f90                                       !!
!!                                                                           !!
!!  Author: Oliver Conway                             Start date: 11/02/2021 !!
!!                                                                           !!
!!  Purpose: Class to read in and store material information from file.      !!
!!                                                                           !!
!!  Revisions:                                                               !!
!!    10/05/2021: Moved some functions from region_class_1d to here          !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE
! Type definition
TYPE,PUBLIC :: region ! This will be the name we instantiate
! Instance variables.
character(80) :: material_id ! material label
type(material),pointer :: materials ! Points to material for region
CONTAINS
! Bound procedures
procedure,public :: set_material_id => set_material_id_sub ! Sets material ID
procedure,public :: associate_material => associate_material_sub ! Points to material
procedure,public :: get_material_id => get_material_id_fn ! Returns material ID
PROCEDURE,PUBLIC :: get_source_flux => get_source_flux_fn ! Returns source flux
PROCEDURE,PUBLIC :: get_absorption => get_absorption_fn ! Returns sigma_a
PROCEDURE,PUBLIC :: get_left_boundary => get_left_boundary_fn ! Returns left BC
PROCEDURE,PUBLIC :: get_right_boundary => get_right_boundary_fn ! Returns right BC
PROCEDURE,PUBLIC :: get_left_albedo => get_left_albedo_fn ! Returns left albedo
PROCEDURE,PUBLIC :: get_right_albedo => get_right_albedo_fn ! Returns right albedo
PROCEDURE,PUBLIC :: get_surface_source => get_surface_source_fn ! Returns surface flux
procedure,public :: get_fission => get_fission_fn ! Returns nu-sigma_f
procedure,public :: get_scatter => get_scatter_fn ! Returns scatter
procedure,public :: get_removal => get_removal_fn ! Returns removal cross sections
procedure,public :: get_probability => get_probability_fn ! Returns fission probability
END TYPE region
! Restrict access to the actual procedure names
PRIVATE :: set_material_id_sub
private :: associate_material_sub
private :: get_material_id_fn
private :: get_source_flux_fn
private :: get_absorption_fn
private :: get_left_boundary_fn
private :: get_right_boundary_fn
private :: get_left_albedo_fn
private :: get_right_albedo_fn
private :: get_surface_source_fn
private :: get_fission_fn
private :: get_scatter_fn
private :: get_removal_fn
private :: get_probability_fn
! Now add methods
CONTAINS

  subroutine set_material_id_sub(this,id)
    !
    ! Subroutine to set material id
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region) :: this ! Region object
    character(len=*),INTENT(IN) :: id
    ! Save data
    this%material_id = id
  end subroutine set_material_id_sub

  subroutine associate_material_sub(this,materials)
    !
    ! Subroutine to point to material
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region) :: this ! Region object
    type(material), target, intent(in) :: materials ! material object
    ! Save data
    this%materials => materials
  end subroutine associate_material_sub

  character(80) FUNCTION get_material_id_fn(this)
    !
    ! Function to return id
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region),INTENT(IN) :: this ! region object
    get_material_id_fn = this%material_id
  END FUNCTION get_material_id_fn

  FUNCTION get_source_flux_fn(this,index) result(value)
    !
    ! Function to return source flux
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region),INTENT(IN) :: this ! region object
    real(dp) :: value
    INTEGER :: index
    if(.not.associated(this%materials)) stop 'Error no material associated with region (source)'
    value = this%materials%get_source_flux(index)
  END FUNCTION get_source_flux_fn

  FUNCTION get_absorption_fn(this,index) result(value)
    !
    ! Function to return macroscopic absorption cross section
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region),INTENT(IN) :: this ! region object
    real(dp) :: value
    INTEGER :: index
    if(.not.associated(this%materials)) stop 'Error no material associated with region (absorption)'
    value = this%materials%get_absorption(index)
  END FUNCTION get_absorption_fn

  character(80) FUNCTION get_left_boundary_fn(this)
    !
    ! Function to return left boundary
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region),INTENT(IN) :: this ! region object
    if(.not.associated(this%materials)) stop 'Error no material associated with region (left boundary)'
    get_left_boundary_fn = this%materials%get_left_boundary()
  END FUNCTION get_left_boundary_fn

  character(80) FUNCTION get_right_boundary_fn(this)
    !
    ! Function to return right boundary
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region),INTENT(IN) :: this ! region object
    if(.not.associated(this%materials)) stop 'Error no material associated with region (right boundary)'
    get_right_boundary_fn = this%materials%get_right_boundary()
  END FUNCTION get_right_boundary_fn

  real(dp) FUNCTION get_left_albedo_fn(this)
    !
    ! Function to return left albedo
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region),INTENT(IN) :: this ! region object
    if(.not.associated(this%materials)) stop 'Error no material associated with region (albedo left)'
    get_left_albedo_fn = this%materials%get_left_albedo()
  END FUNCTION get_left_albedo_fn

  real(dp) FUNCTION get_right_albedo_fn(this)
    !
    ! Function to return right albedo
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region),INTENT(IN) :: this ! region object
    if(.not.associated(this%materials)) stop 'Error no material associated with region (albedo right)'
    get_right_albedo_fn = this%materials%get_right_albedo()
  END FUNCTION get_right_albedo_fn

  real(dp) FUNCTION get_surface_source_fn(this)
    !
    ! Function to return surface flux
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region),INTENT(IN) :: this ! region object
    if(.not.associated(this%materials)) stop 'Error no material associated with region (surface source)'
    get_surface_source_fn = this%materials%get_surface_source()
  END FUNCTION get_surface_source_fn

  FUNCTION get_fission_fn(this,index) result(value)
    !
    ! Function to return nu sigma f
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region),INTENT(IN) :: this ! region object
    real(dp) :: value
    INTEGER,INTENT(IN) :: index
    if(.not.associated(this%materials)) stop 'Error no material associated with region (fission)'
    value = this%materials%get_fission(index)
  END FUNCTION get_fission_fn

  FUNCTION get_scatter_fn(this,row,column) result(value)
    !
    ! Function to return scatter
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region),INTENT(IN) :: this ! region object
    real(dp) :: value
    integer,INTENT(IN) :: row
    integer,OPTIONAL :: column
    if(.not.associated(this%materials)) stop 'Error no material associated with region (scatter)'
    if (present(column)) then
      value = this%materials%get_scatter(row,column)
    else
      value = this%materials%get_scatter(row)
    end if
  END FUNCTION get_scatter_fn

  FUNCTION get_removal_fn(this,group) result(value)
    !
    ! Function to return removal
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region),INTENT(IN) :: this ! region object
    real(dp) :: value
    integer,INTENT(IN) :: group
    if(.not.associated(this%materials)) stop 'Error no material associated with region (removal)'
    value = this%materials%get_removal(group)
  END FUNCTION get_removal_fn

  FUNCTION get_probability_fn(this,index) result(value)
    !
    ! Function to return fission probability
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region),INTENT(IN) :: this ! region object
    real(dp):: value
    INTEGER,INTENT(IN) :: index
    if(.not.associated(this%materials)) stop 'Error no material associated with region (probability)'
    value = this%materials%get_probability(index)
  END FUNCTION get_probability_fn
  
end module region_class
