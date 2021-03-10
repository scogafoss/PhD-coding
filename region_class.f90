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
END TYPE region
! Restrict access to the actual procedure names
PRIVATE :: set_material_id_sub
private :: associate_material_sub
private :: get_material_id_fn
private :: get_source_flux_fn
private :: get_absorption_fn
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
FUNCTION get_source_flux_fn(this) result(value)
    !
    ! Function to return source flux
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region),INTENT(IN) :: this ! region object
    real(dp),allocatable,DIMENSION(:) :: value
    if(.not.associated(this%materials)) stop 'Error no material associated with region'
    value = this%materials%get_source_flux()
  END FUNCTION get_source_flux_fn
FUNCTION get_absorption_fn(this) result(value)
    !
    ! Function to return macroscopic absorption cross section
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(region),INTENT(IN) :: this ! region object
    real(dp),allocatable,DIMENSION(:) :: value
    if(.not.associated(this%materials)) stop 'Error no material associated with region'
    value = this%materials%get_absorption()
  END FUNCTION get_absorption_fn
end module region_class
