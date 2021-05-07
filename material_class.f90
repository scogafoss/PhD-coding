MODULE material_class
  use precision_set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Filename: material_class.f90                                             !!
!!                                                                           !!
!!  Dependant files: precision_set.f90                                       !!
!!                                                                           !!
!!  Author: Oliver Conway                             Start date: 12/01/2021 !!
!!                                                                           !!
!!  Purpose: Class to read in and store material information from file.      !!
!!                                                                           !!
!!  Revisions:                                                               !!
!!    05/02/2021: Updated how material was read in. Added material ID        !!
!!    16/02/2021: Added material macroscopic fission cross section           !!
!!    08/03/2021: Added macroscopic scatter cross section array, changed     !!
!!      source_flux, absorption, fission probability and fission to arrays.  !!
!!    07/05/2021: Added top and bottom boundary conditions and albedo        !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE
! Type definition
TYPE,PUBLIC :: material ! This will be the name we instantiate
! Instance variables.
PRIVATE
character(80) :: left_boundary ! Boundary condition on left side of region
character(80) :: right_boundary ! Boundary condition on right side of region
character(80) :: id = 'default'! Material ID
real(dp) :: left_albedo ! Left albedo value
real(dp) :: right_albedo ! Right albedo value
real(dp) :: surface_source ! Surface flux
real(dp),allocatable,dimension(:) :: source_flux ! Source flux G array
real(dp),allocatable,dimension(:) :: absorption ! Macroscopic absorption cross section G array
real(dp),allocatable,dimension(:) :: fission ! Macroscopic fission cross section G array
real(dp),allocatable,dimension(:) :: fission_prob ! Macroscopic fission cross section G array
real(dp),allocatable,dimension(:,:) :: scatter ! Macroscopic scatter cross section GxG matrix
character(80) :: top_boundary ! Boundary conditon on the top edge of mesh
character(80) :: bottom_boundary ! Boundary at the bottom edge
real(dp) :: top_albedo ! Top albedo value
real(dp) :: bottom_albedo ! Bottom albedo value

CONTAINS
! Bound procedures
PROCEDURE,PUBLIC :: set_variables => set_variables_sub ! Allows user to input desired variables
PROCEDURE,PUBLIC :: read_variables => read_variables_sub ! Reads in variables from file
PROCEDURE,PUBLIC :: get_left_boundary => get_left_boundary_fn ! Returns left BC
PROCEDURE,PUBLIC :: get_right_boundary => get_right_boundary_fn ! Returns right BC
PROCEDURE,PUBLIC :: get_left_albedo => get_left_albedo_fn ! Returns left albedo
PROCEDURE,PUBLIC :: get_right_albedo => get_right_albedo_fn ! Returns right albedo
PROCEDURE,PUBLIC :: get_surface_source => get_surface_source_fn ! Returns surface flux
PROCEDURE,PUBLIC :: get_source_flux => get_source_flux_fn ! Returns source flux
PROCEDURE,PUBLIC :: get_absorption => get_absorption_fn ! Returns sigma_a
procedure,public :: set_id => set_id_sub ! Allows name to be set
procedure,public :: get_id => get_id_fn ! Returns ID
procedure,public :: get_fission => get_fission_fn ! Returns Nu-Sigma_f
procedure,public :: get_removal => get_removal_fn ! Returns (and calculates) the removal cross section
procedure,public :: get_scatter => get_scatter_fn ! Returns scatter matrix
procedure,PUBLIC :: get_probability => get_probability_fn ! Returns fission probability
PROCEDURE,PUBLIC :: get_top_boundary => get_top_boundary_fn ! Returns top BC
PROCEDURE,PUBLIC :: get_bottom_boundary => get_bottom_boundary_fn ! Returns bottom BC
PROCEDURE,PUBLIC :: get_top_albedo => get_top_albedo_fn ! Returns top albedo
PROCEDURE,PUBLIC :: get_bottom_albedo => get_bottom_albedo_fn ! Returns bottom albedo
END TYPE material
! Restrict access to the actual procedure names
PRIVATE :: set_variables_sub
private :: read_variables_sub
private :: get_left_boundary_fn
private :: get_right_boundary_fn
private :: get_left_albedo_fn
private :: get_right_albedo_fn
private :: get_surface_source_fn
private :: get_source_flux_fn
private :: get_absorption_fn
private :: set_id_sub
private :: get_id_fn
private :: get_fission_fn
private :: get_removal_fn
private :: get_scatter_fn
private :: get_probability_fn
private :: get_top_boundary_fn
private :: get_bottom_boundary_fn
private :: get_top_albedo_fn
private :: get_bottom_albedo_fn
! Now add methods
CONTAINS

SUBROUTINE set_variables_sub(this, left_boundary, right_boundary, left_albedo, &
   right_albedo, surface_source, source_flux, absorption, fission, scatter,&
   top_boundary, bottom_boundary)
  !
  ! Subroutine to set the variables
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material) :: this ! Line object
  character(len=*),INTENT(IN) :: left_boundary, right_boundary
  character(len=*),INTENT(IN) :: top_boundary, bottom_boundary
  real(dp),INTENT(IN) :: left_albedo, right_albedo, surface_source
  real(dp),INTENT(IN), dimension(:) :: source_flux, absorption, fission
  real(dp),INTENT(IN), dimension(:,:) :: scatter
  ! Save data
  this%left_boundary = left_boundary
  this%right_boundary = right_boundary
  this%left_albedo = left_albedo
  this%right_albedo = right_albedo
  this%surface_source = surface_source
  this%source_flux = source_flux
  this%absorption = absorption
  this%fission = fission
  this%scatter = scatter
  this%top_boundary = top_boundary
  this%bottom_boundary = bottom_boundary
END SUBROUTINE set_variables_sub

SUBROUTINE read_variables_sub(this, filename)
  !
  ! Subroutine to read in line data from input deck
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material) :: this ! Matrix object
  character(len=*),INTENT(IN) :: filename
  integer :: status ! Checks if end of file.
  integer :: lineskip ! skips desired number of lines.
  character(80) :: line ! Used to test the line.
  real(dp) :: left_albedo, right_albedo, bottom_albedo, top_albedo
  real(dp), allocatable, dimension(:) :: absorption, source_flux, fission, probability
  real(dp), allocatable, dimension(:,:) :: scatter
  real(dp) :: surface_source
  character(80) :: left_boundary, right_boundary, bottom_boundary, top_boundary
  integer groups
  open(10,file=filename,iostat=status)
  ! Read through file.
  do
    read(10,'(A)',iostat=status) line ! (10,'(A)') <-- the 10 indicates write to file 10, the '(A)' indicates read the full line as a string
    ! print *, line
    if (status /= 0) exit ! exit if end of file (or fail).
    ! Check the number of groups
    if (line == 'Groups') then
      read(10,*,iostat=status) line,groups
      ! Allocate the variables to appropriate size.
      allocate(absorption(1:groups))
      allocate(source_flux(1:groups))
      allocate(fission(1:groups))
      allocate(scatter(1:groups,1:groups))
      allocate(probability(1:groups))
    end if
    if (line == "Fission probability") then
      read(10,*,iostat=status) probability
      this%fission_prob = probability
    end if
    ! Check if the material has a defined name
    if (this%id == 'default') then ! If no defined name
      if (line == 'ref (mxx), sigma_a, source flux, nu sigma_f') then
        ! Loop for each of the groups
        do lineskip = 1,groups
          read(10,*,iostat=status) line,absorption(lineskip),source_flux(lineskip),fission(lineskip)
        end do
      ! if (line == 'macroscopic absorption cross section (sigma_a), source flux (S)') then
      !   read(10,*,iostat=status) absorption,source_flux
      !   this%absorption = absorption
      !   this%source_flux = source_flux
      else if (line == '------------SCATTER CROSS SECTIONS------------------') then
        do lineskip = 1,1 ! Skips line.
          read(10,*)
        end do
        do lineskip = 1,groups
          ! Read in columns for "lineskip"th row of the matrix
          read(10,*,iostat=status) scatter(lineskip,:)
        end do
      else if (line == 'Boundaries-') then
        read(10,*,iostat=status) line,left_boundary,right_boundary
        read(10,*,iostat=status) line,left_albedo,right_albedo
        read(10,*,iostat=status) line,bottom_boundary,top_boundary
        read(10,*,iostat=status) line,bottom_albedo,top_albedo
        do lineskip = 1,2 ! Skips two lines.
          read(10,*)
        end do
        read(10,*,iostat=status) line,surface_source
        this%left_boundary = left_boundary
        this%right_boundary = right_boundary
        this%left_albedo = left_albedo
        this%right_albedo = right_albedo
        this%surface_source = surface_source
      end if
    ELSE ! If there is a defined name
      if (line == 'ref (mxx), sigma_a, source flux, nu sigma_f') then
        do ! Loop until the name of the line is found, or the end of the section is reached.
          read(10,'(A)',iostat=status) line
          if (index(line,trim(this%id))/=0) then ! If the line contains the name we are looking for
            backspace (unit = 10) ! Goes back to start of line
            do lineskip = 1,groups
              read(10,*,iostat=status) line,absorption(lineskip),source_flux(lineskip),fission(lineskip)
              ! If any of the lines has the wrong name, implies wrong number of data points
              if(line/=trim(this%id)) stop 'Number of groups does not match data in input deck'
            end do
            exit
          end if
          if (index(line,'---')/=0) stop 'No match for material name in input deck'
        end do
      ! if (line == 'macroscopic absorption cross section (sigma_a), source flux (S)') then
      !   read(10,*,iostat=status) absorption,source_flux
      !   this%absorption = absorption
      !   this%source_flux = source_flux
      else if (line == '------------SCATTER CROSS SECTIONS------------------') then
        do ! Loop until the name of the line is found, or the end of the section is reached.
          read(10,'(A)',iostat=status) line
          if (index(line,trim(this%id))/=0) then ! If the line contains the name we are looking for
            do lineskip = 1,groups
              ! Read in columns for "lineskip"th row of the matrix
              read(10,*,iostat=status) scatter(lineskip,:)
            end do
            exit
          end if
          if (index(line,'---')/=0) stop 'No match for material name in input deck'
        end do
      else if (line == 'Boundaries-') then
        read(10,*,iostat=status) line,left_boundary,right_boundary
        read(10,*,iostat=status) line,left_albedo,right_albedo
        do lineskip = 1,2 ! Skips two lines.
          read(10,*)
        end do
        read(10,*,iostat=status) line,surface_source
        this%left_boundary = left_boundary
        this%right_boundary = right_boundary
        this%left_albedo = left_albedo
        this%right_albedo = right_albedo
        this%surface_source = surface_source
      end if
    end if
  end do
  this%absorption = absorption
  this%source_flux = source_flux
  this%fission = fission
  this%scatter = scatter
  close(10)
END SUBROUTINE read_variables_sub

character(80) FUNCTION get_left_boundary_fn(this)
  !
  ! Function to return left BC
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Line object
  get_left_boundary_fn = this%left_boundary
END FUNCTION get_left_boundary_fn

character(80) FUNCTION get_right_boundary_fn(this)
  !
  ! Function to return right BC
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Line object
  get_right_boundary_fn = this%right_boundary
END FUNCTION get_right_boundary_fn

real(dp) FUNCTION get_left_albedo_fn(this)
  !
  ! Function to return right BC
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Line object
  get_left_albedo_fn = this%left_albedo
END FUNCTION get_left_albedo_fn

real(dp) FUNCTION get_right_albedo_fn(this)
  !
  ! Function to return right BC
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Line object
  get_right_albedo_fn = this%right_albedo
END FUNCTION get_right_albedo_fn

real(dp) FUNCTION get_surface_source_fn(this)
  !
  ! Function to return right BC
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Line object
  get_surface_source_fn = this%surface_source
END FUNCTION get_surface_source_fn

FUNCTION get_source_flux_fn(this,index) result(value)
  !
  ! Function to return right BC
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Line object
  real(dp) :: value
  INTEGER,INTENT(IN) :: index
  value = this%source_flux(index)
END FUNCTION get_source_flux_fn

FUNCTION get_absorption_fn(this,index) result(value)
  !
  ! Function to return right BC
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Line object
  real(dp) :: value
  INTEGER,INTENT(IN) :: index
  value = this%absorption(index)
END FUNCTION get_absorption_fn

subroutine set_id_sub (this, id)
  !
  ! Subroutine to define the name
  !
  implicit none
  ! Declare calling arguments
  class(material) :: this
  character(len=*),INTENT(IN) :: id
  this%id = id
end subroutine set_id_sub

character(80) FUNCTION get_id_fn(this)
  !
  ! Function to return right BC
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Line object
  get_id_fn = this%id
END FUNCTION get_id_fn

FUNCTION get_fission_fn(this,index) result(value)
  !
  ! Function to return nu-sigma_f
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Line object
  INTEGER,INTENT(IN) :: index
  real(dp) :: value
  value = this%fission(index)
END FUNCTION get_fission_fn

function get_removal_fn(this,group) result(value)
  !
  ! Function to return removal cross section
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Line object
  integer,intent(in) :: group
  real(dp) :: value
  integer :: i
  value = 0
  do i = 1, size(this%fission)
    ! Ignore in-group scattering
    if (i/=group) then
      ! Account for all scatterings out of the group
      value = value + this%scatter(group,i)
    end if
  end do
  value = value + this%absorption(group)
end function get_removal_fn

FUNCTION get_scatter_fn(this,row,column) result(value)
  !
  ! Function to return scatter matrix
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Line object
  real(dp) :: value
  integer,INTENT(IN) :: row
  integer,intent(in),optional :: column
  INTEGER :: i
  if (present(column)) then
    value = this%scatter(row,column)
  else
    value = 0
    do i = 1,size(this%absorption)
      value=value+this%scatter(row,i)
    end do
  end if
END FUNCTION get_scatter_fn

FUNCTION get_probability_fn(this,index) result(value)
  !
  ! Function to return fission probabilities of each group
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Line object
  INTEGER,INTENT(IN) :: index
  real(dp) :: value
  value = this%fission_prob(index)
END FUNCTION get_probability_fn

character(80) FUNCTION get_top_boundary_fn(this)
  !
  ! Function to return left BC
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Material object
  get_top_boundary_fn = this%top_boundary
END FUNCTION get_top_boundary_fn

character(80) FUNCTION get_bottom_boundary_fn(this)
  !
  ! Function to return left BC
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Line object
  get_bottom_boundary_fn = this%bottom_boundary
END FUNCTION get_bottom_boundary_fn

real(dp) FUNCTION get_top_albedo_fn(this)
  !
  ! Function to return right BC
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Line object
  get_top_albedo_fn = this%top_albedo
END FUNCTION get_top_albedo_fn

real(dp) FUNCTION get_bottom_albedo_fn(this)
  !
  ! Function to return right BC
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(material),INTENT(IN) :: this ! Line object
  get_bottom_albedo_fn = this%bottom_albedo
END FUNCTION get_bottom_albedo_fn

end module material_class
