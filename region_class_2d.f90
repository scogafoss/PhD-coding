MODULE region_class_2d
    use line_class
    use region_class
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Filename: region_class_2d.f90                                            !!
  !!                                                                           !!
  !!  Dependant files: precision_set.f90, line_class.f90, material_class.f90,  !!
  !!    region_class.f90                                                       !!
  !!                                                                           !!
  !!  Author: Oliver Conway                             Start date: 10/05/2021 !!
  !!                                                                           !!
  !!  Purpose: Class to hold pointers to lines and materials in a 2D region    !!
  !!                                                                           !!
  !!  Revisions:                                                               !!
  !!                                                                           !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  IMPLICIT NONE
  ! Type definition
  TYPE,PUBLIC,extends(region) :: region_2d ! This will be the name we instantiate
  ! Instance variables.
  character(80) :: line_id_x ! line label in x direction(or r)
  character(80) :: line_id_y ! line label in y direction
  type(line),pointer :: line_x ! Points to line for region
  type(line),pointer :: line_y ! Points to line for region
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
  procedure,public :: get_top_boundary => get_top_boundary_fn ! Returns the top boundary condition
  procedure,public :: get_bottom_boundary => get_bottom_boundary_fn ! Returns the bottom boundary condition
  procedure,public :: get_d => get_d_fn ! Function to return diffusion coefficient
  procedure,public :: get_bottom_albedo => get_bottom_albedo_fn ! Returns albedo on bottom edge
  procedure,public :: get_top_albedo => get_top_albedo_fn ! Top edge albedo
  END TYPE region_2d
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
  private :: get_top_boundary_fn
  private :: get_bottom_boundary_fn
  private :: get_bottom_albedo_fn
  private :: get_top_albedo_fn
  ! Now add methods
  CONTAINS
  
    subroutine set_line_id_sub(this,id,dimension)
      !
      ! Subroutine to set material id
      !
      IMPLICIT NONE
      ! Declare calling arguments
      CLASS(region_2d),INTENT(OUT) :: this ! Region object
      integer,INTENT(IN) :: dimension
      character(len=*),INTENT(IN) :: id
      ! Save data
      if (dimension==1) then
        this%line_id_x = id
      elseif (dimension==2) then
        this%line_id_y = id
      else 
        stop 'Line must point in x or y direction'
      endif
    end subroutine set_line_id_sub
  
    subroutine associate_line_sub(this,lines)
      !
      ! Subroutine to point to line
      !
      IMPLICIT NONE
      ! Declare calling arguments
      CLASS(region_2d) :: this ! Region object
      type(line), target, intent(in) :: lines ! line object
      ! Save data
      if(lines%get_dimension()==1) then
        this%line_x => lines
      elseif(lines%get_dimension()==2) then
        this%line_y => lines
      else
        stop 'Dimension of 2D line must be 1 or 2'
      endif
    end subroutine associate_line_sub
  
    real(dp) FUNCTION get_start_fn(this,dimension)
      !
      ! Function to get start
      !
      IMPLICIT NONE
      ! Declare calling arguments
      CLASS(region_2d),INTENT(IN) :: this ! region object
      integer,INTENT(IN) :: dimension ! Dimension of the desired line
      if(dimension==1)then
        if(.not.associated(this%line_x)) stop 'Error no x line associated with region'
        get_start_fn = this%line_x%get_start()
      elseif(dimension==2)then
        if(.not.associated(this%line_y)) stop 'Error no y line associated with region'
        get_start_fn = this%line_y%get_start()
      else
        stop 'Dimension must be 1 or 2'
      endif
    END FUNCTION get_start_fn
  
    real(dp) FUNCTION get_length_fn(this,dimension)
      !
      ! Function to get length of region
      !
      IMPLICIT NONE
      ! Declare calling arguments
      CLASS(region_2d),INTENT(IN) :: this ! region object
      integer,INTENT(IN) :: dimension ! Dimension of the desired line
      if(dimension==1)then
        if(.not.associated(this%line_x)) stop 'Error no x line associated with region'
        get_length_fn = this%line_x%get_length()
      elseif(dimension==2)then
        if(.not.associated(this%line_y)) stop 'Error no y line associated with region'
        get_length_fn = this%line_y%get_length()
      else
        stop 'Dimension must be 1 or 2'
      endif
    END FUNCTION get_length_fn
  
    real(dp) FUNCTION get_delta_fn(this,dimension)
      !
      ! Function to get delta of region
      !
      IMPLICIT NONE
      ! Declare calling arguments
      CLASS(region_2d),INTENT(IN) :: this ! region object
      integer,INTENT(IN) :: dimension ! Dimension of the desired line
      if(dimension==1)then
        if(.not.associated(this%line_x)) stop 'Error no x line associated with region'
        get_delta_fn = this%line_x%get_delta()
      elseif(dimension==2)then
        if(.not.associated(this%line_y)) stop 'Error no y line associated with region'
        get_delta_fn = this%line_y%get_delta()
      else
        stop 'Dimension must be 1 or 2'
      endif
    END FUNCTION get_delta_fn
  
    integer FUNCTION get_steps_fn(this,dimension)
      !
      ! Function to get number of steps
      !
      IMPLICIT NONE
      ! Declare calling arguments
      CLASS(region_2d),INTENT(IN) :: this ! region object
      integer,INTENT(IN) :: dimension ! Dimension of the desired line
      if(dimension==1)then
        if(.not.associated(this%line_x)) stop 'Error no x line associated with region'
        get_steps_fn = this%line_x%get_steps()
      elseif(dimension==2)then
        if(.not.associated(this%line_y)) stop 'Error no y line associated with region'
        get_steps_fn = this%line_y%get_steps()
      else
        stop 'Dimension must be 1 or 2'
      endif
    END FUNCTION get_steps_fn
  
    integer FUNCTION get_geomtype_fn(this,dimension)
      !
      ! Function to return geomtype
      !
      IMPLICIT NONE
      ! Declare calling arguments
      CLASS(region_2d),INTENT(IN) :: this ! region object
      integer,INTENT(IN) :: dimension ! Dimension of the desired line
      if(dimension==1)then
        if(.not.associated(this%line_x)) stop 'Error no x line associated with region'
        get_geomtype_fn = this%line_x%get_geomtype()
      elseif(dimension==2)then
        if(.not.associated(this%line_y)) stop 'Error no y line associated with region'
        get_geomtype_fn = this%line_y%get_geomtype()
      else
        stop 'Dimension must be 1 or 2'
      endif
    END FUNCTION get_geomtype_fn
  
    character(80) FUNCTION get_line_id_fn(this,dimension)
      !
      ! Function to return id
      !
      IMPLICIT NONE
      ! Declare calling arguments
      CLASS(region_2d),INTENT(IN) :: this ! region object
      integer,INTENT(IN) :: dimension ! Dimension of the desired line
      if(dimension==1)then
        if(.not.associated(this%line_x)) stop 'Error no x line associated with region'
        get_line_id_fn = this%line_id_x
      elseif(dimension==2)then
        if(.not.associated(this%line_y)) stop 'Error no y line associated with region'
        get_line_id_fn = this%line_id_y
      else
        stop 'Dimension must be 1 or 2'
      endif
    END FUNCTION get_line_id_fn
  
    subroutine set_steps_sub(this,steps,dimension)
      !
      ! Function to return fission probability
      !
      IMPLICIT NONE
      ! Declare calling arguments
      CLASS(region_2d),INTENT(IN) :: this ! region object
      INTEGER,INTENT(IN) :: steps
      integer,INTENT(IN) :: dimension ! Dimension of the desired line
      if(dimension==1)then
        if(.not.associated(this%line_x)) stop 'Error no x line associated with region'
        call this%line_x%set_steps(steps)
      elseif(dimension==2)then
        if(.not.associated(this%line_y)) stop 'Error no y line associated with region'
        call this%line_y%set_steps(steps)
      else
        stop 'Dimension must be 1 or 2'
      endif
    end subroutine set_steps_sub

    character(80) function get_bottom_boundary_fn(this)
      !
      ! Function to return fission probability
      !
      IMPLICIT NONE
      ! Declare calling arguments
      CLASS(region_2d),INTENT(IN) :: this ! region object
      if(.not.associated(this%materials)) stop 'Error no material associated with region'
      get_bottom_boundary_fn=this%materials%get_bottom_boundary()
    end function get_bottom_boundary_fn

    character(80) function get_top_boundary_fn(this)
      !
      ! Function to return fission probability
      !
      IMPLICIT NONE
      ! Declare calling arguments
      CLASS(region_2d),INTENT(IN) :: this ! region object
      if(.not.associated(this%materials)) stop 'Error no material associated with region'
      get_top_boundary_fn=this%materials%get_top_boundary()
    end function get_top_boundary_fn

    real(dp) FUNCTION get_d_fn(this,group)
      !
      ! Function to get length of region
      !
      IMPLICIT NONE
      ! Declare calling arguments
      CLASS(region_2d),INTENT(IN) :: this ! region object
      INTEGER,INTENT(IN) :: group
      if(.not.associated(this%materials)) stop 'Error no material associated with region'
      get_d_fn = 1/(3*(this%get_absorption(group)+this%get_scatter(group)))
    END FUNCTION get_d_fn

    real(dp) FUNCTION get_bottom_albedo_fn(this)
      !
      ! Function to get bottom albedo of region
      !
      IMPLICIT NONE
      ! Declare calling arguments
      CLASS(region_2d),INTENT(IN) :: this ! region object
      if(.not.associated(this%materials)) stop 'Error no material associated with region'
      get_bottom_albedo_fn = this%materials%get_bottom_albedo()
    END FUNCTION get_bottom_albedo_fn

    real(dp) FUNCTION get_top_albedo_fn(this)
      !
      ! Function to get bottom albedo of region
      !
      IMPLICIT NONE
      ! Declare calling arguments
      CLASS(region_2d),INTENT(IN) :: this ! region object
      if(.not.associated(this%materials)) stop 'Error no material associated with region'
      get_top_albedo_fn = this%materials%get_top_albedo()
    END FUNCTION get_top_albedo_fn
    
  end module region_class_2d
  