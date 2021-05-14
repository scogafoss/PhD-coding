MODULE mesh_class
  use line_class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Filename: mesh_class.f90                                                 !!
!!                                                                           !!
!!  Dependant files: line_class.f90                                          !!
!!                                                                           !!
!!  Author: Oliver Conway                             Start date: 12/01/2021 !!
!!                                                                           !!
!!  Purpose: Class to store x and y values an boundary positions.            !!
!!                                                                           !!
!!  Revisions:                                                               !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE
! Type definition
TYPE,PUBLIC :: mesh ! This will be the name we instantiate
! Instance variables.
PRIVATE
real(dp),allocatable,dimension(:) :: x_coordinates
real(dp),allocatable,dimension(:) :: y_coordinates
integer,allocatable,dimension(:) :: x_boundaries
integer,allocatable,dimension(:) :: y_boundaries

CONTAINS
! Bound procedures
PROCEDURE,PUBLIC :: fill_coordinates => fill_coordinates_sub ! Fills variables using unput line class array
procedure,public :: allocate_coordinates => allocate_coordinates_sub ! Used to allocate size of arrays
procedure,public :: get_x => get_x_fn ! Return x
procedure,public :: get_y => get_y_fn ! Return y
procedure,public :: get_x_boundary => get_x_boundary_fn ! Return x boundary
procedure,public :: get_y_boundary => get_y_boundary_fn ! Return y boundary
procedure,public :: get_x_size => get_x_size_fn ! Return number x nodes
procedure,public :: get_y_size => get_y_size_fn ! Return number y nodes
procedure,public :: number_regions_x => number_regions_x_fn ! Return number x regions
procedure,public :: number_regions_y => number_regions_y_fn ! Return number y regions
procedure,public :: at_boundary_l => at_boundary_l_fn ! Boolean to test if boundary to the left
procedure,public :: at_boundary_r => at_boundary_r_fn
procedure,public :: at_boundary_t => at_boundary_t_fn
procedure,public :: at_boundary_b => at_boundary_b_fn
procedure,public :: at_edge_l => at_edge_l_fn ! Boolean to test if edge of mesh to the left
procedure,public :: at_edge_r => at_edge_r_fn
procedure,public :: at_edge_t => at_edge_t_fn
procedure,public :: at_edge_b => at_edge_b_fn
procedure,public :: r => r_fn ! Returns the index of region which corresponds to desired i,j box
END TYPE mesh
! Restrict access to the actual procedure names
PRIVATE :: fill_coordinates_sub
private :: allocate_coordinates_sub
private :: get_x_fn
private :: get_y_fn
private :: get_x_boundary_fn
private :: get_y_boundary_fn
private :: get_x_size_fn
private :: get_y_size_fn
private :: number_regions_x_fn
private :: number_regions_y_fn
private :: at_boundary_l_fn
private :: at_boundary_r_fn
private :: at_boundary_t_fn
private :: at_boundary_b_fn
private :: at_edge_l_fn
private :: at_edge_r_fn
private :: at_edge_t_fn
private :: at_edge_b_fn
private :: r_fn
! Now add methods
CONTAINS

SUBROUTINE fill_coordinates_sub(this, lines)
  !
  ! Subroutine to set the class variables
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(mesh),intent(out) :: this ! Mesh object
  type(line),dimension(:),INTENT(IN) :: lines
  integer :: nodex,nodey,i,j
  call this%allocate_coordinates(lines)
  nodex = 0
  nodey = 0
  do i = 1, size(lines)
    if (lines(i)%get_dimension() == 1) then
      do j = 1,lines(i)%get_steps()
        nodex=nodex+1
        this%x_coordinates(nodex) = lines(i)%get_start()+(j-0.5)*lines(i)%get_delta()
      enddo
    elseif (lines(i)%get_dimension() == 2) then
      do j =1, lines(i)%get_steps()
        nodey = nodey +1
        this%y_coordinates(nodey) = lines(i)%get_start()+(j-0.5)*lines(i)%get_delta()
      enddo
    else
      Stop 'Need dimension 1 or 2'
    endif
  enddo
END SUBROUTINE fill_coordinates_sub

subroutine allocate_coordinates_sub(this,lines)
  !
  ! Subroutine to allocate the class variable arrays
  !
  implicit NONE
  ! Declare calling arguments
  type(line),intent(in),dimension(:) :: lines
  integer :: i,x_total,y_total,x_element,y_element
  class(mesh),intent(out) :: this
  x_total=0
  y_total=0
  x_element=1
  y_element=1
  !
  ! Allocate number of boundaries
  !
  do i =1,size(lines)
    if (lines(i)%get_dimension()==1) then ! x dimension
      x_total=x_total+1
    elseIf (lines(i)%get_dimension()==2) then ! y dimension
      y_total=y_total+1
    endif
  end do
  allocate(this%x_boundaries(1:x_total))
  allocate(this%y_boundaries(1:y_total))
  !
  ! Allocate the number of x and y coordinates
  !
  x_total=0
  y_total=0
  do i =1,size(lines)
    if (lines(i)%get_dimension()==1) then ! x dimension
      x_total=x_total+lines(i)%get_steps()
      this%x_boundaries(x_element)=x_total ! Boundary at the end of each line
      x_element=x_element+1
    elseIf (lines(i)%get_dimension()==2) then ! y dimension
      y_total=y_total+lines(i)%get_steps()
      this%y_boundaries(y_element)=y_total ! Boundary at the end of each line
      y_element=y_element+1
    endif
  end do
  allocate(this%x_coordinates(1:x_total))
  allocate(this%y_coordinates(1:y_total))
end subroutine allocate_coordinates_sub

real(dp) function get_x_fn(this,index)
  !
  ! function to return x value
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  integer,intent(in) :: index
  get_x_fn = this%x_coordinates(index)
end function get_x_fn

real(dp) function get_y_fn(this,index)
  !
  ! function to return y value
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  integer,intent(in) :: index
  get_y_fn = this%y_coordinates(index)
end function get_y_fn

real(dp) function get_x_boundary_fn(this,index)
  !
  ! function to return x boundary
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  integer,intent(in) :: index
  get_x_boundary_fn = this%x_boundaries(index)
end function get_x_boundary_fn

real(dp) function get_y_boundary_fn(this,index)
  !
  ! function to return y boundary
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  integer,intent(in) :: index
  get_y_boundary_fn = this%y_boundaries(index)
end function get_y_boundary_fn

integer function get_x_size_fn(this)
  !
  ! function to return x value
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  get_x_size_fn = size(this%x_coordinates)
end function get_x_size_fn

integer function get_y_size_fn(this)
  !
  ! function to return y value
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  get_y_size_fn = size(this%y_coordinates)
end function get_y_size_fn

integer function number_regions_x_fn(this)
  !
  ! function to return x value
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  number_regions_x_fn = size(this%x_boundaries)
end function number_regions_x_fn

integer function number_regions_y_fn(this)
  !
  ! function to return x value
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  number_regions_y_fn = size(this%y_boundaries)
end function number_regions_y_fn

logical function at_boundary_l_fn(this,index)
  !
  ! function to check if boundary on left side of box
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  integer,intent(in) :: index ! index for the x_coordinate
  integer :: i
  at_boundary_l_fn =.false. ! False by default
  if (index/=1 .and. index/=size(this%x_coordinates)) then
    do i=1,size(this%x_boundaries)
      if((abs(this%x_coordinates(index-1)-this%x_boundaries(i)))<1e-10_dp) at_boundary_l_fn = .true. ! If box to left is a 'boundary' box then true
    end do
  endif  
end function at_boundary_l_fn

logical function at_boundary_r_fn(this,index)
  !
  ! function to check if boundary on right side of box
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  integer,intent(in) :: index ! index for the x_coordinate
  integer :: i
  at_boundary_r_fn =.false. ! False by default
  if (index/=1 .and. index/=size(this%x_coordinates)) then
    do i=1,size(this%x_boundaries)
      if((abs(this%x_coordinates(index)-this%x_boundaries(i)))<1e-10_dp) at_boundary_r_fn = .true. ! If box is a 'boundary' box then true
    end do
  endif  
end function at_boundary_r_fn

logical function at_boundary_t_fn(this,index)
  !
  ! function to check if boundary on top side of box
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  integer,intent(in) :: index ! index for the y_coordinate
  integer :: i
  at_boundary_t_fn =.false. ! False by default
  if (index/=1 .and. index/=size(this%y_coordinates)) then
    do i=1,size(this%y_boundaries)
      if((abs(this%y_coordinates(index)-this%y_boundaries(i)))<1e-10_dp) at_boundary_t_fn = .true. ! If box is a 'boundary' box then true
    end do
  endif  
end function at_boundary_t_fn

logical function at_boundary_b_fn(this,index)
  !
  ! function to check if boundary on bottom side of box
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  integer,intent(in) :: index ! index for the y_coordinate
  integer :: i
  at_boundary_b_fn =.false. ! False by default
  if (index/=1 .and. index/=size(this%y_coordinates)) then
    do i=1,size(this%y_boundaries)
      if((abs(this%y_coordinates(index-1)-this%y_boundaries(i)))<1e-10_dp) at_boundary_b_fn = .true. ! If box below is a 'boundary' box then true
    end do
  endif  
end function at_boundary_b_fn

logical function at_edge_l_fn(this,index)
  !
  ! function to check if at left edge
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  integer,intent(in) :: index ! index for the x_coordinate
  integer :: i
  at_edge_l_fn =.false. ! False by default
  if (index==1) at_edge_l_fn = .true. ! If first node then true
end function at_edge_l_fn

logical function at_edge_r_fn(this,index)
  !
  ! function to check if at right edge
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  integer,intent(in) :: index ! index for the x_coordinate
  integer :: i
  at_edge_r_fn =.false. ! False by default
  if (index==size(this%x_coordinates)) at_edge_r_fn = .true. ! If last box then true
end function at_edge_r_fn

logical function at_edge_t_fn(this,index)
  !
  ! function to check if at top edge
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  integer,intent(in) :: index ! index for the y_coordinate
  integer :: i
  at_edge_t_fn =.false. ! False by default
  if (index==size(this%y_coordinates)) at_edge_t_fn = .true. ! If last box then true
end function at_edge_t_fn

logical function at_edge_b_fn(this,index)
  !
  ! function to check if at bottom edge
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  integer,intent(in) :: index ! index for the y_coordinate
  integer :: i
  at_edge_b_fn =.false. ! False by default
  if (index==1) at_edge_b_fn = .true. ! If first box then true
end function at_edge_b_fn

integer function r_fn(this,indexx,indexy)
  !
  ! function to search a list-index for the last zero or negative element. If not present returns the first element.
  !
  implicit NONE
  ! Declare calling arguments
  class(mesh),intent(in) :: this
  INTEGER,INTENT(IN) :: indexx,indexy
  INTEGER :: i
  r_fn=1
  do i =size(this%x_boundaries),1,-1
    if(this%x_boundaries(i)-indexx==0) then
      r_fn = i ! If the boudnary - index = 0 implies at edge of a region
      exit
    elseif(this%x_boundaries(i)-indexx<0) then
      r_fn = i+1 ! If the boundary - index is negative this implies index is beyond this region.1
      exit
    endif
  end do
  do i =size(this%y_boundaries),1,-1
    if(this%y_boundaries(i)-indexy==0) then
      r_fn = r_fn + ((i-1)*this%number_regions_x()) ! If the boudnary - index = 0 implies at edge of a region but still in it
      exit
    elseif(this%x_boundaries(i)-indexx<0) then
      r_fn = r_fn + ((i)*this%number_regions_x()) ! If the boudnary - index = 0 implies at edge of a region but still in it
      exit
    endif
  end do
end function r_fn

end module mesh_class
