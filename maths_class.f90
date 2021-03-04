MODULE maths_class
  use precision_set
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Filename: maths_class.f90                                                !!
!!                                                                           !!
!!  Dependant files: precision_set.f90                                       !!
!!                                                                           !!
!!  Author: Oliver Conway                             Start date: 11/01/2021 !!
!!                                                                           !!
!!  Purpose: Class to operate useful mathematics on input data. This class   !!
!!    will interpolate input y values at desired values of x. The desired x  !!
!!    values must be provided with current x and y values.                   !!
!!                                                                           !!
!!  Revisions:                                                               !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE
! Type definition
TYPE,PUBLIC :: maths ! This will be the name we instantiate
! Instance variables.
PRIVATE
real(dp), allocatable, dimension(:) :: x_desired ! Values of x at which y is desired to be known
real(dp), allocatable, dimension(:) :: x_data ! Current known values of x
real(dp), allocatable, dimension(:) :: y_data ! Corresponding known values of y
CONTAINS
! Bound procedures
PROCEDURE,PUBLIC :: set_variables => set_variables_sub
PROCEDURE,PUBLIC :: get_interpolation => get_interpolation_fn
END TYPE maths
! Restrict access to the actual procedure names
PRIVATE :: set_variables_sub
private :: get_interpolation_fn
! Now add methods
CONTAINS
SUBROUTINE set_variables_sub(this, x_desired, x_data, y_data)
  !
  ! Subroutine to set the variables
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(maths) :: this ! Maths object
  real(dp),INTENT(IN),allocatable, dimension(:) :: x_desired
  real(dp),INTENT(IN),allocatable,dimension(:) :: x_data
  real(dp),INTENT(IN),allocatable,dimension(:) :: y_data
  ! Save data
  this%x_desired = x_desired
  this%x_data = x_data
  this%y_data = y_data
END SUBROUTINE set_variables_sub
FUNCTION get_interpolation_fn(this) result(y_desired)
  !
  ! Function to return interpolated y values as an array
  !
  IMPLICIT NONE
  ! Declare calling arguments
  CLASS(maths),INTENT(IN) :: this ! Maths object
  real(dp), dimension(size(this%x_desired)) :: y_desired
  ! Get interpolation
  integer :: x_iterator ! Iterates over the desired x values.
  integer :: data_iterator ! Iterates over the x data valus to interpolate.
  ! Check min and max x_desired are allowed
  if ((minval(this%x_desired) >= minval(this%x_data)) .and. (maxval(this%x_desired) <= maxval(this%x_data))) then
    ! loop over the desired x values.
    do x_iterator = 1,size(this%x_desired)
      ! Increments the data set of x values until reaches point above desired point.
      data_iterator = 1 ! reset to one before iterating again.
      do while (this%x_desired(x_iterator) > this%x_data(data_iterator))
        data_iterator = data_iterator + 1
      end do
      ! Deals with the interpolation.
      if (this%x_desired(x_iterator) == this%x_data(data_iterator)) then ! If the desired point is identical to one provided then immediately set y_desired.
        y_desired(x_iterator) = this%y_data(data_iterator)
      else ! If desired point is between two points then interpolate
        y_desired(x_iterator) = ((this%x_desired(x_iterator)-this%x_data(data_iterator-1))/(this%x_data(data_iterator) &
        -this%x_data(data_iterator-1)))*(this%y_data(data_iterator)-this%y_data(data_iterator-1))+this%y_data(data_iterator-1)
      end if
    end do
  else
    write (*,*) 'x_desired should be within the range of x_data'
  end if
END FUNCTION get_interpolation_fn
END MODULE maths_class
