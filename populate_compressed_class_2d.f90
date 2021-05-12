MODULE populate_compressed_class_2d
    use compressed_matrix_class
    use populate_class
    use region_class_2d
    use mesh_class
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Filename: populate_compressed_class_2d.f90                               !!
  !!                                                                           !!
  !!  Dependant files: compressed_matrix_class.f90, populate_class.f90,        !!
  !!    region_class_2d.f90, mesh_class.f90                                    !!
  !!                                                                           !!
  !!  Author: Oliver Conway                             Start date: 11/05/2021 !!
  !!                                                                           !!
  !!  Purpose: Class to store and discretise compressed nuclear matrix.        !!
  !!                                                                           !!
  !!  Revisions:                                                               !!
  !!                                                                           !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE
  ! Type definition
  TYPE,PUBLIC,extends(populate) :: populate_compressed_2d
  ! Instance variables.
  CONTAINS
  ! Bound procedures
  procedure,public :: populate_matrix => populate_matrix_sub ! Discretises input region array.
  procedure,public :: discretisation => discretisation_sub ! Performs actual discretisation process
END TYPE populate_compressed_2d
  ! Restrict access to the actual procedure names
  private :: populate_matrix_sub
  private :: discretisation_sub
  CONTAINS

  SUBROUTINE populate_matrix_sub(this,in_matrix,regions,group,in_mesh)
    !
    ! Subroutine to perform discretisation.
    !
    implicit none
    class(populate_compressed_2d), intent(inout) :: this ! Nuclear_matrix object
    type(compressed_matrix),INTENT(INOUT) :: in_matrix
    type(region_2d), intent(in), dimension(:) :: regions ! Region object
    integer,dimension(size(regions)) :: boundary_tracker ! Labels the values of i where boundaries between regions are
    type(mesh),intent(in) :: in_mesh ! Mesh to track the points
    integer :: total_steps ! Total steps across regions.
    INTEGER,INTENT(IN) :: group ! Which group's data is required
    !
    ! Allocate the class array sizes.
    !
    call this%allocate_matrix(in_matrix,in_mesh)
    !
    ! Perform the actual discretisation
    !
    call this%discretisation(in_matrix,regions,boundary_tracker,group)
  END SUBROUTINE populate_matrix_sub

  subroutine allocate_matrix_sub(c_matrix,in_mesh)
    !
    ! Subroutine to allocate nuclear matrix
    !
    implicit none
    ! Declare calling arguments
    type(compressed_matrix),INTENT(INOUT) :: c_matrix
    type(mesh),intent(in) :: in_mesh ! Mesh to track the points
    if(in_mesh%number_regions_x()>=in_mesh%number_regions_y()) then
      call c_matrix%set_size(in_mesh%get_x_size(),in_mesh%get_x_size())
    else
      call c_matrix%set_size(in_mesh%get_y_size(),in_mesh%get_y_size())
    endif
  end subroutine allocate_matrix_sub

  subroutine discretisation_sub(this,c_matrix,regions,boundary_tracker,group,in_mesh)
    !
    ! Subroutine to fill elements of compressed matrix using discretisation
    !
    implicit none
    ! Declare calling arguments
    class(populate_compressed), intent(INOUT) :: this ! Nuclear_matrix object
    type(compressed_matrix),INTENT(INOUT) :: c_matrix
    type(mesh),intent(in) :: in_mesh ! Mesh to track the points
    type(region_1d),INTENT(IN), DIMENSION(:) :: regions
    INTEGER,INTENT(IN), DIMENSION(:) :: boundary_tracker
    INTEGER :: region_iterator
    INTEGER :: i
    real(dp) :: delta
    real(dp) :: absorption
    real(dp) :: D
    real(dp) :: temp
    real(dp) :: temp2
    integer,INTENT(IN) :: group
    ! Set the variables for the first region.
    delta = regions(1)%get_delta()
    absorption = regions(1)%get_removal(group)
    D = 1 / (3 * (regions(1)%get_absorption(group)+(regions(1)%get_scatter(group))))
    region_iterator = 1
    DO i = 1 , c_matrix%get_rows()
      ! Check if at a boundary, where average values are needed. No need for last value of i.
      if (i == boundary_tracker(region_iterator) .and. region_iterator < size(regions)) then
        ! D = ((1/(3*absorption))*delta+(1/(3*(regions(region_iterator+1)%get_absorption()))*(regions(region_iterator+1)%get_length()&
        ! /(regions(region_iterator+1)%get_steps()))))/(delta+(regions(region_iterator+1)%get_length()&
        ! /(regions(region_iterator+1)%get_steps())))
        absorption = (((absorption*delta)+(regions(region_iterator+1)%get_removal(group)*regions(region_iterator+1)%get_delta())))/&
        (delta+(regions(region_iterator+1)%get_delta()))
        delta = (delta + (regions(region_iterator+1)%get_delta()))/2
        ! region_iterator=region_iterator+1
      end if
      ! Check if just after a boudary, where values are updated to new region
      if (i==boundary_tracker(region_iterator)+1 .and. region_iterator < size(regions)) then
        region_iterator=region_iterator+1
        delta = regions(region_iterator)%get_delta()
        absorption = regions(region_iterator)%get_removal(group)
        D = 1 / (3 * (regions(region_iterator)%get_absorption(group)+(regions(region_iterator)%get_scatter(group))))
      end if
    !
    !------------------------------------------------------------------------------
    !
      ! At each boundary need more terms
      if(i == boundary_tracker(region_iterator) .and. i /= c_matrix%get_rows()) then
        ! ai,i-1
        temp = -((1/(3*(regions(region_iterator)%get_absorption(group)+&
        (regions(region_iterator)%get_scatter(group)))))/(regions(region_iterator)%get_delta()))
        call c_matrix%add_element(temp,i,i-1)
        ! ai,i
        temp = absorption*delta + ((1/(3*(regions(region_iterator)%get_absorption(group)+&
        (regions(region_iterator)%get_scatter(group)))))/(regions(region_iterator)%get_delta()))+&
        ((1/(3*(regions(region_iterator+1)%get_absorption(group)+&
        (regions(region_iterator+1)%get_scatter(group)))))/(regions(region_iterator+1)%get_delta()))
        call c_matrix%add_element(temp,i,i)
        ! ai,i+1
        temp = -((1/(3*(regions(region_iterator+1)%get_absorption(group)+&
        (regions(region_iterator+1)%get_scatter(group)))))/(regions(region_iterator+1)%get_delta()))
        call c_matrix%add_element(temp,i,i+1)
    !
    !------------------------------------------------------------------------------
    !
      ! Elsewhere need to correct for geometry and position
      ELSE
        ! ai,i-1
        if (i==1) then
          temp= -(1 / (3 * (regions(size(regions))%get_absorption(group)+(regions(size(regions))%get_scatter(group)))))/&
          (regions(size(regions))%get_delta()) ! -D_last/delta_last
          call c_matrix%add_element(temp,i,c_matrix%get_columns())
          call c_matrix%add_element(temp,c_matrix%get_rows(),i) ! Both corners are equal
          temp2= -(1 / (3 * (regions(1)%get_absorption(group)+(regions(1)%get_scatter(group)))))/&
          (regions(size(regions))%get_delta()) !-D_first/delta_first
          call c_matrix%add_element(temp2,i,i+1)
          temp = (-temp)+(-temp2)+(((regions(size(regions))%get_delta()+regions(1)%get_delta())/2)*&
          weighted_average(regions(size(regions))%get_delta(),regions(1)%get_delta(),regions(size(regions))%get_removal(group),&
          regions(1)%get_removal(group)))
          call c_matrix%add_element(temp,i,i)

        else if (i==c_matrix%get_rows()) then
          temp = -(1 / (3 * (regions(size(regions))%get_absorption(group)+(regions(size(regions))%get_scatter(group)))))/&
          (regions(size(regions))%get_delta())
          call c_matrix%add_element(temp,i,i-1)
          temp = (-2*temp)+(regions(size(regions))%get_delta()*regions(size(regions))%get_removal(group))
          call c_matrix%add_element(temp,i,i)

        else
          temp = -(D/(delta))
          ! ai,i-1
          call c_matrix%add_element(temp,i,i-1)
          ! ai,i+1
          call c_matrix%add_element(temp,i,i+1)
          ! ai,i
          temp = (absorption*delta) + (2*D/(delta))
          call c_matrix%add_element(temp,i,i)
        end if
      END IF
    END DO
  end subroutine discretisation_sub

  real(dp) function weighted_average(weight1,weight2,variable1,variable2)
        !
        ! function to calculate weighted average
        !
        IMPLICIT NONE
        real(dp),INTENT(IN) :: weight1
        real(dp),INTENT(IN) :: weight2
        real(dp),INTENT(IN) :: variable1
        real(dp),INTENT(IN) :: variable2
        weighted_average=((variable1*weight1)+(variable2*weight2))/(weight1+weight2)
    end function weighted_average

  END MODULE populate_compressed_class
