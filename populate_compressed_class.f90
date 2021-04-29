MODULE populate_compressed_class
    use compressed_matrix_class
    use populate_class
    use region_class_1d
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Filename: populate_compressed_class.f90                                  !!
  !!                                                                           !!
  !!  Dependant files: compressed_matrix_class.f90, populate_class.f90,        !!
  !!    region_class_1d.f90                                                    !!
  !!                                                                           !!
  !!  Author: Oliver Conway                             Start date: 23/04/2021 !!
  !!                                                                           !!
  !!  Purpose: Class to store and discretise compressed nuclear matrix. This   !!
  !!    matrix class will only be used for the periodic boundary condition.    !!
  !!    This condition only works for slab geometry.                           !!
  !!                                                                           !!
  !!  Revisions:                                                               !!
  !!    28/04/2021: Added to populate class hierarchy                          !!
  !!                                                                           !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE
  ! Type definition
  TYPE,PUBLIC,extends(populate_class) :: populate_compressed ! This will be the name we instantiate
  ! Instance variables.
  CONTAINS
  ! Bound procedures
  procedure,public :: populate => populate_sub ! Discretises input region array.
  procedure,public :: allocate_matrix => allocate_matrix_sub ! Allocates the matrix within discretise_regions_sub
  procedure,public :: discretisation => discretisation_sub ! Performs actual discretisation process
END TYPE populate_compressed
  ! Restrict access to the actual procedure names
  private :: discretise_regions_sub
  private :: allocate_matrix_sub
  private :: discretisation_sub
  CONTAINS

  SUBROUTINE populate_sub(this,c_matrix,regions,group)
    !
    ! Subroutine to perform discretisation.
    !
    implicit none
    class(compressed_nuclear_matrix), intent(out) :: this ! Nuclear_matrix object
    type(compressed_matrix),INTENT(INOUT) :: c_matrix
    type(region_1d), intent(in), dimension(:) :: regions ! Region object
    integer,dimension(size(regions)) :: boundary_tracker ! Labels the values of i where boundaries between regions are
    integer :: total_steps ! Total steps across regions.
    INTEGER,INTENT(IN) :: group ! Which group's data is required
    !
    ! Allocate the class array sizes.
    !
    call this%allocate_matrix(c_matrix,total_steps,regions,boundary_tracker)
    !
    ! Perform the actual discretisation
    !
    call this%discretisation(c_matrix,regions,boundary_tracker,group)
  END SUBROUTINE populate_sub

  subroutine allocate_matrix_sub(this,c_matrix,total_steps,regions,boundary_tracker)
    !
    ! Subroutine to allocate nuclear matrix
    !
    implicit none
    ! Declare calling arguments
    class(compressed_nuclear_matrix), intent(out) :: this ! Nuclear_matrix object
    type(compressed_matrix),INTENT(INOUT) :: c_matrix
    integer,INTENT(OUT) :: total_steps
    type(region_1d),INTENT(IN), DIMENSION(:) :: regions
    INTEGER,INTENT(OUT), DIMENSION(:) :: boundary_tracker
    INTEGER :: region_iterator
    total_steps = 0
    do region_iterator =1,size(regions)
      total_steps = total_steps + regions(region_iterator)%get_steps()
      if(region_iterator==size(regions)) then
        boundary_tracker(region_iterator) = total_steps ! Remove the very last node, store the second to last node (the last node is same as first so no need to repeat)
      else
        boundary_tracker(region_iterator) = total_steps + 1 ! Stores the boundary between regions
      end if
    end do
    call c_matrix%set_size(total_steps,total_steps)
  end subroutine allocate_matrix_sub

  subroutine discretisation_sub(this,c_matrix,regions,boundary_tracker,group)
    !
    ! Subroutine to fill elements of compressed matrix using discretisation
    !
    implicit none
    ! Declare calling arguments
    class(compressed_nuclear_matrix), intent(INOUT) :: this ! Nuclear_matrix object
    type(compressed_matrix),INTENT(INOUT) :: c_matrix
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
    delta = regions(1)%get_length()/regions(1)%get_steps()
    absorption = regions(1)%get_removal(group)
    D = 1 / (3 * (regions(1)%get_absorption(group)+(regions(1)%get_scatter(group))))
    region_iterator = 1
    DO i = 1 , c_matrix%get_rows()
      ! Check if at a boundary, where average values are needed. No need for last value of i.
      if (i == boundary_tracker(region_iterator) .and. region_iterator < size(regions)) then
        ! D = ((1/(3*absorption))*delta+(1/(3*(regions(region_iterator+1)%get_absorption()))*(regions(region_iterator+1)%get_length()&
        ! /(regions(region_iterator+1)%get_steps()))))/(delta+(regions(region_iterator+1)%get_length()&
        ! /(regions(region_iterator+1)%get_steps())))
        absorption = ((absorption*delta)+(regions(region_iterator+1)%get_removal(group)*regions(region_iterator+1)%get_length()&
        /(regions(region_iterator+1)%get_steps())))/(delta+(regions(region_iterator+1)%get_length()&
        /(regions(region_iterator+1)%get_steps())))
        delta = (delta + (regions(region_iterator+1)%get_length()/(regions(region_iterator+1)%get_steps())))/2
        ! region_iterator=region_iterator+1
      end if
      ! Check if just after a boudary, where values are updated to new region
      if (i==boundary_tracker(region_iterator)+1 .and. region_iterator < size(regions)) then
        region_iterator=region_iterator+1
        delta = regions(region_iterator)%get_length()/regions(region_iterator)%get_steps()
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
        (regions(region_iterator)%get_scatter(group)))))/(regions(region_iterator)%get_length()/&
        regions(region_iterator)%get_steps()))
        call c_matrix%add_element(temp,i,i-1)
        ! ai,i
        temp = absorption + ((1/(3*(regions(region_iterator)%get_absorption(group)+&
        (regions(region_iterator)%get_scatter(group)))))/(regions(region_iterator)%get_length()/&
        regions(region_iterator)%get_steps()))+((1/(3*(regions(region_iterator+1)%get_absorption(group)+&
        (regions(region_iterator+1)%get_scatter(group)))))/(regions(region_iterator+1)%get_length()/&
        regions(region_iterator+1)%get_steps()))
        call c_matrix%add_element(temp,i,i)
        ! ai,i+1
        temp = -((1/(3*(regions(region_iterator+1)%get_absorption(group)+&
        (regions(region_iterator+1)%get_scatter(group)))))/(regions(region_iterator+1)%get_length()/&
        regions(region_iterator+1)%get_steps()))
        call c_matrix%add_element(temp,i,i+1)
    !
    !------------------------------------------------------------------------------
    !
      ! Elsewhere need to correct for geometry and position
      ELSE
        ! ai,i-1
        if (i==1) then
          temp= -(1 / (3 * (regions(size(regions))%get_absorption(group)+(regions(size(regions))%get_scatter(group)))))/&
          (regions(size(regions))%get_delta())
          call c_matrix%add_element(temp,i,c_matrix%get_columns())
          call c_matrix%add_element(temp,c_matrix%get_rows(),i) ! Both corners are equal
          temp2= -(1 / (3 * (regions(1)%get_absorption(group)+(regions(1)%get_scatter(group)))))/&
          (regions(size(regions))%get_delta())
          call c_matrix%add_element(temp,i,i)
          temp = (-temp)+(-temp2)+(((regions(size(regions))%get_delta()+regions(1)%get_delta())/2)*&
          weighted_average(regions(size(regions))%get_delta(),regions(1)%get_delta(),regions(size(regions))%get_removal(group),&
          regions(1)%get_removal(group)))
          call c_matrix%add_element(temp,i,i+1)

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
          temp = absorption + (D/(delta))
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

  END MODULE compressed_nuclear_matrix_class
