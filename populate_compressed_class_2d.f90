MODULE populate_compressed_class_2d
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
    call allocate_matrix(in_matrix,in_mesh)
    !
    ! Perform the actual discretisation
    !
    call this%discretisation(in_matrix,regions,group,in_mesh)
  END SUBROUTINE populate_matrix_sub

  subroutine allocate_matrix(c_matrix,in_mesh)
    !
    ! Subroutine to allocate nuclear matrix
    !
    implicit none
    ! Declare calling arguments
    type(compressed_matrix),INTENT(INOUT) :: c_matrix
    type(mesh),intent(in) :: in_mesh ! Mesh to track the points
    call c_matrix%set_size(in_mesh%get_x_size()*in_mesh%get_y_size(),in_mesh%get_x_size()*in_mesh%get_y_size())
  end subroutine allocate_matrix

  subroutine discretisation_sub(this,c_matrix,regions,group,in_mesh)
    !
    ! Subroutine to fill elements of compressed matrix using discretisation
    !
    implicit none
    ! Declare calling arguments
    class(populate_compressed_2d), intent(INOUT) :: this ! Nuclear_matrix object
    type(compressed_matrix),INTENT(INOUT) :: c_matrix
    type(mesh),intent(in) :: in_mesh ! Mesh to track the points
    type(region_2d),INTENT(IN), DIMENSION(:) :: regions
    INTEGER :: i,j
    real(dp) :: al,ar,at,ab,ac
    integer,INTENT(IN) :: group
    ! Set the variables for the first region.
    do j=1,in_mesh%get_y_size() ! loop all y boxes
      do i=1,in_mesh%get_x_size() ! loop all x boxes
        al=get_al(i,j,in_mesh,regions,group)
        ar=get_ar(i,j,in_mesh,regions,group)
        at=get_at(i,j,in_mesh,regions,group)
        ab=get_ab(i,j,in_mesh,regions,group)
        ac=get_ac(i,j,in_mesh,regions,group,al,ar,at,ab)
        call populate_elements(c_matrix,in_mesh,al,ar,at,ab,ac,i,j)
      enddo
    enddo
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

    real(dp) function get_al(i,j,in_mesh,regions,group)
        !
        ! function to calculate matrix element al(i,j)
        !
        IMPLICIT NONE
        integer,INTENT(IN) :: i,j,group
        type(mesh),INTENT(IN) :: in_mesh
        type(region_2d),DIMENSION(:),INTENT(IN) :: regions
        real(dp) :: deltay,deltax1,deltax2,d1,d2,alpha ! the 1 implies the same i and j and input, the 2 is the i or j shifted by +-1.
        !
        ! Define variables used in calculation
        !
        deltay = regions(in_mesh%r(i,j))%get_delta(2)
        !
        ! Perform calculation
        !
        if(.not.in_mesh%at_edge_l(i)) then
          ! Variables
          deltax1 = regions(in_mesh%r(i,j))%get_delta(1)
          deltax2 = regions(in_mesh%r(i-1,j))%get_delta(1)
          d1 = regions(in_mesh%r(i,j))%get_d(group)
          d2 = regions(in_mesh%r(i-1,j))%get_d(group)
          ! Calculation
          get_al=(-2_dp*deltay)/((deltax2/d2)+(deltax1/d1))
        ! If at edge will have logic later in populate_elements to ignore it. 
        elseif(regions(1)%get_left_boundary()=='r') then ! reflective b.c.
          get_al=0 ! Because Jnet = 0
        elseif(regions(1)%get_left_boundary()=='v') then! vacuum b.c.
          get_al=-deltay/2 ! Because Jin = 0, dphi/dx=phi/2D
        elseif(regions(1)%get_left_boundary()=='z') then! zero flux b.c.
          get_al=-1E3_dp ! Very large?
        elseif(regions(1)%get_left_boundary()=='a') then! albedo b.c.
          alpha = regions(1)%get_left_albedo()
          get_al= -1*((1-alpha)/(1+alpha))*deltay/2
        else
          stop 'Unrecognised left boundary condition'
        endif
      end function get_al

      real(dp) function get_ar(i,j,in_mesh,regions,group)
        !
        ! function to calculate matrix element ar(i,j)
        !
        IMPLICIT NONE
        integer,INTENT(IN) :: i,j,group
        type(mesh),INTENT(IN) :: in_mesh
        type(region_2d),DIMENSION(:),INTENT(IN) :: regions
        real(dp) :: deltay,deltax1,deltax2,d1,d2,alpha ! the 1 implies the same i and j and input, the 2 is the i or j shifted by +-1.
        !
        ! Define variables used in calculation
        !
        deltay = regions(in_mesh%r(i,j))%get_delta(2)
        !
        ! Perform calculation
        !
        if(.not.in_mesh%at_edge_r(i)) then
          ! Variables
          deltax1 = regions(in_mesh%r(i,j))%get_delta(1)
          deltax2 = regions(in_mesh%r(i+1,j))%get_delta(1)
          d1 = regions(in_mesh%r(i,j))%get_d(group)
          d2 = regions(in_mesh%r(i+1,j))%get_d(group)
          ! Calculation
          get_ar=(-2*deltay)/((deltax2/d2)+(deltax1/d1))
        ! If at edge will have logic later in populate_elements to ignore it. 
        elseif(regions(1)%get_right_boundary()=='r') then ! reflective b.c.
          get_ar=0 ! Because Jout = 0
        elseif(regions(1)%get_right_boundary()=='v') then ! vacuum b.c.
          get_ar=-deltay/2 ! Because Jin = 0, dphi/dx=-phi/2D
        elseif(regions(1)%get_right_boundary()=='z') then ! zero flux b.c.
          get_ar=-1E3_dp ! Very large?
        elseif(regions(1)%get_right_boundary()=='a') then ! albedo b.c.
          alpha = regions(1)%get_right_albedo()
          get_ar= -1_dp*((1_dp-alpha)/(1_dp+alpha))*deltay/2
        else
          stop 'Unrecognised right boundary condition'
        endif
      end function get_ar

      real(dp) function get_ab(i,j,in_mesh,regions,group)
        !
        ! function to calculate matrix element ab(i,j)
        !
        IMPLICIT NONE
        integer,INTENT(IN) :: i,j,group
        type(mesh),INTENT(IN) :: in_mesh
        type(region_2d),DIMENSION(:),INTENT(IN) :: regions
        real(dp) :: deltax,deltay1,deltay2,d1,d2,alpha ! the 1 implies the same i and j and input, the 2 is the i or j shifted by +-1.
        !
        ! Define variables used in calculation
        !
        deltax = regions(in_mesh%r(i,j))%get_delta(1)
        !
        ! Perform calculation
        !
        if(.not.in_mesh%at_edge_b(j)) then
          ! Variables
          deltay1 = regions(in_mesh%r(i,j))%get_delta(2)
          deltay2 = regions(in_mesh%r(i,j-1))%get_delta(2)
          d1 = regions(in_mesh%r(i,j))%get_d(group)
          d2 = regions(in_mesh%r(i,j-1))%get_d(group)
          ! Calculation
          get_ab=(-2*deltax)/((deltay2/d2)+(deltay1/d1))
        ! If at edge will have logic later in populate_elements to ignore it. 
        elseif(regions(1)%get_bottom_boundary()=='r') then ! reflective b.c.
          get_ab=0 ! Because Jout = 0
        elseif(regions(1)%get_bottom_boundary()=='v') then ! vacuum b.c.
          get_ab=-deltax/2 ! Because Jin = 0, dphi/dy=phi/2D
        elseif(regions(1)%get_bottom_boundary()=='z') then ! zero flux b.c.
          get_ab=-1E3_dp ! Very large?
        elseif(regions(1)%get_bottom_boundary()=='a') then ! albedo b.c.
          alpha = regions(1)%get_bottom_albedo()
          get_ab= -1*((1-alpha)/(1+alpha))*deltax/2
        else
          stop 'Unrecognised bottom boundary condition'
        endif
      end function get_ab

      real(dp) function get_at(i,j,in_mesh,regions,group)
        !
        ! function to calculate matrix element at(i,j)
        !
        IMPLICIT NONE
        integer,INTENT(IN) :: i,j,group
        type(mesh),INTENT(IN) :: in_mesh
        type(region_2d),DIMENSION(:),INTENT(IN) :: regions
        real(dp) :: deltax,deltay1,deltay2,d1,d2,alpha ! the 1 implies the same i and j and input, the 2 is the i or j shifted by +-1.
        !
        ! Define variables used in calculation
        !
        deltax = regions(in_mesh%r(i,j))%get_delta(1)
        !
        ! Perform calculation
        !
        if(.not.in_mesh%at_edge_t(j)) then
          ! Variables  
          deltay1 = regions(in_mesh%r(i,j))%get_delta(2)
          deltay2 = regions(in_mesh%r(i,j+1))%get_delta(2)
          d1 = regions(in_mesh%r(i,j))%get_d(group)
          d2 = regions(in_mesh%r(i,j+1))%get_d(group)
          ! Calculation
          get_at=(-2_dp*deltax)/((deltay2/d2)+(deltay1/d1))
        ! If at edge will have logic later in populate_elements to ignore it. 
        elseif(regions(1)%get_top_boundary()=='r') then ! reflective b.c.
          get_at=0 ! Because Jout = 0
        elseif(regions(1)%get_top_boundary()=='v') then ! vacuum b.c.
          get_at=-deltax/2_dp ! Because Jin = 0, dphi/dy=-phi/2D
        elseif(regions(1)%get_top_boundary()=='z') then ! zero flux b.c.
          get_at=-1E3_dp ! Very large?
        elseif(regions(1)%get_top_boundary()=='a') then ! albedo b.c.
          alpha = regions(1)%get_top_albedo()
          get_at= -1_dp*((1_dp-alpha)/(1_dp+alpha))*deltax/2_dp
        else
          stop 'Unrecognised top boundary condition'
        endif
      end function get_at

      real(dp) function get_ac(i,j,in_mesh,regions,group,al,ar,at,ab)
        !
        ! function to calculate matrix element ac(i,j)
        !
        IMPLICIT NONE
        integer,INTENT(IN) :: i,j,group
        type(mesh),INTENT(IN) :: in_mesh
        type(region_2d),DIMENSION(:),INTENT(IN) :: regions
        real(dp),INTENT(IN) :: al,ar,at,ab
        real(dp) :: deltay,deltax,removal
        !
        ! Define variables used in calculation
        !
        deltay = regions(in_mesh%r(i,j))%get_delta(2)
        deltax = regions(in_mesh%r(i,j))%get_delta(1)
        removal = regions(in_mesh%r(i,j))%get_removal(group)
        !
        ! Perform calculation
        !
        get_ac = (removal*deltax*deltay) - (al + ar + ab + at)
      end function get_ac

      subroutine populate_elements(c_matrix,in_mesh,al,ar,at,ab,ac,i,j)
        !
        ! Subroutine to fill elements of compressed matrix
        !
        implicit none
        ! Declare calling arguments
        type(compressed_matrix),INTENT(INOUT) :: c_matrix
        type(mesh),intent(in) :: in_mesh ! Mesh to track the points
        real(dp),INTENT(IN) :: al,ar,at,ab,ac
        INTEGER,INTENT(IN) :: i,j ! i is column and j is row
        integer :: nodei ! Stores the row and column node in matrix
        !
        ! First populate the central box node
        !
        nodei = i + ((j-1) * (in_mesh%get_x_size()))
        call c_matrix%add_element(ac,nodei,nodei)
        !
        ! Now check for surrounding boxes, if there are any missing , that matrix element will not be added.
        !
        if (.not.in_mesh%at_edge_l(i)) call c_matrix%add_element(al,nodei,nodei-1) ! Left node will always be just to the left of central
        if (.not.in_mesh%at_edge_r(i)) call c_matrix%add_element(ar,nodei,nodei+1) ! Right node will always be just to the right of central
        if (.not.in_mesh%at_edge_b(j)) call c_matrix%add_element(ab,nodei,nodei-in_mesh%get_x_size()) ! Bottom node will always be the row below, so will be central node-(number of x boxes)
        if (.not.in_mesh%at_edge_t(j)) call c_matrix%add_element(at,nodei,nodei+in_mesh%get_x_size()) ! Top node will always be the row above, so will be central node+(number of x boxes)
      end subroutine populate_elements

  END MODULE populate_compressed_class_2d
