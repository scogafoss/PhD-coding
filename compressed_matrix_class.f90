MODULE compressed_matrix_class
    use matrix_class
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Filename: compressed_matrix_class.f90                                    !!
  !!                                                                           !!
  !!  Dependant files: precision_set.f90                                       !!
  !!                                                                           !!
  !!  Author: Oliver Conway                             Start date: 19/04/2021 !!
  !!                                                                           !!
  !!  Purpose: Class to store compressed matrix as three arrays: vlaues,       !!
  !!    column_index, and row_index.                                           !!
  !!                                                                           !!
  !!  Revisions:                                                               !!
  !!    28/04/2021: Added to part of the matrix hierarchy                      !!
  !!                                                                           !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE
  ! Type definition
  TYPE,PUBLIC,extends(matrix) :: compressed_matrix ! This will be the name we instantiate
  ! Instance variables.
  real(dp), allocatable, dimension(:) :: values ! Non-zero values in matrix
  integer, allocatable, dimension(:) :: column_index ! Column index
  integer, allocatable, dimension(:) :: row_index ! Stores the total number of non-zero values summed from previous rows (plus an extra element at the start)
  integer :: rows=0 ! Stores the number of rows
  integer :: columns=0 !Stores the number of columns
  CONTAINS
  ! Bound procedures
  PROCEDURE,PUBLIC :: set_variables => set_variables_sub ! Allows user to input desired matrix
  procedure,public :: vector_multiply => vector_multiply_fn ! Multiply two compressed matrices
  PROCEDURE,PUBLIC :: solve => solve_fn ! Performs CG solve
  PROCEDURE,PUBLIC :: get_element => get_element_fn ! Returns element
  procedure,public :: add_element => add_element_sub ! Adds element to matrix
  procedure,public :: remove_element => remove_element_sub ! Removes element from matrix
  procedure,public :: set_size => set_size_sub ! Set rows and columns dimensions.
  procedure,public :: get_rows => get_rows_fn ! Function to return number of rows
  procedure,public :: get_columns => get_columns_fn ! Function to return number of columns
  procedure,public :: print_all => print_all_sub ! prints the values, then columns, then rows.
  END TYPE compressed_matrix
  ! Restrict access to the actual procedure names
  PRIVATE :: set_variables_sub
  private :: vector_multiply_fn
  private :: solve_fn
  private :: get_element_fn
  private :: add_element_sub
  private :: remove_element_sub
  private :: set_size_sub
  private :: get_rows_fn
  private :: get_columns_fn
  private :: print_all_sub
  ! Now add methods
  CONTAINS

  SUBROUTINE set_variables_sub(this, values, column_index, row_index, rows, columns)
    !
    ! Subroutine to set the variables
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(compressed_matrix) :: this ! Matrix object
    real(dp),INTENT(IN),dimension(:) :: values
    integer,INTENT(IN),dimension(:) :: column_index
    integer,INTENT(IN),dimension(:) :: row_index
    integer,INTENT(IN) :: rows
    integer,INTENT(IN) :: columns
    ! Save data
    this%values = values
    this%column_index = column_index
    this%row_index = row_index
    this%rows = rows
    this%columns = columns
  END SUBROUTINE set_variables_sub

  function vector_multiply_fn(this, vector) result(solution)
    !
    ! Function to return multiplication of two matrices.
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(compressed_matrix),INTENT(IN) :: this ! Matrix object
    real(dp), INTENT(IN), dimension(:) :: vector
    real(dp),DIMENSION(size(vector)) :: solution
    integer :: i
    integer :: j
    integer :: j_index
    if (size(vector)/=this%columns) stop 'Incompatible dimensions between vector and matrix.'
    solution(:)=0_dp
    do i = 1, this%rows
        do j = 1, this%row_index(i+1)-this%row_index(i)
            j_index=this%column_index(j+this%row_index(i)-1)
            solution(i)=solution(i)+(this%get_element(i,j_index)*vector(j_index))
        end do
    end do
END function vector_multiply_fn

  function solve_fn(this,source_flux) result(solution)
    !
    ! Function to return solution to Ax=B matrix equation for a square matrix
    !
    implicit none
    ! Declare calling arguments
    class(compressed_matrix),intent(in) :: this ! Matrix object
    real(dp),INTENT(IN),DIMENSION(:) :: source_flux ! Source B in matrix equation
    real(dp), dimension(size(source_flux)) :: solution ! x in matrix equation
    real(dp),dimension(size(source_flux)) :: residual ! r in equation
    real(dp) :: rsold ! r'r
    real(dp) :: rsnew ! r'r
    real(dp),dimension(size(source_flux)) :: Ap ! matrix * residual, useful constant
    real(dp),dimension(size(source_flux)) :: basis_vector ! p in equation
    real(dp) :: convergence
    real(dp) :: alpha
    integer :: i,j
    ! Check matrix is compatable with CG method
    !!!!!!!!!!!
    ! do i=1,this%rows
    !     do j=1,this%columns
    !         print*,i,j,'=',this%get_element(i,j)
    !     enddo
    ! enddo
    ! print*,source_flux
    !!!!!!!!!!!
    if(this%rows /= this%columns) stop 'Matrix must be square for CG method.'
    do i=1,this%rows
        do j=1, this%columns
            if(this%get_element(i,j)/=this%get_element(j,i)) stop 'Matrix must be symmetric for CG method.'
        end do
    end do
    ! Initial values
    convergence = 1e-8_dp
    solution=0_dp
    residual=0_dp
    residual = source_flux - this%vector_multiply(solution)
    basis_vector = residual
    rsold = dot_product(residual,residual)
    ! Loop
    do i=1,size(source_flux) ! Shouldn't have to loop more than the degrees of freedom
        Ap = this%vector_multiply(basis_vector)
        alpha = rsold / (dot_product(basis_vector,Ap))
        solution = solution + alpha * basis_vector
        residual = residual - alpha * Ap
        rsnew = dot_product(residual,residual)
        if (sqrt(rsnew) < convergence) then
            print*,'CG convergence met'
            exit
        endif
        basis_vector = residual + (rsnew / rsold) * basis_vector
        rsold = rsnew
    end do
  end function solve_fn

  real (dp) function get_element_fn(this,row,column)
    !
    ! Function to return an element of the matrix.
    !
    implicit none
    ! Declare calling arguments
    class(compressed_matrix),intent(in) :: this ! Matrix object
    integer,INTENT(IN) :: row
    integer,INTENT(IN) :: column
    integer :: row_start
    integer :: row_end
    if(row>size(this%row_index)-1) stop 'Row requested is outside of allocated matrix.'
    if(column>size(this%row_index)-1) stop 'Column requested is outside of allocated matrix.' ! I know the matrix will be square so this is correct still.
    row_start = this%row_index(row)
    row_end = this%row_index(row+1)
    if (row_start==row_end) then ! No values in this row
        get_element_fn = 0
    else if (find_location_integer(this%column_index(row_start:row_end-1),column)==0) then ! findloc(array,value) returns zero if value is not in array
        ! There is no non-zero value here - so return zero.
        get_element_fn = 0
    else
        get_element_fn=this%values(row_start+find_location_integer(this%column_index(row_start:row_end-1),column)-1)
    end if
  end function get_element_fn

  subroutine add_element_sub(this,value,row,column)
    !
    ! Subroutine to set the variables
    !
    IMPLICIT NONE
    ! Declare calling arguments
    class(compressed_matrix),INTENT(INOUT) :: this
    real(dp),INTENT(IN) :: value
    integer, INTENT(IN) :: row
    integer, INTENT(IN) :: column
    real(dp),dimension(size(this%values)) :: temp_values
    integer,dimension(size(this%values)) :: temp_columns
    integer :: i
    integer :: new_value_index
    ! Check within matrix
    if(this%rows==0 .or. this%columns==0) stop 'Size of matrix not yet set'
    if(row>this%rows) stop 'Row requested is outside of allocated matrix.'
    if(column>this%columns) stop 'Column requested is outside of allocated matrix.'
    ! If not yet allocated (i.e. only size set) then allocate and add one point.
    if(.not. allocated(this%values)) then
        allocate(this%values(1:1))
        allocate(this%column_index(1:1))
        this%values=value
        this%column_index=column
        this%row_index(row+1:this%rows+1)=2
    ! If already allocated
    else
        if (this%get_element(row,column) /= 0) stop 'Requested row and column already filled by non-zero value.'
        ! Add to size of stored arrays.
        temp_values=this%values
        temp_columns=this%column_index
        deallocate(this%values)
        deallocate(this%column_index)
        allocate(this%values(1:size(temp_values)+1))
        allocate(this%column_index(1:size(temp_columns)+1))
        do i = 1, (this%row_index(row+1)-this%row_index(row)) ! Loop over the number of points in this row
            if (column < this%column_index(i+this%row_index(row))) exit ! If here then the added value needs to go at this index and shift other values along.
        end do
        if (i==(this%row_index(row+1)-this%row_index(row))) i=i+1 ! Added value goes at the end of the row
        this%row_index(row+1:size(this%row_index))=this%row_index(row+1:size(this%row_index))+1 ! Add one to the row index according to location of new value.
        ! Populate the new column and value vectors
        new_value_index=i+this%row_index(row)-1 ! Record position of new value
        ! Shift columns
        this%column_index(1:new_value_index-1)=temp_columns(1:new_value_index-1)
        this%column_index(new_value_index)=column
        this%column_index(new_value_index+1:size(this%column_index))=temp_columns(new_value_index:size(temp_columns))
        ! Shift values
        this%values(1:new_value_index-1)=temp_values(1:new_value_index-1)
        this%values(new_value_index)=value
        this%values(new_value_index+1:size(this%values))=temp_values(new_value_index:size(temp_values))
    end if
  end subroutine add_element_sub

  subroutine remove_element_sub(this,row,column)
    !
    ! Subroutine to set the variables
    !
    IMPLICIT NONE
    ! Declare calling arguments
    class(compressed_matrix),INTENT(INOUT) :: this
    integer, INTENT(IN) :: row
    integer, INTENT(IN) :: column
    real(dp),dimension(size(this%values)) :: temp_values
    integer,dimension(size(this%values)) :: temp_columns
    integer :: i
    integer :: old_value_index
    ! Check within matrix
    if(row>this%rows) stop 'Row requested is outside of allocated matrix.'
    if(column>this%columns) stop 'Column requested is outside of allocated matrix.'
    if (this%get_element(row,column)==0) stop 'No non-zero value present at requested coordinate.'
    ! Add to size of stored arrays.
    temp_values=this%values
    temp_columns=this%column_index
    deallocate(this%values)
    deallocate(this%column_index)
    allocate(this%values(1:size(temp_values)-1))
    allocate(this%column_index(1:size(temp_columns)-1))
    do i = 1, (this%row_index(row+1)-this%row_index(row)) ! Loop over the number of points in this row
        if (temp_columns(i+this%row_index(row)-1)==column) exit ! This is the value to remove
    end do
    this%row_index(row+1:size(this%row_index))=this%row_index(row+1:size(this%row_index))-1 ! Remove one from row index according to location of new value.
    ! Populate the new column and value vectors
    old_value_index=i+this%row_index(row)-1 ! Record position of removed value
    ! Shift columns
    this%column_index(1:old_value_index-1)=temp_columns(1:old_value_index-1)
    this%column_index(old_value_index:size(this%column_index))=temp_columns(old_value_index+1:size(temp_columns))
    ! Shift values
    this%values(1:old_value_index-1)=temp_values(1:old_value_index-1)
    this%values(old_value_index:size(this%values))=temp_values(old_value_index+1:size(temp_values))
end subroutine remove_element_sub

integer function find_location_integer(integer_array,value)
    !
    ! Function to replace the findloc fortran function
    !
    implicit none
    ! Declare calling arguments
    integer,INTENT(IN),DIMENSION(:) :: integer_array
    integer,INTENT(IN) :: value
    integer :: i
    do i = 1, size(integer_array)
        if(integer_array(i)==value) then
            find_location_integer=i
            exit
        else if (i==size(integer_array)) then ! If here it has searched through whole array and not found value.
            find_location_integer=0
        end if
    end do
end function find_location_integer

function transpose_vector(vector) result(solution)
    !
    ! Function to transpose 1D array
    !
    implicit none
    ! Declare calling arguments
    real(dp),INTENT(IN),DIMENSION(:) :: vector
    real(dp),DIMENSION(1,size(vector)) :: solution
    solution(1,:)=vector
end function transpose_vector

subroutine set_size_sub(this,rows,columns)
    !
    ! Sub to set size of matrix
    !
    implicit none
    ! Declare calling arguments
    class(compressed_matrix),INTENT(INOUT) :: this
    INTEGER, INTENT(IN) :: rows
    INTEGER, INTENT(IN) :: columns
    this%rows = rows
    this%columns = columns
    allocate(this%row_index(1:rows+1))
    this%row_index=1
end subroutine set_size_sub

integer function get_rows_fn(this)
    !
    ! Function to return number of rows
    !
    implicit none
    ! Declare calling arguments
    class(compressed_matrix),INTENT(IN) :: this
    get_rows_fn=this%rows
end function get_rows_fn

integer function get_columns_fn(this)
    !
    ! Function to return number of columns
    !
    implicit none
    ! Declare calling arguments
    class(compressed_matrix),INTENT(IN) :: this
    get_columns_fn=this%columns
end function get_columns_fn

subroutine print_all_sub(this)
    !
    ! Sub to set size of matrix
    !
    implicit none
    ! Declare calling arguments
    class(compressed_matrix),INTENT(INOUT) :: this
    print*,this%values
    print*,this%column_index
    print*,this%row_index
end subroutine print_all_sub

  END MODULE compressed_matrix_class
