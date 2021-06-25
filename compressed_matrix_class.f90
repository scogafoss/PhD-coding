MODULE compressed_matrix_class
    use matrix_class
    ! use timer_class
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
  !!    17/06/2021: Added Jacobi preconditioner                                !!
  !!    22/06/2021: Added Cholesky preconditioner                              !!
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
  procedure,public :: jacobi_solve => jacobi_solve_fn ! Performs a PCG solve
  procedure,public :: print_matrix => print_matrix_sub ! Prints the matrix out
  procedure,public :: forward_substitution => forward_substitution_fn ! Performs forward substitution on lower triangular matrix L
  procedure,public :: backward_substitution => backward_substitution_fn ! Performs forward substitution on upper triangular matrix L^T, L is still input into the function
  procedure,public :: calculate_z => calculate_z_fn ! Calculates z in matrix equation LL^Tz=r if L and r are known
  procedure,public :: cholesky_solve => cholesky_solve_fn ! Performs incomplete cholesky PCG solve
  procedure,public :: replace_element => replace_element_sub ! Replaces element of the matrix
  procedure,public :: jocboi_preconditioner => jacobi_preconditioner_sub ! Allocates the jacobi preconditioner matrix (inverse of diagonal)
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
  PRIVATE :: jacobi_solve_fn
  PRIVATE :: print_matrix_sub
  PRIVATE :: forward_substitution_fn
  PRIVATE :: backward_substitution_fn
  PRIVATE :: calculate_z_fn
  PRIVATE :: cholesky_solve_fn
  PRIVATE :: replace_element_sub
  PRIVATE :: jacobi_preconditioner_sub
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
    ! type(timer) :: time
    !
    ! Timer
    !
    ! call time%start_timer()
    ! Check matrix is compatable with CG method
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
        Ap = this%vector_multiply(basis_vector) ! A * p
        alpha = rsold / (dot_product(basis_vector,Ap))
        solution = solution + alpha * basis_vector
        residual = residual - alpha * Ap
        rsnew = dot_product(residual,residual)
        if (sqrt(rsnew) < convergence) then
            print*,'iterations',i
            exit
        elseif (i==size(source_flux)) then
            print*,'CG convergence not met'
        endif
        basis_vector = residual + (rsnew / rsold) * basis_vector
        rsold = rsnew
    end do
    !
    ! Timer
    !
    ! print*, 'Time to complete CG solve in seconds:', time%elapsed_time()
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
        if (this%get_element(row,column) /= 0) then
            call this%replace_element(value,row,column)
        else
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
        endif
    end if
  end subroutine add_element_sub

  subroutine replace_element_sub(this,value,row,column)
    !
    ! Function to replace an element of the matrix.
    !
    implicit none
    ! Declare calling arguments
    class(compressed_matrix),intent(inout) :: this ! Matrix object
    real(dp),INTENT(IN) :: value
    integer,INTENT(IN) :: row
    integer,INTENT(IN) :: column
    integer :: row_start
    integer :: row_end
    if(row>size(this%row_index)-1) stop 'Row requested is outside of allocated matrix.'
    if(column>size(this%row_index)-1) stop 'Column requested is outside of allocated matrix.' ! I know the matrix will be square so this is correct still.
    row_start = this%row_index(row)
    row_end = this%row_index(row+1)
    if (row_start==row_end) then ! No values in this row
        print*,'No element to replace, adding value'
        call this%add_element(value,row,column)
    else if (find_location_integer(this%column_index(row_start:row_end-1),column)==0) then ! findloc(array,value) returns zero if value is not in array
        ! There is no zero value here - so need add element.
        print*,'No element to replace, adding value'
        call this%add_element(value,row,column)
    else
        this%values(row_start+find_location_integer(this%column_index(row_start:row_end-1),column)-1) = value
    end if
  end subroutine replace_element_sub

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
    ! Only make a change if not a zero
    if (this%get_element(row,column)/=0) then
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
    endif
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

subroutine print_matrix_sub(this)
    !
    ! Sub to set size of matrix
    !
    implicit none
    ! Declare calling arguments
    class(compressed_matrix),INTENT(IN) :: this
    INTEGER :: i,j
    real(dp),DIMENSION(:,:),allocatable :: temp
    allocate(temp(this%get_rows(),this%get_columns()))
    do i=1,this%get_rows()
        do j=1,this%get_columns()
            temp(i,j)=this%get_element(i,j)
            ! print*,'row,col',i,j,this%get_element(i,j),temp(i,j)
        enddo
    enddo
    do i=1,this%get_rows()
        print*,temp(i,:)
    enddo
end subroutine print_matrix_sub

function jacobi_solve_fn(this,source_flux) result(solution)
    !
    ! Function to return solution to Ax=b matrix equation for a square matrix
    !
    implicit none
    ! Declare calling arguments
    class(compressed_matrix),intent(in) :: this ! Matrix object
    real(dp),INTENT(IN),DIMENSION(:) :: source_flux ! Source b in matrix equation
    real(dp), dimension(size(source_flux)) :: solution ! x in matrix equation
    real(dp),dimension(size(source_flux)) :: residual,z ! r in equation
    real(dp) :: rsold ! r'r
    real(dp) :: rsnew ! r'r
    real(dp),dimension(size(source_flux)) :: Ap ! matrix * residual, useful constant
    real(dp),dimension(size(source_flux)) :: basis_vector ! p in equation
    real(dp) :: convergence
    real(dp) :: alpha
    integer :: i,j
    type(compressed_matrix) :: preconditioner ! The inverse of the preconditioner, the inverse of the diagonal of A in Ax=B
    ! Check matrix is compatable with CG method and set the preconditioner at the same time
    call preconditioner%jocboi_preconditioner(this)
    ! Initial values
    convergence = 1e-8_dp
    solution=0_dp
    residual=0_dp ! r
    residual = source_flux - this%vector_multiply(solution)
    z=preconditioner%vector_multiply(residual)
    basis_vector = z ! p
    rsold = dot_product(residual,z) ! numerator in alpha
    ! Loop
    do i=1,size(source_flux) ! Shouldn't have to loop more than the degrees of freedom
        Ap = this%vector_multiply(basis_vector) ! A * p
        alpha = rsold / (dot_product(basis_vector,Ap))
        solution = solution + alpha * basis_vector
        residual = residual - alpha * Ap
        z=preconditioner%vector_multiply(residual)
        rsnew = dot_product(residual,z)
        if (sqrt(rsnew) < convergence) then
            print*,'iterations',i
            exit
        elseif (i==size(source_flux)) then
            print*,'CG convergence not met'
        endif
        basis_vector = z + (rsnew / rsold) * basis_vector
        rsold = rsnew
    end do
    !
    ! Timer
    !
    ! print*, 'Time to complete Jacobi PCG solve in seconds:', time%elapsed_time()
  end function jacobi_solve_fn

  subroutine incomplete_cholesky(a,t)
    !
    ! Subroutine to perform incomplete cholesky factorisation
    !
    implicit none
    type(compressed_matrix),INTENT(IN) :: a
    type(compressed_matrix),INTENT(OUT) :: t
    INTEGER :: n,k,i,j
	n = a%get_columns()
    t=a
    ! call t%set_size(n,n)
	do k=1,n
		call t%add_element(sqrt(t%get_element(k,k)),k,k)
		do i=(k+1),n
		    if (t%get_element(i,k)/=0) call t%add_element(t%get_element(i,k)/t%get_element(k,k),i,k)
        enddo
		do j=(k+1),n
		    do i=j,n
		        if (t%get_element(i,j)/=0) call t%add_element(t%get_element(i,j)-t%get_element(i,k)*t%get_element(j,k),i,j)
		    enddo
		enddo
	enddo
    !Remove values for upper triangle
    do i=1,n
        do j=i+1,n
            call t%remove_element(i,j)
        enddo
    enddo
    !test the matrix
    ! print*,'testing cholesky factorisation'
    ! call t%print_matrix()
    ! print*,'the matrix was'
    ! call a%print_matrix()
    ! do i=1,a%get_rows()
    !     do j=1,a%get_columns()
    !         print*,'row,col',i,j,a%get_element(i,j)
    !     enddo
    ! enddo
  end subroutine incomplete_cholesky

  function forward_substitution_fn(this,b) result(x)
    !
    ! Function to perform forward substitution for Lx=b matrix equation, where L is as a lower triangle matrix
    !
    implicit none
    ! Declare calling arguments
    class(compressed_matrix),intent(in) :: this ! Matrix object
    real(dp),INTENT(IN),DIMENSION(:) :: b ! b in matrix equation Lx=b
    real(dp), dimension(size(b)) :: x ! x in matrix equation
    integer :: i,j
    real(dp) :: summation ! summation over the known terms of x
    ! Check matrix is compatable with vector
    if(this%get_columns() /= size(b)) stop 'Matrix must be same size as vector for forward substitution.'
    if(this%get_element(1,2)/=0) stop 'Matrix must be lower triangular.'
    do i=1,this%get_rows()
        x(i)=1.0_dp/this%get_element(i,i)
        summation = 0.0_dp
        do j=1,i-1
            summation = summation + (this%get_element(i,j)*x(j))
        enddo
        x(i)=x(i)*(b(i)-summation)
    end do
    ! ! Test matrix
    ! call this%print_matrix()
    ! print*,'b',b
    ! print*,'x',x
    ! stop 'testing'
  end function forward_substitution_fn

  function backward_substitution_fn(this,b) result(x)
    !
    ! Function to perform backward substitution for Ux=b matrix equation, where U is as a lower triangle matrix
    ! Note any time accessing element of 'this' compressed matrix, need to swap the i and j since the transpose has not been stored
    !
    implicit none
    ! Declare calling arguments
    class(compressed_matrix),intent(in) :: this ! Matrix object
    real(dp),INTENT(IN),DIMENSION(:) :: b ! b in matrix equation Lx=b
    real(dp), dimension(size(b)) :: x ! x in matrix equation
    integer :: i,j,n
    real(dp) :: summation ! summation over the known terms of x
    ! Check matrix is compatable with vector
    if(this%get_columns() /= size(b)) stop 'Matrix must be same size as vector for forward substitution.'
    if(this%get_element(1,2)/=0) stop 'Matrix must be lower triangular.'
    n=size(b)
    do i=1,this%get_rows()
        x(n-i+1)=1.0_dp/this%get_element(n-i+1,n-i+1)
        summation = 0.0_dp
        do j=1,i-1
            summation = summation + (this%get_element(n-j+1,n-i+1)*x(n-j+1))
            ! print*,'n-j+1,n-i+1',n-j+1,n-i+1,this%get_element(n-j+1,n-i+1),x(n-j+1)
        enddo
        x(n-i+1)=x(n-i+1)*(b(n-i+1)-summation)
    end do
    ! ! ! Test matrix
    ! call this%print_matrix()
    ! print*,'b',b
    ! print*,'x',x
    ! stop 'testing'
  end function backward_substitution_fn

  function cholesky_solve_fn(this,source_flux) result(solution)
    !
    ! Function to return solution to Ax=B matrix equation for a square matrix
    !
    implicit none
    ! Declare calling arguments
    class(compressed_matrix),intent(in) :: this ! Matrix object
    real(dp),INTENT(IN),DIMENSION(:) :: source_flux ! Source B in matrix equation
    real(dp), dimension(size(source_flux)) :: solution ! x in matrix equation
    real(dp),dimension(size(source_flux)) :: residual,z ! r in equation
    real(dp) :: rsold ! r'r
    real(dp) :: rsnew ! r'r
    real(dp),dimension(size(source_flux)) :: Ap ! matrix * residual, useful constant
    real(dp),dimension(size(source_flux)) :: basis_vector ! p in equation
    real(dp) :: convergence
    real(dp) :: alpha
    integer :: i,j
    type(compressed_matrix) :: lower ! The inverse of the preconditioner, the inverse of the diagonal of A in Ax=B
    ! Check matrix is compatable with CG method and set the preconditioner at the same time
    if(this%get_rows() /= this%get_columns()) stop 'Matrix must be square for CG method.'
    ! Initial values
    convergence = 1e-8_dp
    solution=0_dp
    residual=0_dp ! r
    residual = source_flux - this%vector_multiply(solution)
    call incomplete_cholesky(this,lower)
    z=lower%calculate_z(residual)
    basis_vector = z ! p
    rsold = dot_product(residual,z) ! numerator in alpha
    ! Loop
    do i=1,size(source_flux) ! Shouldn't have to loop more than the degrees of freedom
        Ap = this%vector_multiply(basis_vector) ! A * p
        alpha = rsold / (dot_product(basis_vector,Ap))
        solution = solution + alpha * basis_vector
        residual = residual - alpha * Ap
        z=lower%calculate_z(residual)
        rsnew = dot_product(residual,z)
        if (sqrt(rsnew) < convergence) then
            print*,'iterations',i
            exit
        elseif (i==size(source_flux)) then
            print*,'CG convergence not met'
        endif
        basis_vector = z + (rsnew / rsold) * basis_vector
        rsold = rsnew
    end do
    !
    ! Timer
    !
    ! print*, 'Time to complete Jacobi PCG solve in seconds:', time%elapsed_time()
  end function cholesky_solve_fn

  function calculate_z_fn(this,r) result(z)
    !
    ! Function to return solution to LL^Tz=r matrix equation, where L is a lower triangular matrix
    !
    implicit none
    ! Declare calling arguments
    class(compressed_matrix),intent(in) :: this ! Matrix object, L in the matrix equation
    real(dp),INTENT(IN),DIMENSION(:) :: r
    real(dp), dimension(size(r)) :: z
    real(dp),dimension(size(r)) :: y ! temporary vector to find z
    ! Check matrix is compatable with r
    if(this%get_columns()/=size(r)) stop 'Matrix must have same number of columns as the size of residual vector'
    ! Substitute to find z
    ! LL^Tz=r
    ! L^Tz=y
    ! Ly=r
    ! Use r to find y, then use y to find z
    y = this%forward_substitution(r)
    z = this%backward_substitution(y)
  end function calculate_z_fn

  subroutine jacobi_preconditioner_sub(this,cmatrix)
    !
    ! Sub to allocate jacobi preconditioner
    !
    implicit none
    ! Declare calling arguments
    class(compressed_matrix),INTENT(OUT) :: this
    class(compressed_matrix),INTENT(IN) :: cmatrix
    INTEGER :: i,j
    call this%set_size(cmatrix%get_rows(),cmatrix%get_columns())
    if(cmatrix%get_rows() /= cmatrix%get_columns()) stop 'Matrix must be square for CG method.'
    do i=1,cmatrix%get_rows()
        ! Inverse of preconditioner is Pii = 1/Aii
        call this%add_element((1/cmatrix%get_element(i,i)),i,i)
        do j=1, cmatrix%get_columns()
            if(cmatrix%get_element(i,j)/=cmatrix%get_element(j,i)) stop 'Matrix must be symmetric for CG method.'
        end do
    end do
end subroutine jacobi_preconditioner_sub

  END MODULE compressed_matrix_class
