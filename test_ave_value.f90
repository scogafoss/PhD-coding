REAL FUNCTION ave_value ( func, first_value, last_value, n )
!
! Purpose:
! To calculate the average value of function "func" over the
! range [first_value, last_value] by taking n evenly-spaced
! samples over the range, and averaging the results. Function
! "func" is passed to this routine via a dummy argument.
!
! Record of revisions:
! Date Programmer Description of change
! ==== ========== =====================
! 11/24/15 S. J. Chapman Original code
!
IMPLICIT NONE
! Data dictionary: declare calling parameter types & definitions
REAL, EXTERNAL :: func ! Function to be evaluated
REAL, INTENT(IN) :: first_value ! First value in range
REAL, INTENT(IN) :: last_value ! Last value in rnage
INTEGER, INTENT(IN) :: n ! Number of samples to average
! Data dictionary: declare local variable types & definitions
REAL :: delta ! Step size between samples
INTEGER :: i ! Index variable
REAL :: sum ! Sum of values to average
! Get step size.
delta = ( last_value - first_value ) / REAL(n-1)
! Accumulate sum.
sum = 0.
DO i = 1, n
sum = sum + func ( REAL(i-1) * delta )
END DO
! Get average.
ave_value = sum / REAL(n)
END FUNCTION ave_value

PROGRAM test_ave_value
!
! Purpose:
! To test function ave_value by calling it with a user-defined
! function my_func.
!
! Record of revisions:
! Date Programmer Description of change
! ==== ========== =====================
! 11/24/15 S. J. Chapman Original code
!
IMPLICIT NONE
! Data dictionary: declare function types
REAL :: ave_value ! Average value of function
REAL, EXTERNAL :: my_function ! Function to evaluate
! Data dictionary: declare local variable types & definitions
REAL :: ave ! Average of my_function
! Call function with func=my_function.
ave = ave_value ( my_function, 0., 1., 3)
WRITE (*,1000) 'my_function', ave
1000 FORMAT ('The average value of ',A,' between 0. and 1. is ', &
F16.6,'.')
END PROGRAM test_ave_value

REAL FUNCTION my_function( x )
IMPLICIT NONE
REAL, INTENT(IN) :: x
my_function = 3. * x
END FUNCTION my_function