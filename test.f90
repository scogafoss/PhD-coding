real function testing(a)
implicit none
real,INTENT(IN) :: a
testing=2.*a
end function testing
program test
    implicit none
    real :: x = 3.143
    real(8) :: y = 2.33
    real :: testing
    print*,testing(x)
    print *, EPSILON(x)
    print *, EPSILON(y)
end program test