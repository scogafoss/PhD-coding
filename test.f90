subroutine testing(a,b)
    implicit none
    integer,INTENT(INOUT) :: a
    integer,INTENT(IN),OPTIONAL :: b
    ! call subtest(a,b)
    if(present(b)) then
        a=b
    else
        a=3
    endif
end subroutine testing

! subroutine subtest(a,b)
!     implicit none
!     INTEGER,INTENT(INOUT) :: a
!     INTEGER,INTENT(IN),OPTIONAL :: b
    ! if(present(b)) then
    !     a=b
    ! else
    !     a=3
    ! endif
! end subroutine subtest

program test_sign
    implicit none
    integer :: a,b
    b=2
    call testing(a,b)
    print *, a
end program test_sign