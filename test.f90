program test
    implicit none
    real, dimension(2,2) :: a
    integer :: i,j
    do i = 1,2
        do j =1,2
            a(i,j)=1
        end do
    end do
    print *, sum(a)
end program test