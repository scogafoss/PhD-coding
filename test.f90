program test
    implicit none
    real, dimension(2,2) :: a
    integer :: i,j
    do i = 1,2
        do j =1,2
            if(i==1 .and. j==1) then
                a(i,j)=2
            else
                a(i,j)=1
            end if
        end do
    end do
    print *, size(a(2,:))
end program test