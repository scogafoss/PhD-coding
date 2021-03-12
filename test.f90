program test
    implicit none
    real, allocatable, dimension(:,:) :: a
    integer :: i,j
    ! do i = 1,2
    !     do j =1,2
    !         if(i==1 .and. j==1) then
    !             a(i,j)=2
    !         else
    !             a(i,j)=1
    !         end if
    !     end do
    ! end do
    allocate(a(1:4,1:4))
    a=1
    print *,a
end program test