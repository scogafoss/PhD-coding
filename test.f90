program test_sign
    use compressed_matrix_class
    implicit none
    type(compressed_matrix) :: c
    real(dp),DIMENSION(9) ::s
    integer :: i,j
    real(dp) :: value
    call c%set_size(9,9)
    do i=1,9
        value=i
        call c%add_element(value,i,i)
    enddo
    do i=3,9,2
        value=i
        call c%add_element(value,1,i)
        call c%add_element(value,i,1)
    enddo    
    do i=1,9
        value=9+i
        s(i)=value
    enddo
    print*,'solve',c%solve(s)

end program test_sign