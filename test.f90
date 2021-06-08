program test_sign
    use compressed_matrix_class
    implicit none
    type(compressed_matrix) :: c
    real(dp),DIMENSION(9) ::s
    integer :: i,j
    real(dp) :: value
    call c%set_size(4,4)
    value=1_dp
    do i=1,4
        do j=1,4
            call c%add_element(value,i,j)
            value=value+1_dp
        enddo
    enddo
    do i=1,4
        do j=1,4
            print*, 'row',i,'column',j,c%get_element(i,j)
        enddo
    enddo

end program test_sign