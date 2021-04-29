program test
    use compressed_matrix_class
    use matrix_class
    implicit none
    ! type(compressed_matrix) :: cm
    ! type(matrix) :: m
    ! real(dp), dimension(4) :: vals
    ! integer,dimension(4) :: cols
    ! integer,dimension(3) :: rows
    ! real(dp),dimension(2)::a
    ! real(dp),dimension(2)::b
    ! real(dp),dimension(2)::c
    ! real(dp),DIMENSION(2) :: d
    ! integer :: i
        
    ! ! do i = 1, 10
    ! !         vals(i) = i
    ! ! end do
    ! ! a(1)=0
    ! ! a(2)=3
    ! ! a(3)=6
    ! ! a(4)=9
    ! ! b(1)=1
    ! ! b(2)=4
    ! ! b(3)=7
    ! ! b(4)=10
    ! ! c(1)=2
    ! ! c(2)=5
    ! ! c(3)=8
    ! ! c(4)=0
    ! ! call m%set_variables(a,b,c)
    ! cols(1)=1
    ! cols(2)=2
    ! cols(3)=1
    ! cols(4)=2
    ! vals(1)=4
    ! vals(2)=1
    ! vals(3)=1
    ! vals(4)=3
    ! rows(1)=1
    ! rows(2)=3
    ! rows(3)=5
    ! d(1) = 1
    ! d(2) = 2
    ! a(1)= 0
    ! a(2)=1
    ! b(1)=4
    ! b(2)=3
    ! c(1)=1
    ! c(2)=0
    ! call m%set_variables(a,b,c)
    ! ! cols(1)=1
    ! ! cols(2)=2
    ! ! cols(3)=1
    ! ! cols(4)=2
    ! ! cols(5)=3
    ! ! cols(6)=2
    ! ! cols(7)=3
    ! ! cols(8)=4
    ! ! cols(9)=3
    ! ! cols(10)=4
    ! ! rows(1)=1
    ! ! rows(2)=3
    ! ! rows(3)=6
    ! ! rows(4)=9
    ! ! rows(5)=11
    ! call cm%set_variables(vals,cols,rows,2,2)
    ! Print*, 'CG',cm%conujgate_gradient_method(d)
    ! Print*, 'Thomas',m%thomas_solve(d)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    type(compressed_matrix) :: cm
    type(matrix) :: m
    real(dp), dimension(10) :: vals
    integer,dimension(10) :: cols
    integer,dimension(5) :: rows
    real(dp),dimension(4)::a
    real(dp),dimension(4)::b
    real(dp),dimension(4)::c
    real(dp),DIMENSION(4) :: d
    integer :: i,j
        
    ! vals=1
    ! a(1)=0
    ! a(2)=1
    ! a(3)=1
    ! a(4)=1
    ! b=1
    ! c(1)=1
    ! c(2)=1
    ! c(3)=1
    ! c(4)=0
    ! call m%set_variables(a,b,c)
    ! cols(1)=1
    ! cols(2)=2
    ! cols(3)=1
    ! cols(4)=2
    ! cols(5)=3
    ! cols(6)=2
    ! cols(7)=3
    ! cols(8)=4
    ! cols(9)=3
    ! cols(10)=4
    ! rows(1)=1
    ! rows(2)=3
    ! rows(3)=6
    ! rows(4)=9
    ! rows(5)=11
    ! d(1)=1
    ! d(2)=2
    ! d(3)=3
    ! d(4)=4
    ! call cm%set_variables(vals,cols,rows,4,4)
    call cm%set_size(4,4)
    call cm%add_element(4.4_dp,2,2)
    call cm%add_element(4.4_dp,3,3)
    ! call cm%add_element(4.4_dp,3,4)
    Print*, 'CG',cm%get_element(4,4)
    do i=1,4
        do j=1,4
            print*,i,j,'=',cm%get_element(i,j)
        end do
    end do
    ! Print*, 'Thomas',m%thomas_solve(d)


end program test