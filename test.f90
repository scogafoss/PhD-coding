subroutine incomplete_cholesky(a)
    implicit none
    real,DIMENSION(4,4),INTENT(INOUT) :: a
    INTEGER :: n,k,i,j
	n = size(a,1)
	do k=1,n
		a(k,k) = sqrt(a(k,k));
		do i=(k+1),n
		    if (a(i,k)/=0) a(i,k) = a(i,k)/a(k,k)
        enddo
		do j=(k+1),n
		    do i=j,n
		        if (a(i,j)/=0) a(i,j) = a(i,j)-a(i,k)*a(j,k)
		    enddo
		enddo
	enddo

    do i=1,n
        do j=i+1,n
            a(i,j) = 0
        enddo
    enddo            
endsubroutine incomplete_cholesky
program test
    implicit none
    INTEGER :: i,j
    real :: t
    real,DIMENSION(4,4) :: a
    t=0.0
    do i=1,4
        do j=1,4
            t=t+1.0
            a(i,j)=3.0/t
        enddo
    enddo
    a(4,3) = 0.0
    print*,'matrix a'
    do i=1,4
        print*,a(i,:)
    enddo
    call incomplete_cholesky(a)
    print*,'matrix L'
    do i=1,4
        print*,a(i,:)
    enddo
end program test