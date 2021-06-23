subroutine incomplete_cholesky(a)
    implicit none
    real(8),DIMENSION(3,3),INTENT(INOUT) :: a
    INTEGER :: n,k,i,j
	n = size(a,1)
	do k=1,n
		a(k,k) = sqrt(a(k,k));
        print*,'a(kk)',k,k,a(k,k)
		do i=(k+1),n
		    if (a(i,k)/=0)then
                 a(i,k) = a(i,k)/a(k,k)
                print*,'a(ik)',i,k,a(i,k)
            endif
        enddo
		do j=(k+1),n
		    do i=j,n
		        if (a(i,j)/=0) then
                    a(i,j) = a(i,j)-a(i,k)*a(j,k)
                    print*,'a(ij)',i,j,a(i,j)
                endif
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
    real(8) :: t
    real(8),DIMENSION(3,3) :: a
    t=0.0_8
    a=0.0_8
    a(1,1)=3.6433691756272402_8
    a(2,2)=1.2561896884477530_8
    a(3,3)=1.6128205128205129_8
    a(1,2)=-0.14336917562724016_8
    a(2,1)=-0.14336917562724016_8
    a(2,3)=-0.51282051282051277_8
    a(3,2)=-0.51282051282051277_8
    print*,'matrix a'
    do i=1,3
        print*,a(i,:)
    enddo
    call incomplete_cholesky(a)
    print*,'matrix L'
    do i=1,3
        print*,a(i,:)
    enddo
end program test