program test
    use precision_set
    use line_class
    use material_class
    use region_class
    use region_class_1d
    use matrix_class
    implicit none
    type(matrix) :: m
    real(dp),allocatable,dimension(:)::a
    real(dp),allocatable,dimension(:)::b
    real(dp),allocatable,dimension(:)::c
    real(dp),allocatable,dimension(:)::s
    allocate(a(1:3))
    allocate(b(1:3))
    allocate(c(1:3))
    allocate(s(1:3))
    a=1
    a(1)=0
    b=2
    c=3
    c(3)=0
    s=4
    call m%set_variables(a,b,c)
    print *,m%thomas_solve(s)
end program test