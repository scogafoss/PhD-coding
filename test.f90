program test
    implicit none
    character(len=23) :: i,x,y,z
    character(len=*),PARAMETER :: FMT = '(I0)'
    character(len=*),PARAMETER :: FMT2 = '(ES21.15 E2)'
    OPEN(12,file='testnewline.vtk')
    write(12,'(A)') '# vtk DataFile Version 2.0'//NEW_LINE('A')//'A mesh'//NEW_LINE('A')//'ASCII'//NEW_LINE('A')//&
    'DATASET UNSTRUCTURED_GRID'//NEW_LINE('A')
    write(i,FMT) 720
    write(x,FMT2) 12.3521
    write(y,FMT2) 1.31e30
    write(z,FMT2) 321.5
    write(12,'(a)') 'POINTS '//adjustl(trim(i))//' double'
    write(12,'(a)')
    write(12,'(a)') adjustl(trim(x))//' '//adjustl(trim(y))//' '//adjustl(trim(z))
    close(12) ! writes a blank line
end program test