MODULE vtk_class
    use mesh_class
    use region_class_2d
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Filename: vtk_class.f90                                                  !!
  !!                                                                           !!
  !!  Dependant files: mesh_class.f90                                          !!
  !!                                                                           !!
  !!  Author: Oliver Conway                             Start date: 16/06/2021 !!
  !!                                                                           !!
  !!  Purpose: Class to write a VTK file to be input to ParaView               !!
  !!                                                                           !!
  !!  Revisions:                                                               !!
  !!                                                                           !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IMPLICIT NONE
  CHARACTER(len=*),PARAMETER :: FMT1 = '(I0)' ! Integer of any length
  CHARACTER(len=*),PARAMETER :: FMT2 = '(ES21.15 E2)' ! Double in scientific form (1.000000000000000E+02)
  ! Type definition
  TYPE,PUBLIC :: vtk ! This will be the name we instantiate
  ! Instance variables.
  CONTAINS
  ! Bound procedures
  PROCEDURE,PUBLIC :: write_vtk => write_vtk_sub ! writes the vtk file
  END TYPE vtk
  ! Restrict access to the actual procedure names
  PRIVATE :: write_vtk_sub
  ! Now add methods
  CONTAINS

  SUBROUTINE write_vtk_sub(this, filename,in_mesh, phi,regions2)
    !
    ! Subroutine to write VTK file
    !
    IMPLICIT NONE
    ! Declare calling arguments
    CLASS(vtk),INTENT(INOUT) :: this ! Matrix object
    CHARACTER(len=*),INTENT(IN) :: filename
    type(mesh),INTENT(IN) :: in_mesh
    real(dp),DIMENSION(:,:),INTENT(IN) :: phi ! 2D array of phi(index,group)
    type(region_2d),DIMENSION(:),INTENT(IN) :: regions2
    INTEGER :: i,j
    !
    ! Open the file for writing
    !
    open(11, file =filename)
    !
    ! Add file header
    !
    call write_header()
    !
    ! Write out coordinates of corners of boxes
    !
    call write_points(in_mesh,regions2)
    !
    ! Write out the cells for the box scheme
    !
    call write_cells(in_mesh)
    !
    ! Write out the cell data
    !
    call write_cell_data(phi)
    !
    ! Close the file
    !
    close(11)
  END SUBROUTINE write_vtk_sub

  SUBROUTINE write_header()
    !
    ! Subroutine to write the header into file
    !
    IMPLICIT NONE
    ! Declare calling arguments
    write(11,'(A)') '# vtk DataFile Version 2.0'//NEW_LINE('A')//'A mesh'//NEW_LINE('A')//'ASCII'//NEW_LINE('A')//'DATASET UNSTRUCTURED_GRID'
    print*,'Added VTK header'
  END SUBROUTINE write_header

  SUBROUTINE write_points(in_mesh,regions2)
    !
    ! Subroutine to write the coordinates into file
    !
    IMPLICIT NONE
    ! Declare calling arguments
    type(mesh),INTENT(IN) :: in_mesh
    type(region_2d),INTENT(IN),DIMENSION(:) :: regions2
    INTEGER :: i,j
    CHARACTER(len=23) :: points,x,y,z
    real(dp) :: deltax,deltay
    write(z,FMT2) 0.0_dp
    ! Write the POINTS header
    write(points,FMT1) (in_mesh%get_x_size()+1)*(in_mesh%get_y_size()+1)
    write(11,'(a)') 'POINTS '//adjustl(trim(points))//' double'
    ! Write the points
    do j=1,in_mesh%get_y_size()
        if(j==1) then ! Bottom corners have extra point (at y-delta/2)
            do i = 1,in_mesh%get_x_size()
                deltay=regions2(in_mesh%r(i,j))%get_delta(2)
                deltax=regions2(in_mesh%r(i,j))%get_delta(1)
                write(y,FMT2) in_mesh%get_y(j)-(deltay/2)
                if(i==1) then ! First corner is before centre of box
                    write(x,FMT2) in_mesh%get_x(i)-(deltax/2)
                    ! Need to write this line now
                    write(11,'(a)') adjustl(trim(x))//' '//adjustl(trim(y))//' '//adjustl(trim(z))
                endif ! All other corners defined after centre
                write(x,FMT2) in_mesh%get_x(i)+(deltax/2)
                write(11,'(a)') adjustl(trim(x))//' '//adjustl(trim(y))//' '//adjustl(trim(z))
            enddo
        endif
        do i = 1,in_mesh%get_x_size()
            deltay=regions2(in_mesh%r(i,j))%get_delta(2)
            deltax=regions2(in_mesh%r(i,j))%get_delta(1)
            write(y,FMT2) in_mesh%get_y(j)+(deltay/2)
            if(i==1) then ! First corner is before centre of box
                write(x,FMT2) in_mesh%get_x(i)-(deltax/2)
                ! Need to write this line now
                write(11,'(a)') adjustl(trim(x))//' '//adjustl(trim(y))//' '//adjustl(trim(z))
            endif ! All other corners defined after centre
            write(x,FMT2) in_mesh%get_x(i)+(deltax/2)
            write(11,'(a)') adjustl(trim(x))//' '//adjustl(trim(y))//' '//adjustl(trim(z))
        enddo
    enddo
    write(11,'(a)') ! Adds a new line at the end of points
    print*,'Added '//adjustl(trim(points))//' VTK box corner points'
  END SUBROUTINE write_points

  SUBROUTINE write_cells(in_mesh)
    !
    ! Subroutine to write the mesh cells indexed from 0
    !
    IMPLICIT NONE
    ! Declare calling arguments
    type(mesh),INTENT(IN) :: in_mesh
    INTEGER :: i,j
    CHARACTER(len=23) :: cells,total,c1,c2,c3,c4,t
    real(dp) :: deltax,deltay
    !
    ! Write the CELLS header
    !
    write(cells,FMT1) in_mesh%get_x_size()*in_mesh%get_y_size()
    write(total,FMT1) in_mesh%get_x_size()*in_mesh%get_y_size()*5
    write(11,'(a)') 'CELLS'//' '//adjustl(trim(cells))//' '//adjustl(trim(total))
    !
    ! Write the points
    !
    write(t,FMT1) 4
    do j=1,in_mesh%get_y_size()
        do i=1,in_mesh%get_x_size()
            write(c1,FMT1) (i+((j-1)*(in_mesh%get_x_size()+1)))-1
            write(c2,FMT1) (i+((j-1)*(in_mesh%get_x_size()+1)))
            write(c3,FMT1) (i+((j)*(in_mesh%get_x_size()+1)))
            write(c4,FMT1) (i+((j)*(in_mesh%get_x_size()+1)))-1
            write(11,'(a)') adjustl(trim(t))//' '//adjustl(trim(c1))//' '//adjustl(trim(c2))//' '//adjustl(trim(c3))//' '//adjustl(trim(c4))
        enddo
    enddo
    write(11,'(a)') ! Add new line at the end
    !
    ! Write the CELL_TYPES
    !
    write(11,'(a)') 'CELL_TYPES'//' '//adjustl(trim(cells))
    write(t,FMT1) 9 ! Type 9 for squares with points anticlockwise
    do i=1,in_mesh%get_x_size()*in_mesh%get_y_size()
        write(11,'(a)') adjustl(trim(t))
    enddo
    write(11,'(a)') ! Add new line at the end
    print*,'Added '//adjustl(trim(cells))//' VTK cells, '//adjustl(trim(total))//' total integers to describe these cells'
  END SUBROUTINE write_cells

  SUBROUTINE write_cell_data(phi)
    !
    ! Subroutine to write the cell data
    !
    IMPLICIT NONE
    ! Declare calling arguments
    real(dp),DIMENSION(:,:),INTENT(IN) :: phi ! 2D array of phi(index,group)
    INTEGER :: i,j
    CHARACTER(len=23) :: cells,group,flux
    !
    ! Write the CELLS_DATA header
    !
    write(cells,FMT1) size(phi(:,1))
    write(11,'(a)') 'CELL_DATA'//' '//adjustl(trim(cells))
    !
    ! Write the data
    !
    ! Whole processs should be repeated for each group
    do i=1,size(phi(1,:))
        ! Write the header for this group
        write(group,FMT1) i
        write(11,'(a)') 'SCALARS FLUX_'//adjustl(trim(group))//' double 1'//NEW_LINE('A')//'LOOKUP_TABLE default'
        do j=1,size(phi(:,1))
            write(flux,FMT2) phi(j,i)
            write(11,'(a)') adjustl(trim(flux))
        enddo
        write(11,'(a)')
        print*,'Added '//adjustl(trim(cells))//' flux values for group '//adjustl(trim(group))
    enddo
  END SUBROUTINE write_cell_data

END MODULE vtk_class
