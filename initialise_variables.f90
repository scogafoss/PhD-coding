module initialise_variables
  use region_class_1d
  use region_class_2d
  use mesh_class
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Filename: initialise_variables.f90                                       !!
!!                                                                           !!
!!  Dependant files: precision_set.f90, line_class.f90, material_class.f90,  !!
!!    region_class_1d.f90                                                    !!
!!                                                                           !!
!!  Author: Oliver Conway                             Start date: 11/02/2021 !!
!!                                                                           !!
!!  Purpose: intialises multiregion set_variables_sub                        !!
!!                                                                           !!
!!  Revisions:                                                               !!
!!    17/02/2021: Updated how input deck is read, due to added nu sigma_f    !!
!!    09/03/2021: Add funcitonality to calculate source flux for each group, !!
!!      and checks if source is fission or volumetric, source flux matrix.   !!
!!    12/03/2021: Corrected material ID reading to not repeat IDs            !!
!!    10/05/2021: Check for 2D conditions                                    !!
!!    21/05/2021: Subdivided code into subroutines                           !!
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE
contains

  subroutine initialise(filename,lines,materials,source_flux,total_groups,regions,regions2,in_mesh)
    character(len=*), intent(in) :: filename
    type(region_1d), intent(inout), allocatable, dimension(:) :: regions
    type(region_2d), intent(inout), allocatable, dimension(:) :: regions2
    type(line),  intent(inout), allocatable, dimension(:) :: lines
    type(material), intent(inout), allocatable, dimension(:) :: materials
    type(mesh),INTENT(INOUT) :: in_mesh
    integer,INTENT(INOUT) :: total_groups
    real(dp),intent(out), allocatable, dimension(:,:) :: source_flux ! Source flux array
    integer,allocatable,dimension(:) :: boundary_tracker
    integer :: i,xlines,ylines
    integer :: j
    INTEGER :: k
    integer :: source_tracker
    integer :: status ! Checks if end of file
    integer :: total_nodes
    integer :: total_steps
    integer :: total_regions
    integer :: total_materials
    integer :: dimension
    character(80) :: file_line
    character(80) :: mid ! Material ID
    character(80) :: lid ! Line ID
    character(80) :: lid2
    character(80) :: fission_or_volumetric ! ID for a volumetric or fission source problem
    !
    ! Read through file.
    !
    call read_input_deck(filename,total_regions,total_materials,total_groups,dimension,regions,regions2,lines,materials,boundary_tracker,fission_or_volumetric)
    ! Allocate the lines, materials and the regions
    ! Also define total number of nodes
    call allocate_variables(regions,regions2,lines,filename,total_nodes,materials,in_mesh,dimension,total_regions,fission_or_volumetric,source_flux,total_groups,boundary_tracker)
  end subroutine initialise

  subroutine correct_source(source,regions,boundary_tracker)
    !
    ! Subroutine to multiply source values by its delta value for periodic BC
    !
    implicit NONE
    ! Declare calling arguments
    real(dp),INTENT(INOUT),DIMENSION(:) :: source
    type(region_1d),INTENT(IN),DIMENSION(:) :: regions
    integer,INTENT(IN),DIMENSION(:) :: boundary_tracker
    integer :: i
    integer :: region_iterator
    region_iterator = 1
    do i = 1, size(source)
      if (i==1) then
        source(1) = source(1) * ((regions(size(regions))%get_delta()+regions(1)%get_delta())/2)
      else if (i == boundary_tracker(region_iterator) .and. i /= size(source)) then
        source(i) = source(i) * ((regions(region_iterator)%get_delta()+regions(region_iterator+1)%get_delta())/2)
        region_iterator = region_iterator +1
      else
        source(i) = source(i)*regions(region_iterator)%get_delta()
      end if
    end do
  end subroutine correct_source

  subroutine read_input_deck(filename,total_regions,total_materials,total_groups,dimension,regions,regions2,lines,materials,boundary_tracker,fission_or_volumetric)
    !
    ! Subroutine to read through file.
    !
    implicit none
    ! Declare calling arguments
    character(len=*),INTENT(IN) :: filename
    integer,INTENT(INOUT) :: total_regions,total_materials,total_groups,dimension
    type(region_1d),INTENT(INOUT),allocatable,dimension(:) :: regions
    type(region_2d),INTENT(INOUT),allocatable,dimension(:) :: regions2
    type(line),INTENT(INOUT),allocatable,dimension(:) :: lines
    type(material),INTENT(INOUT),allocatable,dimension(:) :: materials
    integer :: i,j,status
    character(80) :: file_line
    integer,INTENT(INOUT),allocatable,DIMENSION(:) :: boundary_tracker
    character(80),INTENT(OUT) :: fission_or_volumetric
    open(11,file=filename,iostat=status)
    do
    	read(11,'(A)',iostat=status) file_line ! (10,'(A)') <-- the 10 indicates write to file 10, the '(A)' indicates read the full line as a string
    	if (status /= 0) exit ! exit if end of file (or fail).
        call region_set(file_line,status,total_regions,boundary_tracker)
        call material_set(file_line,status,total_materials,materials)
        call group_set(file_line,status,total_groups)
        call dimension_set(file_line,status,dimension,regions,regions2,total_regions)
        call id_set(file_line,status,dimension,lines,regions,regions2,total_regions)
        call mid_set(file_line,status,total_materials,total_groups,materials)
        call lid_set(file_line,dimension,status,lines)
        call problem_type_set(file_line,status,fission_or_volumetric)
    end do
    close(11)
end subroutine read_input_deck

subroutine allocate_variables(regions,regions2,lines,filename,total_nodes,materials,in_mesh,dimension,total_regions,fission_or_volumetric,source_flux,total_groups,boundary_tracker)
    !
    ! Subroutine to allocate values
    !
    implicit none
    ! Declare calling arguments
    character(len=*),INTENT(IN) :: filename
    character(80),INTENT(IN) :: fission_or_volumetric
    real(dp),INTENT(INOUT),ALLOCATABLE,DIMENSION(:,:) :: source_flux
    integer,INTENT(INOUT) :: total_nodes
    type(region_1d),INTENT(INOUT),allocatable,dimension(:) :: regions
    type(region_2d),INTENT(INOUT),allocatable,dimension(:) :: regions2
    type(line),INTENT(INOUT),allocatable,dimension(:) :: lines
    type(material),INTENT(INOUT),allocatable,dimension(:) :: materials
    type(mesh),INTENT(INOUT) :: in_mesh
    integer :: i,j,total_steps
    integer,INTENT(INOUT),DIMENSION(:) :: boundary_tracker
    integer,intent(in) :: dimension,total_regions,total_groups
    total_nodes=0
    do i = 1, total_regions
        if (dimension == 1) then
            call associate_1d(regions,i,lines,total_nodes,materials,total_regions,filename)
        elseif (dimension ==2) then
            call associate_2d(regions2,i,lines,materials,total_regions,filename)
        else
            stop 'Dimension must be 1 or 2'
        endif
    end do
    call in_mesh%fill_coordinates(lines)
    ! Add one to the nodes so that it is correct number of x values
    total_nodes = total_nodes +1
    ! If the source type is volumetric, allocate the source flux, otherwise don't allocate, use to check later.
    !
    ! If 2D
    !
    if (dimension ==2) then
      call source_2d(fission_or_volumetric,source_flux,in_mesh,total_groups,regions2)
    !
    ! If in 1D
    !
    else
      call source_1d(fission_or_volumetric,source_flux,total_nodes,total_groups,regions)
      ! If periodic correct source
      if (regions(1)%get_left_boundary()=='p') then
        total_nodes = total_nodes -1 ! One less node
        if (regions(size(regions))%get_steps()<2) stop 'Need 2 or more steps per region for periodic boundary'
        call regions(size(regions))%set_steps(regions(size(regions))%get_steps()-1) ! Reduce the number of steps in rightmost region by 1
        do j=1,total_groups
            do i =1,size(regions)
              total_steps = total_steps + regions(i)%get_steps()
              boundary_tracker(i) = total_steps + 1 ! Stores the boundary between regions
            end do
            call correct_source(source_flux(:,j),regions,boundary_tracker)
        end do
    end if
    endif
end subroutine allocate_variables

subroutine region_set(file_line,status,total_regions,boundary_tracker)
    !
    ! Subroutine to allocate values
    !
    implicit none
    ! Declare calling arguments
    character(80),INTENT(INOUT) :: file_line
    integer,INTENT(INOUT) :: status
    integer,intent(out) :: total_regions
    integer,intent(out),allocatable,dimension(:) :: boundary_tracker
    if (file_line == 'Regions') then
        read(11,*,iostat=status) file_line, total_regions
        allocate(boundary_tracker(1:total_regions))
    endif
end subroutine region_set

subroutine material_set(file_line,status,total_materials,materials)
    !
    ! Subroutine to allocate values
    !
    implicit none
    ! Declare calling arguments
    character(80),INTENT(INOUT) :: file_line
    integer,INTENT(INOUT) :: status
    integer,intent(out) :: total_materials
    type(material),intent(inout),allocatable,dimension(:) :: materials
    if (file_line == 'Materials') then
        read(11,*,iostat=status) file_line, total_materials
        allocate(materials(1:total_materials))
    endif
end subroutine material_set

subroutine group_set(file_line,status,total_groups)
    !
    ! Subroutine to allocate values
    !
    implicit none
    ! Declare calling arguments
    character(80),INTENT(INOUT) :: file_line
    integer,INTENT(INOUT) :: status
    integer,intent(out) :: total_groups
    if (file_line == 'Groups') read(11,*,iostat=status) file_line, total_groups
end subroutine group_set

subroutine dimension_set(file_line,status,dimension,regions,regions2,total_regions)
    !
    ! Subroutine to allocate values
    !
    implicit none
    ! Declare calling arguments
    character(80),INTENT(INOUT) :: file_line
    integer,INTENT(INOUT) :: status
    integer,intent(in) :: total_regions
    integer,intent(out) :: dimension
    type(region_1d),INTENT(INOUT),dimension(:),allocatable :: regions
    type(region_2d),INTENT(INOUT),dimension(:),allocatable :: regions2
    if(index(file_line,'Dimension:')/=0) then
        backspace(unit=11) ! Go back to the start of the line
        read(11,*,iostat=status) file_line, dimension
        if(dimension==1) allocate(regions(1:total_regions))
        if(dimension==2) allocate(regions2(1:total_regions))
    endif
end subroutine dimension_set

subroutine id_set(file_line,status,dimension,lines,regions,regions2,total_regions)
    !
    ! Subroutine to allocate values
    !
    implicit none
    ! Declare calling arguments
    character(80),INTENT(INOUT) :: file_line
    integer,INTENT(INOUT) :: status
    integer,intent(in) :: total_regions
    integer,intent(out) :: dimension
    type(region_1d),INTENT(INOUT),dimension(:) :: regions
    type(region_2d),INTENT(INOUT),dimension(:) :: regions2
    type(line),INTENT(INOUT),dimension(:),allocatable :: lines
    character(80) :: lid,lid2,mid
    integer :: i
    if (file_line == 'Region Number, Linex, Liney, Linez, Material ID') then
        if (dimension==1) allocate(lines(1:total_regions))
        do i = 1,total_regions! Loop until no more regions
          if (dimension==1) then
            read(11,*,iostat=status) file_line, lid, mid
            call regions(i)%set_line_id(lid)
            call regions(i)%set_material_id(mid)
            ! Can immediately set the line ID
            call lines(i)%set_id(lid)
          elseif(dimension==2) THEN
            read(11,*,iostat=status) file_line, lid,lid2, mid
            call regions2(i)%set_line_id(lid,1)
            call regions2(i)%set_line_id(lid2,2)
            call regions2(i)%set_material_id(mid)
          endif
            ! Can immediately set the line ID
          ! Set the line and material ID's associated with the region
          ! Add material ID to dummy list of ID's if it has not been recorded already.
        end do
    endif
end subroutine id_set

subroutine mid_set(file_line,status,total_materials,total_groups,materials)
    !
    ! Subroutine to allocate material id
    !
    implicit none
    ! Declare calling arguments
    character(80),INTENT(IN) :: file_line
    integer,INTENT(INOUT) :: status
    integer,intent(in) :: total_materials,total_groups
    type(material),INTENT(INOUT),dimension(:) :: materials
    character(80) :: mid
    integer :: i,j
    if (file_line == 'ref (mxx), sigma_a, source flux, nu sigma_f') then
        do i = 1,total_materials! Loop until no more materials
          do j=1,total_groups ! Make sure to skip over IDs of the same group
            if (j==1)then ! Only read the first group corresponding to each material so no repeats
              read(11,*,iostat=status) mid
              ! Set the material ID
              call materials(i)%set_id(mid)
            else
              read(11,*,iostat=status)
            end if
          end do
          ! Add material ID to dummy list of ID's if it has not been recorded already.
        end do
    endif
end subroutine mid_set

subroutine lid_set(file_line,dimension,status,lines)
    !
    ! Subroutine to allocate material id
    !
    implicit none
    ! Declare calling arguments
    character(80),INTENT(INOUT) :: file_line
    integer,INTENT(INOUT) :: status
    integer,intent(in) :: dimension
    type(line),INTENT(INOUT),allocatable,dimension(:) :: lines
    character(80) :: lid
    integer :: i,total_lines
    total_lines=0
    if (file_line == 'start, length, steps, ref (lxx), dim (x/r,y,z)'.and.dimension==2) then ! only do this for 2D (and later 3D)
        do i =1,2
          read(11,*,iostat=status)
        end do
        read(11,*,iostat=status) file_line, total_lines
        allocate(lines(1:total_lines))
        print*,'allocated lines'
        do i=1,total_lines
          read(11,*,iostat=status) file_line, file_line,file_line,lid,file_line
          call lines(i)%set_id(lid)
        end do
    endif
end subroutine lid_set

subroutine problem_type_set(file_line,status,fission_or_volumetric)
    !
    ! Subroutine to allocate material id
    !
    implicit none
    ! Declare calling arguments
    character(80),INTENT(IN) :: file_line
    integer,INTENT(INOUT) :: status
    character(80),intent(out) :: fission_or_volumetric
    if (file_line == "Source Type - (f)ission or (v)olumetric") read(11,*,iostat=status) fission_or_volumetric
end subroutine problem_type_set

subroutine associate_1d(regions,i,lines,total_nodes,materials,total_regions,filename)
    !
    ! Subroutine to allocate material id
    !
    implicit none
    ! Declare calling arguments
    character(len=*),INTENT(IN) :: filename
    integer,INTENT(INOUT) :: total_nodes
    type(region_1d),intent(inout),dimension(:) :: regions
    type(line),intent(inout),dimension(:) :: lines
    type(material),intent(inout),dimension(:) :: materials
    integer, intent(in) :: i,total_regions
    integer :: j
    call regions(i)%associate_line(lines(i))
    total_nodes=total_nodes+regions(i)%get_steps()! total nodes
    do j = 1, size(materials)
        ! Associate correct rgion with correct material
        if (materials(j)%get_id() == regions(i)%get_material_id()) then
        call regions(i)%associate_material(materials(j))
        print *, 'associated region',i,'with material',j
        end if
        ! Define the materials (at the last region just so it only happens once)
        if (i==total_regions) then
            call materials(j)%read_variables(filename)
        end if
    end do
end subroutine associate_1d

subroutine associate_2d(regions2,i,lines,materials,total_regions,filename)
    !
    ! Subroutine to associate region i with lines and materials
    !
    implicit none
    ! Declare calling arguments
    character(len=*),INTENT(IN) :: filename
    type(region_2d),intent(inout),dimension(:) :: regions2
    type(line),intent(inout),dimension(:) :: lines
    type(material),intent(inout),dimension(:) :: materials
    integer, intent(in) :: i,total_regions
    integer :: j
    do j = 1, size(materials)
        ! Associate correct rgion with correct material
        if (materials(j)%get_id() == regions2(i)%get_material_id()) then
            call regions2(i)%associate_material(materials(j))
            print *, 'associated region',i,'with material ',materials(j)%get_id()
        end if
        ! Define the materials (at the last region just so it only happens once)
        if (i==total_regions) then
            call materials(j)%read_variables(filename)
        end if
    end do
    do j = 1, size(lines)
        ! Define the lines (at the first region just so it only happens once)
        if (i==1) call lines(j)%read_variables(filename)
        ! Associate correct rgion with correct material
        if (lines(j)%get_id() == regions2(i)%get_line_id(1)) then
            call regions2(i)%associate_line(lines(j))
            print *, 'associated region',i,'with x line ',regions2(i)%get_line_id(1)
        elseif (lines(j)%get_id() == regions2(i)%get_line_id(2)) then
            call regions2(i)%associate_line(lines(j))
            print *, 'associated region',i,'with y line ',regions2(i)%get_line_id(2)
        end if
    end do
end subroutine associate_2d

subroutine source_2d(fission_or_volumetric,source_flux,in_mesh,total_groups,regions2)
    !
    ! Subroutine to allocate material id
    !
    implicit none
    ! Declare calling arguments
    character(80),intent(in) :: fission_or_volumetric
    real(dp),intent(out),allocatable,dimension(:,:) :: source_flux
    type(mesh),intent(in) :: in_mesh
    type(region_2d),intent(in),dimension(:) :: regions2
    integer, intent(in) :: total_groups
    integer :: i,j,k
    if (fission_or_volumetric == "v") then
        allocate(source_flux(1:in_mesh%get_x_size()*in_mesh%get_y_size(),1:total_groups))
        do k=1,total_groups
            do j=1,in_mesh%get_y_size()
                do i=1,in_mesh%get_x_size()
                    source_flux(i+((j-1)*in_mesh%get_x_size()),k)=regions2(in_mesh%r(i,j))%get_source_flux(k)
                enddo
            enddo
        enddo
    end if
end subroutine source_2d

subroutine source_1d(fission_or_volumetric,source_flux,total_nodes,total_groups,regions)
    !
    ! Subroutine to allocate material id
    !
    implicit none
    ! Declare calling arguments
    character(80),intent(in) :: fission_or_volumetric
    real(dp),intent(out),allocatable,dimension(:,:) :: source_flux
    integer, intent(in) :: total_groups,total_nodes
    integer :: source_tracker,i,j,k
    type(region_1d),intent(in),dimension(:) :: regions
    if (fission_or_volumetric == "v") then
        allocate(source_flux(1:total_nodes,1:total_groups))
        ! Loop all this for each group
        do k = 1, total_groups
            source_tracker=0 ! Reset to zero for each group
            ! Set the source flux
            do i = 1, size(regions)
                ! First region has extra node at the start
                if (i == 1 .and. i/=size(regions)) then
                    do j = 1, regions(i)%get_steps()+1
                        source_tracker=source_tracker+1
                        ! If not at boundary node
                        if (j/=regions(i)%get_steps()+1) THEN
                            source_flux(source_tracker,k) = regions(i)%get_source_flux(k)
                        ! If at boundary
                        else
                            source_flux(source_tracker,k) = ((regions(i)%get_source_flux(k)*regions(i)%get_length()/regions(i)%get_steps())&
                            +(regions(i+1)%get_source_flux(k)*regions(i+1)%get_length()/regions(i+1)%get_steps()))&
                            /((regions(i)%get_length()/regions(i)%get_steps())+(regions(i+1)%get_length()/regions(i+1)%get_steps()))
                        end if
                    end do
                ! Last region needs no boundary correction
                else if (i == size(regions)) then
                    do j = 1, regions(i)%get_steps()
                        source_tracker = source_tracker + 1
                        if (source_tracker > size(source_flux)) exit ! This is for periodic case
                        source_flux(source_tracker,k) = regions(i)%get_source_flux(k)
                    end do
                ! For internal regions
                else
                    do j = 1, regions(i)%get_steps()
                        source_tracker = source_tracker + 1
                        if (j/=regions(i)%get_steps()) THEN
                            source_flux(source_tracker,k) = regions(i)%get_source_flux(k)
                        ! If at boundary
                        else
                            source_flux(source_tracker,k) = ((regions(i)%get_source_flux(k)*regions(i)%get_length()/regions(i)%get_steps())&
                            +(regions(i+1)%get_source_flux(k)*regions(i+1)%get_length()/regions(i+1)%get_steps()))&
                            /((regions(i)%get_length()/regions(i)%get_steps())+(regions(i+1)%get_length()/regions(i+1)%get_steps()))
                        end if
                    end do
                end if
            end do
        end do
    end if
end subroutine source_1d

END module initialise_variables
