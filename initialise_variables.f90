module initialise_variables
  use region_class_1d
  use region_class_2d
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
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE
contains

  subroutine initialise(filename,lines,materials,source_flux,groups,regions,regions2)
    character(len=*), intent(in) :: filename
    type(region_1d), intent(inout), allocatable, dimension(:) :: regions
    type(region_2d), intent(inout), allocatable, dimension(:) :: regions2
    type(line),  intent(inout), allocatable, dimension(:) :: lines
    type(material), intent(inout), allocatable, dimension(:) :: materials
    integer,INTENT(INOUT) :: groups
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
    integer :: total_groups,dimension
    character(80) :: file_line
    character(80) :: mid ! Material ID
    character(80) :: lid ! Line ID
    character(80) :: lid2
    character(80) :: fission_or_volumetric ! ID for a volumetric or fission source problem
    open(11,file=filename,iostat=status)
    ! Read through file.
    do
    	read(11,'(A)',iostat=status) file_line ! (10,'(A)') <-- the 10 indicates write to file 10, the '(A)' indicates read the full line as a string
    	if (status /= 0) exit ! exit if end of file (or fail).
      if (file_line == 'Regions') then
        read(11,*,iostat=status) file_line, total_regions
        allocate(boundary_tracker(1:total_regions))
      else if (file_line == 'Materials') then
        read(11,*,iostat=status) file_line, total_materials
        allocate(materials(1:total_materials))
      else if (file_line == 'Groups') then
        read(11,*,iostat=status) file_line, total_groups
        groups = total_groups
      else if(index(file_line,'Dimension:')/=0) then
        backspace(unit=11) ! Go back to the start of the line
        read(11,*,iostat=status) file_line, dimension
        if(dimension==1) allocate(regions(1:total_regions))
        if(dimension==2) allocate(regions2(1:total_regions))
      else if (file_line == 'Region Number, Linex, Liney, Linez, Material ID') then
        do i = 1,total_regions! Loop until no more regions
          if (allocated(regions)) then
            read(11,*,iostat=status) file_line, lid, mid
            call regions(i)%set_line_id(lid)
            call regions(i)%set_material_id(mid)
            allocate(lines(1:total_regions))
            ! Can immediately set the line ID
            call lines(i)%set_id(lid)
          elseif(allocated(regions2)) THEN
            read(11,*,iostat=status) file_line, lid,lid2, mid
            call regions2(i)%set_line_id(lid,1)
            call regions2(i)%set_line_id(lid2,2)
            call regions2(i)%set_material_id(mid)
          endif
            ! Can immediately set the line ID
          ! Set the line and material ID's associated with the region
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !lines(i)%read_variables(filename) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DO THIS ONCE THE FILE HAS BEEN CLOSED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Add material ID to dummy list of ID's if it has not been recorded already.
        end do
      else if (file_line == 'ref (mxx), sigma_a, source flux, nu sigma_f') then
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
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !lines(i)%read_variables(filename) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DO THIS ONCE THE FILE HAS BEEN CLOSED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Add material ID to dummy list of ID's if it has not been recorded already.
        end do
      else if (file_line == 'start, length, steps, ref (lxx), dim (x/r,y,z)'.and.dimension==2) then ! only do this for 2D (and later 3D)
        do i =1,3
          read(11,*,iostat=status)
        end do
        read(11,*,iostat=status) file_line, j
        allocate(lines(1:j))
        do i=1,size(lines)
          read(11,*,iostat=status) file_line, file_line,file_line,lid,file_line
          call lines(i)%set_id(lid)
        end do
      else if (file_line == "Source Type - (f)ission or (v)olumetric") then
        read(11,*,iostat=status) fission_or_volumetric
      end if
    end do
    close(11)
    ! Allocate the lines, materials and the regions
    ! Also define total number of nodes
    total_nodes=0
    do i = 1, size(regions)
      call lines(i)%read_variables(filename)
      call regions(i)%associate_line(lines(i))
      total_nodes=total_nodes+regions(i)%get_steps()! total nodes
      do j = 1, size(materials)
        ! Associate correct rgion with correct material
        if (materials(j)%get_id() == regions(i)%get_material_id()) then
          call regions(i)%associate_material(materials(j))
          print *, 'associated region',i,'with material',j
        end if
        ! Define the materials (at the last region just so it only happens once)
        if (i==size(lines)) then
          call materials(j)%read_variables(filename)
        end if
      end do
    end do
    ! Add one to the nodes so that it is correct number of x values
    total_nodes = total_nodes +1
    !
    ! Correct if periodic B.C
    !
    if (regions(1)%get_left_boundary()=='p') then
      total_nodes = total_nodes -1 ! One less node
      if (regions(size(regions))%get_steps()<2) stop 'Need 2 or more steps per region for periodic boundary'
      call regions(size(regions))%set_steps(regions(size(regions))%get_steps()-1) ! Reduce the number of steps in rightmost region by 1
    end if
    ! If the source type is volumetric, allocate the source flux, otherwise don't allocate, use to check later.
    if (fission_or_volumetric == "v") then
      allocate(source_flux(1:total_nodes,1:total_groups))
    end if
    ! If the source flux has been allocated, then fill it.
    if (ALLOCATED(source_flux)) then
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
    if(regions(1)%get_left_boundary()=='p') THEN
      do j=1,total_groups
        do i =1,size(regions)
          total_steps = total_steps + regions(i)%get_steps()
          boundary_tracker(i) = total_steps + 1 ! Stores the boundary between regions
        end do
        call correct_source(source_flux(:,j),regions,boundary_tracker)
      end do
    end if
  end subroutine initialise

  subroutine correct_source(source,regions,boundary_tracker)
    !
    ! Subroutine to multiply source values by its delta value for periodic BC
    !
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

END module initialise_variables
