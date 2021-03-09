module initialise_variables
  use precision_set
  use line_class
  use material_class
  use region_class_1d
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
!!                                                                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE
contains
  subroutine initialise(filename,regions,lines,materials,source_flux)
    character(len=*), intent(in) :: filename
    type(region_1d), intent(inout), allocatable, dimension(:) :: regions
    type(line),  intent(inout), allocatable, dimension(:) :: lines
    type(material), intent(inout), allocatable, dimension(:) :: materials
    real(dp),intent(out), allocatable, dimension(:) :: source_flux ! Source flux array
    integer :: i
    integer :: j
    integer :: source_tracker=0
    integer :: status ! Checks if end of file
    integer :: total_nodes
    integer :: total_regions
    integer :: total_materials
    character(80) :: line
    character(80) :: mid ! Material ID
    character(80) :: lid ! Line ID
    open(11,file=filename,iostat=status)
    ! Read through file.
    do
    	read(11,'(A)',iostat=status) line ! (10,'(A)') <-- the 10 indicates write to file 10, the '(A)' indicates read the full line as a string
    	if (status /= 0) exit ! exit if end of file (or fail).
      if (line == 'Regions') then
        read(11,*,iostat=status) line, total_regions
        allocate(regions(1:total_regions))
        allocate(lines(1:total_regions))
      else if (line == 'Materials') then
        read(11,*,iostat=status) line, total_materials
        allocate(materials(1:total_materials))
      else if (line == 'Region Number, Linex, Liney, Linez, Material ID') then
        do i = 1,total_regions! Loop until no more regions
          read(11,*,iostat=status) line, lid, mid
          ! Set the line and material ID's associated with the region
          call regions(i)%set_line_id(lid)
          call regions(i)%set_material_id(mid)
          ! Can immediately set the line ID
          call lines(i)%set_id(lid)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !lines(i)%read_variables(filename) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DO THIS ONCE THE FILE HAS BEEN CLOSED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Add material ID to dummy list of ID's if it has not been recorded already.
        end do
      else if (line == 'ref (mxx), sigma_a, source flux, nu sigma_f') then
        do i = 1,total_materials! Loop until no more materials
          read(11,*,iostat=status) mid
          ! Set the material ID
          call materials(i)%set_id(mid)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !lines(i)%read_variables(filename) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! DO THIS ONCE THE FILE HAS BEEN CLOSED !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! Add material ID to dummy list of ID's if it has not been recorded already.
        end do
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
        end if
        ! Define the materials (at the last region just so it only happens once)
        if (i==size(lines)) then
          call materials(j)%read_variables(filename)
        end if
      end do
    end do
    ! Add one to the nodes so that it is correct number of x values
    total_nodes = total_nodes +1
    allocate(source_flux(1:total_nodes))
    ! Set the source flux
    do i = 1, size(regions)
      ! First region has extra node at the start
      if (i == 1) then
        do j = 1, regions(i)%get_steps()+1
          source_tracker=source_tracker+1
          ! If not at boundary node
          if (j/=regions(i)%get_steps()+1) THEN
            source_flux(source_tracker) = regions(i)%get_source_flux()
          ! If at boundary
          else
            source_flux(source_tracker) = ((regions(i)%get_source_flux()*regions(i)%get_length()/regions(i)%get_steps())&
            +(regions(i+1)%get_source_flux()*regions(i+1)%get_length()/regions(i+1)%get_steps()))&
            /((regions(i)%get_length()/regions(i)%get_steps())+(regions(i+1)%get_length()/regions(i+1)%get_steps()))
          end if
        end do
      ! Last region needs no boundary correction
    else if (i == size(regions)) then
        do j = 1, regions(i)%get_steps()
          source_tracker = source_tracker + 1
          source_flux(source_tracker) = regions(i)%get_source_flux()
        end do
      ! For internal regions
      else
        do j = 1, regions(i)%get_steps()
          source_tracker = source_tracker + 1
          if (j/=regions(i)%get_steps()) THEN
            source_flux(source_tracker) = regions(i)%get_source_flux()
          ! If at boundary
          else
            source_flux(source_tracker) = ((regions(i)%get_source_flux()*regions(i)%get_length()/regions(i)%get_steps())&
            +(regions(i+1)%get_source_flux()*regions(i+1)%get_length()/regions(i+1)%get_steps()))&
            /((regions(i)%get_length()/regions(i)%get_steps())+(regions(i+1)%get_length()/regions(i+1)%get_steps()))
          end if
        end do
      end if
    end do
  end subroutine initialise
END module initialise_variables
