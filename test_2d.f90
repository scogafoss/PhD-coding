program test_2d
    use read_gem_file
    use error_class
    use populate_compressed_class
    use populate_tridiagonal_class
    use populate_compressed_class_2d
    use initialise_variables
    use solver_class
    implicit none
    character(80) :: filename
    type(region_1d), allocatable, dimension(:) :: regions
    type(region_2d),allocatable,DIMENSION(:) :: regions2
    type(line), target, allocatable, dimension(:) :: lines
    type(material), target, allocatable, dimension(:) :: materials
    type(populate_tridiagonal) :: popt
    type(populate_compressed) :: popc
    type(populate_compressed_2d) :: popc2
    type(tridiagonal_matrix),allocatable,dimension(:) :: matrix_array
    type(compressed_matrix),allocatable,DIMENSION(:) :: c_matrix_array
    type(error) :: error1
    type(mesh) :: in_mesh
    integer :: x_tracker=1
    real(dp), allocatable, dimension(:,:) :: source_flux
    real(dp), allocatable, dimension(:) :: gem_x ! Values of x used by gem
    real(dp), allocatable, dimension(:,:) :: finite_phi ! Values of x used by gem
    REAL(dp), ALLOCATABLE, DIMENSION (:,:) :: gem_phi ! Gem values of phi
    real(dp), allocatable, dimension(:,:) :: finite_x ! All vales of x phi is defined at
    real(dp), allocatable, dimension(:,:) :: interpolated_phi ! stores the phi as an allocatable array because that is needed for some reason
    real(dp) :: keff
    type(solver) :: solve
    character(150) :: gem_file
    integer :: m
    integer :: k
    integer :: i
    integer :: groups
    real(dp),ALLOCATABLE,DIMENSION(:) :: x_coordinate
    !
    ! Define file names
    !
    filename = 'input_deck_2D.dat'
    gem_file = "/mnt/c/Users/scoga/OneDrive - Imperial College &
    London/PhD/Gem events/slab_1d_volum2.event/out/detect/slab_1d_volum2"
    !
    ! Initialise the necessary variables for problem
    !
    call initialise(filename,lines,materials,source_flux,groups,regions,regions2,in_mesh)
    if(allocated(regions2)) then
        allocate(c_matrix_array(1:groups))
    else
        if(regions(1)%get_left_boundary()=='p') then
            allocate(c_matrix_array(1:groups))
        else
            ALLOCATE(matrix_array(1:groups))
        endif
    endif
    !
    ! Populate the necessary matrix array
    !
    do i = 1, groups
        if(allocated(regions2)) then
            call popc2%populate_matrix(c_matrix_array(i),regions2,i,in_mesh)
        else
            if(regions(1)%get_left_boundary()=='p') then
                call popc%populate_matrix(c_matrix_array(i),regions,i)
            else
                call popt%populate_matrix(matrix_array(i),regions,i) ! ith group
            endif
        endif
    end do
    !
    ! Perform the multigroup solve
    !
    if(allocated(regions2)) then
        call solve%multigroup_solver(finite_phi,keff,regions2=regions2,c_matrix=c_matrix_array,source_flux=source_flux,&
        x_coordinate=x_coordinate,in_mesh=in_mesh)
    else
        if(regions(1)%get_left_boundary()=='p') then
            call solve%multigroup_solver(finite_phi,keff,regions=regions,c_matrix=c_matrix_array,source_flux=source_flux,&
            x_coordinate=x_coordinate)
        else
            call solve%multigroup_solver(finite_phi,keff,regions=regions,matrix_array=matrix_array,source_flux=source_flux,&
            x_coordinate=x_coordinate)
        endif
    endif    
    print *, 'Effective neutron multiplication factor:',keff
    !
    ! Compare to gem data
    !
    open(60, file ='x_y1_y2.txt')
    open(70, file ='xgem_ygem1_ygem2.txt')
    do i=1,size(x_coordinate)
        write(60,*) x_coordinate(i), finite_phi(i,:)
    end do
    close (60)
    call read_gem(gem_file,gem_x,gem_phi)
    do i=1,size(gem_x)
        write(70,*) gem_x(i),gem_phi(i,:)
    end do
    close (70)
    do i = 1,groups
        print *,'group',i,'L2 = ',error1%l2_from_fluxes(finite_phi(:,i),x_coordinate,gem_phi(:,i),gem_x), 'delta = ',&
        x_coordinate(2)-x_coordinate(1)
    end do
  end program test_2d
