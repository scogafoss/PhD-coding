program test_multigroup
    use read_gem_file
    use error_class
    use nuclear_matrix_class
    use initialise_variables
    use solver_class
    implicit none
    character(80) :: filename
    type(region_1d), allocatable, dimension(:) :: regions
    type(line), target, allocatable, dimension(:) :: lines
    type(material), target, allocatable, dimension(:) :: materials
    type(nuclear_matrix) :: nuc1
    type(matrix),allocatable,dimension(:) :: matrix_array
    type(error) :: error1
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
    filename = 'input_deck_multigroup.dat'
    gem_file = "/mnt/c/Users/scoga/OneDrive - Imperial College &
    London/PhD/Gem events/slab_1d_volum2.event/out/detect/slab_1d_volum2"
    call initialise(filename,regions,lines,materials,source_flux,groups)
    ALLOCATE(matrix_array(1:groups))
    do i = 1, size(matrix_array) 
        call nuc1%discretise_regions(regions,i) ! ith group
        call matrix_array(i)%set_variables(nuc1%get_a(),nuc1%get_b(),nuc1%get_c())
    end do
    call solve%multigroup_solver(finite_phi,keff, regions, matrix_array,source_flux,x_coordinate)
    print *, 'Effective neutron multiplication factor:',keff
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
        print *,'group',i,'L2 = ',error1%l2_from_fluxes(finite_phi(:,i),x_coordinate,gem_phi(:,i),gem_x), 'delta = ',x_coordinate(2)-x_coordinate(1)
    end do
  end program test_multigroup
  