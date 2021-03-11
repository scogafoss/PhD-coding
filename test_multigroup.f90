program test_multigroup
    use precision_set
    use read_gem_file
    use line_class
    use material_class
    use region_class_1d
    use error_class
    use maths_class
    use nuclear_matrix_class
    use matrix_class
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
    real(dp), allocatable, dimension(:) :: source_flux
    real(dp), allocatable, dimension(:) :: gem_x ! Values of x used by gem
    real(dp), allocatable, dimension(:) :: finite_phi ! Values of x used by gem
    REAL(dp), ALLOCATABLE, DIMENSION (:) :: gem_phi ! Gem values of phi
    real(dp), allocatable, dimension(:) :: finite_x ! All vales of x phi is defined at
    real(dp), allocatable, dimension(:) :: interpolated_phi ! stores the phi as an allocatable array because that is needed for some reason
    real(dp) :: keff
    type(solver) :: solve
    character(150) :: gem_file
    type(maths) :: maths1
    integer :: m
    integer :: k
    integer :: i
    filename = 'input_deck_multigroup.dat'
    call initialise(filename,regions,lines,materials,source_flux)
    ALLOCATE(matrix_array(1:size(regions(1)%get(absorption()))))
    do i = 1, size(matrix_array) ! Does this work???????????????
        call nuc1%discretise_regions(regions)
        call matrix_array(i)%set_variables(nuc1%get_a(),nuc1%get_b(),nuc1%get_c())
    call solve%multigroup_solver(finite_phi,)
    print *, 'Effective neutron multiplication factor:',keff
  end program test_multigroup
  