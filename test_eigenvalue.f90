program test_eigenvalue
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
  implicit none
  character(80) :: filename
  type(region_1d), allocatable, dimension(:) :: regions
  type(line), target, allocatable, dimension(:) :: lines
  type(material), target, allocatable, dimension(:) :: materials
  type(nuclear_matrix) :: nuc1
  type(matrix) :: matrix1
  type(error) :: error1
  integer :: x_tracker=1
  real(dp), allocatable, dimension(:) :: source_flux
  real(dp), allocatable, dimension(:) :: gem_x ! Values of x used by gem
  real(dp), allocatable, dimension(:) :: finite_phi ! Values of x used by gem
	REAL(dp), ALLOCATABLE, DIMENSION (:) :: gem_phi ! Gem values of phi
  real(dp), allocatable, dimension(:) :: finite_x ! All vales of x phi is defined at
  real(dp), allocatable, dimension(:) :: interpolated_phi ! stores the phi as an allocatable array because that is needed for some reason
  real(dp) :: keff
  character(150) :: gem_file
  type(maths) :: maths1
  integer :: m
  integer :: k
  integer :: i
  filename = 'input_deck_new.dat'
  call initialise(filename,regions,lines,materials,source_flux)
  call nuc1%discretise_regions(regions)
  call matrix1%set_variables(nuc1%get_a(),nuc1%get_b(),nuc1%get_c())
  call matrix1%power_iteration(finite_phi, keff, regions)
  print *, 'Effective neutron multiplication factor:',keff
  ! finite_phi = matrix1%thomas_solve(source_flux)
  !
  !--------------------------------------------------------------------------------------------------------------------------------------------
  !
  gem_file = "/mnt/c/Users/scoga/OneDrive - Imperial College &
  London/PhD/Gem events/slab_1d_volum2.event/out/detect/slab_1d_volum2"
  ! Open files to be written to
  open(60, file ='x_y.txt')
  open(70, file ='xgem_ygem.txt')
	! open(2,file='L2_example_plot.txt')

	! Read in the gem data
	call read_gem(gem_file,gem_x,gem_phi)
  allocate(finite_x(1:size(source_flux)))
	do m=1,size(regions)
    if (m==size(regions)) then
      do k=1,regions(m)%get_steps()+1
        if (x_tracker == 1) then
  		    finite_x(x_tracker) = 0
          x_tracker = x_tracker + 1
        else
          finite_x(x_tracker) = finite_x(x_tracker-1)+((regions(m)%get_length()))/regions(m)%get_steps()
          ! print *, 'delta =',((regions(m)%get_length()))/regions(m)%get_steps()
          x_tracker = x_tracker + 1
        end if
      end do
    else
      do k=1,regions(m)%get_steps()
        if (x_tracker == 1) then
  		    finite_x(x_tracker) = 0
          x_tracker = x_tracker + 1
        else
          finite_x(x_tracker) = finite_x(x_tracker-1)+((regions(m)%get_length()))/regions(m)%get_steps()
          ! print *, 'delta =',((regions(m)%get_length()))/regions(m)%get_steps()
          x_tracker = x_tracker + 1
        end if
      end do
    end if
	end do
  print *, finite_x
    ! Populate maths class to perform interpolation
  ! finite_x(size(finite_x))=2
  call maths1%set_variables(gem_x,finite_x,finite_phi)
  interpolated_phi = maths1%get_interpolation()
  do i=1,size(finite_phi)
    write(60,*) finite_x(i),finite_phi(i)
  end do
  do i=1,size(gem_x)
    write(70,*) gem_x(i),gem_phi(i)
  end do
  call error1%set_variables(interpolated_phi,gem_phi)
	PRINT *, 'delta,L2=',finite_x(2)-finite_x(1), error1%get_l2()
  ! print *, 'b=',nuc1%get_b()
  ! Deallocate vectors for next delta x
  close(60)
  close(70)
end program test_eigenvalue
