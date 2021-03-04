program test_classes
  use precision_set
  use timer_class
  use read_gem_file
  use error_class
  use maths_class
  use line_class
  use material_class
  use nuclear_matrix_class
  use matrix_class
  implicit none
  integer :: steps
  real(dp), allocatable, dimension(:) :: source_flux
  type(material) :: m1
  type(line) :: l1
  type(matrix) :: matrix1
  type(nuclear_matrix) :: nuc1
  type(maths) :: maths1
  type(error) :: error1
  TYPE(timer) :: t ! timer object
  character(len = 150) :: gem_file ! Gem file to open.
  real(dp), allocatable, dimension(:) :: gem_x ! Values of x used by gem
  real(dp), allocatable, dimension(:) :: finite_phi ! Values of x used by gem
	REAL(dp), ALLOCATABLE, DIMENSION (:) :: gem_phi ! Gem values of phi
  real(dp), allocatable, dimension(:) :: finite_x ! All vales of x phi is defined at
  real(dp), allocatable, dimension(:) :: interpolated_phi ! stores the phi as an allocatable array because that is needed for some reason
  integer :: m
  integer :: i
  call m1%read_variables('input_deck_new.dat')
  call l1%read_variables('input_deck_new.dat')
  gem_file = "/mnt/c/Users/scoga/OneDrive - Imperial College &
  London/PhD/Gem events/slab_1d_volum2.event/out/detect/slab_1d_volum2"
  ! Open files to be written to
	open(50, file = 'l2_data_gem.txt')  !write each line in file
  open(60, file ='x_y.txt')
  open(70, file ='xgem_ygem.txt')
	! open(2,file='L2_example_plot.txt')

	CALL t%start_timer() !start the timer

	! Read in the gem data
	call read_gem(gem_file,gem_x,gem_phi)
  print*, gem_phi
	! Loop over 1000 different x interval sizes and write to l2 error file.
	DO steps = 10 , 1000
    ! First define the source flux vector
    allocate(source_flux(1:steps+1))
    source_flux = m1%get_source_flux()
    ! Now define the objects for this iteration using steps and source flux.
    call m1%set_variables(m1%get_left_boundary(),m1%get_right_boundary(),&
    m1%get_left_albedo(),m1%get_right_albedo(),m1%get_surface_source(),&
    source_flux(1),m1%get_absorption())
    call l1%set_variables(l1%get_start(),l1%get_length(),steps,l1%get_geomtype())
    call nuc1%discretise(l1,m1)
    call matrix1%set_variables(nuc1%get_a(),nuc1%get_b(),nuc1%get_c())
		! There are steps + 1 values of x (need a value at either boundary even when one step)
    ! Set the x values
    allocate(finite_phi(1:steps+1))
    allocate(finite_x(1:steps+1))
    finite_phi = matrix1%thomas_solve(source_flux)
    if (steps == 12) THEN
      print *, nuc1%get_a()
    end if
		do m=1,size(finite_phi)
			finite_x(m) = l1%get_start() + ((l1%get_length())*(m-1))/steps
		end do
    ! Populate maths class to perform interpolation
    call maths1%set_variables(gem_x,finite_x,finite_phi)
    !
    !
    ! Make sure to update this
		! if (choice == 4) then ! user chose surface - albedo
		! 	source_flux(1) = source_flux(1) - (4*surface*(-steps/length)) !need a correction to the first element of source term
		! end if
    !
    !
		! Make example plot
		! IF (steps == 500) THEN !make a plot when there are 500 nodes
		! 	DO m = 1, size(gem_x)
		! 		write(2,*)  gem_x(m),computational_phi_interpolate(m) , gem_phi(m)
		! 	END DO
		! 	!WRITE(*,*)'phi_c = ',computational_phi,'phi_a = ',analytical_phi,'x = ',x_test
		! END IF
    ! Populate error class to calculate l2
    interpolated_phi = maths1%get_interpolation()
    if (steps == 1000) then
      do i=1,size(finite_phi)
        write(60,*) finite_x(i),finite_phi(i)
      end do
      do i=1,size(gem_x)
        write(70,*) gem_x(i),gem_phi(i)
      end do
    end if
		call error1%set_variables(interpolated_phi,gem_phi)
    ! Write to file
		write(50,*) l1%get_length()/steps, error1%get_l2()
    ! Deallocate vectors for next delta x
		deallocate(source_flux)
		deallocate(finite_phi)
    deallocate(finite_x)
    deallocate(interpolated_phi)
	END DO
	close(50)
  close(60)
  close(70)
	! close(2)
	WRITE (*,*) 'Time =', t%elapsed_time(), ' s'
  ! print *, 'The flux is calculated as:',matrix1%thomas_solve(source_flux)
end program test_classes
