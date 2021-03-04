program test_classes_single
  use precision_set
  use error_class
  use maths_class
  use line_class
  use material_class
  use nuclear_matrix_class
  use matrix_class
  implicit none
  real(dp), allocatable, dimension(:) :: source_flux
  type(material) :: m1
  type(line) :: l1
  type(matrix) :: matrix1
  type(nuclear_matrix) :: nuc1
  call m1%read_variables('input_deck_new.dat')
  call l1%read_variables('input_deck_new.dat')
  call nuc1%discretise(l1,m1)
  call matrix1%set_variables(nuc1%get_a(),nuc1%get_b(),nuc1%get_c())
  allocate(source_flux(1:l1%get_steps()+1))
  source_flux = m1%get_source_flux()
  print *, 'The flux is calculated as:',matrix1%thomas_solve(source_flux)
  print *, 'a:',matrix1%get_a(),'b:',matrix1%get_b(),'c:',matrix1%get_c()
  Print *, 'S:',source_flux
end program test_classes_single
