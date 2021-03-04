module test_compile_module

implicit none

contains

	subroutine print_hello_world()
		write(*,*) 'Hello World!'
	end subroutine

end module
