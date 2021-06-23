program test2
    use ode_class
    implicit none
    type(ode) :: ode1,ode2,ode3
    call ode1%explicit_euler('radiation',2500.0_dp,2.0_dp,6)
    call ode2%implicit_euler('radiation',2500.0_dp,2.0_dp,6)
    call ode3%runge_kutta('radiation',2500.0_dp,2.0_dp,6)
    print*,ode3%get_solution()
end program test2