program x_data
    double precision :: phian
    integer :: i
    double precision :: x
    open(65, file ='x_phi_analytical.txt')
    do i=1, 200
      x=5e-3*(i-1)
      if (x<0.500000000000000) then
        phian = 1.996739782*10.0**(-6.0)*exp(-x*sqrt(3.0))+1.000000001+1.141544254*10.0**(-6.0)* &
        exp(sqrt(3.0)*(x+4.0))-9.693634286*10.0**(-7.0)*exp(-sqrt(3.0)*(x-6.0))+ &
        7.051925394*10.0**(-8.0)*exp((8.0+x)*sqrt(3.0))-1.105760687*10.0**(-8.0)* &
        exp(-(-10.0+x)*sqrt(3.0))-4.077336708*10.0**(-7.0)*exp(sqrt(3.0)*(x+6.0))+ &
        8.398687395*10.0**(-7.0)*exp(x*sqrt(3.0))+1.676554838*10.0**(-7.0)*exp(-(-8.0+x)* &
        sqrt(3.0))+2.713956024*10.0**(-6.0)*exp(-sqrt(3.0)*(x-4.0))-4.651050887*10.0**(-9.0)* &
        exp((10.0+x)*sqrt(3.0))-1.561653034*10.0**(-6.0)*exp(sqrt(3.0)*(x+2.0))- &
        3.712740565*10.0**(-6.0)*exp(-sqrt(3.0)*(x-2.0))
      else
        phian = .4999999966+0.1303085423e-4*exp(2.0*sqrt(3.0)*(x+1.0))+0.4611141268e-4* &
        exp(2.0*x*sqrt(3.0))-0.3107835588e-2*exp(-2.0*x*sqrt(3.0))-0.2109035493e-1 &
        *exp(-2.0*sqrt(3.0)*(x+2.0))-0.1988503424e-4*exp(-2.0*sqrt(3.0)*(x-3.0))- &
        3.518084245*10.0**(-6.0)*exp(2.0*sqrt(3.0)*(x+2.0))-0.8042027998e-2* &
        exp(-2.0*sqrt(3.0)*(x+4.0))+0.3590762399e-2*exp(2.0*sqrt(3.0)*(x-4.0))+ &
        0.2005211854e-2*exp(2.0*sqrt(3.0)*(x-2.0))-0.1422805295e-2* &
        exp(2.0*sqrt(3.0)*(x-5.0))+0.2029582814e-1*exp(-2.0*sqrt(3.0)*(x+3.0))+ &
        1.39400727*10.0**(-6.0)*exp(-2.0*sqrt(3.0)*(x-4.0))-0.5498420208e-3* &
        exp(2.0*sqrt(3.0)*(x-1.0))+2.466294477*10.0**(-7.0)*exp(2.0*sqrt(3.0)*(x+3.0))+ &
        0.736534331e-4*exp(-2.0*sqrt(3.0)*(x-2.0))-0.3731331042e-2* &
        exp(2.0*sqrt(3.0)*(x-3.0))+0.1133392598e-1*exp(-2.0*sqrt(3.0)*(x+1.0))+ &
        0.2606324798e-3*exp(-2.0*sqrt(3.0)*(x-1.0))
      endif
      write(65,*) x,phian
    enddo
    close(65)
end program x_data