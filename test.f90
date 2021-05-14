program test_sign
    implicit none

!   integer, allocatable :: data(:)

!   data = [-3,  3,  0,  2, -3]
!   write(*,*) data

!   write(*,*) findloc(sign(1, data), 1)
!   write(*,*) findloc(sign(1, data), 1, back=.true.)
  print *, sign(-12,1)
  print *, sign(-12,0)
  print *, sign(-12,-1)

  print *, sign(-12.,1.)
  print *, sign(-12.,0.)
  print *, sign(-12.,-1.)
  end program test_sign