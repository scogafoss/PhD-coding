module ShampGordon_wrapper
    use ShampGordon_oli
    ! use ode_store_class
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Filename: ShampGordon_wrapper.f90                                        !!
  !!                                                                           !!
  !!  Dependant files: shampGordon_oli                                         !!
  !!                                                                           !!
  !!  Author: Oliver Conway                             Start date: 20/07/2021 !!
  !!                                                                           !!
  !!  Purpose: wrapper for ShampGordon ODE solver                              !!
  !!                                                                           !!
  !!  Revisions:                                                               !!
  !!                                                                           !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IMPLICIT NONE
  contains
  
    subroutine wrapper(sub,neqn, y, t, tout)
      !
      ! Subroutine for wrapper of ShampGordon ODE solver
      !
      implicit none
      external :: sub
      integer,INTENT(IN) :: neqn
      integer :: iflag
      real(kind=8),INTENT(IN) :: t,tout
      real(kind=8) :: relerr,abserr
      real(kind=8),ALLOCATABLE,DIMENSION(:) :: work
      integer(kind=4),ALLOCATABLE,DIMENSION(:) :: iwork
      real(kind=8),INTENT(INOUT),DIMENSION(:) :: y
      ! neqn = size(y)
      print*,'here'!,size(work)
      iflag=1
      relerr=1.0e-8_8
      abserr=1.0e-8_8
      ALLOCATE(work(1:100+21*NEQN),iwork(1:5))
      call shampODE ( sub, neqn, y, t, tout, relerr, abserr, iflag, work, iwork )
      DEALLOCATE(work,iwork)
    end subroutine wrapper
  
  END module ShampGordon_wrapper
  