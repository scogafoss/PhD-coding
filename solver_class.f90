MODULE solver_class
    use precision_set
    use material_class
    use line_class
    use region_class_1d
    use matrix_class
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!  Filename: solver_class.f90                                               !!
  !!                                                                           !!
  !!  Dependant files: precision_set.f90, region_class_1d.f90                  !!
  !!                                                                           !!
  !!  Author: Oliver Conway                             Start date: 16/03/2021 !!
  !!                                                                           !!
  !!  Purpose: Class to perform multigroup solve                               !!
  !!                                                                           !!
  !!  Revisions:                                                               !!
  !!                                                                           !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    IMPLICIT NONE
    ! Type definition
    TYPE,PUBLIC :: solver ! This will be the name we instantiate
    ! Instance variables.
    ! None needed
    CONTAINS
    ! Bound procedures
    procedure,PUBLIC :: multigroup_solver => multigroup_solver_sub ! Solves multigroup problem
    END TYPE solver
    ! Restrict access to the actual procedure names
    PRIVATE :: multigroup_solver_sub
    ! Now add methods
    contains
    SUBROUTINE multigroup_solver_sub(this, phi, keff, regions, matrix_array,source_flux,x_coordinate)
        !
        ! Subroutine to perfrom power iteration (assumes no volumetric_source)
        !
        IMPLICIT NONE
        ! Declare calling arguments
        CLASS(solver) :: this ! Matrix object
        integer :: phi_iterations,k_iterations
        real(dp),INTENT(inout),allocatable,dimension(:,:) :: phi ! Rows are groups, columns are nodes
        real(dp),allocatable,dimension(:,:) :: phi_temp ! phi from the previous iteration
        real(dp),INTENT(inout) :: keff
        real(dp) :: keff_temp !from the previous iteration
        real(dp),allocatable,dimension(:,:) :: source
        real(dp),ALLOCATABLE,DIMENSION(:) :: source_temp
        type(region_1d),INTENT(in),allocatable,dimension(:) :: regions
        type(matrix),INTENT(IN),dimension(:) :: matrix_array ! array of tridiagonal matrices
        integer,dimension(size(regions)) :: boundary_tracker ! Labels the values of i where boundaries between regions are
        real(dp) :: convergence_criterion
        integer :: region_iterator
        integer :: total_steps
        integer :: i
        integer :: j
        integer :: group
        integer :: groups
        integer :: group_iterator
        real(dp) :: numerator
        real(dp) :: denominator
        real(dp) :: normalisation
        real(dp),intent(inout),allocatable,dimension(:) :: x_coordinate ! Tracks x position in the system
        real(dp),INTENT(IN),allocatable,dimension(:,:) :: source_flux
        ! Convergence condition
        convergence_criterion = 1e-5
        ! Record where boundaries are, so can use correct values of fission and delta
        total_steps = 0
        do region_iterator =1,size(regions)
          total_steps = total_steps + regions(region_iterator)%get_steps()
          boundary_tracker(region_iterator) = total_steps+1 ! tracks the x values where there is a boundary, also the last boundary
        end do
        allocate(source(1:size(matrix_array),1:total_steps+1))
        allocate(source_temp(1:total_steps+1))
        groups = size(matrix_array)
        ! Initial guesses
        keff = 1
        allocate(phi(1:size(matrix_array),1:total_steps+1))
        allocate(phi_temp(1:size(matrix_array),1:total_steps+1))
        phi = 1
        phi_temp=1
        allocate(x_coordinate(1:total_steps+1))
        !
        ! If source flux is not allocated this is an eigenvalue problem
        !
        if (.not.(allocated(source_flux))) then
            !
            ! Do loop to perform iteration on keff (continues until convergence of keff)
            !
            phi_iterations=0
            k_iterations = 0
            do
                !
                ! Do loop to perform iteration on energy groups (continues till convergence of phi)
                !
                do
                    !
                    ! Need to perform calculation for each energy group
                    !
                    do group = 1,groups
                        !
                        ! Correct source flux for fission, but only if there is no volumetric source
                        !                    
                        !need to sum each of the possible sources of fission and scatter
                        source(group,:)=0
                        region_iterator = 1
                        ! Loop for each node
                        do i = 1, total_steps +1
                            ! If at boundary need to average
                            if (i == boundary_tracker(region_iterator) .and. i /= size(source(group,:))) then
                                ! Loop for each group to sum the total fission contribution and scatter
                                do group_iterator = 1,groups ! group iterator -> g' and group -> g
                                    ! If not in the same group then fine
                                    if (group_iterator /= group) then
                                        source(group,i)=source(group,i)+((phi_temp(group_iterator,i)/&
                                        (regions(region_iterator)%get_delta()+regions(region_iterator+1)%get_delta()))*&
                                        (((regions(region_iterator)%get_probability(group)/keff)*&
                                        ((regions(region_iterator)%get_fission(group)*regions(region_iterator)%get_delta())+&
                                        (regions(region_iterator+1)%get_fission(group)*&
                                        regions(region_iterator+1)%get_delta())))+&
                                        ((regions(region_iterator)%get_scatter(group_iterator,group)*&
                                        regions(region_iterator)%get_delta())+&
                                        (regions(region_iterator+1)%get_scatter(group_iterator,group)*&
                                        regions(region_iterator+1)%get_delta()))))
                                    ! If in the same group only fission contributes
                                    else
                                        source(group,i)=source(group,i)+((phi_temp(group_iterator,i)/&
                                        (regions(region_iterator)%get_delta()+&
                                        regions(region_iterator+1)%get_delta()))*&
                                        (regions(region_iterator)%get_probability(group)/keff)*&
                                        ((regions(region_iterator)%get_fission(group)*regions(region_iterator)%get_delta())+&
                                        (regions(region_iterator+1)%get_fission(group)*&
                                        regions(region_iterator+1)%get_delta())))
                                    end if
                                end do
                                region_iterator=region_iterator+1
                            ! If not at boundary
                            else
                                ! Loop for each group to sum the total fission contribution and scatter
                                do group_iterator = 1,groups ! group iterator -> g' and group -> g
                                    ! If not in the same group then fine
                                    if (group_iterator /= group) then
                                        source(group,i) = source(group,i)+(phi_temp(group_iterator,i)*&
                                        (((regions(region_iterator)%get_probability(group)/keff)*&
                                        regions(region_iterator)%get_fission(group))+&
                                        (regions(region_iterator)%get_scatter(group_iterator,group))))
                                        ! If in the same group only fission contributes
                                    else
                                        source(group,i) = source(group,i)+(phi_temp(group_iterator,i)*&
                                        ((regions(region_iterator)%get_probability(group)/keff)*&
                                        regions(region_iterator)%get_fission(group)))
                                    end if
                                end do
                            end if
                        end do
                        !
                        ! Now have the total source so can find the phi iteration for current group
                        !
                        phi_temp(group,:)=phi(group,:)
                        do i=1,size(source_temp)
                            source_temp(i)=source(group,i)
                        end do
                        phi(group,:)=matrix_array(group)%thomas_solve(source_temp)
                        ! This needs to be done for all of the groups, so loop here
                    end do
                    ! Now this iteration will continue until the fluxes converge
                    phi_iterations=phi_iterations+1
                    ! if(abs(sum(phi(1,:))-sum(phi_temp(1,:)))<convergence_criterion) exit
                    if((minval(abs(phi-phi_temp))/maxval(phi))<convergence_criterion) exit
                    ! print *,'phi ',phi(1,:)
                    ! This exits the do loop for flux iterations
                end do
                !
                ! Now calculate next keff and store the previous
                !
                keff_temp=keff
                numerator = 0
                denominator = 0
                region_iterator=1
                ! Summed over nodes
                do i = 1,size(phi(1,:))
                    ! Also each node summed over each group
                    do group = 1,groups
                        ! First find numerator and denominator
                        ! If at boundary use average values
                        if (i == boundary_tracker(region_iterator) .and. i /= size(source(group,:))) then
                            ! Numerator for source
                            numerator = numerator + (((((regions(region_iterator)%get_fission(group)*&
                            (regions(region_iterator)%get_length()/regions(region_iterator)%get_steps()))+&
                            (regions(region_iterator+1)%get_fission(group)*(regions(region_iterator+1)%get_length()/&
                            regions(region_iterator+1)%get_steps())))/((regions(region_iterator)%get_length()/&
                            regions(region_iterator)%get_steps())+(regions(region_iterator+1)%get_length()/&
                            regions(region_iterator+1)%get_steps())))*phi(group,i))*(((regions(region_iterator)%get_length()/&
                            regions(region_iterator)%get_steps())+(regions(region_iterator+1)%get_length()/&
                            regions(region_iterator+1)%get_steps()))/2))
                            ! denominator for source
                            denominator = denominator + (((((regions(region_iterator)%get_fission(group)*&
                            (regions(region_iterator)%get_length()/regions(region_iterator)%get_steps()))+&
                            (regions(region_iterator+1)%get_fission(group)*(regions(region_iterator+1)%get_length()/&
                            regions(region_iterator+1)%get_steps())))/((regions(region_iterator)%get_length()/&
                            regions(region_iterator)%get_steps())+(regions(region_iterator+1)%get_length()/&
                            regions(region_iterator+1)%get_steps())))*phi_temp(group,i))*&
                            (((regions(region_iterator)%get_length()/regions(region_iterator)%get_steps())&
                            +(regions(region_iterator+1)%get_length()/regions(region_iterator+1)%get_steps()))/2))
                            region_iterator=region_iterator+1
                        ! IF not at boundary
                        else
                            numerator = numerator + (regions(region_iterator)%get_fission(group)*phi(group,i)*&
                            (regions(region_iterator)%get_length()/regions(region_iterator)%get_steps()))
                            denominator = denominator + (regions(region_iterator)%get_fission(group)*phi_temp(group,i)*&
                            (regions(region_iterator)%get_length()/regions(region_iterator)%get_steps()))
                        end if
                    end do
                end do
                k_iterations = k_iterations +1
                ! Calculate the new keff
                keff=keff_temp * (numerator/denominator)
                ! If the convergence_criterion has been met, exit the loop
                ! print *,'test converge k1-k0/k0 = ',abs((keff-keff_temp)/keff_temp)
                if (abs((keff-keff_temp)/keff_temp)<convergence_criterion) exit
            end do
        !
        ! If source flux is allocated then no keff iteration needed
        !
        else
            phi_iterations=0
            k_iterations = 0
            source = source_flux
            !
            ! Do loop to perform iteration on energy groups (continues till convergence of phi)
            !
            do
                !
                ! Need to perform calculation for each energy group
                !
                do group = 1,groups
                    !
                    ! Correct source flux for fission, but only if there is no volumetric source
                    !                    
                    !need to sum each of the possible sources of scatter
                    region_iterator = 1
                    ! Loop for each node
                    do i = 1, total_steps +1
                        ! If at boundary need to average
                        if (i == boundary_tracker(region_iterator) .and. i /= size(source(group,:))) then
                            ! Loop for each group to sum the total scatter
                            do group_iterator = 1,groups ! group iterator -> g' and group -> g
                                ! If not in the same group then fine
                                if (group_iterator /= group) then
                                    source(group,i)=source(group,i)+((phi_temp(group_iterator,i)/&
                                    (regions(region_iterator)%get_delta()+&
                                    regions(region_iterator+1)%get_delta()))*&
                                    ((regions(region_iterator)%get_scatter(group_iterator,group)*&
                                    regions(region_iterator)%get_delta())+&
                                    (regions(region_iterator+1)%get_scatter(group_iterator,group)*&
                                    regions(region_iterator+1)%get_delta())))
                                ! If in the same group no additional neutrons
                                else
                                    source(group,i) = source(group,i)
                                end if
                            end do
                            region_iterator=region_iterator+1
                        ! If not at boundary
                        else
                            ! Loop for each group to sum the total scatter
                            do group_iterator = 1,groups ! group iterator -> g' and group -> g
                                ! If not in the same group then fine
                                if (group_iterator /= group) then
                                    ! if(group==2 .and. i==1) then
                                    !     print*,'before: source(2,1)=',source(2,1),'phi_temp(1,1)=',phi_temp(1,1),'regions(1)%get_scatter(1,2)=',regions(1)%get_scatter(1,2)
                                    !     print*,'(in region:',region_iterator,')'
                                    ! end if
                                    source(group,i) = source(group,i)+(phi_temp(group_iterator,i)*&
                                    ((regions(region_iterator)%get_scatter(group_iterator,group))))
                                    ! if(group==2 .and. i==1) then
                                    !     print*,'after: source(2,1)=',source(2,1),'phi_temp(1,1)=',phi_temp(1,1),'regions(1)%get_scatter(1,2)=',regions(1)%get_scatter(1,2)
                                    !     print*,'(in region:',region_iterator,')'
                                    ! end if
                                    ! If in the same group only fission contributes
                                else
                                    source(group,i) = source(group,i)
                                end if
                            end do
                        end if
                    end do
                    !
                    ! Now have the total source so can find the phi iteration for current group
                    !
                    phi_temp(group,:)=phi(group,:)
                    do i=1,size(source_temp)
                        source_temp(i)=source(group,i)
                    end do
                    phi(group,:)=matrix_array(group)%thomas_solve(source_temp)
                    ! This needs to be done for all of the groups, so loop here
                end do
                ! Now this iteration will continue until the fluxes converge
                phi_iterations=phi_iterations+1
                if((minval(abs(phi-phi_temp))/maxval(phi))<convergence_criterion) exit
                ! This exits the do loop for flux iterations
            end do
        end if
        !
        ! Normalise the flux values
        !
        ! First track the x coordinate for each region
        region_iterator=1
        do i = 1, total_steps+1
            if (i==1) then
            x_coordinate(i) = regions(region_iterator)%get_start()
            ! At a boundnary
            else if (i == boundary_tracker(region_iterator) .and. i /= total_steps+1) then
            x_coordinate(i) = x_coordinate(i-1) + regions(region_iterator)%get_length()/regions(region_iterator)%get_steps()
            region_iterator = region_iterator + 1
            else
            x_coordinate(i) = x_coordinate(i-1) + regions(region_iterator)%get_length()/regions(region_iterator)%get_steps()
            end if
        end do
        ! Now calculate the normalisation
        ! If fission source
        if (.not.(allocated(source_flux))) then
            do i =1,total_steps ! note ther are total_steps+1 nodes so this goes to second to last node, which is fine as it integrates over steps.
                do group = 1, groups ! also sum over the total groups
                    ! If at boundary
                    if (i == boundary_tracker(region_iterator) .and. i /= total_steps+1) then
                    region_iterator = region_iterator + 1
                    normalisation=normalisation+(0.5*(regions(region_iterator)%get_length()/regions(region_iterator)%get_steps())*&
                    ((source(group,i)*(x_coordinate(i)**regions(region_iterator)%get_geomtype()))+source(group,i+1)*&
                    (x_coordinate(i+1)**regions(region_iterator)%get_geomtype())))
                    else
                    normalisation=normalisation+(0.5*(regions(region_iterator)%get_length()/regions(region_iterator)%get_steps())*&
                    ((source(group,i)*x_coordinate(i)**regions(region_iterator)%get_geomtype())+source(group,i+1)*x_coordinate(i+1)&
                    **regions(region_iterator)%get_geomtype()))
                    end if
                end do
            end do
            ! Make correction for geometry
            if(regions(region_iterator)%get_geomtype()==1) then! Cylindrical
                normalisation=2*pi_dp*normalisation
            else if (regions(region_iterator)%get_geomtype()==2) then! spherical
                normalisation=4*pi_dp*normalisation
            end if
        ! If volumetric
        else
            normalisation=1
        end if
        ! Now normalise flux
        phi = phi / normalisation
        print *, 'Number of flux iterations = ',phi_iterations
        print *, 'Number of keff iterations = ',k_iterations
        print *,'normalisation = ',normalisation,'sum(s) = ', sum(source)
      END SUBROUTINE multigroup_solver_sub
end module solver_class