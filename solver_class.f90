MODULE solver_class
    use tridiagonal_matrix_class
    use compressed_matrix_class
    use mesh_class
    use region_class_2d
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
  !!    19/03/2021: Added subroutines for the iterations on the flux and keff  !!
  !!      to tidy up the code. Also added a function to return weighted mean.  !!
  !!    22/03/2021: Added subroutine for fixed source to tidy up. Also added   !!
  !!      function to calculate the fission and scatter source. Method of      !!
  !!      flux iteration needed changing, should iterate scatter source not    !!
  !!      changing fission source until keff is changed.                       !!
  !!    23/03/2021: Added function to calculate x_coordinate, and normalisation!!
  !!    24/03/2021: Corrected scatter_source to use correct group's flux       !!
  !!    29/04/2021: Added funciton for compressed matrix solve                 !!
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
    procedure,PUBLIC :: flux_iteration => flux_iteration_sub ! Iteration on flux within multigroup
    procedure,PUBLIC :: weighted_average => weighted_average_fn ! Calcualtes weighted average
    procedure,PUBLIC :: k_iteration => k_iteration_sub ! Performs calculation for next iteration of eigenvalue
    procedure,public :: fixed_source_iteration => fixed_source_iteration_sub ! Performs flux iteration for volumetric (fixed) source
    procedure,public :: fission_reaction_rate => fission_reaction_rate_fn ! Returns the fission source
    procedure,public :: scatter_source => scatter_source_fn ! Returns the scatter source
    procedure,public :: x_coordinates =>  x_coordinates_sub ! Returns the position of each node in the problem
    procedure,public :: normalise => normalise_fn ! Returns normalisation
    procedure,public :: allocate_1d => allocate_1d_sub ! Allocates 1D variables
    END TYPE solver
    ! Restrict access to the actual procedure names
    PRIVATE :: multigroup_solver_sub
    PRIVATE :: flux_iteration_sub
    PRIVATE :: weighted_average_fn
    PRIVATE :: k_iteration_sub
    PRIVATE :: fixed_source_iteration_sub
    PRIVATE :: fission_reaction_rate_fn
    PRIVATE :: scatter_source_fn
    PRIVATE :: x_coordinates_sub
    PRIVATE :: normalise_fn
    PRIVATE :: allocate_1d_sub

    ! Now add methods
    contains

    SUBROUTINE multigroup_solver_sub(this, phi, keff, regions,regions2, matrix_array,c_matrix,source_flux,x_coordinate,in_mesh)
        !
        ! Subroutine to perfrom power iteration (assumes no volumetric_source)
        !
        IMPLICIT NONE
        ! Declare calling arguments
        CLASS(solver),INTENT(INOUT) :: this ! Matrix object
        integer :: phi_iterations,k_iterations
        real(dp),INTENT(inout),allocatable,dimension(:,:) :: phi ! Rows are nodes, columns are groups
        real(dp),allocatable,dimension(:,:) :: phi_temp ! phi from the previous iteration
        real(dp),INTENT(inout) :: keff
        real(dp) :: keff_temp !from the previous iteration
        real(dp),allocatable,dimension(:,:) :: source
        real(dp),ALLOCATABLE,DIMENSION(:) :: source_temp
        type(mesh),INTENT(IN),OPTIONAL :: in_mesh
        type(region_1d),INTENT(in),dimension(:),OPTIONAL :: regions
        type(region_2d),INTENT(in),dimension(:),OPTIONAL :: regions2
        type(tridiagonal_matrix),INTENT(IN),dimension(:),OPTIONAL :: matrix_array ! array of tridiagonal matrices
        type(compressed_matrix),INTENT(IN),DIMENSION(:),OPTIONAL :: c_matrix ! array of compressed matrices.
        integer,dimension(:),ALLOCATABLE :: boundary_tracker ! Labels the values of i where boundaries between regions are
        real(dp) :: convergence_criterion
        integer :: region_iterator
        integer :: total_steps
        integer :: group
        integer :: groups
        real(dp) :: numerator
        real(dp) :: denominator
        real(dp) :: normalisation
        real(dp),intent(inout),allocatable,dimension(:) :: x_coordinate ! Tracks x position in the system
        real(dp),INTENT(IN),allocatable,dimension(:,:) :: source_flux
        real(dp),allocatable,dimension(:) :: f_rate ! fission reaction ratesource
        integer :: max_iteration
        ! Max iteration
        max_iteration = 200
        ! Convergence condition
        convergence_criterion = 1e-5
        ! Record where boundaries are, so can use correct values of fission and delta
        total_steps = 0
        if (present(matrix_array)) then
            call this%allocate_1d(regions,boundary_tracker,total_steps,region_iterator,groups,source,f_rate,source_temp,keff,phi,phi_temp,x_coordinate,matrix_array=matrix_array)
        elseif (present(c_matrix)) then
            if(present(regions)) call this%allocate_1d(regions,boundary_tracker,total_steps,region_iterator,groups,source,f_rate,source_temp,keff,phi,phi_temp,x_coordinate,c_matrix=c_matrix)
            if(present(regions2)) call allocate_2d(regions2,groups,in_mesh,source,f_rate,source_temp,keff,phi,phi_temp,c_matrix)
        else
            stop 'Need a matrix array'
        endif
        !
        ! If source flux is not allocated, then this is an eigenvalue problem.
        !
        if (.not.(allocated(source_flux))) then
            !
            ! Do loop to perform iteration on keff (continues until convergence of keff)
            !
            phi_iterations=0
            k_iterations = 0
            do
                if (phi_iterations>=max_iteration) stop 'Max iterations reached'
                !
                ! Do loop to perform iteration on energy groups (continues till convergence of phi)
                !
                call this%flux_iteration(phi_iterations,phi,phi_temp,keff,source,regions,matrix_array,c_matrix,boundary_tracker,&
                convergence_criterion,total_steps,groups,f_rate,regions2,in_mesh,max_iteration)
                !
                ! Now calculate next keff and store the previous
                !
                keff_temp=keff
                numerator = 0
                denominator = 0
                region_iterator=1
                ! Summed over nodes
                if(present(regions)) call this%k_iteration(phi,phi_temp,source,regions,boundary_tracker,region_iterator,groups,numerator,denominator,f_rate,total_steps,keff,x_coordinate)
                if(present(regions2)) call k_iteration_2d(phi,phi_temp,source,regions2,groups,numerator,denominator,f_rate,keff,in_mesh)
                k_iterations = k_iterations +1
                !
                ! Calculate the new keff
                !
                keff=keff_temp * (numerator/denominator)
                !
                ! If the convergence_criterion has been met, use the new keff and exit the loop
                !
                if (abs((keff-keff_temp)/keff_temp)<convergence_criterion) then
                    if(present(regions))f_rate=this%fission_reaction_rate(groups,phi,regions,boundary_tracker,total_steps)
                    if(present(regions2))f_rate=fission_reaction_rate_2d(groups,phi,regions2,in_mesh)
                    do group=1,groups
                        if(present(regions))source(:,group)=(f_rate*regions(1)%get_probability(group)/keff)+this%scatter_source(group,groups,total_steps,phi,regions,boundary_tracker)
                        if(present(regions2))source(:,group)=(f_rate*regions2(1)%get_probability(group)/keff)+scatter_source_2d(group,groups,phi,regions2,in_mesh)
                    end do
                    do group=1,groups
                        if (present(matrix_array)) phi(:,group)=matrix_array(group)%solve(source(:,group))
                        if (present(c_matrix)) phi(:,group)=c_matrix(group)%solve(source(:,group))
                    end do
                    exit
                end if
            end do
        !
        ! If source flux is allocated then no keff iteration needed
        !
        else
            phi_iterations=0
            k_iterations = 0
            !
            ! Do loop to perform iteration on energy groups (continues till convergence of phi)
            !
            do
                if(present(regions))call this%fixed_source_iteration(groups,source_flux,total_steps,boundary_tracker,phi,phi_temp,regions,matrix_array,c_matrix)
                if(present(regions2))call fixed_source_iteration_2d(in_mesh,groups,source_flux,phi,phi_temp,regions2,c_matrix)
                ! Now this iteration will continue until the fluxes converge
                phi_iterations=phi_iterations+1
                print*,'min=',(minval(abs(phi-phi_temp))/maxval(phi))
                print*,'phi=',(phi(1,1))
                if((minval(abs(phi-phi_temp))/maxval(phi))<convergence_criterion) exit
                ! This exits the do loop for flux iterations
            end do
        end if
        !
        ! Normalise the flux values
        !
        if(present(regions))normalisation = this%normalise(source_flux,f_rate,regions,boundary_tracker,total_steps,x_coordinate,keff)
        if(present(regions2))normalisation=normalise_2d(source_flux,f_rate,regions2,in_mesh)
        ! Now normalise flux
        phi = phi / normalisation
        print *, 'Number of flux iterations = ',phi_iterations
        print *, 'Number of keff iterations = ',k_iterations
        print *,'normalisation = ',normalisation,'sum(s) = ', sum(source)
    END SUBROUTINE multigroup_solver_sub

    subroutine flux_iteration_sub(this,phi_iterations,phi,phi_temp,keff,source,regions,matrix_array,c_matrix,boundary_tracker,&
        convergence_criterion,total_steps,groups,f_rate,regions2,in_mesh,max_iteration)
        !
        ! Subroutine to perfrom power iteration (assumes no volumetric_source)
        !
        IMPLICIT NONE
        ! Declare calling arguments
        CLASS(solver),INTENT(IN) :: this ! Matrix object
        integer,INTENT(INOUT) :: phi_iterations
        real(dp),INTENT(inout),allocatable,dimension(:,:) :: phi ! Rows are groups, columns are nodes
        real(dp),INTENT(INOUT),allocatable,dimension(:,:) :: phi_temp ! phi from the previous iteration
        real(dp),INTENT(in) :: keff
        real(dp),INTENT(INOUT),allocatable,dimension(:,:) :: source
        type(region_1d),INTENT(in),dimension(:),OPTIONAL :: regions
        type(region_2d),INTENT(in),dimension(:),OPTIONAL :: regions2
        type(tridiagonal_matrix),INTENT(IN),dimension(:),optional :: matrix_array ! array of tridiagonal matrices
        type(compressed_matrix),INTENT(IN),dimension(:),optional :: c_matrix
        type(mesh),INTENT(IN) :: in_mesh
        integer,intent(inout),dimension(:) :: boundary_tracker ! Labels the values of i where boundaries between regions are
        real(dp),INTENT(IN) :: convergence_criterion
        integer,INTENT(IN) :: total_steps,max_iteration
        integer :: group
        integer,INTENT(IN) :: groups
        real(dp),intent(inout),DIMENSION(:) :: f_rate ! Fission source
        real(dp),DIMENSION(size(phi(:,1))) :: source_s ! Scatter source
        source=0
        !
        ! Do loop to perform iteration on energy groups (continues till convergence of phi)
        !
        !
        ! First find the fission reaction rate for this iteration
        !
        if(present(regions))f_rate=this%fission_reaction_rate(groups,phi,regions,boundary_tracker,total_steps)
        if(present(regions2))f_rate=fission_reaction_rate_2d(groups,phi,regions2,in_mesh)
        !
        ! Now loop scatter till convergence
        !
        do
            !
            ! Now need to perform scatter calculation for each energy group and iterate
            !
            do group = 1,groups
                !
                ! Correct source flux for fission, but only if there is no volumetric source
                !
                print*,'source s before',source_s
                if(present(regions))source_s=this%scatter_source(group,groups,total_steps,phi,regions,boundary_tracker)
                if(present(regions2))source_s=scatter_source_2d(group,groups,phi,regions2,in_mesh)
                print*,'source s after',source_s
                !
                ! Now have the total source so can find the phi iteration for current group
                !
                phi_temp(:,group)=phi(:,group)
                if(present(regions))source(:,group)=source_s+(f_rate*(regions(1)%get_probability(group)/keff))
                if(present(regions2))source(:,group)=source_s+(f_rate*(regions2(1)%get_probability(group)/keff))
                if(present(matrix_array)) phi(:,group)=matrix_array(group)%solve(source(:,group))
                if(present(c_matrix)) phi(:,group)=c_matrix(group)%solve(source(:,group))
                ! This needs to be done for all of the groups, so loop here
            end do
            ! Now this iteration will continue until the fluxes converge
            phi_iterations=phi_iterations+1
            if (phi_iterations>=max_iteration) stop 'Max flux iterations reached'
            if((minval(abs(phi-phi_temp))/maxval(phi))<convergence_criterion) exit
            ! This exits the do loop for flux iterations
        end do
    end subroutine flux_iteration_sub

    real(dp) function weighted_average_fn(this,weight1,weight2,variable1,variable2)
        !
        ! function to calculate weighted average
        !
        IMPLICIT NONE
        class(solver),intent(in) :: this
        real(dp),INTENT(IN) :: weight1
        real(dp),INTENT(IN) :: weight2
        real(dp),INTENT(IN) :: variable1
        real(dp),INTENT(IN) :: variable2
        weighted_average_fn=((variable1*weight1)+(variable2*weight2))/(weight1+weight2)
    end function weighted_average_fn

    subroutine k_iteration_sub(this,phi,phi_temp,source,regions,boundary_tracker,region_iterator,groups,numerator,denominator,&
        f_rate,total_steps,keff,x_coordinate)
        !
        ! Subroutine to perfrom power iteration (assumes no volumetric_source)
        !
        IMPLICIT NONE
        ! Declare calling arguments
        CLASS(solver),INTENT(IN) :: this ! Matrix object
        real(dp),INTENT(inout),dimension(:,:) :: phi ! Rows are groups, columns are nodes
        real(dp),intent(inout),dimension(:,:) :: phi_temp ! phi from the previous iteration
        real(dp),INTENT(INOUT),dimension(:,:) :: source
        real(dp),INTENT(IN),dimension(:) :: f_rate
        real(dp),dimension(size(f_rate)) :: f_rate_new
        type(region_1d),INTENT(in),dimension(:) :: regions
        integer,intent(inout),dimension(size(regions)) :: boundary_tracker ! Labels the values of i where boundaries between regions are
        integer,INTENT(INOUT) :: region_iterator
        integer :: i
        integer :: group
        integer,INTENT(IN) :: groups
        integer,INTENT(IN) :: total_steps
        real(dp),INTENT(INOUT) :: numerator
        real(dp),INTENT(INOUT) :: denominator
        real(dp),INTENT(in) :: keff
        real(dp),INTENT(IN),DIMENSION(:) :: x_coordinate
        ! First find the new fission source produced by the flux from flux iterations
        f_rate_new=this%fission_reaction_rate(groups,phi,regions,boundary_tracker,total_steps)
        if (sum(phi)<50000.AND.sum(phi)>5000) print*,'f rate =',f_rate,'f rate new=',f_rate_new,'flux =',phi
        ! Summed over nodes
        do i = 1,size(phi(1,:))
            ! First find numerator and denominator
            ! If at boundary use average values
            if (i == boundary_tracker(region_iterator) .and. i /= size(source(:,group))) then
                ! Numerator for source
                numerator=numerator+((((regions(region_iterator)%get_delta()+regions(region_iterator+1)%get_delta())/2)*&
                x_coordinate(i)**regions(region_iterator)%get_geomtype())*(f_rate_new(i)))
                ! denominator for source
                denominator=denominator+((((regions(region_iterator)%get_delta()+regions(region_iterator+1)%get_delta())/2)*&
                x_coordinate(i)**regions(region_iterator)%get_geomtype())*(f_rate(i)))
            ! IF not at boundary
            else
                ! New fission source
                numerator=numerator+(regions(region_iterator)%get_delta()*(f_rate_new(i))*x_coordinate(i)&
                **regions(region_iterator)%get_geomtype())
                ! Previous iteration's fission source
                denominator=denominator+(regions(region_iterator)%get_delta()*(f_rate(i))*x_coordinate(i)&
                **regions(region_iterator)%get_geomtype())
            end if
        end do
        ! print*,'numerator=',numerator,'denominator=',denominator,'phi=',phi
    end subroutine k_iteration_sub

    subroutine k_iteration_2d(phi,phi_temp,source,regions2,groups,numerator,denominator,&
        f_rate,keff,in_mesh)
        !
        ! Subroutine to perfrom power iteration (assumes no volumetric_source)
        !
        IMPLICIT NONE
        ! Declare calling arguments
        real(dp),INTENT(inout),dimension(:,:) :: phi ! Rows are groups, columns are nodes
        real(dp),intent(inout),dimension(:,:) :: phi_temp ! phi from the previous iteration
        real(dp),INTENT(INOUT),dimension(:,:) :: source
        real(dp),INTENT(IN),dimension(:) :: f_rate
        real(dp),dimension(size(f_rate)) :: f_rate_new
        type(region_2d),INTENT(in),dimension(:) :: regions2
        type(mesh),INTENT(IN) :: in_mesh
        integer :: i,j,index
        integer :: group
        integer,INTENT(IN) :: groups
        real(dp),INTENT(INOUT) :: numerator
        real(dp),INTENT(INOUT) :: denominator
        real(dp),INTENT(in) :: keff
        real(dp) :: deltax,deltay
        ! First find the new fission source produced by the flux from flux iterations
        f_rate_new=fission_reaction_rate_2d(groups,phi,regions2,in_mesh)
        if (sum(phi)<50000.AND.sum(phi)>5000) print*,'f rate =',f_rate,'f rate new=',f_rate_new,'flux =',phi
        ! Summed over nodes
        do j=1,in_mesh%get_y_size()
            do i = 1,in_mesh%get_x_size()
                index=in_mesh%index(i,j)
                deltax=regions2(in_mesh%r(i,j))%get_delta(1)
                deltay=regions2(in_mesh%r(i,j))%get_delta(2)
                ! First find numerator and denominator
                ! New fission source
                numerator=numerator+(deltax*(f_rate_new(index))*deltay)
                ! Previous iteration's fission source
                denominator=denominator+(deltax*(f_rate(index))*deltay)
            end do
        enddo
        ! print*,'numerator=',numerator,'denominator=',denominator,'phi=',phi
    end subroutine k_iteration_2d

    subroutine fixed_source_iteration_sub(this,groups,source_flux,total_steps,boundary_tracker,phi,phi_temp,regions,matrix_array,&
        c_matrix)
        !
        ! Subroutine to perfrom flux iteration on fixed source
        !
        IMPLICIT NONE
        ! Declare calling arguments
        CLASS(solver),INTENT(IN) :: this ! Solver object
        real(dp),INTENT(inout),allocatable,dimension(:,:) :: phi ! Rows are groups, columns are nodes
        real(dp),INTENT(INOUT),allocatable,dimension(:,:) :: phi_temp ! phi from the previous iteration
        real(dp),allocatable,dimension(:,:) :: source
        real(dp),INTENT(IN),dimension(:,:) :: source_flux
        type(region_1d),INTENT(in),dimension(:) :: regions
        type(tridiagonal_matrix),INTENT(IN),dimension(:),OPTIONAL :: matrix_array ! array of tridiagonal matrices
        type(compressed_matrix),INTENT(IN),dimension(:),OPTIONAL :: c_matrix ! array of matrices
        integer,intent(inout),dimension(size(regions)) :: boundary_tracker ! Labels the values of i where boundaries between regions are
        integer :: region_iterator
        integer,INTENT(IN) :: total_steps
        integer :: i
        integer :: group
        integer,INTENT(IN) :: groups
        integer :: group_iterator
        !
        ! Need to perform calculation for each energy group
        !
        do group = 1,groups
            !
            ! Correct source flux for fission, but only if there is no volumetric source
            !
            !need to sum each of the possible sources of scatter
            source = source_flux ! volumetric source only
            region_iterator = 1
            ! Loop for each node
            do i = 1, total_steps +1
                ! If at boundary need to average
                if (i == boundary_tracker(region_iterator) .and. i /= size(source(:,group))) then
                    ! Loop for each group to sum the total scatter
                    do group_iterator = 1,groups ! group iterator -> g' and group -> g
                        ! If not in the same group then fine
                        if (group_iterator /= group) then
                            source(i,group)=source(i,group)+(phi(i,group_iterator)*&
                            (this%weighted_average(regions(region_iterator)%get_delta(),&
                            regions(region_iterator+1)%get_delta(),&
                            regions(region_iterator)%get_scatter(group_iterator,group),&
                            regions(region_iterator+1)%get_scatter(group_iterator,group))))
                        ! If in the same group no additional neutrons
                        else
                            source(i,group) = source(i,group)
                        end if
                    end do
                    region_iterator=region_iterator+1
                ! If not at boundary
                else
                    ! Loop for each group to sum the total scatter
                    do group_iterator = 1,groups ! group iterator -> g' and group -> g
                        ! If not in the same group then fine
                        if (group_iterator /= group) then
                            source(i,group) = source(i,group)+(phi(i,group_iterator)*&
                            ((regions(region_iterator)%get_scatter(group_iterator,group))))
                            ! If in the same group only fission contributes
                        else
                            source(i,group) = source(i,group)
                        end if
                    end do
                end if
            end do
            !
            ! Now have the total source so can find the phi iteration for current group
            !
            phi_temp(:,group)=phi(:,group)
            if(present(matrix_array)) phi(:,group)=matrix_array(group)%solve(source(:,group))
            if(present(c_matrix)) phi(:,group)=c_matrix(group)%solve(source(:,group))
            ! print*,'test S', source(:,group)
            ! This needs to be done for all of the groups, so loop here
        end do
    end subroutine fixed_source_iteration_sub

    subroutine fixed_source_iteration_2d(in_mesh,groups,source_flux,phi,phi_temp,regions2,c_matrix)
        !
        ! Subroutine to perfrom flux iteration on fixed source
        !
        IMPLICIT NONE
        ! Declare calling arguments
        class(mesh),INTENT(IN) :: in_mesh
        real(dp),INTENT(inout),allocatable,dimension(:,:) :: phi ! Rows are groups, columns are nodes
        real(dp),INTENT(INOUT),allocatable,dimension(:,:) :: phi_temp ! phi from the previous iteration
        real(dp),allocatable,dimension(:,:) :: source
        real(dp),INTENT(IN),dimension(:,:) :: source_flux
        type(region_2d),INTENT(in),dimension(:) :: regions2
        type(compressed_matrix),INTENT(IN),dimension(:),OPTIONAL :: c_matrix ! array of matrices
        integer :: i,j
        integer :: group
        integer,INTENT(IN) :: groups
        integer :: group_iterator,index
        !
        ! Need to perform calculation for each energy group
        !
        do group = 1,groups
            !
            ! Correct source flux for fission, but only if there is no volumetric source
            !
            !need to sum each of the possible sources of scatter
            source = source_flux ! volumetric source only
            ! Loop for each box
            do j = 1, in_mesh%get_y_size()
                do i = 1, in_mesh%get_x_size()
                    index = in_mesh%index(i,j)
                    ! Loop for each group to sum the total scatter
                    do group_iterator = 1,groups ! group iterator -> g' and group -> g
                        ! If not in the same group then fine
                        if (group_iterator /= group) then
                            source(index,group) = source(index,group)+(phi(index,group_iterator)*&
                            ((regions2(in_mesh%r(i,j))%get_scatter(group_iterator,group))))
                        ! If in the same group no additional in scatter
                        else
                            source(index,group) = source(index,group)
                        end if
                    end do
                enddo
            end do
            ! Now have the total source so can find the phi iteration for current group
            !
            phi_temp(:,group)=phi(:,group)
            phi(:,group)=c_matrix(group)%solve(source(:,group))
            ! print*,'test S', source(:,group)
            ! This needs to be done for all of the groups, so loop here
        end do
    end subroutine fixed_source_iteration_2d

    function fission_reaction_rate_fn(this,groups,phi,regions,boundary_tracker,total_steps) result(fission_source)
        !
        ! function to calculate fission source
        !
        IMPLICIT NONE
        class(solver),intent(in) :: this
        integer,INTENT(IN) :: groups
        integer,INTENT(IN) :: total_steps
        real(dp),INTENT(IN),DIMENSION(:,:) :: phi
        type(region_1d),INTENT(IN),DIMENSION(:) :: regions
        integer,INTENT(IN),DIMENSION(:) :: boundary_tracker
        real(dp),dimension(size(phi(:,1))) :: fission_source ! Fission source
        integer :: group_iterator
        integer :: i
        integer :: region_iterator
        !
        ! Correct source flux for fission, but only if there is no volumetric source
        !
        !need to sum each of the possible sources of fission and scatter
        fission_source=0 ! Reset each source back to zero as new flux will change fission source
        region_iterator = 1
        ! Loop for each node
        do i = 1, total_steps +1
            ! If at boundary need to average
            if (i == boundary_tracker(region_iterator) .and. i /= size(fission_source)) then
                ! Loop for each group to sum the total fission contribution and scatter
                do group_iterator = 1,groups ! group iterator -> g' and group -> g
                    fission_source(i)=fission_source(i)+(phi(i,group_iterator)*&
                    this%weighted_average(regions(region_iterator)%get_delta(),regions(region_iterator+1)%get_delta(),&
                    regions(region_iterator)%get_fission(group_iterator),&
                    regions(region_iterator+1)%get_fission(group_iterator)))
                end do
                region_iterator=region_iterator+1
            ! If not at boundary
            else
                ! Loop for each group to sum the total fission contribution and scatter
                do group_iterator = 1,groups ! group iterator -> g' and group -> g
                    fission_source(i)=fission_source(i)+((phi(i,group_iterator))*&
                    regions(region_iterator)%get_fission(group_iterator))
                end do
            end if
        end do
    end function fission_reaction_rate_fn

    function fission_reaction_rate_2d(groups,phi,regions2,in_mesh) result(fission_source)
        !
        ! function to calculate fission source in 2D
        !
        IMPLICIT NONE
        integer,INTENT(IN) :: groups
        real(dp),INTENT(IN),DIMENSION(:,:) :: phi
        type(region_2d),INTENT(IN),DIMENSION(:) :: regions2
        type(mesh),INTENT(IN) :: in_mesh
        real(dp),dimension(size(phi(:,1))) :: fission_source ! Fission source
        integer :: group_iterator
        integer :: i,j,index
        !
        ! Correct source flux for fission, but only if there is no volumetric source
        !
        !need to sum each of the possible sources of fission and scatter
        fission_source=0 ! Reset each source back to zero as new flux will change fission source
        ! Loop for each node
        do j = 1, in_mesh%get_y_size()
            do i = 1, in_mesh%get_x_size()
                index = i +((j-1)*in_mesh%get_x_size())
                ! Loop for each group to sum the total fission contribution and scatter
                do group_iterator = 1,groups ! group iterator -> g' and group -> g
                    fission_source(index)=fission_source(index)+((phi(index,group_iterator))*&
                    regions2(in_mesh%r(i,j))%get_fission(group_iterator))
                end do
            end do
        enddo
    end function fission_reaction_rate_2d

    function scatter_source_fn(this,group,groups,total_steps,phi,regions,boundary_tracker) result(scatter_source)
        !
        ! function to calculate scatter source
        !
        IMPLICIT NONE
        class(solver),intent(in) :: this
        integer,INTENT(IN) :: group
        integer,INTENT(IN) :: groups
        integer,INTENT(IN) :: total_steps
        real(dp),INTENT(IN),DIMENSION(:,:) :: phi
        type(region_1d),INTENT(IN),DIMENSION(:) :: regions
        integer,INTENT(IN),DIMENSION(:) :: boundary_tracker
        real(dp),dimension(size(phi(:,1))) :: scatter_source ! Fission source
        integer :: group_iterator
        integer :: i
        integer :: region_iterator
        !
        ! Correct source flux for fission, but only if there is no volumetric source
        !
        !need to sum each of the possible sources of fission and scatter
        scatter_source=0 ! Reset each source back to zero as new flux will change fission source
        region_iterator = 1
        ! Loop for each node
        do i = 1, total_steps +1
            ! If at boundary need to average
            if (i == boundary_tracker(region_iterator) .and. i /= size(scatter_source)) then
                ! Loop for each group to sum the total fission contribution and scatter
                do group_iterator = 1,groups ! group iterator -> g' and group -> g
                    ! If not in the same group then fine (ignores ingroup scatter)
                    if (group_iterator /= group) then
                        scatter_source(i)=scatter_source(i)+(phi(i,group_iterator)*&
                        this%weighted_average(regions(region_iterator)%get_delta(),&
                        regions(region_iterator+1)%get_delta(),regions(region_iterator)%get_scatter(group_iterator,group),&
                        regions(region_iterator+1)%get_scatter(group_iterator,group)))
                    end if
                end do
                region_iterator=region_iterator+1
            ! If not at boundary
            else
                ! Loop for each group to sum the total fission contribution and scatter
                do group_iterator = 1,groups ! group iterator -> g' and group -> g
                    ! If not in the same group then fine
                    if (group_iterator /= group) then
                        scatter_source(i)=scatter_source(i)+(phi(i,group_iterator)*&
                        regions(region_iterator)%get_scatter(group_iterator,group))
                    end if
                end do
            end if
        end do
    end function scatter_source_fn

    function scatter_source_2d(group,groups,phi,regions2,in_mesh) result(scatter_source)
        !
        ! function to calculate scatter source
        !
        IMPLICIT NONE
        integer,INTENT(IN) :: group
        integer,INTENT(IN) :: groups
        real(dp),INTENT(IN),DIMENSION(:,:) :: phi
        type(region_2d),INTENT(IN),DIMENSION(:) :: regions2
        type(mesh),INTENT(IN) :: in_mesh
        real(dp),dimension(size(phi(:,1))) :: scatter_source ! Fission source
        integer :: group_iterator
        integer :: i,j,index
        !
        ! Correct source flux for fission, but only if there is no volumetric source
        !
        !need to sum each of the possible sources of fission and scatter
        scatter_source=0 ! Reset each source back to zero as new flux will change fission source
        ! Loop for each node
        do j=1,in_mesh%get_y_size()
            do i = 1, in_mesh%get_x_size()
                index=in_mesh%index(i,j)
                ! Loop for each group to sum the total fission contribution and scatter
                do group_iterator = 1,groups ! group iterator -> g' and group -> g
                    ! If not in the same group then scatter in
                    if (group_iterator /= group) then
                        scatter_source(index)=scatter_source(index)+(phi(index,group_iterator)*&
                        regions2(in_mesh%r(i,j))%get_scatter(group_iterator,group))
                    end if
                end do
            end do
        enddo
    end function scatter_source_2d

    subroutine x_coordinates_sub(this,total_steps,regions,boundary_tracker,x_coordinate)
        !
        ! Subroutine to calculate x_coordinates
        !
        IMPLICIT NONE
        class(solver),intent(in) :: this
        integer,INTENT(IN) :: total_steps
        real(dp),INTENT(INOUT),DIMENSION(:) :: x_coordinate
        type(region_1d),INTENT(IN),DIMENSION(:) :: regions
        integer,INTENT(IN),DIMENSION(:) :: boundary_tracker
        integer :: i
        integer :: region_iterator
        region_iterator=1
        do i = 1, total_steps+1
            if (i==1) then
            x_coordinate(i) = regions(region_iterator)%get_start()
            ! At a boundnary
            else if (i == boundary_tracker(region_iterator) .and. i /= total_steps+1) then
            x_coordinate(i) = x_coordinate(i-1) + regions(region_iterator)%get_delta()
            region_iterator = region_iterator + 1
            else
            x_coordinate(i) = x_coordinate(i-1) + regions(region_iterator)%get_delta()
            end if
        end do
    end subroutine x_coordinates_sub

    real(dp) function normalise_fn(this,source_flux,f_rate,regions,boundary_tracker,total_steps,x_coordinate,keff)
        !
        ! function to calculate normalisation
        !
        IMPLICIT NONE
        class(solver),intent(in) :: this
        real(dp),intent(in),allocatable,DIMENSION(:,:) :: source_flux
        integer,INTENT(IN) :: total_steps
        real(dp),INTENT(INOUT),DIMENSION(:) :: x_coordinate
        type(region_1d),INTENT(IN),DIMENSION(:) :: regions
        integer,INTENT(IN),DIMENSION(:) :: boundary_tracker
        real(dp),INTENT(IN),DIMENSION(:) :: f_rate
        real(dp),intent(in) :: keff
        integer :: i
        integer :: region_iterator
        region_iterator=1
        ! If fission source
        if (.not.(allocated(source_flux))) then
            normalise_fn=0
            ! calculate fission source
            do i =1,total_steps ! note ther are total_steps+1 nodes so this goes to second to last node, which is fine as it integrates over steps.
                ! If at boundary
                if (i == boundary_tracker(region_iterator)) region_iterator = region_iterator + 1
                normalise_fn=normalise_fn+(0.5*(regions(region_iterator)%get_delta())*(((f_rate(i)/keff)*x_coordinate(i)&
                **regions(region_iterator)%get_geomtype())+(f_rate(i+1)/keff)*x_coordinate(i+1)&
                **regions(region_iterator)%get_geomtype()))
            end do
            ! Make correction for geometry
            if(regions(region_iterator)%get_geomtype()==1) then! Cylindrical
                normalise_fn=2*pi_dp*normalise_fn
            else if (regions(region_iterator)%get_geomtype()==2) then! spherical
                normalise_fn=4*pi_dp*normalise_fn
            end if
        ! If volumetric
        else
            normalise_fn=1
        end if
    end function normalise_fn

    subroutine allocate_1d_sub(this,regions,boundary_tracker,total_steps,region_iterator,groups,source,f_rate,source_temp,keff,phi,phi_temp,x_coordinate,matrix_array,c_matrix)
        !
        ! Subroutine to allocate variables for problem
        !
        implicit none
        class(solver),INTENT(INOUT) :: this
        type(region_1d),INTENT(IN),DIMENSION(:) :: regions
        integer,INTENT(OUT),ALLOCATABLE,DIMENSION(:) :: boundary_tracker
        integer,INTENT(INOUT) :: total_steps,region_iterator
        integer,INTENT(OUT) :: groups
        real(dp),INTENT(OUT),ALLOCATABLE,DIMENSION(:,:) :: source,phi,phi_temp
        real(dp),INTENT(OUT),ALLOCATABLE,DIMENSION(:) :: f_rate,source_temp,x_coordinate
        real(dp),INTENT(OUT) :: keff
        type(tridiagonal_matrix),INTENT(IN),DIMENSION(:),OPTIONAL :: matrix_array
        type(compressed_matrix),INTENT(IN),DIMENSION(:),OPTIONAL :: c_matrix
        !
        ! Track where region boundaries are
        !
        allocate(boundary_tracker(1:size(regions)))
        do region_iterator =1,size(regions)
            total_steps = total_steps + regions(region_iterator)%get_steps()
            boundary_tracker(region_iterator) = total_steps+1 ! tracks the x values where there is a boundary, also the edge of last boundary
        end do
        if(present(matrix_array)) then
            groups = size(matrix_array)
        else if(present(c_matrix)) then
            groups = size(c_matrix)
        else
            stop 'need a matrix array for solver'
        end if
        allocate(source(1:total_steps+1,1:groups))
        allocate(f_rate(1:total_steps+1))
        allocate(source_temp(1:total_steps+1))
        ! Initial guesses
        keff = 1
        allocate(phi(1:total_steps+1,1:groups))
        allocate(phi_temp(1:total_steps+1,1:groups))
        phi = 1
        phi_temp=1
        allocate(x_coordinate(1:total_steps+1))
        call this%x_coordinates(total_steps,regions,boundary_tracker,x_coordinate) ! Calculates x coordinates in the problem
    end subroutine allocate_1d_sub

    subroutine allocate_2d(regions2,groups,in_mesh,source,f_rate,source_temp,keff,phi,phi_temp,c_matrix)
        !
        ! Subroutine to allocate variables for problem
        !
        implicit none
        type(mesh),INTENT(IN) :: in_mesh
        type(region_2d),INTENT(IN),DIMENSION(:) :: regions2
        integer,INTENT(OUT) :: groups
        real(dp),INTENT(OUT),ALLOCATABLE,DIMENSION(:,:) :: source,phi,phi_temp
        real(dp),INTENT(OUT),ALLOCATABLE,DIMENSION(:) :: f_rate,source_temp
        real(dp),INTENT(OUT) :: keff
        type(compressed_matrix),INTENT(IN),DIMENSION(:) :: c_matrix
        groups = size(c_matrix)
        allocate(source(1:in_mesh%get_x_size()*in_mesh%get_y_size(),1:groups))
        allocate(f_rate(1:in_mesh%get_x_size()*in_mesh%get_y_size()))
        allocate(source_temp(1:in_mesh%get_x_size()*in_mesh%get_y_size()))
        ! Initial guesses
        keff = 1
        allocate(phi(1:in_mesh%get_x_size()*in_mesh%get_y_size(),1:groups))
        allocate(phi_temp(1:in_mesh%get_x_size()*in_mesh%get_y_size(),1:groups))
        phi = 1
        phi_temp=1
    end subroutine allocate_2d

    real(dp) function trapezoid_2d(f_rate,in_mesh,regions2)
        !
        ! Function to calculate 2d integral using trapezoid rule
        !
        implicit none
        real(dp),INTENT(IN),DIMENSION(:) :: f_rate
        type(region_2d),INTENT(IN),DIMENSION(:) :: regions2
        type(mesh),INTENT(IN) :: in_mesh
        integer :: i,j,index
        real(dp) :: deltax,deltay
        trapezoid_2d=0
        do j=1,in_mesh%get_y_size()
            do i=1,in_mesh%get_x_size()
                index=in_mesh%index(i,j)
                deltax=regions2(in_mesh%r(i,j))%get_delta(1)
                deltay=regions2(in_mesh%r(i,j))%get_delta(2)
                trapezoid_2d=trapezoid_2d+(f_rate(index)*deltax*deltay)
            enddo
        enddo
    end function trapezoid_2d

    real(dp) function normalise_2d(source_flux,f_rate,regions2,in_mesh)
        !
        ! function to calculate normalisation
        !
        IMPLICIT NONE
        real(dp),intent(in),allocatable,DIMENSION(:,:) :: source_flux
        type(region_2d),INTENT(IN),DIMENSION(:) :: regions2
        real(dp),INTENT(IN),DIMENSION(:) :: f_rate
        type(mesh),INTENT(IN) :: in_mesh
        ! If fission source
        if (.not.(allocated(source_flux))) then
            normalise_2d=trapezoid_2d(f_rate,in_mesh,regions2)
        ! If volumetric
        else
            normalise_2d=1
        end if
    end function normalise_2d

end module solver_class
