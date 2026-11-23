module initialize_flow
    use types
    implicit none

    public ::  init_meanflow
    public ::  init_meanflow_freq


    type(meanflow2D) :: flow2D
contains

    !! @brief initilization of the 2D mean flow for the propagation
    !!
    !! @details The function either take the 3D mesh as input and interpolate in the 2D mesh for the propagation
    !! or it read an already interpolated 2D mesh and check that the dimension c orrespond to the parameters of the simulation.
    !! Note that mean flow is interpolated on a coarser grid coresponding to dout in the x direction

    !! @param [in]      kk : number of the angle that is evalutated
    !! @param [in]      simu: parameters of the global simulation
    !! @param [in]      mesh: the 2D mesh coresponding to the kkth angles
    !! @param [in]      nbout: the number of grid step bewteen each recording of the results
    !! @param [in]      nx : number of point in the x direction
    !! @param [in]      nz : number of point in the z direction
    !! @param [in]      c0: the sound speed
    !! @param[inout]    flow3D: the 3D flow before interpolation
    !! @param[inout]    flow2D: the interpolated 2D flow inside the mesh
    !! @param[inout]    dom: parameters of the 2D mesh coresponding to the kk angles
    !! @return
    subroutine init_meanflow(kk, simu, mesh, flow3D, nx, nz, c0)
        use types
        use interpo
        use read_write
        implicit none
        ! input variables
        integer, intent(in)                                 :: kk
        type(simulation), intent(in)                        :: simu
        type(mesh_2D), intent(inout)                        :: mesh
        type(meanflow3D), intent(inout)                     :: flow3D
        integer, intent(inout)                              :: nx, nz
        real, intent(in)                                    :: c0

        !local variables
        integer                                             :: ii
        real                                                :: start, finish
        real                                                :: tau_rad
        real                                                :: pi =3.141592653589793
        ! junk
        real                                                :: dx, dz
        real, allocatable, dimension(:,:)                   :: u_temp, w_temp, c_temp

        tau_rad = simu%tau(kk)*pi/180
        flow2D%ext = simu%external_flow
        flow2D%logarithmic = simu%logarithmic
        flow2D%uniform = simu%uniform
        flow2D%u0 = simu%u0 * cos(tau_rad)

        ! Case where the flow is read from external h5 file and interpolated
        print*, '** init 2D coarse mean flow '
        allocate(flow2D%u(nx,nz))
        allocate(flow2D%w(nx,nz))
        allocate(flow2D%c(nx,nz))
        ! interpolation of 3D mean flow on 2D meshgrid
        ! definition of the 2D grid on which the flow woill be interpolated
        allocate(flow2D%x(nx,nz))
        allocate(flow2D%z(nx,nz))

        ! WARNING the flow is still interpolated on a cartesian grid !!
        ! the flow will be interpolated on the curvilinear grid for the init_meanflow_freq function
        do ii = 1,nz
        flow2D%x(:,ii) = mesh%x_coarse(1:nx)
        enddo
        do ii = 1,nx
        flow2D%z(ii,:) = mesh%z_coarse(1:nz)
        enddo


        if (flow2D%ext) then
            print*,'Interpolation of external flow field flow'
            ! interpolation of the
            CALL CPU_TIME(start)

            flow3D%u_projected = cos(tau_rad)*flow3D%u + sin(tau_rad)*flow3D%v

            ! if flow read is 1D
            if (flow3D%nx==1 .and. flow3D%ny==1) then
                print*, "interpolating 1D profile"
                call  interpolate_from_1D_to_2D(flow3D%u_projected(1,1,:),flow3D%z,u_temp, &
                    & flow2D%x,flow2D%z,tau_rad)
                call  interpolate_from_1D_to_2D(flow3D%w(1,1,:),flow3D%z,w_temp, &
                    & flow2D%x,flow2D%z,tau_rad)
                call  interpolate_from_1D_to_2D(flow3D%c(1,1,:),flow3D%z,c_temp, &
                    & flow2D%x,flow2D%z,tau_rad)

            ! if flow read is 2D 
            else if (flow3D%ny==1) then 
                print*,'Interpolation of 2D flow'
                call interpolate_from_2D_to_2D(flow3D%u_projected(:,1,:),flow3D%x,flow3D%z,u_temp, &
                    & flow3D%Xs,flow2D%x,flow2D%z,tau_rad)
                call interpolate_from_2D_to_2D(flow3D%w(:,1,:),flow3D%x,flow3D%z,w_temp, &
                    & flow3D%Xs,flow2D%x,flow2D%z,tau_rad)
                call interpolate_from_2D_to_2D(flow3D%c(:,1,:),flow3D%x,flow3D%z,c_temp, &
                    & flow3D%Xs,flow2D%x,flow2D%z,tau_rad)

            ! if flow read is 3D 
            else
                print*,'Interpolation of 3D flow'
                call  interpolate_from_3D_to_2D(flow3D%u_projected,flow3D%x,flow3D%y,flow3D%z,u_temp, &
                    & flow3D%Xs,flow3D%Ys, flow2D%x,flow2D%z,tau_rad)
                call  interpolate_from_3D_to_2D(flow3D%w,flow3D%x,flow3D%y,flow3D%z,w_temp, &
                    & flow3D%Xs,flow3D%Ys, flow2D%x,flow2D%z,tau_rad)
                call  interpolate_from_3D_to_2D(flow3D%c,flow3D%x,flow3D%y,flow3D%z,c_temp, &
                    & flow3D%Xs,flow3D%Ys, flow2D%x,flow2D%z,tau_rad)
            end if 
            ! call create_2D_interpoland(u_temp,flow2D%x(:,1),flow2D%z(1,:),flow2D%uKnot,flow2D%utx,flow2D%uty)
            ! call create_2D_interpoland(w_temp,flow2D%x(:,1),flow2D%z(1,:),flow2D%wKnot,flow2D%wtx,flow2D%wty)
            ! call create_2D_interpoland(c_temp,flow2D%x(:,1),flow2D%z(1,:),flow2D%cKnot,flow2D%ctx,flow2D%cty)
            flow2D%u(:,:) = u_temp
            flow2D%w(:,:) = w_temp
            flow2D%c(:,:) = c_temp
            CALL CPU_TIME(finish)
            print*, 'interpolation time :', finish-start, 's'

        else if (flow2D%uniform) then
            print*,'Uniform flow, u0 =', simu%u0
            flow2D%u = flow2D%u0 
            flow2D%c = c0
        else if (flow2D%logarithmic) then
            print*,'Logarithmic flow profile, u0 =', simu%u0
            do ii = 1,nx
            flow2D%u(ii,:) = (flow2D%u0/0.4)*log((mesh%z_coarse+0.1)/0.1)
            enddo
            flow2D%c = c0
        else
            print*,'No flow'
            flow2D%u = 0. 
            flow2D%c = c0
        endif
    end subroutine init_meanflow


    subroutine init_meanflow_freq(pe_model,mesh, nbout, nx, nz, c0, idom)
        use types
        use interpo
        use fast_interp
        implicit none
        ! input variables
        ! integer, intent(in)                                 :: kk
        integer, intent(in) :: pe_model
        ! type(simulation), intent(in)                        :: simu
        type(mesh_2D), intent(in)                        :: mesh
        integer, intent(in)                              :: nbout, nx, nz, idom
        real, intent(in)                                    :: c0

        !local variables
        integer                                             :: ii,jj, nxout
        real                                                :: start, finish
        real                                                :: tau_rad
        real                                                :: pi=3.141592653589793
        ! junk
        real                                                :: dx, dz
        real, allocatable, dimension(:,:)                   :: u_temp, c_temp, w_temp
        real, allocatable, dimension(:,:)                   :: x_temp, z_temp

        ! tau_rad = simu%tau(kk)*pi/180


        nxout = ceiling((1.*nx)/nbout)
        ! Case where the flow is read from external h5 file and interpolated
        allocate(x_temp(nxout,nz))
        allocate(z_temp(nxout,nz))
        allocate(u_temp(nxout,nz))
        allocate(w_temp(nxout,nz))
        allocate(c_temp(nxout,nz))
        x_temp(:,:) = mesh%x_grid(1:nx:nbout, 1:nz)
        z_temp(:,:) = mesh%z_grid(1:nx:nbout, 1:nz)


        if (flow2D%ext) then

            dx  = flow2D%x(2,1) - flow2D%x(1,1)
            dz  = flow2D%z(1,2) - flow2D%z(1,1)
            ! interpolation with spline Keys 
            c_temp = c0
            call bilinear_interp(flow2D%u, mesh%x_coarse(1), mesh%z_coarse(1), dx, &
                            dz, x_temp, z_temp, u_temp, .true.)
            call bilinear_interp(flow2D%w, mesh%x_coarse(1), mesh%z_coarse(1), dx, &
                            dz, x_temp, z_temp, w_temp, .true.)
            call bilinear_interp(flow2D%c, mesh%x_coarse(1), mesh%z_coarse(1), dx, &
                            dz, x_temp, z_temp, c_temp, .true.)
            
            ! interpolation with bspline-fortran
            ! much more costly 
            ! call interpolate_2D(flow2D%uKnot,flow2D%utx,flow2D%uty,x_temp,z_temp,u_temp)
            ! call interpolate_2D(flow2D%wKnot,flow2D%wtx,flow2D%wty,x_temp,z_temp,w_temp)
            ! call interpolate_2D(flow2D%cKnot,flow2D%ctx,flow2D%cty,x_temp,z_temp,c_temp)
            CALL CPU_TIME(finish)
            ! print*, 'interpolation time :', finish -start, 's'
            if (.not. (pe_model==5)) then
                u_temp(:,:) = u_temp(:,:)   * cos(mesh%ac(idom))   & 
                    & + w_temp(:,:) * sin(mesh%ac(idom))
            endif

        else if (flow2D%uniform) then
            u_temp = flow2D%u0
            c_temp = c0

        else if (flow2D%logarithmic) then
            do jj = 1,nxout
            u_temp(jj,:) = (flow2D%u0/0.4)*log((z_temp(jj,:)+0.1)/0.1)
            enddo
            c_temp = c0

        else
            u_temp = 0. 
            c_temp = c0
        endif

        allocate(flow2D%epsilon(nxout,nz))
        allocate(flow2D%epsilon_eff(nxout,nz))
        allocate(flow2D%n(nxout,nz))
        allocate(flow2D%n_eff(nxout,nz))
        allocate(flow2D%m(nxout,nz))
        allocate(flow2D%gamma(nxout,nz))
        allocate(flow2D%tau(nxout,nz))
        allocate(flow2D%h10(nxout,nz))
        allocate(flow2D%h12(nxout,nz))
        allocate(flow2D%h20(nxout,nz))
        allocate(flow2D%h22(nxout,nz))

        flow2D%n                    = c0/c_temp
        flow2D%n_eff                = c0/(c_temp+u_temp)
        flow2D%epsilon              = ((c0/c_temp)**2) - 1
        flow2D%epsilon_eff          = ((c0/(c_temp+u_temp))**2) - 1
        flow2D%m                    = u_temp/c_temp
        flow2D%gamma                = 1/sqrt(1- (flow2D%m**2))
        flow2D%tau                  = flow2D%m * (flow2D%gamma**2) * (flow2D%n - flow2D%m)

        if (pe_model==0) then
            ! WAPE 
            flow2D%h10 =  1.
            flow2D%h12 = 1/(4*(flow2D%n*flow2D%gamma)**2)
            flow2D%h20 = flow2D%n_eff - 1
            flow2D%h22 = 1./(2.*flow2D%n) + (flow2D%n_eff - 1.)/((2*flow2D%n*flow2D%gamma)**2)
        else if (pe_model==1) then 
            ! WAPE_m
            flow2D%h10 =  1. + 0.25*(flow2D%n**2 - 1)
            flow2D%h12 = 1/(4*(flow2D%gamma**2))
            flow2D%h20 = (0.5*(flow2D%gamma**2)*(flow2D%n**2-1)) - (1. + 0.25*(flow2D%n**2 -1))*flow2D%tau
            flow2D%h22 = 0.5 - flow2D%tau/(4*(flow2D%gamma**2))
        else if (pe_model==2) then 
            ! WAPE_eff
            flow2D%h10 =  1. + 0.25*(flow2D%n_eff**2 - 1)
            flow2D%h12 = 0.25
            flow2D%h20 = 0.5*(flow2D%n_eff**2 - 1)
            flow2D%h22 = 0.5
        else if (pe_model==3) then 
            ! NAPE
            flow2D%h10 = 1.
            flow2D%h12 = 0.
            flow2D%h20 = flow2D%n_eff - 1
            flow2D%h22 = 0.5/flow2D%n
        else if (pe_model==4) then 
            ! NAPE_eff
            flow2D%h10 = 1.
            flow2D%h12 = 0.
            flow2D%h20 = 0.5*(flow2D%n_eff**2 - 1)
            flow2D%h22 = 0.5
        else if (pe_model==5) then 

        else
            print *, "the PE model is not defined"
            stop
        end if

        if (pe_model/=5) then 
            deallocate(flow2D%epsilon_eff)
        else
            deallocate(flow2D%m)
            deallocate(flow2D%h10)
            deallocate(flow2D%h12)
            deallocate(flow2D%h20)
            deallocate(flow2D%h22)
        end if 

        deallocate(flow2D%epsilon)
        deallocate(flow2D%gamma)
        deallocate(flow2D%tau)
        deallocate(flow2D%n)
        deallocate(flow2D%n_eff)
        deallocate(x_temp)
        deallocate(z_temp)
        deallocate(u_temp)
        deallocate(w_temp)
        deallocate(c_temp)
    end subroutine init_meanflow_freq

end module
