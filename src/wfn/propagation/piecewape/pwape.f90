module pwape
    implicit none
    public :: pwape_1angle_1freq
contains

    subroutine pwape_1angle_1freq(var0,tau,f,pe_model,source_height,imp,pml,sol,output)
        use types
        use, intrinsic :: ieee_arithmetic
        use utils
        use interpo
        use fast_interp
        use initialize_flow
        ! use interp_2D_spline
        use read_write
        use hdf5
        use splines
        use simple_wape
        use impedance_model
        use create_mesh
        implicit none
        ! input parameters
        ! type(domain), intent(inout)                             :: dom
        type(real), intent(in)                                  :: f
        type(real), intent(in)                                  :: tau
        type(impedance), intent(inout)                          :: imp
        type(variables_0), intent(inout)                        :: var0
        type(perfectly_matched_layer),intent(inout)             :: pml
        real, intent(in)                                        :: source_height
        type(output_parameters),intent(in)                      :: output
        integer, intent(in)                                     :: pe_model
        integer                                                 :: nb_receiver
        ! output
        type(solution), intent(inout)                           :: sol
        type(complex), dimension(:,:), allocatable              :: phi0
        type(complex), dimension(:), allocatable                :: phi0_temp, tx,tz
        type(complex), dimension(:,:), allocatable              :: phi_temp, phiKnot
        type(real), dimension(:,:), allocatable                 :: p_over_p0
        type(real), dimension(:), allocatable                   :: x0, z0

        ! interpolated flow
        ! generel variables
        type(frequency)                                         :: freq
        REAL                                                    :: pi = 3.141592653589793, tau_rad
        real                                                    :: ztest,a, xc, zc, ac
        integer                                                 :: nx, nz
        complex                                                 :: eye = (0.0,1.0), phi_interp
        integer                                                 :: ii,jj, kk, index, idom
        real, dimension(:,:), allocatable   :: x_coarse_grid_rot, z_coarse_grid_rot, x_rec_grid_rot, z_rec_grid_rot
        real                                                    :: abs_source_height, t

         real, parameter :: level_threshold_db = -100.0   ! tune to taste
         real            :: phi0_max_abs, phi0_level_db
         integer :: ixend


        ! time computation 
        real                                                    :: start, finish, t_tot
        real                                                    :: start_propa, finish_propa, t_propa
        real                                                    :: start_dl, finish_dl, t_dl
        real                                                    :: start_interp, finish_interp, t_interp
        real                                                    :: start_starter, finish_starter, t_starter
        real                                                    :: start_solve, finish_solve, t_solve
        real                                                    :: start_flow, finish_flow, t_flow

        t_flow = 0
        t_solve = 0
        t_interp =0 
        t_starter = 0


        ! matrix allocation 
        !----------------------------------------------------------------------
        ! set up the matrix for resolution of each sub domain
        allocate(x_coarse_grid_rot(mesh%nx_coarse, mesh%nz_coarse))
        allocate(z_coarse_grid_rot(mesh%nx_coarse, mesh%nz_coarse))
        allocate(x_rec_grid_rot(mesh%nx_coarse, output%nb_receiver))
        allocate(z_rec_grid_rot(mesh%nx_coarse, output%nb_receiver))
        
        ! set up the starter matric for each subdomain 
        ! (it will be fill with value from the previous sub domain according to the 
        ! connectivity matrix
        allocate(phi0(mesh%nzmax, mesh%ndomain))
        phi0 = 0.0
        allocate(phi0_temp(mesh%nzmax))
        allocate(x0(mesh%nzmax))
        allocate(z0(mesh%nzmax))
        ! sol%deltaL_side = 20000

        print*, 'start vector PE propagation...'
        CALL CPU_TIME(start)
        ! init frequency param
        !----------------------------------------------------------------------

        ! init impedance 
        freq%k0 = 2*pi*f/var0%c
        freq%omega  = 2*pi*f
        call compute_beta(imp%mode, f, imp, var0)

        ! source define from as height with respect to ground
        ! instead of absolute height
        abs_source_height = source_height + mesh%zc(1)

        ! Initalize the starter 
        !----------------------------------------------------------------------
        idom = 1 
        ! salomon start for a point source 
        do kk=1,mesh%nzc(idom)
            z0(kk) = compute_z_grid(1,kk,1)
        enddo

        ! phi0(:,1) = sqrt(eye*freq%k0)*(1.3717-0.3701*(freq%k0*(mesh%z-source_height))**2)* &
        !     & (exp(-((freq%k0*(mesh%z-source_height))**2)/3))

        abs_source_height = -(0.0-mesh%xc(1))*sin(mesh%ac(1)) + (mesh%source_height_abs - mesh%zc(1))*cos(mesh%ac(1))
        phi0(:,1) = sqrt(eye*freq%k0)*(1.3717-0.3701*(freq%k0*(mesh%z-abs_source_height))**2)* &
            & (exp(-((freq%k0*(mesh%z-abs_source_height))**2)/3))

        ! phi0(:,1) = sqrt(eye*freq%k0)*( & 
        !     (1.3717-0.3701*(freq%k0*(mesh%z-abs_source_height))**2)*(exp(-1/3*(freq%k0*(mesh%z-abs_source_height)**2)) + & 
        !     (1-imp%beta)/(1+imp%beta)*(1.3717-0.3701*(freq%k0*(mesh%z+abs_source_height))**2)*(exp(-1/3*(freq%k0*(mesh%z+abs_source_height))**2))))
        ! phi0(:,1) = sqrt(eye*freq%k0)*(1.3717-0.3701*(freq%k0*(z0-abs_source_height))**2)* &
        !     & (exp(-((freq%k0*(z0-abs_source_height))**2)/3))
        ! print*, phi0(:,1)
        ! loop over the domain 
        !----------------------------------------------------------------------
        do idom=1, mesh%ndomain
        ! print*, "idom:", idom
        ! create sub domain mesh 
        nx = mesh%nxc(idom)
        nz = mesh%nzc(idom)
        xc = mesh%xc(idom)
        zc = mesh%zc(idom)
        ac = mesh%ac(idom)
        allocate(mesh%x_grid(nx,nz))
        allocate(mesh%z_grid(nx,nz))
        do ii=1,nx
            do jj=1,nz
            mesh%x_grid(ii,jj) = compute_x_grid(ii,jj,idom)
            mesh%z_grid(ii,jj) = compute_z_grid(ii,jj,idom)
            enddo
        enddo
        CALL CPU_TIME(start_flow)
        ! interpolate flow inside mesh 
        call init_meanflow_freq(pe_model,mesh,output%nbout,nx,nz,var0%c,idom)
        CALL CPU_TIME(finish_flow)
        t_flow = t_flow + finish_flow - start_flow

        CALL CPU_TIME(start_solve)
        allocate(sol%phi(nx, nz))
        allocate(phi_temp(nx, nz))
        allocate(sol%p_over_p0(nx, nz))
        ! damp starter inside the PML 

        ! do ii = nz - pml%size, nz
        !      t = real(ii - (nz - pml%size), kind=8) / real(pml%size, kind=8)
        !      phi0(ii, idom) = phi0(ii, idom) * cos(0.5_8 * pi * t)**2
        ! enddo


        call propa_1angle_1freq_1domain(phi0(1:nz,idom), sol%phi, sol%p_over_p0, &
            f, freq%k0, imp%beta, var0%c,  &
            pe_model,  flow2D%h10, flow2D%h12, flow2D%h20, &
            flow2D%h22, flow2D%m, output%nbout, pml, mesh%dx, mesh%dz,  &
            mesh%x_grid, mesh%z_grid,mesh%source_height_abs)

        deallocate(flow2D%h10)
        deallocate(flow2D%h12)
        deallocate(flow2D%h20)
        deallocate(flow2D%h22)
        deallocate(flow2D%m)
        CALL CPU_TIME(finish_solve)
        t_solve = t_solve + finish_solve - start_solve

        ! fill the starter 
        CALL CPU_TIME(start_starter)
        do ii = idom+1,mesh%ndomain
        if (mesh%connectivity(idom,ii) .eqv. .true.) then
            phi0_temp = ieee_value(0.0, ieee_quiet_nan)
            ! convert absolute coordinate of domain ii 
            ! to local coordinate of domain ii-1 to perform interpolation 
            x0 = 0.0
            z0 = 0.0
            do kk=1,mesh%nzc(ii)
                x0(kk) = (compute_x_grid(1,kk,ii)-xc)*cos(ac) + (compute_z_grid(1,kk,ii)-zc)*sin(ac)
                z0(kk) = -(compute_x_grid(1,kk,ii)-xc)*sin(ac) + (compute_z_grid(1,kk,ii)-zc)*cos(ac)
            enddo
            
            ! clipping z to the ground 
            if ((z0(1)>-1e-10) .and. (z0(1)<0)) then 
                z0(1) = 0.0
            endif 
            call bilinear_interp(sol%phi, &
                0.0, 0.0, mesh%dx, mesh%dz, x0, z0, phi0_temp)

            do kk=1,mesh%nzc(ii)
            if (.not. ieee_is_nan(real(phi0_temp(kk)))) then
                phi0(kk,ii) =  phi0_temp(kk)  * exp(eye * freq%k0 * x0(kk))
                ! phi0(kk,ii) =  phi0_temp(kk) ! * exp(- eye * freq%k0 * x0(1))
            endif 
            end do 
        endif
        enddo
        CALL CPU_TIME(finish_starter)

        t_starter = t_starter + finish_starter-start_starter
        ! interpolation on the receivers 
        ! ---------------------------------------------------------------------
        CALL CPU_TIME(start_interp)
        if (output%side) then
            x_coarse_grid_rot = (mesh%x_coarse_grid-xc)*cos(ac) + (mesh%z_coarse_grid-zc)*sin(ac)
            z_coarse_grid_rot = -(mesh%x_coarse_grid-xc)*sin(ac) +(mesh%z_coarse_grid-zc)*cos(ac)

            call bilinear_interp(sol%p_over_p0(:,1:nz-pml%size), 0.0, 0.0, mesh%dx, mesh%dz, &
            x_coarse_grid_rot, z_coarse_grid_rot, sol%deltaL_side)
            

            ! call spline2D_interp(sol%p_over_p0(1:nx,1:nz,idom), mesh%x, mesh%z, &
            !     x_coarse_grid_rot, z_coarse_grid_rot, sol%deltaL_side)
        end if
        x_rec_grid_rot = (mesh%x_rec_grid-xc)*cos(ac) + (mesh%z_rec_grid-zc)*sin(ac)
        z_rec_grid_rot = -(mesh%x_rec_grid-xc)*sin(ac) +(mesh%z_rec_grid-zc)*cos(ac)
        call bilinear_interp(sol%p_over_p0, 0.0, 0.0, mesh%dx, mesh%dz, &
            x_rec_grid_rot, z_rec_grid_rot, sol%deltaL_rec)


        ! stop the computation if value at the the receiver in the last domain
        ! computed are very small
        ! phi0_max_abs = 20*log10(max(maxval(abs(sol%p_over_p0(10:nx-10, &
        !     & int(floor(20/mesh%dz))))), 1e-10))
        ! if (phi0_max_abs < -50) then
        !     print*, "  [wape] domain", idom, ": field level =", phi0_level_db, &
        !             "dB < threshold (", level_threshold_db, "dB) — stopping propagation"
        !     ! remaining phi0 columns are already 0 → bilinear_interp will return 0
        !     ixend = int(floor(mesh%x_grid(1,1)/mesh%dx_coarse))
        !     sol%deltaL_rec(ixend:mesh%nx_coarse,:) = 1e-10
        !     if (output%side) then
        !     sol%deltaL_side(ixend:mesh%nx_coarse,:) = 1e-10
        !     end if
        !     ! → post-loop 20*log10(max(...,1e-10)) clamps to -200 dB automatically
        !     deallocate(sol%p_over_p0)
        !     deallocate(sol%phi)
        !     deallocate(phi_temp)
        !     deallocate(mesh%x_grid)
        !     deallocate(mesh%z_grid)
        !     t_interp = t_interp + finish_interp- start_interp
        !     exit
        ! endif


        CALL CPU_TIME(finish_interp)
        deallocate(sol%p_over_p0)
        deallocate(sol%phi)
        deallocate(phi_temp)
        deallocate(mesh%x_grid)
        deallocate(mesh%z_grid)
        t_interp = t_interp + finish_interp- start_interp
        enddo

        if (output%side) then
        sol%deltaL_side = 20*log10(max(sol%deltaL_side, 1e-10))
        endif 
        sol%deltaL_rec = 20*log10(max(sol%deltaL_rec, 1e-10))

        deallocate(phi0)
        ! deallocate(sol%phi)
        ! deallocate(sol%p_over_p0)
        CALL CPU_TIME(finish)
        print*,'done.'
        if ((finish - start) > 0.0) then 
        print*, 'x loop time : ', finish-start, 's'
        print*, 'init flow: ', t_flow/(finish-start)*100, '%'
        print*, 'solve : ', t_solve/(finish-start)*100, '%'
        print*, 'fill starter : ', t_starter/(finish-start)*100, '%'
        print*, 'interpolation side : ', t_interp/(finish-start)*100, '%'
        else
            print*, 'init flow: N/A (total time too small to measure)'
            print*, 'solve : N/A'
            print*, 'fill starter : N/A'
            print*, 'interpolation side : N/A'
        endif 

    end subroutine
end module
