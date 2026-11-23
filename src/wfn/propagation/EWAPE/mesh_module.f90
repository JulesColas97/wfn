module mesh_mod
    use kind_mod, only: wp
    use constants, only: PI
    use input_params, only: inputs_t
    use interp_spline_mod, only: interpolate_s
    implicit none
    private :: eval_sinusoidal_surface
    public :: mesh_t,init_coarse_mesh, clean_coarse_mesh,  &
        init_fine_mesh, clean_fine_mesh, mesh, &
        interface, read_interface

    type mesh_t
        integer  :: nx, nz, nx_coarse, nz_coarse
        real(wp) :: dx, dz, dx_coarse, dz_coarse
        real(wp) :: Lx, Lz, Lx_coarse, Lz_coarse
        real(wp) :: x0, z0
        real(wp) :: source_height, source_height_abs

        real(wp), allocatable, dimension(:) :: x, x_coarse
        real(wp), allocatable, dimension(:) :: z, z_coarse
        logical, allocatable, dimension(:, :) :: water_coarse, water_fine
        real(wp), allocatable, dimension(:) :: xi_coarse, xi,xi_r, xi_rr
        integer, allocatable :: js(:), js_coarse(:)
        integer :: nbout
        integer :: nb_bot, jb 
    end type 

    type interface_t 
        integer :: nx, ny
        real(wp), allocatable :: s_2D(:,:), x(:), y(:)
        real(wp), allocatable :: sKnot(:,:), tx(:), ty(:)
        real(wp), allocatable :: sKnot1(:), t1(:)
    end type 

    type(mesh_t) :: mesh
    type(interface_t) :: interface
contains


    subroutine eval_sinusoidal_surface(x, cfg, xi, xi_r, xi_rr)
        real(wp),       intent(in)  :: x(:)
        type(inputs_t), intent(in)  :: cfg
        real(wp),       intent(out) :: xi(:), xi_r(:), xi_rr(:)
        real(wp) :: kw

        kw = 2.0_wp*pi/cfg%lambda_wave
        xi  =   cfg%A_wave * sin(kw*x + cfg%phi_wave)
        xi_r  =  cfg%A_wave * kw    * cos(kw*x + cfg%phi_wave)
        xi_rr = -cfg%A_wave * kw**2 * sin(kw*x + cfg%phi_wave)
    end subroutine


    subroutine read_interface(fname) 
        use hdf5
        use utils_module, only: transpose_2D_matrix
        use interp_spline_mod, only: create_s_interpoland
        implicit none

        ! input variables
        character(len=256) :: fname
        ! HDF5 variables
        integer(HID_T) :: file_space, file_id, dset_id,dspace_id, memspace
        integer :: error,a
        logical :: exists
        integer(HSIZE_T), dimension(2) :: dims_h, maxdims
        integer(HSIZE_T), dimension(1) :: dims_x, dims_y
        real, allocatable, dimension(:,:) :: s_temp
        integer :: ii, jj

        print*, "Reading surface ..."
        ! read of the 3 component of the mesh
        !-------------------------------------------------------------------------
        call h5open_f(error)
        ! the method
        call h5fopen_f (trim(fname), H5F_ACC_RDONLY_F, file_id, error)

        ! Read u 
        call h5lexists_f(file_id, 's', exists, error)
        if (error /= 0) then
            call h5fclose_f(file_id, error)
            print*, "Failed to check existence of dataset s "
            stop
        end if

        if (exists) then
            call h5dopen_f(file_id, 's', dset_id, error)
            call h5dget_space_f(dset_id,dspace_id,error)
            ! Get the filespace (how data is stored in the file)
            call h5dget_space_f(dset_id, file_space, error)
            call h5sget_simple_extent_dims_f(dspace_id, dims_h, maxdims, error)

            interface%ny = int(dims_h(1))
            interface%nx = int(dims_h(2))

                allocate(interface%s_2D(interface%nx, interface%ny))
                allocate(s_temp(interface%ny, interface%nx))

            ! call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_u, error)
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,s_temp, dims_h, error)
            call h5dclose_f(dset_id, error)

            call transpose_2D_matrix(s_temp, interface%s_2D)
            deallocate(s_temp)
        else
            print *, "Dataset h does not exist"
            ! Close the file
            call h5fclose_f(file_id, error)
            stop 
        end if

        ! reading x array 
        call h5lexists_f(file_id, 'x', exists, error)
        if (error /= 0) then
            print *, "Failed to check existence of dataset x " 
            call h5fclose_f(file_id, error)
            stop
        end if

        if (exists) then
            call h5dopen_f(file_id, 'x', dset_id, error)
            call h5dget_space_f(dset_id,dspace_id,error)
            call h5sget_simple_extent_dims_f(dspace_id, dims_x, maxdims, error)
            if (int(dims_x(1)) /= interface%nx) then
                print *, "x array and h does not have the same dimension"
                print *, dims_x(1)
                call h5fclose_f(file_id, error)
                stop 
            end if 

            if (.not. allocated(interface%x)) allocate(interface%x(interface%nx))
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,interface%x, dims_x, error)
            call h5dclose_f(dset_id, error)
        else
            print *, "Dataset x does not exist"
            ! Close the file
            call h5fclose_f(file_id, error)
            stop 
        end if

        ! reading y array 
        call h5lexists_f(file_id, 'y', exists, error)
        if (error /= 0) then
            print *, "Failed to check existence of dataset y " 
            call h5fclose_f(file_id, error)
            stop
        end if

        if (exists) then
            call h5dopen_f(file_id, 'y', dset_id, error)
            call h5dget_space_f(dset_id,dspace_id,error)
            call h5sget_simple_extent_dims_f(dspace_id, dims_y, maxdims, error)
            if (int(dims_y(1)) /= interface%ny) then
                print *, "y array and h does not have the same dimension"
                call h5fclose_f(file_id, error)
                stop 
            end if 

            if (.not. allocated(interface%y)) allocate(interface%y(interface%ny))
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,interface%y, dims_y, error)
            call h5dclose_f(dset_id, error)
        else
            print *, "Dataset y does not exist"
            ! Close the file
            call h5fclose_f(file_id, error)
            stop 
        end if
        ! Close the file
        call h5fclose_f(file_id, error)
        ! Close FORTRAN interface.
        call h5close_f(error)

        
        ! Create 1d interpoland
        if (interface%ny == 1 .and. interface%nx >= 4) then            ! s(x)
            ! interface%is_1d = .true.;  interface%dir = 1
            call create_s_interpoland(interface%s_2D(:,1), interface%x, interface%sKnot1, interface%t1)
        else if (interface%nx >= 4 .and. interface%ny >= 4) then        ! s(x,y)
            call create_s_interpoland(interface%s_2D, interface%x, interface%y, &
                        interface%sKnot, interface%tx, interface%ty)
        else
            print *, "read_interface: need >= 4 points in each non-singleton direction"
            stop
        end if
        ! call create_s_interpoland(interface%s_2D, interface%x, interface%y,  &
        !     interface%sKnot,interface%tx,interface%ty)

            print*, "done"
    end subroutine 

    function surface_at(r, tau, cfg) result(xi)
        real(wp),       intent(in) :: r(:), tau
        type(inputs_t), intent(in) :: cfg
        real(wp) :: xi(size(r))

        if (interface%ny == 1 ) then
            xi = interpolate_s(interface%sKnot1, interface%t1, r*cosd(tau) + cfg%Xs_wave)
        else 
            xi = interpolate_s(interface%sKnot, interface%tx, interface%ty, &
                r*cosd(tau) + cfg%Xs_wave, r*sind(tau) + cfg%Ys_wave)
        end if
    end function

    subroutine init_coarse_mesh(tau,cfg)
        use interp_spline_mod, only: interpolate_s
        implicit none

        type(inputs_t), intent(in)  :: cfg
        real(wp), intent(in) :: tau

        real(wp), allocatable :: dummy1(:), dummy2(:)
        integer :: ii
        real :: tau_rad

        ! conversion in rad
        tau_rad = tau*pi/180.0_wp
        if (tau_rad >= pi) then 
            tau_rad = tau_rad - 2*pi
        end if 

        ! define size 2D domain from 3D dimension and angle
        ! ---------------------------------------------------------------------
        if (tau_rad == 0) then
            mesh%Lx = cfg%Lx1
        elseif (tau_rad == -pi) then
            mesh%Lx = cfg%Lx2
        elseif (tau_rad == pi/2) then
            mesh%Lx = cfg%Ly1     
        elseif (tau_rad == -pi/2) then
            mesh%Lx = cfg%Ly2
        elseif (-pi/2<tau_rad .and.  tau_rad<0) then
            mesh%Lx = min(cfg%Lx1/abs(cos(tau_rad)),cfg%Ly2/abs(sin(tau_rad)))
        elseif (0<tau_rad .and.  tau_rad<pi/2) then
            mesh%Lx = min(cfg%Lx1/abs(cos(tau_rad)),cfg%Ly1/abs(sin(tau_rad)))
        elseif (pi/2<tau_rad .and.  tau_rad<pi) then
            mesh%Lx = min(cfg%Lx2/abs(cos(tau_rad)),cfg%Ly1/abs(sin(tau_rad)))
        elseif (-pi<tau_rad .and.  tau_rad<-pi/2) then
            mesh%Lx = min(cfg%Lx2/abs(cos(tau_rad)),cfg%Ly2/abs(sin(tau_rad)))
        end if

        mesh%dz_coarse = cfg%dz_out
        mesh%dx_coarse = cfg%dx_out

        mesh%nz_coarse = floor((cfg%H + cfg%Lz) / mesh%dz_coarse)
        allocate(mesh%z_coarse(mesh%nz_coarse))
        do ii = 1, mesh%nz_coarse
        mesh%z_coarse(ii) = -cfg%H + mesh%dz_coarse * (ii - 1)
        end do

        mesh%nx_coarse = floor(mesh%Lx/mesh%dx_coarse + 0.5_wp)
        allocate(mesh%x_coarse(mesh%nx_coarse))
        do ii = 1, mesh%nx_coarse
        mesh%x_coarse(ii) = mesh%dx_coarse * ii
        enddo 

        ! interface at the output ranges (slope and curvature are not needed here)
        allocate(mesh%xi_coarse(mesh%nx_coarse), mesh%js_coarse(mesh%nx_coarse))
        if (cfg%external_wave) then 
            mesh%xi_coarse = surface_at(mesh%x_coarse, tau, cfg)

            ! mesh%xi_coarse = interpolate_s(interface%sKnot, interface%tx, interface%ty, &
            !     mesh%x_coarse*cosd(tau) + cfg%Xs_wave, mesh%x_coarse*sind(tau)+cfg%Ys_wave)
        else if (cfg%sinus) then 
            allocate(dummy1(mesh%nx_coarse), dummy2(mesh%nx_coarse))
            call eval_sinusoidal_surface(mesh%x_coarse*cosd(tau), cfg, mesh%xi_coarse, dummy1, dummy2)
        else
            mesh%xi_coarse = 0.0_wp
        endif 

        ! z_coarse(js) <= xi < z_coarse(js+1)   (1-based)
        mesh%js_coarse = floor((mesh%xi_coarse - mesh%z_coarse(1))/mesh%dz_coarse + 1.0e-12_wp) + 1

        allocate(mesh%water_coarse(mesh%nx_coarse, mesh%nz_coarse))
        do ii = 1, mesh%nx_coarse
            mesh%water_coarse(ii,:) = mesh%z_coarse <= mesh%xi_coarse(ii)
        end do
    end subroutine


    subroutine clean_coarse_mesh()
        deallocate(mesh%x_coarse)
        deallocate(mesh%z_coarse)
        deallocate(mesh%js_coarse)
        deallocate(mesh%xi_coarse)
        deallocate(mesh%water_coarse)
    end subroutine



    subroutine init_fine_mesh(tau, freq, cfg)
        use interp_spline_mod, only: interpolate_s
        implicit none

        type(inputs_t), intent(in)  :: cfg
        real(wp), intent(in)  :: tau, freq
        real(wp) :: dx_target
        real(wp) :: kw
        real(wp), allocatable :: xm(:)

        integer :: ii

        mesh%dz = cfg%c_air / freq / cfg%dz_frac


        ! here we choose dx so we can subsample to dx_out without interpolation
        ! interpolation works in a CN scheme where the is 10+points per wavelength
        ! but it is not the case here 
        dx_target = cfg%c_air / freq / cfg%dx_frac
        mesh%nbout = max(1, ceiling(cfg%dx_out/dx_target - 1.0e-9_wp))   ! add nbout to mesh_t
        mesh%dx    = cfg%dx_out / real(mesh%nbout, wp)


        ! add bottom points for bottom pml, or fluid
        if (cfg%seabed_type == 0) then
            mesh%nb_bot = ceiling(cfg%d_pml_bot*cfg%c_water/freq/mesh%dz)
        else if (cfg%seabed_type == 3) then
            mesh%nb_bot = ceiling(cfg%d_pml_bot*cfg%cp_sb/freq/mesh%dz)   ! PML in the sediment
        else
            mesh%nb_bot = 0
        endif

        mesh%nz = floor((cfg%H + cfg%Lz + cfg%pml_size) / mesh%dz) + mesh%nb_bot
        allocate(mesh%z(mesh%nz))
        do ii = 1, mesh%nz
        mesh%z(ii) = -cfg%H - mesh%nb_bot*mesh%dz + mesh%dz * (ii - 1)
        end do

        mesh%nx = floor(mesh%Lx/mesh%dx + 0.5_wp)
        allocate(mesh%x(mesh%nx))
        do ii = 1, mesh%nx
        mesh%x(ii) = mesh%dx * ii
        enddo 


        ! interface bottom
        mesh%jb = 0
        if (cfg%seabed_type == 3) then
            mesh%jb = mesh%nb_bot + 1        ! z(jb) = -H, last sediment node
            mesh%z(mesh%jb) = -cfg%H         ! remove rounding error
        endif

        ! interface
        ! ---------------------------------------------------------------------
        ! interface xi(r), slope and curvature at each fine range  (here: sinusoid)
        allocate(mesh%xi(mesh%nx), &
        mesh%xi_r(mesh%nx), mesh%xi_rr(mesh%nx), mesh%js(mesh%nx))
        
        if (cfg%external_wave) then 
            mesh%xi = surface_at(mesh%x, tau, cfg)

            ! mesh%xi = interpolate_s(interface%sKnot, interface%tx, interface%ty, &
            !     mesh%x*cosd(tau) + cfg%Xs_wave, mesh%x*sind(tau)+cfg%Ys_wave)

            mesh%xi_r(2:mesh%nx-1) = (mesh%xi(3:mesh%nx) - mesh%xi(1:mesh%nx-2)) / (2*mesh%dx) 
            mesh%xi_rr(2:mesh%nx-1) = (mesh%xi(3:mesh%nx) - 2*mesh%xi(2:mesh%nx-1) + & 
                        mesh%xi(1:mesh%nx-2)) / (mesh%dx**2)

            ! mesh%js = floor((mesh%xi - mesh%z(1))/mesh%dz + 1.0e-12_wp) + 1
            !
            ! ! the IIM stencils use z(js-1 : js+2): keep them inside the domain
            ! if (minval(mesh%js) < 2 .or. maxval(mesh%js) > mesh%nz - 2) then
            !     print *, "init_fine_mesh: interface leaves the domain, js range =", minval(mesh%js), maxval(mesh%js)
            !     stop
            ! end if
        else if (cfg%sinus) then 

            ! kw = 2.0_wp*pi/cfg%lambda_wave
            ! mesh%xi  =  cfg%A_wave * sin(kw*mesh%x*cosd(tau) + cfg%phi_wave)
            ! mesh%xi_r  =  cfg%A_wave * kw    * cos(kw*mesh%x*cosd(tau) + cfg%phi_wave)
            ! mesh%xi_rr = -cfg%A_wave * kw**2 * sin(kw*mesh%x*cosd(tau) + cfg%phi_wave)
            allocate(xm(mesh%nx))
            xm    = (mesh%x - 0.5_wp*mesh%dx)*cosd(tau)
            xm(1) = mesh%x(1)*cosd(tau)                       ! starter: interface at x(1)
            kw = 2.0_wp*pi/cfg%lambda_wave
            mesh%xi    =  cfg%A_wave        * sin(kw*xm + cfg%phi_wave)
            mesh%xi_r  =  cfg%A_wave*kw     * cos(kw*xm + cfg%phi_wave)
            mesh%xi_rr = -cfg%A_wave*kw**2  * sin(kw*xm + cfg%phi_wave)
            deallocate(xm)
        else 
            mesh%xi = 0.0_wp
            mesh%xi_r = 0.0_wp
            mesh%xi_rr = 0.0_wp
            ! mesh%js = int(floor(-mesh%z(1)/mesh%dz +1e-9)) + 1
            ! ! clamp 0 
            ! if (abs(mesh%z(mesh%js(1))) < 1e-9_wp*mesh%dz) mesh%z(mesh%js(1)) = 0.0_wp
        endif 

            ! z(js) <= xi < z(js+1)   (1-based)
            mesh%js = floor((mesh%xi - mesh%z(1))/mesh%dz + 1.0e-12_wp) + 1

            ! the IIM stencils use z(js-1 : js+2): keep them inside the domain
            if (minval(mesh%js) < 2 .or. maxval(mesh%js) > mesh%nz - 2) then
                print *, "init_fine_mesh: interface leaves the domain, js range =", minval(mesh%js), maxval(mesh%js)
                stop
            end if
            if (cfg%seabed_type == 3 .and. minval(mesh%js) < mesh%jb + 3) then
                print *, "init_fine_mesh: surface too close to the seabed (H too small)"
                stop
            end if
            allocate(mesh%water_fine(mesh%nx, mesh%nz))
            do ii = 1, mesh%nx
            mesh%water_fine(ii,:) = mesh%z <= mesh%xi(ii)
            end do
        end subroutine


        subroutine clean_fine_mesh()
            deallocate(mesh%x)
            deallocate(mesh%z)
            deallocate(mesh%water_fine)
            deallocate(mesh%xi)
            deallocate(mesh%xi_r)
            deallocate(mesh%xi_rr)
            deallocate(mesh%js)
        end subroutine 


    end module

