module flow_mod
    use kind_mod, only: wp
    use constants, only: PI, UNSET
    use input_params, only: inputs_t
    use mesh_mod, only: mesh_t
    use interp_plane_mod 
    implicit none
    private
    public :: flow3D_t, flow2D_t, init_fine_flow, clean_fine_flow, flow2D
    public :: read_generic_h5flow3D, init_coarse_flow, clean_coarse_flow

    type flow3D_t
        logical                             :: ext                              ! external file
        character(len=100)                  :: fname                            ! file name
        real, allocatable,dimension(:,:,:)  :: u, v,w, u_projected, c             ! 3D fields
        real, allocatable,dimension(:)      :: x, y, z                          ! coordinates
        integer                             :: nx, ny, nz                        ! coordinates
        real                                :: Xs, Ys, Zs, tau                ! source position in flow
    end type flow3D_t

    type flow2D_t
        real, allocatable, dimension (:)    :: utx, uty, ctx, cty, wtx, wty
        real, allocatable, dimension(:,:)   :: uKnot, wKnot
        real, allocatable, dimension(:,:)   :: cKnot
        real, allocatable,dimension(:,:)    :: u, v, w, c, epsilon, gamma, tau, m, n, n_eff ! flow quantities
        complex, allocatable, dimension(:,:) :: epsilon_eff
        real, allocatable,dimension(:,:)    :: h10,h12,h20,h22 ! flow quantities

        real, allocatable,dimension(:,:)    :: x, z                             ! coarse grid
        real                                :: Xs, Ys, Zs                ! WT position
        real                                :: u0
    end type flow2D_t

    type(flow2D_t) :: flow2D
    type(flow3D_t) :: flow3D

contains
    subroutine read_generic_h5flow3D(fname)
        use hdf5
        use utils_module, only: transpose_3D_matrix
        !implicit none

        ! input variables
        character(len=256), intent(in)                         :: fname
        !local variables
        integer                                                 :: nx, ny, nz 
        ! HDF5 variables
        integer(HID_T)                                          :: file_space, file_id, dset_id,dspace_id, memspace
        integer                                                 :: error,a
        logical                                                 :: exists
        character(len=10)                                       :: tau_str
        integer(HSIZE_T), dimension(3)                          :: dims, dims_u, dims_v, dims_w,dims_c, maxdims
        integer(HSIZE_T), dimension(1)                          :: dims_x, dims_y, dims_z
        real, allocatable, dimension(:,:,:)                     :: u_temp
        real, allocatable, dimension(:)                         :: p
        integer                                                 :: ii, jj


        ! open file 
        call h5open_f(error)
        call h5fopen_f (trim(fname), H5F_ACC_RDONLY_F, file_id, error)

        ! reading x array 
        ! ---------------------------------------------------------------------
        call h5lexists_f(file_id, 'x', exists, error)
        if (error /= 0) then
            print *, "Error: Failed to check existence of dataset x " 
            call h5fclose_f(file_id, error)
            stop
        end if

        if (exists) then
            call h5dopen_f(file_id, 'x', dset_id, error)
            call h5dget_space_f(dset_id,dspace_id,error)
            call h5sget_simple_extent_dims_f(dspace_id, dims_x, maxdims, error)
            nx = int(dims_x(1))

            if (.not. allocated(flow3D%x)) allocate(flow3D%x(nx))
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,flow3D%x, dims_x, error)
            call h5dclose_f(dset_id, error)
        else
            print *, "Error: Dataset x does not exist"
            ! Close the file
            call h5fclose_f(file_id, error)
            stop 
        end if

        ! reading y array 
        ! ---------------------------------------------------------------------
        call h5lexists_f(file_id, 'y', exists, error)
        if (error /= 0) then
            print *, "Error: Failed to check existence of dataset y " 
            call h5fclose_f(file_id, error)
            stop
        end if

        if (exists) then
            call h5dopen_f(file_id, 'y', dset_id, error)
            call h5dget_space_f(dset_id,dspace_id,error)
            call h5sget_simple_extent_dims_f(dspace_id, dims_y, maxdims, error)
            ny = int(dims_y(1))

            if (.not. allocated(flow3D%y)) allocate(flow3D%y(ny))
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,flow3D%y, dims_y, error)
            call h5dclose_f(dset_id, error)
        else
            print *, "Error: Dataset y does not exist"
            ! Close the file
            call h5fclose_f(file_id, error)
            stop 
        end if


        ! reading z array 
        ! ---------------------------------------------------------------------
        call h5lexists_f(file_id, 'z', exists, error)
        if (error /= 0) then
            print *, "Error: Failed to check existence of dataset z " 
            call h5fclose_f(file_id, error)
            stop
        end if

        if (exists) then
            call h5dopen_f(file_id, 'z', dset_id, error)
            call h5dget_space_f(dset_id,dspace_id,error)
            call h5sget_simple_extent_dims_f(dspace_id, dims_z, maxdims, error)
            nz = int(dims_z(1))

            if (.not. allocated(flow3D%z)) allocate(flow3D%z(nz))
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,flow3D%z, dims_z, error)
            call h5dclose_f(dset_id, error)
        else
            print *, "Error: Dataset z does not exist"
            ! Close the file
            call h5fclose_f(file_id, error)
            stop 
        end if

        dims = [dims_z(1), dims_y(1), dims_x(1)]

        ! Read not compulsory input (u, v, w, c)
        ! ---------------------------------------------------------------------------
        if (.not. allocated(flow3D%u)) allocate(flow3D%u(nx, ny, nz))
        if (.not. allocated(u_temp)) allocate(u_temp(nz, ny, nx))

        ! Read u 
        ! ---------------------------------------------------------------------------
        call h5lexists_f(file_id, 'u', exists, error)
        if (error /= 0) then
            call h5fclose_f(file_id, error)
            print*, "Failed to check existence of dataset u "
            stop
        end if

        if (exists) then
            call h5dopen_f(file_id, 'u', dset_id, error)
            call h5dget_space_f(dset_id,dspace_id,error)
            ! Get the filespace (how data is stored in the file)
            call h5dget_space_f(dset_id, file_space, error)
            call h5sget_simple_extent_dims_f(dspace_id, dims_u, maxdims, error)
            if (any(dims_u /= dims)) then 
                print *, "u array and x,y,z does not have the same dimension, uset to 0.0"
                print* , dims_u, dims
                flow3D%u = 0.0
            else 

            ! call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_u, error)
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_u, error)
            call h5dclose_f(dset_id, error)
            call transpose_3D_matrix(u_temp, flow3D%u)
            endif 
        else
            print *, "Dataset u does not exist, u set to 0.0"
            flow3D%u = 0.0
            ! Close the file
            call h5fclose_f(file_id, error)
        end if



        ! Read v
        ! ---------------------------------------------------------------------------
        if (.not. allocated(flow3D%v)) allocate(flow3D%v(nx, ny, nz))
        call h5lexists_f(file_id, 'v', exists, error)
        if (error /= 0) then
            print *, "Failed to check existence of dataset v " 
            flow3D%v = 0.0
        end if

        if (exists) then
            call h5dopen_f(file_id, 'v', dset_id, error)
            call h5dget_space_f(dset_id,dspace_id,error)
            call h5sget_simple_extent_dims_f(dspace_id, dims_v, maxdims, error)
            if (any(dims_v /= dims)) then 
                print *, "v array and (x,y,z) does not have the same dimension, vset to 0.0"
                flow3D%v = 0.0
            else 
                ! call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_u, error)
                call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_v, error)
                call h5dclose_f(dset_id, error)

                call transpose_3D_matrix(u_temp, flow3D%v)
            end if 
        else
            print *, "Dataset v does not exist, v set to 0.0"
            flow3D%v = 0.0
        end if


        ! Read W
        ! ---------------------------------------------------------------------------
        if (.not. allocated(flow3D%w)) allocate(flow3D%w(nx, ny, nz))
        call h5lexists_f(file_id, 'w', exists, error)
        if (error /= 0) then
            print *, "Failed to check existence of dataset w " 
            flow3D%w = 0.0
        end if

        if (exists) then
            call h5dopen_f(file_id, 'w', dset_id, error)
            call h5dget_space_f(dset_id,dspace_id,error)
            call h5sget_simple_extent_dims_f(dspace_id, dims_w, maxdims, error)
            if (any(dims_w /= dims)) then 
                print *, "w array and (x,y,z) does not have the same dimension, w set to 0.0"
                flow3D%w = 0.0
            else 
                ! call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_u, error)
                call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_w, error)
                call h5dclose_f(dset_id, error)

                call transpose_3D_matrix(u_temp, flow3D%w)
            end if 
        else
            print *, "Dataset w does not exist, w set to 0.0"
            flow3D%w = 0.0
        end if

        ! Read c
        ! ---------------------------------------------------------------------------
        if (.not. allocated(flow3D%c)) allocate(flow3D%c(nx, ny, nz))
        call h5lexists_f(file_id, 'c', exists, error)
        if (error /= 0) then
            print *, "Failed to check existence of dataset c "
            flow3D%c = UNSET
        end if

        if (exists) then
            call h5dopen_f(file_id, 'c', dset_id, error)
            call h5dget_space_f(dset_id,dspace_id,error)
            call h5sget_simple_extent_dims_f(dspace_id, dims_c, maxdims, error)
            if (any(dims_c /= dims_u)) then
                print *, "u array and c does not have the same dimension, c set to c0"
                flow3D%c = UNSET
            else 

                ! call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_u, error)
                call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_c, error)
                call h5dclose_f(dset_id, error)

                call transpose_3D_matrix(u_temp, flow3D%c)
            end if 
        else
            print *, "Dataset c does not exist, c will be set to c0"
            flow3D%c = UNSET
        end if
        ! Close the file
        call h5fclose_f(file_id, error)
        ! Close FORTRAN interface.
        call h5close_f(error)

        deallocate(u_temp)
        flow3D%nx = nx
        flow3D%ny = ny
        flow3D%nz = nz
    end subroutine read_generic_h5flow3D


    subroutine init_coarse_flow(tau, cfg, mesh)
        use interp_plane_mod 
        use constants, only:pi, kappa,  unset
        implicit none
        real(wp), intent(in) :: tau
        type(inputs_t), intent(in)  :: cfg
        type(mesh_t), intent(in) :: mesh

        real(wp) :: tau_rad 
        real(wp), allocatable :: u_projected(:,:,:)
        integer :: ii, kk

        
        tau_rad = tau * pi / 180.0_wp

        ! initialize 2D fields 
        allocate(flow2D%u(mesh%nx_coarse, mesh%nz_coarse))
        allocate(flow2D%c(mesh%nx_coarse, mesh%nz_coarse)) 
        allocate(u_projected(flow3D%nx, flow3D%ny, flow3D%nz)) 

        ! interpolate 3D fields to 2D planes 
        if (cfg%external_flow) then 
            u_projected = cosd(tau) * flow3D%u  + sind(tau)  * flow3D%v

            if (flow3D%nx==1 .and. flow3D%ny==1) then 
                print*, "Interpolating from 1D profile"
                call interp_to_plane(u_projected(1,1,:), flow3D%z, 2, &
                            mesh%x_coarse, mesh%z_coarse, flow2D%u)
                call interp_to_plane(flow3D%c(1,1,:), flow3D%z, 2, &
                            mesh%x_coarse, mesh%z_coarse, flow2D%c)
            else if (flow3D%ny==1) then 
                print*, "Interpolating from 2D flow field"
                call interp_to_plane(u_projected(:,1,:), flow3D%x, flow3D%z, &
                            mesh%x_coarse, mesh%z_coarse, cfg%Xs_flow, flow2D%u)
                call interp_to_plane(flow3D%c(:,1,:), flow3D%x, flow3D%z, &
                            mesh%x_coarse, mesh%z_coarse, cfg%Xs_flow, flow2D%c)
            else 
                print*, "Interpolating from 3D flow field"
                call interp_to_plane(u_projected, flow3D%x, flow3D%y, flow3D%z, &
                            mesh%x_coarse, mesh%z_coarse, tau, cfg%Xs_flow, cfg%Ys_flow, flow2D%u)

                call interp_to_plane(flow3D%c, flow3D%x, flow3D%y, flow3D%z, &
                            mesh%x_coarse, mesh%z_coarse, tau, cfg%Xs_flow, cfg%Ys_flow, flow2D%c)
            endif 

            ! clip wind speed to 0.0 below interface 
            ! where (mesh%water_coarse) flow2D%u = 0.0


                if (any(flow2D%c == UNSET)) then 
                    print*, "setting c profile to c_air and c_water default"

                    where (mesh%water_coarse)
                        flow2D%c = cfg%c_water
                    elsewhere 
                        flow2D%c = cfg%c_air
                    end where 
                end if
        
        ! compute analytical flow if set in cfg 
        else 
            ! constant c 
            where (mesh%water_coarse)
                flow2D%c = cfg%c_water
            elsewhere 
                flow2D%c = cfg%c_air
            end where 

            ! compute u profile  /
            if (cfg%power) then 
                do kk=1, mesh%nx_coarse
                do ii=mesh%js_coarse(kk)+1, mesh%nz_coarse
                flow2D%c(kk, ii) = flow2D%c(kk, ii) + cfg%u_ref*((mesh%z_coarse(ii)-mesh%xi_coarse(kk))/cfg%z_ref)**cfg%alpha
                enddo
                enddo

            elseif (cfg%logarithmic) then
                do kk=1,mesh%nx_coarse
                do ii=mesh%js_coarse(kk)+1, mesh%nz_coarse
                flow2D%c(kk, ii) = flow2D%c(kk, ii) + cfg%u_ref/kappa * log(mesh%z_coarse(ii)/cfg%z_roughness)
                enddo
                enddo
            elseif (cfg%uniform) then 
                do kk=1,mesh%nx_coarse
                flow2D%c(kk, mesh%js_coarse(kk)+1:mesh%nz_coarse)= flow2D%c(kk, mesh%js_coarse(kk)+1:mesh%nz_coarse) + cfg%u_ref
                enddo
            endif
        endif 
    end subroutine 

    subroutine clean_coarse_flow()
        implicit none
        deallocate(flow2D%c)
        deallocate(flow2D%u)
    end subroutine


    subroutine init_fine_flow(cfg, mesh, f)
        use constants, only: kappa, pi
        use atm_abs_mod
        use interp_plane_mod
        implicit none

        type(inputs_t), intent(in)  :: cfg
        type(mesh_t), intent(in) :: mesh
        real(wp), intent(in) :: f
        integer :: ii, kk
        
        real(wp), allocatable :: u_temp(:,:), c_temp(:,:)
        real(wp) :: alpha_air, eta_air=0.0, alpha_water, eta_water=0.0
        real(wp), allocatable :: eta(:,:)

        allocate(c_temp(mesh%nx, mesh%nz))
        allocate(u_temp(mesh%nx, mesh%nz))
        allocate(flow2D%epsilon_eff(mesh%nz, mesh%nx))
        allocate(eta(mesh%nx, mesh%nz))
        eta = 0.0_wp
        
        ! -----------------------------------------------------------------
        if (cfg%external_flow) then 

            call interp_to_plane(flow2D%u, mesh%x_coarse, mesh%z_coarse, &
                        mesh%x, mesh%z, 0.0_wp, u_temp)
            call interp_to_plane(flow2D%c, mesh%x_coarse, mesh%z_coarse, & 
                        mesh%x, mesh%z, 0.0_wp, c_temp)

            ! where (mesh%water_fine) u_temp = 0.0
            where (mesh%water_fine)
                u_temp = 0.0_wp
                c_temp = cfg%c_water
            elsewhere 
                c_temp = cfg%c_air
            end where

            do ii=1,mesh%nx
                do kk =1, mesh%nz
                flow2D%epsilon_eff(kk,ii) = (cfg%c_ref/(c_temp(ii,kk) + u_temp(ii,kk)))**2 - 1
                enddo
            enddo
        else 
            ! constant c 
            ! -----------------------------------------------------------------
            where (mesh%water_fine)
                c_temp = cfg%c_water
            elsewhere 
                c_temp = cfg%c_air
            end where 

            ! compute u profile  /
            ! -----------------------------------------------------------------
            if (cfg%power) then 
                do kk=1, mesh%nx
                do ii=mesh%js(kk)+1, mesh%nz 
                c_temp(kk, ii) = c_temp(kk, ii) + cfg%u_ref*((mesh%z(ii)-mesh%xi(kk))/cfg%z_ref)**cfg%alpha
                enddo
                enddo

            elseif (cfg%logarithmic) then
                do kk=1, mesh%nx
                do ii=mesh%js(kk)+1, mesh%nz 
                c_temp(kk, ii) = c_temp(kk, ii) + cfg%u_ref/kappa * log(mesh%z(ii)/cfg%z_roughness)
                enddo
                enddo
            elseif (cfg%uniform) then 
                do kk=1,mesh%nx
                c_temp(kk,mesh%js(kk)+1:mesh%nz) = c_temp(kk, mesh%js(kk)+1:mesh%nz) +cfg%u_ref
                enddo
            endif

            ! compute absorption profile 
            ! -----------------------------------------------------------------
            if (cfg%abs_air) then 
                alpha_air = atm_abs_air(cfg%T0, cfg%p0, cfg%rh, f)
                eta_air = alpha_air / (20 * log10(exp(1.0_wp)) * 2* pi * f /cfg%c_air)
                print*, "alpha_air, eta_air = ", alpha_air, eta_air
            endif 

            if (cfg%abs_water) then 
                alpha_water = atm_abs_water( f)
                eta_water = alpha_water / (20 * log10(exp(1.0_wp)) * 2* pi * f /cfg%c_water)
                print*, "alpha_water, eta_water = ", alpha_water, eta_water
            endif 

            where (mesh%water_fine) 
                eta = eta_water
            elsewhere
                eta = eta_air
            end where

            ! compute epsilon_eff 
            ! -----------------------------------------------------------------
            do ii=1,mesh%nx
                do kk =1, mesh%nz
                flow2D%epsilon_eff(kk,ii) = (cfg%c_ref/c_temp(ii,kk)*(1 + cmplx(0,eta(ii,kk))))**2 - 1
                enddo
            enddo

        endif
        deallocate(u_temp)
        deallocate(c_temp)
    end subroutine 


    subroutine clean_fine_flow()
        implicit none
        deallocate(flow2D%epsilon_eff)
    end subroutine
end module
