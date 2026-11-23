program main
    use input_params,  only: inputs_t, read_inputs
    use kind_mod, only: wp
    use utils_module, only: ensure_output_dir
    use mesh_mod
    use flow_mod
    use write_module
    use splitstep_pe_wavy_mod, only: run_splitstep_pe_wavy
    use, intrinsic :: ieee_arithmetic
    implicit none
    type(inputs_t) :: cfg
    integer :: kk, ii
    complex(wp), allocatable, dimension(:,:) :: p_out_side, p_out_rec
    real(wp), allocatable, dimension(:,:) :: TL

    call ieee_set_underflow_mode(gradual=.false.)
    call read_inputs(cfg)
    call ensure_output_dir(cfg%case_name)
    if (cfg%external_flow) call read_generic_h5flow3D(cfg%fname_flow)
    if (cfg%external_wave) call read_interface(cfg%fname_wave)
    ! TODO read surface wave

    ! loop over tau 
    do kk = 1,size(cfg%tau)

    ! build coarse output mesh 
    call init_coarse_mesh(cfg%tau(kk), cfg)

    call init_coarse_flow(cfg%tau(kk), cfg, mesh)



    call init_file(cfg%case_name, cfg%tau(kk), cfg%receiver_heights)
    call write_mesh(cfg%case_name, cfg%tau(kk), mesh)
    call write_flow(cfg%case_name, cfg%tau(kk), flow2D)

    
    ! loop over frequencies 
    do ii = 1, size(cfg%frequencies)
        ! build mesh
        call init_fine_mesh(cfg%tau(kk), cfg%frequencies(ii), cfg)
        ! TODO init fine wave 

        ! interpolate flow 
        call init_fine_flow(cfg, mesh, cfg%frequencies(ii))

        ! run pe 
        call run_splitstep_pe_wavy(cfg%frequencies(ii), cfg, mesh, flow2D, p_out_side, p_out_rec)

        ! save to files
        allocate(TL(mesh%nx_coarse, mesh%nz_coarse)) 
        TL = 20*log10(abs(p_out_side))

        call write_side(cfg%case_name,cfg%tau(kk),cfg%frequencies(ii),TL)
        deallocate(TL)

        allocate(TL(mesh%nx_coarse, size(cfg%receiver_heights))) 
        TL = 20*log10(abs(p_out_rec))
        call write_receiver(cfg%case_name,cfg%tau(kk),cfg%frequencies(ii),TL)
        deallocate(TL)

        call clean_fine_flow()
        call clean_fine_mesh()

    enddo

    call clean_coarse_mesh()
    call clean_coarse_flow()


    enddo 

end program main
