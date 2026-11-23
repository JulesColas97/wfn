module  read_write

implicit none

public :: read_param
public :: write_side
public :: write_side2
public :: write_planes
public :: write_receiver
public :: write_2Dflow
public :: read_generic_h5flow3D

contains

!! @brief Read the inpu parameter file for the simulation by uysing namelist
!!
!! @details this function read the input2.dat files that contaions the parameter for the simulation
!! the input file must be in the same folder that the executable, all vaiables are define like
!! in fortran code (you can use ! to comment line and x=12.5 to define variables )
!!
!!
!! @param [inout]     case_name: the name of the simulation
!! @param [inout]     var0: the medium general parameters
!! @param [inout]     simu: the global 3D parameters for the case
!! @param [inout]     imp: parameters for the ground impedance model
!! @param [inout]     output: parameters for the solution output
!! @param [inout]     pml: Perfectly matched layer boundary condition parameters
!! @param [inout]     flow: mean 3D (or 2D) flow parameters
subroutine read_param(var0,simu,imp,output,pml,flow3D, terrain2D)
    use types

    implicit none
    ! input variables
    !--------------------------------------------------------------------

    type(variables_0), intent(inout)                        :: var0
    type(simulation), intent(inout)                         :: simu
    type(impedance),intent(inout)                           :: imp
    type(output_parameters), intent(inout)                  :: output
    type(perfectly_matched_layer),intent(inout)             :: pml
    type(meanflow3D),intent(inout)                          :: flow3D
    type(terrain_elevation),intent(inout)                   :: terrain2D

    !local variables
    !--------------------------------------------------------------------
    real                                                    :: ratio
    real                                                    :: pi = 3.1415
    complex                                                 :: eye = (0,1)
    integer                                                 :: ii


    ! temporary var0:
    real                                :: c_air, rho_air, gamma0, T0, P0
    ! temporary impedance
    real                                :: impedance_sigmae, impedance_alphae
    integer                             :: impedance_mode
    real :: impedance_Rs, impedance_Om, impedance_q2, impedance_e
    !temporary simu
    character(len=40)                   :: case_name
    integer, parameter                  :: MAX_N = 1000
    real, parameter                     :: UNSET = -1e30
    integer                             :: nb_tau=-1
    integer                             :: nb_frequencies
    real                                :: tau(MAX_N)
    real                                :: frequencies(MAX_N)

    real                                :: Lx1,Lx2,Ly1,Ly2,Lz
    real :: taumin, taumax, dtau 
    real                                :: source_height
    real                                :: dx_frac

    integer                             :: pe_model
    logical                             :: continuation ! not implemented yet, but implemented for GTPE

    !temporary output
    real                                :: dx_out
    integer                             :: nbout
    integer                             :: nb_receiver_heights
    real                                :: receiver_heights(MAX_N)
    real                                :: xplane(MAX_N)
    real                                :: yplane(MAX_N)
    INTEGER                             :: nb_xplane
    INTEGER                             :: nb_yplane
    logical                             :: side, write_flow

    ! temporary pml
    integer                             :: pml_size
    real                                :: pml_param
    real                                :: pml_n

    !temporary flow
    logical                             :: external_flow
    character(len=256)                  :: fname_flow
    real                                :: Xs_flow, Ys_flow
    logical                             :: uniform
    logical                             :: logarithmic
    real                                :: u_ref

    ! temporary terrain 
    logical                             :: external_terrain
    character(len=256)                  :: fname_terrain
    real                                :: Xs_terrain, Ys_terrain
    logical                             :: gaussian_hill
    real                                :: hill_length, hill_height, hill_pos
    logical :: multi_domain
    real :: domain_step
    integer :: i
    character(len=256) :: filename
    character(len=256) :: arg
    integer :: num_args
    logical :: file_exists

    ! define element to be read from input.nml
    namelist/inputs/ case_name, nb_tau, tau,taumin, taumax, frequencies, &
                    Lx1,Lx2,Ly1,Ly2,Lz,dx_frac, &                  
                    continuation,dx_out, nbout, receiver_heights,xplane,yplane, &
                    side,write_flow, &
                    case_name,&
                    c_air, rho_air, gamma0,T0, P0, &
                    impedance_mode, impedance_sigmae, impedance_alphae, &
                    impedance_Rs, impedance_Om, impedance_q2, impedance_e,&
                    pml_size,pml_param, pml_n,&
                    source_height, &
                    external_flow, uniform, logarithmic,u_ref, &
                    fname_flow, Xs_flow, Ys_flow, pe_model, &
                    external_terrain, fname_terrain, Xs_terrain, Ys_terrain, &
                    gaussian_hill, hill_height, hill_length, hill_pos, &
                    multi_domain, domain_step

    taumin = UNSET
    taumax = UNSET
    nb_tau = -1
    tau = UNSET
    frequencies = UNSET
    receiver_heights = UNSET
    xplane = UNSET
    yplane = UNSET

    ! Default values if input not present in the input.nml file
    case_name = 'c0'
    Lx1 = 300.0
    Lx2 = 300.0
    Ly1 = 300.0
    Ly2 = 300.0
    Lz = 100.0
    source_height = 119.0
    dx_frac = 2.0
    u_ref = 0.0
    fname_flow = './u.h5'
    Xs_flow = 0.0
    Ys_flow = 0.0
    fname_terrain = './terrain.h5'
    Xs_terrain = 0.0
    Ys_terrain = 0.0
    dx_out = 1.0
    nbout = 1
    
    pe_model = 0
    uniform = .false.
    external_flow = .false.
    logarithmic = .false.

    external_terrain = .false.
    gaussian_hill = .false.
    hill_height= 0.0
    hill_length = 0.0
    hill_pos = 0.0

    multi_domain = .false.
    domain_step = 0.0
    continuation =.false.
    
    ! default to rigid ground 
    impedance_mode = 0
    ! default parameters for variable porosity model 
    ! grass Attenborough 
    impedance_alphae =  100
    impedance_sigmae = 50000

    ! default parameters for porous model (fro Guibard Thesis) 
    impedance_Rs = 7000
    impedance_Om = 0.6 
    impedance_q2 = 1.3
    impedance_e = 0.13

    ! default PML values 
    pml_size=30
    pml_param=50000.0
    pml_n=2.5


    ! default output
    side = .true.
    write_flow = .true.

    ! default Atmospheric values 
    c_air = 340
    rho_air = 1.2
    gamma0 = 1.4
    ! temperature in Kelvin
    T0 = 290
    P0 = 101325

    filename = 'input.nml'
    ! Count arguments and loop through them
    num_args = command_argument_count()
    
    do i = 1, num_args
        call get_command_argument(i, arg)
        
        ! Check for the -f flag
        if (trim(arg) == '-f') then
            if (i < num_args) then
                call get_command_argument(i + 1, filename)
            else
                print *, "Error: -f flag provided without a filename."
                stop
            end if
        end if
    end do
    
    print *, "Using input file: ", trim(filename)
    inquire(file=trim(filename), exist=file_exists)
    if (.not. file_exists) then
        print*, "File does not exists... stopping"
        stop
    endif

    ! read inout file
    open(unit=1,file=trim(filename))
    read(unit=1,nml=inputs)
    close(1)

    var0%c=c_air
    var0%rho=rho_air
    var0%gamma=gamma0
    var0%T = T0
    var0%P = P0
    
    imp%mode = impedance_mode
    imp%Rs = impedance_Rs
    imp%Om = impedance_Om
    imp%q2 = impedance_q2
    imp%e = impedance_e
    imp%sigmae=impedance_sigmae
    imp%alphae=impedance_alphae

    ! because we are using define type some of elements must be read and the assign to the structure
    simu%case_name = case_name

    if (taumin > UNSET .and. taumax > UNSET .and. nb_tau > 0) then
        ! CASE A: Build from min/max/count
        print *, "Generating tau vector from taumin/taumax..."
        simu%nb_tau = nb_tau
        allocate(simu%tau(simu%nb_tau))
        if (nb_tau == 1) then
            simu%tau(1) = taumin
        else
            dtau = (taumax - taumin) / real(nb_tau)
            do i = 1, nb_tau
                simu%tau(i) = taumin + real(i-1) * dtau
            end do
        end if
        print*, simu%tau
    else
        ! CASE B: Use the vector provided in the NML (default: tau = 0)
        simu%nb_tau = count_set(tau)
        if (simu%nb_tau == 0) then
            tau(1) = 0.0
            simu%nb_tau = 1
        end if
        allocate(simu%tau(simu%nb_tau))
        simu%tau = tau(1:simu%nb_tau)
    end if

    nb_frequencies = count_set(frequencies)
    if (nb_frequencies == 0) then
        frequencies(1) = 100.0
        nb_frequencies = 1
    end if
    simu%nb_freq = nb_frequencies
    allocate(simu%frequencies(simu%nb_freq))
    simu%frequencies= frequencies(1:simu%nb_freq)

    simu%Lx1 = Lx1
    simu%Lx2 = Lx2
    simu%Ly1 = Ly1
    simu%Ly2 = Ly2
    simu%Lz = Lz
    simu%source_height = source_height
    simu%dx_frac = dx_frac
    simu%external_flow = external_flow
    simu%uniform       = uniform
    simu%logarithmic   = logarithmic
    simu%u0            = u_ref
    simu%pe_model      = pe_model 
    simu%continuation  = continuation 
    simu%multi_domain  = multi_domain
    simu%domain_step   = domain_step

    simu%hh            = hill_height
    simu%hl            = hill_length
    simu%hpos          = hill_pos

    if (pe_model==5 .and. multi_domain) then 
        print*, "cant use multi domain with GTPE, choose an other PE model"
        stop 
    end if 

    if (.not.(pe_model==5 .or. multi_domain) .and. (simu%hh/=0 .or. external_terrain)) then 
        print*, "warning you added terrain elevation but without using the GTPE model or multi domain"
        print*, "terrain elevation set to 0.0 every where"
        hill_height = 0 
        hill_length = 0 
    end if 
    
    pml%size = pml_size
    pml%param = pml_param
    pml%n = pml_n

    flow3D%fname = fname_flow
    flow3D%Xs = Xs_flow
    flow3D%Ys = Ys_flow
    
    terrain2D%ext = external_terrain
    terrain2D%fname = fname_terrain
    terrain2D%Xs  = Xs_terrain
    terrain2D%Ys  = Ys_terrain
    terrain2D%gaussian_hill  = gaussian_hill
    terrain2D%hh = hill_height
    terrain2D%hl = hill_length
    terrain2D%hpos = hill_pos

    output%dout = dx_out
    nb_receiver_heights = count_set(receiver_heights)
    if (nb_receiver_heights == 0) then
        receiver_heights(1) = 2.0
        nb_receiver_heights = 1
    end if
    nb_xplane = count_set(xplane)
    nb_yplane = count_set(yplane)
    output%nb_receiver = nb_receiver_heights
    allocate(output%heights(output%nb_receiver))
    output%heights = receiver_heights(1:output%nb_receiver)
    output%side = side
    output%gen_side = side
    output%write_flow = write_flow

    allocate(output%xplane(nb_xplane))
    output%xplane = xplane(1:nb_xplane)    
    allocate(output%yplane(nb_yplane))
    output%yplane = yplane(1:nb_yplane)

contains

    !> number of leading entries of x that differ from the UNSET sentinel
    integer function count_set(x) result(n)
        real, intent(in) :: x(:)
        n = 0
        do while (n < size(x))
            if (x(n+1) <= UNSET) exit
            n = n + 1
        end do
    end function count_set
end subroutine read_param

!! @briefWrite the solution cartographies  in h5 file
subroutine write_side2(case_name,tau,freq,sol,mesh)
    use hdf5
    use utils
    use types
    implicit none

    ! input variables
    character(len=40), intent(in)                           :: case_name
    type(solution), INTENT(IN)                              :: sol
    type(mesh_2D), INTENT(IN)                               :: mesh
    real,intent(in)                                         :: tau, freq

    ! HDF5 variables
    integer(HID_T)                                          :: file_id, group_id
    integer                                                 :: error,a
    character(len=50)                                       :: file_name, group_name
    character(len=10)                                       :: tau_str, freq_str
    integer                                                 :: nx, ny, ndom

    if (int(tau)==tau) then 
        tau_str = trim(adjustl(int_to_str(int(tau))))
    else 
        tau_str = trim(adjustl(float_to_str(tau)))
    end if
    file_name =  trim(case_name)
    file_name=trim('./'//trim(case_name)//'_'//trim(tau_str)//'.h5')

    !file=trim(adjustl(fname))//'_sol_.h5'
    ! Initialize fortran interface
    call h5open_f(error)
    ! create new file
    !call h5fcreate_f(fname, H5F_ACC_TRUNC_F,file_id,error)
    call h5fopen_f (file_name, H5F_ACC_RDWR_F, file_id, error)

    ! create a group (with frequency as name)
    a = int(freq)
    freq_str = adjustl(int_to_str(a))
    group_name = trim('/side/'//freq_str)
    call h5gcreate_f(file_id, group_name, group_id, error)

    ! write p imag and real
    ! if (allocated(sol%p_over_p0)) then 
    !     ndom = size(sol%p_over_p0,3)
    !     nx = size(sol%p_over_p0,1)
    !     ny = size(sol%p_over_p0,2)
    !     call write_3Darray(sol%p_over_p0, nx, ny,ndom, 'deltaL_full',group_id)
    !     ! call write_3Darray(real(sol%phi), nx, ny,ndom, 'phi_r',group_id)
    !     ! call write_3Darray(aimag(sol%phi), nx, ny,ndom, 'phi_i',group_id)
    ! endif 
    ! if (allocated(mesh%x_grid)) then 
    !     ndom = size(mesh%x_grid,3)
    !     nx = size(mesh%x_grid,1)
    !     ny = size(mesh%x_grid,2)
    !     call write_3Darray(mesh%x_grid, nx, ny,ndom, 'x_grid',group_id)
    !     call write_3Darray(mesh%z_grid, nx, ny,ndom, 'z_grid',group_id)
    ! endif
    if (allocated(mesh%Lxc)) then 
        call  write_1Darray(mesh%Lxc,mesh%ndomain,'Lxc',group_id)
        call  write_1Darray(mesh%Lzc,mesh%ndomain,'Lzc',group_id)
        call  write_1Darray(mesh%xc,mesh%ndomain,'xc',group_id)
        call  write_1Darray(mesh%zc,mesh%ndomain,'zc',group_id)
        call  write_1Darray(mesh%ac,mesh%ndomain,'ac',group_id)
    endif 
    nx = size(sol%deltaL_side,1)
    ny = size(sol%deltaL_side,2)
    call write_2Darray(sol%deltaL_side, nx, ny,'TL',group_id)

    ! nx = size(sol%deltaL_rec,1)
    ! ny = size(sol%deltaL_rec,2)
    ! call write_2Darray(sol%deltaL_rec, nx, ny,'deltaL_rec',group_id)

    ! close group
    call h5gclose_f(group_id, error)

    ! close file
    call h5fclose_f(file_id,error)
    !close fortran interface
    call h5close_f(error)
end subroutine write_side2


!! @briefWrite the solution cartographies  in h5 file
subroutine write_side(case_name,tau,freq,sol)
    use hdf5
    use utils
    use types
    implicit none

    ! input variables
    character(len=40), intent(in)                           :: case_name
    type(solution), INTENT(IN)                              :: sol
    real,intent(in)                                         :: tau, freq

    ! HDF5 variables
    integer(HID_T)                                          :: file_id, group_id
    integer                                                 :: error,a
    character(len=50)                                       :: file_name, group_name
    character(len=10)                                       :: tau_str, freq_str
    integer                                                 :: nx, ny

    nx = size(sol%deltaL_side,1)
    ny = size(sol%deltaL_side,2)


    if (int(tau)==tau) then 
        tau_str = trim(adjustl(int_to_str(int(tau))))
    else 
        tau_str = trim(adjustl(float_to_str(tau)))
    end if
    file_name =  trim(case_name)
    file_name=trim('./'//trim(case_name)//'_'//trim(tau_str)//'.h5')

    !file=trim(adjustl(fname))//'_sol_.h5'
    ! Initialize fortran interface
    call h5open_f(error)
    ! create new file
    !call h5fcreate_f(fname, H5F_ACC_TRUNC_F,file_id,error)
    call h5fopen_f (file_name, H5F_ACC_RDWR_F, file_id, error)

    ! create a group (with frequency as name)
    a = int(freq)
    freq_str = adjustl(int_to_str(a))
    group_name = trim('/side/'//freq_str)
    call h5gcreate_f(file_id, group_name, group_id, error)

    ! write p imag and real
  !  call write_2Darray(real(p_vect), ceiling((1.*nx)/nbout), ceiling((1.*ny)/nbout),'p_real',group_id)
  !  call write_2Darray(aimag(p_vect), ceiling((1.*nx)/nbout), ceiling((1.*ny)/nbout),'p_imag',group_id)
    call write_2Darray(sol%deltaL_side, nx, ny,'TL',group_id)


    ! close group
    call h5gclose_f(group_id, error)

    ! close file
    call h5fclose_f(file_id,error)
    !close fortran interface
    call h5close_f(error)
end subroutine write_side



!! @briefWrite the solution cartographies  in h5 file
subroutine write_planes(case_name,tau,freq,sol,output)
    use hdf5
    use utils
    use types
    implicit none

    ! input variables
    character(len=40), intent(in)                           :: case_name
    type(solution), INTENT(IN)                              :: sol
    type(output_parameters),INTENT(IN)                      :: output
    real,intent(in)                                         :: tau, freq

    ! HDF5 variables
    integer(HID_T)                                          :: file_id, group_id
    integer                                                 :: error,a
    character(len=50)                                       :: file_name, group_name
    character(len=10)                                       :: tau_str, freq_str
    integer                                                 :: nx, ny,nz



    if (int(tau)==tau) then 
        tau_str = trim(adjustl(int_to_str(int(tau))))
    else 
        tau_str = trim(adjustl(float_to_str(tau)))
    end if
    file_name =  trim(case_name)
    file_name=trim('./'//trim(case_name)//'_'//trim(tau_str)//'.h5')

    !file=trim(adjustl(fname))//'_sol_.h5'
    ! Initialize fortran interface
    call h5open_f(error)
    ! create new file
    !call h5fcreate_f(fname, H5F_ACC_TRUNC_F,file_id,error)
    call h5fopen_f (file_name, H5F_ACC_RDWR_F, file_id, error)

    ! create a group (with frequency as name)
    a = int(freq)
    freq_str = adjustl(int_to_str(a))
    group_name = trim('/planes/'//freq_str)
    call h5gcreate_f(file_id, group_name, group_id, error)

    ! write p imag and real
    !  call write_2Darray(real(p_vect), ceiling((1.*nx)/nbout), ceiling((1.*ny)/nbout),'p_real',group_id)
    !  call write_2Darray(aimag(p_vect), ceiling((1.*nx)/nbout), ceiling((1.*ny)/nbout),'p_imag',group_id)
    if (size(output%xplane)>0) then
        nx = size(sol%xplane,1)
        ny = size(sol%yplane,1)
        nz = size(sol%yplane,2)
        call write_2Darray(sol%xplane, nx, nz,'xplane',group_id)
        call write_1Darray(sol%ycoord,nx,'xycoord',group_id)
        call write_1Darray(output%xplane,nx,'xxcoord',group_id)
        call write_1Darray(sol%xcount,nx,'xcount',group_id)
    endif 
    if (size(output%yplane)>0) then 
        nx = size(sol%xplane,1)
        ny = size(sol%yplane,1)
        nz = size(sol%yplane,2)
        call write_2Darray(sol%yplane, ny, nz,'yplane',group_id)
        call write_1Darray(sol%xcoord,ny,'yxcoord',group_id)
        call write_1Darray(output%yplane,ny,'yycoord',group_id)
        call write_1Darray(sol%ycount,ny,'ycount',group_id)
    endif 

    ! close group
    call h5gclose_f(group_id, error)

    ! close file
    call h5fclose_f(file_id,error)
    !close fortran interface
    call h5close_f(error)
end subroutine write_planes



!! @brief Write receiver
subroutine write_receiver(case_name,tau,freq, nx,heights, receiver, mesh)
    use hdf5
    use utils
    use types
    implicit none

    ! input variables
    character(len=40), intent(in)                           :: case_name
    integer, intent(in)                                     :: nx
    real,dimension(:,:), intent(in)                          :: receiver
    real,intent(in)                                         :: tau, freq
    real,dimension(:),intent(in)                            :: heights
    type(mesh_2D), INTENT(IN)                               :: mesh

    ! HDF5 variables
    integer(HID_T)                                          :: file_id, group_id
    integer                                                 :: error,a,ii
    character(len=50)                                       :: file_name, group_name
    character(len=10)                                       :: tau_str, freq_str, height_str


    ! compute file name
    !-----------------------------------------------------------------------
    file_name =  trim(case_name)
    ! tau_str = adjustl(float_to_str(tau))
    if (int(tau)==tau) then 
        tau_str = trim(adjustl(int_to_str(int(tau))))
    else 
        tau_str = trim(adjustl(float_to_str(tau)))
    end if
    file_name=trim('./'//trim(case_name)//'_'//trim(tau_str)//'.h5')

    ! Initialize fortran interface
    call h5open_f(error)
    ! create new file
    call h5fopen_f (file_name, H5F_ACC_RDWR_F, file_id, error)

    ! create a group (with frequency as name)
    !-----------------------------------------------------------------------
    ! a = int(freq)
    ! freq_str = adjustl(int_to_str(a))
    if (int(freq)==freq) then 
        freq_str = trim(adjustl(int_to_str(int(freq))))
    else 
        freq_str = trim(adjustl(float_to_str(freq)))
    end if
    group_name = trim('/receiver/'//freq_str)

    call h5gcreate_f(file_id, group_name, group_id, error)
    call write_2Darray(receiver, nx, size(heights),'TL',group_id)
    
    if (allocated(mesh%Lxc)) then 
    call  write_1Darray(mesh%Lxc,mesh%ndomain,'Lxc',group_id)
    call  write_1Darray(mesh%Lzc,mesh%ndomain,'Lzc',group_id)
    call  write_1Darray(mesh%xc,mesh%ndomain,'xc',group_id)
    call  write_1Darray(mesh%zc,mesh%ndomain,'zc',group_id)
    call  write_1Darray(mesh%ac,mesh%ndomain,'ac',group_id)
    endif 
    call h5gclose_f(group_id, error)
    ! close file
    call h5fclose_f(file_id,error)
    !close fortran interface
    call h5close_f(error)
end subroutine write_receiver


subroutine write_mesh(case_name,tau,mesh)
    use hdf5
    use utils
    use types
    implicit none

    ! input variables
    character(len=40), intent(in)                           :: case_name
    type(mesh_2D), INTENT(IN)                               :: mesh
    real,intent(in)                                         :: tau

    ! HDF5 variables
    integer(HID_T)                                          :: file_id, group_id
    integer                                                 :: error
    character(len=50)                                       :: file_name
    character(len=10)                                       :: tau_str, freq_str



    ! tau_str = adjustl(float_to_str(tau))
    if (int(tau)==tau) then 
        tau_str = trim(adjustl(int_to_str(int(tau))))
    else 
        tau_str = trim(adjustl(float_to_str(tau)))
    end if
    file_name =  trim(case_name)
    file_name=trim('./'//trim(case_name)//'_'//trim(tau_str)//'.h5')
    

    ! Initialize fortran interface
    call h5open_f(error)
    ! create new file
    !call h5fcreate_f(fname, H5F_ACC_TRUNC_F,file_id,error)
    call h5fopen_f (file_name, H5F_ACC_RDWR_F, file_id, error)
    ! create group mesh
    call h5gcreate_f(file_id, 'mesh', group_id, error)

    !write mesh
    call  write_1Darray(mesh%x_coarse,size(mesh%x_coarse),'x',group_id)
    call  write_1Darray(mesh%z_coarse,size(mesh%z_coarse),'z',group_id)

    call  write_2Darray(mesh%x_coarse_grid,size(mesh%x_coarse),size(mesh%z_coarse),'x_grid',group_id)
    call  write_2Darray(mesh%z_coarse_grid,size(mesh%x_coarse),size(mesh%z_coarse),'z_grid',group_id)
    call  write_1Darray(mesh%h_coarse,size(mesh%x_coarse),'h',group_id)
    call  write_1real(tau,'tau',group_id)


    ! close group
    call h5gclose_f(group_id, error)

    ! close file
    call h5fclose_f(file_id,error)
    !close fortran interface
    call h5close_f(error)
end subroutine write_mesh

subroutine write_2Dflow(case_name, tau, u, c)
    use hdf5
    use utils
    use types
    implicit none

    ! input variables
    character(len=40), intent(in)                           :: case_name
    real, intent(in)                                        :: tau
    real,dimension(:,:), intent(in)                         :: u, c

    ! HDF5 variables
    integer(HID_T)                                          :: file_id, group_id
    integer                                                 :: error
    character(len=10)                                       :: tau_str
    character(len=50)                                       :: file_name

    if (int(tau)==tau) then
        tau_str = trim(adjustl(int_to_str(int(tau))))
    else
        tau_str = trim(adjustl(float_to_str(tau)))
    end if
    file_name=trim('./'//trim(case_name)//'_'//trim(tau_str)//'.h5')

    ! same file and same layout as EWAPE: group /flow with u and c on the coarse (x, z) mesh
    call h5open_f(error)
    call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, error)
    call h5gcreate_f(file_id, 'flow', group_id, error)
    call write_2Darray(u, size(u,1), size(u,2), 'u',group_id)
    call write_2Darray(c, size(c,1), size(c,2), 'c',group_id)
    call h5gclose_f(group_id, error)
    call h5fclose_f(file_id,error)
    call h5close_f(error)
end subroutine write_2Dflow

subroutine read_generic_h5flow3D(flow3D, c0) !,ratio,Lx, Ly, Lz,z_i, delta, flow3D)
  use hdf5
  use utils
  use types
  !implicit none

  ! input variables
  type(meanflow3D), intent(inout)                         :: flow3D
  real, intent(in)                                        :: c0
  !local variables
  integer                                                 :: nx, ny, nz 
  ! HDF5 variables
  integer(HID_T)                                          :: file_space, file_id, dset_id,dspace_id, memspace
  integer                                                 :: error,a
  logical                                                 :: exists
  character(len=10)                                       :: tau_str
  integer(HSIZE_T), dimension(3)                          :: dims_u, dims_v, dims_w,dims_c, maxdims
  integer(HSIZE_T), dimension(1)                          :: dims_x, dims_y, dims_z
  real, allocatable, dimension(:,:,:)                     :: u_temp
  real, allocatable, dimension(:)                         :: p
  integer                                                 :: ii, jj
  logical                                                 :: file_exists
  integer(hsize_t), dimension(3) :: start, stride, count, block

  ! read of the 3 component of the flow
  !-------------------------------------------------------------------------
  ! read u
  call h5open_f(error)
  ! try multiple time to open the file 
  ! (implementation due to multiple call to the same file when parallelizing the
  ! the method
  error = -1
  do while (error==-1)
    call h5fopen_f (trim(flow3D%fname), H5F_ACC_RDONLY_F, file_id, error)
    call sleep(1)
    ! print*, error
  end do
  
  ! Read u 
  call h5lexists_f(file_id, 'u', exists, error)
  if (error /= 0) then
      call h5fclose_f(file_id, error)
      call stop_error("Failed to check existence of dataset u ")
  end if

  if (exists) then
    call h5dopen_f(file_id, 'u', dset_id, error)
    call h5dget_space_f(dset_id,dspace_id,error)
    ! Get the filespace (how data is stored in the file)
    call h5dget_space_f(dset_id, file_space, error)
    call h5sget_simple_extent_dims_f(dspace_id, dims_u, maxdims, error)

    nz = int(dims_u(1))
    ny = int(dims_u(2))
    nx = int(dims_u(3))
    
    flow3D%nz = int(dims_u(1))
    flow3D%ny = int(dims_u(2))
    flow3D%nx = int(dims_u(3))
    ! Define the hyperslab: read the entire dataset, but reorder the axes

    if (.not. allocated(flow3D%u)) allocate(flow3D%u(nx, ny, nz))
    if (.not. allocated(u_temp)) allocate(u_temp(nz, ny, nx))

    ! call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_u, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_u, error)
    call h5dclose_f(dset_id, error)
    ! flow3D%u = 0.0
    
    call transpose_3D_matrix(u_temp, flow3D%u)
    
    deallocate(u_temp)
  else
      print *, "Dataset u does not exist"
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
    if (int(dims_x(1)) /= nx) then
        print *, "x array and u does not have the same dimension"
        print *, dims_x(1)
        call h5fclose_f(file_id, error)
        stop 
    end if 

    if (.not. allocated(flow3D%x)) allocate(flow3D%x(nx))
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,flow3D%x, dims_x, error)
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
    if (int(dims_y(1)) /= ny) then
        print *, "y array and u does not have the same dimension"
        call h5fclose_f(file_id, error)
        stop 
    end if 

    if (.not. allocated(flow3D%y)) allocate(flow3D%y(ny))
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,flow3D%y, dims_y, error)
    call h5dclose_f(dset_id, error)
  else
      print *, "Dataset y does not exist"
      ! Close the file
      call h5fclose_f(file_id, error)
      stop 
  end if


  ! reading z array 
  call h5lexists_f(file_id, 'z', exists, error)
  if (error /= 0) then
      print *, "Failed to check existence of dataset z " 
      call h5fclose_f(file_id, error)
      stop
  end if

  if (exists) then
    call h5dopen_f(file_id, 'z', dset_id, error)
    call h5dget_space_f(dset_id,dspace_id,error)
    call h5sget_simple_extent_dims_f(dspace_id, dims_z, maxdims, error)
    if (int(dims_z(1)) /= nz) then
        print *, "z array and u does not have the same dimension"
        call h5fclose_f(file_id, error)
        stop 
    end if 

    if (.not. allocated(flow3D%z)) allocate(flow3D%z(nz))
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,flow3D%z, dims_z, error)
    call h5dclose_f(dset_id, error)
  else
      print *, "Dataset z does not exist"
      ! Close the file
      call h5fclose_f(file_id, error)
      stop 
  end if

  ! Read not compulsory input (v, w, c)
  ! ---------------------------------------------------------------------------

  ! Read v
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
    if (any(dims_v /= dims_u)) then 
        print *, "u array and v does not have the same dimension, vset to 0.0"
        flow3D%v = 0.0
    else 
        nz = int(dims_v(1))
        ny = int(dims_v(2))
        nx = int(dims_v(3))
        print *, nx, ny, nz

        if (.not. allocated(flow3D%v)) allocate(flow3D%v(nx, ny, nz))
        if (.not. allocated(u_temp)) allocate(u_temp(nz, ny, nx))

        ! call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_u, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_v, error)
        call h5dclose_f(dset_id, error)
    
        call transpose_3D_matrix(u_temp, flow3D%v)
        deallocate(u_temp)
    end if 
  else
      print *, "Dataset v does not exist, v set to 0.0"
      flow3D%v = 0.0
  end if


  ! Read W
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
    if (any(dims_w /= dims_u)) then 
        print *, "u array and w does not have the same dimension, w set to 0.0"
        flow3D%w = 0.0
    else 
        nz = int(dims_w(1))
        ny = int(dims_w(2))
        nx = int(dims_w(3))

        if (.not. allocated(flow3D%w)) allocate(flow3D%v(nx, ny, nz))
        if (.not. allocated(u_temp)) allocate(u_temp(nz, ny, nx))

        ! call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_u, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_w, error)
        call h5dclose_f(dset_id, error)
    
        call transpose_3D_matrix(u_temp, flow3D%w)
        deallocate(u_temp)
    end if 
  else
      print *, "Dataset w does not exist, w set to 0.0"
      flow3D%w = 0.0
  end if



  ! Read c
  if (.not. allocated(flow3D%c)) allocate(flow3D%c(nx, ny, nz))
  call h5lexists_f(file_id, 'c', exists, error)
  if (error /= 0) then
      print *, "Failed to check existence of dataset c "
      flow3D%c = c0
  end if

  if (exists) then
    call h5dopen_f(file_id, 'c', dset_id, error)
    call h5dget_space_f(dset_id,dspace_id,error)
    call h5sget_simple_extent_dims_f(dspace_id, dims_c, maxdims, error)
    if (any(dims_c /= dims_u)) then
        print *, "u array and c does not have the same dimension, c set to c0"
        flow3D%c = c0
    else 

        nz = int(dims_c(1))
        ny = int(dims_c(2))
        nx = int(dims_c(3))
        print *, nx, ny, nz

        if (.not. allocated(flow3D%c)) allocate(flow3D%c(nx, ny, nz))
        if (.not. allocated(u_temp)) allocate(u_temp(nz, ny, nx))

        ! call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_u, error)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_c, error)
        call h5dclose_f(dset_id, error)
    
        call transpose_3D_matrix(u_temp, flow3D%c)
        deallocate(u_temp)
    end if 
  else
      print *, "Dataset c does not exist, c set to c0"
      flow3D%c = c0
  end if
  ! Close the file
  call h5fclose_f(file_id, error)
  ! Close FORTRAN interface.
  call h5close_f(error)

end subroutine read_generic_h5flow3D


subroutine read_h5_terrain2D(terrain2D)
  use hdf5
  use utils
  use types
  !implicit none

  ! input variables
  type(terrain_elevation), intent(inout)                  :: terrain2D
  ! HDF5 variables
  integer(HID_T)                                          :: file_space, file_id, dset_id,dspace_id, memspace
  integer                                                 :: error,a
  logical                                                 :: exists
  integer(HSIZE_T), dimension(2)                          :: dims_h, maxdims
  integer(HSIZE_T), dimension(1)                          :: dims_x, dims_y
  real, allocatable, dimension(:,:)                       :: h_temp
  integer                                                 :: ii, jj
  logical                                                 :: file_exists
  integer(hsize_t), dimension(3)                          :: start, stride, count, block

  ! read of the 3 component of the flow
  !-------------------------------------------------------------------------
  ! read u
  call h5open_f(error)
  ! try multiple time to open the file 
  ! (implementation due to multiple call to the same file when parallelizing the
  ! the method
  error = -1
  do while (error==-1)
    call h5fopen_f (trim(terrain2D%fname), H5F_ACC_RDONLY_F, file_id, error)
    call sleep(1)
    ! print*, error
  end do
  
  ! Read u 
  call h5lexists_f(file_id, 'h', exists, error)
  if (error /= 0) then
      call h5fclose_f(file_id, error)
      call stop_error("Failed to check existence of dataset h ")
  end if

  if (exists) then
    call h5dopen_f(file_id, 'h', dset_id, error)
    call h5dget_space_f(dset_id,dspace_id,error)
    ! Get the filespace (how data is stored in the file)
    call h5dget_space_f(dset_id, file_space, error)
    call h5sget_simple_extent_dims_f(dspace_id, dims_h, maxdims, error)

    terrain2D%ny = int(dims_h(1))
    terrain2D%nx = int(dims_h(2))
    
    if (.not. allocated(terrain2D%h)) allocate(terrain2D%h(terrain2D%nx, terrain2D%ny))
    if (.not. allocated(h_temp)) allocate(h_temp(terrain2D%ny, terrain2D%nx))

    ! call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,u_temp, dims_u, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,h_temp, dims_h, error)
    call h5dclose_f(dset_id, error)
    
    call transpose_2D_matrix(h_temp, terrain2D%h)
    deallocate(h_temp)
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
    if (int(dims_x(1)) /= terrain2D%nx) then
        print *, "x array and h does not have the same dimension"
        print *, dims_x(1)
        call h5fclose_f(file_id, error)
        stop 
    end if 

    if (.not. allocated(terrain2D%x)) allocate(terrain2D%x(terrain2D%nx))
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,terrain2D%x, dims_x, error)
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
    if (int(dims_y(1)) /= terrain2D%ny) then
        print *, "y array and h does not have the same dimension"
        call h5fclose_f(file_id, error)
        stop 
    end if 

    if (.not. allocated(terrain2D%y)) allocate(terrain2D%y(terrain2D%ny))
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,terrain2D%y, dims_y, error)
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

end subroutine read_h5_terrain2D
end module read_write
