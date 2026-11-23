module write_module

    public :: init_file, write_side, write_mesh, write_flow, write_receiver
    private :: write_1Darray, write_2Darray

contains 


    subroutine init_file(case_name,tau, heights)
        use hdf5
        use utils_module, only: int_to_str, float_to_str
        implicit none

        character(len=40), intent(in)                           :: case_name
        real,dimension(:),intent(in)                            :: heights
        real,intent(in)                                         :: tau 
        integer(HID_T)                                          :: file_id, group_id
        integer                                                 :: error,a

        character(len=10)                                       :: tau_str
        character(len=100)                                      :: file_name

        file_name =  trim(case_name)
        if (int(tau)==tau) then 
            tau_str = trim(adjustl(int_to_str(int(tau))))
        else 
            tau_str = trim(adjustl(float_to_str(tau)))
        end if

        file_name=trim('./'//trim(case_name)//'_'//trim(tau_str)//'.h5')

        call h5open_f(error)
        ! create new file
        call h5fcreate_f(file_name, H5F_ACC_TRUNC_F,file_id,error)
        call h5gcreate_f(file_id, '/side', group_id, error)
        CALL h5gclose_f(group_id, error)
        call h5gcreate_f(file_id, '/receiver', group_id, error)
        call write_1Darray(heights, size(heights) ,'heights',group_id)
        CALL h5gclose_f(group_id, error)
        ! close file
        call h5fclose_f(file_id,error)
        !close fortran interface
        call h5close_f(error)
    end subroutine init_file

    subroutine write_mesh(case_name,tau,mesh)
        use hdf5
        use mesh_mod, only: mesh_t
        use utils_module, only: int_to_str, float_to_str
        implicit none
    
        ! input variables
        character(len=40), intent(in)                           :: case_name
        type(mesh_t), INTENT(IN)                               :: mesh
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
        call  write_1Darray(mesh%xi_coarse,size(mesh%x_coarse),'s',group_id)
    
        ! close group
        call h5gclose_f(group_id, error)
    
        ! close file
        call h5fclose_f(file_id,error)
        !close fortran interface
        call h5close_f(error)
    end subroutine write_mesh

    subroutine write_flow(case_name,tau,flow2D)
        use hdf5
        use kind_mod, only: wp
        use flow_mod, only: flow2D_t
        use utils_module, only: int_to_str, float_to_str
        implicit none
    
        ! input variables
        character(len=40), intent(in)                           :: case_name
        type(flow2D_t), INTENT(IN)                               :: flow2D
        real(wp),intent(in)                                         :: tau
    
        ! HDF5 variables
        integer(HID_T)                                          :: file_id, group_id
        integer                                                 :: error
        integer :: nx, nz
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
        call h5gcreate_f(file_id, 'flow', group_id, error)
    
        !write mesh
        nx = size(flow2D%u,1)
        nz = size(flow2D%u,2)
        call  write_2Darray(flow2D%u,nx,nz,'u',group_id)
        call  write_2Darray(flow2D%c,nx,nz,'c',group_id)
    
        ! close group
        call h5gclose_f(group_id, error)
    
        ! close file
        call h5fclose_f(file_id,error)
        !close fortran interface
        call h5close_f(error)
    end subroutine write_flow


    subroutine write_side(casename,tau,freq,TL)
        use hdf5
        use utils_module, only: int_to_str, float_to_str
        use kind_mod, only: wp
        implicit none

        ! input variables
        character(len=40), intent(in)                           :: casename
        ! complex(wp), intent(in), dimension(:,:)                 :: TL
        real(wp), intent(in), dimension(:,:)                    :: TL
        real(wp),intent(in)                                     :: tau, freq

        ! HDF5 variables
        integer(HID_T)                                          :: file_id, group_id
        integer                                                 :: error,a
        character(len=50)                                       :: file_name, group_name
        character(len=10)                                       :: tau_str, freq_str
        integer                                                 :: nx, nz, ndom

        if (int(tau)==tau) then 
            tau_str = trim(adjustl(int_to_str(int(tau))))
        else 
            tau_str = trim(adjustl(float_to_str(tau)))
        end if
        file_name =  trim(casename)
        file_name=trim('./'//trim(casename)//'_'//trim(tau_str)//'.h5')

        !file=trim(adjustl(fname))//'_sol_.h5'
        ! Initialize fortran interface
        call h5open_f(error)
        ! create new file
        !call h5fcreate_f(fname, H5F_ACC_TRUNC_F,file_id,error)
        call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, error)

        ! create a group (with frequency as name)
        a = int(freq)
        freq_str = adjustl(int_to_str(a))
        group_name = trim('/side/'//freq_str)
        call h5gcreate_f(file_id, group_name, group_id, error)

        nx = size(TL,1)
        nz = size(TL,2)
        call write_2Darray(TL, nx, nz,'TL',group_id)

        ! nx = size(sol%deltaL_rec,1)
        ! ny = size(sol%deltaL_rec,2)
        ! call write_2Darray(sol%deltaL_rec, nx, ny,'deltaL_rec',group_id)

        ! close group
        call h5gclose_f(group_id, error)

        ! close file
        call h5fclose_f(file_id,error)
        !close fortran interface
        call h5close_f(error)

    end subroutine




     subroutine write_receiver(casename, tau, freq, TL)
        use hdf5
        use utils_module, only: int_to_str, float_to_str
        use kind_mod, only: wp
        implicit none

        character(len=40), intent(in)     :: casename
        real(wp), intent(in), dimension(:,:) :: TL          ! (nx_coarse, n_receivers)
        real(wp), intent(in)              :: tau, freq

        integer(HID_T)                    :: file_id, group_id
        integer                           :: error
        character(len=50)                 :: file_name, group_name
        character(len=10)                 :: tau_str, freq_str

        if (int(tau)==tau) then
            tau_str = trim(adjustl(int_to_str(int(tau))))
        else
            tau_str = trim(adjustl(float_to_str(tau)))
        end if
        file_name = trim('./'//trim(casename)//'_'//trim(tau_str)//'.h5')

        call h5open_f(error)
        call h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, error)

        freq_str   = adjustl(int_to_str(int(freq)))
        group_name = trim('/receiver/'//freq_str)
        call h5gcreate_f(file_id, group_name, group_id, error)
        call write_2Darray(TL, size(TL,1), size(TL,2), 'TL', group_id)
        call h5gclose_f(group_id, error)

        call h5fclose_f(file_id, error)
        call h5close_f(error)
    end subroutine


    subroutine write_1Darray(a,n,fname,file_id)
        use hdf5
        implicit none

        character(len=*),intent(in)                             :: fname
        integer,intent(in)                                      :: n
        real,dimension(:),intent(in)                            :: a
        integer(HID_T),intent(in)                               :: file_id

        integer(HSIZE_T),dimension(1)                           :: dim
        integer(HID_T)                                          :: dspace_id, dset_id
        integer                                                 :: error
        dim(1) = n

        call h5screate_simple_f(1,dim,dspace_id,error)
        call h5dcreate_f(file_id,fname, H5T_NATIVE_DOUBLE, dspace_id, dset_id,error)
        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a, dim,error)
        call h5dclose_f(dset_id,error)
        call h5sclose_f(dspace_id,error)

    end subroutine write_1Darray

    subroutine write_2Darray(a,n,m,fname,file_id)
        use hdf5
        implicit none
    
        character(len=*),intent(in)                            :: fname
        integer,intent(in)                                     :: n,m
        real,dimension(:,:),intent(in)                         :: a
        integer(HID_T),intent(in)                              :: file_id
    
        integer(HSIZE_T),dimension(2)                           :: dim
        integer(HID_T)                                          :: dspace_id, dset_id
        integer                                                 :: error
        dim(1) = n
        dim(2) = m
    
        call h5screate_simple_f(2,dim,dspace_id,error)
        call h5dcreate_f(file_id,fname, H5T_NATIVE_DOUBLE, dspace_id, dset_id,error)
        call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a, dim,error)
        call h5dclose_f(dset_id,error)
        call h5sclose_f(dspace_id,error)
    
    end subroutine write_2Darray
end module 

