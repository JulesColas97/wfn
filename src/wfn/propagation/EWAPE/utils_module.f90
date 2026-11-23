module utils_module

    implicit none

    private

    public :: float_to_str, int_to_str, transpose_2D_matrix, transpose_3D_matrix, &
        ensure_output_dir


contains

    !! convert float to string
    character(len=20) function float_to_str(k)
        real, intent(in) :: k
        write (float_to_str, '(F5.1)') k
        float_to_str = adjustl(float_to_str)
    end function float_to_str

    !! convert integer to string
    character(len=20) function int_to_str(k)
        integer, intent(in) :: k
        write (int_to_str, '(I5.4)') k
        int_to_str = adjustl(int_to_str)
    end function int_to_str

    subroutine ensure_output_dir(casename)
        implicit none
        character(len=*), intent(in) :: casename
        integer :: islash, stat

        islash = index(trim(casename), '/', back=.true.)   ! last slash
        if (islash == 0) return                            ! no directory part, e.g. "c0"

        call execute_command_line('mkdir -p "'//casename(1:islash-1)//'"', exitstat=stat)
        if (stat /= 0) then
            print *, "could not create directory: ", casename(1:islash-1)
            stop
        end if
    end subroutine

    subroutine transpose_2D_matrix(a, a_transpose)

        real, intent (in), dimension(:,:)               :: a
        real, intent (out),allocatable, dimension(:,:)  :: a_transpose
        integer, dimension(2)                           :: dims 
        integer                                         :: i,j

        dims = [size(a,1), size(a,2)]

        if (.not. allocated(a_transpose)) allocate(a_transpose(dims(2), dims(1)))
        do i = 1, dims(1)
        do j = 1, dims(2)
        a_transpose(j,i) = a(i,j)
        end do 
        end do 
    end subroutine transpose_2D_matrix


    subroutine transpose_3D_matrix(a, a_transpose)
        real, intent (in), dimension(:,:,:)               :: a
        real, intent (out),allocatable, dimension(:,:,:)  :: a_transpose
        integer, dimension(3)                             :: dims 
        integer                                           :: i,j,k

        dims = [size(a,1), size(a,2), size(a,3)]

        if (.not. allocated(a_transpose)) allocate(a_transpose(dims(3), dims(2), dims(1)))
        do i = 1, dims(1)
        do j = 1, dims(2)
        do k = 1, dims(3)
        a_transpose(k,j,i) = a(i,j,k)
        end do
        end do 
        end do 
    end subroutine transpose_3D_matrix


end module 
