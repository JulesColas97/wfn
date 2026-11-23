module utils

implicit none

private

public :: res_trid
public :: res_trid_2
public :: write_1Darray
public :: write_2Darray
public :: write_3Darray
public :: write_1integer
public :: write_1real
public :: init_file
public :: check_files_exist
public :: linspace
public :: stop_error
public :: float_to_str
public :: int_to_str
public :: transpose_3D_matrix
public :: transpose_2D_matrix
public :: prof_h
public :: prof
! public :: local_to_global
! public :: global_to_local

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

!! @brief Tridiagonal matrix inversion using Thomas alogrithm
!!
!! @details this function return the product between phi and a tri diagonal matruix define by (a,b,c)
!!
!!
!! @param [in]      a: lower diagonal
!! @param [in]      b: middle diagonal
!! @param [in]      c: upper diagonal
!! @param[inout]    phi: a vector
!! @param[inout]    dom: dimension of the matrix and vector phi
!! @return  this function returns the product between the tridiagonal matrix and the vector phi
subroutine res_trid(phi,a,b,c,dms)

    complex, intent(in)                         :: a(:)
    complex, intent(in)                         :: b(:)
    complex, intent(in)                         :: c(:)

    complex, intent(inout)                      :: phi(:)
    integer, intent(in)                         :: dms

    complex, allocatable, dimension(:)          :: u,a2,b2,c2
    integer                                     :: kk

    allocate(u(dms))
    allocate(a2(dms-1))
    allocate(c2(dms-1))
    allocate(b2(dms-1))
    u = phi
    a2 = a
    b2 = b
    c2 = c
    a2(dms-1)=c2(dms-1)/b2(dms)
    u(dms)=u(dms)/b2(dms)

    do kk=dms-1,2,-1
        b2(kk)=b2(kk)-c2(kk)*a2(kk)
        a2(kk-1)=a2(kk-1)/b2(kk)
        u(kk)=(u(kk)-c2(kk)*u(kk+1))/b2(kk)
    enddo

    b2(1)=b2(1)-c2(1)*a2(1)
    u(1)=(u(1)-c2(1)*u(2))/b2(1)

    do kk=2,dms
        u(kk)=u(kk)-a2(kk-1)*u(kk-1)
    enddo
    phi=u
end subroutine res_trid

subroutine res_trid_2(phi,u, a,b,c,dms)

    complex, intent(out)                         :: a(:)
    complex, intent(out)                         :: b(:)
    complex, intent(inout)                         :: c(:)
    complex, intent(out)                         :: u(:)

    complex, intent(out)                           :: phi(:)
    integer, intent(in)                            :: dms

    integer                                     :: kk

    a(dms-1)=c(dms-1)/b(dms)
    u(dms)=u(dms)/b(dms)

    do kk=dms-1,2,-1
        b(kk)=b(kk)-c(kk)*a(kk)
        a(kk-1)=a(kk-1)/b(kk)
        u(kk)=(u(kk)-c(kk)*u(kk+1))/b(kk)
    enddo

    b(1)=b(1)-c(1)*a(1)
    u(1)=(u(1)-c(1)*u(2))/b(1)

    do kk=2,dms
        u(kk)=u(kk)-a(kk-1)*u(kk-1)
    enddo
    phi=u
end subroutine res_trid_2

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

!! @brief Write a integer inside a h5 files
!!
!! @details this function create the data space and the data set  and write an integer inside an h5 file
!!
!!
!! @param [in]      a: the integer to be written
!! @param [in]      fname: the name of the variables (and hence name of the data set)
!! @param [in]      file_id: the file _id which must have be previously open
subroutine write_1integer(a,fname,file_id)

    use hdf5

    character(len=*),intent(in)                             :: fname
    integer,intent(in)                                      :: a
    integer(HID_T),intent(in)                               :: file_id

    integer(HSIZE_T),dimension(1)                           :: dim
    integer(HID_T)                                          :: dspace_id, dset_id
    integer                                                 :: error
    dim(1) = 1

    call h5screate_simple_f(1,dim,dspace_id,error)
    call h5dcreate_f(file_id,fname, H5T_NATIVE_INTEGER, dspace_id, dset_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,a, dim,error)
    call h5dclose_f(dset_id,error)
    call h5sclose_f(dspace_id,error)

end subroutine write_1integer

!! @brief Write a real inside a h5 files
!!
!! @details this function create the data space and the data set  and write a real double inside an h5 file
!!
!!
!! @param [in]      a: the double to be written
!! @param [in]      fname: the name of the variables (and hence name of the data set)
!! @param [in]      file_id: the file _id which must have be previously open
subroutine write_1real(a,fname,file_id)

    use hdf5

    character(len=*),intent(in)                             :: fname
    real,intent(in)                                         :: a
    integer(HID_T),intent(in)                               :: file_id

    integer(HSIZE_T),dimension(1)                           :: dim
    integer(HID_T)                                          :: dspace_id, dset_id
    integer                                                 :: error
    dim(1) = 1
    call h5screate_simple_f(1,dim,dspace_id,error)
    call h5dcreate_f(file_id,fname, H5T_NATIVE_DOUBLE, dspace_id, dset_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a, dim,error)
    call h5dclose_f(dset_id,error)
    call h5sclose_f(dspace_id,error)

end subroutine write_1real

!! @brief Write a 1D double array inside a h5 files
!!
!! @details this function create the data space and the data set  and write a double 1D array inside an h5 file
!!
!!
!! @param [in]      a: the arrayto be written
!! @parma [in]      n: the size of the vector
!! @param [in]      fname: the name of the variable (and hence name of the data set)
!! @param [in]      file_id: the file _id which must have be previously open
subroutine write_1Darray(a,n,fname,file_id)

    use hdf5

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

!! @brief Write a 2D double array inside a h5 files
!!
!! @details this function create the data space and the data set  and write a double 2D array inside an h5 file
!!
!!
!! @param [in]      a: the arrayto be written
!! @parma [in]      n: the size of the first dimension of the array
!! @param [in]      m: the size of the second dimension of the array
!! @param [in]      fname: the name of the variable (and hence name of the data set)
!! @param [in]      file_id: the file _id which must have be previously open
subroutine write_2Darray(a,n,m,fname,file_id)

    use hdf5

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


!! @brief Write a 3D double array inside a h5 files
!!
!! @details this function create the data space and the data set  and write a double 2D array inside an h5 file
!!
!!
!! @param [in]      a: the arrayto be written
!! @parma [in]      n: the size of the first dimension of the array
!! @param [in]      m: the size of the second dimension of the array
!! @param [in]      fname: the name of the variable (and hence name of the data set)
!! @param [in]      file_id: the file _id which must have be previously open
subroutine write_3Darray(a,n,m,p,fname,file_id)

    use hdf5

    character(len=*),intent(in)                            :: fname
    integer,intent(in)                                     :: n,m,p
    real,dimension(:,:,:),intent(in)                         :: a
    integer(HID_T),intent(in)                              :: file_id

    integer(HSIZE_T),dimension(3)                           :: dim
    integer(HID_T)                                          :: dspace_id, dset_id
    integer                                                 :: error
    dim(1) = n
    dim(2) = m
    dim(3) = p

    call h5screate_simple_f(3,dim,dspace_id,error)
    call h5dcreate_f(file_id,fname, H5T_NATIVE_DOUBLE, dspace_id, dset_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a, dim,error)
    call h5dclose_f(dset_id,error)
    call h5sclose_f(dspace_id,error)

end subroutine write_3Darray

!! @brief Init the H5 file where the solution will be written
!!
!! @details this function initialize the h5 files with the relative path given from the execution folder
!! the solution will be written in a file ./case_name_theta.h5
!! the matrix of p_real and imaginary will be stored in a group '/side/$frequency$' and the mesh will be stored
!! in a group 'mesh'
!!
!!
!! @param [in]      case_name: the name of the simulation
!! @parma [in]      theta: the angle of propagation
subroutine init_file(case_name,theta, heights)
    use hdf5
    implicit none

    character(len=40), intent(in)                           :: case_name
    real,dimension(:),intent(in)                            :: heights
    real,intent(in)                                         :: theta 
    integer(HID_T)                                          :: file_id, group_id
    integer                                                 :: error,a

    character(len=10)                                       :: theta_str
    character(len=100)                                      :: file_name

    file_name =  trim(case_name)
    if (int(theta)==theta) then 
        theta_str = trim(adjustl(int_to_str(int(theta))))
    else 
        theta_str = trim(adjustl(float_to_str(theta)))
    end if

    file_name=trim('./'//trim(case_name)//'_'//trim(theta_str)//'.h5')
    
    call h5open_f(error)
    ! create new file
    call h5fcreate_f(file_name, H5F_ACC_TRUNC_F,file_id,error)
    call h5gcreate_f(file_id, '/side', group_id, error)
    CALL h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/planes', group_id, error)
    CALL h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/receiver', group_id, error)
    call write_1Darray(heights, size(heights) ,'heights',group_id)
    CALL h5gclose_f(group_id, error)
    ! close file
    call h5fclose_f(file_id,error)
    !close fortran interface
    call h5close_f(error)
end subroutine init_file

function check_files_exist(case,angles) result(first_angle_index)
    implicit none
    ! Inputs: a list of integer angles
    character(len=32)   :: case
    real, intent(in) :: angles(:)
    integer             :: first_angle_index

    
    character(len=32)   :: filename
    character(len=4)    :: angle_str
    logical             :: file_exists
    integer             :: i

    first_angle_index = 0

    do i = 1, size(angles)
        ! Format the angle into a 4-character string with leading zeros
        ! e.g., 45 becomes "0045"
        write(angle_str, '(I4.4)') int(angles(i))
        
        ! Construct the filename: "c0_0045.h5"
        filename = trim(case) // "_" // angle_str // ".h5"
        
        ! Check if the file exists on the disk
        inquire(file=trim(filename), exist=file_exists)
        
        if (.not. file_exists) then
            print *, "Error: File not found: ", trim(filename)
            first_angle_index = i
            return ! Exit the function immediately on the first failure
        end if
    end do

end function check_files_exist


!! @brief Calculate the height of the ground with a analytical profile
!! depending on the length and the height of the hill and the angle of propagation.
!! The hill is a 2D Hill from the article of Liu (2019)
!!
!!  @param [in]     xi: 1d array of the xi mesh grid
!!  @param [in]     height:  height of the hill
!!  @param [in]     length: length of the hill
!!  @param [in]     theta : angle of propagation in radian
!!  @param [in]     xH    : position of the top of the hill
!!  @param [inout]  H: 1D vector of the height at the ground
!!  @param [inout]  Hp, first derivative along xi of H
!!  @param [inout]  Hpp: second derivative of H
subroutine prof(xi,height,length,theta,xH,H,Hp,Hpp)
  implicit none
  real, dimension(:), intent(in)                          :: xi
  real, intent(in)                                        :: height
  real, intent(in)                                        :: length
  real, intent(in)                                        :: theta
  real, intent(in)                                        :: xH
  real,dimension(:),intent(inout)                         :: H, Hp, Hpp

  logical                                                 :: m1(size(xi))
  real                                                    :: pi = 3.141592653589793

  m1 = logical((((xi*cos(theta))>(xH-length)) .and. ((xi*cos(theta))<(xH+length))))

  where(m1 .eqv. .TRUE.)
    H = height * cos (pi *(xi*cos(theta)-xH)/(2*length))**2
    Hp = -height*pi*cos(theta)*sin(pi*(xi*cos(theta)-xH)/length)/(2*length)
    Hpp = -height*pi*pi*cos(theta)*cos(theta)*cos(pi*(xi*cos(theta)-xH)/length) &
        /(2*length**2)
  elsewhere
    H = 0.
    Hp = 0.
    Hpp = 0.
  end where
end subroutine prof

!! @brief Calculate the height of the ground with a analytical profile
!! depending on the length and the height of the hill and the angle of propagation.
!! The hill is a 2D Hill from the article of Liu (2019)
!!
!!  @param [in]     xi: 1d array of the xi mesh grid
!!  @param [in]     height:  height of the hill
!!  @param [in]     length: length of the hill
!!  @param [in]     theta : angle of propagation in radian
!!  @param [in]     xH    : position of the top of the hill
!!  @result         H: 1D vector of the height at the ground
function prof_h(xi,height,length,theta,xH) result(HH)
    real, dimension(:), intent(in)                          :: xi
    real, intent(in)                                        :: height
    real, intent(in)                                        :: length
    real, intent(in)                                        :: theta
    real, intent(in)                                        :: xH
    real                                                    :: HH(size(xi))

    logical                                                 :: m1(size(xi))
    real                                                    :: pi = 3.141592653589793

    m1 = logical((((xi*cos(theta))>(xH-length)) .and. ((xi*cos(theta))<(xH+length))))
    where(m1 .eqv. .TRUE.)
      HH = height * cos (pi *(xi*cos(theta)-xH)/(2*length))**2
    elsewhere
      HH = 0.
    end where
end function prof_h


! function local_to_global(xp,zp, xc, zc, ac) result(x,z)
!     real, dimension(:), intent(in) :: xp, zp 
!     real, intent(in)               :: xc, zc, ac
!     real, dimension(:)             :: x, z 
!
!     x = xc + xp*cos(ac) - zp*sin(ac)
!     z = zc +xp*sin(ac) + zp*cos(ac)
! end function local_to_global
!
! function global_to_local(x,z, xc, zc, ac) result(xp,zp)
!     real, dimension(:), intent(in) :: x, z 
!     real, intent(in)               :: xc, zc, ac
!     real, dimension(:)             :: xp, zp 
!
!     xp = (x-xc)*cos(ac) + (z-zc)*sin(ac)
!     zp = -(x-xc)*sin(ac) +(z-zc)*cos(ac)
!
! end function global_to_local


!! @brief equivalent to the linspace function in Matlab
!!
!!  @param [in]     from: first value pf the vector
!!  @param [in]     to : last value of the vector
!!  @param [in]     n: number of element of the vector
!!  @param [inout]  array : vector of n elements from 'from' to 'to' with equal
!!                           step
subroutine linspace(from, to, n, array)
    real, intent(in)                                    :: from, to
    integer, intent(in)                                 :: n
    real,allocatable,dimension(:), intent(inout)        :: array

    real :: range
    integer :: i

    allocate(array(n))
    range = to - from

    if (n == 0) return

    if (n == 1) then
        array(1) = from
        return
    end if

    do i=1, n
        array(i) = from + range * (i - 1) / (n - 1)
    end do
end subroutine linspace

subroutine stop_error(msg)
! Aborts the program with nonzero exit code
!
! The statement "stop msg" will return 0 exit code when compiled using
! gfortran. stop_error() uses the statement "stop 1" which returns an exit code
! 1 and a print statement to print the message.
!
! Example
! -------
!
! call stop_error("Invalid argument")

character(len=*) :: msg ! Message to print on stdout
print *, msg
stop 1
end subroutine
end module utils

