program main

use types
use splines
use utils
use interpo
use read_write
use hdf5
use wape
use pwape
use gtpe
use create_mesh
use initialize_flow


implicit none
type(simulation)                                :: simu
type(output_parameters)                         :: output
type(frequency)                                 :: freq
type(impedance)                                 :: imp
type(variables_0)                               :: var0
type(perfectly_matched_layer)                   :: pml
type(solve)                                     :: M
type(solution)                                  :: sol
type(meanflow3D)                                :: flow3D
type(terrain_elevation)                         :: terrain2D

real, allocatable,dimension(:)                  :: frequencies
real, allocatable,dimension(:)                  :: angles
logical :: all_present
integer :: first_angle_index
integer                                         :: ii,kk,jj

! read parameters from input config file
!-------------------------------------------------------------------------------
print*, 'read parameters'
call read_param(var0,simu,imp,output,pml,flow3D,terrain2D)

! read flow field from hdf5 file
!-------------------------------------------------------------------------------
if (simu%external_flow) then
    call read_generic_h5flow3D(flow3D, var0%c)
endif

if (terrain2D%ext) then
    call read_h5_terrain2D(terrain2D)
endif

! if continuation start loop on the first not computed angle
!TODO
if (simu%continuation .eqv. .true.) then 
    first_angle_index = check_files_exist(simu%case_name, simu%tau)
    if (first_angle_index == 0) then 
        print*, "All angles where already computed"
        stop
    else if (first_angle_index>1) then 
        first_angle_index = first_angle_index - 1
    endif
else  
    first_angle_index = 1
endif

! start loop on angles
!##########################################################################
do kk = first_angle_index,size(simu%tau)
    print*, '==================='
    print*, 'init propa angle : ', simu%tau(kk)
    ! init domain and coarse mesh for recording
    ! call init_domain(simu%tau(kk), simu, terrain2D, mesh, dom, output)
    call init_coarse_mesh(simu%tau(kk), simu, terrain2D, pml, output, var0%c)

    ! initiliaze output file
    ! if (simu%continuation .eqv. .false.) then
    call init_file(simu%case_name,simu%tau(kk),  output%heights)
    ! endif 

    ! call init_file(case_name,simu%tau(kk), output%heights)
    ! write mesh in h5 file
    call write_mesh(simu%case_name,simu%tau(kk),mesh)

    ! Initialize 2D Flow 
    !-------------------------------------------------------------------------
    ! interpolation of the meanflow over the coarse grid
    call init_meanflow(kk,simu,mesh,flow3D,mesh%nx_coarse,mesh%nz_coarse,var0%c)


    ! write 2D inetrpolated meanflow in h5 file
    if (output%write_flow .eqv. .true.) then 
    call  write_2Dflow(simu%case_name, simu%tau(kk), flow2D%u(:,:), flow2D%c(:,:))
    print*, 'mean flow written'
    endif 

    ! loop over the frequencies
    !##########################################################################
    do ii = 1, simu%nb_freq
        print*, '-------------------'
        print*, 'start propa for frequency: ', simu%frequencies(ii)
        ! initialize mesh for propa at one frequency
        call init_mesh2(simu%tau(kk), simu%frequencies(ii),&
                        var0%c,simu, terrain2D, pml,output%dout)

        ! nbout for the meanflow interpolation over a coarse grid
        output%nbout = ceiling(output%dout/mesh%dx)
        
        
        ! Init vector for solution
        !-----------------------------------------------------------------------
        print*, "init receiver"
        call init_receiver(simu%tau(kk), output, sol)
        
        ! perform PE resolution
        !-----------------------------------------------------------------------
        if (simu%multi_domain .eqv. .true.) then 
            call  pwape_1angle_1freq(var0,simu%tau(kk),simu%frequencies(ii),&
                            simu%pe_model,simu%source_height,imp,pml, sol, output)
        else if (simu%pe_model /= 5) then 
            call  wape_1angle_1freq(var0,simu%tau(kk),simu%frequencies(ii),&
                            simu%pe_model,simu%source_height,imp,pml, sol, output)

        else if (simu%pe_model == 5) then 
           call init_meanflow_freq(simu%pe_model,mesh,output%nbout,mesh%nx,mesh%nz,var0%c,1)
           call propa_1angle_1freq_gtpe(var0, simu%tau(kk), simu%frequencies(ii),&
                    simu%source_height, imp, pml, flow2D, mesh, sol, output)
        end if
       
        ! writing solution
        !-----------------------------------------------------------------------
        print*, 'Writing solution'
        if (output%side) then
          call write_side2(simu%case_name,simu%tau(kk),simu%frequencies(ii),sol, mesh)
        endif 
        
        call write_receiver(simu%case_name,simu%tau(kk),simu%frequencies(ii),mesh%nx_coarse,&
                      & output%heights, sol%deltaL_rec, mesh)
        ! 
        ! call  write_planes(simu%case_name,simu%tau(kk),simu%frequencies(ii),sol,output)
        ! 
        !dealocate solution
        !-----------------------------------------------------------------------
        if (output%side) then
          deallocate(sol%deltaL_side)
          deallocate(sol%deltaL_side_intermediate)
        endif
        
        deallocate(sol%deltaL_rec)
        deallocate(sol%deltaL_rec_intermediate)
        
        ! deallocate mesh
        call clean_mesh(simu%pe_model)
        
        !deallocate flow
        if (simu%pe_model ==5 ) then 
            deallocate(flow2D%epsilon_eff)
        end if 
        
        
        if (size(output%xplane)>0) then 
          deallocate(sol%xplane)
          deallocate(sol%xcount)
          deallocate(sol%ycoord)
        endif
        
        if (size(output%yplane)>0) then 
          deallocate(sol%yplane)
          deallocate(sol%ycount)
          deallocate(sol%xcoord)
        endif 
    
    enddo
    output%side = output%gen_side
    ! end loop over frequency
    !##########################################################################
    call clean_coarse_mesh(terrain2D%ext)


    ! if (simu%external_flow) then 
        deallocate(flow2D%u)
        deallocate(flow2D%w)
        deallocate(flow2D%c)
        deallocate(flow2D%x)
        deallocate(flow2D%z)
        ! deallocate(flow2D%uKnot)
        ! deallocate(flow2D%utx)
        ! deallocate(flow2D%uty)
        ! deallocate(flow2D%wKnot)
        ! deallocate(flow2D%wtx)
        ! deallocate(flow2D%wty)
        ! deallocate(flow2D%cKnot)
        ! deallocate(flow2D%ctx)
        ! deallocate(flow2D%cty)
    ! endif 
! if delta L true then compute delta L over all frequencies
enddo
! end loop over angles 
!##########################################################################
end program main
