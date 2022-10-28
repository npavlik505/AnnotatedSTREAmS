subroutine wrap_startmpi()
    call startmpi()
end subroutine

subroutine wrap_setup()
    use mod_streams, only: tauw_x
    call setup()
end subroutine

subroutine wrap_init_solver()
    call init_solver()
end subroutine

subroutine wrap_step_solver()
    call step_solver()
end subroutine

subroutine wrap_finalize_solver()
    call finalize_solver()
end subroutine

subroutine wrap_finalize()
    use mod_streams

    call finalize()
    call mpi_finalize(iermpi)
end subroutine

subroutine wrap_copy_gpu_to_cpu()
    ! this follows the calling convention used in 
    ! manage_solver.fF90, so I naive-ly just copy
    ! it here
    call updateghost()
    call prims()
    call copy_gpu_to_cpu()
end subroutine

subroutine wrap_copy_cpu_to_gpu()
    call copy_cpu_to_gpu()
end subroutine

! calculate wall shear stress
! this subtoutine is a HEAVILY chopped down version of `writestatbl`
! in writestatbl.f90
subroutine wrap_tauw_calculate()
    use mod_streams, only: tauw_x, w_av, mykind, y, nx, ny, ncoords, masterproc
    implicit none
!
    integer :: i, j
    real(mykind), dimension(nx, ny) :: ufav, vfav, wfav
    real(mykind) :: dudyw, dy
    real(mykind) :: rmuw
    real(mykind) :: tauw

    ! if this does not get called, w_av is not initialized
    ! correctly and the shear stress is NaN
    call stats2d()

    if (ncoords(3) == 0) then
!
        do j = 1, ny
            do i = 1, nx
                ufav(i, j) = w_av(13, i, j)/w_av(1, i, j)
                vfav(i, j) = w_av(14, i, j)/w_av(1, i, j)
                wfav(i, j) = w_av(15, i, j)/w_av(1, i, j)
            end do
        end do
!
! Mean boundary layer properties
!
        do i = 1, nx
            dudyw = (-22._mykind*ufav(i, 1) + 36._mykind*ufav(i, 2) - 18._mykind*ufav(i, 3) + 4._mykind*ufav(i, 4))/12._mykind
            dy = (-22._mykind*y(1) + 36._mykind*y(2) - 18._mykind*y(3) + 4._mykind*y(4))/12._mykind
            dudyw = dudyw/dy
            rmuw = w_av(20, i, 1)
            tauw = rmuw*dudyw

            ! store shear stress information in the mod_streams array to 
            ! be read in the output file
            tauw_x(i) = tauw
        end do
    end if
!
end subroutine
