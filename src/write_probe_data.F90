module mod_probe 
    real(8), dimension(:,:,:), allocatable :: span_average
    character(len=5) :: nxstr, nystr
end module mod_probe

subroutine to_file(input_filename, i,j)
    use mod_streams
    use mod_sys
    implicit none
    character(len=120) :: input_filename

    integer :: i, j, k
    real(mykind) :: rho, rhou, rhov, rhow, vwrite, uwrite, wwrite

    open(21,file=input_filename)
    write(21, *) "rho, u, v, w"

    do k = 1, nz
        rho =  w(1,i,j,k)
        rhou = w(2,i,j,k)
        rhov = w(3,i,j,k)
        rhow = w(4,i,j,k)

        uwrite = rhou / rho
        vwrite = rhov / rho
        wwrite = rhow / rho

        write(21, "(E36.30, A1, E36.30, A1, E36.30, A1, E36.30)") rho, ",", uwrite, ",", vwrite, ",", wwrite
    enddo

    close(21)
end subroutine to_file

subroutine write_spanwise_probe_helper(probe_number, i, icyc)
    implicit none
    integer :: i, j, probe_number, icyc
    character(len=120) :: filename

    write(filename, "(A20, I1.1, A1, I5.5, A9)") "csv_data/span_probe_", probe_number, "_", icyc, "_visc.csv"
    j = 5
    call to_file(filename, i, j)

    write(filename, "(A20, I1.1, A1, I5.5, A8)") "csv_data/span_probe_", probe_number, "_", icyc, "_log.csv"
    j = 53
    call to_file(filename, i, j)

    write(filename, "(A20, I1.1, A1, I5.5, A9)") "csv_data/span_probe_", probe_number, "_", icyc, "_free.csv"
    j = 168
    call to_file(filename, i, j)

end subroutine write_spanwise_probe_helper

subroutine write_streamwise_probe(probe_number, j)
    use mod_streams
    implicit none

    integer:: probe_number, i,j,k
    real(mykind) :: rho, rhou, rhov, rhow, uwrite, vwrite, wwrite
    character(len=70) :: filename

    k = nz/2

    write(filename, "(A22, I1.1, A1, I5.5, A4)") "csv_data/stream_probe_", probe_number, "_", icyc, ".csv"

    open(299,file=filename)
    write(299, *) "rho, u, v, w"

    do i=1,nx
        rho =  w(1,i,j,k)
        rhou = w(2,i,j,k)
        rhov = w(3,i,j,k)
        rhow = w(4,i,j,k)

        uwrite = rhou / rho
        vwrite = rhov / rho
        wwrite = rhow / rho

        write(299, "(E36.30, A1, E36.30, A1, E36.30, A1, E36.30)") rho, ",", uwrite, ",", vwrite, ",", wwrite

    enddo

end subroutine write_streamwise_probe

subroutine write_probe_data()
    use mod_streams
    use mod_sys

    implicit none
    integer :: i, j, k, dead
    character(len=120) :: filename

    ! The three probe locations 
    !  |____________________|____________________|___________________|____________________|
    !  |                    |                    |                   |                    |    
    !  |                    |                    |                   |                    |    
    !  |         X          |                   X|                   |          X         |    
    !  |                    |                    |                   |                    |    
    !  oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo    
    !  |                    |                    |                   |                    |    
    !  |                    |                    |                   |                    |    
    !  |                    |                    |                   |                    |    
    !  |         X          |                   X|                   |          X         |    
    !  |_________X__________|___________________X|___________________|__________X_________|
    ! 
    !  X is a span-wise probe (into the page)
    !  o is a streamwise probe in the midplane

    ! if we are the main thread
    if (nrank==0) then
        ! write probe data 1/4th in the x direction and half in y direction
        i = nx *2/4
        call write_spanwise_probe_helper(1, i, icyc)

    elseif (nrank==1) then
        ! write probe data 1/2 in the x direction and half in y direction
        i = nx-1
        call write_spanwise_probe_helper(2, i, icyc)

    elseif (nrank ==3) then
        ! write probe data in the middle of the fourth process area, which
        ! comes out to be ~3/4ths of the entire x distance
        write(filename, "(A24, I5.5, A4)") "csv_data/w_probe_data_3_", icyc, ".csv"
        i = nx*2/4
        call write_spanwise_probe_helper(3, i, icyc)
    endif

    call write_streamwise_probe(nrank+1, ny/2)

end subroutine write_probe_data

subroutine init_write_telaps()
    use mod_streams
    implicit none

    if (masterproc) then
        open(995, file="timesteps.csv", action="write", status="replace")
        write(995, *) "telaps"
    endif
end subroutine init_write_telaps

subroutine write_telaps(telaps_in)
    use mod_streams
    use mod_sys
    implicit none

    real(mykind) :: telaps_in

    if (masterproc) then
        open(995, file="timesteps.csv", action="write", position="append")
        write(995, "(E15.10, A1)") telaps_in, ","
    endif
end subroutine write_telaps

! average out all of the values in a span for the vector of conservative variables w
subroutine write_span_averaged
    use mod_streams
    use mod_sys
    use mod_probe

    ! local variables
    character(len=5) :: current_cycle
    character(len=2) :: mpi_process_number
    character(len=70) :: current_filename

    ! works for icyc < 99,999 iterations
    write(current_cycle, "(I5.5)") icyc
    current_filename = "spans/span_average_"//current_cycle//"_average.binary"

    allocate(span_average(nx, ny, 5))
    call helper_average_span(1)
    call helper_average_span(2)
    call helper_average_span(3)
    call helper_average_span(4)
    call helper_average_span(5)

    call write_array_bytes(current_filename)

    ! everything else is deallocated by moving the allocation to span_average 
    ! when writing to vtk
    ! we deallocate here so that we dont error when allocating at the start of 
    ! this subroutine on the second time that it is called
    deallocate(span_average)

end subroutine write_span_averaged

! A helper to average out the values on each of the spans for different variables
! Called from the write_span_averaged subroutine with the name of the file to write to
! and the slice of the variable that we are averaging
subroutine helper_average_span(slice_var)
    use mod_streams
    use mod_sys
    use mod_probe
    implicit none

    ! arguments
    integer:: slice_var

    ! local variables
    real(8) :: current_average
    integer:: i, j, k

    do i = 1,nx
        do j = 1, ny
            current_average = 0.0

            ! sum up all of the values on the z axis
            do k = 1, nz
                ! hopefully the compiler can inline each of these branches in the hot loop
                ! otherwise this is painfully slow
                if (slice_var == 1) then
                    ! if we are dealing with rho
                    current_average = current_average + w(slice_var,i,j,k)
                else if (slice_var ==2 .or. slice_var == 3 .or. slice_var == 4) then
                    ! if we are not dealing with rho then we need to divide by it
                    current_average = current_average + (w(slice_var,i,j,k) / w(1,i,j,k))
                else 
                    current_average = current_average + &
                        (w(2,i,j,k) / w(1,i,j,k))**2 + &
                        (w(3,i,j,k) / w(1,i,j,k))**2 + &
                        (w(4,i,j,k) / w(1,i,j,k))**2 
                end if

            enddo
            ! calculate the mean of the data
            ! TODO: hopefully this is not an overflow somewhere
            current_average = current_average / nz

            span_average(i,j, slice_var) = current_average
        enddo
    enddo

end subroutine helper_average_span

subroutine write_array_bytes(filename)
    use mod_streams
    use mod_probe
    implicit none

    ! the number of points in each of the directions
    character(len=16) :: curr_cycle
    character(len=70) :: filename
    ! the number of variables that will be sent over mpi
    integer :: count
    integer :: tag
    integer :: master_id
    integer :: i 

    count = nx * ny * 5
    master_id = 0

    if (masterproc) then 
        open(23, file=filename, form="unformatted", access="stream")

        ! write all of the current process' span_average information
        call write_current_array_file(23)

        write(*,*) "(master) wrote data for self"

        ! then, recieve all the information from the other MPI procs and write it to the file

        do i = 1,nproc-1
            write(*,*) "(master) wrote data for proc", i
            tag = i
            call MPI_RECV(span_average(1:nx, 1:ny, 1:5), count, mpi_prec, tag, tag, mpi_comm_world, istatus, iermpi)
            call write_current_array_file(23)
        end do
        
        close(23)
    else
        tag = nrank
        write(*,*) "from ", nrank, "sending data now to master"
        call MPI_SEND(span_average(1:nx, 1:ny, 1:5), count, mpi_prec, master_id, tag, mpi_comm_world, iermpi)
    endif


    close(23)
end subroutine write_array_bytes


! write the entire contents of span_average to a file
subroutine write_current_array_file(filehandle)
    use mod_streams
    use mod_probe
    implicit none

    integer :: filehandle
    integer :: arr_idx, i, j

    do i =1,nx
        do j = 1,ny
            do arr_idx = 1,5
                write(filehandle) span_average(i,j,arr_idx)
            enddo
        enddo
    enddo

end subroutine write_current_array_file
