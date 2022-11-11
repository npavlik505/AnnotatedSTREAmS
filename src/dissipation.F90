! calculate the dissipation rate. See lab streams documentation for
! `Calculated Quantities` for the equation being calculated
!
! after calculation, the value is stored in `dissipation_rate` in mod_streams
subroutine dissipation_calculation()
    use mod_streams

    implicit none

    integer :: i, j, k

    real*8 :: uu, vv, ww, rho
    real*8 :: dudx, dvdy, dwdz
    real*8 :: dx, dy, dz
    real*8 :: mpi_dissipation_sum

    ! dissipation_rate is defined in mod_streams, we reset its value here
    dissipation_rate = 0.0

    mpi_dissipation_sum = 0.0

    ! w_gpu(:, :, :, 1) -> rho
    ! w_gpu(:, :, :, 2) -> rho u
    ! w_gpu(:, :, :, 3) -> rho v
    ! w_gpu(:, :, :, 4) -> rho w
    ! w_gpu(:, :, :, 5) -> rho E

    !$cuf kernel do(3) <<<*,*>>> 
    do k = 1,nz
        do j = 1,ny
            do i = 1,nx

                !
                ! du/dx
                !
                dudx = ( &
                    ( &
                        -1 * w_gpu(i+2, j, k, 2) / w_gpu(i+2, j, k, 1) &
                    ) &
                    + &
                    ( &
                        8 * w_gpu(i+1, j, k, 2) / w_gpu(i+1, j, k, 1) &
                    ) &
                    - &
                    ( &
                        8 * w_gpu(i-1, j, k, 2) / w_gpu(i-1, j, k, 1) &
                    ) &
                    + &
                    ( &
                        w_gpu(i-2, j, k, 2) / w_gpu(i-2, j, k, 1) &
                    ) &
                ) &
                ! replace 12 * \Delta x by 6 * (2 \Delta x) since x can vary
                ! with grid size (Without really looking at grid mesh generation
                / ( 6 * ( x_gpu(i+1) - x_gpu(i-1)))

                !
                ! dv/dy
                !
                dvdy = ( &
                    ( &
                        -1 * w_gpu(i, j+2, k, 3) / w_gpu(i, j+2, k, 1) &
                    ) &
                    + &
                    ( &
                        8 * w_gpu(i, j+1, k, 3) / w_gpu(i, j+1, k, 1) &
                    ) &
                    - &
                    ( &
                        8 * w_gpu(i, j-1, k, 3) / w_gpu(i, j-1, k, 1) &
                    ) &
                    + &
                    ( &
                        w_gpu(i, j-2, k, 3) / w_gpu(i, j-2, k, 1) &
                    ) &
                ) &
                ! replace 12 * \Delta y by 6 * (2 \Delta y) since y can vary
                ! with grid size (Without really looking at grid mesh generation
                / ( 6 * ( y_gpu(j+1) - y_gpu(j-1)))

                !
                ! dw/dz
                !
                dwdz = ( &
                    ( &
                        -1 * w_gpu(i, j, k+2, 4) / w_gpu(i, j, k+2, 1) &
                    ) &
                    + &
                    ( &
                        8 * w_gpu(i, j, k+1, 4) / w_gpu(i, j, k+1, 1) &
                    ) &
                    - &
                    ( &
                        8 * w_gpu(i, j, k-1, 4) / w_gpu(i, j, k-1, 1) &
                    ) &
                    + &
                    ( &
                        w_gpu(i, j, k-2, 4) / w_gpu(i, j, k-2, 1) &
                    ) &
                ) &
                ! replace 12 * \Delta z by 6 * (2 \Delta z) since z can vary
                ! with grid size (Without really looking at grid mesh generation
                / ( 6 * ( z_gpu(k+1) - z_gpu(k-1)))

                dx = (x_gpu(i-1) + x_gpu(i+1)) / 2
                dy = (y_gpu(j-1) + y_gpu(j+1)) / 2
                dz = (z_gpu(k-1) + z_gpu(k+1)) / 2

                dissipation_rate = dissipation_rate + &
                    (dudx**2 + dvdy**2 + dwdz**2) * dx * dy * dz
            enddo
        enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()

    dissipation_rate = dissipation_rate / (rlx * rly * rlz)

    ! sum all the values across MPI, store the result in mpi_dissipation_sum
    call MPI_REDUCE(dissipation_rate, mpi_dissipation_sum, 1, mpi_prec,  MPI_SUM, 0, MPI_COMM_WORLD, iermpi)

    ! distribute mpi_dissipation_sum to all MPI procs
    call MPI_BCAST(mpi_dissipation_sum, 1, mpi_prec, 0, MPI_COMM_WORLD, iermpi)

    dissipation_rate = mpi_dissipation_sum
endsubroutine dissipation_calculation
