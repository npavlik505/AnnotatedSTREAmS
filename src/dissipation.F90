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
                ! at i = 1 dx = 0 so the boundaries will be really wonky. 
                ! the same happens for dy at j = 1, and for dz at k = 1
                !
                if (i == 1) then
                    dx = (x_gpu(i+1) - x_gpu(i)) / 2
                else
                    dx = (x_gpu(i+1) - x_gpu(i-1)) / 2
                endif

                if (j == 1) then
                    dy = (y_gpu(j+1) - y_gpu(j)) / 2
                else
                    dy = (y_gpu(j+1) - y_gpu(j-1)) / 2
                endif

                if (k == 1) then
                    dz = (z_gpu(k+1) - z_gpu(k)) / 2
                else
                    dz = (z_gpu(k+1) - z_gpu(k-1)) / 2
                endif

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
                / ( 12 * dx)

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
                / ( 12 * dy)

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
                / ( 12 * dz)

                dissipation_rate = dissipation_rate + &
                    ((dudx**2 + dvdy**2 + dwdz**2) * dx * dy * dz)
            enddo
        enddo
    enddo
    !@cuf iercuda=cudaDeviceSynchronize()

    dissipation_rate = dissipation_rate / (rlx * rly * rlz)

    !write(*, *) "dissipation rate", dissipation_rate

    ! sum all the values across MPI, store the result in mpi_dissipation_sum
    call MPI_REDUCE(dissipation_rate, mpi_dissipation_sum, 1, mpi_prec,  MPI_SUM, 0, MPI_COMM_WORLD, iermpi)

    ! distribute mpi_dissipation_sum to all MPI procs
    call MPI_BCAST(mpi_dissipation_sum, 1, mpi_prec, 0, MPI_COMM_WORLD, iermpi)

    dissipation_rate = mpi_dissipation_sum
endsubroutine dissipation_calculation

subroutine energy_calculation()
    use mod_streams

    implicit none

    integer :: i, j, k

    real*8 :: uu, vv, ww, rho
    real*8 :: dudx, dvdy, dwdz
    real*8 :: dx, dy, dz
    real*8 :: mpi_energy_sum

    ! dissipation_rate is defined in mod_streams, we reset its value here
    energy = 0.0

    mpi_energy_sum = 0.0

    ! w_gpu(:, :, :, 1) -> rho
    ! w_gpu(:, :, :, 2) -> rho u
    ! w_gpu(:, :, :, 3) -> rho v
    ! w_gpu(:, :, :, 4) -> rho w
    ! w_gpu(:, :, :, 5) -> rho E

    !$cuf kernel do(3) <<<*,*>>> 
    do k = 1,nz
        do j = 1,ny
            do i = 1,nx
                rho = w_gpu(i, j, k, 1)
                uu = w_gpu(i, j, k, 2) / rho
                vv = w_gpu(i, j, k, 3) / rho
                ww = w_gpu(i, j, k, 4) / rho

                dx = (x_gpu(i+1) - x_gpu(i-1)) / 2
                dy = (y_gpu(j+1) - y_gpu(j-1)) / 2
                dz = (z_gpu(k+1) - z_gpu(k-1)) / 2

                energy = energy + &
                    (uu**2 + vv**2 + ww**2) * dx * dy * dz
            enddo
        enddo
    enddo

    energy = energy * 0.5

    ! sum all the values across MPI, store the result in mpi_energy_sum 
    call MPI_REDUCE(energy, mpi_energy_sum, 1, mpi_prec,  MPI_SUM, 0, MPI_COMM_WORLD, iermpi)

    ! distribute mpi_energy_sum to all MPI procs
    call MPI_BCAST(mpi_energy_sum, 1, mpi_prec, 0, MPI_COMM_WORLD, iermpi)

    energy = mpi_energy_sum

end subroutine energy_calculation
