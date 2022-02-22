
subroutine bcblow(ilat)
    use mod_streams
    implicit none

    integer :: ilat, start_x
    integer :: i,j,k
    real(mykind) :: jet_velocity, rho

    j = 1
    jet_velocity = 1.0

    ! if on the bottom
    if (ilat == 3) then
        call bcwall_pl(ilat)

        ! TODO: make this selection more dynamic
        if (nrank == 0) then
            start_x = nx / 2
            ! w_gpu(:, :, :, 1) -> rho
            ! w_gpu(:, :, :, 2) -> rho u
            ! w_gpu(:, :, :, 3) -> rho v
            ! w_gpu(:, :, :, 4) -> rho w
            ! w_gpu(:, :, :, 5) -> rho E

            !$cuf kernel do(2) <<<*,*>>>
            do i = start_x,nx
                do k = 1,nz
                    ! do stuff here
                    rho = w_gpu(i,j,k,1)
                    w_gpu(i,j,k,3) = rho * jet_velocity
                enddo
            enddo
            !@cuf iercuda=cudaDeviceSynchronize()
        endif
    else 
        ! TODO: error - we should not be here
    endif


end subroutine bcblow
