subroutine wrap_bc()
    use mod_streams

    call bc(1)
end subroutine

subroutine wrap_startmpi()
    call startmpi()
end subroutine

subroutine wrap_setup()
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
