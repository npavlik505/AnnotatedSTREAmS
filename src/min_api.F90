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
