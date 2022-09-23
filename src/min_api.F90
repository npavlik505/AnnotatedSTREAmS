module othersave
    use mod_streams
    integer :: anotherx
end module

subroutine wrap_bc()
    use mod_streams

    call bc(1)
end subroutine
