!
subroutine hcount_int(n, arr, item, times)
    implicit None
    integer :: n, item, times, i
    integer, dimension(n) :: arr

    times = 0
    do i = 1, n
        if (arr(i) == item) times = times + 1
    end do

    return
end subroutine hcount_int