!
!
!

subroutine Fillbond (na, iza, conn, bondneed, fakebmx)
    implicit None

    integer :: na, i, j, icycle, ineed, bondadd
    integer :: startp, chainp
    integer, dimension(na) :: iza, coord, bondneed, needlist, possible, aval
    integer, dimension(na, na) :: conn, fakebmx
    integer, dimension(4, na) :: neighbor

    coord = 0
    fakebmx = 0

    do i = 1, na
        if (iza(i) > 18) then
            bondneed(i) = 0
            cycle
        end if
        bondneed(i) = min(abs(iza(i) - 10), abs(iza(i)-2), abs(iza(i)-18))
    end do

    do i = 1, na
        do j = i + 1, na
            if (conn(i,j) > 0 .and. bondneed(j) > 0) then
                coord(i) = coord(i) + 1
                coord(j) = coord(j) + 1
                fakebmx(i, j) = fakebmx(i, j) + 1
                fakebmx(j, i) = fakebmx(j, i) + 1
            end if
        end do
        bondneed(i) = bondneed(i) - coord(i)
    end do

    do icycle = 1, 10
        ineed = 0
        possible = 0
        needlist = 0
        aval = 1
        do i = 1, na
            if (bondneed(i) <= 0 ) cycle
            ineed = ineed + 1
            needlist(ineed) = i
            do j = i + 1, na
                if (conn(i, j) > 0 .and. bondneed(j) > 0) then
                    possible(i) = possible(i) + 1
                    possible(j) = possible(j) + 1
                end if
            end do
        end do

        do i = 1, ineed - 1
            do j = i + 1, ineed
                if (possible(needlist(i)) > possible(needlist(j))) then
                    call change_atom(needlist(i), needlist(j))
                end if
            end do
        end do

        chainp = ineed
        do i = 1, ineed
            if (possible(needlist(i)) > 1) then
                chainp = i
                exit
            end if
        end do

        do i = 1, ineed - 1
            if ( i >= chainp) chainp = chainp + 1
            startp = chainp
            do j = startp, ineed
                if (conn(needlist(i), needlist(j)) > 0) then
                    call change_atom(needlist(chainp), needlist(j))
                    chainp = chainp + 1
                end if
            end do
        end do

        bondadd = 0
        do i =1, ineed
            do j = i + 1, ineed
                if ( conn(needlist(i), needlist(j)) > 0 .and. aval(needlist(j)) > 0 .and. aval(needlist(i)) > 0) then
                    bondneed(needlist(i)) = bondneed(needlist(i)) - 1
                    bondneed(needlist(j)) = bondneed(needlist(j)) - 1
                    aval(needlist(i)) = 0
                    aval(needlist(j)) = 0
                    bondadd = bondadd + 1
                    fakebmx(needlist(i), needlist(j)) = fakebmx(needlist(i), needlist(j)) + 1
                    fakebmx(needlist(j), needlist(i)) = fakebmx(needlist(j), needlist(i)) + 1

                    exit
                end if
            end do
        end do

        if (ineed <= 1 .or. bondadd ==0 ) exit
    end do
end subroutine Fillbond

subroutine FillBondSurface (na, iza, conn, bondneed, surface, fakebmx)
    !
    !
    ! Input parameters:
    !   na: number of atoms
    !   iza: atomic number
    !   conn: connection matrix
    !   bondneed: number of bonds needed
    !   surface: surface atoms
    !   fakebmx: fake bond matrix
    !
    ! Other parameters:
    !   icycle: loop index
    !   ineed: number of atoms needed
    !   itail: tail index
    !   startp: start index
    !   chainp: chain index
    !   i, j: loop index
    !   possible: possible bonds
    !   needlist: needed list
    !   coord: coordination number
    !   aval: available bonds
    !   bondadd: number of bonds added
    !

    implicit None

    integer :: na, bondadd
    integer :: i, j, icycle, ineed, itail, startp, chainp
    integer, dimension(na) :: iza, bondneed, surface
    integer, dimension(na) :: needlist, coord, possible, aval
    integer, dimension(na, na) :: conn, fakebmx
    !integer, dimension(4, na) :: neighbor

    ! initialize the bondneed
    ! for atom > 18, bondneed = 0
    ! for atom < 18, bondneed equles to the minimum index difference 
    !   to nearst noble gas
    !   eg, for H, bondneed = 1; for O, bondneed = 2; for C, bondneed = 4.
    coord = 0
    do i = 1, na
        if (iza(i) > 18) then
            bondneed(i) = 0
            cycle
        end if
        bondneed(i) = min(abs(iza(i) - 10), abs(iza(i) - 2), abs(iza(i) - 18))
    end do

    ! initial coordination number, bond matrix, and update bondneed (bondneed - coord)
    !     for atom j, the bondneed should large than 0
    do i = 1, na
        do j = i + 1, na
            if (conn(i, j) > 0 .and. bondneed(j) >0) then
                coord(i) = coord(i) + 1
                coord(j) = coord(j) + 1
                fakebmx(i, j) = fakebmx(i, j) + 1
                fakebmx(j, i) = fakebmx(j, i) + 1
            end if
        end do
        bondneed(i) = bondneed(i) - coord(i)
    end do


    do icycle = 1, 9999
        ineed = 0
        possible = 0
        needlist = 0
        aval = 1
        bondadd = 0

        ! inner loop 1, only consider the atom > 18
        !   possible means the extra bond needed (bondneed-bonded) for each atom
        do i = 1, na
            if (bondneed(i) <= 0) cycle
            ineed = ineed + 1
            needlist(ineed) = i
            do j = i + 1, na
                if (conn(i, j) > 0 .and. bondneed(j) > 0) then
                    possible(i) = possible(i) + 1
                    possible(j) = possible(j) + 1
                end if
            end do

            ! for test:
            ! print *, ineed, needlist
        end do


        ! inner loop 2
        !   resort the needlist to make the surface atoms to the end of the needlist,
        !   because the surface atoms maybe fulfill bondneed
        itail = ineed
        do i = ineed, 1, -1
            if (surface(needlist(i)) == 1) then
                call change_atom(needlist(i), needlist(itail))
                itail = itail -1
            end if
        end do

        ! inner loop 3
        !    resort the need list to make the more bond number with slab to the end
        do i = 1, ineed -1
            do j = i + 1, ineed
                if (possible(needlist(i)) >  possible(needlist(j)) .and. (surface(needlist(i)) + surface(needlist(j)))/=1) then
                    call change_atom(needlist(i), needlist(j))
                end if
            end do
        end do

        ! update the chainp to the atom index needed to be added bond
        chainp = ineed
        do i = 1, ineed
            if (possible(needlist(i)) > 1) then
                chainp = i
                exit
            end if
        end do

        ! inner loop 4
        !    add the bond to the atom in needlist
        if (chainp /= 1) then
            do i = 1, chainp -1
                ! exclude the surface atom, because the surface atom maybe fulfill the requirement
                if (surface(needlist(i)) == 1) cycle
                ! give the i and j add bond and update corresponding matrix
                do j = i + 1, ineed
                    if (conn(needlist(i), needlist(j)) > 0 .and. aval(needlist(j)) > 0 .and. aval(needlist(j)) > 0) then
                        bondneed(needlist(i)) = bondneed(needlist(i)) - 1
                        bondneed(needlist(j)) = bondneed(needlist(j)) - 1
                        aval(needlist(i)) = 0
                        aval(needlist(j)) = 0
                        bondadd = bondadd + 1
                        fakebmx(needlist(i), needlist(j)) = fakebmx(needlist(i), needlist(j)) + 1
                        fakebmx(needlist(j), needlist(i)) = fakebmx(needlist(j), needlist(i)) + 1
                        exit
                    end if
                end do
            end do
            if (bondadd > 0) cycle
        end if

        ! if no bond added before,
        do i = 1, ineed - 1
            if (i >= chainp) chainp = chainp + 1
            startp = chainp
            do j = startp, ineed
                if (conn(needlist(i), needlist(j)) > 0 .and. surface(needlist(j)) == 0) then
                    call change_atom(needlist(chainp), needlist(j))
                    chainp = chainp + 1
                end if
            end do
        end do

        do i = 1, ineed
            if (surface(needlist(i)) == 1) cycle
            do j = i + 1, ineed
                if (conn(needlist(i), needlist(j)) > 0 .and. aval(needlist(j)) > 0 .and. aval(needlist(i)) > 0) then
                    bondneed(needlist(i)) = bondneed(needlist(i)) - 1
                    bondneed(needlist(j)) = bondneed(needlist(j)) - 1
                    aval(needlist(i)) = 0
                    aval(needlist(j)) = 0
                    bondadd = bondadd + 1
                    fakebmx(needlist(i), needlist(j)) = fakebmx(needlist(i), needlist(j)) + 1
                    fakebmx(needlist(j), needlist(i)) = fakebmx(needlist(j), needlist(i)) + 1
                    exit
                end if
            end do
        end do

        if (bondadd == 0) then
            if (ineed <= 1) then
                exit
            else
                do i = 1, ineed
                    do j = i + 1, ineed
                        if (conn(needlist(i), needlist(j)) > 0 .and. aval(needlist(i)) > 0 .and. aval(needlist(j)) > 0) then
                            bondneed(needlist(i)) = bondneed(needlist(i)) - 1
                            bondneed(needlist(j)) = bondneed(needlist(j)) - 1
                            aval(needlist(i)) = 0
                            aval(needlist(j)) = 0
                            bondadd = bondadd + 1
                            fakebmx(needlist(i), needlist(j)) = fakebmx(needlist(i), needlist(j)) + 1
                            fakebmx(needlist(j), needlist(i)) = fakebmx(needlist(j), needlist(i)) + 1
                            exit
                        end if
                    end do
                end do
                if ( bondadd == 0) exit
            end if
        end if
    end do
end subroutine FillBondSurface 

subroutine JudgeBond (na, iza, xa, cell, minstr, fakebmx, bondneed)
    implicit None

    integer :: na
    integer :: i
    integer, dimension(na) :: iza, bondneed
    integer, dimension(na) :: surface
    integer, dimension(na, na) :: conn, fakebmx
    double precision, dimension(3, 3) :: cell
    double precision, dimension(3, na) :: xa
    logical :: minstr, Lsurface

    minstr = .true.
    call get_conn(na, iza, xa, conn, cell)
    call Fillbond(na, iza, conn, bondneed, fakebmx)

    do i = 1, na
        if (bondneed(i) /= 0) then
            minstr = .false.
            exit
        end if
    end do

end subroutine JudgeBond

subroutine JudgeBondSurface(na, iza, xa, cell, minstr, fakebmx, bondneed, surface)
    implicit None

    integer :: na, i
    integer, dimension(na) :: iza, bondneed, surface
    integer, dimension(na, na) :: conn, fakebmx
    double precision, dimension(3, 3) :: cell
    double precision, dimension(3, na) :: xa
    logical :: minstr, Lsurface
    
    !print  *, 'test start'
    minstr = .true.

    call ScreenSurface(na, iza, xa, cell, surface)
    !print *, 'screen surface done, surface is', surface
    call get_conn(na, iza, xa, conn, cell)
    call FillBondSurface(na, iza, conn, bondneed, surface, fakebmx)


    do i = 1, na
        if (bondneed(i) /= 0) then
            if (surface(i) == 1) bondneed(i) = 0
        end if
    end do

    do i = 1, na
        if (bondneed(i) /= 0 ) then
            minstr = .false.
            exit
        end if
    end do
    
end subroutine JudgeBondSurface

subroutine ScreenSurface (na, iza, xa, cell, surface)
    !
    ! Find the surface adsorbed atoms
    !
    ! Input parameters: 
    !   na: number of atoms
    !   iza: atomic number
    !   xa: atomic coordinates
    !   cell: cell parameters
    !   surface: surface atoms
    !
    ! Other parameters:
    !   totbond: total bond length
    !   dist: distance between two atoms
    !   co1, co2: covalent radius
    !   i, j: loop index

    implicit None

    double precision :: totbond, dist, co1, co2
    integer :: na
    integer :: i, j
    integer, dimension(na) :: iza, surface
    double precision, dimension(3, 3) :: cell
    double precision, dimension(3, na) :: xa

    ! exclude the atoms > 18 (only surface adsorbate)
    ! bonded to atoms > 18 (metal atoms) are consider into surface.
    !print *, '    start screen surface'
    surface = 0
    do i = 1, na
        if (iza(i) > 18) cycle
        do j = 1, na
            if (iza(j) > 18) then
                call get_dist(na, dist, i, j, xa, cell)
                !print *, dist

                call radius(iza(i), co1)
                call radius(iza(j), co2)
                totbond = co1 + co2
                !print *, totbond

                if (dist < totbond + 0.4d0) then
                    surface(i) = 1
                    exit
                end if
            end if
        end do
    end do
end subroutine ScreenSurface

subroutine change_atom (n1, n2)
    implicit None
    integer :: t, n1, n2

    t  = n1
    n1 = n2
    n2 = t

    return
end subroutine change_atom

 