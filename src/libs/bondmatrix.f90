! 
! calculate bond matrix
!
subroutine hRadicalOrder (na, iza, bmatrix, radical, fa, lat)
    implicit None
    integer :: na, maxbond, allbond, countatom
    integer, dimension(na) :: iza, radical, bondorder, coord
    integer, dimension(na, na) :: bmatrix
    integer :: i, iatm, junc(2), nc, jatm
    double precision, dimension(3, na) :: fa
    double precision :: lat, theta

    call calcoordnumber (na, bmatrix, coord, bondorder)

    do iatm = 1, na
        ! select different method for different element
        select case(iza(iatm) - 1)
        case (1)
            ! For H atom
            radical(iatm) = abs(coord(iatm) - 1)
        case (6)
            ! For C atom
            maxbond = maxval(bmatrix(iatm, :))
            if (maxbond == 2) then
                radical(iatm) = abs(coord(iatm) - 4)
            else if (maxbond == 2) then
                if (coord(iatm) == 2) then
                    nc = 0
                    junc = (/0, 0/)
                    do i = 1, na
                        if (bmatrix(iatm, i) == 2) then
                            nc = nc + 1
                            junc(nc) = i
                        end if
                    end do
                    if (junc(2) == 0) then
                        radical(iatm) = abs(coord(iatm) - 3)
                    else
                        call hcalbondangle (fa(:, iatm), fa(:, junc(1)), fa(:, junc(2)), lat, theta)

                        if (theta > 150.0) then
                            radical(iatm) = 0
                        else
                            radical(iatm) = abs(coord(iatm) - 3)
                        end if
                    end if
                else if (coord(iatm) == 3) then
                    if (bondorder(iatm) == 5) then
                        countatom = 0
                        do jatm = 1, na
                            if (bmatrix(iatm, jatm) > 0 .and. bondorder(jatm)==5) countatom = countatom + 1
                        end do
                        radical(iatm) = abs(countatom - 2)
                    else
                        radical(iatm) = abs(bondorder(iatm) - 4)
                    end if
                else
                    radical(iatm) = abs(coord(iatm) - 3)
                end if
            else if (maxbond == 3) then
                radical(iatm) = abs(coord(iatm) - 2)
            else
                radical(iatm) = 4
            end if
        case (7)
            ! For N atom
            maxbond = maxval(bmatrix(iatm, :))
            if (maxbond == 1) then
                radical(iatm) = abs(coord(iatm) - 3)
            else if (maxbond == 2) then
                radical(iatm) = abs(coord(iatm) - 2)
            else if (maxbond == 3) then
                radical(iatm) = abs(coord(iatm) - 1)
            else
                radical(iatm) = 3
            end if
        case (8)
            ! For O atom
            maxbond = maxval(bmatrix(iatm, :))
            if (maxbond == 1) then
                radical(iatm) = abs(coord(iatm) - 2)
            else if (maxbond == 2) then
                radical(iatm) = abs(coord(iatm) - 1)
            else
                radical(iatm) = 2
            end if
        end select
    end do
    return
end subroutine hRadicalOrder

subroutine calcoordnumber (na, bmatrix, coord, bondorder)
    implicit None
    integer :: na, i, j
    integer, dimension(na) :: coord, bondorder
    integer, dimension(na, na) :: bmatrix

    do i = 1, na
        coord(i) = 0
        bondorder(i) = sum(bmatrix(i, :))
        do j = 1, na
            if (bmatrix(i, j) > 0) coord(i) = coord(i) + 1
        end do
    end do
    return
end subroutine calcoordnumber

subroutine hbondmatrix (na, fa, iza, lat, bmatrix)
    implicit None
    integer :: i, j, maxbond, idmax, na
    integer, dimension(na) ::iza, conjugate, ctemp
    double precision :: d
    double precision, dimension(3, na) :: fa
    double precision, dimension(3) ::delt
    double precision, dimension(3, 3) :: lat
    integer, dimension(na, na) :: bmatrix
    
    bmatrix = 0
    do i = 1, na-1
        do j = i+1, na
            if (iza(i) > 18 .or. iza(j) > 18) then
                bmatrix(i, j) = 0
                bmatrix(j, i) = 0
                cycle
            end if

            delt = fa(:,i) - fa(:,j)
            delt = delt - dnint(delt)
            delt = matmul(lat, delt)
            d = dsqrt(sum(delt**2.0d0))

            if (iza(i) <= iza(j)) then
                call bond_order(iza(i), iza(j), d, bmatrix(i, j))
                bmatrix(j, i) = bmatrix(i, j)
            else
                call bond_order(iza(j), iza(i), d, bmatrix(i, j))
                bmatrix(j, i) = bmatrix(i, j)
            end if
        end do
    end do
    return
    print *, 'bond matrix'
end subroutine hbondmatrix


subroutine bond_order (n1, n2, d, bij)
    implicit None
    integer :: n1, n2, bij
    double precision :: d
    logical, external :: bsta

    bij = 0
    select case(n1)
        case (1)
            ! For H atom
            select case (n2)
                case (1)
                    if (d < 0.85d0) bij = 1
                case (5)
                    if (d < 1.3d0) bij = 1
                case (6)
                    if (d < 1.2d0) bij = 1
                case (7)
                    if (d < 1.15d0 ) bij = 1
                case(8)
                    if (d < 1.1d0 ) bij = 1
                case(9)
                    if (d < 1.05d0 ) bij = 1
                case(14)
                    if (d < 1.58d0 ) bij = 1
                case(15)
                    if (d < 1.55d0 ) bij = 1
                case(16)
                    if (d < 1.45d0 ) bij = 1
                case(17)
                    if (d < 1.4d0 ) bij = 1
                case(32)
                    if (d < 1.64d0 ) bij = 1
                case(34)
                    if (d < 1.56d0 ) bij = 1
                case(35)
                    if (d < 1.45d0 ) bij = 1
                case(50)
                    if (d < 1.8d0 ) bij = 1
                case(52)
                    if (d < 1.8d0 ) bij = 1
                case(53)
                    if (d < 1.71d0 ) bij = 1
            end select
        case (5)
            ! For B atom
            select case(n2)
                case(5)
                    if (d < 2.0d0 )  bij = 1
                case(6)
                    if (d < 1.66d0 ) bij = 1
                case(17)
                    if (d < 1.85d0 ) bij = 1
            end select
        case(6)
            ! For C atom
            select case(n2)
                case(6)
                    if (d < 1.25d0 ) then
                        bij = 3
                    else if ( d < 1.43d0) then
                        bij = 2
                    else if ( d < 1.64) then
                        bij = 1
                end if
                case(7)
                    if (d < 1.26d0 ) then
                        bij = 3
                    else if ( d < 1.40d0) then
                        bij = 2
                    else if ( d < 1.6) then
                        bij = 1
                    end if
                case(8)
                    if (d<1.15  ) then
                        bij = 3
                    else if ( d < 1.3d0) then
                        bij = 2
                    else if ( d < 1.55d0) then
                        bij = 1
                    end if
                case(9)
                    if (d < 1.48d0 ) bij = 1
                case(14)
                    if (d < 1.96d0 ) bij = 1
                case(15)
                    if (d < 1.97d0 ) bij = 1
                case(16)
                    if (d < 1.7d0 ) then
                        bij = 2
                    else if ( d < 1.92d0) then
                        bij = 1
                    end if
                case(17)
                    if (d < 2.1d0 ) bij = 1
                case(32)
                    if (d < 2.05d0 ) bij = 1
                case(35)
                    if (d < 2.04d0 ) bij = 1
                case(53)
                    if (d < 2.25d0 ) bij = 1
          end select
        case(7)
            ! For N atom
            select case(n2)
                case(7)
                    if (d < 1.2d0 ) then
                        bij = 3
                    else if ( d < 1.35d0) then
                        bij = 2
                    else if ( d < 1.55) then
                        bij = 1
                    end if
                case(8)
                    if (d < 1.15d0 ) then
                        bij = 3
                    else if ( d < 1.3d0) then
                        bij = 2
                    else if ( d < 1.56d0) then
                        bij = 1
                    end if
            end select
        case(8)
            ! For O atom 
            select case(n2)
                case(8)
                    if (d < 1.3d0 ) then
                        bij = 2
                    else if ( d < 1.58d0) then
                        bij = 1
                    end if
                case(14)
                    if (d < 1.73d0 ) bij = 1
                case(15)
                if (d < 1.48d0 ) then
                    bij = 2
                else if ( d < 1.73d0) then
                    bij = 1
                end if
                case(16)
                    if (d < 1.53d0 ) bij = 1
            end select
      end select
    return
end subroutine bond_order

subroutine hsegmentMol (na, fa, iza, lat, group, bmatrix)
    !
    ! This subroutine is used to segment the molecule one by one 
    !
    ! Input parameters:
    !   na : the number of atoms in one str
    !   fa : the fractional coordinates of str
    !   iza : the atomic number of each atom in one str
    !   lat ： the lattice in one str
    !   group : the group number of each atom in one str
    !   bmatrix : the bond matrix of one str
    !
    ! Other parameters:
    !   gid : the group number of the one atom
    !   iatm : the index of the atom

    implicit None
    integer :: gid
    integer :: na
    integer :: iatm
    integer, dimension(na) :: iza
    integer, dimension(na) :: group
    integer, dimension(na, na) :: bmatrix
    double precision, dimension(3, na) :: fa
    double precision, dimension(3, 3) :: lat

    gid = 1
    group = 0

    do iatm =1, na
        ! for non-bond atom or atom with atomic number > 18, the group would labelled one by one
        ! else the start atom would be judgged with all atoms and labeld with gid
        if ((sum(bmatrix(:, iatm)) == 0) .or. iza(iatm) > 18) then
            group(iatm) = gid
            gid = gid + 1
        else
            if (group(iatm) == 0) then
                call hconnectSeg (na, bmatrix, iatm, group, gid)
                gid = gid + 1
            end if
        end if
    end do
end subroutine hsegmentMol

subroutine hsegmentMol2 (na, fa, iza, lat, group, bmatrix)
    !
    ! This subroutine is used to check the group again
    !
    ! Input parameters:
    !   na : the number of atoms in one str
    !   fa : the fractional coordinates of str
    !   iza : the atomic number of each atom in one str
    !   lat ： the lattice in one str
    !   group : the group number of each atom in one str
    !   bmatrix : the bond matrix of one str
    !
    ! Other parameters:
    !   gid : the group number of the one atom
    !  atmid: the index of the atom
    !   iatm: the index of the atom
    implicit None
    integer :: na
    integer :: gid, atmid, iatm
    integer, dimension(na) :: iza
    integer, dimension(na) :: group
    integer, dimension(na, na) :: bmatrix
    double precision, dimension(3, na) :: fa
    double precision, dimension(3, 3) :: lat

    gid = 1
    do iatm = 1, na
        call hjudge_atom_ingroup (na, fa, iza, lat, group, bmatrix, gid, iatm)
    end do

    return
end subroutine hsegmentMol2

recursive subroutine hjudge_atom_ingroup (na, fa, iza, lat, group, bmatrix, gid, atmid)
    !
    ! This subroutine is used to check the group id again
    !
    ! Input parameters:
    !   na : the number of atoms in one str
    !   fa : the fractional coordinates of str
    !   iza : the atomic number of each atom in one str
    !   lat ： the lattice in one str
    !   group : the group number of each atom in one str
    !   bmatrix : the bond matrix of one str
    !   gid : the group number of the one atom
    !   atmid : the index of the atom
    !
    ! Other parameters:
    !   iatm : the index of the atom

    implicit None
    integer :: na, gid, atmid
    integer :: iatm
    integer, dimension(na) :: iza, group
    integer, dimension(na, na) :: bmatrix
    double precision, dimension(3, na) :: fa
    double precision, dimension(3, 3) :: lat

    if (group(atmid) == 0) then
        do iatm = 1, na
            if (bmatrix(iatm, atmid) > 0 .and. group(iatm) > 0) then
                group(atmid) = group(iatm); exit
            end if
        end do
        
        group(atmid) = gid
        gid = gid + 1
    end if
    
    do iatm = 1, na
        if (bmatrix(iatm, atmid) >0 .and. group(iatm)/=group(atmid)) then
            group(iatm) = group(atmid)
            call hjudge_atom_ingroup (na, fa, iza, lat, group, bmatrix, gid, iatm)
        end if
    end do
end subroutine hjudge_atom_ingroup

recursive subroutine hconnectSeg (na, bmatrix, startatom, group, gid)
    !
    !
    ! Input parameters:
    !   na : the number of atoms in one str
    !   bmatrix : the bond matrix of one str
    !   startatom : the start atom to connect
    !   group : the group number of each atom in one str
    !   gid : the group number of the start atom
    !
    ! Other parameters:
    !   i : the index for loop
    !   flag: the flag to judge whether the atom is in the group

    implicit None
    integer :: i
    integer :: na, startatom, gid
    integer, dimension(na, na) :: bmatrix
    integer, dimension(na) :: group
    logical :: flag

    do i = 1, na
        ! if the start atom is bond with atom i and atom i is not in any group
        ! then the atom i would be labeled as gid
        if (bmatrix(startatom, i) > 0 .and. group(i) == 0) then
            group(i) = gid
            call hconnectSeg (na, bmatrix, i, group, gid)
        end if
    end do

    return
end subroutine hconnectSeg



function bsta (d, sd) 
    implicit None
    logical :: bsta
    double precision :: d, sd
    double precision :: tol = 0.1d0

    bsta = .false.
    if (d > (sd-tol) .and. d<(sd+tol)) bsta = .true.

    return
end function bsta