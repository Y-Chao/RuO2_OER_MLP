!
subroutine hcalbondangle (fi, fj, fk, lat, theta)
    implicit None
    double precision, dimension(3) :: fi, fj, fk
    double precision, dimension(3) :: rij, rjk, rik
    double precision :: dij, djk, dik
    double precision, dimension(3,3) :: lat
    double precision :: theta, costheta

    ! fi, fi, fk represents the fractional coordinates
    rij = fj - fi
    rik = fk - fi

    ! the fraction vector 
    rij = rij - dnint(rij) 
    rik = rik - dnint(rik)

    rjk = rik - rij

    dij = dsqrt(sum(matmul (lat, rij))**2.0)
    dik = dsqrt(sum(matmul (lat, rik))**2.0)
    djk = dsqrt(sum(matmul (lat, rjk))**2.0)

    costheta = (dij**2.0 + dik**2.0 - djk**2.0) / (2.0 * dij * dik)
    theta = dacos(costheta) / 3.1415926 * 180.0

    return
end subroutine hcalbondangle
