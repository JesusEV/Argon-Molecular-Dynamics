MODULE random
    IMPLICIT NONE

    INTEGER :: seed

    CONTAINS

!---------------------------------------------------

    REAL(8) FUNCTION Maxwell_Boltzmann_Dist(T, m, kb) RESULT(res)
        REAL(8), INTENT(in) :: T, m, kb
        REAL(8) :: sigma
        REAL(8) :: mu = 0.0
        REAL(8) :: vx, vy, vz

        sigma = SQRT(kb*T/m)

        vx = box_muller_method(sigma,mu)
        vy = box_muller_method(sigma,mu)
        vz = box_muller_method(sigma,mu)

        res = SQRT(vx*vx + vy*vy + vz*vz)

    END FUNCTION Maxwell_Boltzmann_Dist


!---------------------------------------------------

    REAL(8) FUNCTION box_muller_method(sigma, mu) RESULT(res)
        REAL(8), INTENT(in) :: sigma, mu
        REAL(8) :: R, thetha, pi
        REAL(8) :: u, up

        u = RAND()
        up = RAND()

        pi = ACOS(-1.0)        

        R = SQRT(-2*sigma**2 * LOG(1-u))
        thetha = 2*pi*up

        res = R*SIN(thetha) + mu
    END FUNCTION box_muller_method

END MODULE random
