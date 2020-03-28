MODULE kinds
    IMPLICIT NONE
    INTEGER, PARAMETER :: dbl = selected_real_kind(14,200)  ! double precision floating point
    INTEGER, PARAMETER :: sgl = selected_real_kind(6,30)    ! single precision floating point
    INTEGER, PARAMETER :: sln = 200                         ! length of I/O input line
    PRIVATE
    PUBLIC :: sgl, dbl, sln
END MODULE kinds


MODULE utils
    USE kinds
    IMPLICIT NONE
    REAL(kind=dbl), PARAMETER :: pi = 3.14159265359_dbl   ! pi

    PRIVATE
    PUBLIC :: pbc, box_muller_method, pi

    CONTAINS
    
    ! helper function: apply minimum image convention 
    FUNCTION pbc(x, box)
        REAL(kind=dbl), INTENT(IN)  :: x, box
        REAL(kind=dbl) :: pbc
        pbc = x - box*(ANINT(x/box))
    END FUNCTION pbc


    REAL(8) FUNCTION box_muller_method(sigma, mu) RESULT(res)
        REAL(kind=dbl), INTENT(in) :: sigma, mu
        REAL(kind=dbl) :: R, thetha
        REAL(kind=dbl) :: u, up
        u = RAND()
        up = RAND()
        R = SQRT(-2*sigma**2 * LOG(1-u))
        thetha = 2*pi*up
        res = R*SIN(thetha) + mu
    END FUNCTION box_muller_method


END MODULE utils
