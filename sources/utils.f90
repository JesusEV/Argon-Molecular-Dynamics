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

  PRIVATE
  PUBLIC :: azzero, pbc

CONTAINS
! helper function: zero out an array 
  SUBROUTINE azzero(d, n)
    REAL(kind=dbl), DIMENSION(:), INTENT(INOUT) :: d
    INTEGER, INTENT(IN) :: n
    INTEGER :: i

    DO i=1, n
       d(i) = 0.0_dbl
    END DO
  END SUBROUTINE azzero
    
! helper function: apply minimum image convention 
  FUNCTION pbc(x, box)
    REAL(kind=dbl), INTENT(IN)  :: x, box
    REAL(kind=dbl) :: pbc

    pbc = x - box*(ANINT(x/box))
  END FUNCTION pbc
END MODULE utils

! module to hold the complete system information 
MODULE mdsys
  USE kinds
  IMPLICIT NONE
  INTEGER :: natoms,nfi,nsteps
  REAL(kind=dbl) dt, mass, epsilon, sigma, box, rcut
  REAL(kind=dbl) ekin, epot, temp
  REAL(kind=dbl), POINTER, DIMENSION (:) :: rx, ry, rz
  REAL(kind=dbl), POINTER, DIMENSION (:) :: vx, vy, vz
  REAL(kind=dbl), POINTER, DIMENSION (:) :: fx, fy, fz
END MODULE mdsys
