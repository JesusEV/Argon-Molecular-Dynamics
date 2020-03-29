!------------------------------------------------------------------------------
!        MODULE Physical Constants`
!------------------------------------------------------------------------------
! MODULE        : physconst
!
! DESCRIPTION:
!> This module holds the needed physical constants used in the simulation.
!------------------------------------------------------------------------------
MODULE physconst
    USE kinds
    IMPLICIT NONE
    REAL(kind=dbl), PARAMETER :: kboltz =    0.0019872067_dbl   ! boltzman constant in kcal/mol/K
    REAL(kind=dbl), PARAMETER :: mvsq2e = 2390.05736153349_dbl  ! m*v^2 in kcal/mol
    REAL(kind=dbl), PARAMETER :: onesigmavel = 0.001765695_dbl  ! kb T/M in kcal/mol
    PRIVATE
    PUBLIC :: kboltz, mvsq2e, onesigmavel
END MODULE physconst






!------------------------------------------------------------------------------
!        MODULE Molecular Dynamics System`
!------------------------------------------------------------------------------
! MODULE        : mdsys
!
! DESCRIPTION:
!> This module holds the complete system information.
!------------------------------------------------------------------------------
MODULE mdsys
    USE kinds
    IMPLICIT NONE
    INTEGER :: natoms, MD_step,nsteps, thermUpdate
    INTEGER :: pair_num, tracking_size = 10
    REAL(kind=dbl) dt, mass, epsilon, sigma, box, rcut
    REAL(kind=dbl) ekin, epot, temp, res_temp
    REAL(kind=dbl), POINTER, DIMENSION (:) :: rx, ry, rz
    REAL(kind=dbl), POINTER, DIMENSION (:) :: vx, vy, vz
    REAL(kind=dbl), POINTER, DIMENSION (:) :: fx, fy, fz
    REAL(kind=dbl), POINTER, DIMENSION (:) :: dists
    REAL(kind=dbl), POINTER, DIMENSION (:) :: temp_series
END MODULE mdsys





!------------------------------------------------------------------------------
!        MODULE Physics`
!------------------------------------------------------------------------------
! MODULE        : physics
!
! DESCRIPTION:
!> This module contains all physically-related routines.
!------------------------------------------------------------------------------
MODULE physics

    USE utils
    USE kinds
    USE mdsys
    USE physconst

    IMPLICIT NONE

    CONTAINS

!---------------------------------------------------------------------------
!> Get kinetic energy routine. Computes total kinetic energy of the system.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE getekin
        INTEGER :: i
        ekin = 0.0_dbl
        DO i=1, natoms
            ekin = ekin + 0.5_dbl * mvsq2e * mass * (vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i))
        END DO
    END SUBROUTINE getekin



!---------------------------------------------------------------------------
!> Get temperature routine. Computes the temperature of the system.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE gettemp
        temp = 2.0_dbl * ekin/(3.0_dbl*DBLE(natoms-1))/kboltz
        temp_series(MD_step) = temp
    END SUBROUTINE gettemp



!---------------------------------------------------------------------------
!> Force routine. Computes the temperature of the system.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE force
        REAL(kind=dbl) :: r_sq, ffac, dx, dy, dz
        REAL(kind=dbl) :: rcutsq, r6, c6, c12, rinv
        INTEGER :: i, j

        epot=0.0_dbl
        CALL force_to_zero

        rcutsq=rcut*rcut
        c12 = 4.0_dbl*epsilon*sigma**12
        c6  = 4.0_dbl*epsilon*sigma**6

        DO i=1, natoms-1
             DO j=i+1, natoms
                    
                dx=pbc(rx(i) - rx(j), box)
                dy=pbc(ry(i) - ry(j), box)
                dz=pbc(rz(i) - rz(j), box)
                r_sq = dx*dx + dy*dy + dz*dz

                ! compute force and energy if within cutoff */
                IF (r_sq < rcutsq) THEN
                   rinv = 1.0_dbl/r_sq
                   r6 = rinv*rinv*rinv
                   ffac = (12.0_dbl*c12*r6 - 6.0_dbl*c6)*r6*rinv
                   epot = epot + r6*(c12*r6 - c6)

                    fx(i) = fx(i) + dx*ffac
                    fy(i) = fy(i) + dy*ffac
                    fz(i) = fz(i) + dz*ffac

                    fx(j) = fx(j) - dx*ffac
                    fy(j) = fy(j) - dy*ffac
                    fz(j) = fz(j) - dz*ffac
                END IF
             END DO
        END DO
    END SUBROUTINE force



!---------------------------------------------------------------------------
!> Velocity Verlet routine. Updates velocities and positions via Verlet Alg.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE velverlet
        INTEGER :: i

        ! first part: propagate velocities by half and positions by full step
        DO i=1, natoms
            vx(i) = vx(i) + 0.5_dbl * dt / mvsq2e * fx(i) / mass
            vy(i) = vy(i) + 0.5_dbl * dt / mvsq2e * fy(i) / mass
            vz(i) = vz(i) + 0.5_dbl * dt / mvsq2e * fz(i) / mass
            rx(i) = rx(i) + dt*vx(i)
            ry(i) = ry(i) + dt*vy(i)
            rz(i) = rz(i) + dt*vz(i)
        END DO

        ! compute forces and potential energy 
        CALL force

        ! second part: propagate velocities by another half step */
        DO i=1, natoms
            vx(i) = vx(i) + 0.5_dbl * dt / mvsq2e * fx(i) / mass
            vy(i) = vy(i) + 0.5_dbl * dt / mvsq2e * fy(i) / mass
            vz(i) = vz(i) + 0.5_dbl * dt / mvsq2e * fz(i) / mass
        END DO
    END SUBROUTINE velverlet



!---------------------------------------------------------------------------
!> Poors-man Thermostat routine. Rescales temperature of the system.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE thermostat
        USE kinds
        USE mdsys
        USE physconst
        IMPLICIT NONE

        INTEGER :: i

        DO i=1, natoms
            vx(i) = vx(i)*sqrt(res_temp/temp)
            vy(i) = vy(i)*sqrt(res_temp/temp)
            vz(i) = vz(i)*sqrt(res_temp/temp)
        END DO
    END SUBROUTINE thermostat



!---------------------------------------------------------------------------
!> Maxwell Bolzmann velocity initializator. Initialize velocities according
!> to MB velocity distribution at a given temperature.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE MaxBoltz_Dist_vel_init
        INTEGER :: i

        DO i=1, natoms
            vx(i) = box_muller_method(onesigmavel,0.0_dbl)
            vy(i) = box_muller_method(onesigmavel,0.0_dbl)
            vz(i) = box_muller_method(onesigmavel,0.0_dbl)
        END DO
    END SUBROUTINE MaxBoltz_Dist_vel_Init



!---------------------------------------------------------------------------
!> FCC lattice positions iniatializator. Initialize positions in a FCC 
!> lattice.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE fcc_lattice_positions_init
        INTEGER :: i,j,k
        INTEGER :: lattice_size, cells_num, cbrt_cells, cell_idx = 1
        REAL(kind=dbl) :: lattice_spacing

        cells_num = natoms/4
        cbrt_cells = INT(cells_num**(1.0/3))
        lattice_spacing = box/cbrt_cells

        DO i=1, cbrt_cells
            DO j=1, cbrt_cells
                DO k=1, cbrt_cells
                    rx(cell_idx + 0) = (i+0.0_dbl)*lattice_spacing  
                    ry(cell_idx + 0) = (j+0.0_dbl)*lattice_spacing 
                    rz(cell_idx + 0) = (k+0.0_dbl)*lattice_spacing 

                    rx(cell_idx + 1) = (i+0.5_dbl)*lattice_spacing  
                    ry(cell_idx + 1) = (j+0.5_dbl)*lattice_spacing 
                    rz(cell_idx + 1) = (k+0.0_dbl)*lattice_spacing 

                    rx(cell_idx + 2) = (i+0.5_dbl)*lattice_spacing  
                    ry(cell_idx + 2) = (j+0.0_dbl)*lattice_spacing 
                    rz(cell_idx + 2) = (k+0.5_dbl)*lattice_spacing 

                    rx(cell_idx + 3) = (i+0.0_dbl)*lattice_spacing  
                    ry(cell_idx + 3) = (j+0.5_dbl)*lattice_spacing 
                    rz(cell_idx + 3) = (k+0.5_dbl)*lattice_spacing 

                    cell_idx = cell_idx + 4
                END DO
            END DO
        END DO
    END SUBROUTINE fcc_lattice_positions_Init



!---------------------------------------------------------------------------
!> Force to Zero routine. Resets Force values.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE force_to_zero
    INTEGER :: i

    DO i=1, natoms
       fx(i) = 0.0_dbl
       fy(i) = 0.0_dbl
       fz(i) = 0.0_dbl
    END DO
    END SUBROUTINE force_to_zero


!---------------------------------------------------------------------------
!> Get distances routine. Computes distances between the particles in the
!> system.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE get_distances
        INTEGER :: i,j,k
        REAL(kind=dbl) :: dx, dy, dz

        k = 1
        DO i=1, tracking_size
            DO j=i+1, natoms
                dx=pbc(rx(i) - rx(j), box)
                dy=pbc(ry(i) - ry(j), box)
                dz=pbc(rz(i) - rz(j), box)
                dists(k) = SQRT(dx*dx + dy*dy + dz*dz)
                k=k+1
            END DO
        END DO
    END SUBROUTINE get_distances

END MODULE physics

































