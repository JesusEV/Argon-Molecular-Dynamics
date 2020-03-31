!------------------------------------------------------------------------------
!        MODULE Physical Constants`
!------------------------------------------------------------------------------
! MODULE        : Phys_consts
!
! DESCRIPTION:
!> This module holds the needed physical constants used in the simulation.
!------------------------------------------------------------------------------
MODULE Phys_consts
    USE kinds
    IMPLICIT NONE
    REAL(kind=dbl), PARAMETER :: kb =    0.0019872067_dbl       ! kb [kcal k/mol]
    REAL(kind=dbl), PARAMETER :: mv2convfac = 2390.05736153349_dbl  ! mv^2 [kcal / mol]
    REAL(kind=dbl), PARAMETER :: onesigmavel = 0.001765695_dbl  ! vel one sigma [T kb/M]
    PRIVATE
    PUBLIC :: kb, mv2convfac, onesigmavel
END MODULE Phys_consts






!------------------------------------------------------------------------------
!        MODULE Molecular Dynamics System`
!------------------------------------------------------------------------------
! MODULE        : MD_system
!
! DESCRIPTION:
!> This module holds the complete system information.
!------------------------------------------------------------------------------
MODULE MD_system
    USE kinds
    IMPLICIT NONE
    INTEGER :: N, MD_step,MD_steps, thermUpdate
    INTEGER :: pair_num, tracking_size = 10
    REAL(kind=dbl) dt, M, epsilon, sigma, L, r_cut
    REAL(kind=dbl) sq_vel, ekin, epot, temp, res_temp
    REAL(kind=dbl), POINTER, DIMENSION (:) :: rx, ry, rz
    REAL(kind=dbl), POINTER, DIMENSION (:) :: vx, vy, vz
    REAL(kind=dbl), POINTER, DIMENSION (:) :: fx, fy, fz
    REAL(kind=dbl), POINTER, DIMENSION (:) :: dists
    REAL(kind=dbl), POINTER, DIMENSION (:) :: temp_series
    REAL(kind=dbl), POINTER, DIMENSION (:) :: vel_series
END MODULE MD_system





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
    USE MD_system
    USE Phys_consts

    IMPLICIT NONE

    CONTAINS

!---------------------------------------------------------------------------
!> Get kinetic energy routine. Computes total kinetic energy of the system.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE getekin
        INTEGER :: i
        ekin = 0.0_dbl
        DO i=1, N
            sq_vel = (vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i))
            vel_series((MD_step-1)*N+i) = sq_vel
            ekin = ekin + 0.5_dbl * mv2convfac * M * sq_vel
        END DO
    END SUBROUTINE getekin



!---------------------------------------------------------------------------
!> Get temperature routine. Computes the temperature of the system.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE gettemp
        temp = 2.0_dbl * ekin/(3.0_dbl*DBLE(N-1))/kb
        temp_series(MD_step) = temp
    END SUBROUTINE gettemp



!---------------------------------------------------------------------------
!> Force routine. Computes the temperature of the system.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE force
        REAL(kind=dbl) :: r_sq, forcefac, dx, dy, dz
        REAL(kind=dbl) :: r_cutsq, r6, ljconstc6, ljconst12, rinv2
        INTEGER :: i, j

        epot=0.0_dbl
        CALL force_to_zero

        r_cutsq=r_cut*r_cut
        ljconst12 = 4.0_dbl*epsilon*sigma**12
        ljconstc6  = 4.0_dbl*epsilon*sigma**6

        DO i=1, N-1
             DO j=i+1, N
                    
                dx=pbc(rx(i) - rx(j), L)
                dy=pbc(ry(i) - ry(j), L)
                dz=pbc(rz(i) - rz(j), L)
                r_sq = dx*dx + dy*dy + dz*dz

                ! compute force and energy if within cutoff */
                IF (r_sq < r_cutsq) THEN
                   rinv2 = 1.0_dbl/r_sq
                   r6 = rinv2*rinv2*rinv2
                   forcefac = (12.0_dbl*ljconst12*r6 - 6.0_dbl*ljconstc6)*r6*rinv2
                   epot = epot + r6*(ljconst12*r6 - ljconstc6)

                    fx(i) = fx(i) + dx*forcefac
                    fy(i) = fy(i) + dy*forcefac
                    fz(i) = fz(i) + dz*forcefac

                    fx(j) = fx(j) - dx*forcefac
                    fy(j) = fy(j) - dy*forcefac
                    fz(j) = fz(j) - dz*forcefac
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
        DO i=1, N
            vx(i) = vx(i) + 0.5_dbl * dt * fx(i) / (mv2convfac  * M)
            vy(i) = vy(i) + 0.5_dbl * dt * fy(i) / (mv2convfac  * M)
            vz(i) = vz(i) + 0.5_dbl * dt * fz(i) / (mv2convfac  * M)
            rx(i) = rx(i) + dt*vx(i)
            ry(i) = ry(i) + dt*vy(i)
            rz(i) = rz(i) + dt*vz(i)
        END DO

        ! compute forces and potential energy 
        CALL force

        ! second part: propagate velocities by another half step */
        DO i=1, N
            vx(i) = vx(i) + 0.5_dbl * dt * fx(i) / (mv2convfac  * M)
            vy(i) = vy(i) + 0.5_dbl * dt * fy(i) / (mv2convfac  * M)
            vz(i) = vz(i) + 0.5_dbl * dt * fz(i) / (mv2convfac  * M)
        END DO
    END SUBROUTINE velverlet



!---------------------------------------------------------------------------
!> Poors-man Thermostat routine. Rescales temperature of the system.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE thermostat
        USE kinds
        USE MD_system
        USE Phys_consts
        IMPLICIT NONE

        INTEGER :: i

        DO i=1, N
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

        DO i=1, N
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

        cells_num = N/4
        cbrt_cells = INT(cells_num**(1.0/3))
        lattice_spacing = L/cbrt_cells

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

    DO i=1, N
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
            DO j=i+1, N
                dx=pbc(rx(i) - rx(j), L)
                dy=pbc(ry(i) - ry(j), L)
                dz=pbc(rz(i) - rz(j), L)
                dists(k) = SQRT(dx*dx + dy*dy + dz*dz)
                k=k+1
            END DO
        END DO
    END SUBROUTINE get_distances

END MODULE physics

































