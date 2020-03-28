MODULE physconst
  USE kinds
  IMPLICIT NONE
  REAL(kind=dbl), PARAMETER :: kboltz =    0.0019872067_dbl   ! boltzman constant in kcal/mol/K
  REAL(kind=dbl), PARAMETER :: mvsq2e = 2390.05736153349_dbl  ! m*v^2 in kcal/mol
  REAL(kind=dbl), PARAMETER :: onesigmavel = 0.001765695_dbl  ! kb T/M in kcal/mol
  PRIVATE
  PUBLIC :: kboltz, mvsq2e, onesigmavel
END MODULE physconst


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


MODULE physics

    USE kinds
    IMPLICIT NONE

    CONTAINS

        ! compute kinetic energy
        SUBROUTINE getekin
          USE kinds
          USE mdsys, ONLY: natoms, mass, temp, ekin, vx, vy, vz
          USE physconst
          IMPLICIT NONE

          INTEGER :: i

          ekin = 0.0_dbl
          DO i=1, natoms
             ekin = ekin + 0.5_dbl * mvsq2e * mass * (vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i))
          END DO

          temp = 2.0_dbl * ekin/(3.0_dbl*DBLE(natoms-1))/kboltz
        END SUBROUTINE getekin

        ! compute forces 
        SUBROUTINE force
          USE kinds
          USE utils
          USE mdsys
          IMPLICIT NONE

          REAL(kind=dbl) :: r, ffac, dx, dy, dz
          INTEGER :: i, j
          
          epot=0.0_dbl
          CALL azzero(fx,natoms)
          CALL azzero(fy,natoms)
          CALL azzero(fz,natoms)

          DO i=1, natoms
             DO j=1, natoms
                ! particles have no interactions with themselves 
                IF (i==j) CYCLE
                    
                ! get distance between particle i and j 
                !        delta = delta - box*(ANINT(delta/box))
                dx=pbc(rx(i) - rx(j), box)
                dy=pbc(ry(i) - ry(j), box)
                dz=pbc(rz(i) - rz(j), box)
                r = SQRT(dx*dx + dy*dy + dz*dz)
              
                ! compute force and energy if within cutoff */
                IF (r < rcut) THEN
                   ffac = -4.0_dbl*epsilon*(-12.0_dbl*((sigma/r)**12)/r   &
                        +6.0_dbl*(sigma/r)**6/r)
                        
                   epot = epot + 0.5_dbl*4.0_dbl*epsilon*((sigma/r)**12 &
                        -(sigma/r)**6.0)

                   fx(i) = fx(i) + dx/r*ffac
                   fy(i) = fy(i) + dy/r*ffac
                   fz(i) = fz(i) + dz/r*ffac
                END IF
             END DO
          END DO
        END SUBROUTINE force


        ! velocity verlet
        SUBROUTINE velverlet
          USE kinds
          USE mdsys
          USE physconst
          IMPLICIT NONE

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

    SUBROUTINE MaxBoltz_Dist_vel_init

        USE utils
        USE kinds
        USE mdsys
        USE physconst

        IMPLICIT NONE

        INTEGER :: i

        DO i=1, natoms
            vx(i) = box_muller_method(onesigmavel,0.0_dbl)*0.1
            vy(i) = box_muller_method(onesigmavel,0.0_dbl)*0.1
            vz(i) = box_muller_method(onesigmavel,0.0_dbl)*0.1
        END DO

    END SUBROUTINE MaxBoltz_Dist_vel_Init

    SUBROUTINE lattice_positions_init

        USE utils
        USE kinds
        USE mdsys
        USE physconst

        IMPLICIT NONE

        INTEGER :: i,j,k
        INTEGER :: lattice_size, cell_counter = 1
        REAL(kind=dbl) :: cbrt_natoms, lattice_spacing

        cbrt_natoms = natoms**(1.0/3)
        lattice_spacing = box*0.9214/cbrt_natoms
        lattice_size = FLOOR(cbrt_natoms) + 1

        iloop: DO i=1, lattice_size
            DO j=1, lattice_size
                DO k=1, lattice_size
                    rx(cell_counter) = i*lattice_spacing + REAL(RAND())*box/100
                    ry(cell_counter) = j*lattice_spacing + REAL(RAND())*box/100
                    rz(cell_counter) = k*lattice_spacing + REAL(RAND())*box/100
                    cell_counter = cell_counter + 1
                    IF (cell_counter == natoms + 1) EXIT iloop
                END DO
            END DO
        END DO iloop

    END SUBROUTINE lattice_positions_Init


END MODULE physics