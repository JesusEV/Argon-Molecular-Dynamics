!------------------------------------------------------------------------------
!        PROGRAM Molecular Dynamics of Argon`
!------------------------------------------------------------------------------
! Program        : MD
!
! DESCRIPTION:
!> This program simulates Argon gas subject to the Lennard-Jones potential 
!> in a box .
!------------------------------------------------------------------------------
PROGRAM MD
    USE kinds
    USE io
    USE utils
    USE mdsys
    USE physics
    IMPLICIT NONE


!------------------------------------------------------------------------------
! Reading and Initialization of simulation parameters.  
!------------------------------------------------------------------------------  
    INTEGER :: nprint, i
    CHARACTER(len=sln) :: trajfile = './results/positions.dat'
    CHARACTER(len=sln) :: ergfile = './results/energies.dat'
    CHARACTER(len=sln) :: distfile = './results/distances.dat'
    CHARACTER(len=sln) :: tempsfile = './results/temps.dat'

    READ(stdin,*) natoms
    READ(stdin,*) mass
    READ(stdin,*) epsilon
    READ(stdin,*) sigma
    READ(stdin,*) rcut
    READ(stdin,*) box
    READ(stdin,*) nsteps
    READ(stdin,*) dt
    READ(stdin,*) nprint
    READ(stdin,*) res_temp
    READ(stdin,*) thermUpdate
    pair_num = (natoms*(natoms - 1)-(natoms-tracking_size)*(natoms -tracking_size- 1))/2


!------------------------------------------------------------------------------
! Allocation of storage for simulation data.
!------------------------------------------------------------------------------  
    ALLOCATE(rx(natoms),ry(natoms),rz(natoms),&
             vx(natoms),vy(natoms),vz(natoms), &
             fx(natoms),fy(natoms),fz(natoms),&
             dists(pair_num),temp_series(nsteps))


!------------------------------------------------------------------------------
! Initialization of Positions, velocities and Forces.
!------------------------------------------------------------------------------  
    CAll fcc_lattice_positions_init
    CALL MaxBoltz_Dist_vel_init
    CALL force_to_zero


!------------------------------------------------------------------------------
! Files Openning.
!------------------------------------------------------------------------------  
    CALL ioopen(ergfile, trajfile, distfile, tempsfile)



!------------------------------------------------------------------------------
! Molecular Dyanamic Main loop.
!------------------------------------------------------------------------------  
    WRITE(stdout, *) 'Starting simulation with ', natoms, ' atoms for', nsteps, ' steps'
    WRITE(stdout, *) '    MD_step           TEMP                 EKIN                  EPOT&
       &                ETOT'
    DO MD_step=1, nsteps

        ! propagate system and recompute energies
        CALL velverlet
        CALL getekin
        CALL gettemp

        IF (MOD(MD_step, thermUpdate) == 0 ) CALL thermostat

        ! write output, if requested
        IF (mod(MD_step,nprint) == 0) THEN
             CALL get_distances
             CALL output
        END IF
    END DO

    CALL output_temps


!------------------------------------------------------------------------------
! clean up: close files, free memory
!------------------------------------------------------------------------------  
    WRITE(stdout,'(A)') 'Simulation Done.'
    CALL ioclose
    DEALLOCATE(rx,ry,rz,vx,vy,vz,fx,fy,fz,dists,temp_series)

END PROGRAM MD
