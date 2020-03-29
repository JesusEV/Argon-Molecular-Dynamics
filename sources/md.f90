PROGRAM MD
    USE kinds
    USE io
    USE utils
    USE mdsys
    USE physics
    IMPLICIT NONE
  
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


    ! allocate storage for simulation data.
    ALLOCATE(rx(natoms),ry(natoms),rz(natoms),&
             vx(natoms),vy(natoms),vz(natoms), &
             fx(natoms),fy(natoms),fz(natoms),&
             dists(pair_num),temp_series(nsteps))

!-----------------------------------------------
    CAll fcc_lattice_positions_init
    CALL MaxBoltz_Dist_vel_init
    CALL force_to_zero

    CALL ioopen(ergfile, trajfile, distfile, tempsfile)

    WRITE(stdout, *) 'Starting simulation with ', natoms, ' atoms for', nsteps, ' steps'
    WRITE(stdout, *) '    MD_step           TEMP                 EKIN                  EPOT&
       &                ETOT'
    ! main MD loop 
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

    ! clean up: close files, free memory
    CALL output_temps
    WRITE(stdout,'(A)') 'Simulation Done.'
    CALL ioclose

    DEALLOCATE(rx,ry,rz,vx,vy,vz,fx,fy,fz,dists,temp_series)
END PROGRAM MD
