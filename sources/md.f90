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

    ! allocate storage for simulation data.
    ALLOCATE(rx(natoms),ry(natoms),rz(natoms),&
    vx(natoms),vy(natoms),vz(natoms), &
    fx(natoms),fy(natoms),fz(natoms))

!-----------------------------------------------
    CAll fcc_lattice_positions_init
    CALL MaxBoltz_Dist_vel_init
    CALL force_to_zero

    nfi=0
    CALL force
    CALL getekin
    
    CALL ioopen(ergfile, trajfile)

    WRITE(stdout, *) 'Starting simulation with ', natoms, ' atoms for', nsteps, ' steps'
    WRITE(stdout, *) '    NFI           TEMP                 EKIN                  EPOT&
       &                ETOT'
    CALL output

    ! main MD loop 
    DO nfi=1, nsteps
        ! write output, if requested
        IF (mod(nfi,nprint) == 0) CALL output

        ! propagate system and recompute energies
        CALL velverlet
        CALL getekin
        IF (MOD(nfi, thermUpdate) == 0 ) CALL thermostat
    END DO

    ! clean up: close files, free memory
    WRITE(stdout,'(A)') 'Simulation Done.'
    CALL ioclose

    DEALLOCATE(rx,ry,rz,vx,vy,vz,fx,fy,fz)
END PROGRAM MD
