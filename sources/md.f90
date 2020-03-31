!------------------------------------------------------------------------------
!        PROGRAM Molecular Dynamics Liquid Argon`
!------------------------------------------------------------------------------
! Program        : MD
!
! DESCRIPTION:
!> This program simulates Argon Liquid subject to the Lennard-Jones potential 
!> in a box .
!------------------------------------------------------------------------------
PROGRAM MD
    USE kinds
    USE utils
    USE MD_system
    USE physics
    USE io
    IMPLICIT NONE


!------------------------------------------------------------------------------
! Reading and Initialization of simulation parameters.  
!------------------------------------------------------------------------------  
    INTEGER :: nprint, i
    CHARACTER(len=sln) :: trajfile = './results/positions.xyz'
    CHARACTER(len=sln) :: enefile = './results/energies.dat'
    CHARACTER(len=sln) :: distfile = './results/distances.dat'
    CHARACTER(len=sln) :: tempsfile = './results/temps.dat'
    CHARACTER(len=sln) :: sqvelsfile = './results/sq_velocities.dat'

    READ(stdin,*) N
    READ(stdin,*) M
    READ(stdin,*) epsilon
    READ(stdin,*) sigma
    READ(stdin,*) r_cut
    READ(stdin,*) L
    READ(stdin,*) MD_steps
    READ(stdin,*) dt
    READ(stdin,*) nprint
    READ(stdin,*) res_temp
    READ(stdin,*) thermUpdate
    pair_num = (N*(N - 1)-(N-tracking_size)*(N -tracking_size- 1))/2


!------------------------------------------------------------------------------
! Allocation of storage for simulation data.
!------------------------------------------------------------------------------  
    ALLOCATE(rx(N),ry(N),rz(N),&
             vx(N),vy(N),vz(N), &
             fx(N),fy(N),fz(N),&
             dists(pair_num),temp_series(MD_steps),&
             vel_series(MD_steps*N))


!------------------------------------------------------------------------------
! Initialization of Positions, velocities and Forces.
!------------------------------------------------------------------------------  
    CAll fcc_lattice_positions_init
    CALL MaxBoltz_Dist_vel_init
    CALL force_to_zero


!------------------------------------------------------------------------------
! Files Openning.
!------------------------------------------------------------------------------  
    CALL ioopen(enefile, trajfile, distfile, tempsfile, sqvelsfile)



!------------------------------------------------------------------------------
! Molecular Dyanamic Main loop.
!------------------------------------------------------------------------------  
    WRITE(stdout, *) 'MD_step'
    DO MD_step=1, MD_steps

        ! Integrate and compute energies and temperature
        CALL velverlet
        CALL getekin
        CALL gettemp

        IF (MOD(MD_step, thermUpdate) == 0 ) CALL thermostat

        ! write output
        IF (mod(MD_step,nprint) == 0) THEN
             CALL get_distances
             CALL output
        END IF
    END DO

    CALL output_temps_and_vels


!------------------------------------------------------------------------------
! clean up: close files, free memory
!------------------------------------------------------------------------------  
    WRITE(stdout,'(A)') 'Simulation Completed.'
    CALL ioclose
    DEALLOCATE(rx,ry,rz,vx,vy,vz,fx,fy,fz,dists,temp_series, vel_series)

END PROGRAM MD
