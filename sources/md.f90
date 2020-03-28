PROGRAM MD
  USE kinds
  USE io
  USE utils
  USE mdsys
  USE physics
  IMPLICIT NONE
  
  INTEGER :: nprint, i
  CHARACTER(len=sln) :: restfile, trajfile, ergfile

  READ(stdin,*) natoms
  READ(stdin,*) mass
  READ(stdin,*) epsilon
  READ(stdin,*) sigma
  READ(stdin,*) rcut
  READ(stdin,*) box
  CALL getline(stdin,restfile)
  CALL getline(stdin,trajfile)
  CALL getline(stdin,ergfile)
  READ(stdin,*) nsteps
  READ(stdin,*) dt
  READ(stdin,*) nprint

  ! allocate storage for simulation data.
  ALLOCATE(rx(natoms),ry(natoms),rz(natoms),&
       vx(natoms),vy(natoms),vz(natoms), &
       fx(natoms),fy(natoms),fz(natoms))

!-----------------------------------------------
  ! CALL MaxBoltz_Dist_vel_init
  ! CAll lattice_positions_init
! -----------------------------------------------

!-----------------------------------------------
  !read restart 
  OPEN(UNIT=33, FILE=restfile, FORM='FORMATTED', STATUS='OLD')
  DO i=1,natoms
     READ(33,*) rx(i), ry(i), rz(i)
  END DO
  DO i=1,natoms
     READ(33,*) vx(i), vy(i), vz(i)
  END DO
  CLOSE(33)

  ! DO i=1, natoms
  !     WRITE(stdout, *) vx(i), vy(i), vz(i)
  ! END DO
  ! WRITE(stdout, *) '-----------------------------'
  ! DO i=1, natoms
  !     WRITE(stdout, *) rx(i), ry(i), rz(i)
  ! END DO
  ! WRITE(stdout, *) '-----------------------------'
! -----------------------------------------------


  CALL azzero(fx,natoms)
  CALL azzero(fy,natoms)
  CALL azzero(fz,natoms)

  ! initialize forces and energies
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
     IF (mod(nfi,nprint) == 0) THEN
        CALL output
     END IF

        ! propagate system and recompute energies
        CALL velverlet
        CALL getekin
     END DO

     ! clean up: close files, free memory
     WRITE(stdout,'(A)') 'Simulation Done.'
     CALL ioclose

     DEALLOCATE(rx,ry,rz,vx,vy,vz,fx,fy,fz)
END PROGRAM MD
