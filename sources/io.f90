!------------------------------------------------------------------------------
!        MODULE Input/Output
!------------------------------------------------------------------------------
! MODULE        : io
!
! DESCRIPTION:
!> This module is in charge of the reading and writing files with simulation
!> data..
!------------------------------------------------------------------------------
MODULE io
    USE kinds
    USE mdsys
    
    IMPLICIT NONE
    
    PRIVATE 
    INTEGER, PARAMETER :: stdin=5, stdout=6, log=30, xyz=31, dis=32, tem=33
    
    PUBLIC :: ioopen, ioclose, output, output_temps, stdin, stdout

    CONTAINS

!---------------------------------------------------------------------------
!> ioopen routine. Opens all the needed files.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE ioopen(logname, xyzname, disname,temname)
        CHARACTER(LEN=sln) :: logname, xyzname, disname, temname
        OPEN(UNIT=log, FILE=TRIM(logname), STATUS='UNKNOWN', FORM='FORMATTED')
        OPEN(UNIT=xyz, FILE=TRIM(xyzname), STATUS='UNKNOWN', FORM='FORMATTED')
        OPEN(UNIT=dis, FILE=TRIM(disname), STATUS='UNKNOWN', FORM='FORMATTED')
        OPEN(UNIT=tem, FILE=TRIM(temname), STATUS='UNKNOWN', FORM='FORMATTED')
    END SUBROUTINE ioopen



!---------------------------------------------------------------------------
!> ioclose routine. Closes all the needed files.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE ioclose
        CLOSE(UNIT=log)
        CLOSE(UNIT=xyz)
        CLOSE(UNIT=dis)
        CLOSE(UNIT=tem)
    END SUBROUTINE ioclose



!---------------------------------------------------------------------------
!> output routine. Writes simulation data into files.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE output
        INTEGER :: i

        WRITE(log, '(I8,1X,F20.8,1X,F20.8,1X,F20.8,1X,F20.8)') &
             MD_step, temp, ekin, epot, ekin+epot
        WRITE(stdout, '(I8,1X,F20.8,1X,F20.8,1X,F20.8,1X,F20.8)') &
             MD_step, temp, ekin, epot, ekin+epot
        ! WRITE(xyz, '(I8)') natoms
        WRITE(xyz, *)
        WRITE(xyz, '(A,I8,1X,A,F20.8)') 'MD_step=', MD_step, 'etot=', ekin+epot

        DO i=1, natoms
            WRITE(xyz, '(A2,I0.3, 3(1X,F012.8))') &
                'Ar', i, rx(i), ry(i), rz(i)
        END DO

        DO i=1, pair_num
            WRITE(dis, '(F012.8)') dists(i)
        END DO

    END SUBROUTINE output



!---------------------------------------------------------------------------
!> output temps routine. writes temperature data only.
!> @param[in] casefilename
!---------------------------------------------------------------------------
    SUBROUTINE output_temps
        INTEGER :: i
        character(len=23) :: str

        DO i=1, nsteps
            WRITE(str, '(F012.8)') temp_series(i)
            WRITE(tem,'(a)') adjustl(trim(str))    
        END DO
    END SUBROUTINE output_temps

END MODULE io