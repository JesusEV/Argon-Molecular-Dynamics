MODULE io
    USE kinds
    USE mdsys
    
    IMPLICIT NONE
    
    PRIVATE 
    INTEGER, PARAMETER :: stdin=5, stdout=6, log=30, xyz=31
    
    PUBLIC :: ioopen, ioclose, output, stdin, stdout

    CONTAINS

        SUBROUTINE ioopen(logname, xyzname)
            CHARACTER(LEN=sln) :: logname, xyzname
            OPEN(UNIT=log, FILE=TRIM(logname), STATUS='UNKNOWN', FORM='FORMATTED')
            OPEN(UNIT=xyz, FILE=TRIM(xyzname), STATUS='UNKNOWN', FORM='FORMATTED')
        END SUBROUTINE ioopen

        SUBROUTINE ioclose
            CLOSE(UNIT=log)
            CLOSE(UNIT=xyz)
        END SUBROUTINE ioclose

        ! append data to output.
        SUBROUTINE output
            INTEGER :: i
    
            WRITE(log, '(I8,1X,F20.8,1X,F20.8,1X,F20.8,1X,F20.8)') &
                 nfi, temp, ekin, epot, ekin+epot
            WRITE(stdout, '(I8,1X,F20.8,1X,F20.8,1X,F20.8,1X,F20.8)') &
                 nfi, temp, ekin, epot, ekin+epot
            WRITE(xyz, '(I8)') natoms
            WRITE(xyz, '(A,I8,1X,A,F20.8)') 'nfi=', nfi, 'etot=', ekin+epot
    
            DO i=1, natoms
                WRITE(xyz, '(A,1X,F20.8,1X,F20.8,1X,F20.8)') &
                    'Ar ', rx(i), ry(i), rz(i)
            END DO
        END SUBROUTINE output
END MODULE io