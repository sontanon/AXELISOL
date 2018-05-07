PROGRAM ELLSOLVE_F90
	! C BINDINGS MODULE USED FOR PROPER C SIZES.
	USE ISO_C_BINDING
	! SYSTEM MODULES.
	USE IFPORT
	USE IFCORE
	! NO IMPLICIT VARIABLES.
	IMPLICIT NONE 
    ! SOLVER RANGES.
    INTEGER(KIND=C_INT), PARAMETER :: NRINTERIOR_MIN = 32
    INTEGER(KIND=C_INT), PARAMETER :: NRINTERIOR_MAX = 2048
    INTEGER(KIND=C_INT), PARAMETER :: NZINTERIOR_MIN = 32
    INTEGER(KIND=C_INT), PARAMETER :: NZINTERIOR_MAX = 2048
    REAL(KIND=C_DOUBLE), PARAMETER :: DR_MAX = 1.0
    REAL(KIND=C_DOUBLE), PARAMETER :: DR_MIN = 0.000976562
    REAL(KIND=C_DOUBLE), PARAMETER :: DZ_MAX = 1.0
    REAL(KIND=C_DOUBLE), PARAMETER :: DZ_MIN = 0.000976562
	! PARAMETERS.
	INTEGER(KIND=C_INT) :: NRINTERIOR = 256
	INTEGER(KIND=C_INT) :: NZINTERIOR = 64
	REAL(KIND=C_DOUBLE) :: DR = 0.03125
	REAL(KIND=C_DOUBLE) :: DZ = 0.125
	INTEGER(KIND=C_INT) :: NORDER = 2
	CHARACTER(LEN=256) :: DIRNAME = 'output'
	CHARACTER(LEN=256) :: SOLVER = 'flat'
	! GRID VARIABLES.
	INTEGER(KIND=C_INT) :: GHOST, NRTOTAL, NZTOTAL, ARRAY_DIM
	! GRID FUNCTIONS.
	REAL(KIND=C_DOUBLE), DIMENSION(:), ALLOCATABLE :: R, Z, U, F, S, RES
	REAL(KIND=C_DOUBLE), DIMENSION(:), ALLOCATABLE :: A, B, C, D, E
	! AUXILIARY VARIABLES.
	REAL(KIND=C_DOUBLE) :: AUX_R, AUX_Z
	INTEGER(KIND=C_INT) :: I, J, K
	! SOLUTION PARAMETERS.
	INTEGER(KIND=C_INT) :: ROBIN = 1
	REAL(KIND=C_DOUBLE) :: UINF = 1.0
	INTEGER(KIND=C_INT) :: RSYM = 1, ZSYM = 1
	INTEGER(KIND=C_INT) :: PRECOND_USE = 0
	! VARIOUS TIMERS.
	REAL(KIND=C_DOUBLE) :: START_TIME(10), END_TIME(10), TIMES(10)
	! NUMBER OF ARGUMENTS.
	INTEGER(KIND=C_INT) :: NUM_ARG
	! ARGUMENTS.
	CHARACTER(LEN=256) :: ARGS(8), STRING
	! SYSTEM COMMAND STATUS.
	INTEGER(KIND=C_INT) :: ERRNUM
	LOGICAL(KIND=4) :: L_RESULT
	! USER INPUT CHARACTER.
	CHARACTER(KIND=C_CHAR) :: OPT

	! PROCESS PARAMETERS.
	! GET ARGUMENTS FROM COMMAND LINE.
	NUM_ARG = NARGS()

	! SANITY CHECK.
	IF (NUM_ARG .NE. 8) THEN
		PRINT *, 'ELLSOLVEF: WARNING! Usage is $./ELLSOLVEF dirname solver norder NrInterior NzInterior dr dz'
        PRINT *, '           [solver] is the type of solver: flat or general.'
		PRINT *, '           [dirname] is a valid directory string name.'
		PRINT *, '           [norder] = is and integer equal to 2 or 4 corresponding to the finite difference order.'
		PRINT *, '           [NrInterior] and [NzInterior] are itegers equal to the number of interior points in r, z.'
		PRINT *, '           [dr] and [dz] are floating point dobles equal to the spatial step in r, z.'
		PRINT *, 'Press (y/n) to procede with default arguments:'
		OPT = GETCHARQQ()
		IF ((OPT == 'Y').OR.(OPT == 'y')) THEN
			PRINT *, 'ELLSOLVEF: User chose to proceed with default arguments.'
		ELSE
			PRINT *, 'ELLSOLVEF: User chose to abort.'
			CALL EXIT(1)
		END IF
	ELSE
		! GET ARGUMENTS DOING SOME SANITY CHECKS.
		! DIRECTORY NAME.
		CALL GETARG(1, DIRNAME)

        ! SOLVER TYPE.
        CALL GETARG(2, SOLVER)
        IF ((SOLVER .NE. 'flat') .AND. (SOLVER .NE. 'general')) THEN
            PRINT *, 'ELLSOLVEF: ERROR! Unrecongnized solver ', SOLVER, ' Only "flat" or "general" is supported.'
            CALL EXIT(1)
        END IF

		! FINITE DIFFERENCE ORDER.
		CALL GETARG(3, STRING)
		READ (STRING, *) NORDER
        IF ((NORDER .NE. 2) .AND. (NORDER .NE. 4)) THEN
            PRINT *, 'ELLSOLVEF: ERROR! Finite difference ', NORDER, ' is not supported, only 2 or 4.'
            CALL EXIT(1)
        END IF

		! NUMBER OF INTERIOR POINTS.
		CALL GETARG(4, STRING)
		READ (STRING, *) NRINTERIOR
        IF ((NRINTERIOR < NRINTERIOR_MIN) .OR. (NRINTERIOR > NRINTERIOR_MAX) .OR. (NRINTERIOR <= 0)) THEN
            PRINT *, 'ELLSOLVEF: ERROR! NrInterior = ', NRINTERIOR, ' is out of range.'
            CALL EXIT(1)
        END IF

		CALL GETARG(5, STRING)
		READ (STRING, *) NZINTERIOR
        IF ((NZINTERIOR < NZINTERIOR_MIN) .OR. (NZINTERIOR > NZINTERIOR_MAX) .OR. (NZINTERIOR <= 0)) THEN
            PRINT *, 'ELLSOLVEF: ERROR! NzInterior = ', NZINTERIOR, ' is out of range.'
            CALL EXIT(1)
        END IF

		! SPATIAL STEPS.
		CALL GETARG(6, STRING)
		READ (STRING, *) DR
        IF ((DR < DR_MIN) .OR. (DR > DR_MAX)) THEN
            PRINT *, 'ELLSOLVEF: ERROR! dr = ', DR, ' is out of range.'
            CALL EXIT(1)
        END IF

		CALL GETARG(7, STRING)
		READ (STRING, *) DZ
        IF ((DZ < DZ_MIN) .OR. (DZ > DZ_MAX)) THEN
            PRINT *, 'ELLSOLVEF: ERROR! dz = ', DZ, ' is out of range.'
            CALL EXIT(1)
        END IF
	END IF

	! DO I/O ON OUTPUT DIRECTORY.
	! CREATE DIRECTORY.
	L_RESULT = MAKEDIRQQ(DIRNAME)
	! CHECK FOR FAILURE.
	IF (.NOT.L_RESULT) THEN
		ERRNUM = GETLASTERRORQQ()
		! PERMISSION DENIED.
		IF (ERRNUM == ERR$ACCES) THEN
			PRINT *, 'ERROR: Permision denied on MKDIR ', DIRNAME
			CALL EXIT(1)
		ELSE IF (ERRNUM == ERR$EXIST) THEN
			PRINT *, 'ELLSOLVEF: WARNING! Directory ', DIRNAME, ' already existes.'
			PRINT *, 'Press (y/n) to procede and possibly overwrite files:'
			OPT = GETCHARQQ()
			IF ((OPT == 'Y').OR.(OPT == 'y')) THEN
				PRINT *, 'ELLSOLVEF: User chose to proceed.'
			ELSE
				PRINT *, 'ELLSOLVEF: User chose to abort.'
				CALL EXIT(1)
			END IF
		ELSE IF (ERRNUM == ERR$NOENT) THEN
			PRINT *, 'ERROR: Path specified was not found on MKDIR ', DIRNAME
			CALL EXIT(1)
		END IF
	ENDIF

	! CD TO DIRECTORY.
	L_RESULT = CHANGEDIRQQ(DIRNAME)
	! CHECK FOR FAILURE.
	IF (.NOT.L_RESULT) THEN
		PRINT *, 'ELLSOLVEF: WARNIN! Could no CD to directory ', DIRNAME
		PRINT *, 'Press (y/n) to procede and write in current directory:'
		OPT = GETCHARQQ()
		IF ((OPT == 'Y').OR.(OPT == 'y')) THEN
			PRINT *, 'ELLSOLVEF: User chose to proceed.'
		ELSE
			PRINT *, 'ELLSOLVEF: User chose to abort.'
			CALL EXIT(1)
		END IF
	ENDIF

	! GET NUMBER OF GHOST ZONES ACCORDING TO FINITE DIFFERENCE ORDER.
	IF (NORDER == 2) THEN
		GHOST = 2
	ELSEIF (NORDER == 4) THEN
		GHOST = 3
	END IF
	! CALCULATE LONGITUDINAL DIMENSIONS.
	NRTOTAL = GHOST + NRINTERIOR + 1
	NZTOTAL = GHOST + NZINTERIOR + 1
	! TOTAL GRID SIZE.
	ARRAY_DIM = NRTOTAL * NZTOTAL

	! PRINT INFO TO SCREEN.
	PRINT *, 'ELLSOLVEF: System parameters are:'
    PRINT *, '      dirname     = ', DIRNAME
    PRINT *, '      solver      = ', SOLVER
	PRINT *, '		order	    = ', NORDER
	PRINT *, '		NrInterior	= ', NRINTERIOR
	PRINT *, '		NzInterior	= ', NZINTERIOR
	PRINT *, '		dr	        = ', DR
	PRINT *, '		dz	        = ', DZ

	! ALLOCATE GRID FUNCTIONS, INDEX LIKE FORTRAN.
	PRINT *, 'ELLSOLVEF: Allocating memory...'
	ALLOCATE(R(1:ARRAY_DIM), Z(1:ARRAY_DIM), U(1:ARRAY_DIM),&
		 F(1:ARRAY_DIM), S(1:ARRAY_DIM), RES(1:ARRAY_DIM))
	ALLOCATE(A(1:ARRAY_DIM), B(1:ARRAY_DIM), C(1:ARRAY_DIM),&
		 D(1:ARRAY_DIM), E(1:ARRAY_DIM))
	PRINT *, 'ELLSOLVEF: Allocated memory.'

	! FILL GRID FUNCTIONS AND LINEAR SOURCE.
	!$OMP PARALLEL SHARED(R, Z, U, F, S, RES, A, B, C, D, E) PRIVATE(AUX_R, AUX_Z, J, K)
	!$OMP DO SCHEDULE(GUIDED)
	DO I = 1, NRTOTAL
		AUX_R = (DBLE(I - GHOST) - 0.5) * DR
		DO J = 1, NZTOTAL
			AUX_Z = (DBLE(J - GHOST) - 0.5) * DZ
			! LINEAR INDEX.
			K = (I - 1) * NZTOTAL + J
			R(K) = AUX_R
			Z(K) = AUX_Z
            ! SET EVERYTHING ELSE TO ZERO.
			U(K) = 0.0
			F(K) = 0.0
            S(K) = 0.0
			RES(K) = 0.0
            A(K) = 0.0
            B(K) = 0.0
            C(K) = 0.0
            D(K) = 0.0
            E(K) = 0.0
		END DO
	END DO
	!$OMP END DO
	!$OMP END PARALLEL
	PRINT *, 'ELLSOLVEF: Filled grid functions.'

	! CALL C TO PRINT ARRAYS: NOTICE HOW VARIABLES ARE PASSED BY VALUE.
	CALL WRITE_SINGLE_FILE(R, 'r.asc'//C_NULL_CHAR, %VAL(NRTOTAL), %VAL(NZTOTAL))
	CALL WRITE_SINGLE_FILE(Z, 'z.asc'//C_NULL_CHAR, %VAL(NRTOTAL), %VAL(NZTOTAL))
	PRINT *, 'ELLSOLVEF: Printed grid functions.'

	! INTIALIZE MEMORY AND PARAMETERS.
	CALL PARDISO_START(%VAL(NRINTERIOR), %VAL(NZINTERIOR))

    ! CHOOSE BETWEEN FLAT AND GENERAL SOLVER.
    ! FLAT SOLVER.
    IF (SOLVER == 'flat') THEN
        ! FILLL LINEAR SOURCE AND RHS.
        !$OMP PARALLEL SHARED(S, F) PRIVATE(AUX_R, AUX_Z)
            !$OMP DO SCHEDULE(GUIDED)
            DO K = 1, ARRAY_DIM
                AUX_R = R(K)
                AUX_Z = Z(K)
                ! LINEAR SOURCE.
			    S(K) = EXP(-AUX_R**2-AUX_Z**2) * (0.5 + AUX_R**2 * (-3.0 + AUX_R**2 + AUX_Z**2))
                ! RHS.
                F(K) = 0.0
            END DO
            !$OMP END DO
        !$OMP END PARALLEL 

        ! WRITE LINEAR SOURCE AND RHS.
        CALL WRITE_SINGLE_FILE(S, 's.asc'//C_NULL_CHAR, %VAL(NRTOTAL), %VAL(NZTOTAL))
        CALL WRITE_SINGLE_FILE(F, 'f.asc'//C_NULL_CHAR, %VAL(NRTOTAL), %VAL(NZTOTAL))

        ! CALL VANILLA SOLVER.
        ! ONCE AGAIN, NOTICE THAT VARIABLES ARE PASSED BY VALUE.
        PRINT *, 'ELLSOLVEF: Calling normal solver.'
        CALL CPU_TIME(START_TIME(1))
        CALL FLAT_LAPLACIAN(U, RES, S, F, %VAL(UINF), %VAL(ROBIN), %VAL(RSYM), %VAL(ZSYM),&
            %VAL(NRINTERIOR), %VAL(NZINTERIOR), %VAL(GHOST),&
            %VAL(DR), %VAL(DZ), %VAL(NORDER),&
            %VAL(PRECOND_USE))
        CALL CPU_TIME(END_TIME(1))
        TIMES(1) = END_TIME(1) - START_TIME(1)

        ! CALCULATE OTHER TYPES OF SOLVER.
        PRINT *, 'ELLSOLVEF: Solving with CGS.'
        PRECOND_USE = 6
        CALL CPU_TIME(START_TIME(4))
        CALL FLAT_LAPLACIAN(U, RES, S, F, %VAL(UINF), %VAL(ROBIN), %VAL(RSYM), %VAL(ZSYM),&
            %VAL(NRINTERIOR), %VAL(NZINTERIOR), %VAL(GHOST),&
            %VAL(DR), %VAL(DZ), %VAL(NORDER),&
            %VAL(PRECOND_USE))
        CALL CPU_TIME(END_TIME(4))
        TIMES(4) = END_TIME(4) - START_TIME(4)
    ELSE IF (SOLVER == 'general') THEN
        ! FILLL COEFFICIENTS AND RHS.
        !$OMP PARALLEL SHARED(A, B, C, D, E, S, F) PRIVATE(AUX_R, AUX_Z)
            !$OMP DO SCHEDULE(GUIDED)
            DO K = 1, ARRAY_DIM
                AUX_R = R(K)
                AUX_Z = Z(K)
                ! LINEAR SOURCE.
                A(K) = AUX_R
                B(K) = 0.0
                C(K) = AUX_R
                D(K) = 1.0
                E(K) = 0.0
			    S(K) = AUX_R * EXP(-AUX_R**2-AUX_Z**2) * (0.5 + AUX_R**2 * (-3.0 + AUX_R**2 + AUX_Z**2))
                ! RHS.
                F(K) = 0.0
            END DO
            !$OMP END DO
        !$OMP END PARALLEL 

        ! WRITE COEFFICCIENTS AND RHS.
        CALL WRITE_SINGLE_FILE(A, 'a.asc'//C_NULL_CHAR, %VAL(NRTOTAL), %VAL(NZTOTAL))
        CALL WRITE_SINGLE_FILE(B, 'b.asc'//C_NULL_CHAR, %VAL(NRTOTAL), %VAL(NZTOTAL))
        CALL WRITE_SINGLE_FILE(C, 'c.asc'//C_NULL_CHAR, %VAL(NRTOTAL), %VAL(NZTOTAL))
        CALL WRITE_SINGLE_FILE(D, 'd.asc'//C_NULL_CHAR, %VAL(NRTOTAL), %VAL(NZTOTAL))
        CALL WRITE_SINGLE_FILE(E, 'e.asc'//C_NULL_CHAR, %VAL(NRTOTAL), %VAL(NZTOTAL))
        CALL WRITE_SINGLE_FILE(S, 's.asc'//C_NULL_CHAR, %VAL(NRTOTAL), %VAL(NZTOTAL))
        CALL WRITE_SINGLE_FILE(F, 'f.asc'//C_NULL_CHAR, %VAL(NRTOTAL), %VAL(NZTOTAL))

        ! CALL VANILLA SOLVER.
        ! ONCE AGAIN, NOTICE THAT VARIABLES ARE PASSED BY VALUE.
        PRINT *, 'ELLSOLVEF: Calling normal solver.'
        CALL CPU_TIME(START_TIME(1))
        CALL GENERAL_ELLIPTIC(U, RES, A, B, C, D, E, S, F, %VAL(UINF), %VAL(ROBIN), %VAL(RSYM), %VAL(ZSYM),&
            %VAL(NRINTERIOR), %VAL(NZINTERIOR), %VAL(GHOST),&
            %VAL(DR), %VAL(DZ), %VAL(NORDER),&
            %VAL(PRECOND_USE))
        CALL CPU_TIME(END_TIME(1))
        TIMES(1) = END_TIME(1) - START_TIME(1)

        ! CALCULATE OTHER TYPES OF SOLVER.
        PRINT *, 'ELLSOLVEF: Solving with CGS.'
        PRECOND_USE = 6
        CALL CPU_TIME(START_TIME(4))
        CALL GENERAL_ELLIPTIC(U, RES, A, B, C, D, E, S, F, %VAL(UINF), %VAL(ROBIN), %VAL(RSYM), %VAL(ZSYM),&
            %VAL(NRINTERIOR), %VAL(NZINTERIOR), %VAL(GHOST),&
            %VAL(DR), %VAL(DZ), %VAL(NORDER),&
            %VAL(PRECOND_USE))
        CALL CPU_TIME(END_TIME(4))
        TIMES(4) = END_TIME(4) - START_TIME(4)
    END IF

	! PRINT EXECUTION TIMES.
	PRINT *, 'ELLSOLVEF: Normal solver took ', TIMES(1), ' seconds.'
	PRINT *, 'ELLSOLVEF: Solver with CGS took ', TIMES(4), ' seconds.'
	
	! WRITE SOLUTION AND RESIDUAL.
	CALL WRITE_SINGLE_FILE(U, 'u.asc'//C_NULL_CHAR, %VAL(NRTOTAL), %VAL(NZTOTAL))
	CALL WRITE_SINGLE_FILE(RES, 'res.asc'//C_NULL_CHAR, %VAL(NRTOTAL), %VAL(NZTOTAL))

	! CLEAR MEMORY AND PARAMETERS FOR PARDISO.
	CALL PARDISO_STOP()

	! DEALLOCATE MEMORY.
	PRINT *, 'ELLSOLVEF: Deallocating memory...'
	DEALLOCATE(R, Z, U, F, S, RES, A, B, C, D, E)
	PRINT *, 'ELLSOLVEF: Deallocated memory.'
	! ALL DONE.
END PROGRAM ELLSOLVE_F90
