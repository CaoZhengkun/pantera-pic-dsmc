! Copyright (C) 2025 von Karman Institute for Fluid Dynamics (VKI)
!
! This file is part of PANTERA PIC-DSMC, a software for the simulation
! of rarefied gases and plasmas using particles.
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.PANTERA PIC-DSMC

MODULE grid_and_partition 

   USE global
   USE mpi_common
   USE screen
   USE tools
   USE periodic_mesh_utils, ONLY: MATCH_PERIODIC_TRIANGLE, PERIODIC_MESH_OK, &
                                  PERIODIC_MESH_NO_MATCH, PERIODIC_MESH_AMBIGUOUS
   USE, INTRINSIC :: ieee_arithmetic, ONLY: ieee_is_finite

   IMPLICIT NONE

   INTEGER, PARAMETER, PUBLIC :: PERIODIC_MAP_OK = 0
   INTEGER, PARAMETER, PUBLIC :: PERIODIC_MAP_INVALID_GRID = 1
   INTEGER, PARAMETER, PUBLIC :: PERIODIC_MAP_INVALID_BOUNDARY = 2
   INTEGER, PARAMETER, PUBLIC :: PERIODIC_MAP_BAD_TRANSLATION = 3
   INTEGER, PARAMETER, PUBLIC :: PERIODIC_MAP_NO_PARTNER = 4
   INTEGER, PARAMETER, PUBLIC :: PERIODIC_MAP_AMBIGUOUS = 5
   INTEGER, PARAMETER, PUBLIC :: PERIODIC_MAP_BAD_GEOMETRY = 6
   INTEGER, PARAMETER, PUBLIC :: PERIODIC_MAP_GROUP_CONFLICT = 7
   INTEGER, PARAMETER, PUBLIC :: PERIODIC_MAP_BAD_FIELD_TOPOLOGY = 8
 
   CONTAINS

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE CELL_FROM_POSITION -> finds the ID of a grid cell from particle position      !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   SUBROUTINE CELL_FROM_POSITION(XP,YP,  IDCELL)

      ! Note: first cell index is 1.

      IMPLICIT NONE

      REAL(KIND=8), INTENT(IN) :: XP, YP ! Location of particle
      INTEGER, INTENT(OUT)     :: IDCELL ! ID of cell to which the particle belongs

      INTEGER      :: XCELL, YCELL
      REAL(KIND=8) :: DX, DY

      IF (GRID_TYPE == RECTILINEAR_UNIFORM) THEN
         ! Cartesian grid with equally spaced cells
         
         DX = (XMAX - XMIN)/NX
         DY = (YMAX - YMIN)/NY

         !WRITE(*,*) 'XP = ', XP, 'XMIN = ', XMIN, 'DX = ', DX
         XCELL = INT((XP-XMIN)/DX)
         YCELL = INT((YP-YMIN)/DY)

         IF (XCELL .GT. (NX-1)) THEN
            XCELL = NX-1
         ELSE IF (XCELL .LT. 0) THEN
            XCELL = 0
         END IF

         IF (YCELL .GT. (NY-1)) THEN
            YCELL = NY-1
         ELSE IF (YCELL .LT. 0) THEN
            YCELL = 0
         END IF

         IDCELL = XCELL + NX*YCELL + 1
      ELSE IF (GRID_TYPE == RECTILINEAR_NONUNIFORM) THEN
         XCELL = BINARY_SEARCH(XP, XCOORD)
         YCELL = BINARY_SEARCH(YP, YCOORD)

         IF (XCELL .GT. NX) THEN
            XCELL = NX
         ELSE IF (XCELL .LT. 1) THEN
            XCELL = 1
         END IF

         IF (YCELL .GT. NY) THEN
            YCELL = NY
         ELSE IF (YCELL .LT. 1) THEN
            YCELL = 1
         END IF

         IDCELL = XCELL + NX*(YCELL-1)
      END IF

   END SUBROUTINE CELL_FROM_POSITION


   SUBROUTINE PROC_FROM_CELL(IDCELL, IDPROC)

      ! Note: processes go from 0 (usually termed the Master) to MPI_N_THREADS - 1
      ! If the number of MPI processes is only 1, then IDPROC is 0, the one and only process.
      ! Otherwise, check the partition style (variable "DOMPART_TYPE")

      IMPLICIT NONE

      INTEGER, INTENT(IN)     :: IDCELL
      INTEGER, INTENT(OUT)     :: IDPROC
      
      IF (N_MPI_THREADS == 1) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Serial operation

         IDPROC = 0
   
      ELSE

         IDPROC = CELL_PROCS(IDCELL)

      END IF

   END SUBROUTINE PROC_FROM_CELL


   SUBROUTINE ASSIGN_CELLS_TO_PROCS


      !INTEGER :: IDCELL
      !INTEGER :: NCELLSPP
      INTEGER :: I, IPROC, IP, IC
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CENTROID
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: WEIGHT
      INTEGER, DIMENSION(:), ALLOCATABLE :: NP_CELLS
      REAL(KIND=8) :: COORDMAX, COORDMIN, WEIGHT_PER_PROC, CUMULATIVE_WEIGHT
      INTEGER, DIMENSION(:), ALLOCATABLE :: ORDER
      REAL(KIND=8), DIMENSION(3) :: CELLCENTROID

      ! Weight is a generic measure that is ideally directly proportional to the computational cost for each cell.

      IF (.NOT. ALLOCATED(CELL_PROCS)) ALLOCATE(CELL_PROCS(NCELLS))

      IF (N_MPI_THREADS == 1) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Serial operation

         CELL_PROCS = 0
      
      ELSE IF (DIMS == 0 .AND. (N_MPI_THREADS .GT. 1)) THEN

         CALL ERROR_ABORT('ERROR! Simulation in 0d can be run only in series.')

      ELSE IF (GRID_TYPE == UNSTRUCTURED) THEN
         IF (.NOT. LOAD_BALANCE) THEN

            ALLOCATE(CENTROID(NCELLS))


            IF (DIMS == 1) THEN

               DO I = 1, NCELLS
                  CELLCENTROID = (U1D_GRID%NODE_COORDS(:, U1D_GRID%CELL_NODES(1,I)) &
                               +  U1D_GRID%NODE_COORDS(:, U1D_GRID%CELL_NODES(2,I))) / 2.
                  IF (PARTITION_STYLE == STRIPSX) THEN
                     CENTROID(I) = CELLCENTROID(1)
                  ELSE
                     CALL ERROR_ABORT('The specified partition style is not supported. Aborting!')
                  END IF
               END DO
            
            ELSE IF (DIMS == 2) THEN

               DO I = 1, NCELLS
                  CELLCENTROID = (U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(1,I)) &
                               +  U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(2,I)) &
                               +  U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(3,I))) / 3.
                  IF (PARTITION_STYLE == STRIPSX) THEN
                     CENTROID(I) = CELLCENTROID(1)
                  ELSE IF (PARTITION_STYLE == STRIPSY) THEN
                     CENTROID(I) = CELLCENTROID(2)
                  ELSE IF (PARTITION_STYLE == SLICESZ) THEN
                     CENTROID(I) = ATAN2(CELLCENTROID(2), CELLCENTROID(1))
                  ELSE
                     CALL ERROR_ABORT('The specified partition style is not supported. Aborting!')
                  END IF
               END DO

            ELSE IF (DIMS == 3) THEN

               DO I = 1, NCELLS
                  CELLCENTROID = (U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(1,I)) &
                               +  U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(2,I)) &
                               +  U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(3,I)) &
                               +  U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(4,I))) / 4.
                  IF (PARTITION_STYLE == STRIPSX) THEN
                     CENTROID(I) = CELLCENTROID(1)
                  ELSE IF (PARTITION_STYLE == STRIPSY) THEN
                     CENTROID(I) = CELLCENTROID(2)
                  ELSE IF (PARTITION_STYLE == STRIPSZ) THEN
                     CENTROID(I) = CELLCENTROID(3)
                  ELSE IF (PARTITION_STYLE == SLICESX) THEN
                     CENTROID(I) = ATAN2(CELLCENTROID(2), CELLCENTROID(3))
                  ELSE
                     CALL ERROR_ABORT('The specified partition style is not supported. Aborting!')
                  END IF
               END DO
            
            END IF

            COORDMAX = MAXVAL(CENTROID)
            COORDMIN = MINVAL(CENTROID)

            !WRITE(*,*) 'COORDMAX = ', COORDMAX, 'COORDMIN = ', COORDMIN

            DO I = 1, NCELLS
               CELL_PROCS(I) = INT((CENTROID(I)-COORDMIN)/(COORDMAX-COORDMIN)*REAL(N_MPI_THREADS))
               IF (CELL_PROCS(I) < 0) CELL_PROCS(I) = 0
               IF (CELL_PROCS(I) >= N_MPI_THREADS) CELL_PROCS(I) = N_MPI_THREADS - 1
            END DO

            DEALLOCATE(CENTROID)

         ELSE

            ALLOCATE(CENTROID(NCELLS))
            ALLOCATE(WEIGHT(NCELLS))
            ALLOCATE(ORDER(NCELLS))
            ALLOCATE(NP_CELLS(NCELLS))
            ! Assign weight
            NP_CELLS = 0
            DO IP = 1, NP_PROC
               NP_CELLS(particles(IP)%IC) = NP_CELLS(particles(IP)%IC) + 1
            END DO

            IF (PROC_ID .EQ. 0) THEN
               CALL MPI_REDUCE(MPI_IN_PLACE, NP_CELLS, NCELLS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            ELSE
               CALL MPI_REDUCE(NP_CELLS,     NP_CELLS, NCELLS, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            END IF
      
            CALL MPI_BCAST(NP_CELLS, NCELLS, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

            WEIGHT = DBLE(NP_CELLS)
            
            WEIGHT = WEIGHT + 0.1*SUM(WEIGHT) / DBLE(NCELLS)
            
            WEIGHT_PER_PROC = SUM(WEIGHT) / DBLE(N_MPI_THREADS)


            IF (DIMS == 1) THEN

               DO I = 1, NCELLS
                  CELLCENTROID = (U1D_GRID%NODE_COORDS(:, U1D_GRID%CELL_NODES(1,I)) &
                               +  U1D_GRID%NODE_COORDS(:, U1D_GRID%CELL_NODES(2,I))) / 2.
                  IF (PARTITION_STYLE == STRIPSX) THEN
                     CENTROID(I) = CELLCENTROID(1)
                  ELSE
                     CALL ERROR_ABORT('The specified partition style is not supported. Aborting!')
                  END IF

                  ORDER(I) = I
               END DO
             
            ELSE IF (DIMS == 2) THEN

               DO I = 1, NCELLS
                  CELLCENTROID = (U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(1,I)) &
                               +  U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(2,I)) &
                               +  U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(3,I))) / 3.
                  IF (PARTITION_STYLE == STRIPSX) THEN
                     CENTROID(I) = CELLCENTROID(1)
                  ELSE IF (PARTITION_STYLE == STRIPSY) THEN
                     CENTROID(I) = CELLCENTROID(2)
                  ELSE IF (PARTITION_STYLE == STRIPSZ) THEN
                     CALL ERROR_ABORT('The specified partition style is not supported. Aborting!')
                  END IF

                  ORDER(I) = I
               END DO
         
            ELSE IF (DIMS == 3) THEN

               DO I = 1, NCELLS
                  CELLCENTROID = (U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(1,I)) &
                               +  U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(2,I)) &
                               +  U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(3,I)) &
                               +  U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(4,I))) / 4.
                  IF (PARTITION_STYLE == STRIPSX) THEN
                     CENTROID(I) = CELLCENTROID(1)
                  ELSE IF (PARTITION_STYLE == STRIPSY) THEN
                     CENTROID(I) = CELLCENTROID(2)
                  ELSE IF (PARTITION_STYLE == STRIPSZ) THEN
                     CENTROID(I) = CELLCENTROID(3)
                  END IF

                  ORDER(I) = I
               END DO
            END IF

            CALL QUICKSORT(CENTROID, ORDER, 1, NCELLS)

            IPROC = 0
            CUMULATIVE_WEIGHT = 0
            DO I = 1, NCELLS
               IC = ORDER(I)
               CUMULATIVE_WEIGHT = CUMULATIVE_WEIGHT + WEIGHT(IC)
               CELL_PROCS(IC) = IPROC

               IF (CUMULATIVE_WEIGHT > WEIGHT_PER_PROC) THEN
                  IPROC = IPROC + 1
                  CUMULATIVE_WEIGHT = 0
               END IF

            END DO

            DEALLOCATE(CENTROID)
            DEALLOCATE(WEIGHT)
            DEALLOCATE(ORDER)
            DEALLOCATE(NP_CELLS)

         END IF
      ELSE
         CALL ERROR_ABORT('The specified partition style is not supported. Aborting!')
      END IF

      ! ELSE IF (DOMPART_TYPE == 0) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ S T R I P S
      ! !
      ! ! Domain partition type: strips
      ! !
      ! ! This domain partition assigns the same number of grid cells to every process.
      ! ! The domain partition is thus in strips along the "x" axis, as the cells are assumed 
      ! ! to be numbered along x.

      !    NCELLSPP = CEILING(REAL(NCELLS)/REAL(N_MPI_THREADS))

      !    DO IDCELL = 1, NCELLS
      !       ! 3) Here is the process ID 
      !       CELL_PROCS(IDCELL) = INT((IDCELL-1)/NCELLSPP) ! Before was giving wrong result for peculiar combinations of NCELLS and NCELLSPP
      !    END DO

      ! ELSE IF (DOMPART_TYPE == 1) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ B L O C K S
      ! !
      ! ! Domain partition type: blocks
      ! !
      ! ! This domain partition style divides the domain into blocks
      ! ! Note that blocks in this definition are independent from cells! IT IS POSSIBLE TO
      ! ! GENERALIZE IT.

      !    DO IDCELL = 1, NCELLS
      !       I = MOD(IDCELL - 1, NX)
      !       J = (IDCELL - 1)/NX
      !       CELL_PROCS(IDCELL) = I*N_BLOCKS_X/NX + N_BLOCKS_X*(J*N_BLOCKS_Y/NY)
      !    END DO

      ! ELSE ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ E R R O R

      !    WRITE(*,*) 'ATTENTION! Value for variable DOMPART_TYPE = ', DOMPART_TYPE
      !    WRITE(*,*) ' not recognized! Check input file!! In: PROC_FROM_POSITION() ABORTING!'
      !    STOP
      ! END IF

   END SUBROUTINE ASSIGN_CELLS_TO_PROCS


   SUBROUTINE VALIDATE_PERIODIC_BOUNDARY_GROUPS

      IMPLICIT NONE

      INTEGER :: IPG, NUM_GROUP_FACES, ERROR_CODE
      LOGICAL :: HAS_PERIODIC_MASTER
      CHARACTER(LEN=512) :: ERROR_MESSAGE
      REAL(KIND=8) :: MATCH_TOLERANCE

      IF (.NOT. ALLOCATED(GRID_BC)) RETURN

      HAS_PERIODIC_MASTER = .FALSE.
      DO IPG = 1, N_GRID_BC
         IF (.NOT. ANY(GRID_BC(IPG)%PARTICLE_BC == PERIODIC_MASTER)) CYCLE
         HAS_PERIODIC_MASTER = .TRUE.

         NUM_GROUP_FACES = 0
         IF (GRID_TYPE /= UNSTRUCTURED) THEN
            CALL ERROR_ABORT('Periodic master boundary groups require an unstructured mesh.')
            RETURN
         END IF

         SELECT CASE (DIMS)
         CASE (1)
            IF (.NOT. ALLOCATED(U1D_GRID%CELL_EDGES_PG)) THEN
               CALL ERROR_ABORT('Periodic master group has no 1D boundary-element map.')
               RETURN
            END IF
            NUM_GROUP_FACES = COUNT(U1D_GRID%CELL_EDGES_PG == IPG)
         CASE (2)
            IF (.NOT. ALLOCATED(U2D_GRID%CELL_EDGES_PG)) THEN
               CALL ERROR_ABORT('Periodic master group has no 2D boundary-element map.')
               RETURN
            END IF
            NUM_GROUP_FACES = COUNT(U2D_GRID%CELL_EDGES_PG == IPG)
         CASE (3)
            IF (.NOT. ALLOCATED(U3D_GRID%CELL_FACES_PG)) THEN
               CALL ERROR_ABORT('Periodic master group has no 3D boundary-face map.')
               RETURN
            END IF
            NUM_GROUP_FACES = COUNT(U3D_GRID%CELL_FACES_PG == IPG)
         CASE DEFAULT
            CALL ERROR_ABORT('Periodic master boundary groups require Dimensions 1, 2 or 3.')
            RETURN
         END SELECT

         IF (NUM_GROUP_FACES == 0) THEN
            WRITE(*,'(A,A)') 'Periodic master group has no mapped boundary elements: ', &
                             TRIM(GRID_BC(IPG)%PHYSICAL_GROUP_NAME)
            CALL ERROR_ABORT('Check the physical group name and mesh boundary elements.')
            RETURN
         END IF
      END DO

      IF (HAS_PERIODIC_MASTER .AND. DIMS == 3) THEN
         CALL BUILD_PERIODIC_BOUNDARY_FACE_MAP(ERROR_CODE, ERROR_MESSAGE, MATCH_TOLERANCE)
         IF (ERROR_CODE /= PERIODIC_MAP_OK) THEN
            IF (PROC_ID == 0) THEN
               WRITE(*,'(A)') TRIM(ERROR_MESSAGE)
               IF (MATCH_TOLERANCE > 0.d0) &
                  WRITE(*,'(A,ES14.6)') 'Periodic face-matching tolerance: ', MATCH_TOLERANCE
               IF (ERROR_CODE == PERIODIC_MAP_NO_PARTNER) &
                  WRITE(*,'(A)') 'Periodic pairing stopped at the first unmatched face (count at least 1).'
               IF (ERROR_CODE == PERIODIC_MAP_AMBIGUOUS) &
                  WRITE(*,'(A)') 'Periodic pairing stopped at the first ambiguous/reused face (count at least 1).'
            END IF
            CALL ERROR_ABORT('3D periodic boundary-face pairing failed during initialization.')
         END IF
         IF (PROC_ID == 0) THEN
            WRITE(*,'(A,ES14.6)') 'Periodic face-matching tolerance: ', MATCH_TOLERANCE
            WRITE(*,'(A,I0,A)') '3D periodic pairing summary: pairs=', &
               U3D_GRID%NUM_PERIODIC_FACE_PAIRS, ', unmatched=0, ambiguous=0.'
         END IF
      END IF

      ! Some geometry-only unit tests intentionally construct U3D_GRID without
      ! populating the global NNODES counter; field mapping is deferred until
      ! the production initialization path has both values available.
      IF (DIMS == 3 .AND. GRID_TYPE == UNSTRUCTURED .AND. NNODES > 0) THEN
         CALL BUILD_PERIODIC_FIELD_DOF_MAP(ERROR_CODE, ERROR_MESSAGE)
         IF (ERROR_CODE /= PERIODIC_MAP_OK) THEN
            IF (PROC_ID == 0) WRITE(*,'(A)') TRIM(ERROR_MESSAGE)
            CALL ERROR_ABORT('3D periodic field-degree-of-freedom mapping failed during initialization.')
         END IF
      END IF

   END SUBROUTINE VALIDATE_PERIODIC_BOUNDARY_GROUPS


   SUBROUTINE BUILD_PERIODIC_FIELD_DOF_MAP(ERROR_CODE, ERROR_MESSAGE)

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: ERROR_CODE
      CHARACTER(LEN=*), INTENT(OUT) :: ERROR_MESSAGE

      INTEGER :: I, J, IC, IFACE, PARTNER_CELL, PARTNER_FACE
      INTEGER :: SOURCE_NODE, TARGET_NODE, TARGET_VERTEX, ROOT_I, ROOT_J
      INTEGER :: NUMBER_OF_ROOTS, CURRENT_DOF
      INTEGER, DIMENSION(:), ALLOCATABLE :: PARENT, ROOT_TO_DOF

      ERROR_CODE = PERIODIC_MAP_OK
      ERROR_MESSAGE = ''

      IF (NNODES < 1) THEN
         ERROR_CODE = PERIODIC_MAP_INVALID_GRID
         ERROR_MESSAGE = 'Periodic field mapping requires at least one mesh node.'
         RETURN
      END IF
      IF (DIMS /= 3 .OR. GRID_TYPE /= UNSTRUCTURED) THEN
         IF (ALLOCATED(NODE_TO_FIELD_DOF)) DEALLOCATE(NODE_TO_FIELD_DOF)
         IF (ALLOCATED(FIELD_DOF_REPRESENTATIVE)) DEALLOCATE(FIELD_DOF_REPRESENTATIVE)
         FIELD_DOF_COUNT = NNODES
         ALLOCATE(NODE_TO_FIELD_DOF(0:NNODES-1))
         ALLOCATE(FIELD_DOF_REPRESENTATIVE(0:NNODES-1))
         DO I = 0, NNODES-1
            NODE_TO_FIELD_DOF(I) = I
            FIELD_DOF_REPRESENTATIVE(I) = I
         END DO
         RETURN
      END IF
      IF (.NOT. ALLOCATED(U3D_GRID%NODE_COORDS) .OR. &
          .NOT. ALLOCATED(U3D_GRID%FACE_NODES) .OR. &
          U3D_GRID%NUM_NODES /= NNODES) THEN
         ERROR_CODE = PERIODIC_MAP_INVALID_GRID
         ERROR_MESSAGE = 'Periodic field mapping found inconsistent 3D node storage.'
         RETURN
      END IF

      IF (ALLOCATED(NODE_TO_FIELD_DOF)) DEALLOCATE(NODE_TO_FIELD_DOF)
      IF (ALLOCATED(FIELD_DOF_REPRESENTATIVE)) DEALLOCATE(FIELD_DOF_REPRESENTATIVE)
      ALLOCATE(NODE_TO_FIELD_DOF(0:NNODES-1))

      ALLOCATE(PARENT(0:NNODES-1))
      DO I = 0, NNODES-1
         PARENT(I) = I
      END DO

      IF (ALLOCATED(U3D_GRID%PERIODIC_PARTNER_CELL) .AND. &
          ALLOCATED(U3D_GRID%PERIODIC_PARTNER_FACE) .AND. &
          ALLOCATED(U3D_GRID%PERIODIC_VERTEX_PERM)) THEN
         DO IC = 1, U3D_GRID%NUM_CELLS
            DO IFACE = 1, 4
               PARTNER_CELL = U3D_GRID%PERIODIC_PARTNER_CELL(IFACE,IC)
               IF (PARTNER_CELL < 1) CYCLE
               PARTNER_FACE = U3D_GRID%PERIODIC_PARTNER_FACE(IFACE,IC)
               IF (PARTNER_CELL > U3D_GRID%NUM_CELLS .OR. PARTNER_FACE < 1 .OR. PARTNER_FACE > 4) THEN
                  ERROR_CODE = PERIODIC_MAP_BAD_FIELD_TOPOLOGY
                  ERROR_MESSAGE = 'Periodic field mapping found an invalid partner cell or face.'
                  DEALLOCATE(PARENT)
                  RETURN
               END IF
               DO J = 1, 3
                  TARGET_VERTEX = U3D_GRID%PERIODIC_VERTEX_PERM(J,IFACE,IC)
                  IF (TARGET_VERTEX < 1 .OR. TARGET_VERTEX > 3) THEN
                     ERROR_CODE = PERIODIC_MAP_BAD_FIELD_TOPOLOGY
                     ERROR_MESSAGE = 'Periodic field mapping found an invalid face-vertex permutation.'
                     DEALLOCATE(PARENT)
                     RETURN
                  END IF
                  SOURCE_NODE = U3D_GRID%FACE_NODES(J,IFACE,IC)
                  TARGET_NODE = U3D_GRID%FACE_NODES(TARGET_VERTEX,PARTNER_FACE,PARTNER_CELL)
                  IF (SOURCE_NODE < 1 .OR. SOURCE_NODE > NNODES .OR. &
                      TARGET_NODE < 1 .OR. TARGET_NODE > NNODES) THEN
                     ERROR_CODE = PERIODIC_MAP_BAD_FIELD_TOPOLOGY
                     ERROR_MESSAGE = 'Periodic field mapping found a partner node outside the mesh.'
                     DEALLOCATE(PARENT)
                     RETURN
                  END IF

                  ROOT_I = FIND_PERIODIC_NODE_ROOT(PARENT,SOURCE_NODE-1)
                  ROOT_J = FIND_PERIODIC_NODE_ROOT(PARENT,TARGET_NODE-1)
                  IF (ROOT_I /= ROOT_J) PARENT(MAX(ROOT_I,ROOT_J)) = MIN(ROOT_I,ROOT_J)
               END DO
            END DO
         END DO
      END IF

      DO I = 0, NNODES-1
         PARENT(I) = FIND_PERIODIC_NODE_ROOT(PARENT,I)
      END DO

      NUMBER_OF_ROOTS = COUNT([(PARENT(I) == I, I=0,NNODES-1)])
      FIELD_DOF_COUNT = NUMBER_OF_ROOTS
      ALLOCATE(FIELD_DOF_REPRESENTATIVE(0:FIELD_DOF_COUNT-1))
      ALLOCATE(ROOT_TO_DOF(0:NNODES-1))
      ROOT_TO_DOF = -1

      CURRENT_DOF = 0
      DO I = 0, NNODES-1
         IF (PARENT(I) /= I) CYCLE
         FIELD_DOF_REPRESENTATIVE(CURRENT_DOF) = I
         ROOT_TO_DOF(I) = CURRENT_DOF
         CURRENT_DOF = CURRENT_DOF + 1
      END DO

      DO I = 0, NNODES-1
         NODE_TO_FIELD_DOF(I) = ROOT_TO_DOF(PARENT(I))
      END DO

      IF (ANY(NODE_TO_FIELD_DOF < 0) .OR. ANY(NODE_TO_FIELD_DOF >= FIELD_DOF_COUNT)) THEN
         ERROR_CODE = PERIODIC_MAP_BAD_FIELD_TOPOLOGY
         ERROR_MESSAGE = 'Periodic field mapping produced an invalid compressed degree of freedom.'
      END IF

      DEALLOCATE(ROOT_TO_DOF)
      DEALLOCATE(PARENT)

   CONTAINS

      RECURSIVE INTEGER FUNCTION FIND_PERIODIC_NODE_ROOT(PARENT_ARRAY, NODE) RESULT(ROOT)
         INTEGER, DIMENSION(0:), INTENT(INOUT) :: PARENT_ARRAY
         INTEGER, INTENT(IN) :: NODE
         IF (PARENT_ARRAY(NODE) == NODE) THEN
            ROOT = NODE
         ELSE
            PARENT_ARRAY(NODE) = FIND_PERIODIC_NODE_ROOT(PARENT_ARRAY,PARENT_ARRAY(NODE))
            ROOT = PARENT_ARRAY(NODE)
         END IF
      END FUNCTION FIND_PERIODIC_NODE_ROOT

   END SUBROUTINE BUILD_PERIODIC_FIELD_DOF_MAP


   SUBROUTINE BUILD_PERIODIC_BOUNDARY_FACE_MAP(ERROR_CODE, ERROR_MESSAGE, MATCH_TOLERANCE_USED)

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: ERROR_CODE
      CHARACTER(LEN=*), INTENT(OUT) :: ERROR_MESSAGE
      REAL(KIND=8), INTENT(OUT), OPTIONAL :: MATCH_TOLERANCE_USED

      INTEGER :: IC, IFACE, IPG, PG, I, AXIS, N_BOUNDARY_FACES, N_CANDIDATES
      INTEGER :: MASTER_FACE_COUNT, SLAVE_FACE_COUNT, SLAVE_GROUP, MATCH_COUNT
      INTEGER :: MASTER_CELL, MASTER_FACE, CANDIDATE_INDEX, CANDIDATE_CELL, CANDIDATE_FACE
      INTEGER :: CANDIDATE_GROUP, MATCHED_CELL, MATCHED_FACE, MATCHED_GROUP
      INTEGER :: LOWER, UPPER, MIDPOINT, REVERSE_PERMUTATION(3), GEOMETRY_STATUS
      REAL(KIND=8) :: L_REF, COORDINATE_SCALE, ABS_TOL, REL_TOL, TOLERANCE
      REAL(KIND=8) :: EXPECTED_CENTROID_AXIS, AREA_TOLERANCE
      REAL(KIND=8), DIMENSION(3) :: EXTENT, MASTER_CENTROID, CANDIDATE_CENTROID
      REAL(KIND=8), DIMENSION(3) :: MASTER_NORMAL, CANDIDATE_NORMAL, TRANSLATION
      REAL(KIND=8), DIMENSION(3,3) :: MASTER_COORDS, CANDIDATE_COORDS
      REAL(KIND=8), DIMENSION(3) :: NORMALIZED_MASTER, NORMALIZED_CANDIDATE
      INTEGER, DIMENSION(3) :: VERTEX_PERMUTATION, MATCHED_PERMUTATION
      INTEGER, DIMENSION(:), ALLOCATABLE :: CANDIDATE_CELLS, CANDIDATE_FACES
      INTEGER, DIMENSION(:), ALLOCATABLE :: CANDIDATE_GROUPS, CANDIDATE_ORDER
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CANDIDATE_AXIS_COORDS
      LOGICAL, DIMENSION(:), ALLOCATABLE :: IS_MASTER_GROUP, IS_SLAVE_GROUP
      LOGICAL, DIMENSION(:,:), ALLOCATABLE :: USED_PARTNER_FACE
      LOGICAL :: FACE_METRICS_VALID

      ERROR_CODE = PERIODIC_MAP_OK
      ERROR_MESSAGE = ''
      IF (PRESENT(MATCH_TOLERANCE_USED)) MATCH_TOLERANCE_USED = 0.d0
      U3D_GRID%NUM_PERIODIC_FACE_PAIRS = 0

      IF (DIMS /= 3 .OR. GRID_TYPE /= UNSTRUCTURED .OR. .NOT. ALLOCATED(GRID_BC) .OR. &
          .NOT. ALLOCATED(U3D_GRID%NODE_COORDS) .OR. .NOT. ALLOCATED(U3D_GRID%FACE_NODES) .OR. &
          .NOT. ALLOCATED(U3D_GRID%FACE_NORMAL) .OR. .NOT. ALLOCATED(U3D_GRID%FACE_AREA) .OR. &
          .NOT. ALLOCATED(U3D_GRID%CELL_FACES_PG)) THEN
         ERROR_CODE = PERIODIC_MAP_INVALID_GRID
         ERROR_MESSAGE = '3D periodic pairing requires a loaded unstructured tetrahedral boundary map.'
         RETURN
      END IF

      IF (BOOL_X_PERIODIC .OR. BOOL_Y_PERIODIC .OR. BOOL_Z_PERIODIC .OR. ANY(BOOL_PERIODIC)) THEN
         ERROR_CODE = PERIODIC_MAP_GROUP_CONFLICT
         ERROR_MESSAGE = 'Domain_periodicity cannot be combined with 3D unstructured periodic face groups.'
         RETURN
      END IF

      IF (U3D_GRID%NUM_NODES < 1 .OR. U3D_GRID%NUM_CELLS < 1 .OR. N_GRID_BC < 1) THEN
         ERROR_CODE = PERIODIC_MAP_INVALID_GRID
         ERROR_MESSAGE = '3D periodic pairing found an empty mesh or boundary-group table.'
         RETURN
      END IF

      EXTENT = MAXVAL(U3D_GRID%NODE_COORDS, DIM=2) - MINVAL(U3D_GRID%NODE_COORDS, DIM=2)
      L_REF = NORM2(EXTENT)
      COORDINATE_SCALE = MAX(1.d0, MAXVAL(ABS(U3D_GRID%NODE_COORDS)))
      ABS_TOL = 64.d0 * EPSILON(1.d0) * COORDINATE_SCALE
      REL_TOL = 1.d-10
      TOLERANCE = ABS_TOL + REL_TOL * L_REF
      IF (.NOT. IEEE_IS_FINITE(L_REF) .OR. L_REF <= 0.d0 .OR. &
          .NOT. IEEE_IS_FINITE(TOLERANCE) .OR. TOLERANCE <= 0.d0) THEN
         ERROR_CODE = PERIODIC_MAP_INVALID_GRID
         ERROR_MESSAGE = '3D periodic pairing could not derive a finite mesh-scale tolerance.'
         RETURN
      END IF
      IF (PRESENT(MATCH_TOLERANCE_USED)) MATCH_TOLERANCE_USED = TOLERANCE

      ALLOCATE(IS_MASTER_GROUP(N_GRID_BC), IS_SLAVE_GROUP(N_GRID_BC))
      IS_MASTER_GROUP = .FALSE.
      IS_SLAVE_GROUP = .FALSE.
      DO IPG = 1, N_GRID_BC
         IS_MASTER_GROUP(IPG) = ANY(GRID_BC(IPG)%PARTICLE_BC == PERIODIC_MASTER)
      END DO
      IF (.NOT. ANY(IS_MASTER_GROUP)) RETURN

      IF (ALLOCATED(U3D_GRID%PERIODIC_PARTNER_CELL)) DEALLOCATE(U3D_GRID%PERIODIC_PARTNER_CELL)
      IF (ALLOCATED(U3D_GRID%PERIODIC_PARTNER_FACE)) DEALLOCATE(U3D_GRID%PERIODIC_PARTNER_FACE)
      IF (ALLOCATED(U3D_GRID%PERIODIC_PARTNER_GROUP)) DEALLOCATE(U3D_GRID%PERIODIC_PARTNER_GROUP)
      IF (ALLOCATED(U3D_GRID%PERIODIC_VERTEX_PERM)) DEALLOCATE(U3D_GRID%PERIODIC_VERTEX_PERM)
      IF (ALLOCATED(U3D_GRID%PERIODIC_TRANSLATION)) DEALLOCATE(U3D_GRID%PERIODIC_TRANSLATION)
      ALLOCATE(U3D_GRID%PERIODIC_PARTNER_CELL(4,U3D_GRID%NUM_CELLS))
      ALLOCATE(U3D_GRID%PERIODIC_PARTNER_FACE(4,U3D_GRID%NUM_CELLS))
      ALLOCATE(U3D_GRID%PERIODIC_PARTNER_GROUP(4,U3D_GRID%NUM_CELLS))
      ALLOCATE(U3D_GRID%PERIODIC_VERTEX_PERM(3,4,U3D_GRID%NUM_CELLS))
      ALLOCATE(U3D_GRID%PERIODIC_TRANSLATION(3,4,U3D_GRID%NUM_CELLS))
      U3D_GRID%PERIODIC_PARTNER_CELL = -1
      U3D_GRID%PERIODIC_PARTNER_FACE = -1
      U3D_GRID%PERIODIC_PARTNER_GROUP = -1
      U3D_GRID%PERIODIC_VERTEX_PERM = 0
      U3D_GRID%PERIODIC_TRANSLATION = 0.d0

      N_BOUNDARY_FACES = COUNT(U3D_GRID%CELL_FACES_PG > 0)
      IF (N_BOUNDARY_FACES == 0) THEN
         ERROR_CODE = PERIODIC_MAP_INVALID_GRID
         ERROR_MESSAGE = '3D periodic pairing found no tagged boundary triangles.'
         RETURN
      END IF
      ALLOCATE(CANDIDATE_CELLS(N_BOUNDARY_FACES), CANDIDATE_FACES(N_BOUNDARY_FACES))
      ALLOCATE(CANDIDATE_GROUPS(N_BOUNDARY_FACES), CANDIDATE_ORDER(N_BOUNDARY_FACES))
      ALLOCATE(CANDIDATE_AXIS_COORDS(N_BOUNDARY_FACES))
      ALLOCATE(USED_PARTNER_FACE(4,U3D_GRID%NUM_CELLS))
      USED_PARTNER_FACE = .FALSE.

      DO IPG = 1, N_GRID_BC
         IF (.NOT. IS_MASTER_GROUP(IPG)) CYCLE

         IF (.NOT. ALL(GRID_BC(IPG)%PARTICLE_BC == PERIODIC_MASTER) .OR. &
             GRID_BC(IPG)%FIELD_BC /= PERIODIC_MASTER_BC) THEN
            ERROR_CODE = PERIODIC_MAP_INVALID_BOUNDARY
            WRITE(ERROR_MESSAGE,'(A,A,A)') 'Periodic master group has inconsistent particle/field BCs: ', &
                                           TRIM(GRID_BC(IPG)%PHYSICAL_GROUP_NAME), '.'
            RETURN
         END IF

         TRANSLATION = GRID_BC(IPG)%TRANSLATEVEC
         IF (.NOT. ALL(IEEE_IS_FINITE(TRANSLATION)) .OR. NORM2(TRANSLATION) <= TOLERANCE) THEN
            ERROR_CODE = PERIODIC_MAP_BAD_TRANSLATION
            WRITE(ERROR_MESSAGE,'(A,A,A)') 'Periodic master group has a zero or non-finite translation: ', &
                                           TRIM(GRID_BC(IPG)%PHYSICAL_GROUP_NAME), '.'
            RETURN
         END IF

         MASTER_FACE_COUNT = COUNT(U3D_GRID%CELL_FACES_PG == IPG)
         IF (MASTER_FACE_COUNT == 0) THEN
            ERROR_CODE = PERIODIC_MAP_INVALID_BOUNDARY
            WRITE(ERROR_MESSAGE,'(A,A,A)') 'Periodic master group has no mapped faces: ', &
                                           TRIM(GRID_BC(IPG)%PHYSICAL_GROUP_NAME), '.'
            RETURN
         END IF

         AXIS = MAXLOC(ABS(TRANSLATION), DIM=1)
         N_CANDIDATES = 0
         DO IC = 1, U3D_GRID%NUM_CELLS
            DO IFACE = 1, 4
               PG = U3D_GRID%CELL_FACES_PG(IFACE,IC)
               IF (PG <= 0) CYCLE
               IF (PG > N_GRID_BC) THEN
                  ERROR_CODE = PERIODIC_MAP_INVALID_GRID
                  ERROR_MESSAGE = '3D boundary-face map references a physical group outside GRID_BC.'
                  RETURN
               END IF
               IF (IS_MASTER_GROUP(PG)) CYCLE
               N_CANDIDATES = N_CANDIDATES + 1
               CANDIDATE_CELLS(N_CANDIDATES) = IC
               CANDIDATE_FACES(N_CANDIDATES) = IFACE
               CANDIDATE_GROUPS(N_CANDIDATES) = PG
               CANDIDATE_CENTROID = SUM(U3D_GRID%NODE_COORDS(:, &
                  U3D_GRID%FACE_NODES(:,IFACE,IC)), DIM=2) / 3.d0
               CANDIDATE_AXIS_COORDS(N_CANDIDATES) = CANDIDATE_CENTROID(AXIS)
               CANDIDATE_ORDER(N_CANDIDATES) = N_CANDIDATES
            END DO
         END DO

         IF (N_CANDIDATES > 1) CALL SORT_PERIODIC_CANDIDATES(CANDIDATE_AXIS_COORDS, &
                                                              CANDIDATE_ORDER, 1, N_CANDIDATES)
         SLAVE_GROUP = -1

         DO MASTER_CELL = 1, U3D_GRID%NUM_CELLS
            DO MASTER_FACE = 1, 4
               IF (U3D_GRID%CELL_FACES_PG(MASTER_FACE,MASTER_CELL) /= IPG) CYCLE

               MASTER_COORDS = U3D_GRID%NODE_COORDS(:, &
                  U3D_GRID%FACE_NODES(:,MASTER_FACE,MASTER_CELL))
               MASTER_CENTROID = SUM(MASTER_COORDS, DIM=2) / 3.d0
               EXPECTED_CENTROID_AXIS = MASTER_CENTROID(AXIS) + TRANSLATION(AXIS)
               LOWER = 1
               UPPER = N_CANDIDATES + 1
               DO WHILE (LOWER < UPPER)
                  MIDPOINT = (LOWER + UPPER) / 2
                  IF (MIDPOINT <= N_CANDIDATES) THEN
                     IF (CANDIDATE_AXIS_COORDS(CANDIDATE_ORDER(MIDPOINT)) < &
                         EXPECTED_CENTROID_AXIS - TOLERANCE) THEN
                        LOWER = MIDPOINT + 1
                     ELSE
                        UPPER = MIDPOINT
                     END IF
                  ELSE
                     UPPER = MIDPOINT
                  END IF
               END DO

               MATCH_COUNT = 0
               MATCHED_CELL = -1
               MATCHED_FACE = -1
               MATCHED_GROUP = -1
               MATCHED_PERMUTATION = 0
               DO I = LOWER, N_CANDIDATES
                  CANDIDATE_INDEX = CANDIDATE_ORDER(I)
                  IF (CANDIDATE_AXIS_COORDS(CANDIDATE_INDEX) > &
                      EXPECTED_CENTROID_AXIS + TOLERANCE) EXIT

                  CANDIDATE_CELL = CANDIDATE_CELLS(CANDIDATE_INDEX)
                  CANDIDATE_FACE = CANDIDATE_FACES(CANDIDATE_INDEX)
                  CANDIDATE_GROUP = CANDIDATE_GROUPS(CANDIDATE_INDEX)
                  CANDIDATE_COORDS = U3D_GRID%NODE_COORDS(:, &
                     U3D_GRID%FACE_NODES(:,CANDIDATE_FACE,CANDIDATE_CELL))
                  CANDIDATE_CENTROID = SUM(CANDIDATE_COORDS, DIM=2) / 3.d0
                  IF (NORM2(CANDIDATE_CENTROID - MASTER_CENTROID - TRANSLATION) > TOLERANCE) CYCLE

                  CALL MATCH_PERIODIC_TRIANGLE(MASTER_COORDS, CANDIDATE_COORDS, TRANSLATION, &
                                               TOLERANCE, VERTEX_PERMUTATION, GEOMETRY_STATUS)
                  IF (GEOMETRY_STATUS == PERIODIC_MESH_NO_MATCH) CYCLE
                  IF (GEOMETRY_STATUS == PERIODIC_MESH_AMBIGUOUS) THEN
                     ERROR_CODE = PERIODIC_MAP_AMBIGUOUS
                     WRITE(ERROR_MESSAGE,'(A,A,A,I0,A,I0,A)') 'Ambiguous vertices on candidate periodic face for group ', &
                        TRIM(GRID_BC(IPG)%PHYSICAL_GROUP_NAME), ', cell=', CANDIDATE_CELL, ', face=', CANDIDATE_FACE, '.'
                     RETURN
                  ELSE IF (GEOMETRY_STATUS /= PERIODIC_MESH_OK) THEN
                     ERROR_CODE = PERIODIC_MAP_BAD_GEOMETRY
                     ERROR_MESSAGE = 'Invalid coordinates or tolerance while matching periodic triangles.'
                     RETURN
                  END IF

                  CANDIDATE_NORMAL = U3D_GRID%FACE_NORMAL(:,CANDIDATE_FACE,CANDIDATE_CELL)
                  MASTER_NORMAL = U3D_GRID%FACE_NORMAL(:,MASTER_FACE,MASTER_CELL)
                  IF (.NOT. IEEE_IS_FINITE(U3D_GRID%FACE_AREA(MASTER_FACE,MASTER_CELL)) .OR. &
                      .NOT. IEEE_IS_FINITE(U3D_GRID%FACE_AREA(CANDIDATE_FACE,CANDIDATE_CELL)) .OR. &
                      U3D_GRID%FACE_AREA(MASTER_FACE,MASTER_CELL) <= 0.d0 .OR. &
                      U3D_GRID%FACE_AREA(CANDIDATE_FACE,CANDIDATE_CELL) <= 0.d0) THEN
                     ERROR_CODE = PERIODIC_MAP_BAD_GEOMETRY
                     ERROR_MESSAGE = 'Periodic face has a non-finite or non-positive area.'
                     RETURN
                  END IF
                  AREA_TOLERANCE = ABS_TOL * L_REF + REL_TOL * MAX( &
                     U3D_GRID%FACE_AREA(MASTER_FACE,MASTER_CELL), &
                     U3D_GRID%FACE_AREA(CANDIDATE_FACE,CANDIDATE_CELL))
                  IF (ABS(U3D_GRID%FACE_AREA(MASTER_FACE,MASTER_CELL) - &
                          U3D_GRID%FACE_AREA(CANDIDATE_FACE,CANDIDATE_CELL)) > AREA_TOLERANCE) THEN
                     ERROR_CODE = PERIODIC_MAP_BAD_GEOMETRY
                     WRITE(ERROR_MESSAGE,'(A,A,A,I0,A,I0,A)') 'Periodic face area mismatch for group ', &
                        TRIM(GRID_BC(IPG)%PHYSICAL_GROUP_NAME), ', cell=', MASTER_CELL, ', face=', MASTER_FACE, '.'
                     RETURN
                  END IF

                  IF (.NOT. ALL(IEEE_IS_FINITE(MASTER_NORMAL)) .OR. &
                      .NOT. ALL(IEEE_IS_FINITE(CANDIDATE_NORMAL)) .OR. &
                      NORM2(MASTER_NORMAL) <= 0.d0 .OR. NORM2(CANDIDATE_NORMAL) <= 0.d0) THEN
                     ERROR_CODE = PERIODIC_MAP_BAD_GEOMETRY
                     ERROR_MESSAGE = 'Periodic face has an invalid outward normal.'
                     RETURN
                  END IF
                  NORMALIZED_MASTER = MASTER_NORMAL / NORM2(MASTER_NORMAL)
                  NORMALIZED_CANDIDATE = CANDIDATE_NORMAL / NORM2(CANDIDATE_NORMAL)
                  IF (DOT_PRODUCT(NORMALIZED_MASTER,NORMALIZED_CANDIDATE) > -1.d0 + 1.d-8) THEN
                     ERROR_CODE = PERIODIC_MAP_BAD_GEOMETRY
                     WRITE(ERROR_MESSAGE,'(A,A,A,I0,A,I0,A)') 'Periodic face normals are not opposite for group ', &
                        TRIM(GRID_BC(IPG)%PHYSICAL_GROUP_NAME), ', cell=', MASTER_CELL, ', face=', MASTER_FACE, '.'
                     RETURN
                  END IF

                  MATCH_COUNT = MATCH_COUNT + 1
                  MATCHED_CELL = CANDIDATE_CELL
                  MATCHED_FACE = CANDIDATE_FACE
                  MATCHED_GROUP = CANDIDATE_GROUP
                  MATCHED_PERMUTATION = VERTEX_PERMUTATION
                  IF (MATCH_COUNT > 1) THEN
                     ERROR_CODE = PERIODIC_MAP_AMBIGUOUS
                     WRITE(ERROR_MESSAGE,'(A,A,A,I0,A,I0,A)') 'Master face has multiple periodic partners in group ', &
                        TRIM(GRID_BC(IPG)%PHYSICAL_GROUP_NAME), ', cell=', MASTER_CELL, ', face=', MASTER_FACE, '.'
                     RETURN
                  END IF
               END DO

               IF (MATCH_COUNT == 0) THEN
                  ERROR_CODE = PERIODIC_MAP_NO_PARTNER
                  WRITE(ERROR_MESSAGE,'(A,A,A,I0,A,I0,A)') 'No translated partner for periodic group ', &
                     TRIM(GRID_BC(IPG)%PHYSICAL_GROUP_NAME), ', cell=', MASTER_CELL, ', face=', MASTER_FACE, '.'
                  RETURN
               END IF

               IF (SLAVE_GROUP == -1) THEN
                  SLAVE_GROUP = MATCHED_GROUP
               ELSE IF (SLAVE_GROUP /= MATCHED_GROUP) THEN
                  ERROR_CODE = PERIODIC_MAP_GROUP_CONFLICT
                  WRITE(ERROR_MESSAGE,'(A,A,A)') 'One periodic master group maps to more than one slave group: ', &
                                                 TRIM(GRID_BC(IPG)%PHYSICAL_GROUP_NAME), '.'
                  RETURN
               END IF

               IF ((IS_SLAVE_GROUP(MATCHED_GROUP) .AND. SLAVE_GROUP /= MATCHED_GROUP) .OR. &
                   USED_PARTNER_FACE(MATCHED_FACE,MATCHED_CELL)) THEN
                  ERROR_CODE = PERIODIC_MAP_AMBIGUOUS
                  WRITE(ERROR_MESSAGE,'(A,A,A,I0,A,I0,A)') 'Periodic partner face is reused for master group ', &
                     TRIM(GRID_BC(IPG)%PHYSICAL_GROUP_NAME), ', cell=', MASTER_CELL, ', face=', MASTER_FACE, '.'
                  RETURN
               END IF

               IF (.NOT. (ALL(GRID_BC(MATCHED_GROUP)%PARTICLE_BC == VACUUM) .OR. &
                          ALL(GRID_BC(MATCHED_GROUP)%PARTICLE_BC == PERIODIC_SLAVE)) .OR. &
                   .NOT. (GRID_BC(MATCHED_GROUP)%FIELD_BC == NO_BC .OR. &
                          GRID_BC(MATCHED_GROUP)%FIELD_BC == PERIODIC_SLAVE_BC) .OR. &
                   GRID_BC(MATCHED_GROUP)%REACT) THEN
                  ERROR_CODE = PERIODIC_MAP_GROUP_CONFLICT
                  WRITE(ERROR_MESSAGE,'(A,A,A)') 'Geometric periodic partner has a conflicting boundary condition: ', &
                                                 TRIM(GRID_BC(MATCHED_GROUP)%PHYSICAL_GROUP_NAME), '.'
                  RETURN
               END IF

               USED_PARTNER_FACE(MATCHED_FACE,MATCHED_CELL) = .TRUE.
               IS_SLAVE_GROUP(MATCHED_GROUP) = .TRUE.
               U3D_GRID%PERIODIC_PARTNER_CELL(MASTER_FACE,MASTER_CELL) = MATCHED_CELL
               U3D_GRID%PERIODIC_PARTNER_FACE(MASTER_FACE,MASTER_CELL) = MATCHED_FACE
               U3D_GRID%PERIODIC_PARTNER_GROUP(MASTER_FACE,MASTER_CELL) = MATCHED_GROUP
               U3D_GRID%PERIODIC_VERTEX_PERM(:,MASTER_FACE,MASTER_CELL) = MATCHED_PERMUTATION
               U3D_GRID%PERIODIC_TRANSLATION(:,MASTER_FACE,MASTER_CELL) = TRANSLATION

               DO I = 1, 3
                  REVERSE_PERMUTATION(MATCHED_PERMUTATION(I)) = I
               END DO
               U3D_GRID%PERIODIC_PARTNER_CELL(MATCHED_FACE,MATCHED_CELL) = MASTER_CELL
               U3D_GRID%PERIODIC_PARTNER_FACE(MATCHED_FACE,MATCHED_CELL) = MASTER_FACE
               U3D_GRID%PERIODIC_PARTNER_GROUP(MATCHED_FACE,MATCHED_CELL) = IPG
               U3D_GRID%PERIODIC_VERTEX_PERM(:,MATCHED_FACE,MATCHED_CELL) = REVERSE_PERMUTATION
               U3D_GRID%PERIODIC_TRANSLATION(:,MATCHED_FACE,MATCHED_CELL) = -TRANSLATION
               U3D_GRID%NUM_PERIODIC_FACE_PAIRS = U3D_GRID%NUM_PERIODIC_FACE_PAIRS + 1
            END DO
         END DO

         SLAVE_FACE_COUNT = COUNT(U3D_GRID%CELL_FACES_PG == SLAVE_GROUP)
         IF (SLAVE_FACE_COUNT /= MASTER_FACE_COUNT) THEN
            ERROR_CODE = PERIODIC_MAP_NO_PARTNER
            WRITE(ERROR_MESSAGE,'(A,A,A,A,A,I0,A,I0)') 'Periodic group face counts differ: master ', &
               TRIM(GRID_BC(IPG)%PHYSICAL_GROUP_NAME), ', slave ', &
               TRIM(GRID_BC(SLAVE_GROUP)%PHYSICAL_GROUP_NAME), ', master faces=', MASTER_FACE_COUNT, &
               ', slave faces=', SLAVE_FACE_COUNT
            RETURN
         END IF
      END DO

      DO IPG = 1, N_GRID_BC
         IF (.NOT. IS_SLAVE_GROUP(IPG)) CYCLE
         GRID_BC(IPG)%PARTICLE_BC = PERIODIC_SLAVE
         GRID_BC(IPG)%FIELD_BC = PERIODIC_SLAVE_BC
      END DO

   END SUBROUTINE BUILD_PERIODIC_BOUNDARY_FACE_MAP


   RECURSIVE SUBROUTINE SORT_PERIODIC_CANDIDATES(KEYS, ORDER, LEFT, RIGHT)

      IMPLICIT NONE

      REAL(KIND=8), DIMENSION(:), INTENT(IN) :: KEYS
      INTEGER, DIMENSION(:), INTENT(INOUT) :: ORDER
      INTEGER, INTENT(IN) :: LEFT, RIGHT

      INTEGER :: I, J, TEMP
      REAL(KIND=8) :: PIVOT

      I = LEFT
      J = RIGHT
      PIVOT = KEYS(ORDER((LEFT + RIGHT) / 2))
      DO WHILE (I <= J)
         DO WHILE (KEYS(ORDER(I)) < PIVOT)
            I = I + 1
         END DO
         DO WHILE (KEYS(ORDER(J)) > PIVOT)
            J = J - 1
         END DO
         IF (I <= J) THEN
            TEMP = ORDER(I)
            ORDER(I) = ORDER(J)
            ORDER(J) = TEMP
            I = I + 1
            J = J - 1
         END IF
      END DO

      IF (LEFT < J) CALL SORT_PERIODIC_CANDIDATES(KEYS, ORDER, LEFT, J)
      IF (I < RIGHT) CALL SORT_PERIODIC_CANDIDATES(KEYS, ORDER, I, RIGHT)

   END SUBROUTINE SORT_PERIODIC_CANDIDATES

   ! SUBROUTINE ASSIGN_CELLS_TO_PROCS

   !    INTEGER, DIMENSION(:), ALLOCATABLE :: NPC, NPCMOD, NPP
   !    INTEGER :: JP, IC, NPPDESIRED, IPROC
      
   !    ALLOCATE(NPC(NCELLS))
   !    NPC = 0

   !    DO JP = 1, NP_PROC
   !       IC = particles(JP)%IC
   !       NPC(IC) = NPC(IC) + 1
   !    END DO

   !    IF (PROC_ID .EQ. 0) THEN
   !       CALL MPI_REDUCE(MPI_IN_PLACE, NPC, NCELLS, MPI_INTEGER,  MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   !    ELSE
   !       CALL MPI_REDUCE(NPC,          NPC, NCELLS, MPI_INTEGER,  MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   !    END IF

   !    CALL MPI_BCAST(NPC, NCELLS, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)


   !    NPCMOD = NPC + SUM(NPC)/NCELLS/10 + 1

   !    NPPDESIRED = SUM(NPCMOD) / N_MPI_THREADS

   !    ALLOCATE(NPP(N_MPI_THREADS))
   !    NPP = 0

   !    IF (.NOT. ALLOCATED(CELL_PROCS)) ALLOCATE(CELL_PROCS(NCELLS))
   !    CELL_PROCS = 0
   !    IPROC = 0
   !    DO IC = 1, NCELLS
   !       IF (NPP(IPROC) < NPPDESIRED .OR. IPROC == N_MPI_THREADS - 1) THEN
   !          NPP(IPROC) = NPP(IPROC) + NPCMOD(IC)
   !          CELL_PROCS(IC) = IPROC
   !       ELSE
   !          IPROC = IPROC + 1
   !          NPP(IPROC) = NPP(IPROC) + NPCMOD(IC)
   !          CELL_PROCS(IC) = IPROC
   !       END IF
   !    END DO

   !    DEALLOCATE(NPC)
   !    DEALLOCATE(NPP)


   ! END SUBROUTINE ASSIGN_CELLS_TO_PROCS



   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! ! SUBROUTINE PROC_FROM_POSITION -> finds the process ID from particle position              !
   ! ! Note this depends on the parallelization strategy                                         !
   ! ! This should never be called when an unstructured grid is used.                            !
   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   ! SUBROUTINE PROC_FROM_POSITION(XP,YP,  IDPROC)

   !    ! Note: processes go from 0 (usually termed the Master) to MPI_N_THREADS - 1
   !    ! If the number of MPI processes is only 1, then IDPROC is 0, the one and only process.
   !    ! Otherwise, check the partition style (variable "DOMPART_TYPE")

   !    IMPLICIT NONE

   !    REAL(KIND=8), INTENT(IN) :: XP, YP
   !    INTEGER, INTENT(OUT)     :: IDPROC

   !    INTEGER :: NCELLSPP
   !    INTEGER :: IDCELL

   !    IF (N_MPI_THREADS == 1) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Serial operation

   !      IDPROC = 0      

   !    ELSE IF (DOMPART_TYPE == 0) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ S T R I P S
   !    !
   !    ! Domain partition type: strips
   !    !
   !    ! This domain partition assigns the same number of grid cells to every process.
   !    ! The domain partition is thus in strips along the "x" axis, as the cells are assumed 
   !    ! to be numbered along x.
 
   !       ! 1) Find ID of cell where the particle is
   !       CALL CELL_FROM_POSITION(XP, YP, IDCELL)
   !       IF (IDCELL .GT. NX*NY .OR. IDCELL .LT. 1) THEN
   !          WRITE(*,*) 'Error! CELL_FROM_POSITION returned cell:', IDCELL, 'Particle position: ', XP, ', ', YP
   !       END IF

   !       ! 2) Find number of cells for each process (NX*NY*NZ = number of cells)
   !       !    Exceed a bit, so the last processor has slightly less cells, if number
   !       !    of cells is not divisible by the MPI_threads
   !       NCELLSPP = CEILING(REAL(NX*NY)/REAL(N_MPI_THREADS))

   !       ! 3) Here is the process ID 
   !       IDPROC   = INT((IDCELL-1)/NCELLSPP) ! Before was giving wrong result for peculiar combinations of NCELLS and NCELLSPP
   !       IF (IDPROC .GT. N_MPI_THREADS-1 .OR. IDPROC .LT. 0) THEN
   !          WRITE(*,*) 'Error! PROC_FROM_POSITION returned proc:', IDPROC
   !       END IF


   !    ELSE IF (DOMPART_TYPE == 1) THEN ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ B L O C K S
   !    !
   !    ! Domain partition type: blocks
   !    !
   !    ! This domain partition style divides the domain into blocks
   !    ! Note that blocks in this definition are independent from cells! IT IS POSSIBLE TO
   !    ! GENERALIZE IT (and we could and should.. but for PIC I don't really care).

   !       IDPROC = INT((XP-XMIN)/DX_BLOCKS) + N_BLOCKS_X*INT((YP-YMIN)/DY_BLOCKS)

   !    ELSE ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ E R R O R

   !       WRITE(*,*) 'ATTENTION! Value for variable DOMPART_TYPE = ', DOMPART_TYPE
   !       WRITE(*,*) ' not recognized! Check input file!! In: PROC_FROM_POSITION() ABORTING!'
   !       STOP

   !    END IF

   ! END SUBROUTINE PROC_FROM_POSITION


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE COUNTING_SORT -> Sorts particles using the "counting sort"  !!
   ! sorting algorithm. Requires "PROC_FROM_POSITION" subroutine.           !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   SUBROUTINE COUNTING_SORT(PARTICLES_ARRAY, PARTICLES_COUNT)
   
      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), INTENT(INOUT) :: PARTICLES_ARRAY ! array
      INTEGER                      , DIMENSION(:), INTENT(IN)    :: PARTICLES_COUNT ! freq
      INTEGER                      , DIMENSION(:), ALLOCATABLE   :: disp1
      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE   :: sorted
      INTEGER                                                    :: i1
      INTEGER                                                    :: IPROC
   
      ALLOCATE(sorted(SIZE(PARTICLES_ARRAY)))
      ALLOCATE(disp1(0:SIZE(PARTICLES_COUNT)-1))
   
      disp1 = PARTICLES_COUNT

      DO i1 = 1, N_MPI_THREADS - 1
         disp1(i1) = disp1(i1) + disp1(i1-1)
      END DO
   
      DO i1 = SIZE(PARTICLES_ARRAY), 1, -1
   
         ! Find position of particle
         !CALL PROC_FROM_POSITION(PARTICLES_ARRAY(i1)%X, PARTICLES_ARRAY(i1)%Y,  IPROC)
         CALL PROC_FROM_CELL(particles_array(i1)%IC, IPROC)
   
         sorted(disp1(IPROC)) = PARTICLES_ARRAY(i1)
         disp1(IPROC)         = disp1(IPROC) - 1
   
      END DO
   
      PARTICLES_ARRAY = sorted
   
      DEALLOCATE(sorted)
      DEALLOCATE(disp1)
   
      RETURN
   
   END SUBROUTINE COUNTING_SORT


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! SUBROUTINE EXCHANGE -> Exchanges particles between processes !!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE EXCHANGE

   ! This subroutine exchanges particles among processes.
   ! 1) It counts how many particles, in the vector of particles, are to be sent
   ! 2) Copies them in the sendbuf vector and removes them from the "particles" array
   ! 3) Sends around the information on how many they are
   ! 4) Sends them with MPI_ALLTOALLV
   ! 5) Copies them at the bottom of the "particles" array

      IMPLICIT NONE

      TYPE(PARTICLE_DATA_STRUCTURE), DIMENSION(:), ALLOCATABLE :: sendbuf, recvbuf

      INTEGER, DIMENSION(:), ALLOCATABLE :: sendcount, recvcount
      INTEGER, DIMENSION(:), ALLOCATABLE :: senddispl, recvdispl

      INTEGER :: i, IP, JP, JPROC, NP_RECV

      ! Keep a valid base address even on ranks with zero particles.  MPI
      ! implementations may still inspect buffer arguments for zero counts.
      ALLOCATE(sendbuf(MAX(1,NP_PROC))) ! At most NP_PROC particles are exchanged.

      ALLOCATE(sendcount(N_MPI_THREADS)) ! Allocate other vectors
      ALLOCATE(recvcount(N_MPI_THREADS))
      ALLOCATE(senddispl(N_MPI_THREADS))
      ALLOCATE(recvdispl(N_MPI_THREADS))

      DO i = 1,N_MPI_THREADS
         sendcount(i) = 0
         recvcount(i) = 0
         senddispl(i) = 0
         recvdispl(i) = 0
      END DO

      ! Loop on particles, from the last one to the first one and check if they belong to current process
      
      IP = NP_PROC ! Init 
      JP = 0       ! Init

      DO WHILE( IP .GE. 1 )

         !CALL PROC_FROM_POSITION(particles(IP)%X, particles(IP)%Y, JPROC) ! Find which processor the particle belongs to
         CALL PROC_FROM_CELL(particles(IP)%IC, JPROC)

         IF (JPROC .NE. PROC_ID) THEN ! I shall send it to processor JPROC

            ! Increment the number of particles that shall be sent to processor JPROC 
            sendcount(JPROC+1) = sendcount(JPROC+1) + 1 ! Note: processors ID start from 0. Index ID from 1

            ! Copy particle in the send buffer
            JP = JP + 1
            sendbuf(JP) = particles(IP) 

            ! Remove particle from current array
            CALL REMOVE_PARTICLE_ARRAY(IP, particles, NP_PROC)

         END IF

         IP = IP - 1

      END DO

      CALL COUNTING_SORT(sendbuf(1:JP), sendcount) ! Reorder send buffer by process ID

      ! ~~~~~~ At this point, exchange particles among processes ~~~~~~

      ! Say how many particles are to be received
      CALL MPI_ALLTOALL(sendcount, 1, MPI_INTEGER, recvcount, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

      NP_RECV = SUM(recvcount)
      ALLOCATE(recvbuf(MAX(1,NP_RECV)))

      ! Compute position of particle chunks to be sent & received (incremental)
      DO i = 1, N_MPI_THREADS-1
         senddispl(i+1) = senddispl(i) + sendcount(i)
         recvdispl(i+1) = recvdispl(i) + recvcount(i)
      END DO

      ! Compute position of particle chunks to be sent & received (incremental), and send them 
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
      CALL MPI_ALLTOALLV(sendbuf, sendcount, senddispl, MPI_PARTICLE_DATA_STRUCTURE, &
                         recvbuf, recvcount, recvdispl, MPI_PARTICLE_DATA_STRUCTURE, &
                         MPI_COMM_WORLD, ierr)

      ! At this point, write the particles received in the particles vector. Note that among these received,
      ! the processor also has its own particles.

      DO IP = 1, NP_RECV
         CALL ADD_PARTICLE_ARRAY(recvbuf(IP), NP_PROC, particles)
      END DO

      ! ~~~~~~~~ Done. Now deallocate stuff that would waste memory ~~~~~~~
      DEALLOCATE(sendbuf) 
      DEALLOCATE(recvbuf) 
      DEALLOCATE(sendcount)
      DEALLOCATE(recvcount)
      DEALLOCATE(senddispl)
      DEALLOCATE(recvdispl)
      
   END SUBROUTINE EXCHANGE




   SUBROUTINE READ_1D_UNSTRUCTURED_GRID_SU2(FILENAME)

      IMPLICIT NONE

      CHARACTER*256, INTENT(IN) :: FILENAME

      CHARACTER*256 :: LINE, GROUPNAME, DUMMYLINE

      INTEGER, PARAMETER :: in5 = 2385
      INTEGER            :: ios
      INTEGER            :: ReasonEOF

      INTEGER            :: NUM, I, J, FOUND, V1, V2, ELEM_TYPE, NUMELEMS
      REAL(KIND=8)       :: X1, X2
      REAL(KIND=8), DIMENSION(3) :: XYZ, A, B

      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: TEMP_CELL_NEIGHBORS

      INTEGER, DIMENSION(2) :: VLIST2

      INTEGER, DIMENSION(:), ALLOCATABLE :: N_CELLS_WITH_NODE, CELL_WITH_NODE, IOF
      INTEGER :: IDX, JN, JC1, JC2

      LOGICAL, DIMENSION(:), ALLOCATABLE :: NODE_ON_BOUNDARY
      INTEGER :: NUM_BOUNDARY_NODES, NUM_BOUNDARY_ELEM

      ! Open input file for reading
      OPEN(UNIT=in5,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios .NE. 0) THEN
         CALL ERROR_ABORT('Attention, mesh file not found! ABORTING.')
      ENDIF

      ! Read the mesh file. SU2 file format (*.su2)
      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Reading grid file in SU2 format.'
         WRITE(*,*) '==========================================='
      END IF
      
      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE, NUM
         IF (ReasonEOF < 0) EXIT 
         !WRITE(*,*) 'Read line:', LINE, ' number ', NUM
         
         IF (LINE == 'NPOIN=') THEN
            ALLOCATE(U1D_GRID%NODE_COORDS(3,NUM))
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) XYZ(1)
               XYZ(2) = 0.d0
               XYZ(3) = 0.d0 ! Stay on the x axis.
               U1D_GRID%NODE_COORDS(:,I) = XYZ
            END DO
            U1D_GRID%NUM_NODES = NUM
         ELSE IF (LINE == 'NELEM=') THEN
            ALLOCATE(U1D_GRID%CELL_NODES(2,NUM))

            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) ELEM_TYPE, U1D_GRID%CELL_NODES(:,I)
               !WRITE(*,*) 'I read element ', I, ' has nodes ', U1D_GRID%CELL_NODES(:,I)
               IF (ELEM_TYPE .NE. 3) WRITE(*,*) 'Element type was not 3!'
            END DO
            U1D_GRID%CELL_NODES = U1D_GRID%CELL_NODES + 1 ! Start indexing from 1.

            U1D_GRID%NUM_CELLS = NUM

            ALLOCATE(U1D_GRID%CELL_EDGES_PG(2, U1D_GRID%NUM_CELLS))
            U1D_GRID%CELL_EDGES_PG = -1

            ALLOCATE(U1D_GRID%CELL_PG(U1D_GRID%NUM_CELLS))
            U1D_GRID%CELL_PG = -1
      
         ELSE IF (LINE == 'NMARK=') THEN

            ! Assign physical groups to cell edges.
            DO I = 1, NUM
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, GROUPNAME
               IF (LINE .NE. 'MARKER_TAG=') THEN
                  WRITE(*,*) 'Error! did not find marker name.'
               ELSE
                  !WRITE(*,*) 'Found marker tag, with groupname: ', GROUPNAME
               END IF
         
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, NUMELEMS
               IF (LINE .NE. 'MARKER_ELEMS=') THEN
                  WRITE(*,*) 'Error! did not find marker elements.'
               ELSE
                  !WRITE(*,*) 'Found marker elements, with number of elements: ', NUMELEMS
               END IF

               DO J = 1, NUMELEMS
                  READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
               END DO
            END DO

         END IF
      END DO

      REWIND(in5)

      ALLOCATE(N_CELLS_WITH_NODE(U1D_GRID%NUM_NODES))
      ALLOCATE(IOF(U1D_GRID%NUM_NODES))

      N_CELLS_WITH_NODE = 0
      DO I = 1, U1D_GRID%NUM_CELLS
         DO V1 = 1, 2
            JN = U1D_GRID%CELL_NODES(V1,I)
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO
   
      IOF = -1
      IDX = 1
      DO JN = 1, U1D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) .NE. 0) THEN
            IOF(JN) = IDX
            IDX = IDX + N_CELLS_WITH_NODE(JN)
         END IF
      END DO
   
      ALLOCATE(CELL_WITH_NODE(IDX))
      
      N_CELLS_WITH_NODE = 0
      DO I = 1, U1D_GRID%NUM_CELLS
         DO V1 = 1, 2
            JN = U1D_GRID%CELL_NODES(V1,I)
            CELL_WITH_NODE(IOF(JN) + N_CELLS_WITH_NODE(JN)) = I
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO


      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE, NUM
         IF (ReasonEOF < 0) EXIT 
         !WRITE(*,*) 'Read line:', LINE, ' number ', NUM
         
         IF (LINE == 'NPOIN=') THEN
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
            END DO
         ELSE IF (LINE == 'NELEM=') THEN
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
            END DO
         ELSE IF (LINE == 'NMARK=') THEN

            ! Assign physical groups to cells/cell edges.

            ALLOCATE(GRID_BC(NUM)) ! Append the physical group to the list
            N_GRID_BC = NUM

            DO I = 1, NUM
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, GROUPNAME
      
               GRID_BC(I)%PHYSICAL_GROUP_NAME = GROUPNAME
         
               READ(in5,*, IOSTAT=ReasonEOF) LINE, NUMELEMS

               DO J = 1, NUMELEMS
                  READ(in5,'(A)', IOSTAT=ReasonEOF) LINE

                  READ(LINE,*) ELEM_TYPE

                  IF (ELEM_TYPE == 1) THEN ! element in physical group is a vertex.

                     READ(LINE,*) ELEM_TYPE, JN

                     JN = JN + 1
                     
                     IF (N_CELLS_WITH_NODE(JN) > 0) THEN
                        DO IDX = 0, N_CELLS_WITH_NODE(JN) - 1
                           JC1 = CELL_WITH_NODE(IOF(JN) + IDX)
                           FOUND = 0
                           DO V1 = 1, 2
                              IF (U1D_GRID%CELL_NODES(V1,JC1) == JN) THEN
                                 U1D_GRID%CELL_EDGES_PG(V1, JC1) = I
                              END IF
                           END DO
                        END DO
                     END IF

                  ELSE IF (ELEM_TYPE == 3) THEN ! element in physical group is a line.
                     READ(LINE,*) ELEM_TYPE, VLIST2
                     
                     VLIST2 = VLIST2 + 1

                     JN = VLIST2(1)
                     IF (N_CELLS_WITH_NODE(JN) > 0) THEN
                        DO IDX = 0, N_CELLS_WITH_NODE(JN) - 1
                           JC1 = CELL_WITH_NODE(IOF(JN) + IDX)
                           FOUND = 0
                           DO V1 = 1, 2
                              IF (ANY(VLIST2 == U1D_GRID%CELL_NODES(V1,JC1))) THEN
                                 FOUND = FOUND + 1
                              END IF
                           END DO
            
                           IF (FOUND == 2) THEN
                              U1D_GRID%CELL_PG(JC1) = I
                           END IF
                        END DO
                     END IF

                  ELSE
                     WRITE(*,*) 'Error! element type was not point or line.'
                  END IF

               END DO
            END DO

         END IF
      END DO

      ! Done reading
      CLOSE(in5)

      !WRITE(*,*) 'Read grid file. It contains ', U1D_GRID%NUM_NODES, &
      !           'points, and ', U1D_GRID%NUM_CELLS, 'cells.'

      ! Process the mesh: generate connectivity, normals and such...
      !XMIN, XMAX,...

      !DO I = 1, U1D_GRID%NUM_CELLS
      !   WRITE(*,*) U1D_GRID%CELL_NODES(:,I)
      !END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing cell volumes.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute cell volumes
      ALLOCATE(U1D_GRID%SEGMENT_LENGTHS(U1D_GRID%NUM_CELLS))
      ALLOCATE(U1D_GRID%CELL_VOLUMES(U1D_GRID%NUM_CELLS))
      DO I = 1, U1D_GRID%NUM_CELLS
         A = U1D_GRID%NODE_COORDS(:, U1D_GRID%CELL_NODES(1,I))
         B = U1D_GRID%NODE_COORDS(:, U1D_GRID%CELL_NODES(2,I))

         U1D_GRID%SEGMENT_LENGTHS(I) = ABS(A(1)-B(1))
         U1D_GRID%CELL_VOLUMES(I) = U1D_GRID%SEGMENT_LENGTHS(I) * (YMAX-YMIN) * (ZMAX-ZMIN)
      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing grid connectivity.'
         WRITE(*,*) '==========================================='
      END IF

      ! Find cell connectivity
      ALLOCATE(TEMP_CELL_NEIGHBORS(2, U1D_GRID%NUM_CELLS))
      TEMP_CELL_NEIGHBORS = -1



      DO JN = 1, U1D_GRID%NUM_NODES
         !IF (PROC_ID == 0) WRITE(*,*) 'Checking node ', JN, ' of ',  U1D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) > 1) THEN
            DO I = 0, N_CELLS_WITH_NODE(JN) - 1
               DO J = I, N_CELLS_WITH_NODE(JN) - 1
                  IF (I == J) CYCLE
                  JC1 = CELL_WITH_NODE(IOF(JN) + I)
                  JC2 = CELL_WITH_NODE(IOF(JN) + J)

                  FOUND = 0
                  DO V1 = 1, 2
                     DO V2 = 1, 2
                        IF (U1D_GRID%CELL_NODES(V1,JC1) == U1D_GRID%CELL_NODES(V2,JC2)) THEN
                           FOUND = FOUND + 1
                           IF (FOUND .GT. 1) CALL ERROR_ABORT('Error! Found duplicate cells in the mesh!')
                           TEMP_CELL_NEIGHBORS(V1, JC1) = JC2
                           TEMP_CELL_NEIGHBORS(V2, JC2) = JC1
                        END IF
                     END DO
                  END DO

               END DO
            END DO
         END IF
      END DO

      U1D_GRID%CELL_NEIGHBORS = TEMP_CELL_NEIGHBORS



      !WRITE(*,*) 'Generated grid connectivity. '
      !DO I = 1, U1D_GRID%NUM_CELLS
      !   WRITE(*,*) 'Cell ', I, ' neighbors cells ', TEMP_CELL_NEIGHBORS(:, I)
      !END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing face normals.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute segment edge normals
      ALLOCATE(U1D_GRID%EDGE_NORMAL(3, 2, U1D_GRID%NUM_CELLS))
      DO I = 1, U1D_GRID%NUM_CELLS

         U1D_GRID%EDGE_NORMAL(1,1,I) = -1.d0
         U1D_GRID%EDGE_NORMAL(2,1,I) =  0.d0
         U1D_GRID%EDGE_NORMAL(3,1,I) =  0.d0

         U1D_GRID%EDGE_NORMAL(1,2,I) =  1.d0
         U1D_GRID%EDGE_NORMAL(2,2,I) =  0.d0
         U1D_GRID%EDGE_NORMAL(3,2,I) =  0.d0

      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Checking ordering.'
         WRITE(*,*) '==========================================='
      END IF

      DO I = 1, U1D_GRID%NUM_CELLS
         X1 = U1D_GRID%NODE_COORDS(1, U1D_GRID%CELL_NODES(2,I)) &
            - U1D_GRID%NODE_COORDS(1, U1D_GRID%CELL_NODES(1,I))

         IF (X1 < 0) CALL ERROR_ABORT('1D mesh segment are reversed.')
      END DO

      NCELLS = U1D_GRID%NUM_CELLS
      NNODES = U1D_GRID%NUM_NODES



      ALLOCATE(U1D_GRID%BASIS_COEFFS(2,2,NCELLS))

      DO I = 1, NCELLS
         V1 = U1D_GRID%CELL_NODES(1,I)
         V2 = U1D_GRID%CELL_NODES(2,I)

         X1 = U1D_GRID%NODE_COORDS(1, V1)
         X2 = U1D_GRID%NODE_COORDS(1, V2)

         ! These are such that PSI_i = x * BASIS_COEFFS(1,i,IC) + BASIS_COEFFS(2,i,IC)

         U1D_GRID%BASIS_COEFFS(1,1,I) = -1.d0
         U1D_GRID%BASIS_COEFFS(2,1,I) =  X2

         U1D_GRID%BASIS_COEFFS(1,2,I) =  1.d0
         U1D_GRID%BASIS_COEFFS(2,2,I) = -X1
         

         U1D_GRID%BASIS_COEFFS(:,:,I) = U1D_GRID%BASIS_COEFFS(:,:,I)/U1D_GRID%SEGMENT_LENGTHS(I)

      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Creating boundary grid.'
         WRITE(*,*) '==========================================='
      END IF

      ALLOCATE(U1D_GRID%SEGMENT_NODES_BOUNDARY_INDEX(2,NCELLS))
      U1D_GRID%SEGMENT_NODES_BOUNDARY_INDEX = -1
      ALLOCATE(NODE_ON_BOUNDARY(NNODES))
      NODE_ON_BOUNDARY = .FALSE.
      ALLOCATE(U1D_GRID%NODES_BOUNDARY_INDEX(NNODES))
      U1D_GRID%NODES_BOUNDARY_INDEX = -1
      NUM_BOUNDARY_NODES = 0
      NUM_BOUNDARY_ELEM = 0
      DO I = 1, NCELLS
         DO J = 1, 2
            ! If the vertex belongs to any physical group, it should be part of the boundary grid
            ! Later, we may want to filter this further
            IF (U1D_GRID%CELL_EDGES_PG(J,I) .NE. -1) THEN
               NUM_BOUNDARY_ELEM = NUM_BOUNDARY_ELEM + 1
               V1 = U1D_GRID%CELL_NODES(J, I)
               IF (.NOT. NODE_ON_BOUNDARY(V1)) THEN
                  NUM_BOUNDARY_NODES = NUM_BOUNDARY_NODES + 1
                  U1D_GRID%NODES_BOUNDARY_INDEX(V1) = NUM_BOUNDARY_NODES
                  NODE_ON_BOUNDARY(V1) = .TRUE.
               END IF
            END IF
         END DO
      END DO

      U0D_GRID%NUM_POINTS = NUM_BOUNDARY_ELEM
      U0D_GRID%NUM_NODES = NUM_BOUNDARY_NODES
      ALLOCATE(U0D_GRID%POINT_NODES(NUM_BOUNDARY_ELEM))
      ALLOCATE(U0D_GRID%POINT_PG(NUM_BOUNDARY_ELEM))
      ALLOCATE(U0D_GRID%NODE_COORDS(3, NUM_BOUNDARY_NODES))

      DO I = 1, NNODES
         IF (NODE_ON_BOUNDARY(I)) THEN
            U0D_GRID%NODE_COORDS(:,U1D_GRID%NODES_BOUNDARY_INDEX(I)) = U1D_GRID%NODE_COORDS(:,I)
         END IF
      END DO

      NUM_BOUNDARY_ELEM = 0

      DO I = 1, NCELLS
         DO J = 1, 2
            IF (U1D_GRID%CELL_EDGES_PG(J,I) .NE. -1) THEN
               NUM_BOUNDARY_ELEM = NUM_BOUNDARY_ELEM + 1
               U0D_GRID%POINT_PG(NUM_BOUNDARY_ELEM) = U1D_GRID%CELL_EDGES_PG(J,I)
               U1D_GRID%SEGMENT_NODES_BOUNDARY_INDEX(J,I) = NUM_BOUNDARY_ELEM

               V1 = U1D_GRID%CELL_NODES(J, I)
               U0D_GRID%POINT_NODES(NUM_BOUNDARY_ELEM) = U1D_GRID%NODES_BOUNDARY_INDEX(V1)

            END IF
         END DO
      END DO
      
      DEALLOCATE(NODE_ON_BOUNDARY)

      NBOUNDCELLS = NUM_BOUNDARY_ELEM
      NBOUNDNODES = NUM_BOUNDARY_NODES

      ! Compute areas and lengths of boundary mesh
      ALLOCATE(U0D_GRID%VERTEX_AREAS(U0D_GRID%NUM_POINTS))
      U0D_GRID%VERTEX_AREAS = (YMAX-YMIN) * (ZMAX-ZMIN)

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '============================================================='
         WRITE(*,*) 'Done reading grid file.'
         WRITE(*,*) 'It contains ', NNODES, ' nodes and ', NCELLS, ' cells.'
         WRITE(*,*) 'The boundary grid contains ', NBOUNDCELLS, ' lines and ', NBOUNDNODES, ' nodes.'
         WRITE(*,*) '============================================================='
      END IF

   END SUBROUTINE READ_1D_UNSTRUCTURED_GRID_SU2


   SUBROUTINE READ_2D_UNSTRUCTURED_GRID_SU2(FILENAME)

      IMPLICIT NONE

      CHARACTER*256, INTENT(IN) :: FILENAME

      CHARACTER*256 :: LINE, GROUPNAME, DUMMYLINE

      INTEGER, PARAMETER :: in5 = 2385
      INTEGER            :: ios
      INTEGER            :: ReasonEOF

      INTEGER            :: NUM, I, J, FOUND, V1, V2, V3, ELEM_TYPE, NUMELEMS
      INTEGER, DIMENSION(3,2) :: IND
      REAL(KIND=8)       :: X1, X2, X3, Y1, Y2, Y3, LEN, RAD
      REAL(KIND=8), DIMENSION(3) :: XYZ, A, B, C

      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: TEMP_CELL_NEIGHBORS

      INTEGER, DIMENSION(2) :: VLIST2, WHICH1, WHICH2
      INTEGER, DIMENSION(3) :: VLIST3

      INTEGER, DIMENSION(:), ALLOCATABLE :: N_CELLS_WITH_NODE, CELL_WITH_NODE, IOF
      INTEGER :: IDX, JN, JC1, JC2

      LOGICAL, DIMENSION(:), ALLOCATABLE :: NODE_ON_BOUNDARY
      INTEGER :: NUM_BOUNDARY_NODES, NUM_BOUNDARY_ELEM

      ! Open input file for reading
      OPEN(UNIT=in5,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios .NE. 0) THEN
         CALL ERROR_ABORT('Attention, mesh file not found! ABORTING.')
      ENDIF

      ! Read the mesh file. SU2 file format (*.su2)
      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Reading grid file in SU2 format.'
         WRITE(*,*) '==========================================='
      END IF
      
      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE, NUM
         IF (ReasonEOF < 0) EXIT 
         !WRITE(*,*) 'Read line:', LINE, ' number ', NUM
         
         IF (LINE == 'NPOIN=') THEN
            ALLOCATE(U2D_GRID%NODE_COORDS(3,NUM))
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) XYZ(1:2)
               XYZ(3) = 0.d0 ! Stay in the x-y plane.
               U2D_GRID%NODE_COORDS(:,I) = XYZ
            END DO
            U2D_GRID%NUM_NODES = NUM
         ELSE IF (LINE == 'NELEM=') THEN
            ALLOCATE(U2D_GRID%CELL_NODES(3,NUM))

            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) ELEM_TYPE, U2D_GRID%CELL_NODES(:,I)
               !WRITE(*,*) 'I read element ', I, ' has nodes ', U2D_GRID%CELL_NODES(:,I)
               IF (ELEM_TYPE .NE. 5) WRITE(*,*) 'Element type was not 5!'
            END DO
            U2D_GRID%CELL_NODES = U2D_GRID%CELL_NODES + 1 ! Start indexing from 1.

            U2D_GRID%NUM_CELLS = NUM

            ALLOCATE(U2D_GRID%CELL_EDGES_PG(3, U2D_GRID%NUM_CELLS))
            U2D_GRID%CELL_EDGES_PG = -1

            ALLOCATE(U2D_GRID%CELL_PG(U2D_GRID%NUM_CELLS))
            U2D_GRID%CELL_PG = -1
      
         ELSE IF (LINE == 'NMARK=') THEN

            ! Assign physical groups to cell edges.
            DO I = 1, NUM
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, GROUPNAME
               IF (LINE .NE. 'MARKER_TAG=') THEN
                  WRITE(*,*) 'Error! did not find marker name.'
               ELSE
                  !WRITE(*,*) 'Found marker tag, with groupname: ', GROUPNAME
               END IF
         
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, NUMELEMS
               IF (LINE .NE. 'MARKER_ELEMS=') THEN
                  WRITE(*,*) 'Error! did not find marker elements.'
               ELSE
                  !WRITE(*,*) 'Found marker elements, with number of elements: ', NUMELEMS
               END IF

               DO J = 1, NUMELEMS
                  READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
               END DO
            END DO

         END IF
      END DO

      REWIND(in5)

      ALLOCATE(N_CELLS_WITH_NODE(U2D_GRID%NUM_NODES))
      ALLOCATE(IOF(U2D_GRID%NUM_NODES))

      N_CELLS_WITH_NODE = 0
      DO I = 1, U2D_GRID%NUM_CELLS
         DO V1 = 1, 3
            JN = U2D_GRID%CELL_NODES(V1,I)
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO
   
      IOF = -1
      IDX = 1
      DO JN = 1, U2D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) .NE. 0) THEN
            IOF(JN) = IDX
            IDX = IDX + N_CELLS_WITH_NODE(JN)
         END IF
      END DO
   
      ALLOCATE(CELL_WITH_NODE(IDX))
      
      N_CELLS_WITH_NODE = 0
      DO I = 1, U2D_GRID%NUM_CELLS
         DO V1 = 1, 3
            JN = U2D_GRID%CELL_NODES(V1,I)
            CELL_WITH_NODE(IOF(JN) + N_CELLS_WITH_NODE(JN)) = I
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO


      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE, NUM
         IF (ReasonEOF < 0) EXIT 
         !WRITE(*,*) 'Read line:', LINE, ' number ', NUM
         
         IF (LINE == 'NPOIN=') THEN
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
            END DO
         ELSE IF (LINE == 'NELEM=') THEN
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
            END DO
         ELSE IF (LINE == 'NMARK=') THEN

            ! Assign physical groups to cells/cell edges.

            ALLOCATE(GRID_BC(NUM)) ! Append the physical group to the list
            N_GRID_BC = NUM

            DO I = 1, NUM
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, GROUPNAME
      
               GRID_BC(I)%PHYSICAL_GROUP_NAME = GROUPNAME
         
               READ(in5,*, IOSTAT=ReasonEOF) LINE, NUMELEMS

               DO J = 1, NUMELEMS
                  READ(in5,'(A)', IOSTAT=ReasonEOF) LINE

                  READ(LINE,*) ELEM_TYPE

                  IF (ELEM_TYPE == 3) THEN ! element in physical group is a cell edge (segment).

                     READ(LINE,*) ELEM_TYPE, VLIST2
                     
                     VLIST2 = VLIST2 + 1

                     JN = VLIST2(1)
                     IF (N_CELLS_WITH_NODE(JN) > 0) THEN
                        DO IDX = 0, N_CELLS_WITH_NODE(JN) - 1
                           JC1 = CELL_WITH_NODE(IOF(JN) + IDX)
                           FOUND = 0
                           DO V1 = 1, 3
                              IF (ANY(VLIST2 == U2D_GRID%CELL_NODES(V1,JC1))) THEN
                                 FOUND = FOUND + 1
                                 WHICH1(FOUND) = V1
                              END IF
                           END DO
            
                           IF (FOUND == 2) THEN
                              IF (ANY(WHICH1 == 1) .AND. ANY(WHICH1 == 2)) THEN
                                 U2D_GRID%CELL_EDGES_PG(1, JC1) = I
                              ELSE IF (ANY(WHICH1 == 2) .AND. ANY(WHICH1 == 3)) THEN
                                 U2D_GRID%CELL_EDGES_PG(2, JC1) = I
                              ELSE IF (ANY(WHICH1 == 3) .AND. ANY(WHICH1 == 1)) THEN
                                 U2D_GRID%CELL_EDGES_PG(3, JC1) = I
                              END IF
                           END IF
                        END DO
                     END IF

                  ELSE IF (ELEM_TYPE == 5) THEN ! element in physical group is a cell (simplex).
                     READ(LINE,*) ELEM_TYPE, VLIST3
                     
                     VLIST3 = VLIST3 + 1

                     JN = VLIST3(1)
                     IF (N_CELLS_WITH_NODE(JN) > 0) THEN
                        DO IDX = 0, N_CELLS_WITH_NODE(JN) - 1
                           JC1 = CELL_WITH_NODE(IOF(JN) + IDX)
                           FOUND = 0
                           DO V1 = 1, 3
                              IF (ANY(VLIST3 == U2D_GRID%CELL_NODES(V1,JC1))) THEN
                                 FOUND = FOUND + 1
                              END IF
                           END DO
            
                           IF (FOUND == 3) THEN
                              U2D_GRID%CELL_PG(JC1) = I
                           END IF
                        END DO
                     END IF

                  ELSE
                     WRITE(*,*) 'Error! element type was not line or triangle.'
                  END IF

               END DO
            END DO

         END IF
      END DO

      ! Done reading
      CLOSE(in5)

      !WRITE(*,*) 'Read grid file. It contains ', U2D_GRID%NUM_NODES, &
      !           'points, and ', U2D_GRID%NUM_CELLS, 'cells.'

      ! Process the mesh: generate connectivity, normals and such...
      !XMIN, XMAX,...

      !DO I = 1, U2D_GRID%NUM_CELLS
      !   WRITE(*,*) U2D_GRID%CELL_NODES(:,I)
      !END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing cell volumes.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute cell volumes
      ALLOCATE(U2D_GRID%CELL_AREAS(U2D_GRID%NUM_CELLS))
      ALLOCATE(U2D_GRID%CELL_VOLUMES(U2D_GRID%NUM_CELLS))
      DO I = 1, U2D_GRID%NUM_CELLS
         A = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(1,I))
         B = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(2,I))
         C = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(3,I))

         U2D_GRID%CELL_AREAS(I) = 0.5*ABS(A(1)*(B(2)-C(2)) + B(1)*(C(2)-A(2)) + C(1)*(A(2)-B(2)))
         IF (DIMS == 2 .AND. .NOT. AXI) THEN
            U2D_GRID%CELL_VOLUMES(I) = U2D_GRID%CELL_AREAS(I) * (ZMAX-ZMIN)
            !WRITE(*,*) U2D_GRID%CELL_VOLUMES(I)
         END IF
         IF (DIMS == 2 .AND. AXI) THEN
            RAD = (A(2)+B(2)+C(2))/3.
            U2D_GRID%CELL_VOLUMES(I) = U2D_GRID%CELL_AREAS(I) * (ZMAX-ZMIN)*RAD
         END IF
      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing grid connectivity.'
         WRITE(*,*) '==========================================='
      END IF

      ! Find cell connectivity
      ALLOCATE(TEMP_CELL_NEIGHBORS(3, U2D_GRID%NUM_CELLS))
      TEMP_CELL_NEIGHBORS = -1



      DO JN = 1, U2D_GRID%NUM_NODES
         !IF (PROC_ID == 0) WRITE(*,*) 'Checking node ', JN, ' of ',  U2D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) > 1) THEN
            DO I = 0, N_CELLS_WITH_NODE(JN) - 1
               DO J = I, N_CELLS_WITH_NODE(JN) - 1
                  IF (I == J) CYCLE
                  JC1 = CELL_WITH_NODE(IOF(JN) + I)
                  JC2 = CELL_WITH_NODE(IOF(JN) + J)


                  FOUND = 0
                  DO V1 = 1, 3
                     DO V2 = 1, 3
                        IF (U2D_GRID%CELL_NODES(V1,JC1) == U2D_GRID%CELL_NODES(V2,JC2)) THEN
                           FOUND = FOUND + 1
                           IF (FOUND .GT. 2) CALL ERROR_ABORT('Error! Found duplicate cells in the mesh!')
                           WHICH1(FOUND) = V1
                           WHICH2(FOUND) = V2
                        END IF
                     END DO
                  END DO

                  IF (FOUND == 2) THEN
      
                     IF (ANY(WHICH1 == 1) .AND. ANY(WHICH1 == 2)) THEN
                        TEMP_CELL_NEIGHBORS(1, JC1) = JC2
                     ELSE IF (ANY(WHICH1 == 2) .AND. ANY(WHICH1 == 3)) THEN
                        TEMP_CELL_NEIGHBORS(2, JC1) = JC2
                     ELSE IF (ANY(WHICH1 == 3) .AND. ANY(WHICH1 == 1)) THEN
                        TEMP_CELL_NEIGHBORS(3, JC1) = JC2
                     END IF

                     IF (ANY(WHICH2 == 1) .AND. ANY(WHICH2 == 2)) THEN
                        TEMP_CELL_NEIGHBORS(1, JC2) = JC1
                     ELSE IF (ANY(WHICH2 == 2) .AND. ANY(WHICH2 == 3)) THEN
                        TEMP_CELL_NEIGHBORS(2, JC2) = JC1
                     ELSE IF (ANY(WHICH2 == 3) .AND. ANY(WHICH2 == 1)) THEN
                        TEMP_CELL_NEIGHBORS(3, JC2) = JC1
                     END IF
      
                  END IF


               END DO
            END DO
         END IF
      END DO

      U2D_GRID%CELL_NEIGHBORS = TEMP_CELL_NEIGHBORS



      !WRITE(*,*) 'Generated grid connectivity. '
      !DO I = 1, U2D_GRID%NUM_CELLS
      !   WRITE(*,*) 'Cell ', I, ' neighbors cells ', TEMP_CELL_NEIGHBORS(:, I)
      !END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing face normals.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute cell edge normals
      IND(1,:) = [1,2]
      IND(2,:) = [2,3]
      IND(3,:) = [3,1]
      ALLOCATE(U2D_GRID%EDGE_NORMAL(3, 3, U2D_GRID%NUM_CELLS))
      ALLOCATE(U2D_GRID%CELL_EDGES_LEN(3, U2D_GRID%NUM_CELLS))
      DO I = 1, U2D_GRID%NUM_CELLS
         DO J = 1, 3
            X1 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(IND(J,1),I))
            X2 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(IND(J,2),I))
            Y1 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(IND(J,1),I))
            Y2 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(IND(J,2),I))
            LEN = SQRT((Y2-Y1)*(Y2-Y1) + (X2-X1)*(X2-X1))
            U2D_GRID%CELL_EDGES_LEN(J,I) = LEN
            U2D_GRID%EDGE_NORMAL(1,J,I) = (Y2-Y1)/LEN
            U2D_GRID%EDGE_NORMAL(2,J,I) = (X1-X2)/LEN
            U2D_GRID%EDGE_NORMAL(3,J,I) = 0.d0
         END DO
      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Checking ordering.'
         WRITE(*,*) '==========================================='
      END IF

      DO I = 1, U2D_GRID%NUM_CELLS
         X1 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(2,I)) &
            - U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(1,I))
         X2 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(3,I)) &
            - U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(1,I))
         Y1 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(2,I)) &
            - U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(1,I))
         Y2 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(3,I)) &
            - U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(1,I))

         IF (X1*Y2-X2*Y1 < 0) CALL ERROR_ABORT('2D mesh triangles have negative z-normal.')

      END DO

      NCELLS = U2D_GRID%NUM_CELLS
      NNODES = U2D_GRID%NUM_NODES



      ALLOCATE(U2D_GRID%BASIS_COEFFS(3,3,NCELLS))

      DO I = 1, NCELLS
         V1 = U2D_GRID%CELL_NODES(1,I)
         V2 = U2D_GRID%CELL_NODES(2,I)
         V3 = U2D_GRID%CELL_NODES(3,I)

         X1 = U2D_GRID%NODE_COORDS(1, V1)
         X2 = U2D_GRID%NODE_COORDS(1, V2)
         X3 = U2D_GRID%NODE_COORDS(1, V3)
         Y1 = U2D_GRID%NODE_COORDS(2, V1)
         Y2 = U2D_GRID%NODE_COORDS(2, V2)
         Y3 = U2D_GRID%NODE_COORDS(2, V3)

         ! These are such that PSI_i = SUM_j [ x_j * BASIS_COEFFS(j,i,IC) ] + BASIS_COEFFS(3,i,IC)

         U2D_GRID%BASIS_COEFFS(1,1,I) =  Y2-Y3
         U2D_GRID%BASIS_COEFFS(2,1,I) = -(X2-X3)
         U2D_GRID%BASIS_COEFFS(3,1,I) =  X2*Y3 - X3*Y2

         U2D_GRID%BASIS_COEFFS(1,2,I) = -(Y1-Y3)
         U2D_GRID%BASIS_COEFFS(2,2,I) =  X1-X3
         U2D_GRID%BASIS_COEFFS(3,2,I) =  X3*Y1 - X1*Y3

         U2D_GRID%BASIS_COEFFS(1,3,I) = -(Y2-Y1)
         U2D_GRID%BASIS_COEFFS(2,3,I) =  X2-X1
         U2D_GRID%BASIS_COEFFS(3,3,I) =  X1*Y2 - X2*Y1

         U2D_GRID%BASIS_COEFFS(:,:,I) = 0.5*U2D_GRID%BASIS_COEFFS(:,:,I)/U2D_GRID%CELL_AREAS(I)

      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Creating boundary grid.'
         WRITE(*,*) '==========================================='
      END IF

      ALLOCATE(U2D_GRID%CELL_EDGES_BOUNDARY_INDEX(3,NCELLS))
      U2D_GRID%CELL_EDGES_BOUNDARY_INDEX = -1
      ALLOCATE(NODE_ON_BOUNDARY(NNODES))
      NODE_ON_BOUNDARY = .FALSE.
      ALLOCATE(U2D_GRID%NODES_BOUNDARY_INDEX(NNODES))
      U2D_GRID%NODES_BOUNDARY_INDEX = -1
      NUM_BOUNDARY_NODES = 0
      NUM_BOUNDARY_ELEM = 0
      DO I = 1, NCELLS
         DO J = 1, 3
            ! If the edge belongs to any physical group, it should be part of the boundary grid
            ! Later, we may want to filter this further
            IF (U2D_GRID%CELL_EDGES_PG(J,I) .NE. -1) THEN
               NUM_BOUNDARY_ELEM = NUM_BOUNDARY_ELEM + 1
               V1 = U2D_GRID%CELL_NODES(J, I)
               IF (J == 3) THEN
                  V2 = U2D_GRID%CELL_NODES(1, I)
               ELSE
                  V2 = U2D_GRID%CELL_NODES(J+1, I)
               END IF
               IF (.NOT. NODE_ON_BOUNDARY(V1)) THEN
                  NUM_BOUNDARY_NODES = NUM_BOUNDARY_NODES + 1
                  U2D_GRID%NODES_BOUNDARY_INDEX(V1) = NUM_BOUNDARY_NODES
                  NODE_ON_BOUNDARY(V1) = .TRUE.
               END IF
               IF (.NOT. NODE_ON_BOUNDARY(V2)) THEN
                  NUM_BOUNDARY_NODES = NUM_BOUNDARY_NODES + 1
                  U2D_GRID%NODES_BOUNDARY_INDEX(V2) = NUM_BOUNDARY_NODES
                  NODE_ON_BOUNDARY(V2) = .TRUE.
               END IF
               
            END IF
         END DO
      END DO      

      U1D_GRID%NUM_CELLS = NUM_BOUNDARY_ELEM
      U1D_GRID%NUM_NODES = NUM_BOUNDARY_NODES
      ALLOCATE(U1D_GRID%CELL_NODES(2, NUM_BOUNDARY_ELEM))
      ALLOCATE(U1D_GRID%CELL_PG(NUM_BOUNDARY_ELEM))
      ALLOCATE(U1D_GRID%NODE_COORDS(3, NUM_BOUNDARY_NODES))

      DO I = 1, NNODES
         IF (NODE_ON_BOUNDARY(I)) THEN
            U1D_GRID%NODE_COORDS(:,U2D_GRID%NODES_BOUNDARY_INDEX(I)) = U2D_GRID%NODE_COORDS(:,I)
         END IF
      END DO

      NUM_BOUNDARY_ELEM = 0

      DO I = 1, NCELLS
         DO J = 1, 3
            IF (U2D_GRID%CELL_EDGES_PG(J,I) .NE. -1) THEN
               NUM_BOUNDARY_ELEM = NUM_BOUNDARY_ELEM + 1
               U1D_GRID%CELL_PG(NUM_BOUNDARY_ELEM) = U2D_GRID%CELL_EDGES_PG(J,I)
               U2D_GRID%CELL_EDGES_BOUNDARY_INDEX(J,I) = NUM_BOUNDARY_ELEM

               V1 = U2D_GRID%CELL_NODES(J, I)
               IF (J == 3) THEN
                  V2 = U2D_GRID%CELL_NODES(1, I)
               ELSE
                  V2 = U2D_GRID%CELL_NODES(J+1, I)
               END IF
               U1D_GRID%CELL_NODES(1, NUM_BOUNDARY_ELEM) = U2D_GRID%NODES_BOUNDARY_INDEX(V1)
               U1D_GRID%CELL_NODES(2, NUM_BOUNDARY_ELEM) = U2D_GRID%NODES_BOUNDARY_INDEX(V2)

            END IF
         END DO
      END DO
      
      DEALLOCATE(NODE_ON_BOUNDARY)

      NBOUNDCELLS = NUM_BOUNDARY_ELEM
      NBOUNDNODES = NUM_BOUNDARY_NODES

      ! Compute areas and lengths of boundary mesh
      ALLOCATE(U1D_GRID%SEGMENT_LENGTHS(U1D_GRID%NUM_CELLS))
      ALLOCATE(U1D_GRID%SEGMENT_AREAS(U1D_GRID%NUM_CELLS))
      DO I = 1, U1D_GRID%NUM_CELLS
         A = U1D_GRID%NODE_COORDS(:, U1D_GRID%CELL_NODES(1,I))
         B = U1D_GRID%NODE_COORDS(:, U1D_GRID%CELL_NODES(2,I))

         U1D_GRID%SEGMENT_LENGTHS(I) = SQRT((B(1) - A(1))**2 + (B(2) - A(2))**2)
         IF (.NOT. AXI) THEN
            U1D_GRID%SEGMENT_AREAS(I) = U1D_GRID%SEGMENT_LENGTHS(I) * (ZMAX-ZMIN)
         ELSE
            RAD = 0.5*(A(2)+B(2))
            U1D_GRID%SEGMENT_AREAS(I) = U1D_GRID%SEGMENT_LENGTHS(I) * (ZMAX-ZMIN)*RAD
         END IF
      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '============================================================='
         WRITE(*,*) 'Done reading grid file.'
         WRITE(*,*) 'It contains ', NNODES, ' nodes and ', NCELLS, ' cells.'
         WRITE(*,*) 'The boundary grid contains ', NBOUNDCELLS, ' lines and ', NBOUNDNODES, ' nodes.'
         WRITE(*,*) '============================================================='
      END IF

   END SUBROUTINE READ_2D_UNSTRUCTURED_GRID_SU2








   SUBROUTINE READ_3D_UNSTRUCTURED_GRID_SU2(FILENAME)

      IMPLICIT NONE

      CHARACTER*256, INTENT(IN) :: FILENAME

      CHARACTER*256 :: LINE, GROUPNAME, DUMMYLINE

      INTEGER, PARAMETER :: in5 = 2385
      INTEGER            :: ios
      INTEGER            :: ReasonEOF

      INTEGER            :: NUM, I, J, FOUND, V1, V2, V3, V4, ELEM_TYPE, NUMELEMS
      INTEGER, DIMENSION(4,3) :: IND
      REAL(KIND=8), DIMENSION(3) :: XYZ, A, B, C, CROSSP

      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: TEMP_CELL_NEIGHBORS

      INTEGER, DIMENSION(3) :: VLIST3, WHICH1, WHICH2
      INTEGER, DIMENSION(4) :: VLIST4

      INTEGER, DIMENSION(:), ALLOCATABLE :: N_CELLS_WITH_NODE, CELL_WITH_NODE, IOF
      INTEGER :: IDX, JN, JC1, JC2

      REAL(KIND=8) :: VOLUME, X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4

      LOGICAL, DIMENSION(:), ALLOCATABLE :: NODE_ON_BOUNDARY
      INTEGER :: NUM_BOUNDARY_NODES, NUM_BOUNDARY_ELEM

      ! Open input file for reading
      OPEN(UNIT=in5,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios .NE. 0) THEN
         CALL ERROR_ABORT('Attention, mesh file not found! ABORTING.')
      ENDIF

      ! Read the mesh file. SU2 file format (*.su2)
      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Reading grid file in SU2 format.'
         WRITE(*,*) '==========================================='
      END IF
      
      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE, NUM
         IF (ReasonEOF < 0) EXIT 
         !WRITE(*,*) 'Read line:', LINE, ' number ', NUM
         
         IF (LINE == 'NPOIN=') THEN
            ALLOCATE(U3D_GRID%NODE_COORDS(3,NUM))
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) XYZ
               U3D_GRID%NODE_COORDS(:,I) = XYZ
            END DO
            U3D_GRID%NUM_NODES = NUM
         ELSE IF (LINE == 'NELEM=') THEN
            ALLOCATE(U3D_GRID%CELL_NODES(4,NUM))

            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) ELEM_TYPE, U3D_GRID%CELL_NODES(:,I)
               !WRITE(*,*) 'I read element ', I, ' has nodes ', U3D_GRID%CELL_NODES(:,I)
               IF (ELEM_TYPE .NE. 10) CALL ERROR_ABORT('Reading 3D grid found element type was not tetrahedron (type 10).')
            END DO
            U3D_GRID%CELL_NODES = U3D_GRID%CELL_NODES + 1 ! Start indexing from 1.

            U3D_GRID%NUM_CELLS = NUM

            ALLOCATE(U3D_GRID%CELL_FACES_PG(4, U3D_GRID%NUM_CELLS))
            U3D_GRID%CELL_FACES_PG = -1
            ALLOCATE(U3D_GRID%CELL_PG(U3D_GRID%NUM_CELLS))
            U3D_GRID%CELL_PG = -1

         ELSE IF (LINE == 'NMARK=') THEN
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) LINE, GROUPNAME
               IF (LINE .NE. 'MARKER_TAG=') THEN
                  CALL ERROR_ABORT('Error! did not find marker name.')
               ELSE
                  !IF (PROC_ID == 0) WRITE(*,*) 'Found marker tag, with groupname: ', GROUPNAME
               END IF

               READ(in5,*, IOSTAT=ReasonEOF) LINE, NUMELEMS
               IF (LINE .NE. 'MARKER_ELEMS=') THEN
                  CALL ERROR_ABORT('Error! did not find marker elements.')
               ELSE
                  !IF (PROC_ID == 0) WRITE(*,*) 'Found marker elements, with number of elements: ', NUMELEMS
               END IF

               DO J = 1, NUMELEMS
                  READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
               END DO
            END DO
         END IF
      END DO

      REWIND(in5)



      ALLOCATE(N_CELLS_WITH_NODE(U3D_GRID%NUM_NODES))
      ALLOCATE(IOF(U3D_GRID%NUM_NODES))

      N_CELLS_WITH_NODE = 0
      DO I = 1, U3D_GRID%NUM_CELLS
         DO V1 = 1, 4
            JN = U3D_GRID%CELL_NODES(V1,I)
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO
   
      IOF = -1
      IDX = 1
      DO JN = 1, U3D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) .NE. 0) THEN
            IOF(JN) = IDX
            IDX = IDX + N_CELLS_WITH_NODE(JN)
         END IF
      END DO
   
      ALLOCATE(CELL_WITH_NODE(IDX))
      
      N_CELLS_WITH_NODE = 0
      DO I = 1, U3D_GRID%NUM_CELLS
         DO V1 = 1, 4
            JN = U3D_GRID%CELL_NODES(V1,I)
            CELL_WITH_NODE(IOF(JN) + N_CELLS_WITH_NODE(JN)) = I
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO


      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE, NUM
         IF (ReasonEOF < 0) EXIT 
         !WRITE(*,*) 'Read line:', LINE, ' number ', NUM
         
         IF (LINE == 'NPOIN=') THEN
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
            END DO
         ELSE IF (LINE == 'NELEM=') THEN
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) DUMMYLINE
            END DO
         ELSE IF (LINE == 'NMARK=') THEN

            ! Assign physical groups to cell edges.

            ALLOCATE(GRID_BC(NUM)) ! Append the physical group to the list
            N_GRID_BC = NUM

            DO I = 1, NUM
               
               READ(in5,*, IOSTAT=ReasonEOF) LINE, GROUPNAME
      
               GRID_BC(I)%PHYSICAL_GROUP_NAME = GROUPNAME
         
               READ(in5,*, IOSTAT=ReasonEOF) LINE, NUMELEMS

               DO J = 1, NUMELEMS
                  READ(in5,'(A)', IOSTAT=ReasonEOF) LINE

                  READ(LINE,*) ELEM_TYPE

                  IF (ELEM_TYPE == 5) THEN ! element in physical group is a cell cell (simplex).

                     READ(LINE,*) ELEM_TYPE, VLIST3

                     VLIST3 = VLIST3 + 1

                     JN = VLIST3(1)
                     IF (N_CELLS_WITH_NODE(JN) > 0) THEN
                        DO IDX = 0, N_CELLS_WITH_NODE(JN) - 1
                           JC1 = CELL_WITH_NODE(IOF(JN) + IDX)
                           FOUND = 0
                           DO V1 = 1, 4
                              IF (ANY(VLIST3 == U3D_GRID%CELL_NODES(V1,JC1))) THEN
                                 FOUND = FOUND + 1
                                 WHICH1(FOUND) = V1
                              END IF
                           END DO
            
                           IF (FOUND == 3) THEN
               
                              IF (ANY(WHICH1 == 1)) THEN
                                 IF (ANY(WHICH1 == 2)) THEN
                                    IF (ANY(WHICH1 == 3))  THEN
                                       U3D_GRID%CELL_FACES_PG(1, JC1) = I
                                    ELSE IF (ANY(WHICH1 == 4)) THEN
                                       U3D_GRID%CELL_FACES_PG(2, JC1) = I
                                    END IF
                                 ELSE IF (ANY(WHICH1 == 3)) THEN
                                    IF (ANY(WHICH1 == 4)) U3D_GRID%CELL_FACES_PG(4, JC1) = I
                                 END IF
                              ELSE IF (ANY(WHICH1 == 2)) THEN
                                 IF (ANY(WHICH1 == 3) .AND. ANY(WHICH1 == 4)) U3D_GRID%CELL_FACES_PG(3, JC1) = I
                              END IF

                           END IF
                        END DO
                     END IF

                  ELSE IF (ELEM_TYPE == 10) THEN ! element in physical group is a tetrahedron.
                     READ(LINE,*) ELEM_TYPE, VLIST4
                     
                     VLIST4 = VLIST4 + 1

                     JN = VLIST4(1)
                     IF (N_CELLS_WITH_NODE(JN) > 0) THEN
                        DO IDX = 0, N_CELLS_WITH_NODE(JN) - 1
                           JC1 = CELL_WITH_NODE(IOF(JN) + IDX)
                           FOUND = 0
                           DO V1 = 1, 4
                              IF (ANY(VLIST4 == U3D_GRID%CELL_NODES(V1,JC1))) THEN
                                 FOUND = FOUND + 1
                              END IF
                           END DO
            
                           IF (FOUND == 4) THEN
                              U3D_GRID%CELL_PG(JC1) = I
                           END IF
                        END DO
                     END IF

                  ELSE
                     WRITE(*,*) 'Error! element type was not triangle or prism.'
                  END IF

               END DO
            END DO

         END IF
      END DO

      ! Done reading
      CLOSE(in5)

      !WRITE(*,*) 'Read grid file. It contains ', U3D_GRID%NUM_NODES, &
      !           'points, and ', U3D_GRID%NUM_CELLS, 'cells.'

      ! Process the mesh: generate connectivity, normals and such...
      !XMIN, XMAX,...

      !DO I = 1, U3D_GRID%NUM_CELLS
      !   WRITE(*,*) U3D_GRID%CELL_NODES(I,:)
      !END DO
      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing cell volumes.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute cell volumes
      ALLOCATE(U3D_GRID%CELL_VOLUMES(U3D_GRID%NUM_CELLS))
      DO I = 1, U3D_GRID%NUM_CELLS
         A = U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(2,I)) - U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(1,I))
         B = U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(3,I)) - U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(1,I))
         C = U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(4,I)) - U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(1,I))

         U3D_GRID%CELL_VOLUMES(I) = ABS(C(1)*(A(2)*B(3)-A(3)*B(2)) + C(2)*(A(3)*B(1)-A(1)*B(3)) + C(3)*(A(1)*B(2)-A(2)*B(1))) / 6.
      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing grid connectivity.'
         WRITE(*,*) '==========================================='
      END IF


      ALLOCATE(TEMP_CELL_NEIGHBORS(4, U3D_GRID%NUM_CELLS))
      TEMP_CELL_NEIGHBORS = -1

      DO JN = 1, U3D_GRID%NUM_NODES
         !IF (PROC_ID == 0) WRITE(*,*) 'Checking node ', JN, ' of ',  U3D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) > 1) THEN
            DO I = 0, N_CELLS_WITH_NODE(JN) - 1
               DO J = I, N_CELLS_WITH_NODE(JN) - 1
                  IF (I == J) CYCLE
                  JC1 = CELL_WITH_NODE(IOF(JN) + I)
                  JC2 = CELL_WITH_NODE(IOF(JN) + J)


                  FOUND = 0
                  DO V1 = 1, 4
                     DO V2 = 1, 4
                        IF (U3D_GRID%CELL_NODES(V1,JC1) == U3D_GRID%CELL_NODES(V2,JC2)) THEN
                           FOUND = FOUND + 1
                           IF (FOUND .GT. 3) CALL ERROR_ABORT('Error! Found duplicate cells in the mesh!')
                           WHICH1(FOUND) = V1
                           WHICH2(FOUND) = V2
                        END IF
                     END DO
                  END DO

                  IF (FOUND == 3) THEN
      
                     IF (ANY(WHICH1 == 1)) THEN
                        IF (ANY(WHICH1 == 2)) THEN
                           IF (ANY(WHICH1 == 3))  THEN
                              TEMP_CELL_NEIGHBORS(1, JC1) = JC2
                           ELSE IF (ANY(WHICH1 == 4)) THEN
                              TEMP_CELL_NEIGHBORS(2, JC1) = JC2
                           END IF
                        ELSE IF (ANY(WHICH1 == 3)) THEN
                           IF (ANY(WHICH1 == 4)) TEMP_CELL_NEIGHBORS(4, JC1) = JC2
                        END IF
                     ELSE IF (ANY(WHICH1 == 2)) THEN
                        IF (ANY(WHICH1 == 3) .AND. ANY(WHICH1 == 4)) TEMP_CELL_NEIGHBORS(3, JC1) = JC2
                     END IF
      
      
                     IF (ANY(WHICH2 == 1)) THEN
                        IF (ANY(WHICH2 == 2)) THEN
                           IF (ANY(WHICH2 == 3))  THEN
                              TEMP_CELL_NEIGHBORS(1, JC2) = JC1
                           ELSE IF (ANY(WHICH2 == 4)) THEN
                              TEMP_CELL_NEIGHBORS(2, JC2) = JC1
                           END IF
                        ELSE IF (ANY(WHICH2 == 3)) THEN
                           IF (ANY(WHICH2 == 4)) TEMP_CELL_NEIGHBORS(4, JC2) = JC1
                        END IF
                     ELSE IF (ANY(WHICH2 == 2)) THEN
                        IF (ANY(WHICH2 == 3) .AND. ANY(WHICH2 == 4)) TEMP_CELL_NEIGHBORS(3, JC2) = JC1
                     END IF
      
                  END IF


               END DO
            END DO
         END IF
      END DO

      U3D_GRID%CELL_NEIGHBORS = TEMP_CELL_NEIGHBORS

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing face normals.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute cell edge normals
      IND(1,:) = [1,3,2]
      IND(2,:) = [1,2,4]
      IND(3,:) = [2,3,4]
      IND(4,:) = [1,4,3]
      ALLOCATE(U3D_GRID%FACE_NORMAL(3, 4, U3D_GRID%NUM_CELLS))
      ALLOCATE(U3D_GRID%FACE_TANG1(3, 4, U3D_GRID%NUM_CELLS))
      ALLOCATE(U3D_GRID%FACE_TANG2(3, 4, U3D_GRID%NUM_CELLS))
      ALLOCATE(U3D_GRID%FACE_NODES(3, 4, U3D_GRID%NUM_CELLS))
      ALLOCATE(U3D_GRID%CELL_FACES_COEFFS(4, 4, U3D_GRID%NUM_CELLS))
      ALLOCATE(U3D_GRID%FACE_AREA(4, U3D_GRID%NUM_CELLS))

      DO I = 1, U3D_GRID%NUM_CELLS
         DO J = 1, 4
            V1 = U3D_GRID%CELL_NODES(IND(J,1),I)
            V2 = U3D_GRID%CELL_NODES(IND(J,2),I)
            V3 = U3D_GRID%CELL_NODES(IND(J,3),I)

            U3D_GRID%FACE_NODES(1,J,I) = V1
            U3D_GRID%FACE_NODES(2,J,I) = V2
            U3D_GRID%FACE_NODES(3,J,I) = V3

            A = U3D_GRID%NODE_COORDS(:,V1)
            B = U3D_GRID%NODE_COORDS(:,V2)
            C = U3D_GRID%NODE_COORDS(:,V3)
            
            CROSSP = CROSS(B-A,C-A)
            U3D_GRID%FACE_NORMAL(:,J,I) = CROSSP/NORM2(CROSSP)
            U3D_GRID%FACE_AREA(J,I) = 0.5*NORM2(CROSSP)
            
            U3D_GRID%FACE_TANG1(:,J,I) = (B-A)/NORM2(B-A)
            U3D_GRID%FACE_TANG2(:,J,I) = CROSS(U3D_GRID%FACE_NORMAL(:,J,I), U3D_GRID%FACE_TANG1(:,J,I))


            ! The coefficients (a,b,c,d) of a*x + b*y + c*z + d = 0
            U3D_GRID%CELL_FACES_COEFFS(1,J,I) =  A(2)*B(3)-B(2)*A(3) &
                                                +B(2)*C(3)-C(2)*B(3) &
                                                +C(2)*A(3)-A(2)*C(3)
            U3D_GRID%CELL_FACES_COEFFS(2,J,I) = -A(1)*B(3)+B(1)*A(3) &
                                                -B(1)*C(3)+C(1)*B(3) &
                                                -C(1)*A(3)+A(1)*C(3)
            U3D_GRID%CELL_FACES_COEFFS(3,J,I) =  A(1)*B(2)-B(1)*A(2) &
                                                +B(1)*C(2)-C(1)*B(2) &
                                                +C(1)*A(2)-A(1)*C(2)
            U3D_GRID%CELL_FACES_COEFFS(4,J,I) = -A(1)*B(2)*C(3) &
                                                +A(1)*C(2)*B(3) &
                                                +B(1)*A(2)*C(3) &
                                                -B(1)*C(2)*A(3) &
                                                -C(1)*A(2)*B(3) &
                                                +C(1)*B(2)*A(3)

         END DO
      END DO



      NCELLS = U3D_GRID%NUM_CELLS
      NNODES = U3D_GRID%NUM_NODES



      ALLOCATE(U3D_GRID%BASIS_COEFFS(4,4,NCELLS))

      DO I = 1, NCELLS
         VOLUME = U3D_GRID%CELL_VOLUMES(I)
         V1 = U3D_GRID%CELL_NODES(1,I)
         V2 = U3D_GRID%CELL_NODES(2,I)
         V3 = U3D_GRID%CELL_NODES(3,I)
         V4 = U3D_GRID%CELL_NODES(4,I)

         X1 = U3D_GRID%NODE_COORDS(1, V1)
         X2 = U3D_GRID%NODE_COORDS(1, V2)
         X3 = U3D_GRID%NODE_COORDS(1, V3)
         X4 = U3D_GRID%NODE_COORDS(1, V4)
         Y1 = U3D_GRID%NODE_COORDS(2, V1)
         Y2 = U3D_GRID%NODE_COORDS(2, V2)
         Y3 = U3D_GRID%NODE_COORDS(2, V3)
         Y4 = U3D_GRID%NODE_COORDS(2, V4)
         Z1 = U3D_GRID%NODE_COORDS(3, V1)
         Z2 = U3D_GRID%NODE_COORDS(3, V2)
         Z3 = U3D_GRID%NODE_COORDS(3, V3)
         Z4 = U3D_GRID%NODE_COORDS(3, V4)


         ! These are such that PSI_i = SUM_j [ x_j * BASIS_COEFFS(j,i,IC) ] + BASIS_COEFFS(4,i,IC)

         U3D_GRID%BASIS_COEFFS(1,1,I) =  Y2*Z3-Y3*Z2 -Y2*Z4+Y4*Z2 +Y3*Z4-Y4*Z3
         U3D_GRID%BASIS_COEFFS(2,1,I) = -X2*Z3+X3*Z2 +X2*Z4-X4*Z2 -X3*Z4+X4*Z3
         U3D_GRID%BASIS_COEFFS(3,1,I) =  X2*Y3-X3*Y2 -X2*Y4+X4*Y2 +X3*Y4-X4*Y3
         U3D_GRID%BASIS_COEFFS(4,1,I) = -X2*Y3*Z4 +X3*Y2*Z4 +X2*Y4*Z3 -X4*Y2*Z3 -X3*Y4*Z2 +X4*Y3*Z2

         U3D_GRID%BASIS_COEFFS(1,2,I) = -Y1*Z3+Y3*Z1 +Y1*Z4-Y4*Z1 -Y3*Z4+Y4*Z3
         U3D_GRID%BASIS_COEFFS(2,2,I) =  X1*Z3-X3*Z1 -X1*Z4+X4*Z1 +X3*Z4-X4*Z3
         U3D_GRID%BASIS_COEFFS(3,2,I) = -X1*Y3+X3*Y1 +X1*Y4-X4*Y1 -X3*Y4+X4*Y3
         U3D_GRID%BASIS_COEFFS(4,2,I) =  X1*Y3*Z4 -X3*Y1*Z4 -X1*Y4*Z3 +X4*Y1*Z3 +X3*Y4*Z1 -X4*Y3*Z1

         U3D_GRID%BASIS_COEFFS(1,3,I) =  Y1*Z2-Y2*Z1 -Y1*Z4+Y4*Z1 +Y2*Z4-Y4*Z2
         U3D_GRID%BASIS_COEFFS(2,3,I) = -X1*Z2+X2*Z1 +X1*Z4-X4*Z1 -X2*Z4+X4*Z2
         U3D_GRID%BASIS_COEFFS(3,3,I) =  X1*Y2-X2*Y1 -X1*Y4+X4*Y1 +X2*Y4-X4*Y2
         U3D_GRID%BASIS_COEFFS(4,3,I) = -X1*Y2*Z4 +X2*Y1*Z4 +X1*Y4*Z2 -X4*Y1*Z2 -X2*Y4*Z1 +X4*Y2*Z1

         U3D_GRID%BASIS_COEFFS(1,4,I) = -Y1*Z2+Y2*Z1 +Y1*Z3-Y3*Z1 -Y2*Z3+Y3*Z2
         U3D_GRID%BASIS_COEFFS(2,4,I) =  X1*Z2-X2*Z1 -X1*Z3+X3*Z1 +X2*Z3-X3*Z2
         U3D_GRID%BASIS_COEFFS(3,4,I) = -X1*Y2+X2*Y1 +X1*Y3-X3*Y1 -X2*Y3+X3*Y2
         U3D_GRID%BASIS_COEFFS(4,4,I) =  X1*Y2*Z3 -X2*Y1*Z3 -X1*Y3*Z2 +X3*Y1*Z2 +X2*Y3*Z1 -X3*Y2*Z1

         U3D_GRID%BASIS_COEFFS(:,:,I) = -U3D_GRID%BASIS_COEFFS(:,:,I)/6./VOLUME

      END DO





      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Creating boundary grid.'
         WRITE(*,*) '==========================================='
      END IF

      ALLOCATE(U3D_GRID%CELL_FACES_BOUNDARY_INDEX(4,NCELLS))
      U3D_GRID%CELL_FACES_BOUNDARY_INDEX = -1
      ALLOCATE(NODE_ON_BOUNDARY(NNODES))
      NODE_ON_BOUNDARY = .FALSE.
      ALLOCATE(U3D_GRID%NODES_BOUNDARY_INDEX(NNODES))
      U3D_GRID%NODES_BOUNDARY_INDEX = -1
      NUM_BOUNDARY_NODES = 0
      NUM_BOUNDARY_ELEM = 0
      DO I = 1, NCELLS
         DO J = 1, 4
            ! If the face belongs to any physical group, it should be part of the boundary grid
            ! Later, we may want to filter this further
            IF (U3D_GRID%CELL_FACES_PG(J,I) .NE. -1) THEN
               NUM_BOUNDARY_ELEM = NUM_BOUNDARY_ELEM + 1

               V1 = U3D_GRID%CELL_NODES(IND(J,1),I)
               V2 = U3D_GRID%CELL_NODES(IND(J,2),I)
               V3 = U3D_GRID%CELL_NODES(IND(J,3),I)

               IF (.NOT. NODE_ON_BOUNDARY(V1)) THEN
                  NUM_BOUNDARY_NODES = NUM_BOUNDARY_NODES + 1
                  U3D_GRID%NODES_BOUNDARY_INDEX(V1) = NUM_BOUNDARY_NODES
                  NODE_ON_BOUNDARY(V1) = .TRUE.
               END IF
               IF (.NOT. NODE_ON_BOUNDARY(V2)) THEN
                  NUM_BOUNDARY_NODES = NUM_BOUNDARY_NODES + 1
                  U3D_GRID%NODES_BOUNDARY_INDEX(V2) = NUM_BOUNDARY_NODES
                  NODE_ON_BOUNDARY(V2) = .TRUE.
               END IF
               IF (.NOT. NODE_ON_BOUNDARY(V3)) THEN
                  NUM_BOUNDARY_NODES = NUM_BOUNDARY_NODES + 1
                  U3D_GRID%NODES_BOUNDARY_INDEX(V3) = NUM_BOUNDARY_NODES
                  NODE_ON_BOUNDARY(V3) = .TRUE.
               END IF
               
            END IF
         END DO
      END DO      

      U2D_GRID%NUM_CELLS = NUM_BOUNDARY_ELEM
      U2D_GRID%NUM_NODES = NUM_BOUNDARY_NODES
      ALLOCATE(U2D_GRID%CELL_NODES(3, NUM_BOUNDARY_ELEM))
      ALLOCATE(U2D_GRID%CELL_PG(NUM_BOUNDARY_ELEM))
      ALLOCATE(U2D_GRID%NODE_COORDS(3, NUM_BOUNDARY_NODES))

      DO I = 1, NNODES
         IF (NODE_ON_BOUNDARY(I)) THEN
            U2D_GRID%NODE_COORDS(:,U3D_GRID%NODES_BOUNDARY_INDEX(I)) = U3D_GRID%NODE_COORDS(:,I)
         END IF
      END DO

      NUM_BOUNDARY_ELEM = 0

      DO I = 1, NCELLS
         DO J = 1, 4
            IF (U3D_GRID%CELL_FACES_PG(J,I) .NE. -1) THEN
               NUM_BOUNDARY_ELEM = NUM_BOUNDARY_ELEM + 1
               U2D_GRID%CELL_PG(NUM_BOUNDARY_ELEM) = U3D_GRID%CELL_FACES_PG(J,I)
               U3D_GRID%CELL_FACES_BOUNDARY_INDEX(J,I) = NUM_BOUNDARY_ELEM

               V1 = U3D_GRID%CELL_NODES(IND(J,1),I)
               V2 = U3D_GRID%CELL_NODES(IND(J,2),I)
               V3 = U3D_GRID%CELL_NODES(IND(J,3),I)
               U2D_GRID%CELL_NODES(1, NUM_BOUNDARY_ELEM) = U3D_GRID%NODES_BOUNDARY_INDEX(V1)
               U2D_GRID%CELL_NODES(2, NUM_BOUNDARY_ELEM) = U3D_GRID%NODES_BOUNDARY_INDEX(V2)
               U2D_GRID%CELL_NODES(3, NUM_BOUNDARY_ELEM) = U3D_GRID%NODES_BOUNDARY_INDEX(V3)

            END IF
         END DO
      END DO
      
      DEALLOCATE(NODE_ON_BOUNDARY)

      NBOUNDCELLS = NUM_BOUNDARY_ELEM
      NBOUNDNODES = NUM_BOUNDARY_NODES

      ! Compute areas and lengths of boundary mesh
      ALLOCATE(U2D_GRID%CELL_AREAS(U2D_GRID%NUM_CELLS))
      DO I = 1, U2D_GRID%NUM_CELLS
         A = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(1,I))
         B = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(2,I))
         C = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(3,I))

         U2D_GRID%CELL_AREAS(I) = 0.5*NORM2(CROSS(B-A, C-A))
      END DO




      IF (PROC_ID == 0) THEN
         WRITE(*,*) '============================================================='
         WRITE(*,*) 'Done reading grid file.'
         WRITE(*,*) 'It contains ', NNODES, ' nodes and ', NCELLS, ' cells.'
         WRITE(*,*) 'The boundary grid contains ', NBOUNDCELLS, ' surfaces and ', NBOUNDNODES, ' nodes.'
         WRITE(*,*) '============================================================='
      END IF

   END SUBROUTINE READ_3D_UNSTRUCTURED_GRID_SU2



   SUBROUTINE READ_1D_UNSTRUCTURED_GRID_MSH(FILENAME)

      IMPLICIT NONE

      CHARACTER*256, INTENT(IN) :: FILENAME

      CHARACTER*256 :: LINE, GROUPNAME, DUMMYLINE

      INTEGER, PARAMETER :: in5 = 2385
      INTEGER            :: ios
      INTEGER            :: ReasonEOF

      REAL :: MESH_VERSION
      INTEGER :: FILE_TYPE, DATA_SIZE, DUMMY, CELL_COUNT
      INTEGER, DIMENSION(:), ALLOCATABLE :: PG_MAP, PG_PREMAP, TEMP_CELL_PG
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: TEMP_CELL_NODES

      INTEGER            :: NUM, I, J, FOUND, V1, V2, ELEM_TYPE, NUMELEMS
      REAL(KIND=8)       :: X1, X2
      REAL(KIND=8), DIMENSION(3) :: XYZ, A, B

      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: TEMP_CELL_NEIGHBORS

      INTEGER, DIMENSION(2) :: VLIST2

      INTEGER, DIMENSION(:), ALLOCATABLE :: N_CELLS_WITH_NODE, CELL_WITH_NODE, IOF
      INTEGER :: IDX, JN, JC1, JC2, IPG

      LOGICAL, DIMENSION(:), ALLOCATABLE :: NODE_ON_BOUNDARY
      INTEGER :: NUM_BOUNDARY_NODES, NUM_BOUNDARY_ELEM

      ! Open input file for reading
      OPEN(UNIT=in5,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios .NE. 0) THEN
         CALL ERROR_ABORT('Attention, mesh file not found! ABORTING.')
      ENDIF

      ! Read the mesh file. MSH file format (*.msh), verson 2, ASCII only
      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Reading grid file in MSH2 format.'
         WRITE(*,*) '==========================================='
      END IF
      
      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE
         IF (ReasonEOF < 0) EXIT 

         IF (TRIM(LINE) == '$MeshFormat') THEN
            READ(in5,*, IOSTAT=ReasonEOF) MESH_VERSION, FILE_TYPE, DATA_SIZE

            IF (INT(MESH_VERSION) .NE. 2) THEN
               CALL ERROR_ABORT('Attention, only MSH version 2 is supported! ABORTING.')
            END IF
            IF (FILE_TYPE .NE. 0) THEN
               CALL ERROR_ABORT('Attention, only ASCII MSH files are supported! ABORTING.')
            END IF
         END IF

         IF (TRIM(LINE) == '$PhysicalNames') THEN
            READ(in5,*, IOSTAT=ReasonEOF) NUM

            N_GRID_BC = NUM
            ALLOCATE(GRID_BC(NUM))

            ALLOCATE(PG_PREMAP(NUM))
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) ELEM_TYPE, IPG, GRID_BC(I)%PHYSICAL_GROUP_NAME
               PG_PREMAP(I) = IPG
            END DO
            ALLOCATE(PG_MAP(MINVAL(PG_PREMAP):MAXVAL(PG_PREMAP)))
            PG_MAP = -1
            DO I = 1, NUM
               PG_MAP(PG_PREMAP(I)) = I
            END DO

            DEALLOCATE(PG_PREMAP)
         END IF

         IF (TRIM(LINE) == '$Nodes') THEN
            READ(in5,*, IOSTAT=ReasonEOF) NUM

            ALLOCATE(U1D_GRID%NODE_COORDS(3,NUM))
            U1D_GRID%NUM_NODES = NUM

            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) IDX, XYZ
               U1D_GRID%NODE_COORDS(:,I) = XYZ
            END DO
         END IF


         IF (TRIM(LINE) == '$Elements') THEN
            READ(in5,*, IOSTAT=ReasonEOF) NUM

            ALLOCATE(TEMP_CELL_NODES(2, NUM))
            ALLOCATE(TEMP_CELL_PG(NUM))
            CELL_COUNT = 0

            DO I = 1, NUM
               READ(in5,'(A)', IOSTAT=ReasonEOF) LINE
               READ(LINE,*) IDX, ELEM_TYPE

               IF (ELEM_TYPE .EQ. 1) THEN
                  CELL_COUNT = CELL_COUNT + 1
                  READ(LINE,*) IDX, ELEM_TYPE, DUMMY, IPG, DUMMY, TEMP_CELL_NODES(:,CELL_COUNT)
                  TEMP_CELL_PG(CELL_COUNT) = PG_MAP(IPG)
               END IF
            END DO



            ALLOCATE(U1D_GRID%CELL_NODES(2,CELL_COUNT))
            ALLOCATE(U1D_GRID%CELL_EDGES_PG(2, CELL_COUNT))
            ALLOCATE(U1D_GRID%CELL_PG(CELL_COUNT))
            U1D_GRID%NUM_CELLS = CELL_COUNT
            U1D_GRID%CELL_EDGES_PG = -1
            U1D_GRID%CELL_PG = -1

            U1D_GRID%CELL_NODES(:, :) = TEMP_CELL_NODES(:, 1:CELL_COUNT)
            U1D_GRID%CELL_PG(:)       = TEMP_CELL_PG(1:CELL_COUNT)

            DEALLOCATE(TEMP_CELL_NODES)
            DEALLOCATE(TEMP_CELL_PG)
         END IF
      END DO

      REWIND(in5)

      ALLOCATE(N_CELLS_WITH_NODE(U1D_GRID%NUM_NODES))
      ALLOCATE(IOF(U1D_GRID%NUM_NODES))

      N_CELLS_WITH_NODE = 0
      DO I = 1, U1D_GRID%NUM_CELLS
         DO V1 = 1, 2
            JN = U1D_GRID%CELL_NODES(V1,I)
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO
   
      IOF = -1
      IDX = 1
      DO JN = 1, U1D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) .NE. 0) THEN
            IOF(JN) = IDX
            IDX = IDX + N_CELLS_WITH_NODE(JN)
         END IF
      END DO
   
      ALLOCATE(CELL_WITH_NODE(IDX))
      
      N_CELLS_WITH_NODE = 0
      DO I = 1, U1D_GRID%NUM_CELLS
         DO V1 = 1, 2
            JN = U1D_GRID%CELL_NODES(V1,I)
            CELL_WITH_NODE(IOF(JN) + N_CELLS_WITH_NODE(JN)) = I
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO


      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE
         IF (ReasonEOF < 0) EXIT

         IF (TRIM(LINE) == '$Elements') THEN
            READ(in5,*, IOSTAT=ReasonEOF) NUM

            DO I = 1, NUM
               
               READ(in5,'(A)', IOSTAT=ReasonEOF) LINE
               READ(LINE,*) IDX, ELEM_TYPE

               IF (ELEM_TYPE == 15) THEN ! element in physical group is a vertex.

                  READ(LINE,*) IDX, ELEM_TYPE, DUMMY, IPG, DUMMY, JN
                  
                  IF (N_CELLS_WITH_NODE(JN) > 0) THEN
                     DO IDX = 0, N_CELLS_WITH_NODE(JN) - 1
                        JC1 = CELL_WITH_NODE(IOF(JN) + IDX)
                        FOUND = 0
                        DO V1 = 1, 2
                           IF (U1D_GRID%CELL_NODES(V1,JC1) == JN) THEN
                              U1D_GRID%CELL_EDGES_PG(V1, JC1) = PG_MAP(IPG)
                           END IF
                        END DO
                     END DO
                  END IF
               ELSE IF (ELEM_TYPE == 1) THEN
                  ! Element in the physical group is a line.
               ELSE
                  WRITE(*,*) 'Error! element type was not point or line.'
               END IF


            END DO

         END IF
      END DO

      ! Done reading
      CLOSE(in5)
      DEALLOCATE(PG_MAP)


      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing cell volumes.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute cell volumes
      ALLOCATE(U1D_GRID%SEGMENT_LENGTHS(U1D_GRID%NUM_CELLS))
      ALLOCATE(U1D_GRID%CELL_VOLUMES(U1D_GRID%NUM_CELLS))
      DO I = 1, U1D_GRID%NUM_CELLS
         A = U1D_GRID%NODE_COORDS(:, U1D_GRID%CELL_NODES(1,I))
         B = U1D_GRID%NODE_COORDS(:, U1D_GRID%CELL_NODES(2,I))

         U1D_GRID%SEGMENT_LENGTHS(I) = ABS(A(1)-B(1))
         U1D_GRID%CELL_VOLUMES(I) = U1D_GRID%SEGMENT_LENGTHS(I) * (YMAX-YMIN) * (ZMAX-ZMIN)
      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing grid connectivity.'
         WRITE(*,*) '==========================================='
      END IF

      ! Find cell connectivity
      ALLOCATE(TEMP_CELL_NEIGHBORS(2, U1D_GRID%NUM_CELLS))
      TEMP_CELL_NEIGHBORS = -1



      DO JN = 1, U1D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) > 1) THEN
            DO I = 0, N_CELLS_WITH_NODE(JN) - 1
               DO J = I, N_CELLS_WITH_NODE(JN) - 1
                  IF (I == J) CYCLE
                  JC1 = CELL_WITH_NODE(IOF(JN) + I)
                  JC2 = CELL_WITH_NODE(IOF(JN) + J)

                  FOUND = 0
                  DO V1 = 1, 2
                     DO V2 = 1, 2
                        IF (U1D_GRID%CELL_NODES(V1,JC1) == U1D_GRID%CELL_NODES(V2,JC2)) THEN
                           FOUND = FOUND + 1
                           IF (FOUND .GT. 1) CALL ERROR_ABORT('Error! Found duplicate cells in the mesh!')
                           TEMP_CELL_NEIGHBORS(V1, JC1) = JC2
                           TEMP_CELL_NEIGHBORS(V2, JC2) = JC1
                        END IF
                     END DO
                  END DO

               END DO
            END DO
         END IF
      END DO

      U1D_GRID%CELL_NEIGHBORS = TEMP_CELL_NEIGHBORS



      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing face normals.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute segment edge normals
      ALLOCATE(U1D_GRID%EDGE_NORMAL(3, 2, U1D_GRID%NUM_CELLS))
      DO I = 1, U1D_GRID%NUM_CELLS

         U1D_GRID%EDGE_NORMAL(1,1,I) = -1.d0
         U1D_GRID%EDGE_NORMAL(2,1,I) =  0.d0
         U1D_GRID%EDGE_NORMAL(3,1,I) =  0.d0

         U1D_GRID%EDGE_NORMAL(1,2,I) =  1.d0
         U1D_GRID%EDGE_NORMAL(2,2,I) =  0.d0
         U1D_GRID%EDGE_NORMAL(3,2,I) =  0.d0

      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Checking ordering.'
         WRITE(*,*) '==========================================='
      END IF

      DO I = 1, U1D_GRID%NUM_CELLS
         X1 = U1D_GRID%NODE_COORDS(1, U1D_GRID%CELL_NODES(2,I)) &
            - U1D_GRID%NODE_COORDS(1, U1D_GRID%CELL_NODES(1,I))

         IF (X1 < 0) CALL ERROR_ABORT('1D mesh segment are reversed.')
      END DO

      NCELLS = U1D_GRID%NUM_CELLS
      NNODES = U1D_GRID%NUM_NODES



      ALLOCATE(U1D_GRID%BASIS_COEFFS(2,2,NCELLS))

      DO I = 1, NCELLS
         V1 = U1D_GRID%CELL_NODES(1,I)
         V2 = U1D_GRID%CELL_NODES(2,I)

         X1 = U1D_GRID%NODE_COORDS(1, V1)
         X2 = U1D_GRID%NODE_COORDS(1, V2)

         ! These are such that PSI_i = x * BASIS_COEFFS(1,i,IC) + BASIS_COEFFS(2,i,IC)

         U1D_GRID%BASIS_COEFFS(1,1,I) = -1.d0
         U1D_GRID%BASIS_COEFFS(2,1,I) =  X2

         U1D_GRID%BASIS_COEFFS(1,2,I) =  1.d0
         U1D_GRID%BASIS_COEFFS(2,2,I) = -X1
         

         U1D_GRID%BASIS_COEFFS(:,:,I) = U1D_GRID%BASIS_COEFFS(:,:,I)/U1D_GRID%SEGMENT_LENGTHS(I)

      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Creating boundary grid.'
         WRITE(*,*) '==========================================='
      END IF

      ALLOCATE(U1D_GRID%SEGMENT_NODES_BOUNDARY_INDEX(2,NCELLS))
      U1D_GRID%SEGMENT_NODES_BOUNDARY_INDEX = -1
      ALLOCATE(NODE_ON_BOUNDARY(NNODES))
      NODE_ON_BOUNDARY = .FALSE.
      ALLOCATE(U1D_GRID%NODES_BOUNDARY_INDEX(NNODES))
      U1D_GRID%NODES_BOUNDARY_INDEX = -1
      NUM_BOUNDARY_NODES = 0
      NUM_BOUNDARY_ELEM = 0
      DO I = 1, NCELLS
         DO J = 1, 2
            ! If the vertex belongs to any physical group, it should be part of the boundary grid
            ! Later, we may want to filter this further
            IF (U1D_GRID%CELL_EDGES_PG(J,I) .NE. -1) THEN
               NUM_BOUNDARY_ELEM = NUM_BOUNDARY_ELEM + 1
               V1 = U1D_GRID%CELL_NODES(J, I)
               IF (.NOT. NODE_ON_BOUNDARY(V1)) THEN
                  NUM_BOUNDARY_NODES = NUM_BOUNDARY_NODES + 1
                  U1D_GRID%NODES_BOUNDARY_INDEX(V1) = NUM_BOUNDARY_NODES
                  NODE_ON_BOUNDARY(V1) = .TRUE.
               END IF
            END IF
         END DO
      END DO

      U0D_GRID%NUM_POINTS = NUM_BOUNDARY_ELEM
      U0D_GRID%NUM_NODES = NUM_BOUNDARY_NODES
      ALLOCATE(U0D_GRID%POINT_NODES(NUM_BOUNDARY_ELEM))
      ALLOCATE(U0D_GRID%POINT_PG(NUM_BOUNDARY_ELEM))
      ALLOCATE(U0D_GRID%NODE_COORDS(3, NUM_BOUNDARY_NODES))

      DO I = 1, NNODES
         IF (NODE_ON_BOUNDARY(I)) THEN
            U0D_GRID%NODE_COORDS(:,U1D_GRID%NODES_BOUNDARY_INDEX(I)) = U1D_GRID%NODE_COORDS(:,I)
         END IF
      END DO

      NUM_BOUNDARY_ELEM = 0

      DO I = 1, NCELLS
         DO J = 1, 2
            IF (U1D_GRID%CELL_EDGES_PG(J,I) .NE. -1) THEN
               NUM_BOUNDARY_ELEM = NUM_BOUNDARY_ELEM + 1
               U0D_GRID%POINT_PG(NUM_BOUNDARY_ELEM) = U1D_GRID%CELL_EDGES_PG(J,I)
               U1D_GRID%SEGMENT_NODES_BOUNDARY_INDEX(J,I) = NUM_BOUNDARY_ELEM

               V1 = U1D_GRID%CELL_NODES(J, I)
               U0D_GRID%POINT_NODES(NUM_BOUNDARY_ELEM) = U1D_GRID%NODES_BOUNDARY_INDEX(V1)

            END IF
         END DO
      END DO
      
      DEALLOCATE(NODE_ON_BOUNDARY)

      NBOUNDCELLS = NUM_BOUNDARY_ELEM
      NBOUNDNODES = NUM_BOUNDARY_NODES

      ! Compute areas and lengths of boundary mesh
      ALLOCATE(U0D_GRID%VERTEX_AREAS(U0D_GRID%NUM_POINTS))
      U0D_GRID%VERTEX_AREAS = (YMAX-YMIN) * (ZMAX-ZMIN)

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '============================================================='
         WRITE(*,*) 'Done reading grid file.'
         WRITE(*,*) 'It contains ', NNODES, ' nodes and ', NCELLS, ' cells.'
         WRITE(*,*) 'The boundary grid contains ', NBOUNDCELLS, ' lines and ', NBOUNDNODES, ' nodes.'
         WRITE(*,*) '============================================================='
      END IF

   END SUBROUTINE READ_1D_UNSTRUCTURED_GRID_MSH


   SUBROUTINE READ_2D_UNSTRUCTURED_GRID_MSH(FILENAME)

      IMPLICIT NONE

      CHARACTER*256, INTENT(IN) :: FILENAME

      CHARACTER*256 :: LINE, GROUPNAME, DUMMYLINE

      INTEGER, PARAMETER :: in5 = 2385
      INTEGER            :: ios
      INTEGER            :: ReasonEOF

      REAL :: MESH_VERSION
      INTEGER :: FILE_TYPE, DATA_SIZE, DUMMY, CELL_COUNT
      INTEGER, DIMENSION(:), ALLOCATABLE :: PG_MAP, PG_PREMAP, TEMP_CELL_PG
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: TEMP_CELL_NODES

      INTEGER            :: NUM, I, J, FOUND, V1, V2, V3, ELEM_TYPE, NUMELEMS
      INTEGER, DIMENSION(3,2) :: IND
      REAL(KIND=8)       :: X1, X2, X3, Y1, Y2, Y3, LEN, RAD
      REAL(KIND=8), DIMENSION(3) :: XYZ, A, B, C

      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: TEMP_CELL_NEIGHBORS

      INTEGER, DIMENSION(2) :: VLIST2, WHICH1, WHICH2
      INTEGER, DIMENSION(3) :: VLIST3

      INTEGER, DIMENSION(:), ALLOCATABLE :: N_CELLS_WITH_NODE, CELL_WITH_NODE, IOF
      INTEGER :: IDX, JN, JC1, JC2, IPG

      LOGICAL, DIMENSION(:), ALLOCATABLE :: NODE_ON_BOUNDARY
      INTEGER :: NUM_BOUNDARY_NODES, NUM_BOUNDARY_ELEM

      ! Open input file for reading
      OPEN(UNIT=in5,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios .NE. 0) THEN
         CALL ERROR_ABORT('Attention, mesh file not found! ABORTING.')
      ENDIF

      ! Read the mesh file. MSH file format (*.msh), verson 2, ASCII only
      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Reading grid file in MSH2 format.'
         WRITE(*,*) '==========================================='
      END IF

      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE
         IF (ReasonEOF < 0) EXIT 

         IF (TRIM(LINE) == '$MeshFormat') THEN
            READ(in5,*, IOSTAT=ReasonEOF) MESH_VERSION, FILE_TYPE, DATA_SIZE

            IF (INT(MESH_VERSION) .NE. 2) THEN
               CALL ERROR_ABORT('Attention, only MSH version 2 is supported! ABORTING.')
            END IF
            IF (FILE_TYPE .NE. 0) THEN
               CALL ERROR_ABORT('Attention, only ASCII MSH files are supported! ABORTING.')
            END IF
         END IF

         IF (TRIM(LINE) == '$PhysicalNames') THEN
            READ(in5,*, IOSTAT=ReasonEOF) NUM

            N_GRID_BC = NUM
            ALLOCATE(GRID_BC(NUM))

            ALLOCATE(PG_PREMAP(NUM))
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) ELEM_TYPE, IPG, GRID_BC(I)%PHYSICAL_GROUP_NAME
               PG_PREMAP(I) = IPG
            END DO
            ALLOCATE(PG_MAP(MINVAL(PG_PREMAP):MAXVAL(PG_PREMAP)))
            PG_MAP = -1
            DO I = 1, NUM
               PG_MAP(PG_PREMAP(I)) = I
            END DO

            DEALLOCATE(PG_PREMAP)
         END IF

         IF (TRIM(LINE) == '$Nodes') THEN
            READ(in5,*, IOSTAT=ReasonEOF) NUM

            ALLOCATE(U2D_GRID%NODE_COORDS(3,NUM))
            U2D_GRID%NUM_NODES = NUM

            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) IDX, XYZ
               U2D_GRID%NODE_COORDS(:,I) = XYZ
            END DO
         END IF


         IF (TRIM(LINE) == '$Elements') THEN
            READ(in5,*, IOSTAT=ReasonEOF) NUM

            ALLOCATE(TEMP_CELL_NODES(3, NUM))
            ALLOCATE(TEMP_CELL_PG(NUM))
            CELL_COUNT = 0

            DO I = 1, NUM
               READ(in5,'(A)', IOSTAT=ReasonEOF) LINE
               READ(LINE,*) IDX, ELEM_TYPE

               IF (ELEM_TYPE .EQ. 2) THEN
                  CELL_COUNT = CELL_COUNT + 1
                  READ(LINE,*) IDX, ELEM_TYPE, DUMMY, IPG, DUMMY, TEMP_CELL_NODES(:,CELL_COUNT)
                  TEMP_CELL_PG(CELL_COUNT) = PG_MAP(IPG)
               END IF
            END DO



            ALLOCATE(U2D_GRID%CELL_NODES(3,CELL_COUNT))
            ALLOCATE(U2D_GRID%CELL_EDGES_PG(3, CELL_COUNT))
            ALLOCATE(U2D_GRID%CELL_PG(CELL_COUNT))
            U2D_GRID%NUM_CELLS = CELL_COUNT
            U2D_GRID%CELL_EDGES_PG = -1
            U2D_GRID%CELL_PG = -1

            U2D_GRID%CELL_NODES(:, :) = TEMP_CELL_NODES(:, 1:CELL_COUNT)
            U2D_GRID%CELL_PG(:)       = TEMP_CELL_PG(1:CELL_COUNT)

            DEALLOCATE(TEMP_CELL_NODES)
            DEALLOCATE(TEMP_CELL_PG)
         END IF
      END DO

      REWIND(in5)

      ALLOCATE(N_CELLS_WITH_NODE(U2D_GRID%NUM_NODES))
      ALLOCATE(IOF(U2D_GRID%NUM_NODES))

      N_CELLS_WITH_NODE = 0
      DO I = 1, U2D_GRID%NUM_CELLS
         DO V1 = 1, 3
            JN = U2D_GRID%CELL_NODES(V1,I)
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO
   
      IOF = -1
      IDX = 1
      DO JN = 1, U2D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) .NE. 0) THEN
            IOF(JN) = IDX
            IDX = IDX + N_CELLS_WITH_NODE(JN)
         END IF
      END DO
   
      ALLOCATE(CELL_WITH_NODE(IDX))
      
      N_CELLS_WITH_NODE = 0
      DO I = 1, U2D_GRID%NUM_CELLS
         DO V1 = 1, 3
            JN = U2D_GRID%CELL_NODES(V1,I)
            CELL_WITH_NODE(IOF(JN) + N_CELLS_WITH_NODE(JN)) = I
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO


      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE
         IF (ReasonEOF < 0) EXIT

         IF (TRIM(LINE) == '$Elements') THEN
            READ(in5,*, IOSTAT=ReasonEOF) NUM

            ! Assign physical groups to cell edges.
            DO I = 1, NUM

               READ(in5,'(A)', IOSTAT=ReasonEOF) LINE
               READ(LINE,*) IDX, ELEM_TYPE

               IF (ELEM_TYPE == 1) THEN ! element in physical group is a cell edge (segment).

                  READ(LINE,*) IDX, ELEM_TYPE, DUMMY, IPG, DUMMY, VLIST2
                  
                  JN = VLIST2(1)
                  IF (N_CELLS_WITH_NODE(JN) > 0) THEN
                     DO IDX = 0, N_CELLS_WITH_NODE(JN) - 1
                        JC1 = CELL_WITH_NODE(IOF(JN) + IDX)
                        FOUND = 0
                        DO V1 = 1, 3
                           IF (ANY(VLIST2 == U2D_GRID%CELL_NODES(V1,JC1))) THEN
                              FOUND = FOUND + 1
                              WHICH1(FOUND) = V1
                           END IF
                        END DO
         
                        IF (FOUND == 2) THEN
                           IF (ANY(WHICH1 == 1) .AND. ANY(WHICH1 == 2)) THEN
                              U2D_GRID%CELL_EDGES_PG(1, JC1) = PG_MAP(IPG)
                           ELSE IF (ANY(WHICH1 == 2) .AND. ANY(WHICH1 == 3)) THEN
                              U2D_GRID%CELL_EDGES_PG(2, JC1) = PG_MAP(IPG)
                           ELSE IF (ANY(WHICH1 == 3) .AND. ANY(WHICH1 == 1)) THEN
                              U2D_GRID%CELL_EDGES_PG(3, JC1) = PG_MAP(IPG)
                           END IF
                        END IF
                     END DO
                  END IF
               ELSE IF (ELEM_TYPE == 2) THEN
                  ! Element in the physical group is a triangle
               ELSE
                  WRITE(*,*) 'Error! element type was not line or triangle.'
               END IF

            END DO

         END IF
      END DO

      ! Done reading
      CLOSE(in5)
      DEALLOCATE(PG_MAP)


      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing cell volumes.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute cell volumes
      ALLOCATE(U2D_GRID%CELL_AREAS(U2D_GRID%NUM_CELLS))
      ALLOCATE(U2D_GRID%CELL_VOLUMES(U2D_GRID%NUM_CELLS))
      DO I = 1, U2D_GRID%NUM_CELLS
         A = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(1,I))
         B = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(2,I))
         C = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(3,I))

         U2D_GRID%CELL_AREAS(I) = 0.5*ABS(A(1)*(B(2)-C(2)) + B(1)*(C(2)-A(2)) + C(1)*(A(2)-B(2)))
         IF (DIMS == 2 .AND. .NOT. AXI) THEN
            U2D_GRID%CELL_VOLUMES(I) = U2D_GRID%CELL_AREAS(I) * (ZMAX-ZMIN)
         END IF
         IF (DIMS == 2 .AND. AXI) THEN
            RAD = (A(2)+B(2)+C(2))/3.
            U2D_GRID%CELL_VOLUMES(I) = U2D_GRID%CELL_AREAS(I) * (ZMAX-ZMIN)*RAD
         END IF
      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing grid connectivity.'
         WRITE(*,*) '==========================================='
      END IF

      ! Find cell connectivity
      ALLOCATE(TEMP_CELL_NEIGHBORS(3, U2D_GRID%NUM_CELLS))
      TEMP_CELL_NEIGHBORS = -1



      DO JN = 1, U2D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) > 1) THEN
            DO I = 0, N_CELLS_WITH_NODE(JN) - 1
               DO J = I, N_CELLS_WITH_NODE(JN) - 1
                  IF (I == J) CYCLE
                  JC1 = CELL_WITH_NODE(IOF(JN) + I)
                  JC2 = CELL_WITH_NODE(IOF(JN) + J)


                  FOUND = 0
                  DO V1 = 1, 3
                     DO V2 = 1, 3
                        IF (U2D_GRID%CELL_NODES(V1,JC1) == U2D_GRID%CELL_NODES(V2,JC2)) THEN
                           FOUND = FOUND + 1
                           IF (FOUND .GT. 2) CALL ERROR_ABORT('Error! Found duplicate cells in the mesh!')
                           WHICH1(FOUND) = V1
                           WHICH2(FOUND) = V2
                        END IF
                     END DO
                  END DO

                  IF (FOUND == 2) THEN
      
                     IF (ANY(WHICH1 == 1) .AND. ANY(WHICH1 == 2)) THEN
                        TEMP_CELL_NEIGHBORS(1, JC1) = JC2
                     ELSE IF (ANY(WHICH1 == 2) .AND. ANY(WHICH1 == 3)) THEN
                        TEMP_CELL_NEIGHBORS(2, JC1) = JC2
                     ELSE IF (ANY(WHICH1 == 3) .AND. ANY(WHICH1 == 1)) THEN
                        TEMP_CELL_NEIGHBORS(3, JC1) = JC2
                     END IF

                     IF (ANY(WHICH2 == 1) .AND. ANY(WHICH2 == 2)) THEN
                        TEMP_CELL_NEIGHBORS(1, JC2) = JC1
                     ELSE IF (ANY(WHICH2 == 2) .AND. ANY(WHICH2 == 3)) THEN
                        TEMP_CELL_NEIGHBORS(2, JC2) = JC1
                     ELSE IF (ANY(WHICH2 == 3) .AND. ANY(WHICH2 == 1)) THEN
                        TEMP_CELL_NEIGHBORS(3, JC2) = JC1
                     END IF
      
                  END IF


               END DO
            END DO
         END IF
      END DO

      U2D_GRID%CELL_NEIGHBORS = TEMP_CELL_NEIGHBORS



      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing face normals.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute cell edge normals
      IND(1,:) = [1,2]
      IND(2,:) = [2,3]
      IND(3,:) = [3,1]
      ALLOCATE(U2D_GRID%EDGE_NORMAL(3, 3, U2D_GRID%NUM_CELLS))
      ALLOCATE(U2D_GRID%CELL_EDGES_LEN(3, U2D_GRID%NUM_CELLS))
      DO I = 1, U2D_GRID%NUM_CELLS
         DO J = 1, 3
            X1 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(IND(J,1),I))
            X2 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(IND(J,2),I))
            Y1 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(IND(J,1),I))
            Y2 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(IND(J,2),I))
            LEN = SQRT((Y2-Y1)*(Y2-Y1) + (X2-X1)*(X2-X1))
            U2D_GRID%CELL_EDGES_LEN(J,I) = LEN
            U2D_GRID%EDGE_NORMAL(1,J,I) = (Y2-Y1)/LEN
            U2D_GRID%EDGE_NORMAL(2,J,I) = (X1-X2)/LEN
            U2D_GRID%EDGE_NORMAL(3,J,I) = 0.d0
         END DO
      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Checking ordering.'
         WRITE(*,*) '==========================================='
      END IF

      DO I = 1, U2D_GRID%NUM_CELLS
         X1 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(2,I)) &
            - U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(1,I))
         X2 = U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(3,I)) &
            - U2D_GRID%NODE_COORDS(1, U2D_GRID%CELL_NODES(1,I))
         Y1 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(2,I)) &
            - U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(1,I))
         Y2 = U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(3,I)) &
            - U2D_GRID%NODE_COORDS(2, U2D_GRID%CELL_NODES(1,I))

         IF (X1*Y2-X2*Y1 < 0) CALL ERROR_ABORT('2D mesh triangles have negative z-normal.')

      END DO

      NCELLS = U2D_GRID%NUM_CELLS
      NNODES = U2D_GRID%NUM_NODES



      ALLOCATE(U2D_GRID%BASIS_COEFFS(3,3,NCELLS))

      DO I = 1, NCELLS
         V1 = U2D_GRID%CELL_NODES(1,I)
         V2 = U2D_GRID%CELL_NODES(2,I)
         V3 = U2D_GRID%CELL_NODES(3,I)

         X1 = U2D_GRID%NODE_COORDS(1, V1)
         X2 = U2D_GRID%NODE_COORDS(1, V2)
         X3 = U2D_GRID%NODE_COORDS(1, V3)
         Y1 = U2D_GRID%NODE_COORDS(2, V1)
         Y2 = U2D_GRID%NODE_COORDS(2, V2)
         Y3 = U2D_GRID%NODE_COORDS(2, V3)

         ! These are such that PSI_i = SUM_j [ x_j * BASIS_COEFFS(j,i,IC) ] + BASIS_COEFFS(3,i,IC)

         U2D_GRID%BASIS_COEFFS(1,1,I) =  Y2-Y3
         U2D_GRID%BASIS_COEFFS(2,1,I) = -(X2-X3)
         U2D_GRID%BASIS_COEFFS(3,1,I) =  X2*Y3 - X3*Y2

         U2D_GRID%BASIS_COEFFS(1,2,I) = -(Y1-Y3)
         U2D_GRID%BASIS_COEFFS(2,2,I) =  X1-X3
         U2D_GRID%BASIS_COEFFS(3,2,I) =  X3*Y1 - X1*Y3

         U2D_GRID%BASIS_COEFFS(1,3,I) = -(Y2-Y1)
         U2D_GRID%BASIS_COEFFS(2,3,I) =  X2-X1
         U2D_GRID%BASIS_COEFFS(3,3,I) =  X1*Y2 - X2*Y1

         U2D_GRID%BASIS_COEFFS(:,:,I) = 0.5*U2D_GRID%BASIS_COEFFS(:,:,I)/U2D_GRID%CELL_AREAS(I)

      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Creating boundary grid.'
         WRITE(*,*) '==========================================='
      END IF

      ALLOCATE(U2D_GRID%CELL_EDGES_BOUNDARY_INDEX(3,NCELLS))
      U2D_GRID%CELL_EDGES_BOUNDARY_INDEX = -1
      ALLOCATE(NODE_ON_BOUNDARY(NNODES))
      NODE_ON_BOUNDARY = .FALSE.
      ALLOCATE(U2D_GRID%NODES_BOUNDARY_INDEX(NNODES))
      U2D_GRID%NODES_BOUNDARY_INDEX = -1
      NUM_BOUNDARY_NODES = 0
      NUM_BOUNDARY_ELEM = 0
      DO I = 1, NCELLS
         DO J = 1, 3
            ! If the edge belongs to any physical group, it should be part of the boundary grid
            ! Later, we may want to filter this further
            IF (U2D_GRID%CELL_EDGES_PG(J,I) .NE. -1) THEN
               NUM_BOUNDARY_ELEM = NUM_BOUNDARY_ELEM + 1
               V1 = U2D_GRID%CELL_NODES(J, I)
               IF (J == 3) THEN
                  V2 = U2D_GRID%CELL_NODES(1, I)
               ELSE
                  V2 = U2D_GRID%CELL_NODES(J+1, I)
               END IF
               IF (.NOT. NODE_ON_BOUNDARY(V1)) THEN
                  NUM_BOUNDARY_NODES = NUM_BOUNDARY_NODES + 1
                  U2D_GRID%NODES_BOUNDARY_INDEX(V1) = NUM_BOUNDARY_NODES
                  NODE_ON_BOUNDARY(V1) = .TRUE.
               END IF
               IF (.NOT. NODE_ON_BOUNDARY(V2)) THEN
                  NUM_BOUNDARY_NODES = NUM_BOUNDARY_NODES + 1
                  U2D_GRID%NODES_BOUNDARY_INDEX(V2) = NUM_BOUNDARY_NODES
                  NODE_ON_BOUNDARY(V2) = .TRUE.
               END IF
               
            END IF
         END DO
      END DO      

      U1D_GRID%NUM_CELLS = NUM_BOUNDARY_ELEM
      U1D_GRID%NUM_NODES = NUM_BOUNDARY_NODES
      ALLOCATE(U1D_GRID%CELL_NODES(2, NUM_BOUNDARY_ELEM))
      ALLOCATE(U1D_GRID%CELL_PG(NUM_BOUNDARY_ELEM))
      ALLOCATE(U1D_GRID%NODE_COORDS(3, NUM_BOUNDARY_NODES))

      DO I = 1, NNODES
         IF (NODE_ON_BOUNDARY(I)) THEN
            U1D_GRID%NODE_COORDS(:,U2D_GRID%NODES_BOUNDARY_INDEX(I)) = U2D_GRID%NODE_COORDS(:,I)
         END IF
      END DO

      NUM_BOUNDARY_ELEM = 0

      DO I = 1, NCELLS
         DO J = 1, 3
            IF (U2D_GRID%CELL_EDGES_PG(J,I) .NE. -1) THEN
               NUM_BOUNDARY_ELEM = NUM_BOUNDARY_ELEM + 1
               U1D_GRID%CELL_PG(NUM_BOUNDARY_ELEM) = U2D_GRID%CELL_EDGES_PG(J,I)
               U2D_GRID%CELL_EDGES_BOUNDARY_INDEX(J,I) = NUM_BOUNDARY_ELEM

               V1 = U2D_GRID%CELL_NODES(J, I)
               IF (J == 3) THEN
                  V2 = U2D_GRID%CELL_NODES(1, I)
               ELSE
                  V2 = U2D_GRID%CELL_NODES(J+1, I)
               END IF
               U1D_GRID%CELL_NODES(1, NUM_BOUNDARY_ELEM) = U2D_GRID%NODES_BOUNDARY_INDEX(V1)
               U1D_GRID%CELL_NODES(2, NUM_BOUNDARY_ELEM) = U2D_GRID%NODES_BOUNDARY_INDEX(V2)

            END IF
         END DO
      END DO
      
      DEALLOCATE(NODE_ON_BOUNDARY)

      NBOUNDCELLS = NUM_BOUNDARY_ELEM
      NBOUNDNODES = NUM_BOUNDARY_NODES

      ! Compute areas and lengths of boundary mesh
      ALLOCATE(U1D_GRID%SEGMENT_LENGTHS(U1D_GRID%NUM_CELLS))
      ALLOCATE(U1D_GRID%SEGMENT_AREAS(U1D_GRID%NUM_CELLS))
      DO I = 1, U1D_GRID%NUM_CELLS
         A = U1D_GRID%NODE_COORDS(:, U1D_GRID%CELL_NODES(1,I))
         B = U1D_GRID%NODE_COORDS(:, U1D_GRID%CELL_NODES(2,I))

         U1D_GRID%SEGMENT_LENGTHS(I) = SQRT((B(1) - A(1))**2 + (B(2) - A(2))**2)
         IF (.NOT. AXI) THEN
            U1D_GRID%SEGMENT_AREAS(I) = U1D_GRID%SEGMENT_LENGTHS(I) * (ZMAX-ZMIN)
         ELSE
            RAD = 0.5*(A(2)+B(2))
            U1D_GRID%SEGMENT_AREAS(I) = U1D_GRID%SEGMENT_LENGTHS(I) * (ZMAX-ZMIN)*RAD
         END IF
      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '============================================================='
         WRITE(*,*) 'Done reading grid file.'
         WRITE(*,*) 'It contains ', NNODES, ' nodes and ', NCELLS, ' cells.'
         WRITE(*,*) 'The boundary grid contains ', NBOUNDCELLS, ' lines and ', NBOUNDNODES, ' nodes.'
         WRITE(*,*) '============================================================='
      END IF
      
   END SUBROUTINE READ_2D_UNSTRUCTURED_GRID_MSH


   SUBROUTINE READ_3D_UNSTRUCTURED_GRID_MSH(FILENAME)

      IMPLICIT NONE

      CHARACTER*256, INTENT(IN) :: FILENAME

      CHARACTER*256 :: LINE, GROUPNAME, DUMMYLINE

      INTEGER, PARAMETER :: in5 = 2385
      INTEGER            :: ios
      INTEGER            :: ReasonEOF

      REAL :: MESH_VERSION
      INTEGER :: FILE_TYPE, DATA_SIZE, DUMMY, CELL_COUNT
      INTEGER, DIMENSION(:), ALLOCATABLE :: PG_MAP, PG_PREMAP, TEMP_CELL_PG
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: TEMP_CELL_NODES

      INTEGER            :: NUM, I, J, FOUND, V1, V2, V3, V4, ELEM_TYPE, NUMELEMS
      INTEGER, DIMENSION(4,3) :: IND
      REAL(KIND=8), DIMENSION(3) :: XYZ, A, B, C, CROSSP

      INTEGER, DIMENSION(:,:), ALLOCATABLE      :: TEMP_CELL_NEIGHBORS

      INTEGER, DIMENSION(3) :: VLIST3, WHICH1, WHICH2
      INTEGER, DIMENSION(4) :: VLIST4

      INTEGER, DIMENSION(:), ALLOCATABLE :: N_CELLS_WITH_NODE, CELL_WITH_NODE, IOF
      INTEGER :: IDX, JN, JC1, JC2, IPG

      REAL(KIND=8) :: VOLUME, X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4

      LOGICAL, DIMENSION(:), ALLOCATABLE :: NODE_ON_BOUNDARY
      INTEGER :: NUM_BOUNDARY_NODES, NUM_BOUNDARY_ELEM

      ! Open input file for reading
      OPEN(UNIT=in5,FILE=FILENAME, STATUS='old',IOSTAT=ios)

      IF (ios .NE. 0) THEN
         CALL ERROR_ABORT('Attention, mesh file not found! ABORTING.')
      ENDIF

      ! Read the mesh file. MSH file format (*.msh), verson 2, ASCII only
      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Reading grid file in MSH format.'
         WRITE(*,*) '==========================================='
      END IF


      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE
         IF (ReasonEOF < 0) EXIT 

         IF (TRIM(LINE) == '$MeshFormat') THEN
            READ(in5,*, IOSTAT=ReasonEOF) MESH_VERSION, FILE_TYPE, DATA_SIZE

            IF (INT(MESH_VERSION) .NE. 2) THEN
               CALL ERROR_ABORT('Attention, only MSH version 2 is supported! ABORTING.')
            END IF
            IF (FILE_TYPE .NE. 0) THEN
               CALL ERROR_ABORT('Attention, only ASCII MSH files are supported! ABORTING.')
            END IF
         END IF

         IF (TRIM(LINE) == '$PhysicalNames') THEN
            READ(in5,*, IOSTAT=ReasonEOF) NUM

            N_GRID_BC = NUM
            ALLOCATE(GRID_BC(NUM))

            ALLOCATE(PG_PREMAP(NUM))
            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) ELEM_TYPE, IPG, GRID_BC(I)%PHYSICAL_GROUP_NAME
               PG_PREMAP(I) = IPG
            END DO
            ALLOCATE(PG_MAP(MINVAL(PG_PREMAP):MAXVAL(PG_PREMAP)))
            PG_MAP = -1
            DO I = 1, NUM
               PG_MAP(PG_PREMAP(I)) = I
            END DO


            DEALLOCATE(PG_PREMAP)
         END IF

         IF (TRIM(LINE) == '$Nodes') THEN
            READ(in5,*, IOSTAT=ReasonEOF) NUM

            ALLOCATE(U3D_GRID%NODE_COORDS(3,NUM))
            U3D_GRID%NUM_NODES = NUM

            DO I = 1, NUM
               READ(in5,*, IOSTAT=ReasonEOF) IDX, XYZ
               U3D_GRID%NODE_COORDS(:,I) = XYZ
            END DO
         END IF


         IF (TRIM(LINE) == '$Elements') THEN
            READ(in5,*, IOSTAT=ReasonEOF) NUM

            ALLOCATE(TEMP_CELL_NODES(4, NUM))
            ALLOCATE(TEMP_CELL_PG(NUM))
            CELL_COUNT = 0

            DO I = 1, NUM
               READ(in5,'(A)', IOSTAT=ReasonEOF) LINE
               READ(LINE,*) IDX, ELEM_TYPE

               IF (ELEM_TYPE .EQ. 4) THEN
                  CELL_COUNT = CELL_COUNT + 1
                  READ(LINE,*) IDX, ELEM_TYPE, DUMMY, IPG, DUMMY, TEMP_CELL_NODES(:,CELL_COUNT)
                  TEMP_CELL_PG(CELL_COUNT) = PG_MAP(IPG)
               END IF
            END DO



            ALLOCATE(U3D_GRID%CELL_NODES(4,CELL_COUNT))
            ALLOCATE(U3D_GRID%CELL_FACES_PG(4, CELL_COUNT))
            ALLOCATE(U3D_GRID%CELL_PG(CELL_COUNT))
            U3D_GRID%NUM_CELLS = CELL_COUNT
            U3D_GRID%CELL_FACES_PG = -1
            U3D_GRID%CELL_PG = -1

            U3D_GRID%CELL_NODES(:, :) = TEMP_CELL_NODES(:, 1:CELL_COUNT)
            U3D_GRID%CELL_PG(:)       = TEMP_CELL_PG(1:CELL_COUNT)

            DEALLOCATE(TEMP_CELL_NODES)
            DEALLOCATE(TEMP_CELL_PG)
         END IF
      END DO

      REWIND(in5)

      
      ALLOCATE(N_CELLS_WITH_NODE(U3D_GRID%NUM_NODES))
      ALLOCATE(IOF(U3D_GRID%NUM_NODES))

      N_CELLS_WITH_NODE = 0
      DO I = 1, U3D_GRID%NUM_CELLS
         DO V1 = 1, 4
            JN = U3D_GRID%CELL_NODES(V1,I)
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO
   
      IOF = -1
      IDX = 1
      DO JN = 1, U3D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) .NE. 0) THEN
            IOF(JN) = IDX
            IDX = IDX + N_CELLS_WITH_NODE(JN)
         END IF
      END DO
   
      ALLOCATE(CELL_WITH_NODE(IDX))
      
      N_CELLS_WITH_NODE = 0
      DO I = 1, U3D_GRID%NUM_CELLS
         DO V1 = 1, 4
            JN = U3D_GRID%CELL_NODES(V1,I)
            CELL_WITH_NODE(IOF(JN) + N_CELLS_WITH_NODE(JN)) = I
            N_CELLS_WITH_NODE(JN) = N_CELLS_WITH_NODE(JN) + 1
         END DO
      END DO


      DO
         READ(in5,*, IOSTAT=ReasonEOF) LINE
         IF (ReasonEOF < 0) EXIT
         
         IF (TRIM(LINE) == '$Elements') THEN
            READ(in5,*, IOSTAT=ReasonEOF) NUM

            ! Assign physical groups to cell edges.
            DO I = 1, NUM
               
               READ(in5,'(A)', IOSTAT=ReasonEOF) LINE
               READ(LINE,*) IDX, ELEM_TYPE

               IF (ELEM_TYPE == 2) THEN ! Element in physical group is a cell (simplex).

                  READ(LINE,*) IDX, ELEM_TYPE, DUMMY, IPG, DUMMY, VLIST3

                  JN = VLIST3(1)
                  IF (N_CELLS_WITH_NODE(JN) > 0) THEN
                     DO IDX = 0, N_CELLS_WITH_NODE(JN) - 1
                        JC1 = CELL_WITH_NODE(IOF(JN) + IDX)
                        FOUND = 0
                        DO V1 = 1, 4
                           IF (ANY(VLIST3 == U3D_GRID%CELL_NODES(V1,JC1))) THEN
                              FOUND = FOUND + 1
                              WHICH1(FOUND) = V1
                           END IF
                        END DO
         
                        IF (FOUND == 3) THEN
            
                           IF (ANY(WHICH1 == 1)) THEN
                              IF (ANY(WHICH1 == 2)) THEN
                                 IF (ANY(WHICH1 == 3))  THEN
                                    U3D_GRID%CELL_FACES_PG(1, JC1) = PG_MAP(IPG)
                                 ELSE IF (ANY(WHICH1 == 4)) THEN
                                    U3D_GRID%CELL_FACES_PG(2, JC1) = PG_MAP(IPG)
                                 END IF
                              ELSE IF (ANY(WHICH1 == 3)) THEN
                                 IF (ANY(WHICH1 == 4)) U3D_GRID%CELL_FACES_PG(4, JC1) = PG_MAP(IPG)
                              END IF
                           ELSE IF (ANY(WHICH1 == 2)) THEN
                              IF (ANY(WHICH1 == 3) .AND. ANY(WHICH1 == 4)) U3D_GRID%CELL_FACES_PG(3, JC1) = PG_MAP(IPG)
                           END IF

                        END IF
                     END DO
                  END IF

               ELSE IF (ELEM_TYPE == 4) THEN 
                  ! element in physical group is a tetrahedron.
               ELSE
                  WRITE(*,*) 'Error! element type was not triangle or prism.'
               END IF
            END DO

         END IF
      END DO

      ! Done reading
      CLOSE(in5)
      DEALLOCATE(PG_MAP)


      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing cell volumes.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute cell volumes
      ALLOCATE(U3D_GRID%CELL_VOLUMES(U3D_GRID%NUM_CELLS))
      DO I = 1, U3D_GRID%NUM_CELLS
         A = U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(2,I)) - U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(1,I))
         B = U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(3,I)) - U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(1,I))
         C = U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(4,I)) - U3D_GRID%NODE_COORDS(:, U3D_GRID%CELL_NODES(1,I))

         U3D_GRID%CELL_VOLUMES(I) = ABS(C(1)*(A(2)*B(3)-A(3)*B(2)) + C(2)*(A(3)*B(1)-A(1)*B(3)) + C(3)*(A(1)*B(2)-A(2)*B(1))) / 6.
      END DO

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing grid connectivity.'
         WRITE(*,*) '==========================================='
      END IF


      ALLOCATE(TEMP_CELL_NEIGHBORS(4, U3D_GRID%NUM_CELLS))
      TEMP_CELL_NEIGHBORS = -1

      DO JN = 1, U3D_GRID%NUM_NODES
         !IF (PROC_ID == 0) WRITE(*,*) 'Checking node ', JN, ' of ',  U3D_GRID%NUM_NODES
         IF (N_CELLS_WITH_NODE(JN) > 1) THEN
            DO I = 0, N_CELLS_WITH_NODE(JN) - 1
               DO J = I, N_CELLS_WITH_NODE(JN) - 1
                  IF (I == J) CYCLE
                  JC1 = CELL_WITH_NODE(IOF(JN) + I)
                  JC2 = CELL_WITH_NODE(IOF(JN) + J)


                  FOUND = 0
                  DO V1 = 1, 4
                     DO V2 = 1, 4
                        IF (U3D_GRID%CELL_NODES(V1,JC1) == U3D_GRID%CELL_NODES(V2,JC2)) THEN
                           FOUND = FOUND + 1
                           IF (FOUND .GT. 3) CALL ERROR_ABORT('Error! Found duplicate cells in the mesh!')
                           WHICH1(FOUND) = V1
                           WHICH2(FOUND) = V2
                        END IF
                     END DO
                  END DO

                  IF (FOUND == 3) THEN
      
                     IF (ANY(WHICH1 == 1)) THEN
                        IF (ANY(WHICH1 == 2)) THEN
                           IF (ANY(WHICH1 == 3))  THEN
                              TEMP_CELL_NEIGHBORS(1, JC1) = JC2
                           ELSE IF (ANY(WHICH1 == 4)) THEN
                              TEMP_CELL_NEIGHBORS(2, JC1) = JC2
                           END IF
                        ELSE IF (ANY(WHICH1 == 3)) THEN
                           IF (ANY(WHICH1 == 4)) TEMP_CELL_NEIGHBORS(4, JC1) = JC2
                        END IF
                     ELSE IF (ANY(WHICH1 == 2)) THEN
                        IF (ANY(WHICH1 == 3) .AND. ANY(WHICH1 == 4)) TEMP_CELL_NEIGHBORS(3, JC1) = JC2
                     END IF
      
      
                     IF (ANY(WHICH2 == 1)) THEN
                        IF (ANY(WHICH2 == 2)) THEN
                           IF (ANY(WHICH2 == 3))  THEN
                              TEMP_CELL_NEIGHBORS(1, JC2) = JC1
                           ELSE IF (ANY(WHICH2 == 4)) THEN
                              TEMP_CELL_NEIGHBORS(2, JC2) = JC1
                           END IF
                        ELSE IF (ANY(WHICH2 == 3)) THEN
                           IF (ANY(WHICH2 == 4)) TEMP_CELL_NEIGHBORS(4, JC2) = JC1
                        END IF
                     ELSE IF (ANY(WHICH2 == 2)) THEN
                        IF (ANY(WHICH2 == 3) .AND. ANY(WHICH2 == 4)) TEMP_CELL_NEIGHBORS(3, JC2) = JC1
                     END IF
      
                  END IF


               END DO
            END DO
         END IF
      END DO

      U3D_GRID%CELL_NEIGHBORS = TEMP_CELL_NEIGHBORS

      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Computing face normals.'
         WRITE(*,*) '==========================================='
      END IF

      ! Compute cell edge normals
      IND(1,:) = [1,3,2]
      IND(2,:) = [1,2,4]
      IND(3,:) = [2,3,4]
      IND(4,:) = [1,4,3]
      ALLOCATE(U3D_GRID%FACE_NORMAL(3, 4, U3D_GRID%NUM_CELLS))
      ALLOCATE(U3D_GRID%FACE_TANG1(3, 4, U3D_GRID%NUM_CELLS))
      ALLOCATE(U3D_GRID%FACE_TANG2(3, 4, U3D_GRID%NUM_CELLS))
      ALLOCATE(U3D_GRID%FACE_NODES(3, 4, U3D_GRID%NUM_CELLS))
      ALLOCATE(U3D_GRID%CELL_FACES_COEFFS(4, 4, U3D_GRID%NUM_CELLS))
      ALLOCATE(U3D_GRID%FACE_AREA(4, U3D_GRID%NUM_CELLS))

      DO I = 1, U3D_GRID%NUM_CELLS
         DO J = 1, 4
            V1 = U3D_GRID%CELL_NODES(IND(J,1),I)
            V2 = U3D_GRID%CELL_NODES(IND(J,2),I)
            V3 = U3D_GRID%CELL_NODES(IND(J,3),I)

            U3D_GRID%FACE_NODES(1,J,I) = V1
            U3D_GRID%FACE_NODES(2,J,I) = V2
            U3D_GRID%FACE_NODES(3,J,I) = V3

            A = U3D_GRID%NODE_COORDS(:,V1)
            B = U3D_GRID%NODE_COORDS(:,V2)
            C = U3D_GRID%NODE_COORDS(:,V3)
            
            CROSSP = CROSS(B-A,C-A)
            U3D_GRID%FACE_NORMAL(:,J,I) = CROSSP/NORM2(CROSSP)
            U3D_GRID%FACE_AREA(J,I) = 0.5*NORM2(CROSSP)
            
            U3D_GRID%FACE_TANG1(:,J,I) = (B-A)/NORM2(B-A)
            U3D_GRID%FACE_TANG2(:,J,I) = CROSS(U3D_GRID%FACE_NORMAL(:,J,I), U3D_GRID%FACE_TANG1(:,J,I))


            ! The coefficients (a,b,c,d) of a*x + b*y + c*z + d = 0
            U3D_GRID%CELL_FACES_COEFFS(1,J,I) =  A(2)*B(3)-B(2)*A(3) &
                                                +B(2)*C(3)-C(2)*B(3) &
                                                +C(2)*A(3)-A(2)*C(3)
            U3D_GRID%CELL_FACES_COEFFS(2,J,I) = -A(1)*B(3)+B(1)*A(3) &
                                                -B(1)*C(3)+C(1)*B(3) &
                                                -C(1)*A(3)+A(1)*C(3)
            U3D_GRID%CELL_FACES_COEFFS(3,J,I) =  A(1)*B(2)-B(1)*A(2) &
                                                +B(1)*C(2)-C(1)*B(2) &
                                                +C(1)*A(2)-A(1)*C(2)
            U3D_GRID%CELL_FACES_COEFFS(4,J,I) = -A(1)*B(2)*C(3) &
                                                +A(1)*C(2)*B(3) &
                                                +B(1)*A(2)*C(3) &
                                                -B(1)*C(2)*A(3) &
                                                -C(1)*A(2)*B(3) &
                                                +C(1)*B(2)*A(3)

         END DO
      END DO



      NCELLS = U3D_GRID%NUM_CELLS
      NNODES = U3D_GRID%NUM_NODES



      ALLOCATE(U3D_GRID%BASIS_COEFFS(4,4,NCELLS))

      DO I = 1, NCELLS
         VOLUME = U3D_GRID%CELL_VOLUMES(I)
         V1 = U3D_GRID%CELL_NODES(1,I)
         V2 = U3D_GRID%CELL_NODES(2,I)
         V3 = U3D_GRID%CELL_NODES(3,I)
         V4 = U3D_GRID%CELL_NODES(4,I)

         X1 = U3D_GRID%NODE_COORDS(1, V1)
         X2 = U3D_GRID%NODE_COORDS(1, V2)
         X3 = U3D_GRID%NODE_COORDS(1, V3)
         X4 = U3D_GRID%NODE_COORDS(1, V4)
         Y1 = U3D_GRID%NODE_COORDS(2, V1)
         Y2 = U3D_GRID%NODE_COORDS(2, V2)
         Y3 = U3D_GRID%NODE_COORDS(2, V3)
         Y4 = U3D_GRID%NODE_COORDS(2, V4)
         Z1 = U3D_GRID%NODE_COORDS(3, V1)
         Z2 = U3D_GRID%NODE_COORDS(3, V2)
         Z3 = U3D_GRID%NODE_COORDS(3, V3)
         Z4 = U3D_GRID%NODE_COORDS(3, V4)


         ! These are such that PSI_i = SUM_j [ x_j * BASIS_COEFFS(j,i,IC) ] + BASIS_COEFFS(4,i,IC)

         U3D_GRID%BASIS_COEFFS(1,1,I) =  Y2*Z3-Y3*Z2 -Y2*Z4+Y4*Z2 +Y3*Z4-Y4*Z3
         U3D_GRID%BASIS_COEFFS(2,1,I) = -X2*Z3+X3*Z2 +X2*Z4-X4*Z2 -X3*Z4+X4*Z3
         U3D_GRID%BASIS_COEFFS(3,1,I) =  X2*Y3-X3*Y2 -X2*Y4+X4*Y2 +X3*Y4-X4*Y3
         U3D_GRID%BASIS_COEFFS(4,1,I) = -X2*Y3*Z4 +X3*Y2*Z4 +X2*Y4*Z3 -X4*Y2*Z3 -X3*Y4*Z2 +X4*Y3*Z2

         U3D_GRID%BASIS_COEFFS(1,2,I) = -Y1*Z3+Y3*Z1 +Y1*Z4-Y4*Z1 -Y3*Z4+Y4*Z3
         U3D_GRID%BASIS_COEFFS(2,2,I) =  X1*Z3-X3*Z1 -X1*Z4+X4*Z1 +X3*Z4-X4*Z3
         U3D_GRID%BASIS_COEFFS(3,2,I) = -X1*Y3+X3*Y1 +X1*Y4-X4*Y1 -X3*Y4+X4*Y3
         U3D_GRID%BASIS_COEFFS(4,2,I) =  X1*Y3*Z4 -X3*Y1*Z4 -X1*Y4*Z3 +X4*Y1*Z3 +X3*Y4*Z1 -X4*Y3*Z1

         U3D_GRID%BASIS_COEFFS(1,3,I) =  Y1*Z2-Y2*Z1 -Y1*Z4+Y4*Z1 +Y2*Z4-Y4*Z2
         U3D_GRID%BASIS_COEFFS(2,3,I) = -X1*Z2+X2*Z1 +X1*Z4-X4*Z1 -X2*Z4+X4*Z2
         U3D_GRID%BASIS_COEFFS(3,3,I) =  X1*Y2-X2*Y1 -X1*Y4+X4*Y1 +X2*Y4-X4*Y2
         U3D_GRID%BASIS_COEFFS(4,3,I) = -X1*Y2*Z4 +X2*Y1*Z4 +X1*Y4*Z2 -X4*Y1*Z2 -X2*Y4*Z1 +X4*Y2*Z1

         U3D_GRID%BASIS_COEFFS(1,4,I) = -Y1*Z2+Y2*Z1 +Y1*Z3-Y3*Z1 -Y2*Z3+Y3*Z2
         U3D_GRID%BASIS_COEFFS(2,4,I) =  X1*Z2-X2*Z1 -X1*Z3+X3*Z1 +X2*Z3-X3*Z2
         U3D_GRID%BASIS_COEFFS(3,4,I) = -X1*Y2+X2*Y1 +X1*Y3-X3*Y1 -X2*Y3+X3*Y2
         U3D_GRID%BASIS_COEFFS(4,4,I) =  X1*Y2*Z3 -X2*Y1*Z3 -X1*Y3*Z2 +X3*Y1*Z2 +X2*Y3*Z1 -X3*Y2*Z1

         U3D_GRID%BASIS_COEFFS(:,:,I) = -U3D_GRID%BASIS_COEFFS(:,:,I)/6./VOLUME

      END DO





      IF (PROC_ID == 0) THEN
         WRITE(*,*) '==========================================='
         WRITE(*,*) 'Creating boundary grid.'
         WRITE(*,*) '==========================================='
      END IF

      ALLOCATE(U3D_GRID%CELL_FACES_BOUNDARY_INDEX(4,NCELLS))
      U3D_GRID%CELL_FACES_BOUNDARY_INDEX = -1
      ALLOCATE(NODE_ON_BOUNDARY(NNODES))
      NODE_ON_BOUNDARY = .FALSE.
      ALLOCATE(U3D_GRID%NODES_BOUNDARY_INDEX(NNODES))
      U3D_GRID%NODES_BOUNDARY_INDEX = -1
      NUM_BOUNDARY_NODES = 0
      NUM_BOUNDARY_ELEM = 0
      DO I = 1, NCELLS
         DO J = 1, 4
            ! If the face belongs to any physical group, it should be part of the boundary grid
            ! Later, we may want to filter this further
            IF (U3D_GRID%CELL_FACES_PG(J,I) .NE. -1) THEN
               NUM_BOUNDARY_ELEM = NUM_BOUNDARY_ELEM + 1

               V1 = U3D_GRID%CELL_NODES(IND(J,1),I)
               V2 = U3D_GRID%CELL_NODES(IND(J,2),I)
               V3 = U3D_GRID%CELL_NODES(IND(J,3),I)

               IF (.NOT. NODE_ON_BOUNDARY(V1)) THEN
                  NUM_BOUNDARY_NODES = NUM_BOUNDARY_NODES + 1
                  U3D_GRID%NODES_BOUNDARY_INDEX(V1) = NUM_BOUNDARY_NODES
                  NODE_ON_BOUNDARY(V1) = .TRUE.
               END IF
               IF (.NOT. NODE_ON_BOUNDARY(V2)) THEN
                  NUM_BOUNDARY_NODES = NUM_BOUNDARY_NODES + 1
                  U3D_GRID%NODES_BOUNDARY_INDEX(V2) = NUM_BOUNDARY_NODES
                  NODE_ON_BOUNDARY(V2) = .TRUE.
               END IF
               IF (.NOT. NODE_ON_BOUNDARY(V3)) THEN
                  NUM_BOUNDARY_NODES = NUM_BOUNDARY_NODES + 1
                  U3D_GRID%NODES_BOUNDARY_INDEX(V3) = NUM_BOUNDARY_NODES
                  NODE_ON_BOUNDARY(V3) = .TRUE.
               END IF
               
            END IF
         END DO
      END DO      

      U2D_GRID%NUM_CELLS = NUM_BOUNDARY_ELEM
      U2D_GRID%NUM_NODES = NUM_BOUNDARY_NODES
      ALLOCATE(U2D_GRID%CELL_NODES(3, NUM_BOUNDARY_ELEM))
      ALLOCATE(U2D_GRID%CELL_PG(NUM_BOUNDARY_ELEM))
      ALLOCATE(U2D_GRID%NODE_COORDS(3, NUM_BOUNDARY_NODES))

      DO I = 1, NNODES
         IF (NODE_ON_BOUNDARY(I)) THEN
            U2D_GRID%NODE_COORDS(:,U3D_GRID%NODES_BOUNDARY_INDEX(I)) = U3D_GRID%NODE_COORDS(:,I)
         END IF
      END DO

      NUM_BOUNDARY_ELEM = 0

      DO I = 1, NCELLS
         DO J = 1, 4
            IF (U3D_GRID%CELL_FACES_PG(J,I) .NE. -1) THEN
               NUM_BOUNDARY_ELEM = NUM_BOUNDARY_ELEM + 1
               U2D_GRID%CELL_PG(NUM_BOUNDARY_ELEM) = U3D_GRID%CELL_FACES_PG(J,I)
               U3D_GRID%CELL_FACES_BOUNDARY_INDEX(J,I) = NUM_BOUNDARY_ELEM

               V1 = U3D_GRID%CELL_NODES(IND(J,1),I)
               V2 = U3D_GRID%CELL_NODES(IND(J,2),I)
               V3 = U3D_GRID%CELL_NODES(IND(J,3),I)
               U2D_GRID%CELL_NODES(1, NUM_BOUNDARY_ELEM) = U3D_GRID%NODES_BOUNDARY_INDEX(V1)
               U2D_GRID%CELL_NODES(2, NUM_BOUNDARY_ELEM) = U3D_GRID%NODES_BOUNDARY_INDEX(V2)
               U2D_GRID%CELL_NODES(3, NUM_BOUNDARY_ELEM) = U3D_GRID%NODES_BOUNDARY_INDEX(V3)

            END IF
         END DO
      END DO
      
      DEALLOCATE(NODE_ON_BOUNDARY)

      NBOUNDCELLS = NUM_BOUNDARY_ELEM
      NBOUNDNODES = NUM_BOUNDARY_NODES

      ! Compute areas and lengths of boundary mesh
      ALLOCATE(U2D_GRID%CELL_AREAS(U2D_GRID%NUM_CELLS))
      DO I = 1, U2D_GRID%NUM_CELLS
         A = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(1,I))
         B = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(2,I))
         C = U2D_GRID%NODE_COORDS(:, U2D_GRID%CELL_NODES(3,I))

         U2D_GRID%CELL_AREAS(I) = 0.5*NORM2(CROSS(B-A, C-A))
      END DO




      IF (PROC_ID == 0) THEN
         WRITE(*,*) '============================================================='
         WRITE(*,*) 'Done reading grid file.'
         WRITE(*,*) 'It contains ', NNODES, ' nodes and ', NCELLS, ' cells.'
         WRITE(*,*) 'The boundary grid contains ', NBOUNDCELLS, ' surfaces and ', NBOUNDNODES, ' nodes.'
         WRITE(*,*) '============================================================='
      END IF

   END SUBROUTINE READ_3D_UNSTRUCTURED_GRID_MSH


END MODULE grid_and_partition
