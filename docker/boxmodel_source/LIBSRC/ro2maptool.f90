MODULE ro2maptool

! structure for storing RO2 species in the corresponding classes 
  TYPE :: clasro2
    INTEGER             :: nro2          ! number of RO2 in the class
    INTEGER,ALLOCATABLE :: idro2(:)      ! ID of the RO2 species
  END TYPE
  TYPE(clasro2),ALLOCATABLE :: ro2map(:) ! RO2 species map (size: # of pero class)

CONTAINS  
!=======================================================================
! PURPOSE: RO2+RO2 reactions use "counters", i.e. the sum of the 
! concentrations of RO2 species belonging to a given "class". This 
! routine read the files providing the list of species in each class,
! allocate the required memory for each class and store the species 
! that belong to each class.
!=======================================================================
SUBROUTINE readro2(lout,chrsp,ncpero)
  USE boxtool, ONLY: stoperr
  IMPLICIT NONE

  INTEGER, INTENT(IN)          :: lout      ! file unit for information outputs
  CHARACTER(LEN=*), INTENT(IN) :: chrsp(:)  ! list (names) of the species
  INTEGER, INTENT(IN)          :: ncpero    ! # of peroxy (RO2) classes

  CHARACTER(LEN=LEN(chrsp(1))),ALLOCATABLE :: ro2sp(:)
  CHARACTER(LEN=LEN(chrsp(1))) :: liner
  CHARACTER(LEN=20)            :: filnam
  INTEGER :: i,j,nclas,nspe,ispe,ios,iloc
  INTEGER :: flg(SIZE(chrsp))
  LOGICAL :: loerr
  
  INTEGER,PARAMETER :: filu=12                          ! unit of the files to be read 
  CHARACTER(LEN=10),PARAMETER :: progname='readro2'
  CHARACTER(LEN=80) :: mesg1, mesg2

  flg(:) = 0

! Minor change required in this routine if ncpero > 9 (see open below)
  IF (ncpero  > 9) THEN
    mesg1="# of RO2 classes > 9."
    mesg2="Change the program to read more classes." 
    CALL stoperr(progname,mesg1,mesg2)
  ENDIF

  WRITE(lout,*) ' '
  WRITE(lout,*) '---- RO2+RO2 reactions -----'

! LOOP OVER THE RO2 CLASS
! -----------------------
  ro2loop: &
  DO nclas=1,ncpero
    PRINT'("   read RO2 counter file:",i2)',nclas

! open file to be read
! --------------------
    filnam(1:5)='indat'
    WRITE(filnam(6:6),'(i1)') nclas
    filnam(7:)='.ro2 '
    WRITE(lout,*) '  - open: ',TRIM(filnam)
    OPEN(filu, FILE=filnam, STATUS='OLD', IOSTAT=ios)
    IF (ios /= 0) THEN
      mesg1="file not found/open: "//TRIM(filnam)
      CALL stoperr(progname,mesg1)
    ENDIF

! read RO2 species that belong to the class
! -----------------------------------------
    READ(filu,*,IOSTAT=ios) nspe
    IF (ios/=0) THEN 
      mesg1="cannot read # species in file: "//TRIM(filnam)
      CALL stoperr(progname,mesg1)
    ENDIF

    ! allocate required memory to the ro2map
    ro2map(nclas)%nro2 = nspe
    ALLOCATE(ro2map(nclas)%idro2(nspe+1))  ! add 1 slot (just in case)
    ro2map(nclas)%idro2(:)=0
    
    ALLOCATE(ro2sp(nspe+1))     ! not used out of this subroutine
     
    ! loop over the species in the file
    ispe=0
    speloop: & 
    DO                  
      READ(filu,'(a)',IOSTAT=ios) liner
      IF (ios/=0) THEN 
        mesg1="nothing left to read in file: "//TRIM(filnam)
        mesg2="keyword 'END' missing?" 
        CALL stoperr(progname,mesg1,mesg2)
      ENDIF

      ! escape the loop if keyword 'END' is found
      IF (liner(1:3)=='END') THEN   
        IF (ispe/=nspe) THEN
          mesg1="number of species mismatches expected record in: "//TRIM(filnam)
          CALL stoperr(progname,mesg1)
        ENDIF
        EXIT speloop
      ENDIF 

      ! check memory (just in case)
      ispe=ispe+1
      IF (ispe > nspe) THEN
        mesg1="# of RO2 species > expected record in: "//TRIM(filnam)
        CALL stoperr(progname,mesg1)
      ENDIF

      ! check empty line
      IF (LEN_TRIM(liner) == 0) THEN
        mesg1="Unexpected empty line read in: "//TRIM(filnam)
        CALL stoperr(progname,mesg1)
      ENDIF
      
      ! store (temporarily) the species 
      ro2sp(ispe)=liner
    ENDDO speloop
    CLOSE(filu)
  
! get the ID for the RO2 species & store ID in idchemro2
! ------------------------------------------------------
    loerr=.FALSE. ;  iloc=1

! fast seek - require same order in tables (OK in gecko)    
    seekloop1:&    
    DO i=1,ro2map(nclas)%nro2
      DO j=iloc,SIZE(chrsp)
        IF (ro2sp(i)==chrsp(j)) THEN

          ! check for duplicate
          IF (flg(j)/= 0) THEN          
            mesg1="The following species was provided twice: "//TRIM(ro2sp(i))
            mesg2="last found in file: "//TRIM(filnam)
            CALL stoperr(progname,mesg1,mesg2)
          ENDIF

          ! store ID
          ro2map(nclas)%idro2(i)=j
          flg(j)=1 ; iloc=j+1 
          CYCLE seekloop1
        ENDIF
      ENDDO
      
      IF (ro2sp(i)=='GCH3O2') CYCLE ! /!\ patch for 1 RO2 class
      
      ! if that point is reached then species not found => try unsorted loop
      loerr=.TRUE. 
      EXIT seekloop1
    ENDDO seekloop1

! if fast seek did not succeed, then try slow seek (no order required - OK for MCM)
    IF (loerr) THEN  
      PRINT'("   => try slow seek for:",i2)',nclas
      flg(:) = 0

      seekloop2:&  
      DO i=1,ro2map(nclas)%nro2
        DO j=1,SIZE(chrsp)
          IF (ro2sp(i)==chrsp(j)) THEN

            ! check for duplicate
            IF (flg(j)/= 0) THEN
              mesg1="The following species was provided twice: "//TRIM(ro2sp(i))
              mesg2="last found in file: "//TRIM(filnam)
              CALL stoperr(progname,mesg1,mesg2)
            ENDIF

            ! store ID             
            ro2map(nclas)%idro2(i)=j
            flg(j) = 1
            CYCLE seekloop2
          ENDIF
        ENDDO

        ! if that point is reached then species not found
        mesg1="The following species not identified: "//TRIM(ro2sp(i))
        mesg2="but provided in file: "//TRIM(filnam)
        CALL stoperr(progname,mesg1,mesg2)
      ENDDO seekloop2
    ENDIF
    DEALLOCATE(ro2sp)
  ENDDO ro2loop
      
END SUBROUTINE readro2

END MODULE ro2maptool
