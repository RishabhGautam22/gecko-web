MODULE initnbound
IMPLICIT NONE
REAL    :: temp0
REAL    :: rh0
INTEGER :: iyear     ! date - year
INTEGER :: imonth    ! date - month
INTEGER :: iday      ! date - day
REAL    :: sla       ! latitude,
REAL    :: slo       ! longitude
REAL    :: tz        ! time zone

CONTAINS
!=======================================================================  
! Purpose: return T, RH, water for the time provided as input. User 
! can set prescribed function, as needed.
!=======================================================================  
SUBROUTINE getconstrain(time,temp,sumc,rh,water)
  USE parameter_mod, ONLY: avogadro,Rgas
  IMPLICIT NONE

  REAL, INTENT(IN)    :: time        ! time (warning: modulo 1 day)
  REAL, INTENT(OUT)   :: temp        ! temperature
  REAL, INTENT(OUT)   :: sumc        ! M (in molec/cm3)
  REAL, INTENT(OUT)   :: rh          ! relative humidity
  REAL, INTENT(OUT)   :: water       ! water concentration (molec/cm3)
  
! "Magnus formula": psat(H2O) = c1 * EXP(c2*T/(c3+T)) (units Pa)
!  where: T is expressed in degrees ! : T(K)-273.16
!         c1 = 610.94, c2 = 17.625, c3 = 243.04
! Ref: Alduchov & Eskridge 1996, J.Appl.Meteorol. 35 601-609.
! quoted by: Lawrence, 2005, BAMS DOI:10.1175/BAMS-86-2-225
  REAL,PARAMETER :: c1 = 610.94, c2 = 17.625,  c3 = 243.04
  REAL           :: TC, psat_H2O, p_H2O

! set temperature
! --------------------
  ! As provided in the input file:
  temp=temp0

! set pressure (sumc in molec/cm3) 
! --------------------------------
  ! Use perfect gas law :sumc=(pres*6.022E+23)/(8.32*temp)
  ! Constant value (forced by user):
  sumc=2.5E19

! set relative humidity (in %) and the water content 
! -------------------------------------
  ! As provided in the input file:
  rh=rh0

! Water concentration (molec/cm3)
! -------------------------------
  TC = temp-273.16                   ! temperature in Celsius
  psat_H2O = c1*EXP(c2*TC/(c3+TC))   ! saturation water pressure
  p_H2O = psat_H2O*rh/100.           ! water pressure
  water = p_H2O / (Rgas*temp*1.e6/avogadro)

END SUBROUTINE getconstrain

!=======================================================================  
! PURPOSE: Read the initial concentrations of some species
!=======================================================================  
SUBROUTINE readinitial(lout,numsp,chrsp,cbox,nfix,idfix,cfix,noxfix,cnox)
  USE stringdata, ONLY: stringarg
  USE boxtool, ONLY: spe2id,stoperr
  USE parameter_mod, ONLY: smallc
  IMPLICIT NONE

  INTEGER,INTENT(IN)  :: lout              ! file unit (report)
  INTEGER,INTENT(IN)  :: numsp             ! # of species
  CHARACTER(LEN=*),INTENT(IN) :: chrsp(:)  ! name of the species
  REAL,INTENT(OUT)    :: cbox(:)           ! species concentration
  INTEGER,INTENT(OUT) :: nfix              ! number of species having fixed concentration
  INTEGER,INTENT(OUT) :: idfix(:)          ! id of the species having fixed concentration
  REAL,INTENT(OUT)    :: cfix(:)           ! value of the fixed concentration
  INTEGER,INTENT(OUT) :: noxfix            ! flag to fix NOx (NO+O2) concentration 
  REAL,INTENT(OUT)    :: cnox              ! value of the fixed NOx concentration

  CHARACTER(LEN=LEN(chrsp(1))) :: species
  LOGICAL            :: lostop
  CHARACTER(LEN=76)  :: line
  CHARACTER(LEN=4)   :: keyword
  INTEGER            :: i,isp,filu,ios

! parameter for subroutine stringarg (see stringdata module)   
  INTEGER, PARAMETER :: mxarg=10      ! maximum # of arguments in a line
  INTEGER :: posi(mxarg,2)
  INTEGER :: narg

  CHARACTER(LEN=15),PARAMETER :: progname='readinitial'
  CHARACTER(LEN=80) :: mesg1, mesg2

! -----------------------------
! INITIALIZE               
! -----------------------------
  lostop=.FALSE.

! set default value (likely overwritten next)
  cbox(:) = smallc                    ! lower limit for a concentration (avoid 0.) 
  nfix=0  ; idfix(:)=0 ; cfix(:)=0.
  noxfix=0 ; cnox=0.
  
! -----------------------------
! OPEN THE FILE 
! -----------------------------
  filu=12
  OPEN(filu,FILE='indat.key',FORM='formatted',STATUS='old',IOSTAT=ios)
  IF (ios /= 0) THEN
    mesg1="file indat.key not found"
    CALL stoperr(progname,mesg1)
  ENDIF

  WRITE(lout,*)
  WRITE(lout,*)'--------- Keyword input -----------'

! -----------------------------
! READ NEXT LINE
! -----------------------------
  rdline: & 
  DO                 ! escape the loop if keyword 'END ' is found 
    READ(filu,'(a4,(a))',IOSTAT=ios) keyword, line
    IF (ios/=0) THEN 
      mesg1="End of file reached - keyword 'END' missing ?"
      CALL stoperr(progname,mesg1)
    ENDIF

! comment line - read next
    IF (keyword(1:1)=='/' .OR. keyword(1:1)=='!') CYCLE rdline

! CASE SELECTOR FOR KEYWORDS
! -----------------------------
    SELECT CASE (keyword)   

     ! End of file - escape
     CASE ('END ') 
       EXIT rdline

     ! Initial concentrations
     CASE('REAC')
       CALL stringarg(line,narg,posi)
       IF (narg /= 2) THEN            ! species + c*inbox
         WRITE(lout,*) " --error--, in readinitial while reading reac."
         WRITE(lout,*) "            # of arg is not appropriate:"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       species=line(posi(1,1):posi(1,2))
       isp=spe2id(lout,chrsp,numsp,species,csub='readinitial.f90',loerr=.FALSE.)
       IF (isp == 0) THEN
         WRITE(lout,*) " --error--, in readinitial while reading reac."
         WRITE(lout,*) "            species is unkown in line :"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       READ(line(posi(2,1):),*,IOSTAT=ios) cbox(isp)
       IF (ios/=0) THEN
         WRITE(lout,*) "--error--, in readinitial while reading reac:"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       CYCLE rdline

     ! fixed concentrations
     CASE('CFIX')
       nfix = nfix+1
       IF (nfix > SIZE(idfix)) THEN
         WRITE(lout,*) " --error--, in readinitial while reading cfix."
         WRITE(lout,*) "            # of species exceed max allowed: ",SIZE(idfix)
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF

       CALL stringarg(line,narg,posi)
       IF (narg /= 2) THEN            
         WRITE(lout,*) " --error--, in readinitial while reading cfix."
         WRITE(lout,*) "            # of arg is not appropriate:"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       species=line(posi(1,1):posi(1,2))
       isp=spe2id(lout,chrsp,numsp,species,csub='readinitial.f90',loerr=.FALSE.)
       IF (isp == 0) THEN
         WRITE(lout,*) " --error--, in readinitial while reading reac."
         WRITE(lout,*) "            species is unkown in line :"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       idfix(nfix)=isp
       READ(line(posi(2,1):),*,IOSTAT=ios) cfix(nfix)
       IF (ios/=0) THEN
         WRITE(lout,*) "--error--, in readinitial while reading reac:"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
         CYCLE rdline
       ENDIF
       CYCLE rdline

     ! Fixed NOx concentration
     CASE('CNOX')
       READ(line,*,IOSTAT=ios) cnox
       IF (ios/=0) THEN
         WRITE(lout,*) "--error--, in readinitial while reading CNOX:"
         WRITE(lout,*) keyword,' ',TRIM(line) 
         lostop=.TRUE.
       ENDIF
       noxfix=1                       ! raise the "NOx fixed" flag
       CYCLE rdline

    ! Unidentified keyword
    CASE DEFAULT
      WRITE(lout,*)'--error--, in readinitial. Keyword unknown:'
      WRITE(lout,*) keyword,' ',TRIM(line) 
      lostop=.TRUE.
      CYCLE rdline

    END SELECT  
  ENDDO rdline
  CLOSE(filu)

! Check inconsistencies in inputs
! -------------------------------

  ! if NOx is fixed then NO and NO2 must be "free" 
  IF (noxfix/=0) THEN
    IF (nfix/=0) THEN
      DO i=1,nfix
        IF (chrsp(idfix(i))=='GNO ') THEN
          WRITE(lout,*) "--error--, in readinitial. Inconsistency in"
          WRITE(lout,*) " input data: both NO and NOx are fixed."
          lostop=.TRUE.
        ENDIF
        IF (chrsp(idfix(i))=='GNO2 ') THEN
          WRITE(lout,*) "--error--, in readinitial. Inconsistency in"
          WRITE(lout,*) " input data: both NO2 and NOx are fixed."
          lostop=.TRUE.
        ENDIF
      ENDDO
    ENDIF
  ENDIF

! -------------------
! STOP if error found
! -------------------
  IF (lostop) THEN
    mesg1="Miscellaneous errors found while reading indat.key (see *.out file)"
    mesg2="See details in the *.out file"
    CALL stoperr(progname,mesg1,mesg2)
  ENDIF

END SUBROUTINE readinitial

!=======================================================================  
! PURPOSE: Read the initial conditions (temperature, initial C ...) 
! for the simulation and parameters for the time integration (start and 
! stop time, # of time steps, solver tolerance ...) in the namelist
!=======================================================================  
SUBROUTINE rd_nml
  USE parameter_mod, ONLY: flag_ZA,fixedZA,tstart,tstop,ntstep,nskip,&
                           rtol,atol,numit,dtmin,&
                           nvoc,nvic,Rp0,denssoa,densseed,Mseed
  USE mechdata,      ONLY: kdilu, kaerloss
  IMPLICIT NONE
  
  CHARACTER(LEN=9), PARAMETER :: filein='simu.nml'
  INTEGER,PARAMETER           :: nmlu=45
  INTEGER                     :: ierr
  
  NAMELIST /simu_conditions/ tstart,tstop,ntstep,nskip,rtol,atol,numit,dtmin
  NAMELIST /env_conditions/  temp0,rh0,flag_ZA,fixedZA, kdilu, kaerloss,&
                             iday,imonth,iyear,sla,slo,tz
  NAMELIST /mass_transfert_parameters/ Rp0,nvoc,nvic,denssoa,densseed,Mseed

  OPEN(nmlu, file=filein, iostat=ierr)
    IF (ierr/=0) then
      PRINT *,'--error--, problem reading namelist file'
      STOP
    ENDIF

    READ(nmlu,nml=simu_conditions, iostat=ierr)
    IF (ierr/=0) then
      PRINT *,'--error--, problem reading "simu_conditions" arguments in namelist'
      STOP
    ENDIF

    READ(nmlu,nml=env_conditions, iostat=ierr)
    IF (ierr/=0) then
      PRINT *,'--error--, problem reading "env_conditions" arguments in namelist'
      STOP
    ENDIF

    READ(nmlu,nml=mass_transfert_parameters, iostat=ierr)
    IF (ierr/=0) then
      PRINT *,'--error--, problem reading "env_conditions" arguments in namelist'
      STOP
    ENDIF
  CLOSE(nmlu)

  IF ( (imonth <1 .OR. imonth > 12) .OR. &
       (iday  < 1 .OR. iday   > 31)       ) THEN 
    PRINT*, 'problem with the date: ',iday,'/',imonth,'/',iyear
    STOP
  ENDIF
  
END SUBROUTINE rd_nml

!-----------------------------------
! default  parameters for bomxmodel simulations
! => values to use if not supplied by user
!-----------------------------------
SUBROUTINE define_defaults
  USE parameter_mod, ONLY: tstart,tstop,ntstep,nskip,rtol,atol,numit,dtmin,&
                           Rp0,nvoc,nvic,denssoa,densseed,Mseed
  USE mechdata,      ONLY: kdilu, kaerloss
  IMPLICIT NONE

! simulation parameters
    tstart = 0. ; tstop  = 3600.
    ntstep = 12 ; nskip  = 1         
    rtol=0.01   ; atol=100. 
    numit=10    ; dtmin=0.1

    iyear=2000 ; imonth=3 ; iday=21       
    sla=45.    ; slo=0.   ; tz=0.     

    rh0      = 3.
    temp0    = 298.
    kdilu    = 0.      ! dilution rate constant (apply to all species in all phase)
    kaerloss = 0.      ! particle loss rate constant (apply to species in particles only) 

    nvoc     = 0.      ! Non Volatile Organic Compound (concentration)
    nvic     = 0.      ! Non Volatile Inorganic Compound (concentration)
    Rp0      = 1E-5    ! Initial radius of the seed particles 
    denssoa  = 1.06    ! soa density (g.cm-3)    
    densseed = 1.06    ! seed density (g.cm-3)    
    Mseed    = 427.    ! seed molar mass (DOS=427 g/mol)
    
END SUBROUTINE define_defaults

END MODULE initnbound
