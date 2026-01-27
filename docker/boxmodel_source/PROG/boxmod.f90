!-----------------------------------------------------------------------------
! BOX MODEL: 
!-----------------------------------------------------------------------------
PROGRAM boxmod
  USE boxtool, ONLY: getnrec,stoperr
  USE vaporpressure, ONLY: readpvap,soapartition,mtrat,soapartition
  USE transfertool, ONLY: readtransfer,rev2forward,particleloss,get1stpnc
  USE viscositytool, ONLY: readtg,ntg,tg,alpha
  USE diffusivitytool, ONLY: readdv,chge_dvol,ndv,dvol
  USE parameter_mod
  USE lpmaptool, ONLY: lpmap,make_lpmap
  USE phototool, ONLY: mxsza,njlab,jmap,make_jmap,cscj,ratj
  USE ro2maptool, ONLY: readro2,ro2map
  USE solver, ONLY: twostep
  USE mechdata, ONLY: lenspe,numsp,numre,numself, &
                      rxrgnt,rxpdct,idselfreac,kdilu,kaerloss 
  USE frozenstuff, ONLY: resetc,resetnox,nfix,idfix,cfix,noxfix,cnox,&
                         storeid,idch3o2
  USE readchem, ONLY: readmech
  USE phototool, ONLY: readjvalue, getjvalue
  USE initnbound
  USE ratetool, ONLY: ratecoef
  USE ncdf_saving
  IMPLICIT NONE
  
! species concentration 
  REAL,ALLOCATABLE  :: cbox(:)             ! concentration of species i

! variable linked to the chemical scheme
  CHARACTER(LEN=mxlsp),ALLOCATABLE :: chrsp(:)  ! list (names) of the species
  INTEGER :: mxleft                        ! max # of reagents in a reaction
  INTEGER :: mxright                       ! max # of products in a reaction
  INTEGER :: mxaux                         ! max # of auxiliary (extra) data in a rx
  INTEGER :: num_n                         ! # of "simple" thermal rx
  INTEGER :: num_m                         ! # of 'TBODY' rx
  INTEGER :: numfo                         ! # of 'FALLOFF' rx
  INTEGER :: numhv                         ! # of 'HV' rx
  INTEGER :: numextra                      ! # of 'EXTRA' rx
  INTEGER :: numo2                         ! # of 'OXYGEN' rx
  INTEGER :: numiso                        ! # of 'ISOM' rx
  INTEGER,ALLOCATABLE :: id_m(:)           ! mechanism rx index of the ith 'TBODY' reaction
  INTEGER,ALLOCATABLE :: idfo(:)           ! mechanism rx index of the ith 'FALLOFF' reaction
  INTEGER,ALLOCATABLE :: idhv(:)           ! mechanism rx index of the ith 'HV' reaction
  INTEGER,ALLOCATABLE :: idextra(:)        ! mechanism rx index of the ith 'EXTRA' reaction
  INTEGER,ALLOCATABLE :: ido2(:)           ! mechanism rx index of the ith 'OXYGEN' reaction
  INTEGER,ALLOCATABLE :: idiso(:)          ! mechanism rx index of the ith 'ISOM' reaction
  INTEGER,ALLOCATABLE :: rxhvtag(:)        ! tag (label) number for J(xs,qy) of ith hv reaction
                  
  INTEGER :: nummeo2                       ! # of 'MEPERO' (CH3O2) rx 
  INTEGER :: ncpero                        ! # of peroxy classes (for RO2+RO2 reactions)
  INTEGER,ALLOCATABLE :: idmeo2(:)         ! mechanism rx index of the ith 'MEPERO' reaction
  INTEGER,ALLOCATABLE :: numreacro2(:)     ! # of rx with "PEROj" in class [j]
  INTEGER,ALLOCATABLE :: idreacro2(:,:)    ! [i,j] mechanism rx index of the ith 'PEROj' reaction
  REAL,ALLOCATABLE    :: cro2(:)           ! RO2 concentration for 'PERO' in class [j]  
  REAL    :: cmeo2                         ! CH3O2 concentration

  REAL,ALLOCATABLE :: arrhcf(:,:)          ! [i,3] arrhenius parameter for the ith mech. rx index
  REAL,ALLOCATABLE :: focf(:,:)            ! [i,j] j coef. (k0, Fc...) for the ith falloff rx
  REAL,ALLOCATABLE :: hvfact(:)            ! scaling factor applied to reference J
  REAL,ALLOCATABLE :: extracf(:,:)         ! [i,j] j coef. (k...) for the ith "extra" rx)
  REAL,ALLOCATABLE :: isocf(:,:)           ! [i,j] j extra coef. for the ith 'ISOM' rx
  REAL,ALLOCATABLE :: wmol(:)              ! molar mass of the ith species
  REAL,ALLOCATABLE :: qfor(:)              ! rate coefficient for the ith reaction

! index for starting/ending of phase contributors in the list of species
  INTEGER :: nspg                          ! # of gas phase species
  INTEGER :: nspp                          ! # of part. phase species
  INTEGER :: nspw                          ! # of wall phase species
  INTEGER :: ptrgas(2)                     ! pointer (1st and last) gas phase species
  INTEGER :: ptrpart(2)                    ! pointer (1st and last) part. phase species
  INTEGER :: ptrwall(2)                    ! pointer (1st and last) wall phase species
                  
! mass transfer parameters                      
  INTEGER :: numain                        ! # of gas -> part. rxn
  INTEGER :: numaou                        ! # of part. -> gas rxn
  INTEGER :: numwin                        ! # of gas -> wall rxn
  INTEGER :: numwou                        ! # of wall -> gas rxn
  REAL    :: pnc                           
  INTEGER,ALLOCATABLE :: idain(:)          ! mechanism rx index for ith g->p (AIN)
  INTEGER,ALLOCATABLE :: idaou(:)          ! mechanism rx index for ith p->g (AOU)
  INTEGER,ALLOCATABLE :: idwin(:)          ! mechanism rx index for ith g->w (WIN)
  INTEGER,ALLOCATABLE :: idwou(:)          ! mechanism rx index for ith w->g (WOU)
  REAL,ALLOCATABLE    :: woucf(:,:)        ! [i,j] j extra coef. for the ith w->g rx (WOU)
  REAL,ALLOCATABLE    :: wincf(:,:)        ! [i,j] j extra coef. for the ith g->w rx (WIN)
  INTEGER,ALLOCATABLE :: idgsat(:)         ! ID of partitioning species - gas phase
  INTEGER,ALLOCATABLE :: idasat(:)         ! ID of partitioning species - part. phase
  INTEGER,ALLOCATABLE :: idwsat(:)         ! ID of partitioning species - wall phase
  INTEGER,ALLOCATABLE :: gairefor(:)       ! rxn ID of forward g->p matching the reverse p->g
  INTEGER,ALLOCATABLE :: gwirefor(:)       ! rxn ID of forward g->w matching the reverse w->g

! SOA thermodynamic
  CHARACTER(LEN=mxlsp),ALLOCATABLE :: namsat(:)  ! name of the partitioning species (with 'G')
  INTEGER             :: soa_fg            ! flag to solve SOA (0: no SOA, 1: equilib. , 2: dynamic)
  INTEGER             :: nsat              ! # of partitioning species
  REAL,ALLOCATABLE    :: caer(:)           ! particle phase concentration (molec/cm3)
  REAL                :: maer              ! SOA mass concentration
  REAL                :: rt4psat           ! reference T for vapor pressure
  REAL,ALLOCATABLE    :: rpsat(:)          ! vapor pressure at reference T
  REAL,ALLOCATABLE    :: eheat(:)          ! latent heat of evaporation (@ ref T in K-1)
  INTEGER,ALLOCATABLE :: mech2sat(:)       ! [j] index of mech species j in the pvap table

! parameter for the simulation
  REAL    :: time                          ! current time
  REAL    :: deltat                        ! time step (to 'exit' the solver)
  REAL    :: tout                          ! time out (to 'exit' the solver)
  REAL    :: dtmax                         ! maximum delta time (solver)
  REAL    :: timemod                       ! time - modulo 1 day
  REAL    :: temp                          ! temperature
  REAL    :: sumc                          ! M (in molec/cm3)
  REAL    :: rh                            ! relative humidity
  REAL    :: water                         ! water concentration (molec/cm3)

! ouput file numbers
  INTEGER,PARAMETER :: lout=11             ! file unit for information outputs 
  INTEGER,PARAMETER :: lmaer=14            ! output file (formatted) - SOA mass  
  INTEGER,PARAMETER :: lro2=15             ! output file (formatted) - RO2 concentration
  
  INTEGER :: i,k,iskip
  INTEGER :: idum1

  CHARACTER(LEN=15),PARAMETER :: progname='main_box'
  CHARACTER(LEN=80) :: mesg1


! -------------------------------------------
! open log/info file and initialize
! -------------------------------------------
  OPEN(lout,FILE='outdat.out',FORM='formatted')   ! open log/info file
  soa_fg = 1                                      ! assume SOA is computed (overwritten next as needed)

! -------------------------------------------
! READ THE SIZE OF TABLES AND ALLOCATE MEMORY
! -------------------------------------------
  OPEN(12,FILE='indat.li', FORM='unformatted', STATUS='old')
  READ(12) lenspe,mxleft,mxright,mxaux
  READ(12) numsp,numre,num_n,numhv,numo2,num_m,numfo,numextra, &
           nummeo2,idum1,numain,numaou,numwin,numwou,numiso,&
           numself
  READ(12) ncpero

  ALLOCATE(numreacro2(ncpero)) ; numreacro2(:)=0      ! # of rx with RO2 in class [i]
  READ(12) (numreacro2(i),i=1,ncpero)
  CLOSE(12)
  
  ! Allocate memory for the "chemistry" tables
  ALLOCATE(chrsp(numsp)) ;  chrsp(:)=' '              ! list (names) of the species
  ALLOCATE(cbox(numsp))  ;  cbox(:)=smallc            ! concentration of species i
  ALLOCATE(qfor(numre))                               ! rate coefficient for the ith reaction
  ALLOCATE(rxpdct(numre))                             ! product list of rx[i]
  ALLOCATE(rxrgnt(numre))                             ! reagent list of rx[i]
  ALLOCATE(idselfreac(numself,2))                     ! [i,1] mech. rx index and species ID [i,2] for self
  ALLOCATE(arrhcf(numre,3)) ; arrhcf(:,:)=0.          ! [i,3] arrhenius coef. for the ith rx
  ALLOCATE(wmol(numsp)) ; wmol(:)=0.                  ! molar mass of the ith species
  ALLOCATE(id_m(num_m)) ; id_m(:)=0                   ! mechanism rx index of the ith 'TBODY' rxn
  ALLOCATE(idfo(numfo)) ; idfo(:)=0                   ! mechanism rx index of the ith 'FALLOFF' rxn
  ALLOCATE(focf(numfo,mxaux)) ; focf(:,:)=0.          ! [i,j] j coef. for the ith falloff rxn
  ALLOCATE(idhv(numhv)) ; idhv(:)=0                   ! mechanism rx index of the ith 'HV' reaction
  ALLOCATE(rxhvtag(numhv)) ; rxhvtag(:)=0             ! tag (label) number for J(xs,qy) of ith hv reaction
  ALLOCATE(hvfact(numhv)) ; hvfact(:)=0.              ! scaling factor applied to reference J
  ALLOCATE(idextra(numextra)) ; idextra(:)=0          ! mechanism rx index of the ith 'EXTRA' reaction
  ALLOCATE(extracf(numextra,mxaux)) ; extracf(:,:)=0. ! [i,j] j coef. (k...) for the ith "extra" rx
  ALLOCATE(ido2(numo2)) ; ido2(:)=0                   ! mechanism rx index of the ith 'OXYGEN' reaction
  ALLOCATE(idiso(numiso)) ; idiso(:)=0                ! mechanism rx index of the ith 'ISOM' reaction
  ALLOCATE(isocf(numiso,mxaux)) ; isocf(:,:)=0.       ! [i,j] j extra coef. for the ith 'ISOM' rx

  ALLOCATE(idmeo2(nummeo2)) ; idmeo2(:)=0                  ! (mxrpero)
  ALLOCATE(idreacro2(MAXVAL(numreacro2),ncpero))  ; idreacro2(:,:)=0   ! (mxrpero,mxcpero)

  ALLOCATE(ro2map(ncpero))
  ALLOCATE(cro2(ncpero)) ; cro2(:)=0.

  ALLOCATE(idain(numain)) ; idain(:)=0                ! ID # for g->p transfert
  ALLOCATE(idaou(numaou)) ; idaou(:)=0                ! ID # for p->g rxn (mxt)
  ALLOCATE(idwin(numwin)) ; idwin(:)=0                ! ID # for g->w rxn (mxt)
  ALLOCATE(idwou(numwou)) ; idwou(:)=0                ! ID # for w->g rxn (mxt)
  ALLOCATE(woucf(numwou,1)) ; woucf(:,:)=0.           ! extra coef. w->g rxn (mxt,1)
  ALLOCATE(wincf(numwin,1)) ; wincf(:,:)=0.           ! extra coef. g->w rxn (mxt,1)

  ALLOCATE(lpmap(numsp))                              ! loss_production map (see lpmaptool)

  ! Allocate memory for the "thermodynamic" tables and condensed phases
  IF (numain > 0 .OR. numwin > 0) soa_fg=2
  IF (soa_fg > 0) THEN
    nsat=getnrec(lout,'pvap.sat')
    PRINT*, 'nsat=',nsat
    ALLOCATE(caer(nsat)) ; caer(:)=0.               ! particle phase concentration
    ALLOCATE(rpsat(nsat)) ; rpsat(:)=0.             ! vapor pressure at reference T
    ALLOCATE(eheat(nsat)) ; eheat(:)=0.             ! latent heat of evaporation (@ ref T in K-1)
    ALLOCATE(namsat(nsat)) ; namsat(:)=' '          ! name of the partitioning species (with 'G')
    ALLOCATE(idgsat(nsat)) ; idgsat(:)=0            ! ID of partitioning species - gas phase 
    ALLOCATE(idasat(nsat)) ; idasat(:)=0            ! ID of partitioning species - part. phase
    ALLOCATE(idwsat(nsat)) ; idwsat(:)=0            ! ID of partitioning species - wall phase
    ALLOCATE(gairefor(nsat)) ; gairefor(:)=0        ! rxn ID of forward g->p matching the reverse p->g
    ALLOCATE(gwirefor(nsat)) ; gwirefor(:)=0        ! rxn ID of forward g->w matching the reverse w->g
    ALLOCATE(mech2sat(numsp)) ; mech2sat(:)=0       ! [j] index of mech species j in the pvap table

    ntg=getnrec(lout,'Tg.dat')
    IF (ntg/=nsat) THEN
      mesg1="# record in Tg and Pvap are not identical. "
      CALL stoperr(progname,mesg1)
    ENDIF
    PRINT*, '# of Tg=',ntg
    ALLOCATE(Tg(ntg)) ; Tg(:)=0.
    ALLOCATE(alpha(ntg)) ; alpha(:)=0.

    IF (soa_fg==2) THEN 
      ndv=getnrec(lout,'mdv.dat')
      IF (ndv/=nsat) THEN
        mesg1="# record in mdv and Pvap are not identical. "
        CALL stoperr(progname,mesg1)
      ENDIF
      ALLOCATE(dvol(ndv)) ; dvol(:)=0.
    ENDIF

  ENDIF

! -------------------------------------------
! READ CHEMICAL MECHANISM DATA
! -------------------------------------------
  WRITE(6,*) 'reading the chemical scheme ...'
  CALL readmech(lout,numreacro2, &
                nspg,nspp,nspw,ptrgas,ptrpart,ptrwall, &
                idhv,ido2,id_m,idfo,idextra,idmeo2,idreacro2,&
                idain,idaou,idwin,idwou,idiso,&
                arrhcf,rxhvtag,hvfact,focf,extracf,&
                woucf,wincf,isocf,&
                idselfreac,wmol,chrsp,njlab)
  
  ! Allocate memory for J tables and set values
  PRINT'("number of distinct J label identified: ",i4)',njlab
  ALLOCATE(jmap(njlab))                ! map of all J labels
  ALLOCATE(ratj(njlab,mxsza))          ! values of tabulated J values
  ALLOCATE(cscj(njlab,mxsza*3))        ! cubic spline coef. for J values
  CALL make_jmap(numhv,rxhvtag)

  ! Make lpmap structure for indexing prod/loss rxn of each species
  PRINT'("starting index map for chemical loss/production")'
  CALL make_lpmap(numsp,numre)

  ! read data for RO2 counting species (fill ro2map)
  PRINT'("reading the RO2 database")'
  CALL readro2(lout,chrsp,ncpero)

! -------------------------------------------
! READ DATA FOR PHASE PARTITIONING 
! -------------------------------------------
  IF (soa_fg > 0) THEN
    PRINT*, 'reading vapor pressure and Tg data'
    CALL readpvap(chrsp,numsp,nsat,namsat,idgsat,rt4psat,rpsat,eheat,mech2sat)
    CALL readTg(chrsp,numsp,namsat,idgsat)
    IF (soa_fg==2) CALL readdv(chrsp,numsp,namsat,idgsat,idasat,wmol)
  ENDIF 

  ! read ID for gas, particles, wall for a dynamical mass transport
  IF (soa_fg == 2) THEN
    PRINT*, 'calling transfer stuff ...'
    CALL readtransfer(chrsp,numsp,numain,numwin,&
                      nsat,namsat,idasat,idwsat,mech2sat)
    CALL rev2forward(lout,chrsp,numaou,numwin,numwou,nsat,namsat, &
                     idaou,idwin,idwou,gairefor,gwirefor)
    ! pre-calculate parameters required to compute diffusivity
    CALL chge_dvol(idasat,wmol)
  ENDIF
  
! -------------------------------------------
! SET SIMULATION DATA
! -------------------------------------------
  ! define default values
  CALL define_defaults()
  ! read namelist
  CALL rd_nml()
  ! read key file for precursors concentration
  CALL readinitial(lout,numsp,chrsp,cbox,nfix,idfix,cfix,noxfix,cnox)

  ! set particle number concentration (at initial time)         
  CALL get1stpnc(lout,pnc) 

  ! read photolyses frequencies
  PRINT*, 'reading the photolysis data ...'
  IF (numhv > 0)  CALL readjvalue(lout)

  ! store some species ID and reset frozen concentration (if any)
  CALL storeid(lout,numsp,chrsp)
  IF (nfix/=0) CALL resetc(cbox)
  IF (noxfix==1) CALL resetnox(cbox)

  ! open files to save results
  OPEN(lro2,FILE='outdat.ro2',FORM='formatted')  
  IF (soa_fg>0) OPEN(lmaer, FILE='outdat.maer',FORM='formatted')

! open the ncdf file 
! write header and species to the ncdf file
  IF (soa_fg==1) THEN 
    CALL opennc('outdat.nc',numsp,mxlsp,nsat)
    CALL wrtnc_const(chrsp,namsat)
  ELSE
    CALL opennc('outdat.nc',numsp,mxlsp)
    CALL wrtnc_const(chrsp) 
  ENDIF

! -------------------
! INITIALIZE THE RUN
! -------------------
  DO i=1,numsp
    IF (cbox(i) < 1.0E-40) cbox(i)=1.0E-40  ! set a minimum concentration
  ENDDO
  PRINT*, 'Initialize the simuation ...'
  time = tstart
  deltat = (tstop-tstart)/REAL(ntstep)        ! time step 
  dtmax=deltat                                ! max delta time - see solver (2step)
  timemod = MODULO(time,86400.)               ! get time modulo 24h 
  iskip = 0                                   ! initialize the skip counter 

  ! set conditions for current time (temperature, rh, etc) 
  CALL getconstrain(timemod,temp,sumc,rh,water)

  ! compute the partitioning and mass concentration of SOA (in µg/m3)
  IF (soa_fg == 1) THEN       ! Equilibrium approach 
    caer(:)=0.
    CALL soapartition(lout,nsat,idgsat,temp,nvoc,rt4psat,rpsat, & 
                      eheat,cbox,caer)
    maer = SUM( wmol(idgsat(1:nsat)) * caer(1:nsat) ) * multimass
   
  ELSE IF (soa_fg == 2) THEN  ! mass transfer approach
    maer = SUM(wmol(idasat(1:nsat)) * cbox(idasat(1:nsat)) ) * multimass
    WRITE(lmaer,*) time/3600., maer 
  ENDIF  

  ! write data to the ncdf file
  IF (soa_fg==1) THEN 
    CALL wrtnc_tchge(time,cbox,caer)
  ELSE 
    CALL wrtnc_tchge(time,cbox) 
  ENDIF

! --------------------------------------
! SOLVE THE EQUATIONS - INTEGRATION LOOP
! --------------------------------------
  timeloop: &
  DO
    IF (time >= tstop) EXIT timeloop       ! end of the simulation
    tout = time + deltat                   ! time out for the current time step   
    iskip = iskip + 1                      ! record the number of skip so far
    timemod = MODULO(time,86400.)          ! get time modulo 24h

    ! set conditions for current time (temperature, rh, pressure ...)
    CALL getconstrain(timemod,temp,sumc,rh,water)
     
    ! compute photolysis rate constant for current time
    CALL getjvalue(numhv,idhv,hvfact,time,arrhcf)

    ! set CH3O2 and compute RO2 concentration for current time
    cmeo2=cbox(idch3o2)
    IF (cmeo2 < 1E-13) cmeo2=1E-13

    cro2(:)=0.
    DO i=1,ncpero ; DO k=1,ro2map(i)%nro2
        cro2(i)=cro2(i)+cbox(ro2map(i)%idro2(k))
    ENDDO ; ENDDO

    ! compute the rate constant
    CALL ratecoef(numre,numo2,num_m,numfo,numextra,nummeo2,    &
                numiso,numreacro2,ncpero, &
                ido2,id_m,idfo,idextra,idmeo2,idreacro2,idiso, &
                arrhcf,focf,extracf,isocf,         &
                cmeo2,cro2,temp,sumc,water,qfor)
    
    ! compute the rate constant for mass transfer 
    IF (soa_fg == 2) THEN
      CALL mtrat(numain,numaou,numwin,numwou,idain,idaou,idwin,idwou,  &
                 woucf,temp,wmol,pnc,nsat,rt4psat,rpsat,eheat,         &
                 mech2sat,gairefor,gwirefor,nvoc,cbox,idasat,     &
                 rh,qfor)
    ENDIF

    ! call solver
    PRINT'("time: ",f9.1)',time
    CALL twostep(numsp,time,tout,dtmin,dtmax,qfor,numit,atol,rtol, &
                 soa_fg,nsat,idasat,numaou,idaou,nvoc,cbox)

    ! reset current time, frozen stuff and processes not acccounted 
    time=tout    
    IF (nfix /= 0) CALL resetc(cbox)
    IF (noxfix == 1) CALL resetnox(cbox)
    IF (soa_fg == 1) CALL particleloss(soa_fg,kdilu,kaerloss,nsat,deltat,caer)   

    ! compute SOA mass concentration (in µg/m3)
    IF (soa_fg == 1) THEN          
      CALL soapartition(lout,nsat,idgsat,temp,nvoc,rt4psat,rpsat, & 
                        eheat,cbox,caer)
      maer = SUM( wmol(idgsat(1:nsat))*caer(1:nsat) ) * multimass
    ELSE IF (soa_fg == 2) THEN
      maer = SUM( wmol(idasat(1:nsat)) * cbox(idasat(1:nsat)) ) * multimass
    ENDIF  
    
    ! write solution to save-file 
    IF (iskip == nskip) THEN
      WRITE(lmaer,*) time/3600., maer 
      WRITE(lro2,'(10(E13.3))') time/(3600.),(cro2(i),i=1,ncpero) 
      iskip = 0 ! reset iskip
      
      ! save record to ncdf file
      IF (soa_fg==1) THEN 
        CALL wrtnc_tchge(time,cbox,caer)
      ELSE 
        CALL wrtnc_tchge(time,cbox) 
      ENDIF

    ENDIF

  ENDDO timeloop

! -------------------------------------
! END OF THE SIMULATION
! -------------------------------------
  CLOSE(lout) ; CLOSE(lro2) ; CLOSE(lmaer) 

! close the ncdf file
  CALL ncdf_close()
  
  PRINT*, 'no error - end of simulation reached' 
END PROGRAM boxmod
