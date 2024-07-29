	PROGRAM SOUND
C
C   ***   This program accepts input from the sounding contained   ***
C   ***     in the file <<sounding.txt>> and calculates various    ***
C   ***       properties of samples of air raised or lowered       ***
C   ***                   to different levels.                     ***
C
	PARAMETER (NA=2000,NK=20,NS=2000, NAT=2000)
C
	REAL T(NA), P(NA), R(NA), TLR(NK,NA), TLP(NK,NA)
	REAL PTEM(NAT), RTEM(NAT), TTEM(NAT),EVTEM(NAT),ESTEM(NAT)
	REAL EV(NA), ES(NA), TVRDIF(NK,NA), TVPDIF(NK,NA)
	REAL TLVR(NK,NA), TLVP(NK,NA), LW(NK,NA), TVD(NA)
	REAL CAPER(NS,NK), CAPEP(NS,NK), DCAPE(NS,NK), PAR(NS,NK)
	REAL CAPERAV(NK), CAPEPAV(NK), STDCAPEP(NK),STDCAPER(NK)
	REAL PAP(NS,NK), NAR(NS,NK), NAP(NS,NK), PL(NK)
	REAL ATVRD(NS,NK,NA),ATVPD(NS,NK,NA)
	REAL TRDBAR(NK,NA),TPDBAR(NK,NA),STDR(NK,NA),STDP(NK,NA)
	INTEGER JJIR(NK,NA),JJIP(NK,NA)
	CHARACTER*12 AB
	CHARACTER*16 AA
C	CHARACTER NAME*15, OBS*15, PREP*2, OTIME*3, MONTH*3, CNAM*4
	CHARACTER*5 CNAM
	CHARACTER*15 DUMMY1, DUMMY2
	CHARACTER*11 DATETIME
        character*20 soundfile,outfile
C
C   ***   ASSIGN VALUES OF THERMODYNAMIC CONSTANTS     ***
C
        CPD=1005.7
	  CPV=1870.0
        CL=2500.0
        CPVMCL=2320.0
        RV=461.5
        RD=287.04
        EPS=RD/RV
        ALV0=2.501E6
C
	  P1BAR=0.0
	  NPC=0
	  NMAX=0
	DO 8 I=1,NK
	 DO 7 J=1,NA
	  JJIR(I,J)=0
	  JJIP(I,J)=0
	  TRDBAR(I,J)=0.0
	  TPDBAR(I,J)=0.0
	  STDR(I,J)=0.0
	  STDP(I,J)=0.0
	  DO 6 K=1,NS
	   ATVRD(K,I,J)=0.0
	   ATVPD(K,I,J)=0.0
    6   CONTINUE
    7	 CONTINUE
    8 CONTINUE
	NSTOT=0
C
	DO 900 JJ=1,1
C
C
C   ***   Read in the sounding from the file <<sounding>>          ***
c
        write(*,*)'input file name'
        read(*,'(a)')soundfile
        write(*,*)'output file name'
        read(*,'(a)')outfile
        OPEN(UNIT=11,FILE=soundfile,STATUS='OLD')
c       OPEN(UNIT=11,FILE='idealsound.txt',STATUS='OLD')
c	OPEN(UNIT=11,FILE='sounding.txt',STATUS='OLD')
c	READ(11,*)
c	READ(11,*)
c	READ(11,*)
c	READ(11,*)
c	READ(11,*)
c	READ(11,*)
c	READ(11,*)
c	READ(11,*)
c	READ(11,*)
c	READ(11,*)
c	READ(11,*)
c	READ(11,*)
c	READ(11,*)
c	READ(11,*)
C
	KK=1800
	I=0
	OPEN(UNIT=10,FILE=outfile,STATUS='UNKNOWN')
c	OPEN(UNIT=10,FILE='modsound.txt',STATUS='UNKNOWN')
	DO 20 KI=1,KK
c	 READ(11,15,ERR=21,END=21)R1,J1,J2,R2,PK,TK,RK
   15	 FORMAT(1X,F6.1,9X,F5.1,2X,F5.1)
c	 READ(11,*,ERR=21,END=21)R1,J1,J2,R2,PK,TK,RK
         read(11,*,err=21,end=21)PK,TK,RK
         IF(PK.LT.78.0.OR.TK.GT.100.0.OR.RK.GT.100.0)GOTO 20
         IF(PK.lt.-99..or.TK.LT.-99..OR.RK.LT.-99.0) goto 20
          I=I+1
          PTEM(I)=PK
          TTEM(I)=TK
          RTEM(I)=RK
	    IF(TTEM(I).LT.-200.0)TTEM(I)=-200.0
          EVTEM(I)=6.112*EXP(17.67*RTEM(I)/(243.5+RTEM(I)))
          ESTEM(I)=6.112*EXP(17.67*TTEM(I)/(243.5+TTEM(I)))
	    EVTEM(I)=MIN(EVTEM(I),ESTEM(I))
          RTEM(I)=0.622*EVTEM(I)/(PTEM(I)-EVTEM(I))
          TTEM(I)=TTEM(I)+273.15
	    WRITE(10,*)PK,TK,EVTEM(I)/ESTEM(I)
	    WRITE(*,*)PK,TK,EVTEM(I)/ESTEM(I)
   20	CONTINUE
   21	CONTINUE
c	READ(11,*)DUMMY1,DUMMY2,CNAM
c	if(dummy2.eq.'identifier:    ')then
c	 READ(11,*)DUMMY1,DUMMY2,ID
c	 READ(11,*)DUMMY1,DUMMY2,idate,itime
c	else
c	 READ(11,*)DUMMY1,DUMMY2,idate,itime
c	end if
	N=I
	CLOSE(UNIT=11)
	CLOSE(10)
C
C   *** Write header  ***
C
        write(*,*)"enter year,month,day, time "
        read(*,*)IYEAR,IMONTH,IDAY,ITIME
        CNAM = "Mirai"
c	IYEAR=IDATE/10000
c	IY1=IYEAR
c	IF(IYEAR.LT.50)THEN
c	 IYEAR=IYEAR+2000
c	ELSE
c	 IYEAR=IYEAR+1900
c	END IF
c	IMONTH=IDATE/100-100*IY1
c	IDAY=IDATE-10000*IY1-100*IMONTH
	OPEN(UNIT=12,FILE='header.txt',STATUS='UNKNOWN')
      WRITE(12,*)CNAM,' ',ITIME,' ',IDAY,' ',IMONTH,' ',IYEAR
      CLOSE(12)
C
c	IF(PTEM(1).LT.995.0.OR.TTEM(1).LT.273.0.OR.TTEM(1).GT.306.0.OR.
c     1    RTEM(1).LT.0.002.OR.PTEM(N).GE.500.0.OR.PTEM(1).GT.1050.0)
c     2    GOTO 900
	NSTOT=NSTOT+1
C
C   ***  Interpolate soundings to 5 mb intervals
C
	P(1)=PTEM(1)
	T(1)=TTEM(1)
	R(1)=RTEM(1)
	EV(1)=EVTEM(1)
	ES(1)=ESTEM(1)
	NL=(1005.-PTEM(N))/5.0+0.001
	DO I=2,NL
	 P(I)=1005.0-5.0*FLOAT(I)
	 DO J=2,N
	  IF(PTEM(J).LE.P(I))THEN
	   T(I)=(TTEM(J)*(PTEM(J-1)-P(I))+TTEM(J-1)*(P(I)-PTEM(J)))/
     1     (PTEM(J-1)-PTEM(J))
	   R(I)=(RTEM(J)*(PTEM(J-1)-P(I))+RTEM(J-1)*(P(I)-PTEM(J)))/
     1     (PTEM(J-1)-PTEM(J))
	   EV(I)=(EVTEM(J)*(PTEM(J-1)-P(I))+EVTEM(J-1)*(P(I)-PTEM(J)))/
     1     (PTEM(J-1)-PTEM(J))
	   ES(I)=(ESTEM(J)*(PTEM(J-1)-P(I))+ESTEM(J-1)*(P(I)-PTEM(J)))/
     1     (PTEM(J-1)-PTEM(J))
	   GOTO 23
	  END IF
	 END DO
   23 CONTINUE
	END DO
c	print*, p(nl),t(nl),r(nl),p(1),t(1),r(1)
	N=NL
      NMAX=MAX(N,NMAX)
C
C   ***  Begin outer loop, which cycles through parcel origin levels I ***
C  
	DO 500 I=1,NK
C
C   ***  Define various conserved parcel quantities: reversible   ***
C   ***        entropy, S, pseudo-adiabatic entropy, SP,          *** 
C   ***                   and enthalpy, AH                        ***
C
	PL(I)=P(I)
	RS=EPS*ES(I)/(P(I)-ES(I))
	ALV=ALV0-CPVMCL*(T(I)-273.15)
	EM=MAX(EV(I),1.0E-6) 
	S=(CPD+R(I)*CL)*LOG(T(I))-RD*LOG(P(I)-EV(I))+
     1    ALV*R(I)/T(I)-R(I)*RV*LOG(EM/ES(I))
	SP=CPD*LOG(T(I))-RD*LOG(P(I)-EV(I))+
     1    ALV*R(I)/T(I)-R(I)*RV*LOG(EM/ES(I))
	AH=(CPD+R(I)*CL)*T(I)+ALV*R(I)
C      
C   ***  Find the temperature and mixing ratio of the parcel at   ***
C   ***    level I saturated by a wet bulb process                ***
C
	SLOPE=CPD+ALV*ALV*RS/(RV*T(I)*T(I))
	TG=T(I)
	RG=RS  
	DO 100 J=1,20 
	 ALV1=ALV0-CPVMCL*(TG-273.15)
	 AHG=(CPD+CL*RG)*TG+ALV1*RG
	 TG=TG+(AH-AHG)/SLOPE
	 TC=TG-273.15
	 ENEW=6.112*EXP(17.67*TC/(243.5+TC))
	 RG=EPS*ENEW/(P(I)-ENEW)
  100	CONTINUE
C   
C   ***  Calculate conserved variable at top of downdraft   ***
C
	EG=RG*P(I)/(EPS+RG)
	SPD=CPD*LOG(TG)-RD*LOG(P(I)-EG)+
     1    ALV1*RG/TG
	TVD(I)=TG*(1.+RG/EPS)/(1.+RG)-T(I)*(1.+R(I)/EPS)/
     1    (1.+R(I))
	IF(P(I).LT.100.0)TVD(I)=0.0
	RGD0=RG
	TGD0=TG
C
C   ***   Find lifted condensation pressure     ***
C
	RH=R(I)/RS
	RH=MIN(RH,1.0)
	CHI=T(I)/(1669.0-122.0*RH-T(I))
	PLCL=1.0
	IF(RH.GT.0.0)THEN
	 PLCL=P(I)*(RH**CHI)
	END IF
C
C   ***  Begin updraft loop   ***
C
	SUM=0.0
	RG0=R(I)
	TG0=T(I)
	DO 200 J=I,N
C
C   ***  Calculate estimates of the rates of change of the entropies  ***
C   ***           with temperature at constant pressure               ***
C  
	 RS=EPS*ES(J)/(P(J)-ES(J))
	 ALV=ALV0-CPVMCL*(T(J)-273.15)
	 SL=(CPD+R(I)*CL+ALV*ALV*RS/(RV*T(J)*T(J)))/T(J)
	 SLP=(CPD+RS*CL+ALV*ALV*RS/(RV*T(J)*T(J)))/T(J)
C   
C   ***  Calculate lifted parcel temperature below its LCL   ***
C
	 IF(P(J).GE.PLCL)THEN
	  TLR(I,J)=T(I)*(P(J)/P(I))**(RD/CPD)
	  TLP(I,J)=TLR(I,J) 
	  LW(I,J)=0.0
	  TLVR(I,J)=TLR(I,J)*(1.+R(I)/EPS)/(1.+R(I))
	  TLVP(I,J)=TLVR(I,J)
	  TVRDIF(I,J)=TLVR(I,J)-T(J)*(1.+R(J)/EPS)/(1.+R(J))
	  TVPDIF(I,J)=TVRDIF(I,J)
	 ELSE
C
C   ***  Iteratively calculate lifted parcel temperature and mixing   ***
C   ***    ratios for both reversible and pseudo-adiabatic ascent     ***
C
	 TG=T(J)
	 RG=RS
	 DO 150 K=1,20
	  EM=RG*P(J)/(EPS+RG)
	  ALV=ALV0-CPVMCL*(TG-273.15)
	  SG=(CPD+R(I)*CL)*LOG(TG)-RD*LOG(P(J)-EM)+
     1      ALV*RG/TG
	  TG=TG+(S-SG)/SL  
	  TC=TG-273.15
	  ENEW=6.112*EXP(17.67*TC/(243.5+TC))
	  RG=EPS*ENEW/(P(J)-ENEW)           
  150	 CONTINUE
	 TLR(I,J)=TG
	 TLVR(I,J)=TG*(1.+RG/EPS)/(1.+R(I))
	 LW(I,J)=R(I)-RG
	 LW(I,J)=MAX(0.0,LW(I,J))
	 TVRDIF(I,J)=TLVR(I,J)-T(J)*(1.+R(J)/EPS)/(1.+R(J))
c	 IF(TVRDIF(I,J).GE.7.0)GOTO 900
C
C   ***   Now do pseudo-adiabatic ascent   ***
C
	 TG=T(J)
	 RG=RS
	 DO 180 K=1,20 
	  CPW=0.0
	  IF(J.GT.1)THEN
	   CPW=SUM+CL*0.5*(RG0+RG)*(LOG(TG)-LOG(TG0))
	  END IF
	  EM=RG*P(J)/(EPS+RG)
	  ALV=ALV0-CPVMCL*(TG-273.15)
	  SPG=CPD*LOG(TG)-RD*LOG(P(J)-EM)+CPW+
     1      ALV*RG/TG
	  TG=TG+(SP-SPG)/SLP  
	  TC=TG-273.15
	  ENEW=6.112*EXP(17.67*TC/(243.5+TC))
	  RG=EPS*ENEW/(P(J)-ENEW)           
  180	 CONTINUE
	 TLP(I,J)=TG
	 TLVP(I,J)=TG*(1.+RG/EPS)/(1.+RG)
	 TVPDIF(I,J)=TLVP(I,J)-T(J)*(1.+R(J)/EPS)/(1.+R(J))
	 RG0=RG
	 TG0=TG
	 SUM=CPW
 	END IF
  200	CONTINUE
	IF(I.EQ.1)GOTO 500
C
C   ***  Begin downdraft loop   ***
C
	SUM2=0.0
	DO 300 J=I-1,1,-1
C
C   ***  Calculate estimate of the rate of change of entropy          ***
C   ***           with temperature at constant pressure               ***
C  
	 RS=EPS*ES(J)/(P(J)-ES(J))
	 ALV=ALV0-CPVMCL*(T(J)-273.15)
	 SLP=(CPD+RS*CL+ALV*ALV*RS/(RV*T(J)*T(J)))/T(J)
	 TG=T(J)
	 RG=RS
C
C   ***  Do iteration to find downdraft temperature   ***
C
	 DO 250 K=1,20
	  CPW=SUM2+CL*0.5*(RGD0+RG)*(LOG(TG)-LOG(TGD0))
	  EM=RG*P(J)/(EPS+RG)
	  ALV=ALV0-CPVMCL*(TG-273.15)
	  SPG=CPD*LOG(TG)-RD*LOG(P(J)-EM)+CPW+
     1      ALV*RG/TG
	  TG=TG+(SPD-SPG)/SLP  
	  TC=TG-273.15
	  ENEW=6.112*EXP(17.67*TC/(243.5+TC))
	  RG=EPS*ENEW/(P(J)-ENEW)           
  250	 CONTINUE
	 SUM2=CPW
	 TGD0=TG
	 RGD0=RG
	 TLP(I,J)=TG
	 TLVP(I,J)=TG*(1.+RG/EPS)/(1.+RG)
	 TVPDIF(I,J)=TLVP(I,J)-T(J)*(1.+R(J)/EPS)/(1.+R(J))
	 IF(P(I).LT.100.0)TVPDIF(I,J)=0.0
	 TVPDIF(I,J)=MIN(TVPDIF(I,J),0.0)
	 TLR(I,J)=T(J)
	 TLVR(I,J)=T(J)
	 TVRDIF(I,J)=0.0
	 LW(I,J)=0.0
  300	CONTINUE
  500 CONTINUE
C
C  ***  Begin loop to find CAPE, PA, and NA from reversible and ***
C  ***            pseudo-adiabatic ascent, and DCAPE            ***
C
	DO 800 I=1,NK
	 CAPER(JJ,I)=0.0
	 CAPEP(JJ,I)=0.0
	 DCAPE(JJ,I)=0.0
	 PAP(JJ,I)=0.0
	 PAR(JJ,I)=0.0
	 NAP(JJ,I)=0.0
	 NAR(JJ,I)=0.0
C
C   ***   Find lifted condensation pressure     ***
C
	RS=EPS*ES(I)/(P(I)-ES(I))
	RH=R(I)/RS
	RH=MIN(RH,1.0)
	CHI=T(I)/(1669.0-122.0*RH-T(I))
	PLCL=1.0
	IF(RH.GT.0.0)THEN
	 PLCL=P(I)*(RH**CHI)
	END IF
C
C   ***  Find lifted condensation level and maximum level   ***
C   ***               of positive buoyancy                  ***
C
	 ICB=N
	 INBR=1
	 INBP=1
	 DO 550 J=N,I,-1
	  IF(P(J).LT.PLCL)ICB=MIN(ICB,J)
	  IF(TVRDIF(I,J).GT.0.0)INBR=MAX(INBR,J)
	  IF(TVPDIF(I,J).GT.0.0)INBP=MAX(INBP,J)
  550	 CONTINUE
C
C   ***  Do updraft loops        ***
C
	 IF(INBR.GT.I)THEN
	  DO 600 J=I+1,INBR
	   TVM=0.5*(TVRDIF(I,J)+TVRDIF(I,J-1))
	   PM=0.5*(P(J)+P(J-1))
	   IF(TVM.LE.0.0)THEN
	    NAR(JJ,I)=NAR(JJ,I)-RD*TVM*(P(J-1)-P(J))/PM
	   ELSE
	    PAR(JJ,I)=PAR(JJ,I)+RD*TVM*(P(J-1)-P(J))/PM
	   END IF
  600	  CONTINUE
	  CAPER(JJ,I)=PAR(JJ,I)-NAR(JJ,I)
	 END IF
	 IF(INBP.GT.I)THEN
	  DO 650 J=I+1,INBP
	   TVM=0.5*(TVPDIF(I,J)+TVPDIF(I,J-1))
	   PM=0.5*(P(J)+P(J-1))
	   IF(TVM.LE.0.0)THEN
	    NAP(JJ,I)=NAP(JJ,I)-RD*TVM*(P(J-1)-P(J))/PM
	   ELSE
	    PAP(JJ,I)=PAP(JJ,I)+RD*TVM*(P(J-1)-P(J))/PM
	   END IF
  650	  CONTINUE
	  CAPEP(JJ,I)=PAP(JJ,I)-NAP(JJ,I)
	 END IF
C
C  ***       Find DCAPE     ***
C
	 IF(I.EQ.1)GOTO 800
	 DO 700 J=I-1,1,-1
	  TVDIFM=TVPDIF(I,J+1)
	  IF(I.EQ.(J+1))TVDIFM=TVD(I)
	  TVM=0.5*(TVPDIF(I,J)+TVDIFM)
	  PM=0.5*(P(J)+P(J+1))
	  IF(TVM.LT.0.0)THEN
	   DCAPE(JJ,I)=DCAPE(JJ,I)-RD*TVM*(P(J)-P(J+1))/PM
	  END IF
  700	 CONTINUE	  	
  800	CONTINUE
C
	DO 520 I=1,NK
	 INBR=1
	 INBP=1
	 DO 505 J=N,I,-1
	  IF(TVRDIF(I,J).GT.0.0)INBR=MAX(INBR,J)
	  IF(TVPDIF(I,J).GT.0.0)INBP=MAX(INBP,J)
  505	 CONTINUE
	  IMAX=MAX(INBR,I)
	  IMAX=N
	IF(IMAX.GT.I)THEN
	DO 508 J=I,IMAX
	 IF(TVRDIF(I,J).NE.0.0)THEN
	  JJIR(I,J)=JJIR(I,J)+1
	  ATVRD(JJIR(I,J),I,J)=TVRDIF(I,J)
	 END IF
  508	CONTINUE
	END IF
	IMAX=MAX(I,INBP)
	IMAX=N
	DO 510 J=I,IMAX
	 IF(TVPDIF(I,J).NE.0.0)THEN
	  JJIP(I,J)=JJIP(I,J)+1
	  ATVPD(JJIP(I,J),I,J)=TVPDIF(I,J)
	 END IF
  510	CONTINUE
  520	CONTINUE	
	NPC=NPC+1
	P1BAR=P1BAR+P(1)
  900	CONTINUE
  901 CONTINUE
C
	CLOSE(UNIT=15)
C
	P1BAR=P1BAR/FLOAT(NPC)
	P(1)=P1BAR
	PL(1)=P1BAR
	DO 920 I=1,NK     
	PARS=0.0
	PAPS=0.0
	ANAPS=0.0
	ANARS=0.0
	CAPERS=0.0
	CAPEPS=0.0
	DCAPES=0.0
	STDCAPEP(I)=0.0
	STDCAPER(I)=0.0
	JPAR=0
	JPAP=0
	JNAR=0
	JNAP=0
	JCAPER=0
	JCAPEP=0
	JDCAPE=0
	DO 903 JJ=1,NSTOT
	 IF(PAP(JJ,I).NE.0.0)THEN
	  PAPS=PAPS+PAP(JJ,I)
	  JPAP=JPAP+1
	 END IF
	 IF(PAR(JJ,I).NE.0.0)THEN
	  PARS=PARS+PAR(JJ,I)
	  JPAR=JPAR+1
	 END IF
	 IF(NAP(JJ,I).NE.0.0)THEN
	  ANAPS=ANAPS+NAP(JJ,I)
	  JNAP=JNAP+1
	 END IF
	 IF(NAR(JJ,I).NE.0.0)THEN
	  ANARS=ANARS+NAR(JJ,I)
	  JNAR=JNAR+1
	 END IF
	 IF(CAPER(JJ,I).NE.0.0)THEN
	  CAPERS=CAPERS+CAPER(JJ,I)
	  JCAPER=JCAPER+1
	 END IF
	 IF(CAPEP(JJ,I).NE.0.0)THEN
	  CAPEPS=CAPEPS+CAPEP(JJ,I)
	  JCAPEP=JCAPEP+1
	 END IF
	 IF(DCAPE(JJ,I).NE.0.0)THEN
	  DCAPES=DCAPES+DCAPE(JJ,I)
	  JDCAPE=JDCAPE+1
	 END IF
  903	CONTINUE 
	JNAP=MAX(JNAP,1)
	JNAR=MAX(JNAR,1)
	JPAP=MAX(JPAP,1)
	JPAR=MAX(JPAR,1)
	JCAPER=MAX(JCAPER,1)
	JCAPEP=MAX(JCAPEP,1)
	JDCAPE=MAX(JDCAPE,1)
	PAP(1,I)=PAPS/FLOAT(JPAP)
	PAR(1,I)=PARS/FLOAT(JPAR)
	NAP(1,I)=ANAPS/FLOAT(JNAP)
	NAR(1,I)=ANARS/FLOAT(JNAR)
	CAPEPAV(I)=CAPEPS/FLOAT(JCAPEP)
	CAPERAV(I)=CAPERS/FLOAT(JCAPER)
	DCAPE(1,I)=DCAPES/FLOAT(JDCAPE)
	DO 1010 JJ=1,JCAPEP
	 IF(CAPEP(JJ,I).NE.0.0)THEN
	  STDCAPEP(I)=STDCAPEP(I)+(CAPEP(JJ,I)-CAPEPAV(I))**2
	 END IF
 1010	CONTINUE
	STDCAPEP(I)=SQRT(STDCAPEP(I)/FLOAT(JCAPEP))
	DO 1020 JJ=1,JCAPER
	 IF(CAPER(JJ,I).NE.0.0)THEN
	  STDCAPER(I)=STDCAPER(I)+(CAPER(JJ,I)-CAPERAV(I))**2
	 END IF
 1020	CONTINUE
	STDCAPER(I)=SQRT(STDCAPER(I)/FLOAT(JCAPER))
	DO 915 J=1,NMAX
	 IF(JJIR(I,J).GT.0)THEN
	 DO 905 JJ=1,JJIR(I,J)
	  TRDBAR(I,J)=TRDBAR(I,J)+ATVRD(JJ,I,J)
  905	 CONTINUE
	 TRDBAR(I,J)=TRDBAR(I,J)/FLOAT(JJIR(I,J))
	 TRDBAR(I,J)=MAX(TRDBAR(I,J),-4.0)
	 DO 907 JJ=1,JJIR(I,J)
	  STDR(I,J)=STDR(I,J)+(ATVRD(JJ,I,J)-TRDBAR(I,J))**2
  907	 CONTINUE
	 STDR(I,J)=SQRT(STDR(I,J)/FLOAT(JJIR(I,J)))
	 END IF
	 IF(JJIP(I,J).GT.0)THEN
	 DO 908 JJ=1,JJIP(I,J)
	  TPDBAR(I,J)=TPDBAR(I,J)+ATVPD(JJ,I,J)
  908	 CONTINUE
	 TPDBAR(I,J)=TPDBAR(I,J)/FLOAT(JJIP(I,J))
	 TPDBAR(I,J)=MAX(TPDBAR(I,J),-4.0)
	 DO 909 JJ=1,JJIP(I,J)
	  STDP(I,J)=STDP(I,J)+(ATVPD(JJ,I,J)-TPDBAR(I,J))**2
  909	 CONTINUE
	 STDP(I,J)=SQRT(STDP(I,J)/FLOAT(JJIP(I,J)))
	 END IF
  915	CONTINUE
  920	CONTINUE
C
C  ***  Write values of PA, NA, CAPE and DCAPE   ***
C
	OPEN(UNIT=12, FILE='cape.out', STATUS='unknown')
C
	WRITE(12,921)NPC
  921	FORMAT(5X,'Number of soundings = ',I4)
	WRITE(12,801)
  801   FORMAT(20X,'ALL AREAS IN UNITS OF J/kg',//)
	WRITE(12,802)
  802	FORMAT(1X,'Origin',4X,'Rev.',4X,'P.A.',4X,'Rev.',4X,'P.A.',
     1   4X,'Rev.',4X,'P.A.',4x,'Rev.',4x,'P.A.')
	WRITE(12,804)
  804	FORMAT(1X,'p (mb)',4X,' PA ',4X,' PA ',4X,' NA ',4X,' NA ',
     1   4X,'CAPE',4X,'CAPE',4X,'DCAPE',1X,'STDCAPE',4X,'STDCAPE')
	WRITE(12,806)
  806	FORMAT(1X,'------',4X,'----',4X,'----',4X,'----',4X,'----',
     1   4X,'----',4X,'----',4X,'-----',4X,'-------',
     2    4X,'-------')
	DO 820 I=1,NK
	 WRITE(12,810)PL(I),PAR(1,I),PAP(1,I),NAR(1,I),
     1     NAP(1,I),CAPERAV(I),CAPEPAV(I),DCAPE(1,I),STDCAPEP(I),
     2	  STDCAPER(I)
  810	 FORMAT(1X,F6.1,9(F8.1))
  820	CONTINUE
	CLOSE(12)
C
C  ***  Output files for contour plots (e.g. for MATLAB)  ***
C
	OPEN(UNIT=12,FILE='p.out',STATUS='UNKNOWN')
	DO I=1,NMAX
	 WRITE(12,*)P(I)
	END DO	
	CLOSE(12)
C
	OPEN(UNIT=12,FILE='porig.out',STATUS='UNKNOWN')
	DO I=1,NK
	 WRITE(12,*)PL(I)
	END DO
	CLOSE(12)
C
	OPEN(UNIT=12,FILE='tdifrev.out',STATUS='UNKNOWN')
	DO I=1,NMAX
	 WRITE(12,1000)(TRDBAR(J,I),J=1,NK)
	END DO
	CLOSE(12)
C
	OPEN(UNIT=12,FILE='tstdrev.out',STATUS='UNKNOWN')
	DO I=1,NMAX
	 WRITE(12,1000)(STDR(J,I),J=1,NK)
	END DO
	CLOSE(12)
C
	OPEN(UNIT=12,FILE='tdifpseudo.out',STATUS='UNKNOWN')
	DO I=1,NMAX
	 WRITE(12,1000)(TPDBAR(J,I),J=1,NK)
	END DO
	CLOSE(12)
C
	OPEN(UNIT=12,FILE='tstdpseudo.out',STATUS='UNKNOWN')
	DO I=1,NMAX
	 WRITE(12,1000)(STDP(J,I),J=1,NK)
	END DO
	CLOSE(12)
C
 1000 FORMAT(21(1X,F7.3))
C
	STOP
	END

