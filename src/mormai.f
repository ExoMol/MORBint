      PROGRAM MORINT
      IMPLICIT REAL*8 (A-H,O-Z)
C*********************************************************************
C
C       PROGRAM FOR CALCULATING THE ROTATION-VIBRATION INTENSITIES
C       AND LINE STRENGTHS FOR A TRIATOMIC MOLECULE ACCORDING TO
C       THE MORBID SCHEME OF PER JENSEN
C
C*********************************************************************
C
C   THE FOLLOWING PARAMETERS ARE INITIALLY SET:
C
C       LWORK      :  LENGTH OF WORK ARRAY IN UNITS OF REALS
C       NRECUN     :  NUMBER OF UNITS FOR RECORD LENGTH (FOR DIRECT
C                  :  ACCESS FILES) NEEDED TO OBTAIN A RECORD LENGTH
C                  :  OF ONE REAL*8
C       NINPRE     :  NUMBER OF INTEGERS OCCUPYING THE SAME STORAGE
C                  :  AS ONE REAL*8
C       MEPPOW     :  THE SMALLEST NUMBER THAT CAN BE ADDED TO 1.0
C                  :  SO THAT YOUR MACHINE RECOGNIZES THE RESULT
C                  :  AS BEING > 1.0 SHOULD BE 2**(-MEPPOW)
C       LSTPOW     :  THE SMALLEST NUMBER REPRESENTABLE BY YOUR
C                  :  MACHINE SHOULD BE 2**(-LSTPOW)
C
C*********************************************************************
      PARAMETER ( LWORK = 40000000 ,  NRECUN = 8 , NINPRE = 2 ,
     1            KWORK = LWORK*NINPRE )
      PARAMETER ( ABIT = 1.0D-03 )
      PARAMETER ( MEPPOW = 55 , LSTPOW = 78 )
      REAL*8 M1,M2,M3,M,U1,U3,U13,V,RHO,EPS,
     1      CR,SR,CSE,SNE,
     2      CRE,SRE,CORO,EPSP,EPSPP,EPSPPP
C
      REAL*8 AMASS(3,10),GNS(2,2,10)
C
      REAL*8 RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1       AA1,AA3,
     2       F1A1,F2A1,F3A1,F4A1,F1A3,F2A3,F3A3,F4A3,
     3       F11,F1A11,F2A11,F3A11,F33,F1A33,F2A33,F3A33,
     4       F13,F1A13,F2A13,F3A13,
     5       F111,F1A111,F2A111,F333,F1A333,F2A333,
     6       F113,F1A113,F2A113,F133,F1A133,F2A133,
     7       F1111,FA1111,F3333,FA3333,F1113,FA1113,
     8       F1333,FA1333,F1133,FA1133,
     8       RE12 , RE32 , RHOREF , VMIN
C
      REAL*8 ETRIAL , RHOMAX , PNM1 , HBASE , HSTEP , EGUESS ,
     1      PREC
C
      REAL*8 THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1      THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2      VELLGT , PLANCK , AVOGNO , DEGRAD , RADDEG ,
     3      PI
C
      REAL*8   B11,B13,B111,B133,B113,
     2      B1111,B1333,B1113,B1133,
     3      B11111,B13333,B11113,B11333,B11133,
     4      B31,B33,B311,B333,B313,
     5      B3111,B3333,B3113,B3133,
     6      B31111,B33333,B31113,B31333,B31133
C
      REAL*8   CR1,CR3,CR11,CR33,CR13,
     2      CR111,CR333,CR113,CR133,
     3      CR1111,CR3333,CR1113,CR1333,CR1133
C
      INTEGER NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 , NFIL6 ,
     5       NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6       NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7       NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6       ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7       IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8       NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER,
     9       ISOMAX , NATTS  , V0TYPE , IVAR(55) , PARMAX ,
     1       NUMPAR , PRTINT , V2A , V2B
C
      INTEGER V1 ,V2, V3 , V2MXP1, V2P1,
     1       NSTINT , NSERIN , NSERP , NSERQ , KQUA , NTEST ,
     2       NSEPP2 , NSEQP1 , MBASIS ,
     3       MDIM , NFSYM0, NFASY0, NFSYMJ, NFASYJ,
     4       KSTYPA(2) , LSTYPA(2) , JMAX , V2MAX , JMAXP1
C
      INTEGER IQUANT(5,10)
C
      LOGICAL SYMM,POTSYM
C
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW,JDIMST(4),ISTART(4),
     1        LINDEX(4),V2VALU,V2VALL
      REAL*8 WORK(LWORK)
      CHARACTER*90 ELINE
      CHARACTER*4 SYMSYM(4),ASYSYM(4)
      INTEGER IWRK(KWORK),IWRKST,IWRKEN,INTWRK,
     1       NRECS,NRSYM0,NRASY0,NRSYMJ,NRASYJ,
     2       ISSYM0,ISASY0,ISSYMJ,ISASYJ,IENRGY,
     3       IF1ST,IF2ST,IV0ST,IRRST,IGFST,IPHIST,IDERST,
     4       IJACST,IRHSST,IXTXST,IXTYST,ICSSST,IASSST,IIA1ST,
     5       IIA2ST,IRTIRR,LENREC,       IVALST,IPHIL,IPHIR,
     6       IDERR,IDERL,ISTOST,NPHI,LENPHI,LENVAL,IIW1ST,IIW3ST,
     7       IOVEST,IVMAST,IHMAST,NOVA,NOVB,IAMAST,IEVMST,
     8       IUMAST,IEWKST,I,NFCT,LENIAR,LENINT,IRTIST,
     9       ICOMST,IEVIST,IDOMST,IDIM,NORECS,
     1       ITMAST,IESTST,ISMAST,IEWRST,INDBST,NAMDIM,
     1       IHTRST,IWMAST,TAUL,TAUU
      REAL*8 RECS
C
      COMMON /ELECTR/ EC,ZC1,ZC2,ZC3,ZELEC
C
      COMMON /ISOTOP/ IQUANT,AMASS
C
      COMMON /VALUES/ RHO,EPS,EPSP,EPSPP,EPSPPP,
     1               CR,SR,CSE,SNE,CRE,SRE,CORO
C
      COMMON /MOLCUL/ RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1               AA1,AA3,
     2               F1A1,F2A1,F3A1,F4A1,F1A3,F2A3,F3A3,F4A3,
     3               F11,F1A11,F2A11,F3A11,F33,F1A33,F2A33,F3A33,
     4               F13,F1A13,F2A13,F3A13,
     5               F111,F1A111,F2A111,F333,F1A333,F2A333,
     6               F113,F1A113,F2A113,F133,F1A133,F2A133,
     7               F1111,FA1111,F3333,FA3333,F1113,FA1113,
     8               F1333,FA1333,F1133,FA1133,
     8               RE12 , RE32 , M1 , M2 , M3 , M ,
     9               U1 , U3 , U13 , V ,
     1               SYMM
C
      COMMON /BMSC/  BOLTZ
C
      COMMON /INTEG/  ETRIAL , RHOMAX , PNM1 , HBASE , HSTEP , EGUESS ,
     1               RHOREF , VMIN , V0TYPE ,
     1               NSTINT , NSERIN , NSERP , NSERQ , KQUA , NTEST ,
     2               NSEPP2 , NSEQP1 , KSTYPA , LSTYPA
C
      COMMON /BOOKK/ NZERSY , NZERAS , NONESY , NONEAS , NTWOSY ,
     1               NTWOAS , KZERSY , KZERAS , KONESY , KONEAS ,
     2               KTWOSY , KTWOAS
C
      COMMON /DIMEN/  MBASIS ,  V2MAX  , V2MXP1 ,
     1               JMAX   , JMAXP1 , MDIM   , NFSYM0 , NFASY0 ,
     2               NFSYMJ , NFASYJ
C
      COMMON /LSFIT/  PARMAX , NUMPAR , ISOMAX , IVAR
C
      COMMON /SYS/    THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1                     THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2                     VELLGT , PLANCK , AVOGNO , PI , PREC ,
     4                     NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 ,
     5                     NFIL6 , NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6                     NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7                     NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6                     ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7                     IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8                     NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER ,
     9                       NATTS , PRTINT
C
      COMMON/BCOEFF/
     1      B11,B13,B111,B133,B113,
     2      B1111,B1333,B1113,B1133,
     3      B11111,B13333,B11113,B11333,B11133,
     4      B31,B33,B311,B333,B313,
     5      B3111,B3333,B3113,B3133,
     6      B31111,B33333,B31113,B31333,B31133
C
      COMMON/CRCOEF/
     1      CR1,CR3,CR11,CR33,CR13,
     2      CR111,CR333,CR113,CR133,
     3      CR1111,CR3333,CR1113,CR1333,CR1133
C
      COMMON /AMPDI/ IDPMA
C
      COMMON /DIPOLE/P0,PA1,PA2,PA3,PA4,PA5,PA6,PA7,PA8,
     2               P1,P1A1,P2A1,P3A1,P4A1,
     2               P3,P1A3,P2A3,P3A3,P4A3,
     3               P11,P1A11,P2A11,P3A11,P33,P1A33,P2A33,P3A33,
     4               P13,P1A13,P2A13,P3A13,
     5               P111,P1A111,P2A111,P333,P1A333,P2A333,
     6               P113,P1A113,P2A113,P133,P1A133,P2A133,
     7               P1111,PA1111,P3333,PA3333,P1113,PA1113,
     8               P1333,PA1333,P1133,PA1133,
     1               O0,OA1,OA2,OA3,OA4,OA5,OA6,OA7,OA8,
     2               O1,O1A1,O2A1,O3A1,O4A1,
     2               O3,O1A3,O2A3,O3A3,O4A3,
     3               O11,O1A11,O2A11,O3A11,O33,O1A33,O2A33,O3A33,
     4               O13,O1A13,O2A13,O3A13,
     5               O111,O1A111,O2A111,O333,O1A333,O2A333,
     6               O113,O1A113,O2A113,O133,O1A133,O2A133,
     7               O1111,OA1111,O3333,OA3333,O1113,OA1113,
     8               O1333,OA1333,O1133,OA1133
C
C
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
      EQUIVALENCE (WORK(1),IWRK(1))
      DATA SYMSYM/'A1  ','B2  ','B1  ','A2  '/
      DATA ASYSYM/'A''  ','    ','A"  ','    '/
      CALL CLOCKV ( OLDVEC , OLDTIM , 1 , 2 )
C
C
C
C ***** ELEMENTARY CHARGE IN COULOMBS
C
      EC=1.60217733D-19
C
C ***** VACUUM PERMITTIVITY IN FARAD/METER
C
      VCPRM=8.854187817D-12
C
C ***** ELEMENTARY CHARGE IN ELECTROSTATIC UNITS (E.S.U.)
C
      EC=EC/SQRT(4.0D-9*PI*VCPRM)
C
C ***** ELEMENTARY CHARGE IN DEBYE/ANGSTROM =
C
C       1.0E8 * DEBYE/CM = 1.0E-10 * E.S.U. * CM/CM
C
      EC=1.0D10*EC
C*********************************************************************
C
C       PRINT TITLE BANNER
C
C*********************************************************************
C
      CALL HEADER
C
C*********************************************************************
C
C       READ INPUT
C
C*********************************************************************
C
      CALL READIN(POTSYM , FRQLO , FRQHI , XINLIM , TRTLIM ,
     1            GNS    , ABSTMP ,
     1  V2A , NVA , IABA , V2B , NVB , IABB , IOPT ,
     1  ENGLIM )
C
C*********************************************************************
C
C       OPEN DIRECT ACCESS FILE FOR STORING BENDING WAVEFUNCTIONS
C
C*********************************************************************
C
      OPEN(UNIT=NFIL7,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=2*NSTINT*NRECUN,IOSTAT=IERR1)
      IF (IERR1 .NE. 0) THEN
              IERR2=NFIL7
              GOTO 1009
      ENDIF
C
C
C*********************************************************************
C
C       SET UP PARAMETERS FOR DIFFERENT POTENTIAL EXPANSION
C
C*********************************************************************
C
      CALL POTTYP
C
C*********************************************************************
C
C       START LOOP ON ISOTOPES
C
C*********************************************************************
C
      DO 10 I=1,NISOT
      REWIND NFIL1
      REWIND NFIL2
      REWIND NFIL3
      REWIND NFIL4
      REWIND NFIL8
      REWIND NFIL9
      REWIND NFIL10
      REWIND NFIL11
      REWIND NFIL12
      REWIND NFIL13
      REWIND NFIL14
      REWIND NFIL15
      REWIND NFIL20
      M1=AMASS(1,I)
      M2=AMASS(2,I)
      M3=AMASS(3,I)
      M=M1+M2+M3
      SYMM=(POTSYM .AND. ABS(M1-M3) .LT. 1.0E-30)
      JMAX=IQUANT(1,I)
      V2MAX=IQUANT(2,I)
      MBASIS=IQUANT(3,I)
      LSTYPA(1)=IQUANT(4,I)
      LSTYPA(2)=IQUANT(5,I)
      V2MXP1=V2MAX+1
      JMAXP1=JMAX+1
C
C*********************************************************************
C
C       OPEN DIRECT ACCESS FILE FOR STORING CALCULATED ENERGIES
C
C*********************************************************************
C
      IF (MOD(JMAX,2) .EQ. 0) THEN
            IROTE=JMAX/2+1
      ELSE
            IROTE=(JMAX+1)/2
      ENDIF
C
      IF (.NOT. SYMM) THEN
            LRODIM=2*IROTE*V2MXP1*LSTYPA(1)+2
      ELSE
            LRODIM=IROTE*(LSTYPA(1)+LSTYPA(2))*V2MXP1+2
      ENDIF
C
      OPEN(UNIT=NFIL17,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=LRODIM*NRECUN,IOSTAT=IERR1)
      IF (IERR1 .NE. 0) THEN
              IERR2=NFIL17
              GOTO 1009
      ENDIF
C
      IF (SYMM) THEN
            IDIMRO=V2MXP1*LSTYPA(1)
            IDIMCO=V2MXP1*LSTYPA(2)
      ELSE
            IDIMRO=V2MXP1*LSTYPA(1)
            IDIMCO=IDIMRO
      ENDIF
C
C*********************************************************************
C
C       OPEN DIRECT ACCESS FILE FOR STORING THE MATRIX ELEMENTS OF
C       MUZ
C
C*********************************************************************
C
      OPEN(UNIT=NFIL18,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=2*IDIMRO*IDIMCO*NRECUN,IOSTAT=IERR1)
      IF (IERR1 .NE. 0) THEN
              IERR2=NFIL18
              GOTO 1009
      ENDIF
C
C
C*********************************************************************
C
C       OPEN DIRECT ACCESS FILE FOR STORING THE MATRIX ELEMENTS OF
C       MUY
C
C*********************************************************************
C
      IF (LSTYPA(1) .GT. LSTYPA(2)) THEN
           LRECL=(V2MXP1*LSTYPA(1))**2
      ELSE
           LRECL=(V2MXP1*LSTYPA(2))**2
      ENDIF
      OPEN(UNIT=NFIL19,STATUS='SCRATCH',ACCESS='DIRECT',
     1     RECL=LRECL*NRECUN,IOSTAT=IERR1)
      IF (IERR1 .NE. 0) THEN
              IERR2=NFIL19
              GOTO 1009
      ENDIF
C
C
C*********************************************************************
C
C       SEQUENCE FOR CALCULATING THE RHO DEPENDENT FUNCTIONS
C       NECESSARY FOR THE MORBID CALCULATION.
C
C*********************************************************************
C
      IWRKST=1
      INTWRK=LWORK-IWRKST+1
C
C*********************************************************************
C
C THE RHO DEPENDENT FUNCTIONS ARE STORED ON DISK AS FOLLOWS
C
C NFIL1: FUNCTIONS NECESSARY TO DO A J=0 CALCULATION FOR AN ABA
C        MOLECULE.
C NFIL2: ADDITIONAL FUNCTIONS NECESSARY TO DO A J=0 CALCULATION
C        FOR AN ABC MOLECULE.
C NFIL3: FUNCTIONS NECESSARY TO DO A J>0 CALCULATION FOR AN ABA
C        MOLECULE.
C NFIL4: ADDITIONAL FUNCTIONS NECESSARY TO DO A J>0 CALCULATION
C        FOR AN ABC MOLECULE.
C
C*********************************************************************
C
      NFSYM0=40
      NFASY0=29
      NFSYMJ=61
      NFASYJ=46
      NZERSY=67
      NZERAS=47
      NONESY=34
      NONEAS=28
      NTWOSY=18
      NTWOAS=12
      KFSYM0=NFSYM0*NSTINT
      KFASY0=NFASY0*NSTINT
      KFSYMJ=NFSYMJ*NSTINT
      KFASYJ=NFASYJ*NSTINT
      INTSYM=70
      INTASY=50
      LNTSYM=2*INTSYM*V2MXP1*V2MXP1*JMAXP1
      LNTASY=2*INTASY*V2MXP1*V2MXP1*JMAXP1
C
C*********************************************************************
C
C THE NF VARIABLES CONTAIN THE NUMBER OF FUNCTIONS WHOSE VALUES
C ARE TO BE STORED IN NFIL1, NFIL2, NFIL3, AND NFIL4, RESPECTIVELY.
C
C WE NOW DIVIDE UP THE WORK ARRAY TO HOLD THE VARIOUS VARIABLES
C
C*********************************************************************
C
      ISSYM0=IWRKST
      ISASY0=ISSYM0+KFSYM0
      ISSYMJ=ISASY0+KFASY0
      ISASYJ=ISSYMJ+KFSYMJ
      IWRKEN=ISASYJ+KFASYJ
      IF (IWRKEN .GT. LWORK) GOTO 1001
C
C*********************************************************************
C
C RHOFCT CALCULATES THE APPROPRIATE RHO DEPENDENT FUNCTIONS AND
C STORES THEM ON UNITS NFIL1, NFIL2, NFIL3, AND NFIL4.
C
C*********************************************************************
C
      CALL RHOFCT ( WORK(ISSYM0) , WORK(ISASY0) , WORK(ISSYMJ) ,
     1             WORK(ISASYJ) ,
     1             KFSYM0 , KFASY0 , KFSYMJ , KFASYJ , POTSYM )
      CALL PRTIME( 'RHOFCT' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C WE NOW RESERVE SPACE IN THE WORK ARRAY FOR THE PURE BENDING
C ENERGIES AND FOR THE FUNCTIONS NECESSARY FOR THE NUMEROV-COOLEY
C NUMERICAL INTEGRATION. THESE FUNCTIONS ARE: F1, F2, V0+UBEND,
C IRR0, AND G. FOR DEFINITIONS SEE P. JENSEN, COMP. PHYS. REP. 1,
C 1-55 (1983).
C
C*********************************************************************
C
      IENRGY=IWRKST
      IF1ST=IENRGY+V2MXP1*JMAXP1
      IF2ST=IF1ST+NSTINT
      IV0ST=IF2ST+NSTINT
      IRRST=IV0ST+NSTINT
      IGFST=IRRST+NSTINT
      IWRKEN=IGFST+NSTINT
      IF (IWRKEN .GE. LWORK) GOTO 1002
C
C*********************************************************************
C
C NUMFCT CALCULATES THE VALUES OF THESE FUNCTIONS
C
C*********************************************************************
C
      CALL NUMFCT ( WORK(IF1ST) , WORK(IF2ST) , WORK(IV0ST) ,
     1             WORK(IRRST) , WORK(IGFST) )
      CALL PRTIME( 'NUMFCT' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C APART FROM THE SPACE ALREADY RESERVED IN THE WORK ARRAY
C WE NOW RESERVE SPACE FOR ONE BENDING WAVEFUNCTION (THE ONE
C THAT IS CURRENTLY BEING CALCULATED), FOR ITS DERIVATIVE
C AND FOR VARIOUS ARRAYS USED IN CARRYING OUT THE SERIES
C SOLUTION AROUND RHO=0.
C
C*********************************************************************
C
      IRTIST=IWRKEN
      IPHIST=IRTIST+NSTINT
      IDERST=IPHIST+NSTINT
      NSEPP2=NSERP+2
      NSEQP1=NSERQ+1
      IJACST=IDERST+NSTINT
      IRHSST=IJACST+NSERIN*NSEPP2
      IXTXST=IRHSST+NSERIN
      IXTYST=IXTXST+NSEPP2*NSEPP2
      ICSSST=IXTYST+NSEPP2
      IASSST=ICSSST+NSEPP2
      IIA1ST=IASSST+NSEQP1
      LIA1ST=ILCONV(IIA1ST,NINPRE)
      LENIAR=NSEPP2
      IIA2ST=IIA1ST+LENIAR
      LIA2ST=ILCONV(IIA2ST,NINPRE)
      IWRKEN=IIA2ST+LENIAR
      IF (IWRKEN .GT. LWORK) GOTO 1003
C
C*********************************************************************
C
C WAVFUN CALCULATES THE K=0 BENDING FUNCTIONS AND STORES THEM ON
C UNIT NFIL7.
C
C*********************************************************************
C
      CALL WAVFUN ( WORK(IENRGY) , WORK(IF1ST)  , WORK(IF2ST)  ,
     1             WORK(IV0ST)   , WORK(IRRST)  , WORK(IGFST)  ,
     2             WORK(IRTIST)  , WORK(IPHIST) , WORK(IDERST) ,
     3             WORK(IJACST)  , WORK(IRHSST) , WORK(IXTXST) ,
     4             WORK(IXTYST)  , WORK(ICSSST) , WORK(IASSST) ,
     5             IWRK(LIA1ST)  , IWRK(LIA2ST) ,
     6             LENIAR )
      CALL PRTIME( 'WAVFUN' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C WE NOW SET UP FOR CALCULATING THE INTEGRALS OVER THE BENDING
C COORDINATE NECESSARY FOR DOING A J=0 CALCULATION.
C WE SET UP SPACE FOR
C                       - THE INTEGRALS.
C                       - AS MANY VALUES OF THE RHO DEPENDENT
C                         THAT ONE NEEDS AT ANY GIVEN TIME.
C                       - A LEFT AND A RIGHT WAVEFUNCTION AND
C                         THEIR DERIVATIVES.
C
C*********************************************************************
C
      MBASP1=MBASIS+1
      NOBAS=MBASP1*(MBASP1+1)/2
      IESTST=IPHIST
      IVALST=IESTST+NFSYM0*V2MXP1*V2MXP1
      IF (.NOT.SYMM) IVALST=IVALST+NFASY0*V2MXP1*V2MXP1
      IVSYST=IVALST+NFSYM0*NSTINT
      IF (NFSYM0.LT.NFASY0 .AND. .NOT.SYMM)
     1                     IVSYST=IVALST+NFASY0*NSTINT
      IPHIL=IVSYST+NFSYMJ*NSTINT
      IF (NFSYMJ .LT. NFASYJ .AND. .NOT. SYMM)
     1                       IPHIL=IVSYST+NFASYJ*NSTINT
      IDERL=IPHIL+NSTINT
      IPHIR=IDERL+NSTINT
      IDERR=IPHIR+NSTINT
      IWRKEN=IDERR+NSTINT
      LENVAL=IVSYST-IVALST
      LENVSY=IPHIL-IVSYST
      IF (IWRKEN .GT. LWORK) GOTO 1004
C
C*********************************************************************
C
C FIRST WE DO THE DELTA K=0 BENDING MATRIX ELEMENTS
C
C*********************************************************************
C
      CALL INTZER ( WORK(IGFST) , WORK(IVALST) , WORK(IVSYST) ,
     1             WORK(IPHIL)  , WORK(IDERL)  , WORK(IPHIR)  ,
     2             WORK(IDERR)  ,
     4                      LENVAL , LENVSY ,
     1             KFSYM0 , KFASY0 , KFSYMJ , KFASYJ )
      CALL PRTIME( 'INTZER' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C THEN WE DO THE DELTA K=1 BENDING MATRIX ELEMENTS
C
C*********************************************************************
C
      CALL INTONE ( WORK(IVSYST) ,
     1             WORK(IPHIL)  , WORK(IDERL)  , WORK(IPHIR)  ,
     2             WORK(IDERR)  ,
     4             LENVSY ,
     1             KFSYM0 , KFASY0 , KFSYMJ , KFASYJ )
      CALL PRTIME( 'INTONE' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C FINALLY WE DO THE DELTA K=2 BENDING MATRIX ELEMENTS
C
C*********************************************************************
C
      CALL INTTWO ( WORK(IVSYST) ,
     1             WORK(IPHIL)  , WORK(IDERL)  , WORK(IPHIR)  ,
     2             WORK(IDERR)  ,
     4             LENVSY ,
     1             KFSYM0 , KFASY0 , KFSYMJ , KFASYJ )
      CALL PRTIME( 'INTTWO' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
C SET UP PARAMETERS FOR THE CALCULATION OF MORSE OSCILLATOR MATRIX
C ELEMENTS
C
C*********************************************************************
C
      LENIAR=NOBAS
      MDIM=V2MXP1*(LSTYPA(1)+LSTYPA(2))
      IIW1ST=IVALST
      LIW1ST=ILCONV(IIW1ST,NINPRE)
      IIW3ST=IIW1ST+LENIAR
      LIW3ST=ILCONV(IIW3ST,NINPRE)
      IOVEST=IIW3ST+LENIAR
      LENIW=LENIAR
      IWRKEN=IOVEST+2*MBASP1*MBASP1
      IF (IWRKEN .GE. LWORK) GOTO 1005
C
C*********************************************************************
C
      CALL MOPARM ( IWRK(LIW1ST) , IWRK(LIW3ST) , WORK(IOVEST) )
      CALL PRTIME( 'MOPARM' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
      NOPERS=38
      ITMAST=IWRKEN
      ISMAST=ITMAST+NOBAS*NOBAS
      IVMAST=ISMAST+NOBAS*NOBAS*NOPERS
      IUMAST=IVMAST+NOBAS*NOBAS
      IEWRST=IUMAST+NOBAS*2
      IHTRST=IEWRST+NOBAS
      ICOMST=IHTRST+NOBAS
      LCOMST=ILCONV(ICOMST,NINPRE)
      INDBST=ICOMST+LENIAR
      LNDBST=ILCONV(INDBST,NINPRE)
      IDOMST=INDBST+LENIAR
      LDOMST=ILCONV(IDOMST,NINPRE)
      IWRKEN=IDOMST+LENIAR
      IF (IWRKEN .GE. LWORK) GOTO 1006
C
C*********************************************************************
C
      CALL STRTCH ( IWRK(LIW1ST) , IWRK(LIW3ST) , WORK(IOVEST) ,
     1             WORK(ITMAST) , WORK(IESTST) , WORK(ISMAST) ,
     1             WORK(IVMAST) , WORK(IUMAST) , WORK(IEWRST) ,
     1             WORK(IHTRST) , IWRK(LCOMST) , IWRK(LNDBST) ,
     1             IWRK(LDOMST) )
      CALL PRTIME( 'STRTCH' , OLDTIM , OLDVEC )
C
C*********************************************************************
C
      IWMAST=IUMAST
      IWRKEN=IWMAST+NOBAS*NOBAS
      IF (IWRKEN .GE. LWORK) GOTO 1007
C
C*********************************************************************
C
      CALL GENSTR( IWRK(LIW1ST) , IWRK(LIW3ST) , WORK(IOVEST) ,
     1             WORK(ITMAST) , WORK(IESTST) , WORK(ISMAST) ,
     1             WORK(IVMAST) , WORK(IWMAST) , NOPERS )
      CALL PRTIME( 'GENSTR' , OLDTIM , OLDVEC )
C
      REWIND NFIL15
      WRITE (NFIL15) (WORK(III), III=ISMAST,IVMAST-1)
C
C*********************************************************************
C
C  GENERATE MATRIX ELEMENTS OF THE DIPOLE MOMENT OPERATORS
C
C*********************************************************************
C
      ISYMST=IVALST
      IASYST=ISYMST+LNTSYM
      ISMAST=IASYST+LNTASY
      IHMAST=ISMAST+NOBAS*NOBAS*NOPERS
C
      REWIND NFIL15
      READ (NFIL15) (WORK(III), III=ISMAST,IHMAST-1)
C
      IF (SYMM) THEN
            ISTOFF=0
            JSTOFF=KSTYPA(1)
            NDIMRO=LSTYPA(1)
            NDIMCO=LSTYPA(2)
            IDIMRO=V2MXP1*LSTYPA(1)
            IDIMCO=V2MXP1*LSTYPA(2)
      ELSE
            ISTOFF=0
            JSTOFF=0
            NDIMRO=LSTYPA(1)
            NDIMCO=LSTYPA(1)
            IDIMRO=V2MXP1*LSTYPA(1)
            IDIMCO=IDIMRO
      ENDIF
      IWRKEN=IHMAST+IDIMRO*IDIMCO
      IF (IWRKEN .GT. LWORK) GOTO 1012
C
      CALL MUEZ ( WORK(ISYMST) , WORK (IASYST) , WORK(ISMAST) ,
     1            WORK(IHMAST) , NDIMRO        , NDIMCO       ,
     2            IDIMRO       , IDIMCO        , LNTSYM       ,
     3            LNTASY       , ISTOFF        , JSTOFF       ,
     4            NOPERS       )
      CALL PRTIME( 'MUEZ  ' , OLDTIM , OLDVEC )
C
      ISTOFF=0
      NDIMRO=LSTYPA(1)
      IDIMRO=V2MXP1*LSTYPA(1)
      IWRKEN=IHMAST+IDIMRO*IDIMRO
      IF (IWRKEN .GT. LWORK) GOTO 1014
C
      IRECOF=0
      CALL MUEY ( WORK(ISYMST) , WORK (IASYST) , WORK(ISMAST) ,
     1            WORK(IHMAST) , NDIMRO        ,
     2            IDIMRO       ,                 LNTSYM       ,
     3            LNTASY       , ISTOFF        , IRECOF       ,
     4            NOPERS       )
      CALL PRTIME( 'MUEY  ' , OLDTIM , OLDVEC )
C
      IF (SYMM) THEN
      ISTOFF=KSTYPA(1)
      NDIMRO=LSTYPA(2)
      IDIMRO=V2MXP1*LSTYPA(2)
      IWRKEN=IHMAST+IDIMRO*IDIMRO
      IRECOF=2*(JMAX+1)
      IF (IWRKEN .GT. LWORK) GOTO 1014
C
      CALL MUEY ( WORK(ISYMST) , WORK (IASYST) , WORK(ISMAST) ,
     1            WORK(IHMAST) , NDIMRO        ,
     2            IDIMRO       ,                 LNTSYM       ,
     3            LNTASY       , ISTOFF        , IRECOF       ,
     4            NOPERS       )
      CALL PRTIME( 'MUEY  ' , OLDTIM , OLDVEC )
       ENDIF
C
C*********************************************************************
C
C  START LOOP OVER J = 0 THROUGH JMAX
C
C*********************************************************************
C
      ISYMST=IVALST
      IASYST=ISYMST+LNTSYM
      ISMAST=IASYST+LNTASY
      ISTORE=ISMAST+NOBAS*NOBAS*NOPERS
C
      REWIND NFIL15
      READ (NFIL15) (WORK(III), III=ISMAST,ISTORE-1)
C
      IROTOT=(2*JMAX+1)*V2MXP1*(LSTYPA(1)+LSTYPA(2))
      IHMAST=ISTORE+IROTOT
C
      ABSINT=8.0E-36*PI**3*AVOGNO/(3.0D+00*PLANCK*VELLGT)
      TEMPCO=-PLANCK*VELLGT/(BOLTZ*ABSTMP)
      PRTFCT=0.0D+00
      NTRANS=0
C
      DO 200 JP1=1,JMAXP1
C
      J=JP1-1
      FACJ=2.0D+00*FLOAT(J)+1.0D+00
      INDEXJ=0
      INDEXD=1
      DO 160 II=1,4
160   JDIMST(II)=0
C
      DO 180 ITAUP1=1,2
      ITAU=ITAUP1-1
C
      IF (MOD(J,2) .EQ. 0) THEN
            IROTE=J/2+1
            IROTO=J/2
      ELSE
            IROTE=(J+1)/2
            IROTO=IROTE
      ENDIF
      IF (MOD(J+ITAU,2) .NE. 0) IROTE=IROTE-1
C
C*********************************************************************
C
      IF (.NOT. SYMM) THEN
            IRODIM=(IROTE+IROTO)*V2MXP1*LSTYPA(1)
      ELSE
            IRODIM=(IROTE*LSTYPA(1)+IROTO*LSTYPA(2))*V2MXP1
      ENDIF
      IWRKEN=IHMAST+IRODIM*IRODIM
      IF (IWRKEN .GE. LWORK) GOTO 1010
C
C*********************************************************************
C
C  CLEAR ALL ELEMENTS IN THE HAMILTONIAN MATRIX-TO-BE
C
C*********************************************************************
C
      DO 100 IH=IHMAST,IWRKEN
100   WORK(IH)=0.0D+00
C
C*********************************************************************
C
C  SET UP THE EVEN K CORNER OF THE HAMILTONIAN MATRIX BLOCK
C
C*********************************************************************
C
      IF (IROTE .EQ. 0) GOTO 320
      IEO=0
      NOFSET=0
      NDIMST=LSTYPA(1)
      ISTOFF=0
      CALL DELK0 ( WORK(IENRGY) , WORK(IESTST) , WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSET , NDIMST , JP1    , IRODIM , IEO    ,
     3             ITAU   , LNTSYM , LNTASY , ISTOFF ,
     4             NOPERS       )
      NOFSRO=0
      NOFSCO=0
      NDIMRO=LSTYPA(1)
      NDIMCO=LSTYPA(1)
      ISTOFF=0
      JSTOFF=0
      CALL DELK2 (                               WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             IRODIM , IEO    , ITAU   , LNTSYM ,
     4             LNTASY , ISTOFF , JSTOFF ,
     4             NOPERS       )
C
C*********************************************************************
C
C  SET UP THE ODD K CORNER OF THE HAMILTONIAN MATRIX BLOCK
C
C*********************************************************************
C
320   CONTINUE
      IF (IROTO .EQ. 0) GOTO 340
      IEO=1
      NOFSET=IROTE*V2MXP1*NDIMST
      IF (SYMM) THEN
            NDIMST=LSTYPA(2)
            ISTOFF=KSTYPA(1)
      ENDIF
      CALL DELK0 ( WORK(IENRGY) , WORK(IESTST) , WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSET , NDIMST , JP1    , IRODIM , IEO    ,
     3             ITAU   , LNTSYM , LNTASY , ISTOFF ,
     4             NOPERS       )
      NOFSRO=NOFSET
      NOFSCO=NOFSET
      IF (SYMM) THEN
            NDIMRO=LSTYPA(2)
            NDIMCO=LSTYPA(2)
            ISTOFF=KSTYPA(1)
            JSTOFF=KSTYPA(1)
      ENDIF
      CALL DELK2 (                               WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             IRODIM , IEO    , ITAU   , LNTSYM ,
     4             LNTASY , ISTOFF , JSTOFF ,
     4             NOPERS       )
C
C*********************************************************************
C
C  SET UP DELTA K = 1 MATRIX ELEMENTS CONNECTING THE TWO CORNERS
C
C*********************************************************************
C
340   CONTINUE
      IF (IROTO .EQ. 0 .OR. IROTE .EQ. 0) GOTO 360
      IEO=0
      NOFSRO=0
      NOFSCO=NOFSET
      IF (SYMM) THEN
            NDIMRO=LSTYPA(1)
            NDIMCO=LSTYPA(2)
            ISTOFF=0
            JSTOFF=KSTYPA(1)
      ENDIF
      CALL DELK1 (                               WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             IRODIM , IEO    , ITAU   , LNTSYM ,
     4             LNTASY , ISTOFF , JSTOFF ,
     4             NOPERS       )
360   CONTINUE
C
C*********************************************************************
C
C  DIAGONALIZE THE HAMILTONIAN MATRIX BLOCK AND PRINT OUT THE
C  EIGENVALUES
C
C*********************************************************************
C
      IF (IRODIM .EQ. 0) GOTO 380
      IEIGST=IWRKEN
      IEIGVC=IEIGST+IRODIM
      IEWKST=IEIGVC+IRODIM*IRODIM
      IQNUMS=IEWKST+8*IRODIM
      IIWORK=IQNUMS+IRODIM
      IIFAIL=IIWORK+5*IRODIM
      IWRKEN=IIFAIL+IRODIM
      IF (IWRKEN .GE. LWORK) GOTO 1011
C
      CALL DIAROT ( IRODIM , WORK(IHMAST) , WORK(IEIGST) , 
     1 WORK(IEIGVC) ,
     1 WORK(IEWKST) , WORK(IQNUMS) ,
     1 J  ,         ITAU   , 1   ,
     1 WORK(IIWORK) , WORK(IIFAIL) )
C
      IF (JP1 .EQ. 1 .AND. ITAU .EQ. 0) THEN
            E00=WORK(IEIGST)
            IF (IOPT .EQ. 1) WRITE (NFIL6,6600) E00
      ENDIF
      IF (SYMM) THEN
      FACT=FACJ*GNS(ITAU+1,1,I)
      ELSE
      FACT=FACJ*GNS(1,1,I)
      ENDIF
      DO 163 II=1,IRODIM
      IRECN=(4*MOD(J,2)+2*ITAU)*LRODIM+II
      EOUT=WORK(IEIGST+II-1)-E00
      PRTFCT=PRTFCT+FACT*EXP(EOUT*TEMPCO)
163   WRITE (NFIL17,REC=IRECN) EOUT,
     .    (WORK(IHMAST+(II-1)*IRODIM+JJ-1), JJ=1,IRODIM)
      DO 165 II=1,IRODIM
      WORK(ISTORE+INDEXJ)=WORK(IEIGST+II-1)
165   INDEXJ=INDEXJ+1
      JDIMST(INDEXD)=IRODIM
C
380   INDEXD=INDEXD+1
C
C*********************************************************************
C
C
      IF (.NOT. SYMM) GOTO 460
C
C*********************************************************************
C
      IRODIM=(IROTE*LSTYPA(2)+IROTO*LSTYPA(1))*V2MXP1
      IWRKEN=IHMAST+IRODIM*IRODIM
      IF (IWRKEN .GE. LWORK) GOTO 1010
C
C*********************************************************************
C
      DO 120 IH=IHMAST,IWRKEN
120   WORK(IH)=0.0D+00
C
C*********************************************************************
C
C  SET UP THE EVEN K CORNER OF THE HAMILTONIAN MATRIX BLOCK
C
C*********************************************************************
C
      IF (IROTE .EQ. 0) GOTO 400
      IEO=0
      NOFSET=0
      NDIMST=LSTYPA(2)
      ISTOFF=KSTYPA(1)
      CALL DELK0 ( WORK(IENRGY) , WORK(IESTST) , WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSET , NDIMST , JP1    , IRODIM , IEO    ,
     3             ITAU   , LNTSYM , LNTASY , ISTOFF ,
     4             NOPERS       )
      NOFSRO=0
      NOFSCO=0
      NDIMRO=LSTYPA(2)
      NDIMCO=LSTYPA(2)
      ISTOFF=KSTYPA(1)
      JSTOFF=KSTYPA(1)
      CALL DELK2 (                               WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             IRODIM , IEO    , ITAU   , LNTSYM ,
     4             LNTASY , ISTOFF , JSTOFF ,
     4             NOPERS       )
C
C*********************************************************************
C
C  SET UP THE ODD K CORNER OF THE HAMILTONIAN MATRIX BLOCK
C
C*********************************************************************
C
400   CONTINUE
      IF (IROTO .EQ. 0) GOTO 420
      IEO=1
      NOFSET=IROTE*V2MXP1*NDIMST
      NDIMST=LSTYPA(1)
      ISTOFF=0
      CALL DELK0 ( WORK(IENRGY) , WORK(IESTST) , WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSET , NDIMST , JP1    , IRODIM , IEO    ,
     3             ITAU   , LNTSYM , LNTASY , ISTOFF ,
     4             NOPERS       )
      NOFSRO=NOFSET
      NOFSCO=NOFSET
      NDIMRO=LSTYPA(1)
      NDIMCO=LSTYPA(1)
            ISTOFF=0
            JSTOFF=0
      CALL DELK2 (                               WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             IRODIM , IEO    , ITAU   , LNTSYM ,
     4             LNTASY , ISTOFF , JSTOFF ,
     4             NOPERS       )
C
C*********************************************************************
C
C  SET UP DELTA K = 1 MATRIX ELEMENTS CONNECTING THE TWO CORNERS
C
C*********************************************************************
C
420   CONTINUE
      IF (IROTE .EQ. 0 .OR. IROTO .EQ. 0) GOTO 440
      IEO=0
      NOFSRO=0
      NOFSCO=NOFSET
      NDIMRO=LSTYPA(2)
      NDIMCO=LSTYPA(1)
      ISTOFF=KSTYPA(1)
      JSTOFF=0
      CALL DELK1 (                               WORK(ISYMST) ,
     1             WORK(IASYST) , WORK(ISMAST) , WORK(IHMAST) ,
     2             NOFSRO , NOFSCO , NDIMRO , NDIMCO , JP1    ,
     3             IRODIM , IEO    , ITAU   , LNTSYM ,
     4             LNTASY , ISTOFF , JSTOFF ,
     4             NOPERS       )
440   CONTINUE
C
C*********************************************************************
C
      IF (IRODIM .EQ. 0) GOTO 460
      IEIGST=IWRKEN
      IEIGVC=IEIGST+IRODIM
      IEWKST=IEIGVC+IRODIM*IRODIM
      IQNUMS=IEWKST+8*IRODIM
      IIWORK=IQNUMS+IRODIM
      IIFAIL=IIWORK+5*IRODIM
      IWRKEN=IIFAIL+IRODIM
      IF (IWRKEN .GE. LWORK) GOTO 1011
C
      CALL DIAROT ( IRODIM , WORK(IHMAST) , WORK(IEIGST) , 
     1 WORK(IEIGVC) ,
     1 WORK(IEWKST) , WORK(IQNUMS) ,
     1 J  ,         ITAU   , 2   ,
     1 WORK(IIWORK) , WORK(IIFAIL) )
C
      IF (SYMM) THEN
      FACT=FACJ*GNS(ITAU+1,2,I)
      ELSE
      FACT=FACJ*GNS(1,1,I)
      ENDIF
      DO 168 II=1,IRODIM
      IRECN=(4*MOD(J,2)+2*ITAU+1)*LRODIM+II
      EOUT=WORK(IEIGST+II-1)-E00
      PRTFCT=PRTFCT+FACT*EXP(EOUT*TEMPCO)
168   WRITE (NFIL17,REC=IRECN) EOUT,
     .    (WORK(IHMAST+(II-1)*IRODIM+JJ-1), JJ=1,IRODIM)
      DO 170 II=1,IRODIM
      WORK(ISTORE+INDEXJ)=WORK(IEIGST+II-1)
170   INDEXJ=INDEXJ+1
      JDIMST(INDEXD)=IRODIM
460   INDEXD=INDEXD+1
C
C*********************************************************************
C
180   CONTINUE
C
      IF (IOPT .EQ. 1) THEN
C
      WRITE (NFIL6,6500) J
C
      ISTART(1)=ISTORE
      ISTART(2)=ISTART(1)+JDIMST(1)
      ISTART(3)=ISTART(2)+JDIMST(2)
      ISTART(4)=ISTART(3)+JDIMST(3)
C
      DO 182 II=1,4
182   LINDEX(II)=0
C
      II=0
      DO 184 JJ=1,4
      IF (JDIMST(JJ) .NE. 0) THEN
            II=II+1
            LINDEX(II)=JJ
      ENDIF
184   CONTINUE
      NCOL=II
C
      DO 290 II=1,90
290   ELINE(II:II)=' '
      DO 300 II=1,11
300   ELINE(II:II)='-'
      DO 620 JJ=1,NCOL
      IPLLO=12+(JJ-1)*19
      IPLHI=IPLLO+18
      DO 310 II=IPLLO,IPLHI
310   ELINE(II:II)='-'
620   CONTINUE
      WRITE (NFIL6,6300) ELINE
      DO 330 II=1,90
330   ELINE(II:II)=' '
      ELINE(1:1)=':'
      ELINE(11:11)=':'
      DO 640 II=1,NCOL
      JJ=30+(II-1)*19
640   ELINE(JJ:JJ)=':'
      WRITE (NFIL6,6300) ELINE
      DO 350 II=1,NCOL
      JJ=20+(II-1)*19
      KINDEX=LINDEX(II)
      IF (SYMM) THEN
            WRITE (ELINE(JJ:JJ+3),6400) SYMSYM(KINDEX)
      ELSE
            WRITE (ELINE(JJ:JJ+3),6400) ASYSYM(KINDEX)
      ENDIF
350   CONTINUE
      WRITE (NFIL6,6300) ELINE
      DO 660 II=1,90
660   ELINE(II:II)=' '
      DO 370 II=1,90
370   ELINE(II:II)=' '
      ELINE(1:1)=':'
      ELINE(11:11)=':'
      DO 680 II=1,NCOL
      JJ=30+(II-1)*19
680   ELINE(JJ:JJ)=':'
      WRITE (NFIL6,6300) ELINE
      DO 390 II=1,11
390   ELINE(II:II)='-'
      DO 410 JJ=1,NCOL
      IPLLO=12+(JJ-1)*19
      IPLHI=IPLLO+18
      DO 700 II=IPLLO,IPLHI
700   ELINE(II:II)='-'
410   CONTINUE
      WRITE (NFIL6,6300) ELINE
C
      IDIMMX=MAX0(JDIMST(1),JDIMST(2),JDIMST(3),JDIMST(4))
C
      DO 412 II=1,90
412   ELINE(II:II)=' '
      DO 190 II=1,IDIMMX
      ELINE(1:1)=':'
      ELINE(11:11)=':'
      KINDEX=II-1
      WRITE (ELINE(4:6),6000) KINDEX
C
      DO 188 JJ=1,NCOL
      KINDEX=LINDEX(JJ)
      IPLLO=15+(JJ-1)*19
      IPLHI=IPLLO+11
      IF (JDIMST(KINDEX) .GE. II) THEN
           EOUT = WORK(ISTART(KINDEX)+II-1) - E00
           WRITE (ELINE(IPLLO:IPLHI),6100) EOUT
      ELSE
           WRITE (ELINE(IPLLO:IPLHI),6200)
      ENDIF
188   ELINE(IPLHI+4:IPLHI+4)=':'
190   WRITE (NFIL6,6300) ELINE
      DO 720 II=1,11
720   ELINE(II:II)='-'
      DO 740 JJ=1,NCOL
      IPLLO=12+(JJ-1)*19
      IPLHI=IPLLO+18
      DO 430 II=IPLLO,IPLHI
430   ELINE(II:II)='-'
740   CONTINUE
      WRITE (NFIL6,6300) ELINE
C
      ENDIF
C
      IF (.NOT. SYMM) THEN
           IDMRZ=V2MXP1*LSTYPA(1)
           IDMCZ=V2MXP1*LSTYPA(1)
           IDMRY=V2MXP1*LSTYPA(1)
      ELSE
           IDMRZ=V2MXP1*LSTYPA(1)
           IDMCZ=V2MXP1*LSTYPA(2)
           IDMRY=V2MXP1*MAX0(LSTYPA(1),LSTYPA(2))
      ENDIF
      ICLOST=IHMAST+LRODIM
      IZMUST=ICLOST+LRODIM
      IYMUST=IZMUST+IDMRZ*IDMCZ
      IWRKEN=IYMUST+IDMRY*IDMRY
      IF (IWRKEN .GT. LWORK) GOTO 1015
C
      CALL CINTNS ( J     , LRODIM , IDMRZ  , IDMCZ  ,
     .              IDMRY ,
     .              WORK(IHMAST) , WORK(ICLOST) , WORK(IZMUST) ,
     .              WORK(IYMUST) , NTRANS       ,
     .              FRQLO        , FRQHI        , XINLIM       ,
     .              V2A , NVA , IABA , V2B , NVB , IABB ,
     .              TRTLIM , ABSINT , GNS , I ,
     .              ENGLIM )
C
200   CONTINUE
      CLOSE(NFIL17)
      CLOSE(NFIL18)
      CLOSE(NFIL19)
      REWIND NFIL20
C
      ISPACE=INT(FLOAT(LWORK-IWRKST)/6.0D+00)
      NACC=MIN0(NTRANS,ISPACE)
      IF (NTRANS .GT. ISPACE .AND. IOPT .LT. 3) THEN
           WRITE (NFIL6,2200) ISPACE,NTRANS
2200       FORMAT(1H0,'MORBID.CNT.WRN  FIRST ',I8,' TRANSITIONS ',
     .     'PRINTED OUT',/,1H ,'(TOTAL NUMBER: ',I12,')',//)
      ENDIF
      INDEX=IWRKST-1
      DO 4000 JJ=1,NACC
      READ (NFIL20) (WORK(INDEX+II), II=1,6)
4000  INDEX=INDEX+6
C
      IF (IOPT .LT. 4) THEN
C
C*********************************************************************
C
C  SET UP ASCENDING BUBBLE SORT OF THE CALCULATED FREQUENCIES
C
C*********************************************************************
C
      NM1=NACC-1
      DO 4200 II=1,NM1
      IOFF=6*(II-1)
      KK=II
      P=WORK(IWRKST+IOFF+3)
      IP1=II+1
      DO 4100 JJ=IP1,NACC
      JOFF=6*(JJ-1)
      IF (WORK(IWRKST+JOFF+3) .GT. P) GOTO 4100
      KK=JJ
      P=WORK(IWRKST+JOFF+3)
4100  CONTINUE
      IF (KK .EQ. II) GOTO 4200
      KOFF=6*(KK-1)
      DO 4150 JJ=1,6
      P=WORK(IWRKST+KOFF+JJ-1)
      WORK(IWRKST+KOFF+JJ-1)=WORK(IWRKST+IOFF+JJ-1)
4150  WORK(IWRKST+IOFF+JJ-1)=P
4200  CONTINUE
      ENDIF
C
C**********************************************************
C
C     PRINT OUT CALCULATED TRANSITIONS
C
C**********************************************************
C
      IF (IOPT .LT. 3) THEN
      WRITE (NFIL6,5300)
C
      NPROUT=0
      DO 4300 II=1,NACC
      IOFF=6*(II-1)
      XPACKL=WORK(IWRKST+IOFF  )
      XPACKU=WORK(IWRKST+IOFF+1)
      EL    =WORK(IWRKST+IOFF+2)
      FREQ  =WORK(IWRKST+IOFF+3)
      TRMOM =WORK(IWRKST+IOFF+4)
      EU=EL+FREQ
C
      XZ=XPACKL/1.0D+09+ABIT
      JL=INT(XZ)
      XPACKL=XPACKL-JL*1.0D+09
      XZ=XPACKL/1.0D+08+ABIT
      TAUL=INT(XZ)
      XPACKL=XPACKL-TAUL*1.0D+08
      XZ=XPACKL/1.0D+07+ABIT
      IEOL=INT(XZ)
      XPACKL=XPACKL-IEOL*1.0D+07
      XZ=XPACKL/1.0D+05+ABIT
      KVALL=INT(XZ)
      XPACKL=XPACKL-KVALL*1.0D+05
      XZ=XPACKL/1.0D+03+ABIT
      NVALL=INT(XZ)
      XPACKL=XPACKL-NVALL*1.0D+03
      XZ=XPACKL/1.0D+01+ABIT
      V2VALL=INT(XZ)
      XPACKL=XPACKL-V2VALL*1.0D+01
      IABL=NINT(XPACKL)
C
      XZ=XPACKU/1.0D+09+ABIT
      JU=INT(XZ)
      XPACKU=XPACKU-JU*1.0D+09+ABIT
      XZ=XPACKU/1.0D+08
      TAUU=INT(XZ)
      XPACKU=XPACKU-TAUU*1.0D+08
      XZ=XPACKU/1.0D+07+ABIT
      IEOU=INT(XZ)
      XPACKU=XPACKU-IEOU*1.0D+07
      XZ=XPACKU/1.0D+05+ABIT
      KVALU=INT(XZ)
      XPACKU=XPACKU-KVALU*1.0D+05
      XZ=XPACKU/1.0D+03+ABIT
      NVALU=INT(XZ)
      XPACKU=XPACKU-NVALU*1.0D+03
      XZ=XPACKU/1.0D+01+ABIT
      V2VALU=INT(XZ)
      XPACKU=XPACKU-V2VALU*1.0D+01
      IABU=NINT(XPACKU)
C
      IF (SYMM) THEN
      TRINT=ABSINT*FREQ*GNS(TAUL+1,IEOL+1,I)*EXP(EL*TEMPCO)*
     .      (1.0D+00-EXP(FREQ*TEMPCO))*TRMOM/PRTFCT
      ELSE
      TRINT=ABSINT*FREQ*GNS(1,1,I)*EXP(EL*TEMPCO)*
     .      (1.0D+00-EXP(FREQ*TEMPCO))*TRMOM/PRTFCT
      ENDIF
C
      IF (TRINT .GE. TRTLIM) THEN
      NPROUT=NPROUT+1
      IF (SYMM) THEN
      WRITE (NFIL6,5400) JU,SYMSYM(IEOU+1+TAUU*2),KVALU,V2VALU,
     .      SYMSYM(IABU),NVALU,EU
      WRITE (NFIL6,5450) JL,SYMSYM(IEOL+1+TAUL*2),KVALL,V2VALL,
     .      SYMSYM(IABL),NVALL,EL,FREQ,TRMOM,TRINT
      ELSE
      WRITE (NFIL6,5400) JU,ASYSYM(TAUU*2+1),KVALU,V2VALU,
     .      ASYSYM(1),NVALU,EU
      WRITE (NFIL6,5450) JL,ASYSYM(TAUL*2+1),KVALL,V2VALL,
     .      ASYSYM(1),NVALL,EL,FREQ,TRMOM,TRINT
      ENDIF
      ENDIF
4300  CONTINUE
C
      WRITE (NFIL6,5460) XINLIM,NTRANS,NACC,TRTLIM,NPROUT
      GOTO 10
      ENDIF
C
      IF (IOPT .EQ. 3) WRITE (NFIL6,3333) FRQLO,FRQHI
      DO 4320 II=1,NACC
      IOFF=6*(II-1)
      XPACKL=WORK(IWRKST+IOFF  )
      XPACKU=WORK(IWRKST+IOFF+1)
      EL    =WORK(IWRKST+IOFF+2)
      FREQ  =WORK(IWRKST+IOFF+3)
      TRMOM =WORK(IWRKST+IOFF+4)
      EU=EL+FREQ
C
      XZ=XPACKL/1.0D+09+ABIT
      JL=INT(XZ)
      XPACKL=XPACKL-JL*1.0D+09
      XZ=XPACKL/1.0D+08+ABIT
      TAUL=INT(XZ)
      XPACKL=XPACKL-TAUL*1.0D+08
      XZ=XPACKL/1.0D+07+ABIT
      IEOL=INT(XZ)
      XPACKL=XPACKL-IEOL*1.0D+07
      XZ=XPACKL/1.0D+05+ABIT
      KVALL=INT(XZ)
      XPACKL=XPACKL-KVALL*1.0D+05
      XZ=XPACKL/1.0D+03+ABIT
      NVALL=INT(XZ)
      XPACKL=XPACKL-NVALL*1.0D+03
      XZ=XPACKL/1.0D+01+ABIT
      V2VALL=INT(XZ)
      XPACKL=XPACKL-V2VALL*1.0D+01
      IABL=NINT(XPACKL)
C
      XZ=XPACKU/1.0D+09+ABIT
      JU=INT(XZ)
      XPACKU=XPACKU-JU*1.0D+09
      XZ=XPACKU/1.0D+08+ABIT
      TAUU=INT(XZ)
      XPACKU=XPACKU-TAUU*1.0D+08
      XZ=XPACKU/1.0D+07+ABIT
      IEOU=INT(XZ)
      XPACKU=XPACKU-IEOU*1.0D+07
      XZ=XPACKU/1.0D+05+ABIT
      KVALU=INT(XZ)
      XPACKU=XPACKU-KVALU*1.0D+05
      XZ=XPACKU/1.0D+03+ABIT
      NVALU=INT(XZ)
      XPACKU=XPACKU-NVALU*1.0D+03
      XZ=XPACKU/1.0D+01+ABIT
      V2VALU=INT(XZ)
      XPACKU=XPACKU-V2VALU*1.0D+01
      IABU=NINT(XPACKU)
C
      IF (SYMM) THEN
      TRINT=ABSINT*FREQ*GNS(TAUL+1,IEOL+1,I)*EXP(EL*TEMPCO)*
     .      (1.0D+00-EXP(FREQ*TEMPCO))*TRMOM/PRTFCT
      ELSE
      TRINT=ABSINT*FREQ*GNS(1,1,I)*EXP(EL*TEMPCO)*
     .      (1.0D+00-EXP(FREQ*TEMPCO))*TRMOM/PRTFCT
      ENDIF
      IF (IOPT .EQ. 3 .AND. TRINT .GT. TRTLIM) WRITE (NFIL6,3333)
     1                            FREQ,TRINT
      WORK(IWRKST+IOFF+5)=TRINT
4320  CONTINUE
C
      IF (IOPT .EQ. 3) GOTO 10
C
C
C*********************************************************************
C
C  FIND THE 100 STRONGEST TRANSITIONS
C
C*********************************************************************
C
      NM1=MIN0(NACC-1,100)
      DO 4222 II=1,NM1
      IOFF=6*(II-1)
      KK=II
      P=WORK(IWRKST+IOFF+5)
      IP1=II+1
      DO 4122 JJ=IP1,NACC
      JOFF=6*(JJ-1)
      IF (WORK(IWRKST+JOFF+5) .LT. P) GOTO 4122
      KK=JJ
      P=WORK(IWRKST+JOFF+5)
4122  CONTINUE
      IF (KK .EQ. II) GOTO 4222
      KOFF=6*(KK-1)
      DO 4152 JJ=1,6
      P=WORK(IWRKST+KOFF+JJ-1)
      WORK(IWRKST+KOFF+JJ-1)=WORK(IWRKST+IOFF+JJ-1)
4152  WORK(IWRKST+IOFF+JJ-1)=P
4222  CONTINUE
C
C*********************************************************************
C
C  PRINT OUT THE STRONGEST TRANSITION
C
C*********************************************************************
C
      XPACKL=WORK(IWRKST  )
      XPACKU=WORK(IWRKST+1)
      EL    =WORK(IWRKST+2)
      FREQ  =WORK(IWRKST+3)
      TRMOM =WORK(IWRKST+4)
      TRINT =WORK(IWRKST+5)
      EU=EL+FREQ
C
      XZ=XPACKL/1.0D+09+ABIT
      JL=INT(XZ)
      XPACKL=XPACKL-JL*1.0D+09
      XZ=XPACKL/1.0D+08+ABIT
      TAUL=INT(XZ)
      XPACKL=XPACKL-TAUL*1.0D+08
      XZ=XPACKL/1.0D+07+ABIT
      IEOL=INT(XZ)
      XPACKL=XPACKL-IEOL*1.0D+07
      XZ=XPACKL/1.0D+05+ABIT
      KVALL=INT(XZ)
      XPACKL=XPACKL-KVALL*1.0D+05
      XZ=XPACKL/1.0D+03+ABIT
      NVALL=INT(XZ)
      XPACKL=XPACKL-NVALL*1.0D+03
      XZ=XPACKL/1.0D+01+ABIT
      V2VALL=INT(XZ)
      XPACKL=XPACKL-V2VALL*1.0D+01
      IABL=NINT(XPACKL)
C
      XZ=XPACKU/1.0D+09+ABIT
      JU=INT(XZ)
      XPACKU=XPACKU-JU*1.0D+09
      XZ=XPACKU/1.0D+08+ABIT
      TAUU=INT(XZ)
      XPACKU=XPACKU-TAUU*1.0D+08
      XZ=XPACKU/1.0D+07+ABIT
      IEOU=INT(XZ)
      XPACKU=XPACKU-IEOU*1.0D+07
      XZ=XPACKU/1.0D+05+ABIT
      KVALU=INT(XZ)
      XPACKU=XPACKU-KVALU*1.0D+05
      XZ=XPACKU/1.0D+03+ABIT
      NVALU=INT(XZ)
      XPACKU=XPACKU-NVALU*1.0D+03
      XZ=XPACKU/1.0D+01+ABIT
      V2VALU=INT(XZ)
      XPACKU=XPACKU-V2VALU*1.0D+01
      IABU=NINT(XPACKU)
C
      WRITE (NFIL6,5399)
      IF (SYMM) THEN
      WRITE (NFIL6,5400) JU,SYMSYM(IEOU+1+TAUU*2),KVALU,V2VALU,
     .      SYMSYM(IABU),NVALU,EU
      WRITE (NFIL6,5450) JL,SYMSYM(IEOL+1+TAUL*2),KVALL,V2VALL,
     .      SYMSYM(IABL),NVALL,EL,FREQ,TRMOM,TRINT
      ELSE
      WRITE (NFIL6,5400) JU,ASYSYM(TAUU*2+1),KVALU,V2VALU,
     .      ASYSYM(1),NVALU,EU
      WRITE (NFIL6,5450) JL,ASYSYM(TAUL*2+1),KVALL,V2VALL,
     .      ASYSYM(1),NVALL,EL,FREQ,TRMOM,TRINT
      ENDIF
      WRITE (NFIL6,5398)
C
C*********************************************************************
C
C  SORT THE HUNDRED STRONGEST TRANSITIONS ACCORDING TO FREQUENCY
C
C*********************************************************************
C
      NM1P1=NM1
      NM1=NM1-1
      DO 5222 II=1,NM1
      IOFF=6*(II-1)
      KK=II
      P=WORK(IWRKST+IOFF+3)
      IP1=II+1
      DO 5122 JJ=IP1,NM1P1
      JOFF=6*(JJ-1)
      IF (WORK(IWRKST+JOFF+3) .GT. P) GOTO 5122
      KK=JJ
      P=WORK(IWRKST+JOFF+3)
5122  CONTINUE
      IF (KK .EQ. II) GOTO 5222
      KOFF=6*(KK-1)
      DO 5152 JJ=1,6
      P=WORK(IWRKST+KOFF+JJ-1)
      WORK(IWRKST+KOFF+JJ-1)=WORK(IWRKST+IOFF+JJ-1)
5152  WORK(IWRKST+IOFF+JJ-1)=P
5222  CONTINUE
C
C*********************************************************************
C
C  MAKE A TEX INPUT TABLE OF THE 100 STRONGEST TRANSITIONS
C
C*********************************************************************
C
C
      DO 6222 II=1,NM1+1
      IOFF=6*(II-1)
      XPACKL=WORK(IWRKST+IOFF  )
      XPACKU=WORK(IWRKST+IOFF+1)
      EL    =WORK(IWRKST+IOFF+2)
      FREQ  =WORK(IWRKST+IOFF+3)
      TRMOM =WORK(IWRKST+IOFF+4)
      TRINT =WORK(IWRKST+IOFF+5)
      EU=EL+FREQ
C
      XZ=XPACKL/1.0D+09+ABIT
      JL=INT(XZ)
      XPACKL=XPACKL-JL*1.0D+09
      XZ=XPACKL/1.0D+08+ABIT
      TAUL=INT(XZ)
      XPACKL=XPACKL-TAUL*1.0D+08
      XZ=XPACKL/1.0D+07+ABIT
      IEOL=INT(XZ)
      XPACKL=XPACKL-IEOL*1.0D+07
      XZ=XPACKL/1.0D+05+ABIT
      KVALL=INT(XZ)
      XPACKL=XPACKL-KVALL*1.0D+05
      XZ=XPACKL/1.0D+03+ABIT
      NVALL=INT(XZ)
      XPACKL=XPACKL-NVALL*1.0D+03
      XZ=XPACKL/1.0D+01+ABIT
      V2VALL=INT(XZ)
      XPACKL=XPACKL-V2VALL*1.0D+01
      IABL=NINT(XPACKL)
C
      XZ=XPACKU/1.0D+09+ABIT
      JU=INT(XZ)
      XPACKU=XPACKU-JU*1.0D+09
      XZ=XPACKU/1.0D+08+ABIT
      TAUU=INT(XZ)
      XPACKU=XPACKU-TAUU*1.0D+08
      XZ=XPACKU/1.0D+07+ABIT
      IEOU=INT(XZ)
      XPACKU=XPACKU-IEOU*1.0D+07
      XZ=XPACKU/1.0D+05+ABIT
      KVALU=INT(XZ)
      XPACKU=XPACKU-KVALU*1.0D+05
      XZ=XPACKU/1.0D+03+ABIT
      NVALU=INT(XZ)
      XPACKU=XPACKU-NVALU*1.0D+03
      XZ=XPACKU/1.0D+01+ABIT
      V2VALU=INT(XZ)
      XPACKU=XPACKU-V2VALU*1.0D+01
      IABU=NINT(XPACKU)
      IF (SYMM) THEN
      WRITE (NFIL6,5444) JU,KVALU,SYMSYM(IEOU+1+TAUU*2),
     .                   JL,KVALL,SYMSYM(IEOL+1+TAUL*2),
     .                   EL,FREQ,TRMOM,TRINT
      ELSE
      WRITE (NFIL6,5444) JU,KVALU,ASYSYM(TAUU*2+1),
     .                   JL,KVALL,ASYSYM(TAUL*2+1),
     .                   EL,FREQ,TRMOM,TRINT
      ENDIF
6222  CONTINUE
C
3333  FORMAT(' ',F18.6,' ',E18.6)
C
10    CONTINUE
C
5300  FORMAT('1',5X,12('*'),' CALCULATED TRANSITIONS ',12('*')//
     1            '0',18X,'J',1X,'SYM.',1X,'KV',1X,'V2',1X,
     1            'V.SYM.',3X,
     1      'N',10X,'ENERGY',8X,'NUCALC',6X,'LINE STRENGTH',
     1      3X,'INTENSITY'//'0',54X,'-1',12X,'-1',16X,'2',
     1      10X,'-1'/
     1           ' ',52X,'CM',12X,'CM',13X,'DEBYE',8X,
     1      'MOL  CM'//' ',128('-')/)
5398  FORMAT(' ',//)
5399  FORMAT('1',5X,12('*'),' STRONGEST TRANSITION ',12('*')//)
5400  FORMAT('0 UPPER STATE: ',3X,I2,1X,A4,1X,I2,1X,I2,2X,A4,
     1      2X,I2,2X,F15.5)
5450  FORMAT('  LOWER STATE: ',3X,I2,1X,A4,1X,I2,1X,I2,2X,A4,
     1          2X,I2,2X,2F15.5,2E15.5)
5444  FORMAT(' ',I3,' &',I3,' & ',A4,' & \\leftarrow & ',
     .           I3,' &',I3,' & ',A4,' & ',F12.4,/,
     .       '  & ',F12.4,
     .           ' & ????.????',2(' & ',E12.3),' \\cx')
5460  FORMAT('0',4X,'NUMBER OF TRANSITIONS CALCULATED ',
     1            'WITH A LINE STRENGTH LARGER THAN ',
     1            E12.4,' DEBYE**2   : ',I6,/,
     1       '0',4X,'NUMBER OF TRANSITIONS RETAINED FO',
     1            'R SORTING                        ',
     1             12X ,'            : ',I6,/,
     1       '0',4X,'NUMBER OF TRANSITIONS CALCULATED ',
     1            'WITH AN   INTENSITY  LARGER THAN ',
     1            E12.4,' CM/MOL     : ',I6)
6000  FORMAT(I3)
6100  FORMAT(F12.5)
6200  FORMAT(12X)
6300  FORMAT(1H ,4X,A90)
6400  FORMAT(A4)
6500  FORMAT(1H0,12('*'),'  J = ',I3,' ENERGIES  ',12('*'),/)
6600  FORMAT(1H0,12('*'),'  ZERO POINT ENERGY  ',12('*'),//,
     1      12X,'E0 = ',F12.5//)
      CALL PRTITO
      STOP
1001  WRITE (NFIL6,2001)
2001  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' INITIAL FUNCTION CALCULATION')
      STOP
1002  WRITE (NFIL6,2002)
2002  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' CALCULATION OF F1, F2 ETC.')
      STOP
1003  WRITE (NFIL6,2003)
2003  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' NUMEROV-COOLEY INTEGRATION')
      STOP
1004  WRITE (NFIL6,2004)
2004  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' BENDING INTEGRAL CALCULATION (J=0)')
      STOP
1005  WRITE (NFIL6,2005)
2005  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' SETTING UP MORSE OSCILLATOR PARAMETERS')
      STOP
1006  WRITE (NFIL6,2006)
2006  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' DIAGONALIZING THE STRETCHING MATRIX')
      STOP
1007  WRITE (NFIL6,2007)
2007  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' CALCULATING STRETCHING MATRIX ELEMENTS')
      STOP
1009  WRITE (NFIL6,2009) IERR2,IERR1
2009  FORMAT(1H0,'MORBID.CNT.ERR  SCRATCH FILE COULD NOT BE O',
     1          'PENED FOR UNIT ',I3,/,1H ,'IOSTAT = ',I5)
      STOP
1010  WRITE (NFIL6,2010) J
2010  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' ROTATIONAL CALCULATION AT J = ',I3)
      STOP
1011  WRITE (NFIL6,2011) J
2011  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' DIAGONALIZATION AT J = ',I3)
      STOP
1012  WRITE (NFIL6,2012)
2012  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE TO',
     1          ' CALCULATE MUZ MATRIX ELEMENTS')
      STOP
1014  WRITE (NFIL6,2014)
2014  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE TO',
     1          ' CALCULATE MUY MATRIX ELEMENTS')
      STOP
1015  WRITE (NFIL6,2010) J
2015  FORMAT(1H0,'MORBID.CNT.ERR  INSUFFICIENT WORK SPACE FOR',
     1          ' INTENSITY CALCULATION AT J = ',I3)
      END
C
C
      INTEGER FUNCTION ILCONV (I,N)
      ILCONV=N*(I-1)+1
      RETURN
      END
C
C
      SUBROUTINE PRTIME( ROUTIN , OLDTIM , OLDVEC )
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION OLDTIM,NEWTIM,DELTIM,TIM1S,TIM2S
      CHARACTER*6 ROUTIN
C
C
      COMMON /SYS/    THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1                     THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2                     VELLGT , PLANCK , AVOGNO , PI , PREC ,
     4                     NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 ,
     5                     NFIL6 , NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6                     NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7                     NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6                     ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7                     IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8                     NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER ,
     9                       NATTS , PRTINT
C
      IF (IPRINT .EQ. 0) RETURN
      CALL CLOCKV ( VECTIM , NEWTIM , 1 , 2 )
      DELTIM = NEWTIM - OLDTIM
      DELVEC = VECTIM - OLDVEC
      TIM1S=DELTIM/1.0D+03
      TIM2S=NEWTIM/1.0D+03
      TIM3S=DELVEC/1.0D+03
      TIM4S=VECTIM/1.0D+03
      WRITE (NFIL6,2007) ROUTIN, TIM1S, TIM3S, TIM2S, TIM4S
2007  FORMAT(1H0,'INTENS.TIM.INF  FOR ROUTINE ',A6,' : '/
     1   1H ,'                CPU TIME           ',F20.8,' SECONDS'/
     1   1H ,'                VECTOR TIME        ',F20.8,' SECONDS'/
     1   1H ,'                TOTAL CPU TIME     ',F20.8,' SECONDS'/
     1   1H ,'                TOTAL VECTOR TIME  ',F20.8,' SECONDS')
      OLDTIM=NEWTIM
      OLDVEC=VECTIM
      RETURN
      END
C
C
      SUBROUTINE PRTIMJ( J , OLDTIM , OLDVEC )
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION OLDTIM,NEWTIM,DELTIM,TIM1S,TIM2S
C
C
      COMMON /SYS/    THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1                     THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2                     VELLGT , PLANCK , AVOGNO , PI , PREC ,
     4                     NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 ,
     5                     NFIL6 , NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6                     NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7                     NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6                     ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7                     IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8                     NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER ,
     9                       NATTS , PRTINT
C
      IF (IPRINT .EQ. 0) RETURN
      CALL CLOCKV ( VECTIM , NEWTIM , 1 , 2 )
      DELTIM = NEWTIM - OLDTIM
      DELVEC = VECTIM - OLDVEC
      TIM1S=DELTIM/1.0D+03
      TIM2S=NEWTIM/1.0D+03
      TIM3S=DELVEC/1.0D+03
      TIM4S=VECTIM/1.0D+03
      WRITE (NFIL6,2007) J, TIM1S, TIM3S, TIM2S, TIM4S
2007  FORMAT(1H0,'INTENS.TIM.INF  FOR CALCULATING J = ',I3,' : '/
     1   1H ,'                CPU TIME           ',F20.8,' SECONDS'/
     1   1H ,'                VECTOR TIME        ',F20.8,' SECONDS'/
     1   1H ,'                TOTAL CPU TIME     ',F20.8,' SECONDS'/
     1   1H ,'                TOTAL VECTOR TIME  ',F20.8,' SECONDS')
      OLDTIM=NEWTIM
      OLDVEC=VECTIM
      RETURN
      END
C
C
      SUBROUTINE PRTITO
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION NEWTIM,TIM1S
      CHARACTER*6 ROUTIN
C
C
      COMMON /SYS/    THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1                     THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2                     VELLGT , PLANCK , AVOGNO , PI , PREC ,
     4                     NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 ,
     5                     NFIL6 , NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6                     NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7                     NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6                     ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7                     IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8                     NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER ,
     9                       NATTS , PRTINT
C
      CALL CLOCKV ( VECTIM , NEWTIM , 1 , 2 )
      TIM2S=NEWTIM/1.0D+03
      TIM4S=VECTIM/1.0D+03
      WRITE (NFIL6,2007) TIM2S, TIM4S
2007  FORMAT(1H0,'INTENS.CNT.INF  PROGRAM TERMINATED NORMALLY ',/,
     1       1H ,'                ENTIRE CALCULATION USED : '/
     1   1H ,'                TOTAL CPU TIME     ',F20.8,' SECONDS'/
     1   1H ,'                TOTAL VECTOR TIME  ',F20.8,' SECONDS')
      RETURN
      END
C
C
      BLOCK DATA INIVAL
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL SYMM
      REAL*8 M1,M2,M3,M
      INTEGER ISOMAX , IVAR(55) , PARMAX
C
C
      COMMON /MOLCUL/ RHOE,FA1,FA2,FA3,FA4,FA5,FA6,FA7,FA8,
     1               AA1,AA3,
     2               F1A1,F2A1,F3A1,F4A1,F1A3,F2A3,F3A3,F4A3,
     3               F11,F1A11,F2A11,F3A11,F33,F1A33,F2A33,F3A33,
     4               F13,F1A13,F2A13,F3A13,
     5               F111,F1A111,F2A111,F333,F1A333,F2A333,
     6               F113,F1A113,F2A113,F133,F1A133,F2A133,
     7               F1111,FA1111,F3333,FA3333,F1113,FA1113,
     8               F1333,FA1333,F1133,FA1133,
     8               RE12 , RE32 , M1 , M2 , M3 , M ,
     9               U1 , U3 , U13 , V ,
     1               SYMM
C
      COMMON /BMSC/  BOLTZ
C
      COMMON /SYS/    THRSH1 , THRSH2 , THRSH3 , THRSH4 , THRSH5 ,
     1                     THRSH6 , THRSH7 , THRSH8 , THRSH9 , THRSHX ,
     2                     VELLGT , PLANCK , AVOGNO , PI , PREC ,
     4                     NFIL1 , NFIL2 , NFIL3 , NFIL4 , NFIL5 ,
     5                     NFIL6 , NFIL7 , NFIL8 , NFIL9 , NFIL10 ,
     6                     NFIL11 , NFIL12 , NFIL13 , NFIL14 , NFIL15 ,
     7                     NFIL16 , NFIL17 , NFIL18 , NFIL19 , NFIL20 ,
     6                     ITEST  , IPRINT , NSTNR , NSTNIN , IREST ,
     7                     IISOT , IQUAS , ISYMS , NISOT , NQUAS ,
     8                     NUMQUA , NOPTIT , NOPTIM , IOBSER , NOBSER ,
     9                       NATTS , PRTINT
C
      COMMON /LSFIT/  PARMAX , NUMPAR , ISOMAX , IVAR
C
      DATA SYMM/.FALSE./
      DATA VELLGT/2.99792458E+10/,PLANCK/6.6260755E-27/,
     1    AVOGNO/6.0221367E+23/, PI/3.14159265359D+00/,
     2    PREC/0.5D-30/,BOLTZ/1.380658D-16/
      DATA NFIL1/1/,NFIL2/2/,NFIL3/3/,NFIL4/4/,NFIL5/5/,
     1    NFIL6/6/,NFIL7/7/,NFIL8/8/,NFIL9/9/,NFIL10/10/,
     2    NFIL11/11/,NFIL12/12/,NFIL13/13/,NFIL14/14/,NFIL15/15/,
     3    NFIL16/16/,NFIL17/17/,NFIL18/18/,NFIL19/19/,NFIL20/20/
      DATA ITEST/0/,IPRINT/1/,ISOMAX/10/
      END
