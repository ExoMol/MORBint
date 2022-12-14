      SUBROUTINE DIAROT(IDIMTO,HMAT,E,EIGVEC,WSPACE,QNUMBS,JU,ITAU,
     1                  IDENT,IWORK,IFAIL)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 HMAT(IDIMTO,IDIMTO)
      REAL*8 EIGVEC(IDIMTO,IDIMTO),E(IDIMTO)
      REAL*8 WSPACE(8*IDIMTO),QNUMBS(IDIMTO)
      INTEGER IWORK(5*IDIMTO)
      INTEGER IFAIL(IDIMTO)
C
      CHARACTER*1 JOBZ,RANGE,UPLO
      REAL*8 VL,VU
      INTEGER IL,IU
      REAL*8 ABSTOL
      INTEGER MVAL
      INTEGER INFO
C
      CHARACTER*4 SYMSYM(2,2),ASYSYM(2),PRTSYM
      REAL*8 M1,M2,M3,M,U1,U3,U13,V,RHO,EPS,
     1      CR,SR,CSE,SNE,
     2      CRE,SRE,CORO,EPSP,EPSPP,EPSPPP
C
      REAL*8 AMASS(3,10)
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
     1       NUMPAR , PRTINT
C
      INTEGER V1 ,V2, V3 , V2MXP1, V2P1,
     1       NSTINT , NSERIN , NSERP , NSERQ , KQUA , NTEST ,
     2       NSEPP2 , NSEQP1 , MBASIS ,
     3       MDIM , NFSYM0, NFASY0, NFSYMJ, NFASYJ,
     4       KSTYPA(2) , LSTYPA(2) , JMAX , V2MAX , JMAXP1
C
      INTEGER IQUANT(5,10)
C
      LOGICAL SYMM
C
      REAL*8 RMK(2),AAS(2)
      INTEGER NOBAS,MBASP1,LENIW
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
      COMMON /INTEG/  ETRIAL , RHOMAX , PNM1 , HBASE , HSTEP , EGUESS ,
     1               RHOREF , VMIN , V0TYPE ,
     1               NSTINT , NSERIN , NSERP , NSERQ , KQUA , NTEST ,
     2               NSEPP2 , NSEQP1 , KSTYPA , LSTYPA
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
C
      COMMON/MORSE/RMK,AAS
      COMMON/MODIM/NOBAS,MBASP1,LENIW
C
      DATA SYMSYM/'A1  ','B1  ','B2  ','A2  '/
      DATA ASYSYM/'A''  ','A"  '/
C
      IDMM1=IDIMTO-1
      DO 200 II=1,IDMM1
      IIP1=II+1
      DO 100 JJ=IIP1,IDIMTO
100   HMAT(JJ,II)=HMAT(II,JJ)
200   CONTINUE
C
      IF (ITEST .GT. 0) THEN
            WRITE (NFIL6,9000)
            DO 210 II=1,IDIMTO
            WRITE (NFIL6,9010)
            WRITE (NFIL6,9020) (HMAT(JJ,II),JJ=1,IDIMTO)
210         CONTINUE
      ENDIF
C
        JOBZ='V'
        RANGE='A'
        UPLO='U'
        VL=0.0
        VU=0.0
        IL=1
        IU=1
        ABSTOL=0
        MVAL=0
        CALL DSYEVX(JOBZ,RANGE,UPLO,
     1              IDIMTO,HMAT,IDIMTO,
     2              VL,VU,IL,IU,
     3              ABSTOL,
     4              MVAL,E,EIGVEC,IDIMTO,
     5              WSPACE,8*IDIMTO,IWORK,
     6              IFAIL,INFO)
        IF (INFO .NE. 0) THEN
          WRITE (NFIL6,6100)
          STOP
        ENDIF
C
      DO 300 II=1,IDIMTO
      DO 300 JJ=1,IDIMTO
300   HMAT(JJ,II)=EIGVEC(JJ,II)
C
      RETURN
6100  FORMAT(1H0,'  MORBID.MRO.ERR  MATRIX DIAGONALIZATION  ',
     1             'FAILED IN DIAROT'/)
9000  FORMAT('0','  SYMMETRIZED MATRIX BLOCK FROM DIAROT',//)
9010  FORMAT('0')
9020  FORMAT(1H ,5F12.3)
      END
