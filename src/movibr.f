      SUBROUTINE VIBRAM ( CHX11 , CHX13 , CHX31 , CHX33 ,
     1                   CHR11 , CHR13 , CHR31 , CHR33 ,
     2                   PX1   , PX3   , PR1   , PR3   )
      IMPLICIT REAL*8 (A-H,O-Z)
C            :
C DATE       : 29.04.1986
C AUTHOR     : PER JENSEN
C UPDATES    :
C LANGUAGE   : FORTRAN
C PURPOSE    : CALCULATE EXPANSION COEFFICIENTS OF THE TWO VIBRATIONAL
C            : ANGULAR MOMENTA, P SUB X AND P SUB RHO.
C            : THE CODE IS GENERATED BY REDUCE.
C            :
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
      REAL*8 CHX11,CHX13,CHX31,CHX33,CHR11,CHR31,CHR13,CHR33
      REAL*8 PX1(15),PX3(15),PR1(15),PR3(15)
      REAL*8 PX11,PX13,PX111,PX133,PX113,PX1111,PX1333,PX1113,PX1133,
     1       PX11111,PX13333,PX11113,PX11333,PX11133,
     2       PX31,PX33,PX311,PX333,PX313,PX3111,PX3333,PX3113,PX3133,
     3       PX31111,PX33333,PX31113,PX31333,PX31133,
     4       PR11,PR13,PR111,PR133,PR113,PR1111,PR1333,PR1113,PR1133,
     5       PR11111,PR13333,PR11113,PR11333,PR11133,
     6       PR31,PR33,PR311,PR333,PR313,PR3111,PR3333,PR3113,PR3133,
     7       PR31111,PR33333,PR31113,PR31333,PR31133
      REAL*8 A11,A13,A31,A33,A111,A133,A311,A333,A113,A313,
     1      D11,D13,D31,D33,
     2      D111,D131,D311,D331,D113,D133,D313,D333,
     3      D1111,D1311,D3111,D3311,D1133,D1333,D3133,D3333,
     3      D1113,D1313,D3113,D3313
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
      A11=2.E0*RE12
      A13=0.E0
      A111=(V**2-SR**2*U13**2+2.E0*SR**2*U3*U1)/(V**2+2.E0*V*SR**2*
     1      U13+SR**4*U13**2)
      A133=(SR**2*U3**2*RE12**2)/((V**2+2.E0*V*SR**2*U13+SR**4*U13
     1      **2)*RE32**2)
      A113=(2.E0*SR**2*U13*U3*CR*RE12)/((V**2+2.E0*V*SR**2*U13+SR**4
     1      *U13**2)*RE32)
      A31=0.E0
      A33=2.E0*RE32
      A311=(SR**2*U1**2*RE32**2)/((V**2+2.E0*V*SR**2*U13+SR**4*U13
     1      **2)*RE12**2)
      A333=(V**2-SR**2*U13**2+2.E0*SR**2*U3*U1)/(V**2+2.E0*V*SR**2*
     1      U13+SR**4*U13**2)
      A313=(2.E0*SR**2*U13*U1*CR*RE32)/((V**2+2.E0*V*SR**2*U13+SR**4
     1      *U13**2)*RE12)
      D11=1.E0
      D13=0.E0
      D31=0.E0
      D33=1.E0
      D111=(A113*B31+2.E0*A111*B11-2.E0)/(2.E0*RE12)
      D113=(A113*B33+2.E0*A111*B13)/(2.E0*RE12)
      D1111=(A113*RE12*B311-A113*B31+2.E0*A111*RE12*B111-2.E0*A111*
     1      B11+2.E0)/(2.E0*RE12**2)
      D1133=(A113*B333+2.E0*A111*B133)/(2.E0*RE12)
      D1113=(A113*RE12*B313-A113*B33+2.E0*A111*RE12*B113-2.E0*A111*
     1      B13)/(2.E0*RE12**2)
      D311=(B31*A313+2.E0*B11*A311)/(2.E0*RE32)
      D313=(B33*A313+2.E0*B13*A311)/(2.E0*RE32)
      D3111=(B311*A313+2.E0*B111*A311)/(2.E0*RE32)
      D3133=(B333*A313*RE32-B33*A313+2.E0*B133*A311*RE32-2.E0*B13*
     1      A311)/(2.E0*RE32**2)
      D3113=(B313*A313*RE32-B31*A313+2.E0*B113*A311*RE32-2.E0*B11*
     1      A311)/(2.E0*RE32**2)
      D131=(A113*B11+2.E0*A133*B31)/(2.E0*RE12)
      D133=(A113*B13+2.E0*A133*B33)/(2.E0*RE12)
      D1311=(A113*RE12*B111-A113*B11+2.E0*A133*RE12*B311-2.E0*A133*
     1      B31)/(2.E0*RE12**2)
      D1333=(A113*B133+2.E0*A133*B333)/(2.E0*RE12)
      D1313=(A113*RE12*B113-A113*B13+2.E0*A133*RE12*B313-2.E0*A133*
     1      B33)/(2.E0*RE12**2)
      D331=(2.E0*B31*A333+B11*A313)/(2.E0*RE32)
      D333=(2.E0*B33*A333+B13*A313-2.E0)/(2.E0*RE32)
      D3311=(2.E0*B311*A333+B111*A313)/(2.E0*RE32)
      D3333=(2.E0*B333*A333*RE32-2.E0*B33*A333+B133*A313*RE32-B13*
     1      A313+2.E0)/(2.E0*RE32**2)
      D3313=(2.E0*B313*A333*RE32-2.E0*B31*A333+B113*A313*RE32-B11*
     1      A313)/(2.E0*RE32**2)
      PX11=CHX13*D13*B11+CHX31*D11*B31+CHX33*D13*B31+CHX11*D11*
     1      B11
      PX13=CHX13*D13*B13+CHX31*D11*B33+CHX33*D13*B33+CHX11*D11*
     1      B13
      PX111=CHX13*D131*B11+CHX13*D13*B111+CHX31*D111*B31+CHX31*
     1      D11*B311+CHX33*D131*B31+CHX33*D13*B311+CHX11*D111*B11+
     1      CHX11*D11*B111
      PX133=CHX13*D133*B13+CHX13*D13*B133+CHX31*D113*B33+CHX31*
     1      D11*B333+CHX33*D133*B33+CHX33*D13*B333+CHX11*D113*B13+
     1      CHX11*D11*B133
      PX113=CHX13*D133*B11+CHX13*D131*B13+CHX13*D13*B113+CHX31*
     1      D113*B31+CHX31*D111*B33+CHX31*D11*B313+CHX33*D133*B31+
     1      CHX33*D131*B33+CHX33*D13*B313+CHX11*D113*B11+CHX11*D111*
     1      B13+CHX11*D11*B113
      PX1111=CHX13*D1311*B11+CHX13*D131*B111+CHX13*D13*B1111+
     1      CHX31*D1111*B31+CHX31*D111*B311+CHX31*D11*B3111+CHX33*
     1      D1311*B31+CHX33*D131*B311+CHX33*D13*B3111+CHX11*D1111*B11
     1      +CHX11*D111*B111+CHX11*D11*B1111
      PX1333=CHX13*D1333*B13+CHX13*D133*B133+CHX13*D13*B1333+
     1      CHX31*D1133*B33+CHX31*D113*B333+CHX31*D11*B3333+CHX33*
     1      D1333*B33+CHX33*D133*B333+CHX33*D13*B3333+CHX11*D1133*B13
     1      +CHX11*D113*B133+CHX11*D11*B1333
      PX1113=CHX13*D1313*B11+CHX13*D1311*B13+CHX13*D133*B111+
     1      CHX13*D131*B113+CHX13*D13*B1113+CHX31*D1113*B31+CHX31*
     1      D1111*B33+CHX31*D113*B311+CHX31*D111*B313+CHX31*D11*B3113
     1      +CHX33*D1313*B31+CHX33*D1311*B33+CHX33*D133*B311+CHX33*
     1      D131*B313+CHX33*D13*B3113+CHX11*D1113*B11+CHX11*D1111*B13
     1      +CHX11*D113*B111+CHX11*D111*B113+CHX11*D11*B1113
      PX1133=CHX13*D1313*B13+CHX13*D1333*B11+CHX13*D133*B113+
     1      CHX13*D131*B133+CHX13*D13*B1133+CHX31*D1113*B33+CHX31*
     1      D1133*B31+CHX31*D113*B313+CHX31*D111*B333+CHX31*D11*B3133
     1      +CHX33*D1313*B33+CHX33*D1333*B31+CHX33*D133*B313+CHX33*
     1      D131*B333+CHX33*D13*B3133+CHX11*D1113*B13+CHX11*D1133*B11
     1      +CHX11*D113*B113+CHX11*D111*B133+CHX11*D11*B1133
      PX31=CHX13*D33*B11+CHX31*D31*B31+CHX33*D33*B31+CHX11*D31*
     1      B11
      PX33=CHX13*D33*B13+CHX31*D31*B33+CHX33*D33*B33+CHX11*D31*
     1      B13
      PX311=CHX13*D331*B11+CHX13*D33*B111+CHX31*D311*B31+CHX31*
     1      D31*B311+CHX33*D331*B31+CHX33*D33*B311+CHX11*D311*B11+
     1      CHX11*D31*B111
      PX333=CHX13*D333*B13+CHX13*D33*B133+CHX31*D313*B33+CHX31*
     1      D31*B333+CHX33*D333*B33+CHX33*D33*B333+CHX11*D313*B13+
     1      CHX11*D31*B133
      PX313=CHX13*D333*B11+CHX13*D331*B13+CHX13*D33*B113+CHX31*
     1      D313*B31+CHX31*D311*B33+CHX31*D31*B313+CHX33*D333*B31+
     1      CHX33*D331*B33+CHX33*D33*B313+CHX11*D313*B11+CHX11*D311*
     1      B13+CHX11*D31*B113
      PX3111=CHX13*D3311*B11+CHX13*D331*B111+CHX13*D33*B1111+
     1      CHX31*D3111*B31+CHX31*D311*B311+CHX31*D31*B3111+CHX33*
     1      D3311*B31+CHX33*D331*B311+CHX33*D33*B3111+CHX11*D3111*B11
     1      +CHX11*D311*B111+CHX11*D31*B1111
      PX3333=CHX13*D3333*B13+CHX13*D333*B133+CHX13*D33*B1333+
     1      CHX31*D3133*B33+CHX31*D313*B333+CHX31*D31*B3333+CHX33*
     1      D3333*B33+CHX33*D333*B333+CHX33*D33*B3333+CHX11*D3133*B13
     1      +CHX11*D313*B133+CHX11*D31*B1333
      PX3113=CHX13*D3313*B11+CHX13*D3311*B13+CHX13*D333*B111+
     1      CHX13*D331*B113+CHX13*D33*B1113+CHX31*D3113*B31+CHX31*
     1      D3111*B33+CHX31*D313*B311+CHX31*D311*B313+CHX31*D31*B3113
     1      +CHX33*D3313*B31+CHX33*D3311*B33+CHX33*D333*B311+CHX33*
     1      D331*B313+CHX33*D33*B3113+CHX11*D3113*B11+CHX11*D3111*B13
     1      +CHX11*D313*B111+CHX11*D311*B113+CHX11*D31*B1113
      PX3133=CHX13*D3313*B13+CHX13*D3333*B11+CHX13*D333*B113+
     1      CHX13*D331*B133+CHX13*D33*B1133+CHX31*D3113*B33+CHX31*
     1      D3133*B31+CHX31*D313*B313+CHX31*D311*B333+CHX31*D31*B3133
     1      +CHX33*D3313*B33+CHX33*D3333*B31+CHX33*D333*B313+CHX33*
     1      D331*B333+CHX33*D33*B3133+CHX11*D3113*B13+CHX11*D3133*B11
     1      +CHX11*D313*B113+CHX11*D311*B133+CHX11*D31*B1133
      PR11=CHR13*D13*B11+CHR31*D11*B31+CHR33*D13*B31+CHR11*D11*
     1      B11
      PR13=CHR13*D13*B13+CHR31*D11*B33+CHR33*D13*B33+CHR11*D11*
     1      B13
      PR111=CHR13*D131*B11+CHR13*D13*B111+CHR31*D111*B31+CHR31*
     1      D11*B311+CHR33*D131*B31+CHR33*D13*B311+CHR11*D111*B11+
     1      CHR11*D11*B111
      PR133=CHR13*D133*B13+CHR13*D13*B133+CHR31*D113*B33+CHR31*
     1      D11*B333+CHR33*D133*B33+CHR33*D13*B333+CHR11*D113*B13+
     1      CHR11*D11*B133
      PR113=CHR13*D133*B11+CHR13*D131*B13+CHR13*D13*B113+CHR31*
     1      D113*B31+CHR31*D111*B33+CHR31*D11*B313+CHR33*D133*B31+
     1      CHR33*D131*B33+CHR33*D13*B313+CHR11*D113*B11+CHR11*D111*
     1      B13+CHR11*D11*B113
      PR1111=CHR13*D1311*B11+CHR13*D131*B111+CHR13*D13*B1111+
     1      CHR31*D1111*B31+CHR31*D111*B311+CHR31*D11*B3111+CHR33*
     1      D1311*B31+CHR33*D131*B311+CHR33*D13*B3111+CHR11*D1111*B11
     1      +CHR11*D111*B111+CHR11*D11*B1111
      PR1333=CHR13*D1333*B13+CHR13*D133*B133+CHR13*D13*B1333+
     1      CHR31*D1133*B33+CHR31*D113*B333+CHR31*D11*B3333+CHR33*
     1      D1333*B33+CHR33*D133*B333+CHR33*D13*B3333+CHR11*D1133*B13
     1      +CHR11*D113*B133+CHR11*D11*B1333
      PR1113=CHR13*D1313*B11+CHR13*D1311*B13+CHR13*D133*B111+
     1      CHR13*D131*B113+CHR13*D13*B1113+CHR31*D1113*B31+CHR31*
     1      D1111*B33+CHR31*D113*B311+CHR31*D111*B313+CHR31*D11*B3113
     1      +CHR33*D1313*B31+CHR33*D1311*B33+CHR33*D133*B311+CHR33*
     1      D131*B313+CHR33*D13*B3113+CHR11*D1113*B11+CHR11*D1111*B13
     1      +CHR11*D113*B111+CHR11*D111*B113+CHR11*D11*B1113
      PR1133=CHR13*D1313*B13+CHR13*D1333*B11+CHR13*D133*B113+
     1      CHR13*D131*B133+CHR13*D13*B1133+CHR31*D1113*B33+CHR31*
     1      D1133*B31+CHR31*D113*B313+CHR31*D111*B333+CHR31*D11*B3133
     1      +CHR33*D1313*B33+CHR33*D1333*B31+CHR33*D133*B313+CHR33*
     1      D131*B333+CHR33*D13*B3133+CHR11*D1113*B13+CHR11*D1133*B11
     1      +CHR11*D113*B113+CHR11*D111*B133+CHR11*D11*B1133
      PR31=CHR13*D33*B11+CHR31*D31*B31+CHR33*D33*B31+CHR11*D31*
     1      B11
      PR33=CHR13*D33*B13+CHR31*D31*B33+CHR33*D33*B33+CHR11*D31*
     1      B13
      PR311=CHR13*D331*B11+CHR13*D33*B111+CHR31*D311*B31+CHR31*
     1      D31*B311+CHR33*D331*B31+CHR33*D33*B311+CHR11*D311*B11+
     1      CHR11*D31*B111
      PR333=CHR13*D333*B13+CHR13*D33*B133+CHR31*D313*B33+CHR31*
     1      D31*B333+CHR33*D333*B33+CHR33*D33*B333+CHR11*D313*B13+
     1      CHR11*D31*B133
      PR313=CHR13*D333*B11+CHR13*D331*B13+CHR13*D33*B113+CHR31*
     1      D313*B31+CHR31*D311*B33+CHR31*D31*B313+CHR33*D333*B31+
     1      CHR33*D331*B33+CHR33*D33*B313+CHR11*D313*B11+CHR11*D311*
     1      B13+CHR11*D31*B113
      PR3111=CHR13*D3311*B11+CHR13*D331*B111+CHR13*D33*B1111+
     1      CHR31*D3111*B31+CHR31*D311*B311+CHR31*D31*B3111+CHR33*
     1      D3311*B31+CHR33*D331*B311+CHR33*D33*B3111+CHR11*D3111*B11
     1      +CHR11*D311*B111+CHR11*D31*B1111
      PR3333=CHR13*D3333*B13+CHR13*D333*B133+CHR13*D33*B1333+
     1      CHR31*D3133*B33+CHR31*D313*B333+CHR31*D31*B3333+CHR33*
     1      D3333*B33+CHR33*D333*B333+CHR33*D33*B3333+CHR11*D3133*B13
     1      +CHR11*D313*B133+CHR11*D31*B1333
      PR3113=CHR13*D3313*B11+CHR13*D3311*B13+CHR13*D333*B111+
     1      CHR13*D331*B113+CHR13*D33*B1113+CHR31*D3113*B31+CHR31*
     1      D3111*B33+CHR31*D313*B311+CHR31*D311*B313+CHR31*D31*B3113
     1      +CHR33*D3313*B31+CHR33*D3311*B33+CHR33*D333*B311+CHR33*
     1      D331*B313+CHR33*D33*B3113+CHR11*D3113*B11+CHR11*D3111*B13
     1      +CHR11*D313*B111+CHR11*D311*B113+CHR11*D31*B1113
      PR3133=CHR13*D3313*B13+CHR13*D3333*B11+CHR13*D333*B113+
     1      CHR13*D331*B133+CHR13*D33*B1133+CHR31*D3113*B33+CHR31*
     1      D3133*B31+CHR31*D313*B313+CHR31*D311*B333+CHR31*D31*B3133
     1      +CHR33*D3313*B33+CHR33*D3333*B31+CHR33*D333*B313+CHR33*
     1      D331*B333+CHR33*D33*B3133+CHR11*D3113*B13+CHR11*D3133*B11
     1      +CHR11*D313*B113+CHR11*D311*B133+CHR11*D31*B1133
      PX1(1)=0.0D0
      PX1(2)=PX11
      PX1(3)=PX13
      PX1(4)=PX111
      PX1(5)=PX133
      PX1(6)=PX113
      PX1(7)=PX1111
      PX1(8)=PX1333
      PX1(9)=PX1113
      PX1(10)=PX1133
      PX1(11)=0.0D+00
      PX1(12)=0.0D+00
      PX1(13)=0.0D+00
      PX1(14)=0.0D+00
      PX1(15)=0.0D+00
      PX3(1)=0.0D0
      PX3(2)=PX31
      PX3(3)=PX33
      PX3(4)=PX311
      PX3(5)=PX333
      PX3(6)=PX313
      PX3(7)=PX3111
      PX3(8)=PX3333
      PX3(9)=PX3113
      PX3(10)=PX3133
      PX3(11)=0.0D+00
      PX3(12)=0.0D+00
      PX3(13)=0.0D+00
      PX3(14)=0.0D+00
      PX3(15)=0.0D+00
      PR1(1)=0.0D0
      PR1(2)=PR11
      PR1(3)=PR13
      PR1(4)=PR111
      PR1(5)=PR133
      PR1(6)=PR113
      PR1(7)=PR1111
      PR1(8)=PR1333
      PR1(9)=PR1113
      PR1(10)=PR1133
      PR1(11)=0.0D+00
      PR1(12)=0.0D+00
      PR1(13)=0.0D+00
      PR1(14)=0.0D+00
      PR1(15)=0.0D+00
      PR3(1)=0.0D0
      PR3(2)=PR31
      PR3(3)=PR33
      PR3(4)=PR311
      PR3(5)=PR333
      PR3(6)=PR313
      PR3(7)=PR3111
      PR3(8)=PR3333
      PR3(9)=PR3113
      PR3(10)=PR3133
      PR3(11)=0.0D+00
      PR3(12)=0.0D+00
      PR3(13)=0.0D+00
      PR3(14)=0.0D+00
      PR3(15)=0.0D+00
      RETURN
      END
