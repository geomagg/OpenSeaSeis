c-----------------------------------------------------------------------
c     Modelagem 2D acustica no dominio do tempo do programa XDVSSYN2DRS
c     (SEISPAK, "syn3d"), para o modulo SPK_XDVSSYN2D.
c
c     SPKRS  : programa principal (fonte) + syn3 (laco no tempo), com:
c              - velocidade e densidade vindas da memoria (antes:
c                arquivos ONLVELS / ONLRHOS lidos por velcom4/rhocom4)
c              - registro na profundidade IDP e snapshots devolvidos em
c                OUT (antes: arquivos '3dtime', 'snapshot', ONLINE1)
c              - wavelet opcional vinda da memoria (antes: arquivo WAVELET)
c              Diferencas marcadas com SEASEIS:.
c
c     Operadores (APX), como no original:
c       sem densidade: 0,2,3 xrhsfft (pseudo-espectral)
c                      1     rhsfft  (pseudo-espectral 2D)
c                      4..7  xrhscon (x espectral, z diferencas finitas
c                            de meia-largura APX, coeficientes de Holberg)
c       com densidade: 0,2,4..7 test2dr, 1 rhsfftr,
c                      3 xtest3drs (grade intercalada, meia amostra)
c     KUTTA termos da expansao de cos(dt*sqrt(-L)) por passo de tempo.
c-----------------------------------------------------------------------
      SUBROUTINE SPKRS(VEL,RHOIN,IDEN,NX,NZ,NT,DT,DX,DZ,APX,KUTTA,
     :                 RAA,NX0,NZ0,IDP,GAM1,FS,ITYPE,ITZR,WAVE,NTW,
     :                 OUT,NOUT,IERR)
c     VEL(NZ,NX)    velocidade intervalar
c     RHOIN(NZ,NX)  densidade (usada se IDEN=1)
c     ITYPE         2 = secao registrada em IDP: OUT(NT,NX)
c                   1 = snapshots a cada ITZR passos: OUT(NZ,NX,NSNAP)
c     WAVE(NTW)     wavelet injetada na fonte (NTW=0: sem wavelet)
c     IERR          0 ok, 2 sem memoria
      REAL VEL(*),RHOIN(*),WAVE(*),OUT(*)
      REAL, ALLOCATABLE :: A(:)
      COMMON/SPKRSIT/ITABS,ITDZZ
      IERR=0
      NY=1
      NXZ=NX*NZ
      NXYZ=NXZ
c     SEASEIS: u(2*nxyz) | v | rho | s(7*nxyz+folga)
      NA=2*NXYZ+NXYZ+NXYZ+7*NXYZ+100000
      ALLOCATE(A(NA),STAT=IST)
      IF(IST.NE.0)THEN
        IERR=2
        RETURN
      ENDIF
      DO J=1,NA
        A(J)=0.
      ENDDO
c     SEASEIS: reinicia os contadores 'itime' de absorb1 e ddzcon2
      ITABS=0
      ITDZZ=0
      I1=1
      I2=I1+2*NXYZ
      I3=I2+NXYZ
      I4=I3+NXYZ
c   ***   test spike (programa principal syn3d)
      SCALE=+150000.
      NY0=1
      II=(NY0-1)*NXZ+(NX0-1)*NZ+NZ0
      A(II)=SCALE
      IF(RAA.EQ.0.)CALL XYZFILT(A,A(I2),NX0,NY0,NZ0,NX,NY,NZ)
      IF(RAA .EQ. 1.)THEN
        A(II-2-NZ)=.25*SCALE
        A(II-2)   = .5*SCALE
        A(II-2+NZ)=.25*SCALE
        A(II-1-NZ)= .5*SCALE
        A(II-1)   = 1.*SCALE
        A(II-1+NZ)= .5*SCALE
        A(II)     =0.
        A(II+NZ)  =0.
        A(II-NZ)  =0.
        A(II+1-NZ)= -.5*SCALE
        A(II+1)   =   -SCALE
        A(II+1+NZ)= -.5*SCALE
        A(II+2-NZ)=-.25*SCALE
        A(II+2)   = -.5*SCALE
        A(II+2+NZ)=-.25*SCALE
      ENDIF
      CALL SPKRSW(A(I1),A(I2),A(I3),A(I4),VEL,RHOIN,IDEN,NX,NZ,NT,
     :            DT,DX,DZ,APX,KUTTA,RAA,NX0,NZ0,IDP,GAM1,FS,ITYPE,
     :            ITZR,WAVE,NTW,OUT,NOUT)
      DEALLOCATE(A)
      RETURN
      END
c
      SUBROUTINE SPKRSW(U,V,RHO,S,VEL,RHOIN,IDEN,NX,NZ,NT,DT,DX,DZ,APX,
     :                  KUTTA,RAA,NX0,NZ0,IDP,GAM1,FS,ITYPE,ITZR,WAVE,
     :                  NTW,OUT,NOUT)
c     syn3 de xdvssyn2drs.f
      REAL U(*),V(*),RHO(*),S(*),VEL(*),RHOIN(*),WAVE(*),OUT(*)
      LOGICAL EXIST
      NY=1
      NXY=NX
      NXZ=NX*NZ
      NXYZ=NXZ
      NTT=NT
      DDZ=DZ
      DY=DX
      I1=1
      I2=I1+NXYZ
      IS2=NXYZ+1
      IS3=IS2+NXYZ
c   ***   bring in velocities
c     SEASEIS: velocidade e densidade da memoria (antes velcom4/rhocom4)
      DO J=1,NXYZ
        V(J)=VEL(J)
      ENDDO
      IF(IDEN.EQ.1)THEN
        DO J=1,NXYZ
          RHO(J)=RHOIN(J)
        ENDDO
      ELSE
        DO J=1,NXYZ
          RHO(J)=1.
        ENDDO
      ENDIF
c   set surface velocity into bottom absorbing boundary region
      DO K=1,NX
      DO J=NZ-2*ABS(GAM1),NZ
      V((K-1)*NZ+J)=V((K-1)*NZ+1)
      RHO((K-1)*NZ+J)=RHO((K-1)*NZ+1)
      ENDDO
      ENDDO
c     SEASEIS: ITYPE 5 devolve a velocidade usada (antes: arquivo veltest
c              com rho*(v*dt)**2)
      IF(ITYPE.GE.5)THEN
        DO J=1,NXYZ
          OUT(J)=V(J)
        ENDDO
        RETURN
      ENDIF
c   rho*(v*dt)**2 -> v   :   1/rho -> rho
      DO J=1,NXYZ
      V(J)=RHO(J)*(V(J)*DT)**2
      RHO(J)=1./RHO(J)
      ENDDO
c   ***   bring sources in
      IF(RAA .EQ. 2)THEN
        DO J=1,NXYZ
          U(J)=0.
        ENDDO
        U((NX0-1)*NZ+NZ0)=1000.
        CALL FFTFLD2(U,S,V,DZ,DX,NZ,NX)
      ENDIF
      CALL MOVLEV(U(I1),U(I2),NXYZ)
      HFL=6
      S(4)=HFL
c     SEASEIS: 'scale' era variavel local sem valor quando RAA=0 (so e
c              usada em RAA 3/4)
      SCALE=0.
      IF(RAA.GT.0.)SCALE=10000
      IF(RAA.LT.0.)SCALE=10000
      RAA1=ABS(RAA)
      ISWS=1
      IF(RAA1.EQ.3.)CALL SPGAUSSR(NX0,NZ0,4.*SCALE,U,DZ,DX,NZ,NX,ISWS,S)
      IF(RAA1.EQ.4.)CALL SPGAUSSU(NX0,NZ0,4.*SCALE,U,DZ,DX,NZ,NX,ISWS,S)
      CALL CLEAR(U(I1),NXYZ)
c  ***   wavelet (antes: arquivo WAVELET)
      EXIST=(NTW.GT.1)
      IIU=1
      NSNAP=0
      DO 10 I=1,NTT
c  ***   read wavelet into field if exist
      IF(EXIST.AND. I.LT.NTW)THEN
        IF(I.EQ.1)THEN
          U(I1+(NX0-1)*NZ+NZ0-1)=1000.
          CALL MOVLEV(U(I1),U(I2),NXYZ)
          DO J=1,NXYZ
            U(I1+J-1)=WAVE(1)*U(I1+J-1)
            U(I2+J-1)=WAVE(2)*U(I2+J-1)
          ENDDO
        ELSE
          U(I2+(NX0-1)*NZ+NZ0-1)=U(I2+(NX0-1)*NZ+NZ0-1)+1000.*WAVE(I+1)
        ENDIF
      ENDIF
      DO 410 J=1,NXYZ
      U(I1+J-1)=2.*U(I2+J-1)-U(I1+J-1)
410   CONTINUE
      FACT=2.
      CALL MOVLEV(U(I2),S,NXYZ)
      DO K=1,KUTTA
      IF(APX.EQ.0.)THEN
      IF(IDEN.EQ.0)CALL XRHSFFT(S,S,S(IS2),S(IS3),V,DDZ,DX,NZ,NX)
      IF(IDEN.EQ.1)CALL TEST2DR(S,V,RHO,S(IS2),NZ,NX,DDZ,DX)
      ENDIF
      IF(APX.EQ.1.)THEN
      IF(IDEN.EQ.0)CALL RHSFFT(S,S,S(IS3),V,DDZ,DX,NZ,NX)
      IF(IDEN.EQ.1)CALL RHSFFTR(S,S,S(IS3),RHO,V,DDZ,DX,NZ,NX)
      ENDIF
      IF(APX.EQ.2.)THEN
      IF(IDEN.EQ.0)CALL XRHSFFT(S,S,S(IS2),S(IS3),V,DDZ,DX,NZ,NX)
      IF(IDEN.EQ.1)CALL TEST2DR(S,V,RHO,S(IS2),NZ,NX,DDZ,DX)
      ENDIF
      IF(APX.EQ.3.)THEN
      IF(IDEN.EQ.0)CALL XRHSFFT(S,S,S(IS2),S(IS3),V,DDZ,DX,NZ,NX)
      IF(IDEN.EQ.1)CALL XTEST3DRS(S,S(IS2),V,RHO,NZ,NX,NY,DDZ,DX,DY)
      ENDIF
      IF(APX.GE.4.)THEN
      IAPX=NINT(APX)
      IF(IDEN.EQ.0)CALL XRHSCON(S,S,S(IS2),S(IS3),V,DDZ,DX,NZ,NX,IAPX)
      IF(IDEN.EQ.1)CALL TEST2DR(S,V,RHO,S(IS2),NZ,NX,DDZ,DX)
      ENDIF
      FACT=FACT/(2*K*(2*K-1))
      DO M=1,NX*NY
      K1=(M-1)*NZ+I1-1
      K2=(M-1)*NZ-1
      DO J=1,NZ
      U(K1+J)=U(K1+J)+FACT*S(K2+J+1)
      ENDDO
      ENDDO
      ENDDO
c   ***   apply absorbing boundary
      TAP1=.002
      TAP2=.002
      CALL ABSORB1(U(I1),U(I2),V,NZ,NX,NY,DZ,DX,DY,GAM1,FS,TAP1,TAP2)
c   ***   switch indicies
      II=I1
      I1=I2
      I2=II
c   ***   put out surface recording at this time level
c     SEASEIS: direto para OUT(NT,NX) (antes: pi2 '3dtime' + trni2da)
      IF(ITYPE.EQ.2)THEN
        DO J=1,NX
          JJ=(J-1)*NZ+IDP
          OUT((J-1)*NT+I)=U(I1+JJ-1)
        ENDDO
      ENDIF
c     SEASEIS: MOD(I-1,ITZR).EQ.0 (o original, MOD(I,ITZR).EQ.1, nao
c              gravava nada com ITZR=1); snapshots direto para OUT
      IF(ITYPE.EQ.1 .AND. MOD(I-1,ITZR).EQ.0 .AND. IIU+NXZ-1.LE.NOUT)
     :  THEN
        CALL MOVLEV(U(I1),OUT(IIU),NXZ)
        IIU=IIU+NXZ
      ENDIF
10    CONTINUE
      RETURN
      END
