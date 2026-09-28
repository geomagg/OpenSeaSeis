c-----------------------------------------------------------------------
c     Deconvolucao preditiva multi-janela do processo JHHDCON (SEISPAK),
c     adaptada para uso traco a traco no SeaSeis.
c
c     SPKDCON = corpo de JHHDCON (a partir de "NUMBER CRUNCHING"), com:
c       - parametros recebidos como argumentos (antes vinham da TABLE1)
c       - traco sem header (LHDR=0); o offset (RANG) chega em RNG
c       - VEL(I) passado a DCONSUB para cada janela (antes ia VEL(1))
c     DCONSUB = rotina original, com as correcoes marcadas "SEASEIS:"
c-----------------------------------------------------------------------
      SUBROUTINE SPKDCON(S,NT,SMP,NSTM,STM,ETM,VEL,FLNGTH,GAP,WND
     :                  ,MODE,FR,RNG,IRANG,R,Y,X,F,AIGAP,AIGSQ,NCNT
     :                  ,IPFLG)
      REAL S(*),R(*),Y(*),X(NT,*),F(*)
      REAL WND(2,*),GAP(*),FLNGTH(*),STM(*),ETM(*),VEL(*)
      REAL AIGAP(*),AIGSQ(*)
      INTEGER NCNT(*)
      LHDR=0
      DO 100 I=1,NSTM
        ISTRT0=1
        LNGF=nint(FLNGTH(I)/SMP)
        if(MODE.EQ.2)then
           JGAP=nint(GAP(I)/SMP)
        else
           JGAP=nint(GAP(I))
        endif
        TMW=0.001*WND(1,I)
        IW10=nint(WND(1,I)/SMP)
        IW20=nint(WND(2,I)/SMP)
        IF(IW20.GT.NT)IW20=NT
        CALL DCONSUB(LHDR,NT,SMP,S,R,Y,X(1,I),F,WND(1,I),MODE,ISTRT0
     :              ,VEL(I),FR,LNGF,JGAP,TMW,IW10,IW20,RNG
     :              ,AIGAP(I),AIGSQ(I),NCNT(I),IPFLG,IRANG)
  100 CONTINUE
      IPFLG=0
      IW10=nint(WND(1,1)/SMP)
      TMW=0.001*WND(1,1)
      if(vel(1).gt.0.0)then
         xov=rng/vel(1)
         jbias=nint(1000.0*sqrt(tmw*tmw+xov*xov)/smp)-IW10
      else
         jbias=0
      endif
      KTME=0
      DO 200 I=1,NSTM
         ibias=jbias
         IF(I+1.LE.NSTM)THEN
            IW10=nint(WND(1,I+1)/SMP)
            TMW=0.001*WND(1,I+1)
            if(vel(I).gt.0.0)then
               xov=rng/vel(I)
               jbias=nint(1000.0*sqrt(tmw*tmw+xov*xov)/smp)-IW10
            else
               jbias=0
            endif
         ELSE
            jbias=ibias
         ENDIF
         ITMS=1+NINT(STM(I)/SMP)+ibias
         ITME=1+NINT(ETM(I)/SMP)+jbias
         IF(ITME.GT.NT)ITME=NT
         IF(ITME.LT. 1)ITME=1
         IF(ITMS.GT.NT)ITMS=NT
         IF(ITMS.LT. 1)ITMS=1
         IF(I.EQ.1)THEN
            DO 140 J=ITMS,ITME
               S(LHDR+J)=X(J,I)
  140       CONTINUE
         ELSE
            IF(KTME.GT.ITMS)THEN
               WT=1.0/(KTME-ITMS)
               Z=0.0
               DO 160 J=ITMS,KTME
                  S(LHDR+J)=S(LHDR+J)*(1.0-Z) + X(J,I)*Z
                  Z=Z+WT
  160          CONTINUE
               KTME2USE=KTME
            ELSE
               KTME2USE=ITMS
            ENDIF
            DO 180 J=KTME2USE+1,ITME
               S(LHDR+J)=X(J,I)
  180       CONTINUE
         ENDIF
         KTME=ITME
  200 CONTINUE
      RETURN
      END
c
      SUBROUTINE DCONSUB(LHDR,NT,SMP,S,R,Y,X,F,WND,MODE,ISTRT0,VEL,FR
     :                  ,LNGF,JGAP,TMW,IW10,IW20,rng,aigap,aigsq,ncnt
     :                  ,IPFLG,IRANG)
      DIMENSION WND(2),S(*),R(*),Y(*),X(*),F(*)
      NWD=LHDR+NT
      ISTRT=ISTRT0
      if(IRANG.GT.0 .AND. VEL.GT.0.0)THEN
         xov=rng/vel
         ibias=nint(1000.0*sqrt(tmw*tmw+xov*xov)/smp)-IW10
         iw1=iw10+ibias
         iw2=iw20+ibias
         if(iw2.gt.nt)then
            iw2=nt
            iw1=nt-(iw20-iw10)
         endif
         ISTRT=ISTRT+ibias
         if(istrt.gt.nt)istrt=nt
      else
         iw1=iw10
         iw2=iw20
      endif
      IF(IW1.LT.1)IW1=1
      IF(IW2.GT.NT)IW2=NT
      NDX=LHDR+IW1
      LNGW=IW2-IW1+1
c     SEASEIS: modos 0, 1 e 3 leem R() alem de LNGF; calcular ate LNGW
      NLAG=LNGF
      IF(MODE.NE.2)NLAG=MAX(LNGF,LNGW)
c						compute autocorrelation
      DO 20 I=1,NLAG
         CALL XCORL(R(I),S(NDX),S(NDX),LNGW,I-1)
         IF(I.EQ.1)THEN
            IF(R(1).GT.0.0)GO TO 20
c						autocorrelation zero
c						=> trace is dead in window
            GO TO 401
         ENDIF
   20 CONTINUE
      if(mode.eq.0)then
c     SEASEIS: saida em X (antes em S, e X ficava sem valor)
         DO 30 I=1,MIN(LNGW,NT)
            X(I)=R(I)
   30    CONTINUE
         DO 35 I=LNGW+1,NT
            X(I)=0.0
   35    CONTINUE
         return
      endif
c     SEASEIS: valor inicial de igap (modo 3 podia nao definir)
      igap=1
      if(mode.eq.3)then
c						mode=3 causes igap to be
c						set to the location of max
c						autocorrelation between
c						zero crossing -2*jgap and
c						1-2*jgap
         kgap=jgap+jgap+1
         rmax=0
         do 42 i=2,lngw/2
            if(kgap.eq.1 .and. R(i).gt.rmax)then
              rmax=R(i)
              igap=i
            endif
            if(R(i-1)*R(i).ge.0.0) go to 42
            kgap=kgap-1
            if(kgap.eq.0)go to 43
   42    continue
   43    continue
c						mode=1 causes igap to be
c						set to the location of
c						the jgapth zero crossing
      elseif(mode.eq.1)then
         kgap=jgap
         do 52 i=2,lngw/2
            igap=i
            if(R(i-1)*R(i).le.0.0)kgap=kgap-1
            if(kgap.eq.0)go to 53
   52    continue
   53    continue
      else
        igap=jgap
      endif
      ncnt=ncnt+1
      aigap=aigap+igap
      aigsq=aigsq+igap*igap
c						add small stabilizing
c						value to diagonal (noise)
c     SEASEIS: filtro sem coeficientes (gap >= comprimento): nao altera
      IF(R(1).eq.0.0 .OR. LNGF-IGAP.LT.1)THEN
         DO 70 I=1,NT
            X(I)=S(LHDR+I)
   70    CONTINUE
         RETURN
      ENDIF
      R(1)=FR*R(1)
c						set up "right hand side"
c						to predict wavelet at a
c						distance of "GAP" beyond
c						front end of the wavelet
      DO 250 I=1,LNGF-IGAP
         Y(I)=R(I+IGAP)
  250 CONTINUE
c						solve toeplitz matrix
c						using levinson recursion
      CALL TPLTZR(R,Y,F,LNGF-IGAP)
c						convolve prediction filter
c						with data
      IF(IPFLG.EQ.1)THEN
         WRITE(*,*)
         WRITE(*,*)'Decon filter for first trace:'
         WRITE(*,2001)(F(I),I=1,LNGF-IGAP)
         WRITE(*,*)
 2001    FORMAT(10(1X,F7.4))
      ENDIF
      CALL CONV(F,LNGF-IGAP,S(LHDR+1),NT,R)
c						subtract predicted part of
c						energy from trace
      I1=1+IGAP
      I2=ISTRT
      I0=MAX(I1,I2)
      DO 290 I=1,I0-1
         X(I)=S(LHDR+I)
  290 CONTINUE
      DO 300 I=I0,NT
         X(I)=S(LHDR+I)-R(I-IGAP)
  300 CONTINUE
      GO TO 999
  401 CONTINUE
      DO 410 I=1,NT
         X(I)=0.0
  410 CONTINUE
  999 CONTINUE
      RETURN
      END
