c-----------------------------------------------------------------------
c     XCORL, TPLTZR, CONV: rotinas auxiliares chamadas por DCONSUB.
c     As originais da biblioteca SEISPAK nao estavam disponiveis; estas
c     foram escritas a partir da forma como DCONSUB as usa.
c     Se as originais forem encontradas, podem substituir este arquivo.
c-----------------------------------------------------------------------
c
c     XCORL: correlacao de A com B no atraso LAG, N amostras
c            R = soma(i=1..N-LAG) A(i+LAG)*B(i)
c
      SUBROUTINE XCORL(R,A,B,N,LAG)
      REAL A(*),B(*)
      R=0.0
      DO 10 I=1,N-LAG
         R=R+A(I+LAG)*B(I)
   10 CONTINUE
      RETURN
      END
c
c     TPLTZR: resolve o sistema de Toeplitz simetrico T*F = Y,
c             T(i,j) = R(|i-j|+1), i,j=1..N, pela recursao de Levinson
c
      SUBROUTINE TPLTZR(R,Y,F,N)
      REAL R(*),Y(*),F(*)
      PARAMETER (MAXN=4000)
      REAL A(MAXN),B(MAXN)
      IF(N.LT.1)RETURN
      IF(N.GT.MAXN)THEN
         WRITE(*,*)'TPLTZR: filter too long, N =',N,' > ',MAXN
         STOP
      ENDIF
      F(1)=Y(1)/R(1)
      IF(N.EQ.1)RETURN
      A(1)=1.0
      V=R(1)
      DO 100 M=2,N
c                                   erro de predicao (filtro de Levinson)
         E=0.0
         DO 10 J=1,M-1
            E=E+A(J)*R(M-J+1)
   10    CONTINUE
         C=-E/V
         B(1)=1.0
         DO 20 J=2,M-1
            B(J)=A(J)+C*A(M-J+1)
   20    CONTINUE
         B(M)=C
         DO 30 J=1,M
            A(J)=B(J)
   30    CONTINUE
         V=V*(1.0-C*C)
c                                   atualiza a solucao F
         G=0.0
         DO 40 J=1,M-1
            G=G+F(J)*R(M-J+1)
   40    CONTINUE
         Q=(Y(M)-G)/V
         DO 50 J=1,M-1
            F(J)=F(J)+Q*A(M-J+1)
   50    CONTINUE
         F(M)=Q
  100 CONTINUE
      RETURN
      END
c
c     CONV: convolucao de F (NF) com S, primeiras NT amostras
c           R(k) = soma(j=1..min(k,NF)) F(j)*S(k-j+1)
c
      SUBROUTINE CONV(F,NF,S,NT,R)
      REAL F(*),S(*),R(*)
      DO 20 K=1,NT
         SUM=0.0
         DO 10 J=1,MIN(K,NF)
            SUM=SUM+F(J)*S(K-J+1)
   10    CONTINUE
         R(K)=SUM
   20 CONTINUE
      RETURN
      END
