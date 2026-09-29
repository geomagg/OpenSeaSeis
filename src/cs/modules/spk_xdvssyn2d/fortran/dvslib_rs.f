c-----------------------------------------------------------------------
c     Rotinas copiadas sem alteracao de dvslib.f (SEISPAK).
c     Usadas pelo modulo SPK_XDVSSYN2D.
c-----------------------------------------------------------------------
      subroutine clear(a,n) 
      dimension a(*)
      do 10 i=1,n 
10     a(i)=0.
      return
      end 
      subroutine movlev(a,b,n)
      dimension a(*),b(*) 
      do 1 i=1,n
1     b(i)=a(i) 
      return
      end 
      subroutine rtran(a,b,nr,nc) 
      real a(*),b(*) 
      do 10 k=1,nc
      do 10 i=1,nr
      b((i-1)*nc+k)=a((k-1)*nr+i)
10    continue
      return
      end 
      subroutine fftl2d(y,a,nx,ny,isgn)
      complex y(*),a(*)
      if(nx .gt. 1)call fftc(y,a,2,2*nx,nx,ny,-isgn)
      if(ny .gt. 1)call fftc(y,a,2*nx,2,ny,nx,-isgn)
      scale=1./(real(nx)*real(ny))
      if(isgn.eq.-1)then
         do 1 j=1,nx*ny
             y(j)=scale*y(j)
    1    continue
      end if
      return
      end
      SUBROUTINE FFTC(DATA,WORK,INC,JUMP,LEN,LOT,IDIR)
C
C ABSTRACT *************************************************************
C
C     COMPLEX TO COMPLEX, PERIODIC, MULTIPLE FFTS, POWERS OF 2,3,5
C     DATA HAS REAL AND IMAGINARY PARTS IN SINGLE ARRAY, WITH EITHER
C     REAL/IMAG PARTS IN EVERY OTHER ROW, OR EVERY OTHER COLUMN
C
C KEYWORDS
C  FFT
C
C PURPOSE
C
C     COMPLEX TO COMPLEX, PERIODIC, MULTIPLE FFTS, POWERS OF 2,3,5
C     DATA HAS REAL AND IMAGINARY PARTS IN SINGLE ARRAY, WITH EITHER
C     REAL/IMAG PARTS IN EVERY OTHER ROW, OR EVERY OTHER COLUMN
C
C     DATA       REAL ARRAY      IN/OUT      DATA ARRAY OF LENGTH
C                                            2*LEN*LOT (REAL)
C     WORK       REAL ARRAY      SCRATCH     WORK ARRAY OF LENGTH
C         LENGTH SHOULD BE FOUND WITH ROUTINE FFTCL
C     INC        INTEGER         IN    INCREMENT BETWEEN SUCCESSIVE
C                                      REAL ELEMENTS OF DATA
C     JUMP       INTEGER         IN    INCREMENT BETWEEN FIRST ELEMENT
C                                      OF SUCCESSIVE VECTORS
C     LEN        INTEGER         IN    LENGTH OF EACH FFT ( PRIME
C                                      FACTORS MUST BE 2, 3, 5 )
C     LOT        INTEGER         IN    NUMBER OF TRANSFORMS
C     IDIR       INTEGER         IN    SIGN IN EXPONENT; NORMALIZATION
C                                      IS 1 INDEPENDENT OF IDIR
C
C  X(K) = SUM(N=0,LEN-1) X(N)EXP(2*PI*I*IDIR*K*N/LEN)
C
C
C
C  NORMAL USE IS EITHER
C     1)  FFT DOWN COLUMNS (ROW 1,3,... ARE REAL PARTS; 2,4,... ARE IMAG)
C         JUMP ROWS * LOT COLUMNS, INC = 2, JUMP = 2*LEN
C  OR
C     2)  FFT ACROSS ROWS (COL 1,3,... ARE REAL PARTS; 2,4,... ARE IMAG)
C         INC ROWS * LEN COLUMNS, JUMP = 2, INC = 2*LOT
C
C  EXAMPLES
C     1)  INPUT IS A REAL MATRIX WITH 2*LEN ROWS (ROW 1 REAL, ROW 2 IMAG)
C         AND LOT COLUMNS  -  FFT DOWN EACH OF LOT COLUMNS
C         CALL FFTC(DATA,WORK,INC=2,JUMP=2*LEN,LEN,LOT,-1)
C     2)  INPUT IS A MATRIX WITH LOT ROWS AND 2*LEN COLUMNS (COL 1 REAL,
C         COL 2 IMAG )
C         CALL FFTC(DATA,WORK,INC=2*LOT,JUMP=2,LEN,LOT,-1)
C
C     CONVERSION FOR CRAY FFTS IS: (REPLACE C999 WITH FFT)
C     C999(DATAR,DATAI,WORK,TRIGS,IFAX,INC=2,JUMP=2*LEN,LEN,LOT,IDIR)
C     FFTC(DATAR,DATAI,WORK,INC=2,JUMP=2*LEN,LEN,LOT,IDIR)
C
C HISTORY
C
C   MP01     ,  09/20/ 89,  M PORTNEY,  INITIAL RELEASE.
C   SS02     ,  29/SEP/89,  SL SMITH ,  GENERIC VERSION
C   MP03     ,  10/NOV/90,  M PORTNEY   MAKE VECTOR LENGTH PARAMETER STMNT,
C                                       USE FFTCL TO FIND WORK ARRAY,
C                                       REMOVE SEPARATE ARRAY FOR TRIG FACTS
C
* END ******************************************************************
*     PORTABLE VERSION
*     LENV MUST MATCH IN FFTC1L,FFTC1,FFTCL,FFTC,VC999
      PARAMETER (lenv = 8)
      DIMENSION DATA(*),WORK(*),IFAX(13)
      IT = 1 + 2*LEN
      CALL VFTFAX(LEN,IFAX,WORK)
      CALL VC999(DATA,WORK(IT),WORK,IFAX,INC,JUMP,LEN,LOT,IDIR)
C     CALL GENERC('FFTC','WORKING')
      RETURN
      END
      SUBROUTINE FACTR ( N, IP, NP)
      INTEGER ISCR(300), IP(6), NP(6), IPP(300), NPP(300)
      NMAX = SQRT(REAL(N))
      DO 1 I = 1,6
            IP(I) = 0     ! PRIME FACTORS
            NP(I) = 0     ! NUMBER OF OCCURRENCES
    1 CONTINUE
      CALL SIEVE ( NMAX, ISCR, IPP, NPR )
      DO 2 I = 1,NPR
            NPP(I) = 0
    2 CONTINUE
      NREDUC = N
      LFACT = 0
    5 CONTINUE
      DO 20 I = 1,NPR
            IMOD = MOD(NREDUC,IPP(I))
            IF (IMOD .EQ. 0.) THEN
                 NPP(I) = NPP(I) + 1     ! NUMBER OF OCCURS
                 NREDUC = NREDUC/IPP(I)          ! NEW NUMBER
                 IF (NREDUC .NE. 1) THEN
                      GO TO 5            ! FIND NEXT FACTOR
                 ELSE
                      GO TO 999          ! DONE
                 END IF
            END IF
            IF (I .EQ. NPR) THEN
                LFACT = NREDUC              ! NO MORE FACTORS
            END IF
   20 CONTINUE
  999 CONTINUE
      J = 0
      DO 40 I = 1,NPR
            IF (NPP(I) .NE. 0) THEN
                  J = J + 1          ! NUMBER OF PRIMES
                  NP(J) = NPP(I)     ! NUMBER OF OCCURS
                  IP(J) = IPP(I)          ! PRIMES
            END IF
   40 CONTINUE
      NFACT = J                         ! NUMBER OF PRIMES
      IF (LFACT .NE. 0) THEN
           NFACT = NFACT + 1
           NP(NFACT) = 1
           IP(NFACT) = LFACT
      END IF
      RETURN
      END
      SUBROUTINE FFTFLY(F,T,N,CS,NW,IP,LF,LFOLD,INC,JUMP,LOT)
*     COMPLEX F(*),T(*),CS(*)
*     COMPLEX COEF,COEF1,COEF2,COEF3
      real f(*),t(*),cs(*),coef,coef1,coef2,coef3,coef4
      IT=1
      LT0=1
      DO 400 JF=1,NW
         JT=IT
         NDX=LT0
         LT=LT0
         IF(IP.EQ.2)THEN
*           COEF=CS(NDX)
            coefr = cs(2*ndx-1)
            coefi = cs(2*ndx  )
            DO 210 K=JT,JT+LF-1
               I1 = (K+LT0-JT-1)*INC
               I2 = (K-1)*INC
               I3 = (K+LF-1)*INC
               DO 211 I=1,LOT
*                  T(I1+I)=F(I2+I)+COEF*F(I3+I)
                   t(2*(i1+i)-1) = f(2*(i2+i)-1) + coefr*f(2*(i3+i)-1)
     :                                           - coefi*f(2*(i3+i)  )
                   t(2*(i1+i)  ) = f(2*(i2+i)  ) + coefi*f(2*(i3+i)-1)
     :                                           + coefr*f(2*(i3+i)  )
  211          CONTINUE
  210       CONTINUE
         ELSEIF(IP.EQ.3)THEN
*           COEF=CS(NDX)
*           COEF1=COEF*COEF
            coefr = cs(2*ndx-1)
            coefi = cs(2*ndx  )
            coef1r = coefr*coefr - coefi*coefi
            coef1i = coefr*coefi + coefi*coefr
               DO 220 K=JT,JT+LF-1
               I1 = (K+LT0-JT-1)*INC
               I2 = (K-1)*INC
               I3 = (K+LF-1)*INC
               I4 = (K+LF+LF-1)*INC
                  DO 221 I=1,LOT
*                    T(I1+I)=F(I2+I)+COEF*F(I3+I)+COEF1*F(I4+I)
                    t(2*(i1+i)-1) = f(2*(i2+i)-1) + coefr*f(2*(i3+i)-1)
     :                                            - coefi*f(2*(i3+i)  )
     :                                            + coef1r*f(2*(i4+i)-1)
     :                                            - coef1i*f(2*(i4+i)  )
                    t(2*(i1+i)  ) = f(2*(i2+i)  ) + coefi*f(2*(i3+i)-1)
     :                                            + coefr*f(2*(i3+i)  )
     :                                            + coef1i*f(2*(i4+i)-1)
     :                                            + coef1r*f(2*(i4+i)  )
  221             CONTINUE
  220          CONTINUE
         ELSEIF(IP.EQ.5)THEN
*           COEF=CS(NDX)
*           COEF1=COEF*COEF
*           COEF2=COEF*COEF1
*           COEF3=COEF*COEF2
            coefr = cs(2*ndx-1)
            coefi = cs(2*ndx  )
            coef1r = coefr*coefr - coefi*coefi
            coef1i = coefr*coefi + coefi*coefr
            coef2r = coefr*coef1r - coefi*coef1i
            coef2i = coefr*coef1i + coefi*coef1r
            coef3r = coefr*coef2r - coefi*coef2i
            coef3i = coefr*coef2i + coefi*coef2r
            DO 530 K=JT,JT+LF-1
               I1 = (K+LT0-JT-1)*INC
               I2 = (K-1)*INC
               I3 = (K+LF-1)*INC
               I4 = (K+LF+LF-1)*INC
               I5 = (K+LF+LF+LF-1)*INC
               I6 = (K+LF+LF+LF+LF-1)*INC
               DO 531 I=1,LOT
*                 T(I1+I)=F(I2+I)+COEF*F(I3+I)+COEF1*F(I4+I)+
*    :                    COEF2*F(I5+I)+COEF3*F(I6+I)
                  t(2*(i1+i)-1) = f(2*(i2+i)-1) + coefr*f(2*(i3+i)-1)
     :                                          - coefi*f(2*(i3+i)  )
     :                                          + coef1r*f(2*(i4+i)-1)
     :                                          - coef1i*f(2*(i4+i)  )
     :                                          + coef2r*f(2*(i5+i)-1)
     :                                          - coef2i*f(2*(i5+i)  )
     :                                          + coef3r*f(2*(i6+i)-1)
     :                                          - coef3i*f(2*(i6+i)  )
                  t(2*(i1+i)  ) = f(2*(i2+i)  ) + coefi*f(2*(i3+i)-1)
     :                                          + coefr*f(2*(i3+i)  )
     :                                          + coef1i*f(2*(i4+i)-1)
     :                                          + coef1r*f(2*(i4+i)  )
     :                                          + coef2i*f(2*(i5+i)-1)
     :                                          + coef2r*f(2*(i5+i)  )
     :                                          + coef3i*f(2*(i6+i)-1)
     :                                          + coef3r*f(2*(i6+i)  )
  531          CONTINUE
  530       CONTINUE
         ENDIF
         IT=IT+LFOLD
         IF(IT.GT.N)IT=1
         LT0=LT0+LF
  400 CONTINUE
      RETURN
      END
      SUBROUTINE SIEVE (N, ISCR, IP, NP)
*
* Abstract *************************************************************
*
*   prime number generator
*
* Keywords
*
*   Fourier transform, math.
*
* Purpose
*
*   Finds prime numbers using sieve of Eratosthenes
*
* Arguments
*
*   N       INPUT         INTEGER           SCALAR
*     Input number
*
*   ISCR    SCRATCH       INTEGER           ARRAY(*)
*     Scratch array, length is N
*
*   IP      OUTPUT        INTEGER           ARRAY(*)
*     Array of primes, increasing order
*
*   NP      OUTPUT        INTEGER           ARRAY(*)
*     Number of primes found
*
* Common
*
*   None
*
* Errors
*
*   None
*
* Notes
*
*   None
*
* History
*
*   MNP      ,  04/12/88,  Mark Portney,  Initial release.
*
* End ******************************************************************
*
      INTEGER ISCR(*),IP(*)
      N = MAX0(N,3)
      DO 10 I = 1,N
            ISCR(I) = 0
   10 CONTINUE
      DO 20 I = 2,N-1
            DO 30 J =2*I,N,I
                  ISCR(J) = ISCR(J) + 1  ! NUMBER OF TIMES J IS FOUND
C                                    ! IN SIEVE
   30       CONTINUE
   20 CONTINUE
      J = 0
      DO 40 I = 2,N
            IF (ISCR(I) .EQ. 0) THEN
                  J = J + 1          ! NUMBER OF PRIMES
                  IP(J) = I          ! PRIMES
            END IF
   40 CONTINUE
      NP = J                         ! NUMBER OF PRIMES
      RETURN
      END
      SUBROUTINE VC64(F,T,CS,IFAX,INC,JUMP,N,LOT,IDIR)
      INTEGER IPRM(6),NPRM(6),IFAX(*)
      COMPLEX F(*),T(*),CS(*)
      SAVE
C     JUMP = 1 FOR EFFICIENCY
C     INC NOT EQ. 1 FOR EFFICIENCY
      DO 100 I=1,6
         IPRM(I)=IFAX(1+I)
         NPRM(I)=IFAX(7+I)
  100 CONTINUE
      JNVFWD=IDIR
      TEST = JNVFWD*AIMAG(CS(2))
      IF (TEST .GT. 0.) THEN
         DO 50 I = 1,N
            CS(I) = CONJG(CS(I))
   50    CONTINUE
      END IF
      IPC=0
      LF=N
      NW=1
      ISWT=1
  101 CONTINUE
         IPC=IPC+1
         IP=IPRM(IPC)
         NP=NPRM(IPC)
         DO 800 JP=1,NP
            LFOLD=LF
            LF=LF/IP
            NWOLD=NW
            NW=NW*IP
            IF(ISWT.GT.0)THEN
               CALL FFTFLY(F,T,N,CS,NW,IP,LF,LFOLD,INC,JUMP,LOT)
            ELSEIF(ISWT.LT.0)THEN
               CALL FFTFLY(T,F,N,CS,NW,IP,LF,LFOLD,INC,JUMP,LOT)
            ENDIF
            ISWT=-ISWT
  800    CONTINUE
      IF (LF.GT.1) GO TO 101
C
      IF (ISWT.LT.0) THEN
            LOT1 = 2*((LOT+1)/2)
            DO 951 I=1,N*LOT1
               F(I)=T(I)
  951       CONTINUE
      END IF
      RETURN
      END
      SUBROUTINE VC999(DATA,WORK,TRIG,IFAX,INC,JUMP,LEN,LOT,IDIR)
*     LENV MUST MATCH IN FFTC1L,FFTC1,FFTCL,FFTC,VC999
      PARAMETER (lenv = 8)
      DIMENSION DATA(*),WORK(*),TRIG(*)
      IW = 1+2*LEN*MIN0(LOT,LENV)
      LDIR = - IDIR
      DO 10 I = 1,LOT,LENV
         NLOT = MIN0(LOT-I+1,LENV)
         DO 20 J = 1,NLOT
            I1 = 1+(I+J-2)*JUMP
            I2 = 1+2*(J-1)
            DO 30 K = 1,LEN
               WORK(I2+(K-1)*NLOT*2) = DATA(I1+(K-1)*INC)
   30       CONTINUE
            DO 31 K = 1,LEN
               WORK(1+I2+(K-1)*NLOT*2) = DATA(1+I1+(K-1)*INC)
   31       CONTINUE
   20    CONTINUE
         CALL VC64(WORK(1),WORK(IW),TRIG,IFAX,NLOT,1,LEN,NLOT,LDIR)
         DO 40 J = 1,NLOT
            I1 = 1+(I+J-2)*JUMP
            I2 = 1+2*(J-1)
            DO 50 K = 1,LEN
               DATA(I1+(K-1)*INC) = WORK(I2+(K-1)*NLOT*2)
   50       CONTINUE
            DO 51 K = 1,LEN
               DATA(1+I1+(K-1)*INC) = WORK(1+I2+(K-1)*NLOT*2)
   51       CONTINUE
   40    CONTINUE
   10 CONTINUE
      RETURN
      END
      SUBROUTINE VFTFAX(N,IFAX,CS)
      INTEGER IFAX(*),IPRM(6),NPRM(6)
      COMPLEX*8 CS(*), COEF
      PI=3.14159265378
      PI2=PI*2.0
      DTH=-PI2/N
      C=1.
      S=0.
      CD=COS(DTH)
      SD=SIN(DTH)
      COEF=CMPLX(CD,SD)
      CS(1)=CMPLX(C,S)
      NSTEP=10
      DO 30 J=2,N,NSTEP
         IEND=J+NSTEP-1
         IF(IEND.GT.N)IEND=N
         DO 20 I=J,IEND-1
            CS(I)=COEF*CS(I-1)
   20    CONTINUE
         TH=-PI2*(IEND-1)/N
         CS(IEND)=CMPLX(COS(TH),SIN(TH))
   30 CONTINUE
      CALL FACTR ( N, IPRM, NPRM)
      DO 100 I=1,6
         IFAX(1+I)=IPRM(I)
         IFAX(7+I)=NPRM(I)
  100 CONTINUE
      RETURN
      END
