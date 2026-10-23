c     SEASEIS: copias thread-safe (OpenMP) das rotinas do laco de
c     frequencias do SPKBRW: xxmvphs2d, xpspivels1, bndry, fftl2d,
c     dwnn, fftc e auxiliares, movlev, clear. Mesmo codigo, com o
c     prefixo "o", sem SAVE; compilar com -fopenmp (implica
c     -frecursive) e SEM -fno-automatic.
      SUBROUTINE offtc(DATA,WORK,INC,JUMP,LEN,LOT,IDIR)
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
C         CALL offtc(DATA,WORK,INC=2,JUMP=2*LEN,LEN,LOT,-1)
C     2)  INPUT IS A MATRIX WITH LOT ROWS AND 2*LEN COLUMNS (COL 1 REAL,
C         COL 2 IMAG )
C         CALL offtc(DATA,WORK,INC=2*LOT,JUMP=2,LEN,LOT,-1)
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
      CALL ovftfax(LEN,IFAX,WORK)
      CALL ovc999(DATA,WORK(IT),WORK,IFAX,INC,JUMP,LEN,LOT,IDIR)
C     CALL GENERC('FFTC','WORKING')
      RETURN
      END
      SUBROUTINE ofactr ( N, IP, NP)
      INTEGER ISCR(300), IP(6), NP(6), IPP(300), NPP(300)
      NMAX = SQRT(REAL(N))
      DO 1 I = 1,6
            IP(I) = 0     ! PRIME FACTORS
            NP(I) = 0     ! NUMBER OF OCCURRENCES
    1 CONTINUE
      CALL osieve ( NMAX, ISCR, IPP, NPR )
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
      SUBROUTINE offtfly(F,T,N,CS,NW,IP,LF,LFOLD,INC,JUMP,LOT)
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
      SUBROUTINE osieve (N, ISCR, IP, NP)
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
      SUBROUTINE ovc64(F,T,CS,IFAX,INC,JUMP,N,LOT,IDIR)
      INTEGER IPRM(6),NPRM(6),IFAX(*)
      COMPLEX F(*),T(*),CS(*)
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
               CALL offtfly(F,T,N,CS,NW,IP,LF,LFOLD,INC,JUMP,LOT)
            ELSEIF(ISWT.LT.0)THEN
               CALL offtfly(T,F,N,CS,NW,IP,LF,LFOLD,INC,JUMP,LOT)
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
      SUBROUTINE ovc999(DATA,WORK,TRIG,IFAX,INC,JUMP,LEN,LOT,IDIR)
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
         CALL ovc64(WORK(1),WORK(IW),TRIG,IFAX,NLOT,1,LEN,NLOT,LDIR)
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
      SUBROUTINE ovftfax(N,IFAX,CS)
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
      CALL ofactr ( N, IPRM, NPRM)
      DO 100 I=1,6
         IFAX(1+I)=IPRM(I)
         IFAX(7+I)=NPRM(I)
  100 CONTINUE
      RETURN
      END
      subroutine offtl2d(y,a,nx,ny,isgn)
      complex y(*),a(*)
      if(nx .gt. 1)call offtc(y,a,2,2*nx,nx,ny,-isgn)
      if(ny .gt. 1)call offtc(y,a,2*nx,2,ny,nx,-isgn)
      scale=1./(real(nx)*real(ny))
      if(isgn.eq.-1)then
         do 1 j=1,nx*ny
             y(j)=scale*y(j)
    1    continue
      end if
      return
      end
      subroutine odwnn(b,nx,ny,dx,dy,w,shft,vv,sgn)
*     b has length 2*nx*ny
      dimension b(*)
      if (w .eq. 0.) then
         do 1 i  =  1,2*nx*ny
            b(i)  =  0.
    1    continue
         return
      end if
      dkx = 2.*3.141593/nx
      dky = 2.*3.141593/ny
      dkx = dkx*dy/dx
      t = sgn*w*shft
      vvrw2 = vv/w**2 
         k = 1
         ind1=3
         ind2=2*nx-1
         do 21 i = 2,nx/2+1
            akx=(i-1)*dkx
            xk2 = akx**2
            r = 1.-vvrw2*xk2
            if (r .le. 0.) then 
               b(ind1) = 0.
               b(ind1+1) = 0.
               b(ind2) = 0.
               b(ind2+1) = 0.
            else
               a = t*(sqrt(r)-1.)
               er = cos(a)
               ei = sin(a)
               t1  = er*b(ind1)-ei*b(ind1+1)
               b(ind1+1) = er*b(ind1+1)+ei*b(ind1)
               b(ind1) = t1
               t1  = er*b(ind2)-ei*b(ind2+1)
               b(ind2+1) = er*b(ind2+1)+ei*b(ind2)
               b(ind2) = t1
            end if
            ind1 = ind1 + 2
            ind2 = ind2 - 2
   21    continue
         if(ny .eq. 1)return
         k = ny/2+1
         ind = 1 + (k-1)*2*nx
         aky = min0(k-1,ny-k+1)*dky
         aky2 = aky**2
         do 22 i = 1,nx
            akx = min0(i-1,nx-i+1)*dkx
            xk2 = akx**2+aky2
            r = 1.-vvrw2*xk2
            if (r .le. 0.) then 
               b(ind) = 0.
               b(ind+1) = 0.
            else
               a = t*(sqrt(r)-1.)
               er = cos(a)
               ei = sin(a)
               t1  = er*b(ind)-ei*b(ind+1)
               b(ind+1) = er*b(ind+1)+ei*b(ind)
               b(ind) = t1
            end if
            ind = ind + 2
   22    continue
         do 24 k = 2,ny/2
            kk = ny-k+2
            aky = (k-1)*dky
            aky2 = aky**2
            i = 1
            ind = i + (k-1)*2*nx
            indd = i + (kk-1)*2*nx
            akx = min0(i-1,nx-i+1)*dkx
            xk2 = akx**2+aky2
            r = 1.-vvrw2*xk2
            if (r .le. 0.) then 
               b(ind) = 0.
               b(ind+1) = 0.
               b(indd) = 0.
               b(indd+1) = 0.
            else
               a = t*(sqrt(r)-1.)
               er = cos(a)
               ei = sin(a)
               t1  = er*b(ind)-ei*b(ind+1)
               b(ind+1) = er*b(ind+1)+ei*b(ind)
               b(ind) = t1
               t2 = er*b(indd)-ei*b(indd+1)
               b(indd+1) = er*b(indd+1)+ei*b(indd)
               b(indd) = t2
            end if
            i = nx/2+1
            ind = 2*i-1 + (k-1)*2*nx
            indd = 2*i-1 + (kk-1)*2*nx
            akx = min0(i-1,nx-i+1)*dkx
            xk2 = akx**2+aky2
            r = 1.-vvrw2*xk2
            if (r .le. 0.) then 
               b(ind) = 0.
               b(ind+1) = 0.
               b(indd) = 0.
               b(indd+1) = 0.
            else
               a = t*(sqrt(r)-1.)
               er = cos(a)
               ei = sin(a)
               t1  = er*b(ind)-ei*b(ind+1)
               b(ind+1) = er*b(ind+1)+ei*b(ind)
               b(ind) = t1
               t2 = er*b(indd)-ei*b(indd+1)
               b(indd+1) = er*b(indd+1)+ei*b(indd)
               b(indd) = t2
            end if
            i = 2
            ind = 2*i-1 + (k-1)*2*nx
            indd = 2*i-1 + (kk-1)*2*nx
            inb = 2*(nx-i+2)-1 + (k-1)*2*nx
            inbb = 2*(nx-i+2)-1 + (kk-1)*2*nx
CDIR$ IVDEP
         do 43 i = 2,nx/2
            akx = min0(i-1,nx-i+1)*dkx
            xk2 = akx**2+aky2
            r = 1.-vvrw2*xk2
            if (r .le. 0.) then 
               b(ind) = 0.
               b(ind+1) = 0.
               b(indd) = 0.
               b(indd+1) = 0.
               b(inb) = 0.
               b(inb+1) = 0.
               b(inbb) = 0.
               b(inbb+1) = 0.
            else
               a = t*(sqrt(r)-1.)
               er = cos(a)
               ei = sin(a)
               t1  = er*b(ind)-ei*b(ind+1)
               b(ind+1) = er*b(ind+1)+ei*b(ind)
               b(ind) = t1
               t2 = er*b(indd)-ei*b(indd+1)
               b(indd+1) = er*b(indd+1)+ei*b(indd)
               b(indd) = t2
               t3  = er*b(inb)-ei*b(inb+1)
               b(inb+1) = er*b(inb+1)+ei*b(inb)
               b(inb) = t3
               t3 = er*b(inbb)-ei*b(inbb+1)
               b(inbb+1) = er*b(inbb+1)+ei*b(inbb)
               b(inbb) = t3
            end if
            ind = ind + 2
            indd = indd + 2
            inb = inb - 2
            inbb = inbb - 2
   43    continue
   24    continue
      return
      end
      subroutine omovlev(a,b,n)
      dimension a(*),b(*) 
      do 1 i=1,n
1     b(i)=a(i) 
      return
      end
      subroutine oclear(a,n) 
      dimension a(*)
      do 10 i=1,n 
10     a(i)=0.
      return
      end
      subroutine oxxmvphs2d(y,a,b,rs,nx,ny,dx,dy,dt,w,shft,vv,sgn,nvels)
      complex y(*),a(*),b(*) 
      real shft(*),vv(*),rs(*)
      real s1(20)
c   ***   find min and max shifts 
      nxy=nx*ny
      v1=0.
      v2=1000000.
      shft1=1000000.
      shft2=0.
      do 5 i=1,nxy
      v1=amax1(v1,vv(i))
      v2=amin1(v2,vv(i))
      shft1=amin1(shft1,shft(i))
      shft2=amax1(shft2,shft(i))
5     continue
      vmax=1./sqrt(v2)
      vmin=1./sqrt(v1)
      shftmin=shft1
      shftmax=shft2
c   ***   fft, downward continue, and inverse fft for vmin
c   ***   fft, downward continue, and inverse fft for shift=vmin
      call offtl2d(y,b,nx,ny,+1)
c   ***   copy field to scratch array 
      do 7 i=1,nxy
7     a(i)=y(i) 
c   ***   if nvels = 1 the v1 = ( .5*(1./vmin+1./vmax) )**2
c   ***                shft1 = .5*(shftmin+shftmax)
      if(nvels .eq. 1)then
      v1=(.5*(1./vmin + 1./vmax) )**2
      shft1=.5*(shftmin+shftmax)
      end if
      call odwnn(y,nx,ny,dx,dy,w,shft1,v1,sgn)
      call offtl2d(y,b,nx,ny,-1)
c   ***   fft, downward continue, and inverse fft for shift=vmax
      if(vmax-vmin.lt..02*vmax .or. nvels .eq. 1)return
c   ***   velocity non-constant in this interval
c   ***   apply multi-velocity pspi
c   ***   make vmax slightly  larger than max of 1./sqrt(vv)
c   ***   make vmin slightly smaller than min of 1./sqrt(vv)
         dxz=vmin/shftmin
        vmax=1.0001*vmax
        vmin=.9999*vmin
      ichk=0
      dz=nint( shft1*dx*sqrt(v1))
      call oxpspivels1(shft,rs,dz,dt,shftmin,shftmax,nxy,nvels1,ichk)
      shftinc=(shftmax-shftmin)/(nvels-1)
       if(nvels1.le.nvels)then
       do j=1,nvels1
       s1(j)=rs(j)
       end do
       else
       nvels1=nvels
       do j=1,nvels1
       s1(j)=shftmin+(j-1)*shftinc
       end do
       end if
      nsteps=nvels1-1
      call omovlev(y,b,2*nxy)
      call oclear(y,2*nxy)
c   ***   start multiple velocity loop
        vmin=.999*dxz*shftmin
        vmax=1.001*dxz*shftmax
      do istep=1,nsteps
      xs1=(s1(istep)-s1(1))/(shftmax-shftmin)
      xs2=(s1(istep+1)-s1(1))/(shftmax-shftmin)
      shft1=s1(istep)
      shft2=s1(istep+1)
      vv1=vmin+xs1*(vmax-vmin)
      vv2=vmin+xs2*(vmax-vmin)
      v1 = (1./vv1)**2
      v2 = (1./vv2)**2
c   ***   compute contribution of field downward continued using vv1
      do k=1,nxy
      rs(k) = ( vv2-1./sqrt(vv(k)) )/(vv2-vv1)
      rs(k) = amax1(0.,rs(k))
      if(rs(k) .gt. 1.)rs(k)=0.
      end do
c   ***   add contrubution of field downward continued using vv1
      do k=1,nxy
      y(k) = y(k) + rs(k)*b(k)
      end do
c   ***
      call omovlev(a,b,2*nxy)
      call odwnn(b,nx,ny,dx,dy,w,shft2,v2,sgn)
c   ***   inverse fft data and scale
      call offtl2d(b,rs,nx,ny,-1)
c   ***   compute contribution of field downward continued using vv2
      do k=1,nxy
      rs(k)=(1./sqrt(vv(k))-vv1)/(vv2-vv1)
      rs(k)=amax1(0.,rs(k))
      if(rs(k) .ge. 1.)rs(k)=0.
      end do
c   ***   add contribution of field downward continued using vv2
      do k=1,nxy
      y(k) = y(k) + rs(k)*b(k)
      end do
c   ***
      end do
c   ***
      return
      end 
      subroutine oxpspivels1(s1,s2,dz,dt,smin,smax,nx,nvels,icheck)
      real s1(*),s2(*)
      real vlist(2001)
      integer is(2001)
c   ***   clear  index array
c   ***   if smin nearly equal smax then set
c   ***   nvels=1    and   s2(1)=.5*(smin+smax)
         if(smax-smin .lt. .01*smax )then
         nvels=1
         s2(1)=.5*(smin+smax)
         return
         end if
c   ***
      ns=2001
      do j=1,ns
      is(j)=0
      vlist(j)=20000.-10.*(j-1)
      end do
       dzt=dz/dt
c   ***   set index array for histogram of velocity profile
      do j=1,nx
      vv=dzt/s1(j)
       ii=1+nint( (20000.-vv)*.1 )
       is(ii)=is(ii)+1
      end do
c   ***   compute histogram of velocity profile
      jj=0
       imax=0
      do j=2,ns
      imax=max0(is(j),imax)
      end do
      npct=imax/10
      vmn=10.*nint(.1*dzt/smax)
      vmx=10.*nint(.1*dzt/smin)
      jvmn=1+nint( (20000.-vmn)*.1)
      jvmx=1+nint( (20000.-vmx)*.1)
      do j=2,ns
      if(is(j).ge.npct .or. j.eq.jvmn .or. j.eq.jvmx)then
      jj=jj+1
      s2(jj)=dzt/vlist(j)
      end if
      end do
      nvels=jj
      s2(1)=smin
      s2(nvels)=smax
102   format(2i4,8f8.0)
103   format(8i10)
      return
      end
      subroutine obndry(y,nx,ny,lbx,lby)
      complex y(*)
c   ***   apply attenuating absorbing boundary
c   ***   y     nx by ny complex field for one frequency 
c   ***   nx    # xlines
c   ***   ny    # ilines
c   ***   lbx   # of traces to attenuate at each end of inlines
c   ***   lby   # of traces to attenuate at each end of xlines
c   ***   lby trace attenuation region
c   ***   attenuate x boundary regions
      do iy=1,ny
         ii=(iy-1)*nx
         do j=1,lbx
            xx=(99.-lbx+j)/100.
            y(ii+j)=xx*y(ii+j)
         end do
         ii=ii+nx-lbx
         do j=1,lbx
            xx=(100.-j)/100.
            y(ii+j)=xx*y(ii+j)
         end do
      end do
c   ***   end x boundary attenuation
c   ***   attenuate y boundary region
      if(ny .gt. 30)then
      do ix=1,nx
         ii=ix
         do j=1,lby
            xx=(99.-lby+j)/100.
            y(ii+(j-1)*nx)=xx*y(ii+(j-1)*nx)
         end do
         ii=(ny-lby)*nx+ix
         do j=1,lby
            xx=(100.-j)/100
            y(ii+(j-1)*nx)=xx*y(ii+(j-1)*nx)
         end do
      end do
      end if
      return
      end
      subroutine osgimrg2r(z,u,y,shft,t1,t2,nl,nh,nt1,nx,itzr) 
c     SEASEIS: sgimrg2r com OpenMP. Os pontos i sao processados em
c     blocos de NB (laco de frequencias por fora do laco em i, para
c     usar a cache); para cada i as contas e a ordem da soma sao as
c     mesmas do original, logo o resultado e identico.
      parameter (nb=128)
      real u(*),y(*),z(*),shft(*),t1(*),t2(*)
      real q(12,nb)
      pi=3.1415927
      rz=1./itzr
      nws=nh-nl+1 
      delw=2.*pi/nt1
      wf=(nl-1)*delw
!$omp parallel do default(shared) schedule(static)
!$omp& private(i0,ie,i,m,q,k,ar,ai,br,bi,cr,ci,ll,jr,ji,l,f1,f2,tmp)
      do 100 i0=1,nx,nb
      ie=min0(nx,i0+nb-1)
      do l=0,itzr-1
      do i=i0,ie
      z(l*nx+i)=0.
      end do
      end do
c    ***    initialize frequencies
      do i=i0,ie
      m=i-i0+1
      q(1,m)=cos(delw*t1(i))
      q(2,m)=cos(wf*t1(i)) 
      q(3,m)=cos(delw*(t1(i)-shft(i))) 
      q(4,m)=cos(wf*(t1(i)-shft(i))) 
      q(5,m)=cos(delw*rz*(t2(i)-t1(i)+shft(i)))  
      q(6,m)=cos(wf*rz*(t2(i)-t1(i)+shft(i))) 
      q(7,m)=sin(delw*t1(i))
      q(8,m)=sin(wf*t1(i)) 
      q(9,m)=sin(delw*(t1(i)-shft(i))) 
      q(10,m)=sin(wf*(t1(i)-shft(i))) 
      q(11,m)=sin(delw*rz*(t2(i)-t1(i)+shft(i)))  
      q(12,m)=sin(wf*rz*(t2(i)-t1(i)+shft(i))) 
      end do
c     q: 1 sr1  2 sr2  3 sr3  4 sr4  5 sr5  6 sr6
c        7 si1  8 si2  9 si3 10 si4 11 si5 12 si6
      do 20 k=1,nws
      do 30 i=i0,ie
      m=i-i0+1
      ar=q(2,m)
      ai=q(8,m)
      br=q(4,m)
      bi=q(10,m)
      cr=q(6,m)
      ci=q(12,m)
      ll=0
      jr=2*((k-1)*nx+i)-1
      ji=jr+1 
      do 10 l=0,itzr-1
      f1=(itzr-l)*rz
      f2=1.-f1
      z(ll+i)=z(ll+i)+f1*(ar*u(jr)-ai*u(ji))+f2*(br*y(jr)-bi*y(ji))
      tmp=ar*cr-ai*ci
      ai=ar*ci+ai*cr
      ar=tmp
      tmp=br*cr-bi*ci
      bi=br*ci+bi*cr
      br=tmp
      ll=ll+nx
10    continue
c   ***   update next frequency  
      tmp=q(2,m)*q(1,m)-q(8,m)*q(7,m)
      q(8,m)=q(2,m)*q(7,m)+q(8,m)*q(1,m)
      q(2,m)=tmp
      tmp=q(4,m)*q(3,m)-q(10,m)*q(9,m)
      q(10,m)=q(4,m)*q(9,m)+q(10,m)*q(3,m)
      q(4,m)=tmp
      tmp=q(6,m)*q(5,m)-q(12,m)*q(11,m)
      q(12,m)=q(6,m)*q(11,m)+q(5,m)*q(12,m)
      q(6,m)=tmp
30    continue
20    continue
100   continue  
!$omp end parallel do
      return
      end
      subroutine osgivshift1(y,shft,dw,nl,nh,nx,sgn)
      real y(*)
      real shft(*)
c   ***   efficient  complex exponential table computation
      w1=sgn*(nl-1)*dw
      dws=sgn*dw
      nws=nh-nl+1
      nx2=2*nx
c     SEASEIS: OpenMP nos pontos i
!$omp parallel do default(shared) schedule(static)
!$omp& private(i,sr,si,ar,ai,ii,j,t1)
      do 10 i=1,nx
c     s1=cexp(cmplx(0.,dws*shft(i)))
c     a=cexp(cmplx(0.,w1*shft(i)))
      sr=cos(dws*shft(i))
      si=sin(dws*shft(i))
      ar=cos( w1*shft(i))
      ai=sin( w1*shft(i))
      ii=2*(i-1)+1
      do 20 j=1,nws
c     y(ii+i)=a*y(ii+i)
c     a=s1*a
      t1=ar*y(ii)-ai*y(ii+1)
      y(ii+1)=ar*y(ii+1)+ai*y(ii)
      y(ii)=t1
      t1=ar*sr-ai*si
      ai=ar*si+ai*sr
      ar=t1
      ii=ii+nx2
20    continue
10    continue
!$omp end parallel do
      return
      end
