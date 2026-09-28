c-----------------------------------------------------------------------
c     Rotinas copiadas sem alteracao de dvslib.f (SEISPAK), exceto onde
c     marcado com "SEASEIS:". Usadas pelo modulo SPK_DVSSYN3D.
c-----------------------------------------------------------------------
      subroutine clear(a,n) 
      dimension a(*)
      do 10 i=1,n 
10     a(i)=0.
      return
      end
      subroutine consts(theta,vv,vz,alfa,gamma,beta,aa,bb)
      theta2=theta
      pi=3.14159
      theta2=theta2*pi/180. 
      akxdx=sin(theta2)/sqrt(vv)
      if(akxdx.lt..8*pi)go to 1 
      arcs=asin(.8*pi*sqrt(vv)) 
      theta2=amin1(theta2,arcs) 
1     continue
      theta1=.5*theta2
      sqrtvv=sqrt(vv) 
      rsqrtvv=1./sqrtvv 
      akxdx1=sin(theta1)*rsqrtvv
      akxdx2=sin(theta2)*rsqrtvv
      t1=2.*(1.-cos(akxdx1))
      t2=2.*(1.-cos(akxdx2))
      chi1=.5*rsqrtvv/tan(.5*akxdx1)
      chi2=.5*rsqrtvv/tan(.5*akxdx2)
      eta1=vz/tan(vz*(cos(theta1)-1.))
      eta2=vz/tan(vz*(cos(theta2)-1.))
      aa=-(chi2-chi1)/(eta2-eta1) 
      bb=eta1*aa+chi1 
      sinsq1=sin(theta1)**2 
      sinsq2=sin(theta2)**2 
c     beta=(sinsq2-vv*t2)/(sinsq2*t2) 
      rho1=1./(vv*t1/(1.-beta*t1))
      rho2=1./(vv*t2/(1.-beta*t2))
      gamma=-(rho2-rho1)/(eta2-eta1)
      alfa=eta1*gamma+rho1
      return
      end
      subroutine ctran(a,b,nr,nc) 
      complex a(*),b(*) 
      do 10 k=1,nc
      do 10 i=1,nr
      b((i-1)*nc+k)=a((k-1)*nr+i)
10    continue
      return
      end
      subroutine depth2w(y,u,a,b,shft,nl,nh,nt1,itzr,nx) 
      complex a(*),b(*),y(*) 
      dimension u(*),shft(*) 
      nws=nh-nl+1  
      dw=2.*3.141593/nt1 
      w1=(nl-1)*dw 
      ritzr=1./itzr  
      do 10 l=1,itzr 
      ll=(l-1)*nx  
      r=(itzr-l+1)*ritzr
      do 20 k=1,nx 
      a(k)=cexp(cmplx(0.,r*w1*shft(k)))  
20    b(k)=cexp(cmplx(0.,r*dw*shft(k)))  
      do 30 i=1,nws  
      ik=(i-1)*nx  
      do 30 k=1,nx 
      y(ik+k)=y(ik+k)+u(ll+k)*a(k) 
30    a(k)=a(k)*b(k) 
10    continue 
      return 
      end
      subroutine dwn3d(y,yy,s,al,bl,cl,w,v,shft,nx,ny,sgn) 
c   ***   single step of split step 3-d downward/upward continuation.  
c   ***   ny by nx plate at frequency w (data sequential in y) 
c   ***   sgn=+1. for downward/sgn=-1. for upward   continuation.  
      real v(*),shft(*) 
      complex s(*),al(*),bl(*),cl(*),y(*),yy(*)  
      raa=1./6.  
      theta=60.  
      nx1=nx-1 
      rw=0.  
      rw2=0. 
      aa=.25 
      cc=.5  
      aaa=1./12. 
      bbb=4./3.  
      if(w.ne.0.)then  
      rw=1./w  
      rw2=rw**2  
      con1=rw2 
      con2=.5*w  
      ii=nx*(ny/2)+nx/2  
      vv=v(ii)*con1  
      vz=shft(ii)*con2 
      if(v(ii).lt..001)then  
      vv=con1  
      vz=con2  
      end if 
c     SEASEIS: removida a chamada 'call fftl2d(y,b,nx,ny,+1)' (sobra de
c     copia do phs2d: 'b' nao declarado, FFT sem inversa no meio do passo
c     de diferencas finitas em x-y; travava ou embaralhava o campo)
      call consts(theta,vv,vz,aa,cc,raa,aaa,bbb) 
      end if 
      bbb=sgn*bbb  
      cc=sgn*cc  
c   ***   compute coefficients for boundary. 
c   ***   compute interior matrix coefficients.  
      c1=aa*rw2  
      c2=.5*cc*rw  
c   ***   compute coefficients for boundary. 
      do 101 j=1,ny  
      al(j)=cmplx(1.-aaa*sqrt(v(j))*shft(j),-2.*bbb*rw*sqrt(v(j)))
      bl(j)=2.-al(j) 
101   cl(j)=al(j)  
c   ***   compute interior matrix coefficients.  
      do 102 j=ny+1,nx1*ny 
      al(j)=raa+cmplx(c1*v(j),-c2*v(j)*shft(j))  
      bl(j)=1.-2.*al(j)  
102   cl(j)=al(j)  
c   ***   compute coefficients for boundary. 
      do 103 j= nx1*ny+1,nx*ny 
      al(j)=cmplx(1.-aaa*sqrt(v(j))*shft(j),-2.*bbb*rw)
      bl(j)=2.-al(j) 
103   cl(j)=al(j)  
c   ***   rhs = conjg of coef matrix times input field 
      do 10 j=1,ny 
10    s(j)=conjg(al(j))*y(j)+conjg(bl(j))*y(ny+j)  
      do 15 j=ny+1,nx1*ny  
15    s(j)=conjg(bl(j))*y(j)+conjg(al(j))*(y(-ny+j)+y(ny+j)) 
      do 20 j=nx1*ny+1,nx*ny 
20    s(j)=conjg(al(j))*y(j)+conjg(bl(j))*y(-ny+j) 
c   ***   simultaneous solution of systems by gaussian elimination 
      do 200 j=1,ny  
200   bl(j)=1./bl(j) 
      do 210 j=1,ny  
      cl(j)=bl(j)*cl(j)  
210   y(j)=bl(j)*s(j)  
      do 220 i=2,nx  
      m=(i-1)*ny 
      mi=m-ny  
      do 230 j=m+1,m+ny  
230   bl(j)=1./(bl(j)-al(j)*cl(-ny+j)) 
      do 240 j=m+1,m+ny  
      cl(j)=bl(j)*cl(j)  
240   y(j)=bl(j)*(s(j)-al(j)*yy(-ny+j))  
220   continue 
      m=nx1*ny 
      do 250 i=nx1,1,-1  
      m=(i-1)*ny 
      mi=m+ny  
      do 250 j=m+1,m+ny  
250   y(j)=y(j)-cl(j)*yy(ny+j) 
      return 
      end
      subroutine dwn45b(y,yy,s,al,bl,cl,w,w2,v,shft,nx,nl,nh,dw,sgn)
      dimension v(*),shft(*),w(*),w2(*) 
      complex s(*),al(*),bl(*),cl(*),y(*),yy(*),ar1,arn,br1,brn 
      raa=1./6. 
      theta=60.
      nx1=nx-1
      nws=nh-nl+1 
      rw=0. 
      rw2=0.
      if(nl.gt.1)then 
      rw=1./((nl-1)*dw) 
      rw2=rw**2 
      end if
      w(1)=rw 
      w2(1)=rw2 
      do 5 j=2,nws
5     w(j)=1./((nl+j-2)*dw) 
      do 6 j=2,nws
6     w2(j)=w(j)**2 
c   ***   compute boundary coefficients and accuracy constants
      ii=nx1*nws
      vv1=sqrt(v(1))
      vvn=sqrt(v(nx))
      do 221 j=1,nws
      vv=v(nx/2)*w2(j)
      vz=.5*shft(nx/2)/w(j)
      call consts(theta,vv,vz,aa,cc,raa,aaa,bbb)
      bbb=sgn*bbb
      cc=sgn*cc
      al(j)=cmplx(1.-aaa*vv1*shft(1),-2.*bbb*vv1*w(j))
      al(ii+j)=cmplx(1.-aaa*vvn*shft(nx),-2.*bbb*vvn*w(j))
      w(j)=.5*w(j)*cc
      w2(j)=w2(j)*aa
221   continue
      do 223 j=1,nws
      bl(j)=cmplx(2.-real(al(j)),-aimag(al(j)))
      cl(j)=al(j)
223   continue
      do 224 j=ii+1,ii+nws
      bl(j)=cmplx(2.-real(al(j)),-aimag(al(j)))
      cl(j)=al(j)
224   continue
c   ***   compute interior matrix coefficients
      do 100 i=2,nx-1
      ii=(i-1)*nws
      do 100 j=1,nws
      al(ii+j)=cmplx(raa+v(i)*w2(j),-w(j)*v(i)*shft(i))
      bl(ii+j)=cmplx(1.-2.*real(al(ii+j)),-2.*aimag(al(ii+j)))
      cl(ii+j)=al(ii+j)
100   continue
c   ***   rhs = conjg of coef matrix times input field
      ii=nx1*nws
      do 10 j=1,nws 
      s(j)=conjg(al(j))*y(j)+conjg(bl(j))*y(nws+j)
10    continue
      do 9 j=1,nws 
      s(ii+j)=conjg(al(ii+j))*y(ii+j)+conjg(bl(ii+j))*y(ii-nws+j)
9      continue
      do 15 i=2,nx1
      ii=(i-1)*nws
      do 15 j=ii+1,ii+nws
      s(j)=conjg(bl(j))*y(j)+conjg(al(j))*(y(j-nws)+y(j+nws))
15    continue
c   ***   simultaneous solution of systems by gaussian elimination
      do 200 j=1,nws
      w(j)=1./(real(bl(j))**2+aimag(bl(j))**2)
      bl(j)=cmplx(w(j)*real(bl(j)),-w(j)*aimag(bl(j)))
200   continue
      do 210 j=1,nws
      cl(j)=bl(j)*cl(j) 
      y(j)=bl(j)*s(j)
210   continue 
      do 220 i=2,nx 
      ii=(i-1)*nws
      iii=ii-nws
      do 230 j=1,nws
      bl(ii+j)=bl(ii+j)-al(ii+j)*cl(iii+j)
      w(j)=1./(real(bl(ii+j))**2+aimag(bl(ii+j))**2)
      bl(ii+j)=cmplx(w(j)*real(bl(ii+j)),-w(j)*aimag(bl(ii+j)))
230   continue
      do 240 j=1,nws
      cl(ii+j)=bl(ii+j)*cl(ii+j)
240   y(ii+j)=bl(ii+j)*(s(ii+j)-al(ii+j)*yy(iii+j)) 
220   continue
      ii=nx1*nws
      do 250 i=nx1,1,-1 
      ii=(i-1)*nws
      iii=ii+nws
      do 250 j=1,nws
250   y(ii+j)=y(ii+j)-cl(ii+j)*yy(iii+j)
500   continue     
c   ***   apply shift 
      call vshift(y,al,al,s,shft,dw,nl,nh,nx,sgn)
      return
      end
      subroutine dwnn(b,nx,ny,dx,dy,w,shft,vv,sgn)
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
      subroutine fac235(nx,ifax)
c   ***   subroutine fac235 returns in ifax(1) the nearest even integer 
c   ***   >= nx whose only nontrivial factors are 2 3 or 5
c   ***   ifax(2), ifax(3), and ifax(4) contain the number of factors of 2
c   ***   3, and 5 respectively
      integer ifax(*)
      nxx=((nx-1)/2)*2
200   continue
      nxx=nxx+2
      i2=0
      i3=0
      i5=0
      n1=nxx
      n3=n1
300   continue
      n2=n1
      if(mod(n1,2).eq.0)then
      i2=i2+1
      n1=n1/2
      n3=2*n3
      end if
      if(mod(n1,3).eq.0)then
      i3=i3+1
      n1=n1/3
      n3=3*n3
      end if
      if(mod(n1,5).eq.0)then
      i5=i5+1
      n1=n1/5
      n3=5*n3
      end if
c   ***   if n1=1 then finished
      if(n1.eq.1)goto 1000
c   ***   if n1 > 1 but n1=n2 then nxx not factorable by 2,3,5
      if(n1.gt.1 .and. n1.eq.n2)goto 200
c   ***   if n1 > 1 and n1 < n2 then a 2 3 or 5 factor was found
c   ***   make another pass with this nxx
      if(n1.gt.1 .and. n1.lt.n2)goto 300
1000   continue
      ifax(1)=nxx
      ifax(2)=i2
      ifax(3)=i3
      ifax(4)=i5
      return
      end
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
*
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
*
*
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
      subroutine joehint(x,y,xx,yy,dy,mm,nn,mode)  
      dimension x(*),y(*),xx(*),yy(*),dy(*)  
      mm1=mm-1 
c   ***   compute slopes at end points 
      mmm=mm1  
      if(mode.eq.-2)mmm=1  
      do 101 i=1,mmm 
      dy(i)=(y(i+1)-y(i))/(x(i+1)-x(i))  
101   continue 
      dy(mm)=(y(mm)-y(mm-1))/(x(mm)-x(mm-1)) 
      if(mm.lt.3)go to 10  
      if(mode.eq.-3)goto 10  
      if(mode.eq.0)goto 10 
c   ***   compute slopes at interior points  
      do 100 i=2,mm1 
      x0=x(i)  
      xm=x(i-1)-x0 
      xp=x(i+1)-x0 
      fm=y(i-1)  
      f0=y(i)  
      fp=y(i+1)  
      dy(i)=0. 
      pm=fm-f0 
      pp=fp-f0 
      wm=xm**2+pm**2 
      wp=xp**2+pp**2 
      xm2=xm**2/wm 
      xm3=xm*xm2 
      xm4=xm*xm3 
      xm5=xm*xm4 
      xp2=xp**2/wm 
      xp3=xp*xp2 
      xp4=xp*xp3 
      xp5=xp*xp4 
      a1=.25*(xm4-xp4) 
      b1=.333333*(xm3-xp3) 
      g1=.333333*(pm*xm2-pp*xp2) 
      a2=.20*(xm5-xp5) 
      b2=a1  
      g2=.25*(pm*xm3-pp*xp3) 
      den=a2*b1-a1*b2  
      top=a2*g1-a1*g2  
      if(den*1000000..le.top)go to 11  
      dy(i)=top/den  
11    continue 
100   continue 
10    continue 
c   ***   interpolate  
      j=1  
      ii=1 
      do 30 i=1,nn 
3     continue 
      if(ii.ne.i)go to 31  
      delx=x(j+1)-x(j) 
      eps=.000001*delx 
c   ***   compute polynomial coefficients  
      dxr=1./delx**2 
      a1=y(j)  
      a2=dy(j) 
      a3=((y(j+1)-a1)-(x(j+1)-x(j))*a2)*dxr  
      a4=((dy(j+1)-a2)-2.*(x(j+1)-x(j))*a3)*dxr  
      if(mode.eq.-2) goto  31  
      a3=0.  
      a4=0.  
31    continue 
      del1=xx(i)-x(j)  
      del2=x(j+1)-xx(i)  
      if(del1+eps.ge.0..and.del2+eps.ge.0.)go to 2 
      if(del1.lt.0..and.j.eq.1)go to 2 
      if(j.eq.mm1)go to 2  
      j=j+1  
      ii=i 
      go to 3  
2     continue 
      yy(i)=a1+del1*(a2+del1*(a3-del2*a4)) 
30    continue 
4     continue 
      return 
      end
      subroutine movlev(a,b,n)
      dimension a(*),b(*) 
      do 1 i=1,n
1     b(i)=a(i) 
      return
      end
      subroutine mvphs2d(y,a,b,rs,nx,ny,dx,dy,w,shft,vv,sgn,nvels)
      complex y(*),a(*),b(*) 
      real shft(*),vv(*),rs(*)
c   ***   find min and max shifts 
      nxy=nx*ny
      v1=0.
      v2=1000000.
      shft1=1000000.
      shft2=-1000000.
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
      call fftl2d(y,b,nx,ny,+1)
c   ***   copy field to scratch array 
      do 7 i=1,nxy
7     a(i)=y(i) 
c   ***   if nvels = 1 the v1 = ( .5*(1./vmin+1./vmax) )**2
c   ***                shft1 = .5*(shftmin+shftmax)
      if(vmax-vmin.lt..005*vmax .or. nvels .eq. 1)then
      v1=(.5*(1./vmin + 1./vmax) )**2
      shft1=.5*(shftmin+shftmax)
      end if
      call dwnn(y,nx,ny,dx,dy,w,shft1,v1,sgn)
      call fftl2d(y,b,nx,ny,-1)
      if(vmax-vmin .lt. .005*vmax .or. nvels .eq. 1)return
c   ***   fft, downward continue, and inverse fft for shift=vmax
c   ***   velocity non-constant in this interval
c   ***   apply multi-velocity pspi
c   ***   make vmax slightly  larger than max of 1./sqrt(vv)
c   ***   make vmin slightly smaller than min of 1./sqrt(vv)
        vmax=1.0001*vmax
        vmin=.9999*vmin
      nsteps=(vmax-vmin)/(.2*vmax) +1
      nsteps=nvels-1
      vvinc= (vmax-vmin)/nsteps
      shftinc=(shftmax-shftmin)/nsteps
      call movlev(y,b,2*nxy)
      call clear(y,2*nxy)
c   ***   start multiple velocity loop
        do k=1,nxy
        rs(k)=1./sqrt(vv(k))
        end do
      do istep=1,nsteps
      vv1=vmin+(istep-1)*vvinc
      vv2=vv1+vvinc
      v1 = (1./vv1)**2
      v2 = (1./vv2)**2
      shft1=shftmin+(istep-1)*shftinc
      shft2=shft1+shftinc
      rdv=1./(vv2-vv1)
c   ***   compute contribution of field downward continued using vv1
      do k=1,nxy
      rss = (vv2-rs(k))*rdv
      if(0. .lt. rss .and. rss .le. 1.)then
      xr=real(y(k))+rss*real(b(k))
      xi=aimag(y(k))+rss*aimag(b(k))
      y(k)=cmplx(xr,xi)
      end if
      end do
c   ***
      call movlev(a,b,2*nxy)
      call dwnn(b,nx,ny,dx,dy,w,shft2,v2,sgn)
c   ***   inverse fft data and scale
      call fftl2d(b,rs(nxy+1),nx,ny,-1)
c   ***   compute contribution of field downward continued using vv2
      do k=1,nxy
      rss = (rs(k)-vv1)*rdv
      if(0. .lt. rss .and. rss .lt. 1.)then
      xr=real(y(k))+rss*real(b(k))
      xi=aimag(y(k))+rss*aimag(b(k))
      y(k)=cmplx(xr,xi)
      end if
      end do
c   ***
      end do
c   ***
      return
      end
      subroutine phs2d(y,a,b,rs,nx,ny,dx,dy,w,shft,vv,sgn)
      complex y(*),a(*),b(*) 
      real shft(*),vv(*),rs(*)
c   ***   find min and max shifts 
      nxy=nx*ny
      v1=0.
      v2=1000000.
      shft1=1000000.
      shft2=-1000000.
      do 5 i=1,nxy
      v1=amax1(v1,vv(i))
      v2=amin1(v2,vv(i))
      shft1=amin1(shft1,shft(i))
      shft2=amax1(shft2,shft(i))
5     continue
      vmax=1./sqrt(v2)
      vmin=1./sqrt(v1)
c   ***   fft, downward continue, and inverse fft for vmin
c   ***   fft, downward continue, and inverse fft for shift=vmin
      call fftl2d(y,b,nx,ny,+1)
c   ***   copy field to scratch array 
      do 7 i=1,nxy
7     a(i)=y(i) 
      call dwnn(y,nx,ny,dx,dy,w,shft1,v1,sgn)
      call fftl2d(y,b,nx,ny,-1)
c   ***   fft, downward continue, and inverse fft for shift=vmax
      if(vmax-vmin.lt..02*vmax)return
      call dwnn(a,nx,ny,dx,dy,w,shft2,v2,sgn)
      call fftl2d(a,b,nx,ny,-1)
c   ***   interpolate for shifts between vmin and vmax
      r=1./(vmax-vmin)
      eps=1.e-8
      do 41 i=1,nxy
      tmp=(vmax-1./sqrt(vv(i)))*r 
      y(i)=tmp*y(i)+(1.-tmp)*a(i)
41    continue
      return
      end
      subroutine pspi(y,a,b,s,rs,nl,nh,nt1,nx,shft,vv,sgn) 
      complex y(1),a(1),b(1),s(1) 
      real shft(1),vv(1),rs(1)
      nws=nh-nl+1 
      dw=2.*3.1415927/nt1 
c     SEASEIS: dx nao era definido (dwnn usa dy/dx); igual a xmultpspi
      dx=1.
      f1=1./nx
      f2=-1.
c   ***   find min and max shifts 
      v1=0.
      v2=1000000.
      shft1=1000000.
      shft2=0.
      do 5 i=1,nx
      v1=amax1(v1,vv(i))
      v2=amin1(v2,vv(i))
      shft1=amin1(shft1,shft(i))
      shft2=amax1(shft2,shft(i))
5     continue
      vmax=1./sqrt(v2)
      vmin=1./sqrt(v1)
c   ***   fft, downward continue, and inverse fft for shift=vmin
c     call cffts(y,nx,1,nws,nx,+1,ier)
      call fftc(y,s,2,2*nx,nx,nws,-1)
c   ***   copy kx-w field to scratch array 
      do 7 i=1,nx*nws 
7     a(i)=y(i) 
c     call dwnv(y,y,b,rs,rs(nx+1),nx,nl,nh,dw,shft1,v1,sgn)
c   ***   downward continue by phase-shift for shft1,v1
      do iw=nl,nh
      w=(iw-1)*dw
      ii=(iw-nl)*nx+1
      call dwnn(y(ii),nx,1,dx,dx,w,shft1,v1,sgn)
      enddo
c     call cffts(y,nx,1,nws,nx,-1,ier)
c   ***   inverse fft data and scale
      call fftc(y,s,2,2*nx,nx,nws,+1)
      do j=1,nx*nws
      y(j)=f1*y(j)
      end do
      if(vmax-vmin .lt. .01*vmax)then
c     call vshift1(y,a,a,b,shft,dw,nl,nh,nx,sgn)
      call sgivshift1(y,shft,dw,nl,nh,nx,sgn)
      return
      end if
c   ***   fft, downward continue, and inverse fft for shift=vmax
c     call dwnv(a,a,b,rs,rs(nx+1),nx,nl,nh,dw,shft2,v2,sgn)
c   ***   downward continue by phase-shift for shft2,v2
      do iw=nl,nh
      w=(iw-1)*dw
      ii=(iw-nl)*nx+1
      call dwnn(a(ii),nx,1,dx,dx,w,shft2,v2,sgn)
      enddo
c   ***   inverse fft data and scale
      call fftc(a,s,2,2*nx,nx,nws,+1)
      do j=1,nx*nws
      a(j)=f1*a(j)
      end do
c   ***   interpolate for shifts between vmin and vmax
      eps=1.e-8
      r=1./(vmax-vmin+eps)
      do 35 i=1,nx
      rs(i)=(vmax-1./sqrt(vv(i)))*r
35    continue
      l=1
      do 40 j=1,nws 
      do 41 i=1,nx
      y(l)=rs(i)*y(l)+(1.-rs(i))*a(l)
c     y(l)=(rs(i)*cabs(y(l))+rs(nx+i)*cabs(a(l)) )*b(l)/cabs(b(l))
      l=l+1
41    continue
40    continue
c     call vshift1(y,a,a,b,shft,dw,nl,nh,nx,sgn)
      call sgivshift1(y,shft,dw,nl,nh,nx,sgn)
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
      subroutine sgivshift1(y,shft,dw,nl,nh,nx,sgn)
      real y(*)
      real shft(*)
c   ***   efficient  complex exponential table computation
      w1=sgn*(nl-1)*dw
      dws=sgn*dw
      nws=nh-nl+1
      nx2=2*nx
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
      return
      end
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
      subroutine splyn(x,xx,y,f,a,b,c,d,n,m) 
      dimension x(*),xx(*),y(*),f(*),a(*),b(*),c(*),d(*) 
c       cubic spline interpolation. input points are (x(i),y(i)),i=1,n.  
c     output points are (xx(i),f(i)),i=1,m. output independent variables 
c     xx(i) must safisfy x(1)@xx(i). if xx(i)>x(n)then f(i)=0. 
      n1=n-1 
      n2=n-2 
      n3=n-3 
      if(n1.gt.1)go to 10  
      den=1./(x(2)-x(1)) 
      do 9 i=1,m 
      cl=(x(2)-xx(i))*den  
      cr=(xx(i)-x(1))*den  
      f(i)=cl*y(1)+cr*y(2) 
9     continue 
      go to 15 
10    continue 
c       compute lhs and rhs of tridiagonal matrix system 
      do 1 i=1,n2  
      h1=x(i+1)-x(i) 
      h2=x(i+2)-x(i+1) 
      a(i)=h1*h1*h2  
      c(i)=h1*h2*h2  
      b(i)=2.*(a(i)+c(i))  
      d(i+1)=6.*(h1*(y(i+2)-y(i+1))-h2*(y(i+1)-y(i)))  
1     continue 
c        solve tridiagonal matrix system 
      c6=1./6. 
      b(1)=1./b(1) 
      c(1)=c(1)*b(1) 
      d(2)=d(2)*b(1) 
      if(n3.eq.1)go to 6 
      if(n3.eq.0)go to 5 
      do 2 i=2,n3  
      b(i)=1./(b(i)-c(i-1)*a(i-1)) 
      c(i)=b(i)*c(i) 
      d(i+1)=b(i)*(d(i+1)-a(i-1)*d(i)) 
2     continue 
6     continue 
      b(n2)=b(n2)-c(n3)*a(n3)  
      d(n2+1)=(d(n2+1)-a(n3)*d(n2))/b(n2)  
      do 3 k=1,n3  
      i=n2-k 
      d(i+1)=d(i+1)-c(i)*d(i+2)  
3     continue 
5     continue 
      d(1)=0.  
      d(n)=0.  
c       spline coeficients are in d  
c       computation of cubic spline polynomial coefficients  
      do 4 i=1,n1  
      hi=x(i+1)-x(i) 
      hib=1./hi  
      b(i)=.5*d(i) 
      a(i)=(d(i+1)-d(i))*c6*hib  
      c(i)=hib*(y(i+1)-y(i))-c6*(hi*(2.*d(i)+d(i+1)))  
4     continue 
      eps=(x(2)-x(1))*.0001  
      j=1  
      do 11 i=1,m  
13     continue  
      del1=xx(i)-x(j)  
      del2=x(j+1)-xx(i)  
      if(del1+eps.ge.0..and.del2+eps.ge.0.)go to 12  
      if(del1.lt.0..and.j.eq.1)go to 12  
      if(j.eq.n1)go to 12  
      j=j+1  
      go to 13 
12     continue  
      f(i)=y(j)+del1*(c(j)+del1*(b(j)+del1*a(j)))  
11     continue  
14     continue  
15    continue 
      return 
      end
      subroutine trans(b,a,nr,nc) 
      dimension a(*),b(*) 
      jj=1
      do 1 i=1,nc 
       ii=i 
      do 1 j=1,nr 
      b(ii)=a(jj) 
      ii=ii+nc
1     jj=jj+1 
      return
      end
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
*
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
*
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
*
      subroutine vshift(y,a,b,s,shft,dw,nl,nh,nx,sgn)
      complex y(*),s(*),a(*),b(*)
      real shft(*)
c   ***   efficient  complex exponential table computation
c   ***   a and b are equivalent.
      w1=sgn*(nl-1)*dw
      dws=sgn*dw
      nws=nh-nl+1
      do 10 i=1,nx
      s(i)=cexp(cmplx(0.,dws*shft(i)))
10    a(i)=cexp(cmplx(0.,w1*shft(i)))
      do 20 j=2,nws
      ii=(j-2)*nx
      iii=ii+nx
      do 20 i=1,nx
20    b(iii+i)=s(i)*a(ii+i)
c   ***   rearrange data into same order as y
      call ctran(b,s,nx,nws)
      do 30 i=1,nx*nws
30    y(i)=s(i)*y(i)
      return
      end
      subroutine vshift1(y,a,b,s,shft,dw,nl,nh,nx,sgn)
      complex y(*),s(*),a(*),b(*)
      real shft(*)
c   ***   efficient  complex exponential table computation
c   ***   a and b are equivalent.
      w1=sgn*(nl-1)*dw
      dws=sgn*dw
      nws=nh-nl+1
      if(nx.gt.1)then
      do 10 i=1,nx
      s(i)=cexp(cmplx(0.,dws*shft(i)))
10    a(i)=cexp(cmplx(0.,w1*shft(i)))
      do 20 j=1,nws
      ii=(j-1)*nx
      do 20 i=1,nx
      y(ii+i)=a(i)*y(ii+i)
      a(i)=s(i)*a(i)
20    continue
      else
      log2=alog(real(nws)-1.)/alog(2.)
      s(1)=cexp(cmplx(0.,dws*shft(1)))
      s(2)=cexp(cmplx(0.,w1*shft(1)))
      a(1)=cmplx(1.,0.)
      a(2)=s(1)
      i0=1
      do l=1,log2
      i0=2*i0
      s(1)=s(1)*s(1)
      i1=min0(i0,nws-i0)
      do i=1,i1
      b(i0+i)=s(1)*a(i)
      end do
      end do
      do i=1,nws
      y(i)=s(2)*b(i)*y(i)
      end do
      end if
      return
      end
      SUBROUTINE FFTFLZ(F,T,N,CS,NW,IP,LF,LFOLD,INC,JUMP,LOT)
      COMPLEX F(*),T(*),CS(*)
      COMPLEX COEF,COEF1,COEF2,COEF3
      IT=1
      LT0=1
      DO 400 JF=1,NW
         JT=IT
         NDX=LT0
         LT=LT0
         IF(IP.EQ.2)THEN
            COEF=CS(NDX)
            DO 210 K=JT,JT+LF-1
               I1 = (K+LT0-JT-1)*INC
               I2 = (K-1)*INC
               I3 = (K+LF-1)*INC
               DO 211 I=1,LOT
                  T(I1+I)=F(I2+I)+COEF*F(I3+I)
  211          CONTINUE
  210       CONTINUE
         ELSEIF(IP.EQ.3)THEN
            COEF=CS(NDX)
            COEF1=COEF*COEF
               DO 220 K=JT,JT+LF-1
               I1 = (K+LT0-JT-1)*INC
               I2 = (K-1)*INC
               I3 = (K+LF-1)*INC
               I4 = (K+LF+LF-1)*INC
                  DO 221 I=1,LOT
                     T(I1+I)=F(I2+I)+COEF*F(I3+I)+COEF1*F(I4+I)
  221             CONTINUE
  220          CONTINUE
         ELSEIF(IP.EQ.4)THEN
            COEF=CS(NDX)
            COEF1=COEF*COEF
            COEF2=COEF*COEF1
            DO 230 K=JT,JT+LF-1
               I1 = (K+LT0-JT-1)*INC
               I2 = (K-1)*INC
               I3 = (K+LF-1)*INC
               I4 = (K+LF+LF-1)*INC
               I5 = (K+LF+LF+LF-1)*INC
               DO 231 I=1,LOT
                  T(I1+I)=F(I2+I)+COEF*F(I3+I)+COEF1*F(I4+I)+
     :                    COEF2*F(I5+I)
  231          CONTINUE
  230       CONTINUE
         ELSEIF(IP.EQ.5)THEN
            COEF=CS(NDX)
            COEF1=COEF*COEF
            COEF2=COEF*COEF1
            COEF3=COEF*COEF2
            DO 530 K=JT,JT+LF-1
               I1 = (K+LT0-JT-1)*INC
               I2 = (K-1)*INC
               I3 = (K+LF-1)*INC
               I4 = (K+LF+LF-1)*INC
               I5 = (K+LF+LF+LF-1)*INC
               I6 = (K+LF+LF+LF+LF-1)*INC
               DO 531 I=1,LOT
                  T(I1+I)=F(I2+I)+COEF*F(I3+I)+COEF1*F(I4+I)+
     :                    COEF2*F(I5+I)+COEF3*F(I6+I)
  531          CONTINUE
  530       CONTINUE
         ELSE
            DO 320 K=JT,JT+LF-1
               DO 321 I=1,LOT
                  T(1+(K+LT0-JT-1)*INC+(I-1)*JUMP) =
     :               F(1+(K-1)*INC+(I-1)*JUMP)
  321          CONTINUE
  320       CONTINUE
            JT=JT+LF
            DO 350 IFE=2,IP
               COEF=CS(NDX)
               LT=LT0
               DO 345 K=JT,JT+LF-1
                  INDX=1+(K+LT0-JT-1)*INC
                  DO 346 I=1,LOT
C                     T(1+(K+LT0-JT-1)*INC+(I-1)*JUMP) =
C     :                  T(1+(K+LT0-JT-1)*INC+(I-1)*JUMP) +
C     :                  COEF*F(1+(K-1)*INC+(I-1)*JUMP)
                     T(INDX+(I-1)*JUMP) =
     :                  T(INDX+(I-1)*JUMP) +
     :                  COEF*F(1+(K-1)*INC+(I-1)*JUMP)
  346             CONTINUE
  345          CONTINUE
               JT=JT+LF
               NDX=NDX+LT0-1
               IF(NDX.GT.N)NDX=NDX-N
  350       CONTINUE
         ENDIF
         IT=IT+LFOLD
         IF(IT.GT.N)IT=1
         LT0=LT0+LF
  400 CONTINUE
      RETURN
      END
*
