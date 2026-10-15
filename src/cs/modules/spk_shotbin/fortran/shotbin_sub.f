c   ***   SEISPAK xxshotbin: subroutines SPL and INT3, unchanged.
c   ***   The main program BIN (jhh I/O, TABLE1) is replaced by mod_spk_shotbin.cc
      subroutine spl(y,x,f,s,n,m)  
c   ***   resample n equally spaced samples  of y having unit increment  
c   ***   at (x(i),i=1,m) and return results in (f(i),i=1,m).  
      real x(1),f(1),y(1),s(1) 
      data b1/.25/,b2/.2667/,b3/.2679/ 
c   compute spline coefficients  
      hy=1. 
      r=6./(hy*hy) 
      s(1)=0.  
      s(2)=b1*r*(y(3)-2*y(2)+y(1)) 
      s(3)=b2*(r*(y(4)-2*y(3)+y(2))-s(2))  
      do 10 i=4,n-1  
10    s(i)=b3*(r*(y(i+1)-2*y(i)+y(i-1)) -s(i-1)) 
      s(n)=0.  
      do 20 i=n-2,4,-1 
20    s(i)=s(i)-b3*s(i+1)  
      s(3)=s(3)-b2*s(4)  
      s(2)=s(2)-b1*s(3)  
c   compute polynomial coefficients and interpolate  
      do 21 i=1,m  
21    f(i)=0.  
      r1=1./hy 
      r2=1./(6.*hy)  
      r3=hy/6. 
      j=1  
      xy=1.  
      a=s(2)*r2  
      b=0. 
      c=(y(2)-y(1))*r1-s(2)*r3 
      do 30 i=1,m  
13    d1=x(i)-xy 
      d2=hy-d1 
      if(d1.ge.0..and.d2.ge.0.)goto 12 
      j=j+1  
      xy=xy+hy 
      if(j.eq.n)goto 40  
      a=(s(j+1)-s(j))*r2 
      b=.5*s(j)  
      c=(y(j+1)-y(j))*r1-(2.*s(j)+s(j+1))*r3 
      goto 13  
12    continue 
      f(i)=y(j)+d1*(c+d1*(b+d1*a)) 
30    continue 
40    return 
      end  
      subroutine int3(f,x,y,x1,dx,n,m) 
c   interpolate m points (x(i),y(i))  for x = x1 to x1+(n-1)*dx  
c   x1 <= x(1) < x(i) < x(i+1) must be satisfied.  
c   f(x)=0 if x < x1 or x > x1+(n-1)*dx  
c   interpolated results in f  
      double precision d1,d2,a,r1  
      dimension f(*),x(*),y(*) 
      j=1  
      xf=x1-dx 
      a=0. 
      eps=.0001  
      if(x(2)-x(1).gt.eps)a=(y(2)-y(1))/(x(2)-x(1))  
      do 30 i=1,n  
      xf=xf+dx 
13    d1=xf-x(j) 
      d2=x(j+1)-xf 
      if(d1.lt.0.)then 
      f(i)=0.  
      goto 29  
      end if 
      if(d2.lt.0.)then 
      j=j+1  
      if(j.eq.m)goto 40  
      r1=x(j+1)-x(j) 
      a=y(j+1)-y(j)  
      if(r1.gt.eps)then  
      a=a/r1 
      else 
      a=0. 
      end if 
      goto 13  
      else 
      e=d1*a 
      f(i)=y(j)+e  
      end if 
29    continue 
30    continue 
40    do 50 j=i,n  
50    f(j)=0.  
      return 
      end  
