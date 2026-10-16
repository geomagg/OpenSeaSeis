c   ***   SEISPAK xscale: subroutine XINT2, unchanged.
c   ***   The main program (jhh I/O, values read from standard input) is replaced by mod_spk_xscale.cc
      subroutine xint2(f,x,y,x1,dx,n,m) 
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
