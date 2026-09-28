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
