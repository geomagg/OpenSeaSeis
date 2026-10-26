c   ***   SEISPAK tau2vfix: rotinas originais xtau2v, factm, factp
c   ***   (fonte tau2vfix.f; programa principal em mod_spk_tau2vfix.cc).
c   ***   Mudancas marcadas SEASEIS.
      subroutine xtau2v(tdata,vdata,vel,vmin,vmax,nt,dt,jnt)
c  ***  tdata     real   input (nt x jnt) tau deviation data
c  ***  vdata     real   output (nt x jnt) estimated velocity
c  ***  vel       real   input rms velocity sampled at (j-1)*dt
c  ***  vmin      real   minimum velocity to map
c  ***  vmax      real   maximum velocity to map
c  ***  nt        int    number of time samples
c  ***  dt        real   time sample interval in seconds
c  ***  jnt       int    number of tau deviation data traces
      real vdata(nt,jnt), tdata(nt,jnt),vel(nt)
c   ***   map tau-deviation data (dtau) to corrected rms velocity (vc)
c   ***   Corrected rms velocity is estimated from tau deviation using:
c   ***   VC = VM * (1. - DTAU/TAU)
c   ***   DTAU = TAU * (1. - (VC/VM)**2 )
      factor=vdata(1,1)
      ffactor=factor
      sr=1000.*dt
      dv=(vmax-vmin)/(jnt-1)
c   ***   dtau=0 at tau-deviaton trace  rjnt = jnt/2+1
      rjnt=.5*(jnt-1.)+1.
      jnt2=(jnt-1)/2+1
      do n=1,nt
      do j=1,jnt
      vm=vel(n)
      vc=vmin+(j-1)*dv
      tau=(n-1)*dt
      tt=1000.*tau
       if(ffactor.lt..1)then
       factor=factm(tt)
       else
       factor=ffactor
       endif
       dtau=factor*tau*(1.-(vc/vm))
c   ***  set rj to trace # for this dtau
      rj=dtau/dt+rjnt
      jj=int(rj)
      if(1 .le. jj .and. jj .lt. jnt)then
      r1=rj-jj
      r2=1.-r1
      vdata(n,j)=r1*tdata(n,jj+1)+r2*tdata(n,jj)
      else
      vdata(n,j)=0.
      end if
      end do
c   ***   set dtau=0 location to 0
      ii=(vel(n)-vmin)/dv+1.5
c   ***   SEASEIS: so dentro do painel (o original escrevia fora se vel < vmin ou > vmax)
      if(ii.ge.1 .and. ii.le.jnt) vdata(n,ii)=0.
      end do
      end
      function factm(ttx)
c  calibration factors for xtau2v
c  factors are for times 0 - 10000 mls at a 500 mls interval
c  input ttx is time in mls at which factor returned in fact is applied
c  factors are linearly interpolated between values given if f
      real f(21)
      data f/1., .9, .8, .7, .65, .6, .6, .6, .6, .55, .55, .5, .5, .5,
     : .5, .5, .5, .5, .5, .5, .5/
      tt=.001*amax1(0.,amin1(ttx,10000.))
      tt=.5*tt+1.
      it=tt
      x=tt-it
      factm=x*f(it+1)+(1.-x)*f(it)
      return
      end
      function factp(ttx)
c  calibration factors for xtau2v
c  factors are for times 0 - 10000 mls at a 500 mls interval
c  input ttx is time in mls at which factor returned in fact is applied
c  factors are linearly interpolated between values given if f
      real f(21)
      data f/2., .7, .7, .55, .55, .6, .6, .5, .5, .5, .5, .5, .5, .5,
     : .5, .5, .5, .5, .5, .5, .5/
      tt=.001*amax1(0.,amin1(ttx,10000.))
      tt=tt+1.
      it=tt
      x=tt-it
      factp=x*f(it+1)+(1.-x)*f(it)
      return
      end
