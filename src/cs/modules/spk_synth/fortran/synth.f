c-----------------------------------------------------------------------
c     Processo SYNTH do SEISPAK (dados sinteticos), para o modulo SPK_SYNTH.
c
c     SPKSYNTH substitui SUB2RUN + SYNTH do original:
c       - parametros recebidos como argumentos (antes: tabela TABLE1
c         escrita por ESYNTH e lida com GH/GC)
c       - o conjunto de dados de saida (antes: arquivo ONLINE1 criado com
c         PGLOB/TGPINL) e mantido na memoria pelo modulo C++, que fornece
c         GETTRC/PUTTRC(IR,IT,S)
c     As demais rotinas (RDVEL, GAUS, DMOSUB, BN3D, POINTS, LINES, HYPBS,
c     HYPBS3D) sao as originais. Diferencas marcadas com SEASEIS:.
c
c     Cabecalho de cada traco (LHDR=7), como no ONLINE1 original:
c       S(1) RNUM  S(2) TNUM  S(3) RANG  S(4) SX  S(5) GX  S(6) ILIN  S(7) XLIN
c-----------------------------------------------------------------------
      SUBROUTINE SPKSYNTH(S,VFILE,GSSN,GGD,SPJ,SPI,ZRNGTRC,TPWR,LRM,SMP,
     :                    IHMR,IHMT,MPOINT,POINT,MLINE,RLINE,
     :                    MHYPB,HYPB,MHYPB3D,HYPB3D,MBIN3D,BIN3D,
     :                    NREC,IERR)
c     IERR  0 ok
c           1 BIN3D: mais registros do que ranges
c           2 BIN3D: range nao encontrado
c     NREC  numero de registros gerados (HMR, ou o lido de VFILE)
      CHARACTER VFILE*240
      DIMENSION S(*)
      DIMENSION POINT(*), RLINE(*), HYPB(*), HYPB3D(*), BIN3D(*)
      COMMON/SPKSYNE/IERRS
      IERR=0
      IERRS=0
      NREC=IHMR
c     SEASEIS: os vetores ja chegam com o zero final (POINT(MPOINT+1)=0 etc.)
      izrt=zrngtrc+0.4999
      write(*,2001)
      write(*,2002)4, ggd, izrt, lrm, smp, ihmr, ihmt
 2002 format(15x,i5,2x,F6.1,2x,i4,2x,i5,2x,f5.1,2x,I5,2x,I5)
 2001 format(15x,'nbyte',5x,'ggd',2x,'izrt',2x,'  lrm'
     :      ,2x,'  smp',2x,' ihmr',2x,' ihmt')
c                                       indices dos cabecalhos (PGLOB)
      RANG=3.0
      LHDR=7
      NT=1+(LRM/SMP)
      IF(VFILE.EQ.'NONE')THEN
         NNN=1
         DO 90 IR=1,IHMR
           if(mbin3d.ne.0)then
              RNG=bin3d(NNN+1)
              NNN=NNN+bin3d(NNN)+1
              if(nnn.gt.mbin3d .and. IR.LT.IHMR) then
                write(*,*)'More records than ranges.'
c     SEASEIS: erro devolvido ao modulo em vez de STOP
                IERR=1
                RETURN
              endif
           endif
           DO 80 IT=1,IHMT
             do 70 i=1,lhdr+nt
                s(i)=0.0
   70        continue
             s(1)=IR
             s(2)=IT
             if(mbin3d.eq.0)then
                s(3)=(IT-ZRNGTRC)*GGD
                s(4)=SPI+(IR-1)*SPJ*GGD
                s(5)=s(4)+s(3)
                s(6)=0.
                s(7)=0.
             else
                s(3)=RNG
                s(4)=0.0
                s(5)=(IT-1)*GGD
                s(6)=0.
                s(7)=0.
             endif
             call puttrc(IR,IT,S)
   80      CONTINUE
   90    CONTINUE
c
c                                               call the POINTS subroutine
         CALL POINTS(s,point,IHMR,IHMT,NT,LHDR,SMP,GSSN)
c                                               call the LINES subroutine
         CALL LINES(s,rline,IHMR,IHMT,NT,LHDR,SMP,GSSN)
c                                               call the HYPBS subroutine
         CALL HYPBS(s,hypb,IHMR,IHMT,NT,LHDR,SMP,GSSN,TPWR,GGD)
c                                               call the HYPBS3d subroutine
         CALL HYPBS3D(s,hypb3d,IHMR,IHMT,NT,LHDR,SMP,GSSN,TPWR,GGD)
c                                               call the BIN3D subroutine
         CALL BN3D(s,bin3d,IHMR,IHMT,NT,LHDR,SMP,GSSN,RANG,GGD)
         IERR=IERRS
c
      ELSE
         CALL RDVEL(VFILE,S,NREC,LHDR+NT,LHDR,SMP)
      ENDIF
      RETURN
      END
      SUBROUTINE RDVEL(VFILE,S,NREC,NS,LHDR,SMP)
      character vfile*240
      dimension S(NS)
      NT=NS-LHDR
      IRNG=3
      ISX =4
      ISY =5
      ILI =6
      ILX =7
      luin=3
      open(luin,FILE=VFILE,ERR=45)
      LIo=-1000000000
      LXo=-1000000000
      IT=1
      IR=0
   10 continue
      read(luin,*,end=50, err=40) LI, LX, X, Y, TM, VL
c...  write(*,*)'LI,LX,VL,TM,IR:',LI,LX,VL,TM,IR
      IF(LI.NE.LIo .OR. LX.NE.LXo)THEN
         LIo=LI
         LXo=LX
         IF(IR.NE.0) call puttrc(IR,IT,S)
         IR=IR+1
         S(1)=IR
         S(2)=IT
         S(IRNG)=0.
         S(ISX )=X
         S(ISY )=Y
         S(ILI )=LI
         S(ILX )=LX
         N=INT(Tm/SMP)
         DO 20 I=1,N
            S(LHDR+I)=VL
   20    CONTINUE
         VLo=VL
         TMo=TM
c     SEASEIS: apos o ultimo ponto (ou com um so ponto) vale a ultima
c              velocidade lida (antes: V do traco ou trecho anterior)
         V=VL
      ELSE
         M=nint(Tm/SMP)
         ODTM=1.0/(TM-TMo)
         DO 30 I=N+1,M
            IF(I.GT.NT) GO TO 30
c     SEASEIS: tempo da amostra I (antes: T=(N-1)*SMP, o inicio do trecho,
c              o que deixava cada trecho constante em vez de interpolado)
            T=(I-1)*SMP
            F=ODTM*(T-TMo)
            V=F*VL+(1.0-F)*VLo
            S(LHDR+I)=V
   30    CONTINUE
         VLo=VL
         TMo=TM
         N=M
c     SEASEIS: extrapolacao constante com a ultima velocidade lida
         V=VL
      ENDIF
      DO 35 I=N+1,NT
         S(LHDR+I)=V
   35 CONTINUE
      GO TO 10
   40 CONTINUE
      WRITE(*,*)'==================================================='
      WRITE(*,*)'Error reading velocity file. Last rec read ecord IR'
      WRITE(*,*)'==================================================='
      GO TO 50
   45 CONTINUE
      WRITE(*,*)'==================================================='
      WRITE(*,*)'Error opening velocity file.'
      WRITE(*,'(''VFILE:'',A240)') VFILE
      WRITE(*,*)'==================================================='
   50 CONTINUE
      IF(IR.NE.0) call puttrc(IR,IT,S)
      NREC=IR
      return
      end
      subroutine gaus(s,lhdr,nt,smp,tm,amp0,gssn)
      dimension s(1)
      data const /0.693147/
      w=0.5*gssn/smp
      a=const/(w*w)
      iw=2.577*w + 0.49999
      itm= tm + 0.49999
      do 40 i=itm-iw,itm+iw
       if(i.lt.1 .or. i.gt.nt)go to 40
       x=i-tm
       amp=amp0 * exp(-a*x*x)
       s(lhdr+i)=s(lhdr+i) + amp
   40 continue
      return
      end
      subroutine dmosub(v,sx,sy,gx,gy,bx,by,br0,dp0,th,ph
     :                 , t0,t,tn,tr,x,sind,dd,ifl)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  INPUT
c	sx	shot X coordinate
c	sy	shot Y coordinate
c	gx	geophone X coordinate
c       gy      geophone Y coordindate
c	bx	bin center X coordinate
c	by	bin center Y coordinate
c       br0     bin radius
c	dp0	perpendicular distance from bin to plane reflector (> 0).
c	th      angle between z axis and reflector unit normal
c       ph      angle between x axis and projection of unit normal on X-Y plane
c  OUTPUT
c       t0      zero range time at midpoint
c	t       arrival time
c	tn      arrival time after NMO with RMS (here constant) velocity
c	tr	dmo impulse response time at bin center
c       x       distance from bin location to midpoint
c    sind       sin of dip angle
c      dd       perpendicular distance from bin to S-G line
c     ifl       0 if valid 
c               1 if dd > br
c		2 if xoh> 1
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      d2r=3.14159265/180.0
      ntime=1
      br=abs(br0)
   10 continue
c             				find distance from bin to S-G line
      den=sx-gx
      if(abs(den) .gt. 0.1)then
         sl = (sy-gy)/(sx-gx)
         osl=-1.0/sl
         b  = sy-sl*sx
         bp = by+bx*osl
         if(abs(sl).gt.0.00001)then
            xx = (bp-b)/(sl-osl)
            yy = sl*xx+b
         else
            yy=0.5*(sy+gy)
            xx=bx
         endif
      else
         sl=1.0E12
         xx= 0.5*(sx+gx)
         yy= by
      endif
      xbx=xx-bx
      yby=yy-by
      dd = sqrt(xbx*xbx+yby*yby)
 3001 format(1x,F10.3,2x,F10.3)
c					compute unit vector from bin
c					location to shot-geophone line
      if(dd.gt.br)then
         ifl=1
         return
      elseif(dd.gt.0.0)then
         uxx=xbx/dd
         uyy=yby/dd
      else
         uxx=0.0
         uyy=0.0
      endif
c					compute unit normal to surface
      unz=cos(d2r*th)
      uxy=sin(d2r*th)
      unx=uxy*cos(d2r*ph)
      uny=uxy*sin(d2r*ph)
c					compute range and unit vector
c					from geophone to shot
      xgs=sx-gx
      ygs=sy-gy
      rng=sqrt(xgs*xgs+ygs*ygs)
      h= 0.5 * rng
      umx=xgs/rng
      umy=ygs/rng
      umz=0.0
c					compute sin of dip angle which is
c					the cos of the angle between
c					um and un.
      sind=umx*unx+umy*uny
c					if br0 < 0 then set for constructive
c					interference at bin location
      if(br0.lt.0.0 .and. ntime.eq.1)then
	 ntime=2
         symby=sy-by
         sxmbx=sx-bx
         if(sxmbx.gt.0.000001)then
            sa=atan(symby/sxmbx)
         elseif(sxmbx.lt.-0.000001)then 
            sa=atan(symby/sxmbx)+3.14159265
         elseif(symby.gt.0.0)then
            sa=3.14159265/2.
         else
            sa=-3.14159265/2.
         endif   
         sd=(sx-bx + sy-by)/(cos(sa)+sin(sa))
         b=abs(2.0*h*h*sind/(v*tb+4.0*h*h*sind*sind/(v*tb)))
         sd=h-b
         SX=BX+SD*COS(SA)
         SY=BY+SD*SIN(SA)
         GX=SX-RNG*COS(SA)
         GY=SY-RNG*SIN(SA)
	 write(*,*)'New SD and SA:',SD,SA
	 GO TO 10
      endif
c					compute "u" geophone to projected
c					bin distance
      xgb=xx-gx
      ygb=yy-gy
      u=  sqrt(xgb*xgb+ygb*ygb)
c					compute "sd" shot to projected bin
c					distance
      xsb=xx-sx
      ysb=yy-sy
      sd =sqrt(xsb*xsb+ysb*ysb)
c					find midpoint coordinates
      dx =0.5*(sx+gx)
      dy =0.5*(sy+gy)
c					compute "x" midpoint to bin
c					distance
      xbd=xx-dx
      ybd=yy-dy
      x=  sqrt(xbd*xbd+ybd*ybd)
c					compute the cos of the angle between
c					surface unit normal and the line
c					from the bin to S-G line
      cosb=uxx*unx+uyy*uny
c					correct the given depth for bin
c					offset from S-G line
      dp=dp0 + cosb * dd
c					compute midpoint to plane distance
      d0=dp + (h-u)*sind
c					compute bin to plane time
      tb =2.0*dp/v
c					compute midpoint to plane time
      t0 =2.0*d0/v
      if((tb.gt.t0) .or. (x.gt.2.0*h*h/(v*t0)) )x=-x
      vsq=v*v
c					find the vertical RMS velocity
c					NMO corrected time squared
      tn2=(t0*t0-4.0*h*h*sind*sind/vsq)
      tn =sqrt(tn2)
c					find the arrival time associated
c					with the midpoint
      t  =sqrt(tn2+4.0*h*h/vsq)
      xoh=x/h
      if(abs(xoh).gt.1.0)then
         ifl=2
      else
c					compute the dmo corrected time
         ifl=0
         tr=tn*sqrt(1.0-xoh*xoh)
      endif
      return
      end
      subroutine bn3d(s,bin3d,mr,mt,nt,lhdr,smp,gssn,rang,dx)
      dimension bin3d(1), s(1)
c     SEASEIS: codigo de erro para o modulo
      COMMON/SPKSYNE/IERRS
      d2r=3.14159265/180.0
      n=bin3d(1)+0.1
      if(n.le.0) return
      k=1
      ir=1
      irng=rang+0.1
      call gettrc(ir,1,s)
      r0=s(irng)
      write(*,2001)
   10 continue
      n=bin3d(k) + 0.1
      kk=k+n+1
      if(n.le.0) return
      gr=bin3d(k+1)
      if(abs(gr-r0).gt.1.0)then
   15    continue
         ir=ir+1
         if(ir.gt.mr) then
           write(*,*)'Range not found in SYNTH, routine bin3d.'
           write(*,*)'Looking for range:',gr
           write(*,*)'Last range found was:',r0
c     SEASEIS: erro devolvido ao modulo em vez de STOP
           IERRS=2
           return
         endif
         call gettrc(ir,1,s)
         r0=s(irng)
         if(abs(gr-r0).gt.1.0)go to 15
      endif
      k=k+2
   20 continue
      itr=bin3d(k)
      vel=bin3d(k+1)
      dp =bin3d(k+2)
      th =bin3d(k+3)
      ph =bin3d(k+4)
      sd =bin3d(k+5)
      sa =bin3d(k+6)
      bx=0.
      by=0.
      k=k+7
      BR=1.0
      SX=BX+abs(SD)*COS(D2R*SA)
      SY=BY+abs(SD)*SIN(D2R*SA)
      GX=SX-GR*COS(D2R*SA)
      GY=SY-GR*SIN(D2R*SA)
      if(sd.lt.0.0) br=-br
      call dmosub(vel,sx,sy,gx,gy,bx,by,br,dp,th,ph
     :           , t0,t,tn,tr,x,sind,dd,ifl)
      if(ifl.ne.0)go to 30
      it0m=1000.0*t0+0.499
      tm  =1000.0* t+0.499
      itm =tm
      itnm=1000.0*tn+0.499
      itrm=1000.0*tr+0.499
      ix  = x +0.499
      write(*,2002) vel,sd,sa,gr,dp,th,ph,it0m,itm,itnm,itrm,ix,sind
 2002 format(2(F6.0,1x),F4.0,2(1x,F6.0),2(1x,F4.0),4(1x,I5)
     :       ,1x,I6,1x,F6.4)
 2001 format(3x,'VEL',5x,'SD',3x,'SA',5x,'GR',5x,'DP',3x,'TH',3x
     :      ,'PH',4X,'T0',5X,'T',4X,'TN',4X,'TR',6X,'X',2x,'SIN d')
      it0=itr+x/dx +0.49
      if(it0.lt.1 .or. it0.gt.mt) go to 30
      trc=itr+x/dx
      pi=3.14159265
      sm=1.0+tm/smp
      amp0=100.0
      DO 25 i=0,10
         it=it0+i
         xx=it-trc
         if(abs(xx).gt.0.001)then
            amp=amp0*sin(xx*pi)/(xx*pi)
         else
            amp=amp0
         endif
         if(it.gt.mt)go to 25
         call gettrc(ir,it,s)
         call gaus(s,lhdr,nt,smp,sm,amp,gssn)
         call puttrc(ir,it,s)
         if(i.eq.0)go to 25
         it=it0-i
         xx=it-trc
         if(abs(xx).gt.0.001)then
            amp=amp0*sin(xx*pi)/(xx*pi)
         else
            amp=amp0
         endif
         if(it.lt.1)go to 25
         call gettrc(ir,it,s)
         call gaus(s,lhdr,nt,smp,sm,amp,gssn)
         call puttrc(ir,it,s)
   25 continue
   30 continue
      if(k.lt.kk) go to 20
      go to 10
      end
      subroutine points(s,pnt,mr,mt,nt,lhdr,smp,gssn)
      dimension pnt(1), s(1)
      k=1
   10 continue
      n=pnt(k) + 0.1
      k=k+1
      if(n.le.0)return
      ir1=pnt(k)
      it1=pnt(k+1)
      tm1=pnt(k+2)/smp
      a1=pnt(k+3)
      k=k+4
      n=n-4
      ibias=0
      istep=1
   15 continue
      if(n.eq.0)then
       ir2=ir1
       it2=it1
       tm2=tm1
       a2=a1
      else
       ir2=pnt(k)
       it2=pnt(k+1)
       tm2=pnt(k+2)/smp
       a2=pnt(k+3)
       k=k+4
       n=n-4
      endif
      if(ir1.gt.mr)ir1=mr
      if(ir2.gt.mr)ir2=mr
      if(ir1.le.ir2)then
       isgn=1
      else
       isgn=-1
      endif
      do 90 ir=ir1+isgn*ibias,ir2,isgn*istep
       if(ir1.eq.ir2)then
	F=0.0
       else
	D=ir2-ir1
	F=(ir-ir1)/D
       endif
       it=(1.0-F)*it1 + F*it2  +0.49999
       if(it.lt.1 .or. it.gt.mt)go to 90
       amp=(1.0-F)*a1 + f*a2
       tm=1+(1.0-F)*tm1 + F*tm2
       itm= tm + 0.49999
       if((itm.lt.1 .or. itm.gt.nt) .and. gssn.le.0.0)go to 90
       call gettrc(ir,it,s)
       if(gssn.le.0.0)then
         s(lhdr+itm)=s(lhdr+itm) + amp
       else
         call gaus(s,lhdr,nt,smp,tm,amp,gssn)
       endif
       call puttrc(ir,it,s)
   90 continue
      ir1=ir2
      it1=it2
      tm1=tm2
      a1=a2
      ibias=1
      if(n.gt.0)go to 15
      go to 10
      end
      subroutine lines(s,rlne,mr,mt,nt,lhdr,smp,gssn)
      dimension rlne(1), s(1)
      k=1
   10 continue
      n=rlne(k) + 0.1
      k=k+1
      if(n.le.0)return
      ir1=rlne(k)
      it1a=rlne(k+1)
      tm1a=1.0+rlne(k+2)/smp
      a1a=rlne(k+3)
      it1b=rlne(k+4)
      tm1b=1.0+rlne(k+5)/smp
      a1b=rlne(k+6)
      k=k+7
      n=n-7
      ibias=0
      istep=1
   15 continue
      if(n.eq.0)then
       ir2=ir1
       it2a=it1a
       tm2a=tm1a
       a2a=a1a
       it2b=it1b
       tm2b=tm1b
       a2b=a1b
      else
       ir2=rlne(k)
       it2a=rlne(k+1)
       tm2a=1.0+rlne(k+2)/smp
       a2a=rlne(k+3)
       it2b=rlne(k+4)
       tm2b=1.0+rlne(k+5)/smp
       a2b=rlne(k+6)
       k=k+7
       n=n-7
      endif
      if(ir1.le.ir2)then
       isgn=1
      else
       isgn=-1
      endif
      do 90 ir=ir1+isgn*ibias,ir2,isgn*istep
       if(ir.lt.1 .or. ir.gt.mr)go to 90
       if(ir1.eq.ir2)then
	F=0.0
       else
	D=ir2-ir1
	F=(ir-ir1)/D
       endif
       ita=(1.0-F)*it1a + F*it2a  +0.49999
       itb=(1.0-F)*it1b + F*it2b  +0.49999
       aa = (1.0-F)*a1a + F*a2a
       ab = (1.0-F)*a1b + F*a2b
       tma=(1.0-F)*tm1a + F*tm2a
       tmb=(1.0-F)*tm1b + F*tm2b
       if(ita.le.itb)then
	istp=1
       else
	istp=-1
       endif
       do 60 it=ita,itb,istp
	  if(it.lt.1 .or. it.gt.mt)go to 60
	  if(itb.eq.ita)then
	   F=0.0
	  else
	   D=itb-ita
	   F=(it-ita)/D
	  endif
	  tm=(1.0-F)*tma + F*tmb
	  itm= tm + 0.49999
	  amp=(1.0-F)*aa + F*ab
	  if((itm.lt.1 .or. itm.gt.nt) .and. gssn.le.0.0) go to 60
          call gettrc(ir,it,s)
	  if(gssn.le.0.0)then
            s(lhdr+itm)=s(lhdr+itm) + amp
	  else
	    call gaus(s,lhdr,nt,smp,tm,amp,gssn)
	  endif
          call puttrc(ir,it,s)
   60  continue
   90 continue
      ir1=ir2
      it1a=it2a
      it1b=it2b
      tm1a=tm2a
      tm1b=tm2b
      a1a=a2a
      a1b=a2b
      ibias=1
      if(n.gt.0)go to 15
      go to 10
      end
      subroutine hypbs(s,hyp,mr,mt,nt,lhdr,smp,gssn,tpwr,ggd)
      dimension hyp(1), s(1)
      k=1
   10 continue
      n=hyp(k) + 0.1
      k=k+1
      if(n.le.0)return
      ir1=hyp(k)
      it1=hyp(k+1)
      tm1=0.001 * hyp(k+2)
      a1=hyp(k+3)
      v1=hyp(k+4)
      k=k+5
      n=n-5
      ibias=0
      istep=1
   15 continue
      if(n.eq.0)then
       ir2=ir1
       it2=it1
       tm2=tm1
       a2=a1
       v2=v1
      else
       ir2=hyp(k)
       it2=hyp(k+1)
       tm2=0.001 * hyp(k+2)
       a2=hyp(k+3)
       v2=hyp(k+4)
       k=k+5
       n=n-5
      endif
      if(ir1.le.ir2)then
       isgn=1
      else
       isgn=-1
      endif
      do 90 ir=ir1+isgn*ibias,ir2,isgn*istep
       if(ir.lt.1 .or. ir.gt.mr) go to 90
       if(ir1.eq.ir2)then
	F=0.0
       else
	D=ir2-ir1
	F=(ir-ir1)/D
       endif
       it0=(1.0-F)*it1 + F*it2  +0.49999
       amp0=(1.0-F)*a1 + F*a2
       v=(1.0-F)*v1 + F*v2
       tm0=(1.0-F)*tm1 + F*tm2
       do 60 it=1,mt
	  x=ggd*(it-it0)
	  t=sqrt(tm0*tm0+x*x/(v*v))
	  tm=1.0+1000.0*t/smp
	  itm=tm + 0.49999
	  ratio=(tm0/t)**tpwr
	  if((itm.lt.1 .or. itm.gt.nt) .and. gssn.le.0.0) go to 60
          call gettrc(ir,it,s)
	  if(gssn.le.0.0)then
            s(lhdr+itm)=s(lhdr+itm) + ratio * amp0
	  else
	    call gaus(s,lhdr,nt,smp,tm,ratio*amp0,gssn)
	  endif
          call puttrc(ir,it,s)
   60  continue
   90 continue
      ir1=ir2
      it1=it2
      tm1=tm2
      a1=a2
      v1=v2
      ibias=1
      if(n.gt.0)go to 15
      go to 10
      end
      subroutine hypbs3d(s,hyp3d,mr,mt,nt,lhdr,smp,gssn,tpwr,ggd)
      dimension hyp3d(1), s(1)
      k=1
   10 continue
      n=hyp3d(k) + 0.1
      k=k+1
      if(n.le.0)return
      ir1=hyp3d(k)
      it1=hyp3d(k+1)
      tm1=0.001 * hyp3d(k+2)
      a1=hyp3d(k+3)
      v1=hyp3d(k+4)
      k=k+5
      n=n-5
      ibias=0
      istep=1
   15 continue
      if(n.eq.0)then
       ir2=ir1
       it2=it1
       tm2=tm1
       a2=a1
       v2=v1
      else
       ir2=hyp3d(k)
       it2=hyp3d(k+1)
       tm2=0.001 * hyp3d(k+2)
       a2=hyp3d(k+3)
       v2=hyp3d(k+4)
       k=k+5
       n=n-5
      endif
      if(ir1.le.ir2)then
       isgn=1
      else
       isgn=-1
      endif
      do 90 ir=ir1+isgn*ibias,ir2,isgn*istep
       if(ir.lt.1 .or. ir.gt.mr) go to 90
       if(ir1.eq.ir2)then
	F=0.0
       else
	D=ir2-ir1
	F=(ir-ir1)/D
       endif
       it0=(1.0-F)*it1 + F*it2  +0.49999
       amp0=(1.0-F)*a1 + F*a2
       v=(1.0-F)*v1 + F*v2
       time=(1.0-F)*tm1 + F*tm2
       do 70 jr=1,mr
	y=ggd*(jr-ir)
	tm0=sqrt(time*time+y*y/(v*v))
        do 60 it=1,mt
	  x=ggd*(it-it0)
	  t=sqrt(tm0*tm0+x*x/(v*v))
	  tm=1.0+1000.0*t/smp
	  itm=tm + 0.49999
	  ratio=(tm0/t)**tpwr
	  if((itm.lt.1 .or. itm.gt.nt) .and. gssn.le.0.0) go to 60
          call gettrc(jr,it,s)
	  if(gssn.le.0.0)then
            s(lhdr+itm)=s(lhdr+itm) + ratio * amp0
	  else
	    call gaus(s,lhdr,nt,smp,tm,ratio*amp0,gssn)
	  endif
          call puttrc(jr,it,s)
   60   continue
   70  continue
   90 continue
      ir1=ir2
      it1=it2
      tm1=tm2
      a1=a2
      v1=v2
      ibias=1
      if(n.gt.0)go to 15
      go to 10
      end
