c-----------------------------------------------------------------------
c     Migracao 3D pre-empilhamento por registro de no (OBN) do programa
c     test3dbcpmig (SEISPAK, "bcpmig"), para o modulo SPK_3DBCPMIG.
c
c     SPKBSET : tamanhos (fac235 em nx, ny, nt), banda de frequencias,
c               funcao de mute ZMUT e funcao de frequencia HZ
c     SPKBV   : xvel3d com VXT = -3 (ONLVELS) e RFC = 'YES': velocidade
c               media por passo (VELS, grade de saida) e funcao rms VRMS
c     SPKBR   : um registro = xxbin6 (binagem local centrada no no, NMO
c               parcial, FFT) + laco de profundidade de bcpmig para
c               ITYPE 1 (ascendente) e 2 (descendente, espelho), soma no
c               cubo de saida (antes: arquivo 'indata')
c     SPKBT   : xd2t3d com VXT = -3: cubo em profundidade -> tempo
c
c     Diferencas marcadas com SEASEIS:
c       - dado, velocidade, DBASE e cubo de saida em memoria (antes:
c         ONLINE1, ONLVELS, DBASE/shotloc, WSCR1/WINF, velscr, indata)
c       - todas as frequencias num passo (antes: nstepw blocos conforme
c         a memoria; a imagem e linear na frequencia: mesmo resultado)
c       - sem ITYPE 3/4/5 e sem a analise de coerencia (TDAT, OFFLINE4/5)
c-----------------------------------------------------------------------
      SUBROUTINE SPKBSET(NX,NY,NT,APX,NNZ,ITZR,DDZ,DT,FRQ1,FRQ2,AMUTE,
     :                   NXP,NYP,NT1,NLOW,NHIGH,NZ,ZMUT,HZ)
      REAL ZMUT(*),HZ(*)
      INTEGER IFAX(5)
      NXP=NX
      NYP=NY
c   ***   set up nx for pspi (nxx = (2**p1)*(3**p2)*(5**p3) )
      if(apx .gt. 5.)then
      call fac235(nx,ifax)
      nxp=ifax(1)
      end if
c   ***   set up ny for pspi (nyy = (2**p1)*(3**p2)*(5**p3) )
      if(ny .gt. 4 .and. apx .gt. 5.)then
      call fac235(ny,ifax)
      nyp=ifax(1)
      end if
      nz=(nnz-1)/itzr+1
      dz=itzr*ddz
      call fac235(nt,ifax)
      nt1=ifax(1)
      nt2=nt1/2+1
c  ***  determine frequencies to be processed
      ifrq1=2.*dt*nt2*frq1+1
      ifrq2=2.*dt*nt2*frq2+1
      nlow=max0(ifrq1,1)
      nhigh=min0(nt2,ifrq2)
c   ***   zmut function (no 'zmutes' file)
      amute1=amax1(5.,  amin1( 85., amute) )
      tanmute=tan(3.14159*amute1/180.)
      do j=1,nz+1
      zmut(j)=0.
      enddo
      do j=1,nz
       zmute=(j-1)*dz
       rmute=sqrt(1000.*zmute/tanmute)
      zmut(j)=rmute
      enddo
c   ***   frequency function (no 'frequency' file)
      do j=1,nz+1
      hz(j)=frq2
      end do
      RETURN
      END
c
      SUBROUTINE SPKBV(VDAT,NTV,NXV,NYV,DZV,NXO,NYO,NZ,ITZR,DDZ,DT,
     :                 NT,BETA,VELS,VRMS,MAXV,E,X,XX)
c     xvel3d, ramo vel(1) = -3 (ONLVELS) com rfc = 'YES'
      REAL VDAT(NTV,*),VELS(NXO*NYO,*),VRMS(*),E(*),X(*),XX(*)
      dz=itzr*ddz
      np=ntv
      nx=nxo
      ny=nyo
      idpn1=1.
      nxy=nx*ny
      rnxy=1./nxy
      do n=1,2*nz
      vrms(n)=0.
      enddo
      do 10 i=1,ny
      do 10 k=1,nx
      ix=k-idpn1+beta+1
      ix=max0(1, min0(nxv,ix) )
      iy=min0(i,nyv)
c     SEASEIS: gettrc1(iy,ix) do ONLVELS -> VDAT
      itr=(iy-1)*nxv+ix
      xx(1)=ix
      xx(2)=iy
      do j=1,ntv
      xx(2*j+1)=vdat(j,itr)
      xx(2*j+2)=(j-1)*dzv
      enddo
      xx(2*ntv+3)=0
      ii=3
      eps=.1
111   continue
      if(xx(ii) .lt. eps)then
      ii=ii+2
c     SEASEIS: traco de velocidade todo nulo (antes: laco infinito)
      if(ii .gt. 2*np+1)goto 112
      goto 111
      end if
112   continue
      v1=xx(ii)
      v2=v1
      z=0.
      z1=0.
      z2=xx(ii+1)
      eps=.1
      do 101 n=1,nz
      rvbar=0.
      do 104 l=1,itzr
102   continue
      if(ii.gt.2*np+2.or.z.le.z2)goto 103
      if(xx(ii).lt.eps .or. (z .gt. xx(ii+1).and.ii .lt. 2*np))then
      ii=ii+2
      goto 102
      end if
      v1=v2
      v2=xx(ii)
      z1=z2
      z2=xx(ii+1)
103   continue
      vvel=v2
c     rfc = 'YES'
      z=z+ddz
      rvbar=rvbar+1./vvel
104   continue
      e(n)=itzr/rvbar
      vrms(n)=vrms(n)+rnxy*e(n)
101   continue
      irec=(i-1)*nx+k
c     SEASEIS: 'scr' + trni2da -> 'velscr' : VELS(xy, nivel)
      do n=1,nz
      vels(irec,n)=e(n)
      enddo
10    continue
c   ***   return average rms function in vel
      vrms2=vrms(1)*vrms(1)
      tt=0.
      do k=1,nz
      delt=2.*dz/vrms(k)
      vrms2 = ( tt*vrms2 + delt*vrms(k)*vrms(k) )/(tt+delt)
      e(k)=sqrt(vrms2)
      x(k)=tt
      tt=tt+delt
      end do
      do k=1,nt
      e(nt+k)=(k-1)*dt
      end do
      ntout=x(nz)/dt
c     SEASEIS: ntout limitado ao tamanho de VRMS
      ntout=min0(ntout,maxv)
      call spl3(x,e,e(nt+1),vrms,nz,ntout)
      if(ntout .lt. nt)then
      do k=ntout+1,nt
      vrms(k)=vrms(ntout)
      end do
      end if
      return
      end
c
      SUBROUTINE SPKBR(D,NTD,NTR,SXH,SYH,IGX,IGY,IGZ,IP,RP,VELS,VRMS,
     :                 ZMUT,HZ,CUBE,WK,DATA,TR1,YT,XNDX,XNDX1,RR,INFO,
     :                 IFOC,TDAT)
c     IFOC: analise de focalizacao (ngsnx ngsny iforx ifory irnix
c           irniy jnt), ngsnx*ngsny = 0: sem analise; TDAT acumulado
c     IP:  nx ny nt nt1 nlow nhigh nz itzr nxout nyout itype
c     RP:  dx dy dt ddz apx gam1 xll0 yll0 xlr0 ylr0 eps
c     INFO(1..4) = ixout iyout ok(1)/fora(0) nlive ; RINFO via INFO(5..)
      REAL D(NTD,*),SXH(*),SYH(*),VELS(*),VRMS(*),ZMUT(*),HZ(*)
      REAL CUBE(*),WK(*),DATA(*),TR1(*),XNDX(*),XNDX1(*),RR(*)
      INTEGER IP(*),INFO(*),IFOC(*)
      REAL RP(*),TDAT(*)
      COMPLEX YT(*)
      nx=ip(1)
      ny=ip(2)
      nt=ip(3)
      nt1=ip(4)
      nlow=ip(5)
      nhigh=ip(6)
      nz=ip(7)
      itzr=ip(8)
      nxout=ip(9)
      nyout=ip(10)
      itype=ip(11)
      dx=rp(1)
      dy=rp(2)
      dt=rp(3)
      ddz=rp(4)
      apx=rp(5)
      gam1=rp(6)
      xll0=rp(7)
      yll0=rp(8)
      xlr0=rp(9)
      ylr0=rp(10)
      eps=rp(11)
c     planos de frequencia em y(iy), iy = 1+2*nxy (complexo) do SPKBRW
      CALL SPKBB(D,NTD,NTR,SXH,SYH,IGX,IGY,IGZ,WK(4*nx*ny+1),nx,ny,nt,
     : nt1,nlow,nhigh,xll0,yll0,xlr0,ylr0,dx,dy,dt,xll1,yll1,xlr1,ylr1,
     : sx,sy,sz,vrms,eps,TR1,YT,XNDX,XNDX1,RR,nlive)
      INFO(4)=nlive
      CALL SPKBRW(WK,WK,DATA,VELS,ZMUT,HZ,CUBE,nx,ny,nt,nt1,nlow,nhigh,
     : nz,itzr,nxout,nyout,itype,dx,dy,dt,ddz,apx,gam1,xll0,yll0,
     : xll1,yll1,xlr1,ylr1,sx,sy,sz,INFO,IFOC,TDAT)
      RETURN
      END
c
      SUBROUTINE SPKBB(D,NTD,NTR,XHDR,YHDR,IGX,IGY,IGZ,Y,nx,ny,ntt,nt1,
     : nlow,nhigh,xll0,yll0,xlr0,ylr0,dx,dy,dt,xll,yll,xlr,ylr,gx,gy,gz,
     : vel,eps,TR1,YT,XNDX,XNDX1,RR,NLIVE)
c     xxbin6: binagem do registro na grade local nx x ny centrada no no
c     (DBASE do bcdbasefix: xmid = gx, ymid = gy, gz = GZEL), NMO parcial
c     para o centro do bin, FFT; plano de cada frequencia em Y
c     (antes: WSCR1 + ctrni2da -> WINF)
      REAL D(NTD,*),XHDR(*),YHDR(*),VEL(*),TR1(*),XNDX(*),XNDX1(*),RR(*)
      COMPLEX Y(*),YT(*)
c   ***   DBASE (bcdbasefix): GX GY GZEL do traco 1, xmid = gx, ymid = gy
      gx=igx
      gy=igy
      gz=igz
      xmid1=igx
      ymid1=igy
      nt=min0(ntd,ntt)
      nn=nx*ny
      nws=nhigh-nlow+1
c   ***   set up directional unit vectors
      r=sqrt( (xlr0-xll0)**2 + (ylr0-yll0)**2)
      ux=(xlr0-xll0)/r
      uy=(ylr0-yll0)/r
      dxx=ux*dx
      dxy=-uy*dy
      dyx=uy*dx
      dyy=ux*dy
      xll=xmid1-(nx/2)*dxx-(ny/2)*dxy
      yll=ymid1-(nx/2)*dyx-(ny/2)*dyy
c     SEASEIS: sem divisao por zero quando a grade e N-S ou E-W
      ixx=0
      iyy=0
      if(abs(dxx).gt.1.e-3)ixx=nint( (xll-xll0)/dxx )
      if(abs(dyy).gt.1.e-3)iyy=nint( (yll-yll0)/dyy )
      xll=xll0+ixx*dxx
      yll=yll0+iyy*dyy
      xlr=xll+(nx-1)*dxx
      ylr=yll+(nx-1)*dyx
c   ***   first pass to set close bin locations
      factor=2./nt1
      do j=1,nx*ny
      rr(j)=1000000.
      xndx(j)=0
      xndx1(j)=0
      end do
c.  begin binning process
      do  j=1,ntr
      sx=xhdr(j)
      sy=yhdr(j)
      xx=(xhdr(j)-xll)*ux+(yhdr(j)-yll)*uy
      yy=(xhdr(j)-xll)*(-uy)+(yhdr(j)-yll)*ux
      xx=xx/dx
      yy=yy/dy
      ix=nint(xx)+1
      iy=nint(yy)+1
      ii=(iy-1)*nx+ix
       if(1.le.ix .and. ix.le.nx .and. 1 .le.iy .and. iy.le.ny)then
      xxx=xll+ (ix-1)*dxx + (iy-1)*dxy
      yyy=yll+ (ix-1)*dyx + (iy-1)*dyy
      xndx1(ii)=j
      xndx(ii)=j
      dr=sqrt((xxx-sx)**2+(yyy-sy)*2)
      rr(ii)=dr
      endif
      enddo
      do j=1,nx*ny
      xndx(j)=xndx1(j)
      enddo
c   ***   fill in empty bins with traces from neighboring bins
c   ***   make nkk passes nkk=eps/dx+1
      nkk=eps/dx+.5
      do kk=1,nkk
      do iy=1,ny
      do ix=1,nx
       ii=(iy-1)*nx+ix
      if(xndx1(ii) .eq. 0) then
c   ***   find nearest neighboring trace
      do l1=-1,1
      do l2=-1,1
      if(l1 .ne. 0 .or. l2 .ne. 0)then
      ixy=ii+l1*nx+l2
      if(1.le.ixy .and. ixy.le.nx*ny)then
      iii=xndx1(ixy)
      if(iii .gt. 0)then
      xxx=xll+(ix-1)*dxx+(iy-1)*dxy
      yyy=yll+(ix-1)*dyx+(iy-1)*dyy
      sx= xhdr(iii)
      sy= yhdr(iii)
      dr= sqrt( (xxx-sx)**2 + (yyy-sy)**2 )
      if(dr.lt.rr(ii) .and. dr.lt.eps)then
      rr(ii)=dr
      xndx(ii)=iii
      endif
      end if
      end if
      end if
      end do
      end do
      end if
      end do
      end do
      do j=1,nx*ny
      xndx1(j)=xndx(j)
      enddo
      enddo
c   ***   traces of the bins: partial moveout, fft
      nlive=0
      do iy=1,ny
      do ix=1,nx
      ixy=(iy-1)*nx+ix
      xx=xll+(ix-1)*dxx+(iy-1)*dxy
      yy=yll+(ix-1)*dyx+(iy-1)*dyy
      do jj=1,4*ntt+20
      tr1(jj)=0.
      enddo
      if(xndx(ixy).gt.0)then
      nlive=nlive+1
      itrc=xndx(ixy)
      do jj=1,nt
      tr1(jj)=d(jj,itrc)
      enddo
c   ***   find max value for trace
      amax=0
      do jj=1,nt
      amax=amax1(amax,abs(tr1(jj)))
      enddo
c   ***   apply partial movout to bin center
      sx=xhdr(itrc)
      sy=yhdr(itrc)
      rng1=sqrt((gx-sx)**2+(gy-sy)**2)
      rng2=sqrt((gx-xx)**2+(gy-yy)**2)
      call xxnmo(tr1,tr1(ntt+10),ntt,dt,rng1,rng2,vel)
c   ***   fft data; ffted data is scaled by 2/nt1
      factor=2./nt1
      call clear(yt,2*nt1)
c   ***   check for and mute bad values
      do jj=1,nt
      if(abs(tr1(jj)).gt.amax)tr1(jj)=0.
      enddo
      do jj=1,nt
      yt(jj)=factor*tr1(jj)
      end do
      else
      call clear(yt,2*nt1)
      end if
      call fftc(yt,yt(nt1+1),2,2*nt1,nt1,1,-1)
c     SEASEIS: frequencias nlow..nhigh direto no plano de cada frequencia
      do iw=1,nws
      y((iw-1)*nn+ixy)=yt(nlow+iw-1)
      enddo
      enddo
      enddo
      return
      end
c
      SUBROUTINE SPKBRW(Y,S,DATA,VELS,ZMUT,HZ,CUBE,nx,ny,nt,nt1,nlow,
     : nhigh,nz,itzr,nxout,nyout,itype,dx,dy,dt,ddz,apx,gam1,xll,yll,
     : xll1,yll1,xlr1,ylr1,sx,sy,sz,INFO,IFOC,TDAT)
c     laco de profundidade de bcpmig (ITYPE 1 e 2); Y (complexo) e S
c     (real) sao o mesmo vetor, como a EQUIVALENCE (s(1),y(1))
      COMPLEX Y(*)
      REAL S(*),DATA(*),VELS(NXOUT*NYOUT,*),ZMUT(*),HZ(*)
      REAL CUBE(NXOUT*NYOUT,*)
      INTEGER INFO(*),IFOC(*)
      REAL TDAT(*)
c     SEASEIS: analise de velocidade por focalizacao (test3dbcpmig,
c     FORX LORX FORY LORY RNIX RNIY WINDOW): FDAT = 'data' do original
c     (um nivel), WF = y(ia) (espectro + trabalho da fftc), WB = y(ib)
      REAL, ALLOCATABLE :: FDAT(:)
      COMPLEX, ALLOCATABLE :: WF(:),WB(:)
c     SEASEIS: OpenMP no laco de frequencias. Cada thread tem seu
c     proprio scratch a,b,c (antes y(ia),y(ib),y(ic)) em WTH.
      COMPLEX, ALLOCATABLE :: WTH(:)
!$    INTEGER omp_get_thread_num, omp_get_max_threads
      nxy=nx*ny
      nthr=1
!$    nthr=omp_get_max_threads()
      nsc=4*nxy+17*max(nx,ny)+64
      allocate(wth(nsc*nthr))
      ith=0
      nvels=nint(apx-10.)
      ngsnx=ifoc(1)
      ngsny=ifoc(2)
      iforx=ifoc(3)
      ifory=ifoc(4)
      irnix=ifoc(5)
      irniy=ifoc(6)
      jnt=ifoc(7)
      ngsns=ngsnx*ngsny
      sr=1000.*dt
      twindow=(jnt-1)*sr
      ilorx=iforx+(ngsnx-1)*irnix
      ilory=ifory+(ngsny-1)*irniy
      if(ngsns.gt.0)then
      allocate(fdat(jnt*ngsns),wf(22*nt1+10),wb(nhigh-nlow+12))
      endif
      nwdst=nx*ny*itzr
      dz=itzr*ddz
      delw=2.*3.1415927/nt1
      nt2=nt1/2+1
      if(dy.eq.0.)dy=dx
      con2=(dt/dy)**2
      nl=nlow
      nh=nhigh
      nws=nhigh-nlow+1
      nxyws=nxy*nws
      r=sqrt( (xlr1-xll1)**2 + (ylr1-yll1)**2)
      ux=(xlr1-xll1)/r
      uy=(ylr1-yll1)/r
      xxsx=(sx-xll1)*ux+(sy-yll1)*uy
      xxsx=xxsx/dx+1.
      yysy=(sx-xll1)*(-uy)+(sy-yll1)*ux
      yysy=yysy/dy+1.
      alph1=xxsx
      alph2=yysy
      alph3=abs(sz)
      info(3)=0
      if(nint(alph1) .lt. 1 .or. nint(alph1) .gt. nx)return
      if(nint(alph2) .lt. 1 .or. nint(alph2) .gt. ny)return
      info(3)=1
      ixout= nint( ((xll1-xll)*ux + (yll1-yll)*uy )/dx + 1.)
      iyout= nint( (-(xll1-xll)*uy + (yll1-yll)*ux )/dy + 1.)
      info(1)=ixout
      info(2)=iyout
c   ***   compute y (complex) indicies
      itu1=1
      itu2=itu1+nxy
      itd1=itu2+nxy
      itd2=itd1+nxy
       iy=1+2*nxy
      iz=iy+nxyws
      ia=iz+nxyws
      ib=ia+nxy
      ic=ib+nxy
      id=ic+nxy
c   ***   s (real) indicies start at end of complex
      isx=2*(id+nxy)+1
      isy=isx+nxy
      ivx=isy+nxy
      ivy=ivx+nxy
      is=ivy+nxy
c   ***   begin downward continuation loop
      sgn=+1.
      xkshot=+alph3/dz
      kshot=xkshot
      call movlev(y(iy),y(iz),2*nxyws)
      kstart=gam1/dz+1
      kstartu=amax1(gam1/dz+1,alph3/dz+1)
      info(5)=kstart
      info(6)=kstartu
      info(7)=kshot
c.  new3dvspmig code
      call clear(s,4*nxy)
      do 2000 k=1,nz
      call clear(data,nwdst)
c   ***   get velocity
      kk0=k
      con1=dz/dt
c     SEASEIS: velscr -> VELS
      call extractv(vels(1,kk0),s(ivy),nxout,nyout,nx,ny,ixout,iyout)
c   ***   get mute velocity
       if(k .eq. 1)then
       ii=(alph2-1)*nx+alph1
       vmute=s(ivy+ii)
       vstart=vmute
       end if
c   ***
       if(kk0 .le. kstart)then
       xx=con1/vstart
       vv=con2*vstart**2
       do j=1,nxy
       s(ivy+j-1)=vstart
       end do
       end if
c   ***   compute arrival times for direct reflected energy
      xkk=kk0-xkshot
      if(itype.eq.1 .or. itype .eq. 3)then
      call xpgsarr3d2(s(itu1),s(itu2),s(is),s(ivy),xkk,nx,itzr,
     :                dt,ddz,dx,alph1,ny,dy,alph2)
      endif
c   ***   compute arrival times for mirror reflected energy
c   ***   first compute arrival times from location of receiver to
c   ***   surface of water
      if(itype.eq.2 .or. itype.eq.3)then
      if(k .eq. 1)then
      do kk=1,kshot
      xkk=kk+xkshot-kshot
      call xpgsarr3d1(s(itd1),s(itd2),s(is),s(ivy),xkk,nx,itzr,
     :                dt,ddz,dx,alph1,ny,dy,alph2)
      enddo
      end if
c   ***   continue from surface of water downward
      xkk=kk0+xkshot
      call xpgsarr3d1(s(itd1),s(itd2),s(is),s(ivy),xkk,nx,itzr,
     :                dt,ddz,dx,alph1,ny,dy,alph2)
      endif
c   ***
c   ***   downward continue wavefield from recording surface
c   ***
      ii=(alph2-1)*nx+alph1
      do 300 j=0,nxy-1
      s(isy+j)=con1/s(ivy+j)
300   s(ivy+j)=con2*s(ivy+j)**2
      if(kk0 .eq. kstart)then
      xx=kk0*con1/vstart
      vv=con2*vstart**2
      do j=1,nxy
      s(ivx+j-1)=s(ivy+j-1)
      s(isx+j-1)=s(isy+j-1)
      s(ivy+j-1)=vv
      s(isy+j-1)=xx
      nh1=2*dt*nt2*hz(k)+1
      end do
c     SEASEIS: "goto 991" (sair quando i > nh1-nl) virou um IF,
c     para o laco poder ser paralelo; frequencias sao independentes
!$omp parallel do default(shared) private(i,w,i1,ith,ja)
!$omp& schedule(dynamic,1)
      do 100 i=1,nws
      if(i.le.nh1-nl .and. apx .ge. 11.)then
      w=(nl+i-2)*delw
      i1=(i-1)*nxy
!$    ith=omp_get_thread_num()
      ja=1+ith*nsc
c   ***   velocity term s(isy)= (v*dt/dy)**2
      call oxxmvphs2d(y(iy+i1),wth(ja),wth(ja+nxy),wth(ja+2*nxy),
     *      nx,ny,dx,dy,dt,w,s(isy),s(ivy),sgn,nvels)
      end if
100   continue
!$omp end parallel do
      call vshift1(y(iy),y(ia),y(ib),y(ic),s(isy),delw,nl,nh1,nxy,sgn)
c   ***   reset velocity coefficients for this level
      do j=1,nxy
      s(ivy+j-1)=s(ivx+j-1)
      s(isy+j-1)=s(isx+j-1)
      end do
      end if
c   ***   if kk0 .gt. kstart then procede with downward continuation
      if(kk0 .gt. kstart)then
      nh1=2*dt*nt2*hz(k)+1
!$omp parallel do default(shared) private(i,w,i1,ith,ja)
!$omp& schedule(dynamic,1)
      do 200 i=1,nws
      if(i.le.nh1-nl .and. apx .ge. 11.)then
      w=(nl+i-2)*delw
      i1=(i-1)*nxy
!$    ith=omp_get_thread_num()
      ja=1+ith*nsc
      call oxxmvphs2d(y(iy+i1),wth(ja),wth(ja+nxy),wth(ja+2*nxy),
     *      nx,ny,dx,dy,dt,w,s(isy),s(ivy),sgn,nvels)
c   ***   apply attenuating absorbing boundary
      call obndry(y(iy+i1),nx,ny,20,20)
      end if
200   continue
!$omp end parallel do
c   ***   apply thin lens term (refraction correction)
      call osgivshift1(y(iy),s(isy),delw,nl,nh1,nxy,sgn)
      end if
c   ***   image data for direct reflected energies
      if(itype.eq.1 .and.kk0.gt.kstartu)then
      call osgimrg2r(data,y(iz),y(iy),s(isy),s(itu1),s(itu2),nl,nh1,
     *  nt1,nxy,itzr)
c   ***   mute data at mute angle
       zz=(k-1)*dz
       call zmutes(data,zmut,nx,ny,dx,dy,ddz,itzr,zz,alph1,alph2,k)
c     SEASEIS: soma no cubo de saida (antes: arquivo 'indata')
      do jj=1,itzr
      i1=(jj-1)*nxy+1
      i2=(kk0-1)*itzr+jj
      call embed(data(i1),cube(1,i2),nx,ny,nxout,nyout,ixout,iyout)
      end do
      endif
c   ***   image data for mirror reflected energies
       if(itype.eq.2 .and. kk0.gt.kstart)then
      call osgimrg2r(data,y(iz),y(iy),s(isy),s(itd1),s(itd2),nl,nh1,
     *  nt1,nxy,itzr)
c   ***   mute data at mute angle
       zz=(k-1)*dz
       call zmutes(data,zmut,nx,ny,dx,dy,ddz,itzr,zz,alph1,alph2,k)
c     SEASEIS: soma no cubo de saida (antes: arquivo 'indata')
      do jj=1,itzr
      i1=(jj-1)*nxy+1
      i2=(kk0-1)*itzr+jj
      call embed(data(i1),cube(1,i2),nx,ny,nxout,nyout,ixout,iyout)
      end do
      endif
c   ***   SEASEIS: analise de focalizacao (velocity focussing) de
c   ***   test3dbcpmig, ramo itype 1/2: em cada ponto de analise da
c   ***   grade de saida, o campo continuado y(iz) e deslocado pelo
c   ***   tempo de chegada menos meia janela e levado ao tempo; a
c   ***   janela (jnt amostras) e somada em TDAT (nivel kk0)
      if(ngsns .eq. 0)goto 771
      call clear(fdat,jnt*ngsns)
      do jy=1,ny
      do jx=1,nx
      jxx=jx+ixout-1
      jyy=jy+iyout-1
      modx=mod(jxx-iforx,irnix)
      mody=mod(jyy-ifory,irniy)
      if( (jxx.ge.iforx .and. jxx.le.ilorx .and. modx.eq.0) .and.
     *    (jyy.ge.ifory .and. jyy.le.ilory .and. mody.eq.0) )then
      jgsny=(jyy-ifory)/irniy +1
      jgsnx=(jxx-iforx)/irnix +1
      jgsn=(jgsny-1)*ngsnx +jgsnx
      itdat=(jgsn-1)*jnt+1
      if(kk0.gt.kstart)then
      call clear(wf,2*nt1)
      ii=nl
      jxy= (jy-1)*nx + jx
      do jj=0,nws-1
      wf(ii+jj)=y(iz+ jxy-1 + jj*nxy)
      end do
c   ***   select arrival time for this gsn location
      if(itype.eq.1) shf=s(itu1+jxy-1)-.5*twindow/sr
      if(itype.eq.2) shf=s(itd1+jxy-1)-.5*twindow/sr
c   ***   shift complex data by this amount
      call vshift1(wf(ii),wb,wb,wb(nws+3),shf,delw,nl,nh,1,1.)
c   ***   inverse transform
      call fftc(wf,wf(nt1+1),2,2*nt1,nt1,1,+1)
      do jj=1,jnt
      fdat(itdat+jj-1)= fdat(itdat+jj-1)-xx*real(wf(jj))
      end do
      endif
      end if
      end do
      end do
c   ***   sum this level of velocity coherence info to TDAT
      i0=(kk0-1)*jnt*ngsns
      do j=1,ngsns
      itdat=(j-1)*jnt
      do jj=1,jnt
      tdat(i0+itdat+jj)=tdat(i0+itdat+jj)+fdat(itdat+jj)
      enddo
      tdat(i0+itdat+jnt)=fdat(itdat+jnt)
      enddo
771   continue
      call movlev(y(iy),y(iz),2*nxyws)
2000   continue
      RETURN
      END
c
      SUBROUTINE SPKBT(X,U,E,S,VDAT,NTV,NXV,NYV,DZV,NT,SR,NNZ,DDZ,
     :                 DPNX,DPNY,V)
c     xd2t3d, ramo vtype = -3 (ONLVELS): traco em profundidade X -> tempo U
      REAL X(*),U(*),E(*),S(*),VDAT(NTV,*),V(*)
      np=ntv
      iy=max0(1,min0(nint(dpny),nyv))
      ix=max0(1,min0(nint(dpnx),nxv))
c     SEASEIS: gettrc3 do ONLVELS -> VDAT
      itr=(iy-1)*nxv+ix
      v(1)=ix
      v(2)=iy
      do j=1,ntv
      v(2*j+1)=vdat(j,itr)
      v(2*j+2)=(j-1)*dzv
      enddo
      v(2*ntv+3)=0
      ii=3
      eps=.1
      ii0=ii+2*np
      rddz=1./ddz
      dt2=.0005*sr
111   continue
      if(v(ii) .lt. eps)then
      ii=ii+2
c     SEASEIS: traco de velocidade todo nulo (antes: laco infinito)
      if(ii .gt. ii0)goto 112
      goto 111
      end if
112   continue
      v1=v(ii)
      v2=v1
      z=0.
      z1=0.
      z2=v(ii+1)
      t=0.
      eps=.1
      do 101 n=1,nt
102   continue
      if(ii.gt.ii0.or.z.le.z2)goto 103
      if(v(ii) .lt. eps .or. z .gt. v(ii+1))then
      ii=ii+2
      goto 102
      end if
      v1=v2
      v2=v(ii)
      z1=z2
      z2=v(ii+1)
      goto 102
103   continue
c     vint nao e 'yes'
      vel=v2
      e(n)=z*rddz
      z=z+dt2*vel
      t=t+sr
101   continue
c   ***   spl1 resamples data trace in x at sample points defined by e
      call spl1(x,1.,u,e,s,nnz,nt,1.)
      return
      end
      SUBROUTINE SPKBTHR(N)
c     SEASEIS: numero de threads OpenMP (N <= 0: OMP_NUM_THREADS ou
c     todos os nucleos); devolve em N o numero usado
!$    INTEGER omp_get_max_threads
!$    if(n .gt. 0)call omp_set_num_threads(n)
      n=1
!$    n=omp_get_max_threads()
      return
      end
