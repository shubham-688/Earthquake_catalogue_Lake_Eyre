       character*40 ifile,pfile
       character*60 ctext
       dimension itext(10)
       common /p3d/ xpr,ypr,caxt,saxt,cayt,sayt
       real xlat(1000),xlon(1000),xdep(1000),xsec(1000),xmis(1000)
       real alat,alon,adep,asec,amisf
       real qlat,qlon,qdep,qsec
       real xfmin,xfmax,xfbas
       real vlln(12),vlle(12),vlls(12)
       real vadn(12),vadd(12),vads(12)
       real vode(12),vodd(12),vods(12)
       equivalence (itext(1),ctext)  
c 
       dgrad = atan(1.0)/45.
       sq2 = 1.0/sqrt(2.)
c
       dsc = 0.5
       tsc = 2.5
       write(6,*) "Location data file"
       read(5,*) ifile
       open(7,file=ifile)
       write(6,*) "event params - lat, long, dep, sec"
       read(5,*)  olat,olon,odep,osec
       write(6,*) olat,olon,odep,osec
       write(6,*) "ranges km/sec - N, E, dep, sec"
       read(5,*)  rn,re,rd,rs
       write(6,*) rn,re,rd,rs
       write(6,*) "number of iterations"
       read(5,*) niter
       write(6,*) niter
       write(6,*) "style of weight Boltzman/FD/BE (0/1/2)"
       read(5,*) nwe
       write(6,*) nwe
       write(6,*) "beta"
       read(5,*) beta
       write(6,*) beta
       write(6,*) "Reference misfit factor"
       read(5,*) xfref
       write(6,*) xfref
       write(6,*) "Threshold factor"
       read(5,*) xthr
       write(6,*) xthr
       write(6,*) "plot file"
       read(5,*) pfile
       write(6,*) pfile
c
       sn = 4.0/rn
       se = 4.0/re
       dl = 4.0*rd/re
       asc = sn*111.19
       bsc = se*111.19*cos(olat*dgrad)
       scl = cos(olat*dgrad)
       dsc = sn
       sd = dsc
       tsc = 5.0/rs
       d1 = 0.5*dl
c
       read(7,fmt='(a)') ctext
       l = 1
       DO k=0,niter
         read(7,*)
         read(7,*) 
         read(7,*) 
         read(7,*)
         DO j=1,9 
           read(7,*) n,alat,alon,adep,asec,amisf
           xlat(l) = alat
           xlon(l) = alon
           xdep(l) = adep
           xsec(l) = asec
           xmis(l) = amisf
*           write(9,*) n,alat,alon,adep,asec,amisf
*           write(9,*) l,xlat(l),xlon(l),xdep(l),xsec(l),xmis(l)
           l = l+1      
         ENDDO
       ENDDO 
c 
       lmod = l-1
       flmod = 1.0/float(lmod)
c
        xfmax = xmis(l)
        xfmin = xmis(1)
        mopt = 1
        do j=1,lmod
          if(xfmin .gt. xmis(j)) then
            xfmin = xmis(j)
            mopt = j
          end if
          if(xfmax .lt. xmis(j)) xfmax = xmis(j)
        end do
        xfbas = xfmin*0.9999
*        write(9,*) mopt, xfbas, xfmin
        qlat = xlat(mopt)
        qlon = xlon(mopt)
        qdep = xdep(mopt)
        qsec = xsec(mopt)
        xfu = xfmax
        if(xfu.gt.5.) xfu =5.
        ama = xfbas
*        amb = (xfu-xfmin)*0.25+xfbas
        amb = xfbas*xfref
        amc = amb-ama
        write(6,*) "misfit bounds"
        write(6,*) ama,amb,amc
c
c       weighted variances
        vlat = 0.
        vlon = 0.
        vdep = 0.
        vsec = 0.
        ewet = 0.
        l1 = lmod/2.
        do j=l1,lmod
          if(nwe.eq.0) then
             we = exp(-beta*(xmis(j)-xfbas)/xfbas)
          elseif(nwe.eq.1) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas*xfref)/xfbas)+1.0)
          elseif(nwe.eq.2) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas/xfref)/xfbas)-1.0)
          endif
          vlat = vlat + we*(xlat(j)-elat)**2
          vlon = vlon + we*(xlon(j)-elon)**2
          vdep = vdep + we*(xdep(j)-edep)**2
          vsec = vsec + we*(xsec(j)-esec)**2
          ewet = ewet + we
        enddo
        xw = 1.0/ewet
        slat = sqrt(vlat)*xw*111.19
        slon = sqrt(vlon)*xw*111.19*cos(olat*dgrad)
        sdep = sqrt(vdep)*xw
        ssec = sqrt(vsec)*xw
c
c       measure of fit?
        if(nwe.eq.0) then
           wef = exp(-beta*(xfmin-xfbas)/xfbas)
        elseif(nwe.eq.1) then
           wef = 1.0/(exp(beta*(xfmin-xfbas*xfref)/xfbas)+1.0)
        elseif(nwe.eq.2) then
           wef = 1.0/(exp(beta*(xfmin-xfbas/xfref)/xfbas)-1.0)
        endif
        wthres = 0.5+xthr*(wef-0.5)
        write(6,*) "wthres", wthres
        do k=1,12
          ang = float(k-6.5)*30.0*dgrad
          cang = cos(ang)
          sang = sin(ang)
          vlln(k) = 0.25*cang
          vlle(k) = 0.25*sang
          vlls(k) = 0.25
          vadn(k) = 0.25*cang
          vadd(k) = 0.25*sang
          vads(k) = 0.25
          vode(k) = 0.25*cang
          vodd(k) = 0.25*sang
          vods(k) = 0.25
        enddo
        usec = 0.
        dsec = 0.
        do j =l1,lmod
          if(nwe.eq.0) then
             we = exp(-beta*(xmis(j)-xfbas)/xfbas)
          elseif(nwe.eq.1) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas*xfref)/xfbas)+1.0)
          elseif(nwe.eq.2) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas/xfref)/xfbas)-1.0)
          endif
          if(we.gt.wthres) then
            if(qsec-xsec(j) .gt.dsec) dsec =  qsec-xsec(j)     
            if(xsec(j)-qsec .gt.usec) usec =  xsec(j)-qsec                       
            sx = (xlon(j)-qlon)*scl*111.19
            sy = (xlat(j)-qlat)*111.19
            sz = xdep(j)-qdep
            ss = sqrt(sx*sx+sy*sy)
            sr = sqrt(sx*sx+sz*sz)
            st = sqrt(sy*sy+sz*sz)
            su = sqrt(sx*sx+sy*sy+sz*sz)
            rll = atan2(sy,sx)
            rad = atan2(sz,sy)
            rod = atan2(sz,sx)
            write(9,*) j,sx,sy,sz,rll/dgrad,rad/dgrad,rod/dgrad
            do k=1,12
              ta = float(k-7)*30.0*dgrad
              tb = ta+30.0*dgrad
              if(rll.ge.ta .and. rll.lt.tb) then
                if(ss. gt. vlls(k)) then
                  vlls(k) = ss
                  vlle(k) = sx
                  vlln(k) = sy
                endif
              endif
               if(rod.ge.ta .and. rod.lt.tb) then
                if(sr. gt. vods(k)) then
                  vods(k) = sr
                  vode(k) = sx
                  vodd(k) = sz
                endif
              endif
              if(rad.ge.ta .and. rad.lt.tb) then
                if(st. gt. vads(k)) then
                  vads(k) = st
                  vadn(k) = sy
                  vadd(k) = sz
                endif
              endif             
            enddo
          endif
        enddo
        dlat = dlat*111.19
        ulat = ulat*111.19
        dlon = dlon*111.19*scl
        ulon = ulon*111.19*scl
c
       open(8,file=pfile)
       call Hplots(1,0,8,0)
       call typset(0.0,0.0)
       call zpick(3,0,is)
       x1r = 10.0
       y1r = 18.0
       x2r =  9.0
       y2r = 12.0
       x3r =  9.0
       y3r =  5.0
       x4r = 16.0
       y4r = 12.0
c
       call pen(2,0)
       call plot(x1r-6.,y1r-0.25,3)
       call edgert(12.,0.5,0)
       call typnum(x1r-6.9,y1r-0.1,0.30,-rs,0.,-1)
       call typnum(x1r+6.1,y1r-0.1,0.30, rs,0.,-1)
       call typstr(x1r+7.5,y1r-0.1,0.30,"[s]",0.,3)
       call pen(15,0)
       call plot(x1r-3.0,y1r-0.25,3)
       call plot(x1r-3.0,y1r+0.25,2)
       call plot(x1r    ,y1r-0.25,3)
       call plot(x1r    ,y1r+0.25,2)
       call plot(x1r+3.0,y1r-0.25,3)
       call plot(x1r+3.0,y1r+0.25,2)
       xpr = 12.0
       ypr = 10.0
       caxt = cos(25.*dgrad)
       saxt = sin(25.*dgrad)
       cayt = cos(75.*dgrad)
       sayt = sin(75.*dgrad)
       xt = xpr
       yt = ypr+12.
c
c  plot  consistency patch
        call pen(2,0)
        xsa = asc*(qlat-olat)
        ysa = bsc*(qlon-olon)
        zsa = dsc*(qdep-odep)
        xsb = tsc*(qsec-osec)
c
        pdse = dsec*tsc
        puse = usec*tsc
c
        call pen(31,0)
         xpa = xsa+1.05*vlle(1)*se
         ypa = ysa+1.05*vlln(1)*sn
         zpa = zsa
         call plotned(-xpa,ypa,zpa,3)
        do k=1,12
         xpa = xsa+1.05*vlle(k)*se
         ypa = ysa+1.05*vlln(k)*sn
         zpa = zsa
         call plotned(-xpa,ypa,zpa,2)
        enddo        
         xpa = xsa+1.05*vlle(1)*se
         ypa = ysa+1.05*vlln(1)*sn
         zpa = zsa
         call plotned(-xpa,ypa,zpa,2)
        write(8,*) "cf"   
c     
        call pen(31,0)
         xpa = xsa
         ypa = ysa+vode(1)*se
         zpa = zsa+vodd(1)*sd
         call plotned(-xpa,ypa,zpa,3)
        do k=1,12
         xpa = xsa
         ypa = ysa+vode(k)*se
         zpa = zsa+vodd(k)*sd   
         call plotned(-xpa,ypa,zpa,2)
        enddo     
         xpa = xsa
         ypa = ysa+vode(1)*se
         zpa = zsa+vodd(1)*sd   
         call plotned(-xpa,ypa,zpa,2)
        write(8,*) "cf"        
c        
        call pen(31,0)
         xpa = xsa+vadn(1)*sn
         ypa = ysa
         zpa = zsa+vadd(1)*sd
         call plotned(-xpa,ypa,zpa,3)
        do k=1,12
         xpa = xsa+vadn(k)*sn
         ypa = ysa
         zpa = zsa+vadd(k)*sd
         call plotned(xpa,ypa,zpa,2)
        enddo   
         xpa = xsa+vadn(1)*sn
         ypa = ysa
         zpa = zsa+vadd(1)*sd     
         call plotned(-xpa,ypa,zpa,2)
        write(8,*) "cf" 
c
       call pen(19,0)
       call box3ned( 4., 4., 4.,"N","E","Z",0)
       call box3ned(-4., 4., 4.," "," "," ",0)
       call box3ned( 4.,-4., 4.," "," "," ",0)
       call box3ned(-4.,-4., 4.," "," "," ",0)
       call box3ned( 4., 4.,-4.," "," "," ",0)
       call box3ned(-4., 4.,-4.," "," "," ",0)
       call box3ned( 4.,-4.,-4.," "," "," ",0)
       call box3ned(-4.,-4.,-4.," "," "," ",0)
       call numbned( 4.2,0.,-4.,0.3, rn,0.,-1)
       call numbned(-4.7,0.,-4.,0.3,-rn,0.,-1)
       call numbned(0., 4.2,-4.,0.3, re,0.,-1)
       call numbned(0.,-4.7,-4.,0.3,-re,0.,-1)
       call numbned(0.,0.,-4.2,0.3,-rd,0.,-1)
       call numbned(0.,0., 4.2,0.3, rd,0.,-1)       
c
       call pen(2,0)
       call typstr(5.,19.,0.35,itext,0,60)
       DO j=1,lmod
          if(nwe.eq.0) then
             we = exp(-beta*(xmis(j)-xfbas)/xfbas)
          elseif(nwe.eq.1) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas*xfref)/xfbas)+1.0)
          elseif(nwe.eq.2) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas/xfref)/xfbas)-1.0)
          endif
         ysa = asc*(xlat(j)-olat)
         xsa = bsc*(xlon(j)-olon)
         zsb = dsc*(xdep(j)-odep)
         tsb = tsc*(xsec(j)-osec)
         amisf = xmis(j)
         mj = mod(j,9)+1
         ysc = float(mj-5)*0.05
         xm = amc
         if(amisf.lt.amb) then
           xm = amisf-ama
*           sk = 0.1+0.7*xm/amc
*           call greyton(51,sk)
*           call pen(51,0)
           sk = 0.8*xm/amc
           call hsbcol(51,sk,0.96,0.96)
           call pen(51,0)
         else
*           call pen(26,0)
           call pen(10,0)
         endif
         isy = 16
         ist = 13
         nfill = 0
         if(we.gt.wthres) then
           call pen(2,0)
           ist = 23
           nfill = 1
         endif
*         write(6,*) we,xsa,ysa,zsb
c           
           if(abs(xsa).lt.4. .and. abs(ysa).lt.4. 
     &                       .and. abs(zsb).lt.4.) 
     &         call csymned(xsa,ysa,zsb,3,0.15,isy,nfill)          
           if(abs(tsb).lt.6.) 
     &           call csymbl(x1r+tsb,y1r+ysc,3,0.08,ist)            
       ENDDO
c
       call Hplots(0,0,8,1)
       stop
       end 
        subroutine plotned(x,y,z,kp)
c-------------------------------------------------------------
c       plotned - rearrangement of 3D plotting to suit
c                 North, East, Down type of system
c-------------------------------------------------------------
C
        common /p3d/ xpr,ypr,caxt,saxt,cayt,sayt
c
        xx = xpr+x*caxt+y*sayt
        yy = ypr+x*saxt-y*cayt-z
        call plot(xx,yy,kp)
        return
        end
        subroutine csymned(x,y,z,kp,size,it,nfill)
c-------------------------------------------------------------
c       csymned - centred symbol plotting in
c                 North, East, Down type of system
c-------------------------------------------------------------
C
        common /p3d/ xpr,ypr,caxt,saxt,cayt,sayt
c
        nsid = it-10
        xo = xpr+x*caxt+y*sayt
        yo = ypr+x*saxt-y*cayt-z
        cx = xo
        cy = yo
        call plot(cx,cy,kp)
        if(nsid.EQ.0) nsid = 16
        rpa = size*0.5
        ANG = 6.283185308/float(nsid)
        xv = rpa+cx
        yv = cy
        call plot(xv,yv,3)
        sta=0.0
        do 30 i=1,nsid
          sta = sta+ang
          xq = rpa*cos(sta)
          yq = rpa*sin(sta)
          xv = cx+xq*caxt+yq*sayt
          yv = cy+xq*saxt-yq*cayt
          call plot(xv,yv,2)          
 30     continue
        if(nfill.eq.0) write(8,*) "closepath"
        if(nfill.eq.1) write(8,*) "closepath fill"
        xo = cx
        yo = cy
        return
        end
        subroutine typned(x,y,z,size,M,ang,nc)
c-------------------------------------------------------------
c       symbned - symbol plotting in
c                 North, East, Down type of system
c-------------------------------------------------------------
C
        dimension M(1)
        common /p3d/ xpr,ypr,caxt,saxt,cayt,sayt
c
        xx = xpr+x*caxt+y*sayt
        yy = ypr+x*saxt-y*cayt-z
        call typstr(xx,yy,size,M,ang,nc)
        return
        end
        subroutine  numbned(x,y,z,size,ff,ang,nd)
c-------------------------------------------------------------
c       numbned - number plotting in
c                 North, East, Down type of system
c-------------------------------------------------------------
C
       common /p3d/ xpr,ypr,caxt,saxt,cayt,sayt
c
        xx = xpr+x*caxt+y*sayt
        yy = ypr+x*saxt-y*cayt-z
        call typnum(xx,yy,size,ff,ang,nd)
        return
        end
       subroutine box3ned(xl,yl,zl,labx,laby,labz,idash)
C---------------------------------------------------------------------
C       box3ned: draws perspective view of box using
C                North,East,Down axis system
C
C                xl - length of x-side
C                yl - length of y-side
C                zl - length of z-side
C
C       B.L.N. Kennett   RSES, ANU, July 1987
C
C---------------------------------------------------------------------
C
       common /p3d/ xpr,ypr,caxt,saxt,cayt,sayt
        character*1 labx,laby
c
        call plotned(0. ,0. ,0. ,3)
        call plotned(xl+1.25 ,0. ,0. ,2)
        call typned(xl+1.50 ,0.,0.,0.35,labx,0.0,1)
        call plotned(xl ,0. ,0. ,3)
        call plotned(-xl ,0. ,0. ,2)
        call plotned(xl ,0. ,0. ,3)
        call plotned(xl ,0. ,zl ,2)
c       call dashln(2,5)
        if(idash.ne.1)then
          call dashln(2,0.)
          call plotned(0. ,0. , zl, 2)
        end if
        call plotned(0. ,0. ,0. ,3)
c       call dashln(-1,5)
        call dashln(-1,0.)
        call plotned(0. ,yl+1.25 ,0. ,2)
        call typned(0. ,yl+1.50 ,0.,0.35,laby,0.0,1)
        call plotned(0. , yl,0. ,3)
        call plotned(0. ,-yl,0. ,2)
        call plotned(0. ,yl ,0. ,3)
        call plotned(0. ,yl ,zl ,2)
c       call dashln(2,5)
        if(idash.ne.1)then
          call dashln(2,0.)
          call plotned(0. ,0. ,zl ,2)
        end if
        call typned(0. ,0. ,-zl-0.5,0.35,labyz,0.0,1)
c       call dashln(-1,5)
        call dashln(-1,0.)
        call plotned(0. ,yl ,0. ,3)
        call plotned(xl ,yl ,0. ,2)
        call plotned(xl ,yl ,zl ,2)
        call plotned(0. ,yl ,zl ,2)
        call plotned(xl ,yl ,0. ,3)
        call plotned(xl ,0. ,0. ,2)
        call plotned(xl ,yl ,zl ,3)
        call plotned(xl ,0. ,zl ,2)

        return
        end
