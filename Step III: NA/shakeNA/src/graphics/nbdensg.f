       character*40 ifile,pfile
       character*60 ctext
       dimension itext(10)
       common /p3d/ xpr,ypr,caxt,saxt,cayt,sayt
       real xlat(1000),xlon(1000),xdep(1000),xsec(1000),xmis(1000)
       real alat,alon,adep,asec,amisf
       real elat,elon,edep,esec
       real vlat,vlon,vdep,vsec
       real slat,slon,sdep,ssec
       real qlat,qlon,qdep,qsec
       real xfmin,xfmax,xfbas
       real vlln(8),vlle(8)
       real vadn(8),vadd(8)
       real vode(8),vodd(8)
       equivalence (itext(1),ctext)   
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
       sn = 5.0/rn
       se = 5.0/re
       dl = 5.0*rd/re
       asc = sn*111.19
       bsc = se*111.19*cos(olat*dgrad)
       scl = cos(olat*dgrad)
       dsc = 0.5*sn
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
*       write(9,*) lmod,flmod
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
c
        write(6,*) "best fit"
        write(6,*) qlat,qlon,qdep,qsec
c
        xfu = xfmax
        if(xfu.gt.5.) xfu =5.
        ama = xfbas
*        amb = (xfu-xfmin)*0.25+xfbas
        amb = xfbas*xfref
        amc = amb-ama
        write(6,*) "misfit bounds"
        write(6,*) ama,amb,amc
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
c
c       weighted model estimate
        elat = 0.
        elon = 0.
        edep = 0.
        esec = 0.
        ewet = 0.
c        l1 = lmod/2
        l1 = 1
        do j=l1,lmod
          if(nwe.eq.0) then
             we = exp(-beta*(xmis(j)-xfbas)/xfbas)
          elseif(nwe.eq.1) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas*xfref)/xfbas)+1)
          elseif(nwe.eq.2) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas/xfref)/xfbas)-1)
          endif
*          write(9,*) beta,j,xmis(j),we
          if(we.gt.wthres) then
           elat = elat + we*xlat(j)
           elon = elon + we*xlon(j)
           edep = edep + we*xdep(j)
           esec = esec + we*xsec(j)
           ewet = ewet + we
          endif
        enddo
        xw = 1.0/ewet
        elat = elat*xw
        elon = elon*xw
        edep = edep*xw
        esec = esec*xw
c
c       weighted variances
        vlat = 0.
        vlon = 0.
        vdep = 0.
        vsec = 0.
        ewet = 0.
        do j=1,lmod
          if(nwe.eq.0) then
             we = exp(-beta*(xmis(j)-xfbas)/xfbas)
          elseif(nwe.eq.1) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas*xfref)/xfbas)+1.0)
          elseif(nwe.eq.2) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas/xfref)/xfbas)-1.0)
          endif
          if(we.gt.wthres) then
           vlat = vlat + we*(xlat(j)-elat)**2
           vlon = vlon + we*(xlon(j)-elon)**2
           vdep = vdep + we*(xdep(j)-edep)**2
           vsec = vsec + we*(xsec(j)-esec)**2
           ewet = ewet + we
          endif
        enddo
        xw = 1.0/ewet
        slat = sqrt(vlat)*xw*111.19
        slon = sqrt(vlon)*xw*111.19*cos(olat*dgrad)
        sdep = sqrt(vdep)*xw
        ssec = sqrt(vsec)*xw
c
        ulat = 0.
        ulon = 0.
        udep = 0.
        usec = 0.
        dlat = 0.
        dlat = 0.
        ddep = 0.
        dsec = 0.
        ulla = 0.
        ullb = 0.
        ullc = 0.
        ulld = 0.
        dlla = 0.
        dllb = 0.
        dllc = 0.
        dlld = 0.
        udla = 0.
        udlb = 0.
        udlo = 0.
        udlp = 0.        
        ddla = 0.
        ddlb = 0.
        ddlo = 0.
        ddlp = 0.
        do j =l1,lmod
          if(nwe.eq.0) then
             we = exp(-beta*(xmis(j)-xfbas)/xfbas)
          elseif(nwe.eq.1) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas*xfref)/xfbas)+1.0)
          elseif(nwe.eq.2) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas/xfref)/xfbas)-1.0)
          endif
          if(we.gt.wthres) then
            sx = (xlon(j)-qlon)*scl*111.19
            sy = (xlat(j)-qlat)*111.19
            sz = xdep(j)-qdep
            ss = sqrt(sx*sx+sy*sy)
            sr = sqrt(sx*sx+sz*sz)
            st = sqrt(sy*sy+sz*sz)
            su = sqrt(sx*sx+sy*sy+sz*sz)
*            write(9,fmt='(i4,4f12.5)') j,we,sx,sy,sz
            if(qlat-xlat(j) .gt.dlat) dlat =  qlat-xlat(j)     
            if(xlat(j)-qlat .gt.ulat) ulat =  xlat(j)-qlat
            if(qlon-xlon(j) .gt.dlon) dlon =  qlon-xlon(j)     
            if(xlon(j)-qlon .gt.ulon) ulon =  xlon(j)-qlon
            if(qdep-xdep(j) .gt.ddep) ddep =  qdep-xdep(j)     
            if(xdep(j)-qdep .gt.udep) udep =  xdep(j)-qdep
            if(qsec-xsec(j) .gt.dsec) dsec =  qsec-xsec(j)     
            if(xsec(j)-qsec .gt.usec) usec =  xsec(j)-qsec
            if(sx.gt.0. .and. sy.gt.0.) then
               if(ss.gt.ulla) then
                 ulla  = ss
                 ullae = sx
                 ullan = sy              
               endif
            elseif(sx.lt.0. and. sy.lt.0.) then
               if(ss.gt.dlla) then
                 dlla  = ss 
                 dllae = sx
                 dllan = sy              
               endif           
            elseif(sx.lt.0. and. sy.gt.0.) then
               if(ss.gt.ullb) then
                 ullb = ss
                 ullbe = sx
                 ullbn = sy              
               endif          
            elseif(sx.gt.0. and. sy.lt.0.) then
               if(ss.gt.dllb) then
                 dllb = ss 
                 dllbe = sx
                 dllbn = sy
               endif        
            endif
            if(sx.gt.0. .and. sz.gt.0.) then
               if(sr.gt.udlo) udlo = sr
            elseif(sx.lt.0. and. sz.lt.0.) then
               if(sr.gt.ddlo) ddlo = sr            
            elseif(sx.lt.0. and. sz.gt.0.) then
               if(sr.gt.udlp) udlp = sr         
            elseif(sx.gt.0. and. sz.lt.0.) then
               if(sr.gt.ddlp) ddlp = sr         
            endif
            if(sy.gt.0. .and. sz.gt.0.) then
               if(st.gt.udla) udla = st
            elseif(sy.lt.0. and. sz.lt.0.) then
               if(st.gt.ddla) ddla = st            
            elseif(sy.lt.0. and. sz.gt.0.) then
               if(st.gt.udlb) udlb = st         
            elseif(sy.gt.0. and. sz.lt.0.) then
               if(st.gt.ddlb) ddlb = st         
            endif
          endif
        enddo
        dlat = dlat*111.19
        ulat = ulat*111.19
        dlon = dlon*111.19*scl
        ulon = ulon*111.19*scl
c            
        write(6,*) "spreads"
        write(6,*) dlat,ulat
        write(6,*) dlon,ulon
        write(6,*) ddep,udep
        write(6,*) dsec,usec
        write(6,*) dlla,ulla
        write(6,*) dllb,ullb
        write(6,*) ddla,udla
        write(6,*) ddlb,udlb
        write(6,*) ddlo,udlo
        write(6,*) ddlp,udlp       
c
        write(6,*) "ensemble location"
        write(6,*) elat,elon,edep,esec
        write(6,*) "ensemble errors"
        write(6,*) slat,slon,sdep,ssec
c

cc
       open(8,file=pfile)
       call Hplots(1,0,8,0)
       call typset(0.0,0.0)
       call zpick(3,0,is)
       x1r =  8.0
       y1r = 18.0
       x2r =  8.0
       y2r = 12.0
       x3r =  8.0
       y3r =  7.0-0.85-d1
       x4r = 13.0+0.90+d1
       y4r = 12.0
c
       call pen(2,0)
       call plot(x1r-5.,y1r-0.25,3)
       call edgert(10.,0.5,0)
       call typnum(x1r-5.9,y1r-0.1,0.30,-rs,0.,-1)
       call typnum(x1r+5.1,y1r-0.1,0.30, rs,0.,-1)
       call typstr(x1r+5.7,y1r-0.1,0.30,"[s] Origin Time",0.,15)
       call pen(15,0)
       call plot(x1r-3.0,y1r-0.25,3)
       call plot(x1r-3.0,y1r+0.25,2)
       call plot(x1r    ,y1r-0.25,3)
       call plot(x1r    ,y1r+0.25,2)
       call plot(x1r+3.0,y1r-0.25,3)
       call plot(x1r+3.0,y1r+0.25,2)
c
       call pen(2,0)
       call plot(x2r-5.0,y2r-5.0,3)
       call edgert(10.,10.,0)
       call typnum(x2r-5.7,y2r-0.1,0.30,-re,0.,-1)
       call typnum(x2r+5.1,y2r-0.1,0.30, re,0.,-1)
       call typnum(x2r    ,y2r-5.4,0.30,-rn,0.,-1)
       call typnum(x2r    ,y2r+5.1,0.30, rn,0.,-1)
       call typstr(x2r-4.95,y2r+0.05,0.22,"W",0.,1)
       call typstr(x2r+4.73,y2r+0.05,0.22,"E",0.,1)
       call typstr(x2r+0.05,y2r-4.95,0.22,"S",0.,1)
       call typstr(x2r+0.05,y2r+4.73,0.22,"N",0.,1)
       call typstr(x2r-6.0,y2r+1.0,0.30,"[km]",0.,4)
       call pen(15,0)
       call plot(x2r-2.5,y2r-5.0,3)
       call plot(x2r-2.5,y2r+5.0,2)
       call plot(x2r    ,y2r-5.0,3)
       call plot(x2r    ,y2r+5.0,2)
       call plot(x2r+2.5,y2r-5.0,3)
       call plot(x2r+2.5,y2r+5.0,2)
       call plot(x2r-5.0,y2r-2.5,3)
       call plot(x2r+5.0,y2r-2.5,2)
       call plot(x2r-5.0,y2r    ,3)
       call plot(x2r+5.0,y2r    ,2)
       call plot(x2r-5.0,y2r+2.5,3)
       call plot(x2r+5.0,y2r+2.5,2)
c
       d1 = 0.5*dl
       call pen(2,0)
       call plot(x3r-5.0,y3r-d1,3)
       call edgert(10.,dl,0)
       call typstr(x3r-4.95,y3r+0.05,0.22,"W",0.,1)
       call typstr(x3r+4.73,y3r+0.05,0.22,"E",0.,1)
       call typnum(x3r-5.0,y3r-d1-0.4,0.30,-rd,0.,-1)
       call typnum(x3r-5.0,y3r+d1+0.1,0.30, rd,0.,-1)
       call typstr(x3r+5.1,y3r+0.1,0.30,"Depth",0.,5)
       call typstr(x3r+5.1,y3r-0.3,0.30," [km]",0.,5)
       call pen(15,0)
       call plot(x3r-2.5,y3r-d1,3)
       call plot(x3r-2.5,y3r+d1,2)
       call plot(x3r    ,y3r-d1,3)
       call plot(x3r    ,y3r+d1,2)
       call plot(x3r+2.5,y3r-d1,3)
       call plot(x3r+2.5,y3r+d1,2)
       call plot(x3r-5.0,y3r,3)
       call plot(x3r+5.0,y3r,2)
c
       call pen(2,0)
       call plot(x4r-d1,y4r-5.0,3)
       call edgert(dl,10.,0)
       call typnum(x4r-d1-0.7,y4r-5.1,0.30,-rd,0.,-1)
       call typnum(x4r+d1+0.1,y4r-5.1,0.30, rd,0.,-1)
       call typstr(x4r+0.05,y4r-4.95,0.22,"S",0.,1)
       call typstr(x4r+0.05,y4r+4.73,0.22,"N",0.,1)
       call pen(15,0)
       call plot(x4r-d1,y4r-2.5,3)
       call plot(x4r+d1,y4r-2.5,2)
       call plot(x4r-d1,y4r    ,3)
       call plot(x4r+d1,y4r    ,2)
       call plot(x4r-d1,y4r+2.5,3)
       call plot(x4r+d1,y4r+2.5,2)
       call plot(x4r,y4r-5.0,3)
       call plot(x4r,y4r+5.0,2)
c
       call pen(2,0)
       call zpick(6,0,is)
       call typstr(2.08,18.5,0.33,itext,0,60)
c      
c  plot  consistency patch
        call pen(2,0)
        ysa = asc*(qlat-olat)
        xsa = bsc*(qlon-olon)
        ysb = dsc*(qdep-odep)
        xsb = tsc*(qsec-osec)
        write(6,*) "opt",ysa,xsa,ysb,xsb
c
        vlle(1) =  ulon
        vlle(2) =  ulla*sq2
        vlle(3) = 0.
        vlle(4) = -ullb*sq2
        vlle(5) = -dlon
        vlle(6) = -dlla*sq2
        vlle(7) = 0.
        vlle(8) =  dllb*sq2
        vlln(1) = 0.
        vlln(2) =  ulla*sq2
        vlln(3) =  ulat
        vlln(4) =  ullb*sq2
        vlln(5) = 0.
        vlln(6) = -dlla*sq2
        vlln(7) = -dlat
        vlln(8) = -dllb*sq2
c
        vadn(1) =  ulat
        vadn(2) =  udla*sq2
        vadn(3) = 0.
        vadn(4) = -udlb*sq2
        vadn(5) = -dlat
        vadn(6) = -ddla*sq2
        vadn(7) = 0.
        vadn(8) =  ddlb*sq2
        vadd(1) = 0.
        vadd(2) =  udla*sq2
        vadd(3) =  udep
        vadd(4) =  udlb*sq2
        vadd(5) = 0.
        vadd(6) = -ddla*sq2
        vadd(7) = -ddep
        vadd(8) = -ddlb*sq2
c
        vode(1) =  ulon
        vode(2) =  udlo*sq2
        vode(3) = 0.
        vode(4) = -udlp*sq2
        vode(5) = -dlon
        vode(6) = -ddlo*sq2
        vode(7) = 0.
        vode(8) =  ddlp*sq2
        vodd(1) = 0.
        vodd(2) =  udlo*sq2
        vodd(3) =  udep
        vodd(4) =  udlp*sq2
        vodd(5) = 0.
        vodd(6) = -ddlo*sq2
        vodd(7) = -ddep
        vodd(8) = -ddlp*sq2       
c
        pdla = dlat*sn
        pula = ulat*sn
        pdlo = dlon*se
        pulo = ulon*se
        pdde = ddep*dsc
        pude = udep*dsc
        pdse = dsec*tsc
        puse = usec*tsc

        call pen(29,0)
        call plot(x2r+xsa+1.05*vlle(1)*se,y2r+ysa+1.05*vlln(1)*sn,3)
        do k=1,8
          call plot(x2r+xsa+1.05*vlle(k)*se,y2r+ysa+1.05*vlln(k)*sn,2)
        enddo        
        call plot(x2r+xsa+1.05*vlle(1)*se,y2r+ysa+1.05*vlln(1)*sn,2)
        write(8,*) "cf"   
c     
*        call pen(13,0)
*        call plot(x2r+xsa,y2r+ysa,3)
*        call arc(x2r+xsa,y2r+ysa,ulla*sn,0.,90.)
*        call arc(x2r+xsa,y2r+ysa,ullb*sn,90.,180.)
*        call arc(x2r+xsa,y2r+ysa,dlla*sn,180.,270.)
*        call arc(x2r+xsa,y2r+ysa,dllb*sn,270.,360.)
*        write(8,*) "cf"
c
        call pen(29,0)
        call plot(x3r+xsa+vode(1)*se,y3r-ysb-vodd(1)*sd,3)
        do k=1,8
         call plot(x3r+xsa+vode(k)*se,y3r-ysb-vodd(k)*sd,2)
        enddo        
        call plot(x3r+xsa+vode(1)*se,y3r-ysb-vodd(1)*sd,2)
        write(8,*) "cf"        
c        
        call pen(29,0)
        call plot(x4r+ysb+vadd(1)*sd,y4r+ysa+vadn(1)*sn,3)
        do k=1,8
         call plot(x4r+ysb+vadd(k)*sd,y4r+ysa+vadn(k)*sn,2)
        enddo        
        call plot(x4r+ysb+vadd(1)*sd,y4r+ysa+vadn(1)*sn,2)
        write(8,*) "cf" 
c
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
         ysb = dsc*(xdep(j)-odep)
         xsb = tsc*(xsec(j)-osec)
         amisf = xmis(j)
         mj = mod(j,9)+1
         ysc = float(mj-5)*0.05
         xm = amc
         if(amisf.lt.amb) then
           xm = amisf-ama
             sk = 0.1+0.6*xm/amc
             call greyton(51,sk)
             call pen(51,0)
         else
           call pen(26,0)
         endif
         isy = 13
         if(we.gt.wthres) then
           call pen(10,0)
           isy = 24
         endif
c           
         if(abs(xsa).lt.5. .and. abs(ysa).lt.5.) 
     &           call csymbl(x2r+xsa,y2r+ysa,3,0.12,isy)
         if(abs(xsa).lt.5. .and. abs(ysb).lt.d1) 
     &           call csymbl(x3r+xsa,y3r-ysb,3,0.12,isy)
         if(abs(ysa).lt.5. .and. abs(ysb).lt.d1) 
     &           call csymbl(x4r+ysb,y4r+ysa,3,0.12,isy)
         if(abs(xsb).lt.6.) 
     &           call csymbl(x1r+xsb,y1r+ysc,3,0.12,isy)           
       ENDDO
c 
c  plot weighted model
        call pen(27,0)
        ysa = asc*(elat-olat)
        xsa = bsc*(elon-olon)
        ysb = dsc*(edep-odep)
        xsb = tsc*(esec-osec)
        write(6,*) "we",ysa,xsa,ysb,xsb
        ysc = -0.35
        pla = 2.*slat*sn
        plo = 2.*slon*se
        pde = 2.*sdep*dsc
        pse = 2.*ssec*tsc
        if(abs(xsa).lt.5. .and. abs(ysa).lt.5.) 
     &           call csymbl(x2r+xsa,y2r+ysa,3,0.15,13)
        if(abs(xsa).lt.5. .and. abs(ysb).lt.d1) 
     &           call csymbl(x3r+xsa,y3r-ysb,3,0.15,13)
        if(abs(ysa).lt.5. .and. abs(ysb).lt.d1) 
     &           call csymbl(x4r+ysb,y4r+ysa,3,0.15,13)
        if(abs(xsb).lt.6.) 
     &           call csymbl(x1r+xsb,y1r+ysc,3,0.15,13)

        call plot(x2r+xsa,y2r+ysa+pla,3) 
        call plot(x2r+xsa,y2r+ysa-pla,2) 
        call plot(x2r+xsa+plo,y2r+ysa,3) 
        call plot(x2r+xsa-plo,y2r+ysa,2)
        call plot(x3r+xsa+plo,y3r-ysb,3)
        call plot(x3r+xsa-plo,y3r-ysb,2)
        call plot(x3r+xsa,y3r-ysb-pde,3)
        call plot(x3r+xsa,y3r-ysb+pde,2)        
        call plot(x4r+ysb+pde,y4r+ysa,3)
        call plot(x4r+ysb-pde,y4r+ysa,2)
        call plot(x4r+ysb,y4r+ysa-pla,3)
        call plot(x4r+ysb,y4r+ysa+pla,2)
        call plot(x1r+xsb+pse,y1r+ysc,3)
        call plot(x1r+xsb-pse,y1r+ysc,2) 
c
c  plot "optimum model"
        call pen(2,0)
        ysa = asc*(qlat-olat)
        xsa = bsc*(qlon-olon)
        ysb = dsc*(qdep-odep)
        xsb = tsc*(qsec-osec)
        ysc = -0.40
        call pen(29,0)
        if(abs(xsa).lt.5. .and. abs(ysa).lt.5.) 
     &           call csymbl(x2r+xsa,y2r+ysa,3,0.125,23)
        if(abs(xsa).lt.5. .and. abs(ysb).lt.d1) 
     &           call csymbl(x3r+xsa,y3r-ysb,3,0.125,23)
        if(abs(ysa).lt.5. .and. abs(ysb).lt.d1) 
     &           call csymbl(x4r+ysb,y4r+ysa,3,0.125,23)
        if(abs(xsb).lt.6.) 
     &           call csymbl(x1r+xsb,y1r+ysc,3,0.125,23)
        call pen(1,0)
        if(abs(xsa).lt.5. .and. abs(ysa).lt.5.) 
     &           call csymbl(x2r+xsa,y2r+ysa,3,0.125,13)
        if(abs(xsa).lt.5. .and. abs(ysb).lt.d1) 
     &           call csymbl(x3r+xsa,y3r-ysb,3,0.125,13)
        if(abs(ysa).lt.5. .and. abs(ysb).lt.d1) 
     &           call csymbl(x4r+ysb,y4r+ysa,3,0.125,13)
        if(abs(xsb).lt.6.) 
     &           call csymbl(x1r+xsb,y1r+ysc,3,0.125,13)
        call pen(1,0)
        call plot(x2r+xsa,y2r+ysa-pdla,3) 
        call plot(x2r+xsa,y2r+ysa+pula,2)
        call plot(x2r+xsa-pdlo,y2r+ysa,3) 
        call plot(x2r+xsa+pulo,y2r+ysa,2)
        call pen(1,0)
        call plot(x3r+xsa-pdlo,y3r-ysb,3)
        call plot(x3r+xsa+pulo,y3r-ysb,2)
        call plot(x3r+xsa,y3r-ysb+pdde,3)
        call plot(x3r+xsa,y3r-ysb-pude,2)
        call pen(1,0)                
        call plot(x4r+ysb-pdde,y4r+ysa,3)
        call plot(x4r+ysb+pude,y4r+ysa,2)
        call plot(x4r+ysb,y4r+ysa-pdla,3)
        call plot(x4r+ysb,y4r+ysa+pula,2)
        call pen(1,0)
        call plot(x1r+xsb+puse,y1r+ysc,3)
        call plot(x1r+xsb-pdse,y1r+ysc,2) 
c                                
c
c   plot weighting function
       call typstr(25.00,15.6,0.33,"Weighting",0,9)
       call zpick(3,0,is)
       call typstr(23.00,5.0,0.30,"Model",0,5)
       call pen(2,0)
       call zpick(3,0,is)
       xo = 18.
       yo = 6.
       psc = 10.0*flmod
       csc = 25.0/xfu
       call plot(xo,yo,3)
       call plot(xo+10.0,yo,2)
       call plot(xo,yo,3)
       do j=1,lmod,50
         xj = float(j-1)
         call plot(xo+xj*psc,yo,3)
         call plot(xo+xj*psc,yo+0.3,2)
         call typnum(xo+xj*psc-0.1,yo-0.5,0.25,xj,0.0,-3)
       enddo
       call plot(xo,yo,3) 
       call plot(xo,yo+2.*csc,2)
       call plot(xo+10.0,yo,3) 
       call plot(xo+10.0,yo+2.*csc,2)
       do c=0.,1.0,0.25
         call plot(xo    ,yo+c*2.*csc,3)
         call plot(xo+0.3,yo+c*2.*csc,2)
         call typnum(xo-0.9,yo+c*2*csc-0.1,0.28,c,0.0,2)
         call plot(xo+10.,yo+c*2.*csc,3)
         call plot(xo+9.7,yo+c*2.*csc,2)
         call typnum(xo+10.1,yo+c*2*csc-0.1,0.28,c,0.0,2)
       enddo         
       do j=1,lmod
          if(nwe.eq.0) then
             we = exp(-beta*(xmis(j)-xfbas)/xfbas)
          elseif(nwe.eq.1) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas*xfref)/xfbas)+1.0)
          elseif(nwe.eq.2) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas/xfref)/xfbas)-1.0)
          endif
          wf = we
          if(wf.gt.xfmax) wf=wfmax
          isy = 13
         amisf = xmis(j)
         if(amisf.lt.amb) then
           xm = amisf-ama
             sk = 0.1+0.6*xm/amc
             call greyton(51,sk)
             call pen(51,0)
         else
           call pen(26,0)
         endif
          if(wf.gt.wthres) then
            call pen(10,0)
            isy = 24
          endif
         call csymbl(xo+float(j)*psc,yo+wf*csc*2,3,0.08,isy)
       enddo        
c       
       call Hplots(0,0,8,1)
       stop
       end 
