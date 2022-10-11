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
       real vlln(12),vlle(12),vlls(12)
       real vadn(12),vadd(12),vads(12)
       real vode(12),vodd(12),vods(12)
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
        xfu = xfmax
        if(xfu.gt.5.) xfu =5.
        ama = xfbas
*        amb = (xfu-xfmin)*0.25+xfbas
        amb = xfbas*xfref
        amc = amb-ama
        write(6,*) "misfit bounds"
        write(6,*) ama,amb,amc
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
c       weighted model estimate and variances
        elat = 0.
        elon = 0.
        edep = 0.
        esec = 0.
        ewet = 0.
        xbm = 0.
         do j=1,lmod
          if(nwe.eq.0) then
             we = exp(-beta*(xmis(j)-xfbas)/xfbas)
          elseif(nwe.eq.1) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas*xfref)/xfbas)+1)
          elseif(nwe.eq.2) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas/xfref)/xfbas)-1)
          endif
          if(we.gt.wthres) then
            if(xmis(j).gt.xbm) xbm = xmis(j)
            elat = elat + we*xlat(j)
            elon = elon + we*xlon(j)
            edep = edep + we*xdep(j)
            esec = esec + we*xsec(j)
            ewet = ewet + we
          endif
        enddo
        write(6,*) "misfit base", xbm
        xw = 1.0/ewet
        elat = elat*xw
        elon = elon*xw
        edep = edep*xw
        esec = esec*xw
c
        vlat = 0.
        vlon = 0.
        vdep = 0.
        vsec = 0.
        ewet = 0.
         do j=1,lmod
          if(nwe.eq.0) then
             we = exp(-beta*(xmis(j)-xfbas)/xfbas)
          elseif(nwe.eq.1) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas*xfref)/xfbas)+1)
          elseif(nwe.eq.2) then
             we = 1.0/(exp(beta*(xmis(j)-xfbas/xfref)/xfbas)-1)
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
        slat = sqrt(vlat*xw)*111.19
        slon = sqrt(vlon*xw)*111.19*cos(olat*dgrad)
        sdep = sqrt(vdep*xw)
        ssec = sqrt(vsec*xw)
c
c       measure of fit?

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
        do j =1,lmod
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
            sx = (xlon(j)-elon)*scl*111.19
            sy = (xlat(j)-elat)*111.19
            sz = xdep(j)-edep
            ss = sqrt(sx*sx+sy*sy)
            sr = sqrt(sx*sx+sz*sz)
            st = sqrt(sy*sy+sz*sz)
            su = sqrt(sx*sx+sy*sy+sz*sz)
            rll = atan2(sy,sx)
            rad = atan2(sz,sy)
            rod = atan2(sz,sx)
            write(9,*) j,sx,sy,sz,rll/dgrad,rad/dgrad,rod/dgrad
            do k=1,12
              ta = float(k-7.20)*30.0*dgrad
              tb = ta+42.0*dgrad
              if(rll.ge.ta .and. rll.lt.tb) then
                if(ss. gt. vlls(k)) then
                  vlls(k) = ss
                  vlle(k) = xlon(j)
                  vlln(k) = xlat(j)
                endif
              endif
               if(rod.ge.ta .and. rod.lt.tb) then
                if(sr. gt. vods(k)) then
                  vods(k) = sr
                  vode(k) = xlon(j)
                  vodd(k) = xdep(j)
                endif
              endif
              if(rad.ge.ta .and. rad.lt.tb) then
                if(st. gt. vads(k)) then
                  vads(k) = st
                  vadn(k) = xlat(j)
                  vadd(k) = xdep(j)
                endif
              endif             
            enddo
          endif
        enddo
        do k=1,12
          kl = k-1
          if(kl.eq.0) kl = 12
          ang = float(k-6.5)*30.0*dgrad
          cang = cos(ang)
          sang = sin(ang)
          if(vlln(k).eq.0.25*cang .and. vlle(k).eq.0.25*sang) then
             vlln(k) = vlln(kl)
             vlle(k) = vlle(kl)
          endif
          if(vadn(k).eq.0.25*cang .and. vadd(k).eq.0.25*sang) then
             vadn(k) = vadn(kl)
             vadd(k) = vadd(kl)
          endif
          if(vode(k).eq.0.25*cang .and. vodd(k).eq.0.25*sang) then
             vode(k) = vode(kl)
             vodd(k) = vodd(kl)
          endif
        enddo

        dlat = dlat*111.19
        ulat = ulat*111.19
        dlon = dlon*111.19*scl
        ulon = ulon*111.19*scl
c
        write(6,*) "best fit location"
        write(6,*) qlat,qlon,qdep,qsec
        write(6,*) "ensemble location"
        write(6,*) elat,elon,edep,esec
        write(6,*) "ensemble errors"
        write(6,*) slat,slon,sdep,ssec
        write(6,*) "spread - time"
        write(6,*) dsec,usec
        write(6,*) "lat, long" 
        do i=1,12
         write(6,*) vlln(i),vlle(i)
        enddo
        write(6,*) "long,depth"
        do i =1,12
          write(6,*) vode(i),vodd(i)
        enddo
        write(6,*) "lat,depth"
        do i =1,12
          write(6,*) vadn(i),vadd(i)
        enddo
c

c
       stop
       end 
