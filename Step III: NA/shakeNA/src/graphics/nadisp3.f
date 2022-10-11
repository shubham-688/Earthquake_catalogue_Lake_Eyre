       character*40 ifile
       character*60 ctext
       dimension itext(10)
       common /p3d/ xpr,ypr,caxt,saxt,cayt,sayt
       equivalence (itext(1),ctext)   
       dgrad = atan(1.0)/45.
c
       dsc = 0.5
       tsc = 2.5
       write(6,*) "Location data file"
       read(5,*) ifile
       open(7,file=ifile)
       write(6,*) "Model sequence/Misfit (0/1)"
       read(5,*) msf
       write(6,*) "msf:",msf
       write(6,*) "event params - lat, long, dep, sec"
       read(5,*)  olat,olon,odep,osec
       write(6,*) olat,olon,odep,osec
       write(6,*) "ranges km/sec - N, E, dep, sec"
       read(5,*)  rn,re,rd,rs
       write(6,*) rn,re,rd,rs
       write(6,*) "misfit bounds"
       read(5,*)  ama,amb
       write(6,*) ama,amb
       amc = amb-ama
c
       sn = 5.0/rn
       se = 5.0/re
       asc = se*111.19
       bsc = sn*111.19*cos(olat*dgrad)
       dsc = 1.5/rd
       tsc = 6.0/rs      
c
       open(8,file="nadisp3.ps")
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
       caxt = cos(15.*dgrad)
       saxt = sin(15.*dgrad)
       cayt = cos(50.*dgrad)
       sayt = sin(50.*dgrad)
       xt = xpr
       yt = ypr+12.
c
       call pen(15,0)
       call box3ned( 5., 5., 1.5,"S","W","Z",0)
       call box3ned(-5., 5., 1.5," "," "," ",0)
       call box3ned( 5.,-5., 1.5," "," "," ",0)
       call box3ned(-5.,-5., 1.5," "," "," ",0)
       call box3ned( 5., 5.,-1.5," "," "," ",0)
       call box3ned(-5., 5.,-1.5," "," "," ",0)
       call box3ned( 5.,-5.,-1.5," "," "," ",0)
       call box3ned(-5.,-5.,-1.5," "," "," ",0)
c
       call pen(2,0)
       read(7,fmt='(a)') ctext
       call typstr(5.,19.,0.35,itext,0,60)
       DO k=0,20
         read(7,*)
         read(7,*) 
         read(7,*) 
         read(7,*)
         IF(msf.eq.0) then
           sk = 0.04*(k+1)
           call hsbcol(k+31,0.5,sk,0.52+0.5*sk)
           call pen(k+31,0)
         ENDIF
         DO j=1,9 
           read(7,*) n,alat,alon,adep,asec,amisf
c           write(6,*) n,alat,alon,adep,asec,amisf
           ysa = asc*(alat-olat)
           xsa = bsc*(alon-olon)
           zsb = dsc*(adep-odep)
           tsb = tsc*(asec-osec)
           ysc = float(j-5)*0.05
           IF(msf.eq.1) then
            xm = amc
             if(amisf.lt.amb) xm = amisf-ama
             sk = 0.8*xm/amc
             call hsbcol(k+31,sk,0.96,0.96)
             call pen(k+31,0)
           ENDIF
           if(abs(xsa).lt.5. .and. abs(ysa).lt.5. 
     &                       .and. abs(zsb).lt.1.5) 
     &         call csymned(xsa,ysa,zsb,3,0.1,13)          
           if(abs(tsb).lt.6.) 
     &           call csymbl(x1r+tsb,y1r+ysc,3,0.08,23)           
         ENDDO
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
        xx = xpr+x*caxt-y*sayt
        yy = ypr-x*saxt-y*cayt-z
        call plot(xx,yy,kp)
        return
        end
        subroutine csymned(x,y,z,kp,size,it)
c-------------------------------------------------------------
c       csymned - centred symbol plotting in
c                 North, East, Down type of system
c-------------------------------------------------------------
C
        common /p3d/ xpr,ypr,caxt,saxt,cayt,sayt
c
        nsid = it-10
        xo = xpr+x*caxt-y*sayt
        yo = ypr-x*saxt-y*cayt-z
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
          yv = cy-xq*saxt+yq*cayt
          call plot(xv,yv,2)          
 30     continue
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
        xx = xpr+x*caxt-y*sayt
        yy = ypr-x*saxt-y*cayt-z     
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
        xx = xpr+x*caxt-y*sayt
        yy = ypr-x*saxt-y*cayt-z
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
        call typned(xl+1.00 ,-1.,0.,0.35,labx,0.0,1)
        call plotned(xl ,0. ,0. ,3)
        call plotned(xl ,0. ,zl ,2)
c       call dashln(2,5)
        if(idash.ne.1)then
          call dashln(2,0.)
          call plotned(0. ,0. ,zl ,2)
          call plotned(0. ,0. ,0. ,2)
        end if
        call plotned(0. ,0. ,0. ,3)
c       call dashln(-1,5)
        call dashln(-1,0.)
        call plotned(0. ,yl+1.25 ,0. ,2)
        call typned(0. ,yl+2.50 ,-1.,0.35,laby,0.0,1)
        call plotned(0. ,yl ,0. ,3)
        call plotned(0. ,yl ,zl ,2)
c       call dashln(2,5)
        if(idash.ne.1)then
          call dashln(2,0.)
          call plotned(0. ,0. ,zl ,2)
        end if
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
        call typned(xl ,yl ,zl,0.35,labyz,0.0,1)
        return
        end
