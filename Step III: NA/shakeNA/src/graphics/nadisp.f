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
       sn = 5.0/rn
       se = 5.0/re
       dl = 5.0*rd/re
       asc = sn*111.19
       bsc = se*111.19*cos(olat*dgrad)
       dsc = 0.5*sn
       tsc = 5.0/rs
       d1 = 0.5*dl
c
       open(8,file="nadisp.ps")
       call Hplots(1,0,8,0)
       call typset(0.0,0.0)
       call zpick(3,0,is)
       x1r =  9.0
       y1r = 18.0
       x2r =  9.0
       y2r = 12.0
       x3r =  9.0
       y3r =  7.0-0.85-d1
       x4r = 14.0+0.90+d1
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
       call typnum(x2r-5.9,y2r-0.1,0.30,-re,0.,-1)
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
       call typnum(x4r-d1-0.6,y4r-5.1,0.30,-rd,0.,-1)
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
       read(7,fmt='(a)') ctext
       call zpick(6,0,is)
       call typstr(4.08,18.5,0.33,itext,0,60)
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
           ysb = dsc*(adep-odep)
           xsb = tsc*(asec-osec)
           ysc = float(j-5)*0.05
           IF(msf.eq.1) then
            xm = amc
             if(amisf.lt.amb) xm = amisf-ama
             sk = 0.8*xm/amc
             call hsbcol(k+31,sk,0.96,0.96)
             call pen(k+31,0)
           ENDIF
           if(abs(xsa).lt.5. .and. abs(ysa).lt.5.) 
     &           call csymbl(x2r+xsa,y2r+ysa,3,0.08,23)
           if(abs(xsa).lt.5. .and. abs(ysb).lt.d1) 
     &           call csymbl(x3r+xsa,y3r-ysb,3,0.08,23)
           if(abs(ysa).lt.5. .and. abs(ysb).lt.d1) 
     &           call csymbl(x4r+ysb,y4r+ysa,3,0.08,23)
           if(abs(xsb).lt.6.) 
     &           call csymbl(x1r+xsb,y1r+ysc,3,0.08,23)           
         ENDDO
       ENDDO
c      
       call Hplots(0,0,8,1)
       stop
       end 
