       character*24 modnam
       character*40 ctext
       dimension itext(10)
       common /p3d/ xpr,ypr,caxt,saxt,cayt,sayt
       equivalence (itext(1),ctext)   
       dgrad = atan(1.0)/45.
c
       open(8,file="evdk.ps")
       call Hplots(1,0,8,1)
       call typset(0.0,0.0)
       call zpick(3,0,is)
       call plot(24.0,4.0,3)
c      call edgert(8.0,-5.0)
       read(5,*) ctext
       call typstr(24.5,7.5,0.55,itext,0,40)
       read(5,*) ctext
       call typstr(24.5,6.,0.5,itext,0,40)
       read(5,*) ctext
       call typstr(24.5,5.0,0.5,itext,0,40)
c      
       xpr = 14.0
       ypr = 14.0
       caxt = cos(15.*dgrad)
       saxt = sin(15.*dgrad)
       cayt = cos(50.*dgrad)
       sayt = sin(50.*dgrad)
       xt = xpr
       yt = ypr+12.
       call plot(xt-11.0,yt    ,3)
       call plot(xt+11.0,yt    ,2)
       call plot(xt+11.0,yt+0.5,2)
       call plot(xt-11.0,yt+0.5,2)
       call plot(xt-11.0,yt    ,2)
       call typstr(xt+12.,yt,0.6,'origin time [s]',0.0,15)
c
       call typned(13. , 0. , 0. ,0.6,'latitude [km]',0.0,13)
       call typned( 1. ,-12.0, 0. ,0.6,'longitude [km]',0.0,14)
       call typned(11.5, 0. ,0.1 ,0.5,'N',0.0,1)
       call typned(-11.5, 0. ,0.1,0.5,'S',0.0,1)
       call typned( 0. ,11.7 ,0. ,0.5,'E',0.0,1)
       call typned( 0. ,-11.3,0., 0.5,'W',0.0,1)
       call typned(-10.1,-0.3, 0.,0.475,'-20',0.0,4)
       call typned( -5.1,-0.3, 0.,0.475,'-10',0.0,3)
       call typned(  5.1,-0.3, 0.,0.475,'10',0.0,2)
       call typned( 10.1,-0.3, 0.,0.475,'20',0.0,3)
       call typned( -1.1,-10., 0.,0.475,'-20',0.0,4)
       call typned( -1.1, -5., 0.,0.475,'-10',0.0,4)
       call typned( -1.1,  5., 0.,0.475,' 10',0.0,4)
       call typned( -1.1, 10., 0.,0.475,' 20',0.0,4)
       call typned( -1.1, 0.,-10.,0.475,'-20',0.0,4)
       call typned( -1.1, 0., -5.,0.475,'-10',0.0,3)
       call typned( -1.1, 0.,  5.,0.475,' 10',0.0,3)
       call typned( -1.1, 0., 10.,0.475,' 20',0.0,4)
       call typned( 0.5 , 0. ,11.2 ,0.6,'depth [km]',0.0,10)
       call plotned( 0. , 0. , 0. ,3)
       call plotned(-11., 0. , 0. ,2)
       call plotned(-10., 0. , 0. ,2)
       call plotned(-10., 0.3, 0. ,2)
       call plotned(-10., 0. , 0. ,2)
       call plotned( -5., 0. , 0. ,2)
       call plotned( -5., 0.3, 0. ,2)
       call plotned( -5., 0. , 0. ,2)
       call plotned(  5., 0. , 0. ,2)
       call plotned(  5., 0.3, 0. ,2)      
       call plotned(  5., 0. , 0. ,2)
       call plotned( 10., 0. , 0. ,2)
       call plotned( 10., 0.3, 0. ,2)
       call plotned( 10., 0. , 0. ,2)
       call plotned( 11., 0. , 0. ,2)
       call plotned( 0. , 0. , 0. ,3)
       call plotned( 0. ,-11., 0. ,2)
       call plotned( 0. ,-10., 0. ,2)
       call plotned( 0.3,-10., 0. ,2)
       call plotned( 0. ,-10., 0. ,2)
       call plotned( 0. , -5., 0. ,2)
       call plotned( 0.3, -5., 0. ,2)
       call plotned( 0. , -5., 0. ,2)
       call plotned( 0. ,  5., 0. ,2)
       call plotned( 0.3,  5., 0. ,2)
       call plotned( 0. ,  5., 0. ,2)
       call plotned( 0. , 10., 0. ,2)
       call plotned( 0.3, 10., 0. ,2)
       call plotned( 0. , 10., 0. ,2)
       call plotned( 0. , 11., 0. ,2)
       call plotned( 0. , 0. , 0. ,3)
       call plotned( 0. , 0. ,-11.,2) 
       call plotned( 0. , 0. ,-10.,2)
       call plotned( 0.3, 0. ,-10.,2)
       call plotned( 0. , 0. ,-10.,2)
       call plotned( 0. , 0. , -5.,2)
       call plotned( 0.3, 0. , -5.,2)
       call plotned( 0. , 0. , -5.,2)
       call plotned( 0. , 0. ,  5.,2)
       call plotned( 0.3, 0. ,  5.,2)
       call plotned( 0. , 0. ,  5.,2)
       call plotned( 0. , 0. , 10.,2)
       call plotned( 0.3, 0. , 10.,2)
       call plotned( 0. , 0. , 10.,2)
       call plotned( 0. , 0. , 0. ,3)
c	
       call pen(1,0)
       call csymned(0.,0.,0.,3,0.6,15)
       call pen(1,0)
       call plot(xt-10.,yt    ,3)
       call plot(xt-10.,yt+0.5,2)
       call plot(xt- 5.,yt    ,3)
       call plot(xt- 5.,yt+0.5,2)
       call plot(xt+ 5.,yt    ,3)
       call plot(xt+ 5.,yt+0.5,2)
       call plot(xt+10.,yt    ,3)
       call plot(xt+10.,yt+0.5,2)
       call typstr(xt-10.3,yt+0.6,0.475,'-2',0.0,2)
       call typstr(xt-5.3,yt+0.6,0.475,'-1',0.0,2)
       call typstr(xt-0.2,yt+0.6,0.475,'0',0.0,1)
       call typstr(xt+4.8,yt+0.6,0.475,'1',0.0,1)
       call typstr(xt+9.8,yt+0.6,0.475,'2',0.0,1)
       call plot(xt-0.02,yt    ,3)
       call plot(xt-0.02,yt+0.5,2)
       call plot(xt+0.02,yt+0.5,2)
       call plot(xt+0.02,yt    ,2)
       call plotned(0.0,0.0,0.0,3)
       call zpick(1,0,is)
c
       dsc = 0.5
       tsc = 5.0
       write(6,*) "event params - lat, long, dep, sec"
       read(5,*)  olat,olon,odep,osec
       write(6,*) olat,olon,odep,osec
       write(6,*) "nsol"
       asc = 111.19*0.5
       bsc = 111.19*0.5*cos(olat*dgrad)
       read(5,*) nsol
       do 100 k=1,nsol
         read(5,*)  modnam
         write(6,*) modnam
         write(6,*) "pen colour, # sides for symbol"
         read(5,*)  ipen,isym
         write(6,*) ipen,isym
         isym = isym+10
         xk = k
         write(6,*) "event est. - lat,long,depth, sec"
         read(5,*)  clat4,clon4,cdep4,csec4
         write(6,*) clat4,clon4,cdep4,csec4
         call pen(ipen,0)
         xsa = asc*(clat4-olat)
         ysa = bsc*(clon4-olon)
         ysb = dsc*(cdep4-odep)
         xsb = tsc*(csec4-osec)
         call plotned(0.,ysa,0.,3)
         call plotned(xsa,ysa,0.,2)
         call plotned(xsa,0.,0.,2)
         call csymned(xsa,ysa, 0.,3,0.6,isym)
         call plotned(xsa,ysa, 0.,3)
         call csymned(xsa,ysa,ysb,2,0.6,isym)
         call pen(ipen,2)
         call csymned(xsa,ysa,ysb,3,0.4,isym)
         call csymned(xsa,ysa,ysb,3,0.2,isym)
         call pen(2,0)
         call numbned(xsa+0.3,ysa-0.3,ysb+0.1,0.4,xk,0.0,-1)
         call pen(ipen,2)
         call plot(xt+xsb,yt,3)
         call plot(xt+xsb,yt+0.5,2)
         call plot(xt+xsb+0.1,yt+0.25,2)
         call plot(xt+xsb,yt,2)
         call plot(xt+xsb-0.1,yt+0.25,2)
         call plot(xt+xsb,yt+0.5,2)
         call pen(2,0)
         call typnum(xt+xsb-0.12,yt-0.55,0.4,xk,0.0,-1)
         call plot(xpr,ypr,3)
100    continue
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
        nsid = it-9
        nfill = 0
        xo = xpr+x*caxt-y*sayt
        yo = ypr-x*saxt-y*cayt-z
        cx = xo
        cy = yo
        call plot(cx,cy,kp)
        if(nsid.gt.10) then
           nsid = nsid-10
           nfill = 1
        endif
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
        write(8,*) "closepath" 
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
