       character*40 ifile,pfile
       character*60 ctext
       dimension itext(10)
       common /p3d/ xpr,ypr,caxt,saxt,cayt,sayt
       real xfmin,xfmax,xfbas
       equivalence (itext(1),ctext)   
       dgrad = atan(1.0)/45.
       sq2 = 1.0/sqrt(2.)
c
c      Fermi Dirac distribution
       xfbas = 1.5*0.999
       xfref = 2.0
       xthr = 0.975
       beta = 2.0
c
       xfmin = 1.5
c
       wef = 1.0/(exp(beta*(xfmin-xfbas*xfref)/xfbas)+1.0)
       wthres = 0.5+xthr*(wef-0.5)
       write(6,*) "wthres", wthres
c
cc
       open(8,file="fdplot.ps")
       call Hplots(1,0,8,0)
       call typset(0.0,0.0)
       call zpick(3,0,is)
c
c   plot misfit distribution
       call typstr(4.8,7.3,0.33,"Weight",90.,6)
       call typstr(12.0,3.0,0.33,"Misfit",0.,6)
       call pen(1,0)
       call zpick(3,0,is)
       xo = 6.
       yo = 4.
       psc = 2.5
       csc = 10.0
c
       call plot(xo,yo,3)
       call plot(xo+5.0*psc,yo,2)
       call plot(xo,yo,3)
       call plot(xo,yo+csc,2)
c?       do j=1,501,50
c?         xj = float(j-1)*0.01
c?         call plot(xo+xj*psc,yo,3)
c?         call plot(xo+xj*psc,yo+0.3,2)
c?         call typnum(xo+xj*psc-0.1,yo-0.5,0.3,xj,0.0,2)
c?       enddo
       do c=0.,1.0,0.25
         call plot(xo,yo+c*csc,3)
         call plot(xo+0.3,yo+c*csc,2)
         call typnum(xo-1.0,yo+c*csc-0.1,0.3,c,0.0,2)
       enddo 
       kp = 3        
       do j=1,501
         xmis = float(j-1)*0.01
         we = 1.0/(exp(beta*(xmis-xfbas*xfref)/xfbas)+1.0)
         call plot(xo+xmis*psc,yo+we*csc,kp)
         kp = 2
       enddo
c
       xmas = 1.5
       wf = 1.0/(exp(beta*(xmas-xfbas*xfref)/xfbas)+1.0)
       xmbs = 3.0
       wg = 1.0/(exp(beta*(xmbs-xfbas*xfref)/xfbas)+1.0)
       xmis = 1.56558
       we = 1.0/(exp(beta*(xmis-xfbas*xfref)/xfbas)+1.0)
       write(6,*) we
c
       call pen(15,0)
       call plot(xo,yo+we*csc,3)
       call plot(xo+xmis*psc,yo+we*csc,2)
       call plot(xo+xmas*psc,yo+wf*csc,2)
       call plot(xo,yo+wf*csc,2)
       write(8,*) "cf"                     
       call pen(3,0)
       call plot(xo+xmis*psc,yo,3) 
       call plot(xo+xmis*psc,yo+we*csc,2) 
       call plot(xo,yo+we*csc,2)
       call typstr(xo+xmis*psc+0.1,yo+0.2,0.25,"thres",0.,5)       
       call pen(2,0)
       call plot(xo+xmas*psc,yo,3) 
       call plot(xo+xmas*psc,yo+wf*csc,2) 
       call plot(xo,yo+wf*csc,2) 
       call plot(xo+xmbs*psc,yo,3) 
       call plot(xo+xmbs*psc,yo+wg*csc,2) 
       call plot(xo,yo+wg*csc,2)
       call typstr(xo+xmas*psc-0.35,yo-0.4,0.25,"min",0.,3)
       call typstr(xo+xmbs*psc-0.35,yo-0.4,0.25,"ref",0.,3)             
c       
       call Hplots(0,0,8,1)
       stop
       end 
