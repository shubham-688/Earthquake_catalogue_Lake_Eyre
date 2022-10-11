c
c            Program NALIN : Used to plot curves of data from
c                            Neighbourhood Algorithm Location procedure
c                              data points read in from lu 9
c
c

      Program NAlin
      dimension it(60),x1(60),x2(60),x3(60)
      character c1*13,c2*17,c3*9,c4*9
      dimension xar(250),yar(250)     
      character*20 dfile
      DIMENSION ITX(18),ITY(18),ITEXT(18)
      write(6,*)" Line plotting program "
      write(6,*)" "
*     write(6,*)" How many data points are there "
      read(5,*)igen
      write(6,*) "no of generations: ", igen
      call hplots(1,0,8,1)
      call pen(1,0)
*     write(6,*)" Start number for curves"
      read(5,*) ig1
*     write(6,*)" End number for curves"
      read(5,*) ig2
*     write(6,*) "enter data display mode:"
*     write(6,*) "linear: 1, semi-log: 2, log-log :3 "
      read(5,*)istyle
*     write(6,*) "name of data file"
      read(5,*) dfile
      open(9,file=dfile)
c
c
*       Plotting of annotated frames
*
*
*
*                               Input of parameters
 1    continue
*     WRITE(6,217)
*217  FORMAT(//" Enter X and Y origin (cm)    : ")
      READ(5,*) XOR,YOR
      Xdist = Xor
      Ydist = Yor
*     WRITE(6,218)
*218  FORMAT(/," Enter X-axis and Y-axis length (cm) : ")
      READ(5,*) XASL,YASL
*     WRITE(6,703)
*703  FORMAT(/," Enter Xmin and Xmax as function values:")
      READ(5,*) FX0,FXM
*     WRITE(6,704)
*704  FORMAT(/," Enter Ymin and Ymax as function values:")
      READ(5,*) FY0,FYM
*     WRITE(6,705)
*705  FORMAT(
*    #/," Enter step in fnct. values between 2 LARGE tics for X and Y:")
      READ(5,*) DFX,DFY
*     WRITE(6,712)
*712  FORMAT(
*    #/," Enter step in fnct. values between 2 SMALL tics for X and Y:")
      READ(5,*) DDFX,DDFY
*     WRITE(6,706)
*706  FORMAT(/," Enter number of decimals in label for X and Y:")
      READ(5,*) NXDEC,NYDEC
*     write(6,720)
*720  format(/," Enter character size for text, title (cm): ")
      read(5,*) sizl,sizt
*     WRITE(6,707)
*707  FORMAT(/," Enter text X-axis (max. 72 char.):")
      read(5,903)ITX
*     write(6,708)
*708  format(/," Enter text Y-axis (max. 72 char.):")
      read(5,903) ITY
*     write(6,709)
*709  format(/," Enter plot  TITLE (max. 72 char.)")
      read(5,903) ITEXT
 903  format(18A4)
*     write(6,*) " Enter number of font to be used - default 3 "
      read(5,*) ifont
      if(ifont.eq.0) ifont=3
      write(6,*) "Enter pen numbers for frame box, ticks, text"
      read(5,*) ipen1,ipen2,ipen3
C
      XTICL=abs(XASL/(FXM-FX0)*abs(DFX))
      YTICL=abs(YASL/(FYM-FY0)*abs(DFY))
      NXSUB=abs(DFX/DDFX)
      NYSUB=abs(DFY/DDFY)
C
C                                       Check
      WRITE(6,710)
     %  XOR,YOR,XASL,YASL,FX0,FXM,FY0,FYM,DFX,DFY,DDFX,DDFY,
     %  NXDEC,NYDEC,SIZL,SIZT,ITX,ITY,ITEXT,ifont,ipen1,ipen2,ipen3
710   FORMAT(//
     #" Frame parameters:",//
     #" Xor,Yor         :"f8.3,x,f8.3,"  cm",/,
     #" Xasl,Yasl       :"f8.3,x,f8.3,"  cm",/,
     #" Xmin,xmax       :"f8.3,x,f8.3,"  in function values",/,
     #" Ymin,Ymax       :"f8.3,x,f8.3,"  in function values",/,
     #" Large tic space :"f8.3,x,f8.3,"  in function values",/,
     #" Small tic space :"f8.3,x,f8.3,"  in function values",/,
     #" # dec. in label :"3X,I5,3X,I5,/,
     #" char size text  :",f8.3,"  cm",/,
     #" char size title :",f8.3,"  cm",/,
     #" X-txt:"18a4,/,
     #" Y-txt:"18a4,/,
     #" Title:"18a4,/,
     #" Font     :",i5,/,
     #" Frame pen:",i5,/,
     #" Ticks pen:",i5,/,
     #" Text  pen:",i5)
c
C                                       initialize plotting
      call typset(0.0,0.0)
      call zpick(ifont,0,is)
c
C
      hscale = xasl/(fxm-fx0)
      vscale = yasl/(fym-fy0)
 5    continue
      ldash = 12
c
c                      Plot curves
c
      read(9,*)
      read(9,*)
      call pen(ipen,0)
      do 450 i=1,igen
        read(9,91) c1,it(i),c2,x1(i),c3,x3(i),c4,x2(i)
91      format(a,i2,a,f8.5,a,f8.5,a,f8.5)
        write(6,*) c1,it(i),x1(i),x3(i),x2(i)
        do j=1,12
          read(9,*)
        enddo
450   continue
20      call dashln(ldash,0)
      do 520 jq=ig1,ig2
        if(jq.eq.1) then
         do 451 i=1,igen
           xar(i) = it(i)
           yar(i) = x1(i)
           ipn = 3
 451     continue
        elseif(jq.eq.2) then
         do 452 i=1,igen
           xar(i) = it(i)
           yar(i) = x2(i)
           ipn = 4
 452     continue
        elseif(jq.eq.3) then
         do 453 i=1,igen
           xar(i) = it(i)
           yar(i) = x3(i)
           ipn = 5
 453     continue
        endif
        call pen(ipn,0)
          var = fym
          if(yar(1).lt.var) var = yar(1)
          if(istyle.eq.1) then
            xl = xar(1)
            yl = var
          else if(istyle.eq.2) then
            xl = alog10(xar(1))
            yl = var
          else if (istyle.eq.3) then
           xl = alog10(xar(1))
           yl = alog10(var)
          end if
          xu = xdist
          xv = xdist + hscale*(xl+0.5 -fx0)
          yv = ydist + vscale*(yl -fy0)
          call plot(xu,yv,3)
          if(ldash.eq.0)then
           call plot(xv,yv,3)
           call plot(xv,yv,2)
          else
           call plot(xv,yv,2)
         end if
         xu = xv
        do 501 i = 2,igen
          call plot(xu,yv,3)
          var = fym
          if(yar(i).lt.var) var = yar(i)
          if(istyle.eq.1) then
            xl = xar(i)
            yl = var
          else if(istyle.eq.2) then
            xl = alog10(xar(i))
            yl = var
          else if (istyle.eq.3) then
           xl = alog10(xar(i))
           yl = alog10(var)
          end if
          xv = xdist + hscale*(xl+0.5 -fx0)
          yv = ydist + vscale*(yl -fy0)
          if(ldash.eq.0)then
           call plot(xu,yv,3)
           call plot(xv,yv,2)
         else
          call plot(xu,yv,2)
          call plot(xv,yv,2)
         end if
         xu = xv
501    continue
520   continue
c
c     draw frame with labels and text
c
      write(6,*)" Calling frame"
      CALL aframe
     $(XOR,YOR,XASL,YASL,FX0,XTICL,DFX,NXSUB,NXDEC,FY0,YTICL,
     $DFY,NYSUB,NYDEC,ITX,ITY,ITEXT,SIZL,SIZT,ipen1,ipen2,ipen3)
      write(6,*)" Frame done"
10    continue
      call hplots(0,0,8,1)
      stop
      end
