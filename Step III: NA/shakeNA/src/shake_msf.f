c*************************************************************************
c
c      program shake_msf
c
c      Seismology Group, Research School of Earth Sciences
c      The Australian National University
c
c     Calculation of misfit or an externally provided hyypocentre 
c     for comparison with results of non-linear inversion 
c
c     Builds on 
c       Shake       - using iasp software for traveltime calculation
c
c

c----------------------------------------------------------------------*------*
c  SHAKE component:
c
c     Uses Buland tau-p based travel time procedures
c     to improve computational efficiency
c
c
c     AUTHORS:  Brian L.N. Kennett (RSES,ANU)
c               March 1990, January 1999 for shakeNA
c
c     The routines build on the work of Genet M. Edmonson in January 1987     
c     who modified the program JANET by M.S. Sambridge designed for
c     locating local earthquakes for use in a spherical geometry with 
c     teleseismic arrivals.
c
c
c
c*************************************************************************
c
	Program shake_msf
c----------------------------------------------------------------------*------*
c
c     CONSTANTS
c
      integer EVT, STAT, RUN, STOUT, RES, TABL
c
      parameter (RUN=15, STOUT=6, EVT=17, STAT=18, RES=19, TABL=11)
c
c     STOUT   unit number of shake output file
c     EVT     unit number of file containing data about quake
c     STAT    unit number of file containing information about stations
c     RUN     unit number of file containing run information
c     RES     unit number of file for auxiliary results
c
*----------------------------------------------------------------------*------*
c
c     COMMON BLOCK DECLARATIONS
      include 'shakcom.inc'
c
c     cjgl              used for choice of residual statistic
c                       (GAussian,JEffreys,L1)
c
*----------------------------------------------------------------------*------*
c
c     VARIABLES
c
      character*20 modnam
      character*40 cfile1
      character*2  cjgl
c
c
c     cfile             used for file names
c     modnam            name for travel time tables 
c                       (modnam.hed,modnam.tbl)  
c
*----------------------------------------------------------------------*------*
c
c
        character*8 evida
        character*4 wmo1
        integer wyr,wday,whr,wmin,ncev,charcn8
c
        modnam = "ak135"
        call assign(20,2,'ttim1.lis',9)
        call tabin(TABL,modnam)
        phlst(1) = 'basic'
        phlst(2) = 'PP   '
        phlst(3) = '   '
        cjgl = 'L1'
        if(cjgl.eq.'JE') then
          jgl = 2
        elseif(cjgl.eq.'L1') then
          jgl = 3
        else
          jgl = 1
        endif 
        write(6,*) 'Residual Statistic:',cjgl
c
        open(UNIT=RUN,file="epshyp.in")
        open(UNIT=RES,FILE="res.out")
c							
        open(2,file="iwref.lis")
        do m=1,157      
         read(2,1010) evida,wlat,wlng,wdep,
     &                      wyr,wmo1,wday,whr,wmin,wsec 
         write(7,*) evida
         ncev = charcn8(evida)
  
         cfile1 = "Events-PS/"//evida(1:ncev)//".evd"
         open(UNIT=EVT, FILE=cfile1)
         rewind(EVT)
         cfile1 = "Events-PS/"//evida(1:ncev)//".sta"
         open(UNIT=STAT, FILE=cfile1)
         rewind(STAT)
c      
c     OBTAIN QUAKE, STATION & RUN DATA
c
        call input(STOUT, EVT, STAT, RUN, RES)
        prnt(1) = .false.
        prnt(2) = .false.
        prnt(3) = .false.
        call brnset(2,phlst,prnt)
c
c
c       Open direct access tables for ellipticity corrections
        open(21,file='elcordir.tbl',access='direct',
     &     form='formatted',recl=80) 
c
*----------------------------------------------------------------------*------*
c
c     FORMAT STATEMENTS
c
1000  format (1X, A)
1003  format (1X, A)
1004  format ()
1005  format (1X, T18, A, T28, A)
1006  format (1X, T3, A, T13, F8.3, T23, F8.3)
c
*----------------------------------------------------------------------*------*
c
c
        open(1,file="resid.msf") 
        x = wlat
        y = wlng
        z = wdep
        t = 3600.0*float(whr-shr)+60.0*float(wmin-smin)+wsec
        smo = wmo1
c
        call displmsf(1,jgl, x, y, z, t) 
        close(EVT)
        close(STAT)
        rewind(RES)
        rewind(RUN)
      enddo
c
1010  format(a,"| ",3f10.3,"   ",i4," ",a," ",i2," ",i2,":",i2,":",f6.3)
	stop
999   write(stout,1003) 'ERROR on input (shake_msf)'       
	end

*--------------------------------------------------------------------- + ---- +
c DISPLAY ROUTINES
*--------------------------------------------------------------------- + ---- +
c                                                                     displn
      subroutine displmsf (stout, jgl, lt, lg, dp, sc)
c
c     AUTHOR:  Genet M. Edmondson  RSES, ANU
c     DATE:    January 1987 
c     MODIFIED: B.L.N. Kennett March 1990  
c               to allow for changes to residual calculation
c               B.L.N. Kennett May 1991
c               to include slowness, azimuth data
c               B.L.N. Kennett Jan 1999 
c               shakeNA display
c     PURPOSE:
c     * displays results from shake program
c
c     Subroutines and functions required:
c       cstaq - returns residual statistic for given location
c       (ydist, trtm, obst, month - used by cstat)
c
*----------------------------------------------------------------------*------*
c
c     PARAMETERS
c
      integer stout, res
c
      real    lt, lg, dp, sc
c
c     stout   unit number of standard output file
c     res     unit number of file to output results to
c     jgl     index for Residual Statistic
c     LA, LN, D, S
c             bounds on hypocentre location
c     lt, lg, dp, sc
c             hypocentre location
c
*----------------------------------------------------------------------*------*
c
c     COMMON BLOCK DECLARATIONS
      include 'shakcom.inc'
c
*----------------------------------------------------------------------*------*
c
c     VARIABLES
c
      integer  i, jgl
      integer  tmp, tmpmin, tmphrs
      integer  carmin, carhr
      integer  nc, charcn
c
      real     calct, delta, bazim
      real     tmpsec
      real     carsec
      real     resqkp, ressrc
      character*2 cjgl
c
c     i       loop counter
c     jgl     index for residual statistic
c     tmp etc temporary variables
c     carsec, carmin, carhr
c             calculated arrival times at stations
c     calct   theoretical travel time from source to receiver
c     delta   angular distance from source to receiver
c     trash   dummy variable
c     resqkp  residual statistic for results from shake program
c     ressrc  residual statistic for approx source
c
*----------------------------------------------------------------------*------*
c
c     FUNCTIONS
c
      real    cstaq
c
c     cstaq - returns residual statistic
c
*----------------------------------------------------------------------*------*
c
c
c                                    Find residual statistics
      resqkp = cstaq(stout, jgl, sc, lt, lg, dp)
      ressrc = cstaq(stout, jgl, ssec, slat, slng, sdep)
      write(7,*) resqkp,ressrc
c     
      If(jgl.eq.3) cjgl = "LI"      
      If(jgl.eq.2) cjgl = "JE"      
      If(jgl.eq.1) cjgl = "GA"      
c
      write(stout,1000) 'SOLUTION',cjgl
      write(stout,1001) source
      write(stout,1002) 'DATE:'
      write(stout,1003) 'Year', syr, syr
      write(stout,1004) 'Month', smo, smo
      write(stout,1005) 'Day', sday, sday
      write(stout,1005) 'Hour', shr, shr
      write(stout,1005) 'Minute', smin, smin
      write(stout,1006) 'Second', sc, ssec
      write(stout,1002) 'LOCATION:'
      write(stout,1008) 'Latitude', lt, slat
      write(stout,1008) 'Longitude', lg, slng
      write(stout,1008) 'Depth', dp, sdep
      write(stout,1009) 'RESIDUAL:', resqkp, ressrc
   
c
c     CALCULATE & ECHO ARRIVAL TIMES AT STATIONS FOR ORIGIN
c     AT CALCULATED SOLUTION
c
      rslat = (90.-lt)*0.017453292
      call ellref(rslat)
      write(stout,1012) 'Arrival Times for estimated hypocentre'
      write(stout,1015)
c                                 all main phases 
      call depset(dp,usrc)
      numtaz = numtim+numazi
      do 10, i = 1, numevt
      nc = charcn(evphcd(i))
c
c       CALCULATE TRAVEL TIME FROM SOLUTION TO STATION
c
        call ydist(lt, lg, stlat(evind(i)), stlong(evind(i)),
     ^             delta, cazim, bazim)
        call trtm(delta, MAXTT, n, tt,dtdd,dtdh,dddp,phcd)
        cazim = cazim+180.
        if(cazim.gt.360.) cazim = cazim-360.
        do 23 m=1,n
 23      continue
        if(nc.eq.1.and.evphcd(i)(1:1).eq.'P') then
          calct = tt(1)
          calcp = dtdd(1)
c          write(res,*) i,' fast P'
        else
          do 21 l=1,n
            if(evphcd(i)(1:nc).eq.phcd(l)(1:nc)) then
              calct = tt(l)
              calcp = dtdd(l)
c              write(res,*) i,nc,'  ',evphcd(i)
              go to 22
            endif
 21       continue
        endif
c
c       CALCULATE ACTUAL ARRIVAL TIME AT STATION
c       NB: this only takes into account crossing day, month and
c           year boundaries
c
c       tmpsec = ssec + calct
 22     continue
        edist = delta*0.017453292
        bazr = bazim*0.017453292
        call ellcor(edist, bazr, z, evphcd(i)(1:nc), etcor)
        calct = calct+etcor
        tmpsec = sc + calct
        tmpmin = tmpsec / 60
        carsec = tmpsec - (tmpmin * 60)
        tmp = smin + tmpmin
        tmphrs = tmp / 60
        carmin = tmp - (tmphrs * 60)
        carhr = shr + tmphrs
        if (i.le.numtim) then
          residf = obst(stout,sc,i)-calct
          write(stout,1016) evstcd(i),evphcd(i),evtmhr(i),evtmmi(i),
     &          evtmsc(i),residf,carhr,carmin,carsec,delta,bazim,calcp
        elseif (i.gt.numtim .and. i.le.numtaz) then
          resaz =  evtazi(i)-cazim         
          write(stout,1017) evstcd(i),evphcd(i),
     &         evtazi(i),resaz,cazim,delta,bazim
        elseif (i.gt.numtaz) then
          resslo =  evtslo(i)-calcp
          write(stout,1017) evstcd(i),evphcd(i),
     &         evtslo(i),resslo,calcp,delta,bazim
        endif
c
10      continue
c
*----------------------------------------------------------------------*------*
c
c     FORMAT STATEMENTS
c
1000  format (1X, T15, 'HYPOCENTRE ', A, T40, A)
1001  format (1X, T15, 'IWREF', T40, A)
1002  format (1X, T3, A, A)
1003  format (1X, T5, A, T15, I4, T40, I4)
1004  format (1X, T5, A, T15, A, T40, A)
1005  format (1X, T5, A, T15, I2, T40, I2)
1006  format (1X, T5, A, T15, F8.3, T40, F8.3)
1007  format (1X, T5, A, T15, F8.3, T25, F8.3, T40, F8.3, T50, F8.3)
1008  format (1X, T5, A, T15, F8.3, T40, F8.3)
1009  format (1X, T3, A, T15, F12.3, T40, F12.3)
1010  format (1X, T9, 'CODE', T14, 'YEAR', T19, 'MONTH',
     & T25, 'DAY', T29, 'HOUR', T34, 'MIN', T38, 'SEC')
1011  format (1X, T3, I3, T9, A, T14, I4, T19, A, T25, I2, T29, I2, T34,
     & I2, T38, F5.2)
1012  format (/1X, A)
1015  format (1X, T2, 'CODE', T15, 'HOUR', T20, 'MIN',
     & T25, 'SEC', T33, 'RESID', T43, 'chr' , T48, 'cmn',
     & T53, 'csec', T60, 'delta', T70, 'baz', T77,'  slow')
1016  format (1X, T2, A, T7, A, T15, I2, T20, I2, T25,
     & F5.2, T33,  F5.2, T43, I2, T48, I2, T53, F5.2, 
     & T60, F6.2, T70, F5.1, T77, F8.4)
1017  format (1X, T2, A, T7, A, T20, F10.3, T32, F6.2, T48, F10.3, 
     & T60, F6.2, T70, F5.1)
c
*----------------------------------------------------------------------*------*
c
      return
      end
C
c                                                                      charcn
      integer function charcn8 (c)
c
c     AUTHOR:  Brian L.N. Kennett RSES ANU
c     DATE:    December 1985
c     PURPOSE:
c     * This routine returns the number of non-blank
c       characters in a string  ( up to first blank)
c
*----------------------------------------------------------------------*------*
c
c     PARAMETERS
c
      character*8 c
c
c     c is the character string to be counted
c
*----------------------------------------------------------------------*------*
c
c     VARIABLES
c
c     integer charcn8
c
*----------------------------------------------------------------------*------*
c
c
      do 10 k=1,8
        if(c(k:k).eq.' ') then
            charcn8 = k-1
            return
        endif
 10   continue
      charcn8 = 8 
c
*----------------------------------------------------------------------*------*
      return
      end
c
