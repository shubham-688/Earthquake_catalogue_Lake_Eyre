*-----        Location Routines      -------------------------------- + ---- +
c                                                                     ydist
      subroutine ydist
     ^        (dlats, dlons, dlatr, dlonr, delta, cazim, bazim)
c
c     AUTHOR:  Brian L.N. Kennett  RSES, ANU
c     DATE:    January 1985
c     PURPOSE:
c             YDIST        Calculates distance and azimuth
c                          for spheroidal earth between
c                          specified geographic source and
c                          receiver station coordinates
c
*----------------------------------------------------------------------*------*
c     PARAMETERS
c
      real    dlats, dlons, dlatr, dlonr, delta, cazim, bazim
c
c     dlats  latitude of source
c     dlons  longitude of source
c     dlatr  latitude of receiver
c     dlonr  longitude of receiver
c     delta  angular distance
c     cazim  apparent azimuth at an array 
c     bazim   azimuth from epicentre to receiver
c
*----------------------------------------------------------------------*------*
c
*      implicit real*8 (a-h,o-z)
      real   azima
      real*8 gra,ecc,re,ec1,pi,pib2,degr
      real*8 rlats,rlons,rlatr,rlonr
      real*8 glats,glatr,sps,cps,spr,cpr,rs,rr
      real*8 trs,prs,trr,prr
      real*8 AS,BS,CS,DS,ES,GS,HS,KS
      real*8 AR,BR,CR,DR,ER,GR,HR,KR
      real*8 cosdr,deltar,sindr,deltak
      real*8 szs,czs,szr,czr
      real*8 x,y,e
c
c                          radius on spheroid
      gra(x,y,e) = dsqrt( (1.0d0-e)**2 /
     &                   ((1.0d0-e*y)**2 + e*e*x*y ) )
      ecc = 0.003367
      re = 6378.388
      ec1 = (1.0d0-ecc)**2
      pi = 3.141592653589793
      pib2 = pi/2.0
      degr = pi/180.0
      rlats = dlats*degr
      rlons = dlons*degr
      rlatr = dlatr*degr
      rlonr = dlonr*degr
c                          geocentric coordinates
      glats = datan2 ( ec1*dsin(rlats) ,dcos(rlats) )
      glatr = datan2 ( ec1*dsin(rlatr) ,dcos(rlatr) )
      sps = dsin(glats)**2
      cps = dcos(glats)**2
      spr = dsin(glatr)**2
      cpr = dcos(glatr)**2
c                          radii at source,receiver
      rs = re*gra(sps,cps,ecc)
      rr = re*gra(spr,cpr,ecc)
c
      trs = pib2 - glats
      prs = dlons*degr
      trr = pib2 - glatr
      prr = dlonr*degr
c                          direction cosines for source
      AS = dsin(trs)*dcos(prs)
      BS = dsin(trs)*dsin(prs)
      CS = dcos(trs)
      DS = dsin(prs)
      ES = -dcos(prs)
      GS = dcos(trs)*dcos(prs)
      HS = dcos(trs)*dsin(prs)
      KS = -dsin(trs)
c                          direction cosines for receiver
      AR = dsin(trr)*dcos(prr)
      BR = dsin(trr)*dsin(prr)
      CR = dcos(trr)
      DR = dsin(prr)
      ER = -dcos(prr)
      GR = dcos(trr)*dcos(prr)
      HR = dcos(trr)*dsin(prr)
      KR = -dsin(trr)
c                          distance
      cosdr = AS*AR + BS*BR + CS*CR
      deltar = dacos(cosdr)
      sindr = dsin(deltar)
c
      deltak = deltar*0.5d0*(rr+rs)
      delta = deltar/degr
c                          azimuth
      szs = DS*AR + ES*BR
      czs = GS*AR + HS*BR + KS*CR
      szr = DR*AS + ER*BS
      czr = GR*AS + HR*BS + KR*CS
c                          azima - azimuth to source
c                          bazim - backazimuth from source
c                          cazim - apparent azimuth at an array
      if (szr.eq.0.0) then
        bazim = 0.0
        azima = 180.0
      else
        bazim = datan2(-szs ,-czs ) /degr
        azima = datan2(-szr ,-czr ) /degr
      end if
      if( bazim .lt. 0.0) bazim = bazim + 360.0
      cazim = azima + 180.0
      if( azima .lt. 0.0) azima = azima + 360.0
c
      if( cazim.lt. 0.0) cazim = cazim + 360.0
c
      return
      end
*--------------------------------------------------------------------- + ---- +
c                                                                      cstaq
      real   function cstaq(stout, jgl, t, x, y, z)
c
c     AUTHOR:  Genet M. Edmondson  Brian L.N. Kennett RSES, ANU
c     DATE:    January 1987
c     MODIFIED: March 1990 to include Buland technique for travel time
c               calculation
c               May 1991 to include azimuth and slowness information
c
c     PURPOSE:
c     * returns squared residual statistic C for origin time t and source
c       at x, y, z
c
c     Subroutines and functions required:
c       ydist - calculates angular distance from source to receiver
c       depset - sets up trtm calculation for given depth
c       trtm - calculates theoretical traveltime from source to receiver
c       obst - returns observed arrival time in seconds
c       (month - used by obst)
c     Original form of ellipticity corrections
*----------------------------------------------------------------------*------*
c
c     PARAMETERS
c
      integer   stout, jgl
c
      real    t, x, y, z
c
c     stout     unit  number of (standard) output file
c     jgl       index for Residual Statistic
c     t         dummy variable
c     x, y, z   lat, long, depth of source
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
      integer i
      integer charcn, nc
c
      real    delta, bazim, calct
c
c     i       loop counter over phase events
c     delta   angular distance from source to receiver
c     bazim   back azimuth from source
c     cazim   apparent azimuth at array
c     calct   theoretical travel time from source to receiver
c     calcp   theoretical slowness
c
*----------------------------------------------------------------------*------*
c
c     FUNCTIONS
c
      real    obst
c
c     obst - returns observed travel time from source to receiver
c
*----------------------------------------------------------------------*------*
       data r2pi,f1,f2 /0.3989428,0.05,0.95/
c
c                                   stage 1 - ellipicity correction for
c                                   centre of grid
        rslat = (90.-x)*0.017453292
        call ellref(rslat)
c                                    initialise travel time calculation
c                                    (all main phases 1,14: P 1,2)
      call depset(z,usrc)
      edepth = z
c
      cstatm = 0.0
      cstata = 0.0
      cstats = 0.0
      cstaq = 0.0
      n1 = numtim+1
      n2 = numtim+numazi
      n3 = n2+1
      edepth = z
c
      do 10, i = 1, numevt
        nc = charcn(evphcd(i))
        call ydist(x, y, stlat(evind(i)), stlong(evind(i)),
     ^             delta, cazim, bazim)
        cazim = cazim+180.
        if(cazim.gt.360.) cazim = cazim-360.
        call trtm(delta, MAXTT, n, tt,dtdd,dtdh,dddp,phcd)
        if(nc.eq.1.and.evphcd(i)(1:1).eq.'P') then
          calct = tt(1)
          calcp = dtdd(1)
        else
          do 21 l=1,n
            if(evphcd(i)(1:nc).eq.phcd(l)(1:nc)) then
              calct = tt(l)
              calcp = dtdd(l)
              go to 22
            endif
 21       continue
        endif
c
 22     edist = delta*0.017453292
        bazr = bazim*0.017453292
        call ellcor(edist, bazr, z, evphcd(i)(1:nc), etcor)
        calct = calct+etcor
c
       if(i.le.numtim) then
         re = (obst(stout,t,i) - calct)
         si = sigma(i)
         If(jgl.eq.1) Then
           cstatm = cstatm + (re/si)**2
         Elseif(jgl.eq.2) Then
           cstatm = cstatm + (abs(re/si))**1.25
         Elseif(jgl.eq.3) Then
           cstatm = cstatm + abs(re/si)
         Else
         Endif
c
       elseif(i.ge.n1 .and. i.lt.n2) then
         re = (evtazi(i) - cazim)
         si = sigma(i)
         If(jgl.eq.1) Then
           cstata = cstata + (re/si)**2
         Elseif(jgl.eq.2) Then
           cstata = cstata + (abs(re/si))**1.25
         Elseif(jgl.eq.3) Then
           cstata = cstata + abs(re/si)
         Else
         Endif
c
       elseif(i.ge.n3) then
         re = (evtslo(i) - calcp)
         si = sigma(i)
         If(jgl.eq.1) Then
           cstats = cstats + (re/si)**2
         Elseif(jgl.eq.2) Then
           cstats = cstats + (abs(re/si))**1.25
         Elseif(jgl.eq.3) Then
           cstats = cstats + abs(re/si)
         Else
         Endif
       endif
c 
 10   continue
      cstaq = (cstatm + cstata + cstats)/numevt
c%     write(stout,*) 'cstaq = ', cstaq
c
*----------------------------------------------------------------------*------*
c
      return
      end
*--------------------------------------------------------------------- + ---- +
c                                                                      cstaq
      real   function cstaqke(stout, jgl, t, x, y, z)
c
c     AUTHOR:  Genet M. Edmondson  Brian L.N. Kennett RSES, ANU
c     DATE:    January 1987
c     MODIFIED: March 1990 to include Buland technique for travel time
c               calculation
c               May 1991 to include azimuth and slowness information
c
c     PURPOSE:
c     * returns squared residual statistic C for origin time t and source
c       at x, y, z
c
c     Subroutines and functions required:
c       ydist - calculates angular distance from source to receiver
c       depset - sets up trtm calculation for given depth
c       trtm - calculates theoretical traveltime from source to receiver
c       obst - returns observed arrival time in seconds
c       (month - used by obst)
c     Full set of ellipticity corrections
*----------------------------------------------------------------------*------*
c
c     PARAMETERS
c
      integer   stout, jgl
c
      real    t, x, y, z
c
c     stout     unit  number of (standard) output file
c     jgl       index for Residual Statistic
c     t         dummy variable
c     x, y, z   lat, long, depth of source
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
      integer i
      integer charcn, nc
c
      real    delta, bazim, calct
c
c     i       loop counter over phase events
c     delta   angular distance from source to receiver
c     bazim   back azimuth from source
c     cazim   apparent azimuth at array
c     calct   theoretical travel time from source to receiver
c     calcp   theoretical slowness
c
*----------------------------------------------------------------------*------*
c
c     FUNCTIONS
c
      real    obst
c
c     obst - returns observed travel time from source to receiver
c
*----------------------------------------------------------------------*------*
       data r2pi,f1,f2 /0.3989428,0.05,0.95/
       logical lfa
       character*8 phase
       real edepth
c
c                                   stage 1 - ellipicity correction for
c                                   centre of grid
        rslat = (90.-x)*0.017453292
        call kellref(rslat)
c                                    initialise travel time calculation
c                                    (all main phases 1,14: P 1,2)
      call depset(z,usrc)
      edepth = z
c
      cstatm = 0.0
      cstata = 0.0
      cstats = 0.0
      cstaq = 0.0
      n1 = numtim+1
      n2 = numtim+numazi
      n3 = n2+1
c
      do 10, i = 1, numevt
        nc = charcn(evphcd(i))
        call ydist(x, y, stlat(evind(i)), stlong(evind(i)),
     ^             delta, cazim, bazim)
        cazim = cazim+180.
        if(cazim.gt.360.) cazim = cazim-360.
        call trtm(delta, MAXTT, n, tt,dtdd,dtdh,dddp,phcd)
        if(nc.eq.1.and.evphcd(i)(1:1).eq.'P') then
          calct = tt(1)
          calcp = dtdd(1)
        else
          do 21 l=1,n
            if(evphcd(i)(1:nc).eq.phcd(l)(1:nc)) then
              calct = tt(l)
              calcp = dtdd(l)
              go to 22
            endif
 21       continue
        endif
c
 22     edist = delta*0.017453292
        bazr = bazim*0.017453292
        phase = evphcd(i)(1:nc)
        call kellcor(phase,delta,edepth,rslat,bazim,etcor,lfa)
        calct = calct+etcor
c
       if(i.le.numtim) then
         re = (obst(stout,t,i) - calct)
         si = sigma(i)
         If(jgl.eq.1) Then
           cstatm = cstatm + (re/si)**2
         Elseif(jgl.eq.2) Then
           cstatm = cstatm + (abs(re/si))**1.25
         Elseif(jgl.eq.3) Then
           cstatm = cstatm + abs(re/si)
         Else
         Endif
c
       elseif(i.ge.n1 .and. i.lt.n2) then
         re = (evtazi(i) - cazim)
         si = sigma(i)
         If(jgl.eq.1) Then
           cstata = cstata + (re/si)**2
         Elseif(jgl.eq.2) Then
           cstata = cstata + (abs(re/si))**1.25
         Elseif(jgl.eq.3) Then
           cstata = cstata + abs(re/si)
         Else
         Endif
c
       elseif(i.ge.n3) then
         re = (evtslo(i) - calcp)
         si = sigma(i)
         If(jgl.eq.1) Then
           cstats = cstats + (re/si)**2
         Elseif(jgl.eq.2) Then
           cstats = cstats + (abs(re/si))**1.25
         Elseif(jgl.eq.3) Then
           cstats = cstats + abs(re/si)
         Else
         Endif
       endif
c 
 10   continue
      cstaq = (cstatm + cstata + cstats)/numevt
c%     write(stout,*) 'cstaq = ', cstaq
c
*----------------------------------------------------------------------*------*
c
      return
      end
*--------------------------------------------------------------------- + ---- +
c                                                                      obst
      real   function obst (stout, t, i)
c
c     AUTHOR:  Genet M. Edmondson  RSES, ANU
c     DATE:    January 1987
c     PURPOSE:
c     * returns observed travel time (in seconds)
c       with origin time t (in seconds relative to smin)
c       and observed arrrival time given by evtm values
c     * takes into account the possibility that origin and arrival times
c       may cross minute, hour, ...year boundaries
c     * treats all years divisible by 4 as leap years
c
c     Subroutines and functions required:
c       month - converts a character string of roman numerals to an integer
c
*----------------------------------------------------------------------*------*
c
c     PARAMETERS
c
      integer   stout, i
c
      real    t
c
c     stout   unit number of (standard) output file
c     i       number of reading for arrival times
c     t       origin time (seconds relative to smin)
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
      integer  j
      integer  moevt, mosrc
      integer  MONLEN(12)
c
c     j           loop counter
c     moevt       month of arrival time as an integer
c     mosrc       month of origin time as an integer
c     MONLEN(i)   gives the number of days in month i
c
*----------------------------------------------------------------------*------*
c
c     FUNCTIONS
c
      integer   month
c
c     month - converts a character string of roman numerals to an integer
c
*----------------------------------------------------------------------*------*
c
c
c     INITIALISE MONLEN ARRAY
c
      DATA  (MONLEN(j), j = 1, 12) / 31, 28, 31, 30, 31, 30, 31, 31,
     &  30, 31, 30, 31 /
c
      moevt = month(evtmmo(i))
      mosrc = month(smo)
c
      if (evtmyr(i) .eq. syr) then
        if (moevt .eq. mosrc) then
          if (evtmdy(i) .eq. sday) then
            if (evtmhr(i) .eq. shr) then
              if (evtmmi(i) .eq. smin) then
                if (evtmsc(i) .lt. t) then
c
c                 ERROR: ARRIVAL TIME BEFORE ORIGIN TIME
c
                  write(stout,1001) 'Arrival time before origin time'
                  goto 999
                else
                  obst = evtmsc(i) - t
                endif
c
c             MINUTES DIFFER
c
              else if (evtmmi(i) .lt. smin) then
c
c               ERROR: ARRIVAL TIME BEFORE ORIGIN TIME
c
                write(stout,1001) 'Arrival time before origin time'
                goto 999
              else
                obst = (evtmmi(i) - smin)*60 + evtmsc(i) - t
              endif
c
c           HOURS DIFFER
c
            else if (evtmhr(i) .lt. shr) then
c
c             ERROR: ARRIVAL TIME BEFORE ORIGIN TIME
c
              write(stout,1001) 'Arrival time before origin time'
              goto 999
            else
              obst = (evtmhr(i) - shr)*3600 + (evtmmi(i) - smin)*60
     &               + evtmsc(i) - t
            endif
c
c         DAYS DIFFER
c
          else if (evtmdy(i) .lt. sday) then
c
c           ERROR: ARRIVAL TIME BEFORE ORIGIN TIME
c
            write(stout,1001) 'Arrival time before origin time'
            goto 999
          else
            obst = (evtmdy(i) - sday)*24*3600 + (evtmhr(i) - shr)*3600
     &             + (evtmmi(i) - smin)*60 + evtmsc(i) - t
          endif
c
c       MONTHS DIFFER
c
        else if (moevt .lt. mosrc) then
c
c         ERROR: ARRIVAL TIME BEFORE ORIGIN TIME
c
          write(stout,1001) 'Arrival time before origin time'
          goto 999
        else if ((moevt - mosrc) .ne. 1) then
c
c         ERROR: IDIOTIC INPUT WITH MORE THAN A MONTH BETWEEN
c         ORIGIN AND ARRIVAL TIMES
c
          write(stout,1001) 'More than a month between origin and arriva
     &l times'
          goto 999
        else if ((mosrc .eq. 2) .and. (mod(syr,4) .eq. 0)) then
c
c         ITS A LEAP YEAR
c
          obst = (moevt - mosrc)*29*24*3600
     &     + (evtmdy(i) - sday)*24*3600 + (evtmhr(i) - shr)*3600
     &     + (evtmmi(i) - smin)*60 + evtmsc(i) - t
        else
          obst = (moevt - mosrc)*MONLEN(mosrc)*24*3600
     &     + (evtmdy(i) - sday)*24*3600 + (evtmhr(i) - shr)*3600
     &     + (evtmmi(i) - smin)*60 + evtmsc(i) - t
        endif
c
c     YEARS DIFFER
c
      else if (evtmyr(i) .lt. syr) then
c
c       ERROR: ARRIVAL TIME BEFORE ORIGIN TIME
c
        write(stout,1001) 'Arrival time before origin time'
        goto 999
      else if (((evtmyr(i) - syr) .ne. 1) .or. (mosrc .ne. 12)
     & .or. (moevt .ne. 1)) then
c
c       ERROR: IDIOTIC INPUT AS MORE THAN A MONTH BETWEEN ORIGIN AND
c       ARRIVAL TIMES
c
        write(stout,1001) 'More than a month between origin and arrival
     &times'
        goto 999
      else
        obst = MONLEN(12)*24*3600 + (evtmdy(i) - sday)*24*3600
     &   + (evtmhr(i) -shr)*3600 + (evtmmi(i) - smin)*60
     &   + evtmsc(i) - t
      endif
c
c     WRITE RESULTS TO OUTPUT FILE
c
c      write(stout,1002)
c      write(stout,1003) 'ORIGIN', syr, smo, sday, shr, smin, t
c      write(stout,1003) 'ARRIVAL', evtmyr(i), evtmmo(i), evtmdy(i),
c     & evtmhr(i), evtmmi(i), evtmsc(i)
c      write(stout,1004) 'OBSERVED TRAVEL TIME:', obst
c
*----------------------------------------------------------------------*------*
c
c     FORMAT STATEMENTS
c
1000  format (1X, A)
1001  format (1X, 'ERROR: ', A)
1002  format (1X, T3, 'TIME', T11, 'YEAR', T16, 'MONTH', T22, 'DAY',
     & T26, 'HOUR', T31, 'MIN', T35, 'SEC')
1003  format (1X, T3, A, T11, I4, T16, A, T22, I2, T26, I2, T31, I2,
     & T35, F5.2)
1004  format (1X, T3, A, T25, F10.2)
c
*----------------------------------------------------------------------*------*
c
      return
c
999   write(stout,1002)
      write(stout,1003) 'ORIGIN', syr, smo, sday, shr, smin, t
      write(stout,1003) 'ARRIVAL', evtmyr(i), evtmmo(i), evtmdy(i),
     & evtmhr(i), evtmmi(i), evtmsc(i)
c
c      write(res,1002)
c      write(res,1003) 'ORIGIN', syr, smo, sday, shr, smin, t
c      write(res,1003) 'ARRIVAL', evtmyr(i), evtmmo(i), evtmdy(i),
c     & evtmhr(i), evtmmi(i), evtmsc(i)
c
      return
      end
*--------------------------------------------------------------------- + ---- +
c                                                                      month
      integer function month (cmon)
c
c     AUTHOR:  Genet M. Edmondson  RSES, ANU
c     DATE:    February 1987
c     PURPOSE:
c     * converts a string of 4 characters representing roman numerals
c       from 1 to 12 into an integer
c     * assumes the input is one of I, II, III, IV, V, VI, VII, VIII,
c       IX, X, XI, XII
c
*----------------------------------------------------------------------*------*
c
c     PARAMETERS
c
      character*4   cmon
c
c     cmon  character string for roman numerals
c
*----------------------------------------------------------------------*------*
c
c
      if (cmon .eq. 'I') then
        month = 1
      else if (cmon .eq. 'II') then
        month = 2
      else if (cmon .eq. 'III') then
        month = 3
      else if (cmon .eq. 'IV') then
        month = 4
      else if (cmon .eq. 'V') then
        month = 5
      else if (cmon .eq. 'VI') then
        month = 6
      else if (cmon .eq. 'VII') then
        month = 7
      else if (cmon .eq. 'VIII') then
        month = 8
      else if (cmon .eq. 'IX') then
        month = 9
      else if (cmon .eq. 'X') then
        month = 10
      else if (cmon .eq. 'XI') then
        month = 11
      else if (cmon .eq. 'XII') then
        month = 12
      endif
c
*----------------------------------------------------------------------*------*
c
      return
      end
