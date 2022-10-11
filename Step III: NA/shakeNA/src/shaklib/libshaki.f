*-----        Input & Output routines         ------------------------ + ---- +
c                                                                      input
      subroutine input (stout, evt, stat, run, res)
c
c     AUTHOR:  Genet M. Edmondson  RSES, ANU
c     DATE:    January 1987
c     MODIFIED: Brian L.N. Kennett  March 1990
c               B.L.N. Kennett May 1991
c               to include slowness, azimuth data
c
c     PURPOSE:
c     * obtains information for program shake about quake location,
c       stations and run parameters
c     * calculation of indices assumes that there is only one reading
c       from each station
c
*----------------------------------------------------------------------*------*
c
c     PARAMETERS
c
      integer stout, evt, stat, run, res
c
c     stout   unit number of standard output file
c     evt     unit number of file with event information
c     stat    unit number of file with information about stations
c     run     unit number of file with run information
c     res     unit number of file to store results in
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
      character*1   slatdr, slngdr
c
      integer i, j
      integer numst
c
      logical   found
c
c     slatdr  direction (N/S) of latitude
c     slngdr  direction (E/W) of longitude
c     i, j    loop counters etc.
c     numst   number of stations that data is read in about
c     found   used to indicate if matching codes have been found when
c             calculating indices
c
*----------------------------------------------------------------------*------*
c
c
c     INITIALISE VARIABLES
c
      DATA (stind(i), i = 1, MAXST) / MAXST * 0.0 /
      DATA (evind(i), i = 1, MAXEVT) / MAXEVT * 0.0 /
c
c     OBTAIN INFORMATION IDENTIFYING QUAKE
c
      read(evt,1024,ERR=999) region
      read(evt,1024,ERR=999) source
      read(evt,1021,ERR=999) syr, smo, sday, shr, smin, ssec, sersec
      read(evt,1022,ERR=999) slat, serlat
      read(evt,1022,ERR=999) slng, serlng
      read(evt,1022,ERR=999) sdep, serdep
c
c     ECHO QUAKE IDENTIFYING INFORMATION TO OUTPUT FILE
c
      write(stout,1000)
      write(stout,1002) 'REGION', region
      write(stout,1002) 'DATA', source
      write(stout,1003) 'DATE'
      write(stout,1004) 'Year', syr
      write(stout,1005) 'Month', smo
      write(stout,1006) 'Day', sday
      write(stout,1006) 'Hour', shr
      write(stout,1006) 'Minute', smin
      write(stout,1007) 'Second', ssec, sersec
      write(stout,1003) 'LOCATION'
      write(stout,1007) 'Latitude', slat, serlat
      write(stout,1007) 'Longitude', slng, serlng
      write(stout,1007) 'Depth', sdep, serdep
c
c     WRITE QUAKE IDENTIFYING INFORMATION TO RESULTS FILE
c
      write(res,1000)
      write(res,1002) 'REGION', region
      write(res,1002) 'DATA', source
      write(res,1003) 'DATE'
      write(res,1004) 'Year', syr
      write(res,1005) 'Month', smo
      write(res,1006) 'Day', sday
      write(res,1006) 'Hour', shr
      write(res,1006) 'Minute', smin
      write(res,1007) 'Second', ssec, sersec
      write(res,1003) 'LOCATION'
      write(res,1007) 'Latitude', slat, serlat
      write(res,1007) 'Longitude', slng, serlng
      write(res,1007) 'Depth', sdep, serdep
c
c     OBTAIN EVENT DATA
c
      write(stout,1003) 'EVENTS read'
c     write(stout,1011)
c
      write(res,1003) 'EVENTS read'
      write(res,1031)
c                                     read travel times
      i = 1
      numevt = 0
      numtim = 0
      numazi = 0
      numslo = 0
c
210   continue
c                                        modified to read phase code
        read(evt,1030,END=220,ERR=999) evstcd(i), evphcd(i),
     &             evtmhr(i), evtmmi(i), evtmsc(i), evtmer(i)
        if(evstcd(i).eq.'azim') go to 213
        if(evstcd(i).eq.'slow') go to 217
        if(evstcd(i)(1:3).eq.'end') go to 220
        evtmyr(i) = syr
        evtmmo(i) = smo
        evtmdy(i) = sday
        numtim = numtim + 1
        numevt = numevt + 1
c
c       write event data to auxiliary results file
c
        write(res,1032) numevt, evstcd(i), evphcd(i),
     &             evtmhr(i), evtmmi(i), evtmsc(i), evtmer(i)
c
        i = i + 1
        goto 210
c                                   read azimuth information
213   continue
        read(evt,1020,END=220,ERR=999) evstcd(i), evphcd(i),
     &             evtazi(i),evtaer(i)
        if(evstcd(i).eq.'slow') go to 217
        if(evstcd(i)(1:3).eq.'end') go to 220
        numazi = numazi + 1
        numevt = numevt + 1
        write(res,1033) numevt, evstcd(i), evphcd(i),
     &              evtazi(i), evtaer(i)
        i= i+1
        goto 213
c                                   read slowness information
217   continue
        read(evt,1020,END=220,ERR=999) evstcd(i), evphcd(i),
     &             evtslo(i),evtser(i)
        if(evstcd(i)(1:3).eq.'end') go to 220
        numslo = numslo + 1
        numevt = numevt + 1
        write(res,1033) numevt, evstcd(i), evphcd(i),
     &              evtslo(i), evtser(i)
        i= i+1
        goto 217
c
220     continue
      write(stout,*) 'num evt = ', numevt
      write(stout,*) 'numtim,numazi,numslo:',numtim,numazi,numslo
      write(res,*) 'num evt = ', numevt
      write(res,*) 'numtim,numazi,numslo:',numtim,numazi,numslo
c
c     OBTAIN RUN INFORMATION
c
      read(run,*,ERR=999) epslat
      read(run,*,ERR=999) epslng
      read(run,*,ERR=999) epsdep
      read(run,*,ERR=999) epstim
c
c     ECHO RUN INFORMATION TO OUTPUT
c
      write(stout,1001) 'ACCURACY REQUIRED IN SOLUTION:'
      write(stout,1015) 'Latitude', epslat
      write(stout,1015) 'Longitude', epslng
      write(stout,1015) 'Depth', epsdep
      write(stout,1015) 'Time', epstim
c
c     WRITE RUN INFORMATION TO AUXILIARY RESULTS FILE
c
      write(res,1001) 'ACCURACY REQUIRED IN SOLUTION:'
      write(res,1015) 'Latitude', epslat
      write(res,1015) 'Longitude', epslng
      write(res,1015) 'Depth', epsdep
      write(res,1015) 'Time', epstim
c
c     OBTAIN STATION DATA
c
      write(stout,1003) 'STATIONS read'
c     write(stout,1013)
c
      write(res,1003) 'STATIONS read'
      write(res,1013)
c
      read(stat,1001)
      read(stat,1001)
      i = 1
310   continue
        read(stat,1023,END=320,ERR=999) stcode(i), stlat(i), stlong(i),
     &  stelev(i), stwt(i)
        if(stwt(i).eq.0.0) stwt(i) = 1.0
c
c       ECHO STATION DATA & WRITE TO RESULTS FILE
c
c       write(stout,1014) i, stcode(i), stlat(i), stlong(i), stelev(i),
c    &  stwt(i)
        write(res,1014) i, stcode(i), stlat(i), stlong(i), stelev(i),
     &  stwt(i)
c
        i = i + 1
        goto 310
320     continue
      numst = i - 1
      write(stout,*) 'num stat = ', numst
c
c     CALCULATE CROSS-REFERENCE INDICES FOR STATION & EVENT ARRAYS
c
      do 400, j = 1, numevt
c       found = .false.
        write(res,*) evstcd(j)
        do 410 i= 1, numst
          if (evstcd(j) .eq. stcode(i)) then
c           stind(i) = j
            evind(j) = i
          endif
410     continue
400   continue
c
c     Indices of stations, phases events written to unit RES
c     write(res,1001) 'Station indices'
c     do 420, i = 1, numst
c       write(res,1016) i, stind(i)
c420     continue
      write(res,1001) 'Event indices'
      do 430, i = 1, numevt
        write(res,1016) i, evind(i)
430     continue
c
c     CALCULATE WEIGHTING VALUES - sigma
c                            (factor of 5. for teleseisms)
c
      numtaz = numtim+numazi
      do 500, i = 1, numevt
        if(i.le.numtim) then
          sigma(i) = 5.* evtmer(i) * stwt(evind(i))
          write(res,*) i,sigma(i)
        elseif(i.gt.numtim .and. i.le.numtaz) then
          sigma(i) = evtaer(i)
          write(res,*) i,sigma(i)
        elseif(i.gt.numtaz) then
          sigma(i) = evtser(i)
          write(res,*) i,sigma(i)
        endif
500     continue
c
c     ECHO WEIGHTING VALUES
c
      write(res,1001) 'sigma values'
      write(res,*) (sigma(i), i = 1, numtim)
      write(res,*) (sigma(i), i = numtim+1,numtaz)
      write(res,*) (sigma(i), i = numtaz+1,numevt)
c
*----------------------------------------------------------------------*------*
c
c     FORMAT STATEMENTS
c
1000  format ()
1001  format (1X, A)
1002  format (1X, A, ':', T10, A)
1003  format (1X, A, ':')
1004  format (1X, T3, A, T13, I4)
1005  format (1X, T3, A, T13, A)
1006  format (1X, T3, A, T13, I2)
1007  format (1X, T3, A, T13, F7.3, T22, 'Error', T29, F7.3)
1011  format (1X, T9, 'CODE', T14, 'YEAR', T19, 'MONTH',
     & T25, 'DAY', T29, 'HOUR', T34, 'MIN', T38, 'SEC', T44, 'ERROR')
1012  format (1X, T3, I3, T9, A, T14, I4, T19, A, T25, I2, T29, I2, T34,
     & I2, T38, F5.2, T44, F4.2)
1013  format (1X, T9, 'CODE', T14, 'LAT', T25, 'LONG', T36, 'ELEV',
     & T43, 'WEIGHT')
1014  format (1X, T3, I3, T9, A, T14, F10.5, T25, F10.5, T36, F6.3,
     & T43, F5.3)
1015  format (1X, T3, A, T14, F7.3)
1016  format (1X, T3, I3, T7, I3)
1020  format (T1, A, T6, A, T16, F10.3, T30, F10.3)
1021  format (I4, T7, A, T12, I2, T15, I2, T18, I2, T21, F7.3, T29,
     & F7.3)
1022  format (T1, F9.3, T10, F9.3)
1023  format (T1, A, T9, F15.5, T24, F15.5, T39, F15.3, T54, F15.3)
1024  format (A)
1030  format (T1, A, T6,  A, T16, I2, T19, I2, T22,
     & F5.2, T28, F4.2)
1031  format (1X, T9, 'CODE', T14, 'PHASE', 
     &      T25, 'HOUR', T30, 'MIN', T35, 'SEC', T42, 'ERROR')
1032  format (1X, T3, I3, T9, A, T14, A, T25, I2, T30,
     & I2, T35, F5.2, T42, F4.2)
1033  format (1X, T3, I3, T9, A, T14, A, T25, F10.3, T42, F10.3)
c
*----------------------------------------------------------------------*------*
c
      return
c
c     ERROR OCCURRED WHILE READING IN INFORMATION
c
999   write(stout,1001) 'ERROR: reading input (subroutine input)'
      stop
c
      return
      end
*--------------------------------------------------------------------- + ---- +
c                                                                     disply
      subroutine disply (stout, res, jgl,
     &                   LA, LN, D, S, lt, lg, dp, sc)
c
c     AUTHOR:  Genet M. Edmondson  RSES, ANU
c     DATE:    January 1987 
c     MODIFIED: B.L.N. Kennett March 1990  
c               to allow for changes to residual calculation
c               B.L.N. Kennett May 1991
c               to include slowness, azimuth data
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
      real    LA(*), LN(*), D(*), S(*)
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
      real     resqkp, ressrc, resdif, residb
      real     del
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
c     residf  station residual for quake location
c     residb  station residual for approx source 
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
      write(stout,1000) 'BOUNDS'
      write(stout,1001) source
      write(stout,1002) 'DATE:'
      write(stout,1003) 'Year', syr, syr
      write(stout,1004) 'Month', smo, smo
      write(stout,1005) 'Day', sday, sday
      write(stout,1005) 'Hour', shr, shr
      write(stout,1005) 'Minute', smin, smin
      write(stout,1007) 'Second', S(1), S(2), ssec - sersec,
     & ssec + sersec
      write(stout,1002) 'LOCATION:'
      write(stout,1007) 'Latitude', LA(1), LA(2), slat - serlat,
     & slat + serlat
      write(stout,1007) 'Longitude', LN(1), LN(2), slng - serlng,
     & slng + serlng
      write(stout, 1007) 'Depth', D(1), D(2), sdep - serdep,
     & sdep + serdep
c                                    Find residual statistics
c                                    a) best grid point
      qla = 0.5*(LA(1)+LA(2))
      qlg = 0.5*(LN(1)+LN(2))
      qdp = 0.5*(D(1)+D(2))
      qt = 0.5*(S(1)+S(2))
      resmid = cstaq(stout, jgl, qt, qla, qlg, qdp)
c                                    b) quadratic estimate
      resqkp = cstaq(stout, jgl, sc, lt, lg, dp)
c                                    c) approx location
      ressrc = cstaq(stout, jgl, ssec, slat, slng, sdep)
c                                    use best grid point if
c                                    quadratic estimate not an
c                                    improvement
      if(resmid.lt.resqkp) then
        lt = qla
        lg = qlg
        dp = qdp
        sc = qt
        resqkp = resmid
        write (stout,*) "GRID solution used"
      endif
       if(jgl.ne.1) then
           cstarf = cstaq(stout, 1, sc, lt, lg, dp)
           cstar0 = cstaq(stout, 1, ssec, slat, slng, sdep)
       endif
            
c
      write(stout,1000) 'SOLUTION'
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
      if(jgl.ne.1) then
        write(stout,1009) 'Quad Res:', cstarf, cstar0
      endif        
c
c     CALCULATE & ECHO ARRIVAL TIMES AT STATIONS FOR ORIGIN
c     AT CALCULATED SOLUTION
c
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
        do 23 m=1,n
 23      continue
        if(nc.eq.1.and.evphcd(i)(1:1).eq.'P') then
          calct = tt(1)
          calcp = dtdd(1)/111.19
          write(res,*) i,' fast P'
        else
          do 21 l=1,n
            if(evphcd(i)(1:nc).eq.phcd(l)(1:nc)) then
              calct = tt(l)
              calcp = dtdd(l)/111.19
              write(res,*) i,nc,'  ',evphcd(i)
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
     &          evtmsc(i),residf,carhr,carmin,carsec,delta,bazim
        elseif (i.gt.numtim .and. i.le.numtaz) then
          write(stout,1017) evstcd(i),evphcd(i),
     &         evtazi(i),cazim,delta,bazim
        elseif (i.gt.numtaz) then
          write(stout,1017) evstcd(i),evphcd(i),
     &         evtslo(i),calcp,delta,bazim
        endif
c
10      continue
c
*----------------------------------------------------------------------*------*
c
c     FORMAT STATEMENTS
c
1000  format (1X, T15, 'HYPOCENTRE ', A)
1001  format (1X, T15, 'Quake program', T40, A)
1002  format (1X, T3, A)
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
     & T53, 'csec', T60, 'delta', T70, 'baz')
1016  format (1X, T2, A, T7, A, T15, I2, T20, I2, T25,
     & F5.2, T33,  F5.2, T43, I2, T48, I2, T53, F5.2, 
     & T60, F5.2, T70, F5.1)
1017  format (1X, T2, A, T7, A, T20, F10.3, T48, F10.3, 
     & T60, F5.2, T70, F5.1)
c
*----------------------------------------------------------------------*------*
c
      return
      end
*-----        Utility functions         ------------------------------ + ---- +
c                                                                      ready
      logical function ready (stinp, stout)
c
c     AUTHOR:  Genet M. Edmondson  RSES, ANU
c     DATE:    December 1986
c     PURPOSE:
c     * This routine reads a character which is one of y, Y, n and N.
c       If the character read is y or Y the function is true,
c       otherwise it returns false.
c
*----------------------------------------------------------------------*------*
c
c     PARAMETERS
c
      integer stinp, stout
c
c     stinp   unit number of standard input file
c     stout   unit number of standard output file
c
*----------------------------------------------------------------------*------*
c
c     VARIABLES
c
      character*1 c
c
c     c  temporary variable
c
*----------------------------------------------------------------------*------*
c
c
10    continue
        read(stinp, 1010) c
        if (.not.((c .eq. 'y') .or. (c .eq. 'Y') .or. (c .eq. 'n') .or.
     &    (c .eq. 'N'))) then
          write(stout,1000) 'ERROR: incorrect input - re-enter'
          goto 10
        endif
c
      if ((c .eq. 'y') .or. (c .eq. 'Y')) then
        ready = .true.
      else
        ready = .false.
      endif
c
*----------------------------------------------------------------------*------*
c
c     FORMAT STATEMENTS
c
1000  format (1X, A)
1010  format (A)
c
*----------------------------------------------------------------------*------*
c
      return
      end
*--------------------------------------------------------------------- + ---- +
c                                                                      ltou
      subroutine ltou (c)
c
c     AUTHOR:  Genet M. Edmondson  RSES, ANU
c     DATE:    December 1986
c     PURPOSE:
c     * This routine returns the argument as an upper case letter.
c     * This routine is intended for use only where the upper and lower
c       case letters are not interspersed in the collating sequence
c       (the sequence of characters arranged according to the value
c       returned by ichar).
c
*----------------------------------------------------------------------*------*
c
c     PARAMETERS
c
      character*1 c
c
c     c is the character to be converted to upper case
c       (if it is not alreaddy upper case)
c
*----------------------------------------------------------------------*------*
c
c     VARIABLES
c
      integer numc
      integer smalla, biga, smallz
c
c     numc is the ordinal value of c
c     smalla is the ordinal value of 'a'
c     biga is the ordinal value of 'A'
c     smallz is the ordinal value of 'z'
c
*----------------------------------------------------------------------*------*
c
c
c     INITIALISE VARIABLES
c
      smalla = ichar('a')
      biga = ichar('A')
      smallz = ichar('z')
      numc = ichar(c)
c
      if ((numc .ge. smalla) .and. (numc .le. smallz)) then
c
c       c IS LOWER CASE SO CONVERT TO UPPER
c
        c = char(biga - smalla + numc)
      endif
c
*----------------------------------------------------------------------*------*
c                     
      return
      end
*--------------------------------------------------------------------- + ---- +
c                                                                      charcn
      integer function charcn (c)
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
c     integer charcn
c
*----------------------------------------------------------------------*------*
c
c
      do 10 k=1,8
        if(c(k:k).eq.' ') then
            charcn = k-1
            return
        endif
 10   continue
      charcn = 8 
c
*----------------------------------------------------------------------*------*
      return
      end
c
*--------------------------------------------------------------------- + ---- +
c                                                                      todec
      real function todec (d, m, s)
c
c     AUTHOR:  Genet M. Edmondson  RSES, ANU
c     DATE:    January 1987
c     PURPOSE:
c     * converts degrees, minutes and seconds to decimal degrees
c
*----------------------------------------------------------------------*------*
c
c     PARAMETERS
c
      real  d, m, s
c
c     d degrees
c     m minutes
c     s seconds
c
*----------------------------------------------------------------------*------*
c
c
      todec = d + m/60 + s/3600
c
*----------------------------------------------------------------------*------*
c
      return
      end
*--------------------------------------------------------------------- + ---- +
c                                                                      setup
      subroutine setup (stinp, stout, store, datfil)
c
c     AUTHOR:  Genet M. Edmondson  RSES, ANU
c     DATE:    January 1987
c     PURPOSE:
c     * sets up a file for storing data
c     * returns the unit number of the file in 'store' and
c       its name in 'datfil'
c
c     Subroutines and functions required:
c       ready - returns true if either 'y' or 'Y' is read
c
*----------------------------------------------------------------------*------*
c
c     PARAMETERS
c
      integer stinp, stout, store
c
      character*20  datfil
c
c     stinp   unit number of standard input
c     stout   unit number of standard output
c     store   unit number of storage file
c     datfil  name of storage file
c
*----------------------------------------------------------------------*------*
c
c     FUNCTIONS
c
      logical ready
c
c     ready  returns true if 'y' or 'Y' is read
c
*----------------------------------------------------------------------*------*
c
c
      write(stout,1010)
c
10    continue
        write(stout,1000) 'Enter name of file for data storage (<20 char
     &acters):'
        read(stinp,1011,ERR=10) datfil
        write(stout,1001) 'Name of data file is ', datfil
        write(stout,1003)
        if (ready(stinp,stout)) goto 10
c
20    continue
        write(stout,1002) 'Enter unit number of file ', datfil, ':'
        read(stinp,*,ERR=20) store
        write(stout,1004) 'Unit number of file ', datfil, ' is ', store
        write(stout,1003)
        if (ready(stinp,stout)) goto 20
c                                     Masscomp form 
c     open(UNIT=store, FILE=datfil, form='noprint')
c                                     generic form 
      open(UNIT=store, FILE=datfil)
      rewind(store)
c
*----------------------------------------------------------------------*------*
c
c     FORMAT STATEMENTS
c
1000  format (1X, A)
1001  format (1X, A, A)
1002  format (1X, A, A, A)
1003  format (1X, 'Do you want to change this (y/n)?')
1004  format (1X, A, A, A, I2)
1010  format ()
1011  format (A)
c
*----------------------------------------------------------------------*------*
c
      return
      end
