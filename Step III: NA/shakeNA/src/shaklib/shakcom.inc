*----------------------------------------------------------------------*------*
c
c     CONSTANTS FOR COMMON BLOCK DECLARATIONS - Program SHAKE etc
c     (kept in shakcom.inc)
c
      integer MAXST, MAXEVT, MAXTT
c
      parameter (MAXST=1000, MAXEVT=1000, MAXTT=60 )
c
c     MAXST   maximum number of stations used
c     MAXEVT  maximum number of readings used
c     MAXTT   maximum number of calculated phase times
c             at a given range (for trtm)
c
*----------------------------------------------------------------------*------*
c
c     COMMON BLOCK DECLARATIONS
c     (kept in common.inc)
c     (needed in shake.ftn, input.ftn, apploc.ftn, obst.ftn)
c
      character*4   source, smo
      character*48  region
c
      integer   syr, sday, shr, smin
c
      real    ssec
      real    slat, slng, sdep
      real    sersec, serlat, serlng, serdep
c
      COMMON / CQAKID / source, region, syr, smo, sday, shr, smin, ssec,
     & sersec, slat, serlat, slng, serlng, sdep, serdep
c
c     source  name of source of information about quake location (ISC etc.)
c     region  name of area where quake occurred (or a title for quake)
c     syr, smo, sday, shr, smin, ssec
c             year, month, day, hour, minute and second when quake
c             happened according to source
c       syr   4 digit integer
c       sday, shr, smin  2 digit integers
c       smo   roman numerals
c       ssec  real   number with <= 3 decimal places
c     slat, slng, sdep
c             latitude, longitude and depth of quake according to source
c             real  s with <= 3 decimal places
c     sersec, serlat, serlng, serdep
c             error in second, lat., long., and depth given by source
c             real  s with <= 3 decimal places
c
*----------------------------------------------------------------------*------*
c
      character*4 evstcd(MAXEVT), evtmmo(MAXEVT)
      character*8 evphcd(MAXEVT)
c
      integer   evtmyr(MAXEVT), evtmdy(MAXEVT)
      integer   evtmhr(MAXEVT), evtmmi(MAXEVT)
c
      real    evtmsc(MAXEVT), evtmer(MAXEVT)
      real    evtazi(MAXEVT), evtaer(MAXEVT)
      real    evtslo(MAXEVT), evtser(MAXEVT)
c
      COMMON / CEVENT / evphcd, evstcd, evtmyr, evtmmo, evtmdy, 
     &                  evtmhr, evtmmi, evtmsc, evtmer
      COMMON / CAZSLO / evtazi, evtaer, evtslo, evtser
c
c     evstcd(i)   code of station where i-th recording was made
c     evphcd(i)   phase code of arrival for i-th reading
c     evtmyr(i)   year of i-th reading
c     evtmmo(i)   month of i-th reading
c     evtmdy(i)   day of i-th reading
c     evtmhr(i)   hour of i-th reading
c     evtmmi(i)   minute of i-th reading
c     evtmsc(i)   second of i-th reading
c     evtmer(i)   error (seconds) in i-th reading
c                 (currently set to 0.2; used to calculate weighting
c                 function sigma for stations and readings)
c
c     evtazi(i)   azimuth for i-th reading (degrees)
c     evtaer(i)   error in aziumuth (degrees)
c     evtslo(i)   slowness for i-th reading (s/km)
c     evtser(i)   error in slowness (s/km)
c
*----------------------------------------------------------------------*------*
c
      character*4 stcode(MAXST)
c
      real    stlat(MAXST), stlong(MAXST), stelev(MAXST), stwt(MAXST)
c
      COMMON / CSTATN / stcode, stlat, stlong, stelev, stwt
c
c     stcode(i)   code for i-th station
c     stlat(i)    latitude of i-th station
c     stlong(i)   longitude of i-th station
c     stelev(i)   elevation of i-th station
c     stwt(i)     relative station error for i-th station
c                 indicates accuracy of station readings
c                 typically between 1.0 (fairly accurate) and 2.0
c
*----------------------------------------------------------------------*------*
      character*4 refst
c
      integer   stind(MAXST), evind(MAXEVT)
      integer   numevt, numtim, numazi, numslo
c
      real    epslat, epslng, epsdep, epstim, sigma(MAXEVT), rfstcr
c
      COMMON / CRUN / epslat, epslng, epsdep, epstim, sigma, rfstcr,
     & stind, evind, numevt, numtim, numazi, numslo, refst
c
c     epslat, epslng, epsdep, epstim
c       accuracy required in location (latitude, longitude, depth, and
c       time) of quake
c     sigma    weighting function used in calculating residual statistic
c     evstcd(stind(i)) = stcode(i)
c     stcode(evind(i)) = evstcd(i)
c     numevt   number of event readings
c       numtim   number of time readings
c       numazi   number of azimuth readings
c       numslo   number of slowness readings
c     refst    code  of reference station
c     rfstcr   correction to observed arrival time at reference station
c
*----------------------------------------------------------------------*------*
c
      real calt(MAXEVT,7,7,7), calp(MAXEVT,7,7,7)
      real delt(MAXEVT,7,7), cazi(MAXEVT,7,7), bazi(MAXEVT,7,7)
c
      COMMON /COTIM/ calt, calp, delt, cazi, bazi
c
c     common used to store calculated travel times, slowness
c     distances and azimuths for current search grid
c
c     calt   calculated travel time on search grid for each phase event
c     calp   calculated slowness on search grid for each phase event
c     delt   angular distance from stations to epicentral grid
c     cazi   azimuth of wavefront at an array station
c
*----------------------------------------------------------------------*------*
      logical prnt
c
      character*8 phcd(MAXTT),phlst(10)
c
      dimension tt(MAXTT),dtdd(MAXTT),dtdh(MAXTT),dddp(MAXTT)
      dimension usrc(2)
c
      COMMON /PRTFLC/ prnt(3)
      COMMON /TRAVEL/ tt,dtdd,dtdh,dddp,usrc,phcd,phlst
c
c     commons used in calculation of travel times using Buland tau-p
c     procedure
c     prnt     printflags (.true. for extensive debug output)
c     For a given DELTA:
c     tt       calculated travel times for various phases 
c     dtdd     slowness (s/deg)
c     dtdh     depth derivative
c     dddp     d2t/dd2 - amplitude measure
c
*----------------------------------------------------------------------*------*
*
