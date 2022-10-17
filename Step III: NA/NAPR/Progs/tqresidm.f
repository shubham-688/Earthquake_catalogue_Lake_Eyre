c      program tqresid
       subroutine tqresidm(cname)
c                            comparison of empirical and
c                            model times for surface focus
c      the primary output from tqresid is a set of individual phase misfits
c
c ph #						weight
c
c	Ja01
c(1)	P              0.013    73     425.5	[5.0]
c(2)	PP             0.266   128     419.4	[1.5]
c(3)	PcP            0.108    45     325.2	[2.0]
c(4)	S              0.047    56     603.2	[3.0]
c(5)    SS             1.006    95     102.1	[1.0]
c(6)	ScS            0.055    40     139.2	[1.0]
c(7)	SP             0.774    46     130.7	[1.0]
c(8)	ScP            0.101    38     179.6	[2.0]
c(9)	PKPdf          0.123    75     542.7	[4.0]
c(10)	PKPbc          0.495     3     436.0	[4.0]
c(11)	PKPab          0.510    22     448.4	[4.0]
c(12)	PKKPbc         0.285    40      68.6	[1.5]
c(13)	PKKPab         0.081    12      83.1	[1.5]
c(14)	SKSac          0.046    33     301.2	[3.0]
c(15)	SKKSac         0.210   113      40.7	[1.5]
c(16)	SKPdf          0.548    62      46.8	[1.0]
c(17)	SKPbc          0.355     7     131.0	[1.0]
c(18)	P'P'           1.541    37      45.7	[1.0]
c
c	P+PcP          0.121
c	S+ScS          0.102
c	PKP            1.129
c	SKS+SKKS+SKP   1.159
c
c	PP+SS+SP       2.046
c	PKP+PKKP       1.495
c	SP+ScP+SKP     1.778
c
c	P+PcP+PKP      1.249
c	S+ScS+ScP      0.203
c	major P,S      1.452
c
c	mantle ph      2.370
c	core ph.       4.194
c
c	all ISC        6.564
c	rel PKP              1.03 0.25
c
c                    from these entries for the individual phase
c                    we can build different measures of the
c                    travel time fit which can be weighted by
c                    phase importance - indicated in square brackets above]
c
c	"P1"
c		P + PcP + PKP			{1,3,9,10,11}
c	"S1"
c		S + ScS + SKS			{4,6,14}
c	"M1"
c		P + PcP + S + ScS + ScP 	{1,3,4,6,8}
c	"M2"
c		P + PcP + PP + S + ScS + SS + ScP + SP
c						{1,3,2,4,6,5,8,7}
c	"CP"
c		PKP + PKKP + P'P'		{9,10,11,12,13,18}
c	"CS"
c		SKS + SKKS + SKP		{14,15,16,17}
c
c       "A1"    
c               P + PcP + PKP + S + ScS + SKS + ScP
c                                              {1,3,9,10,11,4,6,14,8}
c       "A2"    
c               P + PcP + PKP + S + ScS + SKS + ScP +
c               PKKP + SKKS + SKP
c                                               {1,3,9,10,11,4,6,14,8
c                                               {12,13,15,16,17}
c       "AL"
c		P + PcP + PP + S + ScS + SS + ScP + SP+
c		PKP + PKKP + P'P'+ SKS + SKKS + SKP
c						{1-18}
c
      save
      parameter (max=60)
      logical prnt(3)
      real psi(20),aht(20),wt(20)
      character*8 phcd(max),phlst(10),ophcd
      character*8 cbr(20)
      character*8 modnam,cname
      integer nps(20)
      dimension tt(max),dtdd(max),dtdh(max),dddp(max)
      dimension usrc(2)
      data in/1/,phlst(1)/'all'/,prnt(3)/.false./
      data cbr /"P","PP","PcP","S","SS","ScS","SP","ScP",
     &          "PKPdf","PKPbc","PKPab","PKKPbc","PKKPab",
     &          "SKSac","SKKSac","SKPdf","SKPbc","P'P'",
     &          "Sg","Sb"/
       data wt/5.0,1.5,2.0,3.0,1.0,1.0,1.0,2.0,4.0,4.0,4.0,
     &          1.5,1.5,3.0,1.5,1.0,1.0,1.0,0.0,0.0/
c
      modnam=cname
c
c
      write(6,*) 'This routine for calculating residuals for'
      write(6,*) 'uses a set of precalculated tau-p tables'
      write(6,*) modnam,'.hed  ',modnam,'.tbl'
      write(6,*)
      prnt(1) = .false.
      prnt(2) = .false.
      call assign(10,2,'ttim1.lis',9)
      write(6,*) "ttimes-assign: ttim1.lis"
      call tabin(in,modnam)
      call brnset(1,phlst,prnt)
      nlen = lenstr(modnam)
c                                    choose source depth (0.)
      zs = 0.
      call depset(zs,usrc)
c
      nemp = 1
      nout = 2
      open(nemp,file='surftim.dat')
      open(nout,file='P'//modnam(1:nlen)//'.out')
      write(6,*) 'Output of residuals summary to file'
      write(6,*) 'P',modnam(1:nlen),'.out '
      write(nout,fmt='(a)') modnam
c                                    set information by branch
      do 100 k=1,20
        psi(k) = 0.0
        aht(k) = 0.0
        nps(k) = 0
 100   continue
c                                    loop on delta
      nd = 1
 10   read(nemp,*,end=13) delta,nph
**      write(6,*) delta,nph
      call trtm(delta,max,n,tt,dtdd,dtdh,dddp,phcd)
c
      do 30 i=1,nph
        read(nemp,fmt='(a,2x,2f10.2,f10.0)') ophcd,obst,ovar,ohit
        nc = icharcn(ophcd)
*        write(6,*) ophcd,nc
        do 21 l=1,n
*           write(6,*) ophcd(1:nc),phcd(l)(1:nc),tt(l)
          if(ophcd(1:nc).eq.phcd(l)(1:nc)) then
            resid = tt(l)-obst
            ttc = tt(l)
            if(delta.lt.25.0) goto 22
            call phjset(ophcd,ncd)
            resv = resid/sqrt(ovar)
            psi(ncd) = psi(ncd) + resv*resv
            nps(ncd) = nps(ncd) + 1
            aht(ncd) = aht(ncd) + ohit
*            write(6,*) ncd,psi(ncd),nps(ncd),aht(ncd)
            go to 22
          endif
 21     continue
c
 22   continue
 30   continue
      go to 10
c                                    end delta loop
 13   continue
c>                             residuals    
      do 300 j=1,18
        psi(j) = psi(j)/float(nps(j))
        aht(j) = aht(j)/float(nps(j))
        write(nout,fmt='(a,2x,f10.3,i6,f10.1)') 
     ^      cbr(j),psi(j),nps(j),aht(j)
 300  continue
      ypp  = psi(1)+psi(3)
      yss  = psi(4)+psi(6)
      ycc  = psi(7)+psi(8)+psi(16)+psi(17)
      psk  = psi(9)+psi(10)+psi(11)     
      psp  = psi(1)+psi(3)+psk
      pss  = psi(4)+psi(6)+psi(8)
      psr  = psi(2)+psi(5)+psi(7)
      pscp = psk+psi(12)+psi(13)
      pscs = psi(14)+psi(15)+psi(16)+psi(17)
      psps = psp+pss
      psma = psi(1)+psi(3)+pss+psr
      psco = pscp+pscs+psi(18)
      psto = psma+psco
      yk1  = psi(11)/psi(10)
      yk2  = psi(9)/psi(10)
      write(nout,*)
      write(nout,16) "P+PcP     ", ypp
      write(nout,16) "S+ScS     ", yss
      write(nout,16) "PKP       ", psk
      write(nout,16) "SKS+SKKS  ", pscs
      write(nout,*)
      write(nout,16) "PP+SS+SP  ", psr
      write(nout,16) "PKP+PKKP  ", pscp
      write(nout,16) "SP+ScP+SKP", ycc
      write(nout,*)
      write(nout,16) "P+PcP+PKP ", psp
      write(nout,16) "S+ScS+ScP ", pss
      write(nout,16) "major P,S ", psps
      write(nout,*) 
      write(nout,16) "mantle ph ", psma
      write(nout,16) "core ph.  ", psco
      write(nout,*)
      write(nout,16) "all ISC   ", psto
      write(nout,17) "rel PKP   ", yk1,yk2
      write(nout,*)
c
      yp1 = psi(1)*wt(1)+psi(3)*wt(3)
     ^     +psi(9)*wt(9)+psi(10)*wt(10)+psi(11)*wt(11)
      ys1 = psi(4)*wt(4)+psi(6)*wt(6)+psi(14)*wt(14)
      ym1 = psi(1)*wt(1)+psi(3)*wt(3)
     ^     +psi(4)*wt(4)+psi(6)*wt(6)+psi(8)*wt(8)
      ym2 = psi(1)*wt(1)+psi(3)*wt(3)
     ^     +psi(4)*wt(4)+psi(6)*wt(6)+psi(8)*wt(8)
     ^     +psi(2)*wt(2)+psi(5)*wt(5)+psi(7)*wt(7)
      ycp = psi(9)*wt(9)+psi(10)*wt(10)+psi(11)*wt(11)
     ^     +psi(12)*wt(12)+psi(13)*wt(13)+psi(18)*wt(18)
      ycs = psi(14)*wt(14)+psi(15)*wt(15)
     ^     +psi(16)*wt(16)+psi(17)*psi(17)
      ya1 = psi(1)*wt(1)+psi(3)*wt(3)
     ^     +psi(9)*wt(9)+psi(10)*wt(10)+psi(11)*wt(11)
     ^     +psi(4)*wt(4)+psi(6)*wt(6)+psi(14)*wt(14)
     ^     +psi(8)*wt(8)
      ya2 = psi(1)*wt(1)+psi(3)*wt(3)
     ^     +psi(9)*wt(9)+psi(10)*wt(10)+psi(11)*wt(11)
     ^     +psi(4)*wt(4)+psi(6)*wt(6)+psi(14)*wt(14)
     ^     +psi(8)*wt(8)
     ^     +psi(12)*wt(12)+psi(13)*wt(13)+psi(15)*wt(15)
     ^     +psi(16)*wt(16)+psi(17)*psi(17)
       yal = 0.0
       DO j=1,18
         yal = yal + psi(j)*wt(j)
       ENDDO
      write(nout,*)
      write(nout,*) "weighted measures"
      write(nout,16) "P1  ", yp1
      write(nout,16) "S1  ", ys1
      write(nout,16) "M1  ", ym1
      write(nout,16) "M2  ", ym2
      write(nout,16) "CP  ", ycp
      write(nout,16) "CS  ", ycs
      write(nout,*)
      write(nout,16) "A1  ", ya1
      write(nout,16) "A2  ", ya2
      write(nout,16) "AL  ", yal
c
      call retrns(in)
      call retrns(10)
c      call exit(0) 
      return
16    format(a10,f10.3)
17    format(a10,10x,2f5.2)
      end
*
      subroutine phjset(phcd,ncd)
      character*8 phcd
c
      if (phcd.eq.'Pn' .or. phcd.eq.'P') then
        ncd = 1
      elseif (phcd.eq.'PP') then 
        ncd = 2
      elseif (phcd.eq.'PcP') then 
        ncd = 3
      elseif (phcd.eq.'Sn' .or. phcd.eq.'S') then
        ncd = 4
      elseif (phcd.eq.'SS') then 
        ncd = 5
      elseif (phcd.eq.'ScS') then 
       ncd = 6 
      elseif (phcd.eq.'SPn'.or. phcd.eq.'SP') then 
       ncd = 7
      elseif (phcd.eq.'ScP') then 
       ncd = 8
      elseif (phcd.eq.'PKiKP' .or. phcd.eq.'PKPdf') then 
       ncd = 9
      elseif (phcd.eq.'PKPbc') then 
       ncd = 10
      elseif (phcd.eq.'PKPab') then 
       ncd = 11
      elseif (phcd.eq.'PKKPbc') then 
       ncd = 12
      elseif (phcd.eq.'PKKPab') then 
       ncd = 13
      elseif (phcd.eq.'SKSac') then
       ncd = 14
      elseif (phcd.eq.'SKKSac') then
       ncd = 15
      elseif (phcd.eq.'SKiKP' .or. phcd.eq.'SKPdf') then
       ncd = 16
      elseif (phcd.eq.'SKPbc') then
       ncd = 17
      elseif (phcd.eq."P'P'bc" .or. phcd.eq."P'P'df") then
       ncd = 18
      else
       ncd = 20
      endif
      return
      end
*
      integer function icharcn (c)
c
c     AUTHOR:  Brian L.N. Kennett RSES ANU
c     DATE:    December 1985
c     PURPOSE:
c     * This routine returns the number of non-blank
c       characters in a string  ( up to first blank)
      character*8 c
c     integer charcn
      do 10 k=1,8
        if(c(k:k).eq.' ') then
*          icharcn = k-1
            icharcn = k
            return
        endif
 10   continue
      icharcn = 8 
      return
      end
      integer function lenstr (cs)
c
c     AUTHOR:  Brian L.N. Kennett RSES ANU
c     DATE:    December 1985
c     PURPOSE:
c     * This routine returns the number of non-blank
c       characters in a string  ( up to first blank)
      character*20 cs
c     integer lenstr
      do 10 k=1,20
        if(cs(k:k).eq.' ') then
            lenstr = k-1
*            lenstr = k
            return
        endif
 10   continue
      lenstr = 20 
      return
      end
