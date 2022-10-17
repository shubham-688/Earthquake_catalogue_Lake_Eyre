C-----------------------------------------------------------
c       mesrank
c
c                    routine for ranking misfit outcomes
c                    across a suite of randomly generated models
c                    the data is read from summary Q*.out files
c                    which contain an entry for each model in the form
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
c                    the measures are then ranked across the assemblage of 
c                    models and then best 50 printed out in order
c                    for the particular measure
c
c--------------------------------------------------------------------------
c                    B.L.N. Kennett, RSES ANU
c                    June 1994
c--------------------------------------------------------------------------	
c
	character*10 modnam(990),comment1
        character*20 arg1,filnam1,arg2,filnam2
        character*2 cph,cwt 
        real psi(18),wt(18)
        real eta(990)
        integer ind(990),irk(990)
        data wt/5.0,1.5,2.0,3.0,1.0,1.0,1.0,2.0,4.0,4.0,4.0,
     &          1.5,1.5,3.0,1.5,1.0,1.0,1.0/
c
        write(6,*) "Enter choice of phase option"
        read (5,*) cph
        write(6,*) "Applying weighting? (YE/NO)"
        read (5,*) cwt
c
        callgetarg(1,arg1)
        read(arg1,*) filnam1
        open(2,file=filnam1,status='old')
        callgetarg(2,arg2)
        read(arg2,*) filnam2
        open(8,file=filnam2,status='unknown')
        do 50 k=1,990
          eta(k) = 100.0
          ind(k) = 0
          irk(k) = 990
 50     continue
c
        do 100 k=1,990
          read(2,*,end=101) modnam(k)
	  do 10 j=1,18
            read(2,*) comment1,psi(j),ij,rj
 10	  continue
          do 11 j=19,39
            read(2,*)
 11       continue
          if    (cph.eq."P1") then
            if(cwt.eq."YE") then
              eta(k) = psi(1)*wt(1)+psi(3)*wt(3)
     ^                +psi(9)*wt(9)+psi(10)*wt(10)+psi(11)*wt(11)
            else
              eta(k) = psi(1)+psi(3)
     ^                +psi(9)+psi(10)+psi(11)
            endif
          elseif(cph.eq."S1") then
            if(cwt.eq."YE") then
              eta(k) = psi(4)*wt(4)+psi(6)*wt(6)+psi(14)*wt(14)
            else
              eta(k) = psi(4)+psi(6)+psi(14)
            endif
          elseif(cph.eq."M1") then
            if(cwt.eq."YE") then
              eta(k) = psi(1)*wt(1)+psi(3)*wt(3)
     ^                +psi(4)*wt(4)+psi(6)*wt(6)+psi(8)*wt(8)
            else
              eta(k) = psi(1)+psi(3)
     ^                +psi(4)+psi(6)+psi(8)
            endif
          elseif(cph.eq."M2") then
            if(cwt.eq."YE") then
              eta(k) = psi(1)*wt(1)+psi(3)*wt(3)
     ^                +psi(4)*wt(4)+psi(6)*wt(6)+psi(8)*wt(8)
     ^                +psi(2)*wt(2)+psi(5)*wt(5)+psi(7)*wt(7)
            else
            endif
          elseif(cph.eq."CP") then
            if(cwt.eq."YE") then
              eta(k) = psi(9)*wt(9)+psi(10)*wt(10)+psi(11)*wt(11)
     ^                +psi(12)*wt(12)+psi(13)*wt(13)+psi(18)*wt(18)
            else
              eta(k) = psi(9)+psi(10)+psi(11)
     ^                +psi(12)+psi(13)+psi(18)
           endif
          elseif(cph.eq."CS") then
            if(cwt.eq."YE") then
              eta(k) = psi(14)*wt(14)+psi(15)*wt(15)
     ^                +psi(16)*wt(16)+psi(17)*psi(17)
            else
              eta(k) = psi(14)+psi(15)
     ^                +psi(16)+psi(17)
            endif
          elseif(cph.eq."A1") then
            if(cwt.eq."YE") then
              eta(k) = psi(1)*wt(1)+psi(3)*wt(3)
     ^                +psi(9)*wt(9)+psi(10)*wt(10)+psi(11)*wt(11)
     ^                +psi(4)*wt(4)+psi(6)*wt(6)+psi(14)*wt(14)
     ^                +psi(8)*wt(8)
            else
              eta(k) = psi(1)+psi(3)
     ^                +psi(9)+psi(10)+psi(11)
     ^                +psi(4)+psi(6)+psi(14)
     ^                +psi(8)
            endif
          elseif(cph.eq."A2") then
            if(cwt.eq."YE") then
              eta(k) = psi(1)*wt(1)+psi(3)*wt(3)
     ^                +psi(9)*wt(9)+psi(10)*wt(10)+psi(11)*wt(11)
     ^                +psi(4)*wt(4)+psi(6)*wt(6)+psi(14)*wt(14)
     ^                +psi(8)*wt(8)
     ^                +psi(12)*wt(12)+psi(13)*wt(13)+psi(15)*wt(15)
     ^                +psi(16)*wt(16)+psi(17)*psi(17)
            else
              eta(k) = psi(1)+psi(3)
     ^                +psi(9)+psi(10)+psi(11)
     ^                +psi(4)+psi(6)+psi(14)
     ^                +psi(8)
     ^                +psi(12)+psi(13)+psi(15)
     ^                +psi(16)+psi(17)*psi(17)
            endif
          elseif(cph.eq."AL") then
            if(cwt.eq."YE") then
              eta(k) = 0.0
              do 21 j=1,18
                eta(k) = eta(k) + psi(j)*wt(j)
 21           continue
            else
              eta(k) = 0.0
              do 22 j=1,18
                eta(k) = eta(k) + psi(j)
 22           continue
            endif
          endif
 100	continue
c
 101    close(2)
        kr = k
        call indexx(990,eta,ind)
        call rank(990,ind,irk)
        write(8,*) cph
        do 200 l=1,50
          il = ind(l)
          write(8,*) modnam(il),"   ",eta(il)
200     continue
c
        close(8)
        stop
        end
c++
      subroutine indexx(n,arrin,indx)
      dimension arrin(n),indx(n)
      do 11 j=1,n
        indx(j)=j
11    continue
      l=n/2+1
      ir=n
10    continue
        if(l.gt.1)then
          l=l-1
          indxt=indx(l)
          q=arrin(indxt)
        else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir.eq.1)then
            indx(1)=indxt
            return
          endif
        endif
        i=l
        j=l+l
20      if(j.le.ir)then
          if(j.lt.ir)then
            if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1
          endif
          if(q.lt.arrin(indx(j)))then
            indx(i)=indx(j)
            i=j
            j=j+j
          else
            j=ir+1
          endif
        go to 20
        endif
        indx(i)=indxt
      go to 10
      end
      subroutine rank(n,indx,irank)
      dimension indx(n),irank(n)
      do 11 j=1,n
        irank(indx(j))=n-j+1
11    continue
      return
      end
