      real*8 z1(250),ro1(250),al1(250),be1(250)
      real*8 z2(250),ro2(250),al2(250),be2(250)
      character*20 char1,char2
      character*40 cfile1,cfile2,cname
      
	write(6,*) "tvel file 1"
        read(5,*) cfile1
	write(6,*) "tvel file 2"
	read(5,*) cfile2
c
	open(3,file=cfile1)
	open(4,file=cfile2)
c
c
      read(3,fmt='(a)') char1
      read(4,fmt='(a)') char2
      read(3,*) 
      read(4,*) 

      DO i=1,138
        read(3,*) z1(i),al1(i),be1(i),ro1(i)
        read(4,*) z2(i),al2(i),be2(i),ro2(i)
        dzz = z2(i) - z1(i)
        dal = al2(i)-al1(i)
        dbe = be2(i)-be1(i)
        write(6,fmt = "(f12.2,2f12.4,f12.2)") z1(i), dal,dbe, dzz
      ENDDO
      
      close (3)
      close (4)
      stop
      end
