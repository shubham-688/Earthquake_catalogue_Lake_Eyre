#
IASQ=emtvelin.o
TAUL=libtau.o libsun.o
POST= npost.o
#
#
#  path of executables
BIN = ../bin/
#
RST = rst.o remodl.o setbrn.o tqresidm.o
rst: $(RST) $(TAUL) $(IASQ)
	gfortran -o  $(BIN)rst $(RST) $(IASQ) $(TAUL)
	
#  resampling of model 
REM = remodl.o
remodlv: $(REM) $(TAUL) $(IASQ)
	gfortran -o $(BIN)remodlv $(REM) $(IASQ) $(TAUL)
#
#  set branches  program
SBR = setbrn.o
setbrn: $(SBR) $(TAUL) 
	gfortran -o $(BIN)setbrn $(SBR) $(TAUL)
#
# travel time branches program
TBR =  ttimes.o       
ttimes: $(TBR) $(TAUL)
	gfortran -o $(BIN)ttimes $(TBR) $(TAUL)
#
