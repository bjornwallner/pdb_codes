LGSCORE=lgscore.o molecule.o main.o src/nrutil.c src/jacobi.c  src/eigsrt.c src/nets.c

#-fbranch-probabilities
GOAL=Sscore
SCORE=$(GOAL)
SSCORE =-D$(SCORE)
ARCH = $(shell uname -m)
LFLAG=-O3 -funroll-loops -Isrc/ $(SSCORE)

#-fbranch-probabilities #-fprofile-arcs #-finline-functions #-fomit-frame-pointer -ffast-math -funroll-loops # -fprofile-arcs

LGSCORE: $(LGSCORE)
	$(CC) $(LFLAG) -o $(SCORE)_devel2.$(ARCH) $(LGSCORE) -lm

molecule.o: src/molecule.c src/molecule.h
	$(CC) $(LFLAG) $(CCFLAG) -c src/molecule.c -lm

lgscore.o: lgscore.c lgscore.h
	$(CC) $(LFLAG) $(CCFLAG) -c lgscore.c -lm


.c.o:   
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c -lm

clean:
	rm -rf *.o 

#molecule.o: /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h
#	$(CC) $(LFLAG) $(CCFLAG) -c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.c
#
#lgscore.o: /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.h
#
#.c.o:	
#	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c
#
#
##endif