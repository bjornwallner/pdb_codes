LGSCORE=lgscore.o molecule.o main.o /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/nrutil.c /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/jacobi.c  /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/eigsrt.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/nets/nets.c

#-fbranch-probabilities
SCORE=$(GOAL)

SSCORE =-D$(SCORE)

LFLAG=-lm -O3 -funroll-loops $(SSCORE)#-fbranch-probabilities #-fprofile-arcs #-finline-functions #-fomit-frame-pointer -ffast-math -funroll-loops # -fprofile-arcs

LGSCORE: $(LGSCORE)
	$(CC) $(LFLAG) -o $(SCORE)_devel.$(ARCH) $(LGSCORE)

molecule.o: /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h
	$(CC) $(LFLAG) $(CCFLAG) -c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.c

lgscore.o: /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.h

.c.o:	
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c


#endif