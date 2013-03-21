MAIN=main.o lgscore.o molecule.o /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/nrutil.c /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/jacobi.c  /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/eigsrt.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/nets/nets.c

#-fbranch-probabilities
LFLAG=-lm -O3 -funroll-loops #-fbranch-probabilities #-fprofile-arcs #-finline-functions #-fomit-frame-pointer -ffast-math -funroll-loops # -fprofile-arcs
MAIN: $(MAIN)
	$(CC) $(LFLAG) -o superimpose $(MAIN)

molecule.o: /afs/pdc.kth.se/home/b/bjornw/source/c/pdb/molecule.c /afs/pdc.kth.se/home/b/bjornw/source/c/pdb/molecule.h
	$(CC) $(LFLAG) $(CCFLAG) -c /afs/pdc.kth.se/home/b/bjornw/source/c/pdb/molecule.c

lgscore.o: /afs/pdc.kth.se/home/b/bjornw/source/c/pdb/quality_measure/lgscore.c /afs/pdc.kth.se/home/b/bjornw/source/c/pdb/quality_measure/lgscore.h
	$(CC) $(LFLAG) $(CCFLAG) -c /afs/pdc.kth.se/home/b/bjornw/source/c/pdb/quality_measure/lgscore.c

.c.o:	
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c
