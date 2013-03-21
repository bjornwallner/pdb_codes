MAIN=lgscore.o molecule.o pcons.o /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/nrutil.c /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/jacobi.c  /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/eigsrt.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/nets/nets.c


LFLAG=-lm -O3 -funroll-loops
MAIN: $(MAIN)
	$(CC) $(LFLAG) -o pcons $(MAIN)

molecule.o: /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h
	$(CC) $(LFLAG) $(CCFLAG) -c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.c

lgscore.o: /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.h
	$(CC) $(LFLAG) $(CCFLAG) -c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.c


.c.o:	
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c
