PCONS5=molecule.o lgscore.o pcons5.o /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/nrutil.c /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/jacobi.c  /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/eigsrt.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/nets/nets.c

CCFLAG=-I /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/ -static
LFLAG=-lm -g -O3 
PCONS5: $(PCONS5)
	$(CC) $(LFLAG) -o pcons5 $(PCONS5)

molecule.o: /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h
	$(CC) $(LFLAG) $(CCFLAG) -c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.c

lgscore.o: /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.h
	$(CC) $(LFLAG) $(CCFLAG) -c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.c


.c.o:	
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c
