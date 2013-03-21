MAIN=molecule.o ProQres.o /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/nrutil.c /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/jacobi.c  /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/eigsrt.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/nets/nets.c

CCFLAG=-static -static-libgcc
LFLAG=-lm -g
MAIN: $(MAIN)
	$(CC) $(LFLAG) -o ProQres $(MAIN)

molecule.o: molecule.c molecule.h
	$(CC) $(LFLAG) $(CCFLAG) -c molecule.c

.c.o:	
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c
