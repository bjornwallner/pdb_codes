CONTACT_SURF=molecule.o contact_surf.o /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/nrutil.c /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/jacobi.c  /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/eigsrt.c

CCFLAG=-I /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/
LFLAG=-lm
contact_surf: $(CONTACT_SURF)
	$(CC) $(LFLAG) -o contact_surf $(CONTACT_SURF)
.c.o:	
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c
