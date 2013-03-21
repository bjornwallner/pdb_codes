ProQ_contacts=molecule.o ProQ_contacts.o /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/nrutil.c /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/jacobi.c  /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/eigsrt.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/nets/nets.c



CCFLAG=-I /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/
LFLAG=-lm
contact_surf: $(ProQ_contacts)
	$(CC) $(LFLAG) -o ProQ_contacts $(ProQ_contacts)

molecule.o: /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h 
	$(CC) $(LFLAG) $(CCFLAG) -c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.c

.c.o:	
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c
