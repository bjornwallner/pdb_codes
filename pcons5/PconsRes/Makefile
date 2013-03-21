PconsRes=molecule.o lgscore.o PconsRes.o /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/nrutil.c /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/jacobi.c  /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/eigsrt.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/nets/nets.c
#PconsRes=lgscore.o PconsRes.o 

CCFLAG=-I /afs/pdc.kth.se/home/b/bjornw/bjorn/modules/numrecipies/ 
LFLAG=-lm -g -O3
#LFLAG=-lm -g
PconsRes: $(PconsRes)
	$(CC) $(LFLAG) -o PconsRes $(PconsRes)

molecule.o: /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h
	$(CC) $(LFLAG) $(CCFLAG) -c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.c

lgscore.o: /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.h
	$(CC) $(LFLAG) $(CCFLAG) -c /afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.c


.c.o:	
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c
