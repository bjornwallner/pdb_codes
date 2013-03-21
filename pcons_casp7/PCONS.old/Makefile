MAIN=lgscore.o molecule.o pcons.o src/nrutil.c src/jacobi.c src/eigsrt.c src/nets.c


LFLAG=-lm -O3 -funroll-loops -Isrc/ 
MAIN: $(MAIN)
	$(CC) $(LFLAG) -o bin/pcons $(MAIN)


pcons.o: src/pcons.c
	$(CC) $(LFLAG) $(CCFLAG) -c src/pcons.c

molecule.o: src/molecule.c src/molecule.h
	$(CC) $(LFLAG) $(CCFLAG) -c src/molecule.c

lgscore.o: src/lgscore.c src/lgscore.h
	$(CC) $(LFLAG) $(CCFLAG) -c src/lgscore.c


.c.o:	
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c
