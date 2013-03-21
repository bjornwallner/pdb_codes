MAIN=molecule.o ProQ_contacts.o src/nrutil.c src/jacobi.c src/eigsrt.c src/nets.c

CCFLAG=-static -static-libgcc
LFLAG=-lm -g -O3 -funroll-loops -Isrc/ 
MAIN: $(MAIN)
	$(CC) $(LFLAG) -o bin/ProQ_contacts $(MAIN)


ProQ_contacts.o: src/ProQ_contacts.c 
	$(CC) $(LFLAG) $(CCFLAG) -c src/ProQ_contacts.c

molecule.o: src/molecule.c src/molecule.h
	$(CC) $(LFLAG) $(CCFLAG) -c src/molecule.c

.c.o:	
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c
