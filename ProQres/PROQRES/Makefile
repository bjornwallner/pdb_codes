MAIN=molecule.o ProQres.o src/nrutil.c src/jacobi.c src/eigsrt.c src/nets.c

#CCFLAG=-static -static-libgcc -mdynamic-no-pic
CCFLAG=
LFLAG=-lm -g -O3 -funroll-loops -Isrc/ 
MAIN: $(MAIN)
#	$(CC) $(LFLAG) -o bin/ProQres $(MAIN)
	$(CC) $(LFLAG) -o bin/ProQres.Darwin $(MAIN)


ProQres.o: src/ProQres.c 
	$(CC) $(LFLAG) $(CCFLAG) -c src/ProQres.c

molecule.o: src/molecule.c src/molecule.h
	$(CC) $(LFLAG) $(CCFLAG) -c src/molecule.c

.c.o:	
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c
