CONTACT_MAP=molecule.o contact_map.o  src/nrutil.c src/jacobi.c src/eigsrt.c src/nets.c

#CCFLAG=-static -static-libgcc 
#CCFLAG=
LFLAG=-O3 -funroll-loops -Isrc/


contact_map: $(CONTACT_MAP)
	$(CC) $(LFLAG) -o contact_map_chain $(CONTACT_MAP) -lm 

contact_map.o: src/contact_map.c
	$(CC) $(LFLAG) $(CCFLAG) -c src/contact_map.c -lm 
molecule.o: src/molecule.c src/molecule.h 
	$(CC) $(LFLAG) $(CCFLAG) -c src/molecule.c -lm 

.c.o:	
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c -lm 


