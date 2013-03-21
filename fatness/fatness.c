#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h"
//#ifdef  NOZLIB
//#else
//#include <zlib.h> 
//#endif
//#include <sys/types.h>
//#include <sys/times.h>
//#include <sys/param.h>
//#include <sys/time.h> 

#define PI		3.14159265	/* Useful constant */

main(argc,argv)		/* Main routine */
     int argc;
     char *argv[];
{
  molecule	m[1];		/* Molecule to be read*/
  //dyn_molecule	*m;
  //  int *ignore_res;
  //int           ignore_res[1000]={0};
  // char          ignore_res_file[PATH_MAX]="undef";
  //double        fat;
  /* Parse command line for PDB filename */
  if(argc==2)
    {
      strcpy(m[0].filename,argv[1]);
    }  
  else
    {
      printf("Usage: read_pdb [pdb_file]\n");
      exit(1);
    }
  if(read_molecules(m,'a')==0)
    {
      //fat=fatness(m);
      //      printf("%5.3f",fatness2(m,1)); 
      fatness2(m,1);
    }
        
}

