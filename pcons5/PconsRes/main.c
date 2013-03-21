#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h"
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.h"

//#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_molecule.h"

main(int argc,char *argv[])             /* Main routine */
{
  char file1[1000];
  char file2[1000];
  double* Sstr;
  int     i;
  
  if(argc==3)
    {
      strcpy(file1,argv[1]);
      strcpy(file2,argv[2]);
      LGscore_res(file1,file2);
      //for(i=0;i < 1000;i++)
      //	{
	  //printf("%lf %d\n",Sstr[i],sizeof(Sstr[i]));
      //	}
      //free(Sstr);
      //printf("%d\n",sizeof(double));
    }
}


