#include <stdio.h>
#include <stdlib.h>
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h"
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.h"

//#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_molecule.h"

main(int argc,char *argv[])             /* Main routine */
{
  char file1[100];
  char file2[100];
  double LG;
  
  if(argc==3)
    {
      strcpy(file1,argv[1]);
      strcpy(file2,argv[2]);
      LG=LGscore(file1,file2);
      printf("%lf\n",LG);
    }


}


