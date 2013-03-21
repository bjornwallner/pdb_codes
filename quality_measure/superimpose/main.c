#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <sys/times.h>
#include <time.h>
#include "/afs/pdc.kth.se/home/b/bjornw/source/c/pdb/molecule.h"
#include "/afs/pdc.kth.se/home/b/bjornw/source/c/pdb/quality_measure/lgscore.h"
#include <unistd.h>


//#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_molecule.h"

main(int argc,char *argv[])             /* Main routine */
{
  FILE     *fp;
  double   minsim=121.0;
  double   d0=sqrt(5);
  int      i;
  int      L=4;
  double   factor=0.5;
  int      files=0;
  char     filelist[1000]="undef";
  int      step=1;
  double cpu_time_used;
  static clock_t st_time;
  static clock_t en_time;
  static struct tms st_cpu;
  static struct tms en_cpu;
  
  char file1[2000]="undef";
  char file2[2000]="undef";
  char file3[2000]="undef";
 
  
    
  //char     models[1000][1000];
  //char     pdbs[1000][1000];

  //molecule models[1000];
  //molecule pdbs[1000];

  lgscore LG[1];

  /* Parse command line for PDB files and options */
  /*    temp=1; */
  if(argc == 4)
    {
      strcpy(file1,argv[1]);
      strcpy(file2,argv[2]);
      strcpy(file3,argv[3]);
      //printf("%s\n%s\n%s\n",file1,file2,file3);
      superimpose(file1,file2,file3,minsim,L,factor);

    }
}


