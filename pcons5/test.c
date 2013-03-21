#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h"
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.h"

//#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_molecule.h"

main(int argc,char *argv[])             /* Main routine */
{
  molecule m[1];
  char file[200]="/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/pcons5/testdir/1fl9A.2reb.dali.1.bb";
  
  strcpy(m[0].filename,file);
  read_molecules(&m[0],'c');

}
