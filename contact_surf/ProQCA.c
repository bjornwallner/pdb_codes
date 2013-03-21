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

//double		Levitt_Gerstein();
//double		superimpose_molecules();
//void            strncpy_NULL(char *dest, char *src, size_t n);
//int             read_molecules();
//int             get_type(char *name, char *res);
//double          distance(int atomno1, int atomno2);
//void            print_type(int type_no, FILE *fp);
//int             get_res6(char *res);
//int             get_res(char *res);
//void            print_res(int res,FILE *fp);
//double          crd(int atomno1, int atomno2);   /*closest residue distance */


main(int argc,char *argv[])		/* Main routine */
     //int argc;
     //char *argv[];
{
  int           i=0,j=0;
  network       net1;net2;net3;net4;net5;
  FILE          *fp;
  FILE          *fp2;
  char	        psipredfiles[1000];
  char          psipred[2000];
  char          pdbfiles[1000];
  char          pdbfile[1000];
  double        *features[1000];
  if(argc==3)
    {
      strcpy(pdbfiles,argv[1]);
      strcpy(psipredfiles,argv[2]);
      printf("%s\n%s\n",pdbfiles,psipredfiles);
      //exit(1);
    }
  else
    {
      printf("Usage: read_pdb [pdb_files] [psipredfiles]\n");
      exit(1);
    }
  


  fp=fopen(pdbfiles,"r");
  fp2=fopen(psipredfiles,"r");
  if (fp!=NULL && fp2!=NULL)
    {
      while(fscanf(fp,"%s",pdbfile)!=EOF && fscanf(fp2,"%s",psipred))
	{
	  printf("%s\n",pdbfile);
	  printf("%s\n",psipred);
	  features[i]=calculate_parameters(pdbfile,psipred);
	  i++;
	}
    }
  fclose(fp);
  fclose(fp2);

  printf("\n");
  for(j=0;j<35;j++)
    {
      
      printf("%lf ",features[0][j]);
    }
  printf("\n");

  
  for(i=i-1;i>=0;i--)
    free(features[i]);
}



