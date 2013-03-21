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
  //double LG;
  double minsim=49.0;
  double d0=sqrt(5);
  int     i;
  int L=4;
  double factor=sqrt(5);
  int files=0;

  lgscore LG[1];

  /* Parse command line for PDB files and options */
  i=1;
  /*    temp=1; */
  if(argc > 2)
    {
      while (i<argc)
	{
	  if (strcmp(argv[i],"-L")==0)
	    {
	      i++;
	      L=atoi(argv[i]);
	    }
	  else if (strcmp(argv[i],"-minsim")==0)
	    {
	      i++;
	      minsim=atof(argv[i]);
	    }
	  else if (strcmp(argv[i],"-factor")==0)
	    {
	      i++;
	      factor=atof(argv[i]);
	    }
	  else if (strcmp(argv[i],"-d0")==0)
	    {
	      i++;
	      d0=atof(argv[i]);
	    }
	  else
	    {
	      if(files==0)
		{
		  strcpy(file1,argv[i]);
		}
	      else
		{
		  strcpy(file2,argv[i]);
		}
	      files++;
	    }
	  
	  i++;
	}
      //     printf("%lf\n",d0);
      // LG=LGscore(file1,file2);
      //S=(double*)malloc(sizeof(double)*(MAXRES));
      //TM=(double*)malloc(sizeof(double)*(MAXRES));
   


      LGscore_res(file1,file2,LG,d0,minsim,L,factor);

      for(i=0;i<LG[0].residues;i++)
	{
	  printf("%d %c %lf %lf\n",LG[0].resnum[i],LG[0].residue[i],LG[0].S[i],LG[0].TM[i]);
	}
      printf("%lf %lf %lf \n",LG[0].Ssum,LG[0].TMsum,LG[0].LGscore);


       
      
      printf("L: %d\n",L);
      printf("d0: %lf\n",d0);
      //
      printf("minsim: %lf\n",minsim);
      printf("factor: %lf\n",factor);
      

      printf("residues: %d %d %lf\n",LG[0].residues,LG[0].length,LG[0].maxrms);
      //for(i=0;i < 1000;i++)
      //	{
      //	  printf("%lf %d\n",Sstr[i],sizeof(Sstr[i]));
      //	}
      //tfree(Sstr);
      //printf("%d\n",sizeof(double));
      //printf("hej\n");
      //free(S);
      exit(0);
    }
  else
    {
      fprintf(stderr,"lgscore [file1] [file2] -d0 [default=sqrt(5)] -L [4] -minsim [25.0] -factor [2.23]\n");
    }
}


