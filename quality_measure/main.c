#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <sys/times.h>
#include <time.h>
#include "src/molecule.h"
#include "lgscore.h"
#include <unistd.h>

//#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_molecule.h"


main(int argc,char *argv[])             /* Main routine */
{
  FILE     *fp;
  char     file1[1000]="undef";
  char     file2[1000]="undef";
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
  

  
 
  
    
  //char     models[1000][1000];
  //char     pdbs[1000][1000];

  //molecule models[1000];
  //molecule pdbs[1000];

  lgscore LG[1];

  /* Parse command line for PDB files and options */
  i=1;
  /*    temp=1; */
  if(argc > 2)
    {
      while (i<argc)
	{
	  if(strcmp(argv[i],"-L")==0)
	    {
	      i++;
	      L=atoi(argv[i]);
	    }
	  else if(strcmp(argv[i],"-s")==0)
	    {
	      i++;
	      step=atoi(argv[i]);
	    }
	  else if(strcmp(argv[i],"-minsim")==0)
	    {
	      i++;
	      minsim=atof(argv[i]);
	    }
	  else if(strcmp(argv[i],"-factor")==0)
	    {
	      i++;
	      factor=atof(argv[i]);
	    }
	  else if(strcmp(argv[i],"-d0")==0)
	    {
	      i++;
	      d0=atof(argv[i]);
	    }
	  else if(strcmp(argv[i],"-list")==0)
	    {
	      i++;
	      strcpy(filelist,argv[i]);
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
   
      st_time=times(&st_cpu);
      if(strcmp(filelist,"undef")==0)
	{
	  
	      LGscore_res(file1,file2,LG,d0,minsim,L,factor,step);

      //for(i=0;i<LG[0].residues;i++)
      //	{
      //	  printf("%d %c %lf %lf\n",LG[0].resnum[i],LG[0].residue[i],LG[0].S[i],LG[0].TM[i]);
      //	}
	      //printf("SCORE: %lf %lf %lf \n", LG[0].Ssum,LG[0].TMsum,LG[0].LGscore);
	      for(i=0;i<LG[0].residues;i++)
		{
		  printf("%d %c %lf %lf\n",LG[0].resnum[i],LG[0].residue[i],LG[0].S[i],LG[0].TM[i]);
		}
	      printf("SCORE[Ssum TMsum -log(LGscore)]: %lf %lf %lf \n", LG[0].Ssum,LG[0].TMsum,LG[0].LGscore);
	}
      else
	{
	  fp=fopen(filelist,"r");
	  files=0;
	  if (fp!=NULL)	/* If yes, read in coordinates */
	    {
	      while(fscanf(fp,"%s %s",file1,file2)!=EOF)
		{
		  

		  //strcpy(models[files].filename,file1);
		  //strcpy(pdbs[files].filename,file2);
		  
		   LGscore_res(file1,file2,LG,d0,minsim,L,factor,step);
      //for(i=0;i<LG[0].residues;i++)
      //	{
      //	  printf("%d %c %lf %lf\n",LG[0].resnum[i],LG[0].residue[i],LG[0].S[i],LG[0].TM[i]);
      //	}
		   //printf("%s %s SCORE: %lf %lf %lf \n",file1,file2, LG[0].Ssum,LG[0].TMsum,LG[0].LGscore);
		   printf("SCORE: %s %s %lf %lf %lf %ld \n",file1,file2,LG[0].Ssum,LG[0].TMsum,LG[0].LGscore); //,clock());
		  
		  files++;
		}

	    }
	  else
	    {
	      fprintf(stderr,"Cannot open %s\n",filelist);
	    }


	}


      
     //en_time=times(&en_cpu);
     ////cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
     ////printf("Real Time: %ld, User Time: %ld, System Time: %ld\n",
     ////	     ((en_time - st_time)),
     ////	     ((en_cpu.tms_utime - st_cpu.tms_utime)),
     ////	     ((en_cpu.tms_stime - st_cpu.tms_stime)));
     //
     ////printf("%ld\n",sysconf(_SC_CLK_TCK));
     //
     ////LGscore_res(file1,file2,LG,d0,minsim,L,factor);
     //
     //for(i=0;i<LG[0].residues;i++)
     //	{
     //	  printf("%d %c %lf %lf\n",LG[0].resnum[i],LG[0].residue[i],LG[0].S[i],LG[0].TM[i]);
     //	}
     //printf("SCORE: %lf %lf %lf \n", LG[0].Ssum,LG[0].TMsum,LG[0].LGscore);


       
      
      //printf("L: %d\n",L);
      //printf("d0: %lf\n",d0);
      //
      //printf("minsim: %lf\n",minsim);
      //printf("factor: %lf\n",factor);
      

      //printf("residues: %d %d %lf\n",LG[0].residues,LG[0].length,LG[0].maxrms);
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
#ifdef Sscore
      fprintf(stderr,"Sscore [file1] [file2] -d0 [default=sqrt(5)] -L [4] -minsim [141.0] -factor [0.5] -s [step=2]\n");
#else
      fprintf(stderr,"lgscore [file1] [file2] -d0 [default=sqrt(5)] -L [4] -minsim [141.0] -factor [0.5] -s [step=2]\n");
#endif     
      fprintf(stderr,"\tDefault values are optimized for speed and performance\n");
      fprintf(stderr,"\tOutput\n");
      fprintf(stderr,"\t\tThe following values are reported for each residue:\n");
      fprintf(stderr,"\t\t[resnum] [residue] [S-score] [TM-score]\n");
      
      fprintf(stderr,"\t\tFinally three scores are reported corresponding to:\n");
      fprintf(stderr,"\t\tSsum, TMsum, -log(LGscore) respectively\n");
      
    }
}


