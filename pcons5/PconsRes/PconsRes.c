#include <stdio.h>
#include <stdlib.h>
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h"
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.h"
#define MAXFILE         1000            /* Maximum allowable models */
//#define MAXRES          10000
//#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_molecule.h"

main(int argc,char *argv[])             /* Main routine */
{
  FILE *fp;
  double sim[MAXFILE][MAXRES]={{0}};
  int i=0,j=0,k=0;
  int files=0;
  char filenames[MAXFILE][1000];
  char filename[1000];
  char outfile[1000];
  int number_of_comparison[MAXRES]={0};
  double *Sstr;
  int len;

   if(argc==4)
      {
	strcpy(filename,argv[1]);
	len=atoi(argv[2]);
	strcpy(outfile,argv[3]);
	//	printf("%d %s\n",len,outfile);
	fp=fopen(filename,"r");
	if(fp!=NULL)
	  {
	    while(fscanf(fp,"%s",filenames[i])!=EOF)
	      {
			//		strcpy(filenames[i],buff);
		//tmp=ftell(fp);
		//printf("%s %d %d\n",filenames[i],i,tmp);
		i++;
	      }
	    files=i;
	  }
	else
	  {
	    printf("Cannot open %s",filename);
	    exit(1);
	  }
	fclose(fp);
      }
    else
      {
	printf("\nCalculate average residue similarity\n\nUsage:    PconsRes infile (file with full path to all models) seq-len outfile\n");
	exit(1);
      }
   //exit(1);
   //for(i=0;i<=files;i++)
   //   {
   //	strcpy(m[0].filename,filenames[i]);
   //	m[0].filename[strlen(filenames[i])]='\0';
   //  }
    for(i=0;i<files;i++)
      {
	//printf("%s ",filenames[i]);
	for(j=i+1;j<files;j++)
	   {
	     //printf("%s %s\n",filenames[i],filenames[j]);
	     //if(strcmp(method[i],method[j])!=0)
	       {
		 Sstr=LGscore_res(filenames[i],filenames[j]);
		 k=0;
		 while(Sstr[k] != -1)
		   {
		     sim[i][k]+=Sstr[k];
		     sim[j][k]+=Sstr[k];
		     //printf("%d %f\n",k,Sstr[k]);
		     k++;
		   }
		 number_of_comparison[i]++;
		 number_of_comparison[j]++;
		 free(Sstr);
		 // exit(1);
	       }
	   }
      }
    fp=fopen(strcat(outfile,""),"w");
    if(fp!=NULL)
      {
	for(j=0;j<len;j++)
	  {
	    //printf("%4d ",j+1);
	    for(i=0;i<files;i++)
	      {
		//printf("%d %d %d\n",i,number_of_comparison[i],files);
		sim[i][j]/=number_of_comparison[i];
		//		fprintf(fp,"%d %8.5lf %d %8.5lf",j,sim[i][j],number_of_comparison[i],sim[i][j]/number_of_comparison[i]);
		fprintf(fp,"%8.5lf ",sim[i][j]);
	      }
	    fprintf(fp,"\n");
	  }
	fclose(fp);
      }
    else
      {
	fprintf(stderr,"Cannot open %s\n",strcat(outfile,""));
      }
}
