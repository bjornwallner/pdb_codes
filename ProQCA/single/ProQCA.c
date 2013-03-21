#include <stdio.h>
#include <stdlib.h>
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h"


main(int argc,char *argv[])		/* Main routine */
{
  int           i=0,j=0;
  FILE          *fp_pdb;
  FILE          *fp_psi;
  char	        psipredfiles[1000];
  char          psipred[2000];
  char          pdbfiles[1000];
  char          pdbfile[1000];
  double        *features;
  double        *quality;
  
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

  fp_pdb=fopen(pdbfiles,"r");
  fp_psi=fopen(psipredfiles,"r");
  if (fp_pdb!=NULL && fp_psi!=NULL) 
    {
      while(fscanf(fp_pdb,"%s",pdbfile)!=EOF && fscanf(fp_psi,"%s",psipred))
	{
	  //printf("%s",psipred);
	  quality=(double*)ProQCA(pdbfile,psipred);
	    if(quality[0]!=-1)
	      {
		printf("%-30s %6.4lf %6.4lf\n",pdbfile,quality[0],quality[1]);
	      }
	    else
	      {
		printf("%-30s Error\n",pdbfile);
	      }
	    free(quality);
	}
    }
  fclose(fp_pdb);
  fclose(fp_psi);

  //exit(1);
  //
  //
  //strcpy(lg1,getenv("PROQDIR"));
  //strcpy(lg2,lg1);
  //strcpy(lg3,lg1);
  //strcpy(lg4,lg1);
  //strcpy(lg5,lg1);
  //strcpy(mx1,lg1);
  //strcpy(mx2,lg1);
  //strcpy(mx3,lg1);
  //strcpy(mx4,lg1);
  //strcpy(mx5,lg1);
  //
  //read_net((char*)strcat(lg1,"lg1"),&net_lg[0]);
  //read_net((char*)strcat(lg2,"lg2"),&net_lg[1]);
  //read_net((char*)strcat(lg3,"lg3"),&net_lg[2]);
  //read_net((char*)strcat(lg4,"lg4"),&net_lg[3]);
  //read_net((char*)strcat(lg5,"lg5"),&net_lg[4]);
  //
  //read_net((char*)strcat(mx1,"mx1"),&net_mx[0]);
  //read_net((char*)strcat(mx2,"mx2"),&net_mx[1]);
  //read_net((char*)strcat(mx3,"mx3"),&net_mx[2]);
  //read_net((char*)strcat(mx4,"mx4"),&net_mx[3]);
  //read_net((char*)strcat(mx5,"mx5"),&net_mx[4]);
  //
  //
  //fp_pdb=fopen(pdbfiles,"r");
  //fp_psi=fopen(psipredfiles,"r");
  //if (fp_pdb!=NULL && fp_psi!=NULL) 
  //  {
  //	while(fscanf(fp_pdb,"%s",pdbfile)!=EOF && fscanf(fp_psi,"%s",psipred))
  //	  {
  //	    //printf("%s\n",pdbfile);
  //	    //printf("%s\n",psipred);
  //	    //  printf("%d %d\n",ftell(fp_pdb),ftell(fp_psi));
  //	    features=calculate_parameters(pdbfile,psipred);
  //	    // for(i=0;i<35;i++)
  //	    //  printf("%lf ",features[i]);
  //	    pred_lg=0;
  //	    pred_mx=0;
  //	    for(i=0;i<5;i++)
  //	      pred_lg+=netfwd(features,&net_lg[i]);
  //	    for(i=0;i<5;i++)
  //	      pred_mx+=netfwd(features,&net_mx[i]);
  //	    printf("%-30s %6.4lf %6.4lf\n",pdbfile,pred_lg/5,pred_mx/5);
  //	    free(features);
  //	  }
  //  }
  //fclose(fp_pdb);
  //fclose(fp_psi);
}



