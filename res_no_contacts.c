#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "molecule.h"

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
  molecule      m[1];
  int           i,j,k;
  double        dist;
  double        mean=0;
  int           atoms=0;
  int           residues=0;
  int           res_contacts[20][20]={{}};
  int           res_contact[10000]={};
  int           restype[20]={};
  int           tot_res_contacts=0;
  int           type_i,type_j;
  int           current_res_i=0;
  int           current_res_j=0;
  double        cutoff;
  FILE          *fp;
  double temp;
  int tmp=0;
  //float        sum=0;
  /* Parse command line for PDB filename */
  if(argc==3)
    {
      strcpy(m[0].filename,argv[1]);
      cutoff=strtod(argv[2],(char**)(argv[2]+strlen(argv[2])));
      cutoff=cutoff*cutoff;
    }
  else
    {
      printf("Usage: read_pdb [pdb_file] cutoff\n");
      exit(1);
    }
  if(read_molecules(m)==0)
    {  
      for(i=0;i<m[0].atoms;i++)  
	{
	  if(m[0].atm[i].rescount!=current_res_i) // && res_contact[current_res_i]!=1)
	    {
	      for(j=0,current_res_j=0;j<m[0].atoms;j++)
		{
		  if(m[0].atm[j].rescount!=current_res_j)// && res_contact[current_res_j]!=1)
		    {
		      // printf("%d %d %d %d %f %f\n",current_res_i+1,current_res_j+1,i,j,crd(m,i,j),crd(m,j,i));
		      if((res_contact[current_res_i]==0 || 
			  res_contact[current_res_j]==0 ) && 
			 abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>5 &&
			 crd(m,i,j)<cutoff)
			{
			  
			  res_contact[current_res_i]=1;
			  res_contact[current_res_j]=1;
			  
			}
		      current_res_j++;
		    }
		}				     
	      current_res_i++;
	    }
	}
      residues=0;
      for(i=0;i<current_res_i;i++)
	{
	  //printf("%d %d\n", res_contact[i],i);
	  if(res_contact[i]==0)
	    {
	      tot_res_contacts++;
	    }
	  residues++;
	}
      //temp=(double)tot_res_contacts/residues;
      printf("%f",(double)tot_res_contacts/residues);
      
      //printf("%f %f\n",crd(0,16),crd(16,0));
      //fp=fopen("test","w");	/* Does file exist? */
      //if (fp!=NULL)	/* If yes, write output */
      //	{
      //	  fprintf(fp,"\tC\tN\tO\tCA\tCH3\tCH/CH2\tC(OO-)\tNH\tNH2\t(C)00-\t=0\tOH\tS\tOXT\n");
      //	  for(i=0;i<14;i++)
      //	    {
      //	      print_type(i,fp);
//	      fprintf(fp,"\t");
//	      for(j=0;j<14;j++)
//		{
//		  fprintf(fp,"%d\t",atom_contacts[i][j]);
      //	}
	//      fprintf(fp,"\t");
	//      print_type(i,fp);
	  //    fprintf(fp,"\n");
	//      
//	    }
//	}
      
	  //fprintf(fp,"%d\n",m[0].residues);
      //}
      // fclose(fp);
      
   
//	/*RESIDUE CONTACTS for each residue type */
//	printf("RES: ");
//	for(i=0;i<20;i++)
//	  {
//	    for(j=i;j<20;j++)
//	      {
//		if(i!=j)
//		  {
//		    res_contacts[i][j]=2*res_contacts[i][j];
//		  }
//		
//		//restype[i]=restype[i]+res_contacts[i][j];
//		temp=0;
//		if(tot_res_contacts != 0)
//		  {
//		    temp=(double)res_contacts[i][j]/tot_res_contacts;
//		  }
//		printf("%f ",temp);
//		//      sum+=temp;
//		//printf("%d ",res_contacts[i][j]);
//	      }
//	  }
//	//  printf("\n%f\n",sum);
//    }
    }
}



