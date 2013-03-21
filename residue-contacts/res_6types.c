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
  int           res_contacts[6][6]={{}};
  int           restype[6]={};
  int           tot_res_contacts=0;
  int           type_i,type_j;
  int           current_res_i=0;
  int           current_res_j=0;
  double        cutoff;
  FILE          *fp;
  double temp;
  int tmp=0;
  float         sum=0;
  int           sum2=0;
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
	  // printf("%d %d %d %d %d %d %f \n",m[0].atm[i].rescount,current_res_i,m[0].atm[j].rescount,current_res_j,m[0].atm[i].resnum,m[0].atm[j].resnum,crd(m,i,j));
	  //printf("%d %d\n",m[0].atm[i].rescount,current_res_i);
	  if(m[0].atm[i].rescount!=current_res_i)
	    {
	      for(j=0,current_res_j=0;j<m[0].atoms;j++)
		{
		  // printf("%d %d %d %d %d %d %f \n",m[0].atm[i].rescount,current_res_i,m[0].atm[j].rescount,current_res_j,m[0].atm[i].resnum,m[0].atm[j].resnum,crd(m,i,j));
		  if(m[0].atm[j].rescount!=current_res_j) //current_res_j)
		    {
		      //printf("%d %d %d %d\n",m[0].atm[i].rescount,m[0].atm[j].rescount,m[0].atm[i].resnum,m[0].atm[j].resnum);
		      //printf("%d %d %d %d %d %d %f \n",m[0].atm[i].rescount,current_res_i,m[0].atm[j].rescount,current_res_j,m[0].atm[i].resnum,m[0].atm[j].resnum,crd(m,i,j));
		      //printf("%d %d %d %d %f %f\n",current_res_i+1,current_res_j+1,i,j,crd(i,j),crd(j,i));
		      if(abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>5 &&
			 crd(m,i,j)<cutoff)
			{
			  res_contacts[get_res6(m[0].atm[i].residue)][get_res6(m[0].atm[j].residue)]++;
			  tot_res_contacts++;
			}
		      current_res_j++;
		    }
		}				     
	      current_res_i++;
	    }
	}

      for(i=0;i<6;i++)
	{
	  for(j=0;j<6;j++)
	    {
	      printf("%3d ",res_contacts[i][j]);
	    }
	  printf("\n");
	}


      /*RESIDUE CONTACTS for each residue type */
      printf("RES: ");
      for(i=0;i<6;i++)
	{
	  for(j=i;j<6;j++)
	    {
	      if(i!=j)
		{
		  res_contacts[i][j]=2*res_contacts[i][j];
		}
	      //restype[i]=restype[i]+res_contacts[i][j];
	      temp=0;
	      if(tot_res_contacts != 0)
		{
		  temp=(double)res_contacts[i][j]/tot_res_contacts;
		}
	      printf("%f ",temp);
	      sum+=temp;
	    }
	}
      printf("\n%d\n",sum);
  
    }
  
}



