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
  int           contacts8[2000]={};
  int           contacts10[2000]={};
  int           contacts12[2000]={};
  int           contacts14[2000]={};
  int           contacts16[2000]={};
  int           restype[6]={};
  int           tot_res_contacts=0;
  double        cutoff;
  FILE          *fp;
  double temp;
  int tmp=0;
  float         sum=0;
  int           sum2=0;
  /* Parse command line for PDB filename */
  if(argc==2)
    {
      strcpy(m[0].filename,argv[1]);
      //cutoff=strtod(argv[2],(char**)(argv[2]+strlen(argv[2])));
      //cutoff=cutoff*cutoff;
    }
  else
    {
      printf("Usage: read_pdb [pdb_file] cutoff\n");
     exit(1);
   }
  if(read_molecules_CA(m)==0)
    {  
      //printf("%d\n",m[0].residues);
      for(i=0;i<m[0].atoms;i++)  
	{
	  for(j=i+1;j<m[0].atoms;j++)
	    {
	      //printf("%d %d %d %d %f \n",m[0].atm[i].rescount,m[0].atm[j].rescount,m[0].atm[i].resnum,m[0].atm[j].resnum,crd(m,i,j));
	      // //printf("%d %d %d %d %f %f\n",current_res_i+1,current_res_j+1,i,j,crd(m,i,j),crd(m,j,i));
	      //   if(abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>5 &&
	      // distance(m,i,j)<cutoff)
	      //if(abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>5)
	      //{
	      dist=distance(m,i,j);
	      if(dist<256)
		{
		  contacts16[m[0].atm[i].rescount-1]++;
		  contacts16[m[0].atm[j].rescount-1]++;
		  if(dist<196)
		    {
		      contacts14[m[0].atm[i].rescount-1]++;
		      contacts14[m[0].atm[j].rescount-1]++;
		      if(dist<144)
			{
			  contacts12[m[0].atm[i].rescount-1]++;
			  contacts12[m[0].atm[j].rescount-1]++;
			  if(dist<100)
			    {
			      contacts10[m[0].atm[i].rescount-1]++;
			      contacts10[m[0].atm[j].rescount-1]++;
			      if(dist<64)
				{
				  contacts8[m[0].atm[i].rescount-1]++;
				  contacts8[m[0].atm[j].rescount-1]++;
				}
			    }
			  
			}
		    }
		}
	      
	      //res_contacts[get_res6(m[0].atm[i].residue)][get_res6(m[0].atm[j].residue)]++;
	      //tot_res_contacts++;
	      //}
	    }				     
	}
      
      
      /*RESIDUE CONTACTS for each residue type */
      // printf("RES: ");
      
      for(i=0;i<m[0].atoms;i++)
	{
	  printf("%-5s %-5d %-5d %-5d %-5d %-5d %-5d\n",m[0].atm[i].residue,m[0].atm[i].resnum,contacts8[i],contacts10[i],contacts12[i],contacts14[i],contacts16[i]);
	}
    }
	  //for(j=i;j<6;j++)
	  //  {
	  //	if(i!=j)
	  //	   {
	  //	     res_contacts[i][j]=2*res_contacts[i][j];
	  //	   }
	  //	//restype[i]=restype[i]+res_contacts[i][j];
	  //	temp=0;
	  //	if(tot_res_contacts != 0)
	  //	   {
	  //	     temp=(double)res_contacts[i][j]/tot_res_contacts;
	  //	   }
	  //	printf("%f ",temp);
	  //	sum+=temp;
	  //  }
}



