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
  molecule	m[1];		/* Molecule to be read*/
  int           i,j,k;
  int           atom_contacts[14][14]={{}};
  int           atomtype[14]={};
  int           tot_atom_contacts=0;
  int           res_contacts[6][6]={{}};
  int           restype[6]={};
  int           tot_res_contacts=0;
  int           type_i,type_j;
  int           current_res_i=0;
  int           current_res_j=0;
  double        atom_cutoff=25;  // 5A
  double        res_cutoff=49;   // 7A
  FILE          *fp;
  double temp;
  int tmp=0;
  int           index_i;
  int           index_j;
  /* Parse command line for PDB filename */
  if(argc==2)
    {
      strcpy(m[0].filename,argv[1]);
      //cutoff=strtod(argv[2],(char**)(argv[2]+strlen(argv[2])));
      //cutoff=cutoff*cutoff;
      
    }  
  else
    {
      printf("Usage: read_pdb [pdb_file]\n");
      exit(1);
    }
  if(read_molecules(m)==0)
    {
      for(i=0;i<=m[0].atoms;i++)  
	{
	  //printf("%d %s %s %d\n" ,m[0].atm[i].rescount,m[0].atm[i].name,m[0].atm[i].residue,m[0].atm[i].resnum);
	  
	  for(j=0;j<=m[0].atoms;j++)
	    {
	      if(distance(m,i,j)<atom_cutoff && 
		 abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>1)   /* atom contacts */
		{
		  index_i=get_atomtype(m[0].atm[i].name,m[0].atm[i].residue);
		  index_j=get_atomtype(m[0].atm[j].name,m[0].atm[j].residue);
		  atom_contacts[index_i][index_j]++;
		  //if(index_i<3 && index_j<3)
		  tot_atom_contacts++;
		} 
	    }
	  if(m[0].atm[i].rescount!=current_res_i)
	    {
	      for(j=0,current_res_j=0;j<m[0].atoms;j++)
		{
		  if(m[0].atm[j].rescount!=current_res_j)
		    {
		      //printf("%d %d %d %d %f %f\n",current_res_i+1,current_res_j+1,i,j,crd(i,j),crd(j,i));
		      if(abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>5 &&
			 crd(m,i,j)<res_cutoff)
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
      

      
    
      //printf("ATOM: ");
      for(i=0;i<13;i++)
	{
	  for(j=i;j<13;j++)   // fixing the fact that the diagonal has been counted twice
	    {
	      if(i!=j)
		atom_contacts[i][j]=2*atom_contacts[i][j];
	      temp=0;
	      if(tot_atom_contacts != 0)
		{
		  temp=(double)atom_contacts[i][j]/tot_atom_contacts;
		}
	      printf("%f ",temp);
	      //printf("%d ",atom_contacts[i][j]);
	      // printf("%d ",tot_atom_contacts);
	    }
	}
      printf("\n");
      //printf("RES: ");
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
	      //sum+=temp;
	      //printf("%d ",res_contacts[i][j]);
	    }
	}
      printf("\n");
      printf("%5.3f\n",fatness(m));
      
    }
        
}

