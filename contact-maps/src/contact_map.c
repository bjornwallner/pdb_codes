#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "molecule.h"

//
// Any comments and suggestion may be sent to:
// Author: Björn Wallner
// E-mail: bjorn@sbc.su.se
//   
//


main(argc,argv)		/* Main routine */
     int argc;
     char *argv[];
{
  molecule      m[1];
  int           i,j,k;
  double        d;
  double        mean=0;
  int           atoms=0;
  int           residues=0;
  //  int           contacts[MAXRES][MAXRES]={{}};
  int **contacts;
  double **dist;
  
  //  int           res_contacts[20][20]={{}};
  //  int           restype[20]={};
  //  int           tot_res_contacts=0;
  /// int           type_i,type_j;
  int           current_res_i=0;
  int           current_res_j=0;
  long           selected_resnum;
  int           first=0;
  double        cutoff=0;
  FILE          *fp;
  int error=0;
  double temp;
  int tmp=0;
  int binary=1;
  char chain_name;
  //float        sum=0;
  /* Parse command line for PDB filename */
  if(argc==2)
    {
      strcpy(m[0].filename,argv[1]);
      binary=0;
      //      cutoff=strtod(argv[2],(char**)(argv[2]+strlen(argv[2])));
      //     cutoff=cutoff*cutoff;
    }
  else if(argc==3)
    {
      strcpy(m[0].filename,argv[1]);
      cutoff=strtod(argv[2],(char**)(argv[2]+strlen(argv[2])));
      cutoff=cutoff*cutoff;
    }
  else
    {
      printf("Usage: contact_map [pdb_file] cutoff(for binary otherwise distance)\n");
      exit(1);
    }


  error=read_molecules(m,'a');
  if(error==0)
    {  
      residues=m[0].residues;
      contacts = malloc(residues*sizeof(int *));
      dist = malloc(residues*sizeof(double *));
      for(i=0;i<residues;i++) {
	contacts[i] = malloc(residues*sizeof(int));
	dist[i] = malloc(residues*sizeof(double));
      }
      for(i=0;i<m[0].atoms;i++)  
	{
	  if(m[0].atm[i].rescount!=current_res_i)
	    {
	      for(j=i,current_res_j=current_res_i;j<m[0].atoms;j++)
		{

		  if(m[0].atm[j].rescount!=current_res_j)
		    {
		      //printf("%d %d %d %d %f %f\n",current_res_j,current_res_i,i,j,crd(m,i,j),crd(m,j,i));
		      //	     {
		      //    printf("%d %d %d %d %f %f\n",m[0].atm[j].rescount,m[0].atm[i].rescount,i,j,crd(m,i,j),crd(m,j,i));
		      //printf("%d %s %s %i\n" ,m[0].atm[i].rescount,m[0].atm[i].name,m[0].atm[i].residue,m[0].atm[i].resnum);
		      //printf("%d %s %s %i\n" ,m[0].atm[j].rescount,m[0].atm[j].name,m[0].atm[j].residue,m[0].atm[j].resnum);
		      //printf("%d %d %d %d %f %f\n",m[0].atm[j].rescount,m[0].atm[i].rescount,i,j,crd(m,i,j),crd(m,j,i));
		      //if(abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>5 &&
		      d=crd(m,i,j);
		      if(contacts[current_res_i][current_res_j]==0 && d<cutoff)
			{
			  //printf("%d %d %d %d %f %f\n",m[0].atm[j].rescount,m[0].atm[i].rescount,i,j,crd(m,i,j),crd(m,j,i));
			  contacts[current_res_i][current_res_j]=1;
			  contacts[current_res_j][current_res_i]=1;
			  //res_contacts[get_res(m[0].atm[i].residue)][get_res(m[0].atm[j].residue)]++;
			  //tot_res_contacts++;
			}
		      dist[current_res_i][current_res_j]=d;
		      dist[current_res_j][current_res_i]=d;
		      current_res_j++;
		    }
		}				     
	      current_res_i++;
	    }
	}
      
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
      
    //for(i=0;i<m[0].residues;i++)
    //{
    //	printf("%d: ",i+1); 
    //	for(j=0;j<m[0].residues;j++)
    //	  {
    //	    if(contacts[i][j]==1)
    //	  printf("%d ",j+1);
    //	  }
    //	printf("\n");
    //}
      
      //      printf("contacts:=[");
      if(binary) {
	printf("CUTOFF %lf\n",sqrt(cutoff));
      } else {
	printf("CUTOFF running in distance mode\n");
      }
      printf("FORMAT: RES <counter> <resname> <residue> <chain>: <dist11> <dist12> ...\n");
      //      printf("SEQ %s\n",m[0].sequence);
     for(i=0;i<m[0].residues;i++)
	{
	  // if(m[0].atm[m[0].CA_ref[i]].resnum == selected_resnum)
	    {
	      //printf("RES %d (%d) %c ",m[0].atm[m[0].CA_ref[i]].resnum,i+1,m[0].sequence[i]);
	      chain_name=m[0].atm[m[0].CA_ref[i]].chain[0];
	      if(chain_name == ' ') {
		chain_name='_';
	      }
	      
	      printf("RES %d %s %c %c: ",i+1,m[0].atm[m[0].CA_ref[i]].resname,m[0].sequence[i],chain_name);
	      first=0;
	      //printf("(");
	      for(j=0;j<m[0].residues;j++)
		{
		  
		  if(binary) 
		    {
		      if(contacts[i][j]==1)
			{
			  //  printf("%d (%d) ",m[0].atm[m[0].CA_ref[j]].resnum,j+1);
			  //		      printf("%d ",j+1);
			  printf("%d ",m[0].atm[m[0].CA_ref[j]].resnum);
			  // if(first==0)
			  //   {
			  //     printf("[%d",j+1);
			  //     first=1;
			  //   }
			  // else
			  //   {
			  //     printf(", %d",j+1);
			  //   }
			}
		    } else {
		    printf("%lf ",sqrt(dist[i][j]));
		  }
		}
	      printf("\n");
	    }
	  //if(i+1!=m[0].residues)
	  //  {
	  //    //printf("],\n");
	  //    //printf("\n%d:\n ",m[0].residues);
	  //  }
	  //else
	  //  {
	  //    //printf("]");
	  //  }
	}
     //printf("]:");
     

    }
}

