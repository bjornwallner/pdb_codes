#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h"

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
  int           contacts[MAXRES][MAXRES]={{}};
  //  int           res_contacts[20][20]={{}};
  //  int           restype[20]={};
  //  int           tot_res_contacts=0;
  /// int           type_i,type_j;
  int           current_res_i=0;
  int           current_res_j=0;
  int           first=0;
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
      printf("Usage: contact_map [pdb_file] cutoff\n");
      exit(1);
    }
  if(read_molecules(m)==0)
    {  
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
		      if(contacts[current_res_i][current_res_j]==0 && crd(m,i,j)<cutoff)
			{
			  //printf("%d %d %d %d %f %f\n",m[0].atm[j].rescount,m[0].atm[i].rescount,i,j,crd(m,i,j),crd(m,j,i));
			  contacts[current_res_i][current_res_j]=1;
			  contacts[current_res_j][current_res_i]=1;
			  //res_contacts[get_res(m[0].atm[i].residue)][get_res(m[0].atm[j].residue)]++;
			  //tot_res_contacts++;
			}
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
      
      printf("contacts:=[");
     for(i=0;i<m[0].residues;i++)
	{
	  //printf("\n%d: ",i+1);
	  first=0;
	  //printf("(");
	  for(j=0;j<m[0].residues;j++)
	    {
	      if(contacts[i][j]==1)
		{
		  if(first==0)
		    {
		      printf("[%d",j+1);
		      first=1;
		    }
		  else
		    {
		      printf(", %d",j+1);
		    }
		}
	    }
	  if(i+1!=m[0].residues)
	    {
	      printf("],\n");
	      //printf("\n%d:\n ",m[0].residues);
	    }
	  else
	    {
	      printf("]");
	    }
	}
     printf("]:");
     
     printf("\nPDBseq:=\'%s\':\n",m[0].sequence);
    }
}

