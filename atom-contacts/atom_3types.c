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

#define	MAXATMS		10000		/* Maximum allowable atoms */
#define PI		3.14159265	/* Useful constant */


//molecule	m[1];		/* Molecule to be read*/
//double		LG_pvalue();
//double		Levitt_Gerstein();
//double		superimpose_molecules();

//int             get_type3(char *name, char *res);
//int             get_type(char *name, char *res);
//double          distance(int atomno1, int atomno2);
//void            print_type(int type_no, FILE *fp);
//int             get_res(char *res);
//void            print_res(int res,FILE *fp);
//double          crd(int atomno1, int atomno2);   /*closest residue distance */


main(argc,argv)		/* Main routine */
     int argc;
     char *argv[];
{
  molecule	m[1];		/* Molecule to be read*/
  int           i,j,k;
  double        dist;
  double        mean=0;
  int           atoms=0;
  int           residues=0;
  int           atom_contacts[4][4]={{}};
  int           res_contacts[6][6]={{}};
  int           atomtype[3]={};
  int           restype[6]={};
  int           tot_atom_contacts=0;
  int           type_i,type_j;
  int           current_res_i=0;
  int           current_res_j=0;
  double        cutoff;
  FILE          *fp;
  double temp;
  int tmp=0;
  int           index_i;
  int           index_j;
  
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
	  //printf("%d %s %s %d\n" ,m[0].atm[i].rescount,m[0].atm[i].name,m[0].atm[i].residue,m[0].atm[i].resnum);
	  
	  for(j=0;j<m[0].atoms;j++)
	    {
	      if(distance(m,i,j)<cutoff && 
		 abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>1)   /* atom contacts */
		{
		  index_i=get_atomtype3(m[0].atm[i].name,m[0].atm[i].residue);
		  index_j=get_atomtype3(m[0].atm[j].name,m[0].atm[j].residue);
		  atom_contacts[index_i][index_j]++;
		  if(index_i<3 && index_j<3)
		    tot_atom_contacts++;
		} 
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
      
      //printf("%d\n", tot_res_contacts);
      //printf("\tALA    ARG    ASN    ASP    CYS    GLN    GLU    GLY    HIS    ILE    LEU    LYS    MET    PHE    PRO    SER    THR    TRP    TYR    VAL\n");
      //ATOMKONTAKTER
      //k=0;
      //Calculate totalt number of contacts per atom type
      //printf("ATOMTYPE: ");
      printf("ATOM: ");
      for(i=0;i<3;i++)
	{
	  for(j=i;j<3;j++)   // fixing the fact that the diagonal has been counted twice
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
    }
}

