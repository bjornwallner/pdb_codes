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
  int           arg_number_of_types=13;
  int           atom_contacts[14][14]={{}};
  int           atomtype[14]={};
  int           tot_atom_contacts=0,tot_atom_contacts2=0;
  int           type_i,type_j;
  int           current_res_i=0;
  int           current_res_j=0;
  double        cutoff=25;
  FILE          *fp;
  double temp;
  int tmp=0;
  int           index_i;
  int           index_j;
  char          outfile[10000]="stdout";
  //char**        endptr;
  /* Parse command line for PDB filename */
  i=1;
/*    temp=1; */
  while (i<argc)
    {
      //printf("%s\n",argv[i]);
      if (strcmp(argv[i],"-pdb")==0)
	{
	  strcpy(m[0].filename,argv[i+1]);
	  i++;
	}
      else if (strcmp(argv[i],"-c")==0)
	{
	  //printf("test: %s %d\n",argv[i+1],strlen(argv[i+1]));
	  //cutoff=strtod(argv[i+1],(char**)(argv[i+1]+strlen(argv[i+1])));
	  cutoff=strtod(argv[i+1],NULL); //argv[i+1]+strlen(argv[i+1]));
	  cutoff=cutoff*cutoff;
	  //  printf("test: %s %d %f\n",argv[i+1],strlen(argv[i+1]),cutoff);
	  i++;
	}
      else if (strcmp(argv[i],"-t")==0)
	{
	  
	  arg_number_of_types=atoi(argv[i+1]);
	  i++;
	  //cutoff=strtod(argv[i+1],(char**)(argv[i+1]+strlen(argv[i+1])));
	  //cutoff=cutoff*cutoff;
       	}
      else if (strcmp(argv[i],"-o")==0)
	{
	  strcpy(outfile,argv[i+1]);
	  i++;
	  //arg_number_of_types=atoi(argv[i+1]);
	  //cutoff=strtod(argv[i+1],(char**)(argv[i+1]+strlen(argv[i+1])));
	  //cutoff=cutoff*cutoff;
       	}
      i++;
    }
  
  fp=fopen(m[0].filename,"r");
  if(fp==NULL || 
     (arg_number_of_types !=3 && arg_number_of_types!=13))
    {
      if(arg_number_of_types !=3 && arg_number_of_types!=13)
	printf("Error -t %d is not a valid number of types (choose 3 or 13)\n",arg_number_of_types);
      if(fp==NULL)
	fprintf(stderr,"Cannot open %s\n",m[0].filename);
      printf("Usage: atom_contacts -pdb [pdb_file] -c [cutoff default (5Å)]  -t [types 3/13 (default)]\n");
      exit(1);
    }
  fclose(fp);
  printf("%-11s: %s\n","pdbfile",m[0].filename);
  printf("%-11s: %d\n","atom types",arg_number_of_types);
  printf("%-11s: %8.5lf Å\n","cutoff",sqrt(cutoff));
  printf("%-11s: %s\n","outfile",outfile);
  if(read_molecules(m,'a')==0)
    {

      if(arg_number_of_types==13)
	{
	  for(i=0;i<=m[0].atoms;i++)  
	    {
	      //printf("%d %s %s %d\n" ,m[0].atm[i].rescount,m[0].atm[i].name,m[0].atm[i].residue,m[0].atm[i].resnum);
	      
	      for(j=0;j<=m[0].atoms;j++)
		{
		  if(distance(m,i,j)<cutoff && 
		     abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>1)   /* atom contacts */
		    {
		      index_i=get_atomtype(m[0].atm[i].name,m[0].atm[i].residue);
		      index_j=get_atomtype(m[0].atm[j].name,m[0].atm[j].residue);
		      atom_contacts[index_i][index_j]++;
		      if(index_i<arg_number_of_types && index_j<arg_number_of_types)
			tot_atom_contacts++;
		    } 
		}
	    }   
	  // printf("%f %f\n",crd(0,16),crd(16,0));
	  //fp=fopen("test","w");	/* Does file exist? */
	  if(strcmp(outfile,"stdout")==0)
	    {
	      fp=stdout;
	    }
	  else
	    {
	      fp=fopen(outfile,"w");
	      if (fp==NULL)	
		{
		  fprintf(stderr,"Cannot open outfile: %s\n",outfile);
		  exit(1);
		}
	    }
	  if (fp!=NULL)	/* If yes, write output */
	    {
	      fprintf(fp,"CUTOFF %-8.5f\n",sqrt(cutoff));
	      fprintf(fp,"LENGTH %-8d\n",m[0].residues);
	      // temp=(double)tot_atom_contacts/2;
	      //tot_atom_contacts2=(int)temp;
	      fprintf(fp,"CONTACTS %-8d\n",(int)(double)tot_atom_contacts/2);
	      // fprintf(fp,"\tC\tN\tO\tCA\tCH3\tCH/CH2\tC(OO-)\tNH\tNH2\t(C)00-\t=0\tOH\tS\n");
	      fprintf(fp,"%-6s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s","COUNT", "C","N","O","CA","CH3","CH/CH2","C(OO-)","NH","NH2","(C)OO-","=O","OH","S");
	      fprintf(fp,"\n");
	      k=0;
	      for(i=0;i<13;i++)
		{
		  fprintf(fp,"COUNT ");
		  // print_type(i,fp);
		  //fprintf(fp,"\t");
		  for(j=0;j<13;j++)
		    {
		      temp=0;
		      if(j>=i)
			{
			  if(i==j)
			    {
			      temp=(double)0.5*atom_contacts[i][j];
			    }
			  else
			    {
			      temp=(double)atom_contacts[i][j];
			    }
			}
		      fprintf(fp,"%8.0lf",temp);
		      k=k+atom_contacts[i][j];
		    }
		  fprintf(fp,"    ");
		  print_type(i,fp);
		  fprintf(fp,"\n");
		  
		}
	      printf("%d %d %f\n",k,tot_atom_contacts,0.5*tot_atom_contacts);
	      fprintf(fp,"%-5s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s","FRAC", "C","N","O","CA","CH3","CH/CH2","C(OO-)","NH","NH2","(C)OO-","=O","OH","S");
	      fprintf(fp,"\n");
	      for(i=0;i<13;i++)
		{
		  fprintf(fp,"FRAC ");
		  //print_type(i,fp);
		  //fprintf(fp,"\t");
		  for(j=0;j<13;j++)
		    {
		      temp=0;
		      if(tot_atom_contacts != 0 && j>=i)
			{
			  if(i!=j)
			    {
			      temp=(double)2*atom_contacts[i][j]/tot_atom_contacts;
			    }
			  else
			    {
			      temp=(double)atom_contacts[i][j]/tot_atom_contacts;
			    }
			}
		      fprintf(fp,"%8.5f",temp);
		    }
		  fprintf(fp,"    ");
		  print_type(i,fp);
		  fprintf(fp,"\n");
		  
		}
	    }
	  
	  
	  // fclose(fp);
	  
	  //printf("%d\n", tot_res_contacts);
	  //printf("\tALA    ARG    ASN    ASP    CYS    GLN    GLU    GLY    HIS    ILE    LEU    LYS    MET    PHE    PRO    SER    THR    TRP    TYR    VAL\n");
	  //ATOMKONTAKTER
	  //k=0;
	  //Calculate totalt number of contacts per atom type
	  //printf("ATOMTYPE: ");
	  //exit(1);
	  fprintf(fp,"ATOM: ");
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
		  fprintf(fp,"%f ",temp);
		  //printf("%d ",atom_contacts[i][j]);
		  // printf("%d ",tot_atom_contacts);
		}
	    }
	  fprintf(fp,"\n");
	  fclose(fp);
	}
      else
	{
	  for(i=0;i<=m[0].atoms;i++)  
	    {
	      //printf("%d %s %s %d\n" ,m[0].atm[i].rescount,m[0].atm[i].name,m[0].atm[i].residue,m[0].atm[i].resnum);
	      
	      for(j=0;j<=m[0].atoms;j++)
		{
		  if(distance(m,i,j)<cutoff && 
		     abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>1)   /* atom contacts */
		    {
		      index_i=get_atomtype3(m[0].atm[i].name,m[0].atm[i].residue);
		      index_j=get_atomtype3(m[0].atm[j].name,m[0].atm[j].residue);
		      atom_contacts[index_i][index_j]++;
		      if(index_i<arg_number_of_types && index_j<arg_number_of_types)
			tot_atom_contacts++;
		    } 
		}
	    }   
	  // printf("%f %f\n",crd(0,16),crd(16,0));
	  //fp=fopen("test","w");	/* Does file exist? */
	  if(strcmp(outfile,"stdout")==0)
	    {
	      fp=stdout;
	    }
	  else
	    {
	      fp=fopen(outfile,"w");
	      if (fp==NULL)	
		{
		  fprintf(stderr,"Cannot open outfile: %s\n",outfile);
		  exit(1);
		}
	    }
	  if (fp!=NULL)	/* If yes, write output */
	    {
	      fprintf(fp,"CUTOFF %-8.5f\n",sqrt(cutoff));
	      fprintf(fp,"LENGTH %-8d\n",m[0].residues);
	      // temp=(double)tot_atom_contacts/2;
	      //tot_atom_contacts2=(int)temp;
	      fprintf(fp,"CONTACTS %-8d\n",(int)(double)tot_atom_contacts/2);
	      // fprintf(fp,"\tC\tN\tO\tCA\tCH3\tCH/CH2\tC(OO-)\tNH\tNH2\t(C)00-\t=0\tOH\tS\n");
	      fprintf(fp,"%-5s%8s%8s%8s","COUNT", "C","N","O");
	      fprintf(fp,"\n");
	      k=0;
	      for(i=0;i<3;i++)
		{
		  fprintf(fp,"COUNT ");
		  // print_type(i,fp);
		  //fprintf(fp,"\t");
		  for(j=0;j<3;j++)
		    {
		      temp=0;
		      if(j>=i)
			{
			  if(i==j)
			    {
			      temp=(double)0.5*atom_contacts[i][j];
			    }
			  else
			    {
			      temp=(double)atom_contacts[i][j];
			    }
			}
		      fprintf(fp,"%8.0lf",temp);
		      k=k+atom_contacts[i][j];
		    }
		  fprintf(fp,"    ");
		  if(i==0)
		    fprintf(fp,"C");
		  if(i==1)
		    fprintf(fp,"N");
		  if(i==2)
		    fprintf(fp,"O");
		  fprintf(fp,"\n");
		  
		}
	      //printf("%d %d %f\n",k,tot_atom_contacts,0.5*tot_atom_contacts);
	      fprintf(fp,"%-5s%8s%8s%8s","FRAC", "C","N","O");
	      fprintf(fp,"\n");
	      for(i=0;i<3;i++)
		{
		  
		  fprintf(fp,"FRAC ");
		  //print_type(i,fp);
		  //fprintf(fp,"\t");
		  for(j=0;j<3;j++)
		    {
		      temp=0;
		      if(tot_atom_contacts != 0 && j>=i)
			{
			  if(i!=j)
			    {
			      temp=(double)2*atom_contacts[i][j]/tot_atom_contacts;
			    }
			  else
			    {
			      temp=(double)atom_contacts[i][j]/tot_atom_contacts;
			    }
			}
		      fprintf(fp,"%8.5f",temp);
		    }
		  fprintf(fp,"    ");
		  if(i==0)
		    fprintf(fp,"C");
		  if(i==1)
		    fprintf(fp,"N");
		  if(i==2)
		    fprintf(fp,"O");
		  fprintf(fp,"\n");
		  
		}
	      fprintf(fp,"ATOM: ");
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
		      fprintf(fp,"%f ",temp);
		  //printf("%d ",atom_contacts[i][j]);
		  // printf("%d ",tot_atom_contacts);
		    }
		}
	      fprintf(fp,"\n");
	      fclose(fp);
	    }
	}
    }
}
        
    
    

