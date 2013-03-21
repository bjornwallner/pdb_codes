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
  molecule      m[1];
  int           i,j,k;
  double        dist;
  double        mean=0;
  int           atoms=0;
  int           residues=0;
  int           arg_number_of_types=20;
  int           res_contacts[20][20]={{}};
  int           contacts[3000][3000]={{}};
char          outfile[10000]="stdout";
  int           restype[20]={};
  int           tot_res_contacts=0;
  int           type_i,type_j;
  int           current_res_i=0;
  int           current_res_j=0;
  double        cutoff=49;
  FILE          *fp;
   double temp;
  int tmp=0;
  //float        sum=0;
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
  
  if(strlen(m[0].filename)==0)
    exit(1);
  fp=fopen(m[0].filename,"r");
  if(fp==NULL || 
     (arg_number_of_types !=6 && arg_number_of_types!=20))
    {
      if(arg_number_of_types !=6 && arg_number_of_types!=20)
	fprintf(stderr,"Error -t %d is not a valid number of types (choose 3 or 13)\n",arg_number_of_types);
      if(fp==NULL)
	fprintf(stderr,"Cannot open %s\n",m[0].filename);
      fprintf(stderr,"Usage: res_contacts -pdb [pdb_file] -c [cutoff default (7Å)]  -t [types 6/20 (default)] -o [outfile]\n");
      exit(1);
    }
  fclose(fp);
  printf("%-11s: %s\n","pdbfile",m[0].filename);
  printf("%-11s: %d\n","residue types",arg_number_of_types);
  printf("%-11s: %8.5lf Å\n","cutoff",sqrt(cutoff));
  printf("%-11s: %s\n","outfile",outfile);

  //if(argc==3)
  //  {
  //    strcpy(m[0].filename,argv[1]);
  //    cutoff=strtod(argv[2],(char**)(argv[2]+strlen(argv[2])));
  //    cutoff=cutoff*cutoff;
  //  }
  //else
  //  {
  //    printf("Usage: read_pdb [pdb_file] cutoff\n");
  //    exit(1);
  //  }
  if(read_molecules(m,'a')==0)
    {  
      if(arg_number_of_types==20)
	{
	  for(i=0;i<m[0].atoms;i++)  
	    {
	      if(m[0].atm[i].rescount!=current_res_i)
		{
		  
		  for(j=0,current_res_j=0;j<m[0].atoms;j++)
		    {
		      if(m[0].atm[j].rescount!=current_res_j)
			{
			  //printf("%d %d %d %d %f %f\n",current_res_i,current_res_j,i,j,crd(m,i,j),crd(m,j,i));
			  if(crd(m,i,j)<cutoff)
			    {
			      
			      contacts[current_res_i][current_res_j]=1;
			      contacts[current_res_j][current_res_i]=1;
			      if(abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>5)
				{
				  res_contacts[get_res(m[0].atm[i].residue)][get_res(m[0].atm[j].residue)]++;
				  tot_res_contacts++;
				}
			    }
			  current_res_j++;
			}
		    }				     
		  current_res_i++;
		}
	    }
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
	      // temp=(double)tot_res_contacts/2;
	      //tot_res_contacts2=(int)temp;
	      fprintf(fp,"CONTACTS %-8d\n",(int)(double)tot_res_contacts/2);
	      fprintf(fp,"CONTACT MAP [seq_nr pdb_nr res [contacting seq_nr]]\n");
	      for(i=0;i<m[0].residues;i++)
		{
		  fprintf(fp,"MAP %4d %4d %c ",i+1,m[0].atm[m[0].CA_ref[i]].resnum,m[0].sequence[i]);
		  for(j=0;j<m[0].residues;j++)
		    {
		      if(contacts[i][j]==1)
			{
			  //printf("%d ",m[0].atm[m[0].CA_ref[j]].resnum);
			  fprintf(fp,"%4d ",j+1); 
			}
		    }
		  fprintf(fp,"\n");
		}
	      // fprintf(fp,"\tC\tN\tO\tCA\tCH3\tCH/CH2\tC(OO-)\tNH\tNH2\t(C)00-\t=0\tOH\tS\n");
	      //printf("\tALA    ARG    ASN    ASP    CYS    GLN    GLU    GLY    HIS    ILE    LEU    LYS    MET    PHE    PRO    SER    THR    TRP    TYR    VAL\n");
	      fprintf(fp,"%-6s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s","COUNT", "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL");
	      fprintf(fp,"\n");
	      k=0;
	      for(i=0;i<20;i++)
		{
		  fprintf(fp,"COUNT ");
		  // print_type(i,fp);
		  //fprintf(fp,"\t");
		  for(j=0;j<20;j++)
		    {
		      temp=0;
		      if(j>=i)
			{
			  if(i==j)
			    {
			      temp=(double)0.5*res_contacts[i][j];
			    }
			  else
			    {
			      temp=(double)res_contacts[i][j];
			    }
			}
		      fprintf(fp,"%8.0lf",temp);
		      //k=k+res_contacts[i][j];
		    }
		  fprintf(fp,"    ");
		  print_res(i,fp);
		  fprintf(fp,"\n");
		  
		}
	      //printf("%d %d %f\n",k,tot_res_contacts,0.5*tot_res_contacts);
	      fprintf(fp,"%-5s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s","FRAC", "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL");
	      fprintf(fp,"\n");
	      for(i=0;i<20;i++)
		{
		  fprintf(fp,"FRAC ");
		  //print_type(i,fp);
		  //fprintf(fp,"\t");
		  for(j=0;j<20;j++)
		    {
		      temp=0;
		      if(tot_res_contacts != 0 && j>=i)
			{
			  if(i!=j)
			    {
			      temp=(double)2*res_contacts[i][j]/tot_res_contacts;
			    }
			  else
			    {
			      temp=(double)res_contacts[i][j]/tot_res_contacts;
			    }
			}
		      fprintf(fp,"%9.6f",temp);
		    }
		  fprintf(fp,"    ");
		  print_res(i,fp);
		  fprintf(fp,"\n");
		  //printf("Hello!\n");
		}
	      fprintf(fp,"RES: ");
	      for(i=0;i<20;i++)
		{
		  for(j=i;j<20;j++)
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
		      fprintf(fp,"%f ",temp);
		      //      sum+=temp;
		    }
		}
	      fprintf(fp,"\n");
	    }
	  /*RESIDUE CONTACTS for each residue type */
	  //intf("RES: ");
	  //r(i=0;i<20;i++)
	  //{
	  //  for(j=i;j<20;j++)
	  //	{
	  //	  if(i!=j)
	  //	    {
	  //	      res_contacts[i][j]=2*res_contacts[i][j];
	  //	    }
	  //	  
	  //	  //restype[i]=restype[i]+res_contacts[i][j];
	  //	  temp=0;
	  //	  if(tot_res_contacts != 0)
	  //	    {
	  //	      temp=(double)res_contacts[i][j]/tot_res_contacts;
	  //	    }
	  //	  printf("%f ",temp);
	  //	  //      sum+=temp;
	  //	  //printf("%d ",res_contacts[i][j]);
	  //	}
	  //}
	  //  printf("\n%f\n",sum);
	}
      else
	{
	for(i=0;i<m[0].atoms;i++)  
	    {
	      if(m[0].atm[i].rescount!=current_res_i)
		{
		  
		  for(j=0,current_res_j=0;j<m[0].atoms;j++)
		    {
		      if(m[0].atm[j].rescount!=current_res_j)
			{
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
	      // temp=(double)tot_res_contacts/2;
	      //tot_res_contacts2=(int)temp;
	      fprintf(fp,"CONTACTS %-8d\n",(int)(double)tot_res_contacts/2);
	      // fprintf(fp,"\tC\tN\tO\tCA\tCH3\tCH/CH2\tC(OO-)\tNH\tNH2\t(C)00-\t=0\tOH\tS\n");
	      //printf("\tALA    ARG    ASN    ASP    CYS    GLN    GLU    GLY    HIS    ILE    LEU    LYS    MET    PHE    PRO    SER    THR    TRP    TYR    VAL\n");
	      fprintf(fp,"%-6s%8s%8s%8s%8s%8s%8s","COUNT", "+","-","aromat","polar","hydrop","special");
	      fprintf(fp,"\n");
	      k=0;
	      for(i=0;i<6;i++)
		{
		  fprintf(fp,"COUNT ");
		  // print_type(i,fp);
		  //fprintf(fp,"\t");
		  for(j=0;j<6;j++)
		    {
		      temp=0;
		      if(j>=i)
			{
			  if(i==j)
			    {
			      temp=(double)0.5*res_contacts[i][j];
			    }
			  else
			    {
			      temp=(double)res_contacts[i][j];
			    }
			}
		      fprintf(fp,"%8.0lf",temp);
		      //k=k+res_contacts[i][j];
		    }
		  fprintf(fp,"    ");
		  //print_res(i,fp);
		  fprintf(fp,"\n");
		  
		}
	      //printf("%d %d %f\n",k,tot_res_contacts,0.5*tot_res_contacts);
	      fprintf(fp,"%-5s%9s%9s%9s%9s%9s%9s","FRAC", "+","-","aromat","polar","hydrop","special");
	      fprintf(fp,"\n");
	      for(i=0;i<6;i++)
		{
		  fprintf(fp,"FRAC ");
		  //print_type(i,fp);
		  //fprintf(fp,"\t");
		  for(j=0;j<6;j++)
		    {
		      temp=0;
		      if(tot_res_contacts != 0 && j>=i)
			{
			  if(i!=j)
			    {
			      temp=(double)2*res_contacts[i][j]/tot_res_contacts;
			    }
			  else
			    {
			      temp=(double)res_contacts[i][j]/tot_res_contacts;
			    }
			}
		      fprintf(fp,"%9.6f",temp);
		    }
		  fprintf(fp,"    ");
		  //print_res(i,fp);
		  fprintf(fp,"\n");
		  //printf("Hello!\n");
		}
	      fprintf(fp,"RES: ");
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
		      fprintf(fp,"%f ",temp);
		      //      sum+=temp;
		    }
		}
	       fprintf(fp,"\n");
	      //printf("\n%d\n",sum);
  
	    }
	}
    }
}



