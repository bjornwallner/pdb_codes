#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define	MAXRES		5000		/* Maximum allowable res */
#define	MAXATMS		20000		/* Maximum allowable atoms */
#define TRUE		1		/* Boolean definitions */
#define FALSE		0
#define PI		3.14159265	/* Useful constant */


// Any comments and suggestion may be sent to:
// Author: Björn Wallner
// E-mail: bjorn@sbc.su.se
//   
//


typedef struct {
  struct
  {
    double x,y,z;		/* Atomic coordinates */
    double rms;		/* RMS deviation */
    char residue[8];	/* PDB info for output */
    char name[8];
    int number;
    int resnum;
    int rescount;
    int selected;
  } atm[MAXATMS];
  int CA_ref[MAXRES];
  int res_ref[MAXRES];
  double xcen,ycen,zcen;
  int	atoms;			/* Current # of atoms */
  int   residues;
  char  sequence[MAXRES];
  char	filename[1000];		/* filename to read molecule from */
} molecule;


void strncpy_NULL(char *dest, char *src, size_t n);
int read_molecules(molecule *m);
double crd(molecule *m,int atomno1, int atomno2); 
char aa321(char* res);
double distance(molecule *m,int atomno1,int atomno2);

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
  int contacts2[MAXRES][MAXRES]={{}};
  int tmp=0;
  //float        sum=0;
  /* Parse command line for PDB filename */
  printf("Starting...\n");
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
      //printf("\nPDBseq:=\'%s\':\n",m[0].sequence);
      //exit(1);
      printf("contacts:=[");
     for(i=0;i<m[0].residues;i++)
	{
	  //printf("\n%d: ",i+1);
	  first=0;
	  //printf("(");
	  for(j=0;j<m[0].residues;j++)
	    {
	      if(i==j) //Ensures that the residue is in contact with itself :-)
		contacts[i][j]=1;

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



int read_molecules(molecule *m)	/* Reads in molecules to be superimposed */
{
  int	i,j,k,atoms,residues;	/* Counter variables */
  char	buff[512];	/* Input string */
  char	junk[512];	/* Throwaway string */
  char	line_flag[11];	/* PDB file line mode */
  char	residue[4];	/* PDB atom info for output */
  char	name[4];
  int   resnum;
  char  resname[9];
  char  old_resname[9]="noname";
  char x_temp[11];
  char y_temp[11];
  char z_temp[11];
  //char	number[8];
  int number;
  char chain[2];
  char inscode[2]=" ";
  char alt_loc[2]=" ";
  double	x,y,z;		/* Temporary coordinates values */
  FILE *fp;
  char temp_number[11];
  char temp_resnum[11];
  
  i=0; /* It is only one molecule to be read this was done for the moment instead of changing all "i" to 0 */
  

  //for (i=0;i<2;i++)
  //{
  //#ifdef  NOZLIB*/
  //fp=fopen(m[i].filename,"r");	/* Does file exist? */
  //#else
  //fp=gzopen(m[i].filename,"r");	/* Does file exist? */
  //#endif
  //printf("%s\n",m[0].filename);
  fp=fopen(m[0].filename,"r");	/* Does file exist? */
  if (fp!=NULL)	/* If yes, read in coordinates */
    {
      /* Initialize things */
      //m[0].xcen=m[0].ycen=m[0].zcen=0;
      atoms=0;
      residues=0;
      //#ifdef NOZLIB
      //while(fgets(buff,255,fp)!=NULL)
      //#else
      //while(gzgets(fp,buff,255)!=NULL)
      //#endif
      while(fgets(buff,255,fp)!=NULL)
	{

	  
	  //sscanf(buff,"%-6s %4d  %-3s%1s%3s %1s%5s    %7.3lf %7.3lf %7.3lf",line_flag,&number,name,residue,chain,&resnum,inscode,&x,&y,&z);
	  //printf("%s",buff);
	  //sscanf(buff,"%s %d %s %s %s",line_flag,&number,name,residue,resname);
	  sscanf(buff,"%s",line_flag);
	  if(strcmp("TER",line_flag)==0)
	    break;
	  if (strcmp("ATOM",line_flag)==0)	/* Is it an ATOM entry? */
	    {
	      //printf("%s",&buff[6]);
	      strncpy_NULL(temp_number,&buff[6],5);
	      strncpy_NULL(name,&buff[13],3);
	      strncpy_NULL(alt_loc,&buff[16],1);
	      strncpy_NULL(residue,&buff[17],3);
	      strncpy_NULL(chain,&buff[21],1);
	      strncpy_NULL(temp_resnum,&buff[22],4);
	      strncpy_NULL(resname,&buff[22],5);
	      strncpy_NULL(x_temp,&buff[30],8);
	      strncpy_NULL(y_temp,&buff[38],8);
	      strncpy_NULL(z_temp,&buff[46],8);

	      number=atoi(temp_number);
	      resnum=atoi(temp_resnum);
	      x=atof(x_temp);
	      y=atof(y_temp);
	      z=atof(z_temp);
	      
	      //printf("test: %s %d %s %s %s %s %d %s %lf %lf %lf\n",line_flag,number,name,alt_loc,residue,chain,resnum,resname,x,y,z);
	      //if (strcmp("N",name)==0)	   /*  Is it an N atom => new residue? */
	      //printf("%s %s\n",old_resname,resname);
	      if(strcmp(old_resname,resname)!=0)
		{
		  m[0].sequence[residues]=aa321(residue);
		  residues++;
		  //printf("%s %s\n",resname,residue);
		}
	      //sscanf(&buff[22],"%d %lf %lf %lf",&resnum,&x,&y,&z);
	      //printf("test %d %s %s %d %lf %lf %lf\n",number,name,residue,resnum,x,y,z);
	      m[0].atm[atoms].x=x;
	      m[0].atm[atoms].y=y;
	      m[0].atm[atoms].z=z;
	      m[0].atm[atoms].resnum=resnum;
	      m[0].atm[atoms].number=number; //atoi(number);
	      m[0].atm[atoms].rescount=residues;
	      m[0].atm[atoms].selected=TRUE;
	      strcpy(m[0].atm[atoms].name,name);
	      strcpy(m[0].atm[atoms].residue,residue);
	      //if(strcmp(old_resname,resname)!=0)
	      //	{
	      //	  m[0].sequence[residues-1]=aa321(residue);
	      //	}

	      m[0].xcen+=x;
	      m[0].ycen+=y;
	      m[0].zcen+=z;
	      atoms++;
	      strcpy(old_resname,resname);
	      //}
	    }
	}
      m[0].sequence[residues]='\0';
      m[0].atoms=atoms;
      m[0].residues=residues;
      m[0].xcen=m[0].xcen/atoms;
      m[0].ycen=m[0].ycen/atoms;
      m[0].zcen=m[0].zcen/atoms;


      fclose(fp);
      //#ifdef NOZLIB
      //  fclose(fp);
      //#else
      //	  gzclose(fp);
      //#endif
	  //	  if (atoms!=m[0].atoms)		/* Are file sizes indentical? */
	  //{
	  //  printf("Inconsistent number of atoms in file %s\n",m[0].filename);
	  //  return(1);
	  //}
       
    }
  else
    {
      printf("Couldn't open file %s\n",m[0].filename);
      exit(1);
    }
  return(0);
}

double crd(molecule *m,int atomno1, int atomno2)       /* Closed Residues Distance atomnoX is the first atom of the residue */
{
  int i,j;
  double dist,lowest_dist;
  lowest_dist=999999;
  for(i=atomno1;m[0].atm[i].rescount == m[0].atm[atomno1].rescount;i++)  
    {
      if(strcmp("C  ",m[0].atm[i].name)!=0 && 
	 strcmp("O  ",m[0].atm[i].name)!=0 &&
	 strcmp("N  ",m[0].atm[i].name)!=0)
	{
	  for(j=atomno2;m[0].atm[j].rescount == m[0].atm[atomno2].rescount;j++)
	    {
	      if(strcmp("C  ",m[0].atm[j].name)!=0 && 
		 strcmp("O  ",m[0].atm[j].name)!=0 &&
		 strcmp("N  ",m[0].atm[j].name)!=0)
		{
		  dist=distance(m,i,j);
		  //printf("%f ",dist);
		  if(dist<lowest_dist)
		    lowest_dist=dist;
		  //printf("%s %s %f\n",m[0].atm[j].residue,m[0].atm[j].name,lowest_dist);
		}
	    }
	}
    }
  //printf("%s %s %d %d\n",m[0].atm[i].residue,m[0].atm[i].name,m[0].atm[i].rescount,m[0].atm[atomno1].rescount);
  //printf("%s %s %d %d\n",m[0].atm[j].residue,m[0].atm[j].name,m[0].atm[j].rescount,m[0].atm[atomno2].rescount);
  return lowest_dist;
}

double distance(molecule *m,int atomno1,int atomno2)
{
  //printf("%s %s\n",m[0].atm[atomno1-1].name,m[0].atm[atomno2-1].name);

  return (m[0].atm[atomno1].x-m[0].atm[atomno2].x)*(m[0].atm[atomno1].x-m[0].atm[atomno2].x)+(m[0].atm[atomno1].y-m[0].atm[atomno2].y)*(m[0].atm[atomno1].y-m[0].atm[atomno2].y)+(m[0].atm[atomno1].z-m[0].atm[atomno2].z)*(m[0].atm[atomno1].z-m[0].atm[atomno2].z);
}

char aa321(char* res)
{
  
  if(strcmp("ALA",res)==0)
    return 'A';
  if(strcmp("ARG",res)==0)
    return 'R';
  if(strcmp("ASN",res)==0)
    return 'N';
  if(strcmp("ASP",res)==0)
    return 'D';
  if(strcmp("CYS",res)==0)
    return 'C';
  if(strcmp("GLN",res)==0)
    return 'Q';
  if(strcmp("GLU",res)==0)
    return 'E';
  if(strcmp("GLY",res)==0)
    return 'G';
  if(strcmp("HIS",res)==0)
    return 'H';
  if(strcmp("ILE",res)==0)
    return 'I';
  if(strcmp("LEU",res)==0)
    return 'L';
  if(strcmp("LYS",res)==0)
    return 'K';
  if(strcmp("MET",res)==0)
    return 'M';
  if(strcmp("PHE",res)==0)
    return 'F';
  if(strcmp("PRO",res)==0)
    return 'P';
  if(strcmp("SER",res)==0)
    return 'S';
  if(strcmp("THR",res)==0)
    return 'T';
  if(strcmp("TRP",res)==0)
    return 'W';
  if(strcmp("TYR",res)==0)
    return 'Y';
  if(strcmp("VAL",res)==0)
    return 'V';
}


void strncpy_NULL(char *dest, char *src, size_t n)
{
  strncpy(dest, src, n);
  dest[n]='\0';
}
