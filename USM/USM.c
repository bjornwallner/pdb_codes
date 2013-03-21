#include <stdio.h>
#include <math.h>
#include <stdlib.h>
//#include "molecule.h"
#define	MAXRES		1000		/* Maximum allowable res */
#define	MAXATMS		10000		/* Maximum allowable atoms */
#define TRUE		1		/* Boolean definitions */
#define FALSE		0
#define PI		3.14159265	/* Useful constant */
//#ifdef  NOZLIB
//#else
#include <zlib.h> 
//#endif
//#include <sys/types.h>
//#include <sys/times.h>
//#include <sys/param.h>
//#include <sys/time.h> 
//
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
  molecule      m[2];
  int           i,j,k;
  double        dist;
  double        mean=0;
  int           atoms=0;
  int           residues=0;
  int           contacts[MAXRES][MAXRES]={{}};
  int           contacts2[MAXRES][MAXRES]={{}};
  //  int           res_contacts[20][20]={{}};
  //  int           restype[20]={};
  //  int           tot_res_contacts=0;
  /// int           type_i,type_j;
  int           current_res_i=0;
  int           current_res_j=0;
  long           selected_resnum;
  int           first=0;
  double        cutoff;
  FILE          *fp1;
  FILE          *fp2;
  int    sep=2;
  double temp;
  int tmp=0;
  int random=0;
  char tmpfile1[1000];
  char tmpfile2[1000];
  char tmpfile3[1000];
  char outfile_base[1000];
  //float        sum=0;
  /* Parse command line for PDB filename */
  if(argc==6)
    {
      strcpy(m[0].filename,argv[1]);
      strcpy(m[1].filename,argv[2]);
      strcpy(outfile_base,argv[3]);
      sep=atoi(argv[4]);
      cutoff=strtod(argv[5],(char**)(argv[5]+strlen(argv[5])));
      cutoff=cutoff*cutoff;
      //      printf("%s %d\n",argv[3],sep);
    }
  else
    {
      printf("Usage: USM [pdb_file1] [pdb_file2] [outfile_base] separation cutoff\n");
      exit(1);
    }

  strcat(tmpfile1,outfile_base);
  strcat(tmpfile2,outfile_base);
  strcat(tmpfile3,outfile_base);
  strcat(tmpfile1,"1");
  strcat(tmpfile2,"2");
  strcat(tmpfile3,"3");
  
  //random=rand();
  //printf("%d\n",random);
  //strcat(tmpfile1,random);
  //strcat(tmpfile2,random);
  //  strcat(tmpfile3,random);
  
  //  printf("%s\n%s\n%s\n",tmpfile1,tmpfile2,tmpfile3);
  if(read_molecules(&m[0])==0) {  
      for(i=0;i<m[0].atoms;i++)	{
	  if(m[0].atm[i].rescount!=current_res_i) {
	    for(j=i,current_res_j=current_res_i;j<m[0].atoms;j++) {

	      //printf("%d %d %d %d %f %f\n",current_res_j,current_res_i,i,j,crd(m,i,j),crd(m,j,i));
		  if(m[0].atm[j].rescount!=current_res_j) {
		    if(contacts[current_res_i][current_res_j]==0 && crd(&m[0],i,j)<cutoff) {
		      contacts[current_res_i][current_res_j]=1;
		      contacts[current_res_j][current_res_i]=1;
		    }
		    current_res_j++;
		  }
	    }				     
	    current_res_i++;
	  }
      }
      current_res_i=0;
      current_res_j=0;
      if(read_molecules(&m[1])==0) {  
	for(i=0;i<m[1].atoms;i++) {
	  if(m[1].atm[i].rescount!=current_res_i) {
	    for(j=i,current_res_j=current_res_i;j<m[1].atoms;j++)	{
	      
	      if(m[1].atm[j].rescount!=current_res_j) {
		if(contacts2[current_res_i][current_res_j]==0 && crd(&m[1],i,j)<cutoff) {
		  contacts2[current_res_i][current_res_j]=1;
		  contacts2[current_res_j][current_res_i]=1;
		}
		current_res_j++;
	      }
	    }				     
	    current_res_i++;
	  }
	}
      }

      fp1=fopen(tmpfile1,"w");
      //  fp2=fopen(tmpfile3,"w");
      for(i=0;i<m[0].residues;i++) {
	for(j=0;j<m[0].residues;j++) {
	  if(abs(i-j)>=sep && contacts[i][j]==1) {
	    fprintf(fp1,"%d %d\n",m[0].atm[m[0].CA_ref[i]].resnum,m[0].atm[m[0].CA_ref[j]].resnum);
	    //  fprintf(fp2,"%d %d\n",m[0].atm[m[0].CA_ref[i]].resnum,m[0].atm[m[0].CA_ref[j]].resnum);
	  }
	}
      }
      fclose(fp1);
      fp1=fopen(tmpfile2,"w");

      for(i=0;i<m[1].residues;i++) {
	for(j=0;j<m[1].residues;j++) { 
	  if(abs(i-j)>=sep && contacts[i][j]==1) {
	    fprintf(fp1,"%d %d\n",m[1].atm[m[1].CA_ref[i]].resnum,m[1].atm[m[1].CA_ref[j]].resnum);
	    //    fprintf(fp2,"%d %d\n",m[1].atm[m[1].CA_ref[i]].resnum,m[1].atm[m[1].CA_ref[j]].resnum);
	  }
	}
      }
      fclose(fp1);
      // fclose(fp2);
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
	  //sscanf(buff,"%s %d %s %s %s",line_flag,&number,name,residue,resname);
	  sscanf(buff,"%s",line_flag);
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
	      if(strcmp("CA ",name) == 0)
		{
		  //printf("%s %s %s %d\n",resname,residue,alt_loc,residues);
		  m[0].CA_ref[residues-1]=atoms;
		}
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
  lowest_dist=999999999;
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
