/*
  Program RMS   written by Scott M. Le Grand
  Modified by Arne Elofsson
  Copyright 1992, Scott M. Le Grand
  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h> 
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
#include <sys/time.h> 

#define	MAXATMS		15000		/* Maximum allowable atoms */
#define	MAXRES       	5000		/* Maximum allowable atoms */
#define PI		3.14159265	/* Useful constant */

#define TRUE		1		/* Boolean definitions */
#define FALSE		0
typedef struct {
  struct
  {
    double x,y,z;		/* Atomic coordinates */
    double rms;		/* RMS deviation */
    char residue[8];	/* PDB info for output */
    char name[8];
    char resname[8];
    int number;
    int resnum;
    int rescount;
    int selected;
   int owner;
  } atm[MAXATMS];
  int CA_ref[MAXRES];
  int res_ref[MAXRES];
  double xcen,ycen,zcen;
  int	atoms;			/* Current # of atoms */
  int   residues;
  char  sequence[MAXRES];
  char  ss[MAXRES];
  double score;
  char method[200];
  int  rank;
  char	filename[1000];		/* filename to read molecule from */
  
} molecule;


molecule	m[2];		/* Molecules to be compared */

double		superimpose_molecules();
double          Sscore();
double		rms;			/* final RMS error */
int		verbose_flag;		/* Verbose output flag */
int             matrix_flag;
int             rotation_flag;
int             output_flag;            /* Output flag */
char            output_name[1000];            /* Output name */
char            rot_output_name[1000];            /* Output name */
int             output_flag_ca;            /* Output flag */
char            output_name_ca[1000];            /* Output name */
int             read_molecules();
int             read_molecules2();
int             read_molecules99();
int             read_molecules_ca99();
int atom_exist(char *name,char *resname,molecule *m);
int center_molecule(molecule *m);
void strncpy_NULL(char *dest, char *src, size_t n);
char aa321(char* res);

/*int             temp;*/

main(argc,argv)		/* Main routine */
    int argc;
    char *argv[];
{
 int	files;	/* Number of files parsed */
 int	i;	/* Counter variable */
 double	s[3][3];	/* Final transformation matrix */
 double s_all[3][3];
 double s_ca[3][3];
 /* Initialize */
 files=0;
 verbose_flag=FALSE;
 output_flag=FALSE;
 matrix_flag=FALSE;
 rotation_flag=FALSE;
 /* Parse command line for PDB files and options */
 i=1;
/*    temp=1; */
 while (i<argc)
  {
   if (strcmp(argv[i],"-ca")==0)
    {
     strcpy(output_name_ca,argv[i+1]);
     output_flag_ca=TRUE;
     i++;
    }
   else if (strcmp(argv[i],"-o")==0)
     {
       strcpy(output_name,argv[i+1]);
       output_flag=TRUE;
       i++;
     }
   else if (strcmp(argv[i],"-v")==0)
     {
       verbose_flag=TRUE;
     }
   else if (strcmp(argv[i],"-m")==0)
     {
       matrix_flag=TRUE;
     }
   else if (strcmp(argv[i],"-rot")==0)
     {
       strcpy(rot_output_name,argv[i+1]);
       rotation_flag=TRUE;
       i++;
     }
   else
     {
       strcpy(m[files].filename,argv[i]);
       files++;
     }
   i++;
  }
 if (files==2)	/* Valid input line? */
   {
     /*
     if (read_molecules99()==0)
       {
	 //printf("hej\n");
	 check_molecules(&m[0],&m[1]); // this deletes and change the order some atoms
	 // printf("hej\n");
	 //exit(1);
	 center_molecule(&m[0]);
	 center_molecule(&m[1]);
	 rms=superimpose_molecules(&m[0],&m[1],s_all);
	 //rms=superimpose_molecules(&m[0],&m[1],s);
	 if(matrix_flag)
	   {
	     //  printf("ALL: ");
	     //printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",s_all[0][0],s_all[0][1],s_all[0][2],s_all[1][0],s_all[1][1],s_all[1][2],s_all[2][0],s_all[2][1],s_all[2][2]);
	   }
	 else
	   {
	     output_results();
	     if (output_flag)
	       output_file(output_name);
	   }
	   }

       */
     //     printf("test0\n");
     if (read_molecules_ca99()==0)
       {
	 check_molecules(&m[0],&m[1]); // this deletes some atom records
	 //	 printf("test1\n");
	 //center_molecule(&m[0]);
	 //center_molecule(&m[1]);
	 //	 rms=superimpose_molecules(&m[0],&m[1],s_ca);
	 rms=Sscore(&m[0],&m[1]);
	 //printf("%lf %lf %lf\n%lf %lf %lf\n%lf %lf %lf\n",s_all[0][0],s_all[0][1],s_all[0][2],s_all[1][0],s_all[1][1],s_all[1][2],s_all[2][0],s_all[2][1],s_all[2][2]);
	 //exit(1);

		 //	 printf("test2\n");
	 strcpy(output_name,m[0].filename);
	 strcat(output_name,".S2");
	 output_file(output_name);
	 
       }
   }
 else
   printf("Usage: rmsd [-m prints the rotation matrices] [-rot filename for rotation of file1 to file2] [ -v ] [ -o outfile] [ -ca outfile] file1 file2\n");
}
int read_molecules99(temp2)	/* Reads in molecules to be superimposed */
     int temp2;
{
  int	i,j,k,atoms,residues;	/* Counter variables */
  char	buff[1200];	/* Input string */
  char	temp[1000];	/* Throwaway string */
  char	line_flag[11];	/* PDB file line mode */
  char  desc_flag[20];
  char	residue[4];	/* PDB atom info for output */
  char	name[4];
  int   resnum;
  char  resname[9];
  char  old_resname[9]="noname";
  char alt_loc_check[2]=" ";
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
  int ss_flag=0;

  //i=0; /* It is only one molecule to be read this was done for the moment instead of changing all "i" to 0 */
  

  //for (i=0;i<2;i++)
  //{
  //#ifdef  NOZLIB*/
  //fp=fopen(m[i].filename,"r");	/* Does file exist? */
  //#else
  //fp=gzopen(m[i].filename,"r");	/* Does file exist? */
  //#endif
  //printf("%s\n",m[0].filename);
  for (i=0;i<2;i++)
    {
      strcpy(m[i].method,"undef");
      m[i].rank=-1;
      m[i].score=-9999;
      
      fp=fopen(m[i].filename,"r");	/* Does file exist? */
      if (fp!=NULL)	/* If yes, read in coordinates */
	{
	  /* Initialize things */
	  //m[i].xcen=m[i].ycen=m[i].zcen=0;
	  atoms=0;
	  residues=0;
	  //#ifdef NOZLIB
      //while(fgets(buff,255,fp)!=NULL)
      //#else
      //while(gzgets(fp,buff,255)!=NULL)
      //#endif
	  while(fgets(buff,1000,fp)!=NULL)
	    {
	      //sscanf(buff,"%-6s %4d  %-3s%1s%3s %1s%5s    %7.3lf %7.3lf %7.3lf",line_flag,&number,name,residue,chain,&resnum,inscode,&x,&y,&z);
	      //sscanf(buff,"%s %d %s %s %s",line_flag,&number,name,residue,resname);
	      sscanf(buff,"%s",line_flag);
	      if(strcmp("REMARK",line_flag)==0)
		{
		  sscanf(buff,"%s %s %s",line_flag,desc_flag,temp);
		  if(strcmp("SS",desc_flag)==0)
		    {
		  strcpy(m[i].ss,temp);
		  ss_flag=1;
		  //printf("%s\n",m[i].ss);
		    }
		  if(strcmp("METHOD",desc_flag)==0)
		    {
		      strcpy(m[i].method,temp);
		      //printf("%s\n",m[i].method);
		    }
		  if(strcmp("SCORE",desc_flag)==0)
		    {
		      m[i].score=atof(temp);
		      //strcpy(m[i].method,temp);
		      //printf("%lf\n",m[i].score);
		    }
		}
	      if(strcmp("MODEL",line_flag)==0)
		{
		  sscanf(buff,"%s %s",line_flag,temp);
		  m[i].rank=atoi(temp);
		}
	      
	      if (strcmp("ATOM",line_flag)==0 && buff[13] != 'H')	/* Is it an ATOM entry? */
		{
		  //printf("%s",&buff[6]);
		  strncpy_NULL(temp_number,&buff[6],5);
		  strncpy_NULL(name,&buff[13],3);
		  strncpy_NULL(alt_loc,&buff[16],1);
		  strncpy_NULL(residue,&buff[17],3);
		  strncpy_NULL(chain,&buff[21],1);
		  strncpy_NULL(temp_resnum,&buff[22],5);
		  strncpy_NULL(resname,&buff[22],5);
		  strncpy_NULL(x_temp,&buff[30],8);
		  strncpy_NULL(y_temp,&buff[38],8);
		  strncpy_NULL(z_temp,&buff[46],8);
		  
		  number=atoi(temp_number);
		  resnum=atoi(temp_resnum);
		  x=atof(x_temp);
		  y=atof(y_temp);
		  z=atof(z_temp);
		  
		  //printf("test: %d %s %d %s %s %s %s %d %s %lf %lf %lf\n",i,line_flag,number,name,alt_loc,residue,chain,resnum,resname,x,y,z);
	      //if (strcmp("N",name)==0)	   /*  Is it an N atom => new residue? */
	      //printf("%s %s\n",old_resname,resname);
		  if(strcmp(old_resname,resname)!=0)
		    {
		      m[i].sequence[residues]=aa321(residue);
		      residues++;
		      strcpy(alt_loc_check,alt_loc);
		      //printf("%s %s\n",resname,residue);
		    }
		  //sscanf(&buff[22],"%d %lf %lf %lf",&resnum,&x,&y,&z);
		  // printf("test %d %s %s %d %lf %lf %lf\n",number,name,residue,resnum,x,y,z);
		  if(strcmp(alt_loc_check,alt_loc)==0)
		    {
		      m[i].atm[atoms].x=x;
		      m[i].atm[atoms].y=y;
		      m[i].atm[atoms].z=z;
		      m[i].atm[atoms].resnum=resnum;
		      strcpy(m[i].atm[atoms].resname,resname);
		      m[i].atm[atoms].number=number; //atoi(number);
		      m[i].atm[atoms].rescount=residues;
		      m[i].atm[atoms].selected=TRUE;
		      strcpy(m[i].atm[atoms].name,name);
		      strcpy(m[i].atm[atoms].residue,residue);
		      //if(strcmp(old_resname,resname)!=0)
		      //	{
		      //	  m[i].sequence[residues-1]=aa321(residue);
		      //	}
		      
		      m[i].xcen+=x;
		      m[i].ycen+=y;
		      m[i].zcen+=z;
		      atoms++;
		      strcpy(old_resname,resname);
		    }
		}
	    }
	  m[i].sequence[residues]='\0';
	  m[i].atoms=atoms;
	  m[i].residues=residues;
	  m[i].xcen=m[i].xcen/atoms;
	  m[i].ycen=m[i].ycen/atoms;
	  m[i].zcen=m[i].zcen/atoms;
	  if(ss_flag == 0)
	    {
	      //fprintf(stderr,"No secondary structure information in file!\n");
	    }
	  
	  fclose(fp);
	  //#ifdef NOZLIB
	  //  fclose(fp);
	  //#else
	  //	  gzclose(fp);
	  //#endif
	  //	  if (atoms!=m[i].atoms)		/* Are file sizes indentical? */
	  //{
	  //  printf("Inconsistent number of atoms in file %s\n",m[i].filename);
	  //  return(1);
	  //}
	  
	}
      else
	{
	  printf("Couldn't open file %s\n",m[i].filename);
	  exit(1);
	}
    }
  return(0);
}
int read_molecules_ca99(temp2)	/* Reads in molecules to be superimposed */
     int temp2;
{
  int	i,j,k,atoms,residues;	/* Counter variables */
  char	buff[1200];	/* Input string */
  char	temp[1000];	/* Throwaway string */
  char	line_flag[11];	/* PDB file line mode */
  char  desc_flag[20];
  char	residue[4];	/* PDB atom info for output */
  char	name[4];
  int   resnum;
  char  resname[9];
  char  old_resname[9]="noname";
  char alt_loc_check[2]=" ";
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
  int ss_flag=0;

  //i=0; /* It is only one molecule to be read this was done for the moment instead of changing all "i" to 0 */
  

  //for (i=0;i<2;i++)
  //{
  //#ifdef  NOZLIB*/
  //fp=fopen(m[i].filename,"r");	/* Does file exist? */
  //#else
  //fp=gzopen(m[i].filename,"r");	/* Does file exist? */
  //#endif
  //printf("%s\n",m[0].filename);
  for (i=0;i<2;i++)
    {
      strcpy(m[i].method,"undef");
      m[i].rank=-1;
      m[i].score=-9999;
      
      fp=fopen(m[i].filename,"r");	/* Does file exist? */
      if (fp!=NULL)	/* If yes, read in coordinates */
	{
	  /* Initialize things */
	  //m[i].xcen=m[i].ycen=m[i].zcen=0;
	  atoms=0;
	  residues=0;
	  //#ifdef NOZLIB
      //while(fgets(buff,255,fp)!=NULL)
      //#else
      //while(gzgets(fp,buff,255)!=NULL)
      //#endif
	  while(fgets(buff,1000,fp)!=NULL)
	    {
	      //sscanf(buff,"%-6s %4d  %-3s%1s%3s %1s%5s    %7.3lf %7.3lf %7.3lf",line_flag,&number,name,residue,chain,&resnum,inscode,&x,&y,&z);
	      //sscanf(buff,"%s %d %s %s %s",line_flag,&number,name,residue,resname);
	      sscanf(buff,"%s",line_flag);
	      if(strcmp("REMARK",line_flag)==0)
		{
		  sscanf(buff,"%s %s %s",line_flag,desc_flag,temp);
		  if(strcmp("SS",desc_flag)==0)
		    {
		  strcpy(m[i].ss,temp);
		  ss_flag=1;
		  //printf("%s\n",m[i].ss);
		    }
		  if(strcmp("METHOD",desc_flag)==0)
		    {
		      strcpy(m[i].method,temp);
		      //printf("%s\n",m[i].method);
		    }
		  if(strcmp("SCORE",desc_flag)==0)
		    {
		      m[i].score=atof(temp);
		      //strcpy(m[i].method,temp);
		      //printf("%lf\n",m[i].score);
		    }
		}
	      if(strcmp("MODEL",line_flag)==0)
		{
		  sscanf(buff,"%s %s",line_flag,temp);
		  m[i].rank=atoi(temp);
		}
	      
	      if (strcmp("ATOM",line_flag)==0)	/* Is it an ATOM entry? */
		{

		  strncpy_NULL(name,&buff[13],3);
		  if(strcmp("CA ",name) == 0)
		    {
		      //printf("%s",&buff[6]);
		      strncpy_NULL(temp_number,&buff[6],5);
		      //strncpy_NULL(name,&buff[13],3);
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
		      
		      //printf("test: %d %s %d %s %s %s %s %d %s %lf %lf %lf\n",i,line_flag,number,name,alt_loc,residue,chain,resnum,resname,x,y,z);
		      //if (strcmp("N",name)==0)	   /*  Is it an N atom => new residue? */
		      //printf("%s %s\n",old_resname,resname);
		      if(strcmp(old_resname,resname)!=0)
			{
			  m[i].sequence[residues]=aa321(residue);
			  residues++;
			  strcpy(alt_loc_check,alt_loc);
			  //printf("%s %s\n",resname,residue);
			}
		      if(strcmp(alt_loc_check,alt_loc)==0)
			{
		      //sscanf(&buff[22],"%d %lf %lf %lf",&resnum,&x,&y,&z);
			  // printf("test %d %s %s %d %lf %lf %lf\n",number,name,residue,resnum,x,y,z);
			  m[i].atm[atoms].x=x;
			  m[i].atm[atoms].y=y;
			  m[i].atm[atoms].z=z;
			  m[i].atm[atoms].resnum=resnum;
			  strcpy(m[i].atm[atoms].resname,resname);
			  m[i].atm[atoms].number=number; //atoi(number);
			  m[i].atm[atoms].rescount=residues;
			  m[i].atm[atoms].selected=TRUE;
			  strcpy(m[i].atm[atoms].name,name);
			  strcpy(m[i].atm[atoms].residue,residue);
			  //if(strcmp(old_resname,resname)!=0)
		      //	{
		      //	  m[i].sequence[residues-1]=aa321(residue);
		      //	}
			  
			  m[i].xcen+=x;
			  m[i].ycen+=y;
			  m[i].zcen+=z;
			  atoms++;
			  strcpy(old_resname,resname);
			}
		      //}
		    }
		}
	    }
	  m[i].sequence[residues]='\0';
	  m[i].atoms=atoms;
	  m[i].residues=residues;
	  m[i].xcen=m[i].xcen/atoms;
	  m[i].ycen=m[i].ycen/atoms;
	  m[i].zcen=m[i].zcen/atoms;
	  if(ss_flag == 0)
	    {
	      //fprintf(stderr,"No secondary structure information in file!\n");
	    }
	  
	  fclose(fp);
	  //#ifdef NOZLIB
	  //  fclose(fp);
	  //#else
	  //	  gzclose(fp);
	  //#endif
	  //	  if (atoms!=m[i].atoms)		/* Are file sizes indentical? */
	  //{
	  //  printf("Inconsistent number of atoms in file %s\n",m[i].filename);
	  //  return(1);
	  //}
	  
	}
      else
	{
	  printf("Couldn't open file %s\n",m[i].filename);
	  exit(1);
	}
    }
  return(0);
}

int read_molecules(temp)	/* Reads in molecules to be superimposed */
    int temp;
{
  int	i,j,k,atoms;	/* Counter variables */
  char	buff[512];	/* Input string */
  char	junk[512];	/* Throwaway string */
  char	line_flag[512];	/* PDB file line mode */
  char	residue[8];	/* PDB atom info for output */
  char	name[8];
  int	resnum;
  char  tmpname[2];
  int	owner;
  int	number;
  double	x,y,z;		/* Temporary coordinates values */
  FILE *fp;
    
  for (i=0;i<2;i++)
    {
      fp=fopen(m[i].filename,"r");	/* Does file exist? */
      if (fp!=NULL)	/* If yes, read in coordinates */
	{
	  /* Initialize things */
	  m[i].xcen=m[i].ycen=m[i].zcen=0;
	  atoms=0;
	  while(fgets(buff,255,fp)!=NULL)
	    {
	      buff[26]=' ';
	      sscanf(buff,"%s %d %s %s",line_flag,&number,name,residue);
	      if (strcmp("ATOM",line_flag)==0)	/* Is it an ATOM entry? */
		{
		  tmpname[0]=name[0];
		  tmpname[1]=0;
		  if (strcmp("H",tmpname)!=0)	   /*  Is it an H atom ? */
		    {
		      sscanf(&buff[22],"%d %lf %lf %lf",&resnum,&x,&y,&z);
		      //printf("read1: %s %d %lf %lf %lf\n",m[i].filename,resnum,x,y,z);
		      m[i].atm[atoms].x=x;
		      m[i].atm[atoms].y=y;
		      m[i].atm[atoms].z=z;
		      m[i].atm[atoms].owner=resnum;
		      m[i].atm[atoms].resnum=resnum;
		      m[i].atm[atoms].number=number;
		      strcpy(m[i].atm[atoms].name,name);
		      strcpy(m[i].atm[atoms].residue,residue);
		      m[i].xcen+=x;
		      m[i].ycen+=y;
		      m[i].zcen+=z;
		      atoms++;
		    }
		}
	    }
	  m[i].atoms=atoms;
	  fclose(fp);
	  if (atoms!=m[0].atoms)		/* Are file sizes indentical? */
	    {
	      //printf("Inconsistent number of atoms in file %s %d %d\n",m[i].filename,atoms, m[0].atoms);
	      //return(1);
	    }
	  /* Now center molecule */
	  m[i].xcen/=(double)atoms;
	  m[i].ycen/=(double)atoms;
	  m[i].zcen/=(double)atoms;
	  for (j=0;j<atoms;j++)
	    {
	      m[i].atm[j].x-=m[i].xcen;
	      m[i].atm[j].y-=m[i].ycen;
	      m[i].atm[j].z-=m[i].zcen;
	    }
	}
      else
	{
	  printf("Couldn't open file %s\n",m[i].filename);
	  exit(1);
	}
    }
  return(0);
}
int read_molecules2(temp)	/* Reads in molecules to be superimposed */
    int temp;
{
  int	i,j,k,atoms;	/* Counter variables */
  char	buff[512];	/* Input string */
  char	junk[512];	/* Throwaway string */
  char	line_flag[512];	/* PDB file line mode */
  char	residue[8];	/* PDB atom info for output */
  char	name[8];
  int	owner;
  int	number;
  double	x,y,z;		/* Temporary coordinates values */
  FILE *fp;
    
  for (i=0;i<2;i++)
    {
      fp=fopen(m[i].filename,"r");	/* Does file exist? */
      if (fp!=NULL)	/* If yes, read in coordinates */
	{
	  /* Initialize things */
	  m[i].xcen=m[i].ycen=m[i].zcen=0;
	  atoms=0;
	  while(fgets(buff,255,fp)!=NULL)
	    {
	      buff[26]=' ';
	      sscanf(buff,"%s %d %s %s",line_flag,&number,name,residue);
	      if (strcmp("ATOM",line_flag)==0)	/* Is it an ATOM entry? */
		{
		  if (strcmp("CA",name)==0)	   /*  Is it an CA atom ? */
		    {
		      sscanf(&buff[22],"%d %lf %lf %lf",&owner,&x,&y,&z);
		      //printf("read2: %s %d %lf %lf %lf\n",m[i].filename,owner,x,y,z);
		      m[i].atm[atoms].x=x;
		      m[i].atm[atoms].y=y;
		      m[i].atm[atoms].z=z;
		      m[i].atm[atoms].owner=owner;
		      m[i].atm[atoms].resnum=owner;
		      m[i].atm[atoms].number=number;
		      strcpy(m[i].atm[atoms].name,name);
		      strcpy(m[i].atm[atoms].residue,residue);
		      m[i].xcen+=x;
		      m[i].ycen+=y;
		      m[i].zcen+=z;
		      atoms++;
		    }
		}
	    }
	  m[i].atoms=atoms;
	  fclose(fp);
	  if (atoms!=m[0].atoms)		/* Are file sizes indentical? */
	    {
	      //printf("Inconsistent number of atoms in file %s\n",m[i].filename);
	      //return(1);
	    }
	  /* Now center molecule */
	  m[i].xcen/=(double)atoms;
	  m[i].ycen/=(double)atoms;
	  m[i].zcen/=(double)atoms;
	  for (j=0;j<atoms;j++)
	    {
	      m[i].atm[j].x-=m[i].xcen;
	      m[i].atm[j].y-=m[i].ycen;
	      m[i].atm[j].z-=m[i].zcen;
	    }
	}
      else
	{
	  printf("Couldn't open file %s\n",m[i].filename);
	  exit(1);
	}
    }
  return(0);
}
multiply_matrix(a,b,c)  /* computes C=AB */
    double a[3][3],b[3][3],c[3][3];
{
    int i,j,k;
    for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	    {
		c[i][j]=0.0;
		for (k=0;k<3;k++)
		    c[i][j]+=a[i][k]*b[k][j];
	    }
}
copy_matrix(f,t) /* copy matrix f into matrix t */
    double f[3][3],t[3][3];
{
    int i,j;
    for (i=0;i<3;i++)
	for (j=0;j<3;j++)
	    t[i][j]=f[i][j];
}
transpose_matrix(m) /* Transpose a 3x3 matrix */
    double m[3][3];
{
    double dummy;
    dummy=m[0][1]; m[0][1]=m[1][0]; m[1][0]=dummy;
    dummy=m[0][2]; m[0][2]=m[2][0]; m[2][0]=dummy;
    dummy=m[1][2]; m[1][2]=m[2][1]; m[2][1]=dummy;
}
double Sscore(m1,m2)	/* Find RMS superimposition of m1 on m2 */
    molecule 	*m1,*m2;		/* Molecules to be superimposed */
{
  int i;
  double 		S;		/* Final superimposition score */
  double		x,y,z;		/* Temporary coordinate variables */
  double                d0=5;
  double                rms2=0;
  for (i=0;i<m1->atoms;i++)
	{

	  x=m1->atm[i].x-m2->atm[i].x;
	  y=m1->atm[i].y-m2->atm[i].y;
	  z=m1->atm[i].z-m2->atm[i].z;
	  rms2=x*x+y*y+z*z;
	  S+=1/(1+rms2/d0);
	  m1->atm[i].rms=1/(1+rms2/d0);
	  //	  printf("RMSD %f\n",sqrt(m1->atm[i].rms));
	  // printf("1: %d %d %s %s %d %lf %lf %lf --- ",i,m1->atm[i].number,m1->atm[i].name,m1->atm[i].residue,m1->atm[i].resnum,m1->atm[i].x,m1->atm[i].y,m1->atm[i].z);
	  //printf("2: %d %d %s %s %d %lf %lf %lf\n",i,m2->atm[i].number,m2->atm[i].name,m2->atm[i].residue,m2->atm[i].resnum,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
	  //printf("%f\n",sqrt(m1->atm[i].rms));
	}
    return(S);

}
double superimpose_molecules(m1,m2,s)	/* Find RMS superimposition of m1 on m2 */
    molecule 	*m1,*m2;		/* Molecules to be superimposed */
    double 		s[3][3];	/* Final transformation matrix */
{
    int 		i,j,k;		/* Counter variables */
    double		u[3][3];	/* direct product matrix */
    double 		t[3][3];	/* Temporary storage matrix */
    double		ma[3][3];	/* x axis rotation matrix */
    double		mb[3][3];	/* y axis rotation matrix */
    double		mg[3][3];	/* z axis rotation matrix */
    double 		*d1,*d2;	/* usefule pointers */
    double 		error;		/* Final superimposition error */
    double		alpha=0.0; 	/* Angle of rotation around x axis */
    double		beta=0.0;	/* Angle of rotation around y axis */
    double		gamma=0.0;	/* Angle of rotation around z axis */
    double		x,y,z;		/* Temporary coordinate variables */
    for (i=0;i<3;i++)  /* Initialize matrices */
	for (j=0;j<3;j++)
	    s[i][j]=u[i][j]=0.0;
    s[0][0]=s[1][1]=s[2][2]=1.0;  /* Initialize S matrix to I */
    for (i=0;i<3;i++)  /* Initialize rotation matrices to I */
	for (j=0;j<3;j++)
	    ma[i][j]=mb[i][j]=mg[i][j]=s[i][j];
    for (i=0;i<m1->atoms;i++)  /* Construct U matrix */
	{
	    d1= &(m1->atm[i].x);
	    for (j=0;j<3;j++)
		{
		    d2= &(m2->atm[i].x);
		    for (k=0;k<3;k++)
			{
			    u[j][k]+=(*d1)*(*d2);
			    d2++;
			}
		    d1++;
		}
	}
    
    do
	{
	  
	    error=0.0;
	    /* Calculate x axis rotation */
	    alpha=atan((u[2][1]-u[1][2])/(u[1][1]+u[2][2]));
	    /* Insure we are heading for a minimum, not a maximum */
	    if (cos(alpha)*(u[1][1]+u[2][2])+sin(alpha)*(u[2][1]-u[1][2])<0.0)
		alpha+=PI;
	    ma[1][1]=ma[2][2]=cos(alpha);
	    ma[2][1]=sin(alpha); ma[1][2]= -ma[2][1];
	    transpose_matrix(ma);
	    multiply_matrix(u,ma,t);
	    transpose_matrix(ma);
	    copy_matrix(t,u);
	    multiply_matrix(ma,s,t);
	    copy_matrix(t,s);
	    /* Calculate y axis rotation */
	    beta=atan((u[0][2]-u[2][0])/(u[0][0]+u[2][2]));
	    /* Insure we are heading for a minimum, not a maximum */
	    if (cos(beta)*(u[0][0]+u[2][2])+sin(beta)*(u[0][2]-u[2][0])<0.0)
		beta+=PI;
	    mb[0][0]=mb[2][2]=cos(beta);
	    mb[0][2]=sin(beta); mb[2][0]= -mb[0][2];
	    transpose_matrix(mb); 
	    multiply_matrix(u,mb,t);
	    transpose_matrix(mb);
	    copy_matrix(t,u);
	    multiply_matrix(mb,s,t);
	    copy_matrix(t,s); 
	    /* Calculate z axis rotation */
	    gamma=atan((u[1][0]-u[0][1])/(u[0][0]+u[1][1]));
	    /* Insure we are heading for a minimum, not a maximum */
	    if (cos(gamma)*(u[0][0]+u[1][1])+sin(gamma)*(u[1][0]-u[0][1])<0.0)
		gamma+=PI;
	    mg[0][0]=mg[1][1]=cos(gamma);
	    mg[1][0]=sin(gamma); mg[0][1]= -mg[1][0];
	    transpose_matrix(mg);
	    multiply_matrix(u,mg,t);
	    transpose_matrix(mg);
	    copy_matrix(t,u);
	    multiply_matrix(mg,s,t);
	    copy_matrix(t,s);
	    error=fabs(alpha)+fabs(beta)+fabs(gamma); 
	    // printf("Error: %lf\n",error);
	}
    while (error>0.0001);	/* Is error low enough to stop? */
    /* Now calculate final RMS superimposition */
    error=0.0;
    for (i=0;i<m1->atoms;i++)
	{

	  x=s[0][0]*m2->atm[i].x+s[0][1]*m2->atm[i].y+s[0][2]*m2->atm[i].z;
	  y=s[1][0]*m2->atm[i].x+s[1][1]*m2->atm[i].y+s[1][2]*m2->atm[i].z;
	  z=s[2][0]*m2->atm[i].x+s[2][1]*m2->atm[i].y+s[2][2]*m2->atm[i].z;
	  m2->atm[i].x=x;
	  m2->atm[i].y=y;
	  m2->atm[i].z=z;
	  x=m1->atm[i].x-x;
	  y=m1->atm[i].y-y;
	  z=m1->atm[i].z-z;
	  m1->atm[i].rms=x*x+y*y+z*z;
	  error+=m1->atm[i].rms;
	  //	  printf("RMSD %f\n",sqrt(m1->atm[i].rms));
	  // printf("1: %d %d %s %s %d %lf %lf %lf --- ",i,m1->atm[i].number,m1->atm[i].name,m1->atm[i].residue,m1->atm[i].resnum,m1->atm[i].x,m1->atm[i].y,m1->atm[i].z);
	  //printf("2: %d %d %s %s %d %lf %lf %lf\n",i,m2->atm[i].number,m2->atm[i].name,m2->atm[i].residue,m2->atm[i].resnum,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
	  //printf("%f\n",sqrt(m1->atm[i].rms));
	}
    error/=(double)(m1->atoms);
    return(sqrt(error));
}
output_results()
{
    int	i,j;		/* Counter variable */
    double	sum;		/* Average RMS counter by residue */
    double	resatoms;	/* Number of atoms in each residue */
    char	name[8];	/* Temporary residue name variable */
    printf("Overall average RMS deviation is %lf Angstroms over %d atoms.\n",rms,m[0].atoms);
	    if (verbose_flag)
		{
		    /* First output RMS deviations by atom */
		    printf("\n");
		    printf("RMS deviations by atom\n");
		    printf("Atom \t Name \t  Owner \t RMS deviation\n");
		    for (i=0;i<m[0].atoms;i++)
			printf("%d \t %s \t %3s %4d \t %8.4lf\n",
			       m[0].atm[i].number,m[0].atm[i].name,
			       m[0].atm[i].residue,m[0].atm[i].resnum,
			       sqrt(m[0].atm[i].rms));
		    /* Now output average RMS deviations by residue */
		    printf("\n");
		    printf("Average RMS deviations by residue\n");
		    printf("Residue \t RMS deviation\n"); 
		    j= -1;
		    for (i=0;i<=m[0].atoms;i++)
			{
			    if (m[0].atm[i].owner!=j || i==m[0].atoms)
				{
				    if (j!=-1)
					printf("%3s %4d\t %8.4lf\n",name,j,sqrt(sum/resatoms));
				    sum=0;
				    resatoms=0;
				    j=m[0].atm[i].owner;
				    strcpy(name,m[0].atm[i].residue);		
				}
			    sum+=m[0].atm[i].rms;
			    resatoms++;
			}
		}
}
output_results2()
{
    int	i,j;		/* Counter variable */
    double	sum;		/* Average RMS counter by residue */
    double	resatoms;	/* Number of atoms in each residue */
    char	name[8];	/* Temporary residue name variable */
    printf("Overall average RMS deviation for CA is %lf Angstroms over %d atoms.\n",rms,m[0].atoms);
    if (verbose_flag)
	{
	    /* First output RMS deviations by atom */
	    printf("\n");
	    printf("RMS deviations by CA atom\n");
	    printf("Atom \t Name \t  Owner \t RMS deviation\n");
	    for (i=0;i<m[0].atoms;i++)
		printf("%d \t %s \t %3s %4d \t %8.4lf\n",
		       m[0].atm[i].number,m[0].atm[i].name,
		       m[0].atm[i].residue,m[0].atm[i].resnum,
		       sqrt(m[0].atm[i].rms));
	    
	}
}
output_dme(m1,m2)
    molecule 	*m1,*m2;	 /* Molecules to be superimposed */
{
    double	sum;		/* Average RMS counter by residue */
    double	resatoms;	/* Number of atoms in each residue */
    double      dist1,x1,y1,z1;
    double      dist2,x2,y2,z2;
    char	name[8];	/* Temporary residue name variable */
    int         i,j,count;
    
    sum=0.0;
    count=0;
    for (i=0;i<m1->atoms;i++)
	{
	    for (j=i+1;i<m1->atoms;i++)
		{
		    x1= ( (m1->atm[i].x)-(m1->atm[j].x))*( (m1->atm[i].x)-(m1->atm[j].x));
		    y1= ( (m1->atm[i].y)-(m1->atm[j].y))*( (m1->atm[i].y)-(m1->atm[j].y));
		    z1= ( (m1->atm[i].z)-(m1->atm[j].z))*( (m1->atm[i].z)-(m1->atm[j].z));
		    dist1=sqrt(x1+y1+z1);
		    
		    x2= ( (m2->atm[i].x)-(m2->atm[j].x))*( (m2->atm[i].x)-(m2->atm[j].x));
		    y2= ( (m2->atm[i].y)-(m2->atm[j].y))*( (m2->atm[i].y)-(m2->atm[j].y));
		    z2= ( (m2->atm[i].z)-(m2->atm[j].z))*( (m2->atm[i].z)-(m2->atm[j].z));
		    dist2=sqrt(x2+y2+z2);
		}
	    sum+=(dist1-dist2)*(dist1-dist2);
	    count+=1;
	}
    {
	sum/=(double)(count);
    }
    
    printf("Overall average DME deviation for CA is %lf Angstroms.\n",sqrt(sum));
}
output_file(filename)
    char filename[80];
{
    int i;
    FILE *fp;
    fp=fopen(filename,"w");
    if (fp!=NULL)
	{
	  /*  for (i=0;i<m[1].atoms;i++)
	      {
		fprintf(fp,"ATOM   %4d  %-4s%3s  %5s   %8.3lf%8.3lf%8.3lf  1.0 %5.2lf      RMSD\n",
			m[1].atm[i].number,m[1].atm[i].name,
			m[1].atm[i].residue,m[1].atm[i].resname,
			m[1].atm[i].x,m[1].atm[i].y,m[1].atm[i].z,
			sqrt(m[0].atm[i].rms));
	      
	      }
	      fprintf(fp,"TER \n");*/
/*	    printf("test %d  %d \n",m[1].atoms,m[0].atoms);*/
	    for (i=0;i<m[0].atoms;i++)
	      {
		fprintf(fp,"%3s%5s%8.6lf %8.6lf\n",
			m[0].atm[i].residue,m[0].atm[i].resname,
			m[0].atm[i].rms,5*sqrt(1/m[0].atm[i].rms-1));
	      }
	    //fprintf(fp,"END \n");
	}
    fclose(fp);
}

read_second_molecule(m1,m2,s)	/* Reads in molecules to be superimposed */
  molecule 	*m1,*m2;		/* Molecules to be superimposed */
  double	s[3][3];	/* Final transformation matrix */
{
    int	i,j,k,atoms;	/* Counter variables */
    char	buff[512];	/* Input string */
    char	junk[512];	/* Throwaway string */
    char	line_flag[512];	/* PDB file line mode */
    char	residue[8];	/* PDB atom info for output */
    char	name[8];
    double		u[3][3];	/* direct product matrix */
    double 		t[3][3];	/* Temporary storage matrix */
    double		ma[3][3];	/* x axis rotation matrix */
    double		mb[3][3];	/* y axis rotation matrix */
    double		mg[3][3];	/* z axis rotation matrix */
    double 		*d1,*d2;	/* usefule pointers */
    double 		error;		/* Final superimposition error */
    double		alpha=0.0; 	/* Angle of rotation around x axis */
    double		beta=0.0;	/* Angle of rotation around y axis */
    double		gamma=0.0;	/* Angle of rotation around z axis */
    int	owner;
    int	number;
    double	x,y,z;		/* Temporary coordinates values */
    FILE *fp;
    
    for (i=0;i<2;i++)
	{
	  fp=fopen(m[i].filename,"r");	/* Does file exist? */
	    if (fp!=NULL)	/* If yes, read in coordinates */
		{
		  /* Initialize things */
		  m[i].xcen=m[i].ycen=m[i].zcen=0;
		  atoms=0;
		  while(fgets(buff,255,fp)!=NULL)
		    {
		      sscanf(buff,"%s %d %s %s",line_flag,&number,name,residue);
		      if (strcmp("ATOM",line_flag)==0 && name[0] != 'H')	/* Is it an ATOM entry? */
			{
			  sscanf(&buff[22],"%d %lf %lf %lf",&owner,&x,&y,&z);
			  m[i].atm[atoms].x=x;
			  m[i].atm[atoms].y=y;
			  m[i].atm[atoms].z=z;
			  m[i].atm[atoms].owner=owner;
			  m[i].atm[atoms].number=number;
			  strcpy(m[i].atm[atoms].name,name);
			  strcpy(m[i].atm[atoms].residue,residue);
			  m[i].xcen+=x;
			  m[i].ycen+=y;
			  m[i].zcen+=z;
			  atoms++;
			}
		    }
		  m[i].atoms=atoms;
		  fclose(fp);
		  /* Now center molecule */
		  m[i].xcen/=(double)atoms;
		  m[i].ycen/=(double)atoms;
		  m[i].zcen/=(double)atoms;
		  for (j=0;j<atoms;j++)
		    {
		      m[i].atm[j].x-=m[i].xcen;
		      m[i].atm[j].y-=m[i].ycen;
		      m[i].atm[j].z-=m[i].zcen;
		    }
		}
	    else
		{
		  printf("Couldn't open file %s\n",m[i].filename);
		  exit(1);
		}
	}
    //check_molecules(&m[0],&m[1]);
    for (i=0;i<m1->atoms;i++)
	{

	    x=s[0][0]*m2->atm[i].x+s[0][1]*m2->atm[i].y+s[0][2]*m2->atm[i].z;
	    y=s[1][0]*m2->atm[i].x+s[1][1]*m2->atm[i].y+s[1][2]*m2->atm[i].z;
	    z=s[2][0]*m2->atm[i].x+s[2][1]*m2->atm[i].y+s[2][2]*m2->atm[i].z;
	    m2->atm[i].x=x;
	    m2->atm[i].y=y;
	    m2->atm[i].z=z;
	    x=m1->atm[i].x-x;
	    y=m1->atm[i].y-y;
	    z=m1->atm[i].z-z;
	    m1->atm[i].rms=x*x+y*y+z*z;
	}
    
    return(0);
}

check_molecules(m1,m2)	
    molecule 	*m1,*m2;		
{
  /* this will delete all atoms that are not in both molecules */
  /* Now we should extract the part of the molecules that exists in both */
  int          i,j,length1,length2;
  int          found;
  int          current_resnum;
  char         current_atom[4];
  int          number_removed1=0,number_removed2=0;

  length1=m1->atoms;
  length2=m2->atoms;
  if(length1==0)
    {
      fprintf(stderr,"Length of first protein is zero!\n");
      exit(1);
    }
  if(length2==0)
    {
      fprintf(stderr,"Length of second protein is zero!\n");
      exit(1);
    }
  //  printf("HERE\n");
  //if (minlength>m2->atoms) minlength=m2->atoms;
  i=0;
  //printf("%d\n",m1->atoms);
  while(i<m1->atoms)
    {
      //  printf("hej\n");
      if(!atom_exist(m1->atm[i].name,m1->atm[i].resname,m2))
	{
	  //printf("NO m1: %d %s %s %s\n",i,m1->atm[i].name,m1->atm[i].residue,m1->atm[i].resname);
	  delete_atom(m1,i);
	  number_removed1++;
	  //printf("NO m1: %s %s %d %d %d\n",m1->atm[i].name,m1->atm[i].residue,i,m1->atm[i].resnum,m2->atm[i].resnum);
	}
      else
	{
	  i++;
	  //printf("NO\n");
	  //delete_atom(m1,i);
	}
    }
  i=0;
  while(i<m2->atoms)
    {
      if(!atom_exist(m2->atm[i].name,m2->atm[i].resname,m1))
	{
	  delete_atom(m2,i);
	  number_removed2++;
	  //  printf("NO m2: %s %s %d %s\n",m2->atm[i].name,m2->atm[i].residue,i,m2->atm[i].resname);
	  //printf("YES m2: %s %s %d %d %d\n",m2->atm[i].name,m2->atm[i].residue,i,m2->atm[i].resnum,m1->atm[i].resnum);
	}
      else
	{
	  i++;
	}
    }
//   for(i=0;i<m2->atoms;i++)
//   {
//     printf("1: %d %d %s %s %s %lf %lf %lf --- ",i,m1->atm[i].number,m1->atm[i].name,m1->atm[i].residue,m1->atm[i].resname,m1->atm[i].x,m1->atm[i].y,m1->atm[i].z);
//     printf("2: %d %d %s %s %s %lf %lf %lf\n",i,m2->atm[i].number,m2->atm[i].name,m2->atm[i].residue,m2->atm[i].resname,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
//   }
  // printf("hej\n");
  /*The following will put all atom types in m2 the same order as in m1 */
  //printf("%d %d\n",m1->atoms,m2->atoms);
  for(i=0;i<m1->atoms;i++)
    {
      if(strcmp(m1->atm[i].name,m2->atm[i].name)!=0) //inconsistency do a move to correct
	{
	  j=i;
	  //printf("BEGIN\n");
	  //printf("%5s %5s %5d %5s %5s %5d\n",m1->atm[i].name,m1->atm[i].residue,m1->atm[i].resnum,m2->atm[i].name,m2->atm[i].residue,m2->atm[i].resnum);
	  //exit(1);
	  while(strcmp(m1->atm[i].name,m2->atm[j].name)!=0)
	    {
	      //printf("BEFORE %d %d %3s %3s %5d %3s %3s\n",i,j,m1->atm[i].name,m1->atm[i].residue,m1->atm[i].resnum,m2->atm[i].name,m2->atm[i].residue);
	      j++;
	    }
	 	  
	  //printf("BEFORE %d %d %3s %3s %3s %3s\n",i,j,m1->atm[i].name,m1->atm[i].residue,m2->atm[i].name,m2->atm[i].residue);
	  //printf("BEFORE: %3s %5d %5d %8.3lf%8.3lf%8.3lf\n",m2->atm[i].name,m2->atm[i].resnum,m2->atm[i].number,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
	  
	  move_atom(m2,i,j);
	  //printf("AFTER %d %d %3s %3s %3s %3s\n",i,j,m1->atm[i].name,m1->atm[i].residue,m2->atm[i].name,m2->atm[i].residue);
	  //printf("AFTER: %3s %5d %5d %8.3lf%8.3lf%8.3lf\n",m2->atm[i].name,m2->atm[i].resnum,m2->atm[i].number,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
	  //printf("%3s %5d %5d %8.3lf%8.3lf%8.3lf\n",m2->atm[i].name,m2->atm[i].resnum,m2->atm[i].number,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
	}
    }

  


//  while (i<=m1->atoms &&  i<=m2->atoms ){
//    //printf("%d %d %d\n",i,m1->atm[i].resnum,m2->atm[i].resnum);
//    printf("%d %d < %d\n",i,m1->atm[i].resnum,m2->atm[i].resnum);
//    while ( m1->atm[i].resnum < m2->atm[i].resnum && i< m1->atoms )
//      {
//	printf("Deleting m1: %s %s %d %d %d\n",m1->atm[i].name,m1->atm[i].residue,i,m1->atm[i].resnum,m2->atm[i].resnum);
//      delete_atom(m1,i);	
//      }
//    
//    while (m1->atm[i].resnum > m2->atm[i].resnum &&  i< m2->atoms )
//      {
//	printf("Deleting m2: %s %s %d %d %d\n",m2->atm[i].name,m2->atm[i].residue,i,m1->atm[i].resnum,m2->atm[i].resnum);
//	delete_atom(m2,i);	
//      }
//    i++;
//  }
//  while (m1->atoms < m2->atoms){    delete_atom(m2,m2->atoms)  ;}
//  while (m1->atoms > m2->atoms){    delete_atom(m1,m1->atoms)  ;}

//  printf("sizes: %d %d\n",m1->atoms,m2->atoms);
  ///for(i=0;i<m1->atoms;i++)
  //  {
       //    printf("1: %d %d %s %s %d %lf %lf %lf --- ",i,m1->atm[i].number,m1->atm[i].name,m1->atm[i].residue,m1->atm[i].resnum,m1->atm[i].x,m1->atm[i].y,m1->atm[i].z);
       //	 printf("2: %d %d %s %s %d %lf %lf %lf\n",i,m2->atm[i].number,m2->atm[i].name,m2->atm[i].residue,m2->atm[i].resnum,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
       //   }
//  for(i=0;i<m2->atoms;i++)
//    {
//      printf("2: %d %d %s %s %d %lf %lf %lf\n",i,m2->atm[i].number,m2->atm[i].name,m2->atm[i].residue,m2->atm[i].resnum,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
//     }
//  printf("FIRST\n");
//   for(i=0;i<m1->atoms;i++)
//    {
//      printf("%d %3s %3s %3s %3s %8.3lf%8.3lf%8.3lf\n",i,m1->atm[i].name,m1->atm[i].residue,m2->atm[i].name,m2->atm[i].residue,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
//
//    }
  
  if(number_removed1/length1>0.05 || number_removed2/length2>0.05)
    {
      fprintf(stderr,"Warning: More than 5% of the atoms removed 1 (%5.3f) 2 (%5.3)\n",number_removed1/length1,number_removed2/length2);
    }

  if (m1->atoms > 0){return(0);}
  else{
    printf("no identical  atoms in files %s\n",m[i].filename);
    exit(1);
  }
}

int atom_exist(char *name,char *resname,molecule *m)
{
  int i;
  for(i=0;i<m->atoms;i++)
    {
      if(strcmp(name,m->atm[i].name)==0 && strcmp(m->atm[i].resname,resname)==0)
	return TRUE;
    }
  return FALSE;
}


delete_atom(m1,num)	
    molecule 	*m1;		
    int         num;
{
  int          i,j,k;
  for (i=num;i<m1->atoms;i++)
    {
      j=i+1;
      strcpy(m1->atm[i].name,m1->atm[j].name);
      strcpy(m1->atm[i].residue,m1->atm[j].residue);
      strcpy(m1->atm[i].resname,m1->atm[j].resname);
      m1->atm[i].x=m1->atm[j].x;
      m1->atm[i].y=m1->atm[j].y;
      m1->atm[i].z=m1->atm[j].z;
      /*      m1->atm[i].rms=m1->atm[j].rms; */
      m1->atm[i].number=m1->atm[j].number;
      m1->atm[i].resnum=m1->atm[j].resnum;

    }
  m1->atoms--;
}


/* move atom from num1 to num2*/
move_atom(m1,num1,num2)	
    molecule 	*m1;		
    int         num1;
    int         num2;
{
  int          i,j,k;
  double x,y,z;
  int resnum,number,rescount;
  char name[8],residue[8],resname[8];
  
  //printf("Move from %d %d\n",num1,num2);
  //Save in temporary variables  
  //printf("%3s %5d %5d %8.3lf%8.3lf%8.3lf\n",m1->atm[num1].name,m1->atm[num1].resnum,m1->atm[num1].number,m1->atm[num1].x,m1->atm[num1].y,m1->atm[num1].z);
  //printf("%3s %3s\n",m1->atm[num1].name,m1->atm[num2].name);
  x=m1->atm[num2].x;
  y=m1->atm[num2].y;
  z=m1->atm[num2].z;
  resnum=m1->atm[num2].resnum;
  number=m1->atm[num2].number; //atoi(number);
  //residues=m[i].atm[num2].rescount=residues;
  //m[i].atm[num2].selected=TRUE;
  strcpy(resname,m1->atm[num2].resname);
  strcpy(name,m1->atm[num2].name);
  //printf("%3s %5d %5d %8.3lf%8.3lf%8.3lf\n",name,resnum,number,x,y,z);
  //Move num1 to num2
  m1->atm[num2].x=m1->atm[num1].x;
  m1->atm[num2].y=m1->atm[num1].y;
  m1->atm[num2].z=m1->atm[num1].z;
  m1->atm[num2].resnum=m1->atm[num1].resnum;
  m1->atm[num2].number=m1->atm[num1].number; 
  strcpy(m1->atm[num2].resname,m1->atm[num1].resname);
  strcpy(m1->atm[num2].name,m1->atm[num1].name);
  
  //Complete the move by copying num2 from temp variables.

  m1->atm[num1].x=x;
  m1->atm[num1].y=y;
  m1->atm[num1].z=z;
  m1->atm[num1].resnum=resnum;
  m1->atm[num1].number=number; 
  strcpy(m1->atm[num1].name,name);
  strcpy(m1->atm[num1].resname,resname);
  //printf("%3s %5d %5d %8.3lf%8.3lf%8.3lf\n",m1->atm[num1].name,m1->atm[num1].resnum,m1->atm[num1].number,m1->atm[num1].x,m1->atm[num1].y,m1->atm[num1].z);
  

  //if(strcmp(old_resname,resname)!=0)
  //	{
  //	  m[i].sequence[residues-1]=aa321(residue);
  //	}
  
}

void strncpy_NULL(char *dest, char *src, size_t n)
{
  strncpy(dest, src, n);
  dest[n]='\0';
}
char aa321(char* res)
{
  //char check_res[4];

  //check_res[0]=res[0];
  //check_res[1]=res[1];
  //check_res[2]=res[2];
  //res[3]='\0';
  //printf("%s %s\n",res,res);
  
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
  return 'X';
}
int center_molecule(molecule *m)	/* Reads in molecules to be superimposed */
     //    molecule *m;
{
  int	i;	/* Counter variables */
  int   natoms;  // Number of selected atoms
  double	xcen,ycen,zcen;		/* Temporary coordinates values */
  natoms=0;
  xcen=ycen=zcen=0;
  for (i=0;i<m->atoms;i++){
    //if (m->atm[i].selected){
      xcen+=m->atm[i].x;
      ycen+=m->atm[i].y;
      zcen+=m->atm[i].z;
      natoms++;
      //}
  }
  /* Now center molecule */
  xcen/=(double)natoms;
  ycen/=(double)natoms;
  zcen/=(double)natoms;
  //printf("TEST natoms: %d %f %f %f \n",natoms,xcen,ycen,zcen);
  for (i=0;i<m->atoms;i++)
    {
      m->atm[i].x-=xcen;
      m->atm[i].y-=ycen;
      m->atm[i].z-=zcen;
    }
  return(0);
}
