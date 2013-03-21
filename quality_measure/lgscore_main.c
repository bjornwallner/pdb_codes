/*
  Program LGSCORE   written by Scott M. Le Grand
  Even more modified and completely changed in 1999 by Arne Elofsson 
  Modified by Arne Elofsson in 1994
  Copyright 1992, Scott M. Le Grand and Arne Elofsson (1999)
  The program is available under the Gnu public license (see www.fsf.org)

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; version 2 of the License. With the
 exception that if you use this program in any scientific work you have
 to explicitly state that you have used LGSCORE and cite the relevant
 publication (dependent on what you have used LGSCORE for). My
 publicationlist can be found at http://www.sbc.su.se/~arne/papers/.
 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details. You should have received a
 copy of the GNU General Public License along with this program (in 
 the file gpl.txt); if not, write to the Free Software Foundation, 
 Inc., 59 Temple Place - Suite 330, Boston, MA 02111-13
  */

#include <stdio.h>
#include <math.h> 
#include <stdlib.h>

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

#define	MAXATMS		5000		/* Maximum allowable atoms */
#define PI		3.14159265	/* Useful constant */

#define TRUE		1		/* Boolean definitions */
#define FALSE		0

typedef struct {
  struct {
    double x,y,z;		/* Atomic coordinates */
    double rms;		/* RMS deviation */
    char residue[8];	/* PDB info for output */
    char name[8];
    int number;
    int resnum;
    int selected;
  } atm[MAXATMS];
  double xcen,ycen,zcen;
  int	atoms;			/* Current # of atoms */
  char	filename[800];		/* filename to read molecule from */
} molecule;

molecule	m[2];		/* Molecules to be compared */
int		out_flag;		
int		b_flag;		
int		align_flag;		
int		select1_flag;		
int		select2_flag;		
int		verbose_flag;		/* Verbose output flag */
int             fraction_flag;            /* Output flag */
int             short_flag;
int             arglist;
char            outfile[800];            /* Output name */
double		rms;			/* final RMS error */
int 		ma[2];
/*int             temp;*/
double		LG_pvalue(int N,double MS);
double		Levitt_Gerstein(molecule *m1,molecule *m2);
double		superimpose_molecules();
double          LG_pvalueF(int N,double score);//added by Fang
int read_molecules_ca();	/* Reads in molecules to be superimposed */
int read_molecule(molecule *m);	/* Reads in molecules to be superimposed */
int center_molecule(molecule *m);	
int copy_matrix(); /* copy matrix f into matrix t */
int transpose_matrix(); /* Transpose a 3x3 matrix */
int delete_atom(molecule *m1,int num);	
int check_molecules(molecule *m1,molecule *m2);	
int output_results();
int output_file(char *filename);
void strncpy_NULL(char *dest, char *src, size_t n);
void copymolecule(molecule *m1,molecule *m2);


main(argc,argv)		/* Main routine */
     int argc;
     char *argv[];
{
  //molecule      m0,m1;
  int	files=0;	/* Number of files parsed */
  int	i,j,k,a,b;	/* Counter variable */
  int	L=4;	// Minimum length of fragment for FindMaxSub matching
  int    minatoms=19; /* Variable for min no of atoms for P-value*/
  double	s[3][3];	/* Final transformation matrix */
  double          maxrmsd,minrmsd;
  double          fraction=1;
  double          numfrac=0;
  double          maxpvalue=1,maxpvalue1=1.;
  double          score,maxrms,maxrms1,maxscore,maxscore1;
  double          pvalue;
  int             maxatoms,maxpos,minpos,atoms;
  int             length=0;
  int             maxsub_15=0;
  int             maxsub_20=0;
  int             maxsub_25=0;
  int             maxsub_30=0;
  int             maxsub_35=0;
  int             maxselect1_0[MAXATMS],maxselect1_1[MAXATMS];
  int             maxselect2_0[MAXATMS],maxselect2_1[MAXATMS];
  double          bestsizescore[MAXATMS];
  int             maxatoms1;
  int             all_read=1;
  double          cutoffrmsd=0.;
  double          cutoffrms=3.;
  double          cutoff;
  double ss;//Added by Fang
  int             count;
  int             maxcount=10;
  char            filename[400][800];
  char            target[400][10];
  char            template[400][10];
  char            method[400][20];
  int             rank[400];
  char            filelist[800];
  char            bb[3];
  int             len_filelist;
  int             dotcount=0;
  int             lastchar=' ';
  FILE            *fp;  

  /* Initialize */
  files=0;
  verbose_flag=FALSE;
  out_flag=FALSE;
  select2_flag=FALSE;
  select1_flag=FALSE;
  short_flag=FALSE;
  arglist=FALSE;
      
  /* Parse command line for PDB files and options */
  i=1;
  /*    temp=1; */
  while (i<argc)
    {
      if (strcmp(argv[i],"-v")==0)
	{
	  verbose_flag=TRUE;
	  select2_flag=TRUE;
	  select1_flag=TRUE;
	}
      else if (strcmp(argv[i],"-f")==0)
	{
	  numfrac=0.2;
	  fraction_flag=TRUE;
	}
      else if (strcmp(argv[i],"-1")==0)
	{
	  select1_flag=TRUE;
	}
      else if (strcmp(argv[i],"-2")==0)
	{
	  select2_flag=TRUE;
	}
      else if (strcmp(argv[i],"-o")==0)
	{
	  align_flag=TRUE;
	}
      else if (strcmp(argv[i],"-b")==0)
	{
	  b_flag=TRUE;
	}
      else if (strcmp(argv[i],"-pdb")==0)
	{
	  out_flag=TRUE;
	  i++;
	  strcpy(outfile,argv[i]);
	}
      else if (strcmp(argv[i],"-min")==0)
	{
	  i++;
	  minatoms=atoi(argv[i]);
	}
      else if (strcmp(argv[i],"-s")==0)
	{
	  short_flag=TRUE;
	}
      else if (strcmp(argv[i],"-list")==0)
	{
	  arglist=TRUE;
	  i++;
	  strcpy(filelist,argv[i]);
	}
      else if(!arglist)
	{
	  strcpy(m[files].filename,argv[i]);
	  //printf("%s\n",m[files].filename);
	  files++;
	  len_filelist=2;
	}
      i++;
    }
 
  if(arglist)
    {
      //printf("Reading files from: %s\n",filelist);
      i=0;
      fp=fopen(filelist,"r");	/* Does file exist? */
      if (fp!=NULL)	/* If yes, read in coordinates */
	{
	  while(fscanf(fp,"%s",filename[i])!=EOF)
	    {
	      for(j=0,k=0,dotcount=0;j<strlen(filename[i]);j++)
	      {
		if(filename[i][j] == '.')
		  {
		    //printf("%c\n",filename[i][j]);
		    dotcount++;
		  }
		else if(filename[i][j] == '/')
		  {
		    dotcount=0;
		  }
		else if(dotcount == 2)
		  {
		    method[i][k]=filename[i][j];
		    k++;
		  }
		//  j=sscanf(filename[i],"%*s.%*s.%s",&method[i]);
		  //printf("%d %s %d %s\n",j,filename[i],rank[i],method[i]);
		//lastchar=filename[i][j];
		
	      }
	      method[i][k]='\0';
	      // printf("%s %s\n",filename[i],method[i]);
	      i++;
	    }
	  files=i;
	}
      else
	{
	  printf("Cannot open %s\n",filelist);
	  exit(1);
	}
      //for(i=0;i<len_filelist;i++)
      //	{
      //	  printf("Reading: %s\n",m_list[i].filename);
      //	  read_molecule(&m_list[i]);
      //	}
    }
  else if(files!=2)
    {
      printf("Usage: lgscore [ -a -f -v -1 -2 ] file1 file2\n");
      exit(1);
    }
  for(a=0;a<files;a++)
    {
      for(b=a+1;b<files;b++)
	{
	  //printf("%d %d\n",a,b);
	  if(strcmp(method[a],method[b]))
	    {
	      maxpvalue=1;
	      maxpvalue1=1.;
	      if(arglist)
		{
		  strcpy(m[0].filename,filename[a]);
		  strcpy(m[1].filename,filename[b]);
		  printf("%-40s %-40s\t",filename[a],filename[b]);
		}
	      if(read_molecules_ca()==0)
		{
		  //printf("Hello\n");
		  //if(arglist)
		  //	{
		  //	  molecule m[2];
		  ///	  copymolecule(&m[0],&m_list[a]);
		  //  copymolecule(&m[1],&m_list[b]);
		  //	  printf("%-30s %-30s\t",m[0].filename,m[1].filename);
		  //	}
		  length=max(m[0].atoms,m[1].atoms);
		  check_molecules(&m[0],&m[1]); // this deletes some atoms
		  //printf("testing %d %d\n",m[0].atoms,m[1].atoms);
		  atoms=m[0].atoms;
		  
		  for (i=0;i<m[0].atoms;i++){bestsizescore[i]=0;}
		  for (i=0;i<m[0].atoms-L;i=i+1){
		    //printf("%d %d",i,L);
		    for (j=0;j<m[0].atoms;j++){ // deselect all atoms
		      m[0].atm[j].selected=FALSE;
		      m[1].atm[j].selected=FALSE;
		    }
		    for (j=i;j<i+L;j++){ // select some
		      //printf("%d %d %d\n",j,L,i);
		      m[0].atm[j].selected=TRUE;
		      m[1].atm[j].selected=TRUE;
		    }
		    j=L;
		    center_molecule(&m[0]);
		    center_molecule(&m[1]);
		    //printf("calling superimpose 1\n");
		    rms=superimpose_molecules(&m[0],&m[1],s);	   
		    count=0;
		    while (j<m[0].atoms && rms<(double)(5*(j+225)/225)){
		      count++;
		      minrmsd=9999.0;
		      for (k=0;k<m[0].atoms;k++) {
			if(!m[0].atm[k].selected)
			  {
			    if(m[0].atm[k].rms < minrmsd)
			      {
				minpos=k;
				minrmsd=m[0].atm[k].rms;
			      }
			    else if(m[0].atm[k].rms < cutoffrms)
			      {
				//printf("Selecting %d %f\n",k,m[0].atm[k].rms);
				m[0].atm[k].selected=TRUE;
				m[1].atm[k].selected=TRUE;
				j++;
			      }
			  }
		      }
		      j++;
		      //printf("Selecting %d %f %d\n",minpos,minrmsd,j);
		      m[0].atm[minpos].selected=TRUE;
		      m[1].atm[minpos].selected=TRUE;
		      
		      //	     printf("testing2a %f\n",m[1].atm[i].x);
		      center_molecule(&m[0]);
		      center_molecule(&m[1]);
		      //printf("calling superimpose 2\t %d %d\n",i,j);
		      rms=superimpose_molecules(&m[0],&m[1],s);
		      score=Levitt_Gerstein(&m[0],&m[1]);
		      //pvalue=LG_pvalue(j,score);
		      pvalue=LG_pvalueF(j,score);
		      if (b_flag)
			{
			  if (score > bestsizescore[j]){bestsizescore[j]=score;}
			}
		      if (fraction_flag)
			{
			  if ((double)j/atoms <= fraction )
			    {
			      printf("FRAC1:\t%.2f\t%d\t%.1f\t%.1f\t%e\n",
				     (double)j/atoms,j,rms,score,
				     pvalue);
			      fraction-=numfrac;
			      //printf("TEST %f %f\n",fraction,numfrac);
			    } 
			}
		      if (verbose_flag) output_results();
		      if ((pvalue <= maxpvalue) && (j>minatoms))
			{
			  if (fraction_flag)
			    {
			      printf("GOOD1:\t%.2f\t%d\t%.1f\t%.1f\t%e\n",
				     (double)j/atoms,j,rms,score,
				     pvalue);
			    }
			  //printf ("TEST:\t%d\t%f\t%f\t%f\n",m[0].atoms,pvalue,rms,score);
			  maxpvalue=pvalue;
			  maxatoms=j;
			  maxrms=rms;
			  maxscore=score;
			  for (k=0;k<m[0].atoms;k++)
			    {
			      maxselect1_0[k]=m[0].atm[k].selected;
			    }
			  for (k=0;k<m[1].atoms;k++)
			    {
			      maxselect1_1[k]=m[1].atm[k].selected;
			    }
			}
		    }
		    //	     printf("testing5 %f\n",m[1].atm[i].x);
		  }
		  //printf("%d %d\n",count,j);
	      //	 printf("BEST1:\t%.2f\t%d\t%.1f\t%.1f\t%e\n",
		  // (double)maxatoms/(double)atoms,maxatoms,maxrms,maxscore,maxpvalue);
		  maxatoms1=maxatoms;
		  maxrms1=maxrms;
		  maxscore1=maxscore;
		  maxpvalue1=maxpvalue;
		  //(double)maxatoms/(double)m[0].atoms,maxatoms,maxrms,maxscore,maxpvalue,maxpvalue1);
		  if (maxatoms1<minatoms){maxatoms1=maxrms1=maxscore1=0;maxpvalue1=1;}
		  //           ss=(log10(maxpvalue)+0.6108)/ma1;
		  
		  //lg[a][b]=(log10(1/maxpvalue1));
		  //lg[b][a]=lg[a][b];
		  if(short_flag)
		    {
		      printf("%7.5lf\n",(log10(1/maxpvalue1)));
		    }
		  else
		    {
		      ss=-1000*(log10(maxpvalue1)+0.6108)/ma[0]; //Why ma[0]????? dependent on order
		      printf("SCORE:\t%d\t%d\t%.1f\t%.1f\t%e\t%e\t%7.5lf\n",
			     length,maxatoms1,maxrms1,maxscore1,maxpvalue1,ss,(log10(1/maxpvalue1)));
		    }
		  if (align_flag){
		    for (i=0;i<m[0].atoms;i++){m[0].atm[i].selected=maxselect1_0[i];}
		    for (i=0;i<m[1].atoms;i++){m[1].atm[i].selected=maxselect1_1[i];}
		    center_molecule(&m[0]);
		    center_molecule(&m[1]);
		    //printf("calling superimpose 3\n");
		    rms=superimpose_molecules(&m[0],&m[1],s);	   
		    output_results();
		  }
		  //MAXSUB:   \t1.5\t2.0\t2.5\t3.0\t3.5\n
		  if (b_flag){
		    for (i=0;i<m[0].atoms;i++){
		      //printf("SCORE-TEST\t%d\t%f\t%e\n",i,bestsizescore[i],LG_pvalue(i,bestsizescore[i]));
		      printf("SCORE-TEST\t%d\t%f\t%e\n",i,bestsizescore[i],LG_pvalueF(i,bestsizescore[i]));
		    }
		  }
		  if (out_flag){
		    for (k=0;k<m[0].atoms;k++){m[0].atm[k].selected=maxselect1_0[k];
		    //	    printf("TEST %d %d\n",k,m[0].atm[k].selected);
		    }
		    for (k=0;k<m[1].atoms;k++){m[1].atm[k].selected=maxselect1_1[k];}
		    center_molecule(&m[0]);
		    center_molecule(&m[1]);
		    rms=superimpose_molecules(&m[0],&m[1],s);	   
		    //printf("TEST-out: %s\n",outfile);
		    output_file(outfile);
		  }
		}
	      else
		{
		  printf("Unable to do what you told me... My fault... I'm stupied\n");
		}
	    }
	}
      
    }
}



int read_molecule(molecule *m)	/* Read in molecule store it in m */
{
  int	i,j,k,atoms;	/* Counter variables */
  char	buff[512];	/* Input string */
  char	junk[512];	/* Throwaway string */
  char	line_flag[512];	/* PDB file line mode */
  char	residue[8];	/* PDB atom info for output */
  char	name[8];
  int	resnum;
  int	number;
  double	x,y,z;		/* Temporary coordinates values */
  char  old_resname[9]="noname";
  char  resname[9];
  char x_temp[11];
  char y_temp[11];
  char z_temp[11];
  char chain[2];
  char inscode[2]=" ";
  char alt_loc[2]=" ";
  char temp_number[11];
  char temp_resnum[11];
  FILE *fp;
    
  for (i=0;i<2;i++)
    {
      //#ifdef  NOZLIB
      fp=fopen(m->filename,"r");	/* Does file exist? */
      //#else
      //fp=gzopen(m.filename,"r");	/* Does file exist? */
      //#endif
      if (fp!=NULL)	/* If yes, read in coordinates */
	{
	  /* Initialize things */
	  m->xcen=m->ycen=m->zcen=0;
	  atoms=0;
	  //#ifdef NOZLIB
	  while(fgets(buff,255,fp)!=NULL)
	    //#else
	    //while(gzgets(fp,buff,255)!=NULL)
	    //#endif
	    {
	      sscanf(buff,"%s",line_flag);
	      //sscanf(buff,"%s %d %s %s",line_flag,&number,name,residue);
	      if (strcmp("ATOM",line_flag)==0)	/* Is it an ATOM entry? */
		{
		  strncpy_NULL(name,&buff[13],3);
		  name[3]='\0';
		  if (strcmp("CA ",name)==0)	   /*  Is it an CA atom ? */
		    {
		      //sscanf(&buff[22],"%d %lf %lf %lf",&resnum,&x,&y,&z);
		      // printf("%s %s",m->filename,buff);
		      //		      printf("test %d\n",resnum);
		      //strncpy_NULL(temp_number,&buff[6],5);
		      //strncpy_NULL(name,&buff[13],3);
		      //strncpy_NULL(alt_loc,&buff[16],1);
		    
		      strncpy_NULL(residue,&buff[17],3);
		      strncpy(residue,&buff[17],3);
		      strncpy_NULL(chain,&buff[21],1);
		      strncpy_NULL(temp_resnum,&buff[22],4);
		      strncpy_NULL(resname,&buff[22],5);
		      strncpy_NULL(x_temp,&buff[30],8);
		      strncpy_NULL(y_temp,&buff[38],8);
		      strncpy_NULL(z_temp,&buff[46],8);
		      number=atoi(temp_number);
		      resnum=atoi(temp_resnum);
		      //x=atof(x_temp);
		      //y=atof(y_temp);
		      //z=atof(z_temp);
		      //printf("%s %s %s\n",x_temp,y_temp,z_temp);
		      //printf("%d %d %f %f %f\n",number,resnum,x,y,z);
		      
		      m->atm[atoms].x=atof(x_temp);
		      m->atm[atoms].y=atof(y_temp);
		      m->atm[atoms].z=atof(z_temp);
		      m->atm[atoms].resnum=resnum;
		      m->atm[atoms].number=number;
		      m->atm[atoms].selected=TRUE;
		      strcpy(m->atm[atoms].name,name);
		      strcpy(m->atm[atoms].residue,residue);
		      m->xcen+=x;
		      m->ycen+=y;
		      m->zcen+=z;
		      atoms++;
		    }
		}
	    }
	  //	  printf("test %d %d\n",i,atoms);
	  m->atoms=atoms;
	  //ma=atoms;
	  //#ifdef NOZLIB
	  fclose(fp);
	  //#else
	  //gzclose(fp);
	  //#endif
	  //	  if (atoms!=m[0]->atoms)		/* Are file sizes indentical? */
	  //{
	  //  printf("Inconsistent number of atoms in file %s\n",m->filename);
	  //  return(1);
	  //}
	
	}
      else
	{
	  printf("Couldn't open file %s\n",m->filename);
	  exit(1);
	}
    }
  return(0);
}

int read_molecules_ca()	/* Reads in molecules to be superimposed */
{
  int	i,j,k,atoms;	/* Counter variables */
  char	buff[512];	/* Input string */
  char	junk[512];	/* Throwaway string */
  char	line_flag[512];	/* PDB file line mode */
  char	residue[8];	/* PDB atom info for output */
  char	name[8];
  int	resnum;
  int	number;
  double	x,y,z;		/* Temporary coordinates values */
  char  old_resname[9]="noname";
  char  resname[9];
  char x_temp[11];
  char y_temp[11];
  char z_temp[11];
  char chain[2];
  char inscode[2]=" ";
  char alt_loc[2]=" ";
  char temp_number[11];
  char temp_resnum[11];
  FILE *fp;
    
  for (i=0;i<2;i++)
    {
      //#ifdef  NOZLIB
      fp=fopen(m[i].filename,"r");	/* Does file exist? */
      //#else
      //fp=gzopen(m[i].filename,"r");	/* Does file exist? */
      //#endif
      if (fp!=NULL)	/* If yes, read in coordinates */
	{
	  /* Initialize things */
	  m[i].xcen=m[i].ycen=m[i].zcen=0;
	  atoms=0;
	  //#ifdef NOZLIB
	  while(fgets(buff,255,fp)!=NULL)
	    //#else
	    //while(gzgets(fp,buff,255)!=NULL)
	    //#endif
	    {
	      sscanf(buff,"%s",line_flag);
	      //sscanf(buff,"%s %d %s %s",line_flag,&number,name,residue);
	      if (strcmp("ATOM",line_flag)==0)	/* Is it an ATOM entry? */
		{
		  strncpy_NULL(name,&buff[13],3);
		  name[3]='\0';
		  if (strcmp("CA ",name)==0)	   /*  Is it an CA atom ? */
		    {
		      //sscanf(&buff[22],"%d %lf %lf %lf",&resnum,&x,&y,&z);
		      // printf("%s",buff);
		      //		      printf("test %d\n",resnum);
		      //strncpy_NULL(temp_number,&buff[6],5);
		      //strncpy_NULL(name,&buff[13],3);
		      //strncpy_NULL(alt_loc,&buff[16],1);
		    
		      strncpy_NULL(residue,&buff[17],3);
		      strncpy(residue,&buff[17],3);
		      strncpy_NULL(chain,&buff[21],1);
		      strncpy_NULL(temp_resnum,&buff[22],4);
		      strncpy_NULL(resname,&buff[22],5);
		      strncpy_NULL(x_temp,&buff[30],8);
		      strncpy_NULL(y_temp,&buff[38],8);
		      strncpy_NULL(z_temp,&buff[46],8);
		      number=atoi(temp_number);
		      resnum=atoi(temp_resnum);
		      //x=atof(x_temp);
		      //y=atof(y_temp);
		      //z=atof(z_temp);
		      //printf("%s %s %s\n",x_temp,y_temp,z_temp);
		      //printf("%d %d %f %f %f\n",number,resnum,x,y,z);
		      
		      m[i].atm[atoms].x=atof(x_temp);
		      m[i].atm[atoms].y=atof(y_temp);
		      m[i].atm[atoms].z=atof(z_temp);
		      m[i].atm[atoms].resnum=resnum;
		      m[i].atm[atoms].number=number;
		      m[i].atm[atoms].selected=TRUE;
		      strcpy(m[i].atm[atoms].name,name);
		      strcpy(m[i].atm[atoms].residue,residue);
		      m[i].xcen+=x;
		      m[i].ycen+=y;
		      m[i].zcen+=z;
		      atoms++;
		    }
		}
	    }
	  //	  printf("test %d %d\n",i,atoms);
	  m[i].atoms=atoms;
	  ma[i]=atoms;
	  //#ifdef NOZLIB
	  fclose(fp);
	  //#else
	  //gzclose(fp);
	  //#endif
	  //	  if (atoms!=m[0].atoms)		/* Are file sizes indentical? */
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

int center_molecule(molecule *m)	/* Reads in molecules to be superimposed */
     //    molecule *m;
{
  int	i;	/* Counter variables */
  int   natoms;  // Number of selected atoms
  double	xcen,ycen,zcen;		/* Temporary coordinates values */
  natoms=0;
  xcen=ycen=zcen=0;
  for (i=0;i<m->atoms;i++){
    if (m->atm[i].selected){
      xcen+=m->atm[i].x;
      ycen+=m->atm[i].y;
      zcen+=m->atm[i].z;
      natoms++;
    }
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
int multiply_matrix(a,b,c)  /* computes C=AB */
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
int copy_matrix(f,t) /* copy matrix f into matrix t */
     double f[3][3],t[3][3];
{
  int i,j;
  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      t[i][j]=f[i][j];
}
int transpose_matrix(m) /* Transpose a 3x3 matrix */
     double m[3][3];
{
  double dummy;
  dummy=m[0][1]; m[0][1]=m[1][0]; m[1][0]=dummy;
  dummy=m[0][2]; m[0][2]=m[2][0]; m[2][0]=dummy;
  dummy=m[1][2]; m[1][2]=m[2][1]; m[2][1]=dummy;
}
int delete_atom(molecule *m1,int num)	
{
  int          i,j,k;
  for (i=num;i<m1->atoms;i++)
    {
      j=i+1;
      strcpy(m1->atm[i].name,m1->atm[j].name);
      strcpy(m1->atm[i].residue,m1->atm[j].residue);
      m1->atm[i].x=m1->atm[j].x;
      m1->atm[i].y=m1->atm[j].y;
      m1->atm[i].z=m1->atm[j].z;
      /*      m1->atm[i].rms=m1->atm[j].rms; */
      m1->atm[i].number=m1->atm[j].number;
      m1->atm[i].resnum=m1->atm[j].resnum;

    }
  m1->atoms--;
}

int check_molecules(molecule *m1,molecule *m2)	
{
  /* this will delete all atoms that are not in both molecules */
  /* Now we should extract the part of the molecules that exists in both */
  int          i,j,minlength;
  int          found;

  minlength=m1->atoms;
  if (minlength>m1->atoms) minlength=m2->atoms;
  i=0;
  while (i<=m1->atoms &&  i<=m2->atoms ){
    while ( m1->atm[i].resnum < m2->atm[i].resnum && i< m1->atoms ){
      //printf("Deleting m1: %d %d %d\n",i,m1->atm[i].resnum,m2->atm[i].resnum);
      delete_atom(m1,i);	
    }
    while (m1->atm[i].resnum > m2->atm[i].resnum &&  i< m2->atoms ){
      //	printf("Deleting m2: %d %d %d\n",i,m1->atm[i].resnum,m2->atm[i].resnum);
      delete_atom(m2,i);	
    }
    i++;
  }
  while (m1->atoms < m2->atoms){    delete_atom(m2,m2->atoms)  ;}
  while (m1->atoms > m2->atoms){    delete_atom(m1,m1->atoms)  ;}
}


double Levitt_Gerstein(molecule *m1,molecule *m2)	
{
  int          i,j,last1,last2;
  int          numgap=0;
  double       sum=0.;
  //  double       d0=25.; // 5**2
  double       d0=5.;  // 2.24**2
  double       M=20.;
  //last1=m1->atm[0].resnum-1;
  //last2=m2->atm[0].resnum-1;
  j=0;
  for (i=0;i<m1->atoms;i++)
    {
      if (m1->atm[i].selected){
	//printf("%d\n",m1->atm[i].resnum);
	//last1++;
	//last2++;
	//if (m1->atm[i].resnum > last1 || m2->atm[i].resnum > last2 )
	//  {
	//    	  numgap++;
	//  }
	sum+=(1/(1+m1->atm[i].rms/d0));
	//	printf("Test>\t%d\t%f\t%f\n",numgap,sum,m1->atm[i].rms);
	//last1=m1->atm[i].resnum;
	//last2=m2->atm[i].resnum;
	//j++;
	//	printf("Test> %d  \t%d\t%d\t%f\n",i,m1->atm[i].selected,sum);
      }
    }
  //printf("TEST-LG>\t%d\t%d\t%f\t%f\n",j,numgap,sum,M*(sum-numgap/2));
  return(M*(sum-numgap/2));
}

double LG_pvalue(int N,double MS)
{
  double min_SD=1.e-20;
  double ln,Mean_MS,SD_MS,Z,expZ;
  double ln60=log(120.);
  double b=2.0*ln60*18.411424-4.501719;
  double a=(ln60*ln60)*18.411424+ln60*(-4.501719)+2.637163-b*ln60;
  if (N<6) {
    expZ=1.;
    return(expZ);
  }
  ln=log((double)N);
  if (N<120){
    Mean_MS=ln*ln*18.411424+ln*(-4.501719)+ 2.637163;
    SD_MS=ln*21.351857 -37.521269;
  }else{
    Mean_MS=ln*b+a;
    SD_MS=ln60*21.351857 -37.521269;
  }
  //  if (SD_MS<min_SD) {SD_MS=min_SD;}
  Z=(MS-Mean_MS)/SD_MS;

  expZ=exp((double)-1.*Z);
  //printf("TEST-P> %d\t%f\t%f\t%e\t%e\t%e\t%e\n",N,MS,Mean_MS,SD_MS,Z,expZ,1-exp(-1*expZ));
  if (Z>20){
    return(expZ);
    //  }else if (Z<-100){
    //expZ=1.;
    //return(expZ);
  }else{
    return(1-exp(-expZ));
  }
}

double superimpose_molecules(m1,m2,s)	/* Find RMS superimposition of m1 on m2 */
     molecule 	*m1,*m2;		/* Molecules to be superimposed */
     double 		s[3][3];	/* Final transformation matrix */
{
  int 		i,j,k;		/* Counter variables */
  int                 natoms=0;
  double		u[3][3];	/* direct product matrix */
  double 		t[3][3];	/* Temporary storage matrix */
  double		ma[3][3];	/* x axis rotation matrix */
  double		mb[3][3];	/* y axis rotation matrix */
  double		mg[3][3];	/* z axis rotation matrix */
  double 		*d1,*d2;	/* usefule pointers */
  double 		error;		/* Final superimposition error */
  double 		error2;		/* Final superimposition error */
  double		alpha=0.0; 	/* Angle of rotation around x axis */
  double		beta=0.0;	/* Angle of rotation around y axis */
  double		gamma=0.0;	/* Angle of rotation around z axis */
  double		dist,dist2,x,y,z;		/* Temporary coordinate variables */
  for (i=0;i<3;i++)  /* Initialize matrices */
    for (j=0;j<3;j++)
      s[i][j]=u[i][j]=0.0;
  s[0][0]=s[1][1]=s[2][2]=1.0;  /* Initialize S matrix to I */
  for (i=0;i<3;i++)  /* Initialize rotation matrices to I */
    for (j=0;j<3;j++)
      ma[i][j]=mb[i][j]=mg[i][j]=s[i][j];
  for (i=0;i<m1->atoms;i++)  /* Construct U matrix */
    {
      if (m1->atm[i].selected){
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
    }
  while (error>0.0001);	/* Is error low enough to stop? */
  /* Now calculate final RMS superimposition */
  error=0.0;
  error2=0.0;
  natoms=0;
  /*    printf("testing7b %f %f\n",m2->atm[1].x,error);
	printf ("s1:\t%f\t%f\t%f\n",s[0][0],s[0][1],s[0][2]);
	printf ("s2:\t%f\t%f\t%f\n",s[1][0],s[1][1],s[1][2]);
	printf ("s3:\t%f\t%f\t%f\n",s[2][0],s[2][1],s[2][2]);*/
  if (! (s[0][0]>0 || s[0][0]<0))  {
    //      printf ("bar\n");
    for (i=0;i<3;i++)  /* Initialize matrices */
      for (j=0;j<3;j++)
	s[i][j]=u[i][j]=0.0;
    s[0][0]=s[1][1]=s[2][2]=1.0;  /* Initialize S matrix to I */
  }
  for (i=0;i<m2->atoms;i++)
    {

      //dist=m2->atm[i].x*m2->atm[i].x+m2->atm[i].y*m2->atm[i].y+m2->atm[i].z*m2->atm[i].z;
      //printf("TEST1-m2\t %d\t%6.3f %6.3f %6.3f\t%e\n",i,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z,dist);
      x=s[0][0]*m2->atm[i].x+s[0][1]*m2->atm[i].y+s[0][2]*m2->atm[i].z;
      y=s[1][0]*m2->atm[i].x+s[1][1]*m2->atm[i].y+s[1][2]*m2->atm[i].z;
      z=s[2][0]*m2->atm[i].x+s[2][1]*m2->atm[i].y+s[2][2]*m2->atm[i].z;
      m2->atm[i].x=x;
      m2->atm[i].y=y;
      m2->atm[i].z=z;
      //dist2=m2->atm[i].x*m2->atm[i].x+m2->atm[i].y*m2->atm[i].y+m2->atm[i].z*m2->atm[i].z;
      //printf("TEST2-m2\t %d\t%6.3f %6.3f %6.3f\t%e\t%e\n",i,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z,dist2,dist2-dist);
      x=m1->atm[i].x-x;
      y=m1->atm[i].y-y;
      z=m1->atm[i].z-z;
      m2->atm[i].rms=m1->atm[i].rms=x*x+y*y+z*z;
      error+=m1->atm[i].rms;
      //printf("TEST:\t %d\t%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f \n",i,m1->atm[i].x,m2->atm[i].x,m1->atm[i].y,m2->atm[i].y,m1->atm[i].z,m2->atm[i].z);
      //printf("TEST1:\t %d\t%f %f %f\n",i,x*x,y*y,z*z);
      //printf("TEST2:\t %d %d\t%f\n",i,m1->atm[i].selected,m1->atm[i].rms);
      if (m1->atm[i].selected){
	error2+=m1->atm[i].rms;
	natoms++;
      }
      /*printf ("TEST1: %f %f %f\n",m1->atm[i].x,m1->atm[i].y,m1->atm[i].z);
	printf ("TEST2: %f %f %f\n",m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
	printf ("TEST3: %d %f %f\n",i,error,m1->atm[i].rms); */
    }
  //    printf("testing8 %f\n",m2->atm[1].x);
  error/=(double)(m1->atoms);
  error2/=(double)(natoms);
  return(sqrt(error2));
}
int output_results()
{
  int	i,j;		/* Counter variable */
  double	sum;		/* Average RMS counter by residue */
  double	resatoms;	/* Number of atoms in each residue */
  char	name[8];	/* Temporary residue name variable */
  printf("Overall average RMS deviation is %lf Angstroms.\n",rms);
  /* First output RMS deviations by atom */
  printf("\n");
  printf("RMS deviations by atom\n");
  printf("Atom \t Name \t  residue \t RMS deviation\tAtom \t Name \t  residue \t RMS deviation\n");
  for (i=0;i<m[0].atoms;i++){
    /*	  printf ("Test1: %f %f %f\n",m[0].atm[i].x,m[0].atm[i].y,m[0].atm[i].z);
	  printf ("Test2: %f %f %f\n",m[1].atm[i].x,m[1].atm[i].y,m[1].atm[i].z); */
    if ( m[0].atm[i].selected){
      printf("%d \t %s \t %3s %4d \t %8.4lf\t%d \t %s \t %3s %4d \t %8.4lf\n",
	     m[0].atm[i].number,m[0].atm[i].name,
	     m[0].atm[i].residue,m[0].atm[i].resnum,
	     sqrt(m[0].atm[i].rms),
	     m[1].atm[i].number,m[1].atm[i].name,
	     m[1].atm[i].residue,m[1].atm[i].resnum,
	     sqrt(m[1].atm[i].rms));
    }
  }
}
int output_file(char *filename)
{
  int i;
  FILE *fp;
  fp=fopen(filename,"w");
  if (fp!=NULL)
    {
      for (i=0;i<m[1].atoms;i++)
	{
	  fprintf(fp,"ATOM   %-4d  %-4s%-3s  %-4d     %7.3lf %7.3lf %7.3lf  1.0   %-1d.00     \n",
		  m[1].atm[i].number,m[1].atm[i].name,
		  m[1].atm[i].residue,m[1].atm[i].resnum,
		  m[1].atm[i].x,m[1].atm[i].y,m[1].atm[i].z,
		  m[0].atm[i].selected);
	}
      fprintf(fp,"END \n");
      /*	    printf("test %d  %d \n",m[1].atoms,m[0].atoms);*/
      for (i=0;i<m[0].atoms;i++)
	{
	  fprintf(fp,"ATOM   %-4d  %-4s%-3s  %-4d     %7.3lf %7.3lf %7.3lf  1.0   %-1d.50     \n",
		  m[0].atm[i].number,m[0].atm[i].name,
		  m[0].atm[i].residue,m[0].atm[i].resnum,
		  m[0].atm[i].x,m[0].atm[i].y,m[0].atm[i].z,
		  m[0].atm[i].selected);
	}
      fprintf(fp,"END \n");
    }
  fclose(fp);
}


////////////////////////
//The following is added by Fang
double LG_pvalueF(N,score)
     int N;
     double score;
{
  double Z,s1,mean,SD,sc,expZ;
  s1=(double)N;
  mean=pow(s1,0.3264)/(0.0437*pow(s1,0.0003)+0.0790);
  SD=s1/(0.0417*s1+3.3700);
  Z=(score/10-mean)/SD;
  expZ=exp((double)-1.*Z);
  if (Z>20)s1=expZ;
  else s1=1-exp(-expZ);
  return s1;
}


void strncpy_NULL(char *dest, char *src, size_t n)
{
  strncpy(dest, src, n);
  dest[n]='\0';
}

void copymolecule(molecule *m1,molecule *m2)
{
  /* Copy m2 to m1 */
  int i;
  for(i=0;i<m2->atoms;i++)
    {
      m1->atm[i].x=m2->atm[i].x;
      m1->atm[i].y=m2->atm[i].y;
      m1->atm[i].z=m2->atm[i].z;
      //m1->atm[i].rms=0; //m2->atm[i].rms;
      m1->atm[i].resnum=m2->atm[i].resnum;
      m1->atm[i].number=m2->atm[i].number;
      m1->atm[i].selected=TRUE;
      strcpy(m1->atm[i].name,m2->atm[i].name);
      strcpy(m1->atm[i].residue,m2->atm[i].residue);
      //printf("%-30s %lf %lf %lf == %lf %lf %lf\n",m2->filename,m1->atm[i].x,m1->atm[i].y,m1->atm[i].z,m2->atm[i].x,m2->atm[i].y,m2->atm[i].z);
    }
  m1->xcen=m2->xcen;
  m1->ycen=m2->ycen;
  m1->zcen=m2->zcen;
  m1->atoms=m2->atoms;
  printf("%d %d",m1->atoms,m2->atoms);
  strcpy(m1->filename,m2->filename);
}



//End
///////////////////////////
