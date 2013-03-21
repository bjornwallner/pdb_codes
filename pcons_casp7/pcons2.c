#include <stdio.h>
#include <stdlib.h>
//#include <stddef.h>
#include <sys/types.h>
#include <dirent.h>
#include <sys/times.h>
#include <time.h>
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h"
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.h"
//#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_molecule.h"
#define MAXFILE 2000
#define STRING_BUFFER 1000

void usage();
int SameMethod(char *file1, char *file2);
int ispdb(char *filename);
char *basename(char *fname);

main(int argc,char *argv[])             /* Main routine */
{

  
  DIR            *dip;
  FILE           *fp;
  struct dirent  *dit;
  char           dir[STRING_BUFFER]="";
  char           file2[STRING_BUFFER];
  char           tempfilename[STRING_BUFFER];
  char           infilename[STRING_BUFFER]="";
  char           filenames[MAXFILE][100];
  char           filenames_with_path[MAXFILE][STRING_BUFFER];
  char          target[6]="T0XXX";

  int            i,j,k;
  int            files=0;
  int            L=4;
  int            step=2; 
  int            number_of_comparisons[MAXFILE]={0};
  int            total_number_of_comparisons=0;
  int            maxlen=0;
  int           fastmode=0;
  int           memorymode=1;
  molecule_ca      m[MAXFILE];
  static clock_t st_time;
  static clock_t en_time;
  static struct tms st_cpu;
  static struct tms en_cpu;


  double         minsim=121.0;
  double         d0=sqrt(5);
  double         factor=0.5;
  char           temp[STRING_BUFFER];
  lgscore        LG[1];
  double         LG_average[MAXFILE]={0};
  double         S_average[MAXFILE]={0};
  double         Sstr[MAXFILE][MAXRES]={{0}};
  double         rmsd=0;
  double        sqrt5=sqrt(5);
 
}


void usage()
{
  
  fprintf(stderr,"********************************************************************************\n");
  fprintf(stderr,"*                               PCONS                                          *\n");
  fprintf(stderr,"* Calculates structral consensus for all models in a specified directory       *\n");
  fprintf(stderr,"*                                                                              *\n");
  fprintf(stderr,"* Reference: Bjorn Wallner and Arne Elofsson, Protein Sci., 2006 15(4):900-913 *\n");
  fprintf(stderr,"* For comments, please email to: bjorn@sbc.su.se                               *\n");
  fprintf(stderr,"*                                                                              *\n");
  fprintf(stderr,"*                                                                              *\n");
  fprintf(stderr,"********************************************************************************\n");
  fprintf(stderr,"\nUsage: pcons\t-d <directory w/ models>\n");
  fprintf(stderr,"\t\t-i <inputfile w/ full path to models>\n");
  fprintf(stderr,"\t\t-L <target sequence length, only needed if ALL models are shorter than the target sequence length>\n");
  //fprintf(stderr,"\t\t-f <fast mode, more sloppy heuristics for finding optimal sub alignments>\n");
  //fprintf(stderr,"\t\t-m <memory options\n");
  //fprintf(stderr,"\t\t    0 = re-read from disk\n");
  //fprintf(stderr,"\t\t    1 = read all coordinates in to memory (default)>\n");

  exit(0);

}




int SameMethod(char *file1, char *file2)
{
  int ident=0;
  int diffpos=0;
  int i=0;
  int len1,len2;
  len1=strlen(file1);
  len2=strlen(file2);

  if(len1 != len2)
    {
      return 0;
    }
  
  //All files here will have the same length
  //Check for identical matches
  
  for(i=0;i<len1;i++)
    {
      if(file1[i] == file2[i])
	{
	  ident++;
	}
      else
	{
	  diffpos=i+1;
	}
    }

  // printf("id: %d len %d pos %d ext %s\n",ident,len1,diffpos,&file1[len1-3]);
     
  if(len1-ident == 1 && (diffpos==len1 || diffpos==len1-4 && strcmp(&file1[len1-2],"pdb")))
    {
      //   printf("%s same as %s\n",file1,file2);

      return 1;

    }
  else
    {
      return 0;
    }
}


int ispdb(char *filename)
{
  molecule m[1];
  strcpy(m[0].filename,filename);
  //printf("ispdb: %s\n",filename);
  if(read_molecules(&m[0],'c')==0)
    {
      if(m[0].residues>0)
	{
	  return 1;
	}
    }

  return 0;
}
char *basename(char *fname)
{
  int sl;
    char *s;
    // if (strlen(fname) > STRING_BUFFER-4)
    //  error(NIL,"filename is too long",NIL);
    s = strrchr(fname,'/');
    if (s)
      fname = s+1;
    /* Process suffix */
    return fname;

}

