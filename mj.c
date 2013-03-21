#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//#ifdef  NOZLIB
//#else
//#include <zlib.h> 
//#endif
//#include <sys/types.h>
//#include <sys/times.h>
//#include <sys/param.h>
//#include <sys/time.h> 

#define	MAXATMS		10000		/* Maximum allowable atoms */
#define MAXRES          1000
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
    int number;
    int resnum;
    int rescount;
    int selected;
  } atm[MAXATMS];
  //double xcen,ycen,zcen;
  double res_coord_x[MAXRES];
  double res_coord_y[MAXRES];
  double res_coord_z[MAXRES];
  int	atoms;			/* Current # of atoms */
  int   residues;
  char	filename[1000];		/* filename to read molecule from */
} molecule;

molecule	m[1];		/* Molecule to be read*/
//double		LG_pvalue();
//double		Levitt_Gerstein();
//double		superimpose_molecules();
int             read_molecules();
int             get_type(char *name, char *res);
double          distance(int atomno1, int atomno2);
void            print_type(int type_no, FILE *fp);
int             get_res6(char *res);
int             get_res(char *res);
void            print_res(int res,FILE *fp);
double          residue_distance(int resno1, int resno2);   



main(argc,argv)		/* Main routine */
     int argc;
     char *argv[];
{
  
  int           i,j,k;
  double        dist;
  double        mean=0;
  int           atoms=0;
  int           residues=0;
  int           res_contacts[20][20]={{}};
  int           restype[20]={};
  int           tot_res_contacts=0;
  int           type_i,type_j;
  int           current_res_i=0;
  int           current_res_j=0;
  double        cutoff=42.25;
  FILE          *fp;
  double temp;
  char	buff[512];
  double E_c[20][20]={{}};
  float e[20][20]={{-5.44,-4.99,-5.80,-5.50,-5.83,-4.96,-4.95,-4.16,-3.57,-3.16,-3.11,-2.86,-2.59,-2.85,-2.41,-2.27,-3.60,-2.57,-1.95,-3.07},
		   {0.46,-5.46,-6.56,-6.02,-6.41,-5.32,-5.55,-4.91,-3.94,-3.39,-3.51,-3.03,-2.95,-3.30,-2.57,-2.89,-3.98,-3.12,-2.48,-3.45},
		   {0.54,-0.20,-7.26,-6.84,-7.28,-6.29,-6.16,-5.66,-4.81,-4.13,-4.28,-4.02,-3.75,-4.10,-3.48,-3.56,-4.77,-3.98,-3.36,-4.25},
		   {0.49,-0.01,0.06,-6.54,-7.04,-6.05,-5.78,-5.25,-4.58,-3.78,-4.03,-3.52,-3.24,-3.67,-3.17,-3.27,-4.14,-3.63,-3.01,-3.76},
		   {0.57,0.01,0.03,-0.08,-7.37,-6.48,-6.14,-5.67,-4.91,-4.16,-4.34,-3.92,-3.74,-4.04,-3.40,-3.59,-4.54,-4.03,-3.37,-4.20},
		   {0.52,0.18,0.10,-0.01,-0.04,-5.52,-5.18,-4.62,-4.04,-3.38,-3.46,-3.05,-2.83,-3.07,-2.48,-2.67,-3.58,-3.07,-2.49,-3.32},
		   {0.30,-0.29,0.00,0.02,0.08,0.11,-5.06,-4.66,-3.82,-3.42,-3.22,-2.99,-3.07,-3.11,-2.84,-2.99,-3.98,-3.41,-2.69,-3.73},
		   {0.64,-0.10,0.05,0.11,0.10,0.23,-0.04,-4.17,-3.36,-3.01,-3.01,-2.78,-2.76,-2.97,-2.76,-2.79,-3.52,-3.16,-2.60,-3.19},
		   {0.51,0.15,0.17,0.05,0.13,0.08,0.07,0.09,-2.72,-2.31,-2.32,-2.01,-1.84,-1.89,-1.70,-1.51,-2.41,-1.83,-1.31,-2.03},
		   {0.68,0.46,062,0.62,0.65,0.51,0.21,0.20,0.18,-2.24,-2.08,-1.82,-1.74,-1.66,-1.59,-1.22,-2.15,-1.72,-1.15,-1.87},
		   {0.67,0.28,0.41,0.30,0.40,0.36,0.37,0.13,0.10,0.10,-2.12,-1.96,-1.88,-1.90,-1.80,-1.74,-2.42,-1.90,-1.31,-1.90},
		   {0.69,0.53,0.44,0.59,0.60,0.55,0.38,0.14,0.18,0.14,-0.06,-1.67,-1.58,-1.49,-1.63,-1.48,-2.11,-1.62,-1.05,-1.57},
		   {0.97,0.62,0.72,0.87,0.79,0.77,0.30,0.17,0.36,0.22,0.02,0.10,-1.68,-1.71,-1.68,-1.51,-2.08,-1.64,-1.21,-1.53},
		   {0.64,0.20,0.30,0.37,0.42,0.46,0.19,-0.12,0.24,0.24,0.08,0.11,-0.10,-1.54,-1.46,-1.42,-1.98,-1.80,-1.29,-1.73},
		   {0.91,0.77,0.75,0.71,0.89,0.89,0.30,-0.07,0.26,0.13,-0.14,-0.19,-0.24,-0.09,-1.21,-1.02,-2.32,-2.29,-1.68,-1.33},
		   {0.91,0.30,0.52,0.46,0.55,0.55,0.00,-0.25,0.30,0.36,-0.22,-0.19,-0.21,-0.19,0.05,-0.91,-2.15,-2.27,-1.80,-1.26},
		   {0.65,0.28,0.39,0.66,0.67,0.70,0.08,0.09,0.47,0.50,0.16,0.26,0.29,0.31,-0.19,-0.16,-3.05,-2.16,-1.35,-2.25},
		   {0.93,0.38,0.42,0.41,0.43,0.47,-0.11,-0.30,0.30,0.18,-0.07,-0.01,-0.02,-0.26,-0.91,-1.04,0.14,-1.55,-0.59,-1.70},
		   {0.83,031,033,032,037,033,-0.10,-0.46,0.11,0.03,-0.19,-0.15,-0.30,-0.46,-1.01,-1.28,0.23,0.24,-0.12,-0.97},
		   {0.53,0.16,0.25,0.39,0.35,0.31,-0.33,-0.23,0.20,0.13,0.04,0.14,0.18,-0.08,0.14,0.07,0.15,-0.05,-0.04,-1.75}};
  float q[20]={6.646,6.137,5.870,6.042,6.087,6.155,5.793,6.037,6.334,6.284,6.486,6.582,6.574,6.469,6.487,6.235,6.241,6.318,6.569,5.858};
  float N[8][20]={{212.1,167.2,415.4,527.5,720.4,654.2,185.6,425.0,754.2,677.2,485.1,470.1,333.5,243.9,364.2,306.2,204.6,338.7,251.8,320.8},
		  {242.1,209.7,461.4,619.4,947.3,804.8,172.7,402.6,780.6,499.2,366.8,378.5,201.1,172.1,204.5,204.6,152.0,219.3,125.4,208.4},
		  {178.7,194.4,368.8,549.7,921.2,686.4,95.5,258,494.4,288,235.1,237.9,119.9,90.8,143.9,108.4,95.4,98.5,53.1,116.8},
		  {86.4,107.9,213.0,318.6,553.6,392.8,53.3,121.2,213.0,115.8,116.1,88.5,47.8,52.0,70.0,40.7,55.4,33.8,17.3,50.4},
		  {20.3,42.7,69.2,110.0,205.8,124.7,14.9,36.4,57.8,35.9,45.2,27.3,6.5,13.7,16.2,12.4,14.2,11.4,4.9,13.8},
		  {4.5,6.1,23.7,28.8,59.9,30.0,2.9,11.8,11.6,16.0,9.5,7.1,1.9,1.1,4.9,0.7,2.9,0.6,1.7,1.1},
		  {0.5,0.0,2.5,4.9,5.6,4.7,0,0.8,2.8,1.9,1.8,0.1,0,0,0.1,0,2.0,0,0,0.8},
		  {0,0,0.2,0,0,0.2,0,0,0.2,0.2,0,0,0,0,0,0,0,0,0,0}};
    
  float N_q[20]={201.1,207.6,455.4,616.5,945.0,785.9,175.3,397.2,685,428.7,302.8,296.7,154.5,134,175,182,138.4,180.9,84.26,224.4};
  double E=0;

  /* Parse command line for PDB filename */
  if(argc==2)
    {

      strcpy(m[0].filename,argv[1]);
      //cutoff=strtod(argv[2],(char**)(argv[2]+strlen(argv[2])));
      //cutoff=cutoff*cutoff;
    }
  else
    {
      printf("Usage: mj [pdb_file]\n");
      exit(1);
    }
  
  if(read_molecules()==0)
    {  
      for(i=0;i<m[0].atoms;i++)  
	{
	  if(strcmp("CA ",m[0].atm[i].name)==0)
	  {
	    for(j=0;j<m[0].atoms;j++)
	      {
		if(strcmp("CA ",m[0].atm[j].name)==0)
		  {
		    if(abs(m[0].atm[i].rescount-m[0].atm[j].rescount)>1 &&
		       residue_distance(m[0].atm[i].rescount,m[0].atm[j].rescount)<=cutoff)
		      {
			/*printf("%d %d %s %d %d %s %f %f %f %f %f %f %f %f\n",
			       m[0].atm[i].resnum,
			       m[0].atm[i].rescount,
			        m[0].atm[i].residue,
			        m[0].atm[j].resnum,
			        m[0].atm[j].rescount,
			        m[0].atm[j].residue,
			        residue_distance(m[0].atm[i].rescount,m[0].atm[j].rescount),
			        residue_distance(m[0].atm[j].rescount,m[0].atm[i].rescount),
			        m[0].res_coord_x[m[0].atm[i].rescount],
			        m[0].res_coord_y[m[0].atm[i].rescount],
			        m[0].res_coord_z[m[0].atm[i].rescount],
			        m[0].res_coord_x[m[0].atm[j].rescount],
			        m[0].res_coord_y[m[0].atm[j].rescount],
			        m[0].res_coord_z[m[0].atm[j].rescount]);*/
			res_contacts[get_res(m[0].atm[i].residue)][get_res(m[0].atm[j].residue)]++;
			tot_res_contacts++;
		      }
		  }
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
      /*printf("      ");
       for(i=0;i<20;i++)
 	{
	   print_res(i,stdout);
	   printf("  ");
	  
	 }
	   printf("\n");
      
	  for(i=0;i<20;i++)
	{ 
	   print_res(i,stdout);
	  for(j=0;j<20;j++)
	     {
	       if(i!=j)
		 {
		   res_contacts[i][j]=2*res_contacts[i][j];
		} 
	      printf("%5d",res_contacts[i][j]);
	    } 
	   printf("\n");
	 }
	printf("\n\n");*/
      for(i=0;i<20;i++)
	{
	  for(j=i;j<20;j++)
	    {
	      if(i!=j)
		{
		  res_contacts[i][j]=2*res_contacts[i][j];
		}
	      E_c[i][j]=(double)0.5*e[i][j]*res_contacts[i][j];
	      E+=E_c[i][j];
	      //	      printf("%4.3lf  ",E_c[i][j]);
	    }
	  //printf("\n");
	}
      printf("%lf",E);



   
      /*RESIDUE CONTACTS for each residue type */
      //      printf("RES: ");
      //for(i=0;i<6;i++)
      //	{
      //	  for(j=i;j<6;j++)
      //	    {
      //	      if(i!=j)
      //		{
      //		  res_contacts[i][j]=2*res_contacts[i][j];
      //		}
      //	    
      //	      
      //	      //restype[i]=restype[i]+res_contacts[i][j];
      //	      temp=0;
      //	      if(tot_res_contacts != 0)
      //		{
      //		  temp=(double)res_contacts[i][j]/tot_res_contacts;
      //		}
      //	      printf("%f ",temp);
      //	      //printf("%d ",res_contacts[i][j]);
      //	    }
      //	 }
    }
}




int read_molecules(temp)	/* Reads in molecules to be superimposed */
   int temp;
{
  int	i,j,k,atoms,residues,atom_count; /* Counter variables */
  char	buff[512];	/* Input string */
  char	junk[512];	/* Throwaway string */
  char	line_flag[10];	/* PDB file line mode */
  char	residue[8];	/* PDB atom info for output */
  char	name[8];
  int   resnum;
  char  resname[8];
  char  old_resname[8]="noname";
  char x_temp[10];
  char y_temp[10];
  char z_temp[10];
  int number;
  char chain[1];
  char inscode[1]=" ";
  char alt_loc[1]=" ";
  double	x,y,z;		/* Temporary coordinates values */
  double xcen_sidechain,ycen_sidechain,zcen_sidechain;

  FILE *fp;
  char temp_number[10];
  char temp_resnum[10];
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
      xcen_sidechain=ycen_sidechain=zcen_sidechain=atom_count=0;

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
	      //printf("%s\n",&buff[6]);
	      strncpy(temp_number,&buff[6],5);
	      strncpy(name,&buff[13],3);
	      strncpy(alt_loc,&buff[16],1);
	      strncpy(residue,&buff[17],3);
	      strncpy(chain,&buff[21],1);
	      strncpy(temp_resnum,&buff[22],4);
	      strncpy(resname,&buff[22],5);
	      strncpy(x_temp,&buff[30],8);
	      strncpy(y_temp,&buff[38],8);
	      strncpy(z_temp,&buff[46],8);

	      number=atoi(temp_number);
	      resnum=atoi(temp_resnum);
	      x=atof(x_temp);
	      y=atof(y_temp);
	      z=atof(z_temp);
	      
	      //printf("test: %s %d %s %s %s %s %s %s %lf %lf %lf\n",line_flag,number,name,alt_loc,residue,chain,resnum,resname,x,y,z);
	      //if (strcmp("N",name)==0)	   /*  Is it an N atom => new residue? */
	      //printf("%s %s\n",old_resname,resname);
	      if(strcmp(old_resname,resname)!=0 || strcmp("OXT",name)==0)
		{
		  residues++;
		  if(residues>1)
		    {
		      //Save coordinate for center of sidechain
		      m[0].res_coord_x[residues-1]=xcen_sidechain/atom_count;
		      m[0].res_coord_y[residues-1]=ycen_sidechain/atom_count;
		      m[0].res_coord_z[residues-1]=zcen_sidechain/atom_count;
		      xcen_sidechain=ycen_sidechain=zcen_sidechain=atom_count=0;
		      //printf("%s %d %lf %lf %lf\n",old_resname,residues-1,m[0].res_coord_x[residues-1],m[0].res_coord_y[residues-1],m[0].res_coord_z[residues-1]);
		    }

		  //printf("%s %s\n",resname,residue);
		}
	      //sscanf(&buff[22],"%d %lf %lf %lf",&resnum,&x,&y,&z);
	      //printf("test %d %s %s %d %d %lf %lf %lf\n",number,name,residue,resnum,residues,x,y,z);
	      m[0].atm[atoms].x=x;
	      m[0].atm[atoms].y=y;
	      m[0].atm[atoms].z=z;
	      m[0].atm[atoms].resnum=resnum;
	      m[0].atm[atoms].number=number; //atoi(number);
	      m[0].atm[atoms].rescount=residues;
	      m[0].atm[atoms].selected=TRUE;
	      strcpy(m[0].atm[atoms].name,name);
	      strcpy(m[0].atm[atoms].residue,residue);
	      //m[0].xcen+=x;
	      //m[0].ycen+=y;
	      //m[0].zcen+=z;
	      atoms++;
	      strcpy(old_resname,resname);
	
	      if((strcmp("C  ",name)!=0 && 
		  strcmp("O  ",name)!=0 &&
		  strcmp("N  ",name)!=0 &&
		  strcmp("CA ",name)!=0 &&
		  strcmp("OXT",name)!=0) ||
		 (strcmp("GLY",residue)==0 &&
		  strcmp("CA ",name)==0))
		{
		  	
		  //printf("test %d %s %s %d %lf %lf %lf\n",number,name,residue,resnum,x,y,z);
		  xcen_sidechain+=x;
		  ycen_sidechain+=y;
		  zcen_sidechain+=z;
		  atom_count++;
		}

	    }
	}
      m[0].atoms=atoms;
      m[0].residues=residues;
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
     
int get_type(char *name, char *res)        /* Function which takes a name and residue and return type */
{
  /* Types  - Return number
   *
   * C         - 0
   * N         - 1
   * O         - 2
   * CA        - 3
   * CH3       - 4
   * CH/CH2    - 5
   * C(OO-)    - 6
   * NH        - 7
   * NH2       - 8
   * (C)OO-    - 9
   * =O        - 10
   * OH        - 11
   * S         - 12    
   * OXT       - 13
   */
  //printf("%s %s\n",name,res);
  if(strcmp("C  ",name)==0)
    return 0;
  if(strcmp("N  ",name)==0)
    return 1;
  if(strcmp("O  ",name)==0)
    return 2;
  if(strcmp("CA  ",name)==0)
    return 3;
  if((strcmp("ILE",res)==0 && (strcmp("CD1",name)==0 || strcmp("CG2",name)==0)) ||
     (strcmp("LEU",res)==0 && (strcmp("CD1",name)==0 || strcmp("CD2",name)==0)) ||
     (strcmp("MET",res)==0 && strcmp("CE ",name)==0) ||
     (strcmp("THR",res)==0 && strcmp("CG2",name)==0) ||
     (strcmp("ALA",res)==0 && strcmp("CB ",name)==0) ||
     (strcmp("VAL",res)==0 && (strcmp("CG1",name)==0 || strcmp("CG2",name)==0)))
    return 4; 
  if((strcmp("ASP",res)==0 && strcmp("CG ",name)==0) ||
     (strcmp("GLU",res)==0 && strcmp("CG ",name)==0))
    return 6;
  if(strcmp("NE ",name)==0 || strcmp("ND1",name)==0 || strcmp("NE1",name)==0)
    return 7;
  if(strcmp("NH1",name)==0 || strcmp("NH2",name)==0 || strcmp("ND2",name)==0
     || strcmp("NE2",name)==0 || strcmp("NZ ",name)==0)
    return 8;
  if((strcmp("ASP",res)==0 && (strcmp("OD2",name)==0 || strcmp("OD1",name)==0)) ||
     (strcmp("GLU",res)==0 && (strcmp("OE1",name)==0 || strcmp("OE2",name)==0)))
    return 9;
  if((strcmp("ASN",res)==0 && strcmp("OD1",name)==0) ||
     (strcmp("GLN",res)==0 && strcmp("OE1",name)==0))
    return 10; 
  if(strcmp("OG ",name)==0 || strcmp("OH ",name)==0 || strcmp("OG1",name)==0)
    return 11;
  if(strcmp("SG ",name)==0 || strcmp("SD ",name)==0)
    return 12;
  if(strcmp("OXT",name)==0)
    return 13;
  return 5;
}

void print_type(int type_no, FILE *fp)
{
  if(type_no==0)
    fprintf(fp,"C");
  if(type_no==1)
    fprintf(fp,"N");
  if(type_no==2)
    fprintf(fp,"O");
  if(type_no==3)
    fprintf(fp,"CA");
  if(type_no==4)
    fprintf(fp,"CH3");
  if(type_no==5)
    fprintf(fp,"CH/CH2");
  if(type_no==6)
    fprintf(fp,"C(OO-)");
  if(type_no==7)
    fprintf(fp,"NH");
  if(type_no==8)
    fprintf(fp,"NH2");
  if(type_no==9)
    fprintf(fp,"(C)OO-");
  if(type_no==10)
    fprintf(fp,"=0");
  if(type_no==11)
    fprintf(fp,"OH");
  if(type_no==12)
    fprintf(fp,"S");
  if(type_no==13)
    fprintf(fp,"OXT");
}

int get_res6(char* res)
{
     if(strcmp("ARG",res)==0 || strcmp("LYS",res)==0)
       return 0;
     if(strcmp("ASP",res)==0 || strcmp("GLU",res)==0)
       return 1;
     if(strcmp("HIS",res)==0 || strcmp("PHE",res)==0 || strcmp("TRP",res)==0 || strcmp("TYR",res)==0)
       return 2;
     if(strcmp("ASN",res)==0 || strcmp("GLN",res)==0 || strcmp("SER",res)==0 || strcmp("THR",res)==0)
       return 3;
     if(strcmp("ALA",res)==0 || strcmp("ILE",res)==0 || strcmp("LEU",res)==0 || strcmp("MET",res)==0 || strcmp("VAL",res)==0 || strcmp("CYS",res)==0)
       return 4;
     if(strcmp("GLY",res)==0 || strcmp("PRO",res)==0)
       return 5;
           
}
int get_res(char* res)
{

  if(strcmp("CYS",res)==0)
    return 0;
  if(strcmp("MET",res)==0)
    return 1;
  if(strcmp("PHE",res)==0)
    return 2;
  if(strcmp("ILE",res)==0)
    return 3;
  if(strcmp("LEU",res)==0)
    return 4;
  if(strcmp("VAL",res)==0)
    return 5;
  if(strcmp("TRP",res)==0)
    return 6;
  if(strcmp("TYR",res)==0)
    return 7;
  if(strcmp("ALA",res)==0)
    return 8;
  if(strcmp("GLY",res)==0)
    return 9;
  if(strcmp("THR",res)==0)
    return 10;
  if(strcmp("SER",res)==0)
    return 11;
  if(strcmp("ASN",res)==0)
    return 12;
  if(strcmp("GLN",res)==0)
    return 13;
  if(strcmp("ASP",res)==0)
    return 14;
  if(strcmp("GLU",res)==0)
    return 15;
  if(strcmp("HIS",res)==0)
    return 16;
  if(strcmp("ARG",res)==0)
    return 17;
  if(strcmp("LYS",res)==0)
    return 18;
  if(strcmp("PRO",res)==0)
    return 19;
}
void print_res(int res,FILE *fp)
{
  if(res==0)
    fprintf(fp,"CYS");
  if(res == 1)
    fprintf(fp,"MET");
  if(res == 2)
    fprintf(fp,"PHE");
  if(res == 3)
    fprintf(fp,"ILE");
  if(res == 4)
    fprintf(fp,"LEU");
  if(res == 5)
    fprintf(fp,"VAL");
  if(res == 6)
    fprintf(fp,"TRP");
  if(res == 7)
    fprintf(fp,"TYR");
  if(res == 8)
    fprintf(fp,"ALA");
  if(res == 9)
    fprintf(fp,"GLY");
  if(res == 10)
    fprintf(fp,"THR");
  if(res == 11)
    fprintf(fp,"SER");
  if(res == 12)
    fprintf(fp,"ASN");
  if(res == 13)
    fprintf(fp,"GLN");
  if(res == 14)
    fprintf(fp,"ASP");
  if(res == 15)
    fprintf(fp,"GLU");
  if(res == 16)
    fprintf(fp,"HIS");
  if(res == 17)
    fprintf(fp,"ARG");
  if(res == 18)
    fprintf(fp,"LYS");
  if(res == 19)
    fprintf(fp,"PRO");

}
double residue_distance(int resno1, int resno2)       /* atomnoX is the first atom of the residue */
{
  return (m[0].res_coord_x[resno1]-m[0].res_coord_x[resno2])*(m[0].res_coord_x[resno1]-m[0].res_coord_x[resno2])+
	 (m[0].res_coord_y[resno1]-m[0].res_coord_y[resno2])*(m[0].res_coord_y[resno1]-m[0].res_coord_y[resno2])+
	 (m[0].res_coord_z[resno1]-m[0].res_coord_z[resno2])*(m[0].res_coord_z[resno1]-m[0].res_coord_z[resno2]);

}

double distance(int atomno1,int atomno2)
{
  //printf("%s %s\n",m[0].atm[atomno1-1].name,m[0].atm[atomno2-1].name);

  return (m[0].atm[atomno1].x-m[0].atm[atomno2].x)*(m[0].atm[atomno1].x-m[0].atm[atomno2].x)+(m[0].atm[atomno1].y-m[0].atm[atomno2].y)*(m[0].atm[atomno1].y-m[0].atm[atomno2].y)+(m[0].atm[atomno1].z-m[0].atm[atomno2].z)*(m[0].atm[atomno1].z-m[0].atm[atomno2].z);
}





