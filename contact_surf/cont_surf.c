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

//double		Levitt_Gerstein();
//double		superimpose_molecules();
//void            strncpy_NULL(char *dest, char *src, size_t n);
//int             read_molecules();
//int             get_type(char *name, char *res);
//double          distance(int atomno1, int atomno2);
//void            print_type(int type_no, FILE *fp);
//int             get_res6(char *res);
//int             get_res(char *res);
//void            print_res(int res,FILE *fp);
//double          crd(int atomno1, int atomno2);   /*closest residue distance */


main(int argc,char *argv[])		/* Main routine */
     //int argc;
     //char *argv[];
{
  molecule      m[1];
  int           i,j,k;
  double        dist;
  double        mean=0;
  int           atoms=0;
  int           residues=0;
  int           res_contacts[6][6]={{}};
  int           atom_i=0;
  int           atom_j=0;
  int           contacts[2000]={};
  int           exposed_residues[6]={};
  int           exposed_cut;
  int           buried_cut;
  int           buried_residues[6]={};
  int           total_exposed=0;
  int           total_buried=0;
  int           restype[6]={};
  int           tot_res_contacts=0;
  double        frac_res_contacts[21]={};   // 6*7/2=21
  double        frac_exposed_residues[6]={};
  double        frac_buried_residues[6]={};
  double        fat=0;
  double        cutoff=196;
  double         hcut;
  double         hangle;
  FILE          *fp;
  double        temp;
  int           tmp=0;
  float         sum=0;
  int           sum2=0;
  char          *ss;
  char	        psipredfiles[1000];
  char          psipred[2000];
  char          pdbfiles[1000];
  //char          *psipred;
  int           ss_correct=0;
  double        Q3=0;
  int           ss_pred=0;
  int           ss_real=0;

  /* Parse command line for PDB filename */
//  if(argc==4)
//    {
//	strcpy(m[0].filename,argv[1]);
//	//cutoff=strtod(argv[2],(char**)(argv[2]+strlen(argv[2])));
//	//cutoff=cutoff*cutoff;
//	hcut=atof(argv[2]); //,(char**)(argv[2]+strlen(argv[2]))); //,(char**)(argv[2]+strlen(argv[2])));
//	hangle=atof(argv[3]); //,(char**)(argv[3]+strlen(argv[3])));
//	//printf("%lf %lf\n",hcut,hangle);
//	//for(i=1;i<argc;i++)
//	// 	printf("%s\n",argv[i]);
//    }
//  else
  //    {
  if(argc==3)
    {
      //strcpy(m[0].filename,argv[1]);
      strcpy(pdbfiles,argv[1]);
      //exposed_cut=atoi(argv[2]);
      //buried_cut=atoi(argv[3]);
      //cutoff=cutoff*cutoff;
      strcpy(psipredfiles,argv[2]);
      
      //cutoff=atof(argv[3]);
      //cutoff=cutoff*cutoff;
      //printf("%d %d\n", exposed_cut,buried_cut);
      hcut=3.6;
      hangle=1.2;
      //exit(1);
      printf("%s\n%s\n",pdbfiles,psipredfiles);
    }
  else
    {
      printf("Usage: read_pdb [pdb_file] [psipredfile]\n");
      exit(1);
    }
  //printf("Reading...\n");
  if(read_molecules_backbone(m)==0)
    {  
      //printf("%d\n",m[0].residues);
      for(i=0;i<m[0].residues;i++)  
	{
	  for(j=i+1;j<m[0].residues;j++)
	    {
	      atom_i=m[0].CA_ref[i];
	      atom_j=m[0].CA_ref[j];
	      dist=distance(m,atom_i,atom_j);
	      if(dist<100)
		{
		  contacts[m[0].atm[atom_i].rescount-1]++;
		  contacts[m[0].atm[atom_j].rescount-1]++;
		}
	      if(abs(m[0].atm[atom_i].rescount-m[0].atm[atom_j].rescount)>5 &&
		 dist<cutoff)
		{
		  res_contacts[get_res6(m[0].atm[atom_i].residue)][get_res6(m[0].atm[atom_j].residue)]++;
		  res_contacts[get_res6(m[0].atm[atom_j].residue)][get_res6(m[0].atm[atom_i].residue)]++;
		  tot_res_contacts=tot_res_contacts+2;
		}
	    }
	}
      for(i=0;i<m[0].residues;i++)
	{
	  if(contacts[i]<16)  //exposed
	    {
	      exposed_residues[get_res6(m[0].atm[m[0].CA_ref[i]].residue)]++;
	      total_exposed++;
			      
	    }
	  if(contacts[i]>20)  //buried
	    {
	      
	      buried_residues[get_res6(m[0].atm[m[0].CA_ref[i]].residue)]++;
	      total_buried++;
	    }
	}
      //printf("Exposed Buried\n");
      for(i=0;i<6;i++)
	{
	  //printf("%d/%d %d/%d\n",exposed_residues[i],total_exposed,buried_residues[i],total_buried);
	  if(total_exposed!=0)
	    frac_exposed_residues[i]=(double)exposed_residues[i]/total_exposed;
	  if(total_buried!=0)
	    frac_buried_residues[i]=(double)buried_residues[i]/total_buried;
	  //printf("%5.4f %5.4f\n",frac_exposed_residues[i],frac_buried_residues[i]);
	}
      
//	for(i=0;i<6;i++)
//	  {
//	    printf("%f ",frac_exposed_residues[i]);
//	  }
//	for(i=0;i<6;i++)
//	  {
//	    printf("%f ",frac_buried_residues[i]);
//	  }
//	exit(0);

      for(i=0;i<6;i++)
	{
	  for(j=0;j<6;j++)
	    {
	      if(i!=j)
		{
		  res_contacts[i][j]=2*res_contacts[i][j];
		}
	      //  printf("%3d ",res_contacts[i][j]);
	    }
	  //  printf("\n");
	}
      for(i=0;i<6;i++)
	{ 
	  for(j=i;j<6;j++) 
	    {
	      //if(i!=j)
	      //		{
	      //	  res_contacts[i][j]=2*res_contacts[i][j];
	      //	}
	      //restype[i]=restype[i]+res_contacts[i][j];
	      // temp=0;
	      if(tot_res_contacts != 0)
		{
		  frac_res_contacts[k]=(double)res_contacts[i][j]/tot_res_contacts;
		  //temp=(double)res_contacts[i][j]/tot_res_contacts;
		}
	      //printf("%f ",temp);
	      //sum+=temp;
	      k++;
	    }
	}
      //printf("\n");
    }
  fat=fatness(m);
  //printf("%7.6f\n",fat);

  ss=assign_ss(m,hcut,hangle);
  //printf("SEQ: %s\nSS: %s\n",m[0].sequence,ss);
  //psipred=read_psipred(psipredfile);
  //printf("PSI %s\n",psipred);


  if(strlen(psipred)==strlen(ss))
    {
      for(i=0;i<strlen(psipred);i++)
      {
	//printf("%d %d\n",psipred[i],ss[i]);
	//printf("%d\n",(psipred[i]==ss[i]));
	if(psipred[i]==ss[i])
	  {
	    ss_correct++;
	  }
      }
      //printf("Q3: %d %5.3lf\n",ss_correct,(double)ss_correct/strlen(psipred));
      Q3=(double)ss_correct/strlen(psipred);
    
      for(i=0;i<21;i++)
	{
	  printf("%f ",frac_res_contacts[i]);
	}
      for(i=0;i<6;i++)
	{
	  printf("%f ",frac_exposed_residues[i]);
	}
      for(i=0;i<6;i++)
	{
	  printf("%f ",frac_buried_residues[i]);
	}
      printf("%f %f",fat,Q3);
      //printf("\n");
     


    }
     else
    {
      printf("Different strlen!!!! %d %d\n",strlen(psipred),strlen(ss));
     
    }
  free(ss);
}



