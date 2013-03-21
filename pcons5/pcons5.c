#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/molecule.h"
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_measure/lgscore.h"
#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/pcons5/pcons5.h"

//#include "/afs/pdc.kth.se/home/b/bjornw/bjorn/source/c/pdb/quality_molecule.h"

main(int argc,char *argv[])             /* Main routine */
{
  molecule m[1];
   FILE *fp;
   double score[MAXFILE];
   double proq_lg[MAXFILE];
   double proq_mx[MAXFILE];
   double LG;
   double c1,c2;
   double        *quality;
   double jury_all=0;
   double jury_first=0;
   double pcons5=0;
   double sum_all[MAXFILE]={0};
   double sum_first[MAXFILE]={0};
   double pscore1_5[MAXFILE]={0};
   double pscore3[MAXFILE]={0};
   int count_all[MAXFILE]={0};
   int count_first[MAXFILE]={0};
   int i=0,j=0;
   int files=0;
   int  number_of_methods=0;
   int tmp;
   int rank[MAXFILE];
   char filenames[MAXFILE][200];
   char filename[200];
   char method[MAXFILE][30];
   char cutofffile[200]="/afs/pdc.kth.se/home/b/bjornw/www/bin/cutoff.init";
   

   fp=fopen(cutofffile,"r");
   if(fp==NULL)
     {
       if(getenv("PCONS5DIR")==NULL)
	 {
	   fprintf(stderr,"ENV variable PCONS5DIR needs to be set properly\n");
	   exit(1);
	 }
       strcpy(cutofffile,getenv("PCONS5DIR"));
       strcat(cutofffile,"cutoff.init");
     }

   //printf("file used: %s\n",cutofffile);


//   char *methods[]={"inbgu",
//		    "burnham",
//		    "pdbblast",
//		    "fugue",
//		    "fugu3",
//		    "fugsa",
//		    "mgenthreader",
//		    "samt99",
//		    "foldfit",
//		    "orfeus",
//		    "superfamily",
//		    "orfblast",
//		    "supfampp",
//		    "samt02",
//		    "raptor",
//		    "ffas03",
//		    '\0'};
//
////my @methods_to_use=('burnham',
////		    'ffas03', 
////		    'foldfit',
////		    'forte1', 
////		    'fugsa',
////		    'fugu3',
////		    'inbgu',
////		    'mgenthreader',
////		    'orfblast',
////		    'orfeus',
////		    'pdbblast',
////		    'prospect2',
////		    'raptor',
////		    'supfampp',
////		    'samt99',
////		    'samt02');
//
//   double cutoffs90_1_5[]={23.6,
//			   6.35,
//			   2.301,
//			   5.29,
//			   8.38,
//			   8.27,
//			   0.588,
//			   18.15,
//			   1.3883,
//			   11.75,
//			   1.1979,
//			   3,
//                           10.1733,
//			   0.61995,
//			   6.33,
//			   10.4};
//
//
//   double cutoffs50_1_5[]={4.1,
//			   4.64,
//			   -0.63347,
//			   3.24,
//			   3.32,
//			   3.31,
//			   0.493,
//			   7.68,
//			   -0.50106,
//			   2.68,
//			   -2.2648,
//                           -0.041393,
//                           -1.8319,
//			   -1.2372,
//			   4.6,
//			   5.3};
//
//
//   double cutoffs90_3[]={140.6,
//			 999,
//			 999,
//			 999,
//			 23.65,
//			 28.31,
//			 0.977,
//			 121.9,
//			 9.9508,
//			 999,
//			 29.8239,
//                         54.301,
//			 41.082,
//			 18.3591,
//			 38.91,
//			 51};
//
//
//   double cutoffs50_3[]={27,
//			 12.28,
//			 5.3979,
//			 7.54,
//			 14,
//			 13.49,
//			 0.916,
//			 15.09,
//			 0.96257,
//			 12.78,
//			 1.5272,
//                         16.5229,
//			 10.6289,
//			 2.1386,
//			 9.57,
//			 21.6};
   cutoffs cut;
   
   readcutoffs(&cut,cutofffile);

   //exit(1);
   // Count number of methods
   //   for(number_of_methods=0;methods[number_of_methods]!=NULL;number_of_methods++);
    for(number_of_methods=0;cut.methods[number_of_methods].method[0]!=NULL;number_of_methods++)
     {
       // printf("%s %lf \n",cut.methods[number_of_methods].method,cut.methods[number_of_methods].cutoffs50_3);
     }
    //exit(1);
   //  {
   //    printf("%s %lf %lf %lf  %lf\n",methods[number_of_methods],cutoffs50_1_5[number_of_methods],cutoffs50_3[number_of_methods],cutoffs90_1_5[number_of_methods],cutoffs90_3[number_of_methods]);

   //  }
   if(argc==2)
      {
	strcpy(filename,argv[1]);
	fp=fopen(filename,"r");
	if(fp!=NULL)
	  {
	    while(fscanf(fp,"%s",filenames[i])!=EOF)
	      {
			//		strcpy(filenames[i],buff);
		//tmp=ftell(fp);
		//printf("%s %d %d\n",filenames[i],i,tmp);
		i++;
	      }
	    files=i-1;
	  }
	else
	  {
	    printf("Cannot open %s",filename);
	    exit(1);
	  }
	fclose(fp);
      }
    else
      {
	printf("\nPCONS5 CONSENSUS PREDICTOR\n\nUsage:    pcons5 infile (file with full path to all models)\n\nOBS the PDB should contain a field SS with the secondary structure assignment!\nEx: SS    HHHHHHCCCCCCCCEEEEE\n\n");
	exit(1);
      }
   //printf("Calculating ProQ and Pscores ....\n");
   //printf("%d",files);
   //exit(1);
    for(i=0;i<=files;i++)
      {
	strcpy(m[0].filename,filenames[i]);
	m[0].filename[strlen(filenames[i])]='\0';
	//strncpy_NULL(m[0].filename,filenames[i],strlen(filenames[i]));
	strncpy_NULL(m[0].ss,"",1);
	
	//TETRA	read_molecules_backbone(&m[0]);
	read_molecules(&m[0],'b');

	//printf("%s %s %d\n",filenames[i],m[0].filename,strlen(filenames[i]));
	//	printf("%s\t%d\t%s\t%lf\n",filenames[i],m[0].rank,m[0].method,m[0].score);
	strcpy(method[i],m[0].method);
	rank[i]=m[0].rank;
	score[i]=m[0].score;
	//printf("\n%s\t%d\t%s\t%lf\n",m[0].filename,m[0].rank,m[0].method,m[0].score);
	if(strcmp(method[i],"ffas03")==0)
	  {
	    score[i]=-score[i];
	  }

	if(m[0].score!=0 &&
	   (strcmp(method[i],"blast")==0 ||
	    strcmp(method[i],"pdbblast")==0 ||
	    strcmp(method[i],"foldfit")==0 ||
	    strcmp(method[i],"superfamily")==0 ||
	    strcmp(method[i],"supfampp")==0 ||
	    strcmp(method[i],"samt02")==0 ||
	    strcmp(method[i],"orfblast")==0||
	    strcmp(method[i],"prob_score")==0))
	  {
	    //	  printf("%lf %lf\n",score[i],log10(score[i]));
	    score[i]=log10(1/score[i]);
	    //printf("%lf %lf\n",m[0].score,log(m[0].score));
	    // printf("%lf %lf %lf %lf\n",m[0].score,log(m[0].score),0,0);
	  }


	//printf("%s\n",method[i]);
	//Find method set it index j
	for(j=0;j<number_of_methods;j++)
	  {
	    if(strcmp(method[i],cut.methods[j].method)==0)
	      break;
	  }
	//printf("%s\n",cut.methods[j].method);
	//exit(1);
	pscore1_5[i]=0;
	pscore3[i]=0;

	if(j!=number_of_methods)   //method exist in methods list and j is its index.
	  {
	    if(score[i]>cut.methods[j].cutoffs50_1_5 && cut.methods[j].cutoffs90_1_5 != 999)
	      {
		c1=cut.methods[j].cutoffs50_1_5;
		c2=cut.methods[j].cutoffs90_1_5;

		pscore1_5[i]=(score[i]-c1)/(c2-c1);
		if(pscore1_5[i]>1)
		  pscore1_5[i]=1;
		
		if(score[i]>cut.methods[j].cutoffs50_3 && cut.methods[j].cutoffs90_3 != 999)
		  {
		    c1=cut.methods[j].cutoffs50_3;
		    c2=cut.methods[j].cutoffs90_3;
		    pscore3[i]=(score[i]-c1)/(c2-c1);
		    if(pscore3[i]>1)
		      pscore3[i]=1;
		  }
	      }
	  }

//	if(j!=number_of_methods)   //method exist in methods list and j is its index.
//	  {
//	    if(score[i]>cutoffs50_1_5[j])
//	      {
//		c1=cutoffs50_1_5[j];
//		c2=cutoffs90_1_5[j];
//
//		pscore1_5[i]=(score[i]-c1)/(c2-c1);
//		if(pscore1_5[i]>1)
//		  pscore1_5[i]=1;
//		
//		if(score[i]>cutoffs50_3[j] && cutoffs90_3[j] != 999)
//		  {
//		    c1=cutoffs50_3[j];
//		    c2=cutoffs90_3[j];
//		    pscore3[i]=(score[i]-c1)/(c2-c1);
//		    if(pscore3[i]>1)
//		      pscore3[i]=1;
//		  }
//	      }
//	  }
		
	//printf("HERE %s %lf %lf\n",filenames[i],pscore1_5[i],pscore3[i]);
	//quality=ProQCA(filenames[i],m[0].ss);
	//printf("PARA %s ",filenames[i]);
	quality=ProQCA(filenames[i]);
	//printf("HERE2 %s %lf %lf\n",filenames[i],pscore1_5[i],pscore3[i]);
	proq_lg[i]=quality[0];
	proq_mx[i]=quality[1];
	//printf("%s %lf %lf\n",filenames[i],quality[0],quality[1]);
	free(quality);
	
	//printf("%-90s %7.4lf %7.4lf %7.4lf\n",filenames[i],score[i],proq_lg[i],proq_mx[i]);
      }
    //exit(1);
    // printf("Pairwise comparison\n");
    for(i=0;i<=files;i++)
      {
	//printf("%s ",filenames[i]);
	for(j=i+1;j<=files;j++)
	   {
	     if(strcmp(method[i],method[j])!=0)
	       {
		 LG=LGscore(filenames[i],filenames[j]);
		 //printf("%s %s %lf\n",filenames[i],filenames[j],LG);
		 //pairwise_LG_all[i][j]=LG;
		 //pairwise_LG_all[j][i]=LG;
		 sum_all[i]+=LG;
		 sum_all[j]+=LG;
		 count_all[i]++;
		 count_all[j]++;

		 if(rank[j]==1)
		  {
		    //printf("%d %d %s %s %d %d %lf\n",i,j,method[i],method[j],rank[i],rank[j],LG);
		    //pairwise_LG_first[i][j]=LG;
		    //pairwise_LG_first[j][i]=LG;
		    sum_first[i]+=LG;
		    //sum_first[j]+=LG;
		    count_first[i]++;
		    //count_first[j]++;
		  }
		 if(rank[i]==1)
		   {
		     //  printf("%d %d %s %s %d %d %lf\n",i,j,method[i],method[j],rank[i],rank[j],LG);
		     sum_first[j]+=LG;
		     count_first[j]++;
		   }


		    //}
		 //printf("%s %s %lf\n",filenames[i],filenames[j],LG);
	      }
	   }
      }
   
    //printf("%-90s  %-7s %-7s %-7s %-7s %-7s %-7s %-7s\n","Name","Pcons5","All","First","PROQmx","PROQlg","P1_5","P3");
    //printf("=======================================\n");
   for(i=0;i<=files;i++)
     {
       //printf("%d\n",i);
       jury_all=0;
       jury_first=0;
       if(count_all[i]!=0)
	 jury_all=sum_all[i]/count_all[i];
       if(count_first[i]!=0)
	 jury_first=sum_first[i]/count_first[i];
       
       if(jury_all>1) // only do the proq evaluation if the consensus is ok, lesson learned from a high scoring model in LB with high ProQ but low consensus
	 {
	   pcons5=0.5273*jury_all+0.1995*jury_first+0.6421*proq_mx[i]+0.2654*proq_lg[i]+0.3205*pscore1_5[i]+0.1346*pscore3[i];
	 }
       else
	 {
	   pcons5=jury_all;
	 }
       //printf("%-90s %7.4lf
       printf("%-100s %7.4lf %7.4lf %7.4lf %7.4lf %7.4lf %7.4lf %7.4lf\n",filenames[i],pcons5,jury_all,jury_first,proq_mx[i],proq_lg[i],pscore1_5[i],pscore3[i]);
     }
   //printf("=======================================\n");
}


void readcutoffs(cutoffs *cut,char *file)
{

  FILE *fp;
  int i=0;
  char temp1[100],temp2[100],temp3[100],temp4[100],temp5[100];
  // printf("File: %s\n",file);
  fp=fopen(file,"r");
  if(fp!=NULL)
    {
      //while(fscanf(fp,"%s %lf %lf %lf %lf",cut->methods[i].method,cut->methods[i].cutoffs90_3,cut->methods[i].cutoffs50_3,cut->methods[i].cutoffs50_1_5,cut->methods[i].cutoffs50_1_5)!=EOF)
      while(fscanf(fp,"%s %s %s %s %s",cut->methods[i].method,temp2,temp3,temp4,temp5)!=EOF)
      {
	//printf("%d %s %s %s %s %s\n",i,cut->methods[i].method,temp2,temp3,temp4,temp5);
	//	atof(temp1
	cut->methods[i].cutoffs90_3=atof(temp2);
	cut->methods[i].cutoffs50_3=atof(temp3);
	cut->methods[i].cutoffs90_1_5=atof(temp4);
	cut->methods[i].cutoffs50_1_5=atof(temp5);
      //printf("%d %s\n",i,cut->methods[i].method);
	//printf("%s %lf %lf %lf %lf\n",cut->methods[i].method,cut->methods[i].cutoffs90_3,cut->methods[i].cutoffs50_3,cut->methods[i].cutoffs50_1_5,cut->methods[i].cutoffs50_1_5);
	i++;
      }
    }
  cut->methods[i].method[0]='\0';
}
